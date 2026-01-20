suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DoubletFinder)
  library(data.table)
  library(Seurat)
  library(ggplot2)
})

join_layers_if_needed <- function(obj, verbose = TRUE) {
  stopifnot(inherits(obj, "Seurat"))
  # Seurat v5 multi-layer objects often need JoinLayers before QC/export
  if ("RNA" %in% names(obj@assays)) {
    ok <- FALSE
    n_layers <- NA_integer_
    try({
      n_layers <- length(Layers(obj[["RNA"]]))
      ok <- TRUE
    }, silent = TRUE)
    
    if (ok && !is.na(n_layers) && n_layers > 1) {
      if (verbose) message("JoinLayers: ", n_layers, " layers -> 1")
      obj <- JoinLayers(obj)
    }
  }
  obj
}

calc_qc_metrics <- function(obj, species = "human", do_cellcycle = TRUE, verbose = TRUE) {
  stopifnot(inherits(obj, "Seurat"))
  obj <- join_layers_if_needed(obj, verbose = verbose)
  
  DefaultAssay(obj) <- "RNA"
  
  mt_pattern <- if (species == "human") "^MT-" else "^mt-"
  rb_pattern <- if (species == "human") "^RP[SL]" else "^Rp[sl]"
  hb_pattern <- if (species == "human") "^HB[^(P)]" else "^Hb[^(p)]"
  
  if (verbose) message("Calculating percent.mt / percent.rb / percent.hb")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = rb_pattern)
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = hb_pattern)
  
  if (do_cellcycle) {
    if (verbose) message("Cell cycle scoring (NormalizeData -> CellCycleScoring)")
    obj <- NormalizeData(obj, verbose = FALSE)
    
    # Use updated cell cycle genes if available
    s.genes <- Seurat::cc.genes.updated.2019$s.genes
    g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
    if (species != "human") {
      s.genes <- stringr::str_to_title(s.genes)
      g2m.genes <- stringr::str_to_title(g2m.genes)
    }
    s.genes <- intersect(s.genes, rownames(obj))
    g2m.genes <- intersect(g2m.genes, rownames(obj))
    
    if (length(s.genes) > 10 && length(g2m.genes) > 10) {
      obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
    } else {
      if (verbose) message("Cell cycle genes not found enough; skip CellCycleScoring")
      obj$S.Score <- 0
      obj$G2M.Score <- 0
      obj$Phase <- "Unknown"
    }
  }
  
  obj
}
write_qc_report <- function(obj, out_dir, group_col = NULL) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  if (is.null(group_col)) {
    if ("Sample_Name" %in% colnames(obj@meta.data)) {
      group_col <- "Sample_Name"
    } else if ("GSM_ID" %in% colnames(obj@meta.data)) {
      group_col <- "GSM_ID"
    } else if ("orig.ident" %in% colnames(obj@meta.data)) {
      group_col <- "orig.ident"
    }
  }

  if (!is.null(group_col) && !group_col %in% colnames(obj@meta.data)) {
    group_col <- NULL
  }

  n_groups <- 1
  if (!is.null(group_col)) {
    n_groups <- length(unique(obj@meta.data[[group_col]]))
  }

  calc_width <- 5 + (n_groups * 0.4)
  calc_width <- max(8, min(calc_width, 50))

  features_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.hb")
  features_to_plot <- intersect(features_to_plot, colnames(obj@meta.data))

  pdf(file.path(out_dir, "qc_violin.pdf"), width = calc_width, height = 8)
  
  if (!is.null(group_col)) {
    p <- VlnPlot(obj, features = features_to_plot,
                 group.by = group_col, 
                 pt.size = 0, 
                 ncol = length(features_to_plot))
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + NoLegend()
    print(p)
  } else {
    print(VlnPlot(obj, features = features_to_plot, pt.size = 0, ncol = length(features_to_plot)))
  }
  dev.off()

  scatter_width <- max(10, calc_width * 0.8)
  pdf(file.path(out_dir, "qc_scatter.pdf"), width = scatter_width, height = 6)
  
  group_param <- if(!is.null(group_col)) group_col else NULL
  p1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = group_param)
  p2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = group_param)
  print(p1 + p2)
  dev.off()

  if (!is.null(group_col)) {
    tab <- obj@meta.data %>%
      group_by(across(all_of(group_col))) %>%
      summarise(
        nCells = n(),
        across(all_of(features_to_plot), median, .names = "med_{.col}"),
        .groups = "drop"
      )
  } else {
    tab <- obj@meta.data %>%
      summarise(
        nCells = n(),
        across(all_of(features_to_plot), median, .names = "med_{.col}")
      )
  }

  fwrite(as.data.table(tab), file.path(out_dir, "qc_summary.tsv"), sep = "\t")
  invisible(TRUE)
}
run_doubletfinder_one <- function(sub,
                                  npcs = 20,
                                  expected_rate = 0.075,
                                  resolution = 0.4,
                                  verbose = TRUE) {
  stopifnot(inherits(sub, "Seurat"))
  DefaultAssay(sub) <- "RNA"

  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, verbose = FALSE)
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, npcs = npcs, verbose = FALSE)
  sub <- FindNeighbors(sub, dims = 1:npcs, verbose = FALSE)
  sub <- FindClusters(sub, resolution = resolution, verbose = FALSE)

  df_ns <- asNamespace("DoubletFinder")

  sweep_fun <- NULL
  if (exists("paramSweep_v3", where = df_ns, inherits = FALSE)) {
    sweep_fun <- get("paramSweep_v3", envir = df_ns)
  } else if (exists("paramSweep", where = df_ns, inherits = FALSE)) {
    sweep_fun <- get("paramSweep", envir = df_ns)
  } else {
    stop("DoubletFinder paramSweep function not found in your installed package.")
  }

  sweep.res <- sweep_fun(sub, PCs = 1:npcs, sct = FALSE)

  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)

  best_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  best_pK <- as.numeric(as.character(best_pK))
  if (verbose) message("Best pK = ", best_pK)

  n_cells <- ncol(sub)
  nExp <- round(expected_rate * n_cells)
  homotypic.prop <- DoubletFinder::modelHomotypic(sub$seurat_clusters)
  nExp.adj <- round(nExp * (1 - homotypic.prop))
  if (verbose) message("nExp = ", nExp, " ; homotypic.prop = ", round(homotypic.prop, 3), " ; nExp.adj = ", nExp.adj)

  sub <- DoubletFinder::doubletFinder(
    seu = sub,
    PCs = 1:npcs,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp.adj,
    reuse.pANN = NULL,
    sct = FALSE
  )

  meta_cols <- colnames(sub@meta.data)
  pANN_col <- grep("^pANN", meta_cols, value = TRUE)
  cls_col  <- grep("^DF.classifications", meta_cols, value = TRUE)

  if (length(pANN_col) != 1 || length(cls_col) != 1) {
    stop("DoubletFinder output columns not found. pANN_col=", paste(pANN_col, collapse=","), " cls_col=", paste(cls_col, collapse=","))
  }

  sub$DF_pANN <- as.numeric(sub@meta.data[[pANN_col]])
  sub$DF_class <- as.character(sub@meta.data[[cls_col]])

  sub@meta.data[[pANN_col]] <- NULL
  sub@meta.data[[cls_col]] <- NULL

  sub
}
run_doubletfinder_by_sample <- function(obj,
                                        sample_col = "GSM_ID",
                                        out_dir,
                                        npcs = 20,
                                        expected_rate = 0.075,
                                        min_cells = 200,
                                        verbose = TRUE,
                                        write_global_table = TRUE) {
  stopifnot(inherits(obj, "Seurat"))
  if (!sample_col %in% colnames(obj@meta.data)) stop(paste0("sample_col not found: ", sample_col))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cell_dir <- file.path(out_dir, "cell_level")
  sum_dir  <- file.path(out_dir, "summary")
  dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sum_dir,  recursive = TRUE, showWarnings = FALSE)

  if (verbose) message("Joining layers globally...")
  obj <- JoinLayers(obj)

  obj$DF_class <- NA_character_
  obj$DF_pANN  <- NA_real_

  global_file <- file.path(cell_dir, "doublet_cells_all.tsv")
  if (write_global_table) {
    if (file.exists(global_file)) file.remove(global_file)
  }

  samples <- sort(unique(obj@meta.data[[sample_col]]))

  for (sid in samples) {
    if (verbose) message("Running DoubletFinder: ", sid)

    cells_use <- rownames(obj@meta.data)[obj@meta.data[[sample_col]] == sid]

    if (length(cells_use) < min_cells) {
      if (verbose) message("Skip (too small): ", sid, " n=", length(cells_use))

      obj@meta.data[cells_use, "DF_class"] <- "Singlet"
      obj@meta.data[cells_use, "DF_pANN"]  <- 0

      res_df <- data.frame(
        Barcode = cells_use,
        sample_id = sid,
        DF_class = "Singlet",
        DF_pANN = 0,
        nCount_RNA = obj@meta.data[cells_use, "nCount_RNA", drop = TRUE],
        nFeature_RNA = obj@meta.data[cells_use, "nFeature_RNA", drop = TRUE],
        percent.mt = if ("percent.mt" %in% colnames(obj@meta.data)) obj@meta.data[cells_use, "percent.mt", drop = TRUE] else NA_real_,
        stringsAsFactors = FALSE
      )

      out_one <- file.path(cell_dir, paste0("doublet_cells_", sid, ".tsv.gz"))
      data.table::fwrite(res_df, out_one, sep = "\t")

      if (write_global_table) {
        data.table::fwrite(res_df, global_file, sep = "\t", append = file.exists(global_file))
      }

      tab <- as.data.frame(table(res_df$DF_class, useNA = "ifany"))
      colnames(tab) <- c("DF_class", "nCells")
      tab$sample_id <- sid
      tab$expected_rate <- expected_rate
      data.table::fwrite(tab, file.path(sum_dir, paste0("doublet_table_", sid, ".tsv")), sep = "\t")

      rm(res_df, tab); gc()
      next
    }

    sub <- subset(obj, cells = cells_use)

    sub2 <- tryCatch(
      run_doubletfinder_one(sub, npcs = npcs, expected_rate = expected_rate, verbose = verbose),
      error = function(e) {
        message("DoubletFinder failed for ", sid, " : ", e$message)
        sub$DF_class <- "Unsure"
        sub$DF_pANN  <- NA_real_
        sub
      }
    )

    cells_back <- Cells(sub2)

    obj@meta.data[cells_back, "DF_class"] <- sub2$DF_class
    obj@meta.data[cells_back, "DF_pANN"]  <- sub2$DF_pANN

    res_df <- data.frame(
      Barcode = cells_back,
      sample_id = sid,
      DF_class = as.character(sub2$DF_class),
      DF_pANN = as.numeric(sub2$DF_pANN),
      nCount_RNA = sub2@meta.data[cells_back, "nCount_RNA", drop = TRUE],
      nFeature_RNA = sub2@meta.data[cells_back, "nFeature_RNA", drop = TRUE],
      percent.mt = if ("percent.mt" %in% colnames(sub2@meta.data)) sub2@meta.data[cells_back, "percent.mt", drop = TRUE] else NA_real_,
      seurat_clusters = if ("seurat_clusters" %in% colnames(sub2@meta.data)) as.character(sub2@meta.data[cells_back, "seurat_clusters", drop = TRUE]) else NA_character_,
      stringsAsFactors = FALSE
    )

    out_one <- file.path(cell_dir, paste0("doublet_cells_", sid, ".tsv.gz"))
    data.table::fwrite(res_df, out_one, sep = "\t")

    if (write_global_table) {
      data.table::fwrite(res_df, global_file, sep = "\t", append = file.exists(global_file))
    }

    tab <- as.data.frame(table(res_df$DF_class, useNA = "ifany"))
    colnames(tab) <- c("DF_class", "nCells")
    tab$sample_id <- sid
    tab$expected_rate <- expected_rate
    tab$nDoublet <- sum(res_df$DF_class == "Doublet", na.rm = TRUE)
    tab$doublet_frac <- tab$nDoublet / sum(tab$nCells)
    data.table::fwrite(tab, file.path(sum_dir, paste0("doublet_table_", sid, ".tsv")), sep = "\t")

    rm(sub, sub2, res_df, tab); gc()
  }

  if (verbose) {
    message("Done. DF columns in meta.data: ", paste(grep("^DF_", colnames(obj@meta.data), value = TRUE), collapse = ","))
  }

  obj
}


filter_singlets <- function(obj, keep_unsure = TRUE) {
  stopifnot(inherits(obj, "Seurat"))
  if (!"DF_class" %in% colnames(obj@meta.data)) stop("DF_class not found in meta.data")
  if (keep_unsure) {
    cells_keep <- rownames(obj@meta.data)[obj@meta.data$DF_class %in% c("Singlet", "Unsure")]
  } else {
    cells_keep <- rownames(obj@meta.data)[obj@meta.data$DF_class == "Singlet"]
  }
  subset(obj, cells = cells_keep)
}

export_mtx_bundle <- function(obj,
                              out_dir,
                              assay = "RNA",
                              layer = "counts",
                              join_layers = TRUE,
                              write_gz = FALSE) {
  stopifnot(inherits(obj, "Seurat"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  DefaultAssay(obj) <- assay

  obj_use <- obj

  is_v5 <- FALSE
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    is_v5 <- SeuratObject::Version(obj_use) >= "5.0"
  }

  if (is_v5 && join_layers) {
    obj_use <- JoinLayers(obj_use)
  }

  if (is_v5) {
    mat <- SeuratObject::LayerData(obj_use, assay = assay, layer = layer)
  } else {
    slot_use <- if (layer == "counts") "counts" else "data"
    mat <- Seurat::GetAssayData(obj_use, assay = assay, slot = slot_use)
  }

  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }

  mtx_path <- file.path(out_dir, "matrix.mtx")
  feat_path <- file.path(out_dir, "features.tsv")
  bc_path <- file.path(out_dir, "barcodes.tsv")
  meta_path <- file.path(out_dir, "metadata.csv")

  Matrix::writeMM(mat, file = mtx_path)

  write.table(
    data.frame(gene_id = rownames(mat), gene_name = rownames(mat)),
    file = feat_path,
    row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
  )

  write.table(
    colnames(mat),
    file = bc_path,
    row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
  )

  meta_out <- obj_use@meta.data
  meta_out <- meta_out[colnames(mat), , drop = FALSE]
  write.csv(meta_out, file = meta_path, row.names = TRUE)

  if (isTRUE(write_gz)) {
    gz1 <- gzfile(paste0(mtx_path, ".gz"), "wb"); close(gz1)
    gz2 <- gzfile(paste0(feat_path, ".gz"), "wb"); close(gz2)
    gz3 <- gzfile(paste0(bc_path, ".gz"), "wb"); close(gz3)

    con_in <- file(mtx_path, "rb"); con_out <- gzfile(paste0(mtx_path, ".gz"), "wb")
    writeBin(readBin(con_in, what = "raw", n = 1e9), con_out); close(con_in); close(con_out)
    file.remove(mtx_path)

    con_in <- file(feat_path, "rb"); con_out <- gzfile(paste0(feat_path, ".gz"), "wb")
    writeBin(readBin(con_in, what = "raw", n = 1e9), con_out); close(con_in); close(con_out)
    file.remove(feat_path)

    con_in <- file(bc_path, "rb"); con_out <- gzfile(paste0(bc_path, ".gz"), "wb")
    writeBin(readBin(con_in, what = "raw", n = 1e9), con_out); close(con_in); close(con_out)
    file.remove(bc_path)
  }

  invisible(TRUE)
}

