suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DoubletFinder)
  library(data.table)
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

run_doubletfinder_one <- function(sub,
                                  npcs = 20,
                                  expected_rate = 0.075,
                                  resolution = 0.4,
                                  verbose = TRUE) {
  stopifnot(inherits(sub, "Seurat"))
  DefaultAssay(sub) <- "RNA"
  
  # Minimal preprocessing required by DoubletFinder
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, verbose = FALSE)
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, npcs = npcs, verbose = FALSE)
  sub <- FindNeighbors(sub, dims = 1:npcs, verbose = FALSE)
  sub <- FindClusters(sub, resolution = resolution, verbose = FALSE)
  
  # Parameter sweep (compat across versions)
  if (exists("paramSweep_v3", where = asNamespace("DoubletFinder"), inherits = FALSE)) {
    sweep.res <- DoubletFinder::paramSweep_v3(sub, PCs = 1:npcs, sct = FALSE)
  } else {
    sweep.res <- DoubletFinder::paramSweep(sub, PCs = 1:npcs, sct = FALSE)
  }
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
  
  sub <- DoubletFinder::doubletFinder(sub,
                                      PCs = 1:npcs,
                                      pN = 0.25,
                                      pK = best_pK,
                                      nExp = nExp.adj,
                                      reuse.pANN = FALSE,
                                      sct = FALSE)
  
  # Extract DF columns (names vary by run)
  meta_cols <- colnames(sub@meta.data)
  pANN_col <- grep("^pANN", meta_cols, value = TRUE)
  cls_col  <- grep("^DF.classifications", meta_cols, value = TRUE)
  
  if (length(pANN_col) != 1 || length(cls_col) != 1) {
    stop("DoubletFinder output columns not found. pANN_col=", paste(pANN_col, collapse=","), " cls_col=", paste(cls_col, collapse=","))
  }
  
  sub$DF_pANN <- sub@meta.data[[pANN_col]]
  sub$DF_class <- sub@meta.data[[cls_col]]
  
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
                                        verbose = TRUE) {
  stopifnot(inherits(obj, "Seurat"))
  if (!sample_col %in% colnames(obj@meta.data)) stop("sample_col not found: ", sample_col)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialize result columns
  obj$DF_class <- NA_character_
  obj$DF_pANN <- NA_real_
  
  samples <- sort(unique(obj@meta.data[[sample_col]]))
  
  for (sid in samples) {
    if (verbose) message("Running DoubletFinder: ", sid)
    
    cells_use <- rownames(obj@meta.data)[obj@meta.data[[sample_col]] == sid]
    sub <- subset(obj, cells = cells_use)
    
    if (ncol(sub) < min_cells) {
      if (verbose) message("Skip (too small): ", sid, " n=", ncol(sub))
      obj@meta.data[cells_use, "DF_class"] <- "Singlet"
      obj@meta.data[cells_use, "DF_pANN"] <- 0
      next
    }
    
    sub2 <- tryCatch(
      run_doubletfinder_one(sub, npcs = npcs, expected_rate = expected_rate, verbose = verbose),
      error = function(e) {
        message("DoubletFinder failed for ", sid, " : ", e$message)
        sub$DF_class <- "Unsure"
        sub$DF_pANN <- NA_real_
        sub
      }
    )
    
    obj@meta.data[colnames(sub2), "DF_class"] <- sub2$DF_class
    obj@meta.data[colnames(sub2), "DF_pANN"] <- sub2$DF_pANN
    
    tab <- as.data.table(table(sub2$DF_class, useNA = "ifany"))
    fwrite(tab, file.path(out_dir, paste0("doublet_table_", sid, ".tsv")), sep = "\t")
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


export_mtx_bundle <- function(obj, out_dir, assay = "RNA", layer = "counts") {
  ensure_dirs(out_dir)
  DefaultAssay(obj) <- assay
  
  if ("SeuratObject" %in% loadedNamespaces()) {
    v <- as.character(utils::packageVersion("Seurat"))
  } else {
    v <- "unknown"
  }
  
  if (utils::packageVersion("Seurat") >= "5.0.0") {
    obj <- JoinLayers(obj)
    mat <- LayerData(obj, assay = assay, layer = layer)
  } else {
    mat <- GetAssayData(obj, assay = assay, slot = "counts")
  }
  
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  
  Matrix::writeMM(mat, file = file.path(out_dir, "matrix.mtx"))
  
  write.table(
    data.frame(rownames(mat), rownames(mat)),
    file = file.path(out_dir, "features.tsv"),
    row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
  )
  
  write.table(
    colnames(mat),
    file = file.path(out_dir, "barcodes.tsv"),
    row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
  )
  
  write.csv(obj@meta.data, file = file.path(out_dir, "metadata.csv"), row.names = TRUE)
  
  invisible(TRUE)
}
