#AuxFunctionscellbin
#zhaiyue xu
#


GenerateCellbinSampleData <- function(path_outs) {
  cat("--- 正在准备Cellbin单细胞级元数据 (v12, 修正面积计算) ---\n")
  
  path_outs_norm <- normalizePath(path_outs, winslash = "/", mustWork = TRUE)
  
  image_dir <- file.path(path_outs_norm, "spatial")
  lowres_image_path <- file.path(image_dir, "tissue_lowres_image.png")
  if (!file.exists(lowres_image_path)) stop("找不到低分辨率图像文件:", lowres_image_path)
  image_matrix <- png::readPNG(lowres_image_path)
  
  mappings_path   <- file.path(path_outs_norm, "barcode_mappings.parquet")
  coords_2um_path <- file.path(path_outs_norm, "binned_outputs/square_002um/spatial/tissue_positions.parquet")
  scales_path     <- file.path(path_outs_norm, "segmented_outputs", "spatial", "scalefactors_json.json")
  if (!all(file.exists(mappings_path, coords_2um_path, scales_path))) {
    stop("找不到必需的映射/坐标/比例因子文件。")
  }
  
  barcode_mappings <- arrow::read_parquet(mappings_path)
  coords_2um_df    <- arrow::read_parquet(coords_2um_path)
  colnames(coords_2um_df) <- c("barcode","tissue","row","col","imagerow","imagecol")
  scales <- rjson::fromJSON(file = scales_path)
  
  # 解析 bin 的物理边长（单位：µm），例如 square_002um -> 2
  bin_um <- as.numeric(gsub(".*square_([0-9]+)um.*", "\\1", coords_2um_path))
  if (is.na(bin_um) || bin_um <= 0) {
    warning("未能从路径解析 bin 大小，默认使用 2 µm。")
    bin_um <- 2
  }
  bin_area_um2 <- bin_um^2  # 2µm bin -> 4 µm²
  
  cat("  - 正在计算细胞元数据（包含面积与小格数）...\n")
  # 只取位于细胞内部、且有 cell_id 的 2µm 小格
  cell_bins <- barcode_mappings %>% dplyr::filter(in_cell == TRUE & !is.na(cell_id))
  
  # 关联 2µm 小格的空间坐标
  # 注意：barcode_mappings 里用于 join 的列名可能是 'square_002um'，根据你的文件结构保持一致
  join_key <- "square_002um"
  if (!(join_key %in% colnames(cell_bins))) stop("barcode_mappings 缺少列: ", join_key)
  cell_bins <- dplyr::left_join(
    cell_bins, coords_2um_df,
    by = setNames("barcode", join_key)  # square_002um = barcode
  )
  
  bcs <- cell_bins %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      imagerow = mean(imagerow, na.rm = TRUE),
      imagecol = mean(imagecol, na.rm = TRUE),
      tissue   = max(tissue, na.rm = TRUE),
      num_barcodes = dplyr::n(),        # 该细胞包含多少个 2µm 小格/条形码
      .groups = "drop"
    ) %>%
    dplyr::rename(barcode = cell_id) %>%
    dplyr::mutate(
      # 物理面积（µm²）= 小格数 × 每个小格面积
      area_um2 = num_barcodes * bin_area_um2,
      # 在低清图像坐标系中的坐标（用于和低清 png 对齐可视化）
      imagerow_scaled = imagerow * scales$tissue_lowres_scalef,
      imagecol_scaled = imagecol * scales$tissue_lowres_scalef,
      
      # —— 向后兼容的别名（如果你下游已有依赖）——
      # *如果你真的想看“像素面积”，应从分割 mask 统计像素数；这里不给出像素面积，避免误导*
      cell_area_bins = num_barcodes,    # 用于展示“按 bin 计”的面积规模
      cell_area_pixels = cell_area_bins # 兼容老代码：=bin 数（请在图注中说明含义）
    )
  
  cat("  - 正在加载 Space Ranger 预计算的分析结果...\n")
  cluster_path <- file.path(path_outs_norm, "segmented_outputs/analysis/clustering/gene_expression_graphclust/clusters.csv")
  umap_path    <- file.path(path_outs_norm, "segmented_outputs/analysis/umap/gene_expression_2_components/projection.csv")
  
  if (file.exists(cluster_path)) {
    clusters <- read.csv(cluster_path)
    if ("Cell.ID" %in% colnames(clusters)) clusters <- dplyr::rename(clusters, Barcode = Cell.ID)
    bcs <- dplyr::left_join(bcs, clusters, by = c("barcode" = "Barcode"))
  }
  if (file.exists(umap_path)) {
    umap_coords <- read.csv(umap_path)
    bcs <- dplyr::left_join(bcs, umap_coords, by = c("barcode" = "Barcode"))
  }
  
  bcs <- as.data.frame(bcs)
  bcs$height <- nrow(image_matrix)
  bcs$width  <- ncol(image_matrix)
  rownames(bcs) <- bcs$barcode
  
  return(list(bcs = bcs))
}

GenerateCellbinSeurat <- function(path_outs, sample_id) {
  tmp <- GenerateCellbinSampleData(path_outs = path_outs)
  bcs <- tmp$bcs
  bcs <- bcs %>% dplyr::filter(tissue == 1)
  
  h5_path <- file.path(path_outs, "segmented_outputs", "filtered_feature_cell_matrix.h5")
  counts  <- Read10X_h5(h5_path)
  
  common_cells   <- intersect(colnames(counts), rownames(bcs))
  counts_aligned <- counts[, common_cells, drop = FALSE]
  meta_aligned   <- bcs[common_cells, , drop = FALSE]
  
  sobj <- CreateSeuratObject(counts = counts_aligned, meta.data = meta_aligned,
                             assay = "Spatial", project = sample_id)
  cat("--- Cellbin Seurat对象创建成功！---\n")
  return(sobj)
}


# add_category_pattern
add_category_pattern <- function(
    obj,
    base_col      = "banksy_celltype_detailed",   # metadata 中的细胞类型列
    base_vec      = NULL,                         # 若非 Seurat，可直接传向量作为基础标签
    include_patterns = c("^Epithelial", "^Acinar"), # 命中规则（正则，或置 NULL）
    exclude_patterns = NULL,                      # 排除规则（正则）
    include_levels   = NULL,                      # 精确匹配的类型名（字符向量）
    exclude_levels   = NULL,                      # 精确排除的类型名
    include_cells    = NULL,                      # 直接按条形码加入命中集
    exclude_cells    = NULL,                      # 直接按条形码移出命中集
    label_in   = "Epithelial",
    label_out  = "Non-epithelial",
    new_col    = "EpiSimple",
    ignore_case = TRUE,
    overwrite   = TRUE,
    verbose     = TRUE
){
  # -------- 取条形码与基础标签 --------
  if (!is.null(base_vec)) {
    barcodes <- names(base_vec)
    if (is.null(barcodes)) barcodes <- colnames(obj)
    base <- as.character(base_vec)
  } else {
    stopifnot("Seurat" %in% class(obj))
    stopifnot(base_col %in% colnames(obj@meta.data))
    barcodes <- colnames(obj)
    base <- as.character(obj@meta.data[[base_col]])
  }
  if (length(base) != length(barcodes)) {
    stop("base 向量长度需与细胞数相同。")
  }
  
  hit <- rep(FALSE, length(base))
  
  # 正则包含
  if (!is.null(include_patterns) && length(include_patterns) > 0) {
    for (pat in include_patterns) {
      hit <- hit | grepl(pat, base, ignore.case = ignore_case)
    }
  }
  
  # 精确包含
  if (!is.null(include_levels) && length(include_levels) > 0) {
    hit <- hit | base %in% include_levels
  }
  
  # 条形码直接包含
  if (!is.null(include_cells) && length(include_cells) > 0) {
    hit <- hit | (barcodes %in% include_cells)
  }
  
  # 正则排除
  if (!is.null(exclude_patterns) && length(exclude_patterns) > 0) {
    for (pat in exclude_patterns) {
      hit <- hit & !grepl(pat, base, ignore.case = ignore_case)
    }
  }
  
  # 精确排除
  if (!is.null(exclude_levels) && length(exclude_levels) > 0) {
    hit <- hit & !(base %in% exclude_levels)
  }
  
  # 条形码直接排除
  if (!is.null(exclude_cells) && length(exclude_cells) > 0) {
    hit <- hit & !(barcodes %in% exclude_cells)
  }
  
  # -------- 生成新列 --------
  new_label <- ifelse(hit, label_in, label_out)
  new_label <- factor(new_label, levels = c(label_in, label_out))
  
  # 写回对象或返回向量
  if (!is.null(base_vec)) {
    # 只返回向量
    out <- new_label
    names(out) <- barcodes
  } else {
    if (!overwrite && new_col %in% colnames(obj@meta.data)) {
      stop(sprintf("列 '%s' 已存在，且 overwrite=FALSE。", new_col))
    }
    obj@meta.data[[new_col]] <- new_label
    out <- obj
  }
  
  if (verbose) {
    cat(sprintf("[add_category_pattern] %s 生成完成：\n", new_col))
    print(table(new_label))
  }
  
  return(out)
}
# =========================================================
# ScoreModules: single, smart scorer (SCT-first, else normalize)
# Methods: "auto","Seurat","UCell","zmean","ssGSEA","gsva","AUCell"
# Auto priority (fast-first): Seurat > UCell > zmean
# - If assay is NULL: pick SCT > Spatial > DefaultAssay(object)
# - If using SCT: by default do NOT normalize
# - If NOT SCT: by default normalize only when 'data' slot is empty
# - No built-in summary; use SummarizeScoresByGroup() below
# =========================================================

# ---- helpers ----
.pick_BPPARAM <- function(ncores = 1){
  if (!requireNamespace("BiocParallel", quietly = TRUE)) return(NULL)
  if (is.null(ncores) || ncores < 2) return(BiocParallel::SerialParam())
  if (.Platform$OS.type == "windows") {
    BiocParallel::SnowParam(workers = ncores, type = "SOCK", progressbar = FALSE)
  } else {
    BiocParallel::MulticoreParam(workers = ncores, progressbar = FALSE)
  }
}
.mm_scale <- function(x){
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

ScoreModules <- function(
    object,
    features_list,                           # named list: list(SigA=c(...), SigB=c(...))
    method        = c("auto","Seurat","UCell","zmean","ssGSEA","gsva","AUCell"),
    assay         = NULL,                    # auto-pick: SCT > Spatial > DefaultAssay
    slot_use      = "data",                  # usually "data" for RNA/SCT
    
    # normalization policy (kept simple):
    # - pre_normalize="auto": if SCT -> no normalize; else normalize only if 'data' empty
    # - "always": always NormalizeData(); "skip": never NormalizeData()
    pre_normalize     = c("auto","always","skip"),
    set_default_assay = TRUE,                # set DefaultAssay(object) <- assay
    
    # outputs & post-processing
    name_prefix   = "Score",
    add_to_meta   = TRUE,
    return_matrix = TRUE,
    post_scale    = c("none","z","minmax"),
    ncores        = 1,
    verbose       = TRUE,
    
    # auto order (fast-first)
    auto_order    = c("Seurat","UCell","zmean"),
    
    # overlap thresholds
    min_overlap_prop = 0.40,
    min_genes_abs    = 3
){
  stopifnot(inherits(object, "Seurat"))
  
  # ---- assay auto-pick (SCT first) ----
  if (is.null(assay)) {
    if ("SCT" %in% names(object@assays))       assay <- "SCT"
    else if ("Spatial" %in% names(object@assays)) assay <- "Spatial"
    else assay <- Seurat::DefaultAssay(object)
  }
  
  method        <- match.arg(method)
  post_scale    <- match.arg(post_scale)
  pre_normalize <- match.arg(pre_normalize)
  
  # ---- sanitize gene sets ----
  if (is.null(names(features_list)) || any(names(features_list) == "")) {
    names(features_list) <- paste0("Sig", seq_along(features_list))
  }
  if (any(duplicated(names(features_list)))) names(features_list) <- make.unique(names(features_list))
  features_list <- lapply(features_list, function(v) unique(as.character(v)))
  keep <- lengths(features_list) > 0
  if (!all(keep)) {
    warning("[ScoreModules] Dropping empty sets: ", paste(names(features_list)[!keep], collapse=", "))
    features_list <- features_list[keep]
  }
  if (length(features_list) == 0) stop("[ScoreModules] No valid gene sets supplied.")
  
  # ---- set assay / existence ----
  if (set_default_assay && !is.null(assay)) Seurat::DefaultAssay(object) <- assay
  if (!assay %in% names(object@assays)) stop("[ScoreModules] Assay not found: ", assay)
  
  # ---- normalization (SCT-safe) ----
  if (pre_normalize == "always") {
    object <- Seurat::NormalizeData(object, assay = assay, verbose = FALSE)
  } else if (pre_normalize == "auto") {
    if (!grepl("^SCT", assay, ignore.case = TRUE)) {
      # not SCT: normalize only if 'data' slot appears empty
      Etest <- tryCatch(Seurat::GetAssayData(object, assay = assay, slot = "data")[1:1, 1:1, drop = FALSE],
                        error = function(e) NULL)
      if (is.null(Etest) || all(Etest == 0)) {
        if (verbose) message("[ScoreModules] Normalizing assay '", assay, "' for 'data' slot...")
        object <- Seurat::NormalizeData(object, assay = assay, verbose = FALSE)
      }
    } # if SCT: do nothing
  } # if "skip": do nothing
  
  # ---- expression matrix ----
  E_data    <- Seurat::GetAssayData(object, assay = assay, slot = slot_use)  # genes x cells
  gene_pool <- rownames(E_data)
  
  # ---- overlap filtering ----
  .overlap <- function(feats, nm){
    min_found <- max(min_genes_abs, ceiling(length(feats) * min_overlap_prop))
    f <- intersect(feats, gene_pool)
    if (length(f) < min_found) {
      warning(sprintf("[ScoreModules] Too few genes for %s: %d/%d (min %d)",
                      nm, length(f), length(feats), min_found))
    }
    f
  }
  features_found <- lapply(names(features_list), function(nm) .overlap(features_list[[nm]], nm))
  names(features_found) <- names(features_list)
  nz <- vapply(features_found, length, 1L)
  if (any(nz == 0L)) {
    warning("[ScoreModules] Dropping 0-overlap sets: ", paste(names(features_found)[nz==0], collapse=", "))
    features_found <- features_found[nz > 0]
  }
  if (length(features_found) == 0) stop("[ScoreModules] No overlapped genes with assay matrix.")
  
  # ---- choose method (auto) ----
  if (method == "auto") {
    pick <- NULL
    for (m in auto_order) {
      if (m == "Seurat") { pick <- "Seurat"; break }
      if (m == "UCell"  && requireNamespace("UCell", quietly=TRUE)) { pick <- "UCell"; break }
      if (m == "zmean") { pick <- "zmean"; break }
    }
    if (is.null(pick)) pick <- "zmean"
    method <- pick
    if (verbose) message("[ScoreModules] auto -> using method: ", method)
  }
  
  ScoreMat <- NULL
  
  # ---- Seurat::AddModuleScore ----
  if (method == "Seurat") {
    if (verbose) message("[ScoreModules] Seurat::AddModuleScore (", length(features_found), " sets)")
    object <- Seurat::AddModuleScore(
      object, features = features_found, name = paste0(name_prefix, "_"),
      search = TRUE, assay = assay
    )
    base <- paste0(name_prefix, "_", names(features_found))
    base1 <- paste0(base, "1")
    cols <- if (all(base1 %in% colnames(object@meta.data))) base1 else base
    if (!all(cols %in% colnames(object@meta.data)))
      stop("[ScoreModules] Cannot locate AddModuleScore outputs.")
    ScoreMat <- as.matrix(object@meta.data[, cols, drop = FALSE])
    colnames(ScoreMat) <- base
    colnames(object@meta.data)[match(cols, colnames(object@meta.data))] <- base
  }
  
  # ---- UCell ----
  if (method == "UCell") {
    if (!requireNamespace("UCell", quietly = TRUE))
      stop("[ScoreModules] UCell not installed.")
    if (verbose) message("[ScoreModules] UCell (", length(features_found), " sets)")
    object <- UCell::AddModuleScore_UCell(
      object, features = features_found, name = NULL, assay = assay, ncores = ncores
    )
    uc <- names(features_found)
    if (!all(uc %in% colnames(object@meta.data))) {
      uc_alt <- paste0(uc, "_UCell")
      if (!all(uc_alt %in% colnames(object@meta.data)))
        stop("[ScoreModules] UCell output columns not found.")
      uc <- uc_alt
    }
    out <- paste0(name_prefix, "_", names(features_found))
    colnames(object@meta.data)[match(uc, colnames(object@meta.data))] <- out
    ScoreMat <- as.matrix(object@meta.data[, out, drop = FALSE])
  }
  
  # ---- z-mean (per-gene z then mean) ----
  if (method == "zmean") {
    if (verbose) message("[ScoreModules] z-mean (", length(features_found), " sets)")
    ScoreMat <- sapply(names(features_found), function(nm) {
      g  <- features_found[[nm]]
      mu <- Matrix::rowMeans(E_data[g, , drop = FALSE])
      sd <- matrixStats::rowSds(as.matrix(E_data[g, , drop = FALSE]))
      Z  <- (E_data[g, , drop = FALSE] - mu) / pmax(sd, 1e-6)
      Matrix::colMeans(Z)
    })
    if (is.null(dim(ScoreMat))) ScoreMat <- matrix(ScoreMat, ncol = 1)
    colnames(ScoreMat) <- paste0(name_prefix, "_", names(features_found))
  }
  
  # ---- GSVA / ssGSEA (new API) ----
  if (method %in% c("ssGSEA","gsva")) {
    if (!requireNamespace("GSVA", quietly = TRUE))
      stop("[ScoreModules] GSVA not installed.")
    if (verbose) message("[ScoreModules] GSVA::gsva(method='", method, "')")
    gset <- features_found
    names(gset) <- paste0(name_prefix, "_", names(gset))
    E_dense <- as.matrix(E_data)
    bpp <- .pick_BPPARAM(ncores)
    ScoreMat <- GSVA::gsva(
      expr           = E_dense,
      gset.idx.list  = gset,
      method         = if (method == "ssGSEA") "ssgsea" else "gsva",
      kcdf           = "Gaussian",
      min.sz         = 1, max.sz = Inf,
      BPPARAM        = bpp,
      ssgsea.norm    = TRUE,    # only for ssgsea
      mx.diff        = TRUE,    # only for gsva
      abs.ranking    = FALSE,
      verbose        = FALSE
    )
    ScoreMat <- t(ScoreMat)
    if (add_to_meta) for (j in colnames(ScoreMat)) object[[j]] <- as.numeric(ScoreMat[, j])
  }
  
  # ---- AUCell ----
  if (method == "AUCell") {
    if (!requireNamespace("AUCell", quietly = TRUE))
      stop("[ScoreModules] AUCell not installed.")
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
      stop("[ScoreModules] SummarizedExperiment not installed.")
    if (verbose) message("[ScoreModules] AUCell")
    ranks <- AUCell::buildRankings(as.matrix(E_data), nCores = ncores, plotStats = FALSE, verbose = FALSE)
    auc   <- AUCell::calcAUC(features_found, ranks, nCores = ncores, aucMaxRank = floor(0.05 * nrow(ranks)))
    ScoreMat <- t(SummarizedExperiment::assay(auc))
    colnames(ScoreMat) <- paste0(name_prefix, "_", colnames(ScoreMat))
  }
  
  # ---- post-hoc scaling ----
  if (post_scale != "none") {
    if (post_scale == "z")      ScoreMat <- scale(ScoreMat)
    else if (post_scale == "minmax") ScoreMat <- apply(ScoreMat, 2, .mm_scale)
    ScoreMat <- as.matrix(ScoreMat)
  }
  
  # ---- write to meta if needed (for branches that didn't already) ----
  if (add_to_meta && !(method %in% c("Seurat","UCell","ssGSEA","gsva"))) {
    for (j in colnames(ScoreMat)) object[[j]] <- as.numeric(ScoreMat[, j])
  }
  
  invisible(list(object = object, scores = if (return_matrix) ScoreMat else NULL,
                 method = method, assay = assay))
}


#' @title Integrate Hierarchical Annotations in a Seurat Object
#' @description Creates a final, unified annotation column by prioritizing a fine-grained
#'              annotation column and backfilling missing values with a coarse-grained one.
#'
#' @param object A Seurat object.
#' @param fine_col The name of the metadata column containing the detailed, manually 
#'                 curated annotations (may contain NAs).
#' @param coarse_col The name of the metadata column containing the initial, complete 
#'                   annotations that will be used for backfilling.
#' @param final_col The name of the new metadata column to be created with the 
#'                  integrated annotations.
#'
#' @return The Seurat object with the new, integrated annotation column.
#' @export
#'
#' @examples
#' # Assume 'my_seurat_object' has 'Seurat_L2_relabel' (with NAs) and 'Seurat_L2_label' (complete).
#' # my_seurat_object <- IntegrateAnnotations(
#' #   object = my_seurat_object,
#' #   fine_col = "Seurat_L2_relabel",
#' #   coarse_col = "Seurat_L2_label",
#' #   final_col = "cell_type_final"
#' # )
#' # table(my_seurat_object$cell_type_final)

IntegrateAnnotations <- function(
    object,
    fine_col,
    coarse_col,
    final_col = "cell_type_final"
) {
  
  # --- 1. Input Validation ---
  if (!inherits(object, "Seurat")) {
    stop("Error: 'object' must be a Seurat object.")
  }
  
  meta <- object@meta.data
  
  if (!fine_col %in% colnames(meta)) {
    stop("Error: 'fine_col' (", fine_col, ") not found in object metadata.")
  }
  if (!coarse_col %in% colnames(meta)) {
    stop("Error: 'coarse_col' (", coarse_col, ") not found in object metadata.")
  }
  
  # --- 2. Integration using dplyr::coalesce (Efficient and Clean) ---
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install the 'dplyr' package to use this function.")
  }
  
  # Coalesce returns the first non-NA value from the provided columns.
  # This perfectly matches the logic of "prioritize fine, backfill with coarse".
  final_annotations <- dplyr::coalesce(
    as.character(meta[[fine_col]]), 
    as.character(meta[[coarse_col]])
  )
  
  # --- 3. Add the new column back to the Seurat object ---
  
  # First, check if there are any remaining NAs (would happen if both columns had NAs)
  na_count <- sum(is.na(final_annotations))
  if (na_count > 0) {
    warning(
      "Warning: ", na_count, " cells still have NA annotations because they were NA in both '", 
      fine_col, "' and '", coarse_col, "'. Consider checking your input data."
    )
  }
  
  # Add as a factor for easier plotting and analysis
  object[[final_col]] <- factor(final_annotations)
  
  message("Successfully created new annotation column '", final_col, "'.")
  message("Summary of final annotations:")
  print(table(object[[final_col]], useNA = "ifany"))
  
  return(object)
}



# =========================================================
# SummarizeScoresByGroup (upgraded)
# - Accepts 1, 2, 3+ grouping columns (vector 'group_by')
# - Stable output names: score_mean, score_median, score_sd, score_min, score_max, score_n
# - Optional sorting on any output metric (e.g., sort_by = "score_mean")
# =========================================================
SummarizeScoresByGroup <- function(
    object,
    features,                         # character vector of meta score columns, e.g. c("Score_ECM","Score_T_EX")
    group_by,                         # character vector of 1..N grouping columns
    funs    = c("mean","median"),     # choose from: mean, median, sd, min, max, n
    na_rm   = TRUE,
    sort_by = NULL,                   # e.g. "score_mean", "score_median"; NULL = no sorting
    desc    = TRUE                    # sort direction if sort_by is set
){
  # --- checks ---
  stopifnot(inherits(object, "Seurat"))
  stopifnot(is.character(features), length(features) >= 1)
  stopifnot(is.character(group_by), length(group_by) >= 1)
  miss_g <- setdiff(group_by, colnames(object@meta.data))
  if (length(miss_g) > 0) stop("group_by not in meta.data: ", paste(miss_g, collapse=", "))
  miss_f <- setdiff(features,   colnames(object@meta.data))
  if (length(miss_f) > 0) stop("features not in meta.data: ", paste(miss_f, collapse=", "))
  
  if (!requireNamespace("dplyr", quietly = TRUE) || !requireNamespace("tidyr", quietly = TRUE))
    stop("Requires 'dplyr' and 'tidyr'.")
  
  # --- gather required columns only ---
  df <- cbind(object@meta.data[, group_by, drop = FALSE],
              object@meta.data[, features,   drop = FALSE])
  
  # --- long format: one row per (group..., signature, score) ---
  df_long <- tidyr::pivot_longer(
    df,
    cols      = tidyselect::all_of(features),
    names_to  = "signature",
    values_to = "score"
  )
  
  # --- build aggregation set with stable names ---
  agg <- list()
  if ("mean"   %in% funs) agg$mean   <- ~mean(.x, na.rm = na_rm)
  if ("median" %in% funs) agg$median <- ~stats::median(.x, na.rm = na_rm)
  if ("sd"     %in% funs) agg$sd     <- ~stats::sd(.x, na.rm = na_rm)
  if ("min"    %in% funs) agg$min    <- ~min(.x, na.rm = na_rm)
  if ("max"    %in% funs) agg$max    <- ~max(.x, na.rm = na_rm)
  if ("n"      %in% funs) agg$n      <- ~sum(!is.na(.x))
  if (length(agg) == 0) stop("No valid 'funs' selected.")
  
  # --- summarize by arbitrary number of grouping columns + signature ---
  out <- df_long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by)), .data$signature) |>
    dplyr::summarise(
      dplyr::across(
        .cols  = "score",
        .fns   = agg,
        .names = "score_{.fn}"     # ensures stable column names
      ),
      .groups = "drop"
    )
  
  # --- optional sorting ---
  if (!is.null(sort_by)) {
    if (!sort_by %in% colnames(out)) stop("sort_by not in output: ", sort_by)
    out <- if (desc) dplyr::arrange(out, dplyr::desc(.data[[sort_by]]))
    else       dplyr::arrange(out, .data[[sort_by]])
  }
  
  out
}

# 2.4 若缺少 first_w / second_w / margin：从权重矩阵补全
#     - 需要 Pass-1 保存的权重（行名 = cellbin，列名 = 参考大类）
fill_weights_if_missing <- function(object, weights_rds = weights_rds_path) {
  needed <- c("RCTD_first_w", "RCTD_second_w", "RCTD_margin")
  if (all(needed %in% colnames(object@meta.data))) return(object)
  if (!file.exists(weights_rds)) {
    message("未找到权重 RDS：", weights_rds, "；将仅基于 spot_class 做筛选。")
    object$RCTD_first_w  <- NA_real_
    object$RCTD_second_w <- NA_real_
    object$RCTD_margin   <- NA_real_
    return(object)
  }
  W <- readRDS(weights_rds)
  common <- intersect(rownames(object@meta.data), rownames(W))
  if (length(common) == 0) {
    message("权重矩阵与对象行名无交集；跳过补全。")
    object$RCTD_first_w  <- NA_real_
    object$RCTD_second_w <- NA_real_
    object$RCTD_margin   <- NA_real_
    return(object)
  }
  # 逐细胞拿 first/second 类型的权重
  first  <- as.character(object$RCTD_first[common])
  second <- as.character(object$RCTD_second[common])
  fw <- mapply(function(i, ct) if (!is.na(ct) && ct %in% colnames(W)) W[i, ct] else NA_real_,
               i = common, ct = first, SIMPLIFY = TRUE, USE.NAMES = FALSE)
  sw <- mapply(function(i, ct) if (!is.na(ct) && ct %in% colnames(W)) W[i, ct] else NA_real_,
               i = common, ct = second, SIMPLIFY = TRUE, USE.NAMES = FALSE)
  object$RCTD_first_w  <- NA_real_
  object$RCTD_second_w <- NA_real_
  object$RCTD_margin   <- NA_real_
  object$RCTD_first_w[common]  <- as.numeric(fw)
  object$RCTD_second_w[common] <- as.numeric(sw)
  object$RCTD_margin[common]   <- object$RCTD_first_w[common] - object$RCTD_second_w[common]
  object
}


CellbinQC <- function(object,
                      assay = "Spatial",
                      area_col = "area_um2",
                      mt_pattern = "^MT-",
                      ribo_pattern = "^RP[SL]",
                      # MAD 模式参数（仅在 use_mad=TRUE 时生效）
                      k_mad = 3,
                      cap_mt = 30,
                      # 面积分位（启用 use_area 时）
                      area_q = 0.05,
                      # 可选剔除
                      remove_RBC = FALSE,
                      remove_doublet = FALSE,
                      # 输出
                      make_plots = FALSE,
                      result_path = NULL,
                      sample_id = NULL,
                      verbose = TRUE,
                      # 开关
                      use_ribo = NULL,         # NULL=自动；TRUE/FALSE=强制
                      use_area = NULL,         # NULL=自动；TRUE/FALSE=强制
                      use_mt   = FALSE,        # ★ 默认不按线粒体比例过滤
                      # 极简模式：绝对地板/上限（推荐给 Visium HD/FFPE）
                      use_mad = FALSE,         # ★ 默认极简：不用 MAD，只用 floor / cap
                      floor_nCount = 40,       # FFPE/HD 常用 30–50（原始计数）
                      floor_nFeature = 20,     # FFPE/HD 常用 15–30（原始计数）
                      cap_ribo = 60,           # 极简模式下的核糖体上限（通常 use_ribo=FALSE）
                      # 诊断导出
                      write_diagnostics = TRUE # 输出规则统计 & 被删条码 TSV
){
  stopifnot(inherits(object, "Seurat"))
  old_assay <- Seurat::DefaultAssay(object)
  on.exit(Seurat::DefaultAssay(object) <- old_assay, add = TRUE)
  Seurat::DefaultAssay(object) <- assay
  
  # ------- counts（Seurat v5 兼容） -------
  counts <- tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = "counts"),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = "counts")
  )
  
  # ------- 基础指标（基于 counts 计算） -------
  sum_counts <- Matrix::colSums(counts)
  mt_idx     <- grep(mt_pattern,   rownames(counts))
  ribo_idx   <- grep(ribo_pattern, rownames(counts))
  object[["percent.mt"]]   <- if (length(mt_idx))   100 * Matrix::colSums(counts[mt_idx,  , drop=FALSE]) / pmax(1, sum_counts) else 0
  object[["percent.ribo"]] <- if (length(ribo_idx)) 100 * Matrix::colSums(counts[ribo_idx,, drop=FALSE]) / pmax(1, sum_counts) else 0
  
  # 直接从 counts 现算 nCount/nFeature，并回写 meta
  nCount   <- Matrix::colSums(counts)
  nFeature <- Matrix::colSums(counts > 0)
  object[[paste0("nCount_", assay)]]   <- nCount
  object[[paste0("nFeature_", assay)]] <- nFeature
  
  object$log10_nCount_   <- log10(nCount + 1)
  object$log10_nFeature_ <- log10(nFeature + 1)
  
  # ------- 实用函数 -------
  .mad_cut <- function(x, side = c("lower","upper"), k = 3) {
    side <- match.arg(side)
    med <- stats::median(x, na.rm = TRUE)
    md  <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
    if (is.na(md) || md == 0) md <- stats::IQR(x, na.rm = TRUE)/1.349
    if (side == "lower") med - k*md else med + k*md
  }
  .present_flag <- function(genes, mat) {
    g <- intersect(genes, rownames(mat))
    if (length(g) == 0) return(rep(FALSE, ncol(mat)))
    Matrix::colSums(mat[g, , drop = FALSE] > 0) > 0
  }
  
  # ------- RBC / Doublet 标记（如无则生成） -------
  if (!"rbc_flag" %in% colnames(object@meta.data)) {
    rbc_genes <- c("HBB","HBA1","HBA2")
    object$rbc_flag <- .present_flag(rbc_genes, counts)
  }
  if (!"doublet_flag" %in% colnames(object@meta.data)) {
    epi_genes <- c("EPCAM","KRT8","KRT18","KRT19","KRT7")
    imm_genes <- c("PTPRC","LYZ","LST1")
    epi_sig   <- .present_flag(epi_genes, counts)
    imm_sig   <- .present_flag(imm_genes, counts)
    object$doublet_flag <- epi_sig & imm_sig
  }
  
  mt   <- object$percent.mt
  ribo <- object$percent.ribo
  area <- if (area_col %in% colnames(object@meta.data)) object@meta.data[[area_col]] else rep(Inf, ncol(object))
  
  # ------- 阈值（MAD 模式 vs 极简模式） -------
  if (isTRUE(use_mad)) {
    thr_nCount_low   <- max(0, 10^(.mad_cut(log10(nCount   + 1), "lower", k_mad)))
    thr_nFeature_low <- max(0, 10^(.mad_cut(log10(nFeature + 1), "lower", k_mad)))
    thr_mt_high      <- max(0, min(cap_mt, .mad_cut(mt,   "upper", k_mad)))
    thr_ribo_high    <- max(0,          .mad_cut(ribo, "upper", k_mad))
  } else {
    thr_nCount_low   <- floor_nCount
    thr_nFeature_low <- floor_nFeature
    thr_mt_high      <- cap_mt
    thr_ribo_high    <- cap_ribo
  }
  # 应用“绝对地板”（即便 MAD 模式也不会低于地板）
  if (!is.null(floor_nCount))   thr_nCount_low   <- max(thr_nCount_low,   floor_nCount)
  if (!is.null(floor_nFeature)) thr_nFeature_low <- max(thr_nFeature_low, floor_nFeature)
  
  # —— 是否启用 ribo/area 规则（尊重用户；NULL 才自动）——
  if (is.null(use_ribo)) {
    use_ribo <- !(all(is.na(ribo)) ||
                    (stats::sd(ribo, na.rm=TRUE) < 1e-8) ||
                    (max(ribo, na.rm=TRUE) == 0))
  } else {
    use_ribo <- isTRUE(use_ribo)
  }
  if (is.null(use_area)) {
    use_area <- area_col %in% colnames(object@meta.data)
  } else {
    use_area <- isTRUE(use_area)
  }
  
  thr_area_low <- if (use_area) as.numeric(stats::quantile(area, probs = area_q, na.rm = TRUE)) else -Inf
  q_nCount10   <- as.numeric(stats::quantile(nCount,   0.10, na.rm = TRUE))
  q_nFeat10    <- as.numeric(stats::quantile(nFeature, 0.10, na.rm = TRUE))
  
  # ------- 规则（逐条）-------
  crit1 <- nCount   > thr_nCount_low
  crit2 <- nFeature > thr_nFeature_low
  crit3 <- if (isTRUE(use_mt)) (is.na(mt) | (mt < thr_mt_high)) else TRUE   # ★ 线粒体：默认跳过
  crit4 <- if (use_ribo) (is.na(ribo) | (ribo <= thr_ribo_high)) else TRUE
  crit5 <- if (use_area) !( (area < thr_area_low) & ( (nCount < q_nCount10) | (nFeature < q_nFeat10) ) ) else TRUE
  crit6 <- if (remove_RBC)      !object$rbc_flag      else TRUE
  crit7 <- if (remove_doublet)  !object$doublet_flag  else TRUE
  
  keep <- crit1 & crit2 & crit3 & crit4 & crit5 & crit6 & crit7
  names(keep) <- colnames(object)
  
  before_n <- ncol(object); kept_n <- sum(keep)
  
  # ------- 日志 -------
  if (verbose) {
    mode_tag <- if (isTRUE(use_mad)) "MAD" else "FLOOR"
    msg <- sprintf("[QC mode: %s] thresholds → nCount_%s>%.0f | nFeature_%s>%.0f | %s%s%s%s",
                   mode_tag, assay, thr_nCount_low, assay, thr_nFeature_low,
                   if (isTRUE(use_mt)) sprintf("mt<%.1f%%", thr_mt_high) else "mt:SKIPPED",
                   if (use_ribo) sprintf(" | ribo<=%.1f%%", thr_ribo_high) else " | ribo:SKIPPED",
                   if (use_area) sprintf(" | %s≥P%.0f", area_col, area_q*100) else "",
                   if (remove_RBC) " | -RBC" else "")
    message(msg)
    message(sprintf("QC: %d → %d (%.1f%% kept)", before_n, kept_n, 100*kept_n/before_n))
  }
  
  cells_keep <- names(keep)[keep]
  if (!length(cells_keep)) stop("QC 过严或指标退化：请放宽阈值/关闭 ribo/area/RBC/doublet（mt 已默认跳过）后重试。")
  object_qc <- subset(object, cells = cells_keep)
  
  # ------- 可选绘图 -------
  if (make_plots) {
    if (is.null(result_path)) stop("make_plots=TRUE 需要提供 result_path")
    if (!dir.exists(result_path)) dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
    sid <- if (is.null(sample_id)) "Sample" else sample_id
    if (exists("PlotQC_Summary")) {
      p_before <- PlotQC_Summary(object)
      ggplot2::ggsave(file.path(result_path, sprintf("%s_QC_before.png", sid)),
                      p_before, width = 15, height = 8, dpi = 300, bg = "white")
      p_after  <- PlotQC_Summary(object_qc)
      ggplot2::ggsave(file.path(result_path, sprintf("%s_QC_after.png",  sid)),
                      p_after,  width = 15, height = 8, dpi = 300, bg = "white")
    } else if (verbose) {
      message("未找到 PlotQC_Summary()，跳过 QC 汇总图。")
    }
  }
  
  # ------- 稳健统计（内置） -------
  # 1) 逐规则“贡献度”（不互斥）
  mask_detail <- data.frame(
    crit1_nCount   = crit1,
    crit2_nFeature = crit2,
    crit3_mt       = crit3,
    crit4_ribo     = crit4,
    crit5_area     = if (use_area) crit5 else TRUE,
    crit6_noRBC    = if (remove_RBC) crit6 else TRUE,
    crit7_noDoublet= if (remove_doublet) crit7 else TRUE,
    keep           = keep,
    row.names      = colnames(object),
    check.names    = FALSE
  )
  rule_cols <- c("crit1_nCount","crit2_nFeature","crit3_mt","crit4_ribo","crit5_area","crit6_noRBC","crit7_noDoublet")
  rule_cols <- intersect(rule_cols, colnames(mask_detail))
  rule_drop_counts <- colSums(!as.matrix(mask_detail[, rule_cols, drop=FALSE]))
  
  # 2) 唯一归因：按顺序找第一个失败的规则
  order_rules <- c("crit1_nCount","crit2_nFeature","crit3_mt","crit4_ribo","crit5_area","crit6_noRBC","crit7_noDoublet")
  order_rules <- intersect(order_rules, rule_cols)
  fail_reason <- apply(as.matrix(mask_detail[, order_rules, drop=FALSE]), 1, function(v){
    i <- which(!v)[1]; if (length(i)) order_rules[i] else "kept"
  })
  fail_reason_table <- sort(table(fail_reason), decreasing = TRUE)
  
  # 3) 被删条码清单（含关键指标）
  drop_idx <- which(!mask_detail$keep)
  barcodes_drop <- names(keep)[drop_idx]
  df_drop <- data.frame(
    barcode     = barcodes_drop,
    nCount      = nCount[barcodes_drop],
    nFeature    = nFeature[barcodes_drop],
    percent.mt  = mt[barcodes_drop],
    area_um2    = if ("area_um2" %in% colnames(object@meta.data)) object$area_um2[barcodes_drop] else NA,
    rbc_flag    = if ("rbc_flag" %in% colnames(object@meta.data)) object$rbc_flag[barcodes_drop] else NA,
    doublet_flag= if ("doublet_flag" %in% colnames(object@meta.data)) object$doublet_flag[barcodes_drop] else NA,
    fail_reason = as.character(fail_reason[barcodes_drop]),
    check.names = FALSE
  )
  
  # 可选落盘
  if (isTRUE(write_diagnostics) && !is.null(result_path)) {
    if (!dir.exists(result_path)) dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
    sid <- if (is.null(sample_id)) "Sample" else sample_id
    readr::write_tsv(as.data.frame(rule_drop_counts),
                     file.path(result_path, sprintf("%s_QC_rule_drop_counts.tsv", sid)),
                     col_names = FALSE)
    readr::write_tsv(as.data.frame(fail_reason_table),
                     file.path(result_path, sprintf("%s_QC_fail_reason_table.tsv", sid)),
                     col_names = FALSE)
    if (length(drop_idx)) {
      readr::write_tsv(df_drop,
                       file.path(result_path, sprintf("%s_QC_dropped.tsv", sid)))
    }
  }
  
  # ------- 汇总返回 -------
  thresholds <- c(
    nCount_low        = thr_nCount_low,
    nFeature_low      = thr_nFeature_low,
    percent.mt_high   = thr_mt_high,
    percent.ribo_high = thr_ribo_high,
    area_low          = if (use_area) thr_area_low else NA_real_
  )
  stats <- data.frame(
    before = before_n,
    after  = kept_n,
    kept_frac = kept_n / before_n,
    row.names = NULL
  )
  qc_summary <- list(
    rule_drop_counts  = rule_drop_counts,
    fail_reason_table = fail_reason_table
  )
  
  return(list(
    object_qc   = object_qc,
    keep        = keep,
    thresholds  = thresholds,
    stats       = stats,
    mask_detail = mask_detail,
    diagnostics = qc_summary,
    dropped     = df_drop
  ))
}


# ================================================================
# 子集重归一化→降维→分辨率扫描→选最佳 → 导出 Top markers
# （无打分，仅用每簇Top markers做注释证据）
# 依赖：Seurat, dplyr, ggplot2, readr, writexl(可选), patchwork
# 需要：PlotSpatialDistribution() 已在你的代码库中
# ================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(writexl)
})

# ---------- 小工具 ----------
.safe_genes <- function(genes, obj, assay="SCT"){
  intersect(genes, rownames(obj[[assay]]))
}

.present_flag_counts <- function(object, genes, assay="Spatial"){
  m <- Seurat::GetAssayData(object, assay=assay, slot="counts")
  g <- intersect(genes, rownames(m))
  if (!length(g)) return(rep(FALSE, ncol(object)))
  Matrix::colSums(m[g, , drop=FALSE] > 0) > 0
}

# ---------- 分辨率扫描&选择 ----------
# 扫描多个 resolution，记录 K 和小簇个数；返回“最合适”的 resolution
.pick_resolution_by_k <- function(obj, res_vec, k_target=c(4,10), min_cells=50, ident_prefix="SCT_snn_res."){
  tab <- lapply(res_vec, function(r){
    col <- paste0(ident_prefix, r)
    if (!col %in% colnames(obj@meta.data)) return(NULL)
    ks  <- table(obj@meta.data[[col]])
    data.frame(res=r,
               K=length(ks),
               min_size=min(ks),
               n_small=sum(ks < min_cells))
  })
  sweep <- bind_rows(tab)
  if (is.null(sweep) || !nrow(sweep)) stop("没有可用的分辨率扫描结果。")
  # 评分：优先无小簇 & K 在目标范围内；否则尽量逼近
  kmin <- k_target[1]; kmax <- k_target[2]; kmid <- mean(k_target)
  sweep <- sweep %>%
    mutate(
      fits_k = (K >= kmin & K <= kmax),
      # 越小越好
      score = 1e6 * (n_small > 0) +
        1e4 * (!fits_k) * pmin(abs(K - kmin), abs(K - kmax)) +
        1e2 * abs(K - kmid) +
        abs(res - 0.8)
    ) %>%
    arrange(score, res)
  list(best_res = sweep$res[1], sweep = sweep)
}
# ---- 升级版 subset_recluster：聚类 + （可选）自动合并相似簇 + 可选导出 ----
subset_recluster <- function(
    sub_obj,
    # 维度/特征
    dims_use = 1:30,
    variable.features.n = 6000,
    # 目标簇数/分辨率
    target_k = c(4, 10),
    min_cells_per_cluster = 50,
    res_candidates = seq(0.2, 2.0, by = 0.1),
    res_grid = NULL,
    ident_prefix = "SCT_snn_res.",
    # 邻居图
    k_param = 30,
    prune_snn = 1/15,
    # 回归项
    regress_vars = c("log10_nCount_Spatial","area_um2"),
    assay_in = "Spatial",
    # HVG/PC 钩子
    hvgs_transform_fn = NULL,
    exclude_pc_fn     = NULL,
    # —— 新增：是否执行合并 + 传参（会调用上面的 MergeSimilarClusters）
    do_merge   = FALSE,
    merge_params = list(                     # 这些只是默认，可在调用时覆盖
      out_col      = "merged_clusters",
      dims_use     = 1:20,
      p_val_thresh = 0.01,
      logfc_thresh = 0.25,
      min_pct      = 0.10,
      max_iter     = 20,
      assay        = "SCT",
      de_test      = "presto",
      features_use = NULL,                   # 如果 NULL，会在函数里自动用 head(HVG,2000)
      max_cells_per_ident = 3000,
      prioritize_nearest  = TRUE,
      verbose      = TRUE
    ),
    # 导出
    result_path = NULL,
    sample_id   = "Sample",
    export_loupe = TRUE,
    write_topN   = 50,
    do_umap = TRUE,
    do_neighbors = TRUE,
    verbose = TRUE
){
  stopifnot(inherits(sub_obj, "Seurat"))
  
  # ---------- 回归变量 ----------
  if ("nCount_Spatial" %in% colnames(sub_obj@meta.data)) {
    sub_obj$log10_nCount_Spatial <- log10(sub_obj$nCount_Spatial + 1)
  }
  if ("percent.mt" %in% regress_vars && !"percent.mt" %in% colnames(sub_obj@meta.data)) {
    lyr <- tryCatch({ if ("counts" %in% Layers(sub_obj[[assay_in]])) "counts" else "data" },
                    error = function(e) "counts")
    cts <- tryCatch(GetAssayData(sub_obj, assay = assay_in, layer = lyr),
                    error = function(e) GetAssayData(sub_obj, assay = assay_in, slot = "counts"))
    if (inherits(cts, "dgCMatrix")) {
      mt_idx <- grep("^MT-", rownames(cts))
      tot    <- Matrix::colSums(cts)
      sub_obj$percent.mt <- if (length(mt_idx)) 100 * Matrix::colSums(cts[mt_idx,, drop=FALSE]) / pmax(1, tot) else 0
    }
  }
  vars_reg <- intersect(regress_vars, colnames(sub_obj@meta.data))
  
  # ---------- SCT ----------
  use_glmGamPoi <- requireNamespace("glmGamPoi", quietly = TRUE)
  sub_obj <- SCTransform(
    sub_obj, assay = assay_in,
    vst.flavor = "v2",
    method = if (use_glmGamPoi) "glmGamPoi" else "poisson",
    vars.to.regress = vars_reg,
    variable.features.n = variable.features.n,
    return.only.var.genes = TRUE,
    verbose = FALSE
  )
  DefaultAssay(sub_obj) <- "SCT"
  
  # 可选：HVG 限额
  if (!is.null(hvgs_transform_fn)) {
    v0 <- VariableFeatures(sub_obj)
    v1 <- tryCatch(hvgs_transform_fn(v0), error = function(e) v0)
    if (is.character(v1) && length(v1) >= 200) VariableFeatures(sub_obj) <- unique(v1)
  }
  
  # ---------- PCA ----------
  sub_obj <- RunPCA(sub_obj, verbose = FALSE, features = VariableFeatures(sub_obj))
  
  # 可选：剔除“坏 PC”
  dims <- dims_use[dims_use <= ncol(Embeddings(sub_obj, "pca"))]
  bad_pcs <- integer(0)
  if (!is.null(exclude_pc_fn)) {
    bad_pcs <- tryCatch(exclude_pc_fn(sub_obj), error = function(e) integer(0))
    dims <- setdiff(dims, bad_pcs)
    if (length(dims) < 10) dims <- dims_use[1:min(10, length(dims_use))]
  }
  if (isTRUE(verbose)) message("dims used: ", paste(dims, collapse = ","), 
                               if (length(bad_pcs)) paste0(" | excluded PC: ", paste(bad_pcs, collapse=",")) else "")
  
  # ---------- 图 ----------
  if (isTRUE(do_umap)) {
    sub_obj <- RunUMAP(sub_obj, reduction = "pca", dims = dims, verbose = FALSE)
  }
  if (isTRUE(do_neighbors)) {
    sub_obj <- FindNeighbors(sub_obj, reduction = "pca", dims = dims,
                             k.param = k_param, prune.SNN = prune_snn, verbose = FALSE)
  }
  
  # ---------- 聚类 ----------
  res_vec <- if (!is.null(res_grid)) res_grid else res_candidates
  sub_obj <- FindClusters(sub_obj, resolution = res_vec, verbose = FALSE)
  
  picked <- .pick_resolution_by_k(
    sub_obj,
    res_vec   = res_vec,
    k_target  = target_k,
    min_cells = min_cells_per_cluster,
    ident_prefix = ident_prefix
  )
  best_res <- picked$best_res
  ident_col <- paste0(ident_prefix, best_res)
  Idents(sub_obj) <- ident_col
  
  if (isTRUE(verbose)){
    ks <- table(Idents(sub_obj))
    message(sprintf("Best res = %.1f | K=%d | min_size=%d",
                    best_res, length(ks), min(ks)))
  }
  
  out_files <- list()
  
  # ---------- 可选：导出（聚类原始版本） ----------
  if (!is.null(result_path)) {
    dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
    if (isTRUE(export_loupe)) {
      loupe_df <- data.frame(Barcode=colnames(sub_obj), Category=as.character(Idents(sub_obj)))
      fn <- file.path(result_path, sprintf("%s_res%.1f_loupe.csv", sample_id, best_res))
      readr::write_csv(loupe_df, fn); out_files$loupe_premerge <- fn
      message("[export_loupe_categories] 写出：", fn)
    }
    if (is.numeric(write_topN) && write_topN > 0) {
      mk_all <- FindAllMarkers(
        sub_obj, assay = DefaultAssay(sub_obj),
        only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.10, test.use = "wilcox"
      )
      mk_top <- mk_all %>%
        dplyr::group_by(cluster) %>%
        dplyr::arrange(dplyr::desc(avg_log2FC), dplyr::desc(pct.1 - pct.2)) %>%
        dplyr::slice_head(n = write_topN) %>%
        dplyr::ungroup()
      # 导出为 CSV 文件
      fn_csv <- file.path(result_path, sprintf("%s_res%.1f_loupe.csv", sample_id, best_res))
      readr::write_csv(loupe_df, fn_csv)
      message("[export_top_markers] 写出：", fn_csv)
      # 你可以选择将 markers_premerge 指向 CSV 或 XLSX 文件
      markers_premerge <- fn_csv 
    }
  }
  
  # ---------- （可选）自动合并相似簇 ----------
  merged_info <- list(enabled = FALSE, out_col = NA, merged_k = NA)
  if (isTRUE(do_merge)) {
    # 组装默认参数（自动填 features_use）
    mp <- modifyList(
      list(
        out_col      = "merged_clusters",
        dims_use     = 1:20,
        p_val_thresh = 0.01,
        logfc_thresh = 0.25,
        min_pct      = 0.10,
        max_iter     = 20,
        assay        = "SCT",
        de_test      = "presto",
        features_use = head(VariableFeatures(sub_obj), 2000),
        max_cells_per_ident = 3000,
        prioritize_nearest  = TRUE,
        verbose      = TRUE
      ),
      merge_params
    )
    sub_obj[[mp$out_col]] <- sub_obj[[ident_col]]
    # 调用合并
    sub_obj <- MergeSimilarClusters(
      obj        = sub_obj,
      id_col     = ident_col,
      out_col    = mp$out_col,
      dims_use   = mp$dims_use,
      p_val_thresh = mp$p_val_thresh,
      logfc_thresh = mp$logfc_thresh,
      min_pct      = mp$min_pct,
      max_iter     = mp$max_iter,
      assay        = mp$assay,
      de_test      = mp$de_test,
      features_use = mp$features_use,
      max_cells_per_ident = mp$max_cells_per_ident,
      prioritize_nearest  = mp$prioritize_nearest,
      verbose      = mp$verbose
    )
    Idents(sub_obj) <- mp$out_col
    merged_info$enabled  <- TRUE
    merged_info$out_col  <- mp$out_col
    merged_info$merged_k <- length(levels(sub_obj[[mp$out_col]][,1]))
    
    # 导出（合并后）
    if (!is.null(result_path)) {
      if (isTRUE(export_loupe)) {
        loupe_df2 <- data.frame(Barcode=colnames(sub_obj), Category=as.character(Idents(sub_obj)))
        fn2 <- file.path(result_path, sprintf("%s_merged_loupe.csv", sample_id))
        readr::write_csv(loupe_df2, fn2); out_files$loupe_merged <- fn2
        message("[export_loupe_categories] 写出：", fn2)
      }
      if (is.numeric(write_topN) && write_topN > 0) {
        mk_all2 <- FindAllMarkers(
          sub_obj, assay = DefaultAssay(sub_obj),
          only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.10, test.use = "wilcox"
        )
        mk_top2 <- mk_all2 %>%
          dplyr::group_by(cluster) %>%
          dplyr::arrange(dplyr::desc(avg_log2FC), dplyr::desc(pct.1 - pct.2)) %>%
          dplyr::slice_head(n = write_topN) %>%
          dplyr::ungroup()
        fn2 <- file.path(result_path, sprintf("%s_top%d_markers_merged.csv", sample_id, write_topN))
        readr::write_csv(mk_top2, fn2); 
        library(openxlsx)
        fn_xlsx <- file.path(result_path, sprintf("%s_res%d_markers_merged.xlsx", sample_id, write_topN))
        write.xlsx(mk_top2, fn2, rowNames = FALSE)
        out_files$markers_merged <- fn2 
        message("[export_top_markers] 写出：", fn_xlsx)
        message("[export_top_markers] 写出：", fn2)
      }
    }
  }
  
  return(list(
    obj = sub_obj,
    best_res = best_res,
    sweep = picked$sweep,
    excluded_pcs = bad_pcs,
    hvgs_n = length(VariableFeatures(sub_obj)),
    merged = merged_info,
    out_files = out_files
  ))
}


# ---------- 导出每簇 TopN markers ----------
export_top_markers <- function(obj, assay="SCT", n_top=50,
                               test.use="wilcox",
                               min.pct=0.05,
                               logfc.threshold=0.15){
  Idents(obj) <- Idents(obj) # ensure set
  mk <- FindAllMarkers(
    obj, assay=assay, only.pos=TRUE,
    test.use=test.use,
    min.pct=min.pct, logfc.threshold=logfc.threshold
  )
  topN <- mk %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group=TRUE) %>%
    slice_head(n = n_top) %>%
    ungroup()
  list(markers = mk, topN = topN)
}

# ---------- 运行所有大类子集 ----------
# 需要：object$coarse_byDEG 已存在（大类标签，来自你用 Top50 DEGs 做的粗分）
run_all_subsets <- function(object, outdir, sample_id="Sample",
                            params = list(
                              Epithelial   = list(dims=1:35, v=7000, k=c(6,12), min_cells=60),
                              CAF          = list(dims=1:30, v=6000, k=c(4,8),  min_cells=50),
                              Immune       = list(dims=1:25, v=5000, k=c(5,10), min_cells=20),
                              Endothelial  = list(dims=1:20, v=4000, k=c(2,5),  min_cells=20),
                              Pericyte_SMC = list(dims=1:20, v=4000, k=c(2,4),  min_cells=20),
                              AcinarCore   = list(dims=1:25, v=5000, k=c(3,6),  min_cells=40),
                              Acinar_IsletInterface = list(dims=1:25, v=5000, k=c(3,6),  min_cells=40),
                              Islet        = list(dims=1:20, v=4000, k=c(2,4),  min_cells=20)
                            ),
                            res_candidates = seq(0.2, 2.0, by=0.1),
                            save_excel = TRUE){
  
  stopifnot("coarse_byDEG" %in% colnames(object@meta.data))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  classes <- intersect(names(params), unique(object$coarse_byDEG))
  
  out_info <- list()
  
  for (cls in classes){
    cells <- colnames(object)[object$coarse_byDEG %in% cls]
    if (length(cells) < 30){
      message(sprintf("[%s] 细胞太少（%d），跳过。", cls, length(cells)))
      next
    }
    sub <- subset(object, cells = cells)
    p <- params[[cls]]
    
    rr <- subset_recluster(
      sub_obj = sub,
      dims_use = p$dims,
      variable.features.n = p$v,
      target_k = p$k,
      min_cells_per_cluster = p$min_cells,
      res_candidates = res_candidates,
      verbose = TRUE
    )
    sub <- rr$obj; best_res <- rr$best_res
    
    # 计算 Top markers（FFPE 友好阈值）
    mk <- export_top_markers(sub, assay="SCT",
                             n_top=50, test.use="wilcox",
                             min.pct=0.05, logfc.threshold=0.15)
    
    # 保存
    prefix <- file.path(outdir, paste0(sample_id, "_subset_", cls))
    readr::write_tsv(rr$sweep, paste0(prefix, "_resolution_sweep.tsv"))
    readr::write_tsv(mk$topN,  paste0(prefix, "_Top50.tsv"))
    saveRDS(sub,               paste0(prefix, ".rds"))
    if (isTRUE(save_excel)) {
      writexl::write_xlsx(mk$topN, paste0(prefix, "_Top50.xlsx"))
    }
    
    # 图（修复 ggsave 参数位置）
    p_u <- DimPlot(sub, reduction="umap", label=TRUE) + NoLegend()
    ggsave(filename=paste0(prefix, "_UMAP.png"), plot=p_u,
           width=7, height=6, dpi=300, bg="white")
    
    # 空间分布，用最佳分辨率列作为分组
    grp_col <- paste0("SCT_snn_res.", best_res)
    p_s <- PlotSpatialDistribution(sub, group_by = grp_col)
    ggsave(filename=paste0(prefix, "_Spatial.png"), plot=p_s,
           width=8, height=7, dpi=300, bg="white")
    
    message(sprintf("[%s] best resolution = %.1f; nClusters = %d",
                    cls, best_res, length(unique(Idents(sub)))))
    
    out_info[[cls]] <- list(
      best_res = best_res,
      n_clusters = length(unique(Idents(sub))),
      top50 = mk$topN,
      sweep = rr$sweep
    )
  }
  
  return(out_info)
}

# ----------（可选）按标记拉取稀有类再分析 ----------
# 例：把非浆细胞免疫拉出来再重聚类；或 Endothelial/Pericyte 太少时单独拉出
extract_by_markers <- function(object, markers, assay="Spatial", any_one=TRUE){
  flags <- lapply(markers, function(gs) .present_flag_counts(object, gs, assay=assay))
  if (any_one){
    keep <- Reduce("|", flags)
  } else {
    keep <- Reduce("&", flags)
  }
  subset(object, cells = colnames(object)[keep])
}

# 示例：拉取“非浆细胞免疫”再精细化
# imm_markers <- list(
#   Tcell = c("CD3D","CD3E","CD2"),
#   NK    = c("NKG7","GNLY","KLRD1"),
#   Mono  = c("LYZ","LST1","C1QA","C1QB","C1QC"),
#   Bcell = c("MS4A1","CD79A"),
#   DC    = c("CLEC9A","XCR1","CD1C","FCER1A","LILRA4","TCF4")
# )
# imm_sub <- extract_by_markers(object, imm_markers, assay="Spatial", any_one=TRUE)
# 然后把 imm_sub 当作 run_all_subsets 的一个输入类或单独跑 subset_recluster()




#export_loupe_categories
export_loupe_categories <- function(
    obj,
    meta_col,                           # 要导出的 metadata 列
    out_dir   = ".",                    # 导出文件夹（会自动创建）
    filename  = NULL,                   # 可自定义文件名；默认 "<meta_col>.loupe_category.csv"
    na_label  = "Unlabeled",            # 把 NA 替换为此标签，避免 Loupe 导入报错
    gzip      = FALSE,                  # TRUE 则写成 .csv.gz
    overwrite = TRUE,                   # FALSE 且目标已存在时将报错
    verbose   = TRUE
){
  stopifnot("Seurat" %in% class(obj))
  stopifnot(meta_col %in% colnames(obj@meta.data))
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # 默认文件名
  if (is.null(filename) || !nzchar(filename)) {
    filename <- sprintf("%s.loupe_category.csv", meta_col)
  }
  # 扩展名与 gzip 协调
  if (gzip) {
    if (!grepl("\\.csv(\\.gz)?$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".csv.gz")
    } else if (!grepl("\\.gz$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".gz")
    }
  } else {
    if (!grepl("\\.csv$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".csv")
    }
  }
  
  out_path <- file.path(out_dir, filename)
  if (file.exists(out_path) && !overwrite) {
    stop(sprintf("目标已存在且 overwrite=FALSE：%s", out_path))
  }
  
  # 组织数据
  vec <- as.character(obj@meta.data[[meta_col]])
  vec[is.na(vec)] <- na_label
  df <- data.frame(
    Barcode  = colnames(obj),
    Category = vec,
    check.names = FALSE
  )
  
  # 写文件
  if (gzip) {
    con <- gzfile(out_path, open = "wt")
    on.exit(close(con), add = TRUE)
    write.csv(df, con, row.names = FALSE)
  } else {
    write.csv(df, out_path, row.names = FALSE)
  }
  
  if (verbose) {
    msg <- sprintf("[export_loupe_categories] 已写出：%s（列：Barcode,Category）",
                   tryCatch(normalizePath(out_path), error = function(e) out_path))
    message(msg)
  }
  invisible(out_path)
}



## PlotSpatialCustom

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(forcats)
  library(rlang)
  ok_scattermore <- requireNamespace("scattermore", quietly = TRUE)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

PlotSpatialCustom <- function(
    metadata_df,
    color_by = NULL,
    coord_cols = c("imagecol", "imagerow"),
    spot_size = 0.1,
    alpha = 0.7,
    cmap = "viridis",
    title = NULL,
    palette = NULL
) {
  coord_candidates <- list(
    c("imagecol_scaled","imagerow_scaled"),
    c("imagecol","imagerow")
  )
  hit <- NULL
  if (all(coord_cols %in% colnames(metadata_df))) {
    hit <- coord_cols
  } else {
    for (cc in coord_candidates) {
      if (all(cc %in% colnames(metadata_df))) { hit <- cc; break }
    }
  }
  if (is.null(hit)) stop("找不到坐标列（尝试 imagecol_scaled/imagerow_scaled 与 imagecol/imagerow）")
  x_coord <- hit[1]; y_coord <- hit[2]
  
  aes_base <- aes(x = .data[[x_coord]], y = .data[[y_coord]])
  p <- ggplot(metadata_df, aes_base)
  
  draw_points <- function(p, mapping = NULL) {
    if (ok_scattermore) {
      if (is.null(mapping)) {
        p + scattermore::geom_scattermore(pointsize = spot_size * 20, alpha = alpha)
      } else {
        p + scattermore::geom_scattermore(mapping = mapping, pointsize = spot_size * 20, alpha = alpha)
      }
    } else {
      if (is.null(mapping)) {
        p + geom_point(size = spot_size, alpha = alpha, stroke = 0)
      } else {
        p + geom_point(mapping = mapping, size = spot_size, alpha = alpha, stroke = 0)
      }
    }
  }
  
  if (is.null(color_by)) {
    p <- draw_points(p)
  } else {
    if (!color_by %in% colnames(metadata_df)) stop("找不到 '", color_by, "' 列")
    is_discrete <- is.factor(metadata_df[[color_by]]) || is.character(metadata_df[[color_by]])
    
    p <- draw_points(p, mapping = aes(color = .data[[color_by]]))
    
    if (is_discrete) {
      if (!is.null(palette)) {
        p <- p + scale_color_manual(values = palette, name = color_by)
      } else {
        p <- p + viridis::scale_color_viridis(discrete = TRUE, name = color_by)
      }
    } else {
      p <- p + viridis::scale_color_viridis(option = cmap, name = color_by)
    }
  }
  
  p +
    coord_fixed() +
    scale_y_reverse() +                    # 与真实图像坐标方向一致
    theme_void() +
    labs(title = if (!is.null(title)) title else color_by) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

PlotClusterHighlight <- function(meta_df, cluster_col = "Cluster", cluster_id,
                                 highlight_color = "#E31A1C", bg_color = "grey85",
                                 spot_size = 0.08, title = NULL) {
  if (!(cluster_col %in% colnames(meta_df))) stop("缺少列：", cluster_col)
  df <- meta_df
  label_high <- paste0("Cluster ", cluster_id)
  
  df$.__cluster_hl__ <- ifelse(as.character(df[[cluster_col]]) == as.character(cluster_id),
                               label_high, "Others")
  df$.__cluster_hl__ <- factor(df$.__cluster_hl__, levels = c("Others", label_high))
  
  pal <- setNames(c(bg_color, highlight_color), c("Others", label_high))
  
  PlotSpatialCustom(
    metadata_df = df,
    color_by = ".__cluster_hl__",
    palette = pal,
    spot_size = spot_size,
    alpha = 0.9,
    title = if (is.null(title)) paste0(label_high, " (n=", sum(df[[cluster_col]] == cluster_id), ")") else title
  )
}

#PlotMetricDensityCompare
PlotMetricDensityCompare <- function(meta_df, cluster_col, cluster_id, metric, title = NULL) {
  if (!(metric %in% colnames(meta_df))) return(NULL)
  
  label_high <- paste0("Cluster ", cluster_id)
  
  df <- meta_df %>%
    dplyr::mutate(
      .group = ifelse(as.character(.data[[cluster_col]]) == as.character(cluster_id), label_high, "Others")
    ) %>%
    dplyr::mutate(.group = factor(.group, levels = c("Others", label_high)))
  
  fill_vals <- setNames(c("grey70", "#1F78B4"), c("Others", label_high))
  
  ggplot(df, aes(x = .data[[metric]], fill = .group)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = fill_vals) +
    theme_minimal() +
    labs(title = if (is.null(title)) metric else title, x = NULL, y = "Density") +
    theme(legend.position = "top")
}

PlotClusterComposition <- function(meta_df, cluster_col, cluster_id, celltype_col = "celltype_manual", top_n = 15) {
  if (!(celltype_col %in% colnames(meta_df))) return(NULL)
  df <- meta_df %>%
    dplyr::filter(as.character(.data[[cluster_col]]) == as.character(cluster_id)) %>%
    dplyr::mutate(!!celltype_col := forcats::fct_lump_n(as.factor(.data[[celltype_col]]), n = top_n)) %>%
    dplyr::count(.data[[celltype_col]], name = "n") %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::arrange(dplyr::desc(n))
  if (nrow(df) == 0) return(NULL)
  
  ggplot(df, aes(x = n, y = reorder(!!sym(celltype_col), n), fill = !!sym(celltype_col))) +
    geom_col(show.legend = FALSE) +
    viridis::scale_fill_viridis(discrete = TRUE) +
    theme_minimal() +
    labs(title = paste0("Celltype Composition — Cluster ", cluster_id),
         x = "Cells", y = NULL)
}



.sanitize_filename <- function(filename, replacement = "_") {
  # 移除非法字符: \ / : * ? " < > |
  sanitized <- gsub('[\\/:"*?<>|]', replacement, filename)
  # (可选) 将多个连续的替换符压缩为一个
  sanitized <- gsub(paste0(replacement, "{2,}"), replacement, sanitized)
  # (可选) 移除开头或结尾的替换符
  sanitized <- gsub(paste0("^", replacement, "|", replacement, "$"), "", sanitized)
  return(sanitized)
}


# --- 辅助函数：绘制单个簇的细胞类型组成（无修改）---
PlotClusterComposition <- function(meta_df, cluster_col, cluster_id, celltype_col = "celltype_manual", top_n = 15) {
  # 检查指定的细胞类型列是否存在
  if (!(celltype_col %in% colnames(meta_df))) {
    message(paste0("Warning: Column '", celltype_col, "' not found. Skipping composition plot."))
    return(NULL)
  }
  df <- meta_df %>%
    dplyr::filter(as.character(.data[[cluster_col]]) == as.character(cluster_id)) %>%
    dplyr::mutate(!!sym(celltype_col) := forcats::fct_lump_n(as.factor(.data[[celltype_col]]), n = top_n)) %>%
    dplyr::count(.data[[celltype_col]], name = "n") %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::arrange(dplyr::desc(n))
  if (nrow(df) == 0) return(NULL)
  ggplot(df, aes(x = n, y = reorder(!!sym(celltype_col), n), fill = !!sym(celltype_col))) +
    geom_col(show.legend = FALSE) +
    viridis::scale_fill_viridis(discrete = TRUE) +
    theme_minimal() +
    labs(title = paste0("Celltype Composition (Top ", top_n, ") — Cluster ", cluster_id),
         x = "Number of Cells", y = NULL)
}

# --- 辅助函数：密度对比图（无修改） ---
PlotMetricDensityCompare <- function(meta_df, cluster_col, cluster_id, metric, title = NULL) {
  if (!(metric %in% colnames(meta_df))) return(NULL)
  label_high <- paste0("Cluster ", cluster_id)
  df <- meta_df %>%
    dplyr::mutate(
      .group = ifelse(as.character(.data[[cluster_col]]) == as.character(cluster_id), label_high, "Others")
    ) %>%
    dplyr::mutate(.group = factor(.group, levels = c("Others", label_high)))
  fill_vals <- setNames(c("grey70", "#1F78B4"), c("Others", label_high))
  ggplot(df, aes(x = .data[[metric]], fill = .group)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = fill_vals) +
    theme_minimal() +
    labs(title = if (is.null(title)) metric else title, x = NULL, y = "Density") +
    theme(legend.position = "top")
}


#CreateCellbinOverview
CreateCellbinOverview <- function(object, result_path = "Analysis_Results/",
                                  cluster_col = "Cluster",
                                  celltype_col = "celltype_manual", 
                                  per_cluster_max_metrics = c("nCount_Spatial","nFeature_Spatial","area_um2","num_barcodes"),
                                  small_multiples_ncol = 5,
                                  paginate_n = 25) {
  
  cat("=== 创建Cellbin数据概览 (R版本) ===\n")
  fig_dir <- file.path(result_path, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  meta_df <- object@meta.data
  
  # ... (前面的代码无修改) ...
  # 为了简洁，这里省略
  if (!"nCount_Spatial" %in% colnames(meta_df) && !is.null(object$nCount_Spatial)) {
    meta_df$nCount_Spatial <- object$nCount_Spatial
  }
  if (!"nFeature_Spatial" %in% colnames(meta_df) && !is.null(object$nFeature_Spatial)) {
    meta_df$nFeature_Spatial <- object$nFeature_Spatial
  }
  
  # ---------- 图表1：基本指标空间分布 ----------
  cat("--- 正在创建基本指标空间分布图...\n")
  metrics_to_plot <- list(
    list(metric = "in_tissue",        title = "Tissue Coverage",     cmap = "plasma"),
    list(metric = "area_um2",         title = "CellBin Area (µm²)",  cmap = "viridis"),
    list(metric = "num_barcodes",     title = "# Sub-barcodes",      cmap = "cividis"),
    list(metric = "nCount_Spatial",   title = "Total Counts",        cmap = "plasma"),
    list(metric = "nFeature_Spatial", title = "Number of Genes",     cmap = "cividis")
  )
  
  overview_plots <- lapply(metrics_to_plot, function(m) {
    if (!(m$metric %in% colnames(meta_df))) return(NULL)
    PlotSpatialCustom(
      metadata_df = meta_df, color_by = m$metric,
      title = m$title, cmap = m$cmap, spot_size = 0.05
    )
  })
  overview_plots <- Filter(Negate(is.null), overview_plots)
  if (length(overview_plots) > 0) {
    p_overview <- wrap_plots(overview_plots, ncol = 2)
    # print(p_overview)  # 可以注释掉，加快速度
    ggsave(file.path(fig_dir, "cellbin_basic_metrics_R.png"),
           p_overview, width = 12, height = 10, dpi = 300, bg = "white")
  } else {
    message("提示：未找到可用的基本指标列，跳过基本指标空间图。")
  }
  
  # ---------- 图表2：GraphCluster 空间分布与计数 ----------
  if (cluster_col %in% colnames(meta_df)) {
    cat("--- 正在创建 GraphCluster 分析图...\n")
    tmp_df <- meta_df
    tmp_df[[cluster_col]] <- as.factor(tmp_df[[cluster_col]])
    
    p_cluster_spatial <- PlotSpatialCustom(
      tmp_df, color_by = cluster_col, title = "GraphCluster Spatial Distribution",
      spot_size = 0.05
    )
    
    p_cluster_counts <- tmp_df %>%
      count(.data[[cluster_col]]) %>%
      ggplot(aes(x = as.factor(.data[[cluster_col]]), y = n, fill = as.factor(.data[[cluster_col]]))) +
      geom_bar(stat = "identity", show.legend = FALSE) +
      geom_text(aes(label = n), vjust = -0.3) +
      viridis::scale_fill_viridis(discrete = TRUE) +
      theme_minimal() +
      labs(title = "GraphCluster Cell Counts", x = "Cluster ID", y = "Cell Count") +
      theme(axis.text.x = element_text(angle=45, hjust=1)) # 避免长名字重叠
    
    p_graphcluster <- p_cluster_spatial | p_cluster_counts
    # print(p_graphcluster) # 可以注释掉
    ggsave(file.path(fig_dir, "cellbin_graphcluster_analysis_R.png"),
           p_graphcluster, width = 14, height = 8, dpi = 300, bg = "white")
    
    # ---------- 图表2b：所有簇的小图矩阵（逐簇高亮；分页保存） ----------
    cat("--- 正在创建 GraphCluster 小图矩阵（逐簇高亮）...\n")
    clusters <- levels(as.factor(tmp_df[[cluster_col]]))
    
    hl_plots <- lapply(clusters, function(cl) {
      safe_cl_name <- .sanitize_filename(as.character(cl)) # 清洗标题中的名字
      PlotClusterHighlight(tmp_df, cluster_col = cluster_col, cluster_id = cl,
                           spot_size = 0.06,
                           title = paste0(safe_cl_name, " (n=", sum(tmp_df[[cluster_col]] == cl), ")"))
    })
    
    hl_plots <- Filter(Negate(is.null), hl_plots)
    if (length(hl_plots) > 0) {
      pages <- split(seq_along(hl_plots), ceiling(seq_along(hl_plots)/paginate_n))
      for (i in seq_along(pages)) {
        p_page <- wrap_plots(hl_plots[pages[[i]]], ncol = small_multiples_ncol) +
          plot_annotation(title = paste0("GraphCluster per-cluster highlights — page ", i))
        # print(p_page) # 可以注释掉
        ggsave(file.path(fig_dir, sprintf("graphcluster_highlights_page%02d.png", i)),
               p_page, width = 16, height = 10, dpi = 300, bg = "white")
      }
    }
    
    # ---------- 图表2c：逐簇详细分布（逐簇保存） ----------
    cat("--- 正在创建 GraphCluster 逐簇详细分布图（逐簇保存）...\n")
    metrics_exist <- intersect(per_cluster_max_metrics, colnames(meta_df))
    
    for (cl in clusters) {
      
      # --- 【【【 这是关键的修改点 】】】 ---
      # 1. 对cluster名字进行清洗，以创建安全的文件名
      safe_cl_name <- .sanitize_filename(as.character(cl))
      
      # 2. 构建最终的文件路径
      detail_plot_path <- file.path(fig_dir, sprintf("graphcluster_detail_cluster_%s.png", safe_cl_name))
      
      # (可选) 如果文件已存在，可以跳过以节省时间
      # if (file.exists(detail_plot_path)) {
      #   cat("  - Skipping Cluster '", as.character(cl), "' (plot already exists).\n")
      #   next
      # }
      cat("  - Processing Cluster '", as.character(cl), "' -> ", basename(detail_plot_path), "\n")
      
      p_left <- PlotClusterHighlight(tmp_df, cluster_col = cluster_col, cluster_id = cl,
                                     spot_size = 0.08)
      
      dens_plots <- lapply(metrics_exist, function(m) {
        PlotMetricDensityCompare(tmp_df, cluster_col = cluster_col, cluster_id = cl, metric = m)
      })
      dens_plots <- Filter(Negate(is.null), dens_plots)
      
      right_block <- if (length(dens_plots) > 0) {
        wrap_plots(dens_plots, ncol = 2)
      } else {
        ggplot() + theme_void() + labs(title = "No QC metrics available")
      }
      
      comp_plot <- PlotClusterComposition(tmp_df, cluster_col = cluster_col, cluster_id = cl,
                                          celltype_col = celltype_col, top_n = 15)
      
      if (!is.null(comp_plot)) {
        right_block <- right_block / comp_plot + plot_layout(heights = c(2,1))
      }
      
      p_detail <- p_left | right_block +
        plot_annotation(title = paste0("GraphCluster detail — Cluster ", cl))
      
      # print(p_detail) # 可以注释掉
      
      # 3. 使用安全的文件路径进行保存
      ggsave(detail_plot_path, p_detail, width = 16, height = 9, dpi = 300, bg = "white")
    }
    
  } else {
    message("提示：未找到列 '", cluster_col, "'，跳过 GraphCluster 相关图。")
  }
  
  # ... (图表3和打印摘要部分无修改) ...
  # 为了简洁，这里省略
  cat("--- 正在创建数据质量分布图...\n")
  hist_metrics <- c("nCount_Spatial","nFeature_Spatial",
                    "area_um2","num_barcodes","imagecol","imagerow")
  
  hist_plots <- lapply(hist_metrics, function(m) {
    if (!(m %in% colnames(meta_df))) return(NULL)
    ggplot(meta_df, aes(x = .data[[m]])) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      geom_vline(xintercept = stats::median(meta_df[[m]], na.rm = TRUE),
                 color = "red", linetype = "dashed", linewidth = 1) +
      labs(title = gsub("_", " ", m), x = NULL, y = "Frequency") +
      theme_minimal()
  })
  hist_plots <- Filter(Negate(is.null), hist_plots)
  
  if (length(hist_plots) > 0) {
    p_quality <- wrap_plots(hist_plots, ncol = 3)
    # print(p_quality) # 可以注释掉
    ggsave(file.path(fig_dir, "cellbin_data_quality_R.png"),
           p_quality, width = 15, height = 8, dpi = 300, bg = "white")
  } else {
    message("提示：未找到可用的质量指标列，跳过直方图组合。")
  }
  
  cat("\n=== Cellbin数据分析总结 (R版本) ===\n")
  cat(sprintf("总细胞数: %d\n", ncol(object)))
  if (cluster_col %in% colnames(meta_df)) {
    cat(sprintf("GraphCluster数量: %d\n", nlevels(as.factor(meta_df[[cluster_col]]))))
  } else {
    cat("GraphCluster数量: (未检测到 Cluster 列)\n")
  }
  
  invisible(object)
}


SummarizeNiches_For_CellBin <- function(
    object,
    niche_col,
    celltype_col,
    output_dir = "niche_summary",
    plots_per_page = 9,
    small_multiples_ncol = 3
){
  
  # --- 0. 载入依赖并进行参数检查 ---
  suppressPackageStartupMessages({
    library(Seurat); library(dplyr); library(patchwork); library(ggplot2)
  })
  
  required_cols <- c(niche_col, celltype_col, "imagecol", "imagerow")
  if (!all(required_cols %in% colnames(object@meta.data))) {
    stop("一个或多个必需的列在meta.data中找不到。需要: ", paste(required_cols, collapse=", "))
  }
  
  # 检查绘图函数是否存在于环境中
  required_funcs <- c("PlotClusterHighlight", "PlotClusterComposition")
  if (!all(sapply(required_funcs, exists, mode = "function"))) {
    stop("错误: 绘图所需的辅助函数 (PlotClusterHighlight, PlotClusterComposition) 未加载。请先 source('AuxFunctionscellbin.R')。")
  }
  
  message("--- 开始为 '", niche_col, "' 列进行系统性的Niche注解与总结 ---")
  
  # --- 1. 准备输出目录和数据 ---
  fig_dir <- file.path(output_dir, "figures_per_niche")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  
  meta <- object@meta.data
  meta[[niche_col]] <- factor(meta[[niche_col]]) # 确保是因子
  
  # --- 2. 计算核心指标 (您要求的数据驱动证据) ---
  message("  - (A, B) 正在计算细胞类型组成和香农多样性...")
  
  # A. 细胞类型组成表 (Composition)
  composition_table <- table(meta[[niche_col]], meta[[celltype_col]])
  
  # B. 香农多样性指数 (Diversity / Entropy)
  composition_prop <- prop.table(composition_table, margin = 1)
  shannon_entropy <- apply(composition_prop, 1, function(p) {
    p_filtered <- p[p > 0]
    -sum(p_filtered * log2(p_filtered))
  })
  
  # C. 整合为一个摘要数据框
  # 提取Top 3细胞类型及其比例
  top3_celltypes <- apply(composition_prop, 1, function(p) {
    top_indices <- order(p, decreasing = TRUE)[1:3]
    top_names <- names(p)[top_indices]
    top_props <- round(p[top_indices] * 100, 1)
    valid <- !is.na(top_names)
    paste(paste0(top_names[valid], " (", top_props[valid], "%)"), collapse = "; ")
  })
  
  # 创建最终的摘要数据框
  summary_df <- data.frame(
    niche_id = levels(meta[[niche_col]]),
    total_cells = as.numeric(table(meta[[niche_col]]))
  ) %>%
    mutate(
      shannon_entropy = shannon_entropy[niche_id],
      top_3_celltypes = top3_celltypes[niche_id]
    ) %>%
    arrange(desc(total_cells)) %>%
    relocate(niche_id, total_cells, shannon_entropy, top_3_celltypes)
  
  # 保存摘要表格
  summary_file_path <- file.path(output_dir, "niche_summary_table.csv")
  write.csv(summary_df, summary_file_path, row.names = FALSE)
  message("  - 摘要表格已成功保存至: ", normalizePath(summary_file_path, winslash="/"))
  
  # --- 3. 为每个Niche生成详细的可视化图 ---
  message("  - (C) 正在为每个Niche生成空间分布和组成图...")
  
  niche_ids <- levels(meta[[niche_col]])
  all_niche_plots <- list()
  
  for (id in niche_ids) {
    
    # 使用您自己的函数来创建高亮图
    p_spatial <- PlotClusterHighlight(
      meta_df = meta, 
      cluster_col = niche_col, 
      cluster_id = id,
      title = paste0("Niche ", id) # 保持标题简洁
    )
    
    # 使用您工作流中的函数来创建组成图
    p_composition <- PlotClusterComposition(
      meta_df = meta,
      cluster_col = niche_col,
      cluster_id = id,
      celltype_col = celltype_col,
      top_n = 10 # 只显示Top 10，保持图形清晰
    )
    
    # 组合图
    # 我们用一个空白图来占位，以确保对齐
    p_blank <- ggplot() + theme_void()
    
    # 使用 patchwork 将两个图组合起来
    if (!is.null(p_composition)) {
      # 用 patchwork 的 '|' 操作符将图左右并排
      combined_plot <- p_spatial | p_composition
    } else {
      # 如果没有组成图，就用空白图占位
      combined_plot <- p_spatial | p_blank
    }
    
    all_niche_plots[[as.character(id)]] <- combined_plot
  }
  
  # --- 4. 将所有Niche的图分页保存 ---
  message("  - 正在将所有摘要图分页保存...")
  
  if (length(all_niche_plots) > 0) {
    # 按照摘要表格的顺序对图进行排序 (从大到小)
    sorted_plots <- all_niche_plots[summary_df$niche_id]
    
    pages <- split(seq_along(sorted_plots), ceiling(seq_along(sorted_plots) / plots_per_page))
    
    for (i in seq_along(pages)) {
      p_page <- wrap_plots(sorted_plots[pages[[i]]], ncol = 1) + # 每行一个Niche，更清晰
        plot_annotation(title = paste0("Niche Summary — Page ", i))
      
      page_file_path <- file.path(fig_dir, sprintf("niche_summary_page%02d.png", i))
      
      # 动态调整页面高度
      n_rows_on_page <- length(pages[[i]])
      page_height <- max(6, 4 * n_rows_on_page) # 每个Niche分配4英寸高度
      
      ggsave(page_file_path, p_page, width = 12, height = page_height, 
             dpi = 150, bg = "white", limitsize = FALSE)
      message("    - 已保存: ", normalizePath(page_file_path, winslash="/"))
    }
  }
  
  message("--- Niche注解与总结流程完成！---")
  return(summary_df)
}


suppressPackageStartupMessages({
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
  if (!requireNamespace("stringr", quietly = TRUE))      install.packages("stringr")
  library(RColorBrewer); library(stringr)
})

suppressPackageStartupMessages({
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
  if (!requireNamespace("stringr", quietly = TRUE))      install.packages("stringr")
  library(dplyr); library(tibble); library(RColorBrewer); library(stringr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------ 0) 基础清洗（沿用并略增强） ------------
tidy_pdac_name <- function(x){
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "–|—", "-")
  x <- stringr::str_replace_all(x, "\\s+", " ")
  x
}

# ------------ 1) 解析“家族”（按前缀），并映射到 PDAC 宏类 ------------
# family: 标签第一个 "_" 之前的部分（若无 "_"，取整个词）
# macro : 将 family 归并到宏类（Epithelial/Immune/CAF/Acinar/Vascular/Endocrine/Neural/Interface/Other/Ambient）
parse_pdac_family_macro <- function(labels){
  lbl <- tidy_pdac_name(labels)
  
  # family by prefix
  fam <- ifelse(grepl("_", lbl, fixed = TRUE), sub("_.*$", "", lbl), lbl)
  fam_lc <- tolower(fam)
  
  # normalize families to canonical names
  fam_norm <- dplyr::case_when(
    fam_lc %in% c("tumor","tumour","cancer","epithelial") ~ "Tumor",
    fam_lc %in% c("adm","ductal","duct")                  ~ "ADM",
    grepl("^acinar", fam_lc)                               ~ "Acinar",
    fam_lc %in% c("caf","icaf","mycaf","fibroblast")      ~ "CAF",
    fam_lc %in% c("tam","macrophage","mono","myeloid")    ~ "TAM",
    fam_lc %in% c("plasma","b")                           ~ "Plasma",
    fam_lc %in% c("endocrine","islet")                    ~ "Endocrine",
    fam_lc %in% c("endothelial","vascular")               ~ "Endothelial",
    grepl("^perivascular", fam_lc)                         ~ "Perivascular",
    fam_lc %in% c("smoothmuscle","smoothmuscle ")         ~ "SmoothMuscle",
    fam_lc %in% c("schwann","neural","neuron")            ~ "Schwann",
    grepl("^interface", fam_lc)                            ~ "Interface",
    grepl("ambient", tolower(lbl))                         ~ "Ambient",
    TRUE                                                   ~ stringr::str_to_title(fam) # fallback
  )
  
  # macro mapping (families -> macro classes)
  macro <- dplyr::case_when(
    fam_norm %in% c("Tumor","ADM","Epithelial")        ~ "Epithelial",
    fam_norm %in% c("TAM","Plasma")                    ~ "Immune",
    fam_norm %in% c("CAF","Fibroblast")                ~ "CAF",
    fam_norm %in% c("Acinar")                          ~ "Acinar",
    fam_norm %in% c("Endocrine","Islet")               ~ "Endocrine",
    fam_norm %in% c("Endothelial","Perivascular")      ~ "Vascular",
    fam_norm %in% c("SmoothMuscle")                    ~ "Vascular",
    fam_norm %in% c("Schwann","Neural")                ~ "Neural",
    fam_norm %in% c("Interface")                       ~ "Interface",
    fam_norm %in% c("Ambient")                         ~ "Ambient",
    TRUE                                               ~ "Other"
  )
  
  tibble(label = labels, family = fam_norm, macro = macro)
}

# ------------ 2) 生成顺序调色函数（Brewer 连续梯度 + 截断 + 交错） ------------
# 安全的 zigzag 取色：1, n, 2, n-1, 3, ...
.zigzag <- function(n){
  left <- 1; right <- n
  idx <- integer(0)
  while (left <= right) {
    idx <- c(idx, left)
    if (left < right) idx <- c(idx, right)
    left  <- left + 1
    right <- right - 1
  }
  idx
}

.brewer_seq_fun <- function(pal_name, trim_light = 0.20, trim_dark = 0.85){
  force(pal_name); force(trim_light); force(trim_dark)
  function(n){
    n <- max(1L, as.integer(n))
    maxc <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
    grad <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(maxc, pal_name))(100)
    a <- max(1, floor(100*trim_light)); b <- min(100, ceiling(100*trim_dark))
    grad <- grad[a:b]
    if (n == 1) return(grad[round(length(grad)*0.7)])
    
    # 均匀抽 n 个点后用 zigzag 分散邻近色阶
    idx_base <- round(seq(1, length(grad), length.out = n))
    grad[ .zigzag(length(idx_base)) ] <- grad[idx_base]  # 先放置，再按 zigzag 顺序返回
    grad[seq_len(n)]
  }
}


# ------------ 3) PDAC 宏类到色系的默认映射（可自定义） ------------
pdac_macro_palettes <- function(){
  list(
    Epithelial = .brewer_seq_fun("Reds"),
    Immune     = .brewer_seq_fun("Purples"),
    CAF        = .brewer_seq_fun("Greens"),
    Acinar     = .brewer_seq_fun("Blues"),
    Vascular   = .brewer_seq_fun("PuBuGn"),
    Endocrine  = .brewer_seq_fun("YlOrBr"),
    Neural     = .brewer_seq_fun("YlGnBu"),
    Interface  = .brewer_seq_fun("Oranges"),
    Ambient    = .brewer_seq_fun("Greys"),
    Other      = .brewer_seq_fun("Greys")
  )
}

# ------------ 4) 主调色构建器：支持 mode = "macro" 或 "prefix" ------------
# mode="macro": 同一宏类共享一套色系；mode="prefix": 每个 family 一套色系
build_pdac_palette <- function(
    object,
    group_var = "Category_merged",
    mode      = c("macro","prefix"),
    trim_light = 0.20,
    trim_dark  = 0.85,
    locked_map = NULL,         # 命名向量：指定某些标签的固定颜色
    na_color   = "#BDBDBD",
    order_within = c("freq","alpha"),  # 组内排序：按频数或字母
    light_to_dark = TRUE,
    macro_pal = NULL           # 自定义 list(macro = function(n) colors)
){
  mode <- match.arg(mode)
  order_within <- match.arg(order_within)
  
  md <- object@meta.data
  stopifnot(group_var %in% colnames(md))
  lv <- levels(factor(md[[group_var]]))
  
  info <- parse_pdac_family_macro(lv)
  # ambient override: 若某标签含 Ambient，强行宏类=Ambient
  info$macro <- ifelse(grepl("Ambient", info$label, ignore.case = TRUE), "Ambient", info$macro)
  
  # 计算每个标签大小，用于 shades 排序
  size_tab <- as.data.frame(table(md[[group_var]]), stringsAsFactors = FALSE)
  colnames(size_tab) <- c("label","size")
  info <- dplyr::left_join(info, size_tab, by = "label")
  
  palette_bank <- macro_pal %||% pdac_macro_palettes()
  # 如果 mode="prefix"，用 family 作为“宏类”的键，从而每个家族独立用一套色系
  group_key <- if (mode == "macro") info$macro else info$family
  
  out <- setNames(rep(NA_character_, length(lv)), lv)
  for (g in unique(group_key)) {
    idx <- which(group_key == g)
    sub <- info[idx, , drop = FALSE]
    
    # 组内排序
    if (order_within == "freq") {
      sub <- sub %>% arrange(desc(size), label)
    } else {
      sub <- sub %>% arrange(label)
    }
    
    # 取对应色系生成 n 个色阶
    pal_fun <- palette_bank[[g]]
    if (is.null(pal_fun)) { # 若没有该键，用 Greys 兜底
      pal_fun <- .brewer_seq_fun("Greys", trim_light, trim_dark)
    } else {
      # 包装一下以传入 trim
      base_pal <- pal_fun; pal_fun <- function(n) base_pal(n)
    }
    
    cols <- pal_fun(nrow(sub))
    if (!light_to_dark) cols <- rev(cols)
    names(cols) <- sub$label
    out[sub$label] <- cols
  }
  
  # 固定色覆盖
  if (!is.null(locked_map)) {
    nm <- intersect(names(locked_map), names(out))
    out[nm] <- locked_map[nm]
  }
  
  # 处理 NA 组
  if (any(is.na(md[[group_var]]))) out <- c(out, setNames(na_color, NA))
  
  # 返回 mapping 与 key（方便检查与出图图例）
  key <- info %>% select(label, family, macro, size) %>% arrange(macro, family, desc(size), label)
  list(map = out, key = key)
}

# ------------ 5) 一行画图封装 ------------
PlotSpatialDistribution_pdac <- function(
    object,
    group_var   = "Category_merged",
    mode        = c("macro","prefix"),
    ptsize      = 1.6,
    show_legend = TRUE,
    ...
){
  pal <- build_pdac_palette(object, group_var = group_var, mode = mode)
  PlotSpatialDistribution(object,
                          group_by = group_var,
                          palette  = pal$map,
                          ptsize   = ptsize,
                          show_legend = show_legend,
                          ...)
}


# —— 3) 固定四大类的 Brewer 顺序色系；其余为辅助 —— 
# 用截断梯度 + 均匀取点 + 交错索引提升相邻可分度
make_niche_palette_fixed_core <- function(
    niche_levels,
    trim_light = 0.20,   # 建议 0.15–0.25：更亮更分明
    trim_dark  = 0.85,   # 裁掉最深端，避免发黑
    locked_map = NULL,
    custom_rules = NULL
){
  stopifnot(trim_light >= 0 && trim_light < trim_dark && trim_dark <= 1)
  
  lv_raw <- niche_levels
  lv <- tidy_niche_names(lv_raw)
  df <- assign_niche_groups(lv, custom_rules)
  
  base_palettes <- c(
    # —— 核心四类（固定，绝不改变）——
    Epithelial = "Reds",
    Acinar     = "Blues",
    CAF        = "Greens",
    Immune     = "Purples",
    # —— 辅助类（你不关心时就保持默认）——
    Interface  = "Oranges",
    Endocrine  = "YlOrBr",
    Neural     = "YlGnBu",
    Vascular   = "PuBuGn",
    Other      = "Greys"
  )
  
  brewer_seq <- function(n, pal, trim_light, trim_dark){
    maxc <- RColorBrewer::brewer.pal.info[pal, "maxcolors"]
    grad <- colorRampPalette(RColorBrewer::brewer.pal(maxc, pal))(100)
    a <- max(1, floor(100*trim_light))
    b <- min(100, ceiling(100*trim_dark))
    grad <- grad[a:b]
    if (n == 1) return(grad[round(length(grad)*0.7)])
    idx <- round(seq(1, length(grad), length.out = n))
    # 交错，防止相邻类别拿到相邻色阶
    interleave <- function(v){ as.vector(rbind(v[seq(1, length(v), 2)], v[seq(2, length(v), 2)]))[1:length(v)] }
    grad[interleave(idx)]
  }
  
  # 组内语义排序：Core/Interface/Inflamed 三类优先分开
  order_key <- function(s){
    s <- tolower(s)
    s <- str_replace_all(s, "homogeneous|core", "0_core")
    s <- str_replace_all(s, "interface|border|front|transition", "1_interface")
    s <- str_replace_all(s, "inflamed|stress|stressed", "2_inflamed")
    s
  }
  
  out <- setNames(rep(NA_character_, length(lv)), lv)
  for (g in unique(df$group)) {
    idx <- which(df$group == g)
    sub_lv <- lv[idx]
    palname <- base_palettes[[g]]; if (is.null(palname)) palname <- "Greys"
    sub_lv <- sub_lv[order(order_key(sub_lv))]
    cols   <- brewer_seq(length(sub_lv), palname, trim_light, trim_dark)
    names(cols) <- sub_lv
    out[sub_lv] <- cols
  }
  
  if (!is.null(locked_map)) {
    nm <- intersect(names(locked_map), names(out))
    out[nm] <- locked_map[nm]
  }
  
  pal <- out[ tidy_niche_names(niche_levels) ]
  names(pal) <- niche_levels
  pal
}



#PlotSpatialQC
PlotSpatialQC <- function(
    object,
    qc_feature,
    ptsize = 1.5
) {
  
  # 检查对象中是否存在该QC指标
  if (!qc_feature %in% colnames(object@meta.data)) {
    stop("在Seurat对象的meta.data中找不到指定的QC指标: '", qc_feature, "'")
  }
  
  # 准备绘图数据框
  plot_df <- object@meta.data[, c("imagecol_scaled", "imagerow_scaled", qc_feature)]
  
  # 绘图
  p <- ggplot(plot_df, aes(x = imagecol_scaled, y = -imagerow_scaled, color = .data[[qc_feature]])) +
    geom_scattermore(pointsize = ptsize, pixels = c(2000, 2000)) +
    scale_color_viridis_c() + # 使用viridis色盘
    coord_fixed() +
    theme_void() +
    labs(title = qc_feature, color = "Value")
  
  return(p)
}



# --------------------------------------------------------------------
# 函数名: PlotQC_Summary
# 功能: 一键生成包含所有核心QC指标的组合图。
# --------------------------------------------------------------------
PlotQC_Summary <- function(object) {
  
  # 检查是否已计算percent.mt
  if (!"percent.mt" %in% colnames(object@meta.data)) {
    cat("  - 注意: 'percent.mt' 未找到，正在计算...\n")
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  }
  
  # QC指标列表
  qc_features <- c("nCount_Spatial", "nFeature_Spatial", "percent.mt")
  
  # 生成小提琴图
  vln_plots <- lapply(qc_features, function(feat) PlotViolinQC(object, feat))
  
  # 生成空间热图
  spatial_plots <- lapply(qc_features, function(feat) PlotSpatialQC(object, feat))
  
  # 使用patchwork组合所有图像
  layout <- "
    AAABBBCCC
    DDDEEEFFF
  "
  # 我们只用一行小提琴图和一行空间图
  final_plot <- wrap_plots(vln_plots, ncol = 3) / wrap_plots(spatial_plots, ncol = 3) +
    plot_annotation(
      title = paste("Data Quality Control Summary for Sample:", unique(object$orig.ident)),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
  
  # 打印图像
  print(final_plot)
  
  # (可选) 返回组合图对象
  return(final_plot)
}




# Function used to generate images_tibble 
make_images_tibble<-function (PATH) 
{
  image <- get_spatial_files(PATH, "tissue_lowres_image")
  grobs <- grid::rasterGrob(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
  images_tibble <- tibble(Path = factor(PATH), grob = list(grobs), height = nrow(image), width = ncol(image))
  return(images_tibble)
}




# PlotViolinQC

PlotViolinQC <- function(
    object,
    qc_feature,
    group_by = "orig.ident" # 默认按样本分组
) {
  
  # 检查QC指标是否存在于meta.data中
  if (!qc_feature %in% colnames(object@meta.data)) {
    # 如果不存在，可能是nCount_Spatial等可以直接从assay计算的
    # 但在这里，我们假设它应该在meta.data里
    stop("QC指标 '", qc_feature, "' 在meta.data中未找到。")
  }
  
  p <- VlnPlot(
    object,
    features = qc_feature,
    group.by = group_by,
    pt.size = 0, # 不显示单个点
    
    # **核心修正：明确指定数据层**
    # 对于nCount/nFeature/percent.mt这类元数据特性，它们不依赖于layer
    # 但为了代码的统一性和未来的Seurat更新，我们可以养成好习惯
    # 当VlnPlot的feature是基因时，指定layer很重要
    # layer = "counts" 
    # 对于元数据列，Seurat v5的VlnPlot会自动从@meta.data获取，无需layer
  ) +
    NoLegend() +
    labs(title = qc_feature, x = "")
  
  return(p)
}




#PlotSpatialQC
PlotSpatialQC <- function(
    object,
    qc_feature,
    ptsize = 1.5
) {
  
  # 检查对象中是否存在该QC指标
  if (!qc_feature %in% colnames(object@meta.data)) {
    stop("在Seurat对象的meta.data中找不到指定的QC指标: '", qc_feature, "'")
  }
  
  # 准备绘图数据框
  plot_df <- object@meta.data[, c("imagecol_scaled", "imagerow_scaled", qc_feature)]
  
  # 绘图
  p <- ggplot(plot_df, aes(x = imagecol_scaled, y = -imagerow_scaled, color = .data[[qc_feature]])) +
    geom_scattermore(pointsize = ptsize, pixels = c(2000, 2000)) +
    scale_color_viridis_c() + # 使用viridis色盘
    coord_fixed() +
    theme_void() +
    labs(title = qc_feature, color = "Value")
  
  return(p)
}



# --------------------------------------------------------------------
# 函数名: PlotQC_Summary
# 功能: 一键生成包含所有核心QC指标的组合图。
# --------------------------------------------------------------------
PlotQC_Summary <- function(object) {
  
  # 检查是否已计算percent.mt
  if (!"percent.mt" %in% colnames(object@meta.data)) {
    cat("  - 注意: 'percent.mt' 未找到，正在计算...\n")
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  }
  
  # QC指标列表
  qc_features <- c("nCount_Spatial", "nFeature_Spatial", "percent.mt")
  
  # 生成小提琴图
  vln_plots <- lapply(qc_features, function(feat) PlotViolinQC(object, feat))
  
  # 生成空间热图
  spatial_plots <- lapply(qc_features, function(feat) PlotSpatialQC(object, feat))
  
  # 使用patchwork组合所有图像
  layout <- "
    AAABBBCCC
    DDDEEEFFF
  "
  # 我们只用一行小提琴图和一行空间图
  final_plot <- wrap_plots(vln_plots, ncol = 3) / wrap_plots(spatial_plots, ncol = 3) +
    plot_annotation(
      title = paste("Data Quality Control Summary for Sample:", unique(object$orig.ident)),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
  
  # 打印图像
  print(final_plot)
  
  # (可选) 返回组合图对象
  return(final_plot)
}




# Function used to generate images_tibble 
make_images_tibble<-function (PATH) 
{
  image <- get_spatial_files(PATH, "tissue_lowres_image")
  grobs <- grid::rasterGrob(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
  images_tibble <- tibble(Path = factor(PATH), grob = list(grobs), height = nrow(image), width = ncol(image))
  return(images_tibble)
}

#FindSpatialOutliers_DBSCAN

FindSpatialOutliers_DBSCAN <- function(
    object, 
    coord_cols = c("imagecol_scaled", "imagerow_scaled"),
    eps = 5,
    minPts = 150,
    new_col_name = "is_spatial_outlier",
    plot = TRUE # 新增参数，默认为TRUE，即自动绘图
) {
  library(scattermore)
  # --- 1. 参数检查 (不变) ---
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("请先安装 'dbscan' 包: install.packages('dbscan')")
  }
  if (!all(coord_cols %in% colnames(object@meta.data))) {
    stop("指定的坐标列 '", paste(coord_cols, collapse=", "), "' 在Seurat对象的meta.data中找不到。")
  }
  
  cat("--- 开始使用DBSCAN进行空间离群点检测 ---\n")
  cat("  - 邻域半径 (eps):", eps, "\n")
  cat("  - 核心点最小邻居数 (minPts):", minPts, "\n")
  
  # --- 2. 提取坐标并运行DBSCAN (不变) ---
  coords <- object@meta.data[, coord_cols]
  db_results <- dbscan::dbscan(coords, eps = eps, minPts = minPts)
  
  # --- 3. 解析结果并添加到Seurat对象 (不变) ---
  is_outlier_vec <- ifelse(db_results$cluster == 0, TRUE, FALSE)
  object <- AddMetaData(object, metadata = is_outlier_vec, col.name = new_col_name)
  
  n_outliers <- sum(is_outlier_vec)
  percent_outliers <- round(100 * n_outliers / ncol(object), 2)
  
  cat("--- 检测完成！---\n")
  cat("  - 找到", n_outliers, "个离群点 (占总数的", percent_outliers, "%)\n")
  cat("  - 结果已添加到meta.data的 '", new_col_name, "' 列中。\n")
  
  # --- **4. 新增的绘图模块** ---
  if (plot) {
    cat("--- 正在生成可视化图像... ---\n")
    
    # 准备绘图数据框
    plot_df <- data.frame(
      x_coord = coords[[coord_cols[1]]], # 使用变量，而不是写死的列名
      y_coord = coords[[coord_cols[2]]],
      status = ifelse(object@meta.data[[new_col_name]], "Outlier", "Inlier")
    )
    
    # 创建ggplot对象
    pDBSCAN <- ggplot(plot_df, aes(x = x_coord, y = y_coord)) +
      geom_scattermore(aes(color = status), pointsize = 1.5) + # 使用scattermore更快
      scale_color_manual(values = c("Inlier" = "grey70", "Outlier" = "firebrick")) +
      scale_y_reverse() +
      coord_fixed() +
      theme_void() +
      labs(
        title = "Spatial Outlier Detection (DBSCAN)",
        subtitle = paste(n_outliers, "outliers identified")
      ) +
      theme(legend.position = "bottom")
    
    # **直接打印图像**
    print(pDBSCAN)
  }
  
  # --- 5. 返回更新后的对象 (不变) ---
  return(object)
}




PlotSpatialDistribution <- function(
    object,
    group_by,
    coord_cols   = c("imagecol_scaled", "imagerow_scaled"),
    ptsize       = 1.2,
    palette      = "Set2",
    show_legend  = TRUE,
    spot_shape   = c("square","circle"),
    size_mode    = c("auto","constant"),
    area_col     = NULL,
    bin_um       = NULL,
    microns_per_pixel = NULL,
    scale_lowres      = NULL,
    # --- 新增参数 ---
    title = NULL,                 # 自定义主标题
    legend_title = NULL,          # 自定义图例标题
    legend_spacing_y = 0.3,       # ← 新增：图例项之间的“上下间距”（单位 lines）
    legend_key_height = 0.6,      # ← 新增：图例键（色块/点）高度（lines）
    legend_key_width  = 0.8,      # ← 新增：图例键（色块/点）宽度（lines）
    ...
){
  spot_shape <- match.arg(spot_shape)
  size_mode  <- match.arg(size_mode)
  
  # ---------- 1) 取绘图数据 ----------
  if (inherits(object, "Seurat")) {
    plot_df <- object@meta.data
    # 自动取 scales（优先 @misc，再 Images；若都没有尝试从坐标估算）
    mpp <- microns_per_pixel %||%
      object@misc$microns_per_pixel %||%
      object@misc$spatial_scales$microns_per_pixel %||%
      object@misc$scales$microns_per_pixel %||%
      tryCatch(Images(object)[[1]]@scale.factors$microns_per_pixel, error=function(e) NULL)
    
    s_low <- scale_lowres %||%
      object@misc$spatial_scales$tissue_lowres_scalef %||%
      object@misc$scales$tissue_lowres_scalef %||%
      tryCatch(Images(object)[[1]]@scale.factors$tissue_lowres_scalef, error=function(e) NULL)
    
    # 如仍缺失，基于 hires/lowres 坐标估算 s_low
    if (is.null(s_low) && all(c("imagecol_scaled","imagecol_hires") %in% colnames(plot_df))) {
      ratio <- plot_df$imagecol_scaled / plot_df$imagecol_hires
      s_low <- stats::median(ratio[is.finite(ratio)], na.rm = TRUE)
    }
    if (is.null(mpp)) stop("无法自动获取 microns_per_pixel；可手动传 microns_per_pixel= 数值。")
    if (is.null(s_low)) stop("无法自动获取 tissue_lowres_scalef；可手动传 scale_lowres= 数值。")
    
    # 如未指定 area_col，自动识别 Cellbin 的 area_um2
    if (is.null(area_col) && "area_um2" %in% names(plot_df)) area_col <- "area_um2"
    # 如未指定 bin_um，尝试从 misc 取（Visium HD）
    if (is.null(bin_um)) {
      bin_um <- object@misc$bin_um %||%
        object@misc$spatial_scales$bin_size_um %||% NULL
    }
  } else if (is.data.frame(object)) {
    plot_df <- object
    mpp  <- microns_per_pixel
    s_low<- scale_lowres
    if (is.null(mpp) || is.null(s_low))
      stop("传 data.frame 时请提供 microns_per_pixel 与 scale_lowres（标量或列名）。")
    if (is.null(area_col) && "area_um2" %in% names(plot_df)) area_col <- "area_um2"
  } else {
    stop("'object' 必须是 Seurat 或 data.frame")
  }
  
  # 坐标列检查
  if (!all(coord_cols %in% names(plot_df)))
    stop("缺少坐标列：", paste(coord_cols, collapse=", "))
  x_col <- coord_cols[1]; y_col <- coord_cols[2]
  
  # 分组列
  if (!group_by %in% names(plot_df)) stop("分组列不存在：", group_by)
  plot_df[[group_by]] <- as.factor(plot_df[[group_by]])
  levs <- levels(plot_df[[group_by]])
  
  # ---------- 2) 调色 ----------
  # 先确定图例标题
  legend_name <- if (is.null(legend_title)) group_by else legend_title
  
  if (is.character(palette) && length(palette)==1) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) suppressWarnings(install.packages("RColorBrewer"))
    pal_name <- palette
    if (pal_name %in% rownames(RColorBrewer::brewer.pal.info)) {
      maxc <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
      cols <- RColorBrewer::brewer.pal(min(maxc, max(length(levs), 3)), pal_name)
      plot_colors <- rep(cols, length.out = length(levs)); names(plot_colors) <- levs
    } else {
      stop("未知调色板：", palette)
    }
  } else if (is.character(palette)) {
    plot_colors <- palette; names(plot_colors) <- levs
  } else if (is.function(palette)) {
    tmp <- palette(); plot_colors <- rep(tmp, length.out=length(levs)); names(plot_colors) <- levs
  } else stop("无效 palette")
  
  # ---------- 3) 选择绘制分支 ----------
  library(ggplot2)
  flip_y <- function(y) -y
  
  use_real_tiles <- (spot_shape=="square" && size_mode=="auto" &&
                       ( (!is.null(area_col) && area_col %in% names(plot_df)) || !is.null(bin_um) ))
  
  if (use_real_tiles) {
    # 真实大小：计算每个点在低清像素中的边长
    if (!is.null(area_col) && area_col %in% names(plot_df)) {
      # Cellbin: area_um2 → 低清像素边长
      if (is.numeric(mpp))  mpp_vec <- rep(mpp,  nrow(plot_df)) else mpp_vec <- plot_df[[mpp]]
      if (is.numeric(s_low))sL_vec  <- rep(s_low, nrow(plot_df)) else sL_vec  <- plot_df[[s_low]]
      side <- sqrt(pmax(plot_df[[area_col]],0)) * (sL_vec / mpp_vec)
    } else {
      # 固定 bin（Visium HD）
      if (is.null(bin_um)) stop("size_mode='auto' 需要 area_col 或 bin_um")
      if (is.numeric(mpp))  mpp_vec <- rep(mpp,  nrow(plot_df)) else mpp_vec <- plot_df[[mpp]]
      if (is.numeric(s_low))sL_vec  <- rep(s_low, nrow(plot_df)) else sL_vec  <- plot_df[[s_low]]
      side <- (bin_um / mpp_vec) * s_low
    }
    
    # 计算矩形边框（注意 y 轴翻转）
    cx <- plot_df[[x_col]]; cy <- flip_y(plot_df[[y_col]])
    plot_df$.xmin <- cx - 0.5*side; plot_df$.xmax <- cx + 0.5*side
    plot_df$.ymin <- cy - 0.5*side; plot_df$.ymax <- cy + 0.5*side
    
    p <- ggplot(plot_df, aes(xmin=.xmin, xmax=.xmax, ymin=.ymin, ymax=.ymax,
                             fill=.data[[group_by]])) +
      geom_rect(color=NA, ...) +
      scale_fill_manual(values = plot_colors, name = legend_name)
    
  } else {
    # 符号大小模式（不按物理尺寸）
    if (requireNamespace("scattermore", quietly = TRUE) && spot_shape=="circle") {
      p <- ggplot(plot_df, aes(x=.data[[x_col]], y=flip_y(.data[[y_col]]),
                               color=.data[[group_by]])) +
        scattermore::geom_scattermore(pointsize = ptsize, pixels = c(2000,2000), ...)
    } else {
      shape_val <- if (spot_shape=="circle") 16 else 15
      p <- ggplot(plot_df, aes(x=.data[[x_col]], y=flip_y(.data[[y_col]]),
                               color=.data[[group_by]])) +
        geom_point(size=ptsize, shape=shape_val, ...)
    }
    p <- p + scale_color_manual(values = plot_colors, name = legend_name)
  }
  
  # ---------- 4) 主题 ----------
  plot_title <- if (is.null(title)) paste("Spatial Distribution of", group_by) else title
  p +
    coord_fixed() + theme_void() +
    labs(title = plot_title) +
    theme(
      plot.title      = ggplot2::element_text(hjust=0.5, face="bold", size=14),
      legend.text     = ggplot2::element_text(size=8),
      # ↓↓↓ 新增的三项可控图例布局参数 ↓↓↓
      legend.key.height = grid::unit(legend_key_height, "lines"),
      legend.key.width  = grid::unit(legend_key_width,  "lines"),
      legend.spacing.y  = grid::unit(legend_spacing_y,  "lines")
      # ↑↑↑ 新增的三项可控图例布局参数 ↑↑↑
    ) +
    { if (!show_legend) ggplot2::guides(color="none", fill="none") else NULL }
}



suppressPackageStartupMessages({
  library(Seurat); library(FNN); library(Matrix); library(dplyr)
})

# --------- helpers ----------
.coalesce_num <- function(a, b) as.numeric(ifelse(!is.na(a), a, b))

.get_xy <- function(obj){
  md <- obj@meta.data
  if (!all(c("imagecol","imagerow") %in% colnames(md))) {
    stop("meta.data 缺少 imagecol/imagerow")
  }
  x <- md$imagecol_scaled; y <- md$imagerow_scaled
  if (is.null(x) || all(is.na(x))) x <- md$imagecol
  if (is.null(y) || all(is.na(y))) y <- md$imagerow
  x <- as.numeric(x); y <- as.numeric(y)
  if (any(!is.finite(x) | !is.finite(y))) stop("坐标含 NA/非数值，请检查 imagecol/_scaled 与 imagerow/_scaled")
  cbind(x=x, y=y)
}

.has_data_slot <- function(obj, assay){
  M <- tryCatch(GetAssayData(obj, assay=assay, slot="data"), error=function(e) NULL)
  if (is.null(M)) M <- tryCatch(GetAssayData(obj, assay=assay, layer="data"), error=function(e) NULL)
  !is.null(M) && nrow(M) > 0 && ncol(M) > 0
}

.choose_assay_slot <- function(object,
                               assay=NULL, slot=NULL,
                               prefer="Spatial"){  # 默认更“显色、直观”
  # 显式指定 -> 完全尊重（不做 Normalize）
  if (!is.null(assay) || !is.null(slot)) {
    return(list(
      assay = if (is.null(assay)) DefaultAssay(object) else assay,
      slot  = if (is.null(slot))  "data" else slot,
      normalize_if_needed = FALSE
    ))
  }
  
  prefer <- tolower(as.character(prefer)[1])
  if (!prefer %in% c("spatial","sct","auto")) prefer <- "spatial"
  
  assay_names <- names(object@assays)
  has_spatial <- "Spatial" %in% assay_names
  has_sct     <- "SCT"     %in% assay_names
  
  sel <- NULL
  if (prefer == "spatial") {
    if (has_spatial) sel <- list(assay="Spatial", slot="data", normalize_if_needed=FALSE)
    else if (has_sct) sel <- list(assay="SCT", slot="data", normalize_if_needed=FALSE)
  } else if (prefer == "sct") {
    if (has_sct) sel <- list(assay="SCT", slot="data", normalize_if_needed=FALSE)
    else if (has_spatial) sel <- list(assay="Spatial", slot="data", normalize_if_needed=FALSE)
  } else { # auto
    if (has_sct) sel <- list(assay="SCT", slot="data", normalize_if_needed=FALSE)
    else if (has_spatial) sel <- list(assay="Spatial", slot="data", normalize_if_needed=FALSE)
  }
  
  if (!is.null(sel)) return(sel)
  
  # 兜底：用默认 assay 的 data（不 Normalize）
  a0 <- DefaultAssay(object)
  return(list(assay=a0, slot="data", normalize_if_needed=FALSE))
}


.get_mat <- function(obj, assay, slot="data", normalize_if_needed=FALSE){
  DefaultAssay(obj) <- assay
  if (slot == "data" && normalize_if_needed) {
    if (assay == "SCT") stop("SCT 不应 NormalizeData；选择逻辑错误。")
    obj <- NormalizeData(obj, assay=assay, verbose=FALSE)
  }
  M <- tryCatch(GetAssayData(obj, assay=assay, slot=slot), error=function(e) NULL)
  if (is.null(M)) M <- tryCatch(GetAssayData(obj, assay=assay, layer=slot), error=function(e) NULL)
  if (is.null(M)) stop("无法读取矩阵：assay=", assay, " slot/layer=", slot)
  if (!inherits(M, "dgCMatrix")) M <- as(M, "dgCMatrix")
  list(obj=obj, M=M)
}

# --------- upgraded knn_smooth ----------
knn_smooth <- function(
    object,
    meta_cols = NULL,     # 例如 c("Score_ECM")
    genes     = NULL,     # 例如 c("ITGAV","ITGB1")
    k         = 8,
    method    = "gaussian",  # "mean" 或 "gaussian"
    sigma     = NULL,        # gaussian 时可空：自动用 1-NN 中位距离
    assay     = NULL,        # 默认自动选择；若显式给出将严格使用且不 Normalize
    slot      = NULL,        # 一般为 "data"
    prefer    = "Spatial",   # 默认更显色：Spatial 优先；可改 "SCT" / "auto"
    out_suffix = NULL        # 默认 "_knn{k}" 或 "_wknn{k}"
){
  if (!requireNamespace("FNN", quietly = TRUE)) stop("需要 FNN 包：install.packages('FNN')")
  # method 归一
  method <- tolower(as.character(method)[1])
  if (!method %in% c("mean","gaussian")) method <- "gaussian"
  # 选择 assay/slot
  choice <- .choose_assay_slot(object, assay=assay, slot=slot, prefer=prefer)
  mm <- .get_mat(object, assay=choice$assay, slot=choice$slot,
                 normalize_if_needed = choice$normalize_if_needed)
  object <- mm$obj; M <- mm$M
  
  # kNN
  xy <- .get_xy(object)
  nn <- FNN::get.knn(xy, k=k)
  if (is.null(sigma) && method == "gaussian") {
    sigma <- stats::median(nn$nn.dist[,1], na.rm=TRUE)
    if (!is.finite(sigma) || sigma <= 0) sigma <- 1
  }
  get_weights <- if (method == "mean") {
    function(d) { w <- rep(1, length(d)); w/sum(w) }
  } else {
    function(d) { w <- exp(-(d^2)/(2*sigma^2)); w/sum(w) }
  }
  if (is.null(out_suffix)) {
    out_suffix <- if (method == "mean") paste0("_knn", k) else paste0("_wknn", k)
  }
  
  # 平滑 meta 列
  if (!is.null(meta_cols) && length(meta_cols)) {
    for (nm in meta_cols) {
      v <- object@meta.data[[nm]]
      if (is.null(v)) { message("[skip] meta 列不存在：", nm); next }
      v <- as.numeric(v)
      if (all(!is.finite(v))) { message("[skip] meta 列非数值：", nm); next }
      sm <- v
      for (i in seq_along(v)) {
        idx <- c(i, nn$nn.index[i,]); d <- c(0, nn$nn.dist[i,]); w <- get_weights(d)
        sm[i] <- sum(v[idx] * w, na.rm=TRUE)
      }
      object[[paste0(nm, out_suffix)]] <- sm
    }
  }
  
  # 平滑基因（从“已选好的 assay$data”取）
  if (!is.null(genes) && length(genes)) {
    genes2 <- intersect(genes, rownames(M))
    if (length(genes2) < length(genes)) {
      miss <- setdiff(genes, genes2)
      if (length(miss)) message("[WARN] 基因不存在于 ", choice$assay, "/", choice$slot, "：", paste(miss, collapse=", "))
    }
    for (g in genes2) {
      v <- as.numeric(M[g, ])
      sm <- v
      for (i in seq_along(v)) {
        idx <- c(i, nn$nn.index[i,]); d <- c(0, nn$nn.dist[i,]); w <- get_weights(d)
        sm[i] <- sum(v[idx] * w, na.rm=TRUE)
      }
      object[[paste0(g, out_suffix)]] <- sm
    }
  }
  
  attr(object, "knn_smooth_source") <- list(
    assay=choice$assay, slot=choice$slot,
    normalized = choice$normalize_if_needed,
    method=method, k=k, sigma=sigma
  )
  object
}



get_markers_by_celltype <- function(sobj, assay_for_DE = "SCT",
                                    downsample_per_cluster = 3000L) {
  Idents(sobj) <- factor(sobj$banksy_celltype)
  
  ## 保证 data 层存在
  if (ncol(suppressWarnings(GetAssayData(sobj, assay_for_DE, slot="data"))) == 0) {
    sobj <- NormalizeData(sobj, assay = assay_for_DE)
  }
  
  set.seed(1)
  cells_use <- unlist(tapply(Cells(sobj), Idents(sobj), function(v){
    if (length(v) > downsample_per_cluster) sample(v, downsample_per_cluster) else v
  }))
  sobj_ds <- subset(sobj, cells = cells_use)
  
  mat <- as.matrix(GetAssayData(sobj_ds, assay = assay_for_DE, slot = "data"))
  y   <- as.character(Idents(sobj_ds))
  
  stopifnot(ncol(mat) == length(y))  ## 关键断言，防止再出错
  
  if (requireNamespace("presto", quietly = TRUE)) {
    res <- presto::wilcoxauc(X = mat, y = y, seurat = FALSE)
    colnames(res)[colnames(res)=="group"] <- "cluster"
    res <- res |>
      dplyr::arrange(cluster, dplyr::desc(auc), dplyr::desc(logFC))
  } else {
    res <- FindAllMarkers(
      sobj_ds, assay = assay_for_DE,
      only.pos = TRUE, test.use = "wilcox",
      logfc.threshold = 0.25, min.pct = 0.1
    ) |>
      dplyr::arrange(cluster, dplyr::desc(avg_log2FC))
  }
  res
}

## 调用





# 替代方案：更安全的版本，不修改原始对象
PlotFeatureViolin <- function(
    object,
    features,
    group_by,
    assay = NULL,
    layer = "data",
    sort_by_value = TRUE,
    add_stats = FALSE
) {
  
  # --- 1. 参数和数据检查 ---
  if (!group_by %in% colnames(object@meta.data)) {
    stop("在Seurat对象的meta.data中找不到指定的分组列: '", group_by, "'")
  }
  
  # 检查指定的assay是否存在
  if (!is.null(assay) && !assay %in% Assays(object)) {
    stop("指定的assay '", assay, "' 不存在于Seurat对象中。")
  }
  
  # 如果没有指定assay，且feature不是meta.data列，则使用默认assay
  if (is.null(assay)) {
    is_meta <- all(features %in% colnames(object@meta.data))
    if (!is_meta) {
      assay <- DefaultAssay(object)
    }
  }
  
  # --- 2. 安全方法：创建ggplot数据并手动排序 ---
  if (sort_by_value && length(features) == 1) {
    # 提取数据
    plot_data <- FetchData(object, vars = c(group_by, features), assay = assay, layer = layer)
    
    # 计算中位数并排序
    median_scores <- plot_data %>%
      group_by(.data[[group_by]]) %>%
      summarise(median_val = median(.data[[features]], na.rm = TRUE)) %>%
      arrange(desc(median_val))
    
    # 重新排列数据中的因子水平
    plot_data[[group_by]] <- factor(plot_data[[group_by]], levels = median_scores[[group_by]])
    
    # 手动创建小提琴图
    p <- ggplot(plot_data, aes_string(x = group_by, y = features, fill = group_by)) +
      geom_violin(scale = "width", alpha = 0.7) +
      geom_boxplot(width = 0.1, alpha = 0.8, outlier.shape = NA) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
      ) +
      labs(title = features, y = "Expression")
    
    # 添加统计检验
    if (add_stats && requireNamespace("ggpubr", quietly = TRUE)) {
      p <- p + ggpubr::stat_compare_means(
        label.y.npc = 0.9, 
        method = "wilcox.test", 
        ref.group = ".all.",
        label = "p.signif"
      )
    }
    
  } else {
    # 使用原始的VlnPlot（不排序）
    p <- VlnPlot(
      object,
      features = features,
      group.by = group_by,
      assay = assay,
      layer = layer,
      pt.size = 0
    )
    
    # 添加统计检验
    if (add_stats && requireNamespace("ggpubr", quietly = TRUE)) {
      p <- p + ggpubr::stat_compare_means(
        label.y.npc = 0.9, 
        method = "wilcox.test", 
        ref.group = ".all.",
        label = "p.signif"
      )
    }
    
    # 美化
    p <- p + NoLegend()
    p <- p & theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  }
  
  return(p)
}



# =========================
# 极简 & 稳定版 GSEA（重构）
# =========================
# 依赖：dplyr, ggplot2, msigdbr, fgsea
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# 将数据框中的 list 列统一折叠为字符串，避免写 CSV 报错
collapse_list_cols <- function(df, sep = ";"){
  is_list <- vapply(df, is.list, logical(1))
  if (any(is_list)){
    for (nm in names(df)[is_list]){
      df[[nm]] <- vapply(df[[nm]], function(x){
        if (is.null(x)) return(NA_character_)
        if (length(x) == 0) return(NA_character_)
        paste(as.character(x), collapse = sep)
      }, character(1))
    }
  }
  df
}

# 单组 GSEA（无并行，避免 serialize 警告）
run_fgsea_one <- function(stats_named, pathways, nperm = 2000,
                          minSize = 10, maxSize = 2000, scoreType = "pos"){
  # 清理：去 NA / 去 0 / 名称去重
  ok <- !is.na(stats_named) & is.finite(stats_named)
  stats_named <- stats_named[ok]
  stats_named <- stats_named[stats_named != 0]
  # fgsea（simple 单线程）
  fgsea::fgseaSimple(
    pathways = pathways,
    stats    = stats_named,
    nperm    = nperm,
    minSize  = minSize,
    maxSize  = maxSize,
    scoreType = scoreType,
    nproc    = 1
  )
}

# 主函数：按给定分组字段做 GSEA；导出 CSV（全通路）与 PDF（每组TopN，可选全通路）
RunGSEA_Reboot <- function(
    markers_df,                   # 必含：gene, avg_log2FC, group_by列
    group_by          = "celltype_manual",
    groups_to_analyze = NULL,     # e.g. c("Tumor_Epi_Diff", "CAF_myo")
    msigdb_category   = "H",
    msigdb_subcategory = NULL,
    species           = "Homo sapiens",
    nperm             = 2000,
    minSize           = 10,
    maxSize           = 2000,
    scoreType         = c("auto","pos","std"),
    out_dir           = "GSEA_reboot",
    plot_topN         = 30,
    make_all_pdf      = FALSE     # TRUE 时也画“全通路长图”
){
  req <- c("gene","avg_log2FC", group_by)
  stopifnot(all(req %in% colnames(markers_df)))
  
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("需要安装 msigdbr：install.packages('msigdbr')")
  }
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("需要安装 fgsea：install.packages('fgsea')")
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  df <- markers_df %>%
    dplyr::select(all_of(req)) %>%
    dplyr::rename(.group = !!group_by, .gene = gene, .fc = avg_log2FC)
  
  if (!is.null(groups_to_analyze)){
    df <- df %>% dplyr::filter(.group %in% groups_to_analyze)
  }
  groups <- unique(as.character(df$.group))
  if (!length(groups)) stop("筛选后无可分析的组。")
  
  # 载入基因集一次
  msig <- msigdbr::msigdbr(species = species,
                           category = msigdb_category,
                           subcategory = msigdb_subcategory)
  pathways <- split(msig$gene_symbol, msig$gs_name)
  
  # 逐组运行
  all_list <- vector("list", length(groups))
  names(all_list) <- groups
  
  for (g in groups){
    sub <- df %>% dplyr::filter(.group == g)
    # 去重：同一基因只保留最大 logFC（或用 mean/median 均可）
    sub <- sub %>%
      dplyr::group_by(.gene) %>%
      dplyr::slice_max(order_by = .fc, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup()
    
    stats <- sub$.fc
    names(stats) <- sub$.gene
    
    # 自动判断得分类型：若存在负值，用 "std"，否则 "pos"
    st <- match.arg(scoreType)
    if (st == "auto"){
      if (any(stats < 0, na.rm = TRUE)) st <- "std" else st <- "pos"
    }
    
    fg <- run_fgsea_one(stats, pathways,
                        nperm = nperm, minSize = minSize, maxSize = maxSize, scoreType = st)
    if (nrow(fg)){
      fg$group <- g
      all_list[[g]] <- as.data.frame(fg)
    }
  }
  
  enr <- dplyr::bind_rows(all_list, .id = NULL)
  if (!nrow(enr)) {
    warning("没有任何富集结果。")
    return(invisible(NULL))
  }
  
  # 写 CSV：先折叠 list 列
  enr_out <- collapse_list_cols(enr)
  csv_path <- file.path(out_dir, "GSEA_results_all_groups.csv")
  utils::write.csv(enr_out, csv_path, row.names = FALSE)
  message("CSV 已导出：", csv_path)
  
  # 简洁画图（每组一张 TopN 点图；可选再画“全通路长图”）
  plot_one_group <- function(df_g, title, outfile_top, outfile_all = NULL, topN = 30){
    df_g <- df_g %>% arrange(padj, desc(NES))
    df_g$pathway_clean <- sub("^HALLMARK_", "", df_g$pathway)
    
    # TopN
    df_top <- head(df_g, topN)
    h <- max(4, 0.22 * nrow(df_top) + 2)
    p_top <- ggplot(df_top, aes(x = NES, y = reorder(pathway_clean, NES))) +
      geom_point(aes(size = -log10(padj), color = NES)) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red") +
      scale_size(range = c(1.5, 5)) +
      theme_bw(base_size = 9) +
      labs(title = paste0(title, " (Top ", nrow(df_top), ")"),
           x = "NES", y = NULL)
    ggsave(outfile_top, p_top, width = 7, height = h, dpi = 300, bg = "white")
    
    # 全通路（可选）
    if (!is.null(outfile_all)){
      h2 <- max(4, 0.22 * nrow(df_g) + 2)
      p_all <- ggplot(df_g, aes(x = NES, y = reorder(pathway_clean, NES))) +
        geom_point(aes(size = -log10(padj), color = NES)) +
        scale_color_gradient2(low = "blue", mid = "white", high = "red") +
        scale_size(range = c(1.5, 5)) +
        theme_bw(base_size = 9) +
        labs(title = title, x = "NES", y = NULL)
      ggsave(outfile_all, p_all, width = 7, height = h2, dpi = 300, bg = "white")
    }
  }
  
  for (g in groups){
    df_g <- enr[enr$group == g, , drop = FALSE]
    if (!nrow(df_g)) next
    base <- gsub("[^A-Za-z0-9._-]+", "_", g)
    top_pdf <- file.path(out_dir, paste0("GSEA_", base, "_Top", plot_topN, ".pdf"))
    if (isTRUE(make_all_pdf)){
      all_pdf <- file.path(out_dir, paste0("GSEA_", base, "_ALL.pdf"))
    } else {
      all_pdf <- NULL
    }
    plot_one_group(df_g, title = paste0("MSigDB ", msigdb_category, " - ", g),
                   outfile_top = top_pdf, outfile_all = all_pdf, topN = plot_topN)
  }
  
  invisible(enr_out)
}




PlotSpatialHighlight <- function(
    bcs_df,
    barcodes_to_highlight,
    ptsize = 2.5,
    color_highlight = "red",
    color_background = "grey80"
) {
  
  # 1. 在数据框中创建一个用于标记的列
  bcs_df$highlight_group <- ifelse(bcs_df$barcode %in% barcodes_to_highlight, 
                                   "Highlight", 
                                   "Background")
  
  # 2. 创建一个命名的颜色向量，用于手动指定颜色
  plot_colors <- c(
    "Highlight" = color_highlight,
    "Background" = color_background
  )
  
  # 3. 使用ggplot2和geom_scattermore进行分层绘图
  p <- ggplot(bcs_df, aes(x = imagecol_scaled, y = -imagerow_scaled, color = highlight_group)) +
    
    # 使用scale_color_manual来确保即使只有一个组别，颜色也是正确的
    scale_color_manual(values = plot_colors, name = "Region") +
    
    # 直接用一个geom_scattermore，ggplot会自动按颜色分组绘制
    geom_scattermore(pointsize = ptsize) +
    
    # 美化
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none") # 通常高亮图不需要图例
  
  return(p)
}

#FindSpatialNichesFromScores
library(scattermore)
FindSpatialNichesFromScores <- function(
    bcs_df,
    score_col,
    PATH,
    N_hotspots = 3,
    top_quantile = 0.90,
    niche_size_microns = 200
) {
  
  # --- 1. 参数和包检查 ---
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("请先安装 'MASS' 包: install.packages('MASS')")
  }
  if (!score_col %in% colnames(bcs_df)) {
    stop("在数据框中找不到指定的分数列: '", score_col, "'")
  }
  # 确保所有必需的辅助函数都已加载
  required_funcs <- c("findPeaks3D", "GetSlice")
  if (!all(sapply(required_funcs, exists, where = .GlobalEnv, mode = "function"))) {
    stop("一个或多个必需的辅助函数 (findPeaks3D, GetSlice) 未找到。",
         "请确保完整的 Auxiliary Functions.R 已被 source。")
  }
  
  cat("--- 开始从 '", score_col, "' 分数中识别空间Niche ---\n")
  
  # --- 2. 筛选高分点用于密度估计 ---
  # 移除NA值以确保quantile能正常工作
  valid_scores <- bcs_df[[score_col]][!is.na(bcs_df[[score_col]])]
  score_threshold <- quantile(valid_scores, probs = top_quantile, na.rm = TRUE)
  
  high_score_df <- bcs_df %>% filter(.data[[score_col]] > score_threshold)
  
  if (nrow(high_score_df) < 20) { # 增加一个更稳健的下限
    stop("高分点数量过少 (少于20个)，无法进行可靠的密度估计。请尝试降低 'top_quantile'。")
  }
  
  cat("  - 使用", nrow(high_score_df), "个高分点 (top", (1-top_quantile)*100, "%) 进行密度估计。\n")
  
  # --- 3. 计算二维核密度并找到峰值 ---
  # 使用 scaled 坐标进行密度计算，以获得更好的尺度和一致性
  Kernel <- MASS::kde2d(high_score_df$imagecol_scaled, high_score_df$imagerow_scaled, n = 200)
  
  Peaks <- findPeaks3D(Kernel$z, N_hotspots)
  
  # 将峰值索引转换回 scaled 坐标
  Peaks_df <- data.frame(
    X_scaled = Kernel$x[sapply(Peaks, function(p) p$x)],
    Y_scaled = Kernel$y[sapply(Peaks, function(p) p$y)]
  )
  cat("  - 找到了", nrow(Peaks_df), "个潜在的热点中心。\n")
  
  # --- 4. 以每个热点为中心，提取Niche区域 ---
  niche_regions_barcodes <- vector("list", length = nrow(Peaks_df))
  names(niche_regions_barcodes) <- paste0("Niche_", 1:nrow(Peaks_df))
  
  for (i in 1:nrow(Peaks_df)) {
    peak_xy_scaled <- c(Peaks_df$X_scaled[i], Peaks_df$Y_scaled[i])
    
    # 在scaled坐标系下，计算所有点到理论峰值的欧氏距离
    distances_to_peak <- sqrt((bcs_df$imagecol_scaled - peak_xy_scaled[1])^2 + 
                                (bcs_df$imagerow_scaled - peak_xy_scaled[2])^2)
    
    # 找到距离最近的那个真实存在的条形码作为“种子点”
    center_barcode <- bcs_df$barcode[which.min(distances_to_peak)]
    
    cat("  - Niche", i, "的中心种子点是:", center_barcode, "\n")
    
    # 使用 GetSlice 提取该种子点周围的区域
    # 确保传递给GetSlice的数据框有原始像素坐标
    niche_regions_barcodes[[i]] <- GetSlice(
      Spot = center_barcode,
      SizeMicrons = niche_size_microns,
      BarcodeDF = bcs_df, # GetSlice内部依赖'imagerow'和'imagecol'
      PATH = PATH
    )
  }
  
  cat("--- Niche识别完成！---\n")
  return(niche_regions_barcodes)
}


# Create paths and read different spatial files found in output directory
get_spatial_files<-function (PATH, type) 
{
  if (type == "tissue_lowres_image") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/tissue_lowres_image.png", sep = ""))
  }
  if (type == "tissue_hires_image") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/tissue_hires_image.png", sep = ""))
  }
  if (type == "tissue_positions_list") {
    x <- read.csv(paste(PATH, "/spatial/tissue_positions_list.csv", sep = ""), col.names = c("barcode", "tissue", "row",  "col", "imagerow", "imagecol"), header = F)
  }
  if (type == "aligned_fiducials") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/aligned_fiducials.jpg", sep = ""))
  }
  if (type == "detected_tissue_image") {
    x <- readbitmap::read.bitmap(paste(PATH, "/spatial/detected_tissue_image.jpg", sep = ""))
  }
  if (type == "scales") {
    path_scales <- paste(PATH, "/spatial/scalefactors_json.json", sep = "")
    x <- rjson::fromJSON(file = path_scales)
  }
  return(x)
}

# Define color palettes used throughout the manuscript
ColorPalette<-function()
{
  Colors<-c(paletteer::paletteer_d("ggsci::default_igv")[1:39],"black","azure4")
  names(Colors)<-c('Tumor III','Plasma','Macrophage','CD4 T cell','CAF','vSM','Mature B','Endothelial','Tumor I','CD8 T cell',
                   'Enterocyte','Neutrophil','Proliferating Immune II','Pericytes','Smooth Muscle','Myofibroblast',
                   'Tumor II','Fibroblast','Goblet','Lymphatic Endothelial','Tumor V','Proliferating Macrophages','SM Stress Response',
                   'NK','cDC I','Tumor IV','Proliferating Fibroblast','Epithelial','Tuft','Mast','Unknown III (SM)',
                   'Adipocyte','mRegDC','Enteric Glial','pDC','Vascular Fibroblast','Neuroendocrine','Memory B','Unknwon I (Immune)',"Undetermined","cell")
  
  return(Colors)
  
}

# Read deconvolution results. Wrapped in a function as large RCTD objects slow down the R session
readRCTD<-function(PATH)
{
  Object<-readRDS(PATH)
  Info <- Object@results$results_df
  weights <- get_doublet_weights_modified(Info,Object@results$weights_doublet, Object@cell_type_info$info[[2]])
  Results <- list(DF = Info, Weights = t(weights))
  
  return(Results)
}

# ================================================================
# AuxFunctionscellbin.R — RCTD Pass-1（doublet）工具集
#   - 兼容 Seurat v5 抽取 Cellbin counts/coords
#   - 稳健稀疏矩阵转换：先 CsparseMatrix 再 dgCMatrix（+整数化）
#   - 读取 SCE 参考，做基因交集、全零过滤
#   - RCTD 运行（doublet）
#   - 写回 Seurat + 导出 Loupe 分类 + 三种 PlotSpatialDistribution
# ================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SummarizedExperiment)  # assay()
  library(SingleCellExperiment)
  library(spacexr)
  library(dplyr)
})

# -------------------- 稀疏矩阵保障：Csparse -> dgCMatrix -> integer --------------------
ensure_csparse_integer <- function(mat) {
  # 先尽量变为 Csparse（虚类，覆盖 dgC/lsC/lgC 等），再收敛为 dgCMatrix
  if (!inherits(mat, "CsparseMatrix")) {
    mat <- tryCatch(as(mat, "CsparseMatrix"),
                    error = function(e) as(as.matrix(mat), "CsparseMatrix"))
  }
  if (!inherits(mat, "dgCMatrix")) {
    mat <- tryCatch(as(mat, "dgCMatrix"),
                    error = function(e) {
                      # 若是逻辑/整型/长整型等，强转为双精度
                      as(as(mat, "dgCMatrix"), "dgCMatrix")
                    })
  }
  # RCTD 需要整数 UMI
  if (!is.integer(mat@x)) mat@x <- round(mat@x)
  mat
}

# -------------------- 从 Seurat Cellbin 抽取 counts + coords --------------------
get_cellbin_counts_coords <- function(obj, assay = "Spatial",
                                      umi_min = 200, use_scaled = TRUE) {
  # counts：Seurat v5 优先 layer="counts"，否则回退 slot="counts"
  counts <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = "counts")
  )
  counts <- ensure_csparse_integer(counts)
  
  # coords：优先 meta.data，其次 image@coordinates
  md <- obj@meta.data
  cand_pairs <- if (use_scaled) {
    list(c("imagecol_scaled","imagerow_scaled"),
         c("imagecol","imagerow"),
         c("x","y"), c("X","Y"))
  } else {
    list(c("imagecol","imagerow"),
         c("imagecol_scaled","imagerow_scaled"),
         c("x","y"), c("X","Y"))
  }
  pick_cols <- NULL
  for (p in cand_pairs) if (all(p %in% colnames(md))) { pick_cols <- p; break }
  if (!is.null(pick_cols)) {
    coords <- md[, pick_cols, drop = FALSE]
    colnames(coords) <- c("x","y")
    rownames(coords) <- rownames(md)
  } else {
    if (length(Seurat::Images(obj)) < 1) stop("No image in Seurat object.")
    img <- obj@images[[ Seurat::Images(obj)[1] ]]
    C <- img@coordinates
    cand_pairs_img <- list(c("pxl_col_in_fullres","pxl_row_in_fullres"),
                           c("imagecol","imagerow"),
                           c("x","y"), c("X","Y"))
    for (p in cand_pairs_img) if (all(p %in% colnames(C))) { pick_cols <- p; break }
    if (is.null(pick_cols)) stop("Cannot find coordinate columns.")
    coords <- C[, pick_cols, drop = FALSE]
    colnames(coords) <- c("x","y")
  }
  
  # 对齐细胞
  common <- intersect(colnames(counts), rownames(coords))
  if (length(common) < 100) stop("Too few common cells between counts and coords.")
  counts <- counts[, common, drop = FALSE]
  coords <- coords[common, , drop = FALSE]
  
  # UMI 过滤
  nUMI <- Matrix::colSums(counts)
  keep <- nUMI >= umi_min
  counts <- counts[, keep, drop = FALSE]
  coords <- coords[keep, , drop = FALSE]
  
  # RCTD 偏好 data.frame 坐标；行名与 counts 列名一致
  coords <- data.frame(
    x = as.numeric(coords[,1]),
    y = as.numeric(coords[,2]),
    row.names = rownames(coords),
    check.names = FALSE
  )
  
  list(counts = counts, coords = coords)
}

# -------------------- 结果抽取（兼容不同 spacexr 版本） --------------------
extract_rctd_meta <- function(rctd_obj){
  res <- rctd_obj@results
  if (!is.null(res$results_df) && is.data.frame(res$results_df) && nrow(res$results_df) > 0) {
    df <- as.data.frame(res$results_df)
  } else if (is.data.frame(res) && nrow(res) > 0) {
    df <- as.data.frame(res)
  } else {
    stop("RCTD results appear empty; run.RCTD 可能未跑完或失败。")
  }
  rn <- rownames(df); if (is.null(rn) || !length(rn)) stop("results_df 缺少行名（cell_id）。")
  
  spot_class <- if ("spot_class" %in% names(df)) as.character(df$spot_class)
  else if ("class" %in% names(df)) as.character(df$class)
  else if (!is.null(res$spot_class) && length(res$spot_class) == nrow(df)) as.character(res$spot_class)
  else rep(NA_character_, nrow(df))
  
  first_type  <- if ("first_type"  %in% names(df)) as.character(df$first_type)
  else if ("cell_type" %in% names(df)) as.character(df$cell_type)
  else rep(NA_character_, nrow(df))
  second_type <- if ("second_type" %in% names(df)) as.character(df$second_type)
  else rep(NA_character_, nrow(df))
  
  out <- data.frame(
    cell_id     = rn,
    spot_class  = spot_class,
    first_type  = first_type,
    second_type = second_type,
    stringsAsFactors = FALSE
  )
  rownames(out) <- rn
  out
}

extract_rctd_weights <- function(rctd_obj){
  res <- rctd_obj@results
  cand <- list(res$weights, res$norm_weights, res$deconvolution_weights)
  if (!is.null(res$results_df) && "weights" %in% names(res$results_df)) {
    cand <- c(list(res$results_df$weights), cand)
  }
  for (w in cand) {
    if (!is.null(w)) {
      W <- tryCatch(as.matrix(w), error = function(e) NULL)
      if (!is.null(W) && nrow(W) > 0 && ncol(W) > 0) return(W)
    }
  }
  NULL
}

# -------------------- 写回 Seurat + 导出 Loupe + 画图 --------------------
write_back_rctd_to_seurat <- function(object, meta_df,
                                      result_path,
                                      do_export_loupe = TRUE,
                                      do_plots = TRUE,
                                      sample_id = NULL) {
  stopifnot(is.data.frame(meta_df), "cell_id" %in% colnames(meta_df))
  # 初始化列
  object$RCTD_spot_class <- NA_character_
  object$RCTD_first      <- NA_character_
  object$RCTD_second     <- NA_character_
  # 对齐并写回
  idx <- match(rownames(object@meta.data), meta_df$cell_id)
  keep <- !is.na(idx)
  object$RCTD_spot_class[keep] <- meta_df$spot_class[idx[keep]]
  object$RCTD_first[keep]      <- meta_df$first_type[idx[keep]]
  object$RCTD_second[keep]     <- meta_df$second_type[idx[keep]]
  
  # Loupe 类别导出
  if (do_export_loupe && exists("export_loupe_categories")) {
    export_loupe_categories(
      object,
      meta_col = "RCTD_first",
      out_dir  = result_path,
      filename = NULL,
      na_label = "Unlabeled",
      gzip = FALSE,
      overwrite = TRUE,
      verbose = TRUE
    )
  }
  
  # 空间图（若函数可用）
  if (do_plots && exists("PlotSpatialDistribution")) {
    # 三个图：class / first / second
    suppressWarnings({
      try(PlotSpatialDistribution(object, meta_col = "RCTD_spot_class",
                                  point_size = 0.5,
                                  title = paste0(sample_id %||% "", " RCTD spot_class"),
                                  out_dir = file.path(result_path, "plots")), silent = TRUE)
      try(PlotSpatialDistribution(object, meta_col = "RCTD_first",
                                  point_size = 0.5,
                                  title = paste0(sample_id %||% "", " RCTD first"),
                                  out_dir = file.path(result_path, "plots")), silent = TRUE)
      try(PlotSpatialDistribution(object, meta_col = "RCTD_second",
                                  point_size = 0.5,
                                  title = paste0(sample_id %||% "", " RCTD second"),
                                  out_dir = file.path(result_path, "plots")), silent = TRUE)
    })
  }
  object
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -------------------- 主流程：RunCellbinRCTD_Pass1 --------------------
# 输入：
#   object         : Seurat（Cellbin）
#   ref_rds        : 参考 SCE rds（例：SCP1089 untreated）
#   result_path    : 输出目录（会建子目录 "RCTD_pass1"）
#   umi_min        : UMI 阈值（Cellbin 80–300；8um spot ~100）
#   cell_min_inst  : RCTD CELL_MIN_INSTANCE（参考每大类最少细胞数，20–50）
RunCellbinRCTD_Pass1 <- function(object,
                                 ref_rds,
                                 result_path,
                                 umi_min = 200,
                                 cell_min_inst = 20,
                                 sample_id = NULL) {
  out_dir <- file.path(result_path, "RCTD_pass1")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1) 抽取 Cellbin 数据
  cb <- get_cellbin_counts_coords(object, assay = "Spatial", umi_min = umi_min)
  counts_cb <- cb$counts
  coords_cb <- cb$coords
  rm(cb)
  message("Cellbin after QC: genes=", nrow(counts_cb), " cells=", ncol(counts_cb))
  
  # 2) 加载参考
  sce_ref <- readRDS(ref_rds)
  message("Loaded reference: ", ref_rds)
  
  # 3) 基因交集 + 清洗
  genes_common <- intersect(rownames(counts_cb), rownames(sce_ref))
  if (length(genes_common) < 500) stop("Too few intersecting genes between Cellbin and reference.")
  
  ref_counts <- SummarizedExperiment::assay(sce_ref, "counts")
  ref_counts <- ensure_csparse_integer(ref_counts)
  ref_counts <- ref_counts[genes_common, , drop = FALSE]
  
  counts_use <- counts_cb[genes_common, , drop = FALSE]
  counts_use <- ensure_csparse_integer(counts_use)
  
  keep_g <- Matrix::rowSums(ref_counts) > 0 & Matrix::rowSums(counts_use) > 0
  ref_counts <- ref_counts[keep_g, , drop = FALSE]
  counts_use <- counts_use[keep_g, , drop = FALSE]
  
  ref_nonzero <- Matrix::colSums(ref_counts) > 0
  ref_counts  <- ref_counts[, ref_nonzero, drop = FALSE]
  
  ref_types <- as.character(sce_ref$cell_type_major)[ref_nonzero]
  ok <- !is.na(ref_types)
  ref_counts <- ref_counts[, ok, drop = FALSE]
  ref_types  <- factor(ref_types[ok])
  names(ref_types) <- colnames(ref_counts)
  
  stopifnot(ncol(ref_counts) == length(ref_types))
  if (anyDuplicated(names(ref_types)) > 0) stop("Duplicated reference cell barcodes.")
  
  # counts 与 coords 一一对应
  stopifnot(identical(colnames(counts_use), rownames(coords_cb)))
  
  # 4) 构建并运行 RCTD（doublet）
  ref_obj <- spacexr::Reference(ref_counts, ref_types)
  sp_obj  <- spacexr::SpatialRNA(
    counts = counts_use,
    coords = coords_cb,                       # data.frame (x,y)
    nUMI   = Matrix::colSums(counts_use)
  )
  rctd1 <- spacexr::create.RCTD(sp_obj, ref_obj, CELL_MIN_INSTANCE = cell_min_inst)
  rctd1 <- spacexr::run.RCTD(rctd1, doublet_mode = "doublet")
  
  # 5) 导出结果
  meta_df <- extract_rctd_meta(rctd1)
  write.csv(meta_df, file.path(out_dir, "RCTD_pass1_meta.csv"), row.names = FALSE)
  
  W <- extract_rctd_weights(rctd1)
  if (!is.null(W)) saveRDS(W, file.path(out_dir, "RCTD_pass1_weights.rds"))
  
  # 6) 写回 Seurat + 导出 Loupe + 空间图
  object <- write_back_rctd_to_seurat(object, meta_df,
                                      result_path = out_dir,
                                      do_export_loupe = TRUE,
                                      do_plots = TRUE,
                                      sample_id = sample_id)
  
  # 7) 保存 Seurat
  out_rds <- file.path(out_dir, paste0(sample_id %||% "sample", "_Seurat_with_RCTD_pass1.rds"))
  saveRDS(object, out_rds)
  message("Saved: ", out_rds)
  
  invisible(list(object = object,
                 meta = meta_df,
                 weights = W,
                 rctd = rctd1,
                 out_dir = out_dir))
}

# -------------------- 多样本便捷跑法（各自独立目录） --------------------
# cfg_list: list(
#   list(SampleID="HDA1", OutsPath="I:/HDA1_count/outs/", ResultPath="I:/HDA1_count/outs/Analysis_Results_Cellbin_noqc/"),
#   ...
# )
RunCellbinRCTD_Pass1_ForSamples <- function(cfg_list,
                                            ref_rds,
                                            umi_min = 200,
                                            cell_min_inst = 20) {
  out <- list()
  for (cfg in cfg_list) {
    sample_id <- cfg$SampleID
    message("=== Running RCTD Pass-1 for: ", sample_id, " ===")
    # 用户已有的函数：GenerateCellbinSeurat()
    obj <- GenerateCellbinSeurat(path_outs = cfg$OutsPath, sample_id = sample_id)
    res <- RunCellbinRCTD_Pass1(object = obj,
                                ref_rds = ref_rds,
                                result_path = cfg$ResultPath,
                                umi_min = umi_min,
                                cell_min_inst = cell_min_inst,
                                sample_id = sample_id)
    out[[sample_id]] <- res
  }
  invisible(out)
}




# From spaceXR outputs create weights matrix
get_doublet_weights_modified <- function(ResultDF,Weights,CellTypes) {
  
  barcodes <- rownames(ResultDF)
  my_beta <- matrix(0, nrow = length(barcodes), ncol = length(CellTypes))
  rownames(my_beta) <- barcodes
  colnames(my_beta) <- CellTypes
  
  indexRow_Certain<-which(ResultDF$spot_class %in% c('singlet', 'doublet_certain'))
  indexCol_Certain<-match(ResultDF[indexRow_Certain,'first_type'],colnames(my_beta))
  my_beta[cbind(indexRow_Certain,indexCol_Certain)] <- Weights[indexRow_Certain,1]
  
  indexRow_Doublet<-which(ResultDF$spot_class == "doublet_certain")
  indexCol_Doublet<-match(ResultDF[indexRow_Doublet,'second_type'],colnames(my_beta))
  my_beta[cbind(indexCol_Doublet)] <- Weights[indexRow_Doublet,2]
  
  return(my_beta)
  
  
}

# Add deconvolution results to data.frame created with GenerateSampleData
AddDeconvolutionInfo<-function(BCS,Results,AddWeights=FALSE)
{
  ResultsDF<-Results$DF
  
  index<- match(rownames(ResultsDF),BCS$barcode)
  
  BCS$DeconvolutionClass<-NA
  BCS$DeconvolutionClass[index]<-as.vector(ResultsDF$spot_class)
  
  BCS$DeconvolutionLabel1<-NA
  BCS$DeconvolutionLabel1[index]<-ResultsDF$first_type
  
  BCS$DeconvolutionLabel2<-NA
  BCS$DeconvolutionLabel2[index]<-ResultsDF$scond_type
  
  if(AddWeights)
  {
    Weights<-Results$Weights
    
    index<- match(colnames(Weights),BCS$barcode)
    
    Names<-gsub(" ","",rownames(Weights))
    
    for(jj in 1:nrow(Weights))
    {
      BCS[,Names[jj]]<-NA
      BCS[index,Names[jj]]<-Weights[jj,]
    }
  }
  
  
  return(BCS)
  
}

# Add expression of genes to a data.frame generated by GenerateSampleData from a Seurat Object
AddExpression<-function(Barcodes,Seurat,Genes)
{
  
  Exp<-FetchData(Seurat,Genes)
  
  for(Gx in Genes)
  {
    Barcodes[,Gx]<-NA
    Barcodes[match(rownames(Exp),Barcodes$barcode),Gx]<-Exp[,Gx]
  }
  
  return(Barcodes)
  
}

# Create a square of a given size (in microns) whose center is a given barcode
GetSquare<-function(Spot,SizeMicrons,BarcodeDF,binsize=8)
{
  Xcenter<-BarcodeDF$col[match(Spot,BarcodeDF$barcode)]
  Ycenter<-BarcodeDF$row[match(Spot,BarcodeDF$barcode)]
  
  AddFactor<-round(SizeMicrons/(2*binsize))
  
  Xmin<-Xcenter-AddFactor
  Xmax<-Xcenter+AddFactor
  
  Ymin<-Ycenter-AddFactor
  Ymax<-Ycenter+AddFactor
  
  SquareSection<-BarcodeDF %>% filter(col >= Xmin & col <= Xmax & row >= Ymin & row <= Ymax) %>% pull(barcode)
  
  return(SquareSection)
  
}
TransferBanksyAnno <- function(target, source,
                               col    = "banksy_domain_anno",  # 源列名
                               to_col = "banksy_domain_anno",  # 目标写入列名
                               trim   = TRUE, keep_levels = TRUE) {
  stopifnot(inherits(target, "Seurat"), inherits(source, "Seurat"))
  if (!col %in% colnames(source@meta.data))
    stop("Column '", col, "' not found in source@meta.data")
  
  # 取源向量并按细胞名命名
  v <- source@meta.data[[col]]
  names(v) <- rownames(source@meta.data)
  
  # 对齐细胞名
  tgt_cells <- Cells(target)
  src_cells <- names(v)
  common    <- intersect(tgt_cells, src_cells)
  
  # 写入（未匹配的置 NA）
  out <- rep(NA_character_, length(tgt_cells)); names(out) <- tgt_cells
  out[common] <- as.character(v[common])
  if (trim) out <- trimws(out)
  
  # 设置因子 levels 与源一致（如果源是因子且 keep_levels=TRUE）
  if (is.factor(v) && keep_levels) {
    target[[to_col]] <- factor(out, levels = levels(v))
  } else {
    target[[to_col]] <- out
  }
  
  message(sprintf("Transferred '%s' -> '%s': matched %d / %d cells; %d NA.",
                  col, to_col, length(common), length(tgt_cells),
                  sum(is.na(out))))
  return(target)
}

# 使用：


# Plot gene expression from a data.frame 
PlotExpression<-function(barcodes,Gene,ptsize=2,shape="circle")
{
  barcodes$Expression<-barcodes %>% pull(Gene)
  
  if(shape=="circle")
  {
    Plot<-barcodes %>%
      ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=Expression)) +
      geom_scattermore(pointsize = ptsize,pixels = rep(2000,2))+
      coord_cartesian(expand = FALSE) +
      xlab("") +
      ylab("") +
      theme_set(theme_bw(base_size = 10))+
      theme_minimal() +
      theme(axis.text = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())+
      scale_color_gradient(low="lightgray",high = "red")+
      labs(color=paste0(Gene))
    
  }else if(shape=="square")
  {
    Plot<-barcodes %>%
      ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,fill=Expression)) +
      geom_point(shape=22,size=ptsize,color=alpha("black",0),stroke=0.25)+
      coord_cartesian(expand = FALSE) +
      xlab("") +
      ylab("") +
      theme_set(theme_bw(base_size = 10))+
      theme_minimal() +
      theme(axis.text = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank())+
      scale_fill_gradient(low="lightgray",high = "red")+
      labs(fill=paste0(Gene))
    
  }else{
    stop("Wrong Shape")
  }
}

#PlotExpressionV2

PlotExpressionV2 <- function(
    barcodes,
    Gene,
    ptsize = 2.5,
    shape = "circle",
    colors = c("grey80", "yellow", "red"), # 默认使用一个更好看的色阶
    legend_title = Gene # 默认图例标题就是基因/分数名
) {
  
  # 检查'Gene'列是否存在
  if (!Gene %in% colnames(barcodes)) {
    stop("在数据框中找不到指定的 'Gene' 列: ", Gene)
  }
  
  # 将要绘制的列复制到一个通用列名'Expression'
  barcodes$Expression <- barcodes[[Gene]]
  
  # 根据形状选择美学 (color for circles, fill for squares)
  if (shape == "circle") {
    p <- ggplot(barcodes, aes(x = imagecol_scaled, y = -imagerow_scaled, color = Expression)) +
      geom_scattermore(pointsize = ptsize, pixels = c(2000, 2000)) +
      scale_color_gradientn(
        colors = colors,
        name = legend_title,
        na.value = "lightgrey" # 为NA值指定一个颜色
      )
  } else if (shape == "square") {
    p <- ggplot(barcodes, aes(x = imagecol_scaled, y = -imagerow_scaled, fill = Expression)) +
      geom_point(shape = 22, size = ptsize, color = alpha("black", 0), stroke = 0.25) +
      scale_fill_gradientn(
        colors = colors,
        name = legend_title,
        na.value = "lightgrey"
      )
  } else {
    stop("无效的形状。请选择 'circle' 或 'square'。")
  }
  
  # 添加通用的主题和美化
  p <- p +
    coord_fixed() +
    xlab("") +
    ylab("") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}



CalculateBivariateScore_KDE <- function(
    object,
    score_cols,
    new_score_name = "Coexpression_Score_KDE",
    rescale_input = TRUE
) {
  
  # --- 1. 参数和包检查 ---
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("请先安装 'MASS' 包: install.packages('MASS')")
  }
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("请先安装 'fields' 包: install.packages('fields')")
  }
  if (length(score_cols) != 2 || !all(score_cols %in% colnames(object@meta.data))) {
    stop("'score_cols' 必须是一个包含两个、且都存在于meta.data中的列名的向量。")
  }
  
  cat("--- 开始基于2D KDE计算共表达分数 ---\n")
  cat("  - 输入分数:", paste(score_cols, collapse=" 和 "), "\n")
  
  # --- 2. 提取并准备分数数据 ---
  scores_df <- object@meta.data[, score_cols]
  
  if (rescale_input) {
    cat("  - 正在将输入分数缩放到 0-1 范围...\n")
    scores_df[, 1] <- scales::rescale(scores_df[, 1], to = c(0, 1))
    scores_df[, 2] <- scales::rescale(scores_df[, 2], to = c(0, 1))
  }
  
  # --- 3. 计算二维核密度估计 (KDE) ---
  cat("  - 正在计算二维核密度...\n")
  # 使用 MASS::kde2d 计算密度地形图
  # n = 100 是一个比较均衡的网格精度选择
  kde_result <- MASS::kde2d(scores_df[[1]], scores_df[[2]], n = 100)
  
  # --- 4. 将密度值插值回每个原始点 ---
  cat("  - 正在将密度值插值回每个点...\n")
  # 使用 fields::interp.surface 来根据每个点的(x,y)坐标，查找其在密度图上的“海拔”
  interpolated_density <- fields::interp.surface(kde_result, as.matrix(scores_df))
  
  # --- 5. 将新的KDE分数添加回Seurat对象 ---
  object <- AddMetaData(object, metadata = interpolated_density, col.name = new_score_name)
  
  cat("--- 计算完成！新的分数已添加到 '", new_score_name, "' 列中。\n")
  
  return(object)
}





# Create an additional color palette to be used
ColorsClusters<-function()
{
  x <- c("128 128 128",
         "89 106 55",
         "150 86 53",
         "116 20 12",
         "42 98 24",
         "128 128 38",
         "69 61 134",
         "96 177 119",
         "55 126 127",
         "85 128 176",
         "4 0 123",
         "165 204 79",
         "210 167 65",
         "127 23 134",
         "234 51 35",
         "93 203 207",
         "240 146 53",
         "254 255 84",
         "183 45 130",
         "12 0 197",
         "117 251 76",
         "117 251 141",
         "202 49 66",
         "86 188 249",
         "232 167 108",
         "147 44 231",
         "190 253 91",
         "234 51 247",
         "69 142 247",
         "205 118 147",
         "238 230 151",
         "234 134 119",
         "212 162 217",
         "182 215 228",
         "120 105 230",
         "224 135 232",
         "175 249 162",
         "160 251 214",
         "250 228 200")
  
  
  
  Cols<-sapply(strsplit(x, " "), function(x)
    rgb(x[1], x[2], x[3], maxColorValue=255))
  
  return(Cols)
  
}

# Function to select the barcodes that are within a given distance from a given cell type (cluster)
SelectPeripheryDiscrete<-function(bcs,CellType,distance=50,PATH)
{
  
  SelectedBCs<-bcs %>% filter(DeconvolutionLabel1==CellType)
  
  Result<-GetSlice(SelectedBCs$barcode,distance,bcs,PATH,CellT=CellType)
  
  if(length(distance)>1)
  {
    for(jj in 1:length(Result))
    {
      Result[[jj]]<-Result[[jj]][Result[[jj]]%!in%SelectedBCs$barcode]
    }
    
    return(Result)
    
  }else{
    Result<-Result[Result%!in%SelectedBCs$barcode]
    
    return(Result)
  }
  
  
}

# Equivalent to GetSquare but for a circle instead.
GetSlice<-function(Spot,SizeMicrons,BarcodeDF,PATH,CellT=NA,size="008um")
{
  
  path_scales <- paste0(PATH, "/binned_outputs/square_",size,"/spatial/scalefactors_json.json")
  scales <- rjson::fromJSON(file = path_scales)
  Scale<-(SizeMicrons*scales$spot_diameter_fullres)/as.numeric(unlist(strsplit(size,"um"))[1])
  
  Index<-match(Spot,BarcodeDF$barcode)
  Result<-vector("list",length=length(Index))
  
  for(jj in 1:length(Index))
  {
    Distance<-sqrt(((BarcodeDF$imagecol-BarcodeDF$imagecol[Index[jj]])^2) + ((BarcodeDF$imagerow-BarcodeDF$imagerow[Index[jj]])^2))
    BarcodeDF$Distance<-Distance
    
    if(!is.na(CellT))
    {
      ValTh <- sum(BarcodeDF$DeconvolutionLabel1[BarcodeDF$Distance<min(Scale)]==CellT,na.rm = T)
      if(ValTh < 25)
      {
        next
      }
    }
    
    if(length(Scale)>1)
    {
      Result[[jj]]<-lapply(Scale,function(X){return(BarcodeDF$barcode[BarcodeDF$Distance < X])})
    }else{
      
      Result[[jj]]<-BarcodeDF$barcode[BarcodeDF$Distance < Scale]
    }
    
  }
  
  if(length(Scale)>1)
  {
    Rxx<-vector("list",length=length(Scale))
    names(Rxx)<-as.character(SizeMicrons)
    
    for(ii in 1:length(Scale))
    {
      Rxx[[ii]]<-lapply(Result,function(X){return(X[[ii]])})
      Rxx[[ii]]<-unique(unlist(Rxx[[ii]]))
    }
    
    return(Rxx)
    
  }else{
    Result<-unique(unlist(Result))
    return(Result)
  }
  
}

# Function to plot enrichR results as a barplot
EnrichRBarPlot<-function(Markers,DataBase,TermsX=10,PTh=1e-3,GO=F,colsT=c("firebrick1","dodgerblue"))
{
  
  Genes<-split(Markers$gene,Markers$cluster)
  ResA<-enrichr(Genes[[1]],DataBase)[[1]]
  ResA$Cluster<-names(Genes)[1]
  ResB<-enrichr(Genes[[2]],DataBase)[[1]]
  ResB$Cluster<-names(Genes)[2]
  
  Result<-rbind(ResA,ResB)
  
  if(GO)
  {
    Result$Term<-trimws(sapply(strsplit(Result$Term,"[(]"),function(X){return(X[1])}))
  }
  
  Result<-Result[Result$P.value<PTh,]
  
  Result<-Result %>% group_by(Cluster) %>% slice_max(order_by = P.value, n = TermsX)
  Result$FDR<--log(Result$Adjusted.P.value)
  Result$FDR<-Result$FDR*ifelse(as.numeric(as.factor(Result$Cluster))==1,1,-1)
  Result<-Result[order(Result$FDR),]
  Result$Term<-factor(Result$Term,levels=unique(Result$Term))
  
  Plot<-ggplot(Result,aes(x=Term,y=FDR,fill=Cluster))+geom_bar(stat="identity")+
    theme_classic()+scale_fill_manual(values=colsT)+
    xlab("")+ylab("log(Adj. pvalue)")+theme(axis.title.y=element_blank(),
                                            axis.text.y=element_blank(),
                                            axis.ticks.y=element_blank(),
                                            axis.line.y = element_blank(),
                                            axis.text.x = element_text(face="bold"))+
    geom_text(aes(label = Term,hjust = ifelse(FDR < 0, 0, 1),vjust = 0.5),y=0,size=3)+coord_flip()
  
  return(Plot)
  
}

# Function to create density map and identifies enriched region of a given cell type.
# BarcodeSet is used to restrict the area to a given collection of barcodes, if missing then
# the whole section is used.
PlotDensity<-function(DF,CellType,nBins=3,ptsize=3,Tumor=NA,BarcodeSet=NA)
{
  require(wesanderson)
  
  if(length(CellType)>1)
  {
    stop("Pass only 1 cell type to CellType argument")
  }
  
  # Create DF for density2d
  if(all(!is.na(BarcodeSet)))
  {
    DF2<-DF[DF$tissue==1 & DF$DeconvolutionLabel1%in%CellType & DF$barcode %in% BarcodeSet,]
    DF2$Grouping<-DF2$DeconvolutionLabel1
  }else{
    DF2<-DF[DF$tissue==1 & DF$DeconvolutionLabel1%in%CellType,]
    DF2$Grouping<-DF2$DeconvolutionLabel1
  }
  
  
  DF3 <- DF[DF$tissue==1  & DF$DeconvolutionLabel1==Tumor,]
  DF3 <- DF3 %>% na.omit()
  DF3$XX <- DF3$imagecol_scaled
  DF3$YY <- DF3$imagerow_scaled
  
  DF4 <- DF[DF$tissue==1  & DF$DeconvolutionLabel1%in%CellType,]
  DF4 <- DF4 %>% na.omit()
  DF4$XX <- DF4$imagecol_scaled
  DF4$YY <- DF4$imagerow_scaled
  DF4$Grouping<-DF4$DeconvolutionLabel1
  
  PlotX<-DF %>% filter(tissue == "1") %>% na.omit() %>% 
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled)) +
    geom_scattermore(pointsize = ptsize,pixels = rep(2000,2),col="lightgray")+
    geom_scattermore(data=DF3,pointsize=ptsize,pixels = rep(2000,2),col="gray65")+
    geom_scattermore(data=DF4,pointsize=ptsize,pixels = rep(2000,2),col="red")+
    geom_density_2d(data=DF2,bins=nBins,linewidth=1.2,linetype=1,aes(colour=after_stat(level)),
                    contour_var = "ndensity")+
    scale_colour_gradientn(colours=wes_palette("Zissou1", 20, type = "continuous"),limit=c(0,1))+
    #scale_colour_gradient2(low="dodgerblue",mid = "dodgerblue2" ,high="dodgerblue4",limit=c(0,1))+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    ggtitle(CellType)+
    labs(color="Scaled Density")
  
  return(PlotX)
  
}

# Collect the barcodes that are within an enriched for a given celltype/cluster
SelectEnrichedRegion<-function(CellType,bcs,PATH,N=5,Area=200)
{
  bcs <- bcs %>% filter(tissue==1)
  
  bcsCT<-bcs %>% filter (DeconvolutionLabel1%in%CellType)
  
  
  Kernel<-MASS::kde2d(bcsCT %>% pull(imagecol),
                      bcsCT %>% pull(imagerow),
                      n = 200)
  
  Peaks<-findPeaks3D(Kernel$z,N)
  Peaks<-data.frame(X=sapply(Peaks,function(X){return(X$x)}),
                    Y=sapply(Peaks,function(X){return(X$y)}),
                    Z=sapply(Peaks,function(X){return(X$z)}))
  
  Peaks$X<-Kernel$x[Peaks$X]
  Peaks$Y<-Kernel$y[Peaks$Y]
  
  Result<-vector("list",length=nrow(Peaks))
  
  for(jj in 1:nrow(Peaks))
  {
    BCTmp<-bcs
    BCTmp<-GetDistance(BCTmp,PATH,XY=c(Peaks$X[jj],Peaks$Y[jj]))
    BCTmp<-BCTmp[order(BCTmp$DistanceMicrons),]
    BarcodeX<-BCTmp$barcode[order(BCTmp$DistanceMicrons)[1]]
    Result[[jj]]<-GetSlice(BarcodeX,Area,BCTmp,PATH)
    
  }
  
  return(Result)
  
}

# Find peaks of a 3d function. Used by SelectEnrichedRegion
findPeaks3D <- function(matrix, N) {
  # Function to check if a point is a local maximum
  isLocalMaximum <- function(mat, x, y) {
    neighbors <- c(
      mat[x-1, y], mat[x+1, y], mat[x, y-1], mat[x, y+1],
      mat[x-1, y-1], mat[x-1, y+1], mat[x+1, y-1], mat[x+1, y+1]
    )
    # Remove NA values (edges)
    neighbors <- neighbors[!is.na(neighbors)]
    return(all(mat[x, y] > neighbors))
  }
  
  numRows <- nrow(matrix)
  numCols <- ncol(matrix)
  
  # List to store peaks
  peaks <- list()
  
  for (x in 2:(numRows - 1)) {
    for (y in 2:(numCols - 1)) {
      if (isLocalMaximum(matrix, x, y)) {
        peaks <- c(peaks, list(list(x = x, y = y, z = matrix[x, y])))
      }
    }
  }
  
  # Sort peaks by height and select top N
  peaks <- peaks[order(sapply(peaks, function(peak) -peak$z))]
  if (length(peaks) > N) {
    peaks <- peaks[1:N]
  }
  
  return(peaks)
}

# Function to get the distance to a given barcode or XY coordinate
GetDistance<-function(BCS,PATH,barcode=NA,XY=NA,size="008um")
{
  
  # Transform Scale
  path_scales <- paste0(PATH, "/binned_outputs/square_",size,"/spatial/scalefactors_json.json")
  scales <- rjson::fromJSON(file = path_scales)
  
  # Select Either BC or XY coordinate
  
  if(all(!is.na(barcode),all(!is.na(XY))))
  {
    stop("Either barcode or XY coordinate")
  }else if(!is.na(barcode) & all(is.na(XY)))
  {
    BC_Center<-barcode
    Index<-match(BC_Center,BCS$barcode)
    Distance<-sqrt(((BCS$imagecol-BCS$imagecol[Index])^2) + ((BCS$imagerow-BCS$imagerow[Index])^2))
    DistanceMicrons<-(Distance*8)/scales$spot_diameter_fullres
    BCS$DistanceMicrons<-DistanceMicrons
    
  }else if(is.na(barcode) & all(!is.na(XY)))
  {
    BC_Center<-barcode
    Index<-match(BC_Center,BCS$barcode)
    Distance<-sqrt(((BCS$imagecol-XY[1])^2) + ((BCS$imagerow-XY[2])^2))
    DistanceMicrons<-(Distance*8)/scales$spot_diameter_fullres
    BCS$DistanceMicrons<-DistanceMicrons
  }
  
  return(BCS)
  
}

# Function to find the Tumor (cluster) bins that are within a given distance from a selection of barcodes.
SelectTumor<-function(bcDF,Tumor,barcodes,PATH,DistVal=50,size="008um")
{
  # Get scale factors
  path_scales <- paste0(PATH, "/binned_outputs/square_",size,"/spatial/scalefactors_json.json")
  scales <- rjson::fromJSON(file = path_scales)
  
  # Get the center spot from the given barcodes
  RegionDF <- bcDF %>% filter(barcode%in%barcodes) %>% dplyr::select(barcode,imagerow,imagecol,DeconvolutionLabel1)
  
  MeanSpots<-RegionDF %>% summarise(Xval=(max(imagecol)+min(imagecol))/2,
                                    Yval=(max(imagerow)+min(imagerow))/2)
  
  Distance<-sqrt(((bcDF$imagecol-MeanSpots$Xval)^2) + ((bcDF$imagerow-MeanSpots$Yval)^2))
  CenterSpot<-bcDF$barcode[which.min(Distance)]
  
  # Select Slice 200 microns and Subset the data.frame
  SliceRegion<-GetSlice(CenterSpot,350,bcDF,PATH)
  
  RegionDF <- bcDF %>% filter(barcode%in%SliceRegion | barcode %in%barcodes) %>% dplyr::select(barcode,imagerow,imagecol,DeconvolutionLabel1)
  
  #Get Distance between given barcodes and Tumor Spots in the RegionDF
  RegionDF<- RegionDF %>% filter(DeconvolutionLabel1==Tumor | barcode %in% barcodes)
  TumBC<-RegionDF %>% filter(DeconvolutionLabel1==Tumor) %>% pull(barcode)
  
  XY_Data<-RegionDF[,c("imagecol","imagerow")]
  DistMat<-distances::distances(as.matrix(XY_Data))
  DistMat<-DistMat[match(barcodes,RegionDF$barcode),match(TumBC,RegionDF$barcode)]
  rownames(DistMat)<-barcodes
  colnames(DistMat)<-TumBC
  
  DistMat<-(DistMat*8)/scales$spot_diameter_fullres
  
  # Select for each Region spot the closest tumor Spot
  iix<-apply(DistMat,2,function(X){any(X<DistVal)})
  ClosestTumorSpot<-colnames(DistMat)[iix]
  
  return(ClosestTumorSpot)
  
}




SelectNeighboringCells <- function(
    object, # 输入现在是一个Seurat对象
    source_barcodes,
    target_celltype,
    group_by_col, # 用于识别靶细胞类型的列名
    PATH,
    DistVal = 50,
    search_radius_microns = 350 # 内部搜索半径，用于优化
) {
  
  # --- 1. 参数和数据检查 ---
  if (!inherits(object, "Seurat")) stop("'object' 必须是一个Seurat对象。")
  if (!group_by_col %in% colnames(object@meta.data)) {
    stop("在meta.data中找不到指定的分组列: '", group_by_col, "'")
  }
  
  cat("--- 正在寻找 '", target_celltype, "' 细胞，它们邻近于指定的源区域 ---\n")
  bcDF <- object@meta.data
  bcDF$barcode <- rownames(bcDF)
  
  # --- 2. 找到源区域的几何中心 (与之前类似) ---
  RegionDF <- bcDF %>% filter(barcode %in% source_barcodes)
  if(nrow(RegionDF) == 0) {
    warning("提供的'source_barcodes'在对象中一个都找不到。")
    return(character(0))
  }
  
  MeanSpots <- RegionDF %>% summarise(
    Xval = mean(imagecol, na.rm = TRUE), # 使用平均质心，更准确
    Yval = mean(imagerow, na.rm = TRUE)
  )
  
  # --- 3. 以区域中心为圆心，切出一个更大的“搜索范围” (优化步骤) ---
  # 我们需要一个包含所有点坐标的数据框来调用GetSlice
  full_bcs_df <- GenerateCellbinSampleData_Final(PATH)$bcs
  
  Distance <- sqrt(((full_bcs_df$imagecol - MeanSpots$Xval)^2) + ((full_bcs_df$imagerow - MeanSpots$Yval)^2))
  CenterSpot <- full_bcs_df$barcode[which.min(Distance)]
  
  SliceRegion_barcodes <- GetSlice(CenterSpot, search_radius_microns, full_bcs_df, PATH)
  
  # --- 4. 准备精确距离计算的数据 (核心修改) ---
  # 我们只在搜索范围内进行操作
  search_df <- bcDF %>% filter(barcode %in% SliceRegion_barcodes)
  
  # 筛选出我们的两种目标群体：源细胞 和 靶细胞
  target_barcodes_in_slice <- search_df$barcode[search_df[[group_by_col]] == target_celltype]
  source_barcodes_in_slice <- intersect(search_df$barcode, source_barcodes)
  
  # 准备坐标矩阵
  combined_barcodes <- unique(c(source_barcodes_in_slice, target_barcodes_in_slice))
  XY_Data <- search_df[combined_barcodes, c("imagecol", "imagerow")]
  
  if (length(target_barcodes_in_slice) == 0) {
    cat("  - 在搜索范围内没有找到任何 '", target_celltype, "' 细胞。\n")
    return(character(0))
  }
  
  # --- 5. 计算精确的距离矩阵 (与之前类似) ---
  cat("  - 正在计算", length(source_barcodes_in_slice), "个源细胞与", 
      length(target_barcodes_in_slice), "个靶细胞之间的距离...\n")
  
  DistMat <- as.matrix(dist(XY_Data))
  # 子集化距离矩阵
  DistMat <- DistMat[source_barcodes_in_slice, target_barcodes_in_slice, drop=FALSE]
  
  # 单位转换 (需要scalefactors)
  scales_path <- file.path(PATH, "segmented_outputs/spatial/scalefactors_json.json")
  scales <- rjson::fromJSON(file = scales_path)
  # 1像素的物理尺寸（微米）
  microns_per_pixel <- scales$microns_per_pixel 
  DistMat_microns <- DistMat * microns_per_pixel
  
  # --- 6. 筛选出邻近的靶细胞 (与之前类似) ---
  # 找出哪些靶细胞(列)，与源细胞(行)的最小距离小于阈值
  is_close_enough <- apply(DistMat_microns, 2, function(col_distances) {
    any(col_distances < DistVal)
  })
  
  ClosestTargetCells <- names(is_close_enough[is_close_enough])
  
  cat("--- 找到了", length(ClosestTargetCells), "个邻近的 '", target_celltype, "' 细胞。\n")
  return(ClosestTargetCells)
}



# Plot colocalization of two cell types in  a sample (used for CD4 and CD8 t cells)
PlotColocalization<-function(DF,CellType1,CelltypesID,nBins=3,ptsize=3,Tumor=NA,BarcodeSet=NA,colX=c("darkorchid","gold"),option="CD8")
{
  require(wesanderson)
  require(ggnewscale)
  
  # Create DF for density2d
  if(all(!is.na(BarcodeSet)))
  {
    DF2<-DF[DF$tissue==1 & DF$DeconvolutionLabel1==CellType1 & DF$barcode %in% BarcodeSet,]
    DF2$Grouping<-DF2$DeconvolutionLabel1
    
  }else{
    DF2<-DF[DF$tissue==1 & DF$DeconvolutionLabel1==CellType1,]
    DF2$Grouping<-DF2$DeconvolutionLabel1
  }
  
  
  DF3 <- DF[DF$tissue==1  & DF$DeconvolutionLabel1==Tumor,]
  DF3 <- DF3 %>% na.omit()
  DF3$XX <- DF3$imagecol_scaled
  DF3$YY <- DF3$imagerow_scaled
  
  DF4A <- DF[DF$tissue==1  & DF$DeconvolutionLabel1==CelltypesID[1],]
  DF4A <- DF4A %>% na.omit()
  DF4A$XX <- DF4A$imagecol_scaled
  DF4A$YY <- DF4A$imagerow_scaled
  DF4A$Grouping<-DF4A$DeconvolutionLabel1
  
  DF4B <- DF[DF$tissue==1  & DF$DeconvolutionLabel1==CelltypesID[2],]
  DF4B <- DF4B %>% na.omit()
  DF4B$XX <- DF4B$imagecol_scaled
  DF4B$YY <- DF4B$imagerow_scaled
  DF4B$Grouping<-DF4B$DeconvolutionLabel1
  
  if(option=="CD8")
  {
    ColsGrad<-brewer.pal(9,"Blues")
  }else{
    ColsGrad<-brewer.pal(9,"Greens")
  }
  
  PlotX<-DF %>% filter(tissue == "1") %>% na.omit() %>% 
    ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled)) +
    geom_scattermore(pointsize = ptsize,pixels = rep(2000,2),col="lightgray")+
    geom_scattermore(data=DF3,pointsize=ptsize,pixels = rep(2000,2),col="gray65")+
    geom_scattermore(data=DF4A,pointsize=ptsize+2,pixels = rep(2000,2),col=colX[1])+
    geom_scattermore(data=DF4B,pointsize=ptsize+2,pixels = rep(2000,2),col=colX[2])+
    geom_scattermore(data=DF2,pointsize=ptsize+2,pixels = rep(2000,2),col=ifelse(option=="CD8","dodgerblue","forestgreen"))+
    geom_density_2d(data=DF2,bins=nBins,linewidth=1.2,linetype=1,aes(colour=after_stat(level)),contour_var = "ndensity")+
    scale_colour_gradientn(colours=ColsGrad,limit=c(0,1))+
    labs(color="Scaled Density")+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  
  return(PlotX)
  
}

# Transform barcode names between sizes
TransformBarcodes<-function(Barcodes,SizeOriginal,SizeNew)
{
  SizeO<-str_pad(paste0(SizeOriginal,"um"),5,pad="0")
  SizeN<-str_pad(paste0(SizeNew,"um"),5,pad="0")
  
  Barcodes<-strsplit(Barcodes,"[_-]")
  Barcodes<-data.frame(s="s",size=SizeO,
                       X=as.numeric(sapply(Barcodes,function(X){return(X[3])})),
                       Y=as.numeric(sapply(Barcodes,function(X){return(X[4])})),
                       end="-1")
  
  
  nBinsSide <- SizeOriginal / SizeNew
  
  Xmins <- Barcodes$X * nBinsSide
  Xmaxs <- (Barcodes$X * nBinsSide) + (nBinsSide - 1)
  Ymins <- Barcodes$Y * nBinsSide
  Ymaxs <- (Barcodes$Y * nBinsSide) + (nBinsSide - 1)
  
  Result<-lapply(1:nrow(Barcodes),function(jj)
  {
    Rx <- expand.grid(x = Xmins[jj]:Xmaxs[jj], y = Ymins[jj]:Ymaxs[jj])
    Rx$X <- Barcodes$X[jj]
    Rx$Y <- Barcodes$Y[jj]
    return(Rx)
  })
  
  Result <- do.call(rbind, Result)
  
  OldBC<-paste0("s_",SizeO,"_",str_pad(Result$X,5,pad="0"),"_",str_pad(Result$Y,5,pad="0"),"-1")
  NewBC<-paste0("s_",SizeN,"_",str_pad(Result$x,5,pad="0"),"_",str_pad(Result$y,5,pad="0"),"-1")
  
  Result<-data.frame(Original=OldBC,Transformed=NewBC)
  #Result<-split(Result$Transformed,Result$Original)
  
  return(Result)
  
}

# Used to provide parameters for the segmentation script
NucleiSegmentationScript<-function(BarcodeDF,TransformedDF)
{
  # Get Row and Column image limits for Nuclei Segmentation
  Rx<-data.frame(C1=round(min(BarcodeDF$imagecol[match(unique(TransformedDF$Original),BarcodeDF$barcode)])),
                 C2=round(max(BarcodeDF$imagecol[match(unique(TransformedDF$Original),BarcodeDF$barcode)])),
                 R1=round(min(BarcodeDF$imagerow[match(unique(TransformedDF$Original),BarcodeDF$barcode)])),
                 R2=round(max(BarcodeDF$imagerow[match(unique(TransformedDF$Original),BarcodeDF$barcode)])))
  
  Script<-paste0(" ./NucleiSegmentation.py -i PATH/image.btf -r1 ",Rx$R1," -r2 ",Rx$R2," -c1 ",Rx$C1," -c2 ",Rx$C2," -x /PATH/TO/square_00Xum/ -o PATH/TO/OUTPUT")
  
  message(Script)
}

# Plotting function to at image
geom_spatial<-function (mapping = NULL, data = NULL, stat = "identity", position = "identity", 
                        na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, ...) 
{
  GeomCustom <- ggproto("GeomCustom", Geom, setup_data = function(self, 
                                                                  data, params) {
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    data
  }, draw_group = function(data, panel_scales, coord) {
    vp <- grid::viewport(x = data$x, y = data$y)
    g <- grid::editGrob(data$grob[[1]], vp = vp)
    ggplot2:::ggname("geom_spatial", g)
  }, required_aes = c("grob", "x", "y"))
  layer(geom = GeomCustom, mapping = mapping, data = data, 
        stat = stat, position = position, show.legend = show.legend, 
        inherit.aes = inherit.aes, params = list(na.rm = na.rm, 
                                                 ...))
}

# Negate %in% operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# Plot C-C Communication 
PlotInteraction<-function(LianaRes,SourceCTs=NA,TargetCTs=NA,N=2,Gap=5,ColorsUser=NA,scale=FALSE,alpha=0.2)
{
  # Define Colors if not given
  if(all(is.na(ColorsUser)))
  {
    ColorsTracks<-paletteer::paletteer_d("ggsci::default_igv")
  }else{
    
    ColorsTracks<-ColorsUser
  }
  
  # Generate Frequncy Matrix to be used
  Freqs<-table(LianaRes$source,LianaRes$target)
  Freqs<-matrix(Freqs, ncol = ncol(Freqs), dimnames = dimnames(Freqs))
  
  # Filter source to the selected cell types if any
  if(all(!is.na(SourceCTs) & SourceCTs %in% rownames(Freqs)))
  {
    Freqs<-Freqs[SourceCTs,,drop=F]
  }else{
    
    SourceCTs<-rownames(Freqs)
  }
  
  # Filter interaction with at least N entries
  Freqs[Freqs<N]<-0
  
  if(all(Freqs==0))
  {
    stop("No interactions with given parameters")
  }
  
  # Keep targets with at least 1 interaction
  Freqs<-Freqs[,colSums(Freqs)>0,drop=F]
  
  if(all(Freqs==0))
  {
    stop("No interactions with given parameters")
  }
  
  # Order the CTs (source and target)
  orderTracks <- c(SourceCTs,sort(colnames(Freqs)[colnames(Freqs) %!in% SourceCTs]))
  
  # Define Gaps to visualize better source and sinks
  if(scale)
  {
    GapsTracks<-c(rep(5,nrow(Freqs)-1),50,rep(5,(length(orderTracks)-length(SourceCTs))-1),50)
  }else{
    GapsTracks<-c(rep(5,nrow(Freqs)-1),20,rep(5,(length(orderTracks)-length(SourceCTs))-1),20)
  }
  
  # Generate Color Matrix
  
  ColLinks<-ColorsTracks[match(rownames(Freqs),names(ColorsTracks))]
  ColsAlpha<-adjustcolor(ColLinks,alpha.f = alpha)
  ColorMatrix<-matrix(rep(ColsAlpha,each=ncol(Freqs)),nrow=nrow(Freqs),ncol=ncol(Freqs),dimnames = dimnames(Freqs),byrow = T)
  
  # Used only if we want to highlight specific relationships
  if(all(!is.na(TargetCTs)))
  {
    RowsH<-expand.grid(SourceCTs,TargetCTs)
    RowsH$Col<-ColorsTracks[as.vector(RowsH$Var1)]
    
    Pos<-cbind(match(RowsH$Var1,rownames(Freqs)),match(RowsH$Var2,colnames(Freqs)))
    index<-!is.na(Pos[,1]) & !is.na(Pos[,2])
    Pos<-Pos[index,]
    RowsH<-RowsH[index,]
    
    ColorMatrix[Pos]<-RowsH$Col
    
  }
  
  circos.par(gap.after = GapsTracks,start.degree = -90)
  
  chordDiagram(Freqs, order = orderTracks, annotationTrack = c("grid","name"),grid.col=ColorsTracks,
               direction.type = c("diffHeight", "arrows"),directional = 1,link.arr.type = "big.arrow",
               scale=scale,col = ColorMatrix)
  
  circos.clear()
  
}

## Extra color Palette
ColorsExtra<-function()
{
  Cols<-c("#5580B0","#56BCF9","steelblue","#7F1786","#A5CC4F","#D2A741",
          "#E087E8","#B6D7E4","#932CE7","darkred", "salmon","#458EF7", 
          "#5DCBCF","#E8A76C","#A0FBD6","#75FB4C",  "black","#EA8677", 
          "#965635","#CA3142","#75FB8D","#7869E6","#AFF9A2","#60B177", 
          "#FEFF54","#808080","#EA33F7","#EA3323","#0C00C5","#377E7F", 
          "#2A6218","#F09235","#EEE697","#453D86","#CD7693","#74140C", 
          "#808026","#FAE4C8","#BEFD5B","#D4A2D9","#B72D82","#596A37","#04007B") 
  
  names(Cols)<-c('Bcells-0','Bcells-1','Bcells-2','Bcells-3','Bcells-4','Endothelial-0','Endothelial-1','Endothelial-2','Endothelial-3',
                 'Fibroblast-0','Fibroblast-1','Fibroblast-2','IntestinalEpithelial-0','IntestinalEpithelial-1','IntestinalEpithelial-2',
                 'IntestinalEpithelial-3','IntestinalEpithelial-4','IntestinalEpithelial-5','Myeloid-0','Myeloid-1','Myeloid-2',
                 'Neuronal-0','Neuronal-1','SmoothMuscle-0','SmoothMuscle-1','SmoothMuscle-2','SmoothMuscle-3','SmoothMuscle-4',
                 'Tcells-0','Tcells-1','Tcells-2','Tcells-3','Tumor-0','Tumor-1','Tumor-2','Tumor-3','Unknown-0','Unknown-1',
                 'Unknown-2','Unknown-3','Unknown-4','Unknown-5','Unknown-6')
  
  return(Cols)
}

EnrichRDotPlot<-function(Markers,Database,TermsX=10,N=5)
{
  Cluss<-as.vector(unique(Markers$cluster))
  ResultTerms<-vector("list",length = length(Cluss))
  for(jj in 1:length(Cluss))
  {
    
    Res<-enrichr(Markers[Markers$cluster==Cluss[jj],"gene"],Database)
    RxTmp<-Res[[1]]
    RxTmp$Cluster<-Cluss[jj]
    ResultTerms[[jj]]<-RxTmp
    
  }
  
  ResultTerms<-do.call(rbind,ResultTerms)
  ResultTerms<-ResultTerms[ResultTerms$P.value<1e-3,]
  
  Terms_Result<-ResultTerms %>%  arrange(Adjusted.P.value) %>% group_by(Cluster) %>% slice(1:N) %>% pull(Term)
  Terms_Result<-unique(Terms_Result)
  
  iix<-c()
  for(jj in 1:length(Terms_Result))
  {
    
    iix<-c(iix,which(ResultTerms[,1]==Terms_Result[jj]))
    
  }
  
  Terms_Reduced<-ResultTerms[iix,]
  Terms_Reduced<-Terms_Reduced[!is.na(Terms_Reduced$Term),]
  Terms_Reduced$Term<-factor(Terms_Reduced$Term,levels = unique(Terms_Reduced$Term))
  
  
  PlotX<-ggplot(Terms_Reduced,aes(x=Cluster,y=Term,colour=-log10(P.value)))+geom_point(aes(size=Odds.Ratio))+scale_color_viridis()+theme_classic()
  
  return(PlotX)
  
  
  
}

# Plot C-C Communication 
PlotInteractionGraph<-function(LianaRes,CellTypes=NA,Colors=NA,TitlePlot=NULL)
{
  # Taken from CellChat [netVisual_circle]
  # https://github.com/sqjin/CellChat/blob/e4f68625b074247d619c2e488d33970cc531e17c/R/visualization.R#L1240
  
  # Parameters
  vertex.weight <- 10
  edge.width.max <- 8
  arrow.width <- 1
  arrow.size <- 0.6
  edge.curved <- 0.2
  shape <- 'circle'
  
  
  # Define Colors if not given
  if(all(is.na(Colors)))
  {
    ColorsTracks<-paletteer::paletteer_d("ggsci::default_igv")
  }else{
    
    ColorsTracks<-Colors
  }
  
  if(all(!is.na(CellTypes)))
  {
    LianaRes <- LianaRes %>% filter(source %in% CellTypes & target %in% CellTypes)
  }
  
  NetworkDF<-as.matrix(table(LianaRes$source,LianaRes$target))
  cells.level <- unique(c(rownames(NetworkDF),colnames(NetworkDF)))
  df.net <- reshape2::melt(NetworkDF, value.name = "value")
  colnames(df.net)[1:2] <- c("source","target")
  
  df.net$source <- factor(df.net$source, levels = cells.level)
  df.net$target <- factor(df.net$target, levels = cells.level)
  df.net$value[is.na(df.net$value)] <- 0
  
  NetworkDF <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  
  NetworkDF[is.na(NetworkDF)] <- 0
  
  g <- graph_from_adjacency_matrix(NetworkDF, mode = "directed", weighted = T)
  
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  
  coords<-layout_(g,in_circle())
  
  if(nrow(coords)!=1)
  {
    coords_scale=scale(coords)
    
  }else{
    
    coords_scale<-coords
  }
  
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-ColorsTracks[names(igraph::V(g))]
  igraph::V(g)$frame.color <- ColorsTracks[names(igraph::V(g))]
  igraph::V(g)$label.color <- "black"
  igraph::V(g)$label.cex<-1
  
  edge.weight.max <- max(igraph::E(g)$weight)
  igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-"black"
  igraph::E(g)$label.cex<-0.8
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],0.6)
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))
  
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=0.2, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  
  if (!is.null(TitlePlot)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  
  gg <- recordPlot()
  
  return(gg)
}

# Function to calculate Euclidean distance between two RGB colors
color_distance <- function(color1, color2) {
  sqrt(sum((color1 - color2)^2))
}

# Function to get the N most different colors
most_different_colors <- function(hex_colors, N) {
  # Convert HEX to RGB
  rgb_colors <- t(col2rgb(hex_colors))
  
  # Initialize the first color as the one farthest from the mean color
  mean_color <- colMeans(rgb_colors)
  distances_to_mean <- apply(rgb_colors, 1, function(x) color_distance(x, mean_color))
  selected_colors <- rgb_colors[which.max(distances_to_mean), , drop = FALSE]
  selected_indices <- which.max(distances_to_mean)
  
  # Select the remaining N-1 most different colors
  for (i in 2:N) {
    distances <- apply(rgb_colors, 1, function(x) min(sapply(1:nrow(selected_colors), function(j) color_distance(x, selected_colors[j, ]))))
    next_color_index <- which.max(distances)
    selected_colors <- rbind(selected_colors, rgb_colors[next_color_index, , drop = FALSE])
    selected_indices <- c(selected_indices, next_color_index)
  }
  
  # Return the HEX codes of the selected colors
  return(hex_colors[selected_indices])
}


SpatialAccuracy<-function(barcodes,Path,Genes)
{
  # Generate BC data.frame
  BCS<-GenerateSampleData(Path)$bcs
  BCS<- BCS %>% filter(tissue==1)
  
  # Read full H5 matrix
  Mat<-Read10X_h5(paste0(Path,"/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"))[,BCS$barcode]
  
  # Add nUMI column to data.frame
  BCS$nUMI<-colSums(Mat)
  
  # Generate full section Plot and add squares
  SqDF<-c()
  
  for(index in 1:length(barcodes))
  {
    Selection<-GetSquare(barcodes[index],500,BCS)
    BCS$IsSelection<-BCS$barcode%in%Selection
    SqDF<-rbind(SqDF,BCS %>% filter(IsSelection) %>% summarise(Xmin=min(imagecol_scaled),Xmax=max(imagecol_scaled),Ymin=min(-imagerow_scaled),Ymax=max(-imagerow_scaled),Group="Square",Label=LETTERS[index]))
  }
  
  RowTitle<-sapply(1:length(barcodes), function(x) paste(rep("i", x), collapse = ""))
  
  LabelDF<-data.frame(imagecol_scaled=apply(SqDF[,1:2],1,mean),
                      imagerow_scaled=apply(SqDF[,3:4],1,mean),
                      Label=RowTitle)
  
  PlotAll<-BCS %>% ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=log(nUMI+1))) + 
    geom_scattermore(pointsize = 2,pixels = rep(2000,2))+
    scale_color_viridis(option="B")+
    coord_cartesian(expand = FALSE) +
    xlab("") +
    ylab("") +
    theme_set(theme_bw(base_size = 10))+
    theme_minimal() +
    theme(axis.text = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position="bottom")+
    geom_rect(data=SqDF,aes(xmin=Xmin, xmax=Xmax, ymin=Ymax, ymax=Ymin),
              fill=NA,color="black",linewidth=1,inherit.aes = FALSE)+
    labs(color="log UMI")+geom_text(data=LabelDF,aes(x=imagecol_scaled,y=imagerow_scaled,label=Label,fontface=2),size=6,color="black")
  
  
  PlotResult<-vector("list",length=length(barcodes)*(length(Genes)))
  saveindex<-1
  
  MatrixIndex<-matrix(1:(length(Genes)*length(barcodes)),nrow=length(barcodes),ncol = length(Genes),byrow = T)
  LetterIndex<-MatrixIndex[,1]
  TitleIndex<-MatrixIndex[1,]
  
  for(index in 1:length(barcodes))
  {
    message(barcodes[index])
    Selection<-GetSquare(barcodes[index],500,BCS)
    BCS$IsSelection<-BCS$barcode%in%Selection
    
    for(listI in 1:length(Genes))
    {
      print(names(Genes)[listI])
      BCS$Markers<-colSums(Mat[Genes[[listI]],])
      
      PlotResult[[saveindex]]<-BCS %>% filter(IsSelection) %>%  ggplot(aes(x = imagecol_scaled, y = -imagerow_scaled,color=log(Markers+1)))+
        geom_point(shape=15,size=8)+scale_color_viridis(option="B",limits=c(0,3))+coord_cartesian(expand = FALSE) +
        xlab("") +
        ylab(ifelse(saveindex%in%LetterIndex,RowTitle[index],"")) +
        theme_set(theme_bw(base_size = 10))+
        theme_minimal() +
        theme(axis.text = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              axis.title.y = element_text(angle=0,size=14))+ggtitle(ifelse(saveindex %in% TitleIndex,names(Genes)[saveindex],""),
                                                                    subtitle = ifelse(saveindex %in% TitleIndex,paste(Genes[[saveindex]],collapse=","),""))+
        labs(color="log UMI")
      
      saveindex<-saveindex+1
      
    }
    
    
  }
  
  P2<-Reduce(`+`, PlotResult)+plot_layout(guides = 'collect')
  
  return(PlotAll+P2)
}





# ==== RunBANKSY_Cellbin
RunBANKSY_Cellbin <- function(
    sobj,
    batch_col       = "orig.ident",
    lambda_celltype = 0.2,
    lambda_domain   = 0.8,
    k_geom          = c(15, 30),
    npcs            = 30,
    res_celltype    = 1.2,
    res_domain      = 1.0,
    k_neighbors     = 50,
    agf_mode        = c("auto","on","off"),
    hvg_n           = NULL,
    row_limit_factor= 0.5
){
  agf_mode <- match.arg(agf_mode)
  stopifnot(inherits(sobj, "Seurat"))
  if (!all(c("imagecol","imagerow") %in% colnames(sobj@meta.data))) {
    stop("缺少 imagecol/imagerow；请先运行 GenerateCellbinSeurat()/GenerateCellbinSampleData()。")
  }
  
  suppressPackageStartupMessages({
    library(Banksy); library(SpatialExperiment)
    library(SummarizedExperiment); library(Seurat)
  })
  
  counts <- tryCatch(
    Seurat::GetAssayData(sobj, assay = "Spatial", slot = "counts"),
    error = function(e) Seurat::GetAssayData(sobj, assay = "Spatial", layer = "counts")
  )
  if (!is.null(hvg_n)) {
    sobj_tmp <- sobj; Seurat::DefaultAssay(sobj_tmp) <- "Spatial"
    sobj_tmp <- Seurat::NormalizeData(sobj_tmp)
    sobj_tmp <- Seurat::FindVariableFeatures(sobj_tmp, selection.method = "vst", nfeatures = hvg_n)
    hvgs <- intersect(Seurat::VariableFeatures(sobj_tmp), rownames(counts))
    counts <- counts[hvgs, , drop = FALSE]
    message("HVG 子集：", length(hvgs), " genes")
  }
  n_cells <- ncol(counts)
  if (!batch_col %in% colnames(sobj@meta.data)) batch_col <- "orig.ident"
  
  spe <- SpatialExperiment::SpatialExperiment(
    assays        = list(counts = counts),
    spatialCoords = cbind(x = sobj@meta.data$imagecol, y = sobj@meta.data$imagerow),
    colData       = S4Vectors::DataFrame(sobj@meta.data[, batch_col, drop = FALSE])
  )
  
  fm_cb <- tryCatch(formals(Banksy::computeBanksy), error = function(e) NULL)
  has_compute_agf <- !is.null(fm_cb) && "compute_agf" %in% names(fm_cb)
  has_nharm       <- !is.null(fm_cb) && "n_harmonics" %in% names(fm_cb)
  has_chunk       <- !is.null(fm_cb) && "chunk_size"  %in% names(fm_cb)
  has_rowlim      <- !is.null(fm_cb) && "row_limit_factor" %in% names(fm_cb)
  
  k_vec <- as.numeric(k_geom)
  safe_chunk <- max(1000L, min(nrow(counts),
                               floor((.Machine$integer.max * row_limit_factor) / max(1L, n_cells))))
  
  use_agf_flag <- switch(agf_mode,
                         auto = length(k_vec) > 1,
                         on   = TRUE,
                         off  = FALSE)
  
  cb_args <- list(assay_name = "counts", k_geom = k_vec)
  if (has_compute_agf) cb_args$compute_agf <- use_agf_flag
  if (has_nharm)       cb_args$n_harmonics <- length(k_vec)
  if (has_chunk)       cb_args$chunk_size  <- as.integer(safe_chunk)
  if (has_rowlim)      cb_args$row_limit_factor <- row_limit_factor
  
  ## ---- 这里改：把 x = spe 改为 se = spe（或直接放 spe 作为第一个匿名参数）----
  spe <- do.call(Banksy::computeBanksy, c(list(se = spe), cb_args))
  
  fm_pca <- tryCatch(formals(Banksy::runBanksyPCA), error = function(e) NULL)
  fm_clu <- tryCatch(formals(Banksy::clusterBanksy), error = function(e) NULL)
  has_use_agf_pca <- !is.null(fm_pca) && "use_agf" %in% names(fm_pca)
  has_use_agf_clu <- !is.null(fm_clu) && "use_agf" %in% names(fm_clu)
  
  pca_args <- list(se = spe, lambda = c(lambda_celltype, lambda_domain),
                   npcs = npcs, group = batch_col)   # <- 这里也用 se=
  if (has_use_agf_pca) pca_args$use_agf <- use_agf_flag
  spe <- do.call(Banksy::runBanksyPCA, pca_args)
  
  clu1_args <- list(se = spe, lambda = lambda_celltype,   # <- se=
                    algo = "leiden", k_neighbors = k_neighbors, resolution = res_celltype)
  if (has_use_agf_clu) clu1_args$use_agf <- use_agf_flag
  spe <- do.call(Banksy::clusterBanksy, clu1_args)
  
  clu2_args <- list(se = spe, lambda = lambda_domain,     # <- se=
                    algo = "leiden", k_neighbors = k_neighbors, resolution = res_domain)
  if (has_use_agf_clu) clu2_args$use_agf <- use_agf_flag
  spe <- do.call(Banksy::clusterBanksy, clu2_args)
  
  cn <- Banksy::clusterNames(spe)
  pick_by_lambda <- function(lam_str) {
    cand <- grep(paste0("lam", gsub("\\.", "\\\\.", lam_str)), cn, value = TRUE)
    if (!length(cand)) stop("未找到匹配列（", lam_str, "）：可用列为：", paste(cn, collapse=", "))
    sizes <- sapply(cand, function(cc) length(unique(SummarizedExperiment::colData(spe)[[cc]])))
    cand <- cand[sizes > 1]
    if (!length(cand)) stop("匹配到 lam", lam_str, " 的列，但簇数都 <= 1。")
    cand[which.max(sizes[cand])]
  }
  ct_col <- pick_by_lambda(gsub("\\.?0+$","", as.character(lambda_celltype)))
  dm_col <- pick_by_lambda(gsub("\\.?0+$","", as.character(lambda_domain)))
  
  ord <- match(colnames(sobj), colnames(spe))
  if (anyNA(ord)) {
    warning("Seurat 与 SPE 列名不完全一致，将按 SPE 顺序写回；请确认列名。")
    ord <- seq_len(ncol(spe))
  }
  ct_vec <- as.character(SummarizedExperiment::colData(spe)[[ct_col]][ord])
  dm_vec <- as.character(SummarizedExperiment::colData(spe)[[dm_col]][ord])
  
  sobj$banksy_celltype <- factor(ct_vec)
  sobj$banksy_domain   <- factor(dm_vec)
  
  sobj@misc$banksy <- list(
    version = tryCatch(as.character(utils::packageVersion("Banksy")), error=function(e) NA),
    k_geom  = k_geom, npcs = npcs, k_neighbors = k_neighbors,
    lambda_celltype = lambda_celltype, res_celltype = res_celltype,
    lambda_domain   = lambda_domain,   res_domain   = res_domain,
    batch_col = batch_col,
    agf_used = use_agf_flag,
    chunk_size = if (has_chunk) safe_chunk else NA,
    row_limit_factor = if (has_rowlim) row_limit_factor else NA,
    n_genes = nrow(counts), n_cells = n_cells,
    spe = spe
  )
  
  n_ct <- length(levels(sobj$banksy_celltype))
  n_dm <- length(levels(sobj$banksy_domain))
  message("BANKSY 完成：banksy_celltype(", n_ct, " clusters) / banksy_domain(", n_dm, " clusters) 已写回（原子因子列）。")
  invisible(sobj)
}
# ==========================================================
# CreateCellbinOverview — Banksy/Seurat 安全版（可直接替换原函数）
# 依赖：ggplot2、dplyr、patchwork（可选）
# ==========================================================
CreateCellbinOverviewV2 <- function(
    object,
    result_path,
    cluster_col = NULL,              # 建议传 "banksy_domain" 或 "banksy_celltype"
    sample_col  = "orig.ident",      # 每个点的样本/切片列
    point_size  = 0.1,
    alpha       = 0.9,
    invert_y    = TRUE,              # Visium/HD 常需翻转y
    max_points  = 2e5,               # 大数据绘图抽样上限
    width       = 10, height = 8, dpi = 300
){
  stopifnot(inherits(object, "Seurat"))
  dir.create(result_path, showWarnings = FALSE, recursive = TRUE)
  
  suppressPackageStartupMessages({
    library(ggplot2); library(dplyr)
  })
  
  # ---------- 1) 取 meta 并做基本检查 ----------
  meta <- object@meta.data
  need_cols <- c("imagecol","imagerow")
  if (!all(need_cols %in% colnames(meta))) {
    stop("meta.data 缺少 imagecol / imagerow，请先在 Cellbin 流程中生成这些列。")
  }
  
  # 自动选择 cluster_col
  if (is.null(cluster_col)) {
    cand <- c("banksy_domain","banksy_celltype","seurat_clusters")
    cluster_col <- cand[cand %in% colnames(meta)][1]
    if (is.na(cluster_col)) stop("未找到可用的聚类列；请手动指定 cluster_col。")
    message("未指定 cluster_col，自动使用：", cluster_col)
  } else if (!cluster_col %in% colnames(meta)) {
    stop("cluster_col = '", cluster_col, "' 不在 meta.data 中。")
  }
  
  # ---------- 2) 把聚类列“原子化”为因子（杜绝 list 列） ----------
  vecify <- function(x) {
    if (is.list(x)) x <- unlist(x, use.names = FALSE)
    as.character(x)
  }
  cluster_vec <- vecify(meta[[cluster_col]])
  # 兜底：如果全 NA，报错
  if (all(is.na(cluster_vec))) stop("'", cluster_col, "' 全为 NA，无法绘图。")
  cluster_fac <- factor(cluster_vec)
  
  # 样本列（如不存在则给一个统一标签）
  if (!sample_col %in% colnames(meta)) {
    meta[[sample_col]] <- "sample"
  }
  
  # 组装绘图数据框（只保留必要列，减内存）
  df <- data.frame(
    cell     = colnames(object),
    x        = meta$imagecol,
    y        = meta$imagerow,
    cluster  = cluster_fac,
    sample   = vecify(meta[[sample_col]]),
    stringsAsFactors = FALSE
  )
  
  # ---------- 3) 汇总表 & 导出 ----------
  # 簇大小
  tab_cluster <- as.data.frame(sort(table(df$cluster), decreasing = TRUE))
  colnames(tab_cluster) <- c("cluster","n")
  tab_cluster$prop <- round(tab_cluster$n / nrow(df), 6)
  
  # 分样本 * 簇
  tab_cs <- as.data.frame.matrix(table(df$sample, df$cluster))
  write.csv(tab_cluster, file.path(result_path, "cluster_sizes.csv"), row.names = FALSE)
  write.csv(tab_cs,      file.path(result_path, "sample_by_cluster.csv"))
  
  message("簇数：", nlevels(df$cluster), "；细胞数：", nrow(df))
  
  # ---------- 4) 画图：空间散点 / 簇大小条形 / 分样本堆叠 ----------
  # 抽样（大于上限就随机抽样作图，保证速度；汇总仍用全量）
  if (nrow(df) > max_points) {
    set.seed(1)
    idx <- sample.int(nrow(df), max_points)
    dplot <- df[idx, , drop = FALSE]
    message("点数 ", nrow(df), " 较大，作图抽样 ", nrow(dplot), " 个点（不影响汇总表）。")
  } else {
    dplot <- df
  }
  
  p_sp <- ggplot(dplot, aes(x = x, y = y, color = cluster)) +
    geom_point(size = point_size, alpha = alpha) +
    coord_equal() +
    (if (invert_y) scale_y_reverse() else NULL) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = paste0("Spatial map — ", cluster_col),
         subtitle = paste0("n=", nrow(df), ", clusters=", nlevels(df$cluster))) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.position = "right")
  
  # 簇大小条形图
  p_bar <- ggplot(tab_cluster, aes(x = reorder(cluster, -n), y = n, fill = cluster)) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.3, size = 3) +
    labs(title = "Cluster sizes", x = "cluster", y = "n cells") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.position = "none")
  
  # 分样本堆叠（如果样本只有一个就跳过）
  if (length(unique(df$sample)) > 1) {
    tab_long <- df %>%
      count(sample, cluster, name = "n") %>%
      group_by(sample) %>%
      mutate(prop = n / sum(n)) %>%
      ungroup()
    p_stack <- ggplot(tab_long, aes(x = sample, y = prop, fill = cluster)) +
      geom_col() +
      labs(title = "Cluster composition by sample", x = sample_col, y = "proportion") +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    p_stack <- NULL
  }
  
  # ---------- 5) 保存图像 ----------
  ggsave(file.path(result_path, paste0("spatial_", cluster_col, ".png")),
         plot = p_sp, width = width, height = height, dpi = dpi)
  ggsave(file.path(result_path, "cluster_sizes.png"),
         plot = p_bar, width = width, height = height*0.7, dpi = dpi)
  if (!is.null(p_stack)) {
    ggsave(file.path(result_path, "sample_composition.png"),
           plot = p_stack, width = width, height = height*0.8, dpi = dpi)
  }
  
  # 可选：合并版（若你装了 patchwork）
  if (requireNamespace("patchwork", quietly = TRUE)) {
    if (is.null(p_stack)) {
      combo <- p_sp + p_bar + patchwork::plot_layout(heights = c(3,2), ncol = 1)
    } else {
      combo <- p_sp + p_bar + p_stack + patchwork::plot_layout(heights = c(3,2,2), ncol = 1)
    }
    ggsave(file.path(result_path, paste0("overview_", cluster_col, ".png")),
           plot = combo, width = width, height = height*1.6, dpi = dpi)
  }
  
  message("=== 概览完成：已输出到 ", normalizePath(result_path, winslash = "/"), " ===")
  
  invisible(list(
    cluster_sizes = tab_cluster,
    sample_by_cluster = tab_cs,
    plots = list(spatial = p_sp, sizes = p_bar, by_sample = p_stack)
  ))
}
# ======================================================================
# auto_pick_resolution_plus(): multi-criteria resolution selection
# Scores per resolution:
#   - ASW (silhouette on embeddings; subsampled)
#   - Stability (bootstrap clustering consistency; ARI if available)
#   - Spatial purity (fraction of same-cluster labels among spatial kNN)
#
# Returns: list(best_res, summary, obj)  -- compatible with your V3 caller
# ======================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

auto_pick_resolution_plus <- function(
    obj,
    # graph/embedding configs
    reduction = "pca",      # e.g., "pca" or "harmony"
    dims_use  = 1:30,
    reuse_graph = TRUE,     # reuse existing SNN if already built on the given reduction/dims
    k.param   = 20,
    prune_snn = 1/15,
    
    # resolution grid + constraints
    res_grid  = seq(0.2, 1.2, by = 0.1),
    target_k  = c(9, 30),
    min_cells_per_cluster = 200,
    
    # scoring weights
    w_asw      = 0.4,
    w_stability= 0.4,
    w_spatial  = 0.2,
    
    # ASW options
    asw_max_cells   = 4000,   # subsample size to compute silhouette
    asw_distance    = "euclidean",
    
    # Stability options
    n_boot          = 4,      # number of bootstrap replicates per resolution
    subsample_frac  = 0.8,    # fraction of cells in each bootstrap
    max_cells_stab  = 20000,  # if >, sample down before bootstrapping
    seed            = 1234,
    
    # Spatial purity options
    spatial_k       = 8,      # k for spatial neighbors
    spatial_coords  = c("imagerow","imagecol"),
    
    verbose = TRUE
){
  set.seed(seed)
  
  .msg <- function(...) if (isTRUE(verbose)) message(...)
  .nz  <- function(x) ifelse(is.finite(x), x, NA_real_)
  need_neighbors <- TRUE
  if (isTRUE(reuse_graph)) {
    # If a neighbor graph exists, reuse it (closest behavior to V3)
    need_neighbors <- length(Graphs(obj)) == 0
  }
  if (need_neighbors) {
    .msg(sprintf("Building neighbors on reduction='%s', dims=%s ...",
                 reduction, paste0(range(dims_use), collapse = ":")))
    obj <- FindNeighbors(obj, reduction = reduction, dims = dims_use,
                         k.param = k.param, prune.SNN = prune_snn, verbose = FALSE)
  }
  
  # Embedding for ASW
  emb <- Embeddings(obj, reduction = reduction)[, dims_use, drop = FALSE]
  n   <- nrow(emb)
  has_coords <- all(spatial_coords %in% colnames(obj@meta.data))
  knn_idx <- NULL
  if (has_coords) {
    xy <- as.matrix(obj@meta.data[, spatial_coords])
    if (requireNamespace("FNN", quietly = TRUE)) {
      kn <- FNN::get.knn(xy, k = spatial_k)
      knn_idx <- kn$nn.index
    } else {
      .msg("Package 'FNN' not installed; spatial purity will be NA.")
    }
  } else {
    .msg("Spatial coordinates not found in metadata; spatial purity will be NA.")
  }
  .asw_score <- function(labels){
    if (!requireNamespace("cluster", quietly = TRUE)) return(NA_real_)
    # sample cells to avoid O(n^2)
    m <- min(asw_max_cells, n)
    idx <- if (n > m) sample(n, m) else seq_len(n)
    labs <- as.integer(as.factor(labels[idx]))
    if (length(unique(labs)) < 2) return(NA_real_)
    # dist on embedding subset
    d <- dist(emb[idx, , drop = FALSE], method = asw_distance)
    sil <- tryCatch(cluster::silhouette(labs, d), error = function(e) NULL)
    if (is.null(sil)) return(NA_real_)
    mean(sil[, "sil_width"], na.rm = TRUE)
  }
  
  # Stability via bootstrap; ARI preferred
  .part_sim <- function(lab1, lab2){
    # ARI if possible
    if (requireNamespace("aricode", quietly = TRUE)) {
      return(as.numeric(aricode::ARI(lab1, lab2)))
    } else if (requireNamespace("mclust", quietly = TRUE)) {
      return(as.numeric(mclust::adjustedRandIndex(lab1, lab2)))
    }
    # Fallback: symmetric purity (0..1)
    tab <- table(lab1, lab2)
    pur1 <- sum(apply(tab, 1, max)) / sum(tab)
    pur2 <- sum(apply(tab, 2, max)) / sum(tab)
    (pur1 + pur2)/2
  }
  
  .stability_score <- function(resolution){
    # subsample universe to speed
    pool <- colnames(obj)
    if (length(pool) > max_cells_stab) pool <- sample(pool, max_cells_stab)
    reps <- vector("list", n_boot)
    for (b in seq_len(n_boot)) {
      sub_cells <- sample(pool, max(2, floor(length(pool) * subsample_frac)))
      sub <- subset(obj, cells = sub_cells)
      sub <- FindNeighbors(sub, reduction = reduction, dims = dims_use,
                           k.param = k.param, prune.SNN = prune_snn, verbose = FALSE)
      sub <- FindClusters(sub, resolution = resolution, verbose = FALSE)
      reps[[b]] <- Idents(sub)
    }
    # pairwise similarity on intersections
    sims <- c()
    for (i in 1:(n_boot-1)) {
      for (j in (i+1):n_boot) {
        common <- intersect(names(reps[[i]]), names(reps[[j]]))
        if (length(common) < 50) next
        sims <- c(sims, .part_sim(reps[[i]][common], reps[[j]][common]))
      }
    }
    if (!length(sims)) return(NA_real_)
    mean(sims, na.rm = TRUE)
  }
  
  # Spatial purity averaged over cells
  .spatial_purity <- function(labels){
    if (is.null(knn_idx)) return(NA_real_)
    labs <- as.character(labels)
    # fraction of neighbors with same label for each cell
    frac <- vapply(seq_len(nrow(knn_idx)), function(i){
      nb <- knn_idx[i, ]
      mean(labs[nb] == labs[i], na.rm = TRUE)
    }, numeric(1))
    mean(frac, na.rm = TRUE)
  }
  out <- vector("list", length(res_grid))
  names(out) <- as.character(res_grid)
  
  # ensure reproducibility across FindClusters
  set.seed(seed)
  first_feasible <- NA_real_
  
  for (r in res_grid) {
    # cluster on existing SNN (built above)
    obj <- FindClusters(obj, resolution = r, verbose = FALSE)
    labs <- Idents(obj)
    tb <- table(labs)
    k <- length(tb)
    min_sz <- min(tb)
    
    # base feasibility
    feasible <- (k >= target_k[1] && k <= target_k[2] && min_sz >= min_cells_per_cluster)
    if (isTRUE(feasible) && is.na(first_feasible)) first_feasible <- r
    
    # scores
    asw  <- .asw_score(labs)
    stab <- .stability_score(r)
    spat <- .spatial_purity(labs)
    
    out[[as.character(r)]] <- data.frame(
      resolution = r,
      k = k,
      min_size = as.numeric(min_sz),
      ASW = .nz(asw),
      Stability = .nz(stab),
      SpatialPurity = .nz(spat),
      feasible = feasible
    )
    .msg(sprintf("res=%.2f  k=%d  min=%d  ASW=%.3f  STAB=%.3f  SPAT=%.3f%s",
                 r, k, min_sz, asw %||% NA, stab %||% NA, spat %||% NA,
                 if (feasible) "  [✓ feasible]" else ""))
  }
  
  summary <- do.call(rbind, out)
  mm <- function(x){
    if (all(is.na(x))) return(rep(NA_real_, length(x)))
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) < 1e-12) return(rep(1, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }
  summary$ASW_n      <- mm(summary$ASW)
  summary$Stability_n<- mm(summary$Stability)
  summary$Spatial_n  <- mm(summary$SpatialPurity)
  
  # overall score; if all three NAs at a row, score NA
  summary$score_raw <- w_asw*summary$ASW_n + w_stability*summary$Stability_n + w_spatial*summary$Spatial_n
  
  # Penalize infeasible rows slightly to prefer feasible ones
  penalty <- ifelse(summary$feasible, 0, -0.05)
  summary$score <- summary$score_raw + penalty
  
  # -----------------------------
  # 5) Choose best resolution
  # -----------------------------
  # Prefer feasible max score; else global max
  if (any(summary$feasible, na.rm = TRUE)) {
    cand <- summary[summary$feasible %in% TRUE, ]
    best_idx <- which.max(cand$score)
    best_res <- cand$resolution[best_idx]
  } else {
    best_idx <- which.max(summary$score)
    best_res <- summary$resolution[best_idx]
  }
  
  # Re-cluster obj at best_res so caller receives that state
  obj <- FindClusters(obj, resolution = best_res, verbose = FALSE)
  
  list(best_res = best_res, summary = summary, obj = obj)
}

# little helper for printing
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a


# -------------------- 小工具：自动选择“合适”的分辨率 --------------------
auto_pick_resolution <- function(
    obj,
    dims_use = 1:30,
    res_grid = seq(0.2, 1.2, by = 0.2),
    target_k = c(12, 30),
    min_cells_per_cluster = 200,
    reduction = "pca",           # <--- NEW: can be "harmony"
    reuse_graph = TRUE,          # reuse the last SNN if already built
    k.param = NULL,              # optional: only used when recomputing
    prune_snn = NULL,            # optional: only used when recomputing
    verbose = FALSE
){
  # 1) reuse existing SNN if possible; else (re)build with requested reduction
  need_neighbors <- !reuse_graph || length(Graphs(obj)) == 0
  if (need_neighbors) {
    if (isTRUE(verbose)) message("auto_pick_resolution(): building neighbors on ", reduction)
    obj <- FindNeighbors(
      obj,
      reduction = reduction,
      dims      = dims_use,
      k.param   = if (is.null(k.param)) 20 else k.param,
      prune.SNN = if (is.null(prune_snn)) 1/15 else prune_snn,
      verbose   = FALSE
    )
  }
  
  summary <- list()
  best_res <- tail(res_grid, 1)
  
  # 2) sweep resolutions using the *current* SNN
  for (r in res_grid) {
    obj <- FindClusters(obj, resolution = r, verbose = FALSE)
    # After FindClusters, Idents(obj) already reflects the current resolution
    tb <- table(Idents(obj))
    k <- length(tb); min_sz <- min(tb)
    summary[[as.character(r)]] <- c(k = k, min_size = min_sz)
    if (k >= target_k[1] && k <= target_k[2] && min_sz >= min_cells_per_cluster) {
      best_res <- r
      break
    }
  }
  
  res_df <- do.call(rbind, summary) |> as.data.frame()
  res_df$resolution <- as.numeric(rownames(res_df))
  res_df <- res_df[, c("resolution","k","min_size")]
  
  list(best_res = best_res, summary = res_df, obj = obj)
}

# ------------------------------------------------
# MergeSimilarClusters (v4)
# - 读取 BuildClusterTree 产生的树，逐轮只测试/合并“最近的姐妹对”
# - 仅在 HVGs 上做 DE；每簇最多抽 max_cells_per_ident 个细胞
# - de_test = "presto" 可极大提速（有则用，无则回退到 wilcox）
# ------------------------------------------------
MergeSimilarClusters <- function(
    obj,
    id_col       = "seurat_clusters",
    out_col      = "coarse_clusters",
    dims_use     = 1:20,
    p_val_thresh = 0.01,
    logfc_thresh = 0.25,
    min_pct      = 0.10,
    max_iter     = 20,
    assay        = "SCT",
    de_test      = c("presto","wilcox"),
    features_use = NULL,              # 默认取 HVG（前 2000）
    max_cells_per_ident = 3000,
    prioritize_nearest  = TRUE,       # 优先合并 PCA 空间里最近的姐妹对
    verbose      = TRUE
){
  de_test <- match.arg(de_test)
  stopifnot(id_col %in% colnames(obj@meta.data))
  obj[[out_col]] <- obj[[id_col]]
  Idents(obj) <- out_col
  DefaultAssay(obj) <- assay
  
  # ---- 选基因：优先 HVG，其次全基因 ----
  if (is.null(features_use)) {
    fv <- tryCatch(VariableFeatures(obj), error = function(e) NULL)
    features_use <- if (length(fv)) head(fv, 2000) else rownames(obj)
  } else {
    features_use <- intersect(features_use, rownames(obj))
  }
  
  # ---- helpers ----
  .get_hclust <- function(tr){
    if (inherits(tr, "hclust")) return(tr)
    if (inherits(tr, "dendrogram")) return(stats::as.hclust(tr))
    if (inherits(tr, "phylo")) {
      if (!requireNamespace("ape", quietly = TRUE))
        stop("Need 'ape' to convert 'phylo' to 'hclust'.")
      return(stats::as.hclust(tr))
    }
    stop("Unknown cluster tree class: ", paste(class(tr), collapse=", "))
  }
  .sister_pairs_from_hclust <- function(hc){
    m <- hc$merge; if (is.null(m) || is.null(hc$labels)) return(list())
    idx <- which(m[,1] < 0 & m[,2] < 0)
    lapply(idx, function(i) c(hc$labels[-m[i,1]], hc$labels[-m[i,2]]))
  }
  .take <- function(v, k) if (length(v) > k) sample(v, k) else v
  
  .get_de <- function(obj, c1, c2){
    # 子集 + 等量抽样
    sub <- subset(obj, idents = c(c1, c2))
    Idents(sub) <- factor(Idents(sub))
    .take <- function(v, k) if (length(v) > k) sample(v, k) else v
    f1 <- WhichCells(sub, idents = c1)
    f2 <- WhichCells(sub, idents = c2)
    keep_cells <- c(.take(f1, max_cells_per_ident), .take(f2, max_cells_per_ident))
    sub <- subset(sub, cells = keep_cells)
    
    # 取数据：SCT@data = 基因 × 细胞
    DefaultAssay(sub) <- assay
    M <- GetAssayData(sub, assay = assay, slot = "data")
    feats <- intersect(features_use, rownames(M))
    if (length(feats) < 10L) {
      warning("Too few features after intersection; falling back to all rownames.")
      feats <- rownames(M)
    }
    M <- M[feats, , drop = FALSE]
    
    # 标签顺序必须与列顺序一致
    lab <- factor(Idents(sub)[colnames(sub)])
    
    if (de_test == "presto" && requireNamespace("presto", quietly = TRUE)) {
      # 关键点：X 必须是 基因×细胞；不要转置！
      X <- as(M, "dgCMatrix")
      w <- presto::wilcoxauc(X = X, y = lab)
      
      out <- data.frame(
        gene       = w$feature,
        avg_log2FC = w$logFC,
        p_val      = w$pval,
        stringsAsFactors = FALSE
      )
      out$p_val_adj <- p.adjust(out$p_val, method = "BH")
      out <- out[abs(out$avg_log2FC) >= logfc_thresh & !is.na(out$p_val_adj), , drop = FALSE]
      rownames(out) <- out$gene
      return(out)
    } else {
      mk <- FindMarkers(
        object = obj, ident.1 = c1, ident.2 = c2, assay = assay, test.use = "wilcox",
        features = feats, max.cells.per.ident = max_cells_per_ident,
        logfc.threshold = logfc_thresh, min.pct = min_pct,
        return.thresh = 1, verbose = FALSE
      )
      if (!"p_val_adj" %in% colnames(mk) && "p_val" %in% colnames(mk))
        mk$p_val_adj <- p.adjust(mk$p_val, method = "BH")
      return(mk)
    }
  }

  
  it <- 0
  repeat {
    it <- it + 1
    if (it > max_iter) { warning("Reached max_iter; stop merging."); break }
    
    if (isTRUE(verbose)) message(">>> Build cluster tree & sister pairs (iter ", it, ")")
    obj <- BuildClusterTree(obj, dims = dims_use, verbose = FALSE)
    tr  <- tryCatch(Tool(obj, "BuildClusterTree"), error = function(e) obj@tools$BuildClusterTree)
    if (is.null(tr)) { warning("No cluster tree; abort merging."); break }
    hc  <- .get_hclust(tr)
    sis <- .sister_pairs_from_hclust(hc)
    if (!length(sis)) { if (verbose) message("--- no sister leaf pairs; done."); break }
    
    # 优先“最近”的姐妹对
    if (isTRUE(prioritize_nearest)) {
      emb <- Embeddings(obj, "pca")[, dims_use, drop = FALSE]
      cent <- aggregate(emb, by = list(cluster = Idents(obj)), FUN = mean)
      rownames(cent) <- cent$cluster; cent$cluster <- NULL
      dist_of <- function(a,b){
        if (!(a %in% rownames(cent) && b %in% rownames(cent))) return(Inf)
        pa <- as.numeric(cent[a,,drop=FALSE]); pb <- as.numeric(cent[b,,drop=FALSE])
        sqrt(sum((pa - pb)^2))
      }
      ord <- order(sapply(sis, function(p) dist_of(p[1], p[2])))
      sis <- sis[ord]
    }
    
    merged_once <- FALSE
    for (pair in sis) {
      c1 <- pair[1]; c2 <- pair[2]
      labs_now <- levels(obj[[out_col]][,1])
      if (!(c1 %in% labs_now && c2 %in% labs_now)) next
      if (isTRUE(verbose)) message("    - testing ", c1, " vs ", c2)
      
      mk <- .get_de(obj, c1, c2)
      if (!"p_val_adj" %in% colnames(mk)) mk$p_val_adj <- 1
      sig <- tryCatch(subset(mk, p_val_adj < p_val_thresh), error = function(e) mk[0, , drop = FALSE])
      
      if (nrow(sig) == 0) {
        if (isTRUE(verbose)) message("    >>> merge: ", c1, " ← ", c2, " (no sig DE)")
        labs <- as.character(obj[[out_col]][,1])
        nums <- suppressWarnings(as.numeric(c(c1, c2)))
        if (all(!is.na(nums))) { keep <- as.character(min(nums)); drop <- as.character(max(nums)) }
        else { keep <- min(c(c1, c2)); drop <- max(c(c1, c2)) }
        labs[labs == drop] <- keep
        obj[[out_col]] <- droplevels(factor(labs))
        Idents(obj) <- out_col
        merged_once <- TRUE
        break
      }
    }
    if (!merged_once) { if (verbose) message("--- no merge this round; done."); break }
  }
  if (isTRUE(verbose)) message(">>> merging done. final K = ", length(levels(obj[[out_col]][,1])))
  return(obj)
}


CellbinClusterCore <- function(
    obj,
    # I/O
    result_path = NULL,
    sample_id   = NULL,
    # Normalization
    use_sct = TRUE,                        # TRUE=SCT, FALSE=LogNormalize
    assay_in = c("Spatial","RNA"),
    regress_vars = c("log10_nCount_Spatial","area_um2"), # set NULL for conservative
    sct_features_n     = 3000,
    lognorm_features_n = 2000,
    # Dimensionality & graph
    npcs = 30,
    dims_use = 1:30,
    k_param = 30,
    # Resolution sweep
    res_grid = seq(0.2, 1.2, by = 0.2),
    target_k = c(12, 30),
    min_cells_per_cluster = 200,
    id_col_name = "seurat_clusters",
    # Optional DE-based merging (coarse)
    use_merge = FALSE,
    merge_out_col    = "coarse_clusters",
    merge_p_val_thresh = 0.05,
    merge_logfc_thresh = 0.25,
    merge_min_pct      = 0.1,
    # NEW: outputs
    make_plots = TRUE,
    export_loupe = TRUE,
    write_topN = 50,                  # 每簇TopN(上调)
    write_coords = FALSE,             # 是否导出坐标/UMAP到TSV
    de_test = "wilcox",               # FindAllMarkers test
    de_min_pct = 0.05,
    de_logfc   = 0.15,
    # plotting
    plot_bg = "white",
    umap_w = 7, umap_h = 6,
    spatial_w = 8, spatial_h = 7,
    # misc
    seed = 1234,
    verbose = TRUE
){
  stopifnot(inherits(obj, "Seurat"))
  set.seed(seed)
  if (!is.null(result_path) && !dir.exists(result_path)) dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
  sid <- if (is.null(sample_id)) "Sample" else sample_id
  
  # —— 新增：确保依赖函数已加载（你已说会单独加载，因此只做断言）——
  if (!exists("auto_pick_resolution", mode = "function"))
    stop("auto_pick_resolution() 未加载：请先 source 进来。")
  if (isTRUE(use_merge) && ! (exists("MergeSimilarClusters_v3", mode="function") ||
                              exists("MergeSimilarClusters",    mode="function")))
    stop("请先加载 MergeSimilarClusters_v3() 或 MergeSimilarClusters()。")
  
  # ---------- prepare regressors ----------
  if ("nCount_Spatial" %in% colnames(obj@meta.data) &&
      !"log10_nCount_Spatial" %in% colnames(obj@meta.data)) {
    obj$log10_nCount_Spatial <- log10(obj$nCount_Spatial + 1)
  }
  vars_reg <- if (is.null(regress_vars)) NULL else intersect(regress_vars, colnames(obj@meta.data))
  if (verbose && length(vars_reg)) message("Regress vars: ", paste(vars_reg, collapse=", "))
  
  # ---------- normalization ----------
  assay_in <- match.arg(assay_in)
  if (!assay_in %in% Assays(obj)) {
    assay_in <- if ("Spatial" %in% Assays(obj)) "Spatial" else {
      if ("RNA" %in% Assays(obj)) "RNA" else stop("No 'Spatial' or 'RNA' assay in object.")
    }
  }
  if (use_sct) {
    DefaultAssay(obj) <- assay_in
    obj <- SCTransform(
      obj, assay = assay_in, vst.flavor = "v2",
      method = if (requireNamespace("glmGamPoi", quietly=TRUE)) "glmGamPoi" else "poisson",
      vars.to.regress = vars_reg,
      variable.features.n = sct_features_n,
      return.only.var.genes = TRUE, verbose = FALSE
    )
    DefaultAssay(obj) <- "SCT"
    assay_norm <- "SCT"
  } else {
    DefaultAssay(obj) <- assay_in
    tmp_assay <- paste0(assay_in, "_log")
    obj[[tmp_assay]] <- CreateAssayObject(
      counts = GetAssayData(obj, assay = assay_in, slot = "counts")
    )
    DefaultAssay(obj) <- tmp_assay
    obj <- NormalizeData(obj, normalization.method="LogNormalize", scale.factor=1e4, verbose=FALSE)
    obj <- FindVariableFeatures(obj, selection.method="vst", nfeatures=lognorm_features_n, verbose=FALSE)
    obj <- ScaleData(obj, features=VariableFeatures(obj), verbose=FALSE)
    assay_norm <- tmp_assay
  }
  
  # ---------- PCA / UMAP / neighbors ----------
  obj <- RunPCA(obj, npcs = max(npcs, max(dims_use)), features = VariableFeatures(obj), verbose=FALSE)
  obj <- RunUMAP(obj, dims = dims_use, verbose=FALSE)
  obj <- FindNeighbors(obj, dims = dims_use, k.param = k_param, verbose=FALSE)
  
  # ---------- auto-pick resolution ----------
  pick <- auto_pick_resolution(obj, dims_use=dims_use, res_grid=res_grid,
                               target_k=target_k, min_cells_per_cluster=min_cells_per_cluster)
  obj <- pick$obj
  best_res <- pick$best_res
  if (verbose) message(">>> best resolution: ", best_res)
  if (!is.null(result_path) && !is.null(sample_id)) {
    readr::write_tsv(pick$summary, file.path(result_path, paste0(sid, "_resolution_sweep.tsv")))
  }
  
  # ---------- finalize cluster column ----------
  obj <- FindClusters(obj, resolution=best_res, verbose=FALSE)
  Idents(obj) <- "seurat_clusters"
  obj[[id_col_name]] <- obj$seurat_clusters
  
  # ---------- optional merging ----------
  if (isTRUE(use_merge)) {
    obj <- MergeSimilarClusters(
      obj,
      id_col  = id_col_name,
      out_col = merge_out_col,
      dims_use = dims_use,
      p_val_thresh = merge_p_val_thresh,
      logfc_thresh = merge_logfc_thresh,
      min_pct = merge_min_pct,
      verbose = verbose
    )
  }
  
  # ---------- choose final cluster column ----------
  final_id_col <- if (isTRUE(use_merge)) merge_out_col else id_col_name
  Idents(obj) <- final_id_col
  
  # ---------- outputs: tables & plots ----------
  if (!is.null(result_path) && !is.null(sample_id)) {
    # cluster sizes
    clu_tab <- as.data.frame(table(obj[[final_id_col]][,1])) |>
      dplyr::arrange(desc(Freq))
    readr::write_tsv(clu_tab, file.path(result_path, paste0(sid, "_ClusterSizes_", final_id_col, ".tsv")))
    
    # TopN per cluster (one-vs-all)
    markers_final <- FindAllMarkers(
      obj, assay = assay_norm,
      only.pos = TRUE, test.use = de_test,
      min.pct = de_min_pct, logfc.threshold = de_logfc
    )
    saveRDS(markers_final, file = file.path(result_path, paste0(sid, "_", final_id_col, "_markers.rds")))
    if (!is.null(write_topN) && write_topN > 0) {
      topN <- markers_final |>
        dplyr::group_by(cluster) |>
        dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) |>
        dplyr::slice_head(n = write_topN) |>
        dplyr::ungroup()
      readr::write_tsv(topN, file.path(result_path, paste0(sid, "_", final_id_col, "_Top", write_topN, ".tsv")))
      if (!requireNamespace("writexl", quietly = TRUE)) suppressWarnings(try(install.packages("writexl"), silent=TRUE))
      if (requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(topN, path = file.path(result_path, paste0(sid, "_", final_id_col, "_Top", write_topN, ".xlsx")))
      }
    }
    
    # Plots
    if (isTRUE(make_plots)) {
      p_umap <- DimPlot(obj, reduction="umap", group.by=final_id_col, label=TRUE) + NoLegend()
      ggsave(file.path(result_path, paste0(sid, "_UMAP_", final_id_col, ".png")),
             p_umap, width = umap_w, height = umap_h, dpi = 300, bg = plot_bg)
      
      if (exists("PlotSpatialDistribution")) {
        p_sp <- PlotSpatialDistribution(obj, group_by = final_id_col)
        ggsave(file.path(result_path, paste0(sid, "_Spatial_", final_id_col, ".png")),
               p_sp, width = spatial_w, height = spatial_h, dpi = 300, bg = plot_bg)
      }
    }
    
    # Loupe categories (barcode -> final cluster)
    if (isTRUE(export_loupe) && exists("export_loupe_categories")) {
      export_loupe_categories(obj, meta_col = final_id_col, out_dir = result_path)
    }
    
    # Optional: coordinates/UMAP TSV
    if (isTRUE(write_coords)) {
      coords_umap <- data.frame(
        barcode = colnames(obj),
        UMAP_1  = Embeddings(obj, "umap")[,1],
        UMAP_2  = Embeddings(obj, "umap")[,2],
        imagerow = obj$imagerow,
        imagecol = obj$imagecol,
        cluster  = obj[[final_id_col]][,1]
      )
      readr::write_tsv(coords_umap, file.path(result_path, paste0(sid, "_", final_id_col, "_coords_umap.tsv")))
    }
  }
  
  return(obj)
}

# ===========================================================
# CellbinClusterCore_V3 (patched, cleanup default = FALSE)
# - Regressors supported (LogNorm: ScaleData vars.to.regress; SCT: SCTransform)
# - Seurat v5 friendly (layers with fallback to slots)
# - Robust auto_pick_resolution() call (detect arg names; pass only supported)
# - Robust MergeSimilarClusters() call (detect arg names; pass only supported)
# - Optional radius-based spatial cleanup with CHANGED RATE report & TSV
# ===========================================================
CellbinClusterCore_V3 <- function(
    obj,
    # I/O
    result_path = NULL,
    sample_id   = NULL,
    
    # Normalization & Scaling (Loupe-like defaults)
    use_sct    = FALSE,
    assay_in   = c("Spatial", "RNA"),
    
    # Technical covariates to regress (metadata columns)
    # "nCount_Spatial" will be mapped to "log10_nCount_Spatial" internally
    regress_vars = c("nCount_Spatial"),
    
    lognorm_features_n = 1500,
    exclude_hvgs_regex = "^(MT-|RPS|RPL|HBA|HBB)",
    
    # Dimensionality Reduction & Optional Harmony
    npcs            = 30,
    use_harmony     = F,                # set FALSE for single-slice coarse clustering
    harmony_group_by= "orig.ident",
    harmony_theta   = 2,
    
    # Graph & Clustering
    dims_use  = 1:30,                      # will become 1:npcs if Harmony is used
    k_param   = 20,
    prune_snn = 1/15,
    
    # Resolution sweep
    res_grid = seq(0.2, 1.2, by = 0.1),
    target_k = c(9, 30),
    min_cells_per_cluster = 200,
    id_col_name = "seurat_clusters",
    
    # Optional merging (coarse)
    use_merge        = F,
    merge_out_col    = "coarse_clusters",
    merge_p_val_thresh = 0.001,
    merge_logfc_thresh = 0.5,
    merge_min_pct      = 0.25,   # passed if MergeSimilarClusters() supports 'min_pct'
    
    # Spatial cleanup (radius-based)  —— 默认关闭
    spatial_majority_cleanup = FALSE,
    cleanup_iters  = 1,
    cleanup_radius = 9,                 # pixels/bins
    cleanup_dist_metric = "euclidean",  # "euclidean"|"manhattan"
    
    # Outputs
    make_plots   = TRUE,
    export_loupe = TRUE,
    write_topN   = 50,
    write_coords = TRUE,
    de_test      = "wilcox",
    de_min_pct   = 0.05,
    de_logfc     = 0.15,
    
    # plotting & misc
    plot_bg = "white",
    umap_w = 7, umap_h = 6,
    spatial_w = 8, spatial_h = 7,
    seed = 1234,
    verbose = TRUE
){
  stopifnot(inherits(obj, "Seurat"))
  set.seed(seed)
  if (!is.null(result_path) && !dir.exists(result_path)) dir.create(result_path, recursive = TRUE, showWarnings = FALSE)
  sid <- if (is.null(sample_id)) "Sample" else sample_id
  
  # ---- helper: get counts layer/slot safely (v5/v4) ----
  .get_counts <- function(object, assay){
    lyr <- tryCatch({ if ("counts" %in% Layers(object[[assay]])) "counts" else "data" },
                    error = function(e) "counts")
    tryCatch(GetAssayData(object, assay = assay, layer = lyr),
             error = function(e) GetAssayData(object, assay = assay, slot = "counts"))
  }
  
  # ---- dependencies ----
  if (isTRUE(use_harmony) && !requireNamespace("harmony", quietly = TRUE)) {
    stop("'use_harmony' is TRUE, but package 'harmony' is not installed.")
  }
  if (!exists("auto_pick_resolution", mode = "function"))
    stop("Helper function 'auto_pick_resolution' is not loaded. Please source it.")
  has_merge_fun <- exists("MergeSimilarClusters", mode = "function")
  
  # ---- prepare regressors ----
  if ("nCount_Spatial" %in% colnames(obj@meta.data) &&
      !"log10_nCount_Spatial" %in% colnames(obj@meta.data)) {
    obj$log10_nCount_Spatial <- log10(obj$nCount_Spatial + 1)
  }
  vars_reg <- regress_vars
  if (!is.null(vars_reg) && length(vars_reg)) {
    vars_reg[vars_reg == "nCount_Spatial"] <- "log10_nCount_Spatial"
  }
  if (!is.null(vars_reg) && "percent.mt" %in% vars_reg && !"percent.mt" %in% colnames(obj@meta.data)) {
    assay_guess <- if ("Spatial" %in% Assays(obj)) "Spatial" else DefaultAssay(obj)
    cts <- .get_counts(obj, assay_guess)
    mt_idx <- grep("^MT-", rownames(cts))
    total  <- Matrix::colSums(cts)
    obj$percent.mt <- if (length(mt_idx)) 100 * Matrix::colSums(cts[mt_idx, , drop = FALSE]) / pmax(1, total) else 0
    if (verbose) message("Computed 'percent.mt' from counts via '^MT-' pattern.")
  }
  vars_reg <- intersect(vars_reg, colnames(obj@meta.data))
  if (verbose) {
    if (length(vars_reg)) message("Variables to regress out: ", paste(vars_reg, collapse = ", "))
    else message("No valid regressors found in metadata. Proceeding without regression.")
  }
  
  # ---- choose assay ----
  assay_in <- match.arg(assay_in)
  DefaultAssay(obj) <- assay_in
  
  # ---- normalization & HVG ----
  if (use_sct) {
    if (verbose) message("Step 1: SCTransform workflow...")
    obj <- SCTransform(
      obj, assay = assay_in, vst.flavor = "v2",
      method = if (requireNamespace("glmGamPoi", quietly=TRUE)) "glmGamPoi" else "poisson",
      vars.to.regress = vars_reg,
      variable.features.n = 3000,
      return.only.var.genes = TRUE, verbose = FALSE
    )
    DefaultAssay(obj) <- "SCT"
    assay_norm <- "SCT"
  } else {
    if (verbose) message("Step 1: LogNormalize + HVG + ScaleData (with regression if specified)...")
    tmp_assay <- paste0(assay_in, "_log")
    obj[[tmp_assay]] <- CreateAssayObject(counts = .get_counts(obj, assay_in))
    DefaultAssay(obj) <- tmp_assay
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = lognorm_features_n, verbose = FALSE)
    if (!is.null(exclude_hvgs_regex)) {
      hv <- VariableFeatures(obj)
      bad <- grep(exclude_hvgs_regex, hv, value = TRUE)
      VariableFeatures(obj) <- setdiff(hv, bad)
      if (verbose) message("  - HVGs after exclusion: ", length(VariableFeatures(obj)))
    }
    obj <- ScaleData(obj, features = VariableFeatures(obj), vars.to.regress = vars_reg, verbose = FALSE)
    assay_norm <- tmp_assay
  }
  
  # ---- PCA & optional Harmony ----
  if (verbose) message("Step 2: PCA and optional Harmony correction...")
  obj <- RunPCA(obj, npcs = npcs, features = VariableFeatures(obj), verbose = FALSE)
  
  reduction_to_use <- "pca"
  dims_to_use_clustering <- dims_use
  
  if (isTRUE(use_harmony)) {
    if (!harmony_group_by %in% colnames(obj@meta.data)) {
      warning("`harmony_group_by` not found; falling back to 'orig.ident'.")
      if (!"orig.ident" %in% colnames(obj@meta.data)) obj$orig.ident <- sample_id
      harmony_group_by <- "orig.ident"
    }
    if (verbose) message("  - Applying Harmony (grouping by '", harmony_group_by, "')...")
    obj <- harmony::RunHarmony(
      obj,
      group.by.vars = harmony_group_by,
      reduction = "pca",
      assay.use = assay_norm,
      reduction.save = "harmony",
      theta = harmony_theta,
      verbose = FALSE
    )
    reduction_to_use <- "harmony"
    dims_to_use_clustering <- 1:npcs
  }
  
  # ---- UMAP & neighbors ----
  if (verbose) message("Step 3: UMAP and FindNeighbors (", reduction_to_use, ")...")
  obj <- RunUMAP(obj, dims = dims_to_use_clustering, reduction = reduction_to_use, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims_to_use_clustering, reduction = reduction_to_use,
                       k.param = k_param, prune.SNN = prune_snn, verbose = FALSE)
  
  # ---- auto-pick resolution (robust arg matching) ----
  if (verbose) message("Step 4: Auto-picking clustering resolution...")
  fm <- tryCatch(names(formals(auto_pick_resolution)), error = function(e) character())
  args <- list(obj = obj)
  if ("dims_use" %in% fm)       args$dims_use <- dims_to_use_clustering
  else if ("dims" %in% fm)      args$dims     <- dims_to_use_clustering
  else if ("pc_dims" %in% fm)   args$pc_dims  <- dims_to_use_clustering
  if ("reduction" %in% fm)      args$reduction <- reduction_to_use
  if ("res_grid" %in% fm)       args$res_grid <- res_grid
  if ("target_k" %in% fm)       args$target_k <- target_k
  if ("min_cells_per_cluster" %in% fm) args$min_cells_per_cluster <- min_cells_per_cluster
  pick <- do.call(auto_pick_resolution, args)
  
  obj <- pick$obj
  best_res <- pick$best_res
  if (verbose) message("  - Best resolution: ", best_res)
  
  obj <- FindClusters(obj, resolution = best_res, verbose = FALSE)
  obj[[id_col_name]] <- obj$seurat_clusters
  Idents(obj) <- id_col_name
  
  # ---- optional merging (robust arg matching) ----
  final_id_col <- id_col_name
  if (isTRUE(use_merge)) {
    if (has_merge_fun) {
      if (verbose) message("Step 5: Merging similar clusters...")
      fm2 <- tryCatch(names(formals(MergeSimilarClusters)), error = function(e) character())
      m_args <- list(
        obj          = obj,
        id_col       = id_col_name,
        out_col      = merge_out_col,
        dims_use     = dims_to_use_clustering,
        p_val_thresh = merge_p_val_thresh,
        logfc_thresh = merge_logfc_thresh
      )
      if ("reduction" %in% fm2) m_args$reduction <- reduction_to_use
      if ("min_pct"   %in% fm2) m_args$min_pct   <- merge_min_pct
      if ("verbose"   %in% fm2) m_args$verbose   <- verbose
      obj <- do.call(MergeSimilarClusters, m_args)
      final_id_col <- merge_out_col
    } else {
      warning("MergeSimilarClusters() not found; skipping merging.")
    }
  }
  Idents(obj) <- final_id_col
  
  # ---- spatial majority cleanup (radius-based, with changed-rate) ----
  if (isTRUE(spatial_majority_cleanup) && all(c("imagerow","imagecol") %in% colnames(obj@meta.data))) {
    if (!is.null(cleanup_radius) && cleanup_radius > 0) {
      if (verbose) message("Step 6: Spatial cleanup by radius (", cleanup_dist_metric,
                           "): r=", cleanup_radius, ", iters=", cleanup_iters, "...")
      xy <- obj@meta.data[, c("imagerow", "imagecol")]
      lab_before <- as.character(obj[[final_id_col]][, 1])
      lab <- lab_before
      
      for (it in seq_len(cleanup_iters)) {
        if (verbose) message("  - Cleanup iteration ", it, "...")
        new_labels <- sapply(seq_len(nrow(xy)), function(i) {
          # bounding-box prefilter
          cand <- which(
            abs(xy$imagerow - xy$imagerow[i]) <= cleanup_radius &
              abs(xy$imagecol - xy$imagecol[i]) <= cleanup_radius
          )
          # exact neighborhood
          if (cleanup_dist_metric == "euclidean") {
            d2 <- (xy$imagerow[cand] - xy$imagerow[i])^2 + (xy$imagecol[cand] - xy$imagecol[i])^2
            nbrs <- cand[d2 <= cleanup_radius^2 & d2 > 1e-9]
          } else { # manhattan
            nbrs <- setdiff(cand, i)
          }
          # majority vote
          if (length(nbrs) > 0) {
            tb <- table(lab[nbrs]); names(tb)[which.max(tb)]
          } else {
            lab[i]
          }
        })
        lab <- new_labels
      }
      
      cleaned_col <- paste0(final_id_col, "_cleaned")
      obj[[cleaned_col]] <- factor(lab, levels = sort(unique(lab)))
      Idents(obj) <- cleaned_col
      final_id_col <- cleaned_col
      
      # changed rate
      changed_idx  <- lab_before != lab
      changed_n    <- sum(changed_idx, na.rm = TRUE)
      total_n      <- length(lab_before)
      changed_rate <- if (total_n > 0) changed_n / total_n else NA_real_
      if (verbose) message(sprintf("  - Cleanup changed %d/%d (%.2f%%) labels.",
                                   changed_n, total_n, 100 * changed_rate))
      
      # attach metrics & optional TSV
      attr(obj, "spatial_cleanup_metrics") <- list(
        changed_n    = changed_n,
        total_n      = total_n,
        changed_rate = changed_rate,
        radius       = cleanup_radius,
        iters        = cleanup_iters,
        dist_metric  = cleanup_dist_metric
      )
      if (!is.null(result_path) && !is.null(sample_id)) {
        metadf <- data.frame(
          sample_id    = sample_id,
          radius       = cleanup_radius,
          iters        = cleanup_iters,
          dist_metric  = cleanup_dist_metric,
          changed_n    = changed_n,
          total_n      = total_n,
          changed_rate = changed_rate
        )
        readr::write_tsv(metadf, file.path(result_path, paste0(sid, "_SpatialCleanupMetrics.tsv")))
      }
      if (verbose) message("  - Cleanup finished. Final clusters in '", final_id_col, "'.")
    }
  }
  
  # ---- outputs ----
  if (!is.null(result_path) && !is.null(sample_id)) {
    # cluster sizes
    clu_tab <- as.data.frame(table(obj[[final_id_col]][,1]))
    clu_tab <- clu_tab[order(-clu_tab$Freq), ]
    readr::write_tsv(clu_tab, file.path(result_path, paste0(sid, "_ClusterSizes_", final_id_col, ".tsv")))
    
    # markers
    markers_final <- FindAllMarkers(
      obj, assay = assay_norm,
      only.pos = TRUE, test.use = de_test,
      min.pct = de_min_pct, logfc.threshold = de_logfc
    )
    saveRDS(markers_final, file = file.path(result_path, paste0(sid, "_", final_id_col, "_markers.rds")))
    
    if (!is.null(write_topN) && write_topN > 0) {
      topN <- markers_final |>
        dplyr::group_by(cluster) |>
        dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) |>
        dplyr::slice_head(n = write_topN) |>
        dplyr::ungroup()
      readr::write_tsv(topN, file.path(result_path, paste0(sid, "_", final_id_col, "_Top", write_topN, ".tsv")))
      if (!requireNamespace("writexl", quietly = TRUE)) suppressWarnings(try(install.packages("writexl"), silent = TRUE))
      if (requireNamespace("writexl", quietly = TRUE)) {
        writexl::write_xlsx(topN, path = file.path(result_path, paste0(sid, "_", final_id_col, "_Top", write_topN, ".xlsx")))
      }
    }
    
    # plots
    if (isTRUE(make_plots)) {
      p_umap <- DimPlot(obj, reduction = "umap", group.by = final_id_col, label = TRUE) + NoLegend()
      ggsave(file.path(result_path, paste0(sid, "_UMAP_", final_id_col, ".png")),
             p_umap, width = umap_w, height = umap_h, dpi = 300, bg = plot_bg)
      
      if (exists("PlotSpatialDistribution")) {
        p_sp <- PlotSpatialDistribution(obj, group_by = final_id_col)
      } else if ("Spatial" %in% names(obj@images)) {
        p_sp <- SpatialDimPlot(obj, group.by = final_id_col, label = TRUE)
      } else {
        p_sp <- NULL
      }
      if (!is.null(p_sp)) {
        ggsave(file.path(result_path, paste0(sid, "_Spatial_", final_id_col, ".png")),
               p_sp, width = spatial_w, height = spatial_h, dpi = 300, bg = plot_bg)
      }
    }
    
    # Loupe categories
    if (isTRUE(export_loupe) && exists("export_loupe_categories")) {
      export_loupe_categories(obj, meta_col = final_id_col, out_dir = result_path)
    }
    
    # coords/UMAP TSV
    if (isTRUE(write_coords)) {
      coords_umap <- data.frame(
        barcode = colnames(obj),
        UMAP_1  = Embeddings(obj, "umap")[,1],
        UMAP_2  = Embeddings(obj, "umap")[,2],
        imagerow = obj$imagerow,
        imagecol = obj$imagecol,
        cluster  = obj[[final_id_col]][,1]
      )
      readr::write_tsv(coords_umap, file.path(result_path, paste0(sid, "_", final_id_col, "_coords_umap.tsv")))
    }
  }
  
  return(obj)
}
# =====================================================================
# cellbin_preflight_report(): Preflight diagnostics + V3-ready hints
# Embedded utilities:
#   - PC loadings (dominant genes per PC; pos/neg)
#   - PC ~ cluster separation (eta^2, ANOVA F, p/padj)
#   - PC cluster centroids (cluster means on each PC)
#
# Notes:
# - Safe: works on a temporary copy; original object is not modified.
# - Seurat v4/v5 compatible (uses layers when available; falls back to slots).
# - English-only comments for team sharing & reproducibility.
# =====================================================================
cellbin_preflight_report <- function(
    obj,
    # --- general controls ---
    assay_priority      = c("Spatial","RNA"),
    tmp_hvgs            = 1500,
    tmp_npcs            = 30,
    corr_rho_thresh     = 0.35,      # threshold to flag technical dominance
    spatial_knn_k       = 8,         # for estimating spatial scale (kNN), kept for API
    outdir              = NULL,      # if not NULL, save TSVs into this directory
    prefix              = "Preflight",
    verbose             = TRUE,
    
    # --- embedded PC loadings ---
    do_pc_loadings      = TRUE,
    pc_loadings_pcs     = 1:20,
    pc_loadings_top_n   = 30,
    
    # --- embedded PC~cluster metrics ---
    # If cluster labels are not available (NULL), these two will be skipped.
    cluster_id_col      = NULL,      # e.g., "coarse_clusters" or Idents(obj) if NULL
    do_cluster_metrics  = TRUE,      # compute separation & centroids if labels exist
    pc_cluster_pcs      = 1:20
){
  stopifnot(inherits(obj, "Seurat"))
  
  # ----------------------------
  # Helper: quiet message
  .msg <- function(...) if (isTRUE(verbose)) message(...)
  
  # ----------------------------
  # Helper: write TSV if outdir is provided
  .write_tsv <- function(df, name){
    if (is.null(outdir)) return(invisible(NULL))
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    fn <- file.path(outdir, paste0(prefix, "_", name, ".tsv"))
    readr::write_tsv(df, fn)
    .msg("Saved: ", fn)
    invisible(fn)
  }
  
  # ----------------------------
  # Helper: get counts layer/slot safely (Seurat v5/v4 compatible)
  .get_counts <- function(object, assay){
    layer <- tryCatch({
      if ("counts" %in% Layers(object[[assay]])) "counts" else "data"
    }, error = function(e) "counts")
    x <- tryCatch(
      GetAssayData(object, assay = assay, layer = layer),
      error = function(e) GetAssayData(object, assay = assay, slot = "counts")
    )
    x
  }
  
  # ----------------------------
  # Helper: PC loadings table (balanced pos/neg by top |loading|)
  .pc_loading_table <- function(object, reduction = "pca", pcs = 1:20, top_n = 30){
    loadings <- tryCatch(
      Loadings(object, reduction = reduction),
      error = function(e) object@reductions[[reduction]]@feature.loadings
    )
    pcs <- pcs[pcs <= ncol(loadings)]
    tab_list <- lapply(pcs, function(i){
      v <- loadings[, i]
      ord <- order(abs(v), decreasing = TRUE)
      genes <- rownames(loadings)[ord]
      vals  <- v[ord]
      df <- data.frame(
        PC = paste0("PC_", i),
        gene = genes,
        loading = vals,
        abs_loading = abs(vals),
        direction = ifelse(vals >= 0, "pos", "neg"),
        rank = seq_along(vals),
        stringsAsFactors = FALSE
      )
      keep <- unlist(tapply(seq_len(nrow(df)), df$direction, function(ii) head(ii, top_n)))
      df[sort(keep), ]
    })
    do.call(rbind, tab_list)
  }
  
  # ----------------------------
  # Helper: which PCs separate clusters best (ANOVA-like heuristic)
  .pc_cluster_separation <- function(object, labels, reduction = "pca", pcs = 1:20){
    emb <- Embeddings(object, reduction = reduction)
    pcs <- pcs[pcs <= ncol(emb)]
    lab <- as.factor(labels)
    res <- lapply(pcs, function(i){
      x <- emb[, i]
      groups <- split(x, lab)
      n <- length(x); k <- length(groups)
      m_all <- mean(x)
      bss <- sum(sapply(groups, function(g) length(g) * (mean(g) - m_all)^2))
      wss <- sum(sapply(groups, function(g) sum((g - mean(g))^2)))
      eta2 <- bss / (bss + wss + 1e-12)
      Fval <- (bss / max(1, (k - 1))) / (wss / max(1, n - k))
      pval <- stats::pf(Fval, df1 = max(1, k - 1), df2 = max(1, n - k), lower.tail = FALSE)
      data.frame(PC = paste0("PC_", i), eta2 = eta2, F = Fval, p = pval, stringsAsFactors = FALSE)
    })
    df <- do.call(rbind, res)
    df$padj <- p.adjust(df$p, method = "BH")
    df[order(-df$eta2), ]
  }
  
  # ----------------------------
  # Helper: cluster centroids on selected PCs
  .pc_cluster_centroids <- function(object, labels, reduction = "pca", pcs = 1:10){
    lab <- as.factor(labels)
    emb <- Embeddings(object, reduction = reduction)[, pcs, drop = FALSE]
    means <- aggregate(emb, by = list(cluster = lab), FUN = mean)
    means
  }
  
  # ============================
  # 1) Choose assay
  # ============================
  assay_in <- intersect(assay_priority, Assays(obj))
  assay_in <- if (length(assay_in)) assay_in[1] else DefaultAssay(obj)
  DefaultAssay(obj) <- assay_in
  
  # ============================
  # 2) Metadata & coordinates
  # ============================
  has_coords <- all(c("imagerow","imagecol") %in% colnames(obj@meta.data))
  n_cells <- ncol(obj)
  n_feats <- nrow(obj[[assay_in]])
  imgs <- tryCatch(names(obj@images), error = function(e) character(0))
  
  # Derive log10_nCount_Spatial if missing
  if ("nCount_Spatial" %in% colnames(obj@meta.data) &&
      !"log10_nCount_Spatial" %in% colnames(obj@meta.data)) {
    obj$log10_nCount_Spatial <- log10(obj$nCount_Spatial + 1)
  }
  # Derive percent.mt if missing (best-effort)
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    assay_guess <- if ("Spatial" %in% Assays(obj)) "Spatial" else DefaultAssay(obj)
    cts <- .get_counts(obj, assay_guess)
    if (is(cts, "dgCMatrix")) {
      mt_idx <- grep("^MT-", rownames(cts))
      total  <- Matrix::colSums(cts)
      obj$percent.mt <- if (length(mt_idx)) 100 * Matrix::colSums(cts[mt_idx, , drop = FALSE]) / pmax(1, total) else 0
    }
  }
  
  # ============================
  # 3) Temporary PCA pipeline (lightweight; LogNormalize only)
  # ============================
  obj_tmp <- obj
  DefaultAssay(obj_tmp) <- assay_in
  tmp_assay <- paste0(assay_in, "_preflight")
  obj_tmp[[tmp_assay]] <- CreateAssayObject(counts = .get_counts(obj_tmp, assay_in))
  DefaultAssay(obj_tmp) <- tmp_assay
  obj_tmp <- NormalizeData(obj_tmp, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  obj_tmp <- FindVariableFeatures(obj_tmp, selection.method = "vst", nfeatures = tmp_hvgs, verbose = FALSE)
  obj_tmp <- ScaleData(obj_tmp, features = VariableFeatures(obj_tmp), verbose = FALSE)
  obj_tmp <- RunPCA(obj_tmp, npcs = tmp_npcs, features = VariableFeatures(obj_tmp), verbose = FALSE)
  
  # ============================
  # 4) PC ~ technical correlations
  # ============================
  have_covs <- intersect(
    c("log10_nCount_Spatial","percent.mt","nCount_Spatial","area_um2"),  # ← include area_um2
    colnames(obj_tmp@meta.data)
  )
  pc_mat <- Embeddings(obj_tmp, "pca")[, 1:tmp_npcs, drop = FALSE]
  
  # Variance explained per PC
  pc_var <- tryCatch(Stdev(obj_tmp, "pca")^2,
                     error = function(e) apply(pc_mat, 2, stats::var))
  var_explained <- pc_var / sum(pc_var)
  
  pc_cov_corr <- NULL
  if (length(have_covs)) {
    pc_cov_corr <- sapply(have_covs, function(v) {
      x <- obj_tmp[[v]][,1]
      apply(pc_mat, 2, function(pc) suppressWarnings(cor(pc, x, method = "spearman")))
    })
  }
  
  # ============================
  # 5) Heuristic suggestions
  # ============================
  cumvar <- cumsum(var_explained)
  dims_suggest <- min(max(12, which(cumvar >= 0.85)[1]), 30)
  if (is.na(dims_suggest)) dims_suggest <- 20
  
  # k.param: scale with n; clip to [15, 40]
  k_param_suggest <- max(15, min(40, round(0.0015 * n_cells)))
  prune_snn_suggest <- 1/15
  
  # cleanup radius via spatial NN scale (median 1-NN distance -> radius ~ [5,10])
  cleanup_radius_suggest <- NA_integer_
  spatial_scale <- NA_real_
  if (has_coords) {
    xy <- as.matrix(obj@meta.data[, c("imagerow","imagecol")])
    if (requireNamespace("FNN", quietly = TRUE)) {
      kn1 <- FNN::get.knn(xy, k = 1)
      spatial_scale <- stats::median(kn1$nn.dist, na.rm = TRUE)
      cleanup_radius_suggest <- max(5, min(10, round(spatial_scale * 0.6)))
    } else {
      cleanup_radius_suggest <- 6
    }
  }
  
  # regressors suggestion (now area-aware)
  regress_vars_suggest <- "nCount_Spatial"  # will be mapped to log10... inside V3
  mt_flag   <- FALSE
  area_flag <- FALSE
  rho_area_umi <- NA_real_
  
  if (!is.null(pc_cov_corr) && "percent.mt" %in% colnames(pc_cov_corr)) {
    if (any(abs(pc_cov_corr[, "percent.mt"]) > (corr_rho_thresh + 0.05), na.rm = TRUE)) {
      mt_flag <- TRUE
    }
  }
  
  # (A) global Spearman between area and log10_nCount
  if ("area_um2" %in% colnames(obj_tmp@meta.data) && "log10_nCount_Spatial" %in% colnames(obj_tmp@meta.data)) {
    rho_area_umi <- suppressWarnings(
      cor(obj_tmp$log10_nCount_Spatial, obj_tmp$area_um2, method = "spearman")
    )
    if (!is.na(rho_area_umi) && abs(rho_area_umi) > corr_rho_thresh) area_flag <- TRUE
  }
  # (B) PC-wise area correlations
  if (!is.null(pc_cov_corr) && "area_um2" %in% colnames(pc_cov_corr)) {
    if (any(abs(pc_cov_corr[, "area_um2"]) > corr_rho_thresh, na.rm = TRUE)) area_flag <- TRUE
  }
  
  if (mt_flag)   regress_vars_suggest <- c(regress_vars_suggest, "percent.mt")
  if (area_flag) regress_vars_suggest <- c(regress_vars_suggest, "area_um2")
  
  # harmony: suggest TRUE only if multiple batches (orig.ident/slide_id)
  batch_var <- if ("slide_id" %in% colnames(obj@meta.data)) "slide_id" else "orig.ident"
  use_harmony_suggest <- FALSE
  if (batch_var %in% colnames(obj@meta.data)) {
    if (length(unique(obj@meta.data[[batch_var]])) >= 2) use_harmony_suggest <- TRUE
  }
  
  # ============================
  # 6) Embedded analytics (optional exports)
  # ============================
  pc_loadings_df <- NULL
  sep_df <- NULL
  centroids_df <- NULL
  
  # PC loadings (top |loading| per direction)
  if (isTRUE(do_pc_loadings)) {
    pc_loadings_df <- .pc_loading_table(obj_tmp, reduction = "pca",
                                        pcs = pc_loadings_pcs, top_n = pc_loadings_top_n)
    .write_tsv(pc_loadings_df, "PC_LoadingsTop")
  }
  
  # PC~cluster metrics (if labels are available)
  if (isTRUE(do_cluster_metrics)) {
    labels <- NULL
    if (is.null(cluster_id_col)) {
      ids <- tryCatch(Idents(obj), error = function(e) NULL)
      if (!is.null(ids) && length(ids) == ncol(obj)) labels <- ids
    } else {
      if (cluster_id_col %in% colnames(obj@meta.data)) {
        labels <- obj[[cluster_id_col]][,1]
      } else {
        .msg(sprintf("Cluster column '%s' not found; skipping PC~cluster metrics.", cluster_id_col))
      }
    }
    
    if (!is.null(labels)) {
      pcs_for_cluster <- pc_cluster_pcs[pc_cluster_pcs <= ncol(pc_mat)]
      sep_df <- .pc_cluster_separation(obj_tmp, labels, reduction = "pca", pcs = pcs_for_cluster)
      centroids_df <- .pc_cluster_centroids(obj_tmp, labels,
                                            reduction = "pca",
                                            pcs = pcs_for_cluster[1:min(10, length(pcs_for_cluster))])
      .write_tsv(sep_df,       "PC_ClusterSeparation")
      .write_tsv(centroids_df, "PC_ClusterMeans")
    } else {
      .msg("No cluster labels provided/found; skip PC~cluster metrics.")
    }
  }
  
  # ============================
  # 7) Console report
  # ============================
  if (isTRUE(verbose)) {
    cat("\n========== Cellbin Preflight ==========\n")
    cat("* Assay in use       :", assay_in, "\n")
    cat("* Cells x Features   :", n_cells, "x", n_feats, "\n")
    cat("* Images             :", if (length(imgs)) paste(imgs, collapse=", ") else "none", "\n")
    cat("* Has coordinates    :", has_coords, "\n")
    if (has_coords) {
      rng_row <- range(obj$imagerow, na.rm = TRUE); rng_col <- range(obj$imagecol, na.rm = TRUE)
      cat(sprintf("  - imagerow range   : [%.1f, %.1f]\n", rng_row[1], rng_row[2]))
      cat(sprintf("  - imagecol range   : [%.1f, %.1f]\n", rng_col[1], rng_col[2]))
      if (!is.na(spatial_scale)) cat(sprintf("  - spatial NN scale : %.2f (median 1-NN dist)\n", spatial_scale))
    }
    # QC summaries (quantiles)
    qc_cols <- intersect(c("nCount_Spatial","nFeature_Spatial","percent.mt","log10_nCount_Spatial"),
                         colnames(obj@meta.data))
    for (v in qc_cols) {
      qv <- stats::quantile(obj@meta.data[[v]], probs = c(.05,.25,.5,.75,.95), na.rm = TRUE)
      cat(sprintf("  - %s quantiles     : %s\n", v, paste(signif(qv, 5), collapse=", ")))
    }
    # PC-tech correlations (top 5 PCs)
    if (!is.null(pc_cov_corr)) {
      cat("\nPC ~ technical Spearman rho (top 5 PCs):\n")
      top5 <- 1:min(5, ncol(pc_mat))
      print(round(pc_cov_corr[top5, , drop = FALSE], 3))
    }
    # area–UMI global correlation
    if ("area_um2" %in% colnames(obj_tmp@meta.data)) {
      cat(sprintf("\nSpearman(area_um2, log10_nCount_Spatial) = %.3f\n", rho_area_umi))
    }
    
    # Compute HVG suggestion to print (median feature count aware)
    nf_cols <- intersect(c("nFeature_Spatial","nFeature_RNA"), colnames(obj@meta.data))
    hv_suggest_print <- 1500
    if (length(nf_cols) > 0) {
      med_feat <- median(obj@meta.data[[nf_cols[1]]], na.rm = TRUE)
      hv_suggest_print <- if (is.finite(med_feat) && med_feat < 1000) 1200 else 1500
    }
    
    cat("\n-- Suggested parameters for V3 (coarse) --\n")
    cat("* use_harmony        :", use_harmony_suggest, "\n")
    cat("* regress_vars       :", paste(regress_vars_suggest, collapse = ", "), "\n")
    cat("* lognorm_features_n :", hv_suggest_print, "\n")
    cat("* dims_use           :", paste0("1:", dims_suggest), "\n")
    cat("* k_param            :", k_param_suggest, "\n")
    cat("* prune_snn          :", prune_snn_suggest, "\n")
    if (has_coords) cat("* cleanup_radius     :", cleanup_radius_suggest, "\n")
    cat("========================================\n\n")
  }
  
  # ============================
  # 8) Return bundle (added rho_area_umi)
  # ============================
  hv_suggest <- 1500
  nf_cols <- intersect(c("nFeature_Spatial","nFeature_RNA"), colnames(obj@meta.data))
  if (length(nf_cols) > 0) {
    med_feat <- median(obj@meta.data[[nf_cols[1]]], na.rm = TRUE)
    hv_suggest <- if (is.finite(med_feat) && med_feat < 1000) 1200 else 1500
  }
  
  invisible(list(
    assay_in              = assay_in,
    has_coords            = has_coords,
    pc_cov_corr           = pc_cov_corr,
    dims_use              = 1:dims_suggest,
    k_param               = k_param_suggest,
    prune_snn             = prune_snn_suggest,
    cleanup_radius        = cleanup_radius_suggest,
    lognorm_features_n    = hv_suggest,
    regress_vars          = regress_vars_suggest,
    use_harmony           = use_harmony_suggest,
    rho_area_umi          = rho_area_umi,         # ← NEW: global area–UMI correlation
    pc_loadings           = pc_loadings_df,
    pc_cluster_separation = sep_df,
    pc_cluster_centroids  = centroids_df
  ))
}



# ---- 1) 工具：安全取列为向量（避免 list-col 引发的 dplyr::count 错误）----
.get_vec <- function(df, col){
  x <- df[[col]]
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  as.vector(x)
}

# ---- 2) 保存原始标签，避免覆盖 ----
SaveBanksyRaw <- function(object, cols = c("banksy_domain","banksy_celltype")){
  for (co in cols) {
    if (!is.null(object@meta.data[[co]]) && is.null(object@meta.data[[paste0(co,"_raw")]])){
      object@meta.data[[paste0(co,"_raw")]] <- object@meta.data[[co]]
    }
  }
  object
}

# ---- 3) 空间平滑（多数投票），按批次/切片独立进行 ----
# iter: 平滑轮数；k: 每个点参照的邻居数；by: 批次列（如 "orig.ident"）
SmoothSpatialLabels <- function(object, label_col="banksy_domain", k=6, iter=1, by="orig.ident"){
  stopifnot(all(c("imagecol","imagerow") %in% colnames(object@meta.data)))
  lab <- .get_vec(object@meta.data, label_col)
  if (by %in% colnames(object@meta.data)) {
    batches <- as.character(object@meta.data[[by]])
  } else {
    batches <- rep("all", nrow(object@meta.data))
  }
  coords_all <- as.matrix(object@meta.data[, c("imagecol","imagerow")])
  
  suppressPackageStartupMessages(library(FNN))
  for (t in seq_len(iter)){
    new_lab <- lab
    for (b in unique(batches)){
      idx_b <- which(batches == b)
      coords <- coords_all[idx_b, , drop=FALSE]
      labs_b <- lab[idx_b]
      # KNN 仅在本批次内部查找
      nn <- FNN::get.knn(coords, k = k)$nn.index
      # 多数投票（含自身，稳定边界）
      for (i in seq_along(idx_b)){
        neigh <- labs_b[nn[i,]]
        cand  <- c(labs_b[i], neigh)
        new_lab[idx_b[i]] <- names(sort(table(cand), decreasing = TRUE))[1]
      }
    }
    lab <- new_lab
  }
  object@meta.data[[paste0(label_col,"_smoothed")]] <- factor(lab)
  object
}

# ---- 4) 合并超小簇（噪声岛），同样按批次进行 ----
PruneTinyDomains <- function(object, label_col="banksy_domain_smoothed", min_size=200, by="orig.ident"){
  stopifnot(all(c("imagecol","imagerow") %in% colnames(object@meta.data)))
  lab <- .get_vec(object@meta.data, label_col)
  if (by %in% colnames(object@meta.data)) {
    batches <- as.character(object@meta.data[[by]])
  } else {
    batches <- rep("all", nrow(object@meta.data))
  }
  coords_all <- as.matrix(object@meta.data[, c("imagecol","imagerow")])
  suppressPackageStartupMessages(library(FNN))
  
  new_lab <- lab
  for (b in unique(batches)){
    idx_b <- which(batches == b)
    labs_b <- lab[idx_b]
    tab_b  <- table(labs_b)
    small  <- names(tab_b[tab_b < min_size])
    if (length(small) == 0) next
    
    coords <- coords_all[idx_b, , drop=FALSE]
    nn <- FNN::get.knn(coords, k = 10)$nn.index
    
    for (s in small){
      where <- which(labs_b == s)
      for (ii in where){
        neigh <- labs_b[nn[ii,]]
        neigh <- neigh[neigh != s]
        if (length(neigh)) {
          new_lab[idx_b[ii]] <- names(sort(table(neigh), decreasing=TRUE))[1]
        } else {
          # 极端情况：如果全是同类邻居，就保持不变
          new_lab[idx_b[ii]] <- labs_b[ii]
        }
      }
    }
  }
  object@meta.data[[paste0(gsub("_smoothed$","",label_col), "_refined")]] <- factor(new_lab)
  object
}

# ---- 5) 一键封装：保存原始 → 平滑 → 合并小簇 ----
RefineBanksyDomains <- function(object,
                                label_col="banksy_domain",
                                by="orig.ident",
                                smooth_k=6, smooth_iter=1,
                                min_size=200){
  object <- SaveBanksyRaw(object, cols = c(label_col))
  object <- SmoothSpatialLabels(object, label_col=label_col, k=smooth_k, iter=smooth_iter, by=by)
  object <- PruneTinyDomains(object, label_col=paste0(label_col,"_smoothed"),
                             min_size=min_size, by=by)
  object
}



SmallFragmentFraction <- function(object, label_col="banksy_domain_refined",
                                  k=6, cutoff=50,
                                  coord_cols=c("imagecol","imagerow")){
  meta <- object@meta.data
  xy <- as.matrix(meta[, coord_cols])
  nn <- get.knn(xy, k=k)$nn.index
  g  <- graph_from_edgelist(cbind(rep(1:nrow(nn), ncol(nn)), as.vector(nn)), directed=FALSE)
  
  labs <- meta[[label_col]]
  comp_size <- integer(nrow(meta)); comp_size[] <- NA_integer_
  
  for (lv in unique(labs)) {
    nodes <- which(labs == lv)
    if (!length(nodes)) next
    comps <- components(induced_subgraph(g, nodes))
    sizes <- comps$csize
    comp_size[nodes] <- sizes[comps$membership]
  }
  mean(comp_size <= cutoff, na.rm=TRUE)
}



DomainCoherence <- function(obj, label="banksy_domain_refined", k=6, coords=c("imagecol","imagerow")){
  lab <- obj[[label]][,1]; xy <- as.matrix(obj@meta.data[,coords])
  nn  <- get.knn(xy, k=k)$nn.index
  frac_same <- rowMeans(matrix(lab[nn], ncol=k) == lab)
  by(lab, lab, function(ix) mean(frac_same[which(lab==ix)], na.rm=TRUE))
}




suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(Matrix); library(FNN); library(ggplot2)
})

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(Matrix); library(FNN)
})

BuildNicheAssay_HD <- function(
    object,
    group.by            = "banksy_celltype_detailed",
    neighbors.k         = 30,         # 基本K数
    radius              = NULL,       # 若给半径就用半径邻域
    include.self        = FALSE,
    min_label_fraction  = 0.005,      # 低频标签合并为 "Other"
    restrict.within     = NULL,       # domain 列名（用于判定“跨域”）
    alpha_cross         = 0.3,        # ★ 跨域邻居权重（0~1，0=硬限制，1=不降权）
    dist_kernel         = c("none","gaussian","epanechnikov"), # 距离核
    sigma               = NULL,       # 高斯核宽度（缺省自动取邻距中位数）
    assay.name          = "niche_soft",
    dims.niche          = 1:20,
    niche.resolution    = 0.4,
    seed                = 123
){
  stopifnot(group.by %in% colnames(object@meta.data))
  set.seed(seed)
  dist_kernel <- match.arg(dist_kernel)
  
  ## 1) 坐标
  coords <- tryCatch({
    GetTissueCoordinates(object)[, c("imagecol", "imagerow"), drop = FALSE]
  }, error = function(e){
    md <- object@meta.data
    xs <- grep("imagecol|x|col", colnames(md), ignore.case=TRUE, value=TRUE)[1]
    ys <- grep("imagerow|y|row", colnames(md), ignore.case=TRUE, value=TRUE)[1]
    if (is.na(xs) || is.na(ys)) stop("找不到空间坐标列。")
    md[, c(xs, ys), drop=FALSE]
  })
  colnames(coords) <- c("x","y")
  coords <- as.matrix(coords)
  n <- nrow(coords)
  
  dom <- if (!is.null(restrict.within)) as.character(object[[restrict.within]][,1]) else rep(NA, n)
  labels <- as.character(object[[group.by]][,1]); names(labels) <- colnames(object)
  
  ## 2) 合并低频标签为 "Other"
  lab_tab <- sort(table(labels), decreasing = TRUE)
  keep_labs <- names(lab_tab[ lab_tab/sum(lab_tab) >= min_label_fraction ])
  labels2 <- ifelse(labels %in% keep_labs, labels, "Other")
  labs <- sort(unique(labels2))
  
  ## 3) 建邻域（KNN 或半径）
  # 为了拿到每个邻居的距离，用 FNN 的距离输出
  if (is.null(radius)) {
    k_eff <- min(neighbors.k + as.integer(!include.self), n-1)
    knn <- get.knn(coords, k = k_eff)
    nn_index <- knn$nn.index
    nn_dist  <- knn$nn.dist
  } else {
    # 先取较大的 K，再按半径筛
    k_cap <- min(max(50, neighbors.k*3), n-1)
    knn <- get.knn(coords, k = k_cap)
    nn_index <- knn$nn.index
    nn_dist  <- knn$nn.dist
  }
  
  # 自动 sigma（用于高斯核）
  if (dist_kernel == "gaussian" && (is.null(sigma) || !is.finite(sigma))) {
    sigma <- stats::median(nn_dist, na.rm = TRUE)
    if (!is.finite(sigma) || sigma <= 0) sigma <- 1
  }
  
  ## 4) 计算“加权计数”
  mat_counts <- matrix(0, nrow = length(labs), ncol = n,
                       dimnames = list(labs, colnames(object)))
  
  for (i in seq_len(n)) {
    ids <- nn_index[i, ]; d  <- nn_dist[i, ]
    ids <- ids[!is.na(ids) & ids > 0]
    d   <- d[seq_along(ids)]
    
    # 半径模式：筛选
    if (!is.null(radius)) {
      hit <- which(d <= radius)
      ids <- ids[hit]; d <- d[hit]
    }
    
    # 去掉自身
    if (!include.self) {
      keep <- ids != i
      ids <- ids[keep]; d <- d[keep]
    }
    
    # 若不足 neighbors.k，则截断/补齐
    if (length(ids) > neighbors.k) {
      ids <- ids[1:neighbors.k]; d <- d[1:neighbors.k]
    } else if (length(ids) == 0) {
      next
    }
    
    # 领域内/跨域权重
    w_dom <- rep(1, length(ids))
    if (!all(is.na(dom))) {
      w_dom[ dom[ids] != dom[i] ] <- alpha_cross
    }
    
    # 距离权重
    if (dist_kernel == "gaussian") {
      w_dist <- exp( - (d^2) / (2 * sigma^2) )
    } else if (dist_kernel == "epanechnikov") {
      # h 用 2*sigma 作为带宽
      h <- if (is.null(sigma)) stats::median(d, na.rm = TRUE) * 2 else 2*sigma
      w_dist <- pmax(0, 1 - (d/h)^2)
    } else {
      w_dist <- rep(1, length(ids))
    }
    
    w <- w_dom * w_dist
    if (!any(is.finite(w)) || sum(w) == 0) w <- rep(1, length(ids))  # 回退
    
    # 按标签累加“加权计数”
    labs_nb <- labels2[ids]
    # 用 tapply 汇总到标签
    w_sum <- tapply(w, labs_nb, sum)
    mat_counts[names(w_sum), i] <- mat_counts[names(w_sum), i] + as.numeric(w_sum)
  }
  
  # 归一化为比例（合计=1），避免极少邻居导致的尺度差异
  col_sums <- colSums(mat_counts)
  col_sums[col_sums == 0] <- 1
  mat_prop <- sweep(mat_counts, 2, col_sums, "/")
  
  ## 5) 写入新 assay（counts=加权计数，data=比例）
  object[[assay.name]] <- CreateAssayObject(counts = mat_counts)
  object <- SetAssayData(object, assay = assay.name, slot = "data", new.data = mat_prop)
  
  # 映射保存（下划线→连字符）
  object@misc[[paste0(assay.name, "_feature_map")]] <- data.frame(
    original = sort(unique(as.character(object[[group.by]][,1]))),
    feature  = rownames(object[[assay.name]])
  )
  
  ## 6) niche 空间降维聚类（自适应PC）
  DefaultAssay(object) <- assay.name
  VariableFeatures(object) <- rownames(object[[assay.name]])
  
  object <- ScaleData(object, assay = assay.name, verbose = FALSE)
  num_features <- nrow(object[[assay.name]])
  npcs_to_calc <- max(2, min(max(dims.niche), num_features - 1))
  red_name <- paste0("pca.", assay.name)
  
  object <- RunPCA(object, assay = assay.name, npcs = npcs_to_calc,
                   verbose = FALSE, reduction.name = red_name)
  
  actual_pcs <- ncol(Embeddings(object, reduction = red_name))
  dims_to_use <- 1:min(max(dims.niche), actual_pcs)
  
  object <- FindNeighbors(object, reduction = red_name, dims = dims_to_use)
  object <- FindClusters(object, resolution = niche.resolution,
                         cluster.name = paste0(assay.name, "_cluster"))
  object <- RunUMAP(object, reduction = red_name, dims = dims_to_use,
                    reduction.name = paste0("umap.", assay.name))
  
  # 参数记录，便于复现
  object@misc[[paste0(assay.name, "_params")]] <- list(
    group.by = group.by, neighbors.k = neighbors.k, radius = radius,
    include.self = include.self, min_label_fraction = min_label_fraction,
    restrict.within = restrict.within, alpha_cross = alpha_cross,
    dist_kernel = dist_kernel, sigma = sigma,
    dims.niche = dims.niche, niche.resolution = niche.resolution
  )
  object
}


BuildNicheAssay_Advanced <- function(
    object,
    group.by,
    k_values            = c(15, 50),
    dist_kernels        = c("gaussian", "none"),
    alpha_cross         = 0.4,            # ← 默认开启软限制
    restrict.within     = NULL,
    min_label_fraction  = 0.005,
    assay.name          = "niche_adv",
    normalization       = "CLR"
){
  stopifnot(group.by %in% colnames(object@meta.data))
  stopifnot(length(k_values) == length(dist_kernels))
  
  coords <- tryCatch({
    GetTissueCoordinates(object)[, c("imagecol","imagerow"), drop=FALSE]
  }, error=function(e){
    md <- object@meta.data
    xs <- grep("imagecol|x|col", colnames(md), ignore.case=TRUE, value=TRUE)[1]
    ys <- grep("imagerow|y|row", colnames(md), ignore.case=TRUE, value=TRUE)[1]
    if (is.na(xs) || is.na(ys)) stop("Spatial coordinate columns not found.")
    md[, c(xs, ys), drop=FALSE]
  })
  colnames(coords) <- c("x","y"); coords <- as.matrix(coords); n <- nrow(coords)
  
  dom    <- if (!is.null(restrict.within)) setNames(object[[restrict.within, drop=TRUE]], colnames(object)) else setNames(rep(NA, n), colnames(object))
  labels <- setNames(object[[group.by, drop=TRUE]], colnames(object))
  
  tab <- sort(table(labels), decreasing=TRUE)
  keep <- names(tab[ tab/sum(tab) >= min_label_fraction ])
  labels2 <- setNames(ifelse(labels %in% keep, labels, "Other"), names(labels))
  labs <- sort(unique(labels2))
  
  # 预先一次性取最大K的KNN（自包含对角）
  max_k <- max(k_values)
  knn_max <- FNN::get.knnx(coords, coords, k = max_k + 1)  # +1 for self
  all_counts_list <- list()
  
  # —— σ 的稳健估计：仅用非自邻、正距离 ——
  pos_d <- knn_max$nn.dist[, -1, drop=FALSE]
  pos_d <- pos_d[pos_d > 0]
  sigma_global <- ifelse(length(pos_d) > 0, stats::median(pos_d, na.rm=TRUE), 1)
  
  for (j in seq_along(k_values)) {
    k <- k_values[j]; kernel <- dist_kernels[j]
    mat_counts_k <- matrix(0, nrow=length(labs), ncol=n, dimnames=list(labs, colnames(object)))
    sigma <- if (kernel == "gaussian") sigma_global else 1
    
    for (i in seq_len(n)) {
      nei_idx <- knn_max$nn.index[i, -1][1:k]
      nei_dis <- knn_max$nn.dist[i,  -1][1:k]
      nei_nm  <- colnames(object)[nei_idx]
      
      # 距离权重
      w_dist <- if (kernel == "gaussian") exp(-(nei_dis^2)/(2*sigma^2)) else rep(1, k)
      # 软限制（若提供域信息）
      if (!is.null(restrict.within)) {
        w_dom <- ifelse(dom[nei_nm] == dom[i], 1, alpha_cross)
      } else w_dom <- rep(1, k)
      
      w <- w_dist * w_dom
      w_sum <- tapply(w, labels2[nei_nm], sum, na.rm = TRUE)
      if (length(w_sum) > 0) mat_counts_k[names(w_sum), i] <- as.numeric(w_sum)
    }
    rownames(mat_counts_k) <- paste0("k", k, "_", rownames(mat_counts_k))
    all_counts_list[[j]] <- mat_counts_k
  }
  
  final_counts <- do.call(rbind, all_counts_list)
  
  # 写入 assay：counts=加权计数（稀疏），data=CLR（先比例）
  library(Matrix)
  object[[assay.name]] <- CreateAssayObject(counts = Matrix(final_counts, sparse=TRUE))
  props <- sweep(final_counts, 2, pmax(1, colSums(final_counts)), "/")
  
  if (normalization == "CLR") {
    eps <- 1e-4
    logP <- log(pmax(props, eps))
    clr  <- sweep(logP, 2, colMeans(logP), "-")
    object <- SetAssayData(object, assay=assay.name, slot="data", new.data = clr)
  } else {
    object <- SetAssayData(object, assay=assay.name, slot="data", new.data = log1p(props))
  }
  
  # 轻度降维（只准备好 PCA；聚类放到下一步的一键函数里做）
  DefaultAssay(object) <- assay.name
  VariableFeatures(object) <- rownames(final_counts)
  object <- ScaleData(object, assay=assay.name, verbose=FALSE)
  npcs_to_calc <- min(30, nrow(object[[assay.name]]) - 1)
  object <- RunPCA(object, assay=assay.name, npcs=npcs_to_calc, verbose=FALSE,
                   reduction.name = paste0("pca.", assay.name))
  object
}


# 如果你还没加载我给的后端函数，复制这三行：
# library(Seurat)
ComputeNicheMap_Final <- function(
    object, assay.name="niche_adv",
    vf_top=60, pca_npcs=12,
    graph_k_param=15, prune_snn=0.1, resolution=0.07,
    min_cluster_size=300, rare_threshold=150, target_k=NULL
){
  DefaultAssay(object) <- assay.name
  
  # 选前 vf_top 特征（按 counts 总和排序）
  sums <- rowSums(GetAssayData(object, assay=assay.name, slot="counts"))
  vf <- names(sort(sums, decreasing=TRUE))[1:min(vf_top, length(sums))]
  VariableFeatures(object) <- vf
  
  # 用现成的 CLR（data 槽）做 PCA（再来一次更稳）
  object <- ScaleData(object, assay=assay.name, do.center=FALSE, do.scale=FALSE, verbose=FALSE)
  object <- RunPCA(object, assay=assay.name, npcs=pca_npcs, verbose=FALSE,
                   reduction.name=paste0("pca.",assay.name), reduction.key="NICHEPC_")
  
  if (is.numeric(target_k)) {
    emb <- Embeddings(object, paste0("pca.",assay.name))[, 1:pca_npcs]
    set.seed(1); km <- kmeans(emb, centers=target_k, nstart=25)
    object[[paste0(assay.name,"_final")]] <- factor(km$cluster)
  } else {
    object <- FindNeighbors(object, reduction=paste0("pca.",assay.name), dims=1:pca_npcs,
                            k.param=graph_k_param, nn.method="annoy",
                            prune.SNN=prune_snn, cache.index=TRUE,
                            graph.name=paste0(assay.name,"_snn"))
    object <- FindClusters(object, graph.name=paste0(assay.name,"_snn"),
                           algorithm=4, resolution=resolution,
                           cluster.name=paste0(assay.name,"_final"))
  }
  
  # 小簇并入
  lab <- object[[paste0(assay.name,"_final")]][,1]
  tab <- sort(table(lab), decreasing=TRUE)
  small <- names(tab[tab < min_cluster_size])
  if (length(small)>0) {
    emb  <- Embeddings(object, paste0("pca.",assay.name))[, 1:pca_npcs]
    cent <- rowsum(emb, lab) / as.numeric(tab)
    big  <- names(tab[tab >= min_cluster_size])
    for (s in small) {
      d <- colSums((t(cent[big,,drop=FALSE]) - cent[s,])^2)
      lab[ lab == s ] <- big[ which.min(d) ]
    }
    object[[paste0(assay.name,"_final")]] <- factor(lab)
  }
  
  # 稀有类汇总
  lab2 <- as.character(object[[paste0(assay.name,"_final")]][,1])
  tab2 <- table(lab2); rare <- names(tab2[tab2 < rare_threshold])
  if (length(rare)>0) {
    lab2[ lab2 %in% rare ] <- "niche_rare"
    object[[paste0(assay.name,"_final_coarse")]] <- factor(lab2)
  }
  object
}




RunInferCNV_HD <- function(
    object,
    celltype_col,
    observation_pattern = "^Epithelial", # 用于识别观测组(通常是上皮/肿瘤)的正则匹配模式
    ref_subsample_cap = 30000,           # 参考组下采样上限，防止内存问题
    use_spatial_smoothing = TRUE,        # 是否进行上游 KNN 空间平滑
    knn_k = 10,                          # KNN 平滑的邻居数
    gene_order_file,                     # inferCNV必需的基因顺序文件路径
    result_dir,                          # 所有inferCNV结果的输出目录
    num_threads = 8,                     # inferCNV运行的线程数
    HMM = TRUE,                          # 是否运行HMM来预测离散的CNV状态
    denoise = TRUE                       # 是否进行去噪
){
  
  # --- 0. 载入依赖并进行参数检查 ---
  suppressPackageStartupMessages({
    library(Seurat); library(dplyr); library(Matrix); library(FNN); library(infercnv)
  })
  
  stopifnot(inherits(object, "Seurat"))
  stopifnot(celltype_col %in% colnames(object@meta.data))
  stopifnot(file.exists(gene_order_file))
  
  message("--- 开始为高密度空间数据运行 inferCNV 流程 (v2.0) ---")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- 1. 准备输入数据：细胞选择与Counts矩阵 ---
  message("\nStep 1: 准备输入数据...")
  
  is_obs <- grepl(observation_pattern, object[[celltype_col, drop=TRUE]], ignore.case = TRUE)
  obs_ids <- colnames(object)[is_obs]
  ref_ids <- colnames(object)[!is_obs]
  
  message(sprintf("  - 观测组 (%s): %d 个细胞", observation_pattern, length(obs_ids)))
  message(sprintf("  - 参考组: %d 个细胞", length(ref_ids)))
  
  if (length(ref_ids) > ref_subsample_cap) {
    set.seed(1)
    ref_ids <- sample(ref_ids, ref_subsample_cap)
    message(sprintf("  - 参考组已下采样至 %d 个细胞", ref_subsample_cap))
  }
  
  sel_cells <- c(obs_ids, ref_ids)
  object_sub <- subset(object, cells = sel_cells)
  
  # 获取原始counts
  DefaultAssay(object_sub) <- "Spatial"
  counts <- GetAssayData(object_sub, slot = "counts")
  
  if (isTRUE(use_spatial_smoothing)) {
    message("  - 正在对 counts 矩阵进行 KNN 空间平滑 (k=", knn_k, ")...")
    coords <- GetTissueCoordinates(object_sub, cols = c("imagecol", "imagerow"))
    
    knn <- FNN::get.knn(as.matrix(coords), k = knn_k)
    i <- rep(1:nrow(knn$nn.index), each = knn_k)
    j <- as.vector(t(knn$nn.index)) # 修正：FNN返回的index需要转置再拉平
    
    # 使用 ijx 格式创建权重矩阵，更高效
    W <- sparseMatrix(i = i, j = j, x = 1, dims = c(ncol(counts), ncol(counts)))
    W <- Matrix::rowSums(W)
    W@x <- 1 / W@x[W@i]
    
    counts_smooth <- counts %*% W
    colnames(counts_smooth) <- colnames(counts)
    counts <- counts_smooth
    message("  - 空间平滑完成。")
    rm(counts_smooth, W, knn, coords); gc() # 清理内存
  }
  
  # --- 2. 准备 inferCNV 的输入文件 ---
  message("\nStep 2: 准备 inferCNV 注释和基因顺序文件...")
  
  # 与输入 counts 矩阵对齐基因顺序
  gene_order_df <- readr::read_tsv(gene_order_file, col_names = c("gene", "chr", "start", "end"), col_types = "ccii")
  keep_genes <- intersect(rownames(counts), gene_order_df$gene)
  
  counts_ord <- counts[keep_genes, ]
  gene_order_df_ord <- gene_order_df[match(keep_genes, gene_order_df$gene), ]
  
  gene_order_final_path <- file.path(result_dir, "infercnv.gene_order.txt")
  readr::write_tsv(gene_order_df_ord, gene_order_final_path, col_names = FALSE)
  
  # 创建注释文件
  ann_df <- data.frame(
    cell = colnames(counts_ord),
    group = ifelse(colnames(counts_ord) %in% obs_ids, "Observation", "Reference")
  )
  ann_final_path <- file.path(result_dir, "infercnv.annotations.txt")
  readr::write_tsv(ann_df, ann_final_path, col_names = FALSE)
  
  message("  - 文件已生成并保存至: ", normalizePath(result_dir, winslash="/"))
  
  # --- 3. 创建并运行 inferCNV 对象 ---
  message("\nStep 3: 创建并运行 inferCNV 对象 (这会非常耗时)...")
  
  infer_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts_ord,
    annotations_file  = ann_final_path,
    delim             = "\t",
    gene_order_file   = gene_order_final_path,
    ref_group_names   = c("Reference")
  )
  
  # 运行 inferCNV
  infer_obj_run <- infercnv::run(
    infercnv_obj        = infer_obj,
    cutoff              = 0.1,
    out_dir             = result_dir,
    cluster_by_groups   = TRUE,
    denoise             = denoise,
    HMM                 = HMM,
    HMM_type            = "i6",
    num_threads         = num_threads,
    write_expr_matrix   = TRUE # 确保输出结果矩阵
  )
  
  message("  - inferCNV 运行完成。")
  
  # --- 4. 后处理：计算CNV分数并添加回Seurat对象 ---
  message("\nStep 4: 后处理，计算 CNV 分数并更新 Seurat 对象...")
  
  # 读取去噪后的表达矩阵
  expr_path <- file.path(result_dir, "infercnv.denoised.expr.matrix")
  if(!file.exists(expr_path)) { # HMM 模式下的备用路径
    expr_path <- file.path(result_dir, "infercnv.17_HMM_predHMMi6.hmm_mode-samples.Pnorm_0.5.repr_intensities.observations.txt")
  }
  if(!file.exists(expr_path)) stop("找不到 inferCNV 的输出表达矩阵。请检查运行是否成功。")
  
  E <- as.matrix(read.table(expr_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE))
  
  # 计算 CNV 负荷分数（中位绝对偏移）
  cnv_score <- apply(abs(E - 1), 2, median) # inferCNV去噪后，正常表达中心为1
  
  # 定义阈值并进行肿瘤状态调用
  ref_scores <- cnv_score[intersect(names(cnv_score), ref_ids)]
  thr <- quantile(ref_scores, 0.95, na.rm = TRUE)
  message("  - CNV 分数阈值 (参考组95分位): ", signif(thr, 3))
  
  pred_label <- ifelse(cnv_score > thr, "aneuploid_like", "diploid_like")
  
  # 创建结果数据框
  result_df <- data.frame(
    infercnv_cnv_score = cnv_score[colnames(object)],
    infercnv_call = pred_label[colnames(object)]
  )
  rownames(result_df) <- colnames(object)
  
  # 添加回主对象
  object <- AddMetaData(object, result_df)
  
  # 最终的、结合细胞类型的肿瘤上皮标签
  is_obs_full <- grepl(observation_pattern, object[[celltype_col, drop=TRUE]], ignore.case = TRUE)
  
  final_call <- case_when(
    is_obs_full & object$infercnv_call == "aneuploid_like" ~ "tumor_epithelial",
    is_obs_full & object$infercnv_call == "diploid_like" ~ "normal_epithelial",
    TRUE ~ "non_epithelial"
  )
  
  object$epi_tumor_infercnv <- factor(final_call,
                                      levels=c("tumor_epithelial","normal_epithelial","non_epithelial"))
  
  message("  - CNV 结果已成功添加回 Seurat 对象。")
  
  saveRDS(object, file.path(result_dir, "object_with_infercnv_calls.rds"))
  
  message("--- inferCNV 流程全部完成！ ---")
  return(object)
}




# 创建一个列表来存储两个基因集
pdac_subtypes_gene_list <- list(
  
  # ===================================================================
  # Classical Subtype Gene Set
  # 特征：上皮/腺体分化、细胞黏附、黏蛋白分泌。
  # 核心转录因子：GATA6, HNF4A, HNF1A, KLF4
  # ===================================================================
  "Classical" = c(
    # --- 核心转录调控 ---
    "GATA6",     # <<< Classical 亚型的“主开关” (Master Regulator)
    "HNF4A",     # 肝细胞核因子4α, GATA6 的上游
    "HNF1A",     # 肝细胞核因子1α, GATA6 的下游
    "PDX1",      # 胰腺和十二指肠同源异形盒基因1, 在分化良好的肿瘤中表达
    "KLF4",      # Krüppel样因子4
    
    # --- 角蛋白 (Keratins) & 黏附 ---
    "KRT8",      # 经典的腺上皮角蛋白
    "KRT18",     # 经典的腺上皮角蛋白, 与KRT8配对
    "KRT19",     # 经典的腺上皮角蛋白
    "KRT7",      # 常在胰腺导管和肿瘤中表达
    "KRT20",     # 经典的肠道分化标记
    "EPCAM",     # <<< 上皮细胞黏附分子, 最经典的泛上皮/癌细胞标志物
    "CLDN4",     # Claudin 4, 紧密连接蛋白
    "CLDN18",    # Claudin 18
    "TJP1",      # ZO-1, 紧密连接蛋白
    
    # --- 分泌蛋白 (黏蛋白 & 三叶因子) ---
    "MUC1",      # 跨膜黏蛋白, 经典PDAC标记
    "MUC5AC",    # 分泌型黏蛋白, 胃源性化生 (gastric metaplasia) 的标志
    "MUC13",     # 黏蛋白
    "MUC20",     # 黏蛋白
    "TFF1",      # 三叶因子1, 常与MUC5AC共表达, 胃源性
    "TFF3",      # 三叶因子3, 肠道/杯状细胞标记
    "AGR2",      # 前生长因子2, 参与黏液生成和ER应激
    "AGR3",      # AGR2的同源物
    "REG4",      # 再生岛衍生蛋白4, 肠道干/祖细胞标记
    
    # --- 其他经典型标记 ---
    "MSLN",      # Mesothelin, 经典的PDAC表面抗原
    "CEACAM6",   # 癌胚抗原相关细胞黏附分子6
    "SLC4A4",    # 碳酸氢根转运蛋白
    "ERBB2",     # HER2, 在部分Classical亚型中扩增/高表达
    "KIAA1324"   # (也称为 EVISECT)
  ),
  
  # ===================================================================
  # Basal-like / Squamous Subtype Gene Set
  # 特征：基底细胞/鳞状分化、侵袭、迁移、炎症、糖酵解。
  # 核心转录因子：TP63 (ΔNp63)
  # ===================================================================
  "Basal_like" = c(
    # --- 核心转录调控 ---
    "TP63",      # <<< Basal-like 亚型的“主开关” (Master Regulator), 通常是ΔNp63亚型
    
    # --- 基底/鳞状角蛋白 (Basal/Squamous Keratins) ---
    "KRT17",     # <<< 最经典的 Basal-like 标记之一
    "KRT5",      # 基底细胞角蛋白
    "KRT6A",     # 基底/增殖性角蛋白
    "KRT6B",     # 基底/增殖性角蛋白
    "KRT14",     # 基底细胞角蛋白
    "KRT15",     # 基底细胞角蛋白
    "KRT16",     # 增殖/应激相关角蛋白
    
    # --- 细胞黏附 & 迁移 (EMT相关) ---
    "LAMC2",     # Laminin C2, 基底膜组分, 与侵袭相关
    "LAMB3",     # Laminin B3
    "ITGA6",     # Integrin α6, Laminin受体
    "ITGB4",     # Integrin β4, 与ITGA6配对
    "ITGA2",     # Integrin α2
    "ITGA3",     # Integrin α3
    "COL17A1",   # 胶原蛋白XVII, 半桥粒组分
    "DSG3",      # Desmoglein-3, 桥粒组分, 鳞状上皮标记
    "VIM",       # Vimentin, 经典的间充质标记 (注意可能与CAF信号重叠)
    "FN1",       # Fibronectin 1, 经典的间充质标记 (注意可能与CAF信号重叠)
    "PDPN",      # Podoplanin, 淋巴管和鳞状癌标记
    
    # --- 炎症信号 ---
    "CXCL8",     # IL-8, 强力的中性粒细胞趋化因子
    "IL1A",      # 白细胞介素1α
    "IL1B",      # 白细胞介素1β
    
    # --- 蛋白酶 & 抑制剂 ---
    "SERPINB5",  # Maspin, 丝氨酸蛋白酶抑制剂
    "MMP1",      # 基质金属蛋白酶1
    "MMP10",     # 基质金属蛋白酶10
    "MMP12",     # 基质金属蛋白酶12
    
    # --- 其他 Basal-like 标记 ---
    "S100A2",    # S100 蛋白家族, 经典的 Basal-like 标记
    "S100A11",   # S100 蛋白家族
    "ANXA1",     # Annexin A1
    "SLC2A1",    # GLUT1, 葡萄糖转运蛋白, 与糖酵解/缺氧相关
    "PLEK2"      # Pleckstrin 2
  )
)
emt_gene_showcase_list <- list(
  
  # ==========================================================
  # Epithelial Program (EXPECTED to be HIGH in Classical, LOW in EMT clusters)
  # ==========================================================
  "Epithelial Identity" = c(
    "EPCAM",   # The most canonical epithelial cell adhesion molecule.
    "KRT19",   # A classic simple-epithelial keratin.
    "KRT8",    # Pairs with KRT18 in glandular epithelia.
    "MUC1",    # A key mucin indicating glandular function.
    "CDH1"     # E-Cadherin, the cornerstone of epithelial cell-cell adhesion.
  ),
  
  # ==========================================================
  # Mesenchymal Program (EXPECTED to be LOW in Classical, HIGH in EMT clusters)
  # ==========================================================
  "Mesenchymal Identity" = c(
    # --- Core EMT Transcription Factors (The "Master Switches") ---
    "ZEB1",    # Zinc Finger E-Box Binding Homeobox 1 (The most important one in PDAC)
    "SNAI1",   # Snail
    "TWIST1",  # Twist
    
    # --- Mesenchymal Structural & ECM Genes ---
    "VIM",     # Vimentin (Classic mesenchymal intermediate filament)
    "FN1",     # Fibronectin 1 (Key extracellular matrix protein)
    "COL1A1",  # Collagen Type I Alpha 1 (Hallmark of fibroblasts/mesenchymal cells)
    "SPARC",   # Secreted Protein Acidic and Cysteine Rich
    
    # --- Myofibroblastic / Contractile Genes ---
    "ACTA2",   # Alpha-Smooth Muscle Actin (a-SMA)
    "TAGLN"    # Transgelin
  )
)


## =========================== GMT 工具集 ===========================
`%||%` <- function(a,b) if (is.null(a)) b else a

# 1) 读取 GMT：返回 named list，每个元素是该基因集的基因向量
read_gmt <- function(file) {
  if (!file.exists(file)) stop("GMT not found: ", file)
  ln <- readLines(file, warn = FALSE)
  parsed <- lapply(ln, function(x){
    p <- strsplit(x, "\t")[[1]]
    name <- p[1]
    desc <- if (length(p) >= 2) p[2] else NA_character_
    genes <- unique(p[-c(1,2)])
    list(name=name, desc=desc, genes=genes)
  })
  sets <- setNames(lapply(parsed, `[[`, "genes"), vapply(parsed, `[[`, "", "name"))
  attr(sets, "desc") <- setNames(vapply(parsed, function(x) x$desc %||% "", ""), names(sets))
  sets
}

# 2) 列出/筛选通路名（名字或描述里匹配）
find_pathways <- function(sets, query, mode = c("auto","exact","contains","regex"), ignore_case = TRUE) {
  mode <- match.arg(mode)
  nm  <- names(sets)
  ds  <- attr(sets, "desc"); ds <- ds[nm]
  .match <- function(pat, x) grepl(pat, x, ignore.case = ignore_case, perl = TRUE)
  if (mode == "exact") {
    hits <- nm[ tolower(nm) == tolower(query) ]
  } else if (mode == "contains" || (mode=="auto" && nchar(query) <= 4)) {
    hits <- nm[ .match(query, nm) | .match(query, ds) ]
  } else if (mode == "regex" || mode=="auto") {
    hits <- nm[ .match(query, nm) | .match(query, ds) ]
  }
  hits
}

# 3) 常用别名（可扩展）
default_alias <- function() {
  list(
    EMT = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "EPITHELIAL_MESENCHYMAL_TRANSITION","EMT")
  )
}

# 4) 主函数：从目录/文件提取通路 → features_list（可选与对象交集）
extract_pathways_from_dir <- function(
    dir        = "E:/modifiedcodeR/",
    gmt_file   = "h.all.v2024.1.Hs.symbols.gmt",
    queries    = "EMT",
    mode       = c("auto","exact","contains","regex"),
    ignore_case= TRUE,
    alias_map  = default_alias(),
    filter_present = TRUE,           # 仅保留对象里出现的基因
    object     = NULL,               # 若 filter_present=TRUE，建议提供 Seurat 对象
    assay      = "Spatial",          # 提取对象里基因时用
    min_genes  = 5                   # 过滤后至少保留这么多基因
){
  mode <- match.arg(mode)
  file <- file.path(dir, gmt_file)
  sets <- read_gmt(file)
  
  # 构造查询名列表（考虑别名）
  expand_queries <- function(q){
    if (!is.null(alias_map[[q]])) alias_map[[q]] else q
  }
  q_all <- unique(unlist(lapply(queries, expand_queries), use.names = FALSE))
  
  # 找到匹配到的基因集名
  hit_names <- unique(unlist(lapply(q_all, function(q) find_pathways(sets, q, mode, ignore_case))))
  if (length(hit_names) == 0) stop("No pathway matched: ", paste(queries, collapse=", "))
  
  # 取出并可选与对象交集
  out <- lapply(hit_names, function(nm) sets[[nm]])
  names(out) <- hit_names
  
  if (isTRUE(filter_present) && !is.null(object)) {
    # 尽量从 data 层拿行名；退而求其次 counts
    present <- tryCatch(rownames(Seurat::GetAssayData(object, assay=assay, layer="data")),
                        error=function(e) tryCatch(rownames(Seurat::GetAssayData(object, assay=assay, layer="counts")),
                                                   error=function(e2) rownames(object)))
    out <- lapply(out, function(v) intersect(unique(v), present))
  }
  
  # 去掉过短的 gene set
  n_gene <- vapply(out, length, integer(1))
  keep   <- n_gene >= min_genes
  if (!all(keep)) {
    message("Dropped too-small sets: ", paste(names(out)[!keep], collapse = ", "))
    out <- out[keep]
  }
  if (length(out) == 0) stop("All matched pathways became too small after filtering.")
  
  out
}


#' @title Create a Robust Color Palette for Counted Labels
#' @description Takes a set of counted labels and a base color palette. It intelligently
#'              handles complex names in the base palette (e.g., names containing parentheses)
#'              to create a new palette keyed by the counted labels.
#'
#' @param levels_counted A character vector of the "counted" factor levels 
#'                       (e.g., "CAF (22452)", "Myeloid (Mono/Macro/DC) (12149)").
#' @param base_palette A named character vector where names are the original labels,
#'                     which may themselves contain special characters.
#'                     (e.g., c("CAF" = "red", "Myeloid (Mono/Macro/DC)" = "blue")).
#'
#' @return A named character vector representing the new palette, correctly matched.
#'
CreateMatchingPalette <- function(levels_counted, base_palette) {
  
  # --- 1. 创建用于匹配的“干净键” ---
  
  # a) 从带计数的标签中提取原始标签 (移除最后的括号)
  #    e.g., "Myeloid (Mono/Macro/DC) (12149)" -> "Myeloid (Mono/Macro/DC)"
  keys_from_counted <- sub(" \\([0-9,]+\\)$", "", levels_counted)
  
  # b) 从基础调色板的名字中也提取“干净键”
  #    这一步是新增的！它确保基础调色板的键也是干净的。
  #    e.g., "Myeloid (Mono/Macro/DC)" -> "Myeloid (Mono/Macro/DC)" (保持不变)
  #    e.g., "Perivascular/Smooth muscle" -> "Perivascular/Smooth muscle"
  keys_from_base <- names(base_palette)
  
  
  # --- 2. 建立 干净键 -> 颜色 的映射 ---
  
  # 创建一个以干净键为名字，颜色为值的向量。
  # 这一步非常重要，因为它标准化了我们的查找表。
  color_map <- setNames(unname(base_palette), keys_from_base)
  
  
  # --- 3. 进行匹配 ---
  
  # 使用我们从 counted_labels 中提取的干净键，去 color_map 中查找颜色
  matching_colors <- color_map[keys_from_counted]
  
  # 创建一个新的、以“带计数的标签”为名字的最终调色板
  new_palette <- setNames(matching_colors, levels_counted)
  
  
  # --- 4. 健壮性检查 (与之前相同) ---
  
  unmatched_labels <- names(new_palette[is.na(new_palette)])
  if (length(unmatched_labels) > 0) {
    warning("Some labels did not have a matching color in 'base_palette' and will be gray. Missing: ", 
            paste(unmatched_labels, collapse = ", "))
    new_palette[is.na(new_palette)] <- "grey80" # 为未匹配的项提供一个默认颜色
  }
  
  return(new_palette)
}



#' @title Update or Add Metadata Columns to a Seurat Object
#' @description Safely updates specified metadata columns in a Seurat object 
#'              from an external data frame, ensuring cell barcodes are perfectly aligned.
#'
#' @param seurat_obj A Seurat object to be updated.
#' @param new_metadata A data frame containing the new metadata. Its rownames 
#'                     must be cell barcodes that correspond to the Seurat object.
#' @param cols_to_update A character vector of column names from `new_metadata` 
#'                       that you want to update or add to the Seurat object's metadata. 
#'                       If NULL (default), all columns from `new_metadata` will be used.
#'
#' @return An updated Seurat object with the new metadata.
#' @export
#'
#' @examples
#' # object <- UpdateSeuratMetadata(
#' #   seurat_obj = object,
#' #   new_metadata = RTCDmetadata,
#' #   cols_to_update = c("RCTD_spot_class", "RCTD_first", "RCTD_second")
#' # )
UpdateSeuratMetadata <- function(seurat_obj, new_metadata, cols_to_update = NULL) {
  
  # --- 1. 参数验证 ---
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  if (!is.data.frame(new_metadata)) {
    stop("`new_metadata` must be a data frame.")
  }
  
  cat(">>> Starting metadata update process...\n")
  
  # --- 2. 对齐和验证 Barcodes ---
  barcodes_seurat <- rownames(seurat_obj@meta.data)
  barcodes_new <- rownames(new_metadata)
  
  # 检查barcode集合是否完全一致
  if (length(barcodes_seurat) != length(barcodes_new) || 
      length(intersect(barcodes_seurat, barcodes_new)) != length(barcodes_seurat)) {
    
    stop("Cell barcode sets do not match between the Seurat object and the new metadata. Update aborted.")
  }
  
  # 检查顺序是否一致，如果不一致则重排
  if (!identical(barcodes_seurat, barcodes_new)) {
    cat("... Barcode order differs. Reordering new metadata to match Seurat object...\n")
    new_metadata <- new_metadata[barcodes_seurat, , drop = FALSE]
  }
  
  cat("✅ Barcodes are perfectly aligned.\n")
  
  # --- 3. 确定要更新的列 ---
  if (is.null(cols_to_update)) {
    # 如果用户没有指定，则使用new_metadata中的所有列
    cols_to_update <- colnames(new_metadata)
    cat("... `cols_to_update` not specified. Using all", length(cols_to_update), "columns from new metadata.\n")
  } else {
    # 检查指定的列是否存在于new_metadata中
    missing_cols <- setdiff(cols_to_update, colnames(new_metadata))
    if (length(missing_cols) > 0) {
      stop("The following columns specified in `cols_to_update` were not found in `new_metadata`: ", 
           paste(missing_cols, collapse = ", "))
    }
  }
  
  # --- 4. 执行替换/添加操作 ---
  cat(">>> Updating", length(cols_to_update), "metadata columns...\n")
  
  # 直接用列名索引进行赋值，这会自动处理更新现有列和添加新列
  seurat_obj@meta.data[, cols_to_update] <- new_metadata[, cols_to_update]
  
  cat("✅ Metadata update complete!\n")
  
  # --- 5. 返回更新后的Seurat对象 ---
  return(seurat_obj)
}
AddCountedLabels <- function(seurat_obj, col_to_count, new_col_name,
                             drop_empty_levels = TRUE, big.mark = NULL) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  if (!col_to_count %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", col_to_count, "' not found in Seurat object metadata.")
  }
  
  # --- v5-safe: read as a VECTOR, not a data.frame ---
  x <- seurat_obj@meta.data[[col_to_count]]
  
  # drop unused levels (optional but recommended)
  if (is.factor(x) && drop_empty_levels) x <- droplevels(x)
  
  # ensure factor, capture levels
  f <- factor(x)
  lvl <- levels(f)
  
  # counts aligned to all levels (no NA), even if some have 0
  counts <- table(factor(x, levels = lvl))
  cnt_vec <- as.integer(counts)
  
  # optional nice formatting for counts: big.mark = "," or " "
  if (!is.null(big.mark)) {
    cnt_lab <- format(cnt_vec, big.mark = big.mark, trim = TRUE, scientific = FALSE)
  } else {
    cnt_lab <- as.character(cnt_vec)
  }
  
  counted_levels <- sprintf("%s (%s)", lvl, cnt_lab)
  map <- setNames(counted_levels, lvl)
  counted_values <- map[as.character(f)]
  
  # write back (v5-safe)
  seurat_obj@meta.data[[new_col_name]] <- factor(counted_values, levels = counted_levels)
  return(seurat_obj)
}

#' @title Curate and Consolidate Top-Level Cell Type Annotations
#' @description Merge heterogeneous labels into a consistent top-level annotation,
#'              and order the final factor levels by a biological hierarchy.
#'
#' @param object A Seurat object.
#' @param source_col Metadata column to start from. If NULL, auto-detects from a priority list.
#' @param write_to Name of the new metadata column to create.
#' @param level_order Character vector specifying the desired biological order.
#'                    If NULL, a PDAC-specific default is used. Levels not present
#'                    in the data are ignored; unseen levels are appended alphabetically.
#' @param min_count_to_keep Merge rare cell types into "Other (lowN)" if n < this threshold (0 to disable).
#' @param split_B_and_Plasma Logical, whether to keep B cells and Plasma cells separate.
#' @param split_ductal Logical, whether to keep Ductal, Ductal/ADM, Ductal(PanIN-like) separate.
#' @param drop_qc_from_plots Logical flag returned for downstream plotting.
#'
#' @return list(object, palette, table, drop_qc_from_plots)
CurateCelltypesTop <- function(object,
                               source_col = NULL,
                               write_to   = "cell_type_merged",
                               level_order = NULL,
                               min_count_to_keep = 0,
                               split_B_and_Plasma = TRUE, # 默认拆分 B 和 Plasma
                               split_ductal = TRUE,       # 默认拆分 Ductal 类型
                               drop_qc_from_plots = TRUE) {
  
  stopifnot(inherits(object, "Seurat"))
  
  # --- 1. 自动检测或使用指定的来源列 ---
  if (is.null(source_col)) {
    # 优先使用经过免疫细胞校正的列，其次是Seurat注释，最后是RCTD
    source_col <- dplyr::case_when(
      "cell_type_merged_immuneFirst" %in% colnames(object@meta.data) ~ "cell_type_merged_immuneFirst",
      "Seurat_L1_label"              %in% colnames(object@meta.data) ~ "Seurat_L1_label",
      "RCTD_first"                   %in% colnames(object@meta.data) ~ "RCTD_first",
      TRUE ~ "not_found"
    )
    if (source_col == "not_found") {
      stop("Could not find a valid source column. Please provide one via 'source_col'.")
    }
  }
  
  # --- 2. 强大的同义词归一化 (Synonym Normalization) ---
  # 这是本次修改的核心部分
  
  src_vector <- as.character(object@meta.data[[source_col]])
  
  # 创建一个全面的同义词列表。键是各种可能的写法，值是标准化的名称。
  synonym_map <- c(
    # B细胞谱系
    "B cell" = "B cells", "B_cells" = "B cells",
    # 浆细胞谱系
    "Plasma cell" = "Plasma cells", "Plasma" = "Plasma cells", "Plasmablast" = "Plasma cells",
    # T细胞谱系
    "T cell" = "T cells", "T_cells" = "T cells",
    # 髓系细胞谱系
    "Myeloid" = "Myeloid", "Macrophage" = "Myeloid", "cDC" = "Myeloid", "Mast" = "Myeloid", # Mast可单独或合并
    # 导管细胞谱系 (关键！)
    "Ductal-mucinous (PanIN-like)" = "Ductal(PanIN-like)", "Metaplastic ductal (PanIN-like)" = "Ductal(PanIN-like)",
    "PanIN-like" = "Ductal(PanIN-like)", "ADM" = "Ductal/ADM",
    # 内皮细胞
    "Vascular" = "Endothelial", "NascentEndothelial" = "Endothelial",
    # 周细胞/平滑肌
    "SmoothMuscle" = "Perivascular/SMC", "Perivascular/Smooth muscle" = "Perivascular/SMC",
    # 神经/雪旺细胞
    "Schwann/Glia" = "Neural/Schwann", "Schwann" = "Neural/Schwann",
    # 其他
    "Tumor" = "Tumor/Epithelial", "Fibroblast" = "CAF"
  )
  
  # 使用 recode 进行高效替换
  normalized_vector <- dplyr::recode(src_vector, !!!synonym_map)
  
  # --- 3. 基于逻辑规则的最终类别合并 (Final Category Merging) ---
  
  # a) Myeloid 细胞合并
  myeloid_pattern <- "Myeloid|Macrophage|cDC"
  normalized_vector[grepl(myeloid_pattern, normalized_vector, ignore.case = TRUE)] <- "Myeloid (Mono/Macro/DC)"
  
  # b) 导管细胞按需拆分或合并
  if (!isTRUE(split_ductal)) {
    ductal_pattern <- "Ductal|ADM|PanIN"
    normalized_vector[grepl(ductal_pattern, normalized_vector, ignore.case = TRUE)] <- "Ductal (Combined)"
  }
  
  # c) B细胞和浆细胞按需拆分或合并
  if (!isTRUE(split_B_and_Plasma)) {
    b_plasma_pattern <- "B cells|Plasma cells"
    normalized_vector[grepl(b_plasma_pattern, normalized_vector, ignore.case = TRUE)] <- "B/Plasma cells"
  }
  
  # --- 4. 处理低丰度细胞类型 ---
  if (min_count_to_keep > 0) {
    tbl <- table(normalized_vector)
    low_abundance_types <- names(tbl)[tbl < min_count_to_keep]
    if (length(low_abundance_types) > 0) {
      normalized_vector[normalized_vector %in% low_abundance_types] <- "Other (low abundance)"
    }
  }
  
  # --- 5. 应用生物学顺序和创建因子 ---
  if (is.null(level_order)) {
    level_order <- c(
      # 肿瘤与上皮
      "Tumor/Epithelial", "Ductal(PanIN-like)", "Ductal/ADM", "Ductal", "Ductal (Combined)", "Acinar", "Endocrine",
      # 基质
      "CAF", "Endothelial", "Perivascular/SMC", "Neural/Schwann",
      # 免疫
      "Myeloid (Mono/Macro/DC)", "T cells", "B cells", "Plasma cells", "B/Plasma cells", "Mast",
      # 其他
      "Ambiguous/Doublet","LowQC/Ambiguous/Doublet", "Other (low abundance)"
    )
  }
  
  current_levels <- unique(normalized_vector)
  final_order <- level_order[level_order %in% current_levels]
  missing_from_order <- setdiff(current_levels, final_order)
  if (length(missing_from_order) > 0) {
    final_order <- c(final_order, sort(missing_from_order))
  }
  
  merged_factor <- factor(normalized_vector, levels = final_order)
  
  # --- 6. 生成匹配的调色板 ---
  pal_top <- c(
    "Tumor/Epithelial" = "#d73027", "Ductal(PanIN-like)" = "#6a51a3",  # 最深
    "Ductal/ADM"         = "#7570b3",  # 中等
    "Ductal"             = "#cbc9e2",  # 最浅
    
    
    # "Ductal(PanIN-like)" = "#984ea3", "Ductal/ADM" = "#e78ac3",
    # "Ductal" = "#7570b3", "Ductal (Combined)" = "#7570b3",
    "Acinar" = "#4daf4a", "Endocrine" = "#ffd700",
    "CAF" = "#a65628", "Endothelial" = "#74add1", "Perivascular/SMC" = "#ff7f00", "Neural/Schwann" = "black",
    "Myeloid (Mono/Macro/DC)" = "#483D8B", "T cells" = "#377eb8", "B cells" = "#4292c6",
    "Plasma cells" = "#6baed6", "B/Plasma cells" = "#4292c6", "Mast" = "#e41a1c",
    "Ambiguous/Doublet" = "grey70", "LowQC/Ambiguous/Doublet"= "grey70","Other (low abundance)" = "grey50"
  )
  
  final_palette <- pal_top[final_order]
  names(final_palette) <- final_order
  
  # --- 7. 将结果写回对象并返回 ---
  object[[write_to]] <- merged_factor
  summary_table <- as.data.frame(table(object[[write_to]]))
  colnames(summary_table) <- c("cell_type", "count")
  
  invisible(list(object = object,
                 palette = final_palette,
                 table = summary_table,
                 drop_qc_from_plots = drop_qc_from_plots))
}


theme_slide_light <- function(
    base_size      = 11,
    panel_fill     = "#f7f9fb",
    panel_alpha    = 0.8,
    border_color   = "#d9dfe7",
    legend_position= "right",
    glow           = TRUE,
    glow_fill1     = "#ffffff",
    glow_alpha1    = 0.10,
    glow_fill2     = "#eff5ff",
    glow_alpha2    = 0.05,
    # 新增：legend 尺寸控制
    legend_text_size  = 8,
    legend_title_size = 9,
    legend_key_hw     = c(0.35, 0.45)  # c(height_cm, width_cm)
){
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("请安装 'scales' 包: install.packages('scales')")
  }
  thm <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      # 背景与边框
      panel.background = ggplot2::element_rect(
        fill   = scales::alpha(panel_fill, panel_alpha),
        colour = NA
      ),
      panel.border = ggplot2::element_rect(
        colour = border_color, fill = NA, linewidth = 0.8
      ),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      
      # 无网格 & 无坐标元素
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line        = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_blank(),
      axis.title       = ggplot2::element_blank(),
      
      # 图例
      legend.position      = legend_position,
      legend.background    = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.key           = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.box.background= ggplot2::element_blank(),
      # ↓↓↓ 新增：文字和键尺寸 ↓↓↓
      legend.text  = ggplot2::element_text(size = legend_text_size),
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.key.height = grid::unit(legend_key_hw[1], "cm"),
      legend.key.width  = grid::unit(legend_key_hw[2], "cm"),
      strip.background = ggplot2::element_blank()
    )
  
  if (isTRUE(glow)) {
    return(list(
      thm,
      ggplot2::annotate("rect",
                        xmin = -Inf, xmax =  Inf, ymin = -Inf, ymax =  Inf,
                        fill = glow_fill1, alpha = glow_alpha1
      ),
      ggplot2::annotate("rect",
                        xmin = -Inf, xmax =  Inf, ymin = -Inf, ymax =  Inf,
                        fill = glow_fill2, alpha = glow_alpha2
      )
    ))
  } else {
    return(thm)
  }
}

# ----- Clean journal theme: no gridlines, no background, thin axes -----
theme_ns_clean <- function(base_size = 11, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold", margin = ggplot2::margin(b = 6)),
      axis.title.x       = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y       = ggplot2::element_text(margin = ggplot2::margin(r = 6)),
      panel.grid.major   = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.background   = ggplot2::element_rect(fill = NA, colour = NA),
      plot.background    = ggplot2::element_rect(fill = NA, colour = NA),
      axis.line          = ggplot2::element_line(linewidth = 0.4, colour = "black"),
      axis.ticks         = ggplot2::element_line(linewidth = 0.4, colour = "black"),
      axis.text.x        = ggplot2::element_text(size = base_size - 1, angle = 0, hjust = 0.5),
      legend.title       = ggplot2::element_text(face = "bold"),
      legend.background  = ggplot2::element_blank(),
      legend.key         = ggplot2::element_blank(),
      plot.margin        = grid::unit(c(6, 10, 6, 6), "pt")
    )
}


BandLinePlot <- function(
    object,
    base_colors,
    band_col   = "PERI_band",
    type_col   = "Seurat_L1_label",
    bands      = c("0–50 µm","50–100 µm","100–150 µm",">150 µm"),
    tumor_type = "Tumor/Epithelial",
    title      = "Non-tumor composition trends across periphery",
    ylab       = "Fraction within TME",
    
    types_to_plot = NULL,        # optional whitelist BEFORE ranking
    top_n_types   = NULL,        # NEW: select top-N by mean(frac) across bands; NULL = all
    keep_highlight = TRUE,       # NEW: always keep highlight_types even if not in top-N
    
    line_size  = 1.1,
    point_size = 2.5,
    legend_ncol = 1,
    legend_y_spacing = 0.5,     # << from 0.3 to 0.5
    legend_key_h = 0.6,
    legend_key_w = 0.8,
    legend_text_size = 12,      # << NEW
    legend_title_size = 13,     # << NEW
    highlight_types = NULL,
    highlight_size_mult = 1.5,
    non_highlight_alpha = 0.4,
    width  = 7, height = 5, dpi = 300,
    out_file = NULL
){
  # --- Data prep ---
  normalize_band <- function(x){ gsub("-", "–", x, fixed=TRUE) }
  
  df <- object@meta.data %>%
    tibble::as_tibble(rownames = "cell_id") %>%
    dplyr::transmute(
      band = normalize_band(as.character(.data[[band_col]])),
      type = as.character(.data[[type_col]])
    ) %>%
    dplyr::filter(
      grepl("µm", band),
      type != tumor_type
    )
  
  if (!is.null(types_to_plot)) {
    df <- df %>% dplyr::filter(type %in% types_to_plot)
  }
  
  df <- df %>% dplyr::mutate(band = factor(band, levels = bands[bands %in% unique(df$band)]))
  
  # Count per band/type and compute within-band fractions
  tab_ct <- df %>%
    dplyr::count(band, type, name = "n") %>%
    dplyr::group_by(band) %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::ungroup()
  
  # --- NEW: Top-N selection by mean(frac) across bands ------------------------
  if (!is.null(top_n_types)) {
    top_n_types <- max(1L, as.integer(top_n_types))  # safety clamp
    rank_tbl <- tab_ct %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(mean_frac = mean(frac), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(mean_frac))
    
    keep_types <- utils::head(rank_tbl$type, top_n_types)
    
    if (!is.null(highlight_types) && keep_highlight) {
      keep_types <- union(keep_types, intersect(highlight_types, unique(tab_ct$type)))
    }
    
    tab_ct <- tab_ct %>% dplyr::filter(type %in% keep_types)
  }
  
  # Order legend by mean(frac) of the (possibly filtered) set
  type_order <- tab_ct %>%
    dplyr::group_by(type) %>%
    dplyr::summarise(mean_frac = mean(frac), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_frac)) %>%
    dplyr::pull(type)
  tab_ct$type <- factor(tab_ct$type, levels = type_order)
  
  # Palette: extend base_colors if needed
  pal <- base_colors
  missing <- setdiff(levels(tab_ct$type), names(pal))
  if (length(missing) > 0) {
    hues <- seq(15, 375, length.out = length(missing) + 1)[1:length(missing)]
    extra <- grDevices::hcl(h = hues, c = 60, l = 65)
    names(extra) <- missing
    pal <- c(pal, extra)
  }
  pal_used <- pal[levels(tab_ct$type)]
  
  # Build plot (highlight-aware)
  p <- ggplot2::ggplot(tab_ct, ggplot2::aes(x = band, y = frac, group = type, color = type))
  
  if (!is.null(highlight_types)) {
    df_highlight <- tab_ct %>% dplyr::filter(type %in% highlight_types)
    df_non_highlight <- tab_ct %>% dplyr::filter(!type %in% highlight_types)
    
    p <- p +
      ggplot2::geom_line(data = df_non_highlight, linewidth = line_size, alpha = non_highlight_alpha) +
      ggplot2::geom_point(data = df_non_highlight, size = point_size, alpha = non_highlight_alpha) +
      ggplot2::geom_line(data = df_highlight, linewidth = line_size * highlight_size_mult) +
      ggplot2::geom_point(data = df_highlight, size = point_size * highlight_size_mult)
  } else {
    p <- p +
      ggplot2::geom_line(linewidth = line_size) +
      ggplot2::geom_point(size = point_size)
  }
  
  p <- p +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                limits = c(0, NA),
                                expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::scale_color_manual(values = pal_used, drop = FALSE, name = "Cell type") +
    ggplot2::labs(x = "Distance from tumor", y = ylab, title = title) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.line = ggplot2::element_line(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.4),
      axis.text.x  = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
      legend.position = "right",
      legend.title    = ggplot2::element_text(face = "bold", size = legend_title_size),  # << MOD
      legend.text     = ggplot2::element_text(size = legend_text_size),                   # << MOD
      legend.key.height = grid::unit(legend_key_h, "lines"),
      legend.key.width  = grid::unit(legend_key_w, "lines"),
      legend.spacing.y  = grid::unit(legend_y_spacing, "lines")
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE))
  
  
  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, p, width = width, height = height, dpi = dpi)
  }
  
  invisible(list(plot = p, table = tab_ct, palette = pal_used))
}





# =========================
# PlotspotPieChart (patched)
# ============================================================
# Nature-style Pie/Donut with ggforce explode + smart labeling
# ============================================================
PlotspotPieChart_Nature <- function(
    object,                           # Seurat object or data.frame (metadata)
    celltype_col,                     # e.g. "Seurat_L2_refined"
    subset_expr = NULL,               # e.g. "InvasiveFront_flag == TRUE"
    drop_na = TRUE,                   # drop NA in celltype_col
    order = c("freq_desc","freq_asc","alphabetical","custom"),
    order_levels = NULL,              # when order="custom"
    palette = NULL,                   # named vector or plain vector; will extend
    r = 2,                            # outer radius
    r0 = 0,                           # inner radius (0=pie; 1~1.2=donut)
    explode_by = NULL,                # character vector of groups to explode
    explode_top = NULL,               # integer: explode top-N groups
    explode_amt = 0.15,               # explode distance (0~0.3 typical)
    inside_cutoff = 8,                # % ≥ this ⇒ annotate on pie (inside ring)
    annotate_names = NULL,            # force-annotate these groups on pie
    label_format = c("name_percent","percent","name","both"),
    legend_for_small = TRUE,          # send small slices to legend
    legend_ncol = 1,
    title = NULL, subtitle = NULL, caption = NULL,
    show_ring_for = NULL,             # character vec: draw outer ring for emphasis
    ring_expand = 0.18,               # how much ring expands (outer highlight)
    width = 7, height = 5, dpi = 300,
    out_file = NULL,                  # file path for combined plot (patchwork)
    return_what = c("combined","pie","legend","table") # what to return primarily
) {
  # -------- deps check --------
  pkgs_needed <- c("ggplot2","dplyr","rlang","ggforce","scales")
  for (p in pkgs_needed) if (!requireNamespace(p, quietly = TRUE))
    stop("Please install package: ", p)
  patchwork_ok <- requireNamespace("patchwork", quietly = TRUE)
  ggpubr_ok    <- requireNamespace("ggpubr", quietly = TRUE)
  export_ok    <- requireNamespace("export", quietly = TRUE)
  
  `%>%` <- dplyr::`%>%`
  
  # -------- theme (Nature-ish minimal) --------
  theme_natureish <- function(base_size = 12) {
    ggplot2::theme_void(base_size = base_size) %+replace%
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face="bold", size=base_size+2, hjust=0.5),
        plot.subtitle = ggplot2::element_text(size=base_size, hjust=0.5, margin=ggplot2::margin(b=6)),
        plot.caption  = ggplot2::element_text(size=base_size-1, color="grey30", hjust=1),
        legend.title  = ggplot2::element_blank(),
        legend.text   = ggplot2::element_text(size=base_size),
        plot.margin   = ggplot2::margin(6, 10, 6, 6)
      )
  }
  
  # -------- ingest metadata --------
  if ("Seurat" %in% class(object)) {
    meta <- object@meta.data
  } else if (is.data.frame(object)) {
    meta <- object
  } else stop("`object` must be a Seurat object or a data.frame.")
  
  if (!celltype_col %in% colnames(meta)) {
    stop("Column `", celltype_col, "` not found.")
  }
  if (!is.null(subset_expr)) {
    meta <- dplyr::filter(meta, !!rlang::parse_expr(subset_expr))
  }
  if (drop_na) {
    meta <- dplyr::filter(meta, !is.na(.data[[celltype_col]]))
  }
  if (nrow(meta) == 0) stop("No rows after filtering.")
  
  # coerce to character (avoid list-column surprises)
  if (is.list(meta[[celltype_col]])) {
    meta[[celltype_col]] <- vapply(meta[[celltype_col]], function(x) paste(as.character(x), collapse="; "), character(1))
  }
  meta[[celltype_col]] <- as.character(meta[[celltype_col]])
  
  # -------- composition table --------
  comp <- meta %>%
    dplyr::count(.data[[celltype_col]], name = "n", sort = TRUE) %>%
    dplyr::rename(name = !!rlang::sym(celltype_col)) %>%
    dplyr::mutate(percentage = n / sum(n) * 100)
  
  # -------- ordering --------
  order <- match.arg(order)
  if (order == "freq_desc") comp <- dplyr::arrange(comp, dplyr::desc(percentage))
  if (order == "freq_asc")  comp <- dplyr::arrange(comp, percentage)
  if (order == "alphabetical") comp <- dplyr::arrange(comp, name)
  if (order == "custom") {
    if (is.null(order_levels)) stop("When order='custom', provide `order_levels`.")
    comp <- comp %>% dplyr::mutate(.ord = match(name, order_levels)) %>%
      dplyr::arrange(.ord) %>% dplyr::select(-.ord)
  }
  comp$name <- factor(comp$name, levels = unique(comp$name))
  
  # -------- explode selection --------
  # build explode (focus) vector per slice (same order as comp)
  focus <- rep(0, nrow(comp))
  names(focus) <- as.character(comp$name)
  if (!is.null(explode_by)) {
    focus[names(focus) %in% explode_by] <- explode_amt
  }
  if (!is.null(explode_top)) {
    top_ids <- seq_len(min(explode_top, nrow(comp)))
    focus[top_ids] <- explode_amt
  }
  
  # -------- angles for arcs (needed for ring/labels/legend text) --------
  comp2 <- comp %>%
    dplyr::mutate(
      end_angle   = 2*pi*cumsum(percentage)/100,
      start_angle = dplyr::lag(end_angle, default = 0),
      mid_angle   = 0.5*(start_angle + end_angle),
      legend_lab  = paste0(name, " (", scales::number(percentage, accuracy = 0.1), "%)"),
      focus       = focus[as.character(name)]
    )
  
  # -------- color palette --------
  # Okabe-Ito default; extend as needed
  pal_default <- c("#E69F00","#56B4E9","#009E73","#F0E442",
                   "#0072B2","#D55E00","#CC79A7","#999999")
  labels_present <- as.character(comp2$name)
  if (is.null(palette)) {
    if (length(labels_present) > length(pal_default)) {
      palette <- grDevices::colorRampPalette(pal_default)(length(labels_present))
    } else palette <- pal_default[seq_along(labels_present)]
    names(palette) <- labels_present
  } else {
    if (is.null(names(palette))) {
      if (length(palette) < length(labels_present))
        palette <- grDevices::colorRampPalette(palette)(length(labels_present))
      names(palette) <- labels_present
    } else {
      # reorder and auto-extend if missing
      miss <- setdiff(labels_present, names(palette))
      if (length(miss) > 0) {
        add_cols <- grDevices::colorRampPalette(pal_default)(length(miss))
        names(add_cols) <- miss
        palette <- c(palette, add_cols)
      }
      palette <- palette[labels_present]
    }
  }
  
  # -------- which slices get on-pie annotations --------
  label_format <- match.arg(label_format)
  fmt <- switch(label_format,
                "name_percent" = function(n, p) paste0(n, "\n(", scales::number(p, accuracy=0.1), "%)"),
                "percent"      = function(n, p) paste0(scales::number(p, accuracy=0.1), "%"),
                "name"         = function(n, p) n,
                "both"         = function(n, p) paste0(n, "  ", scales::number(p, accuracy=0.1), "%"))
  
  force_names <- if (is.null(annotate_names)) character(0) else annotate_names
  comp2 <- comp2 %>%
    dplyr::mutate(
      on_pie = (percentage >= inside_cutoff) | (name %in% force_names),
      label_text = fmt(as.character(name), percentage)
    )
  
  # -------- base pie with ggforce (no coord_polar; true geometry) --------
  pie <- ggplot2::ggplot() +
    ggforce::geom_arc_bar(
      data = comp2,
      ggplot2::aes(
        x0 = 0, y0 = 0, r0 = r0, r = r,
        amount = percentage, fill = name, color = name,
        explode = focus
      ),
      stat = "pie", show.legend = FALSE, linewidth = 0.3
    )
  
  # optional outer ring highlight for specific groups
  if (!is.null(show_ring_for)) {
    ring_df <- comp2 %>% dplyr::filter(name %in% show_ring_for)
    if (nrow(ring_df) > 0) {
      pie <- pie +
        ggforce::geom_arc(
          data = ring_df,
          ggplot2::aes(
            x0 = 0, y0 = 0, r = r + ring_expand,
            start = start_angle, end = end_angle
          ),
          linewidth = 0.9, color = "black", inherit.aes = FALSE
        )
    }
  }
  
  # add a thin outline for the whole pie (optional, subtle)
  pie <- pie + ggplot2::coord_fixed()
  
  # fill/color scales
  pie <- pie +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_color_manual(values = palette, drop = FALSE)
  
  # -------- annotate large slices ON the pie --------
  # annotate near the middle radius between r0 and r for visibility
  r_mid <- if (r0 > 0) (r0 + r)/2 else 0.65*r
  ann_df <- comp2 %>% dplyr::filter(on_pie)
  if (nrow(ann_df) > 0) {
    ann_df <- ann_df %>%
      dplyr::mutate(
        x = r_mid * sin(mid_angle),
        y = r_mid * cos(mid_angle)
      )
    pie <- pie +
      ggplot2::annotate(
        "text",
        x = ann_df$x, y = ann_df$y,
        label = ann_df$label_text,
        size = 4, fontface = "plain", lineheight = 1.05
      )
  }
  
  # -------- build legend that contains name + percentage for small slices --------
  small_df <- comp2 %>% dplyr::filter(!on_pie)
  leg_plot <- NULL
  
  if (legend_for_small && nrow(small_df) > 0) {
    # Build a dummy plot to extract a clean legend with our labels+colors
    dummy <- ggplot2::ggplot(small_df, ggplot2::aes(x = 1, y = legend_lab, fill = legend_lab)) +
      ggplot2::geom_tile(ggplot2::aes(width = 0.5, height = 0.5)) +
      ggplot2::scale_fill_manual(values = setNames(palette[as.character(small_df$name)],
                                                   small_df$legend_lab), drop = FALSE) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        legend.position = "right",
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 12),
        legend.key.width  = ggplot2::unit(0.45,"cm"),
        legend.key.height = ggplot2::unit(0.45,"cm"),
        legend.spacing.y  = ggplot2::unit(2.5, "mm")
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol))
    if (ggpubr_ok) {
      leg_plot <- ggpubr::as_ggplot(ggpubr::get_legend(dummy))
    } else {
      # Fallback: keep normal legend on the pie instead
      pie <- pie + ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol))
    }
  }
  
  # -------- titles & theme --------
  pie <- pie +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption) +
    theme_natureish(12)
  
  # -------- combine with patchwork if we have an external legend --------
  combined <- pie
  if (!is.null(leg_plot) && patchwork_ok) {
    combined <- pie + leg_plot + patchwork::plot_layout(widths = c(1.25, 1))
  }
  
  # -------- export combined if path provided --------
  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, plot = combined, width = width, height = height,
                    dpi = dpi, limitsize = FALSE)
  }
  
  # -------- returns --------
  out <- list(
    pie      = pie,
    legend   = leg_plot,
    combined = combined,
    table    = as.data.frame(comp2[, c("name","n","percentage")])
  )
  main <- match.arg(return_what)
  return(out[[main]])
}



#' 绘制 Step 2 分析结果的可视化图表
#'
#' @param results_table 来自 `step2_sender_signals_and_forest` 函数的输出表格 (`$table`)。
#' @param plot_type 要绘制的图形类型，可以是 "volcano" (默认) 或 "dual_forest"。
#' @param top_n 对于 "dual_forest"，显示 Top N 个差异最显著的基因。
#' @param p_val_cutoff 在火山图中标记为显著的 p 值阈值。
#' @param logfc_cutoff 在火山图中标记为显著的 log2FC 阈值。
#'
#' @return 一个 ggplot 对象。
#'
plot_step2_results <- function(
    results_table, 
    plot_type = "volcano",
    top_n = 15,
    p_val_cutoff = 0.05,
    logfc_cutoff = 0.25
) {
  
  # --- 依赖检查 ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'")
  if (!requireNamespace("ggrepel", quietly = TRUE) & plot_type == "volcano") {
    warning("For better labels, please install 'ggrepel'.")
  }
  
  # 移除无效行 (例如之前 Classical 的 NA 行) 和 SIGNATURE 行
  plot_data <- results_table %>%
    dplyr::filter(!is.na(periphery_of), feature != "SIGNATURE")
  
  # ============================================
  #  方案A: 火山-棒棒糖图
  # ============================================
  if (plot_type == "volcano") {
    
    plot_data_volcano <- plot_data %>%
      mutate(
        logP = -log10(pmax(p_value, 1e-10)), # 避免 log10(0)
        label = feature,
        significance_group = case_when(
          p_value < p_val_cutoff & log2_near_over_far > logfc_cutoff ~ "Upregulated",
          p_value < p_val_cutoff & log2_near_over_far < -logfc_cutoff ~ "Downregulated",
          TRUE ~ "Not Significant"
        )
      )
    
    p <- ggplot(plot_data_volcano, aes(x = log2_near_over_far, y = logP)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(p_val_cutoff), linetype = "dotted", color = "red") +
      geom_segment(aes(xend = log2_near_over_far, yend = 0, color = sender),
                   linewidth = 0.8, alpha = 0.6) +
      geom_point(aes(color = sender, shape = significance_group), size = 4, alpha = 0.9) +
      theme_bw(base_size = 14) +
      labs(
        title = "Gene Expression Changes in Stromal Cells at Invasive Front",
        subtitle = "Near vs. Far around Basal-like tumors",
        x = "log2(Mean Expression Near / Far)",
        y = "-log10(p-value)",
        color = "Sender Cell Type",
        shape = "Significance"
      ) +
      scale_color_manual(values = c("myCAF (COL11A1+/desmoplastic)" = "#E41A1C", "SPP1 TAM" = "#377EB8"),
                         guide = guide_legend(override.aes = list(linewidth = 0))) + # 隐藏图例中的线
      scale_shape_manual(values = c("Upregulated" = 16, "Not Significant" = 1, "Downregulated" = 6)) +
      scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
      coord_cartesian(xlim = c(-1.5, 1.5))
    
    # 只有安装了 ggrepel 才添加标签，避免报错
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = subset(plot_data_volcano, p_value < p_val_cutoff & abs(log2_near_over_far) > 0.5),
        aes(label = label),
        box.padding = 0.5,
        max.overlaps = Inf,
        size = 3.5
      )
    }
    
    return(p)
  }
  
  # ============================================
  #  方案B: 双指标森林图
  # ============================================
  if (plot_type == "dual_forest") {
    
    plot_data_forest <- plot_data %>%
      dplyr::arrange(p_value) %>%
      dplyr::slice_head(n = top_n) %>%
      tidyr::pivot_longer(
        cols = c("log2_near_over_far", "log2DR"),
        names_to = "metric_type",
        values_to = "value"
      ) %>%
      mutate(
        label = feature,
        metric_label = ifelse(metric_type == "log2_near_over_far", "Mean Expression (log2FC)", "Detection Rate (log2DR)")
      )
    
    p <- ggplot(plot_data_forest, aes(x = value, y = forcats::fct_reorder(label, value, .fun = max))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(aes(color = sender), size = 4) +
      facet_wrap(~ metric_label) +
      theme_bw(base_size = 12) +
      labs(
        title = paste("Top", top_n, "DEGs in Stromal Cells at Invasive Front"),
        subtitle = "Change in Mean Expression vs. Detection Rate",
        x = "log2(Near / Far)",
        y = "Gene",
        color = "Sender Cell Type"
      ) +
      scale_color_manual(values = c("myCAF (COL11A1+/desmoplastic)" = "#E41A1C", "SPP1 TAM" = "#377EB8"))
    
    return(p)
  }
  
  stop("Invalid 'plot_type'. Please choose 'volcano' or 'dual_forest'.")
}

#' @title 更新母集Seurat对象的细胞注释 (版本2 - 修正因子问题)
#' @description 将子集对象中新的细胞注释回填到母集对象的指定列中。
#'              此版本能自动处理目标列为因子(factor)类型的情况。
#'
#' @param object_master 母集Seurat对象。
#' @param object_subset 子集Seurat对象。
#' @param source_col_subset 子集对象中包含新注释的列名。
#' @param target_col_master 母集对象中需要被更新的目标列名。
#'
#' @return 返回一个更新了注释的母集Seurat对象。
#'
update_master_annotations <- function(object_master,
                                      object_subset,
                                      source_col_subset,
                                      target_col_master) {
  # --- 1. 安全性检查 ---
  if (!inherits(object_master, "Seurat") || !inherits(object_subset, "Seurat")) {
    stop("输入的对象必须是Seurat对象。")
  }
  if (!source_col_subset %in% colnames(object_subset@meta.data)) {
    stop("来源列 '", source_col_subset, "' 在子集对象中不存在。")
  }
  if (!target_col_master %in% colnames(object_master@meta.data)) {
    stop("目标列 '", target_col_master, "' 在母集对象中不存在。")
  }
  
  # --- 2. 创建一个“注释查找图” ---
  meta_subset <- object_subset@meta.data
  annotation_map <- setNames(as.character(meta_subset[[source_col_subset]]), rownames(meta_subset))
  
  # --- 3. 找到母集中需要更新的细胞并执行更新 ---
  master_cells <- rownames(object_master@meta.data)
  cells_to_update <- intersect(master_cells, names(annotation_map))
  
  if (length(cells_to_update) == 0) {
    warning("在母集和子集对象中没有找到匹配的细胞条码，返回原始母集对象。")
    return(object_master)
  }
  
  # 【关键修正】: 在赋值前，确保目标列是字符类型
  if (is.factor(object_master@meta.data[[target_col_master]])) {
    message(paste("目标列 '", target_col_master, "' 是因子类型，将转换为字符类型以接受新注释。"))
    object_master@meta.data[[target_col_master]] <- as.character(object_master@meta.data[[target_col_master]])
  }
  
  new_labels <- annotation_map[cells_to_update]
  
  object_master@meta.data[cells_to_update, target_col_master] <- new_labels
  
  # --- 4. （可选）将更新后的列重新转换为因子类型 ---
  object_master@meta.data[[target_col_master]] <- factor(object_master@meta.data[[target_col_master]])
  
  # --- 5. 返回更新后的母集对象 ---
  message(paste("成功更新了", length(cells_to_update), "个细胞在 '", target_col_master, "' 列中的注释。"))
  return(object_master)
}


#' @title 为Seurat v5 创建一个真正独立的子集对象
#' @description 解决Seurat v5中 subset() 函数的惰性求值问题。
#'              该函数会提取指定细胞的原始counts数据和元数据，
#'              然后从头创建一个新的、干净的Seurat对象。
#'
#' @param object 原始的Seurat v5对象。
#' @param cells_to_keep 一个包含要保留的细胞barcode的字符向量。
#' @param assay 要提取counts数据的assay名称，默认为"Spatial"。
#'
#' @return 一个全新的、独立的Seurat对象。
#'
safe_subset_v5 <- function(object, cells_to_keep, assay = "Spatial") {
  
  message("为Seurat v5创建一个独立的子集对象...")
  
  # --- 1. 从原始对象中提取所需的数据 ---
  
  # 提取原始的counts矩阵
  counts_subset <- GetAssayData(object, assay = assay, layer = "counts")[, cells_to_keep, drop = FALSE]
  
  # 提取元数据
  meta_subset <- object@meta.data[cells_to_keep, , drop = FALSE]
  
  # --- 2. 使用提取的数据创建全新的Seurat对象 ---
  new_sobj <- CreateSeuratObject(
    counts = counts_subset,
    assay = assay,
    meta.data = meta_subset
  )
  
  message("新对象创建完成，包含了 ", ncol(new_sobj), " 个细胞。")
  return(new_sobj)
}



# ==== 独立图例生成函数（默认透明背景）====
StandaloneLegend <- function(
    palette,                         # 命名向量：名称=颜色
    shape = c("circle","square"),    # 圆/方可选（默认圆点）
    title = NULL,                    # 图例标题（NULL=不显示）
    ncol = 2,                        # 图例分成几列显示
    key_size = 6,                    # 点的大小
    text_size = 9,                   # 文字字号
    save_path = NULL,                # 如 "legend.png"/"legend.pdf"；NULL=不保存
    width_in = 3,                    # 保存宽度（英寸）
    dpi = 300                        # 保存分辨率
){
  stopifnot(is.character(palette), !is.null(names(palette)))
  shape <- match.arg(shape)
  
  suppressPackageStartupMessages({
    library(ggplot2); library(grid); library(gtable)
  })
  
  # 构造虚拟数据，仅用于生成图例
  df <- data.frame(lbl = factor(names(palette), levels = names(palette)), x = 1, y = 1)
  
  p <- ggplot(df, aes(x, y, color = lbl)) +
    geom_point(
      size  = key_size,
      shape = if (shape == "circle") 16 else 15,  # 16=实心圆；15=实心方块
      stroke = 0
    ) +
    scale_color_manual(values = palette, name = title) +
    guides(color = guide_legend(
      ncol = ncol, byrow = TRUE,
      override.aes = list(size = key_size, shape = if (shape == "circle") 16 else 15)
    )) +
    theme_void() +
    theme(
      legend.position      = c(0.5, 0.5),   # 面板中心
      legend.justification = c(0.5, 0.5),
      legend.key        = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA),
      plot.background   = element_rect(fill = NA, color = NA),
      legend.text       = element_text(size = text_size),
      legend.title      = element_text(size = text_size + 1, face = "bold")
    )
  
  # 提取图例为 grob
  gt <- ggplotGrob(p)
  legend_grob <- gtable::gtable_filter(gt, "guide-box", trim = TRUE)
  
  # 直接绘制在设备上（RStudio Viewer）
  grid::grid.newpage(); grid::grid.draw(legend_grob)
  
  # 可选保存：PNG 透明底 或 PDF
  if (!is.null(save_path)) {
    rows <- ceiling(length(palette) / ncol)
    height_in <- 0.35 + 0.32 * rows   # 经验高度；可按需微调
    if (grepl("\\.png$", save_path, ignore.case = TRUE)) {
      png(save_path, width = width_in, height = height_in, units = "in", res = dpi, bg = "transparent")
      grid::grid.draw(legend_grob); dev.off()
    } else if (grepl("\\.pdf$", save_path, ignore.case = TRUE)) {
      pdf(save_path, width = width_in, height = height_in, useDingbats = FALSE)
      grid::grid.draw(legend_grob); dev.off()
    } else {
      warning("save_path 未识别（仅支持 .png/.pdf），未保存：", save_path)
    }
  }
  
  invisible(legend_grob)
}


#' @title Create a Minimalist Violin Plot
#' @description A wrapper around Seurat's VlnPlot to create clean, publication-ready violin plots.
#'
#' @param seurat_object A Seurat object.
#' @param feature The name of the gene or metadata column to plot.
#' @param group.by The metadata column to group cells by.
#' @param palette An optional named vector of colors for the groups.
#' @param pt.size Size of the points to display. Set to 0 to hide points. Default is 0.
#' @param y.lab An optional label for the y-axis. Default is "Expression Level".
#' @param y_breaks_n Number of breaks to aim for on the y-axis. Default is 4.
#' @param ... Additional arguments passed to Seurat::VlnPlot.
#'
#' @return A ggplot object.
#'
#' @title Create a Minimalist Violin Plot (v2)
#' @description A wrapper around Seurat's VlnPlot with enhanced controls for theme, legend, and axes.
#'
#' @param seurat_object A Seurat object.
#' @param feature The name of the gene or metadata column to plot.
#' @param group.by The metadata column to group cells by.
#' @param palette An optional named vector of colors for the groups.
#' @param pt.size Size of the points to display. Set to 0 to hide points. Default is 0.
#' @param show.legend Logical. Whether to display the legend. Default is FALSE.
#' @param x.angle Numeric. The angle of the x-axis text labels. Default is 45.
#' @param y.lab An optional label for the y-axis. Default is "Expression Level".
#' @param y_breaks_n Number of breaks to aim for on the y-axis. Default is 4.
#' @param ... Additional arguments passed to Seurat::VlnPlot.
#'
#' @return A ggplot object.
#'
VlnPlot_minimal <- function(
    seurat_object,
    feature,
    group.by,
    palette = NULL,
    pt.size = 0,
    show.legend = FALSE, # <-- 新增参数
    x.angle = 45,        # <-- 新增参数
    y.lab = "Expression Level",
    y_breaks_n = 4,
    ...
) {
  
  # 1. Create the base violin plot
  p <- Seurat::VlnPlot(
    seurat_object,
    features = feature,
    group.by = group.by,
    pt.size = pt.size,
    ... 
  )
  
  # 2. Apply the minimalist theme customizations
  p <- p + theme(
    axis.title.x = element_blank(),
    # --- 使用 x.angle 参数 ---
    axis.text.x = element_text(angle = x.angle, hjust = 1, size = 12, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    
    # --- 使用 show.legend 参数 ---
    legend.position = if (isTRUE(show.legend)) "right" else "none",
    
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold.italic")
  ) +
    # 3. Beautify the y-axis and set labels
    scale_y_continuous(breaks = scales::pretty_breaks(n = y_breaks_n)) +
    labs(title = feature, y = y.lab)
  
  # 4. Apply custom color palette if provided
  if (!is.null(palette)) {
    p <- p + scale_fill_manual(
      values = palette, 
      name = group.by # Use the group.by column name as the legend title
    )
  }
  
  return(p)
}




CreateNamedPalette <- function(
    object,
    group_by,
    palette_fun = viridisLite::turbo,
    grey_keywords = c("LowQC", "Ambiguous", "Doublet", "Uncertain"),
    grey_color = "grey80"
) {
  
  # --- 1. 检查输入 ---
  if (!group_by %in% colnames(object@meta.data)) {
    stop(paste("Column", group_by, "not found in object's meta.data."))
  }
  
  # --- 2. 获取所有唯一的细胞类别/标签 ---
  all_labels <- levels(object@meta.data[[group_by]])
  if (is.null(all_labels)) {
    all_labels <- sort(unique(as.character(object@meta.data[[group_by]])))
  }
  
  # --- 3. 识别需要设为灰色的类别 ---
  # 使用 grepl 进行不区分大小写的关键词匹配
  is_grey_label <- sapply(all_labels, function(lab) {
    any(sapply(grey_keywords, function(key) {
      grepl(key, lab, ignore.case = TRUE)
    }))
  })
  
  grey_labels <- all_labels[is_grey_label]
  color_labels <- all_labels[!is_grey_label]
  
  # --- 4. 为非灰色类别生成颜色 ---
  num_colors_needed <- length(color_labels)
  if (num_colors_needed > 0) {
    colors <- palette_fun(num_colors_needed)
    color_palette_part <- setNames(colors, color_labels)
  } else {
    color_palette_part <- character(0)
  }
  
  # --- 5. 为灰色类别分配灰色 ---
  if (length(grey_labels) > 0) {
    grey_palette_part <- setNames(rep(grey_color, length(grey_labels)), grey_labels)
  } else {
    grey_palette_part <- character(0)
  }
  
  # --- 6. 合并并按原始顺序排序 ---
  full_palette <- c(color_palette_part, grey_palette_part)
  final_palette <- full_palette[all_labels] # 确保输出顺序与因子水平一致
  
  return(final_palette)
}
# 假设您的 Seurat 对象名为 'object'








DimPlotNice2 <- function(
    object,
    reduction,
    group_by,
    palette = c("viridis","turbo","okabe-ito","scico-roma","glasbey"),
    sort_levels = c("none","alpha","numeric","custom"),
    custom_order = NULL,
    sort_desc   = FALSE,
    pt.size     = 0.15,
    label       = TRUE,
    label_size  = 3.5,
    legend_title = NULL,
    legend_ncol  = 1,
    legend_position = "right",
    legend_pt_size = 3,
    title = NULL,
    # axes / scalebar
    show_axes   = c("none", "simple", "scale"),
    axis_title_x = NULL,
    axis_title_y = NULL,
    axis_text_size = 10,
    axis_title_size = 11,
    scale_bar_length = NULL,
    scale_bar_units = "units",
    scale_bar_position = c("bottom_right", "bottom_left", "top_right", "top_left"),
    scale_bar_height = 0.5,
    scale_bar_text_size = 3,
    # NEW: visual style preset
    style = c("soft","crisp")   # <- two modes merged into this function
){
  # ---------- sanity checks ----------
  if (!reduction %in% names(object@reductions))
    stop("Reduction '", reduction, "' not found in object@reductions.")
  if (!group_by %in% colnames(object@meta.data))
    stop("group_by '", group_by, "' not found in meta.data.")
  show_axes <- match.arg(show_axes)
  scale_bar_position <- match.arg(scale_bar_position)
  sort_levels <- match.arg(sort_levels)
  style <- match.arg(style)
  
  # ---------- level ordering ----------
  vals_chr  <- as.character(object@meta.data[[group_by]])
  uniq_vals <- unique(vals_chr)
  if (sort_levels == "custom") {
    if (is.null(custom_order)) stop("Provide custom_order when sort_levels='custom'.")
    breaks <- as.character(custom_order)
  } else if (sort_levels == "alpha") {
    breaks <- sort(uniq_vals, decreasing = sort_desc)
  } else if (sort_levels == "numeric") {
    num <- suppressWarnings(as.numeric(uniq_vals))
    if (sum(is.finite(num)) >= ceiling(length(uniq_vals)/2)) {
      ord <- order(num, decreasing = sort_desc, na.last = TRUE)
      breaks <- uniq_vals[ord]
    } else {
      breaks <- sort(uniq_vals, decreasing = sort_desc)
    }
  } else {
    breaks <- uniq_vals
  }
  object@meta.data[[group_by]] <- factor(vals_chr, levels = breaks)
  
  # ---------- palette helpers ----------
  are_valid_colors <- function(x) {
    all(vapply(x, function(ci) {
      tryCatch({ grDevices::col2rgb(ci); TRUE }, error=function(e) FALSE)
    }, logical(1)))
  }
  okabe_ito_vec <- c("#E69F00","#56B4E9","#009E73","#F0E442",
                     "#0072B2","#D55E00","#CC79A7","#000000")
  get_palette_fun <- function(spec) {
    if (is.function(spec)) return(spec)
    if (is.character(spec)) {
      if (are_valid_colors(spec)) {
        pal_vec <- spec
        return(function(n) {
          if (length(pal_vec) >= n) pal_vec[seq_len(n)] else {
            if (requireNamespace("Polychrome", quietly = TRUE)) {
              c(pal_vec, Polychrome::glasbey.colors(n - length(pal_vec)))
            } else {
              c(pal_vec, grDevices::hcl.colors(n - length(pal_vec), "Spectral"))
            }
          }
        })
      }
      tokens <- spec
      for (tok in tokens) {
        tok <- tolower(tok)
        if (startsWith(tok, "scico-") && requireNamespace("scico", quietly = TRUE)) {
          pal <- sub("^scico-","", tok)
          return(function(n) scico::scico(n, palette = pal))
        }
        if (requireNamespace("viridisLite", quietly = TRUE)) {
          if (tok == "viridis") return(viridisLite::viridis)
          if (tok == "magma")   return(viridisLite::magma)
          if (tok == "plasma")  return(viridisLite::plasma)
          if (tok == "inferno") return(viridisLite::inferno)
          if (tok == "cividis") return(viridisLite::cividis)
          if (tok == "rocket")  return(viridisLite::rocket)
          if (tok == "mako")    return(viridisLite::mako)
          if (tok == "turbo" && "turbo" %in% ls(getNamespace("viridisLite"))) {
            return(viridisLite::turbo)
          }
        }
        if (tok == "glasbey" && requireNamespace("Polychrome", quietly = TRUE)) {
          return(Polychrome::glasbey.colors)
        }
        if (tok %in% c("okabe-ito","okabe_ito","okabeito","okabe")) {
          return(function(n) rep_len(okabe_ito_vec, n))
        }
      }
    }
    return(function(n) grDevices::hcl.colors(n, "Spectral"))
  }
  
  # ---------- build color vector ----------
  if (is.character(palette) && are_valid_colors(palette) && !is.null(names(palette))) {
    cols_vec <- unname(palette[breaks])
    if (anyNA(cols_vec)) {
      miss <- is.na(cols_vec)
      fill_n <- sum(miss)
      fill_cols <- if (requireNamespace("Polychrome", quietly = TRUE)) {
        Polychrome::glasbey.colors(fill_n)
      } else {
        grDevices::hcl.colors(fill_n, "Spectral")
      }
      cols_vec[miss] <- fill_cols
    }
    names(cols_vec) <- breaks
  } else {
    pal_fun <- get_palette_fun(palette)
    cols_vec <- pal_fun(length(breaks))
    names(cols_vec) <- breaks
  }
  
  # ---------- extract embeddings & df ----------
  emb <- Seurat::Embeddings(object, reduction)
  df  <- data.frame(
    x = emb[,1], y = emb[,2],
    group = as.character(object@meta.data[[group_by]]),
    stringsAsFactors = FALSE
  )
  
  # ---------- ordering per style ----------
  # soft: density ordering; crisp: group ordering
  if (style == "soft") {
    grid_n <- 200L
    gx <- cut(df$x, breaks = grid_n, labels = FALSE)
    gy <- cut(df$y, breaks = grid_n, labels = FALSE)
    idx <- (gx - 1L) * grid_n + gy
    counts <- tabulate(idx, nbins = grid_n*grid_n)
    dens <- counts[idx]
    ord <- order(dens, decreasing = FALSE)  # low density on top
    df <- df[ord, , drop = FALSE]
  } else {
    df <- df[order(df$group), , drop = FALSE]
  }
  
  # ---------- aesthetics per style ----------
  if (style == "soft") {
    alpha_pt   <- 0.6
    halo_on    <- TRUE
    halo_mult  <- 2.0
    halo_col   <- "white"
    halo_alpha <- 1.0
  } else { # crisp
    alpha_pt   <- 1.0
    halo_on    <- TRUE
    halo_mult  <- 1.4
    halo_col   <- "white"
    halo_alpha <- 1.0
  }
  
  # ---------- choose backend (auto scattermore for big N) ----------
  use_scattermore <- (nrow(df) > 50000L) && requireNamespace("scattermore", quietly = TRUE)
  
  # ---------- centroids for labels ----------
  cen <- as.data.frame(stats::aggregate(df[,c("x","y")], list(group=df$group), mean))
  names(cen) <- c("group","x","y")
  
  # ---------- base ggplot canvas ----------
  p <- ggplot2::ggplot() +
    ggplot2::coord_equal() +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.title  = ggplot2::element_blank(),
      axis.text   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      panel.grid  = ggplot2::element_blank(),
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.title= ggplot2::element_text(face = "bold"),
      legend.position = legend_position,
      legend.key.height = grid::unit(0.45, "lines"),
      legend.key.width  = grid::unit(0.8,  "lines")
    )
  
  # ---------- draw points (halo + color) ----------
  if (use_scattermore) {
    # halo layer
    if (halo_on) {
      p <- p + scattermore::geom_scattermore(
        data = df, ggplot2::aes(x = x, y = y),
        pointsize = pt.size * halo_mult,
        alpha = halo_alpha, color = halo_col
      )
    }
    # color layer
    p <- p + scattermore::geom_scattermore(
      data = df, ggplot2::aes(x = x, y = y, color = group),
      pointsize = pt.size, alpha = alpha_pt
    )
  } else {
    if (halo_on) {
      p <- p + ggplot2::geom_point(
        data = df, ggplot2::aes(x = x, y = y),
        shape = 16, size = pt.size * halo_mult,
        alpha = halo_alpha, color = halo_col, stroke = 0
      )
    }
    p <- p + ggplot2::geom_point(
      data = df, ggplot2::aes(x = x, y = y, color = group),
      shape = 16, size = pt.size, alpha = alpha_pt, stroke = 0
    )
  }
  
  # ---------- colors, titles, legend ----------
  p <- p + ggplot2::scale_color_manual(
    values = cols_vec[unique(df$group)],
    drop = FALSE,
    guide = ggplot2::guide_legend(
      ncol = legend_ncol, byrow = TRUE,
      override.aes = list(size = legend_pt_size, alpha = 1)
    )
  )
  if (is.null(legend_title)) legend_title <- group_by
  if (is.null(title))        title <- paste0("UMAP — ", group_by)
  p <- p + ggplot2::labs(color = legend_title, title = title)
  
  # ---------- labels (use ggrepel when available) ----------
  if (isTRUE(label)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = cen, ggplot2::aes(x = x, y = y, label = group),
        size = label_size, fontface = "bold",
        box.padding = 0.25, point.padding = 0.3,
        min.segment.length = 0
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = cen, ggplot2::aes(x = x, y = y, label = group),
        size = label_size, fontface = "bold", vjust = -0.8
      )
    }
  }
  
  # ---------- optional axes / scalebar ----------
  if (show_axes == "simple") {
    dims <- colnames(emb)
    xlab <- if (is.null(axis_title_x)) dims[1] else axis_title_x
    ylab <- if (is.null(axis_title_y)) dims[2] else axis_title_y
    p <- p +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = axis_title_size),
        axis.text  = ggplot2::element_text(size = axis_text_size),
        axis.ticks = ggplot2::element_line()
      )
  } else if (show_axes == "scale") {
    x_range <- range(emb[,1]); y_range <- range(emb[,2])
    if (is.null(scale_bar_length)) {
      x_span <- diff(x_range)
      sbl <- round(x_span / 5, 1)
      if (sbl >= 2) sbl <- round(sbl) else if (sbl >= 1) sbl <- round(sbl*2)/2 else sbl <- round(sbl*10)/10
      scale_bar_length <- sbl
    }
    x_margin <- diff(x_range) * 0.05
    y_margin <- diff(y_range) * 0.05
    bar_pos <- switch(
      scale_bar_position,
      "bottom_right" = list(x = x_range[2] - x_margin - scale_bar_length, y = y_range[1] + y_margin),
      "bottom_left"  = list(x = x_range[1] + x_margin,                         y = y_range[1] + y_margin),
      "top_right"    = list(x = x_range[2] - x_margin - scale_bar_length, y = y_range[2] - y_margin),
      "top_left"     = list(x = x_range[1] + x_margin,                         y = y_range[2] - y_margin)
    )
    scale_data <- data.frame(
      x = bar_pos$x, xend = bar_pos$x + scale_bar_length,
      y = bar_pos$y, yend = bar_pos$y
    )
    text_data <- data.frame(
      x = bar_pos$x + scale_bar_length/2,
      y = bar_pos$y + diff(y_range) * 0.02,
      label = paste0(scale_bar_length, " ", scale_bar_units)
    )
    p <- p +
      ggplot2::geom_segment(
        data = scale_data,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "black", linewidth = scale_bar_height, inherit.aes = FALSE
      ) +
      ggplot2::geom_segment(
        data = scale_data,
        ggplot2::aes(x = x, y = y - diff(y_range) * 0.005, xend = x, yend = y + diff(y_range) * 0.005),
        color = "black", linewidth = scale_bar_height, inherit.aes = FALSE
      ) +
      ggplot2::geom_segment(
        data = scale_data,
        ggplot2::aes(x = xend, y = y - diff(y_range) * 0.005, xend = xend, yend = y + diff(y_range) * 0.005),
        color = "black", linewidth = scale_bar_height, inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = text_data,
        ggplot2::aes(x = x, y = y, label = label),
        size = scale_bar_text_size, inherit.aes = FALSE, fontface = "bold"
      )
  }
  
  return(p)
}

# Nicely formatted LR axis × direction bubble plot
PlotLRaxisBubble <- function(
    df,
    # column names (keep defaults if you follow current pipeline)
    col_direction = "direction",
    col_axis      = "axis",
    col_n_pairs   = "n_pairs",
    col_delta     = "sum_delta",
    
    # labels
    x_label       = "Direction (Sender\n→ Receiver)",
    y_label       = "Pathway axis (myCAF ↔ SPP1+ TAM)",
    color_label   = "Σ Δscore\n(Hotspot - Control)",
    size_label    = "No. of LR pairs",
    
    # appearance
    base_size     = 14,
    x_angle       = 45,
    y_wrap_width  = 30,
    point_alpha   = 0.9,
    size_range    = c(3, 9),
    color_low     = "grey85",
    color_high    = "firebrick",
    
    # save options
    save_path     = NULL,
    width         = 13,
    height        = 8,
    dpi           = 300
){
  # ---- basic checks ----
  stopifnot(is.data.frame(df))
  needed <- c(col_direction, col_axis, col_n_pairs, col_delta)
  if (!all(needed %in% colnames(df))) {
    stop("Data frame must contain columns: ", paste(needed, collapse = ", "))
  }
  
  # ---- data prep ----
  plot_data <- df %>%
    dplyr::mutate(
      # normalize arrows to ASCII then reformat for plotting
      direction_plot = stringr::str_replace_all(.data[[col_direction]], "→", "->"),
      axis_plot      = stringr::str_replace_all(.data[[col_axis]],      "→", "->"),
      # reorder axis by sum of delta (larger on top)
      axis_plot      = forcats::fct_reorder(axis_plot, .data[[col_delta]])
    )
  
  # ---- build plot ----
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$direction_plot,
      y = .data$axis_plot
    )
  ) +
    # bubbles
    ggplot2::geom_point(
      ggplot2::aes(
        size  = .data[[col_n_pairs]],
        color = .data[[col_delta]]
      ),
      alpha = point_alpha
    ) +
    
    # color scale (NES-like)
    ggplot2::scale_color_gradient(
      low  = color_low,
      high = color_high,
      name = color_label
    ) +
    
    # size scale
    ggplot2::scale_size_continuous(
      range = size_range,
      name  = size_label
    ) +
    
    # X labels: put newline before arrow to stack sender/receiver
    ggplot2::scale_x_discrete(
      labels = function(x) {
        # "Basal-like -> myCAF" -> "Basal-like\n→ myCAF"
        x %>%
          stringr::str_replace_all("->", "→") %>%
          stringr::str_replace(" → ", "\n→ ")
      }
    ) +
    
    # Y labels: wrap long axis names
    ggplot2::scale_y_discrete(
      labels = function(x) stringr::str_wrap(x, width = y_wrap_width)
    ) +
    
    ggplot2::labs(
      x = x_label,
      y = y_label
    ) +
    
    # theme
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      # grid: keep horizontal lines, drop vertical
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90"),
      panel.grid.minor   = ggplot2::element_blank(),
      # panel border
      panel.border       = ggplot2::element_rect(
        color = "black",
        fill  = NA,
        linewidth = 0.6
      ),
      # axis text
      axis.text.x = ggplot2::element_text(
        color = "black",
        size  = base_size * 0.8,
        angle = x_angle,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = ggplot2::element_text(
        color = "black",
        size  = base_size * 0.8
      ),
      axis.title  = ggplot2::element_text(face = "bold"),
      # legend
      legend.position   = "right",
      legend.title      = ggplot2::element_text(face = "bold"),
      legend.key.height = grid::unit(0.6, "cm"),
      legend.key.width  = grid::unit(0.4, "cm"),
      # margins: extra space bottom for tilted x labels, left for long y labels
      plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 18, l = 20, unit = "mm")
    )
  
  # ---- optional save ----
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot     = p,
      width    = width,
      height   = height,
      dpi      = dpi,
      useDingbats = FALSE
    )
  }
  
  return(p)
}


library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)

# --- 定义通用气泡图函数 ---
PlotLRaxisBubble_Fixed <- function(plot_data, save_path = NULL, width = 9, height = 7) {
  
  # 1. 安全检查
  if (nrow(plot_data) == 0) {
    warning("No data to plot! Check your filtering thresholds.")
    return(NULL)
  }
  
  # 2. 智能排序 Y 轴 (Axis)
  # 按照这些互作在核心细胞间的总强度排序，让强的排上面
  axis_order <- plot_data %>%
    group_by(axis) %>%
    summarise(total = sum(sum_delta, na.rm = TRUE)) %>%
    arrange(total) %>% 
    pull(axis)
  
  plot_data$axis <- factor(plot_data$axis, levels = axis_order)
  
  # 3. 智能排序 X 轴 (Direction)
  # 让发送者聚在一起：先排 myCAF发出的，再排 TAM发出的...
  # 这里做一个简单的字母排序，或者你可以自定义 factor levels
  plot_data$direction <- factor(plot_data$direction, levels = sort(unique(plot_data$direction)))
  
  # 4. 绘图
  p <- ggplot(plot_data, aes(x = direction, y = axis)) +
    # 画气泡
    geom_point(aes(size = sum_delta, color = sum_delta), alpha = 0.9) +
    
    # 颜色：使用 Magma 配色 (深色代表强)，反转让深色对应高值
    scale_color_viridis_c(option = "magma", direction = -1, 
                          name = "Enrichment Score\n(Delta)") +
    
    # 大小：控制气泡范围
    scale_size_continuous(range = c(3, 8), name = "Interaction\nStrength") +
    
    # 主题美化
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black", face = "italic"), # 基因名斜体
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.position = "right"
    ) +
    
    # 标签
    labs(
      title = "Core Signaling Axes in ECM Hotspot",
      subtitle = "Restricted to: myCAF, SPP1+ TAM, Basal-like",
      x = NULL, 
      y = NULL
    )
  
  # 5. 保存
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = width, height = height)
    message(paste0("Plot saved to: ", save_path))
  }
  
  return(p)
}


library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)

# 美化版 + 字符安全版
PlotLRaxisBubble_Fixed <- function(
    plot_data,
    save_path = NULL,
    width  = 10,
    height = 8,
    title    = "Core signaling axes in ECM Hotspot",
    subtitle = "Sender -> Receiver (Hotspot-enriched LR axes)" # 改为普通箭头
){
  # -------- 0. basic check --------
  if (nrow(plot_data) == 0) {
    warning("No data to plot! Check your filtering thresholds.")
    return(NULL)
  }
  
  # make sure we work with characters first
  plot_data <- plot_data %>%
    mutate(
      axis      = as.character(axis),
      direction = as.character(direction)
    )
  
  # -------- 1. order Y-axis (axis) by total strength --------
  axis_order <- plot_data %>%
    group_by(axis) %>%
    summarise(total = sum(sum_delta, na.rm = TRUE), .groups = "drop") %>%
    arrange(total) %>% # ggplot Y轴是从下到上，所以这里从小到大排
    pull(axis)
  
  plot_data$axis <- factor(plot_data$axis, levels = axis_order)
  
  # -------- 2. parse source & target --------
  # 假设方向是用 "->" 或 "->" 连接的，这里统一处理
  # 使用简单的 str_split 更加稳健
  plot_data$source <- sapply(str_split(plot_data$direction, "->"), `[`, 1)
  plot_data$target <- sapply(str_split(plot_data$direction, "->"), `[`, 2)
  
  # 去除空格
  plot_data$source <- str_trim(plot_data$source)
  plot_data$target <- str_trim(plot_data$target)
  
  # preferred biological sender order
  preferred_sources <- c("Basal-like", "myCAF", "SPP1+ TAM") # 核心三角顺序
  source_order <- intersect(preferred_sources, unique(plot_data$source))
  source_order <- c(source_order, setdiff(sort(unique(plot_data$source)), source_order))
  
  # order directions
  plot_data <- plot_data %>%
    arrange(factor(source, levels = source_order), target)
  dir_levels <- unique(plot_data$direction)
  plot_data$direction <- factor(plot_data$direction, levels = dir_levels)
  
  # -------- 3. build plot --------
  p <- ggplot(plot_data, aes(x = direction, y = axis)) +
    # bubble
    geom_point(aes(size = sum_delta, colour = sum_delta), alpha = 0.9) +
    
    # color: viridis magma
    scale_colour_viridis_c(
      option    = "magma",
      direction = -1,
      name      = "Delta Score" # 改为英文
    ) +
    
    # size
    scale_size_continuous(
      range = c(3, 9),
      name  = "Total Delta" # 改为英文
    ) +
    
    # x-axis labels: 换行显示箭头
    scale_x_discrete(labels = function(x) gsub("->", "\n->\n", x)) +
    
    # theme
    theme_bw(base_size = 14) +
    theme(
      axis.title      = element_blank(),
      axis.text.x     = element_text(
        angle = 45, hjust = 1, vjust = 1,
        colour = "black", face = "bold",
        lineheight = 0.8
      ),
      axis.text.y     = element_text(
        colour = "black",
        face   = "italic" 
      ),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title    = title,
      subtitle = subtitle
    )
  
  # -------- 4. save --------
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = width, height = height)
    message("Plot saved to: ", save_path)
  }
  
  return(p)
}



RunLR_Merged_Hotspot_vs_Control<- function(
    object,
    group_col        = "ECMhot_core",
    group_hotspot    = TRUE,
    group_control    = FALSE,
    label_col        = "final_label",
    focus_cells      = c("Basal-like","myCAF","SPP1+ TAM"),
    min_cells        = 20,
    delta_thresh     = 0.1,
    out_dir          = NULL,
    liana_resource   = "OmniPath",
    liana_methods    = c("natmi","connectome","logfc","sca","cellphonedb"),
    # ★ 新增参数 ★
    renormalize      = TRUE,
    norm_assay       = "RNA",
    norm_method      = "LogNormalize",
    norm_scale_factor = 1e4
){
  # --- 1. 输出目录 ---
  if (is.null(out_dir)) {
    base_dir <- if (exists("Config") && !is.null(Config$ECMFront)) Config$ECMFront else "."
    out_dir  <- file.path(base_dir, "LIANA_Merged_Analysis")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  message(">>> Starting LIANA analysis on merged object (V3, with optional renormalization)...")
  
  # --- 2. 处理 Seurat v5 layers ---
  if (any(grepl("\\.", Layers(object)))) {
    message("   Detected split layers (Seurat v5). Joining layers...")
    object <- JoinLayers(object)
  }
  
  # --- 3. 统一归一化（关键修复） ---
  if (renormalize) {
    message("   Renormalizing merged object before subsetting...")
    DefaultAssay(object) <- norm_assay
    object <- NormalizeData(
      object,
      normalization.method = norm_method,
      scale.factor         = norm_scale_factor,
      verbose              = FALSE
    )
  } else {
    message("   [Warning] renormalize = FALSE, will use existing data as-is.")
  }
  
  message(">>> Comparison: ", group_col, " ", group_hotspot, " (Hotspot) vs ", group_control, " (Control)")
  
  # --- 4. 取 Hotspot / Control 细胞 ---
  md <- object@meta.data
  if (!group_col %in% colnames(md)) {
    stop("group_col not found in meta.data: ", group_col)
  }
  
  cells_hotspot <- colnames(object)[md[[group_col]] == group_hotspot]
  cells_control <- colnames(object)[md[[group_col]] == group_control]
  cells_hotspot <- stats::na.omit(cells_hotspot)
  cells_control <- stats::na.omit(cells_control)
  
  message("   Hotspot cells: ", length(cells_hotspot))
  message("   Control  cells: ", length(cells_control))
  if (length(cells_hotspot) == 0 || length(cells_control) == 0) {
    stop("No cells found for one of the groups based on ", group_col)
  }
  
  # --- 5. 辅助函数：过滤稀有细胞类型 ---
  filter_rare_celltypes <- function(seu, min_n) {
    if (is.null(seu) || ncol(seu) == 0) return(NULL)
    ct_counts <- table(Idents(seu))
    keep      <- names(ct_counts[ct_counts >= min_n])
    if (length(keep) == 0) return(NULL)
    subset(seu, idents = keep)
  }
  
  # --- 6. 构建子对象并设置 Idents ---
  obj_hotspot <- subset(object, cells = cells_hotspot)
  obj_control <- subset(object, cells = cells_control)
  
  Idents(obj_hotspot) <- label_col
  Idents(obj_control) <- label_col
  
  obj_hotspot <- filter_rare_celltypes(obj_hotspot, min_cells)
  obj_control <- filter_rare_celltypes(obj_control, min_cells)
  
  if (is.null(obj_hotspot) || is.null(obj_control)) {
    stop("After filtering for min_cells, one group became empty.")
  }
  
  # --- 7. 跑 LIANA ---
  run_liana_one <- function(seu) {
    message("   Running LIANA on subset: ", ncol(seu), " cells...")
    DefaultAssay(seu) <- norm_assay    # 确保使用刚刚归一化的 assay
    res <- liana_wrap(seu, resource = liana_resource, method = liana_methods)
    agg <- liana_aggregate(res)
    agg$score <- -log10(pmax(agg$aggregate_rank, 1e-6))
    return(agg)
  }
  
  message(">>> Processing Hotspot group...")
  liana_hotspot <- run_liana_one(obj_hotspot)
  liana_hotspot$region <- "Hotspot"
  
  message(">>> Processing Control group...")
  liana_control <- run_liana_one(obj_control)
  liana_control$region <- "Control"
  
  # --- 8. 计算 Delta ---
  message(">>> Calculating delta scores...")
  liana_combined <- dplyr::bind_rows(liana_hotspot, liana_control)
  
  pair_deltas_all <- liana_combined %>%
    dplyr::select(source, target, ligand.complex, receptor.complex, region, score) %>%
    tidyr::pivot_wider(names_from = region, values_from = score, values_fill = 0) %>%
    dplyr::mutate(
      delta = Hotspot - Control,
      up_in = dplyr::case_when(
        delta >  delta_thresh ~ "Enriched in Hotspot",
        delta < -delta_thresh ~ "Enriched in Control",
        TRUE                  ~ "No Change"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(delta))
  
  readr::write_csv(pair_deltas_all, file.path(out_dir, "LIANA_Delta_AllPairs.csv"))
  
  # --- 9. 只看 focus_cells & Hotspot↑ ---
  pair_deltas_focus <- pair_deltas_all %>%
    dplyr::filter(
      source %in% focus_cells,
      target %in% focus_cells,
      up_in  == "Enriched in Hotspot"
    )
  
  if (nrow(pair_deltas_focus) == 0) {
    warning("No enriched pairs found among focus_cells.")
    return(NULL)
  }
  
  pair_deltas_focus <- pair_deltas_focus %>%
    dplyr::mutate(axis = paste0(ligand.complex, "-", receptor.complex))
  
  axis_summary <- pair_deltas_focus %>%
    dplyr::mutate(delta_pos = pmax(delta, 0)) %>%
    dplyr::group_by(axis) %>%
    dplyr::summarise(
      n_pairs   = dplyr::n(),
      sum_delta = sum(delta_pos, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(sum_delta))
  
  flow_axis <- pair_deltas_focus %>%
    dplyr::mutate(
      flow      = paste(source, "->", target),
      delta_pos = pmax(delta, 0)
    ) %>%
    dplyr::group_by(flow, axis) %>%
    dplyr::summarise(
      sum_delta = sum(delta_pos, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(sum_delta))
  
  readr::write_csv(pair_deltas_focus, file.path(out_dir, "LIANA_Delta_FocusCells.csv"))
  readr::write_csv(axis_summary,      file.path(out_dir, "LIANA_AxisSummary.csv"))
  
  message(">>> LIANA merged analysis completed. Results at: ", out_dir)
  
  return(list(
    pair_deltas_all   = pair_deltas_all,
    pair_deltas_focus = pair_deltas_focus,
    axis_summary      = axis_summary,
    flow_axis         = flow_axis
  ))
}




library(Seurat)
library(dplyr)
library(msigdbr)
library(clusterProfiler)

#' Run DEG Enrichment (V2 - Two Directions)
#' 
#' 计算每个样本 Hotspot vs Background 的差异基因，并调用 enrichment_dotplot_generic 绘图。
#' 
run_deg_enrichment_v2 <- function(
    object,
    output_dir,
    group_col = "ECM_Hotspot",   # 分组列
    split_col = "orig.ident",    # 样本列
    species   = "Homo sapiens",
    category  = "H",             # MSigDB Category
    subcategory = NULL,
    logfc_thr = 0.25,
    pval_thr  = 0.05,
    
    # 绘图参数
    top_terms = 8,               # 每个方向展示多少个通路
    fig_width = 14, 
    fig_height = 10
) {
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # --- 1. 计算差异基因 (DE Analysis) ---
  samples <- unique(object[[split_col]][,1])
  de_list <- list()
  
  message(">>> [Step 1] Calculating DEGs for each sample (Hotspot vs Background)...")
  
  for (sam in samples) {
    message(paste0("   Processing: ", sam))
    
    # 提取子集
    sub_obj <- subset(object, cells = colnames(object)[object[[split_col]] == sam])
    Idents(sub_obj) <- group_col
    
    # 检查是否有两组
    if (length(unique(Idents(sub_obj))) < 2) {
      warning(paste("Sample", sam, "missing groups. Skipping."))
      next
    }
    
    # 运行 FindMarkers (双向)
    # ident.1 = TRUE (Hotspot), ident.2 = FALSE (Background)
    # LogFC > 0 是 Hotspot, LogFC < 0 是 Background
    tryCatch({
      degs <- FindMarkers(sub_obj, ident.1 = "TRUE", ident.2 = "FALSE", 
                          only.pos = FALSE,  # ★ 关键：同时找上调和下调
                          logfc.threshold = logfc_thr,
                          min.pct = 0.1,
                          verbose = FALSE)
      
      # 整理格式以适配通用绘图函数
      degs$gene <- rownames(degs)
      degs$sample <- sam
      
      # 定义方向 (Direction)
      degs$direction_label <- ifelse(degs$avg_log2FC > 0, 
                                     "Hotspot Core",   # TRUE
                                     "IF Background")  # FALSE
      
      de_list[[sam]] <- degs
      
    }, error = function(e) {
      message(paste("Error in sample", sam, ":", e$message))
    })
  }
  
  # 合并所有结果
  de_all <- bind_rows(de_list)
  
  if (nrow(de_all) == 0) stop("No DEGs found.")
  
  # 保存原始 DE 表
  write.csv(de_all, file.path(output_dir, "DEGs_Hotspot_vs_Background_All.csv"), row.names = FALSE)
  
  # --- 2. 准备绘图参数 ---
  message(">>> [Step 2] Running Enrichment & Plotting...")
  
  # 定义颜色 (Hotspot 用红色系，Background 用蓝色系)
  my_dir_cols <- c(
    "Hotspot Core"  = "#E41A1C", # Red
    "IF Background" = "#377EB8"  # Blue
  )
  
  # 定义 MSigDB 源
  msig_info <- list(species = species, collection = category, subcategory = subcategory)
  
  # --- 3. 调用通用绘图函数 ---
  # 这里直接利用 enrichment_dotplot_generic 的强大功能
  res <- enrichment_dotplot_generic(
    de_all = de_all,
    direction_col = "direction_label", # 我们刚才定义的列名
    dir_levels = c("Hotspot Core", "IF Background"), # 指定显示顺序
    dir_cols = my_dir_cols,            # 指定颜色
    
    msig_source = msig_info,           # 通路数据库
    
    q_cut = pval_thr,                  # P值阈值
    fc_cut = logfc_thr,                # FC 阈值
    top_terms_per_group = top_terms,   # 展示数量
    
    facet_by_sample = TRUE,            # 按样本分面
    shade_direction = TRUE,            # 开启背景色块
    direction_to_legend = TRUE,        # 方向放在图例里，X轴干净
    
    save_prefix = file.path(output_dir, "Enrichment_Comparison"),
    fig_width = fig_width,
    fig_height = fig_height
  )
  
  return(res)
}



