# ==== README / 目录 ====
# VISIUMHDcode 函数分区（按出现顺序）：
# 1. 基础坐标与QC工具
# 2. 数据导入与生成
# 3. 调色与空间分布
# 4. 表达可视化与QC报告
# 5. 标签导出与分布变体
# 6. Seurat流水线与平滑
# 7. RCTD/Cellbin 处理
# 8. 一致性评估与锚定注释
# 9. BANKSY/可视化工具集
# 10. BANKSY/三重生态/边界
# 11. 子集/标签吸收/SoftGate
# 12. Hotspot 与 LR 分析

# ==== 基础坐标与QC工具 ====




suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(ggplot2)
})
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(Matrix)
})
.pick_umap <- function(obj){
  rn <- Reductions(obj)
  if ("seurat_umap" %in% rn) return("seurat_umap")
  if ("umap" %in% rn) return("umap")
  return("umap")
}









suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Please install 'Seurat'")
  if (!requireNamespace("arrow", quietly = TRUE))  stop("Please install 'arrow'")
  if (!requireNamespace("rjson", quietly = TRUE))  stop("Please install 'rjson'")
  if (!requireNamespace("dplyr", quietly = TRUE))  stop("Please install 'dplyr'")
  if (!requireNamespace("png", quietly = TRUE))    stop("Please install 'png'")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Please install 'tibble'")
  if (!requireNamespace("grid", quietly = TRUE))   stop("Please install 'grid'")
  library(Seurat); library(dplyr); library(tibble)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

.parse_bin_um <- function(size_str) {
  x <- suppressWarnings(as.numeric(sub("^0*([0-9]+)um$", "\\1", tolower(size_str))))
  if (is.na(x) || x <= 0) stop("Invalid bin 'size': ", size_str)
  x
}





# 1) 确保有可用的 log-normalized data（仅对非 SCT 进行 NormalizeData 以填充 data layer）
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
get_xy_um <- function(
    object,
    prefer = c("auto","scaled","raw","um"),   # “auto”=自动优先 scaled
    verbose = TRUE
){
  stopifnot(inherits(object, "Seurat"))
  prefer <- match.arg(prefer)
  md <- object@meta.data
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  
  # ---------- 0) 已有“微米坐标”优先使用 ----------
  # 常见命名：x_um/y_um 或 imagecol_um/imagerow_um
  um_x_candidates <- intersect(c("x_um","imagecol_um","x_um_sp","x_um_hd"), colnames(md))
  um_y_candidates <- intersect(c("y_um","imagerow_um","y_um_sp","y_um_hd"), colnames(md))
  has_um_cols <- length(um_x_candidates) == 1 && length(um_y_candidates) == 1
  
  # ---------- 1) 选择像素坐标列 ----------
  # scaled 优先（与 RedrawPeripheryWithPSD 的默认一致）
  x_scaled <- "imagecol_scaled"; y_scaled <- "imagerow_scaled"
  x_raw    <- "imagecol";        y_raw    <- "imagerow"
  
  has_scaled <- all(c(x_scaled, y_scaled) %in% colnames(md))
  has_raw    <- all(c(x_raw,    y_raw)    %in% colnames(md))
  
  src <- list()
  if (prefer == "um" && has_um_cols) {
    src$mode <- "um"
    xcol <- um_x_candidates[1]; ycol <- um_y_candidates[1]
    XY_um <- as.matrix(md[, c(xcol, ycol)])
    colnames(XY_um) <- c("x_um","y_um")
    info <- list(
      already_um = TRUE,
      used_cols = c(xcol, ycol),
      px2um = 1,
      microns_per_pixel = NA,
      tissue_lowres_scalef = NA,
      source = "existing *_um columns"
    )
    if (verbose) message("[get_xy_um] Using existing micrometer columns: ", paste(info$used_cols, collapse=", "))
    return(list(XY_um = XY_um, info = info))
  }
  
  if (prefer %in% c("auto","scaled")) {
    if (has_scaled) {
      xcol <- x_scaled; ycol <- y_scaled; src$mode <- "scaled"
    } else if (has_raw) {
      xcol <- x_raw;    ycol <- y_raw;    src$mode <- "raw"
    } else if (has_um_cols) {
      # 没有像素列但有 *_um，直接返回
      src$mode <- "um"
      xcol <- um_x_candidates[1]; ycol <- um_y_candidates[1]
      XY_um <- as.matrix(md[, c(xcol, ycol)])
      colnames(XY_um) <- c("x_um","y_um")
      info <- list(
        already_um = TRUE,
        used_cols = c(xcol, ycol),
        px2um = 1,
        microns_per_pixel = NA,
        tissue_lowres_scalef = NA,
        source = "existing *_um columns"
      )
      if (verbose) message("[get_xy_um] Using existing micrometer columns: ", paste(info$used_cols, collapse=", "))
      return(list(XY_um = XY_um, info = info))
    } else {
      stop("[get_xy_um] No coordinate columns found (expected imagecol(_scaled)/imagerow(_scaled) or *_um).")
    }
  } else if (prefer == "raw") {
    if (!has_raw) stop("[get_xy_um] Requested raw coords but imagecol/imagerow not found.")
    xcol <- x_raw; ycol <- y_raw; src$mode <- "raw"
  }
  
  # ---------- 2) 读取尺度：microns_per_pixel（mpp）与 tissue_lowres_scalef（s_low） ----------
  mpp <- object@misc$microns_per_pixel %||%
    object@misc$spatial_scales$microns_per_pixel %||%
    object@misc$scales$microns_per_pixel %||%
    tryCatch(Images(object)[[1]]@scale.factors$microns_per_pixel, error = function(e) NULL)
  
  s_low <- object@misc$spatial_scales$tissue_lowres_scalef %||%
    object@misc$scales$tissue_lowres_scalef %||%
    tryCatch(Images(object)[[1]]@scale.factors$tissue_lowres_scalef, error = function(e) NULL)
  
  # 若 s_low 缺失但同时有 scaled 与 hires 列，可用两列比值估算
  if (is.null(s_low) && all(c("imagecol_scaled","imagecol_hires") %in% colnames(md))) {
    ratio <- md$imagecol_scaled / md$imagecol_hires
    s_low <- stats::median(ratio[is.finite(ratio)], na.rm = TRUE)
    if (verbose) message(sprintf("[get_xy_um] tissue_lowres_scalef estimated from scaled/hires = %.6f", s_low))
  }
  
  if (is.null(mpp) || is.null(s_low)) {
    # 若使用的是 *_um 列则不需要 mpp/s_low；其余情况必须有
    stop("[get_xy_um] Missing scale factors: microns_per_pixel (", mpp, ") or tissue_lowres_scalef (", s_low, ").")
  }
  
  # ---------- 3) 计算像素→微米换算系数 ----------
  # 注意：这里使用与 RedrawPeripheryWithPSD 相同的定义：
  #   低清像素（imagecol[_scaled]）× (mpp / s_low) = 微米
  px2um <- as.numeric(mpp) / as.numeric(s_low)
  if (!is.finite(px2um) || px2um <= 0) stop("[get_xy_um] Invalid px2um computed.")
  
  # ---------- 4) 换算 ----------
  XY_px <- as.matrix(md[, c(xcol, ycol)])
  XY_um <- XY_px * px2um
  colnames(XY_um) <- c("x_um","y_um")
  
  info <- list(
    already_um = FALSE,
    used_cols = c(xcol, ycol),
    px2um = px2um,
    microns_per_pixel = mpp,
    tissue_lowres_scalef = s_low,
    source = paste0(src$mode, " pixels → microns")
  )
  
  if (verbose) {
    message(sprintf("[get_xy_um] Using %s columns: %s, %s", src$mode, xcol, ycol))
    message(sprintf("[get_xy_um] microns_per_pixel = %.6f; tissue_lowres_scalef = %.6f; px2um = %.6f µm/px",
                    mpp, s_low, px2um))
  }
  
  return(list(XY_um = XY_um, info = info))
}




# Create or repair the 'data' layer WITHOUT relying on LayerNames(), Assays(), etc.
ensure_log_data <- function(obj, assay = "Spatial", scale.factor = 1e4) {
  # (A) verify assay exists
  assay_ok <- tryCatch(!is.null(obj[[assay]]), error = function(e) FALSE)
  if (!assay_ok) stop(sprintf("Assay '%s' not found on object.", assay))
  
  # (B) probe 'data' by attempting to read it
  has_data <- tryCatch({
    tmp <- GetAssayData(obj, assay = assay, layer = "data")
    !(is.null(tmp) ||
        (inherits(tmp, "dgCMatrix") && Matrix::nnzero(tmp) == 0) ||
        (is.matrix(tmp) && all(tmp == 0)))
  }, error = function(e) FALSE)
  
  if (!has_data) {
    # (C) build log-normalized 'data' from 'counts' in a sparse-safe way
    C <- GetAssayData(obj, assay = assay, layer = "counts")  # dgCMatrix expected
    if (!inherits(C, "dgCMatrix")) C <- as(C, "dgCMatrix")
    
    lib <- Matrix::colSums(C)
    lib[!is.finite(lib) | lib <= 0] <- 1
    
    # column-wise normalize sparsely: t(t(C)/lib) keeps sparsity in Matrix
    N <- t(t(C) / lib) * scale.factor       # still dgCMatrix
    N@x <- log1p(N@x)                       # apply log1p only to non-zeros
    
    obj <- SetAssayData(obj, assay = assay, layer = "data", new.data = N)
  }
  
  return(obj)
}






# 3) 空间 QC 可视化（任意 meta 列都能画）
plot_spatial_meta <- function(obj, col, title = NULL, ptsize = 0.6, assay = "Spatial"){
  obj <- ensure_log_data(obj, assay = assay)
  df <- obj@meta.data %>%
    mutate(x = dplyr::coalesce(imagecol_scaled, imagecol),
           y = -dplyr::coalesce(imagerow_scaled, imagerow))
  pal <- c("#FFFFD4","#FEE391","#FEC44F","#FE9929","#D95F0E")  # YlOrBr
  ggplot(df, aes(x = x, y = y, color = .data[[col]])) +
    geom_point(size = ptsize) +
    coord_fixed() + theme_void() +
    scale_color_gradientn(colors = pal) +
    labs(title = title %||% col, color = col)
}





# 
plot_spatial_meta <- function(obj, value_col, title=""){
  md <- obj@meta.data
  if (!value_col %in% colnames(md)) stop("meta 无列: ", value_col)
  md$imagecol_scaled <- if ("imagecol_scaled" %in% names(md)) md$imagecol_scaled else md$imagecol
  md$imagerow_scaled <- if ("imagerow_scaled" %in% names(md)) md$imagerow_scaled else md$imagerow
  
  if (exists("PlotExpressionV2")) {
    df <- md[, c("barcode","imagecol_scaled","imagerow_scaled", value_col)]
    colnames(df) <- c("barcode","imagecol_scaled","imagerow_scaled","val")
    p <- PlotExpressionV2(
      barcodes = df, Gene = "val",
      ptsize = 1.6, shape = "circle",
      colors = c("#FFFFE5","#FFF7BC","#FEC44F","#D95F0E"),  # 黄系更显色
      legend_title = title
    ) + ggtitle(title)
  } else {
    # fallback：ggplot + scattermore(如有)
    has_scattermore <- requireNamespace("scattermore", quietly = TRUE)
    gg <- ggplot(md, aes(x=imagecol_scaled, y=-imagerow_scaled, color=.data[[value_col]])) +
      (if (has_scattermore) scattermore::geom_scattermore(pointsize=2) else geom_point(size=.2, alpha=.8)) +
      scale_color_viridis_c(option="B") + coord_fixed() + theme_void() + ggtitle(title) +
      theme(plot.title = element_text(hjust=.5, face="bold"))
    p <- gg
  }
  p
}






quick_summary <- function(obj, label){
  cat("\n=== ", label, " ===\n", sep = "")
  x <- obj@meta.data
  fmt <- function(v) sprintf("median=%.0f [IQR %.0f–%.0f]", median(v,na.rm=TRUE), quantile(v,.25,na.rm=TRUE), quantile(v,.75,na.rm=TRUE))
  cat("Spots:", ncol(obj), " | Genes:", nrow(obj), "\n")
  cat("nUMI   : ", fmt(x$nUMI),  "\n")
  cat("nGene  : ", fmt(x$nGene), "\n")
  cat("Detect%:  median=", sprintf("%.1f%%", 100*median(x$detected_rate, na.rm=TRUE)), "\n")
}





# PlotViolinQC

suppressPackageStartupMessages({
  library(Seurat); library(dplyr)
  # patchwork 只是排版用；没装也不影响单图输出
  has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
})

# ——— 统一的小工具：坐标兜底（优先 scaled）
.get_xy_df <- function(md) {
  stopifnot(all(c("imagecol","imagerow") %in% colnames(md)))
  dplyr::mutate(md,
                imagecol_scaled = dplyr::coalesce(.data$imagecol_scaled, .data$imagecol),
                imagerow_scaled = dplyr::coalesce(.data$imagerow_scaled, .data$imagerow)
  )
}





# ——— 1) 小提琴：完全走 Seurat::VlnPlot（meta 列，不涉及 layer）
PlotViolinQC <- function(object, qc_feature, group_by = "orig.ident") {
  if (!qc_feature %in% colnames(object@meta.data)) {
    stop("QC指标 '", qc_feature, "' 在 meta.data 中未找到。")
  }
  VlnPlot(
    object, features = qc_feature, group.by = group_by, pt.size = 0
  ) + NoLegend() + labs(title = qc_feature, x = "")
}





# ——— 2) 空间 QC：改为走你 pipeline 的 PlotExpressionV2（和基因表达同风格）
#      颜色默认黄→橙（和你之前 Score_ECM/ITGAV 的调色一致）
PlotSpatialQC <- function(object, qc_feature, ptsize = 1.6,
                          palette = c("#FFFFD4","#FEE391","#FEC44F","#FE9929","#D95F0E"),
                          legend_title = NULL) {
  
  if (!qc_feature %in% colnames(object@meta.data)) {
    stop("在 meta.data 找不到 QC 指标 '", qc_feature, "'")
  }
  df <- .get_xy_df(object@meta.data) |>
    dplyr::select(barcode, imagecol_scaled, imagerow_scaled, !!qc_feature)
  
  # 用你现有的 PlotExpressionV2 出图（把 QC 指标当作“Gene”列即可）
  p <- PlotExpressionV2(
    barcodes     = df,
    Gene         = qc_feature,
    ptsize       = ptsize,
    shape        = "circle",
    colors       = palette,
    legend_title = legend_title %||% qc_feature
  )
  p
}





# ——— 3) 一键 QC 总览：优先使用你 pipeline；自动补 percent.mt（人类基因名默认 ^MT-）
PlotQC_Summary <- function(object, group_by = "orig.ident",
                           qc_features = c("nCount_Spatial","nFeature_Spatial","percent.mt")) {
  
  # 自动补 percent.mt（如果没有）
  if ("percent.mt" %in% qc_features && !"percent.mt" %in% colnames(object@meta.data)) {
    message("• 计算 percent.mt …")
    # FFPE/人类默认 ^MT-；如是小鼠请自己改成 ^mt-
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  }
  
  # 过滤掉确实不存在的列（以免报错）
  qc_features <- qc_features[qc_features %in% colnames(object@meta.data)]
  if (!length(qc_features)) stop("没有可用的 QC 指标。")
  
  # 小提琴（Seurat 原生）
  vln_list <- lapply(qc_features, function(f) PlotViolinQC(object, f, group_by = group_by))
  
  # 空间图（走 PlotExpressionV2）
  sp_list  <- lapply(qc_features, function(f) PlotSpatialQC(object, f))
  
  if (has_patchwork) {
    final_plot <- patchwork::wrap_plots(vln_list, ncol = length(qc_features)) /
      patchwork::wrap_plots(sp_list,  ncol = length(qc_features)) +
      patchwork::plot_annotation(
        title = paste("QC summary —", paste(unique(object[[group_by]]), collapse = ", ")),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      )
    print(final_plot)
    return(invisible(final_plot))
  } else {
    # 没有 patchwork 就顺序打印
    message("（未安装 patchwork，按序逐图打印）")
    for (p in vln_list) print(p)
    for (p in sp_list)  print(p)
    invisible(NULL)
  }
}





suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(readr); library(tibble)
  has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
})

# —— 若这三个函数已在环境中，可忽略这段；否则自动定义一份 —— #
if (!exists("PlotViolinQC", mode = "function")) {
  .get_xy_df <- function(md) {
    stopifnot(all(c("imagecol","imagerow") %in% colnames(md)))
    dplyr::mutate(md,
                  imagecol_scaled = dplyr::coalesce(.data$imagecol_scaled, .data$imagecol),
                  imagerow_scaled = dplyr::coalesce(.data$imagerow_scaled, .data$imagerow)
    )
  }
  PlotViolinQC <- function(object, qc_feature, group_by = "orig.ident") {
    if (!qc_feature %in% colnames(object@meta.data)) {
      stop("QC指标 '", qc_feature, "' 在 meta.data 中未找到。")
    }
    VlnPlot(object, features = qc_feature, group.by = group_by, pt.size = 0) +
      NoLegend() + labs(title = qc_feature, x = "")
  }
  PlotSpatialQC <- function(object, qc_feature, ptsize = 1.6,
                            palette = c("#FFFFD4","#FEE391","#FEC44F","#FE9929","#D95F0E"),
                            legend_title = NULL) {
    if (!qc_feature %in% colnames(object@meta.data)) {
      stop("在 meta.data 找不到 QC 指标 '", qc_feature, "'")
    }
    df <- .get_xy_df(object@meta.data) |>
      dplyr::select(barcode, imagecol_scaled, imagerow_scaled, !!qc_feature)
    if (!exists("PlotExpressionV2", mode = "function"))
      stop("需要你管线中的 PlotExpressionV2()。")
    PlotExpressionV2(
      barcodes     = df,
      Gene         = qc_feature,
      ptsize       = ptsize,
      shape        = "circle",
      colors       = palette,
      legend_title = legend_title %||% qc_feature
    )
  }
  PlotQC_Summary <- function(object, group_by = "orig.ident",
                             qc_features = c("nCount_Spatial","nFeature_Spatial","percent.mt")) {
    if ("percent.mt" %in% qc_features && !"percent.mt" %in% colnames(object@meta.data)) {
      message("• 计算 percent.mt …"); object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
    }
    qc_features <- qc_features[qc_features %in% colnames(object@meta.data)]
    stopifnot(length(qc_features) > 0)
    vln_list <- lapply(qc_features, function(f) PlotViolinQC(object, f, group_by = group_by))
    sp_list  <- lapply(qc_features, function(f) PlotSpatialQC(object, f))
    if (has_patchwork) {
      p <- patchwork::wrap_plots(vln_list, ncol = length(qc_features)) /
        patchwork::wrap_plots(sp_list,  ncol = length(qc_features)) +
        patchwork::plot_annotation(
          title = paste("QC summary —", paste(unique(object[[group_by]]), collapse = ", ")),
          theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
        )
      print(p); invisible(p)
    } else {
      for (pp in vln_list) print(pp); for (pp in sp_list) print(pp); invisible(NULL)
    }
  }
}

# —— 统一计算 QC（稳健，优先 counts；补 detected_rate & percent.mt）—— #
add_basic_qc <- function(obj, assay = "Spatial", mito_pattern = "^MT-") {
  DefaultAssay(obj) <- assay
  # counts 层最稳（v5 优先 layer=；兼容老版 slot=）
  C <- tryCatch(GetAssayData(obj, assay = assay, layer = "counts"),
                error = function(e) tryCatch(GetAssayData(obj, assay = assay, slot = "counts"), error = function(e2) NULL))
  if (is.null(C)) stop("无法获取 ", assay, " 的 counts。")
  if (!inherits(C, "dgCMatrix")) C <- as(C, "dgCMatrix")
  obj$nCount_Spatial  <- Matrix::colSums(C)
  obj$nFeature_Spatial<- Matrix::colSums(C > 0)
  obj$detected_rate   <- as.numeric(Matrix::colSums(C > 0)) / nrow(C)
  if (!"percent.mt" %in% colnames(obj@meta.data)) {
    mito_idx <- grepl(mito_pattern, rownames(C))
    if (any(mito_idx)) {
      obj$percent.mt <- 100 * Matrix::colSums(C[mito_idx, , drop = FALSE]) / pmax(Matrix::colSums(C), 1)
    } else {
      obj$percent.mt <- 0
    }
  }
  obj
}





# —— 文本摘要（横向对比表使用）—— #
.quick_summary_row <- function(obj, tag) {
  tibble(
    set = tag,
    n_spots   = ncol(obj),
    med_nUMI  = median(obj$nCount_Spatial,   na.rm = TRUE),
    p90_nUMI  = quantile(obj$nCount_Spatial, .90,   na.rm = TRUE),
    med_nGene = median(obj$nFeature_Spatial, na.rm = TRUE),
    p90_nGene = quantile(obj$nFeature_Spatial,.90,  na.rm = TRUE),
    med_detect= median(obj$detected_rate,    na.rm = TRUE)
  )
}




run_hd_qc <- function(path_outs,
                      sample_id  = "sample",
                      modes      = c("all"),      # "cellbin","bin8","bin16","all"
                      outdir     = file.path(path_outs, "Analysis_Results_HD_QC"),
                      make_plots = TRUE,
                      save_pdf   = TRUE) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  if ("all" %in% modes) modes <- c("cellbin","bin8","bin16")
  
  objs <- list()
  
  # 1) 读入
  if ("cellbin" %in% modes) {
    cat("\n[Load] Cell-bin …\n")
    obj_cb <- GenerateCellbinSeurat(path_outs = path_outs, sample_id = paste0(sample_id, "_cellbin"))
    obj_cb <- add_basic_qc(obj_cb)  # 补齐 detected_rate/percent.mt
    objs$cellbin <- obj_cb
  }
  if ("bin8" %in% modes) {
    cat("\n[Load] Bin 8µm …\n")
    obj_b8 <- GenerateHDSeurat(path_outs = path_outs, bin_size = "008um")
    obj_b8 <- add_basic_qc(obj_b8)
    objs$bin8 <- obj_b8
  }
  if ("bin16" %in% modes) {
    cat("\n[Load] Bin 16µm …\n")
    obj_b16 <- GenerateHDSeurat(path_outs = path_outs, bin_size = "016um")
    obj_b16 <- add_basic_qc(obj_b16)
    objs$bin16 <- obj_b16
  }
  
  # 2) 控台摘要 + 汇总表
  for (nm in names(objs)) {
    o <- objs[[nm]]
    cat("\n=== ", toupper(nm), " ===\n", sep = "")
    cat("Spots:", ncol(o), " | Genes:", nrow(o), "\n")
    cat("nUMI   :  median=", round(median(o$nCount_Spatial,   na.rm = TRUE)),
        " [IQR ", round(quantile(o$nCount_Spatial,   .25, na.rm = TRUE)), "–",
        round(quantile(o$nCount_Spatial,   .75, na.rm = TRUE)), "]\n", sep = "")
    cat("nGene  :  median=", round(median(o$nFeature_Spatial, na.rm = TRUE)),
        " [IQR ", round(quantile(o$nFeature_Spatial, .25, na.rm = TRUE)), "–",
        round(quantile(o$nFeature_Spatial, .75, na.rm = TRUE)), "]\n", sep = "")
    cat("Detect%:  median=", round(100*median(o$detected_rate, na.rm = TRUE), 1), "%\n", sep = "")
  }
  summary_tbl <- dplyr::bind_rows(lapply(names(objs), function(nm) .quick_summary_row(objs[[nm]], nm)))
  readr::write_csv(summary_tbl, file.path(outdir, "qc_summary_across_resolutions.csv"))
  print(summary_tbl)
  
  # 3) QC 图：改为使用 PlotQC_Summary（你管线风格：Vln + PlotExpressionV2 空间）
  if (make_plots) {
    if (save_pdf) pdf(file.path(outdir, "qc_summary_plots.pdf"), width = 10, height = 8)
    for (nm in names(objs)) {
      cat("\n[Plot] ", nm, " …\n", sep = "")
      o <- objs[[nm]]
      # 三个核心指标：nCount_Spatial / nFeature_Spatial / percent.mt
      PlotQC_Summary(o, group_by = "orig.ident",
                     qc_features = c("nCount_Spatial","nFeature_Spatial","percent.mt"))
      # 如需检测率的空间图，再加一句：
      print( PlotSpatialQC(o, "detected_rate", legend_title = "Detected rate") )
    }
    if (save_pdf) dev.off()
  }
  
  # 4) 自动建议（阈值基于 FFPE/HD 的经验线；你可按需微调）
  thr <- list(
    bin8    = list(nGene = 500, detect = 0.25),
    bin16   = list(nGene = 350, detect = 0.20),
    cellbin = list(nGene = 200, detect = 0.12)
  )
  adv <- function(tbl) {
    has <- function(name) name %in% tbl$set
    msg <- c()
    if (has("bin8")) {
      r <- tbl[tbl$set=="bin8",]
      if (r$med_nGene >= thr$bin8$nGene && r$med_detect >= thr$bin8$detect) {
        msg <- c(msg, "Bin8：可用于亚群/状态（marker/GSVA/Niche/Banksy）精细分析。")
      } else {
        msg <- c(msg, sprintf("Bin8：偏稀（median nGene=%.0f, detected=%.1f%%），更适合辅助，不建议直接细颗粒度。", r$med_nGene, 100*r$med_detect))
      }
    }
    if (has("bin16")) {
      r16 <- tbl[tbl$set=="bin16",]
      line <- sprintf("Bin16：median nGene=%.0f, detected=%.1f%%。", r16$med_nGene, 100*r16$med_detect)
      if (r16$med_nGene >= thr$bin16$nGene && r16$med_detect >= thr$bin16$detect) {
        if (has("bin8")) {
          r8 <- tbl[tbl$set=="bin8",]
          if (r16$med_nGene >= 1.2 * r8$med_nGene) {
            line <- paste0(line, "（显著高于 Bin8）→ 推荐先用 Bin16 做**大群**映射，再用 Bin8 做细分。")
          } else line <- paste0(line, "→ 可与 Bin8/Cell-bin 搭配使用。")
        } else line <- paste0(line, "→ 可做大群参考/映射。")
      } else line <- paste0(line, "→ 偏稀，谨慎用于大群。")
      msg <- c(msg, line)
    }
    if (has("cellbin")) {
      rc <- tbl[tbl$set=="cellbin",]
      msg <- c(msg, sprintf(
        "Cell-bin：median nGene=%.0f, detected=%.1f%% → 推荐用于**边界/侵袭带**（可叠加轻度 KNN 平滑），配合 Niche/BANKSY。",
        rc$med_nGene, 100*rc$med_detect))
    }
    paste(msg, collapse = "\n")
  }
  recommendation <- adv(summary_tbl)
  writeLines(recommendation)
  readr::write_file(recommendation, file.path(outdir, "resolution_recommendation.txt"))
  
  invisible(list(objects = objs, summary = summary_tbl,
                 recommendation = recommendation, outdir = outdir))
}

# ==== 数据导入与生成 ====






make_images_tibble <- function(PATH) {
  img_path <- file.path(PATH, "spatial", "tissue_lowres_image.png")
  if (!file.exists(img_path)) stop("Low resolution image not found at: ", img_path)
  image <- png::readPNG(img_path)
  grobs <- grid::rasterGrob(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
  tibble::tibble(
    Path   = as.character(PATH),
    grob   = list(grobs),
    height = nrow(image),
    width  = ncol(image)
  )
}





#
GenerateSampleData <- function(PATH, size = "008um") {
  images_tibble <- make_images_tibble(PATH)
  
  bin_um <- .parse_bin_um(size)
  bdir   <- file.path(PATH, "binned_outputs", paste0("square_", size))
  if (!dir.exists(bdir)) stop("Binned dir not found: ", bdir)
  
  tpos_path   <- file.path(bdir, "spatial", "tissue_positions.parquet")
  scales_path <- file.path(bdir, "spatial", "scalefactors_json.json")
  if (!file.exists(tpos_path))   stop("Missing: ", tpos_path)
  if (!file.exists(scales_path)) stop("Missing: ", scales_path)
  
  tpos   <- arrow::read_parquet(tpos_path)
  if (ncol(tpos) >= 6) {
    colnames(tpos)[1:6] <- c("barcode","tissue","row","col","imagerow","imagecol")
    tpos <- tpos[, 1:6]
  } else stop("Unexpected columns in tissue_positions.parquet")
  
  scales <- rjson::fromJSON(file = scales_path)
  s_low   <- as.numeric(scales$tissue_lowres_scalef %||% NA_real_)
  s_hires <- as.numeric(scales$tissue_hires_scalef  %||% NA_real_)
  
  img_row <- images_tibble$height[1]
  img_col <- images_tibble$width[1]
  
  bcs <- tpos %>%
    mutate(
      imagerow_scaled       = imagerow * s_low,
      imagecol_scaled       = imagecol * s_low,
      imagerow_scaled_round = round(imagerow * s_low),
      imagecol_scaled_round = round(imagecol * s_low),
      imagerow_hires        = imagerow * s_hires,
      imagecol_hires        = imagecol * s_hires,
      tissue                = as.factor(tissue),
      height                = img_row,
      width                 = img_col
    )

  clus_csv <- file.path(bdir, "analysis", "clustering", "gene_expression_graphclust", "clusters.csv")
  umap_csv <- file.path(bdir, "analysis", "umap", "gene_expression_2_components", "projection.csv")
  
  if (file.exists(clus_csv)) {
    clusters <- read.csv(clus_csv, stringsAsFactors = FALSE)
    # 1) 识别并改名条码列 -> "Barcode"
    bc_col <- intersect(c("Barcode","barcode","Cell.ID","cell_id"), names(clusters))
    if (length(bc_col) == 0) stop("clusters.csv: cannot find Barcode-like column. Got: ", paste(names(clusters), collapse=", "))
    names(clusters)[names(clusters) == bc_col[1]] <- "Barcode"
    clusters$Barcode <- as.character(clusters$Barcode)
    # 2) 识别聚类列 -> "GraphCluster"
    clus_col <- intersect(c("GraphCluster","graphclust","Cluster","cluster"), names(clusters))
    if (length(clus_col) > 0) {
      names(clusters)[names(clusters) == clus_col[1]] <- "GraphCluster"
    }
    # 去重（防止一对多）
    clusters <- clusters[, unique(c("Barcode","GraphCluster")), drop = FALSE] %>% distinct()
    bcs <- left_join(bcs, clusters, by = c("barcode" = "Barcode"))
  }
  
  if (file.exists(umap_csv)) {
    umap <- read.csv(umap_csv, stringsAsFactors = FALSE)
    # 条码列 -> "Barcode"
    bc_col <- intersect(c("Barcode","barcode","Cell.ID","cell_id"), names(umap))
    if (length(bc_col) == 0) stop("projection.csv: cannot find Barcode-like column. Got: ", paste(names(umap), collapse=", "))
    names(umap)[names(umap) == bc_col[1]] <- "Barcode"
    # UMAP 两列 -> "UMAP_1","UMAP_2"
    if (!all(c("UMAP_1","UMAP_2") %in% names(umap))) {
      num_cols <- names(umap)[sapply(umap, is.numeric)]
      if (length(num_cols) < 2) stop("projection.csv lacks two numeric columns for UMAP")
      names(umap)[names(umap) == num_cols[1]] <- "UMAP_1"
      names(umap)[names(umap) == num_cols[2]] <- "UMAP_2"
    }
    umap <- umap[, c("Barcode","UMAP_1","UMAP_2"), drop = FALSE] %>% distinct()
    bcs  <- left_join(bcs, umap, by = c("barcode" = "Barcode"))
  }
  
  bcs <- as.data.frame(bcs)
  if (anyDuplicated(bcs$barcode)) warning("Duplicated barcodes in bcs.")
  rownames(bcs) <- bcs$barcode
  
  list(
    images_tibble = images_tibble,
    bcs           = bcs,
    scales        = scales,
    bin_um        = bin_um,
    bin_area_um2  = bin_um^2
  )
}






GenerateHDSeurat <- function(path_outs, bin_size = "008um") {
  tmp    <- GenerateSampleData(PATH = path_outs, size = bin_size)
  bcs    <- tmp$bcs
  scales <- tmp$scales
  bin_um <- tmp$bin_um
  
  # 读该 bin 的表达矩阵
  h5_path <- file.path(path_outs, "binned_outputs", paste0("square_", bin_size), "filtered_feature_bc_matrix.h5")
  if (!file.exists(h5_path)) stop("Could not find counts h5 at: ", h5_path)
  counts <- Seurat::Read10X_h5(h5_path)
  
  sobj <- Seurat::CreateSeuratObject(counts, assay = "Spatial")
  
  # 对齐 meta
  common <- intersect(colnames(sobj), rownames(bcs))
  if (length(common) == 0) {
    stop("counts 与元数据无交集：\n- counts 例子：", paste(head(colnames(sobj), 3), collapse=", "),
         "\n- bcs 例子：", paste(head(rownames(bcs), 3), collapse=", "),
         "\n请确认二者条码来源一致（都是 ", file.path("binned_outputs", paste0("square_", bin_size)), " 下的文件）。")
  }
  if (length(common) < ncol(sobj)) {
    warning("Dropping ", ncol(sobj) - length(common), " cells missing coordinates in bcs.")
    sobj <- subset(sobj, cells = common)
  }
  sobj <- Seurat::AddMetaData(sobj, bcs[colnames(sobj), , drop = FALSE])
  
  # 将 scales 与 bin 信息内置到 @misc，并冗余 microns_per_pixel 到 meta
  sobj@misc$spatial_scales     <- scales
  sobj@misc$microns_per_pixel  <- as.numeric(scales$microns_per_pixel %||% NA_real_)
  sobj@misc$bin_um             <- bin_um
  sobj@misc$bin_area_um2       <- bin_um^2
  
  if (!"microns_per_pixel" %in% colnames(sobj@meta.data))
    sobj$microns_per_pixel <- sobj@misc$microns_per_pixel
  
  cat(sprintf("--- HD Binned (%s) Seurat object created. Scales embedded. ---\n", bin_size))
  sobj
}






suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(arrow)
  library(rjson)
  library(png)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b
GenerateCellbinSampleData <- function(path_outs) {
  cat("--- 准备 Cellbin 单细胞元数据 (v15) ---\n")
  
  path_outs_norm <- normalizePath(path_outs, winslash = "/", mustWork = TRUE)
  image_dir      <- file.path(path_outs_norm, "spatial")
  lowres_image   <- file.path(image_dir, "tissue_lowres_image.png")
  
  mappings_path   <- file.path(path_outs_norm, "barcode_mappings.parquet")
  coords_2um_path <- file.path(path_outs_norm, "binned_outputs", "square_002um", "spatial", "tissue_positions.parquet")
  scales_path     <- file.path(path_outs_norm, "segmented_outputs", "spatial", "scalefactors_json.json")
  
  if (!file.exists(lowres_image)) stop("缺少低清图像: ", lowres_image)
  if (!all(file.exists(mappings_path, coords_2um_path, scales_path))) {
    stop("缺少必要文件：\n- ", mappings_path, "\n- ", coords_2um_path, "\n- ", scales_path)
  }
  
  # 低清图尺寸
  img_low    <- png::readPNG(lowres_image)
  height_low <- nrow(img_low); width_low <- ncol(img_low)
  
  # 比例因子
  scales <- rjson::fromJSON(file = scales_path)
  pp_full  <- as.numeric(scales$microns_per_pixel %||% NA_real_)   # full-res: µm/px
  s_low    <- as.numeric(scales$tissue_lowres_scalef %||% NA_real_)
  s_hires  <- as.numeric(scales$tissue_hires_scalef  %||% NA_real_)
  # 不同分辨率下像素的 µm/px（注意：坐标= full-res * scale，因此 px 尺寸 = pp_full / scale）
  pp_low   <- pp_full / s_low
  pp_hires <- pp_full / s_hires
  
  # 读映射与坐标
  barcode_mappings <- arrow::read_parquet(mappings_path)
  coords_2um_df    <- arrow::read_parquet(coords_2um_path)
  colnames(coords_2um_df) <- c("barcode","tissue","row","col","imagerow","imagecol")
  
  # 自动识别 join 键：square_XXXum
  join_col <- grep("^square_\\d+um$", colnames(barcode_mappings), value = TRUE)
  if (length(join_col) == 0) stop("barcode_mappings 中未找到 'square_XXXum' 连接键。")
  join_col <- join_col[1]
  
  # 解析 bin 大小（µm）
  bin_um_from_col  <- suppressWarnings(as.numeric(sub("^square_(\\d+)um$", "\\1", join_col)))
  bin_um_from_path <- suppressWarnings(as.numeric(gsub(".*square_([0-9]+)um.*", "\\1", coords_2um_path)))
  bin_um <- dplyr::coalesce(bin_um_from_col, bin_um_from_path, 2)
  if (is.na(bin_um) || bin_um <= 0) { warning("未能解析 bin 大小，默认 2 µm。"); bin_um <- 2 }
  bin_area_um2 <- bin_um^2
  cat(sprintf("  - bin: %d µm (单格 %.1f µm²)\n", bin_um, bin_area_um2))
  
  # 仅细胞内小格
  cell_bins <- barcode_mappings %>% dplyr::filter(in_cell == TRUE, !is.na(cell_id))
  
  # 连接 2µm 坐标（square_XXXum -> barcode）
  cell_bins <- dplyr::left_join(
    cell_bins, coords_2um_df,
    by = setNames("barcode", join_col)
  )
  
  # 聚合至 cell_id
  bcs <- cell_bins %>%
    dplyr::group_by(cell_id) %>%
    dplyr::summarise(
      imagerow     = mean(imagerow, na.rm = TRUE),
      imagecol     = mean(imagecol, na.rm = TRUE),
      tissue       = max(tissue, na.rm = TRUE),
      num_barcodes = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(barcode = cell_id) %>%
    dplyr::mutate(
      # 物理面积（µm²）与 bin 口径
      area_um2       = num_barcodes * bin_area_um2,
      cell_area_bins = num_barcodes,
      
      # —— 像素等效个数（不是逐像素统计，但物理等效，可复现）——
      px_equiv_fullres = round(area_um2 / (pp_full^2)),
      px_equiv_hires   = round(area_um2 / (pp_hires^2)),
      px_equiv_lowres  = round(area_um2 / (pp_low^2)),
      
      # 坐标缩放（与低清/高清图一致）
      imagerow_scaled = imagerow * s_low,
      imagecol_scaled = imagecol * s_low,
      imagerow_hires  = imagerow * s_hires,
      imagecol_hires  = imagecol * s_hires,
      
      # 便捷字段
      microns_per_pixel = pp_full
    ) %>% as.data.frame()

  # 附加画布尺寸（低清）
  bcs$height_lowres <- height_low
  bcs$width_lowres  <- width_low
  rownames(bcs) <- bcs$barcode

  cat("  - 合并 Space Ranger analysis（clusters / UMAP）...\n")
  cluster_path <- file.path(path_outs_norm, "segmented_outputs", "analysis",
                            "clustering", "gene_expression_graphclust", "clusters.csv")
  umap_path    <- file.path(path_outs_norm, "segmented_outputs", "analysis",
                            "umap", "gene_expression_2_components", "projection.csv")
  
  # Clusters
  if (file.exists(cluster_path)) {
    clusters <- read.csv(cluster_path, stringsAsFactors = FALSE)
    bc_col <- intersect(c("Barcode","barcode","Cell.ID","cell_id"), colnames(clusters))
    if (length(bc_col) == 0) stop("clusters.csv 未找到条形码列。")
    clusters <- clusters %>%
      dplyr::rename(Barcode = !!bc_col[1]) %>%
      dplyr::mutate(Barcode = as.character(Barcode))
    clus_col <- intersect(c("Cluster","cluster","GraphCluster","graphclust"), colnames(clusters))
    if (length(clus_col) > 0) clusters <- dplyr::rename(clusters, GraphCluster = !!clus_col[1])
    bcs <- dplyr::left_join(bcs, clusters, by = c("barcode" = "Barcode"))
  } else {
    message("  · 未找到 clusters.csv，跳过聚类合并。")
  }
  
  # UMAP
  if (file.exists(umap_path)) {
    umap_df <- read.csv(umap_path, stringsAsFactors = FALSE)
    bc_col <- intersect(c("Barcode","barcode","Cell.ID","cell_id"), colnames(umap_df))
    if (length(bc_col) == 0) stop("projection.csv 未找到条形码列。")
    num_cols <- names(umap_df)[sapply(umap_df, is.numeric)]
    if (length(intersect(c("UMAP_1","UMAP_2"), colnames(umap_df))) >= 2) {
      # 已有标准名
    } else if (length(num_cols) >= 2) {
      umap_df <- dplyr::rename(umap_df, UMAP_1 = !!num_cols[1], UMAP_2 = !!num_cols[2])
    } else {
      stop("projection.csv 未识别到 2 列数值型 UMAP 坐标。")
    }
    umap_df <- umap_df %>% dplyr::rename(Barcode = !!bc_col[1])
    bcs <- dplyr::left_join(bcs, umap_df[, c("Barcode","UMAP_1","UMAP_2")],
                            by = c("barcode" = "Barcode"))
    rownames(bcs) <- bcs$barcode 
  } else {
    message("  · 未找到 projection.csv，跳过 UMAP 合并。")
  }
  
  list(
    bcs = bcs,
    scales = scales,
    bin_um = bin_um,
    bin_area_um2 = bin_area_um2
  )
}





GenerateCellbinSeurat <- function(path_outs, sample_id) {
  tmp <- GenerateCellbinSampleData(path_outs = path_outs)
  
  bcs <- tmp$bcs %>% dplyr::filter(tissue == 1)
  
  h5 <- file.path(path_outs, "segmented_outputs", "filtered_feature_cell_matrix.h5")
  if (!file.exists(h5)) stop("缺少表达矩阵: ", h5)
  counts <- Read10X_h5(h5)
  rownames(bcs) <- bcs$barcode 
  common_cells <- intersect(colnames(counts), rownames(bcs))
  if (length(common_cells) == 0) stop("counts 与元数据无交集，请检查条形码命名。")
  
  counts_aligned <- counts[, common_cells, drop = FALSE]
  meta_aligned   <- bcs[common_cells, , drop = FALSE]
  
  sobj <- CreateSeuratObject(
    counts   = counts_aligned,
    meta.data= meta_aligned,
    assay    = "Spatial",
    project  = sample_id
  )
  
  # 附加比例与网格信息（安全放在 @misc）
  sobj@misc$scales             <- tmp$scales
  sobj@misc$microns_per_pixel  <- as.numeric(tmp$scales$microns_per_pixel %||% NA_real_)
  sobj@misc$bin_um             <- tmp$bin_um
  sobj@misc$bin_area_um2       <- tmp$bin_area_um2
  
  # 冗余一列（便于在 meta 里直接使用）
  if (!"microns_per_pixel" %in% colnames(sobj@meta.data))
    sobj$microns_per_pixel <- sobj@misc$microns_per_pixel
  
  cat("--- Cellbin Seurat 对象创建成功（含 analysis 合并、物理面积与像素等效） ---\n")
  sobj
}






suppressPackageStartupMessages({
  library(dplyr); library(sf); library(concaveman)
  library(ggplot2); library(png); library(rjson)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

suppressPackageStartupMessages({
  library(dplyr); library(sf); library(concaveman)
  library(ggplot2); library(png); library(rjson)
})





`%||%` <- function(a,b) if (!is.null(a)) a else b

# ==== 调色与空间分布 ====
# 最终修改版: draw_10x_tissue_outline
# (兼容 Cellbin 路径并重新加入比例尺功能)
# ==== 表达可视化与QC报告 ====

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(sf); library(concaveman)
  library(ggplot2); library(png); library(rjson)
})





`%||%` <- function(a, b) if (!is.null(a)) a else b
`%||%` <- function(a,b) if (!is.null(a)) a else b
draw_10x_tissue_outline <- function(
    object,
    outs_path,
    save_prefix = "HE_tissue_outline",
    resolution  = c("lowres", "hires"),
    # 原参数
    concavity   = 2.2,
    smooth_px   = 1.4,
    simplify_px = 0.8,
    keep        = c("largest", "all"),
    line_halo   = TRUE,
    line_width  = 1.0,
    line_color  = "black",
    halo_color  = "white",
    # --- 比例尺 ---
    add_scale_bar        = FALSE,
    scale_bar_length_mm  = 1,
    scale_bar_n_segments = 4,
    scale_bar_color      = "#FFBF00",
    scale_bar_thickness_px = 8,
    scale_bar_thickness_um = NULL,              # 物理厚度（µm），会覆盖 *_px
    scale_bar_position   = c("bottomright","bottomleft","topright","topleft"),
    scale_bar_margin_frac = c(0.04, 0.05),
    scale_bar_margin_px   = NULL,               # NEW：绝对像素边距（优先于 frac；可与 frac 叠加）
    scale_bar_text       = NULL,
    scale_bar_text_size  = 3.2,
    scale_bar_text_color = "black",
    scale_bar_text_face  = "bold",
    scale_bar_text_offset_px = 10,
    scale_bar_gap_px     = 1,
    # --- 视野裁剪 ---
    clip_to_data   = FALSE,
    clip_margin_frac = c(0.02, 0.02),
    clip_margin_px   = NULL,                    # NEW：绝对像素边距（优先于 frac；可与 frac 叠加）
    # --- 轮廓整体微调（像素级）---
    outline_bias_px = 0,                        # NEW：>0 外扩；<0 内缩；0 不调整
    # --- 线条风格（更柔和） ---
    line_join = c("round","mitre","bevel"),     # NEW
    line_end  = c("round","butt","square"),     # NEW
    # 保存
    save_plots = TRUE,
    width_in = 8, height_in = 7, dpi = 300
){
  resolution <- match.arg(resolution)
  keep       <- match.arg(keep)
  scale_bar_position <- match.arg(scale_bar_position)
  line_join <- match.arg(line_join)
  line_end  <- match.arg(line_end)
  
  `%||%` <- function(a,b) if (is.null(a)) b else a
  get_scales <- function(object, outs_path){
    sc <- object@misc$scales %||% object@misc$spatial_scales %||% NULL
    if (is.null(sc)) {
      cand_paths <- c(
        file.path(outs_path, "segmented_outputs", "spatial", "scalefactors_json.json"),
        file.path(outs_path, "spatial", "scalefactors_json.json")
      )
      sf_json <- cand_paths[file.exists(cand_paths)][1]
      if (is.na(sf_json)) stop("No scalefactors_json.json in spatial/ or segmented_outputs/spatial/.")
      message("Loading scalefactors: ", sf_json)
      sc <- rjson::fromJSON(file = sf_json)
    }
    list(
      mpp_full = as.numeric(sc$microns_per_pixel %||% NA_real_),
      s_low    = as.numeric(sc$tissue_lowres_scalef %||% NA_real_),
      s_hires  = as.numeric(sc$tissue_hires_scalef  %||% NA_real_)
    )
  }
  scales <- get_scales(object, outs_path)
  if (any(is.na(unlist(scales)))) stop("Missing microns_per_pixel / tissue_lowres_scalef / tissue_hires_scalef.")
  
  img_file <- if (resolution=="lowres")
    file.path(outs_path, "spatial", "tissue_lowres_image.png") else
      file.path(outs_path, "spatial", "tissue_hires_image.png")
  if (!file.exists(img_file)) stop("Missing HE image: ", img_file)
  img <- png::readPNG(img_file); W <- ncol(img); H <- nrow(img)
  
  meta <- object@meta.data
  has_scaled <- all(c("imagecol_scaled","imagerow_scaled") %in% colnames(meta))
  has_hires  <- all(c("imagecol_hires","imagerow_hires") %in% colnames(meta))
  has_raw    <- all(c("imagecol","imagerow") %in% colnames(meta))
  
  if (resolution == "lowres") {
    if (has_scaled) { meta$X_plot <- meta$imagecol_scaled; meta$Y_plot <- H - meta$imagerow_scaled
    } else if (has_raw) { meta$X_plot <- meta$imagecol * scales$s_low; meta$Y_plot <- H - (meta$imagerow * scales$s_low)
    } else stop("No columns for lowres mapping.")
  } else {
    if (has_hires) { meta$X_plot <- meta$imagecol_hires; meta$Y_plot <- H - meta$imagerow_hires
    } else if (has_scaled) { fac <- scales$s_hires / scales$s_low; meta$X_plot <- meta$imagecol_scaled * fac; meta$Y_plot <- H - (meta$imagerow_scaled * fac)
    } else if (has_raw) { meta$X_plot <- meta$imagecol * scales$s_hires; meta$Y_plot <- H - (meta$imagerow * scales$s_hires)
    } else stop("No columns for hires mapping.")
  }
  
  if (!"tissue" %in% colnames(meta)) { warning("No 'tissue' column; treat all as in-tissue."); meta$tissue <- 1L }
  
  md_filtered <- meta |> dplyr::filter(tissue==1L, is.finite(X_plot), is.finite(Y_plot)) |> dplyr::select(X_plot, Y_plot)
  if (nrow(md_filtered) < 10) stop("Too few in-tissue points to build outline.")
  
  pts  <- sf::st_as_sf(md_filtered, coords=c("X_plot","Y_plot"), crs=sf::NA_crs_)
  hull <- concaveman::concaveman(pts, concavity=concavity, length_threshold=0)
  poly <- sf::st_cast(sf::st_make_valid(sf::st_as_sf(hull)), "POLYGON", do_split=TRUE)
  if (keep=="largest" && nrow(poly)>1) { poly$area <- as.numeric(sf::st_area(poly)); poly <- poly[which.max(poly$area), , drop=FALSE] }
  poly <- poly |>
    sf::st_buffer(smooth_px) |>
    sf::st_buffer(-smooth_px) |>
    sf::st_simplify(dTolerance=simplify_px) |>
    sf::st_make_valid()
  
  # NEW: 整体外扩/内缩若干像素，修正你看到的“偏移”或让轮廓更贴视觉边缘
  if (is.finite(outline_bias_px) && outline_bias_px != 0) {
    poly <- sf::st_make_valid(sf::st_buffer(poly, outline_bias_px))
  }
  
  # 视窗范围（支持绝对像素边距）
  if (isTRUE(clip_to_data)) {
    bb <- sf::st_bbox(poly)
    # frac 边距
    xpad_f <- (bb$xmax - bb$xmin) * (if (length(clip_margin_frac)>=1) clip_margin_frac[1] else 0.02)
    ypad_f <- (bb$ymax - bb$ymin) * (if (length(clip_margin_frac)>=2) clip_margin_frac[2] else xpad_f)
    # px 边距
    if (!is.null(clip_margin_px)) {
      px_x <- if (length(clip_margin_px)>=1) clip_margin_px[1] else clip_margin_px
      px_y <- if (length(clip_margin_px)>=2) clip_margin_px[2] else px_x
    } else { px_x <- 0; px_y <- 0 }
    xlim <- c(max(0, bb$xmin - xpad_f - px_x), min(W, bb$xmax + xpad_f + px_x))
    ylim <- c(max(0, bb$ymin - ypad_f - px_y), min(H, bb$ymax + ypad_f + px_y))
  } else {
    xlim <- c(0, W); ylim <- c(0, H)
  }
  Wv <- diff(xlim); Hv <- diff(ylim)
  
  p <- ggplot2::ggplot() +
    ggplot2::annotation_raster(img, xmin=0, xmax=W, ymin=0, ymax=H, interpolate=TRUE) +
    { if (isTRUE(line_halo)) ggplot2::geom_sf(data=poly, fill=NA, color=halo_color, linewidth=line_width*1.6,
                                              linejoin=line_join, lineend=line_end) else NULL } +
    ggplot2::geom_sf(data=poly, fill=NA, color=line_color, linewidth=line_width,
                     linejoin=line_join, lineend=line_end) +
    ggplot2::coord_sf(xlim=xlim, ylim=ylim, expand=FALSE, default_crs=sf::NA_crs_) +
    ggplot2::theme_void()
  
  if (isTRUE(add_scale_bar)) {
    s_use <- if (resolution=="hires") scales$s_hires else scales$s_low
    px_per_um <- s_use / scales$mpp_full
    bar_px <- as.integer(max(10, round(scale_bar_length_mm * 1000 * px_per_um)))
    thickness_px <- if (!is.null(scale_bar_thickness_um) && is.finite(scale_bar_thickness_um) && scale_bar_thickness_um>0) {
      as.integer(max(1, round(scale_bar_thickness_um * px_per_um)))
    } else as.integer(max(1, round(scale_bar_thickness_px)))
    
    mx_f <- if (length(scale_bar_margin_frac)>=1) scale_bar_margin_frac[1] else 0.04
    my_f <- if (length(scale_bar_margin_frac)>=2) scale_bar_margin_frac[2] else mx_f
    mx_px <- as.integer(round(mx_f * Wv))
    my_px <- as.integer(round(my_f * Hv))
    # NEW：绝对像素边距（在 frac 基础上叠加）
    if (!is.null(scale_bar_margin_px)) {
      mx_px <- mx_px + as.integer(if (length(scale_bar_margin_px)>=1) scale_bar_margin_px[1] else scale_bar_margin_px)
      my_px <- my_px + as.integer(if (length(scale_bar_margin_px)>=2) scale_bar_margin_px[2] else (scale_bar_margin_px[1] %||% scale_bar_margin_px))
    }
    
    if (scale_bar_position %in% c("bottomright","topright")) x0 <- xlim[2] - mx_px - bar_px else x0 <- xlim[1] + mx_px
    if (scale_bar_position %in% c("bottomright","bottomleft")) y0 <- ylim[1] + my_px else y0 <- ylim[2] - my_px - thickness_px
    x0 <- as.integer(round(x0)); y0 <- as.integer(round(y0)); y1 <- y0 + thickness_px; x1 <- x0 + bar_px
    
    nseg  <- max(1L, as.integer(round(scale_bar_n_segments)))
    gappx <- max(0L, as.integer(round(scale_bar_gap_px))); if (nseg == 1L) gappx <- 0L
    total_gap_px  <- gappx * (nseg - 1L)
    fill_px_total <- max(bar_px - total_gap_px, nseg)
    seg_fill_base <- floor(fill_px_total / nseg)
    rem <- fill_px_total - seg_fill_base * nseg
    seg_fill <- rep(seg_fill_base, nseg); if (rem>0L) seg_fill[seq_len(rem)] <- seg_fill[seq_len(rem)] + 1L
    
    left <- x0; seg_list <- vector("list", nseg)
    for (i in seq_len(nseg)) {
      right <- left + seg_fill[i]
      seg_list[[i]] <- data.frame(xmin=left, xmax=right, ymin=y0, ymax=y1)
      left <- right + if (i<nseg) gappx else 0L
    }
    sb_df <- do.call(rbind, seg_list)
    p <- p + ggplot2::geom_rect(data=sb_df,
                                ggplot2::aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                                inherit.aes=FALSE, fill=scale_bar_color, color=NA)
    if (!is.null(scale_bar_text) && nzchar(scale_bar_text)) {
      text_y <- if (scale_bar_position %in% c("bottomright","bottomleft")) y1 + scale_bar_text_offset_px else y0 - scale_bar_text_offset_px
      p <- p + ggplot2::geom_text(data=data.frame(x = x0 + bar_px/2, y = text_y, lab = scale_bar_text),
                                  ggplot2::aes(x=x, y=y, label=lab),
                                  size=scale_bar_text_size, color=scale_bar_text_color,
                                  fontface=scale_bar_text_face,
                                  vjust = if (scale_bar_position %in% c("bottomright","bottomleft")) 0 else 1)
    }
  }
  
  print(p)
  if (isTRUE(save_plots)) {
    tag <- paste0("_", resolution)
    base <- tools::file_path_sans_ext(save_prefix)
    pdf_out <- paste0(base, tag, ".pdf")
    png_out <- paste0(base, tag, ".png")
    if (dir.exists(dirname(save_prefix))) {
      pdf_out <- file.path(dirname(save_prefix), paste0(basename(base), tag, ".pdf"))
      png_out <- file.path(dirname(save_prefix), paste0(basename(base), tag, ".png"))
    }
    ggplot2::ggsave(pdf_out, p, width=width_in, height=height_in, dpi=dpi, useDingbats=FALSE)
    ggplot2::ggsave(png_out, p, width=width_in, height=height_in, dpi=dpi)
    message("Saved:\n- ", pdf_out, "\n- ", png_out)
  }
  
  list(plot=p, outline=poly, canvas=list(width=W, height=H, resolution=resolution),
       view_bbox=list(xlim=xlim, ylim=ylim))
}

# ==== 标签导出与分布变体 ====








# ==== Seurat流水线与平滑 ====
#  A) Numeric renaming
# ==== RCTD/Cellbin 处理 ====

# Formats a vector of group labels into numeric-named labels (e.g., "C01_Basal").
# It orders by the first integer in each label (if found), else alphabetical.
autonumber_groups <- function(
    x,                            # character or factor vector of group labels (length = cells)
    prefix = "C",                 # leading prefix
    width  = 2,                   # zero-padding width
    fmt    = "{prefix}{num}_{label}",  # final label format
    order_by = c("numeric_then_alpha","alpha"),  # primary ordering rule
    keep_original_levels = TRUE
){
  x_chr <- as.character(x)
  order_by <- match.arg(order_by)
  
  # extract first integer per unique label
  uniq <- unique(x_chr)
  m <- regexpr("-?[0-9]+", uniq)
  key <- ifelse(m > 0, suppressWarnings(as.numeric(regmatches(uniq, m))), NA_real_)
  
  if (order_by == "numeric_then_alpha") {
    ord <- order(is.na(key), key, tolower(uniq))
  } else {
    ord <- order(tolower(uniq))
  }
  
  uniq_ord <- uniq[ord]
  # assign running numbers 1..k by order
  num <- seq_along(uniq_ord)
  
  # formatted numbers
  num_fmt <- sprintf(paste0("%0", width, "d"), num)
  
  # build new labels
  format_one <- function(prefix, num_str, label) {
    sub("\\{prefix\\}", prefix,
        sub("\\{num\\}",    num_str,
            sub("\\{label\\}", label, fmt, fixed = TRUE),
            fixed = TRUE),
        fixed = TRUE)
  }
  new_labels_ord <- mapply(format_one, prefix, num_fmt, uniq_ord, USE.NAMES = FALSE)
  
  # map back to all cells
  map_df <- data.frame(
    old = uniq_ord,
    number = num,
    number_str = num_fmt,
    new = new_labels_ord,
    stringsAsFactors = FALSE
  )
  map <- setNames(map_df$new, map_df$old)
  x_new <- unname(map[x_chr])
  
  # levels for factors
  levels_new <- if (keep_original_levels) new_labels_ord else sort(unique(x_new))
  
  list(
    x_new = x_new,
    map_df = map_df,         # data.frame: old → new (+ numbers)
    rename_map = map,        # named character: old -> new
    levels_new = levels_new  # recommended factor level order
  )
}





# ==== 一致性评估与锚定注释 ====
#  B) High-contrast palette
# ==== BANKSY/可视化工具集 ====

# Returns a *named* color vector aligned to the provided labels (vector or levels).
# "glasbey" (from pals) is default: unlimited, crisp, color-blind friendly.
palette_for_groups <- function(
    x_or_levels,
    scheme = c("glasbey","polychrome","hcl","viridis"),
    tone   = c("vivid","pastel"),   # mild post-processing choice
    seed   = 42L
){
  scheme <- match.arg(scheme)
  tone   <- match.arg(tone)
  
  # accept a vector of labels (factor/character) or a character vector of levels
  if (is.factor(x_or_levels)) {
    lv <- levels(x_or_levels)
    if (is.null(lv) || length(lv) == 0L) lv <- unique(as.character(x_or_levels))
  } else {
    lv <- unique(as.character(x_or_levels))
  }
  lv <- as.character(lv)
  n  <- length(lv)
  if (n == 0L) stop("No group levels found: please pass a non-empty vector.")
  
  # generator
  gen_cols <- function(n){
    if (scheme == "glasbey" && requireNamespace("pals", quietly = TRUE)) {
      # unlimited qualitative; excellent separability
      pals::glasbey(n)
    } else if (scheme == "polychrome" && requireNamespace("Polychrome", quietly = TRUE)) {
      set.seed(seed)
      Polychrome::createPalette(n, seedcolors = c("#ff0000","#00ff00","#0000ff"))
    } else if (scheme == "hcl" && requireNamespace("colorspace", quietly = TRUE)) {
      colorspace::qualitative_hcl(n, palette = "Dark 3")
    } else if (requireNamespace("viridisLite", quietly = TRUE)) {
      viridisLite::viridis(n) # sequential; used as robust fallback
    } else {
      grDevices::hcl(h = seq(15, 375, length.out = n + 1L)[1L:n], c = 60, l = 65)
    }
  }
  
  cols <- gen_cols(n)
  
  # optional toning (gentle, keeps contrasts)
  if (tone == "pastel") {
    if (requireNamespace("colorspace", quietly = TRUE)) {
      cols <- colorspace::lighten(cols, amount = 0.18)
      cols <- colorspace::desaturate(cols, amount = 0.10)
    }
  } else {
    # vivid: slightly deepen to print nicely
    if (requireNamespace("colorspace", quietly = TRUE)) {
      cols <- colorspace::darken(cols, amount = 0.05)
    }
  }
  
  setNames(cols, lv)
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







# ---------- helpers ----------
.fmt_num <- function(x) sub("\\.?0+$", "", format(x, trim = TRUE))
.pick_col_safe <- function(df, col) {
  v <- df[, col, drop = TRUE]
  if (is.data.frame(v)) v <- v[,1, drop = TRUE]
  v
}





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
# family: 标签第一个 "_" 或 "/" 之前的部分（若无，取整个词）
# macro : 将 family 归并到宏类（Epithelial/Immune/CAF/Acinar/Vascular/Endocrine/Neural/Interface/Other/Ambient）
parse_pdac_family_macro <- function(labels){
  lbl <- tidy_pdac_name(labels)
  
  # Normalize / to _ for consistent splitting
  lbl_norm <- sub("/", "_", lbl)
  
  # family by prefix
  fam <- ifelse(grepl("_", lbl_norm, fixed = TRUE), sub("_.*$", "", lbl_norm), lbl_norm)
  fam_lc <- tolower(fam)
  
  # normalize families to canonical names (enhanced with grepl for robustness)
  fam_norm <- dplyr::case_when(
    grepl("tumor|epithelial|cancer", fam_lc)      ~ "Tumor",
    grepl("adm|ductal|duct", fam_lc)              ~ "ADM",
    grepl("^acinar", fam_lc)                      ~ "Acinar",
    grepl("caf|icaf|mycaf|fibroblast", fam_lc)    ~ "CAF",
    grepl("tam|macrophage|mono|myeloid", fam_lc)  ~ "TAM",
    grepl("plasma|b", fam_lc)                     ~ "Plasma",
    grepl("endocrine|islet", fam_lc)              ~ "Endocrine",
    grepl("endothelial|vascular", fam_lc)         ~ "Endothelial",
    grepl("^perivascular", fam_lc)                ~ "Perivascular",
    grepl("smoothmuscle", fam_lc)                 ~ "SmoothMuscle",
    grepl("schwann|neural|neuron|glia", fam_lc)   ~ "Schwann",
    grepl("^interface", fam_lc)                   ~ "Interface",
    grepl("ambient", tolower(lbl))                ~ "Ambient",
    TRUE                                          ~ stringr::str_to_title(fam) # fallback
  )
  
  # macro mapping (enhanced with broader patterns)
  macro <- dplyr::case_when(
    grepl("Tumor|ADM|Epithelial", fam_norm, ignore.case = TRUE) ~ "Epithelial",
    grepl("Immune|TAM|Plasma|B|T|Nk|Mast|Dendritic", fam_norm, ignore.case = TRUE) ~ "Immune",
    grepl("CAF|Fibroblast", fam_norm, ignore.case = TRUE) ~ "CAF",
    grepl("Acinar", fam_norm, ignore.case = TRUE) ~ "Acinar",
    grepl("Endocrine|Islet", fam_norm, ignore.case = TRUE) ~ "Endocrine",
    grepl("Endothelial|Perivascular|Vascular|SmoothMuscle", fam_norm, ignore.case = TRUE) ~ "Vascular",
    grepl("Schwann|Neural|Glia", fam_norm, ignore.case = TRUE) ~ "Neural",
    grepl("Interface", fam_norm, ignore.case = TRUE) ~ "Interface",
    grepl("Ambient", fam_norm, ignore.case = TRUE) ~ "Ambient",
    TRUE ~ "Other"
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

# ==== BANKSY/三重生态/边界 ====










suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

PlotExpressionV2 <- function(
    object,                        # Seurat 或 data.frame
    feature,                       # 基因 或 meta.data 里的分数列名
    assay        = "Spatial",
    slot         = "data",
    source       = c("auto","meta","matrix"),
    # 显示样式
    shape        = c("square","circle"),
    spot_size    = c("auto","constant"),     # square+auto -> 真实大小
    colors       = c("grey80","yellow","red"),
    legend_title = NULL,
    # —— KNN 平滑（可选）——
    knn_smooth   = FALSE,
    k            = 8,
    method       = c("gaussian","mean"),
    sigma        = NULL,
    # —— 真实大小所需（自动从 misc 取；也可显式传标量或列名）——
    area_col     = NULL,                     # 优先：Cellbin 的 "area_um2"
    area_unit    = c("um2","lowres_px"),
    bin_um       = NULL,                     # 其次：HD 的 bin 大小（如 8）
    microns_per_pixel = NULL,                # 标量或列名；若传 Seurat 则自动取
    scale_lowres      = NULL,                # 标量或列名；若传 Seurat 则自动取
    ...
){
  # ---- 0) 参数整理 ----
  source    <- match.arg(source)
  shape     <- match.arg(shape)
  spot_size <- match.arg(spot_size)
  method    <- match.arg(method)
  area_unit <- match.arg(area_unit)
  
  # ---- 1) 取绘图数据 + 尺子 ----
  if (inherits(object, "Seurat")) {
    md <- object@meta.data
    
    # 坐标（优先 lowres）
    if (all(c("imagecol_scaled","imagerow_scaled") %in% colnames(md))) {
      x <- md$imagecol_scaled; y <- md$imagerow_scaled
    } else if (all(c("imagecol","imagerow") %in% colnames(md))) {
      x <- md$imagecol; y <- md$imagerow
    } else stop("meta.data 中找不到低清或原始坐标：imagecol(_scaled)/imagerow(_scaled).")
    
    plot_df <- data.frame(
      barcode = rownames(md),
      imagecol_scaled = x,
      imagerow_scaled = y,
      stringsAsFactors = FALSE
    )
    
    # 自动取尺子（mpp / lowres scale / bin_um）
    mpp <- microns_per_pixel %||%
      object@misc$microns_per_pixel %||%
      object@misc$spatial_scales$microns_per_pixel %||%
      object@misc$scales$microns_per_pixel %||%
      tryCatch(Images(object)[[1]]@scale.factors$microns_per_pixel, error=function(e) NULL)
    
    s_low <- scale_lowres %||%
      object@misc$spatial_scales$tissue_lowres_scalef %||%
      object@misc$scales$tissue_lowres_scalef %||%
      tryCatch(Images(object)[[1]]@scale.factors$tissue_lowres_scalef, error=function(e) NULL)
    
    if (is.null(bin_um)) {
      bin_um <- object@misc$bin_um %||%
        object@misc$spatial_scales$bin_size_um %||% NULL
    }
    
    # 自动识别 Cellbin 面积列
    if (is.null(area_col) && "area_um2" %in% colnames(md)) area_col <- "area_um2"
    
    # 构建 feature 向量（自动判定：meta优先）
    if (source == "auto") {
      if (feature %in% colnames(md)) {
        val <- as.numeric(md[[feature]]); src <- "meta"
      } else {
        DefaultAssay(object) <- assay
        M <- Seurat::GetAssayData(object, assay=assay, slot=slot)
        if (feature %in% rownames(M)) {
          val <- as.numeric(M[feature, colnames(object)])
          src <- "matrix"
        } else {
          stop("在 meta.data 和 ", assay,"/",slot," 都找不到：", feature)
        }
      }
    } else if (source == "meta") {
      stopifnot(feature %in% colnames(md))
      val <- as.numeric(md[[feature]]); src <- "meta"
    } else { # source == "matrix"
      DefaultAssay(object) <- assay
      M <- Seurat::GetAssayData(object, assay=assay, slot=slot)
      stopifnot(feature %in% rownames(M))
      val <- as.numeric(M[feature, colnames(object)]); src <- "matrix"
    }
    
    # 把可能有用的列带进 df（area_um2 等）
    if (!is.null(area_col) && area_col %in% colnames(md))
      plot_df[[area_col]] <- md[[area_col]]
    
  } else if (is.data.frame(object)) {
    plot_df <- object
    if (!all(c("imagecol_scaled","imagerow_scaled") %in% colnames(plot_df)))
      stop("data.frame 需要包含列：imagecol_scaled, imagerow_scaled")
    # 从 df 取 feature 列，若没有再报错（此分支不自动查 assay）
    if (feature %in% colnames(plot_df)) {
      val <- as.numeric(plot_df[[feature]]); src <- "meta"
    } else {
      stop("传入的是 data.frame：请先把 '", feature, "' 列加入其中（或传 Seurat 以便从矩阵读取）。")
    }
    # 尺子
    mpp  <- microns_per_pixel
    s_low<- scale_lowres
    if (is.null(mpp) || is.null(s_low))
      stop("data.frame 模式需要提供 microns_per_pixel 与 scale_lowres（标量或列名）。")
  } else {
    stop("object 必须是 Seurat 或 data.frame")
  }
  
  # ---- 2) 可选：KNN 平滑 ----
  if (knn_smooth) {
    if (!requireNamespace("FNN", quietly = TRUE)) stop("需要 FNN 包：install.packages('FNN')")
    coords <- cbind(plot_df$imagecol_scaled, plot_df$imagerow_scaled)
    nn <- FNN::get.knn(coords, k = min(k, nrow(coords)-1))
    if (is.null(sigma) && method == "gaussian") {
      sigma <- stats::median(nn$nn.dist[,1], na.rm=TRUE)
      if (!is.finite(sigma) || sigma <= 0) sigma <- 1
    }
    # 权重
    dist_mat <- cbind(0, nn$nn.dist)
    if (method == "gaussian") {
      Wrow <- exp(-(dist_mat^2)/(2*sigma^2))
    } else {
      Wrow <- matrix(1, nrow(dist_mat), ncol(dist_mat))
    }
    Wrow <- Wrow / pmax(rowSums(Wrow), 1e-12)
    # 稀疏乘法
    n <- nrow(coords)
    idx_rows <- as.vector(row(Wrow))
    idx_cols <- as.vector(cbind(matrix(seq_len(n), n, 1), nn$nn.index))
    W <- Matrix::sparseMatrix(i=idx_rows, j=idx_cols, x=as.vector(Wrow), dims=c(n,n))
    val <- as.numeric(W %*% val)
    if (is.null(legend_title)) legend_title <- paste0(feature, if (method=="gaussian") "_wknn" else "_knn", k)
  }
  
  # ---- 3) 组装绘图 df ----
  plot_df$Expression <- val
  if (is.null(legend_title)) legend_title <- feature
  
  # ---- 4) 真实大小 or 符号大小 ----
  flip_y <- function(y) -y
  
  use_real_tiles <- (shape=="square" && spot_size=="auto" &&
                       ( (!is.null(area_col) && area_col %in% names(plot_df)) || !is.null(bin_um) ) )
  
  if (use_real_tiles) {
    # 尺子向量化（既支持标量，也支持传列名）
    if (is.character(mpp))  mpp_vec <- plot_df[[mpp]]  else mpp_vec <- rep(mpp,  nrow(plot_df))
    if (is.character(s_low))sL_vec  <- plot_df[[s_low]] else sL_vec  <- rep(s_low, nrow(plot_df))
    
    if (!is.null(area_col) && area_col %in% names(plot_df)) {
      # Cellbin：面积 -> 低清像素边长
      if (area_unit=="um2") {
        side <- sqrt(pmax(plot_df[[area_col]], 0)) * (sL_vec / mpp_vec)
      } else {
        side <- sqrt(pmax(plot_df[[area_col]], 0))
      }
    } else {
      # HD：固定 bin
      if (is.null(bin_um)) stop("size='auto'：需要 area_col 或 bin_um。")
      side <- (bin_um / mpp_vec) * sL_vec
    }
    
    cx <- plot_df$imagecol_scaled; cy <- flip_y(plot_df$imagerow_scaled)
    plot_df$.xmin <- cx - 0.5*side; plot_df$.xmax <- cx + 0.5*side
    plot_df$.ymin <- cy - 0.5*side; plot_df$.ymax <- cy + 0.5*side
    
    p <- ggplot(plot_df, aes(xmin=.xmin, xmax=.xmax, ymin=.ymin, ymax=.ymax, fill=Expression)) +
      geom_rect(color = NA, ...) +
      scale_fill_gradientn(colors = colors, name = legend_title, na.value = "lightgrey")
    
  } else {
    # 符号大小（与坐标无关）
    have_scattermore <- requireNamespace("scattermore", quietly = TRUE)
    if (shape=="circle") {
      if (have_scattermore) {
        p <- ggplot(plot_df, aes(x=imagecol_scaled, y=flip_y(imagerow_scaled), color=Expression)) +
          scattermore::geom_scattermore(pointsize = 2.5, pixels = c(2000,2000), ...) +
          scale_color_gradientn(colors = colors, name = legend_title, na.value = "lightgrey")
      } else {
        p <- ggplot(plot_df, aes(x=imagecol_scaled, y=flip_y(imagerow_scaled), color=Expression)) +
          geom_point(size = 2.5, shape = 16, stroke = 0, ...) +
          scale_color_gradientn(colors = colors, name = legend_title, na.value = "lightgrey")
      }
    } else {
      p <- ggplot(plot_df, aes(x=imagecol_scaled, y=flip_y(imagerow_scaled), fill=Expression)) +
        geom_point(shape = 22, size = 2.5, color = ggplot2::alpha("black", 0), stroke = 0.25, ...) +
        scale_fill_gradientn(colors = colors, name = legend_title, na.value = "lightgrey")
    }
  }
  
  p +
    coord_fixed() +
    xlab("") + ylab("") +
    theme_void() +
    theme(
      plot.title      = element_text(hjust=0.5, face="bold"),
      legend.position = "right"
    )
}





hd_rctd_qc_report <- function(
    object,
    cluster_col      = "Seurat_subcluster_Immune_Macrophage",
    spot_class_col   = "RCTD_spot_class",
    first_calc_col   = "RCTD_first_calc",
    second_calc_col  = "RCTD_second_calc",
    margin_col       = "RCTD_margin",
    first_w_col      = "RCTD_first_w",
    second_w_col     = "RCTD_second_w",
    top_n_pairs      = 12,
    margin_highconf  = 0.20,
    min_doublet_cells= 20,
    doublet_rate_flag= 0.05,
    xy_cols          = c("imagecol_scaled","imagerow_scaled"),
    show_spatial     = TRUE,
    output_dir       = NULL  # <--- 新增参数：用于指定输出目录的路径
){
  # --- 函数核心逻辑部分（未做任何修改）---
  pkgs <- c("dplyr","tidyr","ggplot2","scales","patchwork","RColorBrewer","ggrepel")
  sapply(pkgs, function(p) if (!requireNamespace(p, quietly = TRUE)) stop("Please install package: ", p))
  stopifnot(inherits(object, "Seurat"))
  md <- object@meta.data
  need_cols <- c(cluster_col, spot_class_col, first_calc_col, second_calc_col, margin_col)
  miss <- setdiff(need_cols, colnames(md))
  if (length(miss) > 0) stop("Missing columns in meta.data: ", paste(miss, collapse=", "))
  
  md <- md |>
    dplyr::mutate(
      .cluster = as.character(.data[[cluster_col]]),
      .class   = dplyr::case_when(is.na(.data[[spot_class_col]]) ~ "reject", TRUE ~ as.character(.data[[spot_class_col]])),
      .class   = factor(.class, levels = c("singlet","doublet_certain","doublet_uncertain","reject")),
      .margin  = suppressWarnings(as.numeric(.data[[margin_col]])),
      .first   = as.character(.data[[first_calc_col]]),
      .second  = as.character(.data[[second_calc_col]])
    )
  
  tab_long <- md |>
    dplyr::filter(!is.na(.cluster)) |>
    dplyr::count(.cluster, .class, name = "n") |>
    dplyr::group_by(.cluster) |>
    dplyr::mutate(total = sum(n), prop = n/total) |>
    dplyr::ungroup()
  
  tab_wide <- tidyr::pivot_wider(tab_long, names_from = .class, values_from = c(n, prop), values_fill = 0)
  add0 <- function(df, cols) { for (cc in cols) if (!cc %in% colnames(df)) df[[cc]] <- 0; df }
  tab_wide <- add0(tab_wide, c("n_singlet","n_doublet_certain","n_doublet_uncertain","n_reject",
                               "prop_singlet","prop_doublet_certain","prop_doublet_uncertain","prop_reject"))
  
  tab_sum <- tab_wide |>
    dplyr::mutate(
      doublet_total      = n_doublet_certain + n_doublet_uncertain,
      nonreject_total    = n_singlet + doublet_total,
      doublet_rate_nr    = dplyr::if_else(nonreject_total > 0, doublet_total/nonreject_total, NA_real_),
      reject_rate        = dplyr::if_else((n_reject + nonreject_total) > 0, n_reject/(n_reject + nonreject_total), NA_real_),
      flag_has_doublet   = (doublet_total >= min_doublet_cells) | (!is.na(doublet_rate_nr) & doublet_rate_nr >= doublet_rate_flag)
    ) |>
    dplyr::arrange(dplyr::desc(doublet_rate_nr))
  
  dbl_df <- md |>
    dplyr::filter(.class %in% c("doublet_certain","doublet_uncertain")) |>
    dplyr::mutate(
      .pair = ifelse(is.na(.first) | is.na(.second), NA_character_,
                     vapply(seq_len(n()), function(i) paste(sort(c(.first[i], .second[i])), collapse = " + "), character(1)))
    ) |>
    dplyr::filter(!is.na(.pair))
  
  pair_counts <- dbl_df |>
    dplyr::count(.cluster, .pair, name = "Freq")
  
  top_pairs <- pair_counts |>
    dplyr::group_by(.pair) |>
    dplyr::summarise(Freq=sum(Freq), .groups="drop") |>
    dplyr::arrange(dplyr::desc(Freq)) |>
    dplyr::slice_head(n = top_n_pairs) |>
    dplyr::pull(.pair)
  
  pair_mat <- pair_counts |>
    dplyr::filter(.pair %in% top_pairs) |>
    dplyr::mutate(
      .cluster = factor(.cluster, levels = tab_sum$.cluster),
      .pair    = factor(.pair, levels = rev(top_pairs))
    )
  
  md$highconf_singlet <- md$.class == "singlet" & !is.na(md$.margin) & md$.margin >= margin_highconf
  hc_rate <- md |>
    dplyr::group_by(.cluster) |>
    dplyr::summarise(highconf_singlet_rate = mean(highconf_singlet, na.rm=TRUE), .groups="drop")
  tab_sum <- dplyr::left_join(tab_sum, hc_rate, by = ".cluster")
  
  theme_ns <- function(base_size = 11){
    ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = base_size+1, face="bold"),
        axis.text  = ggplot2::element_text(size = base_size, color = "black"),
        plot.title = ggplot2::element_text(size = base_size+4, face="bold", hjust = 0),
        plot.subtitle = ggplot2::element_text(size = base_size+1, color = "grey30"),
        legend.position = "right",
        legend.title = ggplot2::element_text(face="bold")
      )
  }
  pal_class <- c(singlet="#4C78A8", doublet_certain="#E45756", doublet_uncertain="#F58518", reject="#B0B0B0")
  
  p1 <- tab_long |>
    dplyr::mutate(.cluster = factor(.cluster, levels = tab_sum$.cluster)) |>
    ggplot2::ggplot(ggplot2::aes(x = .cluster, y = prop, fill = .class)) +
    ggplot2::geom_col(width = 0.9, color = NA) +
    ggplot2::scale_fill_manual(values = pal_class, name = "RCTD class") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(x = NULL, y = "Within-cluster fraction", title = "RCTD composition") +
    ggplot2::coord_cartesian(ylim=c(0,1)) +
    theme_ns(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  tab_sum2 <- tab_sum |>
    dplyr::mutate(
      clus_fac = factor(.cluster, levels = rev(.cluster)),
      dr = dplyr::if_else(is.na(doublet_rate_nr), 0, doublet_rate_nr)
    )
  xmax <- max(0.05, max(tab_sum2$dr, na.rm = TRUE) + 0.05)
  
  p2 <- tab_sum2 |>
    ggplot2::ggplot(ggplot2::aes(x = dr, y = clus_fac)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = dr, y = clus_fac, yend = clus_fac),
                          color = "#9E9E9E", linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(color = flag_has_doublet), size = 2.8) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#E45756", `FALSE` = "#4C78A8"), guide = "none") +
    ggplot2::geom_text(ggplot2::aes(label = scales::percent(dr, accuracy = 0.1)),
                       nudge_x = 0.02, size = 3, color = "black") +
    ggplot2::geom_vline(xintercept = doublet_rate_flag, linetype = "dashed", color = "#E45756") +
    ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, xmax)) +
    ggplot2::labs(x = "Doublet rate (non-reject)", y = NULL, title = "Doublet rate") +
    theme_ns(base_size = 11)
  
  p3 <- pair_mat |>
    ggplot2::ggplot(ggplot2::aes(x = .cluster, y = .pair, fill = Freq)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradientn(colors = c("#F1EEF6","#BDC9E1","#74A9CF","#2B8CBE","#045A8D"),
                                  name = "count") +
    ggplot2::labs(x = NULL, y = "Top doublet pairs (global)", title = "Doublet pair heatmap") +
    theme_ns(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  p4 <- md |>
    dplyr::filter(!is.na(.margin)) |>
    ggplot2::ggplot(ggplot2::aes(x = .class, y = .margin, fill = .class)) +
    ggplot2::geom_violin(color = NA, alpha = .85) +
    ggplot2::geom_boxplot(width = .15, outlier.size = 0.5, fill = "white", alpha = .6) +
    ggplot2::scale_fill_manual(values = pal_class, guide = "none") +
    ggplot2::geom_hline(yintercept = margin_highconf, linetype = "dashed", color = "#4C78A8") +
    ggplot2::labs(x = NULL, y = "RCTD margin", title = "Margin distribution") +
    theme_ns(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1
      )
    )
  
  has_xy <- all(xy_cols %in% colnames(md))
  if (isTRUE(show_spatial) && has_xy) {
    p5 <- md |>
      dplyr::filter(!is.na(.data[[xy_cols[1]]]), !is.na(.data[[xy_cols[2]]])) |>
      ggplot2::ggplot(ggplot2::aes(x = .data[[xy_cols[1]]], y = .data[[xy_cols[2]]], color = .class)) +
      ggplot2::geom_point(size = 0.6, alpha = 0.8) +
      ggplot2::scale_color_manual(values = pal_class, name = "RCTD class") +
      ggplot2::coord_fixed() +
      ggplot2::scale_y_reverse() +
      ggplot2::labs(x = NULL, y = NULL, title = "Spatial view") +
      theme_ns(base_size = 11) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())
  } else { p5 <- NULL }
  
  panel <- if (is.null(p5)) {
    (p1 + p2) / (p3 + p4) + patchwork::plot_layout(guides = "collect")
  } else {
    p5_cleaned <- p5 + 
      theme(
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
    panel <- ((p1 + p2) / (p3 + p4)) | p5_cleaned + 
      patchwork::plot_layout(widths = c(2, 1.6), guides = "collect")
  }
  
  any_doublet_flag <- any(tab_sum$flag_has_doublet %in% TRUE, na.rm = TRUE)
  
  # --- 新增的保存模块 ---
  if (!is.null(output_dir)) {
    # 1. 创建文件夹（如果不存在）
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # 2. 保存关键的数据表格为 CSV 文件
    write.csv(tab_sum, file.path(output_dir, "summary_table.csv"), row.names = FALSE)
    write.csv(pair_counts, file.path(output_dir, "doublet_pair_counts.csv"), row.names = FALSE)
    
    # 3. 保存组合好的图面板为 PNG 文件
    ggplot2::ggsave(
      filename = file.path(output_dir, "qc_report_panel.png"),
      plot = panel,
      width = 14,
      height = 8,
      dpi = 300,
      bg = "white"
    )
    
    # 4. 在控制台打印消息提示用户
    message("[OK] QC report and plots saved to: ", normalizePath(output_dir))
  }
  # --------------------
  
  out <- list(
    summary_table   = tab_sum,
    composition_long= tab_long,
    pair_table      = pair_counts,
    highconf_rate   = hc_rate,
    plots = list(p_comp = p1, p_rate = p2, p_pairs = p3, p_margin = p4, p_spatial = p5, panel = panel),
    exists_doublet  = any_doublet_flag,
    params = list(top_n_pairs=top_n_pairs, margin_highconf=margin_highconf,
                  min_doublet_cells=min_doublet_cells, doublet_rate_flag=doublet_rate_flag)
  )
  class(out) <- c("hd_rctd_qc_report", class(out))
  return(out)
}

# ==== 子集/标签吸收/SoftGate ====











suppressPackageStartupMessages({ library(dplyr); library(readr) })

# 识别条形码层级的正则（也可手动指定）
.detect_level_by_barcode <- function(barcodes) {
  b <- barcodes[seq_len(min(length(barcodes), 50))]
  if (any(grepl("^s_002um_", b))) return("square_002um")
  if (any(grepl("^s_008um_", b))) return("square_008um")
  if (any(grepl("^s_016um_", b))) return("square_016um")
  if (any(grepl("^cellid_",  b))) return("cell_id")
  stop("无法从条形码模式推断层级，请手动指定。")
}





# 从 Seurat 对象抽取：条形码向量 + 标签向量
.extract_ids_labels <- function(sobj, label_col) {
  stopifnot(label_col %in% colnames(sobj@meta.data))
  ids <- colnames(sobj)
  labs <- as.character(sobj@meta.data[[label_col]])
  data.frame(id = ids, label = labs, stringsAsFactors = FALSE)
}




## ============================================================
## 8) 子集独立导出 + markers（封装，避免重复）
## ============================================================
subset_and_export <- function(obj, keep_expr, name, logfc.threshold = NULL){
  cells <- WhichCells(obj, expression = !!rlang::parse_expr(keep_expr))
  if (length(cells) < 20) { message("[skip] too few cells for ", name); return(invisible(NULL)) }
  counts_sub <- GetAssayData(obj, assay="Spatial", layer="counts")[, cells, drop=FALSE]
  meta_sub   <- obj@meta.data[cells, , drop=FALSE]
  sobj <- CreateSeuratObject(counts = counts_sub, assay = "Spatial", meta.data = meta_sub)
  sobj <- ensure_log_data(sobj, assay = "Spatial")
  saveRDS(sobj, file.path(Config$Bankspath, paste0(name, "_object.rds")))
  # 导出 markers（若有 L2）
  if ("Seurat_L2" %in% colnames(sobj@meta.data)) {
    ExportTopMarkers(
      sobj,
      out_xlsx = file.path(Config$Bankspath, paste0(name, "_L2_cluster_top30_markers.xlsx")),
      group_col = "Seurat_L2",
      assay = "Spatial",
      n_top = 40,
      only.pos = TRUE,
      logfc.threshold = logfc.threshold %||% 0
    )
  }
  invisible(sobj)
}





# 2) 一键补齐 HVG→PCA→邻接图→UMAP（只缩放 HVG，内存友好）
compute_subset_umap <- function(obj, assay="Spatial",
                                n_hvg=3000, dims_use=1:30, k_nn=30,
                                umap_name="umap", verbose=FALSE){
  DefaultAssay(obj) <- assay
  obj <- NormalizeData(obj, assay=assay, normalization.method="LogNormalize", verbose=verbose)
  obj <- FindVariableFeatures(obj, assay=assay, nfeatures=n_hvg,
                              selection.method="vst", verbose=verbose)
  VariableFeatures(obj, assay=assay) <- .drop_tech_genes(VariableFeatures(obj, assay=assay))
  obj <- ScaleData(obj, assay=assay, features=VariableFeatures(obj, assay=assay),
                   block.size=2000, verbose=verbose)
  obj <- RunPCA(obj, assay=assay, features=VariableFeatures(obj, assay=assay),
                npcs=max(dims_use), verbose=verbose)
  obj <- FindNeighbors(obj, reduction="pca", dims=dims_use, k.param=k_nn, verbose=verbose)
  obj <- RunUMAP(obj, reduction="pca", dims=dims_use,
                 reduction.name=umap_name, verbose=verbose)
  obj
}





transfer_labels_via_2um <- function(mapping,
                                    source_obj, target_obj,
                                    source_label_col,
                                    source_level = NULL,  # 为空则自动根据条码推断
                                    target_level = NULL,  # 为空则自动根据条码推断
                                    filter_in_cell = FALSE,
                                    filter_in_nucleus = FALSE,
                                    min_prop = 0.5,
                                    unknown_label = "Unassigned",
                                    out_col = NULL,
                                    use_datatable = TRUE) {
  # 1) 推断层级
  if (is.null(source_level)) source_level <- .detect_level_by_barcode(colnames(source_obj))
  if (is.null(target_level)) target_level <- .detect_level_by_barcode(colnames(target_obj))
  stopifnot(all(c("square_002um", source_level, target_level) %in% colnames(mapping)))
  
  # 2) 只保留会用到的列，尽量早过滤以加速
  keep_cols <- c("square_002um", source_level, target_level, "in_cell", "in_nucleus")
  map_use <- mapping[, intersect(keep_cols, colnames(mapping))]
  rm(mapping); invisible(gc())
  
  # 3) 过滤 2µm（可选：只计 in_cell / in_nucleus）
  if ("in_cell" %in% colnames(map_use) && isTRUE(filter_in_cell)) {
    map_use <- dplyr::filter(map_use, .data$in_cell)
  }
  if ("in_nucleus" %in% colnames(map_use) && isTRUE(filter_in_nucleus)) {
    map_use <- dplyr::filter(map_use, .data$in_nucleus)
  }
  
  # 4) 只保留来源/目标在“当前对象”内的 2µm 行
  src_ids <- colnames(source_obj)
  tgt_ids <- colnames(target_obj)
  map_use <- dplyr::filter(map_use, .data[[source_level]] %in% src_ids,
                           .data[[target_level]] %in% tgt_ids)
  
  # 5) 把“来源层级标签”贴到 2µm 上
  src_tbl <- .extract_ids_labels(source_obj, source_label_col)
  names(src_tbl) <- c(source_level, "src_label")
  map_lab <- dplyr::inner_join(map_use, src_tbl, by = source_level)
  rm(map_use, src_tbl); invisible(gc())
  
  # 6) 对 2µm -> 目标层级 聚合：统计每个目标单元内各标签所占 2µm 比例
  #    用 data.table 会快很多
  if (use_datatable && requireNamespace("data.table", quietly = TRUE)) {
    DT <- data.table::as.data.table(map_lab)[, .N, by = c(target_level, "src_label")]
    data.table::setnames(DT, "N", "n2um")
    DT_all <- DT[, .(ntot = sum(n2um)), by = target_level]
    DT <- DT[DT_all, on = target_level]
    DT[, prop := n2um / ntot]
    # 每个目标取 prop 最大的标签
    setord <- c(target_level, "-prop")
    DT <- DT[order(get(target_level), -prop)]
    top <- DT[, .SD[1], by = target_level]
    winner <- data.frame(id = top[[target_level]],
                         label = top$src_label,
                         prop  = top$prop, stringsAsFactors = FALSE)
  } else {
    winner <- map_lab %>%
      dplyr::count(.data[[target_level]], src_label, name = "n2um") %>%
      dplyr::group_by(.data[[target_level]]) %>%
      dplyr::mutate(prop = n2um / sum(n2um)) %>%
      dplyr::arrange(.data[[target_level]], dplyr::desc(prop)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(id = .data[[target_level]],
                       label = src_label,
                       prop = prop)
  }
  
  # 7) 低于阈值的目标标记为 unknown
  winner$label_final <- ifelse(winner$prop >= min_prop, winner$label, unknown_label)
  
  # 8) 写回 target meta
  if (is.null(out_col)) {
    out_col <- paste0(source_label_col, "_via2um_from_", sub("^square_", "", source_level))
  }
  target_obj@meta.data[[out_col]]          <- winner$label_final[match(colnames(target_obj), winner$id)]
  target_obj@meta.data[[paste0(out_col,"_maxprop")]] <- winner$prop[match(colnames(target_obj), winner$id)]
  
  # 返回结果与可选的长表（用于审计/画图）
  list(
    target = target_obj,
    winner = winner  # id, label(赢家), prop(最大占比)
  )
}





export_labels_csv <- function(obj, label_col, out_csv) {
  stopifnot(label_col %in% colnames(obj@meta.data))
  df <- data.frame(barcode = colnames(obj),
                   label   = obj@meta.data[[label_col]],
                   check.names = FALSE)
  readr::write_csv(df, out_csv)
  message("Saved: ", out_csv)
}





# ==== Hotspot 与 LR 分析 ====
#  A) Numeric renaming
# =========================

# Formats a vector of group labels into numeric-named labels (e.g., "C01_Basal").
# It orders by the first integer in each label (if found), else alphabetical.
autonumber_groups <- function(
    x,                            # character or factor vector of group labels (length = cells)
    prefix = "C",                 # leading prefix
    width  = 2,                   # zero-padding width
    fmt    = "{prefix}{num}_{label}",  # final label format
    order_by = c("numeric_then_alpha","alpha"),  # primary ordering rule
    keep_original_levels = TRUE
){
  x_chr <- as.character(x)
  order_by <- match.arg(order_by)
  
  # extract first integer per unique label
  uniq <- unique(x_chr)
  m <- regexpr("-?[0-9]+", uniq)
  key <- ifelse(m > 0, suppressWarnings(as.numeric(regmatches(uniq, m))), NA_real_)
  
  if (order_by == "numeric_then_alpha") {
    ord <- order(is.na(key), key, tolower(uniq))
  } else {
    ord <- order(tolower(uniq))
  }
  
  uniq_ord <- uniq[ord]
  # assign running numbers 1..k by order
  num <- seq_along(uniq_ord)
  
  # formatted numbers
  num_fmt <- sprintf(paste0("%0", width, "d"), num)
  
  # build new labels
  format_one <- function(prefix, num_str, label) {
    sub("\\{prefix\\}", prefix,
        sub("\\{num\\}",    num_str,
            sub("\\{label\\}", label, fmt, fixed = TRUE),
            fixed = TRUE),
        fixed = TRUE)
  }
  new_labels_ord <- mapply(format_one, prefix, num_fmt, uniq_ord, USE.NAMES = FALSE)
  
  # map back to all cells
  map_df <- data.frame(
    old = uniq_ord,
    number = num,
    number_str = num_fmt,
    new = new_labels_ord,
    stringsAsFactors = FALSE
  )
  map <- setNames(map_df$new, map_df$old)
  x_new <- unname(map[x_chr])
  
  # levels for factors
  levels_new <- if (keep_original_levels) new_labels_ord else sort(unique(x_new))
  
  list(
    x_new = x_new,
    map_df = map_df,         # data.frame: old → new (+ numbers)
    rename_map = map,        # named character: old -> new
    levels_new = levels_new  # recommended factor level order
  )
}





# =========================
#  B) High-contrast palette
# =========================

# Returns a *named* color vector aligned to the provided labels (vector or levels).
# "glasbey" (from pals) is default: unlimited, crisp, color-blind friendly.
palette_for_groups <- function(
    x_or_levels,
    scheme = c("glasbey","polychrome","hcl","viridis"),
    tone   = c("vivid","pastel"),   # mild post-processing choice
    seed   = 42L
){
  scheme <- match.arg(scheme)
  tone   <- match.arg(tone)
  
  # accept a vector of labels (factor/character) or a character vector of levels
  if (is.factor(x_or_levels)) {
    lv <- levels(x_or_levels)
    if (is.null(lv) || length(lv) == 0L) lv <- unique(as.character(x_or_levels))
  } else {
    lv <- unique(as.character(x_or_levels))
  }
  lv <- as.character(lv)
  n  <- length(lv)
  if (n == 0L) stop("No group levels found: please pass a non-empty vector.")
  
  # generator
  gen_cols <- function(n){
    if (scheme == "glasbey" && requireNamespace("pals", quietly = TRUE)) {
      # unlimited qualitative; excellent separability
      pals::glasbey(n)
    } else if (scheme == "polychrome" && requireNamespace("Polychrome", quietly = TRUE)) {
      set.seed(seed)
      Polychrome::createPalette(n, seedcolors = c("#ff0000","#00ff00","#0000ff"))
    } else if (scheme == "hcl" && requireNamespace("colorspace", quietly = TRUE)) {
      colorspace::qualitative_hcl(n, palette = "Dark 3")
    } else if (requireNamespace("viridisLite", quietly = TRUE)) {
      viridisLite::viridis(n) # sequential; used as robust fallback
    } else {
      grDevices::hcl(h = seq(15, 375, length.out = n + 1L)[1L:n], c = 60, l = 65)
    }
  }
  
  cols <- gen_cols(n)
  
  # optional toning (gentle, keeps contrasts)
  if (tone == "pastel") {
    if (requireNamespace("colorspace", quietly = TRUE)) {
      cols <- colorspace::lighten(cols, amount = 0.18)
      cols <- colorspace::desaturate(cols, amount = 0.10)
    }
  } else {
    # vivid: slightly deepen to print nicely
    if (requireNamespace("colorspace", quietly = TRUE)) {
      cols <- colorspace::darken(cols, amount = 0.05)
    }
  }
  
  setNames(cols, lv)
}





`%||%` <- function(a,b) if (!is.null(a)) a else b
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





`%||%` <- function(a,b) if (!is.null(a)) a else b

PlotSpatialDistribution_bar <- function(
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
    title = NULL,
    legend_title = NULL,
    legend_spacing_y = 0.3,
    legend_key_height = 0.6,
    legend_key_width  = 0.8,
    add_scale_bar        = FALSE,
    scale_bar_length_mm  = 1,
    scale_bar_n_segments = 0,
    scale_bar_color      = "navyblue",
    scale_bar_thickness_px = 4,
    scale_bar_position   = c("bottomright","bottomleft","topright","topleft"),
    scale_bar_margin_frac = c(0.04, 0.05),
    scale_bar_text       = NULL,
    scale_bar_text_size  = 3.2,
    scale_bar_text_color = "black",
    scale_bar_text_face  = "bold",
    scale_bar_text_offset_px = 10,
    scale_bar_gap_px     = 1,
    ...
){
  spot_shape <- match.arg(spot_shape)
  size_mode  <- match.arg(size_mode)
  scale_bar_position <- match.arg(scale_bar_position)
  
  # -----------------------------
  # 0) Helper: extract scale factors (minimal + robust)
  # -----------------------------
  get_scale_factors <- function(obj) {
    mpp <- NULL
    s_low <- NULL
    s_hi  <- NULL
    
    # (1) Try Seurat image slot first
    if (inherits(obj, "Seurat")) {
      mpp  <- tryCatch(Images(obj)[[1]]@scale.factors$microns_per_pixel, error=function(e) NULL)
      s_low <- tryCatch(Images(obj)[[1]]@scale.factors$tissue_lowres_scalef, error=function(e) NULL)
      s_hi  <- tryCatch(Images(obj)[[1]]@scale.factors$tissue_hires_scalef, error=function(e) NULL)
    }
    
    # (2) Fallback: your object stores them in obj@misc$scales / obj@misc
    if (is.null(mpp))  mpp  <- obj@misc$scales$microns_per_pixel %||% obj@misc$microns_per_pixel %||% obj@misc$microns_per_pixel
    if (is.null(s_low)) s_low <- obj@misc$scales$tissue_lowres_scalef %||% obj@misc$scales$tissue_lowres_scalef
    if (is.null(s_hi))  s_hi  <- obj@misc$scales$tissue_hires_scalef %||% obj@misc$scales$tissue_hires_scalef
    
    # (3) Final fallback
    if (is.null(mpp))  mpp  <- 1
    if (is.null(s_low)) s_low <- 1
    if (is.null(s_hi))  s_hi  <- s_low
    
    list(mpp = as.numeric(mpp), s_low = as.numeric(s_low), s_hi = as.numeric(s_hi))
  }
  
  # Decide which scale factor matches coord_cols
  infer_scale_for_coords <- function(coord_cols, s_low, s_hi) {
    cc <- paste(coord_cols, collapse="|")
    # Default: your coord_cols are lowres scaled ("*_scaled")
    if (grepl("hires", cc, ignore.case = TRUE)) return(s_hi)
    if (grepl("scaled", cc, ignore.case = TRUE)) return(s_low)
    # imagecol / imagerow (fullres pixel space)
    return(1)
  }
  
  # -----------------------------
  # 1) Data Prep
  # -----------------------------
  if (inherits(object, "Seurat")) {
    plot_df <- object@meta.data
    
    sc <- get_scale_factors(object)
    mpp <- microns_per_pixel %||% sc$mpp
    
    # IMPORTANT: scale_lowres here is treated as "scale factor for coord_cols".
    # If not provided, auto-infer it from coord_cols (lowres/hires/fullres).
    s_auto <- infer_scale_for_coords(coord_cols, sc$s_low, sc$s_hi)
    s_use  <- scale_lowres %||% s_auto
    
    if (is.null(area_col) && "area_um2" %in% names(plot_df)) area_col <- "area_um2"
    if (is.null(bin_um)) bin_um <- object@misc$bin_um %||% 16
  } else {
    plot_df <- object
    mpp  <- microns_per_pixel %||% 1
    s_use <- scale_lowres %||% 1
    if (is.null(bin_um)) bin_um <- 16
  }
  
  if (!group_by %in% names(plot_df)) stop("Grouping column not found.")
  if (!all(coord_cols %in% names(plot_df))) stop("coord_cols not found in data.")
  
  plot_df[[group_by]] <- as.character(plot_df[[group_by]])
  plot_df[[group_by]][is.na(plot_df[[group_by]])] <- "Unassigned"
  plot_df[[group_by]] <- as.factor(plot_df[[group_by]])
  levs <- levels(plot_df[[group_by]])
  
  # -----------------------------
  # 2) Color logic (keep your fixed behavior)
  # -----------------------------
  legend_name <- if (is.null(legend_title)) group_by else legend_title
  
  if (is.character(palette) && length(palette) == 1 && !is.null(names(palette))) {
    plot_colors <- palette
  } else if (is.character(palette) && length(palette) == 1) {
    if (requireNamespace("RColorBrewer", quietly = TRUE) &&
        palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      cols <- RColorBrewer::brewer.pal(max(3, min(length(levs), 8)), palette)
      plot_colors <- rep(cols, length.out = length(levs))
      names(plot_colors) <- levs
    } else {
      plot_colors <- rep(palette, length.out = length(levs))
      names(plot_colors) <- levs
    }
  } else if (is.character(palette)) {
    if (!is.null(names(palette))) {
      # Critical: keep names if provided
      plot_colors <- palette
    } else {
      plot_colors <- rep(palette, length.out = length(levs))
      names(plot_colors) <- levs
    }
  } else {
    plot_colors <- scales::hue_pal()(length(levs))
    names(plot_colors) <- levs
  }
  
  # Reorder palette strictly by levels (prevents any mismatch)
  if (!is.null(names(plot_colors))) {
    plot_colors <- plot_colors[levs]
  } else {
    names(plot_colors) <- levs
  }
  plot_colors[is.na(plot_colors)] <- "grey90"
  if ("Unassigned" %in% levs) plot_colors["Unassigned"] <- "#FFFFFF"
  
  # -----------------------------
  # 3) Geometry
  # -----------------------------
  library(ggplot2)
  flip_y <- function(y) -y
  
  use_real_tiles <- (spot_shape == "square" && size_mode == "auto")
  
  if (use_real_tiles) {
    mpp_val  <- if (is.numeric(mpp)) mpp else 1
    s_val    <- if (is.numeric(s_use)) s_use else 1
    
    if (!is.null(area_col) && area_col %in% names(plot_df)) {
      side <- sqrt(plot_df[[area_col]]) * (s_val / mpp_val)
    } else {
      side <- (bin_um / mpp_val) * s_val
    }
    
    cx <- plot_df[[coord_cols[1]]]
    cy <- flip_y(plot_df[[coord_cols[2]]])
    plot_df$.xmin <- cx - 0.5*side
    plot_df$.xmax <- cx + 0.5*side
    plot_df$.ymin <- cy - 0.5*side
    plot_df$.ymax <- cy + 0.5*side
    
    p <- ggplot(plot_df, aes(xmin=.xmin, xmax=.xmax, ymin=.ymin, ymax=.ymax, fill=.data[[group_by]])) +
      geom_rect(color=NA, ...) +
      scale_fill_manual(values = plot_colors, name = legend_name, breaks = levs, drop = FALSE, na.value="grey90")
  } else {
    p <- ggplot(plot_df, aes(x=.data[[coord_cols[1]]], y=flip_y(.data[[coord_cols[2]]]), color=.data[[group_by]])) +
      geom_point(size=ptsize, shape=if(spot_shape=="circle") 16 else 15, ...) +
      scale_color_manual(values = plot_colors, name = legend_name, breaks = levs, drop = FALSE, na.value="grey90")
  }
  
  # -----------------------------
  # 4) Theme
  # -----------------------------
  p <- p + coord_fixed() + theme_void() +
    labs(title = title %||% paste("Spatial:", group_by)) +
    theme(
      plot.title = element_text(hjust=0.5, face="bold", size=14),
      legend.position = if(show_legend) "right" else "none",
      legend.spacing.y = unit(legend_spacing_y, "lines"),
      legend.key.height = unit(legend_key_height, "lines"),
      legend.key.width  = unit(legend_key_width,  "lines")
    )
  
  # -----------------------------
  # 5) Scale bar (minimal fix + no per-row mapping)
  # -----------------------------
  if (isTRUE(add_scale_bar)) {
    xr <- range(plot_df[[coord_cols[1]]], na.rm=TRUE)
    yr <- range(flip_y(plot_df[[coord_cols[2]]]), na.rm=TRUE)
    
    mpp_val <- if (is.numeric(mpp)) mpp else 1
    s_val   <- if (is.numeric(s_use)) s_use else 1
    
    # Correct conversion for the coordinate system used by coord_cols:
    # - fullres coords: s_val = 1
    # - lowres coords:  s_val = tissue_lowres_scalef
    # - hires coords:   s_val = tissue_hires_scalef
    # 1 mm = (1000 um) * (s_val / mpp_val) pixels in that coordinate space
    bar_px <- scale_bar_length_mm * 1000 * (s_val / mpp_val)
    
    mx <- scale_bar_margin_frac[1] * diff(xr)
    my <- scale_bar_margin_frac[2] * diff(yr)
    
    x0 <- if (scale_bar_position %in% c("bottomright","topright")) xr[2] - mx - bar_px else xr[1] + mx
    y0 <- if (scale_bar_position %in% c("bottomright","bottomleft")) yr[1] + my else yr[2] - my - scale_bar_thickness_px
    
    label_txt <- scale_bar_text %||% paste0(scale_bar_length_mm, " mm")
    
    p <- p +
      annotate("rect",
               xmin = x0, xmax = x0 + bar_px,
               ymin = y0, ymax = y0 + scale_bar_thickness_px,
               fill = scale_bar_color, color = NA) +
      annotate("text",
               x = x0 + bar_px/2,
               y = y0 + scale_bar_thickness_px + scale_bar_text_offset_px,
               label = label_txt,
               size = scale_bar_text_size,
               colour = scale_bar_text_color,
               fontface = scale_bar_text_face)
  }
  
  return(p)
}





`%||%` <- function(a,b) if (!is.null(a)) a else b

# Minimal-change fix:
# - Auto-pick the proper scalefactor based on coord_cols:
#     * hires coords  -> use tissue_hires_scalef (s_hires)
#     * lowres coords -> use tissue_lowres_scalef (s_low)
# - Use s_use consistently for:
#     * real tile sizing (area_um2 / bin_um)
#     * scale bar px-per-um conversion
# - HE overlay logic remains unchanged; it already aligns to the flipped-y frame.
PlotSpatialDistribution_bar_he <- function(
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
    title = NULL,
    legend_title = NULL,
    legend_spacing_y = 0.3,
    legend_key_height = 0.6,
    legend_key_width  = 0.8,
    add_scale_bar        = FALSE,
    scale_bar_length_mm  = 1,
    scale_bar_n_segments = 4,
    scale_bar_color      = "navyblue",
    scale_bar_thickness_px = 4,
    scale_bar_position   = c("bottomright","bottomleft","topright","topleft"),
    scale_bar_margin_frac = c(0.04, 0.05),
    scale_bar_text       = NULL,
    scale_bar_text_size  = 3.2,
    scale_bar_text_color = "black",
    scale_bar_text_face  = "bold",
    scale_bar_text_offset_px = 10,
    scale_bar_gap_px     = 1,
    # ---- NEW: physical thickness override (μm) ----
    scale_bar_thickness_um = NULL,
    
    # ---- HE background ----
    show_he        = FALSE,
    he_outs_path   = NULL,
    he_resolution  = c("auto","hires","lowres"),
    # ---- NEW: clip view to data extent when HE is shown ----
    clip_to_data   = FALSE,
    ...
){
  spot_shape <- match.arg(spot_shape)
  size_mode  <- match.arg(size_mode)
  scale_bar_position <- match.arg(scale_bar_position)
  he_resolution <- match.arg(he_resolution)
  
  # ---------- 1) Data & spatial scales ----------
  if (inherits(object, "Seurat")) {
    plot_df <- object@meta.data
    
    mpp <- microns_per_pixel %||%
      object@misc$microns_per_pixel %||%
      object@misc$spatial_scales$microns_per_pixel %||%
      object@misc$scales$microns_per_pixel %||%
      tryCatch(Images(object)[[1]]@scale.factors$microns_per_pixel, error=function(e) NULL)
    
    s_low <- scale_lowres %||%
      object@misc$spatial_scales$tissue_lowres_scalef %||%
      object@misc$scales$tissue_lowres_scalef %||%
      tryCatch(Images(object)[[1]]@scale.factors$tissue_lowres_scalef, error=function(e) NULL)
    
    # also try to fetch hires scalefactor
    s_hires <- object@misc$spatial_scales$tissue_hires_scalef %||%
      object@misc$scales$tissue_hires_scalef %||%
      tryCatch(Images(object)[[1]]@scale.factors$tissue_hires_scalef, error=function(e) NULL)
    
    # fallback estimate for s_hires if missing but both coord systems exist
    if (is.null(s_hires) && !is.null(s_low) &&
        all(c("imagecol_scaled","imagecol_hires") %in% colnames(plot_df))) {
      ratio <- plot_df$imagecol_scaled / plot_df$imagecol_hires
      ratio <- stats::median(ratio[is.finite(ratio)], na.rm = TRUE)
      if (is.finite(ratio) && ratio > 0) s_hires <- s_low / ratio
    }
    
    if (is.null(mpp)) stop("microns_per_pixel not found; pass microns_per_pixel = numeric.")
    if (is.null(s_low)) stop("tissue_lowres_scalef not found; pass scale_lowres = numeric.")
    
    # choose scalefactor by coord system (KEY)
    is_hires_coord <- any(grepl("hires", coord_cols, ignore.case = TRUE))
    s_use <- if (isTRUE(is_hires_coord)) {
      if (is.null(s_hires)) stop("tissue_hires_scalef missing while using hires coords.")
      s_hires
    } else {
      s_low
    }
    
    if (is.null(area_col) && "area_um2" %in% names(plot_df)) area_col <- "area_um2"
    if (is.null(bin_um)) {
      bin_um <- object@misc$bin_um %||%
        object@misc$spatial_scales$bin_size_um %||% NULL
    }
  } else if (is.data.frame(object)) {
    plot_df <- object
    mpp  <- microns_per_pixel
    s_low<- scale_lowres
    if (is.null(mpp) || is.null(s_low))
      stop("When passing a data.frame, provide both microns_per_pixel and scale_lowres.")
    s_use <- s_low
    if (is.null(area_col) && "area_um2" %in% names(plot_df)) area_col <- "area_um2"
  } else {
    stop("'object' must be a Seurat object or a data.frame.")
  }
  
  if (!all(coord_cols %in% names(plot_df)))
    stop("Missing coordinate columns: ", paste(coord_cols, collapse=", "))
  x_col <- coord_cols[1]; y_col <- coord_cols[2]
  
  if (!group_by %in% names(plot_df)) stop("Grouping column not found: ", group_by)
  plot_df[[group_by]] <- as.factor(plot_df[[group_by]])
  levs <- levels(plot_df[[group_by]])
  
  legend_name <- if (is.null(legend_title)) group_by else legend_title
  if (is.character(palette) && length(palette)==1) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) suppressWarnings(install.packages("RColorBrewer"))
    pal_name <- palette
    if (pal_name %in% rownames(RColorBrewer::brewer.pal.info)) {
      maxc <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
      cols <- RColorBrewer::brewer.pal(min(maxc, max(length(levs), 3)), pal_name)
      plot_colors <- rep(cols, length.out = length(levs)); names(plot_colors) <- levs
    } else {
      stop("Unknown palette: ", palette)
    }
  } else if (is.character(palette)) {
    plot_colors <- palette; names(plot_colors) <- levs
  } else if (is.function(palette)) {
    tmp <- palette(); plot_colors <- rep(tmp, length.out=length(levs)); names(plot_colors) <- levs
  } else stop("Invalid 'palette'.")
  

  he_layer <- NULL
  if (isTRUE(show_he)) {
    if (!requireNamespace("png", quietly = TRUE)) suppressWarnings(install.packages("png"))
    he_res <- he_resolution
    if (he_res == "auto") {
      he_res <- if (any(grepl("hires", coord_cols, ignore.case = TRUE))) "hires" else "lowres"
    }
    outs_path <- he_outs_path %||% object@misc$OutsPath %||% object@misc$outs_path %||% object@misc$outs
    if (is.null(outs_path)) stop("show_he=TRUE but outs_path is unknown; pass he_outs_path=...")
    img_file <- file.path(outs_path, "spatial",
                          if (he_res=="lowres") "tissue_lowres_image.png" else "tissue_hires_image.png")
    if (!file.exists(img_file)) stop("HE image not found: ", img_file)
    img <- png::readPNG(img_file); Wimg <- ncol(img); Himg <- nrow(img)
    
    # main plot uses y = -y_raw; so HE goes to [ymin=-Himg, ymax=0]
    he_layer <- ggplot2::annotation_raster(img, xmin = 0, xmax = Wimg, ymin = -Himg, ymax = 0, interpolate = TRUE)
  }
  
  library(ggplot2)
  flip_y <- function(y) -y
  
  use_real_tiles <- (spot_shape=="square" && size_mode=="auto" &&
                       ( (!is.null(area_col) && area_col %in% names(plot_df)) || !is.null(bin_um) ))
  
  if (use_real_tiles) {
    mpp_vec <- if (is.numeric(mpp)) rep(mpp, nrow(plot_df)) else plot_df[[mpp]]
    if (!is.null(area_col) && area_col %in% names(plot_df)) {
      side <- sqrt(pmax(plot_df[[area_col]], 0)) * (s_use / mpp_vec)  # use s_use
    } else {
      if (is.null(bin_um)) stop("size_mode='auto' needs area_col or bin_um.")
      side <- (bin_um / mpp_vec) * s_use                                 # use s_use
    }
    cx <- plot_df[[x_col]]; cy <- flip_y(plot_df[[y_col]])
    plot_df$.xmin <- cx - 0.5*side; plot_df$.xmax <- cx + 0.5*side
    plot_df$.ymin <- cy - 0.5*side; plot_df$.ymax <- cy + 0.5*side
    
    p <- ggplot(plot_df, aes(xmin=.xmin, xmax=.xmax, ymin=.ymin, ymax=.ymax,
                             fill=.data[[group_by]])) +
      { if (!is.null(he_layer)) he_layer else NULL } +
      geom_rect(color=NA, ...) +
      scale_fill_manual(values = plot_colors, name = legend_name)
  } else {
    if (requireNamespace("scattermore", quietly = TRUE) && spot_shape=="circle") {
      p <- ggplot(plot_df, aes(x=.data[[x_col]], y=flip_y(.data[[y_col]]),
                               color=.data[[group_by]])) +
        { if (!is.null(he_layer)) he_layer else NULL } +
        scattermore::geom_scattermore(pointsize = ptsize, pixels = c(2000,2000), ...)
    } else {
      shape_val <- if (spot_shape=="circle") 16 else 15
      p <- ggplot(plot_df, aes(x=.data[[x_col]], y=flip_y(.data[[y_col]]),
                               color=.data[[group_by]])) +
        { if (!is.null(he_layer)) he_layer else NULL } +
        geom_point(size=ptsize, shape=shape_val, ...)
    }
    p <- p + scale_color_manual(values = plot_colors, name = legend_name)
  }

  plot_title <- if (is.null(title)) paste("Spatial Distribution of", group_by) else title
  
  # compute data-window for optional clipping (x in raw coords; y in flipped coords)
  xr  <- range(plot_df[[x_col]], na.rm = TRUE)
  yrp <- range(flip_y(plot_df[[y_col]]), na.rm = TRUE)
  
  # if clip_to_data is TRUE and HE is shown, crop the canvas to data extent
  if (isTRUE(clip_to_data) && isTRUE(show_he)) {
    p <- p + coord_fixed(xlim = xr, ylim = yrp, expand = FALSE)
  } else {
    p <- p + coord_fixed()
  }
  
  p <- p +
    theme_void() +
    labs(title = plot_title) +
    theme(
      plot.title        = ggplot2::element_text(hjust=0.5, face="bold", size=14),
      legend.text       = ggplot2::element_text(size=8),
      legend.key.height = grid::unit(legend_key_height, "lines"),
      legend.key.width  = grid::unit(legend_key_width,  "lines"),
      legend.spacing.y  = grid::unit(legend_spacing_y,  "lines")
    ) +
    { if (!show_legend) ggplot2::guides(color="none", fill="none") else NULL }
  

  if (isTRUE(add_scale_bar)) {
    xs <- plot_df[[x_col]]
    ys <- flip_y(plot_df[[y_col]])
    xr <- range(xs, na.rm = TRUE)
    yr <- range(ys, na.rm = TRUE)
    dx <- diff(xr); dy <- diff(yr)
    
    mpp_use  <- if (is.numeric(mpp)) mpp else stats::median(plot_df[[mpp]], na.rm = TRUE)
    px_per_um <- s_use / mpp_use
    
    bar_um <- scale_bar_length_mm * 1000
    bar_px <- as.integer(max(10, round(bar_um * px_per_um)))
    
    # thickness in px (override by μm if provided)
    th_px <- as.integer(round(scale_bar_thickness_px))
    if (!is.null(scale_bar_thickness_um)) {
      th_px <- max(1L, as.integer(round(scale_bar_thickness_um * px_per_um)))
    }
    
    mx <- if (length(scale_bar_margin_frac) >= 1) scale_bar_margin_frac[1] else 0.04
    my <- if (length(scale_bar_margin_frac) >= 2) scale_bar_margin_frac[2] else mx
    
    x0 <- if (scale_bar_position %in% c("bottomright","topright")) xr[2] - mx*dx - bar_px else xr[1] + mx*dx
    y0 <- if (scale_bar_position %in% c("bottomright","bottomleft")) yr[1] + my*dy else yr[2] - my*dy - th_px
    
    x0 <- as.integer(round(x0)); y0 <- as.integer(round(y0))
    y1 <- y0 + as.integer(round(th_px))
    x1 <- x0 + bar_px
    
    nseg  <- max(1L, as.integer(round(scale_bar_n_segments)))
    gappx <- max(0L, as.integer(round(scale_bar_gap_px)))
    if (nseg == 1L) gappx <- 0L
    
    total_gap_px  <- gappx * (nseg - 1L)
    fill_px_total <- max(bar_px - total_gap_px, nseg)
    seg_fill_base <- floor(fill_px_total / nseg)
    rem <- fill_px_total - seg_fill_base * nseg
    seg_fill <- rep(seg_fill_base, nseg)
    if (rem > 0L) seg_fill[seq_len(rem)] <- seg_fill[seq_len(rem)] + 1L
    
    left <- x0
    seg_list <- vector("list", nseg)
    for (i in seq_len(nseg)) {
      right <- left + seg_fill[i]
      seg_list[[i]] <- data.frame(xmin=left, xmax=right, ymin=y0, ymax=y1)
      left <- right + if (i < nseg) gappx else 0L
    }
    sb_df <- do.call(rbind, seg_list)
    
    p <- p + ggplot2::geom_rect(
      data = sb_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = scale_bar_color, color = NA
    )
    
    if (!is.null(scale_bar_text) && nzchar(scale_bar_text)) {
      text_y <- if (scale_bar_position %in% c("bottomright","bottomleft")) y1 + scale_bar_text_offset_px else y0 - scale_bar_text_offset_px
      p <- p + ggplot2::geom_text(
        data = data.frame(x = x0 + bar_px/2, y = text_y, lab = scale_bar_text),
        ggplot2::aes(x = x, y = y, label = lab),
        inherit.aes = FALSE,
        size = scale_bar_text_size, color = scale_bar_text_color,
        fontface = scale_bar_text_face,
        vjust = if (scale_bar_position %in% c("bottomright","bottomleft")) 0 else 1
      )
    }
  }
  
  return(p)
}





# PlotSpatialDistribution_axes + (optional) scale bar + (optional) HE background.
# Core plotting/axes/ROI logic is unchanged; only minimal additive features.
`%||%` <- function(a,b) if (!is.null(a)) a else b

# Small-fix version:
# - Auto-pick the proper spatial scalefactor based on coord_cols:
#     * hires coords  -> use tissue_hires_scalef
#     * lowres coords -> use tissue_lowres_scalef
# - Use the chosen scalefactor (s_use) consistently for:
#     * real tile sizing (area_um2 / bin_um)
#     * ROI um<->px conversion
#     * secondary axes (px -> um/mm)
#     * scale bar px conversion
# Everything else stays the same. Comments are English only.
PlotSpatialDistribution_axes_bar <- function(
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
    title = NULL,
    legend_title = NULL,
    legend_spacing_y = 0.3,
    legend_key_height = 0.6,
    legend_key_width  = 0.8,
    with_axes_grid      = FALSE,
    axes_secondary_unit = c("mm","um","none"),
    axes_breaks_n       = 8,
    grid_minor          = TRUE,
    axes_text_size      = 10,
    grid_major_col      = "grey90",
    grid_minor_col      = "grey95",
    roi_rect_px         = NULL,      # c(xmin, xmax, ymin, ymax) in lowres/hires px (match coord_cols)
    roi_rect_um         = NULL,      # c(xmin_um, xmax_um, ymin_um, ymax_um)
    roi_center_px       = NULL,      # c(x_px, y_px)
    roi_center_um       = NULL,      # c(x_um, y_um)
    roi_radius_um       = NULL,      # numeric, for circle/square
    roi_shape           = c("circle","square"),
    roi_line_color      = "black",
    roi_linewidth       = 0.8,
    roi_fill            = NA,        # e.g. "#00000010" for translucent
    return_roi_barcodes = TRUE,
    add_scale_bar        = FALSE,
    scale_bar_length_mm  = 1,
    scale_bar_n_segments = 0,
    scale_bar_color      = "navyblue",
    scale_bar_thickness_px = 4,
    scale_bar_position   = c("bottomright","bottomleft","topright","topleft"),
    scale_bar_margin_frac = c(0.04, 0.05),
    scale_bar_text       = NULL,    # NULL = no label
    scale_bar_text_size  = 3.2,
    scale_bar_text_color = "black",
    scale_bar_text_face  = "bold",
    scale_bar_text_offset_px = 10,
    scale_bar_gap_px     = 1,
    scale_bar_thickness_um = NULL,  # if set, overrides thickness_px using physical µm
    show_he        = FALSE,               # overlay HE png
    he_outs_path   = NULL,                # 10x outs path
    he_resolution  = c("auto","lowres","hires"),
    clip_to_data   = FALSE,
    ...
){
  spot_shape <- match.arg(spot_shape)
  size_mode  <- match.arg(size_mode)
  axes_secondary_unit <- match.arg(axes_secondary_unit)
  roi_shape  <- match.arg(roi_shape)
  scale_bar_position <- match.arg(scale_bar_position)
  he_resolution <- match.arg(he_resolution)
  
  # ---------- 1) Data & spatial scales ----------
  if (inherits(object, "Seurat")) {
    plot_df <- object@meta.data
    
    mpp <- microns_per_pixel %||%
      object@misc$microns_per_pixel %||%
      object@misc$spatial_scales$microns_per_pixel %||%
      object@misc$scales$microns_per_pixel %||%
      tryCatch(Images(object)[[1]]@scale.factors$microns_per_pixel, error=function(e) NULL)
    
    s_low <- scale_lowres %||%
      object@misc$spatial_scales$tissue_lowres_scalef %||%
      object@misc$scales$tissue_lowres_scalef %||%
      tryCatch(Images(object)[[1]]@scale.factors$tissue_lowres_scalef, error=function(e) NULL)
    
    s_hires <- object@misc$spatial_scales$tissue_hires_scalef %||%
      object@misc$scales$tissue_hires_scalef %||%
      tryCatch(Images(object)[[1]]@scale.factors$tissue_hires_scalef, error=function(e) NULL)
    
    if (is.null(s_hires) && !is.null(s_low) &&
        all(c("imagecol_scaled","imagecol_hires") %in% colnames(plot_df))) {
      ratio <- plot_df$imagecol_scaled / plot_df$imagecol_hires
      ratio <- stats::median(ratio[is.finite(ratio)], na.rm = TRUE)
      if (is.finite(ratio) && ratio > 0) s_hires <- s_low / ratio
    }
    
    if (is.null(mpp)) stop("microns_per_pixel missing; provide microns_per_pixel = numeric.")
    if (is.null(s_low)) stop("tissue_lowres_scalef missing; provide scale_lowres = numeric.")
    
    # choose canvas scale based on coord system
    is_hires_coord <- any(grepl("hires", coord_cols, ignore.case = TRUE))
    s_use <- if (isTRUE(is_hires_coord)) {
      if (is.null(s_hires)) stop("tissue_hires_scalef missing while using hires coords.")
      s_hires
    } else {
      s_low
    }
    
    if (is.null(area_col) && "area_um2" %in% names(plot_df)) area_col <- "area_um2"
    if (is.null(bin_um)) {
      bin_um <- object@misc$bin_um %||%
        object@misc$spatial_scales$bin_size_um %||% NULL
    }
  } else if (is.data.frame(object)) {
    plot_df <- object
    mpp  <- microns_per_pixel
    s_low<- scale_lowres
    if (is.null(mpp) || is.null(s_low))
      stop("For data.frame, provide microns_per_pixel and scale_lowres.")
    s_use <- s_low
    if (is.null(area_col) && "area_um2" %in% names(plot_df)) area_col <- "area_um2"
  } else {
    stop("'object' must be Seurat or data.frame")
  }
  
  if (!all(coord_cols %in% names(plot_df)))
    stop("Missing coord columns: ", paste(coord_cols, collapse=", "))
  x_col <- coord_cols[1]; y_col <- coord_cols[2]
  
  if (!group_by %in% names(plot_df)) stop("Grouping column not found: ", group_by)
  plot_df[[group_by]] <- as.factor(plot_df[[group_by]])
  levs <- levels(plot_df[[group_by]])
  
  # ---------- 2) Colors ----------
  legend_name <- if (is.null(legend_title)) group_by else legend_title
  if (is.character(palette) && length(palette)==1) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) suppressWarnings(install.packages("RColorBrewer"))
    pal_name <- palette
    if (pal_name %in% rownames(RColorBrewer::brewer.pal.info)) {
      maxc <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
      cols <- RColorBrewer::brewer.pal(min(maxc, max(length(levs), 3)), pal_name)
      plot_colors <- rep(cols, length.out = length(levs)); names(plot_colors) <- levs
    } else stop("Unknown palette: ", palette)
  } else if (is.character(palette)) {
    plot_colors <- palette; names(plot_colors) <- levs
  } else if (is.function(palette)) {
    tmp <- palette(); plot_colors <- rep(tmp, length.out=length(levs)); names(plot_colors) <- levs
  } else stop("Invalid palette")
  
  # ---------- 3) Optional HE background ----------
  he_layer <- NULL
  if (isTRUE(show_he)) {
    if (!requireNamespace("png", quietly = TRUE)) suppressWarnings(install.packages("png"))
    he_res <- he_resolution
    if (he_res == "auto") {
      he_res <- if (any(grepl("hires", coord_cols, ignore.case = TRUE))) "hires" else "lowres"
    }
    outs_path <- he_outs_path %||% object@misc$OutsPath %||% object@misc$outs_path %||% object@misc$outs
    if (is.null(outs_path)) stop("show_he=TRUE but outs_path is unknown; pass he_outs_path=...")
    img_file <- file.path(outs_path, "spatial",
                          if (he_res=="lowres") "tissue_lowres_image.png" else "tissue_hires_image.png")
    if (!file.exists(img_file)) stop("HE image not found: ", img_file)
    img <- png::readPNG(img_file); Wimg <- ncol(img); Himg <- nrow(img)
    he_layer <- if (isTRUE(with_axes_grid)) {
      ggplot2::annotation_raster(img, xmin=0, xmax=Wimg, ymin=0, ymax=Himg, interpolate=TRUE)
    } else {
      ggplot2::annotation_raster(img, xmin=0, xmax=Wimg, ymin=-Himg, ymax=0, interpolate=TRUE)
    }
  }
  
  # ---------- 4) Geometry ----------
  library(ggplot2)
  flip_y <- function(y) -y
  y_map  <- function(v) if (isTRUE(with_axes_grid)) v else flip_y(v)
  
  use_real_tiles <- (spot_shape=="square" && size_mode=="auto" &&
                       ((!is.null(area_col) && area_col %in% names(plot_df)) || !is.null(bin_um)))
  
  if (use_real_tiles) {
    mpp_vec <- if (is.numeric(mpp)) rep(mpp, nrow(plot_df)) else plot_df[[mpp]]
    if (!is.null(area_col) && area_col %in% names(plot_df)) {
      side <- sqrt(pmax(plot_df[[area_col]], 0)) * (s_use / mpp_vec)
    } else {
      if (is.null(bin_um)) stop("size_mode='auto' needs area_col or bin_um")
      side <- (bin_um / mpp_vec) * s_use
    }
    cx <- plot_df[[x_col]]; cy <- plot_df[[y_col]]
    plot_df$.xmin <- cx - 0.5*side; plot_df$.xmax <- cx + 0.5*side
    plot_df$.ymin <- y_map(cy - 0.5*side); plot_df$.ymax <- y_map(cy + 0.5*side)
    
    p <- ggplot(plot_df, aes(xmin=.xmin, xmax=.xmax, ymin=.ymin, ymax=.ymax,
                             fill=.data[[group_by]])) +
      { if (!is.null(he_layer)) he_layer else NULL } +
      geom_rect(color=NA, ...) +
      scale_fill_manual(values = plot_colors, name = legend_name)
  } else {
    shape_val <- if (spot_shape=="circle") 16 else 15
    p <- ggplot(plot_df, aes(x=.data[[x_col]], y=y_map(.data[[y_col]]),
                             color=.data[[group_by]])) +
      { if (!is.null(he_layer)) he_layer else NULL } +
      {
        if (requireNamespace("scattermore", quietly = TRUE) && spot_shape=="circle")
          scattermore::geom_scattermore(pointsize = ptsize, pixels = c(2000,2000), ...)
        else
          geom_point(size=ptsize, shape=shape_val, ...)
      } +
      scale_color_manual(values = plot_colors, name = legend_name)
  }
  
  # ---------- 5) Axes / theme ----------
  plot_title <- if (is.null(title)) paste("Spatial Distribution of", group_by) else title
  px_to_um <- function(x) (x * (mpp / s_use))
  px_to_mm <- function(x) (px_to_um(x) / 1000)
  
  if (!isTRUE(with_axes_grid)) {
    # optional clipping to data extent to avoid "large empty HE" look
    if (isTRUE(clip_to_data) && isTRUE(show_he)) {
      xr  <- range(plot_df[[x_col]], na.rm = TRUE)
      yrp <- range(-plot_df[[y_col]], na.rm = TRUE)
      p <- p + coord_fixed(xlim = xr, ylim = yrp, expand = FALSE)
    } else {
      p <- p + coord_fixed()
    }
    p <- p + theme_void() +
      labs(title = plot_title) +
      theme(
        plot.title      = element_text(hjust=0.5, face="bold", size=14),
        legend.text     = element_text(size=8),
        legend.key.height = grid::unit(legend_key_height, "lines"),
        legend.key.width  = grid::unit(legend_key_width,  "lines"),
        legend.spacing.y  = grid::unit(legend_spacing_y,  "lines")
      ) +
      { if (!show_legend) guides(color="none", fill="none") else NULL }
  } else {
    xr <- range(plot_df[[x_col]], na.rm = TRUE)
    yr <- range(plot_df[[y_col]], na.rm = TRUE)
    
    sec_x <- switch(axes_secondary_unit,
                    "mm"  = sec_axis(~ px_to_mm(.), name = "mm (x)"),
                    "um"  = sec_axis(~ px_to_um(.), name = "µm (x)"),
                    "none"= waiver())
    sec_y <- switch(axes_secondary_unit,
                    "mm"  = sec_axis(~ px_to_mm(.), name = "mm (y)"),
                    "um"  = sec_axis(~ px_to_um(.), name = "µm (y)"),
                    "none"= waiver())
    
    p <- p +
      coord_fixed(expand = FALSE) +
      scale_x_continuous(
        name   = "Pixels (x)",
        limits = xr,
        breaks = scales::pretty_breaks(n = axes_breaks_n),
        sec.axis = sec_x
      ) +
      scale_y_reverse(
        name   = "Pixels (y)",
        limits = rev(yr),
        breaks = scales::pretty_breaks(n = axes_breaks_n),
        sec.axis = sec_y
      ) +
      theme_minimal(base_size = axes_text_size) +
      labs(title = plot_title) +
      theme(
        plot.title = element_text(hjust=0.5, face="bold", size=14),
        legend.text = element_text(size=8),
        panel.grid.major = element_line(color = grid_major_col, linewidth = 0.3),
        panel.grid.minor = element_line(color = if (grid_minor) grid_minor_col else NA, linewidth = 0.2)
      ) +
      { if (!show_legend) guides(color="none", fill="none") else NULL }
  }
  
  # ---------- 6) ROI overlay & selection ----------
  any_roi <- FALSE
  sel_idx <- rep(FALSE, nrow(plot_df))
  x_raw <- plot_df[[x_col]]
  y_raw <- plot_df[[y_col]]
  
  .circle_df <- function(cx, cy, rpx, n=240){
    t <- seq(0, 2*pi, length.out = n)
    data.frame(x = cx + rpx*cos(t), y = cy + rpx*sin(t))
  }
  um_to_px <- function(um) (um / mpp) * s_use
  
  .add_rect <- function(p, xmin, xmax, ymin, ymax) {
    ymin_plot <- if (with_axes_grid) ymin else -ymax
    ymax_plot <- if (with_axes_grid) ymax else -ymin
    if (is.na(ymin_plot) || is.na(ymax_plot)) return(p)
    if (ymin_plot > ymax_plot) { tmp <- ymin_plot; ymin_plot <- ymax_plot; ymax_plot <- tmp }
    p + annotate("rect",
                 xmin = xmin, xmax = xmax,
                 ymin = ymin_plot, ymax = ymax_plot,
                 fill = roi_fill, color = roi_line_color, linewidth = roi_linewidth)
  }
  
  if (!is.null(roi_rect_px) && length(roi_rect_px)==4) {
    any_roi <- TRUE
    xmin <- roi_rect_px[1]; xmax <- roi_rect_px[2]
    ymin <- roi_rect_px[3]; ymax <- roi_rect_px[4]
    sel_idx <- sel_idx | (x_raw >= xmin & x_raw <= xmax & y_raw >= ymin & y_raw <= ymax)
    p <- .add_rect(p, xmin, xmax, ymin, ymax)
  }
  if (!is.null(roi_rect_um) && length(roi_rect_um)==4) {
    any_roi <- TRUE
    rx <- um_to_px(roi_rect_um[c(1,2)])
    ry <- um_to_px(roi_rect_um[c(3,4)])
    xmin <- min(rx); xmax <- max(rx); ymin <- min(ry); ymax <- max(ry)
    sel_idx <- sel_idx | (x_raw >= xmin & x_raw <= xmax & y_raw >= ymin & y_raw <= ymax)
    p <- .add_rect(p, xmin, xmax, ymin, ymax)
  }
  if (!is.null(roi_radius_um) && (!is.null(roi_center_px) || !is.null(roi_center_um))) {
    any_roi <- TRUE
    cx <- if (!is.null(roi_center_px)) roi_center_px[1] else um_to_px(roi_center_um[1])
    cy <- if (!is.null(roi_center_px)) roi_center_px[2] else um_to_px(roi_center_um[2])
    rpx <- um_to_px(roi_radius_um)
    
    if (roi_shape == "circle") {
      sel_idx <- sel_idx | ((x_raw - cx)^2 + (y_raw - cy)^2 <= rpx^2)
      cir <- .circle_df(cx, cy, rpx)
      p <- p + geom_path(data = transform(cir, y = if (with_axes_grid) y else -y),
                         aes(x = x, y = y),
                         inherit.aes = FALSE, linewidth = roi_linewidth, color = roi_line_color)
      if (!is.na(roi_fill)) {
        p <- p + geom_polygon(data = transform(cir, y = if (with_axes_grid) y else -y),
                              aes(x = x, y = y),
                              inherit.aes = FALSE, fill = roi_fill, color = NA)
      }
    } else {
      xmin <- cx - rpx; xmax <- cx + rpx; ymin <- cy - rpx; ymax <- cy + rpx
      sel_idx <- sel_idx | (x_raw >= xmin & x_raw <= xmax & y_raw >= ymin & y_raw <= ymax)
      p <- .add_rect(p, xmin, xmax, ymin, ymax)
    }
  }
  
  roi_barcodes <- NULL
  if (any_roi && isTRUE(return_roi_barcodes)) {
    bc <- rownames(plot_df)
    if (is.null(bc) && "barcode" %in% names(plot_df)) bc <- plot_df$barcode
    roi_barcodes <- bc[which(sel_idx)]
  }
  
  # ---------- 7) Scale bar (uses s_use; supports µm thickness) ----------
  if (isTRUE(add_scale_bar)) {
    xs <- plot_df[[x_col]]
    xr <- range(xs, na.rm = TRUE); dx <- diff(xr)
    
    px_per_um <- (s_use / mpp)
    bar_um <- scale_bar_length_mm * 1000
    bar_px <- as.integer(max(10, round(bar_um * px_per_um)))
    
    # choose thickness in px (override by µm if provided)
    th_px <- as.integer(round(scale_bar_thickness_px))
    if (!is.null(scale_bar_thickness_um)) {
      th_px <- max(1L, as.integer(round(scale_bar_thickness_um * px_per_um)))
    }
    
    mx <- if (length(scale_bar_margin_frac) >= 1) scale_bar_margin_frac[1] else 0.04
    my <- if (length(scale_bar_margin_frac) >= 2) scale_bar_margin_frac[2] else mx
    
    if (isTRUE(with_axes_grid)) {
      yr <- range(plot_df[[y_col]], na.rm = TRUE); dy <- diff(yr)
      x0 <- if (scale_bar_position %in% c("bottomright","topright")) xr[2] - mx*dx - bar_px else xr[1] + mx*dx
      y0 <- if (scale_bar_position %in% c("bottomright","bottomleft"))
        yr[2] - my*dy - th_px  # bottom → large y (image coords)
      else
        yr[1] + my*dy
      y1 <- y0 + th_px
    } else {
      ys_plot <- -plot_df[[y_col]]
      yrp <- range(ys_plot, na.rm = TRUE); dyp <- diff(yrp)
      x0 <- if (scale_bar_position %in% c("bottomright","topright")) xr[2] - mx*dx - bar_px else xr[1] + mx*dx
      y0 <- if (scale_bar_position %in% c("bottomright","bottomleft"))
        yrp[1] + my*dyp
      else
        yrp[2] - my*dyp - th_px
      y1 <- y0 + th_px
    }
    
    nseg  <- max(1L, as.integer(round(scale_bar_n_segments)))
    gappx <- max(0L, as.integer(round(scale_bar_gap_px)))
    if (nseg == 1L) gappx <- 0L
    
    total_gap_px  <- gappx * (nseg - 1L)
    fill_px_total <- max(bar_px - total_gap_px, nseg)
    seg_fill_base <- floor(fill_px_total / nseg)
    rem <- fill_px_total - seg_fill_base * nseg
    seg_fill <- rep(seg_fill_base, nseg)
    if (rem > 0L) seg_fill[seq_len(rem)] <- seg_fill[seq_len(rem)] + 1L
    
    left <- x0
    seg_list <- vector("list", nseg)
    for (i in seq_len(nseg)) {
      right <- left + seg_fill[i]
      seg_list[[i]] <- data.frame(xmin=left, xmax=right, ymin=y0, ymax=y1)
      left <- right + if (i < nseg) gappx else 0L
    }
    sb_df <- do.call(rbind, seg_list)
    
    p <- p + ggplot2::geom_rect(
      data = sb_df,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = scale_bar_color, color = NA
    )
    
    if (!is.null(scale_bar_text) && nzchar(scale_bar_text)) {
      text_y <- if (isTRUE(with_axes_grid)) {
        if (scale_bar_position %in% c("bottomright","bottomleft")) y1 + scale_bar_text_offset_px else y0 - scale_bar_text_offset_px
      } else {
        if (scale_bar_position %in% c("bottomright","bottomleft")) y1 + scale_bar_text_offset_px else y0 - scale_bar_text_offset_px
      }
      p <- p + ggplot2::geom_text(
        data = data.frame(x = x0 + bar_px/2, y = text_y, lab = scale_bar_text),
        ggplot2::aes(x = x, y = y, label = lab),
        inherit.aes = FALSE,
        size = scale_bar_text_size, color = scale_bar_text_color,
        fontface = scale_bar_text_face,
        vjust = if (scale_bar_position %in% c("bottomright","bottomleft")) 0 else 1
      )
    }
  }
  
  invisible(list(plot = p, roi_barcodes = roi_barcodes, data = plot_df))
}






make_hier_palette <- function(
    groups,
    base_colors,
    major_sort   = c("alpha","as_is","custom"),
    custom_major_order = NULL,   # 当 major_sort="custom" 时提供
    within_sort  = c("numeric","alpha","none"),
    lighten_amt  = 0.40,
    darken_amt   = 0.40,
    fallback_scheme = c("Glasbey","Dark3"),
    seed = 42
){
  stopifnot(is.character(groups), length(groups) > 0)
  stopifnot(is.character(base_colors), !is.null(names(base_colors)))
  if (!requireNamespace("colorspace", quietly = TRUE))
    stop("Please install 'colorspace'. (install.packages('colorspace'))")
  
  major_sort   <- match.arg(major_sort)
  within_sort  <- match.arg(within_sort)
  fallback_scheme <- match.arg(fallback_scheme)
  
  groups <- as.character(groups)
  uniq_groups <- unique(groups)
  
  # --- helpers ---
  .rx_escape <- function(s) gsub("([.|()\\^{}+$*?\\\\\\[\\]-])", "\\\\\\1", s)
  
  # 1) 大类顺序
  majors <- names(base_colors)
  if (major_sort == "alpha") {
    majors <- sort(majors)
  } else if (major_sort == "custom") {
    stopifnot(!is.null(custom_major_order))
    majors <- as.character(custom_major_order)
  } # "as_is" 保持传入顺序
  
  pal <- character(0)
  assigned <- character(0)
  
  for (maj in majors) {
    maj_rx <- paste0("^", .rx_escape(maj))
    sub <- uniq_groups[grepl(maj_rx, uniq_groups)]
    if (length(sub) == 0) next
    
    # 2) 类内顺序
    if (within_sort == "numeric") {
      # 取末尾数字；无数字的为 NA
      m <- regexpr("(\\d+)$", sub)
      num <- suppressWarnings(as.numeric(ifelse(m > 0, regmatches(sub, m), NA)))
      # 先有数字者按数字升序；再无数字者按字母序
      ord <- order(is.na(num), num, sub, na.last = TRUE)
      sub <- sub[ord]
    } else if (within_sort == "alpha") {
      sub <- sort(sub)
    } # "none" 则保持原有出现顺序
    
    # 3) 为该大类生成亮→暗梯度并按顺序赋色
    base <- base_colors[[maj]]
    col_fun <- grDevices::colorRampPalette(c(
      colorspace::lighten(base, lighten_amt),
      colorspace::darken(base,  darken_amt)
    ))
    cols <- col_fun(length(sub))
    names(cols) <- sub
    
    pal <- c(pal, cols)
    assigned <- c(assigned, sub)
  }
  
  # 4) 未匹配到任何大类前缀的分组（fallback），按字母序再补色
  rest <- setdiff(uniq_groups, assigned)
  if (length(rest) > 0) {
    rest <- sort(rest)
    set.seed(seed)
    fb <- if (fallback_scheme == "Glasbey") {
      grDevices::hcl.colors(length(rest), palette = "Glasbey")
    } else {
      colorspace::qualitative_hcl(length(rest), palette = "Dark 3")
    }
    names(fb) <- rest
    pal <- c(pal, fb)
  }
  
  # 返回按“先大类、后类内”的整体顺序的具名向量
  pal
}

# ==== Seurat?????? ====












# 如果没有安装 Polychrome，请先运行 install.packages("Polychrome")
library(Polychrome)

# 如果没有安装 Polychrome，请先运行：install.packages("Polychrome")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

DimPlotNice <- function(
    object,
    reduction,
    group_by,
    palette = c("viridis","turbo","okabe-ito","scico-roma","glasbey"),
    sort_levels = c("none","alpha","numeric","custom"),
    custom_order = NULL,
    sort_desc   = FALSE,
    pt.size     = 0.15,
    raster      = TRUE,
    label       = TRUE,
    repel       = TRUE,
    label_size  = 3.5,
    legend_title = NULL,
    legend_ncol  = 1,
    legend_position = "right",
    legend_pt_size = 3,        # 默认固定为3（独立于 pt.size）
    title = NULL,
    # === 新增坐标轴参数 ===
    show_axes   = c("none", "simple", "scale"),  # none=隐藏, simple=简单轴, scale=带标尺
    axis_title_x = NULL,
    axis_title_y = NULL,
    axis_text_size = 10,       # 坐标轴文本大小
    axis_title_size = 11,      # 坐标轴标题大小
    scale_bar_length = NULL,   # 标尺长度（NULL=自动）
    scale_bar_units = "units", # 标尺单位文本
    scale_bar_position = c("bottom_right", "bottom_left", "top_right", "top_left"),
    scale_bar_height = 0.5,    # 标尺线粗细
    scale_bar_text_size = 3    # 标尺文本大小
){
  # ---------- 基本检查 ----------
  if (!reduction %in% names(object@reductions))
    stop("Reduction '", reduction, "' not found in object@reductions.")
  if (!group_by %in% colnames(object@meta.data))
    stop("group_by '", group_by, "' not found in meta.data.")
  
  show_axes <- match.arg(show_axes)
  scale_bar_position <- match.arg(scale_bar_position)
  sort_levels <- match.arg(sort_levels)
  
  vals_chr  <- as.character(object@meta.data[[group_by]])
  uniq_vals <- unique(vals_chr)
  
  # ---------- level 顺序 ----------
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
  
  # ---------- palette 解析 ----------
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
      # 颜色向量（不一定带名字）
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
      # 偏好列表 tokens
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
    # 兜底：HCL
    return(function(n) grDevices::hcl.colors(n, "Spectral"))
  }
  
  # 1) 如果 palette 是带名字的颜色向量：按 breaks 名字重排（缺的自动补色）
  if (is.character(palette) && are_valid_colors(palette) && !is.null(names(palette))) {
    cols_vec <- unname(palette[breaks])
    # 对缺失类补色
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
    # 2) 否则按原策略：函数 / 偏好列表 / 无名颜色向量
    pal_fun <- get_palette_fun(palette)
    cols_vec <- pal_fun(length(breaks))
    names(cols_vec) <- breaks
  }
  
  # ---------- 固定因子顺序 ----------
  object@meta.data[[group_by]] <- factor(vals_chr, levels = breaks)
  
  # ---------- 文案 ----------
  if (is.null(legend_title)) legend_title <- group_by
  if (is.null(title))        title <- paste0("UMAP — ", group_by)
  
  # ---------- 绘图 ----------
  p <- Seurat::DimPlot(
    object,
    reduction = reduction,
    group.by  = group_by,
    label     = label,
    repel     = repel,
    label.size= label_size,
    pt.size   = pt.size,
    raster    = raster,
    cols      = unname(cols_vec[breaks])
  ) +
    ggplot2::coord_equal() +
    ggplot2::ggtitle(title) +
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
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        ncol = legend_ncol, 
        byrow = TRUE,
        override.aes = list(
          size = legend_pt_size,  # 固定图例点大小，不受 pt.size 影响
          alpha = 1
        )
      )
    ) +
    ggplot2::labs(color = legend_title)
  
  # ---------- 坐标轴处理 ----------
  if (show_axes == "simple") {
    # 简单坐标轴：显示坐标轴和标题
    dims <- colnames(Seurat::Embeddings(object, reduction))
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
    # 带标尺的坐标轴
    emb <- Seurat::Embeddings(object, reduction)
    x_range <- range(emb[, 1])
    y_range <- range(emb[, 2])
    
    # 自动计算标尺长度（如果未指定）
    if (is.null(scale_bar_length)) {
      x_span <- diff(x_range)
      # 标尺长度为坐标范围的 1/5（可调整）
      scale_bar_length <- round(x_span / 5, 1)
      # 美化：调整到整数或 0.5 的倍数
      if (scale_bar_length >= 2) {
        scale_bar_length <- round(scale_bar_length)
      } else if (scale_bar_length >= 1) {
        scale_bar_length <- round(scale_bar_length * 2) / 2
      } else {
        scale_bar_length <- round(scale_bar_length * 10) / 10
      }
    }
    
    # 计算标尺位置
    x_margin <- diff(x_range) * 0.05
    y_margin <- diff(y_range) * 0.05
    
    bar_pos <- switch(
      scale_bar_position,
      "bottom_right" = list(
        x = x_range[2] - x_margin - scale_bar_length,
        y = y_range[1] + y_margin,
        hjust = 0
      ),
      "bottom_left" = list(
        x = x_range[1] + x_margin,
        y = y_range[1] + y_margin,
        hjust = 0
      ),
      "top_right" = list(
        x = x_range[2] - x_margin - scale_bar_length,
        y = y_range[2] - y_margin,
        hjust = 0
      ),
      "top_left" = list(
        x = x_range[1] + x_margin,
        y = y_range[2] - y_margin,
        hjust = 0
      )
    )
    
    # 标尺数据
    scale_data <- data.frame(
      x = bar_pos$x,
      xend = bar_pos$x + scale_bar_length,
      y = bar_pos$y,
      yend = bar_pos$y
    )
    
    # 标尺文本位置（居中）
    text_data <- data.frame(
      x = bar_pos$x + scale_bar_length / 2,
      y = bar_pos$y + diff(y_range) * 0.02,  # 文本在标尺上方
      label = paste0(scale_bar_length, " ", scale_bar_units)
    )
    
    # 添加标尺到图上
    p <- p +
      ggplot2::geom_segment(
        data = scale_data,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
        color = "black",
        linewidth = scale_bar_height,
        inherit.aes = FALSE
      ) +
      # 标尺两端的竖线
      ggplot2::geom_segment(
        data = scale_data,
        ggplot2::aes(x = x, y = y - diff(y_range) * 0.005, 
                     xend = x, yend = y + diff(y_range) * 0.005),
        color = "black",
        linewidth = scale_bar_height,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_segment(
        data = scale_data,
        ggplot2::aes(x = xend, y = y - diff(y_range) * 0.005, 
                     xend = xend, yend = y + diff(y_range) * 0.005),
        color = "black",
        linewidth = scale_bar_height,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_text(
        data = text_data,
        ggplot2::aes(x = x, y = y, label = label),
        size = scale_bar_text_size,
        inherit.aes = FALSE,
        fontface = "bold"
      )
  }
  # show_axes == "none" 时不做任何修改，保持默认（隐藏坐标轴）
  
  return(p)
}






suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(FNN)      # neighbor-consistency
})

.get_counts <- function(obj, assay) {
  tryCatch(GetAssayData(obj, assay = assay, layer = "counts"),
           error = function(e) tryCatch(GetAssayData(obj, assay = assay, slot = "counts"),
                                        error = function(e2) stop("Cannot get counts for assay=", assay)))
}





# 剔除易主导HVG的技术基因: 线粒体/核糖体/血红蛋白
.drop_tech_genes <- function(genes) {
  bad <- grepl("^MT-|^mt-|^RPL|^RPS|^Rp[sl]|^HBA|^HBB", genes, perl = TRUE)
  genes[!bad]
}





# 空间坐标优先（scaled>raw），否则退化到 UMAP
.pick_coords_for_qc <- function(obj, reduction = "seurat_umap") {
  md <- obj@meta.data
  if (all(c("imagecol_scaled","imagerow_scaled") %in% colnames(md))) {
    as.matrix(md[, c("imagecol_scaled","imagerow_scaled")])
  } else if (all(c("imagecol","imagerow") %in% colnames(md))) {
    as.matrix(md[, c("imagecol","imagerow")])
  } else {
    Embeddings(obj, reduction = reduction)[, 1:2, drop = FALSE]
  }
}





# 邻域一致性：每个点的KNN里，多数派标签的比例，最后取均值
neighbor_consistency <- function(labels, coords, k = 20) {
  idx <- FNN::get.knn(coords, k = k)$nn.index
  frac <- vapply(seq_len(nrow(coords)), function(i) {
    nnlab <- labels[idx[i, ]]
    tab <- sort(table(nnlab), decreasing = TRUE)
    as.numeric(tab[1]) / length(nnlab)
  }, numeric(1))
  mean(frac, na.rm = TRUE)
}




# -------- v4/v5 通吃：释放 scale.data 的小工具 --------
free_scale_data <- function(obj, assay_name) {
  # 推荐用 [[ ]] 取 assay，避免 $/@ 混淆
  aa <- obj[[assay_name]]
  
  # v5: Assay5 有 layers
  if (inherits(aa, "Assay5")) {
    # 安全检查是否存在名为 "scale.data" 的 layer
    has_layer <- FALSE
    lyr_names <- tryCatch(SeuratObject::Layers(aa), error = function(e) character())
    if ("scale.data" %in% lyr_names) has_layer <- TRUE
    
    if (has_layer) {
      # 从 layers 列表里移除该层
      aa@layers <- aa@layers[setdiff(names(aa@layers), "scale.data")]
      obj[[assay_name]] <- aa
    }
  } else {
    # v3/v4: 可能有 slot scale.data
    if ("scale.data" %in% slotNames(aa)) {
      aa@scale.data <- matrix(0, 0, 0)
      obj[[assay_name]] <- aa
    }
  }
  invisible(obj)
}





# -----------------------------
# Level-1：全局表达-only 聚类（可扫分辨率、自动挑最佳）
# -----------------------------
# 返回：更新后的对象（写入 `Seurat_L1_cluster`）、分辨率扫描表、被选中的分辨率与一致性分数
# 依赖：.get_counts, .drop_tech_genes, .pick_coords_for_qc, neighbor_consistency, free_scale_data
# - free_scale_data(obj, assay_name): 之前给过 v4/v5 通吃版本

seurat_hd_level1 <- function(
    object,
    assay              = "Spatial",
    method             = c("lognorm","sct"),
    n_hvg              = 3000,
    dims_use           = 1:30,
    k_nn               = 30,
    res_main           = 0.6,
    do_resolution_sweep= TRUE,
    res_grid           = c(0.4, 0.6, 0.8),
    pick_by            = c("neighbor_consistency"),
    qc_k               = 20,
    regress_vars       = c("percent.mt"),
    run_umap           = TRUE,
    umap_name          = "seurat_umap",
    seed               = 1,
    verbose            = FALSE
){
  stopifnot(inherits(object, "Seurat"))
  set.seed(seed)
  method <- match.arg(method)
  

  DefaultAssay(object) <- assay
  C <- .get_counts(object, assay)
  if (!"percent.mt" %in% colnames(object@meta.data)) {
    mito_idx <- grepl("^MT-|^mt-", rownames(C))
    object$percent.mt <- if (any(mito_idx))
      100 * Matrix::colSums(C[mito_idx, , drop = FALSE]) / pmax(Matrix::colSums(C), 1) else 0
  }
  
  # ---- 归一化 + HVG ----
  if (method == "lognorm") {
    object <- NormalizeData(object, assay = assay,
                            normalization.method = "LogNormalize",
                            scale.factor = 1e4, verbose = verbose)
    object <- FindVariableFeatures(object, assay = assay,
                                   nfeatures = n_hvg, selection.method = "vst",
                                   verbose = verbose)
    VariableFeatures(object, assay = assay) <- .drop_tech_genes(VariableFeatures(object, assay = assay))
    
    # 仅缩放 HVG，避免内存爆炸
    object <- ScaleData(object, assay = assay,
                        features = VariableFeatures(object, assay = assay),
                        vars.to.regress = regress_vars,
                        block.size = 2000, verbose = verbose)
    
    # 使用哪个 assay 进入 PCA
    assay_used <- assay
    
  } else { # method == "sct"
    suppressPackageStartupMessages(requireNamespace("glmGamPoi"))
    object <- SCTransform(object, assay = assay, verbose = verbose,
                          variable.features.n = n_hvg, method = "glmGamPoi")
    DefaultAssay(object) <- "SCT"
    assay_used <- "SCT"
  }
  
  # ---- 确保 PCA 维度 >= max(dims_use)，只跑一次 ----
  .ensure_pca_dims <- function(obj, reduction = "pca", npcs_target = 30, assay_name = DefaultAssay(obj)) {
    need_rerun <- !(reduction %in% names(obj@reductions))
    if (!need_rerun) {
      cur <- ncol(Embeddings(obj, reduction))
      if (is.null(cur) || cur < npcs_target) need_rerun <- TRUE
    }
    if (need_rerun) {
      # 若 HVG 太少，再找一次（防守式）
      if (length(VariableFeatures(obj, assay = assay_name)) < 50) {
        obj <- FindVariableFeatures(obj, assay = assay_name, nfeatures = max(2000, npcs_target*50), verbose = FALSE)
      }
      # 仅缩放 HVG（SCT 分支也安全：ScaleData 只影响缩放矩阵，不会破坏 SCT 归一化）
      obj <- ScaleData(obj, assay = assay_name,
                       features = VariableFeatures(obj, assay = assay_name),
                       block.size = 2000, verbose = FALSE)
      obj <- RunPCA(obj, assay = assay_name,
                    features = VariableFeatures(obj, assay = assay_name),
                    npcs = npcs_target, reduction.name = reduction, verbose = FALSE)
    }
    obj
  }
  
  object <- .ensure_pca_dims(object,
                             reduction   = "pca",
                             npcs_target = max(dims_use),
                             assay_name  = assay_used)
  
  # 释放 scale.data（lognorm 分支才需要；SCT 本身不存巨大 scale.data）
  if (method == "lognorm") {
    object <- free_scale_data(object, assay_name = assay_used)
    gc()
  }
  
  # ---- KNN/SNN 图 & 聚类扫描 ----
  object <- FindNeighbors(object, reduction = "pca", dims = dims_use,
                          k.param = k_nn, verbose = verbose)
  
  .cluster_and_score <- function(obj, res) {
    obj <- FindClusters(obj, resolution = res, verbose = verbose)
    labs <- Idents(obj)
    if (run_umap && !(umap_name %in% Reductions(obj))) {
      obj <- RunUMAP(obj, reduction = "pca", dims = dims_use,
                     reduction.name = umap_name, verbose = verbose, seed.use = seed)
    }
    coords <- .pick_coords_for_qc(obj, reduction = umap_name)
    score  <- neighbor_consistency(labs, coords, k = qc_k)
    list(obj = obj, score = score, labels = labs)
  }
  
  sweep_tbl <- NULL
  best_res  <- res_main
  best_obj  <- NULL
  best_score<- NA_real_
  
  if (isTRUE(do_resolution_sweep)) {
    res_vec <- unique(sort(res_grid))
    tmp_rows <- vector("list", length(res_vec))
    for (i in seq_along(res_vec)) {
      r <- res_vec[i]
      tmp <- .cluster_and_score(object, r)
      tmp_rows[[i]] <- data.frame(resolution = r,
                                  K = length(table(tmp$labels)),
                                  neighbor_consistency = tmp$score)
      # 选分辨率：以 neighbor_consistency 最大为先，并列取更靠近中位的
      if (is.na(best_score) || tmp$score > best_score ||
          (abs(tmp$score - best_score) < 1e-6 &&
           abs(r - median(res_vec)) < abs(best_res - median(res_vec)))) {
        best_score <- tmp$score; best_res <- r; best_obj <- tmp$obj
      }
    }
    sweep_tbl <- dplyr::bind_rows(tmp_rows)
  } else {
    out <- .cluster_and_score(object, res_main)
    best_obj   <- out$obj
    best_score <- out$score
    sweep_tbl  <- data.frame(resolution = res_main,
                             K = length(table(out$labels)),
                             neighbor_consistency = out$score)
  }
  
  object <- best_obj
  object$Seurat_L1_cluster <- factor(as.character(Idents(object)))
  object@misc$seurat_l1 <- list(
    best_resolution = best_res,
    neighbor_consistency = best_score,
    sweep = sweep_tbl,
    params = list(method = method, dims_use = dims_use, k_nn = k_nn,
                  regress_vars = regress_vars, qc_k = qc_k, umap = umap_name,
                  assay_used = assay_used)
  )
  message(sprintf("[L1] chosen resolution = %.2f | neighbor-consistency = %.3f",
                  best_res, best_score))
  return(object)
}




suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(FNN)
})

.check_pick_col <- function(x){
  # returns a plain vector (character/numeric) from vector/factor/data.frame
  if (is.data.frame(x)) x <- x[,1, drop=TRUE]
  if (is.factor(x)) x <- as.character(x)
  x
}





check_clusters_quality <- function(
    object,
    assay              = "Spatial",
    cluster_col        = "seurat_clusters",
    k_spatial          = 6,         # neighbors for spatial continuity
    mito_prefix        = "^MT-",    # mitochondrial gene prefix
    compute_complexity = FALSE,     # top1/entropy (slow on very large data) -> OFF by default
    compute_top1_rbc   = FALSE,     # whether to check "top gene is RBC" (can be slow) -> OFF
    rbc_floor_hi       = 0.10,      # minimum "high RBC" fraction floor
    mito_floor_hi      = 0.10,      # minimum "high mito" fraction floor
    verbose            = TRUE
){
  stopifnot(assay %in% names(object@assays))
  if (!cluster_col %in% colnames(object@meta.data))
    stop("cluster_col not found in meta.data: ", cluster_col)
  
  DefaultAssay(object) <- assay
  
  # ---- counts matrix ----
  mat <- GetAssayData(object, assay = assay, slot = "counts")
  nCount <- if ("nCount_Spatial" %in% colnames(object@meta.data)) .check_pick_col(object$nCount_Spatial) else Matrix::colSums(mat)
  nGene  <- if ("nFeature_Spatial" %in% colnames(object@meta.data)) .check_pick_col(object$nFeature_Spatial) else Matrix::colSums(mat > 0)
  
  # ---- gene sets ----
  hb_genes   <- intersect(c("HBA1","HBA2","HBB","HBD","HBM","HBQ1"), rownames(mat))
  mito_genes <- rownames(mat)[grepl(mito_prefix, rownames(mat))]
  ribo_genes <- rownames(mat)[grepl("^RPL|^RPS", rownames(mat))]
  
  frac <- function(g) if (length(g) == 0) rep(0, ncol(mat)) else Matrix::colSums(mat[g,,drop=FALSE]) / pmax(nCount, 1)
  
  frac_rbc  <- frac(hb_genes)
  frac_mito <- frac(mito_genes)
  frac_ribo <- frac(ribo_genes)
  
  # ---- optional complexity metrics ----
  if (compute_complexity) {
    top1_frac <- apply(mat, 2, function(x){ s <- sum(x); if (s == 0) 0 else max(x)/s })
    entropy   <- apply(mat, 2, function(x){
      s <- sum(x); if (s == 0) return(0)
      p <- x[x>0]/s; -sum(p*log2(p))
    })
  } else {
    top1_frac <- rep(NA_real_, ncol(mat))
    entropy   <- rep(NA_real_, ncol(mat))
  }
  genes_per_1kUMI <- nGene / (nCount/1000 + 1e-9)
  
  # ---- spatial continuity (image coords) ----
  md <- object@meta.data
  cc <- c("imagecol_scaled","imagerow_scaled")
  if (!all(cc %in% colnames(md))) cc <- c("imagecol","imagerow")
  if (!all(cc %in% colnames(md)))
    stop("Missing image coordinates: need imagecol[_scaled] and imagerow[_scaled].")
  
  XY <- as.matrix(md[, cc]); colnames(XY) <- c("x","y")
  
  # ---- cluster labels (safe) ----
  cl <- .check_pick_col(md[[cluster_col]])
  names(cl) <- colnames(object)
  
  # ---- auto-adjust k for safety ----
  k_eff <- max(1L, min(as.integer(k_spatial), nrow(XY) - 1L))
  
  kn <- get.knn(XY, k = k_eff)
  nn_index <- kn$nn.index
  same_label_share <- sapply(seq_len(nrow(nn_index)), function(i){
    mean(cl[nn_index[i,]] == cl[i])
  })
  names(same_label_share) <- colnames(object)
  
  # ---- assemble per-cell QC ----
  df <- data.frame(
    cell = colnames(object),
    cluster = cl,
    nCount = nCount,
    nGene  = nGene,
    frac_rbc = frac_rbc,
    frac_mito = frac_mito,
    frac_ribo = frac_ribo,
    top1_frac = top1_frac,
    entropy = entropy,
    genes_per_1kUMI = genes_per_1kUMI,
    spatial_continuity = same_label_share,
    stringsAsFactors = FALSE
  )
  
  # ---- robust IQR thresholds (global) ----
  thr_low <- function(v) quantile(v, 0.25, na.rm=TRUE) - 1.5*IQR(v, na.rm=TRUE)
  thr_hi  <- function(v) quantile(v, 0.75, na.rm=TRUE) + 1.5*IQR(v, na.rm=TRUE)
  
  thr <- list(
    nCount_low  = thr_low(df$nCount),
    nGene_low   = thr_low(df$nGene),
    rbc_hi      = max(thr_hi(df$frac_rbc),  rbc_floor_hi),
    mito_hi     = max(thr_hi(df$frac_mito), mito_floor_hi),
    top1_hi     = if (all(is.na(df$top1_frac))) NA_real_ else thr_hi(df$top1_frac),
    entropy_low = if (all(is.na(df$entropy)))   NA_real_ else thr_low(df$entropy),
    gpr1k_low   = thr_low(df$genes_per_1kUMI),
    cont_low    = thr_low(df$spatial_continuity)
  )
  
  # ---- per-cell flags ----
  df <- df %>%
    mutate(
      flag_low_counts   = nCount <= thr$nCount_low,
      flag_low_genes    = nGene  <= thr$nGene_low,
      flag_high_rbc     = frac_rbc >= thr$rbc_hi,
      flag_high_mito    = frac_mito >= thr$mito_hi,
      flag_high_top1    = if (!is.na(thr$top1_hi)) top1_frac >= thr$top1_hi else FALSE,
      flag_low_entropy  = if (!is.na(thr$entropy_low)) entropy <= thr$entropy_low else FALSE,
      flag_low_cont     = spatial_continuity <= thr$cont_low
    ) %>%
    mutate(flag_any = flag_low_counts | flag_low_genes | flag_high_rbc |
             flag_high_mito | flag_high_top1 | flag_low_entropy | flag_low_cont)
  
  # ---- RBC top1 dominance (optional) ----
  if (compute_top1_rbc && length(hb_genes) > 0) {
    top_ix <- apply(mat, 2, function(x){ if (sum(x) == 0) NA_integer_ else which.max(x) })
    top_is_rbc <- !is.na(top_ix) & rownames(mat)[top_ix] %in% hb_genes
  } else {
    top_is_rbc <- rep(NA, ncol(mat))
  }
  
  # ---- per-cluster summary ----
  sum_by_cl <- df %>%
    group_by(cluster) %>%
    summarise(
      n = n(),
      med_nCount = median(nCount), med_nGene = median(nGene),
      med_rbc = median(frac_rbc), med_mito = median(frac_mito), med_ribo = median(frac_ribo),
      med_top1 = median(top1_frac, na.rm=TRUE),
      med_entropy = median(entropy, na.rm=TRUE),
      med_gpr1k = median(genes_per_1kUMI),
      med_cont  = median(spatial_continuity),
      bad_rate  = mean(flag_any),
      rbc_hi_rate = mean(flag_high_rbc),
      mito_hi_rate = mean(flag_high_mito),
      low_cont_rate = mean(flag_low_cont),
      rbc_top1_rate = if (all(is.na(top_is_rbc))) NA_real_ else mean(top_is_rbc, na.rm=TRUE)
    ) %>% ungroup()
  
  # ---- risk scoring & decision ----
  add_point <- function(cond, pts) ifelse(cond, pts, 0)
  sum_by_cl <- sum_by_cl %>%
    mutate(
      sc = 0 +
        add_point(med_nGene  <= thr$nGene_low, 1.5) +
        add_point(med_nCount <= thr$nCount_low, 1.0) +
        add_point(rbc_hi_rate >= 0.25,         1.0) +
        add_point(!is.na(rbc_top1_rate) & rbc_top1_rate >= 0.25, 1.0) +
        add_point(mito_hi_rate >= 0.25,        0.5) +
        add_point(!is.na(thr$top1_hi) & med_top1 >= thr$top1_hi, 1.0) +
        add_point(!is.na(thr$entropy_low) & med_entropy <= thr$entropy_low, 1.0) +
        add_point(med_gpr1k <= thr$gpr1k_low,  0.5) +
        add_point(med_cont  <= thr$cont_low,   1.0)
    ) %>%
    mutate(
      decision = dplyr::case_when(
        sc >= 3.0 ~ "LowQuality",
        sc >= 1.5 ~ "Suspect",
        TRUE      ~ "OK"
      ),
      reason = paste0(
        ifelse(med_nGene  <= thr$nGene_low,  "[low nGene] ", ""),
        ifelse(med_nCount <= thr$nCount_low, "[low nCount] ", ""),
        ifelse(rbc_hi_rate >= 0.25,          "[RBC-high many] ", ""),
        ifelse(!is.na(rbc_top1_rate) & rbc_top1_rate >= 0.25, "[RBC-top1 many] ", ""),
        ifelse(mito_hi_rate >= 0.25,         "[mito-high many] ", ""),
        ifelse(!is.na(thr$top1_hi) & med_top1 >= thr$top1_hi, "[high top1] ", ""),
        ifelse(!is.na(thr$entropy_low) & med_entropy <= thr$entropy_low, "[low entropy] ", ""),
        ifelse(med_gpr1k <= thr$gpr1k_low,   "[low genes/1kUMI] ", ""),
        ifelse(med_cont  <= thr$cont_low,    "[low spatial continuity] ", "")
      )
    ) %>%
    arrange(desc(sc), cluster)
  
  if (verbose) {
    message("\n[QC thresholds] IQR-based (global):")
    print(thr)
    message("\n[Cluster-level QC summary ranked by score]")
    print(sum_by_cl %>% select(cluster, n, sc, decision, reason,
                               med_nCount, med_nGene, med_rbc, rbc_hi_rate,
                               med_cont, low_cont_rate))
  }
  
  list(per_cell = df, per_cluster = sum_by_cl, thresholds = thr)
}





# ------------------- helpers -------------------
.sanitize <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

.drop_tech_genes <- function(genes){
  bad <- grepl("^MT-", genes, ignore.case=TRUE) |
    grepl("^RPL|^RPS", genes) |
    genes %in% c("MALAT1","XIST","RN7SL1","RN7SL2","KCNQ1OT1")
  genes[!bad]
}





# 简易近邻一致性评分（用于扫分辨率可选）
.neighbor_consistency <- function(labels, coords, k = 15){
  if (!requireNamespace("FNN", quietly = TRUE)) {
    warning("FNN 未安装，跳过一致性评分（返回 0.5）。install.packages('FNN') 可加速/启用。")
    return(0.5)
  }
  kn <- FNN::get.knn(coords, k = k)$nn.index
  lab <- as.integer(factor(labels))
  same <- vapply(seq_len(nrow(kn)), function(i){
    mean(lab[kn[i,]] == lab[i])
  }, 0.0)
  mean(same)
}





# ------------------- main function -------------------
suppressPackageStartupMessages({ library(Seurat) })

# 仅用于 HVG：屏蔽技术/外泌家族（表达仍可用于打分/可视化）
.drop_tech_genes <- function(genes) {
  genes[!grepl("^MT-|^mt-|^RPL|^RPS|^HBA|^HBB|^IGH|^IGK|^IGL", genes)]
}




.sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# 可选：PCA坐标邻域一致性（需要 FNN）
.neighbor_consistency <- function(labels, coords, k = 15) {
  if (!requireNamespace("FNN", quietly = TRUE)) stop("FNN is required for sweep.")
  idx <- FNN::get.knn(coords, k = k)$nn.index
  frac <- vapply(seq_len(nrow(coords)), function(i){
    nnlab <- labels[idx[i,]]; tab <- sort(table(nnlab), TRUE)
    as.numeric(tab[1]) / length(nnlab)
  }, numeric(1))
  mean(frac, na.rm = TRUE)
}





# 统一、安全地往 meta.data 建列（长度对齐）
.ensure_meta_cols <- function(obj, cols) {
  for (cn in cols) {
    if (!cn %in% colnames(obj@meta.data)) {
      obj@meta.data[[cn]] <- rep(NA_character_, ncol(obj))
    }
  }
  obj
}





suppressPackageStartupMessages({ library(Seurat) })

# 可选：技术/外泌家族从 HVG 里去掉（表达仍保留）
.drop_tech_genes <- function(genes){
  genes[!grepl("^MT-|^RPL|^RPS|^HBA|^HBB|^IGH|^IGK|^IGL", genes)]
}





# 小工具
.sanitize <- function(x) gsub("[^A-Za-z0-9]+","_", x)
.ensure_meta_cols_if_missing <- function(obj, cols){
  for (nm in cols) if (!nm %in% colnames(obj@meta.data)) obj@meta.data[[nm]] <- NA_character_
  obj
}





suppressPackageStartupMessages({ library(Seurat) })

.drop_tech_genes <- function(genes){
  genes[!grepl("^MT-|^RPL|^RPS|^HBA|^HBB|^IGH|^IGK|^IGL", genes)]
}




.sanitize <- function(x) gsub("[^A-Za-z0-9]+","_", x)
.ensure_meta_cols_if_missing <- function(obj, cols){
  for (nm in cols) if (!nm %in% colnames(obj@meta.data)) obj@meta.data[[nm]] <- NA_character_
  obj
}





seurat_hd_level2 <- function(
    object,
    assay          = "Spatial",
    coarse_col,                                # e.g. "Seurat_L1_label"
    targets       = c("CAF"),
    n_hvg         = 3000,
    dims_use      = 1:30,
    k_nn          = 30,
    res_main      = 0.6,
    do_resolution_sweep = FALSE,
    res_grid      = c(0.4, 0.6, 0.8),
    min_cells     = 200,
    regress_vars  = c("percent.mt"),
    run_umap      = TRUE,
    seed          = 1,
    verbose       = FALSE,
    per_class_params = NULL,
    return_subobjects = FALSE,
    cluster_algorithm = NULL,                  # 1=SLM, 3=Leiden（按你的 Seurat 版本）
    keep_existing = TRUE,                      # TRUE=增量（保留旧 L2）；FALSE=覆盖
    save_umap_reduction = TRUE,                # 保存为 reductions：umap_<类>
    save_umap_meta      = TRUE,                # 同步写 meta 列：L2UMAP_<类>_{1,2}
    umap_meta_prefix    = "L2UMAP"
){
  stopifnot(inherits(object, "Seurat"))
  stopifnot(coarse_col %in% colnames(object@meta.data))
  set.seed(seed)
  DefaultAssay(object) <- assay
  
  # 细胞名与 meta 行名一致性
  stopifnot(identical(colnames(object), rownames(object@meta.data)))
  
  # 读取 coarse 标签
  coarse_vec <- as.character(object@meta.data[[coarse_col]])
  
  # 合并默认参数 + per-class 覆盖
  .get_params_for <- function(class_name){
    base <- list(n_hvg=n_hvg, dims_use=dims_use, k_nn=k_nn,
                 res_main=res_main, min_cells=min_cells,
                 regress_vars=regress_vars)
    if (!is.null(per_class_params) && class_name %in% names(per_class_params)) {
      over <- per_class_params[[class_name]]
      for (nm in names(over)) base[[nm]] <- over[[nm]]
    }
    base
  }
  
  # 子对象细化（仅在子对象 meta 计算 percent.mt；表达流程在子对象里跑）
  .refine_one <- function(obj_sub, class_name, params){
    DefaultAssay(obj_sub) <- assay
    
    # 仅当需要回归时计算 percent.mt（写入子对象 meta）
    if ("percent.mt" %in% params$regress_vars && !"percent.mt" %in% colnames(obj_sub@meta.data)) {
      pct <- Seurat::PercentageFeatureSet(obj_sub, pattern="^MT-", assay=assay)
      obj_sub@meta.data$percent.mt <- as.numeric(pct)
      rownames(obj_sub@meta.data) <- colnames(obj_sub) # 强制对齐
    }
    
    obj_sub <- NormalizeData(obj_sub, normalization.method="LogNormalize",
                             scale.factor=1e4, verbose=verbose)
    obj_sub <- FindVariableFeatures(obj_sub, nfeatures = params$n_hvg,
                                    selection.method = "vst", verbose = verbose)
    if (exists(".drop_tech_genes", mode="function"))
      VariableFeatures(obj_sub) <- .drop_tech_genes(VariableFeatures(obj_sub))
    
    obj_sub <- ScaleData(obj_sub, features = rownames(obj_sub),
                         vars.to.regress = params$regress_vars, verbose = verbose)
    
    npcs <- max(params$dims_use)
    obj_sub <- RunPCA(obj_sub, features = VariableFeatures(obj_sub),
                      npcs = npcs, verbose = verbose)
    obj_sub <- FindNeighbors(obj_sub, dims = params$dims_use,
                             k.param = params$k_nn, verbose = verbose)
    
    chosen_res <- params$res_main
    if (isTRUE(do_resolution_sweep)) {
      res_vec <- unique(sort(res_grid))
      scores  <- numeric(length(res_vec))
      coords  <- Embeddings(obj_sub, "pca")[, params$dims_use, drop=FALSE]
      for (i in seq_along(res_vec)) {
        obj_tmp <- if (is.null(cluster_algorithm)) {
          FindClusters(obj_sub, resolution = res_vec[i], verbose = verbose)
        } else {
          FindClusters(obj_sub, resolution = res_vec[i],
                       algorithm = cluster_algorithm, verbose = verbose)
        }
        if (exists(".neighbor_consistency", mode="function")) {
          scores[i] <- .neighbor_consistency(Idents(obj_tmp), coords, k = 15)
        } else scores[i] <- 0
      }
      if (sum(scores) != 0) chosen_res <- res_vec[which.max(scores)]
      if (verbose) message("[", class_name, "] best resolution=", chosen_res)
    }
    
    if (is.null(cluster_algorithm)) {
      obj_sub <- FindClusters(obj_sub, resolution = chosen_res, verbose = verbose)
    } else {
      obj_sub <- FindClusters(obj_sub, resolution = chosen_res,
                              algorithm = cluster_algorithm, verbose = verbose)
    }
    
    # 给该类命名一个 UMAP reduction
    if (run_umap) {
      red_name <- paste0("umap_", .sanitize(class_name))
      obj_sub <- RunUMAP(obj_sub, dims = params$dims_use,
                         reduction.name = red_name, verbose = verbose)
    } else {
      red_name <- NULL
    }
    list(obj = obj_sub, best_res = chosen_res, umap_name = red_name)
  }
  
  # ---------- 关键修复：字符写入，最后再 factor 化 ----------
  object <- .ensure_meta_cols_if_missing(object, "Seurat_L2")
  
  # keep_existing：FALSE 时清空；TRUE 时保留旧值
  if (!keep_existing) object@meta.data$Seurat_L2 <- NA_character_
  
  # 若当前 L2 是 factor，先转字符，避免写入时“invalid factor level, NA generated”
  if (is.factor(object@meta.data$Seurat_L2)) {
    object@meta.data$Seurat_L2 <- as.character(object@meta.data$Seurat_L2)
  }
  old_L2_chr <- as.character(object@meta.data$Seurat_L2)  # 备份（字符）
  
  per_class_meta <- list()
  subobjects <- list()
  
  for (tg in targets) {
    # 精确匹配或以 "<tg>_" 开头
    sel <- which(coarse_vec == tg | startsWith(coarse_vec, paste0(tg, "_")))
    if (!length(sel)) next
    
    params <- .get_params_for(tg)
    if (length(sel) < params$min_cells) {
      if (verbose) message("[", tg, "] skipped: n=", length(sel), " < min_cells")
      next
    }
    
    cells_use <- colnames(object)[sel]
    
    # —— 用 counts + meta 重建“最小子对象”（更稳） —— #
    counts_all <- tryCatch(
      GetAssayData(object, assay = assay, layer = "counts"),
      error = function(e) GetAssayData(object, assay = assay, slot = "counts")
    )
    counts_sub <- counts_all[, cells_use, drop = FALSE]
    meta_sub   <- object@meta.data[cells_use, , drop = FALSE]
    stopifnot(identical(colnames(counts_sub), rownames(meta_sub)))
    
    obj_sub <- CreateSeuratObject(counts = counts_sub, assay = assay, meta.data = meta_sub,
                                  min.cells = 0,      # <-- 新增：不因基因过滤细胞
                                  min.features = 0    # <-- 新增：不因细胞表达量过滤细胞
                                  )
    res <- .refine_one(obj_sub, tg, params)
    obj_sub <- res$obj
    
    # —— 回填前：子簇列若为 factor，转字符 —— #
    lab_col <- paste0("Seurat_subcluster_", .sanitize(tg))
    object <- .ensure_meta_cols_if_missing(object, lab_col)
    if (is.factor(object@meta.data[[lab_col]])) {
      object@meta.data[[lab_col]] <- as.character(object@meta.data[[lab_col]])
    }
    
    # 回填（全用字符写入，避免 NA）
    object@meta.data[cells_use, lab_col] <- as.character(Idents(obj_sub))
    object@meta.data[cells_use, "Seurat_L2"] <-
      paste0(tg, "_", object@meta.data[cells_use, lab_col])
    

    if (run_umap && !is.null(res$umap_name)) {
      umap_name <- res$umap_name
      emb <- Embeddings(obj_sub, umap_name)  # 行名 = cells_use
      
      # 1) 保存为 reduction（允许只包含该类细胞）
      if (isTRUE(save_umap_reduction)) {
        dr <- CreateDimReducObject(
          embeddings = emb,
          key       = paste0("UMAP", substr(.sanitize(tg),1,3), "_"),
          assay     = assay
        )
        object[[umap_name]] <- dr  # 覆盖/复跑可直接更新
      }
      
      # 2) 同步写入 meta（便于导出/快速绘图）
      if (isTRUE(save_umap_meta)) {
        col1 <- paste0(umap_meta_prefix, "_", .sanitize(tg), "_1")
        col2 <- paste0(umap_meta_prefix, "_", .sanitize(tg), "_2")
        object <- .ensure_meta_cols_if_missing(object, c(col1, col2))
        object@meta.data[cells_use, col1] <- emb[,1]
        object@meta.data[cells_use, col2] <- emb[,2]
      }
    }
    
    per_class_meta[[tg]] <- list(
      best_resolution = res$best_res,
      n_cells = length(cells_use),
      n_clusters = length(table(Idents(obj_sub))),
      params_used = params,
      umap_saved = if (run_umap) paste0("umap_", .sanitize(tg)) else NA
    )
    if (isTRUE(return_subobjects)) subobjects[[tg]] <- obj_sub
  }
  
  # —— 循环结束后：统一 factor 化（合并新旧 level） —— #
  new_L2_chr <- as.character(object@meta.data$Seurat_L2)
  lvls <- unique(na.omit(c(old_L2_chr, new_L2_chr)))
  object@meta.data$Seurat_L2 <- factor(new_L2_chr, levels = lvls)
  
  object@misc$seurat_l2 <- list(
    coarse_col = coarse_col,
    classes = per_class_meta,
    options = list(run_umap=run_umap, sweep=do_resolution_sweep,
                   res_grid=res_grid, cluster_algorithm=cluster_algorithm,
                   keep_existing=keep_existing,
                   save_umap_reduction=save_umap_reduction,
                   save_umap_meta=save_umap_meta,
                   umap_meta_prefix=umap_meta_prefix)
  )
  
  message("[L2] finished. Updated: Seurat_L2 (incremental), per-class Seurat_subcluster_*, and class-specific UMAPs")
  if (isTRUE(return_subobjects)) return(list(object=object, subobjects=subobjects))
  object
}






subset_and_export <- function(obj, keep_expr, name, logfc.threshold = NULL){
  cells <- WhichCells(obj, expression = !!rlang::parse_expr(keep_expr))
  if (length(cells) < 20) { message("[skip] too few cells for ", name); return(invisible(NULL)) }
  counts_sub <- GetAssayData(obj, assay="Spatial", layer="counts")[, cells, drop=FALSE]
  meta_sub   <- obj@meta.data[cells, , drop=FALSE]
  sobj <- CreateSeuratObject(counts = counts_sub, assay = "Spatial", meta.data = meta_sub)
  sobj <- ensure_log_data(sobj, assay = "Spatial")
  saveRDS(sobj, file.path(Config$Bankspath, paste0(name, "_object.rds")))
  # 导出 markers（若有 L2）
  if ("Seurat_L2" %in% colnames(sobj@meta.data)) {
    ExportTopMarkers(
      sobj,
      out_xlsx = file.path(Config$Bankspath, paste0(name, "_L2_cluster_top30_markers.xlsx")),
      group_col = "Seurat_L2",
      assay = "Spatial",
      n_top = 40,
      only.pos = TRUE,
      logfc.threshold = logfc.threshold %||% 0
    )
  }
  invisible(sobj)
}





# =========================
# Rank Top-N markers per cluster by a composite score that favors:
# - higher avg_log2FC (intensity)
# - higher (pct.1 - pct.2) (prevalence delta)
# - lower p_val_adj (higher -log10 p)
#
# New args:
#   rank_weights: numeric(3) weights for (fc, delta, p), default c(0.4,0.4,0.2)
#   require_delta_positive: if TRUE, keep genes with (pct.1 - pct.2) > min_delta_pct
#   min_delta_pct: minimal delta (pct.1 - pct.2) required when require_delta_positive=TRUE
#   n_top: default 50 (as requested)
#
# Output:
#   - Single-sheet xlsx with added columns: delta_pct, r_fc, r_delta, r_p, combo_score
#   - Returns the TopN tibble (invisible)
ExportTopMarkers <- function(
    object,
    out_xlsx,
    group_col = "banksy_celltype",
    assay = DefaultAssay(object),
    n_top = 60,
    only.pos = TRUE,
    min.pct = 0.05,
    logfc.threshold = 0.10,
    test.use = "wilcox",
    return.thresh = 0.25,
    normalize_if_missing = TRUE,
    rank_weights = c(fc = 0.4, delta = 0.4, p = 0.2),
    require_delta_positive = TRUE,
    min_delta_pct = 0
){
  stopifnot(inherits(object, "Seurat"))
  if (!group_col %in% colnames(object@meta.data)) {
    stop("`group_col` '", group_col, "' not found in object@meta.data.")
  }
  suppressPackageStartupMessages({ require(Seurat); require(dplyr) })
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Please install 'openxlsx' (install.packages('openxlsx')).")
  }
  
  # --- helper: ensure log-normalized data exists ---
  ensure_log_data <- function(obj, assay){
    # If 'data' slot is empty or all zeros, run NormalizeData quickly.
    mat <- tryCatch(GetAssayData(obj, assay = assay, slot = "data"),
                    error = function(e) NULL)
    need_norm <- is.null(mat) || nrow(mat) == 0 || sum(mat) == 0
    if (isTRUE(need_norm)) {
      message("[info] 'data' slot empty for assay '", assay, "'. Running NormalizeData()...")
      obj <- NormalizeData(obj, assay = assay, normalization.method = "LogNormalize")
    }
    obj
  }
  
  # Use chosen assay
  DefaultAssay(object) <- assay
  
  # Normalize if needed
  if (isTRUE(normalize_if_missing)) {
    object <- ensure_log_data(object, assay = assay)
  }
  
  # Idents
  Idents(object) <- object[[group_col]][, 1]
  
  # DE
  markers <- Seurat::FindAllMarkers(
    object = object, assay = assay,
    only.pos = only.pos,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    test.use = test.use,
    return.thresh = return.thresh
  )
  
  # Basic checks
  need_cols <- c("avg_log2FC","p_val_adj","pct.1","pct.2","cluster")
  missing_cols <- setdiff(need_cols, colnames(markers))
  if (length(missing_cols)) {
    stop("Missing columns in FindAllMarkers() output: ", paste(missing_cols, collapse = ", "))
  }
  if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
  
  # Compute delta pct
  markers <- markers %>%
    mutate(delta_pct = pct.1 - pct.2)
  
  # Optional filter on delta_pct
  if (isTRUE(require_delta_positive)) {
    markers <- markers %>% filter(delta_pct > min_delta_pct)
  }
  
  # Normalize rank weights
  if (is.null(names(rank_weights)) || !all(c("fc","delta","p") %in% names(rank_weights))) {
    names(rank_weights) <- c("fc","delta","p")
  }
  rank_weights <- rank_weights / sum(rank_weights)
  
  # --- robust ranking without using n() in mutate ---
  eps <- 1e-300
  safe_percent_rank <- function(x){
    # Handle edge cases: length 0/1 or all NA -> return 0.5 constants
    if (length(x) <= 1 || all(is.na(x))) return(rep(0.5, length(x)))
    dplyr::percent_rank(x)
  }
  
  markers <- markers %>%
    group_by(cluster) %>%
    mutate(
      r_fc    = safe_percent_rank(avg_log2FC),
      r_delta = safe_percent_rank(delta_pct),
      r_p     = safe_percent_rank(-log10(p_val_adj + eps)),
      combo_score = r_fc * rank_weights["fc"] +
        r_delta * rank_weights["delta"] +
        r_p * rank_weights["p"]
    ) %>%
    ungroup() %>%
    arrange(cluster,
            desc(combo_score),
            desc(avg_log2FC),
            desc(delta_pct),
            p_val_adj)
  
  # Take Top-N per cluster (guard for small groups)
  topn <- markers %>%
    group_by(cluster) %>%
    slice_head(n = n_top) %>%
    ungroup()
  
  # Write xlsx
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "TopN_Scored")
  openxlsx::writeData(wb, "TopN_Scored", topn)
  openxlsx::addFilter(wb, "TopN_Scored", 1, 1)
  openxlsx::freezePane(wb, "TopN_Scored", firstActiveRow = 2)
  openxlsx::setColWidths(wb, "TopN_Scored", cols = 1:ncol(topn), widths = "auto")
  openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
  
  message("[OK] Wrote scored top-", n_top, " markers per cluster to: ", out_xlsx)
  invisible(topn)
}






suppressPackageStartupMessages({
  library(Seurat); library(FNN); library(Matrix); library(dplyr)
})

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





suppressPackageStartupMessages({ library(Matrix); library(FNN); library(Seurat) })

knn_smooth <- function(
    object,
    meta_cols = NULL,          # e.g., c("Score_ECM")
    genes     = NULL,          # e.g., c("ITGAV","ITGB1","ITGB5")
    k         = 8,
    method    = "gaussian",    # "mean" or "gaussian"
    sigma     = NULL,          # if gaussian, auto from 1-NN median if NULL
    assay     = NULL,          # if NULL, choose prefer if exists
    slot      = "data",        # typically "data"
    prefer    = "Spatial",     # preferred assay name
    out_suffix = NULL          # default "_knn{k}" or "_wknn{k}"
){
  stopifnot(inherits(object, "Seurat"))
  # ---- choose assay/slot ----
  if (is.null(assay)) {
    assay <- if (prefer %in% names(object@assays)) prefer else DefaultAssay(object)
  }
  if (is.null(slot)) slot <- "data"
  DefaultAssay(object) <- assay
  M <- GetAssayData(object, assay = assay, slot = slot)  # genes x cells (dgCMatrix)
  
  # ---- coords: prefer scaled; fallback to raw ----
  md <- object@meta.data
  coord <- if (all(c("imagecol_scaled","imagerow_scaled") %in% colnames(md))) {
    as.matrix(md[, c("imagecol_scaled","imagerow_scaled")])
  } else if (all(c("imagecol","imagerow") %in% colnames(md))) {
    as.matrix(md[, c("imagecol","imagerow")])
  } else stop("Cannot find spatial coordinate columns.")
  
  mode <- tolower(method[1]); if (!mode %in% c("mean","gaussian")) mode <- "gaussian"
  
  n <- nrow(coord)
  k_eff <- max(1, min(k, n - 1))  # safety
  nn <- FNN::get.knn(coord, k = k_eff)
  
  # ---- build weights per target cell (row-wise) ----
  # for each i, neighbors are nn.index[i,], distances nn.dist[i,]
  # include self (j=i, d=0). We will build a sparse W so that v_smooth = W %*% v.
  dist_mat <- cbind(0, nn$nn.dist)           # n x (k+1)
  if (mode == "gaussian") {
    if (is.null(sigma)) {
      sig <- stats::median(nn$nn.dist[,1], na.rm = TRUE)
      if (!is.finite(sig) || sig <= 0) sig <- 1
      sigma <- sig
    }
    W_row <- exp(-(dist_mat^2) / (2 * sigma^2))
  } else {
    W_row <- matrix(1, nrow = n, ncol = k_eff + 1)
  }
  # normalize rows to sum 1
  rs <- rowSums(W_row)
  rs[!is.finite(rs) | rs <= 0] <- 1
  W_row <- W_row / rs
  
  # ---- assemble sparse W (n x n), row i has weights on {i ∪ N(i)}
  idx_rows <- as.vector( row( W_row ) )       # 1..n repeated (k+1) times
  # build corresponding columns [self, neighbors]
  self_col <- matrix(rep(seq_len(n), each = 1), nrow = n, ncol = 1)
  neigh_col <- nn$nn.index                    # n x k
  idx_cols <- as.vector( cbind(self_col, neigh_col) )
  wvals    <- as.vector( W_row )
  
  W <- sparseMatrix(i = idx_rows, j = idx_cols, x = wvals, dims = c(n, n))  # row-stochastic
  
  # ---- suffix ----
  if (is.null(out_suffix)) {
    out_suffix <- if (mode == "mean") paste0("_knn", k_eff) else paste0("_wknn", k_eff)
  }
  
  # ---- smooth meta columns (numeric) ----
  if (!is.null(meta_cols) && length(meta_cols)) {
    for (nm in meta_cols) {
      if (!nm %in% colnames(md)) { message("[skip] meta col not found: ", nm); next }
      v <- as.numeric(md[[nm]])
      if (all(!is.finite(v))) { message("[skip] meta col non-numeric: ", nm); next }
      v[!is.finite(v)] <- 0
      v_sm <- as.numeric(W %*% v)
      object[[paste0(nm, out_suffix)]] <- v_sm
    }
  }
  
  # ---- smooth selected genes from assay/slot ----
  if (!is.null(genes) && length(genes)) {
    genes2 <- intersect(genes, rownames(M))
    miss <- setdiff(genes, genes2)
    if (length(miss)) message("[WARN] genes not in ", assay, "/", slot, ": ", paste(miss, collapse = ", "))
    if (length(genes2)) {
      # extract submatrix (genes x cells) and smooth columns via W
      X <- M[genes2, , drop = FALSE]                 # dgCMatrix
      X_sm <- X %*% W                                # genes x cells
      # write back as meta columns (one per gene) for PlotExpressionV2
      for (i in seq_along(genes2)) {
        object[[paste0(genes2[i], out_suffix)]] <- as.numeric(X_sm[i, ])
      }
    }
  }
  
  attr(object, "knn_smooth_source") <- list(
    assay = assay, slot = slot, method = mode, k = k_eff, sigma = sigma
  )
  object
}






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

# ==== RCTD/Cellbin ?? ====








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

#
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





RunRCTD_Pass1 <- function(
    object,                          # 已创建好的 Cellbin Seurat 对象
    ref_rds,                         # 参考 SCE 的 RDS 文件（含 counts、cell_type_major）
    out_dir,                         # 输出文件夹
    assay              = "Spatial",
    umi_min            = 200,        # Cellbin 建议 150-300
    cell_min_instance  = 20,         # RCTD 参数：每类至少多少细胞
    n_cores            = 8,          # 并行核数（spacexr用到）
    use_scaled_coords  = TRUE,       # 坐标优先用 imagecol_scaled/imagerow_scaled
    assign_global      = TRUE,       # 是否把 rctd1 放入 .GlobalEnv
    save_seurat        = TRUE,       # 是否保存带RCTD结果的Seurat对象
    export_loupe_col   = "RCTD_first", # 有 export_loupe_categories() 时导出此列
    min_intersect_genes= 500,        # 与参考相交基因的最低阈值
    verbose            = TRUE
){
  stopifnot(inherits(object, "Seurat"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # —— 屏蔽 Matrix 的“弃用强制转换”告警，并在退出时还原 ——
  old_opt <- getOption("Matrix.warnDeprecatedCoerce")
  options(Matrix.warnDeprecatedCoerce = FALSE)
  on.exit(options(Matrix.warnDeprecatedCoerce = old_opt), add = TRUE)
  
  .msg <- function(...) if (verbose) message(sprintf(...))
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(Matrix)
    require(SummarizedExperiment)
    require(SingleCellExperiment)
    require(spacexr)
    require(dplyr)
  })
  options(mc.cores = n_cores)
  
  # ---- helpers ----
  as_Csparse_dgc <- function(mat) {
    if (!inherits(mat, "CsparseMatrix")) mat <- methods::as(mat, "CsparseMatrix")
    if (!inherits(mat, "dgCMatrix"))     mat <- methods::as(mat, "dgCMatrix")
    mat
  }
  
  ensure_integer_counts <- function(mat){
    mat <- as_Csparse_dgc(mat)
    if (!is.integer(mat@x)) mat@x <- round(mat@x)
    mat
  }
  
  get_cellbin_counts_coords <- function(obj, assay = "Spatial", umi_min = 100, use_scaled = TRUE){
    counts <- tryCatch(
      Seurat::GetAssayData(obj, assay = assay, layer = "counts"),
      error = function(e) Seurat::GetAssayData(obj, assay = assay, slot = "counts")
    )
    counts <- as_Csparse_dgc(counts)  # ← 修正点①：统一到 Csparse/dgC
    
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
      if (length(Seurat::Images(obj)) < 1) stop("Seurat对象里找不到影像坐标。")
      img <- obj@images[[ Seurat::Images(obj)[1] ]]
      C <- img@coordinates
      cand_pairs_img <- list(c("pxl_col_in_fullres","pxl_row_in_fullres"),
                             c("imagecol","imagerow"),
                             c("x","y"), c("X","Y"))
      for (p in cand_pairs_img) if (all(p %in% colnames(C))) { pick_cols <- p; break }
      if (is.null(pick_cols)) stop("在 image@coordinates 中也找不到坐标列。")
      coords <- C[, pick_cols, drop = FALSE]
      colnames(coords) <- c("x","y")
    }
    
    common <- intersect(colnames(counts), rownames(coords))
    if (length(common) < 100) stop("counts 与 coords 交集太小：", length(common))
    counts <- counts[, common, drop = FALSE]
    coords <- coords[common, , drop = FALSE]
    
    nUMI <- Matrix::colSums(counts)
    keep <- nUMI >= umi_min
    counts <- counts[, keep, drop = FALSE]
    coords <- coords[keep, , drop = FALSE]
    
    coords <- data.frame(
      x = as.numeric(coords[,1]),
      y = as.numeric(coords[,2]),
      row.names = rownames(coords),
      check.names = FALSE
    )
    list(counts = counts, coords = coords)
  }
  
  extract_rctd_meta <- function(rctd_obj){
    res <- rctd_obj@results
    if (!is.null(res$results_df) && is.data.frame(res$results_df) && nrow(res$results_df) > 0) {
      df <- as.data.frame(res$results_df)
    } else if (is.data.frame(res) && nrow(res) > 0) {
      df <- as.data.frame(res)
    } else {
      stop("RCTD results 为空；run.RCTD 可能失败。")
    }
    rn <- rownames(df)
    if (is.null(rn) || length(rn) == 0) stop("results_df 没有行名（cell_id）。")
    
    spot_class <- if ("spot_class" %in% names(df)) as.character(df$spot_class) else
      if ("class" %in% names(df)) as.character(df$class) else
        if (!is.null(res$spot_class) && length(res$spot_class) == nrow(df)) as.character(res$spot_class) else
          rep(NA_character_, nrow(df))
    
    first_type  <- if ("first_type"  %in% names(df)) as.character(df$first_type)  else
      if ("cell_type"  %in% names(df)) as.character(df$cell_type)   else NA_character_
    if (length(first_type)  == 1L) first_type  <- rep(first_type,  nrow(df))
    
    second_type <- if ("second_type" %in% names(df)) as.character(df$second_type) else NA_character_
    if (length(second_type) == 1L) second_type <- rep(second_type, nrow(df))
    
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
    return(NULL)
  }
  
  # ---- load reference ----
  .msg(">>> 加载参考：%s", ref_rds)
  sce_ref <- readRDS(ref_rds)
  if (!"counts" %in% SummarizedExperiment::assayNames(sce_ref)) {
    # 兜底：取第一个 assay
    SummarizedExperiment::assay(sce_ref, "counts") <- SummarizedExperiment::assay(sce_ref, SummarizedExperiment::assayNames(sce_ref)[1])
  }
  if (is.null(sce_ref$cell_type_major)) stop("参考中缺少列 'cell_type_major'。")
  
  # ---- extract counts & coords from Cellbin Seurat ----
  .msg(">>> 抽取 Cellbin counts/coords（umi_min=%d）", umi_min)
  cb <- get_cellbin_counts_coords(object, assay = assay, umi_min = umi_min, use_scaled = use_scaled_coords)
  counts_cb <- cb$counts; coords_cb <- cb$coords
  .msg("Cellbin after QC: genes=%d cells=%d", nrow(counts_cb), ncol(counts_cb))
  
  # ---- gene intersection & cleaning ----
  ref_counts <- SummarizedExperiment::assay(sce_ref, "counts")
  ref_counts <- as_Csparse_dgc(ref_counts)   # ← 修正点②：统一到 Csparse/dgC
  ref_counts <- ensure_integer_counts(ref_counts)
  
  genes_common <- intersect(rownames(counts_cb), rownames(ref_counts))
  if (length(genes_common) < min_intersect_genes) {
    stop(sprintf("参考与样本交集基因过少：%d < %d", length(genes_common), min_intersect_genes))
  }
  ref_counts <- ref_counts[genes_common, , drop = FALSE]
  counts_use <- counts_cb[genes_common, , drop = FALSE]
  counts_use <- ensure_integer_counts(counts_use)
  
  keep_g <- Matrix::rowSums(ref_counts) > 0 & Matrix::rowSums(counts_use) > 0
  ref_counts <- ref_counts[keep_g, , drop = FALSE]
  counts_use <- counts_use[keep_g, , drop = FALSE]
  
  ref_nonzero <- Matrix::colSums(ref_counts) > 0
  ref_counts  <- ref_counts[, ref_nonzero, drop = FALSE]
  ref_types   <- as.character(sce_ref$cell_type_major)[ref_nonzero]
  ok <- !is.na(ref_types)
  ref_counts <- ref_counts[, ok, drop = FALSE]
  ref_types  <- factor(ref_types[ok])
  names(ref_types) <- colnames(ref_counts)
  
  if (anyDuplicated(names(ref_types)) > 0) stop("参考条码重复，请唯一化参考列名。")
  stopifnot(identical(colnames(counts_use), rownames(coords_cb)))
  
  # ---- build RCTD and run ----
  .msg(">>> 构建 RCTD 对象（CELL_MIN_INSTANCE=%d）", cell_min_instance)
  ref_obj <- spacexr::Reference(ref_counts, ref_types)
  sp_obj  <- spacexr::SpatialRNA(
    counts = counts_use,
    coords = coords_cb,
    nUMI   = Matrix::colSums(counts_use)
  )
  
  rctd1 <- spacexr::create.RCTD(sp_obj, ref_obj, CELL_MIN_INSTANCE = cell_min_instance)
  .msg(">>> 运行 run.RCTD(doublet)…")
  rctd1 <- spacexr::run.RCTD(rctd1, doublet_mode = "doublet")
  
  # 保存 rctd1
  rctd_rds <- file.path(out_dir, "RCTD_pass1_object.rds")
  saveRDS(rctd1, rctd_rds)
  .msg("已保存 rctd1 到：%s", rctd_rds)
  if (assign_global) {
    assign("rctd1", rctd1, envir = .GlobalEnv)
    .msg("已把 rctd1 放入全局环境。")
  }
  
  # ---- extract meta & weights ----
  meta_df <- extract_rctd_meta(rctd1)
  write.csv(meta_df, file.path(out_dir, "RCTD_pass1_meta.csv"), row.names = FALSE)
  
  W <- extract_rctd_weights(rctd1)
  if (is.null(W)) {
    .msg("未从结果里拿到 weights，构造一个保守fallback（singlet=1/0；doublet=0.5/0.5）…")
    types <- sort(unique(c(meta_df$first_type, meta_df$second_type)))
    types <- types[!is.na(types)]
    W <- matrix(0, nrow = nrow(meta_df), ncol = length(types),
                dimnames = list(meta_df$cell_id, types))
    is_doublet <- grepl("doublet", meta_df$spot_class %||% "", ignore.case = TRUE)
    # singlet
    idx_sing <- which(!is_doublet & !is.na(meta_df$first_type))
    if (length(idx_sing) > 0) {
      W[cbind(meta_df$cell_id[idx_sing], meta_df$first_type[idx_sing])] <- 1
    }
    # doublet
    idx_dou <- which(is_doublet & !is.na(meta_df$first_type))
    if (length(idx_dou) > 0) {
      W[cbind(meta_df$cell_id[idx_dou], meta_df$first_type[idx_dou])] <- 0.5
      has_second <- idx_dou[!is.na(meta_df$second_type[idx_dou])]
      if (length(has_second) > 0) {
        W[cbind(meta_df$cell_id[has_second], meta_df$second_type[has_second])] <-
          W[cbind(meta_df$cell_id[has_second], meta_df$second_type[has_second])] + 0.5
      }
    }
  }
  saveRDS(W, file.path(out_dir, "RCTD_pass1_weights.rds"))
  
  # ---- write back to Seurat meta ----
  object$RCTD_spot_class <- NA_character_
  object$RCTD_first      <- NA_character_
  object$RCTD_second     <- NA_character_
  
  idx <- match(rownames(object@meta.data), meta_df$cell_id)
  keep <- !is.na(idx)
  object$RCTD_spot_class[keep] <- meta_df$spot_class[idx[keep]]
  object$RCTD_first[keep]      <- meta_df$first_type[idx[keep]]
  object$RCTD_second[keep]     <- meta_df$second_type[idx[keep]]
  
  # 归一化 W 并计算 top1/2 & margin
  cells_use <- intersect(rownames(object@meta.data), rownames(W))
  W_sub <- W[cells_use, , drop = FALSE]
  rs <- rowSums(W_sub); nz <- rs > 0
  if (any(nz)) W_sub[nz, ] <- W_sub[nz, , drop = FALSE] / rs[nz]
  
  get_top12 <- function(v){
    if (all(!is.finite(v))) return(c(NA, NA, NA_real_, NA_real_))
    ord <- order(v, decreasing = TRUE)
    t1  <- names(v)[ord[1]]; w1 <- v[ord[1]]
    if (length(v) >= 2) { t2 <- names(v)[ord[2]]; w2 <- v[ord[2]] } else { t2 <- NA; w2 <- 0 }
    c(t1, t2, w1, w2)
  }
  topmat <- t(apply(W_sub, 1, get_top12))
  colnames(topmat) <- c("first_calc","second_calc","first_w","second_w")
  
  object$RCTD_first_calc  <- NA_character_
  object$RCTD_second_calc <- NA_character_
  object$RCTD_first_w     <- NA_real_
  object$RCTD_second_w    <- NA_real_
  object$RCTD_margin      <- NA_real_
  
  m <- match(rownames(object@meta.data), rownames(topmat))
  keep2 <- !is.na(m)
  object$RCTD_first_calc [keep2] <- topmat[m[keep2], "first_calc"]
  object$RCTD_second_calc[keep2] <- topmat[m[keep2], "second_calc"]
  object$RCTD_first_w    [keep2] <- as.numeric(topmat[m[keep2], "first_w"])
  object$RCTD_second_w   [keep2] <- as.numeric(topmat[m[keep2], "second_w"])
  object$RCTD_margin     [keep2] <- object$RCTD_first_w[keep2] - object$RCTD_second_w[keep2]
  
  if (!is.null(export_loupe_col) && exists("export_loupe_categories", mode = "function")) {
    try({
      export_loupe_categories(
        object,
        meta_col = export_loupe_col,
        out_dir  = out_dir,
        filename = NULL,
        na_label = "Unlabeled",
        gzip     = FALSE,
        overwrite= TRUE,
        verbose  = TRUE
      )
    }, silent = !verbose)
  }
  
  if (save_seurat) {
    sid <- object@project.name %||% "Sample"
    saveRDS(object, file.path(out_dir, sprintf("%s_Seurat_with_RCTD.rds", sid)))
    .msg("已保存带 RCTD 结果的 Seurat：%s", file.path(out_dir, sprintf("%s_Seurat_with_RCTD.rds", sid)))
  }
  
  .msg("表：RCTD_first")
  .msg(capture.output(print(table(object$RCTD_first, useNA="ifany"))))
  .msg("表：RCTD_spot_class")
  .msg(capture.output(print(table(object$RCTD_spot_class, useNA="ifany"))))
  .msg("RCTD_margin 摘要：")
  .msg(capture.output(print(summary(object$RCTD_margin))))
  
  invisible(list(
    seurat  = object,
    meta    = meta_df,
    weights = W,
    rctd_rds= rctd_rds
  ))
}






suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(FNN)
})

run_cellbin_ffpe <- function(
    obj,
    assay = "Spatial",
    subset = NULL, cells_use = NULL,
    n_hvg = 3000,
    # --- Smoothing ---
    smooth_mode = c("data","pca"),
    k_smooth = 6,
    bilateral_markers = NULL,        # 显式护边基因；NULL=不护边
    max_bilateral_markers = 64,      # 并行/串行都生效（你可拉大到 150/200）
    # --- Graphs & clustering ---
    dims_use = 1:20,
    k_expr = 20,
    k_sp = 6,
    alpha = 0.8,
    resolution = 1.0,
    do_umap = TRUE,
    # --- Extras ---
    sig_list = NULL,
    seed = 123, verbose = TRUE,
    progress_every = 100,
    save_smoothed_to = NULL,
    # --- Back-compat convenience: bilateral=TRUE -> capped HVGs
    bilateral = NULL,
    # --- NEW: 并行参数 ---
    parallel = FALSE,                # TRUE：并行双边护边平滑
    workers  = NULL,                 # NULL=自动(cores-1)，或指定整数
    chunk_size = 32,                 # 按“基因块”并行，单任务最多多少个基因
    manage_plan = TRUE,              # TRUE：在函数内临时设 plan()
    set_blas_threads = FALSE,        # TRUE：把 OMP/MKL/OPENBLAS 线程数限为1（防过订阅）
    force_uwot = TRUE                # TRUE：强制 R 版 UMAP，绕开 Python/reticulate
){
  set.seed(seed)
  smooth_mode <- match.arg(smooth_mode)
  
  # back-compat: bilateral=TRUE -> 用 capped 顶级 HVG 护边
  if (!is.null(bilateral) && isTRUE(bilateral)) {
    warning("[deprecated] 'bilateral=TRUE' detected; will edge-preserve a CAPPED set of top HVGs.")
    bilateral_markers <- ".__USE_TOP_HVGS__"
  }
  
  # 限制底层 BLAS/OpenMP 线程，避免与 R 层并行冲突
  if (isTRUE(parallel) && isTRUE(set_blas_threads)) {
    Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1)
    }
  }
  
  # helpers ----
  ensure_spatial_logdata <- function(o, assay="Spatial"){
    has <- tryCatch({
      M <- GetAssayData(o, assay=assay, layer="data")
      inherits(M,"dgCMatrix") && ncol(M)>0
    }, error=function(e) FALSE)
    if (has) return(o)
    counts <- tryCatch(GetAssayData(o, assay=assay, layer="counts"),
                       error=function(e) NULL)
    if (is.null(counts) || ncol(counts)==0)
      counts <- tryCatch(GetAssayData(o, assay=assay, slot="counts"),
                         error=function(e) NULL)
    if (is.null(counts) || ncol(counts)==0) stop("No 'counts' to normalize.")
    lib  <- pmax(Matrix::colSums(counts), 1)
    norm <- log1p(t(t(counts) * (1e4/lib)))
    ok <- TRUE
    o2 <- tryCatch(SetAssayData(o, assay=assay, layer="data", new.data=norm),
                   error=function(e){ ok <<- FALSE; NULL })
    if (!ok) o2 <- SetAssayData(o, assay=assay, slot="data", new.data=norm)
    o2
  }
  get_xy <- function(o){
    md <- o@meta.data
    if (all(c("imagecol_scaled","imagerow_scaled") %in% names(md))) {
      cbind(x=as.numeric(md$imagecol_scaled), y=as.numeric(md$imagerow_scaled))
    } else if (all(c("imagecol","imagerow") %in% names(md))) {
      cbind(x=as.numeric(md$imagecol), y=as.numeric(md$imagerow))
    } else stop("No spatial coordinates found: need imagecol*/imagerow*.")
  }
  get_hvgs_ffpe <- function(o, n=3000, assay="Spatial"){
    if (requireNamespace("scran", quietly=TRUE) &&
        requireNamespace("SingleCellExperiment", quietly=TRUE)) {
      if (verbose) message("  - HVG via scran::modelGeneVar …")
      sce <- as.SingleCellExperiment(o, assay=assay)
      fit <- scran::modelGeneVar(sce)
      return(scran::getTopHVGs(fit, n=n))
    }
    if (verbose) message("  - HVG via Seurat::FindVariableFeatures(vst) …")
    o <- FindVariableFeatures(o, assay=assay, selection.method="vst", nfeatures=n, verbose=FALSE)
    VariableFeatures(o)
  }
  set_data_layer <- function(o, assay, mat){
    ok <- TRUE
    o2 <- tryCatch(SetAssayData(o, assay=assay, layer="data", new.data=mat),
                   error=function(e){ ok <<- FALSE; NULL })
    if (!ok) o2 <- SetAssayData(o, assay=assay, slot="data", new.data=mat)
    o2
  }
  ensure_dimnames_square <- function(M, cells){
    if (is.null(rownames(M)) || is.null(colnames(M))) dimnames(M) <- list(cells, cells)
    M
  }
  put_custom_graph <- function(o, W){
    W@x[!is.finite(W@x)] <- 0
    W   <- (W + Matrix::t(W))/2
    diag(W) <- 0
    W   <- Matrix::drop0(W)
    if (requireNamespace("SeuratObject", quietly=TRUE)) {
      g <- SeuratObject::as.Graph(W)
      o@graphs[["custom_graph"]] <- g
    } else {
      o@graphs$custom_graph <- as(W, "dgCMatrix")
    }
    o
  }
  
  knn_smooth_data <- function(o, k=6, hvgs=NULL, bilateral_markers=NULL,
                              assay="Spatial", layer="data",
                              verbose=TRUE, progress_every=100,
                              max_bilat = 64,
                              parallel=FALSE, workers=NULL, chunk_size=32, manage_plan=TRUE){
    if (verbose) message(sprintf("--- KNN smoothing (data): k=%d ---", k))
    E  <- GetAssayData(o, assay=assay, layer=layer)  # genes x cells (sparse ok)
    xy <- get_xy(o); n <- ncol(E)
    # KNN 一次
    nn <- FNN::get.knn(xy, k=k)
    J  <- nn$nn.index         # n x k
    D  <- nn$nn.dist          # n x k
    
    # 空间高斯核
    d1 <- D[,1]
    sigma <- stats::median(d1[d1>0], na.rm=TRUE); if (!is.finite(sigma) || sigma<=0) sigma <- 1
    Wsp_block <- exp(-(D^2)/(2*sigma^2))   # n x k
    
    # 行归一 Wfull（给 HVG 的高斯平滑）
    i_idx <- rep.int(seq_len(n), times=k)
    j_idx <- as.vector(t(J))
    w_nb  <- as.vector(t(Wsp_block))
    Wfull <- Matrix::sparseMatrix(i=c(i_idx, seq_len(n)),
                                  j=c(j_idx, seq_len(n)),
                                  x=c(w_nb, rep(1,n)), dims=c(n,n))
    rs <- Matrix::rowSums(Wfull); rs[!is.finite(rs) | rs==0] <- 1
    Wfull  <- Matrix::Diagonal(x=1/rs) %*% Wfull
    
    # HVG 高斯（一次稀疏乘就够了）
    gsel <- if (!is.null(hvgs)) intersect(hvgs, rownames(E)) else rownames(E)
    E_sm <- E
    if (length(gsel)) {
      if (verbose) message(sprintf("  [genes] Gaussian smoothing on %d HVGs", length(gsel)))
      E_sm[gsel, ] <- E[gsel, ] %*% Matrix::t(Wfull)
    }
    
    # 处理护边基因列表
    if (identical(bilateral_markers, ".__USE_TOP_HVGS__")) {
      markers <- head(gsel, max_bilat)
    } else if (is.null(bilateral_markers)) {
      markers <- character(0)
    } else {
      mk <- intersect(bilateral_markers, rownames(E))
      if (length(mk) > max_bilat) mk <- head(mk, max_bilat)
      markers <- mk
    }
    
    if (!length(markers)) {
      return(set_data_layer(o, assay=assay, mat=E_sm))
    }
    
    # 每基因的稳健尺度（MAD）
    if (verbose) message(sprintf("  [markers] Bilateral smoothing on %d marker(s)%s",
                                 length(markers),
                                 if (max_bilat<Inf) sprintf(" (cap=%d)", max_bilat) else ""))
    sgs <- vapply(markers, function(g){
      v <- as.numeric(E[g, ])
      sg <- stats::mad(v, center = stats::median(v), constant = 1.4826, na.rm = TRUE)
      if (!is.finite(sg) || sg <= 0) {
        sg <- max(1e-6, stats::sd(v, na.rm=TRUE)); if (!is.finite(sg) || sg<=0) sg <- 1
      }
      sg
    }, numeric(1))
    names(sgs) <- markers
    
    # —— 串行路径（稳）——
    bilateral_seq <- function(genes){
      res <- matrix(NA_real_, nrow=length(genes), ncol=n,
                    dimnames=list(genes, colnames(E)))
      t0 <- Sys.time()
      for (ii in seq_along(genes)) {
        g  <- genes[ii]; v <- as.numeric(E[g, ]); sg <- sgs[g]
        Vnb <- v[J]; dim(Vnb) <- dim(J)      # n x k
        Delta <- Vnb - v                     # 按行广播
        Sim <- exp(-(Delta*Delta)/(2*sg*sg))
        Wbl <- Sim * Wsp_block
        rsb <- rowSums(Wbl); rsb[!is.finite(rsb) | rsb==0] <- 1
        Wbl <- Wbl / rsb
        res[ii, ] <- rowSums(Wbl * Vnb)
        if (ii==1 || ii%%progress_every==0 || ii==length(genes)) {
          el <- as.numeric(difftime(Sys.time(), t0, units="secs"))
          rate <- ii / max(el,1e-9); eta <- (length(genes)-ii)/max(rate,1e-9)
          if (verbose) message(sprintf("    %d/%d (%.1f%%) | elapsed %.1fs | ETA %.1fs",
                                       ii, length(genes), 100*ii/length(genes), el, eta))
        }
      }
      res
    }
    
    # —— 并行路径（快）——
    bilateral_par <- function(genes, workers, chunk_size){
      if (!requireNamespace("future.apply", quietly = TRUE)) {
        warning("[parallel=TRUE] future.apply not installed; falling back to sequential.")
        return(bilateral_seq(genes))
      }
      # 切分成块，避免一次性返回超大列表
      chunks <- split(genes, ceiling(seq_along(genes)/max(1L, chunk_size)))
      if (is.null(workers)) {
        workers <- max(1L, parallel::detectCores() - 1L)
      }
      # 临时设置 plan
      if (isTRUE(manage_plan)) {
        if (.Platform$OS.type == "windows") {
          future::plan(future::multisession, workers=workers)
        } else {
          # RStudio/macOS 下也建议 multisession，防止 fork 与外部库冲突
          future::plan(future::multisession, workers=workers)
        }
        on.exit(try(future::plan(future::sequential), silent=TRUE), add=TRUE)
      }
      t0 <- Sys.time()
      res_list <- future.apply::future_lapply(
        chunks,
        FUN = function(chunk_g){
          out <- matrix(NA_real_, nrow=length(chunk_g), ncol=n,
                        dimnames=list(chunk_g, colnames(E)))
          for (ii in seq_along(chunk_g)) {
            g  <- chunk_g[ii]; v <- as.numeric(E[g, ]); sg <- sgs[g]
            Vnb <- v[J]; dim(Vnb) <- dim(J)
            Delta <- Vnb - v
            Sim <- exp(-(Delta*Delta)/(2*sg*sg))
            Wbl <- Sim * Wsp_block
            rsb <- rowSums(Wbl); rsb[!is.finite(rsb) | rsb==0] <- 1
            Wbl <- Wbl / rsb
            out[ii, ] <- rowSums(Wbl * Vnb)
          }
          out
        },
        future.seed = TRUE
      )
      res <- do.call(rbind, res_list)
      if (verbose) {
        el <- as.numeric(difftime(Sys.time(), t0, units="secs"))
        message(sprintf("    [parallel] %d chunks done in %.1fs (workers=%d)",
                        length(chunks), el, workers))
      }
      res
    }
    
    res_mat <- if (isTRUE(parallel)) {
      bilateral_par(markers, workers=workers, chunk_size=chunk_size)
    } else {
      bilateral_seq(markers)
    }
    
    # 回填到 E_sm
    for (g in rownames(res_mat)) E_sm[g, ] <- res_mat[g, ]
    set_data_layer(o, assay=assay, mat=E_sm)
  }
  if (!is.null(cells_use)) {
    keep <- intersect(colnames(obj), cells_use); obj <- obj[, keep]
    if (verbose) message(sprintf("[subset] cells_use: %d spots", ncol(obj)))
  } else if (!is.null(subset)) {
    md <- obj@meta.data
    keep <- with(md, eval(parse(text=subset)))
    if (!is.logical(keep) || length(keep)!=nrow(md))
      stop("subset must evaluate to logical(nrow(meta.data)).")
    obj <- obj[, which(keep)]
    if (verbose) message(sprintf("[subset] expr '%s': %d spots", subset, ncol(obj)))
  }

  obj <- ensure_spatial_logdata(obj, assay=assay)
  DefaultAssay(obj) <- assay

  hvgs <- get_hvgs_ffpe(obj, n=n_hvg, assay=assay)
  if (!length(hvgs)) stop("No HVGs selected. Increase n_hvg or check normalization.")
  VariableFeatures(obj) <- hvgs
  if (verbose) message(sprintf("[HVG] %d features", length(hvgs)))
  if (identical(bilateral_markers, ".__USE_TOP_HVGS__")) {
    bilateral_markers <- head(hvgs, max_bilateral_markers)
  }

  if (ncol(obj) <= max(k_smooth, k_expr, k_sp))
    stop("Too few cells/spots for chosen k. Reduce k_* or increase subset size.")

  obj_smoothed <- obj
  if (smooth_mode == "data") {
    if (verbose) message("[smooth_mode] data: Gaussian on HVGs ± bilateral on markers")
    obj <- knn_smooth_data(
      obj, k=k_smooth, hvgs=hvgs,
      bilateral_markers=bilateral_markers,
      assay=assay, verbose=verbose, progress_every=progress_every,
      max_bilat = max_bilateral_markers,
      parallel = parallel, workers = workers, chunk_size = chunk_size,
      manage_plan = manage_plan
    )
    obj_smoothed <- obj
    if (!is.null(save_smoothed_to)) {
      assign(save_smoothed_to, obj_smoothed, envir = .GlobalEnv)
      if (verbose) message(sprintf("[saved] smoothed object -> %s", save_smoothed_to))
    }
    obj <- ScaleData(obj, assay=assay, features=hvgs, verbose=FALSE)
    obj <- RunPCA(obj, assay=assay, features=hvgs, npcs=max(dims_use), verbose=FALSE)
  } else {
    if (verbose) message("[smooth_mode] pca: smooth PCA embeddings only; data unchanged")
    obj <- ScaleData(obj, assay=assay, features=hvgs, verbose=FALSE)
    obj <- RunPCA(obj, assay=assay, features=hvgs, npcs=max(dims_use), verbose=FALSE)
    xy <- get_xy(obj); n <- nrow(xy)
    nn <- FNN::get.knn(xy, k=k_sp)
    d1 <- nn$nn.dist[,1]
    sigma_sp <- stats::median(d1[d1>0], na.rm=TRUE); if (!is.finite(sigma_sp)||sigma_sp<=0) sigma_sp <- 1
    i <- rep(seq_len(n), each=k_sp); j <- as.vector(nn$nn.index)
    w <- exp(-(as.vector(nn$nn.dist)^2)/(2*sigma_sp^2))
    Wsp_pca <- Matrix::sparseMatrix(i=i, j=j, x=w, dims=c(n,n))
    rW <- 1/Matrix::rowSums(Wsp_pca); rW[!is.finite(rW)] <- 0
    Wsp_pca <- Matrix::Diagonal(x=rW) %*% Wsp_pca
    pcs <- Embeddings(obj,"pca")[, seq_len(max(dims_use)), drop=FALSE]
    pcs_sm <- as.matrix(Wsp_pca %*% pcs)
    obj@reductions$pca@cell.embeddings[, seq_len(ncol(pcs_sm))] <- pcs_sm
  }
   out <- tryCatch({
    # 表达 SNN
    obj <- FindNeighbors(obj, dims=dims_use, k.param=k_expr, graph.name=NULL, verbose=FALSE)
    snn_name <- grep("_snn$", names(obj@graphs), value=TRUE)[1]
    SNN <- as(obj@graphs[[snn_name]], "dgCMatrix")
    cells <- colnames(obj)
    SNN <- ensure_dimnames_square(Matrix::drop0(SNN), cells)
    # 空间图
    xy <- get_xy(obj)
    nn <- FNN::get.knn(xy, k=k_sp)
    i <- rep(seq_len(nrow(xy)), each=k_sp)
    j <- as.vector(nn$nn.index)
    d <- as.vector(nn$nn.dist)
    sigma_sp <- stats::median(nn$nn.dist[,1][nn$nn.dist[,1] > 0], na.rm=TRUE)
    if (!is.finite(sigma_sp) || sigma_sp <= 0) sigma_sp <- 1
    w <- exp(-(d^2)/(2*sigma_sp^2))
    Wsp <- Matrix::sparseMatrix(i=i, j=j, x=w, dims=c(nrow(xy), nrow(xy)))
    Wsp <- (Wsp + Matrix::t(Wsp))/2; diag(Wsp) <- 0
    Wsp <- ensure_dimnames_square(Matrix::drop0(Wsp), cells)
    # 归一 & 融合
    smax <- if (length(SNN@x)) max(SNN@x) else 1
    wmax <- if (length(Wsp@x)) max(Wsp@x) else 1
    if (!is.finite(smax) || smax<=0) smax <- 1
    if (!is.finite(wmax) || wmax<=0) wmax <- 1
    SNNs  <- SNN / smax
    Wsp_s <- Wsp / wmax
    Wmix <- alpha*SNNs + (1-alpha)*Wsp_s
    Wmix <- ensure_dimnames_square(Matrix::drop0(Wmix), cells)
    Wmix <- (Wmix + Matrix::t(Wmix))/2; diag(Wmix) <- 0; Wmix <- Matrix::drop0(Wmix)
    obj  <- put_custom_graph(obj, Wmix)
    # Leiden
    obj <- FindClusters(obj, graph.name="custom_graph", resolution=resolution,
                        algorithm=4, random.seed=seed, verbose=FALSE)
    # UMAP（默认强制 uwot，避免 Python 崩）
    if (do_umap) {
      if (isTRUE(force_uwot)) {
        obj <- RunUMAP(
          obj, reduction="pca", dims=dims_use,
          umap.method="uwot", n.neighbors=max(15, k_expr),
          min.dist=0.35, verbose=FALSE
        )
      } else {
        obj <- RunUMAP(
          obj,
          graph       = "custom_graph",
          umap.method = "umap-learn",
          n.neighbors = max(15, k_expr),
          min.dist    = 0.35,
          verbose     = FALSE
        )
      }
    }
    # UCell
    if (!is.null(sig_list) && length(sig_list)) {
      if (requireNamespace("UCell", quietly=TRUE)) {
        obj <- UCell::AddModuleScore_UCell(obj, features=sig_list, name=NULL, assay=assay)
        for (nm in names(sig_list)) {
          if (nm %in% colnames(obj@meta.data))
            colnames(obj@meta.data)[colnames(obj@meta.data)==nm] <- paste0("Score_", nm)
        }
      } else if (verbose) message("[UCell] not installed; skip.")
    }
    if (verbose) {
      p <- tryCatch(DimPlot(obj, group.by="seurat_clusters", label=TRUE) +
                      ggtitle(sprintf("Leiden (custom graph) | smooth_mode=%s", smooth_mode)),
                    error=function(e) NULL)
      if (!is.null(p)) print(p)
    }
    list(
      object = obj,
      hvgs   = hvgs,
      params = list(n_hvg=n_hvg, smooth_mode=smooth_mode, k_smooth=k_smooth,
                    bilateral_markers=bilateral_markers, max_bilateral_markers=max_bilateral_markers,
                    dims_use=dims_use, k_expr=k_expr, k_sp=k_sp, alpha=alpha,
                    resolution=resolution, do_umap=do_umap, parallel=parallel,
                    workers=workers, chunk_size=chunk_size),
      graphs = NULL,
      smoothed_object = obj_smoothed,
      error = NULL
    )
  }, error = function(e){
    warning(sprintf("[run_cellbin_ffpe] downstream failed: %s\nReturning SMOOTHED object.", conditionMessage(e)))
    list(
      object = obj_smoothed,
      hvgs   = hvgs,
      params = list(n_hvg=n_hvg, smooth_mode=smooth_mode, k_smooth=k_smooth,
                    bilateral_markers=bilateral_markers, max_bilateral_markers=max_bilateral_markers,
                    dims_use=dims_use, k_expr=k_expr, k_sp=k_sp, alpha=alpha,
                    resolution=resolution, do_umap=do_umap, parallel=parallel,
                    workers=workers, chunk_size=chunk_size),
      graphs = NULL,
      smoothed_object = obj_smoothed,
      error = conditionMessage(e)
    )
  })
  
  invisible(out)
}

# ==== ?????????? ====















suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(igraph)
})

score_lineages_get_margin <- function(object, panels, assay = NULL) {
  if (!is.null(assay)) DefaultAssay(object) <- assay
  if (requireNamespace("UCell", quietly = TRUE)) {
    object <- UCell::AddModuleScore_UCell(object, features = panels, name = NULL, ncores = 1)
    S <- as.data.frame(object@meta.data[, names(panels), drop = FALSE])
  } else {
    E <- GetAssayData(object, slot = "data")
    zmean <- function(g) {
      g2 <- intersect(g, rownames(E)); if (length(g2) < 3) return(rep(NA_real_, ncol(E)))
      M <- as.matrix(E[g2, , drop = FALSE]); mu <- rowMeans(M); sd <- apply(M,1,sd); sd[sd<1e-6] <- 1e-6
      colMeans( sweep(M,1,mu,"-")/sd )
    }
    S <- sapply(panels, zmean) %>% as.data.frame()
    colnames(S) <- names(panels); rownames(S) <- colnames(object)
  }
  ord <- t(apply(as.matrix(S), 1, function(v){
    o <- order(v, decreasing = TRUE); c(o[1], o[2])
  }))
  top1  <- mapply(function(i,j) names(S)[i], ord[,1], ord[,2])
  top2  <- mapply(function(i,j) names(S)[j], ord[,1], ord[,2])
  v1    <- S[cbind(seq_len(nrow(S)), ord[,1])]
  v2    <- S[cbind(seq_len(nrow(S)), ord[,2])]
  margin <- v1 - v2
  tibble::tibble(
    barcode = rownames(S),
    lineage_top = top1,
    lineage_second = top2,
    lineage_margin = as.numeric(margin)
  ) %>% dplyr::left_join(S %>% tibble::rownames_to_column("barcode"), by = "barcode")
}






local_consistency <- function(object, graph.name = "custom", idents = NULL){
  A <- object@graphs[[graph.name]]
  if (is.null(A)) stop("Graph not found: ", graph.name)
  labs <- idents %||% Idents(object)
  labs <- as.factor(labs)
  idx_list <- lapply(seq_len(nrow(A)), function(i){ which(A[i, ] != 0) })
  cons <- vapply(seq_len(nrow(A)), function(i){
    nb <- idx_list[[i]]; if (length(nb) == 0) return(NA_real_); mean(labs[nb] == labs[i])
  }, numeric(1))
  tibble::tibble(barcode = colnames(A), local_consistency = cons)
}




summarize_and_flag_clusters <- function(object, lc_df, lin_df, thresholds, cluster_col = NULL){
  labs <- if (is.null(cluster_col)) Idents(object) else factor(object@meta.data[[cluster_col]])
  df <- lc_df %>% inner_join(lin_df[, c("barcode","lineage_margin")], by="barcode") %>%
    mutate(cluster = as.character(labs[barcode]))
  sum_tab <- df %>% group_by(cluster) %>%
    summarise(n = n(),
              med_margin = median(lineage_margin, na.rm = TRUE),
              med_cons   = median(local_consistency, na.rm = TRUE),
              .groups = "drop")
  # 规则
  sum_tab <- sum_tab %>% mutate(
    flag = case_when(
      med_cons >= thresholds$consistency_interface_min &
        med_margin <= thresholds$margin_interface_max ~ "interface_candidate",
      med_cons <= thresholds$consistency_mixed_max &
        med_margin <= thresholds$margin_mixed_max ~ "technical_mixed_candidate",
      TRUE ~ "clean"
    )
  )
  sum_tab
}




neighbor_vote_reassign <- function(object, target_clusters, vote_majority = 0.6, graph.name = "custom"){
  A <- object@graphs[[graph.name]]; if (is.null(A)) stop("Graph not found: ", graph.name)
  cur <- as.character(Idents(object))
  changed <- rep(FALSE, length(cur))
  names(changed) <- names(cur)
  for (cl in target_clusters) {
    cells <- names(cur)[cur == cl]
    if (length(cells) == 0) next
    idx   <- match(cells, colnames(A))
    for (ii in seq_along(idx)) {
      i <- idx[ii]
      nb <- which(A[i, ] != 0)
      if (length(nb) == 0) next
      nb_lab <- cur[colnames(A)[nb]]
      # 排除自身簇的票
      nb_lab2 <- nb_lab[nb_lab != cl]
      if (length(nb_lab2) == 0) next
      tab <- sort(table(nb_lab2), decreasing = TRUE)
      winner <- names(tab)[1]; frac <- as.numeric(tab[1]) / length(nb)
      if (frac >= vote_majority) {
        cur[cells[ii]] <- winner
        changed[cells[ii]] <- TRUE
      }
    }
  }
  list(new_labels = cur, changed = changed)
}




component_merge <- function(object, target_clusters, comp_majority = 0.6, graph.name = "custom"){
  A <- object@graphs[[graph.name]]; if (is.null(A)) stop("Graph not found: ", graph.name)
  cur <- as.character(Idents(object))
  for (cl in target_clusters) {
    cells <- names(cur)[cur == cl]
    if (length(cells) < 2) next
    sub <- A[cells, cells]; sub[sub != 0] <- 1
    G <- igraph::graph_from_adjacency_matrix(as.matrix(sub), mode="undirected", diag=FALSE)
    comp <- igraph::components(G)$membership
    comps <- split(cells, comp)
    for (comp_cells in comps) {
      # 该组件与外界的边
      idx <- match(comp_cells, colnames(A))
      nb_out <- unique(unlist(lapply(idx, function(i) setdiff(which(A[i, ] != 0), match(comp_cells, colnames(A))))))
      if (length(nb_out) == 0) next
      nb_lab <- cur[colnames(A)[nb_out]]
      tab <- sort(table(nb_lab), decreasing=TRUE)
      winner <- names(tab)[1]; frac <- as.numeric(tab[1]) / length(nb_out)
      if (winner != cl && frac >= comp_majority) {
        cur[comp_cells] <- winner
      }
    }
  }
  cur
}





# cluster_col: 若你不想用 Idents(object) 的簇列，可指定 meta.data 的列名
second_pass_refine <- function(
    object,
    panels,
    assay_for_scores = NULL,     # NULL=DefaultAssay（建议 "SpatialSmoothed"）
    graph.name = "custom",
    cluster_col = NULL,
    thresholds = list(
      margin_interface_max   = 0.20,
      consistency_interface_min = 0.60,
      margin_mixed_max       = 0.15,
      consistency_mixed_max  = 0.60,
      vote_majority          = 0.60,
      comp_majority          = 0.60,
      unknown_margin_max     = 0.15
    ),
    interface_label = "Interface",     # 界面簇追加的前缀
    unknown_label   = "Unknown_or_Ambient"
){
  # 1) 评分 & 一致率
  lin_df <- score_lineages_get_margin(object, panels = panels, assay = assay_for_scores)
  lc_df  <- local_consistency(object, graph.name = graph.name,
                              idents = if (is.null(cluster_col)) Idents(object) else object@meta.data[[cluster_col]])
  
  # 2) 簇级汇总与打标
  sum_tab <- summarize_and_flag_clusters(object, lc_df, lin_df, thresholds, cluster_col)
  interface_candidates <- sum_tab$cluster[sum_tab$flag == "interface_candidate"]
  mixed_candidates     <- sum_tab$cluster[sum_tab$flag == "technical_mixed_candidate"]
  
  # 3) 备份一列初始簇
  if (!"cluster_round1" %in% colnames(object@meta.data)) {
    object$cluster_round1 <- as.character(if (is.null(cluster_col)) Idents(object) else object@meta.data[[cluster_col]])
  }
  new_labels <- object$cluster_round1
  
  # 4) 对“技术混杂候选”先做邻居投票
  if (length(mixed_candidates) > 0) {
    nv <- neighbor_vote_reassign(object, mixed_candidates, vote_majority = thresholds$vote_majority, graph.name = graph.name)
    new_labels <- nv$new_labels
    Idents(object) <- factor(new_labels)
    # 5) 组件级合并（在投票之后）
    new_labels <- component_merge(object, mixed_candidates, comp_majority = thresholds$comp_majority, graph.name = graph.name)
  }
  
  # 6) Unknown 放弃阈值（margin 很低且仍未被回收）
  low_margin_cells <- lin_df$barcode[lin_df$lineage_margin <= thresholds$unknown_margin_max]
  still_mixed <- names(new_labels)[new_labels %in% mixed_candidates & names(new_labels) %in% low_margin_cells]
  if (length(still_mixed) > 0) new_labels[still_mixed] <- unknown_label
  
  # 7) 界面簇直接保留并加前缀（可选）
  if (length(interface_candidates) > 0) {
    mask <- new_labels %in% interface_candidates
    new_labels[mask] <- paste0(interface_label, "_", new_labels[mask])
  }
  
  # 8) 写回对象
  object$cluster_round2 <- new_labels
  Idents(object) <- factor(object$cluster_round2)
  
  # 9) 汇总报告
  final_summary <- data.frame(
    cluster = levels(factor(object$cluster_round2)),
    n = as.numeric(table(object$cluster_round2)[levels(factor(object$cluster_round2))]),
    stringsAsFactors = FALSE
  ) %>% left_join(sum_tab, by = c("cluster" = "cluster"))
  
  list(
    object = object,
    cluster_summary_round1 = sum_tab,
    cluster_summary_round2 = final_summary
  )
}







suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })


.find_rctd_cols <- function(md,
                            first_w_cands  = c("first_w","RCTD_first_w","first_weight","first.score","best_score"),
                            second_w_cands = c("second_w","RCTD_second_w","second_weight","second.score"),
                            margin_cands   = c("margin","RCTD_margin","first_margin")){
  list(
    first_w  = first_w_cands [first_w_cands  %in% colnames(md)][1],
    second_w = second_w_cands[second_w_cands %in% colnames(md)][1],
    margin   = margin_cands  [margin_cands   %in% colnames(md)][1]
  )
}





.to_macro_from_patterns <- function(x, patterns, include_other = TRUE){
  # NA-safe mapping by regex pattern list: patterns is a named list, e.g. list(Tumor=c("tumor","epit"), Acinar=c("acinar", ...))
  out <- rep(NA_character_, length(x))
  ok  <- !is.na(x)
  s   <- tolower(as.character(x[ok]))
  hit_any <- rep(FALSE, sum(ok))
  for (nm in names(patterns)){
    rgx <- paste0("(", paste(patterns[[nm]], collapse = "|"), ")")
    m   <- grepl(rgx, s, ignore.case = TRUE, perl = TRUE)
    out[ok][m & !hit_any] <- nm
    hit_any <- hit_any | m
  }
  if (include_other){
    out[ok][!hit_any] <- "Other"
  }
  out
}





# ---------- main unified function ----------
Generate_Anchor_Annotations <- function(
    object,
    label_col          = "RCTD_first",
    spotclass_col      = "RCTD_spot_class",
    weight_col_candidates = list(
      first_w  = c("first_w","RCTD_first_w","first_weight","first.score","best_score"),
      second_w = c("second_w","RCTD_second_w","second_weight","second.score"),
      margin   = c("margin","RCTD_margin","first_margin")
    ),
    macro_patterns = list(
      Tumor   = c("tumou?r","epit","duct","classical","basal","squamoid","emt","proliferative","ifn"),
      Acinar  = c("acinar","adm","exocrine","reg","spink1","reg1[ab]?","reg3[ab]?"),
      CAF     = c("fibro","\\bcaf\\b","strom","mycaf","icaf","perivascular","pericyte","thbs1","ccn1"),
      Immune  = c("immune","tam","macroph","\\bt cell\\b","\\bb cell\\b","plasma","dendritic","nk","mono","neutro"),
      Endothelial = c("endothel","pecam1","cldn5","vwf","kdr","ca4"),
      Endocrine   = c("endocr","islet","alpha","beta","delta","pp","chga","ins","gcg","sst"),
      SmoothMyo   = c("smooth","myo","acta2","tagln","myh11","cspg4"),
      Schwann     = c("schwann","mpz","plp1","mbp")
    ),
    macro_thresholds = list(
      w = c(Tumor=0.55, Acinar=0.55, CAF=0.60, Immune=0.60, Endothelial=0.60, Endocrine=0.60, SmoothMyo=0.60, Schwann=0.60),
      m = c(Tumor=0.15, Acinar=0.15, CAF=0.20, Immune=0.20, Endothelial=0.20, Endocrine=0.20, SmoothMyo=0.20, Schwann=0.20)
    ),
    recover_thresholds = list(first_w = 0.70, margin = 0.20),
    keepable_spot_classes = c("singlet","doublet_uncertain"),
    reject_spot_classes   = c("reject","doublet_certain"),
    unknown_margin_max    = 0.15,
    out_prefix            = "Type0",
    palette = c("Tumor"="#d73027","Acinar"="#1f78b4","CAF"="#33a02c","Immune"="#6a3d9a",
                "Endothelial"="#1b9e77","Endocrine"="#e6ab02","SmoothMyo"="#a6761d","Schwann"="#66a61e",
                "Other"="#a6cee3","Uncertain"="#bdbdbd"),
    plot    = TRUE,
    verbose = TRUE
){
  stopifnot(label_col %in% colnames(object@meta.data),
            spotclass_col %in% colnames(object@meta.data))
  
  md   <- object@meta.data
  cols <- .find_rctd_cols(md,
                          first_w_cands  = weight_col_candidates$first_w,
                          second_w_cands = weight_col_candidates$second_w,
                          margin_cands   = weight_col_candidates$margin)
  if (is.na(cols$first_w)) stop("Cannot find 'first_w' column (try adjusting weight_col_candidates).")
  
  first_w <- md[[cols$first_w]]
  # compute margin if missing (needs second_w)
  margin  <- if (!is.na(cols$margin)) {
    md[[cols$margin]]
  } else {
    if (is.na(cols$second_w)) stop("Neither 'margin' nor 'second_w' found to compute margin.")
    md[[cols$first_w]] - md[[cols$second_w]]
  }
  
  rctd_label <- md[[label_col]]
  has_label  <- !is.na(rctd_label) & nzchar(as.character(rctd_label))
  
  # NA-safe macro mapping
  macro <- .to_macro_from_patterns(rctd_label, macro_patterns, include_other = TRUE)
  
  # spot classes
  spotc <- as.character(md[[spotclass_col]])
  keepable <- ifelse(is.na(spotc), FALSE, spotc %in% keepable_spot_classes)
  rejected <- ifelse(is.na(spotc), FALSE, spotc %in% reject_spot_classes)
  is_singlet <- !is.na(spotc) & spotc == "singlet"
  is_doublet_uncertain <- !is.na(spotc) & spotc == "doublet_uncertain"
  
  # pass thresholds per macro (w & m)
  n    <- nrow(md)
  pass <- logical(n)
  in_set <- !is.na(macro) & macro %in% union(names(macro_thresholds$w), names(macro_thresholds$m))
  pass[in_set] <- ( first_w[in_set] >= macro_thresholds$w[ macro[in_set] ] ) &
    ( margin [in_set] >= macro_thresholds$m[ macro[in_set] ] )
  pass[is.na(pass)] <- FALSE
  
  # recovery for doublet_uncertain
  can_recover <- is_doublet_uncertain &
    !is.na(first_w) & (first_w >= recover_thresholds$first_w) &
    !is.na(margin ) & (margin  >= recover_thresholds$margin)
  can_recover[is.na(can_recover)] <- FALSE
  
  # final accept
  accept <- has_label & keepable & ( (is_singlet & pass) | can_recover )
  accept[is.na(accept)] <- FALSE
  
  # reasoning
  reason <- rep("below_threshold", n)
  reason[!has_label]                 <- "no_rctd_label"
  reason[rejected]                   <- spotc[rejected]
  reason[is_doublet_uncertain & !can_recover] <- "uncertain_not_recovered"
  reason[accept]                     <- "accepted"
  
  # output columns
  out_label  <- rep("Uncertain", n)
  out_label[accept] <- as.character(rctd_label[accept])
  
  md[[paste0(out_prefix, "_macro")]]   <- macro
  md[[paste0(out_prefix, "_first_w")]] <- first_w
  md[[paste0(out_prefix, "_margin")]]  <- margin
  md[[paste0(out_prefix, "_reason")]]  <- reason
  md[[out_prefix]]                     <- out_label
  object@meta.data <- md
  
  # quick summaries
  sum_tab   <- as.data.frame(table(object@meta.data[[out_prefix]], useNA = "ifany"))
  reason_tb <- as.data.frame(table(object@meta.data[[paste0(out_prefix, "_reason")]], useNA = "ifany"))
  
  if (isTRUE(verbose)) {
    message(sprintf("[Generate_Anchor_Annotations] accepted = %d / %d (%.1f%%)",
                    sum(accept), n, 100*mean(accept)))
  }
  
  if (isTRUE(plot)) {
    # spatial plot if user function exists
    if ("PlotSpatialDistribution" %in% ls(envir = .GlobalEnv)) {
      try(print(PlotSpatialDistribution(object, group_by = out_prefix, palette = palette, ptsize = 1.8, show_legend = TRUE)), silent = TRUE)
    }
    print(ggplot(object@meta.data, aes_string(x = paste0(out_prefix, "_first_w"))) +
            geom_histogram(bins = 60) +
            labs(x = "first_w", y = "Count") + theme_classic())
    print(ggplot(object@meta.data, aes_string(x = paste0(out_prefix, "_margin"))) +
            geom_histogram(bins = 60) +
            labs(x = "margin (first-second)", y = "Count") + theme_classic())
  }
  
  return(list(
    object           = object,
    summary          = sum_tab,
    reason           = reason_tb,
    thresholds_used  = list(macro_thresholds = macro_thresholds,
                            recover_thresholds = recover_thresholds,
                            unknown_margin_max = unknown_margin_max,
                            keepable_spot_classes = keepable_spot_classes,
                            reject_spot_classes = reject_spot_classes),
    columns_used     = list(label_col = label_col, spotclass_col = spotclass_col,
                            first_w = cols$first_w, second_w = cols$second_w, margin = cols$margin),
    mapping_preview  = head(data.frame(label = rctd_label, macro = macro))
  ))
}

# ==== BANKSY/?????? ====








# =============================
# audit_pool(): audit an expanded pool
#   - Summarizes how much of a pool comes from anchors vs neighborhood
#   - Checks contamination (RBC/Platelet, epithelial, acinar)
#   - Checks pool-specific positive markers (immune / CAF / tumor / acinar)
#   - Optionally inspects neighborhood alignment and local purity if available
# Requirements:
#   object@meta.data has: Type0_locked (anchors), optional: nbr_major_label/prop, local_purity_major_k
#   DefaultAssay(object) points to the count assay used for marker presence
# All comments in English (as requested).
# =============================

# robust presence helper (auto-ignores missing genes)
if (!exists("safe_presence_vec")) {
  safe_presence_vec <- function(obj, cells, genes, assay="Spatial", min_umi=1L){
    M <- GetAssayData(obj, assay=assay, slot="counts")
    if (!length(cells)) return(logical(0))
    genes_in <- intersect(genes, rownames(M))
    if (length(genes_in) == 0) return(rep(FALSE, length(cells)))
    as.vector(Matrix::colSums(M[genes_in, cells, drop=FALSE] >= min_umi) > 0)
  }
}

audit_pool <- function(object,
                       pool_cells,
                       anchor_label,
                       assay      = "Spatial",
                       min_umi    = 1L,     # marker presence threshold (UMI)
                       purity_col = "local_purity_major_k",
                       nbr_labcol = "nbr_major_label",
                       nbr_pcol   = "nbr_major_prop") {
  md  <- object@meta.data
  all_cells <- colnames(object)
  
  # sanity
  if (!"Type0_locked" %in% colnames(md)) stop("meta.data missing 'Type0_locked'.")
  pool_cells <- intersect(pool_cells, all_cells)
  if (!length(pool_cells)) {
    return(list(
      summary = data.frame(pool=anchor_label, total=0, anchors=0, expanded=0,
                           expanded_pct=NA, pos_mark_pct=NA, rbcplt_pct=NA,
                           epi_pct=NA, acinar_pct=NA, nbr_match_pct=NA,
                           nbr_prop_mean=NA, purity_mean=NA),
      used_markers = list(),
      expanded_cells = character(0)
    ))
  }
  
  # anchors vs expanded
  anchors  <- WhichCells(object, expression = Type0_locked == anchor_label)
  anchors  <- intersect(anchors, pool_cells)
  expanded <- setdiff(pool_cells, anchors)
  
  pct <- function(x) if (!length(x)) NA_real_ else mean(x, na.rm=TRUE)
  
  # marker panels (with fallbacks)
  panel <- list(
    immune_pos = list( # union of these defines "immune"
      ptprc   = c("PTPRC"),
      Tcell   = c("CD3D","CD3E","TRAC"),
      Bcell   = c("MS4A1","CD79A","CD79B"),
      myeloid = c("LST1","LYZ")
    ),
    caf_pos    = c("COL1A1","DCN","LUM"),
    tumor_pos  = c("EPCAM","KRT19","KRT8","KRT18","KRT7","MUC1"),
    acinar_pos = c("CPA1","PRSS1","CTRB1","REG1A"),
    rbc        = c("HBB","HBA1","HBA2","HBD","ALAS2","SLC4A1","GYPA"),
    platelet   = c("PPBP","PF4","TUBB1","ITGA2B","ITGB3","GP9","SPARC"),
    epithelial = c("EPCAM","KRT19","KRT8","KRT18","KRT7"),
    acinar     = c("CPA1","PRSS1","CTRB1","REG1A")
  )
  
  # compute presence on expanded cells
  exp_cells <- expanded
  # If there is no expansion, we will still compute markers on an empty set.
  pos_immune <- if (length(exp_cells)) {
    p_ptprc   <- safe_presence_vec(object, exp_cells, panel$immune_pos$ptprc,   assay, min_umi)
    p_T       <- safe_presence_vec(object, exp_cells, panel$immune_pos$Tcell,   assay, min_umi)
    p_B       <- safe_presence_vec(object, exp_cells, panel$immune_pos$Bcell,   assay, min_umi)
    p_my      <- safe_presence_vec(object, exp_cells, panel$immune_pos$myeloid, assay, min_umi)
    (p_ptprc | p_T | p_B | p_my)
  } else logical(0)
  
  pos_caf     <- safe_presence_vec(object, exp_cells, panel$caf_pos,    assay, min_umi)
  pos_tumor   <- safe_presence_vec(object, exp_cells, panel$tumor_pos,  assay, min_umi)
  pos_acinar  <- safe_presence_vec(object, exp_cells, panel$acinar_pos, assay, min_umi)
  pos_rbc     <- safe_presence_vec(object, exp_cells, panel$rbc,        assay, min_umi)
  pos_plate   <- safe_presence_vec(object, exp_cells, panel$platelet,   assay, min_umi)
  pos_rbcplt  <- (pos_rbc | pos_plate)
  pos_epi     <- safe_presence_vec(object, exp_cells, panel$epithelial, assay, min_umi)
  pos_acin    <- safe_presence_vec(object, exp_cells, panel$acinar,     assay, min_umi)
  
  # choose pool-specific positive marker rate
  pos_rate <- switch(anchor_label,
                     "Immune" = pct(pos_immune),
                     "CAF"    = pct(pos_caf),
                     "Tumor"  = pct(pos_tumor),
                     "Acinar" = pct(pos_acinar),
                     pct(logical(0))
  )
  
  # neighborhood alignment (optional)
  nbr_match <- nbr_prop <- purity <- NA_real_
  if (all(c(nbr_labcol, nbr_pcol) %in% colnames(md)) && length(exp_cells)) {
    nbr_labs <- md[exp_cells, nbr_labcol, drop=TRUE]
    nbr_props<- md[exp_cells, nbr_pcol,   drop=TRUE]
    nbr_match <- pct(nbr_labs == anchor_label)
    nbr_prop  <- mean(nbr_props, na.rm=TRUE)
  }
  if (purity_col %in% colnames(md) && length(exp_cells)) {
    purity <- mean(md[exp_cells, purity_col, drop=TRUE], na.rm=TRUE)
  }
  
  # assemble summary
  sm <- data.frame(
    pool          = anchor_label,
    total         = length(pool_cells),
    anchors       = length(anchors),
    expanded      = length(expanded),
    expanded_pct  = round(100*length(expanded)/max(1, length(pool_cells)), 1),
    pos_mark_pct  = round(100*pos_rate, 1),              # pool-specific positive rate in expanded
    rbcplt_pct    = round(100*pct(pos_rbcplt), 1),       # RBC/Platelet in expanded
    epi_pct       = round(100*pct(pos_epi), 1),          # epithelial in expanded
    acinar_pct    = round(100*pct(pos_acin), 1),         # acinar in expanded
    nbr_match_pct = round(100*nbr_match, 1),
    nbr_prop_mean = ifelse(is.na(nbr_prop), NA, round(nbr_prop, 3)),
    purity_mean   = ifelse(is.na(purity), NA, round(purity, 3))
  )
  
  # record which markers actually existed in matrix
  M <- GetAssayData(object, assay=assay, slot="counts")
  used_markers <- list(
    immune_ptprc = intersect(panel$immune_pos$ptprc,   rownames(M)),
    immune_T     = intersect(panel$immune_pos$Tcell,   rownames(M)),
    immune_B     = intersect(panel$immune_pos$Bcell,   rownames(M)),
    immune_myeloid = intersect(panel$immune_pos$myeloid, rownames(M)),
    caf_pos      = intersect(panel$caf_pos,    rownames(M)),
    tumor_pos    = intersect(panel$tumor_pos,  rownames(M)),
    acinar_pos   = intersect(panel$acinar_pos, rownames(M)),
    rbc          = intersect(panel$rbc,        rownames(M)),
    platelet     = intersect(panel$platelet,   rownames(M))
  )
  
  list(summary = sm, used_markers = used_markers, expanded_cells = expanded)
}





# =============================
# (optional) convenience wrapper to audit multiple pools at once
# pools is a named list: list(Immune=Immune_pool, CAF=CAF_pool, Tumor=Tumor_pool, Acinar=Acinar_pool)
audit_pools <- function(object, pools, assay="Spatial", min_umi=1L,
                        purity_col="local_purity_major_k",
                        nbr_labcol="nbr_major_label", nbr_pcol="nbr_major_prop"){
  out <- lapply(names(pools), function(lbl){
    res <- audit_pool(object, pools[[lbl]], lbl, assay=assay, min_umi=min_umi,
                      purity_col=purity_col, nbr_labcol=nbr_labcol, nbr_pcol=nbr_pcol)
    res$summary
  })
  do.call(rbind, out)
}







suppressPackageStartupMessages({
  library(Seurat); library(Matrix)
  library(SpatialExperiment); library(SingleCellExperiment)
  library(Banksy); library(FNN); library(dplyr)
})

# ===== BANKSY: prepare once + sweep grid + pick best + write back =====
banksy_pool_sweep <- function(
    obj, cells_use, pool_name,
    assay        = "Spatial",
    k_geom       = c(6L,12L),
    M            = 1L,
    use_agf      = TRUE,
    lambdas      = c(0.20, 0.25, 0.30),
    resolutions  = c(0.2, 0.4, 0.6),
    k_neighbors  = 30L,
    npcs         = 30L,
    k_for_score  = 10L,     # for spatial continuity score
    seed         = 1L,
    add_umap     = TRUE,
    attach_best_col = TRUE, # write best column back to obj meta
    write_all_cols  = FALSE # set TRUE if you want all grid labels written back
){
  stopifnot(length(cells_use) > 0)
  set.seed(seed)
  
  # ---------- helpers ----------
  .pick_banksy_umap_name <- function(spe, M, lambda){
    rdn  <- reducedDimNames(spe)
    cand <- grep(paste0("^UMAP_M", M, "_lam"), rdn, value=TRUE)
    if (!length(cand)) return(NULL)
    get_lam <- function(nm) suppressWarnings(as.numeric(sub(".*_lam([0-9.]+).*","\\1", nm)))
    lam_vals <- vapply(cand, get_lam, numeric(1)); lam_vals[is.na(lam_vals)] <- Inf
    cand[ which.min(abs(lam_vals - lambda)) ]
  }
  spatial_continuity_score <- function(coords, labels, k=10){
    idx <- FNN::get.knn(coords, k=k)$nn.index
    same <- sapply(seq_len(nrow(idx)), function(i){
      mean(labels[idx[i,]] == labels[i], na.rm=TRUE)
    })
    mean(same, na.rm=TRUE)
  }
  
  # ---------- subset & build SpatialExperiment ----------
  obj_sub <- subset(obj, cells = cells_use)
  Mx <- GetAssayData(obj_sub, assay=assay, slot="counts")
  if (!inherits(Mx, "dgCMatrix")) Mx <- as(Mx, "dgCMatrix")
  cd <- obj_sub@meta.data
  
  coord_cols <- if (all(c("imagecol_scaled","imagerow_scaled") %in% colnames(cd))) {
    c("imagecol_scaled","imagerow_scaled")
  } else c("imagecol","imagerow")
  coords <- as.matrix(cd[, coord_cols]); colnames(coords) <- c("x","y")
  
  spe <- SpatialExperiment(assays=list(counts=Mx), colData=cd)
  spatialCoords(spe) <- coords
  
  # ---------- compute BANKSY features once ----------
  spe <- computeBanksy(
    spe, assay_name="counts", coord_names=c("x","y"),
    compute_agf=use_agf, k_geom=k_geom, M=M
  )
  
  # ---------- sweep the grid ----------
  grid <- expand.grid(lambda = lambdas, res = resolutions)
  summaries <- vector("list", nrow(grid))
  label_cols <- character(nrow(grid))  # optional: store names to write back later
  
  for (i in seq_len(nrow(grid))) {
    lam <- grid$lambda[i]; rs <- grid$res[i]
    
    # PCA / cluster on prepared BANKSY features
    spe_i <- runBanksyPCA(spe, M=M, lambda=lam, npcs=npcs)
    spe_i <- clusterBanksy(spe_i, use_agf=use_agf, lambda=lam,
                           algo="leiden", k_neighbors=k_neighbors,
                           resolution=rs, seed=seed)
    
    # Optional UMAP
    umap <- NULL
    if (isTRUE(add_umap)) {
      spe_i <- runBanksyUMAP(spe_i, M=M, lambda=lam, npcs=npcs,
                             n_neighbors=k_neighbors, min_dist=0.3, seed=seed)
      un <- .pick_banksy_umap_name(spe_i, M, lam)
      if (!is.null(un)) umap <- as.matrix(reducedDim(spe_i, un))
    }
    
    # fetch labels
    cl_name <- clusterNames(spe_i)[1]
    labs <- as.character(colData(spe_i)[[cl_name]])
    
    # score
    scs <- spatial_continuity_score(spatialCoords(spe_i), labs, k=k_for_score)
    tab <- table(labs)
    summ <- tibble(
      pool               = pool_name,
      lambda             = lam,
      resolution         = rs,
      n_clusters         = length(tab),
      spatial_continuity = scs,
      min_cluster_size   = min(tab),
      median_cluster_size= median(tab)
    )
    summaries[[i]] <- summ
    
    # write back this column if requested
    lab_col <- sprintf("banksy_%s_lam%.2f_res%.2f", pool_name, lam, rs)
    label_cols[i] <- lab_col
    if (isTRUE(write_all_cols)) {
      obj[[lab_col]] <- NA_character_
      obj@meta.data[cells_use, lab_col] <- labs
    }
    
    # keep UMAP only for best later; not attached now to save memory
  }
  
  summary_df <- bind_rows(summaries)
  
  # ---------- pick best (maximize spatial_continuity; tie-break by larger median cluster) ----------
  ord <- order(-summary_df$spatial_continuity, -summary_df$median_cluster_size, summary_df$n_clusters)
  best <- summary_df[ord[1], , drop=FALSE]
  best_lab_col <- sprintf("banksy_%s_lam%.2f_res%.2f",
                          best$pool[1], best$lambda[1], best$resolution[1])
  
  # Recompute best run to get labels + UMAP to attach to subset
  spe_best <- runBanksyPCA(spe, M=M, lambda=best$lambda[1], npcs=npcs)
  spe_best <- clusterBanksy(spe_best, use_agf=use_agf, lambda=best$lambda[1],
                            algo="leiden", k_neighbors=k_neighbors,
                            resolution=best$resolution[1], seed=seed)
  umap_best <- NULL
  if (isTRUE(add_umap)) {
    spe_best <- runBanksyUMAP(spe_best, M=M, lambda=best$lambda[1], npcs=npcs,
                              n_neighbors=k_neighbors, min_dist=0.3, seed=seed)
    un <- .pick_banksy_umap_name(spe_best, M, best$lambda[1])
    if (!is.null(un)) umap_best <- as.matrix(reducedDim(spe_best, un))
  }
  cl_name <- clusterNames(spe_best)[1]
  best_labels <- as.character(colData(spe_best)[[cl_name]])
  
  # attach best column to full obj (if not already)
  if (isTRUE(attach_best_col) && !isTRUE(write_all_cols)) {
    obj[[best_lab_col]] <- NA_character_
    obj@meta.data[cells_use, best_lab_col] <- best_labels
  }
  # ...已写回到 obj 之后，紧接着加这一行：
  obj_sub[[best_lab_col]] <- best_labels
  # attach UMAP to subset for plotting
  if (!is.null(umap_best)) {
    obj_sub[["banksy_umap"]] <- CreateDimReducObject(
      embeddings = umap_best, key = "BK_", assay = assay
    )
  }
  
  list(
    object        = obj,                 # full Seurat with labels column(s)
    object_subset = obj_sub,             # subset with UMAP for visualization
    spe_best      = spe_best,            # SpatialExperiment for best setting
    summary       = summary_df,          # grid summary table
    best          = list(
      lambda   = best$lambda[1],
      resolution = best$resolution[1],
      lab_col  = best_lab_col
    )
  )
}





# ---------- helpers ----------
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





# ========================= ScoreModules (UCell 优先 + Seurat 命名修复) =========================
ScoreModules <- function(
    object,
    features_list,                           # named list: list(SigA=c(...), SigB=c(...))
    method        = c("auto","Seurat","UCell","zmean","ssGSEA","gsva","AUCell"),
    assay         = NULL,                    # auto-pick: SCT > Spatial > DefaultAssay
    slot_use      = "data",
    pre_normalize     = c("auto","always","skip"),
    set_default_assay = TRUE,
    name_prefix   = "Score",
    add_to_meta   = TRUE,
    return_matrix = TRUE,
    post_scale    = c("none","z","minmax"),
    ncores        = 1,
    verbose       = TRUE,
    # <<< changed: UCell highest priority
    auto_order    = c("UCell","Seurat","zmean"),
    min_overlap_prop = 0.40,
    min_genes_abs    = 3
){
  stopifnot(inherits(object, "Seurat"))
  
  # ---- assay auto-pick ----
  if (is.null(assay)) {
    if ("SCT" %in% names(object@assays)) assay <- "SCT"
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
  
  # ---- set assay & normalization policy ----
  if (set_default_assay && !is.null(assay)) Seurat::DefaultAssay(object) <- assay
  if (!assay %in% names(object@assays)) stop("[ScoreModules] Assay not found: ", assay)
  
  if (pre_normalize == "always") {
    object <- Seurat::NormalizeData(object, assay = assay, verbose = FALSE)
  } else if (pre_normalize == "auto") {
    if (!grepl("^SCT", assay, ignore.case = TRUE)) {
      Etest <- tryCatch(Seurat::GetAssayData(object, assay = assay, slot = "data")[1:1, 1:1, drop = FALSE],
                        error = function(e) NULL)
      if (is.null(Etest) || all(Etest == 0)) {
        if (verbose) message("[ScoreModules] Normalizing assay '", assay, "' for 'data' slot...")
        object <- Seurat::NormalizeData(object, assay = assay, verbose = FALSE)
      }
    }
  }
  
  # ---- expression matrix & overlap ----
  E_data    <- Seurat::GetAssayData(object, assay = assay, slot = slot_use)  # genes x cells
  gene_pool <- rownames(E_data)
  
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
  
  # ---- choose method (auto, UCell first) ----
  `%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a
  if (method == "auto") {
    avail <- function(m){
      switch(m,
             "UCell" = requireNamespace("UCell", quietly = TRUE),
             "Seurat"= TRUE,
             "zmean" = TRUE,
             "ssGSEA"= requireNamespace("GSVA", quietly = TRUE),
             "gsva"  = requireNamespace("GSVA", quietly = TRUE),
             "AUCell"= requireNamespace("AUCell", quietly = TRUE),
             FALSE
      )
    }
    pick <- NULL
    for (m in auto_order) if (avail(m)) { pick <- m; break }
    method <- pick %||% "zmean"
    if (verbose) message("[ScoreModules] auto -> using method: ", method)
  }
  
  ScoreMat <- NULL
  
  # ---- Seurat::AddModuleScore (with naming fix) ----
  if (method == "Seurat") {
    if (verbose) message("[ScoreModules] Seurat::AddModuleScore (", length(features_found), " sets)")
    object <- Seurat::AddModuleScore(
      object, features = features_found, name = paste0(name_prefix, "_"),
      search = TRUE, assay = assay
    )
    base <- paste0(name_prefix, "_", names(features_found))        # target names: Score_<name>
    cand_named1 <- paste0(base, "1")                               # Score_<name>1
    cand_named  <- base                                            # Score_<name>
    cand_num    <- paste0(name_prefix, "_", seq_along(features_found))  # Score_1..n
    
    if      (all(cand_named1 %in% colnames(object@meta.data))) cols_in <- cand_named1
    else if (all(cand_named  %in% colnames(object@meta.data))) cols_in <- cand_named
    else if (all(cand_num    %in% colnames(object@meta.data))) cols_in <- cand_num
    else stop("[ScoreModules] Cannot locate AddModuleScore outputs (tried name, name1, numeric).")
    
    # unify names to Score_<name>
    colnames(object@meta.data)[ match(cols_in, colnames(object@meta.data)) ] <- base
    ScoreMat <- as.matrix(object@meta.data[, base, drop = FALSE])
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
    colnames(object@meta.data)[ match(uc, colnames(object@meta.data)) ] <- out
    ScoreMat <- as.matrix(object@meta.data[, out, drop = FALSE])
  }
  
  # ---- z-mean (per-gene z then mean) ----
  if (method == "zmean") {
    if (!requireNamespace("matrixStats", quietly = TRUE))
      stop("[ScoreModules] Requires 'matrixStats' for zmean.")
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
  
  # ---- GSVA / ssGSEA ----
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
      ssgsea.norm    = TRUE,
      mx.diff        = TRUE,
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
  
  # ---- write to meta if needed (branches that didn't already) ----
  if (add_to_meta && !(method %in% c("Seurat","UCell","ssGSEA","gsva"))) {
    for (j in colnames(ScoreMat)) object[[j]] <- as.numeric(ScoreMat[, j])
  }
  
  invisible(list(object = object, scores = if (return_matrix) ScoreMat else NULL,
                 method = method, assay = assay))
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





pretty_violin_nature <- function(
    obj,
    features,                 # 需要绘制的分数列名，如 c("Score_Tumor5", ...)
    group_by,                 # 分组列名
    order        = NULL,
    title        = NULL,
    ylab         = "z-score",
    pt_size      = 0,
    add_box      = TRUE,
    add_median   = TRUE,
    add_zero_line= TRUE,
    zero_line    = 0,
    palette      = NULL,      # 命名向量：c(level="#color", ...)
    same_ylim    = FALSE,
    ylim         = NULL,
    ncol         = NULL,
    nrow         = NULL,
    design       = NULL,
    base_size    = 10,
    axis_lwd     = 0.6,
    # --- 新增：控制图例与合并 ---
    legend_position = "top",
    legend_title    = NULL,   # 默认用分组列名
    collect_legend  = TRUE
){
  stopifnot(inherits(obj, "Seurat"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 not available")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork not available")
  if (!requireNamespace("scales", quietly = TRUE)) stop("scales not available")
  
  `%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a
  
  # -- 分组列为因子 --
  grp <- obj@meta.data[[group_by]]
  if (is.data.frame(grp)) grp <- grp[[1, drop=TRUE]]
  grp <- as.vector(grp)
  if (!is.null(order)) grp <- factor(grp, levels = order) else grp <- factor(grp, levels = unique(grp))
  obj@meta.data[[group_by]] <- grp
  
  # -- 调色板 --
  levs <- levels(grp)
  pal  <- palette %||% setNames(scales::hue_pal()(length(levs)), levs)
  leg_title <- legend_title %||% group_by
  
  # -- 统一 y 轴（可选） --
  if (isTRUE(same_ylim) && is.null(ylim)) {
    vals <- Seurat::FetchData(obj, vars = features)
    r <- range(as.matrix(vals), na.rm = TRUE); pad <- 0.05 * diff(r)
    ylim <- c(r[1] - pad, r[2] + pad)
  }
  
  # -- 逐特征的小提琴图（保留图例，用于后面合并） --
  plist <- Seurat::VlnPlot(
    obj, features = features, group.by = group_by,
    pt.size = pt_size, combine = FALSE
  )
  
  pretty_list <- lapply(seq_along(plist), function(i){
    nm <- features[i]
    nm_clean <- sub("^Score_", "", nm)
    g <- plist[[i]] +
      ggplot2::scale_fill_manual(values = pal, name = leg_title, guide = "legend") +
      ggplot2::labs(title = nm_clean, x = NULL, y = ylab) +
      ggplot2::theme_classic(base_size = base_size) +
      ggplot2::theme(
        # —— X 轴完全隐藏（按你要求）
        axis.title.x = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        # —— 无网格线
        panel.grid = ggplot2::element_blank(),
        # —— 轴线
        axis.line  = ggplot2::element_line(color = "black", linewidth = axis_lwd),
        axis.ticks = ggplot2::element_line(color = "black", linewidth = axis_lwd*0.8),
        # —— 标题与图例
        plot.title = ggplot2::element_text(face = "bold", size = base_size, hjust = 0),
        legend.position = legend_position
      )
    if (isTRUE(add_box)) {
      g <- g + ggplot2::geom_boxplot(
        width = 0.12, outlier.shape = NA, alpha = 0.35,
        color = "black", position = ggplot2::position_dodge(width = 0.9), show.legend = FALSE
      )
    }
    if (isTRUE(add_median)) {
      g <- g + ggplot2::stat_summary(fun = median, geom = "point", size = 1, color = "black", show.legend = FALSE)
    }
    if (isTRUE(add_zero_line)) {
      g <- g + ggplot2::geom_hline(yintercept = zero_line, linetype = "dashed",
                                   linewidth = 0.3, color = "grey30")
    }
    if (!is.null(ylim)) g <- g + ggplot2::coord_cartesian(ylim = ylim)
    g
  })
  
  # -- 布局 & 合并图例 --
  if (is.null(design)) {
    ncol <- ncol %||% length(features)
    nrow <- nrow %||% 1
    p <- patchwork::wrap_plots(pretty_list, ncol = ncol, nrow = nrow)
  } else {
    p <- patchwork::wrap_plots(pretty_list, design = design)
  }
  if (isTRUE(collect_legend)) {
    p <- p + patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = legend_position)
  }
  if (!is.null(title) && nchar(title) > 0) {
    p <- p + patchwork::plot_annotation(title = title)
  }
  return(p)
}






axis_arrow_theme <- function(len_cm = 0.10, linewidth = 0.3, colour = "black",
                             ends = "last", type = "open") {
  theme(
    axis.line.x = element_line(
      colour = colour, linewidth = linewidth,
      arrow  = arrow(length = unit(len_cm, "cm"), ends = ends, type = type)
    ),
    axis.line.y = element_line(
      colour = colour, linewidth = linewidth,
      arrow  = arrow(length = unit(len_cm, "cm"), ends = ends, type = type)
    )
  )
}





suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(ggh4x); library(dplyr); library(scales)
})

plot_umap <- function(
    obj,
    cluster_col,
    reduction    = "banksy_umap",
    facet_col    = NULL,
    colors       = NULL,
    out_prefix   = "UMAP",
    out_dir      = getwd(),
    width        = 2, height = 2, dpi = 300, bg = "white",
    legend_ncol  = 2,
    strip_fill   = "#e6bac5",
    # 箭头参数
    arrow_len_cm = 0.10, arrow_lwd = 0.3, arrow_color = "black",
    arrow_type   = "open", arrow_ends = "last"
){
  stopifnot(reduction %in% names(obj@reductions))
  emb <- Embeddings(obj, reduction); stopifnot(ncol(emb) >= 2)
  
  md <- obj@meta.data
  stopifnot(cluster_col %in% colnames(md))
  Cluster <- md[[cluster_col]]; if (is.data.frame(Cluster)) Cluster <- Cluster[[1, drop=TRUE]]
  Cluster <- factor(Cluster, levels = unique(Cluster))
  
  if (!is.null(facet_col)) {
    stopifnot(facet_col %in% colnames(md))
    label <- md[[facet_col]]; if (is.data.frame(label)) label <- label[[1, drop=TRUE]]
    label <- factor(label, levels = unique(label))
  } else {
    label <- factor("ALL")
  }
  
  df <- tibble::tibble(
    UMAP_1  = emb[,1],
    UMAP_2  = emb[,2],
    Cluster = Cluster,
    label   = label
  ) |> dplyr::filter(!is.na(Cluster), !is.na(UMAP_1), !is.na(UMAP_2))
  
  levs <- levels(df$Cluster)
  if (is.null(colors)) {
    colors <- setNames(hue_pal()(length(levs)), levs)
  } else {
    if (is.null(names(colors))) names(colors) <- levs
    if (!all(levs %in% names(colors))) {
      missing <- setdiff(levs, names(colors))
      colors  <- c(colors, setNames(hue_pal()(length(missing)), missing))
    }
    colors <- colors[levs]
  }
  
  xr <- range(df$UMAP_1, na.rm = TRUE)
  yr <- range(df$UMAP_2, na.rm = TRUE)
  
  # 主图（无图例）
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = 0.03, shape = 16, stroke = 0) +
    scale_color_manual(values = colors) +
    facet_grid(~label) +
    theme_classic() +
    theme(
      plot.background  = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(0,0,0,0),
      axis.title.x     = element_blank(),
      axis.title.y     = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.line.x      = element_line(
        colour = arrow_color, linewidth = arrow_lwd,
        arrow  = arrow(length = unit(arrow_len_cm, "cm"),
                       type = arrow_type, ends = arrow_ends)
      ),
      axis.line.y      = element_line(
        colour = arrow_color, linewidth = arrow_lwd,
        arrow  = arrow(length = unit(arrow_len_cm, "cm"),
                       type = arrow_type, ends = arrow_ends)
      ),
      strip.background = element_rect(fill = strip_fill, color = NA),
      strip.placement  = "outside",
      strip.text       = element_text(size = 8),
      legend.position  = "none",
      aspect.ratio     = 1
    ) +
    scale_x_continuous(limits = xr) +
    scale_y_continuous(limits = yr)
  
  # 单独图例
  p_leg <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = 0, shape = 16, stroke = 0) +
    theme_void() +
    scale_color_manual(values = colors, name = "") +
    theme(
      aspect.ratio      = 1,
      legend.position   = "bottom",
      plot.margin       = margin(0,0,0,0),
      legend.text       = element_text(size = 5),
      legend.spacing.y  = unit(0, "cm"),
      legend.key.height = unit(0, "cm"),
      legend.box.spacing= unit(0, "cm")
    ) +
    guides(color = guide_legend(ncol = legend_ncol,
                                override.aes = list(size = 2, alpha = 1)))
  
  # 主图 + 图例
  p_with_leg <- p +
    theme(
      legend.position   = "bottom",
      legend.text       = element_text(size = 5),
      legend.spacing.y  = unit(0, "cm"),
      legend.key.height = unit(0, "cm"),
      legend.box.spacing= unit(0, "cm")
    ) +
    guides(color = guide_legend(ncol = legend_ncol,
                                override.aes = list(size = 2, alpha = 1)))
  
  # 保存
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(out_dir, paste0(out_prefix, ".png")),
         p, width = width, height = height, dpi = dpi, bg = bg)
  ggsave(file.path(out_dir, paste0(out_prefix, "_legend.pdf")),
         p_leg, width = 3.5, height = 6, dpi = 500, bg = "white")
  ggsave(file.path(out_dir, paste0(out_prefix, "_withLegend.png")),
         p_with_leg, width = width, height = height + 0.6, dpi = dpi, bg = bg)
  
  invisible(list(main = p, legend = p_leg, main_with_legend = p_with_leg,
                 data = df, palette = colors))
}




suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr)
  library(viridis)      # 色标
  # 可选：library(mascarade) # 画圈需要；没装也能出星云主图
})

# —— Nature风黑底主题（精简版，和论文风一致）——
galaxyTheme_black <- function(base_size = 12) {
  theme_void(base_size = base_size) %+replace%
    theme(
      plot.background  = element_rect(color = "black", fill = "black"),
      panel.background = element_rect(fill  = "black", color = NA),
      panel.border     = element_blank(),
      legend.position  = "none",
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text       = element_text(size = base_size*0.8, color = "white"),
      plot.margin      = margin(5.5, 5.5, 5.5, 5.5)
    )
}





# ------------------------------------------------------------
# galaxy + mask（圈边界）UMAP
# ------------------------------------------------------------
plot_galaxy_umap <- function(
    obj,
    cluster_col,                     # 用于加圈的分组列（如 celltype / cluster）
    reduction     = "banksy_umap",   # 也可 "umap" / "tsne"
    facet_col     = NULL,            # 可选：上下分面（如 "Field1": tumor/NL）；NULL=单面板
    facet_ncol    = 1,               # 分面列数（论文例是上下 => 1 列）
    annotate_df   = NULL,            # 可选：数据框 data.frame(x, y, label) 添加文字
    show_points   = TRUE,            # 叠加白色微点
    point_size    = 0.02,
    point_alpha   = 1.0,
    # 星云密度层
    density_geom  = c("raster","polygon"), # "raster" 论文风；"polygon"也可
    density_bins  = 0,               # raster时忽略；polygon时可设等值面数
    density_bw    = NULL,            # NULL=自动；可传 numeric 控制核密度带宽
    density_option= "magma",         # viridis 方案：magma/inferno/plasma…
    # 圈（mask）相关
    add_mask      = TRUE,            # TRUE 需安装 mascarade 包
    minDensity    = 15,              # 圈的松紧
    smoothSigma   = 0.10,            # 圈的平滑度
    mask_lwd      = 0.6,
    mask_linetype = 2,
    mask_color    = "white",
    # 输出
    out_dir       = getwd(),
    out_prefix    = "GalaxyUMAP",
    width         = 8, height = 10, dpi = 300, bg = "black"  # 论文推荐大图导出
){
  stopifnot(inherits(obj, "Seurat"))
  if (!reduction %in% names(obj@reductions))
    stop("Reduction '", reduction, "' not found in obj@reductions.")
  
  emb <- Embeddings(obj, reduction)
  stopifnot(ncol(emb) >= 2)
  
  md <- obj@meta.data
  if (!cluster_col %in% colnames(md))
    stop("cluster_col '", cluster_col, "' not in meta.data.")
  
  Cluster <- md[[cluster_col]]
  if (is.data.frame(Cluster)) Cluster <- Cluster[[1, drop = TRUE]]
  Cluster <- factor(Cluster, levels = unique(Cluster))
  
  if (!is.null(facet_col)) {
    if (!facet_col %in% colnames(md))
      stop("facet_col '", facet_col, "' not in meta.data.")
    facet_lab <- md[[facet_col]]
    if (is.data.frame(facet_lab)) facet_lab <- facet_lab[[1, drop = TRUE]]
    facet_lab <- factor(facet_lab, levels = unique(facet_lab))
  } else {
    facet_lab <- factor("ALL")
  }
  
  df <- tibble::tibble(
    UMAP_1  = emb[,1],
    UMAP_2  = emb[,2],
    Cluster = Cluster,
    label   = facet_lab
  ) |>
    dplyr::filter(!is.na(UMAP_1), !is.na(UMAP_2), !is.na(Cluster))
  
  # —— 星云底图（与论文一致：stat_density_2d + viridis）——
  density_geom <- match.arg(density_geom)
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
    {
      if (density_geom == "raster")
        stat_density_2d(aes(fill = after_stat(density)),
                        geom = "raster", contour = FALSE, bw = density_bw)
      else
        stat_density_2d_filled(aes(fill = after_stat(level)),
                               bins = ifelse(density_bins > 0, density_bins, 60))
    } +
    scale_fill_viridis(option = density_option, guide = "none")
  
  # 白色微点（增强星尘感）
  if (isTRUE(show_points)) {
    p <- p + geom_point(color = "white", size = point_size, alpha = point_alpha)
  }
  
  # 添加注释文字（可选）
  if (!is.null(annotate_df)) {
    stopifnot(all(c("x","y","label") %in% colnames(annotate_df)))
    p <- p + geom_text(
      data = annotate_df, aes(x = x, y = y, label = label),
      inherit.aes = FALSE, color = "white", size = 6
    )
  }
  
  # 分面：论文示例上下两面板
  p <- p + facet_wrap(~label, ncol = facet_ncol)
  
  # 主题：纯黑背景 + 去网格 + 去坐标轴（与论文一致）
  p <- p + galaxyTheme_black()
  
  # —— 可选：圈（mask）——
  maskTable <- NULL
  if (isTRUE(add_mask)) {
    if (!requireNamespace("mascarade", quietly = TRUE)) {
      message("[plot_galaxy_umap] Package 'mascarade' not installed; skip mask.")
    } else {
      # 注意：mask 以“全局”坐标和 cluster 计算；若你要按分面分别圈，可按 label 分组分别生成后 rbind
      maskTable <- mascarade::generateMask(
        dims       = as.matrix(df[, c("UMAP_1","UMAP_2")]),
        cluster    = df$Cluster,
        minDensity = minDensity,
        smoothSigma= smoothSigma
      )
      p <- p + geom_path(
        data = maskTable,
        aes(x = x, y = y, group = group),
        inherit.aes = FALSE,
        linewidth = mask_lwd, linetype = mask_linetype, colour = mask_color
      )
    }
  }
  
  # —— 保存 ——（论文建议大尺寸导出）
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_png <- file.path(out_dir, paste0(out_prefix, ".png"))
  out_pdf <- file.path(out_dir, paste0(out_prefix, ".pdf"))
  ggsave(out_png, p, width = width, height = height, dpi = dpi, bg = bg)
  ggsave(out_pdf, p, width = width, height = height, dpi = dpi, bg = bg)
  
  invisible(list(plot = p, data = df, mask = maskTable,
                 files = list(png = out_png, pdf = out_pdf)))
}





# ============================================================
# Paper-style FeaturePlot grid (updated)
# - Adds panel_tag (corner labels), flexible titles, robust checks
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# Paper palette from the example (continuous)
paper_continuous_colors <- c(
  "#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080"
)

# ---- helper: feature existence check ----
.safe_features <- function(obj, feats){
  feats <- unique(as.character(feats))
  ok   <- feats[feats %in% rownames(obj)]
  miss <- setdiff(feats, ok)
  if (length(miss) > 0) {
    message(sprintf("[WARN] %d features not found: %s",
                    length(miss), paste(miss, collapse = ", ")))
  }
  ok
}





# ---- helper: single FeaturePlot with paper theme ----
.single_featureplot_paper <- function(
    obj, feature, reduction = "umap",
    cols = paper_continuous_colors,
    pt.size = 0.5, max.cutoff = 1.5,
    title = NULL, order_points = TRUE
){
  p <- FeaturePlot(
    obj, features = feature, reduction = reduction,
    pt.size = pt.size, max.cutoff = max.cutoff,
    cols = cols, order = order_points
  ) +
    scale_x_continuous("") + scale_y_continuous("") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks       = element_blank(),
      axis.text        = element_blank(),
      legend.position  = "none",
      plot.title       = element_text(hjust = 0.5, size = 10)
    )
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}





# ---- main: featureplot_grid_paper (UPDATED) ----
# markers: named character vector or named list
#   - names are panel titles (cell types), values are gene symbols
#   - if no names provided, titles fall back to gene symbols
# ===============================
# Paper-style FeaturePlot grid
# ===============================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# ---- 0) Default palette from the paper snippet ----
paper_continuous_colors <- c("#FFEFD5","#E6E6FA","#87CEFA","#6495ED","#4169E1","#0000CD","#000080")

# ---- 1) Helper: safe feature check ----
.safe_features <- function(obj, feats){
  feats <- unique(feats)
  ok <- feats[feats %in% rownames(obj)]
  miss <- setdiff(feats, ok)
  if (length(miss) > 0) {
    message(sprintf("[WARN] %d features not found: %s",
                    length(miss), paste(miss, collapse = ", ")))
  }
  ok
}





# ---- 2) Helper: a single FeaturePlot with the paper theme ----
.single_featureplot_paper <- function(
    obj, feature, reduction = "umap",
    cols = paper_continuous_colors,
    pt.size = 0.5, max.cutoff = 1.5,
    title = NULL
){
  p <- FeaturePlot(
    obj, features = feature, reduction = reduction,
    pt.size = pt.size, max.cutoff = max.cutoff,
    cols = cols, order = TRUE,label = F
  ) +
    scale_x_continuous("") + scale_y_continuous("") +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks       = element_blank(),
      axis.text        = element_blank(),
      legend.position  = "none",
      plot.title       = element_text(hjust = 0.5, size = 10)
    )
  
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}





# ---- 3) Main: featureplot_grid_paper() ----
# markers: named character vector or a named list.
# - If named vector: names are "panel titles", values are gene symbols.
# - If list: names are "panel titles", values are length-1 vectors with gene symbols (for clarity).
# Updated featureplot_grid_paper() with `title_fmt`
featureplot_grid_paper <- function(
    obj,
    markers,
    reduction = "umap",
    cols = paper_continuous_colors,
    pt.size = 0.5,
    max.cutoff = 1.5,
    ncol = 4,
    strip_titles = TRUE,           # show titles or not
    title_fmt   = "%s (%s)",       # sprintf format: first %s = panel title, second %s = gene
    width = 12, height = 6,
    outfile = NULL
){
  # Normalize "markers" to a named character vector
  if (is.list(markers)) {
    stopifnot(!is.null(names(markers)))
    markers_vec <- vapply(markers, function(x) x[[1]], FUN.VALUE = character(1))
    names(markers_vec) <- names(markers)
  } else {
    markers_vec <- markers
  }
  if (is.null(names(markers_vec))) {
    # if not named, use gene symbol itself as panel title
    names(markers_vec) <- markers_vec
  }
  
  # Check existence
  ok_genes <- .safe_features(obj, unname(markers_vec))
  if (length(ok_genes) == 0L) stop("None of the provided features were found in the object.")
  
  # Drop missing features but keep order of panels
  keep_panels <- names(markers_vec)[markers_vec %in% ok_genes]
  markers_vec <- markers_vec[keep_panels]
  
  # Build plots
  plots <- vector("list", length(markers_vec))
  for (i in seq_along(markers_vec)) {
    gene  <- markers_vec[[i]]
    label <- names(markers_vec)[i]
    ttl   <- if (strip_titles) sprintf(title_fmt, label, gene) else NULL
    
    plots[[i]] <- .single_featureplot_paper(
      obj, feature = gene, reduction = reduction,
      cols = cols, pt.size = pt.size, max.cutoff = max.cutoff,
      title = ttl
    )
  }
  
  # Assemble grid
  grid <- wrap_plots(plots, ncol = ncol)
  
  # Export if needed
  if (!is.null(outfile)) {
    ggsave(filename = outfile, plot = grid,
           width = width, height = height, units = "in", dpi = 300, bg = "white")
  }
  grid
}





# ---- 4) Optional: legend-only panel for the continuous scale ----
# Works by drawing an "empty" FeaturePlot to extract the continuous color bar if needed later.
featureplot_continuous_legend <- function(
    obj, feature,
    reduction = "umap",
    cols = paper_continuous_colors,
    pt.size = 0.0, max.cutoff = 1.5
){
  FeaturePlot(
    obj, features = feature, reduction = reduction,
    pt.size = pt.size, max.cutoff = max.cutoff,
    cols = cols
  ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
}







suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
})

# ---------- helper: 清洗基因集（大小写不敏感，缺失基因给建议） ----------
.sanitize_markers <- function(obj, marker_sets, case_insensitive = TRUE) {
  stopifnot(inherits(obj, "Seurat"))
  gp <- rownames(obj)
  gp_up <- toupper(gp)
  
  suggest_one <- function(sym) {
    # 用 agrep 给个相似建议（最多1个）
    idx <- tryCatch(agrep(sym, gp, max.distance = 1, value = FALSE), error = function(e) integer(0))
    if (length(idx)) gp[idx[1]] else NA_character_
  }
  
  out_sets <- list()
  report <- list()
  for (nm in names(marker_sets)) {
    v <- unique(as.character(marker_sets[[nm]]))
    found <- character(0); miss <- character(0); map <- character(0)
    for (g in v) {
      if (g %in% gp) {
        found <- c(found, g); map <- c(map, g)
      } else if (isTRUE(case_insensitive) && toupper(g) %in% gp_up) {
        g_mapped <- gp[ match(toupper(g), gp_up) ]
        found <- c(found, g_mapped); map <- c(map, g_mapped)
      } else {
        miss <- c(miss, g); map <- c(map, NA_character_)
      }
    }
    out_sets[[nm]] <- found
    # 生成报告
    rep_df <- tibble::tibble(
      set = nm,
      input_gene = v,
      mapped_gene = map,
      status = ifelse(!is.na(map), "found", "missing"),
      suggestion = ifelse(is.na(map), vapply(v, suggest_one, character(1)), NA_character_)
    )
    report[[nm]] <- rep_df
  }
  list(marker_sets = out_sets, report = dplyr::bind_rows(report))
}






library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

plot_marker_dot <- function(
    obj,
    marker_sets,
    group_by,
    # --- 分组参数 ---
    grouping_map = NULL,
    group_order  = NULL,
    group_style  = c("color_legend","label","none"),
    # --- 排序参数 ---
    group_sort   = c("none","alpha","numeric","freq","custom"),
    # --- 颜色 ---
    group_colors = NULL,
    # --- 图例 ---
    legend_position = "right",
    legend_box      = "vertical",
    # --- 原有参数 ---
    assay      = NULL,
    cols       = c("lightgrey","red"),
    dot_scale  = 6,
    scale_min  = 0,
    scale_max  = NA,
    strip_fill = "#e6bac5",
    strip_size = 9,
    gene_text_size  = 8,
    group_text_size = 9,
    x_axis_angle    = 45,
    case_insensitive = TRUE,
    drop_empty_sets  = TRUE
){
  group_style <- match.arg(group_style)
  group_sort  <- match.arg(group_sort)
  stopifnot(inherits(obj, "Seurat"))
  if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)
  if (!group_by %in% colnames(obj@meta.data))
    stop("group_by '", group_by, "' not in meta.data")
  
  # ---------- 1) 规范 group levels ----------
  gb <- obj@meta.data[[group_by]]
  if (!is.factor(gb)) gb <- factor(gb)
  
  if (is.null(grouping_map)) {
    # —— 无上层分组：对 group_by levels 排序/重排
    current_levels <- levels(gb)
    if (!is.null(group_order) && length(group_order) > 0) {
      valid <- intersect(group_order, current_levels)
      missing <- setdiff(current_levels, valid)
      gb <- factor(gb, levels = c(valid, missing))
    } else {
      if (group_sort == "alpha") {
        gb <- factor(gb, levels = sort(levels(gb)))
      } else if (group_sort == "numeric") {
        num <- suppressWarnings(as.numeric(levels(gb)))
        if (all(is.na(num))) stop("Cannot sort by numeric: levels not numeric.")
        gb <- factor(gb, levels = levels(gb)[order(num)])
      } else if (group_sort == "freq") {
        fr <- table(gb)
        gb <- factor(gb, levels = names(sort(fr, decreasing = TRUE)))
      } # "none": 保持不变
    }
    obj@meta.data[[group_by]] <- gb
    if (group_style != "none") group_style <- "none"
  } else {
    # —— 有上层分组：构造 group 列
    grp_col <- paste0(group_by, "_group")
    gm <- as.character(grouping_map[as.character(gb)])
    gm[is.na(gm)] <- "Unassigned"
    obj@meta.data[[grp_col]] <- gm
    if (is.null(group_order)) {
      group_order <- sort(unique(obj@meta.data[[grp_col]]))
    } else {
      group_order <- intersect(group_order, unique(obj@meta.data[[grp_col]]))
      if (length(group_order) == 0) group_order <- sort(unique(obj@meta.data[[grp_col]]))
    }
    meta_df <- obj@meta.data |>
      dplyr::mutate(.level = as.character(.data[[group_by]])) |>
      dplyr::distinct(.level, .keep_all = TRUE) |>
      dplyr::arrange(factor(.data[[grp_col]], levels = group_order),
                     .level)
    obj@meta.data[[group_by]] <- factor(obj@meta.data[[group_by]], levels = meta_df$.level)
    gb <- obj@meta.data[[group_by]]
  }
  
  # ---------- 2) 处理基因集 ----------
  all_genes <- rownames(obj[[assay]])
  if (isTRUE(case_insensitive)) {
    lu <- setNames(all_genes, tolower(all_genes))
    ms_clean <- lapply(marker_sets, function(g) {
      out <- unname(lu[tolower(g)])
      out[!is.na(out)]
    })
  } else {
    ms_clean <- lapply(marker_sets, function(g) intersect(g, all_genes))
  }
  if (isTRUE(drop_empty_sets)) {
    keep <- vapply(ms_clean, length, integer(1)) > 0
    if (!all(keep)) message("Dropped empty gene sets: ", paste(names(ms_clean)[!keep], collapse = ", "))
    ms_clean <- ms_clean[keep]
  }
  if (length(ms_clean) == 0) stop("No valid marker genes found in provided sets.")
  
  set_names <- names(ms_clean)
  gene2set  <- stack(ms_clean) |> dplyr::rename(gene=values, set=ind)
  gene2set$set <- factor(gene2set$set, levels = set_names)
  
  # ---------- 3) 取表达矩阵并统计 avg / pct ----------
  M <- Seurat::GetAssayData(obj, assay = assay, slot = "data")
  feats <- unique(gene2set$gene)
  M <- M[feats, , drop = FALSE]
  
  groups <- obj@meta.data[[group_by]]
  idx_by_group <- split(seq_len(ncol(M)), groups)
  
  avg_list <- lapply(names(idx_by_group), function(gid){
    ix <- idx_by_group[[gid]]
    if (length(ix) == 0) rep(NA_real_, nrow(M)) else Matrix::rowMeans(M[, ix, drop = FALSE])
  })
  names(avg_list) <- names(idx_by_group)
  
  pct_list <- lapply(names(idx_by_group), function(gid){
    ix <- idx_by_group[[gid]]
    if (length(ix) == 0) rep(NA_real_, nrow(M)) else Matrix::rowMeans(M[, ix, drop = FALSE] > 0) * 100
  })
  names(pct_list) <- names(idx_by_group)
  
  df <- dplyr::bind_rows(lapply(names(idx_by_group), function(gid){
    tibble::tibble(
      gene  = rownames(M),
      id    = gid,                    # 细胞/亚群
      avg   = as.numeric(avg_list[[gid]]),
      pct   = as.numeric(pct_list[[gid]])
    )
  }))
  
  # ---------- 4) 合并 set 信息并缩放 avg ----------
  df <- df |>
    dplyr::left_join(gene2set, by = "gene") |>
    dplyr::group_by(gene) |>
    dplyr::mutate(avg_scaled = as.numeric(scale(avg))) |>
    dplyr::ungroup()
  
  if (!is.na(scale_min)) df$avg_scaled <- pmax(df$avg_scaled, scale_min)
  if (!is.na(scale_max)) df$avg_scaled <- pmin(df$avg_scaled, scale_max)
  
  # 基因按输入顺序；Y 轴为细胞/亚群
  gene_order_by_set <- unlist(ms_clean[set_names], use.names = FALSE)
  df$gene <- factor(df$gene, levels = gene_order_by_set)           # —— X 轴 = 基因
  df$id   <- factor(df$id,   levels = levels(groups))               # —— Y 轴 = 细胞/亚群
  
  # ---------- 5) 作图：X=基因, Y=细胞；无网格线 ----------
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene, y = id)) +
    ggplot2::geom_point(ggplot2::aes(size = pct, colour = avg_scaled)) +
    ggplot2::scale_size(range = c(0, dot_scale), name = "% expressed") +
    ggplot2::scale_colour_gradient(low = cols[1], high = cols[2], name = "Scaled avg") +
    ggplot2::facet_grid(cols = ggplot2::vars(set), scales = "free_x", space = "free_x") +
    ggplot2::theme_classic() +                                    # —— 无网格线
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_text(size = gene_text_size, angle = x_axis_angle, hjust = 1, vjust = 1),
      axis.text.y  = ggplot2::element_text(size = group_text_size),
      strip.background = ggplot2::element_rect(fill = strip_fill, color = NA),
      strip.text = ggplot2::element_text(size = strip_size, face = "bold", colour = "white"),
      panel.grid.major = ggplot2::element_blank(),                # 双保险去网格
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = legend_position,
      legend.box       = legend_box
    )
  
  # ---------- 6) 分组风格（可选信息展示） ----------
  if (!is.null(grouping_map)) {
    if (group_style == "label") {
      grp_col <- paste0(group_by, "_group")
      gtxt <- paste("Groups:", paste(unique(obj@meta.data[[grp_col]]), collapse = " · "))
      p <- p + ggplot2::labs(subtitle = gtxt)
    }
    # 注：逐个 y 轴标签上色需要 grob 操作，易碎；此处保持简洁稳健
  }
  
  return(p)
}




# ===============================================================
# Raghavan Cell 2021 - Supplemental Table S3 (top30)
# Tumor-only: scBasal / scClassical / IC
# Uses your existing plot_marker_dot()
# ===============================================================

# ---- 跨集合去重：避免 DotPlot 因重复基因名报错 ----
.hd_dedup_sets <- function(marker_sets, priority = names(marker_sets)) {
  stopifnot(is.list(marker_sets))
  seen <- character(0); out <- list()
  for (nm in priority) {
    feats <- unique(marker_sets[[nm]])
    keep  <- setdiff(feats, seen)
    out[[nm]] <- keep
    seen <- c(seen, keep)
  }
  out
}




# ===============================================================
# HD PDAC Marker Panel (All Major Classes)
# Integrates Raghavan Cell 2021 S3 (Tumor subtypes) + Comprehensive PDAC Markers
# From high-impact refs: Peng Cell 2019, Steele NatMed 2020, Hwang NatGen 2022, etc.
# Uses your existing plot_marker_dot()
# ===============================================================

# ===============================================================
# HD PDAC Marker Panel (All Major Classes)
# Integrates Raghavan Cell 2021 S3 (Tumor subtypes) + Comprehensive PDAC Markers
# From high-impact refs: Peng Cell 2019, Steele NatMed 2020, Hwang NatGen 2022, etc.
# ===============================================================

# ---- 跨集合去重：避免 DotPlot 因重复基因名报错 ----
.hd_dedup_sets <- function(marker_sets, priority = names(marker_sets)) {
  stopifnot(is.list(marker_sets))
  seen <- character(0); out <- list()
  for (nm in priority) {
    feats <- unique(marker_sets[[nm]])
    keep  <- setdiff(feats, seen)
    out[[nm]] <- keep
    seen <- c(seen, keep)
  }
  out
}





# ---- 只返回 PDAC 所有大类的签名（core/extended）----
# ===============================================================
# HD PDAC Marker Sets (All Classes)
# panel = "core" (gold anchors, 3-5 genes/class) or "full" (extended, 10-20 genes/class)
# ===============================================================
hd_get_marker_sets_PDAC <- function(
    panel = c("core", "full"),
    dedup_across_sets = TRUE,
    add_registered_extras = TRUE,
    include_emt = TRUE,                 # EMT anchor can be added to either panel
    emt_name   = "EMT_anchor"           # display name of EMT anchor set
) {
  panel <- match.arg(panel)
  
  # --- 1. Tumor (Classical/Basal-like) ---
  Tumor_core <- c("EPCAM", "KRT8", "KRT18", "CEACAM5")
  Tumor_full <- c("KRT19","CLDN4","TACSTD2","LCN2","MSLN","AGR2","MUC1","MUC5AC","TFF1","REG4",
                  "S100A2","KRT7","SOX9","MMP7","SLC4A4","GATA6","TP63","VIM","FN1")
  
  # Subtype panels
  scBasal_core <- c("KRT17","S100A2","KRT6A","KRT5","LY6D","SPRR1B","LAMC2","HMGA2")
  scBasal_full <- c("KRT6A","S100A2","KRT13","KRT17","LY6D","KRT7","KRT5","FAM83A","CD109","GAPDH",
                    "CSTA","PTPN13","MTSS1","SPRR1B","SLC2A1","CAV1","EMP1","SEMA3C","DHRS9","CAST",
                    "FLNA","SLITRK6","COL7A1","AQP3","MT2A","AHNAK2","KRT19","PKP1","YBX3","GJB5")
  scClassical_core <- c("GATA6","CTSE","TFF1","AGR2","CEACAM6","REG4","MUC1","LGALS4","TSPAN8")
  scClassical_full <- c("LGALS4","CTSE","TFF1","AGR2","TSPAN8","CEACAM6","CDH17","TFF2","CLDN18","REG4",
                        "SPINK1","SLC44A4","GPX2","ANXA10","FAM3D","TFF3","MUC13","ANXA13","CEACAM5",
                        "VSIG2","PRSS3","ALDH2","IL1R2","AGR3","CYSTM1","CAPN5","S100A6","ONECUT2","BCAS1")
  IC_core  <- c("VCAM1","IGFBP3","TNFAIP2","SERPINA1","CXCL14","KLF6","TIMP2")
  IC_full  <- c("VCAM1","MALAT1","IGFBP3","SYNE2","TIMP2","TNFAIP2","NFAT5","SERPINA1","CXCL14",
                "CX3CL1","ALPK3","ATHL1","DPYSL2","TMEM139","KLF6","KRT23","MEGF6","HOXB2","PADI1",
                "MFSD4","INPP4B","BAIAP3","HMHA1","ECE1","CAMK2N1","SNORD89","WFDC2","TSHZ2","RILPL2","LAMB2")
  
  # --- 2. Ductal/ADM-like ---
  Ductal_core <- c("KRT19","SOX9","CFTR","HNF1B")
  Ductal_full <- c("KLF5","MMP7","SLC4A4","SOX9","KRT7","ANXA4","PIGR","CITED4","TSPAN8","CLDN18",
                   "PROM1","SPP1","OLFM4","LYZ","AGR3","REG1A","TFF3","KRT17","EPCAM","COL10A1","COL11A1","COL1A1")
  
  # --- 3. Acinar ---
  Acinar_core <- c("PRSS1","CPA1")
  Acinar_full <- c("CPB1","GP2","SYCN","PNLIP","CEL","CLPS","REG3A","AMY2A","PRSS2","CTRC",
                   "SPINK1","REG1B","BICC1","CELA3A","CELA2A","TRY4","CELA1","PRSS2","BHLHA15","PTF1A","RBPJL")
  
  # --- 4. Endocrine ---
  Endocrine_core <- c("CHGA","CHGB","INS")
  Endocrine_full <- c("GCG","SST","IAPP","PCSK1","PCSK2","NEUROD1","NKX6-1","MAFA","PDX1","PAX6",
                      "ARX","GHRL","PPY","SCG5","SCGN","SYT13")
  
  # --- 5. Endothelial ---
  Endothelial_core <- c("PECAM1","VWF","KDR")
  Endothelial_full <- c("FLT1","CLDN5","RAMP2","PLVAP","EGFL7","CDH5","ESAM","TIE1","ENG","ICAM2",
                        "ACKR1","FABP5","EMCN","SOX18","SOX17","CAV1","CD34")
  
  # --- 6. CAF (general) ---
  CAF_core <- c("COL1A1","DCN")
  CAF_full <- c("COL1A2","COL3A1","THY1","THBS2","POSTN","ITGA11","MMP11","FAP","PDGFRA","PDGFRB",
                "ACTA2","INHBA","MFAP5","SPARC","FN1","VCAN","COL6A1","COL6A2","COL6A3")
  
  # --- 6.1 CAF subtypes (NEW): iCAF / myCAF / apCAF ---
  # Core markers kept minimal; full adds supporting markers commonly reported in PDAC
  iCAF_core <- c(
    "CXCL12","PDGFRA","IL6","LIF","CCL2","CCL19","IGFBP3","MFAP5","COL14A1"
  )
  
  iCAF_full <- unique(c(
    iCAF_core,
    "CXCL14","HAS1","HAS2","PTGS2",
    "IGFBP4","IGFBP5","IGFBP7",
    "MFAP4","SLC40A1","APOD","OGN",
    "C1R","C1S","C3"
  ))
  
  ## ii) myCAF: contractile + ECM/fibrotic axis
  myCAF_core <- c(
    "ACTA2","TAGLN","MYL9","MYH11","CNN1","TPM2","COL1A1","COL1A2","THBS2","POSTN"
  )
  
  myCAF_full <- unique(c(
    myCAF_core,
    "SPARC","CCN2","MMP11","MMP14","ITGA11","LOXL1","LOXL2",
    "COL3A1","COL5A1","COL5A2","COL6A1","COL6A2",
    "TPM1","PALLD","CALD1","ACTN1","CSRP1","DPYSL3","ENAH"
  ))
  
  ## iii) apCAF: antigen-presenting CAF (MHC-II + CD74/CIITA)
  apCAF_core <- c(
    "CD74",     # invariant chain，不变链
    "CIITA",    # MHC-II 主调控
    "RFXAP",    # 同上
    "RFXANK",   # 同上
    "IFI30"     # GILT，MHC-II 抗原加工
  )
  
  # --- FULL：在 CORE 基础上扩展“抗原处理/IFNγ反应”轴；DR/DP/DQ 仅作可选项 ---
  apCAF_full <- unique(c(
    apCAF_core,   "RFX5",  "CTSS","CTSL","CTSB","LGMN",    # MHC-II 转录复合体 # 溶酶体蛋白酶/Legumain：辅助 Ii/CD74 加工与肽加载
    "IRF1","STAT1",               # IFNγ 通路：上游强化 CIITA/MHC-II
    # 可选：若你的对象里碰巧存在则会自动命中；没有也不影响
    "HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1"
  ))
  
  
  # --- 7. Perivascular / Smooth Muscle ---
  Perivascular_core <- c("RGS5","NOTCH3","PDGFRB")
  Perivascular_full <- c("MCAM","CSPG4","ACTA2","TAGLN","MYH11","MYL9","CNN1","TPM2","COL4A1","HIGD1B",
                         "ABCC9","KCNJ8","ADIRF","PTGIR","MYLK","DES","SORBS1")
  
  # --- 8. Immune (general) ---
  Immune_core <- c("PTPRC","HLA-DRA","AIF1")
  Immune_full <- c("LST1","HLA-DPB1","CTSS","AIF1","TYROBP","FCGR3A","LAPTM5","PSAP",
                   "C1QA","C1QB","C1QC","CD52","CD53","CD37","IFI30","IFITM3")
  
  # --- 9. T cell ---
  T_Cell_core <- c("CD3D","CD3E","TRAC")
  T_Cell_full <- c("CD3G","CD8A","CD4","IL7R","LEF1","SELL","CCR7","GZMA","GZMK","FOXP3","CTLA4",
                   "PDCD1","CXCL13","IFNG","GNLY","NKG7")
  
  # --- 10. B cell ---
  B_Cell_core <- c("MS4A1","CD79A","CD79B")
  B_Cell_full <- c("CD19","BANK1","CD22","FCRL1","IGHD","PAX5","CR2","CD40","BLNK","BTLA",
                   "CXCR5","STAP1","TCL1A","CD27","CD38","AICDA","BCL6")
  
  # --- 11. Plasma cell ---
  Plasma_core <- c("JCHAIN","IGKC","MZB1")
  Plasma_full <- c("IGHG1","XBP1","SDC1","DERL3","FKBP11","PRDM1","IRF4","SSR4","TXNDC5","SEL1L3",
                   "IGHA1","IGHG3","IGHM","CD27","CD38","SLAMF7")
  
  # --- 12. Macrophage (general) ---
  Macrophage_core <- c("CD68","TYROBP","CSF1R")
  Macrophage_full <- c("CD163","ITGAX","GPNMB","SPP1","MRC1","HLA-DQB1","C1QA","C1QB","APOE",
                       "TREM2","LILRB5","SIGLEC1","VSIG4","MS4A7","FCGR1A","ADGRE1","MARCO")
  # ... (在 Macrophage_full 定义之后)
  
  # --- 13. Neutrophil (中性粒细胞) ---  <-- (我给它加了个编号)
  # 【核心修正】: 添加 Neutrophil_core 的定义
  Neutrophil_core <- c("FCGR3B", "CXCR2", "MMP9") # 例如，使用这三个经典的中性粒细胞marker
  
  # 您可能已经有一个 Neutrophil_full 的定义，如果没有，也需要添加
  # 如果您还没有，可以先让它和 core 保持一致
  Neutrophil_full <- unique(c(
    Neutrophil_core,
    "CEACAM8", "CSF3R", "S100A8", "S100A9", "S100A12", "ELANE", "MPO"
  ))
  # --- 12.1 TAM subtypes (NEW): FOLR2+ vs SPP1+ ---
  # Enforce FOLR2 at first position even after extras merge
  FOLR2_TAM_core <- c("FOLR2","SELENOP","SLC40A1","F13A1","PDK4","STAB1","MS4A6A","OLFML2B","IGFBP4","C1R")
  FOLR2_TAM_full <- c("FOLR2","SELENOP","SLC40A1","F13A1","PDK4","STAB1","MS4A6A","OLFML2B","IGFBP4","IGFBP7",
                      "C1R","C1S","MRC1","CD163","LYVE1","MS4A4A","MERTK","GAS6","ENPP2","SLCO2B1","TREM2","GPX3")
  SPP1_TAM_core  <- c("SPP1","CXCL8","MMP9","MMP12","TREM1","BCL2A1","SERPINE1","CLEC5A","FBP1","SLC2A1")
  SPP1_TAM_full  <- c("SPP1","CXCL8","MMP9","MMP12","TREM1","BCL2A1","SERPINE1","CLEC5A","IL4I1","CA2",
                      "DDIT4","NUPR1","VCAN","LGALS3","GPNMB","TGM2","PLA2G7","SCD","MARCO","SLC11A1","CTSB","FTH1","VEGFA")
  # --- 12.2 TAM polarization (NEW): M1 vs M2 ---
  # M1: IFNγ/TLR4–STAT1/NF-κB inflammatory axis (anti-tumor-leaning)
  M1_TAM_core <- c(
    "IL1B","TNF","IL6","IL12A","IL12B",
    "CXCL9","CXCL10","CXCL11","CCL5",
    "STAT1","IRF1","IRF5",
    "CD80","CD86","CD40",
    "HLA-DRA","TAP1","B2M",
    "FCGR1A","ICAM1"
  )
  # full: add optional/IFN-induced/AP components (keep compact for ST sparsity)
  M1_TAM_full <- unique(c(
    M1_TAM_core,
    "NOS2",          # optional in human (often low), strong in mouse
    "GBP1","CIITA","HLA-DPA1","HLA-DPB1","HLA-DRB1","PSMB9"
  ))
  
  # M2: IL-4/IL-13–STAT6/PPARγ immunoregulatory/repair axis (pro-tumor-leaning)
  M2_TAM_core <- c(
    "MRC1","CD163","MSR1","MERTK",
    "IL10","TGFB1",
    "CCL17","CCL22","CCL18",
    "MMP9","MMP14",
    "VEGFA","PDGFB",
    "PPARG","KLF4","TREM2","APOE","FABP5",
    "CTSB","CTSD","STAB1"
  )
  # full: add context/ECM/scavenger options (kept moderate)
  M2_TAM_full <- unique(c(
    M2_TAM_core,
    "ARG1",          # optional in human (often low), strong in mouse
    "MARCO","SPARC","FN1","VCAN","ITGAV","ITGB5","CCR2"
  ))
  
  # --- 14. Neural/Schwann ---
  Neural_core <- c("SOX10","MPZ","PLP1")
  Neural_full <- c("MBP","PMP22","NGFR","EGR2","S100B","GFAP","NRXN1","CADM3","NCAM1","L1CAM",
                   "ERBB3","GJC3","MAL","MAG","CNP")
  
  # EMT anchor (tumor-intrinsic EMT/invasion axis)
  emt_anchor <- c("VIM","ZEB1","SNAI2","TWIST1","FN1","ITGA5","ITGB1","TGFBI","PDPN","SERPINE1","ITGB5")
  
  # ---- Merge user extras (also supports NEW sets) ----
  if (isTRUE(add_registered_extras)) {
    extra <- getOption("hd.scCBSI.extra",
                       default = list(
                         Tumor=character(0), Ductal=character(0), Acinar=character(0), Endocrine=character(0),
                         Endothelial=character(0), CAF=character(0), Perivascular=character(0), Immune=character(0),
                         T_Cell=character(0), B_Cell=character(0), Plasma=character(0), Macrophage=character(0),
                         Neutrophil=character(0), Neural=character(0),
                         Classical=character(0), Basal=character(0), IC=character(0), EMT=character(0),
                         # NEW extras:
                         # NEW extras:
                         iCAF=character(0), myCAF=character(0), apCAF=character(0),
                         FOLR2_TAM=character(0), SPP1_TAM=character(0),
                         M1_TAM=character(0),  M2_TAM=character(0)  # <--- add these
                         
                       )
    )
    # Merge extras
    Tumor_core <- unique(c(Tumor_core, extra$Tumor));   Tumor_full <- unique(c(Tumor_full, extra$Tumor))
    Ductal_core <- unique(c(Ductal_core, extra$Ductal)); Ductal_full <- unique(c(Ductal_full, extra$Ductal))
    Acinar_core <- unique(c(Acinar_core, extra$Acinar)); Acinar_full <- unique(c(Acinar_full, extra$Acinar))
    Endocrine_core <- unique(c(Endocrine_core, extra$Endocrine)); Endocrine_full <- unique(c(Endocrine_full, extra$Endocrine))
    Endothelial_core <- unique(c(Endothelial_core, extra$Endothelial)); Endothelial_full <- unique(c(Endothelial_full, extra$Endothelial))
    CAF_core <- unique(c(CAF_core, extra$CAF)); CAF_full <- unique(c(CAF_full, extra$CAF))
    iCAF_core <- unique(c(iCAF_core, extra$iCAF)); iCAF_full <- unique(c(iCAF_full, extra$iCAF))
    myCAF_core <- unique(c(myCAF_core, extra$myCAF)); myCAF_full <- unique(c(myCAF_full, extra$myCAF))
    apCAF_core <- unique(c(apCAF_core, extra$apCAF)); apCAF_full <- unique(c(apCAF_full, extra$apCAF))
    Perivascular_core <- unique(c(Perivascular_core, extra$Perivascular)); Perivascular_full <- unique(c(Perivascular_full, extra$Perivascular))
    Immune_core <- unique(c(Immune_core, extra$Immune)); Immune_full <- unique(c(Immune_full, extra$Immune))
    T_Cell_core <- unique(c(T_Cell_core, extra$T_Cell)); T_Cell_full <- unique(c(T_Cell_full, extra$T_Cell))
    B_Cell_core <- unique(c(B_Cell_core, extra$B_Cell)); B_Cell_full <- unique(c(B_Cell_full, extra$B_Cell))
    Plasma_core <- unique(c(Plasma_core, extra$Plasma)); Plasma_full <- unique(c(Plasma_full, extra$Plasma))
    Macrophage_core <- unique(c(Macrophage_core, extra$Macrophage)); Macrophage_full <- unique(c(Macrophage_full, extra$Macrophage))
    Neutrophil_core <- unique(c(Neutrophil_core, extra$Neutrophil)); Neutrophil_full <- unique(c(Neutrophil_full, extra$Neutrophil))
    Neural_core <- unique(c(Neural_core, extra$Neural)); Neural_full <- unique(c(Neural_full, extra$Neural))
    scClassical_core <- unique(c(scClassical_core, extra$Classical)); scClassical_full <- unique(c(scClassical_full, extra$Classical))
    scBasal_core <- unique(c(scBasal_core, extra$Basal)); scBasal_full <- unique(c(scBasal_full, extra$Basal))
    IC_core <- unique(c(IC_core, extra$IC)); IC_full <- unique(c(IC_full, extra$IC))
    emt_anchor <- unique(c(emt_anchor, extra$EMT))
    FOLR2_TAM_core <- unique(c(FOLR2_TAM_core, extra$FOLR2_TAM))
    FOLR2_TAM_full <- unique(c(FOLR2_TAM_full, extra$FOLR2_TAM))
    SPP1_TAM_core  <- unique(c(SPP1_TAM_core,  extra$SPP1_TAM))
    SPP1_TAM_full  <- unique(c(SPP1_TAM_full,  extra$SPP1_TAM))
    
    # ---- NEW: M1/M2 TAM extras
    M1_TAM_core <- unique(c(M1_TAM_core, extra$M1_TAM))
    M1_TAM_full <- unique(c(M1_TAM_full, extra$M1_TAM))
    M2_TAM_core <- unique(c(M2_TAM_core, extra$M2_TAM))
    M2_TAM_full <- unique(c(M2_TAM_full, extra$M2_TAM))
    
  }
  
  # ---- Assemble return sets ----
  # ---- Assemble return sets ----
  if (panel == "core") {
    sets <- list(
      Tumor = Tumor_core,
      scClassical = scClassical_core,
      scBasal = scBasal_core,
      IC = IC_core,
      Ductal = Ductal_core,
      Acinar = Acinar_core,
      Endocrine = Endocrine_core,
      Endothelial = Endothelial_core,
      iCAF = iCAF_core,
      myCAF = myCAF_core,
      apCAF = apCAF_core,
      CAF = CAF_core,
      Perivascular = Perivascular_core,
      Immune = Immune_core,
      T_Cell = T_Cell_core,           # <<< 新增
      B_Cell = B_Cell_core,           # <<< 新增
      Plasma = Plasma_core,           # <<< 新增
      # --- TAM subtypes (keep together)
      FOLR2_TAM = FOLR2_TAM_core,
      SPP1_TAM  = SPP1_TAM_core,
      M1_TAM    = M1_TAM_core,
      M2_TAM    = M2_TAM_core,
      Macrophage = Macrophage_core,
      Neutrophil = Neutrophil_core,
      Neural = Neural_core
    )
  } else {
    sets <- list(
      Tumor = Tumor_full,
      scClassical = scClassical_full,
      scBasal = scBasal_full,
      IC = IC_full,
      Ductal = Ductal_full,
      Acinar = Acinar_full,
      Endocrine = Endocrine_full,
      Endothelial = Endothelial_full,
      iCAF = iCAF_full,
      myCAF = myCAF_full,
      apCAF = apCAF_full,
      CAF = CAF_full,
      Perivascular = Perivascular_full,
      Immune = Immune_full,
      T_Cell = T_Cell_full,           # <<< 新增
      B_Cell = B_Cell_full,           # <<< 新增
      Plasma = Plasma_full,           # <<< 新增
      # --- TAM subtypes (keep together)
      FOLR2_TAM = FOLR2_TAM_full,
      SPP1_TAM  = SPP1_TAM_full,
      M1_TAM    = M1_TAM_full,
      M2_TAM    = M2_TAM_full,
      Macrophage = Macrophage_full,
      Neutrophil = Neutrophil_full,
      Neural = Neural_full
    )
  }
  
  
  # Optionally add EMT
  if (isTRUE(include_emt)) {
    sets[[emt_name]] <- emt_anchor
  }
  
  # Cross-set de-duplication with priority:
  # Give subtype panels higher priority than their parent class (CAF/TAM).
  if (isTRUE(dedup_across_sets)) {
    prio <- intersect(
      c("Tumor","scClassical","scBasal","IC","Ductal","Acinar","Endocrine","Endothelial",
        "iCAF","myCAF","apCAF","CAF",
        "Perivascular","Immune",
        # --- TAM subtype before general:
        "FOLR2_TAM","SPP1_TAM","M1_TAM","M2_TAM","Macrophage",
        "Neutrophil","Neural", emt_name),
      names(sets)
    )
    
    sets <- .hd_dedup_sets(sets, priority = prio)
  }
  sets
}







# ---- 仅保留对象中存在的基因，并打印缺失 ----
.hd_filter_present_markers2 <- function(obj, marker_sets, case_insensitive = TRUE, strip_suffix = TRUE) {
  # 取当前 assay 的所有特征名
  aa <- Seurat::DefaultAssay(obj)
  feats <- tryCatch(rownames(obj[[aa]]), error = function(e) rownames(obj))
  if (is.null(feats)) feats <- rownames(obj)
  
  norm <- function(x) {
    y <- x
    if (strip_suffix) y <- sub("\\.[0-9]+$", "", y)    # 去掉 .1/.2 之类 VEP 风格后缀
    if (strip_suffix) y <- sub("-[0-9]+$", "", y)      # 以及 -1/-2
    if (case_insensitive) y <- toupper(y)
    y
  }
  
  feats_key <- norm(feats)
  names(feats_key) <- feats  # 映回原名
  
  out <- lapply(marker_sets, function(gs) {
    kk <- norm(gs)
    hit <- kk %in% feats_key
    if (!any(hit)) return(character(0))
    # 把命中的 canonical（kk）映射回真实特征名
    remap <- vapply(kk[hit], function(k) {
      # 取第一个匹配的真实名
      names(feats_key)[match(k, feats_key)]
    }, FUN.VALUE = character(1))
    unname(unique(remap))
  })
  out
}






# ---- 安全子集：默认只保留含 "PDAC" 的群（和你之前一致）----
.hd_subset_by_groups <- function(object, group_by, include_groups = NULL) {
  stopifnot(group_by %in% colnames(object@meta.data))
  grp <- as.character(object@meta.data[[group_by]])
  keep <- if (is.null(include_groups)) {
    f <- grepl("PDAC", grp, ignore.case = TRUE); if (!any(f)) rep(TRUE, length(grp)) else f
  } else grp %in% include_groups
  if (!any(keep)) stop("[hd] No cells kept by include_groups/pattern in: ", group_by)
  cells <- rownames(object@meta.data)[keep]
  object[, cells, drop = FALSE]
}





# ---- 主函数：PDAC 所有类别的 DotPlot（复用你的 plot_marker_dot）----
hd_plot_PDAC_dot <- function(
    object,
    group_by       = "Seurat_L2_label",
    include_groups = NULL,
    assay          = NULL,
    which          = "all",
    # 透传到 plot_marker_dot():
    cols           = c("lightgrey","red"),
    dot_scale      = 6,
    scale_min      = 0,
    scale_max      = NA,
    x_axis_angle   = 45,
    legend_position= "right",
    legend_box     = "vertical",
    case_insensitive = TRUE,
    drop_empty_sets  = TRUE,
    grouping_map     = NULL,
    group_order      = NULL,
    panel            = c("core", "full"),
    include_emt      = TRUE,
    emt_name         = "EMT_anchor",
    group_style      = c("none","color_legend","label"),
    group_colors     = NULL,
    # --- 新增 ---
    dedup_across_sets = FALSE,   # <<< 关键：默认不跨集合去重，保住 CAF 的 COL1A1 等
    which_alias        = c(
      "Perivascular/Smooth" = "Perivascular",
      "Perivascular/Smooth muscle" = "Perivascular"
    )
) {
  group_style <- match.arg(group_style)
  panel <- match.arg(panel)
  if (!is.null(assay)) Seurat::DefaultAssay(object) <- assay
  
  obj_sub <- .hd_subset_by_groups(object, group_by = group_by, include_groups = include_groups)
  
  # 取全集（这里把去重开关透传进去）
  sets_all <- hd_get_marker_sets_PDAC(
    panel = panel,
    dedup_across_sets   = dedup_across_sets,  # <<< 透传
    add_registered_extras = TRUE,
    include_emt         = include_emt,
    emt_name            = emt_name
  )
  
  # 别名映射：把 which 中的用户口径转为集合名
  if (!identical(which, "all") && length(which_alias)) {
    which <- unname(ifelse(which %in% names(which_alias), which_alias[which], which))
  }
  
  targets <- if (identical(which, "all")) names(sets_all) else intersect(as.character(which), names(sets_all))
  if (length(targets) == 0) stop("`which` none matched; available: ", paste(names(sets_all), collapse = ", "))
  
  # 更稳健的基因存在性过滤（大小写不敏感、去掉“基因名.数字”后缀）
  marker_sets <- .hd_filter_present_markers2(obj_sub, sets_all[targets],
                                             case_insensitive = TRUE, strip_suffix = TRUE)
  marker_sets <- marker_sets[vapply(marker_sets, length, integer(1)) > 0]
  if (length(marker_sets) == 0) stop("No markers found in object for the requested sets.")
  
  stopifnot(is.function(get0("plot_marker_dot", mode = "function")))
  p <- plot_marker_dot(
    obj              = obj_sub,
    marker_sets      = marker_sets,
    group_by         = group_by,
    grouping_map     = grouping_map,
    group_order      = group_order,
    group_style      = group_style,
    group_colors     = group_colors,
    legend_position  = legend_position,
    legend_box       = legend_box,
    assay            = Seurat::DefaultAssay(obj_sub),
    cols             = cols,
    dot_scale        = dot_scale,
    scale_min        = scale_min,
    scale_max        = scale_max,
    strip_fill       = "#e6bac5",
    strip_size       = 9,
    gene_text_size   = 8,
    group_text_size  = 9,
    x_axis_angle     = x_axis_angle,
    case_insensitive = case_insensitive,
    drop_empty_sets  = drop_empty_sets
  )
  invisible(list(plot = p, markers_used = marker_sets, object_subset = obj_sub))
}

# ==== BANKSY/????/?? ====








# ---------- helpers ----------
.fmt_num <- function(x) sub("\\.?0+$", "", format(x, trim = TRUE))
.pick_col_safe <- function(df, col) {
  v <- df[, col, drop = TRUE]
  if (is.data.frame(v)) v <- v[,1, drop = TRUE]
  v
}




# ---------- helpers ----------
.fmt_num <- function(x) sub("\\.?0+$", "", format(x, trim = TRUE))
.pick_col_safe <- function(df, col) {
  v <- df[, col, drop = TRUE]
  if (is.data.frame(v)) v <- v[,1, drop = TRUE]
  v
}





suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)            # <- for CreateDimReducObject (v5 lives here)
  library(Matrix)
  library(SpatialExperiment)
  library(Banksy)
  library(SummarizedExperiment)
  library(dplyr)
})

# ---------- helpers ----------
.fmt_num <- function(x) sub("\\.?0+$", "", format(x, trim = TRUE))

.pick_col_safe <- function(df, col) {
  v <- df[, col, drop = TRUE]
  if (is.data.frame(v)) v <- v[, 1, drop = TRUE]
  v
}





.as_numeric_matrix_2col <- function(m) {
  m <- as.matrix(m)
  if (!is.numeric(m[,1])) m[,1] <- suppressWarnings(as.numeric(m[,1]))
  if (!is.numeric(m[,2])) m[,2] <- suppressWarnings(as.numeric(m[,2]))
  if (any(!is.finite(m))) stop("[BANKSY] Non-finite coords detected; check your coordinate columns.")
  colnames(m) <- c("sdimx","sdimy")
  m
}





if (!requireNamespace("SeuratObject", quietly = TRUE)) 
  stop("SeuratObject is required for writing UMAP back to Seurat.")
.create_reduction <- function(emb, key, assay){
  SeuratObject::CreateDimReducObject(embeddings = emb, key = key, assay = assay)
}






# ===========================================================
#                         MAIN
# ===========================================================
RunBANKSY_Cellbin <- function(
    sobj,
    assay             = "Spatial",
    batch_col         = "orig.ident",
    # coordinates: prefer scaled; fallback to raw
    coord_cols        = c("imagecol_scaled","imagerow_scaled"),
    coord_cols_fallback = c("imagecol","imagerow"),
    
    # geometry & embedding
    k_geom            = c(15L, 30L),
    npcs              = 30L,
    k_neighbors       = 50L,
    agf_mode          = c("auto","on","off"),
    seed              = 1L,
    
    # two lambdas & resolutions
    lambda_celltype   = 0.20,
    lambda_domain     = 0.80,
    res_celltype      = 1.2,
    res_domain        = 1.0,
    
    # optional: HVG subset (selected on normalized data; subset rows on counts)
    hvg_n             = NULL,
    
    # Banksy memory/compat
    row_limit_factor  = 0.5,
    
    # outputs
    add_umap          = TRUE,
    umap_min_dist     = 0.30,
    store_spe         = FALSE
){
  stopifnot(inherits(sobj, "Seurat"))
  agf_mode <- match.arg(agf_mode)
  
  # ----- counts (Seurat v5-friendly: prefer layer, fallback slot) -----
  DefaultAssay(sobj) <- assay
  counts <- tryCatch(
    GetAssayData(sobj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(sobj, assay = assay, slot = "counts")
  )
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  
  # ----- HVG subset with safeguards -----
  hvgs <- NULL
  if (!is.null(hvg_n)) {
    sobj_tmp <- sobj
    DefaultAssay(sobj_tmp) <- assay
    sobj_tmp <- NormalizeData(sobj_tmp, verbose = FALSE)
    sobj_tmp <- FindVariableFeatures(
      sobj_tmp, selection.method = "vst", nfeatures = hvg_n, verbose = FALSE
    )
    hvgs <- intersect(VariableFeatures(sobj_tmp), rownames(counts))
    if (length(hvgs) >= 50) {
      counts <- counts[hvgs, , drop = FALSE]
      message("[BANKSY] HVG subset: ", length(hvgs), " genes")
    } else {
      message("[BANKSY] HVG<50; skip HVG subset (using all genes).")
      hvgs <- NULL
    }
  }
  
  # ----- coordinates -----
  md <- sobj@meta.data
  cc <- if (all(coord_cols %in% colnames(md))) coord_cols else coord_cols_fallback
  if (!all(cc %in% colnames(md))) {
    stop("[BANKSY] Missing coordinate columns: ", paste(c(coord_cols, coord_cols_fallback), collapse=", "))
  }
  coords <- cbind(x = .pick_col_safe(md, cc[1]), y = .pick_col_safe(md, cc[2]))
  coords <- .as_numeric_matrix_2col(coords)
  
  # ----- batch column -----
  if (!batch_col %in% colnames(md)) {
    warning("[BANKSY] batch_col='", batch_col, "' not in meta.data; fallback to 'orig.ident'")
    batch_col <- "orig.ident"
  }
  
  # ----- assemble SpatialExperiment -----
  spe <- SpatialExperiment(
    assays        = list(counts = counts),
    spatialCoords = coords,
    colData       = S4Vectors::DataFrame(md[, batch_col, drop = FALSE])
  )
  
  # ----- detect function signatures / capabilities -----
  fm_cb     <- tryCatch(formals(Banksy::computeBanksy), error = function(e) NULL)
  fm_pca    <- tryCatch(formals(Banksy::runBanksyPCA),   error = function(e) NULL)
  fm_clu    <- tryCatch(formals(Banksy::clusterBanksy),  error = function(e) NULL)
  fm_umap   <- tryCatch(formals(Banksy::runBanksyUMAP),  error = function(e) NULL)
  
  has_compute_agf  <- !is.null(fm_cb)  && "compute_agf" %in% names(fm_cb)
  has_M_cb         <- !is.null(fm_cb)  && "M"           %in% names(fm_cb)
  has_chunk_cb     <- !is.null(fm_cb)  && "chunk_size"  %in% names(fm_cb)
  has_rowlim_cb    <- !is.null(fm_cb)  && "row_limit_factor" %in% names(fm_cb)
  
  has_use_agf_pca  <- !is.null(fm_pca) && "use_agf"     %in% names(fm_pca)
  has_use_agf_clu  <- !is.null(fm_clu) && "use_agf"     %in% names(fm_clu)
  has_use_agf_umap <- !is.null(fm_umap)&& "use_agf"     %in% names(fm_umap)
  
  # ----- AGF request & k_geom sanity -----
  k_geom <- as.integer(k_geom)
  if (length(k_geom) == 1L && agf_mode != "off") {
    k_geom <- c(k_geom[1], max(30L, as.integer(2L * k_geom[1])))
    message("[BANKSY] AGF requested but single k_geom provided; expanded to c(",
            paste(k_geom, collapse=", "), ").")
  }
  use_agf_requested <- switch(agf_mode,
                              auto = length(k_geom) > 1L,
                              on   = TRUE,
                              off  = FALSE)
  # Effective AGF requires PCA support
  use_agf_effective <- (use_agf_requested && has_use_agf_pca)
  
  # ----- capability log -----
  message(sprintf(
    "[BANKSY] Capabilities: compute_agf=%s | use_agf(PCA)=%s | use_agf(UMAP)=%s | use_agf(CLU)=%s | requested_AGF=%s | effective_AGF=%s",
    has_compute_agf, has_use_agf_pca, has_use_agf_umap, has_use_agf_clu,
    use_agf_requested, use_agf_effective
  ))
  
  # ----- computeBanksy (chunk sizing + guarded args) -----
  n_cells <- ncol(spe)
  safe_chunk <- max(
    1000L,
    min(nrow(counts),
        floor((.Machine$integer.max * row_limit_factor) / max(1L, n_cells)))
  )
  
  cb_args <- list(se = spe, assay_name = "counts", k_geom = k_geom)
  if (has_chunk_cb)     cb_args$chunk_size       <- as.integer(safe_chunk)
  if (has_rowlim_cb)    cb_args$row_limit_factor <- row_limit_factor
  if (has_compute_agf)  cb_args$compute_agf      <- use_agf_requested
  if (has_M_cb && use_agf_requested) cb_args$M   <- 1L  # include gradient (m=1)
  
  set.seed(seed)
  spe <- do.call(Banksy::computeBanksy, cb_args)
  
  # ===== PCA: Scheme-B (always M0; add M1 only if effective) =====
  pca_once <- function(se, use_agf_val) {
    args <- list(
      se     = se,
      lambda = c(lambda_celltype, lambda_domain),
      npcs   = as.integer(npcs),
      group  = batch_col
    )
    if (has_use_agf_pca) args$use_agf <- use_agf_val
    set.seed(seed)
    do.call(Banksy::runBanksyPCA, args)
  }
  # Always compute non-AGF PCA
  spe <- pca_once(spe, FALSE)
  # Compute AGF PCA only when effective
  if (use_agf_effective) {
    spe <- pca_once(spe, TRUE)
  }
  
  # ----- clustering (both lambdas) -----
  k_neighbors <- max(5L, min(as.integer(k_neighbors), ncol(spe) - 1L))
  clu_once <- function(se, lam, res) {
    args <- list(
      se = se, lambda = lam, algo = "leiden",
      k_neighbors = k_neighbors, resolution = res, seed = seed
    )
    if (has_use_agf_clu && use_agf_effective) args$use_agf <- TRUE
    do.call(Banksy::clusterBanksy, args)
  }
  spe <- clu_once(spe, lambda_celltype, res_celltype)
  spe <- clu_once(spe, lambda_domain,   res_domain)
  
  # ----- UMAP (optional; both lambdas) -----
  if (isTRUE(add_umap)) {
    run_umap_once <- function(se, lam) {
      args <- list(
        se = se, lambda = lam, npcs = as.integer(npcs),
        n_neighbors = k_neighbors, min_dist = umap_min_dist, seed = seed
      )
      if (has_use_agf_umap && use_agf_effective) args$use_agf <- TRUE
      do.call(Banksy::runBanksyUMAP, args)
    }
    spe <- run_umap_once(spe, lambda_celltype)
    spe <- run_umap_once(spe, lambda_domain)
  }
  
  # ----- pick cluster columns (robust matching on lambda/res; prefer M1 if effective) -----
  lam_str_ct <- .fmt_num(lambda_celltype)
  lam_str_dm <- .fmt_num(lambda_domain)
  cn <- Banksy::clusterNames(spe)
  
  pick_best_cluster <- function(cnames, lam_str, res_num, prefer_M = c("M1","M0")) {
    prefer_M <- match.arg(prefer_M)
    res_str  <- .fmt_num(res_num)
    pat1 <- paste0("clust_", prefer_M, ".*lam", gsub("\\.", "\\\\.", lam_str), ".*res", res_str)
    hit  <- grep(pat1, cnames, value = TRUE)
    if (length(hit)) return(hit[1])
    pat2 <- paste0("clust_", prefer_M, ".*lam", gsub("\\.", "\\\\.", lam_str))
    hit  <- grep(pat2, cnames, value = TRUE)
    if (length(hit)) return(hit[1])
    pat3 <- paste0("clust_.*lam", gsub("\\.", "\\\\.", lam_str), ".*res", res_str)
    hit  <- grep(pat3, cnames, value = TRUE)
    if (length(hit)) return(hit[1])
    pat4 <- paste0("clust_.*lam", gsub("\\.", "\\\\.", lam_str))
    hit  <- grep(pat4, cnames, value = TRUE)
    if (length(hit)) return(hit[1])
    stop("[BANKSY] No cluster column found for lambda=", lam_str, ". Available: ",
         paste(cnames, collapse=", "))
  }
  
  prefer_M <- if (use_agf_effective) "M1" else "M0"
  ct_col <- pick_best_cluster(cn, lam_str_ct, res_celltype, prefer_M = prefer_M)
  dm_col <- pick_best_cluster(cn, lam_str_dm, res_domain,   prefer_M = prefer_M)
  
  # ----- write back to Seurat meta -----
  if (!setequal(colnames(sobj), colnames(spe))) {
    warning("[BANKSY] Seurat and SPE cell name sets differ; writing back on intersection only.")
  }
  common <- intersect(colnames(sobj), colnames(spe))
  ord_s  <- match(common, colnames(sobj))
  ord_e  <- match(common, colnames(spe))
  
  ct_vec <- as.character(SummarizedExperiment::colData(spe)[[ct_col]])[ord_e]
  dm_vec <- as.character(SummarizedExperiment::colData(spe)[[dm_col]])[ord_e]
  
  ct_name <- paste0("banksy_celltype_lam", lam_str_ct, "_res", .fmt_num(res_celltype))
  dm_name <- paste0("banksy_domain_lam",   lam_str_dm, "_res", .fmt_num(res_domain))
  sobj[[ct_name]] <- NA
  sobj[[dm_name]] <- NA
  sobj@meta.data[ord_s, ct_name] <- ct_vec
  sobj@meta.data[ord_s, dm_name] <- dm_vec
  
  sobj$banksy_celltype <- sobj[[ct_name]]
  sobj$banksy_domain   <- sobj[[dm_name]]
  
  # ----- add UMAP reductions back to Seurat (prefer effective branch; fallback if missing) -----
  if (isTRUE(add_umap)) {
    rdn_names <- reducedDimNames(spe)
    umap_name_for <- function(lam_str, prefer_M = c("M1","M0")) {
      prefer_M <- match.arg(prefer_M)
      prefer <- grep(paste0("^UMAP_", prefer_M, ".*lam", gsub("\\.", "\\\\.", lam_str), "$"),
                     rdn_names, value = TRUE)
      if (length(prefer)) return(prefer[1])
      alt <- grep(paste0("^UMAP_", if (prefer_M=="M1") "M0" else "M1",
                         ".*lam", gsub("\\.", "\\\\.", lam_str), "$"),
                  rdn_names, value = TRUE)
      if (length(alt)) return(alt[1])
      any_hit <- grep(paste0("UMAP.*lam", gsub("\\.", "\\\\.", lam_str)), rdn_names, value = TRUE)
      if (length(any_hit)) return(any_hit[1])
      return(NA_character_)
    }
    umap_name_ct <- umap_name_for(lam_str_ct, prefer_M = prefer_M)
    umap_name_dm <- umap_name_for(lam_str_dm, prefer_M = prefer_M)
    
    add_umap_to_seurat <- function(umap_name, red_name) {
      if (!is.na(umap_name)) {
        emb <- as.matrix(SingleCellExperiment::reducedDim(spe, umap_name))[ord_e, , drop = FALSE]
        key_here <- if (identical(red_name, "banksy_umap_domain")) "BYD_" else "BY_"
        dr  <- .create_reduction(emb, key = key_here, assay = assay)
        sobj@reductions[[red_name]] <- dr
      }
    }
    
    add_umap_to_seurat(umap_name_ct, "banksy_umap")
    add_umap_to_seurat(umap_name_dm, "banksy_umap_domain")
  }
  
  # ===== 写回校验 & 自动回退（从 SPE 补挂） =====
  wrote_ct <- "banksy_umap" %in% names(sobj@reductions) &&
    ncol(sobj@reductions$banksy_umap@cell.embeddings) >= 2
  wrote_dm <- "banksy_umap_domain" %in% names(sobj@reductions) &&
    ncol(sobj@reductions$banksy_umap_domain@cell.embeddings) >= 2
  
  if ((!wrote_ct || !wrote_dm)) {
    # 尝试从 SPE 兜底补挂
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      warning("[BANKSY] UMAP computed but not written; and package 'SingleCellExperiment' not available for fallback.")
    } else {
      tryCatch({
        # 再次确认我们有可用的 dimred 名
        dn <- SingleCellExperiment::reducedDimNames(spe)
        
        if (!wrote_ct && !is.na(umap_name_ct) && umap_name_ct %in% dn) {
          emb_ct <- as.matrix(SingleCellExperiment::reducedDim(spe, umap_name_ct))[ord_e, , drop = FALSE]
          sobj@reductions$banksy_umap <- SeuratObject::CreateDimReducObject(
            embeddings = emb_ct, key = "BY_", assay = assay
          )
          wrote_ct <- TRUE
        }
        if (!wrote_dm && !is.na(umap_name_dm) && umap_name_dm %in% dn) {
          emb_dm <- as.matrix(SingleCellExperiment::reducedDim(spe, umap_name_dm))[ord_e, , drop = FALSE]
          sobj@reductions$banksy_umap_domain <- SeuratObject::CreateDimReducObject(
            embeddings = emb_dm, key = "BYD_", assay = assay
          )
          wrote_dm <- TRUE
        }
      }, error = function(e) {
        warning("[BANKSY] Fallback UMAP write failed: ", conditionMessage(e))
      })
    }
  }
  
  # 最终仍未写回则报警
  if (!wrote_ct) warning("[BANKSY] UMAP (celltype) computed but not written into Seurat reductions.")
  if (!wrote_dm) warning("[BANKSY] UMAP (domain)   computed but not written into Seurat reductions.")
  # ===== 结束 =====
  
  
  # ----- record run info -----
  info <- list(
    package_version      = tryCatch(as.character(utils::packageVersion("Banksy")), error=function(e) NA),
    assay                = assay,
    coord_cols_used      = cc,
    k_geom_used          = k_geom,
    npcs                 = npcs,
    k_neighbors          = k_neighbors,
    seed                 = seed,
    hvg_n_requested      = hvg_n,
    hvg_effective        = if (!is.null(hvgs)) length(hvgs) else NA_integer_,
    lambda_celltype      = lambda_celltype,
    res_celltype         = res_celltype,
    lambda_domain        = lambda_domain,
    res_domain           = res_domain,
    batch_col            = batch_col,
    agf_requested        = use_agf_requested,
    agf_effective        = use_agf_effective,
    has_compute_agf      = has_compute_agf,
    has_use_agf_pca      = has_use_agf_pca,
    has_use_agf_umap     = has_use_agf_umap,
    has_use_agf_clu      = has_use_agf_clu,
    ct_colname_in_spe    = ct_col,
    dm_colname_in_spe    = dm_col,
    umap_ct_dimred       = if (isTRUE(add_umap)) umap_name_ct else NA_character_,
    umap_dm_dimred       = if (isTRUE(add_umap)) umap_name_dm else NA_character_,
    row_limit_factor     = row_limit_factor,
    umap_min_dist        = umap_min_dist
  )
  sobj@misc$banksy <- info
  if (isTRUE(store_spe)) sobj@misc$banksy_spe <- spe  # WARNING: can be large
  
  message("[BANKSY] Done. celltype(", length(unique(na.omit(sobj$banksy_celltype))),
          ") | domain(", length(unique(na.omit(sobj$banksy_domain))), ") written",
          if (isTRUE(add_umap)) " + UMAP reductions" else "")
  
  invisible(sobj)
}







# ============================================================
# Triple-Niche detection & standardized visualization suite
#   - Single-sample Visium HD / Cellbin
#   - Targets example: c("Tumor_EMT","myCAF","TAM_SPP1")
# Deps: Matrix, igraph, dbscan, FNN, ggplot2 (patchwork optional)
# ============================================================

# ---------- helpers: scalars from Seurat object ----------
.get_lowres_mpp <- function(object) {
  mpp_full <- object@misc$microns_per_pixel
  if (is.null(mpp_full)) {
    mpp_full <- tryCatch(Images(object)[[1]]@scale.factors$microns_per_pixel, error=function(e) NA_real_)
  }
  s_low <- object@misc$spatial_scales$tissue_lowres_scalef
  if (is.null(s_low)) s_low <- object@misc$scales$tissue_lowres_scalef
  if (is.null(s_low)) s_low <- tryCatch(Images(object)[[1]]@scale.factors$tissue_lowres_scalef, error=function(e) NA_real_)
  mpp_low <- if (is.finite(mpp_full) && is.finite(s_low) && s_low > 0) mpp_full / s_low else NA_real_
  list(mpp_full = mpp_full, s_low = s_low, mpp_low = mpp_low)
}





# ---------- 1) Detect triple-cell co-occurrence niche ----------
DetectNiche_Triple <- function(
    object,
    type_col,
    targets,
    mode        = c("radius","knn"),
    radius_um   = 50,
    k_neighbors = 20,
    min_each    = 1L,
    n_perm      = 200L,
    write_meta  = TRUE,
    meta_prefix = "TriNiche",
    verbose     = TRUE
){
  stopifnot(inherits(object, "Seurat"))
  mode <- match.arg(mode)
  need_pkg <- c("Matrix","igraph")
  for (p in need_pkg) if (!requireNamespace(p, quietly = TRUE)) stop("Please install '", p, "'.")
  has_dbscan <- requireNamespace("dbscan", quietly = TRUE)
  has_FNN    <- requireNamespace("FNN",    quietly = TRUE)
  
  md <- object@meta.data
  # coordinates
  if (!all(c("imagecol_scaled","imagerow_scaled") %in% names(md))) {
    if (!all(c("imagecol","imagerow") %in% names(md)))
      stop("Need 'imagecol_scaled/imagerow_scaled' or 'imagecol/imagerow' in meta.")
    x <- md$imagecol; y <- md$imagerow
  } else { x <- md$imagecol_scaled; y <- md$imagerow_scaled }
  coords <- cbind(as.numeric(x), as.numeric(y)); rownames(coords) <- rownames(md)
  
  if (!type_col %in% names(md)) stop("type_col '", type_col, "' not found in meta.data.")
  types <- as.character(md[[type_col]])
  lv_all <- sort(unique(types))
  lv_targets <- unique(as.character(targets))
  miss <- setdiff(lv_targets, lv_all); if (length(miss)) stop("targets not found: ", paste(miss, collapse=", "))
  
  # radius->pixels
  sc <- .get_lowres_mpp(object)
  if (mode == "radius" && !is.finite(sc$mpp_low)) {
    warning("[DetectNiche_Triple] Missing microns-to-lowres-pixel conversion; falling back to KNN.")
    mode <- "knn"
  }
  
  # neighborhoods
  neighbor_list <- NULL
  if (mode == "radius") {
    if (!has_dbscan) stop("Please install 'dbscan' (for radius neighborhoods).")
    r_px <- radius_um / sc$mpp_low
    if (verbose) message(sprintf("• Neighborhood: radius = %.1f µm (~%.1f px)", radius_um, r_px))
    fr <- dbscan::frNN(coords, eps = r_px, sort = FALSE)
    neighbor_list <- fr$id
  } else {
    if (!has_FNN) stop("Please install 'FNN' (for KNN).")
    k <- max(1, min(k_neighbors, nrow(coords)-1))
    if (verbose) message(sprintf("• Neighborhood: KNN (k=%d)", k))
    kn <- FNN::get.knn(coords, k = k)
    neighbor_list <- lapply(seq_len(nrow(coords)), function(i) c(i, kn$nn.index[i, ]))
  }
  
  # local composition & triple co-occurrence
  .count_vec <- function(idx, pool) {
    if (length(idx) == 0) return(setNames(rep(0L, length(pool)), pool))
    tab <- table(types[idx]); out <- setNames(rep(0L, length(pool)), pool)
    inter <- intersect(names(tab), pool); out[inter] <- as.integer(tab[inter]); out
  }
  
  n <- nrow(coords)
  counts_targets <- matrix(0L, nrow=n, ncol=length(lv_targets), dimnames=list(rownames(coords), lv_targets))
  comp_all       <- matrix(0,   nrow=n, ncol=length(lv_all),     dimnames=list(rownames(coords), lv_all))
  cooccur_flag   <- logical(n)
  
  for (i in seq_len(n)) {
    nb <- neighbor_list[[i]]
    ct <- .count_vec(nb, lv_targets); counts_targets[i, ] <- ct
    ca <- .count_vec(nb, lv_all); comp_all[i, ] <- if (sum(ca) > 0) ca / sum(ca) else ca
    cooccur_flag[i] <- all(ct >= min_each)
  }
  
  # connected components on co-occurrence subgraph
  comp_id <- integer(n); comp_id[] <- 0L
  if (any(cooccur_flag)) {
    idx_keep <- which(cooccur_flag)
    map_old2new <- integer(n); map_old2new[idx_keep] <- seq_along(idx_keep)
    edges <- cbind(integer(0), integer(0))
    for (ii in idx_keep) {
      nb <- intersect(neighbor_list[[ii]], idx_keep)
      if (length(nb)) {
        a <- rep(map_old2new[ii], length(nb)); b <- map_old2new[nb]
        edges <- rbind(edges, cbind(a, b))
      }
    }
    if (nrow(edges) > 0) {
      g <- igraph::graph_from_edgelist(unique(t(apply(edges,1,sort))), directed = FALSE)
      comp <- igraph::components(g)$membership
      comp_id[idx_keep] <- comp
    }
  }
  
  # permutation enrichment (global)
  p_value <- NA_real_
  if (n_perm > 0L) {
    set.seed(1); obs <- sum(cooccur_flag); perm_counts <- integer(n_perm)
    for (b in seq_len(n_perm)) {
      types_shuf <- sample(types, length(types), replace = FALSE)
      flag_b <- logical(n)
      for (i in seq_len(n)) {
        nb <- neighbor_list[[i]]; if (length(nb) == 0) next
        tt <- table(types_shuf[nb])
        ok <- all(vapply(lv_targets, function(tg) (as.integer(tt[tg]) %||% 0L) >= min_each, logical(1)))
        flag_b[i] <- ok
      }
      perm_counts[b] <- sum(flag_b)
    }
    p_value <- mean(perm_counts >= obs)
    if (verbose) message(sprintf("• Global enrichment: observed=%d; perm p=%.4g (n=%d)", obs, p_value, n_perm))
  }
  
  # write-back
  if (write_meta) {
    object[[paste0(meta_prefix, "_flag")]]   <- cooccur_flag
    object[[paste0(meta_prefix, "_compID")]] <- comp_id
    for (tg in lv_targets) {
      object[[paste0(meta_prefix, "_cnt_",  tg)]] <- counts_targets[, tg]
      # write local proportion of target among all types (from comp_all)
      if (tg %in% colnames(comp_all))
        object[[paste0(meta_prefix, "_prop_", tg)]] <- comp_all[, tg]
    }
  }
  
  list(
    cooccur_flag   = cooccur_flag,
    present_counts = as.data.frame(counts_targets),
    comp           = as.data.frame(comp_all),
    comp_targets   = as.data.frame(comp_all[, lv_targets, drop=FALSE]),
    components     = comp_id,
    p_value        = p_value,
    params         = list(mode=mode, radius_um=if (mode=="radius") radius_um else NA_real_,
                          k_neighbors=if (mode=="knn") k_neighbors else NA_integer_,
                          min_each=min_each, n_perm=n_perm, type_col=type_col, targets=lv_targets,
                          mpp_full=sc$mpp_full, s_low=sc$s_low, mpp_low=sc$mpp_low),
    object         = if (write_meta) object else NULL
  )
}





# ---------- 2) Density score (hotspot) ----------
ComputeTripleNicheDensity <- function(object, flag_col = "TriNiche_flag",
                                      radius_um = 80, meta_name = "TriNiche_density") {
  if (!flag_col %in% colnames(object@meta.data)) stop("Missing flag column: ", flag_col)
  md <- object@meta.data
  x <- if ("imagecol_scaled" %in% names(md)) md$imagecol_scaled else md$imagecol
  y <- if ("imagerow_scaled" %in% names(md)) md$imagerow_scaled else md$imagerow
  coords <- cbind(as.numeric(x), as.numeric(y))
  sc <- .get_lowres_mpp(object)
  if (!is.finite(sc$mpp_low)) stop("Cannot compute density: missing microns-per-pixel conversion.")
  r_px <- radius_um / sc$mpp_low
  if (!requireNamespace("dbscan", quietly = TRUE)) stop("Please install 'dbscan'.")
  fr <- dbscan::frNN(coords, eps = r_px, sort = FALSE)
  flag <- as.logical(md[[flag_col]])
  dens <- vapply(seq_along(fr$id), function(i) {
    nb <- fr$id[[i]]; if (length(nb) == 0) return(0)
    mean(flag[nb])
  }, numeric(1))
  object[[meta_name]] <- dens
  object
}





# ---------- 3) Standardized plots ----------
# 3a) Hotspot map (continuous density)
PlotTripleNicheHotspots <- function(object, density_col = "TriNiche_density",
                                    title = "Triple Niche Hotspots",
                                    ptsize = 1.6) {
  if (!exists("PlotExpressionV2", mode = "function"))
    stop("PlotExpressionV2() is required in your pipeline.")
  PlotExpressionV2(
    object = object,
    feature = density_col,
    source  = "meta",
    shape   = "square",
    spot_size = "auto",
    knn_smooth = FALSE,
    legend_title = title
  ) + ggplot2::ggtitle(title)
}





# 3b) Component outlines (convex hulls over niche-positive nodes)
PlotTripleNicheComponentsOutline <- function(object,
                                             flag_col = "TriNiche_flag",
                                             comp_col = "TriNiche_compID",
                                             min_points = 5L,
                                             line_size = 0.6,
                                             line_color = "#2C7FB8",
                                             fill_alpha = 0.08,
                                             ptsize = 0.6) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  md <- object@meta.data
  if (!all(c(flag_col, comp_col) %in% colnames(md))) stop("Missing columns for outline plotting.")
  df <- data.frame(
    x = if ("imagecol_scaled" %in% names(md)) md$imagecol_scaled else md$imagecol,
    y = if ("imagerow_scaled" %in% names(md)) -md$imagerow_scaled else -md$imagerow,
    flag = md[[flag_col]],
    comp = md[[comp_col]]
  )
  base <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(data = subset(df, !flag), color = "grey85", size = ptsize) +
    ggplot2::geom_point(data = subset(df, flag),   color = "#E41A1C", size = ptsize) +
    ggplot2::coord_fixed() + ggplot2::theme_void() +
    ggplot2::ggtitle("Triple Niche Components (Convex Hull Outlines)")
  
  comps <- unique(df$comp[df$flag & df$comp > 0])
  if (length(comps)) {
    for (cid in comps) {
      sub <- df[df$comp == cid & df$flag, c("x","y")]
      if (nrow(sub) >= min_points) {
        idx <- grDevices::chull(sub$x, sub$y)
        hull <- sub[idx, ]
        base <- base +
          ggplot2::geom_polygon(data = hull, ggplot2::aes(x = x, y = y),
                                inherit.aes = FALSE,
                                color = line_color, fill = line_color,
                                linewidth = line_size, alpha = fill_alpha)
      }
    }
  }
  base
}





# 3c) Composition bars per component (top-N by size)
PlotTripleNicheCompositionBars <- function(object,
                                           type_col = "FinalCoarse",
                                           comp_col = "TriNiche_compID",
                                           flag_col = "TriNiche_flag",
                                           top_n = 5L,
                                           targets_focus = NULL,
                                           normalize = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  md <- object@meta.data
  if (!all(c(type_col, comp_col, flag_col) %in% colnames(md))) stop("Missing columns in meta.")
  df <- md[md[[flag_col]] & md[[comp_col]] > 0, c(type_col, comp_col)]
  colnames(df) <- c("type","comp")
  if (nrow(df) == 0) stop("No niche-positive components to plot.")
  comp_sizes <- sort(table(df$comp), decreasing = TRUE)
  keep <- as.integer(names(head(comp_sizes, top_n)))
  df <- df[df$comp %in% keep, , drop = FALSE]
  tab <- as.data.frame.matrix(table(df$comp, df$type))
  tab$comp <- as.integer(rownames(tab))
  tall <- tidyr::pivot_longer(tab, cols = -comp, names_to = "type", values_to = "count")
  if (!is.null(targets_focus)) tall <- tall[tall$type %in% targets_focus, , drop = FALSE]
  if (normalize) {
    tall <- dplyr::group_by(tall, comp) |>
      dplyr::mutate(prop = count / pmax(sum(count), 1)) |>
      dplyr::ungroup()
  } else { tall$prop <- tall$count }
  ggplot2::ggplot(tall, ggplot2::aes(x = factor(comp), y = prop, fill = type)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::labs(x = "Component (by size)", y = if (normalize) "Proportion" else "Count",
                  title = "Triple Niche Composition (Top Components)") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}





# ---------- 4) One-click pipeline & save ----------
RunTripleNicheSuite <- function(
    object,
    type_col,
    targets,
    outdir,
    radius_um_detect = 50,
    min_each = 1L,
    n_perm = 200L,
    radius_um_density = 80,
    mode = "radius",
    k_neighbors = 20,
    meta_prefix = "TriNiche",
    ptsize = 1.6
){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  res <- DetectNiche_Triple(
    object      = object,
    type_col    = type_col,
    targets     = targets,
    mode        = mode,
    radius_um   = radius_um_detect,
    k_neighbors = k_neighbors,
    min_each    = min_each,
    n_perm      = n_perm,
    write_meta  = TRUE,
    meta_prefix = meta_prefix
  )
  obj2 <- res$object
  
  # density
  obj2 <- ComputeTripleNicheDensity(obj2, flag_col = paste0(meta_prefix,"_flag"),
                                    radius_um = radius_um_density,
                                    meta_name = paste0(meta_prefix,"_density"))
  
  # plots
  p_hot  <- PlotTripleNicheHotspots(obj2, density_col = paste0(meta_prefix,"_density"),
                                    title = "Triple Niche Hotspots", ptsize = ptsize)
  p_out  <- PlotTripleNicheComponentsOutline(obj2,
                                             flag_col = paste0(meta_prefix,"_flag"),
                                             comp_col = paste0(meta_prefix,"_compID"),
                                             ptsize = 0.5)
  p_bars <- PlotTripleNicheCompositionBars(obj2,
                                           type_col = type_col,
                                           comp_col = paste0(meta_prefix,"_compID"),
                                           flag_col = paste0(meta_prefix,"_flag"),
                                           top_n = 6L,
                                           targets_focus = targets,
                                           normalize = TRUE)
  
  # save
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  ggplot2::ggsave(file.path(outdir, "niche_hotspots.pdf"), p_hot, width = 6.5, height = 5.5)
  ggplot2::ggsave(file.path(outdir, "niche_components_outline.pdf"), p_out, width = 6.5, height = 5.5)
  ggplot2::ggsave(file.path(outdir, "niche_composition_bars.pdf"), p_bars, width = 7.5, height = 4.8)
  
  # also write CSVs
  utils::write.csv(res$present_counts, file.path(outdir, "triple_niche_target_counts_per_spot.csv"), row.names = TRUE)
  utils::write.csv(res$comp,           file.path(outdir, "local_composition_alltypes_per_spot.csv"), row.names = TRUE)
  
  invisible(list(object = obj2, results = res,
                 plots = list(hotspots = p_hot, outline = p_out, bars = p_bars),
                 outdir = outdir))
}





library(FNN)  # 用于快速最近邻查找

FindInvasiveFront <- function(
    object, 
    celltype_col = "celltype",        # 元数据中细胞类型列名
    tumor_types = "Tumor",            # 肿瘤细胞类型标识（向量可指定多个类别）
    stroma_types = c("CAF","Immune"), # 基质/免疫细胞类型标识集合
    radius = NULL,                    # 半径阈值（优先使用半径邻域）
    k = 5,                            # 近邻数（radius为空时使用KNN邻域）
    score_col = "tumor_stroma_score", # 新元数据列名：相互作用分数
    flag_col = "is_invasive_front",   # 新元数据列名：是否属于前沿热点
    hotspot_thresh = 0.5              # 判定前沿的分数阈值（此处举例0.5，可调）
) {
  # 提取坐标矩阵
  coords <- as.matrix(object@meta.data[, c("imagecol_scaled","imagerow_scaled")])
  n <- nrow(coords)
  # 确定邻居列表
  if (!is.null(radius)) {
    # 基于半径：先找上限KNN，再筛选距离radius内的
    k_max <- min(n - 1, max(50, k * 3))
    knn_res <- get.knn(coords, k = k_max)
    neighbor_index <- apply(knn_res$nn.index, 1, function(idx_vec) {
      # 去除自身0和NA，然后筛选距离半径内
      d_vec <- knn_res$nn.dist[ idx_vec[1], ]  # FNN::get.knn按行输出邻居距离
      hits <- which(d_vec <= radius)
      # 返回有效邻居索引（排除自身）
      idx <- idx_vec[hits]
      idx[idx > 0 & idx != idx_vec[1]]
    })
  } else {
    # 基于KNN：直接获取每个点的k近邻索引（包含自身，需要去除）
    knn_res <- get.knn(coords, k = k + 1)  # +1留出自身
    neighbor_index <- apply(knn_res$nn.index, 1, function(idx_vec) {
      idx <- idx_vec[idx_vec > 0]        # 去除无效索引（0可能代表自身或无邻居）
      idx <- idx[idx != idx_vec[1]]      # 排除自身
      head(idx, k)                       # 截取前k个邻居
    })
  }
  # 确保neighbor_index按行对应每个spot邻居向量
  # 计算肿瘤-基质相互作用分数
  cell_types <- object@meta.data[[celltype_col]]
  scores <- numeric(n)
  for (i in 1:n) {
    if (cell_types[i] %in% tumor_types) {
      nb_idx <- neighbor_index[[i]]
      if (length(nb_idx) == 0) {
        scores[i] <- 0  # 无邻居则分数记为0
      } else {
        # 统计邻居中非肿瘤细胞的比例
        nb_types <- cell_types[nb_idx]
        n_stroma <- sum(nb_types %in% stroma_types)
        scores[i] <- n_stroma / length(nb_idx)
      }
    } else {
      scores[i] <- NA  # 非肿瘤细胞不计算分数（可设NA）
    }
  }
  # 加入元数据
  object@meta.data[[score_col]] <- scores
  # 判定侵袭前沿热点（根据阈值）
  hotspot_flag <- ifelse(!is.na(scores) & scores >= hotspot_thresh, TRUE, FALSE)
  object@meta.data[[flag_col]] <- hotspot_flag
  return(object)
}








.get_lowres_mpp <- function(object) {
  # 1) mpp (full-res)
  mpp_full <- tryCatch({
    if (!is.null(object@misc$microns_per_pixel)) object@misc$microns_per_pixel else NA_real_
  }, error=function(e) NA_real_)
  if (!is.finite(mpp_full)) {
    mpp_full <- tryCatch(Images(object)[[1]]@scale.factors$microns_per_pixel, error=function(e) NA_real_)
  }
  if (!is.finite(mpp_full) && !is.null(object@misc$scales$microns_per_pixel)) {
    mpp_full <- object@misc$scales$microns_per_pixel
  }
  
  # 2) lowres scale factor (tissue_lowres_scalef)
  s_low <- tryCatch(object@misc$spatial_scales$tissue_lowres_scalef, error=function(e) NULL)
  if (is.null(s_low)) s_low <- tryCatch(object@misc$scales$tissue_lowres_scalef, error=function(e) NULL)
  if (is.null(s_low)) s_low <- tryCatch(Images(object)[[1]]@scale.factors$lowres, error=function(e) NA_real_)
  
  mpp_low <- if (is.finite(mpp_full) && is.finite(s_low) && s_low > 0) mpp_full / s_low else NA_real_
  list(mpp_full = mpp_full, s_low = s_low, mpp_low = mpp_low)
}






DetectInvasionHotspot_BanksyUCell <- function(
    object,
    domain_col,                        # e.g. "BANKSY_Domain" or "Seurat_L2_label"
    tumor_domains,                     # character vec of domain labels considered "Tumor"
    stroma_domains,                    # character vec of domain labels considered "Stroma" (CAF/immune/ECM/vasc)
    module_cols,                       # named list of module score columns (>=1): e.g.
    # list(EMT="ucell_EMT", myCAF="ucell_myCAF", M2="ucell_M2TAM",
    #      Angio="ucell_Angiogenesis", ECM="ucell_ECMremodel")
    weights       = NULL,              # optional numeric weights, same length as module_cols
    mode          = c("radius","knn"),
    radius_um_nb  = 150,               # neighborhood radius for smoothing/graph (µm)  [↑ 默认放大，避免邻域过稀]
    radius_um_den = 200,               # radius for density score (µm)
    q_hot         = 0.90,              # quantile threshold within candidate set
    interface_tau = 0.10,              # minimal interface score to enter candidate set
    min_comp_size = 30L,               # minimal component size (#spots)
    write_meta    = TRUE,
    meta_prefix   = "INVA5",
    verbose       = TRUE,
    rescale_modules = FALSE,  # ← 新增：默认不再二次标准化
    concave_outline = FALSE,   # ← 新增：是否用凹壳描边
    concavity       = 2.0,     # 可选：凹度（越大越贴合，注意噪点）
    length_thresh   = 0.0    # 可选：最短边过滤
    
    
){
  stopifnot(inherits(object, "Seurat"))
  mode <- match.arg(mode)
  
  # deps
  if (!requireNamespace("dbscan", quietly = TRUE)) stop("Please install 'dbscan'")
  if (!requireNamespace("FNN",    quietly = TRUE)) stop("Please install 'FNN'")
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Please install 'igraph'")
  if (!requireNamespace("ggplot2",quietly = TRUE)) stop("Please install 'ggplot2'")
  
  md <- object@meta.data
  
  # ---- coordinates (prefer scaled; fallback raw) ----
  x <- if ("imagecol_scaled" %in% names(md)) md$imagecol_scaled else md$imagecol
  y <- if ("imagerow_scaled" %in% names(md)) md$imagerow_scaled else md$imagerow
  if (any(is.na(x)) || any(is.na(y))) stop("Missing spatial coordinates in meta.")
  coords <- cbind(as.numeric(x), as.numeric(y)); rownames(coords) <- rownames(md)
  

  if (!domain_col %in% names(md)) stop("domain_col not found in meta: ", domain_col)
  dom <- trimws(as.character(md[[domain_col]]))
  tumor_set  <- trimws(as.character(tumor_domains))
  stroma_set <- trimws(as.character(stroma_domains))
  
  # warn if unknown labels provided
  all_labels <- unique(dom)
  unk_t <- setdiff(tumor_set,  all_labels); if (length(unk_t)) warning("Unknown tumor_domains: ", paste(unk_t, collapse=", "))
  unk_s <- setdiff(stroma_set, all_labels); if (length(unk_s)) warning("Unknown stroma_domains: ", paste(unk_s, collapse=", "))
  
  in_tumor  <- dom %in% tumor_set
  in_stroma <- dom %in% stroma_set
  

  mod_names <- names(module_cols)
  if (is.null(mod_names) || any(mod_names == "")) stop("module_cols must be a *named* list.")
  miss <- setdiff(unname(unlist(module_cols)), colnames(md))
  if (length(miss)) stop("Missing module score columns in meta: ", paste(miss, collapse=", "))
  
  M <- matrix(NA_real_, nrow=nrow(coords), ncol=length(module_cols),
              dimnames = list(rownames(coords), names(module_cols)))
  for (j in seq_along(module_cols)) M[, j] <- as.numeric(md[[ module_cols[[j]] ]])
  

  sc <- .get_lowres_mpp(object)
  if (verbose && is.finite(sc$mpp_low)) {
    message(sprintf("Scale: mpp_full=%.4f, lowres_scalef=%.6f -> mpp_low=%.4f µm/px",
                    sc$mpp_full, sc$s_low, sc$mpp_low))
  }
  if (!is.finite(sc$mpp_low)) {
    warning("[INVA5] Missing microns->lowres pixel scale. Falling back to KNN neighbors.")
    mode <- "knn"
  }
  r_nb_px  <- if (mode=="radius") radius_um_nb  / sc$mpp_low else NA_real_
  r_den_px <- if (mode=="radius") radius_um_den / sc$mpp_low else NA_real_
  
  # ---- 1) build neighborhood (smoothing & graph) ----
  if (mode == "radius") {
    fr_nb  <- dbscan::frNN(coords, eps = r_nb_px,  sort = FALSE)
    fr_den <- dbscan::frNN(coords, eps = r_den_px, sort = FALSE)
    
    # include self; empty -> self-only
    nb_list  <- lapply(seq_len(nrow(coords)), function(i) {
      nb <- fr_nb$id[[i]]
      if (length(nb) == 0) i else unique(c(i, nb))
    })
    den_list <- lapply(seq_len(nrow(coords)), function(i) {
      nb <- fr_den$id[[i]]
      if (length(nb) == 0) i else unique(c(i, nb))
    })
    
    deg <- vapply(nb_list, length, 1L)
    if (verbose) message(sprintf("• Neighborhood by radius: %.1f µm (~%.1f px), median degree=%d",
                                 radius_um_nb, r_nb_px, stats::median(deg)))
    # auto fallback if too sparse
    if (stats::median(deg) < 6L) {
      warning("[INVA5] Neighborhood too sparse; switching to KNN fallback (k_nb=20, k_den=40).")
      k_nb  <- 20L; k_den <- 40L
      kn_nb  <- FNN::get.knn(coords, k = k_nb)
      kn_den <- FNN::get.knn(coords, k = k_den)
      nb_list  <- lapply(seq_len(nrow(coords)), function(i) unique(c(i, kn_nb$nn.index[i, ])))
      den_list <- lapply(seq_len(nrow(coords)), function(i) unique(c(i, kn_den$nn.index[i, ])))
      mode <- "knn"
    }
  } else {
    k_nb  <- 20L; k_den <- 40L  # defaults
    kn_nb  <- FNN::get.knn(coords, k = k_nb)
    kn_den <- FNN::get.knn(coords, k = k_den)
    nb_list  <- lapply(seq_len(nrow(coords)), function(i) unique(c(i, kn_nb$nn.index[i, ])))
    den_list <- lapply(seq_len(nrow(coords)), function(i) unique(c(i, kn_den$nn.index[i, ])))
    if (verbose) message(sprintf("• Neighborhood by KNN: k_nb=%d, k_den=%d", k_nb, k_den))
  }
  
  # ---- 2) interface score (Tumor-Stroma mixing) ----
  interface_score <- numeric(nrow(coords))
  for (i in seq_len(nrow(coords))) {
    nb <- nb_list[[i]]
    pT <- mean(in_tumor[nb])
    pS <- mean(in_stroma[nb])
    interface_score[i] <- 4 * pT * pS  # in [0,1], max at pT=pS=0.5
  }
  
  # ---- 3) composite invasive score from modules ----
  # ---- 3) composite invasive score from modules ----
  # (你的 M 已经构建好了)
  if (isTRUE(rescale_modules)) {
    Mz <- scale(M)                    # 需要时才 z-score
  } else {
    Mz <- M                           # 你已在上游标准化，直接用
  }
  if (is.null(weights)) {
    comp_raw <- rowMeans(Mz, na.rm=TRUE)
  } else {
    w <- as.numeric(weights); w <- w / sum(w)
    comp_raw <- as.vector(Mz %*% w)
  }
  
  # spatial smoothing (mean over neighbors) with NaN guard
  comp_smooth <- comp_raw
  for (i in seq_len(nrow(coords))) {
    nb  <- nb_list[[i]]
    val <- mean(comp_raw[nb], na.rm = TRUE)
    if (is.nan(val)) val <- comp_raw[i]  # guard
    comp_smooth[i] <- val
  }
  
  # ---- 4) candidate boundary set & thresholding ----
  is_candidate <- interface_score >= interface_tau
  if (!any(is_candidate)) {
    warning("[INVA5] No candidates passed interface_tau; thresholding on all spots instead.")
    is_candidate <- rep(TRUE, length(interface_score))
  }
  thr <- stats::quantile(comp_smooth[is_candidate], probs = q_hot, na.rm = TRUE, names = FALSE)
  if (!is.finite(thr)) {
    warning("[INVA5] Threshold is not finite; using overall quantile.")
    thr <- stats::quantile(comp_smooth, probs = q_hot, na.rm = TRUE, names = FALSE)
  }
  hot_flag <- is_candidate & (comp_smooth >= thr)
  
  # ---- 5) connected components over hotspot graph ----
  comp_id <- integer(nrow(coords)); comp_id[] <- 0L
  if (any(hot_flag)) {
    map_old2new <- integer(nrow(coords)); map_old2new[which(hot_flag)] <- seq_len(sum(hot_flag))
    edges <- cbind(integer(0), integer(0))
    idx_hot <- which(hot_flag)
    for (ii in idx_hot) {
      nb <- nb_list[[ii]]
      nb_hot <- nb[hot_flag[nb]]
      if (length(nb_hot)) {
        a <- rep(map_old2new[ii], length(nb_hot))
        b <- map_old2new[nb_hot]
        edges <- rbind(edges, cbind(a, b))
      }
    }
    if (nrow(edges) > 0) {
      g  <- igraph::graph_from_edgelist(unique(t(apply(edges,1,sort))), directed = FALSE)
      cc <- igraph::components(g)$membership
      comp_id[hot_flag] <- cc
      # size filter
      keep <- names(which(table(cc) >= min_comp_size))
      keep <- as.integer(keep)
      comp_id[!(comp_id %in% keep)] <- 0L
      hot_flag[comp_id == 0L] <- FALSE
    }
  }
  
  # ---- 6) density score (fraction of hotspot nodes in larger radius) ----
  dens <- numeric(nrow(coords))
  for (i in seq_len(nrow(coords))) {
    dn <- den_list[[i]]
    p  <- mean(hot_flag[dn])
    if (is.nan(p)) p <- 0
    dens[i] <- p
  }
  
  # ---- 7) write-back ----
  if (write_meta) {
    object[[paste0(meta_prefix,"_iface")]]   <- interface_score
    object[[paste0(meta_prefix,"_comp")]]    <- comp_raw
    object[[paste0(meta_prefix,"_smooth")]]  <- comp_smooth
    object[[paste0(meta_prefix,"_flag")]]    <- hot_flag
    object[[paste0(meta_prefix,"_compID")]]  <- comp_id
    object[[paste0(meta_prefix,"_density")]] <- dens
  }
  
  # ---- 8) quick plots (safe, no .data in subset) ----
  make_hotspot_plot <- function(obj, density_col, title="Invasive Hotspot (INVA5 density)") {
    if (exists("PlotExpressionV2", mode="function")) {
      PlotExpressionV2(
        object = obj, feature = density_col, source="meta",
        shape="square", spot_size="auto", knn_smooth=FALSE,
        legend_title = title
      ) + ggplot2::ggtitle(title)
    } else {
      df <- obj@meta.data
      df$x <- if ("imagecol_scaled" %in% names(df)) df$imagecol_scaled else df$imagecol
      df$y <- if ("imagerow_scaled" %in% names(df)) -df$imagerow_scaled else -df$imagerow
      df$val <- df[[density_col]]
      ggplot2::ggplot(df, ggplot2::aes(x=x, y=y, fill=val)) +
        ggplot2::geom_point(shape=22, size=0.6, color=NA) +
        ggplot2::coord_fixed() + ggplot2::theme_void() +
        ggplot2::scale_fill_viridis_c() + ggplot2::ggtitle(title)
    }
  }
  make_outline_plot <- function(obj, flag_col, comp_col, title="Hotspot Components (outline)") {
    df <- obj@meta.data
    df$x <- if ("imagecol_scaled" %in% names(df)) df$imagecol_scaled else df$imagecol
    df$y <- if ("imagerow_scaled" %in% names(df)) -df$imagerow_scaled else -df$imagerow
    p <- ggplot2::ggplot(df, ggplot2::aes(x=x, y=y)) +
      ggplot2::geom_point(data = df[ !df[[flag_col]], ], color="grey85", size=0.4) +
      ggplot2::geom_point(data = df[  df[[flag_col]], ], color="#E41A1C", size=0.4) +
      ggplot2::coord_fixed() + ggplot2::theme_void() + ggplot2::ggtitle(title)
    
    comps <- unique(df[[comp_col]][ df[[flag_col]] & df[[comp_col]]>0 ])
    for (cid in comps) {
      sub <- df[df[[comp_col]]==cid & df[[flag_col]], c("x","y")]
      if (nrow(sub) < 5) next
      
      if (isTRUE(concave_outline) && requireNamespace("concaveman", quietly=TRUE)) {
        # ---- concave hull ----
        hull <- tryCatch({
          hh <- concaveman::concaveman(sub, concavity = concavity, length_threshold = length_thresh)
          colnames(hh) <- c("x","y"); as.data.frame(hh)
        }, error=function(e) NULL)
        if (!is.null(hull) && nrow(hull) >= 3) {
          p <- p + ggplot2::geom_polygon(data=hull, ggplot2::aes(x=x, y=y),
                                         inherit.aes = FALSE,
                                         color="#2C7FB8", fill="#2C7FB8",
                                         alpha=0.10, linewidth=0.6)
          next
        }
      }
      # ---- fallback: convex hull ----
      idx  <- grDevices::chull(sub$x, sub$y)
      hull <- sub[idx, ]
      p <- p + ggplot2::geom_polygon(data=hull, ggplot2::aes(x=x, y=y),
                                     inherit.aes = FALSE,
                                     color="#2C7FB8", fill="#2C7FB8",
                                     alpha=0.10, linewidth=0.6)
    }
    p
  }
  
  
  plist <- list(
    density   = make_hotspot_plot(object, paste0(meta_prefix,"_density")),
    outline   = make_outline_plot(object, paste0(meta_prefix,"_flag"), paste0(meta_prefix,"_compID")),
    composite = if (exists("PlotExpressionV2", mode="function")) {
      PlotExpressionV2(object, feature=paste0(meta_prefix,"_smooth"), source="meta",
                       shape="square", spot_size="auto", knn_smooth=FALSE,
                       legend_title="INVA5 composite (smoothed)") +
        ggplot2::ggtitle("INVA5 composite (smoothed)")
    } else NULL
  )
  
  invisible(list(
    object = object,
    params = list(radius_um_nb=radius_um_nb, radius_um_den=radius_um_den,
                  q_hot=q_hot, interface_tau=interface_tau, min_comp_size=min_comp_size,
                  modules=module_cols, weights=weights, mode=mode,
                  mpp_low = sc$mpp_low, r_nb_px = if (is.finite(sc$mpp_low)) radius_um_nb/sc$mpp_low else NA_real_),
    threshold = thr,
    plots   = plist
  ))
}







# 顶刊风格 100% 堆叠柱：非肿瘤带 × 非肿瘤细胞类型 组成
BandStack100 <- function(
    object,
    base_colors,
    band_col   = "PERI_band",
    type_col   = "Seurat_L1_label",
    bands      = c("0–50 µm","50–100 µm","100–150 µm",">150 µm"),  # 仅非肿瘤带
    tumor_type = "Tumor/Epithelial",
    title      = "Non-tumor composition across periphery bands",
    legend_ncol = 2,
    legend_y_spacing = 0.3,
    legend_key_h = 0.6,
    legend_key_w = 0.8,
    show_N = TRUE,                  # 是否在 x 轴下方标注 N
    N_y    = -0.035,                # N 标签的 y 位置（比例坐标）
    width  = 6.6, height = 4.2, dpi = 300,
    out_file = NULL                 # 如需保存，给文件名（pdf/png/svg）
){
  stopifnot(inherits(object, "Seurat"))
  stopifnot(band_col %in% colnames(object@meta.data))
  stopifnot(type_col %in% colnames(object@meta.data))
  
  # 统一破折号（避免 0-50/0–50 混用）
  normalize_band <- function(x){
    x <- gsub("0-50",   "0–50", x, fixed = TRUE)
    x <- gsub("50-100", "50–100", x, fixed = TRUE)
    x <- gsub("100-150","100–150", x, fixed = TRUE)
    x
  }
  
  df <- object@meta.data |>
    tibble::as_tibble(rownames = "cell_id") |>
    dplyr::transmute(
      band = normalize_band(as.character(.data[[band_col]])),
      type = as.character(.data[[type_col]])
    ) |>
    dplyr::filter(band %in% bands, type != tumor_type) |>
    dplyr::mutate(band = factor(band, levels = bands))
  
  # 计数与百分比
  tab_ct <- df |>
    dplyr::count(band, type, name = "n") |>
    dplyr::group_by(band) |>
    dplyr::mutate(frac = n / sum(n)) |>
    dplyr::ungroup()
  
  # 图例顺序按全局平均占比排序（阅读友好）
  type_order <- tab_ct |>
    dplyr::group_by(type) |>
    dplyr::summarise(mean_frac = mean(frac), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(mean_frac)) |>
    dplyr::pull(type)
  tab_ct$type <- factor(tab_ct$type, levels = type_order)
  
  # 颜色：用你的 base_colors，不足则补 HCL 颜色（恒亮度/饱和度）
  pal <- base_colors
  missing <- setdiff(levels(tab_ct$type), names(pal))
  if (length(missing) > 0) {
    hues <- seq(15, 375, length.out = length(missing) + 1)[1:length(missing)]
    extra <- grDevices::hcl(h = hues, c = 60, l = 65)
    names(extra) <- missing
    pal <- c(pal, extra)
  }
  pal_used <- pal[levels(tab_ct$type)]
  
  # 每带总数（x 轴下方的 N）
  band_totals <- df |>
    dplyr::count(band, name = "N") |>
    dplyr::mutate(label = paste0("N=", scales::comma(N)))
  
  # —— 绘图（无条内百分比标签）——
  p <- ggplot2::ggplot(tab_ct, ggplot2::aes(x = band, y = frac, fill = type)) +
    ggplot2::geom_col(width = 0.88, color = NA) +
    # x 轴下方 N（可关）
    { if (isTRUE(show_N))
      ggplot2::geom_text(data = band_totals, ggplot2::aes(x = band, y = N_y, label = label),
                         inherit.aes = FALSE, vjust = 1, size = 3.1) } +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                expand = ggplot2::expansion(mult = c(0.06, 0.02))) +
    ggplot2::scale_fill_manual(values = pal_used, drop = FALSE, name = "Cell type") +
    ggplot2::coord_cartesian(ylim = c(ifelse(show_N, N_y - 0.01, 0), 1)) +
    ggplot2::labs(x = NULL, y = "Composition (%)", title = title) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.line.y = ggplot2::element_line(linewidth = 0.3),
      axis.ticks.y = ggplot2::element_line(linewidth = 0.3),
      axis.text.x  = ggplot2::element_text(size = 10),
      legend.position = "right",
      legend.title    = ggplot2::element_text(face = "bold"),
      legend.key.height = grid::unit(legend_key_h, "lines"),
      legend.key.width  = grid::unit(legend_key_w, "lines"),
      legend.spacing.y  = grid::unit(legend_y_spacing, "lines")
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = legend_ncol, byrow = TRUE))
  
  if (!is.null(out_file)) {
    ggplot2::ggsave(out_file, p, width = width, height = height, dpi = dpi)
  }
  invisible(list(plot = p, table = tab_ct, totals = band_totals, palette = pal_used))
}




# --- Redraw tumor periphery with your PlotSpatialDistribution -----------------
# Computes distance-to-tumor (µm), band labels, and plots with PlotSpatialDistribution

# Computes distance-to-tumor (µm), band labels, and plots with PlotSpatialDistribution
RedrawPeripheryWithPSD <- function(
    object,
    domain_col     = "Seurat_L1_label",
    tumor_labels   = c("Tumor/Epithelial"),
    bands_um       = c(50, 100, 150),     # outward bands
    radius_um      = 50,                  # for local tumor-neighbor filter (filter only)
    min_neighbors  = 25,                  # avoid isolated cells
    use_scaled_coords    = TRUE,
    add_counts_in_legend = TRUE,
    
    # ---- Tumor core coloring rules ----
    core_subtype  = NULL,                                   # "Classical" / "Basal-like" to force; NULL = auto/grey fallback
    tumor_core_by = NULL,                                   # a meta column to infer subtype (e.g., "Tumor_Subtype")
    core_colors   = c("Classical" = "#1f78b4", "Basal-like" = "#e31a1c"),
    core_grey_if_unknown = "black",
    
    # ---- Fixed band palette with larger gradient (constant across subtypes) ----
    band_palette_fixed = c(
      "0–50 µm"    = "#fee08b",
      "50–100 µm"  = "#fe9929",
      "100–150 µm" = "#d95f02",
      ">150 µm"    = "grey70"
    ),
    
    
    # Optional full manual override; if provided, it wins entirely
    pal_peri = NULL,
    
    title = NULL,
    legend_title = "Distance to tumor"
){
  stopifnot(inherits(object, "Seurat"))
  md <- object@meta.data
  stopifnot(domain_col %in% colnames(md))
  
  # 1) Coordinates & scale: convert px -> µm
  xcol <- if (use_scaled_coords && "imagecol_scaled" %in% names(md)) "imagecol_scaled" else "imagecol"
  ycol <- if (use_scaled_coords && "imagerow_scaled" %in% names(md)) "imagerow_scaled" else "imagerow"
  stopifnot(all(c(xcol, ycol) %in% names(md)))
  
  mpp <- object@misc$microns_per_pixel %||%
    object@misc$spatial_scales$microns_per_pixel %||%
    object@misc$scales$microns_per_pixel %||%
    tryCatch(Images(object)[[1]]@scale.factors$microns_per_pixel, error=function(e) NULL)
  s_low <- object@misc$spatial_scales$tissue_lowres_scalef %||%
    object@misc$scales$tissue_lowres_scalef %||%
    tryCatch(Images(object)[[1]]@scale.factors$tissue_lowres_scalef, error=function(e) NULL)
  if (is.null(s_low) && all(c("imagecol_scaled","imagecol_hires") %in% colnames(md))) {
    ratio <- md$imagecol_scaled / md$imagecol_hires
    s_low <- stats::median(ratio[is.finite(ratio)], na.rm=TRUE)
  }
  if (is.null(mpp) || is.null(s_low)) stop("Missing spatial scales (microns_per_pixel / tissue_lowres_scalef).")
  px2um <- mpp / s_low
  
  XY_all   <- as.matrix(md[, c(xcol, ycol)])
  is_tumor <- md[[domain_col]] %in% tumor_labels
  XY_tumor <- XY_all[is_tumor, , drop = FALSE]
  if (nrow(XY_tumor) < 5) stop("Too few tumor points to compute periphery.")
  
  # 2) STRICT distance-based bands: nearest Euclidean distance to tumor (px) -> µm
  kn_all <- FNN::get.knnx(data = XY_tumor, query = XY_all, k = min(200, nrow(XY_tumor)))
  dist_um_nearest <- kn_all$nn.dist[, 1] * px2um
  
  # 3) Local tumor-neighbor count within radius_um (for filtering only; no effect on distance)
  radius_px <- radius_um / px2um
  if (requireNamespace("dbscan", quietly = TRUE)) {
    fr <- dbscan::frNN(x = XY_tumor, eps = radius_px, query = XY_all, sort = FALSE)
    neighborN <- vapply(fr$id, length, 1L)
  } else {
    k_adapt <- max(200, 5L * min_neighbors)
    k_use   <- min(nrow(XY_tumor), k_adapt)
    kn_ad   <- FNN::get.knnx(data = XY_tumor, query = XY_all, k = k_use)
    neighborN <- rowSums((kn_ad$nn.dist <= radius_px), na.rm = TRUE)
  }
  
  # 4) Band labeling for NON-tumor cells — solely by PERI_dist_um
  bands_um <- sort(unique(bands_um))
  brks  <- c(0, bands_um, Inf)
  labs  <- paste0(brks[-length(brks)], "–", brks[-1], " µm")
  labs[length(labs)] <- paste0(">", bands_um[length(bands_um)], " µm")
  
  band_non_tumor <- rep(labs[length(labs)], nrow(md))
  nt_idx <- which(!is_tumor)
  if (length(nt_idx)) {
    b <- cut(dist_um_nearest[nt_idx], breaks = brks, labels = labs,
             include.lowest = TRUE, right = TRUE)
    band_non_tumor[nt_idx] <- as.character(b)
    # isolate filter: downgrade sparse non-tumor to far band
    band_non_tumor[nt_idx][neighborN[nt_idx] < min_neighbors] <- labs[length(labs)]
  }
  
  band <- band_non_tumor
  band[is_tumor] <- "Tumor core"
  band_levels <- c("Tumor core", labs)
  band <- factor(band, levels = band_levels)
  
  # 5) Write meta
  object$PERI_dist_um    <- dist_um_nearest
  object$PERI_neighborN  <- neighborN
  object$PERI_band       <- band
  
  # 6) Build palette
  if (!is.null(pal_peri)) {
    pal_use <- pal_peri
  } else {
    # decide Tumor core color
    core_col <- core_grey_if_unknown
    if (!is.null(core_subtype) && core_subtype %in% names(core_colors)) {
      core_col <- core_colors[[core_subtype]]
    } else if (!is.null(tumor_core_by) && tumor_core_by %in% colnames(md)) {
      major <- suppressWarnings(names(sort(table(md[[tumor_core_by]][is_tumor]), decreasing = TRUE))[1])
      if (!is.na(major) && major %in% names(core_colors)) core_col <- core_colors[[major]]
    }
    # merge fixed bands with core color; respect level order
    pal_use <- c("Tumor core" = core_col, band_palette_fixed)
    # ensure keys cover all levels
    pal_use <- pal_use[band_levels]
  }
  
  # 7) Optional: add counts into legend labels
  band_to_use <- "PERI_band"
  if (add_counts_in_legend) {
    tab <- table(object$PERI_band)
    lvl <- levels(object$PERI_band)
    lvl_count <- sprintf("%s (%d)", lvl, as.integer(tab[lvl]))
    object$PERI_band_counted <- factor(
      sprintf("%s (%d)", as.character(object$PERI_band), as.integer(tab[as.character(object$PERI_band)])),
      levels = lvl_count
    )
    band_to_use <- "PERI_band_counted"
    pal_use <- setNames(unname(pal_use[lvl]), lvl_count)
  }
  
  # 8) Plot
  if (is.null(title)) {
    title <- sprintf("%s — Tumor periphery (bands: %s µm)",
                     if (!is.null(object@project.name)) object@project.name else "Sample",
                     paste(bands_um, collapse = "/"))
  }
  p <- PlotSpatialDistribution_bar(
    object,
    group_by   = band_to_use,
    spot_shape = "square",
    size_mode  = "auto",
    palette    = pal_use,
    title      = title,
    legend_title = legend_title,
    legend_spacing_y = 0.7, add_scale_bar        = T,scale_bar_length_mm  = 1,  
    legend_key_height = 0.9,
    legend_key_width  = 1.2
  )
  
  res_tab <- as.data.frame(table(object[[band_to_use]]), stringsAsFactors = FALSE)
  colnames(res_tab) <- c("band", "n")
  
  list(object = object, plot = p, table = res_tab)
}








## 先确保这些包已加载
suppressPackageStartupMessages({
  library(dplyr)
})

# ---- 把 RCTD 免疫类型优先写回的函数 ----
PinImmuneByRCTD <- function(
    object,
    seurat_col              = "Seurat_L1_label",   # 你的 Seurat 大类列
    rctd_type_col           = "RCTD_first",        # RCTD 一型结果
    rctd_second_col         = "RCTD_second",       # RCTD 二型结果（用于宽松纳入 doublet）
    rctd_class_col          = "RCTD_spot_class",   # singlet / doublet_certain / doublet_uncertain
    rctd_w1_col             = "RCTD_first_w",      # 第一类型权重
    rctd_margin_col         = "RCTD_margin",       # 第一-第二权重差
    immune_set              = c("B cells","T cells","Plasma cells","Myeloid","cDC","Mast"),
    # 置信阈值（可按需微调）
    w1_thresh_singlet       = 0.60,
    margin_thresh_singlet   = 0.25,
    w1_thresh_doublet       = 0.70,
    allow_doublet_if_second = c("Tumor","Fibroblast"),  # doublet 时，第二类型是这些且 w1 足够高也允许钉
    new_col                 = "cell_type_merged_immuneFirst",
    verbose                 = TRUE
){
  stopifnot(inherits(object, "Seurat"))
  md <- object@meta.data
  
  # 基础检查
  need_cols <- c(seurat_col, rctd_type_col, rctd_second_col, rctd_class_col, rctd_w1_col, rctd_margin_col)
  miss <- setdiff(need_cols, colnames(md))
  if (length(miss) > 0) stop("对象缺少必要列：", paste(miss, collapse=", "))
  
  # 取值并标准化为字符/数值
  type1   <- as.character(md[[rctd_type_col]])
  type2   <- as.character(md[[rctd_second_col]])
  sclass  <- as.character(md[[rctd_class_col]])
  w1      <- suppressWarnings(as.numeric(md[[rctd_w1_col]]))
  margin  <- suppressWarnings(as.numeric(md[[rctd_margin_col]]))
  fallback<- as.character(md[[seurat_col]])
  
  # 高置信条件
  is_singlet <- !is.na(sclass) & sclass == "singlet"
  is_doublet <- !is.na(sclass) & grepl("^doublet", sclass)
  
  high_conf_singlet <- is_singlet & is.finite(w1) & is.finite(margin) &
    (w1 >= w1_thresh_singlet) & (margin >= margin_thresh_singlet)
  
  high_conf_doublet <- is_doublet & is.finite(w1) &
    (w1 >= w1_thresh_doublet) & (!is.na(type2)) & (type2 %in% allow_doublet_if_second)
  
  high_conf <- high_conf_singlet | high_conf_doublet
  
  # 免疫优先：只有当 RCTD 一型是免疫集合 且 高置信 时，才覆盖
  immune_pin <- high_conf & !is.na(type1) & (type1 %in% immune_set)
  
  final <- fallback
  final[immune_pin] <- type1[immune_pin]
  
  # 写回对象
  object[[new_col]] <- factor(final)
  
  if (verbose) {
    cat("Pinned immune (n) by class:\n")
    print(table(type1[immune_pin], useNA = "ifany"))
    cat("\nSummary of", new_col, ":\n")
    print(sort(table(object[[new_col]]), decreasing = TRUE))
  }
  
  invisible(object)
}

# ==== ??/????/SoftGate ====







SubsetSeuratV5 <- function(object, cells, assays = NULL,
                           copy_reductions = TRUE,
                           copy_graphs = FALSE,
                           copy_images = TRUE,
                           keep_idents = TRUE) {
  cat(">>> Starting safe subsetting for Seurat v5...\n")
  stopifnot(inherits(object, "Seurat"))
  
  # cells 保留
  cells <- intersect(cells, colnames(object))
  if (length(cells) == 0) stop("No cells remaining after intersecting with object.")
  
  # assays 名称（字符向量）
  if (is.null(assays)) assays <- names(object@assays) else assays <- intersect(assays, names(object@assays))
  if (length(assays) == 0) stop("No assays to keep.")
  
  # counts：优先 layer='counts'，回退 slot='counts'
  counts_list <- lapply(assays, function(an) {
    mat <- tryCatch(
      GetAssayData(object, assay = an, layer = "counts"),
      error = function(e) GetAssayData(object, assay = an, slot = "counts")
    )
    mat[, cells, drop = FALSE]
  })
  names(counts_list) <- assays
  
  # 元数据
  metadata_subset <- object@meta.data[cells, , drop = FALSE]
  
  # 新对象
  new_obj <- CreateSeuratObject(counts = counts_list[[1]],
                                assay  = names(counts_list)[1],
                                meta.data = metadata_subset)
  if (length(counts_list) > 1) {
    for (i in 2:length(counts_list)) {
      new_obj[[names(counts_list)[i]]] <- CreateAssayObject(counts = counts_list[[i]])
    }
  }
  DefaultAssay(new_obj) <- DefaultAssay(object)
  cat("  - Base object with assays and metadata created for", ncol(new_obj), "cells.\n")
  
  # 降维
  if (isTRUE(copy_reductions)) {
    cat("  - Subsetting dimensional reductions...\n")
    for (reduc_name in Reductions(object)) {
      emb <- tryCatch(Embeddings(object, reduc_name), error = function(e) NULL)
      if (is.null(emb) || nrow(emb) == 0) next
      
      # 对齐到 new_obj 的细胞顺序；缺失填 NA
      new_cells <- Cells(new_obj)
      full_emb <- matrix(NA_real_, nrow = length(new_cells), ncol = ncol(emb),
                         dimnames = list(new_cells, colnames(emb)))
      idx <- match(new_cells, rownames(emb))
      hit <- which(!is.na(idx))
      if (length(hit) > 0) full_emb[hit, ] <- emb[idx[hit], , drop = FALSE]
      
      dr <- CreateDimReducObject(
        embeddings = full_emb,
        key   = tryCatch(Key(object[[reduc_name]]), error = function(e) paste0(toupper(reduc_name), "_")),
        assay = DefaultAssay(new_obj)
      )
      new_obj[[reduc_name]] <- dr
    }
  }
  
  # Graphs（默认不复制；若强制复制，做行列交集，并包裹 try）
  if (isTRUE(copy_graphs)) {
    cat("  - Subsetting graphs...\n")
    for (graph_name in Graphs(object)) {
      try({
        g <- object[[graph_name]]
        keep <- intersect(rownames(g), Cells(new_obj))
        if (length(keep) >= 2) {
          g_sub <- g[keep, keep, drop = FALSE]
          new_obj[[graph_name]] <- g_sub
        }
      }, silent = TRUE)
    }
  }
  
  # Images & misc
  if (isTRUE(copy_images) && length(Images(object)) > 0) {
    cat("  - Copying spatial image data...\n")
    new_obj@images <- object@images
  }
  if (length(object@misc) > 0) {
    new_obj@misc <- object@misc
  }
  
  # Idents 对齐
  if (isTRUE(keep_idents)) {
    try({
      ids <- Idents(object)             # factor, names = cell IDs
      ids <- ids[Cells(new_obj)]        # 严格按顺序取
      ids <- droplevels(ids)
      # 确保 names 完整对齐
      nm <- names(ids)
      if (is.null(nm) || !identical(nm, Cells(new_obj))) {
        ids <- setNames(as.vector(ids), Cells(new_obj))
      }
      Idents(new_obj) <- ids
    }, silent = TRUE)
  }
  
  # 校验
  tryCatch({
    stopifnot(identical(colnames(new_obj), rownames(new_obj@meta.data)))
    validObject(new_obj)
  }, error = function(e) warning("Subset object has minor inconsistencies: ", e$message))
  
  cat("✅ Safe subsetting complete.\n")
  new_obj
}





AbsorbRareLabels <- function(object,
                             label_col,
                             min_size = 50,
                             graph = NULL,
                             graph_candidates = c("Spatial_snn","SCT_snn","RNA_snn","integrated_snn","CCA_snn"),
                             k_pca = 30,
                             dims_use = 1:30,
                             reduction = "pca",
                             min_top_frac = 0.5,
                             exclude_rare_as_targets = TRUE,
                             protect_labels = NULL,
                             allowed_targets = NULL,
                             new_col = NULL,
                             overwrite = FALSE,
                             verbose = TRUE) {
  stopifnot(inherits(object, "Seurat"))
  stopifnot(label_col %in% colnames(object@meta.data))
  
  # --- 读取标签并统计 ---
  lab <- as.character(object@meta.data[[label_col]])
  tab <- sort(table(lab), decreasing = TRUE)
  rare_levels <- names(tab)[tab < min_size]
  if (length(rare_levels) == 0) {
    if (verbose) message("[Absorb] No rare labels (<", min_size, ") found. Nothing to do.")
    res_col <- if (isTRUE(overwrite)) label_col else (new_col %||% paste0(label_col, "_absorbed"))
    object@meta.data[[res_col]] <- factor(lab, levels = names(tab))
    return(list(object = object, summary = data.frame()))
  }
  
  # donor 候选池
  donor_levels <- setdiff(names(tab), rare_levels)
  if (!is.null(protect_labels)) donor_levels <- setdiff(donor_levels, protect_labels)
  if (!is.null(allowed_targets)) donor_levels <- intersect(donor_levels, allowed_targets)
  
  if (length(donor_levels) == 0) {
    stop("[Absorb] No valid target labels to absorb into (check protect/allowed settings).")
  }
  
  # --- 选择邻居来源：优先 SNN 图，回退 PCA ---
  use_graph <- NULL
  if (!is.null(graph)) {
    if (graph %in% Graphs(object)) use_graph <- graph else stop("[Absorb] graph '", graph, "' not found in object.")
  } else {
    cand <- intersect(graph_candidates, Graphs(object))
    if (length(cand) > 0) use_graph <- cand[1]
  }
  if (verbose) {
    message(if (!is.null(use_graph)) paste0("[Absorb] Using graph: ", use_graph)
            else paste0("[Absorb] No graph found; fallback to kNN on ", reduction))
  }
  
  # --- 预备 PCA 嵌入（仅在无图时用） ---
  emb <- NULL
  if (is.null(use_graph)) {
    if (!(reduction %in% Reductions(object))) {
      stop("[Absorb] reduction '", reduction, "' not found. Please compute ", reduction, " first or provide a graph.")
    }
    emb <- Embeddings(object, reduction)
    dims_use <- dims_use[dims_use <= ncol(emb)]
    if (length(dims_use) == 0) stop("[Absorb] dims_use out of range for reduction '", reduction, "'.")
    emb <- emb[, dims_use, drop = FALSE]
  }
  
  # --- 逐细胞计算目标标签 ---
  cells_all <- colnames(object)
  cells_rare <- which(lab %in% rare_levels)
  lab_new <- lab
  audit <- data.frame(cell=names(cells_rare), from=lab[cells_rare],
                      to=NA_character_, top_frac=NA_real_, stringsAsFactors = FALSE)
  
  # 预先算 donor 标签的全局大小（tie-break 用）
  donor_sizes <- tab[donor_levels]
  
  # helper: 从邻居标签计算加权多数
  choose_target <- function(nei_cells, weights, self_label) {
    if (length(nei_cells) == 0) return(c(target=NA_character_, frac=NA_real_))
    nei_labels <- lab[nei_cells]
    # 过滤不允许的目标
    ok <- nei_labels %in% donor_levels
    if (exclude_rare_as_targets) ok <- ok & !(nei_labels %in% rare_levels)
    if (!any(ok)) return(c(target=NA_character_, frac=NA_real_))
    
    nei_labels <- nei_labels[ok]
    ww <- if (is.null(weights)) rep(1, length(nei_labels)) else weights[ok]
    # 汇总到各目标标签
    agg <- tapply(ww, nei_labels, sum)
    agg <- agg[order(agg, decreasing = TRUE)]
    top_label <- names(agg)[1]
    top_frac  <- as.numeric(agg[1]) / sum(agg)
    
    # tie-break：若有平票，选全局更大的 donor
    ties <- names(agg)[agg == agg[1]]
    if (length(ties) > 1) {
      sizes <- donor_sizes[ties]; sizes[is.na(sizes)] <- 0
      top_label <- ties[which.max(sizes)]
    }
    c(target=top_label, frac=top_frac)
  }
  
  # 如果使用图
  G <- NULL
  if (!is.null(use_graph)) {
    G <- object[[use_graph]]
    # 确保行名/列名是全体细胞
    if (!all(cells_all %in% rownames(G)) || !all(cells_all %in% colnames(G))) {
      # 部分图（极少见），做对齐填补
      common <- intersect(cells_all, rownames(G))
      tmp <- Matrix::Matrix(0, nrow=length(cells_all), ncol=length(cells_all),
                            sparse = TRUE, dimnames = list(cells_all, cells_all))
      tmp[common, common] <- G[common, common]
      G <- tmp
    }
  }
  
  # 主循环（向量化也可做，这里为清晰）
  for (i in seq_along(cells_rare)) {
    idx <- cells_rare[i]
    cid <- cells_all[idx]
    
    if (!is.null(G)) {
      # 图邻居及权重（非零边）
      nz <- which(G[cid, ] != 0)
      nei <- names(nz)
      wts <- as.numeric(G[cid, nz])
      pick <- choose_target(nei_cells = nei, weights = wts, self_label = lab[idx])
    } else {
      # PCA kNN
      # 找到与该细胞距离最近的 k_pca 个邻居（不含自身）
      v <- emb[cid, , drop = FALSE]
      d2 <- rowSums((emb - matrix(v, nrow=nrow(emb), ncol=ncol(emb), byrow=TRUE))^2)
      ord <- order(d2, decreasing = FALSE)
      ord <- ord[ord != idx]          # 去掉自身
      nei_idx <- head(ord, k_pca)
      nei <- cells_all[nei_idx]
      wts <- NULL                     # 简单多数票
      pick <- choose_target(nei_cells = nei, weights = wts, self_label = lab[idx])
    }
    
    tgt <- unname(pick["target"]); fr <- as.numeric(pick["frac"])
    if (!is.na(tgt) && !is.na(fr) && fr >= min_top_frac) {
      lab_new[idx] <- tgt
      audit$to[i] <- tgt
      audit$top_frac[i] <- fr
    } else {
      # 不满足阈值：保留原标签
      audit$to[i] <- lab[idx]
      audit$top_frac[i] <- fr
    }
  }
  
  # --- 写回 ---
  res_col <- if (isTRUE(overwrite)) label_col else (new_col %||% paste0(label_col, "_absorbed"))
  # 统一 levels：包含旧的全部 + 新的 donor（通常同一集合）
  lvls <- union(names(tab), unique(lab_new))
  object@meta.data[[res_col]] <- factor(lab_new, levels = lvls)
  
  # 汇总表：rare 类被吸收去向
  sum_df <- aggregate(x=list(n=rep(1, nrow(audit))),
                      by=list(from=audit$from, to=audit$to), FUN=sum)
  sum_df <- sum_df[order(sum_df$from, -sum_df$n), ]
  
  if (verbose) {
    message("[Absorb] Rare labels: ", paste(rare_levels, collapse = ", "))
    message("[Absorb] Wrote new labels to: ", res_col)
  }
  list(object = object, summary = sum_df)
}





`%||%` <- function(a,b){ if(!is.null(a)) a else b }




AbsorbRareLabels <- function(object,
                             label_col,
                             min_size = 50,
                             graph = NULL,
                             graph_candidates = c("Spatial_snn","SCT_snn","RNA_snn","integrated_snn","CCA_snn"),
                             k_pca = 30,
                             dims_use = 1:30,
                             reduction = "pca",
                             min_top_frac = 0.5,
                             exclude_rare_as_targets = TRUE,
                             protect_labels = NULL,
                             allowed_targets = NULL,
                             new_col = NULL,
                             overwrite = FALSE,
                             # === 新增参数 ===
                             force_absorb = FALSE,         # 强制吸收，即使不满足 min_top_frac
                             fallback_strategy = c("largest_neighbor", "largest_global", "nearest_by_distance"),
                             verbose = TRUE,
                             return_details = FALSE) {     # 返回详细信息
  
  stopifnot(inherits(object, "Seurat"))
  stopifnot(label_col %in% colnames(object@meta.data))
  
  fallback_strategy <- match.arg(fallback_strategy)
  
  # --- 读取标签并统计 ---
  lab <- as.character(object@meta.data[[label_col]])
  tab <- sort(table(lab), decreasing = TRUE)
  rare_levels <- names(tab)[tab < min_size]
  
  if (length(rare_levels) == 0) {
    if (verbose) message("[Absorb] No rare labels (<", min_size, ") found. Nothing to do.")
    res_col <- if (isTRUE(overwrite)) label_col else (new_col %||% paste0(label_col, "_absorbed"))
    object@meta.data[[res_col]] <- factor(lab, levels = names(tab))
    return(list(object = object, summary = data.frame(), details = NULL))
  }
  
  if (verbose) {
    message("[Absorb] Found ", length(rare_levels), " rare label(s): ", paste(rare_levels, collapse = ", "))
    message("         Total cells in rare labels: ", sum(tab[rare_levels]))
  }
  
  # donor 候选池
  donor_levels <- setdiff(names(tab), rare_levels)
  if (!is.null(protect_labels)) donor_levels <- setdiff(donor_levels, protect_labels)
  if (!is.null(allowed_targets)) donor_levels <- intersect(donor_levels, allowed_targets)
  
  if (length(donor_levels) == 0) {
    stop("[Absorb] No valid target labels to absorb into (check protect/allowed settings).")
  }
  
  # --- 选择邻居来源：优先 SNN 图，回退 PCA ---
  use_graph <- NULL
  if (!is.null(graph)) {
    if (graph %in% Graphs(object)) {
      use_graph <- graph
    } else {
      stop("[Absorb] graph '", graph, "' not found in object.")
    }
  } else {
    cand <- intersect(graph_candidates, Graphs(object))
    if (length(cand) > 0) use_graph <- cand[1]
  }
  
  if (verbose) {
    message(if (!is.null(use_graph)) 
      paste0("[Absorb] Using graph: ", use_graph)
      else 
        paste0("[Absorb] No graph found; fallback to kNN on ", reduction))
  }
  
  # --- 预备 PCA 嵌入（仅在无图时用或需要距离计算时） ---
  emb <- NULL
  need_emb <- is.null(use_graph) || fallback_strategy == "nearest_by_distance"
  
  if (need_emb) {
    if (!(reduction %in% Reductions(object))) {
      stop("[Absorb] reduction '", reduction, "' not found. Please compute ", reduction, " first or provide a graph.")
    }
    emb <- Embeddings(object, reduction)
    dims_use <- dims_use[dims_use <= ncol(emb)]
    if (length(dims_use) == 0) {
      stop("[Absorb] dims_use out of range for reduction '", reduction, "'.")
    }
    emb <- emb[, dims_use, drop = FALSE]
  }
  
  # --- 逐细胞计算目标标签 ---
  cells_all <- colnames(object)
  cells_rare <- which(lab %in% rare_levels)
  lab_new <- lab
  
  audit <- data.frame(
    cell = cells_all[cells_rare],
    from = lab[cells_rare],
    to = NA_character_, 
    top_frac = NA_real_,
    method = NA_character_,  # 记录使用的方法
    stringsAsFactors = FALSE
  )
  
  # 预先算 donor 标签的全局大小（tie-break 用）
  donor_sizes <- tab[donor_levels]
  
  # helper: 从邻居标签计算加权多数
  choose_target <- function(nei_cells, weights, self_label) {
    if (length(nei_cells) == 0) {
      return(list(target = NA_character_, frac = NA_real_, method = "no_neighbors"))
    }
    
    nei_labels <- lab[nei_cells]
    
    # 过滤不允许的目标
    ok <- nei_labels %in% donor_levels
    if (exclude_rare_as_targets) ok <- ok & !(nei_labels %in% rare_levels)
    
    if (!any(ok)) {
      return(list(target = NA_character_, frac = NA_real_, method = "no_valid_neighbors"))
    }
    
    nei_labels <- nei_labels[ok]
    ww <- if (is.null(weights)) rep(1, length(nei_labels)) else weights[ok]
    
    # 汇总到各目标标签
    agg <- tapply(ww, nei_labels, sum)
    agg <- agg[order(agg, decreasing = TRUE)]
    top_label <- names(agg)[1]
    top_frac  <- as.numeric(agg[1]) / sum(agg)
    
    # tie-break：若有平票，选全局更大的 donor
    ties <- names(agg)[agg == agg[1]]
    if (length(ties) > 1) {
      sizes <- donor_sizes[ties]
      sizes[is.na(sizes)] <- 0
      top_label <- ties[which.max(sizes)]
    }
    
    list(target = top_label, frac = top_frac, method = "neighbor_vote")
  }
  
  # helper: fallback 策略
  apply_fallback <- function(idx, cid, self_label) {
    if (fallback_strategy == "largest_global") {
      # 分配给全局最大的 donor
      target <- names(donor_sizes)[which.max(donor_sizes)]
      return(list(target = target, frac = 0, method = "fallback_largest_global"))
      
    } else if (fallback_strategy == "largest_neighbor") {
      # 从邻居中找最常见的 donor（不管投票比例）
      if (!is.null(G)) {
        nz <- which(G[cid, ] != 0)
        nei <- names(nz)
      } else {
        v <- emb[cid, , drop = FALSE]
        d2 <- rowSums((emb - matrix(v, nrow = nrow(emb), ncol = ncol(emb), byrow = TRUE))^2)
        ord <- order(d2, decreasing = FALSE)
        ord <- ord[ord != idx]
        nei_idx <- head(ord, k_pca)
        nei <- cells_all[nei_idx]
      }
      
      if (length(nei) > 0) {
        nei_labels <- lab[nei]
        nei_labels <- nei_labels[nei_labels %in% donor_levels]
        if (length(nei_labels) > 0) {
          tab_nei <- table(nei_labels)
          target <- names(tab_nei)[which.max(tab_nei)]
          return(list(target = target, frac = 0, method = "fallback_largest_neighbor"))
        }
      }
      # 如果邻居中没有 donor，回退到全局最大
      target <- names(donor_sizes)[which.max(donor_sizes)]
      return(list(target = target, frac = 0, method = "fallback_largest_global"))
      
    } else if (fallback_strategy == "nearest_by_distance") {
      # 基于 PCA 距离找最近的非稀有类细胞
      if (is.null(emb)) {
        stop("Need embeddings for 'nearest_by_distance' strategy")
      }
      
      v <- emb[cid, , drop = FALSE]
      # 只考虑 donor 类的细胞
      donor_cells <- which(lab %in% donor_levels)
      if (length(donor_cells) == 0) {
        target <- names(donor_sizes)[which.max(donor_sizes)]
        return(list(target = target, frac = 0, method = "fallback_largest_global"))
      }
      
      emb_donors <- emb[donor_cells, , drop = FALSE]
      d2 <- rowSums((emb_donors - matrix(v, nrow = nrow(emb_donors), ncol = ncol(emb_donors), byrow = TRUE))^2)
      nearest_idx <- donor_cells[which.min(d2)]
      target <- lab[nearest_idx]
      
      return(list(target = target, frac = 0, method = "fallback_nearest_distance"))
    }
  }
  
  # 如果使用图
  G <- NULL
  if (!is.null(use_graph)) {
    G <- object[[use_graph]]
    # 确保行名/列名是全体细胞
    if (!all(cells_all %in% rownames(G)) || !all(cells_all %in% colnames(G))) {
      common <- intersect(cells_all, rownames(G))
      tmp <- Matrix::Matrix(0, nrow = length(cells_all), ncol = length(cells_all),
                            sparse = TRUE, dimnames = list(cells_all, cells_all))
      tmp[common, common] <- G[common, common]
      G <- tmp
    }
  }
  
  # 主循环
  n_forced <- 0
  n_voted <- 0
  
  for (i in seq_along(cells_rare)) {
    idx <- cells_rare[i]
    cid <- cells_all[idx]
    
    # 尝试邻居投票
    if (!is.null(G)) {
      nz <- which(G[cid, ] != 0)
      nei <- names(nz)
      wts <- as.numeric(G[cid, nz])
      pick <- choose_target(nei_cells = nei, weights = wts, self_label = lab[idx])
    } else {
      v <- emb[cid, , drop = FALSE]
      d2 <- rowSums((emb - matrix(v, nrow = nrow(emb), ncol = ncol(emb), byrow = TRUE))^2)
      ord <- order(d2, decreasing = FALSE)
      ord <- ord[ord != idx]
      nei_idx <- head(ord, k_pca)
      nei <- cells_all[nei_idx]
      wts <- NULL
      pick <- choose_target(nei_cells = nei, weights = wts, self_label = lab[idx])
    }
    
    tgt <- pick$target
    fr <- pick$frac
    method <- pick$method
    
    # 判断是否满足阈值
    if (!is.na(tgt) && !is.na(fr) && fr >= min_top_frac) {
      # 满足阈值，直接吸收
      lab_new[idx] <- tgt
      audit$to[i] <- tgt
      audit$top_frac[i] <- fr
      audit$method[i] <- method
      n_voted <- n_voted + 1
      
    } else if (force_absorb) {
      # 不满足阈值但强制吸收
      fallback <- apply_fallback(idx, cid, lab[idx])
      lab_new[idx] <- fallback$target
      audit$to[i] <- fallback$target
      audit$top_frac[i] <- fallback$frac
      audit$method[i] <- fallback$method
      n_forced <- n_forced + 1
      
    } else {
      # 保留原标签
      audit$to[i] <- lab[idx]
      audit$top_frac[i] <- fr
      audit$method[i] <- paste0("kept_", method)
    }
  }
  
  # --- 写回 ---
  res_col <- if (isTRUE(overwrite)) label_col else (new_col %||% paste0(label_col, "_absorbed"))
  lvls <- union(names(tab), unique(lab_new))
  object@meta.data[[res_col]] <- factor(lab_new, levels = lvls)
  
  # 汇总表
  sum_df <- aggregate(
    x = list(n = rep(1, nrow(audit))),
    by = list(from = audit$from, to = audit$to),
    FUN = sum
  )
  sum_df <- sum_df[order(sum_df$from, -sum_df$n), ]
  
  # 添加方法统计
  method_summary <- aggregate(
    x = list(n = rep(1, nrow(audit))),
    by = list(method = audit$method),
    FUN = sum
  )
  
  if (verbose) {
    message("[Absorb] Results:")
    message("         - Absorbed by voting: ", n_voted, " cells")
    if (force_absorb) {
      message("         - Force absorbed: ", n_forced, " cells")
    }
    message("         - Kept original: ", sum(audit$to == audit$from), " cells")
    message("[Absorb] Wrote new labels to: ", res_col)
  }
  
  # 返回
  result <- list(
    object = object, 
    summary = sum_df,
    method_stats = method_summary
  )
  
  if (return_details) {
    result$details <- audit
  }
  
  return(result)
}





step1_enrichment_and_forest <- function(
    object,
    # ==== 你可自由指定以下4个参数（默认给出示例） ====
    band_cols = c(
      "Classical"     = "PERI_band_Classical",
      "Basal-like"    = "PERI_band_Basal_like",
      "pEMT-invasive" = "PERI_band_pEMT_invasive"
    ),
    analyses_to_run = list(
      `Basal-like` = c("myCAF (COL11A1+/desmoplastic)", "SPP1 TAM"),
      `Classical`  = c("iCAF-like (CCL19+)",
                       "myCAF (COL11A1+/desmoplastic)",
                       "Stromal-adjacent")
    ),
    types_keep = c(
      "myCAF (ECM-rich)","myCAF (COL11A1+/desmoplastic)","iCAF-like (CCL19+)",
      "SPP1 TAM","TAM-Resident/Repair","TAM-ECM-Interface","TAM-interacting",
      "Endothelial","Perivascular/Smooth muscle",
      "B cells","Plasma cells","T cells",
      "Schwann/Glia","Stromal-adjacent"
    ),
    types_drop = c(
      "Classical","Basal-like","pEMT-invasive",
      "Ductal","Ductal/ADM","Ductal(PanIN-like)",
      "Acinar","Acinar-contam/doublet",
      "Endocrine","QC_low","LowQ/Ambig","Low-quality/ambiguous"
    ),
    # ==== 统计/绘图选项（可选） ====
    sample_col = NULL,
    n_perm = 5000,          # 置换次数
    B_boot = 5000,          # bootstrap 次数（细胞级）
    test = c("one.sided","two.sided"),  # 单侧(near>far)或双侧
    title = "Near vs Far enrichment (0–50µm vs ≥100µm)",
    x_lab  = "log2(Near/Far)",
    y_lab  = NULL,
    palette = NULL,         # 如 c("Basal-like"="#1f77b4","Classical"="#ff7f0e")
    save_plot_path = NULL,  # e.g. "Enrichment_Forest_Plot.pdf"
    width = 8, height = 6, dpi = 300
){
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(purrr)
  })
  
  test <- match.arg(test)
  normalize_dash  <- function(x) gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", x)
  normalize_label <- function(x) normalize_dash(trimws(as.character(x)))
  
  types_keep_norm <- normalize_label(types_keep)
  types_drop_norm <- normalize_label(types_drop)
  
  # --- 选样本列 ---
  if (is.null(sample_col)) {
    sample_col <- dplyr::case_when(
      "orig.ident" %in% colnames(object@meta.data) ~ "orig.ident",
      "sample"     %in% colnames(object@meta.data) ~ "sample",
      "Sample"     %in% colnames(object@meta.data) ~ "Sample",
      TRUE ~ NA_character_
    )
  }
  if (is.na(sample_col) || !sample_col %in% colnames(object@meta.data)) {
    stop("No valid sample id column found (tried orig.ident/sample/Sample).")
  }
  
  # --- 工具：规范带 ---
  normalize_band <- function(x) {
    x0 <- sub(" \\(.*\\)$", "", x)
    x0[grepl("^Tumor core", x0)] <- NA
    x1 <- ifelse(x0 %in% c("100–150 µm", ">150 µm"), "Far (≥100µm)",
                 ifelse(x0 == "0–50 µm", "Near (0-50µm)",
                        ifelse(x0 == "50–100 µm", NA, x0)))
    factor(x1, levels = c("Near (0-50µm)", "Far (≥100µm)"))
  }
  
  # --- 抽出“细胞级微环境表”（类型已归一化） ---
  get_microenv_df <- function(band_col, periphery_name){
    md <- object@meta.data
    if (!band_col %in% colnames(md)) {
      warning("Band column not found: ", band_col, " (skip ", periphery_name, ")")
      return(tibble())
    }
    md %>%
      tibble::as_tibble(rownames = "cell_id") %>%
      dplyr::select(
        cell_id,
        sample = all_of(sample_col),
        type_raw = Seurat_L2_refined,
        band_raw = all_of(band_col)
      ) %>%
      dplyr::filter(!is.na(band_raw), !grepl("^Tumor core", band_raw)) %>%
      dplyr::mutate(
        band_group = normalize_band(band_raw),
        type = normalize_label(type_raw)
      ) %>%
      dplyr::filter(
        !is.na(band_group),
        type %in% types_keep_norm,
        !type %in% types_drop_norm
      ) %>%
      dplyr::select(sample, band_group, type) %>%
      dplyr::mutate(periphery_of = periphery_name)
  }
  
  # --- 分层置换 + 细胞级 bootstrap 的近/远富集检验（log2 比的 CI 与 p） ---
  enrich_near_far <- function(df_sub, target_type,
                              n_perm = 5000, B_boot = 5000, eps = 1e-9, test = "one.sided"){
    if (nrow(df_sub) == 0) {
      return(tibble(
        target_cell_type = target_type,
        observed_log2fc_near_vs_far = NA_real_,
        ci95_low = NA_real_, ci95_high = NA_real_,
        p_value = NA_real_, n_perm = n_perm,
        n_near = 0, n_far = 0,
        frac_near = NA_real_, frac_far = NA_real_,
        periphery_of = unique(df_sub$periphery_of)[1]
      ))
    }
    tgt <- normalize_label(target_type)
    is_near <- df_sub$band_group == "Near (0-50µm)"
    is_far  <- df_sub$band_group == "Far (≥100µm)"
    is_tgt  <- df_sub$type == tgt
    n_near_tot <- sum(is_near); n_far_tot <- sum(is_far)
    n_near_tgt <- sum(is_near & is_tgt); n_far_tgt <- sum(is_far & is_tgt)
    
    if (n_near_tot == 0 || n_far_tot == 0) {
      return(tibble(
        target_cell_type = target_type,
        observed_log2fc_near_vs_far = NA_real_,
        ci95_low = NA_real_, ci95_high = NA_real_,
        p_value = NA_real_, n_perm = n_perm,
        n_near = n_near_tgt, n_far = n_far_tgt,
        frac_near = NA_real_, frac_far = NA_real_,
        periphery_of = unique(df_sub$periphery_of)[1]
      ))
    }
    
    frac_near <- n_near_tgt / n_near_tot
    frac_far  <- n_far_tgt  / n_far_tot
    obs_log2  <- log2((frac_near + eps)/(frac_far + eps))
    
    # Bootstrap CI（各带内细胞重抽 → 直接在 log2 比上建 CI）
    set.seed(1)
    near_idx <- which(is_near); far_idx <- which(is_far)
    boot_vals <- replicate(B_boot, {
      bn <- sample(near_idx, size = n_near_tot, replace = TRUE)
      bf <- sample(far_idx,  size = n_far_tot,  replace = TRUE)
      fn <- mean(is_tgt[bn]); ff <- mean(is_tgt[bf])
      log2((fn + eps)/(ff + eps))
    })
    ci <- stats::quantile(boot_vals, c(0.025, 0.975))
    
    # 分层置换（样本内打乱 Near/Far）
    set.seed(2)
    perm_vals <- replicate(n_perm, {
      dfp <- df_sub %>% dplyr::group_by(sample) %>% dplyr::mutate(band_group = sample(band_group)) %>% dplyr::ungroup()
      is_near_p <- dfp$band_group == "Near (0-50µm)"
      is_far_p  <- dfp$band_group == "Far (≥100µm)"
      fn <- mean(dfp$type[is_near_p] == tgt)
      ff <- mean(dfp$type[is_far_p]  == tgt)
      log2((fn + eps)/(ff + eps))
    })
    p_val <- if (test == "one.sided") {
      (sum(perm_vals >= obs_log2) + 1) / (n_perm + 1)
    } else {
      (sum(abs(perm_vals) >= abs(obs_log2)) + 1) / (n_perm + 1)
    }
    
    tibble(
      target_cell_type = target_type,
      observed_log2fc_near_vs_far = obs_log2,
      ci95_low = as.numeric(ci[1]),
      ci95_high = as.numeric(ci[2]),
      p_value = p_val,
      n_perm  = n_perm,
      n_near = n_near_tgt, n_far = n_far_tgt,
      frac_near = frac_near, frac_far = frac_far,
      periphery_of = unique(df_sub$periphery_of)[1]
    )
  }

  all_md <- purrr::imap_dfr(band_cols, ~ get_microenv_df(.x, .y))
  if (nrow(all_md) == 0) stop("No microenvironment cells remained after filtering.")
  
  res_list <- list()
  for (subtype in names(analyses_to_run)) {
    df_sub <- all_md %>% dplyr::filter(periphery_of == subtype)
    if (nrow(df_sub) == 0) next
    for (ct in analyses_to_run[[subtype]]) {
      res_list[[paste(subtype, ct, sep = " | ")]] <-
        enrich_near_far(df_sub, ct, n_perm = n_perm, B_boot = B_boot, test = test)
    }
  }
  final_tab <- dplyr::bind_rows(res_list) %>%
    dplyr::relocate(periphery_of, target_cell_type)
  

  plot_df <- final_tab %>%
    dplyr::mutate(
      y_lab_compact = paste0(periphery_of, " • ", target_cell_type),
      significance = dplyr::case_when(
        is.na(p_value)      ~ "na",
        p_value < 0.001     ~ "***",
        p_value < 0.01      ~ "**",
        p_value < 0.05      ~ "*",
        TRUE                ~ "ns"
      ),
      plot_log2fc = observed_log2fc_near_vs_far
    ) %>%
    dplyr::arrange(periphery_of, dplyr::desc(observed_log2fc_near_vs_far))
  
  p <- ggplot(plot_df, aes(x = plot_log2fc, y = reorder(y_lab_compact, plot_log2fc), color = periphery_of)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey60") +
    geom_errorbarh(aes(xmin = ci95_low, xmax = ci95_high),
                   width = 0.18, alpha = 0.7, linewidth = 0.6, show.legend = FALSE) + # ← 用 width
    geom_point(size = 3.2) +
    geom_text(aes(label = significance),
              nudge_y = -0.28, vjust = 1, hjust = 0.5, size = 3, show.legend = FALSE) + # 星标到误差线下方
    theme_minimal(base_size = 13) +
    theme(
      axis.title.y = element_text(margin = margin(r = 6)),
      panel.grid.major.y = element_line(linetype = "dotted", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.title = element_blank()
    ) +
    labs(title = title, x = x_lab, y = y_lab)
  
  if (!is.null(palette)) {
    p <- p + scale_color_manual(values = palette)
  }
  
  # 自适应 x 轴范围：含 CI
  x_min <- min(plot_df$ci95_low, plot_df$plot_log2fc, na.rm = TRUE)
  x_max <- max(plot_df$ci95_high, plot_df$plot_log2fc, na.rm = TRUE)
  pad <- 0.4
  p <- p + coord_cartesian(xlim = c(floor(x_min - pad), ceiling(x_max + pad)))
  
  if (!is.null(save_plot_path)) {
    ggsave(save_plot_path, p, width = width, height = height, dpi = dpi)
  }
  
  list(table = final_tab, plot = p)
}






# Slimmed “核心6 + 模块化扩展” 版（仅保留：myCAF(COL11A1+) 与 SPP1 TAM）
# 用法示例：
# res_step2 <- step2_sender_signals(
#   object,
#   n_perm = 5000, B_boot = 1000,
#   save_table_path = file.path(Config$Bankspath, "Step2_sender_signals_0-100.csv")
# )

step2_sender_signals <- function(
    object,
    band_cols = c(
      "Classical"     = "PERI_band_Classical",
      "Basal-like"    = "PERI_band_Basal_like",
      "pEMT-invasive" = "PERI_band_pEMT_invasive"
    ),
    periphery_senders = list(
      `Basal-like`    = c("myCAF (COL11A1+/desmoplastic)", "SPP1 TAM"),
      `pEMT-invasive` = c("myCAF (COL11A1+/desmoplastic)", "SPP1 TAM")
    ),

    sender_core6 = list(
      `myCAF (COL11A1+/desmoplastic)` = c( "TGFB1", "TGFB2", "TGFB3", # 检验所有的 TGFB 亚型
                                           "LOXL2", "ITGA11", "COL11A1", "POSTN", "SPARC", "FN1" # 检验所有关键的 ECM 基因
        
      ),
      `SPP1 TAM`                      = c( "SPP1", "TREM2", "APOE", "GPNMB", # 检验 CORE6 的成员
                                           "LGALS9", "VEGFA" # 额外检验我们感兴趣的免疫抑制和血管生成基因
                                          )
    ),
    sender_modules = list(
      `myCAF (COL11A1+/desmoplastic)` = list(
        TGFb_send     = c("TGFB1","TGFB2","THBS2","LTBP1","LTBP3"),
        ECM_crosslink = c("LOX","LOXL2","P4HA1","PXDN"),
        Integrin_mech = c("ITGA11","ITGA5","ITGB1","FN1"),
        Desmo_marker  = c("COL11A1","POSTN","SPARC")
      ),
      `SPP1 TAM` = list(
        Inhibitory    = c("IL10","CD274","PDCD1LG2","MRC1","VSIR"),
        Protease      = c("CTSB","CTSD","CTSS","MMP14"),
        Scavenger     = c("MSR1","STAB1","MARCO"),
        Lipid_meta    = c("LPL","LIPA","PLIN2","PPARG","APOC1")
      )
    ),
    types_keep = c(
      "myCAF (ECM-rich)","myCAF (COL11A1+/desmoplastic)",
      "SPP1 TAM","TAM-Resident/Repair","TAM-ECM-Interface","TAM-interacting",
      "Endothelial","Perivascular/Smooth muscle","B cells","Plasma cells","T cells",
      "Schwann/Glia","Stromal-adjacent"
    ),
    types_drop = c(
      "Classical","Basal-like","pEMT-invasive","Ductal","Ductal/ADM","Ductal(PanIN-like)",
      "Acinar","Acinar-contam/doublet","Endocrine","QC_low","LowQ/Ambig","Low-quality/ambiguous"
    ),
    assay = DefaultAssay(object),
    sample_col = NULL,
    # 统计设置
    n_perm = 10000,
    B_boot = 5000,
    det_threshold = 0,                       # >thr 记为检出
    near_bands = c("0–50 µm","0-50 µm"),     # 可改成含“50–100 µm”做 0–100 近带
    far_bands  = c("100–150 µm","100-150 µm",">150 µm"),
    # 输出
    save_table_path = NULL
){
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(tibble); library(purrr)
  })
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  

  if (is.null(sample_col)) {
    sample_col <- dplyr::case_when(
      "orig.ident" %in% colnames(object@meta.data) ~ "orig.ident",
      "sample"     %in% colnames(object@meta.data) ~ "sample",
      "Sample"     %in% colnames(object@meta.data) ~ "Sample",
      TRUE ~ NA_character_
    )
  }
  if (is.na(sample_col) || !sample_col %in% colnames(object@meta.data)) {
    stop("No valid sample id column found (tried orig.ident/sample/Sample).")
  }

  get_expr <- function(obj, assay){
    m <- tryCatch(GetAssayData(obj, assay = assay, layer = "data"), error = function(e) NULL)
    if (is.null(m)) m <- tryCatch(GetAssayData(obj, assay = assay, slot = "data"), error = function(e) NULL)
    if (is.null(m)) {
      m <- tryCatch(GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
      if (is.null(m)) stop("Cannot fetch expression from assay=", assay)
      m <- log1p(m)
    }
    m
  }
  mat <- get_expr(object, assay)  # gene x cell
  if (is.null(colnames(mat))) stop("Expression matrix has no cell column names.")
  

  to_near_far <- function(x, near_bands, far_bands){
    x0 <- sub(" \\(.*\\)$", "", x)
    out <- ifelse(x0 %in% near_bands, "Near", ifelse(x0 %in% far_bands, "Far", NA))
    factor(out, levels = c("Near","Far"))
  }
  

  get_microenv_cells <- function(band_col, periphery_name){
    md <- object@meta.data
    if (!band_col %in% colnames(md)) return(tibble())
    md %>%
      tibble::as_tibble(rownames = "cell_id") %>%
      dplyr::select(cell_id, sample = all_of(sample_col),
                    type = Seurat_L2_refined,
                    band_raw = all_of(band_col)) %>%
      dplyr::filter(!is.na(band_raw)) %>%
      dplyr::mutate(band_group = to_near_far(band_raw, near_bands, far_bands)) %>%
      dplyr::filter(!is.na(band_group),
                    type %in% types_keep, !type %in% types_drop) %>%
      dplyr::mutate(periphery_of = periphery_name) %>%
      dplyr::select(cell_id, sample, type, band_group, periphery_of)
  }

  test_near_far <- function(df, value_col = "value",
                            n_perm = 10000, B_boot = 5000, eps = 1e-9){
    stopifnot(all(c("sample","band_group", value_col) %in% colnames(df)))
    df <- df %>% dplyr::filter(band_group %in% c("Near","Far"))
    if (nrow(df) == 0) {
      return(list(near_mean=NA_real_, far_mean=NA_real_, log2fc=NA_real_, delta=NA_real_,
                  ci_low=NA_real_, ci_high=NA_real_, p_value=NA_real_, mode="NA",
                  n_samp_paired=0))
    }
    
    avg_sf <- df %>%
      dplyr::group_by(sample, band_group) %>%
      dplyr::summarise(mu = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = band_group, values_from = mu)
    
    n_samp <- sum(!is.na(avg_sf$Near) & !is.na(avg_sf$Far))
    
    if (n_samp >= 2) {
      di <- avg_sf$Near - avg_sf$Far; di <- di[!is.na(di)]
      obs_delta <- mean(di)
      
      set.seed(1)
      boot_vals <- replicate(B_boot, { idx <- sample.int(length(di), replace = TRUE); mean(di[idx]) })
      ci <- stats::quantile(boot_vals, c(0.025, 0.975))
      
      set.seed(2)
      perm_vals <- replicate(n_perm, { sgn <- sample(c(-1,1), length(di), replace = TRUE); mean(sgn * di) })
      p <- (sum(perm_vals >= obs_delta) + 1) / (n_perm + 1)
      
      near_mean <- mean(avg_sf$Near, na.rm = TRUE)
      far_mean  <- mean(avg_sf$Far,  na.rm = TRUE)
      list(near_mean=near_mean, far_mean=far_mean,
           log2fc = log2((near_mean + eps)/(far_mean + eps)),
           delta = obs_delta, ci_low=as.numeric(ci[1]), ci_high=as.numeric(ci[2]),
           p_value=p, mode="sample", n_samp_paired = n_samp)
    } else {
      near_mean <- mean(df[[value_col]][df$band_group == "Near"], na.rm = TRUE)
      far_mean  <- mean(df[[value_col]][df$band_group == "Far"],  na.rm = TRUE)
      obs_delta <- near_mean - far_mean
      
      set.seed(11)
      boot_vals <- replicate(B_boot, {
        df %>% dplyr::group_by(sample) %>%
          dplyr::summarise(
            mu_near = { v <- .data[[value_col]][band_group == "Near"]; if (length(v)==0) NA_real_ else mean(sample(v, length(v), replace = TRUE)) },
            mu_far  = { v <- .data[[value_col]][band_group == "Far"];  if (length(v)==0) NA_real_ else mean(sample(v, length(v), replace = TRUE)) },
            .groups = "drop"
          ) %>% dplyr::summarise(delta = mean(mu_near - mu_far, na.rm = TRUE)) %>% dplyr::pull(delta)
      })
      ci <- stats::quantile(boot_vals, c(0.025, 0.975))
      
      set.seed(22)
      perm_vals <- replicate(n_perm, {
        dfp <- df %>% dplyr::group_by(sample) %>% dplyr::mutate(band_group = sample(band_group)) %>% dplyr::ungroup()
        mn <- dfp %>% dplyr::filter(band_group == "Near") %>% dplyr::summarise(m=mean(.data[[value_col]], na.rm = TRUE)) %>% dplyr::pull(m)
        mf <- dfp %>% dplyr::filter(band_group == "Far")  %>% dplyr::summarise(m=mean(.data[[value_col]], na.rm = TRUE)) %>% dplyr::pull(m)
        (mn - mf)[1]
      })
      p <- (sum(perm_vals >= obs_delta) + 1) / (n_perm + 1)
      
      list(near_mean=near_mean, far_mean=far_mean,
           log2fc = log2((near_mean + eps)/(far_mean + eps)),
           delta = obs_delta, ci_low=as.numeric(ci[1]), ci_high=as.numeric(ci[2]),
           p_value=p, mode="cell", n_samp_paired = n_samp)
    }
  }
  

  all_cells <- purrr::imap_dfr(band_cols, ~ get_microenv_cells(.x, .y))
  if (nrow(all_cells) == 0) stop("No microenvironment cells remained after filtering.")
  if (is.null(all_cells$cell_id)) stop("cell_id missing in microenv table; check inputs.")
  

  results <- list()
  
  for (peri in names(periphery_senders)) {
    df_peri <- all_cells %>% dplyr::filter(periphery_of == peri)
    if (nrow(df_peri) == 0) next
    
    for (sender in periphery_senders[[peri]]) {
      df_sender_cells <- df_peri %>% dplyr::filter(type == sender)
      if (nrow(df_sender_cells) == 0) {
        results[[paste(peri, sender, sep=" | ")]] <- tibble(
          periphery_of = peri, sender = sender, feature = NA_character_, feature_type = NA_character_, module = NA_character_,
          near_mean = NA_real_, far_mean = NA_real_, log2_near_over_far = NA_real_,
          delta_near_minus_far = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_,
          p_value = NA_real_, p_adj = NA_real_, test_mode = NA_character_,
          det_near = NA_real_, det_far = NA_real_, log2DR = NA_real_,
          det_ci95_low = NA_real_, det_ci95_high = NA_real_, p_det = NA_real_, p_det_adj = NA_real_,
          n_cells_near = 0L, n_cells_far = 0L, n_samp_paired = 0L, note = "No cells"
        )
        next
      }
      
      # 元信息：Near/Far 细胞数
      n_cells_meta <- df_sender_cells %>% dplyr::count(band_group) %>%
        tidyr::pivot_wider(names_from = band_group, values_from = n, values_fill = 0)
      n_cells_near <- n_cells_meta$Near %||% 0L
      n_cells_far  <- n_cells_meta$Far  %||% 0L
      
      # 该 sender 的“核心6”和“模块列表”
      core6 <- (sender_core6[[sender]] %||% character(0)) %>% intersect(rownames(mat))
      modules <- sender_modules[[sender]] %||% list()
      
      # 当前 sender 的表达子矩阵（仅 sender 的细胞）
      cell_ids <- df_sender_cells$cell_id
      mat_sub  <- mat[, cell_ids, drop = FALSE]
      
      rows <- list()
      
      # —— 1) 核心6（CORE6；先按基因在“sender细胞中”做 z-score，再均值） —— #
      note_core <- NULL
      if (length(core6) > 0) {
        miss <- setdiff(sender_core6[[sender]], core6)
        if (length(miss) > 0) note_core <- paste("CORE6 missing:", paste(miss, collapse = ","))
        m_core <- mat_sub[core6, , drop = FALSE]
        # gene-wise z-score across sender cells
        m_core_z <- t(scale(t(m_core)))  # may produce NA if var=0
        if (anyNA(m_core_z)) {
          # fallback: center only
          m_core_z <- sweep(m_core, 1, rowMeans(m_core, na.rm = TRUE), FUN = "-")
          m_core_z[is.na(m_core_z)] <- 0
        }
        sig_vec <- colMeans(m_core_z)
        sig_tbl <- tibble(cell_id = names(sig_vec), value = as.numeric(sig_vec)) %>%
          dplyr::inner_join(df_sender_cells, by = "cell_id")
        st  <- test_near_far(sig_tbl, value_col = "value", n_perm = n_perm, B_boot = B_boot)
        std <- test_near_far(sig_tbl %>% dplyr::mutate(value = as.numeric(value > det_threshold)),
                             value_col = "value", n_perm = n_perm, B_boot = B_boot)
        rows[["CORE6"]] <- tibble(
          periphery_of = peri, sender = sender, feature = "CORE6", feature_type = "CORE6", module = NA_character_,
          near_mean = st$near_mean, far_mean = st$far_mean, log2_near_over_far = st$log2fc,
          delta_near_minus_far = st$delta, ci95_low = st$ci_low, ci95_high = st$ci_high,
          p_value = st$p_value, test_mode = st$mode, note = note_core %||% "",
          det_near = std$near_mean, det_far = std$far_mean, log2DR = std$log2fc,
          det_ci95_low = std$ci_low, det_ci95_high = std$ci_high, p_det = std$p_value,
          n_cells_near = n_cells_near, n_cells_far = n_cells_far, n_samp_paired = st$n_samp_paired
        )
      }

      if (length(modules) > 0) {
        for (mod_nm in names(modules)) {
          gp <- intersect(modules[[mod_nm]], rownames(mat_sub))
          miss <- setdiff(modules[[mod_nm]], gp)
          note_mod <- if (length(miss) > 0) paste0("MODULE missing:", paste(miss, collapse = ",")) else ""
          if (length(gp) == 0) next
          m_mod <- mat_sub[gp, , drop = FALSE]
          m_mod_z <- t(scale(t(m_mod)))
          if (anyNA(m_mod_z)) {
            m_mod_z <- sweep(m_mod, 1, rowMeans(m_mod, na.rm = TRUE), FUN = "-")
            m_mod_z[is.na(m_mod_z)] <- 0
          }
          sig_vec <- colMeans(m_mod_z)
          sig_tbl <- tibble(cell_id = names(sig_vec), value = as.numeric(sig_vec)) %>%
            dplyr::inner_join(df_sender_cells, by = "cell_id")
          st  <- test_near_far(sig_tbl, value_col = "value", n_perm = n_perm, B_boot = B_boot)
          std <- test_near_far(sig_tbl %>% dplyr::mutate(value = as.numeric(value > det_threshold)),
                               value_col = "value", n_perm = n_perm, B_boot = B_boot)
          rows[[paste0("MOD:", mod_nm)]] <- tibble(
            periphery_of = peri, sender = sender, feature = paste0("MODULE:", mod_nm), feature_type = "MODULE", module = mod_nm,
            near_mean = st$near_mean, far_mean = st$far_mean, log2_near_over_far = st$log2fc,
            delta_near_minus_far = st$delta, ci95_low = st$ci_low, ci95_high = st$ci_high,
            p_value = st$p_value, test_mode = st$mode, note = note_mod,
            det_near = std$near_mean, det_far = std$far_mean, log2DR = std$log2fc,
            det_ci95_low = std$ci_low, det_ci95_high = std$ci_high, p_det = std$p_value,
            n_cells_near = n_cells_near, n_cells_far = n_cells_far, n_samp_paired = st$n_samp_paired
          )
        }
      }

      if (length(core6) > 0) {
        for (g in core6) {
          vec <- as.numeric(mat_sub[g, , drop = TRUE])
          df_g <- tibble(cell_id = colnames(mat_sub), value = vec) %>%
            dplyr::inner_join(df_sender_cells, by = "cell_id")
          st  <- test_near_far(df_g, value_col = "value", n_perm = n_perm, B_boot = B_boot)
          std <- test_near_far(df_g %>% dplyr::mutate(value = as.numeric(value > det_threshold)),
                               value_col = "value", n_perm = n_perm, B_boot = B_boot)
          rows[[paste0("GENE:", g)]] <- tibble(
            periphery_of = peri, sender = sender, feature = g, feature_type = "GENE", module = NA_character_,
            near_mean = st$near_mean, far_mean = st$far_mean, log2_near_over_far = st$log2fc,
            delta_near_minus_far = st$delta, ci95_low = st$ci_low, ci95_high = st$ci_high,
            p_value = st$p_value, test_mode = st$mode, note = "",
            det_near = std$near_mean, det_far = std$far_mean, log2DR = std$log2fc,
            det_ci95_low = std$ci_low, det_ci95_high = std$ci_high, p_det = std$p_value,
            n_cells_near = n_cells_near, n_cells_far = n_cells_far, n_samp_paired = st$n_samp_paired
          )
        }
      }
      
      results[[paste(peri, sender, sep=" | ")]] <- dplyr::bind_rows(rows)
    }
  }
  
  res_tbl <- dplyr::bind_rows(results)
  

  if (nrow(res_tbl) > 0) {
    res_tbl <- res_tbl %>%
      dplyr::group_by(periphery_of, sender) %>%
      dplyr::mutate(
        p_adj     = p.adjust(p_value, method = "BH"),
        p_det_adj = p.adjust(p_det,    method = "BH")
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(periphery_of, sender,
                     dplyr::desc(feature_type == "CORE6"),
                     dplyr::desc(feature_type == "MODULE"),
                     dplyr::desc(log2_near_over_far))
  }
  
  if (!is.null(save_table_path)) {
    tryCatch({ utils::write.csv(res_tbl, file = save_table_path, row.names = FALSE) },
             error = function(e){ message("Failed to save CSV: ", e$message) })
  }
  return(res_tbl)
}







# =========================
# SoftGateHD: unified soft-gating for Visium HD
# =========================
# Dependencies: Seurat, UCell, dplyr

library(Seurat)
library(UCell)
library(dplyr)

# ---- Default robust signatures for Visium HD ----
.get_default_signatures <- function() {
  list(
    CAF = list(
      POS = c("COL1A1","COL1A2","COL3A1","DCN","LUM","THY1"),
      NEG = list(
        Immune = c("PTPRC","TYROBP","LYZ"),
        Endothelial = c("PECAM1","VWF","KDR"),
        Epithelial  = c("EPCAM","KRT8","KRT18","KRT19","KRT7"),
        Acinar      = c("PRSS1","CPA1","CTRB1"),
        Endocrine   = c("INS","GCG","SST","CHGA")
      )
    ),
    TAM = list(
      POS = c("PTPRC","LYZ","TYROBP","C1QA","C1QB","MS4A7","AIF1","FCGR3A","LST1"),
      NEG = list(
        Epithelial  = c("EPCAM","KRT8","KRT18","KRT19","KRT7"),
        Endothelial = c("PECAM1","VWF","KDR"),
        Fibro       = c("COL1A1","DCN"),
        Acinar      = c("PRSS1","CPA1"),
        Endocrine   = c("INS","CHGA")
      )
    ),
    EPI = list(
      POS = c("EPCAM","KRT8","KRT18","KRT19","KRT7","MUC1"),
      NEG = list(
        Immune      = c("PTPRC","TYROBP","LYZ"),
        Endothelial = c("PECAM1","VWF","KDR"),
        Fibro       = c("COL1A1","DCN"),
        Acinar      = c("PRSS1","CPA1","CTRB1"),
        Endocrine   = c("INS","GCG","SST","CHGA")
      )
    )
  )
}





# ---- Compact PDAC programs for tumor-like step ----
.get_pdac_programs <- function() {
  list(
    PDAC_Basal_like = c("KRT17","LAMC2","LAMB3","ITGA6","LGALS3","S100A2","ITGB4","COL17A1","TACSTD2","MSLN"),
    PDAC_Classical  = c("MUC1","KRT8","KRT18","EPCAM","AGR2","CEACAM6","KRT19","TFF1")
  )
}





# ---- Soft gate one class (POS vs NEG lists) ----
# Drop-in: replace your .soft_gate_one() with this version
.soft_gate_one <- function(obj, prefix, POS, NEG_lists,
                           q_low=0.25, q_high=0.75,
                           leak_weights=c(max=1.0, mean=0.0),
                           assay=NULL, slot="counts") {
  
  # 1) Name every panel so UCell creates predictable columns: <prefix>_<PANEL>_UCell
  feats <- c(
    setNames(list(POS), paste0(prefix, "_POS")),
    setNames(NEG_lists, paste0(prefix, "_NEG", seq_along(NEG_lists)))
  )
  
  # 2) Single-core (Windows safe); explicitly pass assay/slot
  obj <- UCell::AddModuleScore_UCell(
    obj, features = feats, ncores = 1, assay = assay, slot = slot
  )
  
  # 3) Resolve columns
  pos_col <- paste0(prefix, "_POS_UCell")
  if (!pos_col %in% colnames(obj@meta.data)) {
    stop("Missing UCell POS column: ", pos_col,
         ". Check gene symbols & assay/slot (e.g., assay='Spatial').")
  }
  neg_cols <- paste0(prefix, "_NEG", seq_along(NEG_lists), "_UCell")
  neg_cols <- neg_cols[neg_cols %in% colnames(obj@meta.data)]
  
  # 4) Leakage penalty robust to empty NEG
  if (length(neg_cols) == 0) {
    leak_max  <- rep(0, nrow(obj@meta.data))
    leak_mean <- rep(0, nrow(obj@meta.data))
  } else {
    neg_mat   <- as.matrix(obj@meta.data[, neg_cols, drop = FALSE])
    leak_max  <- apply(neg_mat, 1, max,  na.rm = TRUE)
    leak_mean <- apply(neg_mat, 1, mean, na.rm = TRUE)
  }
  leak  <- leak_weights["max"] * leak_max + leak_weights["mean"] * leak_mean
  score <- obj@meta.data[[pos_col]] - leak
  obj@meta.data[[paste0(prefix, "_score")]] <- score
  
  # 5) Quantile scaling → probability in [0,1]
  ql <- as.numeric(quantile(score, q_low,  na.rm = TRUE))
  qh <- as.numeric(quantile(score, q_high, na.rm = TRUE))
  den <- qh - ql
  p <- if (is.finite(den) && den > 0) (score - ql)/den else rep(0.5, length(score))
  p <- pmin(pmax(p, 0), 1)
  obj@meta.data[[paste0(prefix, "_p")]] <- p
  
  # 6) Tri-state
  lab <- ifelse(score >= qh, "Pure",
                ifelse(score <= ql, "Impure", "Ambiguous"))
  obj@meta.data[[paste0(prefix, "_soft")]] <- lab
  
  return(obj)
}






# ---- Optional neighborhood smoothing (majority vote for low-confidence only) ----
.smooth_labels_once <- function(obj, label_col, graph_name="RNA_snn", min_neighbors=5,
                                majority_frac=0.5, exclude_regex="^Ambig_", out_col=NULL) {
  if (is.null(out_col)) out_col <- paste0(label_col, "_smoothed")
  
  # If graph missing, quickly build a PCA kNN graph
  if (!(graph_name %in% names(obj@graphs))) {
    obj <- FindVariableFeatures(obj, nfeatures = 1500, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:30, k.param = 30, verbose = FALSE)
  }
  g <- obj@graphs[[graph_name]]
  lab <- obj@meta.data[[label_col]]
  
  # mark low-confidence cells
  low_idx <- which(lab %in% c("Uncertain", grep("^Ambig_", lab, value=TRUE)))
  if (length(low_idx) == 0) {
    obj@meta.data[[out_col]] <- lab
    return(obj)
  }
  
  for (i in low_idx) {
    nb <- names(which(g[i, ] > 0))
    if (length(nb) >= min_neighbors) {
      vote <- table(lab[nb])
      best <- names(which.max(vote))
      # Only adopt if clear majority and best is not an "Ambig_*"
      if (vote[best] >= majority_frac * sum(vote) && !grepl(exclude_regex, best)) {
        lab[i] <- best
      }
    }
  }
  obj@meta.data[[out_col]] <- lab
  return(obj)
}





# ---- (Optional) Resolve conflicts if multiple classes were gated together ----
.resolve_conflicts <- function(obj, class_prefixes, p_cut=0.6, delta_win=0.15, out_col="SoftGate_final") {
  # Choose class with highest p; if top-2 both ≥ p_cut and too close → Ambig_X+Y; if none ≥ p_cut → Uncertain
  P <- do.call(cbind, lapply(class_prefixes, function(px) obj@meta.data[[paste0(px,"_p")]]))
  colnames(P) <- class_prefixes
  res <- apply(P, 1, function(v){
    ord <- order(v, decreasing = TRUE)
    top <- colnames(P)[ord[1]]; p1 <- v[ord[1]]; p2 <- v[ord[2]]
    if (is.na(p1)) return("Uncertain")
    if (p1 >= 0.8 && (p1 - p2) >= delta_win) return(top)
    if (p1 >= p_cut && p2 >= p_cut && (p1 - p2) < delta_win) return(paste0("Ambig_", top, "+", colnames(P)[ord[2]]))
    if (p1 >= p_cut) return(top)
    return("Uncertain")
  })
  obj@meta.data[[out_col]] <- as.character(res)
  return(obj)
}





# ---- Main unified function ----
SoftGateHD <- function(obj,
                       targets = c("CAF","TAM","EPI","Tumor"),
                       use_default_signatures = TRUE,
                       custom_signatures = NULL,          # list like list(CAF=list(POS=..., NEG=list(...)), ...)
                       q_low = 0.25, q_high = 0.75,
                       leak_weights = c(max=1.0, mean=0.0),
                       do_conflict_resolution = FALSE,
                       smoothing = TRUE,
                       smoothing_graph = "RNA_snn",
                       smoothing_min_neighbors = 5,
                       smoothing_majority_frac = 0.5,
                       tumor_epi_p_min = 0.6,
                       tumor_quantiles = c(0.25, 0.75)) {
  
  sigs <- if (use_default_signatures || is.null(custom_signatures)) .get_default_signatures() else custom_signatures
  
  ran <- c()  # record which classes actually ran
  
  # CAF
  if ("CAF" %in% targets) {
    obj <- .soft_gate_one(obj, "CAF", sigs$CAF$POS, sigs$CAF$NEG, q_low, q_high, leak_weights)
    ran <- c(ran, "CAF")
  }
  # TAM
  if ("TAM" %in% targets) {
    obj <- .soft_gate_one(obj, "TAM", sigs$TAM$POS, sigs$TAM$NEG, q_low, q_high, leak_weights)
    ran <- c(ran, "TAM")
  }
  # EPI
  if ("EPI" %in% targets || "Tumor" %in% targets) {
    obj <- .soft_gate_one(obj, "EPI", sigs$EPI$POS, sigs$EPI$NEG, q_low, q_high, leak_weights)
    ran <- unique(c(ran, "EPI"))
  }
  
  # Optional conflict resolution when multiple classes requested
  if (do_conflict_resolution && length(ran) > 1) {
    obj <- .resolve_conflicts(obj, class_prefixes = ran, p_cut = 0.6, delta_win = 0.15, out_col = "SoftGate_final")
    if (smoothing) {
      obj <- .smooth_labels_once(obj, "SoftGate_final",
                                 graph_name = smoothing_graph,
                                 min_neighbors = smoothing_min_neighbors,
                                 majority_frac = smoothing_majority_frac,
                                 out_col = "SoftGate_final_smoothed")
    }
  }
  
  # Tumor two-step inside epithelial
  if ("Tumor" %in% targets) {
    pdac <- .get_pdac_programs()
    obj <- AddModuleScore_UCell(obj, features = pdac, name = names(pdac))
    
    epi_idx <- which(obj@meta.data$EPI_p >= tumor_epi_p_min)
    tum_p <- rep(NA_real_, nrow(obj@meta.data))
    if (length(epi_idx) > 0) {
      basal  <- obj@meta.data$PDAC_Basal_like_UCell[epi_idx]
      classc <- obj@meta.data$PDAC_Classical_UCell[epi_idx]
      maxprog <- pmax(basal, classc, na.rm = TRUE)
      tq <- quantile(maxprog, probs = tumor_quantiles, na.rm = TRUE)
      tum_p_local <- (maxprog - tq[1]) / (tq[2] - tq[1]); tum_p_local <- pmin(pmax(tum_p_local, 0), 1)
      tum_p[epi_idx] <- tum_p_local
    }
    obj$Tumor_like_p <- tum_p
    
    tum_lab <- rep("nonEpithelial_or_lowEPIprob", length(tum_p))
    tum_lab[which(!is.na(tum_p) & tum_p >= 0.75)] <- "Tumor_like_Pure"
    tum_lab[which(!is.na(tum_p) & tum_p <= 0.25)] <- "Epithelial_nonTumor"
    mid <- which(!is.na(tum_p) & tum_p > 0.25 & tum_p < 0.75)
    tum_lab[mid] <- "Tumor_like_Ambiguous"
    obj$Tumor_like_soft <- tum_lab
  }
  
  # Optional per-class smoothing (when not using joint conflict resolution)
  if (smoothing && (!do_conflict_resolution)) {
    for (px in ran) {
      obj <- .smooth_labels_once(obj,
                                 label_col = paste0(px,"_soft"),
                                 graph_name = smoothing_graph,
                                 min_neighbors = smoothing_min_neighbors,
                                 majority_frac = smoothing_majority_frac,
                                 out_col = paste0(px,"_soft_smoothed"))
    }
    if ("Tumor" %in% targets) {
      obj <- .smooth_labels_once(obj,
                                 label_col = "Tumor_like_soft",
                                 graph_name = smoothing_graph,
                                 min_neighbors = smoothing_min_neighbors,
                                 majority_frac = smoothing_majority_frac,
                                 out_col = "Tumor_like_soft_smoothed")
    }
  }
  
  return(obj)
}

# ==== Hotspot ? LR ?? ====












RunLR_Hotspot_vs_IF <- function(
    object,
    label_col        = "final_label",
    if_flag_col      = "InvasiveFront_flag",
    tumor_subtypes   = c("Classical", "Basal-like"),
    focus_cells      = c("Basal-like","myCAF","SPP1+ TAM"),
    unified_band_col = "PERI_band_Tumor",
    near_labels      = c("Tumor core", "0–50 µm", "50–100 µm"),
    bands_um         = c(50, 100, 150),
    radius_um        = 150,
    min_neighbors    = 15,
    min_cells        = 20,
    delta_thresh     = 0.1,
    out_dir          = NULL,
    liana_resource   = "OmniPath",
    liana_methods    = c("natmi","connectome","logfc","sca","cellphonedb")
){
  # Basic dependency checks
  if (!requireNamespace("liana", quietly = TRUE))  stop("Package 'liana' is required.")
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Package 'Seurat' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE))  stop("Package 'dplyr' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE))  stop("Package 'tidyr' is required.")
  if (!requireNamespace("readr", quietly = TRUE))  stop("Package 'readr' is required.")
  
  obj <- object
  md  <- obj@meta.data
  
  stopifnot(all(c(label_col, if_flag_col) %in% colnames(md)))
  if (is.null(out_dir)) {
    base_dir <- if (exists("Config") && !is.null(Config$Bankspath)) Config$Bankspath else "."
    out_dir  <- file.path(base_dir, "LR_Analysis_Hotspot_vs_IF")
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (!unified_band_col %in% colnames(md)) {
    if (!exists("RedrawPeripheryWithPSD")) {
      stop("RedrawPeripheryWithPSD() not found, cannot generate periphery bands.")
    }
    periphery_results <- RedrawPeripheryWithPSD(
      object        = obj,
      domain_col    = label_col,
      tumor_labels  = tumor_subtypes,
      bands_um      = bands_um,
      radius_um     = radius_um,
      min_neighbors = min_neighbors
    )
    obj[[unified_band_col]] <- periphery_results$object$PERI_band
    md <- obj@meta.data
  }
  
  # Clean band labels and define NearInterface
  clean_band <- sub(" \\(.*", "", md[[unified_band_col]])
  obj$PERI_band_Tumor_clean <- clean_band
  obj$NearInterface <- clean_band %in% near_labels
  
  obj[[if_flag_col]] <- as.logical(md[[if_flag_col]])
  obj$is_hotspot <- obj[[if_flag_col]] & obj$NearInterface
  obj$is_control <- obj[[if_flag_col]]            # whole IF as control
  
  # Helper: filter rare cell types
  filter_rare_celltypes <- function(seu, min_cells = 20) {
    if (is.null(seu) || ncol(seu) == 0) return(NULL)
    ct_counts <- table(Seurat::Idents(seu))
    keep      <- names(ct_counts[ct_counts >= min_cells])
    if (length(keep) < length(ct_counts)) {
      seu <- subset(seu, idents = keep)
    }
    if (length(unique(Seurat::Idents(seu))) < 2) return(NULL)
    seu
  }
  # Hotspot
  obj_hotspot <- subset(obj, subset = is_hotspot)
  if (ncol(obj_hotspot) > 0) {
    obj_hotspot$celltype_liana <- obj_hotspot[[label_col]]
    Seurat::Idents(obj_hotspot) <- "celltype_liana"
    obj_hotspot <- filter_rare_celltypes(obj_hotspot, min_cells)
  } else {
    obj_hotspot <- NULL
  }
  
  # Control (whole IF)
  obj_control <- subset(obj, subset = is_control)
  if (ncol(obj_control) > 0) {
    obj_control$celltype_liana <- obj_control[[label_col]]
    Seurat::Idents(obj_control) <- "celltype_liana"
    obj_control <- filter_rare_celltypes(obj_control, min_cells)
  } else {
    obj_control <- NULL
  }
  run_liana_one <- function(seu){
    if (is.null(seu)) return(NULL)
    res <- liana::liana_wrap(
      seu,
      resource = liana_resource,
      method   = liana_methods
    )
    agg <- liana::liana_aggregate(res)
    agg$score <- -log10(pmax(agg$aggregate_rank, 1e-6))
    agg
  }
  
  liana_hotspot <- run_liana_one(obj_hotspot)
  if (!is.null(liana_hotspot)) liana_hotspot$region <- "Hotspot"
  
  liana_control <- run_liana_one(obj_control)
  if (!is.null(liana_control)) liana_control$region <- "Control"
  
  if (is.null(liana_hotspot) || is.null(liana_control)) {
    warning("LIANA failed or too few cell types in hotspot/control. Returning NULL.")
    return(NULL)
  }

  liana_combined <- dplyr::bind_rows(liana_hotspot, liana_control)
  
  pair_deltas_all <- liana_combined %>%
    dplyr::select(source, target, ligand.complex, receptor.complex, region, score) %>%
    tidyr::pivot_wider(
      names_from  = region,
      values_from = score,
      values_fill = 0
    ) %>%
    dplyr::mutate(
      delta = Hotspot - Control,
      up_in = dplyr::case_when(
        delta >  delta_thresh ~ "Enriched in Hotspot Core",
        delta < -delta_thresh ~ "Enriched in IF Background",
        TRUE                  ~ "No Change"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(delta))
  
  readr::write_csv(
    pair_deltas_all,
    file.path(out_dir, "LIANA_Delta_HotspotCore_vs_IF_allPairs.csv")
  )
  
  pair_deltas_focus <- pair_deltas_all %>%
    dplyr::filter(
      source %in% focus_cells,
      target %in% focus_cells,
      up_in == "Enriched in Hotspot Core"
    )
  
  if (nrow(pair_deltas_focus) == 0) {
    warning("No hotspot-enriched LR pairs between requested focus_cells.")
    return(list(
      object            = obj,
      pair_deltas_all   = pair_deltas_all,
      pair_deltas_focus = pair_deltas_focus,
      triad_pairs       = NULL,
      axis_summary      = NULL,
      flow_axis         = NULL,
      out_dir           = out_dir
    ))
  }

  triad_pairs <- pair_deltas_focus %>%
    dplyr::mutate(
      lig      = toupper(ligand.complex),
      rec      = toupper(receptor.complex),
      axis     = axis_tag_plus(lig, rec),
      delta_pos = pmax(delta, 0)
    )
  
  axis_summary <- triad_pairs %>%
    dplyr::group_by(axis) %>%
    dplyr::summarise(
      n_pairs   = dplyr::n(),
      n_hot     = sum(up_in == "Enriched in Hotspot Core"),
      sum_delta = sum(delta_pos, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(sum_delta))
  
  flow_axis <- triad_pairs %>%
    dplyr::mutate(flow = paste(source, "→", target)) %>%
    dplyr::group_by(flow, axis) %>%
    dplyr::summarise(
      n_hot     = sum(up_in == "Enriched in Hotspot Core"),
      sum_delta = sum(delta_pos, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(sum_delta))
  
  readr::write_csv(
    pair_deltas_focus,
    file.path(out_dir, "LIANA_Delta_HotspotCore_vs_IF_focusCells.csv")
  )
  readr::write_csv(
    axis_summary,
    file.path(out_dir, "LIANA_AxisSummary_focusCells.csv")
  )

  invisible(list(
    object            = obj,
    pair_deltas_all   = pair_deltas_all,
    pair_deltas_focus = pair_deltas_focus,
    triad_pairs       = triad_pairs,
    axis_summary      = axis_summary,
    flow_axis         = flow_axis,
    out_dir           = out_dir
  ))
}









library(dplyr)
library(stringr)
library(dplyr)
library(stringr)

axis_tag_plus <- function(lig, rec){
  lig <- toupper(lig); rec <- toupper(rec)
  
  # --- 1. 基础正则定义 (Families) ---
  is_ecm_lig   <- stringr::str_detect(lig, "^(COL|FN1|LAMA|LAMB|LAMC|TNC|THBS\\d*|SPP1|POSTN|SPARC|BGN|VCAN|COMP|DCN|LUM|SPON2|SPON1|FBN\\d*|AGRN|LPL)")
  is_mtx_lig   <- stringr::str_detect(lig, "^(SPP1|THBS\\d*|CCN\\d*|POSTN|SPARC|LGALS\\d*|DCN|VCAN|MFGE8|SPON1|LPL|ANGPTL\\d*)")
  is_damp_lig  <- stringr::str_detect(lig, "^(S100A8|S100A9|HMGB1|BGN|VCAN)")
  is_wnt_lig   <- stringr::str_detect(lig, "^WNT")
  is_sema_lig  <- stringr::str_detect(lig, "^SEMA") | lig == "CD14"
  is_eph_lig   <- stringr::str_detect(lig, "^EFN") | lig == "EPHA2"
  is_apolipo   <- stringr::str_detect(lig, "^(APOA|APOB|APOC|APOE)")
  is_cytokine  <- stringr::str_detect(lig, "^(IL|IFN|TNF|LTA|OSM|LIF|CSF|MIF|MDK|PTN|GRN|CXCL\\d+|CCL\\d+)")
  is_mmp       <- stringr::str_detect(lig, "^MMP")
  is_plasmin   <- stringr::str_detect(lig, "^(PLAU|PLAT|SERPINE1|SERPINA1)")
  is_timp      <- stringr::str_detect(lig, "^TIMP")
  
  # --- 2. 受体定义 (Receptors) ---
  is_integrinR <- stringr::str_detect(rec, "^ITG[AB]")
  is_ddrR      <- rec %in% c("DDR1","DDR2")
  is_cd44_hspg <- rec %in% c("CD44","SDC1","SDC2","SDC3","SDC4","LRP1","LDLR","LRP8","VLDLR","APP","GPC1")
  is_erbbR     <- rec %in% c("EGFR","ERBB2","ERBB3","ERBB4")
  is_wntR      <- stringr::str_detect(rec, "^FZD") | rec %in% c("ROR1","ROR2","RYK","PTPRK","LRP5","LRP6","PTK7","ANTXR1","LDLR")
  is_ephR      <- stringr::str_detect(rec, "^EPH[AB]") | rec == "EFNA1"
  is_vegfR     <- rec %in% c("KDR","FLT1","FLT4","NRP1","NRP2", "FGFR1") 
  is_tnfr      <- stringr::str_detect(rec, "^TNFRSF|^TNFR")
  is_lair      <- rec %in% c("LAIR1", "LAIR2")
  is_dag       <- rec == "DAG1"
  is_tlrR      <- stringr::str_detect(rec, "^TLR") | rec %in% c("CD14","TREM1")
  is_tgfR      <- rec %in% c("TGFBR1","TGFBR2","ENG","ACVRL1","BMPR2") 
  
  # --- 3. ★★★ 新增：针对您 Unassigned 列表的强力抓取器 ★★★ ---
  
  # 四跨膜蛋白 (Tetraspanins) - 无论作配体还是受体
  tspan_list <- c("CD9","CD81","CD151","CD37","CD53","CD63","CD82","TM4SF1","TSPAN1","TSPAN15")
  is_tspan_interaction <- (lig %in% tspan_list | rec %in% tspan_list)
  
  # APP / CD74 / GPC1 轴
  app_list <- c("APP", "APLP2", "CD74", "GPC1", "CTSD")
  is_app_interaction <- (lig %in% app_list | rec %in% app_list)
  
  dplyr::case_when(
    # ============================================================
    # ★★★ 1. 优先处理您贴出的 Specific Pairs ★★★
    # ============================================================
    
    # 1. Tetraspanin Complex (覆盖 CD82, CD53, CD37, TM4SF1, CD151, CD63)
    is_tspan_interaction & (is_integrinR | stringr::str_detect(rec, "^CD") | stringr::str_detect(lig, "^CD") | lig == "DDR1" | lig == "TIMP1" | stringr::str_detect(lig, "^LAM")) ~ "Tetraspanin Complex / Adhesion",
    
    # 2. APP / CD74 / GPC1 Signaling
    is_app_interaction ~ "APP–CD74/HSPG Signaling",
    
    # 3. VEGF -> Integrin (非经典 VEGF)
    stringr::str_detect(lig, "^VEGF") & (is_integrinR | rec == "GPC1") ~ "VEGF–Integrin/HSPG Crosstalk",
    
    # 4. MMP -> Syndecan
    is_mmp & stringr::str_detect(rec, "^SDC") ~ "MMP–HSPG Crosstalk",
    
    # 5. RTK Heterodimerization (MET-ERBB2)
    lig == "MET" & rec == "ERBB2" ~ "RTK Heterodimerization (MET-ERBB2)",
    
    # 6. Receptor-Receptor Crosstalk (LDLR-LRP1)
    lig == "LDLR" & rec == "LRP1" ~ "ApoE/Lipoprotein uptake",
    
    # ============================================================
    # ★★★ 2. 之前的补丁 (保留防漏) ★★★
    # ============================================================
    stringr::str_detect(lig, "^ANGPTL") & stringr::str_detect(rec, "^SDC") ~ "Matricellular / CD44-LRP1-SDC",
    lig == "DCN" & is_erbbR                 ~ "EGFR/ERBB crosstalk",
    is_plasmin & is_integrinR               ~ "Plasminogen activation (PAI-1/uPA)",
    stringr::str_detect(lig, "^ITG") & rec == "THY1" ~ "Cell-Cell Adhesion (Integrin-Thy1)",
    (lig == "CD74" & rec == "CXCR4") | (lig == "CXCR4" & rec == "CD74") ~ "MIF–CD74/CXCR/CD44 axis",
    is_mmp & rec == "FGFR1"                 ~ "MMP–RTK crosstalk",
    lig == "CXCL12" & rec == "CXCR4"        ~ "CXCL12–CXCR4 Chemotaxis",
    lig == "CXCL12" & is_integrinR          ~ "Chemokine–Integrin crosstalk",
    is_timp & stringr::str_detect(rec, "^ADAM") ~ "Protease inhibition (TIMP–ADAM/MMP)",
    lig == "APOE" & rec == "LSR"            ~ "ApoE/Lipoprotein uptake",
    (lig == "CD44" & rec == "EPCAM") | (lig == "EPCAM" & rec == "CD44") ~ "Cell-Cell Adhesion (CD44/EpCAM)",
    rec == "CD47" & (stringr::str_detect(lig, "^THBS") | lig %in% c("COMP", "FN1")) ~ "CD47 checkpoint / Signal-regulatory",
    lig == "FN1" & rec == "PLAUR"             ~ "Plasminogen activation (PAI-1/uPA)",
    lig == "EDN1" & rec == "ECE1"             ~ "Endothelin processing (EDN1–ECE1)",
    (lig %in% c("PTPRF", "MET", "MST1R") & rec %in% c("MET", "MST1R")) ~ "MET/RON RTK Signaling",
    is_mmp & is_ddrR                          ~ "MMP–DDR crosstalk",
    is_mmp & rec == "CD44"                    ~ "MMP–CD44 ECM remodeling",
    lig == "EPHA2" & rec == "EFNA1"           ~ "Ephrin–Eph",
    lig == "POSTN" & rec == "PTK7"            ~ "WNT (non-canonical/PCP)",
    lig == "MET" & rec == "TNFRSF10B"         ~ "MET-TRAIL Apoptosis crosstalk",
    lig == "MDK" & is_integrinR               ~ "MDK–Integrin axis",
    lig %in% c("MET", "VEGFA") & rec == "CD44" ~ "RTK–CD44 Crosstalk",
    lig == "S100A4" & is_erbbR                ~ "S100A4–ERBB2 (Metastasis)",
    lig %in% c("ERBB2", "ERBB3", "ERBB4") & is_erbbR ~ "EGFR/ERBB crosstalk",
    lig == "ERBB3" & rec == "MUC1"            ~ "EGFR/ERBB crosstalk",
    lig == "FAM3C" & rec == "LAMP1"           ~ "FAM3C–LAMP1 trafficking",
    lig == "TGM2" & is_integrinR              ~ "ECM→Integrin",
    lig == "FCGR2A" & rec == "PDGFRB"         ~ "FcγR–PDGFRB crosstalk",
    
    # ============================================================
    # ★★★ 3. 标准逻辑 (Standard Families) ★★★
    # ============================================================
    stringr::str_detect(lig, "^COL") & is_lair                  ~ "Collagen–LAIR1 (Immune inhibition)",
    (stringr::str_detect(lig, "^LAMA") | lig == "AGRN") & is_dag ~ "ECM→Dystroglycan",
    is_ddrR      & stringr::str_detect(lig, "^COL")             ~ "Collagen→DDR (RTK)",
    is_integrinR & (is_ecm_lig | is_mtx_lig)                    ~ "ECM→Integrin",
    is_cd44_hspg & (is_mtx_lig | is_ecm_lig)                    ~ "Matricellular / CD44-LRP1-SDC",
    is_tlrR | is_damp_lig                                       ~ "Innate/DAMP (TLR2–CD14/TREM1)",
    stringr::str_detect(lig, "^MIF$") & rec %in% c("CD74","CXCR4","CXCR2","CD44") ~ "MIF–CD74/CXCR/CD44 axis",
    is_apolipo & (rec %in% c("LDLR","LRP1","VLDLR","LRP8","APP") | rec == "LSR") ~ "ApoE/Lipoprotein uptake",
    is_wnt_lig  & is_wntR                                       ~ "WNT (non-canonical/PCP)",
    is_sema_lig & (stringr::str_detect(rec, "^PLXN") | rec %in% c("NRP1","NRP2") | is_integrinR) ~ "Semaphorin signaling",
    is_eph_lig  & is_ephR                                       ~ "Ephrin–Eph",
    is_erbbR    & (is_cytokine | stringr::str_detect(lig, "^(EGF|EREG|AREG|HBEGF|TGFA|BTC|LAMB|LAMA|LAMB3)$")) ~ "EGFR/ERBB crosstalk",
    is_tgfR                                                     ~ "TGFβ/Angio/BMP co-receptor",
    is_vegfR                                                    ~ "VEGF signaling",
    
    TRUE ~ "Other/Unassigned"
  )
}







DetectHotspot <- function(
    object,
    feature_col,                 # e.g. "Score_ECM"
    mode             = c("knn", "radius"),
    k_smooth         = 25L,      # KNN neighbors for smoothing (if mode == "knn" or fallback)
    radius_um_smooth = 100,      # Smoothing radius in µm (if mode == "radius" and scale is available)
    q_high           = 0.90,     # High-quantile for hotspot threshold (0–1)
    cluster_eps_um   = 80,       # DBSCAN eps in µm (if scale available; otherwise treated as pixels)
    min_cluster_size = 30L,      # Minimal hotspot core size (number of cells)
    use_scaled_coords = TRUE,    # Prefer imagecol_scaled/imagerow_scaled if present
    write_meta       = TRUE,
    meta_prefix      = NULL,     # If NULL, use feature_col as prefix
    verbose          = TRUE
){
  # -------- small helper: null-coalescing --------
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  mode <- match.arg(mode)
  stopifnot(inherits(object, "Seurat"))
  md <- object@meta.data
  
  if (!feature_col %in% colnames(md)) {
    stop("feature_col not found in meta.data: ", feature_col)
  }
  feat_raw <- as.numeric(md[[feature_col]])
  
  # --- coordinates: scaled first, fallback to raw ---
  if (use_scaled_coords && all(c("imagecol_scaled","imagerow_scaled") %in% colnames(md))) {
    x <- md$imagecol_scaled
    y <- md$imagerow_scaled
  } else if (all(c("imagecol","imagerow") %in% colnames(md))) {
    x <- md$imagecol
    y <- md$imagerow
  } else {
    stop("Missing spatial coordinates: need either (imagecol_scaled,imagerow_scaled) or (imagecol,imagerow).")
  }
  coords <- cbind(as.numeric(x), as.numeric(y))
  rownames(coords) <- rownames(md)
  
  # --- get lowres microns-per-pixel if possible ---
  .get_lowres_mpp_local <- function(object){
    sc_misc <- object@misc %||% list()
    mpp_full <- sc_misc$microns_per_pixel %||%
      (sc_misc$spatial_scales$microns_per_pixel %||%
         (sc_misc$scales$microns_per_pixel %||% NA_real_))
    
    s_low <- sc_misc$spatial_scales$tissue_lowres_scalef %||%
      (sc_misc$scales$tissue_lowres_scalef %||% NA_real_)
    
    if (!is.finite(s_low)) {
      md <- object@meta.data
      if (all(c("imagecol_scaled","imagecol_hires") %in% colnames(md))) {
        ratio <- md$imagecol_scaled / md$imagecol_hires
        s_low <- stats::median(ratio[is.finite(ratio)], na.rm = TRUE)
      }
    }
    
    mpp_low <- if (is.finite(mpp_full) && is.finite(s_low)) mpp_full / s_low else NA_real_
    list(mpp_full = mpp_full, s_low = s_low, mpp_low = mpp_low)
  }
  
  sc <- .get_lowres_mpp_local(object)
  if (verbose && is.finite(sc$mpp_low)) {
    message(sprintf("Scale: mpp_full=%.4f, lowres_scalef=%.6f -> mpp_low=%.4f µm/px",
                    sc$mpp_full, sc$s_low, sc$mpp_low))
  }
  
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Please install 'FNN' package.")
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Please install 'dbscan' package.")
  }
  
  n <- nrow(coords)
  if (n < 5) stop("Too few cells for hotspot detection.")
  
  # --- 1) Build neighborhood for smoothing ---
  if (mode == "radius" && is.finite(sc$mpp_low)) {
    # radius mode in µm (converted to lowres pixels)
    r_px <- radius_um_smooth / sc$mpp_low
    if (verbose) {
      message(sprintf("Using radius-based smoothing: radius_um_smooth=%.1f µm (~%.1f px)",
                      radius_um_smooth, r_px))
    }
    fr <- dbscan::frNN(coords, eps = r_px, sort = FALSE)
    nb_list <- lapply(seq_len(n), function(i){
      nb <- fr$id[[i]]
      if (length(nb) == 0) i else unique(c(i, nb))
    })
    
    # automatic fallback if too sparse
    deg <- vapply(nb_list, length, 1L)
    if (stats::median(deg) < 6L) {
      warning("[DetectHotspot] Neighborhood too sparse with radius; switching to KNN smoothing (k_smooth).")
      kn <- FNN::get.knn(coords, k = k_smooth)
      nb_list <- lapply(seq_len(n), function(i) unique(c(i, kn$nn.index[i, ])))
    }
  } else {
    # KNN smoothing
    if (verbose) {
      message(sprintf("Using KNN-based smoothing: k_smooth=%d", k_smooth))
    }
    kn <- FNN::get.knn(coords, k = k_smooth)
    nb_list <- lapply(seq_len(n), function(i) unique(c(i, kn$nn.index[i, ])))
  }
  
  # --- 2) Smooth the feature over neighbors ---
  feat_smooth <- feat_raw
  for (i in seq_len(n)) {
    nb  <- nb_list[[i]]
    val <- mean(feat_raw[nb], na.rm = TRUE)
    if (is.nan(val)) val <- feat_raw[i]
    feat_smooth[i] <- val
  }
  
  # --- 3) Threshold by high quantile ---
  thr <- stats::quantile(feat_smooth, probs = q_high, na.rm = TRUE, names = FALSE)
  if (!is.finite(thr)) stop("Threshold is not finite; check feature values.")
  if (verbose) {
    message(sprintf("High-quantile threshold (q=%.2f): %.4f", q_high, thr))
  }
  
  high_flag <- feat_smooth >= thr & is.finite(feat_smooth)
  
  # --- 4) DBSCAN clustering within high-flag set ---
  idx_high <- which(high_flag)
  cluster_id <- integer(n); cluster_id[] <- 0L
  
  if (length(idx_high) > 0) {
    coords_high <- coords[idx_high, , drop = FALSE]
    
    # eps for DBSCAN: prefer µm if scale available, otherwise treat as pixels
    if (is.finite(sc$mpp_low)) {
      eps_px <- cluster_eps_um / sc$mpp_low
      if (verbose) {
        message(sprintf("DBSCAN clustering: eps=%.1f µm (~%.1f px), minPts=%d",
                        cluster_eps_um, eps_px, min_cluster_size))
      }
    } else {
      eps_px <- cluster_eps_um
      if (verbose) {
        message(sprintf("DBSCAN clustering: eps=%.1f (treated as pixels), minPts=%d",
                        eps_px, min_cluster_size))
      }
    }
    
    cl <- dbscan::dbscan(coords_high, eps = eps_px, minPts = min_cluster_size)
    cluster_id[idx_high] <- cl$cluster
  } else if (verbose) {
    message("No cells passed the high-quantile threshold; no clusters will be formed.")
  }
  
  # cluster > 0 means part of a connected hotspot core
  core_flag <- cluster_id > 0L
  
  # --- 5) Write back to meta.data ---
  if (is.null(meta_prefix)) meta_prefix <- feature_col
  
  if (write_meta) {
    object[[paste0(meta_prefix, "_raw")]]      <- feat_raw
    object[[paste0(meta_prefix, "_smooth")]]   <- feat_smooth
    object[[paste0(meta_prefix, "_highFlag")]] <- high_flag
    object[[paste0(meta_prefix, "_cluster")]]  <- cluster_id
    object[[paste0(meta_prefix, "_core")]]     <- core_flag
    
    misc_key <- paste0("Hotspot_", meta_prefix)
    object@misc[[misc_key]] <- list(
      feature_col      = feature_col,
      mode             = mode,
      threshold        = thr,
      n_cells          = n,
      n_high           = sum(high_flag),
      n_core           = sum(core_flag),
      mpp_full         = sc$mpp_full,
      mpp_low          = sc$mpp_low,
      radius_um_smooth = radius_um_smooth,
      cluster_eps_um   = cluster_eps_um,
      min_cluster_size = min_cluster_size
    )
  }
  
  if (verbose) {
    message(sprintf("Hotspot highFlag: %d / %d cells (%.1f%%)",
                    sum(high_flag), n, 100 * sum(high_flag)/n))
    message(sprintf("Hotspot core (cluster>0): %d / %d cells (%.1f%%)",
                    sum(core_flag), n, 100 * sum(core_flag)/n))
  }
  
  invisible(list(
    object  = object,
    threshold = thr,
    params = list(
      feature_col       = feature_col,
      mode              = mode,
      k_smooth          = k_smooth,
      radius_um_smooth  = radius_um_smooth,
      q_high            = q_high,
      cluster_eps_um    = cluster_eps_um,
      min_cluster_size  = min_cluster_size,
      use_scaled_coords = use_scaled_coords,
      mpp_low           = sc$mpp_low
    ),
    stats = list(
      n_cells      = n,
      n_high       = sum(high_flag),
      n_core       = sum(core_flag),
      frac_high    = sum(high_flag) / n,
      frac_core    = sum(core_flag) / n
    )
  ))
}










PlotExpressionV3<-function(
    object,                        # Seurat or data.frame
    feature,                       # gene or meta.data column
    assay        = "Spatial",
    slot         = "data",
    source       = c("auto","meta","matrix"),
    # display style
    shape        = c("square","circle"),
    spot_size    = c("auto","constant"),     # square+auto -> true tile size
    colors       = c("grey80","yellow","red"),
    legend_title = NULL,
    
    # ---- scale bar controls ----
    add_scale_bar        = FALSE,
    scale_bar_length_mm  = 1,
    scale_bar_n_segments = 4,
    scale_bar_color      = "navyblue",
    scale_bar_thickness_px = 4,
    scale_bar_position   = c("bottomright","bottomleft","topright","topleft"),
    scale_bar_margin_frac = c(0.04, 0.05),
    scale_bar_text       = NULL,            # NULL -> no label; e.g. "1 mm" to show
    scale_bar_text_size  = 3.2,
    scale_bar_text_color = "black",
    scale_bar_text_face  = "bold",
    scale_bar_text_offset_px = 10,
    scale_bar_gap_px     = 1,
    # physical thickness override (µm) to keep visual width comparable across samples
    scale_bar_thickness_um = NULL,
    
    # ---- optional HE background ----
    show_he        = FALSE,
    he_outs_path   = NULL,                  # e.g. Config$OutsPath
    he_resolution  = c("auto","hires","lowres"),
    clip_to_data   = FALSE,                 # clip view to data extent when HE is shown
    
    # ---- KNN smoothing (optional) ----
    knn_smooth   = FALSE,
    k            = 8,
    method       = c("gaussian","mean"),
    sigma        = NULL,
    # ---- physical scale helpers ----
    area_col     = NULL,                     # e.g. Cellbin "area_um2"
    area_unit    = c("um2","lowres_px"),
    bin_um       = NULL,                     # HD bin size (e.g. 8 µm)
    microns_per_pixel = NULL,                # scalar or column name
    scale_lowres      = NULL,                # scalar or column name
    ...
){
  # ---- small helper: null-coalescing ----
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # ---- 0) argument matching ----
  source       <- match.arg(source)
  shape        <- match.arg(shape)
  spot_size    <- match.arg(spot_size)
  method       <- match.arg(method)
  area_unit    <- match.arg(area_unit)
  scale_bar_position <- match.arg(scale_bar_position)
  he_resolution      <- match.arg(he_resolution)
  
  # ---- 1) fetch plotting data + scales ----
  coord_source <- "scaled"  # track whether we used scaled or raw coords
  
  if (inherits(object, "Seurat")) {
    md <- object@meta.data
    
    # coordinates (prefer lowres-scaled)
    if (all(c("imagecol_scaled","imagerow_scaled") %in% colnames(md))) {
      x <- md$imagecol_scaled; y <- md$imagerow_scaled
      coord_source <- "scaled"
    } else if (all(c("imagecol","imagerow") %in% colnames(md))) {
      x <- md$imagecol; y <- md$imagerow
      coord_source <- "raw"
    } else {
      stop("Cannot find spatial coordinates: need imagecol(_scaled)/imagerow(_scaled).")
    }
    
    plot_df <- data.frame(
      barcode = rownames(md),
      imagecol_scaled = as.numeric(x),
      imagerow_scaled = as.numeric(y),
      stringsAsFactors = FALSE
    )
    
    # auto scales (mpp / lowres scale / bin size)
    mpp <- microns_per_pixel %||%
      object@misc$microns_per_pixel %||%
      object@misc$spatial_scales$microns_per_pixel %||%
      object@misc$scales$microns_per_pixel %||%
      tryCatch(Images(object)[[1]]@scale.factors$microns_per_pixel, error=function(e) NULL)
    
    s_low <- scale_lowres %||%
      object@misc$spatial_scales$tissue_lowres_scalef %||%
      object@misc$scales$tissue_lowres_scalef %||%
      tryCatch(Images(object)[[1]]@scale.factors$tissue_lowres_scalef, error=function(e) NULL)
    
    if (is.null(bin_um)) {
      bin_um <- object@misc$bin_um %||%
        object@misc$spatial_scales$bin_size_um %||% NULL
    }
    
    # auto-detect Cellbin area column
    if (is.null(area_col) && "area_um2" %in% colnames(md)) area_col <- "area_um2"
    
    # feature vector
    if (source == "auto") {
      if (feature %in% colnames(md)) {
        val <- as.numeric(md[[feature]]); src <- "meta"
      } else {
        DefaultAssay(object) <- assay
        M <- Seurat::GetAssayData(object, assay=assay, slot=slot)
        if (feature %in% rownames(M)) {
          val <- as.numeric(M[feature, colnames(object)])
          src <- "matrix"
        } else {
          stop("Feature not found in meta.data or ", assay,"/",slot,": ", feature)
        }
      }
    } else if (source == "meta") {
      stopifnot(feature %in% colnames(md))
      val <- as.numeric(md[[feature]]); src <- "meta"
    } else { # matrix
      DefaultAssay(object) <- assay
      M <- Seurat::GetAssayData(object, assay=assay, slot=slot)
      stopifnot(feature %in% rownames(M))
      val <- as.numeric(M[feature, colnames(object)]); src <- "matrix"
    }
    
    # bring useful columns into df (area, etc.)
    if (!is.null(area_col) && area_col %in% colnames(md))
      plot_df[[area_col]] <- md[[area_col]]
    
  } else if (is.data.frame(object)) {
    plot_df <- object
    if (!all(c("imagecol_scaled","imagerow_scaled") %in% colnames(plot_df)))
      stop("data.frame mode: require columns imagecol_scaled, imagerow_scaled.")
    coord_source <- "scaled"
    
    if (feature %in% colnames(plot_df)) {
      val <- as.numeric(plot_df[[feature]]); src <- "meta"
    } else {
      stop("data.frame mode: please put '", feature, "' as a column or pass a Seurat object.")
    }
    
    mpp  <- microns_per_pixel
    s_low<- scale_lowres
    if (is.null(mpp) || is.null(s_low))
      warning("data.frame mode without microns_per_pixel / scale_lowres: scale bar and true tile size may be unavailable.")
  } else {
    stop("object must be a Seurat or data.frame.")
  }
  
  # ---- 2) optional KNN smoothing on Expression ----
  if (knn_smooth) {
    if (!requireNamespace("FNN", quietly = TRUE)) stop("Please install 'FNN' for knn_smooth.")
    if (!requireNamespace("Matrix", quietly = TRUE)) stop("Please install 'Matrix' for knn_smooth.")
    
    coords <- cbind(plot_df$imagecol_scaled, plot_df$imagerow_scaled)
    nn <- FNN::get.knn(coords, k = min(k, nrow(coords)-1))
    if (is.null(sigma) && method == "gaussian") {
      sigma <- stats::median(nn$nn.dist[,1], na.rm=TRUE)
      if (!is.finite(sigma) || sigma <= 0) sigma <- 1
    }
    # weights
    dist_mat <- cbind(0, nn$nn.dist)
    if (method == "gaussian") {
      Wrow <- exp(-(dist_mat^2)/(2*sigma^2))
    } else {
      Wrow <- matrix(1, nrow(dist_mat), ncol(dist_mat))
    }
    Wrow <- Wrow / pmax(rowSums(Wrow), 1e-12)
    
    n <- nrow(coords)
    idx_rows <- as.vector(row(Wrow))
    idx_cols <- as.vector(cbind(matrix(seq_len(n), n, 1), nn$nn.index))
    W <- Matrix::sparseMatrix(i=idx_rows, j=idx_cols, x=as.vector(Wrow), dims=c(n,n))
    val <- as.numeric(W %*% val)
    
    if (is.null(legend_title))
      legend_title <- paste0(feature, if (method=="gaussian") "_wknn" else "_knn", k)
  }
  
  # ---- 3) build plotting df & base extents ----
  plot_df$Expression <- val
  if (is.null(legend_title)) legend_title <- feature
  
  flip_y <- function(y) -y
  plot_df$y_plot <- flip_y(plot_df$imagerow_scaled)
  
  x_vals <- plot_df$imagecol_scaled
  y_vals <- plot_df$y_plot
  
  x_data_lim <- range(x_vals, na.rm = TRUE)
  y_data_lim <- range(y_vals, na.rm = TRUE)
  if (!all(is.finite(x_data_lim))) x_data_lim <- c(0, 1)
  if (!all(is.finite(y_data_lim))) y_data_lim <- c(0, 1)
  
  # default extents: data-driven
  x_lim <- x_data_lim
  y_lim <- y_data_lim
  
  # ---- 4) optional HE background ----
  he_img   <- NULL
  he_xmin  <- NA_real_
  he_xmax  <- NA_real_
  he_ymin  <- NA_real_
  he_ymax  <- NA_real_
  W_img    <- NA_real_
  H_img    <- NA_real_
  
  if (isTRUE(show_he)) {
    if (is.null(he_outs_path) || !dir.exists(he_outs_path)) {
      warning("show_he = TRUE but 'he_outs_path' is NULL or does not exist; skipping HE background.")
    } else {
      # decide resolution to use
      he_res_use <- he_resolution
      if (he_res_use == "auto") {
        he_res_use <- if (coord_source == "scaled") "lowres" else "hires"
      } else {
        # safeguard: if user picks a resolution inconsistent with coords, warn & auto-correct
        if (coord_source == "scaled" && he_res_use == "hires") {
          warning("Coordinates look lowres-scaled but he_resolution='hires'; switching to 'lowres'.")
          he_res_use <- "lowres"
        }
        if (coord_source == "raw" && he_res_use == "lowres") {
          warning("Coordinates look raw/full-res but he_resolution='lowres'; switching to 'hires'.")
          he_res_use <- "hires"
        }
      }
      
      img_file <- if (he_res_use == "lowres") {
        file.path(he_outs_path, "spatial", "tissue_lowres_image.png")
      } else {
        file.path(he_outs_path, "spatial", "tissue_hires_image.png")
      }
      
      if (!file.exists(img_file)) {
        warning("HE image not found at: ", img_file, "; skipping HE background.")
      } else {
        img <- png::readPNG(img_file)
        W_img <- ncol(img); H_img <- nrow(img)
        he_img  <- img
        he_xmin <- 0
        he_xmax <- W_img
        he_ymin <- -H_img
        he_ymax <- 0
        
        if (!isTRUE(clip_to_data)) {
          # full canvas = full HE image extent
          x_lim <- c(0, W_img)
          y_lim <- c(-H_img, 0)
        } else {
          # clip to data extent with small padding
          pad_x <- diff(x_data_lim) * 0.03
          pad_y <- diff(y_data_lim) * 0.03
          x_lim <- c(x_data_lim[1] - pad_x, x_data_lim[2] + pad_x)
          y_lim <- c(y_data_lim[1] - pad_y, y_data_lim[2] + pad_y)
        }
      }
    }
  } else {
    # no HE: add a small padding for nicer view
    pad_x <- diff(x_data_lim) * 0.03
    pad_y <- diff(y_data_lim) * 0.03
    x_lim <- c(x_data_lim[1] - pad_x, x_data_lim[2] + pad_x)
    y_lim <- c(y_data_lim[1] - pad_y, y_data_lim[2] + pad_y)
  }
  
  canvas_width  <- diff(x_lim)
  canvas_height <- diff(y_lim)
  
  # ---- 5) true tiles vs point symbols ----
  use_real_tiles <- (
    shape=="square" &&
      spot_size=="auto" &&
      ( (!is.null(area_col) && area_col %in% names(plot_df)) || !is.null(bin_um) )
  )
  
  if (use_real_tiles) {
    # vectorized scales (support scalar or column name)
    if (is.character(mpp)) {
      mpp_vec <- as.numeric(plot_df[[mpp]])
    } else {
      mpp_vec <- rep(as.numeric(mpp), nrow(plot_df))
    }
    if (is.character(s_low)) {
      sL_vec <- as.numeric(plot_df[[s_low]])
    } else {
      sL_vec <- rep(as.numeric(s_low), nrow(plot_df))
    }
    
    if (!is.null(area_col) && area_col %in% names(plot_df)) {
      # Cellbin: area -> side in lowres pixels
      if (area_unit=="um2") {
        side <- sqrt(pmax(as.numeric(plot_df[[area_col]]), 0)) * (sL_vec / mpp_vec)
      } else {
        side <- sqrt(pmax(as.numeric(plot_df[[area_col]]), 0))
      }
    } else {
      # HD: constant bin size
      if (is.null(bin_um)) stop("spot_size='auto' requires area_col or bin_um.")
      side <- (as.numeric(bin_um) / mpp_vec) * sL_vec
    }
    
    cx <- plot_df$imagecol_scaled
    cy <- plot_df$y_plot
    plot_df$.xmin <- cx - 0.5*side
    plot_df$.xmax <- cx + 0.5*side
    plot_df$.ymin <- cy - 0.5*side
    plot_df$.ymax <- cy + 0.5*side
    
    p <- ggplot2::ggplot(plot_df)
    
    # HE background if available
    if (!is.null(he_img)) {
      p <- p + ggplot2::annotation_raster(
        he_img,
        xmin = he_xmin, xmax = he_xmax,
        ymin = he_ymin, ymax = he_ymax,
        interpolate = TRUE
      )
    }
    
    p <- p +
      ggplot2::geom_rect(
        ggplot2::aes(xmin=.xmin, xmax=.xmax, ymin=.ymin, ymax=.ymax, fill=Expression),
        color = NA,
        ...
      ) +
      ggplot2::scale_fill_gradientn(
        colors = colors,
        name   = legend_title,
        na.value = "lightgrey"
      )
    
  } else {
    have_scattermore <- requireNamespace("scattermore", quietly = TRUE)
    
    p <- ggplot2::ggplot(plot_df)
    
    # HE background if available
    if (!is.null(he_img)) {
      p <- p + ggplot2::annotation_raster(
        he_img,
        xmin = he_xmin, xmax = he_xmax,
        ymin = he_ymin, ymax = he_ymax,
        interpolate = TRUE
      )
    }
    
    if (shape=="circle") {
      if (have_scattermore) {
        p <- p +
          scattermore::geom_scattermore(
            ggplot2::aes(x=imagecol_scaled, y=y_plot, color=Expression),
            pointsize = 2.5,
            pixels = c(2000, 2000),
            ...
          ) +
          ggplot2::scale_color_gradientn(
            colors   = colors,
            name     = legend_title,
            na.value = "lightgrey"
          )
      } else {
        p <- p +
          ggplot2::geom_point(
            ggplot2::aes(x=imagecol_scaled, y=y_plot, color=Expression),
            size = 2.5,
            shape = 16,
            stroke = 0,
            ...
          ) +
          ggplot2::scale_color_gradientn(
            colors   = colors,
            name     = legend_title,
            na.value = "lightgrey"
          )
      }
    } else { # square symbols (no true tiles)
      p <- p +
        ggplot2::geom_point(
          ggplot2::aes(x=imagecol_scaled, y=y_plot, fill=Expression),
          shape = 22,
          size  = 2.5,
          color = ggplot2::alpha("black", 0),
          stroke = 0.25,
          ...
        ) +
        ggplot2::scale_fill_gradientn(
          colors   = colors,
          name     = legend_title,
          na.value = "lightgrey"
        )
    }
  }
  
  # ---- 6) scale bar (in lowres pixels, with physical length/thickness if possible) ----
  if (isTRUE(add_scale_bar)) {
    # try to get scalar mpp and lowres scale
    if (inherits(object, "Seurat")) {
      if (is.character(mpp)) {
        v <- suppressWarnings(as.numeric(object@meta.data[[mpp]]))
        mpp_num <- v[which(is.finite(v))[1]]
      } else {
        mpp_num <- as.numeric(mpp)
      }
      if (is.character(s_low)) {
        v <- suppressWarnings(as.numeric(object@meta.data[[s_low]]))
        s_low_num <- v[which(is.finite(v))[1]]
      } else {
        s_low_num <- as.numeric(s_low)
      }
    } else {
      mpp_num    <- if (is.character(mpp)) NA_real_ else as.numeric(mpp)
      s_low_num  <- if (is.character(s_low)) NA_real_ else as.numeric(s_low)
    }
    
    if (!is.finite(mpp_num) || !is.finite(s_low_num)) {
      warning("Scale bar requested but microns_per_pixel / scale_lowres are not available; skipping scale bar.")
    } else {
      # µm per lowres pixel
      mpp_low   <- mpp_num / s_low_num
      px_per_um <- 1 / mpp_low
      
      bar_um <- scale_bar_length_mm * 1000
      bar_px <- as.integer(max(10, round(bar_um * px_per_um)))
      
      # thickness: prefer physical; fallback to px
      if (!is.null(scale_bar_thickness_um) && is.finite(scale_bar_thickness_um)) {
        th_px <- as.integer(max(1, round(scale_bar_thickness_um * px_per_um)))
      } else {
        th_px <- as.integer(max(1, round(scale_bar_thickness_px)))
      }
      
      mx <- if (length(scale_bar_margin_frac) >= 1) scale_bar_margin_frac[1] else 0.04
      my <- if (length(scale_bar_margin_frac) >= 2) scale_bar_margin_frac[2] else mx
      
      mx_px <- canvas_width  * mx
      my_px <- canvas_height * my
      
      # anchor
      if (scale_bar_position %in% c("bottomright","topright")) {
        x0 <- x_lim[2] - mx_px - bar_px
      } else {
        x0 <- x_lim[1] + mx_px
      }
      if (scale_bar_position %in% c("bottomright","bottomleft")) {
        y0 <- y_lim[1] + my_px
      } else {
        y0 <- y_lim[2] - my_px - th_px
      }
      x0 <- as.numeric(x0); y0 <- as.numeric(y0)
      x1 <- x0 + bar_px
      y1 <- y0 + th_px
      
      # segmentation with gaps
      nseg  <- max(1L, as.integer(round(scale_bar_n_segments)))
      gappx <- max(0L, as.integer(round(scale_bar_gap_px)))
      if (nseg == 1L) gappx <- 0L
      
      total_gap_px  <- gappx * (nseg - 1L)
      fill_px_total <- max(bar_px - total_gap_px, nseg)
      seg_fill_base <- floor(fill_px_total / nseg)
      rem           <- fill_px_total - seg_fill_base * nseg
      seg_fill      <- rep(seg_fill_base, nseg)
      if (rem > 0L) seg_fill[seq_len(rem)] <- seg_fill[seq_len(rem)] + 1L
      
      left <- x0
      seg_list <- vector("list", nseg)
      for (i in seq_len(nseg)) {
        right <- left + seg_fill[i]
        seg_list[[i]] <- data.frame(
          xmin = left, xmax = right,
          ymin = y0,   ymax = y1
        )
        left <- right + if (i < nseg) gappx else 0L
      }
      sb_df <- do.call(rbind, seg_list)
      
      p <- p +
        ggplot2::geom_rect(
          data = sb_df,
          ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill  = scale_bar_color,
          color = NA
        )
      
      # optional label
      if (!is.null(scale_bar_text) && nzchar(scale_bar_text)) {
        text_y <- if (scale_bar_position %in% c("bottomright","bottomleft")) {
          y1 + scale_bar_text_offset_px
        } else {
          y0 - scale_bar_text_offset_px
        }
        txt_df <- data.frame(
          x     = x0 + bar_px/2,
          y     = text_y,
          label = scale_bar_text
        )
        p <- p +
          ggplot2::geom_text(
            data = txt_df,
            ggplot2::aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size     = scale_bar_text_size,
            color    = scale_bar_text_color,
            fontface = scale_bar_text_face,
            vjust    = if (scale_bar_position %in% c("bottomright","bottomleft")) 0 else 1
          )
      }
    }
  }
  
  # ---- 7) final layout ----
  p <- p +
    ggplot2::coord_fixed(xlim = x_lim, ylim = y_lim) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust=0.5, face="bold"),
      legend.position = "right"
    )
  
  return(p)
}






library(Seurat)
library(dplyr)
library(dbscan)
library(FNN)
library(ggplot2)
library(Matrix) # ★★★ 必须加载这个包用于矩阵运算

DetectECMHotspot_MultiSample <- function(
    obj_list,                     
    feature_col = "Score_ECM",    
    celltype_col = "final_label", 
    target_celltypes = c("Basal-like", "Classical", "myCAF", "iCAF", "SPP1+ TAM", "FOLR2+ TAM"),
    
    radius_um_smooth = 120,       
    cluster_eps_um = 80,          
    min_cluster_size = 40,        
    seed = 123,                   
    out_dir = NULL                
) {
  
  # --- 内部辅助函数：获取坐标和分辨率 ---
  get_coords_and_scale <- function(object, obj_name) {
    md <- object@meta.data
    # 1. 坐标获取
    if (all(c("imagecol_scaled","imagerow_scaled") %in% colnames(md))) {
      coords <- cbind(md$imagecol_scaled, md$imagerow_scaled)
    } else if (all(c("imagecol","imagerow") %in% colnames(md))) {
      coords <- cbind(md$imagecol, md$imagerow)
    } else {
      stop(paste(obj_name, ": Missing spatial coordinates (imagecol/imagerow)."))
    }
    # 2. 分辨率获取
    sc_misc <- object@misc %||% list()
    mpp_full <- sc_misc$microns_per_pixel %||% (sc_misc$spatial_scales$microns_per_pixel %||% NA_real_)
    s_low <- sc_misc$spatial_scales$tissue_lowres_scalef %||% NA_real_
    if (!is.finite(s_low) && all(c("imagecol_scaled","imagecol_hires") %in% colnames(md))) {
      ratio <- md$imagecol_scaled / md$imagecol_hires
      s_low <- stats::median(ratio[is.finite(ratio)], na.rm = TRUE)
    }
    mpp_low <- if (is.finite(mpp_full) && is.finite(s_low)) mpp_full / s_low else NA_real_
    return(list(coords = coords, mpp_low = mpp_low))
  }
  
  # =========================================================================
  # Step 1: 预计算平滑评分 (极速矩阵版)
  # =========================================================================
  message(">>> Step 1: Calculating Smoothed Scores (Matrix Acceleration)...")
  
  obj_list <- lapply(names(obj_list), function(nm) {
    message(paste("   Processing:", nm, "..."))
    obj <- obj_list[[nm]]
    if(!feature_col %in% colnames(obj@meta.data)) stop(paste(nm, ": Feature missing."))
    
    raw_vals <- obj[[feature_col]][,1]
    # 将 NA 替换为 0，防止矩阵运算崩溃
    raw_vals[is.na(raw_vals)] <- 0
    
    info <- get_coords_and_scale(obj, nm)
    n_cells <- nrow(info$coords)
    
    # --- 核心优化：使用邻接表构建稀疏矩阵 ---
    ids_list <- list()
    
    if(is.finite(info$mpp_low)) {
      r_px <- radius_um_smooth / info$mpp_low
      # dbscan frNN 极快
      fr <- dbscan::frNN(info$coords, eps = r_px, sort = FALSE)
      ids_list <- fr$id
    } else {
      warning(paste(nm, ": Scale missing, using KNN=25."))
      kn <- FNN::get.knn(info$coords, k=25)
      # 将 KNN 矩阵转为 list
      ids_list <- split(kn$nn.index, row(kn$nn.index))
    }
    
    # ★★★ 矩阵化魔法 (Matrix Magic) ★★★
    # 1. 构建邻接矩阵 (Adjacency Matrix)
    # 这是一个 n x n 的大矩阵，如果 i 和 j 是邻居，则位置为 1
    
    # 展开 list 为 i, j 索引
    lens <- lengths(ids_list)
    # j 是邻居的索引
    j_idx <- unlist(ids_list)
    # i 是当前细胞的索引 (根据邻居数量重复)
    i_idx <- rep(seq_len(n_cells), times = lens)
    
    # 添加"自己"到邻居列表 (Smoothing 包含自己)
    j_idx <- c(j_idx, seq_len(n_cells))
    i_idx <- c(i_idx, seq_len(n_cells))
    
    # 创建稀疏矩阵 (避免内存爆炸)
    # 这里的 1 代表权重。如果是高斯平滑，这里可以放高斯权重。我们用均值，所以是 1。
    Adj <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = 1, dims = c(n_cells, n_cells))
    
    # 2. 归一化 (Row Normalization)
    # 每一行的和代表该细胞有多少个邻居。我们需要除以这个数来求平均。
    row_sums <- Matrix::rowSums(Adj)
    row_sums[row_sums == 0] <- 1 # 防止除以0
    
    # 3. 矩阵乘法求平滑值
    # Smoothed = (A * Raw) / N_neighbors
    # 这一步在 C++ 底层完成，极快
    smooth_vals <- as.numeric((Adj %*% raw_vals) / row_sums)
    
    obj$temp_smooth_score <- smooth_vals
    obj$temp_sample_id <- nm 
    return(obj)
  })
  names(obj_list) <- names(obj_list) 
  
  # =========================================================================
  # Step 2: 提取数据、过滤细胞、平衡抽样
  # =========================================================================
  message(">>> Step 2: Balanced Sampling...")
  
  all_data_list <- lapply(obj_list, function(obj) {
    obj@meta.data %>% 
      dplyr::select(dplyr::all_of(c("temp_smooth_score", celltype_col))) %>%
      dplyr::mutate(Sample = obj$temp_sample_id)
  })
  
  df_all <- dplyr::bind_rows(all_data_list)
  
  # 强制类型转换，防止因子报错
  df_all$Sample <- as.character(df_all$Sample)
  
  # 过滤
  df_target <- df_all %>% dplyr::filter(!!rlang::sym(celltype_col) %in% target_celltypes)
  
  if(nrow(df_target) == 0) stop("No target cells found!")
  
  # 检查
  counts <- table(df_target$Sample)
  min_count <- min(counts)
  message("   Target cells per sample: ")
  print(counts)
  
  # 抽样
  set.seed(seed)
  df_balanced <- df_target %>%
    dplyr::group_by(Sample) %>%
    dplyr::sample_n(min_count) %>%
    dplyr::ungroup()
  
  # =========================================================================
  # Step 3: 全局阈值 (K-means)
  # =========================================================================
  message(">>> Step 3: Determining Global Threshold...")
  
  valid_scores <- df_balanced$temp_smooth_score
  valid_scores <- valid_scores[is.finite(valid_scores)]
  
  km <- kmeans(valid_scores, centers = 2, nstart = 25)
  centers <- sort(km$centers)
  GLOBAL_THRESHOLD <- mean(centers)
  
  message("-----------------------------------------------------------")
  message(sprintf("★ GLOBAL THRESHOLD: %.4f", GLOBAL_THRESHOLD))
  message("-----------------------------------------------------------")
  
  if (!is.null(out_dir)) {
    p <- ggplot(df_balanced, aes(x = temp_smooth_score, fill = Sample)) +
      geom_density(alpha = 0.4) +
      geom_vline(xintercept = GLOBAL_THRESHOLD, color = "red", linetype = "dashed") +
      theme_minimal() +
      labs(title = "Global Threshold (Balanced)", x = "Score")
    try(ggsave(file.path(out_dir, "Global_ECM_Threshold.pdf"), p, width=6, height=4))
  }
  
  # =========================================================================
  # Step 4: 回溯应用 & DBSCAN
  # =========================================================================
  message(">>> Step 4: Final Clustering (DBSCAN)...")
  
  final_obj_list <- lapply(obj_list, function(obj) {
    
    feat_smooth <- obj$temp_smooth_score
    high_flag <- feat_smooth >= GLOBAL_THRESHOLD
    
    info <- get_coords_and_scale(obj, obj$temp_sample_id)
    cluster_id <- integer(length(high_flag)); cluster_id[] <- 0L
    
    idx_high <- which(high_flag)
    if (length(idx_high) > 0) {
      coords_high <- info$coords[idx_high, , drop=FALSE]
      eps_px <- if(is.finite(info$mpp_low)) cluster_eps_um / info$mpp_low else cluster_eps_um
      
      # dbscan 也是 C++ 写的，处理 50万个点通常只需要几秒
      cl <- dbscan::dbscan(coords_high, eps = eps_px, minPts = min_cluster_size)
      cluster_id[idx_high] <- cl$cluster
    }
    
    obj$ECMhot_smooth   <- feat_smooth
    obj$ECMhot_highFlag <- high_flag
    obj$ECMhot_core     <- cluster_id > 0L
    obj$ECMhot_cluster  <- cluster_id
    
    obj$temp_smooth_score <- NULL
    obj$temp_sample_id <- NULL
    
    return(obj)
  })
  
  message(">>> Done.")
  return(final_obj_list)
}







library(liana)
library(tidyverse)
library(Seurat)

RunLR_Merged_Hotspot_vs_Control_v2 <- function(
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
    liana_methods    = c("natmi","connectome","logfc","sca","cellphonedb")
){
  
  # --- 1. 设置输出路径 ---
  if (is.null(out_dir)) {
    # 尝试从 Config 获取路径，如果没有则使用默认
    base_dir <- if (exists("Config") && !is.null(Config$ECMFront)) Config$ECMFront else "."
    out_dir  <- file.path(base_dir, "LIANA_Merged_Analysis")
  }
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  message(">>> Starting LIANA analysis on Merged Object...")
  
  # ★★★ 修复 1: Seurat V5 图层合并 ★★★
  # 检查是否有 split layers, 如果有则合并，否则 LIANA 找不到数据
  if (any(grepl("\\.", Layers(object)))) {
    message("   Detected split layers (Seurat V5). Joining layers now...")
    object <- JoinLayers(object)
  }
  
  message(">>> Comparison: ", group_col, " ", group_hotspot, " (Hotspot) vs ", group_control, " (Control)")
  
  # --- 2. 准备数据与过滤 ---
  
  # ★★★ 修复 2: 更安全的子集提取方法 ★★★
  # 获取 Hotspot 细胞 ID
  cells_hotspot <- colnames(object)[object[[group_col]] == group_hotspot]
  # 获取 Control 细胞 ID
  cells_control <- colnames(object)[object[[group_col]] == group_control]
  
  # 移除 NA
  cells_hotspot <- na.omit(cells_hotspot)
  cells_control <- na.omit(cells_control)
  
  message(paste0("   Hotspot cells found: ", length(cells_hotspot)))
  message(paste0("   Control cells found: ", length(cells_control)))
  
  if (length(cells_hotspot) == 0 || length(cells_control) == 0) {
    stop("Error: No cells found for one of the groups based on ", group_col)
  }
  
  # 辅助函数：过滤稀有细胞
  filter_rare_celltypes <- function(seu, min_n) {
    if (is.null(seu) || ncol(seu) == 0) return(NULL)
    ct_counts <- table(Idents(seu))
    keep <- names(ct_counts[ct_counts >= min_n])
    if (length(keep) == 0) return(NULL)
    subset(seu, idents = keep)
  }
  
  # 提取 Hotspot 对象
  obj_hotspot <- subset(object, cells = cells_hotspot)
  Idents(obj_hotspot) <- label_col
  obj_hotspot <- filter_rare_celltypes(obj_hotspot, min_cells)
  
  # 提取 Control 对象
  obj_control <- subset(object, cells = cells_control)
  Idents(obj_control) <- label_col
  obj_control <- filter_rare_celltypes(obj_control, min_cells)
  
  if (is.null(obj_hotspot) || is.null(obj_control)) {
    stop("Error: After filtering for min_cells, one group became empty.")
  }
  
  # --- 3. 运行 LIANA ---
  run_liana_one <- function(seu){
    message("   Running LIANA on subset: ", length(colnames(seu)), " cells...")
    res <- liana_wrap(seu, resource = liana_resource, method = liana_methods)
    agg <- liana_aggregate(res)
    agg$score <- -log10(pmax(agg$aggregate_rank, 1e-6)) 
    return(agg)
  }
  
  # 分别运行
  message(">>> Processing Hotspot Group...")
  liana_hotspot <- run_liana_one(obj_hotspot)
  liana_hotspot$region <- "Hotspot"
  
  message(">>> Processing Control Group...")
  liana_control <- run_liana_one(obj_control)
  liana_control$region <- "Control"
  
  # --- 4. 计算 Delta (差异) ---
  message(">>> Calculating Delta Scores...")
  liana_combined <- bind_rows(liana_hotspot, liana_control)
  
  # ★★★ 修复了这里的管道操作逻辑 ★★★
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
  
  # --- 5. 聚焦关键细胞 (Focus Cells) ---
  pair_deltas_focus <- pair_deltas_all %>%
    dplyr::filter(
      source %in% focus_cells,
      target %in% focus_cells,
      up_in == "Enriched in Hotspot"
    )
  
  if (nrow(pair_deltas_focus) == 0) {
    warning("No enriched pairs found among focus_cells.")
    return(NULL)
  }
  
  # 添加 Axis 标签
  pair_deltas_focus <- pair_deltas_focus %>%
    dplyr::mutate(axis = paste0(ligand.complex, "-", receptor.complex))
  
  # --- 6. 统计汇总 ---
  # Axis Summary
  axis_summary <- pair_deltas_focus %>%
    dplyr::mutate(delta_pos = pmax(delta, 0)) %>%
    dplyr::group_by(axis) %>%
    dplyr::summarise(
      n_pairs   = dplyr::n(),
      sum_delta = sum(delta_pos, na.rm = TRUE)
    ) %>%
    dplyr::arrange(dplyr::desc(sum_delta))
  
  # Flow Summary
  flow_axis <- pair_deltas_focus %>%
    dplyr::mutate(
      flow = paste(source, "->", target),
      delta_pos = pmax(delta, 0)
    ) %>%
    dplyr::group_by(flow, axis) %>%
    dplyr::summarise(
      sum_delta = sum(delta_pos, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(sum_delta))
  
  readr::write_csv(pair_deltas_focus, file.path(out_dir, "LIANA_Delta_FocusCells.csv"))
  readr::write_csv(axis_summary, file.path(out_dir, "LIANA_AxisSummary.csv"))
  
  message(">>> Analysis Complete. Results saved to: ", out_dir)
  
  return(list(
    pair_deltas_all = pair_deltas_all,
    pair_deltas_focus = pair_deltas_focus,
    axis_summary = axis_summary,
    flow_axis = flow_axis
  ))
}







library(Seurat)
library(Matrix)
library(dbscan)
library(ggplot2)
library(viridis)
library(scales)

#' Plot Spatial Ligand-Receptor Interaction (Nature/Science Publication Quality)
#' 
#' @description
#' 顶刊级空间互作图：
#' - 科学配色：深空背景 + 荧光热图
#' - 多层视觉：轮廓层 + 互作层 + 高亮层
#' - 专业标注：物理比例尺 + 渐变图例 + 统计标注
#' - 出版级分辨率：矢量化 + 高 DPI
#' 
Plot_Spatial_LR_Nature <- function(
    object, 
    ligand, 
    receptor, 
    radius_um = 50,
    
    # --- 视觉风格 ---
    style = c("dark_fire", "dark_ice", "white_minimal", "nature_classic"),
    show_tissue_outline = TRUE,      # 显示组织轮廓
    outline_color = "grey25",
    outline_alpha = 0.3,
    
    # --- 互作热图控制 ---
    min_threshold = 0.01,            # 过滤噪音
    max_percentile = 99,             # 颜色上限(避免极值压缩动态范围)
    use_log_scale = FALSE,           # 对数变换(适合跨度大的数据)
    
    # --- 高亮热点 ---
    highlight_top = 0.05,            # 高亮前5%的高互作区域
    highlight_color = "white",
    highlight_size = 1.2,
    
    # --- 统计标注 ---
    add_statistics = TRUE,           # 显示互作强度分布统计
    stat_position = "topleft",
    
    # --- 物理比例尺 ---
    add_scale_bar = TRUE,
    scale_bar_length_mm = 0.5,
    scale_bar_position = "bottomright",
    scale_bar_color = NULL,          # NULL = 自动根据风格选择
    
    # --- 图例 ---
    legend_title = "Interaction\nPotential",
    legend_position = "right",
    
    # --- 输出 ---
    save_path = NULL,
    width = 10,
    height = 10,
    dpi = 600,                       # 出版级 DPI
    
    ...
) {
  
  style <- match.arg(style)
  
  # ============================================================================
  # 1. 数据准备与平滑计算
  # ============================================================================
  if (!all(c(ligand, receptor) %in% rownames(object))) {
    missing <- setdiff(c(ligand, receptor), rownames(object))
    stop(paste("Gene(s) not found:", paste(missing, collapse=", ")))
  }
  
  message(sprintf(">>> Nature-Quality Plot: %s × %s (Radius: %d µm)", 
                  ligand, receptor, radius_um))
  
  # 比例尺
  mpp <- object@misc$microns_per_pixel %||% 
    object@misc$scales$microns_per_pixel %||% 1
  if ("microns_per_pixel" %in% colnames(object@meta.data)) {
    mpp <- object@meta.data$microns_per_pixel[1]
  }
  mpp <- as.numeric(mpp)
  radius_px <- radius_um / mpp
  
  # 坐标
  coords <- cbind(object@meta.data$imagecol, object@meta.data$imagerow)
  
  # 表达量
  DefaultAssay(object) <- "Spatial"
  mat <- Seurat::GetAssayData(object, layer = "data")
  expr_L <- as.numeric(mat[ligand, ])
  expr_R <- as.numeric(mat[receptor, ])
  
  # 空间平滑
  message("   Spatial smoothing...")
  nn <- dbscan::frNN(coords, eps = radius_px)
  ids <- nn$id
  lens <- lengths(ids)
  
  total_cells <- ncol(object)
  i_idx <- c(rep(which(lens > 0), times = lens[lens > 0]), 1:total_cells)
  j_idx <- c(unlist(ids[lens > 0]), 1:total_cells)
  
  W <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = 1, 
                            dims = c(total_cells, total_cells))
  W <- W / pmax(Matrix::rowSums(W), 1)
  
  L_smooth <- as.numeric(W %*% expr_L)
  R_smooth <- as.numeric(W %*% expr_R)
  lr_score <- L_smooth * R_smooth
  
  # 对数变换(可选)
  if (use_log_scale) {
    lr_score <- log1p(lr_score * 1000) # log(1 + x*1000) 保留动态范围
  }
  
  # ============================================================================
  # 2. 配色方案 (Nature/Science 风格)
  # ============================================================================
  color_schemes <- list(
    dark_fire = list(
      bg = "#000000",
      outline = "grey20",
      gradient = c("#000000", "#1a0033", "#4d0066", "#cc0033", 
                   "#ff6600", "#ffcc00", "#ffffff"),
      text = "white",
      scalebar = "white"
    ),
    dark_ice = list(
      bg = "#000814",
      outline = "grey25",
      gradient = c("#000814", "#001d3d", "#003566", "#0077b6", 
                   "#00b4d8", "#90e0ef", "#ffffff"),
      text = "white",
      scalebar = "cyan"
    ),
    white_minimal = list(
      bg = "white",
      outline = "grey85",
      gradient = c("grey95", "grey70", "#fee0d2", "#fc9272", 
                   "#de2d26", "#a50f15", "#67000d"),
      text = "black",
      scalebar = "black"
    ),
    nature_classic = list(
      bg = "#f8f9fa",
      outline = "grey80",
      gradient = c("#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", 
                   "#1b7837", "#00441b"),
      text = "#2c3e50",
      scalebar = "#2c3e50"
    )
  )
  
  scheme <- color_schemes[[style]]
  if (is.null(scale_bar_color)) scale_bar_color <- scheme$scalebar
  
  # ============================================================================
  # 3. 数据层级构建
  # ============================================================================
  
  # 动态范围调整
  score_max <- quantile(lr_score, max_percentile/100, na.rm = TRUE)
  lr_score_clipped <- pmin(lr_score, score_max)
  
  # 识别高互作区域
  threshold_high <- quantile(lr_score, 1 - highlight_top, na.rm = TRUE)
  
  plot_df <- data.frame(
    x = coords[,1],
    y = -coords[,2],  # 反转Y轴
    score = lr_score_clipped,
    is_high = lr_score > threshold_high,
    has_signal = lr_score > min_threshold
  )
  
  # 分层数据
  df_outline <- plot_df[!plot_df$has_signal, ]  # 组织轮廓层
  df_signal <- plot_df[plot_df$has_signal, ]    # 互作信号层
  df_hotspot <- plot_df[plot_df$is_high, ]      # 热点高亮层
  
  # 按分数排序(让高值画在上面)
  df_signal <- df_signal[order(df_signal$score), ]
  
  # ============================================================================
  # 4. ggplot 绘图 (多层叠加)
  # ============================================================================
  
  p <- ggplot()
  
  # Layer 1: 组织轮廓 (背景层)
  if (show_tissue_outline && nrow(df_outline) > 0) {
    p <- p + geom_point(
      data = df_outline,
      aes(x = x, y = y),
      color = scheme$outline,
      alpha = outline_alpha,
      size = 0.5,
      stroke = 0
    )
  }
  
  # Layer 2: 互作热图 (主层)
  p <- p + geom_point(
    data = df_signal,
    aes(x = x, y = y, color = score),
    size = 0.8,
    stroke = 0
  ) +
    scale_color_gradientn(
      colors = scheme$gradient,
      name = legend_title,
      limits = c(0, score_max),
      oob = scales::squish,  # 处理超出范围的值
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 1,
        barheight = 10,
        frame.colour = scheme$text,
        ticks.colour = scheme$text
      )
    )
  
  # Layer 3: 热点高亮 (顶层)
  if (highlight_top > 0 && nrow(df_hotspot) > 0) {
    p <- p + geom_point(
      data = df_hotspot,
      aes(x = x, y = y),
      color = highlight_color,
      size = highlight_size,
      alpha = 0.6,
      stroke = 0
    )
  }
  
  # ============================================================================
  # 5. 物理比例尺 (矢量化)
  # ============================================================================
  if (add_scale_bar) {
    x_range <- range(plot_df$x, na.rm = TRUE)
    y_range <- range(plot_df$y, na.rm = TRUE)
    
    bar_um <- scale_bar_length_mm * 1000
    bar_px <- bar_um / mpp
    
    # 位置计算
    margin <- 0.05
    if (scale_bar_position == "bottomright") {
      x0 <- x_range[2] - margin * diff(x_range) - bar_px
      y0 <- y_range[1] + margin * diff(y_range)
    } else if (scale_bar_position == "bottomleft") {
      x0 <- x_range[1] + margin * diff(x_range)
      y0 <- y_range[1] + margin * diff(y_range)
    } else {
      x0 <- x_range[1] + margin * diff(x_range)
      y0 <- y_range[2] - margin * diff(y_range) - 50
    }
    
    # 分段比例尺
    n_segments <- 4
    seg_width <- bar_px / n_segments
    
    for (i in 1:n_segments) {
      seg_color <- if (i %% 2 == 1) scale_bar_color else scheme$bg
      p <- p + annotate("rect",
                        xmin = x0 + (i-1)*seg_width,
                        xmax = x0 + i*seg_width,
                        ymin = y0,
                        ymax = y0 + 40,
                        fill = seg_color,
                        color = scale_bar_color,
                        linewidth = 0.5)
    }
    
    # 标签
    p <- p + annotate("text",
                      x = x0 + bar_px/2,
                      y = y0 + 70,
                      label = sprintf("%g mm", scale_bar_length_mm),
                      color = scheme$text,
                      size = 3.5,
                      fontface = "bold")
  }
  
  # ============================================================================
  # 6. 统计标注
  # ============================================================================
  if (add_statistics) {
    stats_text <- sprintf(
      "n = %s cells\nMedian = %.3f\nMax = %.3f",
      format(sum(plot_df$has_signal), big.mark = ","),
      median(df_signal$score, na.rm = TRUE),
      max(df_signal$score, na.rm = TRUE)
    )
    
    x_range <- range(plot_df$x, na.rm = TRUE)
    y_range <- range(plot_df$y, na.rm = TRUE)
    
    stat_x <- x_range[1] + 0.05 * diff(x_range)
    stat_y <- y_range[2] - 0.05 * diff(y_range)
    
    p <- p + annotate("label",
                      x = stat_x,
                      y = stat_y,
                      label = stats_text,
                      hjust = 0,
                      vjust = 1,
                      size = 3,
                      color = scheme$text,
                      fill = alpha(scheme$bg, 0.8),
                      label.size = NA,
                      fontface = "italic")
  }
  
  # ============================================================================
  # 7. 主题与标注
  # ============================================================================
  p <- p +
    coord_fixed() +
    labs(
      title = bquote(bold(.(ligand)) %<->% bold(.(receptor))),
      subtitle = sprintf("Spatially-resolved interaction | Radius: %d µm", radius_um)
    ) +
    theme_void() +
    theme(
      # 背景
      plot.background = element_rect(fill = scheme$bg, color = NA),
      panel.background = element_rect(fill = scheme$bg, color = NA),
      
      # 标题
      plot.title = element_text(
        color = scheme$text,
        size = 18,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        color = scheme$text,
        size = 11,
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      
      # 图例
      legend.position = legend_position,
      legend.background = element_rect(fill = alpha(scheme$bg, 0.9), color = NA),
      legend.text = element_text(color = scheme$text, size = 9),
      legend.title = element_text(color = scheme$text, size = 11, face = "bold"),
      
      # 边距
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # ============================================================================
  # 8. 保存 (出版级质量)
  # ============================================================================
  if (!is.null(save_path)) {
    # 矢量格式(首选)
    pdf_path <- sub("\\.png$", ".pdf", save_path)
    ggsave(pdf_path, p, width = width, height = height, device = cairo_pdf)
    message(sprintf("✓ Saved vector: %s", pdf_path))
    
    # 高分辨率位图(备用)
    ggsave(save_path, p, width = width, height = height, dpi = dpi, 
           bg = scheme$bg)
    message(sprintf("✓ Saved raster: %s (%d DPI)", save_path, dpi))
  }
  
  return(p)
}





# Null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b






library(Seurat)
library(Matrix)
library(dbscan)
library(ggplot2)
library(viridis)

# small helper: null-coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b

# small helper: 0–1 scaling with safety
.scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[2] - rng[1]) || (rng[2] - rng[1]) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}





#' Plot Spatial Ligand-Receptor Interaction (Dark, Nature-style)
#'
#' Compute spatially smoothed LR score and plot with PlotExpressionV3.
library(Seurat)
library(Matrix)
library(dbscan)
library(ggplot2)
library(viridis)

Plot_Spatial_LR <- function(
    object,
    ligand,
    receptor,
    radius_um = 50,
    outline_floor = 0.15,   # 稍微提高一点底噪值，确保轮廓更实
    save_path = NULL,
    ...
) {
  # --- 1. 获取数据 ---
  # 兼容 V5 的 JoinLayers 逻辑，防止数据拿不到
  if ("Spatial" %in% names(object@assays)) {
    DefaultAssay(object) <- "Spatial"
    # 尝试获取 data 层，如果为空则尝试 counts
    mat <- tryCatch(
      Seurat::GetAssayData(object, assay = "Spatial", layer = "data"),
      error = function(e) Seurat::GetAssayData(object, assay = "Spatial", slot = "data")
    )
  } else {
    stop("Assay 'Spatial' not found.")
  }
  
  if (!all(c(ligand, receptor) %in% rownames(mat))) {
    missing <- setdiff(c(ligand, receptor), rownames(mat))
    stop("Gene(s) not found in Spatial data: ", paste(missing, collapse = ", "))
  }
  
  expr_L <- as.numeric(mat[ligand, ])
  expr_R <- as.numeric(mat[receptor, ])
  
  # --- 2. 坐标 & 物理尺度 ---
  md <- object@meta.data
  if (all(c("imagecol_scaled","imagerow_scaled") %in% colnames(md))) {
    coords <- cbind(md$imagecol_scaled, md$imagerow_scaled)
    coord_source <- "scaled"
  } else if (all(c("imagecol","imagerow") %in% colnames(md))) {
    coords <- cbind(md$imagecol, md$imagerow)
    coord_source <- "raw"
  } else {
    stop("Cannot find spatial coordinates: need imagecol(_scaled)/imagerow(_scaled).")
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  .scale01 <- function(x){
    rng <- range(x, na.rm = TRUE)
    if (!is.finite(rng[2] - rng[1]) || (rng[2] - rng[1]) == 0) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }
  
  # 尝试获取物理分辨率
  mpp <- object@misc$microns_per_pixel %||%
    object@misc$spatial_scales$microns_per_pixel %||%
    object@misc$scales$microns_per_pixel %||%
    tryCatch(Images(object)[[1]]@scale.factors$microns_per_pixel, error = function(e) NULL)
  
  s_low <- object@misc$spatial_scales$tissue_lowres_scalef %||%
    object@misc$scales$tissue_lowres_scalef %||%
    tryCatch(Images(object)[[1]]@scale.factors$tissue_lowres_scalef, error = function(e) NULL)
  
  if (is.null(mpp) || !is.finite(as.numeric(mpp))) {
    warning("Scale factor not found; assuming 1 µm per pixel.")
    mpp <- 1
  }
  mpp <- as.numeric(mpp)
  
  # 计算像素半径
  if (identical(coord_source, "scaled") && !is.null(s_low) && is.finite(as.numeric(s_low))) {
    radius_px <- radius_um * as.numeric(s_low) / mpp
  } else {
    radius_px <- radius_um / mpp
  }
  message(sprintf(">>> Calculating LR Score: %s-%s (Radius: %.1f µm ≈ %.1f px)", ligand, receptor, radius_um, radius_px))
  
  # --- 3. 空间平滑 (矩阵加速版) ---
  nn <- dbscan::frNN(coords, eps = radius_px)
  ids  <- nn$id
  lens <- lengths(ids)
  valid <- lens > 0
  if (!any(valid)) stop("No neighbors found; increase radius_um.")
  
  i_idx <- rep(which(valid), times = lens[valid])
  j_idx <- unlist(ids[valid])
  
  total_cells <- ncol(object)
  # 加入自环 (Self-loop)
  i_idx <- c(i_idx, seq_len(total_cells))
  j_idx <- c(j_idx, seq_len(total_cells))
  
  W <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = 1, dims = c(total_cells, total_cells))
  rs <- Matrix::rowSums(W); rs[rs == 0] <- 1
  W <- W / rs
  
  # 计算平滑表达量
  L_smooth <- as.numeric(W %*% expr_L)
  R_smooth <- as.numeric(W %*% expr_R)
  
  # 归一化后相乘
  L_s <- .scale01(L_smooth)
  R_s <- .scale01(R_smooth)
  
  score_raw <- L_s * R_s
  score_raw[!is.finite(score_raw)] <- 0
  score_pos <- pmax(score_raw, 0)
  
  # --- 4. 高值截断 (Top 0.5%) ---
  q_hi <- stats::quantile(score_pos[score_pos > 0], probs = 0.995, na.rm = TRUE)
  if (!is.finite(q_hi) || q_hi <= 0) {
    score_scaled <- score_pos
  } else {
    score_scaled <- pmin(score_pos, q_hi) / q_hi
  }
  score_scaled[!is.finite(score_scaled)] <- 0
  
  # --- 5. 视觉增强处理 ---
  # outline_floor 决定了非热点区域的“底色亮度”
  outline_floor <- max(0, min(0.4, outline_floor)) 
  
  score_vis <- ifelse(
    score_scaled > 0.01, # 只有真正有分数的点才参与渐变
    outline_floor + (1 - outline_floor) * score_scaled,
    outline_floor # 背景点设为地板值
  )
  
  feature_name_vis <- paste0("LR_", ligand, "_", receptor, "_vis")
  object[[feature_name_vis]] <- score_vis
  
  # --- 6. 调色板修改 (关键步骤) ---
  # 第一个颜色改为 grey25，这样在纯黑背景下能看到深灰色的组织轮廓
  dark_palette <- c(
    "grey25",           # <--- 背景轮廓色 (可见的深灰)
    viridis::magma(50)  # <--- 热点渐变色 (紫-红-黄-白)
  )
  
  # --- 7. 绘图 ---
  # 假设 PlotExpressionV3 是您环境中的辅助函数
  p <- PlotExpressionV3(
    object      = object,
    feature     = feature_name_vis,
    source      = "meta",
    shape       = "circle",
    spot_size   = "constant",
    colors      = dark_palette,
    legend_title = paste0(ligand, " × ", receptor, "\nLR hotspot score"),
    ...
  )
  
  # 黑色主题美化
  p <- p +
    theme(
      plot.background   = element_rect(fill = "black", color = "black"),
      panel.background  = element_rect(fill = "black", color = "black"),
      legend.background = element_rect(fill = "black", color = NA),
      legend.key        = element_rect(fill = "black", color = NA),
      legend.text       = element_text(color = "white", size = 10),
      legend.title      = element_text(color = "white", face = "bold", size = 11),
      plot.title        = element_text(color = "white", size = 16, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(color = "grey80", size = 12, hjust = 0.5),
      axis.text         = element_blank(),
      axis.title        = element_blank(),
      axis.ticks        = element_blank(),
      panel.grid        = element_blank()
    ) +
    ggtitle(
      paste0(ligand, " – ", receptor),
      subtitle = "Smoothed Interaction (Dark Mode)"
    )
  
  # 保存
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 6, height = 6, dpi = 600)
    message("Plot saved to: ", save_path)
  }
  
  return(p)
}




library(liana)
library(dplyr)
library(tidyr)
library(readr)
library(Seurat)

RunLR_Merged_Hotspot_vs_Control_v2 <- function(
    object,
    group_col        = "ECM_Hotspot",   # 分组列（Hotspot vs Control）
    group_hotspot    = TRUE,           # Hotspot 标记
    group_control    = FALSE,          # Control 标记
    label_col        = "final_label",  # 细胞类型列
    focus_cells      = NULL,           # 关注的细胞类型向量（可为全谱系）
    min_cells        = 20,             # 每组中每个 cell type 的最小细胞数
    delta_thresh     = 0.1,            # 判定“富集”的 Δ 阈值
    out_dir          = NULL,
    liana_resource   = "OmniPath",
    liana_methods    = c("natmi","connectome","logfc","sca","cellphonedb"),
    assay            = NULL,           # 可选：指定使用哪个 assay（默认用当前 DefaultAssay）
    slot             = "data"          # 可选：传给 LIANA 使用的表达层
){
  message("========== RunLR_Merged_Hotspot_vs_Control_v2 (ultimate) ==========")
  
  # ------------------------------------------------------------------
  # 0. 基本检查 + 准备工作
  # ------------------------------------------------------------------
  if (!label_col %in% colnames(object@meta.data)) {
    stop("label_col '", label_col, "' 不在 meta.data 里，请检查。")
  }
  if (!group_col %in% colnames(object@meta.data)) {
    stop("group_col '", group_col, "' 不在 meta.data 里，请检查。")
  }
  
  # 指定 assay（如果提供）
  if (!is.null(assay)) {
    if (!assay %in% Assays(object)) {
      stop("指定的 assay '", assay, "' 不存在于该 Seurat 对象中。")
    }
    DefaultAssay(object) <- assay
  }
  
  # Seurat v5: 如果有 split layers（e.g. counts.data），先合并
  if (any(grepl("\\.", Layers(object)))) {
    message("   Detected split layers in Seurat v5 -> JoinLayers()...")
    object <- JoinLayers(object)
  }
  
  # 去掉 group 或 label 为 NA 的细胞（保险）
  md  <- object@meta.data
  keep_cells <- rownames(md)[
    !is.na(md[[group_col]]) &
      !is.na(md[[label_col]])
  ]
  object <- subset(object, cells = keep_cells)
  md  <- object@meta.data  # 更新
  
  # 设置 Idents = 细胞类型
  Idents(object) <- md[[label_col]]
  
  # 如果用户提供 focus_cells，则先看哪些真的存在
  if (!is.null(focus_cells)) {
    ct_all   <- sort(unique(as.character(Idents(object))))
    ct_exist <- intersect(focus_cells, ct_all)
    ct_miss  <- setdiff(focus_cells, ct_all)
    
    if (length(ct_exist) == 0) {
      stop("focus_cells 中的细胞类型在对象中一个都找不到，请检查拼写。")
    }
    
    if (length(ct_miss) > 0) {
      message("   [Warning] 以下 focus_cells 在对象中不存在，将被忽略：")
      message("   ", paste(ct_miss, collapse = ", "))
    }
    
    focus_cells_use <- ct_exist
  } else {
    # 如果没给 focus_cells，就默认用所有细胞类型
    focus_cells_use <- sort(unique(as.character(Idents(object))))
  }
  
  # 输出路径
  if (is.null(out_dir)) {
    base_dir <- if (exists("Config") && !is.null(Config$ECMFront)) Config$ECMFront else "."
    out_dir  <- file.path(base_dir, "LIANA_Merged_Analysis")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  message("   Results will be written to: ", out_dir)
  
  message(">>> Comparison design: ", group_col, "  (Hotspot=", group_hotspot,
          "  vs  Control=", group_control, ")")
  
  # ------------------------------------------------------------------
  # 1. 按组提取细胞
  # ------------------------------------------------------------------
  cells_hotspot <- rownames(md)[md[[group_col]] == group_hotspot]
  cells_control <- rownames(md)[md[[group_col]] == group_control]
  
  cells_hotspot <- na.omit(cells_hotspot)
  cells_control <- na.omit(cells_control)
  
  message("   Hotspot cells: ", length(cells_hotspot))
  message("   Control  cells: ", length(cells_control))
  
  if (length(cells_hotspot) == 0 || length(cells_control) == 0) {
    stop("Error: Hotspot 或 Control 组细胞数为 0，请检查 ", group_col, " 的取值。")
  }
  
  # 辅助函数：过滤掉在该组里细胞数 < min_cells 的 cell type
  filter_rare_celltypes <- function(seu, min_n, label_col, focus_cells_use) {
    if (is.null(seu) || ncol(seu) == 0) return(NULL)
    
    Idents(seu) <- seu@meta.data[[label_col]]
    # 如果只关心 focus_cells，用它们来过滤
    seu <- subset(seu, idents = focus_cells_use)
    if (ncol(seu) == 0) return(NULL)
    
    ct_counts <- table(Idents(seu))
    keep_ct   <- names(ct_counts[ct_counts >= min_n])
    if (length(keep_ct) == 0) return(NULL)
    
    message("      Cell types kept (>= ", min_n, " cells):")
    message("      ", paste(paste0(keep_ct, " (n=", ct_counts[keep_ct], ")"), collapse = "; "))
    
    seu_sub <- subset(seu, idents = keep_ct)
    return(seu_sub)
  }
  
  # Hotspot 对象
  message(">>> Building Hotspot object...")
  obj_hotspot <- subset(object, cells = cells_hotspot)
  obj_hotspot <- filter_rare_celltypes(
    seu = obj_hotspot, 
    min_n = min_cells,
    label_col = label_col,
    focus_cells_use = focus_cells_use
  )
  
  # Control 对象
  message(">>> Building Control object...")
  obj_control <- subset(object, cells = cells_control)
  obj_control <- filter_rare_celltypes(
    seu = obj_control, 
    min_n = min_cells,
    label_col = label_col,
    focus_cells_use = focus_cells_use
  )
  
  if (is.null(obj_hotspot) || is.null(obj_control)) {
    stop("Error: 在按 min_cells 过滤后，Hotspot 或 Control 为空，请适当调低 min_cells 或检查细胞数量。")
  }
  
  # ------------------------------------------------------------------
  # 2. 跑 LIANA（按组）
  # ------------------------------------------------------------------
  run_liana_one <- function(seu, liana_resource, liana_methods, slot){
    message("   Running LIANA on ", ncol(seu), " cells ...")
    # 确保 DefaultAssay 和 slot 与整体一致
    DefaultAssay(seu) <- DefaultAssay(object)
    res <- liana_wrap(
      seu,
      resource = liana_resource,
      method   = liana_methods,
      assay    = DefaultAssay(seu),
      slot     = slot
    )
    agg <- liana_aggregate(res)
    # 官方文档中 aggregate_rank 是 RRA-p 类似的值，可以看作 p-value:contentReference[oaicite:1]{index=1}
    # 这里取 score = -log10(aggregate_rank)，方便做差（越大越显著）
    agg$score <- -log10(pmax(agg$aggregate_rank, 1e-10))
    return(agg)
  }
  
  message(">>> Running LIANA: Hotspot group")
  liana_hotspot <- run_liana_one(obj_hotspot, liana_resource, liana_methods, slot)
  liana_hotspot$region <- "Hotspot"
  
  message(">>> Running LIANA: Control group")
  liana_control <- run_liana_one(obj_control, liana_resource, liana_methods, slot)
  liana_control$region <- "Control"
  
  # ------------------------------------------------------------------
  # 3. 计算 Delta（Hotspot vs Control）
  # ------------------------------------------------------------------
  message(">>> Calculating Δscore (Hotspot - Control) ...")
  liana_combined <- bind_rows(liana_hotspot, liana_control)
  
  pair_deltas_all <- liana_combined %>%
    dplyr::select(
      source, target,
      ligand.complex, receptor.complex,
      region, score
    ) %>%
    tidyr::pivot_wider(
      names_from  = region,
      values_from = score,
      values_fill = 0
    ) %>%
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
  
  # ------------------------------------------------------------------
  # 4. 聚焦 focus_cells（全谱系或子集）
  # ------------------------------------------------------------------
  pair_deltas_focus <- pair_deltas_all %>%
    dplyr::filter(
      source %in% focus_cells_use,
      target %in% focus_cells_use,
      up_in == "Enriched in Hotspot"
    )
  
  if (nrow(pair_deltas_focus) == 0) {
    warning("在 focus_cells 中没有找到富集于 Hotspot 的互作。")
    return(list(
      pair_deltas_all   = pair_deltas_all,
      pair_deltas_focus = NULL,
      axis_summary      = NULL,
      flow_axis         = NULL
    ))
  }
  
  pair_deltas_focus <- pair_deltas_focus %>%
    dplyr::mutate(axis = paste0(ligand.complex, "-", receptor.complex))
  
  # ------------------------------------------------------------------
  # 5. 统计 Axis / Flow 汇总
  # ------------------------------------------------------------------
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
  
  message(">>> LIANA contrast analysis finished.")
  
  return(list(
    pair_deltas_all   = pair_deltas_all,
    pair_deltas_focus = pair_deltas_focus,
    axis_summary      = axis_summary,
    flow_axis         = flow_axis
  ))
}




library(Seurat)
library(Matrix)
library(dplyr)
library(dbscan)
library(ComplexHeatmap)
library(circlize)

Identify_HD_Niche <- function(
    object,
    celltype_col,
    feature_col = "Score_ECM",
    radius_um = 50,
    n_niche_clusters = 10,
    prob_quantile = 0.60,
    weight_by_abundance = TRUE,
    exclude_clusters = c("LowQC/Ambiguous/Doublet", "Unassigned", "Unknown", "Doublet/Ambiguous"),
    required_celltypes = c("myCAF", "SPP1+ TAM", "Endothelial"),
    top_n_niches = 2,
    hotspot_core_celltypes = c("myCAF", "SPP1+ TAM", "FOLR2+ TAM", "Endothelial", "Perivascular/SMC")
){
  # small helper
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  message(">>> [Phase 1] Initialization (Final + core)...")
  md <- object@meta.data
  
  if (!all(c("imagecol", "imagerow") %in% colnames(md))) {
    stop("Missing spatial coordinates 'imagecol' / 'imagerow' in meta.data.")
  }
  if (!celltype_col %in% colnames(md)) {
    stop("celltype_col not found in meta.data: ", celltype_col)
  }
  if (!feature_col %in% colnames(md)) {
    stop("feature_col not found in meta.data: ", feature_col)
  }
  
  # ------------------------------------------------------------------
  # 1. Coordinates and scale (已修改：增加对 misc$scales 的检查)
  # ------------------------------------------------------------------
  coords <- cbind(md$imagecol, md$imagerow)
  rownames(coords) <- rownames(md)
  
  # 修改了下面这一段赋值逻辑：
  mpp <- object@misc$microns_per_pixel %||% 
    object@misc$scales$microns_per_pixel %||%  # <--- 新增：兼容 HDA1 结构
    (md$microns_per_pixel[1])
  
  if (is.null(mpp) || is.na(mpp)) {
    warning("microns_per_pixel not found; assuming mpp = 1.")
    mpp <- 1
  }
  radius_px <- radius_um / mpp
  message(sprintf("   Radius: %d µm -> %.2f pixels (mpp = %.4f)", radius_um, radius_px, mpp))
  
  # ------------------------------------------------------------------
  # 2. Neighborhood composition (Result A logic)
  # ------------------------------------------------------------------
  message(">>> [Phase 2] Building spatial composition features...")
  
  # Label handling: keep all rows, but re-label excluded clusters
  raw_labels <- as.character(md[[celltype_col]])
  raw_labels[is.na(raw_labels)] <- "Unassigned"
  
  if (!is.null(exclude_clusters)) {
    mask_excl <- raw_labels %in% exclude_clusters
    raw_labels[mask_excl] <- "Excluded"
  }
  
  ct_factor <- as.factor(raw_labels)
  ct_mat <- model.matrix(~ ct_factor - 1)
  colnames(ct_mat) <- levels(ct_factor)
  
  # Remove "Excluded" column from feature space
  if ("Excluded" %in% colnames(ct_mat)) {
    ct_mat <- ct_mat[, colnames(ct_mat) != "Excluded", drop = FALSE]
  }
  
  # Neighborhood search
  nn_res <- dbscan::frNN(coords, eps = radius_px)
  ids   <- nn_res$id
  lens  <- lengths(ids)
  valid_cells <- lens > 0
  
  if (sum(valid_cells) < 100) {
    stop("Radius too small: very few cells have neighbors. Increase radius_um.")
  }
  
  # Sparse row-normalized kernel W
  i_idx  <- rep(which(valid_cells), times = lens[valid_cells])
  j_idx  <- unlist(ids[valid_cells])
  w_vals <- 1 / rep(lens[valid_cells], times = lens[valid_cells])
  
  W <- Matrix::sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = w_vals,
    dims = c(nrow(coords), nrow(coords))
  )
  
  # Neighborhood composition
  niche_mat <- as.matrix(W %*% ct_mat)
  rownames(niche_mat) <- rownames(coords)
  
  # Abundance weighting
  if (weight_by_abundance) {
    message("   Applying abundance weighting (inverse global proportion)...")
    global_props <- colMeans(ct_mat)
    epsilon      <- 1e-4
    weights_vec  <- 1 / (global_props + epsilon)
    weights_vec  <- pmin(weights_vec, 100)
    niche_mat_weighted <- t(t(niche_mat) * weights_vec)
    niche_mat_for_clustering <- log1p(niche_mat_weighted)
  } else {
    niche_mat_for_clustering <- niche_mat
  }
  
  # ECM smoothing
  ecm_score <- md[[feature_col]]
  ecm_score[is.na(ecm_score)] <- 0
  ecm_smooth <- as.numeric(W %*% matrix(ecm_score, ncol = 1))
  
  # ------------------------------------------------------------------
  # 3. K-means clustering into niches
  # ------------------------------------------------------------------
  message(sprintf(">>> [Phase 3] K-means clustering into %d niches...", n_niche_clusters))
  set.seed(123)
  
  niche_features <- niche_mat_for_clustering[valid_cells, , drop = FALSE]
  km_res <- kmeans(
    x        = niche_features,
    centers  = n_niche_clusters,
    iter.max = 50,
    nstart   = 5
  )
  
  niche_ids <- rep(NA_character_, nrow(coords))
  niche_ids[valid_cells] <- paste0("CN_", km_res$cluster)
  
  object$Niche_ID         <- niche_ids
  object$Score_ECM_Smooth <- ecm_smooth
  
  # Per-niche stats
  niche_stat <- data.frame(
    Niche = niche_ids[valid_cells],
    ECM   = ecm_smooth[valid_cells],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(Niche) %>%
    dplyr::summarise(
      Mean_ECM = mean(ECM, na.rm = TRUE),
      N_cells  = dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(Mean_ECM))
  
  # ------------------------------------------------------------------
  # 4. Select ECM-rich target niches
  # ------------------------------------------------------------------
  message(">>> [Phase 4] Selecting ECM-rich target niches...")
  
  centers <- km_res$centers
  rownames(centers) <- paste0("CN_", 1:n_niche_clusters)
  colnames(centers) <- colnames(ct_mat)
  
  target_niches <- character(0)
  
  if (!is.null(required_celltypes)) {
    valid_req <- intersect(required_celltypes, colnames(centers))
    if (length(valid_req) > 0) {
      message("   Required cell-type enrichment filter: ", paste(valid_req, collapse = ", "))
      z_centers <- scale(centers)
      candidate_niches <- rownames(z_centers)[
        apply(z_centers[, valid_req, drop = FALSE], 1, function(x) any(x > 0.5))
      ]
      
      if (length(candidate_niches) > 0) {
        filtered_stats <- niche_stat %>%
          dplyr::filter(Niche %in% candidate_niches)
        target_niches <- head(filtered_stats$Niche, top_n_niches)
      } else {
        warning("   No niches enriched for required_celltypes; falling back to top ECM niches.")
        target_niches <- head(niche_stat$Niche, top_n_niches)
      }
    } else {
      warning("   None of required_celltypes present in centers; using top ECM niches.")
      target_niches <- head(niche_stat$Niche, top_n_niches)
    }
  } else {
    target_niches <- head(niche_stat$Niche, top_n_niches)
  }
  
  message("   Selected Target Niches: ", paste(target_niches, collapse = ", "))
  
  # ------------------------------------------------------------------
  # 5. Define ECM hotspots
  # ------------------------------------------------------------------
  message(">>> [Phase 5] Defining ECM hotspots (belt + core)...")
  
  is_target_niche <- !is.na(object$Niche_ID) & object$Niche_ID %in% target_niches
  
  is_hotspot <- rep(FALSE, nrow(md))
  names(is_hotspot) <- rownames(md)
  thresh_val <- NA_real_
  
  if (sum(is_target_niche, na.rm = TRUE) > 0) {
    target_vals <- ecm_smooth[is_target_niche]
    target_vals <- target_vals[is.finite(target_vals)]
    
    if (length(target_vals) > 0) {
      thresh_val <- stats::quantile(
        target_vals,
        probs = prob_quantile,
        na.rm = TRUE,
        names = FALSE
      )
      message(sprintf("   Local ECM threshold (q=%.2f): %.4f", prob_quantile, thresh_val))
      is_hotspot[is_target_niche] <- ecm_smooth[is_target_niche] >= thresh_val
    } else {
      warning("   No finite ECM values in target niches; ECM_Hotspot will be all FALSE.")
    }
  } else {
    warning("   No cells in target niches; ECM_Hotspot will be all FALSE.")
  }
  
  object$ECM_Hotspot <- is_hotspot
  
  message(sprintf("   ECM_Hotspot cells: %d (%.2f%% of all cells)",
                  sum(is_hotspot),
                  100 * sum(is_hotspot) / length(is_hotspot)))
  
  # --- core ---
  core_flag <- rep(FALSE, length(is_hotspot))
  names(core_flag) <- names(is_hotspot)
  
  if (!is.null(hotspot_core_celltypes) && length(hotspot_core_celltypes) > 0) {
    lab <- as.character(md[[celltype_col]])
    core_flag <- is_hotspot & (lab %in% hotspot_core_celltypes)
  }
  
  object$ECM_Hotspot_core <- core_flag
  
  # ------------------------------------------------------------------
  # 6. Heatmap
  # ------------------------------------------------------------------
  message(">>> [Phase 6] Building niche composition heatmap...")
  
  if (!requireNamespace("viridis", quietly = TRUE)) {
    warning("Package 'viridis' not installed; using default colors in Heatmap.")
    col_fun <- circlize::colorRamp2(c(min(centers), max(centers)), c("white", "black"))
    hm <- ComplexHeatmap::Heatmap(
      t(centers),
      name = "Enrichment",
      col  = col_fun,
      column_title = "Niche Composition Profile",
      cluster_rows = TRUE,
      cluster_columns = TRUE
    )
  } else {
    hm <- ComplexHeatmap::Heatmap(
      t(centers),
      name = "Enrichment",
      col  = viridis::viridis(100),
      column_title = "Niche Composition Profile",
      cluster_rows = TRUE,
      cluster_columns = TRUE
    )
  }
  
  # ------------------------------------------------------------------
  # 7. Return
  # ------------------------------------------------------------------
  message(">>> [Done] Niche identification and ECM hotspot calling completed.")
  
  return(list(
    object                 = object,
    niche_stats            = niche_stat,
    target_niches          = target_niches,
    hotspot_threshold      = thresh_val,
    hotspot_core_celltypes = hotspot_core_celltypes,
    plot_heatmap           = hm
  ))
}








library(liana)
library(dplyr)
library(tidyr)
library(readr)
library(Seurat)

RunLR_Merged_Hotspot_vs_Control_v3 <- function(
    object,
    group_col        = "ECM_Hotspot",   # column for group labels (Hotspot vs Control)
    group_hotspot    = TRUE,            # value defining Hotspot
    group_control    = FALSE,           # value defining Control
    label_col        = "final_label",   # column with cell-type annotation
    focus_cells      = NULL,            # vector of cell types of interest (character)
    min_cells        = 20,              # minimal cells per cell type per group
    delta_thresh     = 0.1,             # threshold for "enriched" delta
    out_dir          = NULL,
    liana_resource   = "OmniPath",
    liana_methods    = c("natmi","connectome","logfc","sca","cellphonedb"),
    assay            = NULL,            # which assay to use (NULL = keep DefaultAssay)
    slot             = "data"           # which slot to use for expression
){
  message("========== RunLR_Merged_Hotspot_vs_Control_v3 (symmetric design) ==========")
  
  # ------------------------------------------------------------------
  # 0. Basic checks and setup
  # ------------------------------------------------------------------
  if (!label_col %in% colnames(object@meta.data)) {
    stop("label_col '", label_col, "' not found in meta.data.")
  }
  if (!group_col %in% colnames(object@meta.data)) {
    stop("group_col '", group_col, "' not found in meta.data.")
  }
  
  # Debug info about group_col (you can comment these out later)
  message(">>> DEBUG: Checking group_col structure...")
  message("   Class: ", class(object@meta.data[[group_col]]))
  message("   Unique values: ", paste(unique(object@meta.data[[group_col]]), collapse = ", "))
  
  # Use specified assay if provided  ---- (FIXED HERE)
  # ---- 替换 v3 里原来的 assay 检查块 ----
  if (!is.null(assay)) {
    aa <- Seurat::Assays(object)
    
    # 兼容两种情况：
    # 1）aa 是字符向量（来自 SeuratObject） -> 直接用
    # 2）aa 是带名字的 list（标准 Seurat）   -> 用 names()
    if (is.character(aa)) {
      assay_names <- aa
    } else {
      assay_names <- names(aa)
    }
    
    if (is.null(assay_names) || length(assay_names) == 0) {
      warning(
        "Could not infer assay names from Assays(object); ",
        "will still try to set DefaultAssay(object) <- '", assay, "'."
      )
      # 直接尝试设置，如果不存在 Seurat 自己会报错
      Seurat::DefaultAssay(object) <- assay
    } else {
      if (!assay %in% assay_names) {
        stop(
          "Assay '", assay, "' not found in Seurat object. Existing assays: ",
          paste(assay_names, collapse = ", ")
        )
      }
      Seurat::DefaultAssay(object) <- assay
    }
  }
  
  
  # Join Seurat v5 layers if needed
  if (any(grepl("\\.", Layers(object)))) {
    message("   Detected split layers in Seurat v5 -> JoinLayers()...")
    object <- JoinLayers(object)
  }
  
  # Remove cells with NA group or NA label
  md <- object@meta.data
  keep_cells <- rownames(md)[
    !is.na(md[[group_col]]) &
      !is.na(md[[label_col]])
  ]
  object <- subset(object, cells = keep_cells)
  md <- object@meta.data
  
  # Set Idents = cell type
  Idents(object) <- as.character(md[[label_col]])
  
  # ------------------------------------------------------------------
  # 1. Determine focus_cells_use
  # ------------------------------------------------------------------
  ct_all <- sort(unique(as.character(Idents(object))))
  
  if (!is.null(focus_cells)) {
    focus_cells <- as.character(focus_cells)
    ct_exist <- intersect(focus_cells, ct_all)
    ct_miss  <- setdiff(focus_cells, ct_all)
    
    if (length(ct_exist) == 0) {
      stop("None of the focus_cells are present in the object. Please check spelling.")
    }
    if (length(ct_miss) > 0) {
      message("   [Warning] The following focus_cells are not present and will be ignored:")
      message("   ", paste(ct_miss, collapse = ", "))
    }
    focus_cells_use <- sort(ct_exist)
  } else {
    focus_cells_use <- ct_all
  }
  
  # ------------------------------------------------------------------
  # 2. Output directory
  # ------------------------------------------------------------------
  if (is.null(out_dir)) {
    base_dir <- if (exists("Config") && !is.null(Config$ECMFront)) Config$ECMFront else "."
    out_dir  <- file.path(base_dir, "LIANA_Merged_Analysis_v3")
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  message("   Results will be written to: ", out_dir)
  
  message(">>> Comparison design: ", group_col, "  (Hotspot=", group_hotspot,
          "  vs  Control=", group_control, ")")
  
  # ------------------------------------------------------------------
  # 3. Extract cells for Hotspot and Control (robust for logical/character)
  # ------------------------------------------------------------------
  group_vals <- md[[group_col]]
  
  if (is.logical(group_vals)) {
    cells_hotspot <- rownames(md)[which(group_vals == group_hotspot)]
    cells_control <- rownames(md)[which(group_vals == group_control)]
  } else {
    cells_hotspot <- rownames(md)[which(as.character(group_vals) == as.character(group_hotspot))]
    cells_control <- rownames(md)[which(as.character(group_vals) == as.character(group_control))]
  }
  
  cells_hotspot <- cells_hotspot[!is.na(cells_hotspot)]
  cells_control <- cells_control[!is.na(cells_control)]
  
  message("   Hotspot cells: ", length(cells_hotspot))
  message("   Control  cells: ", length(cells_control))
  
  if (length(cells_hotspot) == 0 || length(cells_control) == 0) {
    stop("Error: Hotspot or Control group has 0 cells. Please check ", group_col, " values.")
  }
  
  # Build temporary Seurat objects per group
  message(">>> Building group-specific objects (before symmetry filter)...")
  obj_hotspot_tmp <- subset(object, cells = cells_hotspot)
  obj_control_tmp <- subset(object, cells = cells_control)
  
  obj_hotspot_tmp <- SetIdent(
    obj_hotspot_tmp,
    cells = colnames(obj_hotspot_tmp),
    value = as.character(obj_hotspot_tmp@meta.data[[label_col]])
  )
  obj_control_tmp <- SetIdent(
    obj_control_tmp,
    cells = colnames(obj_control_tmp),
    value = as.character(obj_control_tmp@meta.data[[label_col]])
  )
  
  # Restrict to focus_cells
  obj_hotspot_tmp <- subset(obj_hotspot_tmp, idents = focus_cells_use)
  obj_control_tmp <- subset(obj_control_tmp, idents = focus_cells_use)
  
  if (ncol(obj_hotspot_tmp) == 0 || ncol(obj_control_tmp) == 0) {
    stop("After restricting to focus_cells, one of the groups has 0 cells.")
  }
  
  # ------------------------------------------------------------------
  # 4. Symmetric cell-type design (min_cells in BOTH groups)
  # ------------------------------------------------------------------
  ct_hot_counts <- table(Idents(obj_hotspot_tmp))
  ct_ctl_counts <- table(Idents(obj_control_tmp))
  
  ct_hot_keep <- names(ct_hot_counts[ct_hot_counts >= min_cells])
  ct_ctl_keep <- names(ct_ctl_counts[ct_ctl_counts >= min_cells])
  
  common_ct <- sort(intersect(ct_hot_keep, ct_ctl_keep))
  
  if (length(common_ct) == 0) {
    stop(
      "No cell types have >= ", min_cells,
      " cells in BOTH Hotspot and Control after focus_cells filtering."
    )
  }
  
  message("   Common cell types kept (>= ", min_cells, " cells in both groups):")
  for (ct in common_ct) {
    n_hot <- ct_hot_counts[ct]
    n_ctl <- ct_ctl_counts[ct]
    message("      ", ct, ": Hotspot=", n_hot, ", Control=", n_ctl)
  }
  
  obj_hotspot <- subset(obj_hotspot_tmp, idents = common_ct)
  obj_control <- subset(obj_control_tmp, idents = common_ct)
  
  if (ncol(obj_hotspot) == 0 || ncol(obj_control) == 0) {
    stop("After applying common_ct filter, one of the group objects is empty.")
  }
  
  # ------------------------------------------------------------------
  # 5. Run LIANA in each group
  # ------------------------------------------------------------------
  run_liana_one <- function(seu, liana_resource, liana_methods, slot, group_name){
    message("   Running LIANA on ", group_name, " (", ncol(seu), " cells; ",
            length(unique(Idents(seu))), " cell types ) ...")
    DefaultAssay(seu) <- DefaultAssay(object)
    res <- liana_wrap(
      seu,
      resource = liana_resource,
      method   = liana_methods,
      assay    = DefaultAssay(seu),
      slot     = slot
    )
    agg <- liana_aggregate(res)
    agg$score <- -log10(pmax(agg$aggregate_rank, 1e-10))
    return(agg)
  }
  
  message(">>> Running LIANA: Hotspot group")
  liana_hotspot <- run_liana_one(
    seu            = obj_hotspot,
    liana_resource = liana_resource,
    liana_methods  = liana_methods,
    slot           = slot,
    group_name     = "Hotspot"
  )
  liana_hotspot$region <- "Hotspot"
  
  message(">>> Running LIANA: Control group")
  liana_control <- run_liana_one(
    seu            = obj_control,
    liana_resource = liana_resource,
    liana_methods  = liana_methods,
    slot           = slot,
    group_name     = "Control"
  )
  liana_control$region <- "Control"
  
  # ------------------------------------------------------------------
  # 6. Compute delta scores
  # ------------------------------------------------------------------
  message(">>> Calculating Δscore (Hotspot - Control) ...")
  liana_combined <- dplyr::bind_rows(liana_hotspot, liana_control)
  
  pair_deltas_all <- liana_combined %>%
    dplyr::select(
      source, target,
      ligand.complex, receptor.complex,
      region, score
    ) %>%
    tidyr::pivot_wider(
      names_from  = region,
      values_from = score,
      values_fill = 0
    ) %>%
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
  
  # ------------------------------------------------------------------
  # 7. Focus on enriched interactions among common_ct
  # ------------------------------------------------------------------
  pair_deltas_focus <- pair_deltas_all %>%
    dplyr::filter(
      source %in% common_ct,
      target %in% common_ct,
      up_in == "Enriched in Hotspot"
    )
  
  if (nrow(pair_deltas_focus) == 0) {
    warning("No interactions enriched in Hotspot among common cell types.")
    return(list(
      pair_deltas_all   = pair_deltas_all,
      pair_deltas_focus = NULL,
      axis_summary      = NULL,
      flow_axis         = NULL,
      common_celltypes  = common_ct
    ))
  }
  
  pair_deltas_focus <- pair_deltas_focus %>%
    dplyr::mutate(axis = paste0(ligand.complex, "-", receptor.complex))
  
  # ------------------------------------------------------------------
  # 8. Summaries
  # ------------------------------------------------------------------
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
  
  message(">>> LIANA symmetric contrast analysis finished.")
  
  return(list(
    pair_deltas_all   = pair_deltas_all,
    pair_deltas_focus = pair_deltas_focus,
    axis_summary      = axis_summary,
    flow_axis         = flow_axis,
    common_celltypes  = common_ct
  ))
}





#------------------------------------------------------------
# RunMultiNiche_ECMHotspot
#   Wrapper around multinichenetr to compare ECM_Hotspot+ vs -
#   across multiple samples (here: HDD1 / HDA1 / HDD2).
#------------------------------------------------------------
RunMultiNiche_ECMHotspot <- function(
    seu,                          # Seurat object (merged)
    sample_col  = "orig.ident",   # sample / patient id
    group_col   = "ECM_Hotspot",  # ECM hotspot indicator (logical or factor)
    label_col   = "final_label",  # cell type annotation
    focus_cells = NULL,           # optional vector of cell types to keep
    min_cells   = 50,             # min cells per sample × celltype
    # multinichenetr / nichenetr resources:
    lr_network            = NULL, # data.frame with ligand / receptor columns
    ligand_target_matrix  = NULL, # ligand–target matrix from NicheNet
    out_dir    = NULL,
    contrast_names = c("Hotspot", "Control"),  # names for TRUE / FALSE
    fraction_cutoff  = 0.05,   # fraction cutoff for gene filtering
    min_sample_prop  = 0.5,    # at least this fraction of samples express the gene
    logFC_threshold  = 0.5,    # DE threshold for ligand activity
    p_val_threshold  = 0.05,
    top_n_targets    = 250,
    top_n_lr_pairs   = 50
){
  # ---- 0. Load required packages ----
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    stop("Please install 'SingleCellExperiment'.")
  if (!requireNamespace("multinichenetr", quietly = TRUE))
    stop("Please install 'multinichenetr'.")
  if (!requireNamespace("nichenetr", quietly = TRUE))
    stop("Please install 'nichenetr'.")
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("Please install 'Seurat'.")
  
  library(SingleCellExperiment)
  library(multinichenetr)
  library(nichenetr)
  library(dplyr)
  library(ggplot2)
  
  message("========== RunMultiNiche_ECMHotspot ==========")
  
  # ---- 1. Basic checks on meta.data ----
  md <- seu@meta.data
  
  for (col in c(sample_col, group_col, label_col)) {
    if (!col %in% colnames(md)) {
      stop("Column '", col, "' not found in meta.data.")
    }
  }
  
  # ---- 2. Prepare group variable (Hotspot vs Control) ----
  grp_raw <- md[[group_col]]
  
  # map to character labels used by multinichenetr
  if (is.logical(grp_raw)) {
    if (length(contrast_names) != 2) {
      stop("contrast_names must have length 2 for logical group_col.")
    }
    grp_chr <- ifelse(grp_raw, contrast_names[1], contrast_names[2])
  } else {
    # if already factor/character, just coerce to character
    grp_chr <- as.character(grp_raw)
  }
  
  md[[group_col]] <- grp_chr
  
  # ---- 3. Optionally restrict to focus_cells ----
  if (!is.null(focus_cells)) {
    focus_cells <- as.character(focus_cells)
    keep_ct <- md[[label_col]] %in% focus_cells
    n_before <- nrow(md)
    n_after  <- sum(keep_ct)
    message("   Restricting to focus_cells: kept ", n_after, " / ", n_before, " cells.")
    if (n_after == 0) {
      stop("No cells left after filtering by focus_cells. Check cell type names.")
    }
    seu <- subset(seu, cells = rownames(md)[keep_ct])
    md  <- seu@meta.data
  }
  
  # ---- 4. Convert Seurat -> SingleCellExperiment ----
  message("   Converting Seurat -> SingleCellExperiment ...")
  sce <- Seurat::as.SingleCellExperiment(seu, assay = Seurat::DefaultAssay(seu))
  
  # ensure metadata columns exist in colData
  cd <- as.data.frame(colData(sce))
  
  cd[[sample_col]] <- as.character(md[[sample_col]])
  cd[[group_col]]  <- as.character(md[[group_col]])
  cd[[label_col]]  <- as.character(md[[label_col]])
  
  # use sample as batch (you have 3 samples: HDD1 / HDA1 / HDD2)
  batch_col <- sample_col
  cd[[batch_col]] <- cd[[sample_col]]
  
  colData(sce) <- S4Vectors::DataFrame(cd)
  
  # ---- 5. Prepare arguments for multinichenetr ----
  sample_id   <- sample_col
  group_id    <- group_col
  celltype_id <- label_col
  batches     <- batch_col
  
  # use all cell types present (after optional focus filter)
  senders_oi   <- sort(unique(cd[[celltype_id]]))
  receivers_oi <- sort(unique(cd[[celltype_id]]))
  
  message("   Unique cell types (senders/receivers):")
  message("   ", paste(senders_oi, collapse = ", "))
  
  # ---- 6. Load or check LR resources ----
  if (is.null(lr_network) || is.null(ligand_target_matrix)) {
    stop(
      "Please provide 'lr_network' and 'ligand_target_matrix'.\n",
      "Typical usage:\n",
      "  lr_network_all <- readRDS('lr_network_human_allInfo_30112033.rds')\n",
      "  lr_network <- lr_network_all %>% dplyr::distinct(ligand, receptor)\n",
      "  ligand_target_matrix <- readRDS('ligand_target_matrix_nsga2r_final.rds')"
    )
  }
  
  # ensure lr_network has ligand & receptor columns
  if (!all(c("ligand","receptor") %in% colnames(lr_network))) {
    stop("lr_network must contain columns 'ligand' and 'receptor'.")
  }
  
  # ---- 7. Define contrasts (Hotspot vs Control) ----
  # groups should be something like "Hotspot" and "Control"
  groups_chr <- unique(cd[[group_id]])
  message("   Groups in group_col: ", paste(groups_chr, collapse = ", "))
  
  # Expect exactly two groups
  if (length(groups_chr) != 2) {
    stop("Expected exactly 2 groups in group_col (e.g. Hotspot vs Control).")
  }
  
  # Make sure order matches contrast_names
  # contrast_names[1] = "Hotspot", [2] = "Control" by default
  g1 <- contrast_names[1]
  g2 <- contrast_names[2]
  
  if (!all(c(g1, g2) %in% groups_chr)) {
    stop("contrast_names (", g1, ", ", g2, ") not consistent with groups in group_col.")
  }
  
  # multinichenetr expects a string like "'Hotspot-Control','Control-Hotspot'"
  contrasts_oi <- paste0("'", g1, "-", g2, "','", g2, "-", g1, "'")
  
  # contrast table for prioritization
  contrast_tbl <- tibble::tibble(
    contrast = contrasts_oi,
    group    = c(g1, g2)
  )
  
  # sender–receiver grid
  sender_receiver_tbl <- tibble::tibble(
    sender   = senders_oi,
    receiver = receivers_oi
  )
  
  # grouping table (sample × group × batch)
  grouping_tbl <- cd[, c(sample_id, group_id, batches)] %>%
    dplyr::distinct()
  
  # ---- 8. Abundance info (who can send / receive?) ----
  message(">>> Step 1: get_abundance_info ...")
  abundance_info <- multinichenetr::get_abundance_info(
    sce          = sce,
    sample_id    = sample_id,
    group_id     = group_id,
    celltype_id  = celltype_id,
    min_cells    = min_cells,
    senders_oi   = senders_oi,
    receivers_oi = receivers_oi,
    batches      = batches
  )
  
  # ---- 9. Gene expression fraction (filter genes) ----
  message(">>> Step 2: get_frac_exprs ...")
  frq_list <- multinichenetr::get_frac_exprs(
    sce           = sce,
    sample_id     = sample_id,
    celltype_id   = celltype_id,
    group_id      = group_id,
    batches       = batches,
    min_cells     = min_cells,
    fraction_cutoff = fraction_cutoff,
    min_sample_prop = min_sample_prop
  )
  
  # ---- 10. Pseudobulk abundance + expression per sample × celltype ----
  message(">>> Step 3: process_abundance_expression_info ...")
  abundance_expression_info <- multinichenetr::process_abundance_expression_info(
    sce          = sce,
    sample_id    = sample_id,
    group_id     = group_id,
    celltype_id  = celltype_id,
    min_cells    = min_cells,
    senders_oi   = senders_oi,
    receivers_oi = receivers_oi,
    lr_network   = lr_network,
    batches      = batches,
    frq_list     = frq_list,
    abundance_info = abundance_info
  )
  
  # ---- 11. Differential expression (pseudobulk, muscat) ----
  message(">>> Step 4: get_DE_info (pseudobulk DE) ...")
  DE_info <- multinichenetr::get_DE_info(
    sce          = sce,
    sample_id    = sample_id,
    group_id     = group_id,
    celltype_id  = celltype_id,
    batches      = batches,
    covariates   = NA,
    contrasts_oi = contrasts_oi,
    min_cells    = min_cells,
    expressed_df = frq_list$expressed_df
  )
  
  # ---- 12. Ligand activity (NicheNet score) ----
  message(">>> Step 5: get_ligand_activities_targets_DEgenes ...")
  ligand_activities_targets_DEgenes <- multinichenetr::get_ligand_activities_targets_DEgenes(
    receiver_de          = DE_info$celltype_de$de_output_tidy,
    receivers_oi         = receivers_oi,
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold      = logFC_threshold,
    p_val_threshold      = p_val_threshold,
    p_val_adj            = TRUE,
    top_n_target         = top_n_targets,
    verbose              = TRUE,
    n.cores              = 1
  )
  
  # ---- 13. Integrated prioritization of LR pairs ----
  message(">>> Step 6: generate_prioritization_tables ...")
  prioritization_tables <- multinichenetr::generate_prioritization_tables(
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de   = multinichenetr::combine_sender_receiver_de(
      sender_de   = DE_info$celltype_de$de_output_tidy,
      receiver_de = DE_info$celltype_de$de_output_tidy,
      senders_oi  = senders_oi,
      receivers_oi= receivers_oi,
      lr_network  = lr_network
    ),
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl        = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl        = grouping_tbl,
    scenario            = "regular",
    fraction_cutoff     = fraction_cutoff,
    abundance_data_receiver = abundance_info$abundance_data_receiver,
    abundance_data_sender   = abundance_info$abundance_data_sender,
    ligand_activity_down = FALSE
  )
  
  # ---- 14. Get top LR pairs and default plots ----
  message(">>> Step 7: get_top_n_lr_pairs & make_sample_lr_prod_activity_plots ...")
  prioritized_tbl_oi <- multinichenetr::get_top_n_lr_pairs(
    prioritization_tables,
    top_n         = top_n_lr_pairs,
    rank_per_group = FALSE
  )
  
  plot_oi <- multinichenetr::make_sample_lr_prod_activity_plots(
    prioritization_tables,
    prioritized_tbl_oi
  )
  
  # ---- 15. Save results (optional) ----
  if (is.null(out_dir)) {
    out_dir <- getwd()
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  readr::write_csv(
    prioritized_tbl_oi,
    file.path(out_dir, "MultiNiche_Prioritized_LR_top.csv")
  )
  
  saveRDS(
    prioritization_tables,
    file.path(out_dir, "MultiNiche_prioritization_tables.rds")
  )
  
  message(">>> MultiNicheNet ECM_Hotspot analysis finished.")
  message("    Results written to: ", out_dir)
  
  invisible(list(
    sce                    = sce,
    abundance_info         = abundance_info,
    frq_list               = frq_list,
    DE_info                = DE_info,
    ligand_activities      = ligand_activities_targets_DEgenes,
    prioritization_tables  = prioritization_tables,
    prioritized_tbl_oi     = prioritized_tbl_oi,
    plot_oi                = plot_oi
  ))
}




RunLR_liana <- function(
    object,
    sample_col       = "orig.ident",   
    group_col        = "ECM_Hotspot",  
    target_group     = TRUE,           
    label_col        = "final_label",  
    focus_cells      = NULL,           
    min_cells        = 10,             
    liana_resource   = "OmniPath",
    liana_methods    = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
    out_dir          = NULL,
    assay            = NULL,
    # ===【新增参数】===
    mode             = "sample_wise"   # 选项: "sample_wise" (默认，分样本) 或 "merged" (合并)
){
  message("========== RunLR_liana (Mode: ", mode, ") ==========")
  
  # ---- 1. 基础设置与过滤 (通用步骤) ----
  if (!label_col %in% colnames(object@meta.data)) stop("label_col not found.")
  if (!is.null(group_col) && !group_col %in% colnames(object@meta.data)) stop("group_col not found.")
  
  if (is.null(out_dir)) {
    base_dir <- if (exists("Config") && !is.null(Config$ECMFront)) Config$ECMFront else "."
    out_dir  <- file.path(base_dir, paste0("LIANA_", mode)) # 路径随模式自动变
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  if (!is.null(assay)) DefaultAssay(object) <- assay
  if (any(grepl("\\.", Layers(object)))) object <- JoinLayers(object)
  Idents(object) <- object@meta.data[[label_col]]
  
  # 1.1 按 Group 筛选 (只保留 Hotspot)
  if (!is.null(group_col) && !is.null(target_group)) {
    message(">>> Filtering group: ", group_col, " == ", target_group)
    object <- subset(object, cells = rownames(object@meta.data)[object@meta.data[[group_col]] == target_group])
  }
  
  # 1.2 按 Focus Cells 筛选
  if (!is.null(focus_cells)) {
    valid_focus <- intersect(focus_cells, unique(as.character(Idents(object))))
    if (length(valid_focus) == 0) stop("No focus cells found.")
    object <- subset(object, idents = valid_focus)
  }
  
  # ============================================================================
  # ---- 2. 分支逻辑 (Minimal Modification) ----
  # ============================================================================
  
  if (mode == "merged") {
    # >>>>> 分支 A: 合并分析 (Merged) <<<<<
    message(">>> Running LIANA on ALL samples merged together...")
    
    # 过滤稀有细胞 (全局)
    ct_counts <- table(Idents(object))
    valid_ct  <- names(ct_counts)[ct_counts >= min_cells]
    if (length(valid_ct) < 2) stop("Less than 2 cell types remain after min_cells filter.")
    object <- subset(object, idents = valid_ct)
    
    # 运行 LIANA
    liana_res <- liana_wrap(object, resource = liana_resource, method = liana_methods, verbose = TRUE)
    liana_agg <- liana_aggregate(liana_res) %>%
      mutate(score = -log10(aggregate_rank + 1e-10)) %>%
      arrange(desc(score))
    
    # 保存与返回
    write_csv(liana_agg, file.path(out_dir, "LIANA_Merged_Results.csv"))
    message(">>> Finished. Saved to: ", out_dir)
    return(list(merged_res = liana_agg))
    
  } else {
    # >>>>> 分支 B: 分样本分析 (Sample-wise) <<<<<
    # (保留原本的循环逻辑)
    sample_ids <- unique(as.character(object@meta.data[[sample_col]]))
    message(">>> Iterating over ", length(sample_ids), " samples...")
    
    res_list <- list()
    for (sid in sample_ids) {
      cells_s <- rownames(object@meta.data)[object@meta.data[[sample_col]] == sid]
      if (length(cells_s) < 50) next 
      
      obj_s <- subset(object, cells = cells_s)
      # 样本内过滤
      ct_counts <- table(Idents(obj_s))
      valid_ct  <- names(ct_counts)[ct_counts >= min_cells]
      if (length(valid_ct) < 2) next
      
      tryCatch({
        liana_res <- liana_wrap(subset(obj_s, idents = valid_ct), resource = liana_resource, method = liana_methods, verbose = FALSE)
        agg <- liana_aggregate(liana_res)
        agg$sample_id <- sid
        res_list[[sid]] <- agg
      }, error = function(e) message("   Error in ", sid, ": ", e$message))
    }
    
    if (length(res_list) == 0) stop("No samples successful.")
    all_res <- dplyr::bind_rows(res_list)
    
    # 统计 Frequency
    summary_df <- all_res %>%
      group_by(source, target, ligand.complex, receptor.complex) %>%
      summarise(
        n_samples_detected = n_distinct(sample_id),
        mean_aggregate_rank = mean(aggregate_rank, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        frequency = n_samples_detected / length(res_list),
        consensus_score = frequency * (-log10(mean_aggregate_rank + 1e-10))
      ) %>%
      arrange(desc(consensus_score))
    
    write_csv(summary_df, file.path(out_dir, "LIANA_SampleWise_Summary.csv"))
    message(">>> Finished. Saved to: ", out_dir)
    return(list(raw_results = res_list, combined_df = all_res, summary_df = summary_df))
  }
}
