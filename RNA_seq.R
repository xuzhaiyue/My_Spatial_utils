suppressPackageStartupMessages({
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) stop("Please install ComplexHeatmap.")
  if (!requireNamespace("circlize", quietly = TRUE)) stop("Please install circlize.")
})
library(reshape2)            # 加载包
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggsci) # 顶刊配色包
read_excel_all_sheets <- function(path, id_column = "Sheet_Source") {
  # 1. 获取 Excel 中所有的 sheet 名称
  sheet_names <- readxl::excel_sheets(path)
  
  # 2. 循环读取每个 sheet，并将其存入列表
  list_of_dfs <- lapply(sheet_names, function(sheet) {
    message("正在读取 Sheet: ", sheet)
    
    # 读取当前 sheet
    df <- readxl::read_excel(path, sheet = sheet)
    
    # 在第一列添加 sheet 的名称，方便区分来源
    df[[id_column]] <- sheet
    
    return(df)
  })
  
  # 3. 将列表中的所有数据框合并成一个总的数据框
  # 使用 dplyr::bind_rows 可以自动处理列名不完全匹配的情况
  final_df <- dplyr::bind_rows(list_of_dfs)
  
  return(final_df)
}
Plot_ViolinPub <- function(
    df,                       
    mode = c("panel", "merged", "stratified"),
    sample_levels = c("HDA1","HDD1","HDD2"),
    lr_levels = NULL,         
    point_n = 3000,
    seed = 123,
    fill_colors = c("Non_ECM" = "#8DA0CB", "ECM_Hotspot" = "#FC8D62"),
    ylab = "Ligand Intensity (Smoothed, Normalized)",
    title = NULL,
    facet_scales = "free_y"
) {
  mode <- match.arg(mode)
  
  # ---- 0) Sanity Check ----
  stopifnot(is.data.frame(df))
  need <- c("sample","lr_pair","group","y")
  miss <- setdiff(need, colnames(df))
  if (length(miss) > 0) stop("df missing columns: ", paste(miss, collapse = ", "))
  
  df_plot <- df |>
    dplyr::mutate(
      group = factor(group, levels = c("Non_ECM","ECM_Hotspot")),
      sample = factor(sample, levels = sample_levels)
    ) |>
    dplyr::filter(is.finite(y))
  
  if (!is.null(lr_levels)) {
    df_plot <- df_plot |> dplyr::mutate(lr_pair = factor(lr_pair, levels = lr_levels))
  } else {
    df_plot <- df_plot |> dplyr::mutate(lr_pair = factor(lr_pair))
  }
  
  # ---- 1) Downsample points only (修复了这里的报错) ----
  set.seed(seed)
  
  if (mode == "panel") {
    df_points <- df_plot |>
      dplyr::group_by(sample, lr_pair, group) |>
      # 【修复】直接使用 n = point_n，dplyr 会自动处理行数不足的情况
      dplyr::slice_sample(n = point_n) |> 
      dplyr::ungroup()
  } else {
    df_points <- df_plot |>
      dplyr::group_by(lr_pair, group) |>
      # 【修复】直接使用 n = point_n
      dplyr::slice_sample(n = point_n) |> 
      dplyr::ungroup()
  }
  
  # ---- 2) Stats + annotation (Wilcoxon) ----
  if (mode == "panel") {
    stat_tbl <- df_plot |>
      dplyr::group_by(sample, lr_pair) |>
      dplyr::summarise(
        p = tryCatch(suppressWarnings(stats::wilcox.test(y ~ group)$p.value), error = function(e) NA_real_),
        ymax = stats::quantile(y, 0.995, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        p_adj = stats::p.adjust(p, method = "BH"),
        label = dplyr::case_when(
          p_adj < 1e-4 ~ "****",
          p_adj < 1e-3 ~ "***",
          p_adj < 1e-2 ~ "**",
          p_adj < 5e-2 ~ "*",
          TRUE ~ "ns"
        )
      )
    
    stat_anno <- stat_tbl |>
      dplyr::filter(!is.na(p), !is.na(ymax)) |>
      dplyr::mutate(
        group1 = "Non_ECM",
        group2 = "ECM_Hotspot",
        y.position = ymax + ymax * 0.1,
        p.label = label
      )
    
  } else if (mode == "merged") {
    stat_tbl <- df_plot |>
      dplyr::group_by(lr_pair) |>
      dplyr::summarise(
        p = tryCatch(suppressWarnings(stats::wilcox.test(y ~ group)$p.value), error = function(e) NA_real_),
        ymax = stats::quantile(y, 0.995, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        p_adj = stats::p.adjust(p, method = "BH"),
        label = dplyr::case_when(
          p_adj < 1e-4 ~ "****",
          p_adj < 1e-3 ~ "***",
          p_adj < 1e-2 ~ "**",
          p_adj < 5e-2 ~ "*",
          TRUE ~ "ns"
        )
      )
    
    stat_anno <- stat_tbl |>
      dplyr::filter(!is.na(p), !is.na(ymax)) |>
      dplyr::mutate(
        group1 = "Non_ECM",
        group2 = "ECM_Hotspot",
        y.position = ymax + ymax * 0.1,
        p.label = label
      )
    
  } else {
    # Stratified (coin package)
    if (!requireNamespace("coin", quietly = TRUE)) {
      stop("mode='stratified' requires package 'coin'. install.packages('coin')")
    }
    
    stat_tbl <- lapply(levels(df_plot$lr_pair), function(pair) {
      sub <- df_plot[df_plot$lr_pair == pair, , drop = FALSE]
      sub$group <- factor(sub$group, levels = c("Non_ECM","ECM_Hotspot"))
      # Coin check: if strictly one group level, cannot test
      if (length(unique(sub$group)) < 2) return(NULL)
      
      p <- tryCatch(
        coin::pvalue(coin::wilcox_test(y ~ group | sample, data = sub, distribution = "asymptotic")),
        error = function(e) NA_real_
      )
      
      if (is.na(p)) return(NULL)
      
      data.frame(
        lr_pair = pair,
        p = as.numeric(p),
        ymax = stats::quantile(sub$y, 0.995, na.rm = TRUE)
      )
    }) |>
      dplyr::bind_rows() |>
      dplyr::mutate(
        p_adj = stats::p.adjust(p, method = "BH"),
        label = dplyr::case_when(
          p_adj < 1e-4 ~ "****",
          p_adj < 1e-3 ~ "***",
          p_adj < 1e-2 ~ "**",
          p_adj < 5e-2 ~ "*",
          TRUE ~ "ns"
        ),
        group1 = "Non_ECM",
        group2 = "ECM_Hotspot",
        y.position = ymax + ymax * 0.1,
        p.label = label
      )
    
    stat_anno <- stat_tbl
  }
  
  # ---- 3) Plot ----
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = group, y = y)) +
    ggplot2::geom_jitter(
      data = df_points,
      ggplot2::aes(color = group),
      width = 0.25, height = 0,
      size = 0.1, alpha = 0.3
    ) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = group),
      width = 0.8, trim = TRUE, scale = "width",
      alpha = 0.6, linewidth = 0.4
    ) +
    ggplot2::geom_boxplot(
      width = 0.15, outlier.shape = NA,
      fill = "white", alpha = 0.8, linewidth = 0.3
    ) +
    ggpubr::stat_pvalue_manual(
      stat_anno,
      label = "p.label",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      size = ifelse(mode == "panel", 2.5, 5),
      remove.bracket = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::scale_color_manual(values = fill_colors) +
    ggplot2::labs(x = NULL, y = ylab, title = title) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", color = "black"),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black"),
      legend.position = "none",
      panel.grid = ggplot2::element_blank()
    )
  
  if (mode == "panel") {
    p <- p + ggplot2::facet_grid(lr_pair ~ sample, scales = facet_scales)
  } else {
    p <- p + ggplot2::facet_wrap(~ lr_pair, ncol = 1, scales = facet_scales)
  }
  
  attr(p, "stat_table") <- stat_tbl
  return(p)
}
        
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

#' 绘制 Nature 风格的小提琴图 + 统计显著性标注
#'
#' @param data 绘图用的长格式数据 (dataframe)
#' @param stat_data 统计检验结果 (dataframe, 需包含 p.adj, group1, group2, y.position)
#' @param x_col X轴列名 (字符串, e.g. "TAM_Subtype")
#' @param y_col Y轴数值列名 (字符串, e.g. "dist_um")
#' @param fill_col 填充颜色列名 (字符串, 通常同 x_col)
#' @param facet_formula 分面公式 (formula, e.g. metric ~ sample)
#' @param palette 颜色向量 (named vector)
#' @param y_expand_mult Y轴顶部预留空间的倍数 (默认 0.18，用于放星号)
#' @param save_path (可选) 保存 PDF 的完整路径。如果为 NULL 则不保存。
#' @param width 保存宽度 (默认 6)
#' @param height 保存高度 (默认 5.2)
#'
#' @return ggplot 对象
plot_nature_violin <- function(
    data, 
    stat_data, 
    x_col, 
    y_col, 
    fill_col, 
    facet_formula = NULL, 
    palette = NULL, 
    y_expand_mult = 0.18,
    save_path = NULL,
    width = 6, 
    height = 5.2
) {
  
  # 1. 预处理：给统计表添加星号列 (p_star)
  # 确保 stat_data 只有在该列存在时才处理，避免重复
  stat_data_plot <- stat_data %>%
    mutate(
      p_star = case_when(
        p.adj <= 1e-4 ~ "****",
        p.adj <= 1e-3 ~ "***",
        p.adj <= 1e-2 ~ "**",
        p.adj <= 5e-2 ~ "*",
        TRUE          ~ "ns"
      )
    )
  
  # 2. 确保 x 轴是因子 (保持顺序)
  # 如果传入的数据中该列已经是因子，则保留；否则转为因子
  if (!is.factor(data[[x_col]])) {
    data[[x_col]] <- as.factor(data[[x_col]])
  }
  
  # 3. 绘图核心
  p <- ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], fill = .data[[fill_col]])) +
    # 小提琴层
    geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.90) +
    # 箱线图层 (细长风格)
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.55,
                 color = "grey15", size = 0.35) +
    # 中位数点
    stat_summary(fun = median, geom = "point", size = 1.4, color = "black") +
    
    # 颜色
    scale_fill_manual(values = palette) +
    
    # Y轴扩展 (防止星号被切)
    scale_y_continuous(expand = expansion(mult = c(0.02, y_expand_mult))) +
    
    # 统计显著性层
    stat_pvalue_manual(
      stat_data_plot,
      label = "p_star",
      xmin = "group1", # rstatix 默认列名
      xmax = "group2", # rstatix 默认列名
      y.position = "y.position", # 必须确保 stat_data 里有这列
      tip.length = 0.01,
      bracket.size = 0.35,
      size = 4.2,
      hide.ns = FALSE   # 是否隐藏 ns，可根据需求改为 TRUE
    ) +
    
    # 坐标轴与主题
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(color = "black", face = "bold", size = 10),
      axis.text.y = element_text(color = "black"),
      panel.grid = element_blank()
    )
  
  # 4. 处理分面 (如果有)
  if (!is.null(facet_formula)) {
    p <- p + facet_grid(facet_formula, scales = "free_y")
  }
  
  # 5. 保存逻辑
  if (!is.null(save_path)) {
    # 自动创建目录 (如果不存在)
    dir_name <- dirname(save_path)
    if (!dir.exists(dir_name)) dir.create(dir_name, recursive = TRUE)
    
    ggsave(save_path, plot = p, width = width, height = height, dpi = 600, useDingbats = FALSE)
    message("Plot saved to: ", save_path)
  }
  
  return(p)
}
plot_gene_boxplot <- function(expr_mat, 
                              meta_data, 
                              genes, 
                              group_col, 
                              pair_col = NULL,
                              palette_style = "npg", 
                              test_method = "t.test",
                              y_label = "Normalized Expression (log2)",
                              show_p_value = "p.format",
                              # --- 新增参数 ---
                              show_sample_labels = FALSE, 
                              label_size = 3) {
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(ggpubr)
    library(ggsci)
    library(ggrepel) # 必须加载这个包用于文字排版
  })
  
  # --- 1. 数据清洗 ---
  valid_genes <- intersect(genes, rownames(expr_mat))
  if (length(valid_genes) == 0) stop("❌ 错误：找不到指定基因！")
  
  sub_expr <- as.data.frame(t(expr_mat[valid_genes, , drop = FALSE]))
  sub_expr$Sample <- rownames(sub_expr)
  
  meta_use <- as.data.frame(meta_data)
  meta_use$Sample <- rownames(meta_use)
  df_plot <- merge(sub_expr, meta_use, by = "Sample")
  
  df_long <- df_plot %>%
    tidyr::pivot_longer(
      cols = all_of(valid_genes),
      names_to = "Gene", 
      values_to = "Expression"
    )
  
  if (!is.factor(df_long[[group_col]])) {
    df_long[[group_col]] <- as.factor(df_long[[group_col]])
  }
  
  sym_group <- sym(group_col)
  
  # --- 2. 基础绘图 ---
  p <- ggplot(df_long, aes(x = !!sym_group, y = Expression))
  
  if (!is.null(pair_col)) {
    sym_pair <- sym(pair_col)
    p <- p + geom_line(aes(group = !!sym_pair), color = "gray80", alpha = 0.6, linewidth = 0.3)
  }
  
  p <- p + geom_boxplot(aes(fill = !!sym_group), 
                        alpha = 0.8, 
                        width = 0.5, 
                        outlier.shape = NA,
                        size = 0.3)
  
  # 使用 position_jitter 固定随机散点位置，方便文字对齐
  pos <- position_jitter(width = 0.1, seed = 1)
  
  p <- p + geom_jitter(position = pos, size = 1.5, alpha = 0.7, 
                       color = "black", shape = 21, fill = "white")
  
  # --- 3. 新增：添加样本标签图层 ---
  if (show_sample_labels) {
    p <- p + geom_text_repel(
      aes(label = Sample),
      position = pos,     # 必须与 jitter 的 position 完全一致，文字才能对准点
      size = label_size,
      box.padding = 0.3,   # 文字距离点的距离
      point.padding = 0.3,
      min.segment.length = 0, # 总是显示引线
      segment.color = "grey50",
      force = 2,           # 标签排斥力度
      max.overlaps = Inf   # 强制显示所有标签
    )
  }
  
  # --- 4. 美化与配色 ---
  if (palette_style == "npg") {
    p <- p + scale_fill_npg() + scale_color_npg()
  } else if (palette_style == "aaas") {
    p <- p + scale_fill_aaas() + scale_color_aaas()
  } else if (palette_style == "nejm") {
    p <- p + scale_fill_nejm() + scale_color_nejm()
  } else if (palette_style == "lancet") {
    p <- p + scale_fill_lancet() + scale_color_lancet()
  } else {
    p <- p + scale_fill_brewer(palette = "Set1")
  }
  
  p <- p + theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(color = "black", size = 12, face = "bold"),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.y = element_text(color = "black", size = 13, face = "plain", margin = margin(r = 10)),
      axis.title.x = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(size = 13, face = "bold.italic")
    ) +
    labs(y = y_label)
  
  # --- 5. 统计分析 ---
  groups <- levels(df_long[[group_col]])
  if (length(groups) == 2) {
    comps <- list(c(groups[1], groups[2]))
    p <- p + stat_compare_means(
      comparisons = comps,
      method = test_method,
      label = show_p_value,
      size = 4,
      vjust = 0.5
    )
  }
  
  # --- 6. 分面处理 ---
  if (length(valid_genes) > 1) {
    p <- p + facet_wrap(~Gene, scales = "free_y", ncol = min(4, length(valid_genes)))
  }
  
  return(p)
}

plot_boxplot_v2 <- function(gsva_df,
                            score_vars = c("Classical", "Basal", "Diff_Score"),
                            cell_lines = NULL,
                            contrast_map = list(
                              Acinar = c("0h", "6h", "24h"),
                              Ductal = c("0h", "6h", "24h")
                            ),
                            treatment_col = "Time",
                            cell_line_col = "Cell_Type",
                            replicate_col = "Cell_Line",
                            sample_col = "sample",
                            palette_style = "npg",
                            connect_replicates = TRUE, 
                            test_method = "t.test",
                            p_label = "p.signif", 
                            y_label = "GSVA score",
                            title = NULL) {
  
  # 显式加载依赖，确保函数内可用
  suppressMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ggpubr)
    library(ggsci)
  })
  
  df <- gsva_df
  
  # 1. 自动过滤：根据 contrast_map 提取数据
  keep_idx <- rep(FALSE, nrow(df))
  for (cl in names(contrast_map)) {
    tr <- contrast_map[[cl]]
    # 显式使用 [[]] 获取列内容进行逻辑判断
    keep_idx <- keep_idx | (df[[cell_line_col]] == cl & df[[treatment_col]] %in% tr)
  }
  df <- df[keep_idx, , drop = FALSE]
  
  # 确保因子顺序
  df[[treatment_col]] <- factor(df[[treatment_col]], levels = unique(unlist(contrast_map)))
  
  # 2. 转为长格式 (使用 dplyr:: 前缀防止冲突)
  df_long <- df %>%
    dplyr::select(dplyr::all_of(c(sample_col, cell_line_col, treatment_col, replicate_col, score_vars))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(score_vars), names_to = "Score", values_to = "Value")
  
  # 3. 计算多组对比统计 (Ref vs Others)
  stat_list <- list()
  for (cl in names(contrast_map)) {
    tr <- contrast_map[[cl]]
    ref_group <- tr[1] 
    comp_groups <- tr[-1] 
    
    for (sc in score_vars) {
      for (target in comp_groups) {
        sub_data <- df_long %>% dplyr::filter(Score == sc, .data[[cell_line_col]] == cl)
        
        # 提取数值
        x <- sub_data$Value[sub_data[[treatment_col]] == ref_group]
        y <- sub_data$Value[sub_data[[treatment_col]] == target]
        
        if(length(x) >= 2 && length(y) >= 2) {
          pval <- if(test_method == "t.test") t.test(x, y)$p.value else wilcox.test(x, y)$p.value
          
          # 格式化星号
          lab <- if(p_label == "p.signif") {
            if(pval < 0.001) "***" else if(pval < 0.01) "**" else if(pval < 0.05) "*" else "ns"
          } else paste0("p=", formatC(pval, format = "f", digits = 3))
          
          # 确定 Y 轴位置 (阶梯式升高以防重叠)
          y_base <- max(sub_data$Value, na.rm = TRUE)
          y_range <- max(sub_data$Value, na.rm = TRUE) - min(sub_data$Value, na.rm = TRUE)
          y_pos <- y_base + (which(comp_groups == target) * 0.15 * y_range)
          
          stat_list[[length(stat_list) + 1]] <- data.frame(
            Cell = cl, Score = sc, group1 = ref_group, group2 = target,
            p = pval, label = lab, y.position = y_pos
          )
        }
      }
    }
  }
  
  df_p <- dplyr::bind_rows(stat_list)
  if(nrow(df_p) > 0) colnames(df_p)[1] <- cell_line_col
  
  # 4. 绘图
  # 使用 sym() 转换字符串为符号
  p <- ggplot(df_long, aes(x = .data[[treatment_col]], y = Value))
  
  # 连线层 (放在箱线图下面)
  if (connect_replicates && !is.null(replicate_col)) {
    p <- p + geom_line(aes(group = .data[[replicate_col]]), color = "gray80", alpha = 0.6)
  }
  
  p <- p + 
    geom_boxplot(aes(fill = .data[[treatment_col]]), alpha = 0.7, outlier.shape = NA, width = 0.6, size = 0.4) +
    geom_jitter(width = 0.1, shape = 21, fill = "white", size = 1.8, stroke = 0.5) +
    facet_grid(Score ~ .data[[cell_line_col]], scales = "free_y") +
    theme_classic(base_size = 14) +
    theme(strip.background = element_blank(), 
          strip.text = element_text(face = "bold", size = 12),
          panel.grid.major.y = element_line(color = "gray95"),
          legend.position = "none", 
          axis.title.x = element_blank()) +
    labs(y = y_label, title = title)
  
  # 应用 NPG 配色
  if(palette_style == "npg") p <- p + ggsci::scale_fill_npg()
  
  # 添加统计标签
  if(nrow(df_p) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(df_p, label = "label", xmin = "group1", xmax = "group2", 
                                        tip.length = 0.01, size = 3.5)
  }
  
  return(p)
}
plot_boxplot <- function(gsva_df,
                         score_vars = c("Classical", "Basal"),
                         # Which cell lines to include (NULL = all)
                         cell_lines = NULL,
                         # Per-cell-line contrast definition: list(cell_line = c("Control","Drug"))
                         contrast_map = list(
                           Suit007 = c("Control", "MRTX1133"),
                           PDC145  = c("Control", "RMC6236")
                         ),
                         treatment_col = "treatment",
                         cell_line_col = "cell_line",
                         replicate_col = "replicate",
                         sample_col = "sample",
                         # Plot look
                         palette_style = "npg",
                         show_points = TRUE,
                         point_size = 1.8,
                         point_alpha = 0.75,
                         box_alpha = 0.85,
                         box_width = 0.55,
                         connect_replicates = FALSE,   # draws lines by replicate across treatments (use only if truly paired)
                         line_alpha = 0.45,
                         line_size = 0.35,
                         # Stats
                         test_method = c("wilcox.test", "t.test"),
                         paired = FALSE,
                         p_label = c("p.format", "p.signif"),
                         hide_ns = TRUE,
                         p_y_pad_frac = 0.10,          # how high above max to place p label
                         # Facet
                         facet_grid = TRUE,            # TRUE: Score x CellLine grid; FALSE: facet_wrap by Score
                         facet_ncol = 2,
                         # Background by cell line (soft)
                         add_cellline_bg = F,       # requires ggnewscale; if not installed, auto-fallback to FALSE
                         bg_alpha = 0.06,
                         bg_colors = NULL,             # named vector: c(PDC145="#...", Suit007="#...")
                         # Labels
                         y_label = "GSVA score",
                         title = NULL) {
  # -----------------------------
  # Dependencies
  # -----------------------------
  req_pkgs <- c("ggplot2", "dplyr", "tidyr", "rlang")
  for (p in req_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
  }
  if (!requireNamespace("ggsci", quietly = TRUE)) stop("Missing package: ggsci")
  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Missing package: ggpubr")
  
  test_method <- match.arg(test_method)
  p_label <- match.arg(p_label)
  
  # -----------------------------
  # Input checks
  # -----------------------------
  if (!is.data.frame(gsva_df)) stop("gsva_df must be a data.frame.")
  needed_cols <- c(score_vars, treatment_col, cell_line_col)
  miss <- setdiff(needed_cols, colnames(gsva_df))
  if (length(miss) > 0) stop("Missing columns in gsva_df: ", paste(miss, collapse = ", "))
  
  df <- gsva_df
  
  # Optional filter cell lines
  if (!is.null(cell_lines)) {
    df <- df[df[[cell_line_col]] %in% cell_lines, , drop = FALSE]
  }
  
  # Apply contrast_map filter (Control vs Drug per cell line)
  if (!is.null(contrast_map) && length(contrast_map) > 0) {
    keep_idx <- rep(FALSE, nrow(df))
    for (cl in names(contrast_map)) {
      tr <- contrast_map[[cl]]
      if (length(tr) != 2) stop("contrast_map[[", cl, "]] must be length-2: c('Control','Drug')")
      keep_idx <- keep_idx | (df[[cell_line_col]] == cl & df[[treatment_col]] %in% tr)
    }
    df <- df[keep_idx, , drop = FALSE]
  }
  
  if (nrow(df) == 0) stop("No data left after filtering. Check cell_lines / contrast_map / column names.")
  
  # Ensure factors for ordering (within each facet we free_x, so order mainly affects comparisons)
  df[[treatment_col]] <- as.factor(df[[treatment_col]])
  df[[cell_line_col]] <- as.factor(df[[cell_line_col]])
  if (!is.null(replicate_col) && replicate_col %in% colnames(df)) {
    df[[replicate_col]] <- as.factor(df[[replicate_col]])
  }
  
  # Long format
  df_long <- df |>
    dplyr::select(dplyr::all_of(c(sample_col, cell_line_col, treatment_col, replicate_col, score_vars))) |>
    tidyr::pivot_longer(cols = dplyr::all_of(score_vars),
                        names_to = "Score",
                        values_to = "Value")
  
  # -----------------------------
  # Compute per-panel p-values (robust; avoids missing-group warnings)
  # -----------------------------
  # Helper to format p label
  fmt_p <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p_label == "p.signif") {
      if (p < 0.0001) "****" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
    } else {
      # p.format
      if (p < 0.0001) "p<1e-4" else paste0("p=", formatC(p, format = "f", digits = 3))
    }
  }
  
  # Build comparisons table from contrast_map that actually exist in df_long
  comps <- data.frame()
  if (!is.null(contrast_map) && length(contrast_map) > 0) {
    for (cl in names(contrast_map)) {
      tr <- contrast_map[[cl]]
      comps <- rbind(comps, data.frame(cell_line = cl, group1 = tr[1], group2 = tr[2]))
    }
  } else {
    # If no contrast_map, do nothing (no p-values)
    comps <- data.frame(cell_line = character(), group1 = character(), group2 = character())
  }
  
  # Map column names
  names(comps)[1] <- cell_line_col
  
  # Compute p per (cell_line, Score, comparison)
  df_p <- NULL
  if (nrow(comps) > 0) {
    # Precompute per panel y positions
    panel_max <- df_long |>
      dplyr::group_by(.data[[cell_line_col]], .data[["Score"]]) |>
      dplyr::summarise(ymax = max(Value, na.rm = TRUE),
                       ymin = min(Value, na.rm = TRUE),
                       .groups = "drop") |>
      dplyr::mutate(yr = dplyr::if_else(is.finite(ymax - ymin), ymax - ymin, 0),
                    y_pos = ymax + p_y_pad_frac * dplyr::if_else(yr == 0, 1, yr))
    
    out_list <- list()
    k <- 1
    
    for (i in seq_len(nrow(comps))) {
      cl <- comps[[cell_line_col]][i]
      g1 <- comps$group1[i]
      g2 <- comps$group2[i]
      
      sub <- df_long |>
        dplyr::filter(.data[[cell_line_col]] == cl,
                      .data[[treatment_col]] %in% c(g1, g2))
      
      if (nrow(sub) == 0) next
      
      # Per Score
      for (sc in unique(sub$Score)) {
        sub_sc <- sub[sub$Score == sc, , drop = FALSE]
        
        # Skip if any group missing
        if (!all(c(g1, g2) %in% unique(sub_sc[[treatment_col]]))) next
        
        # Paired test via replicate matching (only if replicate_col exists)
        pval <- NA_real_
        if (paired && !is.null(replicate_col) && replicate_col %in% colnames(sub_sc)) {
          wide <- sub_sc |>
            dplyr::select(dplyr::all_of(c(replicate_col, treatment_col)), Value) |>
            tidyr::pivot_wider(names_from = dplyr::all_of(treatment_col), values_from = Value)
          
          if (all(c(g1, g2) %in% colnames(wide))) {
            x <- wide[[g1]]
            y <- wide[[g2]]
            ok <- is.finite(x) & is.finite(y)
            x <- x[ok]; y <- y[ok]
            if (length(x) >= 2) {
              pval <- if (test_method == "t.test") stats::t.test(x, y, paired = TRUE)$p.value
              else stats::wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value
            }
          }
        } else {
          x <- sub_sc$Value[sub_sc[[treatment_col]] == g1]
          y <- sub_sc$Value[sub_sc[[treatment_col]] == g2]
          x <- x[is.finite(x)]; y <- y[is.finite(y)]
          if (length(x) >= 2 && length(y) >= 2) {
            pval <- if (test_method == "t.test") stats::t.test(x, y)$p.value
            else stats::wilcox.test(x, y, exact = FALSE)$p.value
          }
        }
        
        lab <- fmt_p(pval)
        if (hide_ns && identical(lab, "ns")) next
        
        yp <- panel_max |>
          dplyr::filter(.data[[cell_line_col]] == cl, .data[["Score"]] == sc) |>
          dplyr::pull(y_pos)
        if (length(yp) == 0) yp <- max(sub_sc$Value, na.rm = TRUE)
        
        out_list[[k]] <- data.frame(
          cell_line = cl, Score = sc, group1 = g1, group2 = g2,
          p = pval, label = lab, y.position = yp
        )
        names(out_list[[k]])[1] <- cell_line_col
        k <- k + 1
      }
    }
    
    if (length(out_list) > 0) {
      df_p <- dplyr::bind_rows(out_list)
    }
  }
  
  # -----------------------------
  # Plot
  # -----------------------------
  library(ggplot2)
  
  # Base plot
  p <- ggplot(df_long, aes(x = .data[[treatment_col]], y = Value))
  
  # Optional soft background by cell line (uses ggnewscale to avoid fill-scale conflict)
  has_ggnewscale <- requireNamespace("ggnewscale", quietly = TRUE)
  if (add_cellline_bg && has_ggnewscale) {
    if (is.null(bg_colors)) {
      levs <- levels(df_long[[cell_line_col]])
      # default soft colors (no need to specify; will be set via scale_fill_manual)
      bg_colors <- stats::setNames(rep("#B0B0B0", length(levs)), levs)
      # If you want custom, pass bg_colors explicitly.
    }
    bg_df <- df_long |>
      dplyr::distinct(.data[[cell_line_col]]) |>
      dplyr::mutate(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
    
    p <- p +
      geom_rect(
        data = bg_df,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = .data[[cell_line_col]]),
        inherit.aes = FALSE,
        alpha = bg_alpha,
        color = NA
      ) +
      scale_fill_manual(values = bg_colors, guide = "none") +
      ggnewscale::new_scale_fill()
  }
  
  # Paired lines (optional)
  if (connect_replicates && !is.null(replicate_col) && replicate_col %in% colnames(df_long)) {
    p <- p + geom_line(aes(group = interaction(.data[[cell_line_col]], Score, .data[[replicate_col]])),
                       alpha = line_alpha, linewidth = line_size, color = "gray70")
  }
  
  # Box + points
  p <- p +
    geom_boxplot(aes(fill = .data[[treatment_col]]),
                 alpha = box_alpha, width = box_width,
                 outlier.shape = NA, linewidth = 0.35)
  
  if (show_points) {
    p <- p + geom_jitter(width = 0.10, height = 0,
                         size = point_size, alpha = point_alpha,
                         shape = 21, color = "black", fill = "white")
  }
  
  # Colors
  if (palette_style == "npg") {
    p <- p + ggsci::scale_fill_npg()
  } else if (palette_style == "aaas") {
    p <- p + ggsci::scale_fill_aaas()
  } else if (palette_style == "nejm") {
    p <- p + ggsci::scale_fill_nejm()
  } else if (palette_style == "lancet") {
    p <- p + ggsci::scale_fill_lancet()
  } else {
    p <- p + scale_fill_brewer(palette = "Set1")
  }
  
  # Facet
  if (facet_grid) {
    p <- p + facet_grid(rows = vars(Score),
                        cols = vars(!!rlang::sym(cell_line_col)),
                        scales = "free_x", space = "free_x")
  } else {
    p <- p + facet_wrap(~Score, ncol = facet_ncol, scales = "free_y")
  }
  
  # Stats annotation
  if (!is.null(df_p) && nrow(df_p) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      df_p,
      label = "label",
      y.position = "y.position",
      xmin = "group1",
      xmax = "group2",
      tip.length = 0.01,
      size = 3.8
    )
  }
  
  # Theme
  p <- p +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(color = "black", size = 11, face = "bold"),
      axis.text.y = element_text(color = "black", size = 11),
      axis.title.y = element_text(color = "black", size = 13, margin = margin(r = 8)),
      axis.title.x = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(y = y_label, title = title)
  
  return(p)
}



# =========================== Flexible Heatmap (ComplexHeatmap) ===========================
flexi_heatmap <- function(
    data,
    row_id_col         = NULL,         # if data.frame: column holding row IDs; if NULL and has rownames, use them
    annotation_df      = NULL,         # optional row annotation table (e.g., cell -> Group)
    annotation_id_col  = NULL,         # which column in annotation_df matches row IDs; required if annotation_df is given
    annotation_keep    = NULL,         # which columns from annotation_df to display; if NULL, use all except the id
    # ----- value handling -----
    transform          = c("none","log2","log1p"),  # transform BEFORE scaling
    scale_mode         = c("none","row_z","col_z","row_minmax","col_minmax"), # scaling AFTER transform
    use_original       = FALSE,        # if TRUE, force scale_mode="none" (keep original magnitude)
    cap_z              = 2.5,          # cap for z-scores
    na_strategy        = c("none","impute_row_median","impute_col_median","impute_global_median"),
    # ----- clustering & distance -----
    cluster_rows       = TRUE,
    cluster_cols       = TRUE,
    dist_rows          = c("euclidean","pearson","spearman","cosine","manhattan"),
    dist_cols          = c("euclidean","pearson","spearman","cosine","manhattan"),
    hclust_method      = "ward.D2",
    row_split_by       = NULL,         # a column name from annotation_df to split rows
    col_split_by       = NULL,         # split columns by a vector (length = ncol) or NULL
    # ----- aesthetics -----
    palette            = c("blue-white-red","navy-white-firebrick","viridis","magma","plasma","terrain"),
    legend_name        = NULL,
    na_col             = "#F1F1F1",
    show_row_names     = TRUE,
    show_column_names  = TRUE,
    row_names_gp       = grid::gpar(fontsize = 9),
    column_names_gp    = grid::gpar(fontsize = 10),
    label_wrap         = 28,           # wrap long column labels (approx chars per line)
    annotate_values    = FALSE,        # overlay numbers for small matrices
    value_digits       = 2,
    ann_colors         = NULL,         # NEW: colors for row annotations (named list: list(ColName = c(Level=Color,...)))
    # ----- sizing & export -----
    width              = NULL,         # default auto by data size (mm)
    height             = NULL,         # default auto by data size (mm)
    units              = "mm",
    dpi                = 600,
    save               = TRUE,
    outdir             = "AUC_HC_WardD2_out",
    prefix             = "AUC_Heatmap",
    # ----- restrictions / return -----
    col_only           = FALSE,        # NEW: if TRUE, row_* scaling is auto-converted to col_*
    return_ht          = TRUE          # return the Heatmap object
){
  library(ComplexHeatmap); library(circlize)
  
  # ---- explicit matches (fixes "argument is missing, with no default") ----
  transform   <- match.arg(transform,   c("none","log2","log1p"))
  scale_mode  <- match.arg(scale_mode,  c("none","row_z","col_z","row_minmax","col_minmax"))
  dist_rows   <- match.arg(dist_rows,   c("euclidean","pearson","spearman","cosine","manhattan"))
  dist_cols   <- match.arg(dist_cols,   c("euclidean","pearson","spearman","cosine","manhattan"))
  na_strategy <- match.arg(na_strategy, c("none","impute_row_median","impute_col_median","impute_global_median"))
  palette     <- match.arg(palette,     c("blue-white-red","navy-white-firebrick","viridis","magma","plasma","terrain"))
  
  # use_original overrides any scaling
  if (isTRUE(use_original)) scale_mode <- "none"
  
  # enforce column-only normalization if requested
  if (isTRUE(col_only) && scale_mode %in% c("row_z","row_minmax")) {
    scale_mode <- sub("^row_", "col_", scale_mode)
    message("[flexi_heatmap] col_only=TRUE: switched scale_mode to ", scale_mode)
  }
  
  # ---------- helpers ----------
  wrap_label <- function(x, k=28){
    vapply(as.character(x), function(s){
      if (nchar(s) <= k) return(s)
      s <- gsub("([_/.-])", "\\1 ", s)   # add break opportunities
      parts <- strsplit(s, " +")[[1]]
      out <- ""; cur <- ""
      for (w in parts){
        if (nchar(paste(cur, w)) <= k) { cur <- paste(cur, w) }
        else { out <- paste0(out, trimws(cur), "\n"); cur <- w }
      }
      paste0(out, trimws(cur))
    }, FUN.VALUE = character(1))
  }
  
  cosine_dist <- function(M, by="row"){
    X <- if (by == "row") M else t(M)
    X[is.na(X)] <- 0
    nrm <- sqrt(rowSums(X^2)); nrm[nrm == 0] <- 1
    Xn <- X / nrm
    as.dist(1 - Xn %*% t(Xn))
  }
  corr_dist <- function(M, method=c("pearson","spearman"), by="row"){
    method <- match.arg(method)
    X <- if (by == "row") M else t(M)
    C <- suppressWarnings(cor(t(X), method = method, use = "pairwise.complete.obs"))
    C[is.na(C)] <- 0
    as.dist(1 - C)
  }
  
  to_numeric_matrix <- function(d, row_id_col){
    if (is.matrix(d)) {
      mat <- d
      rn <- rownames(mat)
      if (is.null(rn)) stop("Matrix must have rownames or provide row_id_col.")
    } else {
      df <- as.data.frame(d, check.names = FALSE)
      if (!is.null(row_id_col)) {
        stopifnot(row_id_col %in% names(df))
        rn <- as.character(df[[row_id_col]])
        df[[row_id_col]] <- NULL
      } else if (!is.null(rownames(df))) {
        rn <- rownames(df)
      } else {
        stop("Provide row_id_col when data.frame has no rownames.")
      }
      # coerce numeric
      num <- suppressWarnings(as.data.frame(lapply(df, function(x) as.numeric(as.character(x))), check.names = FALSE))
      keep <- vapply(num, function(v) sum(!is.na(v)) > 0, logical(1))
      num <- num[keep]
      if (ncol(num) < 2) stop("Need at least 2 numeric columns to draw a heatmap.")
      mat <- as.matrix(num)
      rownames(mat) <- rn
    }
    # drop all-NA rows/cols
    row_keep <- rowSums(is.na(mat)) < ncol(mat)
    col_keep <- colSums(is.na(mat)) < nrow(mat)
    mat <- mat[row_keep, col_keep, drop = FALSE]
    if (nrow(mat) == 0 || ncol(mat) == 0) stop("Matrix is empty after removing all-NA rows/cols.")
    mat
  }
  
  impute_matrix <- function(M, how){
    if (how == "none") return(M)
    if (how == "impute_row_median") {
      for (i in seq_len(nrow(M))) {
        mi <- M[i, ]
        if (anyNA(mi)) {
          med <- stats::median(mi, na.rm = TRUE); if (!is.finite(med)) med <- 0
          mi[is.na(mi)] <- med; M[i, ] <- mi
        }
      }
    } else if (how == "impute_col_median") {
      for (j in seq_len(ncol(M))) {
        mj <- M[, j]
        if (anyNA(mj)) {
          med <- stats::median(mj, na.rm = TRUE); if (!is.finite(med)) med <- 0
          mj[is.na(mj)] <- med; M[, j] <- mj
        }
      }
    } else if (how == "impute_global_median") {
      med <- stats::median(M, na.rm = TRUE); if (!is.finite(med)) med <- 0
      M[is.na(M)] <- med
    }
    M
  }
  
  apply_transform <- function(M, mode){
    if (mode == "none") return(M)
    if (mode == "log2")  return(log2(pmax(M, 0) + 1))
    if (mode == "log1p") return(log1p(pmax(M, 0)))
    stop("Unknown transform: ", mode)
  }
  apply_scale <- function(M, mode, cap_z){
    if (mode == "none") return(M)
    if (mode == "row_z") {
      mu <- rowMeans(M, na.rm = TRUE)
      sdv <- matrixStats::rowSds(M, na.rm = TRUE)
      sdv[sdv == 0 | !is.finite(sdv)] <- 1
      Z <- (M - mu) / sdv
      Z[Z >  cap_z] <-  cap_z
      Z[Z < -cap_z] <- -cap_z
      return(Z)
    }
    if (mode == "col_z") {
      mu <- colMeans(M, na.rm = TRUE)
      sdv <- matrixStats::colSds(M, na.rm = TRUE)
      sdv[sdv == 0 | !is.finite(sdv)] <- 1
      Z <- sweep(M, 2, mu, "-")
      Z <- sweep(Z, 2, sdv, "/")
      Z[Z >  cap_z] <-  cap_z
      Z[Z < -cap_z] <- -cap_z
      return(Z)
    }
    if (mode == "row_minmax") {
      rng <- matrixStats::rowRanges(M)
      span <- rng[,2] - rng[,1]; span[span == 0] <- 1
      return((M - rng[,1]) / span)
    }
    if (mode == "col_minmax") {
      mn <- apply(M, 2, min, na.rm = TRUE); mx <- apply(M, 2, max, na.rm = TRUE)
      span <- mx - mn; span[span == 0] <- 1
      Z <- sweep(M, 2, mn, "-")
      Z <- sweep(Z, 2, span, "/")
      return(Z)
    }
    stop("Unknown scale_mode: ", mode)
  }
  
  pick_col_fun <- function(M, palette, cap_z, mode){
    if (mode %in% c("row_z","col_z")) {
      # symmetric around 0
      cols <- switch(
        palette,
        "blue-white-red"       = c("#2C7BB6","white","#D7191C"),
        "navy-white-firebrick" = c("navy","white","firebrick3"),
        "viridis"              = grDevices::colorRampPalette(c("#440154","#21908C","#FDE725"))(3),
        "magma"                = grDevices::colorRampPalette(c("#000004","#B53679","#FCFDBF"))(3),
        "plasma"               = grDevices::colorRampPalette(c("#0D0887","#CC4778","#F0F921"))(3),
        "terrain"              = grDevices::colorRampPalette(c("#0B3C49","#88B365","#FDE725"))(3)
      )
      circlize::colorRamp2(c(-cap_z, 0, cap_z), cols)
    } else {
      rng <- range(M, na.rm = TRUE); if (!is.finite(rng[1])) rng <- c(0,1)
      cols <- switch(
        palette,
        "blue-white-red"       = grDevices::colorRampPalette(c("#2C7BB6","#E0ECF4","#FDD9C4","#D7191C"))(11),
        "navy-white-firebrick" = grDevices::colorRampPalette(c("navy","#F2F2F2","firebrick3"))(11),
        "viridis"              = grDevices::colorRampPalette(c("#440154","#31688E","#35B779","#FDE725"))(11),
        "magma"                = grDevices::colorRampPalette(c("#000004","#51127C","#B63679","#FBFCBF"))(11),
        "plasma"               = grDevices::colorRampPalette(c("#0D0887","#9C179E","#ED7953","#F0F921"))(11),
        "terrain"              = grDevices::colorRampPalette(c("#00441B","#006D2C","#74C476","#F7FCF5"))(11)
      )
      circlize::colorRamp2(seq(rng[1], rng[2], length.out = length(cols)), cols)
    }
  }
  
  build_dend <- function(M, axis=c("row","col"), dist_mode="euclidean", method="ward.D2"){
    axis <- match.arg(axis)
    if (dist_mode %in% c("pearson","spearman")) {
      d <- corr_dist(M, method = dist_mode, by = ifelse(axis=="row","row","col"))
    } else if (dist_mode == "cosine") {
      d <- cosine_dist(M, by = ifelse(axis=="row","row","col"))
    } else {
      X <- if (axis == "row") M else t(M)
      d <- dist(X, method = dist_mode)
    }
    stats::as.dendrogram(stats::hclust(d, method = method))
  }
  
  # -------------------- assemble matrix --------------------
  mat <- to_numeric_matrix(data, row_id_col = row_id_col)
  
  # -------------------- row annotation (optional) --------------------
  row_ha <- NULL; row_split <- NULL
  if (!is.null(annotation_df)) {
    if (is.null(annotation_id_col)) stop("Provide annotation_id_col when annotation_df is given.")
    ann <- as.data.frame(annotation_df, check.names = FALSE)
    stopifnot(annotation_id_col %in% names(ann))
    rid <- rownames(mat)
    mm  <- match(rid, as.character(ann[[annotation_id_col]]))
    ann2 <- ann[mm, , drop = FALSE]
    rownames(ann2) <- rid
    if (is.null(annotation_keep)) {
      annotation_keep <- setdiff(colnames(ann2), annotation_id_col)
    }
    if (length(annotation_keep)) {
      ann_show <- ann2[, annotation_keep, drop = FALSE]
      # ensure factors for discrete columns
      ann_show[] <- lapply(ann_show, function(x) if (is.numeric(x)) x else factor(x))
      
      # build color list from ann_colors if provided
      col_list <- NULL
      if (!is.null(ann_colors)) {
        col_list <- list()
        for (nm in names(ann_colors)) {
          if (nm %in% colnames(ann_show)) col_list[[nm]] <- ann_colors[[nm]]
        }
        if (!length(col_list)) col_list <- NULL
      }
      
      row_ha <- ComplexHeatmap::rowAnnotation(
        df = ann_show,
        col = col_list,  # pass user-defined colors
        annotation_name_gp = grid::gpar(fontface = "bold")
      )
      
      # row split by a column in annotation
      if (!is.null(row_split_by)) {
        stopifnot(row_split_by %in% colnames(ann_show))
        row_split <- ann_show[[row_split_by]]
      }
    }
  }
  
  # -------------------- transform / impute / scale --------------------
  mat <- apply_transform(mat, transform)
  if (na_strategy != "none") mat <- impute_matrix(mat, na_strategy)
  mat <- apply_scale(mat, scale_mode, cap_z)
  
  # -------------------- clustering (robust to split vectors) --------------------
  is_scalar_number <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  
  split_rows_is_vec <- !is.null(row_split)    && !is_scalar_number(row_split)
  split_cols_is_vec <- !is.null(col_split_by) && !is_scalar_number(col_split_by)
  
  if (split_rows_is_vec && length(row_split) != nrow(mat)) {
    stop("row_split is a vector but its length != nrow(mat).")
  }
  if (split_cols_is_vec && length(col_split_by) != ncol(mat)) {
    stop("column_split is a vector but its length != ncol(mat).")
  }
  
  if (isTRUE(cluster_rows)) {
    row_dend <- if (split_rows_is_vec) TRUE else build_dend(mat, "row", dist_rows, hclust_method)
  } else {
    row_dend <- FALSE
  }
  if (isTRUE(cluster_cols)) {
    col_dend <- if (split_cols_is_vec) TRUE else build_dend(mat, "col", dist_cols, hclust_method)
  } else {
    col_dend <- FALSE
  }
  
  # -------------------- colors --------------------
  col_fun <- pick_col_fun(mat, palette, cap_z, scale_mode)
  if (is.null(legend_name)) legend_name <- if (scale_mode %in% c("row_z","col_z")) "Z" else "Value"
  
  # -------------------- column labels (wrapping) --------------------
  cn <- colnames(mat)
  if (!is.null(label_wrap) && isTRUE(label_wrap > 0)) cn <- wrap_label(cn, k = label_wrap)
  
  # -------------------- small-matrix value overlay --------------------
  cell_fun <- NULL
  if (isTRUE(annotate_values) && nrow(mat) <= 60 && ncol(mat) <= 40) {
    cell_fun <- function(j, i, x, y, w, h, fill){
      v <- mat[i, j]
      if (is.na(v)) return(NULL)
      grid::grid.text(sprintf(paste0("%.", value_digits, "f"), v), x, y, gp = grid::gpar(col = "black", fontsize = 7))
    }
  }
  
  # -------------------- auto size (mm) --------------------
  if (is.null(width))  width  <- max(120, 50 + 6 * ncol(mat))  # ~6 mm per column
  if (is.null(height)) height <- max(100, 50 + 4 * nrow(mat))  # ~4 mm per row
  
  # -------------------- build heatmap --------------------
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name                = legend_name,
    col                 = col_fun,
    na_col              = na_col,
    cluster_rows        = row_dend,
    cluster_columns     = col_dend,
    row_split           = row_split,
    column_split        = col_split_by,
    show_row_names      = show_row_names,
    show_column_names   = show_column_names,
    row_names_gp        = row_names_gp,
    column_names_gp     = column_names_gp,
    column_names_side   = "top",
    column_labels       = cn,  # keep using column_labels for custom labels
    heatmap_legend_param = list(title_gp = grid::gpar(fontface = "bold")),
    cell_fun            = cell_fun,
    border              = TRUE
  )
  
  if (!is.null(row_ha)) {
    ht <- row_ha + ht
  }
  
  # -------------------- export --------------------
  if (isTRUE(save)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    mm_to_in <- function(x) as.numeric(x)/25.4
    
    pdf(file.path(outdir, paste0(prefix, "_heatmap.pdf")),
        width = mm_to_in(width), height = mm_to_in(height))
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
    grDevices::dev.off()
    
    png(file.path(outdir, paste0(prefix, "_heatmap.png")),
        width = width, height = height, units = units, res = dpi, type = "cairo")
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
    grDevices::dev.off()
  }
  
  if (isTRUE(return_ht)) return(ht) else invisible(NULL)
}
# ======================== (end) Flexible Heatmap (ComplexHeatmap) ========================

# ============== 工具函数 ==============
pick_best <- function(base){
  cand <- list.files(pattern = paste0("^", base, "\\.(xlsx|csv|xls)$"), ignore.case = TRUE)
  if (length(cand) == 0) stop("找不到文件：", base, ".xlsx/csv/xls")
  ord <- order(match(tolower(tools::file_ext(cand)), c("xlsx","csv","xls")))
  cand[ord][1]
}

# 广谱识别列名（返回实际列名或 NA）
find_col <- function(nms, keywords_exact = character(), keywords_regex = character()){
  nms_low <- tolower(nms)
  # 先精确关键词
  for (k in tolower(keywords_exact)) {
    hit <- which(nms_low == k)
    if (length(hit)) return(nms[hit[1]])
  }
  # 再正则关键词
  for (rgx in keywords_regex) {
    hit <- which(str_detect(nms_low, rgx))
    if (length(hit)) return(nms[hit[1]])
  }
  NA_character_
}

# 判断列名是否明显为注释列（非样本列）
is_annotation_name <- function(colname){
  pat <- "(^len$|length|gene_length|biotype|type$|class$|category$|chr|chromosome|start$|end$|strand$|gc$|desc|description|biotype|tx|transcript|exon|fpkm|tpm|rpkm)"
  str_detect(tolower(colname), pat)
}

# 读取一个 counts 文件为标准化 tibble：Symbol + 样本列（不改样本名）
read_counts_generic <- function(path){
  ext <- tolower(tools::file_ext(path))
  df <- switch(ext,
               "xlsx" = readxl::read_excel(path, .name_repair = "minimal"),
               "xls"  = readxl::read_excel(path, .name_repair = "minimal"),
               "csv"  = readr::read_csv(path, show_col_types = FALSE),
               stop("不支持的扩展名：", ext))
  df <- as_tibble(df)
  names(df) <- trimws(names(df))
  
  # 识别 Symbol / #ID 列
  id_col <- find_col(names(df),
                     keywords_exact = c("#id","gene_id","id"),
                     keywords_regex = c("ensembl", "^ensg", "gene[_\\.]?id"))
  sym_col <- find_col(names(df),
                      keywords_exact = c("symbol","gene_symbol","gene","gene_name","hgnc_symbol"),
                      keywords_regex = c("external[_\\.]?gene[_\\.]?name", "symbol"))
  
  if (is.na(sym_col) && DROP_ROWS_WITHOUT_SYMBOL) {
    message("警告：未找到 Symbol 列；由于 DROP_ROWS_WITHOUT_SYMBOL=TRUE，将丢弃无 Symbol 的行（即使存在 #ID）。")
  }
  if (is.na(sym_col) && !DROP_ROWS_WITHOUT_SYMBOL && is.na(id_col)) {
    stop("既没有可用的 Symbol 也没有 #ID，无法定位基因。")
  }
  
  # 标注注释列
  anno_cols <- names(df)[vapply(names(df), is_annotation_name, logical(1))]
  anno_cols <- union(anno_cols, c(sym_col, id_col))
  sample_cols <- setdiff(names(df), anno_cols)
  
  # 强制样本列为数值；把完全非数值的列剔除（保留 0/整数/小数）
  keep_sample <- logical(length(sample_cols))
  for (i in seq_along(sample_cols)) {
    xi <- suppressWarnings(as.numeric(as.character(df[[ sample_cols[i] ]])))
    # 保留：不全是 NA；且不是全空字符
    keep_sample[i] <- !all(is.na(xi))
    if (keep_sample[i]) df[[ sample_cols[i] ]] <- xi
  }
  sample_cols <- sample_cols[keep_sample]
  
  # 仅保留 Symbol + 样本列（以及必要时 #ID）
  out <- tibble::tibble(Symbol = if (!is.na(sym_col)) as.character(df[[sym_col]]) else NA_character_)
  if (!is.na(id_col)) out <- mutate(out, `#ID` = as.character(df[[id_col]]))
  
  if (length(sample_cols) == 0) stop("没有检测到任何数值型样本列：", path)
  
  out <- bind_cols(out, df[, sample_cols, drop = FALSE])
  
  # 丢弃无 Symbol 的行或用 #ID 顶上
  if (DROP_ROWS_WITHOUT_SYMBOL) {
    out <- filter(out, !is.na(Symbol) & Symbol != "")
  } else {
    if (is.na(sym_col)) {
      out <- mutate(out, Symbol = `#ID`)
    } else {
      out <- mutate(out, Symbol = if_else(is.na(Symbol) | Symbol == "",
                                          if ("#ID" %in% names(out)) `#ID` else Symbol,
                                          Symbol))
    }
    out <- filter(out, !is.na(Symbol) & Symbol != "")
  }
  
  # 处理重复 Symbol（同一文件中）
  fun <- switch(DUPLICATE_SYMBOL_RULE,
                "sum"   = ~ sum(replace_na(.x, 0)),
                "max"   = ~ max(replace_na(.x, 0)),
                "first" = ~ .x[which.max(!is.na(.x))],
                stop("未知的 DUPLICATE_SYMBOL_RULE: ", DUPLICATE_SYMBOL_RULE))
  out <- out %>%
    group_by(Symbol) %>%
    summarise(across(all_of(sample_cols), fun), .groups = "drop")
  
  out
}





suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(stringr)
  library(tibble); library(writexl); library(purrr)
})

merge_counts_for_deseq2 <- function(
    dir = ".",
    files = NULL,
    bases = NULL,
    group_labels = NULL,
    drop_rows_without_symbol = TRUE,
    duplicate_symbol_rule    = c("sum","max","first"),
    drop_all_zero_genes      = TRUE,
    stop_if_many_decimals    = TRUE,
    decimal_tol              = 1e-8,
    decimal_rate_threshold   = 0.001,
    allow_duplicate_sample_names = FALSE,
    duplicate_suffix_format  = "{name}.{group}",
    write_xlsx_preview       = TRUE,
    out_prefix               = ""
){
  duplicate_symbol_rule <- match.arg(duplicate_symbol_rule)
  
  ## ---------- helpers ----------
  pick_best <- function(base, .dir = dir){
    cand <- list.files(.dir, pattern = paste0("^", base, "\\.(xlsx|csv|xls)$"),
                       ignore.case = TRUE, full.names = TRUE)
    if (!length(cand)) stop("找不到文件：", base, ".xlsx/csv/xls 于 ", normalizePath(.dir))
    ord <- order(match(tolower(tools::file_ext(cand)), c("xlsx","csv","xls")))
    cand[ord][1]
  }
  find_col <- function(nms, keywords_exact = character(), keywords_regex = character()){
    nms_low <- tolower(nms)
    for (k in tolower(keywords_exact)) {
      hit <- which(nms_low == k); if (length(hit)) return(nms[hit[1]])
    }
    for (rgx in keywords_regex) {
      hit <- which(str_detect(nms_low, rgx)); if (length(hit)) return(nms[hit[1]])
    }
    NA_character_
  }
  is_annotation_name <- function(colname){
    pat <- "(^len$|length|gene_length|biotype|type$|class$|category$|chr|chromosome|start$|end$|strand$|gc$|desc|description|tx|transcript|exon|fpkm|tpm|rpkm)"
    str_detect(tolower(colname), pat)
  }
  
  ## 智能读取平面文本（csv/tsv/分号），自动选择“最像”的解析结果
  read_table_any <- function(path, guess_max = 100000){
    try_readers <- list(
      function() readr::read_csv(path, show_col_types = FALSE, guess_max = guess_max, progress = FALSE),
      function() readr::read_tsv(path, show_col_types = FALSE, guess_max = guess_max, progress = FALSE),
      function() readr::read_delim(path, delim = ";", show_col_types = FALSE,
                                   guess_max = guess_max, progress = FALSE,
                                   locale = locale(decimal_mark = ",", grouping_mark = "."))
    )
    cand <- list()
    for (rd in try_readers) {
      df <- tryCatch(rd(), error = function(e) NULL)
      if (is.null(df)) next
      df <- as_tibble(df)
      # 评分：列数越多越好；数值潜力越高越好
      ncol_score <- ncol(df)
      # 清理千分位/空格后评估数值潜力
      num_potential <- mean(vapply(df, function(v){
        if (is.numeric(v)) return(TRUE)
        if (!is.character(v)) return(FALSE)
        vv <- str_replace_all(v, "[ ,]", "")
        suppressWarnings(!all(is.na(as.numeric(vv))))
      }, logical(1)))
      cand[[length(cand)+1]] <- list(df=df, score = ncol_score + num_potential)
    }
    if (!length(cand)) stop("无法读取文本表：", basename(path))
    cand[[ which.max(vapply(cand, `[[`, numeric(1), "score")) ]]$df
  }
  
  ## 读取并清洗 counts（智能分隔 + 数字清理）
  read_counts_generic <- function(path){
    ext <- tolower(tools::file_ext(path))
    df <- switch(ext,
                 "xlsx" = readxl::read_excel(path, .name_repair = "minimal"),
                 "xls"  = readxl::read_excel(path, .name_repair = "minimal"),
                 "csv"  = read_table_any(path),
                 stop("不支持的扩展名：", ext))
    df <- as_tibble(df); names(df) <- trimws(names(df))
    
    id_col  <- find_col(names(df),
                        keywords_exact=c("#id","gene_id","id"),
                        keywords_regex=c("ensembl","^ensg","gene[_\\.]?id"))
    sym_col <- find_col(names(df),
                        keywords_exact=c("symbol","gene_symbol","gene","gene_name","hgnc_symbol"),
                        keywords_regex=c("external[_\\.]?gene[_\\.]?name","symbol"))
    
    if (is.na(sym_col) && drop_rows_without_symbol)
      message("警告：未找到 Symbol；将丢弃无 Symbol 的行（", basename(path), "）。")
    if (is.na(sym_col) && !drop_rows_without_symbol && is.na(id_col))
      stop("既无 Symbol 也无 #ID：", basename(path))
    
    anno_cols <- names(df)[vapply(names(df), is_annotation_name, logical(1))]
    anno_cols <- union(anno_cols, c(sym_col, id_col))
    sample_cols <- setdiff(names(df), anno_cols)
    
    ## 清理字符数字：去空格/千分位逗号；支持科学计数；尽最大努力转 numeric
    keep_sample <- logical(length(sample_cols))
    for (i in seq_along(sample_cols)) {
      colv <- df[[ sample_cols[i] ]]
      if (is.numeric(colv)) {
        keep_sample[i] <- TRUE
      } else {
        if (is.character(colv)) {
          colv <- str_replace_all(colv, "[\\s,]", "")  # 去空格/逗号（千分位）
        }
        xi <- suppressWarnings(as.numeric(as.character(colv)))
        keep_sample[i] <- !all(is.na(xi))
        if (keep_sample[i]) df[[ sample_cols[i] ]] <- xi
      }
    }
    sample_cols <- sample_cols[keep_sample]
    if (!length(sample_cols)) stop("未检测到数值型样本列：", basename(path))
    
    out <- tibble(Symbol = if (!is.na(sym_col)) as.character(df[[sym_col]]) else NA_character_)
    if (!is.na(id_col)) out <- mutate(out, `#ID` = as.character(df[[id_col]]))
    out <- bind_cols(out, df[, sample_cols, drop = FALSE])
    
    if (drop_rows_without_symbol) {
      out <- filter(out, !is.na(Symbol) & Symbol != "")
    } else {
      if (is.na(sym_col)) out <- mutate(out, Symbol = `#ID`)
      else out <- mutate(out, Symbol = if_else(is.na(Symbol) | Symbol=="",
                                               if ("#ID"%in%names(out)) `#ID` else Symbol, Symbol))
      out <- filter(out, !is.na(Symbol) & Symbol != "")
    }
    
    agg_fun <- switch(duplicate_symbol_rule,
                      "sum"=~sum(replace_na(.x,0)),
                      "max"=~max(replace_na(.x,0)),
                      "first"=~.x[which.max(!is.na(.x))]
    )
    out |> group_by(Symbol) |> summarise(across(all_of(sample_cols), agg_fun), .groups="drop")
  }
  
  ## ---------- 解析输入 ----------
  if (is.null(files)) {
    if (!is.null(bases)) {
      files <- vapply(bases, pick_best, character(1))
    } else {
      cand <- list.files(dir, pattern="^All_gene_counts.*\\.(xlsx|csv|xls)$",
                         ignore.case=TRUE, full.names=TRUE)
      if (!length(cand)) stop("未提供 files/bases，且目录下未找到 All_gene_counts*.xlsx/csv/xls")
      base <- tools::file_path_sans_ext(basename(cand))
      files <- tapply(cand, base, function(v){
        ord <- order(match(tolower(tools::file_ext(v)), c("xlsx","csv","xls"))); v[ord][1]
      }) |> unname()
    }
  }
  files <- normalizePath(files, mustWork = TRUE)
  message("使用文件：\n", paste(sprintf("  - %s", basename(files)), collapse="\n"))
  
  lst <- lapply(files, read_counts_generic)
  sample_lists <- lapply(lst, function(x) setdiff(names(x), "Symbol"))
  
  dups <- Reduce(intersect, sample_lists, init = sample_lists[[1]])
  if (length(dups) && !allow_duplicate_sample_names) {
    stop("发现跨文件重复样本名：\n  - ", paste(dups, collapse=", "),
         "\n请在源文件改名，或设 allow_duplicate_sample_names=TRUE（自动追加后缀）。")
  }
  if (length(dups) && allow_duplicate_sample_names) {
    if (is.null(group_labels)) group_labels <- paste0("batch", seq_along(files))
    for (i in seq_along(lst)) {
      smp <- setdiff(names(lst[[i]]), "Symbol")
      hit <- intersect(smp, dups)
      if (length(hit)) {
        newn <- vapply(hit, \(nm) str_replace_all(duplicate_suffix_format,
                                                  c("\\{name\\}"=nm, "\\{group\\}"=group_labels[i])), "")
        names(lst[[i]])[match(hit, names(lst[[i]]))] <- newn
      }
    }
    sample_lists <- lapply(lst, function(x) setdiff(names(x), "Symbol"))
  }
  
  merged <- purrr::reduce(
    lst,
    function(x, y) dplyr::full_join(x, y, by = "Symbol")
  )
  merged <- dplyr::mutate(merged, dplyr::across(-Symbol, ~ tidyr::replace_na(.x, 0)))
  
  if (drop_all_zero_genes) merged <- merged |> filter(rowSums(across(-Symbol)) > 0)
  
  num_mat <- as.matrix(merged[, -1, drop=FALSE])
  nonint_rate <- mean(abs(num_mat - round(num_mat)) > decimal_tol)
  if (stop_if_many_decimals && nonint_rate > decimal_rate_threshold) {
    stop(sprintf("检测到 %.2f%% 条目为非整数（疑似 FPKM/TPM）。", 100*nonint_rate))
  }
  
  if (is.null(group_labels)) group_labels <- paste0("batch", seq_along(files))
  if (length(group_labels) != length(files))
    stop("group_labels 长度需要与文件数一致。")
  
  coldata <- tibble(
    sample      = unlist(sample_lists, use.names = FALSE),
    batch       = rep(group_labels, lengths(sample_lists)),
    source_file = rep(basename(files), lengths(sample_lists)),
    condition   = NA_character_
  )
  readr::write_csv(coldata, file.path(dir, paste0(out_prefix, "coldata_group.csv")))
  readr::write_csv(merged,   file.path(dir, paste0(out_prefix, "counts_merged_for_deseq2.csv")))
  
  counts_mat <- merged |> column_to_rownames("Symbol") |> as.matrix()
  if (any(abs(counts_mat - round(counts_mat)) > decimal_tol)) counts_mat <- round(counts_mat)
  storage.mode(counts_mat) <- "integer"
  
  if (!setequal(colnames(counts_mat), coldata$sample))
    stop("coldata_group.csv 的 sample 集合与 counts 列名不一致。")
  counts_mat <- counts_mat[, coldata$sample, drop=FALSE]
  saveRDS(counts_mat, file.path(dir, paste0(out_prefix, "counts_merged_for_deseq2.rds")))
  
  if (write_xlsx_preview) {
    writexl::write_xlsx(list(
      counts_merged = as.data.frame(merged),
      coldata_group = as.data.frame(coldata)
    ), path = file.path(dir, paste0(out_prefix, "counts_and_group.xlsx")))
  }
  
  cat(sprintf("
完成：
  - 文件数：%d
  - 基因数：%d
  - 样本数：%d
  - 输出：%s{counts_merged_for_deseq2.csv,.rds}, %scoldata_group.csv%s
",
              length(files), nrow(merged), ncol(merged)-1,
              out_prefix, out_prefix, if (write_xlsx_preview) paste0(", ", out_prefix, "counts_and_group.xlsx") else ""
  ))
  
  invisible(list(counts = counts_mat, coldata = coldata, merged = merged))
}

score_gene_sets_bulk <- function(
    expr,                      # matrix: genes x samples
    meta,                      # data.frame: sample_id, group, optional batch
    gene_sets,                 # named list of character vectors
    method = c("ssgsea", "singscore", "mean_z"),
    min_genes = 10,
    # ssGSEA options
    ssgsea_alpha = 0.25,
    ssgsea_normalize = FALSE,
    # rescaling
    do_minmax01 = TRUE,
    minmax_by = c("all", "batch"),
    batch_col = "batch",
    verbose = FALSE
) {
  method <- match.arg(method)
  minmax_by <- match.arg(minmax_by)
  
  stopifnot(is.matrix(expr) || is.data.frame(expr))
  expr <- as.matrix(expr)
  storage.mode(expr) <- "double"
  
  if (is.null(rownames(expr)) || is.null(colnames(expr))) {
    stop("expr must have rownames (genes) and colnames (sample IDs).")
  }
  if (!all(c("sample_id", "group") %in% colnames(meta))) {
    stop("meta must contain columns: sample_id and group.")
  }
  if (!all(meta$sample_id %in% colnames(expr))) {
    missing <- setdiff(meta$sample_id, colnames(expr))
    stop("These meta$sample_id are not found in expr colnames: ", paste(missing, collapse = ", "))
  }
  
  # Reorder columns to match meta
  expr <- expr[, meta$sample_id, drop = FALSE]
  
  # Remove duplicated genes if any
  if (any(duplicated(rownames(expr)))) {
    expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
  }
  
  # Intersect genes with gene sets
  gene_sets2 <- lapply(gene_sets, function(gs) intersect(unique(gs), rownames(expr)))
  gene_sets2 <- gene_sets2[lengths(gene_sets2) >= min_genes]
  if (length(gene_sets2) == 0) stop("No gene sets left after intersecting with expr (min_genes too strict?).")
  
  score_mat <- NULL
  
  if (method == "ssgsea") {
    # Check if using new GSVA API (version >= 1.50)
    has_new_api <- "ssgseaParam" %in% getNamespaceExports("GSVA")
    
    if (has_new_api) {
      # --- New API (GSVA >= 1.50) ---
      # 直接构造参数对象并调用，不使用动态参数过滤，防止S4分发错误
      tryCatch({
        # 1. Construct parameter object
        param <- GSVA::ssgseaParam(
          exprData = expr,
          geneSets = gene_sets2,
          minSize = min_genes,
          maxSize = Inf,
          alpha = ssgsea_alpha,
          normalize = ssgsea_normalize
        )
        
        # 2. Run GSVA
        es <- GSVA::gsva(param, verbose = verbose)
        
      }, error = function(e) {
        stop("GSVA (New API) failed: ", e$message)
      })
      
      # Extract matrix from possible return types (SummarizedExperiment or matrix)
      if (inherits(es, "SummarizedExperiment")) {
        an <- SummarizedExperiment::assayNames(es)
        if ("es" %in% an) {
          score_mat <- SummarizedExperiment::assay(es, "es")
        } else {
          score_mat <- SummarizedExperiment::assay(es)
        }
      } else {
        score_mat <- as.matrix(es)
      }
      
    } else {
      # --- Old API (GSVA < 1.50) ---
      tryCatch({
        es <- GSVA::gsva(
          expr = expr,
          gset.idx.list = gene_sets2,
          method = "ssgsea",
          ssgsea.norm = ssgsea_normalize,
          min.sz = min_genes,
          verbose = verbose
        )
        score_mat <- as.matrix(es)
      }, error = function(e) {
        stop("GSVA (Old API) failed: ", e$message)
      })
    }
  }
  
  if (method == "singscore") {
    if (!requireNamespace("singscore", quietly = TRUE)) {
      stop("Package 'singscore' is required. Install via BiocManager::install('singscore').")
    }
    # Singscore usually requires rank data, but simpleScore handles it internally
    out_list <- lapply(names(gene_sets2), function(nm) {
      gs <- gene_sets2[[nm]]
      # simpleScore expects the full expression matrix and a gene set
      sc <- singscore::simpleScore(expr, upSet = gs, centerScore = FALSE)
      sc$TotalScore
    })
    score_mat <- do.call(rbind, out_list)
    rownames(score_mat) <- names(gene_sets2)
    colnames(score_mat) <- colnames(expr)
  }
  
  if (method == "mean_z") {
    z <- t(scale(t(expr)))
    out_list <- lapply(names(gene_sets2), function(nm) {
      gs <- gene_sets2[[nm]]
      colMeans(z[gs, , drop = FALSE], na.rm = TRUE)
    })
    score_mat <- do.call(rbind, out_list)
    rownames(score_mat) <- names(gene_sets2)
    colnames(score_mat) <- colnames(expr)
  }
  
  score_df <- as.data.frame(t(score_mat))
  score_df$sample_id <- rownames(score_df)
  score_df <- dplyr::left_join(meta, score_df, by = "sample_id")
  
  minmax01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (!is.finite(r[1]) || !is.finite(r[2]) || diff(r) == 0) return(rep(0.5, length(x)))
    (x - r[1]) / diff(r)
  }
  
  if (isTRUE(do_minmax01)) {
    gs_names <- names(gene_sets2)
    if (minmax_by == "all") {
      for (nm in gs_names) {
        score_df[[paste0(nm, "_01")]] <- minmax01(score_df[[nm]])
      }
    } else {
      if (!(batch_col %in% colnames(score_df))) {
        stop("minmax_by='batch' requires meta to contain column: ", batch_col)
      }
      score_df <- score_df %>%
        dplyr::group_by(.data[[batch_col]]) %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(gs_names), ~ minmax01(.x), .names = "{.col}_01")) %>%
        dplyr::ungroup()
    }
  }
  
  list(score_df = score_df, used_gene_sets = gene_sets2, score_mat = score_mat)
}

library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggExtra)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggExtra)
library(ggrepel)

# ==============================================================================
# 通用绘图函数 (支持颜色+形状双重分组)
# ==============================================================================
plot_nature_scatter<-  function(data, 
                                x_col, 
                                y_col, 
                                color_col,        
                                shape_col = NULL, 
                                label_col = NULL,
                                title = NULL,
                                subtitle = NULL,
                                x_lab = NULL, 
                                y_lab = NULL,
                                show_marginal = TRUE, 
                                label_size = 3) {
  
  if (is.null(x_lab)) x_lab <- x_col
  if (is.null(y_lab)) y_lab <- y_col
  
  # 1. 全局映射 (只含颜色)
  p <- ggplot(data, aes(x = .data[[x_col]], 
                        y = .data[[y_col]], 
                        color = .data[[color_col]], 
                        fill = .data[[color_col]]))
  
  # 2. 图层添加
  if (!is.null(shape_col)) {
    # 【关键修改】增加更多形状代码，支持最多 10 种分组
    # 16=圆, 17=三角, 15=方, 18=菱, 3=十字, 4=叉, 8=米字, 9=菱形叉, 10=十字圆, 12=方块十
    my_shapes <- c(16, 17, 15, 18, 3, 4, 8, 9, 10, 12)
    
    p <- p + geom_point(aes(shape = .data[[shape_col]]), size = 4, alpha = 0.8, stroke = 1) + 
      scale_shape_manual(values = my_shapes) 
  } else {
    p <- p + geom_point(size = 4, alpha = 0.8, shape = 21, stroke = 0.5, color = "white")
  }
  
  p <- p + stat_ellipse(aes(group = .data[[color_col]]), 
                        geom = "polygon", level = 0.68, alpha = 0.15, segments = 51, show.legend = FALSE) +
    geom_smooth(aes(group = .data[[color_col]]), 
                method = "lm", se = FALSE, linetype = "dashed", size = 0.8, alpha = 0.5) +
    stat_cor(aes(group = .data[[color_col]]), 
             method = "pearson", 
             label.x.npc = "left", label.y.npc = "top", 
             size = 4, show.legend = FALSE) +
    scale_color_npg() + 
    scale_fill_npg() +
    theme_classic(base_size = 16) +
    labs(x = x_lab, y = y_lab, color = "Group", fill = "Group", shape = "Cell Line") + 
    theme(
      axis.line = element_line(linewidth = 0.8, color = "black"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", size = 14),
      legend.position = c(0.85, 0.25), 
      legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
      legend.box = "vertical",
      legend.title = element_text(face="bold", size=10),
      legend.text = element_text(size=9)
    )
  
  # 3. 标签
  if (!is.null(label_col)) {
    p <- p + geom_text_repel(
      aes(label = .data[[label_col]]),
      size = label_size,
      color = "black",
      min.segment.length = 0,
      box.padding = 0.5,
      max.overlaps = Inf,
      show.legend = FALSE
    )
  }
  
  # 4. 边缘图
  if (show_marginal) {
    p_final <- ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE, alpha = 0.6, size = 5)
  } else {
    p_final <- p
  }
  
  # 5. 标题
  if (!is.null(title) || !is.null(subtitle)) {
    full_text <- title
    if(!is.null(subtitle)) full_text <- paste0(full_text, "\n", subtitle)
    p_final <- annotate_figure(p_final, top = text_grob(full_text, face = "bold", size = 16, lineheight = 1.2))
  }
  
  return(p_final)
}



# --- 1) Read a GMT file into a named list ---
read_gmt_as_list <- function(gmt_file, unique_genes = TRUE) {
  lines <- readLines(gmt_file, warn = FALSE)
  gs <- vector("list", length(lines))
  nm <- character(length(lines))
  
  for (i in seq_along(lines)) {
    parts <- strsplit(lines[[i]], "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) next
    nm[i] <- parts[1]
    genes <- parts[-c(1, 2)]
    genes <- genes[genes != ""]
    if (unique_genes) genes <- unique(genes)
    gs[[i]] <- genes
  }
  
  names(gs) <- nm
  gs <- gs[names(gs) != ""]
  gs
}

# --- 2) Pick the latest GMT file by regex pattern (fallback: last in sorted list) ---
pick_latest_gmt <- function(gmt_dir, pattern) {
  files <- list.files(gmt_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No GMT files matched pattern: ", pattern)
  files <- sort(files)
  files[[length(files)]]
}

# --- 3) Extract gene sets by exact names or keyword patterns ---
extract_gene_sets <- function(gs_list, keys, exact = TRUE, ignore_case = TRUE) {
  gs_names <- names(gs_list)
  
  pick_one <- function(key) {
    if (exact) {
      hit <- which(gs_names == key)
    } else {
      hit <- which(grepl(key, gs_names, ignore.case = ignore_case))
    }
    if (length(hit) == 0) stop("Gene set not found: ", key)
    if (length(hit) > 1) {
      warning("Multiple gene sets matched '", key, "'. Using the first: ", gs_names[hit[1]])
    }
    gs_names[hit[1]]
  }
  
  picked <- vapply(keys, pick_one, character(1))
  out <- gs_list[picked]
  names(out) <- picked
  out
}

# --- 4) One-stop wrapper: make list(KRAS=..., EMT=...) from your GMT directory ---
make_gene_sets <- function(
    gmt_dir,
    hallmark_pattern = "^h\\.all\\..*\\.gmt$",
    kras_set_name = "HALLMARK_KRAS_SIGNALING_UP",
    emt_set_name  = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    exact = TRUE,
    force_upper = TRUE
) {
  gmt_file <- pick_latest_gmt(gmt_dir, hallmark_pattern)
  gs_all <- read_gmt_as_list(gmt_file)
  
  picked <- extract_gene_sets(gs_all, keys = c(kras_set_name, emt_set_name), exact = exact)
  
  kras_genes <- picked[[kras_set_name]]
  emt_genes  <- picked[[emt_set_name]]
  
  if (force_upper) {
    kras_genes <- unique(toupper(kras_genes))
    emt_genes  <- unique(toupper(emt_genes))
  }
  
  list(KRAS = kras_genes, EMT = emt_genes)
}



# ===========================
# Modular RNA-seq pipeline (cts + cells)
# ===========================
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(ggprism)
  library(ggrepel)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(pROC)
  library(grid)
})

# ---- Science-style theme & palette ----
library(ggplot2)

theme_science <- function(base_size = 12, base_family = "Helvetica", show_ygrid = TRUE){
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # 文字
      plot.title   = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title   = element_text(face = "bold"),
      axis.text    = element_text(colour = "black"),
      # 轴与刻度
      axis.line    = element_line(colour = "black", linewidth = 0.6),
      axis.ticks   = element_line(colour = "black", linewidth = 0.6),
      axis.ticks.length = grid::unit(2.5, "pt"),
      # 网格：仅保留浅色 y 方向主网格（可选）
      panel.grid.major.y = if (show_ygrid) element_line(colour = "grey90", linewidth = 0.4) else element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      # 背景与边框
      panel.border = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      # 图例与分面
      legend.position = "right",
      legend.title    = element_text(face = "bold"),
      legend.background = element_blank(),
      legend.key      = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
      strip.text      = element_text(face = "bold")
    )
}

# 色盲友好 Okabe–Ito 调色板
pal_science <- function(n){
  cols <- c("#000000","#E69F00","#56B4E9","#009E73",
            "#F0E442","#0072B2","#D55E00","#CC79A7","#999999")
  cols[seq_len(min(n, length(cols)))]
}

# 与你现有代码兼容：把 prism_theme 指向新的主题
prism_theme <- theme_science(base_size = 12)

ensure_integer_counts <- function(
    cts,
    gene_id_col   = NULL,                   # 手动指定基因ID列，如 "#ID" / "SYMBOL" / "ENSEMBL"
    prefer        = c("auto","#ID","ENSEMBL","SYMBOL"),
    dedup_genes   = c("sum","mean","max","first","none"),
    drop_cols     = NULL,
    fill_na       = 0,
    round_to_int  = TRUE,
    allow_negative= FALSE,
    verbose       = TRUE
){
  prefer      <- match.arg(prefer)
  dedup_genes <- match.arg(dedup_genes)
  
  df <- as.data.frame(cts, check.names = FALSE)
  
  # 1) 选择基因ID列：若未指定，则自动选择“唯一值最多”的那列（tie-break: #ID > ENSEMBL > SYMBOL）
  if (is.null(gene_id_col)) {
    cand_map <- list(
      `#ID`     = c("#ID","Id","ID","Geneid","GeneID"),
      ENSEMBL   = c("ENSEMBL","Ensembl","ENSG","ENSEMBL_ID"),
      SYMBOL    = c("SYMBOL","Symbol","Gene","gene","GeneSymbol","GeneName")
    )
    # 展开候选并找出在数据中的列
    hits <- lapply(cand_map, function(v) intersect(v, names(df)))
    # 评估唯一值数量
    uniq_n <- sapply(names(hits), function(k){
      cols <- hits[[k]]
      if (length(cols)==0) return(-1L)
      u <- max(sapply(cols, function(cn) length(unique(df[[cn]]))), na.rm = TRUE)
      as.integer(u)
    })
    # 决策
    if (prefer == "auto") {
      # 选 uniq 值最多者；同分按 #ID > ENSEMBL > SYMBOL
      order_pref <- c("#ID","ENSEMBL","SYMBOL")
      mx <- max(uniq_n)
      if (mx < 0) stop("ensure_integer_counts: 未找到可用的基因ID列，请提供 gene_id_col。")
      winners <- names(uniq_n)[uniq_n == mx]
      pick_family <- order_pref[match(order_pref, winners, nomatch = 0) > 0][1]
    } else {
      pick_family <- prefer
      if (uniq_n[pick_family] < 0) stop("ensure_integer_counts: 未找到偏好列 '", prefer, "'，请检查列名或指定 gene_id_col。")
    }
    gene_id_col <- hits[[pick_family]][1]
    if (verbose) message(sprintf("[ensure_integer_counts] 使用基因ID列: %s", gene_id_col))
  }
  
  if (!gene_id_col %in% names(df))
    stop("ensure_integer_counts: 指定 gene_id_col='", gene_id_col, "' 不在数据里。")
  
  ids <- as.character(df[[gene_id_col]])
  df[[gene_id_col]] <- NULL
  
  # 2) 丢弃指定列
  if (!is.null(drop_cols)) {
    keep <- setdiff(names(df), drop_cols)
    df <- df[, keep, drop = FALSE]
  }
  
  # 3) 仅保留可转数值的列
  num <- suppressWarnings(lapply(df, function(x) as.numeric(as.character(x))))
  is_use <- vapply(num, function(v) any(!is.na(v)), logical(1))
  num <- as.data.frame(num[is_use], check.names = FALSE)
  if (ncol(num) == 0) stop("ensure_integer_counts: 没有可用的数值列。")
  
  # 4) 填 NA
  if (anyNA(num)) {
    if (verbose) message(sprintf("[ensure_integer_counts] 检测到 NA，已用 %s 填充。", as.character(fill_na)))
    num[is.na(num)] <- fill_na
  }
  M <- as.matrix(num)
  
  # 5) 负值检查
  if (!allow_negative && any(M < 0, na.rm = TRUE)) {
    bad <- which(M < 0, arr.ind = TRUE)
    ex  <- bad[seq_len(min(5, nrow(bad))), , drop = FALSE]
    exs <- apply(ex, 1, function(ix) sprintf("%s / %s = %s", ids[ix[1]], colnames(M)[ix[2]], M[ix[1], ix[2]]))
    stop(sprintf("counts 矩阵包含负值（示例：%s）。", paste(exs, collapse="; ")))
  }
  
  # 6) 如需合并重复基因
  if (anyDuplicated(ids)) {
    if (dedup_genes == "none") {
      dupn <- sum(duplicated(ids))
      stop("ensure_integer_counts: gene_id_col='", gene_id_col, "' 存在重复（", dupn, " 条）。设置 dedup_genes='sum/mean/max/first' 之一进行合并，或改用唯一ID列（如 #ID/ENSEMBL）。")
    }
    if (verbose) message(sprintf("[ensure_integer_counts] 发现重复基因ID（n=%d），按 '%s' 合并。", sum(duplicated(ids)), dedup_genes))
    if (dedup_genes == "sum") {
      M <- rowsum(M, group = ids, reorder = FALSE)
      ids <- rownames(M)
    } else {
      split_idx <- split(seq_along(ids), ids)
      combine_one <- switch(
        dedup_genes,
        mean = function(X) matrix(colMeans(X), nrow=1),
        max  = function(X) matrix(apply(X,2,max), nrow=1),
        first= function(X) matrix(X[1,,drop=FALSE], nrow=1)
      )
      out <- lapply(split_idx, function(ii){
        X <- M[ii,,drop=FALSE]
        if (dedup_genes=="sum") matrix(colSums(X), nrow=1) else combine_one(X)
      })
      M <- do.call(rbind, out)
      ids <- names(split_idx)
      rownames(M) <- ids
    }
  } else {
    rownames(M) <- ids
  }
  
  # 7) 取整为整数计数
  if (round_to_int) {
    if (any(abs(M - round(M)) > 1e-8, na.rm = TRUE) && verbose)
      message("[ensure_integer_counts] 检测到非整数计数，已四舍五入。")
    M <- round(M)
  }
  storage.mode(M) <- "integer"
  return(M)
}

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(ggpubr) # 必须加载

plot_gene_response <- function(
    expression_matrix,       # 表达矩阵 (VST 或 Counts)
    cell_info,               # 分组信息
    genes_to_plot,           # 基因列表
    work_dir, 
    out_prefix = "Gene_Expression",
    do_log2 = FALSE,         # 是否做 log2 转换
    
    # --- 新增参数 ---
    test_method = "t.test",  # 可选: "t.test" 或 "wilcox.test"
    # ----------------
    
    colors = c("Sensitive" = "#3B82F6", "Resistant" = "#EF4444")
) {
  
  # 1. 数据准备
  valid_genes <- intersect(genes_to_plot, rownames(expression_matrix))
  if (length(valid_genes) == 0) stop("所选基因在矩阵中不存在。")
  if (length(valid_genes) < length(genes_to_plot)) warning("部分基因未找到: ", paste(setdiff(genes_to_plot, rownames(expression_matrix)), collapse=", "))
  
  sub_mat <- expression_matrix[valid_genes, , drop = FALSE]
  
  if (do_log2) {
    sub_mat <- log2(sub_mat + 1)
    ylab_text <- "Log2(Counts + 1)"
  } else {
    ylab_text <- "Expression Level (Normalized/VST)"
  }
  
  df_values <- as.data.frame(t(sub_mat))
  df_values$cell <- rownames(df_values)
  
  cell_info$cell <- as.character(cell_info$cell)
  if (!"Group" %in% colnames(cell_info)) stop("cell_info 缺少 'Group' 列")
  
  df_plot <- left_join(df_values, cell_info[, c("cell", "Group")], by = "cell")
  df_plot <- na.omit(df_plot)
  
  df_long <- df_plot %>%
    pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Expression")
  
  # 2. 绘制 Boxplot
  p_box <- ggplot(df_long, aes(x = Group, y = Expression)) +
    geom_boxplot(aes(fill = Group), alpha = 0.6, outlier.shape = NA) +
    geom_jitter(aes(color = Group), width = 0.2, size = 2, alpha = 0.8) +
    geom_text_repel(aes(label = cell), size = 3, max.overlaps = 15, box.padding = 0.3) +
    facet_wrap(~Gene, scales = "free_y") +
    theme_bw(base_size = 14) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    
    # --- 关键修改：使用变量 test_method ---
    stat_compare_means(method = test_method, label = "p.signif", 
                       comparisons = list(c("Sensitive", "Resistant"))) +
    # ------------------------------------
  
  labs(title = paste0("Gene Expression (Test: ", test_method, ")"), 
       y = ylab_text, x = "Group")
  
  # 3. 绘制 Scatter (仅当2个基因时)
  p_scatter <- NULL
  if (length(valid_genes) == 2) {
    g1 <- valid_genes[1]; g2 <- valid_genes[2]
    p_scatter <- ggplot(df_plot, aes(x = .data[[g1]], y = .data[[g2]], color = Group)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text_repel(aes(label = cell), size = 3, max.overlaps = 20) +
      geom_smooth(method = "lm", se = FALSE, linetype = "dashed", alpha = 0.5) +
      theme_bw(base_size = 14) +
      scale_color_manual(values = colors) +
      labs(title = paste0("Scatter: ", g1, " vs ", g2), x = g1, y = g2)
  }
  
  # 4. 保存
  ggsave(file.path(work_dir, paste0(out_prefix, "_Boxplot.pdf")), p_box, width = 6 + length(valid_genes)*2, height = 6)
  if (!is.null(p_scatter)) ggsave(file.path(work_dir, paste0(out_prefix, "_Scatter.pdf")), p_scatter, width = 7, height = 6)
  
  message(paste("Plots saved to:", work_dir))
  return(list(p_box = p_box, p_scatter = p_scatter))
}

harmonize_sample_names <- function(
    cts,
    cells,
    cell_col      = NULL,                    # which column in 'cells' holds sample names
    group_col     = NULL,                    # optional group column name
    gene_col      = NULL,                    # which column in 'cts' holds gene IDs if rownames are empty
    fix_target    = c("cells","cts","none"), # rename which side to match the other
    drop_missing  = TRUE,                    # drop samples that appear only on one side
    reorder_cols  = TRUE,                    # order cts columns to follow 'cells' order
    dedup_genes   = c("sum","max","mean","first","none"),
    combine_dupe_samples = c("sum","mean","stop"),
    alias_table   = NULL,                    # optional data.frame(alias, canonical)
    allow_oneoff  = 0,                       # 0 = exact only; >0 = fuzzy distance for agrep/adist
    normalize_fun = NULL,                    # custom normalization; default removes case/space/_/-
    to_integer    = TRUE,                    # cast storage.mode(cts) <- "integer"
    verbose       = TRUE
){
  fix_target <- match.arg(fix_target)
  dedup_genes <- match.arg(dedup_genes)
  combine_dupe_samples <- match.arg(combine_dupe_samples)
  
  .msg <- function(...) if (verbose) message("[harmonize] ", sprintf(...))
  default_normalize <- function(x){
    x <- toupper(as.character(x))
    x <- gsub("\\s+", "", x)   # remove spaces
    x <- gsub("_", "", x)      # remove underscores
    x <- gsub("-", "", x)      # remove dashes
    x <- gsub("\\.", "", x)    # remove dots
    x <- gsub("\\W", "", x)    # remove non-word chars
    x
  }
  norm <- if (is.null(normalize_fun)) default_normalize else normalize_fun
  
  keep_unique_first <- function(x){
    ux <- !duplicated(x) | is.na(x)
    x[ux]
  }
  
  combine_cols <- function(mat, new_names, how = "sum"){
    d <- as.data.frame(mat, check.names = FALSE)
    names(d) <- new_names
    dup_tab <- table(names(d))
    if (any(dup_tab > 1L)) {
      dupe <- names(dup_tab)[dup_tab > 1L]
      for (nm in dupe) {
        idx <- which(names(d) == nm)
        if (how == "sum") {
          d[[nm]] <- rowSums(as.matrix(d[idx]), na.rm = TRUE)
        } else if (how == "mean") {
          d[[nm]] <- rowMeans(as.matrix(d[idx]), na.rm = TRUE)
        } else if (how == "stop") {
          stop("Duplicated samples after renaming: ", nm,
               " (use combine_dupe_samples = 'sum'/'mean').")
        }
        if (length(idx) > 1L) {
          # keep first, drop the rest
          drop_idx <- idx[-1L]
          d[drop_idx] <- NULL
        }
      }
    }
    as.matrix(d)
  }
  
  if (is.null(cell_col)) {
    cand <- intersect(c("Cell","cell","Sample","sample","cell_id","sample_id","SampleID","ID"), colnames(cells))
    if (length(cand) == 0L)
      stop("Cannot find a sample-name column in 'cells'. Please specify 'cell_col'.")
    cell_col <- cand[1L]
  }
  cells <- as.data.frame(cells, stringsAsFactors = FALSE)
  cells[[cell_col]] <- as.character(cells[[cell_col]])
  if (!is.null(group_col) && !group_col %in% names(cells)) {
    .msg("Specified group_col='%s' not found; ignoring.", group_col)
    group_col <- NULL
  }
  gene_candidates <- c("Symbol","SYMBOL","Gene","gene","GeneSymbol","GeneName","ID","#ID",
                       "Ensembl","ENSEMBL","ENSEMBL_ID","GENEID")
  if (!is.null(gene_col)) gene_candidates <- c(gene_col, gene_candidates)
  
  if (is.matrix(cts)) {
    mat <- cts
  } else {
    cts_df <- as.data.frame(cts, check.names = FALSE)
    # Detect gene column if rownames are empty or numeric
    if (is.null(rownames(cts_df)) || all(is.na(rownames(cts_df))) ||
        all(grepl("^[0-9]+$", rownames(cts_df)))) {
      hit <- intersect(gene_candidates, colnames(cts_df))
      if (length(hit) > 0L) {
        gene_col <- hit[1L]
        gene_vec <- cts_df[[gene_col]]
        cts_df[[gene_col]] <- NULL
        # Coerce the rest to numeric
        for (j in seq_along(cts_df)) {
          cts_df[[j]] <- suppressWarnings(as.numeric(cts_df[[j]]))
        }
        mat <- as.matrix(cts_df)
        rownames(mat) <- make.names(gene_vec, unique = FALSE)
      } else {
        stop("Cannot determine gene column in 'cts'. Provide 'gene_col'.")
      }
    } else {
      # Try to coerce all columns to numeric; keep as is if already matrix-like
      mat <- as.matrix(suppressWarnings(apply(cts_df, 2, as.numeric)))
      rownames(mat) <- rownames(cts_df)
    }
  }
  # Drop rows with empty/NA gene ids
  bad_gene <- is.na(rownames(mat)) | rownames(mat) == "" | rownames(mat) == "NA"
  if (any(bad_gene)) mat <- mat[!bad_gene, , drop = FALSE]
  
  # Deduplicate gene rows
  if (dedup_genes != "none") {
    .msg("Deduplicating genes by '%s' ...", dedup_genes)
    u <- unique(rownames(mat))
    if (length(u) < nrow(mat)) {
      if (dedup_genes == "sum") {
        mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
      } else if (dedup_genes == "max") {
        mat <- do.call(rbind, lapply(split(mat, rownames(mat)), function(m){
          m <- as.matrix(m); m[which.max(rowSums(m)), , drop = FALSE]
        }))
      } else if (dedup_genes == "mean") {
        mat <- do.call(rbind, lapply(split(mat, rownames(mat)), function(m){
          m <- as.matrix(m); matrix(colMeans(m), nrow = 1,
                                    dimnames = list(rownames(m)[1], colnames(m)))
        }))
      } else if (dedup_genes == "first") {
        mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
      }
    }
  }
  if (to_integer) storage.mode(mat) <- "integer"
  
  # ---------- 3) build normalized name maps ----------
  cts_names   <- colnames(mat)
  cells_names <- as.character(cells[[cell_col]])
  cts_norm    <- norm(cts_names)
  cells_norm  <- norm(cells_names)
  
  # Optional alias mapping (two cols: alias, canonical)
  alias_map <- NULL
  if (!is.null(alias_table)) {
    stopifnot(all(c("alias","canonical") %in% colnames(alias_table)))
    alias_map <- setNames(as.character(alias_table$canonical), norm(alias_table$alias))
  }
  
  # ---------- 4) match samples ----------
  # unique canonical map from cells norm -> cells canonical
  canon_by_norm <- split(cells_names, cells_norm)
  canon_by_norm <- vapply(canon_by_norm, function(v) v[1L], FUN.VALUE = character(1))
  
  matched <- rep(NA_character_, length(cts_names))
  reason  <- rep("no_match", length(cts_names))
  
  for (i in seq_along(cts_names)) {
    nn <- cts_norm[i]
    # exact norm match to cells
    if (nn %in% names(canon_by_norm)) {
      matched[i] <- canon_by_norm[[nn]]
      reason[i]  <- "cells_exact"
      next
    }
    # alias table
    if (!is.null(alias_map) && nn %in% names(alias_map)) {
      matched[i] <- alias_map[[nn]]
      reason[i]  <- "alias_map"
      next
    }
    # fuzzy (optional)
    if (allow_oneoff > 0) {
      hit <- agrep(nn, cells_norm, max.distance = allow_oneoff, value = FALSE)
      if (length(hit) == 1L) {
        matched[i] <- cells_names[hit]
        reason[i]  <- sprintf("fuzzy_%s", allow_oneoff)
      } else if (length(hit) > 1L) {
        reason[i]  <- sprintf("ambiguous_%s", allow_oneoff)
      }
    }
  }
  
  map <- data.frame(
    cts_old   = cts_names,
    cts_norm  = cts_norm,
    matched   = matched,
    reason    = reason,
    stringsAsFactors = FALSE
  )
  
  # ---------- 5) rename / align according to fix_target ----------
  issues <- character(0)
  
  if (fix_target == "cells") {
    new_names <- ifelse(!is.na(map$matched), map$matched, map$cts_old)
    # Combine duplicates if any
    mat2 <- combine_cols(mat, new_names, how = combine_dupe_samples)
    
    # Reorder and/or drop to match cells
    common <- intersect(colnames(mat2), cells_names)
    cts_only   <- setdiff(colnames(mat2), cells_names)
    cells_only <- setdiff(cells_names, colnames(mat2))
    
    if (drop_missing) {
      if (length(cts_only))   issues <- c(issues, sprintf("Dropped from cts: %s", paste(cts_only, collapse = ", ")))
      if (length(cells_only)) issues <- c(issues, sprintf("Dropped from cells: %s", paste(cells_only, collapse = ", ")))
      mat2   <- mat2[, common, drop = FALSE]
      cells2 <- cells[cells[[cell_col]] %in% common, , drop = FALSE]
    } else {
      cells2 <- cells
    }
    
    if (reorder_cols) {
      ord <- match(cells2[[cell_col]], colnames(mat2))
      ord <- ord[!is.na(ord)]
      if (length(ord)) mat2 <- mat2[, ord, drop = FALSE]
    }
    
  } else if (fix_target == "cts") {
    # Rename cells to cts (less common but supported)
    rename_map <- na.omit(map[, c("matched","cts_old")])
    if (nrow(rename_map)) {
      idx <- match(cells[[cell_col]], rename_map$matched)
      repl <- rename_map$cts_old[idx]
      repl[is.na(repl)] <- cells[[cell_col]][is.na(repl)]
      cells2 <- cells
      cells2[[cell_col]] <- repl
    } else {
      cells2 <- cells
    }
    
    common <- intersect(colnames(mat), cells2[[cell_col]])
    cts_only   <- setdiff(colnames(mat), cells2[[cell_col]])
    cells_only <- setdiff(cells2[[cell_col]], colnames(mat))
    
    if (drop_missing) {
      if (length(cts_only))   issues <- c(issues, sprintf("Dropped from cts: %s", paste(cts_only, collapse = ", ")))
      if (length(cells_only)) issues <- c(issues, sprintf("Dropped from cells: %s", paste(cells_only, collapse = ", ")))
      mat2   <- mat[, common, drop = FALSE]
      cells2 <- cells2[cells2[[cell_col]] %in% common, , drop = FALSE]
    } else {
      mat2 <- mat
    }
    
    if (reorder_cols) {
      ord <- match(cells2[[cell_col]], colnames(mat2))
      ord <- ord[!is.na(ord)]
      if (length(ord)) mat2 <- mat2[, ord, drop = FALSE]
    }
    
  } else { # "none"
    mat2 <- mat
    cells2 <- cells
  }
  
  dropped <- list(
    cts_only   = setdiff(colnames(mat), cells_names),
    cells_only = setdiff(cells_names, colnames(mat))
  )
  
  # Final sanity
  if (ncol(mat2) == 0L)
    stop("After harmonization there are 0 overlapping samples. Check your inputs or set drop_missing=FALSE.")
  
  invisible(list(
    cts = mat2,
    cells = cells2,
    map = map,
    dropped = dropped,
    issues = issues
  ))
}
# ================= (end) Enhanced harmonize_sample_names =================

build_coldata <- function(cells, sample_names, group_col = "Group", group_levels = NULL, batch_col = NULL) {
  # 1. 安全检查
  if (is.null(cells) || nrow(cells) == 0) stop("metadata (cells) is empty.")
  
  # 确保是 data.frame
  cells <- as.data.frame(cells)
  
  # 2. 智能寻找样本名列 (ID Column)
  # 我们需要找到 cells 中哪一列的内容，跟 sample_names (counts的列名) 重合最多
  found_col <- NULL
  max_overlap <- -1
  
  # 遍历每一列进行测试
  for (nm in colnames(cells)) {
    # 计算交集数量
    ov <- sum(as.character(cells[[nm]]) %in% sample_names)
    if (ov > max_overlap) {
      max_overlap <- ov
      found_col <- nm
    }
  }
  
  # 如果列里没找到，看看行名是否匹配？
  if (max_overlap == 0) {
    if (sum(rownames(cells) %in% sample_names) > 0) {
      cells$SampleID_Auto <- rownames(cells)
      found_col <- "SampleID_Auto"
    } else {
      stop("build_coldata Failed: Cannot find any column in 'cells' that matches the counts sample names.")
    }
  }
  
  # 3. 对齐与构建
  # 仅保留共有的样本
  common_samples <- intersect(as.character(cells[[found_col]]), sample_names)
  if (length(common_samples) == 0) stop("No common samples between metadata and counts.")
  
  # 按照 counts (sample_names) 的顺序或者交集顺序排列
  # 这里我们取交集，并以 metadata 中的该列为行名
  rownames(cells) <- as.character(cells[[found_col]])
  
  # 提取子集
  meta_sub <- cells[common_samples, , drop = FALSE]
  
  # 构建最终的 colData (data.frame)
  res <- data.frame(row.names = common_samples, stringsAsFactors = FALSE)
  
  # 填充分组 Group
  if (!is.null(group_col) && group_col %in% colnames(meta_sub)) {
    res$Group <- meta_sub[[group_col]]
    if (!is.null(group_levels)) {
      res$Group <- factor(res$Group, levels = group_levels)
    } else {
      res$Group <- as.factor(res$Group)
    }
  } else {
    stop("Group column '", group_col, "' not found in metadata.")
  }
  
  # 填充批次 Batch
  if (!is.null(batch_col) && batch_col %in% colnames(meta_sub)) {
    res$Batch <- as.factor(meta_sub[[batch_col]])
  }
  
  # [关键] 显式添加 Sample 列，防止后续报错
  res$Sample <- rownames(res)
  
  # 最后按照 sample_names 的顺序重排一下 (如果 counts 里有 metadata 里没有的，会在外层过滤)
  # 这里只返回共有的即可
  return(res)
}
# ========== 3) 基因过滤 ==========
filter_counts <- function(counts, 
                          min_count = 10, 
                          min_samples_prop = 0.10, 
                          remove_pattern = "^NewGene_") {
  
  # 1. (新增) 按名称模式过滤：移除匹配 remove_pattern 的基因
  # 如果不想过滤任何模式，可以将 remove_pattern 设为 NULL
  if (!is.null(remove_pattern) && remove_pattern != "") {
    is_bad_name <- grepl(remove_pattern, rownames(counts))
    if (any(is_bad_name)) {
      message(sprintf("Pre-filtering: Removing %d genes matching pattern '%s'", sum(is_bad_name), remove_pattern))
      counts <- counts[!is_bad_name, , drop = FALSE]
    }
  }
  
  # 2. (原有) 按表达量过滤：至少在 k 个样本中 count > min_count
  k <- max(1, ceiling(ncol(counts) * min_samples_prop))
  keep_expr <- rowSums(counts > min_count) >= k
  
  # 返回过滤后的 counts 和保留的索引逻辑向量（仅针对当前的行）
  list(counts = counts[keep_expr, , drop = FALSE], keep_idx = keep_expr)
}

# ========== 4) QC（总量、相关性、PCA、分布、离群）==========
qc_plots <- function(counts_filtered, colData, output_dir, prefix = "Project") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # —— 强制对齐：行名=样本名，顺序=counts 列顺序 ——
  colData <- as.data.frame(colData)
  if (is.null(rownames(colData)) || !all(rownames(colData) == colData$Sample)) {
    rownames(colData) <- colData$Sample
  }
  colData <- colData[colnames(counts_filtered), , drop = FALSE]
  
  # 1) 总量柱状图
  sample_sums <- colSums(counts_filtered)
  p_total <- ggplot(data.frame(Sample = names(sample_sums), TotalCounts = sample_sums),
                    aes(Sample, TotalCounts)) +
    geom_bar(stat = "identity") +
    prism_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste0(prefix, " - Total Counts per Sample"), x = "", y = "Total Counts")
  ggsave(file.path(output_dir, paste0(prefix, "_Total_Counts_per_Sample.png")),
         p_total, width = 10, height = 6, dpi = 300)
  
  # 2) VST
  dds0 <- DESeq2::DESeqDataSetFromMatrix(counts_filtered, colData, design = ~ 1)
  vs   <- DESeq2::vst(dds0, blind = TRUE)
  vsm  <- SummarizedExperiment::assay(vs)
  
  # 3) 相关性热图（确保注释行名与列名一致）
  cor_mat <- stats::cor(vsm)
  ann_col <- data.frame(Group = colData$Group)
  rownames(ann_col) <- rownames(colData)  # 必须与 colnames(cor_mat) 一致
  pheatmap::pheatmap(
    cor_mat,
    annotation_col = ann_col,
    main = paste0(prefix, " - Sample Correlation (VST)"),
    color = colorRampPalette(c("navy", "white", "firebrick"))(100),
    filename = file.path(output_dir, paste0(prefix, "_Sample_Correlation_Heatmap.png"))
  )
  
  # 4) PCA
  set.seed(2025)
  pca <- prcomp(t(vsm), scale. = TRUE)
  var1 <- round(100 * (pca$sdev[1]^2 / sum(pca$sdev^2)), 1)
  var2 <- round(100 * (pca$sdev[2]^2 / sum(pca$sdev^2)), 1)
  pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], colData[rownames(pca$x), , drop = FALSE])
  
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 3) +
    prism_theme +
    labs(title = paste0(prefix, " - PCA (VST)"),
         x = paste0("PC1 (", var1, "%)"),
         y = paste0("PC2 (", var2, "%)"))
  if ("Batch" %in% colnames(colData)) p_pca <- p_pca + aes(shape = Batch)
  ggsave(file.path(output_dir, paste0(prefix, "_PCA.png")), p_pca, width = 6, height = 6, dpi = 300)
  
  # 5) 表达分布（log2）
  log_counts <- log2(counts_filtered + 1)
  png(file.path(output_dir, paste0(prefix, "_Log2_Expression_Distribution.png")),
      width = 1000, height = 700, res = 130)
  boxplot(log_counts, las = 2, main = paste0(prefix, " - Log2 Expression per Sample"),
          col = "lightblue", border = "navy")
  dev.off()
  
  # 6) 离群（VST 后每样本标准差）
  sample_sd <- apply(vsm, 2, sd)
  outliers <- names(sample_sd)[sample_sd > mean(sample_sd) + 2 * stats::sd(sample_sd)]
  writeLines(paste("Potential outliers:",
                   ifelse(length(outliers) > 0, paste(outliers, collapse = ", "), "none")),
             con = file.path(output_dir, paste0(prefix, "_Outliers.txt")))
  
  list(vst_mat = vsm, dds_vst = vs, outliers = outliers,
       plots = list(total = p_total, pca = p_pca))
}

write_dds_summaries <- function(dds, output_dir, prefix = "Project"){
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # 保存 dds
  saveRDS(dds, file.path(output_dir, paste0(prefix, "_DESeq2_DDS.rds")))
  # 导出样本信息
  write.csv(as.data.frame(SummarizedExperiment::colData(dds)),
            file.path(output_dir, paste0(prefix, "_DESeq2_colData.csv")))
  # 导出 size factors
  sf <- data.frame(sample = colnames(dds), sizeFactor = DESeq2::sizeFactors(dds))
  write.csv(sf, file.path(output_dir, paste0(prefix, "_DESeq2_sizeFactors.csv")), row.names = FALSE)
  # 导出基因离散度
  disp <- data.frame(gene = rownames(dds), dispersion = DESeq2::dispersions(dds))
  write.csv(disp, file.path(output_dir, paste0(prefix, "_DESeq2_geneDispersions.csv")), row.names = FALSE)
}


# ========== 5) 差异分析（DESeq2）==========
run_deseq2 <- function(counts_filtered, colData, design_formula = ~ Group,
                       contrast = c("Group","Resistant","Sensitive"),
                       output_dir, prefix = "Project",
                       # Always computed:
                       lognorm_scale_factor = 1e4,
                       # Optional TPM/FPKM:
                       calc_tpm_fpkm = FALSE,
                       gene_length_bp = NULL) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # ----------------------------
  # 1) DESeq2 (unchanged)
  # ----------------------------
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered,
                                        colData   = colData,
                                        design    = design_formula)
  dds <- DESeq2::DESeq(dds)
  
  res <- DESeq2::results(dds, contrast = contrast)
  res <- as.data.frame(res)
  res$Symbol <- rownames(res)
  res <- res[order(res$padj), ]
  
  base_filename <- paste0(prefix, "_DESeq2_", contrast[2], "_vs_", contrast[3])
  out_csv  <- file.path(output_dir, paste0(base_filename, ".csv"))
  out_xlsx <- file.path(output_dir, paste0(base_filename, ".xlsx"))
  
  write.csv(res, out_csv, row.names = FALSE)
  writexl::write_xlsx(res, out_xlsx)
  
  # Prepare count matrix once
  counts_mat <- as.matrix(counts_filtered)
  storage.mode(counts_mat) <- "double"
  libsize <- colSums(counts_mat)
  
  # ----------------------------
  # 2) Always: CPM
  # CPM = counts / library_size * 1e6
  # ----------------------------
  cpm_mat <- t(t(counts_mat) / libsize * 1e6)
  cpm_mat <- as.data.frame(cpm_mat)
  
  cpm_save <- cpm_mat
  cpm_save$Symbol <- rownames(cpm_save)
  cpm_file <- file.path(output_dir, paste0(prefix, "_CPM_Matrix.xlsx"))
  writexl::write_xlsx(cpm_save, cpm_file)
  
  # ----------------------------
  # 3) Always: LogNorm (Seurat-style log1p(CP10K))
  # log1p( counts / library_size * scale_factor )
  # ----------------------------
  lognorm_mat <- log1p(t(t(counts_mat) / libsize * lognorm_scale_factor))
  lognorm_mat <- as.data.frame(lognorm_mat)
  
  lognorm_save <- lognorm_mat
  lognorm_save$Symbol <- rownames(lognorm_save)
  lognorm_file <- file.path(output_dir, paste0(prefix, "_LogNorm_Matrix_sf", lognorm_scale_factor, ".xlsx"))
  writexl::write_xlsx(lognorm_save, lognorm_file)
  
  # ----------------------------
  # 4) Always: VST (for PCA/heatmap/GSVA)
  # ----------------------------
  vst_obj <- DESeq2::vst(dds, blind = TRUE)
  vst_mat <- SummarizedExperiment::assay(vst_obj) %>% as.data.frame()
  
  vst_save <- vst_mat
  vst_save$Symbol <- rownames(vst_save)
  vst_file <- file.path(output_dir, paste0(prefix, "_VST_Matrix.xlsx"))
  writexl::write_xlsx(vst_save, vst_file)
  
  # ----------------------------
  # 5) Optional: TPM / FPKM
  # ----------------------------
  tpm_mat <- NULL
  fpkm_mat <- NULL
  
  if (isTRUE(calc_tpm_fpkm)) {
    if (is.null(gene_length_bp)) {
      stop("calc_tpm_fpkm=TRUE requires gene_length_bp (named numeric vector; bp).")
    }
    
    len <- gene_length_bp[rownames(counts_mat)]
    if (any(is.na(len))) {
      missing_n <- sum(is.na(len))
      stop("gene_length_bp missing for ", missing_n, " genes (match rownames(counts_filtered)).")
    }
    len_kb <- len / 1000
    
    # FPKM = (counts / gene_length_kb) / (library_size / 1e6)
    fpkm <- (counts_mat / len_kb) / (libsize / 1e6)
    
    # TPM: RPK then per-sample scaling to 1e6
    rpk <- counts_mat / len_kb
    tpm <- t(t(rpk) / colSums(rpk) * 1e6)
    
    fpkm_mat <- as.data.frame(fpkm)
    tpm_mat  <- as.data.frame(tpm)
    
    fpkm_save <- fpkm_mat; fpkm_save$Symbol <- rownames(fpkm_save)
    tpm_save  <- tpm_mat;  tpm_save$Symbol  <- rownames(tpm_save)
    
    fpkm_file <- file.path(output_dir, paste0(prefix, "_FPKM_Matrix.xlsx"))
    tpm_file  <- file.path(output_dir, paste0(prefix, "_TPM_Matrix.xlsx"))
    
    writexl::write_xlsx(fpkm_save, fpkm_file)
    writexl::write_xlsx(tpm_save,  tpm_file)
  }
  
  # ----------------------------
  # Return
  # ----------------------------
  list(
    dds = dds,
    res = res,
    vst_mat = vst_mat,
    lognorm_mat = lognorm_mat,
    cpm_mat = cpm_mat,
    tpm_mat = tpm_mat,
    fpkm_mat = fpkm_mat
  )
}



suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})
summarize_de_results <- function(res_df, 
                                 output_dir, 
                                 prefix = "ProjA", 
                                 padj_cutoff = 0.05, 
                                 log2fc_cutoff = 1,
                                 remove_pattern = "^NewGene_",
                                 filter_metric = "padj") { # 【新增参数】默认用 padj，可选 "pvalue"
  
  # 强制加载 dplyr 防止报错
  if (!requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
  require(ggplot2)
  
  # 1. 过滤 (Filtering)
  # 去掉无意义的基因，并去掉【选定指标】为 NA 的行
  clean_res <- res_df %>%
    dplyr::filter(!grepl(remove_pattern, Symbol)) %>%
    dplyr::filter(!is.na(.data[[filter_metric]])) # 【修改】动态过滤
  
  message(sprintf("Original rows: %d -> Cleaned rows: %d (Metric: %s)", 
                  nrow(res_df), nrow(clean_res), filter_metric))
  
  # 2. 分类 (Categorize)
  # Resistant vs Sensitive
  clean_res <- clean_res %>%
    dplyr::mutate(Category = dplyr::case_when(
      # 【修改】使用 .data[[filter_metric]] 动态引用列
      .data[[filter_metric]] < padj_cutoff & log2FoldChange > log2fc_cutoff ~ "Resistant_Marker",
      .data[[filter_metric]] < padj_cutoff & log2FoldChange < -log2fc_cutoff ~ "Sensitive_Marker",
      TRUE ~ "Not_Significant"
    ))
  
  # 3. 统计数量 (Statistics)
  stats_table <- clean_res %>%
    dplyr::group_by(Category) %>%
    dplyr::summarise(Count = dplyr::n())
  
  print("=== 差异基因统计 (Differential Expression Summary) ===")
  print(stats_table)
  
  # 4. 提取 Top 30 (Extract Top Genes)
  # 【修改】排序和选择列时使用动态指标
  top_resistant <- clean_res %>%
    dplyr::filter(Category == "Resistant_Marker") %>%
    dplyr::arrange(.data[[filter_metric]]) %>% # 按选定指标排序
    head(30) %>%
    dplyr::select(Symbol, log2FoldChange, dplyr::all_of(filter_metric)) # 导出对应指标列
  
  top_sensitive <- clean_res %>%
    dplyr::filter(Category == "Sensitive_Marker") %>%
    dplyr::arrange(.data[[filter_metric]]) %>% # 按选定指标排序
    head(30) %>%
    dplyr::select(Symbol, log2FoldChange, dplyr::all_of(filter_metric)) # 导出对应指标列
  
  # 5. 打印到控制台 (Print to Console)
  cat("\n>>> Top Resistant Markers (Up-regulated in Resistant):\n")
  if(nrow(top_resistant) > 0) print(as.data.frame(top_resistant)) else print("None found.")
  
  cat("\n>>> Top Sensitive Markers (Down-regulated in Resistant):\n")
  if(nrow(top_sensitive) > 0) print(as.data.frame(top_sensitive)) else print("None found.")
  
  # 6. 保存表格 (Save Tables)
  write.csv(clean_res, file.path(output_dir, paste0(prefix, "_Clean_DE_Results.csv")), row.names = FALSE)
  write.csv(top_resistant, file.path(output_dir, paste0(prefix, "_Top30_Resistant_Markers.csv")), row.names = FALSE)
  write.csv(top_sensitive, file.path(output_dir, paste0(prefix, "_Top30_Sensitive_Markers.csv")), row.names = FALSE)
  
  # 7. 绘图 (Visualizations)
  if(nrow(stats_table %>% dplyr::filter(Category != "Not_Significant")) > 0) {
    p_bar <- ggplot(stats_table %>% dplyr::filter(Category != "Not_Significant"), 
                    aes(x = Category, y = Count, fill = Category)) +
      geom_bar(stat = "identity", width = 0.6) +
      geom_text(aes(label = Count), vjust = -0.5) +
      scale_fill_manual(values = c("Resistant_Marker" = "#d73027", "Sensitive_Marker" = "#4575b4")) +
      (if (exists("prism_theme")) prism_theme else theme_classic()) +
      labs(title = paste0(prefix, ": Number of DE Genes (", filter_metric, " < ", padj_cutoff, ")"), 
           x = "", y = "Count") +
      theme(legend.position = "none")
    
    ggsave(file.path(output_dir, paste0(prefix, "_DE_Counts_Barplot.png")), p_bar, width = 5, height = 5)
  }
  
  return(list(
    clean_df = clean_res,
    stats = stats_table,
    top_resistant = top_resistant,
    top_sensitive = top_sensitive
  ))
}


plot_volcano <- function(
    res_df, output_dir, prefix = "Project",
    log2fc_thresh = 1, p_thresh = 0.05, top_n = 20,
    p_col = c("pvalue","padj"),
    label_sig_only = TRUE,        # label only significant genes (pass p & |LFC|)
    label_split    = TRUE,        # split labels to Up/Down halves
    gene_cols = c("Symbol","SYMBOL","Gene","gene","Gene.symbol","X.ID"),
    point_size_range = c(0.5, 2.5),
    alpha_points     = 0.6,
    palette          = c(Up = "#d73027", Down = "#4575b4", None = "grey80"),
    xlim_range       = c(-6, 6),  # NULL = auto; numeric range fixes x-axis
    save_png         = TRUE,
    show_counts      = TRUE,
    counts_labels    = c(Up = "Up", Down = "Down"),
    counts_fmt       = "{label}: {n}",
    counts_text_size = 4,
    counts_bold      = TRUE,
    arrow_length_frac = 0.12,   # arrow length as a fraction of x-range
    arrow_xpad_frac   = 0.04,   # horizontal padding from plot borders
    counts_y_frac     = 0.96,   # y position as fraction of max(-log10 p)
    arrow_size        = 0.6,
    arrow_colors      = c(Up = "#d73027", Down = "#4575b4")
){
  p_col <- match.arg(p_col)
  
  # --- gene symbol column detection ---
  sym_col <- gene_cols[gene_cols %in% names(res_df)]
  sym_col <- if (length(sym_col)) sym_col[1] else "Symbol"
  if (!sym_col %in% names(res_df)) res_df[[sym_col]] <- NA_character_
  
  # --- ensure adj.P exists when chosen ---
  if (p_col == "padj" && !"padj" %in% names(res_df)) {
    if (!"pvalue" %in% names(res_df)) stop("p_col='padj' selected but no 'padj' or 'pvalue' found.")
    res_df$padj <- p.adjust(res_df$pvalue, method = "BH")
  }
  
  # --- prepare data & flags ---
  df <- res_df %>%
    mutate(
      p_used   = .data[[p_col]],
      neglog10 = -log10(pmax(p_used, .Machine$double.xmin)),
      significant = case_when(
        log2FoldChange >  log2fc_thresh & p_used < p_thresh ~ "Up",
        log2FoldChange < -log2fc_thresh & p_used < p_thresh ~ "Down",
        TRUE ~ "None"
      ),
      significant = factor(significant, levels = c("Up","Down","None"))
    )
  
  up_count   <- sum(df$significant == "Up",   na.rm = TRUE)
  down_count <- sum(df$significant == "Down", na.rm = TRUE)
  
  # --- select labels ---
  lab_pool <- if (label_sig_only) df %>% filter(significant != "None") else df
  if (label_split) {
    n_each <- ceiling(top_n/2)
    lab_df <- bind_rows(
      lab_pool %>% filter(significant == "Up")   %>% arrange(p_used, desc(abs(log2FoldChange))) %>% slice_head(n = n_each),
      lab_pool %>% filter(significant == "Down") %>% arrange(p_used, desc(abs(log2FoldChange))) %>% slice_head(n = n_each)
    )
  } else {
    lab_df <- lab_pool %>% arrange(p_used, desc(abs(log2FoldChange))) %>% slice_head(n = top_n)
  }
  lab_df <- distinct(lab_df, .data[[sym_col]], .keep_all = TRUE)
  
  theme_use <- if (exists("prism_theme", inherits = TRUE)) get("prism_theme") else theme_bw()
  
  # --- base plot ---
  p <- ggplot(df, aes(log2FoldChange, neglog10)) +
    geom_point(aes(size = neglog10, color = significant), alpha = alpha_points) +
    geom_point(data = lab_df, aes(size = neglog10), shape = 21, color = "black", fill = "#ff7f00") +
    {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        ggrepel::geom_text_repel(
          data = lab_df, aes(label = .data[[sym_col]]),
          box.padding = 0.5, max.overlaps = Inf
        )
      } else {
        geom_text(data = lab_df, aes(label = .data[[sym_col]]), vjust = -0.4)
      }
    } +
    scale_color_manual(values = palette) +
    scale_size(range = point_size_range) +
    geom_vline(xintercept = c(-log2fc_thresh, log2fc_thresh), lty = 4) +
    geom_hline(yintercept = -log10(p_thresh), lty = 4) +
    {
      if (is.null(xlim_range)) NULL else xlim(xlim_range)
    } +
    labs(
      title = paste0(prefix, " - Volcano Plot (", if (p_col=="padj") "adj.P" else "P", ")"),
      x = "Log2 Fold Change", y = "-Log10 P-value"
    ) +
    theme_use
  
  # --- counts + arrows on top left/right ---
  if (show_counts) {
    # x-range used for positioning (respect fixed xlim if provided)
    xr <- if (!is.null(xlim_range)) {
      range(xlim_range)
    } else {
      range(df$log2FoldChange[is.finite(df$log2FoldChange)], na.rm = TRUE)
    }
    dx <- diff(xr)
    
    # y position near top
    y_top <- counts_y_frac * max(df$neglog10, na.rm = TRUE)
    
    # left (Down) arrow: point to the left
    xL0 <- xr[1] + dx * arrow_xpad_frac
    xL1 <- xL0 - dx * arrow_length_frac
    # right (Up) arrow: point to the right
    xR0 <- xr[2] - dx * arrow_xpad_frac
    xR1 <- xR0 + dx * arrow_length_frac
    
    # labels
    lab_down <- sub("\\{label\\}", counts_labels["Down"], sub("\\{n\\}", down_count, counts_fmt))
    lab_up   <- sub("\\{label\\}", counts_labels["Up"],   sub("\\{n\\}", up_count,   counts_fmt))
    
    p <- p +
      # left arrow + text
      annotate("segment", x = xL0, xend = xL1, y = y_top, yend = y_top,
               colour = arrow_colors["Down"], size = arrow_size,
               arrow = arrow(type = "closed", length = grid::unit(3, "mm"))) +
      annotate("text", x = xL0 + dx*0.01, y = y_top,
               label = lab_down, hjust = 0, vjust = 0.5,
               size = counts_text_size, fontface = if (counts_bold) "bold" else "plain",
               colour = arrow_colors["Down"]) +
      # right arrow + text
      annotate("segment", x = xR0, xend = xR1, y = y_top, yend = y_top,
               colour = arrow_colors["Up"], size = arrow_size,
               arrow = arrow(type = "closed", length = grid::unit(3, "mm"))) +
      annotate("text", x = xR0 - dx*0.01, y = y_top,
               label = lab_up, hjust = 1, vjust = 0.5,
               size = counts_text_size, fontface = if (counts_bold) "bold" else "plain",
               colour = arrow_colors["Up"])
  }
  
  # --- export ---
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, paste0(prefix, "_Volcano.pdf")), p, width = 7, height = 6, dpi = 300)
  if (save_png) ggsave(file.path(output_dir, paste0(prefix, "_Volcano.png")), p, width = 7, height = 6, dpi = 300)
  return(p)
}




run_go_enrich_pretty <- function(
    gene_symbols,
    ont = c("MF","BP","CC"),   # 可传向量或 "ALL"
    top_n = 10,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    universe_symbols = NULL,
    outdir = getwd(),
    out_prefix = NULL          # 作为基底前缀；会自动在后面加 _MF/_BP/_CC
){
  onts <- ont
  if (length(onts) == 1 && toupper(onts) %in% c("ALL","GOALL")) onts <- c("MF","BP","CC")
  valid_onts <- c("MF","BP","CC")
  onts <- toupper(onts)
  onts <- intersect(onts, valid_onts)
  if (length(onts) == 0) stop("ont 需为 'MF','BP','CC' 或 'ALL'。")
  
  # ---- 内部：SYMBOL→ENTREZ（去重，大小写容错）----
  map_sym_to_entrez <- function(symbols){
    symbols <- unique(toupper(symbols))
    if(length(symbols) == 0) return(character(0))
    tb <- suppressMessages(clusterProfiler::bitr(symbols,
                                                 fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
    tb <- tb[!duplicated(tb$SYMBOL), ]
    unname(tb$ENTREZID)
  }
  
  # 输入/背景 映射一次（复用）
  entrez_in <- map_sym_to_entrez(gene_symbols)
  if(length(entrez_in) == 0) stop("没有成功映射到 ENTREZID 的输入基因。")
  
  entrez_universe <- NULL
  if(!is.null(universe_symbols)){
    entrez_universe <- map_sym_to_entrez(universe_symbols)
    if(length(entrez_universe) == 0) entrez_universe <- NULL
  }
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  results <- list()
  combined_sheets <- list()
  
  for (ont_i in onts){
    # ----- ORA -----
    ego <- enrichGO(
      gene          = entrez_in,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = ont_i,
      universe      = entrez_universe,
      pvalueCutoff  = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      qvalueCutoff  = qvalueCutoff,
      minGSSize     = minGSSize,
      maxGSSize     = maxGSSize,
      readable      = TRUE
    )
    df <- as.data.frame(ego)
    if(nrow(df) == 0){
      message("【", ont_i, "】没有富集到条目。")
      next
    }
    
    term_size <- as.numeric(sub("/.*","", df$BgRatio))
    df$RichFactor <- df$Count / term_size
    df$pval_log  <- -log10(df$pvalue)
    
    df_top <- df |>
      dplyr::arrange(p.adjust, pvalue) |>
      dplyr::slice_head(n = min(top_n, nrow(df))) |>
      dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))
    
    p1 <- ggplot(df_top, aes(x = -log10(pvalue), y = Description)) +
      geom_bar(stat = "identity", fill = "#7f4ea850", width = 0.8) +
      geom_text(aes(x = 0.2, label = Description), hjust = 0, vjust = 1, size = 4.5) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        plot.margin  = margin(t = 10, r = 10, b = 10, l = 10),
        axis.line    = element_line(linewidth = 0.75, color = "black"),
        axis.ticks.x = element_line(linewidth = 0.75),
        axis.ticks.y = element_blank(),
        axis.text    = element_text(size = 12, color = "black"),
        axis.title   = element_text(size = 12, color = "black"),
        plot.title   = element_text(size = 12, color = "black"),
        panel.grid   = element_blank()
      ) +
      labs(x = "-log10(pvalue)", y = "", title = paste0("GO ", ont_i, " (Top ", nrow(df_top), ")")) +
      scale_x_continuous(expand = c(0, 0)) +
      coord_cartesian(clip = 'off')
    
    # ----- 图 2：双轴图 -----
    scale_factor <- max(df_top$pval_log, na.rm = TRUE) / max(df_top$RichFactor, na.rm = TRUE)
    if(!is.finite(scale_factor) || scale_factor == 0) scale_factor <- 1
    p2 <- ggplot(df_top, aes(y = reorder(Description, pval_log))) +
      geom_bar(aes(x = pval_log), stat = "identity", fill = "#7f4ea860") +
      geom_path(aes(x = RichFactor * scale_factor, y = reorder(Description, pval_log), group = 1),
                color = "black", linewidth = 0.75) +
      geom_point(aes(x = RichFactor * scale_factor),
                 shape = 21, fill = "#7f4ea8", color = "black", size = 3) +
      scale_x_continuous(
        name = "-log10(pvalue)",
        sec.axis = sec_axis(~ . / scale_factor, name = "Rich Factor")
      ) +
      labs(y = "", title = paste0("GO ", ont_i, " (Top ", nrow(df_top), ")")) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = NA, color = NA),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text  = element_text(color = "black", size = 12),
        axis.title = element_text(size = 12),
        axis.ticks = element_line(colour = "black", linewidth = 1)
      )
    
    # ----- 输出文件名：基底前缀 + _MF/_BP/_CC -----
    base_prefix <- if (is.null(out_prefix)) "GO" else out_prefix
    prefix_i <- paste0(base_prefix, "_", ont_i)
    
    writexl::write_xlsx(df,     file.path(outdir, paste0(prefix_i, "_enrich_all.xlsx")))
    writexl::write_xlsx(df_top, file.path(outdir, paste0(prefix_i, "_enrich_top", nrow(df_top), ".xlsx")))
    ggsave(file.path(outdir, paste0(prefix_i, "_bar.pdf")),      p1, width = 4,   height = 4.5)
    ggsave(file.path(outdir, paste0(prefix_i, "_dualaxis.pdf")), p2, width = 5.2, height = 4.5)
    
    # 合并 Excel 的 sheet
    combined_sheets[[paste0(ont_i, "_all")]] <- df
    combined_sheets[[paste0(ont_i, "_top", nrow(df_top))]] <- df_top
    
    results[[ont_i]] <- list(result = df, top = df_top, plots = list(bar = p1, dual = p2))
  }
  
  # 额外写一个合并 Excel（多 sheet）
  if (length(combined_sheets) > 0) {
    comb_file <- file.path(outdir, paste0(ifelse(is.null(out_prefix), "GO", out_prefix), "_enrich_combined.xlsx"))
    writexl::write_xlsx(combined_sheets, comb_file)
  }
  
  invisible(results)
}

#rounded_enrich_plot
rounded_enrich_plot <- function(
    gene_symbols,                     # 向量：SYMBOL
    signature_name   = "Signature",   # 标题
    db               = c("GO:BP","GO:MF","GO:CC","KEGG","Reactome",
                         "MSigDB:H","MSigDB:C2","MSigDB:C5:BP","MSigDB:C5:MF","MSigDB:C5:CC","MSigDB:C6","MSigDB:C7"),
    top_per_cat      = 5,
    pAdjustMethod    = "BH",
    pvalueCutoff     = 1,
    qvalueCutoff     = 1,
    minGSSize        = 10,
    maxGSSize        = 2000,
    universe_symbols = NULL,          # 背景（SYMBOL
    msigdb_release   = "2024.1.Hs",
    msigdb_local_dir = "E:/modifiedcodeR",
    reactome_fallback = "E:/modifiedcodeR/c2.cp.reactome.v2024.1.Hs.symbols.gmt",
    prefer_online    = TRUE,
    kegg_offline_gmt = NULL,          # 本地 KEGG GMT（无网时）
    # ------- 外观与输出 -------
    pal = c(  # 清新 6 色（命名，严格映射）
      "BP"        = "#c3e1e6",
      "MF"        = "#f3dfb7",
      "CC"        = "#dcc6dc",
      "KEGG"      = "#96c38e",
      "Reactome"  = "#f6c1c7",
      "MSigDB:H"  = "#c1e7e3"
    ),
    show_geneIDs = TRUE,              # 图上显示 geneID（斜体）
    geneID_size = 3.5,
    geneID_vjust = 2.6,
    geneID_x = 0.10,                  # geneID 文本的 x 位置（与 -log10 轴同单位）
    outdir  = getwd(),
    prefix  = "rounded_enrich",
    width   = 7, height = 6, dpi = 300
){
  # 依赖
  suppressPackageStartupMessages({
    library(gground); library(ggprism); library(tidyverse)
    library(clusterProfiler); library(org.Hs.eg.db)
  })
  .sanitize <- function(x) gsub("[^[:alnum:]_.-]+","_", x)
  .ensure_writable_dir <- function(path){
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    testfile <- file.path(path, paste0(".rwtest_", as.integer(runif(1,1e6,9e6)), ".tmp"))
    ok <- tryCatch({ writeBin(as.raw(0), testfile); file.remove(testfile); TRUE }, error=function(e) FALSE)
    if (!ok) {
      alt <- file.path(tempdir(), "rounded_enrich_outputs")
      dir.create(alt, recursive = TRUE, showWarnings = FALSE)
      message("⚠️ 输出目录不可写：", path, "  → 改存到临时目录：", alt)
      return(alt)
    }
    path
  }
  
  outdir_ok <- .ensure_writable_dir(outdir)
  s_prefix <- .sanitize(prefix)
  s_sig    <- .sanitize(signature_name)
  file_base <- file.path(outdir_ok, paste0(s_prefix, "_", s_sig))
  .sym2entrez_df <- function(symbols){
    symbols <- unique(na.omit(toupper(symbols)))
    if(length(symbols) == 0) return(tibble(SYMBOL=character(), ENTREZID=character()))
    tb <- suppressMessages(clusterProfiler::bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db))
    tb <- tb[!duplicated(tb$SYMBOL), , drop=FALSE]
    tibble::as_tibble(tb)
  }
  
  # —— 工具：下载文件（静默）
  .dl_quiet <- function(url, dest){
    tryCatch({
      suppressWarnings(utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE))
      file.exists(dest) && file.info(dest)$size > 0
    }, error = function(e) FALSE)
  }
  
  # —— 工具：MSigDB 文件名
  .msigdb_filename <- function(code, release){
    switch(code,
           "H"           = sprintf("h.all.v%s.symbols.gmt", release),
           "C2"          = sprintf("c2.cp.v%s.symbols.gmt", release),
           "C2:REACTOME" = sprintf("c2.cp.reactome.v%s.symbols.gmt", release),
           "C5:BP"       = sprintf("c5.go.bp.v%s.symbols.gmt", release),
           "C5:MF"       = sprintf("c5.go.mf.v%s.symbols.gmt", release),
           "C5:CC"       = sprintf("c5.go.cc.v%s.symbols.gmt", release),
           "C6"          = sprintf("c6.all.v%s.symbols.gmt", release),
           "C7"          = sprintf("c7.all.v%s.symbols.gmt", release),
           stop("未知 MSigDB 代码：", code)
    )
  }
  
  # —— 工具：获取 GMT（优先在线 → 本地；Reactome 专用回退）
  .get_gmt <- function(collection_code){
    fname <- .msigdb_filename(collection_code, msigdb_release)
    target <- file.path(outdir_ok, fname)  # 下载到可写目录
    if (file.exists(target)) return(target)
    if (isTRUE(prefer_online)) {
      url <- sprintf("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/%s/%s", msigdb_release, fname)
      if (.dl_quiet(url, target)) return(target)
    }
    local_path <- file.path(msigdb_local_dir, fname)
    if (file.exists(local_path)) return(local_path)
    if (collection_code == "C2:REACTOME" && file.exists(reactome_fallback)) return(reactome_fallback)
    stop("无法获取 GMT：", collection_code,
         "\n  尝试在线：", fname, "；本地：", local_path,
         if(collection_code=="C2:REACTOME") paste0("；Reactome回退：", reactome_fallback) else "")
  }
  
  .make_pretty_terms <- function(x){
    if (length(x) == 0) return(x)
    x <- gsub("_", " ", x, fixed = TRUE)
    x <- gsub("-", " ", x, fixed = TRUE)
    x <- trimws(x)
    x <- tools::toTitleCase(tolower(x))
    keep_low <- c("and","or","of","in","to","via","the","for","by","on","with","from","into","a","an")
    prett <- vapply(strsplit(x, "\\s+"), function(w){
      if (length(w) == 0) return("")
      w[-1] <- ifelse(tolower(w[-1]) %in% keep_low, tolower(w[-1]), w[-1])
      paste(w, collapse = " ")
    }, character(1))
    trimws(prett)
  }
  genes_syms <- unique(na.omit(toupper(gene_symbols)))
  map_df <- .sym2entrez_df(genes_syms)
  genes_entrez <- unique(map_df$ENTREZID)
  unmapped_syms <- setdiff(genes_syms, unique(map_df$SYMBOL))
  
  message(sprintf("SYMBOL→ENTREZ 映射：%d / %d (%.2f%%) 成功；未映射 %d",
                  length(genes_entrez), length(genes_syms),
                  if (length(genes_syms)==0) 0 else 100*length(genes_entrez)/length(genes_syms),
                  length(unmapped_syms)))
  readr::write_csv(map_df, paste0(file_base, "_mapping_SYMBOL2ENTREZ.csv"))
  if (length(unmapped_syms) > 0) {
    writeLines(unmapped_syms, paste0(file_base, "_unmapped_SYMBOL.txt"))
  }
  universe_symbols <- if(is.null(universe_symbols)) NULL else unique(toupper(universe_symbols))
  entrez_universe <- if(!is.null(universe_symbols)) {
    uni_map <- .sym2entrez_df(universe_symbols)
    unique(uni_map$ENTREZID)
  } else NULL
  all_tables <- list()
  
  for (d in unique(db)) {
    df <- NULL
    
    if (startsWith(d, "GO:")) {
      ont <- sub("^GO:", "", d)  # BP/MF/CC
      if (length(genes_entrez) == 0) next
      eg <- tryCatch({
        enrichGO(gene = genes_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                 ont = ont, universe = entrez_universe,
                 pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
                 minGSSize = minGSSize, maxGSSize = maxGSSize)
      }, error = function(e) NULL)
      if (!is.null(eg)) eg <- setReadable(eg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      df <- if(!is.null(eg)) as.data.frame(eg) else data.frame()
      if (nrow(df) > 0) df$ONTOLOGY <- ont
      
    } else if (d == "KEGG") {
      if (!isTRUE(prefer_online) && !is.null(kegg_offline_gmt) && file.exists(kegg_offline_gmt)) {
        # 离线 KEGG：GMT + enricher（基于 SYMBOL）
        TERM2GENE <- clusterProfiler::read.gmt(kegg_offline_gmt)
        colnames(TERM2GENE) <- c("term","gene")
        eg <- tryCatch({
          enricher(gene = genes_syms, TERM2GENE = TERM2GENE, universe = universe_symbols,
                   pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, pAdjustMethod = pAdjustMethod,
                   minGSSize = minGSSize, maxGSSize = maxGSSize)
        }, error = function(e) NULL)
        df <- if(!is.null(eg)) as.data.frame(eg) else data.frame()
      } else {
        # 在线 KEGG（ENTREZ）+ setReadable → SYMBOL（静默请求提示）
        if (length(genes_entrez) == 0) next
        eg <- tryCatch({
          suppressMessages(
            enrichKEGG(gene = genes_entrez, organism = "hsa", keyType = "ncbi-geneid",
                       universe = entrez_universe, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod,
                       qvalueCutoff = qvalueCutoff, minGSSize = minGSSize, maxGSSize = maxGSSize)
          )
        }, error = function(e) NULL)
        if (!is.null(eg)) eg <- setReadable(eg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        df <- if(!is.null(eg)) as.data.frame(eg) else data.frame()
      }
      if (nrow(df) > 0) df$ONTOLOGY <- "KEGG"
      
    } else {
      # Reactome / MSigDB：GMT + enricher（基于 SYMBOL）
      coll_code <- switch(d,
                          "Reactome"     = "C2:REACTOME",
                          "MSigDB:H"     = "H",
                          "MSigDB:C2"    = "C2",
                          "MSigDB:C5:BP" = "C5:BP",
                          "MSigDB:C5:MF" = "C5:MF",
                          "MSigDB:C5:CC" = "C5:CC",
                          "MSigDB:C6"    = "C6",
                          "MSigDB:C7"    = "C7",
                          stop("不支持的 db：", d)
      )
      gmt_path <- .get_gmt(coll_code)
      TERM2GENE <- clusterProfiler::read.gmt(gmt_path)
      colnames(TERM2GENE) <- c("term","gene")
      # 去掉前缀（Reactome/Hallmark），其余保持原有集合名
      if (coll_code == "C2:REACTOME") TERM2GENE$term <- sub("^REACTOME_", "", TERM2GENE$term)
      if (coll_code == "H")           TERM2GENE$term <- sub("^HALLMARK_", "", TERM2GENE$term)
      
      eg <- tryCatch({
        enricher(gene = genes_syms, TERM2GENE = TERM2GENE, universe = universe_symbols,
                 pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, pAdjustMethod = pAdjustMethod,
                 minGSSize = minGSSize, maxGSSize = maxGSSize)
      }, error = function(e) NULL)
      df <- if(!is.null(eg)) as.data.frame(eg) else data.frame()
      if (nrow(df) > 0) {
        df$ONTOLOGY <- d
        # 术语美化：对 Reactome / MSigDB（尤其 H）转为可发表标题
        df$Description <- .make_pretty_terms(df$Description)
      }
    }
    
    if (!is.null(df) && nrow(df) > 0) {
      if (!"p.adjust" %in% colnames(df) && "qvalue" %in% colnames(df)) df$p.adjust <- df$qvalue
      all_tables[[d]] <- df
    }
  }
  
  if (length(all_tables) == 0) {
    warning("没有任何数据库富集到条目。")
    return(invisible(NULL))
  }
  
  # —— 汇总与 TopN（修复 slice_head 常数 & Count 兜底）
  GO <- dplyr::bind_rows(all_tables[names(all_tables) %in% c("GO:BP","GO:MF","GO:CC")])
  KEGG <- if("KEGG" %in% names(all_tables)) all_tables[["KEGG"]] else NULL
  OTHERS <- dplyr::bind_rows(all_tables[!(names(all_tables) %in% c("GO:BP","GO:MF","GO:CC","KEGG"))])
  
  bind_all <- dplyr::bind_rows(GO, KEGG, OTHERS) %>%
    dplyr::mutate(
      Description = as.character(Description),
      p.adjust = suppressWarnings(as.numeric(p.adjust)),
      Count = dplyr::coalesce(as.numeric(Count), lengths(strsplit(geneID, "/")))
    ) %>%
    dplyr::filter(is.finite(p.adjust))
  
  use_pathway <- bind_all %>%
    dplyr::group_by(ONTOLOGY) %>%
    dplyr::arrange(p.adjust, dplyr::desc(Count), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_per_cat) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ONTOLOGY = factor(
      ONTOLOGY,
      levels = c(intersect(c("BP","CC","MF"), unique(ONTOLOGY)),
                 intersect("KEGG", unique(ONTOLOGY)),
                 setdiff(unique(ONTOLOGY), c("BP","CC","MF","KEGG")))
    )) %>%
    dplyr::arrange(ONTOLOGY, p.adjust, dplyr::desc(Count)) %>%
    dplyr::mutate(Description = factor(Description, levels = Description)) %>%
    tibble::rowid_to_column('index')
  
  # —— 左侧分类条数据
  width_left <- 0.5
  xaxis_max <- max(-log10(use_pathway$p.adjust), na.rm = TRUE) + 1
  rect.data <- use_pathway %>%
    dplyr::group_by(ONTOLOGY) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      xmin = -3 * width_left, xmax = -2 * width_left,
      ymax = cumsum(n),
      ymin = dplyr::lag(ymax, default = 0) + 0.6,
      ymax = ymax + 0.4
    )
  
  # —— 颜色：按名字精确映射；缺失类别自动补足
  levs <- levels(use_pathway$ONTOLOGY)
  if (!is.null(names(pal))) {
    pal_named <- pal
    missing <- setdiff(levs, names(pal_named))
    if (length(missing) > 0) {
      base_cols <- unname(pal)
      if (length(base_cols) == 0) base_cols <- RColorBrewer::brewer.pal(max(3, length(missing)), "Set2")
      pal_named <- c(pal_named, stats::setNames(rep_len(base_cols, length(missing)), missing))
    }
    pal_use <- pal_named[levs]
  } else {
    pal_use <- rep_len(pal, length(levs))
  }
  
  # —— 画图
  p <- use_pathway %>%
    ggplot2::ggplot(ggplot2::aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    gground::geom_round_col(ggplot2::aes(y = Description), width = 0.6, alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(x = 0.05, label = Description), hjust = 0, size = 5) +
    { if (isTRUE(show_geneIDs) && "geneID" %in% colnames(use_pathway))
      ggplot2::geom_text(ggplot2::aes(x = geneID_x, label = geneID, colour = ONTOLOGY),
                         hjust = 0, vjust = geneID_vjust, size = geneID_size,
                         fontface = 'italic', show.legend = FALSE) } +
    ggplot2::geom_point(ggplot2::aes(x = -width_left, size = Count), shape = 21) +
    ggplot2::geom_text(ggplot2::aes(x = -width_left, label = Count), size = 3.5) +
    ggplot2::scale_size_continuous(name = 'Count', range = c(5, 12)) +
    gground::geom_round_rect(
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
      data = rect.data, radius = grid::unit(2,'mm'), inherit.aes = FALSE, alpha = 0.95
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
      data = rect.data, inherit.aes = FALSE
    ) +
    ggplot2::annotate("segment", x = 0, y = 0, xend = xaxis_max, yend = 0, linewidth = 1.5) +
    ggplot2::labs(
      x = expression(-log[10](adj.~pvalue)), y = NULL,
      title = paste0(signature_name, " | Enrichment (ORA)")
    ) +
    ggplot2::scale_fill_manual(name = 'Category', values = pal_use, drop = FALSE) +
    ggplot2::scale_colour_manual(values = pal_use, guide = "none") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(), expand = ggplot2::expansion(c(0,0))) +
    ggprism::theme_prism() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  # —— 输出（稳健保存：若失败，自动回落到临时目录）
  .safe_save <- function(plot_obj, base, width, height, dpi){
    ok_pdf <- tryCatch({ ggplot2::ggsave(paste0(base, ".pdf"), plot_obj, width = width, height = height); TRUE }, error=function(e) FALSE)
    ok_png <- tryCatch({ ggplot2::ggsave(paste0(base, ".png"), plot_obj, width = width, height = height, dpi = dpi); TRUE }, error=function(e) FALSE)
    if (!(ok_pdf && ok_png)) {
      alt_dir <- .ensure_writable_dir(tempdir())
      alt_base <- file.path(alt_dir, basename(base))
      ggplot2::ggsave(paste0(alt_base, ".pdf"), plot_obj, width = width, height = height)
      ggplot2::ggsave(paste0(alt_base, ".png"), plot_obj, width = width, height = height, dpi = dpi)
      message("⚠️ 已改存到临时目录：", alt_base, ".pdf / .png")
      return(invisible(c(paste0(alt_base, ".pdf"), paste0(alt_base, ".png"))))
    }
    invisible(c(paste0(base, ".pdf"), paste0(base, ".png")))
  }
  
  .safe_write_csv <- function(df, path){
    ok <- tryCatch({ readr::write_csv(df, path); TRUE }, error=function(e) FALSE)
    if (!ok) {
      alt <- file.path(.ensure_writable_dir(tempdir()), basename(path))
      readr::write_csv(df, alt)
      message("⚠️ 表格已改存到临时目录：", alt)
      return(invisible(alt))
    }
    invisible(path)
  }
  
  .safe_save(p, file_base, width, height, dpi)
  .safe_write_csv(use_pathway, paste0(file_base, "_top.csv"))
  
  invisible(list(plot = p, table = use_pathway,
                 mapping = map_df, unmapped = unmapped_syms,
                 out_base = file_base))
}
# 
# 
# pal6 <- c(
#   "BP"        = "#c3e1e6",  # 清新浅蓝
#   "MF"        = "#f3dfb7",  # 柔和米杏
#   "CC"        = "#dcc6dc",  # 淡薰衣草
#   "KEGG"      = "#96c38e",  # 鼠尾草绿
#   "Reactome"  = "#f6c1c7",  # 珊瑚粉
#   "MSigDB:H"  = "#c1e7e3"   # 薄荷青
# )
# 
# res <- rounded_enrich_plot(
#   gene_symbols     = sig_genes,            # ← 建议确认就是 geneVEnn
#   signature_name   = "geneVEnn",
#   db               = c("GO:BP","GO:MF","GO:CC","KEGG","Reactome","MSigDB:H"),
#   top_per_cat      = 5,
#   msigdb_release   = "2024.1.Hs",
#   msigdb_local_dir = "E:/modifiedcodeR",
#   reactome_fallback= "E:/modifiedcodeR/c2.cp.reactome.v2024.1.Hs.symbols.gmt",
#   prefer_online    = TRUE,
#   pal              = pal6,                 # ← 用 pal=...
#   show_geneIDs     = TRUE,
#   outdir           = outdir,
#   prefix           = "rounded_enrich_geneVEnn",
#   width = 12, height = 14
# )

prettify_hallmark_terms <- function(x) {
  # 1. 基础清理
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x, fixed = TRUE)
  
  # 2. Title Case 转换 (介词小写)
  to_title_smart <- function(s) {
    s <- tools::toTitleCase(tolower(s))
    # 定义不需要大写的介词/冠词
    keep_low <- c("and", "or", "of", "in", "to", "via", "the", "for", "by", "on", "with", "from", "into", "a", "an", "at", "as")
    
    # 分割并处理
    w <- strsplit(s, "\\s+")[[1]]
    if (length(w) > 1) {
      # 跳过第一个词，检查其余词
      w[-1] <- ifelse(tolower(w[-1]) %in% keep_low, tolower(w[-1]), w[-1])
    }
    paste(w, collapse = " ")
  }
  
  # 对向量应用 Title Case
  x <- vapply(x, to_title_smart, character(1))
  
  # 3. 常见生物学术语/符号修正
  x <- gsub("\\bG2M\\b", "G2/M", x, ignore.case = TRUE)
  x <- gsub("\\bNFKB\\b", "NF-κB", x, ignore.case = TRUE)
  x <- gsub("\\bTNFA\\b", "TNF-α", x, ignore.case = TRUE)
  x <- gsub("\\bTGF BETA\\b", "TGF-β", x, ignore.case = TRUE)
  x <- gsub("\\bINTERFERON ALPHA\\b", "Interferon-α", x, ignore.case = TRUE)
  x <- gsub("\\bINTERFERON GAMMA\\b", "Interferon-γ", x, ignore.case = TRUE)
  x <- gsub("(?i)\\bIL[- ]?(\\d+)\\b", "IL-\\1", x, perl = TRUE) # IL6 -> IL-6
  x <- gsub("(?i)\\bMTORC1\\b", "mTORC1", x, perl = TRUE)
  x <- gsub("(?i)\\bMTOR\\b",   "mTOR",   x, perl = TRUE)
  x <- gsub("(?i)\\bJAK\\b",   "JAK",   x, perl = TRUE)
  x <- gsub("(?i)\\bSTAT(\\d)\\b", "STAT\\1", x, perl = TRUE)
  x <- gsub("(?i)\\bNOTCH\\b",   "Notch",   x, perl = TRUE)
  
  # 4. 基因名大写修正 (KRAS, MYC, etc.)
  x <- gsub("\\bMyc\\b",  "MYC",  x)
  x <- gsub("\\bKras\\b", "KRAS", x)
  x <- gsub("\\bE2f\\b",  "E2F",  x)
  x <- gsub("\\bWnt\\b",  "WNT",  x)
  x <- gsub("\\bPi3k\\b", "PI3K", x)
  x <- gsub("\\bAkt\\b",  "AKT",  x)
  x <- gsub("\\bMapk\\b", "MAPK", x)
  x <- gsub("\\bRos\\b",  "ROS",  x)
  x <- gsub("\\bDna\\b",  "DNA",  x)
  x <- gsub("\\bUv\\b",   "UV",   x)
  
  # 5. 上下调标记 (可选)
  x <- gsub("(?i)\\bUP\\b$", "(Up)", x, perl = TRUE)
  x <- gsub("(?i)\\bDN\\b$", "(Down)", x, perl = TRUE)
  
  return(x)
}
# --- 把数据库名字印得更友好（用于标题）---
pretty_db_label <- function(db){
  d <- toupper(db)
  if (grepl("^GO_", d)) return(paste("GO", sub("^GO_", "", db)))
  if (d == "KEGG")      return("KEGG")
  if (d == "REACTOME")  return("Reactome")
  if (d == "HALLMARK")  return("MSigDB Hallmark")
  db
}
plot_gsea_bubble <- function(
    data,
    x_var = "Group_Label",   
    y_var = "ID",            
    size_var = "NES",        
    color_var = "NES",       
    
    # 核心排序功能
    order_method = c("mean_asc", "mean_desc", "abs_mean", "cluster"), 
    
    # --- 数据清洗 (接入了你的高级美化函数) ---
    clean_names = TRUE,      
    
    # 视觉参数
    colors = c("#313695", "white", "#A50026"), 
    bubble_range = c(2, 8),  
    title = "Hallmark Pathways Enrichment",
    subtitle = "Comparing Treatment vs Control",
    x_lab = "",
    y_lab = ""
) {
  
  # 0. 检查参数
  order_method <- match.arg(order_method)
  plot_data <- data
  
  # 1. 高级名称美化 (Integration Point)
  if(clean_names) {
    # 这里直接调用我们在外部定义的 prettify_hallmark_terms
    if(exists("prettify_hallmark_terms")) {
      plot_data[[y_var]] <- prettify_hallmark_terms(plot_data[[y_var]])
    } else {
      # 如果函数没定义，回退到简单清洗并警告
      warning("prettify_hallmark_terms function not found. Using basic cleaning.")
      plot_data[[y_var]] <- gsub("HALLMARK_", "", plot_data[[y_var]])
      plot_data[[y_var]] <- gsub("_", " ", plot_data[[y_var]])
    }
  }
  
  # 2. 核心：计算排序逻辑
  # 2.1 准备聚合数据
  agg_data <- plot_data %>%
    group_by(.data[[y_var]]) %>%
    summarise(
      Mean_Val = mean(.data[[color_var]], na.rm = TRUE),
      Abs_Mean_Val = mean(abs(.data[[color_var]]), na.rm = TRUE),
      .groups = "drop"
    )
  
  final_levels <- NULL
  
  if (order_method == "cluster") {
    # --- 聚类模式 ---
    mat_data <- plot_data %>%
      dplyr::select(all_of(c(y_var, x_var, color_var))) %>%
      pivot_wider(names_from = all_of(x_var), values_from = all_of(color_var), values_fill = 0) %>%
      column_to_rownames(var = y_var)
    
    dist_mat <- dist(mat_data)
    hc <- hclust(dist_mat)
    final_levels <- hc$labels[hc$order] 
    
  } else if (order_method == "mean_asc") {
    # --- 平均值从小到大 (从上往下) ---
    final_levels <- agg_data %>% arrange(desc(Mean_Val)) %>% pull(.data[[y_var]])
    
  } else if (order_method == "mean_desc") {
    # --- 平均值从大到小 (从上往下) ---
    final_levels <- agg_data %>% arrange(Mean_Val) %>% pull(.data[[y_var]])
    
  } else if (order_method == "abs_mean") {
    # --- 绝对值从小到大 (从上往下) ---
    final_levels <- agg_data %>% arrange(desc(Abs_Mean_Val)) %>% pull(.data[[y_var]])
  }
  
  # 应用排序
  plot_data[[y_var]] <- factor(plot_data[[y_var]], levels = final_levels)
  
  # 3. 绘图
  p <- ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(size = abs(.data[[size_var]]), color = .data[[color_var]])) +
    scale_color_gradient2(low = colors[1], mid = colors[2], high = colors[3], midpoint = 0) +
    scale_size(range = bubble_range) +
    theme_bw() +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_lab, 
      y = y_lab,
      color = "NES",
      size = "|NES|"
    ) +
    theme(
      axis.text.x = element_text(size = 11, color = "black", face = "bold"),
      axis.text.y = element_text(size = 10, color = "black"), # 字体稍大，因为名字变漂亮了
      panel.grid.major = element_line(linetype = "dashed", color = "grey90"),
      plot.title = element_text(face = "bold", size = 14)
    )
  
  return(p)
}




suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
  library(writexl)
})
RunEnrichmentSuite <- function(
    gene_symbols,
    db_options = c("GO_BP", "KEGG", "Reactome", "Hallmark"),
    top_n = 10,
    pvalueCutoff = 1.0,
    qvalueCutoff = 1.0,
    minGSSize = 10,
    maxGSSize = 500,
    pAdjustMethod = "BH",
    universe_symbols = NULL,
    outdir = getwd(),
    out_prefix = "Enrichment",
    # >>> 新增：标题标签 & 组色 <<<
    title_tag   = NULL,     # 例如 "YourLab · Basal-like"
    plot_color  = "#7f4ea8" # 该组的主题色（来自 custom_palette）
){
  message("--- (1/4) Converting gene SYMBOLs to ENTREZIDs ---")
  n_genes_before <- length(unique(toupper(gene_symbols)))
  entrez_in <- suppressMessages(clusterProfiler::bitr(
    unique(toupper(gene_symbols)), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
  ))
  if(nrow(entrez_in) == 0) { warning("没有成功映射到 ENTREZID 的输入基因。"); return(invisible(NULL)) }
  
  entrez_universe <- NULL
  if(!is.null(universe_symbols)){
    entrez_universe_df <- suppressMessages(clusterProfiler::bitr(
      unique(toupper(universe_symbols)), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db
    ))
    entrez_universe <- entrez_universe_df$ENTREZID
  }
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  final_results <- list(); all_sheets <- list()
  
  message("--- (2/4) Running enrichment analysis for selected databases ---")
  for (db in db_options) {
    message(sprintf("  - Processing database: %s", db))
    ego <- NULL
    if (grepl("^GO_", db, ignore.case = TRUE)) {
      ont_type <- sub("GO_", "", db); if (!ont_type %in% c("BP","MF","CC")) next
      ego <- enrichGO(entrez_in$ENTREZID, universe=entrez_universe, OrgDb=org.Hs.eg.db,
                      keyType="ENTREZID", ont=ont_type, pvalueCutoff=pvalueCutoff,
                      pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff,
                      minGSSize=minGSSize, maxGSSize=maxGSSize, readable=TRUE)
    } else if (toupper(db) == "KEGG") {
      ego <- enrichKEGG(entrez_in$ENTREZID, universe=entrez_universe, organism='hsa', keyType='kegg',
                        pvalueCutoff=pvalueCutoff, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff,
                        minGSSize=minGSSize, maxGSSize=maxGSSize)
      if (!is.null(ego)) ego <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    } else if (toupper(db) == "REACTOME") {
      ego <- ReactomePA::enrichPathway(entrez_in$ENTREZID, universe=entrez_universe, organism="human",
                                       pvalueCutoff=pvalueCutoff, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff,
                                       minGSSize=minGSSize, maxGSSize=maxGSSize, readable=TRUE)
    } else if (toupper(db) == "HALLMARK") {
      m_t2g <- msigdbr(species="Homo sapiens", category="H") %>% dplyr::select(gs_name, entrez_gene)
      ego <- enricher(entrez_in$ENTREZID, universe=entrez_universe, TERM2GENE=m_t2g,
                      pvalueCutoff=pvalueCutoff, pAdjustMethod=pAdjustMethod, qvalueCutoff=qvalueCutoff,
                      minGSSize=minGSSize, maxGSSize=maxGSSize)
      if (!is.null(ego)) ego <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    }
    
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) { message(sprintf("    【%s】无显著条目。", db)); next }
    df <- as.data.frame(ego)
    
    # >>> Hallmark 名称清洗（出版风） <<<
    if (toupper(db) == "HALLMARK" && "Description" %in% names(df)) {
      df$Description <- prettify_hallmark_terms(df$Description)
    }
    
    # 指标
    if ("BgRatio" %in% names(df)) {
      term_size <- as.numeric(sub("/.*", "", df$BgRatio))
      df$RichFactor <- df$Count / term_size
    } else if (all(c("GeneRatio", "BgRatio") %in% names(df))) {
      num_in <- as.numeric(sub("/.*", "", df$GeneRatio))
      term_size <- as.numeric(sub("/.*", "", df$BgRatio))
      df$RichFactor <- num_in / term_size
    } else df$RichFactor <- NA_real_
    df$pval_log <- -log10(df$pvalue)
    
    df_top <- df %>% arrange(p.adjust, pvalue) %>% slice_head(n = min(top_n, nrow(df))) %>%
      mutate(Description = factor(Description, levels = rev(unique(Description))))
    
    # 颜色（来自组色）
    bar_fill  <- grDevices::adjustcolor(plot_color, alpha.f = 0.55)
    bar_fill2 <- grDevices::adjustcolor(plot_color, alpha.f = 0.70)
    
    # 标题：组名标签 + 数据库名
    ttl <- paste0(if (!is.null(title_tag)) paste0(title_tag, " — ") else "",
                  pretty_db_label(db), " (Top ", nrow(df_top), ")")
    
    # 原风格 + 组色
    p1 <- ggplot(df_top, aes(x = pval_log, y = Description)) +
      geom_bar(stat = "identity", fill = bar_fill, width = 0.8) +
      geom_text(aes(x = 0.05 * max(pval_log, na.rm=TRUE), label = Description), hjust = 0, size = 4) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_blank(),
        plot.margin  = margin(t = 10, r = 10, b = 10, l = 10),
        axis.line    = element_line(linewidth = 0.75, color = "black"),
        axis.ticks.x = element_line(linewidth = 0.75),
        axis.ticks.y = element_blank(),
        axis.text    = element_text(size = 12, color = "black"),
        axis.title   = element_text(size = 12, color = "black"),
        panel.grid   = element_blank(),
        plot.title   = element_text(hjust = 0.5, face = "bold", colour = plot_color)
      ) +
      labs(x = "-log10(pvalue)", y = "", title = ttl) +
      scale_x_continuous(expand = c(0, 0.01))
    
    scale_factor <- max(df_top$pval_log, na.rm = TRUE) / max(df_top$RichFactor, na.rm = TRUE)
    if(!is.finite(scale_factor) || scale_factor == 0) scale_factor <- 1
    
    p2 <- ggplot(df_top, aes(y = Description)) +
      geom_bar(aes(x = pval_log), stat = "identity", fill = bar_fill2) +
      geom_path(aes(x = RichFactor * scale_factor, group = 1), color = "black", linewidth = 0.75) +
      geom_point(aes(x = RichFactor * scale_factor), shape = 21, fill = plot_color, color = "black", size = 3) +
      scale_x_continuous(name = "-log10(pvalue)", sec.axis = sec_axis(~ . / scale_factor, name = "Rich Factor")) +
      labs(y = "", title = ttl) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = NA, color = NA),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text  = element_text(color = "black", size = 12),
        axis.title = element_text(size = 12),
        axis.ticks = element_line(colour = "black", linewidth = 1),
        plot.title = element_text(hjust = 0.5, face = "bold", colour = plot_color)
      )
    
    prefix_i <- paste0(out_prefix, "_", db)
    writexl::write_xlsx(list(enrichment_results = df), file.path(outdir, paste0(prefix_i, "_all.xlsx")))
    ggsave(file.path(outdir, paste0(prefix_i, "_bar.pdf")),      p1, width = 4,   height = max(4, 0.25 * nrow(df_top) + 1))
    ggsave(file.path(outdir, paste0(prefix_i, "_dualaxis.pdf")), p2, width = 5.2, height = max(4, 0.25 * nrow(df_top) + 1))
    # --- NEW: save tabular results alongside plots (TSV + RDS) ---
    utils::write.table(
      df,
      file.path(outdir, paste0(prefix_i, "_all.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    utils::write.table(
      df_top,
      file.path(outdir, paste0(prefix_i, "_top", min(top_n, nrow(df_top)), ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    saveRDS(
      list(all = df, top = df_top),
      file.path(outdir, paste0(prefix_i, "_tables.rds"))
    )
    all_sheets[[db]] <- df
    final_results[[db]] <- list(result = df, plot_bar = p1, plot_dual = p2)
  }
  
  message("--- (3/4) Saving combined results to a single Excel file ---")
  if (length(all_sheets) > 0) {
    writexl::write_xlsx(all_sheets, file.path(outdir, paste0(out_prefix, "_ALL_DATABASES_combined.xlsx")))
  }
  message("--- (4/4) Enrichment suite analysis complete! ---")
  invisible(final_results)
}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#        multi_db_text_enrich_plot (v2 - 字体大小优化版)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
multi_db_text_enrich_plot <- function(
    sig_list,
    db = c("Reactome","GO:BP","GO:MF","GO:CC","KEGG",
           "MSigDB:H","MSigDB:C2","MSigDB:C5:BP","MSigDB:C5:MF","MSigDB:C5:CC","MSigDB:C6","MSigDB:C7"),
    top_n = 10,
    universe_symbols = NULL,
    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
    minGSSize = 10, maxGSSize = 2000,
    msigdb_release   = "2024.1.Hs",
    msigdb_local_dir = "E:/modifiedcodeR",
    reactome_fallback = "E:/modifiedcodeR/c2.cp.reactome.v2024.1.Hs.symbols.gmt",
    prefer_online    = TRUE,
    msigdb_manifest  = NULL,
    kegg_offline_gmt = NULL,
    title_map = NULL,
    color_map = NULL,
    default_low_high = c("#98bf92","#006a01"),
    
    # --- 字体大小控制参数 ---
    size_mult = 1.0,          # [修改] 全局基础字号倍数 (现在是乘数，不是加数)
    font_size_range = c(3, 8), # [新增] 将 -log10(padj) 缩放到这个字体大小范围内
    
    show_legends = TRUE,
    outdir = getwd(), prefix = "PathTextEnrich",
    save_single = TRUE, save_combined = TRUE,
    width_each = 5.0, height_each = 3.0, dpi = 300
){
  # ... (前面的依赖包加载和辅助函数保持不变) ...
  pkgs <- c("clusterProfiler","org.Hs.eg.db","ggplot2","dplyr","tibble","patchwork")
  for(p in pkgs){ if(!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE) }
  lapply(pkgs, require, character.only = TRUE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  .msigdb_filename <- function(code, release){
    switch(code,
           "H"           = sprintf("h.all.v%s.symbols.gmt", release),
           "C2"          = sprintf("c2.cp.v%s.symbols.gmt", release),
           "C2:REACTOME" = sprintf("c2.cp.reactome.v%s.symbols.gmt", release),
           "C5:BP"       = sprintf("c5.go.bp.v%s.symbols.gmt", release),
           "C5:MF"       = sprintf("c5.go.mf.v%s.symbols.gmt", release),
           "C5:CC"       = sprintf("c5.go.cc.v%s.symbols.gmt", release),
           "C6"          = sprintf("c6.all.v%s.symbols.gmt", release),
           "C7"          = sprintf("c7.all.v%s.symbols.gmt", release),
           stop("未知 MSigDB 集合代码：", code)
    )
  }
  .dl_quiet <- function(url, dest){
    tryCatch({
      suppressWarnings(utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE))
      file.exists(dest) && file.info(dest)$size > 0
    }, error = function(e) FALSE)
  }
  manifest <- NULL
  if (!is.null(msigdb_manifest) && file.exists(msigdb_manifest)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
    manifest <- tryCatch(jsonlite::read_json(msigdb_manifest, simplifyVector = TRUE), error = function(e) NULL)
  }
  .get_gmt <- function(collection_code){
    if (!is.null(manifest) && !is.null(manifest[[collection_code]]) && !is.na(manifest[[collection_code]])) {
      fp <- manifest[[collection_code]]
      if (file.exists(fp)) return(fp)
    }
    fname <- .msigdb_filename(collection_code, msigdb_release)
    target <- file.path(outdir, fname)
    if (file.exists(target)) return(target)
    if (prefer_online) {
      url <- sprintf("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/%s/%s", msigdb_release, fname)
      if (.dl_quiet(url, target)) return(target)
    }
    local_path <- file.path(msigdb_local_dir, fname)
    if (file.exists(local_path)) return(local_path)
    if (collection_code == "C2:REACTOME" && file.exists(reactome_fallback)) return(reactome_fallback)
    stop("无法获取 GMT：", collection_code,
         "\n  尝试：", sprintf("MSigDB/%s", fname), "（在线）；本地：", local_path,
         if(collection_code=="C2:REACTOME") paste0("；Reactome回退：", reactome_fallback) else "")
  }

  
  # ---------- (核心修改点) 统一绘图面板 ----------
  .panel_plot <- function(df, sig_name, low_high, mult, size_range){
    if (is.null(df) || nrow(df) == 0) return(NULL)
    
    df <- df |>
      dplyr::arrange(p.adjust, pvalue) |>
      dplyr::slice_head(n = min(top_n, nrow(df)))
    
    df$Description <- factor(df$Description, levels = df$Description)
    df$pval_log <- -log10(df$p.adjust)

    # --- 字体大小缩放逻辑 ---
    # 使用 scales::rescale 将 pval_log 线性映射到您定义的 font_size_range
    # `from` 参数设定了 pval_log 的原始范围，这里我们取当前 top_n 个通路中的范围
    min_logp <- min(df$pval_log, na.rm = TRUE)
    max_logp <- max(df$pval_log, na.rm = TRUE)
    
    # 防止只有一个通路时 min 和 max 相等导致除零错误
    if (min_logp == max_logp) {
      df$txt_size <- mean(size_range) * mult
    } else {
      df$txt_size <- scales::rescale(df$pval_log, to = size_range, from = c(min_logp, max_logp)) * mult
    }

    y_top <- nrow(df) + 0.6
    
    ggplot(df, aes(x = 1, y = rev(Description))) +
      # 使用新的 txt_size
      geom_text(aes(label = Description, colour = pval_log, size = txt_size), hjust = 0.5) +
      scale_color_gradient(low = low_high[1], high = low_high[2],
                           name = "Significance\n(-log10 adj. p-val.)") +
      scale_size_identity() + # size_identity() 意味着直接使用 txt_size 列的值作为字号
      scale_x_continuous(expand = c(0,0)) +
      labs(x = "", y = "", title = sig_name) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"), # 标题加粗
        legend.position = if (show_legends) "right" else "none"
      ) +
      annotate("segment", x = 0, xend = 2, y = y_top, yend = y_top, color = "black", linewidth = 1.1)
  }
  
  # ... (后面的 .run_one_db 函数和主逻辑保持不变, 但需要传递新参数) ...
  uni_syms <- if(!is.null(universe_symbols)) unique(toupper(universe_symbols)) else NULL
  .sym2entrez <- function(symbols){
    symbols <- unique(toupper(symbols))
    if(length(symbols) == 0) return(character(0))
    tb <- suppressMessages(clusterProfiler::bitr(symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db))
    tb <- tb[!duplicated(tb$SYMBOL), ]
    unname(tb$ENTREZID)
  }
  entrez_universe <- if(!is.null(uni_syms)) { x <- .sym2entrez(uni_syms); if(length(x)==0) NULL else x } else NULL
  
  .run_one_db <- function(db_code){
    sig_names <- names(sig_list)
    .lowhigh <- function(s) if(!is.null(color_map) && !is.null(color_map[[s]])) color_map[[s]] else default_low_high
    
    panels <- list(); tables <- list()
    
    if (db_code == "KEGG" || startsWith(db_code, "GO:")) {
      for (s in sig_names){
        g_syms <- unique(na.omit(toupper(sig_list[[s]])))
        title_i <- if(!is.null(title_map) && !is.null(title_map[[s]])) title_map[[s]] else s
        df <- NULL
        
        if (db_code == "KEGG") {
          if (!isTRUE(prefer_online) && !is.null(kegg_offline_gmt) && file.exists(kegg_offline_gmt)) {
            TERM2GENE <- clusterProfiler::read.gmt(kegg_offline_gmt)
            colnames(TERM2GENE) <- c("term","gene")
            eg <- tryCatch({
              enricher(gene = g_syms, TERM2GENE = TERM2GENE,
                       universe = uni_syms, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
                       pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize)
            }, error = function(e) NULL)
            df <- if(!is.null(eg)) as.data.frame(eg) else data.frame()
          } else {
            entrez <- .sym2entrez(g_syms)
            eg <- tryCatch({
              enrichKEGG(gene = entrez, organism = "hsa", keyType = "ncbi-geneid",
                         universe = entrez_universe, pvalueCutoff = pvalueCutoff,
                         pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
                         minGSSize = minGSSize, maxGSSize = maxGSSize)
            }, error = function(e) NULL)
            df <- if(!is.null(eg)) as.data.frame(eg) else data.frame()
          }
        } else {
          ont <- sub("^GO:", "", db_code)
          entrez <- .sym2entrez(g_syms)
          eg <- tryCatch({
            enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                     ont = ont, universe = entrez_universe,
                     pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
                     minGSSize = minGSSize, maxGSSize = maxGSSize, readable = TRUE)
          }, error = function(e) NULL)
          df <- if(!is.null(eg)) as.data.frame(eg) else data.frame()
        }
        
        tables[[s]] <- df
        # 传递新参数 font_size_range
        p <- .panel_plot(df, title_i, .lowhigh(s), size_mult, font_size_range)
        panels[[s]] <- p
        if (isTRUE(save_single) && !is.null(p)) {
          ggsave(file.path(outdir, sprintf("%s_%s_%s.png", prefix, gsub(":","-",db_code), s)),
                 plot = p, width = width_each, height = height_each, dpi = dpi)
          ggsave(file.path(outdir, sprintf("%s_%s_%s.pdf",  prefix, gsub(":","-",db_code), s)),
                 plot = p, width = width_each, height = height_each)
          if (nrow(df)>0) write.csv(df, file.path(outdir, sprintf("%s_%s_%s_top%d.csv", prefix, gsub(":","-",db_code), s, min(top_n,nrow(df)))), row.names = FALSE)
        }
      }
      
    } else {
      coll_code <- switch(db_code,
                          "Reactome"     = "C2:REACTOME",
                          "MSigDB:H"     = "H",
                          "MSigDB:C2"    = "C2",
                          "MSigDB:C6"    = "C6",
                          "MSigDB:C7"    = "C7",
                          "MSigDB:C5:BP" = "C5:BP",
                          "MSigDB:C5:MF" = "C5:MF",
                          "MSigDB:C5:CC" = "C5:CC",
                          stop("不支持的 db：", db_code)
      )
      gmt_path <- .get_gmt(coll_code)
      TERM2GENE <- clusterProfiler::read.gmt(gmt_path)
      colnames(TERM2GENE) <- c("term","gene")
      TERM2GENE$term <- gsub("^REACTOME_", "", TERM2GENE$term)
      uni_here <- if(!is.null(uni_syms)) unique(toupper(uni_syms)) else NULL
      
      for (s in sig_names){
        genes <- unique(na.omit(toupper(sig_list[[s]])))
        title_i <- if(!is.null(title_map) && !is.null(title_map[[s]])) title_map[[s]] else s
        eg <- tryCatch({
          enricher(gene = genes, TERM2GENE = TERM2GENE, universe = uni_here,
                   pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff,
                   pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize)
        }, error = function(e) NULL)
        df <- if(!is.null(eg)) as.data.frame(eg) else data.frame()
        tables[[s]] <- df
        # 传递新参数 font_size_range
        p <- .panel_plot(df, title_i, .lowhigh(s), size_mult, font_size_range)
        panels[[s]] <- p
        if (isTRUE(save_single) && !is.null(p)) {
          ggsave(file.path(outdir, sprintf("%s_%s_%s.png", prefix, gsub(":","-",db_code), s)),
                 plot = p, width = width_each, height = height_each, dpi = dpi)
          ggsave(file.path(outdir, sprintf("%s_%s_%s.pdf",  prefix, gsub(":","-",db_code), s)),
                 plot = p, width = width_each, height = height_each)
          if (nrow(df)>0) write.csv(df, file.path(outdir, sprintf("%s_%s_%s_top%d.csv", prefix, gsub(":","-",db_code), s, min(top_n,nrow(df)))), row.names = FALSE)
        }
      }
    }
    
    panels_valid <- Filter(Negate(is.null), panels)
    combined <- NULL
    if (length(panels_valid) >= 1) {
      combined <- patchwork::wrap_plots(panels_valid, nrow = 1)
      if (isTRUE(save_combined)) {
        ggsave(file.path(outdir, sprintf("%s_%s_combined.png", prefix, gsub(":","-",db_code))),
               plot = combined, width = width_each * length(panels_valid), height = height_each, dpi = dpi)
        ggsave(file.path(outdir, sprintf("%s_%s_combined.pdf",  prefix, gsub(":","-",db_code))),
               plot = combined, width = width_each * length(panels_valid), height = height_each)
      }
    }
    list(panels = panels, combined = combined, tables = tables)
  }
  
  results <- list()
  db <- unique(db)
  for (d in db) results[[d]] <- .run_one_db(d)
  invisible(results)
}

run_hallmark_gsea <- function(
    x,
    output_dir,
    prefix           = "Proj",
    species          = "Homo sapiens",       # set to "Homo sapiens" for human, etc.
    collection       = "H",                  # Hallmark
    gene_cols        = c("Symbol","SYMBOL","Gene","gene","X.ID"),
    rank_by          = c("stat","signed_logp","log2FC"),
    lfc_col          = "log2FoldChange",
    pval_col         = "pvalue",
    padj_col         = "padj",
    minSize          = 10,
    maxSize          = 5000,
    method           = c("multilevel","simple"),
    nperm            = 10000,                # only used if method="simple"
    sig_in_prop      = NULL,                 # keep pathways significant in at least this fraction of groups
    sig_padj         = 0.05,                 # significance cutoff for padj
    min_groups       = NULL,                 # alternative to sig_in_prop: absolute min number of groups significant
    padj_cut         = NULL,                 # optional global padj filter for plotting
    nes_cap          = 3,
    theme_obj        = NULL,
    width_h          = NULL, height_h = NULL,   # horizontal plot size
    width_v          = NULL, height_v = NULL,   # vertical plot size
    units            = "mm",
    dpi              = 600,
    cluster_paths    = TRUE,                    # cluster pathways?
    cluster_groups   = FALSE,                   # cluster groups?
    cluster_metric   = c("NES","signed_log10p"),
    dist_method      = "euclidean",
    hclust_method    = "ward.D2",
    path_order_by        = c("cluster","absNES","NES","signed_log10p"),
    path_order_stat      = c("mean","median","max_abs"),
    path_order_group     = NULL,
    path_order_decreasing= TRUE,
    pretty_names     = TRUE,
    save_plots       = TRUE,
    return_plots     = TRUE,
    grey_non_sig     = FALSE,                  # gray-out non-significant points (used only when non-sig are kept)
    color_sig_padj   = NULL,                   # threshold for coloring as significant; defaults to sig_padj
    non_sig_color    = "grey80",
    non_sig_alpha    = 0.55,
    plot_only_sig    = NULL                    # NEW: auto-hide non-sig points for <=2 groups (TRUE/FALSE/NULL=auto)
){
  suppressPackageStartupMessages({
    library(dplyr); library(tibble); library(purrr)
    library(msigdbr); library(fgsea)
    library(ggplot2); library(forcats); library(stringr); library(tidyr)
  })
  
  # ----- theme -----
  .default_theme <- theme_bw(base_size = 12) + theme(
    panel.grid = element_blank(),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
  use_theme <- if (!is.null(theme_obj)) theme_obj else
    if (exists("prism_theme", inherits = TRUE)) get("prism_theme", inherits = TRUE) else
      .default_theme
  
  # ----- output directory -----
  gsea_dir <- file.path(output_dir, paste0(prefix, "_GSEA_Hallmark"))
  dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ----- MSigDB Hallmark gene sets -----
  msig_tbl <- msigdbr::msigdbr(species = species, category = collection) %>%
    dplyr::select(gs_name, gene_symbol) %>% distinct()
  if (nrow(msig_tbl) == 0) stop("No Hallmark gene sets found; check species/collection.")
  hallmark_sets <- split(msig_tbl$gene_symbol, msig_tbl$gs_name)
  
  # ----- pathway pretty labels -----
  .pretty_pathway <- function(x){
    y <- gsub("^HALLMARK_", "", x)
    y <- gsub("_", " ", y)
    repl <- list(
      "TNFA"="TNF-α", "NFKB"="NF-κB", "TGF BETA"="TGF-β",
      "IL6"="IL-6", "IL2"="IL-2", "IL 2"="IL-2", "IL 6"="IL-6",
      "G2M"="G2/M", "MYC"="MYC", "KRAS"="KRAS", "E2F"="E2F",
      "UV RESPONSE UP"="UV response (up)",
      "UV RESPONSE DN"="UV response (down)"
    )
    for(k in names(repl)) y <- gsub(k, repl[[k]], y, fixed = TRUE)
    y <- tools::toTitleCase(tolower(y))
    caps <- c("MYC","KRAS","E2F","DNA","RNA","TNF-α","NF-κB","G2/M","TGF-β","IL-2","IL-6")
    for (c in caps) y <- gsub(tools::toTitleCase(tolower(c)), c, y, fixed = TRUE)
    y
  }
  
  # ----- choose gene column & ranks -----
  rank_by        <- match.arg(rank_by)
  method         <- match.arg(method)
  cluster_metric <- match.arg(cluster_metric)
  path_order_by  <- match.arg(path_order_by)
  path_order_stat<- match.arg(path_order_stat)
  
  .pick_gene_col <- function(df, gene_cols) {
    gc <- gene_cols[gene_cols %in% names(df)]
    if (length(gc) == 0) stop("No gene column found (tried: ", paste(gene_cols, collapse = "/"), ").")
    gc[1]
  }
  
  .make_ranks <- function(df, gene_col){
    # ensure padj exists if only p-value is present
    if (!(padj_col %in% names(df)) && (pval_col %in% names(df))) {
      df[[padj_col]] <- p.adjust(df[[pval_col]], method = "BH")
    }
    if (rank_by == "stat" && ("stat" %in% names(df))) {
      df_rank <- df %>% filter(!is.na(.data[[gene_col]]), !is.na(stat)) %>%
        group_by(.data[[gene_col]]) %>% slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
        ungroup() %>% dplyr::select(!!gene_col, stat)
      sc <- setNames(df_rank$stat, df_rank[[gene_col]])
    } else if (rank_by == "log2FC" && (lfc_col %in% names(df))) {
      df_rank <- df %>% filter(!is.na(.data[[gene_col]]), !is.na(.data[[lfc_col]])) %>%
        group_by(.data[[gene_col]]) %>% slice_max(order_by = abs(.data[[lfc_col]]), n = 1, with_ties = FALSE) %>%
        ungroup() %>% dplyr::select(!!gene_col, !!lfc_col)
      sc <- setNames(df_rank[[lfc_col]], df_rank[[gene_col]])
    } else {
      if (!(pval_col %in% names(df)) || !(lfc_col %in% names(df)))
        stop("Need 'stat' OR both ", lfc_col, " + ", pval_col, " to construct ranking scores.")
      df_rank <- df %>% filter(!is.na(.data[[gene_col]]), !is.na(.data[[pval_col]]), !is.na(.data[[lfc_col]])) %>%
        mutate(signed_logp = sign(.data[[lfc_col]]) * -log10(pmax(.data[[pval_col]], .Machine$double.xmin))) %>%
        group_by(.data[[gene_col]]) %>% slice_max(order_by = abs(signed_logp), n = 1, with_ties = FALSE) %>%
        ungroup() %>% dplyr::select(!!gene_col, signed_logp)
      sc <- setNames(df_rank$signed_logp, df_rank[[gene_col]])
    }
    sort(sc, decreasing = TRUE)
  }
  
  # ----- run fgsea for one group -----
  # ----- run fgsea for one group -----
  .run_one <- function(df, nm){
    gene_col <- .pick_gene_col(df, gene_cols)
    ranks <- .make_ranks(df, gene_col)
    overlap <- sum(names(ranks) %in% unique(unlist(hallmark_sets)))
    cov_pct <- round(100 * overlap / max(1, length(ranks)), 2)
    message(sprintf("[GSEA:%s] ranks=%s, overlap=%s (%.2f%%)", nm, length(ranks), overlap, cov_pct))
    
    set.seed(1)
    if (method == "multilevel") {
      fg <- fgsea(pathways = hallmark_sets, stats = ranks,
                  minSize = minSize, maxSize = maxSize)
    } else {
      fg <- fgsea(pathways = hallmark_sets, stats = ranks,
                  minSize = minSize, maxSize = maxSize, nperm = nperm)
    }
    
    # 修改点 1：使用 transmute 时保留 leadingEdge，并将其转换为字符串格式
    as_tibble(fg) %>%
      transmute(
        Group = nm, 
        Pathway = pathway, 
        NES, pval, padj, size,
        leadingEdge = map_chr(leadingEdge, ~paste(.x, collapse = ";")) # 将基因列表转为 "Gene1;Gene2"
      ) %>%
      arrange(padj, desc(abs(NES)))
  }
  
  # ----- support data.frame or named list -----
  if (is.data.frame(x)) {
    gsea_tbls <- list(ONE = .run_one(x, "ONE"))
  } else if (is.list(x)) {
    if (is.null(names(x))) stop("If 'x' is a list, it must be a *named* list (names are group labels).")
    gsea_tbls <- purrr::imap(x, ~ .run_one(.x, .y))
  } else stop("'x' must be a data.frame or a named list of data.frames.")
  
  # ----- write per-group tables & merge -----
  purrr::iwalk(gsea_tbls, ~ write.csv(.x, file.path(gsea_dir, sprintf("%s_Hallmark_GSEA_%s.csv", prefix, .y)), row.names = FALSE))
  gsea_all <- dplyr::bind_rows(gsea_tbls) %>% dplyr::filter(!is.na(padj))
  if (nrow(gsea_all) == 0) stop("fgsea returned empty results; check ranks and gene set coverage.")
  
  # ----- group count for auto logic -----
  n_groups <- length(unique(gsea_all$Group))
  
  # ----- pathway-level filtering for multi-group context -----
  # Keep pathways that are significant (padj <= sig_padj) in at least a proportion or absolute number of groups.
  need_k <- if (!is.null(min_groups)) min_groups else
    if (!is.null(sig_in_prop)) ceiling(sig_in_prop * n_groups) else NA_integer_
  
  if (!is.na(need_k) && n_groups > 1) {
    keep_paths <- gsea_all %>%
      mutate(sig = padj <= sig_padj) %>%
      group_by(Pathway) %>%
      summarise(n_sig = sum(sig, na.rm = TRUE), .groups = "drop") %>%
      filter(n_sig >= need_k) %>%
      pull(Pathway)
    
    plot_df <- if (length(keep_paths)) {
      gsea_all %>% filter(Pathway %in% keep_paths)
    } else {
      message(sprintf("No pathways met padj<=%.3g in >=%d/%d groups; falling back to all pathways.", sig_padj, need_k, n_groups))
      gsea_all
    }
  } else {
    # Optional global padj filter when not using multi-group proportion
    if (!is.null(padj_cut) && padj_cut < 1) {
      keep_paths <- gsea_all %>% filter(padj <= padj_cut) %>% pull(Pathway) %>% unique()
      plot_df <- if (length(keep_paths)) gsea_all %>% filter(Pathway %in% keep_paths) else gsea_all
      if (length(keep_paths) == 0) message("No significant pathways under 'padj_cut'; plotting all.")
    } else {
      plot_df <- gsea_all
    }
  }
  
  # ----- NEW: auto hide non-significant points when few groups -----
  # If plot_only_sig = NULL (auto): TRUE when n_groups <= 2; FALSE otherwise.
  if (is.null(plot_only_sig)) plot_only_sig <- (n_groups <= 2)
  if (isTRUE(plot_only_sig)) {
    plot_df <- plot_df %>% dplyr::filter(is.finite(padj) & !is.na(padj) & padj <= sig_padj)
    if (nrow(plot_df) == 0) {
      message(sprintf("No significant pathways to plot (padj <= %.3g). Returning tables only.", sig_padj))
    }
  }
  
  # ----- plot helper columns -----
  plot_df <- plot_df %>%
    mutate(Log10_p = -log10(pmax(padj, .Machine$double.xmin)),
           NES_cap = pmax(pmin(NES, nes_cap), -nes_cap))
  if (is.null(color_sig_padj)) color_sig_padj <- sig_padj
  plot_df <- plot_df %>% mutate(is_sig = ifelse(is.finite(padj) & !is.na(padj), padj <= color_sig_padj, FALSE))
  
  # ----- pretty pathway label mapping -----
  if (pretty_names) {
    path_map <- tibble(Pathway = unique(plot_df$Pathway)) %>%
      mutate(Pathway_label = .pretty_pathway(Pathway))
  } else {
    path_map <- tibble(Pathway = unique(plot_df$Pathway),
                       Pathway_label = Pathway)
  }
  plot_df <- plot_df %>% left_join(path_map, by = "Pathway")
  write.csv(path_map, file.path(gsea_dir, paste0(prefix, "_Pathway_Name_Mapping.csv")), row.names = FALSE)
  
  # ----- metric matrix for ordering/clustering -----
  mat_metric <- if (cluster_metric == "NES") {
    plot_df %>% dplyr::select(Pathway, Group, NES) %>%
      tidyr::pivot_wider(names_from = Group, values_from = NES, values_fill = 0)
  } else {
    plot_df %>% mutate(signed = sign(NES) * Log10_p) %>%
      dplyr::select(Pathway, Group, signed) %>%
      tidyr::pivot_wider(names_from = Group, values_from = signed, values_fill = 0)
  }
  mat <- as.matrix(mat_metric[,-1, drop = FALSE])
  rownames(mat) <- mat_metric[[1]]
  
  # ----- default orders (fallbacks) -----
  path_order <- unique(plot_df$Pathway)
  grp_order  <- unique(plot_df$Group)
  
  # ----- pathway ordering -----
  if (path_order_by == "cluster") {
    if (cluster_paths && nrow(mat) >= 2) {
      d  <- dist(scale(mat), method = dist_method)
      hc <- hclust(d, method = hclust_method)
      path_order <- rownames(mat)[hc$order]
    }
  } else {
    metric_df <- plot_df %>%
      mutate(metric = dplyr::case_when(
        path_order_by == "absNES"        ~ abs(NES),
        path_order_by == "NES"           ~ NES,
        path_order_by == "signed_log10p" ~ sign(NES) * Log10_p
      )) %>%
      dplyr::select(Pathway, Group, metric)
    
    if (!is.null(path_order_group)) {
      metric_x <- metric_df %>%
        dplyr::filter(Group == path_order_group) %>%
        dplyr::group_by(Pathway) %>%
        dplyr::summarise(val = dplyr::first(metric), .groups = "drop")
    } else {
      metric_x <- metric_df %>%
        dplyr::group_by(Pathway) %>%
        dplyr::summarise(
          val = dplyr::case_when(
            path_order_stat == "max_abs" ~ max(abs(metric), na.rm = TRUE),
            path_order_stat == "mean"    ~ mean(metric, na.rm = TRUE),
            TRUE                         ~ stats::median(metric, na.rm = TRUE)
          ),
          .groups = "drop"
        ) %>%
        dplyr::arrange(dplyr::desc(val))
      
    }
    
    metric_x <- if (isTRUE(path_order_decreasing)) {
      dplyr::arrange(metric_x, dplyr::desc(val))
    } else {
      dplyr::arrange(metric_x, val)
    }
    path_order <- metric_x$Pathway
  }
  
  # ----- group clustering (optional) -----
  if (cluster_groups && ncol(mat) >= 2) {
    d2  <- dist(scale(t(mat)), method = dist_method)
    hc2 <- hclust(d2, method = hclust_method)
    grp_order <- colnames(mat)[hc2$order]
  }
  
  # ----- factor levels for labels -----
  label_order <- path_map %>%
    mutate(idx = match(Pathway, path_order)) %>%
    arrange(idx) %>% pull(Pathway_label)
  
  plot_df <- plot_df %>%
    mutate(
      Pathway_label = factor(Pathway_label, levels = label_order),
      Group         = factor(Group, levels = grp_order)
    )
  # lock the order once (use your current factor order)
  group_lvls <- levels(droplevels(plot_df$Group))
  # if you want the vertical axis to show top-to-bottom reversed, use rev(group_lvls)
  # group_lvls_rev <- rev(group_lvls)
  
  # ----- plotting -----
  if (grey_non_sig) {
    # gray non-sig and color sig (only effective if non-sig points were not filtered out)
    p_h <- ggplot(plot_df, aes(x = Pathway_label, y = Group)) +
      geom_point(
        data = subset(plot_df, !is_sig),
        aes(size = Log10_p),
        color = non_sig_color, alpha = non_sig_alpha, na.rm = TRUE
      ) +
      geom_point(
        data = subset(plot_df,  is_sig),
        aes(size = Log10_p, color = NES_cap),
        alpha = 0.9, na.rm = TRUE
      ) +
      scale_size_continuous(name = "-log10(padj)", range = c(1.8, 8)) +
      scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0, name = "NES") +
      use_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9)) +
      labs(title = paste0(prefix, " - Hallmark GSEA (horizontal)"), x = NULL, y = NULL) +
      scale_x_discrete(limits = levels(plot_df$Pathway_label), drop = FALSE)  +scale_y_discrete(limits = group_lvls, drop = FALSE)   # <-- lock Group order
    
    p_v <- ggplot(plot_df, aes(x = Group, y = Pathway_label)) +
      geom_point(
        data = subset(plot_df, !is_sig),
        aes(size = Log10_p),
        color = non_sig_color, alpha = non_sig_alpha, na.rm = TRUE
      ) +
      geom_point(
        data = subset(plot_df,  is_sig),
        aes(size = Log10_p, color = NES_cap),
        alpha = 0.9, na.rm = TRUE
      ) +
      scale_size_continuous(name = "-log10(padj)", range = c(1.8, 8)) +
      scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0, name = "NES") +
      use_theme +
      labs(title = paste0(prefix, " - Hallmark GSEA (vertical)"), x = NULL, y = NULL) +
      scale_y_discrete(limits = rev(levels(plot_df$Pathway_label)), drop = FALSE) +scale_x_discrete(limits = group_lvls, drop = FALSE)   # <-- lock Group order
  } else {
    p_h <- ggplot(plot_df, aes(x = Pathway_label, y = Group)) +
      geom_point(aes(size = Log10_p, color = NES_cap), alpha = 0.9, na.rm = TRUE) +
      scale_size_continuous(name = "-log10(padj)", range = c(1.8, 8)) +
      scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0, name = "NES") +
      use_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9)) +
      labs(title = paste0(prefix, " - Hallmark GSEA (horizontal)"), x = NULL, y = NULL) +
      scale_x_discrete(limits = levels(plot_df$Pathway_label), drop = FALSE) +scale_y_discrete(limits = group_lvls, drop = FALSE)   # <-- lock Group order
    
    p_v <- ggplot(plot_df, aes(x = Group, y = Pathway_label)) +
      geom_point(aes(size = Log10_p, color = NES_cap), alpha = 0.9, na.rm = TRUE) +
      scale_size_continuous(name = "-log10(padj)", range = c(1.8, 8)) +
      scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0, name = "NES") +
      use_theme +
      labs(title = paste0(prefix, " - Hallmark GSEA (vertical)"), x = NULL, y = NULL) +
      scale_y_discrete(limits = rev(levels(plot_df$Pathway_label)), drop = FALSE) +scale_x_discrete(limits = group_lvls, drop = FALSE)   # <-- lock Group order
  }
  
  # ----- save plots (skip if empty) -----
  if (save_plots && nrow(plot_df) > 0) {
    if (is.null(width_h))  width_h  <- max(140, 5 + 0.22 * length(levels(plot_df$Pathway_label)) * 10)
    if (is.null(height_h)) height_h <- max( 60, 5 + 0.50 * length(levels(plot_df$Group))          * 10)
    if (is.null(width_v))  width_v  <- max(150, 5 + 0.55 * length(levels(plot_df$Group))          * 12)
    if (is.null(height_v)) height_v <- max(120, 5 + 0.20 * length(levels(plot_df$Pathway_label))  * 16)
    
    ggsave(file.path(gsea_dir, paste0(prefix, "_GSEA_Hallmark_Bubble_Horizontal.pdf")),
           p_h, width = width_h, height = height_h, units = units, dpi = dpi)
    ggsave(file.path(gsea_dir, paste0(prefix, "_GSEA_Hallmark_Bubble_Horizontal.png")),
           p_h, width = width_h, height = height_h, units = units, dpi = dpi)
    ggsave(file.path(gsea_dir, paste0(prefix, "_GSEA_Hallmark_Bubble_Vertical.pdf")),
           p_v, width = width_v, height = height_v, units = units, dpi = dpi)
    ggsave(file.path(gsea_dir, paste0(prefix, "_GSEA_Hallmark_Bubble_Vertical.png")),
           p_v, width = width_v, height = height_v, units = units, dpi = dpi)
  }
  
  # ----- export tables -----
  # 修改点 2：在导出 PlotUsed 表格时，也带上 leadingEdge
  write.csv(gsea_all, file.path(gsea_dir, paste0(prefix, "_Hallmark_GSEA_All_Groups_full.csv")), row.names = FALSE)
  
  # 修改点 3：确保 plot_df 也包含这一列（plot_df 是由 gsea_all 过滤得到的，所以会自动包含）
  # 只是在导出 PlotUsed.csv 时显式写出它
  if ("leadingEdge" %in% names(plot_df)) {
    write.csv(plot_df %>% dplyr::select(Group, Pathway, Pathway_label, NES, padj, Log10_p, NES_cap, size, leadingEdge),
              file.path(gsea_dir, paste0(prefix, "_Hallmark_GSEA_PlotUsed.csv")), row.names = FALSE)
  }

  out <- list(
    per_group   = gsea_tbls,
    combined    = gsea_all,
    plot_df     = plot_df,
    outdir      = gsea_dir,
    plots       = if (return_plots && nrow(plot_df) > 0) list(horizontal = p_h, vertical = p_v) else NULL
  )
  return(out)
}
library(limma)
library(ggplot2)
library(reshape2)
library(ggfortify)


correct_systematic_shift <- function(expr_mat, 
                                     meta_data, 
                                     group_col = "treatment", 
                                     method = "quantile",
                                     outdir = "QC_Correction", 
                                     prefix = "BatchFix") {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # --- 内部绘图函数 ---
  plot_dist <- function(mat, title_suffix) {
    # PCA
    pca_res <- prcomp(t(mat[apply(mat, 1, var) > 0, ]), scale. = TRUE)
    p_pca <- autoplot(pca_res, data = meta_data, colour = group_col, size = 3) +
      theme_bw() +
      labs(title = paste0("PCA: ", title_suffix))
    
    # Boxplot (抽样 5000 基因以提速)
    set.seed(123)
    if(nrow(mat) > 5000) mat_sub <- mat[sample(nrow(mat), 5000), ] else mat_sub <- mat
    
    df_long <- melt(mat_sub)
    colnames(df_long) <- c("Gene", "Sample", "Expression")
    df_long <- merge(df_long, meta_data, by.x = "Sample", by.y = "row.names")
    
    p_box <- ggplot(df_long, aes_string(x = "Sample", y = "Expression", fill = group_col)) +
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.1) +
      theme_bw() +
      labs(title = paste0("Distributions: ", title_suffix), y = "Expression") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
    
    return(list(pca = p_pca, box = p_box))
  }
  
  message(">>> 1. 绘制矫正前 (Pre-Correction) 图像...")
  plots_pre <- plot_dist(expr_mat, "Raw / Pre-Correction")
  ggsave(file.path(outdir, paste0(prefix, "_01_Pre_PCA.pdf")), plots_pre$pca, width = 6, height = 5)
  ggsave(file.path(outdir, paste0(prefix, "_01_Pre_Boxplot.pdf")), plots_pre$box, width = 10, height = 6)
  
  message(paste0(">>> 2. 执行归一化矫正 (Method: ", method, ")..."))
  # 核心：使用 limma 的 normalizeBetweenArrays 进行分位数归一化
  expr_corrected <- limma::normalizeBetweenArrays(expr_mat, method = method)
  
  message(">>> 3. 绘制矫正后 (Post-Correction) 图像...")
  plots_post <- plot_dist(expr_corrected, "Quantile Normalized")
  ggsave(file.path(outdir, paste0(prefix, "_02_Post_PCA.pdf")), plots_post$pca, width = 6, height = 5)
  ggsave(file.path(outdir, paste0(prefix, "_02_Post_Boxplot.pdf")), plots_post$box, width = 10, height = 6)
  
  message(paste0(">>> 完成。结果保存在: ", outdir))
  return(expr_corrected)
}
de_universal <- function(
    expr,
    meta,
    sample_id_col = NULL,               
    group_col = "treatment",
    group_levels = c("control", "drug"),
    pair_col = NULL,                    
    require_pairs = TRUE,
    choose_onepair = TRUE,
    drop_all_zero_samples = TRUE,
    
    method = c("auto", "limma_block", "limma_fixed", "voom", "deseq2"),
    data_type = c("auto", "counts", "log"),  
    log_offset = 1,                   
    transform = c("auto", "none", "log2p1"),
    
    batch_cols = NULL,                  
    batch_mode = c("design", "removeBatchEffect"), 
    extra_covariates = NULL,            
    
    contrast = NULL,                    
    outdir = NULL,
    prefix = "DE",
    
    do_hallmark = FALSE,
    hallmark_args = list(),             
    
    verbose = TRUE
) {
  method      <- match.arg(method)
  data_type   <- match.arg(data_type)
  transform   <- match.arg(transform)
  batch_mode  <- match.arg(batch_mode)
  
  # ---- packages ----
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
  })
  
  # ---- align meta rownames to expr colnames ----
  meta <- as.data.frame(meta)
  if (!is.null(sample_id_col)) {
    stopifnot(sample_id_col %in% colnames(meta))
    rownames(meta) <- as.character(meta[[sample_id_col]])
  }
  stopifnot(!is.null(colnames(expr)))
  stopifnot(!is.null(rownames(meta)))
  
  # enforce same samples
  common_samples <- intersect(colnames(expr), rownames(meta))
  if (length(common_samples) < 2) stop("Too few overlapping samples between expr and meta.")
  expr <- expr[, common_samples, drop = FALSE]
  meta <- meta[common_samples, , drop = FALSE]
  
  # order meta to expr
  meta <- meta[colnames(expr), , drop = FALSE]
  stopifnot(identical(colnames(expr), rownames(meta)))
  
  # ---- basic checks ----
  stopifnot(group_col %in% colnames(meta))
  # Enforce factor levels
  meta[[group_col]] <- factor(meta[[group_col]], levels = group_levels)
  
  if (!is.null(pair_col)) {
    stopifnot(pair_col %in% colnames(meta))
    meta[[pair_col]] <- factor(meta[[pair_col]])
  }
  
  # ---- drop all-zero/invalid samples ----
  if (isTRUE(drop_all_zero_samples)) {
    is_all_zero <- apply(expr, 2, function(x) all(is.na(x)) || sum(abs(x), na.rm = TRUE) == 0)
    if (any(is_all_zero)) {
      if (verbose) message("Dropping all-zero/invalid samples: ", paste(colnames(expr)[is_all_zero], collapse = ", "))
      expr <- expr[, !is_all_zero, drop = FALSE]
      meta <- meta[!is_all_zero, , drop = FALSE]
    }
  }
  
  # ---- require paired cell lines (keep only those with both groups) ----
  if (isTRUE(require_pairs) && !is.null(pair_col)) {
    pair_table <- meta %>%
      mutate(sample = rownames(meta)) %>%
      group_by(.data[[pair_col]]) %>%
      summarise(
        n_g1 = sum(.data[[group_col]] == group_levels[1]),
        n_g2 = sum(.data[[group_col]] == group_levels[2]),
        .groups = "drop"
      )
    
    paired_levels <- pair_table %>%
      filter(n_g1 > 0, n_g2 > 0) %>%
      pull(.data[[pair_col]]) %>%
      as.character()
    
    if (verbose) message("Paired ", pair_col, " kept (", length(paired_levels), "): ", paste(paired_levels, collapse = ", "))
    
    keep_samples <- rownames(meta)[as.character(meta[[pair_col]]) %in% paired_levels]
    expr <- expr[, keep_samples, drop = FALSE]
    meta <- meta[keep_samples, , drop = FALSE]
  }
  
  # ---- optionally keep exactly 1 pair per block ----
  if (isTRUE(choose_onepair) && !is.null(pair_col)) {
    md <- meta %>% mutate(sample = rownames(meta))
    md_keep <- md %>%
      group_by(.data[[pair_col]]) %>%
      slice(which(.data[[group_col]] == group_levels[1])[1]) %>%
      bind_rows(
        md %>% group_by(.data[[pair_col]]) %>% slice(which(.data[[group_col]] == group_levels[2])[1])
      ) %>%
      ungroup()
    
    keep_samples <- md_keep$sample
    expr <- expr[, keep_samples, drop = FALSE]
    meta <- meta[keep_samples, , drop = FALSE]
    meta <- meta[colnames(expr), , drop = FALSE]
  }
  
  # ---- infer data_type if auto ----
  if (data_type == "auto") {
    x <- as.numeric(expr[1:min(nrow(expr), 200), 1:min(ncol(expr), 30)])
    x <- x[is.finite(x)]
    has_neg <- any(x < 0)
    frac_int <- mean(abs(x - round(x)) < 1e-8)
    if (!has_neg && frac_int > 0.9) data_type <- "counts" else data_type <- "log"
    if (verbose) message("Inferred data_type = ", data_type)
  }
  
  # ---- choose method if auto ----
  if (method == "auto") {
    method <- if (data_type == "counts") "deseq2" else "limma_block"
    if (verbose) message("Auto-selected method = ", method)
  }
  
  # ---- transform if needed ----
  expr_use <- as.matrix(expr)
  storage.mode(expr_use) <- "numeric"
  
  if (transform == "auto") {
    transform <- if (data_type == "counts") "none" else "none"
  }
  if (transform == "log2p1") {
    expr_use <- log2(expr_use + log_offset)
  }
  
  # ---- batch handling (recommended: include in design) ----
  batch_terms <- character(0)
  if (!is.null(batch_cols)) {
    stopifnot(all(batch_cols %in% colnames(meta)))
    for (bc in batch_cols) meta[[bc]] <- factor(meta[[bc]])
    batch_terms <- batch_cols
  }
  if (!is.null(extra_covariates)) {
    stopifnot(all(extra_covariates %in% colnames(meta)))
  }
  
  # ============================================================
  # Run DE
  # ============================================================
  res_out <- NULL
  fit_out <- NULL
  used <- list(expr = expr_use, meta = meta)
  
  if (method %in% c("limma_block", "limma_fixed", "voom")) {
    suppressPackageStartupMessages(library(limma))
  }
  if (method == "voom") {
    suppressPackageStartupMessages(library(edgeR))
  }
  if (method == "deseq2") {
    suppressPackageStartupMessages(library(DESeq2))
  }
  
  if (method == "limma_block") {
    stopifnot(data_type == "log")
    stopifnot(!is.null(pair_col))
    
    # design: group + optional batch/covariates
    design_terms <- c(paste0("0 + ", group_col), batch_terms, extra_covariates)
    design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
    design <- model.matrix(design_formula, data = meta)
    
    # FIX: Ensure column names are valid R names to avoid matching issues later
    colnames(design) <- make.names(colnames(design))
    
    # optional removeBatchEffect
    if (!is.null(batch_cols) && batch_mode == "removeBatchEffect") {
      expr_use <- removeBatchEffect(expr_use, batch = meta[[batch_cols[1]]])
    }
    
    block <- meta[[pair_col]]
    corfit <- duplicateCorrelation(expr_use, design, block = block)
    
    fit <- lmFit(expr_use, design, block = block, correlation = corfit$consensus)
    
    # contrast
    g1 <- group_levels[1]; g2 <- group_levels[2]
    if (is.null(contrast)) contrast <- c(g2, g1)
    
    # --- FIXED COLUMN DETECTION LOGIC ---
    # Instead of regex matching, construct the expected column names
    # When using ~0 + group_col, model.matrix creates columns like "group_colLevel"
    # We also sanitize them because we ran make.names() on design colnames
    
    expected_col_g2 <- make.names(paste0(group_col, contrast[1]))
    expected_col_g1 <- make.names(paste0(group_col, contrast[2]))
    
    # Check if they exist in design
    if (!expected_col_g2 %in% colnames(design)) stop("Column '", expected_col_g2, "' not found in design matrix.")
    if (!expected_col_g1 %in% colnames(design)) stop("Column '", expected_col_g1, "' not found in design matrix.")
    
    contrast_formula <- paste(expected_col_g2, "-", expected_col_g1)
    if (verbose) message("Contrasting: ", contrast_formula)
    
    cm <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
    
    tt <- topTable(fit2, number = Inf, sort.by = "P")
    tt$Symbol <- rownames(tt)
    tt <- tt %>%
      as.data.frame() %>%
      transmute(
        Symbol,
        log2FoldChange = logFC,
        AveExpr,
        stat = t,
        pvalue = P.Value,
        padj = adj.P.Val
      ) %>% arrange(padj, desc(abs(log2FoldChange)))
    
    res_out <- tt
    fit_out <- list(fit = fit2, corfit = corfit)
    
  } else if (method == "limma_fixed") {
    stopifnot(data_type == "log")
    
    fixed_terms <- character(0)
    if (!is.null(pair_col)) fixed_terms <- c(fixed_terms, paste0("factor(", pair_col, ")"))
    design_terms <- c(paste0("0 + ", group_col), fixed_terms, batch_terms, extra_covariates)
    design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
    design <- model.matrix(design_formula, data = meta)
    colnames(design) <- make.names(colnames(design))
    
    fit <- lmFit(expr_use, design)
    
    g1 <- group_levels[1]; g2 <- group_levels[2]
    if (is.null(contrast)) contrast <- c(g2, g1)
    
    # Same fixed logic
    expected_col_g2 <- make.names(paste0(group_col, contrast[1]))
    expected_col_g1 <- make.names(paste0(group_col, contrast[2]))
    
    if (!expected_col_g2 %in% colnames(design)) stop("Column '", expected_col_g2, "' not found in design matrix.")
    if (!expected_col_g1 %in% colnames(design)) stop("Column '", expected_col_g1, "' not found in design matrix.")
    
    contrast_formula <- paste(expected_col_g2, "-", expected_col_g1)
    cm <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
    
    tt <- topTable(fit2, number = Inf, sort.by = "P")
    tt$Symbol <- rownames(tt)
    tt <- tt %>%
      as.data.frame() %>%
      transmute(
        Symbol,
        log2FoldChange = logFC,
        AveExpr,
        stat = t,
        pvalue = P.Value,
        padj = adj.P.Val
      ) %>% arrange(padj, desc(abs(log2FoldChange)))
    
    res_out <- tt
    fit_out <- list(fit = fit2)
    
  } else if (method == "voom") {
    stopifnot(data_type == "counts")
    
    y <- edgeR::DGEList(counts = round(expr_use))
    y <- edgeR::calcNormFactors(y)
    
    design_terms <- c(paste0("0 + ", group_col))
    if (!is.null(pair_col)) design_terms <- c(design_terms, paste0("factor(", pair_col, ")"))
    design_terms <- c(design_terms, batch_terms, extra_covariates)
    design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
    
    design <- model.matrix(design_formula, data = meta)
    colnames(design) <- make.names(colnames(design))
    
    v <- limma::voom(y, design, plot = FALSE)
    fit <- lmFit(v, design)
    
    g1 <- group_levels[1]; g2 <- group_levels[2]
    if (is.null(contrast)) contrast <- c(g2, g1)
    
    expected_col_g2 <- make.names(paste0(group_col, contrast[1]))
    expected_col_g1 <- make.names(paste0(group_col, contrast[2]))
    
    contrast_formula <- paste(expected_col_g2, "-", expected_col_g1)
    cm <- makeContrasts(contrasts = contrast_formula, levels = design)
    
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
    
    tt <- topTable(fit2, number = Inf, sort.by = "P")
    tt$Symbol <- rownames(tt)
    tt <- tt %>%
      as.data.frame() %>%
      transmute(
        Symbol,
        log2FoldChange = logFC,
        AveExpr,
        stat = t,
        pvalue = P.Value,
        padj = adj.P.Val
      ) %>% arrange(padj, desc(abs(log2FoldChange)))
    
    res_out <- tt
    fit_out <- list(fit = fit2, voom = v)
    
  } else if (method == "deseq2") {
    stopifnot(data_type == "counts")
    
    counts <- round(expr_use)
    mode(counts) <- "integer"
    
    terms <- character(0)
    if (!is.null(batch_terms)) terms <- c(terms, batch_terms)
    if (!is.null(pair_col))    terms <- c(terms, pair_col)
    terms <- c(terms, group_col)
    design_formula <- as.formula(paste("~", paste(terms, collapse = " + ")))
    
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = meta, design = design_formula)
    dds <- DESeq2::DESeq(dds)
    
    if (is.null(contrast)) {
      contrast <- c(group_col, group_levels[2], group_levels[1])
    } else {
      if (length(contrast) == 2) contrast <- c(group_col, contrast[1], contrast[2])
    }
    res <- as.data.frame(DESeq2::results(dds, contrast = contrast))
    res$Symbol <- rownames(res)
    
    res_out <- res %>%
      transmute(
        Symbol,
        log2FoldChange = log2FoldChange,
        stat = stat,
        pvalue = pvalue,
        padj = padj
      ) %>% arrange(padj, desc(abs(log2FoldChange)))
    
    fit_out <- list(dds = dds)
    
  } else {
    stop("Unsupported method: ", method)
  }
  
  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    
    # 1. Write CSV
    write.csv(res_out, file.path(outdir, paste0(prefix, "_DE.csv")), row.names = FALSE)
    write.csv(meta, file.path(outdir, paste0(prefix, "_MetaUsed.csv")), row.names = TRUE)
    
    # 2. Write Excel (if openxlsx installed)
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      openxlsx::write.xlsx(res_out, file.path(outdir, paste0(prefix, "_DE.xlsx")))
      # For metadata, we preserve rownames by setting rowNames = TRUE
      openxlsx::write.xlsx(meta, file.path(outdir, paste0(prefix, "_MetaUsed.xlsx")), rowNames = TRUE)
      if (verbose) message("Saved Excel files to ", outdir)
    } else {
      warning("Package 'openxlsx' is not installed. Skipping .xlsx export.")
    }
  }
  gsea_out <- NULL
  if (isTRUE(do_hallmark)) {
    if (!exists("run_hallmark_gsea", mode = "function")) {
      warning("run_hallmark_gsea() not found in environment. Did you source RNA_seq.R?")
    } else {
      ha <- hallmark_args
      ha$x <- res_out
      ha$output_dir <- ifelse(is.null(outdir), getwd(), outdir)
      ha$prefix <- prefix
      gsea_out <- do.call(run_hallmark_gsea, ha)
    }
  }
  
  invisible(list(
    res = res_out,
    fit = fit_out,
    used = used,
    gsea = gsea_out
  ))
}

library(ggplot2)
library(ggsci)
library(pheatmap)
library(factoextra)
library(ConsensusClusterPlus)
library(NMF)

#' 通用无监督聚类分析函数 (全套图表版)
#'
#' @param expr_mat 表达矩阵 (行=基因, 列=样本)。建议使用 Log2(FPKM+1) 或 VST 数据。
#' @param methods 聚类方法。可选: "kmeans", "consensus", "nmf", "all"。默认 "all"。
#' @param k_range 尝试的 K 值范围，用于评估 (例如 2:6)。
#' @param best_k 指定的最佳 K 值 (用于最终出图)。若为 NULL，则只做评估，不产出最终结果。
#' @param n_top_genes 使用多少个高变基因进行聚类。默认 2000。
#' @param outdir 输出目录。
#' @param seed 随机种子。
#' @param nmf_nrun NMF 重复次数。测试建议 10，正式建议 30。
#'
#' @return 返回包含分组信息的列表
library(ggplot2)
library(ggsci)
library(pheatmap)
library(factoextra)
library(ConsensusClusterPlus)
library(NMF)
library(RColorBrewer)

#' 顶刊风格无监督聚类分析函数 (V2)
#'
#' @param expr_mat 表达矩阵 (行=基因, 列=样本)。建议使用 Log2 矩阵。
#' @param methods 聚类方法。默认 "all"。
#' @param best_k 指定最佳 K 值。可以是单个数字 (所有方法用同一个 K)，也可以是列表 (如 list(kmeans=3, consensus=2, nmf=2))。
#' @param n_top_genes 使用前多少个高变基因。默认 2000。
#' @param outdir 输出目录。
#' @param seed 随机种子。
#'
#' @return 包含分组信息的列表
run_clustering_analysis<- function(expr_mat, 
                                       methods = "all", 
                                       k_range = 2:6,
                                       best_k = 2, # 默认值，会被下面的逻辑覆盖
                                       n_top_genes = 2000,
                                       outdir = "Final_Clustering_Result",
                                       seed = 123,
                                       nmf_nrun = 30) {
  
  # --- 0. 参数解析与环境准备 ---
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  if ("all" %in% methods) methods <- c("kmeans", "consensus", "nmf")
  
  # 处理 best_k 参数：如果是数字，转为列表；如果是列表，补全缺失项
  if (!is.list(best_k)) {
    k_list <- list(kmeans = best_k, consensus = best_k, nmf = best_k)
  } else {
    k_list <- best_k
    # 补全未指定的为默认值 2，防止报错
    if(is.null(k_list$kmeans)) k_list$kmeans <- 2
    if(is.null(k_list$consensus)) k_list$consensus <- 2
    if(is.null(k_list$nmf)) k_list$nmf <- 2
  }
  
  message(">>> 开始聚类分析流程 (V2 Pro)...")
  message(paste0(">>> 设定 K 值: Kmeans=", k_list$kmeans, ", Consensus=", k_list$consensus, ", NMF=", k_list$nmf))
  
  # --- [内部绘图函数] 顶刊风格 PCA ---
  draw_pca_pub <- function(mat_t, cluster_vec, title, filename) {
    # PCA 计算
    pca_res <- prcomp(mat_t, scale. = TRUE)
    var_explained <- round(pca_res$sdev^2 / sum(pca_res$sdev^2) * 100, 1)
    
    pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], 
                         Sample = rownames(mat_t), 
                         Cluster = factor(cluster_vec))
    
    # 绘图
    p <- ggplot(pca_df, aes(x=PC1, y=PC2, fill=Cluster)) +
      # 添加置信椭圆 (虚线，浅色填充)
      stat_ellipse(geom = "polygon", alpha = 0.15, aes(color = Cluster), linetype = "dashed", show.legend = FALSE) +
      # 核心散点：使用 shape=21 (实心带边框)，白色边框，增加质感
      geom_point(size = 4, shape = 21, color = "white", stroke = 0.8, alpha = 0.9) +
      # 配色：使用 Nature (NPG) 配色
      scale_fill_npg(name = "Subtype") +
      scale_color_npg(name = "Subtype") +
      # 主题微调
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major = element_line(color = "gray90", size = 0.2), # 极简网格
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title = element_text(face = "bold"),
        legend.position = "right",
        aspect.ratio = 1 # 保持正方形比例
      ) +
      labs(title = title, 
           x = paste0("PC1 (", var_explained[1], "%)"), 
           y = paste0("PC2 (", var_explained[2], "%)"))
    
    ggsave(filename, p, width = 6, height = 5, dpi = 300)
  }
  
  # --- [内部绘图函数] 顶刊风格热图 ---
  draw_heatmap_pub <- function(mat, cluster_vec, title, filename) {
    # 排序样本
    ord <- order(cluster_vec)
    mat_ord <- mat[, ord]
    cluster_ord <- cluster_vec[ord]
    
    # 构建注释条
    ann_col <- data.frame(Subtype = factor(cluster_ord))
    rownames(ann_col) <- colnames(mat_ord)
    
    # 颜色配置
    n_cls <- length(unique(cluster_vec))
    # 使用 NPG 配色，如果类别多于 10 个会自动循环，这里假设 <10
    cols <- pal_npg("nrc")(n_cls) 
    names(cols) <- levels(factor(cluster_vec))
    ann_colors <- list(Subtype = cols)
    
    # 数据裁剪（增强对比度）：Z-score 超过 2.5 的拉平
    mat_scaled <- t(scale(t(mat_ord))) # 按行 Scale
    mat_scaled[mat_scaled > 2.5] <- 2.5
    mat_scaled[mat_scaled < -2.5] <- -2.5
    
    # 绘图
    pdf(filename, width = 8, height = 10)
    pheatmap(mat_scaled,
             show_rownames = FALSE, 
             show_colnames = TRUE,
             cluster_cols = FALSE, # 已排序，不聚类
             cluster_rows = TRUE,  # 基因聚类
             annotation_col = ann_col,
             annotation_colors = ann_colors,
             # 经典红蓝配色：蓝(低)-白(中)-红(高)
             color = colorRampPalette(c("#3C5488FF", "white", "#E64B35FF"))(100),
             border_color = NA,    # 去掉格子边框，看起来更干净
             main = title,
             fontsize = 10,
             fontsize_col = 8,
             treeheight_row = 30)
    dev.off()
  }
  
  # --- 1. 数据筛选 ---
  expr_mat <- as.matrix(expr_mat)
  n_genes <- min(n_top_genes, nrow(expr_mat))
  message(paste0(">>> 筛选前 ", n_genes, " 个高变基因..."))
  gene_vars <- apply(expr_mat, 1, var)
  top_genes <- head(order(gene_vars, decreasing = TRUE), n_genes)
  
  mat_gene_row <- expr_mat[top_genes, ]
  mat_sample_row <- t(mat_gene_row)
  
  results_list <- list()
  
  # ============================================================
  # 1. K-means
  # ============================================================
  if ("kmeans" %in% methods) {
    message("\n--- Running K-means ---")
    sub_dir <- file.path(outdir, "Kmeans")
    if (!dir.exists(sub_dir)) dir.create(sub_dir)
    
    # 评估图 (不保存 K，只看趋势)
    pdf(file.path(sub_dir, "00_Evaluation_Elbow_Silhouette.pdf"), width = 10, height = 5)
    p1 <- fviz_nbclust(mat_sample_row, kmeans, method = "wss", k.max = max(k_range)) + labs(subtitle = "Elbow Method")
    p2 <- fviz_nbclust(mat_sample_row, kmeans, method = "silhouette", k.max = max(k_range)) + labs(subtitle = "Silhouette Method")
    print(p1); print(p2)
    dev.off()
    
    # 最终出图 (使用 k_list$kmeans)
    curr_k <- k_list$kmeans
    set.seed(seed)
    km_res <- kmeans(mat_sample_row, centers = curr_k, nstart = 25)
    clusters <- paste0("C", km_res$cluster)
    
    write.csv(data.frame(Sample=names(km_res$cluster), Cluster=clusters), 
              file.path(sub_dir, paste0("Groups_K", curr_k, ".csv")), row.names=F)
    
    draw_pca_pub(mat_sample_row, clusters, paste0("K-means (K=", curr_k, ")"), 
                 file.path(sub_dir, paste0("PCA_K", curr_k, ".pdf")))
    draw_heatmap_pub(mat_gene_row, clusters, paste0("K-means Heatmap (K=", curr_k, ")"), 
                     file.path(sub_dir, paste0("Heatmap_K", curr_k, ".pdf")))
    
    results_list$kmeans <- clusters
  }
  
  # ============================================================
  # 2. Consensus Clustering
  # ============================================================
  if ("consensus" %in% methods) {
    message("\n--- Running Consensus ---")
    sub_dir <- file.path(outdir, "Consensus")
    d_cons <- mat_gene_row - rowMeans(mat_gene_row)
    
    # 评估过程 (生成所有 K 的图)
    res_cons <- ConsensusClusterPlus(d_cons, maxK = max(k_range), reps = 1000, 
                                     pItem = 0.8, pFeature = 1, title = sub_dir, 
                                     clusterAlg = "hc", distance = "pearson", seed = seed, plot = "png")
    
    # 最终出图 (使用 k_list$consensus)
    curr_k <- k_list$consensus
    cls <- res_cons[[curr_k]][["consensusClass"]]
    clusters <- paste0("C", cls)
    names(clusters) <- colnames(mat_gene_row)
    
    write.csv(data.frame(Sample=names(clusters), Cluster=clusters), 
              file.path(sub_dir, paste0("Groups_K", curr_k, ".csv")), row.names=F)
    
    draw_pca_pub(mat_sample_row, clusters, paste0("Consensus (K=", curr_k, ")"), 
                 file.path(sub_dir, paste0("PCA_K", curr_k, ".pdf")))
    draw_heatmap_pub(mat_gene_row, clusters, paste0("Consensus Heatmap (K=", curr_k, ")"), 
                     file.path(sub_dir, paste0("Heatmap_K", curr_k, ".pdf")))
    
    results_list$consensus <- clusters
  }
  
  # ============================================================
  # 3. NMF
  # ============================================================
  if ("nmf" %in% methods) {
    message("\n--- Running NMF ---")
    if (min(mat_gene_row) < 0) {
      warning("⚠️ NMF 跳过：输入包含负数。")
    } else {
      sub_dir <- file.path(outdir, "NMF")
      if (!dir.exists(sub_dir)) dir.create(sub_dir)
      
      # 评估
      estim_res <- nmf(mat_gene_row, rank = k_range, method = "brunet", nrun = 10, seed = seed)
      pdf(file.path(sub_dir, "00_Evaluation_Rank.pdf")); plot(estim_res); dev.off()
      
      # 最终出图 (使用 k_list$nmf)
      curr_k <- k_list$nmf
      message(paste0("   Running NMF Rank=", curr_k, "..."))
      nmf_res <- nmf(mat_gene_row, rank = curr_k, method = "brunet", nrun = nmf_nrun, seed = seed)
      
      cls <- predict(nmf_res)
      clusters <- paste0("C", cls)
      names(clusters) <- names(cls)
      
      write.csv(data.frame(Sample=names(clusters), Cluster=clusters), 
                file.path(sub_dir, paste0("Groups_K", curr_k, ".csv")), row.names=F)
      
      draw_pca_pub(mat_sample_row, clusters, paste0("NMF (K=", curr_k, ")"), 
                   file.path(sub_dir, paste0("PCA_K", curr_k, ".pdf")))
      draw_heatmap_pub(mat_gene_row, clusters, paste0("NMF Heatmap (K=", curr_k, ")"), 
                       file.path(sub_dir, paste0("Heatmap_K", curr_k, ".pdf")))
      
      # NMF 特有的一致性图 (也非常漂亮，建议保留)
      pdf(file.path(sub_dir, paste0("NMF_ConsensusMap_K", curr_k, ".pdf")))
      consensusmap(nmf_res, labRow=NA, labCol=NA, main=paste0("NMF Consensus (K=", curr_k, ")"))
      dev.off()
      
      results_list$nmf <- clusters
    }
  }
  
  message("\n>>> 分析完成！")
  return(results_list)
}

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr); library(purrr)
  library(clusterProfiler); library(ReactomePA); library(msigdbr)
  library(ggplot2); library(grid)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr); library(purrr)
  library(clusterProfiler); library(ReactomePA); library(msigdbr)
  library(ggplot2); library(grid)
})
run_multi_enrich_lollipop <- function(
    de_csv     = "I:/Kras_figure/01.09.2025/RNAseq_Intermediate_vs_Sensitive/DE_Medium_vs_Sensitive_limma_trend.csv",
    outdir     = "I:/Kras_figure/01.09.2025/RNAseq_Intermediate_vs_Sensitive/",
    species    = c("Mus musculus","Homo sapiens"),
    p_for_de   = c("pvalue","padj"),
    p_cut_gene = 0.05,
    lfc_gene   = 1,
    p_for_plot = c("pvalue","padj"),
    term_p_use = c("pvalue","padj"),
    term_p_cut = NULL,       
    db_choices = NULL,
    method     = c("segment","col","bar"), 
    top_up     = 10, top_down = 10,
    up_left    = FALSE,
    
    # ===== 颜色设置 =====
    up_color   = "#f29325",
    down_color = "#007172",
    
    font_size_axis  = 20,   # 坐标轴刻度数字 (theme size)
    font_size_title = 22,   # 标题和轴标签 (theme size)
    font_size_label = 6,    # 通路名称和GeneRatio数字 (geom_text size)
    font_size_group = 10,   # Up/Down 大标签 (annotate size)
    dashed_p   = 0.05,
    dashed_size= 0.6,
    width_mm   = 210, height_mm = 180, dpi = 600,
    title_prefix = NULL
){
  species    <- match.arg(species)
  p_for_de   <- match.arg(p_for_de)
  p_for_plot <- match.arg(p_for_plot)
  term_p_use <- match.arg(term_p_use)
  method     <- match.arg(method)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  title_prefix <- title_prefix %||% tools::file_path_sans_ext(basename(de_csv))
  
  mytheme <- theme(
    # 坐标轴刻度数字大小
    axis.text.x = element_text(hjust = 0.5, size = font_size_axis), 
    axis.text.y = element_text(size = font_size_axis),
    
    # 坐标轴标题大小
    axis.title.x = element_text(size = font_size_axis), 
    axis.title.y = element_text(size = font_size_axis),
    
    # 图例和主标题大小
    plot.title = element_text(hjust = 0.5, size = font_size_title),
    legend.title = element_text(size = font_size_title),
    legend.text  = element_text(size = font_size_axis),
    
    axis.ticks.y = element_blank(),
    axis.line = element_line(size = 1),
    plot.margin = unit(c(1,1,1,1), "cm"),
    legend.position = "right",
    legend.background = element_rect(fill = 'transparent'),
    axis.line.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank()
  )
  
  if (!file.exists(de_csv)) stop("找不到差异文件：", de_csv)
  message(">>> 读取差异文件：", de_csv)
  de <- suppressMessages(readr::read_csv(de_csv, guess_max = 1e6))
  
  pick_first <- function(cands, nm=names(de)) {
    hit <- cands[cands %in% nm]
    if (length(hit)) hit[1] else NA
  }
  gene_col <- pick_first(c("Gene","SYMBOL","Symbol","gene","Gene.symbol","X.ID","GENE","GeneName","Gene.name"))
  lfc_col  <- pick_first(c("log2FoldChange","logFC","log2FC","LFC","log2fc", "avg_log2FC"))
  pval_col <- pick_first(c("pvalue","P.Value","pval","PValue","p_value","p_val"))
  padj_col <- pick_first(c("padj","FDR","adj.P.Val","qvalue","q_value","p.adjust","padj.BH","p_val_adj"))
  
  if (is.na(gene_col) || is.na(lfc_col)) stop("差异文件未识别到 Gene 或 log2FC 列。")
  if (is.na(padj_col) && !is.na(pval_col)) {
    de$padj <- p.adjust(de[[pval_col]], "BH"); padj_col <- "padj"
  }
  
  .pretty_pathway <- function(x){
    y <- gsub("^HALLMARK_", "", x); y <- gsub("_", " ", y)
    repl <- c("TNFA"="TNF-α","NFKB"="NF-κB","TGF BETA"="TGF-β","IL6"="IL-6","IL2"="IL-2",
              "G2M"="G2/M","UV RESPONSE UP"="UV response (up)","UV RESPONSE DN"="UV response (down)")
    for (k in names(repl)) y <- gsub(k, repl[[k]], y, fixed = TRUE)
    y <- tools::toTitleCase(tolower(y))
    caps <- c("DNA","RNA","MYC","KRAS","E2F","TNF-α","NF-κB","G2/M","TGF-β","IL-2","IL-6")
    for (c in caps) y <- gsub(tools::toTitleCase(tolower(c)), c, y, fixed = TRUE)
    y
  }
  
  de_clean <- de %>%
    transmute(
      Gene = .data[[gene_col]],
      LFC  = .data[[lfc_col]],
      pval = if (!is.na(pval_col)) .data[[pval_col]] else NA_real_,
      padj = if (!is.na(padj_col)) .data[[padj_col]] else NA_real_
    ) %>%
    filter(!is.na(Gene), is.finite(LFC))
  
  # ---------- DE 基因 ----------
  if (p_for_de == "pvalue") {
    if (!any(is.finite(de_clean$pval))) stop("文件没有可用的 pvalue。")
    sig <- de_clean %>% filter(is.finite(pval), pval <= p_cut_gene)
  } else {
    if (!any(is.finite(de_clean$padj))) stop("文件没有可用的 padj。")
    sig <- de_clean %>% filter(is.finite(padj), padj <= p_cut_gene)
  }
  up_genes   <- unique(sig$Gene[sig$LFC >=  lfc_gene])
  down_genes <- unique(sig$Gene[sig$LFC <= -lfc_gene])
  if (!length(up_genes) && !length(down_genes))
    stop("阈值后无基因：请放宽 p_cut_gene 或降低 lfc_gene。")
  
  # ---------- 注释包 ----------
  orgdb_pkg <- if (species=="Mus musculus") "org.Mm.eg.db" else "org.Hs.eg.db"
  if (!requireNamespace(orgdb_pkg, quietly = TRUE))
    stop("缺少注释包：", orgdb_pkg, "（BiocManager::install('", orgdb_pkg, "')）")
  OrgDb <- get(orgdb_pkg, envir = asNamespace(orgdb_pkg))
  
  detect_id <- function(x){
    x <- x[!is.na(x)]
    if (length(x) == 0) return("SYMBOL")
    if (mean(grepl("^ENSG|^ENSMUSG|^ENS", x, ignore.case = TRUE)) > 0.7) return("ENSEMBL")
    if (mean(grepl("^[0-9]+$", x)) > 0.7) return("ENTREZID")
    "SYMBOL"
  }
  id_detected <- detect_id(c(up_genes, down_genes))
  message(">>> 基因 ID 识别：", id_detected)
  
  map_to <- function(genes, fromType, toType){
    if (!length(genes)) return(tibble(from=character(), to=character()))
    res <- tryCatch(
      clusterProfiler::bitr(genes, fromType=fromType, toType=toType, OrgDb=OrgDb),
      error=function(e) NULL
    )
    if (is.null(res)) return(tibble(from=character(), to=character()))
    colnames(res)[1:2] <- c("from","to")
    dplyr::distinct(res, from, to)
  }
  
  # ---------- 数据库清单 ----------
  all_db <- c("GO:BP","GO:MF","GO:CC","KEGG","Reactome",
              "MSigDB:H","MSigDB:C2","MSigDB:C5:BP","MSigDB:C5:MF",
              "MSigDB:C5:CC","MSigDB:C6","MSigDB:C7")
  db_vec <- (db_choices %||% all_db) %>% intersect(all_db)
  if (!length(db_vec)) stop("db_choices 为空。")
  
  kegg_code    <- if (species=="Mus musculus") "mmu" else "hsa"
  reactome_org <- if (species=="Mus musculus") "mouse" else "human"
  
  do_one_db <- function(db_tag){
    message(">>> 富集：", db_tag, "  [DE:", p_for_de, "≤", p_cut_gene, "; 作图:", p_for_plot, "]")
    
    if (db_tag %in% c("GO:BP","GO:MF","GO:CC")){
      keyType <- id_detected
      if (!keyType %in% c("SYMBOL","ENSEMBL","ENTREZID")) keyType <- "SYMBOL"
      ont <- sub("GO:", "", db_tag)
      e_up <- if (length(up_genes))
        enrichGO(gene=unique(up_genes), OrgDb=OrgDb, keyType=keyType, ont=ont,
                 pvalueCutoff=1, qvalueCutoff=1, readable=TRUE) else NULL
      e_dn <- if (length(down_genes))
        enrichGO(gene=unique(down_genes), OrgDb=OrgDb, keyType=keyType, ont=ont,
                 pvalueCutoff=1, qvalueCutoff=1, readable=TRUE) else NULL
      
    } else if (db_tag=="KEGG"){
      ent_up   <- map_to(up_genes,   id_detected, "ENTREZID")
      ent_down <- map_to(down_genes, id_detected, "ENTREZID")
      e_up <- if (nrow(ent_up))
        enrichKEGG(unique(ent_up$to), organism = kegg_code,
                   pvalueCutoff=1, qvalueCutoff=1) else NULL
      e_dn <- if (nrow(ent_down))
        enrichKEGG(unique(ent_down$to), organism = kegg_code,
                   pvalueCutoff=1, qvalueCutoff=1) else NULL
      if (!is.null(e_up)) e_up <- setReadable(e_up, OrgDb = OrgDb, keyType = "ENTREZID")
      if (!is.null(e_dn)) e_dn <- setReadable(e_dn, OrgDb = OrgDb, keyType = "ENTREZID")
      
    } else if (db_tag=="Reactome"){
      ent_up   <- map_to(up_genes,   id_detected, "ENTREZID")
      ent_down <- map_to(down_genes, id_detected, "ENTREZID")
      e_up <- if (nrow(ent_up))
        ReactomePA::enrichPathway(ent_up$to, organism=reactome_org,
                                  pvalueCutoff=1, qvalueCutoff=1, readable=TRUE) else NULL
      e_dn <- if (nrow(ent_down))
        ReactomePA::enrichPathway(ent_down$to, organism=reactome_org,
                                  pvalueCutoff=1, qvalueCutoff=1, readable=TRUE) else NULL
      
    } else { # MSigDB
      if (db_tag=="MSigDB:H")      msig <- msigdbr(species=species, category="H")
      else if (db_tag=="MSigDB:C2") msig <- msigdbr(species=species, category="C2")
      else if (db_tag=="MSigDB:C6") msig <- msigdbr(species=species, category="C6")
      else if (db_tag=="MSigDB:C7") msig <- msigdbr(species=species, category="C7")
      else if (db_tag=="MSigDB:C5:BP") msig <- msigdbr(species=species, category="C5", subcategory="BP")
      else if (db_tag=="MSigDB:C5:MF") msig <- msigdbr(species=species, category="C5", subcategory="MF")
      else if (db_tag=="MSigDB:C5:CC") msig <- msigdbr(species=species, category="C5", subcategory="CC")
      
      t2g <- msig %>% dplyr::select(gs_name, gene_symbol) %>% distinct()
      sym_up   <- if (id_detected=="SYMBOL") tibble(to=up_genes) else map_to(up_genes,   id_detected, "SYMBOL")
      sym_down <- if (id_detected=="SYMBOL") tibble(to=down_genes) else map_to(down_genes, id_detected, "SYMBOL")
      e_up <- if (nrow(sym_up))
        enricher(gene = unique(sym_up$to), TERM2GENE = t2g,
                 pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1) else NULL
      e_dn <- if (nrow(sym_down))
        enricher(gene = unique(sym_down$to), TERM2GENE = t2g,
                 pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1) else NULL
      
      if (!is.null(e_up))
        e_up@result$Description <- ifelse(
          grepl("^HALLMARK_", e_up@result$Description),
          .pretty_pathway(e_up@result$Description),
          e_up@result$Description
        )
      if (!is.null(e_dn))
        e_dn@result$Description <- ifelse(
          grepl("^HALLMARK_", e_dn@result$Description),
          .pretty_pathway(e_dn@result$Description),
          e_dn@result$Description
        )
    }
    list(up=e_up, down=e_dn)
  }
  
  to_df <- function(ek, p_for_plot){
    empty <- tibble(
      Description = character(),
      p.adjust    = numeric(),
      pvalue      = numeric(),
      Count       = integer(),
      GeneRatio   = character(),
      minusLog10  = numeric()
    )
    if (is.null(ek)) return(empty)
    tb <- as_tibble(ek@result)
    if (!nrow(tb)) return(empty)
    metric_col <- if (p_for_plot=="padj" && "p.adjust" %in% names(tb)) "p.adjust" else "pvalue"
    tb %>%
      mutate(minusLog10 = -log10(pmax(.data[[metric_col]], .Machine$double.xmin))) %>%
      arrange(.data[[metric_col]], desc(Count))
  }
  
  # ---------- 画图三模式 ----------
  plot_one <- function(db_tag, up_tbl, dn_tbl){
    safe <- gsub("[^A-Za-z0-9_]+", "_", db_tag)
    
    .metric_for_terms <- function(tb){
      if (!nrow(tb)) return(NA_character_)
      if (term_p_use=="padj" && "p.adjust" %in% names(tb)) "p.adjust" else "pvalue"
    }
    eff_cut <- term_p_cut %||% dashed_p
    
    if (nrow(up_tbl)) {
      mc <- .metric_for_terms(up_tbl)
      up_tbl <- up_tbl %>% filter(is.finite(.data[[mc]]), .data[[mc]] <= eff_cut)
    }
    if (nrow(dn_tbl)) {
      mc <- .metric_for_terms(dn_tbl)
      dn_tbl <- dn_tbl %>% filter(is.finite(.data[[mc]]), .data[[mc]] <= eff_cut)
    }
    
    up_top <- up_tbl %>% slice_head(n = top_up)
    dn_top <- dn_tbl %>% slice_head(n = top_down)
    
    if (nrow(up_top)==0 && nrow(dn_top)==0) {
      message("  ⚠️ ", db_tag, " 无显著条目可画。")
      return(NULL)
    }
    
    # ========== segment/col：镜像棒棒糖 ==========
    if (method %in% c("segment","col")){
      sign_up   <- if (up_left) -1 else +1
      sign_down <- if (up_left) +1 else -1
      
      ord_up   <- if (nrow(up_top)) up_top %>% arrange(desc(minusLog10)) %>% pull(Description) else character(0)
      ord_down <- if (nrow(dn_top)) dn_top %>% arrange(desc(minusLog10)) %>% pull(Description) else character(0)
      lvls_raw <- c(ord_down, rev(ord_up))
      lvls     <- unique(lvls_raw)
      
      dat <- bind_rows(
        if (nrow(up_top)) up_top %>% mutate(threshold="Up",   signed_log10 = sign_up   * minusLog10),
        if (nrow(dn_top)) dn_top %>% mutate(threshold="Down", signed_log10 = sign_down * minusLog10)
      ) %>%
        mutate(Description = factor(Description, levels = lvls))
      
      ymax <- max(abs(dat$signed_log10), na.rm = TRUE)
      ylim <- c(-ymax - 0.6, ymax + 0.6)
      
      thr_y <- -log10(pmax(dashed_p, .Machine$double.xmin))
      
      p <- ggplot(dat, aes(x = Description, y = signed_log10)) +
        coord_flip(clip = "off") +
        { if (method=="segment")
          geom_segment(aes(xend = Description, y = 0, yend = signed_log10,
                           linetype = threshold, color = threshold), size = 2)
          else
            geom_col(aes(fill = threshold), width = 0.1) } +
        geom_point(aes(size = Count, color = threshold)) +
        scale_size_continuous(range = c(4,10)) +
        scale_color_manual(values = c("Up"=up_color, "Down"=down_color)) +
        { if (method=="segment")
          scale_linetype_manual(values = c("Up"="dashed","Down"="solid"))
          else guides(linetype="none") } +
        { if (method=="col") scale_fill_manual(values = c("Up"=up_color, "Down"=down_color)) else NULL } +
        geom_hline(yintercept = c(-thr_y, thr_y), color = 'grey60',
                   size = dashed_size, lty='dashed') +
        # 标签 (使用 font_size_label)
        geom_text(data = dat %>% filter(threshold=="Up"),
                  aes(y = if (up_left) 0.1 else -0.1, x=Description, label = Description),
                  hjust = if (up_left) 0 else 1, color = 'black', size = font_size_label,
                  inherit.aes = FALSE) +
        geom_text(data = dat %>% filter(threshold=="Down"),
                  aes(y = if (up_left) -0.1 else 0.1, x=Description, label = Description),
                  hjust = if (up_left) 1 else 0, color = 'black', size = font_size_label,
                  inherit.aes = FALSE) +
        scale_y_continuous(limits = ylim, expand = expansion(add = c(0.1,0.1))) +
        labs(
          x = NULL,
          y = if (p_for_plot=="padj") bquote("-"~Log[10]~"(adj. P)")
          else bquote("-"~Log[10]~"(P value)"),
          title = paste0(title_prefix, " - ", db_tag, " enrichment"),
          color = "Group", linetype = "Group", size = "Count"
        ) +
        theme_bw() + mytheme +
        scale_x_discrete(labels = NULL) +
        # Up/Down 标签 (使用 font_size_group)
        annotate(
          "text", x = max(1, round(length(lvls)*0.25)),
          y = if (up_left)  0.85*max(abs(ylim)) else -0.85*max(abs(ylim)),
          label = "Down", size = font_size_group, color = down_color
        ) +
        annotate(
          "text", x = max(1, round(length(lvls)*0.75)),
          y = if (up_left) -0.85*max(abs(ylim)) else  0.85*max(abs(ylim)),
          label = "Up",   size = font_size_group, color = up_color
        )
      
      ggplot2::ggsave(file.path(outdir, paste0("Enrich_", safe, "_", method, "_", p_for_plot, ".pdf")),
                      p, width = width_mm, height = height_mm, units = "mm", dpi = dpi)
      ggplot2::ggsave(file.path(outdir, paste0("Enrich_", safe, "_", method, "_", p_for_plot, ".png")),
                      p, width = width_mm, height = height_mm, units = "mm", dpi = dpi)
      
      write_csv(dat %>% dplyr::select(Description, threshold, Count, signed_log10),
                file.path(outdir, paste0("Enrich_", safe, "_", method, "_PlotData_", p_for_plot, ".csv")))
      return(p)
    }
    
    # ========== bar：对称条形（p5） ==========
    if (method=="bar"){
      if (!nrow(up_top) && !nrow(dn_top)) {
        message("  ⚠️ ", db_tag, " 无条目可画。")
        return(NULL)
      }
      up_dat <- up_top %>% mutate(threshold="Up",   pvalue_signed =  minusLog10)
      dn_dat <- dn_top %>% mutate(threshold="Down", pvalue_signed = -minusLog10)
      
      up_levels      <- up_dat %>% arrange(pvalue_signed) %>% pull(Description)
      dn_levels      <- dn_dat %>% arrange(pvalue_signed) %>% pull(Description)
      levels_all_raw <- c(dn_levels, up_levels)
      levels_all     <- unique(levels_all_raw) 
      
      dat <- bind_rows(up_dat, dn_dat) %>%
        mutate(
          Description = factor(Description, levels = levels_all),
          Count       = as.integer(Count)
        )
      
      ymax   <- max(abs(dat$pvalue_signed), na.rm = TRUE)
      top_lim <- ceiling((ymax + 0.5)/5)*5
      thr_y   <- -log10(pmax(dashed_p, .Machine$double.xmin))
      
      p4 <- ggplot(dat, aes(x = Description, y = pvalue_signed, fill = threshold)) +
        geom_col() +
        geom_hline(yintercept = c(-thr_y, thr_y), color = 'grey60',
                   size = dashed_size, linetype = 'dashed') +
        coord_flip() +
        # GeneRatio 数字 (使用 font_size_label)
        geom_text(aes(label = GeneRatio), size = font_size_label) + 
        scale_fill_manual(values = c('Up'=up_color, 'Down'=down_color)) +
        labs(
          x = NULL,
          y = bquote("-"~Log[10]~"(P value)"),
          title = paste0(title_prefix, " - ", db_tag, " enrichment")
        )
      
      p5 <- p4 +
        # 通路文字 (使用 font_size_label)
        geom_text(
          data = dat %>% filter(threshold=="Up"),
          aes(y = -0.1, x=Description, label = Description),
          hjust = 1, color = 'black', size = font_size_label, inherit.aes = FALSE
        ) +
        geom_text(
          data = dat %>% filter(threshold=="Down"),
          aes(y = 0.1, x=Description, label = Description),
          hjust = 0, color = 'black', size = font_size_label, inherit.aes = FALSE
        ) +
        scale_x_discrete(labels = NULL) +
        scale_y_continuous(
          expand = expansion(add = c(0.1, 0.1)),
          limits = c(-top_lim, top_lim),
          breaks = seq(-top_lim, top_lim, 5),
          labels = c(seq(top_lim, 0, -5), seq(5, top_lim, 5))
        ) +
        # Up/Down 标签 (使用 font_size_group)
        annotate(
          "text", x = max(2, round(0.8*nrow(dat))), y =  0.67*top_lim,
          label = "Up",   size = font_size_group, color = up_color
        ) +
        annotate(
          "text", x = max(1, round(0.2*nrow(dat))), y = -0.67*top_lim,
          label = "Down", size = font_size_group, color = down_color
        ) +
        theme_bw() + mytheme + theme(legend.position = "none")
      
      ggplot2::ggsave(file.path(outdir, paste0("Enrich_", safe, "_bar_", p_for_plot, ".pdf")),
                      p5, width = width_mm, height = height_mm, units = "mm", dpi = dpi)
      ggplot2::ggsave(file.path(outdir, paste0("Enrich_", safe, "_bar_", p_for_plot, ".png")),
                      p5, width = width_mm, height = height_mm, units = "mm", dpi = dpi)
      write_csv(dat %>% dplyr::select(Description, threshold, Count, pvalue_signed, GeneRatio),
                file.path(outdir, paste0("Enrich_", safe, "_bar_PlotData_", p_for_plot, ".csv")))
      return(p5)
    }
    
    NULL
  }
  
  plots_list  <- list()
  tables_up   <- list()
  tables_down <- list()
  
  for (db in db_vec) {
    er <- do_one_db(db)
    up_tbl <- to_df(er$up, p_for_plot)
    dn_tbl <- to_df(er$down, p_for_plot)
    
    tables_up[[db]]   <- up_tbl
    tables_down[[db]] <- dn_tbl
    
    if (nrow(up_tbl))
      write_csv(up_tbl,
                file.path(outdir, paste0(title_prefix, "_", gsub("[^A-Za-z0-9_]+","_",db),
                                         "_UP_", p_for_plot, ".csv")))
    if (nrow(dn_tbl))
      write_csv(dn_tbl,
                file.path(outdir, paste0(title_prefix, "_", gsub("[^A-Za-z0-9_]+","_",db),
                                         "_DOWN_", p_for_plot, ".csv")))
    
    p_obj <- plot_one(db, up_tbl, dn_tbl)
    plots_list[[db]] <- p_obj
  }
  
  invisible(list(
    up_genes    = up_genes,
    down_genes  = down_genes,
    id_detected = id_detected,
    databases   = db_vec,
    enrich_up   = tables_up, 
    enrich_down = tables_down,  
    plots       = plots_list    
  ))
}


############make_gene_upset_plot################
make_gene_upset_plot <- function(
    de_list,
    outdir,
    prefix              = "DE",
    gene_cols           = c("Symbol","SYMBOL","Gene","gene","X.ID"),
    lfc_col             = "log2FoldChange",
    padj_col            = "padj",
    pval_col            = "pvalue",
    padj_cut            = 0.05,
    lfc_cut             = 1,
    direction           = c("both","up","down","separate"),
    width_mm            = 200,
    height_mm           = 120,
    dpi                 = 600,
    csv_basename        = NULL,
    min_intersection    = NULL
){
  # --- 1. 参数校验和依赖检查 ---
  stopifnot(is.list(de_list), length(de_list) > 0, !is.null(names(de_list)))
  stopifnot(is.numeric(lfc_cut) && lfc_cut >= 0)
  direction <- match.arg(direction)
  
  suppressPackageStartupMessages({
    has_cu <- requireNamespace("ComplexUpset", quietly = TRUE)
    has_ur <- requireNamespace("UpSetR", quietly = TRUE)
    if (!has_cu && !has_ur) {
      stop("请先安装 ComplexUpset 或 UpSetR：\n",
           "install.packages('ComplexUpset') 或 install.packages('UpSetR')")
    }
  })
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # --- 2. 内部辅助函数 ---
  .prep_df <- function(df, df_name) {
    # 查找基因列
    gene_col_found <- gene_cols[gene_cols %in% names(df)]
    if (length(gene_col_found) == 0) {
      stop("在 '", df_name, "' 中找不到任何指定的基因列: ", paste(gene_cols, collapse=", "))
    }
    
    # 确保 padj 列存在
    if (!(padj_col %in% names(df)) && (pval_col %in% names(df))) {
      message("在 '", df_name, "' 中找不到 '", padj_col, "' 列，正在从 '", pval_col, "' 计算...")
      df[[padj_col]] <- p.adjust(df[[pval_col]], method = "BH")
    }
    
    # 确保 lfc 列存在
    if (!(lfc_col %in% names(df))) {
      alt_lfc <- c("log2FC", "logFC", "LFC")
      hit <- alt_lfc[alt_lfc %in% names(df)]
      if (length(hit) == 0) stop("在 '", df_name, "' 中找不到 log2FC 列（尝试了: '", lfc_col, "', ", paste(alt_lfc, collapse=", "), "）")
      lfc_col <- hit[1]
    }
    
    df %>%
      dplyr::rename(Gene = dplyr::all_of(gene_col_found[1])) %>%
      dplyr::filter(!is.na(Gene), !is.na(.data[[lfc_col]]), !is.na(.data[[padj_col]])) %>%
      dplyr::select(Gene, padj = dplyr::all_of(padj_col), lfc = dplyr::all_of(lfc_col))
  }
  
  # --- 3. 数据准备：提取和过滤基因集 ---
  
  # 3.1. 从每个 DE 结果中提取上调/下调基因
  processed_sets <- purrr::map(names(de_list), function(nm) {
    df <- .prep_df(de_list[[nm]], nm)
    
    is_sig <- df$padj <= padj_cut & !is.na(df$padj)
    
    list(
      up   = unique(df$Gene[is_sig & df$lfc >=  lfc_cut]),
      down = unique(df$Gene[is_sig & df$lfc <= -lfc_cut]),
      both = unique(df$Gene[is_sig & abs(df$lfc) >= lfc_cut])
    )
  })
  names(processed_sets) <- names(de_list)
  
  # 3.2. 根据 'direction' 参数构建最终用于绘图的基因集列表
  final_gene_sets <- switch(
    direction,
    both = purrr::map(processed_sets, ~.x$both),
    up = purrr::map(processed_sets, ~.x$up),
    down = purrr::map(processed_sets, ~.x$down),
    separate = {
      up_sets <- purrr::map(processed_sets, ~.x$up)
      names(up_sets) <- stringr::str_c(names(up_sets), "_up")
      down_sets <- purrr::map(processed_sets, ~.x$down)
      names(down_sets) <- stringr::str_c(names(down_sets), "_down")
      c(up_sets, down_sets)
    }
  )
  
  # 3.3. 创建基因-集合逻辑矩阵
  universe <- unique(unlist(final_gene_sets, use.names = FALSE))
  logical_matrix_cols <- purrr::map_dfc(final_gene_sets, ~ universe %in% .x)
  logical_matrix <- dplyr::bind_cols(data.frame(Gene = universe), logical_matrix_cols)
  
  # --- 4. 保存逻辑矩阵 ---
  if (is.null(csv_basename)) {
    csv_basename <- paste0(prefix, "_GeneSet_", direction, "_logical_matrix.csv")
  }
  readr::write_csv(logical_matrix, file.path(outdir, csv_basename))
  
  # --- 5. 绘图 ---
  set_names <- names(final_gene_sets)
  plot_object <- NULL
  
  if (has_cu) {
    message("正在使用 ComplexUpset 生成图像...")
    
    # 核心逻辑：动态检测 ComplexUpset 的 API 版本以保持兼容性。
    # 新版本 (>=1.3.0) 使用 `sets` 和 `base_annotations` 参数。
    # 旧版本使用 `intersect` 和 `annotations` 参数。
    # 这段代码通过检查函数签名来自动适应。
    cu_ns <- asNamespace("ComplexUpset")
    cu_upset <- get("upset", envir = cu_ns)
    fml <- names(formals(cu_upset))
    
    arg_sets_name <- if ("sets" %in% fml) "sets" else if ("intersect" %in% fml) "intersect" else "intersect" # 默认为旧版
    ann_key <- if ("base_annotations" %in% fml) "base_annotations" else "annotations"
    
    base_ann <- list(
      'Intersection size' = ComplexUpset::intersection_size(text = list(size = 3))
    )
    
    # 新版 ComplexUpset 有独立的 set_size 注释函数
    if ("set_size" %in% getNamespaceExports("ComplexUpset")) {
      base_ann[['Set size']] <- ComplexUpset::set_size(text = list(size = 3))
    }
    
    args_list <- list(data = as.data.frame(logical_matrix))
    args_list[[arg_sets_name]] <- set_names
    args_list[[ann_key]] <- base_ann
    if ("min_size" %in% fml && !is.null(min_intersection)) {
      args_list[["min_size"]] <- min_intersection
    }
    
    plot_object <- do.call(cu_upset, args_list)
    
    # 保存图像
    gg_filename_base <- file.path(outdir, paste0(prefix, "_Gene_UpSet_", direction))
    ggplot2::ggsave(paste0(gg_filename_base, ".pdf"), plot_object, 
                    width = width_mm, height = height_mm, units = "mm", dpi = dpi)
    ggplot2::ggsave(paste0(gg_filename_base, ".png"), plot_object, 
                    width = width_mm, height = height_mm, units = "mm", dpi = dpi)
    
  } else if (has_ur) {
    message("ComplexUpset 未安装，正在使用 UpSetR 生成图像...")
    
    m_for_upsetr <- as.data.frame(logical_matrix)
    rownames(m_for_upsetr) <- m_for_upsetr$Gene
    m_for_upsetr$Gene <- NULL
    
    pdf_filename <- file.path(outdir, paste0(prefix, "_Gene_UpSet_", direction, ".pdf"))
    grDevices::pdf(pdf_filename, width = width_mm/25.4, height = height_mm/25.4)
    UpSetR::upset(m_for_upsetr, sets = set_names, nsets = length(set_names), order.by = "freq")
    grDevices::dev.off()
  }
  
  # --- 6. 返回结果 ---
  invisible(
    list(
      plot = plot_object,
      matrix = logical_matrix,
      gene_sets = final_gene_sets
    )
  )
}




#######plot_gsea_bubble_top5#################################################
suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(dplyr); library(readr); library(rlang)
})

plot_gsea_bubble_top5 <- function(
    file,
    use_adjusted  = TRUE,
    top_n         = 5,
    top_up_extra  = 3,
    sig_cutoff    = 0.05,
    out_prefix    = NULL,
    width = 5.5, height = 5.5, dpi = 300,
    show_y0 = FALSE,
    alpha_min = 0.20,
    alpha_max = 0.95
){
  # ---- 读取 ----
  # ---- 读取 ----
  gsea_res <- readr::read_csv(file, show_col_types = FALSE)
  nm <- names(gsea_res)
  
  # ---- term 列兼容：优先匹配这些候选列名 ----
  term_candidates <- c("term","pathway","ID","gs_name","Pathway","Description")
  term_col <- term_candidates[term_candidates %in% nm][1]
  if (is.na(term_col)) stop("找不到 term/pathway/ID/gs_name/Pathway/Description 列: ", file)
  names(gsea_res)[names(gsea_res) == term_col] <- "term"
  
  # ---- NES 必须有 ----
  if (!"NES" %in% names(gsea_res)) stop("找不到 NES 列: ", file)
  
  # ---- 集合大小列兼容 ----
  size_candidates <- c("setSize","size","Size")
  size_col <- size_candidates[size_candidates %in% nm][1]
  if (is.na(size_col)) stop("找不到 setSize/size/Size 列: ", file)
  names(gsea_res)[names(gsea_res) == size_col] <- "setSize"
  
  # ---- p 值列兼容：clusterProfiler 常用 p.adjust（也顺手支持 qvalues）----
  if ("p.adjust" %in% names(gsea_res) && !("padj" %in% names(gsea_res))) {
    gsea_res$padj <- gsea_res[["p.adjust"]]
  }
  if (!("padj" %in% names(gsea_res)) && "qvalues" %in% names(gsea_res)) {
    gsea_res$padj <- gsea_res[["qvalues"]]
  }
  
  
  has_padj <- "padj"   %in% names(gsea_res)
  has_pval <- "pval"   %in% names(gsea_res)
  has_pv   <- "pvalue" %in% names(gsea_res)
  
  if (use_adjusted && has_padj) {
    gsea_res <- gsea_res %>% mutate(p_for_sig = padj)
  } else if (has_pval) {
    gsea_res <- gsea_res %>% mutate(p_for_sig = pval)
  } else if (has_pv) {
    gsea_res <- gsea_res %>% mutate(p_for_sig = pvalue)
  } else if (has_padj) {
    gsea_res <- gsea_res %>% mutate(p_for_sig = padj)
  } else {
    stop("既无 padj 也无 pval/pvalue: ", file)
  }
  
  # ---- 避免 -log10(0) ----
  gsea_res <- gsea_res %>% mutate(p_for_sig = pmax(p_for_sig, 1e-300))
  
  # ---- 整理 & 可视化（与你现有代码相同）----
  plot_df <- gsea_res %>%
    arrange(desc(NES)) %>%
    mutate(
      xlab      = row_number(),
      setSize_1 = setSize / 10,
      logP      = -log10(p_for_sig)
    )
  
  logP_min <- min(plot_df$logP, na.rm = TRUE)
  logP_max <- max(plot_df$logP, na.rm = TRUE)
  rng <- max(logP_max - logP_min, 1e-9)
  plot_df <- plot_df %>%
    mutate(alpha_val = alpha_min + (logP - logP_min) / rng * (alpha_max - alpha_min))
  
  label_global <- plot_df %>%
    arrange(p_for_sig, desc(abs(NES))) %>%
    slice_head(n = min(top_n, nrow(.)))
  
  label_up_extra <- plot_df %>%
    filter(NES > 0, p_for_sig < sig_cutoff, !(term %in% label_global$term)) %>%
    arrange(p_for_sig, desc(NES)) %>%
    slice_head(n = min(top_up_extra, nrow(.)))
  
  label_df <- bind_rows(label_global, label_up_extra) %>% distinct(term, .keep_all = TRUE)
  
  blues_pal <- c("#E8F1FA","#CFE2F3","#9FC4E6","#6AA0CF","#3E7CB1","#1F5A8A","#0C3C6E")
  
  p <- ggplot(plot_df, aes(x = xlab, y = NES)) +
    { if(show_y0) geom_hline(yintercept = 0, color = "#A8A8A8", linewidth = 0.5) } +
    geom_point(
      aes(size = setSize_1, color = logP, alpha = alpha_val),
      shape = 16
    ) +
    scale_alpha(range = c(alpha_min, alpha_max), guide = "none") +
    scale_color_gradientn(colors = blues_pal, limits = c(0, logP_max),
                          name = expression(-log[10]*"(p)")) +
    scale_size_continuous(range = c(3,10), name = "Gene set\nsize (×10)") +
    scale_x_continuous(breaks = pretty(plot_df$xlab, n = 6),
                       labels = pretty(plot_df$xlab, n = 6)) +
    scale_y_continuous(breaks = seq(floor(min(plot_df$NES, na.rm=TRUE)),
                                    ceiling(max(plot_df$NES, na.rm=TRUE)), by = 1)) +
    labs(x = "Hallmark gene sets (sorted by NES)",
         y = "Normalized Enrichment Score (NES)",
         title = if(!is.null(out_prefix)) out_prefix else gsub("_synergy_pathways\\.csv$","",basename(file))) +
    theme_classic(base_size = 14) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.line    = element_line(linewidth = 0.7),
          legend.title = element_text(size = 12, face = "bold"),
          legend.key   = element_blank(),
          plot.title   = element_text(hjust = 0.5, face = "bold"))
  
  p_final <- p +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = term),
      size = 3.2, max.overlaps = 20,
      box.padding = 0.45, point.padding = 0.25,
      segment.color = "grey35", segment.size = 0.35,
      color = "black"
    )
  
  if (is.null(out_prefix)) out_prefix <- gsub("\\.csv$","", basename(file))
  ggsave(paste0(out_prefix, "_GSEA_Hallmark_bubble.pdf"), p_final,
         width = width, height = height, dpi = dpi)
  ggsave(paste0(out_prefix, "_GSEA_Hallmark_bubble.png"), p_final,
         width = width, height = height, dpi = dpi, bg = "white")
  
  message("✓ Saved: ", out_prefix, "_GSEA_Hallmark_bubble.[pdf|png]")
  invisible(p_final)
}

# GSEA 气泡图（按 NES 排序；颜色=NES 三色渐变；大小=-log10(p)） -------------------
plot_gsea_bubble_top5 <- function(
    file,
    use_adjusted  = TRUE,     # TRUE: 用 padj；否则回退 pval/pvalue
    top_n         = 5,        # 全局最显著标注数（先按 p，再按 |NES|）
    top_up_extra  = 3,        # 追加标注：NES>0 且显著（p<sig_cutoff）的前 n 个
    sig_cutoff    = 0.05,     # 显著阈值（padj 或 p）
    nes_cap       = 3,        # NES 折顶范围 [-nes_cap, nes_cap]
    out_prefix    = NULL,
    width = 5.5, height = 5.5, dpi = 300,
    show_y0 = FALSE,          # 显示 NES=0 基线
    alpha_min = 0.25,         # 不显著更透明
    alpha_max = 0.95,         # 显著更实
    pretty_names = TRUE,      # 是否美化通路名
    theme_base   = ggplot2::theme_classic(base_size = 14)
){
  # ---- 读取并统一列名 ----
  gsea_res <- readr::read_csv(file, show_col_types = FALSE)
  nm <- names(gsea_res)
  
  # term / pathway
  term_col <- intersect(nm, c("pathway","term","ID","gs_name","Pathway"))
  if(length(term_col)==0) stop("找不到 pathway/term 列: ", file)
  gsea_res <- gsea_res %>% rename(term = !!term_col[1])
  
  # NES
  if(!"NES" %in% names(gsea_res)) stop("找不到 NES 列: ", file)
  
  # p 列（优先 padj）
  has_padj <- "padj"   %in% names(gsea_res)
  has_pval <- "pval"   %in% names(gsea_res)
  has_pv   <- "pvalue" %in% names(gsea_res)
  
  if(use_adjusted && has_padj){
    gsea_res <- gsea_res %>% mutate(p_for_sig = padj)
  } else if(has_pval){
    gsea_res <- gsea_res %>% mutate(p_for_sig = pval)
  } else if(has_pv){
    gsea_res <- gsea_res %>% mutate(p_for_sig = pvalue)
  } else if(has_padj){
    gsea_res <- gsea_res %>% mutate(p_for_sig = padj)
  } else {
    stop("既无 padj 也无 pval/pvalue: ", file)
  }
  gsea_res <- gsea_res %>% mutate(p_for_sig = pmax(p_for_sig, 1e-300)) # 避免 Inf
  
  # 可选 setSize/size（不强制）
  size_col <- intersect(names(gsea_res), c("setSize","size"))
  if(length(size_col)>0) gsea_res <- gsea_res %>% rename(setSize = !!size_col[1])
  
  # ---- 整理数据：按 NES 排序 + 计算 -log10(p) + NES 折顶 + 标签 ----
  plot_df <- gsea_res %>%
    mutate(
      term_label = if (pretty_names) .pretty_pathway(term) else term,
      NES_cap    = pmax(pmin(NES, nes_cap), -nes_cap),
      logP       = -log10(p_for_sig)
    ) %>%
    arrange(desc(NES)) %>%
    mutate(xlab = row_number(),
           is_sig = p_for_sig < sig_cutoff,
           alpha_val = ifelse(is_sig, alpha_max, alpha_min))
  
  # 标签选择
  label_global <- plot_df %>%
    arrange(p_for_sig, desc(abs(NES))) %>%
    slice_head(n = min(top_n, nrow(.)))
  
  label_up_extra <- plot_df %>%
    filter(NES > 0, is_sig, !(term %in% label_global$term)) %>%
    arrange(p_for_sig, desc(NES)) %>%
    slice_head(n = min(top_up_extra, nrow(.)))
  
  label_df <- bind_rows(label_global, label_up_extra) %>%
    distinct(term, .keep_all = TRUE)
  
  # ---- 颜色：三色渐变（与 run_hallmark_gsea 一致语义）----
  # 使用 gradientn，values 定位到 [-nes_cap, 0, nes_cap] 的 0/0.5/1
  nes_cols   <- c("navy", "white", "firebrick")
  nes_limits <- c(-nes_cap, nes_cap)
  nes_values <- c(0, 0.5, 1)
  
  # ---- 绘图 ----
  p <- ggplot(plot_df, aes(x = xlab, y = NES)) +
    { if (show_y0) geom_hline(yintercept = 0, color = "#A8A8A8", linewidth = 0.5) } +
    geom_point(
      aes(size = logP, color = NES_cap, alpha = alpha_val),
      shape = 16
    ) +
    scale_alpha(range = c(alpha_min, alpha_max), guide = "none") +
    scale_size_continuous(name = expression(-log[10]*"(p)"), range = c(3, 10)) +
    scale_color_gradientn(
      colours = nes_cols,
      values  = nes_values,
      limits  = nes_limits,
      name    = "NES"
    ) +
    scale_x_continuous(breaks = pretty(plot_df$xlab, n = 6),
                       labels = pretty(plot_df$xlab, n = 6)) +
    scale_y_continuous(
      breaks = seq(floor(min(plot_df$NES, na.rm = TRUE)),
                   ceiling(max(plot_df$NES, na.rm = TRUE)), by = 1)
    ) +
    labs(
      x = "Hallmark gene sets (sorted by NES)",
      y = "Normalized Enrichment Score (NES)",
      title = if(!is.null(out_prefix)) out_prefix else gsub("_synergy_pathways\\.csv$","",basename(file))
    ) +
    theme_base +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.line    = element_line(linewidth = 0.7),
      legend.title = element_text(size = 12, face = "bold"),
      legend.key   = element_blank(),
      plot.title   = element_text(hjust = 0.5, face = "bold")
    )
  
  # 强化不遮挡的标签（repel）
  set.seed(1)
  p_final <- p +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = term_label),
      size = 3.2,
      max.overlaps = Inf,              # 尽量都放出来
      box.padding = 0.6,
      point.padding = 0.35,
      min.segment.length = 0,
      segment.size = 0.35,
      segment.color = "grey35",
      color = "black",
      force = 1.5,                      # 适度外推
      force_pull = 0.8
    )
  
  # ---- 保存 ----
  if(is.null(out_prefix)) out_prefix <- gsub("\\.csv$","", basename(file))
  ggsave(paste0(out_prefix, "_GSEA_Hallmark_bubble.pdf"), p_final,
         width = width, height = height, dpi = dpi)
  ggsave(paste0(out_prefix, "_GSEA_Hallmark_bubble.png"), p_final,
         width = width, height = height, dpi = dpi, bg = "white")
  
  message("✓ Saved: ", out_prefix, "_GSEA_Hallmark_bubble.[pdf|png]")
  invisible(p_final)
}

library(tidyverse)
library(stringr)
library(tidyverse)
library(stringr)

make_pub_bubble_multi <- function(df,
                                  cell_lines = c("SUIT007", "PDC145"),
                                  treat_map = list(
                                    SUIT007 = c("Src", "MRTX", "MRTX+SRC"),
                                    PDC145  = c("Src", "RMC",  "RMC+SRC")
                                  ),
                                  title_text = "Hallmark GSEA Enrichment",
                                  top_n = 25,
                                  nes_cap = 3,
                                  
                                  # --- 新增功能开关 ---
                                  remove_insig = TRUE,  # TRUE = 删除在所有组中都不显著的通路
                                  sort_by = "abs_nes",  # 排序方式: "abs_nes", "nes", "padj"
                                  
                                  gap_after_first_block = TRUE,
                                  low_col = "#2166ac",  
                                  mid_col = "white",
                                  high_col = "#b2182b", 
                                  bg_cols = c(SUIT007 = "#F0F4FA", PDC145 = "#FFF5F0"), 
                                  bg_alpha = 1,
                                  out_dir = NULL,
                                  out_prefix = NULL,
                                  dpi = 600) {
  
  cap01 <- function(x, lo = 1e-300) pmax(x, lo)
  
  # --- 1. 数据清洗 ---
  d <- df %>% filter(cell_line %in% cell_lines) %>%
    mutate(
      treat = as.character(treat),
      treat_key = paste0(cell_line, "__", treat),
      pathway_pretty = as.character(pathway_pretty),
      NES = as.numeric(NES),
      padj = as.numeric(padj)
    ) %>%
    filter(!is.na(NES), !is.na(padj), !is.na(treat_key), !is.na(pathway_pretty))
  
  x_levels <- unlist(lapply(cell_lines, function(cl) paste0(cl, "__", treat_map[[cl]])), use.names = FALSE)
  gap_level <- "___GAP___"
  
  if (gap_after_first_block && length(cell_lines) >= 2) {
    n1 <- length(treat_map[[cell_lines[1]]])
    x_levels <- c(x_levels[1:n1], gap_level, x_levels[(n1 + 1):length(x_levels)])
  }
  
  d <- d %>%
    mutate(
      treat_key = factor(treat_key, levels = x_levels),
      mlog10 = pmin(-log10(cap01(padj)), 20), 
      sig = padj <= 0.05,
      NES_plotted = pmax(pmin(NES, nes_cap), -nes_cap)
    ) %>%
    filter(!is.na(treat_key))
  
  # --- 2. 关键过滤步骤：删除完全不显著的通路 ---
  if (remove_insig) {
    # 找出至少在一个组中 padj <= 0.05 的通路
    sig_pathways <- d %>%
      group_by(pathway_pretty) %>%
      summarise(min_padj = min(padj, na.rm = TRUE), .groups = "drop") %>%
      filter(min_padj <= 0.05) %>%
      pull(pathway_pretty)
    
    # 过滤数据
    d <- d %>% filter(pathway_pretty %in% sig_pathways)
    
    if(nrow(d) == 0) stop("No significant pathways found across any groups! Try setting remove_insig = FALSE.")
  }
  
  # --- 3. 筛选 Top N 通路 (在过滤后的数据中选 Top) ---
  top_terms <- d %>%
    group_by(pathway_pretty) %>%
    summarise(score = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(score)) %>%
    slice_head(n = top_n) %>%
    pull(pathway_pretty)
  
  d <- d %>% filter(pathway_pretty %in% top_terms)
  
  # --- 4. Y轴排序 ---
  y_stats <- d %>%
    group_by(pathway_pretty) %>%
    summarise(
      mean_nes = mean(NES, na.rm = TRUE),
      max_abs_nes = max(abs(NES), na.rm = TRUE),
      max_sig = max(mlog10, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (sort_by == "nes") {
    y_levels <- y_stats %>% arrange(mean_nes) %>% pull(pathway_pretty)
  } else if (sort_by == "padj") {
    y_levels <- y_stats %>% arrange(max_sig) %>% pull(pathway_pretty)
  } else {
    y_levels <- y_stats %>% arrange(max_abs_nes) %>% pull(pathway_pretty)
  }
  
  d <- d %>% mutate(pathway_pretty = factor(pathway_pretty, levels = y_levels))
  
  # --- 5. 绘图准备 ---
  n_pathways <- length(y_levels)
  x_pos <- seq_along(x_levels); names(x_pos) <- x_levels
  
  block_df <- purrr::map_dfr(cell_lines, function(cl) {
    levs <- paste0(cl, "__", treat_map[[cl]])
    idx <- x_pos[levs]
    center <- mean(idx)
    tibble(cell_line = cl, xmin = min(idx) - 0.5, xmax = max(idx) + 0.5, xcenter = center)
  })
  
  # --- 6. 绘图 ---
  p <- ggplot(d, aes(x = treat_key, y = pathway_pretty)) +
    
    geom_rect(
      data = block_df, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = 0.5, ymax = n_pathways + 0.5), 
      fill = unname(bg_cols[block_df$cell_line]), 
      alpha = bg_alpha
    ) +
    
    geom_point(
      aes(size = mlog10, fill = NES_plotted, alpha = sig),
      shape = 21, color = "grey30", stroke = 0.3
    ) +
    
    geom_text(
      data = block_df, inherit.aes = FALSE,
      aes(x = xcenter, y = n_pathways + 1.0, label = cell_line),
      fontface = "bold", size = 5, color = "black", vjust = 0
    ) +
    
    scale_x_discrete(
      labels = function(x) {
        x <- as.character(x)
        ifelse(x == gap_level, "", sub("^[^_]+__", "", x))
      },
      drop = FALSE
    ) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.1), guide = "none") +
    scale_size_continuous(
      range = c(1.5, 8), 
      name = expression(-log[10](italic("P")["adj"])),
      guide = guide_legend(override.aes = list(fill = "grey50", color="grey30"))
    ) +
    scale_fill_gradient2(
      low = low_col, mid = mid_col, high = high_col, midpoint = 0, 
      name = "NES",
      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title.vjust = 1)
    ) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b=20)),
      axis.text.x = element_text(color="black", face = "bold", size = 11, angle = 45, hjust = 1),
      axis.text.y = element_text(color="black", size = 11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
      legend.position = "right",
      plot.margin = margin(t = 30, r = 20, b = 10, l = 10)
    ) +
    coord_cartesian(clip = "off")
  
  if (!is.null(out_dir) && !is.null(out_prefix)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(out_dir, paste0(out_prefix, ".pdf")),
           p, width = 11, height = 9, units = "in", device = cairo_pdf)
    ggsave(file.path(out_dir, paste0(out_prefix, ".png")),
           p, width = 11, height = 9, units = "in", dpi = dpi, bg = "white")
  }
  
  return(p)
}


plot_hallmark_bubble <- function(
    files,
    groups = NULL,
    output_dir = getwd(),
    prefix = "Hallmark_from_synergy",
    # Visual + ordering knobs aligned with run_hallmark_gsea
    nes_cap = 3,
    sig_padj = 0.05,
    cluster_paths = TRUE,
    cluster_groups = FALSE,
    cluster_metric = c("NES","signed_log10p"),
    dist_method = "euclidean",
    hclust_method = "ward.D2",
    path_order_by = c("cluster","absNES","NES","signed_log10p"),
    path_order_stat = c("max_abs","mean","median"),
    path_order_group = NULL,
    path_order_decreasing = TRUE,
    pretty_names = TRUE,
    grey_non_sig = TRUE,
    non_sig_color = "grey80",
    non_sig_alpha = 0.55,
    width_h = NULL, height_h = NULL,
    width_v = NULL, height_v = NULL,
    units = "mm", dpi = 600,
    theme_obj = NULL,
    # New: significance filter mode
    filter_sig_mode = c("none","all_groups","at_least_k"),
    min_k = NULL
){
  suppressPackageStartupMessages({
    library(dplyr); library(readr); library(tidyr); library(stringr)
    library(ggplot2); library(forcats)
  })
  
  # Pretty pathway names (same as in run_hallmark_gsea)
  .pretty_pathway <- function(x){
    y <- gsub("^HALLMARK_", "", x)
    y <- gsub("_", " ", y)
    repl <- list(
      "TNFA"="TNF-α", "NFKB"="NF-κB", "TGF BETA"="TGF-β",
      "IL6"="IL-6", "IL2"="IL-2", "IL 2"="IL-2", "IL 6"="IL-6",
      "G2M"="G2/M", "MYC"="MYC", "KRAS"="KRAS", "E2F"="E2F",
      "UV RESPONSE UP"="UV response (up)",
      "UV RESPONSE DN"="UV response (down)"
    )
    for(k in names(repl)) y <- gsub(k, repl[[k]], y, fixed = TRUE)
    y <- tools::toTitleCase(tolower(y))
    caps <- c("MYC","KRAS","E2F","DNA","RNA","TNF-α","NF-κB","G2/M","TGF-β","IL-2","IL-6")
    for (c in caps) y <- gsub(tools::toTitleCase(tolower(c)), c, y, fixed = TRUE)
    y
  }
  
  cluster_metric <- match.arg(cluster_metric)
  path_order_by  <- match.arg(path_order_by)
  path_order_stat<- match.arg(path_order_stat)
  filter_sig_mode<- match.arg(filter_sig_mode)
  
  if (is.null(theme_obj)) {
    theme_obj <- theme_bw(base_size = 12) + theme(
      panel.grid = element_blank(),
      axis.text  = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Reader for one synergy CSV
  read_one <- function(f){
    df <- readr::read_csv(f, show_col_types = FALSE)
    nm <- names(df)
    term_col <- intersect(nm, c("pathway","term","ID","gs_name","Pathway"))
    if(length(term_col)==0) stop("No pathway/term column found: ", f)
    if(!"NES" %in% nm) stop("NES column missing: ", f)
    
    df <- df %>% dplyr::rename(Pathway = !!term_col[1])
    if ("padj" %in% nm) pad <- df$padj
    else if ("pval" %in% nm) pad <- df$pval
    else if ("pvalue" %in% nm) pad <- df$pvalue
    else stop("No padj/pval/pvalue in: ", f)
    
    size_col <- intersect(nm, c("setSize","size"))
    if(length(size_col)>0) df <- df %>% dplyr::rename(size = !!size_col[1])
    
    tibble::tibble(
      Pathway = as.character(df$Pathway),
      NES     = as.numeric(df$NES),
      padj    = as.numeric(pad),
      size    = if ("size" %in% names(df)) as.numeric(df$size) else NA_real_
    )
  }
  
  tbls <- lapply(files, read_one)
  if (is.null(groups)) groups <- tools::file_path_sans_ext(basename(files))
  stopifnot(length(groups) == length(tbls))
  
  gsea_all <- dplyr::bind_rows(lapply(seq_along(tbls), function(i){
    tbls[[i]] %>% mutate(Group = groups[i])
  })) %>%
    mutate(
      padj = pmax(padj, .Machine$double.xmin),
      Log10_p = -log10(padj),
      NES_cap = pmax(pmin(NES, nes_cap), -nes_cap)
    ) %>%
    dplyr::select(Group, Pathway, NES, NES_cap, padj, Log10_p, size)
  
  # Significance intersection filter (applied BEFORE building label map)
  n_groups <- dplyr::n_distinct(gsea_all$Group)
  keep_paths <- NULL
  if (filter_sig_mode == "all_groups") {
    keep_paths <- gsea_all %>%
      mutate(sig = padj <= sig_padj) %>%
      group_by(Pathway) %>%
      summarise(n_sig = sum(sig, na.rm = TRUE), .groups = "drop") %>%
      filter(n_sig == n_groups) %>%
      pull(Pathway)
  } else if (filter_sig_mode == "at_least_k") {
    if (is.null(min_k)) min_k <- max(2L, floor(0.75 * n_groups))
    keep_paths <- gsea_all %>%
      mutate(sig = padj <= sig_padj) %>%
      group_by(Pathway) %>%
      summarise(n_sig = sum(sig, na.rm = TRUE), .groups = "drop") %>%
      filter(n_sig >= min_k) %>%
      pull(Pathway)
  }
  if (!is.null(keep_paths)) {
    gsea_all <- gsea_all %>% filter(Pathway %in% keep_paths)
    if (nrow(gsea_all) == 0) {
      out_csv <- file.path(output_dir, paste0(prefix, "_No_Intersection_found.csv"))
      write.csv(tibble::tibble(info = "No pathways pass the intersection filter."), out_csv, row.names = FALSE)
      stop(sprintf("No pathways are significant (padj ≤ %.3g) in all %d groups. Wrote: %s",
                   sig_padj, n_groups, out_csv))
    }
  }
  
  # Build pretty label map AFTER filtering so only kept pathways appear as levels
  if (pretty_names) {
    path_map <- tibble::tibble(Pathway = unique(gsea_all$Pathway)) %>%
      mutate(Pathway_label = .pretty_pathway(Pathway))
  } else {
    path_map <- tibble::tibble(Pathway = unique(gsea_all$Pathway),
                               Pathway_label = Pathway)
  }
  gsea_all <- gsea_all %>% left_join(path_map, by = "Pathway")
  
  # Ordering (same logic as run_hallmark_gsea)
  mat_metric <- if (cluster_metric == "NES") {
    gsea_all %>% dplyr::select(Pathway, Group, NES) %>%
      tidyr::pivot_wider(names_from = Group, values_from = NES, values_fill = 0)
  } else {
    gsea_all %>% mutate(signed = sign(NES) * Log10_p) %>%
      dplyr::select(Pathway, Group, signed) %>%
      tidyr::pivot_wider(names_from = Group, values_from = signed, values_fill = 0)
  }
  mat <- as.matrix(mat_metric[,-1, drop = FALSE]); rownames(mat) <- mat_metric[[1]]
  
  path_order <- unique(gsea_all$Pathway)
  grp_order  <- unique(gsea_all$Group)
  
  if (path_order_by == "cluster") {
    if (cluster_paths && nrow(mat) >= 2) {
      d  <- dist(scale(mat), method = dist_method)
      hc <- hclust(d, method = hclust_method)
      path_order <- rownames(mat)[hc$order]
    }
  } else {
    metric_df <- gsea_all %>%
      mutate(metric = dplyr::case_when(
        path_order_by == "absNES"        ~ abs(NES),
        path_order_by == "NES"           ~ NES,
        path_order_by == "signed_log10p" ~ sign(NES) * Log10_p
      )) %>% dplyr::select(Pathway, Group, metric)
    
    if (!is.null(path_order_group)) {
      metric_x <- metric_df %>% filter(Group == path_order_group) %>%
        group_by(Pathway) %>% summarise(val = dplyr::first(metric), .groups = "drop")
    } else {
      metric_x <- metric_df %>%
        group_by(Pathway) %>%
        summarise(
          val = dplyr::case_when(
            path_order_stat == "max_abs" ~ max(abs(metric), na.rm = TRUE),
            path_order_stat == "mean"    ~ mean(metric, na.rm = TRUE),
            TRUE                         ~ stats::median(metric, na.rm = TRUE)
          ),
          .groups = "drop"
        )
    }
    metric_x <- if (isTRUE(path_order_decreasing)) arrange(metric_x, desc(val)) else arrange(metric_x, val)
    path_order <- metric_x$Pathway
  }
  
  if (cluster_groups && ncol(mat) >= 2) {
    d2  <- dist(scale(t(mat)), method = dist_method)
    hc2 <- hclust(d2, method = hclust_method)
    grp_order <- colnames(mat)[hc2$order]
  }
  
  # Factor levels consistent with filtered set
  label_order <- path_map %>%
    mutate(idx = match(Pathway, path_order)) %>%
    arrange(idx) %>% pull(Pathway_label)
  
  plot_df <- gsea_all %>%
    mutate(
      Pathway_label = factor(Pathway_label, levels = label_order),
      Group         = factor(Group, levels = grp_order),
      is_sig        = padj <= sig_padj
    ) %>%
    droplevels()
  
  # Horizontal plot (size = -log10(padj), color = NES)
  if (grey_non_sig) {
    p_h <- ggplot(plot_df, aes(x = Pathway_label, y = Group)) +
      geom_point(
        data = subset(plot_df, !is_sig),
        aes(size = Log10_p),
        color = non_sig_color, alpha = non_sig_alpha, na.rm = TRUE
      ) +
      geom_point(
        data = subset(plot_df,  is_sig),
        aes(size = Log10_p, color = NES_cap),
        alpha = 0.9, na.rm = TRUE
      ) +
      scale_size_continuous(name = "-log10(padj)", range = c(1.8, 8)) +
      scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0, name = "NES") +
      theme_obj +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9)) +
      labs(title = paste0(prefix, " - Hallmark GSEA (Horizontal)"), x = NULL, y = NULL) +
      scale_x_discrete(drop = TRUE)
  } else {
    p_h <- ggplot(plot_df, aes(x = Pathway_label, y = Group)) +
      geom_point(aes(size = Log10_p, color = NES_cap), alpha = 0.9, na.rm = TRUE) +
      scale_size_continuous(name = "-log10(padj)", range = c(1.8, 8)) +
      scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0, name = "NES") +
      theme_obj +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9)) +
      labs(title = paste0(prefix, " - Hallmark GSEA (Horizontal)"), x = NULL, y = NULL) +
      scale_x_discrete(drop = TRUE)
  }
  
  # Vertical plot
  if (grey_non_sig) {
    p_v <- ggplot(plot_df, aes(x = Group, y = Pathway_label)) +
      geom_point(
        data = subset(plot_df, !is_sig),
        aes(size = Log10_p),
        color = non_sig_color, alpha = non_sig_alpha, na.rm = TRUE
      ) +
      geom_point(
        data = subset(plot_df,  is_sig),
        aes(size = Log10_p, color = NES_cap),
        alpha = 0.9, na.rm = TRUE
      ) +
      scale_size_continuous(name = "-log10(padj)", range = c(1.8, 8)) +
      scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0, name = "NES") +
      theme_obj +
      labs(title = paste0(prefix, " - Hallmark GSEA (Vertical)"), x = NULL, y = NULL) +
      scale_y_discrete(limits = rev(levels(plot_df$Pathway_label)), drop = TRUE)
  } else {
    p_v <- ggplot(plot_df, aes(x = Group, y = Pathway_label)) +
      geom_point(aes(size = Log10_p, color = NES_cap), alpha = 0.9, na.rm = TRUE) +
      scale_size_continuous(name = "-log10(padj)", range = c(1.8, 8)) +
      scale_color_gradient2(low = "navy", mid = "white", high = "firebrick", midpoint = 0, name = "NES") +
      theme_obj +
      labs(title = paste0(prefix, " - Hallmark GSEA (Vertical)"), x = NULL, y = NULL) +
      scale_y_discrete(limits = rev(levels(plot_df$Pathway_label)), drop = TRUE)
  }
  
  # Size heuristics (mm)
  if (is.null(width_h))  width_h  <- max(140, 5 + 0.22 * length(levels(plot_df$Pathway_label)) * 10)
  if (is.null(height_h)) height_h <- max( 60, 5 + 0.50 * length(levels(plot_df$Group))          * 10)
  if (is.null(width_v))  width_v  <- max( 90, 5 + 0.55 * length(levels(plot_df$Group))          * 10)
  if (is.null(height_v)) height_v <- max(120, 5 + 0.20 * length(levels(plot_df$Pathway_label))  * 10)
  
  ggsave(file.path(output_dir, paste0(prefix, "_Bubble_Horizontal.pdf")),
         p_h, width = width_h, height = height_h, units = units, dpi = dpi)
  ggsave(file.path(output_dir, paste0(prefix, "_Bubble_Horizontal.png")),
         p_h, width = width_h, height = height_h, units = units, dpi = dpi, bg = "white")
  ggsave(file.path(output_dir, paste0(prefix, "_Bubble_Vertical.pdf")),
         p_v, width = width_v, height = height_v, units = units, dpi = dpi)
  ggsave(file.path(output_dir, paste0(prefix, "_Bubble_Vertical.png")),
         p_v, width = width_v, height = height_v, units = units, dpi = dpi, bg = "white")
  
  write.csv(plot_df %>% dplyr::select(Group, Pathway, Pathway_label, NES, padj, Log10_p, NES_cap, size),
            file.path(output_dir, paste0(prefix, "_PlotUsed.csv")), row.names = FALSE)
  
  list(horizontal = p_h, vertical = p_v, plot_df = plot_df, outdir = output_dir)
}




# ========== 8) 单基因 AUC ==========
# 需要：DESeq2, pROC；若做批次校正需 limma
suppressPackageStartupMessages({ library(DESeq2); library(pROC) })

# 生成用于AUC的表达矩阵：VST 或 归一化log
prep_expr_for_auc <- function(counts, colData, method = c("vst","norm_log"),
                              batch_col = NULL) {
  method <- match.arg(method)
  stopifnot(identical(colnames(counts), rownames(colData)))
  
  if (method == "vst") {
    dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ 1)
    vs  <- vst(dds, blind = TRUE)
    mat <- assay(vs)
  } else {
    dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ 1)
    dds <- estimateSizeFactors(dds)
    mat <- counts(dds, normalized = TRUE)
    mat <- log2(mat + 1)
  }
  
  # 可选：批次校正（对表达矩阵做线性去批次）
  if (!is.null(batch_col) && batch_col %in% colnames(colData)) {
    if (!requireNamespace("limma", quietly = TRUE)) {
      warning("未安装 limma，跳过批次校正。"); return(mat)
    }
    batch <- factor(colData[[batch_col]])
    mat   <- limma::removeBatchEffect(mat, batch = batch)
  }
  mat
}

# 基于表达矩阵计算单基因AUC
gene_auc_from_matrix <- function(expr_mat, colData, positive_class = "Resistant",
                                 out_csv = NULL) {
  y <- factor(colData$Group,
              levels = c(setdiff(levels(factor(colData$Group)), positive_class), positive_class))
  
  auc_vals <- apply(expr_mat, 1, function(gx) {
    ro <- tryCatch(pROC::roc(response = y, predictor = as.numeric(gx), quiet = TRUE, direction = "auto"),
                   error = function(e) NULL)
    if (is.null(ro)) return(NA_real_)
    as.numeric(pROC::auc(ro))
  })
  
  auc_tbl <- tibble::tibble(Symbol = rownames(expr_mat), AUC = as.numeric(auc_vals)) |>
    dplyr::arrange(dplyr::desc(AUC))
  
  if (!is.null(out_csv)) utils::write.csv(auc_tbl, out_csv, row.names = FALSE)
  auc_tbl
}

# ------------ Drop-in replacement: gene_auc_performance ------------
# ------------ Drop-in replacement: gene_auc_performance ------------
gene_auc_performance <- function(
    counts_filtered,
    colData,
    method = c("vst","log2cpm","log1p"),
    batch_col = NULL,                  # e.g., "Batch"
    class_col = "Group",               # label column in colData
    positive_class = "Resistant",      # which class is treated as positive for 'auc'
    sample_col = NULL,                 # set if rownames(colData) are not sample IDs (e.g., "Sample")
    allow_case_insensitive = TRUE,
    normalize_fun = NULL,              # custom ID normalizer (function)
    output_dir = NULL,
    prefix = "ProjA"
){
  method <- match.arg(method)
  if (is.null(counts_filtered) || nrow(counts_filtered) == 0 || ncol(counts_filtered) == 0)
    stop("counts_filtered is empty.")
  
  # ---------- helpers ----------
  .norm_default <- function(x) toupper(gsub("[^A-Z0-9]", "", trimws(as.character(x))))
  .norm <- if (is.null(normalize_fun)) .norm_default else normalize_fun
  .numeric_rownames <- function(x) length(x) > 0 && all(grepl("^[0-9]+$", x))
  
  # ---------- force rownames(colData) to sample IDs ----------
  if (!is.null(sample_col) && sample_col %in% colnames(colData)) {
    rownames(colData) <- as.character(colData[[sample_col]])
  } else if (is.null(rownames(colData)) || any(is.na(rownames(colData))) || .numeric_rownames(rownames(colData))) {
    cand <- intersect(c("Sample","Cell","sample","cell","SampleID","ID"), colnames(colData))
    if (length(cand) == 0)
      stop("Could not determine sample IDs in colData. Set sample_col (e.g., 'Sample').")
    rownames(colData) <- as.character(colData[[cand[1]]])
  }
  
  # ---------- align colData rows to counts columns ----------
  cn <- colnames(counts_filtered)
  rn <- rownames(colData)
  
  m <- match(cn, rn)
  if (anyNA(m) && allow_case_insensitive) {
    m2 <- match(toupper(cn), toupper(rn)); m[is.na(m)] <- m2[is.na(m)]
  }
  if (anyNA(m)) {
    m3 <- match(.norm(cn), .norm(rn)); m[is.na(m)] <- m3[is.na(m)]
  }
  if (anyNA(m)) {
    drop_idx <- which(is.na(m))
    warning(sprintf(
      "Dropping %d samples not found in colData after normalization: %s",
      length(drop_idx), paste0(cn[drop_idx], collapse = ", ")
    ))
    keep_idx <- setdiff(seq_along(cn), drop_idx)
    if (length(keep_idx) == 0) stop("No overlapping samples between counts and colData.")
    counts_filtered <- counts_filtered[, keep_idx, drop = FALSE]
    cn <- colnames(counts_filtered)
    m  <- match(.norm(cn), .norm(rownames(colData)))
    if (anyNA(m)) stop("Failed to align counts with colData even after dropping unmatched samples.")
  }
  
  colData_aligned <- colData[m, , drop = FALSE]
  rownames(colData_aligned) <- cn
  stopifnot(identical(colnames(counts_filtered), rownames(colData_aligned)))
  
  # ---------- labels: must be binary ----------
  if (!class_col %in% colnames(colData_aligned))
    stop(sprintf("class_col='%s' not found in colData.", class_col))
  y <- factor(colData_aligned[[class_col]])
  if (nlevels(y) != 2L)
    stop(sprintf("gene_auc_performance supports binary classes only; found %d levels in '%s'.",
                 nlevels(y), class_col))
  if (!positive_class %in% levels(y))
    stop(sprintf("positive_class='%s' not present in '%s'.", positive_class, class_col))
  neg_class <- setdiff(levels(y), positive_class)[1]
  y_bin <- as.integer(y == positive_class)
  n_pos <- sum(y_bin == 1); n_neg <- sum(y_bin == 0)
  if (n_pos == 0 || n_neg == 0)
    stop("Both positive and negative classes must be present to compute AUC.")
  
  # ---------- expression transform ----------
  expr <- switch(
    method,
    vst = {
      if (!requireNamespace("DESeq2", quietly = TRUE))
        stop("DESeq2 is required for method='vst'.")
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_filtered,
                                            colData = colData_aligned,
                                            design = ~ 1)
      vs <- DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
      as.matrix(SummarizedExperiment::assay(vs))
    },
    log2cpm = {
      lib <- colSums(counts_filtered)
      lib[lib == 0] <- 1
      cpm <- t(t(counts_filtered) / lib * 1e6)
      log2(cpm + 1)
    },
    log1p = {
      log1p(counts_filtered)
    }
  )
  
  # optional batch removal
  if (!is.null(batch_col)) {
    if (!batch_col %in% colnames(colData_aligned))
      stop(sprintf("batch_col='%s' not found in colData.", batch_col))
    if (requireNamespace("limma", quietly = TRUE)) {
      expr <- limma::removeBatchEffect(expr, batch = factor(colData_aligned[[batch_col]]))
    } else {
      warning("Package 'limma' not available; skip batch removal.")
    }
  }
  
  # ---------- fast per-gene AUC (Mann–Whitney) ----------
  fast_auc <- function(x, y01){
    r <- rank(x, ties.method = "average", na.last = "keep")
    n1 <- sum(y01 == 1, na.rm = TRUE)
    n0 <- sum(y01 == 0, na.rm = TRUE)
    if (n1 == 0L || n0 == 0L) return(NA_real_)
    (sum(r[y01 == 1], na.rm = TRUE) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  auc_pos <- apply(expr, 1, fast_auc, y01 = y_bin)
  names(auc_pos) <- rownames(expr)
  
  # ---------- assemble: bidirectional + corrected ----------
  auc_resistant <- if (positive_class == "Resistant") auc_pos else (1 - auc_pos)
  auc_sensitive <- 1 - auc_resistant  # inverse for the other class
  best_class    <- ifelse(auc_resistant >= 0.5, "Resistant", "Sensitive")
  auc_corrected <- pmax(auc_resistant, auc_sensitive)
  signed_auc    <- (auc_resistant - 0.5) * 2  # [-1,1], positive means Resistant
  
  res <- tibble::tibble(
    gene          = rownames(expr),
    # legacy outputs (kept for compatibility):
    auc           = if (positive_class == "Resistant") auc_resistant else auc_sensitive,
    direction     = ifelse( (if (positive_class == "Resistant") auc_resistant else auc_sensitive) >= 0.5,
                            positive_class, setdiff(c("Resistant","Sensitive"), positive_class)[1]),
    # new richer outputs:
    auc_resistant = as.numeric(auc_resistant),
    auc_sensitive = as.numeric(auc_sensitive),
    auc_corrected = as.numeric(auc_corrected),
    signed_auc    = as.numeric(signed_auc),
    best_class    = best_class
  ) |>
    dplyr::arrange(dplyr::desc(auc_corrected))
  
  # ---------- export ----------
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(res, file.path(output_dir, paste0(prefix, "_gene_auc_corrected.csv")))
    # also export class-specific full tables for convenience
    res_res <- res |> dplyr::arrange(dplyr::desc(auc_resistant))
    res_sen <- res |> dplyr::arrange(dplyr::desc(auc_sensitive))
    readr::write_csv(res_res, file.path(output_dir, paste0(prefix, "_gene_auc_resistant.csv")))
    readr::write_csv(res_sen, file.path(output_dir, paste0(prefix, "_gene_auc_sensitive.csv")))
  }
  
  return(res)
}
# ------------ end: gene_auc_performance ------------


# ============================ One-click RNA-seq pipeline (FIXED & COMPLETE) ============================
run_rnaseq_pipeline <- function(
    cts, cells, fix_target = c("cells","cts"),
    output_dir = "Results", prefix = "Project",
    # --- filtering ---
    min_count = 10, min_samples_prop = 0.10,
    # --- design / contrast ---
    group_levels = c("Sensitive","Resistant"),
    design_formula = ~ Group,
    contrast = c("Group","Resistant","Sensitive"),
    # --- steps ---
    do_gsea = TRUE,
    do_auc  = TRUE,
    do_heatmap = TRUE,
    # --- counts ---
    cts_gene_id_col  = NULL,                           
    cts_prefer       = c("auto","#ID","ENSEMBL","SYMBOL"),
    cts_dedup_genes  = c("sum","mean","max","first","none"),
    # --- AUC options ---
    auc_method = c("vst","norm_log"),
    auc_batch_col = NULL,
    positive_class = NULL,      
    # --- HARMONIZE ---
    harm_cell_col            = NULL,
    harm_group_col           = NULL,
    harm_gene_col            = NULL,        
    harm_drop_missing        = TRUE,
    harm_reorder_cols        = TRUE,
    harm_dedup_genes         = c("sum","max","mean","first","none"),
    harm_combine_dupe_samples= c("sum","mean","stop"),
    harm_alias_table         = NULL,
    harm_allow_oneoff        = 0,
    harm_to_integer          = TRUE,
    harm_verbose             = TRUE,
    # --- GSEA ---
    gsea_species        = "Homo sapiens",
    gsea_collection     = "H",
    gsea_rank_by        = c("stat","signed_logp","log2FC"),
    gsea_sig_in_prop    = NULL,
    gsea_sig_padj       = 0.05,
    gsea_min_groups     = NULL,
    gsea_padj_cut       = NULL,
    gsea_plot_only_sig  = NULL,
    gsea_cluster_paths  = TRUE,
    gsea_cluster_groups = FALSE,
    gsea_path_order_by  = c("cluster","absNES","NES","signed_log10p"),
    gsea_path_order_stat= c("max_abs","mean","median"),
    gsea_path_order_group = NULL,
    gsea_path_order_decreasing = TRUE,
    gsea_grey_non_sig   = FALSE,
    gsea_color_sig_padj = NULL,
    gsea_non_sig_color  = "grey80",
    gsea_non_sig_alpha  = 0.55,
    # --- Heatmap ---
    heatmap_matrix      = c("vst","counts"),
    heatmap_gene_select = c("mad","padj","stat","lfc","none"), # 确保支持这些选项
    heatmap_top_n       = 100,
    heatmap_padj_max    = 0.05,
    heatmap_scale_mode  = c("col_z","col_minmax","none"),
    heatmap_cap_z       = 2,
    heatmap_transform   = c("none","log2","log1p"),
    heatmap_palette     = c("navy-white-firebrick","blue-white-red","viridis","magma","plasma","terrain"),
    heatmap_label_wrap  = 28,
    heatmap_dist_rows   = c("euclidean","spearman","pearson","cosine","manhattan"),
    heatmap_dist_cols   = c("euclidean","spearman","pearson","cosine","manhattan"),
    heatmap_hclust_method = "ward.D2",
    heatmap_row_split_by  = "Group",
    heatmap_ann_colors    = NULL,
    heatmap_save          = TRUE
){
  suppressPackageStartupMessages({
    library(Matrix); library(matrixStats); library(DESeq2); library(SummarizedExperiment)
  })
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---------- Prebuilt map for Gene ID -> SYMBOL ----------
  sym_map <- NULL
  if (is.data.frame(cts) && !is.null(cts_gene_id_col) &&
      cts_gene_id_col %in% names(cts) && "SYMBOL" %in% names(cts)) {
    key <- as.character(cts[[cts_gene_id_col]])
    val <- as.character(cts[["SYMBOL"]])
    keep <- !is.na(key) & key != "" & !is.na(val) & val != ""
    mp   <- unique(data.frame(key = key[keep], sym = val[keep], stringsAsFactors = FALSE))
    mp   <- mp[!duplicated(mp$key), , drop = FALSE]
    sym_map <- setNames(mp$sym, mp$key)
  }
  
  # ---------- 1) counts robustify ----------
  cts_prefer      <- match.arg(cts_prefer)
  cts_dedup_genes <- match.arg(cts_dedup_genes)
  cts <- ensure_integer_counts(
    cts,
    gene_id_col  = cts_gene_id_col,
    prefer       = cts_prefer,
    dedup_genes  = cts_dedup_genes,
    verbose      = TRUE
  )
  
  # ---------- 2) harmonize sample names ----------
  fix_target <- match.arg(fix_target)
  harm_dedup_genes         <- match.arg(harm_dedup_genes)
  harm_combine_dupe_samples<- match.arg(harm_combine_dupe_samples)
  hz <- harmonize_sample_names(
    cts          = cts,
    cells        = cells,
    cell_col     = harm_cell_col,
    group_col    = harm_group_col,
    gene_col     = harm_gene_col, 
    fix_target   = fix_target,
    drop_missing = harm_drop_missing,
    reorder_cols = harm_reorder_cols,
    dedup_genes  = harm_dedup_genes,
    combine_dupe_samples = harm_combine_dupe_samples,
    alias_table  = harm_alias_table,
    allow_oneoff = harm_allow_oneoff,
    to_integer   = harm_to_integer,
    verbose      = harm_verbose
  )
  cts2   <- hz$cts
  cells2 <- hz$cells
  
  # ---------- 3) build colData ----------
  colData <- build_coldata(cells2, colnames(cts2),
                           group_col = "Group",
                           group_levels = group_levels,
                           batch_col = NULL)
  
  # ---------- 4) filter low-count genes ----------
  filt <- filter_counts(cts2, min_count = min_count, min_samples_prop = min_samples_prop)
  counts_filtered <- filt$counts
  utils::write.csv(counts_filtered, file.path(output_dir, paste0(prefix, "_Filtered_Counts.csv")))
  
  # ---------- 5) ALIGN (ROBUST FIX) ----------
  colData <- as.data.frame(colData)
  if (is.null(rownames(colData))) {
    stop("colData has no rownames. build_coldata failed to set sample IDs.")
  }
  colData <- colData[colnames(counts_filtered), , drop = FALSE]
  if (!"Sample" %in% colnames(colData)) {
    colData$Sample <- rownames(colData)
  }
  if (!identical(colnames(counts_filtered), rownames(colData))) {
    stop(paste("Alignment failed.",
               "Counts cols:", paste(head(colnames(counts_filtered)), collapse=","),
               "colData rows:", paste(head(rownames(colData)), collapse=",")))
  }
  
  # ---------- 6) QC ----------
  qc <- qc_plots(counts_filtered, colData, output_dir, prefix)
  
  # ---------- 7) DESeq2 ----------
  de <- run_deseq2(counts_filtered, colData,
                   design_formula = design_formula,
                   contrast = contrast,
                   output_dir = output_dir, prefix = prefix)
  
  # save VST
  vst_mat <- de$vst_mat
  utils::write.csv(vst_mat, file.path(output_dir, paste0(prefix, "_VST_Matrix.csv")))
  
  # Ensure Symbol column
  de$res <- as.data.frame(de$res)
  if (!("Symbol" %in% names(de$res) || "SYMBOL" %in% names(de$res))) {
    rid <- rownames(de$res)
    sym_vec <- NULL
    if (!is.null(sym_map)) sym_vec <- unname(sym_map[match(rid, names(sym_map))])
    if (is.null(sym_vec) || all(is.na(sym_vec))) sym_vec <- rid 
    de$res[["Symbol"]] <- sym_vec
  }
  
  # ---------- 8) Volcano & Summary ----------
  volcano <- plot_volcano(de$res, output_dir = output_dir, prefix = prefix)
  
  de_summary <- summarize_de_results(
    res_df = de$res, 
    output_dir = output_dir, 
    prefix = prefix, 
    padj_cutoff = 0.05, 
    log2fc_cutoff = 1,
    remove_pattern = "^NewGene_"
  )
  
  # ---------- 9) Hallmark GSEA ----------
  gsea_out <- NULL
  if (isTRUE(do_gsea)) {
    gsea_out <- run_hallmark_gsea(
      x = de$res, output_dir = output_dir, prefix = prefix,
      species = gsea_species, collection = gsea_collection,
      rank_by = match.arg(gsea_rank_by),
      sig_padj = gsea_sig_padj,
      grey_non_sig = gsea_grey_non_sig,
      path_order_by = match.arg(gsea_path_order_by),
      cluster_paths = gsea_cluster_paths
    )
  }
  
  # ---------- 10) AUC ----------
  auc_out <- NULL
  if (isTRUE(do_auc)) {
    if (is.null(positive_class)) positive_class <- as.character(group_levels[2L])
    auc_out <- gene_auc_performance(
      counts_filtered = counts_filtered, colData = colData,
      method = match.arg(auc_method), batch_col = auc_batch_col,
      positive_class = positive_class, output_dir = output_dir, prefix = prefix
    )
  }
  
  # ---------- 11) Heatmap (UPDATED LOGIC) ----------
  heatmap_obj <- NULL
  if (isTRUE(do_heatmap)) {
    heatmap_matrix <- match.arg(heatmap_matrix)
    select_mode    <- match.arg(heatmap_gene_select)
    
    # 1. 准备矩阵
    base_mat <- switch(heatmap_matrix, "vst" = vst_mat, "counts" = as.matrix(counts_filtered))
    
    # 2. 准备筛选依据
    sel_genes <- rownames(base_mat)
    
    # --- ★★★ NEW LOGIC STARTS HERE ★★★ ---
    if (select_mode == "mad") {
      # 变异最大
      mads <- matrixStats::rowMads(base_mat, na.rm = TRUE)
      ord  <- order(mads, decreasing = TRUE)
      sel_genes <- rownames(base_mat)[ord]
    } else if (select_mode %in% c("padj", "lfc", "stat")) {
      # 基于 DESeq2 结果筛选
      res_df <- de$res
      
      # 确保矩阵里的基因都在 DESeq2 结果里
      common_genes <- intersect(rownames(base_mat), rownames(res_df))
      res_df <- res_df[common_genes, , drop = FALSE]
      
      if (select_mode == "padj") {
        # 按 P值 从小到大
        ord <- order(res_df$padj, decreasing = FALSE, na.last = TRUE)
      } else if (select_mode == "lfc") {
        # 按 |log2FC| 从大到小
        ord <- order(abs(res_df$log2FoldChange), decreasing = TRUE, na.last = TRUE)
      } else if (select_mode == "stat") {
        # 按 |stat| 从大到小
        ord <- order(abs(res_df$stat), decreasing = TRUE, na.last = TRUE)
      }
      sel_genes <- rownames(res_df)[ord]
    }
    # --- ★★★ NEW LOGIC ENDS HERE ★★★ ---
    
    # 3. 截取 Top N
    n_keep <- min(heatmap_top_n, length(sel_genes))
    sel_genes <- sel_genes[seq_len(n_keep)]
    
    # 4. 构建绘图矩阵
    hm_mat <- t(base_mat[sel_genes, , drop = FALSE])
    ann_df <- data.frame(cell = rownames(hm_mat),
                         Group = as.character(colData[rownames(hm_mat), "Group"]),
                         stringsAsFactors = FALSE)
    
    # 5. 绘图
    heatmap_obj <- flexi_heatmap(
      data = hm_mat, annotation_df = ann_df, annotation_id_col = "cell",
      annotation_keep = "Group", row_split_by = heatmap_row_split_by,
      transform = match.arg(heatmap_transform), scale_mode = match.arg(heatmap_scale_mode),
      cap_z = heatmap_cap_z, cluster_rows = TRUE, cluster_cols = TRUE,
      palette = match.arg(heatmap_palette), ann_colors = heatmap_ann_colors,
      save = isTRUE(heatmap_save), outdir = output_dir,
      prefix = paste0(prefix, "_Heatmap_", heatmap_matrix, "_", select_mode), # 添加筛选方式到文件名
      return_ht = TRUE
    )
  }
  
  cat(sprintf("\n[run_rnaseq_pipeline] Done. Output: %s\n", normalizePath(output_dir)))
  
  invisible(list(
    cts = cts2, cells = cells2, colData = colData,
    counts_filtered = counts_filtered, vst_matrix = vst_mat,
    qc = qc, dds = de$dds, res = de$res, volcano = volcano,
    summary = de_summary, gsea = gsea_out, auc = auc_out, heatmap = heatmap_obj
  ))
}



# ============================================================
# 分层+总模型 一键跑：按 cell_line 分别做 KO vs NC，再做总体模型
# 依赖：ensure_integer_counts, filter_counts, qc_plots, run_deseq2,
#      plot_volcano, run_hallmark_gsea, gene_auc_performance
# ============================================================
run_rnaseq_pipeline_by_strata <- function(
    counts, colData,
    group_col    = "Group",                 # 分组列（例如 KO / NC）
    group_levels = c("NC","KO"),            # 水平顺序，第二个视为“阳性”/对照2
    strata_col   = "cell_line",             # 分层列（例如 1942/2249/2435/4222）
    per_strata   = TRUE,                    # 是否逐层跑
    overall      = TRUE,                    # 是否总体模型
    interaction_test = FALSE,               # 是否做交互 LRT（检测 cell_line:Group）
    output_dir   = "Results",
    prefix       = "Project",
    # 过滤
    min_count = 10, min_samples_prop = 0.10,
    # GSEA / AUC
    do_gsea = TRUE, do_auc = TRUE,
    auc_method = c("vst","norm_log"), auc_batch_col = NULL
){
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  auc_method <- match.arg(auc_method)
  
  # —— 对齐 & 因子水平 —— 
  counts <- ensure_integer_counts(counts)
  colData <- as.data.frame(colData)
  if (is.null(rownames(colData))) rownames(colData) <- colData[[1]]
  stopifnot(setequal(colnames(counts), rownames(colData)))
  colData <- colData[colnames(counts), , drop = FALSE]
  colData[[group_col]]  <- factor(colData[[group_col]], levels = group_levels)
  
  out <- list()
  
  # ============= 总体模型（控制 cell_line） =============
  if (isTRUE(overall)) {
    design_overall <- if (!is.null(strata_col) && strata_col %in% names(colData)) {
      as.formula(paste("~", strata_col, "+", group_col))
    } else {
      as.formula(paste("~", group_col))
    }
    
    filt <- filter_counts(counts, min_count, min_samples_prop)
    cf   <- filt$counts
    cd   <- colData[colnames(cf), , drop = FALSE]
    cd[[group_col]] <- droplevels(cd[[group_col]])
    pre  <- paste0(prefix, "_overall")
    
    qc <- qc_plots(cf, cd, output_dir, pre)
    
    de <- run_deseq2(cf, cd,
                     design_formula = design_overall,
                     contrast       = c(group_col, group_levels[2], group_levels[1]),
                     output_dir     = output_dir, prefix = pre
    )
    
    # >>> 新增 VST 保存
    vst_mat <- DESeq2::vst(de$dds, blind = TRUE)
    vst_mat <- SummarizedExperiment::assay(vst_mat)
    utils::write.csv(vst_mat, file.path(output_dir, paste0(pre, "_VST_Matrix.csv")))
    
    plot_volcano(de$res, output_dir, pre)
    if (do_gsea) run_hallmark_gsea(de$res, output_dir, pre)
    if (do_auc)  gene_auc_performance(cf, cd, method = auc_method,
                                      batch_col = auc_batch_col,
                                      positive_class = group_levels[2],
                                      output_dir = output_dir, prefix = pre)
    
    # 可选：交互作用 LRT
    if (interaction_test && strata_col %in% names(cd)) {
      full    <- as.formula(paste("~", strata_col, "+", group_col, "+", paste0(strata_col, ":", group_col)))
      reduced <- as.formula(paste("~", strata_col, "+", group_col))
      dds_int <- DESeq2::DESeqDataSetFromMatrix(cf, cd, design = full)
      dds_int <- DESeq2::DESeq(dds_int, test = "LRT", reduced = reduced)
      res_int <- as.data.frame(DESeq2::results(dds_int))
      res_int$Symbol <- rownames(res_int)
      write.csv(res_int, file.path(output_dir, paste0(pre, "_Interaction_LRT.csv")), row.names = FALSE)
      out$overall <- list(qc = qc, de = de, vst = vst_mat, interaction = list(dds = dds_int, res = res_int))
    } else {
      out$overall <- list(qc = qc, de = de, vst = vst_mat)
    }
  }
  
  # ============= 分层（每个 cell_line 内 KO vs NC） =============
  if (isTRUE(per_strata)) {
    lv <- unique(colData[[strata_col]])
    lv <- lv[!is.na(lv)]
    out$per_strata <- vector("list", length(lv))
    names(out$per_strata) <- as.character(lv)
    
    for (s in lv) {
      idx <- which(colData[[strata_col]] == s)
      sub_counts <- counts[, idx, drop = FALSE]
      sub_cd     <- colData[idx, , drop = FALSE]
      if (length(unique(sub_cd[[group_col]])) < 2) {
        message("跳过层 ", s, "（只有一个组别）"); next
      }
      
      filt <- filter_counts(sub_counts, min_count, min_samples_prop)
      cf   <- filt$counts
      cd   <- sub_cd[colnames(cf), , drop = FALSE]
      cd[[group_col]] <- droplevels(cd[[group_col]])
      pre  <- paste0(prefix, "_", strata_col, s)
      
      qc <- qc_plots(cf, cd, output_dir, pre)
      
      de <- run_deseq2(cf, cd,
                       design_formula = as.formula(paste("~", group_col)),
                       contrast       = c(group_col, group_levels[2], group_levels[1]),
                       output_dir     = output_dir, prefix = pre
      )
      
      # >>> 新增 VST 保存
      vst_mat <- DESeq2::vst(de$dds, blind = TRUE)
      vst_mat <- SummarizedExperiment::assay(vst_mat)
      utils::write.csv(vst_mat, file.path(output_dir, paste0(pre, "_VST_Matrix.csv")))
      
      plot_volcano(de$res, output_dir, pre)
      if (do_gsea) run_hallmark_gsea(de$res, output_dir, pre)
      if (do_auc)  gene_auc_performance(cf, cd, method = auc_method,
                                        batch_col = auc_batch_col,
                                        positive_class = group_levels[2],
                                        output_dir = output_dir, prefix = pre)
      
      out$per_strata[[as.character(s)]] <- list(qc = qc, de = de, vst = vst_mat)
    }
  }
  
  invisible(out)
}




suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(ggrepel); library(stringr)})

# —— 多组火山图（Nat Commun 风）——
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(ggrepel); library(stringr)})

# ============= Nat Commun 同款多比较组差异火山图（严格版） =============
natcomm_multi_volcano <- function(de_list,
                                  group_colname   = "cell_line",   # 分面列：与你的 all_DE_results 中的 cell_line 一致
                                  gene_colname    = "Gene",
                                  lfc_colname     = "log2FoldChange",
                                  padj_colname    = "padj",
                                  top_n_per_group = 10,             # 每组挑多少个 Top 基因上色+标注（按 |log2FC|）
                                  fc_thr          = 1,              # 竖向阈值线（原文示例 0.05；RNA-seq常用 1）
                                  palette_colors  = c('#c0322f','#26aee7','#96ca39','#51aa44','#f7bcd0',
                                                      '#f1bd58','#68c2b2','#f74141','#b783b6','#3565ab'),
                                  theme_add       = theme_bw(),
                                  file_prefix     = "NatComm_MultiVolcano",
                                  facet_levels    = NULL,           # 自定义面板顺序（可留空）
                                  width_per_panel = 3.6,            # 每个分面宽度（英寸）
                                  height          = 4.2,            # 图高（英寸）
                                  seed_label      = 233) {
  
  stopifnot(is.list(de_list), length(de_list) >= 1)
  
  # 1) 合并为“原文”格式数据框：dt(gene, neuron_type, avg_logFC, p_val_adj)
  dt <- bind_rows(de_list) %>%
    transmute(
      gene        = .data[[gene_colname]],
      neuron_type = .data[[group_colname]],
      avg_logFC   = as.numeric(.data[[lfc_colname]]),
      p_val_adj   = as.numeric(.data[[padj_colname]])
    ) %>%
    filter(is.finite(avg_logFC)) %>%
    mutate(p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj))
  
  # 因子顺序：遵循原文做法（levels=unique）；或使用用户指定的 facet_levels
  if (is.null(facet_levels)) {
    dt$neuron_type <- factor(dt$neuron_type, levels = unique(dt$neuron_type))
  } else {
    dt$neuron_type <- factor(dt$neuron_type, levels = facet_levels)
  }
  
  # 2) 每组挑 Top 基因（按 |avg_logFC|），并确保基因唯一（原文 distinct 再 top_n）
  sig <- dt %>%
    group_by(neuron_type) %>%
    distinct(gene, .keep_all = TRUE) %>%
    slice_max(order_by = abs(avg_logFC), n = top_n_per_group, with_ties = FALSE) %>%
    ungroup()
  
  # 动态对齐标签所需的 y 上限常数 K（原文用 150；这里按数据自适应）
  y_all   <- -log10(pmax(dt$p_val_adj, 1e-300))
  K_align <- max(y_all, na.rm = TRUE) * 1.08  # 稍微往右留点边距
  
  # 调色板长度不够就循环补齐
  if (length(palette_colors) < nlevels(dt$neuron_type)) {
    palette_colors <- rep(palette_colors, length.out = nlevels(dt$neuron_type))
  }
  names(palette_colors) <- levels(dt$neuron_type)
  
  # 3) 原文式绘图流水线
  # 基图：灰点 + 坐标翻转
  p <- ggplot() +
    geom_point(data = dt,
               aes(x = avg_logFC, y = -log10(p_val_adj)),
               size = 0.8, color = "grey70") +
    coord_flip()
  
  # 分面：一行多列（原文 facet_grid . ~ neuron_type）
  p1 <- p + facet_grid(. ~ neuron_type)
  
  # 只给“目标基因（各组Top）”着色
  p2 <- p1 +
    geom_point(data = sig,
               aes(x = avg_logFC, y = -log10(p_val_adj), color = neuron_type),
               size = 1.2)
  
  # 阈值线 + 主题 + 配色（完全按原文结构）
  p3 <- p2 +
    geom_vline(xintercept = c(-fc_thr, fc_thr),
               linewidth = 0.5, color = "grey50", linetype = "dashed") +
    scale_color_manual(values = palette_colors, drop = FALSE) +
    theme_add +
    theme(
      legend.position = "none",
      panel.grid      = element_blank(),
      axis.text       = element_text(size = 10),
      axis.text.x     = element_text(angle = 45, vjust = 0.8),
      strip.text.x    = element_text(size = 10, face = "bold")
    ) +
    labs(x = "log2 Fold Change", y = expression(-log[10]("adj. P")),
         title = "Nat Commun-style multi-group volcano")
  
  # 4) 标签两种：普通版（p4）与对齐式（p5，严格复刻原文 nudge_y=K - y）
  # 普通标签
  p4 <- p3 +
    geom_text_repel(
      data = sig,
      aes(x = avg_logFC, y = -log10(p_val_adj),
          label = gene, color = neuron_type),
      size = 3, fontface = "bold.italic",
      seed = seed_label
    )
  
  # 对齐式标签（重点）：把标签沿 y 轴推到统一的 K_align
  p5 <- p3 +
    geom_text_repel(
      data = sig,
      aes(x = avg_logFC, y = -log10(p_val_adj),
          label = gene, color = neuron_type),
      fontface = "italic",
      seed = seed_label,
      size = 3,
      min.segment.length = 0,
      force = 12, force_pull = 2,
      box.padding = 0.1,
      max.overlaps = Inf,
      segment.linetype = 3,
      segment.alpha = 0.5,
      nudge_y = K_align - (-log10(sig$p_val_adj)),  # —— 原文思想：对齐到“固定的右侧” —— #
      direction = "y",
      hjust = 0
    )
  
  # 自适应宽度保存（每个分面 width_per_panel 英寸）
  npanels <- nlevels(dt$neuron_type)
  w <- max(8, npanels * width_per_panel)
  
  ggsave(paste0(file_prefix, "_basic.pdf"), p3, width = w, height = height, units = "in", device = cairo_pdf)
  ggsave(paste0(file_prefix, "_basic.png"), p3, width = w, height = height, units = "in", dpi = 300)
  
  ggsave(paste0(file_prefix, "_labels_normal.pdf"), p4, width = w, height = height, units = "in", device = cairo_pdf)
  ggsave(paste0(file_prefix, "_labels_normal.png"), p4, width = w, height = height, units = "in", dpi = 300)
  
  ggsave(paste0(file_prefix, "_labels_aligned.pdf"), p5, width = w, height = height, units = "in", device = cairo_pdf)
  ggsave(paste0(file_prefix, "_labels_aligned.png"), p5, width = w, height = height, units = "in", dpi = 300)
  
  list(dt = dt, sig = sig, p_basic = p3, p_label = p4, p_align = p5,
       K_align = K_align, width = w, height = height)
}


suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel); library(ggnewscale)
})

multi_group_fc_on_y <- function(
    de_list,
    group_colname    = "cell_line",
    gene_colname     = "Gene",
    lfc_colname      = "log2FoldChange",
    padj_colname     = "padj",
    lfc_thr          = 1,
    padj_thr         = 0.05,
    top_n_per_dir    = 5,
    label_genes      = NULL,
    highlight_mode   = c("add","only"),
    user_fill        = "#FFD700",
    user_size        = 3.2,
    user_shape       = 21,
    user_stroke      = 0.8,
    dir_colors       = c(Up="#f46d43", Ns="#bdbdbd", Down="#0080BB"),
    group_band_palette = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462'),
    size_range       = c(0.2, 5),
    alpha_range      = c(0.05, 0.6),
    y_limits         = c(-10, 10),
    # <<< 新增这三个参数 >>>
    group_indicator  = c("tile","legend","none"),
    tile_labels      = FALSE,      # 是否在色块上写字（按你的需求默认不写）
    show_group_legend= TRUE,       # 是否给色块加图例（按你的需求默认显示）
    x_text_angle     = 45,
    tile_height_frac = 0.10,
    file_prefix      = "MultiGroup_FC_on_Y",
    width_per_group  = 2.4,
    height           = 6,
    theme_add        = theme_bw()
){
  highlight_mode <- match.arg(highlight_mode)
  group_indicator <- match.arg(group_indicator)
  
  all_deg <- dplyr::bind_rows(de_list) %>%
    dplyr::transmute(
      group     = .data[[group_colname]],
      gene_name = .data[[gene_colname]],
      log2FoldChange = as.numeric(.data[[lfc_colname]]),
      padj      = as.numeric(.data[[padj_colname]])
    ) %>%
    dplyr::mutate(
      padj = ifelse(is.na(padj), 1, padj),
      direction = dplyr::case_when(
        log2FoldChange >  lfc_thr & padj < padj_thr ~ "Up",
        log2FoldChange < -lfc_thr & padj < padj_thr ~ "Down",
        TRUE ~ "Ns"
      )
    ) %>% dplyr::filter(is.finite(log2FoldChange))
  
  grp_fac <- factor(all_deg$group, levels = unique(all_deg$group))
  all_deg$jittered_group <- jitter(as.numeric(grp_fac), amount = 0.4)
  
  topN <- all_deg %>%
    dplyr::filter(direction %in% c("Up","Down")) %>%
    dplyr::group_by(group, direction) %>%
    dplyr::slice_max(order_by = abs(log2FoldChange), n = top_n_per_dir, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    tidyr::drop_na(gene_name) %>%
    dplyr::left_join(all_deg %>% dplyr::select(group, gene_name, jittered_group),
                     by = c("group","gene_name")) %>%
    dplyr::distinct(group, gene_name, .keep_all = TRUE)
  
  user_df <- NULL
  if (!is.null(label_genes)) {
    tolow <- function(x) tolower(trimws(as.character(x)))
    if (is.list(label_genes)) {
      add_all <- if ("ALL" %in% names(label_genes)) tolow(label_genes[["ALL"]]) else character(0)
      user_df <- purrr::imap_dfr(label_genes, function(v, gname){
        if (gname=="ALL") return(NULL)
        gg <- tolow(v)
        all_deg %>% dplyr::filter(group == gname, tolower(gene_name) %in% gg)
      })
      if (length(add_all) > 0) {
        user_df <- dplyr::bind_rows(user_df, all_deg %>% dplyr::filter(tolower(gene_name) %in% add_all))
      }
    } else {
      gg <- tolow(label_genes)
      user_df <- all_deg %>% dplyr::filter(tolower(gene_name) %in% gg)
    }
    if (!is.null(user_df) && nrow(user_df) > 0 && !"jittered_group" %in% names(user_df)) {
      user_df <- dplyr::left_join(user_df,
                                  all_deg %>% dplyr::select(group, gene_name, jittered_group),
                                  by = c("group","gene_name"))
    }
    if (!is.null(user_df) && nrow(user_df) > 0) {
      user_df <- dplyr::distinct(user_df, group, gene_name, .keep_all = TRUE)
    }
  }
  
  if (!is.null(user_df) && highlight_mode=="only") {
    topN <- user_df
  } else if (!is.null(user_df) && highlight_mode=="add") {
    topN <- dplyr::bind_rows(topN, user_df) %>% dplyr::distinct(group, gene_name, .keep_all = TRUE)
  }
  if (!"jittered_group" %in% names(topN) && nrow(topN) > 0) {
    topN <- dplyr::left_join(topN,
                             all_deg %>% dplyr::select(group, gene_name, jittered_group),
                             by = c("group","gene_name"))
  }
  
  if (is.null(y_limits)) {
    ymax <- max(stats::quantile(abs(all_deg$log2FoldChange), 0.995, na.rm = TRUE), lfc_thr*1.5, 2)
    y_limits <- c(-ymax, ymax)
  }
  y_span <- diff(y_limits)
  
  if (length(group_band_palette) < nlevels(grp_fac)) {
    group_band_palette <- rep(group_band_palette, length.out = nlevels(grp_fac))
  }
  names(group_band_palette) <- levels(grp_fac)
  
  p <- ggplot() +
    geom_jitter(
      data = all_deg %>% dplyr::filter(direction == "Ns"),
      aes(x = jittered_group, y = log2FoldChange, color = direction,
          size = abs(log2FoldChange), alpha = abs(log2FoldChange)),
      width = 0, shape = 16
    ) +
    geom_jitter(
      data = all_deg %>% dplyr::filter(direction != "Ns"),
      aes(x = jittered_group, y = log2FoldChange, color = direction,
          size = abs(log2FoldChange), alpha = abs(log2FoldChange)),
      width = 0, shape = 16
    )
  
  if (group_indicator == "tile") {
    band_y  <- y_limits[1] + 0.08 * y_span
    band_h  <- tile_height_frac * y_span
    band_df <- data.frame(x = seq_len(nlevels(grp_fac)), group = levels(grp_fac), y = band_y)
    
    p <- p + geom_tile(
      data = band_df,
      aes(x = x, y = y, fill = group),
      height = band_h, alpha = 0.65, inherit.aes = FALSE, width = 0.95
    )
    
    # << 不在色块上写字 >>
    if (isTRUE(tile_labels)) {
      p <- p + geom_text(
        data = band_df, aes(x = x, y = y, label = group),
        size = 4.2, inherit.aes = FALSE
      )
    }
    
    # << 给色块一个 legend，用于映射 group >>
    p <- p + scale_fill_manual(values = group_band_palette, name = "Group",
                               guide = if (isTRUE(show_group_legend)) "legend" else "none")
  } else if (group_indicator == "legend") {
    leg_df <- data.frame(group = levels(grp_fac), x = 1, y = y_limits[1])
    p <- p + geom_point(
      data = leg_df, aes(x = x, y = y, fill = group),
      inherit.aes = FALSE, shape = 21, size = 0, alpha = 0, show.legend = TRUE
    ) +
      scale_fill_manual(values = group_band_palette, name = "Group") +
      guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1)))
  }
  ggnewscale::new_scale_fill()
  
  if (!is.null(topN) && nrow(topN) > 0) {
    p <- p + geom_jitter(
      data = topN,
      aes(x = jittered_group, y = log2FoldChange, size = abs(log2FoldChange), fill = direction),
      width = 0, shape = 21, stroke = 1.2, color = "black", alpha = 0.85
    )
  }
  if (!is.null(user_df) && nrow(user_df)>0 && highlight_mode=="add") {
    p <- p + geom_jitter(
      data = user_df,
      aes(x = jittered_group, y = log2FoldChange),
      width = 0, shape = user_shape, size = user_size,
      stroke = user_stroke, color = "black", fill = user_fill, alpha = 0.95
    )
  }
  
  if (!is.null(topN) && nrow(topN) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = topN,
      aes(x = jittered_group, y = log2FoldChange, label = gene_name),
      box.padding = unit(0.4, "lines"),
      point.padding = unit(0.2, "lines"),
      segment.color = "black", force = 1.5, max.overlaps = Inf
    )
  }
  
  p <- p +
    geom_hline(yintercept = c(-lfc_thr, lfc_thr)) +
    scale_y_continuous(limits = y_limits) +
    scale_size(range = size_range) +
    scale_alpha(range = alpha_range) +
    scale_color_manual(values = dir_colors, name = "Direction") +
    scale_fill_manual(values = dir_colors, guide = "none") +
    scale_x_continuous(
      breaks = seq_len(nlevels(grp_fac)),
      labels = levels(grp_fac)
    ) +
    labs(x = "Group", y = "log2 Fold Change", title = "Multi-group volcano (x=Group, y=log2FC)") +
    theme_add +
    theme(
      panel.grid = element_blank(),
      legend.background = element_rect(color = "#969696"),
      axis.text.x = element_text(angle = x_text_angle, hjust = 1),
      plot.margin = margin(10, 10, 18, 10)
    )
  
  w <- max(10, nlevels(grp_fac) * width_per_group)
  ggsave(paste0(file_prefix, ".pdf"), p, width = w, height = height, units = "in", device = cairo_pdf)
  ggsave(paste0(file_prefix, ".png"), p, width = w, height = height, units = "in", dpi = 600)
  p
}
make_upset_input_from_de_list <- function(
    de_list,
    gene_col   = "Gene",
    lfc_col    = "log2FoldChange",
    padj_col   = "padj",
    which_set  = c("all","up","down"),
    lfc_thr    = 1,
    padj_thr   = 0.05,
    min_set_size = 1
){
  which_set <- match.arg(which_set)
  stopifnot(is.list(de_list), length(de_list) > 0)
  
  sets_list <- lapply(names(de_list), function(nm){
    df <- as.data.frame(de_list[[nm]])
    if (!gene_col %in% names(df) && !is.null(rownames(df))) {
      df <- tibble::rownames_to_column(df, gene_col)
    }
    stopifnot(all(c(gene_col, lfc_col, padj_col) %in% names(df)))
    ok <- is.finite(df[[padj_col]]) & df[[padj_col]] < padj_thr
    if (which_set == "up") {
      ok <- ok & (df[[lfc_col]] >  lfc_thr)
    } else if (which_set == "down") {
      ok <- ok & (df[[lfc_col]] < -lfc_thr)
    } else {
      ok <- ok & (abs(df[[lfc_col]]) > lfc_thr)
    }
    unique(df[[gene_col]][ok])
  })
  names(sets_list) <- names(de_list)
  
  # 过滤空集合 & 按集合大小排序
  sets_list <- sets_list[lengths(sets_list) >= min_set_size]
  if (length(sets_list) == 0) stop("没有通过阈值的集合（empty sets after filtering）")
  set_sizes <- sort(lengths(sets_list), decreasing = TRUE)
  set_order <- names(set_sizes)
  
  # 构造 presence/absence 矩阵
  all_genes <- unique(unlist(sets_list))
  upset_df  <- as.data.frame(
    sapply(sets_list, function(s) all_genes %in% s),
    row.names = all_genes
  )
  upset_df$gene <- rownames(upset_df)
  
  list(sets_list = sets_list, upset_df = upset_df, set_order = set_order)
}



# --- 工具：从逻辑矩阵提取 TopN “exclusive” 交集 ---
.compute_top_intersections <- function(upset_df, set_order, top_n = NULL){
  M <- as.data.frame(upset_df[, set_order, drop = FALSE])
  # 全转为逻辑
  M[] <- lapply(M, function(x){
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    if (is.character(x)) return(tolower(x) %in% c("1","t","true","y","yes"))
    if (is.factor(x)) return(as.character(x) %in% c("1","T","TRUE","true","y","yes"))
    as.logical(x)
  })
  if (!ncol(M)) return(list())
  
  # 每行的“精确组合” = 该行为 TRUE 的集合名
  comb_key <- apply(M, 1, function(row) {
    on <- set_order[isTRUE(row)]
    if (length(on) == 0) NA_character_ else paste(on, collapse = "&&")
  })
  comb_key <- comb_key[!is.na(comb_key)]
  if (!length(comb_key)) return(list())
  
  tab <- sort(table(comb_key), decreasing = TRUE)
  if (!is.null(top_n)) tab <- utils::head(tab, top_n)
  combos <- strsplit(names(tab), "&&", fixed = TRUE)
  # 确保纯 character
  combos <- lapply(combos, as.character)
  names(combos) <- NULL
  combos
}

# --- 主函数：优先 ComplexUpset，失败自动回退 UpSetR ---
plot_upset_auto <- function(
    de_input,
    style = c("nature","science"),
    set_palette = NULL,
    top_n_intersections = NULL,
    out_prefix = NULL,
    width = 11, height = 6.2, dpi = 450
){
  style <- match.arg(style)
  stopifnot(is.list(de_input),
            all(c("sets_list","upset_df","set_order") %in% names(de_input)))
  
  sets_list <- de_input$sets_list
  upset_df  <- de_input$upset_df
  set_order <- as.character(de_input$set_order)  # ★ 强制字符
  
  # 统一为 data.frame，避免 tibble 一些奇怪行为
  upset_df <- as.data.frame(upset_df)
  # 仅保留 set_order 里的列（存在的那部分）
  set_order <- intersect(set_order, colnames(upset_df))
  stopifnot(length(set_order) >= 2)
  
  # 列转逻辑（0/1/TRUE/FALSE）
  upset_df[set_order] <- lapply(upset_df[set_order], function(x){
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    if (is.character(x)) return(tolower(x) %in% c("1","t","true","y","yes"))
    if (is.factor(x))   return(as.character(x) %in% c("1","T","TRUE","true","y","yes"))
    as.logical(x)
  })
  
  # 默认色板
  if (is.null(set_palette)) {
    set_palette <- c('#8dd3c7','#ffffb3','#bebada','#fb8072',
                     '#80b1d3','#fdb462','#b3de69','#fccde5',
                     '#d9d9d9','#bc80bd','#ccebc5','#ffed6f')
  }
  set_palette <- rep(set_palette, length.out = length(set_order))
  names(set_palette) <- set_order
  
  # 主题
  nature_theme  <- ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Helvetica", color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.text = ggplot2::element_text(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      strip.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold")
    )
  science_theme <- ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Helvetica", color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
      axis.line = ggplot2::element_line(color = "black"),
      legend.title = ggplot2::element_text(face = "bold")
    )
  
  # 计算 TopN exclusive 组合（如未指定则让包自己排）
  intersections_arg <- NULL
  if (!is.null(top_n_intersections)) {
    intersections_arg <- .compute_top_intersections(upset_df, set_order, top_n = top_n_intersections)
    if (length(intersections_arg) == 0) intersections_arg <- NULL
  }
  
  # ========== 优先 ComplexUpset ==========
  if (requireNamespace("ComplexUpset", quietly = TRUE)) {
    p <- tryCatch({
      ComplexUpset::upset(
        upset_df,
        intersect     = set_order,
        intersections = intersections_arg,  # 可能为 NULL
        base_annotations = list(
          'Intersection size' = ComplexUpset::intersection_size(
            text = list(size = 3, vjust = -0.4),
            mapping = ggplot2::aes(fill = after_stat(degree))
          ) +
            {
              if (style == "nature")
                ggplot2::scale_fill_gradientn(
                  colours = c("#cbd5e8", "#80c5bb", "#1b9e77"),
                  name = "Intersection degree"
                )
              else
                ggplot2::scale_fill_gradientn(
                  colours = c("#e0e0e0", "#9e9e9e", "#212121"),
                  name = "Degree"
                )
            } +
            if (style == "nature") nature_theme else science_theme
        ),
        set_sizes = ComplexUpset::upset_set_size(
          text = list(size = 3, hjust = 1.02),
          mapping = ggplot2::aes(fill = after_stat(set))
        ) + ggplot2::scale_fill_manual(values = set_palette, guide = "none"),
        matrix = (
          ComplexUpset::upset_matrix(
            point_size = 2,
            matrix_dot_alpha = if (style == "nature") 0.9 else 1,
            stripes = ComplexUpset::upset_stripes(
              alpha = if (style == "nature") 0.08 else 0.05,
              colours = if (style == "nature") set_palette else rep("#000000", length(set_order))
            )
          ) + ggplot2::theme(panel.grid = ggplot2::element_blank())
        ),
        width_ratio = 0.25,
        min_size = 1,
        sort_intersections_by = "cardinality",
        wrap = FALSE
      ) + ggplot2::labs(title = "UpSet")
    }, error = function(e){
      message("ComplexUpset 绘图失败：", conditionMessage(e), "\n将自动回退到 UpSetR。")
      NULL
    })
    
    if (!is.null(p)) {
      if (!is.null(out_prefix)) {
        ggplot2::ggsave(paste0(out_prefix, ".pdf"), p, width = width, height = height, device = grDevices::cairo_pdf)
        ggplot2::ggsave(paste0(out_prefix, ".png"), p, width = width, height = height, dpi = dpi)
      }
      return(p)
    }
  }
  
  # ========== 回退 UpSetR ==========
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("既无 ComplexUpset 也无 UpSetR。请先安装：install.packages('ComplexUpset') 或 install.packages('UpSetR')")
  }
  fl <- UpSetR::fromList(sets_list)
  sets.bar.color <- set_palette[match(colnames(fl), names(set_palette))]
  sets.bar.color[is.na(sets.bar.color)] <- "#9e9e9e"
  
  grDevices::pdf(ifelse(is.null(out_prefix), "UpSetR.pdf", paste0(out_prefix, ".pdf")),
                 width = width, height = height)
  p <- UpSetR::upset(
    fl,
    sets = colnames(fl),
    keep.order = TRUE,
    order.by = "freq",
    decreasing = TRUE,
    nintersects = if (is.null(top_n_intersections)) 40 else top_n_intersections,
    main.bar.color = if (style=="nature") "#4F4F4F" else "#212121",
    sets.bar.color = sets.bar.color,
    text.scale = 1.4
  )
  grDevices::dev.off()
  invisible(p)
}



suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr); library(purrr)
  library(ggplot2)
  # 绘图包在用到的分支里再按需加载
})

# ---------- 1) 从 all_DE_results 生成集合 ----------
.guess_gene_col <- function(df){
  cand <- c("Gene","Symbol","gene","gene_name","GENE","SYMBOL","gene_id","ENSEMBL","ensembl")
  col <- cand[cand %in% names(df)]
  if (length(col)==0) stop("在 DE 结果里没找到基因列，请手动传 gene_col。")
  col[1]
}

make_deg_sets_from_all_de <- function(
    all_DE_results,
    groups      = NULL,       # 选择哪些组；NULL=用列表顺序全部
    which_set   = c("up","down","sig","all"), # 提取哪类集合
    lfc_thr     = 1,
    padj_thr    = 0.05,
    gene_col    = NULL,
    lfc_col     = "log2FoldChange",
    padj_col    = "padj",
    unique_genes = TRUE,
    drop_empty   = TRUE
){
  which_set <- match.arg(which_set)
  if (is.null(groups)) groups <- names(all_DE_results)
  if (is.null(groups) || length(groups)<2) stop("至少选择两个组来画韦恩图。")
  
  sets <- list()
  for (g in groups){
    df <- all_DE_results[[g]]
    if (is.null(df)) next
    gc <- if (is.null(gene_col)) .guess_gene_col(df) else gene_col
    if (!lfc_col %in% names(df) && which_set!="all") stop("缺少列：", lfc_col)
    if (!padj_col %in% names(df) && which_set!="all") stop("缺少列：", padj_col)
    
    genes <- dplyr::mutate(df,
                           .gene = .data[[gc]],
                           .lfc  = if (lfc_col %in% names(df)) as.numeric(.data[[lfc_col]]) else NA_real_,
                           .padj = if (padj_col %in% names(df)) as.numeric(.data[[padj_col]]) else NA_real_
    )
    
    sel <- switch(which_set,
                  up   = genes %>% filter(.padj < padj_thr, .lfc >  lfc_thr),
                  down = genes %>% filter(.padj < padj_thr, .lfc < -lfc_thr),
                  sig  = genes %>% filter(.padj < padj_thr, abs(.lfc) > lfc_thr),
                  all  = genes
    ) %>% pull(.gene) %>% as.character()
    
    if (unique_genes) sel <- unique(na.omit(sel))
    if (!(drop_empty && length(sel)==0)) sets[[g]] <- sel
  }
  if (length(sets)<2) stop("有效集合少于2个，无法绘制韦恩图。")
  sets
}

# ---------- 2) 计算所有交并集，写出 CSV ----------
compute_intersections_df <- function(gene_list){
  # gene x set 的逻辑矩阵
  all_genes <- sort(unique(unlist(gene_list)))
  mat <- sapply(gene_list, function(v) all_genes %in% v)
  colnames(mat) <- names(gene_list)
  inc <- as.data.frame(mat)
  inc$gene <- all_genes
  
  # 枚举非空组合
  sets <- names(gene_list); k <- length(sets)
  combos <- unlist(lapply(1:k, function(m) combn(sets, m, simplify = FALSE)), recursive = FALSE)
  
  res <- purrr::map_dfr(combos, function(cb){
    in_all  <- rowSums(inc[, cb, drop=FALSE]) == length(cb)
    out_oth <- if (length(setdiff(sets, cb))==0) rep(TRUE, nrow(inc)) else
      rowSums(inc[, setdiff(sets, cb), drop=FALSE]) == 0
    genes <- inc$gene[in_all & out_oth]
    tibble(combo = paste(cb, collapse = "&"),
           n = length(genes),
           genes = paste(genes, collapse = ";"))
  })
  res %>% arrange(desc(n))
}

# ---------- 3) 各类韦恩图绘制（自动根据集合数选择风格） ----------
plot_venn_suite <- function(
    gene_list,
    out_dir    = "Results/Venn",
    prefix     = "Venn",
    # 经典圆形 VennDiagram 配置
    vd_fill    = c("#80C8CF","#FBB9C0","#FEF392","#9DAFD1","#E790C0"),
    vd_alpha   = 0.5,
    # eulerr 配置（面积近似比例）
    eul_fill   = c("#80C8CF","#FBB9C0","#FEF392","#9DAFD1","#E790C0"),
    # ggVennDiagram 配置（3 组填充深浅）
    ggd_shape  = "301",
    width      = 4, height = 4, dpi = 300
){
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  nset <- length(gene_list)
  if (nset < 2 || nset > 5) stop("仅支持 2–5 个集合；更多集合建议 UpSet 图。")
  
  # 1) 输出集合清单
  for (nm in names(gene_list)){
    writeLines(gene_list[[nm]], file.path(out_dir, paste0(prefix,"_Set_", nm, ".txt")))
  }
  # 2) 输出交集表
  inter_df <- compute_intersections_df(gene_list)
  write.csv(inter_df, file.path(out_dir, paste0(prefix, "_Intersections.csv")), row.names = FALSE)
  
  # 3) 根据 nset 选择风格并保存
  if (nset <= 3){
    # 3a) 经典 VennDiagram
    if (!requireNamespace("VennDiagram", quietly = TRUE))
      stop("需要安装 VennDiagram 包。")
    vd <- VennDiagram::venn.diagram(
      x = gene_list,
      category.names = names(gene_list),
      filename = NULL, output = FALSE,
      disable.logging = TRUE,
      lwd = 2, lty = 'blank',
      fill = vd_fill[seq_len(nset)], alpha = vd_alpha,
      cex = 0.9, fontfamily = "sans",
      cat.cex = 0.9, cat.fontfamily = "sans", rotation = 1
    )
    gr <- cowplot::plot_grid(vd)
    ggsave(file.path(out_dir, paste0(prefix, "_Classic_", nset, "sets.pdf")),
           gr, width = width, height = height, device = cairo_pdf, dpi = dpi)
    ggsave(file.path(out_dir, paste0(prefix, "_Classic_", nset, "sets.png")),
           gr, width = width, height = height, dpi = dpi)
    
    # 3b) eulerr（面积近似比例）
    if (!requireNamespace("eulerr", quietly = TRUE))
      stop("需要安装 eulerr 包。")
    fit <- try(eulerr::euler(gene_list), silent = TRUE)
    if (!inherits(fit, "try-error")){
      p <- plot(fit,
                quantities = list(type = "counts", font = 1, cex = 0.85),
                fill = eul_fill[seq_len(nset)], edges = list(col = "black", lwd = 2), alpha = 0.8,
                labels = list(col = "black", font = 2, cex = 0.9)
      )
      ggplot2::ggsave(file.path(out_dir, paste0(prefix, "_Euler_", nset, "sets.pdf")),
                      plot = p, width = width, height = height, device = "pdf")
      ggplot2::ggsave(file.path(out_dir, paste0(prefix, "_Euler_", nset, "sets.png")),
                      plot = p, width = width, height = height, dpi = dpi)
    }
    
    # 3c) ggVennDiagram（3 组：填充色深浅代表大小）
    if (nset == 3){
      if (!requireNamespace("ggVennDiagram", quietly = TRUE))
        stop("需要安装 ggVennDiagram 包。")
      v   <- ggVennDiagram::Venn(gene_list)
      dat <- ggVennDiagram::process_data(v, shape_id = ggd_shape)
      p <- ggplot() +
        geom_polygon(aes(X, Y, fill = count, group = id), data = ggVennDiagram::venn_regionedge(dat)) +
        geom_path(aes(X, Y, color = id, group = id), data = ggVennDiagram::venn_setedge(dat), linewidth = 2, show.legend = FALSE) +
        geom_text(aes(X, Y, label = name), data = ggVennDiagram::venn_setlabel(dat)) +
        geom_text(aes(X, Y, label = label),
                  data = ggVennDiagram::venn_regionlabel(dat) |> dplyr::mutate(label = paste0(count, "\n", round(count/sum(count),2)*100, "%"))) +
        scale_fill_gradient(name = "Count", low = "white", high = "#C92027") +
        scale_color_manual(values = c("#91C73E", "#00B4C8", "#F16922")) +
        coord_equal() + theme_void()
      ggsave(file.path(out_dir, paste0(prefix, "_ggVennDiagram_3sets.pdf")),
             p, width = width, height = height, device = cairo_pdf, dpi = dpi)
      ggsave(file.path(out_dir, paste0(prefix, "_ggVennDiagram_3sets.png")),
             p, width = width, height = height, dpi = dpi)
    }
    
  } else if (nset == 4){
    # 4 组：ggvenn 清晰度最好
    if (!requireNamespace("ggvenn", quietly = TRUE))
      stop("需要安装 ggvenn 包。")
    p <- ggvenn::ggvenn(
      gene_list,
      show_percentage = FALSE,
      fill_color = vd_fill[1:4],
      fill_alpha = 0.5,
      stroke_linetype = "dashed",
      stroke_size = 1,
      text_size = 5,
      set_name_size = 5
    ) + theme_void()
    ggsave(file.path(out_dir, paste0(prefix, "_ggvenn_4sets.pdf")),
           p, width = width, height = height, device = cairo_pdf, dpi = dpi)
    ggsave(file.path(out_dir, paste0(prefix, "_ggvenn_4sets.png")),
           p, width = width, height = height, dpi = dpi)
    
  } else if (nset == 5){
    # 5 组：eulerr::venn 用逻辑矩阵
    if (!requireNamespace("eulerr", quietly = TRUE))
      stop("需要安装 eulerr 包。")
    # incidence matrix
    allg <- sort(unique(unlist(gene_list)))
    mat  <- as.data.frame(sapply(gene_list, function(v) allg %in% v))
    colnames(mat) <- names(gene_list); rownames(mat) <- allg
    fit <- eulerr::venn(mat)
    pal <- scales::alpha(RColorBrewer::brewer.pal(5, "Set2"), 0.3)
    p <- plot(fit,
              quantities = list(type = "counts", font = 1, cex = 0.7),
              fill = pal, edges = RColorBrewer::brewer.pal(5, "Set2"), alpha = 1,
              labels = list(col = "black", font = 1, cex = 1)
    )
    ggplot2::ggsave(file.path(out_dir, paste0(prefix, "_eulerr_5sets.pdf")),
                    plot = p, width = width, height = height, device = "pdf")
    ggplot2::ggsave(file.path(out_dir, paste0(prefix, "_eulerr_5sets.png")),
                    plot = p, width = width, height = height, dpi = dpi)
  }
  
  invisible(list(intersections = inter_df, sets = gene_list))
}



# ========= helper 1: 行Z分数并截断 =========
row_zscore <- function(mat, cap = 2) {
  mat <- as.matrix(mat)
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z[z >  cap] <-  cap
  z[z < -cap] <- -cap
  z
}

#align_by_meta
align_by_meta <- function(vst_mat, meta_df, dataset_name, order_levels, order_cells = NULL) {
  stopifnot(!is.null(colnames(vst_mat)),
            "SampleID" %in% names(meta_df),
            "Dataset"  %in% names(meta_df),
            "Group"    %in% names(meta_df))
  
  # 过滤到目标数据集并规范因子
  meta_sub <- meta_df %>%
    dplyr::filter(.data$Dataset == dataset_name) %>%
    dplyr::mutate(
      Group = toupper(Group),
      Group = factor(Group, levels = order_levels)
    )
  
  if ("cell_line" %in% names(meta_sub) && !is.null(order_cells)) {
    meta_sub <- meta_sub %>%
      dplyr::mutate(cell_line = factor(cell_line, levels = order_cells))
  } else if (!"cell_line" %in% names(meta_sub)) {
    meta_sub$cell_line <- NA_character_
  }
  
  # 只保留在矩阵中的样本；并提示未匹配的样本
  in_mat   <- meta_sub$SampleID %in% colnames(vst_mat)
  if (!all(in_mat)) {
    missing_ids <- meta_sub$SampleID[!in_mat]
    if (length(missing_ids)) {
      message("[", dataset_name, "] meta 中有 ", length(missing_ids),
              " 个样本在矩阵列名中找不到：",
              paste(utils::head(missing_ids, 5), collapse = ", "),
              if (length(missing_ids) > 5) " ..." else "")
    }
    meta_sub <- meta_sub[in_mat, , drop = FALSE]
  }
  if (!nrow(meta_sub)) stop("[", dataset_name, "] 没有任何样本能与矩阵列名匹配。")
  
  # 排序：cell_line -> Group -> replicate(若有) -> SampleID
  has_rep <- "replicate" %in% names(meta_sub)
  if (has_rep) {
    meta_sub <- meta_sub %>%
      dplyr::arrange(cell_line, Group, replicate, SampleID)
  } else {
    meta_sub <- meta_sub %>%
      dplyr::arrange(cell_line, Group, SampleID)
  }
  
  samp_order <- meta_sub$SampleID
  mat_out <- as.matrix(vst_mat[, samp_order, drop = FALSE])
  
  # 返回：基因 x 样本
  list(mat = mat_out, meta = meta_sub, samples = samp_order)
}


# ===================== 函数：多数据集环状热图（按原始列名绘制法） =====================
plot_multi_dataset_circos <- function(
    vst_list, meta_df, sel_genes,
    order_dataset, order_within,              # 保留入参以兼容，不再用于扇区拆分
    de_list = NULL,
    start_degree = 90, gap_between = 32, open_gap = NULL,
    cap_z = 2, col_fun = NULL,
    dataset_cols = NULL,                      # 每个 dataset 的颜色（可传入，不传用默认柔和色）
    cell_cols    = NULL,                      # 每个 cell 的颜色（可传入，不传用浅色系）
    group_cols   = c("P"="#4E79A7","R"="#E15759","Other"="#999999"),
    ring_order   = c("group","cell","dataset"),   # 默认 Group 最靠近热图，其次 Cell，最外是 Dataset
    paint_dataset = TRUE, paint_cell = TRUE, paint_group = TRUE,
    alpha_group = 0.55, alpha_cell = 0.50, alpha_dataset = 0.5,
    tile_track_height = 0.30,
    sample_cex = 0.55, gene_cex = 0.85,
    gene_tick_mm = 1.0, gene_text_offset_mm = 1.2,
    add_legends = TRUE, venn = TRUE, venn_size = 0.28,
    use_padj = TRUE, padj_cutoff = 0.05, lfc_cutoff = 1,
    outdir = getwd(), outfile = file.path(outdir, "circos_multi_ds_selected_genes_final.pdf"),
    font_family = "Times",
    order_cells = NULL,
    label_all_genes = FALSE,
    highlight_genes = NULL,        # 你手动给一组基因（字符向量）
    auto_highlight  = c("none","var","range"),  # 自动挑选：按总体方差/极差
    n_highlight     = 8,
    highlight_col   = "#000000",
    highlight_cex   = 1.18,
    highlight_tick_lwd = 2.2
){
  stopifnot(length(order_dataset) >= 1, !is.null(order_cells))
  auto_highlight <- match.arg(auto_highlight)
  
  row_zscore <- function(M, cap = 2){
    Z <- t(scale(t(M))); Z[is.na(Z)] <- 0
    Z[Z >  cap] <-  cap; Z[Z < -cap] <- -cap; Z
  }
  align_by_meta <- function(vst_mat, meta_df, dataset_name, order_cells){
    md <- meta_df %>% dplyr::filter(.data$Dataset == dataset_name)
    stopifnot(all(c("SampleID","Group","cell_line") %in% names(md)))
    md <- md %>%
      mutate(
        cell_line = factor(as.character(cell_line), levels = order_cells),
        Group     = factor(toupper(as.character(Group)), levels = c("P","R"))
      ) %>%
      arrange(cell_line, Group, SampleID)
    samp <- md$SampleID
    keep <- intersect(samp, colnames(vst_mat))
    md   <- md[match(keep, md$SampleID), , drop = FALSE]
    mat  <- vst_mat[, keep, drop = FALSE]
    list(mat = mat, meta = md, samples = keep)
  }
  adj_alpha <- function(cols, a) unname(vapply(cols, function(cc) adjustcolor(cc, alpha.f = a), character(1)))
  
  # Dataset：柔和（Okabe–Ito子集）
  if (is.null(dataset_cols)) {
    base_ds <- c("#86B6F2","#F2B6A0","#B6E3A8","#F2D28A")  # 柔和蓝/橙/绿/黄
    dataset_cols <- setNames(rep(base_ds, length.out = length(order_dataset)), order_dataset)
  } else dataset_cols <- dataset_cols[order_dataset]
  
  # Cell：浅色系（Set3/自定义）
  if (is.null(cell_cols)) {
    pal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
    cell_cols <- setNames(pal(length(order_cells)), order_cells)
  } else {
    miss <- setdiff(order_cells, names(cell_cols))
    if (length(miss)) {
      pal2 <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      cell_cols <- c(cell_cols, setNames(pal2(length(miss)), miss))
    }
    cell_cols <- cell_cols[order_cells]
  }
  
  # Group：用 Okabe–Ito 两色，适度降低饱和度由 alpha 控制（更不刺眼）
  group_cols <- group_cols[c("P","R","Other")]
  make_sample_by_gene <- function(vst_mat, meta_df, dataset_name, order_cells, sel_genes, cap=2){
    aln <- align_by_meta(vst_mat, meta_df, dataset_name, order_cells)
    rn <- rownames(aln$mat)
    map_upper_to_orig <- setNames(rn, toupper(rn))
    sel_upper <- toupper(sel_genes)
    sel_matched <- map_upper_to_orig[sel_upper]
    found <- !is.na(sel_matched)
    X <- matrix(NA_real_, nrow = ncol(aln$mat), ncol = length(sel_genes),
                dimnames = list(aln$samples, sel_genes))
    if (any(found)) {
      sub   <- aln$mat[sel_matched[found], aln$samples, drop = FALSE]
      sub_z <- row_zscore(sub, cap = cap)
      X[, sel_genes[found]] <- t(sub_z)
    }
    list(X = X, meta = aln$meta, samples = aln$samples)
  }
  
  pieces <- lapply(order_dataset, function(ds){
    stopifnot(ds %in% names(vst_list))
    make_sample_by_gene(vst_list[[ds]], meta_df, ds, order_cells, sel_genes, cap = cap_z)
  })
  names(pieces) <- order_dataset
  
  data_matrix <- do.call(rbind, lapply(pieces, `[[`, "X"))
  meta_comb   <- do.call(dplyr::bind_rows, lapply(pieces, `[[`, "meta"))
  
  # 扇区：Dataset|Cell；每扇区内部 P→R
  sector_levels <- unlist(lapply(order_dataset, function(ds){
    paste(ds, order_cells, sep = "|")
  }), use.names = FALSE)
  split_comb <- factor(paste(meta_comb$Dataset, as.character(meta_comb$cell_line), sep="|"),
                       levels = sector_levels)
  
  ord <- order(split_comb,
               factor(toupper(as.character(meta_comb$Group)), levels = c("P","R")),
               meta_comb$SampleID)
  data_matrix <- data_matrix[ord, , drop = FALSE]
  meta_comb   <- meta_comb[ord, , drop = FALSE]
  split_comb  <- droplevels(split_comb[ord])
  
  # gaps
  levs <- levels(split_comb)
  gaps <- c()
  for (ds in order_dataset) {
    k <- sum(startsWith(levs, paste0(ds, "|")))
    gaps <- c(gaps, rep(0, k - 1), gap_between)
  }
  if (length(gaps) > 0) {
    gaps[length(gaps)] <- if (is.null(open_gap)) gap_between else open_gap
  }
  
  if (is.null(highlight_genes) && auto_highlight != "none") {
    mfun <- switch(auto_highlight,
                   var   = function(v) stats::var(v, na.rm = TRUE),
                   range = function(v) diff(range(v, na.rm = TRUE)))
    s <- apply(data_matrix, 2, mfun)
    s[!is.finite(s)] <- 0
    highlight_genes <- names(sort(s, decreasing = TRUE))[seq_len(min(n_highlight, length(s)))]
  }
  
  hi_set <- if (is.null(highlight_genes)) character(0) else toupper(unique(highlight_genes))
  if(!is.null(outfile)) pdf(outfile, width = 10, height = 10)
  on.exit({ if(!is.null(outfile)) dev.off() }, add = TRUE)
  par(family = font_family)
  circos.clear()
  circos.par(start.degree = start_degree,
             gap.after    = gaps,
             track.margin = c(0, 0.01),
             cell.padding = c(0, 0, 0, 0))
  
  # 表达色标
  if (is.null(col_fun)) {
    col_fun <- circlize::colorRamp2(
      seq(-cap_z, cap_z, length.out = 100),
      colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100)
    )
  }
  
  # 主热图
  do.call(circos.heatmap, list(
    data_matrix,
    split         = split_comb,
    cluster       = FALSE,
    bg.border     = "black",
    bg.lwd        = 1,
    cell.border   = "white",
    cell.lwd      = 0.5,
    rownames.side = "outside",
    rownames.cex  = sample_cex,
    col           = col_fun,
    track.height  = tile_track_height
  ))
  
  # 基因名 + 高亮
  circos.track(
    track.index = get.current.track.index(),
    bg.border = NA,
    panel.fun = function(x,y){
      if (CELL_META$sector.numeric.index == length(levels(split_comb))) {
        cn <- colnames(data_matrix)
        is_hi <- toupper(cn) %in% hi_set
        n  <- length(cn)
        cell_h <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
        y_coords <- seq(CELL_META$cell.ylim[1] + cell_h/2,
                        CELL_META$cell.ylim[2] - cell_h/2, length.out = n)
        for (i in seq_len(n)) {
          circos.lines(
            c(CELL_META$cell.xlim[2],
              CELL_META$cell.xlim[2] + convert_x(gene_tick_mm, "mm")),
            c(y_coords[i], y_coords[i]),
            col = if (is_hi[i]) highlight_col else "black",
            lwd = if (is_hi[i]) highlight_tick_lwd else 1.3
          )
        }
        circos.text(
          rep(CELL_META$cell.xlim[2], n) + convert_x(gene_text_offset_mm, "mm"),
          y_coords, cn,
          cex = ifelse(is_hi, gene_cex * highlight_cex, gene_cex),
          adj = c(0, 0.5),
          facing = "inside",
          col = ifelse(is_hi, highlight_col, "black"),
          font = ifelse(is_hi, 2, 1)
        )
      }
    }
  )
  
  # ---------- 三个环：按 ring_order 从热图向外叠加 ----------
  add_ring <- function(kind){
    if (kind == "group" && isTRUE(paint_group)) {
      grp_fac <- factor(as.character(meta_comb$Group),
                        levels = intersect(c("P","R","Other"), names(group_cols)))
      levs_grp <- levels(grp_fac)
      # 颜色（同一基调，透明度协调）
      col_map  <- setNames(adj_alpha(group_cols[levs_grp], alpha_group),
                           seq_along(levs_grp))
      grp_mat  <- matrix(as.integer(grp_fac), ncol = 1,
                         dimnames = list(rownames(data_matrix), "Group"))
      do.call(circos.heatmap, list(
        grp_mat,
        split        = split_comb,
        col          = col_map,
        cluster      = FALSE,
        track.height = 0.04,
        cell.border  = NA,
        bg.border    = NA
      ))
    } else if (kind == "cell" && isTRUE(paint_cell)) {
      sector_cells <- sub("^.*\\|", "", levels(split_comb))
      cell_strip_col <- adj_alpha(cell_cols[sector_cells], alpha_cell)
      circos.track(
        ylim = c(0, 1), track.height = 0.045,
        bg.col = cell_strip_col, bg.border = NA,
        panel.fun = function(x, y) {}
      )
    } else if (kind == "dataset" && isTRUE(paint_dataset)) {
      sector_datasets <- sub("\\|.*$", "", levels(split_comb))
      dataset_strip_col <- adj_alpha(dataset_cols[sector_datasets], alpha_dataset)
      circos.track(
        ylim = c(0, 1), track.height = 0.035,
        bg.col = dataset_strip_col, bg.border = NA,
        panel.fun = function(x, y) {}
      )
    }
  }
  for (k in ring_order) add_ring(k)
  
  # 图例（只对启用的环）
  if (isTRUE(add_legends)) {
    lg_list <- list(
      ComplexHeatmap::Legend(
        title = "Expression z-score",
        col_fun = col_fun,
        at = round(seq(-cap_z, cap_z, length.out = 5), 2),
        title_position = "leftcenter-rot",
        title_gp = gpar(fontsize = 13),
        labels_gp = gpar(fontsize = 12)
      )
    )
    if (paint_dataset) {
      lg_list <- c(lg_list, list(
        ComplexHeatmap::Legend(
          title = "Dataset",
          legend_gp = gpar(fill = dataset_cols[order_dataset]),
          labels = order_dataset,
          title_position = "topleft",
          title_gp = gpar(fontsize = 12),
          labels_gp = gpar(fontsize = 11)
        )
      ))
    }
    if (paint_cell) {
      lg_list <- c(lg_list, list(
        ComplexHeatmap::Legend(
          title = "Cell line",
          legend_gp = gpar(fill = cell_cols[order_cells]),
          labels = order_cells,
          nrow = ceiling(length(order_cells)/2),
          title_position = "topleft",
          title_gp = gpar(fontsize = 12),
          labels_gp = gpar(fontsize = 10)
        )
      ))
    }
    if (paint_group) {
      used_grp <- intersect(names(group_cols), unique(as.character(meta_comb$Group)))
      lg_list <- c(lg_list, list(
        ComplexHeatmap::Legend(
          title = "Group (P/R)",
          legend_gp = gpar(fill = adj_alpha(group_cols[used_grp], alpha_group)),
          labels = used_grp,
          title_position = "topleft",
          title_gp = gpar(fontsize = 12),
          labels_gp = gpar(fontsize = 11)
        )
      ))
    }
    y_pos <- 0.90
    for (lg in lg_list) {
      draw(lg, x = unit(0.99,"npc") - unit(5,"mm"), y = unit(y_pos,"npc"), just=c("right","top"))
      y_pos <- y_pos - 0.18
    }
  }
  order_dataset<-order_dataset %>% rev()
  # Venn / Euler
  if (isTRUE(venn) && !is.null(de_list)) {
    de_sets <- lapply(order_dataset, function(ds){
      x <- de_list[[ds]]
      if (is.null(x)) return(character(0))
      if (is.data.frame(x)) {
        cols <- colnames(x)
        pcol <- if (use_padj && "padj" %in% cols) "padj"
        else if (!use_padj && "pval" %in% cols) "pval"
        else if ("padj" %in% cols) { message("[", ds, "] 指定的p列缺失，回退 padj"); "padj" }
        else if ("pval" %in% cols) { message("[", ds, "] 指定的p列缺失，回退 pval"); "pval" }
        else stop("de_list[['", ds, "']] 需包含 'padj' 或 'pval' 列")
        stopifnot(all(c("Gene","log2FC") %in% cols))
        x |>
          dplyr::filter(!is.na(.data[[pcol]]),
                        .data[[pcol]] <= padj_cutoff,
                        is.finite(log2FC), abs(log2FC) >= lfc_cutoff) |>
          dplyr::pull(Gene) |>
          unique()
      } else if (is.character(x)) unique(x)
      else stop("de_list[['", ds, "']] 必须是 data.frame 或 字符向量")
    })
    names(de_sets) <- order_dataset
    fit <- eulerr::euler(de_sets)
    pushViewport(viewport(width = venn_size, height = venn_size,
                          x = unit(0.5,"npc"), y = unit(0.5,"npc"),
                          just = c("center","center")))
    grid.draw(plot(
      fit,
      fills = list(fill = adj_alpha(dataset_cols[order_dataset], 0.38), alpha = 1),
      edges = list(col = "black", lty = 2, lwd = 2),
      labels = FALSE,
      quantities = TRUE
    ))
    grid.text("DE genes overlap", x = unit(0.5, "npc"), y = unit(-0.06, "npc"),
              gp = gpar(fontsize = 12, lineheight = 0.9))
    upViewport()
  }
  
  invisible(list(
    matrix = data_matrix, meta = meta_comb, split = split_comb,
    dataset_cols = dataset_cols, cell_cols = cell_cols, group_cols = group_cols,
    highlight_genes = highlight_genes
  ))
}

# 
# 
# vst_list <- list(Organoid = vst_org, PIK3IP1 = vst_pik)
# order_dataset <- c("Organoid","PIK3IP1")
# order_within  <- list(Organoid = ORDER_ORG, PIK3IP1 = ORDER_PIK)
# de_list <- list(Organoid = de_org, PIK3IP1 = de_pik)
# 
# outfile <- file.path(outdir, "circos_multi_ds_selected_genes_final.pdf")
# plot_multi_dataset_circos(
#   vst_list = vst_list,
#   meta_df  = all_meta,
#   sel_genes = c("ANO1","ANXA3","CALB1","CAV1","CAVIN2","CLIC3","CLTRN","CRIP2","FLNC","HYI","LGALS9","RAB27B","RHOB","VCAM1"),
#   order_dataset = order_dataset,
#   order_within  = order_within,
#   de_list = de_list,
#   gap_between = 32,
#   tile_track_height = 0.30,
#   sample_cex = 0.55,
#   gene_cex = 0.85,
#   venn = TRUE,
#   use_padj = TRUE,  # ← 用校正后的 padj；改成 FALSE 则用 pval
#   # outfile 可以不传，将默认写到 outdir/circos_multi_ds_selected_genes_final.pdf
# )















# ===================== prot_clean_qc_de.R =====================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
  library(purrr); library(ggplot2); library(patchwork); library(limma); library(scales)
})

# ---------- 工具 ----------
median_center <- function(mat){
  med <- apply(mat, 2, median, na.rm = TRUE)
  sweep(mat, 2, med, "-")
}
impute_lowabund_by_group <- function(M, grp, downshift = 1.8, width = 0.3, seed = 20250831){
  set.seed(seed)
  out <- M
  for(g in unique(grp)){
    cols <- which(grp == g)
    obs  <- as.numeric(M[, cols, drop = FALSE]); obs <- obs[is.finite(obs)]
    if(!length(obs)) next
    s <- stats::sd(obs, na.rm = TRUE); if(!is.finite(s) || s == 0) s <- stats::sd(M, na.rm = TRUE)
    mu <- stats::median(obs, na.rm = TRUE) - downshift * s
    sig<- width * s
    miss_idx <- which(is.na(M[, cols, drop = FALSE]), arr.ind = TRUE)
    if(nrow(miss_idx) > 0){
      repl <- stats::rnorm(nrow(miss_idx), mean = mu, sd = sig)
      for(i in seq_len(nrow(miss_idx))){
        out[miss_idx[i,"row"], cols[miss_idx[i,"col"]]] <- repl[i]
      }
    }
  }
  out
}
keep_by_presence <- function(mat, groups, min_rep = 2, both_sides = FALSE){
  idxP <- which(groups == "P"); idxR <- which(groups == "R")
  P_ok <- rowSums(!is.na(mat[, idxP, drop = FALSE])) >= min_rep
  R_ok <- rowSums(!is.na(mat[, idxR, drop = FALSE])) >= min_rep
  keep <- if (both_sides) (P_ok & R_ok) else (P_ok | R_ok)
  keep[is.na(keep)] <- FALSE
  keep
}

# ---------- 读入 & 仅保留成对细胞系 ----------
load_gene_and_meta <- function(gene_file = "proteomics_gene_abundance.tsv",
                               smeta_file = "proteomics_sample_metadata.tsv",
                               restrict_lines = NULL,   # 例：c("HPAF1","P1005","PL45","S007","SKPC1","SU86","Suit2")
                               require_pairs = TRUE){
  gene_df <- read_tsv(gene_file, show_col_types = FALSE) %>%
    filter(GeneSymbol != "") %>% distinct(GeneSymbol, .keep_all = TRUE)
  smeta   <- read_tsv(smeta_file, show_col_types = FALSE) %>%
    mutate(condition = toupper(condition)) %>% filter(condition %in% c("P","R"))
  
  # 只保留 P&R 都存在的细胞系（可选）
  if (require_pairs){
    paired <- smeta %>% count(cell_line, condition) %>%
      tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0) %>%
      filter(P > 0, R > 0) %>% pull(cell_line)
  } else {
    paired <- unique(smeta$cell_line)
  }
  if (!is.null(restrict_lines)) paired <- intersect(paired, restrict_lines)
  
  smeta <- smeta %>% filter(cell_line %in% paired)
  # gene_df 真实存在的样本列
  sample_cols <- intersect(smeta$sample_id, names(gene_df))
  smeta <- smeta %>% filter(sample_id %in% sample_cols)
  
  list(gene_df = gene_df, smeta = smeta, lines = sort(unique(smeta$cell_line)))
}

# ---------- (A) 每条细胞系做 R vs P 并合并（返回 DE_comb） ----------
run_de_per_line <- function(gene_df, smeta,
                            min_rep = 2, both_sides = FALSE,
                            log_offset = 1,
                            normalize = c("median","quantile","none"),
                            do_impute = TRUE, impute_downshift = 1.8, impute_width = 0.3, seed = 20250831,
                            eb_trend = TRUE, eb_robust = TRUE,
                            outdir = "DE_gene_per_line", write_files = TRUE){
  normalize <- match.arg(normalize)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  stopifnot(all(smeta$sample_id %in% names(gene_df)))
  
  per_line <- function(cl){
    submeta <- smeta %>% filter(cell_line == cl) %>% arrange(condition, replicate)
    cols <- submeta$sample_id
    M <- as.matrix(gene_df[, cols, drop = FALSE]); rownames(M) <- gene_df$GeneSymbol
    # 过滤
    keep <- keep_by_presence(M, groups = submeta$condition, min_rep = min_rep, both_sides = both_sides)
    M <- M[keep,, drop = FALSE]
    # 预处理
    Mlog <- log2(M + log_offset)
    if (normalize == "quantile") Mnor <- normalizeBetweenArrays(Mlog, method = "quantile") else
      if (normalize == "median") Mnor <- median_center(Mlog) else Mnor <- Mlog
    Mimp <- if (do_impute) impute_lowabund_by_group(Mnor, submeta$condition,
                                                    downshift = impute_downshift, width = impute_width, seed = seed) else Mnor
    # limma：R - P
    grp <- factor(submeta$condition, levels = c("P","R"))
    design <- model.matrix(~ 0 + grp); colnames(design) <- c("P","R")
    fit <- lmFit(Mimp, design)
    fit2 <- contrasts.fit(fit, makeContrasts(RminusP = R - P, levels = design))
    fit2 <- eBayes(fit2, trend = eb_trend, robust = eb_robust)
    tt  <- topTable(fit2, coef = "RminusP", number = Inf, sort.by = "none")
    out <- tibble(
      GeneSymbol = rownames(Mimp),
      log2FC     = tt$logFC,
      AveExpr    = tt$AveExpr,
      t          = tt$t,
      P.Value    = tt$P.Value,
      adj.P.Val  = p.adjust(tt$P.Value, method = "BH"),
      cell_line  = cl
    ) %>% arrange(adj.P.Val, desc(abs(log2FC)))
    if (write_files) write_tsv(out, file.path(outdir, paste0("DE_gene_", cl, "_R_vs_P.tsv")))
    out
  }
  
  res_list <- map(unique(smeta$cell_line), per_line)
  DE_comb  <- bind_rows(res_list)
  if (write_files) write_tsv(DE_comb, file.path(outdir, "DE_gene_Combined_Summary.tsv"))
  return(DE_comb)
}

# ---------- (B) 样本层面 QC ----------
qc_sample_level <- function(gene_df, smeta, outdir = "QC_sample",
                            log_offset = 1,
                            normalize = c("median","quantile","none"),
                            do_impute = TRUE,
                            impute_downshift = 1.8, impute_width = 0.3, seed = 20250831,
                            write_processed = FALSE,
                            processed_dir = "Processed_Expression"){
  normalize <- match.arg(normalize)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  samples <- smeta$sample_id
  stopifnot(all(samples %in% names(gene_df)))
  M <- as.matrix(gene_df[, samples, drop = FALSE])
  rownames(M) <- gene_df$GeneSymbol
  miss_pre <- is.na(M)
  
  # ---- 预处理 ----
  Mlog <- log2(M + log_offset)
  if (normalize == "quantile") {
    Mnor <- normalizeBetweenArrays(Mlog, method = "quantile")
  } else if (normalize == "median") {
    Mnor <- median_center(Mlog)
  } else {
    Mnor <- Mlog
  }
  grp <- smeta$condition
  if (do_impute) {
    Mimp <- impute_lowabund_by_group(Mnor, grp, downshift = impute_downshift, width = impute_width, seed = seed)
  } else {
    Mimp <- Mnor
  }
  
  # ---- 写出 QC 后的表达矩阵 + 分组（可选）----
  if (write_processed) {
    dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_tsv(tibble::rownames_to_column(as.data.frame(Mlog), "GeneSymbol"),
                     file.path(processed_dir, "expression_log2.tsv"))
    readr::write_tsv(tibble::rownames_to_column(as.data.frame(Mnor), "GeneSymbol"),
                     file.path(processed_dir, "expression_norm.tsv"))
    readr::write_tsv(tibble::rownames_to_column(as.data.frame(Mimp), "GeneSymbol"),
                     file.path(processed_dir, "expression_imputed.tsv"))
    readr::write_tsv(smeta %>% dplyr::select(sample_id, cell_line, condition, replicate),
                     file.path(processed_dir, "sample_metadata_used.tsv"))
  }
  
  # ---- 画 QC 图（保持原逻辑）----
  mk_long <- function(Mm, tag){
    as_tibble(Mm, rownames = "GeneSymbol") |>
      pivot_longer(-GeneSymbol, names_to="sample_id", values_to="value") |>
      left_join(smeta, by = "sample_id") |>
      mutate(stage = tag)
  }
  d_pre  <- mk_long(Mlog, "log2")
  d_post <- mk_long(Mimp, paste0("norm_", normalize, ifelse(do_impute, "_imputed","")))
  
  p_den <- ggplot(bind_rows(d_pre, d_post),
                  aes(x = value, group = sample_id, linetype = stage)) +
    geom_density() + facet_wrap(~stage, ncol = 1, scales = "free_y") +
    labs(title = "Density per sample (pre vs post)", x = "Intensity (log2 scale)") +
    theme_bw(base_size = 12)
  ggsave(file.path(outdir, "density_pre_post.png"), p_den, width = 8, height = 6, dpi = 150)
  
  p_box <- ggplot(bind_rows(d_pre, d_post),
                  aes(x = sample_id, y = value, fill = condition)) +
    geom_boxplot(outlier.shape = NA) + coord_flip() +
    facet_wrap(~stage, ncol = 1, scales = "free_y") +
    labs(title = "Boxplot per sample (pre vs post)", y = "Intensity (log2)") +
    theme_bw(base_size = 12)
  ggsave(file.path(outdir, "boxplot_pre_post.png"), p_box, width = 9, height = 9, dpi = 150)
  
  X <- t(Mimp)  # 样本 × 特征
  pc <- prcomp(X, center = TRUE, scale. = FALSE)
  pca_df <- data.frame(pc$x[,1:2]) %>%
    mutate(sample_id = rownames(.)) %>%
    left_join(smeta, by = "sample_id")
  p_pca <- ggplot(pca_df, aes(PC1, PC2, shape = condition, color = cell_line, label = sample_id)) +
    geom_point(size = 3) + labs(title = "PCA (post)") + theme_bw(base_size = 12)
  ggsave(file.path(outdir, "PCA_post.png"), p_pca, width = 8, height = 6, dpi = 150)
  
  if (requireNamespace("uwot", quietly = TRUE)) {
    set.seed(seed)
    um <- uwot::umap(X, n_neighbors = 15, min_dist = 0.3, metric = "euclidean")
    um_df <- as.data.frame(um) %>% setNames(c("UMAP1","UMAP2")) %>%
      mutate(sample_id = rownames(X)) %>% left_join(smeta, by = "sample_id")
    p_umap <- ggplot(um_df, aes(UMAP1, UMAP2, shape = condition, color = cell_line, label = sample_id)) +
      geom_point(size = 3) + labs(title = "UMAP (post)") + theme_bw(base_size = 12)
    ggsave(file.path(outdir, "UMAP_post.png"), p_umap, width = 8, height = 6, dpi = 150)
  }
  
  miss_rate <- colMeans(miss_pre, na.rm = TRUE)
  miss_df <- tibble(sample_id = names(miss_rate), miss_rate = miss_rate) %>%
    left_join(smeta, by = "sample_id")
  p_miss <- ggplot(miss_df, aes(x = reorder(sample_id, miss_rate), y = miss_rate, fill = condition)) +
    geom_col() + coord_flip() + scale_y_continuous(labels = scales::percent_format()) +
    labs(title = "Missing rate per sample (pre-imputation)", y = "Missing fraction") +
    theme_bw(base_size = 12)
  ggsave(file.path(outdir, "missing_rate_per_sample.png"), p_miss, width = 8, height = 6, dpi = 150)
  
  if (do_impute) {
    imputed_counts <- colSums(is.na(Mnor))
    imputed_frac <- imputed_counts / nrow(Mnor)
    imp_df <- tibble(sample_id = names(imputed_frac), imputed_frac = imputed_frac) %>%
      left_join(smeta, by = "sample_id")
    p_imp <- ggplot(imp_df, aes(x = reorder(sample_id, imputed_frac), y = imputed_frac, fill = condition)) +
      geom_col() + coord_flip() + scale_y_continuous(labels = scales::percent_format()) +
      labs(title = "Imputed fraction per sample", y = "Imputed fraction of features") +
      theme_bw(base_size = 12)
    ggsave(file.path(outdir, "imputed_fraction_per_sample.png"), p_imp, width = 8, height = 6, dpi = 150)
  }
  
  invisible(list(Mlog = Mlog, Mnor = Mnor, Mimp = Mimp))
}


# ---------- (C) 结果层面 QC（需要 DE_comb；本函数也能从文件读） ----------
qc_result_level <- function(DE_comb = NULL, de_file = NULL, outdir = "QC_result",
                            padj_thresh = 0.05, lfc_thresh = 0.58, top_n = 30,
                            lines_focus = NULL){
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(DE_comb)) {
    stopifnot(!is.null(de_file) && file.exists(de_file))
    DE_comb <- read_tsv(de_file, show_col_types = FALSE)
  }
  if (!is.null(lines_focus)) DE_comb <- DE_comb %>% filter(cell_line %in% lines_focus)
  stopifnot(all(c("GeneSymbol","log2FC","AveExpr","P.Value","adj.P.Val","cell_line") %in% names(DE_comb)))
  
  p_hist <- ggplot(DE_comb, aes(x = P.Value)) + geom_histogram(bins = 40) +
    labs(title = "P-value histogram (all tests)", x = "P-value") + theme_bw(12)
  ggsave(file.path(outdir, "pvalue_histogram.png"), p_hist, width = 7, height = 5, dpi = 150)
  
  for (cl in unique(DE_comb$cell_line)) {
    sub <- DE_comb %>% filter(cell_line == cl) %>%
      mutate(sig = adj.P.Val <= padj_thresh & abs(log2FC) >= lfc_thresh)
    p_vol <- ggplot(sub, aes(x = log2FC, y = -log10(P.Value), color = sig)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red3")) +
      labs(title = paste0("Volcano: ", cl), x = "log2FC (R - P)", y = "-log10(P)") + theme_bw(12)
    p_ma  <- ggplot(sub, aes(x = AveExpr, y = log2FC, color = sig)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red3")) +
      labs(title = paste0("MA: ", cl), x = "Average expression (A)", y = "M = log2FC") + theme_bw(12)
    ggsave(file.path(outdir, paste0("volcano_", cl, ".png")), p_vol, width = 6.5, height = 5, dpi = 150)
    ggsave(file.path(outdir, paste0("MA_", cl, ".png")),      p_ma,  width = 6.5, height = 5, dpi = 150)
  }
  
  dircount <- DE_comb %>%
    mutate(direction = case_when(
      adj.P.Val <= padj_thresh & log2FC >  0 ~ "Up",
      adj.P.Val <= padj_thresh & log2FC <= 0 ~ "Down",
      TRUE ~ NA_character_
    )) %>% filter(!is.na(direction)) %>%
    distinct(cell_line, GeneSymbol, direction) %>%
    count(GeneSymbol, direction, name = "n_lines") %>%
    group_by(GeneSymbol) %>% mutate(max_lines = max(n_lines)) %>% ungroup() %>%
    arrange(desc(max_lines), desc(n_lines)) %>%
    slice_head(n = top_n)
  
  p_cons <- ggplot(dircount, aes(x = reorder(GeneSymbol, max_lines), y = n_lines, fill = direction)) +
    geom_col(position = position_dodge(width = 0.8)) + coord_flip() +
    labs(title = "Top genes by cross-line consistency",
         y = "# of lines significant (adj.P ≤ 0.05)", x = "Gene") + theme_bw(12)
  ggsave(file.path(outdir, "top_gene_consistency.png"), p_cons, width = 8, height = 10, dpi = 150)
  
  invisible(TRUE)
}

# ---------- (D) 跨细胞系 主效应 ----------
de_main_effect <- function(gene_df, smeta, outdir = "DE_main_effect",
                           log_offset = 1,
                           normalize = c("median","quantile","none"),
                           do_impute = TRUE, impute_downshift = 1.8, impute_width = 0.3, seed = 20250831,
                           min_rep = 2, both_sides = FALSE,
                           method = c("duplicateCorrelation","fixed_cellline")){
  normalize <- match.arg(normalize)
  method    <- match.arg(method)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  stopifnot(all(smeta$sample_id %in% names(gene_df)))
  
  M <- as.matrix(gene_df[, smeta$sample_id, drop = FALSE]); rownames(M) <- gene_df$GeneSymbol
  
  # 按细胞线做 presence 过滤，取并集
  keep_any <- rep(FALSE, nrow(M))
  for (cl in unique(smeta$cell_line)) {
    idx <- which(smeta$cell_line == cl)
    K <- keep_by_presence(M[, idx, drop = FALSE], groups = smeta$condition[idx],
                          min_rep = min_rep, both_sides = both_sides)
    keep_any <- keep_any | K
  }
  M <- M[keep_any, , drop = FALSE]
  
  Mlog <- log2(M + log_offset)
  if (normalize == "quantile") Mnor <- normalizeBetweenArrays(Mlog, method = "quantile") else
    if (normalize == "median") Mnor <- median_center(Mlog) else Mnor <- Mlog
  grp <- factor(smeta$condition, levels = c("P","R"))
  Mimp <- if (do_impute) impute_lowabund_by_group(Mnor, grp, downshift = impute_downshift, width = impute_width, seed = seed) else Mnor
  
  if (method == "duplicateCorrelation") {
    design <- model.matrix(~ 0 + grp); colnames(design) <- c("P","R")
    block  <- factor(smeta$cell_line)
    corfit <- duplicateCorrelation(Mimp, design, block = block)
    fit <- lmFit(Mimp, design, block = block, correlation = corfit$consensus)
    fit2 <- contrasts.fit(fit, makeContrasts(RminusP = R - P, levels = design))
    fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  } else {
    design <- model.matrix(~ 0 + grp + factor(smeta$cell_line))
    colnames(design)[1:2] <- c("P","R")
    fit <- lmFit(Mimp, design)
    fit2 <- contrasts.fit(fit, makeContrasts(RminusP = R - P, levels = design))
    fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
  }
  
  tt <- topTable(fit2, coef = "RminusP", number = Inf, sort.by = "none")
  res <- tibble(
    GeneSymbol = rownames(Mimp),
    log2FC     = tt$logFC,
    AveExpr    = tt$AveExpr,
    t          = tt$t,
    P.Value    = tt$P.Value,
    adj.P.Val  = p.adjust(tt$P.Value, method = "BH")
  ) %>% arrange(adj.P.Val, desc(abs(log2FC)))
  
  write_tsv(res, file.path(outdir, "DE_main_effect_R_vs_P.tsv"))
  p_hist <- ggplot(res, aes(x = P.Value)) + geom_histogram(bins = 40) +
    labs(title = "Main-effect model: P-value histogram", x = "P-value") + theme_bw(12)
  ggsave(file.path(outdir, "pvalue_hist_main_effect.png"), p_hist, width = 7, height = 5, dpi = 150)
  
  invisible(res)
}

# ---------- (E) 一键总控：从“整理后的蛋白组学数据”出发 ----------
run_all_from_clean <- function(
    gene_file  = "proteomics_gene_abundance.tsv",
    smeta_file = "proteomics_sample_metadata.tsv",
    restrict_lines = NULL, require_pairs = TRUE,
    # 预处理参数
    log_offset = 1, normalize = c("median","quantile","none"),
    do_impute = TRUE, impute_downshift = 1.8, impute_width = 0.3, seed = 20250831,
    min_rep = 2, both_sides = FALSE,
    eb_trend = TRUE, eb_robust = TRUE,
    # 输出目录
    de_dir = "DE_gene_per_line", qc_sample_dir = "QC_sample",
    qc_result_dir = "QC_result", main_effect_dir = "DE_main_effect",
    # 新增：是否写出处理后的矩阵
    write_processed = TRUE, processed_dir = "Processed_Expression"
){
  L <- load_gene_and_meta(gene_file, smeta_file, restrict_lines, require_pairs)
  
  # 每线差异 & 合并
  DE_comb <- run_de_per_line(L$gene_df, L$smeta,
                             min_rep = min_rep, both_sides = both_sides,
                             log_offset = log_offset, normalize = match.arg(normalize),
                             do_impute = do_impute, impute_downshift = impute_downshift, impute_width = impute_width, seed = seed,
                             eb_trend = eb_trend, eb_robust = eb_robust,
                             outdir = de_dir, write_files = TRUE)
  
  # 样本层面 QC + 写出表达矩阵/分组
  qc_objs <- qc_sample_level(L$gene_df, L$smeta, outdir = qc_sample_dir,
                             log_offset = log_offset, normalize = match.arg(normalize),
                             do_impute = do_impute, impute_downshift = impute_downshift, impute_width = impute_width, seed = seed,
                             write_processed = write_processed, processed_dir = processed_dir)
  
  # 结果层面 QC
  qc_result_level(DE_comb, de_file = NULL, outdir = qc_result_dir,
                  padj_thresh = 0.05, lfc_thresh = 0.58, top_n = 30, lines_focus = L$lines)
  
  # 主效应
  res_main <- de_main_effect(L$gene_df, L$smeta, outdir = main_effect_dir,
                             log_offset = log_offset, normalize = match.arg(normalize),
                             do_impute = do_impute, impute_downshift = impute_downshift, impute_width = impute_width, seed = seed,
                             min_rep = min_rep, both_sides = both_sides, method = "duplicateCorrelation")
  
  invisible(list(DE_comb = DE_comb, qc = qc_objs, main_effect = res_main, lines = L$lines))
}

# res <- run_all_from_clean(
#   gene_file  = "proteomics_gene_abundance.tsv",
#   smeta_file = "proteomics_sample_metadata.tsv",
#   restrict_lines = c("HPAF1","P1005","PL45","S007","SKPC1","SU86","Suit2"),
#   require_pairs  = TRUE,
#   normalize = "median", do_impute = TRUE,
#   write_processed = TRUE,                       # ← 打开写出
#   processed_dir   = "Processed_Expression"      # ← 输出目录
# )


## =========================
## 2) Filtering, normalization, voom, ANOVA  (FUNCTION)
## =========================
voom_anova <- function(counts_mat,
                       colData_df,
                       group_col = "treatment",
                       fdr_thresh = 0.05,
                       trend = TRUE,
                       robust = TRUE,
                       use_quality_weights = FALSE,
                       # expression filters
                       filter_min_count = 10,         # counts
                       filter_min_total_count = 15,   # counts
                       filter_cpm = 1,                # CPM threshold
                       filter_min_samples = NULL,     # e.g. 6 (全样本数中至少6个样本 CPM>=1)
                       # effect-size / variability filters
                       min_delta_lfc = NA_real_,      # max(group mean) - min(group mean) on log2 scale
                       min_sd = NA_real_,             # SD across group means
                       # fallback caps
                       min_sig_genes = 200,
                       max_sig_genes = Inf) {
  
  counts_mat <- as.matrix(counts_mat); storage.mode(counts_mat) <- "numeric"
  stopifnot(ncol(counts_mat) == nrow(colData_df))
  group <- factor(colData_df[[group_col]])
  
  # ---- filtering ----
  dge <- edgeR::DGEList(counts = counts_mat)
  keep <- edgeR::filterByExpr(dge, group = group,
                              min.count = filter_min_count,
                              min.total.count = filter_min_total_count)
  if (!is.null(filter_min_samples)) {
    cpm_mat <- edgeR::cpm(dge)
    keep_cpm <- rowSums(cpm_mat >= filter_cpm) >= filter_min_samples
    keep <- keep & keep_cpm
  }
  dge <- dge[keep,, keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  message("After expression filtering: ", nrow(dge), " genes")
  
  # ---- voom + lmFit ----
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  v <- if (isTRUE(use_quality_weights)) {
    limma::voomWithQualityWeights(dge, design, plot = FALSE)
  } else {
    limma::voom(dge, design, plot = FALSE)
  }
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)
  
  # ---- global F (ANOVA-like) ----
  tabF <- limma::topTable(fit, coef = NULL, number = Inf, sort.by = "F", adjust.method = "BH")
  tabF$FDR <- tabF$adj.P.Val
  sig_flag <- tabF$FDR < fdr_thresh
  message("FDR <", fdr_thresh, " : ", sum(sig_flag), " genes")
  
  # ---- effect-size filters on group means ----
  prof_mu <- fit$coefficients[, levels(group), drop = FALSE]
  if (!all(rownames(prof_mu) == rownames(tabF))) {
    prof_mu <- prof_mu[rownames(tabF), , drop = FALSE]
  }
  
  if (!is.na(min_delta_lfc)) {
    delta <- apply(prof_mu, 1, function(x) diff(range(x)))
    sig_flag <- sig_flag & (delta >= min_delta_lfc)
    message("Δ(max-min) ≥ ", min_delta_lfc, " : ", sum(sig_flag), " genes")
  }
  if (!is.na(min_sd)) {
    sdg <- apply(prof_mu, 1, sd)
    sig_flag <- sig_flag & (sdg >= min_sd)
    message("SD ≥ ", min_sd, " : ", sum(sig_flag), " genes")
  }
  
  sig_genes <- rownames(tabF)[sig_flag]
  
  # fallback窗：至少/至多
  if (length(sig_genes) < min_sig_genes) {
    sig_genes <- rownames(tabF)[seq_len(min(min_sig_genes, nrow(tabF)))]
    message("Fallback min_sig_genes → kept ", length(sig_genes))
  }
  if (is.finite(max_sig_genes) && length(sig_genes) > max_sig_genes) {
    sig_genes <- sig_genes[seq_len(max_sig_genes)]
    message("Capped to max_sig_genes → kept ", length(sig_genes))
  }
  
  list(
    group_levels = levels(group),
    group_factor = group,
    dge   = dge,
    voom  = v,
    fit   = fit,
    tabF  = tabF,
    sig_genes = sig_genes
  )
}



## =========================
## 3+4) Build row-Z profiles + clustering (FUNCTION with named methods)
## =========================
profiles_and_cluster <- function(
    fit_obj,
    sig_genes,
    group_levels,
    method = c("knn_jaccard_louvain_merge","hierarchical_wardD2","kmeans","fuzzy_cmeans"),
    K_target = 8,
    # --- A: kNN 图参数 ---
    a_knn_frac = 0.10, a_k_min = 10, a_k_max = 50,
    # --- B: 层次聚类 ---
    use_corr_for_hclust = TRUE,
    # --- C: k-means ---
    km_nstart = 50, km_iter = 1000,
    # --- D: fuzzy c-means ---
    m_fuzz = 1.5,                   # fuzzifier (1.25–2 常用)
    # --- 小簇折叠（方案1：逐基因重分配） ---
    collapse_small = TRUE,          # 是否启用小簇处理
    min_cluster_size = 100,         # 小簇阈值：簇大小 < 该值则视为小簇
    size_ratio_trigger = 8,         # 触发条件之二：max(size)/min(size) >= 该值
    max_collapse_iter = 5,          # 重分配最多迭代次数
    verbose = TRUE
){
  method <- match.arg(method)
  
  ## ---------- 构建 row-Z profiles ----------
  prof <- fit_obj$coefficients[sig_genes, group_levels, drop = TRUE] %>% as.matrix()
  prof_z <- t(scale(t(prof))); prof_z[is.na(prof_z)] <- 0
  
  ## helpers
  compute_centroids <- function(mat, labels) {
    ks <- sort(unique(labels))
    cen <- do.call(rbind, lapply(ks, function(k){
      colMeans(mat[names(labels)[labels == k], , drop = FALSE], na.rm = TRUE)
    }))
    rownames(cen) <- ks
    cen
  }
  relabel_by_size <- function(labels) {
    tab <- sort(table(labels), decreasing = TRUE)
    map <- setNames(seq_along(tab), names(tab))
    setNames(as.integer(map[as.character(labels)]), names(labels))
  }
  collapse_small_clusters_per_gene <- function(labels, prof_z, min_size = 150, max_iter = 5) {
    labels <- as.integer(labels); names(labels) <- rownames(prof_z)
    for (iter in seq_len(max_iter)) {
      tab <- sort(table(labels), decreasing = TRUE)
      small <- as.integer(names(tab[tab < min_size]))
      big   <- as.integer(names(tab[tab >= min_size]))
      if (length(small) == 0L) break
      if (length(big) == 0L) { # 所有簇都很小：并到当前最大簇
        largest <- as.integer(names(tab)[1]); labels[] <- largest; break
      }
      cen <- compute_centroids(prof_z, labels)
      for (k in small) {
        idx <- names(labels)[labels == k]
        if (!length(idx)) next
        cors <- sapply(big, function(bk) {
          as.numeric(cor(t(prof_z[idx, , drop = FALSE]), cen[as.character(bk), ]))
        })
        if (is.null(dim(cors))) cors <- matrix(cors, ncol = 1, dimnames = list(idx, big))
        rownames(cors) <- idx; colnames(cors) <- as.character(big)
        best_big <- colnames(cors)[max.col(cors, ties.method = "first")]
        labels[idx] <- as.integer(best_big)
      }
    }
    labels <- relabel_by_size(labels)
    cen2 <- compute_centroids(prof_z, labels)
    rho2 <- numeric(nrow(prof_z)); names(rho2) <- rownames(prof_z)
    for (k in sort(unique(labels))) {
      ix <- names(labels)[labels == k]
      r  <- as.numeric(cor(t(prof_z[ix, , drop = FALSE]), cen2[as.character(k), ]))
      rho2[ix] <- ifelse(is.na(r), 0, r)
    }
    list(labels = labels, centroids = cen2, rho = rho2, membership01 = (rho2 + 1) / 2)
  }
  
  ## ---------- 初始聚类 ----------
  fuzzy_bonus <- list(fuzzy_membership = NULL, fuzzy_membership_max = NULL)
  
  if (method == "knn_jaccard_louvain_merge") {
    n <- nrow(prof_z)
    k <- min(a_k_max, max(a_k_min, round(a_knn_frac * n)))
    
    # kNN
    nn  <- FNN::get.knn(prof_z, k = k)
    nbr <- lapply(seq_len(n), function(i) nn$nn.index[i, ])
    ids <- rownames(prof_z)
    
    # SNN (Jaccard) 加权图
    edges <- vector("list", length = n * k); ctri <- 0L
    for (i in seq_len(n)) {
      ni <- nbr[[i]]
      for (j in ni) if (i < j) {
        shared <- length(intersect(ni, nbr[[j]]))
        if (shared > 0) {
          uni <- length(union(ni, nbr[[j]]))
          w <- shared / uni
          ctri <- ctri + 1L
          edges[[ctri]] <- c(ids[i], ids[j], w)
        }
      }
    }
    edges <- edges[seq_len(ctri)]
    edges_df <- as.data.frame(do.call(rbind, edges), stringsAsFactors = FALSE)
    colnames(edges_df) <- c("from","to","weight")
    edges_df$weight <- as.numeric(edges_df$weight)
    
    g <- igraph::graph_from_data_frame(edges_df, directed = FALSE, vertices = data.frame(name = ids))
    g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum"))
    
    # Leiden 优先，否则 Louvain
    cl_obj <- NULL
    if ("cluster_leiden" %in% ls("package:igraph")) {
      set.seed(1)
      cl_obj <- tryCatch(igraph::cluster_leiden(g, weights = igraph::E(g)$weight), error = function(e) NULL)
    }
    if (is.null(cl_obj)) {
      set.seed(1)
      cl_obj <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
    }
    cl_init <- igraph::membership(cl_obj) %>% as.integer() %>% setNames(igraph::V(g)$name)
    
    # 合并到目标 K（按质心相关，Ward.D2）
    cent_init <- compute_centroids(prof_z, cl_init)
    d_cent <- as.dist(1 - cor(t(cent_init)))
    hc <- hclust(d_cent, method = "ward.D2")
    lab_merge <- cutree(hc, k = K_target); names(lab_merge) <- rownames(cent_init)
    cl_final <- setNames(lab_merge[as.character(cl_init)], names(cl_init))
    cl_final <- relabel_by_size(cl_final)
    centroids_final <- compute_centroids(prof_z, cl_final)
    method_used <- sprintf("knn_jaccard_louvain_merge[k=%d]", k)
    
  } else if (method == "hierarchical_wardD2") {
    
    d <- if (use_corr_for_hclust) {
      as.dist(1 - cor(t(prof_z), method = "pearson"))
    } else {
      dist(prof_z, method = "euclidean")
    }
    hc <- hclust(d, method = "ward.D2")
    cl_final <- cutree(hc, k = K_target) %>% relabel_by_size()
    centroids_final <- compute_centroids(prof_z, cl_final)
    method_used <- "hierarchical_wardD2"
    
  } else if (method == "kmeans") {
    
    km <- kmeans(prof_z, centers = K_target, nstart = km_nstart, iter.max = km_iter, algorithm = "Lloyd")
    cl_final <- km$cluster %>% setNames(rownames(prof_z)) %>% relabel_by_size()
    centroids_final <- km$centers[sort(unique(cl_final)), , drop = FALSE]
    method_used <- "kmeans"
    
  } else if (method == "fuzzy_cmeans") {
    
    if (!requireNamespace("e1071", quietly = TRUE)) {
      stop("Please install package 'e1071' to use method = 'fuzzy_cmeans'.")
    }
    cm <- e1071::cmeans(prof_z, centers = K_target, m = m_fuzz,
                        iter.max = 100, dist = "euclidean", method = "cmeans")
    # 硬标签用于落盘/可视化分面
    cl_final <- cm$cluster %>% setNames(rownames(prof_z)) %>% relabel_by_size()
    # 质心（cmeans 的中心已在 prof_z 空间）
    centroids_final <- cm$centers
    rownames(centroids_final) <- sort(unique(cl_final))  # 对齐标签
    centroids_final <- centroids_final[as.character(sort(unique(cl_final))), , drop = FALSE]
    
    # 记录 fuzzy 隶属度（不参与 collapse；用于绘图可选）
    fuzzy_bonus$fuzzy_membership <- cm$membership[rownames(prof_z), , drop = FALSE]
    fuzzy_bonus$fuzzy_membership_max <- apply(fuzzy_bonus$fuzzy_membership, 1, max)
    
    method_used <- sprintf("fuzzy_cmeans[m=%.2f]", m_fuzz)
    
  } else {
    stop("Unknown method name.")
  }
  
  ## ---------- 初始 membership（与质心相关；保持与其他方法一致） ----------
  rho <- numeric(nrow(prof_z)); names(rho) <- rownames(prof_z)
  for (k in sort(unique(cl_final))) {
    ix <- names(cl_final)[cl_final == k]
    r  <- as.numeric(cor(t(prof_z[ix, , drop = FALSE]), centroids_final[as.character(k), ]))
    rho[ix] <- ifelse(is.na(r), 0, r)
  }
  membership01 <- (rho + 1) / 2
  
  sizes_initial <- sort(table(cl_final), decreasing = TRUE)
  if (verbose) {
    message("Initial cluster sizes: ", paste(sprintf("%s:%d", names(sizes_initial), sizes_initial), collapse = " "))
  }
  
  ## ---------- 触发条件：明显不均衡则逐基因重分配 ----------
  collapse_triggered <- FALSE
  if (isTRUE(collapse_small)) {
    sz <- as.integer(sizes_initial)
    size_ratio <- if (length(sz) >= 2) max(sz) / max(min(sz), 1L) else 1
    need_collapse <- (min(sz) < min_cluster_size) || (size_ratio >= size_ratio_trigger)
    
    if (need_collapse) {
      collapse_triggered <- TRUE
      if (verbose) {
        message(sprintf("Imbalance detected (min=%d, ratio=%.1f). Collapsing small clusters per gene (min_size=%d)...",
                        min(sz), size_ratio, min_cluster_size))
      }
      fix <- collapse_small_clusters_per_gene(cl_final, prof_z,
                                              min_size = min_cluster_size,
                                              max_iter = max_collapse_iter)
      cl_final       <- fix$labels
      centroids_final<- fix$centroids
      rho            <- fix$rho
      membership01   <- fix$membership01
      method_used    <- paste0(method_used, "+collapseSmall[min=", min_cluster_size, "]")
    }
  }
  
  sizes_final <- sort(table(cl_final), decreasing = TRUE)
  if (verbose && collapse_triggered) {
    message("Final cluster sizes:   ", paste(sprintf("%s:%d", names(sizes_final), sizes_final), collapse = " "))
  }
  
  # 返回
  c(
    list(
      method        = method_used,
      prof          = prof,
      prof_z        = prof_z,
      clusters      = cl_final,            # named integer vector (gene -> cluster)
      centroids     = centroids_final,     # matrix [K x 4]
      rho           = rho,                 # membership (raw Pearson r to centroid)
      membership01  = membership01,        # scaled to [0,1]
      sizes_initial = sizes_initial,
      sizes_final   = sizes_final,
      collapse_applied = collapse_triggered
    ),
    fuzzy_bonus    # 仅 fuzzy_cmeans 时包含；否则为 NULL
  )
}







## =========================
## RUN ENRICHMENT PER CLUSTER (KEGG/REACTOME/GO:BP/HALLMARK)
## compatible with msigdbr >= 10 (collection/subcollection)
## =========================
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr)
  library(clusterProfiler); library(msigdbr)
})

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr)
  library(clusterProfiler); library(msigdbr); library(org.Hs.eg.db)
})

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tidyr)
  library(clusterProfiler)
})

# ---- 用本地 GMT 构建 TERM2GENE，并对每簇做 ORA ----
run_cluster_enrichment_local <- function(
    prof_z, clusters_hc,
    gmt_dir,
    dbs = c("KEGG","REACTOME","GO:BP","HALLMARK"),
    top_per_db = 2,
    q_cutoff = 0.25,
    minGSSize = 10,
    pAdjustMethod = "BH",
    verbose = TRUE
){
  stopifnot(is.matrix(prof_z), all(!is.na(rownames(prof_z))))
  stopifnot(all(c("gene","cluster") %in% colnames(clusters_hc)))
  stopifnot(dir.exists(gmt_dir))
  
  ## gene universe（SYMBOL）与 cluster 向量
  universe <- unique(rownames(prof_z))
  cl_vec   <- setNames(clusters_hc$cluster, clusters_hc$gene)
  cl_vec   <- cl_vec[names(cl_vec) %in% universe]
  levs     <- sort(unique(cl_vec))
  
  ## 小工具：在目录里按正则找单个文件（取第一个匹配）
  find_gmt <- function(pattern) {
    files <- list.files(gmt_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
    if (length(files) == 0) return(NA_character_) else return(files[[1]])
  }
  
  ## 读取 GMT → TERM2GENE（gs_name, gene_symbol）
  read_t2g <- function(path) {
    if (is.na(path)) return(NULL)
    df <- tryCatch(clusterProfiler::read.gmt(path), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(NULL)
    colnames(df)[1:2] <- c("gs_name","gene_symbol")
    df
  }
  
  ## 为各 DB 寻找并读取文件
  t2g_list <- list()
  
  if ("KEGG" %in% dbs) {
    # 先找 *kegg* 的专用 GMT；若不存在，再考虑 c2.cp*.gmt 里筛选 CP:KEGG（这里你已经有专用 KEGG）
    p_kegg <- find_gmt("^c2\\.cp\\.kegg.*Hs\\.symbols\\.gmt$")
    if (verbose) message("[GMT] KEGG: ", ifelse(is.na(p_kegg), "NOT FOUND", basename(p_kegg)))
    t2g_list$KEGG <- read_t2g(p_kegg)
  }
  if ("REACTOME" %in% dbs) {
    p_rea <- find_gmt("^c2\\.cp\\.reactome.*Hs\\.symbols\\.gmt$")
    if (verbose) message("[GMT] REACTOME: ", ifelse(is.na(p_rea), "NOT FOUND", basename(p_rea)))
    t2g_list$REACTOME <- read_t2g(p_rea)
  }
  if ("GO:BP" %in% dbs) {
    p_gobp <- find_gmt("^c5\\.go\\.bp.*Hs\\.symbols\\.gmt$")
    if (verbose) message("[GMT] GO:BP: ", ifelse(is.na(p_gobp), "NOT FOUND", basename(p_gobp)))
    t2g_list[["GO:BP"]] <- read_t2g(p_gobp)
  }
  if ("HALLMARK" %in% dbs) {
    p_hall <- find_gmt("^h\\.all.*Hs\\.symbols\\.gmt$")
    if (verbose) message("[GMT] HALLMARK: ", ifelse(is.na(p_hall), "NOT FOUND", basename(p_hall)))
    t2g_list$HALLMARK <- read_t2g(p_hall)
  }
  
  ## 富集器（通用：用 enricher + TERM2GENE）
  do_enrich_t2g <- function(genes_sym, TERM2GENE, db_label){
    if (is.null(TERM2GENE) || nrow(TERM2GENE) == 0) return(NULL)
    genes_sym <- intersect(unique(genes_sym), universe)
    if (length(genes_sym) < 3) return(NULL)
    e <- tryCatch(
      enricher(gene = genes_sym,
               TERM2GENE = TERM2GENE,
               universe = universe,
               pAdjustMethod = pAdjustMethod,
               qvalueCutoff = q_cutoff,
               minGSSize = minGSSize),
      error = function(e) NULL
    )
    if (is.null(e)) return(NULL)
    df <- as.data.frame(e)
    if (!nrow(df)) return(NULL)
    df$db <- db_label
    df
  }
  
  ## 逐簇跑
  per_cluster_full <- vector("list", length(levs)); names(per_cluster_full) <- as.character(levs)
  per_cluster_top  <- vector("list", length(levs)); names(per_cluster_top)  <- as.character(levs)
  
  for (cl in levs) {
    genes_cl <- names(cl_vec)[cl_vec == cl]
    
    all_df <- bind_rows(lapply(intersect(names(t2g_list), dbs), function(db){
      do_enrich_t2g(genes_cl, TERM2GENE = t2g_list[[db]], db_label = db)
    }))
    
    per_cluster_full[[as.character(cl)]] <- all_df
    
    if (!is.null(all_df) && nrow(all_df)) {
      top_df <- all_df %>%
        filter(!is.na(p.adjust)) %>%
        mutate(Description = str_replace_all(Description, "_", " ")) %>%
        group_by(db) %>%
        arrange(p.adjust, .by_group = TRUE) %>%
        slice_head(n = top_per_db) %>%
        ungroup() %>%
        transmute(
          cluster = as.integer(cl),
          db,
          term = Description,
          Count,
          p.adjust,
          geneID,
          GeneRatio,
          BgRatio
        )
      per_cluster_top[[as.character(cl)]] <- top_df
    } else {
      per_cluster_top[[as.character(cl)]] <- NULL
    }
  }
  
  ## 合并输出（稳健处理空）
  combined_top  <- suppressWarnings(bind_rows(per_cluster_top,  .id = "cluster_chr"))
  if (is.null(combined_top) || nrow(combined_top) == 0) {
    combined_top <- tibble(cluster = integer(), db = character(), term = character(),
                           Count = integer(), p.adjust = numeric(),
                           geneID = character(), GeneRatio = character(), BgRatio = character())
  } else {
    combined_top <- combined_top %>%
      mutate(cluster = as.integer(cluster)) %>%
      select(-cluster_chr) %>%
      arrange(cluster, db, p.adjust)
  }
  
  combined_full <- suppressWarnings(bind_rows(per_cluster_full, .id = "cluster_chr"))
  if (!is.null(combined_full) && nrow(combined_full)) {
    combined_full <- combined_full %>%
      mutate(cluster = as.integer(cluster_chr)) %>%
      select(-cluster_chr)
  } else {
    combined_full <- tibble()
  }
  
  list(
    levs = levs,
    universe = universe,
    term2gene = t2g_list,
    per_cluster_full = per_cluster_full,
    per_cluster_top  = per_cluster_top,
    combined_full = combined_full,
    combined_top  = combined_top,
    source = "GMT(local)"
  )
}









plot_line_heat_enrich <- function(
    prof_z,            # 矩阵：基因 x 条件（行Z）
    clusters_hc,       # data.frame: gene, cluster, membership
    res_enr,           # run_cluster_enrichment*_ 的结果（含 $combined_top）
    outdir = getwd(),
    filename_base = "mfuzz_style_line_heat_dot",
    # 视觉参数
    cap_z = 2,
    max_rows_per_cluster = Inf,     # 为热图/折线展示每簇最多基因数（按 membership 取Top-N）
    size_mode = c("fixed_label","scaled","fixed"),
    fixed_size = 3.5,               # 等大小点半径（"fixed"/"fixed_label"）
    label_text_size = 2.6,          # 圆点内数字的字号（"fixed_label"）
    scaled_max_size = 5.5,          # 连续缩放的最大点径（"scaled"）
    width = 12, height = 9, dpi = 320,
    return_plot = TRUE
){
  size_mode <- match.arg(size_mode)
  
  stopifnot(is.matrix(prof_z), all(!is.na(rownames(prof_z))))
  stopifnot(all(c("gene","cluster","membership") %in% colnames(clusters_hc)))
  stopifnot("combined_top" %in% names(res_enr))
  
  suppressPackageStartupMessages({
    library(tidyverse); library(ggplot2); library(forcats); library(RColorBrewer); library(patchwork)
  })
  
  # ==== 0) 基础对象 ====
  cond_cols <- colnames(prof_z)
  levs <- sort(unique(clusters_hc$cluster))
  # 簇向量 & 隶属度（对齐 prof_z 行序）
  cl_vec   <- setNames(clusters_hc$cluster,    clusters_hc$gene)[rownames(prof_z)]
  memb_vec <- setNames(clusters_hc$membership, clusters_hc$gene)[rownames(prof_z)]
  stopifnot(!any(is.na(cl_vec)), !any(is.na(memb_vec)))
  
  # ==== 1) tidy 数据（折线/热图共用） ====
  exp_df <- as.data.frame(prof_z) %>%
    rownames_to_column("Identifier") %>%
    mutate(
      clusterindex = cl_vec[Identifier],
      maxMembership = memb_vec[Identifier]
    ) %>%
    pivot_longer(all_of(cond_cols), names_to = "Condition", values_to = "expression") %>%
    mutate(
      clusterindex = factor(clusterindex, levels = levs),
      Condition    = factor(Condition, levels = cond_cols),
      Time         = as.numeric(Condition)  # 1..4
    )
  
  # 可选：限制展示行数（Top-N membership）
  if (is.finite(max_rows_per_cluster)) {
    keep <- exp_df %>%
      distinct(Identifier, clusterindex, maxMembership) %>%
      group_by(clusterindex) %>%
      slice_max(order_by = maxMembership, n = max_rows_per_cluster, with_ties = FALSE) %>%
      ungroup()
    exp_df_plot <- exp_df %>% semi_join(keep, by = c("Identifier","clusterindex"))
  } else {
    exp_df_plot <- exp_df
  }
  
  # 簇中心
  centers <- exp_df_plot %>%
    group_by(clusterindex, Time) %>%
    summarise(Centre = mean(expression, na.rm = TRUE), .groups = "drop")
  
  # ==== 2) 调色板 ====
  SpatialColors <- colorRampPalette(rev(brewer.pal(11, "Spectral"))) # membership
  pal_expr <- rev(hcl.colors(11, "RdBu"))                            # 热图
  
  # ==== 3) 左侧折线 ====
  p1 <- ggplot(exp_df_plot, aes(x = Time, y = expression)) +
    geom_line(aes(group = Identifier, colour = maxMembership),
              linewidth = 0.25, alpha = 0.55) +
    geom_line(data = centers, aes(x = Time, y = Centre),
              colour = "grey30", linewidth = 0.9, inherit.aes = FALSE) +
    facet_grid(rows = vars(clusterindex), switch = "y") +
    scale_colour_gradientn("Membership", colours = SpatialColors(100),
                           limits = c(0,1),
                           breaks = c(0.1,0.3,0.5,0.7,0.9)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq_along(cond_cols), labels = cond_cols) +
    scale_y_continuous(limits = c(-cap_z, cap_z), expand = c(0, 0)) +
    theme_bw() +
    theme(
      legend.position = "top", legend.direction = "horizontal",
      legend.title = element_text(hjust = 0.5),
      axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
      strip.placement = "outside", strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(10, 5, 10, 10, "pt")
    ) +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5))
  
  # ==== 4) 中间热图 ====
  p2 <- exp_df_plot %>%
    ggplot(aes(x = Time, y = Identifier, fill = expression)) +
    geom_tile(color = NA) +
    scale_fill_gradientn("Expression",
                         colours = pal_expr,
                         limits  = c(-cap_z, cap_z),
                         breaks  = c(-cap_z, -cap_z/2, 0, cap_z/2, cap_z),
                         labels  = c(paste0("≤", -cap_z), -cap_z/2, 0, cap_z/2, paste0("≥", cap_z))) +
    scale_x_continuous(expand = c(0, 0), breaks = seq_along(cond_cols), labels = cond_cols) +
    facet_grid(rows = vars(clusterindex), scales = "free_y", space = "free_y") +
    theme_bw() +
    theme(
      legend.position = "top", legend.direction = "horizontal",
      legend.title = element_text(hjust = 0.5),
      axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
      strip.text = element_blank(), strip.background = element_blank(),
      panel.grid  = element_blank(),
      panel.border= element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(10, 10, 10, 0, "pt")
    ) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
  
  # ==== 5) 右侧富集点图 ====
  dot_df <- res_enr$combined_top %>%
    transmute(
      cluster = as.integer(cluster),
      cls     = db,
      go_term = term,
      Count   = Count,
      padj    = p.adjust
    ) %>%
    arrange(cluster, cls, padj) %>%
    group_by(cluster, cls) %>%
    slice_head(n = 2) %>%
    ungroup() %>%
    group_by(cluster) %>%
    mutate(
      Position = rev(seq_len(n())),
      ml10 = -log10(padj)
    ) %>%
    ungroup()
  
  # 对齐簇顺序，并把 cluster 前缀编码到术语中（方便 strip 关闭时识别）
  cluster_levels <- levels(exp_df_plot$clusterindex)
  dot_df <- dot_df %>%
    mutate(
      cls     = factor(cls, levels = c("KEGG","REACTOME","GO:BP","HALLMARK")),
      cluster = factor(cluster, levels = cluster_levels),
      go_term_lab = paste0("[C", as.integer(cluster), "] ", go_term)
    )
  
  db_cols <- c(
    KEGG     = "#ff7700",
    REACTOME = "#ac56af",
    "GO:BP"  = "#00b338",
    HALLMARK = "#0091EA"
  )
  
  base_p3 <- ggplot(dot_df,
                    aes(x = 0, y = fct_reorder(go_term_lab, Position),
                        fill = cls)) +
    scale_fill_manual("DB", values = db_cols, drop = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(position = "right") +
    facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top", legend.direction = "horizontal",
      legend.box = "horizontal", legend.justification = "left",
      legend.spacing.x = unit(0, "cm"),
      axis.title = element_blank(), axis.text.x = element_blank(),
      axis.ticks.length = unit(0, "pt"), axis.line = element_blank(),
      axis.text.y.right = element_text(size = 6),
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(10, 0, 0, 0, "pt")
    ) +
    guides(fill = guide_legend(title = "Database", title.position = "top", title.hjust = 0.5, nrow = 1))
  
  if (nrow(dot_df) == 0) {
    p3 <- ggplot() + theme_void() + ggtitle("No enriched terms")
  } else if (size_mode == "scaled") {
    p3 <- base_p3 +
      geom_point(aes(size = Count), shape = 21, colour = "grey20", stroke = 0.25) +
      scale_size_area(max_size = scaled_max_size, name = "Count") +
      guides(size = guide_legend(title.position = "top", title.hjust = 0.5))
  } else if (size_mode == "fixed_label") {
    p3 <- base_p3 +
      geom_point(shape = 21, colour = "grey20", stroke = 0.25,
                 size = fixed_size) +
      geom_text(aes(label = Count), size = label_text_size, fontface = 2) +
      guides(size = "none")
  } else { # "fixed"
    p3 <- base_p3 +
      geom_point(shape = 21, colour = "grey20", stroke = 0.25,
                 size = fixed_size) +
      guides(size = "none")
  }
  
  # ==== 6) 拼图并保存 ====
  combo <- p1 + p2 + p3 + plot_layout(widths = c(1, 1.05, 0.28))
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, paste0(filename_base, ".pdf")), combo, width = width, height = height)
  ggsave(file.path(outdir, paste0(filename_base, ".png")), combo, width = width, height = height, dpi = dpi, bg = "white")
  
  if (return_plot) return(invisible(list(p = combo, p1 = p1, p2 = p2, p3 = p3)))
}



run_line_heat_dot_pipeline <- function(
    # ---- INPUTS ----
    counts_mat, colData_df, group_col = "treatment",
    outdir, gmt_dir,
    seed = 1,
    
    # ---- step2: voom_anova 参数（与你已有函数一致）----
    fdr_thresh = 1e-3, trend = TRUE, robust = TRUE, use_quality_weights = TRUE,
    filter_min_count = 10, filter_min_total_count = 15,
    filter_cpm = 1, filter_min_samples = 8,
    min_delta_lfc = 1.0, min_sd = 0.4,
    min_sig_genes = 200, max_sig_genes = Inf,
    
    # ---- step3+4: 聚类参数 ----
    method = c("knn_jaccard_louvain_merge","hierarchical_wardD2","kmeans","fuzzy_cmeans"),
    K_target = 8,
    a_knn_frac = 0.10, a_k_min = 10, a_k_max = 50,      # kNN
    use_corr_for_hclust = TRUE,                          # hclust
    km_nstart = 50, km_iter = 1000,                      # kmeans
    m_fuzz = 1.6,                                        # fuzzy c-means
    collapse_small = TRUE, min_cluster_size = 150,
    size_ratio_trigger = 8, max_collapse_iter = 5,
    
    # ---- 富集参数（本地 GMT）----
    dbs = c("KEGG","REACTOME","GO:BP","HALLMARK"),
    top_per_db = 2, q_cutoff = 0.25, minGSSize = 10, pAdjustMethod = "BH",
    
    # ---- 绘图参数 ----
    filename_base = "mfuzz_style_line_heat_dot",
    cap_z = 2, max_rows_per_cluster = Inf,
    height_mode = c("equal","proportional"),           # 若你的 plot 函数不支持会自动忽略
    size_mode   = c("fixed_label","fixed","scaled"),
    fixed_size = 3.5, label_text_size = 2.6, scaled_max_size = 5.5,
    width = 14, height = 9, dpi = 320,
    
    # ---- 输出控制 ----
    save_tables = TRUE,
    return_objects = TRUE
){
  set.seed(seed)
  method <- match.arg(method)
  height_mode <- match.arg(height_mode)
  size_mode   <- match.arg(size_mode)
  
  if (missing(outdir) || is.null(outdir)) outdir <- getwd()
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  message("== Step 2: Filtering / voom / ANOVA ==")
  res_ana <- voom_anova(
    counts_mat, colData_df, group_col = group_col,
    fdr_thresh = fdr_thresh, trend = trend, robust = robust, use_quality_weights = use_quality_weights,
    filter_min_count = filter_min_count, filter_min_total_count = filter_min_total_count,
    filter_cpm = filter_cpm, filter_min_samples = filter_min_samples,
    min_delta_lfc = min_delta_lfc, min_sd = min_sd,
    min_sig_genes = min_sig_genes, max_sig_genes = max_sig_genes
  )
  message("Genes kept for profiling: ", length(res_ana$sig_genes))
  
  message("== Step 3+4: Build row-Z + Clustering ==")
  res_clust <- profiles_and_cluster(
    fit_obj = res_ana$fit,
    sig_genes = res_ana$sig_genes,
    group_levels = res_ana$group_levels,
    method = method, K_target = K_target,
    a_knn_frac = a_knn_frac, a_k_min = a_k_min, a_k_max = a_k_max,
    use_corr_for_hclust = use_corr_for_hclust,
    km_nstart = km_nstart, km_iter = km_iter,
    m_fuzz = m_fuzz,
    collapse_small = collapse_small, min_cluster_size = min_cluster_size,
    size_ratio_trigger = size_ratio_trigger, max_collapse_iter = max_collapse_iter,
    verbose = TRUE
  )
  
  # 汇总 cluster 表
  clusters_df <- tibble(
    gene       = rownames(res_clust$prof_z),
    cluster    = as.integer(res_clust$clusters[gene]),
    rho        = res_clust$rho[gene],
    membership = res_clust$membership01[gene]
  ) %>% arrange(cluster, desc(membership), gene)
  
  centroids_df <- as.data.frame(res_clust$centroids) %>%
    rownames_to_column("cluster") %>%
    mutate(cluster = as.integer(cluster))
  
  # ---- 保存表格 ----
  if (isTRUE(save_tables)) {
    write.csv(res_ana$tabF, file.path(outdir, "ANOVA_all_genes_topTableF.csv"), row.names = TRUE)
    write.csv(clusters_df,  file.path(outdir, sprintf("clusters_%s_K%d.csv", res_clust$method, K_target)), row.names = FALSE)
    write.csv(centroids_df, file.path(outdir, sprintf("centroids_%s_K%d.csv", res_clust$method, K_target)), row.names = FALSE)
    write.csv(as.data.frame(table(clusters_df$cluster)),
              file.path(outdir, sprintf("cluster_sizes_%s_K%d.csv", res_clust$method, K_target)), row.names = FALSE)
  }
  
  message("== Step 5: Enrichment with local GMT ==")
  stopifnot(dir.exists(gmt_dir))
  clusters_hc <- clusters_df %>% dplyr::select(gene, cluster, membership)
  
  res_enr <- run_cluster_enrichment_local(
    prof_z      = res_clust$prof_z,
    clusters_hc = clusters_hc,
    gmt_dir     = gmt_dir,
    dbs         = dbs,
    top_per_db  = top_per_db,
    q_cutoff    = q_cutoff,
    minGSSize   = minGSSize,
    pAdjustMethod = pAdjustMethod,
    verbose = TRUE
  )
  
  if (isTRUE(save_tables)) {
    write.csv(res_enr$combined_top,
              file.path(outdir, "Cluster_Enrichment_topN_byDB.csv"),
              row.names = FALSE)
  }
  
  message("== Step 6: Plot line + heatmap + enrichment dots ==")
  # 兼容不同版本的 plot_line_heat_enrich（有的无 height_mode 参数）
  plot_args <- list(
    prof_z = res_clust$prof_z,
    clusters_hc = clusters_hc,
    res_enr = res_enr,
    outdir = outdir,
    filename_base = filename_base,
    cap_z = cap_z,
    max_rows_per_cluster = max_rows_per_cluster,
    size_mode = size_mode,
    fixed_size = fixed_size,
    label_text_size = label_text_size,
    scaled_max_size = scaled_max_size,
    width = width, height = height, dpi = dpi,
    return_plot = FALSE
  )
  if ("height_mode" %in% names(formals(plot_line_heat_enrich))) {
    plot_args$height_mode <- height_mode
  }
  do.call(plot_line_heat_enrich, plot_args)
  
  message("DONE. Method = ", res_clust$method,
          " | genes clustered = ", nrow(clusters_df),
          " | clusters = ", length(unique(clusters_df$cluster)))
  
  if (isTRUE(return_objects)) {
    return(invisible(list(
      res_ana = res_ana,
      res_clust = res_clust,
      clusters_df = clusters_df,
      centroids_df = centroids_df,
      res_enr = res_enr,
      outdir = outdir,
      files = list(
        topF = file.path(outdir, "ANOVA_all_genes_topTableF.csv"),
        clusters = file.path(outdir, sprintf("clusters_%s_K%d.csv", res_clust$method, K_target)),
        centroids = file.path(outdir, sprintf("centroids_%s_K%d.csv", res_clust$method, K_target)),
        sizes = file.path(outdir, sprintf("cluster_sizes_%s_K%d.csv", res_clust$method, K_target)),
        enr_top = file.path(outdir, "Cluster_Enrichment_topN_byDB.csv"),
        fig_pdf = file.path(outdir, paste0(filename_base, ".pdf")),
        fig_png = file.path(outdir, paste0(filename_base, ".png"))
      )
    )))
  }
}


# 加载必要的R包
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

theme_publication <- function(base_size=12, show_grid_lines = FALSE){ # <-- 新增参数
  th <- theme_classic(base_size=base_size) +
    theme(
      plot.title = element_text(face="bold", size=rel(1.1), hjust=0, margin=margin(b=5)),
      axis.title = element_text(face="bold", size=rel(1)),
      axis.text = element_text(size=rel(0.9), color="black"),
      axis.line = element_line(color="black", linewidth=0.5),
      axis.ticks = element_line(color="black", linewidth=0.5),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # --- 根据参数决定是否显示网格线 ---
  if (show_grid_lines) {
    th <- th + theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4)
    )
  } else {
    th <- th + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }
  
  return(th)
}



universal_rank_plot <- function(data, 
                                id_col, 
                                value_col, 
                                label_items = character(), 
                                top_n_each_side = 0,
                                title = "Rank Plot",
                                xlab = "Rank",
                                ylab = "Value",
                                show_connector_line = F,
                                add_threshold_lines = NULL,
                                colors = c(pos = "#B2182B", neg = "#2166AC"),
                                base_point_size = 1.2,
                                highlight_point_size = 3,
                                highlight_fill = "gold",
                                highlight_color = "black",
                                label_font_face = "plain",
                                show_grid_lines = T, # <-- 新增参数
                                base_font_size = 12) { # <-- 新增参数
  
  # ... (数据准备和确定标注点的代码与上一版完全相同) ...
  df_ranked <- data %>%
    select(id = {{id_col}}, value = {{value_col}}) %>%
    filter(is.finite(value)) %>%
    arrange(desc(value)) %>%
    mutate(
      rank = row_number(),
      side = if_else(value >= 0, "pos", "neg")
    )
  
  labs <- df_ranked %>% 
    filter(
      toupper(id) %in% toupper(label_items) |
        rank <= top_n_each_side |
        rank > (nrow(.) - top_n_each_side)
    ) %>% 
    distinct(id, .keep_all = TRUE)
  
  # --- 绘图 ---
  p <- ggplot(df_ranked, aes(x = rank, y = value))
  
  if (show_connector_line) {
    p <- p + geom_line(linewidth = 0.3, colour = "grey80")
  }
  
  p <- p + geom_point(aes(colour = side), size = base_point_size, alpha = 0.7) +
    scale_color_manual(values = colors)
  
  if (!is.null(add_threshold_lines)) {
    p <- p + geom_hline(yintercept = add_threshold_lines, linetype = "dashed",
                        colour = "grey40", linewidth = 0.6)
  }
  
  if (nrow(labs) > 0) {
    p <- p + geom_point(data = labs, shape = 21, size = highlight_point_size, 
                        fill = highlight_fill, colour = highlight_color, stroke = 0.7)
    
    p <- p + geom_text_repel(
      data = labs, aes(label = id), size = 3.5, fontface = label_font_face,
      box.padding = 0.4, point.padding = 0.5, segment.color = 'grey50',
      min.segment.length = 0, max.overlaps = Inf, seed = 7
    )
  }
  
  # --- 应用主题和标签 (将参数传递给主题函数) ---
  p <- p + theme_publication(base_size = base_font_size, 
                             show_grid_lines = show_grid_lines) + # <-- 传递参数
    labs(x = xlab, y = ylab, title = title)
  
  return(p)
}



suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(readr)
  library(writexl)
  library(scales)
})

# ---------- Core function (minimal-change + sample gap + direction legend) ----------
hallmark_enrichment_dotplot <- function(
    de_all,
    obj_list = NULL,                         # optional: list of Seurat objects per sample
    gmt_path = "E:/modifiedcodeR/h.all.v2024.1.Hs.symbols.gmt",
    q_cut = 0.05,
    fc_cut = 0.25,
    top_terms_per_group = 10,
    wrap_width = 40,
    x_order_by_sample = NULL,               # optional: custom order of samples on x-axis
    save_prefix = NULL,                     # optional: if set, save csv/xlsx/pdf with this prefix
    fig_width = 12,
    fig_height = 9,
    # --- layout & background switches ---
    facet_by_sample = TRUE,                 # 分面显示样本，产生明显间隔
    shade_direction = TRUE,                 # 面板内用淡色背景区分方向
    dir_cols = NULL,                        # 必填（开启 shade 或方向图例时）：c("Basal-like"="#e31a1c","Classical"="#1f78b4")
    shade_alpha = 0.06,                     # 背景透明度
    panel_spacing_x = 1.2,                  # 样本面板水平间距
    # --- NEW: 把方向放入图例而不是X轴文本 ---
    direction_to_legend = TRUE              # TRUE: 隐藏X轴方向文字，并添加“Direction”图例
){
  # ---- 0) Basic checks ----
  req_cols <- c("sample","comparison","gene","up_in","avg_log2FC","p_val","p_val_adj","pct.1","pct.2","delta_pct")
  if (!all(req_cols %in% colnames(de_all))) {
    stop("de_all must contain columns: ", paste(req_cols, collapse=", "))
  }
  if (!file.exists(gmt_path)) {
    stop("GMT file not found: ", gmt_path)
  }
  # dir_cols 必须提供（当需要背景或方向图例时）
  if (isTRUE(shade_direction) || isTRUE(direction_to_legend)) {
    if (is.null(dir_cols) || !all(c("Basal-like","Classical") %in% names(dir_cols))) {
      stop("Please provide 'dir_cols' like c('Basal-like'='#e31a1c','Classical'='#1f78b4').")
    }
  }
  
  TERM2GENE <- clusterProfiler::read.gmt(gmt_path)
  
  run_enrich <- function(genes_up, universe = NULL) {
    if (length(genes_up) < 5) return(NULL)
    enricher(
      gene         = genes_up,
      TERM2GENE    = TERM2GENE,
      universe     = universe,
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2
    )
  }
  
  # ---- 1) Build per-sample, per-direction gene lists ----
  samples <- unique(de_all$sample)
  if (!is.null(x_order_by_sample)) {
    samples <- intersect(x_order_by_sample, samples)
  }
  
  enrich_list <- list()
  for (sm in samples) {
    de_sm <- de_all %>% filter(sample == sm)
    
    genes_basal_up <- de_sm %>%
      filter(up_in == "Basal-like", p_val_adj < q_cut, avg_log2FC >  fc_cut) %>%
      pull(gene) %>% unique()
    
    genes_class_up <- de_sm %>%
      filter(up_in == "Classical",  p_val_adj < q_cut, avg_log2FC < -fc_cut) %>%
      pull(gene) %>% unique()
    
    universe <- NULL
    if (!is.null(obj_list) && !is.null(obj_list[[sm]])) {
      universe <- tryCatch(rownames(obj_list[[sm]]), error = function(e) NULL)
    }
    
    ek_basal <- run_enrich(genes_basal_up, universe)
    ek_class <- run_enrich(genes_class_up, universe)
    
    tidy_one <- function(ek, direction) {
      if (is.null(ek)) return(NULL)
      as.data.frame(ek) %>%
        mutate(sample = sm, direction = direction) %>%
        dplyr::select(sample, direction, ID, Description, GeneRatio, BgRatio,
                      pvalue, p.adjust, qvalue, Count, everything())
    }
    
    enrich_list[[sm]] <- bind_rows(
      tidy_one(ek_basal, "Basal-like"),
      tidy_one(ek_class, "Classical")
    )
  }
  
  enrich_all <- bind_rows(enrich_list[!vapply(enrich_list, is.null, logical(1))])
  if (nrow(enrich_all) == 0) {
    stop("No enrichment results. Consider relaxing q_cut/fc_cut or check DE input.")
  }
  
  # ---- 2) Build dotplot data ----
  gene_ratio_num <- function(x) {
    sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
  }
  
  plot_df <- enrich_all %>%
    group_by(sample, direction) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_terms_per_group) %>%
    ungroup() %>%
    mutate(
      Group        = paste(sample, direction, sep = " | "),
      P_adj        = p.adjust,
      GeneRatioNum = gene_ratio_num(GeneRatio),
      Pathway      = pretty_hallmark(Description, width = wrap_width)
    )
  
  dir_order <- c("Basal-like","Classical")
  plot_df$direction <- factor(plot_df$direction, levels = dir_order)
  
  if (!is.null(x_order_by_sample)) {
    plot_df$sample <- factor(plot_df$sample, levels = x_order_by_sample)
  }
  
  plot_df <- plot_df %>%
    arrange(sample, direction, P_adj) %>%
    mutate(Group = factor(Group, levels = unique(Group)))
  
  path_order <- plot_df %>%
    group_by(Pathway) %>%
    summarise(min_padj = min(P_adj, na.rm = TRUE), .groups = "drop") %>%
    arrange(min_padj) %>%
    pull(Pathway)
  plot_df$Pathway <- factor(plot_df$Pathway, levels = rev(path_order))
  
  # ---- 3) Plot ----
  fill_scale <- scale_fill_gradientn(
    colors = c("#B30000", "#F46D43", "#FEE08B", "#ABD9E9", "#4575B4"),
    values = scales::rescale(c(0, 0.01, 0.02, 0.05)),
    limits = c(0, 0.05), oob = scales::squish,
    name = "FDR (adj. p)",
    breaks = c(0.001, 0.01, 0.02, 0.05),
    labels = c("0.001", "0.01", "0.02", "0.05")
  )
  
  size_scale <- scale_size_continuous(
    range = c(3.8, 11),
    breaks = scales::pretty_breaks(n = 3),
    name = "Gene ratio"
  )
  
  base_theme <- theme_minimal(base_size = 14) +
    theme(
      panel.border       = element_rect(fill = NA, color = "black", linewidth = 0.8),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.title.x       = element_blank(),
      axis.title.y       = element_blank(),
      axis.text.x        = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      legend.title       = element_text(size = 12),
      legend.text        = element_text(size = 11),
      plot.margin        = margin(t = 10, r = 12, b = 10, l = 10)
    )
  
  guide_setup <- guides(
    fill  = guide_colourbar(order = 1, reverse = TRUE, ticks = FALSE,
                            barheight = grid::unit(12, "lines"),
                            barwidth  = grid::unit(0.7,  "lines"),
                            frame.colour = "black", frame.linewidth = 0.8),
    size  = guide_legend(order = 2,
                         override.aes = list(shape = 21, fill = NA, color = "black", stroke = 1.2))
    # color (Direction) guide 在下面根据需要再加（order = 3）
  )
  
  if (isTRUE(facet_by_sample)) {
    p <- ggplot(plot_df, aes(x = direction, y = Pathway)) +
      { if (isTRUE(shade_direction)) {
        list(
          # Basal-like 背景
          geom_rect(
            data = data.frame(sample = unique(plot_df$sample)),
            aes(ymin = -Inf, ymax = Inf),
            xmin = 0.5, xmax = 1.5,
            inherit.aes = FALSE,
            fill = unname(dir_cols["Basal-like"]),
            alpha = shade_alpha, color = NA
          ),
          # Classical 背景
          geom_rect(
            data = data.frame(sample = unique(plot_df$sample)),
            aes(ymin = -Inf, ymax = Inf),
            xmin = 1.5, xmax = 2.5,
            inherit.aes = FALSE,
            fill = unname(dir_cols["Classical"]),
            alpha = shade_alpha, color = NA
          )
        )
      } } +
      geom_point(
        aes(fill = P_adj, size = GeneRatioNum),
        shape = 21, color = "white", stroke = 0.7, alpha = 0.95
      ) +
      facet_grid(. ~ sample, scales = "free_x", space = "free_x") +
      fill_scale + size_scale +
      scale_x_discrete(position = "top") +
      scale_y_discrete(labels = scales::label_wrap(35)) +
      base_theme +
      theme(
        strip.background = element_rect(fill = "#f5f5f5", color = NA),
        strip.text.x     = element_text(face = "bold", size = 12),
        panel.spacing.x  = grid::unit(panel_spacing_x, "lines"),
        axis.text.x      = if (direction_to_legend) element_blank() else element_text(angle = 0)
      ) +
      guide_setup
    
    # --- NEW: 把方向放入图例（使用 color 图例，避免与 fill 冲突）---
    if (isTRUE(direction_to_legend)) {
      dir_df <- data.frame(direction = factor(names(dir_cols), levels = dir_order))
      p <- p +
        geom_point(
          data = dir_df,
          inherit.aes = FALSE,
          aes(x = 1, y = 1, color = direction),
          size = 4, alpha = 0  # 不在图上显示，仅用于生成图例
        ) +
        scale_color_manual(
          name = "Direction",
          values = dir_cols,
          breaks = dir_order
        ) +
        guides(color = guide_legend(order = 3, override.aes = list(alpha = 1, size = 5)))
    }
    
  } else {
    # 非分面模式（保持原有后备逻辑）
    p <- ggplot(plot_df, aes(x = Group, y = Pathway)) +
      geom_point(
        aes(fill = P_adj, size = GeneRatioNum),
        shape = 21, color = "white", stroke = 0.7, alpha = 0.95
      ) +
      fill_scale + size_scale +
      scale_x_discrete(position = "top") +
      scale_y_discrete(labels = scales::label_wrap(35)) +
      base_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0)) +
      guide_setup
  }
  
  # ---- 4) Save outputs (optional) ----
  if (!is.null(save_prefix)) {
    readr::write_csv(enrich_all, paste0(save_prefix, "_Hallmark_Enrichment_ALL.csv"))
    readr::write_csv(plot_df,   paste0(save_prefix, "_Hallmark_Dotplot_Data.csv"))
    xl <- plot_df %>% arrange(sample, direction, P_adj)
    split_list <- split(xl, paste0(xl$sample, "_", xl$direction))
    writexl::write_xlsx(split_list, path = paste0(save_prefix, "_Hallmark_Enrichment_bySheet.xlsx"))
    ggsave(paste0(save_prefix, "_Hallmark_Dotplot_AllSamples.pdf"),
           p, width = fig_width, height = fig_height, useDingbats = FALSE)
  }
  
  list(
    plot = p,
    plot_data = plot_df,
    enrichment_tables = enrich_all
  )
}
# Pretty print Hallmark term -> publication style
pretty_hallmark <- function(x, width = 40) {
  out <- x
  # Remove prefix and underscores
  out <- gsub("^HALLMARK_", "", out)
  out <- gsub("_", " ", out)
  
  # Title case baseline
  out <- stringr::str_to_title(out)
  
  # Common tokens (restore correct case/symbols)
  fixes <- list(
    "Tgf Beta" = "TGF-β",
    "Tnfa" = "TNFα",
    "Nf Kb" = "NF-κB",
    "Mtorc1" = "mTORC1",
    "Kras" = "KRAS",
    "Myc" = "MYC",
    "Il " = "IL-",
    "Dn" = "Down",
    "Up" = "Up",
    "G2m" = "G2/M",
    "Uv" = "UV"
  )
  for (k in names(fixes)) out <- gsub(k, fixes[[k]], out, fixed = TRUE)
  
  # Specific long names → concise, publication style
  out <- gsub("Epithelial Mesenchymal Transition", "Epithelial–Mesenchymal Transition (EMT)", out)
  out <- gsub("Tnfa Signaling Via Nf-κb", "TNFα Signaling via NF-κB", out)
  out <- gsub("Tgf-β Signaling", "TGF-β Signaling", out)
  out <- gsub("Mtorc1 Signaling", "mTORC1 Signaling", out)
  out <- gsub("Oxidative Phosphorylation", "Oxidative Phosphorylation", out)
  out <- gsub("Apical Junction", "Apical Junction", out)
  out <- gsub("Myc Targets V1", "MYC Targets v1", out)
  out <- gsub("Myc Targets V2", "MYC Targets v2", out)
  out <- gsub("Kras Signaling Up", "KRAS Signaling (Up)", out)
  out <- gsub("Kras Signaling Down", "KRAS Signaling (Down)", out)
  
  # Final wrap
  out <- stringr::str_wrap(out, width = width)
  out
}


# ------------------------------------------------------------
# Plot violins of signature scores grouped by macro classes
# (修正版 - 删除了重复的 fig_height, 确保所有新参数都在)
# ------------------------------------------------------------
vln_signatures_by_macro <- function(
    obj_list,
    signatures     = c("Score_scBasal","Score_scClassical"),
    group_col      = "Tumor_L2_macro_std",
    include_groups = c("Classical","Basal-like"),
    pal_macro3     = c("Classical"="#1f78b4","Basal-like"="#e31a1c"),
    sample_order   = NULL,
    signature_order = c("Basal score","Classical score"),
    same_y         = TRUE,
    add_jitter     = FALSE,
    jitter_alpha   = 0.15,
    jitter_width   = 0.15,
    one_column     = FALSE,
    save_path      = NULL,
    fig_width      = 12,
    fig_height     = 6.5, # <-- 只保留一个
    
    # --- 新功能参数 ---
    flip_axes      = FALSE,
    signature_bg_colors = NULL
){
  # --- 0) Require namespaces ---
  pkgs <- c("dplyr","tidyr","ggplot2","purrr", "ggnewscale") # ggnewscale 是必需的
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Please install/load packages: ", paste(miss, collapse=", "))
  
  # --- 1) & 2) & 3) & 4) ... (函数的前半部分保持不变) ...
  # --- 1) Basic checks ---
  stopifnot(is.list(obj_list), length(obj_list) > 0)
  need_cols <- c(group_col, signatures)
  
  bad <- purrr::keep(names(obj_list), function(nm){
    !all(need_cols %in% colnames(obj_list[[nm]]@meta.data))
  })
  if (length(bad)) {
    stop("These Seurat objects are missing required columns (",
         paste(need_cols, collapse=", "), "): ",
         paste(bad, collapse=", "))
  }
  
  # --- 2) Bind meta from all samples ---
  df <- purrr::imap(obj_list, function(x, nm){
    x@meta.data |>
      dplyr::select(dplyr::all_of(c(group_col, signatures))) |>
      dplyr::mutate(sample = nm, .before = 1)
  }) |>
    dplyr::bind_rows()
  
  # --- 3) Keep two macro groups only; order factors ---
  df <- df |>
    dplyr::filter(.data[[group_col]] %in% include_groups) |>
    dplyr::mutate(
      macro  = factor(.data[[group_col]], levels = include_groups),
      sample = if (!is.null(sample_order))
        factor(sample, levels = sample_order) else factor(sample)
    )
  
  # --- 4) Long format for facets (signature as rows) ---
  df_long <- df |>
    tidyr::pivot_longer(cols = dplyr::all_of(signatures),
                        names_to = "signature", values_to = "score") |>
    dplyr::mutate(
      signature_label = dplyr::recode(
        signature,
        "Score_scBasal"      = "Basal score",
        "Score_scClassical"  = "Classical score",
        .default = signature
      )
    )
  
  if (!is.null(signature_order)) {
    df_long <- df_long |>
      dplyr::mutate(signature_label = factor(signature_label, levels = signature_order))
  }
  
  pal_use <- pal_macro3[include_groups]
  if (any(is.na(pal_use))) {
    stop("pal_macro3 must define colors for: ", paste(include_groups, collapse=", "))
  }
  
  # --- 5) 构建面板布局 ---
  if (isTRUE(one_column)) {
    sm_levels <- if (!is.null(sample_order)) sample_order else levels(df$sample)
    if (is.null(sm_levels)) sm_levels <- unique(df_long$sample)
    sig_levels <- if (!is.null(signature_order)) signature_order else unique(df_long$signature_label)
    panel_levels <- as.vector(outer(sm_levels, sig_levels, paste, sep = " • "))
    df_long <- df_long |>
      dplyr::mutate(
        panel = paste0(as.character(sample), " • ", as.character(signature_label)),
        panel = factor(panel, levels = panel_levels)
      )
  }
  
  # --- 6) Plot ---
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = macro, y = score))
  
  if (!is.null(signature_bg_colors)) {
    if (isTRUE(one_column)) {
      bg_df <- df_long |>
        dplyr::distinct(panel, signature_label) |>
        dplyr::mutate(bg_col = signature_bg_colors[as.character(signature_label)])
    } else {
      bg_df <- df_long |>
        dplyr::distinct(sample, signature_label) |>
        dplyr::mutate(bg_col = signature_bg_colors[as.character(signature_label)])
    }
    
    p <- p +
      ggplot2::geom_rect(
        data = bg_df,
        ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = bg_col),
        inherit.aes = FALSE, color = NA
      ) +
      ggplot2::scale_fill_identity() +
      ggnewscale::new_scale_fill()
  }
  
  p <- p +
    ggplot2::geom_violin(ggplot2::aes(fill = macro), width = 0.9, trim = TRUE, color = "white", linewidth = 0.3) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", alpha = 0.7,
                          color = "black", linewidth = 0.3,
                          position = ggplot2::position_dodge(width = 0.9)) +
    (if (isTRUE(add_jitter))
      ggplot2::geom_jitter(width = jitter_width, alpha = jitter_alpha, size = 0.6)
     else NULL) +
    ggplot2::scale_fill_manual(values = pal_use, name = NULL) +
    {
      if (isTRUE(one_column)) {
        ggplot2::facet_wrap(~ panel, ncol = 1, scales = if (same_y) "fixed" else "free_y")
      } else {
        ggplot2::facet_grid(rows = ggplot2::vars(signature_label),
                            cols  = ggplot2::vars(sample),
                            scales = if (same_y) "fixed" else "free_y")
      }
    } +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.border       = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.8),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.title.x       = ggplot2::element_blank(),
      axis.title.y       = ggplot2::element_blank(),
      legend.position    = "top",
      strip.text.x       = ggplot2::element_text(face = "bold"),
      strip.text.y       = ggplot2::element_text(face = "bold"),
      plot.margin        = ggplot2::margin(10,12,10,10)
    )
  
  if (isTRUE(flip_axes)) p <- p + ggplot2::coord_flip()
  
  # --- 7) Save (optional) ---
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = fig_width, height = fig_height, useDingbats = FALSE)
    message("Saved violin figure to: ", save_path)
  }
  
  return(p)
}

vln_signatures_by_macro <- function(
    obj_list,
    signatures      = c("Score_scBasal","Score_scClassical"),
    group_col       = "Tumor_L2_macro_std",
    include_groups  = c("Classical","Basal-like"),
    pal_macro3      = c("Classical"="#1f78b4","Basal-like"="#e31a1c"),
    sample_order    = NULL,
    signature_order = c("Basal score","Classical score"),
    same_y          = TRUE,
    add_jitter      = FALSE,
    jitter_alpha    = 0.15,
    jitter_width    = 0.15,
    one_column      = FALSE,
    save_path       = NULL,
    fig_width       = 12,
    fig_height      = 6.5,
    
    # layout helpers
    flip_axes           = FALSE,
    signature_bg_colors = NULL,
    
    # p-value options
    add_p_value     = FALSE,
    p_engine        = c("manual", "ggpubr"),
    p_method        = "wilcox.test",
    p_comparisons   = list(c("Classical","Basal-like")),
    p_label         = "p.signif",
    p_y_position    = NULL,
    p_step_increase = 0.1,
    p_size          = 4
){
  p_engine <- match.arg(p_engine)
  
  # --- Dependencies / checks ---
  pkgs <- c("dplyr","tidyr","ggplot2","purrr","ggnewscale","ggpubr","rstatix")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Please install/load packages: ", paste(miss, collapse=", "))
  
  # --- Data preparation ---
  stopifnot(is.list(obj_list), length(obj_list) > 0)
  need_cols <- c(group_col, signatures)
  bad <- purrr::keep(names(obj_list), function(nm){ 
    !all(need_cols %in% colnames(obj_list[[nm]]@meta.data)) 
  })
  if (length(bad)) { 
    stop("Objects missing columns: ", paste(bad, collapse=", ")) 
  }
  
  df <- purrr::imap(obj_list, ~ .x@meta.data |> 
                      dplyr::select(dplyr::all_of(c(group_col, signatures))) |> 
                      dplyr::mutate(sample = .y, .before = 1)) |> 
    dplyr::bind_rows()
  
  df <- df |> 
    dplyr::filter(.data[[group_col]] %in% include_groups) |> 
    dplyr::mutate(
      macro = factor(.data[[group_col]], levels = include_groups), 
      sample = if (!is.null(sample_order)) 
        factor(sample, levels = sample_order) 
      else 
        factor(sample)
    )
  
  df_long <- df |> 
    tidyr::pivot_longer(
      cols = dplyr::all_of(signatures), 
      names_to = "signature", 
      values_to = "score"
    ) |> 
    dplyr::mutate(
      signature_label = dplyr::recode(
        signature, 
        "Score_scBasal" = "Basal score", 
        "Score_scClassical" = "Classical score", 
        .default = signature
      )
    )
  
  if (!is.null(signature_order)) { 
    df_long <- df_long |> 
      dplyr::mutate(signature_label = factor(signature_label, levels = signature_order)) 
  }
  
  pal_use <- pal_macro3[include_groups]
  if (any(is.na(pal_use))) { 
    stop("pal_macro3 must define colors for: ", paste(include_groups, collapse=", ")) 
  }
  
  # One-column layout
  if (isTRUE(one_column)) {
    sm_levels <- if (!is.null(sample_order)) sample_order else levels(df$sample)
    if (is.null(sm_levels)) sm_levels <- unique(df_long$sample)
    sig_levels <- if (!is.null(signature_order)) signature_order else unique(df_long$signature_label)
    panel_levels <- as.vector(outer(sm_levels, sig_levels, paste, sep = " • "))
    df_long <- df_long |> 
      dplyr::mutate(
        panel = paste0(as.character(sample), " • ", as.character(signature_label)), 
        panel = factor(panel, levels = panel_levels)
      )
  }
  
  # --- Base plot ---
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = macro, y = score))
  
  # Background colors
  if (!is.null(signature_bg_colors)) {
    if (isTRUE(one_column)) { 
      bg_df <- df_long |> 
        dplyr::distinct(panel, signature_label) |> 
        dplyr::mutate(bg_col = signature_bg_colors[as.character(signature_label)]) 
    } else { 
      bg_df <- df_long |> 
        dplyr::distinct(sample, signature_label) |> 
        dplyr::mutate(bg_col = signature_bg_colors[as.character(signature_label)]) 
    }
    p <- p + 
      ggplot2::geom_rect(
        data = bg_df, 
        ggplot2::aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = bg_col), 
        inherit.aes = FALSE, color = NA
      ) + 
      ggplot2::scale_fill_identity() + 
      ggnewscale::new_scale_fill()
  }
  
  # Violin + boxplot + jitter
  p <- p +
    ggplot2::geom_violin(
      ggplot2::aes(fill = macro), 
      width = 0.9, trim = TRUE, 
      color = "white", linewidth = 0.3
    ) +
    ggplot2::geom_boxplot(
      width = 0.18, outlier.shape = NA, 
      fill = "white", alpha = 0.7, 
      color = "black", linewidth = 0.3
    ) +
    {if (isTRUE(add_jitter)) 
      ggplot2::geom_jitter(width = jitter_width, alpha = jitter_alpha, size = 0.6) 
    } +
    ggplot2::scale_fill_manual(values = pal_use, name = NULL)
  
  # Faceting
  if (isTRUE(one_column)) { 
    p <- p + ggplot2::facet_wrap(
      ~ panel, ncol = 1, 
      scales = if (same_y) "fixed" else "free_y"
    )
  } else { 
    p <- p + ggplot2::facet_grid(
      rows = ggplot2::vars(signature_label), 
      cols = ggplot2::vars(sample), 
      scales = if (same_y) "fixed" else "free_y"
    )
  }
  
  # --- P-value layer (FIXED) ---
  if (isTRUE(add_p_value)) {
    
    # Identify valid facets (must have exactly 2 groups)
    if (isTRUE(one_column)) {
      valid_facets <- df_long |> 
        dplyr::group_by(panel) |> 
        dplyr::summarise(n_groups = dplyr::n_distinct(macro), .groups = "drop") |> 
        dplyr::filter(n_groups == 2L)
      
      # Calculate max y for positioning
      max_pos <- df_long |> 
        dplyr::semi_join(valid_facets, by = "panel") |> 
        dplyr::group_by(panel) |> 
        dplyr::summarise(y_max = max(score, na.rm = TRUE), .groups = "drop")
      
      # Compute p-values
      .p_stars <- function(p) { 
        as.character(cut(
          p, 
          breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf), 
          labels = c("****","***","**","*","ns"), 
          right = TRUE
        )) 
      }
      
      cmp <- ggpubr::compare_means(
        score ~ macro, 
        group.by = "panel", 
        data = df_long, 
        method = p_method
      ) |> 
        dplyr::filter(
          (group1 == p_comparisons[[1]][1] & group2 == p_comparisons[[1]][2]) | 
            (group1 == p_comparisons[[1]][2] & group2 == p_comparisons[[1]][1])
        ) |> 
        dplyr::semi_join(valid_facets, by = "panel") |> 
        dplyr::left_join(max_pos, by = "panel") |> 
        dplyr::mutate(
          p.signif = .p_stars(p),
          p.format = if (requireNamespace("rstatix", quietly = TRUE)) 
            rstatix::p_format(p) 
          else 
            formatC(p, format = "e", digits = 2)
        )
      
      # Set y position
      if (is.null(p_y_position)) { 
        cmp <- cmp |> 
          dplyr::mutate(y.position = y_max * (1 + max(0.05, p_step_increase))) 
      } else { 
        cmp <- cmp |> 
          dplyr::mutate(y.position = p_y_position) 
      }
      
      # KEY FIX: Add required columns with proper values
      cmp <- cmp |> 
        dplyr::mutate(
          xmin = 1,  # Position for first group
          xmax = 2,  # Position for second group
          # Remove problematic columns or set to NULL
          step.increase = 0,
          bracket.nudge.y = 0,
          bracket.shorten = 0
        )
      
    } else {
      # For grid layout
      valid_facets <- df_long |> 
        dplyr::group_by(sample, signature_label) |> 
        dplyr::summarise(n_groups = dplyr::n_distinct(macro), .groups = "drop") |> 
        dplyr::filter(n_groups == 2L)
      
      max_pos <- df_long |> 
        dplyr::semi_join(valid_facets, by = c("sample","signature_label")) |> 
        dplyr::group_by(sample, signature_label) |> 
        dplyr::summarise(y_max = max(score, na.rm = TRUE), .groups = "drop")
      
      .p_stars <- function(p) { 
        as.character(cut(
          p, 
          breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf), 
          labels = c("****","***","**","*","ns"), 
          right = TRUE
        )) 
      }
      
      cmp <- ggpubr::compare_means(
        score ~ macro, 
        group.by = c("sample","signature_label"), 
        data = df_long, 
        method = p_method
      ) |> 
        dplyr::filter(
          (group1 == p_comparisons[[1]][1] & group2 == p_comparisons[[1]][2]) | 
            (group1 == p_comparisons[[1]][2] & group2 == p_comparisons[[1]][1])
        ) |> 
        dplyr::semi_join(valid_facets, by = c("sample","signature_label")) |> 
        dplyr::left_join(max_pos, by = c("sample","signature_label")) |> 
        dplyr::mutate(
          p.signif = .p_stars(p),
          p.format = if (requireNamespace("rstatix", quietly = TRUE)) 
            rstatix::p_format(p) 
          else 
            formatC(p, format = "e", digits = 2)
        )
      
      if (is.null(p_y_position)) { 
        cmp <- cmp |> 
          dplyr::mutate(y.position = y_max * (1 + max(0.05, p_step_increase))) 
      } else { 
        cmp <- cmp |> 
          dplyr::mutate(y.position = p_y_position) 
      }
      
      # KEY FIX: Add required columns
      cmp <- cmp |> 
        dplyr::mutate(
          xmin = 1,
          xmax = 2,
          step.increase = 0,
          bracket.nudge.y = 0,
          bracket.shorten = 0
        )
    }
    
    # Add p-values to plot
    p <- p + 
      ggpubr::stat_pvalue_manual(
        cmp, 
        label = if (identical(p_label, "p.signif")) "p.signif" else "p.format", 
        y.position = "y.position",
        xmin = "xmin",
        xmax = "xmax",
        tip.length = 0.01, 
        size = p_size,
        inherit.aes = FALSE
      )
  }
  
  # --- Theme and axis labels ---
  p <- p +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.border       = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.8),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.title.x       = ggplot2::element_blank(),
      axis.title.y       = ggplot2::element_text(face = "bold", size = 12),
      legend.position    = "top",
      strip.text.x       = ggplot2::element_text(face = "bold"),
      strip.text.y       = ggplot2::element_text(face = "bold", angle = 0),
      plot.margin        = ggplot2::margin(10,12,10,10)
    ) +
    ggplot2::labs(y = "Signature Score (Z-score)")
  
  # Axis flip + title swap
  if (isTRUE(flip_axes)) {
    p <- p + 
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Signature Score (Z-score)", y = NULL)
  }
  
  # --- Save (optional) ---
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      save_path, p, 
      width = fig_width, 
      height = fig_height, 
      useDingbats = FALSE
    )
    message("Saved violin figure to: ", save_path)
  }
  
  return(p)
}


# Drop-in replacement: paste this above where pretty_term() is used
pretty_term <- function(
    x,
    wrap_width = 40,   # 新参数名
    width = NULL,      # 兼容老调用: pretty_term(..., width=)
    prettify = c("auto","all","hallmark","go","kegg","reactome","msig","none")
){
  if (!is.null(width)) wrap_width <- width
  prettify <- match.arg(prettify)
  
  if (prettify == "none") return(stringr::str_wrap(x, width = wrap_width))
  
  # 统一向量化处理
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    
    s0 <- s
    
    # ---------- 1) 去前缀 ----------
    # 常见 MSigDB / 集合前缀
    s0 <- gsub("^(HALLMARK|REACTOME|KEGG|BIOCARTA|PID|WP|NABA|HP|MP|GOBP|GOCC|GOMF|GO_BP|GO_CC|GO_MF)_", "", s0)
    # 去掉 GO ID
    s0 <- gsub("\\(GO:\\d{7}\\)", "", s0)         # 例如 "... (GO:0008150)"
    s0 <- gsub("^GO:\\d{7}[ _-]*", "", s0)        # 例如 "GO:0008150_BIOLOGICAL_PROCESS"
    
    # ---------- 2) 基础清洗 ----------
    s0 <- gsub("_+", " ", s0)                     # 下划线 → 空格
    s0 <- gsub("\\s+", " ", s0)                   # 多空格折叠
    s0 <- trimws(s0)
    
    # ---------- 3) 标题化（先统一成 Title Case，再修复缩写） ----------
    s0 <- stringr::str_to_title(s0)
    
    # 保留常见缩写/通用基因通路名的大小写
    fix_acronyms <- function(z) {
      # 一般信号分子/通路缩写
      repl <- c(
        "Tgf"="TGF", "Tnf"="TNF", "Nf"="NF", "Kappa"="Kappa", "Il"="IL", "Ifn"="IFN",
        "Vegf"="VEGF", "Egfr"="EGFR", "Mapk"="MAPK", "Erk"="ERK", "Jnk"="JNK", "P38"="p38",
        "Akt"="AKT", "Pi3k"="PI3K", "MtOr"="mTOR", "Wnt"="WNT", "Tlr"="TLR",
        "Dna"="DNA", "Rna"="RNA", "Atp"="ATP", "Ampk"="AMPK", "Tgf-β"="TGF-β"
      )
      for (k in names(repl)) {
        z <- gsub(paste0("\\b", k, "\\b"), repl[[k]], z)
      }
      z
    }
    s0 <- fix_acronyms(s0)
    
    # ---------- 4) 连接词 & 语气词优化 ----------
    # AND/OR/VIA/OF 等的风格化
    s0 <- gsub("\\bAnd\\b", "&", s0)
    s0 <- gsub("\\bOr\\b", "/", s0)
    s0 <- gsub("\\bVia\\b", "via", s0)
    s0 <- gsub("\\bOf\\b", "of", s0)
    s0 <- gsub("\\bTo\\b", "to", s0)
    s0 <- gsub("\\bIn\\b", "in", s0)
    s0 <- gsub("\\bFor\\b", "for", s0)
    
    # ---------- 5) 希腊字母/特例 ----------
    # 常见 Beta/Alpha/Gamma/Kappa/Delta → 希腊字母
    greek_map <- list(" Alpha"="-α", " Beta"="-β", " Gamma"="-γ", " Kappa"="-κ", " Delta"="-δ")
    for (k in names(greek_map)) {
      s0 <- gsub(k, greek_map[[k]], s0, fixed = TRUE)
    }
    # NF-κB 特例（Title Case后会成"Nf Kappa B"）
    s0 <- gsub("\\bNf[- ]?Kappa[- ]?B\\b", "NF-κB", s0, ignore.case = TRUE)
    # TGF-β 特例
    s0 <- gsub("\\bTgf[- ]?Beta\\b", "TGF-β", s0, ignore.case = TRUE)
    # IL-1β / IL-6 等
    s0 <- gsub("\\bIl[- ]?(\\d+) Beta\\b", "IL-\\1β", s0, ignore.case = TRUE)
    s0 <- gsub("\\bIl[- ]?(\\d+) Alpha\\b", "IL-\\1α", s0, ignore.case = TRUE)
    
    # ---------- 6) 罗马数字（Type II / III / IV 等） ----------
    s0 <- gsub("\\bType Ii\\b", "Type II", s0)
    s0 <- gsub("\\bType Iii\\b", "Type III", s0)
    s0 <- gsub("\\bType Iv\\b", "Type IV", s0)
    
    # ---------- 7) “signaling pathway” 等冗余收敛 ----------
    s0 <- gsub("\\bSignalling\\b", "Signaling", s0)  # 英式→美式统一
    s0 <- gsub("\\bSignaling Pathway\\b", "signaling", s0)
    s0 <- gsub("\\bPathway\\b$", "", s0)            # 末尾单独 Pathway 去掉
    s0 <- gsub("\\s+", " ", s0); s0 <- trimws(s0)
    
    # ---------- 8) 包装换行 ----------
    stringr::str_wrap(s0, width = wrap_width)
  }, FUN.VALUE = character(1))
}


# =====================================================================
# Generic enrichment dot-plot for arbitrary "direction" labels
# - Works with GMT (symbols) or MSigDB (msigdbr, `collection=` API)
# - Robust to either `p_val_adj` or `p.adjust`
# - No ID conversion needed (uses clusterProfiler::enricher with SYMBOLs)
# - Optional per-sample universe from obj_list (rownames assumed SYMBOLs)
# - Clean background shading per direction via annotate(), no ggplot warnings
# - Optional legend for directions (keeps x-axis clean)
# =====================================================================

enrichment_dotplot_generic <- function(
    de_all,
    direction_col        = "up_in",              # column with direction labels (e.g., "Basal-like"/"Classical", "iCAF"/"myCAF")
    dir_levels           = NULL,                 # order of directions on x-axis; default: unique order in data
    dir_cols             = NULL,                 # optional named colors for directions; required if shade_direction or direction_to_legend = TRUE
    gmt_path             = NULL,                 # path to a .gmt file (preferred when provided)
    msig_source          = NULL,                 # list(species="Homo sapiens", collection="H"). `category` also accepted for backward-compat.
    q_cut                = 0.05,                 # adj. p-value threshold on DE table
    fc_cut               = 0.25,                 # abs(log2FC) threshold on DE table
    top_terms_per_group  = 10,                   # top terms per (sample x direction) by adjusted p
    wrap_width           = 40,                   # text wrap width for pathway labels
    x_order_by_sample    = NULL,                 # order of samples on the facet columns
    obj_list             = NULL,                 # optional named list of Seurat objects; if present, use per-sample rownames as universe
    save_prefix          = NULL,                 # if not NULL, save CSV/XLSX/PDF using this prefix
    fig_width            = 12,
    fig_height           = 9,
    facet_by_sample      = TRUE,                 # TRUE: facets by sample; x-axis shows directions only
    shade_direction      = TRUE,                 # add faint background color bands for each direction
    shade_alpha          = 0.06,                 # transparency for background shading
    panel_spacing_x      = 1.2,                  # horizontal spacing between sample facets
    direction_to_legend  = TRUE                  # move direction labels/colors into a legend (hides x text)
){
  # ----- helpers -----
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  .need <- function(pkgs){
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss)) stop("Please install/load packages: ", paste(miss, collapse = ", "))
  }
  
  .need(c("dplyr","tidyr","stringr","ggplot2","clusterProfiler"))
  
  # ----- 0) Basic checks & robust column mapping -----
  if (!all(c("sample","gene", direction_col) %in% colnames(de_all))) {
    stop("`de_all` must contain columns: sample, gene, and ", direction_col)
  }
  p_col <- if ("p_val_adj" %in% names(de_all)) "p_val_adj"
  else if ("p.adjust" %in% names(de_all)) "p.adjust"
  else stop("`de_all` must contain either 'p_val_adj' or 'p.adjust'")
  if (!("avg_log2FC" %in% names(de_all))) stop("`de_all` must contain 'avg_log2FC'")
  
  de_all <- dplyr::mutate(
    de_all,
    PAdj      = .data[[p_col]],
    Direction = .data[[direction_col]]
  )
  
  # direction order
  if (is.null(dir_levels)) {
    dir_levels <- unique(as.character(de_all[[direction_col]]))
  }
  de_all$Direction <- factor(de_all$Direction, levels = dir_levels)
  
  # sample order
  if (!is.null(x_order_by_sample)) {
    de_all$sample <- factor(de_all$sample, levels = x_order_by_sample)
  } else {
    de_all$sample <- factor(de_all$sample)
  }
  
  # ----- 1) TERM2GENE source (precedence: GMT > MSigDB) -----
  if (!is.null(gmt_path)) {
    if (!file.exists(gmt_path)) stop("GMT not found: ", gmt_path)
    TERM2GENE <- clusterProfiler::read.gmt(gmt_path) %>%
      dplyr::rename(term = 1, gene = 2)
  } else if (!is.null(msig_source)) {
    .need("msigdbr")
    # backward-compat: allow `category`; prefer `collection`
    collection <- msig_source$collection %||% msig_source$category %||% "H"
    species    <- msig_source$species    %||% "Homo sapiens"
    TERM2GENE <- msigdbr::msigdbr(species = species, collection = collection) %>%
      dplyr::select(gs_name, gene_symbol) %>%
      dplyr::rename(term = gs_name, gene = gene_symbol)
  } else {
    stop("Provide either `gmt_path` or `msig_source`.")
  }
  
  # fallback pretty term: Hallmark nicer; otherwise wrap
  
  # ----- 2) Build per-sample, per-direction gene lists and run enricher() -----
  samples <- levels(de_all$sample)
  enrich_list <- list()
  
  for (sm in samples) {
    de_sm <- dplyr::filter(de_all, .data$sample == sm)
    
    # per-sample universe (optional)
    uni <- NULL
    if (!is.null(obj_list) && !is.null(obj_list[[as.character(sm)]])) {
      uni <- tryCatch(rownames(obj_list[[as.character(sm)]]), error = function(e) NULL)
    }
    
    for (lev in dir_levels) {
      de_sd <- dplyr::filter(de_sm, .data$Direction == lev)
      genes_up <- de_sd %>%
        dplyr::filter(.data$PAdj < q_cut, abs(.data$avg_log2FC) >= fc_cut) %>%
        dplyr::pull(gene) %>%
        unique()
      
      if (length(genes_up) < 5) next
      
      ek <- tryCatch(
        clusterProfiler::enricher(
          gene         = genes_up,
          TERM2GENE    = TERM2GENE,
          universe     = uni,
          pvalueCutoff = 0.1,
          qvalueCutoff = 0.2
        ),
        error = function(e) NULL
      )
      if (is.null(ek)) next
      
      df <- as.data.frame(ek)
      if (!nrow(df)) next
      
      df2 <- df %>%
        dplyr::mutate(sample = sm, direction = lev) %>%
        dplyr::select(sample, direction, ID, Description, GeneRatio, BgRatio,
                      pvalue, p.adjust, qvalue, Count, dplyr::everything())
      
      enrich_list[[paste(sm, lev, sep = "||")]] <- df2
    }
  }
  
  enrich_all <- dplyr::bind_rows(enrich_list)
  if (is.null(enrich_all) || !nrow(enrich_all)) {
    stop("No enrichment results after filtering. Consider relaxing q_cut/fc_cut or check DE input.")
  }
  
  # ----- 3) Prepare plot data -----
  gene_ratio_num <- function(x) sapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]))
  
  plot_df <- enrich_all %>%
    dplyr::group_by(sample, direction) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_terms_per_group) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Group        = paste(sample, direction, sep = " | "),
      P_adj        = p.adjust,
      GeneRatioNum = gene_ratio_num(GeneRatio),
      Pathway      = pretty_term(Description, width = wrap_width),
      direction    = factor(direction, levels = dir_levels),
      sample       = factor(sample, levels = samples)
    )
  
  # order y by best FDR
  path_order <- plot_df %>%
    dplyr::group_by(Pathway) %>%
    dplyr::summarise(min_padj = min(P_adj, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(min_padj) %>%
    dplyr::pull(Pathway)
  plot_df$Pathway <- factor(plot_df$Pathway, levels = rev(path_order))
  
  # ----- 4) Scales & theme -----
  fill_scale <- ggplot2::scale_fill_gradientn(
    colors = c("#B30000", "#F46D43", "#FEE08B", "#ABD9E9", "#4575B4"),
    values = scales::rescale(c(0, 0.01, 0.02, 0.05)),
    limits = c(0, 0.05), oob = scales::squish,
    name = "FDR (adj. p)",
    breaks = c(0.001, 0.01, 0.02, 0.05),
    labels = c("0.001", "0.01", "0.02", "0.05")
  )
  size_scale <- ggplot2::scale_size_continuous(
    range = c(3.8, 11),
    breaks = scales::pretty_breaks(n = 3),
    name   = "Gene ratio"
  )
  base_theme <- ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.border       = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.8),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      axis.title.x       = ggplot2::element_blank(),
      axis.title.y       = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      legend.title       = ggplot2::element_text(size = 12),
      legend.text        = ggplot2::element_text(size = 11),
      plot.margin        = ggplot2::margin(t = 10, r = 12, b = 10, l = 10)
    )
  
  # ----- 5) Build plot -----
  if (isTRUE(facet_by_sample)) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = direction, y = Pathway))
    
    # background shading per direction (no warnings)
    if (isTRUE(shade_direction)) {
      if (is.null(dir_cols)) stop("Please provide `dir_cols` when `shade_direction = TRUE`.")
      # shade each discrete x band according to dir_cols (replicated to facets)
      for (i in seq_along(dir_levels)) {
        col_i <- unname(dir_cols[ dir_levels[i] ])
        if (is.na(col_i)) col_i <- "#00000010"
        p <- p + ggplot2::annotate(
          "rect",
          xmin = i - 0.5, xmax = i + 0.5, ymin = -Inf, ymax = Inf,
          fill = col_i, alpha = shade_alpha, color = NA
        )
      }
    }
    
    p <- p +
      ggplot2::geom_point(
        ggplot2::aes(fill = P_adj, size = GeneRatioNum),
        shape = 21, color = "white", stroke = 0.7, alpha = 0.95
      ) +
      ggplot2::facet_grid(. ~ sample, scales = "fixed", space = "fixed")  +
      fill_scale + size_scale +
      ggplot2::scale_x_discrete(drop = FALSE, limits = dir_levels, position = "top") +
      ggplot2::scale_y_discrete(labels = scales::label_wrap(35)) +
      base_theme +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "#f5f5f5", color = NA),
        strip.text.x     = ggplot2::element_text(face = "bold", size = 12),
        panel.spacing.x  = grid::unit(panel_spacing_x, "lines"),
        axis.text.x      = if (direction_to_legend) ggplot2::element_blank() else ggplot2::element_text(angle = 0)
      ) +
      ggplot2::guides(
        fill = ggplot2::guide_colourbar(order = 1, reverse = TRUE, ticks = FALSE,
                                        barheight = grid::unit(12, "lines"),
                                        barwidth  = grid::unit(0.7,  "lines"),
                                        frame.colour = "black", frame.linewidth = 0.8),
        size = ggplot2::guide_legend(order = 2,
                                     override.aes = list(shape = 21, fill = NA, color = "black", stroke = 1.2))
      )
    
    # optional direction legend (keep x-labels clean)
    if (isTRUE(direction_to_legend)) {
      if (is.null(dir_cols)) stop("Please provide `dir_cols` when `direction_to_legend = TRUE`.")
      dir_df <- data.frame(direction = factor(dir_levels, levels = dir_levels), x = 1, y = 1)
      p <- p +
        ggplot2::geom_point(
          data = dir_df, inherit.aes = FALSE,
          ggplot2::aes(x = x, y = y, color = direction),
          size = 4, alpha = 0
        ) +
        ggplot2::scale_color_manual(
          name   = "Direction",
          values = dir_cols,
          breaks = dir_levels
        ) +
        ggplot2::guides(color = ggplot2::guide_legend(order = 3, override.aes = list(alpha = 1, size = 5)))
    }
    
  } else {
    # non-faceted layout (x = "sample | direction")
    plot_df <- dplyr::mutate(plot_df, Group = factor(Group, levels = unique(Group)))
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Group, y = Pathway)) +
      ggplot2::geom_point(
        ggplot2::aes(fill = P_adj, size = GeneRatioNum),
        shape = 21, color = "white", stroke = 0.7, alpha = 0.95
      ) +
      fill_scale + size_scale +
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::scale_y_discrete(labels = scales::label_wrap(35)) +
      base_theme +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0)) +
      ggplot2::guides(
        fill = ggplot2::guide_colourbar(order = 1, reverse = TRUE, ticks = FALSE,
                                        barheight = grid::unit(12, "lines"),
                                        barwidth  = grid::unit(0.7,  "lines"),
                                        frame.colour = "black", frame.linewidth = 0.8),
        size = ggplot2::guide_legend(order = 2,
                                     override.aes = list(shape = 21, fill = NA, color = "black", stroke = 1.2))
      )
  }
  
  # ----- 6) Save outputs (optional) -----
  if (!is.null(save_prefix)) {
    .need(c("readr","writexl"))
    readr::write_csv(enrich_all, paste0(save_prefix, "_Enrichment_ALL.csv"))
    readr::write_csv(plot_df,   paste0(save_prefix, "_Dotplot_Data.csv"))
    
    xl <- plot_df %>% dplyr::arrange(sample, direction, P_adj)
    split_list <- split(x = xl, f = paste0(xl$sample, "_", xl$direction))
    writexl::write_xlsx(split_list, path = paste0(save_prefix, "_Enrichment_bySheet.xlsx"))
    
    ggplot2::ggsave(paste0(save_prefix, "_Dotplot.pdf"),
                    p, width = fig_width, height = fig_height, useDingbats = FALSE)
  }
  
  list(
    plot = p,
    plot_data = plot_df,
    enrichment_tables = enrich_all
  )
}

## =========================
voom_anova <- function(counts_mat,
                       colData_df,
                       group_col = "treatment",
                       fdr_thresh = 0.05,
                       trend = TRUE,
                       robust = TRUE,
                       use_quality_weights = FALSE,
                       # expression filters
                       filter_min_count = 10,         # counts
                       filter_min_total_count = 15,   # counts
                       filter_cpm = 1,                # CPM threshold
                       filter_min_samples = NULL,     # e.g. 6 (全样本数中至少6个样本 CPM>=1)
                       # effect-size / variability filters
                       min_delta_lfc = NA_real_,      # max(group mean) - min(group mean) on log2 scale
                       min_sd = NA_real_,             # SD across group means
                       # fallback caps
                       min_sig_genes = 200,
                       max_sig_genes = Inf) {
  
  counts_mat <- as.matrix(counts_mat); storage.mode(counts_mat) <- "numeric"
  stopifnot(ncol(counts_mat) == nrow(colData_df))
  group <- factor(colData_df[[group_col]])
  
  # ---- filtering ----
  dge <- edgeR::DGEList(counts = counts_mat)
  keep <- edgeR::filterByExpr(dge, group = group,
                              min.count = filter_min_count,
                              min.total.count = filter_min_total_count)
  if (!is.null(filter_min_samples)) {
    cpm_mat <- edgeR::cpm(dge)
    keep_cpm <- rowSums(cpm_mat >= filter_cpm) >= filter_min_samples
    keep <- keep & keep_cpm
  }
  dge <- dge[keep,, keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  message("After expression filtering: ", nrow(dge), " genes")
  
  # ---- voom + lmFit ----
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  v <- if (isTRUE(use_quality_weights)) {
    limma::voomWithQualityWeights(dge, design, plot = FALSE)
  } else {
    limma::voom(dge, design, plot = FALSE)
  }
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit, trend = trend, robust = robust)
  
  # ---- global F (ANOVA-like) ----
  tabF <- limma::topTable(fit, coef = NULL, number = Inf, sort.by = "F", adjust.method = "BH")
  tabF$FDR <- tabF$adj.P.Val
  sig_flag <- tabF$FDR < fdr_thresh
  message("FDR <", fdr_thresh, " : ", sum(sig_flag), " genes")
  
  # ---- effect-size filters on group means ----
  prof_mu <- fit$coefficients[, levels(group), drop = FALSE]
  if (!all(rownames(prof_mu) == rownames(tabF))) {
    prof_mu <- prof_mu[rownames(tabF), , drop = FALSE]
  }
  
  if (!is.na(min_delta_lfc)) {
    delta <- apply(prof_mu, 1, function(x) diff(range(x)))
    sig_flag <- sig_flag & (delta >= min_delta_lfc)
    message("Δ(max-min) ≥ ", min_delta_lfc, " : ", sum(sig_flag), " genes")
  }
  if (!is.na(min_sd)) {
    sdg <- apply(prof_mu, 1, sd)
    sig_flag <- sig_flag & (sdg >= min_sd)
    message("SD ≥ ", min_sd, " : ", sum(sig_flag), " genes")
  }
  
  sig_genes <- rownames(tabF)[sig_flag]
  
  # fallback窗：至少/至多
  if (length(sig_genes) < min_sig_genes) {
    sig_genes <- rownames(tabF)[seq_len(min(min_sig_genes, nrow(tabF)))]
    message("Fallback min_sig_genes → kept ", length(sig_genes))
  }
  if (is.finite(max_sig_genes) && length(sig_genes) > max_sig_genes) {
    sig_genes <- sig_genes[seq_len(max_sig_genes)]
    message("Capped to max_sig_genes → kept ", length(sig_genes))
  }
  
  list(
    group_levels = levels(group),
    group_factor = group,
    dge   = dge,
    voom  = v,
    fit   = fit,
    tabF  = tabF,
    sig_genes = sig_genes
  )
}


######################################################################
######################################################################
###################################
###################################
#Drug

score_kras_resistance <- function(expr, signature_df,
                                  gene_col = "Gene", weight_col = "Weight",
                                  log1p = TRUE, zscore_genes = TRUE,
                                  return_01 = TRUE) {
  # expr: matrix/data.frame, rows=genes, cols=samples
  
  expr <- as.matrix(expr)
  sig  <- signature_df %>% dplyr::select(all_of(gene_col), all_of(weight_col)) %>% distinct()
  
  common <- intersect(rownames(expr), sig[[gene_col]])
  if (length(common) < 10) stop("Too few signature genes found in expr (check gene symbols).")
  
  w <- sig[[weight_col]][match(common, sig[[gene_col]])]
  names(w) <- common
  
  X <- expr[common, , drop = FALSE]
  if (log1p) X <- log2(X + 1)
  
  if (zscore_genes) {
    # Z-score each gene across samples to reduce scale differences
    X <- t(scale(t(X)))
    X[is.na(X)] <- 0
  }
  
  raw_score <- as.numeric(crossprod(w, X))
  names(raw_score) <- colnames(X)
  
  out <- data.frame(sample_id = names(raw_score), KRASi_resistance_score = raw_score)
  
  if (return_01) {
    rng <- range(out$KRASi_resistance_score, na.rm = TRUE)
    if (diff(rng) == 0) {
      out$KRASi_resistance_score01 <- 0.5
    } else {
      out$KRASi_resistance_score01 <- (out$KRASi_resistance_score - rng[1]) / diff(rng)
    }
  }
  out
}





