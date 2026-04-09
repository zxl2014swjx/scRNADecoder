################################## START ##################################
#' data_preprocess.R
# 输入: 单细胞转录组表达矩阵，行是基因，列是细胞
# 输出: 1.Data_preprocess文件夹及其所有分析结果
#' 单细胞表达谱矩阵预处理函数
#' 对单细胞表达谱矩阵进行质量控制、过滤和标准化
#' @param input 表达谱矩阵，行是基因，列是细胞
#' @param gene_min_cells 基因过滤条件：一个基因至少要在多少个细胞中表达（或百分比），默认5%（0.05）
#' @param cell_min_genes 细胞过滤条件：一个细胞至少要表达多少个基因，默认200
#' @param cell_max_genes 细胞过滤条件：一个细胞最多表达多少个基因（用于去除doublet），默认NULL（不过滤）
#' @param cell_max_mito 细胞过滤条件：线粒体基因表达百分比最大值，默认20%（0.2）
#' @param normalization 标准化方法："log"（log1p），"zscore"，"cpm"，"scran" 或 "none"
#' @param scale_factor 缩放因子，用于CPM标准化，默认1e6
#' @param do_scaling 是否在log转换后进行z-score缩放，默认TRUE
#' @param remove_mito 是否移除线粒体基因，默认FALSE
#' @param mito_pattern 线粒体基因匹配模式，默认"^MT-|^mt-"
#' @param do_hvg 是否选择高变基因，默认FALSE
#' @param n_hvg 选择的高变基因数量，默认2000
#' @param verbose 是否显示详细信息，默认TRUE
data_preprocess <- function(input,gene_min_cells = 0.05,cell_min_genes = 200,cell_max_genes = NULL,cell_max_mito = 0.2,normalization = "log",scale_factor = 1e6,do_scaling = TRUE,remove_mito = FALSE,mito_pattern = "^MT-|^mt-",do_hvg = FALSE,n_hvg = 2000,verbose = TRUE,output_dir="1.Data_preprocess") {
  expr_matrix<-read.table(input,header=T,sep="\t",row.names=1)
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
  if (verbose) {
    cat(sprintf("输入矩阵维度: %d 基因 × %d 细胞\n", nrow(expr_matrix), ncol(expr_matrix)))
    cat(sprintf("总表达值: %.0f\n", sum(expr_matrix)))
  }
  original_matrix <- expr_matrix
  if (!is.matrix(expr_matrix)) {
    expr_matrix <- as.matrix(expr_matrix)
    if (verbose) cat("  将输入数据转换为矩阵格式\n")
  }
  rows_zero <- rowSums(expr_matrix) == 0
  cols_zero <- colSums(expr_matrix) == 0
  if (any(rows_zero) || any(cols_zero)) {
    if (verbose) {
      cat(sprintf("  移除全零行: %d 个\n", sum(rows_zero)))
      cat(sprintf("  移除全零列: %d 个\n", sum(cols_zero)))
    }
    expr_matrix <- expr_matrix[!rows_zero, !cols_zero]
  }
  if (verbose) cat("\n识别线粒体基因...\n")
  mito_genes <- grep(mito_pattern, rownames(expr_matrix), value = TRUE)
  if (length(mito_genes) > 0) {
    if (verbose) {
      cat(sprintf("  找到 %d 个线粒体基因: %s\n", length(mito_genes), paste(head(mito_genes, 3), collapse = ", ")))
      if (length(mito_genes) > 3) cat(" ...等\n")
    }
  } else {
    if (verbose) cat("  未找到线粒体基因\n")
  }
  if (verbose) cat("\n计算质量控制指标...\n")
  cell_stats <- data.frame(
    cell_id = colnames(expr_matrix),
    n_genes = colSums(expr_matrix > 0),  # 表达基因数
    total_counts = colSums(expr_matrix),  # 总表达量
    mito_percent = 0)
  if (length(mito_genes) > 0) {
    mito_counts <- colSums(expr_matrix[mito_genes, , drop = FALSE])
    cell_stats$mito_percent <- mito_counts / cell_stats$total_counts
    cell_stats$mito_percent[is.na(cell_stats$mito_percent)] <- 0
  }
  gene_stats <- data.frame(
    gene_id = rownames(expr_matrix),
    n_cells = rowSums(expr_matrix > 0),  # 表达细胞数
    total_counts = rowSums(expr_matrix),  # 总表达量
    mean_expression = rowMeans(expr_matrix),
    is_mito = rownames(expr_matrix) %in% mito_genes
  )
  if (verbose) {
    cat("  细胞统计:\n")
    cat(sprintf("    平均每个细胞表达基因数: %.1f (范围: %d-%d)\n",mean(cell_stats$n_genes), min(cell_stats$n_genes), max(cell_stats$n_genes)))
    cat(sprintf("    平均每个细胞总表达量: %.1f (范围: %.0f-%.0f)\n",mean(cell_stats$total_counts),min(cell_stats$total_counts), max(cell_stats$total_counts))) 
    if (length(mito_genes) > 0) {
      cat(sprintf("    平均线粒体基因百分比: %.2f%% (范围: %.2f%%-%.2f%%)\n",mean(cell_stats$mito_percent) * 100,min(cell_stats$mito_percent) * 100,max(cell_stats$mito_percent) * 100))
    }
    cat("  基因统计:\n")
    cat(sprintf("    平均每个基因在细胞中表达: %.1f 个细胞 (范围: %d-%d)\n",mean(gene_stats$n_cells),min(gene_stats$n_cells), max(gene_stats$n_cells)))
  }
  if (verbose) cat("\n过滤低质量基因...\n")
  n_cells_total <- ncol(expr_matrix)
  if (gene_min_cells < 1) {
    gene_min_cells_abs <- ceiling(gene_min_cells * n_cells_total)
  } else {
    gene_min_cells_abs <- gene_min_cells
  }
  genes_to_keep <- gene_stats$n_cells >= gene_min_cells_abs
  if (verbose) {
    cat(sprintf("  基因过滤条件: 至少在 %d 个细胞中表达 (%.1f%%)\n", gene_min_cells_abs, gene_min_cells * 100))
    cat(sprintf("  保留基因: %d / %d (%.1f%%)\n", sum(genes_to_keep), nrow(expr_matrix),sum(genes_to_keep) / nrow(expr_matrix) * 100))
    if (any(!genes_to_keep)) {
      filtered_genes <- gene_stats[!genes_to_keep, ]
      cat(sprintf("  过滤基因中: %.1f%% 是线粒体基因\n",sum(filtered_genes$is_mito) / nrow(filtered_genes) * 100))
    }
  }
  expr_matrix <- expr_matrix[genes_to_keep, ]
  gene_stats <- gene_stats[genes_to_keep, ]
  if (verbose) cat("\n过滤低质量细胞...\n")
  cell_stats$n_genes <- colSums(expr_matrix > 0)
  cell_stats$total_counts <- colSums(expr_matrix)
  if (length(mito_genes) > 0) {
    remaining_mito <- intersect(mito_genes, rownames(expr_matrix))
    if (length(remaining_mito) > 0) {
      mito_counts <- colSums(expr_matrix[remaining_mito, , drop = FALSE])
      cell_stats$mito_percent <- mito_counts / cell_stats$total_counts
      cell_stats$mito_percent[is.na(cell_stats$mito_percent)] <- 0
    }
  }
  cells_to_keep <- cell_stats$n_genes >= cell_min_genes
  if (!is.null(cell_max_genes)) {cells_to_keep <- cells_to_keep & (cell_stats$n_genes <= cell_max_genes)}
  if (!is.null(cell_max_mito) && length(mito_genes) > 0) {cells_to_keep <- cells_to_keep & (cell_stats$mito_percent <= cell_max_mito)}
  if (verbose) {
    if (!is.null(cell_max_genes)) {cat(sprintf("    最大表达基因数: %d (去除doublet)\n", cell_max_genes))}
    if (!is.null(cell_max_mito) && length(mito_genes) > 0) {cat(sprintf("    最大线粒体百分比: %.1f%%\n", cell_max_mito * 100))}
    cat(sprintf("  保留细胞: %d / %d (%.1f%%)\n", sum(cells_to_keep), ncol(expr_matrix),sum(cells_to_keep) / ncol(expr_matrix) * 100))
  }
  expr_matrix <- expr_matrix[, cells_to_keep]
  cell_stats <- cell_stats[cells_to_keep, ]
  if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {stop("错误: 过滤后矩阵为空，请调整过滤参数")}
  if (verbose) {cat(sprintf("  过滤后矩阵维度: %d 基因 × %d 细胞\n", nrow(expr_matrix), ncol(expr_matrix)))}
  if (remove_mito && length(mito_genes) > 0) {
    if (verbose) cat("\n6. 移除线粒体基因...\n")
    mito_in_matrix <- intersect(mito_genes, rownames(expr_matrix))
    if (length(mito_in_matrix) > 0) {
      expr_matrix <- expr_matrix[!rownames(expr_matrix) %in% mito_in_matrix, ]
      if (verbose) {
        cat(sprintf("  移除了 %d 个线粒体基因\n", length(mito_in_matrix)))
        cat(sprintf("  移除后矩阵维度: %d 基因 × %d 细胞\n", nrow(expr_matrix), ncol(expr_matrix)))
      }
    }
  }
  if (verbose) cat(sprintf("\n标准化: %s...\n", normalization))
  count_matrix <- expr_matrix
  if (normalization == "log") {
    expr_norm <- log1p(expr_matrix)
    if (verbose) cat("  应用log1p转换: log(1 + x)\n")
  } else if (normalization == "zscore") {
    expr_norm <- log1p(expr_matrix)
    if (verbose) cat("  应用log1p转换\n")
  } else if (normalization == "cpm") {
    expr_norm <- t(t(expr_matrix) / colSums(expr_matrix)) * scale_factor
    if (verbose) cat(sprintf("  应用CPM标准化 (缩放因子: %.0f)\n", scale_factor))
  } else if (normalization == "scran") {
    if (!requireNamespace("scran", quietly = TRUE)) {
      stop("请安装scran包: BiocManager::install('scran')")
    }
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = expr_matrix))
    clusters <- scran::quickCluster(expr_matrix)
    sce <- scran::computeSumFactors(sce, clusters = clusters)
    sce <- scater::logNormCounts(sce)
    expr_norm <- as.matrix(logcounts(sce))
    if (verbose) cat("  应用scran标准化 (使用去卷积方法)\n")
  } else if (normalization == "none") {
    expr_norm <- expr_matrix
    if (verbose) cat("  跳过标准化步骤\n")
  } else {
    stop(sprintf("不支持的标准化方法: %s", normalization))
  }
  if (do_scaling && normalization %in% c("log", "zscore")) {
    if (verbose) cat("  应用z-score缩放\n")
    gene_means <- rowMeans(expr_norm)
    gene_sds <- apply(expr_norm, 1, sd)
    gene_sds[gene_sds == 0] <- 1
    expr_norm <- (expr_norm - gene_means) / gene_sds
  }
  if (do_hvg) {
    if (verbose) cat(sprintf("\n8. 选择高变基因 (%d个)...\n", n_hvg))
    if (n_hvg >= nrow(expr_norm)) {
      if (verbose) cat("  警告: 指定的高变基因数量超过总基因数，使用所有基因\n")
      hvg_selected <- 1:nrow(expr_norm)
    } else {
      gene_means <- rowMeans(expr_norm)
      gene_vars <- apply(expr_norm, 1, var)
      gene_dispersion <- gene_vars / (gene_means + 1e-8)
      hvg_selected <- order(gene_vars, decreasing = TRUE)[1:min(n_hvg, length(gene_vars))]
      if (verbose) {
        cat(sprintf("  选择前 %d 个高变基因\n", length(hvg_selected)))
        cat(sprintf("  高变基因的方差范围: %.3f - %.3f\n", min(gene_vars[hvg_selected]), max(gene_vars[hvg_selected])))
      }
    }
    expr_norm <- expr_norm[hvg_selected, ]
    count_matrix <- count_matrix[hvg_selected, ]
  }
  
  if (verbose) {
    cat(sprintf("原始矩阵: %d 基因 × %d 细胞\n", nrow(original_matrix), ncol(original_matrix)))
    cat(sprintf("处理后矩阵: %d 基因 × %d 细胞\n", nrow(expr_norm), ncol(expr_norm)))
    cat(sprintf("基因保留率: %.1f%%\n", nrow(expr_norm) / nrow(original_matrix) * 100))
    cat(sprintf("细胞保留率: %.1f%%\n", ncol(expr_norm) / ncol(original_matrix) * 100))
    cat(sprintf("标准化方法: %s%s\n", normalization,ifelse(do_scaling && normalization %in% c("log", "zscore"), " + z-score", "")))
    if (do_hvg) {cat(sprintf("高变基因: %d 个\n", nrow(expr_norm)))}
  }
  result <- list(processed_matrix = expr_norm,count_matrix = count_matrix,gene_stats = gene_stats,cell_stats = cell_stats,
    parameters = list(gene_min_cells = gene_min_cells,cell_min_genes = cell_min_genes,cell_max_genes = cell_max_genes,cell_max_mito = cell_max_mito,normalization = normalization,scale_factor = scale_factor,do_scaling = do_scaling,remove_mito = remove_mito,mito_pattern = mito_pattern,do_hvg = do_hvg,n_hvg = n_hvg),
    filtering_summary = list(original_genes = nrow(original_matrix),original_cells = ncol(original_matrix),filtered_genes = nrow(expr_norm),filtered_cells = ncol(expr_norm),gene_retention = nrow(expr_norm) / nrow(original_matrix) * 100,cell_retention = ncol(expr_norm) / ncol(original_matrix) * 100),
    mito_genes = mito_genes)
  return(result)
}
#' @param result 预处理结果对象
#' @param output_dir 输出目录
#' @export
basic_visualize_preprocessing <- function(result, output_dir = ".") {
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
  cell_stats <- result$cell_stats
  gene_stats <- result$gene_stats
  png_file <- file.path(output_dir, "quality_control_plots.png")
  png(png_file, width = 1200, height = 1600, res = 150)
  # 设置布局：4行2列
  par(mfrow = c(4, 2), mar = c(4, 4, 3, 2), oma = c(0, 0, 2, 0))
  hist(cell_stats$n_genes, breaks = 50,col = rgb(0.3, 0.5, 0.7, 0.7), border = "white",xlab = "Number of Genes", ylab = "Frequency",main = "A. Genes per Cell", cex.main = 1.2)
  box()
  hist(cell_stats$total_counts, breaks = 50,col = rgb(0.1, 0.6, 0.3, 0.7), border = "white",xlab = "Total Counts", ylab = "Frequency",main = "B. Counts per Cell", cex.main = 1.2)
  box()
  if ("mito_percent" %in% colnames(cell_stats) && 
      any(cell_stats$mito_percent > 0, na.rm = TRUE)) {
    hist(cell_stats$mito_percent * 100, breaks = 50,col = rgb(1, 0.6, 0.1, 0.7), border = "white",xlab = "Mitochondrial %", ylab = "Frequency",main = "C. Mitochondrial Percentage", cex.main = 1.2)
    box()
  } else {
    plot(1, 1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "C. No Mitochondrial Genes")
    text(1, 1, "No mitochondrial genes detected", cex = 1.2)
  }
  hist(gene_stats$n_cells, breaks = 50,col = rgb(0.6, 0.3, 0.7, 0.7), border = "white",xlab = "Number of Cells", ylab = "Frequency",main = "D. Cells per Gene", cex.main = 1.2)
  box()
  hist(gene_stats$mean_expression, breaks = 50,col = rgb(0.7, 0.4, 0.3, 0.7), border = "white",xlab = "Mean Expression", ylab = "Frequency",main = "E. Mean Expression per Gene", cex.main = 1.2)
  box()
  smoothScatter(cell_stats$n_genes, cell_stats$total_counts,xlab = "Number of Genes", ylab = "Total Counts",main = "F. Genes vs Counts Correlation",colramp = colorRampPalette(c("white", "blue", "darkblue")))
  abline(lm(total_counts ~ n_genes, data = cell_stats), col = "red", lwd = 2)
  boxplot(cell_stats$n_genes,col = rgb(0.3, 0.5, 0.7, 0.5),ylab = "Number of Genes",main = "G. Genes per Cell (Boxplot)")
  boxplot(gene_stats$n_cells,col = rgb(0.6, 0.3, 0.7, 0.5),ylab = "Number of Cells",main = "H. Cells per Gene (Boxplot)")
  mtext("Single-cell RNA-seq Quality Control Plots", outer = TRUE, cex = 1.5, font = 2, line = -1)
  dev.off()
  cat(sprintf("质量控制图表已保存到: %s\n", png_file))
}
################################## END ##################################