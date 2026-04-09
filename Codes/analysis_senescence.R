################################## START ##################################
#' analysis_senescence.R
# 输入: data_preprocess.txt, 用户自定义features基因列表
# 输出: 3.Senescence_score文件夹及其所有分析结果
# 设置工作目录和加载包
cat("设置工作环境...\n")
setwd(".")  # 设置工作目录
# 检查并安装必要的包
required_packages <- c("reshape2", "pheatmap", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("安装缺失的包:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, dependencies = TRUE)}
# 加载所有包
suppressPackageStartupMessages({
  library(reshape2)
  library(pheatmap)
  library(ggplot2)
})
cat("\n创建输出目录...\n")
dir.create("3.Senescence_score", showWarnings = FALSE)
cat("\n1. 正在加载数据...\n")
expr_mat<-read.table("1.Data_preprocess/data_preprocessed.txt",header=T,sep="\t",row.names=1,check.names=F)
gene<-read.table("Example_Data/senescence_genesets.txt",header=T,sep="\t")
gene_modules <- list(
Senescence = gene[which(gene$Type=="Senescence"),]$Genes,
Stemness = gene[which(gene$Type=="Stemness"),]$Genes)
calculate_module_score <- function(expr_matrix, features, 
  pool = NULL, nbin = 24, ctrl = 100, name = "Module", seed = 1, search = FALSE, verbose = TRUE) {
  if (!is.null(seed)) {set.seed(seed = seed)}
  if (!is.matrix(expr_matrix) && !inherits(expr_matrix, "Matrix")) {
    expr_matrix <- as.matrix(expr_matrix)}
  if (is.null(rownames(expr_matrix))) {
    stop("expr_matrix must have rownames (gene names)")}
  if (is.null(colnames(expr_matrix))) {
    colnames(expr_matrix) <- paste0("Cell", 1:ncol(expr_matrix))}
  if (is.list(features)) {
    features_list <- features
  } else if (is.character(features)) {
    features_list <- list(Module1 = features)
  } else {
    stop("features must be either a character vector or a list of character vectors")
  }
  if (is.null(features_list) || length(features_list) == 0) {
    stop("Missing input feature list")
  }
  features_old <- features_list
  features_list <- lapply(features_list, function(genes) {
    missing_genes <- setdiff(genes, rownames(expr_matrix))
    if (length(missing_genes) > 0 && verbose) {
      message("The following genes are not present in the matrix: ", 
              paste(missing_genes, collapse = ", "))
    }
    return(intersect(genes, rownames(expr_matrix)))
  })
  valid_lengths <- sapply(features_list, length) > 0
  if (!all(valid_lengths)) {
    stop("The following feature lists have no valid genes present in the matrix: ", 
         paste(names(features_list)[!valid_lengths], collapse = ", "))
  }
  if (verbose) {
    message(paste("Processing", length(features_list), "gene modules..."))
    for (i in seq_along(features_list)) {
      module_name <- names(features_list)[i]
      n_genes <- length(features_list[[i]])
      message(paste("  Module", i, "(", module_name, "):", n_genes, "genes"))
    }
  }
  pool <- pool %||% rownames(expr_matrix)
  # 计算每个基因在所有样本中的平均表达量
  if (verbose) message("Calculating gene average expression...")
  gene_means <- Matrix::rowMeans(expr_matrix[pool, , drop = FALSE])
  gene_means <- gene_means[order(gene_means)]
  # 将基因按表达量分箱
  if (verbose) message("Binning genes by expression level...")
  # 添加微小噪声以避免平局
  gene_means_noisy <- gene_means + rnorm(n = length(gene_means)) / 1e30
  # 使用cut_number将基因分成nbin个箱
  gene_bins <- cut(x = gene_means_noisy, breaks = nbin, labels = FALSE, include.lowest = TRUE)
  names(gene_bins) <- names(gene_means)
  n_modules <- length(features_list)
  # 为每个模块选择背景基因
  if (verbose) message("Selecting control genes...")
  ctrl_genes <- vector("list", length = n_modules)
  for (i in 1:n_modules) {
    module_genes <- features_list[[i]]
    ctrl_genes[[i]] <- character(0)
    for (gene in module_genes) {
      if (gene %in% names(gene_bins)) {
        # 找到与目标基因在同一箱内的所有基因
        same_bin <- names(gene_bins)[gene_bins == gene_bins[gene]]
        # 移除目标基因本身
        same_bin <- setdiff(same_bin, gene)
        # 从同一箱中随机选择ctrl个基因作为背景
        if (length(same_bin) >= ctrl) {
          selected <- sample(same_bin, size = ctrl, replace = FALSE)
        } else {
          # 如果同一箱内基因不够，使用所有可用基因
          selected <- same_bin
          if (verbose && length(same_bin) < ctrl) {
            message(paste("  Warning: Only", length(same_bin), "control genes available for", gene,  "(requested", ctrl, ")"))
          }
        }
        ctrl_genes[[i]] <- c(ctrl_genes[[i]], selected)
      }
    }
    ctrl_genes[[i]] <- unique(ctrl_genes[[i]])
  }
  # 计算背景基因的平均表达
  if (verbose) message("Calculating control gene scores...")
  ctrl_scores <- matrix(0, nrow = n_modules, ncol = ncol(expr_matrix))
  for (i in 1:n_modules) {
    ctrl_genes_use <- ctrl_genes[[i]]
    if (length(ctrl_genes_use) > 0) {
      ctrl_scores[i, ] <- Matrix::colMeans(expr_matrix[ctrl_genes_use, , drop = FALSE])
    } else {
      ctrl_scores[i, ] <- 0
    }
  }
  # 计算目标基因模块的平均表达
  if (verbose) message("Calculating module gene scores...")
  module_scores <- matrix(0, nrow = n_modules, ncol = ncol(expr_matrix))
  for (i in 1:n_modules) {
    module_genes <- features_list[[i]]
    module_scores[i, ] <- Matrix::colMeans(expr_matrix[module_genes, , drop = FALSE])
  }
  # 计算调整后的模块得分：模块基因平均表达 - 背景基因平均表达
  if (verbose) message("Calculating adjusted module scores...")
  adjusted_scores <- module_scores - ctrl_scores
  # 设置行列名
  if (is.null(names(features_list))) {
    rownames(adjusted_scores) <- paste0(name, 1:n_modules)
  } else {
    rownames(adjusted_scores) <- names(features_list)
  }
  colnames(adjusted_scores) <- colnames(expr_matrix)
  # 转换为数据框（样本为行，模块为列）
  result_df <- as.data.frame(t(adjusted_scores))
  if (verbose) message("Module score calculation completed!")
  return(result_df)
}
module_scores <- calculate_module_score(expr_matrix = expr_mat,features = gene_modules,nbin = 24,ctrl = 50,verbose = TRUE)
write.table(module_scores,"3.Senescence_score/module_scores.txt",sep="\t",quote=F)

#' 可视化模块得分
#' 绘制模块得分的可视化图形
#' @param module_scores 模块得分
#' @param output_file 输出文件路径
#' @param plot_type 图形类型
visualize_module_scores <- function(module_scores, output_file = NULL, plot_type = "heatmap") {
  if (!require("ggplot2", quietly = TRUE) || !require("pheatmap", quietly = TRUE)) {install.packages(c("ggplot2", "pheatmap"))}
  library(ggplot2)
  library(pheatmap)
  if (!is.null(output_file)) {pdf(output_file, width = 10, height = 8)}
  if (plot_type == "heatmap") {
    pheatmap(t(module_scores),main = "Module Scores Heatmap",color = colorRampPalette(c("blue", "white", "red"))(50),show_colnames = ncol(module_scores) <= 50,fontsize_row = 8,fontsize_col = 8)  
  } else if (plot_type == "boxplot") {
    plot_data <- reshape2::melt(as.matrix(module_scores))
    colnames(plot_data) <- c("Sample", "Module", "Score")
    p <- ggplot(plot_data, aes(x = Module, y = Score, fill = Module)) +geom_boxplot() +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ggtitle("Module Score Distribution")
    print(p)
  } else if (plot_type == "violin") {
    plot_data <- reshape2::melt(as.matrix(module_scores))
    colnames(plot_data) <- c("Sample", "Module", "Score")
    p <- ggplot(plot_data, aes(x = Module, y = Score, fill = Module)) +geom_violin(scale = "width") +geom_boxplot(width = 0.1, fill = "white") +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ggtitle("Module Score Distribution (Violin Plot)")
    print(p)
  }
  if (!is.null(output_file)) {dev.off()
  }
}
pdf("./3.Senescence_score/heatmap.pdf")
visualize_module_scores(module_scores, plot_type = "heatmap")
dev.off()
pdf("./3.Senescence_score/boxplot.pdf")
visualize_module_scores(module_scores, plot_type = "boxplot")
dev.off()
################################## END ##################################