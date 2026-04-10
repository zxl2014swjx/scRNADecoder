################################## START ##################################
# analysis_celltype.R
# 输入: data_preprocess.txt, 以及跨组织人体成纤维细胞亚型参考图谱
# 输出: 2.Identify_celltype文件夹及其所有分析结果
analysis_celltype <- function(query_matrix = NULL,ref_matrix = NULL,output_dir = "./2.Identify_celltype",project_name = "Example",create_plots = TRUE,plot_width = 12,plot_height = 10,plot_dpi = 300) {
required_packages <- c("ggplot2", "RColorBrewer")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(sprintf("请安装以下包: %s", paste(missing_packages, collapse = ", ")))}
  library(ggplot2)
  library(RColorBrewer)
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)} 
annotate_scRNA <- function(query_matrix=NULL, ref_matrix=NULL, method = "pearson",top_genes = 2000,min_confidence = 0.1) {
  query_matrix<-read.table(query_matrix,header=T,sep="\t",row.names=1)
  ref_matrix = read.table(unz("Reference.zip", "Cand_Major8_subtype16_average.txt"),header=T,sep="\t",row.names=1)
  if (!is.matrix(query_matrix)) query_matrix <- as.matrix(query_matrix)
  if (!is.matrix(ref_matrix)) ref_matrix <- as.matrix(ref_matrix)
  cat(sprintf("查询数据: %d 个基因, %d 个细胞\n", nrow(query_matrix), ncol(query_matrix)))
  cat(sprintf("参考数据: %d 个基因, %d 个细胞类型\n", nrow(ref_matrix), ncol(ref_matrix)))
  common_genes <- intersect(rownames(query_matrix), rownames(ref_matrix))
  if (length(common_genes) == 0) {stop("错误: 查询数据和参考数据没有共同的基因名")}
  cat(sprintf("共有 %d 个共同基因\n", length(common_genes)))
  query <- query_matrix[common_genes, , drop = FALSE]
  ref <- ref_matrix[common_genes, , drop = FALSE]
  if (!is.null(top_genes) && top_genes < nrow(query)) {
    gene_vars <- apply(ref, 1, var)
    selected_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_genes, length(gene_vars))]
    query <- query[selected_genes, , drop = FALSE]
    ref <- ref[selected_genes, , drop = FALSE]
    cat(sprintf("已选择 %d 个特征基因\n", length(selected_genes)))
  }
  n_cells <- ncol(query)
  n_celltypes <- ncol(ref)
  similarity_matrix <- matrix(0, nrow = n_cells, ncol = n_celltypes,dimnames = list(colnames(query), colnames(ref)))
  if (method %in% c("spearman", "pearson")) {
    library(doParallel)
    library(foreach)
    n_cores <- min(parallel::detectCores() - 1, 10)
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    cat(sprintf("使用 %d 个核心进行并行计算...\n", n_cores))
    results <- foreach(i = 1:n_cells, .combine = rbind, .packages = c("matrixStats")) %dopar% {
      cell_expr <- query[, i]
      cor_results <- apply(ref, 2, function(ref_expr) {
        if (method == "spearman") {
          cor(cell_expr, ref_expr, method = "spearman")
        } else {
          cor(cell_expr, ref_expr, method = "pearson")
        }
      })
      cor_results
    }
    stopCluster(cl)
    similarity_matrix <- results
    rownames(similarity_matrix) <- colnames(query)
    colnames(similarity_matrix) <- colnames(ref)
  } else if (method == "cosine") {
    for (i in 1:n_cells) {
      cell_expr <- query[, i]
      for (j in 1:n_celltypes) {
        ref_expr <- ref[, j]
        similarity <- sum(cell_expr * ref_expr) / 
          (sqrt(sum(cell_expr^2)) * sqrt(sum(ref_expr^2)) + 1e-8)
        similarity_matrix[i, j] <- similarity
      }
      if (i %% 100 == 0) cat(sprintf("已处理 %d/%d 个细胞\n", i, n_cells))
    }
  } else if (method == "knn") {
    cat("使用KNN方法...\n")
    query_t <- t(query)
    ref_t <- t(ref)
    library(FNN)
    knn_result <- get.knnx(ref_t, query_t, k = 5)
    similarity_matrix <- 1 / (1 + knn_result$nn.dist[, 1])
  } else {
    stop(sprintf("不支持的方法: %s", method))
  }
  predictions <- character(n_cells)
  confidence_scores <- numeric(n_cells)
  for (i in 1:n_cells) {
    cell_scores <- similarity_matrix[i, ]
    max_score <- max(cell_scores, na.rm = TRUE)
    best_type <- names(which.max(cell_scores))
    predictions[i] <- best_type
    confidence_scores[i] <- max_score
  }
  uncertain_cells <- confidence_scores < min_confidence
  if (any(uncertain_cells)) {
    predictions[uncertain_cells] <- "Unknown"
    cat(sprintf("有 %d 个细胞的置信度低于阈值，标记为Unknown\n", sum(uncertain_cells)))
  }
  result <- data.frame(
    cell_id = colnames(query),
    predicted_celltype = predictions,
    confidence = confidence_scores,
    is_uncertain = uncertain_cells,
    stringsAsFactors = FALSE)
  score_df <- as.data.frame(similarity_matrix)
  colnames(score_df) <- paste0("score_", colnames(score_df))
  final_result <- cbind(result, score_df)
  rownames(final_result) <- NULL
  cat("注释完成！\n")
  cat(sprintf("预测了 %d 个细胞，%d 个细胞类型\n", n_cells, length(unique(predictions[!uncertain_cells]))))
  return(list(
    predictions = final_result,
    similarity_matrix = similarity_matrix,
    parameters = list(
      method = method,
      top_genes = ifelse(is.null(top_genes), "all", top_genes),
      min_confidence = min_confidence,
      n_cells = n_cells,
      n_celltypes = n_celltypes
    )
  ))
}
  start_time <- Sys.time()
  annotation_result <- annotate_scRNA(
    query_matrix = query_matrix,
    ref_matrix = ref_matrix,
    method = "pearson",
    top_genes = 2000,
    min_confidence = 1e-5)
  end_time <- Sys.time()
  cat(sprintf("注释完成，耗时: %.1f 分钟\n", as.numeric(difftime(end_time, start_time, units = "mins"))))
create_basic_plots <- function(annotation_result, plot_dir, project_name,plot_width, plot_height, plot_dpi) {
  df <- annotation_result$predictions
  p1 <- ggplot(df, aes(x = confidence)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "black") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
    geom_vline(xintercept = 0.3, linetype = "dashed", color = "orange", size = 1) +
    labs(title = "Distribution of Confidence Scores",subtitle = paste("Mean:", round(mean(df$confidence), 3),"Median:", round(median(df$confidence), 3)),x = "Confidence Score",y = "Number of Cells") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          panel.grid.minor = element_blank())
    ggsave(file.path(plot_dir, paste0(project_name, "_confidence_distribution.png")), p1, width = plot_width, height = plot_height, dpi = plot_dpi)
  p2 <- ggplot(df, aes(x = predicted_celltype, y = confidence, fill = predicted_celltype)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 0.5) +
    labs(title = "Confidence Scores by Cell Type",x = "Cell Type",y = "Confidence Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          panel.grid.major.x = element_blank())
    ggsave(file.path(plot_dir, paste0(project_name, "_confidence_by_celltype.png")), p2, width = plot_width, height = plot_height, dpi = plot_dpi)     
  celltype_counts <- as.data.frame(table(df$predicted_celltype))
  colnames(celltype_counts) <- c("CellType", "Count")
  celltype_counts <- celltype_counts[order(celltype_counts$Count, decreasing = TRUE), ]
  celltype_counts$Percentage <- round(celltype_counts$Count / sum(celltype_counts$Count) * 100, 1)
  celltype_counts$CellType <- factor(celltype_counts$CellType, levels = celltype_counts$CellType)
  p3 <- ggplot(celltype_counts, aes(x = reorder(CellType, -Count), y = Count, fill = CellType)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), vjust = -0.5, size = 3) +
    labs(title = "Cell Type Distribution",x = "Cell Type",y = "Number of Cells") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          panel.grid.major.x = element_blank())
  conf_summary <- aggregate(confidence ~ predicted_celltype, df, FUN = function(x) c(mean = mean(x), sd = sd(x)))
  conf_summary <- do.call(data.frame, conf_summary)
  colnames(conf_summary) <- c("CellType", "mean_confidence", "sd_confidence")
  ggsave(file.path(plot_dir, paste0(project_name, "_celltype_distribution.png")), p3, width = plot_width, height = plot_height, dpi = plot_dpi)
  p4 <- ggplot(conf_summary, aes(x = reorder(CellType, -mean_confidence), y = mean_confidence, fill = mean_confidence)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean_confidence - sd_confidence, ymax = mean_confidence + sd_confidence), width = 0.2) +
    scale_fill_gradient(low = "blue", high = "red") +
    labs(title = "Mean Confidence by Cell Type",x = "Cell Type",y = "Mean Confidence (± SD)",fill = "Confidence") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          panel.grid.major.x = element_blank())
    ggsave(file.path(plot_dir, paste0(project_name, "_mean_confidence_heatmap.png")), p4, width = plot_width, height = plot_height, dpi = plot_dpi)
}
create_advanced_plots <- function(annotation_result, plot_dir, project_name,plot_width, plot_height, plot_dpi) {
  df <- annotation_result$predictions
  sim_matrix <- annotation_result$similarity_matrix
  if (nrow(sim_matrix) > 50) {
    plot_sim <- sim_matrix[1:min(50, nrow(sim_matrix)), ]
  } else {
    plot_sim <- sim_matrix
  }
  plot_sim_long <- reshape2::melt(plot_sim)
  colnames(plot_sim_long) <- c("Cell", "CellType", "Similarity")
  p1 <- ggplot(plot_sim_long, aes(x = CellType, y = Cell, fill = Similarity)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Similarity Matrix (Top 50 Cells)",x = "Reference Cell Types",y = "Query Cells") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    ggsave(file.path(plot_dir, paste0(project_name, "_similarity_heatmap.png")), p1, width = plot_width, height = plot_height, dpi = plot_dpi)
  top_similarities <- apply(sim_matrix, 1, max)
  celltype_top_sim <- data.frame(CellType = df$predicted_celltype,TopSimilarity = top_similarities)
  p2 <- ggplot(celltype_top_sim, aes(x = CellType, y = TopSimilarity, fill = CellType)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", size = 0.5) +
    labs(title = "Top Similarity Score Distribution by Cell Type",x = "Cell Type",y = "Top Similarity Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          panel.grid.major.x = element_blank())
    ggsave(file.path(plot_dir, paste0(project_name, "_top_similarity_violin.png")), p2, width = plot_width, height = plot_height, dpi = plot_dpi)
  celltype_avg_sim <- aggregate(. ~ df$predicted_celltype, as.data.frame(sim_matrix), mean)
  rownames(celltype_avg_sim) <- celltype_avg_sim[, 1]
  celltype_avg_sim <- as.matrix(celltype_avg_sim[, -1])
  celltype_avg_sim_long <- reshape2::melt(celltype_avg_sim)
  colnames(celltype_avg_sim_long) <- c("PredictedType", "ReferenceType", "AvgSimilarity")
  p3 <- ggplot(celltype_avg_sim_long, aes(x = ReferenceType, y = PredictedType, fill = AvgSimilarity)) +
    geom_tile() +
    geom_text(aes(label = round(AvgSimilarity, 2)), color = "black", size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5, limits = c(0, 1)) +
    labs(title = "Average Similarity: Predicted vs Reference Cell Types",x = "Reference Cell Types",y = "Predicted Cell Types",fill = "Avg\nSimilarity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    ggsave(file.path(plot_dir, paste0(project_name, "_avg_similarity_heatmap.png")), p3, width = plot_width, height = plot_height, dpi = plot_dpi)
  plot_data <- data.frame(Cell = df$cell_id,CellType = df$predicted_celltype,Confidence = df$confidence,TopSimilarity = top_similarities)
  p4 <- ggplot(plot_data, aes(x = TopSimilarity, y = Confidence, color = CellType)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = "Confidence vs Top Similarity Score",x = "Top Similarity Score",y = "Confidence Score") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    ggsave(file.path(plot_dir, paste0(project_name, "_confidence_vs_similarity.png")), p4, width = plot_width, height = plot_height, dpi = plot_dpi)
}
save_results <- function(annotation_result, output_dir, project_name) {
  result_file <- file.path(output_dir, sprintf("%s_annotations.csv", project_name))
  write.csv(annotation_result$predictions, result_file, row.names = FALSE)
  sim_file <- file.path(output_dir, sprintf("%s_similarity_matrix.rds", project_name))
  saveRDS(annotation_result$similarity_matrix, sim_file)
  rdata_file <- file.path(output_dir, sprintf("%s_results.RData", project_name))
  save(annotation_result, file = rdata_file)
}
  if (create_plots) {
    plot_dir <- output_dir
    if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
    create_basic_plots(annotation_result, plot_dir, project_name,  plot_width, plot_height, plot_dpi)
    create_advanced_plots(annotation_result, plot_dir, project_name,plot_width, plot_height, plot_dpi)
    save_results(annotation_result, output_dir, project_name)
  return(annotation_result)
}
}
################################## END ##################################