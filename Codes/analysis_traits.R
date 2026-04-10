################################## START ##################################
# analysis_traits.R
# 输入: gene_traits_matrix.csv, normalized_gene_samples.csv, fibroblast_samples.csv
# 输出: 4.Traits_association文件夹及其所有分析结果
required_packages <- c("tidyverse", "reshape2", "pheatmap", "ggplot2", "ggdendro", 
                      "corrplot", "RColorBrewer", "viridis", "ggpubr", "ComplexHeatmap")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("安装缺失的包:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(pheatmap)
  library(ggplot2)
  library(ggdendro)
  library(corrplot)
  library(RColorBrewer)
  library(viridis)
  library(ggpubr)
})
dir.create("4.Traits_association", showWarnings = FALSE)
gene_traits_df<-read.table(unz("Reference.zip", "Result_scDRS_correlations.txt"),header=T,sep="\t",row.names=1)
gene_samples_df<-read.table("1.Data_preprocess/data_preprocessed.txt",header=T,sep="\t",row.names=1)
fibroblast_df<-read.table("2.Identify_celltype/Example_annotations.csv",
  header=T,sep=",",row.names=1)
fibroblast_df<-fibroblast_df[,-c(1:3)]
fibroblast_df<-t(fibroblast_df[,1:16])
common_samples <- intersect(colnames(gene_samples_df), colnames(fibroblast_df))
cat("  基因表达和成纤维细胞数据中的共同样本数: ", length(common_samples), "\n")
if(length(common_samples) == 0) {
  stop("错误: 没有共同的样本!")
}
gene_samples_df <- gene_samples_df[, common_samples, drop = FALSE]
fibroblast_df <- fibroblast_df[, common_samples, drop = FALSE]
common_genes <- intersect(rownames(gene_traits_df), rownames(gene_samples_df))
cat("  基因-性状矩阵和表达矩阵中的共同基因数: ", length(common_genes), "\n")
if(length(common_genes) == 0) {
  stop("错误: 没有共同的基因!")
}
gene_traits_df <- gene_traits_df[common_genes, , drop = FALSE]
gene_samples_df <- gene_samples_df[common_genes, , drop = FALSE]
dim(gene_traits_df)
gene_traits_df[1:6,1:6]
dim(gene_samples_df)
gene_samples_df[1:6,1:6]
dim(fibroblast_df)
fibroblast_df[1:6,1:6]

predict_traits <- function(gene_expression, gene_traits_matrix, fibroblast_data) {
  traits_prediction <- matrix(0, nrow = ncol(gene_expression),ncol = ncol(gene_traits_matrix),dimnames = list(colnames(gene_expression),colnames(gene_traits_matrix)))
  for(trait in colnames(gene_traits_matrix)) {
    trait_weights <- gene_traits_matrix[, trait, drop = FALSE]
    valid_genes <- rownames(trait_weights)[trait_weights[,1] != 0]
    valid_genes <- intersect(valid_genes, rownames(gene_expression))
    if(length(valid_genes) > 0) {
      expression_subset <- gene_expression[valid_genes, , drop = FALSE]
      weights_subset <- trait_weights[valid_genes, 1]
      expression_normalized <- t(apply(expression_subset, 1, function(x) {
        (x - mean(x)) / sd(x)
      }))
      for(sample in colnames(gene_expression)) {
        traits_prediction[sample, trait] <- sum(expression_normalized[, sample] * weights_subset)
      }
    }
  }
  traits_pred_df <- as.data.frame(traits_prediction)
  if(nrow(fibroblast_data) > 0) {
    fibroblast_proportion <- as.numeric(fibroblast_data[1, ])
    names(fibroblast_proportion) <- colnames(fibroblast_data)
  } else {
    fibroblast_proportion <- rep(0, ncol(fibroblast_data))
  }
  write.csv(traits_pred_df, "4.Traits_association/traits_prediction_raw.csv")
  write.csv(data.frame(Sample = names(fibroblast_proportion), Fibroblast_Proportion = fibroblast_proportion),"4.Traits_association/fibroblast_proportion.csv", row.names = FALSE)
  return(list(traits = traits_pred_df, fibroblast = fibroblast_proportion))
}
pre_traits_pred<-predict_traits(gene_samples_df,gene_traits_df,fibroblast_df)
traits_pred<-pre_traits_pred$traits
fibroblast_prop<-pre_traits_pred$fibroblast

classify_traits <- function(traits_df, threshold = 1.0) {
  trait_classes <- data.frame(
    Sample = rownames(traits_df),
    High_Risk_Traits = character(nrow(traits_df)),
    Protective_Traits = character(nrow(traits_df)),
    Neutral_Traits = character(nrow(traits_df)),
    Top_Predicted_Trait = character(nrow(traits_df)),
    Trait_Risk_Score = numeric(nrow(traits_df)),
    stringsAsFactors = FALSE
  )
  for(i in 1:nrow(traits_df)) {
    sample_name <- rownames(traits_df)[i]
    sample_traits <- traits_df[i, ]
    high_risk_idx <- which(sample_traits > threshold)
    protective_idx <- which(sample_traits < -threshold)
    neutral_idx <- which(sample_traits >= -threshold & sample_traits <= threshold)
    trait_classes$High_Risk_Traits[i] <- paste(names(sample_traits)[high_risk_idx], collapse = ",")
    trait_classes$Protective_Traits[i] <- paste(names(sample_traits)[protective_idx], collapse = ",")
    trait_classes$Neutral_Traits[i] <- paste(names(sample_traits)[neutral_idx], collapse = ",")
    trait_classes$Top_Predicted_Trait[i] <- names(which.max(sample_traits))
    trait_classes$Trait_Risk_Score[i] <- mean(as.numeric(sample_traits))
  }
  rownames(trait_classes) <- trait_classes$Sample
  return(trait_classes)
}
trait_classes <- classify_traits(traits_pred)
write.csv(trait_classes, "4.Traits_association/trait_classification.csv", row.names = FALSE)

create_heatmap <- function(traits_df, output_path) {
  plot_data <- as.matrix(t(traits_df)) 
  if(nrow(plot_data) > 50) {
    trait_vars <- apply(plot_data, 1, var)
    top_traits <- names(sort(trait_vars, decreasing = TRUE))[1:50]
    plot_data <- plot_data[top_traits, ]
  }
  if(ncol(plot_data) > 100) {
    set.seed(123)
    if(ncol(plot_data) > 100) {
      sampled_cols <- sample(1:ncol(plot_data), 100)
      plot_data <- plot_data[, sampled_cols]
    }
  }
  pheatmap(plot_data,main = "Predicted Traits Heatmap",color = colorRampPalette(c("blue", "white", "red"))(100),scale = "row", cluster_rows = TRUE,cluster_cols = TRUE,show_colnames = ncol(plot_data) <= 50,show_rownames = TRUE,fontsize_row = 8,fontsize_col = 8,filename = output_path,width = 12,height = 10)
  cat("  热图已保存到:", output_path, "\n")
}
create_heatmap(traits_pred, "4.Traits_association/traits_heatmap.png")

create_trait_distribution <- function(traits_df, output_path) {
  trait_vars <- apply(traits_df, 2, var)
  top_traits <- names(sort(trait_vars, decreasing = TRUE))[1:min(20, ncol(traits_df))]
  plot_data <- traits_df[, top_traits] %>%
    as.data.frame() %>%
    mutate(Sample = rownames(.)) %>%
    pivot_longer(cols = -Sample, names_to = "Trait", values_to = "Score")
  p <- ggplot(plot_data, aes(x = reorder(Trait, Score, FUN = median), y = Score)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
    coord_flip() +
    labs(title = "Distribution of Top Trait Predictions",x = "Trait",y = "Prediction Score") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),panel.grid.major = element_line(color = "gray90"),panel.grid.minor = element_blank())
  ggsave(output_path, p, width = 10, height = 8, dpi = 300)
  cat("  性状分布图已保存到:", output_path, "\n")
}
create_trait_distribution(traits_pred, "4.Traits_association/trait_distribution.png")

create_fibroblast_correlation <- function(traits_df, fibroblast_prop, output_path) {
  correlations <- sapply(traits_df, function(x) {
    cor(x, fibroblast_prop, use = "complete.obs")})
  top_n <- min(15, length(correlations))
  top_corr <- sort(abs(correlations), decreasing = TRUE)[1:top_n]
  top_traits <- names(top_corr)
  plot_data <- data.frame(
    Trait = top_traits,
    Correlation = correlations[top_traits],
    AbsCorrelation = abs(correlations[top_traits])
  )
  plot_data <- plot_data[order(plot_data$Correlation), ]
  plot_data$Trait <- factor(plot_data$Trait, levels = plot_data$Trait)
  p <- ggplot(plot_data, aes(x = Trait, y = Correlation, fill = Correlation > 0)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),name = "Correlation",labels = c("Positive", "Negative")) +
    coord_flip() +
    labs(title = "Top Traits Correlated with Fibroblast Proportion",x = "Trait",y = "Correlation Coefficient") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),legend.position = "right")
  ggsave(output_path, p, width = 10, height = 8, dpi = 300)
  cat("  成纤维细胞相关性图已保存到:", output_path, "\n")
}
create_fibroblast_correlation(traits_pred, fibroblast_prop, "4.Traits_association/fibroblast_correlation.png")

create_risk_score_plot <- function(trait_classes, output_path) {
  plot_data <- trait_classes %>%
    arrange(Trait_Risk_Score) %>%
    mutate(Sample_Index = 1:n(),Risk_Level = ifelse(Trait_Risk_Score > 0, "High Risk", "Low Risk"))
  p <- ggplot(plot_data, aes(x = reorder(Sample, Trait_Risk_Score), y = Trait_Risk_Score, fill = Risk_Level)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    scale_fill_manual(values = c("High Risk" = "red", "Low Risk" = "green")) +
    labs(title = "Trait Risk Scores Across Samples",x = "Samples (sorted by risk score)",y = "Average Trait Risk Score",fill = "Risk Level") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),legend.position = "top",panel.grid.major.x = element_blank())
  ggsave(output_path, p, width = 12, height = 6, dpi = 300)
  cat("  风险评分图已保存到:", output_path, "\n")
}
create_risk_score_plot(trait_classes, "4.Traits_association/risk_scores.png")

create_trait_network <- function(traits_df, output_path, top_n = 30) {
  trait_corr <- cor(traits_df, use = "complete.obs")
  if(ncol(trait_corr) > top_n) {
    mean_abs_corr <- apply(abs(trait_corr), 2, mean)
    top_traits <- names(sort(mean_abs_corr, decreasing = TRUE))[1:top_n]
    corr_matrix <- trait_corr[top_traits, top_traits]
  } else {
    corr_matrix <- trait_corr
  }
  png(output_path, width = 1200, height = 1000, res = 150)
  corrplot(corr_matrix,method = "color",type = "full",order = "hclust",tl.col = "black",tl.srt = 45,addrect = 3,rect.col = "black",rect.lwd = 2,cl.pos = "r",cl.ratio = 0.1,cl.align.text = "l",cl.cex = 0.8,number.cex = 0.7,mar = c(0, 0, 2, 0))
  title(main = "Trait Correlation Network", cex.main = 1.5)
  dev.off()
  cat("  性状相关性网络图已保存到:", output_path, "\n")
}
create_trait_network(traits_pred, "4.Traits_association/trait_correlation_network.png")

generate_analysis_report <- function(traits_df, trait_classes, fibroblast_prop, output_path) {
  report_lines <- c()
  report_lines <- c(report_lines, paste(rep("=", 60), collapse = ""))
  report_lines <- c(report_lines, "TRAIT PREDICTION ANALYSIS REPORT")
  report_lines <- c(report_lines, paste(rep("=", 60), collapse = ""))
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "1. BASIC INFORMATION")
  report_lines <- c(report_lines, sprintf("   Total samples analyzed: %d", nrow(traits_df)))
  report_lines <- c(report_lines, sprintf("   Total traits predicted: %d", ncol(traits_df)))
  report_lines <- c(report_lines, sprintf("   Date of analysis: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "2. TRAIT STATISTICS")
  high_risk_counts <- sapply(strsplit(trait_classes$High_Risk_Traits, ","), function(x) ifelse(x[1] == "", 0, length(x)))
  report_lines <- c(report_lines, sprintf("   Average high-risk traits per sample: %.2f", mean(high_risk_counts)))
  report_lines <- c(report_lines, sprintf("   Samples with high-risk traits: %d", sum(high_risk_counts > 0)))
  protective_counts <- sapply(strsplit(trait_classes$Protective_Traits, ","), function(x) ifelse(x[1] == "", 0, length(x)))
  report_lines <- c(report_lines, sprintf("   Average protective traits per sample: %.2f", mean(protective_counts)))
  report_lines <- c(report_lines, "")
  top_traits <- sort(table(trait_classes$Top_Predicted_Trait), decreasing = TRUE)[1:5]
  report_lines <- c(report_lines, "3. MOST FREQUENTLY PREDICTED TOP TRAITS")
  for(i in 1:length(top_traits)) {
    trait_name <- names(top_traits)[i]
    count <- top_traits[i]
    percentage <- (count / nrow(traits_df)) * 100
    report_lines <- c(report_lines, sprintf("   %s: %d samples (%.1f%%)", trait_name, count, percentage))
  }
  report_lines <- c(report_lines, "")
  high_risk_samples <- trait_classes %>% 
    filter(Trait_Risk_Score > 1) %>%
    arrange(desc(Trait_Risk_Score))
  report_lines <- c(report_lines, "4. HIGH-RISK SAMPLES IDENTIFICATION")
  report_lines <- c(report_lines, sprintf("   Samples with risk score > 1: %d", nrow(high_risk_samples)))
  if(nrow(high_risk_samples) > 0) {
    report_lines <- c(report_lines, "   Top 5 high-risk samples:")
    for(i in 1:min(5, nrow(high_risk_samples))) {
      report_lines <- c(report_lines, sprintf("   %d. %s: Risk score = %.3f", i, high_risk_samples$Sample[i],high_risk_samples$Trait_Risk_Score[i]))
    }
  }
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "5. FIBROBLAST CORRELATION ANALYSIS")
  correlations <- sapply(traits_df, function(x) {
    cor(x, fibroblast_prop, use = "complete.obs")
  })
  top_pos_corr <- sort(correlations, decreasing = TRUE)[1:3]
  top_neg_corr <- sort(correlations)[1:3]
  report_lines <- c(report_lines, "   Traits most positively correlated with fibroblast proportion:")
  for(i in 1:length(top_pos_corr)) {
    report_lines <- c(report_lines, sprintf("   - %s: r = %.3f", names(top_pos_corr)[i], top_pos_corr[i]))
  }
  report_lines <- c(report_lines, "   Traits most negatively correlated with fibroblast proportion:")
  for(i in 1:length(top_neg_corr)) {
    report_lines <- c(report_lines, sprintf("   - %s: r = %.3f", names(top_neg_corr)[i], top_neg_corr[i]))
  }
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "6. OUTPUT FILES")
  report_lines <- c(report_lines, "   The following files have been generated:")
  report_lines <- c(report_lines, "   4.Traits_association/traits_prediction_raw.csv - Raw trait predictions")
  report_lines <- c(report_lines, "   4.Traits_association/trait_classification.csv - Trait classification results")
  report_lines <- c(report_lines, "   4.Traits_association/fibroblast_proportion.csv - Fibroblast proportions")
  report_lines <- c(report_lines, "   4.Traits_association/ - Directory containing all visualizations")
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, paste(rep("=", 60), collapse = ""))
  report_lines <- c(report_lines, "ANALYSIS COMPLETE")
  report_lines <- c(report_lines, paste(rep("=", 60), collapse = ""))
  writeLines(report_lines, con = output_path)
  cat(paste(report_lines, collapse = "\n"), "\n")
  return(report_lines)
}
report <- generate_analysis_report(traits_pred, trait_classes, fibroblast_prop,  "4.Traits_association/analysis_report.txt")
write.csv(traits_pred, "4.Traits_association/traits_predictions_detailed.csv")
top_traits_per_sample <- data.frame(Sample = rownames(traits_pred),Top_1 = NA, Top_2 = NA, Top_3 = NA, Top_4 = NA, Top_5 = NA,stringsAsFactors = FALSE)
for(i in 1:nrow(traits_pred)) {
  sample_name <- rownames(traits_pred)[i]
  top_5 <- sort(traits_pred[i, ], decreasing = TRUE)[1:5]
  for(j in 1:5) {
    if(j <= length(top_5)) {
      trait_name <- names(top_5)[j]
      score <- top_5[j]
      top_traits_per_sample[i, paste0("Top_", j)] <- sprintf("%s (%.3f)", trait_name, score)
    }
  }
}
write.csv(top_traits_per_sample, "4.Traits_association/top_traits_per_sample.csv", row.names = FALSE)
trait_stats <- data.frame(
  Trait = colnames(traits_pred),
  Mean_Score = colMeans(traits_pred, na.rm = TRUE),
  Std_Score = apply(traits_pred, 2, sd, na.rm = TRUE),
  Min_Score = apply(traits_pred, 2, min, na.rm = TRUE),
  Max_Score = apply(traits_pred, 2, max, na.rm = TRUE),
  Positive_Samples = colSums(traits_pred > 0, na.rm = TRUE),
  Negative_Samples = colSums(traits_pred < 0, na.rm = TRUE)
)
write.csv(trait_stats, "4.Traits_association/trait_statistics.csv", row.names = FALSE)

create_summary_plot <- function(traits_df, trait_classes, output_path) {
  summary_data <- trait_classes %>%
    mutate(Risk_Level = ifelse(Trait_Risk_Score > 0, "High Risk", "Low Risk"),
           Num_High_Risk = sapply(strsplit(High_Risk_Traits, ","),  function(x) ifelse(x[1] == "", 0, length(x))),
           Num_Protective = sapply(strsplit(Protective_Traits, ","),  function(x) ifelse(x[1] == "", 0, length(x))))
  p1 <- ggplot(summary_data, aes(x = Num_High_Risk, fill = Risk_Level)) +
    geom_histogram(binwidth = 1, alpha = 0.7, position = "dodge") +
    labs(title = "Distribution of High-Risk Traits",x = "Number of High-Risk Traits",y = "Count") +
    theme_minimal() +
    scale_fill_manual(values = c("High Risk" = "red", "Low Risk" = "green"))
  p2 <- ggplot(summary_data, aes(x = Trait_Risk_Score, fill = Risk_Level)) +
    geom_density(alpha = 0.5) +
    labs(title = "Distribution of Trait Risk Scores",x = "Risk Score", y = "Density") +
    theme_minimal() +
    scale_fill_manual(values = c("High Risk" = "red", "Low Risk" = "green"))
  p3 <- ggplot(summary_data, aes(x = Num_High_Risk, y = Num_Protective, color = Risk_Level)) +
    geom_point(alpha = 0.6, size = 3) +
    labs(title = "High-Risk vs Protective Traits",x = "Number of High-Risk Traits",y = "Number of Protective Traits") +
    theme_minimal() +
    scale_color_manual(values = c("High Risk" = "red", "Low Risk" = "green"))
  summary_plot <- ggarrange(p1, p2, p3, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
  summary_plot <- annotate_figure(summary_plot, top = text_grob("Summary of Trait Predictions", face = "bold", size = 16))
  ggsave(output_path, summary_plot, width = 14, height = 10, dpi = 300)
  cat("  汇总图已保存到:", output_path, "\n")
}
create_summary_plot(traits_pred, trait_classes, "4.Traits_association/summary_plot.png")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("分析完成！\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("所有结果已保存到 4.Traits_association/ 目录\n\n")
cat("生成的主要文件：\n")
cat("1. 4.Traits_association/traits_prediction_raw.csv - 原始性状预测得分\n")
cat("2. 4.Traits_association/trait_classification.csv - 性状分类结果\n")
cat("3. 4.Traits_association/top_traits_per_sample.csv - 每个样本的前5个预测性状\n")
cat("4. 4.Traits_association/trait_statistics.csv - 性状统计摘要\n")
cat("5. 4.Traits_association/analysis_report.txt - 分析报告\n\n")
cat("生成的可视化图表：\n")
cat("1. 4.Traits_association/traits_heatmap.png - 性状预测热图\n")
cat("2. 4.Traits_association/trait_distribution.png - 性状得分分布\n")
cat("3. 4.Traits_association/fibroblast_correlation.png - 成纤维细胞相关性\n")
cat("4. 4.Traits_association/risk_scores.png - 样本风险评分\n")
cat("5. 4.Traits_association/trait_correlation_network.png - 性状相关性网络\n")
cat("6. 4.Traits_association/summary_plot.png - 分析汇总图\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
################################## END ##################################
