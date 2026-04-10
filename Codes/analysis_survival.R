################################## START ##################################
# analysis_survival.R
# 输入: 细胞类型注释，临床生存信息
# 输出: 5.Survival_prediction文件夹及其所有分析结果
required_packages <- c("tidyverse", "survival", "survminer", "ggplot2", "pheatmap", "ggpubr", "forestplot", "reshape2", "RColorBrewer", "viridis", "cowplot", "gridExtra")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("安装缺失的包:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(pheatmap)
  library(ggpubr)
  library(forestplot)
  library(reshape2)
  library(RColorBrewer)
  library(viridis)
  library(cowplot)
  library(gridExtra)
})

dir.create("5.Survival_prediction", showWarnings = FALSE)
new_samples_df<-read.table("2.Identify_celltype/Example_annotations.csv",header=T,sep=",",row.names=1)
new_samples_df<-new_samples_df[,-c(1:3)]
new_samples_df<-t(new_samples_df[,1:16])
if(file.exists(unz("Reference.zip", "fibroblast_TCGA_samples.txt"))) {
  tcga_fibroblast <- read.table(unz("Reference.zip", "fibroblast_TCGA_samples.txt"),header=T,sep="\t")
  tcga_fibroblast_df <- as.data.frame(tcga_fibroblast)
  rownames(tcga_fibroblast_df) <- tcga_fibroblast_df[,1]
  tcga_fibroblast_df <- tcga_fibroblast_df[,-1]
  tcga_fibroblast_df <- t(tcga_fibroblast_df)
  cat("  TCGA成纤维细胞数据: ", nrow(tcga_fibroblast_df), "行 ×", ncol(tcga_fibroblast_df), "个样本\n")
} else {
  stop("错误: 找不到 fibroblast_TCGA_samples.txt 文件")
}
if(file.exists(unz("Reference.zip", "TCGA_clinical_data.txt"))) {
  tcga_clinical <- read.table(unz("Reference.zip", "TCGA_clinical_data.txt"),header=T,sep="\t")
  cat("  TCGA临床数据: ", nrow(tcga_clinical), "个样本 ×", ncol(tcga_clinical), "个变量\n")
} else {
  stop("错误: 找不到 TCGA_clinical_data.txt 文件")
}
cat("\n2. 进行数据预处理...\n")
tcga_samples <- colnames(tcga_fibroblast_df)
cat("  TCGA样本数量: ", length(tcga_samples), "\n")
fibroblast_prop_tcga <- as.numeric(tcga_fibroblast_df[1, ])
names(fibroblast_prop_tcga) <- colnames(tcga_fibroblast_df)
if(!is.null(new_samples_df)) {
  new_samples_names <- colnames(new_samples_df)
  cat("  新样本数量: ", length(new_samples_names), "\n")
  fibroblast_prop_new <- as.numeric(new_samples_df[1, ])
  names(fibroblast_prop_new) <- colnames(new_samples_df)
} else {
  fibroblast_prop_new <- NULL
}
calculate_simple_roc <- function(time, status, marker, time_point) {
  event_at_time <- ifelse(time <= time_point & status == 1, 1, 0)
  valid_idx <- !is.na(event_at_time) & !is.na(marker)
  event_at_time <- event_at_time[valid_idx]
  marker <- marker[valid_idx]
  if(length(unique(event_at_time)) < 2) {
    return(NA) 
  }
  n <- length(event_at_time)
  n_cases <- sum(event_at_time == 1)
  n_controls <- sum(event_at_time == 0)
  
  if(n_cases == 0 | n_controls == 0) {
    return(NA)
  }
  scores_cases <- marker[event_at_time == 1]
  scores_controls <- marker[event_at_time == 0]
  # 计算AUC
  auc <- 0
  for(i in 1:n_cases) {
    auc <- auc + sum(scores_cases[i] > scores_controls)
  }
  auc <- auc / (n_cases * n_controls)
  return(auc)
}
analysis_survival <- function(tcga_fibroblast_data, tcga_clinical_data, new_samples_data = NULL, cancer_type = "all", survival_time = "OS.time", survival_status = "OS") {
  cat(sprintf("\n对癌症类型 '%s' 进行生存分析...\n", cancer_type))
  results <- list()
  # 4.1 准备TCGA数据
  if(cancer_type != "all") {
    if("cancer_type" %in% colnames(tcga_clinical_data)) {
      tcga_clinical_subset <- tcga_clinical_data %>% 
        filter(cancer_type == !!cancer_type)
    } else {
      warning("临床数据中没有cancer_type列，将使用所有数据")
      tcga_clinical_subset <- tcga_clinical_data
    }
  } else {
    tcga_clinical_subset <- tcga_clinical_data
  }
  if(!"sample" %in% colnames(tcga_clinical_subset)) {
    sample_id_col <- colnames(tcga_clinical_subset)[1]
    colnames(tcga_clinical_subset)[colnames(tcga_clinical_subset) == sample_id_col] <- "sample"
  }
  analysis_data <- data.frame(
    sample = tcga_clinical_subset$sample,
    fibroblast_proportion = NA,
    time = as.numeric(tcga_clinical_subset[[survival_time]]),
    status = as.numeric(tcga_clinical_subset[[survival_status]])
  )
  for(i in 1:nrow(analysis_data)) {
    sample_id <- analysis_data$sample[i]
    if(sample_id %in% names(tcga_fibroblast_data)) {
      analysis_data$fibroblast_proportion[i] <- tcga_fibroblast_data[sample_id]
    }
  }
  analysis_data <- analysis_data[complete.cases(analysis_data), ]
  if(nrow(analysis_data) < 10) {
    warning(sprintf("癌症类型 %s 的样本数太少 (%d)，跳过分析", cancer_type, nrow(analysis_data)))
    return(NULL)
  }
  cat(sprintf("  可用样本数: %d\n", nrow(analysis_data)))
  cat(sprintf("  事件数: %d\n", sum(analysis_data$status)))
  fibroblast_median <- median(analysis_data$fibroblast_proportion, na.rm = TRUE)
  analysis_data$fibroblast_group <- ifelse(analysis_data$fibroblast_proportion > fibroblast_median, "High", "Low")
  analysis_data$fibroblast_group <- factor(analysis_data$fibroblast_group, levels = c("Low", "High"))
  surv_obj <- Surv(time = analysis_data$time,  event = analysis_data$status)
  tryCatch({
    cox_model <- coxph(surv_obj ~ fibroblast_proportion, data = analysis_data)
    cox_summary <- summary(cox_model)
    hr <- exp(cox_summary$coefficients[1, 1])
    hr_lower <- exp(cox_summary$coefficients[1, 1] - 1.96 * cox_summary$coefficients[1, 3])
    hr_upper <- exp(cox_summary$coefficients[1, 1] + 1.96 * cox_summary$coefficients[1, 3])
    p_value <- cox_summary$coefficients[1, 5]
    results$cox_model <- cox_model
    results$cox_summary <- cox_summary
    results$hr <- hr
    results$hr_ci <- c(hr_lower, hr_upper)
    results$p_value <- p_value
    cat(sprintf("  Cox模型结果: HR = %.2f (%.2f-%.2f), P = %.4f\n", hr, hr_lower, hr_upper, p_value))
  }, error = function(e) {
    cat(sprintf("  Cox模型拟合失败: %s\n", e$message))
    results$cox_model <- NULL
    results$hr <- NA
    results$hr_ci <- c(NA, NA)
    results$p_value <- NA
  })
  tryCatch({
    logrank_test <- survdiff(surv_obj ~ fibroblast_group, data = analysis_data)
    logrank_p <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
    results$logrank_p <- logrank_p
    cat(sprintf("  Log-rank检验: P = %.4f\n", logrank_p))
  }, error = function(e) {
    cat(sprintf("  Log-rank检验失败: %s\n", e$message))
    results$logrank_p <- NA
  })
  if(nrow(analysis_data) >= 20) {
    tryCatch({
      median_time <- median(analysis_data$time[analysis_data$status == 1], na.rm = TRUE)
      if(!is.na(median_time) && median_time > 0) {
        auc <- calculate_simple_roc(analysis_data$time, analysis_data$status, analysis_data$fibroblast_proportion,median_time)
        results$auc_median <- auc
        cat(sprintf("  中位时间点的AUC: %.3f\n", auc))
      }
    }, error = function(e) {
      results$auc_median <- NA
    })
  }
  if(!is.null(new_samples_data)) {
    new_samples_risk <- data.frame(
      sample = names(new_samples_data),
      fibroblast_proportion = new_samples_data,
      risk_group = ifelse(new_samples_data > fibroblast_median, "High Risk", "Low Risk"),
      risk_score = NA
    )
    if(!is.null(results$cox_model) && !is.na(results$hr)) {
      for(i in 1:nrow(new_samples_risk)) {
        new_samples_risk$risk_score[i] <- exp(log(results$hr) * new_samples_risk$fibroblast_proportion[i])
      }
    }
    results$new_samples_risk <- new_samples_risk
  }
  results$analysis_data <- analysis_data
  results$fibroblast_median <- fibroblast_median
  results$cancer_type <- cancer_type
  results$n_samples <- nrow(analysis_data)
  results$n_events <- sum(analysis_data$status)
  results$surv_obj <- surv_obj
  return(results)
}
cat("\n3. 对所有癌症类型进行生存分析...\n")
# 获取所有癌症类型
if("cancer_type" %in% colnames(tcga_clinical)) {
  cancer_types <- unique(tcga_clinical$cancer_type)
  cat("  检测到", length(cancer_types), "种癌症类型\n")
  if(length(cancer_types) > 33) {
    cat("  注意: 检测到", length(cancer_types), "种癌症类型，将分析前33种\n")
    cancer_types <- cancer_types[1:33]
  }
} else {
  cancer_types <- "All_TCGA"
  tcga_clinical$cancer_type <- "All_TCGA"
}
all_results <- list()
summary_results <- data.frame()
for(cancer in cancer_types) {
  cat(sprintf("  分析 %s... ", cancer))
  result <- analysis_survival(tcga_fibroblast_data = fibroblast_prop_tcga,tcga_clinical_data = tcga_clinical,new_samples_data = fibroblast_prop_new,cancer_type = cancer)
  if(!is.null(result)) {
    all_results[[cancer]] <- result
    summary_results <- rbind(summary_results, data.frame(
      Cancer_Type = cancer,
      N_Samples = result$n_samples,
      N_Events = result$n_events,
      HR = ifelse(is.null(result$hr), NA, result$hr),
      HR_Lower = ifelse(is.null(result$hr_ci), NA, result$hr_ci[1]),
      HR_Upper = ifelse(is.null(result$hr_ci), NA, result$hr_ci[2]),
      P_Value = ifelse(is.null(result$p_value), NA, result$p_value),
      LogRank_P = ifelse(is.null(result$logrank_p), NA, result$logrank_p),
      AUC_Median = ifelse(is.null(result$auc_median), NA, result$auc_median),
      Fibroblast_Median = result$fibroblast_median,
      Significant = ifelse(!is.na(result$p_value) && result$p_value < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    ))
    cat("完成\n")
  } else {
    cat("跳过（样本不足）\n")
  }
}

create_survival_plots <- function(results_list, output_dir) {
  for(cancer in names(results_list)) {
    result <- results_list[[cancer]]    
    if(is.null(result)) next
    if(is.null(result$analysis_data)) next 
    # 准备数据
    surv_data <- result$analysis_data
    surv_data$group <- surv_data$fibroblast_group
    if(length(unique(surv_data$group)) < 2) next
    tryCatch({
      fit <- survfit(Surv(time, status) ~ group, data = surv_data)
      p <- ggsurvplot(fit,data = surv_data,pval = TRUE,pval.method = TRUE,conf.int = FALSE,risk.table = TRUE,risk.table.height = 0.25,ggtheme = theme_minimal(),palette = c("#2E9FDF", "#E7B800"),title = paste("Survival Analysis -", cancer),xlab = "Time (days)",ylab = "Survival Probability",legend = "none",legend.title = "Fibroblast",legend.labs = c("Low", "High"),font.title = 16,font.x = 14,font.y = 14,font.tickslab = 12)
      plt<-arrange_ggsurvplots(list(p),nrow=1,ncol=1)
      dev.off()
      dev.off()
      output_file <- file.path(output_dir, paste0("survival_curve_", gsub("[^A-Za-z0-9]", "_", cancer), ".png"))
      ggsave(output_file, plt, width = 10, height = 8)
    }, error = function(e) {
      cat(sprintf("无法为 %s 生成生存曲线: %s\n", cancer, e$message))
    })
  }
  cat("  生存曲线图已保存\n")
}
create_forest_plot <- function(summary_df, output_path) {
  summary_df_clean <- summary_df[!is.na(summary_df$HR) & !is.na(summary_df$P_Value), ]
  if(nrow(summary_df_clean) == 0) {
    cat("  没有有效的HR数据用于生成森林图\n")
    return(NULL)
  }
  summary_df_clean <- summary_df_clean[order(summary_df_clean$HR), ]
  tabletext <- cbind(
    c("Cancer Type", as.character(summary_df_clean$Cancer_Type)),
    c("N", sprintf("%d", summary_df_clean$N_Samples)),
    c("Events", sprintf("%d", summary_df_clean$N_Events)),
    c("HR (95% CI)", sprintf("%.2f (%.2f-%.2f)", summary_df_clean$HR, summary_df_clean$HR_Lower, summary_df_clean$HR_Upper)),
    c("P Value", sprintf("%.3f", summary_df_clean$P_Value))
  )
  mean <- c(NA, summary_df_clean$HR)
  lower <- c(NA, summary_df_clean$HR_Lower)
  upper <- c(NA, summary_df_clean$HR_Upper)
  tryCatch({
    png(output_path, width = 1200, height=1200)
    forestplot(labeltext = tabletext,mean = mean,lower = lower,upper = upper,is.summary = c(TRUE, rep(FALSE, nrow(summary_df_clean))),xlab = "Hazard Ratio",zero = 1,boxsize = 0.3,col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),txt_gp = fpTxtGp(label = gpar(cex = 0.8),ticks = gpar(cex = 0.8),xlab = gpar(cex = 1)))
    dev.off()
    cat("  森林图已保存到:", output_path, "\n")
  }, error = function(e) {
    cat("  无法生成森林图:", e$message, "\n")
  })
}
create_risk_prediction_plot <- function(all_results, output_path) {
  if(is.null(fibroblast_prop_new)) {
    cat("  没有新样本数据，跳过风险预测图\n")
    return(NULL)
  }
  all_new_predictions <- data.frame()
  for(cancer in names(all_results)) {
    result <- all_results[[cancer]]
    if(!is.null(result) && !is.null(result$new_samples_risk)) {
      temp_df <- result$new_samples_risk
      temp_df$Cancer_Type <- cancer
      temp_df$HR <- ifelse(is.null(result$hr), NA, result$hr)
      temp_df$P_Value <- ifelse(is.null(result$p_value), NA, result$p_value)
      all_new_predictions <- rbind(all_new_predictions, temp_df)
    }
  }
  if(nrow(all_new_predictions) == 0) {
    cat("  没有可用的新样本预测数据\n")
    return(NULL)
  }
  tryCatch({
    heatmap_data <- all_new_predictions %>%
      select(sample, Cancer_Type, risk_score) %>%
      pivot_wider(names_from = Cancer_Type, values_from = risk_score) %>%
      column_to_rownames("sample")
    heatmap_data[is.na(heatmap_data)] <- colMeans(heatmap_data, na.rm = TRUE)
    heatmap_data_scaled <- as.matrix(scale(heatmap_data))
    png(output_path, width = 1000, height = 800)
    pheatmap(heatmap_data_scaled,main = "New Samples Risk Prediction Across Cancer Types",color = colorRampPalette(c("blue", "white", "red"))(100),cluster_rows = TRUE,cluster_cols = TRUE,show_rownames = ifelse(nrow(heatmap_data_scaled) <= 50, TRUE, FALSE),show_colnames = TRUE,fontsize_row = 10,fontsize_col = 10,border_color = NA)
    dev.off()
    cat("  新样本风险预测图已保存到:", output_path, "\n")
  }, error = function(e) {
    cat("  无法生成风险预测热图:", e$message, "\n")
  })
}
create_hr_distribution_plot <- function(summary_df, output_path) {
  plot_data <- summary_df[!is.na(summary_df$HR) & !is.na(summary_df$P_Value), ]
  if(nrow(plot_data) == 0) {
    cat("  没有有效的HR数据用于生成分布图\n")
    return(NULL)
  }
  tryCatch({
    p <- ggplot(plot_data, aes(x = reorder(Cancer_Type, HR), y = HR, color = Significant, fill = Significant)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = HR_Lower, ymax = HR_Upper), width = 0.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.5) +
      scale_color_manual(values = c("Yes" = "red", "No" = "blue")) +
      scale_fill_manual(values = c("Yes" = "red", "No" = "blue")) +
      coord_flip() +
      labs(title = "Hazard Ratios Across Cancer Types",x = "Cancer Type",y = "Hazard Ratio (95% CI)",color = "Significant (P < 0.05)",fill = "Significant (P < 0.05)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),axis.text.y = element_text(size = 10),legend.position = "bottom")
    ggsave(output_path, p, width = 12, height = max(8, nrow(plot_data) * 0.3), dpi = 300)
    cat("  HR分布图已保存到:", output_path, "\n")
  }, error = function(e) {
    cat("  无法生成HR分布图:", e$message, "\n")
  })
}
create_new_samples_risk_barplot <- function(all_results, output_path) {
  if(is.null(fibroblast_prop_new)) {
    cat("  没有新样本数据，跳过风险条形图\n")
    return(NULL)
  }
  valid_cancers <- names(all_results)[sapply(all_results, function(x) !is.null(x) && !is.null(x$new_samples_risk))]
  if(length(valid_cancers) == 0) {
    cat("  没有可用的新样本风险数据\n")
    return(NULL)
  }
  cancer <- valid_cancers[1]
  result <- all_results[[cancer]]
  tryCatch({
    new_risk <- result$new_samples_risk
    new_risk <- new_risk[order(new_risk$fibroblast_proportion, decreasing = TRUE), ]
    new_risk$sample <- factor(new_risk$sample, levels = new_risk$sample)
    p <- ggplot(new_risk, aes(x = sample, y = fibroblast_proportion, fill = risk_group)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = result$fibroblast_median, linetype = "dashed", color = "red", linewidth = 1) +
      annotate("text", x = nrow(new_risk)/2, y = result$fibroblast_median * 1.1,label = paste("Median =", round(result$fibroblast_median, 3)),color = "red", size = 4) +
      scale_fill_manual(values = c("High Risk" = "red", "Low Risk" = "blue")) +
      labs(title = paste("New Samples Fibroblast Proportion and Risk Prediction"),subtitle = paste("Cancer Type:", result$cancer_type, ifelse(!is.na(result$hr), paste("| HR =", round(result$hr, 2), "| P =", sprintf("%.3f", result$p_value)),"")),x = "Sample",y = "Fibroblast Proportion",fill = "Risk Group") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),plot.title = element_text(size = 16, face = "bold"),plot.subtitle = element_text(size = 12),legend.position = "top")
    ggsave(output_path, p, width = min(20, nrow(new_risk) * 0.5), height = 8, dpi = 300)
    cat("  新样本风险条形图已保存到:", output_path, "\n")
  }, error = function(e) {
    cat("  无法生成新样本风险条形图:", e$message, "\n")
  })
}
create_auc_plot <- function(summary_df, output_path) {
  plot_data <- summary_df[!is.na(summary_df$AUC_Median), ]
  if(nrow(plot_data) == 0) {
    cat("  没有AUC数据用于生成性能图\n")
    return(NULL)
  }
  tryCatch({
    p <- ggplot(plot_data, aes(x = reorder(Cancer_Type, AUC_Median), y = AUC_Median, fill = AUC_Median > 0.5)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5) +
      scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "orange"),name = "AUC > 0.5") +
      coord_flip() +
      labs(title = "AUC at Median Survival Time Across Cancer Types",x = "Cancer Type",y = "AUC") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),axis.text.y = element_text(size = 10),legend.position = "bottom")
    ggsave(output_path, p, width = 10, height = max(6, nrow(plot_data) * 0.3), dpi = 300)
    cat("  AUC性能图已保存到:", output_path, "\n")
  }, error = function(e) {
    cat("  无法生成AUC性能图:", e$message, "\n")
  })
}
create_summary_plot <- function(summary_df, output_path) {
  plot_data <- summary_df[!is.na(summary_df$HR) & !is.na(summary_df$P_Value), ]
  if(nrow(plot_data) == 0) {
    cat("  没有有效数据用于生成汇总图\n")
    return(NULL)
  }
  tryCatch({
    p1 <- ggplot(plot_data, aes(x = N_Samples)) +
      geom_histogram(fill = "steelblue", color = "white", bins = 20) +
      labs(title = "Distribution of Sample Sizes",x = "Number of Samples",y = "Count") +theme_minimal()
    p2 <- ggplot(plot_data, aes(x = HR)) +
      geom_histogram(fill = "coral", color = "white", bins = 20) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
      labs(title = "Distribution of Hazard Ratios",x = "Hazard Ratio", y = "Count") +
      theme_minimal()
    p3 <- ggplot(plot_data, aes(x = reorder(Cancer_Type, -log10(P_Value)), y = -log10(P_Value))) +
      geom_bar(stat = "identity", fill = ifelse(plot_data$P_Value < 0.05, "red", "gray")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      labs(title = "Significance of Associations",x = "Cancer Type",y = "-log10(P Value)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # 组合图
    summary_plot <- plot_grid(p1, p2, p3, ncol = 2, nrow = 2)
    title <- ggdraw() + 
      draw_label("Summary of Fibroblast Survival Analysis", fontface = 'bold', size = 20)
    final_plot <- plot_grid(title, summary_plot, ncol = 1, rel_heights = c(0.1, 1))
    ggsave(output_path, final_plot, width = 16, height = 12, dpi = 300)
    cat("  汇总图已保存到:", output_path, "\n")
  }, error = function(e) {
    cat("  无法生成汇总图:", e$message, "\n")
  })
}

create_survival_plots(all_results, "5.Survival_prediction/visualizations")
if(nrow(summary_results) > 0) {
  create_forest_plot(summary_results, "5.Survival_prediction/")
  create_hr_distribution_plot(summary_results,"5.Survival_prediction/hr_distribution.png")
  create_auc_plot(summary_results,"5.Survival_prediction/auc_performance.png")
  create_summary_plot(summary_results,"5.Survival_prediction/summary_plot.png")
  write_csv(summary_results, "5.Survival_prediction/summary_results.csv")
  cat("摘要结果已保存: 5.Survival_prediction/summary_results.csv\n")
}
create_risk_prediction_plot(all_results,"5.Survival_prediction/new_samples_risk_heatmap.png")
create_new_samples_risk_barplot(all_results,"5.Survival_prediction/new_samples_risk_barplot.png")
cat("\n6. 生成综合报告...\n")
generate_comprehensive_report <- function(all_results, summary_df, output_path) {
  report_lines <- c()
  report_lines <- c(report_lines, paste(rep("=", 70), collapse = ""))
  report_lines <- c(report_lines, "FIBROBLAST SURVIVAL ANALYSIS COMPREHENSIVE REPORT")
  report_lines <- c(report_lines, paste(rep("=", 70), collapse = ""))
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, sprintf("Analysis Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "1. ANALYSIS OVERVIEW")
  report_lines <- c(report_lines, paste(rep("-", 50), collapse = ""))
  report_lines <- c(report_lines, sprintf("Total cancer types analyzed: %d", length(all_results)))
  report_lines <- c(report_lines, sprintf("Total TCGA samples: %d", length(fibroblast_prop_tcga)))
  if(!is.null(fibroblast_prop_new)) {
    report_lines <- c(report_lines, sprintf("Total new samples: %d", length(fibroblast_prop_new)))
  }
  valid_results <- summary_df[!is.na(summary_df$HR) & !is.na(summary_df$P_Value), ]
  report_lines <- c(report_lines, sprintf("Valid cancer types with HR results: %d", nrow(valid_results)))
  report_lines <- c(report_lines, "")
  significant_cancers <- valid_results %>% filter(Significant == "Yes")
  report_lines <- c(report_lines, "2. SIGNIFICANT FINDINGS")
  report_lines <- c(report_lines, paste(rep("-", 50), collapse = ""))
  report_lines <- c(report_lines, sprintf("Number of cancer types with significant association (P < 0.05): %d", nrow(significant_cancers)))
  if(nrow(significant_cancers) > 0) {
    report_lines <- c(report_lines, "")
    report_lines <- c(report_lines, "Top 5 most significant associations:")
    top5 <- significant_cancers %>% arrange(P_Value) %>% head(5)
    for(i in 1:min(5, nrow(top5))) {
      report_lines <- c(report_lines, sprintf("  %d. %s: HR = %.2f (%.2f-%.2f), P = %.4f", i, 
        top5$Cancer_Type[i],top5$HR[i],top5$HR_Lower[i],top5$HR_Upper[i],top5$P_Value[i]))
    }
  }
  report_lines <- c(report_lines, "")
  high_risk_cancers <- valid_results %>% filter(HR > 1.5 & Significant == "Yes")
  if(nrow(high_risk_cancers) > 0) {
    report_lines <- c(report_lines, "3. HIGH-RISK CANCER TYPES (HR > 1.5)")
    report_lines <- c(report_lines, paste(rep("-", 50), collapse = ""))
    for(i in 1:nrow(high_risk_cancers)) {
      report_lines <- c(report_lines,sprintf("  %s: High fibroblast proportion associated with poor prognosis (HR = %.2f)",high_risk_cancers$Cancer_Type[i],high_risk_cancers$HR[i]))
    }
    report_lines <- c(report_lines, "")
  }
  protective_cancers <- valid_results %>% filter(HR < 0.67 & Significant == "Yes")
  if(nrow(protective_cancers) > 0) {
    report_lines <- c(report_lines, "4. PROTECTIVE CANCER TYPES (HR < 0.67)")
    report_lines <- c(report_lines, paste(rep("-", 50), collapse = ""))
    for(i in 1:nrow(protective_cancers)) {
      report_lines <- c(report_lines,sprintf("  %s: High fibroblast proportion associated with better prognosis (HR = %.2f)",protective_cancers$Cancer_Type[i],protective_cancers$HR[i]))
    }
    report_lines <- c(report_lines, "")
  }
  auc_cancers <- valid_results[!is.na(valid_results$AUC_Median), ]
  if(nrow(auc_cancers) > 0) {
    report_lines <- c(report_lines, "5. PREDICTIVE PERFORMANCE (AUC)")
    report_lines <- c(report_lines, paste(rep("-", 50), collapse = ""))
    high_perf <- auc_cancers %>% filter(AUC_Median > 0.7) %>% arrange(desc(AUC_Median))
    if(nrow(high_perf) > 0) {
      report_lines <- c(report_lines, "Cancer types with good predictive performance (AUC > 0.7):")
      for(i in 1:min(5, nrow(high_perf))) {
        report_lines <- c(report_lines,sprintf("  %s: AUC = %.3f", high_perf$Cancer_Type[i],high_perf$AUC_Median[i]))
      }
      report_lines <- c(report_lines, "")
    }
  }
  if(!is.null(fibroblast_prop_new)) {
    report_lines <- c(report_lines, "6. NEW SAMPLES PREDICTION")
    report_lines <- c(report_lines, paste(rep("-", 50), collapse = ""))
    if(length(all_results) > 0) {
      first_cancer <- names(all_results)[1]
      result <- all_results[[first_cancer]]
      if(!is.null(result$new_samples_risk)) {
        new_risk <- result$new_samples_risk
        high_risk_samples <- new_risk %>% filter(risk_group == "High Risk")
        low_risk_samples <- new_risk %>% filter(risk_group == "Low Risk")
        report_lines <- c(report_lines, sprintf("Based on %s model:", first_cancer))
        report_lines <- c(report_lines, sprintf("  High-risk samples: %d (%.1f%%)", nrow(high_risk_samples), nrow(high_risk_samples)/nrow(new_risk)*100))
        report_lines <- c(report_lines, sprintf("  Low-risk samples: %d (%.1f%%)", nrow(low_risk_samples),nrow(low_risk_samples)/nrow(new_risk)*100))
        if(nrow(high_risk_samples) > 0) {
          report_lines <- c(report_lines, "")
          report_lines <- c(report_lines, "  Top 5 highest risk samples:")
          top5_high <- high_risk_samples %>% 
            arrange(desc(fibroblast_proportion)) %>% 
            head(5)
          for(i in 1:nrow(top5_high)) {
            report_lines <- c(report_lines,sprintf("    %s: Proportion = %.3f, Risk Score = %.3f",top5_high$sample[i],top5_high$fibroblast_proportion[i],top5_high$risk_score[i]))
          }
        }
      }
    }
  }
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "7. OUTPUT FILES")
  report_lines <- c(report_lines, paste(rep("-", 50), collapse = ""))
  report_lines <- c(report_lines, "The following files have been generated:")
  report_lines <- c(report_lines, "  5.Survival_prediction/summary_results.csv - Summary statistics for all cancer types")
  report_lines <- c(report_lines, "  5.Survival_prediction/ - Directory containing all visualizations")
  report_lines <- c(report_lines, "  5.Survival_prediction/ - Directory containing R data objects")
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, "Visualization files include:")
  report_lines <- c(report_lines, "  - survival_curve_*.png: Survival curves for each cancer type")
  report_lines <- c(report_lines, "  - forest_plot.png: Forest plot of hazard ratios")
  report_lines <- c(report_lines, "  - hr_distribution.png: Distribution of HRs across cancer types")
  report_lines <- c(report_lines, "  - auc_performance.png: AUC performance across cancer types")
  report_lines <- c(report_lines, "  - summary_plot.png: Summary visualization")
  if(!is.null(fibroblast_prop_new)) {
    report_lines <- c(report_lines, "  - new_samples_risk_heatmap.png: Risk prediction heatmap for new samples")
    report_lines <- c(report_lines, "  - new_samples_risk_barplot.png: Bar plot of new samples risk prediction")
  }
  report_lines <- c(report_lines, "")
  report_lines <- c(report_lines, paste(rep("=", 70), collapse = ""))
  report_lines <- c(report_lines, "ANALYSIS COMPLETE")
  report_lines <- c(report_lines, paste(rep("=", 70), collapse = ""))
  writeLines(report_lines, con = output_path)
  cat(paste(report_lines, collapse = "\n"), "\n")
  return(report_lines)
}
report <- generate_comprehensive_report(all_results, summary_results,"5.Survival_prediction/analysis_report.txt")
tryCatch({
  saveRDS(all_results, "5.Survival_prediction/all_results.rds")
  saveRDS(summary_results, "5.Survival_prediction/summary_results.rds")
  cat("  模型和结果已保存到 models/ 目录\n")
}, error = function(e) {
  cat("  无法保存模型文件:", e$message, "\n")
})
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("生存分析完成！\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("所有结果已保存到 5.Survival_prediction/ 目录\n\n")
cat("生成的主要文件：\n")
cat("1. 5.Survival_prediction/summary_results.csv - 各癌症类型摘要统计\n")
cat("2. 5.Survival_prediction/analysis_report.txt - 详细分析报告\n")
cat("3. 5.Survival_prediction/ - 各癌症类型分析模型\n\n")
cat("生成的可视化图表：\n")
cat("1. survival_curve_*.png - 各癌症类型的生存曲线\n")
cat("2. forest_plot.png - HR森林图\n")
cat("3. hr_distribution.png - HR分布图\n")
cat("4. auc_performance.png - AUC性能图\n")
cat("5. summary_plot.png - 汇总分析图\n")
if(!is.null(fibroblast_prop_new)) {
  cat("6. new_samples_risk_heatmap.png - 新样本风险热图\n")
  cat("7. new_samples_risk_barplot.png - 新样本风险条形图\n")
}
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("\n8. 重要发现总结：\n")
if(nrow(summary_results) > 0) {
  valid_results <- summary_results[!is.na(summary_results$HR) & !is.na(summary_results$P_Value), ]
  if(nrow(valid_results) > 0) {
    high_risk <- valid_results %>% 
      filter(Significant == "Yes" & HR > 1.3) %>%
      arrange(desc(HR))
    if(nrow(high_risk) > 0) {
      cat("\n高风险癌症（HR > 1.3）：\n")
      for(i in 1:min(5, nrow(high_risk))) {
        cat(sprintf("  %s: HR = %.2f, P = %.4f\n", high_risk$Cancer_Type[i],high_risk$HR[i],high_risk$P_Value[i]))
      }
    }
    protective <- valid_results %>% 
      filter(Significant == "Yes" & HR < 0.77) %>%
      arrange(HR)
    if(nrow(protective) > 0) {
      cat("\n保护性癌症（HR < 0.77）：\n")
      for(i in 1:min(5, nrow(protective))) {
        cat(sprintf("  %s: HR = %.2f, P = %.4f\n", protective$Cancer_Type[i],protective$HR[i],protective$P_Value[i]))
      }
    }
    most_sig <- valid_results %>% 
      filter(Significant == "Yes") %>%
      arrange(P_Value) %>%
      head(3)
    if(nrow(most_sig) > 0) {
      cat("\n最显著的关联：\n")
      for(i in 1:nrow(most_sig)) {
        cat(sprintf("  %s: HR = %.2f, P = %.4f\n", most_sig$Cancer_Type[i],most_sig$HR[i],most_sig$P_Value[i]))
      }
    }
  } else {
    cat("  没有有效的HR结果\n")
  }
} else {
  cat("  没有可用的分析结果\n")
}
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("分析完成！\n")
################################## END ##################################