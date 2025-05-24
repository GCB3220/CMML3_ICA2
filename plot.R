# This is an R script used to plot a comparison graph of ARI, NMI, silhouette coefficient and runtime.
# Execution environment can be found in the sessioninfo at the end of the script.


# packages
library(ggplot2)
library(dplyr)

# plot function
plot_metric_comparison <- function(seurat_data, citefuse_data, metric_name) {
  df_summary <- data.frame(
    group = rep(c("Seurat", "CiteFuse"), each = length(seurat_data)),
    value = c(seurat_data, citefuse_data)
  ) %>%
    group_by(group) %>%
    summarise(
      mean = mean(value),
      sd = sd(value),
      .groups = "drop"
    )
  
  title_text <- paste("Mean", metric_name, "± Standard Deviation")
  y_label <- paste(metric_name, "Score")
  
  ggplot(df_summary, aes(x = group, y = mean, fill = group)) +
    geom_col(width = 0.6, alpha = 0.8) +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), 
                  width = 0.2, linewidth = 0.6) +
    labs(title = title_text, 
         y = y_label, 
         x = "") + 
    ylim(0, 1) +
    theme_minimal() +
    scale_y_continuous(limits = c(0, max(df_summary$mean + df_summary$sd) * 1.1),
                       expand = expansion(mult = c(0, 0))) +  # 完全从0开始
    theme(
      panel.grid = element_blank(),  # 简写取消所有网格线
      axis.title.y = element_text(face = "bold", size = 12),  # Y轴标签加粗
      axis.title.x = element_text(face = "bold", size = 12),  # X轴标签加粗
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),  # 加粗图例标题
      legend.text = element_text(face = "bold", size = 10),   # 加粗图例文字
      plot.title = element_text(face = "bold", hjust = 0.5)   # 加粗标题
    )
}

# load
seurat_data <- read.csv('combined_seurat_results_wide.csv')
citefuse_data <- read.csv('combined_citefuse_results_wide.csv')

# ARI
ari_seurat <- seurat_data$ARI
ari_citefuse <- citefuse_data$ARI
p_ari <- plot_metric_comparison(ari_seurat, ari_citefuse, 'ARI')
# NMI
nmi_seurat <- seurat_data$NMI
nmi_citefuse <- citefuse_data$NMI
p_nmi <- plot_metric_comparison(nmi_seurat, nmi_citefuse, 'NMI')
# Silhouette
sil_seurat <- seurat_data$Silhouette
sil_citefuse <- citefuse_data$Silhouette
p_sil <- plot_metric_comparison(sil_seurat, sil_citefuse, 'Silhouette')
# time
time_seurat <- seurat_data$Time_min
time_citefuse <- citefuse_data$Time_min
p_time <- plot_metric_comparison(time_seurat, time_citefuse, 'Time(min)')

save.image('plot.RData')


# > sessionInfo()
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26100)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
# [2] LC_CTYPE=Chinese (Simplified)_China.utf8   
# [3] LC_MONETARY=Chinese (Simplified)_China.utf8
# [4] LC_NUMERIC=C                               
# [5] LC_TIME=Chinese (Simplified)_China.utf8    
# 
# time zone: Asia/Shanghai
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods  
# [7] base     
# 
# other attached packages:
#   [1] dplyr_1.1.4   ggplot2_3.5.2
# 
# loaded via a namespace (and not attached):
#   [1] labeling_0.4.3     RColorBrewer_1.1-3 R6_2.6.1          
# [4] tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.3    
# [7] gtable_0.3.6       glue_1.7.0         tibble_3.2.1      
# [10] dichromat_2.0-0.1  pkgconfig_2.0.3    generics_0.1.4    
# [13] lifecycle_1.0.4    cli_3.6.3          scales_1.4.0      
# [16] grid_4.4.1         vctrs_0.6.5        withr_3.0.2       
# [19] compiler_4.4.1     rstudioapi_0.17.1  tools_4.4.1       
# [22] pillar_1.10.2      crayon_1.5.3       rlang_1.1.4    