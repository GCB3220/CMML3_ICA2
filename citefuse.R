# This R script integrates CITE-seq data using CiteFuse.
# The specific process and parameter settings here refer to: https://bioconductor.org/packages/release/bioc/vignettes/CiteFuse/inst/doc/CiteFuse.html  
# All parameters remained unchanged.
# The dataset used is the bmcite dataset from SeuratData. 
# First, the data is randomly sampled to 4000 cells. Then, data preprocessing is performed, followed by integration of RNA and ADT data. Finally, cell clustering is carried out.
# The clustering performance is evaluated using ARI, NMI and Silhouette Score.
# The above process will loop 5 times
# The performance of data integration is measured by taking the average of three indicators.
# All data were retained in the end
# Execution environment can be found in the sessioninfo at the end of the script.


# packages
library(CiteFuse)
library(scater)
library(SingleCellExperiment)
library(DT)
library(patchwork)
library(Seurat)
library(SeuratData)
library(mclust)      # ARI
library(aricode)     # NMI
library(cluster)     # silhouette

# random seeds
random_seeds <- c(123, 234, 345, 456, 567)
# load
data("bmcite")

# loop
for(seed in random_seeds) {
  start_time <- Sys.time()
  # extract 4000 randomly
  set.seed(seed)
  random_cells <- sample(colnames(bmcite), size = 4000, replace = FALSE)
  bm <- subset(bmcite, cells = random_cells)
  true_labels <- bm$celltype.l2
  # preprocessing
  CITEseq <- list(RNA = bm@assays$RNA@counts, ADT = bm@assays$ADT@counts)
  sce_citeseq <- preprocessing(CITEseq)
  # normalize
  sce_citeseq <- scater::logNormCounts(sce_citeseq)
  sce_citeseq <- normaliseExprs(sce_citeseq, altExp_name = "ADT", transform = "log")
  # CiteFuse
  sce_citeseq <- CiteFuse(sce_citeseq)
  # t-SNE
  sce_citeseq <- reducedDimSNF(sce_citeseq,
                               method = "tSNE", 
                               dimNames = "tSNE_joint")
  # Louvain
  set.seed(2024)  
  SNF_W_louvain <- igraphClustering(sce_citeseq, method = "louvain")
  sce_citeseq$SNF_W_louvain <- as.factor(SNF_W_louvain)
  # plot
  sce_citeseq$true_labels <- true_labels
  p <- visualiseDim(sce_citeseq, 
                    dimNames = "tSNE_joint", 
                    colour_by = "true_labels") +
    labs(title = paste("celltype.l2 (Seed:", seed, ")")) +
    theme(legend.position = "right", 
          legend.text = element_text(size = 12))
  
  # benchmarkng
  cluster_labels <- sce_citeseq$SNF_W_louvain
  ari <- ARI(cluster_labels, true_labels)
  nmi <- NMI(cluster_labels, true_labels)
  snf_embedding <- reducedDim(sce_citeseq, "tSNE_joint")
  sil_dist <- dist(snf_embedding) 
  sil <- silhouette(as.numeric(cluster_labels), sil_dist)
  sil_score <- mean(sil[, "sil_width"])
  
  end_time <- Sys.time()
  
  # save
  results <- data.frame(
    Seed = seed,
    Metric = c("ARI", "NMI", "Silhouette", 'Time'),
    Value = c(ari, nmi, sil_score, as.numeric(end_time - start_time, units = "mins")),
    Description = c(
      "Adjusted Rand Index (0-1, higher is better)",
      "Normalized Mutual Information (0-1, higher is better)",
      "Mean Silhouette Width (-1 to 1, higher is better)",
      "Running time (min)"
    )
  )
  result_file <- paste0("citefuse_results_", seed, ".csv")
  write.csv(results, result_file, row.names = FALSE)
  rdata_file <- paste0("citefuse_", seed, ".RData")
  save.image(rdata_file)
  message(paste("Completed CiteFuse analysis with seed", seed))
}

all_results <- lapply(random_seeds, function(seed) {
  df <- read.csv(paste0("citefuse_results_", seed, ".csv"))
  data.frame(
    Seed = seed,
    ARI = df$Value[df$Metric == "ARI"],
    NMI = df$Value[df$Metric == "NMI"],
    Silhouette = df$Value[df$Metric == "Silhouette"],
    Time_min = df$Value[df$Metric == "Time"]
  )
})
combined_results <- do.call(rbind, all_results)
write.csv(combined_results, "combined_citefuse_results_wide.csv", row.names = FALSE)


sessionInfo()
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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cluster_2.1.6                 aricode_1.0.3                
# [3] mclust_6.1.1                  pbmcMultiome.SeuratData_0.1.4
# [5] bmcite.SeuratData_0.3.0       SeuratData_0.2.2.9002        
# [7] Seurat_5.3.0                  SeuratObject_5.1.0           
# [9] sp_2.2-0                      patchwork_1.3.0              
# [11] DT_0.33                       scater_1.32.1                
# [13] ggplot2_3.5.2                 scuttle_1.14.0               
# [15] SingleCellExperiment_1.26.0   SummarizedExperiment_1.34.0  
# [17] Biobase_2.64.0                GenomicRanges_1.56.2         
# [19] GenomeInfoDb_1.40.1           IRanges_2.38.1               
# [21] S4Vectors_0.42.1              BiocGenerics_0.50.0          
# [23] MatrixGenerics_1.16.0         matrixStats_1.5.0            
# [25] CiteFuse_1.16.0              
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22          splines_4.4.1             later_1.3.2              
# [4] tibble_3.2.1              polyclip_1.10-7           fastDummies_1.7.5        
# [7] lifecycle_1.0.4           edgeR_4.2.2               globals_0.18.0           
# [10] lattice_0.22-6            MASS_7.3-61               magrittr_2.0.3           
# [13] limma_3.60.6              plotly_4.10.4             metapod_1.12.0           
# [16] httpuv_1.6.15             sctransform_0.4.2         spam_2.11-1              
# [19] spatstat.sparse_3.1-0     reticulate_1.42.0         cowplot_1.1.3            
# [22] pbapply_1.7-2             bayesm_3.1-6              RColorBrewer_1.1-3       
# [25] abind_1.4-8               zlibbioc_1.50.0           Rtsne_0.17               
# [28] purrr_1.0.2               mixtools_2.0.0.1          ggraph_2.2.1             
# [31] tensorA_0.36.2.1          rappdirs_0.3.3            tweenr_2.0.3             
# [34] GenomeInfoDbData_1.2.12   ggrepel_0.9.6             irlba_2.3.5.1            
# [37] listenv_0.9.1             spatstat.utils_3.1-3      pheatmap_1.0.12          
# [40] goftest_1.2-3             RSpectra_0.16-2           dqrng_0.4.1              
# [43] spatstat.random_3.3-3     fitdistrplus_1.2-2        parallelly_1.43.0        
# [46] DelayedMatrixStats_1.26.0 codetools_0.2-20          DelayedArray_0.30.1      
# [49] ggforce_0.4.2             tidyselect_1.2.1          UCSC.utils_1.0.0         
# [52] farver_2.1.2              ScaledMatrix_1.12.0       viridis_0.6.5            
# [55] spatstat.explore_3.4-2    jsonlite_2.0.0            BiocNeighbors_1.22.0     
# [58] tidygraph_1.3.1           progressr_0.15.1          ggridges_0.5.6           
# [61] survival_3.7-0            dbscan_1.2.2              segmented_2.1-4          
# [64] tools_4.4.1               ica_1.0-3                 Rcpp_1.0.12              
# [67] glue_1.7.0                gridExtra_2.3             SparseArray_1.4.8        
# [70] dplyr_1.1.4               withr_3.0.2               fastmap_1.2.0            
# [73] bluster_1.14.0            rhdf5filters_1.16.0       digest_0.6.36            
# [76] rsvd_1.0.5                R6_2.6.1                  mime_0.12                
# [79] scattermore_1.2           tensor_1.5                dichromat_2.0-0.1        
# [82] spatstat.data_3.1-6       tidyr_1.3.1               generics_0.1.4           
# [85] data.table_1.17.0         robustbase_0.99-4-1       graphlayouts_1.2.2       
# [88] httr_1.4.7                htmlwidgets_1.6.4         S4Arrays_1.4.1           
# [91] uwot_0.2.3                pkgconfig_2.0.3           gtable_0.3.6             
# [94] lmtest_0.9-40             XVector_0.44.0            htmltools_0.5.8.1        
# [97] dotCall64_1.2             scales_1.4.0              png_0.1-8                
# [100] spatstat.univar_3.1-2     scran_1.32.0              rstudioapi_0.17.1        
# [103] reshape2_1.4.4            nlme_3.1-165              rhdf5_2.48.0             
# [106] zoo_1.8-14                cachem_1.1.0              stringr_1.5.1            
# [109] KernSmooth_2.23-24        parallel_4.4.1            miniUI_0.1.2             
# [112] vipor_0.4.7               pillar_1.10.2             grid_4.4.1               
# [115] vctrs_0.6.5               RANN_2.6.2                randomForest_4.7-1.2     
# [118] promises_1.3.2            BiocSingular_1.20.0       beachmat_2.20.0          
# [121] xtable_1.8-4              beeswarm_0.4.0            locfit_1.5-9.12          
# [124] cli_3.6.3                 compiler_4.4.1            rlang_1.1.4              
# [127] crayon_1.5.3              future.apply_1.11.3       labeling_0.4.3           
# [130] plyr_1.8.9                ggbeeswarm_0.7.2          stringi_1.8.7            
# [133] viridisLite_0.4.2         deldir_2.0-4              BiocParallel_1.38.0      
# [136] lazyeval_0.2.2            spatstat.geom_3.3-6       compositions_2.0-8       
# [139] Matrix_1.7-0              RcppHNSW_0.6.0            sparseMatrixStats_1.16.0 
# [142] future_1.33.2             Rhdf5lib_1.26.0           statmod_1.5.0            
# [145] shiny_1.10.0              kernlab_0.9-33            ROCR_1.0-11              
# [148] igraph_2.1.4              memoise_2.0.1             DEoptimR_1.1-3-1    
