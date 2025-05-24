# This R script integrates CITE-seq data using Seurat.
# The specific process and parameter settings here refer to: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html  
# Since each run only samples 4,000 cells, the npcs parameter in the RunPCA function was modified to 10, and the dim.list parameter in the subsequent FindMultiModalNeighbors function was also adjusted accordingly (1:10). All other parameters remained unchanged.
# The dataset used is the bmcite dataset from SeuratData. 
# First, the data is randomly sampled to 4000 cells. Then, data preprocessing is performed, followed by integration of RNA and ADT data. Finally, cell clustering is carried out.
# The clustering performance is evaluated using ARI, NMI and Silhouette Score.
# The above process will loop 5 times.
# The performance of data integration is measured by taking the average of three indicators.
# All data were retained in the end
# Execution environment can be found in the sessioninfo at the end of the script.


# packages
library(Seurat)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(dplyr)
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
  # RNA
  DefaultAssay(bm) <- 'RNA'
  bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 10)
  # ADT
  DefaultAssay(bm) <- 'ADT'
  VariableFeatures(bm) <- rownames(bm[["ADT"]])
  bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
    ScaleData() %>% RunPCA(npcs = 10, reduction.name = 'apca')
  # interagtion
  bm <- FindMultiModalNeighbors(
    bm, reduction.list = list("pca", "apca"), 
    dims.list = list(1:10, 1:10), modality.weight.name = "RNA.weight"
  )
  # cluster
  bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  bm <- FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
  # plot
  p1 <- DimPlot(bm, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
  p2 <- DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = FALSE, repel = TRUE, label.size = 2.5) + 
    theme(legend.position = "right") 
  
  # benchmarking
  cluster_labels <- bm@meta.data$seurat_clusters 
  true_labels <- bm$celltype.l2  
  ari <- ARI(cluster_labels, true_labels)
  nmi <- NMI(cluster_labels, true_labels)
  wnn_embeddings <- Embeddings(bm, reduction = "wnn.umap")
  sil <- silhouette(as.numeric(cluster_labels), dist(wnn_embeddings))
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
  result_file <- paste0("seurat_results_", seed, ".csv")
  write.csv(results, result_file, row.names = FALSE)
  rdata_file <- paste0("seurat_", seed, ".RData")
  save.image(rdata_file)
  message(paste("Completed analysis with seed", seed))
}


all_results <- lapply(random_seeds, function(seed) {
  df <- read.csv(paste0("seurat_results_", seed, ".csv"))
  data.frame(
    Seed = seed,
    ARI = df$Value[df$Metric == "ARI"],
    NMI = df$Value[df$Metric == "NMI"],
    Silhouette = df$Value[df$Metric == "Silhouette"],
    Time_min = df$Value[df$Metric == "Time"]
  )
})
combined_results <- do.call(rbind, all_results)
write.csv(combined_results, "combined_seurat_results_wide.csv", row.names = FALSE)

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cluster_2.1.6                 aricode_1.0.3                
# [3] mclust_6.1.1                  dplyr_1.1.4                  
# [5] cowplot_1.1.3                 ggplot2_3.5.2                
# [7] pbmcMultiome.SeuratData_0.1.4 bmcite.SeuratData_0.3.0      
# [9] SeuratData_0.2.2.9002         Seurat_5.3.0                 
# [11] SeuratObject_5.1.0            sp_2.2-0                     
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3          rstudioapi_0.17.1          
# [3] jsonlite_2.0.0              magrittr_2.0.3             
# [5] ggbeeswarm_0.7.2            spatstat.utils_3.1-3       
# [7] farver_2.1.2                zlibbioc_1.50.0            
# [9] vctrs_0.6.5                 ROCR_1.0-11                
# [11] DelayedMatrixStats_1.26.0   spatstat.explore_3.4-2     
# [13] htmltools_0.5.8.1           S4Arrays_1.4.1             
# [15] BiocNeighbors_1.22.0        SparseArray_1.4.8          
# [17] sctransform_0.4.2           parallelly_1.43.0          
# [19] KernSmooth_2.23-24          htmlwidgets_1.6.4          
# [21] ica_1.0-3                   plyr_1.8.9                 
# [23] plotly_4.10.4               zoo_1.8-14                 
# [25] igraph_2.1.4                mime_0.12                  
# [27] lifecycle_1.0.4             pkgconfig_2.0.3            
# [29] rsvd_1.0.5                  Matrix_1.7-0               
# [31] R6_2.6.1                    fastmap_1.2.0              
# [33] GenomeInfoDbData_1.2.12     MatrixGenerics_1.16.0      
# [35] fitdistrplus_1.2-2          future_1.33.2              
# [37] shiny_1.10.0                digest_0.6.36              
# [39] patchwork_1.3.0             S4Vectors_0.42.1           
# [41] tensor_1.5                  scater_1.32.1              
# [43] RSpectra_0.16-2             irlba_2.3.5.1              
# [45] GenomicRanges_1.56.2        beachmat_2.20.0            
# [47] labeling_0.4.3              progressr_0.15.1           
# [49] spatstat.sparse_3.1-0       httr_1.4.7                 
# [51] polyclip_1.10-7             abind_1.4-8                
# [53] compiler_4.4.1              withr_3.0.2                
# [55] BiocParallel_1.38.0         viridis_0.6.5              
# [57] fastDummies_1.7.5           MASS_7.3-61                
# [59] rappdirs_0.3.3              DelayedArray_0.30.1        
# [61] tools_4.4.1                 vipor_0.4.7                
# [63] lmtest_0.9-40               beeswarm_0.4.0             
# [65] httpuv_1.6.15               future.apply_1.11.3        
# [67] goftest_1.2-3               glue_1.7.0                 
# [69] nlme_3.1-165                promises_1.3.2             
# [71] grid_4.4.1                  Rtsne_0.17                 
# [73] reshape2_1.4.4              generics_0.1.4             
# [75] gtable_0.3.6                spatstat.data_3.1-6        
# [77] tidyr_1.3.1                 data.table_1.17.0          
# [79] ScaledMatrix_1.12.0         BiocSingular_1.20.0        
# [81] XVector_0.44.0              BiocGenerics_0.50.0        
# [83] spatstat.geom_3.3-6         RcppAnnoy_0.0.22           
# [85] ggrepel_0.9.6               RANN_2.6.2                 
# [87] pillar_1.10.2               stringr_1.5.1              
# [89] spam_2.11-1                 RcppHNSW_0.6.0             
# [91] later_1.3.2                 splines_4.4.1              
# [93] lattice_0.22-6              survival_3.7-0             
# [95] deldir_2.0-4                tidyselect_1.2.1           
# [97] SingleCellExperiment_1.26.0 miniUI_0.1.2               
# [99] scuttle_1.14.0              pbapply_1.7-2              
# [101] gridExtra_2.3               IRanges_2.38.1             
# [103] SummarizedExperiment_1.34.0 scattermore_1.2            
# [105] stats4_4.4.1                Biobase_2.64.0             
# [107] matrixStats_1.5.0           stringi_1.8.7              
# [109] UCSC.utils_1.0.0            lazyeval_0.2.2             
# [111] codetools_0.2-20            tibble_3.2.1               
# [113] cli_3.6.3                   uwot_0.2.3                 
# [115] xtable_1.8-4                reticulate_1.42.0          
# [117] dichromat_2.0-0.1           Rcpp_1.0.12                
# [119] GenomeInfoDb_1.40.1         globals_0.18.0             
# [121] spatstat.random_3.3-3       png_0.1-8                  
# [123] spatstat.univar_3.1-2       parallel_4.4.1             
# [125] dotCall64_1.2               sparseMatrixStats_1.16.0   
# [127] listenv_0.9.1               viridisLite_0.4.2          
# [129] scales_1.4.0                ggridges_0.5.6             
# [131] purrr_1.0.2                 crayon_1.5.3               
# [133] rlang_1.1.4    

