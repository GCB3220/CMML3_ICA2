# CMML3 ICA2 

This is a repository storing the CCML3 ICA2 code.

This repository stores three R scripts. Among them, `seurat.r` and `citefuse.r` are used for processing and integrating CITE-seq data using Seurat and CiteFuse respectively, as well as calculating the corresponding ARI, NMI, and Silhouette Score. The `plot.r` script is used to draw bar charts that display these metrics and the running time. All three scripts have been tested and run in the following environment:
```
> sessionInfo()
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8 
[2] LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] patchwork_1.3.0               DT_0.33                      
 [3] scater_1.32.1                 scuttle_1.14.0               
 [5] SingleCellExperiment_1.26.0   SummarizedExperiment_1.34.0  
 [7] Biobase_2.64.0                GenomicRanges_1.56.2         
 [9] GenomeInfoDb_1.40.1           IRanges_2.38.1               
[11] S4Vectors_0.42.1              BiocGenerics_0.50.0          
[13] MatrixGenerics_1.16.0         matrixStats_1.5.0            
[15] CiteFuse_1.16.0               cluster_2.1.6                
[17] aricode_1.0.3                 mclust_6.1.1                 
[19] cowplot_1.1.3                 pbmcMultiome.SeuratData_0.1.4
[21] bmcite.SeuratData_0.3.0       SeuratData_0.2.2.9002        
[23] Seurat_5.3.0                  SeuratObject_5.1.0           
[25] sp_2.2-0                      dplyr_1.1.4                  
[27] ggplot2_3.5.2                

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22          splines_4.4.1             later_1.3.2              
  [4] tibble_3.2.1              polyclip_1.10-7           fastDummies_1.7.5        
  [7] lifecycle_1.0.4           edgeR_4.2.2               globals_0.18.0           
 [10] lattice_0.22-6            MASS_7.3-61               magrittr_2.0.3           
 [13] limma_3.60.6              plotly_4.10.4             metapod_1.12.0           
 [16] httpuv_1.6.15             sctransform_0.4.2         spam_2.11-1              
 [19] spatstat.sparse_3.1-0     reticulate_1.42.0         pbapply_1.7-2            
 [22] bayesm_3.1-6              RColorBrewer_1.1-3        abind_1.4-8              
 [25] zlibbioc_1.50.0           Rtsne_0.17                purrr_1.0.2              
 [28] mixtools_2.0.0.1          ggraph_2.2.1              tensorA_0.36.2.1         
 [31] rappdirs_0.3.3            tweenr_2.0.3              GenomeInfoDbData_1.2.12  
 [34] ggrepel_0.9.6             irlba_2.3.5.1             listenv_0.9.1            
 [37] spatstat.utils_3.1-3      pheatmap_1.0.12           goftest_1.2-3            
 [40] RSpectra_0.16-2           dqrng_0.4.1               spatstat.random_3.3-3    
 [43] fitdistrplus_1.2-2        parallelly_1.43.0         DelayedMatrixStats_1.26.0
 [46] codetools_0.2-20          DelayedArray_0.30.1       ggforce_0.4.2            
 [49] tidyselect_1.2.1          UCSC.utils_1.0.0          farver_2.1.2             
 [52] ScaledMatrix_1.12.0       viridis_0.6.5             spatstat.explore_3.4-2   
 [55] jsonlite_2.0.0            BiocNeighbors_1.22.0      tidygraph_1.3.1          
 [58] progressr_0.15.1          ggridges_0.5.6            survival_3.7-0           
 [61] dbscan_1.2.2              segmented_2.1-4           tools_4.4.1              
 [64] ica_1.0-3                 Rcpp_1.0.12               glue_1.7.0               
 [67] gridExtra_2.3             SparseArray_1.4.8         withr_3.0.2              
 [70] fastmap_1.2.0             bluster_1.14.0            rhdf5filters_1.16.0      
 [73] digest_0.6.36             rsvd_1.0.5                R6_2.6.1                 
 [76] mime_0.12                 scattermore_1.2           tensor_1.5               
 [79] dichromat_2.0-0.1         spatstat.data_3.1-6       tidyr_1.3.1              
 [82] generics_0.1.4            data.table_1.17.0         robustbase_0.99-4-1      
 [85] graphlayouts_1.2.2        httr_1.4.7                htmlwidgets_1.6.4        
 [88] S4Arrays_1.4.1            uwot_0.2.3                pkgconfig_2.0.3          
 [91] gtable_0.3.6              lmtest_0.9-40             XVector_0.44.0           
 [94] htmltools_0.5.8.1         dotCall64_1.2             scales_1.4.0             
 [97] png_0.1-8                 spatstat.univar_3.1-2     scran_1.32.0             
[100] rstudioapi_0.17.1         reshape2_1.4.4            nlme_3.1-165             
[103] rhdf5_2.48.0              zoo_1.8-14                cachem_1.1.0             
[106] stringr_1.5.1             KernSmooth_2.23-24        parallel_4.4.1           
[109] miniUI_0.1.2              vipor_0.4.7               pillar_1.10.2            
[112] grid_4.4.1                vctrs_0.6.5               RANN_2.6.2               
[115] randomForest_4.7-1.2      promises_1.3.2            BiocSingular_1.20.0      
[118] beachmat_2.20.0           xtable_1.8-4              beeswarm_0.4.0           
[121] locfit_1.5-9.12           cli_3.6.3                 compiler_4.4.1           
[124] rlang_1.1.4               crayon_1.5.3              future.apply_1.11.3      
[127] labeling_0.4.3            plyr_1.8.9                ggbeeswarm_0.7.2         
[130] stringi_1.8.7             viridisLite_0.4.2         deldir_2.0-4             
[133] BiocParallel_1.38.0       lazyeval_0.2.2            spatstat.geom_3.3-6      
[136] compositions_2.0-8        Matrix_1.7-0              RcppHNSW_0.6.0           
[139] sparseMatrixStats_1.16.0  future_1.33.2             Rhdf5lib_1.26.0          
[142] statmod_1.5.0             shiny_1.10.0              kernlab_0.9-33           
[145] ROCR_1.0-11               igraph_2.1.4              memoise_2.0.1            
[148] DEoptimR_1.1-3-1
```