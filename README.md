# GSFM
These files are the key script for the GSFM
Code used for the analysis of the paper:

GSFM ：A Genome-Scale Functional Module transformation to represent drug efficacy for in silico drug discovery.

Saisai Tian, Xuyang Liao, Wen Cao, Xinyi Wu, Zexi Chen, Jinyuan Lu, Qun Wang, Jinbo Zhang, Luonan Chen, Weidong Zhang

# STEP 1 Data processing
Description: Process data from TCGA and LINCS databases,and further measured the correlation between each sample and cell line, and removed samples that are not correlated to the cell lines.
## input 
     Cancer Gene Expression Matrix
     Drug-induced Gene Expression Matrix
     CCLE data
     
# STEP 2 Data transformation
Description: Calculation of the four biologically interpretable quantifiers (BIQs): GSFM_Up, GSFM_Down, GSFM_ssGSEA, and GSFM_TF.
## input
     Processing Cancer Gene Expression Matrix (step1)
     Processing Drug-induced Gene Expression Matrix (step1)
     Hallmark genesets
     Database:TF-target pairs by Garcia-Alonso et al.
## Build docker
     cd 2.data transformation/docker
     docker build --tag gsfm-script:latest --file GSFM.dockerfile .
     cd ..
     docker run -it -d  --restart=always --name GSFM-notebook   -p 12101:8888 --log-opt max-size=10m --log-opt max-file=5 -v `pwd`/project:/project   gsfm-script
     docker container exec -it GSFM-notebook bash
     jupyter-notebook --ip=0.0.0.0 --port=8888 --allow-root --no-browser
     run ./project/Notebooks BRCA_GSFM.ipynb

# STEP 3 Drug discovery
Description: Calculation of RS-GSFM, RS-Gene, RS-Viper                                                                                 
 step3.1:Calculation of RS-GSFM                                
 step3.2:Calculation of RS-Gene                                      
 step3.3:Calculation of RS-Viper                              
## input
     Processing Cancer Gene Expression Matrix (step1)
     Processing Drug-induced Gene Expression Matrix (step1) 
     Processing Cancer GSFM Active Matrix (step2)
     Processing Drug-induced GSFM Active Matrix  (step2)
     

# System requirement
R studio 4.1.1 R version 4.1.1 (2021-08-10)

Docker version:24.0.6

Platform: x86_64-pc-linux-gnu (64-bit)

Running under: Ubuntu 20.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3;  LAPACK version 3.9.0

locale:

 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:

[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:

 [1] ggrepel_0.9.1         RColorBrewer_1.1.3    org.Hs.eg.db_3.12.0     
 [4] shinyWidgets_0.5.7    plotly_4.9.3          signatureSearch_1.4.6      
 [7] leaflet_2.1.1         DT_0.17               limma_3.46.0   
[10] shinyBS_0.61          ROCR_1.0.11           stringr_1.5.0        
[13] purrr_1.0.1           readr_1.4.0           tibble_3.1.8         
[16] tidyverse_1.3.0       tidyr_1.2.0           dplyr_1.0.8          
[19] ggplot2_3.3.5         pROC_1.17.0.1         memoise_2.0.1        
[22] cowplot_1.1.1          

loaded via a namespace (and not attached):
 [1] methods_4.3.1     rio_0.5.29        utf8_1.2.3        cellranger_1.1.0  readxl_1.4.2      magrittr_2.0.3    glue_1.6.2       
 [8] tibble_3.2.1      foreign_0.8-85    pkgconfig_2.0.3   lifecycle_1.0.3   utils_4.3.1       cli_3.6.1         zip_2.3.0        
[15] fansi_1.0.4       openxlsx_4.2.5.2  vctrs_0.6.3       graphics_4.3.1    data.table_1.14.8 grDevices_4.3.1   stats_4.3.1      
[22] compiler_4.3.1    forcats_1.0.0     haven_2.5.2       base_4.3.1        rstudioapi_0.14   tools_4.3.1       curl_5.0.2       
[29] hms_1.1.3         pillar_1.9.0      Rcpp_1.0.11       rlang_1.1.1       stringi_1.7.12    datasets_4.3.1 



     
    
