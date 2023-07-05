# elephantsealPhenology
Code and data for: Oosthuizen WC, Pistorius PA, Bester MN, Altwegg, R, de Bruyn PJN. 2023. Reproductive phenology is a repeatable, heritable trait linked to the timing of other life history events in a migratory marine predator. Proceedings of the Royal Society B. DOI: 10.1098/rspb.2023-1170.

This code contains updates (is a live version) of the Dryad version associated with the paper. 

We used mixed models to quantify how reproductive timing varied across 1,755 female southern elephant seals (Mirounga leonina) breeding at Marion Island in the Southern Ocean (1989 - 2019), and to identify the factors that correlate with phenological shifts within- and between individuals. We found strong support for covariation in the timing of breeding arrival dates and the timing of the preceding moult. Breeding arrival dates were more repeatable at the individual-level, as compared to the population-level, even after accounting for individual traits (wean date as a pup, age and breeding experience) associated with phenological variability. Mother-daughter similarities in breeding phenology were also evident, indicating that additive genetic effects may contribute to between-individual variation in breeding phenology. Over 30 years, elephant seal phenology did not change towards earlier or later dates, and we found no correlation between annual fluctuations in phenology and indices of environmental variation. 

## Sharing/Access information

This long-term data will be best used in collaboration with the authors, given their knowledge of how the data were collected. Please contact us to collaborate if you are interested in using this data for further work.  

## Code/Software

Analyses were conducted in R (v4.1.3 and later versions) (R Core Team. 2022). Analysis can be replicated using: OosthuizenWC_ProceedingsB2023-1170_all code.

utils::sessionInfo()
R version 4.2.1 (2022-06-23 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=English_South Africa.utf8  LC_CTYPE=English_South Africa.utf8   
[3] LC_MONETARY=English_South Africa.utf8 LC_NUMERIC=C                         
[5] LC_TIME=English_South Africa.utf8    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] suncalc_0.5.1       performance_0.10.0  boot_1.3-28         rptR_0.9.22        
 [5] scales_1.2.1        optimx_2022-4.30    cowplot_1.1.1       rsoi_0.5.5         
 [9] ggeffects_1.1.3     sjPlot_2.8.11       chron_2.3-58        broom.mixed_0.2.9.4
[13] broom_1.0.4         lme4_1.1-30         Matrix_1.5-1        lubridate_1.8.0    
[17] forcats_1.0.0       stringr_1.5.0       dplyr_1.1.1         purrr_1.0.1        
[21] readr_2.1.3         tidyr_1.3.0         tibble_3.2.1        ggplot2_3.4.2      
[25] tidyverse_2.0.0    

loaded via a namespace (and not attached):
 [1] viridis_0.6.2        viridisLite_0.4.1    splines_4.2.1       
 [4] modelr_0.1.11        datawizard_0.7.1     colorBlindness_0.1.9
 [7] bayestestR_0.13.0    globals_0.16.2       numDeriv_2016.8-1.1 
[10] pillar_1.9.0         backports_1.4.1      lattice_0.20-45     
[13] glue_1.6.2           digest_0.6.31        RColorBrewer_1.1-3  
[16] minqa_1.2.4          colorspace_2.1-0     plyr_1.8.8          
[19] pkgconfig_2.0.3      listenv_0.9.0        xtable_1.8-4        
[22] mvtnorm_1.1-3        tzdb_0.3.0           emmeans_1.8.1-1     
[25] generics_0.1.3       sjlabelled_1.2.0     withr_2.5.0         
[28] furrr_0.3.1          cli_3.4.1            magrittr_2.0.3      
[31] effectsize_0.8.1     estimability_1.4.1   GGally_2.1.2        
[34] future_1.32.0        fansi_1.0.4          parallelly_1.35.0   
[37] nlme_3.1-160         MASS_7.3-57          data.table_1.14.8   
[40] tools_4.2.1          hms_1.1.3            lifecycle_1.0.3     
[43] munsell_0.5.0        compiler_4.2.1       gridGraphics_0.5-1  
[46] rlang_1.1.0          grid_4.2.1           nloptr_2.0.3        
[49] parameters_0.19.0    rstudioapi_0.14      gtable_0.3.3        
[52] codetools_0.2-18     sjstats_0.18.1       reshape_0.8.9       
[55] sjmisc_2.8.9         R6_2.5.1             gridExtra_2.3       
[58] knitr_1.42           utf8_1.2.3           insight_0.19.1      
[61] stringi_1.7.12       parallel_4.2.1       Rcpp_1.0.10         
[64] vctrs_0.6.1          tidyselect_1.2.0     xfun_0.38           
[67] coda_0.19-4         
