
                                                                            ![epbmclogo](analysis/data/derived_data/horse_DNA_image.png#style=centerme)

# Single cell resolution landscape of equine peripheral blood mononuclear cells reveals diverse cell types including T-bet<sup>+</sup> B cells.

This repository contains the data analysis scripts and code for our
paper:

> Roosheel S. Patel*, Joy E. Tomlinson*, Thomas J. Divers, Gerlinde R.
> Van de Walle, Brad R. Rosenberg, (2020). *Single cell resolution
> landscape of equine peripheral blood mononuclear cells reveals diverse
> cell types including T-bet+ B cells*. BMC Biology <DOI pending>

Our pre-print is online here:

> Roosheel S. Patel*, Joy E. Tomlinson*, Thomas J. Divers, Gerlinde R.
> Van de Walle, Brad R. Rosenberg, (2020). *Single cell resolution
> landscape of equine peripheral blood mononuclear cells reveals diverse
> cell types including T-bet+ B cells* biorxiv, Accessed 30 Dec 2020.
> Online at <https://doi.org/10.1101/2020.05.05.077362>

## Contents

The **analysis** directory contains:

  - [:file\_folder: paper](/analysis/paper): R Markdown document which
    includes code to reproduce the figures generated by our analysis.
  - [:file\_folder: data](/analysis/data): Supplementary data files and
    metadata used in this analysis

## Additional files required for analysis

For surface annotation of differentially expressed genes, surface
protein annotations were pulled from the [Cell Surface Protein Atlas
website](https://wlab.ethz.ch/cspa/#abstract).

> Bausch-Fluck D, Hofmann A, Bock T, Frei AP, Cerciello F, et al. (2015)
> *A Mass Spectrometric-Derived Cell Surface Protein Atlas*. PLoS One
> <https://doi.org/10.1371/journal.pone.0121314>

## Licenses

**Code :** MIT License, Copyright (c) 2020 Roosheel Patel

The organization of this repository was adapted from the rrtools
research compendium structure

> Marwick, B. (2017). Computational reproducibility in archaeological
> research: Basic principles and a case study of their
> implementation.*Journal of Archaeological Method and Theory*, 24(2),
> 424-450.<https://doi.org/10.1007/s10816-015-9272-9>

## Session Information

Session information (computer environment and R/package versions) used
to generate the analyses and figures in the manuscript:

<details>

> R version 3.6.2 (2019-12-12) Platform: x86\_64-pc-linux-gnu (64-bit)
> Running under: CentOS Linux 7 (Core) Matrix products: default
> BLAS/LAPACK:
> /hpc/packages/minerva-centos7/intel/parallel\_studio\_xe\_2019/compilers\_and\_libraries\_2019.0.117/linux/mkl/lib/intel64\_lin/libmkl\_gf\_lp64.so
> Random number generation: RNG: Mersenne-Twister Normal: Inversion
> Sample: Rounding locale: \[1\] LC\_CTYPE=en\_US.UTF-8 LC\_NUMERIC=C  
> \[3\] LC\_TIME=en\_US.UTF-8 LC\_COLLATE=en\_US.UTF-8  
> \[5\] LC\_MONETARY=en\_US.UTF-8 LC\_MESSAGES=en\_US.UTF-8  
> \[7\] LC\_PAPER=en\_US.UTF-8 LC\_NAME=C  
> \[9\] LC\_ADDRESS=C LC\_TELEPHONE=C  
> \[11\] LC\_MEASUREMENT=en\_US.UTF-8 LC\_IDENTIFICATION=C  
> other attached packages: \[1\] readr\_1.4.0 modelr\_0.1.8 tidyr\_1.1.2
> reshape2\_1.4.4 here\_1.0.0 circlize\_0.4.11 ComplexHeatmap\_2.2.0
> \[8\] GEOquery\_2.54.1 factoextra\_1.0.7 cluster\_2.1.0
> openxlsx\_4.2.3 legocolors\_0.2.0 ggforce\_0.3.2 ggtree\_2.0.4  
> \[15\] viridis\_0.5.1 viridisLite\_0.3.0 BiocManager\_1.30.10
> jntools\_0.1.0 mgsub\_1.7.2 ape\_5.4-1 patchwork\_1.1.0  
> \[22\] dplyr\_1.0.2 edgeR\_3.28.1 limma\_3.42.2 scales\_1.1.1
> lemon\_0.4.5 ggpubr\_0.4.0 readxl\_1.3.1  
> \[29\] ggrepel\_0.8.2 RColorBrewer\_1.1-2 ggnetwork\_0.5.8
> sctransform\_0.2.0 uwot\_0.1.9 rlang\_0.4.9 rrtools\_0.1.0  
> \[36\] future\_1.20.1 clustree\_0.4.3 ggraph\_2.0.4 ggplot2\_3.3.2
> gdata\_2.18.0 cowplot\_1.1.0 gtools\_3.8.2  
> \[43\] tibble\_3.0.4 AnnotationDbi\_1.48.0 IRanges\_2.20.2
> S4Vectors\_0.24.4 Biobase\_2.46.0 BiocGenerics\_0.32.0
> Matrix\_1.2-18  
> \[50\] renv\_0.12.3 Seurat\_3.1.0  
> loaded via a namespace (and not attached): \[1\] reticulate\_1.18
> tidyselect\_1.1.0 RSQLite\_2.2.1 htmlwidgets\_1.5.2 Rtsne\_0.15
> devtools\_2.3.2  
> \[7\] munsell\_0.5.0 codetools\_0.2-18 ica\_1.0-2 miniUI\_0.1.1.1
> withr\_2.3.0 colorspace\_2.0-0  
> \[13\] knitr\_1.30 rstudioapi\_0.13 ROCR\_1.0-11 ggsignif\_0.6.0
> tensor\_1.5 listenv\_0.8.0  
> \[19\] labeling\_0.4.2 git2r\_0.27.1 polyclip\_1.10-0 bit64\_4.0.5
> farver\_2.0.3 rprojroot\_2.0.2  
> \[25\] treeio\_1.10.0 parallelly\_1.21.0 vctrs\_0.3.5 generics\_0.1.0
> xfun\_0.19 R6\_2.5.0  
> \[31\] clue\_0.3-58 graphlayouts\_0.7.1 rsvd\_1.0.3 locfit\_1.5-9.4
> spatstat.utils\_1.17-0 assertthat\_0.2.1  
> \[37\] promises\_1.1.1 gtable\_0.3.0 globals\_0.14.0 processx\_3.4.5
> goftest\_1.2-2 tidygraph\_1.2.0  
> \[43\] clisymbols\_1.2.0 GlobalOptions\_0.1.2 splines\_3.6.1
> rstatix\_0.6.0 lazyeval\_0.2.2 broom\_0.7.2  
> \[49\] yaml\_2.2.1 abind\_1.4-5 backports\_1.2.1 httpuv\_1.5.4
> tools\_3.6.1 usethis\_1.6.3  
> \[55\] bookdown\_0.21 ellipsis\_0.3.1 sessioninfo\_1.1.1
> ggridges\_0.5.2 Rcpp\_1.0.5 plyr\_1.8.6  
> \[61\] purrr\_0.3.4 ps\_1.5.0 prettyunits\_1.1.1 rpart\_4.1-15
> deldir\_0.2-3 GetoptLong\_1.0.4  
> \[67\] pbapply\_1.4-3 zoo\_1.8-8 haven\_2.3.1 fs\_1.5.0
> magrittr\_2.0.1 data.table\_1.13.4  
> \[73\] lmtest\_0.9-38 RANN\_2.6.1 fitdistrplus\_1.1-3
> matrixStats\_0.57.0 pkgload\_1.1.0 hms\_0.5.3  
> \[79\] mime\_0.9 evaluate\_0.14 xtable\_1.8-4 rio\_0.5.16 shape\_1.4.5
> gridExtra\_2.3  
> \[85\] testthat\_3.0.0 compiler\_3.6.1 KernSmooth\_2.23-18
> crayon\_1.3.4 htmltools\_0.5.0 mgcv\_1.8-33  
> \[91\] later\_1.1.0.1 DBI\_1.1.0 tweenr\_1.0.1 MASS\_7.3-53
> car\_3.0-10 cli\_2.2.0  
> \[97\] igraph\_1.2.6 forcats\_0.5.0 pkgconfig\_2.0.3 rvcheck\_0.1.8
> foreign\_0.8-71 plotly\_4.9.2.1  
> \[103\] xml2\_1.3.2 stringr\_1.4.0 callr\_3.5.1 digest\_0.6.27
> RcppAnnoy\_0.0.17 spatstat.data\_1.5-2  
> \[109\] rmarkdown\_2.5 cellranger\_1.1.0 leiden\_0.3.6 tidytree\_0.3.3
> curl\_4.3 shiny\_1.5.0  
> \[115\] rjson\_0.2.20 lifecycle\_0.2.0 nlme\_3.1-150 jsonlite\_1.7.1
> carData\_3.0-4 desc\_1.2.0  
> \[121\] fansi\_0.4.1 pillar\_1.4.7 lattice\_0.20-41 fastmap\_1.0.1
> httr\_1.4.2 pkgbuild\_1.1.0  
> \[127\] survival\_3.2-7 glue\_1.4.2 remotes\_2.2.0 zip\_2.1.1
> spatstat\_1.64-1 png\_0.1-7  
> \[133\] bit\_4.0.4 stringi\_1.5.3 blob\_1.2.1 memoise\_1.1.0
> irlba\_2.3.3 future.apply\_1.6.0
