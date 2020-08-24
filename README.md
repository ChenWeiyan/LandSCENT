
<!-- README.md is generated from README.Rmd. Please edit that file -->
LandSCENT package
=========

[![DOI](https://zenodo.org/badge/167301603.svg)](https://zenodo.org/badge/latestdoi/167301603)

<!-- badges: start -->
<!-- badges: end -->
`LandSCENT` (Landscape Single Cell Entropy) is a R-package for the analysis of single-cell RNA-Seq data. One important feature of this package is the computation of signaling entropy, which allows single cells to be ordered according to differentiation potency. LandSCENT also integrates cell density with potency distribution to dissect cell types across all potency states and generates high-quality figures to show this. In the latest version, we employ diffusion maps to infer trajectory within the data and also generate figures to show this.

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ChenWeiyan/LandSCENT")
```

Softerware Version
------------------

At the time of writing, we tested on version 0.99.3.

Operating system and dependencies
---------------------------------

1. This package is developed on Linux under version 3.10.0-693.5.2.el7.x86\_64 and Red Hat 4.8.5-16.

2. Using the most recent version of R is strongly recommended (R 3.6 at the time of writing).

3. The following are several packages from CRAN and Bioconductor that `LandSCENT` uses:

``` r
cluster (version >= 2.0.9), corpcor (version >= 1.6.9), igraph (version >= 1.2.4.1), isva (version >= 1.9), mclust (version >= 5.4.3), marray (version >= 1.62.0), scater (version >= 1.12.0), Biobase (version >= 2.44.0), BiocGenerics (version >= 0.30.0), SummarizedExperiment (version >= 1.14.0), SingleCellExperiment (version >= 1.6.0), Rtsne (version >= 0.15), irlba (version >= 2.3.3), plot3D (version >= 1.1.1), MASS(version >= 7.3-51.4), dbscan (version >= 1.1-3), monocle (version >= 2.12.0), DelayedArray (version >= 0.10.0), Matrix (version >= 1.2-17), destiny (version >= 2.14.0), ggplot2 (version >= 3.1.1), ggthemes (version >= 4.2.0), wordspace (version >= 0.2-6)
```
