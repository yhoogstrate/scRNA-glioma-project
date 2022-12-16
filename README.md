scRNA glioma project
================

-   [Introduction](#introduction)
-   [Dataset overview](#dataset-overview)

<!-- README.md is generated from README.Rmd through: devtools::build_readme() . Please edit that file -->
<!-- Nice example package: https://github.com/hadley/babynames -->

------------------------------------------------------------------------

![license](https://img.shields.io/badge/license-GPL--3-blue.svg)
[![GitHub
issues](https://img.shields.io/github/issues/yhoogstrate/scRNA-glioma-project.svg)]()
![rversion](https://img.shields.io/badge/R%20version-%3E4.2.0-lightgrey.svg)

# Introduction

This projects aims to clean-up and integrate publicly available glioma
s\[c/n\]RNA.

# Dataset overview

``` r
devtools::load_all()
```

    ## ℹ Loading scRNAgliomaProject

``` r
library(ggplot2)
```

``` r
ggplot(data_per_sample, ggplot2::aes(x=date_added, y=reorder(sid, date_added), col=IDH.status)) +
  geom_point() +
  labs(y=NULL)
```

``` r
data_per_sample |> 
  dplyr::arrange(date_added) |> 
  dplyr::select(sample_name, IDH.status, date_added, doi) |> 
  dplyr::rename(IDH = IDH.status) |> 
  dplyr::mutate(doi = gsub("^.+\\/","",doi))
```

    ##                 sample_name IDH date_added            doi
    ## 1 van Hijfte - GBM Sample Y  wt 2022-12-12 zenodo.6546712
