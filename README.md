
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Causality Explained

Further details will be provided after the preprint is submitted.

## Installation

To install the current version of Causality Explained:

``` r
  if(!"devtools" %in% installed.packages()[, "Package"]) {
  install.packages("devtools") }

  if(!"CausalityExplained" %in% installed.packages()[, "Package"]) {
  devtools::install_github("martintnr/CausalityExplained") }
```

## Example

To process the coronary heart disease trait from CARDIoGRAMplusC4D using
the Causality Explained framework, with phenome-wide UK Biobank data
serving as predictors:

``` r
library(data.table)
library(parallel)
library(CausalityExplained)

setwd("where/you/want/CausalityExplained")

system("wget -O Data.zip  https://zenodo.org/records/14860110/files/Data.zip?download=1")
system("unzip Data.zip && rm Data.zip")

Results <- CausalityExplained::Causality_Explained(CutOff = 5e-08, DataPath = "Data/",
       Outcomes <- c("CAD_CARDIoGRAMplusC4D.txt"),
       NBCores = 1, BonferroniCorrection = F)
```
