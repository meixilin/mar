
loadMARdependencies_ <- function() {
    # set.seed(1) # because some models are stochastic
    # tidy packages
    library(tidyverse)
    library(knitr)
    library(dplyr)
    # plotting
    library(ggplot2)
    library(cowplot)
    theme_set(theme_cowplot())
    library(RColorBrewer)
    library(latex2exp)
    # spatial
    library(raster)
    library(sp)
    library(sads)
    library(sars)
    # computational
    library(Rcpp)
    library(RcppArmadillo)
    library(adehabitatHR)
}
loadMARdependencies <- function() {
    suppressMessages(loadMARdependencies_())
}
