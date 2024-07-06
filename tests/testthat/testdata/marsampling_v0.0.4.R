# Title: Generate reproducible MARsampling results using Moi's original code
# Author: Meixi Lin
# used to create test files to test backward compatibility before correcting some newly discovered bugs

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

# setwd("tests/testthat/testdata")
library(mar)
sessionInfo()
# v0.0.3[4] c581a7e986e9bf330fd2763901f864e4df318d93

# def functions --------
myMARsampling <- function(genemaps, myscheme) {
    set.seed(7)
    # TOFIX: only when run interactively in Rmarkdown console, MARextinction_sim outputs waring message
    # no non-missing arguments to min; returning Infno non-missing arguments to max; returning -Inf
    outdt = mar:::MARsampling(genemaps, scheme = myscheme)
    return(outdt)
}

# def variables --------
# only check the "random" scheme as the MARPIPELINE did not include inwards or outwards options
# myschemes = c("random", "inwards", "outwards")

# fix the errors on create_gene_maps --------
load(file = 'coords-joshua.rda')
load(file = 'genome-joshua.rda')
load(file = 'genemaps-joshua.rda')
newmaps <- create_gene_maps(coords, genomes, bystep = 0.005)

# confirm that the new genemaps have 290 individuals
cellStats(newmaps[[1]], 'sum')
plot(newmaps[[1]])
plot(genemaps[[1]], add = T)

# save the new genemaps
save(newmaps, file = 'genemaps_new-joshua.rda')

# run marsampling --------
# # before adding the debug option
# yy = myMARsampling(newmaps, "random")
# oyy = yy
# # add debug option
# yy = myMARsampling(newmaps, "random")
# # check previous run without extent gives the same first 7 columns
# all.equal(oyy, yy[,-8]) # T

# formally run and save the output from MARsampling
mares <- myMARsampling(newmaps, "random")
save(mares, file = "mares_new-joshua.rda")

# cleanup --------
date()
closeAllConnections()
