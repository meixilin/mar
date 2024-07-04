# Title: Generate reproducible extinction-sim results using Moi's original code
# Author: Meixi Lin
# Date: Wed Apr 17 20:06:41 2024

# this file intends to become obselete at v0.1.0 when mutdiv calculation corrects for `sum(cells > 0)` bug
# used to create test files to test backward compatibility before correcting some newly discovered bugs

# preparation --------
rm(list = ls())
cat("\014")
options(echo = TRUE, stringsAsFactors = FALSE)

library(mar)
# version: v0.0.1 commit id: 3c576665cedd871ecf773addca8656110e90d9c2

# def functions --------
loadSomeRData <- function(vars, file) {
    E <- new.env()
    load(file=file, envir=E)
    return(mget(vars, envir=E, inherits=F))
}

load_prevdata <- function() {
    testdt = loadSomeRData("outdfl", file = '~/AMoiLab/mar_related/pi_extinct/data/reproduce_extinctionsim/v3c576/reproduce_extinctionsim.RData')[[1]][[1]]
    testdt1 = testdt[testdt$version == "v3c576" & testdt$type == "random", 1:9]
    testdt2 = testdt[testdt$version == "v3c576" & testdt$type == "southnorth", 1:9]
    rownames(testdt1) = NULL
    rownames(testdt2) = NULL
    return(list(testdt1, testdt2))
}


myMARextinction_sim <- function(genemaps, myscheme) {
    set.seed(7)
    # TOFIX: only when run interactively in Rmarkdown console, MARextinction_sim outputs waring message
    # no non-missing arguments to min; returning Infno non-missing arguments to max; returning -Inf
    outdt = mar:::MARextinction_sim(genemaps, scheme = myscheme)
    return(outdt)
}

# def variables --------
# outwards extinction not implemented yet
myschemes = c("random", "inwards", "southnorth", "radial")

# load data --------
species = 'joshua'
load('~/AMoiLab/mar_related/MAR2.0/pi_extinct/ARCHIVE/data-raw/tmpobjects/genemaps-joshua.rda')
# check that they are the same with pi_extinct tests
testdt = load_prevdata()

# main --------
outdts = lapply(myschemes, function(ss) {
    outdt = myMARextinction_sim(genemaps, ss)
    saveRDS(outdt, file = paste0('xsim-', species, '-', ss, '.rds'))
    return(outdt)
})

# check random and southnorth with previous
stopifnot(all.equal(outdts[[1]], testdt[[1]]))
stopifnot(all.equal(outdts[[3]], testdt[[2]]))

# cleanup --------
date()
closeAllConnections()

# the genemaps-joshua is copied to testdata too
# Rscript --vanilla extinction-sim_v0.0.1.R

# update after commit 6289 --------
ss = "radial"
outdt = myMARextinction_sim(genemaps, ss)
saveRDS(outdt, file = paste0('xsim-', species, '-', ss, '.rds'))




