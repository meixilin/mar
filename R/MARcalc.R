
# calculate MAR/EMAR relationship
MARcalc <- function(mardf, Mtype = c("M", "E", "thetaw", "thetapi"), Atype = c("A", "Asq")) {
    Mtype = match.arg(Mtype)
    Atype = match.arg(Atype)
    # remove NA or zero data
    tmpdf = mardf[, c(Atype, Mtype)]
    tmpdf = tmpdf[(tmpdf[,Mtype] > 0 & !is.na(tmpdf[,Mtype])), ]
    stopifnot(nrow(tmpdf) > 0)
    stopifnot(length(unique(tmpdf[,Mtype])) > 3)
    # run MAR analyses
    mar <- sars::sar_power(tmpdf)
    return(mar)
}



