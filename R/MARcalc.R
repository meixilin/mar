
# calculate MAR/EMAR relationship
.Mtype = c("M", "E", "thetaw", "thetapi")
.Atype = c("A", "Asq")

MARcalc <- function(mardf, Mtype = .Mtype, Atype = .Atype) {
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

.marsummary <- function(mar) {
    marsum <- summary(mar)
    outdf <- list(
        model = marsum$Model,
        c = marsum$Parameters[1, "Estimate"],
        z = marsum$Parameters[2, "Estimate"],
        c_p = marsum$Parameters[1, "Pr(>|t|)"],
        z_p = marsum$Parameters[2, "Pr(>|t|)"],
        R2_adj = marsum$R2a
    )
    return(outdf)
}

MARcalc_all <- function(mardf, verbose = TRUE, return_mar = FALSE) {
    MA <- expand.grid(M = .Mtype, A = .Atype, stringsAsFactors = FALSE)
    mars <- apply(MA, 1, function(x) MARcalc(mardf, x[1], x[2]))
    names(mars) <- apply(MA, 1, paste, collapse = "_")
    marsummary <- cbind(MA, do.call(rbind, lapply(mars, .marsummary)))
    if (verbose) {
        print(marsummary)
    }
    if (return_mar) {
        output <- list(mars = mars, marsummary = marsummary)
    } else {
        output <- marsummary
    }
    class(output) <- "marcalc"
    return(output)
}


