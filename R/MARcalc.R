# calculate MAR/EMAR relationship
.Mtype = c("M", "E", "thetaw", "thetapi")
.Atype = c("A", "Asq")

#' Title
#'
#' @param mardf
#' @param Mtype
#' @param Atype
#'
#' @return
#' @export
#'
#' @examples
MARcalc <- function(mardf, Mtype = .Mtype, Atype = .Atype) {
    Mtype = match.arg(Mtype)
    Atype = match.arg(Atype)
    # remove NA or zero data
    tmpdf = mardf[, c(Atype, Mtype)]
    tmpdf = tmpdf[(tmpdf[,Mtype] > 0 & !is.na(tmpdf[,Mtype])), ]
    # previously had a check for length(unique(tmpdf[,Mtype])) < 4
    if (nrow(tmpdf) == 0) {
        warning(paste0("MAR for Mtype = ", Mtype, ", Atype = ", Atype, " cannot be calculated"))
        mar <- NULL
    } else {
        # run MAR analyses
        mar <- sars::sar_power(tmpdf)
    }
    return(mar)
}

.marsummary <- function(mar) {
    # Default output structure
    default_outdf <- list(
        model = NA_character_,
        c = NA_real_,
        z = NA_real_,
        c_p = NA_real_,
        z_p = NA_real_,
        R2_adj = NA_real_
    )

    outdf <- tryCatch({
        marsum <- summary(mar)
        outdf <- list(
            model = marsum$Model,
            c = marsum$Parameters[[1, "Estimate"]],
            z = marsum$Parameters[[2, "Estimate"]],
            c_p = marsum$Parameters[[1, "Pr(>|t|)"]],
            z_p = marsum$Parameters[[2, "Pr(>|t|)"]],
            R2_adj = marsum$R2a
        )
        return(outdf)
    },
    error = function(e) {
        return(default_outdf)
    })

    return(outdf)
}

.pipe_MARcalc <- function(mardf, verbose = TRUE) {
    message(paste0("MAR built from scheme: ", attr(mardf, "scheme")))
    MA <- expand.grid(M = .Mtype, A = .Atype, stringsAsFactors = FALSE)
    mars <- apply(MA, 1, function(x) MARcalc(mardf, x[1], x[2]))
    names(mars) <- apply(MA, 1, paste, collapse = "_")
    marsuml <- lapply(mars, .marsummary)
    output <- cbind(MA, do.call(rbind, lapply(marsuml, as.data.frame, stringsAsFactors = FALSE)))
    if (verbose) {
        print(output)
    }
    class(output) <- c(class(output), "marcalc")
    return(output)
}
