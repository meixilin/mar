# calculate MAR relationship
.Mtype = c("M", "E", "thetaw", "thetapi")
.Atype = c("A", "Asq")

#' Calculate the mutations-area relationship (MAR) using the power-law model
#'
#' This function calculates the MAR using the power-law model from the provided data frame.
#' The data frame should contain columns for the diversity metric (e.g. richness, evenness)
#' and the area metric (e.g. total area, area squared).
#'
#' @param mardf A data frame with columns of diversity metric and area.
#'        Output from \link{MARsampling} or \link{MARextinction}.
#' @param Mtype The diversity metric to use. Default is "M" (richness), allowed values are `r toString(.Mtype)`.
#' @param Atype The area metric to use. Default is "A" (area), allowed values are `r toString(.Atype)`.
#'
#' @return A fitted power model object from sars::sar_power
#' @export
#'
#' @examples
#' area <- 1:10
#' mutations <- c(50, 75, 100, 125, 150, 175, 200, 225, 250, 275)
#' data <- data.frame(A = area, M = mutations)
#' mar <- MARcalc(data, Mtype = "M", Atype = "A")
#' summary(mar)
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
