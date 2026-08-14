# calculate MAR relationship
.Mtype <- c("M", "E", "thetaw", "thetapi")
.Atype <- c("A", "Asq", "N")

#' Calculate the mutations-area relationship (MAR) using the power-law model
#'
#' This function calculates the MAR using the power-law model from the provided data frame.
#' The data frame should contain columns for the diversity metric (e.g. number of mutation, nucleotide diversity)
#' and the area metric (e.g. occupied area, total squared area).
#'
#' @param mardf A data frame with columns of diversity metric and area.
#'        Output from \link{MARsampling}, \link{MARextinction} or \link{MARtheory}.
#' @param Mtype The diversity metric to use. Default is "M" (number of mutation), allowed values are `r toString(.Mtype)`.
#' @param Atype The area metric to use. Default is "A" (area), allowed values are `r toString(.Atype)`.
#'
#' @return A fitted power model object from [sars::sar_power()]
#' @seealso [MARsampling()], [MARextinction()] and [MARtheory()] for generating the input `mardf`.
#' @export
#'
#' @examples
#' mardf <- MARsampling(gm1001g, nrep = 2, xfrac = 0.1, myseed = 42)
#' mar <- MARcalc(mardf, Mtype = "M", Atype = "A")
#' summary(mar)
#'
#' # the same power law model fitted against sample size rather than area
#' MARcalc(mardf, Mtype = "M", Atype = "N")
MARcalc <- function(mardf, Mtype = .Mtype, Atype = .Atype) {
    Mtype <- match.arg(Mtype)
    Atype <- match.arg(Atype)
    # remove NA or zero data
    tmpdf <- mardf[, c(Atype, Mtype)]
    tmpdf <- tmpdf[(tmpdf[, Mtype] > 0 & !is.na(tmpdf[, Mtype])), ]
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

    outdf <- tryCatch(
        {
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
        }
    )

    return(outdf)
}
