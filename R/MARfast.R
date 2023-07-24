
#' Title Faster implementation of `MARsampling`
#'
#' @param coord data.frame with 3 columns,
#' @param geno
#' @param samples
#' @param res
#' @param quorum
#' @param minN
#'
#' @return
#' @export
#'
#' @examples
#' library(mar)
#' data("coords-Arabidopsis_thaliana")
#' data("geno-Arabidopsis_thaliana")
#' ### Run MAR
#' MAR = MARfastGT(coord = coords, geno = genomes, res = 0.01)
#' ### Plot MAR
#' plot(MAR$A_sq, MAR$M)
#' plot(log10(MAR$A_sq), log10(MAR$M))

MARfastGT <- function(coord, geno, samples = 100, res = 0.05, quorum = T, minN = 1) {

    # calculate lon/lat range of study area , add 50% of range to each edge
    lonrange <- range(coord[, 1]) + c(-1, 1) * diff(range(coord[, 1])) * 0.50 # range coords ± 50% of study area
    latrange <- range(coord[, 2]) + c(-1, 1) * diff(range(coord[, 2])) * 0.50 # range coords ± 50% of study area

    # set min-max dimension of squares to compute MAR
    minD <- res # minimal dimension of squares: map resolution
    maxD <- min(diff(lonrange), diff(latrange)) / 2 # maximal dimension of squares: half of smallest side of study area

    # create output container
    out <- c()

    # compute boxes, based on number of samples set as input ()
    i <- 1
    while (i <= samples) {

        # randomly sample a side size to compute the box, from a uniform distribution
        sqd <- runif(1, minD, maxD)

        ### randomly pick a square center, from a random distribution
        x <- runif(1, min(lonrange) + sqd, max(lonrange) - sqd)
        y <- runif(1, min(latrange) + sqd, max(latrange) - sqd)

        ## set square boundaries
        x1 <- x - sqd
        x2 <- x + sqd
        y1 <- y - sqd
        y2 <- y + sqd

        # find samples inside boundaries
        SAM.in <- which(coord[, 1] > x1 & coord[, 1] < x2 & coord[, 2] > y1 & coord[, 2] < y2)

        # calculate number of samples in square
        N <- length(SAM.in)

        # compute mutation count
        M <- sum(apply(geno[SAM.in, , drop = F], 2, function(x) { # for every column of SNP matrix...
            # count major and minor alleles
            ac <- c((sum(x == 0, na.rm = T) * 2) + sum(x == 1, na.rm = T), (sum(x == 2, na.rm = T) * 2) + sum(x == 1, na.rm = T))
            # return: if both alleles are >0, there is a mutation
            return(sum(ac > 0) == 2)
        }))

        ### compute area of square
        A_sq <- (x2 - x1) * (y2 - y1)

        # End iteration based on input : with or without quorum

        # with quorum: all samples must have a value
        if (quorum) {
            if (N >= minN & M > 0) { # if there are samples and mutations
                # add to output
                out <- rbind(out, c("M" = M, "A_sq" = A_sq, "N" = N))
                # update i
                i <- i + 1
            } else { # if there are no samples -> square is empty
                "" # do nothing: new spatial sample to be picked
            }
        }

        # without quorum: if sample is empty, don't keep the sample, but do not collect a new one in replacement
        if (quorum == F) {
            if (N >= minN & M > 0) { # if there are samples
                # add to output
                out <- rbind(out, c("M" = M, "A_sq" = A_sq, "N" = N))
            }
            # update i in any case
            i <- i + 1
        }
    }

    out <- as.data.frame(out)
    return(out)
}
