#' Print Method for margeno Objects
#'
#' Displays a summary of genetic data including sample counts, genomic sites,
#' and a preview of the genotype matrix
#'
#' @param x An object of class "margeno" containing genetic data
#' @param ... Additional arguments passed to print
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' print(gm1001g$geno)

print.margeno <- function(x, ...) {
    cat("margeno object\n")
    cat("    number of samples: ", length(x$sample.id), "\n")
    cat("    number of genomic sites: ", length(x$variant.id), "\n")
    cat("    ploidy: ", x$ploidy, "\n")
    cat("\n")
    cat("head of sampleid, variantid, position, chromosome, genotype:\n")
    AC <- rowSums(head(x$genotype))
    df <- cbind(head(data.frame(variant.id = x$variant.id,
                                position = .null_check(x$position),
                                chromosome = .null_check(x$chromosome))),
                AC,
                x$genotype[1:min(6, nrow(x$genotype)),1:min(6, ncol(x$genotype))])
    colnames(df)[5:length(colnames(df))] <- head(x$sample.id)
    print(df)
    return(invisible())
}

#' Print Method for marmaps Objects
#'
#' Displays a summary of spatial mapping data including sample counts,
#' geographic ranges, and a preview of sample locations
#'
#' @param x An object of class "marmaps" containing mapping information
#' @param ... Additional arguments passed to print
#'
#' @return Invisibly returns NULL
#' @export
#'
#' @examples
#' # Create and print a marmaps object
#' print(gm1001g$maps)
print.marmaps <- function(x, ...) {
    cat("marmaps object\n")
    cat("    number of samples: ", length(x$sample.id), "\n")
    cat("    longitude range: [", min(x$lonlat[,1]), ", ", max(x$lonlat[,1]), "]\n", sep = "")
    cat("    latitude range: [", min(x$lonlat[,2]), ", ", max(x$lonlat[,2]), "]\n", sep = "")
    cat("\n")
    cat("samplemap raster layer:\n")
    print(x$samplemap)
    cat("head of sampleid and lonlat:\n")
    print(head(data.frame(sample.id = x$sample.id, longitude = x$lonlat[,1], latitude = x$lonlat[,2])))
    return(invisible())
}

.null_check <- function(x) {
    if(is.null(x)) {
        return(NA)
    } else {
        return(x)
    }
}

