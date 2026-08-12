#' mar: Fitting Mutations-Area Relationships using population
#' genomic samples of a species
#'
#' mar package provides functions to approximate the genetic diversity
#' loss in a species from habitat reduction.
#'
#' @section Key functions:
#' \describe{
#'   \item{\code{\link{text_parser}}, \code{\link{vcf_parser}}, \code{\link{lonlat_parser}}}{read genotypes and sample coordinates from files}
#'   \item{\code{\link{genomaps}}}{pairs the genotypes and the sample map into the core data object for downstream analyses}
#'   \item{\code{\link{MARsampling}}}{samples genetic diversity across varying spatial scales}
#'   \item{\code{\link{MARcalc}}}{fit power-law model}
#'   \item{\code{\link{MARextinction}}}{projects the expected loss of genetic diversity under hypothetical area reduction scenarios}
#'   \item{\code{\link{MARtheory}}}{predicts the sampling relationship expected from the site frequency spectrum alone}
#'   \item{\code{\link{sfs}}, \code{\link{expsfs}}}{the observed site frequency spectrum and its neutral expectation}
#' }
#'
#' @section Example data:
#' \code{\link{gm1001g}} (Arabidopsis 1001 Genomes) and \code{\link{gmexp}}
#' (simulated population) are ready-made `genomaps` objects used
#' throughout the examples and vignettes.
#'
#' @seealso Useful links:
#' \itemize{
#'   \item{Package repository: \url{https://github.com/meixilin/mar}}
#'   \item{Package paper: \url{https://www.biorxiv.org/content/10.1101/2025.09.09.675155v1.full}}
#' }
#'
#' @docType package
#' @keywords internal
#' @name mar
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
