#' Arabidopsis thaliana 1001 Genomes example data
#'
#' A [genomaps()] object derived from the Arabidopsis 1001 Genomes Project.
#' A total of 1004 geo-referenced samples from the native continental range
#' (excluding the USA, Canary Islands or Japan samples) were kept, and 10,000
#' biallelic polymorphic SNPs were randomly sampled from chromosome 1.
#'
#' The raw inputs used to build this object are shipped in `inst/extdata`
#' (`1001g_genotypes.txt.gz`, `1001g_accessions.txt`, `1001g_chrpos.txt` and
#' `1001g_lonlat.txt`); see the examples of [text_parser()] and
#' [lonlat_parser()] for how to read raw data.
#'
#' @format ## `gm1001g`
#' An S3 `genomaps` object created using [genomaps()] function.
#' \describe{
#'   \item{geno}{An S3 `margeno` object with 10,000 variants x 1,004 diploid
#'     samples. See [margeno()] function.}
#'   \item{maps}{An S3 `marmaps` object with 1,004 samples on a 58 x 95 cell
#'     raster. See [marmaps()] function.}
#' }
#' @source <https://1001genomes.org/data/GMI-MPI/releases/v3.1/SNP_matrix_imputed_hdf5/>
#' @seealso [gmexp] for the simulated example data.
#' @examples
#' print(gm1001g)
"gm1001g"

#' Simulated MARtheory example data
#'
#' A [genomaps()] object derived from a simulated population
#' We simulated 50 diploid individuals sampled from a randomly mating population
#' growing exponentially, using the python package `msprime`. Coordinates are
#' randomly assigned to simulate a randomly distributed population.
#' For this population, the MAR can be predicted theoretically from the SFS.
#'
#' The simulated VCF and coordinates are shipped in `inst/extdata`
#' (`gmexp.vcf.gz` and `gmexp_lonlat.csv`); see the examples of [vcf_parser()]
#' and [lonlat_parser()] for how to read raw data.
#'
#' @format ## `gmexp`
#' An S3 `genomaps` object created using [genomaps()] function.
#' \describe{
#'   \item{geno}{An S3 `margeno` object with 885 variants x 50 diploid samples.
#'     See [margeno()] function.}
#'   \item{maps}{An S3 `marmaps` object with 50 samples. See [marmaps()] function.}
#' }
#' @source data-raw/gmexp_simulate.py
#' @seealso [gm1001g] for the empirical example data.
#' @examples
#' print(gmexp)
"gmexp"
