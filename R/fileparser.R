.guess_delim <- function(firstline) {
    tab_count <- length(grep("\t", firstline, fixed = TRUE))
    comma_count <- length(grep(",", firstline, fixed = TRUE))
    stopifnot("First line must contain a tab or comma delimiter" = any(c(tab_count, comma_count) > 0))
    delim <- ifelse(tab_count >= comma_count, "\t", ",")
    return(delim)
}

# open either txt or txt.gz file
.open_txt <- function(filename) {
    if (grepl(".gz$", filename)) {
        con <- gzfile(filename, "r")
    } else {
        con <- file(filename, "r")
    }
    return(con)
}

.firstline <- function(filename, myheader = NULL) {
    con <- .open_txt(filename)
    firstline <- readLines(con, n = 1)
    close(con)
    # check header
    if (!is.null(myheader)) {
        stopifnot("File header does not match the expected pattern" = grepl(myheader, firstline, ignore.case = TRUE)) # double check expected pattern
    }
    # guess delimiter
    delim <- .guess_delim(firstline)
    return(delim)
}

.read_table <- function(filename, header, sep) {
    con <- .open_txt(filename)
    tryCatch(
        {
            df <- utils::read.table(con, header = header, sep = sep, row.names = NULL, stringsAsFactors = FALSE)
        },
        finally = close(con)
    )
    return(df)
}

# only allow files without headers, just the genotype matrix. row is SNPs, column is samples, value is the amount of alternative alleles.
# Not allowing het2hom as the script has been updated to match ploidy.
# If want to use het2hom, manually convert the genotype matrix and set ploidy to one.
.read_genotype <- function(geno.fn, ploidy) {
    # read first line to guess the format
    delim <- .firstline(geno.fn)
    # not allow row or column names
    df <- .read_table(geno.fn, header = FALSE, sep = delim)
    df <- data.matrix(df)
    dimnames(df) <- NULL
    # check that the data only contains 0,1,...,ploidy
    # .valid_genotype(df, ploidy)
    # if (het2hom) {
    #     # convert all values to 0/1
    #     message("converting heterozygotes to homozygotes")
    #     print(table(as.vector(df)))
    #     df = apply(df, 2, function(xx) ifelse(xx > 1, 1, xx))
    #     print(table(as.vector(df)))
    # }
    .valid_genotype(df, ploidy)
    return(df)
}

# header has to be CHR/CHROM POS
.read_pos <- function(pos.fn) {
    delim <- .firstline(pos.fn, "(CHROM|CHR)\\s*[,\t]\\s*POS")
    df <- .read_table(pos.fn, header = TRUE, sep = delim)
    return(list(df[[1]], df[[2]]))
}

# no header or delimiter allowed for samp.fn (any single column file)
.read_column <- function(fn) {
    df <- .read_table(fn, header = FALSE, sep = "")
    stopifnot("samp.fn must be a single-column file with no header or delimiter" = ncol(df) == 1)
    return(df[[1]])
}

# header has to be id lon lat or id longitude latitude (in that order)
.read_lonlat <- function(lonlat.fn) {
    delim <- .firstline(lonlat.fn, "ID\\s*[,\t]\\s*(LON(GITUDE)?)\\s*[,\t]\\s*(LAT(ITUDE)?)")
    df <- .read_table(lonlat.fn, header = TRUE, sep = delim)
    .valid_lonlat(as.matrix(df[, 2:3]))
    return(df)
}

#' Parse text files containing genotype data
#'
#' Reads genotype data from text files (txt, csv, tsv) along with optional sample IDs and position information
#'
#' @param geno.fn Path to genotype file. The file must be a txt/csv/tsv file (can be gzipped), and contain a matrix where rows are SNPs and columns are samples, with values representing count of alternative alleles. No headers or row names are allowed.
#' @param samp.fn Optional path to sample ID file. The file must be a single column file with no header.
#' @param pos.fn Optional path to chromosome position file. The file must have header with CHR/CHROM and POS columns
#' @param ploidy Integer specifying the ploidy level of the samples (default: 2)
#'
#' @return A [margeno()] object containing the parsed genotype data and associated information
#' @seealso [vcf_parser()] to read genotypes from a VCF instead, and
#'   [lonlat_parser()] to read coordinates.
#' @export
#'
#' @examples
#' # read the genotypes behind the gm1001g example object
#' geno.fn <- system.file("extdata", "1001g_genotypes.txt.gz", package = "mar")
#' samp.fn <- system.file("extdata", "1001g_accessions.txt", package = "mar")
#' pos.fn <- system.file("extdata", "1001g_chrpos.txt", package = "mar")
#'
#' # the genotype matrix alone is enough
#' geno <- text_parser(geno.fn)
#'
#' # sample IDs and positions are optional but recommended
#' geno <- text_parser(geno.fn, samp.fn, pos.fn, ploidy = 2)
#' print(geno)
text_parser <- function(geno.fn, samp.fn = NULL, pos.fn = NULL, ploidy = 2) {
    # check if geno.fn is a valid txt file
    txt.ext <- c(".txt", ".txt.gz", ".csv", ".csv.gz", ".tsv", ".tsv.gz")
    stopifnot("Invalid file type. Must be txt, csv, or tsv" = any(sapply(txt.ext, function(xx) grepl(xx, geno.fn))))
    # read txt file
    genotype <- .read_genotype(geno.fn, ploidy)
    # read sample file if exists
    if (!is.null(samp.fn)) {
        sample.id <- .read_column(samp.fn)
    } else {
        sample.id <- NULL
    }
    # read chromosome and position file if exists
    if (!is.null(pos.fn)) {
        poslist <- .read_pos(pos.fn)
    } else {
        poslist <- list(NULL, NULL)
    }

    # create margeno object
    margeno <- margeno(
        sample.id = sample.id,
        variant.id = NULL, # TODO allow inputs for variant.id
        position = poslist[[2]],
        chromosome = poslist[[1]],
        genotype = genotype,
        ploidy = ploidy
    )

    return(margeno)
}

#' Parse a VCF file into a margeno object
#'
#' Reads the genotypes of a VCF (Variant Call Format) file into a [margeno()]
#' object. Only the `GT` field is used: each call is converted to a count of
#' alternative alleles, and missing calls (`./.`) become `NA`. Sample IDs,
#' chromosomes and positions are taken from the VCF itself. The whole file is
#' read into memory, so this is intended for a modest SNP panels used in MAR
#' analyses rather than for whole-genome all-sites VCFs.
#'
#' @param vcf.fn Path to the input VCF file. Can be either .vcf or .vcf.gz format
#' @param ploidy Integer specifying the ploidy level of the samples (default: 2).
#'
#' @return A [margeno()] object containing the parsed genotypes, sample IDs,
#'   chromosomes and positions.
#' @seealso [text_parser()] to read genotype matrices already in text form, and
#'   [lonlat_parser()] to read coordinates.
#' @export
#'
#' @examples
#' # simulated expanding population shipped with the package
#' vcf.fn <- system.file("extdata", "gmexp.vcf.gz", package = "mar")
#' geno <- vcf_parser(vcf.fn)
#' print(geno)
vcf_parser <- function(vcf.fn, ploidy = 2) {
    con <- if (grepl("\\.gz$", vcf.fn)) gzfile(vcf.fn) else file(vcf.fn)

    lines <- readLines(con)
    close(con)

    header_line <- grep("^#CHROM", lines, value = TRUE)
    data_lines <- lines[!grepl("^#", lines)]

    col_names <- sub("^#", "", strsplit(header_line, "\t")[[1]])
    fixed_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
    sample.id <- setdiff(col_names, fixed_cols)
    stopifnot("VCF can not be empty" = length(sample.id) > 0)

    split_lines <- strsplit(data_lines, "\t", fixed = TRUE)
    n_var <- length(split_lines)
    n_col <- length(col_names)

    fields <- matrix(unlist(split_lines, use.names = FALSE),
        nrow = n_var, ncol = n_col, byrow = TRUE,
        dimnames = list(NULL, col_names)
    )

    chromosome <- fields[, "CHROM"]
    position <- as.integer(fields[, "POS"])

    format_fields <- strsplit(fields[, "FORMAT"], ":")
    gt_idx <- vapply(format_fields, function(f) match("GT", f), integer(1))

    sample_col_idx <- match(sample.id, col_names)
    genotype <- matrix(NA_real_, nrow = n_var, ncol = length(sample.id))

    for (j in seq_along(sample.id)) {
        cell_split <- strsplit(fields[, sample_col_idx[j]], ":", fixed = TRUE)
        gt_strings <- mapply(function(cell, idx) cell[idx], cell_split, gt_idx)
        genotype[, j] <- vapply(gt_strings, .gt_to_dosage, numeric(1))
    }

    .valid_genotype(genotype, ploidy)

    margeno <- margeno(
        sample.id = sample.id,
        variant.id = NULL, # TODO: allow VCF variant id inputs
        position = position,
        chromosome = chromosome,
        genotype = genotype,
        ploidy = ploidy
    )

    return(margeno)
}

# Converts a GT string ("0/1", "1|1", "./.", etc.) to ALT-allele dosage.
.gt_to_dosage <- function(gt) {
    if (is.na(gt) || gt %in% c(".", "./.", ".|.")) {
        return(NA_real_)
    }
    alleles <- strsplit(gt, "[/|]")[[1]]
    if (any(alleles == ".")) {
        return(NA_real_)
    }
    sum(alleles != "0")
}

#' Parse longitude/latitude coordinates from text file
#'
#' Reads a file containing sample IDs with their corresponding longitude and latitude coordinates.
#' The input file must have a header with columns: ID, LON/LONGITUDE, LAT/LATITUDE (in that order).
#' Sample IDs must be unique and in the same order as the Sample IDs provided in the genotype matrix.
#' If Sample IDs were not provided in [text_parser()], sample IDs should be 1,2,3, ... , N.
#' @param lonlat.fn Path to input file (txt/csv/tsv, can be gzipped) containing coordinates. No missing values allowed.
#' @inheritDotParams marmaps mapcrs mapres incrs
#'
#' @return A marmaps object. See [marmaps()].
#'         Returns error if coordinates contain missing values or incorrect number of columns.
#' @seealso [marmaps()], and [text_parser()] / [vcf_parser()] for the matching
#'   genotypes.
#' @export
#'
#' @examples
#' # read the coordinates behind the gm1001g example object in longitude-latitude formats (EPSG:4326)
#' # and reprojected onto the Equal Earth Greenwich (EPSG: 8857)
#' lonlat.fn <- system.file("extdata", "1001g_lonlat.txt", package = "mar")
#' maps <- lonlat_parser(lonlat.fn, mapcrs = "EPSG:8857")
#' print(maps)
#' plot(maps)
#'
#' # simulated locations without a coordinate reference system, so no projection is applied
#' lonlat.fn <- system.file("extdata", "gmexp_lonlat.csv", package = "mar")
#' lonlat_parser(lonlat.fn, incrs = "", mapcrs = "")
lonlat_parser <- function(lonlat.fn, ...) {
    # check if lonlat.fn is a valid txt file
    txt.ext <- c(".txt", ".txt.gz", ".csv", ".csv.gz", ".tsv", ".tsv.gz")
    stopifnot("Invalid file type. Must be txt, csv, or tsv" = any(sapply(txt.ext, function(xx) grepl(xx, lonlat.fn))))
    # read txt file
    lonlatdf <- .read_lonlat(lonlat.fn)
    marmap <- marmaps(lonlatdf, ...)
    return(marmap)
}
