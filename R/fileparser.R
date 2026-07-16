# helper functions
.strip_ext <- function(filename, extensions) {
    bn <- basename(filename)
    matchext <- sapply(extensions, function(ext) grepl(paste0(ext, "$"), bn))
    stopifnot("File extention does not match allowed extintions" = sum(matchext) == 1)
    bn <- sub(paste0(extensions[matchext], "$"), "", bn)
    return(bn)
}

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
#' @param geno.fn Path to genotype file. Must be a txt/csv/tsv file (can be gzipped). Should contain a matrix where rows are SNPs and columns are samples, with values representing count of alternative alleles
#' @param samp.fn Optional path to sample ID file. Must be a single column file with no header
#' @param pos.fn Optional path to position file. Must have header with CHR/CHROM and POS columns
#' @param ploidy Integer specifying the ploidy level of the samples (default: 2)
#'
#' @return A margeno object containing the parsed genotype data and associated information
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with just genotype file
#' geno <- text_parser("genotypes.txt")
#'
#' # With sample IDs and positions
#' geno <- text_parser("genotypes.txt", "samples.txt", "positions.txt")
#'
#' # For haploid data
#' geno <- text_parser("genotypes.txt", ploidy = 1)
#' }
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
        sample.id <- seq_len(ncol(genotype))
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
        variant.id = seq_len(nrow(genotype)), # TODO not allow inputs for variant.id
        position = poslist[[2]],
        chromosome = poslist[[1]],
        genotype = genotype,
        ploidy = ploidy
    )

    return(margeno)
}

#' Convert VCF file to GDS format
#'
#' This function converts a VCF (Variant Call Format) file to GDS (Genomic Data Structure) format
#' using SeqArray package.
#'
#' @param vcf.fn Path to the input VCF file. Can be either .vcf or .vcf.gz format
#' @param gds.fn Optional. Path for the output GDS file. If NULL, will use the same name as vcf file
#'               with .gds extension
#' @param opengds Logical. If TRUE, opens and returns the GDS file handle. If FALSE, returns the
#'                path to the created GDS file. Default is FALSE
#'
#' @return If opengds=TRUE, returns an opened GDS file connection. If opengds=FALSE, returns the
#'         path to the created GDS file as a character string
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert VCF to GDS
#' gds_file <- vcf_parser("input.vcf")
#'
#' # Convert and open GDS file
#' gds_conn <- vcf_parser("input.vcf.gz", opengds = TRUE)
#' }
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
        variant.id = seq_len(n_var),
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
#' While coordinates do not need to be in "+proj=longlat +datum=WGS84" projection, the WGS84 projection was used in testing..
#' Sample IDs must be unique and in the same order as the Sample IDs provided in the genotype matrix.
#'
#' @param lonlat.fn Path to input file (txt/csv/tsv, can be gzipped) containing coordinates. No missing values allowed.
#'
#' @return A data frame containing sample IDs and their corresponding longitude/latitude coordinates.
#'         Returns error if coordinates contain missing values or incorrect number of columns.
#' @export
#'
#' @examples
#' \dontrun{
#' # Read coordinates from file
#' coords <- lonlat_parser("sample_locations.txt")
#' }
lonlat_parser <- function(lonlat.fn) {
    # check if lonlat.fn is a valid txt file
    txt.ext <- c(".txt", ".txt.gz", ".csv", ".csv.gz", ".tsv", ".tsv.gz")
    stopifnot("Invalid file type. Must be txt, csv, or tsv" = any(sapply(txt.ext, function(xx) grepl(xx, lonlat.fn))))
    # read txt file
    lonlatdf <- .read_lonlat(lonlat.fn)
    marmap <- marmaps(lonlatdf, mapres = NULL, mapcrs = "EPSG:8857")
    return(marmap)
}
