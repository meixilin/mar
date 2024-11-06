# helper functions
.strip_ext <- function(filename, extensions) {
    bn <- basename(filename)
    matchext <- sapply(extensions, function(ext) grepl(paste0(ext, "$"), bn))
    stopifnot(sum(matchext) == 1)
    bn <- sub(paste0(extensions[matchext], "$"), "", bn)
    return(bn)
}

.guess_delim <- function(firstline) {
    tab_count <- length(grep("\t", firstline, fixed = TRUE))
    comma_count <- length(grep(",", firstline, fixed = TRUE))
    stopifnot(any(c(tab_count, comma_count) > 0))
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
    if(!is.null(myheader)) {
        stopifnot(grepl(myheader, firstline, ignore.case = TRUE))
    }
    # guess delimiter
    delim <- .guess_delim(firstline)
    return(delim)
}

.read_table <- function(filename, header, sep) {
    con <- .open_txt(filename)
    tryCatch({
        df <- read.table(con, header = header, sep = sep, row.names = NULL, stringsAsFactors = FALSE)
    },
        finally = close(con)
    )
    return(df)
}

# only allow files without headers, just the genotype matrix. row is SNPs, column is samples.
.read_genotype <- function(geno.fn, het2hom) {
    # read first line to guess the format
    delim <- .firstline(geno.fn)
    # not allow row or column names
    df <- .read_table(geno.fn, header = FALSE, sep = delim)
    df <- as.matrix(df)
    dimnames(df) <- NULL
    # check that the data only contains 0,1,2
    .valid_genotype(df)
    # if convert het to hom
    if (het2hom) {
        # convert all values to 0/1
        message("converting heterozygotes to homozygotes")
        print(table(as.vector(df)))
        df = apply(df, 2, function(xx) ifelse(xx == 2, 1, xx))
        print(table(as.vector(df)))
    }
    .valid_genotype(df)
    return(df)
}

# header has to be CHR/CHROM POS
.read_pos <- function(pos.fn) {
    delim <- .firstline(pos.fn, "(CHROM|CHR)\\s*[,\t]\\s*POS")
    df <- .read_table(pos.fn, header = TRUE, sep = delim)
    return(list(df[[1]], df[[2]]))
}

# no header or delimiter allowed for samp.fn
.read_samp <- function(samp.fn) {
    df <- .read_table(samp.fn, header = FALSE, sep = "")
    stopifnot(ncol(df) == 1)
    return(df[[1]])
}

# header has to be id lon lat or id longitude latitude (in that order)
.read_lonlat <- function(lonlat.fn) {
    delim <- .firstline(lonlat.fn, "ID\\s*[,\t]\\s*(LON(GITUDE)?)\\s*[,\t]\\s*(LAT(ITUDE)?)")
    df <- .read_table(lonlat.fn, header = TRUE, sep = delim)
    .valid_lonlat(as.matrix(df[,2:3]))
    return(df)
}

# genotype txt/csv parser, with or without chromosome and position information
# default ploidy is 2
text_parser <- function(geno.fn, samp.fn = NULL, pos.fn = NULL, ploidy = 2, het2hom = TRUE) {
    if (ploidy != 2) {
        stop("ploidy other than 2 is better handled by SeqArray")
    }
    # check if geno.fn is a valid txt file
    txt.ext <- c(".txt", ".txt.gz", ".csv", ".csv.gz", ".tsv", ".tsv.gz")
    stopifnot(any(sapply(txt.ext, function(xx) grepl(xx, geno.fn))))
    # read txt file
    genotype <- .read_genotype(geno.fn)
    # read sample file if exists
    if (!is.null(samp.fn)) {
        sample.id <- .read_samp(samp.fn)
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

# vcf parser
vcf_parser <- function(vcf.fn, gds.fn = NULL, opengds = FALSE) {
    # check if vcf.fn is a valid vcf file
    vcf.ext <- c(".vcf", ".vcf.gz")
    stopifnot(any(sapply(vcf.ext, function(xx) grepl(xx, vcf.fn))))
    # assign name if gds.fn is not provided
    if (is.null(gds.fn)) {
        gds.fn <- paste0(.strip_ext(vcf.fn, vcf.ext), ".gds")
    }
    # create gds file
    SeqArray::seqVCF2GDS(vcf.fn, gds.fn)
    if (opengds) {
        f <- SeqArray::seqOpen(gds.fn)
        return(f)
    } else {
        return(gds.fn)
    }
}

# plink parser
plink_parser <- function(plink.fn, gds.fn = NULL, opengds = FALSE) {
    # add suffix to plink.fn
    bed.fn <- paste0(plink.fn, ".bed")
    fam.fn <- paste0(plink.fn, ".fam")
    bim.fn <- paste0(plink.fn, ".bim")
    # check if inputs exists
    stopifnot(all(sapply(c(bed.fn, fam.fn, bim.fn), file.exists)))
    # assign name if gds.fn is not provided
    if (is.null(gds.fn)) {
        gds.fn <- paste0(.strip_ext(bed.fn, ".bed"), ".gds")
    }
    # create gds file
    SeqArray::seqBED2GDS(bed.fn, fam.fn, bim.fn, gds.fn)
    if (opengds) {
        f <- SeqArray::seqOpen(gds.fn)
        return(f)
    } else {
        return(gds.fn)
    }
}

# lonlat file parser
lonlat_parser <- function(lonlat.fn, mapres = NULL, mapcrs = "+proj=longlat +datum=WGS84") {
    # check if lonlat.fn is a valid txt file
    txt.ext <- c(".txt", ".txt.gz", ".csv", ".csv.gz", ".tsv", ".tsv.gz")
    stopifnot(any(sapply(txt.ext, function(xx) grepl(xx, lonlat.fn))))
    # read txt file
    lonlatdf <- .read_lonlat(lonlat.fn)
    sample.id <- lonlatdf[[1]]
    lonlat <- as.matrix(lonlatdf[,2:3])

    # create marmaps object
    marmaps <- marmaps(
        sample.id = sample.id,
        lonlat = lonlat,
        mapres = mapres,
        mapcrs = mapcrs
    )
    return(marmaps)
}


