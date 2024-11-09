# genofile is of length 1
MARPIPELINE <- function(name,
                        workdir,
                        genofile,
                        lonlatfile,
                        extra_file = list(samplefile = NULL, posfile = NULL, subsample = NULL, subvariant = NULL),
                        filetype = c('text', 'vcf', 'plink'),
                        option_geno = list(ploidy = 2, het2hom = FALSE, maxsnps = 1000000),
                        option_map = list(mapres = NULL, mapcrs = "+proj=longlat +datum=WGS84"),
                        option_mar = list(scheme = .MARsampling_schemes, nrep = 10, quorum = FALSE, animate = FALSE, myseed = NULL),
                        marsteps = c("data", "gm", "sfs", "mar", "ext"),
                        saveobj = FALSE) {
# Define some variables --------------------------------------------------------
    # all potential output files (and objects)
    ofn <- list(
        data = c("genodata", "mapsdata"),
        gm = c("gm"),
        sfs = c("sfs"),
        mar = c("mardflist", "marlist"),
        ext = c("extdflist", "extlist")
    )

# Check and file setup ---------------------------------------------------------
    message(paste0("MARPIPELINE starts at ", Sys.time(), "."))
    print(sessionInfo())
    filetype <- match.arg(filetype)
    marsteps <- match.arg(marsteps, several.ok = TRUE)

    .valid_files(genofile, lonlatfile, extra_file, filetype, marsteps)
    .valid_options(filetype, option_geno, option_map)
    outfile <- .valid_output(name, workdir, marsteps, ofn)

    setwd(workdir) # change to working directory

# Load data --------------------------------------------------------------------
    if ("data" %in% marsteps) {
        message("MARPIPELINE loading genomic and geographic data ...")
        # genodata can be either margeno or a character pointing to gds file
        genodata <- switch(filetype,
            text = text_parser(genofile, samp.fn = extra_file$samplefile, pos.fn = extra_file$posfile,
                               ploidy = option_geno$ploidy, het2hom = option_geno$het2hom),
            vcf = vcf_parser(genofile),
            plink = plink_parser(genofile)
        )

        # mapdata is a list of sample.id and matrix of lonlat
        mapsdata <- lonlat_parser(lonlatfile, mapres = option_map$mapres, mapcrs = option_map$mapcrs)
    }

# Filter data and construct genomaps object ------------------------------------
    if ("gm" %in% marsteps) {
        message("MARPIPELINE constructing genomaps object ...")
        .required_objects("data", ofn, outfile) # requires output from "data" step
        # load SeqArray if needed
        if (filetype != "text") {
            stopifnot(file.exists(genodata))
            tempgeno <- SeqArray::seqOpen(genodata)
        } else {
            tempgeno <- genodata
        }
        # subset samples if needed
        if (!is.null(extra_file$subsample)) {
            message("Subsetting samples ...")
            subsamples <- .read_column(extra_file$subsample)
            tempgeno <- .filter_genosample(tempgeno, subsamples)
            mapsdata <- mapsdata[mapsdata$id %in% subsamples, ]
        }
        # subset variants if needed
        if (!is.null(extra_file$subvariant)) {
            message("Subsetting variants ...")
            subvariants <- .read_column(extra_file$subvariant)
            tempgeno <- .filter_genovariant(tempgeno, subvariants)
        }
        # prune variants if number of variants is too large
        if (option_geno$maxsnps < .get_genodata(tempgeno, "num.variant")) {
            message("Too many variants, pruning ...")
            subvariants <- sort(sample(.get_genodata(tempgeno, "variant.id"), option_geno$maxsnps))
            tempgeno <- .filter_genovariant(tempgeno, subvariants)
        }

        # convert SeqArray object to margeno object after filtering
        if (filetype != "text") {
            tempgeno <- .seqarray2margeno(tempgeno)
            SeqArray::seqClose(genodata)
        }

        # construct genomaps object
        tempmaps <- marmaps(mapsdata, option_map$mapres, option_map$mapcrs)
        gm <- genomaps(tempgeno, tempmaps)

        # remove temporary objects
        rm(tempgeno, tempmaps)
    }

# Build SAR / SFS --------------------------------------------------------------
    if ("sfs" %in% marsteps) {
        message("site frequency spectrum ...")
        # freq <- data.table::fread(freqfile)
        # M <- (ceiling(freq$MAF * freq$NCHROBS))
        # M <- M[!is.na(M) & M != 0]
        # Mrad <- rad(M)

        # models <- c("weibull", "lnorm", "ls", "bs", "geom", "mzsm")
        # fitted_models <- lapply(models, function(m) fitsad(M, m))
        # names(fitted_models) <- models
        # fitted_models$sfs <- predSFS(M)

        # tableaic <- do.call(AICtab, fitted_models[names(fitted_models) != "sfs"])
        # liksfs <- likelihoodSFSfold(M)

        # if (saveobjects) {
        #     save(file = filetable, object = tableaic)
        #     save(file = fileliksfs, object = liksfs)
        # }
    }

# MARsampling ------------------------------------------------------------------
    if ("mar" %in% marsteps) {
        message("MARPIPELINE sampling the MAR distribution ...")
        .required_objects("gm", ofn, outfile) # requires output from "gm" step

        # Create sampling
        mardflist <- lapply(option_mar$scheme, function(scheme) {
            message(paste0("Sampling scheme: ", scheme))
            MARsampling(gm = gm, scheme = scheme, nrep = option_mar$nrep, xfrac = option_mar$xfrac,
                        quorum = option_mar$quorum, animate = option_mar$animate, myseed = option_mar$myseed)
        })
        names(mardflist) <- option_mar$scheme

        # Calculate MAR
        marlist <- lapply(mardflist, MARcalc_all)
        names(marlist) <- option_mar$scheme
    }

# MARextinction simulation -----------------------------------------------------
    if ("ext" %in% marsteps) {
        message("MARPIPELINE simulating extinction of distribution cells ...")
        .required_objects("gm", ofn, outfile) # requires output from "gm" step

        # Create extinction scheme
        extdflist <- lapply(option_mar$scheme, function(scheme) {
            message(paste0("Extinction scheme: ", scheme))
            MARextinction(gm = gm, scheme = scheme, nrep = option_mar$nrep, xfrac = option_mar$xfrac, animate = option_mar$animate, myseed = option_mar$myseed)
        })
        names(extdflist) <- option_mar$scheme

        # Calculate MAR
        extlist <- lapply(extdflist, MARcalc_all)
        names(extlist) <- option_mar$scheme
    }

# Save data and exit -----------------------------------------------------------
    if (saveobj) {
        message("MARPIPELINE saving objects ...")
        .save_objects(marsteps, ofn, outfile)
    }

    message(paste0("MARPIPELINE successfully ends at ", Sys.time(), "."))
    return(invisible())
}

# functions used in MARPIPELINE ------------------------------------------------
# allows NULL in all file inputs and will not check
.valid_files <- function(genofile, lonlatfile, extra_file, filetype, marsteps) {
    if ("data" %in% marsteps) {
        if (filetype == "plink") {
            # not checking for NULL as ".bed" file will not exists
            genofile <- paste0(genofile, c(".bed", ".fam", ".bim"))
        }
        allfiles <- c(genofile, lonlatfile, unlist(extra_file))
        stopifnot(all(file.exists(allfiles)))
    }
    return(invisible())
}

.valid_options <- function(filetype, option_geno, option_map) {
    # check for seqarray
    if (filetype != "text") {
        if (!requireNamespace("SeqArray", quietly = TRUE)) {
            stop("SeqArray package is not installed. Run `BiocManager::install(\"SeqArray\")` first.")
        }
    }
    if (option_geno$het2hom) {
        warning("het2hom = TRUE. Converting heterozygotes to homozygotes.")
    }
    if (option_geno$ploidy != 2) {
        warning("ploidy != 2. Be careful with the interpretation of the results.")
    }
    return(invisible())
}

.valid_output <- function(name, workdir, marsteps, ofn) {
    # check if workdir exists
    stopifnot(dir.exists(workdir))
    # check if required files exist if not running all steps
    if (any(c("sfs", "mar", "ext") %in% marsteps) & !("gm" %in% marsteps)) {
        stopifnot(file.exists(paste0(workdir, "/genomaps_", name, ".rda")))
    }
    # outfile does not have workdir as workdir will be set in pipeline
    outfile <- lapply(ofn, function(x) paste0(x, "_", name, ".rda"))
    return(outfile)
}

.required_objects <- function(marsteps, ofn, outfile) {
    for (ii in seq_along(ofn[[marsteps]])) {
        obj <- ofn[[marsteps]][ii]
        ofile <- outfile[[marsteps]][ii]

        if(!exists(obj, envir = parent.frame())) {
            message(paste0("loaded ", ofile, " as ", obj))
            load(ofile, envir = parent.frame())
        }

    }
    return(invisible())
}

.save_objects <- function(marsteps, ofn, outfile) {
    for (ii in marsteps) {
        for (jj in seq_along(ofn[[ii]])) {
            save(list = ofn[[ii]][jj], file = outfile[[ii]][jj])
        }
    }
}
