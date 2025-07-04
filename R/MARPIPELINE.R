#' Main pipeline function. Run MAR (mutations-area relationship) Analysis Pipeline
#'
#' @description
#' MARPIPELINE performs a complete MAR analysis workflow, including data loading,
#' building `genomaps` object, Species Abundance Distribution (SAD) fitting, Site Frequency
#' Spectrum (SFS) calculation, sampling of geographical areas to build MAR, and extinction simulations given MAR predictions.
#'
#' @param name Required. Character string. Base name for output files.
#' @param workdir Required. Character string. Working directory path where outputs will be saved. Need to be created before running the pipeline.
#' @param genofile Required. Character string. Path to input genomic data file. Full path to the file is required.
#' @param lonlatfile Required. Character string. Path to file containing longitude/latitude coordinates. Full path to the file is required.
#' @param samplefile Optional. Character string. Path to sample metadata file. Only used if filetype is 'text'.
#' @param posfile Optional. Character string. Path to variant position file. Only used if filetype is 'text'.
#' @param subsamplefile Optional. Character string. Path to file with sample IDs to subset.
#' @param subvariantfile Optional. Character string. Path to file with variant IDs to subset.
#' @param filetype Optional. Character string. Input genomic file format: 'text', 'vcf', or 'plink'. Default is 'text'.
#' @param randseed Optional. Integer. Random seed for reproducibility in steps mar and ext. Default is NULL.
#' @param option_geno List of genomic data options. Must provide complete list with all elements:
#'   \itemize{
#'     \item ploidy: Ploidy level of samples (default: 2)
#'     \item maxsnps: Maximum number of SNPs to analyze (default: 1000000)
#'   }
#' @param option_map List of mapping options. Must provide complete list with all elements:
#'   \itemize{
#'     \item mapres: Map resolution
#'     \item mapcrs: Coordinate reference system (default: WGS84)
#'   }
#' @param option_sadsfs List of SAD/SFS analysis options. Must provide complete list with all elements:
#'   \itemize{
#'     \item sad_models: SAD models to fit
#'     \item folded: Whether to use folded spectrum (default: TRUE)
#'   }
#' @param option_marext List of MAR/extinction analysis options. Must provide complete list with all elements:
#'   \itemize{
#'     \item scheme: Sampling schemes to use
#'     \item nrep: Number of replicates (default: 10)
#'     \item xfrac: Fraction for sampling (default: 0.01)
#'     \item quorum: Logical; require all MAR sampling boundaries to have at least one sample (default: TRUE)
#'     \item animate: Logical; create animation (default: FALSE)
#'   }
#' @param marsteps Character vector specifying which analysis steps to run:
#'   "data", "gm", "sfs", "mar", "ext", "plot"
#' @param saveobj Logical. Whether to save intermediate objects as .rda files.
#'
#' @note All option lists (option_geno, option_map, option_sadsfs, option_marext) must be provided
#' with their complete set of elements. Partial lists are not supported.
#'
#' @return Returns invisibly. Creates output files in working directory:
#'   \itemize{
#'     \item PDF plots for selected analysis steps if plot is in marsteps
#'     \item .rda files with analysis objects if saveobj=TRUE
#'   }
#'
#' @examples
#' \dontrun{
#' MARPIPELINE(
#'   name = "myanalysis",
#'   workdir = "path/to/output",
#'   genofile = "path/to/genotype.vcf",
#'   lonlatfile = "path/to/coordinates.txt",
#'   filetype = "vcf",
#'   option_geno = list(ploidy = 2, maxsnps = 1000000),
#'   option_map = list(mapres = NULL, mapcrs = "OGC:CRS84"),
#'   option_sadsfs = list(sad_models = .sad_models, folded = TRUE),
#'   option_marext = list(scheme = .MARsampling_schemes, nrep = 10,
#'                       xfrac = 0.01, quorum = TRUE, animate = FALSE),
#'   marsteps = c("data", "gm", "sfs", "mar", "plot")
#' )
#' }
#'
#' @export
MARPIPELINE <- function(name,
                        workdir,
                        genofile,
                        lonlatfile,
                        samplefile = NULL,
                        posfile = NULL,
                        subsamplefile = NULL,
                        subvariantfile = NULL,
                        filetype = c('text', 'vcf', 'plink'),
                        randseed = NULL,
                        option_geno = list(ploidy = 2, maxsnps = 1000000),
                        option_map = list(mapres = NULL, mapcrs = "OGC:CRS84"),
                        option_sadsfs = list(sad_models = .sad_models, folded = TRUE),
                        option_marext = list(scheme = .MARsampling_schemes, nrep = 10, xfrac = 0.01, quorum = TRUE, animate = FALSE),
                        marsteps = c("data", "gm", "sfs", "mar", "ext", "plot"),
                        saveobj = FALSE) {
# Define some variables --------------------------------------------------------
    options(warn = 1) # print warning as they occur
    # all potential output files (and objects)
    ofn <- list(
        data = c("genodata", "mapsdata"),
        gm = c("gm"),
        sfs = c("genosfs", "marsad", "sadsfs"),
        mar = c("mardflist", "marlist"),
        ext = c("extdflist", "extlist")
    )
    # combine the options into a single extra_file list
    extra_file <- list(samplefile = samplefile, posfile = posfile, subsample = subsamplefile, subvariant = subvariantfile)

# Check and file setup ---------------------------------------------------------
    message(paste0("MARPIPELINE starts at ", Sys.time(), "."))
    print(utils::sessionInfo())
    filetype <- match.arg(filetype)
    marsteps <- match.arg(marsteps, several.ok = TRUE)

    .valid_files(genofile, lonlatfile, extra_file, filetype, marsteps)
    .valid_options(filetype, option_geno, option_map, marsteps)
    outfile <- .valid_output(name, workdir, ofn)

    setwd(workdir) # change to working directory

# Load data --------------------------------------------------------------------
    if ("data" %in% marsteps) {
        message("MARPIPELINE loading genomic and geographic data ...")
        # genodata can be either margeno or a character pointing to gds file
        genodata <- switch(filetype,
            text = text_parser(genofile, samp.fn = extra_file$samplefile, pos.fn = extra_file$posfile,
                               ploidy = option_geno$ploidy),
            vcf = vcf_parser(genofile),
            plink = plink_parser(genofile)
        )

        # mapdata is a list of sample.id and matrix of lonlat
        mapsdata <- lonlat_parser(lonlatfile)
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
            mapsdata <- mapsdata[mapsdata[, 1] %in% subsamples, ]
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
            tempgeno0 <- .seqarray2margeno(tempgeno)
            SeqArray::seqClose(tempgeno)
            tempgeno <- tempgeno0
        }

        # construct genomaps object
        tempmaps <- marmaps(mapsdata, option_map$mapres, option_map$mapcrs)
        gm <- genomaps(tempgeno, tempmaps)

        # remove temporary objects
        rm(tempgeno, tempmaps)
    }

# Build SAR / SFS --------------------------------------------------------------
    if ("sfs" %in% marsteps) {
        message("MARPIPELINE fitting SAD models and calculating SFS ...")
        .required_objects("gm", ofn, outfile)
        genosfs <- sfs(AC = .get_AC(gm$geno),
                       N = length(gm$maps$sample.id),
                       ploidy = gm$geno$ploidy,
                       folded = option_sadsfs$folded)
        # fit SAD models
        marsad <- MARsad(gm = gm, sad_models = option_sadsfs$sad_models, folded = option_sadsfs$folded)
        message("SAD models AIC:")
        print(marsad$AICtabs)
        sadsfs <- .pipe_sadsfs(gm, marsad, genosfs, folded = option_sadsfs$folded)
        message("SAD model predictions compared with data in SFS:")
        print(sadsfs$statdf)
    }

# MARsampling ------------------------------------------------------------------
    if ("mar" %in% marsteps) {
        message("MARPIPELINE sampling the MAR distribution ...")
        .required_objects("gm", ofn, outfile) # requires output from "gm" step

        # Create sampling
        mardflist <- lapply(option_marext$scheme, function(scheme) {
            message(paste0("Sampling scheme: ", scheme))
            MARsampling(gm = gm, scheme = scheme, nrep = option_marext$nrep, xfrac = option_marext$xfrac,
                        quorum = option_marext$quorum, animate = option_marext$animate, myseed = randseed)
        })
        names(mardflist) <- option_marext$scheme

        # Calculate MAR
        marlist <- lapply(mardflist, .pipe_MARcalc)
        names(marlist) <- option_marext$scheme
    }

# MARextinction simulation -----------------------------------------------------
    if ("ext" %in% marsteps) {
        message("MARPIPELINE simulating extinction of distribution cells ...")
        .required_objects("gm", ofn, outfile) # requires output from "gm" step

        # Create extinction scheme
        extdflist <- lapply(option_marext$scheme, function(scheme) {
            message(paste0("Extinction scheme: ", scheme))
            MARextinction(gm = gm, scheme = scheme, nrep = option_marext$nrep, xfrac = option_marext$xfrac, animate = option_marext$animate, myseed = randseed)
        })
        names(extdflist) <- option_marext$scheme

        # Calculate MAR
        extlist <- lapply(extdflist, .pipe_MARcalc)
        names(extlist) <- option_marext$scheme
    }

# Plotting of given steps ------------------------------------------------------
    # To avoid issues with file versions, only plot if relevant step is run
    if ("plot" %in% marsteps) {
        message("MARPIPELINE plotting ...")
        if ("gm" %in% marsteps) {
            .pdf_plot(name, "gm", 6, 6)
            plot.marmaps(gm$maps)
            dev.off()
        }
        if ("sfs" %in% marsteps) {
            .pdf_plot(name, "sfs", 8, 9)
            .pipe_plot.marsadsfs(sadsfs, marsad$AICtabs)
            dev.off()
        }
        if ("mar" %in% marsteps) {
            lapply(option_marext$scheme, function(scheme) {
                .pdf_plot(name, paste0("mar.", scheme), 8, 9)
                .pipe_plot.marsamp(mardf = mardflist[[scheme]], mar = marlist[[scheme]])
                dev.off()
            })
        }
        if ("ext" %in% marsteps) {
            lapply(option_marext$scheme, function(scheme) {
                .pdf_plot(name, paste0("ext.", scheme), 8, 9)
                .pipe_plot.marextinct(extdf = extdflist[[scheme]], ext = extlist[[scheme]])
                dev.off()
            })
        }
    }

# Save data and exit -----------------------------------------------------------
    if (saveobj) {
        message("MARPIPELINE saving objects ...")
        for (ii in setdiff(marsteps, "plot")) {
            for (jj in seq_along(ofn[[ii]])) {
                save(list = ofn[[ii]][jj], file = outfile[[ii]][jj])
            }
        }
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

.valid_options <- function(filetype, option_geno, option_map, marsteps) {
    # check for seqarray
    if (filetype != "text") {
        if (!requireNamespace("SeqArray", quietly = TRUE)) {
            stop("SeqArray package is not installed. Run `BiocManager::install(\"SeqArray\")` first.")
        }
    }
    if (option_geno$ploidy != 2) {
        warning("Ploidy != 2. Be careful with the interpretation of the results.")
    }

    # check for step dependency
    if ("plot" %in% marsteps) {
        if (!any(marsteps %in% c("gm", "sfs", "mar", "ext"))) {
            stop("Plotting requires at least one of the following steps: gm, sfs, mar, ext")
        }
    }

    return(invisible())
}

.valid_output <- function(name, workdir, ofn) {
    # check if workdir exists
    stopifnot(dir.exists(workdir))
    # outfile does not have workdir as workdir will be set in pipeline
    outfile <- lapply(ofn, function(x) paste0(x, "_", name, ".rda"))
    return(outfile)
}

.pdf_plot <- function(name, plotname, ww, hh) {
    outname <- paste0(plotname, "_", name, ".pdf")
    pdf(file = outname, width = ww, height = hh)
}

.required_objects <- function(marstep, ofn, outfile) {
    for (ii in seq_along(ofn[[marstep]])) {
        obj <- ofn[[marstep]][ii]
        ofile <- outfile[[marstep]][ii]

        if(!exists(obj, envir = parent.frame())) {
            stopifnot(file.exists(ofile))
            message(paste0("loaded ", ofile, " as ", obj))
            load(ofile, envir = parent.frame())
        }
    }
    return(invisible())
}
