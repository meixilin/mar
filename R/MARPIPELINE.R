#' MARPIPELINE wrapper function
#'
#' @param name
#' @param workdir
#' @param genofile
#' @param lonlatfile
#' @param extra_file
#' @param filetype
#' @param option_geno
#' @param option_map
#' @param option_marext
#' @param marsteps
#' @param saveobj
#'
#' @return
#' @export
#'
#' @examples
MARPIPELINE <- function(name,
                        workdir,
                        genofile,
                        lonlatfile,
                        extra_file = list(samplefile = NULL, posfile = NULL, subsample = NULL, subvariant = NULL),
                        filetype = c('text', 'vcf', 'plink'),
                        option_geno = list(ploidy = 2, maxsnps = 1000000),
                        option_map = list(mapres = NULL, mapcrs = "+proj=longlat +datum=WGS84"),
                        option_sadsfs = list(sad_models = .sad_models),
                        option_marext = list(scheme = .MARsampling_schemes, nrep = 10, quorum = FALSE, animate = FALSE, myseed = NULL),
                        marsteps = c("data", "gm", "sfs", "mar", "ext", "plot"),
                        saveobj = FALSE) {
# Define some variables --------------------------------------------------------
    options(warn = 1) # print warning as they occur
    # all potential output files (and objects)
    ofn <- list(
        data = c("genodata", "mapsdata"),
        gm = c("gm"),
        sfs = c("sfslist"),
        mar = c("mardflist", "marlist"),
        ext = c("extdflist", "extlist")
    )

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
        message("MARPIPELINE fitting SAD models and calculating SFS ...")
        .required_objects("gm", ofn, outfile)
        AC <- .get_AC(gm$geno) # allele counts used for SAD and SFS
        sfslist <- MARsadsfs(AC, N = length(gm$maps$sample.id), ploidy = option_geno$ploidy, sad_models = option_sadsfs$sad_models)
        message("SAD models AIC:")
        print(sfslist$AICtabs)
        message("SAD model predictions compared with data in SFS:")
        print(sfslist$statdf)
    }

# MARsampling ------------------------------------------------------------------
    if ("mar" %in% marsteps) {
        message("MARPIPELINE sampling the MAR distribution ...")
        .required_objects("gm", ofn, outfile) # requires output from "gm" step

        # Create sampling
        mardflist <- lapply(option_marext$scheme, function(scheme) {
            message(paste0("Sampling scheme: ", scheme))
            MARsampling(gm = gm, scheme = scheme, nrep = option_marext$nrep, xfrac = option_marext$xfrac,
                        quorum = option_marext$quorum, animate = option_marext$animate, myseed = option_marext$myseed)
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
            MARextinction(gm = gm, scheme = scheme, nrep = option_marext$nrep, xfrac = option_marext$xfrac, animate = option_marext$animate, myseed = option_marext$myseed)
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
            .pipe_plot.marsadsfs(sfslist)
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
