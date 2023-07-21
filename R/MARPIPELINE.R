MARPIPELINE <- function(name = "amaranthus",
                        # ****************** #
                        # define parameters
                        nsnps = 10000,
                        maxsnps = 100000,
                        nsamplesMAR = 10,
                        mapresolution = 1,
                        REDO_FREQ = TRUE,
                        REDO_GENOME = TRUE,
                        REDO_GENEMAPS = TRUE,
                        REDO_MARSAMPLING = TRUE,
                        REDO_MAR = TRUE,
                        REDO_EXTINCTION = TRUE,
                        # ****************** #
                        # path to store
                        mypath = "tmpobjects",
                        # files
                        freqfile = "../data/1001g/1001gbi.chr1.frq",
                        famfile = "../data/1001g/1001gbi.chr1.fam",
                        mapfile = "../data/1001g/1001gbi.chr1.map",
                        bedfile = "../data/1001g/1001gbi.chr1.bed",
                        latlonfile = "../data/1001g/1001gbi.chr1.famlatlon.csv",
                        # ****************** #
                        debug = FALSE,
                        saveobjects = TRUE) {
    # load all packages needed
    require(mar)
    name <- name
    message("run name ", name)
    # Create files needed
    filetable <- paste0(mypath, "/tableaic-", name, ".rda")
    fileliksfs <- paste0(mypath, "/liksfs-", name, ".rda")
    filecoords <- paste0(mypath, "/coords-", name, ".rda")
    filegenome <- paste0(mypath, "/genome-", name, ".rda")
    filegenemaps <- paste0(mypath, "/genemaps-", name, ".rda")
    filemares <- paste0(mypath, "/mares-", name, ".rda")
    filemar <- paste0(mypath, "/mar-", name, ".rda")
    fileemar <- paste0(mypath, "/emar-", name, ".rda")
    filemarz <- paste0(mypath, "/marz-", name, ".rda")
    fileemarz <- paste0(mypath, "/emarz-", name, ".rda")
    fileextinctionmar <- paste0(mypath, "/extinctionmar-", name, ".rda")
    fileextinctionsim <- paste0(mypath, "/extinctionsim-", name, ".rda")
    fileextinctionsimR2 <- paste0(mypath, "/extinctionsimR2-", name, ".rda")
    ## freq ##############################################################################
    if (REDO_FREQ) {
        message("site frequency spectrum ...")
        freq <- data.table::fread(freqfile) # %>% head(nsnps) # uncomment for all data
        M <- (ceiling(freq$MAF * freq$NCHROBS))
        M[is.na(M)] <- 0
        M <- M[M != 0]
        Mrad <- rad(M)
        # fit models
        # m.g <- fitsad(M, "gamma",)
        m.w <- fitsad(M, "weibull")
        m.ln <- fitsad(M, "lnorm")
        m.ls <- fitsad(M, "ls")
        m.bs <- fitsad(M, "bs")
        m.geom <- fitsad(M, "geom")
        m.mzsm <- fitsad(M, sad = "mzsm")
        # m.poi <- fitsad(M, "poilog")
        m.sfs <- predSFS(M)

        tableaic <- AICtab(
            m.bs,
            m.mzsm,
            # m.poi,
            m.ln,
            m.ls,
            # m.g,
            m.w,
            m.geom
        )
        # tableR2<-list(
        #   "bs"=pseudoR2(Mrad$abund,radpred(m.bs)$abund),
        #   "mzsm"=pseudoR2(Mrad$abund,radpred(m.mzsm)$abund),
        #   "lognormal"=pseudoR2(Mrad$abund,radpred(m.ln)$abund),
        #   "logseries"=pseudoR2(Mrad$abund,radpred(m.ls)$abund),
        #   "gamma"=pseudoR2(Mrad$abund,radpred(m.g)$abund),
        #   "weibull"=pseudoR2(Mrad$abund,radpred(m.w)$abund),
        #   "geom"=pseudoR2(Mrad$abund,radpred(m.geom)$abund),
        #   "sfs"=pseudoR2(Mrad$abund,rad(m.sfs)$abund)
        # )
        liksfs <- likelihoodSFSfold(M)
        if (saveobjects) save(file = filetable, object = tableaic)
        if (saveobjects) save(file = fileliksfs, object = liksfs)
    }
    ## genome ##############################################################################
    if (REDO_GENOME) {
        message("loading genomic and geographic data ...")
        stopifnot(nsnps < maxsnps)
        ## Genomic data
        fam <- data.table::fread(famfile)
        map <- data.table::fread(mapfile)
        if (nrow(map) < maxsnps) maxsnps <- nrow(map)
        bed <- readbed(bedfile,
            N = nrow(fam), p = nrow(map),
            myrows = 1:nrow(fam), mycols = 1:maxsnps
        )
        ## Latitude and Longitude information
        latlon <- utils::read.csv(latlonfile)
        stopifnot(all(c("id", "latitude", "longitude") %in% colnames(latlon)))

        # Design objects for later functions
        if (nrow(map) < nsnps) nsnps <- nrow(map)
        # set.seed(1)
        genomes <- bed[
            # remove those that do not have long/lat, previously parsed to match fam
            !is.na(latlon$latitude) & !is.na(latlon$longitude),
            # use only limited SNPs for computational efficiency
            sample(1:ncol(bed), nsnps)
        ]
        coords <- cbind(latlon$longitude, latlon$latitude)
        coords <- coords[!is.na(latlon$latitude) & !is.na(latlon$longitude), ]
        # save objects
        if (saveobjects) save(file = filecoords, coords)
        if (saveobjects) save(file = filegenome, genomes)
        message("...done")
    }
    ## genemaps ##############################################################################
    # Create sampling
    if (REDO_GENEMAPS) {
        message("creating gene maps ...")
        load(file = filecoords)
        load(file = filegenome)
        # Create map object
        genemaps <- create_gene_maps(coords, genomes, by = mapresolution)
        # genemaps<-create_gene_maps(coords, genomes, by=1) # for debug
        # save object
        if (saveobjects) save(file = filegenemaps, genemaps)
        message("...done")
    }
    ## mars##############################################################################
    if (REDO_MARSAMPLING) {
        message("sampling the distribution ...")
        load(file = filegenemaps)
        # Create sampling
        mares <- MARsampling(genemaps, scheme = "random", samples = nsamplesMAR)
        # save object
        if (saveobjects) save(file = filemares, mares)
        message("...done")
    }
    #****************************************************************************** #
    if (REDO_MAR) {
        message("building the mutations area relationship ...")
        load(file = filemares)
        # head(mares); tail(mares)
        # Build MAR
        # tmpmar<-dplyr::select(filter(mares,M>0,asub>0), asub,M)
        tmpmar <- dplyr::select(filter(mares, M > 0, a > 0), a, M)
        # changed as of Aug 6 2021 to the total area, it seems less sensitive
        if (length(dim(tmpmar)[1] == 0 | unique(tmpmar$M)) < 3) {
            message("same mutations in every location! cannot calculate MAR")
            mar <- NA
            marz <- NA
        } else {
            mar <- sar_power(tmpmar)
            marz <- marcoef(mar)
            # if(debug) print(mar)
        }
        # Build EMAR
        tmpemar <- dplyr::select(filter(mares, E > 0, asub > 0), asub, E)
        if (dim(tmpemar)[1] == 0 | length(unique(tmpmar$E)) < 3) {
            message("no endemic mutations! cannot calculate EMAR")
            emar <- NA
            emarz <- NA
        } else {
            emar <- sar_power(tmpemar)
            emarz <- marcoef(emar)
        }

        save(file = filemar, mar)
        save(file = fileemar, emar)
        save(file = filemarz, marz)
        save(file = fileemarz, emarz)
        message("...done")
    }
    ## extinction##############################################################################
    if (REDO_EXTINCTION) {
        message("in silico extinction of distribution cells ...")
        load(filegenemaps)
        load(filemar)
        load(fileemar)
        ### Simulate extinctions
        # random losss of cells in grid
        random.X <-
            MARextinction_sim(genemaps,
                scheme = "random",
                xfrac = 0.01
            ) %>%
            mutate(type = "random")
        # extinction from outside to
        inward.X <-
            MARextinction_sim(genemaps,
                scheme = "inwards",
                xfrac = 0.01
            ) %>%
            mutate(type = "inwards")
        # from outside to the midpoint value
        inward.X.center <-
            MARextinction_sim(genemaps,
                scheme = "inwards",
                centerfun = function(x) (min(x) + max(x)) / 2
            ) %>%
            mutate(type = "inwards.center")
        # from south to north
        sn.X <-
            MARextinction_sim(genemaps,
                scheme = "southnorth",
                xfrac = 0.01
            ) %>%
            mutate(type = "southnorth")
        # radial from a central point
        radial.X <-
            MARextinction_sim(genemaps,
                scheme = "radial",
                centerfun = median,
                xfrac = 0.05
            ) %>%
            mutate(type = "radial")
        radial.spain.X <-
            MARextinction_sim(genemaps,
                scheme = "radial",
                centerfun = min,
                xfrac = 0.05
            ) %>%
            mutate(type = "radial.xymin")
        radial.scand.X <-
            MARextinction_sim(genemaps,
                scheme = "radial",
                centerfun = max,
                xfrac = 0.05
            ) %>%
            mutate(type = "radial.xymax")
        # assemble all simulations
        xsim <- rbind(
            random.X,
            inward.X,
            inward.X.center,
            sn.X,
            radial.X,
            radial.spain.X,
            radial.scand.X
        )
        save(file = fileextinctionsim, xsim)

        # Create a mar of the extinction values themselves
        marsim <- sars::sar_power(data = dplyr::select(xsim, asub, M))
        save(file = fileextinctionmar, marsim)

        # quantify how good the fit is
        # R2_randomMAR<-marR2(mar, dplyr::select(xsim,asub,M))
        # R2_randomEmAR<-marR2(emar, dplyr::select(xsim,asub,M))
        # R2sim<-list(R2_randomMAR,R2_randomEmAR)

        # save objects
        # save(file = fileextinctionsimR2, R2sim)

        message("...done")
    }
    ## end##############################################################################
    message("END")
}
