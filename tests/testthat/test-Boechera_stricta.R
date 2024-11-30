rm(list = ls())
library(mar)

genofile="/global/scratch/projects/fc_moilab/kmua/MAR2/Boechera_stricta_full_concatened.filtered.final.1.vcf.gz"
lonlatfile = "/global/scratch/projects/fc_moilab/kmua/MAR2/05_Bstricta_MitchellOlds_wgs/latlong.tsv"
workdir = "/global/scratch/projects/fc_moilab/meixilin/mar/tests/testthat/testdata/bost"

# run MARPIPELINE
devtools::load_all()

MARPIPELINE(name = "bost",
            genofile = genofile,
            lonlatfile = lonlatfile,
            workdir = workdir,
            filetype = 'vcf',
            marsteps = c('ext'),
            saveobj = TRUE)

# marextinct rerun --------
setwd(workdir)
load('gm_bost.rda')
devtools::load_all()

# same as MARsampling ------------------------------------------------------
# match schemes (default to random)
scheme = "random"
# calculate and store raster area in the given gm$maps$samplemap
gmarea = .areaofraster(gm$maps$samplemap)
# the point where most samples are available (for inwards / outwards sampling)
maxids <- raster::which.max(gm$maps$samplemap)
if (length(maxids) > 1) {
    warning('More than one cell with maximum samples')
}
r0c0 <- raster::rowColFromCell(gm$maps$samplemap, maxids[1])
# End same as MARsampling --------------------------------------------------
xfrac = 0.01
nrep = 10
extlist <- .extlist_sample(gm, xfrac, scheme, nrep, r0c0)

# if need to plot
animate = TRUE
if (animate) {
    lapply(extlist, .animate_MARextinction, gm = gm)
}

# calculate area and genetic diversity in each extant cell list
ii = 1
outlist <- lapply(seq_along(extlist), function(ii) {
    outl <- lapply(extlist[[ii]], mutdiv.cells, gm = gm, gmarea = gmarea)
    out <- do.call(rbind, lapply(outl, as.data.frame, stringsAsFactors = FALSE))
    # append end theta (zero in all)
    out[nrow(out)+1, ] <- rep(0, ncol(out))
    out$repid <- ii # replicate id
    return(out)
})
outdf <- do.call(rbind, outlist)

# set outdf as a marsamp class
class(outdf) <- c(class(outdf), "marextinct") # marextinction output class
attr(outdf, 'scheme') <- scheme