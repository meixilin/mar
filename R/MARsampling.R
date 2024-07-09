allowed_schemes = c("random", "inwards", "outwards", "southnorth", "northsouth")

# match with old MARsampling
# TODO: Add not gridded sampling
MARsampling <- function(gm, scheme = allowed_schemes, nrep = 10, gridded = TRUE, plot = FALSE, myseed = NULL) {
    # set seed if specified
    if (!is.null(myseed)) {
        set.seed(myseed)
    }
    # calculate and store raster area in the given gm$samplemap
    gmarea = areaofraster(gm$samplemap)
    # the x and y number of cells in gm$samplemap
    lonrange <- dim(gm$samplemap)[1]
    latrange <- dim(gm$samplemap)[2]
    minrange <- min(lonrange, latrange) # the maximum size box can become
    # match schemes (default to random)
    scheme = match.arg(scheme)
    # differences btw different schemes are in the bounding boxes selected for diversity calculations
    bboxlist = switch(scheme,
                      random = {
                          lapply(1:minrange, function(sidesize) {
                              ll = list()
                              # TODO: did not update for reproducibility. Will update next.
                              for (ii in 1:nrep) {
                                  xmin <- sample(1:(lonrange - sidesize), 1); xmax <- xmin + sidesize
                                  ymin <- sample(1:(latrange - sidesize), 1); ymax <- ymin + sidesize
                                  ll = c(ll, list(c(xmin, xmax, ymin, ymax)))
                              }
                              return(ll)
                          }) %>% unlist(., recursive = FALSE)
                      },
                      inwards = {

                      }
    )
    # calculate area and genetic diversity in each bounding boxes
    outlist = lapply(bboxlist, .mutdiv.gridded, gm = gm, gmarea = gmarea)
    outdf = do.call(dplyr::bind_rows, lapply(outlist, as.data.frame))
    return(outdf)
}
#
# # sampling methods
# .random_sample <- function(gm) {
#     # ranges of samples
#     lonrange <- dim(genemaps[[1]])[1]
#     latrange <- dim(genemaps[[1]])[2]
# }
