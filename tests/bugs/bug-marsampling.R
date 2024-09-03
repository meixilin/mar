# this is a demonstration on how old MARsampling method cannot sample the full extent of the maps
rm(list = ls())
load('../testdata/mares_new-joshua.rda')
load('../testdata/gm-joshua.rda')

get_extents <- function(gm, string) {
    bbox = as.integer(strsplit(string, ';')[[1]])
    out = rowcol_extent(gm, bbox[1], bbox[2], bbox[3], bbox[4])
    return(out)
}

plotter <- function(ii) {
    plot(gm$samplemap)
    for (jj in 1:10) {
        id = (ii-1)*10+jj
        print(id)
        plot(mares_ext[[id]], add = T, col = sample(grDevices::colors(),1))
    }
}

mares_ext = lapply(mares$extent, get_extents, gm = gm)

# ii = 1, all boxes contains 4 cells instead of 1 cell only
png(filename = 'bug_marsampling_4cells.png')
plotter(1)
dev.off()

png(filename = 'bug_marsampling_sidesize28.png')
plotter(28)
dev.off()

png(filename = 'bug_marsampling_sidesize27.png')
plotter(27)
dev.off()

