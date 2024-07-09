

test_that("genemaps works", {
    library(raster)
    myspecies = c('joshua', 'arabidopsis')
    for (ii in myspecies) {
        load(paste0('testdata/coords-', ii, '.rda'))
        load(paste0('testdata/genome-', ii, '.rda'))
        # test for previous one
        gm <- genemaps(coords, genomes)
        # write new gm if it doesn't exists
        if (!file.exists(paste0('testdata/gm-', ii, '.rda'))) {
            save(gm, file = paste0('testdata/gm-', ii, '.rda'))
        }

        # test that the cellid is correctly relinked
        cellids = gm$sample$cellid
        samplemap = gm$samplemap
        rmat = t(as.matrix(samplemap))
        rmatid = which(!is.na(rmat))
        expect_true(all(rmat[rmatid] == table(cellids)))
        ## plot for the first sample
        plot(samplemap)
        plot(rasterFromCells(samplemap, cells = cellids[1]), col = 'red', add = T)
        points(x = coords[1,1], y = coords[1,2])
        # test that the number of samples is correct
        expect_equal(nrow(coords), cellStats(gm$samplemap, 'sum'))
    }
})
