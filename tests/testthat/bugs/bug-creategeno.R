# this is able to reproduce the genemaps from coords and genome data
# therefore, coords and genome can be directly used to generate mutdiv with new sampling methods
# BUG: THE OLD create_gene_map can miss some samples at the edge. Need updates.

load('tests/testthat/testdata/coords-joshua.rda')
load('tests/testthat/testdata/genome-joshua.rda')
load('tests/testthat/testdata/genemaps-joshua.rda')

library(mar)
genemaps_rep = create_gene_maps(coord = coords, geno = genomes, bystep = 0.005)
# all TRUE
compareRaster(genemaps[[1]], genemaps_rep[[1]])
compareRaster(genemaps[[2]], genemaps_rep[[2]])

# BUG note the shrinked extent from the generated raster compared with lonlat dataframe
apply(coords, 2, range)
extent(genemaps[[1]])

# you can see it from this plot:
absmap = genemaps_rep[[1]]
absmap[absmap > 0] = 1
plot(absmap)
abline(h = min(coords[,2]))
abline(h = extent(genemaps[[1]])@ymin)
points(coords[,1], coords[,2])

# BUG: 3 samples were left out
cellStats(genemaps_rep[[1]],'sum')
nrow(coords)



