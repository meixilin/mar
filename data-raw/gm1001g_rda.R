# Generate RDA file for mar package
options(echo = TRUE)
library(mar)

genodata <- text_parser(
    geno.fn = '1001g_genotypes.txt.gz',
    samp.fn = '1001g_accessions.txt',
    pos.fn = '1001g_chrpos.txt',
    ploidy = 2
)

mapsdata <- lonlat_parser('1001g_lonlat.txt', mapcrs = 'EPSG:8857')

gm1001g <- genomaps(genodata, mapsdata)

usethis::use_data(gm1001g)
