# Generate RDA file for mar package
options(echo = TRUE)
library(mar)

genodata <- mar::text_parser(
    geno.fn = '1001g_genotypes.txt.gz', 
    samp.fn = '1001g_accessions.txt', 
    pos.fn = '1001g_chrpos.txt', 
    ploidy = 2
)

mapsdata <- mar::lonlat_parser('1001g_lonlat.txt')
mapsdata <- mar::marmaps(mapsdata, mapres = NULL, mapcrs = "+proj=longlat +datum=WGS84")

gm1001g <- mar::genomaps(genodata, mapsdata)

usethis::use_data(gm1001g)
