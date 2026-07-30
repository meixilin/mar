# Generate RDA file for mar package
options(echo = TRUE)
library(mar)

genodata <- vcf_parser('gmexp.vcf.gz')
# no crs needed for random locations
mapsdata <- lonlat_parser('gmexp_lonlat.csv', incrs = "", mapcrs = "")

gmexp <- genomaps(genodata, mapsdata)

usethis::use_data(gmexp)
