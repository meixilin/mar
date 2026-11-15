# Generate RDA file for mar package
options(echo = TRUE)
library(mar)

genodata <- vcf_parser('gmexp.vcf.gz')
# random locations on an arbitrary plane: "local" gives terra a Cartesian
# (meter) CRS so cellSize() area still works, without a real-world projection
mapsdata <- lonlat_parser('gmexp_lonlat.csv', incrs = "local", mapcrs = "local")

gmexp <- genomaps(genodata, mapsdata)

usethis::use_data(gmexp)
