# test MARPIPELINE
library(mar)

name = "joshua"
workdir = "~/AMoiLab/mar_related/mar/tests/testthat/MARPIPELINE"
genofile = "/Users/linmeixi/AMoiLab/mar_related/mar/tests/testthat/testdata/fileparser/genome.joshua.tsv"
lonlatfile = "/Users/linmeixi/AMoiLab/mar_related/mar/tests/testthat/testdata/fileparser/lonlat.joshua.csv"
MARPIPELINE(name = name, workdir = workdir, genofile = genofile, lonlatfile = lonlatfile, saveobj = TRUE)


option_mar = list(scheme = "inwards", nrep = 10, quorum = FALSE, animate = FALSE, myseed = NULL)
scheme = "inwards"
mardf = MARsampling(gm = gm, scheme = scheme, nrep = option_mar$nrep, xfrac = option_mar$xfrac,
            quorum = option_mar$quorum, animate = option_mar$animate, myseed = option_mar$myseed)

mar:::.MARsampling_schemes

extdf = mar:::MARextinction(gm,scheme = "random")
