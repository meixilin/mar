# test MARPIPELINE
library(mar)

name = "joshua"
workdir = "~/AMoiLab/mar_related/mar/tests/testthat/MARPIPELINE"
genofile = "/Users/linmeixi/AMoiLab/mar_related/mar/tests/testthat/testdata/fileparser/genome.joshua.tsv"
lonlatfile = "/Users/linmeixi/AMoiLab/mar_related/mar/tests/testthat/testdata/fileparser/lonlat.joshua.csv"
MARPIPELINE(name = name, workdir = workdir, genofile = genofile, lonlatfile = lonlatfile, saveobj = TRUE)

obj <- MARsadsfs(AC, N = 290, ploidy = 2)

AC <- mar:::.get_AC(gm$geno)
N = 290
ploidy = 2

library(sads)
mar:::plot.sfs(mar::sfs(moths, N = 500, ploidy = 2))
moths.mzsm <- fitsad(moths, "mzsm")
vars = coef(moths.mzsm)

moths.mzsm.r <- radpred(moths.mzsm)
test = moths.mzsm.r
all.equal(moths.mzsm.r, test)

yy = rad(rmzsm(n = 240, J = vars[1], theta = vars[2]))


# load genomaps
setwd(workdir)
load("gm_joshua2.rda")
load("marsadsfs_joshua3.rda")
load("mardflist_joshua3.rda")
load("marlist_joshua3.rda")
load("extdflist_joshua3.rda")

methods(class = "marmaps")
plot(gm$maps)

option_mar = list(scheme = "inwards", nrep = 10, quorum = FALSE, animate = FALSE, myseed = NULL)
scheme = "inwards"
mardf = MARsampling(gm = gm, scheme = scheme, nrep = option_mar$nrep, xfrac = option_mar$xfrac,
            quorum = option_mar$quorum, animate = option_mar$animate, myseed = option_mar$myseed)

mar:::.MARsampling_schemes

extdf = mar:::MARextinction(gm,scheme = "random")

load("marlist_joshua2.rda")
load("extdflist_joshua3.rda")
load("extlist_joshua2.rda")

extdf = extdflist[[1]]
ext = extlist[[1]]


mardf = mardflist[[1]]
df = marlist[[1]]
str(marlist[[1]])
