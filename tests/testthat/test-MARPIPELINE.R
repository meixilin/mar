# test MARPIPELINE
name = "joshua"
workdir = "~/AMoiLab/mar_related/mar/tests/testthat/MARPIPELINE"
genofile = "/Users/linmeixi/AMoiLab/mar_related/mar/tests/testthat/testdata/fileparser/genome.joshua.tsv"
genofile = "/Users/linmeixi/AMoiLab/mar_related/mar/tests/testthat/testdata/fileparser/genome.joshua.csv" # bad file

lonlatfile = "/Users/linmeixi/AMoiLab/mar_related/mar/tests/testthat/testdata/fileparser/lonlat.joshua.csv"

filetype = 'text'
marsteps = "gm"

extra_file = list(samplefile = NULL, posfile = NULL, subsample = NULL, subvariant = NULL)
option_geno = list(ploidy = 2, het2hom = FALSE, maxsnps = 1000000, nreps = 10)
option_map = list(mapres = NULL, mapcrs = "+proj=longlat +datum=WGS84")


text_parser(genofile, samp.fn = extra_file$samplefile, pos.fn = extra_file$posfile,
                               ploidy = option_geno$ploidy, het2hom = option_geno$het2hom)