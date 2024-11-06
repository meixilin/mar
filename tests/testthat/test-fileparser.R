# Mock data creation functions
gen_mock <- function(myext, content) {
    filename <- tempfile(fileext = myext)
    writeLines(content, filename)
    return(filename)
}

unlink_files <- function() {
    unlink(list.files(path = tempdir(), pattern = "^file", full.names = TRUE))
}

vcf_content <- c(
    "##fileformat=VCFv4.2",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3",
    "1\t100\trs1\tA\tT\t100\tPASS\tAF=0.5\tGT:DP\t0/0:30\t0/1:35\t1/1:32",
    "1\t200\trs2\tG\tC\t95\tPASS\tAF=0.16\tGT:DP\t0/1:28\t0/0:31\t0/0:29",
    "2\t300\trs3\tT\tA\t98\tPASS\tAF=0.16\tGT:DP\t0/0:33\t0/0:30\t0/1:34",
    "2\t400\trs4\tC\tG\t99\tPASS\tAF=0.66\tGT:DP\t0/1:31\t1/1:36\t0/1:32"
)

gds_checks <- function(result, genomat, alleles) {
    expect_s4_class(result, "SeqVarGDSClass")

    # Check the contents of the GDS file
    expect_equal(.get_genodata(result, "sample.id"), c("Sample1", "Sample2", "Sample3"))
    expect_equal(.get_genodata(result, "variant.id"), 1:4) #ID not kept here in GDS
    expect_equal(.get_genodata(result, "chromosome"), c("1", "1", "2", "2"))
    expect_equal(.get_genodata(result, "position"), c(100, 200, 300, 400))
    expect_equal(.get_genodata(result, "genotype"), genomat)
    expect_equal(seqGetData(result, "allele"), alleles)
    seqClose(result)
    return(invisible())
}


# Test cases
test_that(".strip_ext works correctly", {
    expect_equal(.strip_ext("file.txt", c(".txt")), "file")
    expect_equal(.strip_ext("file.txt.gz", c(".txt", ".txt.gz")), "file")
    expect_error(.strip_ext("file.txt", c(".csv", ".gz")))
    expect_error(.strip_ext("file.txt.gz", c(".txt.gz", ".gz")))
})

test_that(".guess_delim works correctly", {
    expect_equal(.guess_delim("A,B,C"), ",")
    expect_equal(.guess_delim("A\tB\tC"), "\t")
    expect_error(.guess_delim("A B C"))
})

test_that(".read_genotype works correctly", {
    # working case
    res <- .read_genotype(gen_mock(".txt", "0\t1\t2\n1\t0\t2\n"), het2hom = FALSE)
    expect_equal(res, matrix(c(0,1,2,1,0,2), nrow=2, byrow=TRUE))
    # het2hom
    res <- .read_genotype(gen_mock(".txt", "0\t1\t2\n1\t0\t2\n"), het2hom = TRUE)
    expect_equal(res, matrix(c(0,1,1,1,0,1), nrow=2, byrow=TRUE))
    # error case
    expect_error(.read_genotype(gen_mock(".txt", "0\t1\t2\n1\t0\t2\n3\t1\t0\t2\n")))
    expect_error(.read_genotype(gen_mock(".txt", "0\t1\t2\n1\t0\tNA\n"), het2hom = FALSE))
    unlink_files()
})

test_that(".read_pos works correctly", {
    # working case
    res <- .read_pos(gen_mock(".txt", "CHROM\tPOS\n1\t100\n2\t200"))
    expect_equal(res[[1]], c(1, 2))
    expect_equal(res[[2]], c(100, 200))
    unlink_files()
})

test_that(".read_samp works correctly", {
    res <- .read_samp(gen_mock(".txt", "sample1\nsample2\nsample3"))
    expect_equal(res, c("sample1", "sample2", "sample3"))
    expect_error(.read_samp(gen_mock(".txt", "sample1,data1\nsample2,data2\nsample3,data3")))
    unlink_files()
})

test_that(".read_lonlat works correctly", {
    # working case with comma delimiter
    res <- .read_lonlat(gen_mock(".txt", "SAMPID,LONGITUDE,LATITUDE\n1,123.45,67.89\n2,-45.67,12.34\n"))
    expect_equal(res$LONGITUDE, c(123.45, -45.67))
    expect_equal(res$LATITUDE, c(67.89, 12.34))

    # working case with tab delimiter and shorter column names
    res <- .read_lonlat(gen_mock(".txt", "ID\tLON\tLAT\n1\t123.45\t67.89\n2\t-45.67\t12.34\n"))
    expect_equal(res$LON, c(123.45, -45.67))
    expect_equal(res$LAT, c(67.89, 12.34))

    # error cases
    expect_error(.read_lonlat(gen_mock(".txt", "ID,LAT,LON\n1,67.89,123.45\n"))) # wrong column order
    expect_error(.read_lonlat(gen_mock(".txt", "ID,LONGITUDE\n1,123.45\n"))) # missing latitude
    unlink_files()
})

test_that("text_parser works correctly", {
    temp_geno <- gen_mock(".txt", "0\t1\t2\n1\t0\t2\n")
    temp_pos <- gen_mock(".txt", "CHROM\tPOS\n1\t100\n2\t200")
    result <- text_parser(temp_geno, pos.fn = temp_pos, het2hom = FALSE)
    expect_s3_class(result, "margeno")
    expect_equal(result$sample.id, 1:3)
    expect_equal(result$variant.id, c(1, 2))
    expect_equal(result$position, c(100, 200))
    expect_equal(result$chromosome, c("1", "2"))
    expect_equal(result$genotype, matrix(c(0,1,2,1,0,2), nrow=2, byrow=TRUE))
    expect_equal(result$ploidy, 2)
    result <- text_parser(temp_geno, pos.fn = temp_pos)
    expect_equal(result$genotype, matrix(c(0,1,1,1,0,1), nrow=2, byrow=TRUE))

    # not input pos.fn
    result <- text_parser(temp_geno)
    expect_equal(result$position, NULL)
    expect_equal(result$chromosome, NULL)
    unlink_files()
})

test_that("vcf_parser and .get_genodata works correctly", {
    temp_vcf <- gen_mock(".vcf", vcf_content)
    result <- vcf_parser(temp_vcf, opengds = TRUE)
    genomat <- matrix(c(0,1,2,1,0,0,0,0,1,1,2,1), nrow=4, byrow=TRUE)
    gds_checks(result, genomat, c("A,T","G,C","T,A","C,G"))
    unlink(c(temp_vcf, result$filename))
    unlink_files()
})

test_that("plink_parser works correctly", {
    # requires PLINK installed
    temp_vcf <- gen_mock(".vcf", vcf_content)
    temp_prefix <- .strip_ext(temp_vcf, ".vcf")
    plink="/Applications/Lab/plink_mac_20221210/plink" # TODO: cleanup
    system(paste0(plink, " --vcf ", temp_vcf, " --make-bed --out ", temp_prefix))
    result <- plink_parser(temp_prefix, opengds = TRUE)
    # plink automatically flips the alleles by MAF
    genomat <- matrix(c(0,1,2,1,0,0,0,0,1,1,0,1), nrow=4, byrow=TRUE)
    gds_checks(result, genomat, c("A,T","G,C","T,A","G,C"))
    unlink(list.files(pattern = temp_prefix))
    unlink_files()
})

test_that("lonlat_parser works correctly", {
    # working case
    result <- lonlat_parser(gen_mock(".txt", "ID\tLONGITUDE\tLATITUDE\nSample1\t-73.935242\t40.730610\nSample2\t-118.243683\t34.052235"))
    lonlat <- matrix(c(-73.935242, 40.730610, -118.243683, 34.052235), nrow=2, byrow=TRUE)
    colnames(lonlat) <- c("LONGITUDE", "LATITUDE")
    # error case
    expect_equal(result[[1]], c("Sample1", "Sample2"))
    expect_equal(result[[2]], lonlat)
    unlink_files()

})

# remove any temp files
unlink_files()
closeAllConnections()
