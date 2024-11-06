
# TODO: Add support for filtering by genomic regions and samples
sfs <- function(gm, folded = TRUE, nonzero = TRUE) {
    table(rowSums(gm$genotype))
}


# expected SFS
predSFS <- function(M, n = 1000) {
    n = 20
    M = 1:1000
    n = 1000
    lenM = 1593510
    # abundance = c 1/q
    # log(snps) = c
    # q=M/max(M,na.rm=TRUE) * 0.5 # because they are oriented to minor
    # q= round(q *100)/ 100
    # countsq=table(q)
    # lm(log(countsq) ~ log(fn(names(countsq))))
    theta = lenM / Hn(n)
    x = 1:n
    expsfs = theta/x
    fsfs = (expsfs + rev(expsfs))[1:(n/2)]
    barplot(fsfs)
    hist(m.sfs, breaks = 500)

    samp <- sample(
        x = seq(1, 49, 0.01) / 100, size = length(M),
        prob = 1 / (seq(1, 49, 0.01) / 100), replace = TRUE
    )
    # samp=samp * max(M,na.rm=TRUE) #* 2 # because I used 0.49 as maximum frequency
    samp <- ceiling(samp * n) # to make this a
    hist(samp)
    return(fn(samp))

}


# Open the GDS file
gds_file <- "path/to/your/file.gds"
gds <- seqOpen(gds_file)

# Get the sample size
sample_size <- seqSummary(gds)$sample.id

# Calculate allele frequencies
allele_freq <- seqAlleleFreq(gds)

# Close the GDS file
seqClose(gds)

# Function to calculate SFS
calculate_sfs <- function(freq, n) {
    sfs <- table(round(freq * (2 * n)))
    sfs <- sfs[-c(1, length(sfs))]  # Remove 0 and 2n
    return(sfs)
}

# Generate SFS
sfs <- calculate_sfs(allele_freq, length(sample_size))

# Plot SFS
barplot(sfs,
        main="Site Frequency Spectrum",
        xlab="Minor Allele Count",
        ylab="Frequency")

# Print SFS data
print(sfs)

# Optionally, save SFS data to a file
write.table(sfs, file="sfs_data.txt", sep="\t", quote=FALSE)


message("site frequency spectrum ...")
freqfile = './testdata/1001g/1001gbi.chr1.frq'
freq <- data.table::fread(freqfile) # %>% head(nsnps) # uncomment for all data
M <- (ceiling(freq$MAF * freq$NCHROBS))
M[is.na(M)] <- 0
M <- rep(1, times = 100)
expsfs0 <- mar:::predSFS(M, n = 10)
hist(expsfs, breaks = 1:5)
M <- M[M != 0]

lenM = 100
theta = lenM / mar:::Hn(n = 10)
n = 10
x = 1:n
expsfs = theta/x
fsfs = (expsfs + rev(expsfs))[1:(n/2)]

hist(M)
Mrad <- rad(M)
# fit models
# m.g <- fitsad(M, "gamma",)
m.w <- fitsad(M, "weibull")
m.ln <- fitsad(M, "lnorm")
m.ls <- fitsad(M, "ls")
m.bs <- fitsad(M, "bs")
m.geom <- fitsad(M, "geom")
m.mzsm <- fitsad(M, sad = "mzsm")
# m.poi <- fitsad(M, "poilog")
m.sfs <- mar:::predSFS(M)
hist(m.sfs)

tableaic <- AICtab(
    m.bs,
    m.mzsm,
    # m.poi,
    m.ln,
    m.ls,
    # m.g,
    m.w,
    m.geom
)
# tableR2<-list(
#   "bs"=pseudoR2(Mrad$abund,radpred(m.bs)$abund),
#   "mzsm"=pseudoR2(Mrad$abund,radpred(m.mzsm)$abund),
#   "lognormal"=pseudoR2(Mrad$abund,radpred(m.ln)$abund),
#   "logseries"=pseudoR2(Mrad$abund,radpred(m.ls)$abund),
#   "gamma"=pseudoR2(Mrad$abund,radpred(m.g)$abund),
#   "weibull"=pseudoR2(Mrad$abund,radpred(m.w)$abund),
#   "geom"=pseudoR2(Mrad$abund,radpred(m.geom)$abund),
#   "sfs"=pseudoR2(Mrad$abund,rad(m.sfs)$abund)
# )
liksfs <- likelihoodSFSfold(M)
if (saveobjects) save(file = filetable, object = tableaic)
if (saveobjects) save(file = fileliksfs, object = liksfs)
