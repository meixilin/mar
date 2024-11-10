
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
