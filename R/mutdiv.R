
# genetic diversity estimator
mutdiv <- function(raster_samples, raster_mutmaps, rest_mutmaps) {
    require(raster)
    # Get the number of samples
    N <- sum(fn(values(raster_samples)), na.rm = T)
    # freqs
    P <- fn(apply(values(raster_mutmaps), 2, function(cells) sum(cells, na.rm = T))) / N
    # at least one cell has to have a presence of a mutation, and sum over
    M_ <- fn(apply(values(raster_mutmaps), 2, function(cells) any(cells > 0)))
    M_[is.na(M_)] <- 0
    M <- sum(M_)
    # find endemisms
    E_ <- apply(values(rest_mutmaps), 2, function(cells) any(cells > 0))
    E_[is.na(E_)] <- 0
    table(M_, E_)
    E <- sum(M_ & !E_)
    # Get the number of SNPs for the sample
    L <- dim(raster_mutmaps)[3]
    # compute diversity, Theta Waterson and Theta Pi (pairwise)
    if (N > 1 & M > 0) {
        theta <- M / (Hn(N-1) * L)
        thetapi <- (N/(N-1)) * sum(2 * P * (1 - P), na.rm = T) / L
    } else {
        theta <- 0
        thetapi <- 0
    }
    # area taking into account only grid cells with data
    # asub= sum(raster_samples[] > 0, na.rm = T) * (res(raster_samples)[1]*res(raster_samples)[2])
    asub <- areaofraster(raster_samples)

    # area based on simple square
    a <- dim(raster_samples)[1] * res(raster_samples)[1] * dim(raster_samples)[2] * res(raster_samples)[2]
    # return
    return(data.frame(thetaw = theta, pi = thetapi, M = M, E = E, N = N, a = a, asub = asub))
}
