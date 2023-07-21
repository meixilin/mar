require(dplyr)
################################################################################
whittakerplot2 <- function(y, color = "grey") {
    myM <- sads::rad(y)
    ggplot(data.frame(myM)) +
        geom_point(aes(x = rank, y = abund), color = color) +
        geom_line(aes(x = rank, y = abund), color = color) +
        # geom_point(aes(x=(rank(-y)), y=log(y*100)+1e-10) ,color=color)+
        # geom_line(aes(x=(rank(-y)), y=log(y*100)+1e-10) ,color=+color)+
        labs(y = "log (% abundance)", x = "rank of abundance") +
        scale_y_log10()
}
whittakerplot <- whittakerplot_base <- function(y, color = "grey") {
    ggplot(data.frame(y)) +
        geom_point(aes(x = (rank(-y)), y = log(y * 100) + 1e-10), color = color) +
        geom_line(aes(x = (rank(-y)), y = log(y * 100) + 1e-10), color = color) +
        labs(y = "log (% relative abundance)", x = "rank of abundance")
}
whittakerplot_add <- function(y, color = "grey", p) {
    p +
        geom_point(data = data.frame(y), aes(x = (rank(-y)), y = log(y * 100) + 1e-10), color = color) +
        geom_line(data = data.frame(y), aes(x = (rank(-y)), y = log(y * 100) + 1e-10), color = color)
}
prestonplot <- prestonplot_base <- function(y, color = "grey") {
    ggplot(data.frame(y)) +
        geom_histogram(aes(x = log(y)), fill = color, color = "white") +
        labs(y = "# of species", x = "bin of abundance (log # counts)")
}
prestonplot_add <- function(y, color = "grey", p) {
    p + geom_histogram(data.frame(y), aes(x = log(y)), fill = color, color = "white")
}
################################################################################

# deterministic geometric series
geometricseries50 <- function(S = 100, N = 10000) {
    # deterministic series taking fixed 50% fraction
    out <- sapply(1:S, function(i) (1 / (2))^i)
    if (!is.null(N)) out <- ceiling(out * N)
    return(out)
}
geometricseriesgeneral <- function(S = 100, k = 0.5, N = 10000) {
    # deterministic series taking fixed k fraction
    # every time a new species arrives
    out <- sapply(1:S, function(i) k * (1 - k)^(i - 1))
    if (!is.null(N)) out <- ceiling(out * N)
    return(out)
}
dominancepreemption <- function(S = 100, N = 10000) {
    # stochastic series taking fractions  0.5-1
    # of niches as species arrive.
    # This should converges to geometric series with k=0.75
    y <- c()
    y[1] <- 1
    for (i in 2:(S + 1)) {
        y[i] <- y[i - 1] * runif(1, min = 0.5, max = 1)
    }
    out <- y[-1]
    if (!is.null(N)) out <- ceiling(out * N)
    return(out)
}
# randomfraction<-function(S, N=NULL){
#   # stochastic series similar to dominance pre-emption
#   # but taking habitat fractions from 0 to 1.
#  y<-c()
#  y[1]=1
#  for(i in 2:(S+1)){
#    y[i] <-  y[i-1] * runif(1,min=0,max=1) <- need to reimplement because any fragment should have same prob being broken
#  }
# out= y[-1]
# if(!is.null(N)) out = ceiling(out * N)
# return(out)
# }
# composite<-function(S){stop("composite model not implemented!")}
# randomassortment<-function(S){stop("random assortment not implemented!")}
# dominancedecay<-function(S){stop("dominancedecay not implemented!")}

brokenstick <- function(S = 100, N = 10000) {
    # Also called MacArthur model 1957 and 1960
    # which is similar to random fraction except
    # for that larger habitat fragments are more
    # likely to be colonized (i.e. brroken by next piece)
    out <- (1 / S) * sapply(1:S, function(i) {
        sum(sapply(i:S, function(x) 1 / x))
    })
    if (!is.null(N)) out <- ceiling(out * N)
    return(out)
}
.normalize <- function(x) {
    (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
}
logseries <- function(S = 100, N = 1000,
                      alpha = 0.8,
                      # x=N/(alpha+N)
                      x = 1) {
    # Fiser 1943 created a new distribution derived from
    # that observations of any given species are Poisson distributed.
    # This function of species of given abundance is the
    # probability density function of the log series (-1/ log(1-k)) * (k^i/i)
    # the x axis is not the rank of species, but the units of abundance
    possibleabundances <- 1:round(N * 0.5) # the bound is that all N indiuviduals are 1 species
    # athough realistically no species is going to have more than 50% of the data
    out <- sapply(possibleabundances, function(i) alpha * x^i / i)
    out <- sample(x = 1:length(out), size = S, prob = .normalize(out) + 1) # no need to do proportaional to obs, as this model already incorporates that
    return(rev(sort(out)))
}
lognormal <- function(S, N = 10000,
                      mu = N / S,
                      sigma = 1000) {
    # S_mode= N/S,
    # R=seq(-S*0.1,S*0.1)
    # Preston 1948.
    out <- dnorm(x = seq(0, N), mean = mu, sd = sigma, log = T)
    out <- sample(log(seq(0, N)), size = S, prob = .normalize(out) + 1, replace = T)
    return(rev(sort(out)))
}
################################################################################
whittakerplot_combine <- function(listabundances = list(y_geometric, y_dominancepreemption, y_broken, y_logseries, y_lognormal),
                                  listlabels = list("geometric", "dominance preemption", "broken stick", "log series", "log normal")) {
    toplot <- do.call(
        rbind,
        lapply(1:length(listabundances), function(i) {
            ytmp <- data.frame(y = listabundances[[i]])
            ytmp$y <- ytmp$y / max(ytmp$y, na.rm = T)
            ytmp$x <- 1:length(listabundances[[i]])
            ytmp$method <- listlabels[[i]]
            return(ytmp)
        })
    )
    toplot$method <- factor(toplot$method, levels = as.character(listlabels))
    ggplot(toplot) +
        geom_point(aes(y = log(y), x = rank(x), color = method, group = method)) +
        geom_line(aes(y = log(y), x = rank(x), color = method, group = method)) +
        scale_color_manual("method", values = brewer.pal(length(listlabels), "Reds"), breaks = as.character(listlabels)) +
        labs(y = "log (% relative abundance)", x = "rank of abundance")
}
