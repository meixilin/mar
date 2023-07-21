

extractheritability <- function(h2fit) {
    # mean(h2fit)
    # plot(h2fit)
    # autocorr(h2fit)
    # formateo<-(c(mean(h2fit),HPDinterval(h2fit)[1],HPDinterval(h2fit)[2]))
    # plot(proph2[,1],type="l")
    # plot(density(proph2[,1]))
    plot(h2fit, type = "l")

    formateo <- c(mean(h2fit), quantile(h2fit, p = 0.025), quantile(h2fit, p = 0.975))
    formateo <- round(formateo, digits = 4)
    toexport <- (paste(formateo[1], " (", formateo[2], " , ", formateo[3], ")", sep = ""))
    return(toexport)
}

exh2 <- function(modelrandom = model$VCV) {
    allrandom <- apply(modelrandom, 1, FUN = sum)
    proph2 <- apply(modelrandom, 2, FUN = function(x) {
        x / allrandom
    })
    apply(proph2, 2, FUN = extractheritability)
}

summaryglmmMCMC <- function(model) {
    fixres <- summary(model)$solutions
    term <- rep("fixed", dim(fixres)[1])
    fixres <- cbind(fixres, term)
    randres <- summary(model)$Gcovariances
    errorres <- summary(model)$Rcovariances

    propVar <- randres[, "post.mean"] / sum(randres[, "post.mean"] + errorres[, "post.mean"])
    propVarerror <- errorres[, "post.mean"] / sum(randres[, "post.mean"] + errorres[, "post.mean"])

    term <- rep("random", length(propVar))
    randres <- cbind(randres, propVar, term)
    term <- "residual"
    errorres <- cbind(errorres, propVarerror, term)

    toreturn <- rbind(fixres, randres, errorres)
    return(toreturn)
}


myglmmMCMC.f <- function(data, formfix, formran, uncorrelated = F) {
    library(MCMCglmm)

    if (uncorrelated == F) {
        lmfull <- MCMCglmm(
            data = data,
            formfix, random = formran,
            pl = T, pr = T,
            nitt = 50000, burnin = 5000,
            ginverse = list(Genotype = Ai)
        )
    }
    if (uncorrelated == T) {
        lmfull <- MCMCglmm(
            data = data,
            formfix, random = formran,
            pl = T, pr = T,
            nitt = 50000, burnin = 5000
        )
    }

    tr <- summaryglmmMCMC(lmfull)
    print(tr)
    print("Deviance Information Criterion:")
    print(lmfull$DIC)
    return(list(summary = tr, DIC = lmfull$DIC))
}
