require(dplyr)
################################################################################
martable <- function(marmod) {
    if (!is.null(marmod$Model)) {
        mytable <- marmod$Model$coefficients
        colnames(mytable)[4] <- "P"
        rownames(mytable)[1] <- "logC"
        rownames(mytable)[2] <- "z"
        mytable[, 4] <- format(mytable[, 4], scientific = TRUE, digits = 3)
    } else {
        mytable <- marmod$sigConf
        colnames(mytable)[4] <- "P"
    }
    return(mytable)
}

marcoef <- function(marmod, dig = 3, conf = 0.95) {
    require(dplyr)
    if (!is.null(marmod$Model)) {
        qs <- (1 - conf) / 2
        mid <- marmod$Model$coefficients[2, 1] %>%
            round(dig)
        se <- marmod$Model$coefficients[2, 2]
        low <- (mid + qnorm(qs) * se) %>%
            round(dig)
        up <- (mid + qnorm(1 - qs) * se) %>%
            round(dig)
        mycoef <- paste0(mid, " (", low, " - ", up, ")")
    } else {
        mytable <- marmod$sigConf
        mid <- mytable[2, 1] %>% round(dig)
        low <- mytable[2, 5] %>% round(dig)
        up <- mytable[2, 6] %>% round(dig)
        mycoef <- paste0(mid, " (", low, " - ", up, ")")
    }
    return(mycoef)
}

marcoefnum <- function(marcoefficient) {
    tmp <- gsub(marcoefficient, pattern = ")", replacement = "", fixed = T)
    tmp <- gsub(tmp, pattern = "(", replacement = "", fixed = T)
    tmp <- gsub(tmp, pattern = " - ", replacement = " ", fixed = T)
    strsplit(tmp, split = " ", fixed = T)
}

marplot <- function(marmod) {
    ggplot(
        filter(marmod$data, S > 0, A > 0),
        aes(y = S, x = A)
    ) +
        geom_point(color = "black") +
        scale_x_log10() +
        scale_y_log10() +
        stat_smooth(method = "glm", fill = "lightgrey", color = "grey") +
        labs(y = "mutations", x = "area") +
        coord_cartesian(ylim = range(marmod$data$S) + 1)
}

marpseudoR2 <- function(marmod, newdata) {
    yi <- (newdata[, 2]) / max(newdata[, 2])
    ymean <- mean(yi)
    yhat <- (newdata[, 1] / max(newdata[, 1]))^marmod$par[2]
    SSreg <- sum((yhat - ymean)^2)
    SStot <- sum((yi - ymean)^2)
    SSres <- sum((yi - yhat)^2)
    R2 <- (SSres / SStot)
    # R2 = 1-(SSreg/SStot)
    R2
}
marR2 <- function(marmod, newdata) {
    yi <- (newdata[, 2])
    ymean <- mean(yi)
    yhat <- (sar_pred(marmod, newdata[, 1])$Prediction)
    # SSreg= sum( (yhat-ymean)^2 )
    # SStot= sum( (yi-ymean)^2 )
    # R2 =(SSreg/SStot)
    # R2
    summary(lm(yi ~ yhat))$r.squared
}

xcoef <- function(xvec, dig = 3) {
    mid <- xvec[1] %>% round(dig)
    low <- xvec[2] %>% round(dig)
    up <- xvec[3] %>% round(dig)
    mycoef <- paste0(mid, " (", low, " - ", up, ")")
}
extinctioncoef <- function(mymarcoef, ax) {
    if (is.character(mymarcoef)) {
        1 - ax^c(fn(unlist(marcoefnum(mymarcoef))))
    } else if (is.numeric(mymarcoef)) {
        1 - ax^mymarcoef
    } else {
        message("unknown data type calculate x coef, giving NA")
        NA
    }
}
