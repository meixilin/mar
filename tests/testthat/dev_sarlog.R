library(mar)

mardf <- MARsampling(gm1001g)


plot(N ~ A, data = mardf)
plot(N ~ Asq, data = mardf)

summary(lm(N ~ A, data = mardf))
summary(lm(N ~ Asq, data = mardf))

Mtype = "thetaw"
Atype = "A"

# remove NA or zero data
tmpdf = mardf[, c(Atype, Mtype)]
tmpdf = tmpdf[(tmpdf[,Mtype] > 0 & !is.na(tmpdf[,Mtype])), ]
# run MAR analyses
mar <- sars::sar_multi(tmpdf, obj = c('power', "linear", "loga", "koba"))
mar <- sars::sar_multi(tmpdf, obj = c('power', "loga"))

lapply(mar, function(x) {summary(x)$AIC})
plot(mar$power)
plot(mar$loga)
plot(mar$koba)
# sars::display_sars_models()
plot(mar[c(1,4)])
summary(mar$power)
summary(mar$koba)

theta = 3.2576e+03
p = 5.130e-04
z = 1.2360e+05
theta/p

theta = 2.5345e+03
p = 0.21787
z = 1.7535e+02
theta/p

c1 = 2534.4560
z1 = 175.3528

c2 = 198.4825188
z2 = 0.4558545

curve(c1 * log(1 + x/z1), from = 0, to = 3000)
curve(c2 * x^z2, col = 'blue', add = T)

curve(c1 * (1/z1/(1 + x/z1)), from = 0, to = 3000)
curve(c2 * (x^(z2 - 1) * z2), col = 'blue', add = T)

D(expression(c*log(1+x/z)), name = 'x')
D(expression(c*x^z), name = 'x')
