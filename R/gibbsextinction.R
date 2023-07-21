RaxIUCN <- function(N) {
    # 1 Extinct (or extinct in wiled), Critically Endangered, Endangered, Vulnerable, no concern
    ps <- c(0.0022539583 + 0.0007759529, 0.0863524673, 0.1587562584, 0.1562805993, 0.5955807638)
    ax_end <- c(100, 99, 79, 49, 29)
    ax_start <- c(99, 80, 50, 30, 0)
    drawcategory <- sample(1:5, N, prob = ps, replace = T)
    runif(N, min = ax_start[drawcategory], max = ax_end[drawcategory])
}
Rz <- function(N, mar_means, mar_ci95) {
    mm <- mean(mar_means, na.rm = T)
    ss <- (1 / sum(!is.na(mar_ci95))) * sqrt(sum(((mar_ci95 - mm) / 1.96)^2, na.rm = T))
    rnorm(N, mm, ss)
}
gibbsextinction <- function(iter = 350000, mar_means, mar_ci95) {
    ax <- RaxIUCN(iter)
    Rz <- Rz(iter, mar_means, mar_ci95)
    mx <- 1 - (1 - ax)^Rz
    return(mx)
}

# Rax<-function(IUCN=NULL){
#   if(is.null(IUCN)){
#     load("tmpobjects/IUCN.rda")
#     ax_start<-c(0,)
#   }
#   pro<-IUCN[,4]/sum(IUCN[,4])
#   iucnsample<- sapply(1:nrow(IUCN), function(i){
#     j=ifelse(i==1,1,i-1)
#     runif(IUCN[i,"Number_plant_spp"],min = IUCN[i,"Area_reduction_criteria"], max=IUCN[j,"Area_reduction_criteria"]) # this prob wrong
#   })
#
# }

# iucnsample<- sapply(1:nrow(IUCN), function(i){
#   j=ifelse(i==1,1,i-1)
#   runif(IUCN[i,"Number_plant_spp"],min = IUCN[i,"Area_reduction_criteria"], max=IUCN[j,"Area_reduction_criteria"]) # this prob wrong
# })
# iucnsample=fn(unlist(iucnsample))
# # hist(fn(unlist(iucnsample)))
# iucn_quantiles<-quantile(iucnsample, probs=c(.25,.5,.75))
