
# x="test_N_10_n_1_s_0.00_mig_0.1.out"
# parseslim<-function(x){
#   tmp<-x
#   tmp<-gsub(tmp,pattern = "test_",replacement = "")
#   tmp<-gsub(tmp,pattern = ".out",replacement = "")
#   tmp<-strsplit(tmp,split = "_",fixed = T)  [[1]]
#   tmp2<-tmp[c(2,4,6,8)]
#   names(tmp2)<-tmp[c(2,4,6,8)-1]
#   tmp2
# }
parseparams <- function(x = "test_N_100_n_1_s_0.1_mig_1e-1_rep_1.out") {
    x <- gsub(x, pattern = "test_", replacement = "")
    x <- gsub(x, pattern = ".out", replacement = "")
    res_ <- strsplit(x, split = "_", fixed = T)[[1]]
    resul <- res_[c(2, 4, 6, 8)]
    names(resul) <- res_[c(2, 4, 6, 8) - 1]
    return(resul)
}
getz <- function(x, path = "../dataslim/") {
    tmp <- read.table(paste0(path, "/", x), header = T)
    sars::sar_power(tmp)
}

ReadLastLinesSlim <- function(x, n) {
    tmp <- system(paste("cat", x), intern = T)
    tmp <- tail(tmp, 15)
    diditrun <- any(grepl("segsites", tmp))
    if (!diditrun) {
        out <- matrix(c(NA, NA, NA, NA, NA, NA), nrow = 2)
        out <- data.frame(out)
        colnames(out) <- c("area", "segsites", "fst")
    } else {
        out <- tmp %>%
            strsplit(., "\t") %>%
            do.call(rbind, .) %>%
            data.frame()
        colnames(out) <- out[1, ]
        out <- out[-1, ]
    }
    return(out)
}
ReadLastLinesSlim2 <- function(x, n) {
    con <- file(x)
    open(con)
    out <- scan(con, n, what = "char(0)", sep = "\n", quiet = TRUE)
    diditrun <- !any(grepl("segsites", out))
    if (incorrect) {
        out <- matrix(c(NA, NA, NA, NA, NA, NA), nrow = 2)
        out <- data.frame(out)
        colnames(out) <- c("area", "segsites", "fst")
    } else {
        while (TRUE) {
            tmp <- scan(con, 1, what = "char(0)", sep = "\n", quiet = TRUE)
            if (length(tmp) == 0) {
                close(con)
                break
            }
            out <- c(out[-1], tmp)
        }
        out <- out %>%
            strsplit(., "\t") %>%
            do.call(rbind, .) %>%
            data.frame()
        if (dim(out)[2] == 3) {
            colnames(out) <- c("area", "segsites", "fst")
            out <- out[-1, ]
        } else {
            colnames(out) <- c("area", "segsites")
            out <- out[-1, ]
            out$fst <- NA
        }
    }
    return(out)
}

# ReadLastLinesSlim <- function(x,n){
#   con <- file(x)
#   open(con)
#   out <- scan(con,n,what="char(0)",sep="\n",quiet=TRUE)
#   incorrect<-any(grepl("initializeGenomicElement", out))
#   if(incorrect)
#     while(TRUE){
#       tmp <- scan(con,1,what="char(0)",sep="\n",quiet=TRUE)
#       if(length(tmp)==0) {close(con) ; break }
#       out <- c(out[-1],tmp)
#     }
#   if(is.null(out)){
#     out<-matrix(c(NA,NA,NA,NA,NA,NA),nrow=2)
#     colnames<-c("area","segsites","fst")
#     out<-data.frame(out)
#   }else{
#     out<-out %>%
#       strsplit(.,"\t") %>%
#       do.call(rbind,.) %>%
#       data.frame
#     colnames(out)<-c("area","segsites","fst")
#     out<-out[-1,]
#   }
#   return(out)
# }
