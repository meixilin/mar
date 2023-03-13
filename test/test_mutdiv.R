# TAKEAWAY: 
# Mon Mar 13 00:12:38 2023
# L in mutdiv is already considering all the SNPs, no need to change. 
# The increase in theta pi is due to the decreased number of samples. 


# run the original functions --------
rm(list = ls())
species = 'joshua'
dtpath = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/mutationarearelationship/mar/tmpobjects/'
load(paste0(dtpath, '/genemaps-', species, '.rda'))
set.seed(7)
yy0 = MARextinction_random(genemaps)

# load data --------
rm(list = ls()) # keey the yy0 though
require(mar)

species = 'joshua'
dtpath = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/mutationarearelationship/mar/tmpobjects/'
load(paste0(dtpath, '/genemaps-', species, '.rda'))

# extract data --------
raster_samples<-genemaps[[1]]
raster_mutmaps<-genemaps[[2]]
rest_mutmaps<-raster_mutmaps

# setup my mutdiv --------
mutdiv<-function(raster_samples, raster_mutmaps,rest_mutmaps, rasterN = integer(0)){
  require(raster)
  # Get the number of samples
  N=sum(fn(values(raster_samples)),na.rm=T)
  if (length(rasterN) == 0) {
    rasterN = N
  } 
  # freqs
  P<-fn(apply(values(raster_mutmaps), 2, function(cells) sum(cells>0,na.rm=T) )) / rasterN
  # at least one cell has to have a presence of a mutation, and sum over
  M_<-fn(apply(values(raster_mutmaps), 2, function(cells) any(cells>0)) )
  M_[is.na(M_) ]<-0
  M<-sum(M_)
  # find endemisms
  E_<-apply(values(rest_mutmaps), 2, function(cells) any(cells >0)) 
  E_[is.na(E_) ]<-0
  table(M_, E_)
  E<-sum(M_ & !E_)
  # Get the number of SNPs for the sample
  L<-dim(raster_mutmaps)[3]
  # compute diversity, Theta Waterson
  if(N>0 & M >0){
    theta<-M/(Hn(rasterN)*L)
    thetapi<-sum(2*P*(1-P),na.rm=T)/L
  }else{
    theta=0
    thetapi=0
  }
  # area taking into account only grid cells with data
  # asub= sum(raster_samples[] > 0, na.rm = T) * (res(raster_samples)[1]*res(raster_samples)[2])
  asub<-areaofraster(raster_samples)
  
  # area based on simple square
  a= dim(raster_samples)[1] * res(raster_samples)[1] * dim(raster_samples)[2] * res(raster_samples)[2] 
  # return
  return(data.frame(thetaw=theta, pi=thetapi,M=M, E=E,N=N, rasterN=rasterN, a=a, asub=asub))
}

mutdiv(raster_samples, raster_mutmaps,rest_mutmaps)
mutdiv(raster_samples, raster_mutmaps,rest_mutmaps, rasterN = 100)

# MARextinction --------
xfrac = 0.01
MARextinction_random2<-function(genemaps,xfrac=0.01,centerfun=median, debug=FALSE){
  require(raster)
  require(dplyr)
  raster_samples<-genemaps[[1]]
  raster_mutmaps<-genemaps[[2]]
  rest_mutmaps<-raster_mutmaps
  # get the  present locations
  gridpresent<-which(apply(values(raster_mutmaps),1,function(x) any(!is.na(x)))==TRUE)
  A=length(gridpresent)
  Astart=A
  xstep=ceiling(xfrac*A)
  # iterate
  listres<-list()
  # calculate original diversity
  listres<-c(listres,
             list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)))
  # extract the original rasterN
  rasterN = listres[[1]]$rasterN
  while(A > 1){
    # extinct some grids
    toextinct<-sample(gridpresent,xstep,replace = TRUE)
    values(raster_mutmaps)[toextinct,]<-NA
    values(raster_samples)[toextinct]<-NA
    values(rest_mutmaps)[toextinct,]<-NA
    # calculate diversity
    listres<-c(listres,
               list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps, rasterN = rasterN)))
    # recalculate area remaining
    gridpresent<-which(apply(values(raster_mutmaps),1,function(x) any(!is.na(x)))==TRUE)
    A=A-xstep
  }
  res<- listres %>% do.call(rbind,.) %>% data.frame %>%
    mutate(ax = 1-(asub/max(asub,na.rm=T)),
           mx = 1-(M/max(M,na.rm=T)))
  return(res)
}

set.seed(7)
yy = MARextinction_random2(genemaps)
yyt = yy[,1:9]
diffdf(yyt, yy0)


plot(yy$thetaw, yy0$thetaw)


