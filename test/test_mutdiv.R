# TAKEAWAY: 
# Mon Mar 13 00:12:38 2023
# L in mutdiv is already considering all the SNPs, no need to change. 
# The increase in theta pi is due to the decreased number of samples. 

# load data --------

rm(list = ls())
require(mar)

species = 'joshua'
dtpath = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/mutationarearelationship/mar/tmpobjects/'
load(paste0(dtpath, '/genemaps-', species, '.rda'))

# run mutdiv
# number of samples for each raster cell
raster_samples<-genemaps[[1]]
# for each mutation, there is a raster layer around the range of the data
raster_mutmaps<-genemaps[[2]]
rest_mutmaps<-raster_mutmaps
xfrac <-  0.1
# the first mutdiv
mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)

# run extinction --------
# get the  present locations (80 cells)
gridpresent<-which(apply(values(raster_mutmaps),1,function(x) any(!is.na(x)))==TRUE)
A=length(gridpresent)
Astart=A
xstep=ceiling(xfrac*A)
# iterate
listres<-list()
# calculate original diversity
listres<-c(listres,
           list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)))
listall <- vector('list', length = 10)
# get original L
L<-dim(raster_mutmaps)[3]
id = 1
while(A > 1) {
  # extinct some grids
  set.seed(7)
  toextinct<-sample(gridpresent,xstep,replace = TRUE)
  values(raster_mutmaps)[toextinct,]<-NA
  values(raster_samples)[toextinct]<-NA
  values(rest_mutmaps)[toextinct,]<-NA
  # calculate diversity
  listres<-c(listres,
             list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps, L)))
  listres
  listall[[id]] <- list(raster_samples, raster_mutmaps, rest_mutmaps)
  # recalculate area remaining
  gridpresent<-which(apply(values(raster_mutmaps),1,function(x) any(!is.na(x)))==TRUE)
  gridpresent
  A=A-xstep
  id = id + 1
}

res<- listres %>% do.call(rbind,.) %>% data.frame %>%
  mutate(ax = 1-(asub/max(asub,na.rm=T)),
         mx = 1-(M/max(M,na.rm=T)))

# list all --------
myid = 6
raster_samples = listall[[myid]][[1]]
plot(raster_samples)
raster_mutmaps = listall[[myid]][[2]]
rest_mutmaps = listall[[myid]][[3]]

# Get the number of samples
N=sum(fn(values(raster_samples)),na.rm=T)
# freqs (the increase is caused because of the decreased number of samples)
P<-fn(apply(values(raster_mutmaps), 2, function(cells) sum(cells>0,na.rm=T) )) / N

# only keep the cells that are not 
# at least one cell has to have a presence of a mutation, and sum over
M_<-fn(apply(values(raster_mutmaps), 2, function(cells) any(cells>0)) )
M_[is.na(M_) ]<-0
M<-sum(M_)
# find endemisms
E_<-apply(values(rest_mutmaps), 2, function(cells) any(cells >0)) 
E_[is.na(E_) ]<-0
table(M_, E_)
E<-sum(M_ & !E_)
# Sum samples across cells
N<-sum(values(raster_samples),na.rm=TRUE)
# Get the number of SNPs for the sample (options for using the same L)
# Problem is that dim(raster_mutmaps)[3] will always be the same. 
if (rasterL == -1) {
  L<-dim(raster_mutmaps)[3]
} else {
  L <- rasterL
}

# compute diversity, Theta Waterson
if(N>0 & M >0){
  theta<-M/(Hn(N)*L)
  thetapi<-sum(2*P*(1-P),na.rm=T)/L
}else{
  theta=0
  thetapi=0
}

# radial extinction --------
species = 'joshua'
centerfun = median
dtpath = '/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/mutationarearelationship/mar/tmpobjects/'
load(paste0(dtpath, '/genemaps-', species, '.rda'))

# get the normal ones first
zz = MARextinction_sim(genemaps,scheme = "radial",
                       centerfun = median,
                       xfrac=0.05)
# debug
require(raster)
require(dplyr)
raster_samples<-genemaps[[1]]
raster_mutmaps<-genemaps[[2]]
rest_mutmaps<-raster_mutmaps
# get most abundant location, or median
tmp<-xyFromCell(genemaps[[1]], which(values(genemaps[[1]]) >0))
startcoordextinction<-fn(apply(tmp,2, centerfun)) %>% t() %>% data.frame %>%
  rename(x=X1, y=X2)
# extinction cell from which we will expand extinction
id.cell <- raster::extract(genemaps[[1]],SpatialPoints(startcoordextinction), cellnumbers=TRUE)[1]
startrowcol<-rowColFromCell(genemaps[[1]], id.cell)
# get the  present locations
gridpresent<-which(apply(values(raster_mutmaps),1,function(x) any(!is.na(x)))==TRUE)
A=length(gridpresent)
Astart=A
xstep=ceiling(xfrac*A)
# get the latlong of every cell in the grid
locs = raster::as.data.frame(raster_samples,xy=TRUE)
# get distance to the extinction point
alldist<-as.matrix(dist(rbind(startcoordextinction,locs[,1:2]),method = 'euclidean'))[1,-1] %>%fn()
# iterate
listres<-list()
# calculate original diversity
listres<-c(listres,
           list(mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)))
if(debug) print(listres)
while(A > 1){ # change 0 for Astop if wanted to stop earlier
  if(debug) message("A ",A)
  # extinct some grids. get the top that are closest in distance
  toextinct<-gridpresent[which( alldist[gridpresent] < sort(alldist[gridpresent])[xstep])]
  # extinct those values
  values(raster_samples)[toextinct]<-NA
  values(raster_mutmaps)[toextinct,]<-NA
  values(rest_mutmaps)[toextinct,]<-NA
  # calculate diversity
  tmpdiv<-mutdiv(raster_samples, raster_mutmaps, rest_mutmaps)
  listres<-c(listres,list(tmpdiv))
  # recalculate area remaining
  gridpresent<-which(apply(values(raster_mutmaps),1,
                           function(x){any(!is.na(x))})
  )
  A=A-xstep
  if(debug) message("asub ",tmpdiv$asub)
  # if(debug) plot(raster_samples>100) # ** for debug# for debug
}
## end function
res<- listres %>% do.call(rbind,.) %>% data.frame %>%
  mutate(ax = 1-(asub/max(asub,na.rm=T)),
         mx = 1-(M/max(M,na.rm=T)))


# SCRATCH ========






















load('/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/pi_extinct/data/extinctionsim-new/extinctionsim-new-joshua.rda') 



zz = xsim %>% 
  dplyr::filter(type == 'inwards')
zz = zz[,1:9]

library(diffdf)
diffdf(yy,zz)

plot(zz$pi, yy$pi)
abline(a = 0, b = 1)
plot(zz$pi)
points(yy$pi, col = 'red')

View(cbind(yy$pi, zz$pi))
