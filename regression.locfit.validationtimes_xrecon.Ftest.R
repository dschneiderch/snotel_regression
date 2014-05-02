##############################
#setwd('.')
library(lubridate)
library(sp)
library(spdep)
library(raster)
library(ggplot2)
library(MASS)
library(R.matlab)
require(locfit)
require(rgdal)

load('modeldata.validationtimes.RData')

sitecoords=readMat('sitecoords.mat')
slt=sitecoords[[1]][,1]
sln=sitecoords[[1]][,2]
snotel.locs=data.frame(lat=slt,long=sln)
coordinates(snotel.locs)<- ~long+lat
proj4string(snotel.locs) <- '+proj=longlat +datum=WGS84'
projstr='+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0'#usgs albers for conus
snotel.locs.usgs=spTransform(snotel.locs,CRS(projstr))
spdll=spDists(snotel.locs,longlat=T)#distances in km


load('ucophv.stack.RData')
#ucophv.stack.scaled=ucophv.stack
#m=cellStats(ucophv.stack,'mean')
#stdev=cellStats(ucophv.stack,'sd')
#for(i in 1:nlayers(ucophv.stack)){  
#     ucophv.stack.scaled[[i]]=(ucophv.stack[[i]]-m[i])/stdev[i]
#   }
    
##*********************
# # Map domain and snotel stations
# #*** For Upper CO River Basin domain
upleft=c(43.75,-112.25)
lowright=c(33,-104.125)
lon.low=-112.25
lon.high=-104.125
lat.low=33
lat.high=43.75
lim=c(lon.low, lon.high, lat.low, lat.high)
##


load('recon.v2.stack.validation.RData')
## #REGRESSION AND STEPAIC
##  snotelrecon=array(NA,dim=c(length(dates2keep),12,237))
##  snotelrecon.sc=snotelrecon
## ## Read in all recon images and create dataframe of snotel pixel values
## rv='v2.senlat.nldas.usace.no_canopy_windcorr'         
## for(i in 1:length(dates2keep)){
##      for(yri in 2000:2011){
##          yr=year(dates2keep[i])       
##          mth=month(dates2keep[i],label=T,abbr=T)
##          dy=sprintf('%02d',day(dates2keep[i]))
##          snodis_geotiffs=paste('/Volumes/WSC/SWE_SNODIS/recon/',rv,'/UpperColoradoRiver/',yr,'/tif',sep='')
## ### Get current recon
##          im.string=paste(snodis_geotiffs,'/',dy,toupper(mth),yr,'.tif',sep='')
##          r=raster(im.string,crs='+proj=longlat +datum=NAD83')
##          recon.v2.stack[[i]]=r
##          names(r)='recon'
##          snotelrecon[i,yri-1999,]=as.numeric(raster::extract(r,snotel.locs))
##          snotelrecon.sc[i,yri-1999,]= scale(snotelrecon[i,yri-1999,])
##      }
##  }
#
## recon.v2.stack.scaled=recon.v2.stack
## m=cellStats(recon.v2.stack,'mean')
## stdev=cellStats(recon.v2.stack,'sd')
## for(i in 1:nlayers(recon.v2.stack)){
##      recon.v2.stack.scaled[[i]]=(recon.v2.stack[[i]]-m[i])/stdev[i]
##  }

## save(file='recon.v2.stack.validation.RData',list=c('recon.v2.stack','recon.v2.stack.scaled','snotelrecon','snotelrecon.sc','rv','dates2keep'))



#x variables will get scaled in locfit call
swe.sc=as.data.frame(scale(snotel2model[,c('Latitude','Longitude','Elev','Aspect','Slope','RegionalSlope','RegionalAspect','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','NWd2ocean','SWd2ocean')]))
swe.sc=cbind(swe.sc,snotel=snotel2model$swe,mth=snotel2model$mth,yr=snotel2model$yr,date=snotel2model$date)
# create dataframe from ucophv.stack with which to predict
data2predict=data.frame(rasterToPoints(ucophv.stack.scaled))


mdl.wrecon=list()
mdl.wphv=list()
pr.phvrcn.surface.stack=list()
pr.phv.surface.stack=list()
rcn.sig.phv=list()
#
i=1
## for (i in 1:length(dates2keep)){
for(i in 17:18){
    yr=year(dates2keep[i])
    mth=month(dates2keep[i],label=T,abbr=T)
    dy=day(dates2keep[i])
    print(paste('iteration - ',i,': ',dy,mth,yr,sep=''))
    
    swe2model=swe.sc[swe.sc$date==dates2keep[i],]
    ## z2keep=which(swe2model$snotel>0)
    ## swe2model=swe2model[z2keep,]
    ## samplesize[j]=nrow(swe2model)
    swe2model$snotel[swe2model$snotel==0]=0.0001
    snotel.sc=scale(swe2model$snotel)
    ## N=nrow(swe2model)
    ## swe2model$snotel=swe2model$snotel+runif(N,0.001,0.002)
    
    ## fn=paste('plots','/','snotel.datadistribution.','01',mth,yr,'.pdf',sep='') #filename
    ## pdf(fn)
    ## x11=swe2model$snotel
    ## hist(x11,probability=T)
    ## sm.density(x11,add=T)
    ## gparam=fitdistr(x11,'gamma')
    ## lines(x11[order(x11)], dgamma(x11[order(x11)], shape=gparam$estimate[[1]],rate=gparam$estimate[[2]]),col='red')
    ## lines(x11[order(x11)], dnorm(x11[order(x11)], mean=mean(swe2model$snotel), sd=sd(swe2model$snotel)),col='blue')
    ## legend('topright',legend=c('red=gamma','blue=gauss','black=density'))
    ## dev.off()
         
    mdl.wphv[[i]]=list()
    mdl.wrecon[[i]]=list()
    pr.phv.surface.stack[[i]]=stack()
    pr.phvrcn.surface.stack[[i]]=stack()
### Base model         
    ## start.param=fitdistr(swe2model$snotel,'gamma')
    mdl=locfit(snotel~lp(Latitude + Elev + RegionalSlope + Wbdiff + NWbdiff + SWbdiff,nn=0.3,deg=1),family='gamma',data=swe2model)
    mdl.wphv[[i]]=mdl
    for (ryr in seq(2000,2011)){
        print(paste('ryr = ',ryr))
        r=recon.v2.stack.scaled[[i]]
        names(r)='recon'
        data2predict$recon=data.frame(rasterToPoints(r))$recon
        ## site.recon=as.numeric(raster::extract(r,snotel.locs))
        ## swe2model$recon=site.recon#[z2keep]
        swe2model$recon=snotelrecon[i,ryr-1999,]
              
        mdl.recon=locfit(snotel~lp(Latitude + Elev + RegionalSlope + Wbdiff + NWbdiff + SWbdiff + recon,nn=0.3,deg=1),data=swe2model,family='gamma')
        mdl.wrecon[[i]][[ryr-1999]]=mdl.recon
              
### Predict surface from regression
          #PHV and RECON surface
        pr.phvrcn.surface=predict(mdl.recon,newdata=data2predict)
        pr.phvrcn.surface[data2predict$recon==0]=0
#              pr.surface=pr.phvrcn.surface*attr(snotel.sc,'scaled:scale')+attr(snotel.sc,'scaled:center')
        pr.phvrcn.surface[pr.phvrcn.surface<0]=0
#                  projection(pr.phvrcn.surface)='+proj=longlat +datum=NAD83'
        pr.phvrcn.surface=data.frame(cbind(x=data2predict$x,y=data2predict$y,swe=pr.phvrcn.surface))
     pr.phvrcn.surface.stack[[i]]=addLayer(pr.phvrcn.surface.stack[[i]],rasterFromXYZ(pr.phvrcn.surface,crs='+proj=longlat +datum=NAD83'))
       ## pr.phvrcn.surface.stack[[i]][[ryr-1999]]=pr.phvrcn.surface
          }

         #PHV only surface 
    pr.phv.surface=predict(mdl,newdata=data2predict)
    pr.phv.surface[data2predict$recon==0]=0
 #   pr.phv.surface=pr.phv.surface*attr(snotel.sc,'scaled:scale')+attr(snotel.sc,'scaled:center')
    pr.phv.surface[pr.phv.surface<0]=0
 #                 projection(pr.phv.surface)='+proj=longlat +datum=NAD83'
    pr.phv.surface=data.frame(cbind(x=data2predict$x,y=data2predict$y,swe=pr.phv.surface))
    pr.phv.surface.stack[[i]]=addLayer(pr.phv.surface.stack[[i]],rasterFromXYZ(pr.phv.surface,crs='+proj=longlat +datum=NAD83'))
    ## pr.phv.surface.stack[[i]][[ryr-1999]]=pr.phv.surface
}



## library(RColorBrewer)
## surface=pr.phv.surface
## surface$cuts=cut(surface$swe,seq(0,max(surface$swe),0.05),include.lowest=F)
## ggplot(surface)+geom_raster(aes(x,y,fill=cuts))+scale_fill_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(surface$cuts))),drop=F)+theme_bw()


## save.image(file='latest.image.RData')
save(file='output/prediction.locfit.surfaces.gamma.snotelunscaled.xrecon.validation.URG.RData',list=c('pr.phvrcn.surface.stack','pr.phv.surface.stack','recon.v2.stack','dates2keep'))


urg.phv.surface.stack=list(pr.phv.surface.stack[[17]],pr.phv.surface.stack[[18]])
urg.phvrcn.surface.stack=list(pr.phvrcn.surface.stack[[17]],pr.phvrcn.surface.stack[[18]])
urg.recon.v2.surface.stack=list(recon.v2.stack[[17]],recon.v2.stack[[18]])
urg.dates=dates2keep[17:18]

## csu.phv.surface.stack=list()
## csu.phvrcn.surface.stack=list()
## csu.recon.v2.surface.stack=list()
## for(i in 1:16){
## 	csu.phv.surface.stack[[i]]=pr.phv.surface.stack[[i]]
## 	csu.phvrcn.surface.stack[[i]]=pr.phvrcn.surface.stack[[i]]
## 	csu.recon.v2.surface.stack[[i]]=recon.v2.stack[[i]]
## }
## csu.dates=dates2keep[1:16]

## #dates2keep has dates for glv surveys in reverse chrono order. fixing it for the glv surface file
## glv.phv.surface.stack=list()
## glv.phvrcn.surface.stack=list()
## glv.recon.v2.surface.stack=list()
## for(i in 30:19){
## 	glv.phv.surface.stack[[i-18]]=pr.phv.surface.stack[[i]]
## 	glv.phvrcn.surface.stack[[i-18]]=pr.phvrcn.surface.stack[[i]]
## 	glv.recon.v2.surface.stack[[i-18]]=recon.v2.stack[[i]]
## }
## glv.dates=dates2keep[30:19]


## save(file='output/urg.surfaces.validation.RData',list=c('urg.phv.surface.stack','urg.phvrcn.surface.stack','urg.recon.v2.surface.stack','urg.dates'))
## save(file='output/csu.surfaces.validation.RData',list=c('csu.phv.surface.stack','csu.phvrcn.surface.stack','csu.recon.v2.surface.stack','csu.dates'))
## save(file='output/glv.surfaces.validation.RData',list=c('glv.phv.surface.stack','glv.phvrcn.surface.stack','glv.recon.v2.surface.stack','glv.dates'))



     
     ## ###### Load mdl.stepaic or mdl.stepaic.gamma from R/snotel_regression
     ## # either gaussian or gamma
     ## load('mdl.stepaic.RData') #
     ## for(i in seq(1:39)){
     ##      qplot(predict(mdl.stepaic[[i]],type='response'),resid(mdl.stepaic[[i]]))
     ##      ggsave(paste('plots/resids/resid.gauss.',i,'.png',sep=''))
     ##      #
     ##      x11=residuals(mdl.stepaic[[i]])
     ##      png(paste('plots/resids/resids.hist.gauss.',i,'.png',sep=''))
     ##      hist(x11,probability=T) #some positive skew. not terrible
     ##      sm.density(x11,add=T) 
     ##      lines(x11[order(x11)], dnorm(x11[order(x11)], mean=mean(x11), sd=sd(x11)), col=       "red")
     ##      dev.off()
     ## }
     
     
##      rm(list=ls())


## source('plot.bestpredictedsurfaces.R')
     
## source('xvalidated.xrecon.regression.skill.R')


## ii=1
## for( ii in seq(1,nlayers(ucophv.stack.scaled))){
##     print(cellStats(ucophv.stack.scaled[[ii]],range) )
##     ## rp=data.frame(rasterToPoints(ucophv.stack.scaled[[ii]]))
##     ## colnames(rp) <- c('lat','long','var')
##     ## rn=names(ucophv.stack.scaled[[ii]])
##     ## g=ggplot()+geom_raster(data=rp,aes(x=long,y=lat,fill=var))+labs(title=rn)
##     ## ggsave(plot=g,filename=paste('plots/ucovarplots/',rn,'.pdf',sep=''))
## }
    
