##############################
setwd('~/Documents/R/snotel_regression')
library(sp)
library(spdep)
library(raster)
library(ggplot2)
library(DMwR)
library(R.matlab)
library(lubridate)

load('~/Documents/R/snotel_regression/snoteldata.RData')

sitecoords=readMat('~/Documents/R/import.process.matlabdata/snotel_attributes/sitecoords.mat')
slt=sitecoords[[1]][,1]
sln=sitecoords[[1]][,2]
snotel.locs=data.frame(lat=slt,long=sln)
coordinates(snotel.locs)<- ~long+lat
proj4string(snotel.locs) <- '+proj=longlat +datum=WGS84'

load('ucophv.RData')
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

 ## REGRESSION AND STEPAIC
mdl.stepaic.phv=list()
mdl.wrecon=list()
mdl.predict=list()
pr.phvrcn.surface.stack=list()
pr.phv.surface.stack=list()
rcn.sig.phv=list()

load('recon.v2.stack.RData')

## ## Read in all recon images and create rasterstack
## recon.v2.stack=stack()
## rv='v2.senlat.nldas.usace.no_canopy_windcorr'         
##  yr=2000
## j=1
##  for (yr in seq(2000,2011)){
##     rswe_maxdoy=as.numeric(strftime(strptime(paste(yr,'0301'),format='%Y%m%d',tz='MST'),'%j'))
##     rswe_enddoy=as.numeric(strftime(strptime(paste(yr,'0831'),format='%Y%m%d',tz='MST'),'%j'))
##  for (dofy in seq(rswe_maxdoy,rswe_enddoy)){ 
##         print(paste0('iteration ',yr,' ',dofy,': '))
##         snodis_geotiffs=paste('/Volumes/WSC/SWE_SNODIS/recon/',rv,'/UpperColoradoRiver/',yr,'/swe_tif',sep='')
## ### Get current recon
##         getdoy=sprintf("%003d",dofy)
##         wdt=strftime(as.POSIXct(paste0(getdoy,yr),'%j%Y',tz='MST'),'%d%b%Y',tz='MST')
##         im.string=paste(snodis_geotiffs,'/',toupper(wdt),'.tif',sep='')
##         r=raster(im.string,crs='+proj=longlat +datum=NAD83')
##         names(r)=toupper(wdt)
##         ## names(recon.v2.stack[[j]])=toupper(wdt)
##         recon.v2.stack=addLayer(recon.v2.stack,r)
##         j=j+1
##     }
## }
#
## Scale values in space and save to a file.
##  recon.v2.stack.scaled=recon.v2.stack
##  m=cellStats(recon.v2.stack,'mean')
##  stdev=cellStats(recon.v2.stack,'sd')
##  for(i in 1:nlayers(recon.v2.stack)){
##      recon.v2.stack.scaled[[i]]=(recon.v2.stack[[i]]-m[i])/stdev[i]
##  }
## ## print(str(recon.v2.stack.scaled))
## save(file='recon.v2.stack.RData',list=c('recon.v2.stack'))#,'recon.v2.stack.scaled'))




#scale the x variables.
phv_scatt=scale(phvrec)
phvrec.sc=as.data.frame(phv_scatt)
swe=data.frame(swe)
swe$doy=as.numeric(strftime(swe$date,'%j'))

samplesize=NULL
yr=2000
for (yr in seq(2000,2011)){
    rswe_maxdoy=as.numeric(strftime(strptime(paste(yr,'0301'),format='%Y%m%d',tz='MST'),'%j'))
    rswe_enddoy=as.numeric(strftime(strptime(paste(yr,'0831'),format='%Y%m%d',tz='MST'),'%j'))
    for (dofy in seq(rswe_maxdoy,rswe_enddoy){ 
        print(paste0('iteration - ',yr,' ',dofy,': '))
        
        swesub=subset(swe,doy==dofy & year(date)==yr)
        
        snotel_scatt=scale(swesub$snotel)
        swesub$snotel=snotel_scatt
        swe2model=cbind(swesub,phvrec.sc)
        z2keep=which(swe2model$snotel>0)
        if(length(z2keep)==0){
            next
        }
        swe2model=swe2model[z2keep,]
        ## samplesize[j]=nrow(swe2model)
        ## swe2model$snotel[swe2model$snotel==0]=0.0001

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

          
          mdl.stepaic.phv=list()
          mdl.wrecon[[j]]=list()
          mdl.predict[[j]]=list()
          pr.phv.surface.stack[[j]]=stack()
          pr.phvrcn.surface.stack[[j]]=stack()
          rcn.sig.phv=matrix(NA,nrow=12,ncol=12)
          sr=list()
          ### Base model         
          ## start.param=fitdistr(swe2model$snotel,'gamma')
          mdl=glm(snotel~Lat + Long+Elev+Eastness+Northness+Slope+RegionalSlope+RegionalEastness+RegionalNorthness+FtprtW+Wbdiff+NWbdiff+SWbdiff+Wd2ocean+NWd2ocean+SWd2ocean,family=gaussian,data=swe2model)

        for (ryr in seq(2000,2011)){
            print(paste('ryr = ',ryr))
            rswe_maxdoy=as.numeric(strftime(strptime(paste(ryr,'0301'),format='%Y%m%d',tz='MST'),'%j'))
            rswemax=subset(swe[,'recon'],year(swe$date)==ryr & swe$doy==rswe_maxdoy)
            rswemax.sc=scale(rswemax)
            
            if(dofy < rswe_maxdoy){
                siterecon=rswemax.sc
            } else {
                siterecon=subset(swe$recon,swe$doy==dofy & year(swe$date)==ryr)
            }
            swe2model$recon=siterecon[z2keep]
                              
### Reduce with stepaic
            ## PHV only
            if(ryr==2000){
                mdl.stepaic.phv[[j]]=stepAIC(mdl, scope=list(upper= ~ Lat + Long + Elev+ Eastness + Northness + Slope + RegionalSlope + RegionalEastness + RegionalNorthness + FtprtW + Wbdiff + NWbdiff + SWbdiff + Wd2ocean + NWd2ocean + SWd2ocean, lower=~1), direction='both',trace=F,k=log(nrow(swe2model)))
            } 
            
            
### Reduce with stepaic
            ## PHV and RECON
            newformula=paste(mdl$formula[2],mdl$formula[1],mdl.stepaic.phv[[j]]$formula[3],'+ recon',sep=' ')
            print(newformula)
            
            mdl.wrecon[[j]][[ryr-1999]]=glm(newformula,data=swe2model,family=gaussian)
            mdl.anova=anova(mdl.stepaic.phv[[j]],mdl.wrecon[[j]][[ryr-1999]],test='F')
              
### Predict surface from regression
              if(mdl.anova$P[2]<0.05){
                  rcn.sig.phv[j,ryr-1999]=1
              } else {
                  rcn.sig.phv[j,ryr-1999]=0
              }


#PHV only surface #[[-nlayers(var.stack)]]
            r=recon.v2.stack[[yr-1999]]#using this to mask so use yr of simulation
            var.stack=ucophv.stack
            pr.phv.surface=predict(var.stack,mdl.stepaic.phv[[j]][[ryr-1999]],type='response')
            pr.phv.surface[r==0]=0
            pr.surface=pr.surface*attr(swe2model$snotel,'scaled:scale')+attr(swe2model$snotel,'scaled:center')

            pr.phv.surface[pr.phv.surface<0]=0
            projection(pr.phv.surface)='+proj=longlat +datum=NAD83'
            pr.phv.surface.stack[[j]]=addLayer(pr.phv.surface.stack[[j]],pr.phv.surface)
            
#PHV and RECON surface
            var.stack=stack(ucophv.stack,recon.v2.stack[[ryr-1999]])#predictor recon surface changes
            pr.phvrcn.surface=predict(var.stack,mdl.wrecon[[j]][[ryr-1999]],type='response')
            pr.phvrcn.surface[r==0]=0
            pr.surface=pr.surface*attr(swe2model$snotel,'scaled:scale')+attr(swe2model$snotel,'scaled:center')
            pr.phvrcn.surface[pr.phvrcn.surface<0]=0
            projection(pr.phvrcn.surface)='+proj=longlat +datum=NAD83'
            pr.phvrcn.surface.stack[[j]]=addLayer(pr.phvrcn.surface.stack[[j]],pr.phvrcn.surface)
          }
          j=j+1
      }
 }

## save.image(file='latest.image.RData')
save(file='latest.model/mdl.stepaic.gauss.snotelunscaled.xrecon.Ftest.RData',list=c('mdl.stepaic.phv','mdl.wrecon','rcn.sig.phv','recon.v2.stack'))
save(file='latest.model/prediction.surfaces.gauss.snotelunscaled.xrecon.Ftest.RData',list=c('pr.phvrcn.surface.stack','pr.phv.surface.stack','recon.v2.stack'))
     
     
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
    
