##############################
setwd('~/Documents/R/snotel_regression')
library(gstat)
library(sp)
library(spdep)
library(raster)
library(ggplot2)
library(MASS)
library(R.matlab)

load('~/Documents/R/import.process.matlabdata/snoteldata.RData')
swe=cbind(swe,phvrec)
swe$mth=factor(swe$mth,c('Mar','Apr','May'))

sitecoords=readMat('~/Documents/R/import.process.matlabdata/snotel_attributes/sitecoords.mat')
slt=sitecoords[[1]][,1]
sln=sitecoords[[1]][,2]
snotel.locs=data.frame(lat=slt,long=sln)
coordinates(snotel.locs)<- ~long+lat
proj4string(snotel.locs) <- '+proj=longlat +datum=WGS84'

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

 ## REGRESSION AND STEPAIC
load('snotel_regression/recon.v2.stack.RData')
mdl.stepaic.phv=list()
mdl.wrecon=list()
mdl.predict=list()
pr.phvrcn.surface.stack=list()
pr.phv.surface.stack=list()
rcn.sig.phv=list()
# ## ## Read in all recon images and create rasterstack
# recon.v2.stack=stack()
# rv='v2.senlat.nldas.usace.no_canopy_windcorr'         
#  yr=2000
#  mth='Mar'
#  for (yr in seq(2000,2011)){
#       for (mth in c('MAR','APR','MAY')){
#            snodis_geotiffs=paste('/Volumes/WSC/SWE_SNODIS/recon/',rv,'/UpperColoradoRiver/',yr,'/tif',sep='')
#            ### Get current recon
#            im.string=paste(snodis_geotiffs,'/','01',toupper(mth),yr,'.tif',sep='')
#            r=raster(im.string,crs='+proj=longlat +datum=NAD83')
#            names(r)='recon'
#            recon.v2.stack=addLayer(recon.v2.stack,r)
#       }
# }
# 
# ## ## print(str(recon.v2.stack))
#  recon.v2.stack.scaled=recon.v2.stack
#  m=cellStats(recon.v2.stack,'mean')
#  stdev=cellStats(recon.v2.stack,'sd')
#  for(i in 1:nlayers(recon.v2.stack)){
#      recon.v2.stack.scaled[[i]]=(recon.v2.stack[[i]]-m[i])/stdev[i]
#  }
## print(str(recon.v2.stack.scaled))

#scale the x variables.
swe.sc=as.data.frame(scale(swe[,c('Latitude','Longitude','Elev','Aspect','Slope','RegionalSlope','RegionalAspect','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','NWd2ocean','SWd2ocean')]))
swe.sc=cbind(swe.sc,snotel=swe$snotel,mth=swe$mth,yr=swe$yr,date=swe$date)
## print(str(swe.sc))

samplesize=NULL
j=1
j2=1
yr=2000
mth='Mar'
for (yr in seq(2000,2011)){
     for (mth in c('Mar','Apr','May')){ 
          print(paste('iteration - ',j,': ',mth,yr,sep=''))

          swe2model=swe.sc[toupper(swe.sc$mth)==toupper(mth) & swe.sc$yr==yr,]
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

          
          mdl.stepaic.phv[[j]]=list()
          mdl.wrecon[[j]]=list()
          mdl.predict[[j]]=list()
          pr.phv.surface.stack[[j]]=stack()
          pr.phvrcn.surface.stack[[j]]=stack()
          rcn.sig.phv[[j]]=list()
          sr=list()
          ### Base model         
          ## start.param=fitdistr(swe2model$snotel,'gamma')
          mdl=glm(snotel~Latitude + Longitude+Elev+Aspect+Slope+RegionalSlope+RegionalAspect+FtprtW+Wbdiff+NWbdiff+SWbdiff+Wd2ocean+NWd2ocean+SWd2ocean,family=gaussian,data=swe2model)
          for (ryr in seq(2000,2011)){
              print(paste('ryr = ',ryr))
              print(paste('j2 = ',j2))
              r=recon.v2.stack.scaled[[j2]]
              names(r)='recon'
              var.stack=addLayer(ucophv.stack.scaled,r)
              site.recon=as.numeric(raster::extract(r,snotel.locs))
              swe2model$recon=site.recon#[z2keep]

             
### Reduce with stepaic
              ## PHV only
              if(ryr==2000){
                  mdl.stepaic.phv[[j]][[ryr-1999]]=stepAIC(mdl, scope=list(upper= ~ Latitude + Longitude + Elev+ Aspect + Slope + RegionalSlope + RegionalAspect + FtprtW + Wbdiff + NWbdiff + SWbdiff + Wd2ocean + NWd2ocean + SWd2ocean, lower=~1), direction='both',trace=F,k=log(nrow(swe2model)))
              } else {
                  mdl.stepaic.phv[[j]][[ryr-1999]]=mdl.stepaic.phv[[j]][[1]]
              }

              
### Reduce with stepaic
              ## PHV and RECON
              newformula=paste(mdl$formula[2],mdl$formula[1],mdl.stepaic.phv[[j]][[ryr-1999]]$formula[3],'+ recon',sep=' ')
              print(newformula)
              
              mdl.wrecon[[j]][[ryr-1999]]=glm(newformula,data=swe2model,family=gaussian)

              mdl.anova=anova(mdl.stepaic.phv[[j]][[ryr-1999]],mdl.wrecon[[j]][[ryr-1999]],test='F')


              
### Predict surface from regression

              if(mdl.anova$P[2]<0.05){
                  rcn.sig.phv[[j]][[ryr-1999]]=1
              } else {
                  rcn.sig.phv[[j]][[ryr-1999]]=0
              }

          #PHV and RECON surface
                  pr.phvrcn.surface=predict(var.stack,mdl.wrecon[[j]][[ryr-1999]],type='response')
                  pr.phvrcn.surface[r==0]=0
                  ## pr.surface=pr.surface*attr(snotel.sc,'scaled:scale')+attr(snotel.sc,'scaled:center')
                  pr.phvrcn.surface[pr.phvrcn.surface<0]=0
                  projection(pr.phvrcn.surface)='+proj=longlat +datum=NAD83'
                  pr.phvrcn.surface.stack[[j]]=addLayer(pr.phvrcn.surface.stack[[j]],pr.phvrcn.surface)

          #PHV only surface #[[-nlayers(var.stack)]]
                  pr.phv.surface=predict(var.stack,mdl.stepaic.phv[[j]][[ryr-1999]],type='response')
                  pr.phv.surface[r==0]=0
                  ## pr.surface=pr.surface*attr(snotel.sc,'scaled:scale')+attr(snotel.sc,'scaled:center')
                  pr.phv.surface[pr.phv.surface<0]=0
                  projection(pr.phv.surface)='+proj=longlat +datum=NAD83'
                  pr.phv.surface.stack[[j]]=addLayer(pr.phv.surface.stack[[j]],pr.phv.surface)
              


              j2=j2+3

          }
          j=j+1
          if(j%%3==1) j2=1
          if(j%%3==2) j2=2
          if(j%%3==0) j2=3
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
    
