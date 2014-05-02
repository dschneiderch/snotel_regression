setwd('~/Documents/R/snotel_regression')
library(sp)
library(spdep)
library(raster)
library(ggplot2)
library(scales)
library(MASS)
library(R.matlab)
library(reshape2)
library(plyr)
library(RColorBrewer)


load('~/Documents/R/import.process.matlabdata/snoteldata.RData')
swe=cbind(swe,phvrec)

sitecoords=readMat('../import.process.matlabdata/snotel_attributes/sitecoords.mat')
slt=sitecoords[[1]][,1]
sln=sitecoords[[1]][,2]
snotel.locs=data.frame(lat=slt,long=sln)
coordinates(snotel.locs)<- ~long+lat
proj4string(snotel.locs) <- '+proj=longlat +datum=WGS84'


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

load('recon.v2.stack.RData')
## REGRESSION AND STEPAIC
## snotelrecon=array(NA,dim=c(12,237,3))
## snotelrecon.sc=snotelrecon
## ## ## Read in all recon images and create dataframe of snotel pixel values
## mth='MAR'
## yr=2000
## mthind=1
## rv='v2.senlat.nldas.usace.no_canopy_windcorr'         
## for (mth in c('MAR','APR','MAY')){
##     for (yr in seq(2000,2011)){
##          snodis_geotiffs=paste('/Volumes/WSC/SWE_SNODIS/recon/',rv,'/UpperColoradoRiver/',yr,'/tif',sep='')
##           ### Get current recon
##          im.string=paste(snodis_geotiffs,'/','01',toupper(mth),yr,'.tif',sep='')
##          r=raster(im.string,crs='+proj=longlat +datum=NAD83')
##          names(r)='recon'
##          snotelrecon[yr-1999,,mthind]=as.numeric(raster::extract(r,snotel.locs))
##          snotelrecon.sc[yr-1999,,mthind]= scale(snotelrecon[yr-1999,,mthind])
##      }
##     mthind=mthind+1
## }



swe$mth=factor(swe$mth,levels=c('Mar','Apr','May'))
swe.sc=as.data.frame(scale(swe[,c('Latitude','Longitude','Elev','Aspect','Slope','RegionalSlope','RegionalAspect','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','NWd2ocean','SWd2ocean')]))
swe.sc=cbind(swe.sc,snotel=swe$snotel,mth=swe$mth,yr=swe$yr,date=swe$date)
## print(str(swe.sc))


j=1
yr=2001
mth='MAY'
mthind=1 
cvitmax=237
pr.phvrcn=NULL
pr.phv=NULL
SSphv=NULL
SSphvrcn=NULL
SSrecon=NULL
swediff.phv=NULL
swediff.phvrcn=NULL
swediffpct.phv=NULL
swediffpct.phvrcn=NULL
#
swediff.phvmodel=array(NA,dim=c(12,3,12,237))
swediff.phvrcnmodel=array(NA,dim=c(12,3,12,237))
swediffpct.phvmodel=array(NA,dim=c(12,3,12,237))
swediffpct.phvrcnmodel=array(NA,dim=c(12,3,12,237))
#
r2.phvmodel=array(NA,dim=c(12,3,12,cvitmax))
r2.phvrcnmodel=array(NA,dim=c(12,3,12,cvitmax))
r2.reconmodel=array(NA,dim=c(12,3,12,cvitmax))
for (mth in c('MAR','APR','MAY')){ #MUST BE CAPITAL, 3 LETTER ABBREVIATIONS 
    for (yr in seq(2000,2011)){
          print(paste('iteration - ',j,': ',mth,yr,sep=''))
          print(paste('month index = ',mthind,sep=''))
          
          swe2model.all=swe.sc[toupper(swe.sc$mth)==mth & swe.sc$yr==yr,]
          swe2model.all$snotel[swe2model.all$snotel==0]=runif(1,0,0.00001)
          for (ryr in seq(2000,2011)){
              print(paste('ryr = ',ryr))
              swe2model.all$recon=as.numeric(snotelrecon.sc[ryr-1999,,mthind])
              for(dropiter in seq(1,237)){
                  z2keep=seq(1,237,1)[-dropiter]
                  swe2model=swe2model.all[z2keep,]
                  snotel.sc=scale(swe2model$snotel)
### Base model         
                ## start.param=fitdistr(swe2model$snotel,'gamma')
                  mdl=glm(snotel~Latitude + Longitude + Elev+Aspect+Slope+RegionalSlope+RegionalAspect+FtprtW+Wbdiff+NWbdiff+SWbdiff+Wd2ocean+NWd2ocean+SWd2ocean,family=gaussian,data=swe2model)
                       
### Reduce with stepaic
                  ## PHV only
                  mdl.stepaic.phv=stepAIC(mdl, scope=list(upper= ~ Latitude + Longitude + Elev+ Aspect + Slope + RegionalSlope + RegionalAspect + FtprtW + Wbdiff + NWbdiff + SWbdiff + Wd2ocean + NWd2ocean + SWd2ocean, lower=~1), direction='both',trace=F,k=log(nrow(swe2model)))
                                  
### Reduce with stepaic
                  ## PHV and RECON
                  newformula=paste(mdl$formula[2],mdl$formula[1],mdl.stepaic.phv$formula[3],'+ recon',sep=' ')
                 # print(newformula)
                  mdl.wrecon=glm(newformula,data=swe2model,family=gaussian)
                  
### Predict surface from regression
                  #PHV and RECON surface
                  pr.phvrcn[dropiter]=predict(mdl.wrecon,type='response',newdata=swe2model.all[-z2keep,])

#scale data back 
                  ## pr.phvrcn=pr.phvrcn*attr(snotel.sc,'scaled:scale')+attr(snotel.sc,'scaled:center')
                  
                  #PHV only surface
                  pr.phv[dropiter]=predict(mdl.stepaic.phv,type='response',newdata=swe2model.all[-z2keep,])

# scale data back
                  ## pr.phv=pr.phv*attr(snotel.sc,'scaled:scale')+attr(snotel.sc,'scaled:center')

              }
### Compute skill
              yobs=swe2model.all$snotel
              
              yrecon=as.numeric(snotelrecon[ryr-1999,,mthind])
              SSrecon = cor(yrecon,yobs)^2

              SSphv = cor(pr.phv,yobs)^2#1 - SSE/SSNull
              pr.phv[yrecon==0]=0
              swediff.phv = pr.phv-yobs
              swediffpct.phv = (pr.phv-yobs)/yobs
              swediffpct.phv[is.na(swediffpct.phv)]=0
              swediffpct.phv[is.infinite(swediffpct.phv)]=1000
                            
              SSphvrcn = cor(pr.phvrcn,yobs)^2 #1 - SSE/SSNull
              pr.phvrcn[yrecon==0]=0
              swediff.phvrcn = pr.phvrcn-yobs
              swediffpct.phvrcn = (pr.phvrcn-yobs)/yobs
              swediffpct.phvrcn[is.na(swediffpct.phvrcn)]=0
              swediffpct.phvrcn[is.infinite(swediffpct.phvrcn)]=1000

              
                  ## if(length(which(is.na(SSrecon)))!=0) print(paste('NA values - recon ',mth,yr,'- recon',ryr,sep=''))
         
              r2.phvmodel[yr-1999,mthind,ryr-1999,]=SSphv
              r2.phvrcnmodel[yr-1999,mthind,ryr-1999,]=SSphvrcn
              r2.reconmodel[yr-1999,mthind,ryr-1999,]=SSrecon
              swediff.phvmodel[yr-1999,mthind,ryr-1999,]=swediff.phv
              swediff.phvrcnmodel[yr-1999,mthind,ryr-1999,]=swediff.phvrcn
              swediffpct.phvmodel[yr-1999,mthind,ryr-1999,]=swediffpct.phv
              swediffpct.phvrcnmodel[yr-1999,mthind,ryr-1999,]=swediffpct.phvrcn
              
          }
          j=j+1
      } 
    mthind=mthind+1
}


## In a much smaller loop, find swediff for recon vs snotel
swediff.reconmodel=array(NA,dim=c(12,3,12,237))
swediffpct.reconmodel=array(NA,dim=c(12,3,12,237))
j=1
mthcnt=1
for (mthi in c('Mar','Apr','May')){ #MUST BE CAPITAL, 3 LETTER ABBREVIATIONS 
    for (yri in seq(2000,2011)){
          print(paste('iteration - ',j,': ',mthi,yri,sep=''))
          print(paste('month index = ',mthcnt,sep=''))
          for(ryr in seq(2000,2011)){                            #  
              swe2model.all=swe.sc[as.character(swe.sc$mth)==mthi & swe.sc$yr==yri,]
              yobs=swe2model.all$snotel
              yrecon=as.numeric(snotelrecon[ryr-1999,,mthcnt])
              print(paste('yobs: ',which(is.na(yobs)),sep=''))
              print(paste('yrecon: ',which(is.na(yrecon)),sep=''))
                                        #   
              swediff.reconmodel[yri-1999,mthcnt,ryr-1999,]=yrecon-yobs
              swediffpct.reconmodel[yri-1999,mthcnt,ryr-1999,]=(yrecon-yobs)/yobs
              swediffpct.reconmodel[yri-1999,mthcnt,ryr-1999,is.na(swediffpct.reconmodel[yri-1999,mthcnt,ryr-1999,])]=0
              swediffpct.reconmodel[yri-1999,mthcnt,ryr-1999,is.infinite(swediffpct.reconmodel[yri-1999,mthcnt,ryr-1999,])]=1000

          }
                                        #
          j=j+1
      }
    mthcnt=mthcnt+1
}


## **************
# Not necessary in drop 1. only if doing drop10
## **************

## ## Find median r2 skill from cv iterations
## ## And find mean swediff from cv iterations
## r2phvmed=array(NA,dim=c(12,3,12))
## r2phvrcnmed=array(NA,dim=c(12,3,12))
## r2reconmed=array(NA,dim=c(12,3,12))
## swediff.phv.avg=array(NA,dim=c(12,3,12,237))
## swediff.phvrcn.avg=array(NA,dim=c(12,3,12,237))
## swediffpct.phv.avg=array(NA,dim=c(12,3,12,237))
## swediffpct.phvrcn.avg=array(NA,dim=c(12,3,12,237))
## for(yr in 1:12){
##     for(mth in 1:3){
##         r2phvmed[yr,mth,]=apply(r2.phvmodel[yr,mth,,],1,median)
##         r2phvrcnmed[yr,mth,]=apply(r2.phvrcnmodel[yr,mth,,],1,median)
##         r2reconmed[yr,mth,]=apply(r2.reconmodel[yr,mth,,],1,median)
##         for(ryr in 1:12){
##             swediff.phv.avg[yr,mth,ryr,]=apply(swediff.phvmodel[yr,mth,ryr,,],2,mean,na.rm=T)
##             swediff.phvrcn.avg[yr,mth,ryr,]=apply(swediff.phvrcnmodel[yr,mth,ryr,,],2,mean,na.rm=T)
##             swediffpct.phv.avg[yr,mth,ryr,]=apply(swediffpct.phvmodel[yr,mth,ryr,,],2,mean,na.rm=T)
##             swediffpct.phvrcn.avg[yr,mth,ryr,]=apply(swediffpct.phvrcnmodel[yr,mth,ryr,,],2,mean,na.rm=T)    
##         }
##     }
## }
## ****************

# ****** Reassign names so we can use existing code.
r2phvmed=r2.phvmodel
r2phvrcnmed=r2.phvrcnmodel
r2reconmed=r2.reconmodel
#
swediff.phv.avg=swediff.phvmodel
swediff.phvrcn.avg=swediff.phvrcnmodel
swediffpct.phv.avg=swediffpct.phvmodel
swediffpct.phvrcn.avg=swediffpct.phvrcnmodel

##### Dataframe of skill scores by yr, mth, and ryr for boxplot
r2df=data.frame(mth=factor(c('Mar','Apr','May'),levels=c('Mar','Apr','May')),yr=c(rep(2000,3),rep(2001,3),rep(2002,3),rep(2003,3),rep(2004,3),rep(2005,3),rep(2006,3),rep(2007,3),rep(2008,3),rep(2009,3),rep(2010,3),rep(2011,3)))
#
r2df=cbind(r2df,data.frame(reconyear=c(rep(2000,nrow(r2df)),rep(2001,nrow(r2df)),rep(2002,nrow(r2df)),rep(2003,nrow(r2df)),rep(2004,nrow(r2df)),rep(2005,nrow(r2df)),rep(2006,nrow(r2df)),rep(2007,nrow(r2df)),rep(2008,nrow(r2df)),rep(2009,nrow(r2df)),rep(2010,nrow(r2df)),rep(2011,nrow(r2df)))))
#
r2df=r2df[order(r2df$reconyear,r2df$mth),]
r2df=cbind(r2df,r2phv=as.vector(r2phvmed),r2phvrcn=as.vector(r2phvrcnmed),r2recon=as.vector(r2reconmed))
#
#Find the mean r2 for phv and recon. Each recon year has cvitmax iterations for each timestep. the mean is used to represent r2 skill for a given timestep since phv and recon modesl don't change with the rcnyear used.
for(iyr in 2000:2011){
    for(imth in c('Mar','Apr','May')){
        r2df[r2df$mth==imth & r2df$yr==iyr,'r2phv']=apply(r2df[imth==r2df$mth & iyr==r2df$yr,c('r2phv','r2recon')],2,mean)[1]
        ## r2df[imth==r2df$mth & iyr==r2df$yr,'r2recon']=apply(r2df[imth==r2df$mth & iyr==r2df$yr,c('r2phv','r2recon')],2,mean)[2]
    }
}

tmp=melt(r2df,id=c('mth','yr','reconyear'))
ddply(tmp,'variable',summarize, mean(value))


rcnyr.phvrcn=matrix(NA,nrow=12,ncol=3)
rcnyr.recon=matrix(NA,nrow=12,ncol=3)
for(yrd in 2000:2011){
    mthi=1
    for(mthd in c('Mar','Apr','May')){
        r2.max=ddply(.data=r2df,c('mth','yr'),summarize,which.max(r2phvrcn),max(r2phvrcn,na.rm=T))
        r2.max.recon=ddply(.data=r2df,c('mth','yr'),summarize,which.max(r2recon),max(r2recon,na.rm=T))
        names(r2.max)=c('mth','yr','modelnum','maxr2')
        names(r2.max.recon)=c('mth','yr','modelnum','maxr2')
        r2.max$reconyear=r2.max$modelnum+1999
        r2.max.recon$reconyear=r2.max.recon$modelnum+1999
        r2.max=r2.max[order(r2.max$yr),] #places it in order mth, yr same as the stack. assumes factor mth is ordered correctly.
        r2.max.recon=r2.max.recon[order(r2.max.recon$yr),] #places it in order mth, yr same as the stack. assumes factor mth is ordered correctly.
        rcnyr.phvrcn[yrd-1999,mthi]=subset(r2.max,mth==mthd & yr==yrd)$reconyear
        rcnyr.recon[yrd-1999,mthi]=subset(r2.max.recon,mth==mthd & yr==yrd)$reconyear
        mthi=mthi+1
    }
}

bestphvrcn=data.frame(rcnyr.phvrcn)
rownames(bestphvrcn) <- 2000:2011
colnames(bestphvrcn) <- c('Mar','Apr','May')
write.table(x=bestphvrcn,file='../swe.validation/best_recon.phvrcn.drop1.txt',sep=' ')

bestrecon=data.frame(rcnyr.recon)
rownames(bestrecon) <- 2000:2011
colnames(bestrecon) <- c('Mar','Apr','May')
write.table(x=bestphvrcn,file='../swe.validation/best_recon.recon.drop1.txt',sep=' ')




## save.image('cvdrop1.xrecon.regression.gauss.snotelunscaled.skill.image.RData')
load('cvdrorp1.xrecon.regression.gauss.snotelunscaled.skill.image.RData')






##############-------------- MAPS ---------------------------

require(lubridate)
states=map_data('state')
#### Map of Avg SWE Difference in drop10 prediction
# and percent swe difference
swediffmap=data.frame()
for(yri in 2000:2011){
    for(mthi in 3:5){
        mths=sprintf('0%d',mthi)
        mths=strftime(ymd(paste(yri,mths,'01',sep=''),tz='MST'),'%b')
        for(ryri in 2000:2011){
            df=data.frame(mth=mths,
                yr=yri,
                reconyear=ryri,
                swediff.phv=swediff.phv.avg[yri-1999,mthi-2,ryri-1999,],
                swediff.phvrcn=swediff.phvrcn.avg[yri-1999,mthi-2,ryri-1999,],
                swediff.recon=swediff.reconmodel[yri-1999,mthi-2,ryri-1999,],
                swediffpct.phv=swediffpct.phv.avg[yri-1999,mthi-2,ryri-1999,],
                swediffpct.phvrcn=swediffpct.phvrcn.avg[yri-1999,mthi-2,ryri-1999,],
                swediffpct.recon=swediffpct.reconmodel[yri-1999,mthi-2,ryri-1999,],
                lat=snotel.locs$lat,
                long=snotel.locs$long)
            swediffmap=rbind(swediffmap,df)
        }
    }
}
#
#

## Get shaded relief for background
relief=raster('shaded.relief.grey/GRAY_HR_SR.tif')
e=extent(x=c(-112.5,-104.25),y=c(32,44))
ucorelief=crop(relief,e)
ucorelief.agg=aggregate(ucorelief,2)
ucorelief.df=data.frame(rasterToPoints(ucorelief.agg))
names(ucorelief.df)=c('long','lat','ter')
#
## Get info about best rcnyr for phvrcn model from drop1 r2
bestmodel=read.table('~/Documents/R/swe.validation/best_recon.phvrcn.drop1.txt')




##--------------- PERCENT DIFFERENCES IN SWE INSTEAD OF SWE DIFF (METERS)
dev.set(2)
## Map regression phv % swe differences with observed facet by year
swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
swediffmapplot$swediffpct.phv=swediffmapplot$swediffpct.phv*100
swediffmapplot$cuts=cut(swediffmapplot$swediffpct.phv,breaks=seq(-200,200,length.out=17))
#
swediffmapplot$extreme=NA
swediffmapplot$extreme[swediffmapplot$swediffpct.phv > 200] ='> 200%'
swediffmapplot$extreme[swediffmapplot$swediffpct.phv < -200] = '< -200%'
#
extremesz=2
szmin=2
szstep=1
szmax=szmin+szstep*(length(levels(swediffmapplot$cuts))-1)/2
#
gphv=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,swediffpct.phv <= 200 | swediffpct.phv >= -200),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
    geom_point(data=subset(swediffmapplot,swediffpct.phv > 200 | swediffpct.phv < -200),aes(x=long,y=lat,shape=extreme),size=3,fill='black',alpha=0.75)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
    scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F)+
    scale_shape_manual(values=c(25,24))+
    guides(shape=guide_legend('Extremes'),
           alpha=F,
           size=guide_legend('SWE Differences (%)'),
           colour=guide_legend('SWE Differences (%)'))+
    coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey80'),legend.background=element_rect(fill='grey80'))+
    ggtitle('April 1 Differences with Observed SNOTEL SWE\nRegression with only Physiographics')+facet_wrap(~yr)
show(gphv)

## ggsave(plot=gphv,file='plots/swediffpct.phv.map.pdf',height=12,width=12)

dev.set(3)
## Map regression with rcn % swe differences with observed facet by year
swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr & yr==2002)
swediffmapplot$swediffpct.phvrcn=swediffmapplot$swediffpct.phvrcn*100
swediffmapplot$cuts=cut(swediffmapplot$swediffpct.phvrcn,breaks=seq(-200,200,length.out=17))
#
swediffmapplot$extreme=NA
swediffmapplot$extreme[swediffmapplot$swediffpct.phvrcn > 200] ='> 200%'
swediffmapplot$extreme[swediffmapplot$swediffpct.phvrcn < -200] = '< -200%'
#
extremesz=2
szmin=2
szstep=1
szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts))-2)/2)
#
gphvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,swediffpct.phvrcn <= 200 | swediffpct.phvrcn >= -200),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
    geom_point(data=subset(swediffmapplot,swediffpct.phvrcn > 200 | swediffpct.phvrcn < -200),aes(x=long,y=lat,shape=extreme),size=3,fill='black',alpha=0.75)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
    scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F)+
    scale_shape_manual(values=c(25,24))+
    guides(shape=guide_legend('Extremes'),
           alpha=F,
           size=guide_legend('SWE Differences (%)'),
           colour=guide_legend('SWE Differences (%)'))+
    coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey80'),legend.background=element_rect(fill='grey80'))+
    ggtitle('April 1 Differences with Observed SNOTEL SWE\nRegression with Physiographics and RCN')+facet_wrap(~yr)
show(gphvrcn)

## ggsave(plot=gphvrcn,file='plots/swediffpct.phvrcn.2001.map.pdf',height=12,width=10)




dev.set(4)
## Map reconstruction % swe differences with observed. facet by year
swediffmapplot=subset(swediffmap,mth=='Apr')
swediffmapplot$swediffpct.recon=swediffmapplot$swediffpct.recon*100
swediffmapplot$cuts=cut(swediffmapplot$swediffpct.recon,breaks=seq(-200,200,length.out=17))
#
swediffmapplot$extreme=NA
swediffmapplot$extreme[swediffmapplot$swediffpct.recon > 200] ='> 200%'
swediffmapplot$extreme[swediffmapplot$swediffpct.recon < -200] = '< -200%'
#
extremesz=2
szmin=2
szstep=1
szmax=szmin+szstep*(length(levels(swediffmapplot$cuts))-1)/2
#
grecon=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,swediffpct.recon <= 200 | swediffpct.recon >= -200),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
    geom_point(data=subset(swediffmapplot,swediffpct.recon > 200 | swediffpct.recon < -200),aes(x=long,y=lat,shape=extreme),size=3,fill='black',alpha=0.75)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
    scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F)+
    scale_shape_manual(values=c(25,24),drop=F)+
    guides(shape=guide_legend('Extremes'),
           alpha=F,
           size=guide_legend('SWE Differences (%)'),
           colour=guide_legend('SWE Differences (%)'))+
    coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey80'),
          legend.background=element_rect(fill='grey80'))+
    ggtitle('April 1 Differences with Observed SNOTEL SWE\nReconstruction')+
    facet_wrap(~yr)
show(grecon)

## ggsave(plot=grecon,file='plots/swediffpct.recon.map.pdf',height=12,width=12)


facet1_names=list(
    'swephvpct.phv'='Regression w/o RCN',
    'swephvpct.phvrcn'='Regression w/ RCN',
    'swephvpct.recon'='Reconstruction')
plot_labeller <- function(variable,value){
    if(variable=='yr'){
        return(value)
    } else {
        return(facet1_names[value])
    }
 }

#### Combined mapping - rows = modeltype, columns 4 years.
numf=11
outlier=250
swediffmapplot=melt(swediffmap,id=c('mth','yr','reconyear','lat','long','swediff.phv','swediff.phvrcn','swediff.recon'))
#
mnth='Apr'
mthcol=which(names(bestmodel)==mnth)
swediffmapplot=subset(swediffmapplot,(yr==2002 & reconyear==bestmodel[2002-1999,mthcol]) | (yr==2000 & reconyear==bestmodel[2000-1999,mthcol]) | (yr==2011 & reconyear==bestmodel[2011-1999,mthcol]))
swediffmapplot=subset(swediffmapplot,as.character(mth)==mnth)
swediffmapplot$value=swediffmapplot$value*100
swediffmapplot$value[swediffmapplot$value > outlier]=outlier-1
swediffmapplot$value[swediffmapplot$value < -outlier]=-outlier+1
swediffmapplot$cuts=cut(swediffmapplot$value,breaks=seq(-outlier,outlier,length.out=numf),right=F,include.lowest=T)
#
legendlabels=levels(swediffmapplot$cuts)
outlabel=outlier*2/(numf-1)*((numf-1)/2-1)
legendlabels[1]=paste('< ',-outlabel,sep='')
legendlabels[length(levels(swediffmapplot$cuts))]=paste('>= ',outlabel,sep='')
#
szmin=2
szstep=.5
szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts)))/2)
#
gphvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,value < 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
    geom_point(data=subset(swediffmapplot,value >= 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
    scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F,labels=legendlabels)+
    scale_shape_manual(values=c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2)),drop=F,labels=legendlabels)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F,labels=legendlabels)+
    labs(x='Longitude',y='Latitude')+
    guides(alpha=F,
           size=guide_legend('SWE Differences (%)'),
           colour=guide_legend('SWE Differences (%)'),
           shape=guide_legend('SWE Differences (%)'))+
    coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey80'),
          legend.key.size=unit(1.5,'lines'),
          legend.text=element_text(size=14),
          legend.title=element_text(size=16),
          ## legend.justification='right',
          legend.background=element_rect(fill='grey80'),
          legend.position='right',
          axis.text=element_text(size=14),
          plot.title=element_text(size=18),
          axis.title=element_text(size=16),
          strip.text=element_text(size=14,face='bold'))+
    ggtitle('April 1 % Differences with Observed SNOTEL SWE\n')+
    facet_grid(variable~yr,labeller=plot_labeller)
show(gphvrcn)
ggsave(plot=gphvrcn,'plots/allmodels.pctdiff.2002.2000.2011.png',dpi=300,width=11,height=10)


#********************************
## Map differences between %differences of the regression models facet by year
mnth='Apr'
mthcol=which(names(bestmodel)==mnth)
swediffmapplot=subset(swediffmap,(yr==2002 & reconyear==bestmodel[2002-1999,mthcol]) | (yr==2000 & reconyear==bestmodel[2000-1999,mthcol]) | (yr==2011 & reconyear==bestmodel[2011-1999,mthcol]))
swediffmapplot=subset(swediffmapplot,as.character(mth)==mnth)
     #
numf=11
#
swediffmapplot$mdldiffpct.phv.phvrcn=(abs(swediffmapplot$swediffpct.phv)-abs(swediffmapplot$swediffpct.phvrcn))*100
swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn > outlier]=outlier-1
swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn < -outlier]=-outlier+1
swediffmapplot$cuts=cut(swediffmapplot$mdldiffpct.phv.phvrcn,breaks=seq(-outlier,outlier,length.out=numf),right=F, include.lowest=T)
#
legendlabels=levels(swediffmapplot$cuts)
outlabel=outlier*2/(numf-1)*(numf-1)/2-outlier*2/(numf-1)
legendlabels[1]=paste('< ',-outlabel,sep='')
legendlabels[length(levels(swediffmapplot$cuts))]=paste('>= ',outlabel,sep='')
#
extremesz=4
szmin=3
szstep=1.5
szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts)))/2) 
#
gdiff.phv.phvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn >= 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
    geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn < 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
    scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F, labels=legendlabels)+
    scale_shape_manual(values=c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2)),drop=F,labels=legendlabels)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F,
                      labels=legendlabels)+
    labs(x='Longitude', y='Latitude')+
    guides(alpha=F,
           size=guide_legend('Error Differences (%)'),
           colour=guide_legend('Error Differences (%)'),
           shape=guide_legend('Error Differences (%)'))+
    coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey80'),
          legend.key.size=unit(1.5,'lines'),
          legend.background=element_rect(fill='grey80'),
          legend.text=element_text(size=14),
          legend.title=element_text(size=16),
          axis.text=element_text(size=14),
          plot.title=element_text(size=18),
          axis.title=element_text(size=16),
          strip.text=element_text(size=14,face='bold'))+
    facet_wrap(~yr)
    ## ggtitle(paste('April 1 Differences between % Error of Regression w/o RCN and Regression w/ RCN\nRelative to Observed SNOTEL SWE\nBlue triangles indicate a smaller error with Regression w/ RCN',sep=''))+
   #
show(gdiff.phv.phvrcn)
 ## ggsave(plot=gdiff.phv.phvrcn,file='plots/mdldiff.phv.phvrcn.2000.2002.2011.png',dpi=300,width=14, height=6.5)
## }




## Map differences between %differences of phvrcn and recon models. facet by year
outliers=100
swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
swediffmapplot$mdldiffpct.recon.phvrcn=(abs(swediffmapplot$swediffpct.recon)-abs(swediffmapplot$swediffpct.phvrcn))*100
swediffmapplot$cuts=cut(swediffmapplot$mdldiffpct.recon.phvrcn,breaks=seq(-outliers,outliers,length.out=21))
#
swediffmapplot$extreme=NA
swediffmapplot$extreme[swediffmapplot$mdldiffpct.recon.phvrcn > outliers] = paste('> ',outliers,'%',sep='')
swediffmapplot$extreme[swediffmapplot$mdldiffpct.recon.phvrcn < -outliers] = paste('< ',-outliers,'%',sep='')
#
extremesz=1
szmin=1
szstep=.5
szmax=szmin+szstep*(length(levels(swediffmapplot$cuts))-1)/2
#
gdiff.recon.phvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,mdldiffpct.recon.phvrcn <= outliers | mdldiffpct.recon.phvrcn >= -outliers),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
    geom_point(data=subset(swediffmapplot,mdldiffpct.recon.phvrcn > outliers | mdldiffpct.recon.phvrcn < -outliers),aes(x=long,y=lat,shape=extreme),size=extremesz,fill='black',alpha=0.75)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
    scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F)+
    scale_shape_manual(values=c(25,24),drop=F)+
    guides(shape=guide_legend('Extremes'),
           alpha=F,
           size=guide_legend('Error Differences (%)'),
           colour=guide_legend('Error Differences (%)'))+
    coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey80'),legend.background=element_rect(fill='grey80'))+
    ggtitle('April 1 Differences between Reconstruction and Regression with PHV/RCN\nRelative to Observed SNOTEL SWE\nBlue points indicate a smaller error with Regression w/ RCN')+
    facet_wrap(~yr)
show(gdiff.recon.phvrcn)

## ggsave(plot=gdiff.recon.phvrcn,file='plots/modeldiffpct.recon.phvrcn.map.pdf',height=12,width=12)

dev.new()
## Map differences between %differences of phv and recon models. facet by year
outliers=100
swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
swediffmapplot$mdldiffpct.recon.phv=(abs(swediffmapplot$swediffpct.recon)-abs(swediffmapplot$swediffpct.phv))*100
swediffmapplot$cuts=cut(swediffmapplot$mdldiffpct.recon.phv,breaks=seq(-outliers,outliers,length.out=21))
#
swediffmapplot$extreme=NA
swediffmapplot$extreme[swediffmapplot$mdldiffpct.recon.phv > outliers] = paste('> ',outliers,'%',sep='')
swediffmapplot$extreme[swediffmapplot$mdldiffpct.recon.phv < -outliers] = paste('< ',-outliers,'%',sep='')
#
extremesz=1
szmin=1
szstep=.5
szmax=szmin+szstep*(length(levels(swediffmapplot$cuts))-1)/2
#
gdiff.recon.phv=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,mdldiffpct.recon.phv <= outliers | mdldiffpct.recon.phv >= -outliers),aes(x=long,y=lat,color=cuts,size=cuts),alpha=0.75)+
    geom_point(data=subset(swediffmapplot,mdldiffpct.recon.phv > outliers | mdldiffpct.recon.phv < -outliers),aes(x=long,y=lat,shape=extreme),size=extremesz,fill='black',alpha=0.75)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F)+
    scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F)+
    scale_shape_manual(values=c(25,24),drop=F)+
    guides(shape=guide_legend('Extremes'),
           alpha=F,
           size=guide_legend('Error Differences (%)'),
           colour=guide_legend('Error Differences (%)'))+
    coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey80'),legend.background=element_rect(fill='grey80'))+
    ggtitle('April 1 Differences between Reconstruction and Regression with only PHV\nRelative to Observed SNOTEL SWE\nBlue points indicate a smaller error with Regression with PHV')+
    facet_wrap(~yr)
show(gdiff.recon.phv)

## ggsave(plot=gdiff.recon.phvrcn,file='plots/modeldiffpct.recon.phvrcn.map.pdf',height=12,width=12)







str(swediffmap)

swehist=swediffmap
swehist[,c('swediffpct.phv','swediffpct.phvrcn','swediffpct.recon')]=swehist[,c('swediffpct.phv','swediffpct.phvrcn','swediffpct.recon')]*100
#
swehist=subset(swehist,yr==reconyear & mth=='Apr')
avgs=ddply(swehist[,c('mth','yr','swediffpct.phv','swediffpct.phvrcn','swediffpct.recon')],.(yr,mth),summarize,avg.phv=mean(swediffpct.phv,na.rm=T),avg.phvrcn=mean(swediffpct.phvrcn,na.rm=T),avg.recon=mean(swediffpct.recon,na.rm=T))
#
avgs=subset(avgs,mth='Apr')
#
avg.phv=NULL
avg.phvrcn=NULL
avg.recon=NULL
for(i in 1:12){
    avg.phv=rbind(avg.phv,rep(avgs$avg.phv[i],237))
    avg.phvrcn=rbind(avg.phvrcn,rep(avgs$avg.phvrcn[i],237))
    avg.recon=rbind(avg.recon,rep(avgs$avg.recon[i],237))
}

swehist=cbind(swehist,avg=c(as.vector(t(avg.phv)),as.vector(t(avg.phvrcn)),as.vector(t(avg.recon))))
swehistplot=melt(swehist,id=c('mth','yr','reconyear','lat','long','swediff.phv','swediff.phvrcn','swediff.recon','avg'))
swediffhistplot1=subset(swehistplot, yr<2004)

ggplot(data=swediffhistplot1)+geom_histogram(aes(x=value),binwidth=1000)
+
    facet_grid(variable~yr)



ggplot()+geom_histogram(data=swediffmap,aes(x=swediffpct.phvrcn),binwidth=0.05)+facet_wrap(~yr)

ggplot()+geom_histogram(data=swediffmap,aes(x=swediffpct.phv),binwidth=0.05)+facet_wrap(~yr)

ggplot()+geom_histogram(data=swediffmap,aes(x=swediffpct.recon),binwidth=0.05)+facet_wrap(~yr)






#Bias
bias.phvrcn=mean(swediffmap$swediff.phvrcn,na.rm=T)
bias.phv=mean(swediffmap$swediff.phv,na.rm=T)
bias.recon=mean(swediffmap$swediff.recon,na.rm=T)

#MAE
mae.phvrcn=sum(abs(swediffmap$swediff.phvrcn),na.rm=T)/nrow(swediffmap)#0.11
mae.phv=sum(abs(swediffmap$swediff.phv),na.rm=T)/nrow(swediffmap)#0.13
mae.recon=sum(abs(swediffmapplot$swediff.recon),na.rm=T)/nrow(swediffmapplot)#0.20

# phv - phvrcn
mdldiff.phv.phvrcn=abs(swediffmap$swediff.phv)-abs(swediffmap$swediff.phvrcn)
mean(mdldiff.phv.phvrcn,na.rm=T)
ggplot()+geom_histogram(aes(x=mdldiff.phv.phvrcn),binwidth=0.05)
#
# recon - phvrcn
mdldiff.recon.phvrcn=abs(swediffmap$swediff.recon)-abs(swediffmap$swediff.phvrcn)
mean(mdldiff.recon.phvrcn,na.rm=T)
ggplot()+geom_histogram(aes(x=mdldiff.recon.phvrcn),binwidth=0.05)
#
# recon - phv
mdldiff.recon.phv=abs(swediffmap$swediff.recon)-abs(swediffmap$swediff.phv)
mean(mdldiff.recon.phv,na.rm=T)
ggplot()+geom_histogram(aes(x=mdldiff.recon.phv),binwidth=0.05)
