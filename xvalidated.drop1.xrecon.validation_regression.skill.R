## setwd('~/Documents/R/snotel_regression')
library(sp)
library(spdep)
library(raster)
library(ggplot2)
library(scales)
library(MASS)
library(R.matlab)
library(lubridate)
library(reshape2)
library(plyr)
library(RColorBrewer)


load('modeldata.validationtimes.RData')

sitecoords=readMat('sitecoords.mat')
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

load('recon.v2.stack.validation.RData')
#REGRESSION AND STEPAIC
 snotelrecon=array(NA,dim=c(length(dates2keep),12,237))
 snotelrecon.sc=snotelrecon
## ## ## Read in all recon images and create dataframe of snotel pixel values
 ## mth='MAR'
## ## yr=2000
rv='v2.senlat.nldas.usace.no_canopy_windcorr'         
for(i in 1:length(dates2keep)){
     for(yr in 2000:2011){
     #yr=year(dates2keep[i])       
     mth=month(dates2keep[i],label=T,abbr=T)
     dy=sprintf('%02d',day(dates2keep[i]))
     snodis_geotiffs=paste('/Volumes/WSC/SWE_SNODIS/recon/',rv,'/UpperColoradoRiver/',yr,'/tif',sep='')
### Get current recon
     im.string=paste(snodis_geotiffs,'/',dy,toupper(mth),yr,'.tif',sep='')
     r=raster(im.string,crs='+proj=longlat +datum=NAD83')
     recon.v2.stack[[i]]=r
     names(r)='recon'
     snotelrecon[i,yr-1999,]=as.numeric(raster::extract(r,snotel.locs))
     snotelrecon.sc[i,yr-1999,]= scale(snotelrecon[i,yr-1999,])
     }
}
#
recon.v2.stack.scaled=recon.v2.stack
m=cellStats(recon.v2.stack,'mean')
stdev=cellStats(recon.v2.stack,'sd')
for(i in 1:nlayers(recon.v2.stack)){
     recon.v2.stack.scaled[[i]]=(recon.v2.stack[[i]]-m[i])/stdev[i]
 }
#
save(file='recon.v2.stack.validation.RData',list=c('recon.v2.stack','recon.v2.stack.scaled','snotelrecon','snotelrecon.sc','rv','dates2keep'))

swe.sc=as.data.frame(scale(snotel2model[,c('Latitude','Longitude','Elev','Aspect','Slope','RegionalSlope','RegionalAspect','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','NWd2ocean','SWd2ocean')]))
swe.sc=cbind(swe.sc,snotel=snotel2model$swe,mth=snotel2model$mth,yr=snotel2model$yr,date=snotel2model$date)
## print(str(swe.sc))



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
swediff.phvmodel=array(NA,dim=c(length(dates2keep),12,237))
swediff.phvrcnmodel=array(NA,dim=c(length(dates2keep),12,237))
swediffpct.phvmodel=array(NA,dim=c(length(dates2keep),12,237))
swediffpct.phvrcnmodel=array(NA,dim=c(length(dates2keep),12,237))
#
r2.phvmodel=array(NA,dim=c(length(dates2keep),12))
r2.phvrcnmodel=array(NA,dim=c(length(dates2keep),12))
r2.reconmodel=array(NA,dim=c(length(dates2keep),12))
#length(dates2keep)
for(j in 1:length(dates2keep)){
    yr=year(dates2keep[j])
    mth=month(dates2keep[j],label=T,abbr=T)
    dy=sprintf('%02d',day(dates2keep[j]))
    print(paste('iteration - ',j,': ',dy,mth,yr,sep=''))
          
          swe2model.all=swe.sc[swe.sc$date==dates2keep[j],]
          swe2model.all$snotel[swe2model.all$snotel==0]=runif(1,0,0.00001)
          for (ryr in seq(2000,2011)){
              print(paste('ryr = ',ryr))
              swe2model.all$recon=as.numeric(snotelrecon.sc[j,ryr-1999,])
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
                  SSphv = cor(pr.phv,yobs)^2#1 - SSE/SSNull
                  swediff.phv = pr.phv-yobs
                  swediffpct.phv = (pr.phv-yobs)/yobs
                  
                  SSphvrcn = cor(pr.phvrcn,yobs)^2 #1 - SSE/SSNull
                  swediff.phvrcn = pr.phvrcn-yobs
                  swediffpct.phvrcn = (pr.phvrcn-yobs)/yobs

                  yrecon=as.numeric(snotelrecon[j,yr-1999,])
                  SSrecon = cor(yrecon,yobs)^2
                  ## if(length(which(is.na(SSrecon)))!=0) print(paste('NA values - recon ',mth,yr,'- recon',ryr,sep=''))
         
              r2.phvmodel[j,ryr-1999]=SSphv
              r2.phvrcnmodel[j,ryr-1999]=SSphvrcn
              r2.reconmodel[j,ryr-1999]=SSrecon
              swediff.phvmodel[j,ryr-1999,]=swediff.phv
              swediff.phvrcnmodel[j,ryr-1999,]=swediff.phvrcn
              swediffpct.phvmodel[j,ryr-1999,]=swediffpct.phv
              swediffpct.phvrcnmodel[j,ryr-1999,]=swediffpct.phvrcn
              
	}
       
} 



## In a much smaller loop, find swediff for recon vs snotel
swediff.reconmodel=array(NA,dim=c(length(dates2keep),237))
swediffpct.reconmodel=array(NA,dim=c(length(dates2keep),237))
for(j in 1:length(dates2keep)){
	yr=year(dates2keep[j])
	mth=month(dates2keep[j],label=T,abbr=T)
	dy=sprintf('%02d',day(dates2keep[j]))
	print(paste('iteration - ',j,': ',dy,mth,yr,sep=''))
                                        #
	swe2model.all=swe.sc[swe.sc$date==dates2keep[j],]
	yobs=swe2model.all$snotel
	yrecon=as.numeric(snotelrecon[j,yr-1999,])
	print(paste('yobs: ',which(is.na(yobs)),sep=''))
	print(paste('yrecon: ',which(is.na(yrecon)),sep=''))
                                        #   
    swediff.reconmodel[j,]=yrecon-yobs
    #swediffpct.reconmodel[j,]=(yrecon-yobs)/yobs
	swediffpct.reconmodel[j,]=ifelse((yrecon-yobs)/yobs=='NaN',0,(yrecon-yobs)/yobs)
          #
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
ndf=data.frame()
for(ryr in 2000:2011){
	date=dates2keep
	site=sites4dates
	reconyear=ryr
	r2phv=r2phvmed[,ryr-1999]
	r2phvrcn=r2phvrcnmed[,ryr-1999]
	r2recon=r2reconmed[,ryr-1999]
	ndf=rbind(ndf,data.frame(date,site,reconyear,r2phv,r2phvrcn,r2recon))
 	}
r2df=ndf


#
##Find the mean r2 for phv and recon. Each recon year has cvitmax iterations for each timestep. the mean is used to represent r2 skill for a given timestep since phv and recon modesl don't change with the rcnyear used.
#for(i in 1:length(dates2keep)){
#     di=dates2keep[i]
#    r2df[r2df$date==di,'r2phv']=apply(r2df[di==r2df$date,c('r2phv','r2recon')],2,mean)[1]
#    r2df[r2df$date==di,'r2recon']=apply(r2df[r2df$date==di,c('r2phv','r2recon')],2,mean)[2]
#}


r2.max=ddply(.data=r2df,.(date,site),summarize,which.max(r2phvrcn),max(r2phvrcn))
names(r2.max)=c('date','site','modelnum','maxr2')
r2.max$reconyear=r2.max$modelnum+1999

write.table(x=r2.max,file='best_recon.phvrcn.validation.drop1.txt',sep=' ',row.names=F)

r2df.rt=subset(r2df,reconyear<year(date))
r2.rt.max=ddply(.data=r2df.rt,.(date,site),summarize,which.max(r2phvrcn),max(r2phvrcn))
names(r2.rt.max)=c('date','site','modelnum','maxr2')
r2.rt.max$reconyear=r2.rt.max$modelnum+1999

write.table(x=r2.rt.max,file='swe.validation/best_recon.rt.phvrcn.validation.drop1.txt',sep=' ',row.names=F)


#save.image('cvdrop1.validationtimesteps.xrecon.regression.gauss.snotelunscaled.skill.image.RData')
l#oad('snotel_regression/cvdrop1.validationtimesteps.xrecon.regression.gauss.snotelunscaled.skill.image.RData')





## #### Boxplot facet by model
## # note: each boxplot is a range of 12*cvitmax values. In calculation of r2, for pr.phv and recon only yobs changes (a different 10% each time). for pr.phvrcn, the recon year changes the prediction in addition to yobs sampling being different.
## model_names <- list(
##    'r2phv'='Regression w/o RCN',
##    'r2phvrcn'='Regression w/ RCN',
##    'r2recon'='Reconstruction')
## model_labeller <- function(variable,value) return(model_names[value])
## r2df.sub=subset(r2df,mth=='Apr')
## r2df.sub.m=melt(r2df.sub,id=c('yr','mth','reconyear'))
## r2sameyr=subset(r2df.sub.m,yr==reconyear)
## ggplot()+geom_boxplot(data=r2df.sub.m,aes(as.factor(yr),y=value),outlier.shape=3)+
##     geom_point(data=r2sameyr,aes(x=as.factor(yr),y=value),col='red')+
##     labs(x='Year Modeled',y=expression(r^2),title='April 1 Model Skill Ensemble')+
##     facet_grid(~variable,labeller=model_labeller)+
##     theme_minimal()


#### Boxplot facet by month
r2df.sameyr=subset(r2df,yr==reconyear)
gr2box=ggplot()+geom_boxplot(data=r2df,aes(as.factor(yr),r2phvrcn,col='Regression with RCN Ensemble'),outlier.shape=3,show_guide=T)+
    geom_boxplot(data=r2df,aes(x=as.factor(yr),y=r2phv,col='Regression w/o RCN'),show_guide=F)+
    geom_boxplot(data=r2df,aes(x=as.factor(yr),y=r2recon,col='Reconstruction'),show_guide=F)+
    geom_point(data=r2df.sameyr,aes(x=as.factor(yr),y=r2phvrcn,col='Regression w/ RCN from Same Year'),show_guide=T)+
    scale_color_manual(values=c('blue','red','green','black'))+
    facet_grid(~mth)+
    xlab('Year Modeled')+
    ylab(expression(r^2))+
    ylim(c(0,0.8))+
    theme_bw()+
    guides(color=guide_legend('Model',override.aes=list(shape=c(3,20,3,0))))+
    theme(legend.position = c(0.85,.18),legend.background=element_rect(color='black'))+
    labs(title=expression(paste('Ensemble of LOO Cross-Validated Skill Scores (',r^2,')',sep='')))
#
show(gr2box)
## ggsave(plot=gr2box,file='plots/r2.xval.boxplot.gauss.pdf',width=14,height=5)


### Boxplot realtime r2 by month
r2df.realtime=subset(r2df,yr>reconyear)
r2df.sameyr=subset(r2df,yr==reconyear)
gr2rtbox=ggplot()+geom_boxplot(data=r2df.realtime,aes(as.factor(yr),r2phvrcn,col='Regression with Real-time RCN Ensemble'),outlier.shape=3,show_guide=F)+
    geom_boxplot(data=r2df,aes(x=as.factor(yr),y=r2phv,col='Regression w/o RCN'),show_guide=F)+
    geom_point(data=r2df,aes(x=as.factor(yr),y=r2recon,col='Reconstruction'),shape=6)+
    geom_point(data=r2df.sameyr,aes(x=as.factor(yr),y=r2phvrcn,col='Regression w/ RCN from Same Year'),shape=6)+
    scale_color_manual(values=c('blue','red','green','black'))+
    facet_grid(~mth)+
    xlab('Year Modeled')+
    ylab(expression(r^2))+
    ylim(c(0,0.8))+
    theme_bw()+
    guides(color=guide_legend('Model',override.aes=list(shape=c(6,6,1,1))))+
    theme(legend.position = c(0.85,.18),legend.background=element_rect(color='black'))+
    labs(title=expression(paste('Real-time Ensemble of LOO Cross-Validated Skill Scores (',r^2,')',sep='')))
#
show(gr2rtbox)
## ggsave(plot=gr2rtbox,file='plots/r2.xval.realtime.boxplot.gauss.ps',width=14,height=5)


## Boxplot r2 mlr
r2df.modelm=melt(r2df,id=c('mth','yr','reconyear'))
ggplot(r2df.modelm)+geom_boxplot(aes(x=variable,y=value))+
    scale_x_discrete(labels=c('Regression\nw/o RCN','Regression\nw/ RCN','Recon'))+
    xlab('')+
    ylab(expression(r^2))+
    theme_minimal()

## ggsave('cv.mlrphv.ps')





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
                swediff.recon=swediff.reconmodel[yri-1999,mthi-2,],
                swediffpct.phv=swediffpct.phv.avg[yri-1999,mthi-2,ryri-1999,],
                swediffpct.phvrcn=swediffpct.phvrcn.avg[yri-1999,mthi-2,ryri-1999,],
                swediffpct.recon=swediffpct.reconmodel[yri-1999,mthi-2,],
                lat=snotel.locs$lat,
                long=snotel.locs$long)
            swediffmap=rbind(swediffmap,df)
        }
    }
}
#
#
#
## Get shaded relief for background
relief=raster('shaded.relief.grey/GRAY_HR_SR.tif')
e=extent(x=c(-112.5,-104.25),y=c(32,44))
ucorelief=crop(relief,e)
ucorelief.agg=aggregate(ucorelief,2)
ucorelief.df=data.frame(rasterToPoints(ucorelief.agg))
names(ucorelief.df)=c('long','lat','ter')

## ## Get cities for map
## cities=readOGR('ne_10m_populated_places','ne_10m_populated_places')
## cities=data.frame(subset(cities,LATITUDE > 32 & LATITUDE < 44 & LONGITUDE > -113 & LONGITUDE < -104))
## str(cities)
## ggplot()+geom_point(data=cities,aes(x=LONGITUDE,y=LATITUDE))+
##     geom_text(data=cities,aes(x = LONGITUDE, y = LATITUDE, label = NAMEALT))


## ## Map regression phv swe differences facet by year
## swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
## gphv=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,fill=ter))+
##     scale_fill_gradient(low='grey5',high='grey70')+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=swediffmapplot,aes(x=long,y=lat,color=swediff.phv,size=abs(swediff.phv)),alpha=0.8)+
##     scale_color_gradient2(low='red',mid='white',high='blue',limits=range(c(swediffmap$swediff.phv,swediffmap$swediff.phvrcn,swediffmap$swediff.recon),na.rm=T),na.value='green')+
##     guides(fill=F,color=guide_legend('SWE Differences (m)'),size=guide_legend('SWE Difference\nMagnitudes (m)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     theme_bw()+
##     ggtitle('April 1 SWE Differences with Observed\nRegression with only Physiographics')+
##     facet_wrap(~yr)
## show(gphv)
## ## ggsave(plot=gphv,file='plots/swediff.phv.map.pdf',height=12,width=12)





## dev.set(3)
## ## Map regression with rcn swe differences facet by year
## swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
## gphvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,fill=ter))+
##     scale_fill_gradient(low='grey5',high='grey70')+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=swediffmapplot,aes(x=long,y=lat,color=swediff.phvrcn,size=abs(swediff.phvrcn)),alpha=0.8)+
##     scale_color_gradient2(low='red',mid='white',high='blue',limits=range(c(swediffmap$swediff.phv,swediffmap$swediff.phvrcn,swediffmap$swediff.recon),na.rm=T),na.value='green')+
##     guides(fill=F,color=guide_legend('SWE Differences (m)'),size=guide_legend('SWE Difference\nMagnitudes(m)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     facet_wrap(~yr)+
##     theme_bw()+
##     ggtitle('April 1 Differences with Observed SWE\nRegression with Physiographics and Reconstructed SWE')
## show(gphvrcn)
## ggsave(plot=gphvrcn,file='plots/swediff.phvrcn.map.pdf',height=12,width=12)

## dev.set(4)
## ## Map differences between regression models facet by year
## swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
## swediffmapplot$mdldiff.phv.phvrcn=abs(swediffmapplot$swediff.phv)-abs(swediffmapplot$swediff.phvrcn)
## gdiff=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,fill=ter))+
##     scale_fill_gradient(low='grey5',high='grey70')+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=swediffmapplot,aes(x=long,y=lat,color=mdldiff,size=abs(mdldiff)),alpha=0.8)+
##     scale_color_gradient2(low='red',mid='white',high='blue',limits=range(c(swediffmapplot$mdldiff),na.rm=T),na.value='green')+
##     guides(fill=F,color=guide_legend('SWE Differences (m)'),size=guide_legend('SWE Difference\nMagnitudes (m)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     facet_wrap(~yr)+
##     theme_bw()+
##     ggtitle('April 1 Differences between Regression Models\nRelative to Observed SNOTEL SWE\nPositive (blue) indicates Regression w/o RCN had a greater difference with observed SWE')
## ## show(gdiff)
## ggsave(plot=gdiff,file='plots/modeldiffpct.phv.phvrcn.map.pdf',height=12,width=12)


## dev.set(5)
## ## Map reconstruction swe differences facet by year
## swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
## grecon=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,fill=ter))+
##     scale_fill_gradient(low='grey5',high='grey70')+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=swediffmapplot,aes(x=long,y=lat,color=swediff.recon,size=abs(swediff.recon)),alpha=0.8)+
##     scale_color_gradient2(low='red',mid='white',high='blue',limits=range(c(swediffmap$swediff.phv,swediffmap$swediff.phvrcn,swediffmap$swediff.recon),na.rm=T),na.value='green')+
##     guides(fill=F,color=guide_legend('SWE Differences (m)'),size=guide_legend('SWE Difference\nMagnitudes(m)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     facet_wrap(~yr)+
##     theme_bw()+
##     ggtitle('April 1 Differences with Observed SWE\nReconstructed SWE')
## ## show(grecon)
## ggsave(plot=grecon,file='plots/swediff.recon.map.pdf',height=12,width=12)


## dev.set(6)
## ## Map differences between phvrcn and recon models facet by year
## swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==yr)
## swediffmapplot$mdldiff.recon.phvrcn=abs(swediffmapplot$swediff.recon)-abs(swediffmapplot$swediff.phvrcn)
## gdiffrecon=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,fill=ter))+
##     scale_fill_gradient(low='grey5',high='grey70')+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=swediffmapplot,aes(x=long,y=lat,color=recon.phvrcndiff,size=abs(recon.phvrcndiff)),alpha=0.8)+
##     scale_color_gradient2(low='red',mid='white',high='blue',limits=range(c(swediffmapplot$mdldiff,swediffmapplot$recon.phvrcndiff),na.rm=T),na.value='green')+
##     guides(fill=F,color=guide_legend('SWE Differences (m)'),size=guide_legend('SWE Difference\nMagnitudes (m)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     facet_wrap(~yr)+
##     theme_bw()+
##     ggtitle('April 1 Differences between Regression with RCN and Recon\nRelative to Observed SNOTEL SWE\nPositive (blue) indicates Recon had a greater difference with observed SWE')
## ## show(gdiffrecon)
## ggsave(plot=gdiffrecon,file='plots/modeldiffpct.recon.phvrcn.map.pdf',height=12,width=12)



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
    theme(legend.key=element_rect(fill='grey80'),legend.background=element_rect(fill='grey80'))+
    ggtitle('April 1 Differences with Observed SNOTEL SWE\nReconstruction')+
    facet_wrap(~yr)
show(grecon)

## ggsave(plot=grecon,file='plots/swediffpct.recon.map.pdf',height=12,width=12)


facet1_names=list(
    'swephvpct.phv'='Regression-PHV',
    'swephvpct.phvrcn'='Regression-PHV/RCN',
    'swephvpct.recon'='Reconstruction')
plot_labeller <- function(variable,value){
    if(variable=='yr'){
        return(value)
    } else {
        return(facet1_names[value])
    }
 }

#### Combined mapping - rows = modeltype, columns 4 years.
outlier=250
swediffmapplot=melt(swediffmap,id=c('mth','yr','reconyear','lat','long','swediff.phv','swediff.phvrcn','swediff.recon'))
#
swediffmapplot=subset(swediffmapplot,mth=='Apr' & (reconyear==rcnyr[9] &  yr==2008) | (reconyear==rcnyr[10] &  yr==2009) | (reconyear==rcnyr[11] &  yr==2010) | (reconyear==rcnyr[12] &  yr==2011))
swediffmapplot$value=swediffmapplot$value*100
swediffmapplot$value[swediffmapplot$value > outliers]=outliers
swediffmapplot$value[swediffmapplot$value < -outliers]=-outliers+1
swediffmapplot$cuts=cut(swediffmapplot$value,breaks=seq(-outlier,outlier,length.out=11))
#
legendlabels=levels(swediffmapplot$cuts)
legendlabels[1]='<= -200'
legendlabels[length(levels(swediffmapplot$cuts))]='> 200'
#
szmin=2
szstep=.5
szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts)))/2)
#
gphvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,value < 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
    geom_point(data=subset(swediffmapplot,value > 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
    scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F,labels=legendlabels)+
    scale_shape_manual(values=c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2)),drop=F,labels=legendlabels)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F,labels=legendlabels)+
    guides(alpha=F,
           size=guide_legend('SWE Differences (%)'),
           colour=guide_legend('SWE Differences (%)'),
           shape=guide_legend('SWE Differences (%)'))+
    coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey80'),
          legend.background=element_rect(fill='grey80'),
          strip.text=element_text(size=12,face='bold'))+
    ggtitle('April 1 Differences with Observed SNOTEL SWE\n')+
    facet_grid(variable~yr,labeller=plot_labeller)
show(gphvrcn)
## ggsave(plot=gphvrcn,'plots/allmodels.pctdiff.2008_2012.pdf',width=15,height=15)


#********************************
## Map differences between %differences of the regression models facet by year
 for(yri in seq(2000,2011,1)){

     if(!exists(rcnyr)){
         rcnyr=read.table('best_recon.phvrcn.drop1.txt')
     }
     rcnyr=rcnyr[yri-1999,2]
     
     outliers=250
swediffmapplot=subset(swediffmap,mth=='Apr' & reconyear==rcnyr & yr==yri)
swediffmapplot$mdldiffpct.phv.phvrcn=(abs(swediffmapplot$swediffpct.phv)-abs(swediffmapplot$swediffpct.phvrcn))*100
swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn > outliers]=outliers
swediffmapplot$mdldiffpct.phv.phvrcn[swediffmapplot$mdldiffpct.phv.phvrcn < -outliers]=-outliers+1
swediffmapplot$cuts=cut(swediffmapplot$mdldiffpct.phv.phvrcn,breaks=seq(-outliers,outliers,length.out=11))
#
legendlabels=levels(swediffmapplot$cuts)
legendlabels[1]='<= -200'
legendlabels[length(levels(swediffmapplot$cuts))]='> 200'
#
extremesz=4
szmin=3
szstep=1.5
szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts)))/2) 
#
gdiff.phv.phvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(1,0.5))+
    geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
    geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn > 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
    geom_point(data=subset(swediffmapplot,mdldiffpct.phv.phvrcn < 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
    scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F, labels=legendlabels)+
    scale_shape_manual(values=c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2)),drop=F,labels=legendlabels)+
    scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F,
                      labels=legendlabels)+
    guides(alpha=F,
           size=guide_legend('Error Differences (%)'),
           colour=guide_legend('Error Differences (%)'),
           shape=guide_legend('Error Differences (%)'))+
    coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
    theme_bw()+
    theme(legend.key=element_rect(fill='grey80'),
          legend.background=element_rect(fill='grey80'))+
    ggtitle(paste('April 1, ',yri,' Differences between Regression-PHV and Regression-PHV/RCN\nRelative to Observed SNOTEL SWE\nBlue triangles indicate a smaller error with Regression-PHV/RCN',sep=''))+
    facet_wrap(~yr)
show(gdiff.phv.phvrcn)
ggsave(plot=gdiff.phv.phvrcn,file=paste('plots/mdldiff.phv.phvrcn.',yri,'.pdf',sep=''),width=10, height=12)
}




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
    ggtitle('April 1 Differences between Reconstruction and Regression with PHV/RCN\nRelative to Observed SNOTEL SWE\nBlue points indicate a smaller error with Regression with PHV/RCN')+
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
