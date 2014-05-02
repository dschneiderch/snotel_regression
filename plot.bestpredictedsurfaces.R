##############################
setwd('~/Documents/R/snotel_regression')
## library(RColorBrewer)
## library(sp)
## library(spdep)
## library(raster)
## library(ggplot2)
## library(scales)
## library(MASS)
## library(rgdal)
## require(plyr)
## require(reshape2)
## require(grDevices)
 require(R.matlab)
## require(ncdf4)

## load('snoteldata.RData')
## load('ucophv.RData')
#load('latest.model/mdl.stepaic.gauss.snotelscaled.xrecon.Ftest.RData')
#load('latest.model/prediction.surfaces.gauss.snotelscaled.xrecon.Ftest.RData')
## load('/Volumes/snowFTP/pub/fromDominik/forBalaji/regress.data.RData')

## load('~/Documents/R/snotel_regression/latest.model/mdl.stepaic.gauss.snotelscaled.xrecon.Ftest.RData')
# load('~/Documents/R/snotel_regression/latest.model/prediction.surfaces.gauss.snotelunscaled.xrecon.Ftest.RData')

## ##*********************
## # # Map domain and snotel stations
## # #*** For Upper CO River Basin domain
## upleft=c(43.75,-112.25)
## lowright=c(33,-104.125)
## lat.low=33
## lat.high=43.75
## lim=c(lon.low, lon.high, lat.low, lat.high)
## ##

## ## ------ Load Snotel coordinates
## sitecoords=readMat('../import.process.matlabdata/snotel_attributes/sitecoords.mat')
## slt=sitecoords[[1]][,1]
## sln=sitecoords[[1]][,2]
## snotel.locs=data.frame(lat=slt,long=sln)
## coordinates(snotel.locs)<- ~long+lat
## proj4string(snotel.locs) <- '+proj=longlat +datum=WGS84'
## ##

## ## --- Load snotel data
## load('~/Documents/R/import.process.matlabdata/snoteldata.RData')
## swe=cbind(swe,phvrec)
## swe$mth=factor(swe$mth,c('Mar','Apr','May'))
## #scale the x variables.
## swe.sc=as.data.frame(scale(swe[,c('Latitude','Longitude','Elev','Aspect','Slope','RegionalSlope','RegionalAspect','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','NWd2ocean','SWd2ocean')]))
## swe.sc=cbind(swe.sc,snotel=swe$snotel,mth=swe$mth,yr=swe$yr,date=swe$date)
## # scale the y variable inside loop


## col.pal=brewer.pal(7, 'RdYlBu') # or use 'RdBu'

## print('starting plotting of surfaces....')




## ## ###--------------


## states <- map_data("state")
## ## ggplot(states)+geom_polygon(data=states,aes(x=long,y=lat,group=group),fill='white',alpha=0.5,color='grey')+coord_cartesian(xlim=c(-112.25,-104.125),ylim=c(33,43.75))+theme_bw()

## source('get_bestmodel.R')
## snotel.locs=data.frame(snotel.locs) 
## statsdf.pr=NULL
## statsdf.r=NULL
## g=list()
## jm=1
## ryr=1
## yrcnt=1
## for( jm in seq(1:nrow(r2.max))){
##     print(paste('Beginning iteration -',jm,'of',nrow(r2.max)))

##     if(jm%%3==1) mth='Mar'
##     if(jm%%3==2) mth='Apr'
##     if(jm%%3==0) mth='May'     

   
## ### ---   Get best regression surface and plot model
##     rcnyr=get_bestmodel(mth,yrcnt+1999)
##     r=recon.v2.stack[[jm]]
##     prswe=pr.surface.stack[[jm]][[rcnyr-1999]]
##     prswe[r==0]=0
##     prswe[prswe<0]=0
##     df=data.frame(rasterToPoints(prswe))#pr.surface.stack[[jm]][[ryr]]))
##     colnames(df)=c('long','lat','swe')

##     r2k=seq(1,nrow(df),1)
##     df.sub=df[r2k,]
##     df.sub$swe[df.sub$swe==0]=NA

##     statsdf.pr$max.m[jm]=max(df$swe,na.rm=T)
##     statsdf.pr$vol.km3[jm]=sum(df$swe,na.rm=T)*500*500/1000^2
##     statsdf.pr$avg.m[jm]=mean(df$swe,na.rm=T)
##     statsdf.pr$med.m[jm]=median(df$swe,na.rm=T)
##     statsdf.pr$mth[jm]=mth
##     statsdf.pr$yr[jm]=yrcnt+1999
##     statsdf.pr$model[jm]='regression'
##     ## df.sub$breaks=cut_interval(x=df.sub$swe,length=0.25)

##     brks=c(0.001,0.5,1,1.5,2,2.5)   
##     gpr=ggplot()+ geom_raster(data=df.sub,aes(x=long,y=lat,fill=swe))+
##         geom_path(data=states,aes(x=long,y=lat,group=group),col='black')+
##             scale_fill_gradientn(colours=brewer.pal(name='RdYlBu',n=11),na.value='grey90',breaks=brks,labels=brks)+
##                 geom_point(data=snotel.locs,aes(x=long,y=lat),size=1)+
##                     theme_minimal()+
##                         coord_cartesian(xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
##                             ggtitle(paste('best regression (gaussian) surface - 01',mth,yrcnt+1999,sep=''))
##   ## show(gpr)


    
##     ## ## longitude=ncdim_def('Longitude','degrees',coordinates(prswe)[,1])
##     ## ## latitude=ncdim_def('Latitude','degrees',coordinates(prswe)[,2])
##     ## ## swevar=ncvar_def('swe',units='meters',dim=list(longitude,latitude),missval=-9999)
##     ## ncfn='~/Documents/R/snotel_regression/plots/netcdf.surfaces/tempsurface.nc'
##     ## ## ncnew <- nc_create(ncfn,vars=swevar)
##     ## ## ncvar_put(ncnew,swevar,values(prswe))
##     ## ## nc_close(ncnew)

##      ncfn='~/Documents/R/snotel_regression/plots/netcdf.surfaces/tempsurface.nc'
##      writeRaster(prswe,ncfn,format='CDF',overwrite=T)
##      ncfnnew=paste('~/Documents/R/snotel_regression/plots/netcdf.surfaces/regression.surface.01',mth,yrcnt+1999,'.nc',sep='')
##      syscmd=paste('gdal_translate -of netcdf -a_srs NAD83 -a_nodata -9999',ncfn,ncfnnew,sep=' ')
##      system(syscmd)
##      system('rm ~/Documents/R/snotel_regression/plots/netcdf.surfaces/tempsurface.nc')


    
## ## ###----- Get recon and plot for comparison
## ##     df=data.frame(rasterToPoints(r))
## ##     names(df)=c('long','lat','swe')

## ##     r2k=seq(1,nrow(df),1)
## ##     df.sub=df[r2k,]
## ##     df.sub$swe[df.sub$swe==0]=NA

## ##     statsdf.r$max.m[jm]=max(df$swe,na.rm=T)
## ##     statsdf.r$vol.km3[jm]=sum(df$swe,na.rm=T)*500*500/1000^2
## ##     statsdf.r$avg.m[jm]=mean(df$swe,na.rm=T)
## ##     statsdf.r$med.m[jm]=median(df$swe,na.rm=T)
## ##     statsdf.r$mth[jm]=mth
## ##     statsdf.r$yr[jm]=yrcnt+1999
## ##     statsdf.r$model[jm]='recon'
## ##     #df.sub$breaks=cut_interval(x=df.sub$swe,length=0.25)

## ##     ## grcn=ggplot()+geom_raster(data=df.sub,aes(x=long,y=lat,fill=swe))+
## ##     ##     scale_fill_gradientn(colours=brewer.pal(name='RdYlBu',n=11),na.value='grey90',breaks=brks,limits=c(0.001,3),labels=brks)+
## ##     ##         geom_path(data=states,aes(x=long,y=lat,group=group),col='black')+
## ##     ##             theme_minimal()+
## ##     ##                 coord_cartesian(xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
## ##     ##                     ggtitle(paste('recon v2 surface - 01',mth,yrcnt+1999,sep=''))
## ##     ## ## show(grcn)

##     fn1=paste('plots/gaussian.surfaces/regression.surface.01',mth,yrcnt+1999,'.pdf',sep='')
## ##     ## fn2=paste('plots/recon.surface.01',mth,yrcnt+1999,'.pdf',sep='')
##     ggsave(plot=gpr,filename=fn1,width=6,height=7)    
## ##     ## ggsave(plot=grcn,filename=fn2,width=6,height=7)

    
##     if(jm%%3==0) yrcnt=yrcnt+1
## }

## statsdf=rbind(data.frame(statsdf.pr),data.frame(statsdf.r))
## write.table(statsdf,'surface.gaussian.snotelunscaled.stats.csv',sep=',',row.names = F)





#### ---------------- Combined mapping - rows = modeltype, columns 4 years.
library(lubridate)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(plyr)
library(reshape2)
library(raster)
library(rasterVis)
library(sp)

## Get shaded relief for background
relief=raster('shaded.relief.grey/GRAY_HR_SR.tif')
e=extent(x=c(-112.5,-104),y=c(32,44))
ucorelief=crop(relief,e)
ucorelief.agg=aggregate(ucorelief,2)
ucorelief.df=as.data.frame(ucorelief.agg,xy=T)
names(ucorelief.df)=c('long','lat','ter')

## Get state lines from maps package
states <- map_data("state")

## Get SNOTEL locations
sitecoords=readMat('../import.process.matlabdata/snotel_attributes/sitecoords.mat')
slt=sitecoords[[1]][,1]
sln=sitecoords[[1]][,2]
snotel.locs=data.frame(lat=slt,long=sln)
coordinates(snotel.locs)<- ~long+lat
proj4string(snotel.locs) <- '+proj=longlat +datum=WGS84'
#

## Assign correct projection to reconstruction
#projection(recon.v2.stack) <- '+proj=longlat +NAD83'

## jm=2 #april 1, 2000
## jm=8 #april 1 2002
## jm=35 #april 1 2011

## Get info for best model (recon year that gives best r2 in drop1 for phvrcn model)
bestmodel=read.table('~/Documents/R/swe.validation/best_recon.phvrcn.drop1.txt',sep=' ',header=T)
#

rast2df <- function(prdate){
yr=year(prdate)
mth=month(prdate)
mtha=month(prdate,label=T,abbr=T)
dy=day(prdate)
jm=((yr-1999)-1)*3+(mth-2)
rcnyr=bestmodel[rownames(bestmodel) %in% as.character(yr),names(bestmodel) %in% mtha]
ryr=rcnyr-1999
#
prswe.recon=recon.v2.stack[[jm]]
prswe.recon=sampleRegular(prswe.recon,size=3e6,asRaster=T)
dat.recon=as.data.frame(prswe.recon,xy=T)
names(dat.recon) <- c('x','y','swe')
dat.recon$model='Reconstruction'
dat.recon$yr=yr
#
prswe.phvrcn=pr.phvrcn.surface.stack[[jm]][[ryr]]
prswe.phvrcn[recon.v2.stack[[jm]]==0]=0
prswe.phvrcn=sampleRegular(prswe.phvrcn,size=3e6,asRaster=T)
dat.phvrcn=as.data.frame(prswe.phvrcn,xy=T)
names(dat.phvrcn) <- c('x','y','swe')
dat.phvrcn$model='Regression.PHVRCN'
dat.phvrcn$yr=yr
#
prswe.phv=pr.phv.surface.stack[[jm]][[ryr]]
prswe.phv[recon.v2.stack[[jm]]==0]=0
prswe.phv=sampleRegular(prswe.phv,size=3e6,asRaster=T)
dat.phv=as.data.frame(prswe.phv,xy=T)
names(dat.phv) <- c('x','y','swe')
dat.phv$model='Regression.PHV'
dat.phv$yr=yr
#
dat=rbind(dat.phv,dat.phvrcn,dat.recon)

return(dat)
}

#prdate=as.Date('April-1-2000','%B-%d-%Y')
#plotswe=rast2df(prdate)
#prdate=as.Date('April-1-2002','%B-%d-%Y')
#plotswe=rbind(plotswe,rast2df(prdate))
#prdate=as.Date('April-1-2011','%B-%d-%Y')
#plotswe=rbind(plotswe,rast2df(prdate))
plotswe$model=factor(plotswe$model,levels=c('Regression.PHV','Regression.PHVRCN','Reconstruction'))


#rm(list=c('recon.v2.stack','pr.phvrcn.surface.stack','pr.phv.surface.stack'))
load('~/Documents/R/snotel_regression/surfdf.image.april.2000.2002.2011.RData')
#
## plotswe$model=ifelse(plotswe$model=='Regression.PHV','Regression w/ PHV',
##        ifelse(plotswe$model=='Regression.PHVRCN','Regression w/ PHVRCN','Reconstruction'))
## plotswe$model=factor(plotswe$model,levels=c('Regression w/ PHV','Regression w/ PHVRCN','Reconstruction'))

facet1_names=list(
    'Regression.PHV'='Regression w/o RCN',
    'Regresion.PHVRCN'='Regression w/ RCN',
    'Reconstruction'='Reconstruction')
plot_labeller <- function(variable,value){
    if(variable=='model'){
        return(facet1_names[value])
    } else {
        return(value)
    }
 }


swe2plot=subset(plotswe,yr==2000)
plotswe$swe[plotswe$swe==0]=NA
inc=0.15
maxval=ceiling((max(swe2plot$swe,na.rm=T))*10)/10
maxval=ceiling(maxval/inc)*inc
swe2plot$cuts <- cut(swe2plot$swe,breaks=seq(0,max(subset(swe2plot,as.character(model)=='Regression.PHV' | as.character(model)=='Regression.PHVRCN')$swe,na.rm=T),inc),include.lowest=F)
#
colramppal1=colorRampPalette(brewer.pal(9,"RdYlBu"))(length(levels(swe2plot$cuts)))
colramppal2=colorRampPalette(brewer.pal(9,"RdYlBu"))(length(levels(swe2plot$cuts))*1.5)
colpal=c(colramppal1[c(1,3,6)],colramppal2[10:18])
                                        #
g=ggplot()+geom_raster(data=swe2plot,aes(x,y,fill=cuts))+
    scale_fill_manual(values = colpal,drop=F,labels=levels(swe2plot$cuts))+
    geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(0.45,0.05))+
    geom_path(data=states,aes(x=long,y=lat,group=group),col='grey20')+
    geom_point(data=data.frame(snotel.locs),aes(x=long,y=lat),colour='black',shape=19,size=1)+
    guides(fill=guide_legend('SWE (m)'),alpha=F)+
    coord_cartesian(xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
    labs(x='Longitude',y='Latitude')+
    facet_grid(model~yr,labeller=plot_labeller)+
    theme_bw()+
    theme(legend.key.size=unit(2,'lines'),
          axis.text=element_text(size=14),
          plot.title=element_text(size=18),
          axis.title=element_text(size=16),
          legend.text=element_text(size=14),
          legend.title=element_text(size=16),
          strip.text=element_text(size=14,face='bold'))        
## show(g)#
ggsave(plot=g,file='plots/swe.surf.2000.take2.png',dpi=300,width=5.5,height=10)


## --- Get difference of phv and phvrcn surfaces
swediff=subset(plotswe,as.character(model)=='Regression w/ PHV' | as.character(model)=='Regression w/ PHVRCN')
swediff.cast=dcast(swediff,x+y+yr~model,value.var='swe')
swediff.cast$diff=(swediff.cast$'Regression w/ PHVRCN' - swediff.cast$'Regression w/ PHV')/swediff.cast$'Regression w/ PHV'*100
#
numf=11
outlier=150
#
## maxval=ceiling((max(abs(range(swediff.cast$diff,na.rm=T)),na.rm=T))*10)/10
## maxval=ceiling(maxval/inc)*inc
#
swe2plot=swediff.cast
swe2plot$value=swe2plot$diff
swe2plot$value[swe2plot$diff > outlier]=outlier
swe2plot$value[swe2plot$diff < -outlier]=-outlier+1
swe2plot$cuts=cut(swe2plot$value,breaks=seq(-outlier,outlier,length.out=numf))
#
legendlabels=levels(swe2plot$cuts)
outlabel=outlier*2/(numf-1)*(numf-1)/2-outlier*2/(numf-1)
legendlabels[1]=paste('<= ',-outlabel,sep='')
legendlabels[length(levels(swe2plot$cuts))]=paste('> ',outlabel,sep='')
#
## swe2plot$cuts <- cut(swe2plot$diff,breaks=seq(-maxval,maxval,length.out=maxval*2/inc),include.lowest=F)
#
gdiff=ggplot()+geom_raster(data=swe2plot,aes(x,y,fill=cuts))+
    scale_fill_manual(values = colorRampPalette(brewer.pal(11,"RdBu"))(length(levels(swe2plot$cuts))),drop=F,labels=legendlabels)+
    geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
    scale_alpha_continuous(range=c(0.45,0.05))+
    geom_path(data=states,aes(x=long,y=lat,group=group),col='grey20')+
    geom_point(data=data.frame(snotel.locs),aes(x=long,y=lat),colour='black',shape=19,size=1)+
    guides(fill=guide_legend('SWE (%)'),alpha=F)+
    coord_cartesian(xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
    labs(x='Longitude',y='Latitude')+#, title='April 1 % Differences in Estimated SWE between the\nRegression w/ PHV and Regression w/ PHVRCN.\nBlue indicates Regression w/ PHVRCN estimated more SWE')+
    facet_wrap(~yr)+
    theme_bw()+
    theme(plot.title=element_text(size=20),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16),
          strip.text=element_text(size=12,face='bold'),
          legend.key.size=unit(1.5,'lines'),
          legend.text=element_text(size=14),
          legend.title=element_text(size=16))
## show(gdiff)
#
## ggsave(plot=gdiff,file='plots/swediff.surf.2000.2002.2011.png',dpi=300, width=10, height=5.3, units='in')


## lon.low=-112.25
## lon.high=-104.125

#
## mnth='Apr'
## mthcol=which(names(bestmodel)==mnth)
## swediffmapplot=subset(swediffmapplot,(yr==2002 & reconyear==bestmodel[2002-1999,mthcol]) | (yr==2000 & reconyear==bestmodel[2000-1999,mthcol]) | (yr==2011 & reconyear==bestmodel[2011-1999,mthcol]))
## swediffmapplot=subset(swediffmapplot,as.character(mth)==mnth)
## swediffmapplot$value=swediffmapplot$value*100
## swediffmapplot$value[swediffmapplot$value > outlier]=outlier
## swediffmapplot$value[swediffmapplot$value < -outlier]=-outlier+1
## swediffmapplot$cuts=cut(swediffmapplot$value,breaks=seq(-outlier,outlier,length.out=11))

## #
## legendlabels=levels(swediffmapplot$cuts)
## legendlabels[1]='<= -200'
## legendlabels[length(levels(swediffmapplot$cuts))]='> 200'
## #
## szmin=2
## szstep=.5
## szmax=ceiling((szmin+szstep*length(levels(swediffmapplot$cuts)))/2)
## #
## gphvrcn=ggplot()+geom_raster(data=ucorelief.df,aes(x=long,y=lat,alpha=ter))+
##     scale_alpha_continuous(range=c(1,0.5))+
##     geom_path(data=states,aes(x=long,y=lat,group=group),color='grey10')+
##     geom_point(data=subset(swediffmapplot,value < 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
##     geom_point(data=subset(swediffmapplot,value > 0),aes(x=long,y=lat,size=cuts,color=cuts,shape=cuts),alpha=0.75)+
##     scale_color_manual(values = colorRampPalette(brewer.pal(11,"RdYlBu"))(length(levels(swediffmapplot$cuts))),drop=F,labels=legendlabels)+
##     scale_shape_manual(values=c(rep(16,length(levels(swediffmapplot$cuts))/2),rep(17,length(levels(swediffmapplot$cuts))/2)),drop=F,labels=legendlabels)+
##     scale_size_manual(values=c(seq(szmax,szmin,-szstep),seq(szmin,szmax,szstep)),drop=F,labels=legendlabels)+
##     guides(alpha=F,
##            size=guide_legend('SWE Differences (%)'),
##            colour=guide_legend('SWE Differences (%)'),
##            shape=guide_legend('SWE Differences (%)'))+
##     coord_cartesian(xlim=c(-112.5,-104.25),ylim=c(33,44))+
##     theme_bw()+
##     theme(legend.key=element_rect(fill='grey80'),
##           legend.background=element_rect(fill='grey80'),
##           strip.text=element_text(size=12,face='bold'))+
##     ggtitle('April 1 Differences with Observed SNOTEL SWE\n')+
##     facet_grid(variable~yr,labeller=plot_labeller)
## ## show(gphvrcn)
## ggsave(plot=gphvrcn,'plots/allmodels.pctdiff.2002.2000.2011.png',width=15,height=15,dpi=300)


## brks=c(0.001,0.5,1,1.5,2,2.5)   
##     gpr=ggplot()+ geom_raster(data=df.sub,aes(x=long,y=lat,fill=swe))+
##         geom_path(data=states,aes(x=long,y=lat,group=group),col='black')+
##             scale_fill_gradientn(colours=brewer.pal(name='RdYlBu',n=11),na.value='grey90',breaks=brks,labels=brks)+
##                 geom_point(data=snotel.locs,aes(x=long,y=lat),size=1)+
##                     theme_minimal()+
##                         coord_cartesian(xlim=c(-112.25,-104.125),ylim=c(33,43.75))+
##                             ggtitle(paste('best regression (gaussian) surface - 01',mth,yrcnt+1999,sep=''))
