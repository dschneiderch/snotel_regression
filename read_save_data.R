library(ggplot2)
library(raster)
library(lubridate)
library(R.matlab)
library(stats)
library(matlab)
library(reshape2)
# Read snotel, recon, phv data from matlab output.
setwd('~/Documents/R/snotel_regression')

data=readMat('~/Documents/MATLAB/recon_sntl_alldays/swe_phv_snotel_matrices.mat')


## PHV DATAFRAME
phvrec=data[[1]]
#Assign appropriate column names for phv
colnames(phvrec)=c('Lat','Long','Elev','Eastness','Northness','Slope','RegionalSlope','RegionalEastness','RegionalNorthness','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','NWd2ocean','SWd2ocean')
phvrec=as.data.frame(phvrec)


#get attribute data from binary matlab files.
sitenames=readMat('snotel_attributes/sitenames.mat')
sn=as.character(sitenames[[1]])
#
sitestate=readMat('snotel_attributes/sitestate.mat')
ss=as.character(sitestate[[1]])
#
sitecoords=readMat('snotel_attributes/sitecoords.mat')
slt=sitecoords[[1]][,1]
sln=sitecoords[[1]][,2]
#
stationid=readMat('snotel_attributes/stationid.mat')
sid=as.character(stationid[[1]])
#
snotelatt=data.frame(sn,sid,ss,slt,sln)

## SNOTEL DATAFRAME
datevec=as.Date('1999-10-1'):as.Date('2012-09-30')
class(datevec) <- 'Date'
sntlrec=data[[2]]
colnames(sntlrec) <- sid
snotelrec=as.data.frame(sntlrec)
snotelrec$date=datevec
snotelrec=melt(snotelrec,id='date')
colnames(snotelrec) <- c('date','stationid','swe.m')
snotelrec$lat=rep(slt,each=length(datevec))
snotelrec$long=rep(sln,each=length(datevec))
snotelrec$state=rep(ss,each=length(datevec))
snotelrec$station=rep(sn,each=length(datevec))

## RECON DATAFRAME
datevec=as.Date('1999-10-1'):as.Date('2013-09-30')
class(datevec) <- 'Date'
rcnrec=data[[3]]
colnames(rcnrec) <- sid
reconrec=as.data.frame(rcnrec)
reconrec$date=datevec
reconrec=melt(reconrec,id='date')
colnames(reconrec) <- c('date','stationid','swe.m')
reconrec$lat=rep(slt,each=length(datevec))
reconrec$long=rep(sln,each=length(datevec))
reconrec$state=rep(ss,each=length(datevec))
reconrec$station=rep(sn,each=length(datevec))


## COMBINED SWE SPATIAL DATAFRAME
swe=subset(reconrec,date<'2012-10-1')
colnames(swe) <- c('date','stationid','recon','lat','long','state','station')
swe$snotel=snotelrec$swe.m
coordinates(swe)=~long+lat
projection(swe)='+proj=longlat +datum=NAD83'
str(swe)
save(list=c('phvrec','swe','reconrec','snotelrec'),file='~/Documents/R/snotel_regression/snoteldata.RData')

### Get UCO domain terrain variables
### LOAD UCO TERRAIN VARIABLES
uco_phv=readMat('uco_variables_MASTER.mat')
names(uco_phv)=c('Lat','Long','RegionalSlope','RegionalAspect','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','Waz','NWd2ocean','NWaz','SWd2ocean','SWaz','Aspect','Slope','Elev')
### Stack each variable
Lat=raster(uco_phv$Lat,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Long=raster(uco_phv$Long,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
RegionalSlope=raster(uco_phv$RegionalSlope,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
RegionalEastness=raster(uco_phv$RegionalAspect,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
RegionalNorthness=sin(RegionalSlope)*RegionalEastness
FtprtW=raster(uco_phv$FtprtW,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Wbdiff=raster(uco_phv$Wbdiff,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
NWbdiff=raster(uco_phv$NWbdiff,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
SWbdiff=raster(uco_phv$SWbdif,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Wd2ocean=raster(uco_phv$Wd2ocean,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
NWd2ocean=raster(uco_phv$NWd2ocean,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
SWd2ocean=raster(uco_phv$SWd2ocean,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Eastness=raster(uco_phv$Aspect,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Slope=raster(uco_phv$Slope,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
Northness=Eastness*sin(Slope)
Elev=raster(uco_phv$Elev,xmn=-112.25,xmx=-104.125,ymn=33,ymx=43.75)
ucophv.stack=stack(list(Lat,Long,RegionalSlope,RegionalEastness,RegionalNorthness,FtprtW,Wbdiff,NWbdiff,SWbdiff,Wd2ocean,NWd2ocean,SWd2ocean,Eastness,Northness,Slope,Elev))
names(ucophv.stack)=c('Lat','Long','RegionalSlope','RegionalEastness','RegionalNorthness','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','NWd2ocean','SWd2ocean','Eastness','Northness','Slope','Elev')
#
#### Combine each variable into a data frame
ucophv=as.numeric(uco_phv[[1]])
for (i in seq(2,length(uco_phv))){
     ucophv=cbind(ucophv,as.numeric(uco_phv[[i]]))
}
ucophv=as.data.frame(ucophv)
colnames(ucophv)=c('Lat','Long','RegionalSlope','RegionalEastness','FtprtW','Wbdiff','NWbdiff','SWbdiff','Wd2ocean','Waz','NWd2ocean','NWaz','SWd2ocean','SWaz','Eastness','Slope','Elev')
ucophv$RegionalNorthness=sin(acos(ucophv$RegionalEastness))
ucophv$Northness=sin(acos(ucophv$Eastness))
str(ucophv)

save(file='ucophv.RData',list=c('ucophv','ucophv.stack'))
