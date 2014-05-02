setwd('~/Documents/R/snotel_regression')
library(gstat)
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
#
load('cvdrop1.xrecon.regression.gauss.snotelunscaled.skill.image.RData')



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
gr2box=ggplot()+geom_boxplot(data=r2df.sameyr,aes(x=as.factor(yr),y=r2recon,col='Reconstruction'),show_guide=F)+
    geom_boxplot(data=r2df,aes(as.factor(yr),r2phvrcn,col='Regression with RCN Ensemble'),alpha=0.5,outlier.shape=3,show_guide=T)+#
    geom_boxplot(data=r2df,aes(x=as.factor(yr),y=r2phv,col='Regression w/o RCN'),show_guide=F)+
    geom_point(data=r2df.sameyr,aes(x=as.factor(yr),y=r2phvrcn,col='Regression w/ RCN from  Same Year'),show_guide=T)+
    scale_color_manual(values=c('blue','red','green','black'))+
    facet_grid(~mth)+
    xlab('Year Modeled')+
    ylab(expression(r^2))+
    ylim(c(0,0.8))+
    theme_bw()+
    guides(color=guide_legend('Model',override.aes=list(shape=c(3,20,3,3))))+
    theme(legend.position = c(0.9,.17),
          legend.background=element_rect(color='black'),
          axis.text.y=element_text(size=12),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=14),
          strip.text=element_text(size=14))#+
    ## labs(title=expression(paste('Ensemble of Cross-Validated Skill Scores (',r^2,')',sep='')))
#
show(gr2box)
 ggsave(plot=gr2box,file='plots/r2.xvaldrop1.boxplot.gauss.png',dpi=300,width=15.5,height=5)


### Boxplot realtime r2 by month
r2df.realtime=subset(r2df,yr>reconyear)
r2df.sameyr=subset(r2df,yr==reconyear)
gr2rtbox=ggplot()+#geom_boxplot(data=r2df.realtime,aes(x=as.factor(yr),y=r2recon,col='Reconstruction'),outlier.shape=3)+
    geom_boxplot(data=r2df.realtime,aes(as.factor(yr),r2phvrcn,col='Regression with Real-time RCN Ensemble'),outlier.shape=3,fill=NA)+
    geom_boxplot(data=r2df,aes(x=as.factor(yr),y=r2phv,col='Regression w/o RCN'))+
  #  geom_point(data=r2df.sameyr,aes(x=as.factor(yr),y=r2phvrcn,col='Regression w/ RCN from Same Year'),shape=6)+
    scale_color_manual(values=c('green','black'))+#'red','green'
    facet_grid(~mth)+
    xlab('Year Modeled')+
    ylab(expression(r^2))+
    ylim(c(0,0.8))+
    theme_bw()+
    guides(color=guide_legend('Model'))+#,override.aes=list(shape=c(6,1,1))))+
    theme(legend.position = c(0.89,.11),
          legend.background=element_rect(color='black'),
          axis.text.y=element_text(size=12),
          axis.text.x=element_text(size=10),
          axis.title=element_text(size=14),
          strip.text=element_text(size=14))
          #+
    ## labs(title=expression(paste('Real-time Ensemble of Cross-Validated Skill Scores (',r^2,')',sep='')))
#
## show(gr2rtbox)
ggsave(plot=gr2rtbox,file='plots/r2.xval.realtime.boxplot.gauss.png',dpi=300,width=15.5,height=5)


## Boxplot r2 mlr
r2df.modelm=melt(r2df,id=c('mth','yr','reconyear'))
r2df.best=data.frame()
rcnyr=read.table('best_recon.phvrcn.drop1.txt')
for(yr in 2000:2011){
    for(mthi in 1:3){
        temp=subset(r2df.modelm,yr==yr & reconyear==rcnyr[yr-1999,mthi])
        r2df.best=rbind(r2df.best,temp)
    }
}
r2df.phvrcnall=subset(r2df.modelm,variable=='r2phvrcn')

#
ggplot()+
    geom_boxplot(data=r2df.best,aes(x=variable,y=value,color='c1',linetype='a'),position='dodge',fill=NA)+
    geom_boxplot(data=r2df.phvrcnall,aes(x=variable,y=value,color='c2',linetype='b'),position='dodge',fill=NA)+
    scale_x_discrete(labels=c('Regression\nw/o RCN','Regression\nw/ RCN','Recon'))+
    scale_colour_manual(values=c('grey10','red'),labels=c('Best Available Model','Full Model Ensemble'))+
    scale_linetype_manual(values=c(1,2),labels=c('Best Available Model','Full Model Ensemble'))+
    guides(color=guide_legend('Model Scenario'),
           linetype=guide_legend('Model Scenario'))+
    xlab('')+
    ylab(expression(r^2))+
    theme_classic()+
    theme(legend.position=c(.15,.9))
 ## ggsave('plots/cv.mlrphv.pdf',width=170,height=100,units='mm')




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
