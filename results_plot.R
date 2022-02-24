library(MASS)
library(stats)
library(lattice)
library(tidyr)

library(ggplot2)
library(gridExtra)
library(data.table)

# to draw figures from processed data
draw_from_figdata=TRUE
# thus setting simple_implement=FALSE
simple_implement=!draw_from_figdata
# or draw figures from simple_implement simulation by setting the opposite

setwd('../Phospho-MV-simulations-main')
cd=getwd()

######  plotting  ######
## read data for Figure S1
## The datasets are from a LUAD study. Refer to: https://doi.org/10.1016/j.cell.2020.05.043 ##
# data.luad.phospho=read.table('../phoMatrixNoImpute.txt',header = T,sep='\t',quote='')
# luad.phos.pro=data.luad.phospho[,-c(2,4:5)]
# luad.phos=luad.phos.pro[,-c(1:2)]
# 
# allNA=which(rowSums(is.na(luad.phos))==158)
# luad.phos=luad.phos[-c(allNA),]
# 
## read data of proteomics
# luad.proteome=read.csv('../LUAD_protein.csv',header = T)
# luad.proteome=luad.proteome[,-c(1:3)]

### read results
if (simple_implement==FALSE){
  if (draw_from_figdata==FALSE){
    # All the raw results data is too large to upload to GitHub; and stored in several files
    result1=read.table("../Simu_results_n=5.txt",sep='\t',header=T)
    result2=read.table("../Simu_results_n=10.txt",sep='\t',header=T)
    result3=read.table("../Simu_results_n=30.txt",sep='\t',header=T)
    result4=read.table("../Simu_results_n=50.txt",sep='\t',header=T)
    result5=read.table("../Simu_results_n=80.txt",sep='\t',header=T)
    result6=read.table("../Simu_results_n=100.txt",sep='\t',header=T)
    result=rbind(result1,result2,result3,result4,result5,result6)
    result.auc=result
}}else{
  result.auc=read.table(paste0(cd,'/Simu_results.txt'),sep='\t',header=T)
}


####### ZI continuous models --- Infer, Conclude (Box plot) ######
source('../results_plot_functions.R')

## change the order of labels in box-plots
order_lab=as.factor(c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin','twoT','twoWilcox','SDA','SDA_robust','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p','MSstatsPTM','MSstatsPTM-AFT'))
fac=function(result,ol=order_lab){
  f <- as.factor(result$model)
  result$model <- factor(f, levels=levels(f)[ol])
  return(result)
}

if (draw_from_figdata==TRUE){
  result_fig2=read.table(paste0(cd,"simulation_results/data_fig2.txt"),sep='\t',header=T)
  result_fig3=read.table(paste0(cd,"simulation_results/data_fig3.txt"),sep='\t',header=T)
  result_fig4=read.table(paste0(cd,"simulation_results/data_fig4.txt"),sep='\t',header=T)
  result_fig6=read.table(paste0(cd,"simulation_results/data_fig6.txt"),sep='\t',header=T)
  result_figs8=read.table(paste0(cd,"simulation_results/data_fig s8.txt"),sep='\t',header=T)
  result_fig2=fac(result_fig2)
  result_fig3=fac(result_fig3)
  result_fig4=fac(result_fig4)
  #result_fig6=fac(result_fig6)
  f <- as.factor(result_fig6$model)
  ol=as.factor(c('T-test','Wilcoxon','Wilcoxon-SampMin','ModT','ModT-bPCA','twoT','twoWilcox','SDA','ZIG_2p','ZILN_2p','MSstatsPTM','MSstatsPTM-AFT'))
  result_fig6$model <- factor(f, levels=levels(f)[ol])
  
  f <- as.factor(result_figs8$model)
  result_figs8$model <- factor(f, levels=levels(f)[ol])
}else{
  result.auc=fac(result.auc)
}

## 1.imputation or not
mlist=c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin','MSstatsPTM','MSstatsPTM-AFT')

if (draw_from_figdata){result.auc=result_fig2}
plotbox(data=result.auc, n=10, zr=.3, mnar=NULL, fc=2, evaluation='pAUROC',modellist=mlist)
ggsave("fig2a_revised.png", width = 9, height = 5, dpi=600)
if (simple_implement==FALSE){
  # because the sample size of 30 not in simple_implement simulation; the following is the same
plotbox(data=result.auc, n=30, zr=.3, mnar=NULL, fc=2, evaluation='pAUROC',modellist=mlist)
ggsave("fig2b_revised.png", width = 9, height = 5, dpi=600)
}
# Fig2: v2
plotbox2(data=result.auc, n=c(10,30), zr=.3, fc=2, evaluation='pAUROC',modellist=mlist,No.fig='fig2',textsizeL = 15)
ggsave("fig2_revised.png", width = 6, height = 5, dpi=600)

# 1.2
if (draw_from_figdata){result.auc=result_fig3}
# Fig3: v2
plotbox2(data=result.auc, n=c(5,10,30,50,80,100), zr=c(.3,.6), mnar=0, fc=1.5, evaluation='pAUROC',modellist=mlist,No.fig='fig3ab',textsizeL = 15)
ggsave("fig3ab_revised.png", width = 6, height = 5, dpi=600)
plotbox2(data=result.auc, n=c(10,30), zr=c(.3,.4,.5,.6,.7,.8), mnar=0, fc=2, evaluation='pAUROC',modellist=mlist,No.fig='fig3cd',textsizeL = 15)
ggsave("fig3cd_revised.png", width = 6, height = 5, dpi=600)


## 2. Zero-Inflated models  vs  Non-Impute
mlist=c('T-test','Wilcoxon','ModT','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p','MSstatsPTM')
if (draw_from_figdata){result.auc=result_fig4}

if (simple_implement==FALSE){
plotbox(data=result.auc, n=30, zr=.7, mnar=NULL, fc=2, evaluation='pAUROC',modellist=mlist)
ggsave("fig4_legend.png", width = 9, height = 5, dpi=600)
}
# Fig4: v2
plotbox2(data=result.auc, n=30, zr=.7, fc=2, evaluation='pAUROC',modellist=mlist,No.fig='fig4',textsizeL = 15)
ggsave("fig4_revised.png", width = 6, height = 3, dpi=600)


## 3. ZI-models (5)  vs  Non-Impute (4) + Impute (3) --- Finals
mlist=c('T-test','Wilcoxon','Wilcoxon-SampMin','ModT','ModT-bPCA','twoT','twoWilcox','SDA','ZIG_2p','ZILN_2p','MSstatsPTM','MSstatsPTM-AFT')
if (draw_from_figdata){result.auc=result_fig6}
# 3.1 ModT-bPCA; 3.2 ModT; 3.3 two-part T
# Fig6: v2
plotbox2(data=result.auc, n=NULL, zr=NULL, fc=NULL, evaluation='pAUROC',modellist=mlist,No.fig='fig6',textsizeL = 15)
ggsave("fig6_revised.png", width = 6, height = 8, dpi=600)
plotbox2(data=result.auc, n=NULL, zr=NULL, fc=NULL, evaluation='pAUROC',modellist=mlist,No.fig='fig6',textsizeL = 15,ylimit = .75)
ggsave("fig6_.75.png", width = 6, height = 8, dpi=600)


# Figure S9
if (draw_from_figdata){result.auc=result_figs8}
plotbox(data=result.auc, n=10, zr=.7, mnar=NULL, fc=1.5, evaluation='pAUROC',modellist=mlist)
ggsave("figs9a_legend.png", width = 9, height = 6, dpi=600)

plotbox2(data=result.auc, n=10, zr=c(0.7,0.6,0.3), fc=1.5, evaluation='pAUROC',modellist=mlist,No.fig='figsi',textsizeL = 15)
ggsave("figs9.png", width = 6, height = 8, dpi=600)
plotbox2(data=result.auc, n=10, zr=c(0.7,0.6,0.3), fc=1.5, evaluation='pAUROC',modellist=mlist,No.fig='figsi',textsizeL = 15,ylimit = .75)
ggsave("figs9_.75.png", width = 6, height = 8, dpi=600)

# verified for simple_implement==F;draw_from_figdata==T
# verified for simple_implement==T;draw_from_figdata==F


###### Best Method Result ######
# calculate mean of results
if (draw_from_figdata==FALSE){
result.mean.all=matrix(NA,0,10)
for (n in c(5,10,30,50,80,100)){
  for (zr in c(.3,.4,.5,.6,.7,.8,.9)){
    for (mnar in c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)){
      for (fc in c(1.2,1.5,1.8,2,3,5,10)) {
        tmp=result.auc[(result.auc[,1]==n)&(result.auc[,3]==zr)&(result.auc[,4]==mnar)&(result.auc[,2]==fc),]
        result.mean.all=rbind(result.mean.all,method_mean(tmp))
      }
    }
  }
}


# already checked---is consistent.
write.table(result.mean.all,"D:/datasets/simu_/Simu_results_mean.txt", quote=FALSE, sep="\t", row.names=F, col.names=T, append=FALSE)

result.mean.all=read.table('D:/datasets/simu_/Simu_results_mean.txt',sep='\t',header=T)

result.mean.all2=result.mean.all[result.mean.all$model!='SDA_robust',]
write.table(result.mean.all2,"D:/datasets/simu_/Simu_results_mean2.txt", quote=FALSE, sep="\t", row.names=F, col.names=T, append=FALSE)

### select best methods

## best methods under 3,234 scenarios
result.mean.best=best_mean(result.mean.all,idx='pAUROC')
# NA row due to pAUROC NA
result.mean.best=na.omit(result.mean.best)
#sum(is.na(result.mean.all$pAUROC))==554
#554+3234==3788
write.table(result.mean.best,"D:/datasets/simu_/Simu_results_mean_best.txt", quote=FALSE, sep="\t", row.names=F, col.names=T, append=FALSE)

result.mean.best2=best_mean(result.mean.all2,idx='pAUROC')
write.table(result.mean.best2,"D:/datasets/simu_/Simu_results_mean_best2.txt", quote=FALSE, sep="\t", row.names=F, col.names=T, append=FALSE)

}else{
  # if draw_from_figdata==TRUE
  result.mean.all2=read.table(paste0(cd,'/simulation_results/Simu_results_MSstatsPTM_mean_all.txt'),sep='\t',header=T)
  result.mean.best3=read.table(paste0(cd,'/simulation_results/Simu_results_mean_best3.txt'),sep='\t',header=T)
}


## best methods across MNAR ratios
result.mm=mean_mean(result.mean.all2,idx='pAUROC')
colnames(result.mm)=c("Nsample","fc","zr","model","E_pAUROC")

result.bb=best_best(result.mm)


####### ZI continuous models --- best method function ######
# color
#colorpalette=c("#DB2F20","#9F2218",'#EE3768',"#E68200","#E6B600","#E6EB00","#A5CF47","#00B159","#1F868F",'#4eb7d9',"#3E5CC5","#7149D7","#AF5CE2","#B981DB","#CCCCCC","#8C8C8C")
#image(x=1:16,y=1,z=as.matrix(1:16),col=colorRampPalette(colorpalette)(16))

## change the order of legends !!! important
order_lab=as.factor(c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p','MSstatsPTM')) #,'MSstatsPTM-AFT'))
fac <- as.factor(result.mean.best3$model)
result.mean.best3$model <- factor(fac, levels=levels(fac)[order_lab])

order_lab2=as.factor(c('T-test','T-bPCA','Wilcoxon','Wilcoxon-bPCA','ModT','ModT-bPCA','twoT','twoWilcox','SDA','ZILN_2p','MSstatsPTM'))
fac2 <- as.factor(result.bb$model)
result.bb$model <- factor(fac2, levels=levels(fac2)[order_lab2])

## get all legends
tmp=result.mean.best3[(result.mean.best3$Nsample==50)&(result.mean.best3$zr==.4),]
test=tmp
test$model[1:17]=c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p','MSstatsPTM')
test=test[,c(2,4,5,6)]
test$fc=as.factor(test$fc)
test$MNAR=as.factor(test$MNAR)
levels(test$model)=c('T-test','T-bPCA','T-SM','Wilcox','W-bPCA','W-SM','ModT','ModT-bPCA','ModT-SM','2part-T','2part-W','SDA','ZIG_DM','ZILN_DM','ZIG_2part','ZILN_2part','MSstatsPTM')
colnames(test)=c('Var2',"Var1",'method','value')
resplot3(test,samplesize=0,zeroratio=0)

## get 10 legends
test=result.bb[(result.bb$Nsample==30),]
test$model[1:3]=c('T-bPCA','Wilcoxon-bPCA','Wilcoxon')
test=test[,c(2:5)]
test$fc=as.factor(test$fc)
test$zr=as.factor(test$zr)
levels(test$model)=c('T-test','T-bPCA','Wilcox','W-bPCA','ModT','ModT-bPCA','2part-T','2part-W','SDA','ZILN_2part','MSstatsPTM')
colnames(test)=c('Var2',"Var1",'method','value')
resplot3b(test,samplesize=0)

## 42 scenarios
# Figure S3-S8
setwd(paste0(cd,'/pAUROC_best_methods'))
for (n in c(5,10,30,50,80,100)){
  for (ZR in c(.3,.4,.5,.6,.7,.8,.9)){
    tmp=result.mean.best3[(result.mean.best3$Nsample==n)&(result.mean.best3$zr==ZR),]
    temp=tmp[,c(2,4,5,6)]
    temp$fc=as.factor(temp$fc)
    temp$MNAR=as.factor(temp$MNAR)
    levels(temp$model)=c('T-test','T-bPCA','T-SM','Wilcox','W-bPCA','W-SM','ModT','ModT-bPCA','ModT-SM','2part-T','2part-W','SDA','ZIG_DM','ZILN_DM','ZIG_2part','ZILN_2part','MSstatsPTM')
    colnames(temp)=c('Var2',"Var1",'method','value')
    resplot3(temp,samplesize=n,zeroratio=ZR)
  }
}

## 6 charts (best across MNAR) 
# Figure 5
setwd(paste0(cd,'/pAUROC_best_across'))
for (n in c(5,10,30,50,80,100)){
  tmp=result.bb[(result.bb$Nsample==n),]
  temp=tmp[,c(2:5)]
  temp$fc=as.factor(temp$fc)
  temp$zr=as.factor(temp$zr)
  levels(temp$model)=c('T-test','T-bPCA','Wilcox','W-bPCA','ModT','ModT-bPCA','2part-T','2part-W','SDA','ZILN_2part','MSstatsPTM')
  colnames(temp)=c('Var2',"Var1",'method','value')
  resplot3b(temp,samplesize=n)
}

####### quantify the decrease in pAUC by AFT ######
aft=result.mean.all2[(result.mean.all2$model %in% c('MSstatsPTM','MSstatsPTM-AFT'))&(result.mean.all2$MNAR<1),]
# if consider mnar==1, sometimes AFT improve by 70%
aft2=c()
for (n in c(5,10,30,50,80,100)){
  for (zr in c(.3,.4,.5,.6,.7,.8,.9)){
    for (fc in c(1.2,1.5,1.8,2,3,5,10)) {
      for (mnar in c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9)){
        tmp1=aft[(aft$Nsample==n)&(aft$zr==zr)&(aft$fc==fc)&(aft$MNAR==mnar)&(aft$model=='MSstatsPTM'),]
        tmp2=aft[(aft$Nsample==n)&(aft$zr==zr)&(aft$fc==fc)&(aft$MNAR==mnar)&(aft$model=='MSstatsPTM-AFT'),]
        aft2=c(aft2,(tmp1$pAUROC-tmp2$pAUROC)/tmp1$pAUROC)
      }
    }
  }
}
max(aft2)
summary(aft2)
sum(aft2<0)/length(aft2)

