library(MASS)
library(stats)
library(lattice)
library(tidyr)

library(ggplot2)
library(gridExtra)
library(data.table)

######  plotting  ######
## read data
## The datasets are from a LUAD study. Refer to: https://doi.org/10.1016/j.cell.2020.05.043 ##
data.luad.phospho=read.table('../phoMatrixNoImpute.txt',header = T,sep='\t',quote='')
luad.phos.pro=data.luad.phospho[,-c(2,4:5)]
luad.phos=luad.phos.pro[,-c(1:2)]

allNA=which(rowSums(is.na(luad.phos))==158)
luad.phos=luad.phos[-c(allNA),]
# 22,564 sites * 79 pairs

## read data of proteomics
luad.proteome=read.csv('../LUAD_protein.csv',header = T)
luad.proteome=luad.proteome[,-c(1:3)]

### read results
result1=read.table("../Simu_results_n=5.txt",sep='\t',header=T)
result2=read.table("../Simu_results_n=10.txt",sep='\t',header=T)
result3=read.table("../Simu_results_n=30.txt",sep='\t',header=T)
result4=read.table("../Simu_results_n=50.txt",sep='\t',header=T)
result5=read.table("../Simu_results_n=80.txt",sep='\t',header=T)
result6=read.table("../Simu_results_n=100.txt",sep='\t',header=T)
result=rbind(result1,result2,result3,result4,result5,result6)
result.auc=result


####### ZI continuous models --- Infer, Conclude (Box plot) ######
source('../results_plot_functions.R')

## change the order of labels in box-plots
order_lab=as.factor(c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin','twoT','twoWilcox','SDA','SDA_robust','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p'))
fac <- as.factor(result.auc$model)
result.auc$model <- factor(fac, levels=levels(fac)[order_lab])


## 1.imputation or not
mlist=c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin')

plotbox(data=result.auc, n=10, zr=.3, mnar=NULL, fc=2, evaluation='pAUROC',modellist=mlist)
ggsave("fig2a.png", width = 15, height = 6, dpi=600)
plotbox(data=result.auc, n=30, zr=.3, mnar=NULL, fc=2, evaluation='pAUROC',modellist=mlist)
ggsave("fig2b.png", width = 15, height = 6, dpi=600)

# 1.2
plotbox(data=result.auc, n=NULL, zr=.3, mnar=0, fc=1.5, evaluation='pAUROC',modellist=mlist,textsizeM = 19,textsizeL = 22)
ggsave("fig3a.png", width = 10, height = 6, dpi=600)
plotbox(data=result.auc, n=NULL, zr=.6, mnar=0, fc=1.5, evaluation='pAUROC',modellist=mlist,textsizeM = 19,textsizeL = 22)
ggsave("fig3b.png", width = 10, height = 6, dpi=600)
# 1.3
plotbox(data=result.auc, n=10, zr=NULL, mnar=0, fc=2, evaluation='pAUROC',modellist=mlist,textsizeM = 19,textsizeL = 22)
ggsave("fig3c.png", width = 11, height = 6, dpi=600)
plotbox(data=result.auc, n=30, zr=NULL, mnar=0, fc=2, evaluation='pAUROC',modellist=mlist,textsizeM = 19,textsizeL = 22)
ggsave("fig3d.png", width = 11, height = 6, dpi=600)


## 2. Zero-Inflated models  vs  Non-Impute
mlist=c('T-test','Wilcoxon','ModT','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p')

plotbox(data=result.auc, n=30, zr=.7, mnar=NULL, fc=2, evaluation='pAUROC',modellist=mlist)
ggsave("fig4.png", width = 15, height = 6, dpi=600)


## 3. ZI-models (5)  vs  Non-Impute (3) + Impute (2) --- Finals
mlist=c('T-test','Wilcoxon','Wilcoxon-SampMin','ModT','ModT-bPCA','twoT','twoWilcox','SDA','ZIG_2p','ZILN_2p')

#3.1 ModT-bPCA
plotbox(data=result.auc, n=5, zr=.3, mnar=NULL, fc=1.5, evaluation='pAUROC',modellist=mlist,ylimit = .75)
ggsave("fig6a.png", width = 15, height = 6, dpi=600)

#3.2 ModT
plotbox(data=result.auc, n=50, zr=.9, mnar=NULL, fc=2, evaluation='pAUROC',modellist=mlist)
ggsave("fig6b.png", width = 15, height = 6, dpi=600)

#plotbox(data=result.auc, n=10, zr=.6, mnar=NULL, fc=1.5, evaluation='pAUROC',modellist=mlist)

#3.3 two-part T
plotbox(data=result.auc, n=80, zr=.7, mnar=NULL, fc=1.5, evaluation='pAUROC',modellist=mlist)
ggsave("fig6c.png", width = 15, height = 6, dpi=600)

# Figure S8
plotbox(data=result.auc, n=10, zr=.7, mnar=NULL, fc=1.5, evaluation='pAUROC',modellist=mlist)
ggsave("figs8a.png", width = 15, height = 6, dpi=600)
plotbox(data=result.auc, n=10, zr=.6, mnar=NULL, fc=1.5, evaluation='pAUROC',modellist=mlist)
ggsave("figs8b.png", width = 15, height = 6, dpi=600)
plotbox(data=result.auc, n=10, zr=.3, mnar=NULL, fc=1.5, evaluation='pAUROC',modellist=mlist)
ggsave("figs8c.png", width = 15, height = 6, dpi=600)


###### Best Method Result ######
# calculate mean of results
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


## best methods across MNAR ratios
result.mm=mean_mean(result.mean.all2,idx='pAUROC')
colnames(result.mm)=c("Nsample","fc","zr","model","E_pAUROC")

result.bb=best_best(result.mm)


####### ZI continuous models --- best method function ######
# color
#colorpalette=c("#DB2F20","#9F2218",'#EE3768',"#E68200","#E6B600","#E6EB00","#A5CF47","#00B159","#1F868F",'#4eb7d9',"#3E5CC5","#7149D7","#AF5CE2","#B981DB","#CCCCCC","#8C8C8C")
#image(x=1:16,y=1,z=as.matrix(1:16),col=colorRampPalette(colorpalette)(16))

## get all legends
tmp=result.mean.best2[(result.mean.best2$Nsample==50)&(result.mean.best2$zr==.4),]
test=tmp
test$model[1:16]=c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p')
test=test[,c(2,4,5,6)]
test$fc=as.factor(test$fc)
test$MNAR=as.factor(test$MNAR)
levels(test$model)=c('T-test','T-bPCA','T-SM','Wilcox','W-bPCA','W-SM','ModT','ModT-bPCA','ModT-SM','2part-T','2part-W','SDA','ZIG_DM','ZILN_DM','ZIG_2part','ZILN_2part')
colnames(test)=c('Var2',"Var1",'method','value')
resplot3(test,samplesize=0,zeroratio=0)

## get 10 legends
test=result.bb[(result.bb$Nsample==30),]
test$model[1:2]=c('T-bPCA','Wilcoxon-bPCA')
test=test[,c(2:5)]
test$fc=as.factor(test$fc)
test$zr=as.factor(test$zr)
levels(test$model)=c('T-test','T-bPCA','Wilcox','W-bPCA','ModT','ModT-bPCA','2part-T','2part-W','SDA','ZILN_2part')
colnames(test)=c('Var2',"Var1",'method','value')
resplot3b(test,samplesize=0)


## change the order of legends !!! important
order_lab=as.factor(c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p'))
fac <- as.factor(result.mean.best2$model)
result.mean.best2$model <- factor(fac, levels=levels(fac)[order_lab])

order_lab2=as.factor(c('T-test','T-bPCA','Wilcoxon','Wilcoxon-bPCA','ModT','ModT-bPCA','twoT','twoWilcox','SDA','ZILN_2p'))
fac2 <- as.factor(result.bb$model)
result.bb$model <- factor(fac2, levels=levels(fac2)[order_lab2])

## 42 scenarios
# Figure S2-S7
#setwd('D:/R/pAUROC_best_methods')
for (n in c(5,10,30,50,80,100)){
  for (ZR in c(.3,.4,.5,.6,.7,.8,.9)){
    tmp=result.mean.best2[(result.mean.best2$Nsample==n)&(result.mean.best2$zr==ZR),]
    temp=tmp[,c(2,4,5,6)]
    temp$fc=as.factor(temp$fc)
    temp$MNAR=as.factor(temp$MNAR)
    levels(temp$model)=c('T-test','T-bPCA','T-SM','Wilcox','W-bPCA','W-SM','ModT','ModT-bPCA','ModT-SM','2part-T','2part-W','SDA','ZIG_DM','ZILN_DM','ZIG_2part','ZILN_2part')
    colnames(temp)=c('Var2',"Var1",'method','value')
    resplot3(temp,samplesize=n,zeroratio=ZR)
  }
}

## 6 charts (best across MNAR) 
# Figure 5
#setwd('D:/R/pAUROC_best_across')
for (n in c(5,10,30,50,80,100)){
  tmp=result.bb[(result.bb$Nsample==n),]
  temp=tmp[,c(2:5)]
  temp$fc=as.factor(temp$fc)
  temp$zr=as.factor(temp$zr)
  levels(temp$model)=c('T-test','T-bPCA','Wilcox','W-bPCA','ModT','ModT-bPCA','2part-T','2part-W','SDA','ZILN_2part')
  colnames(temp)=c('Var2',"Var1",'method','value')
  resplot3b(temp,samplesize=n)
}

