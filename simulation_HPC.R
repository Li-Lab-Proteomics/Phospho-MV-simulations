library(MASS)
library(stats)
if(!require('lattice')) install.packages('lattice', repos="https://cloud.r-project.org/")
library(lattice)
if(!require('tidyr')) install.packages('tidyr', repos="https://cloud.r-project.org/")
library(tidyr)
if(!require('invgamma')) install.packages('invgamma', repos="https://cloud.r-project.org/")
library(invgamma)
if(!require('BiocManager')) install.packages('BiocManager', repos="https://cloud.r-project.org/")
if(!require('limma')) BiocManager::install('limma')
library(limma)
if(!require('SDAMS')) BiocManager::install('SDAMS')
library(SDAMS)
if(!require('pcaMethods')) BiocManager::install('pcaMethods')
library(pcaMethods)
# modified SDA
if(!require('trust')) install.packages('trust', repos="https://cloud.r-project.org/")
library(trust)
# MSstatsPTM
if(!require('reshape2')) install.packages('reshape2', repos="https://cloud.r-project.org/")
library(reshape2)
if(!require('survival')) install.packages('survival', repos="https://cloud.r-project.org/")
library(survival)

setwd('../Phospho-MV-simulations-main')
cd=getwd()
simple_implement=TRUE
# number of simulation rounds
nround=ifelse(simple_implement,5,500)

# modified SDA
## mod_SDA is a modified version of 'SDA' implemented in the package 'SDAMS'. ##
source(paste(cd,'/mod_SDA/SDA_control.R',sep=''))
source(paste(cd,'/mod_SDA/data_clean.R',sep=''))
source(paste(cd,'/mod_SDA/SDA.R',sep=''))
source(paste(cd,'/mod_SDA/result_summary.R',sep=''))
# Mean-Based 2 Group Comparison without Covariate Adjustment
source(paste(cd,'/ZI_DM_tests.R',sep=''))

######  functions  ######

gen_simudata=function(pepsize=10000, Nsample=79, rowmu=18.3, rowmusd=0.2545187, DElist, FC=1){
  if (!is.integer(DElist)){
    print('DE sites index should be integer!')
  }else if (length(DElist)>pepsize){
    print('DE index mismatch.')
  }else{
    control=matrix(NA,nrow=pepsize,ncol=Nsample,byrow=F)
    case=matrix(NA,nrow=pepsize,ncol=Nsample,byrow=F)
    rmu=rnorm(pepsize,rowmu,rowmusd)
    rowsd=rinvgamma(pepsize,shape=7.393877,scale=0.4040921)
    for (m in 1:pepsize) {
      control[m,]=rlnorm(Nsample,meanlog = rmu[m],sdlog = rowsd[m])
      
      if (m %in% DElist){
        case[m,]=rlnorm(Nsample,meanlog = (rmu[m]+log(FC)),sdlog = rowsd[m])
      } else{
        case[m,]=rlnorm(Nsample,meanlog = rmu[m],sdlog = rowsd[m])
      }
      
    }
    return(list(control,case))
  }
}
MNAR_Threshold<-function(dataset,miss_prop=0.5,MNAR_rate=1,sd_q=0.01,log=TRUE){
  dataset<-as.matrix(dataset)
  samplenumber<-round(nrow(dataset)*ncol(dataset)*miss_prop*MNAR_rate)
  if (miss_prop==0) {
    return(dataset)
  }
  if (miss_prop>0) {
    if (log==TRUE){
      dataset.log=log(dataset)
      mean_u<-quantile(dataset.log,na.rm = T,miss_prop)
      threshold<-rnorm(length(dataset.log),mean = mean_u,sd = sd_q)
      id_1=which(dataset.log<threshold)
      if (samplenumber>length(id_1)){
        print('Threshold MNAR pool less than sampling number!')
        id_2=id_1
      } else{
        id_2=sample(id_1,samplenumber)
      }
      dataset[id_2]<-NA
      return(dataset)
    } else{
      mean_u<-quantile(dataset,na.rm = T,miss_prop)
      threshold<-rnorm(length(dataset),mean = mean_u,sd = sd_q)
      id_1=which(dataset<threshold)
      if (samplenumber>length(id_1)){
        print('Threshold MNAR pool less than sampling number!')
        id_2=id_1
      } else{
        id_2=sample(id_1,samplenumber)
      }
      dataset[id_2]<-NA
      return(dataset)
    }
  } else{
    print('MNAR ratio should not be negative!')
  }
}
miss_CN_mean<-function(dataset, MNAR_rate=0.8, missrate=0.1, log=TRUE, method='Thmat', logis=0, sd_q=0.01){
  elements<-nrow(dataset)*ncol(dataset)
  MCAR_rate <- 1-MNAR_rate
  if (method=='logis'){
    dataset1<-MNAR_logistic_mean(dataset,logis_a = logis,miss_prop = missrate*MNAR_rate,log=T)
  } else if (method=='Thmat'){
    dataset1=MNAR_Threshold(dataset, miss_prop=missrate, MNAR_rate = MNAR_rate, sd_q=sd_q, log=T)
  }
  almiss=which(!is.na(dataset1))
  rdm<-sample(almiss,elements*missrate*MCAR_rate)
  dataset1[rdm]<-NA
  return(dataset1) 
}

del_setNArow=function(dataset){
  if (ncol(dataset)%%2!=0) print('Ncol of dataset is not 2x !')
  n=ncol(dataset)/2
  x1=rowSums(is.na(dataset[,1:n]))
  x2=rowSums(is.na(dataset[,(n+1):(2*n)]))
  if (max(x1,x2)==n){
    index=unique(c(which(x1==n),which(x2==n)))
  } else{
    index=NA
  }
  return(index)
}
NAcol=function(dataset){
  n=nrow(dataset)
  y=colSums(is.na(dataset))
  if (sum(y==n)==0){
    return(FALSE)
  } else{
    NAsamp=which(y==n)
    return(list(TRUE,NAsamp,dataset[,-NAsamp]))
  }
}
transfer_na2zero=function(set){
  set.zr=set
  set.zr[is.na(set)]=0
  return(set.zr)
}
SampMin=function(x){
  for (i in 1:dim(x)[2])
    x[,i][is.na(x[,i])]=min(x[,i],na.rm=T)
  return(x)
}

sda_feat_ind=function(featstr=sda$feat.names){
  x=unlist(strsplit(featstr,split='_'))
  ind=as.numeric(x[seq(0,length(x),2)])
  return(ind)
}

modT=function(data,ncontrol,ncase){
  design=model.matrix(~-1+factor(c(rep(1,ncontrol),rep(2,ncase))))
  colnames(design)=c('control','case')
  fit=lmFit(data,design)
  contrast.matrix=makeContrasts(case-control,levels=design)
  fit1=contrasts.fit(fit,contrast.matrix)
  fit2=eBayes(fit1)
  pvalue.modT <- p.adjust(fit2$p.value, method = "BH")
  return(pvalue.modT)
}

## The function of calculating two-part statistics is referred to the paper of Taylor & Pollard (2009). ##
## Refer to: DOI: 10.2202/1544-6115.1425 ##
TwoPart <- function(data, group, test="t.test", point.mass=0){
  Index1 <- c(group==1)
  Group1 <- data[Index1]
  Group0 <- data[!Index1]
  n1 <- length(Group1)
  n2 <- length(Group0)
  obs <- c(n1, n2)
  success <- c(sum(Group1!=point.mass), sum(Group0!=point.mass))
  pointmass <- obs-success
  if (sum(success)==0) {
    T2 <- 0
    B2 <- 0
  }
  else if ((success[1]==0)|(success[2]==0)) {
    T2 <- 0
    B2 <- prop.test(pointmass, obs)$statistic
  }
  else if ((success[1]==1)|(success[2]==1)){
    T2 <- 0
    B2 <- prop.test(pointmass, obs)$statistic
  }
  else {
    uniq1 <- length(unique(Group1[Group1!=point.mass]))
    uniq2 <- length(unique(Group0[Group0!=point.mass]))
    if ((uniq1 < 2) & (uniq2 < 2)){
      T2 <- 0
      if (sum(pointmass)==0)
        B2 <- 0
      else
        B2 <- prop.test(pointmass, obs)$statistic
    }
    else if (sum(pointmass)==0){
      B2 <- 0
      if (test=="t.test")
        T2 <- t.test(data~group)$statistic^2
      if (test=="wilcoxon") {
        W <- wilcox.test(data~group, exact=FALSE)$statistic
        mu <- (n1*n2)/2
        sigma <- sqrt((n1*n2*(n1+n2+1))/12)
        T2 <- ((abs(W-mu)-0.5)/sigma)^2
      }
    }
    else {
      B2 <- prop.test(pointmass, obs)$statistic
      contIndex <- data!=point.mass
      cont <- data[contIndex]
      cGroup <- group[contIndex]
      n1c <- sum(cGroup==1)
      n2c <- sum(cGroup==0)
      if (test=="t.test")
        T2 <- t.test(cont~cGroup)$statistic^2
      if (test=="wilcoxon") {
        W <- wilcox.test(cont~cGroup, exact=FALSE)$statistic
        mu <- (n1c*n2c)/2
        sigma <- sqrt((n1c*n2c*(n1c+n2c+1))/12)
        T2 <- ((abs(W-mu)-0.5)/sigma)^2
      }
    }
  }
  X2 <- B2+T2
  if ((T2==0)|(B2==0)) {
    X2pv <- 1-pchisq(X2,1)
  } else {
    X2pv <- 1-pchisq(X2,2)
  }
  ans <- list(statistic=as.numeric(X2), pvalue=as.numeric(X2pv))
  return(ans)
}

ziDMtest=function(x,y,nit=25,method='zig'){
  xx=matrix(c(rep(1,nsample*2),x),ncol=2,byrow=F)
  if (method=='zig'){
    mod_g=zig.glm(x=xx,y=y,nit=nit)
    zig_test=zig.uni(mod=mod_g,n0=nsample,n1=nsample)
    return(zig_test$DM.p)
  } else if (method=='ziln'){
    mod_ln=ziln.glm(x=xx,y=y,nit=nit)
    ziln_test=ziln.uni(mod=mod_ln,n0=nsample,n1=nsample)
    return(ziln_test$DM.p[1,1])
  }
}
## ZIG/ZILN 2-part tests are based on the model estimation codes from Mills (2013)'s PhD thesis. ##
## Refer to: https://doi.org/10.17077/etd.7v3bafbd ##
ziWaldTest=function(x,y,nit=25,method='zig'){
  xx=matrix(c(rep(1,nsample*2),x),ncol=2,byrow=F)
  
  Index1 <- c(x==1)
  Group1 <- y[Index1]
  Group0 <- y[!Index1]
  
  n1 <- length(Group1)
  n2 <- length(Group0)
  obs <- c(n1, n2)
  success <- c(sum(Group1!=0), sum(Group0!=0))
  pointmass <- obs-success
  if (sum(success)==0) {
    W1 <- 0
    pv <- NA
  }
  else if ((success[1]==0)|(success[2]==0)) {
    W1 <- 0
    pv <- NA
  }
  # else if ((success[1]==1)|(success[2]==1)){
  #   T2 <- 0
  #   glmb=glm.bin(x=xx,y=y,nit=nit)
  #   B2 <- prop.test(pointmass, obs)$statistic
  # }
  else if ((pointmass[1]==0)|(pointmass[2]==0)){
    if (method=='zig'){
      mod_g=glm.gamma(x=xx,y=y,nit=nit)
      tau1=mod_g$b.c[2]
      sd_tau1=mod_g$sd1.c.mle[2]
      W1=tau1^2/sd_tau1^2
      pv <- 1-pchisq(W1,1)
    } else if (method=='ziln'){
      mod_ln=lm.lnorm(x=xx,y=y)
      tau1=mod_ln$b.c[2]
      sd_tau1=mod_ln$sdb.c[2]
      W1=tau1^2/sd_tau1^2
      pv <- 1-pchisq(W1,1)
    }
  }
  else {
    if (method=='zig'){
      mod_zig=zig.glm(x=xx,y=y,nit=nit)
      beta1=mod_zig$b.b[2]
      sd_beta1=mod_zig$sd1.b[2]
      tau1=mod_zig$b.c[2]
      sd_tau1=mod_zig$sd1.c.mle[2]
      W1=beta1^2/sd_beta1^2+tau1^2/sd_tau1^2
      pv <- 1-pchisq(W1,2)
      #return(list(statistic=as.numeric(W1), pvalue=as.numeric(pv)))
    } else if (method=='ziln'){
      mod_ziln=ziln.glm(x=xx,y=y,nit=nit)
      beta1=mod_ziln$b.b[2]
      sd_beta1=mod_ziln$sd1.b[2]
      tau1=mod_ziln$b.c[2]
      sd_tau1=mod_ziln$sd1.c[2]
      W1=beta1^2/sd_beta1^2+tau1^2/sd_tau1^2
      pv <- 1-pchisq(W1,2)
    }
  }
  return(list(statistic=as.numeric(W1), pvalue=as.numeric(pv)))
}
## MSstatsPTM. Refer to: Kohler, D.;  Tsai, T.-H.;  Huang, T.;  Staniak, M.;  Choi, M.; Vitek, O. MSstatsPTM: Statistical Characterization of Post-translational Modifications. 2021. R package version 1.2.4.
msstats=function(input, impute=FALSE){
  n_features=nrow(input)
  p=c()
  if (impute==TRUE){
    input=AFT_impute_summarize(input)
    for (i in seq_len(n_features)){
      input_feature=input[input$FEATURE==i,]
      full_fit = lm(newABUNDANCE ~ GROUP, data = input_feature)
      p[i]=summary(full_fit)$coefficients[2,4]
    }
  }else{
    input=MSstats_input(input)
    for (i in seq_len(n_features)){
      input_feature=input[input$FEATURE==i,]
      full_fit = lm(value ~ GROUP, data = input_feature)
      #df_full = full_fit[["df.residual"]]
      p[i]=summary(full_fit)$coefficients[2,4]
    }
  }
  return(p)
}
MSstats_input=function(raw,logTrans=2){
  if (logTrans==2){
    raw=log2(raw)
  }
  tmp=melt(raw,varnames=c('FEATURE','RUN'))
  tmp$RUN=factor(tmp$RUN)
  tmp$FEATURE=factor(tmp$FEATURE)
  tmp['cen']=ifelse(is.na(tmp$value),0,1)
  tmp['GROUP']=factor(rep(c('Ctrl','Tumor'),each=nrow(tmp)/2))
  return(tmp)
}
AFT_impute=function(dataset){
  input=MSstats_input(dataset)
  missingness_filter = is.finite(input$value)
  n_total = nrow(input[missingness_filter, ])
  n_features = data.table::uniqueN(input[missingness_filter, 'FEATURE'])
  n_runs = data.table::uniqueN(input[missingness_filter, 'RUN'])
  countdf = n_total  < n_features + n_runs - 1
  #set.seed(100)
  if (n_features == 1L) {
    fit = survreg(Surv(value, cen, type = "left") ~ RUN,
                  data = input, dist = "gaussian")    
  } else {
    if (countdf) {
      fit = survreg(Surv(value, cen, type = "left") ~ RUN,
                    data = input, dist = "gaussian")
    } else {
      fit = survreg(Surv(value, cen, type = "left") ~ FEATURE + RUN,
                    data = input, dist = "gaussian")
    }
  }
  input['predicted']=predict(fit, newdata = input)
  input['predicted']=ifelse(input$cen==0, input$predicted, NA)
  input['newABUNDANCE']=ifelse(input$cen==0, input$predicted, input$value)
  return(input)
}
AFT_impute_summarize=function(data,num_fea=1){
  n_features=nrow(data)
  #n_split=n_features%/%num_pro
  input=NULL
  if (num_fea==1){
    for (i in seq_len(n_features)){
      dataset=t(data.frame(data[i,]))
      rownames(dataset)=c(i)
      input=rbind(input,AFT_impute(dataset))
    }
  }
  return(input)
}

######  random seeds  ######
set.seed(12)
seed500=sample(1:100000,500)
seed500=sort(seed500)
HPCseed=data.frame(1:500,seed500)
#write.table(HPCseed, paste(cd, "/simu_/seed.txt", sep=""), quote=FALSE, sep="\t", row.names=F, col.names=F, append=FALSE)

######  start simulation  ######

sitesize=10000
sample.size = c(5,10,30,50,80,100)
zero.ratio <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
mnar.ratio = seq(0,1,0.1)
fclist <- c(1.2, 1.5, 1.8, 2.0, 3.0, 5.0, 10)
proportion=.5
if (simple_implement==TRUE){
  sample.size = c(5,10)
  fclist <- c(1.5, 2.0)
}

# simulation framework designed for MNAR - site level
dir.create(paste(cd, "/simu_results", sep=""))
for (round in 1:nround){
  dir.create(paste(cd, "/simu_results/round", round, sep=""))
  # set seed
  seed500=HPCseed$seed500[round]
  
  for (nsample in sample.size) {
    set.seed(seed500)
    group=c(rep(0,nsample),rep(1,nsample))
    
    for (fc in fclist) {
      DEindex=sample(sitesize,sitesize*proportion)
      truth=rep(0,sitesize)
      truth[DEindex]=1
      
      # generate non-zero intensity for case matrix
      datalist=gen_simudata(pepsize=sitesize, Nsample=nsample, DElist=DEindex, FC=fc)
      control=datalist[[1]]
      case=datalist[[2]]

      data=cbind(control,case)
      
      # define MNAR & MAR cells
      for (zr in zero.ratio) {
        for (mnar in mnar.ratio){
          cat("\nconstruct analyze file under sample size=", nsample, ", fc=", fc, ", zero ratio=", zr, ", MNAR ratio=", mnar, ":\n")
          data.na=miss_CN_mean(dataset=data, MNAR_rate = mnar, missrate=zr, log=T, method='Thmat')
          #index.na.df=del_allNArow(data.na)
          index.na=del_setNArow(data.na)

          truth.na=truth
          if (!sum(is.na(index.na))){
            data.na=data.na[-index.na,]
            truth.na=truth.na[-index.na]
          }
          if (length(truth.na)<1000) print('Too Many NA rows!')
          if (length(truth.na)==0){
            print('All sites are NA (missing)!')
            break
          }
          
          # imputation
          colsamp=NAcol(data.na)
          BPCA=TRUE
          if (!colsamp[[1]]){
            ### No NA col ###
            pc=try(pca(data.na,method='bpca',nPcs=2,verbose=F),silent = T)
            if ('try-error' %in% class(pc)){
              #write.table(as.data.frame(pc[1]), paste(cd, "/simu_results/Error_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
              pcsp=strsplit(pc[1],'computationally')[[1]]
              if (length(pcsp)==2){
                ## mod_bPCA is a modified version of bPCA imputation method implemented in the package 'pcaMethods'. ##
                source(paste(cd,'/R/mod_bpca/pca.R',sep=''))
                source(paste(cd,'/R/mod_bpca/bpca.R',sep=''))
                source(paste(cd,'/R/mod_bpca/BPCA_initmodel.R',sep=''))
                source(paste(cd,'/R/mod_bpca/BPCA_dostep.R',sep=''))
                source(paste(cd,'/R/mod_bpca/checkData.R',sep=''))
                source(paste(cd,'/R/mod_bpca/prep.R',sep=''))
                source(paste(cd,'/R/mod_bpca/repmat.R',sep=''))
                source(paste(cd,'/R/mod_bpca/errorHierarchic.R',sep=''))
                source(paste(cd,'/R/mod_bpca/derrorHierarchic.R',sep=''))
                source(paste(cd,'/R/mod_bpca/AllClasses.R',sep=''))
                source(paste(cd,'/R/mod_bpca/AllGenerics.R',sep=''))
                source(paste(cd,'/R/mod_bpca/methods-pcaRes.R',sep=''))
                my_pc=try(mod_pca(data.na,method='mod_bpca',nPcs=2,verbose=F),silent=T)
                if ('try-error' %in% class(my_pc)){
                  #write.table(as.data.frame(my_pc[1]), paste(cd, "/simu_results/Error_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, "_2.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
                  #write.table(data.na, paste(cd, "/simu_results/DataError_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
                  BPCA=FALSE
                } else{
                  data.bPCA=completeObs(my_pc)
                }
                library(pcaMethods)
              } else{
                #write.table(data.na, paste(cd, "/simu_results/DataError_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
                BPCA=FALSE
              }
            } else{
              data.bPCA=completeObs(pc)
            }
            data.SampMin=SampMin(data.na)
          } else{
            ### NA col ###
            data.Tsm=SampMin(colsamp[[3]])
            # global min  % before %  imputation
            ris <- integer(ncol(data.Tsm)+length(colsamp[[2]]))
            if (length(ris)!=nsample*2) print('Error in bPCA imputation!')
            ris[colsamp[[2]]] <- ncol(data.Tsm)+1L
            ris[-colsamp[[2]]] <- seq_len(ncol(data.Tsm))
            
            pc=try(pca(colsamp[[3]],method='bpca',nPcs=2,verbose=F),silent = T)
            if ('try-error' %in% class(pc)){
              #write.table(as.data.frame(pc[1]), paste(cd, "/simu_results/Error_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
              pcsp=strsplit(pc[1],'computationally')[[1]]
              if (length(pcsp)==2){
                source(paste(cd,'/R/mod_bpca/pca.R',sep=''))
                source(paste(cd,'/R/mod_bpca/bpca.R',sep=''))
                source(paste(cd,'/R/mod_bpca/BPCA_initmodel.R',sep=''))
                source(paste(cd,'/R/mod_bpca/BPCA_dostep.R',sep=''))
                source(paste(cd,'/R/mod_bpca/checkData.R',sep=''))
                source(paste(cd,'/R/mod_bpca/prep.R',sep=''))
                source(paste(cd,'/R/mod_bpca/repmat.R',sep=''))
                source(paste(cd,'/R/mod_bpca/errorHierarchic.R',sep=''))
                source(paste(cd,'/R/mod_bpca/derrorHierarchic.R',sep=''))
                source(paste(cd,'/R/mod_bpca/AllClasses.R',sep=''))
                source(paste(cd,'/R/mod_bpca/AllGenerics.R',sep=''))
                source(paste(cd,'/R/mod_bpca/methods-pcaRes.R',sep=''))
                my_pc=try(mod_pca(colsamp[[3]],method='mod_bpca',nPcs=2,verbose=F),silent=T)
                if ('try-error' %in% class(my_pc)){
                  #write.table(as.data.frame(my_pc[1]), paste(cd, "/simu_results/Error_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, "_2.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
                  #write.table(data.na, paste(cd, "/simu_results/DataError_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
                  BPCA=FALSE
                } else{
                  data.temp=completeObs(my_pc)
                  data.bPCA=cbind(data.temp,min(data.Tsm))[,ris]
                }
                library(pcaMethods)
              } else{
                #write.table(data.na, paste(cd, "/simu_results/DataError_round", round, "_seed", seed500, "_", nsample, "_", fc, "_", zr, "_", mnar, ".txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=FALSE)
                BPCA=FALSE
              }
            } else{
              data.temp=completeObs(pc)
              data.bPCA=cbind(data.temp,min(data.Tsm))[,ris]
            }
            data.SampMin=cbind(data.Tsm,min(data.Tsm))[,ris]
          }
          # imputation done
          
          # split control/case sets
          control.SampMin=data.SampMin[,1:nsample]
          case.SampMin=data.SampMin[,(nsample+1):(nsample*2)]
          control.log.SampMin=log2(control.SampMin)
          case.log.SampMin=log2(case.SampMin)
          
          if (BPCA==TRUE){
            control.bPCA=data.bPCA[,1:nsample]
            case.bPCA=data.bPCA[,(nsample+1):(nsample*2)]
            control.log.bPCA=log2(control.bPCA)
            case.log.bPCA=log2(case.bPCA)
          }
          
          control.na=data.na[,1:nsample]
          case.na=data.na[,(nsample+1):(nsample*2)]
          control.zr=transfer_na2zero(control.na)
          case.zr=transfer_na2zero(case.na)
          # case.zr<-raw 0
          # case.na<-raw NA
          
          data.zr=as.matrix(cbind(control.zr,case.zr))
          rownames(data.zr) <- paste("feat_", 1:nrow(data.zr), sep = '')
          
          control.log.zr=log2(control.zr+1)
          case.log.zr=log2(case.zr+1)
          control.log.na=log2(control.na)
          case.log.na=log2(case.na)
          
          
          #two-sided t-test
          cat("----------T-test----------\n")
          pvalue.log.t.unpair=c()
          for (i in 1:nrow(case.log.na)) {
            p <- try(t.test(as.numeric(case.log.na[i,]), as.numeric(control.log.na[i,]),alternative='two.sided',paired = F,na.omit=T),silent=T)
            if ('try-error' %in% class(p)) {
              pvalue.log.t.unpair[i] <- NA
            } else pvalue.log.t.unpair[i] <- p$p.value
          }
          pvalue.log.t.unpair.BH=p.adjust(pvalue.log.t.unpair,method='BH')
          
          pvalue.t.sampmin=c()
          for (i in 1:nrow(case.log.SampMin)) {
            p <- try(t.test(as.numeric(case.log.SampMin[i,]), as.numeric(control.log.SampMin[i,]),alternative='two.sided',paired = F,na.omit=T),silent=T)
            if ('try-error' %in% class(p)) {
              pvalue.t.sampmin[i] <- NA
            } else pvalue.t.sampmin[i] <- p$p.value
          }
          pvalue.t.sampmin.BH=p.adjust(pvalue.t.sampmin,method='BH')

          if (BPCA==TRUE){
          pvalue.t.bPCA=c()
          for (i in 1:nrow(case.log.bPCA)) {
            p <- try(t.test(as.numeric(case.log.bPCA[i,]), as.numeric(control.log.bPCA[i,]),alternative='two.sided',paired = F,na.omit=T),silent=T)
            if ('try-error' %in% class(p)) {
              pvalue.t.bPCA[i] <- NA
            } else pvalue.t.bPCA[i] <- p$p.value
          }
          pvalue.t.bPCA.BH=p.adjust(pvalue.t.bPCA,method='BH')
          }else{
            pvalue.t.bPCA.BH=rep(NA,nrow(case.log.SampMin))
          }
          
          #wilcoxon test
          cat("----------wilcoxon----------\n")
          pvalue.wilcox.unpair <- c()
          for (i in 1:nrow(case.na)) {
            p <- try(wilcox.test(as.numeric(case.na[i,]), as.numeric(control.na[i,]),paired = F),silent=T)
            if ('try-error' %in% class(p)) {
              pvalue.wilcox.unpair[i] <- NA
            } else pvalue.wilcox.unpair[i] <- p$p.value
          }
          pvalue.wilcox.unpair.BH <- p.adjust(pvalue.wilcox.unpair, method = "BH")
          #names(pvalue.wilcox.BH) <- as.character(1:nrow(case.na))
          #hits.BH <- names(pvalue.wilcox.BH)[which(pvalue.wilcox < 0.05)]
          
          pvalue.wilcox.sampmin <- c()
          for (i in 1:nrow(case.SampMin)) {
            p <- try(wilcox.test(as.numeric(case.SampMin[i,]), as.numeric(control.SampMin[i,]),paired = F),silent=T)
            if ('try-error' %in% class(p)) {
              pvalue.wilcox.sampmin[i] <- NA
            } else pvalue.wilcox.sampmin[i] <- p$p.value
          }
          pvalue.wilcox.sampmin.BH <- p.adjust(pvalue.wilcox.sampmin, method = "BH")

          if (BPCA==TRUE){
          pvalue.wilcox.bPCA <- c()
          for (i in 1:nrow(case.bPCA)) {
            p <- try(wilcox.test(as.numeric(case.bPCA[i,]), as.numeric(control.bPCA[i,]),paired = F),silent=T)
            if ('try-error' %in% class(p)) {
              pvalue.wilcox.bPCA[i] <- NA
            } else pvalue.wilcox.bPCA[i] <- p$p.value
          }
          pvalue.wilcox.bPCA.BH <- p.adjust(pvalue.wilcox.bPCA, method = "BH")
          }else{
            pvalue.wilcox.bPCA.BH=rep(NA,nrow(case.SampMin))
          }
          
          # LIMMA: Moderated T-test (log2-conversion data)
          cat("----------Moderated T-test----------\n")
          pvalue.modT=try(modT(cbind(control.log.na,case.log.na),ncontrol=nsample,ncase=nsample),silent=T)
		      if ('try-error' %in% class(pvalue.modT)){
			      pvalue.modT=rep(NA,nrow(case.log.na))
			    }
          #r=decideTests(fit2,method='global',adjust.method = 'BH',p.value = .05,lfc=log2(1.2))
          #t=topTable(fit2,number=5000,adjust.method = 'BH')
          ## with lfc=log2(1.2), summary(r) <- 2504 up & 4 down, which is almost exactly, while only p<.05 <- 2571 DE.

          pvalue.modT.sampmin=try(modT(cbind(control.log.SampMin,case.log.SampMin),ncontrol=nsample,ncase=nsample),silent=T)
		      if ('try-error' %in% class(pvalue.modT.sampmin)){
			      pvalue.modT.sampmin=rep(NA,nrow(case.log.SampMin))
			    }

          if (BPCA==TRUE){
          pvalue.modT.bPCA=try(modT(cbind(control.log.bPCA,case.log.bPCA),ncontrol=nsample,ncase=nsample),silent=T)
		      if ('try-error' %in% class(pvalue.modT.bPCA)){
			      pvalue.modT.bPCA=rep(NA,nrow(case.log.bPCA))
			      }
          }else{
            pvalue.modT.bPCA=rep(NA,nrow(case.log.SampMin))
          }
          
          # 2-part t-test
          cat("----------2-part T-test----------\n")
          pvalue.2t=c()
          for (i in 1:nrow(case.log.zr)) {
            twoT <- try(TwoPart(as.numeric(c(control.log.zr[i,],case.log.zr[i,])), group, test="t.test", point.mass=0), silent=T)
            if ('try-error' %in% class(twoT)) {
              pvalue.2t[i] <- NA
            } else pvalue.2t[i] <- twoT$pvalue
          }
          pvalue.2t.BH=p.adjust(pvalue.2t,method='BH')
          
          # 2-part wilcoxon test
          cat("----------2-part Wilcoxon----------\n")
          pvalue.2wilcox=c()
          for (i in 1:nrow(case.zr)) {
            twoWILCOX <- try(TwoPart(as.numeric(c(control.zr[i,],case.zr[i,])), group, test="wilcoxon", point.mass=0), silent=T)
            if ('try-error' %in% class(twoWILCOX)) {
              pvalue.2wilcox[i] <- NA
            } else pvalue.2wilcox[i] <- twoWILCOX$pvalue
          }
          pvalue.2wilcox.BH=p.adjust(pvalue.2wilcox,method='BH')
          
          
          # SDA
          cat("----------Semi-parametric Differential Abundance----------\n")
          groupinfo=data.frame(grouping=matrix(rep(0:1,each=nsample),ncol=1))
          exampleSE=createSEFromMatrix(feature = data.zr,colData = groupinfo)
          
          # 1
          sda=try(SDA(sumExp = exampleSE, pi0=1),silent = T)
          pvalue.sda.BH=rep(NA,nrow(case.zr))
          #pvalue.sda.BH=p.adjust(sda$pv_2part,method='BH')
          if (!('try-error' %in% class(sda))){
            pvalue.sda.BH[sda_feat_ind(featstr=sda$feat.names)]=sda$qv_2part
          }
          if (length(pvalue.sda.BH)!=nrow(case.zr)){
            pvalue.sda.BH=rep(0,nrow(case.zr))
          }
          
          # 2
          my.sda=try(my.SDA(sumExp = exampleSE, pi0=1),silent = T)
          if ('try-error' %in% class(my.sda)){
            pvalue.mysda.BH=rep(NA,nrow(case.zr))
          } else{
            pvalue.mysda.BH=my.sda$qv_2part
            if (length(pvalue.mysda.BH)!=nrow(case.zr)){
              pvalue.mysda.BH=rep(0,nrow(case.zr))
            }
          }
          
          
          # ZIG & ZILN test
          cat("----------ZIG & ZILN 2-part test----------\n")
          #explore ZI-continuous, 'BH' (aka 'fdr')
          
          pvalue.zig=c()
          for (i in 1:nrow(case.zr)) {
            zig_test <- try(ziWaldTest(x=group,y=as.numeric(c(control.zr[i,],case.zr[i,])), method="zig", nit=50), silent=T)
            if ('try-error' %in% class(zig_test)) {
              pvalue.zig[i] <- NA
            } else pvalue.zig[i] <- zig_test$pvalue
          }
          pvalue.zig.fdr <- p.adjust(pvalue.zig, method = "BH")
          
          pvalue.ziln=c()
          for (i in 1:nrow(case.zr)) {
            ziln_test <- try(ziWaldTest(x=group,y=as.numeric(c(control.zr[i,],case.zr[i,])), method="ziln", nit=50), silent=T)
            if ('try-error' %in% class(ziln_test)) {
              pvalue.ziln[i] <- NA
            } else pvalue.ziln[i] <- ziln_test$pvalue
          }
          pvalue.ziln.fdr <- p.adjust(pvalue.ziln, method = "BH")
          
          
          cat("----------ZIG & ZILN Mean-based test----------\n")
          #explore ZI-continuous, 'BH' (aka 'fdr')
          
          pvalue.zigDM=c()
          for (i in 1:nrow(case.zr)) {
            zig_test <- try(ziDMtest(x=group,y=as.numeric(c(control.zr[i,],case.zr[i,])), method="zig", nit=50), silent=T)
            if ('try-error' %in% class(zig_test)) {
              pvalue.zigDM[i] <- NA
            } else pvalue.zigDM[i] <- zig_test
          }
          pvalue.zigDM.fdr <- p.adjust(pvalue.zigDM, method = "BH")
          
          pvalue.zilnDM=c()
          for (i in 1:nrow(case.zr)) {
            ziln_test <- try(ziDMtest(x=group,y=as.numeric(c(control.zr[i,],case.zr[i,])), method="ziln", nit=50), silent=T)
            if ('try-error' %in% class(ziln_test)) {
              pvalue.zilnDM[i] <- NA
            } else pvalue.zilnDM[i] <- ziln_test
          }
          pvalue.zilnDM.fdr <- p.adjust(pvalue.zilnDM, method = "BH")
		
          # MSstatsPTM
          cat("----------MSstatsPTM----------\n")
          pvalue.msstats=msstats(data.na,impute=F)
          pvalue.msstats.BH=p.adjust(pvalue.msstats,method='BH')
          
          cat("----------MSstatsPTM_AFT----------\n")
          pvalue.msstats.AFT=msstats(data.na,impute=T)
          pvalue.msstats.AFT.BH=p.adjust(pvalue.msstats.AFT,method='BH')
          
		
          output.data <- data.frame(pvalue.log.t.unpair.BH, pvalue.t.sampmin.BH, pvalue.t.bPCA.BH, pvalue.wilcox.unpair.BH, pvalue.wilcox.sampmin.BH, pvalue.wilcox.bPCA.BH, pvalue.modT, pvalue.modT.sampmin, pvalue.modT.bPCA, pvalue.2t.BH, pvalue.2wilcox.BH, pvalue.sda.BH, pvalue.mysda.BH, pvalue.zig.fdr, pvalue.ziln.fdr, pvalue.zigDM.fdr, pvalue.zilnDM.fdr, pvalue.msstats.BH, pvalue.msstats.AFT.BH, truth.na)
          colnames(output.data) <- c("T-test", 'T-SampMin', 'T-bPCA', "Wilcoxon", 'Wilcoxon-SampMin', 'Wilcoxon-bPCA', "ModT", 'ModT-SampMin', 'ModT-bPCA', 'twoT', 'twoWilcox', 'SDA_robust', 'SDA', 'ZIG_2p', 'ZILN_2p', 'ZIG_DM', 'ZILN_DM', "MSstatsPTM", 'MSstatsPTM-AFT', "TrueDE")
          write.table(output.data, paste(cd, "/simu_results/round", round, "/Simu_", nsample, "_", fc, "_", zr, "_", mnar, "_for_", round, "_times.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, append=FALSE)
        }
      }
    }
  }
}






