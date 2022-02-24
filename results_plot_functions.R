library(ggplot2)
library(gridExtra)

labeldict=data.frame(c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p','MSstatsPTM','MSstatsPTM-AFT'),
                     c('T-test','T-bPCA','T-SampMin','Wilcox','Wilcox-bPCA','Wilcox-SampMin','ModT','ModT-bPCA','ModT-SampMin','2part-T','2part-Wilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2part','ZILN_2part','MSstatsPTM','MSstatsPTM-AFT'))
colnames(labeldict)=c('rawname','legend')

## color
# v2: Figure 2 & 3
colorpalette=c('#FFC1C1','#FF6347','#CD2626',  #'#FF6A6A','#FF3030','#8B1A1A',
               '#90EE90','#32CD32','#228B22',
               '#ADD8E6','#00BFFF','#4876FF',  #'#0000FF','#0000CD',
               '#4eb7d9',"#3E5CC5","#7149D7","#FF99FF","#CC99FF","#CCCCCC","#8C8C8C",
               '#FFBBFF','#E066FF')
# raw (Fig 4)
colorpalette=c("#DB2F20","#9F2218",'#EE3768',"#E68200","#E6B600","#E6EB00","#A5CF47","#00B159","#1F868F",'#4eb7d9',"#3E5CC5","#7149D7","#FF99FF","#CC99FF","#CCCCCC","#8C8C8C",'#FFD700','#CDAD00')
# Fig 6 & 5 & Fig SI
colorpalette=c("#DB2F20","#9F2218",'#EE3768',"#E68200","#E6B600","#E6EB00","#A5CF47","#00B159","#1F868F",'#4eb7d9',"#3E5CC5","#7149D7","#FF99FF","#CC99FF","#CCCCCC","#8C8C8C","#FF99FF","#CC99FF")
# Figs SI
#colorpalette=c("#DB2F20","#9F2218",'#EE3768',"#E68200","#E6B600","#E6EB00","#A5CF47","#00B159","#1F868F",'#4eb7d9',"#3E5CC5","#7149D7","#FF99FF","#CC99FF","#CCCCCC","#8C8C8C","#FFDAB9","#CC99FF")

# distinct color palette
colpal0=data.frame(c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p','MSstatsPTM','MSstatsPTM-AFT'),colorpalette)
colnames(colpal0)=c('model','colorpalette')

colpal=data.frame(c('T-test','T-bPCA','T-SM','Wilcox','W-bPCA','W-SM','ModT','ModT-bPCA','ModT-SM','2part-T','2part-W','SDA','ZIG_DM','ZILN_DM','ZIG_2part','ZILN_2part','MSstatsPTM','MSstatsPTM-AFT'),colorpalette)
colnames(colpal)=c('model','colorpalette')


####### ZI continuous models --- Box plot function ######
## The box-plot function is modified from the codes in Ding et al (2020)'s excellent work. ##
## Refer to: https://github.com/hmsch/proteomics-simulations ##

plotbox=function(data=result.auc,n=5,zr=.3,mnar=0,fc=NULL,evaluation='pAUROC',modellist=c('T-test','T-SampMin','T-bPCA','Wilcoxon','Wilcoxon-SampMin','Wilcoxon-bPCA','ModT','ModT-SampMin','ModT-bPCA','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p'),textsizeL=20,textsizeM=17,ylimit=1.0){
  if(evaluation %in% c('pAUROC','AUROC')){
    palettetmp=colpal0$colorpalette[colpal0$model %in% modellist]
    legendtmp=labeldict$legend[labeldict$rawname %in% modellist]
    if(is.null(fc)){
      # xlab: FC --- fig_version 1
      df=result.auc[(result.auc[,1]==n)&(result.auc[,3]==zr)&(result.auc[,4]==mnar)&(result.auc$model %in% modellist),]
      ggplot(df, aes(x = as.factor(fc), fill=model, y = get(evaluation)))+theme_bw()+
        geom_boxplot(outlier.shape=1, outlier.size=0.1, lwd=0.3) + 
        # Add lines between each group
        geom_vline(xintercept=seq(1.5, length(unique(df$fc))-0.5, 1),
                   lwd=0.5, color="grey") +
        # what confusing command below ?!
        #scale_y_continuous(limits = quantile(df$pAUROC, c(0.01, 1.0),na.rm=T)) +
        scale_fill_manual(values = palettetmp,breaks = modellist,labels=legendtmp) +
        scale_y_continuous(limits = c(0.5,1)) +
        theme( 
          panel.grid=element_blank(),
          # remove the vertical grid lines
          panel.grid.major.x = element_blank() ,
          # explicitly set the horizontal lines (or they will disappear too)
          panel.grid.major.y = element_line(size=.3, color="grey"),
          # Set legend size
          legend.key.size = unit(1.5, "lines"),
          #legend.position = "top",
          axis.text = element_text(size = textsizeM),
          axis.title = element_text(size = textsizeL),
          legend.text = element_text(size = textsizeM)
        ) +
        coord_fixed(ratio=9) +
        labs(
          #title='pAUROC of all methods', 
          x='Fold change', y=evaluation, fill="")
    }else if(is.null(mnar)){
      # xlab: MNAR --- fig_version 2
      df=result.auc[(result.auc[,4] %in% c(0,.2,.4,.6,.8,1))&(result.auc[,1]==n)&(result.auc[,2]==fc)&(result.auc[,3]==zr)&(result.auc$model %in% modellist),]
      fig<-
      ggplot(df, aes(x = as.factor(MNAR), y = get(evaluation), fill=model, color=model))+theme_bw()+
        geom_boxplot(outlier.shape=1, outlier.size=0.1, lwd=0.3,position = position_dodge(.9),width=.8)+
        # position .9 + width .8 / position .8 +width .7
        stat_boxplot(geom='errorbar',position = position_dodge(.9),width=.3)+
        geom_vline(xintercept=seq(1.5, length(unique(df$MNAR))-0.5, 1),
                   lwd=0.5, color="grey") +
        scale_fill_manual(values=palettetmp,breaks=modellist,labels=legendtmp) +
        scale_color_manual(values=palettetmp,breaks=NULL,labels=NULL) +
        ##stat_summary(fun.data=f,geom='crossbar',color='black',width=0.1,aes(x = as.factor(MNAR))) +
        scale_y_continuous(limits = c(0.5, ylimit),expand=c(0,0)) +
        #scale_x_discrete(breaks=c(0,.2,.4,.6,.8,1))+
        theme( 
          panel.grid=element_blank(),
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line(size=.3, color="grey"),
          legend.key.size = unit(1.5, "lines"),
          axis.text = element_text(size = textsizeM),
          axis.title = element_text(size = textsizeL),
          legend.text = element_text(size = textsizeM),
          legend.position = 'bottom' #'none'#
        ) +
        coord_fixed(ratio = 6/(2*ylimit-1))+ #9/(2*ylimit-1))+
        labs(x='MNAR ratio', y=evaluation, fill="")
      #dat <- ggplot_build(fig)$data[[1]]
      #a=data.frame(rep(modellist,dim(dat)[1]/length(modellist)));colnames(a)=c('model')
      #dat=cbind(dat,a)
      #return(fig + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
      #                               y=middle, yend=middle), colour="black", size=0.5))
      return(fig)
    }else if(is.null(zr)){
      # xlab: zero ratio
      df=result.auc[(result.auc[,1]==n)&(result.auc[,2]==fc)&(result.auc[,4]==mnar)&(result.auc$model %in% modellist),]
      fig<-
        ggplot(df, aes(x = as.factor(zr), y = get(evaluation), fill=model, color=model))+theme_bw()+
        geom_boxplot(outlier.shape=1, outlier.size=0.1, lwd=0.3) + 
        geom_vline(xintercept=seq(1.5, length(unique(df$zr))-0.5, 1),
                   lwd=0.5, color="grey") +
        scale_fill_manual(values = palettetmp,breaks = modellist,labels=legendtmp) +
        scale_color_manual(values = palettetmp,breaks=NULL,labels=NULL) +
        scale_y_continuous(limits = c(0.5, 1.0)) +
        theme( 
          panel.grid=element_blank(),
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line(size=.3, color="grey"),
          legend.key.size = unit(1.5, "lines"),
          axis.text = element_text(size = textsizeM),
          axis.title = element_text(size = textsizeL),
          legend.text = element_text(size = textsizeM)
        ) +
        coord_fixed(ratio = 9) +
        labs(x='Zero ratio', y=evaluation, fill="")
      dat <- ggplot_build(fig)$data[[1]]
      a=data.frame(rep(modellist,dim(dat)[1]/length(modellist)));colnames(a)=c('model')
      dat=cbind(dat,a)
      return(fig + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                              y=middle, yend=middle), colour="black", size=0.5))
    }else if(is.null(n)){
      # xlab: Nsample
      df=result.auc[(result.auc[,2]==fc)&(result.auc[,3]==zr)&(result.auc[,4]==mnar)&(result.auc$model %in% modellist),]
      fig<-
        ggplot(df, aes(x = as.factor(Nsample), y = get(evaluation), fill=model, color=model))+theme_bw()+
        geom_boxplot(outlier.shape=1, outlier.size=0.1, lwd=0.3) + 
        geom_vline(xintercept=seq(1.5, length(unique(df$Nsample))-0.5, 1),
                   lwd=0.5, color="grey") +
        scale_fill_manual(values = palettetmp,breaks = modellist,labels=legendtmp) +
        scale_color_manual(values = palettetmp,breaks=NULL,labels=NULL) +
        scale_y_continuous(limits = c(0.5, 1.0)) +
        theme( 
          panel.grid=element_blank(),
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_line(size=.3, color="grey"),
          legend.key.size = unit(1.5, "lines"),
          axis.text = element_text(size = textsizeM),
          axis.title = element_text(size = textsizeL),
          legend.text = element_text(size = textsizeM)
        ) +
        coord_fixed(ratio = 9) +
        labs(x='Sample size', y=evaluation, fill="")
      dat <- ggplot_build(fig)$data[[1]]
      a=data.frame(rep(modellist,dim(dat)[1]/length(modellist)));colnames(a)=c('model')
      dat=cbind(dat,a)
      return(fig + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                              y=middle, yend=middle), colour="black", size=0.5))
    }
  }else print('Error: Eval index should be pAUROC or AUROC')
}
# test
#plotbox(data=result.auc, n=5, zr=.3, mnar=0, fc=NULL, evaluation='pAUROC')
#plotbox(data=result.auc, n=5, zr=.3, mnar=0, fc=NULL, evaluation='AUROC')

plotbox2=function(data=result.auc,n=c(10,30),zr=.3,mnar=c(0,.2,.4,.6,.8,1),fc=2,evaluation='pAUROC',modellist=c('T-test','T-SampMin','T-bPCA','Wilcoxon','Wilcoxon-SampMin','Wilcoxon-bPCA','ModT','ModT-SampMin','ModT-bPCA','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p'),textsizeL=17,textsizeM=14,ylimit=1.0,No.fig='fig2'){
  if(evaluation %in% c('pAUROC','AUROC')){
    palettetmp=colpal0$colorpalette[colpal0$model %in% modellist]
    legendtmp=labeldict$legend[labeldict$rawname %in% modellist]
    df=result.auc[(result.auc[,4] %in% mnar)&(result.auc[,1] %in% n)&(result.auc[,2] %in% fc)&(result.auc[,3] %in% zr)&(result.auc$model %in% modellist),]
    if (No.fig %in% c('fig2','fig3cd')){
      df['SampleSize']=as.factor(paste0('N=',as.character(df$Nsample)))
    }else if (No.fig %in% c('fig3ab','figsi')){
      df['MissingRatio']=as.factor(paste0('Missing % = ',as.character(df$zr*100),'%'))
    }else if (No.fig=='fig6'){
      df=result.auc[(result.auc$model %in% modellist)&(result.auc[,4] %in% mnar),]
      df=rbind(df[(df$Nsample==5)&(df$zr==.3)&(df$fc==1.5),],
               df[(df$Nsample==50)&(df$zr==.9)&(df$fc==2),],
               df[(df$Nsample==80)&(df$zr==.7)&(df$fc==1.5),])
      df['scenario']=as.factor(ifelse(df$Nsample==5,'Scenario 1',ifelse(df$Nsample==50,'Scenario 2','Scenario 3')))
    }
    
    fig<-ggplot(df, aes(y = get(evaluation), fill=model, color=model))+theme_bw()
    if (No.fig=='fig2'){
      fig=fig+facet_grid(SampleSize~MNAR)+
        labs(x='MNAR ratio', y=evaluation, fill="")
    }else if (No.fig=='fig3ab'){
      fig=fig+facet_grid(MissingRatio~Nsample)+
        labs(x='Sample size', y=evaluation, fill="")
    }else if (No.fig=='fig3cd'){
      fig=fig+facet_grid(SampleSize~zr)+
        labs(x='Missing ratio', y=evaluation, fill="")
    }else if (No.fig=='fig4'){
      fig=fig+facet_wrap(~MNAR,nrow=1)+
        labs(x='MNAR ratio', y=evaluation, fill="")
    }else if (No.fig=='fig6'){
      fig=fig+facet_grid(scenario~MNAR)+#,scales='free_y')+
        labs(x='MNAR ratio', y=evaluation, fill="")
    }else if (No.fig=='figsi'){
      fig=fig+facet_grid(MissingRatio~MNAR)+
        labs(x='MNAR ratio', y=evaluation, fill="")
    }
    
    if (No.fig %in% c('fig6','figsi')){
      fig=fig+geom_vline(xintercept=c(-0.412,-0.337,-0.262,-0.188,-0.113,-0.04,0.04,0.113,0.188,0.262,0.337,0.412),lwd=0.5, color="gray91")
    } else {
      fig=fig+geom_vline(xintercept=c(-0.41,-0.327,-0.246,-0.164,-0.082,0,0.082,0.164,0.246,0.327,0.408),lwd=0.5, color="gray91")
        #seq(1.5, length(unique(df$MNAR))-0.5, 1),lwd=0.5, color="grey") +
    }
    fig<-fig+
      geom_boxplot(outlier.shape=1, outlier.size=0.1, lwd=0.3,position = position_dodge(.9),width=.8)+
      stat_boxplot(geom='errorbar',position = position_dodge(.9),width=.3)+
      scale_fill_manual(values=palettetmp,breaks=modellist,labels=legendtmp) +
      scale_color_manual(values=palettetmp,breaks=NULL,labels=NULL) +
      scale_y_continuous(limits = c(0.5, ylimit),expand=c(0,0)) +
      scale_x_continuous(labels=NULL)+
      theme( 
        panel.grid=element_blank(),
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(size=.3, color="gray91"),
        legend.key.size = unit(1.5, "lines"),
        axis.text = element_text(size = textsizeM-6),
        axis.title = element_text(size = textsizeL),
        legend.text = element_text(size = textsizeM),
        legend.position = 'none', #'bottom'
        axis.ticks.x=element_blank(),
        panel.spacing.y=unit(0.5, "lines")
      ) +
      coord_fixed(ratio = 6/(2*ylimit-1)) #9/(2*ylimit-1))+
    return(fig)
  }else print('Error: Eval index should be pAUROC or AUROC')
}


method_mean=function(x){
  temp=matrix(NA,0,10)
  method=names(table(x$model))
  for (i in 1:length(method)){
    result=x[x$model==method[i],]
    n=ncol(result)
    auc=colMeans(result[,(n-4):n],na.rm = T)
    simu.para=result[1,c(1:4,6)]
    mm=cbind(simu.para,t(as.data.frame(auc)),row.names=NULL)
    temp=rbind(temp,mm)
  }
  return(temp)
}

###### Best Method Result ######
## select best methods
best_mean=function(x,idx='pAUROC'){
  res=data.frame()
  if(idx=='pAUROC'){
    idx_col=7
  }else if(idx=='AUROC'){
    idx_col=6
  }else{
    print('Evaluation index should be "AUROC" or "pAUROC"!')
  }
  for (n in c(5,10,30,50,80,100)){
    for (ZR in c(.3,.4,.5,.6,.7,.8,.9)){
      for (mnar in c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)){
        for (FC in c(1.2,1.5,1.8,2,3,5,10)) {
          best=x[(x$Nsample==n)&(x$zr==ZR)&(x$MNAR==mnar)&(x$fc==FC),]
          bes=best[best[,idx_col]==max(best[,idx_col],na.rm=T),]
          if(nrow(res)==0){
            res=bes[,c(1:5,idx_col)]
          }else{
            res=rbind(res,bes[,c(1:5,idx_col)])
          }
        }
      }
    }
  }
  return(res)
}

mean_mean=function(x,idx='pAUROC',mnar_high=1,mnar_low=0){
  res=data.frame()
  if(idx=='pAUROC'){
    idx_col=7
  }else if(idx=='AUROC'){
    idx_col=6
  }else{
    print('Evaluation index should be "AUROC" or "pAUROC"!')
  }
  for (n in c(5,10,30,50,80,100)){
    for (ZR in c(.3,.4,.5,.6,.7,.8,.9)){
      for (FC in c(1.2,1.5,1.8,2,3,5,10)) {
        for (m in c('T-test','T-bPCA','T-SampMin','Wilcoxon','Wilcoxon-bPCA','Wilcoxon-SampMin','ModT','ModT-bPCA','ModT-SampMin','twoT','twoWilcox','SDA','ZIG_DM','ZILN_DM','ZIG_2p','ZILN_2p','MSstatsPTM','MSstatsPTM-AFT')){
          best=x[(x$Nsample==n)&(x$zr==ZR)&(x$model==m)&(x$fc==FC),]
          bes=best[1,c(1:3,5,idx_col)]
          bes[,5]=mean(best[,idx_col],na.rm=T)
          if(nrow(res)==0){
            res=bes
          }else{
            res=rbind(res,bes)
          }
        }
      }
    }
  }
  return(res)
}

best_best=function(x){
  res=data.frame()
  for (n in c(5,10,30,50,80,100)){
    for (ZR in c(.3,.4,.5,.6,.7,.8,.9)){
      for (FC in c(1.2,1.5,1.8,2,3,5,10)) {
        best=x[(x$Nsample==n)&(x$zr==ZR)&(x$fc==FC),]
        bes=best[best$E_pAUROC==max(best$E_pAUROC,na.rm=T),]
        if(nrow(res)==0){
          res=bes
        }else{
          res=rbind(res,bes)
        }
      }
    }
  }
  return(res)
}


####### ZI continuous models --- best method function ######

resplot3<-function(dataset,type='pAUROC',samplesize=5,zeroratio=.3){
  bbb=dataset
  palettetmp=colpal$colorpalette[colpal$model %in% bbb$method]
  pcor<-ggplot(bbb,aes(Var2,Var1))+geom_tile(aes(fill=method),colour="white")
  ppplot<-pcor+
    scale_fill_manual(values=palettetmp)+
    theme(axis.text.x = element_text(vjust = 0.5,hjust = 0.5,angle = 0))+
    coord_fixed(ratio = 1)+
    theme(axis.text = element_text(size = 15,family = "ARL"))+
    #theme(plot.margin = unit(c(0.1,0,0,0),unit ="mm" ))+
    labs(x="Fold Change",y="MNAR Ratio",title=paste0('N=',samplesize,', Missing%=',zeroratio),family="ARL")+
    theme(axis.title = element_text(size = 18))+
    theme(plot.title = element_text(size = 18,hjust = 0.5,family = "ARL",face = 'bold'))+
    theme(legend.key.width = unit(6,'mm'),legend.key.height = unit(5,'mm'))+
    theme(legend.text = element_text(size=15))+
    theme(legend.title = element_text(size = 15))+
    geom_text(aes(Var2,Var1,label=format(value,digits=2,nsmall=2)),color='black',size=5)
  svg(filename = paste(type,'_n',samplesize,'_zr_',zeroratio,'_value',".svg",sep = ""))
  plot(ppplot)
  dev.off()
}

resplot3b<-function(dataset,type='E_pAUROC',samplesize=5){
  bbb=dataset
  palettetmp=colpal$colorpalette[colpal$model %in% bbb$method]
  pcor<-ggplot(bbb,aes(Var2,Var1))+geom_tile(aes(fill=method),colour="white")
  ppplot<-pcor+
    scale_fill_manual(values=palettetmp)+
    theme(axis.text.x = element_text(vjust = 0.5,hjust = 0.5,angle = 0))+
    coord_fixed(ratio = 1)+
    theme(axis.text = element_text(size = 17,family = "ARL"))+
    #theme(plot.margin = unit(c(0.1,0,0,0),unit ="mm" ))+
    labs(x="Fold Change",y="Missing Ratio",title=paste0('N=',samplesize),family="ARL")+
    theme(axis.title = element_text(size = 24))+
    theme(plot.title = element_text(size = 24,hjust = 0.5,family = "ARL",face = 'bold'))+
    theme(legend.key.width = unit(6,'mm'),legend.key.height = unit(5,'mm'))+
    theme(legend.text = element_text(size=17))+
    theme(legend.title = element_text(size = 17))+
    theme(legend.position = 'none')+
    geom_text(aes(Var2,Var1,label=format(value,digits=2,nsmall=2)),color='black',size=7)
  svg(filename = paste(type,'_n',samplesize,'_value.svg',sep = ""))
  plot(ppplot)
  dev.off()
}



