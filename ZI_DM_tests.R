# ZIG & ZILN models test
## The codes in this script are from Mills (2013)'s PhD thesis. ##
## Refer to: https://doi.org/10.17077/etd.7v3bafbd ##

### model estimation code ###
## Fisher scoring algorithms
# Step-Halving
step.halve=function(d,loglik.old,loglik.new,x,y,bhat,bdiff){
  eps=1e-04
  if(loglik.old<=loglik.new+eps){
    step.ok=1
    return(c(step.ok,bhat,loglik.new))
  } else{
    bhat=bhat-bdiff
    for(step in 1:10){
      bdiff=bdiff/2
      bhat=bhat+bdiff
      if(d=='gamma'){
        muhat=exp(x%*%bhat)
        loglik.new=sum(-log(muhat)-y/muhat)
      }
      if(d=='bern'){
        niln=log(1/(1+exp(x %*% bhat)))
        niln[is.na(niln)]=0
        loglik.new=crossprod(x %*% bhat,y)-sum(niln)
      }
      if(loglik.new+eps>=loglik.old){
        step.ok=2
        return(c(step.ok,bhat,loglik.new))
      }
    }
    step.ok=0
    c(step.ok,bhat,loglik.new)
  }
}
# Fisher Scoring for binomial
glm.bin=function(x,y,nit,tol=0.0000001){
  if(sum(y)==0){
    cat("\n\nAll the Ys are zeros!\n\n")
  }else y.b=ifelse(y==0,0,1)
  # Initial estimate of intercept
  b0.b=log(mean(y.b)/(1-mean(y.b)))
  # Setup for iteration routine
  bdiff.b=1
  if(length(x)==length(y)){
    nb=1
    bhat.b=t(b0.b)
  }else{
    nb=dim(x)[2]
    bhat.b=c(b0.b,rep(0,nb-1))
  }
  it=0
  space=""
  niln=log(1/(1+exp(x%*%bhat.b)))
  niln[is.na(niln)]=0
  loglik.null.b=crossprod(x%*%bhat.b,y.b)+sum(niln)
  loglik.old.b=loglik.null.b
  
  # Finding MLEs and standard errors
  while (it <= nit && max(abs(bdiff.b)) > tol){
    it=it+1
    p=exp(x%*%bhat.b)/(1+exp(x%*%bhat.b))
    inf=t(x)%*%(diag(length(p))*as.vector(p*(1-p)))%*% x
    resid.b=y.b-p
    score.b=t(x)%*%(resid.b)
    bhat.b.old=bhat.b
    bhat.b=bhat.b+solve(inf)%*%score.b
    p=exp(x%*%bhat.b)/(1+exp(x%*%bhat.b))
    niln=log(1/(1+exp(x%*%bhat.b)))
    niln[is.na(niln)]=0
    loglik.b=crossprod(x%*%bhat.b,y.b)+sum(niln)
    bdiff.b=bhat.b-bhat.b.old
    
    st.res=step.halve('bern',loglik.old.b,loglik.b,x,y.b,bhat.b,bdiff.b)
    if (st.res[1]==2){
      bhat.b=st.res[2:(nb+1)]
      loglik.b=st.res[(nb+2)]
    }
    if (st.res[1]==0){
      cat("Step halving failed to find better estimates \n")
      break
    }
    loglik.old.b=loglik.b
  }
  # Calculate results
  if (it>nit){
    cat("Failed to converge in",nit,"steps\n")
    iconvg=0
  } else if(st.res[1]==0){
    cat("Step halving failed\n")
  } else{
    niln=log(1/(1+exp(x%*%bhat.b)))
    niln[is.na(niln)]=0
    # log-likelihood at mle
    loglik.b=crossprod(x%*%bhat.b,y.b)+sum(niln)
    inf=t(x)%*%(diag(length(p))*as.vector(p*(1-p)))%*%x
    resid.b=y.b-p
    score.b=t(x)%*%resid.b
    infinv=solve(inf)
    sd1.b=sqrt(diag(infinv))
  }
  return(list(llb=loglik.b,b.b=bhat.b,sd1.b=sd1.b,score.b=score.b,invinf.b=infinv))
}

# Fisher Scoring for gamma
glm.gamma=function(x,y,nit,tol=0.00000001){
  y.c=y[y!=0]
  # Initial estimate of intercept
  b0=log(mean(y.c))
  # Saturated log-likelihood/v (scaled; c(y,v) dropped)
  loglik.sat=sum(-log(y.c)-1)
  # Setup for iteration routine
  bdiff=1
  nc=dim(as.matrix(x))[2]
  #x.c contains the rows in x corresponding to non-zero y's
  x.c=matrix(x[y!=0],nrow=length(y.c),ncol=nc)
  bhat=c(b0,rep(0,nc-1))
  it=0
  space=""
  muhat=exp(x.c%*%bhat)
  # Null log-likelihood/v (scaled version; dropping c(y,v))
  # Note: scaled v cancels out in maximizaion; c(y,v) doesn't impact maximization
  loglik.null=sum(-log(muhat)-y.c/muhat)
  deviance.null=2*(loglik.sat-loglik.null)
  loglik.old=loglik.null
  
  #Finding MLEs and standard errors
  while (it <= nit && max(abs(bdiff)) > tol){
    it=it+1
    muhat=exp(x.c%*%bhat)
    resid=y.c-muhat
    bhat.old=bhat
    inf=crossprod(x.c,x.c)
    score=t(x.c)%*%(resid/muhat)
    bhat=bhat.old+solve(inf)%*%score
    muhat.old=muhat
    muhat=exp(x.c%*%bhat)
    loglik=sum(-log(muhat)-y.c/muhat)
    bdiff=bhat-bhat.old
    # Step-Halving
    st.res=step.halve('gamma',loglik.old,loglik,x.c,y.c,bhat,bdiff)
    if (st.res[1]==2){
      bhat=st.res[2:(nc+1)]
      loglik=st.res[(nc+2)]
    }
    if (st.res[1]==0){
      cat("Step halving failed to find better estimates\n")
      break
    }
    loglik.old=loglik
  }
  # Calculate results
  if (it > nit){
    cat("Failed to converge in",nit,"steps\n")
    iconvg=0
  }else if (st.res[1]==0){
    cat("Step halving failed\n")
  }else{
    iconvg=1
    muhat=exp(x.c%*%bhat)
    loglik=sum(-log(muhat)-y.c/muhat)
    deviance=2*(loglik.sat-loglik)
    vinv.mle=deviance/(length(y.c)-nc)
    covb.mle=vinv.mle*solve(crossprod(x.c,x.c))
    
    sdb.mle=sqrt(diag(covb.mle))
    resid.c=y.c-muhat
    score=t(x.c)%*%(resid.c/muhat)
  }
  return(list(llc=loglik,b.c=bhat,sd1.c.mle=sdb.mle,disp.mle=vinv.mle,invinf.c=covb.mle,score.c=score,dev=deviance,ncont=length(y.c),resid.c=resid.c))
}

# Putting the two parts together
zig.glm=function(x,y,nit,tol=0.00000001){
  glmb=glm.bin(x,y,nit,tol)
  glmg=glm.gamma(x,y,nit,tol)
  llzig=glmb$llb+glmg$llc
  py=exp(x%*%glmb$b.b)/(1+exp(x%*%glmb$b.b))
  muhat=exp(x%*%glmg$b.c)
  pred.y=py*muhat
  resid=y-pred.y
  return(list(llb=glmb$llb,llc=glmg$llc,llzic=llzig,b.b=glmb$b.b,b.c=glmg$b.c,sd1.b=glmb$sd1.b,sd1.c.mle=glmg$sd1.c.mle,disp.mle=glmg$disp.mle,invinf.b=glmb$invinf.b,invinf.c=glmg$invinf.c,score.b=glmb$score.b,score.c=glmg$score.c,dev.c=glmg$dev,ncont=glmg$ncont,resid.c=glmg$resid.c,pred.y=pred.y,resid=resid))
}

## Log-Normal Regression Estimation
lm.lnorm=function(x,y){
  if (sum(y)==0) cat("\n\nAll the Ys are zeros!\n\n")
  y.c=y[y!=0]
  if (length(x)==length(y)){
    x.c=matrix(x[y!=0],nrow=length(y.c),ncol=1)
  } else{
    x.c=matrix(x[y!=0],nrow=length(y.c),ncol=dim(x)[2])
  }
  nc=dim(x.c)[2]
  invinf.unsc=solve(crossprod(x.c,x.c))
  bhat=invinf.unsc%*%crossprod(x.c,log(y.c))
  resid.c=log(y.c)-x.c%*%bhat
  sigma2=crossprod(resid.c)/(length(y.c)-nc)
  #sdb.c=sqrt(sigma2*diag(invinf.unsc))
  sdb.c=sqrt(sigma2[1]*diag(invinf.unsc))
  llc=crossprod(resid.c)
  return(list(llc=llc,b.c=bhat,sdb.c=sdb.c,sigma2=sigma2,invinf.c=as.numeric(sigma2)*invinf.unsc,resid.c=resid.c))
}

# Putting ZILN pieces together
ziln.glm=function(x,y,nit){
  glmb=glm.bin(x,y,nit)
  glmln=lm.lnorm(x,y)
  llziln=glmb$llb+glmln$llc
  py=exp(x%*%glmb$b.b)/(1+exp(x%*%glmb$b.b))
  muhat=exp(x%*%glmln$b.c)
  pred.y=py*muhat
  resid=y-pred.y
  return(list(llb=glmb$llb,llc=glmln$llc,llzic=llziln,b.b=glmb$b.b,b.c=glmln$b.c,sd1.b=glmb$sd1.b,sd1.c=glmln$sdb.c,disp=glmln$sigma2,invinf.b=glmb$invinf.b,invinf.c=glmln$invinf.c,resid.c=glmln$resid.c,pred.y=pred.y,resid=resid))
}

### Mean-based 2 Group Comparison without Covariate Adjustment
## Mean-based Tests for Zero Inflated Gamma
#--- IRM is lRM ---#
zig.uni=function(mod,n0,n1){
  n=n0+n1
  p1=exp(mod$b.b[1]+mod$b.b[2])/(1+exp(mod$b.b[1]+mod$b.b[2]))
  p0=exp(mod$b.b[1])/(1+exp(mod$b.b[1]))
  n0c=n0*p0
  n1c=n1*p1
  nc=n0c+n1c
  M0=exp(mod$b.c[1])*p0
  M1=exp(mod$b.c[1]+mod$b.c[2])*p1
  
  DM=exp(mod$b.c[1]+mod$b.c[2])*p1-exp(mod$b.c[1])*p0
  VDM=exp(2*mod$b.c[1])*p0*(1-p0)/n0+exp(2*mod$b.c[1]+2*mod$b.c[2])*p1*(1-p1)/n1+p0^2*exp(2*mod$b.c[1])*mod$disp.mle/n0c+p1^2*exp(2*(mod$b.c[1]+mod$b.c[2]))*mod$disp.mle/n1c
  # Note: p0^2/n0c=p0/n0
  DMtest=DM^2/VDM
  DM.p=1-pchisq(DMtest,1)
  
  IRM=mod$b.c[2]+mod$b.b[2]+log(1-p1)-log(1-p0)
  VIRM=mod$disp.mle*(nc/(n0c*n1c))+(1-p0)/(n0*p0)+(1-p1)/(n1*p1)
  # Note: nc/(n0c*n1c)=(n0*p0+n1*p1)/(n0*p0*n1*p1)=1/(n0*p0)+1/(n1*p1)
  IRMtest=IRM^2/VIRM
  IRM.p=1-pchisq(IRMtest,1)
  return(list(M0=M0,M1=M1,DM=DM,VDM=VDM,DMtest=DMtest,DM.p=DM.p,IRM=IRM,VIRM=VIRM,IRMtest=IRMtest,IRM.p=IRM.p))
}

wilcoxon=function(yw,xw){
  r=rank(yw)
  s2r=var(r)
  r0=mean(r[xw==0])
  r1=mean(r[xw==1])
  Wlcx=(r0-r1)^2/(s2r*(1/length(r[xw==0])+1/length(r[xw==1])))
  return(Wlcx)
}

## Mean-based Tests for Zero Inflated Log-Normal
ziln.uni=function(mod,n0,n1){
  n=n0+n1
  p1=exp(mod$b.b[1]+mod$b.b[2])/(1+exp(mod$b.b[1]+mod$b.b[2]))
  p0=exp(mod$b.b[1])/(1+exp(mod$b.b[1]))
  n0c=n0*p0
  n1c=n1*p1
  nc=n0c+n1c
  M0=exp(mod$b.c[1]+mod$disp/2)*p0
  M1=exp(mod$b.c[1]+mod$b.c[2]+mod$disp/2)*p1
  
  DM=M1-M0
  VDM=exp(2*mod$b.c[1]+mod$disp)*p0*(1-p0)/n0+exp(2*mod$b.c[1]+2*mod$b.c[2]+mod$disp)*p1*(1-p1)/n1+p0^2*exp(2*mod$b.c[1]+mod$disp)*mod$disp/n0c+p1^2*exp(2*(mod$b.c[1]+mod$b.c[2]+mod$disp/2))*mod$disp/n1c+mod$disp^2/(2*nc)*(exp(mod$b.c[1]+mod$b.c[2]+mod$disp/2)*p1-exp(mod$b.c[1]+mod$disp/2)*p0)^2
  # Note: p0^2/n0c=p0/n0
  
  DMtest=DM^2/VDM
  DM.p=1-pchisq(DMtest,1)
  
  IRM=mod$b.c[2]+mod$b.b[2]+log(1-p1)-log(1-p0)
  VIRM=mod$disp*(nc/(n0c*n1c))+(1-p0)/(n0*p0)+(1-p1)/(n1*p1)
  # Note: nc/(n0c*n1c)=(n0*p0+n1*p1)/(n0*p0*n1*p1)=1/(n0*p0)+1/(n1*p1)
  IRMtest=IRM^2/VIRM
  IRM.p=1-pchisq(IRMtest,1)
  
  return(list(M0=M0,M1=M1,DM=DM,VDM=VDM,DMtest=DMtest,DM.p=DM.p,IRM=IRM,VIRM=VIRM,IRMtest=IRMtest,IRM.p=IRM.p))
}













