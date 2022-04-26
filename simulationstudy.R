library(pracma)
library(cubature)
library(expm)
library("copula")
library(CVTuningCov)
library(JM)
library(joineR)
library("abind")
library("magic")
library("matrixcalc")
library(Deriv)
library(mvtnorm)
library(emdbook)
library("rgl")
library(latex2exp)
library(tikzDevice)
options(tikzMetricPackages = c("\\usepackage{amsmath}",
                               "\\usepackage{xcolor}",
                               "\\usepackage{tikz}",
                               "\\usetikzlibrary{calc}"))


GaussianHermite <- function(r){                     #function for generating Gaussian hermite nodes and their corresponding weights
  esp <- 3.0e-14
  pim4 <- 1/sqrt(sqrt(pi)) 	
  maxit <- 1000
  
  x <- w <- rep(0,r) 
  m <- (r+1 - ((r+1)%%2))/2
  p <- 2*r
  q <- 2*r+1
  
  for(i in 1:m){
    if(i==1) z <- sqrt(q) -1.85575*q^(-.1667)
    if(i==2) z <- z-1.14^((.425)/z)
    if(i==3) z <- 1.86*z - 0.86*x[1]
    if(i==4) z <- 1.91*z - 0.91*x[2]
    if(i!=1 & i!=2 & i!=3 & i!=4)  z <- 2*z - x[i-2]
    
    its <- 1
    dif <- 1
    while(its<maxit & dif>esp){
      p1 <- pim4
      p2 <- 0
      for(j in 1:r){
        p3 <- p2
        p2 <- p1
        p1 <- z*sqrt(2/j)*p2 - sqrt((j-1)/j)*p3
      }
      pp <- sqrt(p)*p2
      z <- z - p1/pp
      dif <- abs(p1/pp)
      its <- its +1
    }
    
    x[i] <- z
    x[r+1-i] <- -z
    w[i] <- 2/(pp^2)
    w[r+1-i] <- w[i]
  }
  return(list(round(w,15),round(x,15)))
}


car1 <- function(timepoint, rho) {                       #function for generating continuous AR1 correlation matrix
  exponent <- abs(matrix(timepoint, nrow =length(timepoint), ncol =length(timepoint), byrow = TRUE) - 
                    (timepoint))
  rho^exponent
}


#Function for generating data from the proposed joint model with up to ni longitudinal measurements between time (0,10) and an observed event time.
Y.lmmsurcondicopulajoineRcorrNOPLOTni=function(beta1,beta2,n,sigma,var.random,lambda,alpha,rate,trprob,genprob,ageprob,
                                               rho1,rho2,surlmmcor,lmmcor,ni)
{
  timepoint=seq(0,10,length=ni)
  if(surlmmcor=="car") R21=rho1^(c(10-timepoint+1))       #Select the correlation structure R_ty between the two sub-models
  if(surlmmcor=="ex") R21=rep(rho1,ni)                    #after conditioning on the random effects b_i
  R12=t(R21)
  R11=1
  if(lmmcor=="ex") R22=matrix(c(rep(c(1,rep(rho2,ni)),(ni-1)),1),ncol=ni)   #Select the correlation structure R_y within the longitudinal process
  if(lmmcor=="car") R22=car1(timepoint,rho2)                                 #after conditioning on the random effects b_i
  R=rbind(cbind(R11,R12),cbind(R21,R22))
  Y=matrix(0,nrow=ni*n,ncol=13)
  Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
  G.eigen=eigen(var.random)
  G.sqrt=G.eigen$vectors%*%diag(sqrt(G.eigen$values))%*%t(G.eigen$vectors)
  for(subj in 1:n)
  {
    treat=rbinom(1,1,trprob)
    gender=rbinom(1,1,genprob)
    age=sample(c(0,1,2),1,prob=ageprob)
    Xi1=matrix(c(c(Zi1),rep(treat,ni),rep(gender,ni),rep(as.numeric(age==1),ni)
                 ,rep(as.numeric(age==2),ni)),ncol=6)
    Xi2=matrix(c(treat,gender,as.numeric(age==1),as.numeric(age==2)),ncol=4)
    ceni=rexp(1,rate)
    bi=G.sqrt%*%rnorm(2)
    Yi=Xi1%*%beta1+Zi1%*%bi
    Ci=lambda*exp(beta2[1]*treat+beta2[2]*gender
                  +beta2[3]*Xi2[3]+beta2[4]*Xi2[4]+alpha[1]*bi[1])
    paramar=list(list(min=0,max=1))
    for(i in 2:(ni+1))
    {
      paramar[[i]]=list(mean=Yi[i-1],sd=sigma)
    }         
    Joint=rMvdc(1,mvdc(copula = ellipCopula(family = "normal", dim=(ni+1),dispstr = "un",
                                            param =c(R[lower.tri(R,diag=F)])),margins = c("unif", rep("norm",ni)),paramMargins=paramar))
    Ui=1-Joint[1]
    if(alpha[2]==0) {
      Ti=-log(Ui)/Ci
    } else {
      if((1-alpha[2]*bi[2]*log(Ui)/Ci)<=0) {
        Ti=max(ceni,50)
      } else {
        Ti=1/(alpha[2]*bi[2])*log(1-alpha[2]*bi[2]*log(Ui)/Ci)
      }}
    Yi=Joint[2:(ni+1)]
    Tiplot=min(timepoint[sum(timepoint<Ti)],timepoint[sum(timepoint<ceni)],10)
    Tiobs=min(Ti,ceni,10)
    indi=as.numeric(Ti<=ceni&&Ti<=10)
    Y[((subj-1)*ni+1):(subj*ni),1]=Yi
    Y[((subj-1)*ni+1):(subj*ni),2]=Xi1[,2]
    Y[((subj-1)*ni+1):(subj*ni),3]=subj
    Y[((subj-1)*ni+1):(subj*ni),4]=Ti 
    Y[((subj-1)*ni+1):(subj*ni),5]=ceni
    Y[((subj-1)*ni+1):(subj*ni),6]=Tiobs
    Y[((subj-1)*ni+1):(subj*ni),7]=Tiplot
    Y[((subj-1)*ni+1):(subj*ni),8]=indi
    Y[((subj-1)*ni+1):(subj*ni),9]=treat
    Y[((subj-1)*ni+1):(subj*ni),10]=gender
    Y[((subj-1)*ni+1):(subj*ni),11]=age
    Y[((subj-1)*ni+1):(subj*ni),12]=bi[1]
    Y[((subj-1)*ni+1):(subj*ni),13]=bi[2]
  }
  Y=data.frame(Y)
  colnames(Y)=c('resp','ti','subject','surti','centi','obsti','ploti','indicator'
                ,'treat','gender','age','bi0','bi1')
  return(Y)
}


#log-likelihood function for the model with extra correlation (R_ty) between the two sub-models and (R_y) within longitudinal process after 
#conditioning on the random effects. 
apploglikObmatrixni=function(theta,data,surlmmcor,lmmcor,m)
{
  beta1=theta[1:6]
  beta2=theta[7:10]
  D11=theta[11]
  D22=theta[12]
  D12=theta[13]
  rho1=theta[14]
  rho2=theta[15]
  lambda=theta[16]
  sigma=theta[17]
  alpha=theta[18]
  timepoint=data$ti[data$subject==1]
  ni=length(timepoint)
  if(surlmmcor=="car") R21=rho1^(10-timepoint+1)
  if(surlmmcor=="ex") R21=rep(rho1,ni)
  R12=t(R21)
  R11=1
  if(lmmcor=="ex") R22=matrix(c(rep(c(1,rep(rho2,ni)),(ni-1)),1),ncol=ni)
  if(lmmcor=="car") R22=car1(timepoint,rho2)
  R=rbind(cbind(R11,R12),cbind(R21,R22))
  Sigma=diag(c(1,rep(sigma,ni)))%*%R%*%diag(c(1,rep(sigma,ni)))
  Sigmay=Sigma[-1,-1]
  n=length(unique(data$subject))
  a=GaussianHermite(m)[[2]]
  w=GaussianHermite(m)[[1]]
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  ll=-1.0e40
  if(is.positive.definite(R)>0&sigma>0&lambda>0&is.positive.definite(D)>0)
  {
    ll=0
    for(sub in 1:n)
    {
      Xi1=matrix(c(rep(1,ni),timepoint,data$treat[data$subject==sub],
                   data$gender[data$subject==sub],as.numeric(data$age[data$subject==sub]==1),
                   as.numeric(data$age[data$subject==sub]==2)),ncol=6)
      Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
      Xi2=c(unique(data$treat[data$subject==sub]),unique(data$gender[data$subject==sub]),
            unique(as.numeric(data$age[data$subject==sub]==1)),
            unique(as.numeric(data$age[data$subject==sub]==2)))
      obsti=unique(data$obsti[data$subject==sub])
      indi=unique(data$indicator[data$subject==sub])
      dimi=sum(timepoint<=obsti)+1
      dimyi=dimi-1
      Sigmai=Sigma[1:dimi,1:dimi]
      Sigmayi=matrix(Sigmay[1:dimyi,1:dimyi],ncol=dimyi)
      V11=matrix(Zi1[1:dimyi,],ncol=2)%*%D%*%t(matrix(Zi1[1:dimyi,],ncol=2))+Sigmayi
      V22=D
      V12=matrix(Zi1[1:dimyi,],ncol=2)%*%D
      V21=t(V12)
      mubicon=V21%*%solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)[1:dimyi]
      Sigmabicon=V22-V21%*%solve(V11)%*%V12
      roti1=matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),ncol=2)
      roti=eigen(Sigmabicon)$vector%*%diag(sqrt(eigen(Sigmabicon)$value))%*%roti1
      bij=sqrt(2)*roti%*%t(as.matrix(expand.grid(a,a)))+c(mubicon)
      wij=expand.grid(w,w)[,1]*expand.grid(w,w)[,2]/pi
      Cij=lambda*exp(c(Xi2%*%beta2)+alpha*bij[1,])
      hij=Cij*exp(alpha*bij[2,]*obsti)
      Sij=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*obsti)-1))
      fij=hij*Sij
      Zyij=(matrix(rep(data$resp[data$subject==sub],m^2),ncol=m^2)-
              matrix(rep(Xi1%*%beta1,m^2),ncol=m^2)-Zi1%*%bij)[1:dimyi,]
      mutijcon=t(Sigmai[1,-1])%*%solve(Sigmayi)%*%Zyij
      sigmatijcon=sqrt(1-t(Sigmai[1,-1])%*%solve(Sigmayi)%*%Sigmai[1,-1])
      ftijcon=indi*dnorm(-qnorm(Sij),mutijcon,sigmatijcon)*fij/dnorm(-qnorm(Sij))+
        (1-indi)*(pnorm((qnorm(Sij)+mutijcon)/c(sigmatijcon)))
      posticondimean=sum(ftijcon*wij)     
      ll=ll+(-dimyi/2*log(2*pi)-0.5*log(det(V11))-0.5*t((data$resp[data$subject==sub]-Xi1%*%beta1)[1:dimyi])%*%
               solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)[1:dimyi]+log(posticondimean))
    }
  }
  if(is.na(ll)) ll=-1.0e40
  if(abs(ll)==Inf) ll=-1.0e40
  return(-ll)
}

#log-likelihood function for the model with R_ty=0 between the two sub-models and R_y=I within longitudinal process after 
#conditioning on the random effects. This is equivalent to a conventional joint model assuming conditional independence.
apploglikObuncormatrixni=function(theta,data,m)
{
  beta1=theta[1:6]
  beta2=theta[7:10]
  D11=theta[11]
  D22=theta[12]
  D12=theta[13]
  lambda=theta[14]
  sigma=theta[15]
  alpha=theta[16]
  timepoint=data$ti[data$subject==1]
  ni=length(timepoint)
  n=length(unique(data$subject))
  a=GaussianHermite(m)[[2]]
  w=GaussianHermite(m)[[1]]
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  ll=-1.0e40
  if(sigma>0&lambda>0&is.positive.definite(D)>0)
  {
    ll=0
    for(sub in 1:n)
    {
      Xi1=matrix(c(rep(1,ni),timepoint,data$treat[data$subject==sub],
                   data$gender[data$subject==sub],as.numeric(data$age[data$subject==sub]==1),
                   as.numeric(data$age[data$subject==sub]==2)),ncol=6)
      Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
      Xi2=c(unique(data$treat[data$subject==sub]),unique(data$gender[data$subject==sub]),
            unique(as.numeric(data$age[data$subject==sub]==1)),
            unique(as.numeric(data$age[data$subject==sub]==2)))
      obsti=unique(data$obsti[data$subject==sub])
      indi=unique(data$indicator[data$subject==sub])
      dimyi=sum(timepoint<=obsti)
      Sigmayi=sigma^2*diag(ni)[1:dimyi,1:dimyi]
      V11=matrix(Zi1[1:dimyi,],ncol=2)%*%D%*%t(matrix(Zi1[1:dimyi,],ncol=2))+Sigmayi
      V22=D
      V12=matrix(Zi1[1:dimyi,],ncol=2)%*%D
      V21=t(V12)
      mubicon=V21%*%solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)[1:dimyi]
      Sigmabicon=V22-V21%*%solve(V11)%*%V12
      roti1=matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),ncol=2)
      roti=eigen(Sigmabicon)$vector%*%diag(sqrt(eigen(Sigmabicon)$value))%*%roti1
      bij=sqrt(2)*roti%*%t(as.matrix(expand.grid(a,a)))+c(mubicon)
      wij=expand.grid(w,w)[,1]*expand.grid(w,w)[,2]/pi
      Cij=lambda*exp(c(Xi2%*%beta2)+alpha*bij[1,])
      hij=Cij*exp(alpha*bij[2,]*obsti)
      Sij=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*obsti)-1))
      fij=hij*Sij
      Zyij=(matrix(rep(data$resp[data$subject==sub],m^2),ncol=m^2)-
              matrix(rep(Xi1%*%beta1,m^2),ncol=m^2)-Zi1%*%bij)[1:dimyi,]
      ftijcon=indi*fij+(1-indi)*Sij
      posticondimean=sum(ftijcon*wij)     #equivalent to denoij above except constant pi^{-p/2} weight
      ll=ll+(-dimyi/2*log(2*pi)-0.5*log(det(V11))-0.5*t((data$resp[data$subject==sub]-Xi1%*%beta1)[1:dimyi])%*%
               solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)[1:dimyi]+log(posticondimean))
    }
  }
  if(is.na(ll)) ll=-1.0e40
  if(abs(ll)==Inf) ll=-1.0e40
  return(-ll)
}

#log-likelihood function for the model with R_ty=0 between the two sub-models after conditioning on the random effects. 
#This is model assuming conditional independence between the two sub-models but not within the longitudinal process.
apploglikObLongmatrixni=function(theta,data,lmmcor,m)
{
  beta1=theta[1:6]
  beta2=theta[7:10]
  D11=theta[11]
  D22=theta[12]
  D12=theta[13]
  rho2=theta[14]
  lambda=theta[15]
  sigma=theta[16]
  alpha=theta[17]
  timepoint=data$ti[data$subject==1]
  ni=length(timepoint)
  R21=rep(0,ni)
  R12=t(R21)
  R11=1
  if(lmmcor=="ex") R22=matrix(c(rep(c(1,rep(rho2,ni)),(ni-1)),1),ncol=ni)
  if(lmmcor=="car") R22=car1(timepoint,rho2)
  R=rbind(cbind(R11,R12),cbind(R21,R22))
  Sigma=diag(c(1,rep(sigma,ni)))%*%R%*%diag(c(1,rep(sigma,ni)))
  Sigmay=Sigma[-1,-1]
  n=length(unique(data$subject))
  a=GaussianHermite(m)[[2]]
  w=GaussianHermite(m)[[1]]
  D=matrix(c(D11,D12,D12,D22),ncol=2)
  ll=-1.0e40
  if(is.positive.definite(R)>0&sigma>0&lambda>0&is.positive.definite(D)>0)
  {
    ll=0
    for(sub in 1:n)
    {
      Xi1=matrix(c(rep(1,ni),timepoint,data$treat[data$subject==sub],
                   data$gender[data$subject==sub],as.numeric(data$age[data$subject==sub]==1),
                   as.numeric(data$age[data$subject==sub]==2)),ncol=6)
      Zi1=matrix(c(rep(1,ni),timepoint),ncol=2)
      Xi2=c(unique(data$treat[data$subject==sub]),unique(data$gender[data$subject==sub]),
            unique(as.numeric(data$age[data$subject==sub]==1)),
            unique(as.numeric(data$age[data$subject==sub]==2)))
      obsti=unique(data$obsti[data$subject==sub])
      indi=unique(data$indicator[data$subject==sub])
      dimi=sum(timepoint<=obsti)+1
      dimyi=dimi-1
      Sigmai=Sigma[1:dimi,1:dimi]
      Sigmayi=matrix(Sigmay[1:dimyi,1:dimyi],ncol=dimyi)
      V11=matrix(Zi1[1:dimyi,],ncol=2)%*%D%*%t(matrix(Zi1[1:dimyi,],ncol=2))+Sigmayi
      V22=D
      V12=matrix(Zi1[1:dimyi,],ncol=2)%*%D
      V21=t(V12)
      mubicon=V21%*%solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)[1:dimyi]
      Sigmabicon=V22-V21%*%solve(V11)%*%V12
      roti1=matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),ncol=2)
      roti=eigen(Sigmabicon)$vector%*%diag(sqrt(eigen(Sigmabicon)$value))%*%roti1
      bij=sqrt(2)*roti%*%t(as.matrix(expand.grid(a,a)))+c(mubicon)
      wij=expand.grid(w,w)[,1]*expand.grid(w,w)[,2]/pi
      Cij=lambda*exp(c(Xi2%*%beta2)+alpha*bij[1,])
      hij=Cij*exp(alpha*bij[2,]*obsti)
      Sij=exp(-Cij/(alpha*bij[2,])*(exp(alpha*bij[2,]*obsti)-1))
      fij=hij*Sij
      Zyij=(matrix(rep(data$resp[data$subject==sub],m^2),ncol=m^2)-
              matrix(rep(Xi1%*%beta1,m^2),ncol=m^2)-Zi1%*%bij)[1:dimyi,]
      ftijcon=indi*fij+(1-indi)*Sij
      posticondimean=sum(ftijcon*wij)     #equivalent to denoij above except constant pi^{-p/2} weight
      ll=ll+(-dimyi/2*log(2*pi)-0.5*log(det(V11))-0.5*t((data$resp[data$subject==sub]-Xi1%*%beta1)[1:dimyi])%*%
               solve(V11)%*%(data$resp[data$subject==sub]-Xi1%*%beta1)[1:dimyi]+log(posticondimean))
    }
  }
  if(is.na(ll)) ll=-1.0e40
  if(abs(ll)==Inf) ll=-1.0e40
  return(-ll)
}


#Estimated parameters from the proposed joint model with R_ty not equal to 0 and R_y not equal to identity matrix I
Y.lmmsurcondicopulajoineRGHObcorest=function(beta1,beta2,sigma,var.random,lambda,alpha,rho1,rho2,N,m,Jointdata,applogli,surlmmcor,lmmcor)
{
  gradmatrix=matrix(0,nrow=N,ncol=18)
  estmatrix=matrix(0,nrow=N,ncol=18)
  samplestmatrix=matrix(0,nrow=1,ncol=18)
  sdmatrix=matrix(0,nrow=N,ncol=18)
  biasmatrix=matrix(0,nrow=N,ncol=18)
  tvaluematrix=matrix(0,nrow=1,ncol=18)
  colnames(estmatrix)=c('beta01','beta11','beta21','beta31','beta41','beta51',
                        'beta12','beta22','beta32','beta42',
                        'D11','D22','D12','rho1','rho2',
                        'lambda','sigma','alpha')
  colnames(gradmatrix)=c('gbeta01','gbeta11','gbeta21','gbeta31','gbeta41','gbeta51',
                         'gbeta12','gbeta22','gbeta32','gbeta42',
                         'gD11','gD22','gD12','grho1','grho2',
                         'glambda','gsigma','galpha')
  colnames(samplestmatrix)=c('ssdbeta01','ssdbeta11','ssdbeta21','ssdbeta31','ssdbeta41','ssdbeta51',
                             'ssdbeta12','ssdbeta22','ssdbeta32','ssdbeta42',
                             'ssdD11','ssdD22','ssdD12','ssdrho1','ssdrho2',
                             'ssdlambda','ssdsigma','ssdalpha')
  colnames(sdmatrix)=c('sdbeta01','sdbeta11','sdbeta21','sdbeta31','sdbeta41','sdbeta51',
                       'sdbeta12','sdbeta22','sdbeta32','sdbeta42',
                       'sdD11','sdD22','sdD12','sdrho1','sdrho2',
                       'sdlambda','sdsigma','sdalpha')
  colnames(biasmatrix)=c('bbeta01','bbeta11','bbeta21','bbeta31','bbeta41','bbeta51',
                         'bbeta12','bbeta22','bbeta32','bbeta42',
                         'bD11','bD22','bD12','brho1','brho2',
                         'blambda','bsigma','balpha')
  colnames(tvaluematrix)=c('tbeta01','tbeta11','tbeta21','tbeta31','tbeta41','tbeta51',
                           'tbeta12','tbeta22','tbeta32','tbeta42',
                           'tD11','tD22','tD12','trho1','trho2',
                           'tlambda','tsigma','talpha')
  i=0
  k=1
  while(k<=N)
  {
    i=i+1
    estGHOb=nlm(f=applogli,p=c(beta1,beta2,c(var.random)[1],c(var.random)[4],c(var.random)[2],rho1,rho2,lambda,
                               sigma,unique(alpha)),data=Jointdata[[i]],surlmmcor=surlmmcor,lmmcor=lmmcor,m=m, hessian=T,iterlim=1000)
    if(estGHOb$code==1)
    {
      estmatrix[k,]=estGHOb$estimate
      sdmatrix[k,]=sqrt(diag(solve(estGHOb$hessian)))
      gradmatrix[k,]=estGHOb$gradient
      biasmatrix[k,]=estmatrix[k,]-c(beta1,beta2,c(var.random)[1],c(var.random)[4],c(var.random)[2],
                                     rho1,rho2,lambda,sigma,unique(alpha))
      k=k+1
    }
  }
  samplestmatrix[1,]=apply(estmatrix, 2, sd)
  tvaluematrix[1,]=colMeans(biasmatrix)/samplestmatrix[1,]
  results=list(gradmatrix,estmatrix,sdmatrix,colMeans(estmatrix),colMeans(biasmatrix),samplestmatrix,tvaluematrix,i)
  return(results)
}


#Estimated parameters from the proposed joint model with R_ty=0 and R_y=I. This is equivalent to a conventional joint model 
#assuming conditional independence.
Y.lmmsurcondicopulajoineRGHObuncorest=function(beta1,beta2,sigma,var.random,lambda,alpha,N,m,Jointdata,applogli)
{
  gradmatrix=matrix(0,nrow=N,ncol=16)
  estmatrix=matrix(0,nrow=N,ncol=16)
  samplestmatrix=matrix(0,nrow=1,ncol=16)
  sdmatrix=matrix(0,nrow=N,ncol=16)
  biasmatrix=matrix(0,nrow=N,ncol=16)
  tvaluematrix=matrix(0,nrow=1,ncol=16)
  colnames(estmatrix)=c('beta01','beta11','beta21','beta31','beta41','beta51',
                        'beta12','beta22','beta32','beta42',
                        'D11','D22','D12',
                        'lambda','sigma','alpha')
  colnames(gradmatrix)=c('gbeta01','gbeta11','gbeta21','gbeta31','gbeta41','gbeta51',
                         'gbeta12','gbeta22','gbeta32','gbeta42',
                         'gD11','gD22','gD12',
                         'glambda','gsigma','galpha')
  colnames(samplestmatrix)=c('ssdbeta01','ssdbeta11','ssdbeta21','ssdbeta31','ssdbeta41','ssdbeta51',
                             'ssdbeta12','ssdbeta22','ssdbeta32','ssdbeta42',
                             'ssdD11','ssdD22','ssdD12',
                             'ssdlambda','ssdsigma','ssdalpha')
  colnames(sdmatrix)=c('sdbeta01','sdbeta11','sdbeta21','sdbeta31','sdbeta41','sdbeta51',
                       'sdbeta12','sdbeta22','sdbeta32','sdbeta42',
                       'sdD11','sdD22','sdD12',
                       'sdlambda','sdsigma','sdalpha')
  colnames(biasmatrix)=c('bbeta01','bbeta11','bbeta21','bbeta31','bbeta41','bbeta51',
                         'bbeta12','bbeta22','bbeta32','bbeta42',
                         'bD11','bD22','bD12',
                         'blambda','bsigma','balpha')
  colnames(tvaluematrix)=c('tbeta01','tbeta11','tbeta21','tbeta31','tbeta41','tbeta51',
                           'tbeta12','tbeta22','tbeta32','tbeta42',
                           'tD11','tD22','tD12',
                           'tlambda','tsigma','talpha')
  i=0
  k=1
  while(k<=N)
  {
    i=i+1
    estGHOb=nlm(f=applogli,p=c(beta1,beta2,c(var.random)[1],c(var.random)[4],c(var.random)[2],lambda,
                               sigma,unique(alpha)),data=Jointdata[[i]],m=m, hessian=T,iterlim=1000)
    if(estGHOb$code==1)
    {
      estmatrix[k,]=estGHOb$estimate
      sdmatrix[k,]=sqrt(diag(solve(estGHOb$hessian)))
      gradmatrix[k,]=estGHOb$gradient
      biasmatrix[k,]=estmatrix[k,]-c(beta1,beta2,c(var.random)[1],c(var.random)[4],c(var.random)[2],
                                     lambda,sigma,unique(alpha))
      k=k+1
    }
  }
  samplestmatrix[1,]=apply(estmatrix, 2, sd)
  tvaluematrix[1,]=colMeans(biasmatrix)/samplestmatrix[1,]
  results=list(gradmatrix,estmatrix,sdmatrix,colMeans(estmatrix),colMeans(biasmatrix),samplestmatrix,tvaluematrix,i)
  return(results)
}


#Estimated parameters from the proposed joint model with R_ty=0 but R_y not equal to I. 
#This is model assuming conditional independence between the two sub-models but not within the longitudinal process.
Y.lmmsurcondicopulajoineRGHObLongest=function(beta1,beta2,sigma,var.random,lambda,alpha,rho2,N,m,Jointdata,applogli,lmmcor)
{
  gradmatrix=matrix(0,nrow=N,ncol=17)
  estmatrix=matrix(0,nrow=N,ncol=17)
  samplestmatrix=matrix(0,nrow=1,ncol=17)
  sdmatrix=matrix(0,nrow=N,ncol=17)
  biasmatrix=matrix(0,nrow=N,ncol=17)
  tvaluematrix=matrix(0,nrow=1,ncol=17)
  colnames(estmatrix)=c('beta01','beta11','beta21','beta31','beta41','beta51',
                        'beta12','beta22','beta32','beta42',
                        'D11','D22','D12','rho2',
                        'lambda','sigma','alpha')
  colnames(gradmatrix)=c('gbeta01','gbeta11','gbeta21','gbeta31','gbeta41','gbeta51',
                         'gbeta12','gbeta22','gbeta32','gbeta42',
                         'gD11','gD22','gD12','grho2',
                         'glambda','gsigma','galpha')
  colnames(samplestmatrix)=c('ssdbeta01','ssdbeta11','ssdbeta21','ssdbeta31','ssdbeta41','ssdbeta51',
                             'ssdbeta12','ssdbeta22','ssdbeta32','ssdbeta42',
                             'ssdD11','ssdD22','ssdD12','ssdrho2',
                             'ssdlambda','ssdsigma','ssdalpha')
  colnames(sdmatrix)=c('sdbeta01','sdbeta11','sdbeta21','sdbeta31','sdbeta41','sdbeta51',
                       'sdbeta12','sdbeta22','sdbeta32','sdbeta42',
                       'sdD11','sdD22','sdD12','sdrho2',
                       'sdlambda','sdsigma','sdalpha')
  colnames(biasmatrix)=c('bbeta01','bbeta11','bbeta21','bbeta31','bbeta41','bbeta51',
                         'bbeta12','bbeta22','bbeta32','bbeta42',
                         'bD11','bD22','bD12','brho2',
                         'blambda','bsigma','balpha')
  colnames(tvaluematrix)=c('tbeta01','tbeta11','tbeta21','tbeta31','tbeta41','tbeta51',
                           'tbeta12','tbeta22','tbeta32','tbeta42',
                           'tD11','tD22','tD12','trho2',
                           'tlambda','tsigma','talpha')
  i=0
  k=1
  while(k<=N)
  {
    i=i+1
    estGHOb=nlm(f=applogli,p=c(beta1,beta2,c(var.random)[1],c(var.random)[4],c(var.random)[2],rho2,lambda,
                               sigma,unique(alpha)),data=Jointdata[[i]],lmmcor=lmmcor,m=m, hessian=T,iterlim=1000)
    if(estGHOb$code==1)
    {
      estmatrix[k,]=estGHOb$estimate
      sdmatrix[k,]=sqrt(diag(solve(estGHOb$hessian)))
      gradmatrix[k,]=estGHOb$gradient
      biasmatrix[k,]=estmatrix[k,]-c(beta1,beta2,c(var.random)[1],c(var.random)[4],c(var.random)[2],rho2,
                                     lambda,sigma,unique(alpha))
      k=k+1
    }
  }
  samplestmatrix[1,]=apply(estmatrix, 2, sd)
  tvaluematrix[1,]=colMeans(biasmatrix)/samplestmatrix[1,]
  results=list(gradmatrix,estmatrix,sdmatrix,colMeans(estmatrix),colMeans(biasmatrix),samplestmatrix,tvaluematrix,i)
  return(results)
}



