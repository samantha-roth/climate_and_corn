#### 
## from PICAR2
plotRF<-function(dat,rangeDat=dat,label="Plot",location,length.out=10,pch=16,cex=1){
  breaks <- seq(range(rangeDat,na.rm = TRUE)[1],
                range(rangeDat,na.rm = TRUE)[2],
                length.out=length.out)
  pal <- tim.colors(length(breaks)-1,alpha = 1)
  fb <- classIntervals(dat, n = length(pal),
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(x=location[,1],y=location[,2],col=col, pch=pch,cex=cex,
       main=label)
}

################################
## Using Ming-Hui Chen's paper in Journal of Computational and Graphical Stats.
hpd <- function(samp,p=0.05){
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],1:2]
  return(hpd)
}

# Exponential Covariance Function
expCov<-function(distMat,phi){
  exp(-distMat/phi)
}

sqeCov<-function(distMat,phi){
  exp(-0.5*(distMat/phi)^2)
}

matCov<-function(distMat,phi){
  (1+(sqrt(5)*(distMat/phi))+((5*distMat^2)/(3*(phi^2))))*exp(-(sqrt(5)*(distMat/phi)))
}


# Matern Cov Function + Acceptance Rate function
Matern <- function(d, param = c(scale = 1, range = 1, smoothness = 2)) {
  scale <- param[1]
  range <- param[2]
  smoothness <- param[3]
  if (any(d < 0))
    stop("distance argument must be nonnegative")
  d <- d / range
  d[d == 0] <- 1e-10
  rootcon<-sqrt(2*smoothness)
  con <- (2^(smoothness - 1)) * gamma(smoothness)
  con <- 1 / con
  return(scale * con * ((rootcon*d)^smoothness) * besselK(rootcon*d, smoothness))
}

accRateFunc<-function(x){
  accRate<-(length(unique(x))-1)/(length(x)-1)
  return(accRate)
}

# Summary 
summaryFunction<-function(mcmcDat,bmseThresh=0.01,time){
  
  # Parameters
  summaryMat<-rbind(apply(mcmcDat,2,mean),
                    apply(mcmcDat,2,hpd),
                    apply(mcmcDat,2,accRateFunc),
                    bmmat(mcmcDat)[,2],
                    abs(apply(mcmcDat,2,mean))*bmseThresh,
                    apply(mcmcDat,2,ess),
                    apply(mcmcDat,2,ess)/time)
  
  rownames(summaryMat)<-c("Mean","95%CI-Low","95%CI-High",
                          "Accept","BMSE",paste(bmseThresh,"x mean"),
                          "ESS","ESS/sec")
  return(summaryMat)
}



# Rank Selection
MLE_FindRank_Poisson<-function(XMat,XMatCV,dimSeq,AMat,AMatCV,obsCV,obsMod,MoransOperatorEig){
  CVMSPE<-vector("numeric")
  betaMatList<-list() # Contains Beta Parameters
  for(i in 1:ncol(XMat)){betaMatList[[i]]<-matrix(NA,nrow=length(dimSeq),ncol=3)}
  for(jk in 1:length(dimSeq)){
    if(jk%%50==0){print(jk)}
    keepM<-1:dimSeq[jk]
    mBase<-(AMat%*%MoransOperatorEig$vectors[,keepM])
    mBaseCV<-(AMatCV%*%MoransOperatorEig$vectors[,keepM])
    lm1<-glm(obsMod~0+cbind(as.matrix(mBase),XMat),family = "poisson")  
    coeffs<-lm1$coefficients
    estMean<-coeffs[-(1:ncol(mBase))]
    lowCI<-estMean-1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))]
    highCI<-estMean+1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))]
    for(k in 1:length(betaMatList)){betaMatList[[k]][jk,]<-rbind(estMean,lowCI,highCI)[,k]}
    predCV<-exp(cbind(as.matrix(mBaseCV),XMatCV)%*%coeffs)  
    CVMSPE[jk]<-mean((predCV-obsCV)^2)
  }
  
  return(list(CVMSPE,betaMatList))
}


MLE_FindRank_Binary<-function(XMat,XMatCV,dimSeq,AMat,AMatCV,obsCV,obsMod,MoransOperatorEig){
  CVMSPE<-vector("numeric")
  betaMatList<-list() # Contains Beta Parameters
  for(i in 1:ncol(XMat)){betaMatList[[i]]<-matrix(NA,nrow=length(dimSeq),ncol=3)}
  for(jk in 1:length(dimSeq)){
    if(jk%%50==0){print(jk)}
    keepM<-1:dimSeq[jk]
    mBase<-(AMat%*%MoransOperatorEig$vectors[,keepM])
    mBaseCV<-(AMatCV%*%MoransOperatorEig$vectors[,keepM])
    lm1<-glm(obsMod~0+cbind(as.matrix(mBase),XMat),family = "binomial")  
    coeffs<-lm1$coefficients
    estMean<-coeffs[-(1:ncol(mBase))]
    lowCI<-estMean-1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))]
    highCI<-estMean+1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))]
    for(k in 1:length(betaMatList)){betaMatList[[k]][jk,]<-rbind(estMean,lowCI,highCI)[,k]}
    foo<-exp(cbind(as.matrix(mBaseCV),XMatCV)%*%coeffs)  
    predCV<-foo/(1+foo)
    predCV<-ifelse(predCV>0.5,1,0)
    CVMSPE[jk]<-mean((predCV-obsCV)^2)
  }
  
  return(list(CVMSPE,betaMatList))
}



MLE_FindRank_Linear<-function(XMat,XMatCV,dimSeq,AMat,AMatCV,obsCV,obsMod,MoransOperatorEig){
  CVMSPE<-vector("numeric")
  betaMatList<-list() # Contains Beta Parameters
  for(i in 1:ncol(XMat)){betaMatList[[i]]<-matrix(NA,nrow=length(dimSeq),ncol=3)}
  for(jk in 1:length(dimSeq)){
    if(jk%%50==0){print(jk)}
    keepM<-1:dimSeq[jk]
    mBase<-(AMat%*%MoransOperatorEig$vectors[,keepM])
    mBaseCV<-(AMatCV%*%MoransOperatorEig$vectors[,keepM])
    lm1<-glm(obsMod~0+cbind(as.matrix(mBase),XMat),family = "gaussian")  
    coeffs<-lm1$coefficients
    estMean<-coeffs[-(1:ncol(mBase))]
    lowCI<-estMean-1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))]
    highCI<-estMean+1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))]
    for(k in 1:length(betaMatList)){betaMatList[[k]][jk,]<-rbind(estMean,lowCI,highCI)[,k]}
    predCV<-cbind(as.matrix(mBaseCV),XMatCV)%*%coeffs
    CVMSPE[jk]<-mean((predCV-obsCV)^2)
  }
  
  return(list(CVMSPE,betaMatList))
}


# Cross- Validation
cvFunction.pois<-function(mcmcDat,AMatCV,mBase,XMatCV,gpWCV,obsCV,burnin=0.5*nrow(mcmcDat[[1]])){
  basisMat<-AMatCV%*%mBase
  cvPred<-tcrossprod(basisMat,mcmcDat[[2]][-(1:burnin),]) # Random Effects Prediction
  betaInd<-grep("beta",colnames(mcmcDat[[1]]))
  cvPreXB<-XMatCV%*%t(mcmcDat[[1]][-(1:burnin),betaInd]) # Mean prediction
  
  LoglambdaPred<-cvPred+cvPreXB
  cvPredW<-exp(apply(LoglambdaPred,1,mean))
  
  predVal<-cvPredW # MCMCResults
  
  predVal2<-pWCVPois # True
  cvSummary<-c(mean((predVal-obsCVPois)^2), #CVMSPE
               truthCVMSPEPois)
  predValuesMat<-cbind(predVal,predVal2,obsCVPois) #Combine
  return(list(cvSummary,predValuesMat))
}

cvFunction.bin<-function(mcmcDat,AMatCV,mBase,XMatCV,obsCV,burnin=0.5*nrow(mcmcDat[[1]])){
  basisMat<-AMatCV%*%mBase
  cvPred<-tcrossprod(basisMat,mcmcDat[[2]][-(1:burnin),]) # Random Effects Prediction
  betaInd<-grep("beta",colnames(mcmcDat[[1]]))
  cvPreXB<-XMatCV%*%t(mcmcDat[[1]][-(1:burnin),betaInd]) # Mean prediction
  foo<-apply(cvPred+cvPreXB,1,mean)
  cvPredW<-exp(foo)/(1+exp(foo))
  
  predVal<-ifelse(cvPredW>0.5,1,0) # MCMCResults
  
  predVal2<-pWCVBin # True
  cvSummary<-c(mean((predVal-obsCV)^2), #CVMSPE
               truthCVMSPEBin)
  predValuesMat<-cbind(predVal,predVal2,obsCV) #Combine
  return(list(cvSummary,predValuesMat))
}


cvFunction.linear<-function(mcmcDat,AMatCV,mBase,XMatCV,gpWCV,obsCV,burnin=0.5*nrow(mcmcDat[[1]])){
  basisMat<-AMatCV%*%mBase
  cvPred<-tcrossprod(basisMat,mcmcDat[[2]][-(1:burnin),]) # Random Effects Prediction
  betaInd<-grep("beta",colnames(mcmcDat[[1]]))
  cvPreXB<-XMatCV%*%t(mcmcDat[[1]][-(1:burnin),betaInd]) # Mean prediction
  predVal<-apply(cvPred+cvPreXB,1,mean)
  predVal2<-pWCVLinear # True
  cvSummary<-c(mean((predVal-obsCV)^2), #CVMSPE
               truthCVMSPELinear)
  predValuesMat<-cbind(predVal,predVal2,obsCV) #Combine
  return(list(cvSummary,predValuesMat))
}

cvFunction.linear.nogpWCV<-function(mcmcDat,AMatCV,mBase,XMatCV,obsCV,burnin=0.5*nrow(mcmcDat[[1]])){
  basisMat<-AMatCV%*%mBase
  cvPred<-tcrossprod(basisMat,mcmcDat[[2]][-(1:burnin),])`` # Random Effects Prediction
  betaInd<-grep("beta",colnames(mcmcDat[[1]]))
  cvPreXB<-XMatCV%*%t(mcmcDat[[1]][-(1:burnin),betaInd]) # Mean prediction
  predVal<-apply(cvPred+cvPreXB,1,mean)
  #predVal2<-pWCVLinear # True
  cvSummary<-mean((predVal-obsCV)^2) #CVMSPE
  predValuesMat<-cbind(predVal,obsCV) #Combine
  return(list(cvSummary,predValuesMat))
}

