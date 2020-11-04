rm(list=ls())

# Initialize
library(nimble);library(mvtnorm);library(fields)
setwd("C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR/run")
source(file = "C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR/source/sharedFunctions.R") # Shared Functions
source(file = "C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR/source/nimbleSource.R") # Nimble Functions
source(file="C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR/source/batchmeans.R")
load(file="C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR/samples/Crop1981_SpatialData.RData") # Spatial Data
#load(file="C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR/samples/mymesh1_newparams.RData") # Mesh
load(file="C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR/samples/mydatamesh2.RData") # Mesh
#Select Rank for PICAR
dimSeq<-seq(2, 200, by=1) # ADjust the endpoint (50) and the resolution (1)
heuristicResults<-MLE_FindRank_Linear(XMat=XMat,XMatCV=XMatCV,dimSeq=dimSeq,
                                      AMat=AMat,AMatCV=AMatCV,obsCV=obsCVLinear,obsMod=obsModLinear,
                                      MoransOperatorEig = MoransOperatorEig)

par(mfrow=c(2,2))
plot(x=dimSeq,y=heuristicResults[[1]],typ="l" , main="CVMSPE") # CVMSPE
for(k in 1:ncol(XMat)){
  plot(x=dimSeq,y=heuristicResults[[2]][[k]][,1],typ="l",ylim=range(heuristicResults[[2]][[k]],na.rm = TRUE),
       main=paste("Beta",k))
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,2],lty=2)
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,3],lty=2)
  
}

# Select Rank of Moran's Basis Functions
# We choose the rank that yields the lowest CVMSPE based on the above
pBase<-dimSeq[which.min(heuristicResults[[1]])] ; print(pBase)
mBase<-MoransOperatorEig$vectors[,1:pBase] # Moran's Basis Function
p<-ncol(mBase) #p is the number of bases functions
#p<- 20
X<-XMat
M<-as.matrix(AMat%*%mBase)
MQM<-diag(p) #For nimble, the identity matrix yields the least numerical issues 

# Preliminaries for NIMBLE
niter=100000
consts   <- list(n=n,p=p)
data     <- list(Z=obsModLinear,X=XMat,M=M,MQM=MQM,mn=rep(0,p))
inits    <- list(beta1=rnorm(1), beta2= rnorm(1), beta3= rnorm(1), beta4= rnorm(1), beta5= rnorm(1), beta6= rnorm(1), beta7= rnorm(1), beta8= rnorm(1), beta9= rnorm(1), tau=2, 
                 delta=rnorm(p))

# Build Model
Rmodel<-nimbleModel(code=linear_PICAR_string, data = data,  constants=consts, inits = inits)
Cpicar <- compileNimble(Rmodel)
picarConf <- configureMCMC(Rmodel, print = TRUE, multivariateNodesAsScalars=TRUE) # original statement: picarConf <- configureMCMC(Rmodel, print = TRUE)
picarConf$addMonitors(c("beta1", "beta2", "beta3","beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "tau", "delta"))
picarMCMC <- buildMCMC(picarConf)

# Run the mcmc algorithm
CpicarMCMC <- compileNimble(picarMCMC)
set.seed(0)
pt<-proc.time()
CpicarMCMC$run(niter)
ptFinal<-proc.time()-pt
samples <- as.matrix(CpicarMCMC$mvSamples)
samples_picar_mesh2_dimseq200<- samples
save(samples_picar_mesh2_dimseq200,file="C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR/output/samples_picar_mymesh2_dimseq200_niter100k.RData")

par(mfrow= c(3,3))
for(i in 1:9) plot(1:niter, samples[1:niter,i], main= paste("beta",i, " mean= ", mean(samples[1:niter,i])))

# samples_picar_mymesh1_dimseq50_niter100k looks like the Markov chains get to around the right place by 40k
# samples_picar_mymesh1_dimseq30_niter100k looks like the Markov chains get to around the right place by 40k
#par(mfrow= c(3,3))
#for(i in 1:9) plot(50001:niter, samples[50001:niter,i], main= paste("beta",i, " mean= ", mean(samples[50001:niter,i])))
# Summary
# Table
deltaInd<-grep("delta",colnames(samples))
summaryMat<-list()
summaryMat[[1]]<-round(summaryFunction(samples[,-deltaInd], time=ptFinal[3]),3)
summaryMat[[2]]<-round(summaryFunction(samples[,deltaInd], time=ptFinal[3]),3)
summaryMat[[1]]
apply(summaryMat[[2]],1,mean) # Mean results for the reparameterized random effects/basis coefficients
save(summaryMat,samples,ptFinal,file="C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR/output/linearPICAR_mesh2_dseq200_niter100k.RData")


.
# trace Plots and Posterior Density
# pdf(file = "../output/BinaryPICAR.pdf",width=11,height=8.5)
par(mfrow=c(5,2),mar=c(2,2,2,2))  
#sampInd<-floor(seq(1,nrow(samples),length.out = 1000))
sampInd<- seq(1,nrow(samples))
for(i in 1:9){
  plot.ts(samples[sampInd,paste("beta",i,sep="")], 
          main=paste("beta",i,sep="")); 
  #abline(h=c(beta,phi)[i],col="red",lwd=2)
  plot(density(samples[sampInd,paste("beta",i,sep="")]),
       main=paste("beta",i,sep="")); 
  #abline(v=beta[i],col="red",lwd=2)
}
# dev.off()

# Get Cross Validation

mcmcDat<-list()
mcmcDat[[1]]<-samples[,-deltaInd]
mcmcDat[[2]]<-samples[,deltaInd]
cvSummary<-cvFunction.linear.nogpWCV(mcmcDat=mcmcDat,AMatCV=AMatCV,mBase=mBase, XMatCV=XMatCV,obsCV=obsCVLinear)
#the original statement was cvSummary<-cvFunction.linear(mcmcDat=mcmcDat,AMatCV=AMatCV,mBase=mBase, 
#XMatCV=XMatCV,gpWCV=gpWCV,obsCV=obsCVLinear), however since our data is real, we don't know gpWCV, so we got rid of it
print(cvSummary[[1]]) #CVMSPE
