rm(list=ls()) ## from PICAR 2
library(fields) ; library(mvtnorm) ; library(classInt)
setwd("~/Dropbox/PICAR/") # Set Directory
source(file = "source/sharedFunctions.R")

#Parameters
set.seed(2020)
n=500 ; ncv=100
beta=c(1,1,1) ; phi=0.2 ; sigma2=1
# Generate Locations 
## Split into Model + Cross-Validation
gridLocation<-cbind(runif(n,min = 0,max = 1),runif(n,min = 0,max = 1))
CVgridLocation<-cbind(runif(ncv,min = 0,max = 1),runif(ncv,min = 0,max = 1))
comboLocation<-rbind(gridLocation,CVgridLocation)
distMatFull<-as.matrix(rdist(comboLocation))
# Create Indices
modInd<-1:n
CVInd<-(n+1):nrow(distMatFull)
# Covariates
XMat<-cbind(runif(n,-1,1),runif(n,-1,1),runif(n,-1,1))
XMatCV<-cbind(runif(ncv,-1,1),runif(ncv,-1,1),runif(ncv,-1,1))
XB<-XMat%*%beta
cvXB<-XMatCV%*%beta
XBFull<-rbind(XB,cvXB)

# Covariance Matrix
# I used a Matern Cov function with parameters phi, sigma2, and fixed nu=2.5. 
# Nu can be adjusted by replacing matCov with expCov (exponential nu=0.5) or sqeCov (squared exponential nu=infinity)
CovMat<-sigma2*matCov(distMatFull,phi)

# Latent Gaussian Random Field
gpWFull <- as.numeric(rmvnorm(n=1,mean=rep(0,nrow(CovMat)),sigma = CovMat,method = "chol"))
pWFullLinear<-gpWFull+XBFull
pWFullPois<-exp(gpWFull+XBFull)
pWFullBin<-exp(gpWFull+XBFull)/(1+exp(gpWFull+XBFull))

# Observations
obsFullLinear<-pWFullLinear
obsFullPois<-sapply(pWFullPois,rpois,n=1)
obsFullBin<-sapply(pWFullBin,rbinom,n=1,size=1)


##################
# Create Model fitting and validation sample sets
##################
## Model Fitting Sample
gpWMod<-gpWFull[modInd]# Latent Process
pWModLinear<-pWFullLinear[modInd] # Expected Value Linear
pWModPois<-pWFullPois[modInd] # Expected Value Poisson
pWModBin<-pWFullBin[modInd] # Expected Value Binary
obsModLinear<-obsFullLinear[modInd] # Observations Linear
obsModPois<-obsFullPois[modInd] # Observations Poisson
obsModBin<-obsFullBin[modInd] # Observations Binary

## CV Sample
# Corresponds to model fitting sample
gpWCV<-gpWFull[CVInd]
pWCVLinear<-pWFullLinear[CVInd]
pWCVPois<-pWFullPois[CVInd]
pWCVBin<-pWFullBin[CVInd]
obsCVLinear<-obsFullLinear[CVInd]
obsCVPois<-obsFullPois[CVInd]
obsCVBin<-obsFullBin[CVInd]

# CVMSPE: Cross validation mean squared prediction error if we knew the truth
truthCVMSPELinear<-mean((pWCVLinear-obsCVLinear)^2) # Linear
truthCVMSPEPois<-mean((pWCVPois-obsCVPois)^2) # Poisson
truthCVMSPEBin<-mean((pWCVBin-obsCVBin)^2) # Binary

# Save Data
rm(CovMat, distMatFull) # Large matrices
save.image(file="samples/SpatialData.RData") # Save Data 


# Plots
par(mfrow=c(2,2))
plotRF(dat=obsFullLinear, location = gridLocation , label="Linear Observations")
plotRF(dat=obsFullPois, location = gridLocation , label="Count Observations")
plot(x=gridLocation[,1], y=gridLocation[,2], col = obsFullBin+1 , pch=16 , main="Binary Observations")

