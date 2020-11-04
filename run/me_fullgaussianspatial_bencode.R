# Refitting the full Gaussian spatial model with Ben's nimble code instead of mine to see
# if PICAR performs any better in estimating the intercept term
#setwd("/gpfs/group/kzk10/default/private/svr5482/PICAR/PICAR")
setwd("C:/Users/saman/Dropbox/svr5482/PICAR/PICAR")
load("samples/Crop1981_SpatialData.RData")
library(batchmeans)
library(nimble,warn.conflicts = FALSE)

## define a function to calculate distances
dists<- function(df= df, n=nrow(df)){
  distmat<- matrix(NA, nrow = n, ncol = n)
  for(i in 1:(n-1)){
    i
    distmat[i,i]<-0
    for(j in (i+1):n){
      j
      distmat[i,j]<- sqrt((df$lat[i]-df$lat[j])^2 + (df$lon[i]-df$lon[j])^2 )
      distmat[j,i]<- distmat[i,j]
    }
  }
  distmat[n,n]<-0
  return(distmat)
}

# Covariance Functions for NIMBLE
# Exponential
# inputs: matrix of distances, a value for phi
expcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- exp(-dists[i,j]/phi)
      }
    }
    
    return(result)
  })
cExpcov <- compileNimble(expcov)


##Full Gaussian Spatial Model: Ben's code
# Nimble Model - Linear Spatial Model
# need initial values for: sigma2, tau2, phi, covMat, beta1 to beta9
# constants: n, X
# data: data <- list(y= train_data$Yield)
linear_model_string <- nimbleCode({
  
  # Data Model
  Z[1:n] ~ dmnorm(mean = mn[1:n], cov = fullCovMat[1:n,1:n])
  
  # Constant and Cov Matrix
  mn[1:n]<-beta1*X[1:n,1] + beta2*X[1:n,2] + beta3*X[1:n,3] +  beta4*X[1:n,4] + beta5*X[1:n,5] + beta6*X[1:n,6] + beta7*X[1:n,7] + beta8*X[1:n,8] + beta9*X[1:n,9]
  covMat[1:n,1:n]<- expcov(dists[1:n,1:n],phi) # Need to use the right Covariance function
  fullCovMat[1:n,1:n]<- sigma2*covMat[1:n,1:n]+tau2*diag(n)
  
  # Process Model
  # None needed because Ws are marginalized out 
  
  # Parameter Model
  # Set different priors
  sigma2   ~  dinvgamma(0.2, 0.2) 
  tau2   ~  dinvgamma(0.2, 0.2) 
  phi   ~  dunif(0,1)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
  beta4 ~  dnorm(0, sd=sqrt(100))
  beta5 ~  dnorm(0, sd=sqrt(100))
  beta6 ~  dnorm(0, sd=sqrt(100))
  beta7 ~  dnorm(0, sd=sqrt(100))
  beta8 ~  dnorm(0, sd=sqrt(100))
  beta9 ~  dnorm(0, sd=sqrt(100))
})

# Set initial values, constants, data
# Preliminaries for NIMBLE
niter=50000
consts   <- list(n= nrow(XMat), dists= distMatMod, X= XMat)
data     <- list(Z=obsModLinear)
inits    <- list(beta1=rnorm(1), beta2= rnorm(1), beta3= rnorm(1), beta4= rnorm(1), beta5= rnorm(1), beta6= rnorm(1), beta7= rnorm(1), beta8= rnorm(1), beta9= rnorm(1), 
                 tau2= 5, sigma2= 5, phi= dunif(1,0,1))
inits$covMat <- cExpcov(distMatMod, inits$phi)


# build model
Rmodel<-nimbleModel(code=linear_model_string, data = data,  constants=consts, inits = inits)
Cpicar <- compileNimble(Rmodel)
picarConf <- configureMCMC(Rmodel, print = TRUE) # original statement: picarConf <- configureMCMC(Rmodel, print = TRUE)
picarConf$addMonitors(c("beta1", "beta2", "beta3","beta4", "beta5", "beta6", "beta7", "beta8", "beta9", "tau2", "sigma2", "phi"))
picarMCMC <- buildMCMC(picarConf)

# Run the mcmc algorithm
CpicarMCMC <- compileNimble(picarMCMC)
set.seed(0)
pt<-proc.time()
CpicarMCMC$run(niter)
ptFinal<-proc.time()-pt
samples_fg <- as.matrix(CpicarMCMC$mvSamples)
save(samples_fg,file="C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR 2/PICAR/output/samples_fg.RData")

par(mfrow= c(3,3))
for(i in 1:9) plot(1:niter, samples_fg[1:niter,i], main= paste("beta",i, " mean= ", mean(samples_fg[1:niter,i])))

#burn in up to wherever the last MC seems to start sticking around the same place
par(mfrow= c(3,3))
for(i in 1:9) plot(25001:niter, samples_fg[25001:niter,i], main= paste("beta",i, " mean= ", mean(samples_fg[25001:niter,i])))

# trace Plots and Posterior Density
# pdf(file = "../output/BinaryPICAR.pdf",width=11,height=8.5)
par(mfrow=c(5,2),mar=c(2,2,2,2))  
#sampInd<-floor(seq(1,nrow(samples),length.out = 1000))
sampInd<- seq(1,nrow(samples_fg))
for(i in 1:9){
  plot.ts(samples_fg[sampInd,paste("beta",i,sep="")], 
          main=paste("beta",i,sep="")); 
  #abline(h=c(beta,phi)[i],col="red",lwd=2)
  plot(density(samples_fg[sampInd,paste("beta",i,sep="")]),
       main=paste("beta",i,sep="")); 
  #abline(v=beta[i],col="red",lwd=2)
}

## Predict on the test set to get the CVMSPE
bmb1<- bm(samples[,"beta1"])
bmb2<- bm(samples[,"beta2"])
bmb3<- bm(samples[,"beta3"])
bmb4<- bm(samples[,"beta4"])
bmb5<- bm(samples[,"beta5"])
bmb6<- bm(samples[,"beta6"])
bmb7<- bm(samples[,"beta7"])
bmb8<- bm(samples[,"beta8"])
bmb9<- bm(samples[,"beta9"])
bmsig2<- bm(samples[,"sigma2"])
bmtau2<- bm(samples[,"tau2"])
bmphi<- bm(samples[,"phi"])

fulldata<- rbind(test_data, train_data)
distmat<- dists(fulldata, n= nrow(fulldata))
n= nrow(fulldata)
est_covmat= matrix(NA, nrow= n, ncol= n)
for(i in 1:n){
  for(j in 1:n){
    if(j==i) est_covmat[i,j] <- bmsig2$est + bmtau2$est
    if(j!=i) est_covmat[i,j] <- bmsig2$est*exp(-distmat[i,j]/bmphi$est)
  }
}

#define the estimated covariance matrix between the training data and the training data
est_covmat_traintrain<- est_covmat[(nrow(test_data)+1):n,(nrow(test_data)+1):n]

#define the estimated covariance matrix between the training data and the testing data
est_covmat_traintest<- est_covmat[1:nrow(test_data),(nrow(test_data)+1):n]

mu_train<- bmb1$est + bmb2$est*train_data$GDD + bmb3$est*train_data$EDD + bmb4$est*train_data$Pr + bmb5$est*((train_data$Pr)^2) + bmb6$est*train_data$Tmax + bmb7$est*train_data$Tmin + bmb8$est*train_data$lastSF + bmb9$est*train_data$firstFF

pred_bayes_sp<- bmb1$est + bmb2$est*test_data$GDD + bmb3$est*test_data$EDD + bmb4$est*test_data$Pr + bmb5$est*((test_data$Pr)^2) + bmb6$est*test_data$Tmax + bmb7$est*test_data$Tmin + bmb8$est*test_data$lastSF + bmb9$est*test_data$firstFF + est_covmat_traintest%*%solve(est_covmat_traintrain)%*%(train_data$Yield - mu_train)

CV_bayes_sp<-mean((pred_bayes_sp - test_data$Yield)^2)
sse_bayes_sp<- sum((pred_bayes_sp - test_data$Yield)^2)
NPE_bayes_sp<-sqrt(mean(CV_bayes_sp))/mean(test_data[ ,"Yield"],na.rm=TRUE)
NPE_bayes_sp