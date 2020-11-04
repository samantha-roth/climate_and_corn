# Covariance Functions for NIMBLE
# Exponential
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

# Matern nu=2.5
matcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- (1+(sqrt(5)*(dists[i,j]/phi))+((5*dists[i,j]^2)/(3*(phi^2))))*exp(-(sqrt(5)*(dists[i,j]/phi)))
      }
    }
    
    return(result)
  })

######################################################################
######################################################################
# Full Hierarchical Spatial models - linear, binary, and counts
######################################################################
######################################################################

# Nimble Model - Linear Spatial Model
linear_model_string <- nimbleCode({
  
  # Data Model
  
  Z[1:n] ~ dmnorm(mean = mn[1:n], cov = fullCovMat[1:n,1:n])
  
  # Constant and Cov Matrix
  fullCovMat[1:n,1:n]<- sigma2*covMat[1:n,1:n]+tau2*diag(n)
  mn[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,3]
  covMat[1:n,1:n]<- expcov(dists[1:n,1:n],phi) # Need to use the right Covariance function
  
  
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
})


# Nimble Model - Binary
binary_model_string <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    Z[i] ~ dbinom(prob = prob[i] , size = 1)
    prob[i] <- exp(W[i]+XB[i])/(1+exp(W[i]+XB[i]))
  }
  
  # Constant and Cov Matrix
  XB[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,2]
  covMat[1:n,1:n]<- expcov(dists[1:n,1:n],phi)
  fullCovMat[1:n,1:n]<- sigma2*covMat[1:n,1:n]
  
  # Process Model
  W[1:n] ~ dmnorm(mean = mn[1:n], cov = fullCovMat[1:n,1:n])
  
  # Parameter Model
  sigma2   ~  dgamma(0.2, 0.2)
  phi   ~  dunif(0,1)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
})


# Nimble Model - Count
poisson_model_string <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    lambda[i] <- exp(W[i]+XB[i])
    Z[i] ~ dpois(lambda[i])
  }
  
  # Constant and Cov Matrix
  XB[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,2]
  covMat[1:n,1:n]<- expcov(dists[1:n,1:n],phi)
  fullCovMat[1:n,1:n]<- sigma2*covMat[1:n,1:n]
  
  # Process Model
  W[1:n] ~ dmnorm(mean = mn[1:n], cov = fullCovMat[1:n,1:n])
  
  # Parameter Model
  sigma2   ~  dinvgamma(0.2, 0.2)
  phi   ~  dunif(0,1)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
})

######################################################################
######################################################################
# PICAR-based Spatial models - linear, binary, and counts
######################################################################
######################################################################


linear_PICAR_string <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    Z[i] ~ dnorm(mean = meanV[i] , var = tau2)
  }
  # Mean Function
  meanV[1:n]<-XB[1:n]+W[1:n]
  XB[1:n]<- beta1*X[,1] + beta2*X[,2] + beta3*X[,3] + beta4*X[,4] + beta5*X[,5] + beta6*X[,6] + beta7*X[,7] + beta8*X[,8] + beta9*X[,9]
  W[1:n]<-M[1:n,1:p]%*%delta[1:p]
  precMat[1:p,1:p]<-tau * MQM[1:p,1:p]
  
  # Process Model
  delta[1:p] ~ dmnorm(mean = mn[1:p], prec = precMat[1:p,1:p])
  # Parameter Model
  tau   ~  dgamma(0.5,2000)
  tau2 ~  dinvgamma(0.2, 0.2) 
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


# Nimble Model - Binary
binary_PICAR_string <- nimbleCode({
  # Data Model
  for(i in 1:n){
    prob[i] <- exp(W[i]+XB[i])/(1+exp(W[i]+XB[i]))
    Z[i] ~ dbinom(prob = prob[i] , size = 1)
  }
  # Constant and Cov Matrix
  XB[1:n]<- beta0*X[,1] +beta1*X[,2] + beta2*X[,3] + beta3*X[,4] + beta4*X[,5] + beta5*X[,6] + beta6*X[,7] + beta7*X[,8] + beta8*X[,9]
  W[1:n]<-M[1:n,1:p]%*%delta[1:p]
  precMat[1:p,1:p]<-tau * MQM[1:p,1:p]
  # Process Model
  delta[1:p] ~ dmnorm(mean = mn[1:p], prec = precMat[1:p,1:p])
  # Parameter Model
  tau   ~  dgamma(0.5,2000)
  beta0 ~  dnorm(0, sd=sqrt(100))
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
  beta4 ~  dnorm(0, sd=sqrt(100))
  beta5 ~  dnorm(0, sd=sqrt(100))
  beta6 ~  dnorm(0, sd=sqrt(100))
  beta7 ~  dnorm(0, sd=sqrt(100))
  beta8 ~  dnorm(0, sd=sqrt(100))
})


# Nimble Model - Poisson
poisson_PICAR_string <- nimbleCode({
  # Data Model
  for(i in 1:n){
    lambda[i] <- exp(W[i]+XB[i])
    Z[i] ~ dpois(lambda[i])
  }
  # Constant and Cov Matrix
  XB[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,2]
  W[1:n]<-M[1:n,1:p]%*%delta[1:p]
  precMat[1:p,1:p]<-tau * MQM[1:p,1:p]
  # Process Model
  delta[1:p] ~ dmnorm(mean = mn[1:p], prec = precMat[1:p,1:p])
  # Parameter Model
  tau   ~  dgamma(0.5,1/2000)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
})
