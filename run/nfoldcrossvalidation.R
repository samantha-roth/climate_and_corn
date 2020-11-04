#Compare the PICAR and Full Gaussian Models using Bayesian Cross-Validation

setwd("C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR")

lpd_picar_mesh0_dimseq100_niter100k<- runCrossValidate(picarConf, k= 1589, foldFunction = "random", lossFunction= "predictive")
