rm(list=ls())
library(fields) ; library(mvtnorm) ; library(INLA)
#setwd("~/Dropbox/PICAR/") # Set Directory
setwd("C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR")
load("samples/Crop1981_SpatialData.RData")
source(file = "source/sharedFunctions.R")

# Make Mesh using INLA
#mesh0
mesh <- inla.mesh.2d(gridLocation, 
                     max.edge=c(0.05,0.5),
                     cutoff = 0.025,
                     offset=c(0.1, 0.1))
# my adjustments to mesh, mesh1
mesh <- inla.mesh.2d(gridLocation, 
                     max.edge=c(0.05,0.8),
                     min.angle= 10, #smaller minimum triangle angle
                     cutoff = 0.01, #smaller minimum distance between points
                     offset=c(0.1, 0.1))
# second adjustments to mesh, mesh2
mesh <- inla.mesh.2d(gridLocation, 
                     max.edge=c(1),
                     min.angle= 5, #smaller minimum triangle angle
                     cutoff = 0.005, #smaller minimum distance between points
                     offset=c(0.1, 0.1))
# third adjustments to mesh, mesh3
mesh <- inla.mesh.2d(gridLocation, 
                     max.edge=c(1),
                     min.angle= 5, #smaller minimum triangle angle
                     cutoff = 0.05, #smaller minimum distance between points
                     offset=c(0.1, 0.1))
# fourth adjustments to mesh, mesh4
mesh <- inla.mesh.2d(gridLocation, 
                     max.edge=c(1),
                     min.angle= 5, #smaller minimum triangle angle
                     cutoff = 0.0025, #smaller minimum distance between points
                     offset=c(0.1, 0.1))
# fifth adjustments to mesh, mesh5
mesh <- inla.mesh.2d(gridLocation, 
                     max.edge=c(1),
                     min.angle= 5, #minimum triangle angle
                     cutoff = 0.04, #minimum distance between points
                     offset=c(0.1, 0.1))

# sixth adjustments to mesh, mesh6
mesh <- inla.mesh.2d(gridLocation, 
                     max.edge=c(1),
                     min.angle= 5, #minimum triangle angle
                     cutoff = 0.03, #minimum distance between points
                     offset=c(0.1, 0.1))
# seventh adjustments to mesh, mesh7
mesh <- inla.mesh.2d(gridLocation, 
                     max.edge=c(0.05,0.8),
                     min.angle= 10, #smaller minimum triangle angle
                     cutoff = 0.02, #smaller minimum distance between points
                     offset=c(0.1, 0.1))
# eighth adjustments to mesh, mesh8
mesh <- inla.mesh.2d(gridLocation, 
                     max.edge=c(0.05,0.8),
                     min.angle= 10, #smaller minimum triangle angle
                     cutoff = 0.015, #smaller minimum distance between points
                     offset=c(0.1, 0.1))


# Tips for mesh construction:
# ADjust max.edge and cutoff parameters in the inla.mesh.2d() function
# This will allow you to adjust the density of the mesh vertices. 

# Projector Matrix 
AMat <- inla.spde.make.A(mesh, loc=gridLocation)  # model-fitting locations
AMatCV <- inla.spde.make.A(mesh, loc=CVgridLocation) # validation locations


# FIgure for Mesh and observation locations 
par(mfrow=c(1,1),mar=c(2,2,2,2))
plot(mesh,main="")
points(x=mesh$loc[,1], y=mesh$loc[,2],col="black",pch=16,cex=0.4)
points(x=gridLocation[,1], y=gridLocation[,2],col="blue",pch=16,cex=.8)
points(x=CVgridLocation[,1], y=CVgridLocation[,2],col="red",pch=18,cex=1.2)
mtext("Mesh: INLA",cex=2)
legend("topright", legend=c("Model Fitting" , "Validation" , "Mesh Vertices"),
       col=c("blue","red","black"), pch=c(16,16,16))

#######################################################
# Generate Moran's Basis Functions
# See Hughes and Haran (2012) or Lee and Haran (2020+) for details
wp<-ncol(AMat)
DMat<-diag(apply(mesh$graph$vv,1,sum))
WeightMat<-mesh$graph$vv
PrecMat<-DMat-WeightMat  
Nnew<-nrow(WeightMat)
OrthSpace<-diag(Nnew)-(rep(1,Nnew)%*%t(rep(1,Nnew)))/Nnew
MoransOperator<-OrthSpace%*%(WeightMat%*%OrthSpace)# Moran's Operator
# Moran's Basis functions
MoransOperatorEig<-eigen(MoransOperator) # Eigenvectors of the Moran's Operator
#######################################################

# Save File
save(mesh, AMat, AMatCV, MoransOperatorEig, file="C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR/samples/mydatamesh2.RData")
