##generateSamples, but with my data
rm(list=ls())
library(fields) ; library(mvtnorm) ; library(classInt)
#setwd("~/Dropbox/PICAR/") # Set Directory #ben
setwd("C:/Users/saman/Dropbox/climate and crop research- summer 2020/PICAR/PICAR") #me
source(file = "source/sharedFunctions.R")
load("Data_prism")

df<- Data_prism[Data_prism$year==1,]
df<- df[df$StateANSI==42,]
#add yield increase column to data. yes 1 no 0
df$Yield_increase = ifelse(df$Yield_anomaly > 0, 1, 0)

#I'll convert the latitude and longitudes into grid locations between 0 and 1 for compatibility with picar
max(df$lat); min(df$lon)
s_lat<- df$lat- min(df$lat)
s_lon<- df$lon - max(df$lon)

scale_lat<- s_lat/max(s_lat)
scale_lon<- 1+ s_lon/abs(min(s_lon))

df$sclat<- scale_lat
df$sclon<- scale_lon

# create training set, test set
totlen=dim(df)[1]
cvind<-c(1:totlen)
set.seed(2020)
groups<-split(cvind,sample(rep(1:10,each = floor(totlen/10))))
test_indx<-groups[[1]] #10% of data used as testing data, other 90% training data
test_data<-df[test_indx,]
train_indx<-cvind[-test_indx]
train_data<- df[train_indx,]
df2<- rbind(train_data, test_data)

n= nrow(train_data); modInd= 1:n; cvInd= (n+1):nrow(df2)


# set the grid locations for training, testing data. Combine
gridLocation<- cbind(train_data$sclon, train_data$sclat)
CVgridLocation<- cbind(test_data$sclon, test_data$sclat)
comboLocation<-rbind(gridLocation,CVgridLocation)

distMatMod<- as.matrix(rdist(gridLocation))
distMatCV<- as.matrix(rdist(CVgridLocation))
distMatFull<-as.matrix(rdist(comboLocation))

# Observations
obsFullBin<- df2$Yield_increase
# Observations that we will use to create the model
obsModBin<-obsFullBin[modInd] 
obsCVBin<- obsFullBin[cvInd]

# Plot the data
plotRF(dat=obsFullBin, location = gridLocation , label="Binary Observations")

# Covariates
XMat<-cbind(rep(1, nrow(train_data)), train_data$GDD, train_data$EDD, train_data$Pr, (train_data$Pr)^2, train_data$Tmax, train_data$Tmin, train_data$lastSF, train_data$firstFF)
XMatCV<-cbind(rep(1, nrow(test_data)), test_data$GDD, test_data$EDD, test_data$Pr, (test_data$Pr)^2, test_data$Tmax, test_data$Tmin, test_data$lastSF, test_data$firstFF)


save.image(file="samples/CropPA1981_SpatialData.RData") # Save Data 