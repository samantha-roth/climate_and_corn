##My version of Prism_S&Rmodel.R
setwd("c:/Users/smr31/advising related/climate and crop research- summer 2020/Climate_CornYield-master-new/Climate_CornYield-master")
#remove all data/variables and plots
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(ncdf4)
library(ggplot2)
library(usmap)
library(precintcon)
library(BMS)
library(adaptMCMC)
library(geoR) #for variogram
library(housingData) #for returning county centroid
library(gstat)
library(fields)
library(binaryLogic)
load("Data_prism")
Data_prism<-Data_prism[complete.cases(Data_prism), ]
Data_prism$StateANSI<-factor(Data_prism$StateANSI)
Data_prism$year<-factor(Data_prism$year)
### reimplement S&R paper's model using GDD/EDD (piecewise linear model)


#quadratic yield trend in each state 
statenum<-length(levels(Data_prism$StateANSI))
countynum<-length(levels(Data_prism$fips))
Year<-length(levels(Data_prism$year))
Yield_stateyear<-rep(NA,statenum*Year) # the quadratic yield of each state each year
for (i in 1:statenum){
  index<-which(Data_prism$StateANSI == levels(Data_prism$StateANSI)[i])
  Yield_state<-Data_prism$Yield[index]
  Year_state<-Data_prism$year[index]
  
  newdata<-data.frame(Yield=log(Yield_state),year=as.integer(Year_state),year_sqr=as.integer(Year_state)^2)
  model<-lm(Yield~poly(year,2),data = newdata) #a quadratic time trend fit for each state
  Yield_predict<-predict(model,newdata) 
  for (j in 1:Year){
    index<-which(Year_state == levels(Data_prism$year)[j])
    Yield_stateyear[(i-1)*Year+j]=mean(Yield_predict[index]) #each county each year has a "trend" value
  }
}

trendindex<-(as.integer(Data_prism$StateANSI)-1)*Year+as.integer(Data_prism$year)  #find the index of detrending yield based on state and year
Data_prism$logDTYield<-log(Data_prism$Yield-Yield_stateyear[trendindex])   #this anomaly is: yield - the quadratic trend
##ME: I took the log after subtracting the time trend because the other way felt wrong

model<-lm(logDTYield~ GDD+EDD+poly(Pr,2)+SPI+year+fips,data= Data_prism, weights= Data_prism$Area) #year and fips are fixed effects
summary(model)
hind<-predict(model,Data_prism)
Data_prism$hindcast<-exp(hind+Yield_stateyear[trendindex]) #the hindcasts
Data_prism$residual<-Data_prism$Yield-Data_prism$hindcast  #the residuals

##ME: I'm testing this model without weights corresponding to Area. 
##The above weighting system implies that the variance that 
##since each weight is the reciprocal of the variance of the measurement,
##the error variance decreases as the Area harvested increases? This seems backwards.
##So I'll fit a SLR without weights and see if there is any trend between area and residuals.
model2<-lm(logDTYield~ GDD+EDD+poly(Pr,2)+SPI+year+fips,data= Data_prism) #year and fips are fixed effects
summary(model2)
resids<- model2$resid 
plot(Data_prism$Area, resids)
hind<-predict(model2,Data_prism)
Data_prism$hindcast<-exp(hind+Yield_stateyear[trendindex]) #the hindcasts
Data_prism$residual<-Data_prism$Yield-Data_prism$hindcast  #the residuals

##ME: END MY ADDITIONAL MODEL

#each year's residual plot
for (i in 1:32){
  index<-which(Data_prism$year==i)
  a=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
               data=Data_prism[index,c(12,26)], values = "residual") + labs(title = paste("Yield residual of each county in ",i+1980,sep=""))+ 
    scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-45,45),name="bu/acre")+theme(plot.title = element_text(size=14))
  plot(a)
}
