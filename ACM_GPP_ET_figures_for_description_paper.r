
#calculate_annual_fluxes_timeseries_grid<-function(var_in,area_in,steps_in_year,time) {

#     # loop over each year first
#     for (t in seq(1,(time-steps_in_year),steps_in_year)) {
#          var_in[,,]

#     } # time loop 


#} # end function calculate_annual_fluxes_timeseries_grid

###
## ACM-GPP-ET discription paper figures
### 

###
## Set working directory

setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")
source("./cardamom_functions/load_all_cardamom_functions.r")
setwd("/home/lsmallma/WORK/GREENHOUSE/models/ACM_GPP_ET/")

###
## Read in each of the data sets

## ACM-GPP-ET runs
load("./outputs/calibration_output.RData")             # Calibration dataset
load("./outputs/validation_nowater_output.RData")      # Out of sample SPA, water held at field capacity
load("./outputs/validation_water_output.RData")        # Out of sample SPA, water allowed to vary
load("./outputs/fluxnet_independent_validation.RData") # Independent FLUXNET GPP and ET derived estimates
load("./outputs/global_1x1_degree_2001_2015.RData")
#load("./outputs/global_1x1_degree_2001_2015_co2_plus100.RData")
#load("./outputs/global_1x1_degree_2001_2015_Tair_plus1.RData")
#load("./outputs/global_1x1_degree_2001_2015_NUE_half.RData")
#load("./outputs/global_1x1_degree_2001_2015_NUE_half_co2_plus100.RData")
#load("./outputs/global_1x1_degree_2001_2015_NUE_half_Tair_plus1.RData")

## SPA calibration / validation sets - what is the impact of dynamic water stocks
nowater = read.csv("/home/lsmallma/gcel/ACM_GPP_ET_RECALIBRATION/output_files/acm_recal_with_spa_200pixels_continuous_timeseries_obs_whole_unfiltered_iWUE_trunk_nowater_copy.csv") 
water = read.csv("/home/lsmallma/gcel/ACM_GPP_ET_RECALIBRATION/output_files/acm_recal_with_spa_200pixels_continuous_timeseries_obs_whole_unfiltered_iWUE_trunk_water_copy.csv")

# Check that all drivers are in agreement
par(mfrow=c(3,4))
plot(water$lai - nowater$lai, main="LAI")
plot(water$avgN - nowater$avgN, main="avgN")
plot(water$co2 - nowater$co2, main="Mean CO2")
plot(water$sat_avg - nowater$sat_avg, main="Mean temperature")
plot(water$vpd_avg - nowater$vpd_avg, main="Mean VPD")
plot(water$swrad_avg - nowater$swrad_avg, main="Mean shortwave radiation")
plot(water$wind_avg - nowater$wind_avg, main="Mean wind")
plot(water$ppt_avg - nowater$ppt_avg, main="Mean rainfall")
plot(water$sfc_pressure - nowater$sfc_pressure, main="Surface pressure")
plot(water$porosity - nowater$porosity, main="Porosity")
par(mfrow=c(3,4))
# Check other fluxes / model diagnositics
plot(water$ground_heat ~ nowater$ground_heat, main="Ground heat") ; abline(0,1,col="red")
plot(water$sensible ~ nowater$sensible, main="Sensible") ; abline(0,1,col="red")
plot(water$energy_balance_residual ~ nowater$energy_balance_residual, main="Energy balance residual") ; abline(0,1,col="red")
plot(water$netrad ~ nowater$netrad, main="Netrad") ; abline(0,1,col="red")
plot(water$canopy_temp ~ nowater$canopy_temp, main="Canopy temperature") ; abline(0,1,col="red")
plot(water$sun_canopy_temp ~ nowater$sun_canopy_temp, main="Sun temperature") ; abline(0,1,col="red")
plot(water$shade_canopy_temp ~ nowater$shade_canopy_temp, main="Shade temperature") ; abline(0,1,col="red")
plot(water$SWP ~ nowater$SWP, main="wSWP") ; abline(0,1,col="red")
plot(water$SWC ~ nowater$SWC, main="SWC") ; abline(0,1,col="red")
plot(water$soil_conductance ~ nowater$soil_conductance, main="Soil conductance") ; abline(0,1,col="red")
# Check GPP delta against drivers
SPA_deltaGPP = water$GPP - nowater$GPP
SPA_deltaET = water$Evap - nowater$Evap
par(mfrow=c(2,4))
plot(SPA_deltaGPP ~ nowater$sat_avg)
plot(SPA_deltaGPP ~ nowater$ppt_avg)
plot(SPA_deltaGPP ~ nowater$vpd_avg)
plot(SPA_deltaGPP ~ nowater$swrad_avg)
plot(SPA_deltaET ~ nowater$sat_avg)
plot(SPA_deltaET ~ nowater$ppt_avg)
plot(SPA_deltaET ~ nowater$vpd_avg)
plot(SPA_deltaET ~ nowater$swrad_avg)

###
## First what did SPA do in response to water change

# what is the mean GPP for fixed and dynamics SPA simulations?
mean(water$GPP) ; mean(nowater$GPP) # water = 4.49 ; nowater = 4.58 gC/m2/day (or a difference of -1.886239 %)
summary(lm(water$GPP~nowater$GPP))$adj.r.squared # 0.9759
# above field capacity 
mean(water$GPP[which(water$SWP > nowater$SWP)]) ; mean(nowater$GPP[which(water$SWP > nowater$SWP)]) # water = 5.812257 ; nowater = 5.812182 gC/m2/day (or a difference of 0.00129 %)
summary(lm(water$GPP[which(water$SWP > nowater$SWP)]~nowater$GPP[which(water$SWP > nowater$SWP)]))$adj.r.squared # 0.9999117
# below field capacity
mean(water$GPP[which(water$SWP < nowater$SWP)]) ; mean(nowater$GPP[which(water$SWP < nowater$SWP)]) # water = 3.7896 ; nowater = 4.082068 gC/m2/day (or a difference of -7.164702 %)
summary(lm(water$GPP[which(water$SWP < nowater$SWP)]~nowater$GPP[which(water$SWP < nowater$SWP)]))$adj.r.squared # 0.8792734

# estimate rmse of difference between SPA simulations
rmse(water$GPP,nowater$GPP) #  0.6651699 gC/m2/day
# now for when above field capacity
rmse(water$GPP[which(water$SWP > nowater$SWP)],nowater$GPP[which(water$SWP > nowater$SWP)]) # 0.08 gC/m2/day
# now for when below field capacity
rmse(water$GPP[which(water$SWP < nowater$SWP)],nowater$GPP[which(water$SWP < nowater$SWP)]) # 1.221998 gC/m2/day
# estimate mean bias across whole database
mean(water$GPP-nowater$GPP) # -0.08642161 gC/m2/day
# now for when above field capacity
mean(water$GPP[which(water$SWP > nowater$SWP)]-nowater$GPP[which(water$SWP > nowater$SWP)]) # 0.0000746537 gC/m2/day
# now for when below field capacity
mean(water$GPP[which(water$SWP < nowater$SWP)]-nowater$GPP[which(water$SWP < nowater$SWP)]) # -0.2924683 gC/m2/day

# how much of the time is spend above field capacity?
length(which(water$SWP > nowater$SWP)) / length(nowater$SWP)  # 58 %
# how much of the time is spend below field capacity?
length(which(water$SWP < nowater$SWP)) / length(nowater$SWP)  # 29 %
# how much of the time is at field capacity?
length(which(water$SWP == nowater$SWP)) / length(nowater$SWP) # 12 %

## Now repeat R2, rmse and bias for each site with calibration / validation SPA set

simulated_pixels = paste(water$lat,water$long,sep="")
simulated_pixels = unique(simulated_pixels)
combined_var = paste(water$lat,water$long,sep="")

SPA_GPP_r2 = array(NA, dim=c(length(simulated_pixels)))
SPA_GPP_rmse = array(NA, dim=c(length(simulated_pixels)))
SPA_GPP_bias = array(NA, dim=c(length(simulated_pixels)))
SPA_ET_r2 = array(NA, dim=c(length(simulated_pixels)))
SPA_ET_rmse = array(NA, dim=c(length(simulated_pixels)))
SPA_ET_bias = array(NA, dim=c(length(simulated_pixels)))
SPA_WUE_r2 = array(NA, dim=c(length(simulated_pixels)))
SPA_WUE_rmse = array(NA, dim=c(length(simulated_pixels)))
SPA_WUE_bias = array(NA, dim=c(length(simulated_pixels)))

# Loop through each site now
for (i in seq(1, length(simulated_pixels))){
     # determine where in the dataset all the time points are for a particular site
     how_many = which(combined_var == simulated_pixels[i])
     # calculate R2, rmse and bias for each site
     # GPP
     SPA_GPP_r2[i] = summary(lm(water$GPP[how_many]~nowater$GPP[how_many]))$adj.r.squared
     SPA_GPP_rmse[i] = rmse(water$GPP[how_many],nowater$GPP[how_many])
     SPA_GPP_bias[i] = mean(water$GPP[how_many]-nowater$GPP[how_many])
     # Evapo-Transpiration
     SPA_ET_r2[i] = summary(lm(water$Evap[how_many]~nowater$Evap[how_many]))$adj.r.squared
     SPA_ET_rmse[i] = rmse(water$Evap[how_many],nowater$Evap[how_many])
     SPA_ET_bias[i] = mean(water$Evap[how_many]-nowater$Evap[how_many])
     # WUE
     water_trans = water$Evap[how_many] - water$soilevap[how_many] - water$wetevap[how_many]
     nowater_trans = nowater$Evap[how_many] - nowater$soilevap[how_many] - nowater$wetevap[how_many]
     filter = which(water$GPP[how_many] > 0 & nowater$GPP[how_many] & water_trans > 0 & nowater_trans > 0)
     SPA_WUE_r2[i] = summary(lm(as.vector(water$GPP[how_many[filter]]/water_trans[filter])~as.vector(nowater$GPP[how_many[filter]]/nowater_trans[filter])))$adj.r.squared
     SPA_WUE_rmse[i] = rmse(as.vector(water$GPP[how_many[filter]]/water_trans[filter]),as.vector(nowater$GPP[how_many[filter]]/nowater_trans[filter]))
     SPA_WUE_bias[i] = mean(as.vector(water$GPP[how_many[filter]]/water_trans[filter])-as.vector(nowater$GPP[how_many[filter]]/nowater_trans[filter]))
}
print("SPA water-nowater output: ")
print(paste("    GPP R2   = ",round(mean(SPA_GPP_r2,na.rm=TRUE),digits=3)," ",sep=""))
print(paste("    GPP RMSE = ",round(mean(SPA_GPP_rmse,na.rm=TRUE),digits=3)," gC/m2/day",sep=""))
print(paste("    GPP BIAS = ",round(mean(SPA_GPP_bias,na.rm=TRUE),digits=3)," gC/m2/day",sep=""))
print(paste("     ET R2   = ",round(mean(SPA_ET_r2,na.rm=TRUE),digits=3)," ",sep=""))
print(paste("     ET RMSE = ",round(mean(SPA_ET_rmse,na.rm=TRUE),digits=3)," kgH2O/m2/day",sep=""))
print(paste("     ET BIAS = ",round(mean(SPA_ET_bias,na.rm=TRUE),digits=3)," kgH2O/m2/day",sep=""))
print(paste("    WUE R2   = ",round(mean(SPA_WUE_r2,na.rm=TRUE),digits=3)," ",sep=""))
print(paste("    WUE RMSE = ",round(mean(SPA_WUE_rmse,na.rm=TRUE),digits=3)," kgH2O/m2/day",sep=""))
print(paste("    WUE BIAS = ",round(mean(SPA_WUE_bias,na.rm=TRUE),digits=3)," kgH2O/m2/day",sep=""))

my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral"))))
#fig_height=6000 ; fig_width=4000
#jpeg(file="./FIGURES/Cal_val_paper_figure_3_heat_map.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(4,2), mar=c(4.4, 4.2, 3.4, 2), omi=c(0.2, 0.2, 0.2, 0.40))
var1 = nowater$GPP
var2 = water$GPP
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,
              main=expression(paste("GPP"," (gC ",m^-2," da",y^-1,")")), ylab="", xlab="",cex=0.5,pch=16,cex.axis=2.0,cex.lab=1.6,cex.main=2.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0055,diff(range(var2,na.rm=TRUE))*0.0055),nbin=128*20) ; abline(0,1,col="red",lwd=3)
#text(0,max(var2*0.82,na.rm=TRUE),labels=bquote(RMSE == .(round(mean(fluxnet_validation_output$gpp_rmse,na.rm=TRUE),2))), cex=1.5, pos=4)
#text(0,max(var2*0.72,na.rm=TRUE),labels=bquote(Bias == .(round(mean(fluxnet_validation_output$gpp_bias,na.rm=TRUE),2))), cex=1.5, pos=4)
#text(0,max(var2*0.92,na.rm=TRUE),labels=bquote(R^2 == .(round(mean(fluxnet_daily_GPP_r2,na.rm=TRUE),2))), cex=1.5, pos=4)
mtext(expression(paste("SPA - dynamic soil water")),side = 2,cex=2.0, padj = -1.5, adj = 0.5)
mtext(expression(paste("SPA - fixed soil water")),side = 1,cex=2.0, padj = 1.75, adj = 1.7)
var1 = nowater$Evap
var2 = water$Evap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main=expression(paste("Evapo-transpiration"," (kg",H[2],"O ",m^-2," da",y^-1,")")), ylab="", xlab="",
              cex=0.5,pch=16,cex.axis=2.0,cex.lab=1.6,cex.main=2.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0055,diff(range(var2,na.rm=TRUE))*0.0055),nbin=128*20) ; abline(0,1,col="red",lwd=3)
#text(max(var1*0.82,na.rm=TRUE),max(var2*0.82,na.rm=TRUE),labels=bquote(RMSE == .(round(mean(fluxnet_validation_output$ET_rmse,na.rm=TRUE),2))), cex=1.5, pos=4)
#text(max(var1*0.72,na.rm=TRUE),max(var2*0.72,na.rm=TRUE),labels=bquote(Bias == .(round(mean(fluxnet_validation_output$ET_bias,na.rm=TRUE),2))), cex=1.5, pos=4)
#text(max(var1*0.92,na.rm=TRUE),max(var2*0.92,na.rm=TRUE),labels=bquote(R^2 == .(round(mean(fluxnet_daily_ET_r2,na.rm=TRUE),2))), cex=1.5, pos=4)

hist(SPA_GPP_r2, main="", xlab=expression(paste("Site specific GPP ",R^2,"",sep="")), ylab="", cex.lab = 2.3,cex.axis = 2, cex=2, col="grey") 
hist(SPA_ET_r2, main="", xlab=expression(paste("Site specific ET ",R^2,"",sep="")), ylab="", cex.lab = 2.3,cex.axis = 2, cex=2, col="grey") 
hist(SPA_GPP_rmse, main="", xlab=expression(paste("Site specific GPP RMSE (gC/",m^2,"/day)",sep="")), ylab="", cex.lab = 2.3,cex.axis = 2, cex=2, col="grey")
hist(SPA_ET_rmse, main="", xlab=expression(paste("Site specific ET RMSE (kg",H[2],"O/",m^2,"/day)",sep="")), ylab="", cex.lab = 2.3,cex.axis = 2, cex=2, col="grey")
hist(SPA_GPP_bias, main="", xlab=expression(paste("Site specific GPP Bias (gC/",m^2,"/day)",sep="")), ylab="", cex.lab = 2.3,cex.axis = 2, cex=2, col="grey")
hist(SPA_ET_bias, main="", xlab=expression(paste("Site specific ET Bias (kg",H[2],"O/",m^2,"/day)",sep="")), ylab="", cex.lab = 2.3,cex.axis = 2, cex=2, col="grey")
#dev.off()

###
## What did ACM do in response to change in water

# what is the mean GPP for fixed and dynamics ACM simulations?
mean(validation_water_output$mean_gpp) ; mean(validation_nowater_output$mean_gpp) # water = 4.415887 ; nowater = 4.465642 gC/m2/day (or a difference of -1.114174 %)
summary(lm(validation_water_output$mean_gpp~validation_nowater_output$mean_gpp))$adj.r.squared # 0.9836108
# above field capacity 
mean(validation_water_output$mean_gpp[which(validation_water_output$mean_wSWP > validation_nowater_output$mean_wSWP)]) ; mean(validation_nowater_output$mean_gpp[which(validation_water_output$mean_wSWP > validation_nowater_output$mean_wSWP)]) # water = 4.98649 ; nowater = 4.985598 gC/m2/day (or a difference of 0.01789153 %)
summary(lm(validation_water_output$mean_gpp[which(validation_water_output$mean_wSWP > validation_nowater_output$mean_wSWP)]~validation_nowater_output$mean_gpp[which(validation_water_output$mean_wSWP > validation_nowater_output$mean_wSWP)]))$adj.r.squared # 0.9994
# below field capacity
mean(validation_water_output$mean_gpp[which(validation_water_output$mean_wSWP < validation_nowater_output$mean_wSWP)]) ; mean(validation_nowater_output$mean_gpp[which(validation_water_output$mean_wSWP < validation_nowater_output$mean_wSWP)]) # water = 3.587599 ; nowater = 3.710879 gC/m2/day (or a difference of -3.322124 %)
summary(lm(validation_water_output$mean_gpp[which(validation_water_output$mean_wSWP < validation_nowater_output$mean_wSWP)]~validation_nowater_output$mean_gpp[which(validation_water_output$mean_wSWP < validation_nowater_output$mean_wSWP)]))$adj.r.squared # 0.9507061

# estimate rmse of difference between ACM simulations
rmse(validation_water_output$mean_gpp,validation_nowater_output$mean_gpp) #  0.51470562 gC/m2/day
# now for when above field capacity
rmse(validation_water_output$mean_gpp[which(validation_water_output$mean_wSWP > validation_nowater_output$mean_wSWP)],validation_nowater_output$mean_gpp[which(validation_water_output$mean_wSWP > validation_nowater_output$mean_wSWP)]) # 0.004212332 gC/m2/day
# now for when below field capacity
rmse(validation_water_output$mean_gpp[which(validation_water_output$mean_wSWP < validation_nowater_output$mean_wSWP)],validation_nowater_output$mean_gpp[which(validation_water_output$mean_wSWP < validation_nowater_output$mean_wSWP)]) # 0.8565555 gC/m2/day
# estimate mean bias across whole database
mean(validation_water_output$mean_gpp-validation_nowater_output$mean_gpp) #  -0.04975568 gC/m2/day
# now for when above field capacity
mean(validation_water_output$mean_gpp[which(validation_water_output$mean_wSWP > validation_nowater_output$mean_wSWP)]-validation_nowater_output$mean_gpp[which(validation_water_output$mean_wSWP > validation_nowater_output$mean_wSWP)]) #  0.0008921788 gC/m2/day
# now for when below field capacity
mean(validation_water_output$mean_gpp[which(validation_water_output$mean_wSWP < validation_nowater_output$mean_wSWP)]-validation_nowater_output$mean_gpp[which(validation_water_output$mean_wSWP < validation_nowater_output$mean_wSWP)]) # -0.1232798 gC/m2/day

# how much of the time is spend above field capacity?
length(which(validation_water_output$mean_wSWP > validation_nowater_output$mean_wSWP)) / length(validation_nowater_output$mean_wSWP)  # 60 %
# how much of the time is spend below field capacity?
length(which(validation_water_output$mean_wSWP < validation_nowater_output$mean_wSWP)) / length(validation_nowater_output$mean_wSWP)  # 40 %
# how much of the time is at field capacity?
length(which(validation_water_output$mean_wSWP == validation_nowater_output$mean_wSWP)) / length(validation_nowater_output$mean_wSWP) # ~0 %

## Now repeat R2, rmse and bias for each site with calibration / validation SPA set

simulated_pixels = paste(validation_water_output$drivers$lat,validation_water_output$drivers$long,sep="")
simulated_pixels = unique(simulated_pixels)
combined_var = paste(validation_water_output$drivers$lat,validation_water_output$drivers$long,sep="")

ACM_GPP_r2 = array(NA, dim=c(length(simulated_pixels)))
ACM_GPP_rmse = array(NA, dim=c(length(simulated_pixels)))
ACM_GPP_bias = array(NA, dim=c(length(simulated_pixels)))
ACM_ET_r2 = array(NA, dim=c(length(simulated_pixels)))
ACM_ET_rmse = array(NA, dim=c(length(simulated_pixels)))
ACM_ET_bias = array(NA, dim=c(length(simulated_pixels)))
ACM_WUE_r2 = array(NA, dim=c(length(simulated_pixels)))
ACM_WUE_rmse = array(NA, dim=c(length(simulated_pixels)))
ACM_WUE_bias = array(NA, dim=c(length(simulated_pixels)))

# Loop through each site now
for (i in seq(1, length(simulated_pixels))){
     # determine where in the dataset all the time points are for a particular site
     how_many = which(combined_var == simulated_pixels[i])
     # calculate R2, rmse and bias for each site
     # GPP
     ACM_GPP_r2[i] = summary(lm(validation_water_output$mean_gpp[how_many]~validation_nowater_output$mean_gpp[how_many]))$adj.r.squared
     ACM_GPP_rmse[i] = rmse(validation_water_output$mean_gpp[how_many],validation_nowater_output$mean_gpp[how_many])
     ACM_GPP_bias[i] = mean(validation_water_output$mean_gpp[how_many]-validation_nowater_output$mean_gpp[how_many])
     # Evapo-Transpiration
     water_Evap = validation_water_output$mean_transpiration + validation_water_output$mean_wetcanopyevap + validation_water_output$mean_soilevaporation
     nowater_Evap = validation_nowater_output$mean_transpiration + validation_nowater_output$mean_wetcanopyevap + validation_nowater_output$mean_soilevaporation
     ACM_ET_r2[i] = summary(lm(water_Evap[how_many]~nowater_Evap[how_many]))$adj.r.squared
     ACM_ET_rmse[i] = rmse(water_Evap[how_many],nowater_Evap[how_many])
     ACM_ET_bias[i] = mean(water_Evap[how_many]-nowater_Evap[how_many])
     # WUE
     water_trans = validation_water_output$mean_transpiration[how_many]
     nowater_trans = validation_nowater_output$mean_transpiration[how_many]
     filter = which(water_trans > 0 & nowater_trans > 0) ; how_many = how_many[filter]
     water_trans = water_trans[filter] ; nowater_trans = nowater_trans[filter]
     ACM_WUE_r2[i] = summary(lm(as.vector(validation_water_output$mean_gpp[how_many]/water_trans)~as.vector(validation_nowater_output$mean_gpp[how_many]/nowater_trans)))$adj.r.squared
     ACM_WUE_rmse[i] = rmse(as.vector(validation_water_output$mean_gpp[how_many]/water_trans),as.vector(validation_nowater_output$mean_gpp[how_many]/nowater_trans))
     ACM_WUE_bias[i] = mean(as.vector(validation_water_output$mean_gpp[how_many]/water_trans)-as.vector(validation_nowater_output$mean_gpp[how_many]/nowater_trans))
}
print("ACM water-nowater output: ")
print(paste("    GPP R2   = ",round(mean(ACM_GPP_r2,na.rm=TRUE),digits=3)," ",sep=""))
print(paste("    GPP RMSE = ",round(mean(ACM_GPP_rmse,na.rm=TRUE),digits=3)," gC/m2/day",sep=""))
print(paste("    GPP BIAS = ",round(mean(ACM_GPP_bias,na.rm=TRUE),digits=3)," gC/m2/day",sep=""))
print(paste("     ET R2   = ",round(mean(ACM_ET_r2,na.rm=TRUE),digits=3)," ",sep=""))
print(paste("     ET RMSE = ",round(mean(ACM_ET_rmse,na.rm=TRUE),digits=3)," kgH2O/m2/day",sep=""))
print(paste("     ET BIAS = ",round(mean(ACM_ET_bias,na.rm=TRUE),digits=3)," kgH2O/m2/day",sep=""))
print(paste("    WUE R2   = ",round(mean(ACM_WUE_r2,na.rm=TRUE),digits=3)," ",sep=""))
print(paste("    WUE RMSE = ",round(mean(ACM_WUE_rmse,na.rm=TRUE),digits=3)," kgH2O/m2/day",sep=""))
print(paste("    WUE BIAS = ",round(mean(ACM_WUE_bias,na.rm=TRUE),digits=3)," kgH2O/m2/day",sep=""))

#ACM
#[1] "    GPP R2   = 0.95 "
#[1] "    GPP RMSE = 0.206 gC/m2/day"
#[1] "    GPP BIAS = -0.086 gC/m2/day"
#[1] "     ET R2   = 0.897 "
#[1] "     ET RMSE = 0.185 kgH2O/m2/day"
#[1] "     ET BIAS = -0.009 kgH2O/m2/day"
#[1] "    WUE R2   = 0.543 "
#[1] "    WUE RMSE = 0.278 kgH2O/m2/day"
#[1] "    WUE BIAS = -0.033 kgH2O/m2/day"

#SPA
#[1] "    GPP R2   = 0.97 "
#[1] "    GPP RMSE = 0.107 gC/m2/day"
#[1] "    GPP BIAS = -0.050 gC/m2/day"
#[1] "     ET R2   = 0.958 "
#[1] "     ET RMSE = 0.139 kgH2O/m2/day"
#[1] "     ET BIAS = -0.018 kgH2O/m2/day"
#[1] "    WUE R2   = 0.885 "
#[1] "    WUE RMSE = 0.139 kgH2O/m2/day"
#[1] "    WUE BIAS = -0.018 kgH2O/m2/day"

## CARDAMOM run used to provide LAI and root C
cardamom = nc_open("/disk/scratch/local.2/lsmallma/Forest2020/C_cycle_analyses/DALEC_GSI_DFOL_CWD_FR_1_2001_2015_NEE_GPP_Rh_Ra_Bio_lit_cwd_som_timeseries.nc")
cardamom_latitude = ncvar_get(cardamom,"lat") ; cardamom_longitude = ncvar_get(cardamom,"long")

# ## Read in FLUXCOM GPP estimates for independent validation
# list_of_files = list.files("/home/lsmallma/gcel/FLUXCOM/", full.name=TRUE)
# list_of_files = list_of_files[grepl("GPP.annual",list_of_files)]
# years_to_do = 2001:2013
# GPP_types = c("ANN","MARS","RF")
# for (y in seq(1,length(years_to_do))) {
#     for (m in seq(1,length(GPP_types))) {
#          tmp_list = list_of_files[grepl(years_to_do[y], list_of_files)]
#          tmp_list = tmp_list[grepl(GPP_types[m], tmp_list)]
#          datain = nc_open(tmp_list)
#          tmp = ncvar_get(datain, "GPP")
#          if (y == 1 & m == 1) {
#              fluxcom_gpp = array(0, dim=c(dim(tmp),length(years_to_do)))
#              fluxcom_gpp[,,y] = tmp * (1/length(GPP_types))
#          } else {
#              fluxcom_gpp[,,y] = fluxcom_gpp[,,y] + (tmp * (1/length(GPP_types)))
#          }
#     } # methods to loop
# } # years to loop
# 
# ## Updates for comparison later on with ACM-GPP-ET
# # Convert units from gC/m2/day -> gC/m2/yr
# fluxcom_gpp = fluxcom_gpp * 365.25
# # Filter out zeros to NA
# fluxcom_gpp[which(fluxcom_gpp == 0)] = NA
# # Fix latitude orientation to match ACM-GPP-ET
# fluxcom_gpp = fluxcom_gpp[,dim(fluxcom_gpp)[2]:1,]
# 
# ## Regrid to 1 x 1 degree to match with ACM-GPP-ET
# tmp_out_array = array(NA, dim=c(dim(global_output$mean_gpp),length(years_to_do)))
# target_grid = raster(global_output$mean_gpp) 
# for (y in seq(1, length(years_to_do))) {
#      tmp_in = raster(fluxcom_gpp[,,y])    
#      tmp_out = resample(tmp_in, target_grid, method="bilinear", filename="")
#      tmp_out_array[,,y] = as.vector(t(tmp_out))
# }
# # now replace the original fluxcom at 0.5 degree to 1 degree
# fluxcom_gpp = tmp_out_array ; rm(tmp_in,tmp_out,tmp_out_array,target_grid)
# fluxcom_gpp_total = array(0,dim=dim(fluxcom_gpp)[1:2])
# for (t in seq(1,dim(fluxcom_gpp)[3])) {fluxcom_gpp_total = fluxcom_gpp_total + fluxcom_gpp[,,t]} ; fluxcom_gpp_total = fluxcom_gpp_total / dim(fluxcom_gpp)[3]

###
## Some statistical / summary values
###

# print("FLUXCOM output: ")
# print(paste("       GPP = ",round(sum(fluxcom_gpp_total*global_output_NUE_half$area*1e-15,na.rm=TRUE),digits=1)," PgC",sep=""))

print("Calibration output: ")
spa_transpiration = (calibration_output$drivers$Evap - calibration_output$drivers$soilevap - calibration_output$drivers$wetevap)
acm_et = calibration_output$mean_transpiration + calibration_output$mean_soilevaporation + calibration_output$mean_wetcanopyevap
tmp = which(calibration_output$drivers$GPP > 0.1 & calibration_output$mean_gpp > 0.1 & calibration_output$mean_transpiration > 0.1 & spa_transpiration > 0.1 & calibration_output$drivers$lai > 0.1)
acm_wue = sum(calibration_output$mean_gpp[tmp]) / sum(calibration_output$mean_transpiration[tmp])
spa_wue = sum(calibration_output$drivers$GPP[tmp]) / sum(spa_transpiration[tmp])
print(paste("       GPP R2     = ",round(calibration_output$gpp_r2,digits=3),sep="")) 
print(paste("       GPP Bias   = ",round(calibration_output$gpp_bias,digits=3),sep=""))
print(paste("       GPP RMSE   = ",round(calibration_output$gpp_rmse,digits=3),sep=""))
print(paste("       Trans R2   = ",round(calibration_output$transpiration_r2,digits=3),sep=""))
print(paste("       Trans Bias = ",round(calibration_output$transpiration_bias,digits=3),sep=""))
print(paste("       Trans RMSE = ",round(calibration_output$transpiration_rmse,digits=3),sep=""))
print(paste("       Esoil R2   = ",round(calibration_output$soilevaporation_r2,digits=3),sep=""))
print(paste("       Esoil Bias = ",round(calibration_output$soilevaporation_bias,digits=3),sep=""))
print(paste("       Esoil RMSE = ",round(calibration_output$soilevaporation_rmse,digits=3),sep=""))
print(paste("       Ewet R2    = ",round(calibration_output$wetcanopyevap_r2,digits=3),sep=""))
print(paste("       Ewet Bias  = ",round(calibration_output$wetcanopyevap_bias,digits=3),sep=""))
print(paste("       Ewet RMSE  = ",round(calibration_output$wetcanopyevap_rmse,digits=3),sep=""))
print(paste("       SPA WUE = ",round(mean(spa_wue),digits=2)," gC/kgH2O",sep=""))
print(paste("       ACM WUE = ",round(mean(acm_wue),digits=2)," gC/kgH2O",sep=""))
print(paste("       WUE R2 = ",round(summary(lm(spa_wue ~ acm_wue))$adj.r.squared,digits=3)," ",sep=""))
print(paste("       WUE BIAS = ",round(mean(acm_wue-spa_wue,na.rm=TRUE),digits=3)," gC/kgH2O",sep=""))
print(paste("       WUE RMSE = ",round(rmse(acm_wue,spa_wue),digits=3)," gC/kgH2O",sep=""))
print(paste("       SPA T/ET   = ",round(sum(spa_transpiration)/sum(calibration_output$drivers$Evap),digits=3),sep=""))
print(paste("       SPA S/ET   = ",round(sum(calibration_output$drivers$soilevap)/sum(calibration_output$drivers$Evap),digits=3),sep=""))
print(paste("       SPA W/ET   = ",round(sum(calibration_output$drivers$wetevap)/sum(calibration_output$drivers$Evap),digits=3),sep=""))
print(paste("       ACM T/ET   = ",round(sum(calibration_output$mean_transpiration)/sum(acm_et),digits=3),sep=""))
print(paste("       ACM S/ET   = ",round(sum(calibration_output$mean_soilevaporation)/sum(acm_et),digits=3),sep=""))
print(paste("       ACM W/ET   = ",round(sum(calibration_output$mean_wetcanopyevap)/sum(acm_et),digits=3),sep=""))

print("Validation output (nowater): ")
spa_transpiration = (validation_nowater_output$drivers$Evap - validation_nowater_output$drivers$soilevap - validation_nowater_output$drivers$wetevap)
acm_et = validation_nowater_output$mean_transpiration + validation_nowater_output$mean_soilevaporation + validation_nowater_output$mean_wetcanopyevap
tmp = which(validation_nowater_output$drivers$GPP > 0.1 & spa_transpiration > 0.1 & validation_nowater_output$mean_gpp > 0.1 & validation_nowater_output$mean_transpiration > 0.1 & validation_nowater_output$drivers$lai > 0.1)
acm_wue = validation_nowater_output$mean_gpp / validation_nowater_output$mean_transpiration
spa_wue = validation_nowater_output$drivers$GPP / spa_transpiration
print(paste("       GPP Mean   = ",round(mean(as.vector(validation_nowater_output$mean_gpp)),digits=2),sep=""))
print(paste("       Trans Mean = ",round(mean(as.vector(validation_nowater_output$mean_transpiration)),digits=2),sep=""))
print(paste("       Esoil Mean = ",round(mean(as.vector(validation_nowater_output$mean_soilevaporation)),digits=2),sep=""))
print(paste("       Ewet Mean  = ",round(mean(as.vector(validation_nowater_output$mean_wetcanopyevap)),digits=2),sep=""))
print(paste("-------------------","-",sep=""))
print(paste("       GPP R2     = ",round(validation_nowater_output$gpp_r2,digits=2),sep=""))
print(paste("       GPP Bias   = ",round(validation_nowater_output$gpp_bias,digits=2),sep=""))
print(paste("       GPP RMSE   = ",round(validation_nowater_output$gpp_rmse,digits=2),sep=""))
print(paste("       Trans R2   = ",round(validation_nowater_output$transpiration_r2,digits=2),sep=""))
print(paste("       Trans Bias = ",round(validation_nowater_output$transpiration_bias,digits=2),sep=""))
print(paste("       Trans RMSE = ",round(validation_nowater_output$transpiration_rmse,digits=2),sep=""))
print(paste("       Esoil R2   = ",round(validation_nowater_output$soilevaporation_r2,digits=2),sep=""))
print(paste("       Esoil Bias = ",round(validation_nowater_output$soilevaporation_bias,digits=2),sep=""))
print(paste("       Esoil RMSE = ",round(validation_nowater_output$soilevaporation_rmse,digits=2),sep=""))
print(paste("       Ewet R2    = ",round(validation_nowater_output$wetcanopyevap_r2,digits=2),sep=""))
print(paste("       Ewet Bias  = ",round(validation_nowater_output$wetcanopyevap_bias,digits=2),sep=""))
print(paste("       Ewet RMSE  = ",round(validation_nowater_output$wetcanopyevap_rmse,digits=2),sep=""))
print(paste("       SPA WUE = ",round(mean(spa_wue[tmp],na.rm=TRUE),digits=2)," gC/kgH2O",sep=""))
print(paste("       ACM WUE = ",round(mean(acm_wue[tmp],na.rm=TRUE),digits=2)," gC/kgH2O",sep=""))
print(paste("       WUE R2     = ",round(summary(lm(spa_wue[tmp] ~ acm_wue[tmp]))$adj.r.squared,digits=2)," ",sep=""))
print(paste("       WUE BIAS   = ",round(mean(acm_wue[tmp]-spa_wue[tmp],na.rm=TRUE),digits=2)," gC/kgH2O",sep=""))
print(paste("       WUE BIAS = ",round(mean(acm_wue[tmp]-spa_wue[tmp],na.rm=TRUE),digits=2)," gC/kgH2O",sep=""))
print(paste("       SPA T/ET   = ",round(sum(spa_transpiration)/sum(validation_nowater_output$drivers$Evap),digits=2),sep=""))
print(paste("       SPA S/ET   = ",round(sum(validation_nowater_output$drivers$soilevap)/sum(validation_nowater_output$drivers$Evap),digits=2),sep=""))
print(paste("       SPA W/ET   = ",round(sum(validation_nowater_output$drivers$wetevap)/sum(validation_nowater_output$drivers$Evap),digits=2),sep=""))
print(paste("       ACM T/ET   = ",round(sum(validation_nowater_output$mean_transpiration)/sum(acm_et),digits=2),sep=""))
print(paste("       ACM S/ET   = ",round(sum(validation_nowater_output$mean_soilevaporation)/sum(acm_et),digits=2),sep=""))
print(paste("       ACM W/ET   = ",round(sum(validation_nowater_output$mean_wetcanopyevap)/sum(acm_et),digits=2),sep=""))

print("Validation output (water): ")
SPA_deltaGPP = validation_water_output$drivers$GPP - validation_nowater_output$drivers$GPP
SPA_deltaET = validation_water_output$drivers$Evap - validation_nowater_output$drivers$Evap

ACM_deltaET = validation_water_output$mean_transpiration + validation_water_output$mean_soilevaporation + validation_water_output$mean_wetcanopyevap
ACM_deltaET = ACM_deltaET - (validation_nowater_output$mean_transpiration + validation_nowater_output$mean_soilevaporation + validation_nowater_output$mean_wetcanopyevap)
ACM_deltaGPP = validation_water_output$mean_gpp - validation_nowater_output$mean_gpp
ACM_deltaGPP_r2 = summary(lm(SPA_deltaGPP ~ ACM_deltaGPP))$adj.r.squared 
ACM_deltaET_r2 = summary(lm(SPA_deltaET ~ ACM_deltaET))$adj.r.squared
par(mfrow=c(2,2))
plot(SPA_deltaGPP ~ ACM_deltaGPP) ; abline(0,1,col="red") ; plot(SPA_deltaET ~ ACM_deltaET) ; abline(0,1,col="red")

#summary(lm(SPA_deltaGPP ~ water$sat_avg))$adj.r.squared ; summary(lm(ACM_deltaGPP ~ water$sat_avg))$adj.r.squared 
#summary(lm(SPA_deltaGPP ~ water$swrad_avg))$adj.r.squared ; summary(lm(ACM_deltaGPP ~ water$swrad_avg))$adj.r.squared 
#summary(lm(SPA_deltaGPP ~ water$ppt_avg))$adj.r.squared ; summary(lm(ACM_deltaGPP ~ water$ppt_avg))$adj.r.squared 
#summary(lm(SPA_deltaGPP ~ water$vpd_avg))$adj.r.squared ; summary(lm(ACM_deltaGPP ~ water$vpd_avg))$adj.r.squared

#summary(lm(SPA_deltaET ~ water$sat_avg))$adj.r.squared ; summary(lm(ACM_deltaET ~ water$sat_avg))$adj.r.squared 
#summary(lm(SPA_deltaET ~ water$swrad_avg))$adj.r.squared ; summary(lm(ACM_deltaET ~ water$swrad_avg))$adj.r.squared 
#summary(lm(SPA_deltaET ~ water$ppt_avg))$adj.r.squared ; summary(lm(ACM_deltaET ~ water$ppt_avg))$adj.r.squared 
#summary(lm(SPA_deltaET ~ water$vpd_avg))$adj.r.squared ; summary(lm(ACM_deltaET ~ water$vpd_avg))$adj.r.squared

spa_transpiration = (validation_water_output$drivers$Evap - validation_water_output$drivers$soilevap - validation_water_output$drivers$wetevap)
acm_et = validation_water_output$mean_transpiration + validation_water_output$mean_soilevaporation + validation_water_output$mean_wetcanopyevap
#tmp = which(validation_water_output$drivers$GPP > 1 & spa_transpiration > 0.5 & validation_water_output$drivers$lai > 1.0)
tmp = which(validation_water_output$drivers$GPP > 0.1 & spa_transpiration > 0.1 & validation_water_output$mean_gpp > 0.1 & validation_water_output$mean_transpiration > 0.1 & validation_water_output$drivers$lai > 0.1)
acm_wue = validation_water_output$mean_gpp / validation_water_output$mean_transpiration
spa_wue = validation_water_output$drivers$GPP / spa_transpiration
print(paste("       GPP Mean   = ",round(mean(as.vector(validation_water_output$mean_gpp)),digits=2),sep=""))
print(paste("       Trans Mean = ",round(mean(as.vector(validation_water_output$mean_transpiration)),digits=2),sep=""))
print(paste("       Esoil Mean = ",round(mean(as.vector(validation_water_output$mean_soilevaporation)),digits=2),sep=""))
print(paste("       Ewet Mean  = ",round(mean(as.vector(validation_water_output$mean_wetcanopyevap)),digits=2),sep=""))
print(paste("-------------------","-",sep=""))
print(paste("       GPP R2     = ",round(validation_water_output$gpp_r2,digits=2),sep=""))
print(paste("       GPP Bias   = ",round(validation_water_output$gpp_bias,digits=2),sep=""))
print(paste("       GPP RMSE   = ",round(validation_water_output$gpp_rmse,digits=2),sep=""))
print(paste("       Trans R2   = ",round(validation_water_output$transpiration_r2,digits=2),sep=""))
print(paste("       Trans Bias = ",round(validation_water_output$transpiration_bias,digits=2),sep=""))
print(paste("       Trans RMSE = ",round(validation_water_output$transpiration_rmse,digits=2),sep=""))
print(paste("       Esoil R2   = ",round(validation_water_output$soilevaporation_r2,digits=2),sep=""))
print(paste("       Esoil Bias = ",round(validation_water_output$soilevaporation_bias,digits=2),sep=""))
print(paste("       Esoil RMSE = ",round(validation_water_output$soilevaporation_rmse,digits=2),sep=""))
print(paste("       Ewet R2    = ",round(validation_water_output$wetcanopyevap_r2,digits=2),sep=""))
print(paste("       Ewet Bias  = ",round(validation_water_output$wetcanopyevap_bias,digits=2),sep=""))
print(paste("       Ewet RMSE  = ",round(validation_water_output$wetcanopyevap_rmse,digits=2),sep=""))
print(paste("       SPA WUE    = ",round(mean(spa_wue[tmp]),digits=1)," gC/kgH2O",sep=""))
print(paste("       ACM WUE    = ",round(mean(acm_wue[tmp]),digits=1)," gC/kgH2O",sep=""))
print(paste("       WUE R2     = ",round(summary(lm(spa_wue[tmp] ~ acm_wue[tmp]))$adj.r.squared,digits=2)," ",sep=""))
print(paste("       WUE BIAS   = ",round(mean(acm_wue[tmp]-spa_wue[tmp],na.rm=TRUE),digits=2)," gC/kgH2O",sep=""))
print(paste("       WUE RMSE   = ",round(rmse(acm_wue[tmp],spa_wue[tmp]),digits=2)," gC/kgH2O",sep=""))
print(paste("       SPA T/ET   = ",round(sum(spa_transpiration)/sum(validation_water_output$drivers$Evap),digits=2),sep=""))
print(paste("       SPA S/ET   = ",round(sum(validation_water_output$drivers$soilevap)/sum(validation_water_output$drivers$Evap),digits=2),sep=""))
print(paste("       SPA W/ET   = ",round(sum(validation_water_output$drivers$wetevap)/sum(validation_water_output$drivers$Evap),digits=2),sep=""))
print(paste("       ACM T/ET   = ",round(sum(validation_water_output$mean_transpiration)/sum(acm_et),digits=2),sep=""))
print(paste("       ACM S/ET   = ",round(sum(validation_water_output$mean_soilevaporation)/sum(acm_et),digits=2),sep=""))
print(paste("       ACM W/ET   = ",round(sum(validation_water_output$mean_wetcanopyevap)/sum(acm_et),digits=2),sep=""))

print("Independent validation (fluxnet2015): ")
nos_years = 15
steps_per_year = ceiling(dim(fluxnet_validation_output$observation_wue)[2]/nos_years)
nos_sites = dim(fluxnet_validation_output$observation_wue)[1]
WUE_r2 = array(NA,dim=c(nos_sites))
GPP_r2 = array(NA,dim=c(nos_sites))
ET_r2 = array(NA,dim=c(nos_sites))
GPP_rmse = array(NA,dim=c(nos_sites))
ET_rmse = array(NA,dim=c(nos_sites))
model_WUE = array(NA,dim=c(nos_sites,nos_years))
obs_WUE = array(NA,dim=c(nos_sites,nos_years))
model_annual_gpp = array(NA,dim=c(nos_sites,nos_years))
model_annual_ET = array(NA,dim=c(nos_sites,nos_years))
fluxnet_annual_gpp = array(NA,dim=c(nos_sites,nos_years))
fluxnet_annual_ET = array(NA,dim=c(nos_sites,nos_years))
model_ET = fluxnet_validation_output$timeseries_transpiration + fluxnet_validation_output$timeseries_soilevaporation + fluxnet_validation_output$timeseries_wetcanopyevap
model_GPP = fluxnet_validation_output$timeseries_gpp
obs_GPP = fluxnet_validation_output$observation_gpp
obs_ET = fluxnet_validation_output$observation_ET
model_vpd = fluxnet_validation_output$timeseries_vpd
obs_vpd = fluxnet_validation_output$observed_vpd
for (i in seq(1, dim(fluxnet_validation_output$observation_wue)[1])) {
    # these fluxes should be prior to WUE restrictions...
    tmp = which(is.na(obs_GPP[i,]) == FALSE)
    mod = daily_mean(model_GPP[i,tmp],steps_per_year,steps_per_year-1) ; obs = daily_mean(obs_GPP[i,tmp],steps_per_year,steps_per_year-1)
    if (length(which(is.na(obs) == FALSE)) > 3 & length(which(is.na(mod) == FALSE)) > 3) {
        GPP_r2[i] = summary(lm(obs  ~ mod ))$adj.r.squared
        GPP_rmse[i] = rmse(obs,mod)
    }
    tmp = which(is.na(obs_ET[i,]) == FALSE)
    mod = daily_mean(model_ET[i,tmp],steps_per_year,steps_per_year-1) ; obs = daily_mean(obs_ET[i,tmp],steps_per_year,steps_per_year-1)
    if (length(which(is.na(obs) == FALSE)) > 3 & length(which(is.na(mod) == FALSE)) > 3) {
        ET_r2[i] = summary(lm(obs  ~ mod ))$adj.r.squared
        ET_rmse[i] = rmse(obs,mod)
    }
    filter = which(model_GPP[i,] < 0.2 | model_ET[i,] < 0.2 | obs_GPP[i,] < 0.2 | obs_ET[i,] < 0.2 | fluxnet_validation_output$timeseries_lai[i,] < 0.3 | is.na(obs_ET[i,]) | is.na(obs_GPP[i,]) & fluxnet_validation_output$timeseries_precipitation[i,] >= 2.386993e-05)
    model_GPP[i,filter] = NA ; model_ET[i,filter] = NA
    obs_GPP[i,filter] = NA ; obs_ET[i,filter] = NA
    model_annual_gpp[i,] = daily_sum(model_GPP[i,]*model_vpd[i,],steps_per_year,steps_per_year-1) 
    model_annual_ET[i,] = daily_sum(model_ET[i,],steps_per_year,steps_per_year-1) 
    fluxnet_annual_gpp[i,] = daily_sum(obs_GPP[i,]*obs_vpd[i,],steps_per_year,steps_per_year-1) 
    fluxnet_annual_ET[i,] = daily_sum(obs_ET[i,],steps_per_year,steps_per_year-1) 
    if (length(which(is.na(model_annual_gpp[i,]) == FALSE)) > 3) {
        model_WUE[i,] = model_annual_gpp[i,] / model_annual_ET[i,]
        obs_WUE[i,] = fluxnet_annual_gpp[i,] / fluxnet_annual_ET[i,]
        WUE_r2[i] = summary(lm(obs_WUE[i,] ~ model_WUE[i,]))$adj.r.squared
    }
}
print(paste("Annual:Obs IWUE = ",round(mean(obs_WUE,na.rm=TRUE),digits=1)," gC/kgH2O",sep=""))
print(paste("       ACM IWUE = ",round(mean(model_WUE,na.rm=TRUE),digits=1)," gC/kgH2O",sep=""))
print(paste("        IWUE R2 = ",round(mean(WUE_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("         GPP R2 = ",round(mean(GPP_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("          ET R2 = ",round(mean(ET_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("       GPP RMSE = ",round(mean(GPP_rmse,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("        ET RMSE = ",round(mean(ET_rmse,na.rm=TRUE),digits=2)," ",sep=""))

steps_per_year = ceiling(dim(fluxnet_validation_output$observation_wue)[2]/nos_years)
days_in_period=steps_per_year/12
nos_sites = dim(fluxnet_validation_output$observation_wue)[1]
WUE_r2 = array(NA,dim=c(nos_sites))
GPP_r2 = array(NA,dim=c(nos_sites))
ET_r2 = array(NA,dim=c(nos_sites))
GPP_rmse = array(NA,dim=c(nos_sites))
ET_rmse = array(NA,dim=c(nos_sites))
model_monthly_gpp = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/days_in_period))))
model_monthly_ET = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/days_in_period))))
fluxnet_monthly_gpp = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/days_in_period))))
fluxnet_monthly_ET = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/days_in_period))))
model_WUE = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/days_in_period))))
obs_WUE = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/days_in_period))))
model_ET = fluxnet_validation_output$timeseries_transpiration + fluxnet_validation_output$timeseries_soilevaporation + fluxnet_validation_output$timeseries_wetcanopyevap
model_GPP = fluxnet_validation_output$timeseries_gpp
obs_GPP = fluxnet_validation_output$observation_gpp
obs_ET = fluxnet_validation_output$observation_ET
for (i in seq(1, dim(fluxnet_validation_output$observation_wue)[1])) {
    # these fluxes should be prior to WUE restrictions...
    tmp = which(is.na(obs_GPP[i,]) == FALSE)
    mod = daily_mean(model_GPP[i,tmp],days_in_period,days_in_period-4) ; obs = daily_mean(obs_GPP[i,tmp],days_in_period,days_in_period-4)
    if (length(which(is.na(obs) == FALSE)) > 3 & length(which(is.na(mod) == FALSE)) > 3) {
        GPP_r2[i] = summary(lm(obs  ~ mod ))$adj.r.squared
        GPP_rmse[i] = rmse(obs,mod)
    }
    tmp = which(is.na(obs_ET[i,]) == FALSE)
    mod = daily_mean(model_ET[i,tmp],days_in_period,days_in_period-4) ; obs = daily_mean(obs_ET[i,tmp],days_in_period,days_in_period-4)
    if (length(which(is.na(obs) == FALSE)) > 3 & length(which(is.na(mod) == FALSE)) > 3) {
        ET_r2[i] = summary(lm(obs  ~ mod ))$adj.r.squared
        ET_rmse[i] = rmse(obs,mod)
     }
    filter = which(model_GPP[i,] < 0.2 | model_ET[i,] < 0.2 | obs_GPP[i,] < 0.2 | obs_ET[i,] < 0.2 | fluxnet_validation_output$timeseries_lai[i,] < 0.3 | is.na(obs_ET[i,]) | is.na(obs_GPP[i,]) & fluxnet_validation_output$timeseries_precipitation[i,] >= 2.386993e-05)
    model_GPP[i,filter] = NA ; model_ET[i,filter] = NA
    obs_GPP[i,filter] = NA ; obs_ET[i,filter] = NA
    model_monthly_gpp[i,] = daily_sum(model_GPP[i,]*model_vpd[i,],days_in_period,days_in_period-4) 
    model_monthly_ET[i,] = daily_sum(model_ET[i,],days_in_period,days_in_period-4) 
    fluxnet_monthly_gpp[i,] = daily_sum(obs_GPP[i,]*obs_vpd[i,],days_in_period,days_in_period-4) 
    fluxnet_monthly_ET[i,] = daily_sum(obs_ET[i,],days_in_period,days_in_period-4) 
    if (length(which(is.na(model_monthly_gpp[i,]) == FALSE)) > 3) {
        model_WUE[i,] = model_monthly_gpp[i,] / model_monthly_ET[i,]
        obs_WUE[i,] = fluxnet_monthly_gpp[i,] / fluxnet_monthly_ET[i,]
        WUE_r2[i] = summary(lm(obs_WUE[i,] ~ model_WUE[i,]))$adj.r.squared
    }
}
print(paste("Monthly:Obs IWUE = ",round(mean(obs_WUE,na.rm=TRUE),digits=1)," gC/kgH2O",sep=""))
print(paste("        ACM IWUE = ",round(mean(model_WUE,na.rm=TRUE),digits=1)," gC/kgH2O",sep=""))
print(paste("       IWUE R2   = ",round(mean(WUE_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("        GPP R2   = ",round(mean(GPP_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("        ET R2    = ",round(mean(ET_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("        GPP RMSE = ",round(mean(GPP_rmse,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("        ET RMSE  = ",round(mean(ET_rmse,na.rm=TRUE),digits=2)," ",sep=""))

days_in_period=7
steps_per_year = ceiling(dim(fluxnet_validation_output$observation_wue)[2]/nos_years)
nos_sites = dim(fluxnet_validation_output$observation_wue)[1]
WUE_r2 = array(NA,dim=c(nos_sites))
GPP_r2 = array(NA,dim=c(nos_sites))
ET_r2 = array(NA,dim=c(nos_sites))
GPP_rmse = array(NA,dim=c(nos_sites))
ET_rmse = array(NA,dim=c(nos_sites))
model_weekly_gpp = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/7))))
model_weekly_ET = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/7))))
fluxnet_weekly_gpp = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/7))))
fluxnet_weekly_ET = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/7))))
model_WUE = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/7))))
obs_WUE = array(NA,dim=c(nos_sites,ceiling(nos_years*(365.25/7))))
model_ET = fluxnet_validation_output$timeseries_transpiration + fluxnet_validation_output$timeseries_soilevaporation + fluxnet_validation_output$timeseries_wetcanopyevap
model_GPP = fluxnet_validation_output$timeseries_gpp
obs_GPP = fluxnet_validation_output$observation_gpp
obs_ET = fluxnet_validation_output$observation_ET
for (i in seq(1, dim(fluxnet_validation_output$observation_wue)[1])) {
    # these fluxes should be prior to WUE restrictions...
    tmp = which(is.na(obs_GPP[i,]) == FALSE)
    mod = daily_mean(model_GPP[i,tmp],days_in_period,days_in_period-1) ; obs = daily_mean(obs_GPP[i,tmp],days_in_period,days_in_period-1)
    if (length(which(is.na(obs) == FALSE)) > 3 & length(which(is.na(mod) == FALSE)) > 3) {
	      GPP_r2[i] = summary(lm(obs  ~ mod ))$adj.r.squared
	      GPP_rmse[i] = rmse(obs,mod)
    }
    tmp = which(is.na(obs_ET[i,]) == FALSE)
    mod = daily_mean(model_ET[i,tmp],days_in_period,days_in_period-1) ; obs = daily_mean(obs_ET[i,tmp],days_in_period,days_in_period-1)
    if (length(which(is.na(obs) == FALSE)) > 3 & length(which(is.na(mod) == FALSE)) > 3) {
	      ET_r2[i] = summary(lm(obs  ~ mod ))$adj.r.squared
    	  ET_rmse[i] = rmse(obs,mod)
    }
    filter = which(model_GPP[i,] < 0.2 | model_ET[i,] < 0.2 | obs_GPP[i,] < 0.2 | obs_ET[i,] < 0.2 | fluxnet_validation_output$timeseries_lai[i,] < 0.3 | is.na(obs_ET[i,]) | is.na(obs_GPP[i,]) & fluxnet_validation_output$timeseries_precipitation[i,] >= 2.386993e-05)
    model_GPP[i,filter] = NA ; model_ET[i,filter] = NA
    obs_GPP[i,filter] = NA ; obs_ET[i,filter] = NA
    model_weekly_gpp[i,] = daily_sum(model_GPP[i,]*model_vpd[i,],days_in_period,days_in_period-1) 
    model_weekly_ET[i,] = daily_sum(model_ET[i,],days_in_period,days_in_period-1) 
    fluxnet_weekly_gpp[i,] = daily_sum(obs_GPP[i,]*obs_vpd[i,],days_in_period,days_in_period-1) 
    fluxnet_weekly_ET[i,] = daily_sum(obs_ET[i,],days_in_period,days_in_period-1) 
    if (length(which(is.na(model_weekly_gpp[i,]) == FALSE)) > 3) {
	      model_WUE[i,] = model_weekly_gpp[i,] / model_weekly_ET[i,]
	      obs_WUE[i,] = fluxnet_weekly_gpp[i,] / fluxnet_weekly_ET[i,]
        WUE_r2[i] = summary(lm(obs_WUE[i,] ~ model_WUE[i,]))$adj.r.squared
    }
}
print(paste("Weekly:Obs IWUE = ",round(mean(obs_WUE,na.rm=TRUE),digits=1)," gC/kgH2O",sep=""))
print(paste("       ACM IWUE = ",round(mean(model_WUE,na.rm=TRUE),digits=1)," gC/kgH2O",sep=""))
print(paste("        IWUE R2 = ",round(mean(WUE_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("         GPP R2 = ",round(mean(GPP_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("          ET R2 = ",round(mean(ET_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("       GPP RMSE = ",round(mean(GPP_rmse,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("        ET RMSE = ",round(mean(ET_rmse,na.rm=TRUE),digits=2)," ",sep=""))

days_in_period=1
steps_per_year = ceiling(dim(fluxnet_validation_output$observation_wue)[2]/nos_years)
nos_sites = dim(fluxnet_validation_output$observation_wue)[1]
WUE_r2 = array(NA,dim=c(nos_sites))
GPP_r2 = array(NA,dim=c(nos_sites))
ET_r2 = array(NA,dim=c(nos_sites))
GPP_rmse = array(NA,dim=c(nos_sites))
ET_rmse = array(NA,dim=c(nos_sites))
model_weekly_gpp = array(NA,dim=c(nos_sites,dim(fluxnet_validation_output$timeseries_soilevaporation)[2]))
model_weekly_ET = array(NA,dim=c(nos_sites,dim(fluxnet_validation_output$timeseries_soilevaporation)[2]))
fluxnet_weekly_gpp = array(NA,dim=c(nos_sites,dim(fluxnet_validation_output$timeseries_soilevaporation)[2]))
fluxnet_weekly_ET = array(NA,dim=c(nos_sites,dim(fluxnet_validation_output$timeseries_soilevaporation)[2]))
model_WUE = array(NA,dim=c(nos_sites,dim(fluxnet_validation_output$timeseries_soilevaporation)[2]))
obs_WUE = array(NA,dim=c(nos_sites,dim(fluxnet_validation_output$timeseries_soilevaporation)[2]))
model_ET = fluxnet_validation_output$timeseries_transpiration + fluxnet_validation_output$timeseries_soilevaporation + fluxnet_validation_output$timeseries_wetcanopyevap
model_GPP = fluxnet_validation_output$timeseries_gpp
obs_GPP = fluxnet_validation_output$observation_gpp
obs_ET = fluxnet_validation_output$observation_ET
for (i in seq(1, dim(fluxnet_validation_output$observation_wue)[1])) {
    # these fluxes should be prior to WUE restrictions...
    tmp = which(is.na(obs_GPP[i,]) == FALSE)
    mod = daily_mean(model_GPP[i,tmp],days_in_period,days_in_period) ; obs = daily_mean(obs_GPP[i,tmp],days_in_period,days_in_period)
    if (length(which(is.na(obs) == FALSE)) > 3 & length(which(is.na(mod) == FALSE)) > 3) {
        GPP_r2[i] = summary(lm(obs  ~ mod ))$adj.r.squared
        GPP_rmse[i] = rmse(obs,mod)
    }
    tmp = which(is.na(obs_ET[i,]) == FALSE)
    mod = daily_mean(model_ET[i,tmp],days_in_period,days_in_period) ; obs = daily_mean(obs_ET[i,tmp],days_in_period,days_in_period)
    if (length(which(is.na(obs) == FALSE)) > 3 & length(which(is.na(mod) == FALSE)) > 3) {
        ET_r2[i] = summary(lm(obs  ~ mod ))$adj.r.squared
        ET_rmse[i] = rmse(obs,mod)
    }
    filter = which(model_GPP[i,] < 0.2 | model_ET[i,] < 0.2 | obs_GPP[i,] < 0.2 | obs_ET[i,] < 0.2 | fluxnet_validation_output$timeseries_lai[i,] < 0.3 | is.na(obs_ET[i,]) | is.na(obs_GPP[i,]) & fluxnet_validation_output$timeseries_precipitation[i,] >= 2.386993e-05)
    model_GPP[i,filter] = NA ; model_ET[i,filter] = NA
    obs_GPP[i,filter] = NA ; obs_ET[i,filter] = NA
    model_weekly_gpp[i,] = daily_sum(model_GPP[i,]*model_vpd[i,],days_in_period,days_in_period) 
    model_weekly_ET[i,] = daily_sum(model_ET[i,],days_in_period,days_in_period) 
    fluxnet_weekly_gpp[i,] = daily_sum(obs_GPP[i,]*obs_vpd[i,],days_in_period,days_in_period) 
    fluxnet_weekly_ET[i,] = daily_sum(obs_ET[i,],days_in_period,days_in_period) 
    if (length(which(is.na(model_weekly_gpp[i,]) == FALSE)) > 3 & length(which(is.na(fluxnet_weekly_gpp[i,]) == FALSE)) > 3) {
        model_WUE[i,] = model_weekly_gpp[i,] / model_weekly_ET[i,]
        obs_WUE[i,] = fluxnet_weekly_gpp[i,] / fluxnet_weekly_ET[i,]
        WUE_r2[i] = summary(lm(obs_WUE[i,] ~ model_WUE[i,]))$adj.r.squared
    }
}
fluxnet_daily_GPP_r2 = GPP_r2
fluxnet_daily_ET_r2 = ET_r2
fluxnet_daily_WUE_r2 = WUE_r2
print(paste("Daily: Obs IWUE = ",round(mean(obs_WUE,na.rm=TRUE),digits=1)," gC/kgH2O",sep=""))
print(paste("       ACM IWUE = ",round(mean(model_WUE,na.rm=TRUE),digits=1)," gC/kgH2O",sep=""))
print(paste("        IWUE R2 = ",round(mean(WUE_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("         GPP R2 = ",round(mean(GPP_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("          ET R2 = ",round(mean(ET_r2,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("       GPP RMSE = ",round(mean(GPP_rmse,na.rm=TRUE),digits=2)," ",sep=""))
print(paste("        ET RMSE = ",round(mean(ET_rmse,na.rm=TRUE),digits=2)," ",sep=""))

print("Daily GPP by IGBP")
print(aggregate(GPP_r2,by=list(fluxnet_validation_output$IGBP),mean, na.rm=TRUE))
print("Daily ET by IGBP")
print(aggregate(ET_r2,by=list(fluxnet_validation_output$IGBP),mean, na.rm=TRUE))
print("Daily WUE by IGBP")
print(aggregate(WUE_r2,by=list(fluxnet_validation_output$IGBP),mean, na.rm=TRUE))

# print("Global output: ")
# print(paste("       GPP = ",round(global_output$global_mean_annual_gpp,digits=1)," PgC",sep=""))
# print(paste("       Transpiration = ",round(global_output$global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
# print(paste("       Soil Evaporation = ",round(global_output$global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
# print(paste("       Wet Canopy Evaporation = ",round(global_output$global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
# print(paste("       ET = ",round(global_output$global_mean_annual_et,digits=1)," PgH2O",sep=""))
# print(paste("       WUE = ",round(global_output$global_mean_annual_wue,digits=2)," gC/kgH2O",sep=""))
# print(paste("       wSWP = ",round(global_output$global_mean_annual_wSWP,digits=1)," MPa",sep=""))
# print(paste("       Water in root zone = ",round(global_output$global_mean_annual_rootwatermm,digits=1),sep=""))
# 
# print("Global output: co2_plus100")
# print(paste("       GPP = ",round(global_output_co2_plus100$global_mean_annual_gpp,digits=1)," PgC",sep=""))
# print(paste("       Transpiration = ",round(global_output_co2_plus100$global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
# print(paste("       Soil Evaporation = ",round(global_output_co2_plus100$global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
# print(paste("       Wet Canopy Evaporation = ",round(global_output_co2_plus100$global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
# print(paste("       ET = ",round(global_output_co2_plus100$global_mean_annual_et,digits=1)," PgH2O",sep=""))
# print(paste("       WUE = ",round(global_output_co2_plus100$global_mean_annual_wue,digits=2)," gC/kgH2O",sep=""))
# print(paste("       wSWP = ",round(global_output_co2_plus100$global_mean_annual_wSWP,digits=1)," MPa",sep=""))
# print(paste("       Water in root zone = ",round(global_output_co2_plus100$global_mean_annual_rootwatermm,digits=1),sep=""))
# 
# print("Global output: NUE_half")
# print(paste("       GPP = ",round(global_output_NUE_half$global_mean_annual_gpp,digits=1)," PgC",sep=""))
# print(paste("       Transpiration = ",round(global_output_NUE_half$global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
# print(paste("       Soil Evaporation = ",round(global_output_NUE_half$global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
# print(paste("       Wet Canopy Evaporation = ",round(global_output_NUE_half$global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
# print(paste("       ET = ",round(global_output_NUE_half$global_mean_annual_et,digits=1)," PgH2O",sep=""))
# print(paste("       WUE = ",round(global_output_NUE_half$global_mean_annual_wue,digits=2)," gC/kgH2O",sep=""))
# print(paste("       wSWP = ",round(global_output_NUE_half$global_mean_annual_wSWP,digits=1)," MPa",sep=""))
# print(paste("       Water in root zone = ",round(global_output_NUE_half$global_mean_annual_rootwatermm,digits=1),sep=""))
# 
# print("Global output: NUE_half_co2_plus100")
# print(paste("       GPP = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_gpp,digits=1)," PgC",sep=""))
# print(paste("       Transpiration = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
# print(paste("       Soil Evaporation = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
# print(paste("       Wet Canopy Evaporation = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
# print(paste("       ET = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_et,digits=1)," PgH2O",sep=""))
# print(paste("       WUE = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_wue,digits=2)," gC/kgH2O",sep=""))
# print(paste("       wSWP = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_wSWP,digits=1)," MPa",sep=""))
# print(paste("       Water in root zone = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_rootwatermm,digits=1),sep=""))
# 
# print("Global output: NUE_half_Tair_plus1")
# print(paste("       GPP = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_gpp,digits=1)," PgC",sep=""))
# print(paste("       Transpiration = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
# print(paste("       Soil Evaporation = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
# print(paste("       Wet Canopy Evaporation = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
# print(paste("       ET = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_et,digits=1)," PgH2O",sep=""))
# print(paste("       WUE = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_wue,digits=2)," gC/kgH2O",sep=""))
# print(paste("       wSWP = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_wSWP,digits=1)," MPa",sep=""))
# print(paste("       Water in root zone = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_rootwatermm,digits=1),sep=""))
# 
# ###
# ## GPP mean annual correlations
# 
# nos_years = 15 ; steps_per_years = dim(global_output_NUE_half$timeseries_gpp)[3] / nos_years
# annual_gpp = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2],nos_years))
# annual_et = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2],nos_years))
# annual_wue = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2],nos_years))
# annual_lai = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2],nos_years))
# annual_root = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2],nos_years))
# 
# r2_gpp_temperature = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_gpp_radiation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_gpp_vpd = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_gpp_precipitation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_et_temperature = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_et_radiation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_et_vpd = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_et_precipitation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_wue_temperature = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_wue_radiation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_wue_vpd = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# r2_wue_precipitation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
# 
# ## global_output_NUE_half
# # Mean annual conditions
# for (slot_i in seq(1,dim(global_output_NUE_half$timeseries_gpp)[1])) {
#      for (slot_j in seq(1,dim(global_output_NUE_half$timeseries_gpp)[2])) {
#           # calculate annual averages needed
#           a = 1 ; b = steps_per_years
#           for (y in seq(1,nos_years)) {
#                annual_gpp[slot_i,slot_j,y] = mean(global_output_NUE_half$timeseries_gpp[slot_i,slot_j,a:b],na.rm=TRUE)
#                annual_et[slot_i,slot_j,y] = mean(global_output_NUE_half$timeseries_transpiration[slot_i,slot_j,a:b] + global_output_NUE_half$timeseries_wetcanopyevap[slot_i,slot_j,a:b] + global_output_NUE_half$timeseries_soilevaporation[slot_i,slot_j,a:b],na.rm=TRUE)
# #               mod_wue = global_output_NUE_half$timeseries_WUE[slot_i,slot_j,a:b] ; mod_wue[which(mod_wue == Inf)] = NA
#                mod_wue = annual_gpp[slot_i,slot_j,y] / annual_et[slot_i,slot_j,y]
#                annual_wue[slot_i,slot_j,y] = mean(mod_wue,na.rm=TRUE)
#                annual_lai[slot_i,slot_j,y] = mean(global_output_NUE_half$timeseries_lai[slot_i,slot_j,a:b],na.rm=TRUE)
#                annual_root[slot_i,slot_j,y] = mean(global_output_NUE_half$timeseries_root[slot_i,slot_j,a:b],na.rm=TRUE)
#                a = a + steps_per_years ; b = b + steps_per_years
#           } # time
# 	  if (length(which(is.na(annual_gpp[slot_i,slot_j,]) == FALSE)) > 3) {
#               ## then directly calculate their correlation with temperature, radiation, VPD and precipitation
#               # GPP
#               r2_gpp_temperature[slot_i,slot_j] = summary(lm(annual_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_temperature[slot_i,slot_j,]))$adj.r.squared
#               r2_gpp_radiation[slot_i,slot_j] = summary(lm(annual_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_radiation[slot_i,slot_j,]))$adj.r.squared
#               r2_gpp_vpd[slot_i,slot_j] = summary(lm(annual_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_vpd[slot_i,slot_j,]))$adj.r.squared
#               r2_gpp_precipitation[slot_i,slot_j] = summary(lm(annual_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_precipitation[slot_i,slot_j,]))$adj.r.squared
#               # ET
#               r2_et_temperature[slot_i,slot_j] = summary(lm(annual_et[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_temperature[slot_i,slot_j,]))$adj.r.squared
#               r2_et_radiation[slot_i,slot_j] = summary(lm(annual_et[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_radiation[slot_i,slot_j,]))$adj.r.squared
#               r2_et_vpd[slot_i,slot_j] = summary(lm(annual_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_vpd[slot_i,slot_j,]))$adj.r.squared
#               r2_et_precipitation[slot_i,slot_j] = summary(lm(annual_et[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_precipitation[slot_i,slot_j,]))$adj.r.squared
#               # ET
#               r2_wue_temperature[slot_i,slot_j] = summary(lm(annual_wue[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_temperature[slot_i,slot_j,]))$adj.r.squared
#               r2_wue_radiation[slot_i,slot_j] = summary(lm(annual_wue[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_radiation[slot_i,slot_j,]))$adj.r.squared
#               r2_wue_vpd[slot_i,slot_j] = summary(lm(annual_wue[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_vpd[slot_i,slot_j,]))$adj.r.squared
#               r2_wue_precipitation[slot_i,slot_j] = summary(lm(annual_wue[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_precipitation[slot_i,slot_j,]))$adj.r.squared
#           }
#      } # j
# } # i
# 
# fig_height=4000 ; fig_width=7200
# my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral")))) 
# colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
# jpeg(file="./FIGURES/Cal_val_paper_figure_S3.jpg", height=fig_height, width=fig_width, res=400, quality=100)
# par(mfrow=c(3,4), mar=c(1.4, 4.2, 2.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
# # Summary figures GPP, ET, WUE
# var1=(r2_gpp_temperature) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main=expression(paste("Mean temperature (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# mtext("GPP",side = 2,cex=1.8, padj = 1.2)
# var1=(r2_gpp_radiation) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main=expression(paste("Mean Radiation (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(r2_gpp_precipitation) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main=expression(paste("Mean precipitation (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(r2_gpp_vpd) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main=expression(paste("Mean VPD (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(r2_et_temperature) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# mtext("ET",side = 2,cex=1.8, padj = 1.2)
# var1=(r2_et_radiation) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(r2_et_precipitation) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(r2_et_vpd) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(r2_wue_temperature) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# mtext("WUE",side = 2,cex=1.8, padj = 1.2)
# var1=(r2_wue_radiation) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(r2_wue_precipitation) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(r2_wue_vpd) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# dev.off()
# 
# fluxcom_r2_gpp_temperature = array(NA,dim=c(dim(fluxcom_gpp)[1:2]))
# fluxcom_r2_gpp_radiation = array(NA,dim=c(dim(fluxcom_gpp)[1:2]))
# fluxcom_r2_gpp_vpd = array(NA,dim=c(dim(fluxcom_gpp)[1:2]))
# fluxcom_r2_gpp_precipitation = array(NA,dim=c(dim(fluxcom_gpp)[1:2]))
# 
# ## FLUXCOM (2001-2013)
# # Mean annual conditions
# for (slot_i in seq(1,dim(global_output_NUE_half$timeseries_gpp)[1])) {
#      for (slot_j in seq(1,dim(global_output_NUE_half$timeseries_gpp)[2])) {
#           if (length(which(is.na(fluxcom_gpp[slot_i,slot_j,]) == FALSE)) > 0 & length(which(is.na(global_output_NUE_half$mean_annual_temperature[slot_i,slot_j,]) == FALSE)) > 0) {
#               ## then directly calculate their correlation with temperature, radiation, VPD and precipitation
#               fluxcom_r2_gpp_temperature[slot_i,slot_j] = summary(lm(fluxcom_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_temperature[slot_i,slot_j,1:13]))$adj.r.squared
#               fluxcom_r2_gpp_radiation[slot_i,slot_j] = summary(lm(fluxcom_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_radiation[slot_i,slot_j,1:13]))$adj.r.squared
#               fluxcom_r2_gpp_vpd[slot_i,slot_j] = summary(lm(fluxcom_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_vpd[slot_i,slot_j,1:13]))$adj.r.squared
#               fluxcom_r2_gpp_precipitation[slot_i,slot_j] = summary(lm(fluxcom_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_precipitation[slot_i,slot_j,1:13]))$adj.r.squared
#           }
#      } # j
# } # i
# 
# fig_height=4000 ; fig_width=7200
# my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral")))) 
# colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
# jpeg(file="./FIGURES/Cal_val_paper_figure_S4.jpg", height=fig_height, width=fig_width, res=400, quality=100)
# par(mfrow=c(2,2), mar=c(1.4, 4.4, 2.8, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
# # Summary figures GPP, ET, WUE
# var1=(fluxcom_r2_gpp_temperature) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main=expression(paste("Mean temperature (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# mtext("FLUXCOM GPP",side = 2,cex=1.8, padj = 1.4)
# var1=(fluxcom_r2_gpp_radiation) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main=expression(paste("Mean Radiation (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(fluxcom_r2_gpp_precipitation) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main=expression(paste("Mean precipitation (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(fluxcom_r2_gpp_vpd) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main=expression(paste("Mean VPD (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# dev.off()
# 
# jpeg(file="./FIGURES/Cal_val_paper_figure_S5.jpg", height=fig_height, width=fig_width*1.5, res=400, quality=100)
# var1=(fluxcom_gpp_total) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main=expression(paste("FLUXCOM GPP"," (gC ",m^-2," y",r^-1,")")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(global_output_NUE_half$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
# zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
# image.plot(var1, main=expression(paste("ACM-GPP-ET GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(global_output_NUE_half$mean_gpp*365.25)-(fluxcom_gpp_total) ; colour_choices=colour_choices_upper(length(var1))
# image.plot(var1, main=expression(paste("Bias (ACM-GPP-ET - FLUXCOM)"," (gC ",m^-2," y",r^-1,")")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# dev.off()

###
## Begin output of figures
###

###
## Remote single site at which drought response is dependent not on model concept but rather the number of root layers

remove = which(validation_water_output$drivers$lat < 29.12275+1e-4 & validation_water_output$drivers$lat > 29.12275-1e-4 & validation_water_output$drivers$long < 81.0+1e-4 & validation_water_output$drivers$long  > 81.0-1e-4)
validation_water_output$mean_lai[remove] = NA
validation_water_output$mean_wSWP[remove] = NA
validation_water_output$mean_lwp[remove] = NA
validation_water_output$mean_soilevaporation[remove] = NA
validation_water_output$mean_drainagemm[remove] = NA
validation_water_output$mean_wetcanopyevap[remove] = NA
validation_water_output$mean_runoffmm[remove] = NA
validation_water_output$mean_transpiration[remove] = NA
validation_water_output$mean_rootwatermm[remove] = NA
validation_water_output$mean_gpp[remove] = NA
validation_water_output$mean_WUE[remove] = NA
validation_water_output$mean_ci[remove] = NA

validation_nowater_output$mean_lai[remove] = NA
validation_nowater_output$mean_wSWP[remove] = NA
validation_nowater_output$mean_lwp[remove] = NA
validation_nowater_output$mean_soilevaporation[remove] = NA
validation_nowater_output$mean_drainagemm[remove] = NA
validation_nowater_output$mean_wetcanopyevap[remove] = NA
validation_nowater_output$mean_runoffmm[remove] = NA
validation_nowater_output$mean_transpiration[remove] = NA
validation_nowater_output$mean_rootwatermm[remove] = NA
validation_nowater_output$mean_gpp[remove] = NA
validation_nowater_output$mean_WUE[remove] = NA
validation_nowater_output$mean_ci[remove] = NA

###
## Paper figures

# extract the lat / long information for the calibration sites
data = read.csv("/home/lsmallma/gcel/ACM_GPP_ET_RECALIBRATION/data_for_ACM_era_200pixels_continuous.csv",header=TRUE)
#latitude = data$Lat ; longitude = data$Lon
latitude = calibration_output$drivers$lat ; longitude = calibration_output$drivers$long
latlong = paste(latitude,longitude)
latlong = unique(latlong)
latlong = unlist(strsplit(latlong," "))
latitude = as.numeric(latlong[seq(1,length(latlong),2)])
longitude = as.numeric(latlong[seq(2,length(latlong),2)])
# extract lat / long information for the validation sites
latitude_fluxnet = fluxnet_validation_output$latitude
longitude_fluxnet = fluxnet_validation_output$longitude
# filter out those which were rejected
latitude_fluxnet = latitude_fluxnet[which(is.na(fluxnet_validation_output$gpp_r2) == FALSE)]
longitude_fluxnet = longitude_fluxnet[which(is.na(fluxnet_validation_output$gpp_r2) == FALSE)]
latlong_fluxnet = paste(latitude_fluxnet,longitude_fluxnet)
latlong_fluxnet = unique(latlong_fluxnet)
latlong_fluxnet = unlist(strsplit(latlong_fluxnet," "))
latitude_fluxnet = as.numeric(latlong_fluxnet[seq(1,length(latlong_fluxnet),2)])
longitude_fluxnet = as.numeric(latlong_fluxnet[seq(2,length(latlong_fluxnet),2)])
# calculate temperature and precipitation means for each site
site_mean_temperature = rep(NA, length(latitude)) ; site_mean_rainfall = rep(NA, length(latitude))
fluxnet_mean_temperature = rep(NA, length(latitude_fluxnet)) ; fluxnet_mean_rainfall = rep(NA, length(latitude_fluxnet))
# Do calibration sites first...
for (n in seq(1,length(latitude))) {
    tmp = which(data$Lat == latitude[n] & data$Lon == longitude[n])
    site_mean_temperature[n] = mean(as.vector(unlist(data[tmp,7:14])))-273.15 # K->C
    site_mean_rainfall[n] = mean(as.vector(unlist(data[tmp,31:38]))) # kg/m2/s
}
# ...then fluxnet validation sites second
for (n in seq(1,length(latitude_fluxnet))) {
    tmp = which(fluxnet_validation_output$latitude == latitude_fluxnet[n] & fluxnet_validation_output$longitude == longitude_fluxnet[n])
    fluxnet_mean_temperature[n] = fluxnet_validation_output$mean_temperature[tmp] # oC
    fluxnet_mean_rainfall[n] = fluxnet_validation_output$mean_precipitation[tmp] # kg/m2/s
}
# Adjust the lat/long from training dataset to that of the current analysis
cardamom_latitude = cardamom_latitude[which(is.na(cardamom_latitude) == FALSE)]
cardamom_longitude = cardamom_longitude[which(is.na(cardamom_longitude) == FALSE)]
for (n in seq(1,length(latitude))) {
    tmp = closest2d(n,cardamom_latitude,cardamom_longitude,latitude,longitude,1)
    latitude[n] = cardamom_latitude[tmp]
    longitude[n] = cardamom_longitude[tmp] 
}
# adjust to 0-360 and 0-180
longitude = longitude + 180
latitude = latitude + 90
# now do the same for the validation set
for (n in seq(1,length(latitude_fluxnet))) {
    tmp = closest2d(n,cardamom_latitude,cardamom_longitude,latitude_fluxnet,longitude_fluxnet,1)
    latitude_fluxnet[n] = cardamom_latitude[tmp]
    longitude_fluxnet[n] = cardamom_longitude[tmp] 
}
# adjust to 0-360 and 0-180
longitude_fluxnet = longitude_fluxnet + 180
latitude_fluxnet = latitude_fluxnet + 90

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_1.jpg", height=fig_height, width=fig_width, res=400, quality=100)
# Mean status of biophysical inputs
my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral")))) 
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(global_output$mean_lai) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
par(fig=c(0.0,0.70,0.0,1.0),mar=c(0.6, 0.8, 2.4, 0.5), omi=c(0.2, 0.2, 0.2, 0.4))
image.plot(var1, main=expression(paste("Mean LAI"," (",m^2,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
par(fig=c(0.0,0.70,0.0,1.0),new=TRUE)
plot(longitude,latitude, xlab="", ylab="", pch=16,cex=1.25,xaxt = "n", yaxt = "n",ylim=c(1,179), xlim=c(1,359), xaxs="i", yaxs="i")
par(new=TRUE)
plot(longitude_fluxnet,latitude_fluxnet, xlab="", ylab="", pch=16, col="white",cex=1.25,xaxt = "n", yaxt = "n",ylim=c(1,179), xlim=c(1,359), xaxs="i", yaxs="i")
par(new=TRUE)
plot(longitude_fluxnet,latitude_fluxnet, xlab="", ylab="", pch=1, col="black",cex=1.25,xaxt = "n", yaxt = "n",ylim=c(1,179), xlim=c(1,359), xaxs="i", yaxs="i")
par(fig=c(0.70,1.0,0.0,1.0),mar=c(4.6, 4.6, 0.4, 0.5), omi=c(0.2, 0.2, 0.2, 0.4), new=TRUE)
var1 = as.vector(global_output$mean_temperature)
var2 = as.vector(global_output$mean_precipitation*(365.25*60*60*24))
smoothScatter(var1,var2, ylim=c(0,5000), 
              pch=16, cex.axis=1.6, cex.lab=1.6, cex.main=1.6,nrpoints=0,colramp=my_colours,  
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.015,diff(range(var2,na.rm=TRUE))*0.0025),nbin=128*10,
              ylab="Mean precipitation (kgH2O/m/yr)", xlab="Mean air temperature (Celcuis)")
points(site_mean_temperature, site_mean_rainfall*(365.25*60*60*24), pch=16, col="black",cex=1.25)
points(fluxnet_mean_temperature, fluxnet_mean_rainfall*(365.25*60*60*24), pch=16, col="white",cex=1.25)
points(fluxnet_mean_temperature, fluxnet_mean_rainfall*(365.25*60*60*24), pch=1, col="black",cex=1.25)
dev.off()

my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral"))))
fig_height=4000 ; fig_width=4200
jpeg(file="./FIGURES/Cal_val_paper_figure_2_heat_map_calibration.jpg", height=fig_height, width=fig_width, res=400, quality=100)

par(mfrow=c(2,2), mar=c(3.9, 4.15, 3.8, 0.4), omi=c(0.5, 0.1, 0.1, 0.1))
var1 = calibration_output$mean_gpp
var2 = calibration_output$drivers$GPP
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,
              main="", ylab="SPA", xlab="",cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
text(0,max(var2*0.82),labels=bquote(RMSE == .(round(calibration_output$gpp_rmse,2))), cex=1.2, pos=4)
text(0,max(var2*0.72),labels=bquote(Bias == .(round(calibration_output$gpp_bias,2))), cex=1.2, pos=4)
text(0,max(var2*0.92),labels=bquote(R^2 == .(round(calibration_output$gpp_r2,2))), cex=1.2, pos=4)
mtext("Calibration",side = 3, cex=1.8, padj = -1.7, adj = 1.40)
mtext(expression(paste("GPP (gC/",m^2,"/day)")), side=3, cex = 1.6, padj = -0.1, adj = 0.5)

var1 = calibration_output$mean_transpiration
var2 = (calibration_output$drivers$Evap-calibration_output$drivers$soilevap-calibration_output$drivers$wetevap)
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
text(0,max(var2*0.82),labels=bquote(RMSE == .(round(calibration_output$transpiration_rmse,2))), cex=1.2, pos=4)
text(0,max(var2*0.72),labels=bquote(Bias == .(round(calibration_output$transpiration_bias,2))), cex=1.2, pos=4)
text(0,max(var2*0.92),labels=bquote(R^2 == .(round(calibration_output$transpiration_r2,2))), cex=1.2, pos=4)
mtext(expression(paste("Transpiration (kg",H[2],"O/",m^2,"/day)")), side=3, cex = 1.6, padj = -0.1, adj = 0.5)

var1 = calibration_output$mean_soilevaporation
var2 = calibration_output$drivers$soilevap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="SPA", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
text(0,max(var2*0.82),labels=bquote(RMSE == .(round(calibration_output$soilevaporation_rmse,2))), cex=1.2, pos=4)
text(0,max(var2*0.72),labels=bquote(Bias == .(round(calibration_output$soilevaporation_bias,2))), cex=1.2, pos=4)
text(0,max(var2*0.92),labels=bquote(R^2 == .(round(calibration_output$soilevaporation_r2,2))), cex=1.2, pos=4)
mtext(expression(paste("Soil evaporation (kg",H[2],"O/",m^2,"/day)")), side=3, cex = 1.6, padj = -0.1, adj = 0.5)

var1 = calibration_output$mean_wetcanopyevap
var2 = calibration_output$drivers$wetevap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
text(-0.75,max(var2*0.82),labels=bquote(RMSE == .(round(calibration_output$wetcanopyevap_rmse,2))), cex=1.2, pos=4)
text(-0.75,max(var2*0.72),labels=bquote(Bias == .(round(calibration_output$wetcanopyevap_bias,2))), cex=1.2, pos=4)
text(-0.75,max(var2*0.92),labels=bquote(R^2 == .(round(calibration_output$wetcanopyevap_r2,2))), cex=1.2, pos=4)
mtext(expression(paste("Wet Canopy evaporation (kg",H[2],"O/",m^2,"/day)")), side=3, cex = 1.6, padj = -0.1, adj = 0.5)

dev.off()

# my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral"))))
# fig_height=6000 ; fig_width=4000
# jpeg(file="./FIGURES/Cal_val_paper_figure_2_heat_map.jpg", height=fig_height, width=fig_width, res=400, quality=100)
# par(mfrow=c(4,3), mar=c(4.4, 4.2, 4.4, 1.0), omi=c(0.2, 0.2, 0.3, 0.40))
# var1 = calibration_output$mean_gpp
# var2 = calibration_output$drivers$GPP
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,
#               main="Calibration", ylab="SPA", xlab="",cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(calibration_output$gpp_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(calibration_output$gpp_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(calibration_output$gpp_r2,2))), cex=1.2, pos=4)
# var1 = validation_nowater_output$mean_gpp
# var2 = validation_nowater_output$drivers$GPP
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,
#               main="Validation-Fixed Soil Water", ylab="", xlab="",cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_nowater_output$gpp_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_nowater_output$gpp_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_nowater_output$gpp_r2,2))), cex=1.2, pos=4)
# mtext(expression(paste("GPP"," (gC ",m^-2," da",y^-1,")")),side = 3,cex=1.8, padj = -1.2, adj = 0.5)
# var1 = validation_water_output$mean_gpp
# var2 = validation_water_output$drivers$GPP
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,
#               main="Validation-Dynamic Soil Water", ylab="", xlab="",cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_water_output$gpp_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_water_output$gpp_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_water_output$gpp_r2,2))), cex=1.2, pos=4)
# 
# var1 = calibration_output$mean_transpiration
# var2 = (calibration_output$drivers$Evap-calibration_output$drivers$soilevap-calibration_output$drivers$wetevap)
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="SPA", xlab="",
#               cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(calibration_output$transpiration_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(calibration_output$transpiration_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(calibration_output$transpiration_r2,2))), cex=1.2, pos=4)
# var1 = validation_nowater_output$mean_transpiration
# var2 = (validation_nowater_output$drivers$Evap-validation_nowater_output$drivers$soilevap-validation_nowater_output$drivers$wetevap)
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="",
#               cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_nowater_output$transpiration_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_nowater_output$transpiration_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_nowater_output$transpiration_r2,2))), cex=1.2, pos=4)
# mtext(expression(paste("Transpiration"," (kgH2O ",m^-2," da",y^-1,")")),side = 3,cex=1.8, padj = -0.6, adj = 0.5)
# var1 = validation_water_output$mean_transpiration
# var2 = (validation_water_output$drivers$Evap-validation_water_output$drivers$soilevap-validation_water_output$drivers$wetevap)
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="",cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_water_output$transpiration_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_water_output$transpiration_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_water_output$transpiration_r2,2))), cex=1.2, pos=4)
# 
# var1 = calibration_output$mean_soilevaporation
# var2 = calibration_output$drivers$soilevap
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="SPA", xlab="",
#               cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(calibration_output$soilevaporation_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(calibration_output$soilevaporation_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(calibration_output$soilevaporation_r2,2))), cex=1.2, pos=4)
# var1 = validation_nowater_output$mean_soilevaporation
# var2 = validation_nowater_output$drivers$soilevap
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="",
#               cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_nowater_output$soilevaporation_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_nowater_output$soilevaporation_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_nowater_output$soilevaporation_r2,2))), cex=1.2, pos=4)
# mtext(expression(paste("Soil evaporation"," (kgH2O ",m^-2," da",y^-1,")")),side = 3,cex=1.8, padj = -0.6, adj = 0.5)
# var1 = validation_water_output$mean_soilevaporation
# var2 = validation_water_output$drivers$soilevap
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="",
#               cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_water_output$soilevaporation_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_water_output$soilevaporation_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_water_output$soilevaporation_r2,2))), cex=1.2, pos=4)
# 
# var1 = calibration_output$mean_wetcanopyevap
# var2 = calibration_output$drivers$wetevap
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="SPA", xlab="ACM-GPP-ET",
#               cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(calibration_output$wetcanopyevap_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(calibration_output$wetcanopyevap_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(calibration_output$wetcanopyevap_r2,2))), cex=1.2, pos=4)
# var1 = validation_nowater_output$mean_wetcanopyevap
# var2 = validation_nowater_output$drivers$wetevap
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="ACM-GPP-ET",
#               cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_nowater_output$wetcanopyevap_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_nowater_output$wetcanopyevap_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_nowater_output$wetcanopyevap_r2,2))), cex=1.2, pos=4)
# mtext(expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," da",y^-1,")")),side = 3,cex=1.8, padj = -0.6, adj = 0.5)
# var1 = validation_water_output$mean_wetcanopyevap
# var2 = validation_water_output$drivers$wetevap
# smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="ACM-GPP-ET",
#               cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
#               transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
# text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_water_output$wetcanopyevap_rmse,2))), cex=1.2, pos=4)
# text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_water_output$wetcanopyevap_bias,2))), cex=1.2, pos=4)
# text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_water_output$wetcanopyevap_r2,2))), cex=1.2, pos=4)
# dev.off()

my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral"))))
fig_height=4000 ; fig_width=6000
jpeg(file="./FIGURES/Cal_val_paper_figure_3_validation.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,3), mar=c(4.4, 4.2, 4.4, 1.0), omi=c(0.2, 0.2, 0.3, 0.40))
var1 = validation_water_output$mean_gpp
var2 = validation_water_output$drivers$GPP
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,
              main=expression(paste("GPP"," (gC ",m^-2," da",y^-1,")")), ylab="SPA", xlab="ACM-GPP-ET",cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_water_output$gpp_rmse,2))), cex=1.2, pos=4)
text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_water_output$gpp_bias,2))), cex=1.2, pos=4)
text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_water_output$gpp_r2,2))), cex=1.2, pos=4)

var1 = validation_water_output$mean_transpiration
var2 = (validation_water_output$drivers$Evap-validation_water_output$drivers$soilevap-validation_water_output$drivers$wetevap)
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main=expression(paste("Transpiration"," (kg",H[2],"O ",m^-2," da",y^-1,")")),xlim=range(var1,na.rm=TRUE),ylim=range(var2,na.rm=TRUE), ylab="SPA", xlab="ACM-GPP-ET",cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_water_output$transpiration_rmse,2))), cex=1.2, pos=4)
text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_water_output$transpiration_bias,2))), cex=1.2, pos=4)
text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_water_output$transpiration_r2,2))), cex=1.2, pos=4)
mtext("Validation-Dynamic Soil Water",side = 3,cex=1.8, padj = -2.2, adj = 0.5)

var1 = validation_water_output$mean_soilevaporation
var2 = validation_water_output$drivers$soilevap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main=expression(paste("Soil evaporation"," (kg",H[2],"O ",m^-2," da",y^-1,")")),ylab="SPA", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_water_output$soilevaporation_rmse,2))), cex=1.2, pos=4)
text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_water_output$soilevaporation_bias,2))), cex=1.2, pos=4)
text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_water_output$soilevaporation_r2,2))), cex=1.2, pos=4)

var1 = validation_water_output$mean_wetcanopyevap
var2 = validation_water_output$drivers$wetevap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main=expression(paste("Wet canopy evaporation"," (kg",H[2],"O ",m^-2," da",y^-1,")")), ylab="SPA", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_water_output$wetcanopyevap_rmse,2))), cex=1.2, pos=4)
text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_water_output$wetcanopyevap_bias,2))), cex=1.2, pos=4)
text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_water_output$wetcanopyevap_r2,2))), cex=1.2, pos=4)

var1 = validation_water_output$mean_rootwatermm
var2 = validation_water_output$drivers$SWC
validation_water_output$rootwatermm_rmse = rmse(var1,var2)
validation_water_output$rootwatermm_bias = mean(var1-var2,na.rm=TRUE)
validation_water_output$rootwatermm_r2 = summary(lm(var2 ~ var1))$adj.r.squared
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main=expression(paste("Soil moisture"," (0-10cm; kg",H[2],"O ",m^-2," day)")), ylab="SPA", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,xlim=c(0,55),ylim=c(0,55),
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_water_output$rootwatermm_rmse,2))), cex=1.2, pos=4)
text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_water_output$rootwatermm_bias,2))), cex=1.2, pos=4)
text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_water_output$rootwatermm_r2,2))), cex=1.2, pos=4)

var1 = validation_water_output$mean_gpp / validation_water_output$mean_transpiration
var2 = validation_water_output$drivers$Evap - validation_water_output$drivers$soilevap - validation_water_output$drivers$wetevap
filter = which(validation_water_output$mean_gpp > 0.5 & validation_water_output$mean_transpiration > 0.5 & validation_water_output$drivers$GPP > 0.5 & var2 > 0.5)
var2 = validation_water_output$drivers$GPP / var2
var1 = var1[filter] ; var2 = var2[filter]
validation_water_output$WUE_rmse = rmse(var1,var2)
validation_water_output$WUE_bias = mean(var1-var2,na.rm=TRUE)
validation_water_output$WUE_r2 = summary(lm(var2 ~ var1))$adj.r.squared
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main=expression(paste("WUE"," (gC / kg",H[2],"O )")), ylab="SPA", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
expression(paste("WUE"," (gC / kg",H[2],"O )"))
text(0,max(var2*0.82),labels=bquote(RMSE == .(round(validation_water_output$WUE_rmse,2))), cex=1.2, pos=4)
text(0,max(var2*0.72),labels=bquote(Bias == .(round(validation_water_output$WUE_bias,2))), cex=1.2, pos=4)
text(0,max(var2*0.92),labels=bquote(R^2 == .(round(validation_water_output$WUE_r2,2))), cex=1.2, pos=4)

dev.off()

my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral"))))
fig_height=4000 ; fig_width=6000
jpeg(file="./FIGURES/Cal_val_paper_figure_4_emergent.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,3), mar=c(3.8, 4.6, 1.0, 1.0), omi=c(0.2, 0.2, 0.3, 0.40))
var2 = validation_water_output$drivers$Evap - validation_water_output$drivers$soilevap - validation_water_output$drivers$wetevap
filter = which(validation_water_output$mean_gpp > 0.5 & validation_water_output$mean_transpiration > 0.5 & validation_water_output$drivers$GPP > 0.5 & var2 > 0.5)
var2 = validation_water_output$mean_gpp / validation_water_output$mean_transpiration
var1 = validation_water_output$mean_lai
var1 = var1[filter] ; var2 = var2[filter]
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab=expression(paste("WUE"," (gC / kg",H[2],"O )")), xlab="",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = validation_water_output$mean_lwp
var1 = var1[filter] 
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = validation_water_output$mean_wSWP
var1 = var1[filter] 
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)

var2 = validation_water_output$mean_ci / validation_water_output$drivers$co2_avg
var1 = validation_water_output$mean_lai
var1 = var1[filter] ; var2 = var2[filter]
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab=expression(paste("ci:ca"," (ppm)")), xlab=expression(paste("LAI"," (",m^2,"/",m^2,")")),
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) 
var1 = validation_water_output$mean_lwp
var1 = var1[filter] 
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab=expression(paste("min LWP"," (MPa)")),
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) 
var1 = validation_water_output$mean_wSWP
var1 = var1[filter] 
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab=expression(paste("wSWP"," (MPa)")),
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.7,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) 

dev.off()

my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral"))))
fig_height=6000 ; fig_width=4500
jpeg(file="./FIGURES/Cal_val_paper_figure_5_heat_map.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,2), mar=c(4.4, 4.2, 3.4, 2), omi=c(0.2, 0.2, 0.2, 0.40))
var1 = fluxnet_validation_output$timeseries_gpp
var2 = fluxnet_validation_output$observation_gpp
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,
              main=expression(paste("GPP"," (gC ",m^-2," da",y^-1,")")), ylab="", xlab="",cex=0.5,pch=16,cex.axis=2.0,cex.lab=1.6,cex.main=2.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0055,diff(range(var2,na.rm=TRUE))*0.0055),nbin=128*20) ; abline(0,1,col="red",lwd=3)
text(0,max(var2*0.82,na.rm=TRUE),labels=bquote(RMSE == .(round(mean(fluxnet_validation_output$gpp_rmse,na.rm=TRUE),2))), cex=1.5, pos=4)
text(0,max(var2*0.72,na.rm=TRUE),labels=bquote(Bias == .(round(mean(fluxnet_validation_output$gpp_bias,na.rm=TRUE),2))), cex=1.5, pos=4)
text(0,max(var2*0.92,na.rm=TRUE),labels=bquote(R^2 == .(round(mean(fluxnet_daily_GPP_r2,na.rm=TRUE),2))), cex=1.5, pos=4)
#text(0,max(var2*0.92,na.rm=TRUE),labels=bquote(R^2 == .(round(mean(fluxnet_validation_output$gpp_r2,na.rm=TRUE),2))), cex=1.5, pos=4)
mtext(expression(paste("Fluxnet2015")),side = 2,cex=2.0, padj = -1.5, adj = 0.5)
mtext(expression(paste("ACM-GPP-ET")),side = 1,cex=2.0, padj = 1.75, adj = 1.7)
var1 = fluxnet_validation_output$timeseries_transpiration+fluxnet_validation_output$timeseries_wetcanopyevap+fluxnet_validation_output$timeseries_soilevaporation
var2 = fluxnet_validation_output$observation_ET
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main=expression(paste("Evapo-transpiration"," (kg",H[2],"O ",m^-2," da",y^-1,")")), ylab="", xlab="",
              cex=0.5,pch=16,cex.axis=2.0,cex.lab=1.6,cex.main=2.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0055,diff(range(var2,na.rm=TRUE))*0.0055),nbin=128*20) ; abline(0,1,col="red",lwd=3)
text(max(var1*0,na.rm=TRUE),max(var2*0.82,na.rm=TRUE),labels=bquote(RMSE == .(round(mean(fluxnet_validation_output$ET_rmse,na.rm=TRUE),2))), cex=1.5, pos=4)
text(max(var1*0,na.rm=TRUE),max(var2*0.72,na.rm=TRUE),labels=bquote(Bias == .(round(mean(fluxnet_validation_output$ET_bias,na.rm=TRUE),2))), cex=1.5, pos=4)
text(max(var1*0,na.rm=TRUE),max(var2*0.92,na.rm=TRUE),labels=bquote(R^2 == .(round(mean(fluxnet_daily_ET_r2,na.rm=TRUE),2))), cex=1.5, pos=4)
#text(max(var1*0.92,na.rm=TRUE),max(var2*0.92,na.rm=TRUE),labels=bquote(R^2 == .(round(mean(fluxnet_validation_output$ET_r2,na.rm=TRUE),2))), cex=1.5, pos=4)

hist(fluxnet_validation_output$gpp_r2, main="", xlab=expression(paste("Site specific GPP ",R^2,"",sep="")), ylab="", cex.lab = 2.3,cex.axis = 2, cex=2, col="grey") 
hist(fluxnet_validation_output$ET_r2, main="", xlab=expression(paste("Site specific ET ",R^2,"",sep="")), ylab="", cex.lab = 2.3,cex.axis = 2, cex=2, col="grey") 
hist(fluxnet_validation_output$gpp_rmse, main="", xlab=expression(paste("Site specific GPP RMSE (gC/",m^2,"/day)",sep="")), ylab="", cex.lab = 2.3,cex.axis = 2, cex=2, col="grey")
hist(fluxnet_validation_output$ET_rmse, main="", xlab=expression(paste("Site specific ET RMSE (kg",H[2],"O/",m^2,"/day)",sep="")), ylab="", cex.lab = 2.3,cex.axis = 2, cex=2, col="grey")
dev.off()

fig_height=4000 ; fig_width=7000
# I reorder the groups order : I change the order of the factor data$names
fluxnet_validation_output$IGBP=factor(fluxnet_validation_output$IGBP , levels=levels(fluxnet_validation_output$IGBP)[c(3,8,6,5,7,1,13,2,9,10)])
jpeg(file="./FIGURES/Boxplot_fluxnet_statistics_by_IGBP.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,1), mar=c(5,6,3,1))
boxplot(fluxnet_daily_GPP_r2~fluxnet_validation_output$IGBP, main="GPP", xlab="", ylab=expression(paste(R^2,"",sep="")), cex.lab = 2.3,cex.axis = 2.3, cex=2.3, cex.main=2.3, lwd=2, pch=16, range = 0)
abline(mean(GPP_r2,na.rm=TRUE),0,col="grey", lwd=3)
boxplot(fluxnet_daily_ET_r2~fluxnet_validation_output$IGBP, main="ET", xlab="Vegetation Type", ylab=expression(paste(R^2,"",sep="")), cex.lab = 2.3,cex.axis = 2.3, cex=2.3, cex.main=2.3, lwd=2, pch=16, range = 0) 
abline(mean(ET_r2,na.rm=TRUE),0,col="grey", lwd=3)
dev.off()

###
## Calibration Figures

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/calibration_XY_flux_components.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,3), mar=c(4.4, 4.2, 2.4, 2.5), omi=c(0.2, 0.2, 0.2, 0.40))
plot(calibration_output$drivers$GPP~calibration_output$mean_gpp,main="GPP", ylab="SPA", xlab="ACM_GPP_ET",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,1,col="red",lwd=3)
plot(calibration_output$drivers$Evap~as.vector(calibration_output$mean_transpiration+calibration_output$mean_soilevaporation+calibration_output$mean_wetcanopyevap),main="Evapotranspiration", ylab="SPA", xlab="ACM_GPP_ET",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,1,col="red",lwd=3)
plot(calibration_output$drivers$soilevap~calibration_output$mean_soilevaporation,main="Soil evaporation", ylab="SPA", xlab="ACM_GPP_ET",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,1,col="red",lwd=3)
plot(calibration_output$drivers$wetevap~calibration_output$mean_wetcanopyevap,main="Wet canopy evaporation", ylab="SPA", xlab="ACM_GPP_ET",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,1,col="red",lwd=3)
transpiration = (calibration_output$drivers$Evap-calibration_output$drivers$soilevap-calibration_output$drivers$wetevap)
plot(transpiration~calibration_output$mean_transpiration,main="Transpiration", ylab="SPA", xlab="ACM_GPP_ET",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,1,col="red",lwd=3)
plot(calibration_output$mean_wSWP,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6, main="Weighted soil water potential", ylab="SPA", xlab="ACM_GPP_ET")
dev.off()

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/calibration_XY_transpiration_residuals_vs_drivers.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,3), mar=c(4.4, 4.2, 4.4, 2.5), omi=c(0.2, 0.2, 0.2, 0.40))
transpiration = (calibration_output$drivers$Evap-calibration_output$drivers$soilevap-calibration_output$drivers$wetevap)
plot(transpiration~calibration_output$mean_transpiration,main="", ylab="SPA", xlab="ACM_GPP_ET",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,1,col="red",lwd=3)
transpiration_residual = transpiration - calibration_output$mean_transpiration
plot(transpiration_residual~calibration_output$drivers$sat_max,main="Transpiration", ylab="Residual (ACM-SPA)", xlab="Temperature",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=2.2) ; abline(0,0,col="red",lwd=3)
plot(transpiration_residual~calibration_output$drivers$swrad_avg,main="", ylab="Residual (ACM-SPA)", xlab="Short-Wave radiation",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,0,col="red",lwd=3)
plot(transpiration_residual~calibration_output$drivers$vpd_avg,main="", ylab="Residual (ACM-SPA)", xlab="VPD",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,0,col="red",lwd=3)
plot(transpiration_residual~calibration_output$drivers$wind_avg,main="", ylab="Residual (ACM-SPA)", xlab="Wind",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,0,col="red",lwd=3)
plot(transpiration_residual~calibration_output$drivers$ppt_avg,main="", ylab="Residual (ACM-SPA)", xlab="Rain",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,0,col="red",lwd=3)
dev.off()

###
## global_propogation

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_S1.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,2), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Mean status of biophysical inputs
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(global_output$mean_lai) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Mean LAI"," (",m^2,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
par(new=TRUE) ; plot(longitude,latitude, pch=16,cex=1.0,xaxt = "n", yaxt = "n",ylim=c(1,179), xlim=c(1,359))
var1=(global_output$mean_root) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Mean Fine Root"," (gC/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
par(new=TRUE) ; plot(longitude,latitude, pch=16,cex=1.0,xaxt = "n", yaxt = "n",ylim=c(1,179), xlim=c(1,359))
# Standard deviation of biophysical inputs
var1=(global_output$sd_lai) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Standard Deviation LAI"," (",m^2,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
par(new=TRUE) ; plot(longitude,latitude, pch=16,cex=1.0,xaxt = "n", yaxt = "n",ylim=c(1,179), xlim=c(1,359))
var1=(global_output$sd_root) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Standard Deviation Fine Root"," (gC/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
par(new=TRUE) ; plot(longitude,latitude, pch=16,cex=1.0,xaxt = "n", yaxt = "n",ylim=c(1,179), xlim=c(1,359))
dev.off()

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/All_fluxers_global_output.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(global_output$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output$mean_transpiration+global_output$mean_wetcanopyevap+global_output$mean_soilevaporation)*365.25 ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("ET"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=global_output$mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("WUE (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# ET component figures transpiration, wetcanopy evaporation, soil evaporation
var1=(global_output$mean_transpiration*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Transpiration"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output$mean_wetcanopyevap*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output$mean_soilevaporation*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# Soil runoff / drainage and status
var1=(global_output$mean_runoffmm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil run-off"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output$mean_drainagemm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil drainage"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output$mean_rootwatermm) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Water in rooting zone"," (",kg,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

fig_height=4500 ; fig_width=3000
jpeg(file="./FIGURES/Cal_val_paper_figure_3.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,1), mar=c(1.2, 1.2, 1.8, 1.5), omi=c(0.1, 0.1, 0.1, 0.20))
# Summary figures GPP, ET, WUE
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
# ET component figures transpiration, wetcanopy evaporation, soil evaporation
mean_ET = global_output$mean_transpiration+global_output$mean_wetcanopyevap+global_output$mean_soilevaporation
var1=(global_output$mean_transpiration/mean_ET) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(0,1)
image.plot(var1, main="Transpiration / ET",zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.2,legend.width=2.0,cex=1.8,axis.args=list(cex.axis=2.2,hadj=0.1))
var1=(global_output$mean_wetcanopyevap/mean_ET) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(0,1)
image.plot(var1, main="Wet canopy evaporation / ET",zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.2,legend.width=2.0,cex=1.8,axis.args=list(cex.axis=2.2,hadj=0.1))
var1=(global_output$mean_soilevaporation/mean_ET) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(0,1)
image.plot(var1, main="Soil evaporation / ET",zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.2,legend.width=2.0,cex=1.8,axis.args=list(cex.axis=2.2,hadj=0.1))
dev.off()

###############################################
###############################################
###############################################
###############################################
# 
# soil = read.table("/home/lsmallma/Desktop/acm_gpp_et_soil_test.txt", sep="", header=FALSE)
# par(mfrow=c(2,3))
# plot(soil[,1])
# plot(soil[,2])
# plot(soil[,3])
# 
# plot(soil_residual[test_site]~soil[,1])
# plot(soil_residual[test_site]~soil[,2])
# plot(soil_residual[test_site]~soil[,3])
# 
# setwd("/home/lsmallma/WORK/GREENHOUSE/models/ACM_GPP_ET/")
# load("./outputs/validation_water_output.RData")
# load("./outputs/validation_nowater_output.RData")
# 
# names(validation_water_output)
# names(validation_water_output$drivers)
# 
# gpp_residual = validation_water_output$mean_gpp-validation_water_output$drivers$GPP
# trans_residual = validation_water_output$mean_transpiration-(validation_water_output$drivers$Evap-validation_water_output$drivers$soilevap-validation_water_output$drivers$wetevap)
# soil_residual = validation_water_output$mean_soilevaporation-validation_water_output$drivers$soilevap
# swp_residual = validation_water_output$mean_wSWP-validation_water_output$drivers$SWP
# 
# spa_transpiration = (validation_water_output$drivers$Evap-validation_water_output$drivers$soilevap-validation_water_output$drivers$wetevap)
# 
# locations = which(abs(soil_residual) > 4)
# validation_water_output$drivers$lat[locations]
# validation_water_output$drivers$long[locations]
# test_site = which(validation_water_output$drivers$lat == validation_water_output$drivers$lat[locations[2]] & validation_water_output$drivers$long == validation_water_output$drivers$long[locations[2]])
# par(mfrow=c(3,3))
# plot(validation_water_output$mean_gpp[test_site])
# plot(validation_water_output$mean_lai[test_site])
# plot(validation_water_output$mean_WUE[test_site])
# plot(validation_water_output$mean_transpiration[test_site])
# plot(validation_water_output$mean_soilevaporation[test_site])
# plot(validation_water_output$mean_drainagemm[test_site])
# plot(validation_water_output$mean_runoffmm[test_site])
# plot(validation_water_output$mean_wSWP[test_site])
# plot(validation_water_output$mean_rootwatermm[test_site])
# 
# par(mfrow=c(2,5))
# plot(validation_water_output$mean_gpp[test_site], type="l", col="red", ylim=c(0,max(validation_water_output$drivers$GPP[test_site])))
# lines(validation_water_output$drivers$GPP[test_site], col="black")
# plot(validation_water_output$mean_soilevaporation[test_site], type="l", col="red")
# lines(validation_water_output$drivers$soilevap[test_site], col="black")
# plot(validation_water_output$mean_wSWP[test_site], type="l", col="red")
# lines(validation_water_output$drivers$SWP[test_site], col="black")
# plot(validation_water_output$mean_rootwatermm[test_site], type="l", col="red")
# lines(validation_water_output$drivers$SWC[test_site], col="black")
# plot(validation_water_output$mean_transpiration[test_site], type="l", col="red", ylim=c(0,max(spa_transpiration[test_site])))
# lines(spa_transpiration[test_site], col="black")
# 
# plot(validation_water_output$mean_gpp[test_site]~validation_water_output$drivers$GPP[test_site])
# abline(0,1,col="red", lwd=3)
# plot(validation_water_output$mean_soilevaporation[test_site]~validation_water_output$drivers$soilevap[test_site])
# abline(0,1,col="red", lwd=3)
# plot(validation_water_output$mean_wSWP[test_site]~validation_water_output$drivers$SWP[test_site])
# abline(0,1,col="red", lwd=3)
# plot(validation_water_output$mean_rootwatermm[test_site]~validation_water_output$drivers$SWC[test_site])
# abline(0,1,col="red", lwd=3)
# plot(validation_water_output$mean_transpiration[test_site]~spa_transpiration[test_site])
# abline(0,1,col="red", lwd=3)
# 
# par(mfrow=c(4,4))
# plot(soil_residual[test_site]~validation_water_output$drivers$doy[test_site], main="", ylab="Residual", xlab="DOYT")
# plot(soil_residual[test_site]~validation_water_output$drivers$lai[test_site], main="", ylab="Residual", xlab="LAI")
# plot(soil_residual[test_site]~validation_water_output$drivers$avgN[test_site], main="", ylab="Residual", xlab="avgN")
# plot(soil_residual[test_site]~validation_water_output$drivers$swrad_avg[test_site], main="", ylab="Residual", xlab="swrad_avg")
# plot(soil_residual[test_site]~validation_water_output$drivers$co2_avg[test_site], main="", ylab="Residual", xlab="CO2")
# plot(soil_residual[test_site]~validation_water_output$drivers$sat_max[test_site], main="", ylab="Residual", xlab="sat_max")
# plot(soil_residual[test_site]~validation_water_output$drivers$vpd_avg[test_site], main="", ylab="Residual", xlab="vpd_avg")
# plot(soil_residual[test_site]~validation_water_output$drivers$wind_avg[test_site], main="", ylab="Residual", xlab="wind_avg")
# plot(soil_residual[test_site]~validation_water_output$drivers$ppt_avg[test_site], main="", ylab="Residual", xlab="rainfall")
# plot(soil_residual[test_site]~validation_water_output$drivers$SWP[test_site], main="", ylab="Residual", xlab="wSWP")
# plot(soil_residual[test_site]~validation_water_output$drivers$SWC[test_site], main="", ylab="Residual", xlab="SWC")
# plot(soil_residual[test_site]~validation_water_output$drivers$Rtot[test_site], main="", ylab="Residual", xlab="Rtot")
# plot(soil_residual[test_site]~validation_water_output$drivers$energy_balance_residual[test_site], main="", ylab="Residual", xlab="energy balance")
# plot(soil_residual[test_site]~validation_water_output$drivers$netrad[test_site], main="", ylab="Residual", xlab="NetRad")
# 
# validation_water_output$drivers$porosity[test_site]
# 
# jpeg("./FIGURES/gpp_residual_investigation_water.jpg", width=7200, height = 4000, res = 400, quality = 100)
# par(mfrow=c(3,4))
# plot(gpp_residual~validation_water_output$drivers$lai, main="", ylab="Residual", xlab="LAI")
# plot(gpp_residual~validation_water_output$drivers$avgN, main="", ylab="Residual", xlab="avgN")
# plot(gpp_residual~validation_water_output$drivers$swrad_avg, main="", ylab="Residual", xlab="swrad_avg")
# plot(gpp_residual~validation_water_output$drivers$co2_avg, main="", ylab="Residual", xlab="CO2")
# plot(gpp_residual~validation_water_output$drivers$sat_max, main="", ylab="Residual", xlab="sat_max")
# plot(gpp_residual~validation_water_output$drivers$vpd_avg, main="", ylab="Residual", xlab="vpd_avg")
# plot(gpp_residual~validation_water_output$drivers$wind_avg, main="", ylab="Residual", xlab="wind_avg")
# plot(gpp_residual~validation_water_output$drivers$ppt_avg, main="", ylab="Residual", xlab="rainfall")
# plot(gpp_residual~validation_water_output$drivers$SWP, main="", ylab="Residual", xlab="sfc_pressure")
# plot(gpp_residual~validation_water_output$drivers$SWC, main="", ylab="Residual", xlab="SWC")
# plot(gpp_residual~validation_water_output$drivers$energy_balance_residual, main="", ylab="Residual", xlab="energy balance")
# plot(gpp_residual~validation_water_output$drivers$netrad, main="", ylab="Residual", xlab="NetRad")
# dev.off()
# 
# jpeg("./FIGURES/trans_residual_investigation_water.jpg", width=7200, height = 4000, res = 400, quality = 100)
# par(mfrow=c(3,4))
# plot(trans_residual~validation_water_output$drivers$lai, main="", ylab="Residual", xlab="LAI")
# plot(trans_residual~validation_water_output$drivers$avgN, main="", ylab="Residual", xlab="avgN")
# plot(trans_residual~validation_water_output$drivers$swrad_avg, main="", ylab="Residual", xlab="swrad_avg")
# plot(trans_residual~validation_water_output$drivers$co2_avg, main="", ylab="Residual", xlab="CO2")
# plot(trans_residual~validation_water_output$drivers$sat_max, main="", ylab="Residual", xlab="sat_max")
# plot(trans_residual~validation_water_output$drivers$vpd_avg, main="", ylab="Residual", xlab="vpd_avg")
# plot(trans_residual~validation_water_output$drivers$wind_avg, main="", ylab="Residual", xlab="wind_avg")
# plot(trans_residual~validation_water_output$drivers$ppt_avg, main="", ylab="Residual", xlab="rainfall")
# plot(trans_residual~validation_water_output$drivers$SWP, main="", ylab="Residual", xlab="sfc_pressure")
# plot(trans_residual~validation_water_output$drivers$SWC, main="", ylab="Residual", xlab="SWC")
# plot(trans_residual~validation_water_output$drivers$energy_balance_residual, main="", ylab="Residual", xlab="energy balance")
# plot(trans_residual~validation_water_output$drivers$netrad, main="", ylab="Residual", xlab="NetRad")
# dev.off()
# 
# jpeg("./FIGURES/soil_residual_investigation_water.jpg", width=7200, height = 4000, res = 400, quality = 100)
# par(mfrow=c(3,4))
# plot(soil_residual~validation_water_output$drivers$lai, main="", ylab="Residual", xlab="LAI")
# plot(soil_residual~validation_water_output$drivers$avgN, main="", ylab="Residual", xlab="avgN")
# plot(soil_residual~validation_water_output$drivers$swrad_avg, main="", ylab="Residual", xlab="swrad_avg")
# plot(soil_residual~validation_water_output$drivers$co2_avg, main="", ylab="Residual", xlab="CO2")
# plot(soil_residual~validation_water_output$drivers$sat_max, main="", ylab="Residual", xlab="sat_max")
# plot(soil_residual~validation_water_output$drivers$vpd_avg, main="", ylab="Residual", xlab="vpd_avg")
# plot(soil_residual~validation_water_output$drivers$wind_avg, main="", ylab="Residual", xlab="wind_avg")
# plot(soil_residual~validation_water_output$drivers$ppt_avg, main="", ylab="Residual", xlab="rainfall")
# plot(soil_residual~validation_water_output$drivers$SWP, main="", ylab="Residual", xlab="wSWP")
# plot(soil_residual~validation_water_output$drivers$SWC, main="", ylab="Residual", xlab="SWC")
# plot(soil_residual~validation_water_output$drivers$energy_balance_residual, main="", ylab="Residual", xlab="energy balance")
# plot(soil_residual~validation_water_output$drivers$netrad, main="", ylab="Residual", xlab="NetRad")
# dev.off()
# 
# jpeg("./FIGURES/wSWP_residual_investigation_water.jpg", width=7200, height = 4000, res = 400, quality = 100)
# par(mfrow=c(3,4))
# plot(swp_residual~validation_water_output$drivers$lai, main="", ylab="Residual", xlab="LAI")
# plot(swp_residual~validation_water_output$drivers$avgN, main="", ylab="Residual", xlab="avgN")
# plot(swp_residual~validation_water_output$drivers$swrad_avg, main="", ylab="Residual", xlab="swrad_avg")
# plot(swp_residual~validation_water_output$drivers$co2_avg, main="", ylab="Residual", xlab="CO2")
# plot(swp_residual~validation_water_output$drivers$sat_max, main="", ylab="Residual", xlab="sat_max")
# plot(swp_residual~validation_water_output$drivers$vpd_avg, main="", ylab="Residual", xlab="vpd_avg")
# plot(swp_residual~validation_water_output$drivers$wind_avg, main="", ylab="Residual", xlab="wind_avg")
# plot(swp_residual~validation_water_output$drivers$ppt_avg, main="", ylab="Residual", xlab="rainfall")
# plot(swp_residual~validation_water_output$drivers$sfc_pressure, main="", ylab="Residual", xlab="sfc_pressure")
# plot(swp_residual~validation_water_output$drivers$SWC, main="", ylab="Residual", xlab="SWC")
# plot(swp_residual~validation_water_output$drivers$energy_balance_residual, main="", ylab="Residual", xlab="energy balance")
# plot(swp_residual~validation_water_output$drivers$netrad, main="", ylab="Residual", xlab="NetRad")
# dev.off()
# 
# par(mfrow=c(2,2))
# plot(soil_residual~validation_water_output$drivers$SWP, main="", ylab="Residual", xlab="SWP")
# plot(soil_residual~validation_water_output$mean_wSWP, main="", ylab="Residual", xlab="SWP")
# plot(soil_residual[which(validation_water_output$mean_wSWP > -0.1)]~validation_water_output$drivers$swrad_avg[which(validation_water_output$mean_wSWP > -0.1)], main="", ylab="Residual", xlab="SW-RAD")
# plot(soil_residual[which(validation_water_output$mean_wSWP > -0.1)]~validation_water_output$drivers$swrad_avg[which(validation_water_output$mean_wSWP > -0.1)], main="", ylab="Residual", xlab="SW-RAD")
# 
# par(mfrow=c(2,2))
# plot(validation_water_output$mean_rootwatermm~validation_water_output$drivers$SWC, main="Rootwatermm", ylab="ACM", xlab="SPA")
# abline(0,1,col="red",lwd=3)
# plot(validation_water_output$mean_wSWP~validation_water_output$drivers$SWP, main="wSWP", ylab="ACM", xlab="SPA")
# abline(0,1,col="red",lwd=3)
# plot(validation_water_output$mean_rootwatermm[which(abs(swp_residual) < 0.2)]~validation_water_output$drivers$SWC[which(abs(swp_residual) < 0.2)], main="Rootwatermm", ylab="ACM", xlab="SPA")
# abline(0,1,col="red",lwd=3)
# plot(validation_water_output$mean_wSWP[which(abs(swp_residual) < 0.2)]~validation_water_output$drivers$SWP[which(abs(swp_residual) < 0.2)], main="wSWP", ylab="ACM", xlab="SPA", ylim=c(-0.5,0))
# abline(0,1,col="red",lwd=3)
# 
# length(which(abs(swp_residual) < 0.2)) / length(swp_residual)
# summary(lm(validation_water_output$mean_soilevaporation[which(abs(swp_residual) < 0.2)]~validation_water_output$drivers$soilevap[which(abs(swp_residual) < 0.2)]))
# 
# soilevap_rad_intercept    = 1.122969e-02 ; soilevap_rad_coef = 1.748044e+00 
# soilevap = soilevap_rad_intercept + (soilevap * soilevap_rad_coef)
# 
# uncorrected_soilevap = (validation_water_output$mean_soilevaporation - soilevap_rad_intercept) / soilevap_rad_coef
# plot(validation_water_output$mean_soilevaporation~uncorrected_soilevap)
# plot(validation_water_output$drivers$soilevap~uncorrected_soilevap)
# 
# summary(lm(validation_water_output$drivers$soilevap~uncorrected_soilevap))
# 
# usoil_residual = uncorrected_soilevap-validation_water_output$drivers$soilevap
# jpeg("./FIGURES/usoil_residual_investigation_water.jpg", width=7200, height = 4000, res = 400, quality = 100)
# par(mfrow=c(3,4))
# plot(usoil_residual~validation_water_output$drivers$lai, main="", ylab="Residual", xlab="LAI")
# plot(usoil_residual~validation_water_output$drivers$avgN, main="", ylab="Residual", xlab="avgN")
# plot(usoil_residual~validation_water_output$drivers$swrad_avg, main="", ylab="Residual", xlab="swrad_avg")
# plot(usoil_residual~validation_water_output$drivers$co2_avg, main="", ylab="Residual", xlab="CO2")
# plot(usoil_residual~validation_water_output$drivers$sat_max, main="", ylab="Residual", xlab="sat_max")
# plot(usoil_residual~validation_water_output$drivers$vpd_avg, main="", ylab="Residual", xlab="vpd_avg")
# plot(usoil_residual~validation_water_output$drivers$wind_avg, main="", ylab="Residual", xlab="wind_avg")
# plot(usoil_residual~validation_water_output$drivers$ppt_avg, main="", ylab="Residual", xlab="rainfall")
# plot(usoil_residual~validation_water_output$drivers$sfc_pressure, main="", ylab="Residual", xlab="sfc_pressure")
# plot(usoil_residual~validation_water_output$drivers$SWC, main="", ylab="Residual", xlab="SWC")
# plot(usoil_residual~validation_water_output$drivers$energy_balance_residual, main="", ylab="Residual", xlab="energy balance")
# plot(usoil_residual~validation_water_output$drivers$netrad, main="", ylab="Residual", xlab="NetRad")
# dev.off()
# 
# uncorrected_soilevap = (validation_nowater_output$mean_soilevaporation - soilevap_rad_intercept) / soilevap_rad_coef
# usoil_residual = uncorrected_soilevap-validation_nowater_output$drivers$soilevap
# jpeg("./FIGURES/usoil_residual_investigation_nowater.jpg", width=7200, height = 4000, res = 400, quality = 100)
# par(mfrow=c(3,4))
# plot(usoil_residual~validation_nowater_output$drivers$lai, main="", ylab="Residual", xlab="LAI")
# plot(usoil_residual~validation_nowater_output$drivers$avgN, main="", ylab="Residual", xlab="avgN")
# plot(usoil_residual~validation_nowater_output$drivers$swrad_avg, main="", ylab="Residual", xlab="swrad_avg")
# plot(usoil_residual~validation_nowater_output$drivers$co2_avg, main="", ylab="Residual", xlab="CO2")
# plot(usoil_residual~validation_nowater_output$drivers$sat_max, main="", ylab="Residual", xlab="sat_max")
# plot(usoil_residual~validation_nowater_output$drivers$vpd_avg, main="", ylab="Residual", xlab="vpd_avg")
# plot(usoil_residual~validation_nowater_output$drivers$wind_avg, main="", ylab="Residual", xlab="wind_avg")
# plot(usoil_residual~validation_nowater_output$drivers$ppt_avg, main="", ylab="Residual", xlab="rainfall")
# plot(usoil_residual~validation_nowater_output$drivers$sfc_pressure, main="", ylab="Residual", xlab="sfc_pressure")
# plot(usoil_residual~validation_nowater_output$drivers$SWC, main="", ylab="Residual", xlab="SWC")
# plot(usoil_residual~validation_nowater_output$drivers$energy_balance_residual, main="", ylab="Residual", xlab="energy balance")
# plot(usoil_residual~validation_nowater_output$drivers$netrad, main="", ylab="Residual", xlab="NetRad")
# dev.off()
# 
# 
# 
# 
# 
# 
# setwd("/Users/lsmallma/WORK/Github/GCEL/ACM_GPP_ET/")
# load("./outputs/global_1x1_degree_2001_2015_NUE_half.RData")
# load("./outputs/global_1x1_degree_2001_2015_NUE_half_co2_plus100.RData")
# load("./outputs/global_1x1_degree_2001_2015_NUE_half_Tair_plus4.RData")
# 
# fig_height=4000 ; fig_width=7200
# jpeg(file="./FIGURES/Cal_val_paper_figure_5.jpg", height=fig_height, width=fig_width, res=400, quality=100)
# par(mfrow=c(2,3), mar=c(1.4, 1.2, 2.4, 8.5), omi=c(0.1, 0.1, 0.1, 0.1))
# ## Summary figures GPP, ET, WUE; CO2 + 100
# colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
# var1=((global_output_NUE_half_co2_plus100$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25))/(global_output_NUE_half$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
# zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
# image.plot(var1, main=expression(paste(Delta,"GPP"," (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# mtext(expression(paste("C",O[2],"+100 ppm")),side = 2,cex=1.8, padj = 1.2)
# var1=((global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25)
# var1 = var1-((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
# var1 = var1 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)  
# colour_choices=colour_choices_upper(length(var1))
# zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
# image.plot(var1, main=expression(paste(Delta,"ET"," (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1=(global_output_NUE_half_co2_plus100$mean_wue - global_output_NUE_half$mean_wue) / global_output_NUE_half$mean_wue ; colour_choices=colour_choices_upper(length(var1))
# zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
# image.plot(var1, main=expression(paste(Delta,"WUE (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# ## Summary figures GPP, ET, WUE; Tair + 1
# var1=((global_output_NUE_half_Tair_plus4$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25))/(global_output_NUE_half$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
# zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
# image.plot(var1, main=expression(paste(Delta,"GPP"," (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# mtext(expression(paste("Tair+",1^o,"C",sep="")),side = 2,cex=1.8, padj = 1.2)
# var1=((global_output_NUE_half_Tair_plus4$mean_transpiration+global_output_NUE_half_Tair_plus4$mean_wetcanopyevap+global_output_NUE_half_Tair_plus4$mean_soilevaporation)*365.25)
# var1 = var1-((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
# var1 = var1 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)  
# colour_choices=colour_choices_upper(length(var1))
# zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
# image.plot(var1, main=expression(paste(Delta,"ET"," (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# var1 = (global_output_NUE_half_Tair_plus4$mean_wue - global_output_NUE_half$mean_wue) / global_output_NUE_half$mean_wue ; colour_choices=colour_choices_upper(length(var1))
# zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
# image.plot(var1, main=expression(paste(Delta,"WUE (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# dev.off()
# 
