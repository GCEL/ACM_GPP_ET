
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
load("./outputs/global_1x1_degree_2001_2015_co2_plus100.RData")
load("./outputs/global_1x1_degree_2001_2015_Tair_plus1.RData")
load("./outputs/global_1x1_degree_2001_2015_NUE_half.RData")
load("./outputs/global_1x1_degree_2001_2015_NUE_half_co2_plus100.RData")
load("./outputs/global_1x1_degree_2001_2015_NUE_half_Tair_plus1.RData")

## CARDAMOM run used to provide LAI and root C
cardamom = nc_open("/disk/scratch/local.2/lsmallma/Forest2020/C_cycle_analyses/DALEC_GSI_DFOL_CWD_FR_1_2001_2015_NEE_GPP_Rh_Ra_Bio_lit_cwd_som_timeseries.nc")
cardamom_latitude = ncvar_get(cardamom,"lat") ; cardamom_longitude = ncvar_get(cardamom,"long")

## Read in FLUXCOM GPP estimates for independent validation
list_of_files = list.files("/home/lsmallma/gcel/FLUXCOM/", full.name=TRUE)
list_of_files = list_of_files[grepl("GPP.annual",list_of_files)]
years_to_do = 2001:2013
GPP_types = c("ANN","MARS","RF")
for (y in seq(1,length(years_to_do))) {
    for (m in seq(1,length(GPP_types))) {
         tmp_list = list_of_files[grepl(years_to_do[y], list_of_files)]
         tmp_list = tmp_list[grepl(GPP_types[m], tmp_list)]
         datain = nc_open(tmp_list)
         tmp = ncvar_get(datain, "GPP")
         if (y == 1 & m == 1) {
             fluxcom_gpp = array(0, dim=c(dim(tmp),length(years_to_do)))
             fluxcom_gpp[,,y] = tmp * (1/length(GPP_types))
         } else {
             fluxcom_gpp[,,y] = fluxcom_gpp[,,y] + (tmp * (1/length(GPP_types)))
         }
    } # methods to loop
} # years to loop

## Updates for comparison later on with ACM-GPP-ET
# Convert units from gC/m2/day -> gC/m2/yr
fluxcom_gpp = fluxcom_gpp * 365.25
# Filter out zeros to NA
fluxcom_gpp[which(fluxcom_gpp == 0)] = NA
# Fix latitude orientation to match ACM-GPP-ET
fluxcom_gpp = fluxcom_gpp[,dim(fluxcom_gpp)[2]:1,]

## Regrid to 1 x 1 degree to match with ACM-GPP-ET
tmp_out_array = array(NA, dim=c(dim(global_output$mean_gpp),length(years_to_do)))
target_grid = raster(global_output$mean_gpp) 
for (y in seq(1, length(years_to_do))) {
     tmp_in = raster(fluxcom_gpp[,,y])    
     tmp_out = resample(tmp_in, target_grid, method="bilinear", filename="")
     tmp_out_array[,,y] = as.vector(t(tmp_out))
}
# now replace the original fluxcom at 0.5 degree to 1 degree
fluxcom_gpp = tmp_out_array ; rm(tmp_in,tmp_out,tmp_out_array,target_grid)
fluxcom_gpp_total = array(0,dim=dim(fluxcom_gpp)[1:2])
for (t in seq(1,dim(fluxcom_gpp)[3])) {fluxcom_gpp_total = fluxcom_gpp_total + fluxcom_gpp[,,t]} ; fluxcom_gpp_total = fluxcom_gpp_total / dim(fluxcom_gpp)[3]

###
## Some statistical / summary values
###

print("FLUXCOM output: ")
print(paste("       GPP = ",round(sum(fluxcom_gpp_total*global_output_NUE_half$area*1e-15,na.rm=TRUE),digits=1)," PgC",sep=""))

print("Global output: ")
print(paste("       GPP = ",round(global_output$global_mean_annual_gpp,digits=1)," PgC",sep=""))
print(paste("       Transpiration = ",round(global_output$global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
print(paste("       Soil Evaporation = ",round(global_output$global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
print(paste("       Wet Canopy Evaporation = ",round(global_output$global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
print(paste("       ET = ",round(global_output$global_mean_annual_et,digits=1)," PgH2O",sep=""))
print(paste("       WUE = ",round(global_output$global_mean_annual_wue,digits=2)," gC/kgH2O",sep=""))
print(paste("       wSWP = ",round(global_output$global_mean_annual_wSWP,digits=1)," MPa",sep=""))
print(paste("       Water in root zone = ",round(global_output$global_mean_annual_rootwatermm,digits=1),sep=""))

print("Global output: co2_plus100")
print(paste("       GPP = ",round(global_output_co2_plus100$global_mean_annual_gpp,digits=1)," PgC",sep=""))
print(paste("       Transpiration = ",round(global_output_co2_plus100$global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
print(paste("       Soil Evaporation = ",round(global_output_co2_plus100$global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
print(paste("       Wet Canopy Evaporation = ",round(global_output_co2_plus100$global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
print(paste("       ET = ",round(global_output_co2_plus100$global_mean_annual_et,digits=1)," PgH2O",sep=""))
print(paste("       WUE = ",round(global_output_co2_plus100$global_mean_annual_wue,digits=2)," gC/kgH2O",sep=""))
print(paste("       wSWP = ",round(global_output_co2_plus100$global_mean_annual_wSWP,digits=1)," MPa",sep=""))
print(paste("       Water in root zone = ",round(global_output_co2_plus100$global_mean_annual_rootwatermm,digits=1),sep=""))

print("Global output: NUE_half")
print(paste("       GPP = ",round(global_output_NUE_half$global_mean_annual_gpp,digits=1)," PgC",sep=""))
print(paste("       Transpiration = ",round(global_output_NUE_half$global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
print(paste("       Soil Evaporation = ",round(global_output_NUE_half$global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
print(paste("       Wet Canopy Evaporation = ",round(global_output_NUE_half$global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
print(paste("       ET = ",round(global_output_NUE_half$global_mean_annual_et,digits=1)," PgH2O",sep=""))
print(paste("       WUE = ",round(global_output_NUE_half$global_mean_annual_wue,digits=2)," gC/kgH2O",sep=""))
print(paste("       wSWP = ",round(global_output_NUE_half$global_mean_annual_wSWP,digits=1)," MPa",sep=""))
print(paste("       Water in root zone = ",round(global_output_NUE_half$global_mean_annual_rootwatermm,digits=1),sep=""))

print("Global output: NUE_half_co2_plus100")
print(paste("       GPP = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_gpp,digits=1)," PgC",sep=""))
print(paste("       Transpiration = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
print(paste("       Soil Evaporation = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
print(paste("       Wet Canopy Evaporation = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
print(paste("       ET = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_et,digits=1)," PgH2O",sep=""))
print(paste("       WUE = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_wue,digits=2)," gC/kgH2O",sep=""))
print(paste("       wSWP = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_wSWP,digits=1)," MPa",sep=""))
print(paste("       Water in root zone = ",round(global_output_NUE_half_co2_plus100$global_mean_annual_rootwatermm,digits=1),sep=""))

print("Global output: NUE_half_Tair_plus1")
print(paste("       GPP = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_gpp,digits=1)," PgC",sep=""))
print(paste("       Transpiration = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
print(paste("       Soil Evaporation = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
print(paste("       Wet Canopy Evaporation = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
print(paste("       ET = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_et,digits=1)," PgH2O",sep=""))
print(paste("       WUE = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_wue,digits=2)," gC/kgH2O",sep=""))
print(paste("       wSWP = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_wSWP,digits=1)," MPa",sep=""))
print(paste("       Water in root zone = ",round(global_output_NUE_half_Tair_plus1$global_mean_annual_rootwatermm,digits=1),sep=""))

###
## GPP mean annual correlations

nos_years = 15 ; steps_per_years = dim(global_output_NUE_half$timeseries_gpp)[3] / nos_years
annual_gpp = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2],nos_years))
annual_et = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2],nos_years))
annual_wue = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2],nos_years))
annual_lai = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2],nos_years))
annual_root = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2],nos_years))

r2_gpp_temperature = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_gpp_radiation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_gpp_vpd = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_gpp_precipitation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_et_temperature = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_et_radiation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_et_vpd = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_et_precipitation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_wue_temperature = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_wue_radiation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_wue_vpd = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))
r2_wue_precipitation = array(NA,dim=c(dim(global_output_NUE_half$timeseries_gpp)[1:2]))

## global_output_NUE_half
# Mean annual conditions
for (slot_i in seq(1,dim(global_output_NUE_half$timeseries_gpp)[1])) {
     for (slot_j in seq(1,dim(global_output_NUE_half$timeseries_gpp)[2])) {
          # calculate annual averages needed
          a = 1 ; b = steps_per_years
          for (y in seq(1,nos_years)) {
               annual_gpp[slot_i,slot_j,y] = mean(global_output_NUE_half$timeseries_gpp[slot_i,slot_j,a:b],na.rm=TRUE)
               annual_et[slot_i,slot_j,y] = mean(global_output_NUE_half$timeseries_transpiration[slot_i,slot_j,a:b] + global_output_NUE_half$timeseries_wetcanopyevap[slot_i,slot_j,a:b] + global_output_NUE_half$timeseries_soilevaporation[slot_i,slot_j,a:b],na.rm=TRUE)
               annual_wue[slot_i,slot_j,y] = mean(global_output_NUE_half$timeseries_WUE[slot_i,slot_j,a:b],na.rm=TRUE)
               annual_lai[slot_i,slot_j,y] = mean(global_output_NUE_half$timeseries_lai[slot_i,slot_j,a:b],na.rm=TRUE)
               annual_root[slot_i,slot_j,y] = mean(global_output_NUE_half$timeseries_root[slot_i,slot_j,a:b],na.rm=TRUE)
               a = a + steps_per_years ; b = b + steps_per_years
          } # time
	  if (length(which(is.na(annual_gpp[slot_i,slot_j,]) == FALSE)) > 0) {
              ## then directly calculate their correlation with temperature, radiation, VPD and precipitation
              # GPP
              r2_gpp_temperature[slot_i,slot_j] = summary(lm(annual_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_temperature[slot_i,slot_j,]))$adj.r.squared
              r2_gpp_radiation[slot_i,slot_j] = summary(lm(annual_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_radiation[slot_i,slot_j,]))$adj.r.squared
              r2_gpp_vpd[slot_i,slot_j] = summary(lm(annual_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_vpd[slot_i,slot_j,]))$adj.r.squared
              r2_gpp_precipitation[slot_i,slot_j] = summary(lm(annual_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_precipitation[slot_i,slot_j,]))$adj.r.squared
              # ET
              r2_et_temperature[slot_i,slot_j] = summary(lm(annual_et[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_temperature[slot_i,slot_j,]))$adj.r.squared
              r2_et_radiation[slot_i,slot_j] = summary(lm(annual_et[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_radiation[slot_i,slot_j,]))$adj.r.squared
              r2_et_vpd[slot_i,slot_j] = summary(lm(annual_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_vpd[slot_i,slot_j,]))$adj.r.squared
              r2_et_precipitation[slot_i,slot_j] = summary(lm(annual_et[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_precipitation[slot_i,slot_j,]))$adj.r.squared
              # ET
              r2_wue_temperature[slot_i,slot_j] = summary(lm(annual_wue[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_temperature[slot_i,slot_j,]))$adj.r.squared
              r2_wue_radiation[slot_i,slot_j] = summary(lm(annual_wue[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_radiation[slot_i,slot_j,]))$adj.r.squared
              r2_wue_vpd[slot_i,slot_j] = summary(lm(annual_wue[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_vpd[slot_i,slot_j,]))$adj.r.squared
              r2_wue_precipitation[slot_i,slot_j] = summary(lm(annual_wue[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_precipitation[slot_i,slot_j,]))$adj.r.squared
          }
     } # j
} # i

fig_height=4000 ; fig_width=7200
my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral")))) 
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
jpeg(file="./FIGURES/Cal_val_paper_figure_S3.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,4), mar=c(1.4, 4.2, 2.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
var1=(r2_gpp_temperature) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main=expression(paste("Mean temperature (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
mtext("GPP",side = 2,cex=1.8, padj = 1.2)
var1=(r2_gpp_radiation) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main=expression(paste("Mean Radiation (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(r2_gpp_precipitation) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main=expression(paste("Mean precipitation (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(r2_gpp_vpd) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main=expression(paste("Mean VPD (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(r2_et_temperature) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
mtext("ET",side = 2,cex=1.8, padj = 1.2)
var1=(r2_et_radiation) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(r2_et_precipitation) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(r2_et_vpd) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(r2_wue_temperature) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
mtext("WUE",side = 2,cex=1.8, padj = 1.2)
var1=(r2_wue_radiation) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(r2_wue_precipitation) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(r2_wue_vpd) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main="", col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

fluxcom_r2_gpp_temperature = array(NA,dim=c(dim(fluxcom_gpp)[1:2]))
fluxcom_r2_gpp_radiation = array(NA,dim=c(dim(fluxcom_gpp)[1:2]))
fluxcom_r2_gpp_vpd = array(NA,dim=c(dim(fluxcom_gpp)[1:2]))
fluxcom_r2_gpp_precipitation = array(NA,dim=c(dim(fluxcom_gpp)[1:2]))

## FLUXCOM (2001-2013)
# Mean annual conditions
for (slot_i in seq(1,dim(global_output_NUE_half$timeseries_gpp)[1])) {
     for (slot_j in seq(1,dim(global_output_NUE_half$timeseries_gpp)[2])) {
          if (length(which(is.na(fluxcom_gpp[slot_i,slot_j,]) == FALSE)) > 0 & length(which(is.na(global_output_NUE_half$mean_annual_temperature[slot_i,slot_j,]) == FALSE)) > 0) {
              ## then directly calculate their correlation with temperature, radiation, VPD and precipitation
              fluxcom_r2_gpp_temperature[slot_i,slot_j] = summary(lm(fluxcom_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_temperature[slot_i,slot_j,1:13]))$adj.r.squared
              fluxcom_r2_gpp_radiation[slot_i,slot_j] = summary(lm(fluxcom_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_radiation[slot_i,slot_j,1:13]))$adj.r.squared
              fluxcom_r2_gpp_vpd[slot_i,slot_j] = summary(lm(fluxcom_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_vpd[slot_i,slot_j,1:13]))$adj.r.squared
              fluxcom_r2_gpp_precipitation[slot_i,slot_j] = summary(lm(fluxcom_gpp[slot_i,slot_j,] ~ global_output_NUE_half$mean_annual_precipitation[slot_i,slot_j,1:13]))$adj.r.squared
          }
     } # j
} # i

fig_height=4000 ; fig_width=7200
my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral")))) 
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
jpeg(file="./FIGURES/Cal_val_paper_figure_S4.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,2), mar=c(1.4, 4.4, 2.8, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
var1=(fluxcom_r2_gpp_temperature) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main=expression(paste("Mean temperature (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
mtext("FLUXCOM GPP",side = 2,cex=1.8, padj = 1.4)
var1=(fluxcom_r2_gpp_radiation) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main=expression(paste("Mean Radiation (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(fluxcom_r2_gpp_precipitation) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main=expression(paste("Mean precipitation (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(fluxcom_r2_gpp_vpd) ; colour_choices=colour_choices_upper(length(var1))
image.plot(var1, main=expression(paste("Mean VPD (",R^2,")",sep="")), col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

###
## Begin output of figures
###

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
# calculate temperature and precipitation means for each site
site_mean_temperature = rep(NA, length(latitude)) ; site_mean_rainfall = rep(NA, length(latitude))
for (n in seq(1,length(latitude))) {
    tmp = which(data$Lat == latitude[n] & data$Lon == longitude[n])
    site_mean_temperature[n] = mean(as.vector(unlist(data[tmp,7:14])))-273.15 # K->C
    site_mean_rainfall[n] = mean(as.vector(unlist(data[tmp,31:38]))) # kg/m2/s
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
plot(longitude,latitude, xlab="", ylab="", pch=16,cex=1.0,xaxt = "n", yaxt = "n",ylim=c(1,179), xlim=c(1,359), xaxs="i", yaxs="i")
par(fig=c(0.70,1.0,0.0,1.0),mar=c(4.6, 4.6, 0.4, 0.5), omi=c(0.2, 0.2, 0.2, 0.4), new=TRUE)
var1 = as.vector(global_output$mean_temperature)
var2 = as.vector(global_output$mean_precipitation*(365.25*60*60*24))
smoothScatter(var1,var2, ylim=c(0,5000), 
              pch=16, cex.axis=1.6, cex.lab=1.6, cex.main=1.6,nrpoints=0,colramp=my_colours,  
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.015,diff(range(var2,na.rm=TRUE))*0.0025),nbin=128*10,
              ylab="Mean precipitation (kgH2O/m/yr)", xlab="Mean air temperature (Celcuis)")
points(site_mean_temperature, site_mean_rainfall*(365.25*60*60*24), pch=16, col="black")
dev.off()

fig_height=4000 ; fig_width=6000
jpeg(file="./FIGURES/Cal_val_paper_figure_2.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,2), mar=c(4.4, 4.2, 4.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
plot(calibration_output$drivers$GPP~calibration_output$mean_gpp,main="GPP", ylab="SPA", xlab="ACM-GPP-ET",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,1,col="red",lwd=3)
transpiration = (calibration_output$drivers$Evap-calibration_output$drivers$soilevap-calibration_output$drivers$wetevap)
plot(transpiration~calibration_output$mean_transpiration,main="Transpiration", ylab="SPA", xlab="ACM-GPP-ET",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,1,col="red",lwd=3)
plot(calibration_output$drivers$soilevap~calibration_output$mean_soilevaporation,main="Soil evaporation", ylab="SPA", xlab="ACM-GPP-ET",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,1,col="red",lwd=3)
plot(calibration_output$drivers$wetevap~calibration_output$mean_wetcanopyevap,main="Wet canopy evaporation", ylab="SPA", xlab="ACM-GPP-ET",pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6) ; abline(0,1,col="red",lwd=3)
dev.off()

my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral")))) 
fig_height=4000 ; fig_width=6000
jpeg(file="./FIGURES/Cal_val_paper_figure_2_heat_map.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,2), mar=c(4.4, 4.2, 2.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
var1 = calibration_output$mean_gpp
var2 = calibration_output$drivers$GPP
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,
              main=expression(paste("GPP"," (gC ",m^-2," da",y^-1,")")), ylab="SPA", xlab="ACM-GPP-ET",cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
var1 = calibration_output$mean_transpiration
var2 = (calibration_output$drivers$Evap-calibration_output$drivers$soilevap-calibration_output$drivers$wetevap)
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main=expression(paste("Transpiration"," (kgH2O ",m^-2," da",y^-1,")")), ylab="SPA", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
var1 = calibration_output$mean_soilevaporation
var2 = calibration_output$drivers$soilevap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main=expression(paste("Soil evaporation"," (kgH2O ",m^-2," da",y^-1,")")), ylab="SPA", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
var1 = calibration_output$mean_wetcanopyevap
var2 = calibration_output$drivers$wetevap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main=expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," da",y^-1,")")), ylab="SPA", xlab="ACM-GPP-ET", 
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
dev.off()

my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral"))))
fig_height=6000 ; fig_width=4000
jpeg(file="./FIGURES/Cal_val_paper_figure_2a_heat_map.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(4,2), mar=c(4.4, 4.2, 2.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
var1 = calibration_output$mean_gpp
var2 = calibration_output$drivers$GPP
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,
              main="Calibration", ylab="SPA", xlab="",cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
mtext(expression(paste("GPP"," (gC ",m^-2," da",y^-1,")")),side = 3,cex=1.8, padj = -0.2, adj = 2.3)
var1 = validation_water_output$mean_gpp
var2 = validation_water_output$drivers$GPP
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,
              main="Validation", ylab="", xlab="",cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)


var1 = calibration_output$mean_transpiration
var2 = (calibration_output$drivers$Evap-calibration_output$drivers$soilevap-calibration_output$drivers$wetevap)
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="SPA", xlab="",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
mtext(expression(paste("Transpiration"," (kgH2O ",m^-2," da",y^-1,")")),side = 3,cex=1.8, padj = -0.2, adj = 2.3)
var1 = validation_water_output$mean_transpiration
var2 = (validation_water_output$drivers$Evap-validation_water_output$drivers$soilevap-validation_water_output$drivers$wetevap)
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)


var1 = calibration_output$mean_soilevaporation
var2 = calibration_output$drivers$soilevap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="SPA", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
mtext(expression(paste("Soil evaporation"," (kgH2O ",m^-2," da",y^-1,")")),side = 3,cex=1.8, padj = -0.2, adj = 2.3)
var1 = validation_water_output$mean_soilevaporation
var2 = validation_water_output$drivers$soilevap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)


var1 = calibration_output$mean_wetcanopyevap
var2 = calibration_output$drivers$wetevap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="SPA", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
mtext(expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," da",y^-1,")")),side = 3,cex=1.8, padj = -0.2, adj = 2.3)
var1 = validation_water_output$mean_wetcanopyevap
var2 = validation_water_output$drivers$wetevap
smoothScatter(var1,var2,nrpoints=0,colramp=my_colours,main="", ylab="SPA", xlab="ACM-GPP-ET",
              cex=0.5,pch=16,cex.axis=1.6,cex.lab=1.6,cex.main=1.6,
              transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10) ; abline(0,1,col="red",lwd=3)
dev.off()

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_3.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(global_output_NUE_half$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25 ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("ET"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=global_output_NUE_half$mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("WUE (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# ET component figures transpiration, wetcanopy evaporation, soil evaporation
var1=(global_output_NUE_half$mean_transpiration*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Transpiration"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_wetcanopyevap*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_soilevaporation*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# Soil runoff / drainage and status
var1=(global_output_NUE_half$mean_runoffmm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil run-off"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_drainagemm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil drainage"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_rootwatermm) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Water in rooting zone"," (",kg,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

#fluxcom_gpp_mean = array(0,dim=dim(fluxcom_gpp)[1:2])
#for (t in seq(1,dim(fluxcom_gpp)[3])) {fluxcom_gpp_mean = fluxcom_gpp_mean + fluxcom_gpp[,,t]} ; fluxcom_gpp_mean = fluxcom_gpp_mean / dim(fluxcom_gpp)[3]
#fluxcom_gpp_difference = (global_output_NUE_half$mean_gpp*365.25) - fluxcom_gpp_mean
#image.plot(fluxcom_gpp_difference)

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_4.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(global_output_NUE_half$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25 ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("ET"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=global_output_NUE_half$mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("WUE (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# ET component figures transpiration, wetcanopy evaporation, soil evaporation
var1=(global_output_NUE_half$mean_transpiration*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Transpiration"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_wetcanopyevap*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_soilevaporation*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# Soil runoff / drainage and status
var1=(global_output_NUE_half$mean_runoffmm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil run-off"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_drainagemm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil drainage"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_rootwatermm) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Water in rooting zone"," (",kg,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_5.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,3), mar=c(1.4, 1.2, 2.4, 8.5), omi=c(0.1, 0.1, 0.1, 0.1))
## Summary figures GPP, ET, WUE; CO2 + 100
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=((global_output_NUE_half_co2_plus100$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25))/(global_output_NUE_half$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
image.plot(var1, main=expression(paste(Delta,"GPP"," (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
mtext(expression(paste("C",O[2],"+100 ppm")),side = 2,cex=1.8, padj = 1.2)
var1=((global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25)
var1 = var1-((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var1 = var1 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)  
colour_choices=colour_choices_upper(length(var1))
zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
image.plot(var1, main=expression(paste(Delta,"ET"," (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_co2_plus100$mean_wue - global_output_NUE_half$mean_wue) / global_output_NUE_half$mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
image.plot(var1, main=expression(paste(Delta,"WUE (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
## Summary figures GPP, ET, WUE; Tair + 1
var1=((global_output_NUE_half_Tair_plus1$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25))/(global_output_NUE_half$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
image.plot(var1, main=expression(paste(Delta,"GPP"," (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
mtext(expression(paste("Tair+",1^o,"C",sep="")),side = 2,cex=1.8, padj = 1.2)
var1=((global_output_NUE_half_Tair_plus1$mean_transpiration+global_output_NUE_half_Tair_plus1$mean_wetcanopyevap+global_output_NUE_half_Tair_plus1$mean_soilevaporation)*365.25)
var1 = var1-((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var1 = var1 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)  
colour_choices=colour_choices_upper(length(var1))
zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
image.plot(var1, main=expression(paste(Delta,"ET"," (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1 = (global_output_NUE_half_Tair_plus1$mean_wue - global_output_NUE_half$mean_wue) / global_output_NUE_half$mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=max(abs(quantile(var1,prob=c(0.001,0.999),na.rm=TRUE))) ; zaxis = c(-1*zaxis,zaxis)
image.plot(var1, main=expression(paste(Delta,"WUE (Fraction)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_6.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,4), mar=c(4.4, 5.2, 2.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
var1 = global_output_NUE_half$mean_temperature
var2 = (((global_output_NUE_half_co2_plus100$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25) 
plot(var1,var2, ylab=expression(paste(Delta,"GPP"," (fraction)")),xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_radiation
var2 = (((global_output_NUE_half_co2_plus100$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
plot(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
mtext(expression(paste("C",O[2],"+100 ppm")),side = 3,cex=2.2, padj = -0.2, adj = 2.3)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = (((global_output_NUE_half_co2_plus100$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
plot(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_vpd
var2 = (((global_output_NUE_half_co2_plus100$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
plot(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_temperature
var2 = ((global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
plot(var1,var2, ylab=expression(paste(Delta,"ET"," (fraction)")), xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_radiation
var2 = ((global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
plot(var1,var2, ylab="", xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = ((global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
plot(var1,var2, ylab="", xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_vpd
var2 = ((global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
plot(var1,var2, ylab="", xlab="Mean VPD (kPa)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_temperature
var2 = (global_output_NUE_half_co2_plus100$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
plot(var1,var2, ylab=expression(paste(Delta,"WUE (fraction)")), xlab="Mean Temperature (Celcius)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_radiation
var2 = (global_output_NUE_half_co2_plus100$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
plot(var1,var2, ylab="", xlab=expression(paste("Mean Radiation (MJ/",m^2,"/day)",sep="")), main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = (global_output_NUE_half_co2_plus100$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
plot(var1,var2, ylab="", xlab=expression(paste("Mean water in root zone (kg",H[2],"O/",m^2,")",sep="")), main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_vpd
var2 = (global_output_NUE_half_co2_plus100$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
plot(var1,var2, ylab="", xlab="Mean VPD (kPa)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
dev.off()

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_7.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,4), mar=c(4.4, 5.2, 2.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
var1 = global_output_NUE_half$mean_temperature
var2 = (((global_output_NUE_half_Tair_plus1$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25) 
plot(var1,var2, ylab=expression(paste(Delta,"GPP"," (fraction)")),xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_radiation
var2 = (((global_output_NUE_half_Tair_plus1$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
plot(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
mtext(expression(paste("Tair+",1^o,"C")),side = 3,cex=2.2, padj = -0.2, adj = 1.5)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = (((global_output_NUE_half_Tair_plus1$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
plot(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_vpd
var2 = (((global_output_NUE_half_Tair_plus1$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
plot(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_temperature
var2 = ((global_output_NUE_half_Tair_plus1$mean_transpiration+global_output_NUE_half_Tair_plus1$mean_wetcanopyevap+global_output_NUE_half_Tair_plus1$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
plot(var1,var2, ylab=expression(paste(Delta,"ET"," (fraction)")), xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_radiation
var2 = ((global_output_NUE_half_Tair_plus1$mean_transpiration+global_output_NUE_half_Tair_plus1$mean_wetcanopyevap+global_output_NUE_half_Tair_plus1$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
plot(var1,var2, ylab="", xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = ((global_output_NUE_half_Tair_plus1$mean_transpiration+global_output_NUE_half_Tair_plus1$mean_wetcanopyevap+global_output_NUE_half_Tair_plus1$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
plot(var1,var2, ylab="", xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_vpd
var2 = ((global_output_NUE_half_Tair_plus1$mean_transpiration+global_output_NUE_half_Tair_plus1$mean_wetcanopyevap+global_output_NUE_half_Tair_plus1$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
plot(var1,var2, ylab="", xlab="Mean VPD (kPa)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_temperature
var2 = (global_output_NUE_half_Tair_plus1$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
plot(var1,var2, ylab=expression(paste(Delta,"WUE (fraction)")), xlab="Mean Temperature (Celcius)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_radiation
var2 = (global_output_NUE_half_Tair_plus1$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
plot(var1,var2, ylab="", xlab=expression(paste("Mean Radiation (MJ/",m^2,"/day)",sep="")), main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = (global_output_NUE_half_Tair_plus1$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
plot(var1,var2, ylab="", xlab=expression(paste("Mean water in root zone (kg",H[2],"O/",m^2,")",sep="")), main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
var1 = global_output_NUE_half$mean_vpd
var2 = (global_output_NUE_half_Tair_plus1$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
plot(var1,var2, ylab="", xlab="Mean VPD (kPa)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8)
dev.off()

my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral")))) 
fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_6_heat_map.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,4), mar=c(4.4, 5.2, 2.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
var1 = global_output_NUE_half$mean_temperature
var2 = (((global_output_NUE_half_co2_plus100$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25) 
smoothScatter(var1,var2, ylab=expression(paste(Delta,"GPP"," (fraction)")),xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_radiation
var2 = (((global_output_NUE_half_co2_plus100$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
smoothScatter(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
mtext(expression(paste("C",O[2],"+100 ppm")),side = 3,cex=2.2, padj = -0.2, adj = 2.3)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = (((global_output_NUE_half_co2_plus100$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
smoothScatter(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_vpd
var2 = (((global_output_NUE_half_co2_plus100$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
smoothScatter(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_temperature
var2 = ((global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
smoothScatter(var1,var2, ylab=expression(paste(Delta,"ET"," (fraction)")), xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_radiation
var2 = ((global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
smoothScatter(var1,var2, ylab="", xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = ((global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
smoothScatter(var1,var2, ylab="", xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_vpd
var2 = ((global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
smoothScatter(var1,var2, ylab="", xlab="Mean VPD (kPa)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_temperature
var2 = (global_output_NUE_half_co2_plus100$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
smoothScatter(var1,var2, ylab=expression(paste(Delta,"WUE (fraction)")), xlab="Mean Temperature (Celcius)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_radiation
var2 = (global_output_NUE_half_co2_plus100$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
smoothScatter(var1,var2, ylab="", xlab=expression(paste("Mean Radiation (MJ/",m^2,"/day)",sep="")), main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = (global_output_NUE_half_co2_plus100$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
smoothScatter(var1,var2, ylab="", xlab=expression(paste("Mean water in root zone (kg",H[2],"O/",m^2,")",sep="")), main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_vpd
var2 = (global_output_NUE_half_co2_plus100$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
smoothScatter(var1,var2, ylab="", xlab="Mean VPD (kPa)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
dev.off()

my_colours=colorRampPalette(c("white",rev(brewer.pal(11,"Spectral")))) 
fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_7_heat_map.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,4), mar=c(4.4, 5.2, 2.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
var1 = global_output_NUE_half$mean_temperature
var2 = (((global_output_NUE_half_Tair_plus1$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25) 
smoothScatter(var1,var2, ylab=expression(paste(Delta,"GPP"," (fraction)")),xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_radiation
var2 = (((global_output_NUE_half_Tair_plus1$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
smoothScatter(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
mtext(expression(paste("Tair+",1^o,"C")),side = 3,cex=2.2, padj = -0.2, adj = 1.5)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = (((global_output_NUE_half_Tair_plus1$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
smoothScatter(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_vpd
var2 = (((global_output_NUE_half_Tair_plus1$mean_gpp*365.25)-(global_output_NUE_half$mean_gpp*365.25)))/(global_output_NUE_half$mean_gpp*365.25)
smoothScatter(var1,var2, ylab="",xlab="", main="", cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_temperature
var2 = ((global_output_NUE_half_Tair_plus1$mean_transpiration+global_output_NUE_half_Tair_plus1$mean_wetcanopyevap+global_output_NUE_half_Tair_plus1$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
smoothScatter(var1,var2, ylab=expression(paste(Delta,"ET"," (fraction)")), xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_radiation
var2 = ((global_output_NUE_half_Tair_plus1$mean_transpiration+global_output_NUE_half_Tair_plus1$mean_wetcanopyevap+global_output_NUE_half_Tair_plus1$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
smoothScatter(var1,var2, ylab="", xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = ((global_output_NUE_half_Tair_plus1$mean_transpiration+global_output_NUE_half_Tair_plus1$mean_wetcanopyevap+global_output_NUE_half_Tair_plus1$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
smoothScatter(var1,var2, ylab="", xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_vpd
var2 = ((global_output_NUE_half_Tair_plus1$mean_transpiration+global_output_NUE_half_Tair_plus1$mean_wetcanopyevap+global_output_NUE_half_Tair_plus1$mean_soilevaporation)*365.25)
var2 = var2 - ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
var2 = var2 / ((global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25)
smoothScatter(var1,var2, ylab="", xlab="Mean VPD (kPa)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_temperature
var2 = (global_output_NUE_half_Tair_plus1$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
smoothScatter(var1,var2, ylab=expression(paste(Delta,"WUE (fraction)")), xlab="Mean Temperature (Celcius)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_radiation
var2 = (global_output_NUE_half_Tair_plus1$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
smoothScatter(var1,var2, ylab="", xlab=expression(paste("Mean Radiation (MJ/",m^2,"/day)",sep="")), main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_rootwatermm #(as.vector(global_output_NUE_half$mean_precipitation)*365.25*60*60*24)
var2 = (global_output_NUE_half_Tair_plus1$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
smoothScatter(var1,var2, ylab="", xlab=expression(paste("Mean water in root zone (kg",H[2],"O/",m^2,")",sep="")), main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_vpd
var2 = (global_output_NUE_half_Tair_plus1$mean_wue - global_output_NUE_half$mean_wue)/global_output_NUE_half$mean_wue 
smoothScatter(var1,var2, ylab="", xlab="Mean VPD (kPa)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16, cex.lab=1.8,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
dev.off()

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_S1.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,4), mar=c(2.4, 4.2, 4.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
var1=global_output_NUE_half$mean_gpp*365.25
plot(var1~global_output_NUE_half$mean_temperature, ylab=expression(paste("GPP"," (gC/m2/yr)")),xlab="Mean Temperature (Celcius)", main="", cex=1.8, cex.axis=2.0, pch=16)
plot(var1~global_output_NUE_half$mean_radiation, ylab=expression(paste("GPP"," (gC/m2/yr)")),xlab="Mean Radiation (MJ/m2/day)", main="", cex=1.8, cex.axis=2.0, pch=16)
plot(var1~global_output_NUE_half$mean_precipitation, ylab=expression(paste("GPP"," (gC/m2/yr)")),xlab="Mean Precipitation (kgH2O/m2/s)", main="", cex=1.8, cex.axis=2.0, pch=16)
plot(var1~global_output_NUE_half$mean_vpd, ylab=expression(paste("GPP"," (gC/m2/yr)")),xlab="Mean VPD (kPa)", main="", cex=1.8, cex.axis=2.0, pch=16)
var1=(global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25
plot(var1~global_output_NUE_half$mean_temperature, ylab=expression(paste("ET"," (kgH2O/m2/yr)")), xlab="Mean Temperature (Celcius)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16)
plot(var1~global_output_NUE_half$mean_radiation, ylab=expression(paste("ET"," (kgH2O/m2/yr")), xlab="Mean Radiation (MJ/m2/day)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16)
plot(var1~global_output_NUE_half$mean_precipitation, ylab=expression(paste("ET"," (kgH2O/m2/yr)")), xlab="Mean Precipitaion (kgH2O/m2/s)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16)
plot(var1~global_output_NUE_half$mean_vpd, ylab=expression(paste("ET"," (kgH2O/m2/yr)")), xlab="Mean VPD (kPa)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16)
var1=global_output_NUE_half$mean_wue
plot(var1~global_output_NUE_half$mean_temperature, ylab=expression(paste("WUE (gC/kgH2O)")), xlab="Mean Temperature (Celcius)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16)
plot(var1~global_output_NUE_half$mean_radiation, ylab=expression(paste("WUE (gC/kgH2O)")), xlab="Mean Radiation (MJ/m2/day)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16)
plot(var1~global_output_NUE_half$mean_precipitation, ylab=expression(paste("WUE (gC/kgH2O)")), xlab="Mean Precipitation (kgH2O/m2/s)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16)
plot(var1~global_output_NUE_half$mean_vpd, ylab=expression(paste("WUE (gC/kgH2O)")), xlab="Mean VPD (kPa)", main="", cex.main=2.4, cex=1.8, cex.axis=2.0, pch=16)
dev.off()

fig_height=4000 ; fig_width=6000
jpeg(file="./FIGURES/Cal_val_paper_figure_S2.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,2), mar=c(4.4, 5.6, 1.4, 1.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
var1 = global_output_NUE_half$mean_lai
var2 = global_output_NUE_half$mean_gpp*365.25
smoothScatter(var1,var2, ylab=expression(paste("GPP"," (gC/m2/yr)")),xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.2, pch=16, cex.lab=2.2, xlim=c(0,6),
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_root
smoothScatter(var1,var2, ylab="",xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.2, pch=16, cex.lab=2.2,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_lai
var2=(global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25
smoothScatter(var1,var2, ylab=expression(paste("ET"," (kgH2O/m2/yr)")), xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.2, pch=16, cex.lab=2.2, xlim=c(0,6),
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_root
smoothScatter(var1,var2, ylab="", xlab="", main="", cex.main=2.4, cex=1.8, cex.axis=2.2, pch=16, cex.lab=2.2,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_lai
var2=global_output_NUE_half$mean_wue 
smoothScatter(var1,var2, ylab=expression(paste("WUE (gC/kgH2O)")), xlab=expression(paste("Mean LAI (",m^2,"/",m^2,")",sep="")), main="", cex.main=2.4, cex=1.8, cex.axis=2.2, pch=16, cex.lab=2.2, xlim=c(0,6),
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
var1 = global_output_NUE_half$mean_root
smoothScatter(var1,var2, ylab="", xlab=expression(paste("Mean Root (gC/",m^2,")",sep="")), main="", cex.main=2.4, cex=1.8, cex.axis=2.2, pch=16, cex.lab=2.2,
              nrpoints=0,colramp=my_colours,transformation = function(x) x^.25, bandwidth=c(diff(range(var1,na.rm=TRUE))*0.0075,diff(range(var2,na.rm=TRUE))*0.0075),nbin=128*10)
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

###
## Global_co2_plus100
###

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/All_fluxers_global_output_co2_plus100.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(global_output_co2_plus100$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_co2_plus100$mean_transpiration+global_output_co2_plus100$mean_wetcanopyevap+global_output_co2_plus100$mean_soilevaporation)*365.25 ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("ET"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=global_output_co2_plus100$mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("WUE (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# ET component figures transpiration, wetcanopy evaporation, soil evaporation
var1=(global_output_co2_plus100$mean_transpiration*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Transpiration"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_co2_plus100$mean_wetcanopyevap*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_co2_plus100$mean_soilevaporation*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# Soil runoff / drainage and status
var1=(global_output_co2_plus100$mean_runoffmm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil run-off"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_co2_plus100$mean_drainagemm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil drainage"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_co2_plus100$mean_rootwatermm) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Water in rooting zone"," (",kg,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/All_fluxes_difference_global_co2_minus_global.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=((global_output_co2_plus100$mean_gpp*365.25)-(global_output$mean_gpp*365.25)) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste(Delta,"GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=((global_output_co2_plus100$mean_transpiration+global_output_co2_plus100$mean_wetcanopyevap+global_output_co2_plus100$mean_soilevaporation)*365.25)
var1 = var1-((global_output$mean_transpiration+global_output$mean_wetcanopyevap+global_output$mean_soilevaporation)*365.25) 
colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste(Delta,"ET"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=global_output_co2_plus100$mean_wue - global_output$mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste(Delta,"WUE (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# ET component figures transpiration, wetcanopy evaporation, soil evaporation
var1=(global_output_co2_plus100$mean_transpiration*365.25)-(global_output$mean_transpiration*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste(Delta,"Transpiration"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_co2_plus100$mean_wetcanopyevap*365.25)-(global_output$mean_wetcanopyevap*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(-0.1,0.1)#quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste(Delta,"Wet canopy evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_co2_plus100$mean_soilevaporation*365.25)-(global_output$mean_soilevaporation*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste(Delta,"Soil evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# Soil runoff / drainage and status
var1=(global_output_co2_plus100$mean_runoffmm*365.25)-(global_output$mean_runoffmm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste(Delta,"Soil run-off"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_co2_plus100$mean_drainagemm*365.25)-(global_output$mean_drainagemm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste(Delta,"Soil drainage"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_co2_plus100$mean_rootwatermm)-(global_output$mean_rootwatermm) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste(Delta,"Water in rooting zone"," (",kg,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

###
## Global output NUE half

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/All_fluxes_global_output_NUE_half.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(global_output_NUE_half$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_transpiration+global_output_NUE_half$mean_wetcanopyevap+global_output_NUE_half$mean_soilevaporation)*365.25 ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("ET"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=global_output_NUE_half$mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("WUE (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# ET component figures transpiration, wetcanopy evaporation, soil evaporation
var1=(global_output_NUE_half$mean_transpiration*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Transpiration"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_wetcanopyevap*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_soilevaporation*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# Soil runoff / drainage and status
var1=(global_output_NUE_half$mean_runoffmm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil run-off"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_drainagemm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil drainage"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half$mean_rootwatermm) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Water in rooting zone"," (",kg,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

###
## NUE half CO2 + 100

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/All_fluxers_global_output_NUE_half_co2_plus100.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(global_output_NUE_half_co2_plus100$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_co2_plus100$mean_transpiration+global_output_NUE_half_co2_plus100$mean_wetcanopyevap+global_output_NUE_half_co2_plus100$mean_soilevaporation)*365.25 ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("ET"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=global_output_NUE_half_co2_plus100$mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("WUE (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# ET component figures transpiration, wetcanopy evaporation, soil evaporation
var1=(global_output_NUE_half_co2_plus100$mean_transpiration*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Transpiration"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_co2_plus100$mean_wetcanopyevap*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_co2_plus100$mean_soilevaporation*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# Soil runoff / drainage and status
var1=(global_output_NUE_half_co2_plus100$mean_runoffmm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil run-off"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_co2_plus100$mean_drainagemm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil drainage"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_co2_plus100$mean_rootwatermm) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Water in rooting zone"," (",kg,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

###
## NUE half Tair + 1 


fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/All_fluxers_global_output_NUE_half_Tair_plus1.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(global_output_NUE_half_Tair_plus1$mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_Tair_plus1$mean_transpiration+global_output_NUE_half_Tair_plus1$mean_wetcanopyevap+global_output_NUE_half_Tair_plus1$mean_soilevaporation)*365.25 ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("ET"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=global_output_NUE_half_Tair_plus1$mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("WUE (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# ET component figures transpiration, wetcanopy evaporation, soil evaporation
var1=(global_output_NUE_half_Tair_plus1$mean_transpiration*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Transpiration"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_Tair_plus1$mean_wetcanopyevap*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_Tair_plus1$mean_soilevaporation*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# Soil runoff / drainage and status
var1=(global_output_NUE_half_Tair_plus1$mean_runoffmm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil run-off"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_Tair_plus1$mean_drainagemm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Soil drainage"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(global_output_NUE_half_Tair_plus1$mean_rootwatermm) ; colour_choices=colour_choices_upper(length(var1))
zaxis=quantile(var1,na.rm=TRUE,prob=c(0.001,0.999))
image.plot(var1, main=expression(paste("Water in rooting zone"," (",kg,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.4,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
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
