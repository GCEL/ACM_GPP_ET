###
## function to return x,y coordinate from an array which is nearest to a provided lat / long value
###

###
## Create needed ACM_GPP_ET shared object

# Packages
library(RColorBrewer)
# set to the working directory that this script should be called from
setwd("/home/lsmallma/WORK/GREENHOUSE/models/ACM_GPP_ET/") ; wkdir = getwd()
# compile the shared object containing ACM_GPP and ACM_ET
system("gfortran ./src/ACM_GPP_ET.f90 ./src/ACM_GPP_ET_R_interface.f90 -o ./src/acm_gpp_et.so -fPIC -shared")
system("mv ./src/acm_gpp_et.so .")

###
## Borrow met data from an existing CARDAMOM analysis

drivers = read.csv("/home/lsmallma/gcel/ACM_GPP_ET_RECALIBRATION/output_files/acm_recal_with_spa_200pixels_continuous_timeseries_obs_iWUE_trunk_nowater_copy.csv")

# determine the mean values for each driver which we will now repeat...
met=array(-9999,dim=c(1,12))
met[,1] = 1  # day of analysis
met[,2] = mean(drivers$sat_min)  # min temperature (oC)
met[,3] = mean(drivers$sat_max)  # max temperature (oC)
met[,4] = mean(drivers$swrad_avg)*(24*60*60*1e-6)  # SW Radiation (MJ.m-2.day-1)
met[,5] = mean(drivers$co2_avg)  # CO2 ppm
met[,6] = mean(drivers$doy)  # day of year
met[,7] = mean(drivers$ppt_avg)/(60*60)  # rainfall (kgH2O.m-2.hr-1 -> kg.m-2.s-1)
met[,8] = -9999 # drivers$sat_avg[n] # avg temperature (oC)
met[,9] = mean(drivers$wind_avg) # avg wind speed (m.s-1)
met[,10]= mean(drivers$vpd_avg)*1000 # avg VPD (kPa->Pa)

###
## Some ACM_GPP_ET parameters

output_dim=11 ; nofluxes = 8 ; nopools = 1 ; nopars = 4 ; nos_iter = 1 ; iterations = 100

lai_iterations = seq(0,10,length.out=iterations) # rep(10.0,times=iterations) # seq(0,10,length.out=iterations)
avN_iterations = rep(mean(drivers$avgN),times=iterations) # seq(0.5,5,length.out=iterations) # rep(mean(drivers$avgN),times=iterations)

###
## Define our output variables based on the grid of the CARDAMOM analysis we are borrowing

mean_lai = array(NA, dim=c(iterations))
mean_gpp = array(NA, dim=c(iterations))
mean_transpiration = array(NA, dim=c(iterations))
mean_wetcanopyevap = array(NA, dim=c(iterations))
mean_soilevaporation = array(NA, dim=c(iterations))
mean_rootwatermm = array(NA, dim=c(iterations))
mean_runoffmm = array(NA, dim=c(iterations))
mean_drainagemm = array(NA, dim=c(iterations))
mean_WUE = array(NA, dim=c(iterations))
mean_wSWP = array(NA, dim=c(iterations))

# iterative process through the years...
for (n in seq(1,iterations)) {

     if (n%%1000 == 0){print(paste("...beginning site:",n," of ",dim(drivers)[1], sep=""))}

     # ecosystem state drivers now rather than meteorology
     met[,11]= lai_iterations[n] # LAI
     met[,12]= 100 #100  # root C stocks
     # parameters
     parameters = array(NA, dim=c(nopars,nos_iter))
     parameters[1,] = avN_iterations[n]  # foliar N (gN.m-2)
     parameters[2,] = -9999 # min leaf water potential (MPa)
     parameters[3,] = 100   # root biomass needed to reach 50 % depth
     parameters[4,] = 2   # max root depth (m)

     # other inputs
     lat = mean(drivers$lat)
     # search location of soils data
#     i1=unlist(closest2d(1,soils_data$lat_wanted,soils_data$long_wanted,drivers$lat[n],drivers$long[n],1))[1]
#     soil_info=c(pmax(1,soils_data$sand_top[i1]),pmax(1,soils_data$sand_bot[i1]),pmax(1,soils_data$clay_top[i1]),pmax(1,soils_data$clay_bot[i1]) )
#     if (soil_info[2] == 1) {soil_info[2] = soil_info[1]}
#     if (soil_info[4] == 1) {soil_info[4] = soil_info[3]}
     soil_info=c(mean(drivers$sand_top),mean(drivers$sand_bot),mean(drivers$clay_top),mean(drivers$clay_bot))
        if (is.loaded("racmgppet") == FALSE) { dyn.load("./acm_gpp_et.so") }
        tmp=.Fortran("racmgppet",output_dim=as.integer(output_dim),met=as.double(t(met)),pars=as.double(parameters)
                                ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                ,lat=as.double(lat),nopars=as.integer(nopars),nomet=as.integer(dim(met)[2])
                                ,nofluxes=as.integer(nofluxes),nopools=as.integer(nopools)
                                ,nodays=as.integer(dim(met)[1])
                                ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                ,soil_frac_clay=as.double(array(c(soil_info[3],soil_info[3],soil_info[4],soil_info[4]),dim=c(4)))
                                ,soil_frac_sand=as.double(array(c(soil_info[1],soil_info[1],soil_info[2],soil_info[2]),dim=c(4))) )
        output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
        if (n == dim(drivers)[1]) {dyn.unload("./acm_gpp_et.so")}
        rm(tmp) ; gc()

     # assign outputs to out final grids
     mean_lai[n] = mean(output[,,1]) # lai
     mean_gpp[n] = mean(output[,,2]) # GPP (gC.m-2.day-1)
     mean_transpiration[n] = mean(output[,,3])   # transpiration (kg.m-2.day-1)
     mean_wetcanopyevap[n] = mean(output[,,4])   # wetcanopy evaporation (kg.m-2.day-1)
     mean_soilevaporation[n] = mean(output[,,5]) # soil evaporation (kg.m-2.day-1)
     mean_wSWP[n] = mean(output[,,6])            # weighted soil water potential (MPa)
     mean_WUE[n] = mean_gpp[n]/mean_transpiration[n] # water use efficiency (gC/kgH2O)
     mean_rootwatermm[n] = mean(output[,,7])     # water in rooting zone (mm)
     mean_runoffmm[n] = mean(output[,,8])        # surface runoff (mm)
     mean_drainagemm[n] = mean(output[,,9])      # drainage from soil column (mm)

} # site loop

###
## Begin output to figures
###

par(mfrow=c(1,1), mar=c(5,5,3,1))
#plot(as.vector(mean_gpp[seq(2,length(mean_gpp),1)]-mean_gpp[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="LAI", ylab="GPP", cex=2, cex.lab=2, cex.main=2,cex.axis=2, col="blue")
#plot(as.vector(mean_WUE[seq(2,length(mean_gpp),1)]-mean_WUE[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="LAI", ylab="WUE", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
plot(as.vector(mean_transpiration[seq(2,length(mean_gpp),1)]-mean_transpiration[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="LAI", ylab="Transpiration", cex=2, cex.lab=2, cex.main=2,cex.axis=2)

#par(mfrow=c(2,4), mar=c(5,5,3,1))
#plot(mean_gpp ~ lai_iterations, type="l", lwd=3, xlab="LAI", ylab="GPP", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_WUE ~ lai_iterations, type="l", lwd=3, xlab="LAI", ylab="WUE", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_transpiration ~ lai_iterations, type="l", lwd=3, xlab="LAI", ylab="Transpiration", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_wetcanopyevap ~ lai_iterations, type="l", lwd=3, xlab="LAI", ylab="Wet Canopy", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_soilevaporation ~ lai_iterations, type="l", lwd=3, xlab="LAI", ylab="Soil Evaporation", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_rootwatermm ~ lai_iterations, type="l", lwd=3, xlab="LAI", ylab="SWC", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_runoffmm ~ lai_iterations, type="l", lwd=3, xlab="LAI", ylab="Runoff", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_drainagemm ~ lai_iterations, type="l", lwd=3, xlab="LAI", ylab="Drainage", cex=2, cex.lab=2, cex.main=2,cex.axis=2)

#par(mfrow=c(2,4), mar=c(5,5,3,1))
#plot(as.vector(mean_gpp[seq(2,length(mean_gpp),1)]-mean_gpp[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="LAI", ylab="GPP", cex=2, #cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_WUE[seq(2,length(mean_gpp),1)]-mean_WUE[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="LAI", ylab="WUE", cex=2, #cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_transpiration[seq(2,length(mean_gpp),1)]-mean_transpiration[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="LAI", #ylab="Transpiration", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_wetcanopyevap[seq(2,length(mean_gpp),1)]-mean_wetcanopyevap[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="LAI", #ylab="Wet Canopy", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_soilevaporation[seq(2,length(mean_gpp),1)]-mean_soilevaporation[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, #lab="LAI", ylab="Soil Evaporation", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_rootwatermm[seq(2,length(mean_gpp),1)]-mean_rootwatermm[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="LAI", #ylab="SWC", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_runoffmm[seq(2,length(mean_gpp),1)]-mean_runoffmm[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="LAI", #ylab="Runoff", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_drainagemm[seq(2,length(mean_gpp),1)]-mean_drainagemm[seq(1,length(mean_gpp)-1,1)]) ~ lai_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="LAI", #ylab="Drainage", cex=2, cex.lab=2, cex.main=2,cex.axis=2)

#par(mfrow=c(2,4), mar=c(5,5,3,1))
#plot(mean_gpp ~ avN_iterations, type="l", lwd=3, xlab="avN", ylab="GPP", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_WUE ~ avN_iterations, type="l", lwd=3, xlab="avN", ylab="WUE", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_transpiration ~ avN_iterations, type="l", lwd=3, xlab="avN", ylab="Transpiration", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_wetcanopyevap ~ avN_iterations, type="l", lwd=3, xlab="avN", ylab="Wet Canopy", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_soilevaporation ~ avN_iterations, type="l", lwd=3, xlab="avN", ylab="Soil Evaporation", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_rootwatermm ~ avN_iterations, type="l", lwd=3, xlab="avN", ylab="SWC", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_runoffmm ~ avN_iterations, type="l", lwd=3, xlab="avN", ylab="Runoff", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(mean_drainagemm ~ avN_iterations, type="l", lwd=3, xlab="avN", ylab="Drainage", cex=2, cex.lab=2, cex.main=2,cex.axis=2)

#par(mfrow=c(2,4), mar=c(5,5,3,1))
#plot(as.vector(mean_gpp[seq(2,length(mean_gpp),1)]-mean_gpp[seq(1,length(mean_gpp)-1,1)]) ~ avN_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="avN", ylab="GPP", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_WUE[seq(2,length(mean_gpp),1)]-mean_WUE[seq(1,length(mean_gpp)-1,1)]) ~ avN_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="avN", ylab="WUE", cex=2, #cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_transpiration[seq(2,length(mean_gpp),1)]-mean_transpiration[seq(1,length(mean_gpp)-1,1)]) ~ avN_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="avN", #ylab="Transpiration", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_wetcanopyevap[seq(2,length(mean_gpp),1)]-mean_wetcanopyevap[seq(1,length(mean_gpp)-1,1)]) ~ avN_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="avN", #ylab="Wet Canopy", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_soilevaporation[seq(2,length(mean_gpp),1)]-mean_soilevaporation[seq(1,length(mean_gpp)-1,1)]) ~ avN_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, #xlab="avN", ylab="Soil Evaporation", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_rootwatermm[seq(2,length(mean_gpp),1)]-mean_rootwatermm[seq(1,length(mean_gpp)-1,1)]) ~ avN_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="avN", ylab="SWC", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_runoffmm[seq(2,length(mean_gpp),1)]-mean_runoffmm[seq(1,length(mean_gpp)-1,1)]) ~ avN_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="avN", ylab="Runoff", cex=2, cex.lab=2, cex.main=2,cex.axis=2)
#plot(as.vector(mean_drainagemm[seq(2,length(mean_gpp),1)]-mean_drainagemm[seq(1,length(mean_gpp)-1,1)]) ~ avN_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="avN", ylab="Drainage", cex=2, cex.lab=2, cex.main=2,cex.axis=2)

#par(mfrow=c(1,1))
#plot(as.vector(mean_gpp[seq(2,length(mean_gpp),1)]-mean_gpp[seq(1,length(mean_gpp)-1,1)]) ~ avN_iterations[seq(2,length(mean_gpp),1)], type="l", lwd=3, xlab="avN", ylab="GPP", cex=2, cex.lab=2, cex.main=2,cex.axis=2, ylim=c(0,0.5))
#lines(as.vector(mean_gpp[seq(2,length(mean_gpp),1)]-mean_gpp[seq(1,length(mean_gpp)-1,1)]) ~ avN_iterations[seq(2,length(mean_gpp),1)], lwd=3)



