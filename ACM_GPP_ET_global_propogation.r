###
## Create needed ACM_GPP_ET shared object

# set to the working directory that this script should be called from
setwd("/home/lsmallma/WORK/GREENHOUSE/models/ACM_GPP_ET/") ; wkdir = getwd()
# compile the shared object containing ACM_GPP and ACM_ET
system("gfortran ./src/ACM_GPP_ET.f90 ./src/ACM_GPP_ET_R_interface.f90 -o ./src/acm_gpp_et.so -fPIC -shared")
system("mv ./src/acm_gpp_et.so .")

###
## Borrow met data from an existing CARDAMOM analysis

# set to the cardamom working directory for the moment
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")
## Load needed libraries and internal functions
source("./cardamom_functions/load_all_cardamom_functions.r")
# define file name for PROJECT file we will be borrowing from
PROJECTfile=paste("./CARDAMOM_OUTPUTS/DALECN_GSI_BUCKET_MHMCMC/global_1x1_new_acm/infofile.RData",sep="")
load(PROJECTfile)
# this information will be used in loop to create the met inputs needed for the emulator, moving on the next task
cardamom = nc_open("/disk/scratch/local.2/lsmallma/Forest2020/C_cycle_analyses/DALEC_GSI_DFOL_CWD_FR_1_2001_2015_NEE_GPP_Rh_Ra_Bio_lit_cwd_som_timeseries.nc")
lai = ncvar_get(cardamom,"lai_median") ; root = ncvar_get(cardamom,"root_median")
## return back to working directory
setwd(wkdir)

###
## Define our output variables based on the grid of the CARDAMOM analysis we are borrowing

mean_lai = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
sd_lai = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
sd_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_gpp = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_transpiration = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_wetcanopyevap = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_soilevaporation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_rootwatermm = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_wue = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_wSWP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
min_wSWP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_runoffmm = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_drainagemm = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

timeseries_lai = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
timeseries_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
timeseries_gpp = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
timeseries_transpiration = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
timeseries_wetcanopyevap = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
timeseries_soilevaporation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
timeseries_rootwatermm = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
timeseries_WUE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
timeseries_wSWP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
timeseries_runoffmm = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
timeseries_drainagemm = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))

# work out area matrix for the pixels in meters
# include adjustment for g-> Tg (*1e-12)
if (PROJECT$grid_type == "UK") {
    area=array(PROJECT$resolution**2, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
} else if (PROJECT$grid_type == "wgs84") {
    # generate the lat / long grid again
    output=generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
    # then generate the area estimates for each pixel
    area=calc_pixel_area(output$lat,output$long,PROJECT$resolution)
    # this output is in vector form and we need matching array shapes so...
    area=array(area, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
} else {
    stop("valid spatial grid option not selected (UK, or wgs84)")
}

###
## Some ACM_GPP_ET parameters

output_dim=9 ; nofluxes = 6 ; nopools = 1 ; nopars = 4 ; nos_iter = 1

# iterative process through the years...
for (n in seq(1,PROJECT$nosites)) {

     if (n%%1000 == 0){print(paste("...beginning site:",n," of ",PROJECT$nosites, sep=""))}

     # locational information
     slot_j=as.numeric(PROJECT$sites[n])/PROJECT$long_dim
     slot_i=as.numeric(PROJECT$sites[n])-(floor(slot_j)*PROJECT$long_dim)
     if(slot_i == 0) {slot_i = PROJECT$long_dim} ; slot_j=ceiling(slot_j)

     # load the met data for each site
     drivers=read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))

     # met note that the dimension here are different to that of drivers$met
     met=array(-9999,dim=c(length(PROJECT$model$timestep_days),12))
     met[,1] = drivers$met[,1]  # day of analysis
     met[,2] = drivers$met[,2]  # min temperature (oC)
     met[,3] = drivers$met[,3]  # max temperature (oC)
     met[,4] = drivers$met[,4]  # SW Radiation (MJ.m-2.day-1)
     met[,5] = drivers$met[,5]  # CO2 ppm
     met[,6] = drivers$met[,6]  # day of year
     met[,7] = drivers$met[,7]  # rainfall (kg.m-2.s-1)
     met[,8] = drivers$met[,14] # avg temperature (oC)
     met[,9] = drivers$met[,15] # avg wind speed (m.s-1)
     met[,10]= drivers$met[,16] # avg VPD (Pa)
     # Extract LAI (m2/m2) and root (gC/m2) from CARDAMOM analysis
     met[,11]= lai[n,]
     met[,12]= root[n,]

     # Assuming I have not LAI and root information we will run the analysis
     if (length(which(is.na(met[,11]) == TRUE)) == 0 & length(which(is.na(met[,12]) == TRUE)) == 0) {

	       # parameters
	       parameters = array(NA, dim=c(nopars,nos_iter))
	       parameters[1,] = 1.89  # foliar N (gN.m-2)
	       parameters[2,] = -9999 # min leaf water potential (MPa)
	       parameters[3,] = 100   # root biomass needed to reach 50 % depth
	       parameters[4,] = 2.0   # max root depth (m)

	       # other inputs
	       lat = drivers$lat
	       soil_info=c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
	       if (length(which(met[,11] > 0)) > 0) {
               # If the shared object has not been loaded yet do so...
	           if (is.loaded("racmgppet") == FALSE) { dyn.load("./acm_gpp_et.so") }
	           tmp=.Fortran("racmgppet",output_dim=as.integer(output_dim),met=as.double(t(met)),pars=as.double(parameters)
                                     ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                     ,lat=as.double(lat),nopars=as.integer(nopars),nomet=as.integer(dim(met)[2])
                                     ,nofluxes=as.integer(nofluxes),nopools=as.integer(nopools),nodays=as.integer(dim(met)[1])
                                     ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                     ,soil_frac_clay=as.double(array(c(soil_info[3],soil_info[4],soil_info[4]),dim=c(3)))
                                     ,soil_frac_sand=as.double(array(c(soil_info[1],soil_info[2],soil_info[2]),dim=c(3))) )
                   # extract output from the analysis
                   output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
                   # If this is the last site in the list best un-load the shared onject now
                   if (n == PROJECT$sites[length(PROJECT$sites)]) {dyn.unload("./acm_gpp_et.so")}
                   rm(tmp) ; gc()
	       } # If have LAI data

	       # assign time series to grid
	       timeseries_lai[slot_i,slot_j,] = (output[,1:dim(met)[1],1])             # lai (m2/m2)
               timeseries_root[slot_i,slot_j,] = root[n,]                              # root (gC/m2)
	       timeseries_gpp[slot_i,slot_j,] = (output[,1:dim(met)[1],2])             # GPP (gC.m-2.day-1)
	       timeseries_transpiration[slot_i,slot_j,] = (output[,1:dim(met)[1],3])   # transpiration (kg.m-2.day-1)
	       timeseries_wetcanopyevap[slot_i,slot_j,] = (output[,1:dim(met)[1],4])   # wetcanopy evaporation (kg.m-2.day-1)
	       timeseries_soilevaporation[slot_i,slot_j,] = (output[,1:dim(met)[1],5]) # soil evaporation (kg.m-2.day-1)
	       timeseries_wSWP[slot_i,slot_j,] = (output[,1:dim(met)[1],6])            # weighted soil water potential (MPa)
	       timeseries_rootwatermm[slot_i,slot_j,] = (output[,1:dim(met)[1],7])     # water in rooting zone (mm)
               timeseries_runoffmm[slot_i,slot_j,] = (output[,1:dim(met)[1],8])        # surface runoff (mm)
               timeseries_drainagemm[slot_i,slot_j,] = (output[,1:dim(met)[1],9])      # drainage / underflow from bottom of soil column (mm)

	       # assign timeseries mean values to grid
	       mean_lai[slot_i,slot_j] = mean(output[,1:dim(met)[1],1])             # lai (m2/m2)
               mean_root[slot_i,slot_j] = mean(root[n,])                            # root (gC/m2)
	       sd_lai[slot_i,slot_j] = sd(output[,1:dim(met)[1],1])                 # lai (m2/m2)
               sd_root[slot_i,slot_j] = sd(root[n,])                                # root (gC/m2)
	       mean_gpp[slot_i,slot_j] = mean(output[,1:dim(met)[1],2])             # GPP (gC.m-2.day-1)
	       mean_transpiration[slot_i,slot_j] = mean(output[,1:dim(met)[1],3])   # transpiration (kg.m-2.day-1)
	       mean_wetcanopyevap[slot_i,slot_j] = mean(output[,1:dim(met)[1],4])   # wetcanopy evaporation (kg.m-2.day-1)
	       mean_soilevaporation[slot_i,slot_j] = mean(output[,1:dim(met)[1],5]) # soil evaporation (kg.m-2.day-1)
	       mean_wSWP[slot_i,slot_j] = mean(output[,1:dim(met)[1],6])            # weighted soil water potential (MPa)
               min_wSWP[slot_i,slot_j] = min(output[,1:dim(met)[1],6])              # weighted soil water potential (MPa)
	       mean_rootwatermm[slot_i,slot_j] = mean(output[,1:dim(met)[1],7])     # water in rooting zone (mm)
               mean_runoffmm[slot_i,slot_j] = mean(output[,1:dim(met)[1],8])        # surface runoff (mm)
               mean_drainagemm[slot_i,slot_j] = mean(output[,1:dim(met)[1],9])      # drainage / underflow from bottom of soil column (mm)

     } # have got LAI and root infromation

} # site loop

###
## Generate some statistics
###

# Calculate the time series and mean values for water use efficiency (gC/kgH2O)
#mean_wue = mean_gpp/mean_transpiration
#timeseries_WUE = timeseries_gpp/timeseries_transpiration
mean_wue = mean_gpp/(mean_transpiration + mean_wetcanopyevap + mean_soilevaporation)
timeseries_wue = timeseries_gpp/(timeseries_transpiration + timeseries_wetcanopyevap + timeseries_soilevaporation)

# Global mean GPP (PgC/yr); note 1e-15 is conversion from gC to PgC
global_mean_annual_gpp = sum(mean_gpp*365.25*area*1e-15,na.rm=TRUE)
# Global mean Transpiration (PgH2O/yr); note 1e-12 is conversion from kgH2O to PgH2O
global_mean_annual_transpiration = sum(mean_transpiration*365.25*area*1e-12,na.rm=TRUE)
# Global mean Soil evaporation (PgH2O/yr); note 1e-12 is conversion from kgH2O to PgH2O
global_mean_annual_soilevaporation = sum(mean_soilevaporation*365.25*area*1e-12,na.rm=TRUE)
# Global mean Wet canopy evaporation (PgH2O/yr); note 1e-12 is conversion from kgH2O to PgH2O
global_mean_annual_wetcanopyevap = sum(mean_wetcanopyevap*365.25*area*1e-12,na.rm=TRUE)
# Global mean evapo-transpiration (PgH2O/yr); note 1e-12 is conversion from kgH2O to PgH2O
global_mean_annual_et = sum((mean_wetcanopyevap+mean_transpiration+mean_soilevaporation)*365.25*area*1e-12,na.rm=TRUE)
# Global mean water use efficiency (gC/kgH2O)
global_mean_annual_wue = mean(mean_wue,na.rm=TRUE)
# Global mean weighted soil water potential (MPa)
global_mean_annual_wSWP = mean(mean_wSWP,na.rm=TRUE)
# Global mean water in rooted zone (kgH2O/m2)
global_mean_annual_rootwatermm = mean(mean_rootwatermm,na.rm=TRUE)
# Global mean LAI (m2/m2)
global_mean_annual_lai = mean(mean_lai,na.rm=TRUE)
# Global mean LAI (gC/m2)
global_mean_annual_root = mean(mean_root,na.rm=TRUE)

###
## Generate figures
###

jpeg(file="./FIGURES/Cal_val_paper_figure_3.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,2), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Mean status of biophysical inputs
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(mean_lai) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("LAI"," (",m^2,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_root) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Fine Root"," (gC/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# Standard deviation of biophysical inputs
var1=(sd_lai) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("LAI"," (",m^2,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(sd_root) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Fine Root"," (gC/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/Cal_val_paper_figure_4.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(3,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
# Summary figures GPP, ET, WUE
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_transpiration+mean_wetcanopyevap+mean_soilevaporation)*365.25 ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("ET"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=mean_wue ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("WUE (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# ET component figures transpiration, wetcanopy evaporation, soil evaporation
var1=(mean_transpiration*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Transpiration"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_wetcanopyevap*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Wet canopy evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_soilevaporation*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Soil evaporation"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
# Soil runoff / drainage and status
var1=(mean_runoffmm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Soil run-off"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_drainagemm*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Soil drainage"," (kgH2O ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_rootwatermm) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Water in rooting zone"," (",kg,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

###
## Save output to files for later use
###

units=c("LAI = m2/m2","Roots = gC/m2","Water use efficiency (WUE) = gC/kgH2O"
       ,"GPP = gC/m2/day, global_mean_annual_gpp = PgC"
       ,"All water fluxes = kgH2O/m2/day except global_mean* = PgH2O"
       ,"mean_rootwatermm = kg/m2","All soil water potentials (SWP) = MPa")
# Save output for later use
global_output = list(        units = units,
            global_mean_annual_gpp = global_mean_annual_gpp,
  global_mean_annual_transpiration = global_mean_annual_transpiration,
global_mean_annual_soilevaporation = global_mean_annual_soilevaporation,
  global_mean_annual_wetcanopyevap = global_mean_annual_wetcanopyevap,
             global_mean_annual_et = global_mean_annual_et,
            global_mean_annual_wue = global_mean_annual_wue,
           global_mean_annual_wSWP = global_mean_annual_wSWP,
    global_mean_annual_rootwatermm = global_mean_annual_rootwatermm,
            global_mean_annual_lai = global_mean_annual_lai,
           global_mean_annual_root = global_mean_annual_root,
                          mean_lai = mean_lai,
                         mean_root = mean_root,
                            sd_lai = sd_lai,
                           sd_root = sd_root,
                          mean_gpp = mean_gpp,
                mean_transpiration = mean_transpiration,
                mean_wetcanopyevap = mean_wetcanopyevap,
              mean_soilevaporation = mean_soilevaporation,
                  mean_rootwatermm = mean_rootwatermm,
                          mean_wue = mean_wue,
                         mean_wSWP = mean_wSWP,
                          min_wSWP = min_wSWP,
                     mean_runoffmm = mean_runoffmm,
                   mean_drainagemm = mean_drainagemm,
                    timeseries_lai = timeseries_lai,
                   timeseries_root = timeseries_root,
                    timeseries_gpp = timeseries_gpp,
          timeseries_transpiration = timeseries_transpiration,
          timeseries_wetcanopyevap = timeseries_wetcanopyevap,
        timeseries_soilevaporation = timeseries_soilevaporation,
            timeseries_rootwatermm = timeseries_rootwatermm,
                    timeseries_WUE = timeseries_WUE,
                   timeseries_wSWP = timeseries_wSWP,
               timeseries_runoffmm = timeseries_runoffmm,
             timeseries_drainagemm = timeseries_drainagemm)

# Now save the file
save(global_output, file="./outputs/global_1x1_degree_2001_2015.RData")

###
## Print some default information to the user
###
print(paste("Global GPP = ",round(global_mean_annual_gpp,digits=1)," PgC",sep=""))
print(paste("Global Transpiration = ",round(global_mean_annual_transpiration,digits=1)," PgH2O",sep=""))
print(paste("Global Soil Evaporation = ",round(global_mean_annual_soilevaporation,digits=1)," PgH2O",sep=""))
print(paste("Global Wet Canopy Evaporation = ",round(global_mean_annual_wetcanopyevap,digits=1)," PgH2O",sep=""))
print(paste("Global ET = ",round(global_mean_annual_et,digits=1)," PgH2O",sep=""))
print(paste("Global WUE = ",round(mean_wue,digits=2)," gC/kgH2O",sep=""))
print(paste("Global wSWP = ",round(global_mean_annual_wSWP,digits=1)," MPa",sep=""))
print(paste("Global Water in root zone = ",round(global_mean_annual_rootwatermm,digits=1),sep=""))

###
## Save output to file for later analysis
###
