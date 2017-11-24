###
## Create needed ACM_GPP_ET shared object

# set to the working directory that this script should be called from
setwd("/home/lsmallma/WORK/R/Scripts/ACM_GPP_ET/") ; wkdir = getwd()
# compile the shared object containing ACM_GPP and ACM_ET
system("gfortran ./src/ACM_GPP_ET.f90 ./src/ACM_GPP_ET_R_interface.f90 -o ./src/acm_gpp_et.so -fPIC -shared")
system("mv ./src/acm_gpp_et.so .")

###
## Borrow met data from an existing CARDAMOM analysis

# set to the cardamom working directory for the moment
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM_R/")
## Load needed libraries and internal functions
source("./cardamom_functions/load_all_cardamom_functions.r")
# define file name for PROJECT file we will be borrowing from
PROJECTfile=paste("./CARDAMOM_OUTPUTS/DALECN_GSI_BUCKET_MHMCMC/global_1x1_new_acm/infofile.RData",sep="")
load(PROJECTfile) 
# this information will be used in loop to create the met inputs needed for the emulator, moving on the next task

## return back to working directory
setwd(wkdir)

###
## Define our output variables based on the grid of the CARDAMOM analysis we are borrowing

mean_lai = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_gpp = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_transpiration = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_wetcanopyevap = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_soilevaporation = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_rootwatermm = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_WUE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
mean_wSWP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

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

output_dim=7 ; nofluxes = 4 ; nopools = 1 ; nopars = 4 ; nos_iter = 1

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
     # ecosystem state drivers now rather than meteorology
     met[,11]= mean(drivers$obs[which(drivers$obs[,2] != -9999),2]) # LAI
     met[,12]= met[,11]*80 #50 #100  # root C stocks
     # parameters
     parameters = array(NA, dim=c(nopars,nos_iter))
     parameters[1,] = 1.89  # foliar N (gN.m-2)
     parameters[2,] = -9999 # min leaf water potential (MPa)
     parameters[3,] = 100   # root biomass needed to reach 50 % depth
     parameters[4,] = 1.5   # max root depth (m)

     # other inputs
     lat = drivers$lat
     soil_info=c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
     if (length(which(met[,11] > 0)) > 0) {
        if (is.loaded("racmgppet") == FALSE) { dyn.load("./acm_gpp_et.so") }
        tmp=.Fortran("racmgppet",output_dim=as.integer(output_dim),met=as.double(t(met)),pars=as.double(parameters)
                                ,out_var=as.double(array(0,dim=c(nos_iter,(dim(met)[1]),output_dim)))
                                ,lat=as.double(lat),nopars=as.integer(nopars),nomet=as.integer(dim(met)[2])
                                ,nofluxes=as.integer(nofluxes),nopools=as.integer(nopools)
                                ,nodays=as.integer(dim(met)[1])
                                ,deltat=as.double(array(0,dim=c(as.integer(dim(met)[1])))),nos_iter=as.integer(nos_iter)
                                ,soil_frac_clay=as.double(array(c(soil_info[3],soil_info[4]),dim=c(3)))
                                ,soil_frac_sand=as.double(array(c(soil_info[1],soil_info[2]),dim=c(3))) )
        output=tmp$out_var ; output=array(output, dim=c(nos_iter,(dim(met)[1]),output_dim))
        if (n == PROJECT$sites[length(PROJECT$sites)]) {dyn.unload("./acm_gpp_et.so")}
        rm(tmp) ; gc()
     }

     # assign outputs to out final grids
     mean_lai[slot_i,slot_j] = mean(output[,,1]) # lai
     mean_gpp[slot_i,slot_j] = mean(output[,,2]) # GPP (gC.m-2.day-1)
     mean_transpiration[slot_i,slot_j] = mean(output[,,3])   # transpiration (kg.m-2.day-1)
     mean_wetcanopyevap[slot_i,slot_j] = mean(output[,,4])   # wetcanopy evaporation (kg.m-2.day-1)
     mean_soilevaporation[slot_i,slot_j] = mean(output[,,5]) # soil evaporation (kg.m-2.day-1)
     mean_wSWP[slot_i,slot_j] = mean(output[,,6])            # weighted soil water potential (MPa)
     mean_WUE[slot_i,slot_j] = mean_gpp[slot_i,slot_j]/mean_transpiration[slot_i,slot_j] # water use efficiency (gC/kgH2O)
     mean_rootwatermm[slot_i,slot_j] = mean(output[,,7])     # water in rooting zone (mm)

} # site loop

fig_height=4000 ; fig_width=7200
jpeg(file="./FIGURES/figure_0.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(mean_gpp) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("GPP"," (gC ",m^-2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_transpiration) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Transpiration"," (kg ",m^-2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_wetcanopyevap) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Wet canopy evaporation"," (kg ",m^-2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_soilevaporation) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("Soil evaporation"," (kg ",m^-2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_WUE) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("WUE (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_wSWP) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("wSWP"," (MPa)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

jpeg(file="./FIGURES/figure_1.jpg", height=fig_height, width=fig_width, res=400, quality=100)
par(mfrow=c(2,3), mar=c(1.4, 1.2, 2.4, 6.5), omi=c(0.2, 0.2, 0.2, 0.40))
colour_choices_upper = colorRampPalette((brewer.pal(11,"Spectral")))
var1=(mean_gpp*365.25) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("GPP"," (gC ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_transpiration+mean_wetcanopyevap+mean_soilevaporation)*365.25 ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("ET"," (kg ",m^-2," y",r^-1,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=mean_gpp/(mean_transpiration+mean_wetcanopyevap+mean_soilevaporation) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("GPP/ET (gC/kgH2O)")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
var1=(mean_lai) ; colour_choices=colour_choices_upper(length(var1))
zaxis=c(min(var1,na.rm=TRUE),max(var1,na.rm=TRUE))
image.plot(var1, main=expression(paste("LAI"," (",m^2,"/",m^2,")")),zlim=zaxis, col=colour_choices,axes=FALSE, cex.main=2.6,legend.width=3.0,cex=1.8,axis.args=list(cex.axis=2.0,hadj=0.1))
dev.off()

print(paste("Global GPP = ",round(sum(mean_gpp*365.25*area*1e-15,na.rm=TRUE),digits=1)," PgC",sep=""))
