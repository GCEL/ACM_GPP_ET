
module CARBON_MODEL_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL        &
         ,opt_max_scaling     &
         ,dble_one,dble_zero  &
         ,wSWP_time           &
         ,nos_soil_layers     &
         ,soil_frac_clay      &
         ,soil_frac_sand

!!!!!!!!!
! Parameters
!!!!!!!!!

! useful technical parameters
logical :: do_iWUE = .true., & ! Use iWUE or WUE for stomatal optimisation
 do_energy_balance = .false.   ! Calculate steady-state energy balance for GPP~Transpiration
double precision, parameter :: dble_zero = 0d0    &
                              ,dble_one = 1d0     &
                              ,vsmall = tiny(0d0)

integer, parameter :: nos_root_layers = 3, nos_soil_layers = nos_root_layers + 1
double precision, parameter :: pi = 3.1415927d0,  &
                             pi_1 = pi**(-dble_one), &
                              pi2 = pi**2,        &
                           two_pi = pi*2d0,       &
                       deg_to_rad = pi/180d0,     &
              sin_dayl_deg_to_rad = sin( 23.45d0 * deg_to_rad ), & ! repeated function in acm
                          gravity = 9.8067d0,     & ! acceleration due to gravity, ms-1
                            boltz = 5.670400d-8,  & ! Boltzmann constant (W.m-2.K-4)
                       emissivity = 0.96d0,       &
                      emiss_boltz = emissivity * boltz, &
                  sw_par_fraction = 0.5d0,        & ! fraction of short-wave radiation which is PAR
                           freeze = 273.15d0,     &
                       gs_H2O_CO2 = 1.646259d0,   & ! The ratio of H20:CO2 diffusion for gs (Jones appendix 2)
                     gs_H2O_CO2_1 = gs_H2O_CO2 ** (-dble_one), &
                       gb_H2O_CO2 = 1.37d0,       & ! The ratio of H20:CO2 diffusion for gb (Jones appendix 2)
          partial_molar_vol_water = 18.05d-6,     & ! partial molar volume of water, m3 mol-1 at 20C
                       umol_to_gC = 1d-6*12d0,    & ! conversion of umolC -> gC
                 mmol_to_kg_water = 1.8d-5,       & ! milli mole conversion to kg
                   mol_to_g_water = 18d0,         & ! molecular mass of water (g)
                     mol_to_g_co2 = 12d0,         & ! molecular mass of CO2 (g)
                     g_to_mol_co2 = 1d0/12d0,     &
!snowscheme       density_of_water = 998.9d0,         & ! density of !water kg.m-3
                   gas_constant_d = 287.04d0,     & ! gas constant for dry air (J.K-1.mol-1)
                             Rcon = 8.3144d0,     & ! Universal gas constant (J.K-1.mol-1)
                        vonkarman = 0.41d0,       & ! von Karman's constant
                      vonkarman_2 = vonkarman**2, & ! von Karman's constant^2
                            cpair = 1004.6d0        ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

! photosynthesis / respiration parameters
double precision, parameter :: &
                    kc_saturation = 310d0,        & ! CO2 half saturation, at reference temperature (298.15 K)
                 kc_half_sat_conc = 23.956d0,     & ! CO2 half sat, sensitivity coefficient
               co2comp_saturation = 36.5d0,       & ! CO2 compensation point, at reference temperature (298.15 K)
            co2comp_half_sat_conc = 9.46d0          ! CO2 comp point, sensitivity coefficient

! hydraulic parameters
double precision, parameter :: &
                       tortuosity = 2.5d0,        & ! tortuosity
                           gplant = 5d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                      root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                      root_radius = 0.00029d0,    & ! root radius (m) Bonen et al 2014 = 0.00029
                                                    ! Williams et al 1996 = 0.0001
                    root_radius_1 = root_radius**(-dble_one), &
              root_cross_sec_area = pi * root_radius**2, & ! root cross sectional area (m2)
                                                           ! = pi * root_radius * root_radius
                     root_density = 0.31d6,       & ! root density (g biomass m-3 root)
                                                    ! 0.5e6 Williams et al 1996
                                                    ! 0.31e6 Bonan et al 2014
          root_mass_length_coef_1 = (root_cross_sec_area * root_density)**(-dble_one), &
               const_sfc_pressure = 101325d0,     & ! (Pa)  Atmospheric surface pressure
                             head = 0.009807d0,   & ! head of pressure (MPa/m)
                           head_1 = 101.968d0       ! inverse head of pressure (m/MPa)

! structural parameters
double precision, parameter :: &
            canopy_height_default = 9d0,          & ! canopy height assumed to be 9 m
             tower_height_default = canopy_height_default + 2d0, & ! tower (observation) height assumed to be 2 m above canopy
                         min_wind = 0.1d0,        & ! minimum wind speed at canopy top
                     min_drythick = 0.01d0,       & ! minimum dry thickness depth (m)
                        min_layer = 0.03d0,       & ! minimum thickness of the second rooting layer (m)
                      soil_roughl = 0.05d0,       & ! soil roughness length (m)
                   top_soil_depth = 0.1d0,        & ! thickness of the top soil layer (m)
                   mid_soil_depth = 0.2d0,        & ! thickness of the second soil layer (m)
                          min_lai = 1.5d0,        & ! minimum LAI assumed for aerodynamic conductance calculations (m2/m2)
                         min_root = 5d0,          & ! minimum root biomass (gBiomass.m-2)
                  min_throughfall = 0.2d0,        & ! minimum fraction of precipitation which
                                                    ! is through fall
                      min_storage = 0.2d0           ! minimum canopy water (surface) storage (mm)

! timing parameters
double precision, parameter :: &
                 seconds_per_hour = 3600d0,       & ! Number of seconds per hour
                  seconds_per_day = 86400d0,      & ! number of seconds per day
                seconds_per_day_1 = 1d0/seconds_per_day ! inverse of seconds per day

!!!!!!!!!
! Module level variables
!!!!!!!!!

! hydraulic model variables
integer :: water_retention_pass, soil_layer, sunrise, sunset
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand ! clay and sand fractions of soil
double precision, dimension(nos_root_layers) :: uptake_fraction, & !
                                                     water_flux, & ! potential transpiration flux (mmol.m-2.s-1)
                                                         demand    ! maximum potential canopy hydraulic demand
double precision, dimension(nos_soil_layers+1) :: SWP, & ! soil water potential (MPa)
                                    soil_conductivity, & ! soil conductivity
                                            waterloss, & ! water loss from specific soil layers (m)
                                            watergain, & ! water gained by specfic soil layers (m)
                                       field_capacity, & ! soil field capacity (m3.m-3)
                                       soil_waterfrac, & ! soil water content (m3.m-3)
                                             porosity, & ! soil layer porosity, (fraction)
                                      layer_thickness, & ! thickness of soil layers (m)
                      cond1, cond2, cond3, potA, potB    ! Saxton equation values

double precision :: root_reach, root_biomass, &
                  drythick, & ! estimate of the thickness of the dry layer at soil surface (m)
                      wSWP, & ! weighted soil water potential (MPa) used in GSI calculate.
                              ! Removes / limits the fact that very low root density in young plants
                              ! give values too large for GSI to handle.
              tower_height, & ! height of drivers above the canopy (m)
             canopy_height, & ! height of canopy (m)
                 max_depth, & ! maximum possible root depth (m)
                    root_k, & ! biomass to reach half max_depth
                    runoff, & ! runoff (kgH2O.m-2.day-1)
                 underflow, & ! drainage from the bottom of soil column (kgH2O.m-2.day-1)
  new_depth,previous_depth, & ! depth of bottom of soil profile
               canopy_wind, & ! wind speed (m.s-1) at canopy top
                     ustar, & ! friction velocity (m.s-1)
                  ustar_Uh, &
            air_density_kg, & ! air density kg/m3
                    roughl, & ! roughness length (m)
              displacement, & ! zero plane displacement (m)
                max_supply, & ! maximum water supply (mmolH2O/m2/day)
                     meant, & ! mean air temperature (oC)
                   meant_K, & ! mean air temperature (K)
                 maxt_lag1, &
                     leafT, & ! canopy temperature (oC)
          mean_annual_temp, &
        canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
          soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
            soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
             sky_lwrad_Wm2, & ! sky absorbed longwave radiation (W.m-2)
                 ci_global, & ! internal CO2 concentration (ppm or umol/mol)
      stomatal_conductance, & ! maximum stomatal conductance (mmolH2O.m-2.s-1)
   aerodynamic_conductance, & ! bulk surface layer conductance (m.s-1)
          soil_conductance, & ! soil surface conductance (m.s-1)
         convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
       air_vapour_pressure, & ! Vapour pressure of the air (kPa)
                    lambda, & ! latent heat of vapourisa/tion (J.kg-1)
                     psych, & ! psychrometric constant (kPa K-1)
                     slope, & ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
              snow_storage, & ! snow storage (kgH2O/m2)
            canopy_storage, & ! water storage on canopy (kgH2O.m-2)
      intercepted_rainfall    ! intercepted rainfall rate equivalent (kg.m-2.s-1)

! Module level variables for ACM_GPP_ET parameters
double precision :: delta_gs, & ! day length corrected gs increment mmolH2O/m2/dayl
                       avN, & ! average foliar N (gN/m2)
                      iWUE, & ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
                       NUE, & ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                              ! ,unlimited by CO2, light and photoperiod (gC/gN/m2leaf/day)
               pn_max_temp, & ! Maximum temperature for photosynthesis (oC)
               pn_opt_temp, & ! Optimum temperature fpr photosynthesis (oC)
               pn_kurtosis, & ! Kurtosis of photosynthesis temperature response
                        e0, & ! Quantum yield gC/MJ/m2/day PAR
              co2_half_sat, & ! CO2 at which photosynthesis is 50 % of maximum (ppm)
            co2_comp_point, & ! CO2 at which photosynthesis > 0 (ppm)
                    minlwp, & ! min leaf water potential (MPa)
 max_lai_lwrad_transmitted, & ! Max fraction of LW from sky transmitted by canopy
lai_half_lwrad_transmitted, & ! LAI at which canopy LW transmittance = 50 %
    max_lai_nir_reflection, & ! Max fraction of NIR reflected by canopy
   lai_half_nir_reflection, & ! LAI at which canopy NIR refection = 50 %
    max_lai_par_reflection, & ! Max fraction of PAR refected by canopy
   lai_half_par_reflection, & ! LAI at which canopy PAR reflection = 50 %
   max_lai_par_transmitted, & ! minimum transmittance = 1-par
  lai_half_par_transmitted, & ! LAI at which 50 %
   max_lai_nir_transmitted, & ! minimum transmittance = 1-par
  lai_half_nir_transmitted, & ! LAI at which 50 %
   max_lai_lwrad_reflected, & !
  lai_half_lwrad_reflected, & ! LAI at which 50 % LW is reflected back to sky
     soil_swrad_absorption, & ! Fraction of SW rad absorbed by soil
     max_lai_lwrad_release, & ! 1-Max fraction of LW emitted from canopy to be
    lai_half_lwrad_release    ! LAI at which LW emitted from canopy to be released at 50 %

! Module level variables for step specific met drivers
double precision :: mint, & ! minimum temperature (oC)
                    maxt, & ! maximum temperature (oC)
      airt_zero_fraction, & ! fraction of air temperature above freezing
                   swrad, & ! incoming short wave radiation (MJ/m2/day)
                     co2, & ! CO2 (ppm)
                     doy, & ! Day of year
                rainfall, & ! rainfall (kgH2O/m2/s)
                snowfall, &
               snow_melt, & ! snow melt (kgH2O/m2/s)
                wind_spd, & ! wind speed (m/s)
                  vpd_pa, & ! Vapour pressure deficit (Pa)
                     lai, & ! leaf area index (m2/m2)
                  root_C    ! fine root stock (gC/m2)

! Module level varoables for step specific timing information
double precision :: cos_solar_zenith_angle, &
                    seconds_per_step, & !
                       days_per_step, & !
                     days_per_step_1, & !
                        dayl_seconds, & ! day length in seconds
                      dayl_seconds_1, &
                          dayl_hours    ! day length in hours

double precision, dimension(:), allocatable ::    deltat_1, & ! inverse of decimal days
                                               soilwatermm, &
                                                 wSWP_time

contains
!
!--------------------------------------------------------------------
!
  subroutine CARBON_MODEL(start,finish,met,pars,deltat,nodays,lat,FLUXES,POOLS &
                         ,nopars,nomet,nopools,nofluxes)

    ! The Carbon model, seen here, is derived from that used within the CARDAMOM
    ! framework to run the DALEC model. Here, all C cycle components have been
    ! removed except the aggregated canopy model (ACM_GPP) for Gross Primary Productivity (GPP),
    ! the aggregated canopy model (ACM_ET) for evapotranspiration and the 3 soil layer bucket model.

    implicit none

    ! declare input variables
    integer, intent(in) :: start    &
                          ,finish   &
                          ,nopars   & ! number of paremeters in vector
                          ,nomet    & ! number of meteorological fields
                          ,nofluxes & ! number of model fluxes
                          ,nopools  & ! number of model pools
                          ,nodays     ! number of days in simulation

    double precision, intent(in) :: met(nomet,nodays) & ! met drivers
                         ,deltat(nodays)    & ! time step in decimal days
                         ,pars(nopars)      & ! number of parameters
                         ,lat                 ! site latitude (degrees)

    double precision, dimension((nodays+1),nopools), intent(inout) :: POOLS ! vector of ecosystem pools

    double precision, dimension(nodays,nofluxes), intent(inout) :: FLUXES ! vector of ecosystem fluxes

    ! declare general local variables
    double precision ::  infi &
                         ,tmp &
               ,act_pot_ratio &
               ,transpiration &
             ,soilevaporation &
              ,wetcanopy_evap &
            ,snow_sublimation &
               ,ET_pot,ET_net &
                     ,deltaWP & ! deltaWP (MPa) minlwp-soilWP
                  ,isothermal &
            ,deltaR,deltaTemp &
                        ,Rtot   ! MPa.s-1.m-2.mmol-1

    integer :: nxp,n

    ! met drivers are:
    ! 1st run day
    ! 2nd min daily temp (oC)
    ! 3rd max daily temp (oC)
    ! 4th Radiation (MJ.m-2.day-1)
    ! 5th CO2 (ppm)
    ! 6th DOY
    ! 7th precipitation (kg.m-2.s-1)
    ! 8th avg daily temperature (oC)
    ! 9th avg daily wind speed (m.s-1)
    ! 10th avg daily VPD (Pa)
    ! 11th LAI
    ! 12th root C

    ! POOLS are:
    ! 1 = water in rooting zone (mm)

    ! FLUXES are:
    ! 1 = GPP (gC.m-2.day-1)
    ! 2 = transpiration (kg.m-2.day-1)
    ! 3 = wetcanopy evaporation (kg.m-2.day-1)
    ! 4 = soil evaporation (kg.m-2.day-1)

    ! PARAMETERS
    ! 3 process parameters;

    ! p(1) foliar nitrogen concentration (gN.m-2)
    ! p(2) min leaf water potential (MPa)
    ! p(3) root biomass to 50 % depth (g.m-2)
    ! p(4) maxmimum rooting depth (m)

! profiling example
!real :: begin, done,f1=0,f2=0,f3=0,f4=0,f5=0,total_time = 0
!real :: Rtot_track_time=0, aero_time=0 , soilwater_time=0 , acm_et_time = 0 , Rm_time = 0
!call cpu_time(begin)
!call cpu_time(done)

    ! infinity check requirement
    infi = 0d0

    ! load ACM-GPP-ET parameters
    NUE                        = 1.182549d+01  ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                               ! ,unlimited by CO2, light and
                                               ! photoperiod (gC/gN/m2leaf/day)
    pn_max_temp                = 5.357174d+01  ! Maximum temperature for photosynthesis (oC)
    pn_opt_temp                = 3.137242d+01  ! Optimum temperature for photosynthesis (oC)
    pn_kurtosis                = 1.927458d-01  ! Kurtosis of photosynthesis temperature response
    e0                         = 5.875662d+00  ! Quantum yield gC/MJ/m2/day PAR
    max_lai_lwrad_transmitted  = 7.626683d-01  ! Max fractional reduction of LW from sky transmitted through canopy
    lai_half_lwrad_transmitted = 7.160363d-01  ! LAI at which canopy LW transmittance reduction = 50 %
    max_lai_nir_reflection     = 4.634860d-01  ! Max fraction of NIR reflected by canopy
    lai_half_nir_reflection    = 1.559148d+00  ! LAI at which canopy NIR reflected = 50 %
    minlwp                     =-1.996830d+00  ! minimum leaf water potential (MPa)
    max_lai_par_reflection     = 1.623013d-01  ! Max fraction of PAR reflected by canopy
    lai_half_par_reflection    = 1.114360d+00  ! LAI at which canopy PAR reflected = 50 %
    lai_half_lwrad_reflected   = 1.126214d+00  ! LAI at which 50 % LW is reflected back to sky
    iWUE                       = 1.602503d-06  ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
    soil_swrad_absorption      = 6.643079d-01  ! Fraction of SW rad absorbed by soil
    max_lai_par_transmitted    = 8.079519d-01  ! Max fractional reduction in PAR transmittance by canopy
    lai_half_par_transmitted   = 9.178784d-01  ! LAI at which PAR transmittance reduction = 50 %
    max_lai_nir_transmitted    = 8.289803d-01  ! Max fractional reduction in NIR transmittance by canopy
    lai_half_nir_transmitted   = 1.961831d+00  ! LAI at which NIR transmittance reduction = 50 %
    max_lai_lwrad_release      = 9.852855d-01  ! Max fraction of LW emitted (1-par) from canopy to be released
    lai_half_lwrad_release     = 7.535450d-01  ! LAI at which LW emitted from canopy to be released at 50 %
    max_lai_lwrad_reflected    = 1.955832d-02  ! LAI at which 50 % LW is reflected back to sky

    ! load ACM-GPP-ET parameters
!    NUE                       = pars(5+1)     ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
!                                              ! ,unlimited by CO2, light and
!                                              ! photoperiod
!                                              ! (gC/gN/m2leaf/day)
!    pn_max_temp               = pars(5+2)  ! Maximum temperature for photosynthesis (oC)
!    pn_opt_temp               = pars(5+3)  ! Optimum temperature for photosynthesis (oC)
!    pn_kurtosis               = pars(5+4)  ! Kurtosis of photosynthesis temperature response
!    e0                        = pars(5+5)  ! Quantum yield gC/MJ/m2/day PAR
!    max_lai_lwrad_absorption  = pars(5+6)  ! Max fraction of LW from sky absorbed by canopy
!    lai_half_lwrad_absorption = pars(5+7)  ! LAI at which canopy LW absorption = 50 %
!    max_lai_nir_absorption    = pars(5+8)  ! Max fraction of NIR absorbed by canopy
!    lai_half_nir_absorption   = pars(5+9)  ! LAI at which canopy NIR absorption = 50 %
!    minlwp                    = pars(5+10) ! minimum leaf water potential (MPa)
!    max_lai_par_absorption    = pars(5+11)  ! Max fraction of PAR absorbed by canopy
!    lai_half_par_absorption   = pars(5+12)  ! LAI at which canopy PAR absorption = 50 %
!    lai_half_lwrad_to_sky     = pars(5+13)  ! LAI at which 50 % LW is reflected back to sky
!    iWUE                      = pars(5+14)  ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
!    soil_swrad_absorption     = pars(5+15)  ! Fraction of SW rad absorbed by soil
!    max_lai_swrad_reflected   = pars(5+16)  ! Max fraction of SW reflected back to sky
!    lai_half_swrad_reflected  = (lai_half_nir_absorption+lai_half_par_absorption) * 0.5d0
!    max_lai_lwrad_release     = pars(5+17)  ! Max fraction of LW emitted from canopy to be released
!    lai_half_lwrad_release    = pars(5+18)  ! LAI at which LW emitted from canopy to be released at 50 %

    ! check loaded parameters
    avN = pars(1)
    ! if input minlwp not empty over-write default value
    if (pars(2) /= -9999d0) minlwp = pars(2)
    deltaWP = minlwp
    ! load root distribution information
    root_k = pars(3) ; max_depth = pars(4)
    ! load canopy information if available
    canopy_height = canopy_height_default
    if (size(pars) > 4 .and. pars(5) > 0d0) canopy_height = pars(5) ! if some cases use input canopy height
    tower_height = canopy_height + 2d0

    ! reset values
    intercepted_rainfall = dble_zero ; canopy_storage = dble_zero ; snow_storage = dble_zero

    ! SHOULD TURN THIS INTO A SUBROUTINE CALL AS COMMON TO BOTH DEFAULT AND CROPS
    if (.not.allocated(deltat_1)) then
       allocate(deltat_1(nodays),wSWP_time(nodays),soilwatermm(nodays))
       deltat_1 = deltat**(-dble_one)
    endif

    ! zero variables not done elsewhere
    water_flux = dble_zero
    ! initialise some time invarient parameters
    call saxton_parameters(soil_frac_clay,soil_frac_sand)
    call initialise_soils(soil_frac_clay,soil_frac_sand)

    ! determine the number of seconds per time step being used
    seconds_per_step = seconds_per_day * deltat(1)
    days_per_step = deltat(1)
    mint = met(2,1)  ! minimum temperature (oC)
    maxt = met(3,1)  ! maximum temperature (oC)
    maxt_lag1 = maxt
    leafT = maxt
    if (met(8,1) /= -9999d0) then
        meant = met(8,1)  ! mean air temperature (oC)
!        mean_annual_temp = sum(met(8,1:nodays)) / dble(nodays)
    else
        meant = (met(3,1)+met(2,1))*0.5d0
!        mean_annual_temp = sum((met(3,1:nodays)+met(2,1:nodays))*0.5d0) / dble(nodays)
    endif
    meant_K = meant + freeze
    lai = met(11,1) ! leaf area index (m2/m2)
    root_C = met(12,1) ! fine root stock (gC/m2)

    ! initialise root reach based on initial conditions
    root_biomass = max(min_root,root_C*2d0)
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! determine initial soil layer thickness
    layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
    layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
    layer_thickness(4) = max_depth - sum(layer_thickness(1:3))
    layer_thickness(5) = top_soil_depth

    previous_depth = max(top_soil_depth,root_reach)
    ! needed to initialise soils
    call calculate_Rtot(Rtot)
    ! used to initialise soils
    call calculate_update_soil_water(dble_zero,dble_zero,dble_zero,FLUXES(1,2)) ! assume no evap or rainfall
    ! store soil water content of the rooting zone (mm)
    POOLS(1,1) = 1d3*soil_waterfrac(1)*layer_thickness(1)

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! set lag information using previous time step value for temperature
      maxt_lag1 = maxt

      ! Incoming drivers
      mint = met(2,n)  ! minimum temperature (oC)
      maxt = met(3,n)  ! maximum temperature (oC)
      leafT = maxt     ! initial canopy temperature (oC)
      swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
      co2 = met(5,n)   ! CO2 (ppm)
      doy = met(6,n)   ! Day of year
      rainfall = max(dble_zero,met(7,n)) ! rainfall (kgH2O/m2/s)
      ! calculate mean air temperature dependending on available information
      if (met(8,n) /= -9999d0) then
          meant = met(8,n)  ! mean air temperature (oC)
      else
          meant = (met(3,n)+met(2,n))*0.5d0
      endif
      meant_K = meant + freeze ! oC -> K
      airt_zero_fraction = (maxt-dble_zero) / (maxt-mint) ! fraction of temperture period above freezing
      wind_spd = met(9,n) ! wind speed (m/s)
      vpd_pa = met(10,n)  ! Vapour pressure deficit (Pa)
      lai = met(11,n)     ! leaf area index (m2/m2)
      root_C = met(12,n)  ! fine root stock (gC/m2)

      ! calculate daylength in hours and seconds
      call calculate_daylength((doy-(deltat(n)*0.5d0)),lat)
      ! plus various time related variables needed thoughout
      dayl_seconds_1 = dayl_seconds ** (-dble_one)
      seconds_per_step = seconds_per_day * deltat(n)
      days_per_step = deltat(n)
      days_per_step_1 = deltat_1(n)

      ! snowing or not...?
      snow_melt = dble_zero ; snowfall = dble_zero
      if (mint < dble_zero .and. maxt > dble_zero) then
          ! if minimum temperature is below freezing point then we weight the
          ! rainfall into snow or rain based on proportion of temperature below
          ! freezing
          snowfall = dble_one - airt_zero_fraction
          snowfall = rainfall * snowfall ; rainfall = rainfall - snowfall
          ! Add rainfall to the snowpack and clear rainfall variable
          snow_storage = snow_storage + (snowfall*seconds_per_step)
          ! Also melt some of the snow
          snow_melt = airt_zero_fraction
          ! otherwise we assume snow is melting at 10 % per day light hour
          snow_melt = min(snow_storage, snow_melt * snow_storage * dayl_hours * 0.1d0 * deltat(n))
          snow_storage = snow_storage - snow_melt
      elseif (maxt < dble_zero) then
          ! if whole day is below freezing then we should assume that all
          ! precipitation is snowfall
          snowfall = rainfall ; rainfall = dble_zero
          ! Add rainfall to the snowpack and clear rainfall variable
          snow_storage = snow_storage + (snowfall*seconds_per_step)
      else if (mint > dble_zero) then
          ! otherwise we assume snow is melting at 10 % per day light hour
          snow_melt = min(snow_storage, snow_storage * dayl_hours * 0.1d0 * deltat(n))
          snow_storage = snow_storage - snow_melt
      end if

      !!!!!!!!!!
      ! calculate soil water potential and total hydraulic resistance
      !!!!!!!!!!

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      root_biomass = max(min_root,root_C*2d0)
      ! estimate drythick for the current step
      drythick = max(min_drythick, top_soil_depth * min(dble_one,dble_one - (soil_waterfrac(1) / porosity(1))))
      call calculate_Rtot(Rtot)
      ! Pass wSWP to output variable and update deltaWP between minlwp and
      ! current weighted soil WP
      wSWP_time(n) = wSWP ; deltaWP = min(dble_zero,minlwp-wSWP)

      !!!!!!!!!!
      ! Calculate surface exchange coefficients
      !!!!!!!!!!

      ! calculate some temperature dependent meteorologial properties
      call acm_meteorological_constants(maxt)
      ! calculate aerodynamic using consistent approach with SPA
      call calculate_aerodynamic_conductance

      !!!!!!!!!!
      ! Determine net shortwave and isothermal longwave energy balance
      !!!!!!!!!!

      call calculate_shortwave_balance
      call calculate_longwave_isothermal(leafT,maxt)

      !!!!!!!!!!
      ! Estimate approximate wet canopy evaporation and impact on energy balance
      !!!!!!!!!!

      ! if desired calculate the steady-state energy balance
      if (do_energy_balance) then
          isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * seconds_per_day_1)
          call update_net_radiation(isothermal,leafT,lai,dble_one &
                                   ,dble_zero,aerodynamic_conductance,vpd_pa &
                                   ,deltaTemp,deltaR)
          ! update long wave and canopy temperature based on potential canopy surface flux
          canopy_lwrad_Wm2 = canopy_lwrad_Wm2 + deltaR
          leafT = leafT + deltaTemp
          ! Canopy intercepted rainfall evaporation (kgH2O/m2/day)
          call calculate_wetcanopy_evaporation(wetcanopy_evap,act_pot_ratio,canopy_storage,dble_zero)
          ! restore temperature and radiation values
          leafT = leafT - deltaTemp ; canopy_lwrad_Wm2 = canopy_lwrad_Wm2 - deltaR
      else 
          ! Canopy intercepted rainfall evaporation (kgH2O/m2/day)
          call calculate_wetcanopy_evaporation(wetcanopy_evap,act_pot_ratio,canopy_storage,dble_zero)
      endif ! do energy balance

      ! calculate radiation absorption and estimate stomatal conductance
      call acm_albedo_gc(abs(deltaWP),Rtot)

      !!!!!!!!!!
      ! Evaptranspiration (kgH2O.m-2.day-1)
      !!!!!!!!!!

      ! Canopy transpiration (kgH2O/m2/day)
      call calculate_transpiration(transpiration)
      ! Soil surface (kgH2O.m-2.day-1)
      call calculate_soil_evaporation(soilevaporation)
      ! restrict transpiration to positive only
      transpiration = max(dble_zero,transpiration)

      ! if snow present assume that soilevaporation is sublimation of soil first
      snow_sublimation = dble_zero
      if (snow_storage > dble_zero) then
          snow_sublimation = soilevaporation
          if (snow_sublimation*deltat(n) > snow_storage) snow_sublimation = snow_storage * deltat_1(n)
          soilevaporation = soilevaporation - snow_sublimation
          snow_storage = snow_storage - (snow_sublimation * deltat(n))
      end if

      !!!!!!!!!!
      ! GPP (gC.m-2.day-1)
      !!!!!!!!!!

      ! reset output variable
      ci_global = dble_zero
      if (stomatal_conductance > vsmall) then
          FLUXES(n,1) = max(dble_zero,acm_gpp(stomatal_conductance))
      else
          FLUXES(n,1) = dble_zero
      endif

      !!!!!!!!!!
      ! Update water balance
      !!!!!!!!!!

      ! determine potential evaporation which sources its water from the soil,
      ! i.e. excluding wet canopy surface
      ET_pot = transpiration + soilevaporation

      ! add any snow melt to the rainfall now that we have already dealt with the canopy interception
      rainfall = rainfall + (snow_melt / seconds_per_step)
      ! do mass balance (i.e. is there enough water to support ET)
      call calculate_update_soil_water(transpiration,soilevaporation,((rainfall-intercepted_rainfall)*seconds_per_day),ET_net)

      !!!!!!!!!!
      ! Assign output variables
      !!!!!!!!!!

      ! store soil water content of the top soil layer 0-10 cm (mm)
      POOLS(n,1) = 1d3*soil_waterfrac(1)*layer_thickness(1)
!      POOLS(n,1) = sum(1d3*soil_waterfrac(1:3)*layer_thickness(1:3))

      ! Restrict transpiration to positive only
      FLUXES(n,2) = transpiration
      ! Wet canopy surface evaporation
      FLUXES(n,3) = wetcanopy_evap
      ! pass soil evaporation to output
      FLUXES(n,4) = soilevaporation + snow_sublimation
      ! pass soil surface runoff and underflow (drainage) from soil column kg/m2/day
      FLUXES(n,5) = runoff
      FLUXES(n,6) = underflow
      ! estimate mean leaf water potential based on ratio of actual vs potential water supply
      FLUXES(n,7) = minlwp * (transpiration / (max_supply * mmol_to_kg_water))
      if (max_supply <= vsmall) FLUXES(n,7) = dble_zero
      ! output diagnositic, internal leaf CO2 concentration
      FLUXES(n,8) = ci_global

      !!!!!!!!!!
      ! Bug checking
      !!!!!!!!!!

!      do nxp = 1, nopools
!         if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < 0d0 .or. &
!             abs(POOLS(n+1,nxp)) == abs(log(infi))) then
!            print*,"step",n,"POOL",nxp
!            print*,"met",met(:,n)
!            print*,"POOLS",POOLS(n,:)
!            print*,"FLUXES",FLUXES(n,:)
!            print*,"POOLS+1",POOLS(n+1,:)
!            print*,"wSWP",wSWP
!            print*,"waterfrac",soil_waterfrac
!            stop
!         endif
!      enddo
! 
!      do nxp = 1, nofluxes
!         if (FLUXES(n,nxp) /= FLUXES(n,nxp) .or. abs(FLUXES(n,nxp)) == abs(log(infi))) then
!            print*,"step",n,"FLUX",nxp
!            print*,"met",met(:,n)
!            print*,"POOLS",POOLS(n,:)
!            print*,"FLUXES",FLUXES(n,:)
!            print*,"POOLS+1",POOLS(n+1,:)
!            print*,"wSWP",wSWP
!            print*,"waterfrac",soil_waterfrac
!            stop
!         endif
!      enddo

    end do ! nodays loop

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm_gpp(gs)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosynthesis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    ! declare input variables
    double precision, intent(in) :: gs

    ! declare local variables
    double precision :: pn, pd, pp, qq, ci, mult, pl &
                       ,gc ,gs_mol, gb_mol

    ! Temperature adjustments for Michaelis-Menten coefficients
    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
    ! See McMurtrie et al., (1992) Australian Journal of Botany, vol 40, 657-677
    co2_half_sat   = arrhenious(kc_saturation,kc_half_sat_conc,leafT)
    co2_comp_point = arrhenious(co2comp_saturation,co2comp_half_sat_conc,leafT)

    !
    ! Metabolic limited photosynthesis
    !

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1)
    pn = lai*avN*NUE*opt_max_scaling(pn_max_temp,pn_opt_temp,pn_kurtosis,leafT)

    !
    ! Diffusion limited photosynthesis
    !

    ! daily canopy conductance (mmolH2O.m-2.s-1-> molCO2.m-2.day-1)
    ! The ratio of H20:CO2 diffusion is 1.646259 (Jones appendix 2).
    ! i.e. gcH2O*1.646259 = gcCO2
    gs_mol = gs * 1d-3 * seconds_per_day * gs_H2O_CO2
    ! canopy level boundary layer conductance unit change
    ! (m.s-1 -> mol.m-2.day-1) assuming sea surface pressure only.
    ! Note the ratio of H20:CO2 diffusion through leaf level boundary layer is
    ! 1.37 (Jones appendix 2).
    gb_mol = aerodynamic_conductance * seconds_per_day * convert_ms1_mol_1 * gb_H2O_CO2
    ! Combining in series the stomatal and boundary layer conductances
    gc = (gs_mol ** (-dble_one) + gb_mol ** (-dble_one)) ** (-dble_one)

    ! pp and qq represent limitation by metabolic (temperature & N) and
    ! diffusion (co2 supply) respectively
    pp = (pn/umol_to_gC)/gc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm or umol/mol)
    mult = co2+qq-pp
    ci = 0.5d0*(mult+sqrt((mult*mult)-4d0*(co2*qq-pp*co2_comp_point)))
    ci = min(ci,co2) ! C3 can't have more CO2 than is in the atmosphere
    ci_global = ci
    ! calculate CO2 limited rate of photosynthesis (gC.m-2.day-1)
    pd = (gc * (co2-ci)) * umol_to_gC
    ! scale to day light period as this is then consistent with the light
    ! capture period (1/24 = 0.04166667)
    pd = pd * dayl_hours * 0.04166667d0

    !
    ! Light limited photosynthesis
    !

    ! calculate light limted rate of photosynthesis (gC.m-2.day-1)
    pl = e0 * canopy_par_MJday

    !
    ! CO2 and light co-limitation
    !

    ! calculate combined light and CO2 limited photosynthesis
    acm_gpp = pl*pd/(pl+pd)

    ! don't forget to return
    return

  end function acm_gpp
  !
  !----------------------------------------------------------------------
  !
  double precision function find_gs(gs_in)

    ! Calculate CO2 limited photosynthesis as a function of metabolic limited
    ! photosynthesis (pn), atmospheric CO2 concentration and stomatal
    ! conductance (gs_in). Photosynthesis is calculated twice to allow for
    ! testing of senstivity to iWUE (iWUE).

    ! arguments
    double precision, intent(in) :: gs_in

    ! local variables
    double precision :: tmp,airt_save,lw_save, &
                        isothermal,deltaTemp,deltaR
    double precision :: gs_high, gs_store, &
                        gpp_high, gpp_low, &
                        evap_high, evap_low

    if (do_iWUE) then

        !!!!!!!!!!
        ! Optimise intrinsic water use efficiency
        !!!!!!!!!!

        ! if desired calculate the steady-state energy balance
        if (do_energy_balance) then
            ! save values which will need to be reset
            airt_save = leafT ; lw_save = canopy_lwrad_Wm2
            ! estimate energy balance without wet evaporation effects
            isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)
            call update_net_radiation(isothermal,leafT,lai,dble_one &
                                     ,gs_in,aerodynamic_conductance,vpd_pa &
                                     ,deltaTemp,deltaR)
            ! note that both the leafT and canopy LW have an implicit day -> day length correction
            canopy_lwrad_Wm2 = canopy_lwrad_Wm2 + deltaR
            leafT = leafT + deltaTemp
        endif
        ! estimate photosynthesis with current estimate of gs
        gpp_low = acm_gpp(gs_in)

        ! Increment gs
        gs_high = gs_in + delta_gs
        ! if desired calculate the steady-state energy balance
        if (do_energy_balance) then
            leafT = airt_save ; canopy_lwrad_Wm2 = lw_save
            ! estimate energy balance without wet evaporation effects
            isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)
            call update_net_radiation(isothermal,leafT,lai,dble_one &
                                     ,gs_in,aerodynamic_conductance,vpd_pa &
                                     ,deltaTemp,deltaR)
            ! note that both the leafT and canopy LW have an implicit day -> day length correction
            canopy_lwrad_Wm2 = canopy_lwrad_Wm2 + deltaR
            leafT = leafT + deltaTemp
        endif
        ! estimate photosynthesis with incremented gs
        gpp_high = acm_gpp(gs_high)

        ! determine impact of gs increment on pd and how far we are from iWUE
        find_gs = iWUE - ((gpp_high - gpp_low)/lai)
!        find_gs = iWUE - (gpp_high - gpp_low)

    else ! iWUE = .true. / .false.

        !!!!!!!!!!
        ! Optimise water use efficiency
        !!!!!!!!!!

        ! Globally stored upper stomatal conductance estimate in memory
        gs_store = stomatal_conductance
        ! now assign the current estimate
        stomatal_conductance = gs_in
        if (do_energy_balance) then
            ! save values which will need to be reset
            airt_save = leafT ; lw_save = canopy_lwrad_Wm2
            ! estimate energy balance without wet evaporation effects
            isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)
            call update_net_radiation(isothermal,leafT,lai,dble_one &
                                     ,gs_in,aerodynamic_conductance,vpd_pa &
                                     ,deltaTemp,deltaR)
            ! note that both the leafT and canopy LW have an implicit day -> day length correction
            canopy_lwrad_Wm2 = canopy_lwrad_Wm2 + deltaR
            leafT = leafT + deltaTemp
        endif
        ! estimate photosynthesis with current estimate of gs
        gpp_low = acm_gpp(gs_in)
        call calculate_transpiration(evap_low)

        ! Increment gs
        gs_high = gs_in + delta_gs
        ! now assign the incremented estimate
        stomatal_conductance = gs_high
        ! if desired calculate the steady-state energy balance
        if (do_energy_balance) then
            leafT = airt_save ; canopy_lwrad_Wm2 = lw_save
            ! estimate energy balance without wet evaporation effects
            isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)
            call update_net_radiation(isothermal,leafT,lai,dble_one &
                                     ,gs_in,aerodynamic_conductance,vpd_pa &
                                     ,deltaTemp,deltaR)
            ! note that both the leafT and canopy LW have an implicit day -> day length correction
            canopy_lwrad_Wm2 = canopy_lwrad_Wm2 + deltaR
            leafT = leafT + deltaTemp
        endif
        ! estimate photosynthesis with incremented gs
        gpp_high = acm_gpp(gs_high)
        call calculate_transpiration(evap_high)

        ! estimate marginal return on GPP for water loss, less water use efficiency criterion (gC.kgH2O-1.m-2.s-1)
        find_gs = ((gpp_high - gpp_low)/(evap_high - evap_low)) / lai
!        find_gs = ((gpp_high - gpp_low)/(evap_high - evap_low))
        find_gs = find_gs - iWUE

        ! return original stomatal value back into memory
        stomatal_conductance = gs_store

    end if ! iWUE = .true. / .false.

    ! now if I have been changing these drivers, best put them back to normal
    if (do_energy_balance) then
        leafT = airt_save ; canopy_lwrad_Wm2 = lw_save
    endif

    ! remember to return back to the user
    return

  end function find_gs
  !
  !------------------------------------------------------------------
  !
  subroutine acm_albedo_gc(deltaWP,Rtot)

    ! Determines 1) an approximation of canopy conductance (gc) mmolH2O.m-2.s-1
    ! based on potential hydraulic flow, air temperature and absorbed radiation.
    ! 2) calculates absorbed shortwave radiation (W.m-2) as function of LAI

    implicit none

    ! arguments
    double precision, intent(in) :: deltaWP, & ! minlwp-wSWP (MPa)
                                       Rtot    ! total hydraulic resistance (MPa.s-1.m-2.mmol-1)

    ! local variables
    double precision :: denom, isothermal, deltaTemp, deltaR
    double precision, parameter :: max_gs = 2500d0, & ! mmolH2O.m-2.s-1
                                   min_gs = 0.0001d0, & !
                                   tol_gs = 4d0!4d0       ! 

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (deltaWP > vsmall) then
       ! Determine potential water flow rate (mmolH2O.m-2.dayl-1)
       max_supply = (deltaWP/Rtot) * seconds_per_day
    else
       ! set minimum (computer) precision level flow
       max_supply = vsmall
    end if

    if (lai > vsmall) then

        ! there is lai therefore we have have stomatal conductance

        ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
        ! maximum possible evaporation for the day.
        ! This will then be reduced based on CO2 limits for diffusion based
        ! photosynthesis
        denom = slope * ((canopy_swrad_MJday * 1d6 * dayl_seconds_1) + canopy_lwrad_Wm2) &
              + (air_density_kg * cpair * vpd_pa * 1d-3 * aerodynamic_conductance)
        denom = (denom / (lambda * max_supply * mmol_to_kg_water * dayl_seconds_1)) - slope
        denom = denom / psych
        stomatal_conductance = aerodynamic_conductance / denom

        ! convert m.s-1 to mmolH2O.m-2.s-1
        stomatal_conductance = stomatal_conductance * 1d3 * convert_ms1_mol_1
        ! if conditions are dew forming then set conductance to maximum as we are not going to be limited by water demand
        if (stomatal_conductance <= dble_zero .or. stomatal_conductance > max_gs) stomatal_conductance = max_gs

        ! if we are potentially limited by stomatal conductance or we are using instrinsic water use efficiency (rather than WUE)
        ! then iterate to find optimum gs otherwise just go with the max...
        if (stomatal_conductance /= max_gs .or. do_iWUE ) then
            ! If there is a positive demand for water then we will solve for photosynthesis limits on gs through iterative solution
            delta_gs = 1d-3*lai ! mmolH2O/m2leaf/day
            stomatal_conductance = zbrent('acm_albedo_gc:find_gs',find_gs,min_gs,stomatal_conductance,tol_gs)
        end if

        ! if desired calculate the steady-state energy balance
        if (do_energy_balance) then
            isothermal = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)
            call update_net_radiation(isothermal,leafT,lai,dble_one &
                                     ,stomatal_conductance,aerodynamic_conductance,vpd_pa &
                                     ,deltaTemp,deltaR)
            ! note that both the leafT and canopy LW have an implicit day -> day length correction
            canopy_lwrad_Wm2 = canopy_lwrad_Wm2 + deltaR
            leafT = leafT + deltaTemp
        endif

    else

        ! if no LAI then there can be no stomatal conductance
        stomatal_conductance = dble_zero    
 
    endif ! if LAI > vsmall

  end subroutine acm_albedo_gc
  !
  !------------------------------------------------------------------
  !
  subroutine acm_meteorological_constants(input_temperature)

    ! Determine some multiple use constants

    implicit none

    ! arguments
    double precision, intent(in) :: input_temperature

    ! local variables
    double precision :: s, mult

    ! Density of air (kg/m3)
    air_density_kg = 353d0/(input_temperature+freeze)
    ! Conversion ratio for m.s-1 -> mol.m-2.s-1
    convert_ms1_mol_1 = const_sfc_pressure / ((input_temperature+freeze)*Rcon)
    ! latent heat of vapourisation,
    ! function of air temperature (J.kg-1)
    if (input_temperature < dble_one) then
        lambda = 2.835d6
    else
        lambda = 2501000d0-2364d0*input_temperature
    endif
    ! psychrometric constant (kPa K-1)
    psych = (0.0646d0*exp(0.00097d0*input_temperature))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult = input_temperature+237.3d0
    ! 2502.935945 = 0.61078*17.269*237.3
    s = 2502.935945d0*exp(17.269d0*input_temperature/mult)
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = s/(mult*mult)

  end subroutine acm_meteorological_constants
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_transpiration(transpiration)

    ! Models leaf cnaopy transpiration based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(out) :: transpiration ! kgH2O.m-2.day-1

    ! local variables
    double precision :: canopy_radiation & ! isothermal net radiation (W/m2)
                           ,water_supply & ! Potential water supply to canopy from soil (kgH2O.m-2.day-1)
                                  ,gs,gb   ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! calculate potential water supply (kgH2O.m-2.day-1)
    ! provided potential upper bound on evaporation
    water_supply = max_supply * mmol_to_kg_water

    ! Change units of potential stomatal conductance
    ! (mmolH2O.m-2.s-1 -> m.s-1).
    ! Note assumption of sea surface pressure only
    gs = stomatal_conductance / (convert_ms1_mol_1 * 1d3)
    ! Combine in series stomatal conductance with boundary layer
    gb = aerodynamic_conductance

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kg.m-2.day-1)
    transpiration = (slope*canopy_radiation) + (air_density_kg*cpair*vpd_pa*1d-3*gb)
    ! Calculate the transpiration flux and restrict by potential water supply
    ! over the day
    transpiration = min(water_supply,(transpiration / (lambda*(slope+(psych*(dble_one+gb/gs)))))*dayl_seconds)

  end subroutine calculate_transpiration
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_wetcanopy_evaporation(wetcanopy_evap,act_pot_ratio,storage,transpiration)

    ! Estimates evaporation of canopy intercepted rainfall based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(in) :: transpiration      ! kgH2O/m2/day
    double precision, intent(inout) :: storage         ! canopy water storage kgH2O/m2
    double precision, intent(out) :: wetcanopy_evap, & ! kgH2O.m-2.day-1
                                      act_pot_ratio    ! Ratio of potential evaporation to actual

    ! local variables
    double precision :: day, night
    double precision :: canopy_radiation, & ! isothermal net radiation (W/m2)
                                      gb   ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! Combine in series stomatal conductance with boundary layer
    gb = aerodynamic_conductance

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * seconds_per_day_1)
!    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * dayl_seconds_1)

    !!!!!!!!!!
    ! Calculate canopy evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! Calculate numerator of Penman Montheith (kgH2O.m-2.day-1)
    wetcanopy_evap = (slope*canopy_radiation) + (air_density_kg*cpair*vpd_pa*1d-3*gb)
    ! Calculate the potential wet canopy evaporation, limited by energy used for
    ! transpiration
    wetcanopy_evap = (wetcanopy_evap / (lambda*(slope+psych))) * seconds_per_day !dayl_seconds
    ! substract transpiration from potential surface evaporation
    wetcanopy_evap = wetcanopy_evap - transpiration

    ! dew is unlikely to occur (if we had energy balance) if mint > 0
    if (wetcanopy_evap < dble_zero .and. mint > dble_zero) wetcanopy_evap = dble_zero
    ! Sublimation is unlikely to occur (if we had energy balance) if maxt < 0
    if (wetcanopy_evap > dble_zero .and. maxt < dble_zero) wetcanopy_evap = dble_zero

    ! Remember potential evaporation to later calculation of the potential
    ! actual ratio
    act_pot_ratio = wetcanopy_evap

    ! assuming there is any rainfall, currently water on the canopy or dew formation
    if (rainfall > dble_zero .or. storage > dble_zero .or. wetcanopy_evap < dble_zero) then
        ! Update based on canopy water storage
        call canopy_interception_and_storage(wetcanopy_evap,storage)
    else
        ! there is no water movement possible
        intercepted_rainfall = dble_zero ; wetcanopy_evap = dble_zero
    endif

    ! now calculate the ratio of potential to actual evaporation
    if (act_pot_ratio == dble_zero) then
        act_pot_ratio = dble_zero
    else
        act_pot_ratio = abs(wetcanopy_evap / act_pot_ratio)
    endif

  end subroutine calculate_wetcanopy_evaporation
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_soil_evaporation(soilevap)

    ! Estimate soil surface evaporation based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).

    implicit none

    ! arguments
    double precision, intent(out) :: soilevap ! kgH2O.m-2.day-1

    ! local variables
    double precision :: soil_radiation & ! isothermal net radiation (W/m2)
                      ,water_diffusion & ! Diffusion of water through soil matrix (m.s-1)
                                ,esurf & ! see code below
                                 ,esat & ! soil air space saturation vapour pressure
                                  ,gws & ! water vapour conductance through soil air space (m.s-1)
                                   ,Qc

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    soil_radiation = soil_lwrad_Wm2 + (soil_swrad_MJday * 1d6 * dayl_seconds_1)
    ! estimate ground heat flux from statistical approximation, positive if energy moving up profile
    ! NOTE: linear coefficient estimates from SPA simulations
    Qc = -0.4108826d0 * (maxt - maxt_lag1)
    soil_radiation = soil_radiation + Qc

    !!!!!!!!!!
    ! Calculate soil evaporative fluxes (kgH2O/m2/day)
    !!!!!!!!!!

    ! calculate saturated vapour pressure (kPa), function of temperature.
    esat = 0.1d0 * exp( 1.80956664d0 + ( 17.2693882d0 * (maxt+freeze) - 4717.306081d0 ) / ( maxt+freeze - 35.86d0 ) )
    air_vapour_pressure = esat - (vpd_pa * 1d-3)

    ! Estimate water diffusion rate (m2.s-1) Jones (2014) appendix 2
    water_diffusion = 24.2d-6 * ( (maxt+freeze) / 293.2d0 )**1.75d0
    ! Soil conductance to water vapour diffusion (m s-1)...
    gws = porosity(1) * water_diffusion / (tortuosity*drythick)

    ! vapour pressure in soil airspace (kPa), dependent on soil water potential
    ! - Jones p.110. partial_molar_vol_water
    esurf = esat * exp( 1d6 * SWP(1) * partial_molar_vol_water / ( Rcon * (maxt+freeze) ) )
    ! now difference in vapour pressure between soil and canopy air spaces
    esurf = esurf - air_vapour_pressure

    ! Estimate potential soil evaporation flux (kgH2O.m-2.day-1)
    soilevap = (slope*soil_radiation) + (air_density_kg*cpair*esurf*soil_conductance)
    soilevap = soilevap / (lambda*(slope+(psych*(dble_one+soil_conductance/gws))))
    soilevap = soilevap * dayl_seconds

    ! dew is unlikely to occur (if we had energy balance) if mint > 0
!    if (soilevap < dble_zero .and. mint > dble_one) soilevap = dble_zero
    ! Sublimation is unlikely to occur (if we had energy balance) if maxt < 0
!    if (soilevap > dble_zero .and. maxt < dble_one) soilevap = dble_zero

  end subroutine calculate_soil_evaporation
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_aerodynamic_conductance

    !
    ! Calculates the aerodynamic or bulk canopy conductance (m.s-1). Here we
    ! assume neutral conditions due to the lack of an energy balance calculation
    ! in either ACM or DALEC. The equations used here are with SPA at the time
    ! of the calibration
    !

    implicit none

    ! local variables
    double precision :: mixing_length_momentum, & ! mixing length parameter for momentum (m)
                         length_scale_momentum    ! length scale parameter for momentum (m)

    ! calculate the zero plane displacement and roughness length
    call z0_displacement(ustar_Uh)
    ! calculate friction velocity at tower height (reference height ) (m.s-1)
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
!    ustar = (wind_spd / log((tower_height-displacement)/roughl)) * vonkarman
    ustar = wind_spd * ustar_Uh
    ! both length scale and mixing length are considered to be constant within
    ! the canopy (under dense canopy conditions) calculate length scale (lc)
    ! for momentum absorption within the canopy; Harman & Finnigan (2007)
    ! and mixing length (lm) for vertical momentum within the canopy Harman & Finnigan (2008)
    if (lai > min_lai) then
        length_scale_momentum = (4d0*canopy_height) / lai
        mixing_length_momentum = max(canopy_height*0.02d0, 2d0*(ustar_Uh**3)*length_scale_momentum)
    else
        length_scale_momentum = vonkarman * tower_height
        mixing_length_momentum = canopy_height * vonkarman
    endif

    ! based on Harman & Finnigan (2008); neutral conditions only
    call log_law_decay

    ! calculate soil surface conductance
    call calculate_soil_conductance(mixing_length_momentum)

    ! now we are interested in the within canopy wind speed,
    ! here we assume that the wind speed just inside of the canopy is most important.
    canopy_wind = canopy_wind*exp((ustar_Uh*((canopy_height*0.75d0)-canopy_height))/mixing_length_momentum)

    ! calculate leaf level conductance (m/s) for water vapour under forced convective conditions
    call average_leaf_conductance(aerodynamic_conductance)

  end subroutine calculate_aerodynamic_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine average_leaf_conductance(gv_forced)

    !
    ! Subroutine calculates the forced conductance of water vapour for non-cylinder within canopy leaves (i.e. broadleaf)
    ! Free convection (i.e. that driven by energy balance) is negelected here due to the lack of an energy balance
    ! calculation in DALEC. Should a energy balance be added then this code could be expanded include free conductance
    ! Follows a simplified approach to that used in SPA (Smallman et al 2013).
    !

    implicit none

    ! arguments
    double precision, intent(out) :: gv_forced ! canopy conductance (m/s) for water vapour under forced convection

    ! local parameters
    double precision, parameter :: leaf_width = 0.08d0   & ! leaf width (m)
                                          ,Pr = 0.72d0     ! Prandtl number
    ! local variables
    double precision :: Dwv & ! Diffusion coefficient of water in air (m2.s-1); air temperature and pressure dependant
                              ! variables for the more advanced
                              ! boundary conditions
            ,nusselt_forced & ! Nusselt value under forced convection
         ,dynamic_viscosity & ! dynamic viscosity (kg.m-2.s-1)
       ,kinematic_viscosity & ! kinematic viscosity (m2.s-1)
                 ,Sh_forced & ! Sherwood number under forced convection
                        ,Re   ! Reynolds number

    ! Determine diffusion coefficient (m2s-1), temperature dependant (pressure dependence neglected). Jones p51;
    ! 0.0000242 = conversion to make diffusion specific for water vapor (um2.s-1)
    Dwv = 0.0000242d0*(((leafT+freeze)/293.15d0)**1.75d0)
    ! Calculate the dynamic viscosity of air
    dynamic_viscosity = (((leafT+freeze)**1.5d0)/((leafT+freeze)+120d0))*1.4963d-6
    kinematic_viscosity = dynamic_viscosity/air_density_kg
    Re = (leaf_width*canopy_wind)/kinematic_viscosity
    ! calculate nusselt value under forced convection conditions
    nusselt_forced = (1.18d0*(Pr**(0.33d0))*(sqrt(Re)))
    ! update specific Sherwood numbers
    Sh_forced = 0.962d0*nusselt_forced
    ! This is the forced conductance of water vapour for the current leaf
    gv_forced = ((Dwv*Sh_forced)/leaf_width)*0.5d0
    ! apply lai correction
    gv_forced = gv_forced * lai

  end subroutine average_leaf_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine log_law_decay

    ! Standard log-law above canopy wind speed (m.s-1) decay under neutral
    ! conditions.
    ! See Harman & Finnigan 2008; Jones 1992 etc for details.

    implicit none

    ! log law decay
    canopy_wind = (ustar / vonkarman) * log((canopy_height-displacement) / roughl)

    ! set minimum value for wind speed at canopy top (m.s-1)
    canopy_wind = max(min_wind,canopy_wind)

  end subroutine log_law_decay
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_field_capacity

    ! field capacity calculations for saxton eqns !

    implicit none

    ! local variables..
    integer        :: i
    double precision :: x1, x2

    x1 = 0.1d0 ; x2 = 0.7d0 ! low/high guess
    do i = 1 , nos_soil_layers+1
       water_retention_pass = i
       ! field capacity is water content at which SWP = -10 kPa
       field_capacity(i) = zbrent('water_retention:water_retention_saxton_eqns', &
                                   water_retention_saxton_eqns , x1 , x2 , 0.001d0 )
    enddo

  end subroutine calculate_field_capacity
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_daylength(doy,lat)

    ! Subroutine uses day of year and latitude (-90 / 90 degrees) as inputs,
    ! combined with trigonomic functions to calculate
    ! 1) day length in hours and seconds
    ! 2) hour of sunrise and sunset
    ! 3) cosine of solar zenith angle to allow scaling of evaporation over course of day

    implicit none

    ! arguments
    double precision, intent(in) :: doy, lat

    ! local variables
    double precision :: dec, mult, sinld, cosld, aob

    !
    ! Estimate solar geometry variables needed
    !

    ! declination
!    dec = - asin( sin( 23.45d0 * deg_to_rad ) * cos( 2d0 * pi * ( doy + 10d0 ) / 365d0 ) )
    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) / 365d0 ) )
    ! latitude in radians
    mult = lat * deg_to_rad
    ! day length is estimated as the ratio of sin and cos of the product of declination an latitude in radiation
    sinld = sin( mult ) * sin( dec )
    cosld = cos( mult ) * cos( dec )
    aob = max(-dble_one,min(dble_one,sinld / cosld))

    ! estimate day length in hours and seconds and upload to module variables
    dayl_hours = 12d0 * ( dble_one + 2d0 * asin( aob ) * pi_1 )
    dayl_seconds = dayl_hours * seconds_per_hour

    ! estimate sun rise and run set hours
    sunrise = 12 - nint(dayl_hours*0.5d0) ; sunset = sunrise + nint(dayl_hours)

    ! estimate the solar cosine zenith angle for 12 noon
    cos_solar_zenith_angle = sinld + cosld

    ! return to user
    return

  end subroutine calculate_daylength
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_longwave_isothermal(canopy_temperature,soil_temperature)

    ! Subroutine estimates the isothermal net longwave radiation (W.m-2) for
    ! the canopy and soil surface. SPA uses a complex multi-layer radiative
    ! transfer scheme including reflectance, transmittance any absorption.
    ! However, for a given canopy vertical profiles, the LAI absorption
    ! relationship is readily predicted via Michaelis-Menten or
    ! non-rectangular hyperbola as done here.

    implicit none

    ! arguments
    double precision, intent(in) :: canopy_temperature, soil_temperature ! oC

    ! local variables
    double precision :: lwrad, & ! downward long wave radiation from sky (W.m-2)
        longwave_release_soil, & ! emission of long wave radiation from surfaces per m2
      longwave_release_canopy, & ! assuming isothermal condition (W.m-2)
            trans_lw_fraction, &
        reflected_lw_fraction, &
         absorbed_lw_fraction, &
      canopy_release_fraction, & ! fraction of longwave emitted from within the canopy to ultimately be released
   canopy_absorption_from_sky, & ! canopy absorbed radiation from downward LW (W.m-2)
  canopy_absorption_from_soil, & ! canopy absorbed radiation from soil surface (W.m-2)
                  canopy_loss, & ! longwave radiation released from canopy surface (W.m-2).
                                 ! i.e. this value is released from the top and the bottom
     soil_absorption_from_sky, & ! soil absorbed radiation from sky (W.m-2)
  soil_absorption_from_canopy    ! soil absorbed radiation emitted from canopy (W.m-2)

    ! estimate long wave radiation from atmosphere (W.m-2)
    lwrad = emiss_boltz * (maxt+freeze-20d0) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_soil = emiss_boltz * (soil_temperature+freeze) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release_canopy = emiss_boltz * (canopy_temperature+freeze) ** 4

    !!!!!!!!!!
    ! Determine fraction of longwave absorbed by canopy and returned to the sky
    !!!!!!!!!!

    ! calculate fraction of longwave radiation coming from the sky to pentrate to the soil surface
    trans_lw_fraction = dble_one - (max_lai_lwrad_transmitted*lai)/(lai+lai_half_lwrad_transmitted)
    ! calculate the fraction of longwave radiation from sky which is reflected back into the sky
    reflected_lw_fraction = (max_lai_lwrad_reflected*lai) / (lai+lai_half_lwrad_reflected)
    ! calculate absorbed longwave radiation coming from the sky
    absorbed_lw_fraction = dble_one - trans_lw_fraction - reflected_lw_fraction
    ! Calculate the potential absorption of longwave radiation lost from the
    ! canopy to soil / sky
    canopy_release_fraction = dble_one - (max_lai_lwrad_release*lai) / (lai+lai_half_lwrad_release)

    !!!!!!!!!!
    ! Distribute longwave from sky
    !!!!!!!!!!

    ! long wave absorbed by the canopy from the sky
    canopy_absorption_from_sky = lwrad * absorbed_lw_fraction
    ! Long wave absorbed by soil from the sky, soil absorption assumed to be equal to emissivity
    soil_absorption_from_sky = trans_lw_fraction * lwrad * emissivity
    ! Long wave reflected directly back into sky
    sky_lwrad_Wm2 = lwrad * reflected_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave from soil
    !!!!!!!!!!

    ! First, calculate longwave radiation coming up from the soil plus the radiation which is reflected
    canopy_absorption_from_soil = longwave_release_soil + (trans_lw_fraction * lwrad * (dble_one-emissivity))
    ! Second, use this total to estimate the longwave returning to the sky
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + (canopy_absorption_from_soil * trans_lw_fraction)
    ! Third, now calculate the longwave from the soil surface absorbed by the canopy
    canopy_absorption_from_soil = canopy_absorption_from_soil * absorbed_lw_fraction

    !!!!!!!!!!
    ! Distribute longwave originating from the canopy itself
    !!!!!!!!!!

    ! calculate two-sided long wave radiation emitted from canopy which is
    ! ultimately lost from to soil or sky (i.e. this value is used twice, once
    ! to soil once to sky)
    canopy_loss = longwave_release_canopy * lai * canopy_release_fraction
    ! Calculate longwave absorbed by soil which is released by the canopy itself
    soil_absorption_from_canopy = canopy_loss * emissivity
    ! Canopy released longwave returned to the sky
    sky_lwrad_Wm2 = sky_lwrad_Wm2 + canopy_loss

    !!!!!!!!!!
    ! Isothermal net long wave canopy and soil balance (W.m-2)
    !!!!!!!!!!

    ! determine isothermal net canopy. Note two canopy_loss used to account for
    ! upwards and downwards emissions
    canopy_lwrad_Wm2 = (canopy_absorption_from_sky + canopy_absorption_from_soil) - (canopy_loss + canopy_loss)
    ! determine isothermal net soil
    soil_lwrad_Wm2 = (soil_absorption_from_sky + soil_absorption_from_canopy) - longwave_release_soil

  end subroutine calculate_longwave_isothermal
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_shortwave_balance

    ! Subroutine estimates the canopy and soil absorbed shortwave radiation (MJ/m2/day).
    ! Radiation absorption is paritioned into NIR and PAR for canopy, and NIR +
    ! PAR for soil.

    ! SPA uses a complex multi-layer radiative transfer scheme including
    ! reflectance, transmittance any absorption. However, for a given
    ! canopy vertical profiles, the LAI absorption relationship is readily
    ! predicted via Michaelis-Menten or non-rectangular hyperbola as done here.

    implicit none

    ! local variables
    double precision :: balance                    &
                       ,absorbed_nir_fraction_soil &
                       ,absorbed_par_fraction_soil &
                       ,fsnow                      &
                       ,soil_par_MJday             &
                       ,soil_nir_MJday             &
                       ,trans_nir_MJday            &
                       ,trans_par_MJday            &
                       ,canopy_nir_MJday           &
                       ,refl_par_MJday             &
                       ,refl_nir_MJday             &
                       ,reflected_nir_fraction     & !
                       ,reflected_par_fraction     & !
                       ,absorbed_nir_fraction      & !
                       ,absorbed_par_fraction      & !
                       ,trans_nir_fraction         & !
                       ,trans_par_fraction

    ! local parameters
    double precision, parameter :: newsnow_nir_abs = 0.27d0 & ! NIR absorption fraction
                                  ,newsnow_par_abs = 0.05d0   ! PAR absorption fraction

    !!!!!!!!!!
    ! Determine canopy absorption / reflectance as function of LAI
    !!!!!!!!!!

    ! Canopy transmitted of PAR & NIR radiation towards the soil
    trans_par_fraction = dble_one - (lai*max_lai_par_transmitted) &
                       / (lai+lai_half_par_transmitted)
    trans_nir_fraction = dble_one - (lai*max_lai_nir_transmitted) &
                       / (lai+lai_half_nir_transmitted)
    ! Canopy reflected of near infrared and photosynthetically active radiation
    reflected_nir_fraction = (lai*max_lai_nir_reflection) &
                          / (lai+lai_half_nir_reflection)
    reflected_par_fraction = (lai*max_lai_par_reflection) &
                          / (lai+lai_half_par_reflection)
    ! Canopy absorption of near infrared and photosynthetically active radiation
    absorbed_nir_fraction = dble_one - reflected_nir_fraction - trans_nir_fraction
    absorbed_par_fraction = dble_one - reflected_par_fraction - trans_par_fraction

    !!!!!!!!!!
    ! Estimate canopy absorption of incoming shortwave radiation
    !!!!!!!!!!

    ! Estimate incoming shortwave radiation absorbed, transmitted and reflected by the canopy (MJ.m-2.day-1)
    canopy_par_MJday = (sw_par_fraction * swrad * absorbed_par_fraction)
    canopy_nir_MJday = ((dble_one - sw_par_fraction) * swrad * absorbed_nir_fraction)
    trans_par_MJday = (sw_par_fraction * swrad * trans_par_fraction)
    trans_nir_MJday = ((dble_one - sw_par_fraction) * swrad * trans_nir_fraction)
    refl_par_MJday = (sw_par_fraction * swrad * reflected_par_fraction)
    refl_nir_MJday = ((dble_one - sw_par_fraction) * swrad * reflected_nir_fraction)

    !!!!!!!!!
    ! Estimate soil absorption of shortwave passing through the canopy
    !!!!!!!!!

    ! Update soil reflectance based on snow cover
    if (snow_storage > dble_zero) then
        fsnow = dble_one - exp( - snow_storage * 1d-2 )  ! fraction of snow cover on the ground
        absorbed_par_fraction_soil = ((dble_one - fsnow) * soil_swrad_absorption) + (fsnow * newsnow_par_abs)
        absorbed_nir_fraction_soil = ((dble_one - fsnow) * soil_swrad_absorption) + (fsnow * newsnow_nir_abs)
    else
        absorbed_par_fraction_soil = soil_swrad_absorption
        absorbed_nir_fraction_soil = soil_swrad_absorption
    endif

    ! Then the radiation incident and ultimately absorbed by the soil surface itself (MJ.m-2.day-1)
    soil_par_MJday = trans_par_MJday * absorbed_par_fraction_soil
    soil_nir_MJday = trans_nir_MJday * absorbed_nir_fraction_soil
    ! combine totals for use is soil evaporation
    soil_swrad_MJday = soil_nir_MJday + soil_par_MJday

    !!!!!!!!!
    ! Estimate canopy absorption of soil reflected shortwave radiation
    ! This additional reflection / absorption cycle is needed to ensure > 0.99
    ! of incoming radiation is explicitly accounted for in the energy balance.
    !!!!!!!!!

    ! Update the canopy radiation absorption based on the reflected radiation (MJ.m-2.day-1)
    canopy_par_MJday = canopy_par_MJday + ((trans_par_MJday-soil_par_MJday) * absorbed_par_fraction)
    canopy_nir_MJday = canopy_nir_MJday + ((trans_nir_MJday-soil_nir_MJday) * absorbed_nir_fraction)
    ! Update the total radiation reflected back into the sky, i.e. that which is now transmitted through the canopy
    refl_par_MJday = refl_par_MJday + ((trans_par_MJday-soil_par_MJday) * trans_par_fraction)
    refl_nir_MJday = refl_nir_MJday + ((trans_nir_MJday-soil_nir_MJday) * trans_nir_fraction)

    ! Combine to estimate total shortwave canopy absorbed radiation
    canopy_swrad_MJday = canopy_par_MJday + canopy_nir_MJday

    ! check energy balance
    balance = swrad - canopy_par_MJday - canopy_nir_MJday - refl_par_MJday - refl_nir_MJday - soil_swrad_MJday
!    if ((balance - swrad) / swrad > 0.01) then
!        print*,"SW residual frac = ",(balance - swrad) / swrad,"SW residual = ",balance,"SW in = ",swrad
!    endif

  end subroutine calculate_shortwave_balance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_Rtot(Rtot)

    ! Purpose of this subroutine is to calculate the minimum soil-root hydraulic
    ! resistance input into ACM. The approach used here is identical to that
    ! found in SPA.

    ! declare inputs
    double precision,intent(inout) :: Rtot ! MPa.s-1.m-2.mmol-1

    ! local variables
    integer :: i
    double precision :: bonus, sum_water_flux, &
                        soilR1,soilR2,transpiration_resistance,root_reach_local, &
                        root_depth_50
    double precision, dimension(nos_root_layers) :: root_mass    &
                                                   ,root_length  &
                                                   ,ratio
    double precision, parameter :: root_depth_frac_50 = 0.25d0 ! fractional soil depth above which 50 %
                                                               ! of the root mass is assumed to be located

    ! reset water flux
    water_flux = dble_zero ; wSWP = dble_zero
    ratio = dble_zero ; ratio(1) = dble_one ; root_mass = dble_zero
    ! calculate soil depth to which roots reach
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
!    transpiration_resistance = (gplant * lai)**(-dble_one)
    transpiration_resistance = canopy_height / (gplant * lai)

    !!!!!!!!!!!
    ! calculate current steps soil hydraulic conductivity
    !!!!!!!!!!!

    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i,soil_waterfrac(i),soil_conductivity(i))
    end do ! soil layers

    !!!!!!!!!!!
    ! Calculate root profile
    !!!!!!!!!!!

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 3 root layer system this can be estimates more simply

    ! top 25 % of root profile
    root_depth_50 = root_reach * root_depth_frac_50
    if (root_depth_50 <= layer_thickness(1)) then
        ! Greater than 50 % of the fine root biomass can be found in the top
        ! soil layer

        ! Start by assigning all 50 % of root biomass to the top soil layer
        root_mass(1) = root_biomass * 0.5d0
        ! Then quantify how much additional root is found in the top soil layer
        ! assuming that the top 25 % depth is found somewhere within the top
        ! layer
        bonus = (root_biomass-root_mass(1)) &
              * (layer_thickness(1)-root_depth_50) / (root_reach - root_depth_50)
        root_mass(1) = root_mass(1) + bonus
        ! partition the remaining root biomass between the seconds and third
        ! soil layers
        if (root_reach > sum(layer_thickness(1:2))) then
            root_mass(2) = (root_biomass - root_mass(1)) &
                         * (layer_thickness(2)/(root_reach-layer_thickness(1)))
            root_mass(3) = root_biomass - sum(root_mass(1:2))
        else
            root_mass(2) = root_biomass - root_mass(1)
        endif
    else if (root_depth_50 > layer_thickness(1) .and. root_depth_50 <= sum(layer_thickness(1:2))) then
        ! Greater than 50 % of fine root biomass found in the top two soil
        ! layers. We will divide the root biomass uniformly based on volume,
        ! plus bonus for the second layer (as done above)
        root_mass(1) = root_biomass * 0.5d0 * (layer_thickness(1)/root_depth_50)
        root_mass(2) = root_biomass * 0.5d0 * ((root_depth_50-layer_thickness(1))/root_depth_50)
        ! determine bonus for the seconds layer
        bonus = (root_biomass-sum(root_mass(1:2))) &
              * ((sum(layer_thickness(1:2))-root_depth_50)/(root_reach-root_depth_50))
        root_mass(2) = root_mass(2) + bonus
        root_mass(3) = root_biomass - sum(root_mass(1:2))
    else
        ! Greater than 50 % of fine root biomass stock spans across all three
        ! layers
        root_mass(1) = root_biomass * 0.5d0 * (layer_thickness(1)/root_depth_50)
        root_mass(2) = root_biomass * 0.5d0 * (layer_thickness(2)/root_depth_50)
        root_mass(3) = root_biomass - sum(root_mass(1:2))
    endif
    ! now convert root mass into lengths
    root_length = root_mass * root_mass_length_coef_1
!    root_length = root_mass / (root_density * root_cross_sec_area)

    !!!!!!!!!!!
    ! Calculate hydraulic properties and each rooted layer
    !!!!!!!!!!!

    ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
    root_reach_local = min(root_reach,layer_thickness(1))
    soilR1=soil_resistance(root_length(1),root_reach_local,soil_conductivity(1)*head_1)
    soilR2=root_resistance(root_mass(1),root_reach_local)
    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    ! also note that the head correction is now added rather than
    ! subtracted in SPA equations because deltaWP is soilWP-minlwp not
    ! soilWP prior to application of minlwp
    demand = abs(minlwp-SWP(1:nos_root_layers))+head*canopy_height
    water_flux(1) = demand(1)/(transpiration_resistance + soilR1 + soilR2)

    ! second root layer
    if (root_mass(2) > dble_zero) then
        root_reach_local = min(root_reach,layer_thickness(2))
        soilR1=soil_resistance(root_length(2),root_reach_local,soil_conductivity(2)*head_1)
        soilR2=root_resistance(root_mass(2),root_reach_local)
        water_flux(2) = demand(2)/(transpiration_resistance + soilR1 + soilR2)
    endif ! roots present in the seconds layer?

    ! Bottom root layer
    if (root_mass(3) > dble_zero ) then
       ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
       soilR1=soil_resistance(root_length(3),layer_thickness(3),soil_conductivity(3)*head_1)
       soilR2=root_resistance(root_mass(3),layer_thickness(3))
       ! calculate and accumulate steady state water flux in mmol.m-2.s-1
       water_flux(3) = demand(3)/(transpiration_resistance + soilR1 + soilR2)
    endif ! roots present in third layer?
    ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))

    ! if freezing then assume soil surface is frozen
    if (meant < dble_one) then
        water_flux(1) = dble_zero
        ratio(1) = dble_zero
        ratio(2:nos_root_layers) = layer_thickness(2:nos_root_layers) / sum(layer_thickness(2:nos_root_layers))
    endif
    ! calculate sum value
    sum_water_flux = sum(water_flux)

    ! calculate weighted SWP and uptake fraction
    wSWP = sum(SWP(1:nos_root_layers) * water_flux(1:nos_root_layers))
    uptake_fraction(1:nos_root_layers) = water_flux(1:nos_root_layers) / sum_water_flux
    wSWP = wSWP / sum_water_flux

    ! sanity check in case of zero flux
    if (sum_water_flux == dble_zero) then
        wSWP = -20d0
        uptake_fraction = dble_zero ; uptake_fraction(1) = dble_one
    endif

    ! determine effective resistance (MPa.s-1.m-2.mmol-1)
    Rtot = sum(demand) / sum(water_flux)

    ! finally convert transpiration flux (mmol.m-2.s-1)
    ! into kg.m-2.step-1 for consistency with ET in "calculate_update_soil_water"
    water_flux = water_flux * mmol_to_kg_water * seconds_per_step

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !-----------------------------------------------------------------
  !
  subroutine canopy_interception_and_storage(potential_evaporation,storage)

    ! Simple daily time step integration of canopy rainfall interception, runoff
    ! and rainfall (kgH2O.m-2.s-1). NOTE: it is possible for intercepted rainfall to be
    ! negative if stored water running off into the soil is greater than
    ! rainfall (i.e. when leaves have died between steps)

    implicit none

    ! arguments
    double precision, intent(inout) :: storage, & ! canopy water storage (kgH2O/m2)
                         potential_evaporation    ! wet canopy evaporation (kgH2O.m-2.day-1),
                                                  ! enters as potential but leaves as water balance adjusted
    ! local variables
    integer :: i, hr
    double precision :: a, through_fall, max_storage, max_storage_1, daily_addition, wetcanopy_evaporation &
                       ,potential_drainage_rate ,drain_rate, evap_rate, initial_canopy, co_mass_balance, dx, tmp1, tmp2, tmp3
    ! local parameters
    double precision, parameter :: CanIntFrac = -0.5d0,  & ! Coefficient scaling rainfall interception fraction with LAI
                                  CanStorFrac = 0.1d0,   & ! Coefficient scaling canopy water storage with LAI
                                 RefDrainRate = 0.002d0, & ! Reference drainage rate (mm/min; Rutter et al 1975)
                                  RefDrainLAI = 1.05d0,  & ! Reference drainage LAI (m2/m2; Rutter et al 1975)
                                 RefDrainCoef = 3.7d0,   & ! Reference drainage Coefficient (Rutter et al 1975)
                               RefDrainCoef_1 = RefDrainCoef ** (-dble_one)

    ! hold initial canopy storage in memory
    initial_canopy = storage
    ! determine maximum canopy storage & through fall fraction
    through_fall = max(min_throughfall,exp(CanIntFrac*lai))
    ! maximum canopy storage (mm); minimum is applied to prevent errors in
    ! drainage calculation. Assume minimum capacity due to wood
    max_storage = max(min_storage,CanStorFrac*lai) ; max_storage_1 = max_storage**(-dble_one)
    ! potential intercepted rainfall (kgH2O.m-2.s-1)
    intercepted_rainfall = rainfall * (dble_one - through_fall)
    
    ! calculate drainage coefficients (Rutter et al 1975); Corsican Pine
    ! 0.002 is canopy specific coefficient modified by 0.002*(max_storage/1.05)
    ! where max_storage is the canopy maximum capacity (mm) (LAI based) and
    ! 1.05 is the original canopy capacitance
    a = log( RefDrainRate * ( max_storage / RefDrainLAI ) ) - RefDrainCoef * max_storage

    ! average rainfall intercepted by canopy (kgH2O.m-2.day-1)
    daily_addition = intercepted_rainfall * seconds_per_day 

    ! reset cumulative variables
    through_fall = dble_zero ; wetcanopy_evaporation = dble_zero
    drain_rate = dble_zero ; evap_rate = dble_zero

    ! deal with rainfall additions first
    do i = 1, int(days_per_step)

       ! add rain to the canopy and overflow as needed
       storage = storage + daily_addition

       if (storage > max_storage) then

           if (potential_evaporation > dble_zero) then

               ! assume co-access to available water above max_storage by both drainage and
               ! evaporation. Water below max_storage is accessable by evaporation only. 

               ! Trapezium rule for approximating integral of drainage rate.
               ! Allows estimation of the mean drainage rate between starting point and max_storage,
               ! thus the time period appropriate for co-access can be quantified
               dx = storage - ((storage + max_storage)*0.5d0)
               tmp1 = exp(a + (RefDrainCoef * storage)) 
               tmp2 = exp(a + (RefDrainCoef * max_storage)) 
               tmp3 = exp(a + (RefDrainCoef * (storage+dx))) 
               potential_drainage_rate = 0.5d0 * dx * ((tmp1 + tmp2) + 2d0 * tmp3)
               potential_drainage_rate = potential_drainage_rate * 1440d0
              
               ! restrict evaporation and drainage to the quantity above max_storage
               evap_rate = potential_evaporation ; drain_rate = min(potential_drainage_rate,storage-max_storage)

               ! limit based on available water if total demand is greater than excess
               co_mass_balance = ((storage-max_storage) / (evap_rate + drain_rate))
               evap_rate = evap_rate * co_mass_balance ; drain_rate = drain_rate * co_mass_balance 

               ! estimate evaporation from remaining water, less that already removed from storage and evaporation energy used
               evap_rate = evap_rate + min(potential_evaporation - evap_rate, storage - evap_rate - drain_rate)

           else 

               ! load dew formation to the current local evap_rate variable
               evap_rate = potential_evaporation
               ! restrict drainage the quantity above max_storage, adding dew formation too
               drain_rate = (storage - evap_rate) - max_storage

           endif

       else
           ! no drainage just apply evaporation / dew formation fluxes directly
           evap_rate = potential_evaporation
           drain_rate = dble_zero
           if (evap_rate > dble_zero) then
               ! evaporation restricted by fraction of surface actually covered
               ! in water
               evap_rate = evap_rate * min(dble_one,storage * max_storage_1)
               ! and the total amount of water
               evap_rate = min(evap_rate,storage)
           else 
               ! then dew formation has occurred, if this pushes storage > max_storage add it to drainage
               drain_rate = max(dble_zero,(storage - evap_rate) - max_storage)
           endif ! evap_rate > 0
       endif ! storage > max_storage

       ! update canopy storage with water flux
       storage = max(dble_zero,storage - evap_rate - drain_rate)
       wetcanopy_evaporation = wetcanopy_evaporation + evap_rate
       through_fall = through_fall + drain_rate

    end do ! days

    ! correct intercepted rainfall rate to kgH2O.m-2.s-1
    intercepted_rainfall = intercepted_rainfall - ((through_fall * days_per_step_1) * seconds_per_day_1)

    ! sanity checks; note 1e-8 prevents precision errors causing flags
    if (intercepted_rainfall > rainfall .or. storage < dble_zero &
   .or. (wetcanopy_evaporation * days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day)) ) then
       print*,"Condition 1",intercepted_rainfall > rainfall
       print*,"Condition 2",storage < dble_zero
       print*,"Condition 3",(wetcanopy_evaporation * days_per_step_1) > (1d-8 + initial_canopy + (rainfall*seconds_per_day))
       print*,"storage (kgH2O/m2)",storage,"max_storage (kgH2O/m2)",max_storage,"initial storage (kgH2O/m2)", initial_canopy
       print*,"rainfall (kgH2O/m2/day)", rainfall*seconds_per_day, "through_fall (kgH2O/m2/day)", (through_fall * days_per_step_1)
       print*,"through_fall_total (kgH2O/m2/step)",through_fall
       print*,"potential_evaporation (kgH2O/m2/day)",potential_evaporation
       print*,"actual evaporation    (kgH2O/m2/day)",wetcanopy_evaporation * days_per_step_1
       stop
    endif

    ! average evaporative flux to daily rate (kgH2O/m2/day)
    potential_evaporation = wetcanopy_evaporation * days_per_step_1

    ! final clearance of canopy storage of version small values at the level of system precision
    if (storage < 10d0*vsmall) storage = dble_zero

  end subroutine canopy_interception_and_storage
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_update_soil_water(ET_leaf,ET_soil,rainfall_in,corrected_ET)

   !
   ! Function updates the soil water status and layer thickness
   ! Soil water profile is updated in turn with evaporative losses,
   ! rainfall infiltration and gravitational drainage
   ! Root layer thickness is updated based on changes in the rooting depth from
   ! the previous step
   !

   implicit none

   ! arguments
   double precision, intent(in) :: ET_leaf,ET_soil & ! evapotranspiration estimate (kgH2O.m-2.day-1)
                                      ,rainfall_in   ! rainfall (kgH2O.m-2.day-1)
   double precision, intent(out) :: corrected_ET     ! water balance corrected evapotranspiration (kgH2O/m2/day)

   ! local variables
   integer :: day
   double precision ::  depth_change, water_change
   double precision, dimension(nos_root_layers) :: avail_flux, evaporation_losses

   ! reset soil water exchanges
   underflow = dble_zero ; runoff = dble_zero ; corrected_ET = dble_zero

   ! to allow for smooth water balance integration carry this out at daily time step
   do day = 1, nint(days_per_step)

      !!!!!!!!!!
      ! Evaporative losses
      !!!!!!!!!!

      ! Assume leaf transpiration is drawn from the soil based on the
      ! update_fraction estimated in calculate_Rtot
      evaporation_losses = ET_leaf * uptake_fraction
      ! Assume all soil evaporation comes from the soil surface only
      evaporation_losses(1) = evaporation_losses(1) + ET_soil
      ! can not evaporate from soil more than is available (m -> mm)
      avail_flux = soil_waterfrac(1:nos_root_layers) * layer_thickness(1:nos_root_layers) * 1d3
      where (evaporation_losses > avail_flux) evaporation_losses = avail_flux * 0.999d0

      ! this will update the ET estimate outside of the function
      ! days_per_step corrections happens outside of the loop below
      corrected_ET = corrected_ET + sum(evaporation_losses)

      ! pass information to waterloss variable and zero watergain
      ! convert kg.m-2 (or mm) -> Mg.m-2 (or m)
      waterloss = dble_zero ; watergain = dble_zero
      waterloss(1:nos_root_layers) = evaporation_losses(1:nos_root_layers)*1d-3
      ! update soil water status with evaporative losses
      soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)) &
                                           + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                        / layer_thickness(1:nos_soil_layers)
      ! reset soil water flux variables
      waterloss = dble_zero ; watergain = dble_zero

      !!!!!!!!!!
      ! Gravitational drainage
      !!!!!!!!!!

      ! determine drainage flux between surface -> sub surface and sub surface
      call gravitational_drainage

      ! update soil water status with drainage
      soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)) &
                                           + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                        / layer_thickness(1:nos_soil_layers)
      ! reset soil water flux variables
      waterloss = dble_zero ; watergain = dble_zero

      !!!!!!!!!!
      ! Rainfall infiltration drainage
      !!!!!!!!!!

      ! determine infiltration from rainfall (kgH2O/m2/step),
      ! if rainfall is probably liquid / soil surface is probably not frozen
      if (rainfall_in > dble_zero) then
          call infiltrate(rainfall_in)
      else
          runoff = runoff + (rainfall_in * days_per_step_1)
      endif ! is there any rain to infiltrate?
      ! update soil profiles. Convert fraction into depth specific values (rather than m3/m3) then update fluxes
      soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)) &
                                           + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                        / layer_thickness(1:nos_soil_layers)
      ! reset soil water flux variables
      waterloss = dble_zero ; watergain = dble_zero

      ! mass balance check, at this point do not try and adjust evaporation to
      ! correct for lack of supply. Simply allow for drought in next time step
      ! instead...
      where (soil_waterfrac <= dble_zero)
             soil_waterfrac = vsmall
      end where

   end do ! days_per_step

   ! apply time step correction kgH2O/m2/step -> kgH2O/m2/day
   corrected_ET = corrected_ET * days_per_step_1
   underflow = underflow * days_per_step_1
   runoff = runoff * days_per_step_1

   !!!!!!!!!!
   ! Update soil layer thickness
   !!!!!!!!!!

   depth_change = dble_zero ; water_change = dble_zero
   ! if roots extent down into the bucket
   if (root_reach > (top_soil_depth+mid_soil_depth) .or. previous_depth > (top_soil_depth+mid_soil_depth)) then
      ! how much has root depth extended since last step?
      depth_change = root_reach - previous_depth

      ! if there has been an increase
      if (depth_change > dble_zero .and. root_reach > sum(layer_thickness(1:2))+min_layer) then

         ! determine how much water is within the new volume of soil
         water_change = soil_waterfrac(nos_soil_layers) * depth_change
         ! now assign that new volume of water to the deep rooting layer
         soil_waterfrac(nos_root_layers) = ((soil_waterfrac(nos_root_layers) * layer_thickness(nos_root_layers)) &
                                            + water_change) / (layer_thickness(nos_root_layers)+depth_change)
         ! explicitly update the soil profile if there has been rooting depth
         ! changes
         layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
         layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
         layer_thickness(4) = max_depth - sum(layer_thickness(1:3))


      elseif (depth_change < dble_zero .and. root_reach > layer_thickness(1)+min_layer) then

         ! determine how much water is lost from the old volume of soil
         water_change = soil_waterfrac(nos_root_layers) * abs(depth_change)
         ! now assign that new volume of water to the deep rooting layer
         soil_waterfrac(nos_soil_layers) = ((soil_waterfrac(nos_soil_layers) * layer_thickness(nos_soil_layers)) &
                                            + water_change) / (layer_thickness(nos_soil_layers)+abs(depth_change))

         ! explicitly update the soil profile if there has been rooting depth
         ! changes
         layer_thickness(1) = top_soil_depth ; layer_thickness(2) = mid_soil_depth
         layer_thickness(3) = max(min_layer,root_reach-sum(layer_thickness(1:2)))
         layer_thickness(4) = max_depth - sum(layer_thickness(1:3))

      else

         ! we don't want to do anything, just recycle the previous depth

      end if ! depth change

   end if ! root reach beyond top layer

   ! in all cases keep track of the previous rooted depth
   previous_depth = root_reach

   ! finally update soil water potential
   call soil_water_potential

!   ! sanity check for catastrophic failure
!   do soil_layer = 1, nos_soil_layers
!      if (soil_waterfrac(soil_layer) < 0d0 .and. soil_waterfrac(soil_layer) > -0.01d0) then
!          soil_waterfrac(soil_layer) = 0d0
!      endif
!      if (soil_waterfrac(soil_layer) < 0d0 .or. soil_waterfrac(soil_layer) /= soil_waterfrac(soil_layer)) then
!         print*,'ET',ET,"rainfall",rainfall_in
!         print*,'evaporation_losses',evaporation_losses
!         print*,"watergain",watergain
!         print*,"waterloss",waterloss
!         print*,'depth_change',depth_change
!         print*,"soil_waterfrac",soil_waterfrac
!         print*,"porosity",porosity
!         print*,"layer_thicknes",layer_thickness
!         print*,"Uptake fraction",uptake_fraction
!         print*,"max_depth",max_depth,"root_k",root_k,"root_reach",root_reach
!         print*,"fail" ; stop
!      endif
!   end do

   ! explicit return needed to ensure that function runs all needed code
   return

  end subroutine calculate_update_soil_water
  !
  !-----------------------------------------------------------------
  !
  subroutine infiltrate(rainfall_in)

    ! Takes surface_watermm and distrubutes it among top !
    ! layers. Assumes total infilatration in timestep.   !

    implicit none

    ! arguments
    double precision, intent(in) :: rainfall_in ! rainfall (kg.m-2.day-1)

    ! local argumemts
    integer :: i
    double precision    :: add   & ! surface water available for infiltration (m)
                          ,wdiff   ! available space in a given soil layer for water to fill (m)

    ! convert rainfall water from mm -> m (or kgH2O.m-2.day-1 -> MgH2O.m-2.day-1)
    add = rainfall_in * 1d-3

    do i = 1 , nos_soil_layers
       ! determine the available pore space in current soil layer
       wdiff = max(dble_zero,(porosity(i)-soil_waterfrac(i))*layer_thickness(i)-watergain(i)+waterloss(i))
       ! is the input of water greater than available space
       ! if so fill and subtract from input and move on to the next
       ! layer
       if (add > wdiff) then
          ! if so fill and subtract from input and move on to the next layer
          watergain(i) = watergain(i) + wdiff
          add = add - wdiff
       else
          ! otherwise infiltate all in the current layer
          watergain(i) = watergain(i) + add
          add = dble_zero
       end if
       ! if we have added all available water we are done
       if (add <= dble_zero) then
           add = dble_zero
           exit
       end if

    end do ! nos_soil_layers

    ! if after all of this we have some water left assume it is runoff (kgH2O.m-2.day-1)
    ! NOTE that runoff is reset outside of the daily soil loop
    runoff = runoff + (add * 1d3)

  end subroutine infiltrate
  !
  !-----------------------------------------------------------------
  !
  subroutine gravitational_drainage

    ! integrator for soil gravitational drainage !

    implicit none

    ! local variables..
    integer :: d, nos_integrate
    double precision  :: liquid & ! liquid water in local soil layer (m3/m3)
                    ,drainlayer & ! field capacity of local soil layer (m3/m3)
                         ,unsat & ! unsaturated pore space in soil_layer below the current (m3/m3) 
                        ,change & ! absolute volume of water drainage in current layer (m3)
                      ,drainage & ! drainage rate of current layer (m/day)
                   ,local_drain & ! drainage of current layer (m/nos_minutes)
      ,iceprop(nos_soil_layers)

    ! local parameters
    integer, parameter :: nos_hours_per_day = 1440, nos_minutes = 360

    ! calculate soil ice proportion; at the moment
    ! assume everything liquid
    iceprop = dble_zero
    ! except the surface layer in the mean daily temperature is < 0oC
    if (meant < dble_one) iceprop(1) = dble_one

    do soil_layer = 1, nos_soil_layers

       ! soil water capacity of the current layer
       drainlayer = field_capacity( soil_layer )
       ! liquid content of the soil layer
       liquid     = soil_waterfrac( soil_layer ) &
                  * ( dble_one - iceprop( soil_layer ) )

       ! initial conditions; i.e. is there liquid water and more water than
       ! layer can hold
       if ( liquid > drainlayer ) then

          ! unsaturated volume of layer below (m3 m-2)..
          unsat = max( dble_zero , ( porosity( soil_layer+1 ) - soil_waterfrac( soil_layer+1 ) ) &
                             * layer_thickness( soil_layer+1 ) / layer_thickness( soil_layer ) )

          d = 1 ; nos_integrate = nos_hours_per_day / nos_minutes
          drainage = dble_zero ; local_drain = dble_zero
          do while (d <= nos_integrate .and. liquid > drainlayer)
              ! estimate drainage rate (m/s) 
              call calculate_soil_conductivity(soil_layer,liquid,local_drain)
              ! scale to total number of seconds in increment
              local_drain = local_drain * dble(nos_minutes * 60)
              local_drain = min(liquid-drainlayer,local_drain)
              liquid = liquid - local_drain
              drainage = drainage + local_drain
              d = d + 1
          end do ! integrate over time

          ! layer below cannot accept more water than unsat
          if ( drainage > unsat ) drainage = unsat
          ! water loss from this layer (m3)
          change = drainage * layer_thickness(soil_layer)
          ! update soil layer below with drained liquid
          watergain( soil_layer + 1 ) = watergain( soil_layer + 1 ) + change
          waterloss( soil_layer     ) = waterloss( soil_layer     ) + change

       end if ! some liquid water and drainage possible

    end do ! soil layers

    ! estimate drainage from bottom of soil column (kgH2O/m2/day)
    ! NOTES: that underflow is reset outside of the daily soil loop
    underflow = underflow + (waterloss(nos_soil_layers) * 1d3)

  end subroutine gravitational_drainage
  !
  !-----------------------------------------------------------------
  !
  subroutine soil_porosity(soil_frac_clay,soil_frac_sand)

   ! Porosity is estimated from Saxton equations. !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand
    ! local variables..
    double precision, parameter :: H = 0.332d0, &
                                   J = -7.251d-4, &
                                   K = 0.1276d0

    ! loop over soil layers..
    porosity(1:nos_soil_layers) = H + J * soil_frac_sand(1:nos_soil_layers) + &
                                  K * log10(soil_frac_clay(1:nos_soil_layers))
    ! then assign same to core layer
    porosity(nos_soil_layers+1) = porosity(nos_soil_layers)

  end subroutine soil_porosity
  !
  !---------------------------------------------------------------------
  !
  subroutine initialise_soils(soil_frac_clay,soil_frac_sand)

    !
    ! Subroutine calculate the soil layers field capacities and sets the initial
    ! soil water potential set to field capacity
    !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    integer :: i

    ! include some hardcoded boundaries for the Saxton equations
    where (soil_frac_sand < 5d0) soil_frac_sand = 5d0
    where (soil_frac_clay < 5d0) soil_frac_clay = 5d0
    where (soil_frac_clay > 60d0) soil_frac_clay = 60d0
    ! calculate soil porosity (m3/m3)
    call soil_porosity(soil_frac_clay,soil_frac_sand)
    ! calculate field capacity (m3/m-3)
    call calculate_field_capacity
    ! calculate initial soil water fraction
    soil_waterfrac = field_capacity
    ! calculate initial soil water potential
    SWP = dble_zero
    call soil_water_potential
    ! seperately calculate the soil conductivity as this applies to each layer
    do i = 1, nos_soil_layers
       call calculate_soil_conductivity(i,soil_waterfrac(i),soil_conductivity(i))
    end do ! soil layers
    ! but apply the lowest soil layer to the core as well in initial conditions
    soil_conductivity(nos_soil_layers+1) = soil_conductivity(nos_soil_layers)

    ! final sanity check for porosity
    where (porosity <= field_capacity) porosity = field_capacity + 0.01d0

  end subroutine initialise_soils
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_soil_conductivity(soil_layer,waterfrac,conductivity)

    ! Calculate the soil conductivity (m s-1) of water based on soil
    ! characteristics and current water content

    implicit none

    ! arguments
    integer, intent(in) :: soil_layer
    double precision, intent(in) :: waterfrac
    double precision, intent(out) :: conductivity

    ! soil conductivity for the dynamic soil layers (i.e. not including core)
    conductivity = cond1(soil_layer) * exp(cond2(soil_layer)+cond3(soil_layer)/waterfrac)

    ! protection against floating point error
    if (waterfrac < 0.05d0) conductivity = 1d-30

  end subroutine calculate_soil_conductivity
  !
  !------------------------------------------------------------------
  !
  subroutine saxton_parameters(soil_frac_clay,soil_frac_sand)

    ! Calculate the key parameters of the Saxton, that is cond1,2,3 !
    ! and potA,B                                                    !

    implicit none

    ! arguments
    double precision, dimension(nos_soil_layers) :: soil_frac_clay &
                                                   ,soil_frac_sand

    ! local variables
    integer :: i
    double precision, parameter :: A = -4.396d0,  B = -0.0715d0,       CC = -4.880d-4, D = -4.285d-5, &
                                   E = -3.140d0,  F = -2.22d-3,         G = -3.484d-5, H = 0.332d0,   &
                                   J = -7.251d-4, K = 0.1276d0,         P = 12.012d0,  Q = -7.551d-2, &
                                   R = -3.895d0,  T = 3.671d-2,         U = -0.1103d0, V = 8.7546d-4, &
                                   mult1 = 100d0, mult2 = 2.778d-6

    ! layed out in this manor to avoid memory management issues in module
    ! variables
    potA(1:nos_soil_layers) = A + (B * soil_frac_clay) + &
                            (CC * soil_frac_sand * soil_frac_sand) + &
                            (D * soil_frac_sand * soil_frac_sand * soil_frac_clay)
    potA(1:nos_soil_layers) = exp(potA(1:nos_soil_layers))
    potA(1:nos_soil_layers) = potA(1:nos_soil_layers) * mult1

    potB(1:nos_soil_layers) = E + (F * soil_frac_clay * soil_frac_clay) + &
                             (G * soil_frac_sand * soil_frac_sand * soil_frac_clay)

    cond1(1:nos_soil_layers) = mult2
    cond2(1:nos_soil_layers) = P + (Q * soil_frac_sand)
    cond3(1:nos_soil_layers) = R + (T * soil_frac_sand) + (U * soil_frac_clay) + &
                              (V * soil_frac_clay * soil_frac_clay)

    ! assign bottom of soil column value to core
    potA(nos_soil_layers+1)  = potA(nos_soil_layers)
    potB(nos_soil_layers+1)  = potB(nos_soil_layers)
    cond1(nos_soil_layers+1) = mult2
    cond2(nos_soil_layers+1) = cond2(nos_soil_layers)
    cond3(nos_soil_layers+1) = cond3(nos_soil_layers)

  end subroutine saxton_parameters
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_soil_conductance(lm)

    ! proceedsure to solve for soil surface resistance based on Monin-Obukov
    ! similarity theory stability correction momentum & heat are integrated
    ! through the under canopy space and canopy air space to the surface layer
    ! references are Nui & Yang 2004; Qin et al 2002
    ! NOTE: conversion to conductance at end

    implicit none

    ! declare arguments
    double precision, intent(in) :: lm

    ! local variables
    double precision :: canopy_decay & ! canopy decay coefficient for soil exchange
                       ,Kh_canht       ! eddy diffusivity at canopy height (m2.s-1)

    ! parameters
    double precision, parameter :: foliage_drag = 0.2d0 ! foliage drag coefficient

    ! calculate eddy diffusivity at the top of the canopy (m2.s-1)
    ! Kaimal & Finnigan 1994; for near canopy approximation
    Kh_canht=vonkarman*ustar*(canopy_height-displacement)

    ! calculate canopy decay coefficient with stability correction
    ! NOTE this is not consistent with canopy momentum decay done by Harman &
    ! Finnigan (2008)
    canopy_decay = (((foliage_drag*canopy_height*max(min_lai,lai))/lm)**0.5d0)

    ! approximation of integral for soil resistance
    soil_conductance = canopy_height/(canopy_decay*Kh_canht) &
                     * (exp(canopy_decay*(dble_one-(soil_roughl/canopy_height)))- &
                        exp(canopy_decay*(dble_one-((roughl+displacement)/canopy_height))))

    ! convert resistance (s.m-1) to conductance (m.s-1)
    soil_conductance = soil_conductance ** (-dble_one)

  end subroutine calculate_soil_conductance
  !
  !----------------------------------------------------------------------
  !
  subroutine soil_water_potential

    ! Find SWP without updating waterfrac yet (we do that in !
    ! waterthermal). Waterfrac is m3 m-3, soilwp is MPa.     !

    implicit none

    ! reformulation aims to remove if statement within loop to hopefully improve
    ! optimisation
    SWP(1:nos_soil_layers) = -0.001d0 * potA(1:nos_soil_layers) &
                           * soil_waterfrac(1:nos_soil_layers)**potB(1:nos_soil_layers)
    where (SWP(1:nos_soil_layers) < -20d0) SWP(1:nos_soil_layers) = -20d0

  end subroutine soil_water_potential
  !
  !------------------------------------------------------------------
  !
  subroutine update_net_radiation(isothermal,tempC,area_scaling,act_pot_ratio &
                                 ,sfc_exchange,aero_exchange,vapour_gradient,deltaTemp,deltaR)

    ! Use steady state solution of evaporation, convective (sensible) and radiative heat loss
    ! to update isothermal net radiation to net.
    ! Area scaling (e.g. lai) is an input to allow for common useage for soil (neglecting ground heat) and canopy

    ! arguments
    double precision, intent(in) ::      tempC, & ! input surface / air temperature (oC)
                                    isothermal, & ! isothermal net radiation (SW+LW; W/m2)
                                  area_scaling, & ! area scaling to apply (m2/m2)
                                 act_pot_ratio, & ! ratio of potential to actual evaporation, i.e. (avail / potenial)
                                  sfc_exchange, & ! surface exchange conductance (m/s; e.g. stomatal conductance)
                                 aero_exchange, & ! aerodynamic exchange conductance (m/s; e.g. aerodynamic conductance)
                               vapour_gradient    ! vapour pressure gradient (Pa; either VPD or between air and soil)
    double precision, intent(out) :: deltaTemp, & ! surface temperature difference (K)
                                     deltaR    ! surface longwave radiation difference (W/m2); subtract from isothermal longwave

    ! local variables
    double precision ::  tempK, & ! ambient temperature as K
          heat_loss_resistance, & ! resistance to heat loss from radiative and convection (s/m)
        aerodynamic_resistance, & ! aerodynamic resistance to water or heat exchangce (s/m)
           stomatal_resistance, & ! stomatal resistance to water exchangce (s/m)
              water_resistance, & ! serial combination of resistances to water evaporation
 thermal_gains, thermal_losses

    ! ambient temperature C -> K
    tempK = tempC + freeze

    !
    ! Calculate resistance to heat loss (s/m)
    !

    ! First estimate radiative loss term, initially calculated as conductance)
    heat_loss_resistance = 4d0 * emissivity * boltz * tempK ** 3 / (air_density_kg * cpair)
    ! Combine in parallel with convective conductances with area correction
    heat_loss_resistance = heat_loss_resistance + (2d0 * aero_exchange / area_scaling)
    ! convert from conductance m/s to s/m
    heat_loss_resistance = heat_loss_resistance ** (-1d0)

    !
    ! Convert aerodynamic and stomatal conductances to reisistance of water flux
    !

    aerodynamic_resistance = (aero_exchange/area_scaling) ** (-1d0)
    if (sfc_exchange == dble_zero) then
        ! if being used for surface water flux
        stomatal_resistance = dble_zero
    else
        ! if used for transpiration
        stomatal_resistance = (sfc_exchange/area_scaling) ** (-1d0)
    endif

    !
    ! Estimate thermal gains and losses (K) to calculate temperature difference
    !

    water_resistance = (aerodynamic_resistance + stomatal_resistance)
    thermal_gains = (heat_loss_resistance * water_resistance * psych * (isothermal/area_scaling)) &
                  / (air_density_kg * cpair * ((psych*water_resistance) + (slope*heat_loss_resistance)))
    thermal_losses = (heat_loss_resistance * vapour_gradient * 1d-3) &
                   / ((psych*water_resistance) + (slope*heat_loss_resistance))
    ! determine surface temperature difference (K); should be added to the canopy temperature
    deltaTemp = thermal_gains - thermal_losses
    ! apply actual potential ratio to scale wet surface evaporation when the
    ! supply of water is limited
    deltaTemp = deltaTemp * act_pot_ratio

    ! estimate update between isothermal to net radiation (W/m2), including area correction
    ! note that this MUST be added from the longwave component outside of this function
    deltaR = -4d0 * emissivity * boltz * tempK ** 3 * ( deltaTemp )
    deltaR = deltaR * area_scaling

    ! return to user
    return

  end subroutine update_net_radiation
  !
  !------------------------------------------------------------------
  !
  subroutine z0_displacement(ustar_Uh)

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    implicit none

    ! arguments
    double precision, intent(out) :: ustar_Uh ! ratio of friction velocity over wind speed at canopy top
    ! local variables
    double precision  sqrt_cd1_lai &
                     ,local_lai &
                     ,phi_h       ! roughness sublayer influence function
    double precision, parameter :: cd1 = 7.5d0,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003d0, & ! Substrate drag coefficient
                                    Cr = 0.3d0,   & ! Roughness element drag coefficient
!                          ustar_Uh_max = 0.3,   & ! Maximum observed ratio of
                                                   ! (friction velocity / canopy top wind speed) (m.s-1)
                          ustar_Uh_max = 1d0, ustar_Uh_min = 0.2d0, &
                                    Cw = 2d0      ! Characterises roughness sublayer depth (m)

    ! assign new value to min_lai to avoid max min calls
    local_lai = max(min_lai,lai)
    sqrt_cd1_lai = sqrt(cd1 * local_lai)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    displacement = (dble_one-((dble_one-exp(-sqrt_cd1_lai))/sqrt_cd1_lai))*canopy_height

    ! calculate estimate of ratio of friction velocity / canopy wind speed; with
    ! max value set at
    ustar_Uh = max(ustar_Uh_min,min(sqrt(Cs+Cr*local_lai*0.5d0),ustar_Uh_max))
    ! calculate roughness sublayer influence function;
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    phi_h = 0.19314718056d0
!    phi_h = log(Cw)-dble_one+Cw**(-dble_one) ! DO NOT FORGET TO UPDATE IF Cw CHANGES

    ! finally calculate roughness length, dependant on displacement, friction
    ! velocity and lai.
    roughl = ((dble_one-displacement/canopy_height)*exp(-vonkarman*ustar_Uh-phi_h))*canopy_height

    ! sanity check
!    if (roughl /= roughl) then
!        write(*,*)"TLS:  ERROR roughness length calculations"
!        write(*,*)"Roughness lenght", roughl, "Displacement", displacement
!        write(*,*)"canopy height", canopy_height, "lai", lai
!    endif

  end subroutine z0_displacement
  !
  !------------------------------------------------------------------
  !
  !------------------------------------------------------------------
  ! Functions other than the primary ACM and ACM ET are stored
  ! below this line.
  !------------------------------------------------------------------
  !
  !------------------------------------------------------------------
  !
  pure function arrhenious( a , b , t )

    ! The equation is simply...                        !
    !    a * exp( b * ( t - 25.0 ) / ( t + 273.15 ) )  !
    ! However, precision in this routine matters as it !
    ! affects many others. To maximise precision, the  !
    ! calculations have been split & d0 has been used. !

    implicit none

    ! arguments..
    double precision,intent(in) :: a , b , t
    double precision            :: arrhenious

    ! local variables..
    double precision :: answer, denominator, numerator

    numerator   = t - 25d0
    denominator = t + freeze
    answer      = a * exp( b * dble_one * numerator / denominator )
    arrhenious  = answer

  end function arrhenious
  !
  !----------------------------------------------------------------------
  !
  double precision function opt_max_scaling( max_val , optimum , kurtosis , current )

    ! estimates a 0-1 scaling based on a skewed guassian distribution with a
    ! given optimum, maximum and kurtosis

    implicit none

    ! arguments..
    double precision,intent(in) :: max_val, optimum, kurtosis, current

    ! local variables..
    double precision :: dummy

    if ( current >= max_val ) then
       opt_max_scaling = dble_zero
    else
       dummy     = ( max_val - current ) / ( max_val - optimum )
       dummy     = exp( log( dummy ) * kurtosis * ( max_val - optimum ) )
       opt_max_scaling = dummy * exp( kurtosis * ( current - optimum ) )
    end if

  end function opt_max_scaling
  !
  !------------------------------------------------------------------
  !
  double precision function root_resistance(root_mass,thickness)

   !
   ! Calculates root hydraulic resistance (MPa m2 s mmol-1) in a soil-root zone
   !

   implicit none

   ! arguments
   double precision :: root_mass, & ! root biomass in layer (gbiomass)
                       thickness    ! thickness of soil zone roots are in

   ! calculate root hydraulic resistance
   root_resistance = root_resist / (root_mass*thickness)

   ! return
   return

  end function root_resistance
  !
  !-----------------------------------------------------------------
  !
  double precision function soil_resistance(root_length,thickness,soilC)

    !
    ! Calculates the soil hydraulic resistance (MPa m2 s mmol-1) for a given
    ! soil-root zone
    !

    implicit none

    ! arguments
    double precision :: root_length, & ! root length in soil layer (m)
                          thickness, & ! thickness of soil layer (m)
                              soilC    ! soil conductivity m2.s-1.MPa-1

    ! local variables
    double precision :: rs, rs2

    ! calculate
    rs  = (root_length*pi)**(-0.5d0)
    rs2 = log( rs * root_radius_1 ) / (two_pi*root_length*thickness*soilC)
    ! soil water resistance
    soil_resistance = rs2*1d-9*mol_to_g_water

    ! return
    return

  end function soil_resistance
  !
  !------------------------------------------------------------------
  !
  double precision function water_retention_saxton_eqns( xin )

    ! field capacity calculations for saxton eqns !

    implicit none

    ! arguments..
    double precision, intent(in) :: xin

    ! local variables..
    double precision ::soil_wp

    ! calculate the soil water potential (MPa)..
    ! note that some modifications to scaling values have been made compared to
    ! SPA src to reduce computational cost
!    soil_wp = -0.001 * potA( water_retention_pass ) * xin**potB( water_retention_pass )
!    water_retention_saxton_eqns = -1000.0 * soil_wp + 10.0    ! 10 kPa represents air-entry swp
    soil_wp = potA( water_retention_pass ) * xin**potB( water_retention_pass )
    water_retention_saxton_eqns = -1d0 * soil_wp + 10d0    ! 10 kPa represents air-entry swp

    return

  end function water_retention_saxton_eqns
  !
  !------------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  ! Generic mathematical functions such as bisection and intergrator proceedures
  ! are stored below here
  !------------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------
  !
  double precision function zbrent( called_from , func , x1 , x2 , tol )

    ! This is a bisection routine. When ZBRENT is called, we provide a    !
    !  reference to a particular function and also two values which bound !
    !  the arguments for the function of interest. ZBRENT finds a root of !
    !  the function (i.e. the point where the function equals zero), that !
    !  lies between the two bounds.                                       !
    ! For a full description see Press et al. (1986).                     !

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    double precision,intent(in) :: tol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
       double precision function func( xval )
         double precision, intent(in) :: xval
       end function func
    end interface

    ! local variables..
    integer            :: iter
    integer,parameter  :: ITMAX = 30
    double precision   :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    double precision,parameter :: EPS = 3d-8

    ! calculations...
    a  = x1
    b  = x2
    fa = func( a )
    fb = func( b )

    ! Check that we haven't (by fluke) already started with the root..
    if ( fa .eq. 0d0 ) then
      zbrent = a
      return
    elseif ( fb .eq. 0d0 ) then
      zbrent = b
      return
    end if
    ! Ensure the supplied x-values give y-values that lie either
    ! side of the root and if not flag an error message...
    if ( sign(1d0,fa) .eq. sign(1d0,fb) ) then
       fa = func( a )
       fb = func( b )
       ! tell me otherwise what is going on
!       print*,"Supplied values must bracket the root of the function.",new_line('x'),  &
!         "     ","You supplied x1:",x1,new_line('x'),                     &
!         "     "," and x2:",x2,new_line('x'),                             &
!         "     "," which give function values of fa :",fa,new_line('x'),  &
!         "     "," and fb:",fb," .",new_line('x'),                        &
!         " zbrent was called by: ",trim(called_from)
!       fa = func( a )
!       fb = func( b )
    end if
    c = b
    fc = fb

    do iter = 1 , ITMAX

       ! If the new value (f(c)) doesn't bracket
       ! the root with f(b) then adjust it..
       if ( sign(1d0,fb) .eq. sign(1d0,fc) ) then
          c  = a
          fc = fa
          d  = b - a
          e  = d
       end if
       if ( abs(fc) .lt. abs(fb) ) then
          a  = b
          b  = c
          c  = a
          fa = fb
          fb = fc
          fc = fa
       end if
       tol1 = 2d0 * EPS * abs(b) + 0.5d0 * tol
       xm   = 0.5d0 * ( c - b )
       if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0d0 ) ) then
          zbrent = b
          return
       end if
       if ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) then
          s = fb / fa
          if ( a .eq. c ) then
             p = 2d0 * xm * s
             q = 1d0 - s
          else
             q = fa / fc
             r = fb / fc
             p = s * ( 2d0 * xm * q * ( q - r ) - ( b - a ) * ( r - 1d0 ) )
             q = ( q - 1d0 ) * ( r - 1d0 ) * ( s - 1d0 )
          end if
          if ( p .gt. 0d0 ) q = -q
          p = abs( p )
          if ( (2d0*p) .lt. min( 3d0*xm*q-abs(tol1*q) , abs(e*q) ) ) then
             e = d
             d = p / q
          else
             d = xm
             e = d
          end if
       else
          d = xm
          e = d
       end if
       a  = b
       fa = fb
       if ( abs(d) .gt. tol1 ) then
          b = b + d
       else
          b = b + sign( tol1 , xm )
       end if
       fb = func(b)
    enddo

!    print*,"zbrent has exceeded maximum iterations",new_line('x'),&
!           "zbrent was called by: ",trim(called_from)

    zbrent = b

  end function zbrent
  !
  !------------------------------------------------------------------
  !
!
!--------------------------------------------------------------------
!
end module CARBON_MODEL_MOD
