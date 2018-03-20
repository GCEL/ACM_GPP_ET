
module CARBON_MODEL_MOD

implicit none

! make all private
private

! explicit publics
public :: CARBON_MODEL        &
         ,dble_one,dble_zero  &
         ,wSWP_time           &
         ,nos_soil_layers     &
         ,soil_frac_clay      &
         ,soil_frac_sand

!!!!!!!!!
! Parameters
!!!!!!!!!

! useful technical parameters
double precision, parameter :: xacc = 1d-4        & ! accuracy parameter for zbrent bisection proceedure ! 0.0001
                              ,dble_zero = 0d0    &
                              ,dble_one = 1d0     &
                              ,vsmall = tiny(0d0)

integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
double precision, parameter :: pi = 3.1415927d0,  &
                             pi_1 = pi**(-dble_one),   &
                              pi2 = pi**2,        &
                           two_pi = pi*2d0,       &
                       deg_to_rad = pi/180d0,     &
              sin_dayl_deg_to_rad = sin( 23.45d0 * deg_to_rad ), & ! repeated function in acm
                          gravity = 9.8067d0,     & ! acceleration due to gravity, ms-1
                            boltz = 5.670400d-8,  & ! Boltzmann constant (W.m-2.K-4)
                       emissivity = 0.96d0,       &
                      emiss_boltz = emissivity * boltz, &
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
                    kc_saturation = 310d0,        & ! CO2 half saturation, saturation value
                 kc_half_sat_conc = 23.956d0,     & ! CO2 half sat, half sat
               co2comp_saturation = 36.5d0,       & ! CO2 compensation point, saturation
            co2comp_half_sat_conc = 9.46d0          ! CO2 comp point, half sat
                                                    ! Each of these are temperature
                                                    ! sensitivty
! hydraulic parameters
double precision, parameter :: &
                       tortuosity = 2.5d0,          & ! tortuosity
                           gplant = 5d0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1)
                      root_resist = 25d0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                      root_radius = 0.00029d0,    & ! root radius (m) Bonen et al 2014 = 0.00029
                                                    ! Williams et al 1996 = 0.0001
                    root_radius_1 = root_radius**(-dble_one), &
              root_cross_sec_area = 3.141593d-08, & ! root cross sectional area (m2)
                                                    ! = pi * root_radius * root_radius
                     root_density = 0.31d6,       & ! root density (g biomass m-3 root)
                                                    ! 0.5e6 Williams et al 1996
                                                    ! 0.31e6 Bonan et al 2014
          root_mass_length_coef_1 = (root_cross_sec_area * root_density)**(-1.d0), &
               const_sfc_pressure = 101325d0,     & ! (Pa)  Atmospheric surface pressure
                             head = 0.009807d0,   & ! head of pressure (MPa/m)
                           head_1 = 101.968d0       ! inverse head of pressure (m/MPa)

! structural parameters
double precision, parameter :: &
                    canopy_height = 9d0,          & ! canopy height assumed to be 9 m
                     tower_height = canopy_height + 2d0, & ! tower (observation) height assumed to be 2 m above canopy
                         min_wind = 0.1d0,        & ! minimum wind speed at canopy top
                     min_drythick = 0.001d0,      & ! minimum dry thickness depth (m)
                        min_layer = 0.01d0,       & ! minimum thickness of the second rooting layer (m)
                      soil_roughl = 0.05d0,       & ! soil roughness length (m)
                   top_soil_depth = 0.3d0,        & ! depth to which we conider the top soil to extend (m)
                         min_root = 5d0,          & ! minimum root biomass (gBiomass.m-2)
                  min_throughfall = 0.2d0,        & ! minimum fraction of precipitation which
                                                    ! is through fall
                     min_storage = 0.2d0            ! minimum canopy water (surface) storage (mm)

! timing parameters
double precision, parameter :: &
                  seconds_per_day = 86400d0,      & ! number of seconds per day
                seconds_per_day_1 = 1d0/seconds_per_day ! inverse of seconds per day

!!!!!!!!!
! Module level variables
!!!!!!!!!

! hydraulic model variables
integer :: water_retention_pass, soil_layer
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand, &
                                                layer_thickness    ! thickness of soil layers (m)
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
                      cond1, cond2, cond3, potA, potB    ! Saxton equation values

double precision :: root_reach, root_biomass, &
                  drythick, & ! estimate of the thickness of the dry layer at soil surface (m)
                    soilRT, &
                      wSWP, & ! weighted soil water potential (MPa) used in GSI calculate.
                              ! Removes / limits the fact that very low root density in young plants
                              ! give values too large for GSI to handle.
                 max_depth, & ! maximum possible root depth (m)
                    root_k, & ! biomass to reach half max_depth
   liquid,drainlayer,unsat, & ! variables used in drainage (m)
                    runoff, & ! runoff (kgH2O.m-2.day-1)
                 underflow, & ! drainage from the bottom of soil column (kgH2O.m-2.day-1)
  new_depth,previous_depth, & ! depth of bottom of soil profile
               canopy_wind, & ! wind speed (m.s-1) at canopy top
                     ustar, & ! friction velocity (m.s-1)
            air_density_kg, & ! air density kg/m3
                    roughl, & ! roughness length (m)
              displacement, & ! zero plane displacement (m)
                max_supply, & ! maximum water supply (mmolH2O/m2/day)
                     meant, & ! mean air temperature (oC)
                   meant_K, & ! mean air temperature (K)
        canopy_swrad_MJday, & ! canopy_absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_par_MJday, & ! canopy_absorbed PAR radiation (MJ.m-2.day-1)
          soil_swrad_MJday, & ! soil absorbed shortwave radiation (MJ.m-2.day-1)
          canopy_lwrad_Wm2, & ! canopy absorbed longwave radiation (W.m-2)
            soil_lwrad_Wm2, & ! soil absorbed longwave radiation (W.m-2)
      stomatal_conductance, & ! maximum stomatal conductance (mmolH2O.m-2.day-1)
   aerodynamic_conductance, & ! bulk surface layer conductance (m.s-1)
          soil_conductance, & ! soil surface conductance (m.s-1)
         convert_ms1_mol_1, & ! Conversion ratio for m.s-1 -> mol.m-2.s-1
                    lambda, & ! latent heat of vapourisa/tion (J.kg-1)
                     psych, & ! psychrometric constant (kPa K-1)
                     slope, & ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
            canopy_storage, & ! water storage on canopy (kg.m-2)
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
 max_lai_lwrad_absorption, & ! Max fraction of LW from sky absorbed by canopy
lai_half_lwrad_absorption, & ! LAI at which canopy LW absorption = 50 %
   max_lai_nir_absorption, & ! Max fraction of NIR absorbed by canopy
  lai_half_nir_absorption, & ! LAI at which canopy NIR absorption = 50 %
   max_lai_par_absorption, & ! Max fraction of PAR absorbed by canopy
  lai_half_par_absorption, & ! LAI at which canopy PAR absorption = 50 %
  max_lai_swrad_reflected, & ! Max fraction of SW rad reflected back to sky
 lai_half_swrad_reflected, & ! LAI at which SW reflection = 50 %
    lai_half_lwrad_to_sky, & ! LAI at which 50 % LW is reflected back to sky
    soil_swrad_absorption, & ! Fraction of SW rad absorbed by soil
    max_lai_lwrad_release, & ! Max fraction of LW emitted from canopy to be
   lai_half_lwrad_release, & ! LAI at which LW emitted from canopy to be released at 50 %
   soilevap_rad_intercept, & ! Intercept (kgH2O/m2/day) on linear adjustment to soil evaporation
                             ! to account for non-calculation of  energy balance
        soilevap_rad_coef    ! Coefficient on linear adjustment to soil evaporation to account for
                             ! non-calculation of energy balance


! Module level variables for step specific met drivers
double precision :: mint, & ! minimum temperature (oC)
                    maxt, & ! maximum temperature (oC)
                   swrad, & ! incoming short wave radiation (MJ/m2/day)
                     co2, & ! CO2 (ppm)
                     doy, & ! Day of year
                rainfall, & ! rainfall (kgH2O/m2/s)
                wind_spd, & ! wind speed (m/s)
                  vpd_pa, & ! Vapour pressure deficit (Pa)
                     lai, & ! leaf area index (m2/m2)
                  root_C    ! fine root stock (gC/m2)

! Module level varoables for step specific timing information
double precision :: seconds_per_step, & !
                       days_per_step, & !
                     days_per_step_1, & !
                        dayl_seconds, & ! day length in seconds
                          dayl_hours    ! day length in hours

double precision, dimension(:), allocatable ::    deltat_1, & ! inverse of decimal days
                                               soilwatermm, &
                                                 wSWP_time
double precision :: gpp_co2
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
    double precision ::  infi,tmp,tmp2, tmp3, tmp1 = 0d0 &
               ,ET_pot,ET_net &
                     ,deltaWP & ! deltaWP (MPa) minlwp-soilWP
                        ,Rtot   ! MPa.s-1.m-2.mmol-1

    integer :: p,f,nxp,n,test,m

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
    NUE                       = 1.850535d+01  ! Photosynthetic nitrogen use efficiency at optimum temperature (oC)
                                              ! ,unlimited by CO2, light and
                                              ! photoperiod
                                              ! (gC/gN/m2leaf/day)
    pn_max_temp               = 6.982614d+01  ! Maximum temperature for photosynthesis (oC)
    pn_opt_temp               = 3.798068d+01  ! Optimum temperature for photosynthesis (oC)
    pn_kurtosis               = 1.723531d-01  ! Kurtosis of photosynthesis temperature response
    e0                        = 4.489652d+00  ! Quantum yield gC/MJ/m2/day PAR
    max_lai_lwrad_absorption  = 9.282892d-01  ! Max fraction of LW from sky absorbed by canopy
    lai_half_lwrad_absorption = 5.941333d-01  ! LAI at which canopy LW absorption = 50 %
    max_lai_nir_absorption    = 8.333743d-01  ! Max fraction of NIR absorbed by canopy
    lai_half_nir_absorption   = 2.148633d+00  ! LAI at which canopy NIR absorption = 50 %
    minlwp                    = -1.990154d+00 ! minimum leaf water potential (MPa)
    max_lai_par_absorption    = 8.737539d-01  ! Max fraction of PAR absorbed by canopy
    lai_half_par_absorption   = 1.804925d+00  ! LAI at which canopy PAR absorption = 50 %
    lai_half_lwrad_to_sky     = 2.489314d+00  ! LAI at which 50 % LW is reflected back to sky
    iWUE                      = 1.722579d-02  ! Intrinsic water use efficiency (gC/m2leaf/day/mmolH2Ogs)
    soil_swrad_absorption     = 7.375071d-01  ! Fraction of SW rad absorbed by soil
    max_lai_swrad_reflected   = 2.796492d-01  ! Max fraction of SW reflected back to sky
    lai_half_swrad_reflected  = (lai_half_nir_absorption+lai_half_par_absorption) * 0.5d0
    max_lai_lwrad_release     = 2.481599d-01  ! Max fraction of LW emitted from canopy to be released
    lai_half_lwrad_release    = 5.020443d-01  ! LAI at which LW emitted from canopy to be released at 50 %
    soilevap_rad_intercept    = 1.122969d-02  ! Intercept (kgH2O/m2/day) on linear adjustment to soil evaporation
                                              ! to account for non-calculation
                                              ! of energy balance
    soilevap_rad_coef         = 1.748044d+00  ! Coefficient on linear adjustment to
                                              ! soil evaporation to account for
                                              ! non-calculation of energy
                                              ! balance
    ! check loaded parameters
    avN = pars(1)
    ! if input minlwp not empty over-write default value
    if (pars(2) /= -9999d0) minlwp = pars(2)
    deltaWP = minlwp
    root_k = pars(3) ; max_depth = pars(4)

    ! reset values
    intercepted_rainfall = dble_zero ; canopy_storage = dble_zero

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
    if (met(8,1) /= -9999d0) then
        meant = met(8,1)  ! mean air temperature (oC)
    else
        meant = (met(3,1)+met(2,1))*0.5d0
    endif
    meant_K = meant + freeze
    lai = met(11,1) ! leaf area index (m2/m2)
    root_C = met(12,1) ! fine root stock (gC/m2)

    ! initialise root reach based on initial conditions
    root_biomass = max(min_root,root_C*2d0)
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! determine initial soil layer thickness
    layer_thickness(1) = top_soil_depth ; layer_thickness(2) = max(min_layer,root_reach-layer_thickness(1))
    layer_thickness(3) = max_depth - sum(layer_thickness(1:2))
    previous_depth = max(top_soil_depth,root_reach)
    ! needed to initialise soils
    call calculate_Rtot(Rtot)
    ! used to initialise soils
    FLUXES(1,2) = calculate_update_soil_water(dble_zero,dble_zero,dble_zero) ! assume no evap or rainfall
    ! store soil water content of the rooting zone (mm)
    POOLS(1,1) = 1d3*sum(soil_waterfrac(1:nos_root_layers)*layer_thickness(1:nos_root_layers))

    !
    ! Begin looping through each time step
    !

    do n = start, finish

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! Incoming drivers
      mint = met(2,n)  ! minimum temperature (oC)
      maxt = met(3,n)  ! maximum temperature (oC)
      swrad = met(4,n) ! incoming short wave radiation (MJ/m2/day)
      co2 = met(5,n)   ! CO2 (ppm)
      doy = met(6,n)   ! Day of year
      rainfall = max(dble_zero,met(7,n)) ! rainfall (kgH2O/m2/s)
      if (met(8,n) /= -9999d0) then
          meant = met(8,n)  ! mean air temperature (oC)
      else
          meant = (met(3,n)+met(2,n))*0.5d0
      endif
      meant_K = meant + freeze
      wind_spd = met(9,n) ! wind speed (m/s)
      vpd_pa = met(10,n)  ! Vapour pressure deficit (Pa)
      lai = met(11,n)     ! leaf area index (m2/m2)
      root_C = met(12,n)  ! fine root stock (gC/m2)

      ! calculate daylength in hours and seconds
      dayl_hours = daylength_hours((doy-(deltat(n)*0.5d0)),lat)
      dayl_seconds = dayl_hours * 3600d0
      seconds_per_step = seconds_per_day * deltat(n)
      days_per_step = deltat(n)
      days_per_step_1 = deltat_1(n)

      !!!!!!!!!!
      ! calculate water balance
      !!!!!!!!!!

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      root_biomass = max(min_root,root_C*2d0)
      ! estimate drythick for the current step
      drythick = max(min_drythick, top_soil_depth * min(dble_one,dble_one - (soil_waterfrac(1) / porosity(1))))
      call calculate_Rtot(Rtot)
      ! pass Rtot to output variable and update deltaWP between minlwp and
      ! current weighted soil WP
      wSWP_time(n) = wSWP ; deltaWP = min(dble_zero,minlwp-wSWP)

      !!!!!!!!!!
      ! Calculate surface exchange coefficients
      !!!!!!!!!!

      ! calculate aerodynamic using consistent approach with SPA
      call calculate_aerodynamic_conductance

      ! calculate variables used commonly between ACM_GPP and ACM_ET
      call acm_albedo_gc(abs(deltaWP),Rtot)

      !!!!!!!!!!
      ! GPP (gC.m-2.day-1)
      !!!!!!!!!!

      if (stomatal_conductance > dble_zero) then
         FLUXES(n,1) = max(dble_zero,acm_gpp(stomatal_conductance))
      else
         FLUXES(n,1) = dble_zero
      endif

      ! Potential latent energy (kg.m-2.day-1)
      call acm_et(FLUXES(n,3),FLUXES(n,2),FLUXES(n,4))

      ! determine potential evaporation which sources its water from the soil,
      ! i.e. excluding wet canopy surface
      ET_pot = FLUXES(n,2) + FLUXES(n,4)
      ! do mass balance (i.e. is there enough water to support ET)
      ET_net = calculate_update_soil_water(FLUXES(n,2)*days_per_step,FLUXES(n,4)*days_per_step, &
                                          ((rainfall-intercepted_rainfall)*seconds_per_step))

      ! store soil water content of the rooting zone (mm)
      POOLS(n,1) = 1d3*sum(soil_waterfrac(1:nos_root_layers)*layer_thickness(1:nos_root_layers))

      ! pass soil surface runoff and underflow (drainage) from soil column kg/m2/day
      FLUXES(n,5) = runoff
      FLUXES(n,6) = underflow

      !!!!!!!!!!
      ! Bug checking
      !!!!!!!!!!

      do nxp = 1, nopools
         if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < 0d0 .or. &
             abs(POOLS(n+1,nxp)) == abs(log(infi))) then
            print*,"step",n,"POOL",nxp
            print*,"met",met(:,n)
            print*,"POOLS",POOLS(n,:)
            print*,"FLUXES",FLUXES(n,:)
            print*,"POOLS+1",POOLS(n+1,:)
            print*,"wSWP",wSWP
            print*,"waterfrac",soil_waterfrac
            stop
         endif
      enddo

      do nxp = 1, nofluxes
         if (FLUXES(n,nxp) /= FLUXES(n,nxp) .or. abs(FLUXES(n,nxp)) == abs(log(infi))) then
            print*,"step",n,"FLUX",nxp
            print*,"met",met(:,n)
            print*,"POOLS",POOLS(n,:)
            print*,"FLUXES",FLUXES(n,:)
            print*,"POOLS+1",POOLS(n+1,:)
            print*,"wSWP",wSWP
            print*,"waterfrac",soil_waterfrac
            stop
         endif
      enddo

    end do ! nodays loop

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm_gpp(gs)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
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
    co2_half_sat   = arrhenious(kc_saturation,kc_half_sat_conc,maxt)
    co2_comp_point = arrhenious(co2comp_saturation,co2comp_half_sat_conc,maxt)

    !
    ! Metabolic limited photosynthesis
    !

    ! maximum rate of temperature and nitrogen (canopy efficiency) limited
    ! photosynthesis (gC.m-2.day-1)
    pn = lai*avN*NUE*opt_max_scaling(pn_max_temp,pn_opt_temp,pn_kurtosis,maxt)

    !
    ! Diffusion limited photosynthesis
    !

    ! daily canopy conductance (mmolH2O.m-2.day-1-> molCO2.m-2.day-1)
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
    ! calculate CO2 limited rate of photosynthesis (gC.m-2.day-1)
    pd = (gc * (co2-ci)) * umol_to_gC
    ! scale to day light period as this is then consistent with the light
    ! capture period
    pd = pd * (dayl_hours / 24d0)

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
    double precision :: gs_high, gpp_high, gpp_low

    ! estimate photosynthesis with current estimate of gs
    gpp_low = acm_gpp(gs_in)

    ! Increment gs
    gs_high = gs_in + delta_gs
    ! estimate photosynthesis with incremented gs
    gpp_high = acm_gpp(gs_high)

    ! determine impact of gs increment on pd and how far we are from iWUE
    find_gs = iWUE - ((gpp_high - gpp_low)/lai)

  end function find_gs
  !
  !------------------------------------------------------------------
  !
  subroutine acm_et(wetcanopy_evap,transpiration,soilevap)

    ! Three response function(s) based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kgH20.m-2.day-1).
    ! Function 1 calculates the potential canopy level evaporation.
    ! Function 2 calculates the transpiration linking to hydraulic limitations.
    ! F1-F2 = potential canopy wet surface evaporation.
    ! Function 3 quantifies the potential soil surface evaporation flux

    implicit none

    ! arguments
    double precision, intent(out) :: wetcanopy_evap, transpiration, soilevap ! kgH2O.m-2.day-1

    ! local variables
    double precision :: &
              canopy_radiation, soil_radiation & ! isothermal net radiation (W/m2)
                              ,water_diffusion & ! Diffusion of water through soil matrix (m.s-1)
                                 ,water_supply & ! Potential water supply to canopy from soil (kgH2O.m-2.day-1)
                                    ,Jm3kPaK_1 &
                                        ,esurf & ! see code below
                                         ,esat & ! soil air space saturation vapour pressure
                                          ,gws & ! water vapour conductance through soil air space (m.s-1)
                                           ,ea & ! water vapour content of air (Pa)
                                        ,gs,gb   ! stomatal and boundary layer conductance (m.s-1)

    !!!!!!!!!!
    ! Estimate energy radiation balance (W.m-2)
    !!!!!!!!!!

    ! Absorbed shortwave radiation MJ.m-2.day-1 -> J.m-2.s-1
    canopy_radiation = canopy_lwrad_Wm2 + (canopy_swrad_MJday * 1d6 * seconds_per_day_1)
    soil_radiation = soil_lwrad_Wm2 + (soil_swrad_MJday * 1d6 * seconds_per_day_1)

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! calculate potential water supply (kgH2O.m-2.day-1)
    ! provided potential upper bound on evaporation
    water_supply = max_supply * mmol_to_kg_water

    ! Change units of potential stomatal conductance
    ! (m.s-1 <-> mmolH2O.m-2.day-1).
    ! Note assumption of sea surface pressure only
    gs = stomatal_conductance / (convert_ms1_mol_1 * 1d3)
    ! Combine in series stomatal conductance with boundary layer
    gb = aerodynamic_conductance

    !!!!!!!!!!
    ! Calculate evaporative fluxes (W.m-2)
    !!!!!!!!!!

    ! Calculate energy change due to VPD per kelvin (multiple use value)
    Jm3kPaK_1 = air_density_kg*cpair*vpd_pa*1d-3
    ! Calculate numerator of Penman Montheith (kg.m-2.day-1)
    wetcanopy_evap = (slope*canopy_radiation) + (Jm3kPaK_1*gb)
    ! Calculate the transpiration flux and restrict by potential water supply
    ! over the day
    transpiration = min(water_supply,(wetcanopy_evap / (lambda*(slope+(psych*(dble_one+gb/gs)))))*seconds_per_day)
    transpiration = max(dble_zero,transpiration)
    ! Calculate the potential wet canopy evaporation, limited by energy used for
    ! transpiration
    wetcanopy_evap = (wetcanopy_evap / (lambda*(slope+psych))) * seconds_per_day
    wetcanopy_evap = wetcanopy_evap - transpiration

    ! Update based on canopy water storage
    call canopy_interception_and_storage(wetcanopy_evap)

    ! Estimate water diffusion rate (m2.s-1) Jones (2014) appendix 2
    water_diffusion = 24.2d-6 * ( (maxt+freeze) / 293.2d0 )**1.75d0
    ! Soil conductance to water vapour diffusion (m s-1)...
    gws = porosity(1) * water_diffusion / (tortuosity*drythick)
    ! apply potential flow restriction at this stage
    gws = min(gws,(soil_waterfrac(1)*top_soil_depth*1d3)/dayl_seconds)

    ! calculate saturated vapour pressure (kPa), function of temperature.
    esat = 0.1d0 * exp( 1.80956664d0 + ( 17.2693882d0 * (maxt+freeze) - 4717.306081d0 ) / ( maxt+freeze - 35.86d0 ) )
    ! vapour pressure of air (kPa)
    ea = esat - (vpd_pa * 1d-3)
    ! vapour pressure in soil airspace (kPa), dependent on soil water potential
    ! - Jones p.110. partial_molar_vol_water
    esurf = esat * exp( 1d6 * SWP(1) * partial_molar_vol_water / ( Rcon * (maxt+freeze) ) )
    ! calculate VPD of the soil surface (kPa)
    esurf = esat - esurf
    ! now difference in VPD between soil air space and canopy
    esurf = (vpd_pa * 1d-3) - esurf

    ! update soil isothermal net radiation to net radiation
!    soil_radiation = soil_radiation + calculate_soil_netrad_adjustment(soil_radiation,gs,gb,esurf)
    ! Estimate potential soil evaporation flux (kg.m-2.day-1)
    soilevap = (slope*soil_radiation) + (air_density_kg*cpair*esurf*soil_conductance)
    soilevap = soilevap / (lambda*(slope+(psych*(dble_one+soil_conductance/gws))))
    soilevap = soilevap * seconds_per_day
    ! apply statistical adjustment to account for energy balance
    soilevap = soilevap_rad_intercept + (soilevap * soilevap_rad_coef)
    soilevap = max(dble_zero,min(soil_waterfrac(1) * top_soil_depth * 1d3, soilevap))

  end subroutine acm_et
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
    double precision :: ustar_Uh

    ! calculate the zero plane displacement and roughness length
    call z0_displacement(ustar_Uh)
    ! calculate friction velocity at tower height (reference height ) (m.s-1)
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
!    ustar = (wind_spd / log((tower_height-displacement)/roughl)) * vonkarman
    ustar = wind_spd * ustar_Uh

    ! based on Harman & Finnigan (2008); neutral conditions only
    call log_law_decay

    ! calculate soil surface conductance
    call calculate_soil_conductance(ustar_Uh)
    ! calculate bulk conductance (Jones p68)
    aerodynamic_conductance = (canopy_wind * vonkarman_2) &
                            / (log((canopy_height-displacement)/roughl))**2

  end subroutine calculate_aerodynamic_conductance
  !
  !------------------------------------------------------------------
  !
  subroutine log_law_decay

    ! Standard log-law above canopy wind speed (m.s-1) decay under neutral
    ! conditions.
    ! See Harman & Finnigan 2008; Jones 1992 etc for details.

    implicit none

    ! local parameters
    double precision, parameter :: min_wind = 0.01d0 ! minimum wind speed at canopy top

    ! log law decay
    canopy_wind = (ustar / vonkarman) * log(displacement / roughl)

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
                                   water_retention_saxton_eqns , x1 , x2 , xacc )
    enddo

  end subroutine calculate_field_capacity
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
    double precision :: s, denom, mult, tmp
    double precision, parameter :: max_gs = 5000d0, & ! mmolH2O.m-2.s-1
                                   min_gs = 5d-12

    !!!!!!!!!!
    ! determine some multiple use constants
    !!!!!!!!!

    ! Density of air (kg/m3)
    air_density_kg = 353d0/(maxt+freeze)
    ! Conversion ratio for m.s-1 -> mol.m-2.s-1
    convert_ms1_mol_1 = const_sfc_pressure / ((maxt+freeze)*Rcon)
    ! latent heat of vapourisation,
    ! function of air temperature (J.kg-1)
    if (maxt < dble_one) then
        lambda = 2.835d6
    else
        lambda = 2501000d0-2364d0*maxt
    endif
    ! psychrometric constant (kPa K-1)
    psych = (0.0646d0*exp(0.00097d0*maxt))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult = maxt+237.3d0
    ! 2502.935945 = 0.61078*17.269*237.3
    s = 2502.935945d0*exp(17.269d0*maxt/mult)
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = s/(mult*mult)

    !!!!!!!!!!
    ! Determine net shortwave and isothermal longwave energy balance
    !!!!!!!!!!

    call calculate_shortwave_balance
    call calculate_longwave_isothermal

    !!!!!!!!!!
    ! Calculate stomatal conductance under H2O and CO2 limitations
    !!!!!!!!!!

    if (deltaWP > vsmall) then
       ! Determine potential water flow rate (mmolH2O.m-2.dayl-1)
       max_supply = (deltaWP/Rtot) * dayl_seconds
    else
       ! set minimum (computer) precision level flow
       max_supply = vsmall
    end if

    ! Invert Penman-Monteith equation to give gs (m.s-1) needed to meet
    ! maximum possible evaporation for the day.
    ! This will then be reduced based on CO2 limits for diffusion based
    ! photosynthesis
    denom = slope * ((canopy_swrad_MJday * 1d6 * seconds_per_day_1) + canopy_lwrad_Wm2) &
            + (air_density_kg*cpair*vpd_pa*1d-3*aerodynamic_conductance)
    denom = (denom / (lambda * max_supply * mmol_to_kg_water * seconds_per_day_1)) - slope
    denom = denom / psych
    stomatal_conductance = aerodynamic_conductance / denom

    ! convert m.s-1 to mmolH2O.m-2.s-1
    stomatal_conductance = stomatal_conductance * 1d3 * convert_ms1_mol_1
    if (stomatal_conductance < dble_zero .or. stomatal_conductance > max_gs) stomatal_conductance = max_gs

    ! solve for photosynthesis limits on gs through iterative solution
    delta_gs = 1d-3*dayl_seconds ! mmolH2O/m2/dayl
    stomatal_conductance = zbrent('acm_albedo_gc:find_gs',find_gs,min_gs,stomatal_conductance,delta_gs)

  end subroutine acm_albedo_gc
  !
  !------------------------------------------------------------------
  !
  subroutine calculate_longwave_isothermal

    ! Subroutine estimates the isothermal net longwave radiation (W.m-2) for
    ! the canopy and soil surface. SPA uses a complex multi-layer radiative
    ! transfer scheme including reflectance, transmittance any absorption.
    ! However, for a given canopy vertical profiles, the LAI absorption
    ! relationship is readily predicted via Michaelis-Menten or
    ! non-rectangular hyperbola as done here.

    implicit none

    ! local variables
    double precision :: lwrad, & ! downward long wave radiation from sky (W.m-2)
             longwave_release, & ! emission of long wave radiation from surfaces per m2
                                 ! assuming isothermal condition (W.m-2)
     canopy_absorbed_fraction, & ! fraction of incoming longwave absorbed by canopy
        sky_returned_fraction, & ! fraction of incoming longwave reflected back into sky
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
    longwave_release = emiss_boltz * (maxt+freeze) ** 4

    !!!!!!!!!!
    ! Determine fraction of longwave absorbed by canopy and returned to the sky
    !!!!!!!!!!

    ! Calculate potential canopy absorbed fraction of longwave radiation
    ! incoming from outside of the canopy
    ! as a function of LAI
    canopy_absorbed_fraction = max_lai_lwrad_absorption*lai/(lai+lai_half_lwrad_absorption)
    ! calculate the fraction of longwave radiation from sky which is reflected
    ! back into the sky
    sky_returned_fraction = (dble_one-emissivity) - &
                             (((dble_one-emissivity)*lai) / (lai+lai_half_lwrad_to_sky))
    ! Calculate the potential absorption of longwave radiation lost from the
    ! canopy to soil / sky
    canopy_release_fraction = max_lai_lwrad_release - &
                             ((max_lai_lwrad_release*lai) / (lai+lai_half_lwrad_release))

    !!!!!!!!!!
    ! Calculate fluxes for canopy long wave absorption
    !!!!!!!!!!

    ! long wave absorbed by the canopy from the sky
    canopy_absorption_from_sky = lwrad * canopy_absorbed_fraction
    ! long wave absorbed by the canopy from soil
    canopy_absorption_from_soil = longwave_release * canopy_absorbed_fraction
    ! calculate two-sided long wave radiation emitted from canopy which is
    ! ultimately lost from to soil or sky (i.e. this value is used twice, once
    ! to soil once to sky)
    canopy_loss = 2d0 * longwave_release * lai * canopy_release_fraction

    !!!!!!!!!!
    ! Calculate fluxes for soil long wave absorption
    !!!!!!!!!!

    ! Long wave absorbed by soil from the sky.
    ! Accounting for fraction absorbed by the canopy and returned to the sky.
    ! We assume that long wave absorption is equivalent the the emissivity of
    ! the surface
    soil_absorption_from_sky = dble_one - (canopy_absorbed_fraction + sky_returned_fraction)
    soil_absorption_from_sky = soil_absorption_from_sky * lwrad * emissivity
    ! Calculate longwave absorbed by soil which si released by the canopy itself
    soil_absorption_from_canopy = canopy_loss * emissivity

    !!!!!!!!!!
    ! Isothermal net long wave canopy and soil balance (W.m-2)
    !!!!!!!!!!

    ! determine isothermal net canopy. Note two canopy_loss used to account for
    ! upwards and downwards emissions
    canopy_lwrad_Wm2 = (canopy_absorption_from_sky + canopy_absorption_from_soil) - (canopy_loss + canopy_loss)
    ! determine isothermal net soil
    soil_lwrad_Wm2 = (soil_absorption_from_sky + soil_absorption_from_canopy) - longwave_release

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

    ! parameters
    double precision, parameter :: sw_diffuse_fraction = 0.5d0

    ! local variables
    double precision :: absorbed_nir_fraction & !
                       ,absorbed_par_fraction & !
                       ,sky_absorped_fraction   !

    !!!!!!!!!!
    ! Determine canopy absorption / reflectance as function of LAI
    !!!!!!!!!!

    ! Canopy absorption of near infrared and photosynthetically active radiation
    absorbed_nir_fraction = max_lai_nir_absorption*lai/(lai+lai_half_nir_absorption)
    absorbed_par_fraction = max_lai_par_absorption*lai/(lai+lai_half_par_absorption)
    ! Canopy reflectance of SW radiation back into the sky
    sky_absorped_fraction = max_lai_swrad_reflected*lai/(lai+lai_half_swrad_reflected)

    !!!!!!!!!!
    ! Determine SW balance
    !!!!!!!!!!

    ! now determine shortwave radiation absorbed by the canopy (MJ.m-2.day-1)
    canopy_par_MJday = (sw_diffuse_fraction * swrad * absorbed_par_fraction)
    canopy_swrad_MJday = canopy_par_MJday + (sw_diffuse_fraction * swrad * absorbed_nir_fraction)
    ! absorption of shortwave radiation by soil (MJ.m-2.day-1)
    soil_swrad_MJday = dble_one - (sw_diffuse_fraction*absorbed_nir_fraction + &
                                   sw_diffuse_fraction*absorbed_par_fraction + sky_absorped_fraction)
    soil_swrad_MJday = swrad * soil_swrad_MJday * soil_swrad_absorption

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
    double precision :: slpa, cumdepth, prev, curr, sum_water_flux, &
                        soilR1,soilR2,transpiration_resistance,root_reach_local
    double precision, dimension(nos_root_layers) :: root_mass    &
                                                   ,soilRT_local &
                                                   ,root_length  &
                                                   ,ratio

    ! reset water flux
    water_flux = dble_zero ; wSWP = dble_zero ; soilRT_local = dble_zero ; soilRT = dble_zero
    ratio = dble_zero ; ratio(1) = dble_one
    ! calculate soil depth to which roots reach
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! calculate the plant hydraulic resistance component. Currently unclear
    ! whether this actually varies with height or whether tall trees have a
    ! xylem architecture which keeps the whole plant conductance (gplant) 1-10 (ish).
!    transpiration_resistance = (gplant * lai)**(-dble_one)
    transpiration_resistance = canopy_height / (gplant * lai)

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 2 root layer system this can be estimates more simply

    ! top 25 % of root profile
    slpa = (root_reach * 0.25d0) - layer_thickness(1)
    if (slpa <= dble_zero) then
        ! > 50 % of root is in top layer
        root_mass(1) = root_biomass * 0.5d0
        root_mass(1) = root_mass(1) + ((root_biomass-root_mass(1)) * (abs(slpa)/root_reach))
    else
        ! < 50 % of root is in bottom layer
        root_mass(1) = root_biomass * 0.5d0 * (layer_thickness(1)/(abs(slpa)+layer_thickness(1)))
    endif
    root_mass(2) = max(dble_zero,root_biomass - root_mass(1))
    root_length = root_mass * root_mass_length_coef_1
!    root_length = root_mass / (root_density * root_cross_sec_area)

    ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
    root_reach_local = min(root_reach,layer_thickness(1))
    soilR1=soil_resistance(root_length(1),root_reach_local,soil_conductivity(1)*head_1)
    soilR2=root_resistance(root_mass(1),root_reach_local)
    soilRT_local(1)=soilR1 + soilR2 + transpiration_resistance
    ! calculate and accumulate steady state water flux in mmol.m-2.s-1
    ! NOTE: Depth correction already accounted for in soil resistance
    ! calculations and this is the maximum potential rate of transpiration
    ! assuming saturated soil and leaves at their minimum water potential.
    ! also note that the head correction is now added rather than
    ! subtracted in SPA equations because deltaWP is soilWP-minlwp not
    ! soilWP prior to application of minlwp
    demand = abs(minlwp-SWP(1:nos_root_layers))+head*canopy_height
    water_flux(1) = demand(1)/(transpiration_resistance + soilR1 + soilR2)
    ! Bottom root layer
    if (root_mass(2) > dble_zero ) then
       ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
       soilR1=soil_resistance(root_length(2),layer_thickness(2),soil_conductivity(2)*head_1)
       soilR2=root_resistance(root_mass(2),layer_thickness(2))
       soilRT_local(2)=soilR1 + soilR2 + transpiration_resistance
       ! calculate and accumulate steady state water flux in mmol.m-2.s-1
!       demand = abs(minlwp-SWP(2))+head*canopy_height
       water_flux(2) = demand(2)/(transpiration_resistance + soilR1 + soilR2)
       ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))
    endif ! roots present in second layer?

    ! if freezing then assume soil surface is frozen
    if (meant < dble_one) then
        water_flux(1) = dble_zero ; ratio(1) = dble_zero ; ratio(2) = dble_one
    endif
    ! calculate sum value
    sum_water_flux = sum(water_flux)

    ! calculate weighted SWP and uptake fraction
    wSWP = sum(SWP(1:nos_root_layers) * water_flux(1:nos_root_layers))
    soilRT = sum(soilRT_local(1:nos_root_layers) * water_flux(1:nos_root_layers))
    uptake_fraction(1:nos_root_layers) = water_flux(1:nos_root_layers) / sum_water_flux
    wSWP = wSWP / sum_water_flux
    soilRT = soilRT / sum_water_flux

    ! sanity check in case of zero flux
    if (sum_water_flux == dble_zero) then
        wSWP = -20d0 ; soilRT = sum(soilRT_local)*0.5d0
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
  subroutine canopy_interception_and_storage(potential_evaporation)

    ! Simple daily time step integration of canopy rainfall interception, runoff
    ! and rainfall (kgH2O.m-2.s-1). NOTE: it is possible for intercepted rainfall to be
    ! negative if stored water running off into the soil is greater than
    ! rainfall (i.e. when leaves have died between steps)

    implicit none

    ! arguments
    double precision, intent(inout) :: potential_evaporation  ! wet canopy evaporation (kg.m-2.day-1),
                                                              ! enters as potential but leaves as water balance adjusted
    ! local variables
    integer :: i
    double precision :: tmp, through_fall, max_storage, max_storage_1 &
                       ,daily_addition, wetcanopy_evaporation

    ! determine maximum canopy storage & through fall fraction
    through_fall=max(min_throughfall,exp(-0.5d0*lai))
    ! maximum canopy storage (mm); minimum is applied to prevent errors in
    ! drainage calculation. Assume minimum capacity due to wood
    max_storage=max(min_storage,0.1d0*lai) ; max_storage_1 = max_storage**(-dble_one)
    ! potential intercepted rainfall (kgH2O.m-2.s-1)
    intercepted_rainfall = rainfall * (dble_one - through_fall)
    ! average rainfall intercepted by canopy (kgH2O.m-2.day-1)
    daily_addition = intercepted_rainfall * seconds_per_day

    tmp = dble_zero ; through_fall = dble_zero ; wetcanopy_evaporation = dble_zero
    ! if we have rainfall or canopy water currently stored
    if (intercepted_rainfall > dble_zero .or. canopy_storage > dble_zero) then

        ! intergrate over canopy for each day
        do i = 1, int(days_per_step)

           ! add new rain to the canopy
           canopy_storage = canopy_storage + daily_addition

           ! how much is over and above the max_storage capacity?
           tmp = max(dble_zero, canopy_storage - max_storage)
           ! add this back to the through fall
           through_fall = through_fall + tmp
           ! remove the difference from the canopy
           canopy_storage = canopy_storage - tmp

           ! scale potential evaporation by ratio of current to max storage
           tmp = min(canopy_storage,potential_evaporation * min(dble_one,canopy_storage * max_storage_1))
           ! add to the running total of wet canopy evaporation
           wetcanopy_evaporation = wetcanopy_evaporation + tmp
           ! now remove from canopy
           canopy_storage = canopy_storage - tmp

           ! in case of dew formation do overflow calculation again
           ! how much is over and above the max_storage capacity?
           tmp = max(dble_zero, canopy_storage - max_storage)
           ! add this back to the through fall
           through_fall = through_fall + tmp
           ! remove the difference from the canopy
           canopy_storage = canopy_storage - tmp

        end do ! day looping

        ! sanity checks
        ! NOTE: addition of 1e-10 is to prevent precision error causing stop
        if (canopy_storage > (max_storage + vsmall) .or. canopy_storage < dble_zero) then
           print*,"Canopy water storage mass balance error!!"
           print*,"through_fall_total",through_fall
           print*,"canopy_storage",canopy_storage,"max_storage",max_storage
           print*,"potential_evaporation",potential_evaporation,"actual",wetcanopy_evaporation * days_per_step_1
           stop
        endif

        ! average evaporative flux to daily rate (kgH2O/m2/day)
        potential_evaporation = wetcanopy_evaporation * days_per_step_1
        ! correct intercepted rainfall rate to kgH2O.m-2.s-1
        intercepted_rainfall = intercepted_rainfall - ((through_fall * days_per_step_1) * seconds_per_day_1)

        ! sanity check
        if (intercepted_rainfall > rainfall) then
            print*,"Canopy intercepted rainfall cannot be greater than rainfall!!"
            print*,"rainfall", rainfall, "through_fall", (through_fall * days_per_step_1 * seconds_per_day_1)
        endif

    end if ! we have some rainfall

  end subroutine canopy_interception_and_storage
  !
  !-----------------------------------------------------------------
  !
  subroutine infiltrate(rainfall_in)

    ! Takes surface_watermm and distrubutes it among top !
    ! layers. Assumes total infilatration in timestep.   !

    implicit none

    ! arguments
    double precision, intent(in) :: rainfall_in ! rainfall (kg.m-2.step-1)

    ! local argumemts
    integer :: i
    double precision    :: add   & ! surface water available for infiltration (m)
                          ,wdiff   ! available space in a given soil layer for water to fill (m)

    ! convert rainfall water from mm -> m (or kgH2O.m-2.step-1 -> MgH2O.m-2.step-1)
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

    end do

    ! if after all of this we have some water left assume it is runoff
    ! converted to kgH2O.m-2.day-1
    runoff = add * 1d3 * days_per_step_1

  end subroutine infiltrate
  !
  !-----------------------------------------------------------------
  !
  subroutine gravitational_drainage

    ! integrator for soil gravitational drainage !

    implicit none

    ! local variables..
    double precision  :: change, drainage, iceprop(nos_soil_layers)

    ! calculate soil ice proportion; at the moment
    ! assume everything liquid
    iceprop = dble_zero
    ! except the surface layer in the mean daily temperature is < 1oC
    if (meant < dble_one) iceprop(1) = dble_one

    do soil_layer = 1, nos_soil_layers

       ! liquid content of the soil layer, i.e. fraction avaiable for drainage
       liquid     = soil_waterfrac( soil_layer ) &
                  * ( dble_one - iceprop( soil_layer ) )
       ! soil water capacity of the current layer
       drainlayer = field_capacity( soil_layer )

       ! initial conditions; i.e. is there liquid water and more water than
       ! layer can hold
       if ( (liquid > dble_zero) .and. (soil_waterfrac( soil_layer ) > drainlayer) ) then

          ! unsaturated volume of layer below (m3 m-2)..
          unsat = max( dble_zero , ( porosity( soil_layer+1 ) - soil_waterfrac( soil_layer+1 ) ) &
                             * layer_thickness( soil_layer+1 ) / layer_thickness( soil_layer ) )

          ! potential drainage over time step
          drainage = soil_conductivity( soil_layer ) * seconds_per_step

!          ! gravitational drainage above field_capacity
!          ! already convered above
!          if ( soil_waterfrac(soil_layer) < drainlayer ) drainage = dble_zero
          ! ice does not drain
          if ( drainage > liquid ) drainage = liquid
          ! layer below cannot accept more water than unsat
          if ( drainage > unsat ) drainage = unsat
          ! water loss from this layer
          change = drainage * layer_thickness(soil_layer)
          ! update soil layer below with drained liquid
          watergain( soil_layer + 1 ) = watergain( soil_layer + 1 ) + change
          waterloss( soil_layer     ) = waterloss( soil_layer     ) + change

       end if ! some liquid water and drainage possible

    end do ! soil layers

    ! estimate drainage from bottom of soil column (kgH2O/m2/day)
    underflow = waterloss(nos_soil_layers) * 1d3 * days_per_step_1

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
    integer :: i
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
    ! seperately calculate the soil conductivity as this applies to each layer
    call calculate_soil_conductivity
    ! but apply the lowest soil layer to the core as well in initial conditions
    soil_conductivity(nos_soil_layers+1) = soil_conductivity(nos_soil_layers)

    ! final sanity check for porosity
    where (porosity <= field_capacity) porosity = field_capacity + 0.01d0

  end subroutine initialise_soils
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_soil_conductivity

    ! Calculate the soil conductivity (m s-1) of water based on soil
    ! characteristics and current water content

    implicit none

    ! soil conductivity for the dynamic soil layers (i.e. not including core)
    soil_conductivity(1:nos_soil_layers) = cond1(1:nos_soil_layers) &
                                        * exp(cond2(1:nos_soil_layers)+cond3(1:nos_soil_layers)/soil_waterfrac(1:nos_soil_layers))

    ! protection against floating point error
    where (soil_waterfrac < 0.05d0)
          soil_conductivity = 1d-30
    end where

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
                                   mult1 = 100d0, mult2 = 2.778d-6, mult3 = 1000d0

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
  subroutine calculate_soil_conductance(ustar_Uh)

    ! proceedsure to solve for soil surface resistance based on Monin-Obukov
    ! similarity theory stability correction momentum & heat are integrated
    ! through the under canopy space and canopy air space to the surface layer
    ! references are Nui & Yang 2004; Qin et al 2002
    ! NOTE: conversion to conductance at end

    implicit none

    ! declare arguments
    double precision, intent(in) :: ustar_Uh

    ! local variables
    double precision :: canopy_decay & ! canopy decay coefficient for soil exchange
                       ,lc,lm &
                       ,mult1,mult2,mult3,mult4 &
                       ,beta &         ! ustar/wind_spd ratio
                       ,Kh_canht       ! eddy diffusivity at canopy height (m2.s-1)

    ! parameters
    double precision, parameter :: foliage_drag = 0.2d0, & ! foliage drag coefficient
                                   beta_max = 1d0, beta_min = 0.2d0, &
                                   min_lai = 1d0, &
                                   most_soil = 1d0 ! Monin-Obukov similarity theory stability correction.
                                                   ! As no sensible heat flux
                                                   ! calculated,
                                                   ! assume neutral conditions
                                                   ! only

    ! ratio of friction velocity and canopy wind speed
    beta = ustar_Uh
!    beta = min(beta_max,max(ustar/canopy_wind,beta_min))
    ! both length scale and mixing length are considered to be constant within
    ! the canopy (under dense canopy conditions) calculate length scale (lc)
    ! for momentum absorption within the canopy; Harman & Finnigan (2007)
    ! and mixing length (lm) for vertical momentum within the canopy Harman & Finnigan (2008)
    if (lai > min_lai) then
        lc = (4d0*canopy_height) / lai
        lm = max(canopy_height*0.02d0, 2d0*(beta**3)*lc)
    else
        lc = vonkarman * tower_height
        lm = canopy_height * vonkarman
    endif

    ! calculate eddy diffusivity at the top of the canopy (m2.s-1)
    ! Kaimal & Finnigan 1994; for near canopy approximation
    Kh_canht=vonkarman*ustar*(canopy_height-displacement)

    ! calculate canopy decay coefficient with stability correction
    ! NOTE this is not consistent with canopy momentum decay done by Harman &
    ! Finnigan (2008)
    canopy_decay = (((foliage_drag*canopy_height*max(min_lai,lai))/lm)**0.5d0)*(most_soil**0.5d0)

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

    ! local variables..
    integer :: i

    ! reformulation aims to remove if statement within loop to hopefully improve
    ! optimisation
    SWP(1:nos_soil_layers) = -0.001d0 * potA(1:nos_soil_layers) &
                           * soil_waterfrac(1:nos_soil_layers)**potB(1:nos_soil_layers)
    where (SWP(1:nos_soil_layers) < -20d0) SWP(1:nos_soil_layers) = -20d0
!    where (soil_waterfrac(1:nos_soil_layers) < 0.005)
!        SWP(1:nos_soil_layers) = -9999.0
!    end where

  end subroutine soil_water_potential
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
                               min_lai = 1d0,   & ! Minimum LAI parameter as height does not vary with growth
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
  !------------------------------------------------------------------
  !
  double precision function daylength_hours(doy,lat)

    ! Function uses day of year and latitude (-90 / 90 degrees) as inputs,
    ! combined with trigonomic functions to calculate
    ! day length in hours

    implicit none

    ! arguments
    double precision, intent(in) :: doy, lat

    ! local variables
    double precision :: dec, mult, sinld, cosld, aob

!    dec = - asin( sin( 23.45 * deg_to_rad ) * cos( 2.0 * pi * ( doy + 10.0 ) / 365.0 ) )
    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10d0 ) / 365d0 ) )
    mult = lat * deg_to_rad
    sinld = sin( mult ) * sin( dec )
    cosld = cos( mult ) * cos( dec )
    aob = max(-dble_one,min(dble_one,sinld / cosld))

    ! define output
    daylength_hours = 12d0 * ( dble_one + 2d0 * asin( aob ) * pi_1 )

    ! return to user
    return

  end function daylength_hours
  !
  !------------------------------------------------------------------
  !
  double precision function calculate_update_soil_water(ET_leaf,ET_soil,rainfall_in)

   !
   ! Function updates the soil water status and layer thickness
   ! Soil water profile is updated in turn with evaporative losses,
   ! rainfall infiltration and gravitational drainage
   ! Root layer thickness is updated based on changes in the rooting depth from
   ! the previous step
   !

   implicit none

   ! arguments
   double precision, intent(in) :: ET_leaf,ET_soil & ! evapotranspiration estimate (kg.m-2.step-1)
                                      ,rainfall_in   ! rainfall (kg.m-2.step-1)

   ! local variables
   double precision ::  depth_change, water_change, tmp
   double precision, dimension(nos_root_layers) :: avail_flux, evaporation_losses

   ! seperately calculate the soil conductivity as this applies to each layer
   call calculate_soil_conductivity

   !!!!!!!!!!
   ! Evaporative losses
   !!!!!!!!!!

   ! Assume leaf transpiration is drawn from the soil based on the
   ! update_fraction estimated in calculate_Rtot
   evaporation_losses = ET_leaf * uptake_fraction
   ! Assume all soil evaporation comes from the soil surface only
   evaporation_losses(1) = evaporation_losses(1) + ET_soil
   ! can not evaporate from soil more than is available
   avail_flux = soil_waterfrac(1:nos_root_layers) * layer_thickness(1:nos_root_layers)
   where (evaporation_losses > avail_flux) evaporation_losses = avail_flux * 0.99d0

   ! this will update the ET estimate outside of the function
   ! unit / time correction also occurs outside of this function
   calculate_update_soil_water = sum(evaporation_losses)

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

   ! update osil water status with drainage
   soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness(1:nos_soil_layers)) &
                                        + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                     / layer_thickness(1:nos_soil_layers)
   ! reset soil water flux variables
   waterloss = dble_zero ; watergain = dble_zero

   !!!!!!!!!!
   ! Rainfal infiltration drainage
   !!!!!!!!!!

   ! Reset runoff variable before assignment 
   runoff = dble_zero

   ! determine infiltration from rainfall,
   ! if rainfall is probably liquid / soil surface is probably not frozen
   if (meant >= dble_zero .and. rainfall_in > dble_zero) then
       call infiltrate(rainfall_in)
   else
       runoff = rainfall_in * days_per_step_1
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

   !!!!!!!!!!
   ! Update soil layer thickness
   !!!!!!!!!!

   depth_change = dble_zero ; water_change = dble_zero
   ! if roots extent down into the bucket
   if (root_reach > top_soil_depth .or. previous_depth > top_soil_depth) then
      ! how much has root depth extended since last step?
      depth_change = root_reach - previous_depth

      ! if there has been an increase
      if (depth_change > dble_zero .and. root_reach > layer_thickness(1)+min_layer) then

         ! determine how much water is within the new volume of soil
         water_change = soil_waterfrac(nos_soil_layers) * depth_change
         ! now assign that new volume of water to the deep rooting layer
         soil_waterfrac(nos_root_layers) = ((soil_waterfrac(nos_root_layers) * layer_thickness(nos_root_layers)) &
                                            + water_change) / (layer_thickness(nos_root_layers)+depth_change)
         ! explicitly update the soil profile if there has been rooting depth
         ! changes
         layer_thickness(1) = top_soil_depth ; layer_thickness(2) = max(min_layer,root_reach-layer_thickness(1))
         layer_thickness(3) = max_depth - sum(layer_thickness(1:2))

      elseif (depth_change < dble_zero .and. root_reach > layer_thickness(1)+min_layer) then

         ! determine how much water is lost from the old volume of soil
         water_change = soil_waterfrac(nos_root_layers) * abs(depth_change)
         ! now assign that new volume of water to the deep rooting layer
         soil_waterfrac(nos_soil_layers) = ((soil_waterfrac(nos_soil_layers) * layer_thickness(nos_soil_layers)) &
                                            + water_change) / (layer_thickness(nos_soil_layers)+abs(depth_change))

         ! explicitly update the soil profile if there has been rooting depth
         ! changes
         layer_thickness(1) = top_soil_depth ; layer_thickness(2) = max(min_layer,root_reach-layer_thickness(1))
         layer_thickness(3) = max_depth - sum(layer_thickness(1:2))

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

  end function calculate_update_soil_water
  !
  !------------------------------------------------------------------
  !
  double precision function calculate_soil_netrad_adjustment(soil_radiation,gs,gb,vpd_in)

    ! function estimates the steady state solution to canopy temperature by
    ! balancing approximate turbulent fluxes

    implicit none

    ! arguments
    double precision, intent(in) ::  &
                           soil_radiation, & ! estimate of net absorbed radiation (SW+LW) J.m-2.day-1
                                   vpd_in, & ! vapour pressure deficit (Pa)
                                    gs,gb   ! canopy conductance (m.s-1)
    ! local variables
    double precision :: thermal_gains, thermal_losses, sm_1_kPaK_1,    &
                        soil_evaporative_resistance,                 &
                        soil_thermal_resistance,radiative_conductance, &
                        delta_temperature

    !!!!!!!!!!
    ! Calculate total thermal resistance in form of radiation (gr) and sensible
    ! heat. Units (m.s-1)
    !!!!!!!!!!

    ! calculate conductance rate for radiative heat exchange (m.s-1)
    radiative_conductance = (4d0 * emiss_boltz * (maxt+freeze) ** 3) / (air_density_kg * cpair)
    ! combine with conductance of sensible heat (aerodynamic_conductance) in
    ! parallel as these processes are indeed occuring in parallel.
    ! NOTE: that aerodynamic_conductance already implicitly contains LAI scaling
    ! from roughness length and displacement height components
    soil_thermal_resistance = radiative_conductance + gb * 0.93d0
    ! NOTE: conversion to resistance (s.m-1)
    soil_thermal_resistance = (soil_thermal_resistance) ** (-dble_one)

    !!!!!!!!!
    ! Calculate total evaporation resistance (s.m-1)
    !!!!!!!!!

    ! Canopy conductance (lai scaled anologue of stomatal conductance) is
    ! combined in series with aerodynamic conductance.
    ! NOTE: conversion of conductances to resistance
    soil_evaporative_resistance = gs ** (-dble_one) + gb ** (-dble_one)

    !!!!!!!!!
    ! Determine energy gains by canopy and then losses
    !!!!!!!!!

    ! calculate common denominator (units: sm-1 kPaK-1)
    sm_1_kPaK_1 = (psych * soil_evaporative_resistance) + (slope * soil_thermal_resistance)
    ! calculate thermal gains to the system (K)
    thermal_gains = soil_thermal_resistance * soil_evaporative_resistance * psych * soil_radiation
    thermal_gains = thermal_gains / (air_density_kg * cpair * sm_1_kPaK_1)
    ! calculate thermal losses to the system (K)
    thermal_losses = (soil_thermal_resistance * vpd_in * 1d-3) / sm_1_kPaK_1
    ! calculate net difference in canopy temperature
    delta_temperature = thermal_gains - thermal_losses

    calculate_soil_netrad_adjustment = - 4d0 * emiss_boltz * ( maxt+freeze ) ** 3 * ( delta_temperature )

    return

  end function calculate_soil_netrad_adjustment
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
    double precision,intent(in)             :: tol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
       double precision function func( xval )
         double precision ,intent(in) :: xval
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
end module CARBON_MODEl_MOD
