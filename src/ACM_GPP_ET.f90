
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
double precision, parameter :: xacc = 0.0001      & ! accuracy parameter for zbrent bisection proceedure ! 0.0001
                              ,dble_zero = 0.0    &
                              ,dble_one = 1.0     &
                              ,vsmall = tiny(0.0)

integer, parameter :: nos_root_layers = 2, nos_soil_layers = nos_root_layers + 1
double precision, parameter :: pi = 3.1415927,    &
                             pi_1 = pi**(-dble_one),   &
                              pi2 = pi**2,        &
                           two_pi = pi*2.0,       &
                       deg_to_rad = pi/180d0,     & 
              sin_dayl_deg_to_rad = sin( 23.45 * deg_to_rad ), & ! repeated function in acm
                          gravity = 9.8067,       & ! acceleration due to gravity, ms-1
                            boltz = 5.670400e-8,  & ! Boltzmann constant (W.m-2.K-4)
                       emissivity = 0.96,         &
                      emiss_boltz = emissivity * boltz, &
                           freeze = 273.15,       &
          partial_molar_vol_water = 18.05e-6,     & ! partial molar volume of water, m3 mol-1 at 20C
                       umol_to_gC = 1e-6*12,      & ! conversion of umolC -> gC
                 mmol_to_kg_water = 1.8e-5,       & ! milli mole conversion to kg
                   mol_to_g_water = 18.0,         & ! molecular mass of water
!snowscheme       density_of_water = 998.9,         & ! density of !water kg.m-3
                   gas_constant_d = 287.04,       & ! gas constant for dry air (J.K-1.mol-1) 
                             Rcon = 8.3144,       & ! Universal gas constant (J.K-1.mol-1)
                        vonkarman = 0.41,         & ! von Karman's constant
                      vonkarman_2 = 0.41**2,      & ! von Karman's constant^2
                            cpair = 1004.6          ! Specific heat capacity of air; used in energy balance J.kg-1.K-1

! photosynthesis / respiration parameters
double precision, parameter :: &
                    kc_saturation = 310.0,        & ! CO2 half saturation, saturation value
                 kc_half_sat_conc = 23.956,       & ! CO2 half sat, half sat
               co2comp_saturation = 36.5,         & ! CO2 compensation point, saturation
            co2comp_half_sat_conc = 9.46            ! CO2 comp point, half sat
                                                    ! Each of these are temperature
                                                    ! sensitivty
! hydraulic parameters
double precision, parameter :: &
                       tortuosity = 2.5,          & ! tortuosity
                           gplant = 5.0,          & ! plant hydraulic conductivity (mmol m-1 s-1 MPa-1) 
                      root_resist = 25.0,         & ! Root resistivity (MPa s g mmol−1 H2O)
                      root_radius = 0.0001,       & ! root radius (m) Bonen et al 2014 = 0.00029
                    root_radius_1 = root_radius**(-dble_one), &
              root_cross_sec_area = 3.141593e-08, & ! root cross sectional area (m2)
                                                    ! = pi * root_radius * root_radius 
                     root_density = 0.5e6,        & ! root density (g biomass m-3 root) 
                                                    ! 0.5e6 Williams et al 1996                                       
                                                    ! 0.31e6 Bonan et al 2014                  
          root_mass_length_coef_1 = (root_cross_sec_area * root_density)**(-1.d0), &
                             head = 0.009807,     & ! head of pressure (MPa/m)
                           head_1 = 101.968,      & ! inverse head of pressure (m/MPa)
                   minlwp_default = -2.060814       ! min leaf water potential (MPa)
! structural parameters
double precision, parameter :: &
                    canopy_height = 9.0,          & ! canopy height assumed to be 9 m
                     tower_height = canopy_height + 2.0, & ! tower (observation) height assumed to be 2 m above canopy
                     min_drythick = 0.001,        & ! minimum dry thickness depth (m)
                      soil_roughl = 0.05,         & ! soil roughness length (m)
                   top_soil_depth = 0.3,          & ! depth to which we conider the top soil to extend (m)
                         min_root = 5.0,          & ! minimum root biomass (gBiomass.m-2)
                  min_throughfall = 0.2,          & ! minimum fraction of precipitation which
                                                    ! is through fall 
                     min_storage = 0.2              ! minimum canopy water (surface) storage (mm)

! timing parameters
double precision, parameter :: &
                  seconds_per_day = 86400,        & ! number of seconds per day
                seconds_per_day_1 = 1d0/86400       ! inverse of seconds per day

!!!!!!!!!
! Module level variables
!!!!!!!!!

! hydraulic model variables 
integer :: water_retention_pass, soil_layer
double precision, dimension(nos_soil_layers) :: soil_frac_clay,soil_frac_sand, &
                                                layer_thickness    ! thickness of soil layers (m)
double precision, dimension(nos_root_layers) :: uptake_fraction, & ! 
                                                     water_flux    ! potential transpiration flux (mmol.m-2.s-1)
double precision, dimension(nos_soil_layers+1) :: SWP, & ! soil water potential (MPa) 
                                    soil_conductivity, & ! soil conductivity
                                            waterloss, & ! water loss from specific soil layers (m)
                                            watergain, & ! water gained by specfic soil layers (m)
                                       field_capacity, & ! soil field capacity (m3.m-3)
                                       soil_waterfrac, & ! soil water content (m3.m-3)
                                             porosity, & ! soil layer porosity, (fraction)
                      cond1, cond2, cond3, potA, potB    ! Saxton equation values

double precision :: root_reach, root_biomass,soil_depth, &
                  drythick, & ! estimate of the thickness of the dry layer at soil surface (m)
                    demand, & ! maximum potential canopy hydraulic demand
                    soilRT, &
                    minlwp, & ! minimum leaf water potential (MPa), either
                              ! parameter input or default value
                      wSWP, & ! weighted soil water potential (MPa) used in GSI calculate. 
                              ! Removes / limits the fact that very low root density in young plants
                              ! give values too large for GSI to handle.
                 max_depth, & ! maximum possible root depth (m)
                    root_k, & ! biomass to reach half max_depth
   liquid,drainlayer,unsat, & ! variables used in drainage (m)
                    runoff, & ! runoff (kg.m-2.step)
          seconds_per_step, & !
                        x1, & ! lower boundary condition for zbrent calculation
                        x2, & ! upper boundary condition for zbrent calculation
  new_depth,previous_depth, & ! depth of bottom of soil profile
   aerodynamic_conductance, & ! bulk surface layer conductance
                    roughl, & ! roughness length (m)
              displacement    ! zero plane displacement (m)

double precision, dimension(:), allocatable ::    deltat_1, & ! inverse of decimal days 
                                               soilwatermm, &
                                                 wSWP_time, &

save 

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
 
    ! The subroutine calls an Aggregated Canopy Model to simulate GPP (ACM_GPP)
    ! The subroutine calls an Aggregated Canopy Model to simulate Evapotranspiration (ACM_ET)

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
              ,rainfall_input &
               ,ET_pot,ET_net &
                       ,meant & ! mean air temperature (oC)
                     ,deltaWP & ! deltaWP (MPa) minlwp-soilWP
              ,canopy_storage & ! water storage on canopy (kg.m-2)
        ,intercepted_rainfall & ! intercepted rainfall rate equivalent (kg.m-2.s-1)
                 ,gpppars(12) & ! ACM inputs (LAI+met)
                 ,evappars(8) & ! ACM_ET parameters
               ,constants(10)   ! parameters for ACM

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

    ! check load parameters
    minlwp = minlwp_default
    if (pars(2) /= -9999) minlwp = pars(2)
    root_k = pars(3) ; max_depth = pars(4)

    ! load some values
    gpppars(4) = pars(1) ! foliar N gN/m2
    gpppars(7) = lat
    gpppars(9) = abs(minlwp) ! leafWP-soilWP (i.e. -2-0) ! p11 from ACM recal
                             ! NOTE: sign is forced positive for use in varous
                             ! equations
    gpppars(10) = dble_one ! totaly hydraulic resistance ! p12 from ACM recal (updated)
    gpppars(11) = pi

    ! assign acm parameters
    constants(1)=1.047076e+01  ! Nitrogen use efficiency (gC/gN per m2) ! 11.606093711
    constants(2)=9.768712e-03  ! Day length correction coefficient 
!    constants(3)=24.09383     ! co2 compensation point (ppm)
!    constants(4)=304.1392     ! co2 half saturation (ppm)
    constants(5)=9.534137e-02  ! Day length correction constant
    constants(6)=1.191852e-01  ! hydraulic adjustment coefficient on Rtot 
    constants(7)=7.396011e+00  ! maximum canopy quantum yield interception satuation (gC/MJ)
    constants(8)=3.426868e-02  ! Amax temperature response exponential coefficient
    constants(9)=2.567800e+00  ! LAI**2 at half saturation quantum yield (m2/m2)
    constants(10)=3.262981e-01  ! hydraulic adjustment exponent

    ! assign acm_et parameters
!    evappars(1) = 0.255475204  ! maximum fraction of radiation absorbed for canopy evaporation (fraction)
!    evappars(2) = 0.472094244  ! LAI**2 at half saturation radiation absorption (m2/m2)
!    evappars(3) = 0.038749412  ! maximum short wave radiation (W.m-2) input for which gc > 0
!    evappars(4) = 147.216593167! optimum short wave radiation (W.m-2) for gc openning
!    evappars(5) = 2.468762476  ! max gradient for gc response to hydraulic resistance (MPa)
!    evappars(6) = 14.598361448 ! half saturation for gc response to hydraulic resistance (MPa)
!    evappars(7) = 0.2496333    ! maxmimum fraction of radiation absorbed for soil evaporation (fraction)
!    evappars(8) = 0.05          ! kurtosis of radiation response on gc

    evappars(1) = 7.929720e-01
    evappars(2) = 2.474387e-01
    evappars(3) = 6.997828e+02
    evappars(4) = 4.042512e+02  
    evappars(5) = 3.802051e-03  
    evappars(6) = 3.349832e-02  
    evappars(7) = 0.7446446 !0.7782783 !pars(5)
    evappars(8) = 2.176627e-03

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
    ! initialise root reach based on initial conditions
    root_biomass = max(min_root,met(12,1)*2)
    root_reach = max_depth * root_biomass / (root_k + root_biomass)
    ! determine initial soil layer thickness
    layer_thickness(1) = top_soil_depth ; layer_thickness(2) = max(0.1,root_reach-layer_thickness(1))
    layer_thickness(3) = max_depth - sum(layer_thickness(1:2))
    previous_depth = max(top_soil_depth,root_reach)
    soil_depth = dble_zero ; previous_depth = dble_zero
    ! needed to initialise soils 
    call calculate_Rtot(gpppars(9),gpppars(10),met(11,1) &
                       ,deltat(1),((met(3,1)+met(2,1))*0.5))
    ! used to initialise soils
    FLUXES(1,2) = calculate_update_soil_water(dble_zero,dble_zero,((met(3,1)+met(2,1))*0.5))
    ! store soil water content of the rooting zone (mm)
    POOLS(1,1) = 1e3*sum(soil_waterfrac(1:nos_root_layers)*layer_thickness(1:nos_root_layers))

    ! 
    ! Begin looping through each time step
    ! 

    do n = start, finish

      !!!!!!!!!!
      ! assign drivers and update some prognostic variables
      !!!!!!!!!!

      ! timing variable
      seconds_per_step = seconds_per_day * deltat(n)
      meant = (met(3,n)+met(2,n))*0.5
      rainfall_input = max(dble_zero,met(7,n))

      ! load next met / lai values for ACM and acm_et
      gpppars(1)=met(11,n) ! lai
      gpppars(2)=met(3,n)  ! max temp
      gpppars(3)=met(2,n)  ! min temp
      gpppars(5)=met(5,n)  ! co2
      gpppars(6)=daylength_hours((met(6,n)-(deltat(n)*0.5)),gpppars(7)) 
      gpppars(8)=met(4,n)  ! radiation

      ! Temperature adjustments for Michaelis-Menten coefficients 
      ! for CO2 (kc) and O2 (ko) and CO2 compensation point
      constants(3) = arrhenious(co2comp_saturation,co2comp_half_sat_conc,gpppars(2))
      constants(4) = arrhenious(kc_saturation,kc_half_sat_conc,gpppars(2))

      !!!!!!!!!!
      ! calculate water balance
      !!!!!!!!!!

      ! calculate the minimum soil & root hydraulic resistance based on total
      ! fine root mass ! *2*2 => *RS*C->Bio
      root_biomass = max(min_root,met(12,n)*2)
      deltaWP = min(dble_zero,minlwp-wSWP)
      gpppars(9) = abs(deltaWP) ! update deltaWP
      ! estimate drythick for the current step
      drythick = max(min_drythick, top_soil_depth * min(dble_one,(soil_waterfrac(1) / field_capacity(1))))
      call calculate_Rtot(gpppars(9),gpppars(10),met(11,n) &
                         ,deltat(n),meant)
      ! pass Rtot to output variable and update deltaWP between minlwp and
      ! current weighted soil WP 
      wSWP_time(n) = wSWP 
      ! calculate aerodynamic resistance (1/conductance) using consistent
      ! approach with SPA
      call calculate_aerodynamic_conductance(met(11,n),met(9,n)) 

      ! Potential latent energy (kg.m-2.day-1)
      if (deltaWP < dble_zero) then
          call acm_et(n,deltat(n),gpppars(6),rainfall_input,meant,met(4,n),met(10,n),met(11,n),met(9,n) &
                     ,evappars,gpppars(9),gpppars(10),canopy_storage &
                     ,intercepted_rainfall,FLUXES(n,3),FLUXES(n,2),FLUXES(n,4))
      else
          FLUXES(n,2) = dble_zero 
          FLUXES(n,3) = dble_zero
          FLUXES(n,4) = dble_zero
      endif
      ! determine potential evaporation (excluding wet canopy surface)
      ET_pot = FLUXES(n,2) + FLUXES(n,4)
      ! do mass balance (i.e. is there enough water to support ET)
      ET_net = calculate_update_soil_water(ET_pot*deltat(n), & 
                                          ((rainfall_input-intercepted_rainfall)*seconds_per_step),meant)
      ! now reverse the time correction (step -> day)
      ET_net = ET_net * deltat_1(n)
      ! update the specific soil and transpiration fluxes
      if (ET_net > dble_zero) then
          FLUXES(n,2) = FLUXES(n,2) * (ET_net/ET_pot)
          FLUXES(n,4) = FLUXES(n,4) * (ET_net/ET_pot)
      endif
      ! store soil water content of the rooting zone (mm)
      POOLS(n,1) = 1e3*sum(soil_waterfrac(1:nos_root_layers)*layer_thickness(1:nos_root_layers))

      !!!!!!!!!!
      ! GPP (gC.m-2.day-1)
      !!!!!!!!!!

      if (FLUXES(n,2) > vsmall .and. met(11,n) > 1e-10 .and. &
          deltaWP < dble_zero .and. abs(gpppars(10)) /= abs(log(infi))) then
         FLUXES(n,1) = acm_gpp(gpppars,constants)
      else
         FLUXES(n,1) = dble_zero
      endif

      !!!!!!!!!!
      ! Bug checking
      !!!!!!!!!!

      do nxp = 1, nopools
         if (POOLS(n+1,nxp) /= POOLS(n+1,nxp) .or. POOLS(n+1,nxp) < 0) then
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
         if (FLUXES(n,nxp) /= FLUXES(n,nxp)) then
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

    end do ! nodays loop

  end subroutine CARBON_MODEL
  !
  !------------------------------------------------------------------
  !
  double precision function acm_gpp(drivers,constants)

    ! the Aggregated Canopy Model, is a Gross Primary Productivity (i.e.
    ! Photosyntheis) emulator which operates at a daily time step. ACM can be
    ! paramaterised to provide reasonable results for most ecosystems.

    implicit none

    ! declare input variables
    double precision, intent(in) :: drivers(12) & ! acm input requirements
                         ,constants(10) ! ACM parameters

    ! declare local variables
    double precision :: gc, pn, pd, pp, qq, ci, e0, dayl, cps, nit &
             ,mult, light_gpp, mint, maxt, radiation, co2, lai    &
             ,deltaWP,Rtot,NUE,temp_exponent,dayl_coef &
             ,dayl_const,hydraulic_exponent,hydraulic_temp_coef &
             ,co2_comp_point,co2_half_sat,lai_coef,lai_const

    ! load driver values to correct local vars
    lai = drivers(1)
    maxt = drivers(2)
    mint = drivers(3)
    nit = drivers(4)   
    co2 = drivers(5)
    dayl = drivers(6)
    radiation = drivers(8)

    ! load parameters into correct local vars
    deltaWP = drivers(9)
    Rtot = drivers(10)
    NUE = constants(1)
    dayl_coef = constants(2)
    co2_comp_point = constants(3) 
    co2_half_sat = constants(4)
    dayl_const = constants(5)
    hydraulic_temp_coef = constants(6)
    lai_coef = constants(7)
    temp_exponent = constants(8)
    lai_const = constants(9)
    hydraulic_exponent = constants(10)

!    ! Temperature adjustments for Michaelis-Menten coefficients 
!    ! for CO2 (kc) and O2 (ko) and CO2 compensation point
!    co2_half_sat   = arrhenious(kc_saturation,kc_half_sat_conc,maxt)
!    co2_comp_point = arrhenious(co2comp_saturation,co2comp_half_sat_conc,maxt)

    ! daily canopy conductance (m.s-1) ; NOTE that ratio of H20:CO2 diffusion is
    ! 1.646259 (Jones appendix 2). i.e. gcCO2/1.646259 = gcH2O
    ! also note that deltaWP has been set outside of this function as absolute
    ! value even through the actual variable would be negative
    gc = deltaWP**(hydraulic_exponent)/((hydraulic_temp_coef*Rtot))
    ! maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
    pn = lai*nit*NUE*exp(temp_exponent*maxt)
    ! pp and qq represent limitation by diffusion and metabolites respecitively
    pp = pn/gc ; qq = co2_comp_point-co2_half_sat
    ! calculate internal CO2 concentration (ppm)
    mult = co2+qq-pp
    ci = 0.5*(mult+sqrt((mult*mult)-4.0*(co2*qq-pp*co2_comp_point)))
    ! Michaelis–Menten limitation of maximum quantium efficiency by leaf area
    ! (i.e. light interception)
    mult = lai*lai
    e0 = lai_coef*mult/(mult+lai_const)
    ! calculate CO2 limited rate of photosynthesis
    pd = gc*(co2-ci)
    ! calculate light limted rate of photosynthesis
    light_gpp = e0 * radiation
    ! calculate combined light and CO2 limited photosynthesis
    cps = light_gpp*pd/(light_gpp+pd)
    ! correct for day length variation
    acm_gpp = cps*(dayl_coef*dayl+dayl_const)
    ! don't forget to return
    return

  end function acm_gpp
  !
  !------------------------------------------------------------------
  !
  subroutine acm_et(n,deltat,dayl,rainfall,meant,swrad_MJ,vpd_pa,lai,wind_spd,pars,deltaWP,Rtot &
                   ,canopy_storage,intercepted_rainfall,wetcanopy_evap,transpiration,soilevap)

    ! Three response function(s) based on the Penman-Monteith model of
    ! evapotranspiration used to estimate SPA's daily evapotranspiration flux
    ! (kg.m-2.day-1). 
    ! Function 1 calculates the potential canopy level evaporation.
    ! Function 2 calculates the transpiration linking to hydraulic limitations. 
    ! F1-F2 = potential canopy wet surface evaporation. 
    ! Function 3 quantifies the potential soil surface evaporation flux

    implicit none

    ! arguments
    integer, intent(in) :: n
    double precision, intent(in) :: deltat   & ! days in time step
                                   ,dayl     & ! day length in hours
                                   ,rainfall & ! rainfall (kg.m-2.s-1)
                                   ,meant    & ! daily mean temperature (oC)
                                   ,swrad_MJ & ! daily short wave rad (MJ.m-2.day-1)
                                   ,vpd_pa   & ! daily mean vapoure pressure deficit (Pa)
                                   ,wind_spd & ! wind speed (m.s-1)
                                   ,lai      & ! daily leaf area index (m2/m2)
                                   ,deltaWP  & ! soil-minleafwp (MPa)
                                   ,Rtot       ! total soil root hydraulic resistance (MPa.s-1.m-2.mmol-1)

    double precision, intent(in), dimension(8) :: pars
    double precision, intent(inout) :: canopy_storage, intercepted_rainfall
    double precision, intent(out) :: wetcanopy_evap, transpiration, soilevap ! kg.m-2.day-1

    ! local variables
    double precision :: mult &
                    ,meant_K &
               ,dayl_seconds &
                      ,swrad & ! daily mean radiation
                      ,lwrad &
           ,longwave_release &
                      ,decay &
     ,total_soil_conductance &
                        ,gws & !
            ,water_diffusion & ! 
                    ,fun_rad &
                   ,fun_Rtot &
                    ,s,slope &
                      ,psych &
                 ,lai_albedo & 
         ,canopy_temperature &
           ,canopy_radiation &
             ,soil_radiation &
           ,soil_conductance & 
                     ,lambda &
                        ,rho &
                  ,Jm3kPaK_1 &
                         ,gc
  
    ! convert meant air temperature to Kelvin
    meant_K = meant + freeze
    dayl_seconds = dayl * 3600.0
    soil_radiation = dble_zero ; canopy_radiation = dble_zero

    !!!!!!!!!!
    ! Estimate isotermal longwave radiation balance (W.m-2)
    !!!!!!!!!!

    ! estimate long wave radiation from atmosphere (W.m-2)
    lwrad = emiss_boltz * (meant_K-20.0) ** 4
    ! estimate isothermal long wave emission per unit area
    longwave_release = emiss_boltz * meant_K ** 4

    ! calculate integral of (Beer's law) exponential decay 
    ! (0.5 = decay coefficient assuming mean angle of interception 30 degrees,
    ! see SPA light.f90)
    decay = min(dble_one,(dble_one-exp(-0.5 * lai)) / 0.5)

    ! Estimate long wave incident on canopy from sky and soil surface.
    ! Canopy absorption of longwave is assumed to be equal to emissivity.
    ! Likewise, include assumption of long wave from soil surface in isothermal
    ! conditions.
    canopy_radiation = emissivity * lwrad * decay + longwave_release * decay * emissivity

    ! assume that longwave emitted from both sides of the leaves are largely
    ! re-absorbed by other parts of the canopy. i.e. losses area limited to the
    ! exposed canopy at top and bottom.
    canopy_radiation = canopy_radiation - &
                     ( 2.0 * longwave_release * min(dble_one,lai) )

    ! assign longwave which makes it through the canopy to soil surface
    ! plus longwave from the canopy itself which we assume to pass down the
    ! canopy
    soil_radiation = emissivity * lwrad * ( dble_one - decay ) + &
                    (emissivity * longwave_release * min(dble_one,lai))
    ! less long wave loss from the surface assuming isothermal conditions
    soil_radiation = soil_radiation - longwave_release

    !!!!!!!!!!
    ! Estimate shortwave radiation balance (W.m-2)
    !!!!!!!!!

    ! convert incoming shortwave radiation from MJ.m-2.day-1 -> J.m-2.day-1
    ! as shortwave radiation will infact be more concentrated (i.e. comes down
    ! in day time only) we average radiation over this period only
    ! calculate day length (seconds)
    swrad = (swrad_MJ*1e6) / dayl_seconds
!    swrad = swrad_MJ*1e6*seconds_per_day_1

    ! Calculate leaf area interaction on radiation balance
    mult = lai*lai
    lai_albedo = pars(1)*mult/(mult+pars(2))
    ! Estimate net radiation balance for canopy 
    canopy_radiation = canopy_radiation + (swrad * lai_albedo)

    ! Assuming the same leaf area sensitivity, but different soil albedo maximum
    ! quantify the effective albedo of soil surface
    soil_radiation = soil_radiation + (swrad * (dble_one-lai_albedo) * pars(7))

    !!!!!!!!!!
    ! Calculate canopy conductance (to water vapour)
    !!!!!!!!!!

    ! calculate canopy conductance of evaporation. Assumes logitistic functions
    ! to radiation and hydraulic resistance
    fun_rad  = opt_max_scaling(pars(3),pars(4),pars(8),swrad)
    fun_Rtot = (dble_one + exp(pars(5)*(Rtot-pars(6))))**(-dble_one)
          gc = max(0.00005,deltaWP * fun_rad * fun_Rtot)

    !!!!!!!!!!
    ! Update canopy energy balance
    !!!!!!!!!!
    if (lai > dble_zero) then
        ! Quantify leaf / canopy temperature difference
        canopy_temperature = calculate_canopy_temperature(meant,meant_K,lai,canopy_radiation,vpd_pa,gc,decay)
        ! update canopy radiation balance
        canopy_radiation = canopy_radiation - (lai*4.0*emiss_boltz*meant_K**3*(canopy_temperature-meant))
    else
        canopy_temperature = meant
    endif

    ! calculate coefficient for Penman Montieth
    ! density of air (kg.m-3)
    rho = 353.0/(canopy_temperature+freeze)
    if (meant < dble_one) then
        lambda = 2.835e6
    else
        ! latent heat of vapourisation (J.kg-1)
        lambda = 2501000.0-2364.0*canopy_temperature
    endif
    ! psychrometric constant (kPa K-1)
    psych = (0.0646*exp(0.00097*canopy_temperature))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult = canopy_temperature+237.3
    ! 2502.935945 = 0.61078*17.269*237.3
    s = 2502.935945*exp(17.269*canopy_temperature/mult)
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = s/(mult*mult)
    ! Calculate energy change due to VPD per kelvin (multiple use value)
    Jm3kPaK_1 = rho*cpair*vpd_pa*1e-3

    ! Calculate numerator of Penman Montheith (kg.m-2.day-1); note vpd Pa->kPa
    wetcanopy_evap = (slope*canopy_radiation*dayl_seconds) + (Jm3kPaK_1*aerodynamic_conductance)
    ! Calculate the transpiration flux 
    transpiration = max(dble_zero,wetcanopy_evap / (lambda*(slope+(psych*(dble_one+aerodynamic_conductance/gc)))))
    ! Calculate the potential wet canopy evaporation
    wetcanopy_evap = max(dble_zero,(wetcanopy_evap / (lambda*(slope+psych))) - transpiration )
    ! Update based on canopy water storage
    call canopy_interception_and_storage(n,deltat,rainfall,lai, &
                                         wetcanopy_evap,canopy_storage,intercepted_rainfall)

    ! Estimate water diffusion rate (m.s-1)
    water_diffusion = 24.2e-6 * ( meant_K / 293.2 )**1.75
    ! Soil conductance to water vapour diffusion (m s-1)...
    gws = porosity(1) * water_diffusion / ( tortuosity * drythick )
    ! Calculate the soil conductance (m.s-1)
    call soil_surface_conductance(lai,wind_spd,soil_conductance)
    ! Combine to total soil conductance 
    total_soil_conductance = ( gws**(-dble_one) + soil_conductance**(-dble_one) )**(-dble_one)
    ! Estimate potential soil evaporation flux (kg.m-2.day-1)
    soilevap = (slope*soil_radiation*dayl_seconds) + (Jm3kPaK_1*soil_conductance)
    soilevap = max(dble_zero,soilevap / (lambda*(slope+(psych*(dble_one+total_soil_conductance)))))

  end subroutine acm_et
  !
  !------------------------------------------------------------------
  !
   subroutine calculate_aerodynamic_conductance(lai,wind_spd)

    ! 
    ! Calculates the aerodynamic or bulk canopy conductance (m.s-1). Here we
    ! assume neutral conditions due to the lack of an energy balance calculation
    ! in either ACM or DALEC. The equations used here are with SPA at the time
    ! of the calibration
    ! 

    implicit none

    ! arguments
    double precision, intent(in) :: lai, wind_spd

    ! calculate the zero plane displacement and roughness length
    call z0_displacement(lai)

    ! calculate bulk conductance (Jones p68)
    aerodynamic_conductance = (wind_spd * vonkarman_2) / (log((canopy_height-displacement)/roughl))**2

  end subroutine calculate_aerodynamic_conductance
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_field_capacity

    ! field capacity calculations for saxton eqns !

    implicit none

    ! local variables..
    integer        :: i

    x1 = 0.1 ; x2 = 0.7 ! low/high guess
    do i = 1 , nos_soil_layers+1
       water_retention_pass = i
       ! field capacity is water content at which SWP = -10 kPa
       field_capacity(i) = zbrent('water_retention:water_retention_saxton_eqns', &
                                   water_retention_saxton_eqns , x1 , x2 , xacc )
    enddo

  end subroutine calculate_field_capacity
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_Rtot(deltaWP,Rtot,lai,deltat,meant)
  
    ! Purpose of this subroutine is to calculate the minimum soil-root hydraulic
    ! resistance input into ACM. The approach used here is identical to that
    ! found in SPA.

    ! declare inputs
    double precision,intent(in) :: deltat,deltaWP,lai,meant
    double precision,intent(inout) :: Rtot

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
    transpiration_resistance = (gplant * lai)**(-dble_one)
!    transpiration_resistance = canopy_height / (gplant * lai)

    ! The original SPA src generates an exponential distribution which aims
    ! to maintain 50 % of root biomass in the top 25 % of the rooting depth.
    ! In a simple 2 root layer system this can be estimates more simply

    ! top 25 % of root profile
    slpa = (root_reach * 0.25) - layer_thickness(1) 
    if (slpa <= dble_zero) then
        ! > 50 % of root is in top layer
        root_mass(1) = root_biomass * 0.5
        root_mass(1) = root_mass(1) + ((root_biomass-root_mass(1)) * (abs(slpa)/root_reach))
    else
        ! < 50 % of root is in bottom layer
        root_mass(1) = root_biomass * 0.5 * (layer_thickness(1)/(abs(slpa)+layer_thickness(1)))
    endif
    root_mass(2) = max(dble_zero,root_biomass - root_mass(1))
    root_length = root_mass * root_mass_length_coef_1
!    root_length = root_mass / (root_density * root_cross_sec_area)
    !! Top root layer.
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
    demand = deltaWP+head*canopy_height
    water_flux(1) = demand/(transpiration_resistance + soilR1 + soilR2)
    ! Bottom root layer
    if (root_mass(2) > dble_zero ) then
       ! soil conductivity converted from m.s-1 -> m2.s-1.MPa-1 by head
       soilR1=soil_resistance(root_length(2),layer_thickness(2),soil_conductivity(2)*head_1)
       soilR2=root_resistance(root_mass(2),layer_thickness(2))
       soilRT_local(2)=soilR1 + soilR2 + transpiration_resistance
       ! calculate and accumulate steady state water flux in mmol.m-2.s-1
       water_flux(2) = demand/(transpiration_resistance + soilR1 + soilR2)
       ratio = layer_thickness(1:nos_root_layers)/sum(layer_thickness(1:nos_root_layers))
    endif ! roots present in second layer?

    ! if freezing then assume soil surface is frozen
    if (meant < dble_one) then 
        water_flux(1) = dble_zero ; ratio(1) = dble_zero
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
        wSWP = -20.0 ; soilRT = sum(soilRT_local)*0.5
        uptake_fraction = dble_zero ; uptake_fraction(1) = dble_one
    endif

    ! determine effective resistance
    Rtot = demand / sum(water_flux*ratio)
    ! finally convert transpiration flux (mmol.m-2.s-1) 
    ! into kg.m-2.step-1 for consistency with ET in "calculate_update_soil_water"
    water_flux = water_flux * mmol_to_kg_water * seconds_per_step

    ! and return
    return

  end subroutine calculate_Rtot
  !
  !-----------------------------------------------------------------
  !
  subroutine canopy_interception_and_storage(n,deltat,rainfall,lai,potential_evaporation &
                                            ,canopy_storage,intercepted_rainfall)

    ! Description

    implicit none

    ! arguments
    integer, intent (in) :: n
    double precision, intent(in) :: deltat, &
                                  rainfall, & ! incoming rainfall (kg.m-2.s-1)
                                       lai    ! current leaf area index (m2/m2) 
    double precision, intent(inout) :: potential_evaporation & ! wet canopy evaporation (kg.m-2.day-1), 
                                                               ! enters as potential but leaves as water balance adjusted 
                                             ,canopy_storage   ! water currently stored on the canopy (kg.m-2)
    double precision, intent(out) :: intercepted_rainfall ! total rainfall intercepted (kg.m-2.s-1)
                                                          ! NOTE: it is possible for this to be
                                                          ! negative if stored water running off
                                                          ! into the soil is greater than
                                                          ! rainfall (i.e. when leaves hahve died between steps) 
    ! local variables
    integer :: i
    double precision :: tmp, through_fall, max_storage, max_storage_1 &
                       ,daily_addition, wetcanopy_evaporation

    ! determine maximum canopy storage & through fall fraction (mm)
    through_fall=max(min_throughfall,exp(-0.5*lai))
    ! maximum canopy storage (mm); minimum is applied to prevent errors in
    ! drainage calculation. Assume minimum capacity due to wood
    max_storage=max(min_storage,0.1*lai) ; max_storage_1 = max_storage**(-dble_one)
    ! potential intercepted rainfall (kg.m-2.s-1)
    intercepted_rainfall = rainfall * (dble_one - through_fall)
    ! average rainfall intercepted by canopy (kg.m-2.day-1)
    daily_addition = intercepted_rainfall * seconds_per_day

    tmp = dble_zero ; through_fall = dble_zero ; wetcanopy_evaporation = dble_zero 
    if (daily_addition > potential_evaporation .or. &
        daily_addition > max_storage .or.           &
        canopy_storage > vsmall) then

        ! intergrate over canopy for each day
        do i = 1, int(deltat)

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

        end do ! day looping

        ! sanity checks 
        ! NOTE: addition of 1e-10 is to prevent precision error causing stop
!        if (canopy_storage > (max_storage*(dble_one+1e-10)) .or. canopy_storage < dble_zero) then
!           print*,"Canopy water storage mass balance error!!"
!           print*,"through_fall_total",through_fall
!           print*,"canopy_storage",canopy_storage,"max_storage",max_storage
!           print*,"potential_evaporation",potential_evaporation,"actual",wetcanopy_evaporation / deltat
!           stop
!        endif

        ! average fluxes to daily rates
        potential_evaporation = wetcanopy_evaporation * deltat_1(n)
        ! correct intercepted rainfall rate to kg.m-2.s-1    
        intercepted_rainfall = intercepted_rainfall - ((through_fall * deltat_1(n)) * seconds_per_day_1)

!        ! sanity check
!        if (intercepted_rainfall > rainfall) then
!            print*,"Canopy intercepted rainfall cannot be greater than rainfall!!"
!            print*,"rainfall", rainfall, "through_fall", ((through_fall / deltat) * seconds_per_day_1)
!        endif

    else

       ! we can assume that all water will be evaporated and daily averages can
       ! be assumed for the whole time step
       
       ! assuming all intercepted rainfall was evaporated
       potential_evaporation = daily_addition
       ! intercepted rainfall therefore remains unchanged...
       
       ! canopy water storage is also assumed to remain empty as all intercepted
       ! water was evaporated away

    end if

  end subroutine canopy_interception_and_storage
  !
  !-----------------------------------------------------------------
  !
  subroutine infiltrate(rainfall)

    ! Takes surface_watermm and distrubutes it among top !
    ! layers. Assumes total infilatration in timestep.   !

    implicit none

    ! arguments 
    double precision, intent(in) :: rainfall ! rainfall (kg.m-2.step-1)

    ! local argumemts
    integer :: i
    double precision    :: add   & ! surface water available for infiltration (m)
                          ,wdiff   ! available space in a given soil layer for water to fill (m)

    ! convert rainfall water from mm -> m (or kg.m-2.step-1 -> Mg.m-2.step-1)
    add = rainfall * 1e-3 
    do i = 1 , nos_soil_layers
       ! determine the available pore space in current soil layer
       wdiff = max(dble_zero,(porosity(i)-soil_waterfrac(i))*layer_thickness(i)-watergain(i)+waterloss(i))
       ! is the input of water greater than available space
       ! if so fill and subtract from input and move on to the next
       ! layer
       if (add > wdiff) then
          ! if so fill and subtract from input and move on to the next layer
          watergain(i) = watergain(i)+wdiff
          add = add-wdiff
       else
          ! otherwise infiltate all in the current layer
          watergain(i) = watergain(i)+add
          add = dble_zero
       end if
       ! if we have added all available water we are done
       if (add <= dble_zero) then
           add = dble_zero
           exit
       end if
    end do

    ! if after all of this we have some water left assume it is runoff
    runoff = add * 1e3

  end subroutine infiltrate
  !
  !-----------------------------------------------------------------
  !
  subroutine gravitational_drainage(meant)

    ! integrator for soil gravitational drainage !

    implicit none

    ! arguments
    double precision, intent(in) :: meant ! daily mean temperature (oC)

    ! local variables..
    double precision  :: change, drainage, iceprop(nos_soil_layers)


    ! calculate soil ice proportion; at the moment
    ! assume everything liquid
    iceprop = dble_zero
    ! except the surface layer in the mean daily temperature is < 1oC
    if (meant < dble_one) iceprop(1) = dble_one

    do soil_layer = 1, nos_soil_layers

       ! liquid content of the soil layer, i.e. fraction avaiable for drainage
       liquid     = soil_waterfrac( soil_layer ) * ( dble_one - iceprop( soil_layer ) )     ! liquid fraction
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
    double precision, parameter :: H = 0.332, &
                                   J = -7.251e-4, &
                                   K = 0.1276

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

  end subroutine initialise_soils
  !
  !-----------------------------------------------------------------
  !
  subroutine calculate_soil_conductivity

    ! Calculate the soil conductivity (m s-1) of water based on soil
    ! characteristics and
    ! current water content

    implicit none

    ! soil conductivity for the dynamic soil layers (i.e. not including core)
    soil_conductivity(1:nos_soil_layers) = cond1(1:nos_soil_layers) &
                                        * exp(cond2(1:nos_soil_layers)+cond3(1:nos_soil_layers)/soil_waterfrac(1:nos_soil_layers))

    ! protection against floating point error
    where (soil_waterfrac < 0.05)
          soil_conductivity = 1e-30
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
    double precision, parameter :: A = -4.396, B = -0.0715, CC = -4.880e-4, D = -4.285e-5, &
                                   E = -3.140, F = -2.22e-3, G = -3.484e-5, H = 0.332,     &
                                   J = -7.251e-4, K = 0.1276, P = 12.012, Q = -7.551e-2,   &
                                   R = -3.895, T = 3.671e-2, U = -0.1103, V = 8.7546e-4,   &
                                   mult1 = 100.0, mult2 = 2.778e-6, mult3 = 1000.0

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
  subroutine soil_surface_conductance(lai,wind_spd,soil_conductance)

    ! proceedsure to solve for soil surface resistance based on Monin-Obukov
    ! similarity theory stability correction momentum & heat are integrated 
    ! through the under canopy space and canopy air space to the surface layer
    ! references are Nui & Yang 2004; Qin et al 2002
    ! NOTE: conversion to conductance at end

    implicit none

    ! declare arguments
    double precision, intent(in) :: lai, wind_spd
    double precision, intent(out) :: soil_conductance ! soil conductance to atmopshere for exchange (m.s-1)

    ! local variables
    double precision :: canopy_decay & ! canopy decay coefficient for soil exchange
                       ,ustar        & ! friction velocity (m.s-1)
                       ,beta         & ! ratio of u*/Uh
                       ,lc           & ! length scale for vertical momentum absorption within the canopy
                       ,lm           & ! mixing length (m)
                       ,Uh           & ! wind speed at canopy top (m.s-1)
                       ,Kh_canht       ! eddy diffusivity at canopy height (m2.s-1)

    ! parameters
    double precision, parameter :: foliage_drag = 0.2, & ! foliage drag coefficient
                                   beta_max = 1.0, beta_min = 0.2, min_lai = 1.0, &
                                   most_soil = 1.0 ! Monin-Obukov similarity theory stability correction.
                                                   ! As no sensible heat flux calculated, 
                                                   ! assume neutral conditions only

    ! calculate friction velocity at tower height (reference height) (m.s-1) 
    ! WARNING neutral conditions only; see WRF module_sf_sfclay.F for 'with
    ! stability versions'
    ustar = (wind_spd / log((tower_height-displacement)/roughl)) * vonkarman

    ! calculate this steps beta coefficient; minimum value see Finnigan & Harman (2007)
    Uh   = wind_spd ! should really be log-law decayed wind speed at top of canopy
    beta = min(beta_max,max(ustar/Uh,beta_min))

    ! both length scale and mixing length are considered to be constant within
    ! the canopy (under dense canopy conditions)
    ! calculate length scale (lc) for momentum absorption within the canopy; Harman & Finnigan (2007)
    ! and mixing length (lm) for vertical momentum within the canopy; Harman & Finnigan (2008)
    if (lai > min_lai) then
        lc = (4.0*canopy_height) / lai
        lm = max(canopy_height*0.02, 2.0*(beta**3)*lc)
    else
!        lc = vonkarman * tower_height
        lm = canopy_height * vonkarman
    endif

    ! calculate eddy diffusivity at the top of the canopy (m2.s-1) 
    ! Kaimal & Finnigan 1994; for near canopy approximation
    Kh_canht=vonkarman*ustar*(canopy_height-displacement)

    ! calculate canopy decay coefficient with stability correction
    ! NOTE this is not consistent with canopy momentum decay done by Harman &
    ! Finnigan (2008)
    canopy_decay = (((foliage_drag*canopy_height*max(min_lai,lai))/lm)**0.5)*(most_soil**0.5)
    ! approximation of integral for soil resistance; maximum soil surface
    ! resistance applied at 200 s.m-1, conductance = 5e-3 m.s-1
    soil_conductance=canopy_height/(canopy_decay*Kh_canht) & 
                   *(exp(canopy_decay*(dble_one-(soil_roughl/canopy_height)))- &
                     exp(canopy_decay*(dble_one-((roughl+displacement)/canopy_height))))

    ! convert resistance (s.m-1) to conductance (m.s-1)
    soil_conductance = soil_conductance ** (-dble_one)

  end subroutine soil_surface_conductance
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
    SWP(1:nos_soil_layers) = -0.001 * potA(1:nos_soil_layers) &
                           * soil_waterfrac(1:nos_soil_layers)**potB(1:nos_soil_layers)
    where (SWP(1:nos_soil_layers) < -20) SWP(1:nos_soil_layers) = -20
!    where (soil_waterfrac(1:nos_soil_layers) < 0.005)
!        SWP(1:nos_soil_layers) = -9999.0
!    end where

  end subroutine soil_water_potential
  ! 
  !------------------------------------------------------------------
  !
  subroutine z0_displacement(lai)

    ! dynamic calculation of roughness length and zero place displacement (m)
    ! based on canopy height and lai. Raupach (1994)

    implicit none

    ! arguments
    double precision, intent(in) :: lai
    ! local variables
    double precision  sqrt_cd1_lai &
                     ,local_lai &
                     ,ustar_Uh  & ! ratio of friction velocity over wind speed at canopy top
                     ,phi_h       ! roughness sublayer influence function
    double precision, parameter :: cd1 = 7.5,   & ! Canopy drag parameter; fitted to data
                                    Cs = 0.003, & ! Substrate drag coefficient
                                    Cr = 0.3,   & ! Roughness element drag coefficient
                          ustar_Uh_max = 0.3,   & ! Maximum observed ratio of (friction velocity / canopy top wind speed) (m.s-1)
                               min_lai = 1.0,   & ! Minimum LAI parameter as height does not vary with growth
                                    Cw = 2.0      ! Characterises roughness sublayer depth (m)

    ! assign new value to min_lai to avoid max min calls
    local_lai = max(min_lai,lai)
    sqrt_cd1_lai = sqrt(cd1 * local_lai)

    ! calculate displacement (m); assume minimum lai 1.0 or 1.5 as height is not
    ! varied
    displacement=(dble_one-((dble_one-exp(-sqrt_cd1_lai))/sqrt_cd1_lai))*canopy_height

    ! calculate estimate of ratio of friction velocity / canopy wind speed; with
    ! max value set at
    ustar_Uh=min(sqrt(Cs+Cr*local_lai*0.5),ustar_Uh_max)
    ! calculate roughness sublayer influence function; 
    ! this describes the departure of the velocity profile from just above the
    ! roughness from the intertial sublayer log law
    phi_h = 0.19314718056
!    phi_h = log(Cw)-dble_one+Cw**(-dble_one) ! DO NOT FORGET TO UPDATE IF Cw CHANGES

    ! finally calculate roughness length, dependant on displacement, friction
    ! velocity and lai.
    roughl=((dble_one-displacement/canopy_height)*exp(-vonkarman*ustar_Uh-phi_h))*canopy_height

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

    numerator   = t - 25.0
    denominator = t + freeze
    answer      = a * exp( b * dble_one * numerator / denominator )
    arrhenious  = answer

  end function arrhenious
  !
  !------------------------------------------------------------------ 
  !
  double precision function calculate_canopy_temperature(meant,meant_K,lai,canopy_radiation,vpd_pa,gc,decay)

    ! function estimates the steady state solution to canopy temperature by
    ! balancing approximate turbulent fluxes

    implicit none

    ! arguments
    double precision, intent(in) :: meant, & ! air temperature (oC)
                                  meant_K, & ! air temperatyre (K)
                                      lai, &
                         canopy_radiation, & ! estimate of net absorbed radiation (SW+LW) J.m-2.day-1
                                   vpd_pa, & ! vapour pressure deficit (Pa)
                                    decay, &
                                       gc    ! canopy conductance (m.s-1)

    ! local variables
    double precision :: rho, lambda, psych, mult, s, slope, &
                        thermal_gains, thermal_losses, sm_1_kPaK_1,    &
                        canopy_evaporative_resistance,                 &
                        canopy_thermal_resistance,radiative_conductance

    ! calculate coefficients for Penman Montieth
    ! density of air (kg.m-3)
    rho = 353.0/meant_K
    if (meant < dble_one) then
        lambda = 2.835e6
    else
        ! latent heat of vapourisation (J.kg-1)
        lambda = 2501000.0-2364.0*meant
    endif
    ! psychrometric constant (kPa K-1)
    psych = (0.0646*exp(0.00097*meant))
    ! Straight line approximation of the true slope; used in determining
    ! relationship slope
    mult = meant+237.3
    ! 2502.935945 = 0.61078*17.269*237.3
    s = 2502.935945*exp(17.269*meant/mult)
    ! Rate of change of saturation vapour pressure with temperature (kPa.K-1)
    slope = s/(mult*mult)

    !!!!!!!!!!
    ! Calculate total thermal resistance in form of radiation (gr) and sensible
    ! heat. Units (m.s-1)
    !!!!!!!!!!
    
    ! calculate conductance rate for radiative heat exchange (m.s-1)
    radiative_conductance = (4.0 * emiss_boltz * meant_K ** 3) / (rho * cpair)
    ! scale for leaf area
    radiative_conductance =  min(dble_one,lai) * radiative_conductance !* (dble_one - decay)
    ! combine with conductance of sensible heat (aerodynamic_conductance) in
    ! parallel as these processes are indeed occuring in parallel. 
    ! NOTE: that aerodynamic_conductance already implicitly contains LAI scaling
    ! from roughness length and displacement height components
    canopy_thermal_resistance = radiative_conductance + aerodynamic_conductance! * 2

    ! NOTE: conversion to resistance (s.m-1)
    canopy_thermal_resistance = (canopy_thermal_resistance) ** (-dble_one)

    ! MAY NEED TO CARRY FORWARD LONG WAVE ASSUMPTION THAT HALF LOST WILL BE
    ! FOUND BY OTHER PARTS OF THE CANOPY...

    !!!!!!!!!
    ! Calculate total evaporation resistance (s.m-1)
    !!!!!!!!!

    ! Canopy conductance (lai scaled anologue of stomatal conductance) is
    ! combined in series with areodynamic conductance. 
    ! NOTE: conversion of conductances to resistance
    if (gc > dble_zero) then
        canopy_evaporative_resistance = gc ** (-dble_one) + aerodynamic_conductance ** (-dble_one)
    else
        canopy_evaporative_resistance = aerodynamic_conductance ** (-dble_one)
    endif

    !!!!!!!!!
    ! Determine energy gains by canopy and then losses
    !!!!!!!!!

    ! calculate common denominator (units: sm-1 kPaK-1)
    sm_1_kPaK_1 = psych * canopy_evaporative_resistance + slope * canopy_thermal_resistance
    ! calculate thermal gains to the system (K)
    thermal_gains = (dble_one + canopy_thermal_resistance * canopy_evaporative_resistance * psych * (canopy_radiation/lai)) 
    thermal_gains = thermal_gains / (rho * cpair * sm_1_kPaK_1)
    ! calculate thermal losses to the system (K)
    thermal_losses = (canopy_thermal_resistance * vpd_pa * 1e-3) / sm_1_kPaK_1
    ! calculate net difference in canopy temperature and thus new canopy
    ! temperature (K)
    calculate_canopy_temperature = meant + (thermal_gains - thermal_losses)

    return

  end function calculate_canopy_temperature
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
    dec = - asin( sin_dayl_deg_to_rad * cos( two_pi * ( doy + 10.0 ) / 365.0 ) )
    mult = lat * deg_to_rad
    sinld = sin( mult ) * sin( dec )
    cosld = cos( mult ) * cos( dec )
    aob = max(-dble_one,min(dble_one,sinld / cosld))
    
    ! define output
    daylength_hours = 12.0 * ( dble_one + 2.0 * asin( aob ) * pi_1 )

    ! return to user
    return

  end function daylength_hours
  !
  !------------------------------------------------------------------
  !
  double precision function calculate_update_soil_water(ET,rainfall,meant)

   !
   ! 1) Limits ET by available water in the soil
   ! 2) Updates soil water balance based on ET, drainage and rainfall
   !

   implicit none

   ! arguments
   double precision, intent(in) :: ET & ! evapotranspiration estimate (kg.m-2.step-1)
                            ,rainfall & ! rainfall (kg.m-2.step-1)
                               ,meant   ! daily mean temperature (oC)

   ! local variables
   double precision ::  depth_change
   double precision, dimension(nos_root_layers) :: avail_flux, evaporation_losses

   ! seperately calculate the soil conductivity as this applies to each layer
   call calculate_soil_conductivity

   ! for simplicity assume that all evaporation occurs in same distribution as
   ! transpiration. As the surface layer will also have the bulk of the roots
   ! too and steady state flux from surface is linked to water availability too.
   evaporation_losses = ET * uptake_fraction

   ! limit water losses from each layer to that available
   avail_flux = soil_waterfrac(1:nos_root_layers) & 
              * layer_thickness(1:nos_root_layers) * 1e3
   do soil_layer = 1, nos_root_layers
      if (evaporation_losses(soil_layer) > avail_flux(soil_layer)) then 
         ! just to give a buffer against numerical precision error make
         ! extraction slightly less than available
         evaporation_losses(soil_layer) = avail_flux(soil_layer) * 0.99
      endif
   end do

!   where (evaporation_losses < dble_zero) evaporation_losses = dble_zero
   ! this will update the ET estimate outside of the function
   ! unit / time correction also occurs outside of this function
   calculate_update_soil_water = sum(evaporation_losses)

   ! pass information to waterloss variable and zero watergain 
   ! convert kg.m-2 (or mm) -> Mg.m-2 (or m)
   waterloss = dble_zero ; watergain = dble_zero
   waterloss(1:nos_root_layers) = evaporation_losses(1:nos_root_layers)*1e-3 
   ! Update soil provides based on evaporative losses to avoid double movement
   ! of water from drainage
   ! Convert fraction into depth specific values (rather than m3/m3) then update soil water fractions
   soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness) &
                                        + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                     / layer_thickness(1:nos_soil_layers)
   ! now re-zero water losses
   waterloss = dble_zero ; watergain = dble_zero
   ! determine drainage flux between surface -> sub surface and sub surface
   call gravitational_drainage(meant)
   ! update soil water profile again
   soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness) &
                                        + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                     / layer_thickness(1:nos_soil_layers)
   ! now re-zero water losses
   waterloss = dble_zero ; watergain = dble_zero
   ! determine infiltration from rainfall,
   ! if rainfall is probably liquid / soil surface is probably not frozen
   if (meant >= dble_one .and. rainfall > dble_zero) then
       call infiltrate(rainfall)
       ! update soil profiles. Convert fraction into depth specific values (rather than m3/m3) then update fluxes
       soil_waterfrac(1:nos_soil_layers) = ((soil_waterfrac(1:nos_soil_layers)*layer_thickness) &
                                            + watergain(1:nos_soil_layers) - waterloss(1:nos_soil_layers)) &
                                         / layer_thickness(1:nos_soil_layers)
   endif ! was there any rain to infiltrate?

   ! if roots extent down into the bucket 
   if (root_reach > layer_thickness(1)) then
      ! how much has root depth extended since last step?
      depth_change = root_reach - previous_depth
      ! if there has been an increase
      if (depth_change > dble_zero) then
          ! calculate weighting between current lowest root layer and new soil 
          depth_change = depth_change / (depth_change+layer_thickness(nos_root_layers))
          ! add to bottom root layer
          soil_waterfrac(nos_root_layers) = (soil_waterfrac(nos_root_layers)*(dble_one-depth_change)) & 
                                          + (soil_waterfrac(nos_soil_layers)*depth_change)
      else
          ! calculate weighting between bottom soil layer and new bit coming from lowest root 
          depth_change = abs(depth_change) / (abs(depth_change)+layer_thickness(nos_root_layers))
          ! and add back to the bottom soil layer
          soil_waterfrac(nos_soil_layers) = (soil_waterfrac(nos_soil_layers)*(dble_one-depth_change)) &
                                          + (soil_waterfrac(nos_root_layers)*depth_change)
      end if ! depth change 

   end if ! root reach beyond top layer
   ! update new soil states
   previous_depth = root_reach
   ! determine soil layer thickness
   layer_thickness(1) = top_soil_depth ; layer_thickness(2) = max(0.1,root_reach-layer_thickness(1))
   layer_thickness(3) = max_depth - sum(layer_thickness(1:2))
   ! finally update soil water potential
   call soil_water_potential

!   ! sanity check for catastrophic failure
!   do soil_layer = 1, nos_soil_layers
!      if (soil_waterfrac(soil_layer) < 0.0 .and. soil_waterfrac(soil_layer) > -0.01) then
!          soil_waterfrac(soil_layer) = 0.0
!      endif
!      if (soil_waterfrac(soil_layer) < 0.0 .or. soil_waterfrac(soil_layer) /= soil_waterfrac(soil_layer)) then
!         print*,'ET',ET,"rainfall",rainfall
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
    rs  = (root_length*pi)**(-0.5) 
    rs2 = log( rs * root_radius_1 ) / (two_pi*root_length*thickness*soilC)
    ! soil water resistance
    soil_resistance = rs2*1e-9*mol_to_g_water

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
    water_retention_saxton_eqns = -1.0 * soil_wp + 10.0    ! 10 kPa represents air-entry swp

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
    double precision,parameter :: EPS = 3e-8

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
       tol1 = 2.0 * EPS * abs(b) + 0.5 * tol
       xm   = 0.5 * ( c - b )
       if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0d0 ) ) then
          zbrent = b
          return
       end if
       if ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) then
          s = fb / fa
          if ( a .eq. c ) then
             p = 2.0 * xm * s
             q = 1.0 - s
          else
             q = fa / fc
             r = fb / fc
             p = s * ( 2.0 * xm * q * ( q - r ) - ( b - a ) * ( r - 1.0 ) )
             q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 )
          end if
          if ( p .gt. 0.0 ) q = -q
          p = abs( p )
          if ( (2.0*p) .lt. min( 3.0*xm*q-abs(tol1*q) , abs(e*q) ) ) then
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
