
subroutine racmgppet(output_dim,met,pars,out_var,lat  &
                    ,nopars,nomet,nofluxes,nopools    &
                    ,nodays,deltat,nos_iter           &
                    ,soil_frac_clay_in,soil_frac_sand_in)

  use CARBON_MODEL_MOD, only: CARBON_MODEL, wSWP_time, nos_soil_layers, &
                              soil_frac_clay,soil_frac_sand

  ! subroutine specificially deals with the calling of the fortran code model by
  ! R

  implicit none

  ! declare input variables
  integer, intent(in) :: nopars         & ! number of paremeters in vector
                        ,output_dim     & !
                        ,nos_iter       & !
                        ,nomet          & ! number of meteorological fields
                        ,nofluxes       & ! number of model fluxes
                        ,nopools        & ! number of model pools
                        ,nodays           ! number of days in simulation

  double precision, intent(in) :: met(nomet,nodays)   & ! met drivers, note reverse of needed
                       ,pars(nopars,nos_iter)         & ! number of parameters
                       ,soil_frac_clay_in(nos_soil_layers) & ! clay in soil (%)
                       ,soil_frac_sand_in(nos_soil_layers) & ! sand in soil (%)
                       ,lat                 ! site latitude (degrees)

  double precision, intent(inout) :: deltat(nodays) ! time step in decimal days

  ! output declaration
  double precision, intent(out), dimension(nos_iter,nodays,output_dim) :: out_var

  ! local variables
  integer i
  ! vector of ecosystem pools
  double precision, dimension((nodays+1),nopools) :: POOLS
  ! vector of ecosystem fluxes
  double precision, dimension(nodays,nofluxes) :: FLUXES

  ! zero initial conditions
  POOLS = 0d0 ; FLUXES = 0d0
  out_var = 0d0

  ! update soil parameters
  soil_frac_clay=soil_frac_clay_in
  soil_frac_sand=soil_frac_sand_in

  ! generate deltat step from input data
  deltat(1) = met(1,1)
  do i = 2, nodays
     deltat(i)=met(1,i)-met(1,(i-1))
  end do

  ! begin iterations
  do i = 1, nos_iter

     ! reset at the beginning of each iteration
     POOLS = 0d0 ; FLUXES = 0d0

     ! call the models
     call CARBON_MODEL(1,nodays,met,pars(1:nopars,i),deltat,nodays &
                      ,lat,FLUXES,POOLS,nopars,nomet,nopools,nofluxes)

!if (i == 1) then
!    open(unit=666,file="/home/lsmallma/out.csv", &
!         status='replace',action='readwrite' )
!write(666,*)"deltat",deltat
!    write(666,*),"GSI",FLUXES(:,14)(1:365)
!    close(666)
!endif

     ! now allocate the output the our 'output' variable
     out_var(i,1:nodays,1)  = met(11,1:nodays)    ! LAI output for consistency check
     out_var(i,1:nodays,2)  = FLUXES(1:nodays,1)  ! GPP (gC.m-2.day-1)
     out_var(i,1:nodays,3)  = FLUXES(1:nodays,2)  ! transpiration (kgH2O.m-2.day-1)
     out_var(i,1:nodays,4)  = FLUXES(1:nodays,3)  ! wet canopy evaporation (kgH2O.m-2.day-1)
     out_var(i,1:nodays,5)  = FLUXES(1:nodays,4)  ! soil evaporation (kgH2O.m-2.day-1)
     out_var(i,1:nodays,6)  = wSWP_time(1:nodays) ! weighted soil water potential (MPa)
     out_var(i,1:nodays,7)  = POOLS(1:nodays,1)   ! Water in rooting zone (mm)
     out_var(i,1:nodays,8)  = FLUXES(1:nodays,5)  ! runoff (kgH2O.m-2.day-1)
     out_var(i,1:nodays,9)  = FLUXES(1:nodays,6)  ! drainage / underflow (kgH2O.m-2.day-1)
     out_var(i,1:nodays,10)  = FLUXES(1:nodays,7) ! estimate mean LWP (MPa)
     out_var(i,1:nodays,11)  = FLUXES(1:nodays,8) ! internal leaf CO2 concentration (umol/mol)

  end do ! nos_iter loop

  ! return back to the subroutine then
  return

end subroutine racmgppet
