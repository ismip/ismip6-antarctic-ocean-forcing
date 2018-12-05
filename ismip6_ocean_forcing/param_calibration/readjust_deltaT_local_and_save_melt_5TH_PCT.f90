program melting

!*************************************************************************************************************
!
! N. Jourdain, IGE-CNRS, Grenoble, Oct. 2018
!
! Program to test the non-local melting parameterization
!
! Will print values on screen and create melt_rates.nc
! 
! To use this script (here with the gfortran compiler) :
!    alias gfortran='gfortran -ffixed-line-length-none -ffree-line-length-none' # to read long lines
!    export NC_INC='-I /usr/local/include'                                      # platform-depedent 
!    export NC_LIB='-L /usr/local/lib -lnetcdf -lnetcdff'                       # platform-depedent
!    gfortran -c $NC_INC example_melt_NON_LOCAL_param.f90 
!    gfortran -o run_example example_melt_NON_LOCAL_param.o $NC_LIB
!    ./run_example 
!
!
!*************************************************************************************************************

USE netcdf

IMPLICIT NONE
 
INTEGER :: fidA, fidB, fidC, status, dimID_nbounds, dimID_z, dimID_y, dimID_x, mnbounds, mz, my, mx, z_ID,   &
&          thermal_forcing_ID, lon_ID, lat_ID, y_ID, x_ID, basinNumber_ID, icemask_shelves_ID, ice_draft_ID, &
&          thickness_ID, surface_ID, fidM, fidK, dimID_pct, kinf, ksup, ii, jj, kk, kbasin, Nbasin, mpct,    &
&          melt_ID, DeltaT_basin_pc05_ID, DeltaT_basin_pc50_ID, DeltaT_basin_pc95_ID, K0_ID, iter
 
CHARACTER(LEN=150) :: file_TF, file_basinNumbers, file_K, file_topo, file_melt_out

INTEGER*4, ALLOCATABLE, DIMENSION(:,:) :: basinNumber
 
REAL*8,ALLOCATABLE,DIMENSION(:) :: z, y, x, total, control, mean_TF, IS_area, K0, dT, dT0
 
REAL*8,ALLOCATABLE,DIMENSION(:,:) :: icemask_shelves, thickness, surface, ice_draft, TF_draft, &
&                                    mesh_area, melt
 
REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: thermal_forcing, DeltaT_basin_pc05, DeltaT_basin_pc50, DeltaT_basin_pc95

REAL*8 :: rhoi_SI, rhosw_SI, rhofw_SI, Lf_SI, cpw_SI, cste

LOGICAL :: ll_z_upward

!---------------------------------------------------------------------
! Constants

rhoi_SI   =  918.d0             ! Ice density (kg/m^3)
rhosw_SI  = 1028.d0             ! Sea water density (kg/m^3)
rhofw_SI  = 1000.d0             ! Fresh water density (kg/m^3)
Lf_SI     =    3.34d5           ! Fusion Latent heat of Ice (J/kg)
cpw_SI    = 3974.d0             ! Specific heat of sea water (J/kg/K)

cste = (rhosw_SI*cpw_SI/(rhoi_SI*Lf_SI))**2  ! in K^(-2)

!---------------------------------------------------------------------
! input files :

file_TF  = 'obs_thermal_forcing_1995-2017_8km_x_60m.nc'   ! 3d thermal forcing 
file_basinNumbers  = 'basinNumbers_8km.nc'                ! basins ID
file_K  = 'coeff_K0_DeltaT_quadratic_local_Rignot_Depoorter_with_TS_uncert.nc'
file_topo  = 'bedmap2_8km.nc'                            ! bedmap2 ice-shelf draft and mask

! output file containing melt rates :
file_melt_out = 'melt_rates_LOCAL_K0_DeltaT.nc'
 
!---------------------------------------------------------------------
! Read thermal forcing and dimensions
 
write(*,*) 'Reading ', TRIM(file_TF)
 
status = NF90_OPEN(TRIM(file_TF),0,fidA); call erreur(status,.TRUE.,"read TF")
 
status = NF90_INQ_DIMID(fidA,"z",dimID_z); call erreur(status,.TRUE.,"inq_dimID_z")
status = NF90_INQ_DIMID(fidA,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidA,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
 
status = NF90_INQUIRE_DIMENSION(fidA,dimID_z,len=mz); call erreur(status,.TRUE.,"inq_dim_z")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=my); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mx); call erreur(status,.TRUE.,"inq_dim_x")
  
ALLOCATE(  z(mz)  ) 
ALLOCATE(  thermal_forcing(my,mx,mz)  ) 
ALLOCATE(  y(my)  ) 
ALLOCATE(  x(mx)  ) 
 
status = NF90_INQ_VARID(fidA,"z",z_ID); call erreur(status,.TRUE.,"inq_z_ID")
status = NF90_INQ_VARID(fidA,"thermal_forcing",thermal_forcing_ID); call erreur(status,.TRUE.,"inq_thermal_forcing_ID")
status = NF90_INQ_VARID(fidA,"y",y_ID); call erreur(status,.TRUE.,"inq_y_ID")
status = NF90_INQ_VARID(fidA,"x",x_ID); call erreur(status,.TRUE.,"inq_x_ID")
 
status = NF90_GET_VAR(fidA,z_ID,z); call erreur(status,.TRUE.,"getvar_z")
status = NF90_GET_VAR(fidA,thermal_forcing_ID,thermal_forcing); call erreur(status,.TRUE.,"getvar_thermal_forcing")
status = NF90_GET_VAR(fidA,y_ID,y); call erreur(status,.TRUE.,"getvar_y")
status = NF90_GET_VAR(fidA,x_ID,x); call erreur(status,.TRUE.,"getvar_x")
 
status = NF90_CLOSE(fidA); call erreur(status,.TRUE.,"close_file")

! To adapt to the ISM grid :
ALLOCATE ( mesh_area(my,mx)  )
mesh_area(:,:) = abs( (x(3)-x(2)) * (y(3)-y(2)) )

!-------------------------------------------------------------------
! Read basin numbers
 
write(*,*) 'Reading ', TRIM(file_basinNumbers)
status = NF90_OPEN(TRIM(file_basinNumbers),0,fidB); call erreur(status,.TRUE.,"read basinNumbers")
ALLOCATE(  basinNumber(mx,my)  ) 
status = NF90_INQ_VARID(fidB,"basinNumber",basinNumber_ID); call erreur(status,.TRUE.,"inq_basinNumber_ID")
status = NF90_GET_VAR(fidB,basinNumber_ID,basinNumber); call erreur(status,.TRUE.,"getvar_basinNumber")
status = NF90_CLOSE(fidB); call erreur(status,.TRUE.,"close_file")

!------------------------------------------------------------------
! Read K coefficient for melting parameterization

write(*,*) 'Reading ', TRIM(file_K) 
status = NF90_OPEN(TRIM(file_K),0,fidK); call erreur(status,.TRUE.,"read")
status = NF90_INQ_DIMID(fidK,"pct",dimID_pct); call erreur(status,.TRUE.,"inq_dimID_pct")
status = NF90_INQUIRE_DIMENSION(fidK,dimID_pct,len=mpct); call erreur(status,.TRUE.,"inq_dim_pct")
ALLOCATE(  DeltaT_basin_pc95(mx,my,mpct)  ) 
ALLOCATE(  DeltaT_basin_pc50(mx,my,mpct)  ) 
ALLOCATE(  DeltaT_basin_pc05(mx,my,mpct)  ) 
ALLOCATE(  K0(mpct)  ) 
status = NF90_INQ_VARID(fidK,"DeltaT_basin_pc95",DeltaT_basin_pc95_ID); call erreur(status,.TRUE.,"inq_DeltaT_basin_pc95_ID")
status = NF90_INQ_VARID(fidK,"DeltaT_basin_pc50",DeltaT_basin_pc50_ID); call erreur(status,.TRUE.,"inq_DeltaT_basin_pc50_ID")
status = NF90_INQ_VARID(fidK,"DeltaT_basin_pc05",DeltaT_basin_pc05_ID); call erreur(status,.TRUE.,"inq_DeltaT_basin_pc05_ID")
status = NF90_INQ_VARID(fidK,"K0",K0_ID); call erreur(status,.TRUE.,"inq_K0_ID")
status = NF90_GET_VAR(fidK,DeltaT_basin_pc95_ID,DeltaT_basin_pc95); call erreur(status,.TRUE.,"getvar_DeltaT_basin_pc95")
status = NF90_GET_VAR(fidK,DeltaT_basin_pc50_ID,DeltaT_basin_pc50); call erreur(status,.TRUE.,"getvar_DeltaT_basin_pc50")
status = NF90_GET_VAR(fidK,DeltaT_basin_pc05_ID,DeltaT_basin_pc05); call erreur(status,.TRUE.,"getvar_DeltaT_basin_pc05")
status = NF90_GET_VAR(fidK,K0_ID,K0); call erreur(status,.TRUE.,"getvar_K0")
status = NF90_CLOSE(fidK); call erreur(status,.TRUE.,"close_file")

!-----------------------------------------------------------------------------------------------------
! Read BEDMAP2 ice-shelf draft and mask (this part should eventually come from the ice-sheet model)
 
write(*,*) 'Reading ', TRIM(file_topo)
 
status = NF90_OPEN(TRIM(file_topo),0,fidC); call erreur(status,.TRUE.,"read")
 
ALLOCATE(  icemask_shelves(mx,my)  ) 
ALLOCATE(  thickness(mx,my)  ) 
ALLOCATE(  surface(mx,my)  ) 
ALLOCATE(  ice_draft(mx,my)  ) 
 
status = NF90_INQ_VARID(fidC,"icemask_shelves",icemask_shelves_ID); call erreur(status,.TRUE.,"inq_icemask_shelves_ID")
status = NF90_INQ_VARID(fidC,"thickness",thickness_ID); call erreur(status,.TRUE.,"inq_thickness_ID")
status = NF90_INQ_VARID(fidC,"surface",surface_ID); call erreur(status,.TRUE.,"inq_surface_ID")
 
status = NF90_GET_VAR(fidC,icemask_shelves_ID,icemask_shelves); call erreur(status,.TRUE.,"getvar_icemask_shelves")
status = NF90_GET_VAR(fidC,thickness_ID,thickness); call erreur(status,.TRUE.,"getvar_thickness")
status = NF90_GET_VAR(fidC,surface_ID,surface); call erreur(status,.TRUE.,"getvar_surface")
 
status = NF90_CLOSE(fidC); call erreur(status,.TRUE.,"close_file")

ice_draft(:,:) = surface(:,:) - thickness(:,:)
DEALLOCATE( surface, thickness )

!-------------------------------------------------------------------
! Melt rate calculation

Nbasin=MAXVAL(basinNumber)+1 ! number of basins

ALLOCATE( TF_draft(mx,my), melt(mx,my)  )
ALLOCATE( mean_TF(Nbasin), IS_area(Nbasin)  )

mean_TF(:) = 0.d0
IS_area(:) = 0.d0

if ( z(2) .gt. z(1) ) then
  ll_z_upward = .true.
else
  ll_z_upward = .false.
endif

do ii=1,mx
do jj=1,my

  if ( icemask_shelves(ii,jj) .ge. 5.d-1 .and. ice_draft(ii,jj) .lt. 0.d0 ) then

    ! 1 -  Linear interpolation of the thermal forcing on the ice draft depth :
    if ( ll_z_upward ) then
      kinf=1
      do kk=2,mz-1
        if ( z(kk) .le. ice_draft(ii,jj) ) kinf = kk
      enddo
      ksup = kinf + 1
    else
      ksup=mz
      do kk=mz-1,2,-1
        if ( z(kk) .le. ice_draft(ii,jj) ) ksup = kk
      enddo
      kinf = ksup - 1
    endif

    TF_draft(ii,jj) = (   (z(ksup)-ice_draft(ii,jj)) * thermal_forcing(ii,jj,kinf)   &
    &                   + (ice_draft(ii,jj)-z(kinf)) * thermal_forcing(ii,jj,ksup) ) / (z(ksup)-z(kinf))

    ! 2 -  Mean Thermal forcing in individual basins (NB: fortran norm while basins start at zero) :
    mean_TF( basinNumber(ii,jj)+1 ) = mean_TF( basinNumber(ii,jj)+1 ) + mesh_area(ii,jj) * TF_draft(ii,jj)
    IS_area( basinNumber(ii,jj)+1 ) = IS_area( basinNumber(ii,jj)+1 ) + mesh_area(ii,jj)

  else

    TF_draft(ii,jj) = -9999.9d0

  endif

enddo
enddo

mean_TF(:) = mean_TF(:) / IS_area(:) ! assuming that there are ice shelves in all basins

! 3  - Calculation of melting rate :

! melt rate in m/yr (meters of pure water per year) :
! [ * rhofw_SI / rhoi_SI to get it in meters of ice per year ]
do ii=1,mx
do jj=1,my
  if ( TF_draft(ii,jj) .gt. -5.d0 ) then
    melt(ii,jj) = K0(1) * cste * max( TF_draft(ii,jj) + DeltaT_basin_pc50(ii,jj,1), 0.d0 )* max( TF_draft(ii,jj) + DeltaT_basin_pc50(ii,jj,1), 0.d0 )
  else
    melt(ii,jj) = 0.d0
  endif
enddo
enddo

!-------------------------------------------------------------------
! Testing total melt rates per basin 
! (so that you can check while you adapt this program to your ice sheet model):

ALLOCATE( dT(Nbasin), dT0(Nbasin), total(Nbasin), control(Nbasin) )

! obtained from true_median in calculate_K0_DeltaT_quadratic.f90 :
control(:) = (/ 57.81, 30.99, 37.24, 111.55, 111.76, 21.36, 17.04, 65.06, 145.44, 252.76, 89.39, 157.00, 22.42, 4.21, 99.27, 22.72 /)
dT(:) = (/ -0.242, 1.026, -0.016, 0.578, 1.000, 0.801, 0.273, -0.210, 0.625, 1.178, -0.459, -0.784, 0.486, -0.027, -0.024, 0.206 /)

dT0(:) = dT(:)

total(:) = 0.d0
do ii=1,mx
do jj=1,my
  total( basinNumber(ii,jj)+1 ) = total( basinNumber(ii,jj)+1 ) + melt(ii,jj) * mesh_area(ii,jj) * 1.d-9  ! in Gt/yr
enddo
enddo

DO iter=1,10000

  do kbasin=1,Nbasin
    if     ( total(kbasin)-control(kbasin) .gt.  0.001*control(kbasin) ) then
      dT(kbasin) = dT(kbasin)-0.001
    elseif ( total(kbasin)-control(kbasin) .lt. -0.001*control(kbasin) ) then
      dT(kbasin) = dT(kbasin)+0.001 
    endif
  enddo

  total(:) = 0.d0
  do ii=1,mx
  do jj=1,my
    if ( TF_draft(ii,jj)+dT(basinNumber(ii,jj)+1) .gt. -5.d0 ) then
      melt(ii,jj) = K0(1) * cste * max( TF_draft(ii,jj) + dT(basinNumber(ii,jj)+1), 0.d0 )* max( TF_draft(ii,jj) + dT(basinNumber(ii,jj)+1), 0.d0 )
    else
      melt(ii,jj) = 0.d0
    endif
    total( basinNumber(ii,jj)+1 ) = total( basinNumber(ii,jj)+1 ) + melt(ii,jj) * mesh_area(ii,jj) * 1.d-9  ! in Gt/yr
  enddo
  enddo

ENDDO

do kbasin=1,Nbasin
  write(*,883) kbasin, dT0(kbasin), dT(kbasin)
enddo
883 FORMAT('deltaT in basin ', i2.2,' modified from ', f7.3, ' to ', f7.3)

write(*,*) '----------------------------------------------------------'
write(*,*) dT(:)
write(*,*) '----------------------------------------------------------'

do kbasin=1,Nbasin
  write(*,100) kbasin, total(kbasin), control(kbasin)
enddo
100 FORMAT(' melt in basin ',i2.2,' is ',f9.2,' Gt/yr   [ Ctrl value = ',f9.2,' ]')

write(*,*) ' Total Antarctic ice-shelf basal mass balance = ', sum(total), ' Gt/yr   [ Ctrl value = 1533.82 ]'

!-------------------------------------------------------------------
! Save melt rate in a netcdf file :

write(*,*) 'Creating ', TRIM(file_melt_out)
 
status = NF90_CREATE(TRIM(file_melt_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidM,"y",my,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"x",mx,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")
  
status = NF90_DEF_VAR(fidM,"y",NF90_DOUBLE,(/dimID_y/),y_ID); call erreur(status,.TRUE.,"def_var_y_ID")
status = NF90_DEF_VAR(fidM,"x",NF90_DOUBLE,(/dimID_x/),x_ID); call erreur(status,.TRUE.,"def_var_x_ID")
status = NF90_DEF_VAR(fidM,"melt_rate",NF90_DOUBLE,(/dimID_x,dimID_y/),melt_ID); call erreur(status,.TRUE.,"def_var_melt_basinK_mean_ID")

status = NF90_PUT_ATT(fidM,melt_ID,"long_name","melt rate"); call erreur(status,.TRUE.,"put_att_melt_basinK_mean_ID")
status = NF90_PUT_ATT(fidM,melt_ID,"units","m.w.e / yr"); call erreur(status,.TRUE.,"put_att_melt_basinK_mean_ID")
 
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using readjust_deltaT_local_and_save_melt.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"K coefficient file",TRIM(file_K)); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 
 
status = NF90_PUT_VAR(fidM,y_ID,y); call erreur(status,.TRUE.,"var_y_ID")
status = NF90_PUT_VAR(fidM,x_ID,x); call erreur(status,.TRUE.,"var_x_ID")
status = NF90_PUT_VAR(fidM,melt_ID,melt); call erreur(status,.TRUE.,"var_melt_basinK_mean_ID")
 
status = NF90_CLOSE(fidM); call erreur(status,.TRUE.,"final")

!------------------

end program melting


!==========================================================

SUBROUTINE erreur(iret, lstop, chaine)
  ! pour les messages d'erreur
  USE netcdf
  INTEGER, INTENT(in)                     :: iret
  LOGICAL, INTENT(in)                     :: lstop
  CHARACTER(LEN=*), INTENT(in)            :: chaine
  !
  CHARACTER(LEN=80)                       :: message
  !
  IF ( iret .NE. 0 ) THEN
    WRITE(*,*) 'ROUTINE: ', TRIM(chaine)
    WRITE(*,*) 'ERREUR: ', iret
    message=NF90_STRERROR(iret)
    WRITE(*,*) 'CA VEUT DIRE:',TRIM(message)
    IF ( lstop ) STOP
  ENDIF
  !
END SUBROUTINE erreur
