program modif

USE netcdf
 
IMPLICIT NONE
 
INTEGER :: fidISFMSK, fidBASIN, status, dimID_y, dimID_x, my, mx, y_ID, x_ID, basinNumber_ID, dimID_Nisf, Nisf, IF_mask_ID, GL_mask_ID, isfmask_ID, &
&          err_melt_isf_ID, melt_isf_ID, front_max_lat_ID, front_min_lat_ID, front_max_lon_ID, front_min_lon_ID, front_ice_dep_avg_ID, fidK, k4, N4,&
&          front_ice_dep_min_ID, front_bot_dep_avg_ID, front_bot_dep_max_ID, lat_ID, lon_ID, name_reg_ID, name_isf_ID, fidM, mx2, my2, fidTF, kk,   &
&          TF_draft_ID, Nbasin, kbasin, Nstat, Nstat2, kstat, kisf, ii, jj, melt_ID, para, nn_data_melt, ipct, Npct, Nmerge

INTEGER :: melt_DeltaT_pc05_ID, melt_DeltaT_pc50_ID, melt_DeltaT_pc95_ID, K0_ID, pct_ID, dimID_pct, Nsmooth, NminR, NmaxR, NminD, NmaxD,    &
&          DeltaT_basin_pc05_ID, DeltaT_basin_pc50_ID, DeltaT_basin_pc95_ID
 
CHARACTER(LEN=150) :: file_basin, file_isf_mask, file_melt_out, file_TF, file_K_out

CHARACTER(LEN=50) :: str_para, str_data_melt, str_TS_uncertainty

INTEGER*1,ALLOCATABLE,DIMENSION(:) :: isthere, istwice

INTEGER*2,ALLOCATABLE,DIMENSION(:,:) :: isfmask, tmpmsk, tmpmsk2

INTEGER*4,ALLOCATABLE,DIMENSION(:) :: indx, indx2, NN

INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: basinNumber
 
REAL*8,ALLOCATABLE,DIMENSION(:) :: y, x, err_melt_isf, melt_isf, xmelt, Kfac, mean_melt_basin, std_melt_basin, median_melt_basin,  &
&                                  thermo, true_median, correc_med_D, correc_med_R, perct, K0_pct

REAL*8,ALLOCATABLE,DIMENSION(:,:) :: TF_draft, thermo2d, DeltaT, DeltaT_pc05, DeltaT_pc50, DeltaT_pc95, melt_basin_D, err_melt_basin_D

REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: DeltaT_basin_pc05, DeltaT_basin_pc50, DeltaT_basin_pc95, melt_DeltaT_pc05, melt_DeltaT_pc50, melt_DeltaT_pc95 

REAL*8 :: dx, dy, rr, r1, r2, pi, cste, rhoi_SI, rhosw_SI, rhofw_SI, Lf_SI, cpw_SI, aa, check_median, err_TF, meanTF, meanTF2, &
&         areaTF, dis, check_max_DeltaT_pc05, check_max_DeltaT_pc50, check_max_DeltaT_pc95

LOGICAL :: positive, ll_to_do, ll_upscaling, ll_TS_uncertainty

!---------------------------------------

file_isf_mask  = 'ice_shelf_mask_bedmap2_8km.nc'
file_basin     = 'basinNumbers_8km.nc'
file_TF        = 'obs_thermal_forcing_on_bedmap2_ice_draft.nc'

para = 1      ! =1 quadratic non-local :  TF.<TF>
              ! =2 quadratic local : TF.TF

nn_data_melt = 3       ! =1 Rignot et al. (2013)
                       ! =2 Depoorter et al. (2013)
                       ! =3 both Rignot et al. (2013) and Depoorter et al. (2013)

Nstat = 100000 !! number of samples to sample observation-based melt rates

Nsmooth = 0  !! size of the running window (running mean between -Nsmooth/2 and +Nsmooth/2)
             !! No smoothing if Nsmooth < 2

ll_upscaling = .false.  ! =.false. if we consider that our resolution only represents the monitored
                        !          glaciers (>~100km2) [ false recommended for 8km resolution] 

ll_TS_uncertainty = .true. ! =.true. to take into account uncertainties on TF from WOA+MEOP

! Basins to merge with kbasin-1 (e.g. NN=5 if 5 and 4 are merged)
Nmerge=0 ! =0 if no basins to merge
if ( Nmerge .gt. 0 ) then
  ALLOCATE( NN(Nmerge) )
  NN=(/ 13 /) ! python numbering (same as basin number in file_basin)
endif

pi = acos(-1.d0)

rhoi_SI   =  918.d0             ! Ice density (kg/m^3)
rhosw_SI  = 1028.d0             ! Sea water density (kg/m^3)
rhofw_SI  = 1000.d0             ! Fresh water density (kg/m^3)
Lf_SI     =    3.34d5           ! Fusion Latent heat of Ice (J/kg)
cpw_SI    = 3974.d0             ! Specific heat of sea water (J/kg/K)

cste = (rhosw_SI*cpw_SI/(rhoi_SI*Lf_SI))**2  ! in K^(-2)
 
CALL RANDOM_SEED

!---------------------------------------
! name of output files :

if     ( para .eq. 1 ) then
  str_para='quadratic_non_local'
elseif ( para .eq. 2 ) then
  str_para='quadratic_local    '
endif

if     ( nn_data_melt .eq. 1 ) then
  str_data_melt='Rignot          '
elseif ( nn_data_melt .eq. 2 ) then
  str_data_melt='Depoorter       '
elseif ( nn_data_melt .eq. 3 ) then
  str_data_melt='Rignot_Depoorter'
endif

if (ll_TS_uncertainty) then
  str_TS_uncertainty='with_TS_uncert'
else
  str_TS_uncertainty='no_TS_uncert  '
endif

write(file_melt_out,751) TRIM(str_para), TRIM(str_data_melt), TRIM(str_TS_uncertainty)
write(file_K_out,752) TRIM(str_para), TRIM(str_data_melt), TRIM(str_TS_uncertainty)

751 FORMAT('melt_K0_DeltaT_',a,'_',a,'_',a,'.nc')
752 FORMAT('coeff_K0_DeltaT_',a,'_',a,'_',a,'.nc')

write(*,*) TRIM(file_melt_out)
write(*,*) TRIM(file_K_out)

!---------------------------------------------------------------------
! Uncertainty (stddev) on temperature profiles for WOA+MEOP.
! ( calculated as the mean standard deviation between 80S and 60S
!   only considering points with more than 3 valid points (TT_DD)
!   and assuming that the uncertainty on TF comes from T and not from S )
! 
! dep_TF_err( 1) = 0    ; TF_err( 1) = 0.8079
! dep_TF_err( 2) = 5    ; TF_err( 2) = 0.8103
! dep_TF_err( 3) = 10   ; TF_err( 3) = 0.8046
! dep_TF_err( 4) = 15   ; TF_err( 4) = 0.7967
! dep_TF_err( 5) = 20   ; TF_err( 5) = 0.7871
! dep_TF_err( 6) = 25   ; TF_err( 6) = 0.7774
! dep_TF_err( 7) = 30   ; TF_err( 7) = 0.7675
! dep_TF_err( 8) = 35   ; TF_err( 8) = 0.7473
! dep_TF_err( 9) = 40   ; TF_err( 9) = 0.7235
! dep_TF_err(10) = 45   ; TF_err(10) = 0.6968
! dep_TF_err(11) = 50   ; TF_err(11) = 0.6660
! dep_TF_err(12) = 55   ; TF_err(12) = 0.6271
! dep_TF_err(13) = 60   ; TF_err(13) = 0.5900
! dep_TF_err(14) = 65   ; TF_err(14) = 0.5575
! dep_TF_err(15) = 70   ; TF_err(15) = 0.5342
! dep_TF_err(16) = 75   ; TF_err(16) = 0.5168
! dep_TF_err(17) = 80   ; TF_err(17) = 0.5035
! dep_TF_err(18) = 85   ; TF_err(18) = 0.4951
! dep_TF_err(19) = 90   ; TF_err(19) = 0.4911
! dep_TF_err(20) = 95   ; TF_err(20) = 0.4889
! dep_TF_err(21) = 100  ; TF_err(21) = 0.4912
! dep_TF_err(22) = 125  ; TF_err(22) = 0.4955
! dep_TF_err(23) = 150  ; TF_err(23) = 0.4619
! dep_TF_err(24) = 175  ; TF_err(24) = 0.4046
! dep_TF_err(25) = 200  ; TF_err(25) = 0.3556
! dep_TF_err(26) = 225  ; TF_err(26) = 0.3219
! dep_TF_err(27) = 250  ; TF_err(27) = 0.2881
! dep_TF_err(28) = 275  ; TF_err(28) = 0.2472
! dep_TF_err(29) = 300  ; TF_err(29) = 0.2337
! dep_TF_err(30) = 325  ; TF_err(30) = 0.2250
! dep_TF_err(31) = 350  ; TF_err(31) = 0.2158
! dep_TF_err(32) = 375  ; TF_err(32) = 0.2089
! dep_TF_err(33) = 400  ; TF_err(33) = 0.2000
! dep_TF_err(34) = 425  ; TF_err(34) = 0.1930
! dep_TF_err(35) = 450  ; TF_err(35) = 0.1861
! dep_TF_err(36) = 475  ; TF_err(36) = 0.1792
! dep_TF_err(37) = 500  ; TF_err(37) = 0.1687
! dep_TF_err(38) = 550  ; TF_err(38) = 0.1571
! dep_TF_err(39) = 600  ; TF_err(39) = 0.1524
! dep_TF_err(40) = 650  ; TF_err(40) = 0.1464
! dep_TF_err(41) = 700  ; TF_err(41) = 0.1441
! dep_TF_err(42) = 750  ; TF_err(42) = 0.1396
! dep_TF_err(43) = 800  ; TF_err(43) = 0.1348
! dep_TF_err(44) = 850  ; TF_err(44) = 0.1304
! dep_TF_err(45) = 900  ; TF_err(45) = 0.1200
! dep_TF_err(46) = 950  ; TF_err(46) = 0.0986
! dep_TF_err(47) = 1000 ; TF_err(47) = 0.0962
! dep_TF_err(48) = 1050 ; TF_err(48) = 0.0869
! dep_TF_err(49) = 1100 ; TF_err(49) = 0.0836
! dep_TF_err(50) = 1150 ; TF_err(50) = 0.0808
! dep_TF_err(51) = 1200 ; TF_err(51) = 0.0766
! dep_TF_err(52) = 1250 ; TF_err(52) = 0.0752
! dep_TF_err(53) = 1300 ; TF_err(53) = 0.0724
! dep_TF_err(54) = 1350 ; TF_err(54) = 0.0707
! dep_TF_err(55) = 1400 ; TF_err(55) = 0.0685
! dep_TF_err(56) = 1450 ; TF_err(56) = 0.0642
! dep_TF_err(57) = 1500 ; TF_err(57) = 0.0614
! dep_TF_err(58) = 1550 ; TF_err(58) = 0.0615
! dep_TF_err(59) = 1600 ; TF_err(59) = 0.0593
! dep_TF_err(60) = 1650 ; TF_err(60) = 0.0565
! dep_TF_err(61) = 1700 ; TF_err(61) = 0.0543
! dep_TF_err(62) = 1750 ; TF_err(62) = 0.0543
! dep_TF_err(63) = 1800 ; TF_err(63) = 0.0518
! dep_TF_err(64) = 1850 ; TF_err(64) = 0.0503
! dep_TF_err(65) = 1900 ; TF_err(65) = 0.0480
! dep_TF_err(66) = 1950 ; TF_err(66) = 0.0445
! dep_TF_err(67) = 2000 ; TF_err(67) = 0.0595
! dep_TF_err(68) = 2100 ; TF_err(68) = 0.0577
! dep_TF_err(69) = 2200 ; TF_err(69) = 0.0556
! dep_TF_err(70) = 2300 ; TF_err(70) = 0.0548
! dep_TF_err(71) = 2400 ; TF_err(71) = 0.0555
! dep_TF_err(72) = 2500 ; TF_err(72) = 0.0546
! dep_TF_err(73) = 2600 ; TF_err(73) = 0.0517
! dep_TF_err(74) = 2700 ; TF_err(74) = 0.0490
! dep_TF_err(75) = 2800 ; TF_err(75) = 0.0485
! dep_TF_err(76) = 2900 ; TF_err(76) = 0.0442
! dep_TF_err(77) = 3000 ; TF_err(77) = 0.0378
! dep_TF_err(78) = 3100 ; TF_err(78) = 0.0380
! dep_TF_err(79) = 3200 ; TF_err(79) = 0.0383
! dep_TF_err(80) = 3300 ; TF_err(80) = 0.0272
! dep_TF_err(81) = 3400 ; TF_err(81) = 0.0269
! dep_TF_err(82) = 3500 ; TF_err(82) = 0.0229
! dep_TF_err(83) = 3600 ; TF_err(83) = 0.0215
! dep_TF_err(84) = 3700 ; TF_err(84) = 0.0209
! dep_TF_err(85) = 3800 ; TF_err(85) = 0.0173
! dep_TF_err(86) = 3900 ; TF_err(86) = 0.0155
! dep_TF_err(87) = 4000 ; TF_err(87) = 0.0161

! Taking mean stddev at 500m depth and apply it everywhere :
! (Taking locally would cancel the uncertainty over large areas
!  so better to apply coherent error )

err_TF = 0.1687 

!---------------------------------------
! Read basin numbers :
 
write(*,*) 'Reading ', TRIM(file_basin)
 
status = NF90_OPEN(TRIM(file_basin),0,fidBASIN); call erreur(status,.TRUE.,"read")
 
status = NF90_INQ_DIMID(fidBASIN,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidBASIN,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
 
status = NF90_INQUIRE_DIMENSION(fidBASIN,dimID_y,len=my); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidBASIN,dimID_x,len=mx); call erreur(status,.TRUE.,"inq_dim_x")
  
ALLOCATE(  y(my)  ) 
ALLOCATE(  x(mx)  ) 
ALLOCATE(  basinNumber(mx,my)  ) 
 
status = NF90_INQ_VARID(fidBASIN,"y",y_ID); call erreur(status,.TRUE.,"inq_y_ID")
status = NF90_INQ_VARID(fidBASIN,"x",x_ID); call erreur(status,.TRUE.,"inq_x_ID")
status = NF90_INQ_VARID(fidBASIN,"basinNumber",basinNumber_ID); call erreur(status,.TRUE.,"inq_basinNumber_ID")
 
status = NF90_GET_VAR(fidBASIN,y_ID,y); call erreur(status,.TRUE.,"getvar_y")
status = NF90_GET_VAR(fidBASIN,x_ID,x); call erreur(status,.TRUE.,"getvar_x")
status = NF90_GET_VAR(fidBASIN,basinNumber_ID,basinNumber); call erreur(status,.TRUE.,"getvar_basinNumber")
 
status = NF90_CLOSE(fidBASIN); call erreur(status,.TRUE.,"close_file")

Nbasin = maxval(basinNumber) + 1 ! python convention for basinNumber
dx = abs(x(3)-x(2)) 
dy = abs(y(3)-y(2)) 
write(*,*) '     Resolution is ', dx, dy
write(*,*) '     Number of basins : ', Nbasin
write(*,*) '     Number of basins to merge : ', Nmerge

! Merge some basins
do kk=1,Nmerge
  WHERE( basinNumber(:,:) .eq. NN(kk) )
    basinNumber(:,:) = NN(kk)-1
  ENDWHERE
enddo

!---------------------------------------
! Read ice shelf mask and melt rates :
 
write(*,*) 'Reading ', TRIM(file_isf_mask)
 
status = NF90_OPEN(TRIM(file_isf_mask),0,fidISFMSK); call erreur(status,.TRUE.,"read")
 
status = NF90_INQ_DIMID(fidISFMSK,"Nisf",dimID_Nisf); call erreur(status,.TRUE.,"inq_dimID_Nisf")
status = NF90_INQ_DIMID(fidISFMSK,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidISFMSK,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
 
status = NF90_INQUIRE_DIMENSION(fidISFMSK,dimID_Nisf,len=Nisf); call erreur(status,.TRUE.,"inq_dim_Nisf")
status = NF90_INQUIRE_DIMENSION(fidISFMSK,dimID_y,len=my2); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidISFMSK,dimID_x,len=mx2); call erreur(status,.TRUE.,"inq_dim_x")

if ( mx2 .ne. mx .or. my2 .ne. my ) then
  write(*,*) '~!@#$%^* ERROR : dimension mismatch >>> CHECK netcdf FILES >>>> STOP'
  stop
endif
  
ALLOCATE(  isfmask(mx,my)  ) 
ALLOCATE(  err_melt_isf(Nisf)  ) 
ALLOCATE(  melt_isf(Nisf)  ) 
 
status = NF90_INQ_VARID(fidISFMSK,"isfmask",isfmask_ID); call erreur(status,.TRUE.,"inq_isfmask_ID")
status = NF90_INQ_VARID(fidISFMSK,"err_melt_isf",err_melt_isf_ID); call erreur(status,.TRUE.,"inq_err_melt_isf_ID")
status = NF90_INQ_VARID(fidISFMSK,"melt_isf",melt_isf_ID); call erreur(status,.TRUE.,"inq_melt_isf_ID")
 
status = NF90_GET_VAR(fidISFMSK,isfmask_ID,isfmask); call erreur(status,.TRUE.,"getvar_isfmask")
status = NF90_GET_VAR(fidISFMSK,err_melt_isf_ID,err_melt_isf); call erreur(status,.TRUE.,"getvar_err_melt_isf")
status = NF90_GET_VAR(fidISFMSK,melt_isf_ID,melt_isf); call erreur(status,.TRUE.,"getvar_melt_isf")
 
status = NF90_CLOSE(fidISFMSK); call erreur(status,.TRUE.,"close_file")

! uniform upscaling for Rignot et al. (2013) if required :
IF ( ll_upscaling ) THEN
  melt_isf(:) = melt_isf(:) * 1500.d0 / 1325.d0
  err_melt_isf(:) = err_melt_isf(:) * 1500.d0 / 1325.d0
ENDIF

!-----------------------------------------------------------------
! Read thermal forcing (previously interpolated on the ice draft)
 
write(*,*) 'Reading ', TRIM(file_TF)
 
status = NF90_OPEN(TRIM(file_TF),0,fidTF); call erreur(status,.TRUE.,"read")
 
ALLOCATE(  TF_draft(my,mx)  ) 
 
status = NF90_INQ_VARID(fidTF,"TF_draft",TF_draft_ID); call erreur(status,.TRUE.,"inq_TF_draft_ID")
 
status = NF90_GET_VAR(fidTF,TF_draft_ID,TF_draft); call erreur(status,.TRUE.,"getvar_TF_draft")
 
status = NF90_CLOSE(fidTF); call erreur(status,.TRUE.,"close_file")

! remove NaNs :
do ii=1,mx
do jj=1,my
  if ( isnan(TF_draft(ii,jj)) ) TF_draft(ii,jj) = 0.d0
enddo
enddo

!-----------------------------------------------------------------
! Calculation...
 
N4=4

ALLOCATE( tmpmsk(mx,my), tmpmsk2(mx,my) )
ALLOCATE( isthere(Nisf), istwice(Nisf) )
ALLOCATE( xmelt(Nstat), Kfac(Nstat), thermo(Nstat) )
ALLOCATE( mean_melt_basin(Nbasin), std_melt_basin(Nbasin), median_melt_basin(Nbasin) )
ALLOCATE( correc_med_R(Nbasin), correc_med_D(Nbasin), true_median(Nbasin) )
ALLOCATE( melt_basin_D(Nbasin,N4), err_melt_basin_D(Nbasin,N4) )

! Melt values per basin (Depoorter et al. 2013) [2nd dimension for multiple groups of ice shelves]
melt_basin_D     (:,:) = 0.d0 
err_melt_basin_D (:,:) = 0.d0

melt_basin_D( 1,1) =  24.d0                   ; err_melt_basin_D( 1,1) =      13.d0                    ! K'A+AA'  - JF+AR
melt_basin_D( 1,2) =        39.d0             ; err_melt_basin_D( 1,2) =            16.d0              !
melt_basin_D( 2,1) =  36.d0                   ; err_melt_basin_D( 2,1) =      16.d0                    ! A'B  - NE
melt_basin_D( 3,1) =  39.d0                   ; err_melt_basin_D( 3,1) =      21.d0                    ! BC   - AIS
melt_basin_D( 4,1) =  26.d0                   ; err_melt_basin_D( 4,1) =      13.d0                    ! CC'  - W+SHA
melt_basin_D( 4,2) =        76.d0             ; err_melt_basin_D( 4,2) =            23.d0              !
melt_basin_D( 5,1) =   5.d0                   ; err_melt_basin_D( 5,1) =       4.d0                    ! C'D  - VAN+TOT+MU+POR
melt_basin_D( 5,2) =        64.d0             ; err_melt_basin_D( 5,2) =            12.d0              !
melt_basin_D( 5,3) =              28.d0       ; err_melt_basin_D( 5,3) =                  7.d0         !
melt_basin_D( 5,4) =                    18.d0 ; err_melt_basin_D( 5,4) =                       10.d0   !
melt_basin_D( 6,1) =  13.d0                   ; err_melt_basin_D( 6,1) =       4.d0                    ! DD'  - ADE+MER+NIN+COO
melt_basin_D( 6,2) =         5.d0             ; err_melt_basin_D( 6,2) =             4.d0              !
melt_basin_D( 6,3) =               0.d0       ; err_melt_basin_D( 6,3) =                  3.d0         !
melt_basin_D( 6,4) =                     3.d0 ; err_melt_basin_D( 6,4) =                       7.d0    !
melt_basin_D( 7,1) =   7.d0                   ; err_melt_basin_D( 7,1) =       2.d0                    ! D'E  - REN + DRY
melt_basin_D( 7,2) =         5.d0             ; err_melt_basin_D( 7,2) =             1.d0              !
melt_basin_D( 8,1) =  34.d0                   ; err_melt_basin_D( 8,1) =      25.d0                    ! EF'  - RIS+SUL
melt_basin_D( 8,2) =        28.d0             ; err_melt_basin_D( 8,2) =             6.d0              ! 
melt_basin_D( 9,1) =   6.d0                   ; err_melt_basin_D( 9,1) =       6.d0                    ! F'G  - LAN+GET
melt_basin_D( 9,2) =       136.d0             ; err_melt_basin_D( 9,2) =            23.d0              !
melt_basin_D(10,1) =  78.d0                   ; err_melt_basin_D(10,1) =       7.d0                    ! GH   - CD+THW+PI+COS
melt_basin_D(10,2) =        69.d0             ; err_melt_basin_D(10,2) =            18.d0              !
melt_basin_D(10,3) =              95.d0       ; err_melt_basin_D(10,3) =                  14.d0        !
melt_basin_D(10,4) =                    11.d0 ; err_melt_basin_D(10,4) =                        3.d0   !
melt_basin_D(11,1) =  86.d0                   ; err_melt_basin_D(11,1) =      22.d0                    ! HH'  - ABB+VEN
melt_basin_D(11,2) =        15.d0             ; err_melt_basin_D(11,2) =            3.d0               !
melt_basin_D(12,1) = 144.d0                   ; err_melt_basin_D(12,1) =      42.d0                    ! H'I+II'  - GEO+WOR
melt_basin_D(12,2) =        10.d0             ; err_melt_basin_D(12,2) =            4.d0               !
melt_basin_D(13,1) =  18.d0                   ; err_melt_basin_D(13,1) =       8.d0                    ! I'I" - LBC
melt_basin_D(14,1) =   4.d0                   ; err_melt_basin_D(14,1) =       4.d0                    ! I"J  - none in Depoorter => attribute Rignot's value and remove from I'K' upscaling
melt_basin_D(15,1) =  50.d0                   ; err_melt_basin_D(15,1) =      40.d0                    ! JK   - FRIS
melt_basin_D(16,1) =  26.d0                   ; err_melt_basin_D(16,1) =      16.d0                    ! KK'  - BRL

! small correction for {Quar+Ekstrom_Atka}: to move from kbasin=1 to kbasin=16 based on surface area in QGIS :
aa =   ( 2131.d0 + 6870.d0 + 1993.d0 ) &
&    / ( 2131.d0 + 6870.d0 + 1993.d0 + 10845.d0 + 40947.d0 + 2096.d0 + 7321.d0 + 8571.d0 + 21615.d0 + 33129.d0 )
melt_basin_D( 1,1) = melt_basin_D( 1,1) - aa * ( melt_basin_D( 1,1) + melt_basin_D( 1,2) )
melt_basin_D(16,1) = melt_basin_D(16,1) + aa * ( melt_basin_D( 1,1) + melt_basin_D( 1,2) )
err_melt_basin_D( 1,1) = err_melt_basin_D( 1,1) - aa * (err_melt_basin_D( 1,1) + err_melt_basin_D( 1,2) )
err_melt_basin_D(16,1) = err_melt_basin_D(16,1) + aa * (err_melt_basin_D( 1,1) + err_melt_basin_D( 1,2) )

!-------------------------------------------------------
IF ( ll_upscaling ) THEN
  ! Weddell Sea Upscaling :
  aa = SUM( melt_basin_D(13,:) + melt_basin_D(14,:) + melt_basin_D(15,:) + melt_basin_D(16,:) )
  melt_basin_D(13,:) = melt_basin_D(13,:) * ( 1.d0 + 9.d0/aa ) ! 9.d0 = 13.d0 - melt_basin_D(14,1)  [see comment above]
  melt_basin_D(14,:) = melt_basin_D(14,:) * ( 1.d0 + 9.d0/aa )
  melt_basin_D(15,:) = melt_basin_D(15,:) * ( 1.d0 + 9.d0/aa )
  melt_basin_D(16,:) = melt_basin_D(16,:) * ( 1.d0 + 9.d0/aa )
  err_melt_basin_D(13,:) = err_melt_basin_D(13,:) * ( 1.d0 + 2.d0/aa )
  err_melt_basin_D(14,:) = err_melt_basin_D(14,:) * ( 1.d0 + 2.d0/aa )
  err_melt_basin_D(15,:) = err_melt_basin_D(15,:) * ( 1.d0 + 2.d0/aa )
  err_melt_basin_D(16,:) = err_melt_basin_D(16,:) * ( 1.d0 + 2.d0/aa )
  ! West Indian Ocean upscaling :
  aa = melt_basin_D( 1,2) + SUM( melt_basin_D( 2,:) + melt_basin_D( 3,:) ) + melt_basin_D( 4,1)
  melt_basin_D( 1,2) = melt_basin_D( 1,2) * ( 1.d0 + 40.d0/aa )
  melt_basin_D( 2,:) = melt_basin_D( 2,:) * ( 1.d0 + 40.d0/aa )
  melt_basin_D( 3,:) = melt_basin_D( 3,:) * ( 1.d0 + 40.d0/aa )
  melt_basin_D( 4,1) = melt_basin_D( 4,1) * ( 1.d0 + 40.d0/aa )
  err_melt_basin_D( 1,2) = err_melt_basin_D( 1,2) * ( 1.d0 +  6.d0/aa )
  err_melt_basin_D( 2,:) = err_melt_basin_D( 2,:) * ( 1.d0 +  6.d0/aa )
  err_melt_basin_D( 3,:) = err_melt_basin_D( 3,:) * ( 1.d0 +  6.d0/aa )
  err_melt_basin_D( 4,1) = err_melt_basin_D( 4,1) * ( 1.d0 +  6.d0/aa )
  ! East Indian Ocean upscaling :
  aa = melt_basin_D( 4,2) + SUM( melt_basin_D( 5,:) + melt_basin_D( 6,:) ) + melt_basin_D( 7,1)
  melt_basin_D( 4,2) = melt_basin_D( 4,2) * ( 1.d0 + 82.d0/aa )
  melt_basin_D( 5,:) = melt_basin_D( 5,:) * ( 1.d0 + 82.d0/aa )
  melt_basin_D( 6,:) = melt_basin_D( 6,:) * ( 1.d0 + 82.d0/aa )
  melt_basin_D( 7,1) = melt_basin_D( 7,1) * ( 1.d0 + 82.d0/aa )
  err_melt_basin_D( 4,2) = err_melt_basin_D( 4,2) * ( 1.d0 + 27.d0/aa )
  err_melt_basin_D( 5,:) = err_melt_basin_D( 5,:) * ( 1.d0 + 27.d0/aa )
  err_melt_basin_D( 6,:) = err_melt_basin_D( 6,:) * ( 1.d0 + 27.d0/aa )
  err_melt_basin_D( 7,1) = err_melt_basin_D( 7,1) * ( 1.d0 + 27.d0/aa )
  ! Ross Sea upscaling :
  aa =  melt_basin_D( 7,2) + SUM ( melt_basin_D( 8,:) )
  melt_basin_D( 7,2) = melt_basin_D( 7,2) * ( 1.d0 + 12.d0/aa )
  melt_basin_D( 8,:) = melt_basin_D( 8,:) * ( 1.d0 + 12.d0/aa )
  err_melt_basin_D( 7,2) = err_melt_basin_D( 7,2) * ( 1.d0 + 4.d0/aa )
  err_melt_basin_D( 8,:) = err_melt_basin_D( 8,:) * ( 1.d0 + 4.d0/aa )
  ! Amundsen Sea upscaling :
  aa = SUM( melt_basin_D( 9,:) + melt_basin_D(10,:) )
  melt_basin_D( 9,:) = melt_basin_D( 9,:) * ( 1.d0 + 89.d0/aa )
  melt_basin_D(10,:) = melt_basin_D(10,:) * ( 1.d0 + 89.d0/aa )
  err_melt_basin_D( 9,:) = err_melt_basin_D( 9,:) * ( 1.d0 + 20.d0/aa )
  err_melt_basin_D(10,:) = err_melt_basin_D(10,:) * ( 1.d0 + 20.d0/aa )
  ! Bellingshausen Sea upscaling :
  aa = SUM( melt_basin_D(11,:) + melt_basin_D(12,:) )
  melt_basin_D(11,:) = melt_basin_D(11,:) * ( 1.d0 + 25.d0/aa )
  melt_basin_D(12,:) = melt_basin_D(12,:) * ( 1.d0 + 25.d0/aa )
  err_melt_basin_D(11,:) = err_melt_basin_D(11,:) * ( 1.d0 + 6.d0/aa )
  err_melt_basin_D(12,:) = err_melt_basin_D(12,:) * ( 1.d0 + 6.d0/aa )
ENDIF

!--------------------------------------------------------------------------------
! Apply corrections to conserve the median melt of each basin [put 1.d0 to start]
! (needed because negative melt not allowed with the quadratic parameterization
!  so only positive extreme melts are taken when sampling the distribution ) :
!correc_med_R(:) = 1.d0
!correc_med_D(:) = 1.d0 
correc_med_R = (/ 1.0000, 1.0000, 0.9834, 1.0000, 1.0000, 1.0000, 1.0000, 0.9545, 1.0000, 1.0000, 1.0000, 1.0000, 0.5228, 0.4525, 1.0000, 0.9377 /)
correc_med_D = (/ 1.0000, 1.0000, 0.9793, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.8667, 0.9202, 0.9471 /)

!-------------------------------------------------------------
! First, finding K0 value (single value for the entire AIS) :

Nstat2 = Nstat/10
if ( nn_data_melt .eq. 1 ) then ! Rignot
  NminR=1
  NmaxR=Nstat2
  NminD=1
  NmaxD=0
elseif ( nn_data_melt .eq. 2 ) then ! Depoorter
  NminR=1
  NmaxR=0
  NminD=1
  NmaxD=Nstat2
elseif ( nn_data_melt .eq. 3 ) then ! Rignot & Depoorter
  NminR=1
  NmaxR=Nstat2/2
  NminD=Nstat2/2 + 1
  NmaxD=Nstat2
endif

ALLOCATE( thermo2d(mx,my) )
xmelt(:) = 0.d0
DO kstat=1,Nstat2 ! less sampling (only one melt rate to match)
   positive = .false.
   do while ( .not. positive ) ! we only consider positive values for K
     thermo2d(:,:) = 0.d0
     xmelt(kstat) = 0.d0
     ! Sampling observation-based melt rate (including uncertainty) :
     CALL RANDOM_NUMBER(r1) ! 1st random number
     CALL RANDOM_NUMBER(r2) ! 2nd random number
     if ( kstat .ge. NminR .and. kstat .le. NmaxR ) then ! Rignot et al. (2013)
       xmelt(kstat) = 1325.d0 + 175.d0 * sqrt ( - 2.d0 * log ( r1 ) ) * cos ( 2.d0 * pi * r2 )  ! total AIS melt in Gt/yr
     elseif ( kstat .ge. NminD .and. kstat .le. NmaxD ) then ! Depoorter et al. (2013)
       xmelt(kstat) = 1193.d0 + 163.d0 * sqrt ( - 2.d0 * log ( r1 ) ) * cos ( 2.d0 * pi * r2 )  ! total AIS melt in Gt/yr
     endif
     ! Sampling observed thermal forcing (including uncertainty if required) :
     do kbasin=1,Nbasin
       if ( ll_TS_uncertainty ) then
         CALL RANDOM_NUMBER(r1)
         CALL RANDOM_NUMBER(r2)
         ! sample normal distribution of mean=0 and stddev=err_TF :
         rr = err_TF * sqrt ( - 2.d0 * log ( r1 ) ) * cos ( 2.d0 * pi * r2 )
       else
         rr = 0.d0
       endif
       ! identify ice-shelves in this basin
       tmpmsk(:,:) = 0
       do ii=1,mx
       do jj=1,my
         if ( isfmask(ii,jj) .ge. 2 .and. basinNumber(ii,jj) .eq. kbasin-1 ) then
           tmpmsk(ii,jj) = 1
         endif
       enddo
       enddo
       meanTF  = sum(sum(tmpmsk(:,:)*TF_draft(:,:),2),1) / sum(sum(tmpmsk(:,:),2),1)
       !
       if ( para .eq. 1 ) then
         thermo2d(:,:) = thermo2d(:,:) + tmpmsk(:,:)*(TF_draft(:,:)+rr)*abs(meanTF+rr)*cste*dx*dy
       elseif ( para .eq. 2 ) then
         thermo2d(:,:) = thermo2d(:,:) + tmpmsk(:,:)*(TF_draft(:,:)+rr)*(TF_draft(:,:)+rr)*cste*dx*dy
       endif
     enddo
     thermo(kstat) = sum(sum(thermo2d(:,:),2),1)
     ! K coefficient of same sign for xmelt(kstat) and thermo(kstat) (i.e. allows melt<0 if TF<0)
     if ( sign(1.d0,xmelt(kstat)*thermo(kstat)) .gt. 0.d0 ) then
      Kfac(kstat) = xmelt(kstat)*1.d9 / thermo(kstat)
      positive = .true.
     endif
   enddo
ENDDO
DEALLOCATE( thermo2d)


!Npct = 15
Npct = 3
ALLOCATE( perct(Npct), K0_pct(Npct), DeltaT(Nstat,Npct) )
ALLOCATE( DeltaT_pc05(Nbasin,Npct) )
ALLOCATE( DeltaT_pc50(Nbasin,Npct) )
ALLOCATE( DeltaT_pc95(Nbasin,Npct) )
!perct(:) = (/0.01, 0.05, 0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75, 0.80, 0.90, 0.95, 0.99 /)
perct(:) = (/0.05, 0.50, 0.95/) !! keep median in the middle
indx2=(/ (ii, ii=1,Nstat2) /)
CALL hpsort_eps_epw(Nstat2,Kfac,indx2,1.d-12) ! this sorts Kfac
do ipct=1,Npct
  K0_pct(ipct) = Kfac(NINT(perct(ipct)*Nstat2)) 
enddo

write(*,969) K0_pct(1), K0_pct(1+Npct/2), K0_pct(Npct)
969 FORMAT('K0 : ',d12.4,' ',d12.4,' ',d12.4)

!----------
! Then, find Delta T values in each basin

true_median          =  0.d0
mean_melt_basin(:)   =  0.d0
median_melt_basin(:) =  0.d0
std_melt_basin(:)    = -1.d36
DeltaT_pc05(:,:)     = -1.d36
DeltaT_pc50(:,:)     = -1.d36
DeltaT_pc95(:,:)     = -1.d36
istwice(:) = 0

DO kbasin=1,Nbasin

  ll_to_do = .true.
  do kk=1,Nmerge
    if ( kbasin-1 .eq. NN(kk) ) ll_to_do = .false. 
  enddo

  IF ( .NOT. ll_to_do ) THEN
    write(*,*) kbasin, 'MERGED WITH PREVIOUS BASIN >>>>>>> NOTHING TO DO'
    CYCLE
  ENDIF

  ! identify ice-shelves with known melt rate in each IMBIE basin :
  isthere(:) = 0
  tmpmsk(:,:) = 0
  DO kisf=10,75
    do ii=1,mx
    do jj=1,my
      ! ice-shelves in this basin (removing ice-shelves covering more than their main basin)
      if ( isfmask(ii,jj) .eq. kisf .and. basinNumber(ii,jj) .eq. kbasin-1  &
      &    .and. .not. ( kisf .eq. 10 .and. kbasin-1 .ne.  7 )              &  ! Ross
      &    .and. .not. ( kisf .eq. 18 .and. kbasin-1 .ne.  5 )              &  ! Cook
      &    .and. .not. ( kisf .eq. 21 .and. kbasin-1 .ne. 14 )              &  ! Ronne
      &    .and. .not. ( kisf .eq. 31 .and. kbasin-1 .ne.  2 )              &  ! Amery
      &    .and. .not. ( kisf .eq. 32 .and. kbasin-1 .ne.  3 )              &  ! Conger/Glenzer 
      &    .and. .not. ( kisf .eq. 39 .and. kbasin-1 .ne. 10 )              &  ! Abbot
      &    .and. .not. ( kisf .eq. 42 .and. kbasin-1 .ne. 12 )              &  ! Larsen C
      &    .and. .not. ( kisf .eq. 45 .and. kbasin-1 .ne.  0 )              &  ! Roi Baudoin
      &    .and. .not. ( kisf .eq. 62 .and. kbasin-1 .ne.  5 )              &  ! Dibble
      &    .and. .not. ( kisf .eq. 70 .and. kbasin-1 .ne.  0 )     ) then      ! Jelbart 
        tmpmsk(ii,jj) = 1
        isthere(kisf) = 1
      endif 
    enddo
    enddo
    if ( isthere(kisf) .eq. 1 ) istwice(kisf) = istwice(kisf) + 1
  ENDDO

  areaTF  = sum(sum(tmpmsk(:,:)*dx*dy,2),1)
  ! mean thermal forcing with zero DeltaT :
  meanTF  = sum(sum(tmpmsk(:,:)*TF_draft(:,:)*dx*dy,2),1) / areaTF
  meanTF2 = sum(sum(tmpmsk(:,:)*TF_draft(:,:)*TF_draft(:,:)*dx*dy,2),1) / areaTF

  if ( nn_data_melt .eq. 1 ) then ! Rignot
    NminR=1
    NmaxR=Nstat
    NminD=1
    NmaxD=0
    do kisf=10,75
      if ( isthere(kisf) .eq. 1 ) then
        true_median(kbasin) = true_median(kbasin) + melt_isf(kisf)
      endif
    enddo
  elseif ( nn_data_melt .eq. 2 ) then ! Depoorter
    NminR=1
    NmaxR=0
    NminD=1
    NmaxD=Nstat
    do k4=1,N4
      true_median(kbasin) = true_median(kbasin) + melt_basin_D(kbasin,k4)
    enddo
  elseif ( nn_data_melt .eq. 3 ) then ! Rignot & Depoorter
    NminR=1
    NmaxR=Nstat/2
    NminD=Nstat/2 + 1
    NmaxD=Nstat
    do kisf=10,75
      if ( isthere(kisf) .eq. 1 ) then
        true_median(kbasin) = true_median(kbasin) + 0.5 * melt_isf(kisf)
      endif
    enddo
    do k4=1,N4
      true_median(kbasin) = true_median(kbasin) + 0.5 * melt_basin_D(kbasin,k4)
    enddo
  endif

  !--- Calculate the K coefficients :
  xmelt(:) = 0.d0
  DO kstat=1,Nstat
   positive = .false.
   do while ( .not. positive ) ! we only consider positive values for K
     xmelt(kstat) = 0.d0 ! reset in case .not. positive after a loop
     ! Sampling observation-based melt rate (including uncertainty) :
     if ( kstat .ge. NminR .and. kstat .le. NmaxR ) then ! Rignot et al. (2013)
       do kisf=10,75 ! all monitored ice shelves
         if ( isthere(kisf) .eq. 1 ) then
           CALL RANDOM_NUMBER(r1) ! 1st random number
           CALL RANDOM_NUMBER(r2) ! 2nd random number
           ! sample normal distribution of mean=melt_isf and stddev=err_melt_isf :
           xmelt (kstat) = xmelt (kstat) + correc_med_R(kbasin) * ( melt_isf(kisf) + err_melt_isf(kisf) * sqrt ( - 2.d0 * log ( r1 ) ) * cos ( 2.d0 * pi * r2 ) )
         endif
       enddo
     elseif ( kstat .ge. NminD .and. kstat .le. NmaxD ) then ! Depoorter et al. (2013)
       do k4=1,N4
         CALL RANDOM_NUMBER(r1)
         CALL RANDOM_NUMBER(r2)
         xmelt(kstat) = xmelt(kstat) + correc_med_D(kbasin) * ( melt_basin_D(kbasin,k4) + err_melt_basin_D(kbasin,k4) * sqrt ( - 2.d0 * log ( r1 ) ) * cos ( 2.d0 * pi * r2 ) ) ! in Gt/yr
       enddo
     endif
     ! Sampling observed thermal forcing (including uncertainty if required) :
     if ( ll_TS_uncertainty ) then
       CALL RANDOM_NUMBER(r1)
       CALL RANDOM_NUMBER(r2)
       ! sample normal distribution of mean=0 and stddev=err_TF :
       rr = err_TF * sqrt ( - 2.d0 * log ( r1 ) ) * cos ( 2.d0 * pi * r2 )
     else
       rr = 0.d0
     endif
     if ( xmelt(kstat) .gt. 0.d0 ) then
       if ( para .eq. 1 ) then
         do ipct=1,Npct
           CALL RANDOM_NUMBER(r1)
           ! resolving 2nd order equation by replacing TF by TF + DeltaT in the melt parameterization :
           dis = xmelt(kstat) / ( K0_pct(ipct) * areaTF * cste * 1.d-9 ) ! discriminent: if K0_pct matches perfectly for this basin, dis = meanTF^2
           if ( dis .ge. 0 ) then
             DeltaT(kstat,ipct) = - meanTF + sqrt(dis)                     ! DeltaT needs to be zero if K0_pct matches perfectly for this basin
           else
             DeltaT(kstat,ipct) = 99999.d0 * SIGN(1.d0,r1-0.5d0) ! randomly either -99999 or + 99999
           endif
         enddo
       elseif ( para .eq. 2 ) then
         do ipct=1,Npct
           CALL RANDOM_NUMBER(r1)
           ! resolving 2nd order equation by replacing TF by TF + DeltaT in the melt parameterization :
           dis = meanTF*meanTF - meanTF2 + xmelt(kstat) / ( K0_pct(ipct) * areaTF * cste * 1.d-9 ) ! discriminent.
           if ( dis .ge. 0 ) then
             DeltaT(kstat,ipct) = - meanTF + sqrt(dis)                    ! DeltaT needs to be zero if K0_pct matches perfectly for this basin
           else
             DeltaT(kstat,ipct) = 99999.d0 * SIGN(1.d0,r1-0.5d0) ! randomly either -99999 or + 99999
           endif
         enddo
       endif
       positive = .true.
     endif
   enddo
  ENDDO

  mean_melt_basin(kbasin) = sum(xmelt) / Nstat   ! in Gt/yr
  std_melt_basin(kbasin)  = sqrt( sum( ( xmelt(:) - mean_melt_basin(kbasin) )**2 ) / Nstat )
  indx=(/ (ii, ii=1,Nstat) /)
  CALL hpsort_eps_epw(Nstat,xmelt,indx,1.d-12) ! this sorts xmelt
  median_melt_basin(kbasin) = xmelt(NINT(0.50*Nstat))
  CALL hpsort_eps_epw(Nstat,Kfac,indx,1.d-12) ! this sorts Kfac
  do ipct=1,Npct
    DeltaT_pc05(kbasin,ipct) = DeltaT(NINT(0.05*Nstat),ipct) !  5th pct of DeltaT matching with either 5th, 50th, or 95th pct of K0
    DeltaT_pc50(kbasin,ipct) = DeltaT(NINT(0.50*Nstat),ipct) ! 50th pct of DeltaT matching with either 5th, 50th, or 95th pct of K0
    DeltaT_pc95(kbasin,ipct) = DeltaT(NINT(0.95*Nstat),ipct) ! 95th pct of DeltaT matching with either 5th, 50th, or 95th pct of K0
  enddo

  write(*,962) kbasin, meanTF, DeltaT_pc05(kbasin,Npct), DeltaT_pc50(kbasin,1+Npct/2), DeltaT_pc95(kbasin,1)

ENDDO ! kbasin=1,Nbasin

962 FORMAT(i3,' ',f7.3,' ',f7.3,' ',f7.3,' ',f7.3)

write(*,*) ' '
write(*,*) 'WARNING : IF FOLLOWING RATIO NOT CLOSE TO 1.0 YOU NEED TO CORRECT THE PRESCRIBED MEAN'
write(*,*) '          BY TUNING correc_med_R OR correc_med_D'
DO kbasin=1,Nbasin
  write(*,*) '          ratio of true mean to statistical mean = ', true_median(kbasin) / median_melt_basin(kbasin)
ENDDO
write(*,*) ' '
write(*,*) true_median(:) / median_melt_basin(:)
write(*,*) ' '

! check if some basins are counted twice :
do kisf=1,Nisf
  if ( istwice(kisf) .gt. 1 ) then
    write(*,*) ' '
    write(*,*) ' !@#$%^ WARNING : ice shelf ', kisf, ' is counted ', istwice(kisf), ' times'
    write(*,*) ' '
  endif
enddo

! assign correct values for the basins that have been merged with kbasin-1
do kk=1,Nmerge
  mean_melt_basin   (NN(kk)+1) = 0.d0
  std_melt_basin    (NN(kk)+1) = 0.d0
  median_melt_basin (NN(kk)+1) = 0.d0
  DeltaT_pc05 (NN(kk)+1,:) = DeltaT_pc05 (NN(kk),:)
  DeltaT_pc50 (NN(kk)+1,:) = DeltaT_pc50 (NN(kk),:)
  DeltaT_pc95 (NN(kk)+1,:) = DeltaT_pc95 (NN(kk),:)
enddo

! Check sum :
write(*,*) ' '
write(*,*) 'STATISTICAL SUM = ', SUM(median_melt_basin), ' Gt/yr'
write(*,*) 'vs EXPECTED SUM = ', SUM(true_median), ' Gt/yr'

!----------------------------------------------------------
! 2d DeltaT :

ALLOCATE( melt_DeltaT_pc05(mx,my,Npct), melt_DeltaT_pc50(mx,my,Npct), melt_DeltaT_pc95(mx,my,Npct) )
ALLOCATE( DeltaT_basin_pc05(mx,my,Npct), DeltaT_basin_pc50(mx,my,Npct), DeltaT_basin_pc95(mx,my,Npct) )

DO kbasin=1,Nbasin
  do ii=1,mx
  do jj=1,my
    if ( basinNumber(ii,jj) .eq. kbasin-1 ) then
      DeltaT_basin_pc05(ii,jj,:) = DeltaT_pc05(kbasin,:)
      DeltaT_basin_pc50(ii,jj,:) = DeltaT_pc50(kbasin,:)
      DeltaT_basin_pc95(ii,jj,:) = DeltaT_pc95(kbasin,:)
    endif
  enddo
  enddo
ENDDO

!-------------------------------
! Smoothing DeltaT transitions :

if ( Nsmooth .ge. 2 ) then
  CALL smooth(mx,my,Npct,Nsmooth,DeltaT_basin_pc05)
  CALL smooth(mx,my,Npct,Nsmooth,DeltaT_basin_pc50)
  CALL smooth(mx,my,Npct,Nsmooth,DeltaT_basin_pc95)
endif

!----------------------------
! Recalculating melt rates :

check_max_DeltaT_pc05 = 0.d0 
check_max_DeltaT_pc50 = 0.d0 
check_max_DeltaT_pc95 = 0.d0 

write(*,*) ' '
DO kbasin=1,Nbasin
  ! identify ice-shelves
  tmpmsk(:,:) = 0
  check_median = 0.d0
  do ii=1,mx
  do jj=1,my
    if ( isfmask(ii,jj) .ge. 2 .and. kbasin-1 .eq. basinNumber(ii,jj) ) then
      tmpmsk(ii,jj) = 1
    endif 
  enddo
  enddo
  if ( para .eq. 1 ) then
    meanTF = sum(sum(tmpmsk(:,:)*TF_draft(:,:),2),1) / sum(sum(tmpmsk(:,:),2),1)
    do ii=1,mx
    do jj=1,my
      if ( tmpmsk(ii,jj) .eq. 1 ) then
        do ipct=1,Npct
          melt_DeltaT_pc05(ii,jj,ipct) = K0_pct(ipct) * cste * (meanTF+DeltaT_pc05(kbasin,ipct)) * (TF_draft(ii,jj)+DeltaT_pc05(kbasin,ipct)) ! in meters of pure water per year
          melt_DeltaT_pc50(ii,jj,ipct) = K0_pct(ipct) * cste * (meanTF+DeltaT_pc50(kbasin,ipct)) * (TF_draft(ii,jj)+DeltaT_pc50(kbasin,ipct)) ! in meters of pure water per year
          melt_DeltaT_pc95(ii,jj,ipct) = K0_pct(ipct) * cste * (meanTF+DeltaT_pc95(kbasin,ipct)) * (TF_draft(ii,jj)+DeltaT_pc95(kbasin,ipct)) ! in meters of pure water per year
          check_max_DeltaT_pc05 = max(check_max_DeltaT_pc05,melt_DeltaT_pc05(ii,jj,ipct))
          check_max_DeltaT_pc50 = max(check_max_DeltaT_pc50,melt_DeltaT_pc50(ii,jj,ipct))
          check_max_DeltaT_pc95 = max(check_max_DeltaT_pc95,melt_DeltaT_pc95(ii,jj,ipct))
        enddo
        check_median = check_median + K0_pct(1+Npct/2) * cste * (meanTF+DeltaT_pc50(kbasin,1+Npct/2)) * (TF_draft(ii,jj)+DeltaT_pc50(kbasin,1+Npct/2)) * dx * dy * 1.d-9  ! Gt/yr
      endif
    enddo
    enddo
  elseif ( para .eq. 2 ) then
    do ii=1,mx
    do jj=1,my
      if ( tmpmsk(ii,jj) .eq. 1 ) then
        do ipct=1,Npct
          melt_DeltaT_pc05(ii,jj,ipct) = K0_pct(ipct) * cste * (TF_draft(ii,jj)+DeltaT_pc05(kbasin,ipct)) * (TF_draft(ii,jj)+DeltaT_pc05(kbasin,ipct)) ! in meters of pure water per year
          melt_DeltaT_pc50(ii,jj,ipct) = K0_pct(ipct) * cste * (TF_draft(ii,jj)+DeltaT_pc50(kbasin,ipct)) * (TF_draft(ii,jj)+DeltaT_pc50(kbasin,ipct)) ! in meters of pure water per year
          melt_DeltaT_pc95(ii,jj,ipct) = K0_pct(ipct) * cste * (TF_draft(ii,jj)+DeltaT_pc95(kbasin,ipct)) * (TF_draft(ii,jj)+DeltaT_pc95(kbasin,ipct)) ! in meters of pure water per year
          check_max_DeltaT_pc05 = max(check_max_DeltaT_pc05,melt_DeltaT_pc05(ii,jj,ipct))
          check_max_DeltaT_pc50 = max(check_max_DeltaT_pc50,melt_DeltaT_pc50(ii,jj,ipct))
          check_max_DeltaT_pc95 = max(check_max_DeltaT_pc95,melt_DeltaT_pc95(ii,jj,ipct))
        enddo
        check_median = check_median + K0_pct(1+Npct/2) * cste * (TF_draft(ii,jj)+DeltaT_pc05(kbasin,1+Npct/2)) * (TF_draft(ii,jj)+DeltaT_pc05(kbasin,1+Npct/2)) * dx * dy * 1.d-9  ! Gt/yr
      endif
    enddo
    enddo
  endif
  write(*,848) kbasin, check_median, true_median(kbasin), median_melt_basin(kbasin)
ENDDO
848 FORMAT('Check median melt for basin ',i3,' : ',f6.2,'   vs true : ',f6.2,'   vs sampled : ',f6.2,'  Gt/yr')

write(*,*) ' '
write(*,*) 'CHECKING maximum melt rates :'
write(*,*) '   check_max_DeltaT_pc05 = ', check_max_DeltaT_pc05
write(*,*) '   check_max_DeltaT_pc50 = ', check_max_DeltaT_pc50
write(*,*) '   check_max_DeltaT_pc95 = ', check_max_DeltaT_pc95

! Checking sum when applied to all ice shelves :
write(*,*) ' '
write(*,*) 'vs REINTEGRATRED SUM = ', sum( sum( melt_DeltaT_pc50(:,:,1+Npct/2), 2 ), 1) * dx * dy * 1.d-9

! Checking ice shelf per ice shelf:
DO kisf=2,Nisf
  ! identify ice-shelves
  tmpmsk(:,:) = 0
  do ii=1,mx
  do jj=1,my
    if ( isfmask(ii,jj) .eq. kisf ) then
      tmpmsk(ii,jj) = 1
    endif 
  enddo
  enddo
  check_median = 0.d0 
  if ( para .eq. 1 ) then
    meanTF = sum(sum(tmpmsk(:,:)*(TF_draft(:,:)+DeltaT_basin_pc50(:,:,1+Npct/2)),2),1) / sum(sum(tmpmsk(:,:),2),1)
    do ii=1,mx
    do jj=1,my
      if ( tmpmsk(ii,jj) .eq. 1 ) then
        check_median = check_median + K0_pct(1+Npct/2) * cste * meanTF * (TF_draft(ii,jj)+DeltaT_basin_pc50(ii,jj,1+Npct/2)) * dx * dy * 1.d-9   ! in Gt/yr
      endif
    enddo
    enddo
  elseif ( para .eq. 2 ) then
    do ii=1,mx
    do jj=1,my
      if ( tmpmsk(ii,jj) .eq. 1 ) then
        check_median = check_median + K0_pct(1+Npct/2) * cste * (TF_draft(ii,jj)+DeltaT_basin_pc50(ii,jj,1+Npct/2)) * (TF_draft(ii,jj)+DeltaT_basin_pc50(ii,jj,1+Npct/2)) * dx * dy * 1.d-9   ! in Gt/yr
      endif
    enddo
    enddo
  endif
  write(*,849) kisf, check_median, melt_isf(kisf)
ENDDO
849 FORMAT('TOTAL MELT FOR ICE SHELF ',i3,' : ',f6.2,'   vs PRESCRIBED : ',f6.2,'   Gt/yr')

!------------------------------------------------
! Writing new netcdf file containing melt rates:

write(*,*) ' ' 
write(*,*) 'Creating ', TRIM(file_melt_out)
 
status = NF90_CREATE(TRIM(file_melt_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidM,"y",my,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"x",mx,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"pct",Npct,dimID_pct); call erreur(status,.TRUE.,"def_dimID_pct")
  
status = NF90_DEF_VAR(fidM,"y",NF90_DOUBLE,(/dimID_y/),y_ID); call erreur(status,.TRUE.,"def_var_y_ID")
status = NF90_DEF_VAR(fidM,"x",NF90_DOUBLE,(/dimID_x/),x_ID); call erreur(status,.TRUE.,"def_var_x_ID")
status = NF90_DEF_VAR(fidM,"pct",NF90_DOUBLE,(/dimID_pct/),pct_ID); call erreur(status,.TRUE.,"def_var_pct_ID")
status = NF90_DEF_VAR(fidM,"K0",NF90_DOUBLE,(/dimID_pct/),K0_ID); call erreur(status,.TRUE.,"def_var_K0_ID")
status = NF90_DEF_VAR(fidM,"melt_DeltaT_pc05",NF90_DOUBLE,(/dimID_x,dimID_y,dimID_pct/),melt_DeltaT_pc05_ID); call erreur(status,.TRUE.,"def_var_melt_DeltaT_pc05_ID")
status = NF90_DEF_VAR(fidM,"melt_DeltaT_pc50",NF90_DOUBLE,(/dimID_x,dimID_y,dimID_pct/),melt_DeltaT_pc50_ID); call erreur(status,.TRUE.,"def_var_melt_DeltaT_pc50_ID")
status = NF90_DEF_VAR(fidM,"melt_DeltaT_pc95",NF90_DOUBLE,(/dimID_x,dimID_y,dimID_pct/),melt_DeltaT_pc95_ID); call erreur(status,.TRUE.,"def_var_melt_DeltaT_pc95_ID")

status = NF90_PUT_ATT(fidM,pct_ID,"long_name","K0 percentile") ; call erreur(status,.TRUE.,"put_att_pct_ID")
status = NF90_PUT_ATT(fidM,pct_ID,"units","-") ; call erreur(status,.TRUE.,"put_att_pct_ID")
status = NF90_PUT_ATT(fidM,K0_ID,"long_name","K0 coefficient") ; call erreur(status,.TRUE.,"put_att_K0_ID")
status = NF90_PUT_ATT(fidM,K0_ID,"units","-") ; call erreur(status,.TRUE.,"put_att_K0_ID")
status = NF90_PUT_ATT(fidM,melt_DeltaT_pc05_ID,"long_name","melt rate obtained with basin-5th-percentile DeltaT values"); call erreur(status,.TRUE.,"put_att_melt_DeltaT_pc05_ID")
status = NF90_PUT_ATT(fidM,melt_DeltaT_pc05_ID,"units","m.w.e / yr"); call erreur(status,.TRUE.,"put_att_melt_DeltaT_pc05_ID")
status = NF90_PUT_ATT(fidM,melt_DeltaT_pc50_ID,"long_name","melt rate obtained with basin-median DeltaT values"); call erreur(status,.TRUE.,"put_att_melt_DeltaT_pc50_ID")
status = NF90_PUT_ATT(fidM,melt_DeltaT_pc50_ID,"units","m.w.e / yr"); call erreur(status,.TRUE.,"put_att_melt_DeltaT_pc50_ID")
status = NF90_PUT_ATT(fidM,melt_DeltaT_pc95_ID,"long_name","melt rate obtained with basin-95th-percentile DeltaT values"); call erreur(status,.TRUE.,"put_att_melt_DeltaT_pc95_ID")
status = NF90_PUT_ATT(fidM,melt_DeltaT_pc95_ID,"units","m.w.e / yr"); call erreur(status,.TRUE.,"put_att_melt_DeltaT_pc95_ID")
 
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using calculate_K_coefficients_quadratic.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
if     ( para .eq. 1 ) then
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"melt parameterization","Quadratic non-local"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( para .eq. 2 ) then
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"melt parameterization","Quadratic local"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
endif
if     ( nn_data_melt .eq. 1 .and. ll_upscaling ) then
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"calibration melt rates","from Rignot et al. (2013) with upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( nn_data_melt .eq. 2 .and. ll_upscaling ) then
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"calibration melt rates","from Depoorter et al. (2013) with upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( nn_data_melt .eq. 3 .and. ll_upscaling ) then
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"calibration melt rates","from both Rignot et al. (2013) and Depoorter et al. (2013) with upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( nn_data_melt .eq. 1 ) then
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"calibration melt rates","from Rignot et al. (2013) without upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( nn_data_melt .eq. 2 ) then
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"calibration melt rates","from Depoorter et al. (2013) without upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( nn_data_melt .eq. 3 ) then
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"calibration melt rates","from both Rignot et al. (2013) and Depoorter et al. (2013) without upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
endif
if ( ll_TS_uncertainty ) then
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"calibration thermal forcing","expanded WOA2013-v2+MEOP and its uncertainty"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
else
  status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"calibration thermal forcing","expanded WOA2013-v2+MEOP (no uncertainty)"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
endif
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"calibration topography","BEDMAP2"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"basins","IMBIE2"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"number of samples",Nstat); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"smoothing window",Nsmooth); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")

status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 
 
status = NF90_PUT_VAR(fidM,y_ID,y); call erreur(status,.TRUE.,"var_y_ID")
status = NF90_PUT_VAR(fidM,x_ID,x); call erreur(status,.TRUE.,"var_x_ID")
status = NF90_PUT_VAR(fidM,pct_ID,perct); call erreur(status,.TRUE.,"var_pct_ID")
status = NF90_PUT_VAR(fidM,K0_ID,K0_pct); call erreur(status,.TRUE.,"var_K0_ID")
status = NF90_PUT_VAR(fidM,melt_DeltaT_pc05_ID,melt_DeltaT_pc05); call erreur(status,.TRUE.,"var_melt_DeltaT_pc05_ID")
status = NF90_PUT_VAR(fidM,melt_DeltaT_pc50_ID,melt_DeltaT_pc50); call erreur(status,.TRUE.,"var_melt_DeltaT_pc50_ID")
status = NF90_PUT_VAR(fidM,melt_DeltaT_pc95_ID,melt_DeltaT_pc95); call erreur(status,.TRUE.,"var_melt_DeltaT_pc95_ID")
 
status = NF90_CLOSE(fidM); call erreur(status,.TRUE.,"final")

!----------------------------------------------------
! Writing new netcdf file containing K coefficients:
 
write(*,*) 'Creating ', TRIM(file_K_out)
 
status = NF90_CREATE(TRIM(file_K_out),NF90_NOCLOBBER,fidK); call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidK,"y",my,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidK,"x",mx,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidK,"pct",Npct,dimID_pct); call erreur(status,.TRUE.,"def_dimID_pct")
  
status = NF90_DEF_VAR(fidK,"y",NF90_DOUBLE,(/dimID_y/),y_ID); call erreur(status,.TRUE.,"def_var_y_ID")
status = NF90_DEF_VAR(fidK,"x",NF90_DOUBLE,(/dimID_x/),x_ID); call erreur(status,.TRUE.,"def_var_x_ID")
status = NF90_DEF_VAR(fidK,"pct",NF90_DOUBLE,(/dimID_pct/),pct_ID); call erreur(status,.TRUE.,"def_var_pct_ID")
status = NF90_DEF_VAR(fidK,"K0",NF90_DOUBLE,(/dimID_pct/),K0_ID); call erreur(status,.TRUE.,"def_var_K0_ID")
status = NF90_DEF_VAR(fidK,"DeltaT_basin_pc05",NF90_DOUBLE,(/dimID_x,dimID_y,dimID_pct/),DeltaT_basin_pc05_ID); call erreur(status,.TRUE.,"def_var_DeltaT_basin_pc05_ID")
status = NF90_DEF_VAR(fidK,"DeltaT_basin_pc50",NF90_DOUBLE,(/dimID_x,dimID_y,dimID_pct/),DeltaT_basin_pc50_ID); call erreur(status,.TRUE.,"def_var_DeltaT_basin_pc50_ID")
status = NF90_DEF_VAR(fidK,"DeltaT_basin_pc95",NF90_DOUBLE,(/dimID_x,dimID_y,dimID_pct/),DeltaT_basin_pc95_ID); call erreur(status,.TRUE.,"def_var_DeltaT_basin_pc95_ID")

status = NF90_PUT_ATT(fidK,pct_ID,"long_name","K0 percentile") ; call erreur(status,.TRUE.,"put_att_pct_ID")
status = NF90_PUT_ATT(fidK,pct_ID,"units","-") ; call erreur(status,.TRUE.,"put_att_pct_ID")
status = NF90_PUT_ATT(fidK,K0_ID,"long_name","K0 coefficient") ; call erreur(status,.TRUE.,"put_att_K0_ID")
status = NF90_PUT_ATT(fidK,K0_ID,"units","-") ; call erreur(status,.TRUE.,"put_att_K0_ID")
status = NF90_PUT_ATT(fidK,DeltaT_basin_pc05_ID,"long_name","basin-5th-percentile DeltaT values"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc05_ID")
status = NF90_PUT_ATT(fidK,DeltaT_basin_pc05_ID,"units","m.w.e / yr"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc05_ID")
status = NF90_PUT_ATT(fidK,DeltaT_basin_pc50_ID,"long_name","basin-median DeltaT values"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc50_ID")
status = NF90_PUT_ATT(fidK,DeltaT_basin_pc50_ID,"units","m.w.e / yr"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc50_ID")
status = NF90_PUT_ATT(fidK,DeltaT_basin_pc95_ID,"long_name","basin-95th-percentile DeltaT values"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc95_ID")
status = NF90_PUT_ATT(fidK,DeltaT_basin_pc95_ID,"units","m.w.e / yr"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc95_ID")
 
status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"history","Created using calculate_K_coefficients_quadratic.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
if     ( para .eq. 1 ) then
  status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"melt parameterization","Quadratic non-local"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( para .eq. 2 ) then
  status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"melt parameterization","Quadratic local"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
endif
if     ( nn_data_melt .eq. 1 .and. ll_upscaling ) then
  status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"calibration melt rates","from Rignot et al. (2013) with upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( nn_data_melt .eq. 2 .and. ll_upscaling ) then
  status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"calibration melt rates","from Depoorter et al. (2013) with upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( nn_data_melt .eq. 3 .and. ll_upscaling ) then
  status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"calibration melt rates","from both Rignot et al. (2013) and Depoorter et al. (2013) with upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( nn_data_melt .eq. 1 ) then
  status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"calibration melt rates","from Rignot et al. (2013) without upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( nn_data_melt .eq. 2 ) then
  status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"calibration melt rates","from Depoorter et al. (2013) without upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
elseif ( nn_data_melt .eq. 3 ) then
  status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"calibration melt rates","from both Rignot et al. (2013) and Depoorter et al. (2013) without upscaling"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
endif
if ( ll_TS_uncertainty ) then
  status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"calibration thermal forcing","expanded WOA2013-v2+MEOP and its uncertainty"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
else
  status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"calibration thermal forcing","expanded WOA2013-v2+MEOP (no uncertainty)"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
endif
status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"calibration topography","BEDMAP2"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"basins","IMBIE2"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"number of samples",Nstat); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
status = NF90_PUT_ATT(fidK,NF90_GLOBAL,"smoothing window",Nsmooth); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
 
status = NF90_ENDDEF(fidK); call erreur(status,.TRUE.,"fin_definition") 
 
status = NF90_PUT_VAR(fidK,y_ID,y); call erreur(status,.TRUE.,"var_y_ID")
status = NF90_PUT_VAR(fidK,x_ID,x); call erreur(status,.TRUE.,"var_x_ID")
status = NF90_PUT_VAR(fidK,pct_ID,perct); call erreur(status,.TRUE.,"var_pct_ID")
status = NF90_PUT_VAR(fidK,K0_ID,K0_pct); call erreur(status,.TRUE.,"var_K0_ID")
status = NF90_PUT_VAR(fidK,DeltaT_basin_pc05_ID,DeltaT_basin_pc05); call erreur(status,.TRUE.,"var_DeltaT_basin_pc05_ID")
status = NF90_PUT_VAR(fidK,DeltaT_basin_pc50_ID,DeltaT_basin_pc50); call erreur(status,.TRUE.,"var_DeltaT_basin_pc50_ID")
status = NF90_PUT_VAR(fidK,DeltaT_basin_pc95_ID,DeltaT_basin_pc95); call erreur(status,.TRUE.,"var_DeltaT_basin_pc95_ID")
 
status = NF90_CLOSE(fidK); call erreur(status,.TRUE.,"final")

end program modif

!==========
SUBROUTINE smooth(nx,ny,nz,Nsm,M)
! Smooting 2d matrix M (dimensions mx,my) with running average of width Nsm

IMPLICIT NONE
INTEGER, INTENT(IN) :: nx, ny, nz, Nsm
REAL(KIND=8), DIMENSION(nx,ny,nz), INTENT(INOUT) :: M
REAL(KIND=8), DIMENSION(nx,ny,nz) :: Ms
INTEGER :: i, j, ims, ips, jms, jps

do i=1,nx
do j=1,ny
  ims = MAX( 1,i-Nsm/2)
  ips = MIN(nx,i+Nsm/2)
  jms = MAX( 1,j-Nsm/2)
  jps = MIN(ny,j+Nsm/2)  
  Ms(i,j,:) = sum( sum( M(ims:ips,jms:jps,:), 2 ), 1 ) / ( (ips-ims+1)*(jps-jms+1) ) 
enddo
enddo

M(:,:,:) = Ms(:,:,:)

END SUBROUTINE smooth


!==========
subroutine hpsort_eps_epw (n, ra, ind, eps)
!========================================================================================
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
!                                                                            
! This file is distributed under the terms of the GNU General Public         
! License. See the file `LICENSE' in the root directory of the               
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
!                                                                            
! Adapted from flib/hpsort_eps
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! sort an array ra(1:n) into ascending order using heapsort algorithm,
! and considering two elements being equal if their values differ
! for less than "eps".
! n is input, ra is replaced on output by its sorted rearrangement.
! create an index table (ind) by making an exchange in the index array
! whenever an exchange is made on the sorted data array (ra).
! in case of equal values in the data array (ra) the values in the
! index array (ind) are used to order the entries.
! if on input ind(1)  = 0 then indices are initialized in the routine,
! if on input ind(1) != 0 then indices are assumed to have been
!                initialized before entering the routine and these
!                indices are carried around during the sorting process
!
! no work space needed !
! free us from machine-dependent sorting-routines !
!
! adapted from Numerical Recipes pg. 329 (new edition)
!
implicit none  
!-input/output variables
integer, intent(in)   :: n  
real(kind=8), intent(in)  :: eps
integer :: ind (n)  
real(kind=8) :: ra (n)
!-local variables
integer :: i, ir, j, l, iind  
real(kind=8) :: rra  
!
! initialize index array
IF (ind (1) .eq.0) then  
   DO i = 1, n  
      ind (i) = i  
   ENDDO
ENDIF
! nothing to order
IF (n.lt.2) return  
! initialize indices for hiring and retirement-promotion phase
l = n / 2 + 1  

ir = n  

sorting: do 

  ! still in hiring phase
  IF ( l .gt. 1 ) then  
     l    = l - 1  
     rra  = ra (l)  
     iind = ind (l)  
     ! in retirement-promotion phase.
  ELSE  
     ! clear a space at the end of the array
     rra  = ra (ir)  
     !
     iind = ind (ir)  
     ! retire the top of the heap into it
     ra (ir) = ra (1)  
     !
     ind (ir) = ind (1)  
     ! decrease the size of the corporation
     ir = ir - 1  
     ! done with the last promotion
     IF ( ir .eq. 1 ) then  
        ! the least competent worker at all !
        ra (1)  = rra  
        !
        ind (1) = iind  
        exit sorting  
     ENDIF
  ENDIF
  ! wheter in hiring or promotion phase, we
  i = l  
  ! set up to place rra in its proper level
  j = l + l  
  !
  DO while ( j .le. ir )  
     IF ( j .lt. ir ) then  
        ! compare to better underling
        IF ( hslt( ra (j),  ra (j + 1) ) ) then  
           j = j + 1  
        !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
           ! this means ra(j) == ra(j+1) within tolerance
         !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
        ENDIF
     ENDIF
     ! demote rra
     IF ( hslt( rra, ra (j) ) ) then  
        ra (i) = ra (j)  
        ind (i) = ind (j)  
        i = j  
        j = j + j  
     !else if ( .not. hslt ( ra(j) , rra ) ) then
        !this means rra == ra(j) within tolerance
        ! demote rra
       ! if (iind.lt.ind (j) ) then
       !    ra (i) = ra (j)
       !    ind (i) = ind (j)
       !    i = j
       !    j = j + j
       ! else
           ! set j to terminate do-while loop
       !    j = ir + 1
       ! endif
        ! this is the right place for rra
     ELSE
        ! set j to terminate do-while loop
        j = ir + 1  
     ENDIF
  ENDDO
  ra (i) = rra  
  ind (i) = iind  

END DO sorting    
contains 

!  internal function 
!  compare two real number and return the result

logical function hslt( a, b )
  REAL(kind=8) :: a, b
  IF( abs(a-b) <  eps ) then
    hslt = .false.
  ELSE
    hslt = ( a < b )
  end if
end function hslt

  !
end subroutine hpsort_eps_epw

!======================================================================
SUBROUTINE linear_regress(nn,x,y,a,b,r)
! get a, b in y = a.x + b
! and optionally correlation r

IMPLICIT NONE

INTEGER, INTENT(IN) :: nn

REAL(KIND=8), DIMENSION(nn), INTENT(IN) :: x, y

REAL(KIND=8), INTENT(OUT) :: a, b

REAL(KIND=8), INTENT(OUT), OPTIONAL :: r

REAL(KIND=8) :: sumx, sumy, sumx2, sumy2, sumxy

sumx  = SUM(x)
sumy  = SUM(y)                                                           
sumx2 = SUM(x*x)                                                             
sumy2 = SUM(y*y)                                                             
sumxy = SUM(x*y)                                                             

a = (nn * sumxy  -  sumx * sumy) / (nn * sumx2 - sumx**2)          ! slope
b = (sumy * sumx2  -  sumx * sumxy) / (nn * sumx2  -  sumx**2)     ! y-intercept
r = (sumxy - sumx * sumy / nn)                       &             ! correlation coefficient
&   / sqrt((sumx2 - sumx**2/nn) * (sumy2 - sumy**2/nn))

END SUBROUTINE linear_regress

!======================================================================

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
