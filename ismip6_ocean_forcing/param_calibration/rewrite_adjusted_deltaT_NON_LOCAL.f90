program modif
 
USE netcdf
 
IMPLICIT NONE
 
INTEGER :: fidA, status, dimID_pct, dimID_x, dimID_y, mpct, mx, my, DeltaT_basin_pc95_ID, DeltaT_basin_pc50_ID, DeltaT_basin_pc05_ID, K0_ID, pct_ID, x_ID, y_ID, fidM, ii, jj, kbasin, Nbasin
 
CHARACTER(LEN=100) :: file_in, file_out
 
REAL*8,ALLOCATABLE,DIMENSION(:) :: K0, pct, x, y, dT_median, dT_5thpct, dT_95thpct
 
REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: DeltaT_basin_pc95, DeltaT_basin_pc50, DeltaT_basin_pc05
 
INTEGER :: fidB, basinNumber_ID
 
CHARACTER(LEN=100) :: file_bn
 
INTEGER*4, ALLOCATABLE,DIMENSION(:,:) :: basinNumber
 
file_in  = 'coeff_K0_DeltaT_quadratic_non_local_Rignot_Depoorter_with_TS_uncert.nc'
file_bn  = 'basinNumbers_8km.nc'
file_out = 'coeff_K0_DeltaT_quadratic_non_local_Rignot_Depoorter_with_TS_uncert_ADJUSTED.nc'
 
!---------------------------------------
! Read netcdf input file :
 
write(*,*) 'Reading ', TRIM(file_in)
 
status = NF90_OPEN(TRIM(file_in),0,fidA); call erreur(status,.TRUE.,"read")
 
status = NF90_INQ_DIMID(fidA,"pct",dimID_pct); call erreur(status,.TRUE.,"inq_dimID_pct")
status = NF90_INQ_DIMID(fidA,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
status = NF90_INQ_DIMID(fidA,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
 
status = NF90_INQUIRE_DIMENSION(fidA,dimID_pct,len=mpct); call erreur(status,.TRUE.,"inq_dim_pct")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mx); call erreur(status,.TRUE.,"inq_dim_x")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=my); call erreur(status,.TRUE.,"inq_dim_y")
  
ALLOCATE(  DeltaT_basin_pc95(mx,my,mpct)  ) 
ALLOCATE(  DeltaT_basin_pc50(mx,my,mpct)  ) 
ALLOCATE(  DeltaT_basin_pc05(mx,my,mpct)  ) 
ALLOCATE(  K0(mpct)  ) 
ALLOCATE(  pct(mpct)  ) 
ALLOCATE(  x(mx)  ) 
ALLOCATE(  y(my)  ) 
 
status = NF90_INQ_VARID(fidA,"DeltaT_basin_pc95",DeltaT_basin_pc95_ID); call erreur(status,.TRUE.,"inq_DeltaT_basin_pc95_ID")
status = NF90_INQ_VARID(fidA,"DeltaT_basin_pc50",DeltaT_basin_pc50_ID); call erreur(status,.TRUE.,"inq_DeltaT_basin_pc50_ID")
status = NF90_INQ_VARID(fidA,"DeltaT_basin_pc05",DeltaT_basin_pc05_ID); call erreur(status,.TRUE.,"inq_DeltaT_basin_pc05_ID")
status = NF90_INQ_VARID(fidA,"K0",K0_ID); call erreur(status,.TRUE.,"inq_K0_ID")
status = NF90_INQ_VARID(fidA,"pct",pct_ID); call erreur(status,.TRUE.,"inq_pct_ID")
status = NF90_INQ_VARID(fidA,"x",x_ID); call erreur(status,.TRUE.,"inq_x_ID")
status = NF90_INQ_VARID(fidA,"y",y_ID); call erreur(status,.TRUE.,"inq_y_ID")
 
status = NF90_GET_VAR(fidA,DeltaT_basin_pc95_ID,DeltaT_basin_pc95); call erreur(status,.TRUE.,"getvar_DeltaT_basin_pc95")
status = NF90_GET_VAR(fidA,DeltaT_basin_pc50_ID,DeltaT_basin_pc50); call erreur(status,.TRUE.,"getvar_DeltaT_basin_pc50")
status = NF90_GET_VAR(fidA,DeltaT_basin_pc05_ID,DeltaT_basin_pc05); call erreur(status,.TRUE.,"getvar_DeltaT_basin_pc05")
status = NF90_GET_VAR(fidA,K0_ID,K0); call erreur(status,.TRUE.,"getvar_K0")
status = NF90_GET_VAR(fidA,pct_ID,pct); call erreur(status,.TRUE.,"getvar_pct")
status = NF90_GET_VAR(fidA,x_ID,x); call erreur(status,.TRUE.,"getvar_x")
status = NF90_GET_VAR(fidA,y_ID,y); call erreur(status,.TRUE.,"getvar_y")
 
status = NF90_CLOSE(fidA); call erreur(status,.TRUE.,"close_file")
 
!---------------------------------------
 
write(*,*) 'Reading ', TRIM(file_bn)
 
status = NF90_OPEN(TRIM(file_bn),0,fidB); call erreur(status,.TRUE.,"read")
 
ALLOCATE(  basinNumber(mx,my)  ) 
 
status = NF90_INQ_VARID(fidB,"basinNumber",basinNumber_ID); call erreur(status,.TRUE.,"inq_basinNumber_ID")
 
status = NF90_GET_VAR(fidB,basinNumber_ID,basinNumber); call erreur(status,.TRUE.,"getvar_basinNumber")
 
status = NF90_CLOSE(fidB); call erreur(status,.TRUE.,"close_file")

!---------------------------------------
! Modification of the variables :

Nbasin=MAXVAL(basinNumber)+1

ALLOCATE( dT_median(Nbasin), dT_5thpct(Nbasin), dT_95thpct(Nbasin) )

! from readjust_deltaT_non_local_and_save_melt_MEDIAN.f90
dT_median(:) = (/ -0.15948466855502003, 0.80184953618220900, 9.6753354313715290e-2, 0.56328612973301229, 0.82906067997372057, 0.39593084278299073, 5.6393672343579837e-3, -0.14920448303891182, 0.40210906578254157, 1.0150054864802700, 8.3537124490397963e-2, -0.61512307466741678, 3.0671006394673017e-3, 1.7114340467726169e-2, -8.7542537953089417e-2, 0.11194347327395332 /)

! from readjust_deltaT_non_local_and_save_melt_5TH_PCT.f90 
dT_5thpct(:) = (/ -6.5043324647154455e-2, 1.0472582447705274, 0.20921525675042540, 0.76178908115001309, 1.1792981050454865, 0.57918066447514505, 0.15517749550279625, -9.8787773672779355e-2, 0.66582502282767730, 1.4957194779281202, 0.30697253392664681, -0.36818897749567703, 0.10526877945684110, 6.5159526911439558e-2, -1.9564068543699964e-2, 0.18252923064859006 /)

! from readjust_deltaT_non_local_and_save_melt_95TH_PCT.f90
dT_95thpct(:) = (/ -0.23630915758124554, 0.60702431713629579, 7.6838002034667618E-003, 0.40596305415191170, 0.55025710122782834, 0.24984633114643806, -0.11431269883131545, -0.18920082728732751, 0.19280805020297576, 0.63262598741837905, -9.3694969056708266E-002, -0.81118798605542719, -7.8252953569171502E-002, -2.0682286102327835E-002, -0.14119038769283254, 5.6703636134596835E-002 /)

DeltaT_basin_pc50(:,:,:) = -99999.9

do kbasin=1,Nbasin
 do ii=1,mx
 do jj=1,my
   if ( basinNumber(ii,jj)+1 .eq. kbasin ) then
     DeltaT_basin_pc50(ii,jj,1) = dT_5thpct(kbasin)
     DeltaT_basin_pc50(ii,jj,2) = dT_median(kbasin)
     DeltaT_basin_pc50(ii,jj,3) = dT_95thpct(kbasin)
   endif
 enddo
 enddo
enddo

!---------------------------------------
! Writing new netcdf file :
 
write(*,*) 'Creating ', TRIM(file_out)
 
status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidM,"pct",mpct,dimID_pct); call erreur(status,.TRUE.,"def_dimID_pct")
status = NF90_DEF_DIM(fidM,"x",mx,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")
status = NF90_DEF_DIM(fidM,"y",my,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
  
status = NF90_DEF_VAR(fidM,"DeltaT_basin_pc95",NF90_DOUBLE,(/dimID_x,dimID_y,dimID_pct/),DeltaT_basin_pc95_ID); call erreur(status,.TRUE.,"def_var_DeltaT_basin_pc95_ID")
status = NF90_DEF_VAR(fidM,"DeltaT_basin_pc50",NF90_DOUBLE,(/dimID_x,dimID_y,dimID_pct/),DeltaT_basin_pc50_ID); call erreur(status,.TRUE.,"def_var_DeltaT_basin_pc50_ID")
status = NF90_DEF_VAR(fidM,"DeltaT_basin_pc05",NF90_DOUBLE,(/dimID_x,dimID_y,dimID_pct/),DeltaT_basin_pc05_ID); call erreur(status,.TRUE.,"def_var_DeltaT_basin_pc05_ID")
status = NF90_DEF_VAR(fidM,"K0",NF90_DOUBLE,(/dimID_pct/),K0_ID); call erreur(status,.TRUE.,"def_var_K0_ID")
status = NF90_DEF_VAR(fidM,"pct",NF90_DOUBLE,(/dimID_pct/),pct_ID); call erreur(status,.TRUE.,"def_var_pct_ID")
status = NF90_DEF_VAR(fidM,"x",NF90_DOUBLE,(/dimID_x/),x_ID); call erreur(status,.TRUE.,"def_var_x_ID")
status = NF90_DEF_VAR(fidM,"y",NF90_DOUBLE,(/dimID_y/),y_ID); call erreur(status,.TRUE.,"def_var_y_ID")
 
status = NF90_PUT_ATT(fidM,DeltaT_basin_pc95_ID,"units","m.w.e / yr"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc95_ID")
status = NF90_PUT_ATT(fidM,DeltaT_basin_pc95_ID,"long_name","basin-95th-percentile DeltaT values"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc95_ID")
status = NF90_PUT_ATT(fidM,DeltaT_basin_pc50_ID,"units","m.w.e / yr"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc50_ID")
status = NF90_PUT_ATT(fidM,DeltaT_basin_pc50_ID,"long_name","basin-median DeltaT values"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc50_ID")
status = NF90_PUT_ATT(fidM,DeltaT_basin_pc05_ID,"units","m.w.e / yr"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc05_ID")
status = NF90_PUT_ATT(fidM,DeltaT_basin_pc05_ID,"long_name","basin-5th-percentile DeltaT values"); call erreur(status,.TRUE.,"put_att_DeltaT_basin_pc05_ID")
status = NF90_PUT_ATT(fidM,K0_ID,"units","-"); call erreur(status,.TRUE.,"put_att_K0_ID")
status = NF90_PUT_ATT(fidM,K0_ID,"long_name","K0 coefficient"); call erreur(status,.TRUE.,"put_att_K0_ID")
status = NF90_PUT_ATT(fidM,pct_ID,"units","-"); call erreur(status,.TRUE.,"put_att_pct_ID")
status = NF90_PUT_ATT(fidM,pct_ID,"long_name","K0 percentile"); call erreur(status,.TRUE.,"put_att_pct_ID")
 
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created on Wed 28 Nov 2018 16:26:22 +07"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
 
status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 
 
status = NF90_PUT_VAR(fidM,DeltaT_basin_pc95_ID,DeltaT_basin_pc95); call erreur(status,.TRUE.,"var_DeltaT_basin_pc95_ID")
status = NF90_PUT_VAR(fidM,DeltaT_basin_pc50_ID,DeltaT_basin_pc50); call erreur(status,.TRUE.,"var_DeltaT_basin_pc50_ID")
status = NF90_PUT_VAR(fidM,DeltaT_basin_pc05_ID,DeltaT_basin_pc05); call erreur(status,.TRUE.,"var_DeltaT_basin_pc05_ID")
status = NF90_PUT_VAR(fidM,K0_ID,K0); call erreur(status,.TRUE.,"var_K0_ID")
status = NF90_PUT_VAR(fidM,pct_ID,pct); call erreur(status,.TRUE.,"var_pct_ID")
status = NF90_PUT_VAR(fidM,x_ID,x); call erreur(status,.TRUE.,"var_x_ID")
status = NF90_PUT_VAR(fidM,y_ID,y); call erreur(status,.TRUE.,"var_y_ID")
 
status = NF90_CLOSE(fidM); call erreur(status,.TRUE.,"final")

end program modif



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
