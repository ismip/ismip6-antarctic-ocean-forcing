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
 
file_in  = 'coeff_K0_DeltaT_quadratic_local_Rignot_Depoorter_with_TS_uncert.nc'
file_bn  = 'basinNumbers_8km.nc'
file_out = 'coeff_K0_DeltaT_quadratic_local_Rignot_Depoorter_with_TS_uncert_ADJUSTED.nc'
 
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

! from readjust_deltaT_local_and_save_melt_MEDIAN.f90
dT_median(:) = (/ -0.43000000761821866, 0.91100001742597669, 0.11000000522471964, 0.67600001371465623, 1.0170000008074567, 0.49299998441711068, -9.0000114869326353e-3, -0.13599998992867768, 0.45999999216292053, 1.2219999753870070, -3.3999971230514348e-2, -0.92199998605065048, 6.1999981291592121e-2, 4.1000002529472113e-2, -8.6000003153458238e-2, 0.13099999667610973 /)

! from readjust_deltaT_local_and_save_melt_5TH_PCT.f90 
dT_5thpct(:) = (/ -0.24199999868869781, 1.1700000297278166, 0.23500001116190106, 0.88200002349913120, 1.3810000180965289, 0.68499999353662133, 0.16599999682512134, -7.9999987268820405e-2, 0.74900000588968396, 1.7309999995632097, 0.24700004211626947, -0.56999996933154762, 0.16799998632632196, 9.1000004904344678e-2, -8.9999994961544871e-3, 0.20800000033341348 /)

! from readjust_deltaT_local_and_save_melt_95TH_PCT.f90
dT_95thpct(:) = (/ -0.61000001616775990, 0.71400000806897879, 1.1000000522471964e-2, 0.52100000635255128, 0.74199998774565756, 0.34799997752998024, -0.14900001813657582, -0.17999999201856554, 0.23399998142849654, 0.83499995700549334, -0.27199998253490776, -1.2079999996349216, -1.8000022508203983e-2, 3.0000007245689631e-3, -0.14800000609830022, 7.1999993873760104e-2 /)

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
