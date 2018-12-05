program modif
 
USE netcdf
 
IMPLICIT NONE
 
INTEGER :: fidWOA, status, dimID_nbounds, dimID_z, dimID_y, dimID_x, mnbounds, mz, my, mx, thermal_forcing_ID, z_bnds_ID, z_ID, y_ID, x_ID, fidDFT, my2, mx2, thickness_ID, surface_ID, fidM, ii, jj, kk, kinf, ksup
 
CHARACTER(LEN=150) :: file_WOA, file_out, file_DFT
 
REAL*8,ALLOCATABLE,DIMENSION(:) :: x, y, z

REAL*8,ALLOCATABLE,DIMENSION(:,:) :: thickness, surface, ice_draft, z_bnds, TF_draft
 
REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: thermal_forcing

LOGICAL :: ll_z_upward
 
file_DFT  = 'bedmap2_8km.nc'
!file_WOA  = 'woa_meop_thermal_forcing_1955-2012_8km_x_60m.nc'
file_WOA  = 'obs_thermal_forcing_1995-2017_8km_x_60m.nc'
file_out  = 'obs_thermal_forcing_on_bedmap2_ice_draft.nc'
 
!---------------------------------------
! Read 3D thermal forcing :
 
write(*,*) 'Reading ', TRIM(file_WOA)
 
status = NF90_OPEN(TRIM(file_WOA),0,fidWOA); call erreur(status,.TRUE.,"read")
 
status = NF90_INQ_DIMID(fidWOA,"nbounds",dimID_nbounds); call erreur(status,.TRUE.,"inq_dimID_nbounds")
status = NF90_INQ_DIMID(fidWOA,"z",dimID_z); call erreur(status,.TRUE.,"inq_dimID_z")
status = NF90_INQ_DIMID(fidWOA,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidWOA,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
 
status = NF90_INQUIRE_DIMENSION(fidWOA,dimID_nbounds,len=mnbounds); call erreur(status,.TRUE.,"inq_dim_nbounds")
status = NF90_INQUIRE_DIMENSION(fidWOA,dimID_z,len=mz); call erreur(status,.TRUE.,"inq_dim_z")
status = NF90_INQUIRE_DIMENSION(fidWOA,dimID_y,len=my); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidWOA,dimID_x,len=mx); call erreur(status,.TRUE.,"inq_dim_x")
  
ALLOCATE(  thermal_forcing(my,mx,mz)  ) 
ALLOCATE(  TF_draft(my,mx)  ) 
ALLOCATE(  z_bnds(mnbounds,mz)  ) 
ALLOCATE(  z(mz)  ) 
ALLOCATE(  y(my)  ) 
ALLOCATE(  x(mx)  ) 
 
status = NF90_INQ_VARID(fidWOA,"thermal_forcing",thermal_forcing_ID); call erreur(status,.TRUE.,"inq_thermal_forcing_ID")
status = NF90_INQ_VARID(fidWOA,"z_bnds",z_bnds_ID); call erreur(status,.TRUE.,"inq_z_bnds_ID")
status = NF90_INQ_VARID(fidWOA,"z",z_ID); call erreur(status,.TRUE.,"inq_z_ID")
status = NF90_INQ_VARID(fidWOA,"y",y_ID); call erreur(status,.TRUE.,"inq_y_ID")
status = NF90_INQ_VARID(fidWOA,"x",x_ID); call erreur(status,.TRUE.,"inq_x_ID")
 
status = NF90_GET_VAR(fidWOA,thermal_forcing_ID,thermal_forcing); call erreur(status,.TRUE.,"getvar_thermal_forcing")
status = NF90_GET_VAR(fidWOA,z_bnds_ID,z_bnds); call erreur(status,.TRUE.,"getvar_z_bnds")
status = NF90_GET_VAR(fidWOA,z_ID,z); call erreur(status,.TRUE.,"getvar_z")
status = NF90_GET_VAR(fidWOA,y_ID,y); call erreur(status,.TRUE.,"getvar_y")
status = NF90_GET_VAR(fidWOA,x_ID,x); call erreur(status,.TRUE.,"getvar_x")
 
status = NF90_CLOSE(fidWOA); call erreur(status,.TRUE.,"close_file")

!---------------------------------------
! Read BEDMAP2's ice draft :
 
write(*,*) 'Reading ', TRIM(file_DFT)
 
status = NF90_OPEN(TRIM(file_DFT),0,fidDFT); call erreur(status,.TRUE.,"read")
 
status = NF90_INQ_DIMID(fidDFT,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidDFT,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
 
status = NF90_INQUIRE_DIMENSION(fidDFT,dimID_y,len=my2); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidDFT,dimID_x,len=mx2); call erreur(status,.TRUE.,"inq_dim_x")
  
if ( mx .ne. mx2 .or. my .ne. my2 ) then
  write(*,*) '~!@#$%^* ERROR : dimension mismatch >>>>> stop'
  stop
endif

ALLOCATE(  thickness(mx,my)  ) 
ALLOCATE(  surface(mx,my)  ) 
ALLOCATE(  ice_draft(mx,my)  ) 
 
status = NF90_INQ_VARID(fidDFT,"thickness",thickness_ID); call erreur(status,.TRUE.,"inq_thickness_ID")
status = NF90_INQ_VARID(fidDFT,"surface",surface_ID); call erreur(status,.TRUE.,"inq_surface_ID")
 
status = NF90_GET_VAR(fidDFT,thickness_ID,thickness); call erreur(status,.TRUE.,"getvar_thickness")
status = NF90_GET_VAR(fidDFT,surface_ID,surface); call erreur(status,.TRUE.,"getvar_surface")
 
status = NF90_CLOSE(fidDFT); call erreur(status,.TRUE.,"close_file")

ice_draft(:,:) = surface(:,:) - thickness(:,:)
DEALLOCATE( surface, thickness )
 
!----------------------------------------------------------------------
! Linear interpolation of the thermal forcing on the ice draft depth :

if ( z(2) .gt. z(1) ) then
  ll_z_upward = .true.
else
  ll_z_upward = .false.
endif


do ii=1,mx
do jj=1,my

  if ( ice_draft(ii,jj) .lt. 0.d0 ) then

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
    write(*,*) 'zinf, ice_draft, zsup = ', z(kinf), ice_draft(ii,jj), z(ksup)

    TF_draft(ii,jj) = (   (z(ksup)-ice_draft(ii,jj)) * thermal_forcing(ii,jj,kinf)   &
    &                   + (ice_draft(ii,jj)-z(kinf)) * thermal_forcing(ii,jj,ksup) ) / (z(ksup)-z(kinf))

    write(*,555) ii, jj, kinf, ksup, z(kinf), z(ksup), ice_draft(ii,jj), thermal_forcing(ii,jj,kinf), thermal_forcing(ii,jj,ksup), TF_draft(ii,jj)

  else

    TF_draft(ii,jj) = 0.d0 ! will ensure zero melt in case of unexpected overlap
    !write(*,*) ii, jj, '........'

  endif

enddo
enddo
 
555 FORMAT(i4,i4,i3,i3,' ',f7.1,' ',f7.1,' ',f7.1,' ',f5.2,' ',f5.2,' ',f5.2)

!---------------------------------------
! Writing new netcdf file :
 
write(*,*) 'Creating ', TRIM(file_out)
 
status = NF90_CREATE(TRIM(file_out),NF90_NOCLOBBER,fidM); call erreur(status,.TRUE.,'create')
 
status = NF90_DEF_DIM(fidM,"y",my,dimID_y); call erreur(status,.TRUE.,"def_dimID_y")
status = NF90_DEF_DIM(fidM,"x",mx,dimID_x); call erreur(status,.TRUE.,"def_dimID_x")
  
status = NF90_DEF_VAR(fidM,"TF_draft",NF90_DOUBLE,(/dimID_y,dimID_x/),thermal_forcing_ID); call erreur(status,.TRUE.,"def_var_thermal_forcing_ID")
status = NF90_DEF_VAR(fidM,"y",NF90_DOUBLE,(/dimID_y/),y_ID); call erreur(status,.TRUE.,"def_var_y_ID")
status = NF90_DEF_VAR(fidM,"x",NF90_DOUBLE,(/dimID_x/),x_ID); call erreur(status,.TRUE.,"def_var_x_ID")
 
status = NF90_PUT_ATT(fidM,thermal_forcing_ID,"long_name","thermal forcing on the ice draft"); call erreur(status,.TRUE.,"put_att_thermal_forcing_ID")
status = NF90_PUT_ATT(fidM,thermal_forcing_ID,"units","degrees_celsius"); call erreur(status,.TRUE.,"put_att_thermal_forcing_ID")
!status = NF90_PUT_ATT(fidM,thermal_forcing_ID,"_FillValue",0.d0); call erreur(status,.TRUE.,"put_att_thermal_forcing_ID")
status = NF90_PUT_ATT(fidM,y_ID,"units","meters"); call erreur(status,.TRUE.,"put_att_y_ID")
!status = NF90_PUT_ATT(fidM,y_ID,"_FillValue",0.d0); call erreur(status,.TRUE.,"put_att_y_ID")
status = NF90_PUT_ATT(fidM,x_ID,"units","meters"); call erreur(status,.TRUE.,"put_att_x_ID")
!status = NF90_PUT_ATT(fidM,x_ID,"_FillValue",0.d0); call erreur(status,.TRUE.,"put_att_x_ID")
 
status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using interpolate_TF_on_bedmap2_ice_draft.f90"); call erreur(status,.TRUE.,"put_att_GLOBAL_ID")
 
status = NF90_ENDDEF(fidM); call erreur(status,.TRUE.,"fin_definition") 
 
status = NF90_PUT_VAR(fidM,thermal_forcing_ID,TF_draft); call erreur(status,.TRUE.,"var_thermal_forcing_ID")
status = NF90_PUT_VAR(fidM,y_ID,y); call erreur(status,.TRUE.,"var_y_ID")
status = NF90_PUT_VAR(fidM,x_ID,x); call erreur(status,.TRUE.,"var_x_ID")
 
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
