program modif                                         

USE netcdf                                            
USE polarstereo

IMPLICIT NONE                                         

INTEGER :: fidA, status, dimID_x, dimID_y, lon_ID, lat_ID, bedrock_topography_ID, x_ID, y_ID, bed_ID, surface_ID,    &
&          fidM, isfmask_ID, fidDFT, ice_base_topography_ID, fidSRF, surface_elevation_ID, ii, jj, Nisf, dimID_Nisf, &
&          name_isf_ID, name_reg_ID, strlen, dimID_strlen, iim1, iip1, jjm1, jjp1, kisf, GL_mask_ID, IF_mask_ID,     &
&          front_bot_dep_avg_ID, front_ice_dep_avg_ID, front_min_lon_ID, front_max_lon_ID, front_min_lat_ID,         &
&          front_max_lat_ID, front_bot_dep_max_ID, front_ice_dep_min_ID, melt_isf_ID, err_melt_isf_ID, mx, my,       &
&          open_ocean_mask_ID, icemask_shelves_ID, icemask_grounded_ID, thickness_ID, rs, imin, imax, jmin, jmax,    &
&          area_isf_ID

CHARACTER(LEN=150) :: file_ice, file_msk

CHARACTER(LEN=19),ALLOCATABLE,DIMENSION(:) :: name_isf, name_reg  !! adapt strlen if (LEN=19) is modified.

INTEGER*2,ALLOCATABLE,DIMENSION(:,:) :: isfmask,  & ! 0 grounded; 1 ocean; >1 ice shelf number
&                                       GLINE,    & ! 0 or ice shelf number if adjacent to grounding line
&                                       FRONT,    & ! 0 or ice shelf number if adjacent to ice shelf front
&                                       tmppin, rmppin

REAL*8,ALLOCATABLE,DIMENSION(:) :: x, y  

REAL*8,ALLOCATABLE,DIMENSION(:,:) :: x2d, y2d, lon, lat, bedrock_topography, ice_base_topography, surface_elevation, thickness, open_ocean_mask, icemask_shelves, icemask_grounded 

REAL*8,ALLOCATABLE,DIMENSION(:) :: front_bot_dep_avg, front_ice_dep_avg

REAL*8 :: rn_ice

REAL*8,ALLOCATABLE,DIMENSION(:) :: front_min_lon, front_max_lon, front_min_lat, front_max_lat, front_bot_dep_max, front_ice_dep_min, &
&                                  err_melt_isf, melt_isf, our_isf_area, area_isf

REAL*8 :: tmplon, itmpmax, itmpmin, jtmpmin, jtmpmax, mesh_area

LOGICAL :: llisf

file_ice = 'bedmap2_8km.nc'
file_msk = 'ice_shelf_mask_bedmap2_8km.nc'

!---------------------------------------                   
! Read bedrock, srf, thickness :

write(*,*) 'Reading ', TRIM(file_ice)
 
status = NF90_OPEN(TRIM(file_ice),0,fidA); call erreur(status,.TRUE.,"read")
 
status = NF90_INQ_DIMID(fidA,"y",dimID_y); call erreur(status,.TRUE.,"inq_dimID_y")
status = NF90_INQ_DIMID(fidA,"x",dimID_x); call erreur(status,.TRUE.,"inq_dimID_x")
 
status = NF90_INQUIRE_DIMENSION(fidA,dimID_y,len=my); call erreur(status,.TRUE.,"inq_dim_y")
status = NF90_INQUIRE_DIMENSION(fidA,dimID_x,len=mx); call erreur(status,.TRUE.,"inq_dim_x")
  
ALLOCATE(  open_ocean_mask(mx,my)  )
ALLOCATE(  icemask_shelves(mx,my)  )
ALLOCATE(  icemask_grounded(mx,my)  )
!ALLOCATE(  thickness(mx,my)  ) 
!ALLOCATE(  surface(mx,my)  ) 
ALLOCATE(  surface_elevation(mx,my)  )
ALLOCATE(  thickness(mx,my)  )
ALLOCATE(  bedrock_topography(mx,my)  ) 
ALLOCATE(  ice_base_topography(mx,my)  ) 
ALLOCATE(  y(my)  ) 
ALLOCATE(  x(mx)  ) 
ALLOCATE( isfmask(mx,my) )
 
status = NF90_INQ_VARID(fidA,"icemask_shelves",icemask_shelves_ID); call erreur(status,.TRUE.,"inq_icemask_shelves_ID")
status = NF90_INQ_VARID(fidA,"open_ocean_mask",open_ocean_mask_ID); call erreur(status,.TRUE.,"inq_open_ocean_mask_ID")
status = NF90_INQ_VARID(fidA,"icemask_grounded",icemask_grounded_ID); call erreur(status,.TRUE.,"inq_icemask_grounded_ID")
status = NF90_INQ_VARID(fidA,"thickness",thickness_ID); call erreur(status,.TRUE.,"inq_thickness_ID")
status = NF90_INQ_VARID(fidA,"surface",surface_ID); call erreur(status,.TRUE.,"inq_surface_ID")
status = NF90_INQ_VARID(fidA,"bed",bed_ID); call erreur(status,.TRUE.,"inq_bed_ID")
status = NF90_INQ_VARID(fidA,"y",y_ID); call erreur(status,.TRUE.,"inq_y_ID")
status = NF90_INQ_VARID(fidA,"x",x_ID); call erreur(status,.TRUE.,"inq_x_ID")
 
status = NF90_GET_VAR(fidA,icemask_shelves_ID,icemask_shelves); call erreur(status,.TRUE.,"getvar_icemask_shelves")
status = NF90_GET_VAR(fidA,icemask_grounded_ID,icemask_grounded); call erreur(status,.TRUE.,"getvar_icemask_grounded")
status = NF90_GET_VAR(fidA,open_ocean_mask_ID,open_ocean_mask); call erreur(status,.TRUE.,"getvar_open_ocean_mask")
status = NF90_GET_VAR(fidA,thickness_ID,thickness); call erreur(status,.TRUE.,"getvar_thickness")
status = NF90_GET_VAR(fidA,surface_ID,surface_elevation); call erreur(status,.TRUE.,"getvar_surface")
status = NF90_GET_VAR(fidA,bed_ID,bedrock_topography); call erreur(status,.TRUE.,"getvar_bed")
status = NF90_GET_VAR(fidA,y_ID,y); call erreur(status,.TRUE.,"getvar_y")
status = NF90_GET_VAR(fidA,x_ID,x); call erreur(status,.TRUE.,"getvar_x")
 
status = NF90_CLOSE(fidA); call erreur(status,.TRUE.,"close_file")

ice_base_topography(:,:) = surface_elevation(:,:) - thickness(:,:)
DEALLOCATE( surface_elevation, thickness )

!---------------------------------------
! Get lon/lat

ALLOCATE(  lon(mx,my), lat(mx,my) )
ALLOCATE(  x2d(mx,my), y2d(mx,my) )
do ii=1,mx
  y2d(ii,:) = y(:)
enddo
do jj=1,my
  x2d(:,jj) = x(:)
enddo
CALL stereo_loc_xy_to_latlon(mx,my,x2d,y2d,lat,lon)
DEALLOCATE( x2d, y2d ) 

mesh_area = abs( (x(3)-x(2))*(y(3)-y(2)) ) * 1.e-6  !! in km^2

lon(381,1:380) = 180.0 ! to avoid interpolation to zero between -180 and 180

!---------------------------------------                      
! Mask calculation

write(*,*) 'Calculating mask...'

do ii=1,mx
do jj=1,my

  if ( icemask_shelves(ii,jj) .ge. 5.d-1 ) then

    isfmask(ii,jj) = 2  ! floating ice

    ! Now attribute an ID per ice shelf :
    if ( lon(ii,jj) .ge.  -62.7 .and. lon(ii,jj) .lt.  -61.75.and. lat(ii,jj) .ge.  -66.16 .and. lat(ii,jj) .lt.  -65.55 )  isfmask(ii,jj) = 26 ! Larsen B
    if ( lon(ii,jj) .ge.  -61.75.and. lon(ii,jj) .lt.  -61.2 .and. lat(ii,jj) .ge.  -66.10 .and. lat(ii,jj) .lt.  -65.49 )  isfmask(ii,jj) = 26 ! Larsen B
    if ( lon(ii,jj) .ge.  -61.2 .and. lon(ii,jj) .lt.  -60.5 .and. lat(ii,jj) .ge.  -66.00 .and. lat(ii,jj) .lt.  -65.49 )  isfmask(ii,jj) = 26 ! Larsen B
    if ( lon(ii,jj) .ge.  -64.0 .and. lon(ii,jj) .lt.  -58.0 .and. lat(ii,jj) .ge.  -73.10 .and. lat(ii,jj) .lt.  -68.43 )  isfmask(ii,jj) = 55 ! Larsen D ##NB## keep before Larsen C
    if ( lon(ii,jj) .ge.  -66.0 .and. lon(ii,jj) .lt.  -65.0 .and. lat(ii,jj) .ge.  -68.90 .and. lat(ii,jj) .lt.  -67.00 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -65.0 .and. lon(ii,jj) .lt.  -64.0 .and. lat(ii,jj) .ge.  -68.90 .and. lat(ii,jj) .lt.  -66.40 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -64.0 .and. lon(ii,jj) .lt.  -63.7 .and. lat(ii,jj) .ge.  -68.52 .and. lat(ii,jj) .lt.  -66.15 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -63.7 .and. lon(ii,jj) .lt.  -63.4 .and. lat(ii,jj) .ge.  -68.47 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -63.4 .and. lon(ii,jj) .lt.  -62.9 .and. lat(ii,jj) .ge.  -68.43 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -62.9 .and. lon(ii,jj) .lt.  -62.8 .and. lat(ii,jj) .ge.  -68.46 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -62.8 .and. lon(ii,jj) .lt.  -62.7 .and. lat(ii,jj) .ge.  -68.47 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -62.7 .and. lon(ii,jj) .lt.  -62.6 .and. lat(ii,jj) .ge.  -68.48 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -62.6 .and. lon(ii,jj) .lt.  -62.5 .and. lat(ii,jj) .ge.  -68.50 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -62.5 .and. lon(ii,jj) .lt.  -62.4 .and. lat(ii,jj) .ge.  -68.51 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -62.4 .and. lon(ii,jj) .lt.  -62.3 .and. lat(ii,jj) .ge.  -68.52 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -62.3 .and. lon(ii,jj) .lt.  -62.2 .and. lat(ii,jj) .ge.  -68.54 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -62.2 .and. lon(ii,jj) .lt.  -62.1 .and. lat(ii,jj) .ge.  -68.55 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -62.1 .and. lon(ii,jj) .lt.  -62.0 .and. lat(ii,jj) .ge.  -68.57 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -62.0 .and. lon(ii,jj) .lt.  -61.9 .and. lat(ii,jj) .ge.  -68.58 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -61.9 .and. lon(ii,jj) .lt.  -61.8 .and. lat(ii,jj) .ge.  -68.59 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -61.8 .and. lon(ii,jj) .lt.  -61.7 .and. lat(ii,jj) .ge.  -68.61 .and. lat(ii,jj) .lt.  -66.16 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -61.7 .and. lon(ii,jj) .lt.  -61.6 .and. lat(ii,jj) .ge.  -68.62 .and. lat(ii,jj) .lt.  -66.10 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -61.6 .and. lon(ii,jj) .lt.  -61.5 .and. lat(ii,jj) .ge.  -68.63 .and. lat(ii,jj) .lt.  -66.10 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -61.5 .and. lon(ii,jj) .lt.  -61.4 .and. lat(ii,jj) .ge.  -68.65 .and. lat(ii,jj) .lt.  -66.10 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -61.4 .and. lon(ii,jj) .lt.  -61.3 .and. lat(ii,jj) .ge.  -68.66 .and. lat(ii,jj) .lt.  -66.10 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -61.3 .and. lon(ii,jj) .lt.  -61.2 .and. lat(ii,jj) .ge.  -68.67 .and. lat(ii,jj) .lt.  -66.10 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -61.2 .and. lon(ii,jj) .lt.  -61.1 .and. lat(ii,jj) .ge.  -68.69 .and. lat(ii,jj) .lt.  -66.10 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -61.1 .and. lon(ii,jj) .lt.  -59.0 .and. lat(ii,jj) .ge.  -68.70 .and. lat(ii,jj) .lt.  -66.00 )  isfmask(ii,jj) = 42 ! Larsen C
    if ( lon(ii,jj) .ge.  -63.0 .and. lon(ii,jj) .lt.  -59.0 .and. lat(ii,jj) .ge.  -73.60 .and. lat(ii,jj) .lt.  -73.10 )  isfmask(ii,jj) = 68 ! Larsen E
    if ( lon(ii,jj) .ge.  -62.0 .and. lon(ii,jj) .lt.  -60.0 .and. lat(ii,jj) .ge.  -74.35 .and. lat(ii,jj) .lt.  -73.80 )  isfmask(ii,jj) = 60 ! Larsen F
    if ( lon(ii,jj) .ge.  -62.5 .and. lon(ii,jj) .lt.  -61.2 .and. lat(ii,jj) .ge.  -74.85 .and. lat(ii,jj) .lt.  -74.45 )  isfmask(ii,jj) = 27 ! Larsen G
    if ( lon(ii,jj) .ge.  -86.0 .and. lon(ii,jj) .lt.  -60.6 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -60.6 .and. lon(ii,jj) .lt.  -60.5 .and. lat(ii,jj) .ge.  -83.26 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -60.5 .and. lon(ii,jj) .lt.  -60.4 .and. lat(ii,jj) .ge.  -83.22 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -60.4 .and. lon(ii,jj) .lt.  -60.3 .and. lat(ii,jj) .ge.  -83.19 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -60.3 .and. lon(ii,jj) .lt.  -60.2 .and. lat(ii,jj) .ge.  -83.15 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -60.2 .and. lon(ii,jj) .lt.  -60.1 .and. lat(ii,jj) .ge.  -83.11 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -60.1 .and. lon(ii,jj) .lt.  -60.0 .and. lat(ii,jj) .ge.  -83.07 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -60.0 .and. lon(ii,jj) .lt.  -59.9 .and. lat(ii,jj) .ge.  -83.03 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -59.9 .and. lon(ii,jj) .lt.  -59.8 .and. lat(ii,jj) .ge.  -82.99 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -59.8 .and. lon(ii,jj) .lt.  -59.7 .and. lat(ii,jj) .ge.  -82.96 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -59.7 .and. lon(ii,jj) .lt.  -59.6 .and. lat(ii,jj) .ge.  -82.92 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -59.6 .and. lon(ii,jj) .lt.  -59.5 .and. lat(ii,jj) .ge.  -82.88 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -59.5 .and. lon(ii,jj) .lt.  -59.4 .and. lat(ii,jj) .ge.  -82.84 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -59.4 .and. lon(ii,jj) .lt.  -59.3 .and. lat(ii,jj) .ge.  -82.80 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -59.3 .and. lon(ii,jj) .lt.  -59.2 .and. lat(ii,jj) .ge.  -82.76 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -59.2 .and. lon(ii,jj) .lt.  -59.1 .and. lat(ii,jj) .ge.  -82.73 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -59.1 .and. lon(ii,jj) .lt.  -59.0 .and. lat(ii,jj) .ge.  -82.69 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -59.0 .and. lon(ii,jj) .lt.  -58.9 .and. lat(ii,jj) .ge.  -82.65 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -58.9 .and. lon(ii,jj) .lt.  -58.8 .and. lat(ii,jj) .ge.  -82.61 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -58.8 .and. lon(ii,jj) .lt.  -58.7 .and. lat(ii,jj) .ge.  -82.57 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -58.7 .and. lon(ii,jj) .lt.  -58.6 .and. lat(ii,jj) .ge.  -82.53 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -58.6 .and. lon(ii,jj) .lt.  -58.5 .and. lat(ii,jj) .ge.  -82.50 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -58.5 .and. lon(ii,jj) .lt.  -58.4 .and. lat(ii,jj) .ge.  -82.46 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -58.4 .and. lon(ii,jj) .lt.  -58.3 .and. lat(ii,jj) .ge.  -82.42 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -58.3 .and. lon(ii,jj) .lt.  -58.2 .and. lat(ii,jj) .ge.  -82.38 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -58.2 .and. lon(ii,jj) .lt.  -58.1 .and. lat(ii,jj) .ge.  -82.34 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -58.1 .and. lon(ii,jj) .lt.  -58.0 .and. lat(ii,jj) .ge.  -82.30 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -58.0 .and. lon(ii,jj) .lt.  -57.9 .and. lat(ii,jj) .ge.  -82.27 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -57.9 .and. lon(ii,jj) .lt.  -57.8 .and. lat(ii,jj) .ge.  -82.23 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -57.8 .and. lon(ii,jj) .lt.  -57.7 .and. lat(ii,jj) .ge.  -82.19 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -57.7 .and. lon(ii,jj) .lt.  -57.6 .and. lat(ii,jj) .ge.  -82.15 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -57.6 .and. lon(ii,jj) .lt.  -57.5 .and. lat(ii,jj) .ge.  -82.11 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -57.5 .and. lon(ii,jj) .lt.  -57.4 .and. lat(ii,jj) .ge.  -82.07 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -57.4 .and. lon(ii,jj) .lt.  -57.3 .and. lat(ii,jj) .ge.  -82.04 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -57.3 .and. lon(ii,jj) .lt.  -57.2 .and. lat(ii,jj) .ge.  -82.00 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -57.2 .and. lon(ii,jj) .lt.  -57.1 .and. lat(ii,jj) .ge.  -81.96 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -57.1 .and. lon(ii,jj) .lt.  -57.0 .and. lat(ii,jj) .ge.  -81.92 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -57.0 .and. lon(ii,jj) .lt.  -56.9 .and. lat(ii,jj) .ge.  -81.88 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -56.9 .and. lon(ii,jj) .lt.  -56.8 .and. lat(ii,jj) .ge.  -81.84 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -56.8 .and. lon(ii,jj) .lt.  -56.7 .and. lat(ii,jj) .ge.  -81.81 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -56.7 .and. lon(ii,jj) .lt.  -56.6 .and. lat(ii,jj) .ge.  -81.77 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -56.6 .and. lon(ii,jj) .lt.  -56.5 .and. lat(ii,jj) .ge.  -81.73 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -56.5 .and. lon(ii,jj) .lt.  -56.4 .and. lat(ii,jj) .ge.  -81.69 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -56.4 .and. lon(ii,jj) .lt.  -56.3 .and. lat(ii,jj) .ge.  -81.65 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -56.3 .and. lon(ii,jj) .lt.  -56.2 .and. lat(ii,jj) .ge.  -81.61 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -56.2 .and. lon(ii,jj) .lt.  -56.1 .and. lat(ii,jj) .ge.  -81.58 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -56.1 .and. lon(ii,jj) .lt.  -56.0 .and. lat(ii,jj) .ge.  -81.54 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -56.0 .and. lon(ii,jj) .lt.  -55.9 .and. lat(ii,jj) .ge.  -81.50 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -55.9 .and. lon(ii,jj) .lt.  -55.8 .and. lat(ii,jj) .ge.  -81.46 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -55.8 .and. lon(ii,jj) .lt.  -55.7 .and. lat(ii,jj) .ge.  -81.42 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -55.7 .and. lon(ii,jj) .lt.  -55.6 .and. lat(ii,jj) .ge.  -81.38 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -55.6 .and. lon(ii,jj) .lt.  -55.5 .and. lat(ii,jj) .ge.  -81.35 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -55.5 .and. lon(ii,jj) .lt.  -55.4 .and. lat(ii,jj) .ge.  -81.31 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -55.4 .and. lon(ii,jj) .lt.  -55.3 .and. lat(ii,jj) .ge.  -81.27 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -55.3 .and. lon(ii,jj) .lt.  -55.2 .and. lat(ii,jj) .ge.  -81.23 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -55.2 .and. lon(ii,jj) .lt.  -55.1 .and. lat(ii,jj) .ge.  -81.19 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -55.1 .and. lon(ii,jj) .lt.  -55.0 .and. lat(ii,jj) .ge.  -81.15 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -55.0 .and. lon(ii,jj) .lt.  -54.9 .and. lat(ii,jj) .ge.  -81.12 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -54.9 .and. lon(ii,jj) .lt.  -54.8 .and. lat(ii,jj) .ge.  -81.08 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -54.8 .and. lon(ii,jj) .lt.  -54.7 .and. lat(ii,jj) .ge.  -81.04 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -54.7 .and. lon(ii,jj) .lt.  -54.2 .and. lat(ii,jj) .ge.  -81.00 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -54.2 .and. lon(ii,jj) .lt.  -52.0 .and. lat(ii,jj) .ge.  -80.70 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -52.0 .and. lon(ii,jj) .lt.  -47.5 .and. lat(ii,jj) .ge.  -80.30 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -64.0 .and. lon(ii,jj) .lt.  -62.5 .and. lat(ii,jj) .ge.  -74.85 .and. lat(ii,jj) .lt.  -74.45 )  isfmask(ii,jj) = 21
    if ( lon(ii,jj) .ge.  -60.6 .and. lon(ii,jj) .lt.  -60.5 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -83.26 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -60.5 .and. lon(ii,jj) .lt.  -60.4 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -83.22 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -60.4 .and. lon(ii,jj) .lt.  -60.3 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -83.19 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -60.3 .and. lon(ii,jj) .lt.  -60.2 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -83.15 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -60.2 .and. lon(ii,jj) .lt.  -60.1 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -83.11 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -60.1 .and. lon(ii,jj) .lt.  -60.0 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -83.07 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -60.0 .and. lon(ii,jj) .lt.  -59.9 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -83.03 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -59.9 .and. lon(ii,jj) .lt.  -59.8 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.99 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -59.8 .and. lon(ii,jj) .lt.  -59.7 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.96 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -59.7 .and. lon(ii,jj) .lt.  -59.6 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.92 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -59.6 .and. lon(ii,jj) .lt.  -59.5 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.88 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -59.5 .and. lon(ii,jj) .lt.  -59.4 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.84 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -59.4 .and. lon(ii,jj) .lt.  -59.3 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.80 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -59.3 .and. lon(ii,jj) .lt.  -59.2 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.76 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -59.2 .and. lon(ii,jj) .lt.  -59.1 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.73 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -59.1 .and. lon(ii,jj) .lt.  -59.0 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.69 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -59.0 .and. lon(ii,jj) .lt.  -58.9 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.65 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -58.9 .and. lon(ii,jj) .lt.  -58.8 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.61 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -58.8 .and. lon(ii,jj) .lt.  -58.7 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.57 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -58.7 .and. lon(ii,jj) .lt.  -58.6 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.53 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -58.6 .and. lon(ii,jj) .lt.  -58.5 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.50 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -58.5 .and. lon(ii,jj) .lt.  -58.4 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.46 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -58.4 .and. lon(ii,jj) .lt.  -58.3 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.42 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -58.3 .and. lon(ii,jj) .lt.  -58.2 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.38 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -58.2 .and. lon(ii,jj) .lt.  -58.1 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.34 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -58.1 .and. lon(ii,jj) .lt.  -58.0 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.30 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -58.0 .and. lon(ii,jj) .lt.  -57.9 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.27 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -57.9 .and. lon(ii,jj) .lt.  -57.8 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.23 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -57.8 .and. lon(ii,jj) .lt.  -57.7 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.19 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -57.7 .and. lon(ii,jj) .lt.  -57.6 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.15 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -57.6 .and. lon(ii,jj) .lt.  -57.5 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.11 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -57.5 .and. lon(ii,jj) .lt.  -57.4 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.07 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -57.4 .and. lon(ii,jj) .lt.  -57.3 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.04 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -57.3 .and. lon(ii,jj) .lt.  -57.2 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -82.00 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -57.2 .and. lon(ii,jj) .lt.  -57.1 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.96 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -57.1 .and. lon(ii,jj) .lt.  -57.0 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.92 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -57.0 .and. lon(ii,jj) .lt.  -56.9 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.88 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -56.9 .and. lon(ii,jj) .lt.  -56.8 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.84 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -56.8 .and. lon(ii,jj) .lt.  -56.7 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.81 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -56.7 .and. lon(ii,jj) .lt.  -56.6 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.77 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -56.6 .and. lon(ii,jj) .lt.  -56.5 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.73 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -56.5 .and. lon(ii,jj) .lt.  -56.4 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.69 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -56.4 .and. lon(ii,jj) .lt.  -56.3 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.65 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -56.3 .and. lon(ii,jj) .lt.  -56.2 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.61 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -56.2 .and. lon(ii,jj) .lt.  -56.1 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.58 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -56.1 .and. lon(ii,jj) .lt.  -56.0 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.54 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -56.0 .and. lon(ii,jj) .lt.  -55.9 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.50 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -55.9 .and. lon(ii,jj) .lt.  -55.8 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.46 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -55.8 .and. lon(ii,jj) .lt.  -55.7 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.42 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -55.7 .and. lon(ii,jj) .lt.  -55.6 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.38 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -55.6 .and. lon(ii,jj) .lt.  -55.5 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.35 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -55.5 .and. lon(ii,jj) .lt.  -55.4 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.31 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -55.4 .and. lon(ii,jj) .lt.  -55.3 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.27 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -55.3 .and. lon(ii,jj) .lt.  -55.2 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.23 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -55.2 .and. lon(ii,jj) .lt.  -55.1 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.19 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -55.1 .and. lon(ii,jj) .lt.  -55.0 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.15 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -55.0 .and. lon(ii,jj) .lt.  -54.9 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.12 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -54.9 .and. lon(ii,jj) .lt.  -54.8 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.08 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -54.8 .and. lon(ii,jj) .lt.  -54.7 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.04 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -54.7 .and. lon(ii,jj) .lt.  -54.2 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -81.00 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -52.0 .and. lon(ii,jj) .lt.  -47.5 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -80.30 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -47.5 .and. lon(ii,jj) .lt.  -28.0 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -78.00 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -54.2 .and. lon(ii,jj) .lt.  -52.0 .and. lat(ii,jj) .ge.  -84.00 .and. lat(ii,jj) .lt.  -80.70 )  isfmask(ii,jj) = 11
    if ( lon(ii,jj) .ge.  -29.0 .and. lon(ii,jj) .lt.  -20.8 .and. lat(ii,jj) .ge.  -76.05 .and. lat(ii,jj) .lt.  -73.70 )  isfmask(ii,jj) = 69 ! Stancomb Brunt
    if ( lon(ii,jj) .ge.  -20.8 .and. lon(ii,jj) .lt.  -16.0 .and. lat(ii,jj) .ge.  -76.05 .and. lat(ii,jj) .lt.  -74.45 )  isfmask(ii,jj) = 69 ! Stancomb Brunt
    if ( lon(ii,jj) .ge.  -20.8 .and. lon(ii,jj) .lt.  -10.0 .and. lat(ii,jj) .ge.  -74.45 .and. lat(ii,jj) .lt.  -71.65 )  isfmask(ii,jj) = 43
    if ( lon(ii,jj) .ge.  -12.0 .and. lon(ii,jj) .lt.  -10.0 .and. lat(ii,jj) .ge.  -71.65 .and. lat(ii,jj) .lt.  -70.00 )  isfmask(ii,jj) = 28
    if ( lon(ii,jj) .ge.  -10.0 .and. lon(ii,jj) .lt.   -7.5 .and. lat(ii,jj) .ge.  -72.00 .and. lat(ii,jj) .lt.  -70.90 )  isfmask(ii,jj) = 12 ! Ekstrom
    if ( lon(ii,jj) .ge.  -10.0 .and. lon(ii,jj) .lt.   -7.8 .and. lat(ii,jj) .ge.  -70.90 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) = 12 ! Ekstrom
    if ( lon(ii,jj) .ge.   -7.8 .and. lon(ii,jj) .lt.   -6.1 .and. lat(ii,jj) .ge.  -70.90 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) = 57 ! Atka
    if ( lon(ii,jj) .ge.   -6.1 .and. lon(ii,jj) .lt.   -2.9 .and. lat(ii,jj) .ge.  -70.90 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) = 70 ! Jelbart
    if ( lon(ii,jj) .ge.   -6.3 .and. lon(ii,jj) .lt.   -2.7 .and. lat(ii,jj) .ge.  -72.50 .and. lat(ii,jj) .lt.  -70.90 )  isfmask(ii,jj) = 70 ! Jelbart
    if ( lon(ii,jj) .ge.   -2.9 .and. lon(ii,jj) .lt.   -2.6 .and. lat(ii,jj) .ge.  -70.90 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) = 44 ! Fimbul
    if ( lon(ii,jj) .ge.   -2.6 .and. lon(ii,jj) .lt.    7.6 .and. lat(ii,jj) .ge.  -72.50 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) = 44 ! Fimbul
    if ( lon(ii,jj) .ge.    7.6 .and. lon(ii,jj) .lt.    9.3 .and. lat(ii,jj) .ge.  -71.50 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) = 29
    if ( lon(ii,jj) .ge.    9.3 .and. lon(ii,jj) .lt.   12.7 .and. lat(ii,jj) .ge.  -70.10 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) = 13 ! Nivl
    if ( lon(ii,jj) .ge.    9.3 .and. lon(ii,jj) .lt.   13.05.and. lat(ii,jj) .ge.  -71.50 .and. lat(ii,jj) .lt.  -70.10 )  isfmask(ii,jj) = 13 ! Nivl
    if ( lon(ii,jj) .ge.   12.7 .and. lon(ii,jj) .lt.   16.15.and. lat(ii,jj) .ge.  -70.10 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) = 58 ! Lazarev
    if ( lon(ii,jj) .ge.   13.05.and. lon(ii,jj) .lt.   16.15.and. lat(ii,jj) .ge.  -71.50 .and. lat(ii,jj) .lt.  -70.10 )  isfmask(ii,jj) = 58 ! Lazarev
    if ( lon(ii,jj) .ge.   16.15.and. lon(ii,jj) .lt.   23.9 .and. lat(ii,jj) .ge.  -71.50 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) = 71
    if ( lon(ii,jj) .ge.   23.9 .and. lon(ii,jj) .lt.   33.5 .and. lat(ii,jj) .ge.  -71.50 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) = 45
    if ( lon(ii,jj) .ge.   33.5 .and. lon(ii,jj) .lt.   38.02.and. lat(ii,jj) .ge.  -70.50 .and. lat(ii,jj) .lt.  -68.60 )  isfmask(ii,jj) = 30
    if ( lon(ii,jj) .ge.   38.02.and. lon(ii,jj) .lt.   39.5 .and. lat(ii,jj) .ge.  -71.00 .and. lat(ii,jj) .lt.  -69.60 )  isfmask(ii,jj) = 14
    if ( lon(ii,jj) .ge.   47.3 .and. lon(ii,jj) .lt.   47.85.and. lat(ii,jj) .ge.  -67.80 .and. lat(ii,jj) .lt.  -67.60 )  isfmask(ii,jj) = 59 ! Rayner/Thyer
    if ( lon(ii,jj) .ge.   47.85.and. lon(ii,jj) .lt.   49.0 .and. lat(ii,jj) .ge.  -68.10 .and. lat(ii,jj) .lt.  -67.50 )  isfmask(ii,jj) = 59 ! Rayner/Thyer
    if ( lon(ii,jj) .ge.   55.5 .and. lon(ii,jj) .lt.   57.5 .and. lat(ii,jj) .ge.  -66.87 .and. lat(ii,jj) .lt.  -66.50 )  isfmask(ii,jj) = 72
    if ( lon(ii,jj) .ge.   55.5 .and. lon(ii,jj) .lt.   57.5 .and. lat(ii,jj) .ge.  -67.60 .and. lat(ii,jj) .lt.  -66.87 )  isfmask(ii,jj) = 46 ! Wilma/Robert/Downer
    if ( lon(ii,jj) .ge.   66.0 .and. lon(ii,jj) .lt.   74.0 .and. lat(ii,jj) .ge.  -75.00 .and. lat(ii,jj) .lt.  -69.60 )  isfmask(ii,jj) = 31 ! Amery
    if ( lon(ii,jj) .ge.   68.0 .and. lon(ii,jj) .lt.   74.4 .and. lat(ii,jj) .ge.  -69.60 .and. lat(ii,jj) .lt.  -68.30 )  isfmask(ii,jj) = 31 ! Amery
    if ( lon(ii,jj) .ge.   74.0 .and. lon(ii,jj) .lt.   76.0 .and. lat(ii,jj) .ge.  -70.20 .and. lat(ii,jj) .lt.  -69.60 )  isfmask(ii,jj) = 15 ! Publications
    if ( lon(ii,jj) .ge.   74.3 .and. lon(ii,jj) .lt.   76.0 .and. lat(ii,jj) .ge.  -69.60 .and. lat(ii,jj) .lt.  -69.30 )  isfmask(ii,jj) = 15 ! Publications
    if ( lon(ii,jj) .ge.   80.5 .and. lon(ii,jj) .lt.   89.0 .and. lat(ii,jj) .ge.  -68.00 .and. lat(ii,jj) .lt.  -66.00 )  isfmask(ii,jj) = 61
    if ( lon(ii,jj) .ge.   94.5 .and. lon(ii,jj) .lt.  101.0 .and. lat(ii,jj) .ge.  -67.50 .and. lat(ii,jj) .lt.  -64.50 )  isfmask(ii,jj) = 73 ! Shackleton **NB** keep before Tracy Tremenchus
    if ( lon(ii,jj) .ge.  100.55.and. lon(ii,jj) .lt.  102.7 .and. lat(ii,jj) .ge.  -66.10 .and. lat(ii,jj) .lt.  -65.00 )  isfmask(ii,jj) = 47 ! Tracy Tremenchus
    if ( lon(ii,jj) .ge.  100.4 .and. lon(ii,jj) .lt.  100.55.and. lat(ii,jj) .ge.  -66.10 .and. lat(ii,jj) .lt.  -65.85 )  isfmask(ii,jj) = 47 ! Tracy Tremenchus
    if ( lon(ii,jj) .ge.  100.35.and. lon(ii,jj) .lt.  100.4 .and. lat(ii,jj) .ge.  -66.12 .and. lat(ii,jj) .lt.  -65.90 )  isfmask(ii,jj) = 47 ! Tracy Tremenchus
    if ( lon(ii,jj) .ge.  100.3 .and. lon(ii,jj) .lt.  100.35.and. lat(ii,jj) .ge.  -66.14 .and. lat(ii,jj) .lt.  -65.95 )  isfmask(ii,jj) = 47 ! Tracy Tremenchus
    if ( lon(ii,jj) .ge.  102.7 .and. lon(ii,jj) .lt.  104.5 .and. lat(ii,jj) .ge.  -66.50 .and. lat(ii,jj) .lt.  -65.30 )  isfmask(ii,jj) = 32
    if ( lon(ii,jj) .ge.  108.8 .and. lon(ii,jj) .lt.  112.5 .and. lat(ii,jj) .ge.  -67.70 .and. lat(ii,jj) .lt.  -66.20 )  isfmask(ii,jj) = 16 ! Vincennes Bay
    if ( lon(ii,jj) .ge.  113.0 .and. lon(ii,jj) .lt.  117.4 .and. lat(ii,jj) .ge.  -67.80 .and. lat(ii,jj) .lt.  -66.52 )  isfmask(ii,jj) = 48 ! Totten
    if ( lon(ii,jj) .ge.  117.4 .and. lon(ii,jj) .lt.  123.0 .and. lat(ii,jj) .ge.  -67.80 .and. lat(ii,jj) .lt.  -66.40 )  isfmask(ii,jj) = 33 ! Moscow Univ.
    if ( lon(ii,jj) .ge.  126.48.and. lon(ii,jj) .lt.  128.05.and. lat(ii,jj) .ge.  -67.20 .and. lat(ii,jj) .lt.  -66.45 )  isfmask(ii,jj) = 17 ! Holmes
    if ( lon(ii,jj) .ge.  133.7 .and. lon(ii,jj) .lt.  135.4 .and. lat(ii,jj) .ge.  -66.70 .and. lat(ii,jj) .lt.  -65.80 )  isfmask(ii,jj) = 62 ! Dibble
    if ( lon(ii,jj) .ge.  144.2 .and. lon(ii,jj) .lt.  147.0 .and. lat(ii,jj) .ge.  -67.48 .and. lat(ii,jj) .lt.  -66.50 )  isfmask(ii,jj) = 49 ! Mertz
    if ( lon(ii,jj) .ge.  143.5 .and. lon(ii,jj) .lt.  145.4 .and. lat(ii,jj) .ge.  -68.00 .and. lat(ii,jj) .lt.  -67.35 )  isfmask(ii,jj) = 49 ! Mertz
    if ( lon(ii,jj) .ge.  146.75.and. lon(ii,jj) .lt.  149.0 .and. lat(ii,jj) .ge.  -68.60 .and. lat(ii,jj) .lt.  -67.90 )  isfmask(ii,jj) = 34 ! Ninnis
    if ( lon(ii,jj) .ge.  146.2 .and. lon(ii,jj) .lt.  146.75.and. lat(ii,jj) .ge.  -68.40 .and. lat(ii,jj) .lt.  -67.98 )  isfmask(ii,jj) = 34 ! Ninnis
    if ( lon(ii,jj) .ge.  150.9 .and. lon(ii,jj) .lt.  154.1 .and. lat(ii,jj) .ge.  -69.50 .and. lat(ii,jj) .lt.  -67.50 )  isfmask(ii,jj) = 18 ! Cook
    if ( lon(ii,jj) .ge.  160.8 .and. lon(ii,jj) .lt.  162.03.and. lat(ii,jj) .ge.  -71.80 .and. lat(ii,jj) .lt.  -69.90 )  isfmask(ii,jj) = 63 ! Rennick
    if ( lon(ii,jj) .ge.  162.03.and. lon(ii,jj) .lt.  162.6 .and. lat(ii,jj) .ge.  -71.80 .and. lat(ii,jj) .lt.  -70.45 )  isfmask(ii,jj) = 63 ! Rennick
    if ( lon(ii,jj) .ge.  163.5 .and. lon(ii,jj) .lt.  164.6 .and. lat(ii,jj) .ge.  -71.25 .and. lat(ii,jj) .lt.  -70.25 )  isfmask(ii,jj) = 74 ! Lillie
    if ( lon(ii,jj) .ge.  167.4 .and. lon(ii,jj) .lt.  169.2 .and. lat(ii,jj) .ge.  -73.65 .and. lat(ii,jj) .lt.  -72.90 )  isfmask(ii,jj) = 50 ! Mariner
    if ( lon(ii,jj) .ge.  167.23.and. lon(ii,jj) .lt.  167.4 .and. lat(ii,jj) .ge.  -73.52 .and. lat(ii,jj) .lt.  -72.90 )  isfmask(ii,jj) = 50 ! Mariner
    if ( lon(ii,jj) .ge.  166.5 .and. lon(ii,jj) .lt.  167.23.and. lat(ii,jj) .ge.  -73.45 .and. lat(ii,jj) .lt.  -72.90 )  isfmask(ii,jj) = 50 ! Mariner
    if ( lon(ii,jj) .ge.  164.6 .and. lon(ii,jj) .lt.  165.5 .and. lat(ii,jj) .ge.  -73.98 .and. lat(ii,jj) .lt.  -73.50 )  isfmask(ii,jj) = 35 ! Aviator
    if ( lon(ii,jj) .ge.  165.15.and. lon(ii,jj) .lt.  165.8 .and. lat(ii,jj) .ge.  -74.06 .and. lat(ii,jj) .lt.  -73.88 )  isfmask(ii,jj) = 35 ! Aviator
    if ( lon(ii,jj) .ge.  165.8 .and. lon(ii,jj) .lt.  166.5 .and. lat(ii,jj) .ge.  -74.20 .and. lat(ii,jj) .lt.  -73.93 )  isfmask(ii,jj) = 35 ! Aviator
    if ( lon(ii,jj) .ge.  161.0 .and. lon(ii,jj) .lt.  164.1 .and. lat(ii,jj) .ge.  -75.20 .and. lat(ii,jj) .lt.  -74.30 )  isfmask(ii,jj) = 19 ! Nansen
    if ( lon(ii,jj) .ge.  155.0 .and. lon(ii,jj) .lt.  162.45.and. lat(ii,jj) .ge.  -75.50 .and. lat(ii,jj) .lt.  -75.20 )  isfmask(ii,jj) = 64 ! Drygalsky
    if ( lon(ii,jj) .ge.  162.45.and. lon(ii,jj) .lt.  166.5 .and. lat(ii,jj) .ge.  -75.45 .and. lat(ii,jj) .lt.  -75.20 )  isfmask(ii,jj) = 64 ! Drygalsky
    if ( lon(ii,jj) .ge.  163.1 .and. lon(ii,jj) .lt.  166.5 .and. lat(ii,jj) .ge.  -75.70 .and. lat(ii,jj) .lt.  -75.20 )  isfmask(ii,jj) = 64 ! Drygalsky
    if ( lon(ii,jj) .ge.  150.0 .and. lon(ii,jj) .le.  180.0 .and. lat(ii,jj) .ge.  -86.00 .and. lat(ii,jj) .lt.  -77.00 )  isfmask(ii,jj) = 10 ! Ross
    if ( lon(ii,jj) .ge. -180.0 .and. lon(ii,jj) .lt. -145.0 .and. lat(ii,jj) .ge.  -86.00 .and. lat(ii,jj) .lt.  -77.70 )  isfmask(ii,jj) = 10 ! Ross
    if ( lon(ii,jj) .ge. -158.0 .and. lon(ii,jj) .lt. -156.0 .and. lat(ii,jj) .ge.  -77.40 .and. lat(ii,jj) .lt.  -76.70 )  isfmask(ii,jj) = 36 ! Withrow
    if ( lon(ii,jj) .ge. -155.0 .and. lon(ii,jj) .lt. -152.35.and. lat(ii,jj) .ge.  -77.60 .and. lat(ii,jj) .lt.  -76.40 )  isfmask(ii,jj) = 20
    if ( lon(ii,jj) .ge. -152.35.and. lon(ii,jj) .lt. -146.5 .and. lat(ii,jj) .ge.  -77.80 .and. lat(ii,jj) .lt.  -76.40 )  isfmask(ii,jj) = 65
    if ( lon(ii,jj) .ge. -146.6 .and. lon(ii,jj) .lt. -142.5 .and. lat(ii,jj) .ge.  -77.80 .and. lat(ii,jj) .lt.  -76.50 )  isfmask(ii,jj) = 65
    if ( lon(ii,jj) .ge. -149.5 .and. lon(ii,jj) .lt. -146.6 .and. lat(ii,jj) .ge.  -76.40 .and. lat(ii,jj) .lt.  -75.00 )  isfmask(ii,jj) = 51
    if ( lon(ii,jj) .ge. -146.6 .and. lon(ii,jj) .lt. -142.3 .and. lat(ii,jj) .ge.  -76.50 .and. lat(ii,jj) .lt.  -75.00 )  isfmask(ii,jj) = 51
    if ( lon(ii,jj) .ge. -142.3 .and. lon(ii,jj) .lt. -140.0 .and. lat(ii,jj) .ge.  -76.00 .and. lat(ii,jj) .lt.  -75.00 )  isfmask(ii,jj) = 37
    if ( lon(ii,jj) .ge. -135.5 .and. lon(ii,jj) .lt. -118.0 .and. lat(ii,jj) .ge.  -76.00 .and. lat(ii,jj) .lt.  -74.60 )  isfmask(ii,jj) = 22 ! Getz
    if ( lon(ii,jj) .ge. -135.5 .and. lon(ii,jj) .lt. -126.1 .and. lat(ii,jj) .ge.  -74.60 .and. lat(ii,jj) .lt.  -72.50 )  isfmask(ii,jj) = 22
    if ( lon(ii,jj) .ge. -124.3 .and. lon(ii,jj) .lt. -114.5 .and. lat(ii,jj) .ge.  -74.60 .and. lat(ii,jj) .lt.  -72.50 )  isfmask(ii,jj) = 22
    if ( lon(ii,jj) .ge. -126.1 .and. lon(ii,jj) .lt. -124.3 .and. lat(ii,jj) .ge.  -74.60 .and. lat(ii,jj) .lt.  -73.90 )  isfmask(ii,jj) = 22
    if ( lon(ii,jj) .ge. -114.5 .and. lon(ii,jj) .lt. -111.2 .and. lat(ii,jj) .ge.  -74.90 .and. lat(ii,jj) .lt.  -74.10 )  isfmask(ii,jj) = 52 ! Dotson
    if ( lon(ii,jj) .ge. -115.5 .and. lon(ii,jj) .lt. -113.0 .and. lat(ii,jj) .ge.  -75.50 .and. lat(ii,jj) .lt.  -74.90 )  isfmask(ii,jj) = 52
    if ( lon(ii,jj) .ge. -111.2 .and. lon(ii,jj) .lt. -108.6 .and. lat(ii,jj) .ge.  -76.00 .and. lat(ii,jj) .lt.  -74.50 )  isfmask(ii,jj) = 38 ! Crosson
    if ( lon(ii,jj) .ge. -113.0 .and. lon(ii,jj) .lt. -112.6 .and. lat(ii,jj) .ge.  -75.50 .and. lat(ii,jj) .lt.  -74.90 )  isfmask(ii,jj) = 38
    if ( lon(ii,jj) .ge. -112.6 .and. lon(ii,jj) .lt. -111.2 .and. lat(ii,jj) .ge.  -75.50 .and. lat(ii,jj) .lt.  -74.87 )  isfmask(ii,jj) = 38
    if ( lon(ii,jj) .ge. -108.6 .and. lon(ii,jj) .lt. -103.5 .and. lat(ii,jj) .ge.  -77.00 .and. lat(ii,jj) .lt.  -74.10 )  isfmask(ii,jj) = 23
    if ( lon(ii,jj) .ge. -103.5 .and. lon(ii,jj) .lt.  -95.0 .and. lat(ii,jj) .ge.  -77.00 .and. lat(ii,jj) .lt.  -74.10 )  isfmask(ii,jj) = 66 ! PIG
    if ( lon(ii,jj) .ge. -101.6 .and. lon(ii,jj) .lt.  -98.0 .and. lat(ii,jj) .ge.  -73.90 .and. lat(ii,jj) .lt.  -73.60 )  isfmask(ii,jj) = 53 ! Cosgrove
    if ( lon(ii,jj) .ge. -105.0 .and. lon(ii,jj) .lt.  -98.0 .and. lat(ii,jj) .ge.  -73.60 .and. lat(ii,jj) .lt.  -73.10 )  isfmask(ii,jj) = 53
    if ( lon(ii,jj) .ge. -105.0 .and. lon(ii,jj) .lt.  -89.0 .and. lat(ii,jj) .ge.  -73.20 .and. lat(ii,jj) .lt.  -72.40 )  isfmask(ii,jj) = 39 ! Abbot
    if ( lon(ii,jj) .ge.  -98.5 .and. lon(ii,jj) .lt.  -89.0 .and. lat(ii,jj) .ge.  -73.70 .and. lat(ii,jj) .lt.  -73.20 )  isfmask(ii,jj) = 39 ! Abbot
    if ( lon(ii,jj) .ge. -105.0 .and. lon(ii,jj) .lt. -100.0 .and. lat(ii,jj) .ge.  -72.40 .and. lat(ii,jj) .lt.  -72.10 )  isfmask(ii,jj) = 39 ! Abbot
    if ( lon(ii,jj) .ge. -100.0 .and. lon(ii,jj) .lt.  -98.0 .and. lat(ii,jj) .ge.  -72.40 .and. lat(ii,jj) .lt.  -72.20 )  isfmask(ii,jj) = 39 ! Abbot
    if ( lon(ii,jj) .ge.  -96.5 .and. lon(ii,jj) .lt.  -94.0 .and. lat(ii,jj) .ge.  -72.40 .and. lat(ii,jj) .lt.  -72.21 )  isfmask(ii,jj) = 39 ! Abbot
    if ( lon(ii,jj) .ge.  -95.8 .and. lon(ii,jj) .lt.  -95.0 .and. lat(ii,jj) .ge.  -72.21 .and. lat(ii,jj) .lt.  -72.07 )  isfmask(ii,jj) = 39 ! Abbot
    if ( lon(ii,jj) .ge.  -89.0 .and. lon(ii,jj) .lt.  -85.7 .and. lat(ii,jj) .ge.  -73.80 .and. lat(ii,jj) .lt.  -72.40 )  isfmask(ii,jj) = 24 ! Venable
    if ( lon(ii,jj) .ge.  -85.0 .and. lon(ii,jj) .lt.  -80.0 .and. lat(ii,jj) .ge.  -74.30 .and. lat(ii,jj) .lt.  -72.50 )  isfmask(ii,jj) = 67 ! Ferrigno
    if ( lon(ii,jj) .ge.  -75.0 .and. lon(ii,jj) .lt.  -74.0 .and. lat(ii,jj) .ge.  -70.50 .and. lat(ii,jj) .lt.  -69.62 )  isfmask(ii,jj) = 75 ! Wilkins
    if ( lon(ii,jj) .ge.  -74.0 .and. lon(ii,jj) .lt.  -69.5 .and. lat(ii,jj) .ge.  -71.15 .and. lat(ii,jj) .lt.  -69.62 )  isfmask(ii,jj) = 75 ! Wilkins
    if ( lon(ii,jj) .ge.  -73.2 .and. lon(ii,jj) .lt.  -69.5 .and. lat(ii,jj) .ge.  -72.50 .and. lat(ii,jj) .lt.  -71.50 )  isfmask(ii,jj) = 54 ! Bach
    if ( lon(ii,jj) .ge.  -74.0 .and. lon(ii,jj) .lt.  -69.5 .and. lat(ii,jj) .ge.  -72.50 .and. lat(ii,jj) .lt.  -71.70 )  isfmask(ii,jj) = 54 ! Bach
    if ( lon(ii,jj) .ge.  -79.0 .and. lon(ii,jj) .lt.  -77.0 .and. lat(ii,jj) .ge.  -74.30 .and. lat(ii,jj) .lt.  -72.54 )  isfmask(ii,jj) = 40 ! Stange
    if ( lon(ii,jj) .ge.  -77.0 .and. lon(ii,jj) .lt.  -75.7 .and. lat(ii,jj) .ge.  -74.30 .and. lat(ii,jj) .lt.  -72.50 )  isfmask(ii,jj) = 40 ! Stange
    if ( lon(ii,jj) .ge.  -75.7 .and. lon(ii,jj) .lt.  -74.8 .and. lat(ii,jj) .ge.  -74.30 .and. lat(ii,jj) .lt.  -72.90 )  isfmask(ii,jj) = 40 ! Stange
    if ( lon(ii,jj) .ge.  -74.8 .and. lon(ii,jj) .lt.  -74.55.and. lat(ii,jj) .ge.  -74.30 .and. lat(ii,jj) .lt.  -73.50 )  isfmask(ii,jj) = 40 ! Stange
    if ( lon(ii,jj) .ge.  -74.55.and. lon(ii,jj) .lt.  -74.5 .and. lat(ii,jj) .ge.  -74.60 .and. lat(ii,jj) .lt.  -73.64 )  isfmask(ii,jj) = 40 ! Stange
    if ( lon(ii,jj) .ge.  -69.5 .and. lon(ii,jj) .lt.  -66.0 .and. lat(ii,jj) .ge.  -74.30 .and. lat(ii,jj) .lt.  -69.50 )  isfmask(ii,jj) = 25 ! Georges VI
    if ( lon(ii,jj) .ge.  -74.55.and. lon(ii,jj) .lt.  -72.4 .and. lat(ii,jj) .ge.  -74.30 .and. lat(ii,jj) .lt.  -72.65 )  isfmask(ii,jj) = 25 ! Georges VI
    if ( lon(ii,jj) .ge.  -72.4 .and. lon(ii,jj) .lt.  -66.0 .and. lat(ii,jj) .ge.  -74.30 .and. lat(ii,jj) .lt.  -72.50 )  isfmask(ii,jj) = 25 ! Georges VI
    if ( lon(ii,jj) .ge.  -75.0 .and. lon(ii,jj) .lt.  -74.55.and. lat(ii,jj) .ge.  -73.50 .and. lat(ii,jj) .lt.  -72.50 )  isfmask(ii,jj) = 25 ! Georges VI
    if ( lon(ii,jj) .ge.  -74.55.and. lon(ii,jj) .lt.  -74.5 .and. lat(ii,jj) .ge.  -73.64 .and. lat(ii,jj) .lt.  -73.50 )  isfmask(ii,jj) = 25 ! Georges VI
    if ( lon(ii,jj) .ge.  -69.0 .and. lon(ii,jj) .lt.  -66.0 .and. lat(ii,jj) .ge.  -69.80 .and. lat(ii,jj) .lt.  -68.50 )  isfmask(ii,jj) = 41 ! Wordie

    !---- small ice shelves (no name?) :
    ! Amundsen (keep 77:85 to spread melt rates):
    if ( lon(ii,jj) .ge. -140.0 .and. lon(ii,jj) .lt. -138.0 .and. lat(ii,jj) .ge.  -75.20 .and. lat(ii,jj) .lt.  -75.10 )  isfmask(ii,jj) = 77
    if ( lon(ii,jj) .ge. -138.0 .and. lon(ii,jj) .lt. -136.0 .and. lat(ii,jj) .ge.  -75.20 .and. lat(ii,jj) .lt.  -74.97 )  isfmask(ii,jj) = 78
    if ( lon(ii,jj) .ge. -138.0 .and. lon(ii,jj) .lt. -136.6 .and. lat(ii,jj) .ge.  -74.97 .and. lat(ii,jj) .lt.  -74.72 )  isfmask(ii,jj) = 79
    if ( lon(ii,jj) .ge. -136.6 .and. lon(ii,jj) .lt. -135.3 .and. lat(ii,jj) .ge.  -74.80 .and. lat(ii,jj) .lt.  -74.60 )  isfmask(ii,jj) = 80
    if ( lon(ii,jj) .ge. -126.1 .and. lon(ii,jj) .lt. -124.5 .and. lat(ii,jj) .ge.  -73.90 .and. lat(ii,jj) .lt.  -73.60 )  isfmask(ii,jj) = 81
    if ( lon(ii,jj) .ge. -114.7 .and. lon(ii,jj) .lt. -113.5 .and. lat(ii,jj) .ge.  -74.10 .and. lat(ii,jj) .lt.  -73.93 )  isfmask(ii,jj) = 82
    if ( lon(ii,jj) .ge. -111.2 .and. lon(ii,jj) .lt. -110.5 .and. lat(ii,jj) .ge.  -74.35 .and. lat(ii,jj) .lt.  -74.25 )  isfmask(ii,jj) = 83
    if ( lon(ii,jj) .ge. -110.8 .and. lon(ii,jj) .lt. -110.0 .and. lat(ii,jj) .ge.  -74.45 .and. lat(ii,jj) .lt.  -74.35 )  isfmask(ii,jj) = 84
    if ( lon(ii,jj) .ge. -102.3 .and. lon(ii,jj) .lt. -101.0 .and. lat(ii,jj) .ge.  -74.05 .and. lat(ii,jj) .lt.  -73.85 )  isfmask(ii,jj) = 85
    ! Bellingshausen / West Peninsula
    if ( lon(ii,jj) .ge.  -96.8 .and. lon(ii,jj) .lt.  -95.0 .and. lat(ii,jj) .ge.  -72.07 .and. lat(ii,jj) .lt.  -71.80 )  isfmask(ii,jj) = 86
    if ( lon(ii,jj) .ge.  -97.1 .and. lon(ii,jj) .lt.  -95.0 .and. lat(ii,jj) .ge.  -72.21 .and. lat(ii,jj) .lt.  -72.07 )  isfmask(ii,jj) = 86
    if ( lon(ii,jj) .ge.  -97.1 .and. lon(ii,jj) .lt.  -96.5 .and. lat(ii,jj) .ge.  -72.30 .and. lat(ii,jj) .lt.  -72.21 )  isfmask(ii,jj) = 86
    if ( lon(ii,jj) .ge.  -97.7 .and. lon(ii,jj) .lt.  -97.1 .and. lat(ii,jj) .ge.  -72.30 .and. lat(ii,jj) .lt.  -72.07 )  isfmask(ii,jj) = 87
    if ( lon(ii,jj) .ge.  -97.7 .and. lon(ii,jj) .lt.  -96.9 .and. lat(ii,jj) .ge.  -72.07 .and. lat(ii,jj) .lt.  -71.80 )  isfmask(ii,jj) = 87
    if ( lon(ii,jj) .ge.  -98.5 .and. lon(ii,jj) .lt.  -97.7 .and. lat(ii,jj) .ge.  -72.20 .and. lat(ii,jj) .lt.  -71.80 )  isfmask(ii,jj) = 88
    if ( lon(ii,jj) .ge.  -99.5 .and. lon(ii,jj) .lt.  -98.5 .and. lat(ii,jj) .ge.  -72.20 .and. lat(ii,jj) .lt.  -71.80 )  isfmask(ii,jj) = 89
    if ( lon(ii,jj) .ge. -101.0 .and. lon(ii,jj) .lt.  -99.5 .and. lat(ii,jj) .ge.  -72.10 .and. lat(ii,jj) .lt.  -71.80 )  isfmask(ii,jj) = 89
    if ( lon(ii,jj) .ge.  -85.7 .and. lon(ii,jj) .lt.  -85.2 .and. lat(ii,jj) .ge.  -73.60 .and. lat(ii,jj) .lt.  -73.10 )  isfmask(ii,jj) = 90
    if ( lon(ii,jj) .ge.  -85.2 .and. lon(ii,jj) .lt.  -84.8 .and. lat(ii,jj) .ge.  -73.60 .and. lat(ii,jj) .lt.  -73.30 )  isfmask(ii,jj) = 90
    if ( lon(ii,jj) .ge.  -81.9 .and. lon(ii,jj) .lt.  -81.0 .and. lat(ii,jj) .ge.  -74.00 .and. lat(ii,jj) .lt.  -73.70 )  isfmask(ii,jj) = 91
    if ( lon(ii,jj) .ge.  -82.5 .and. lon(ii,jj) .lt.  -81.9 .and. lat(ii,jj) .ge.  -74.00 .and. lat(ii,jj) .lt.  -73.70 )  isfmask(ii,jj) = 92
    if ( lon(ii,jj) .ge.  -81.0 .and. lon(ii,jj) .lt.  -80.0 .and. lat(ii,jj) .ge.  -73.55 .and. lat(ii,jj) .lt.  -73.20 )  isfmask(ii,jj) = 93
    if ( lon(ii,jj) .ge.  -77.5 .and. lon(ii,jj) .lt.  -77.0 .and. lat(ii,jj) .ge.  -72.54 .and. lat(ii,jj) .lt.  -72.40 )  isfmask(ii,jj) = 94
    if ( lon(ii,jj) .ge.  -75.7 .and. lon(ii,jj) .lt.  -75.55.and. lat(ii,jj) .ge.  -72.90 .and. lat(ii,jj) .lt.  -72.70 )  isfmask(ii,jj) = 95
    if ( lon(ii,jj) .ge.  -75.55.and. lon(ii,jj) .lt.  -75.4 .and. lat(ii,jj) .ge.  -72.90 .and. lat(ii,jj) .lt.  -72.70 )  isfmask(ii,jj) = 96
    if ( lon(ii,jj) .ge.  -73.0 .and. lon(ii,jj) .lt.  -72.4 .and. lat(ii,jj) .ge.  -72.65 .and. lat(ii,jj) .lt.  -72.50 )  isfmask(ii,jj) = 97
    if ( lon(ii,jj) .ge.  -75.0 .and. lon(ii,jj) .lt.  -74.2 .and. lat(ii,jj) .ge.  -71.75 .and. lat(ii,jj) .lt.  -71.50 )  isfmask(ii,jj) = 98
    if ( lon(ii,jj) .ge.  -74.2 .and. lon(ii,jj) .lt.  -73.0 .and. lat(ii,jj) .ge.  -71.70 .and. lat(ii,jj) .lt.  -71.50 )  isfmask(ii,jj) = 99
    if ( lon(ii,jj) .ge.  -74.2 .and. lon(ii,jj) .lt.  -73.5 .and. lat(ii,jj) .ge.  -71.50 .and. lat(ii,jj) .lt.  -71.20 )  isfmask(ii,jj) = 99
    if ( lon(ii,jj) .ge.  -73.5 .and. lon(ii,jj) .lt.  -72.0 .and. lat(ii,jj) .ge.  -71.50 .and. lat(ii,jj) .lt.  -71.15 )  isfmask(ii,jj) =100
    if ( lon(ii,jj) .ge.  -74.8 .and. lon(ii,jj) .lt.  -74.4 .and. lat(ii,jj) .ge.  -70.90 .and. lat(ii,jj) .lt.  -70.50 )  isfmask(ii,jj) =101
    if ( lon(ii,jj) .ge.  -74.4 .and. lon(ii,jj) .lt.  -74.0 .and. lat(ii,jj) .ge.  -70.90 .and. lat(ii,jj) .lt.  -70.50 )  isfmask(ii,jj) =102
    if ( lon(ii,jj) .ge.  -72.5 .and. lon(ii,jj) .lt.  -72.0 .and. lat(ii,jj) .ge.  -69.62 .and. lat(ii,jj) .lt.  -69.50 )  isfmask(ii,jj) =103
    if ( lon(ii,jj) .ge.  -67.0 .and. lon(ii,jj) .lt.  -66.7 .and. lat(ii,jj) .ge.  -67.35 .and. lat(ii,jj) .lt.  -67.15 )  isfmask(ii,jj) =104
    if ( lon(ii,jj) .ge.  -67.2 .and. lon(ii,jj) .lt.  -66.7 .and. lat(ii,jj) .ge.  -67.60 .and. lat(ii,jj) .lt.  -67.40 )  isfmask(ii,jj) =105
    if ( lon(ii,jj) .ge.  -66.7 .and. lon(ii,jj) .lt.  -66.5 .and. lat(ii,jj) .ge.  -67.70 .and. lat(ii,jj) .lt.  -67.50 )  isfmask(ii,jj) =106
    ! Eastern Peninsula :
    if ( lon(ii,jj) .ge.  -61.1 .and. lon(ii,jj) .lt.  -59.5 .and. lat(ii,jj) .ge.  -65.30 .and. lat(ii,jj) .lt.  -64.70 )  isfmask(ii,jj) =107 ! Larsen A/B
    if ( lon(ii,jj) .ge.  -61.3 .and. lon(ii,jj) .lt.  -61.1 .and. lat(ii,jj) .ge.  -65.10 .and. lat(ii,jj) .lt.  -64.90 )  isfmask(ii,jj) =108
    if ( lon(ii,jj) .ge.  -61.8 .and. lon(ii,jj) .lt.  -61.3 .and. lat(ii,jj) .ge.  -65.15 .and. lat(ii,jj) .lt.  -64.95 )  isfmask(ii,jj) =109
    if ( lon(ii,jj) .ge.  -62.2 .and. lon(ii,jj) .lt.  -61.9 .and. lat(ii,jj) .ge.  -65.27 .and. lat(ii,jj) .lt.  -65.12 )  isfmask(ii,jj) =110
    if ( lon(ii,jj) .ge.  -62.3 .and. lon(ii,jj) .lt.  -62.0 .and. lat(ii,jj) .ge.  -65.37 .and. lat(ii,jj) .lt.  -65.27 )  isfmask(ii,jj) =111
    if ( lon(ii,jj) .ge.  -62.3 .and. lon(ii,jj) .lt.  -62.0 .and. lat(ii,jj) .ge.  -65.43 .and. lat(ii,jj) .lt.  -65.37 )  isfmask(ii,jj) =112
    if ( lon(ii,jj) .ge.  -62.3 .and. lon(ii,jj) .lt.  -62.0 .and. lat(ii,jj) .ge.  -65.49 .and. lat(ii,jj) .lt.  -65.43 )  isfmask(ii,jj) =113
    if ( lon(ii,jj) .ge.  -62.0 .and. lon(ii,jj) .lt.  -61.8 .and. lat(ii,jj) .ge.  -65.55 .and. lat(ii,jj) .lt.  -65.49 )  isfmask(ii,jj) =114
    if ( lon(ii,jj) .ge.  -61.5 .and. lon(ii,jj) .lt.  -60.0 .and. lat(ii,jj) .ge.  -73.80 .and. lat(ii,jj) .lt.  -73.60 )  isfmask(ii,jj) =115
    if ( lon(ii,jj) .ge.  -61.2 .and. lon(ii,jj) .lt.  -60.8 .and. lat(ii,jj) .ge.  -74.50 .and. lat(ii,jj) .lt.  -74.45 )  isfmask(ii,jj) =116
    ! Dronning Maud Land 
    if ( lon(ii,jj) .ge.  -36.1 .and. lon(ii,jj) .lt.  -35.4 .and. lat(ii,jj) .ge.  -78.20 .and. lat(ii,jj) .lt.  -78.00 )  isfmask(ii,jj) =117
    if ( lon(ii,jj) .ge.  -35.4 .and. lon(ii,jj) .lt.  -34.2 .and. lat(ii,jj) .ge.  -78.00 .and. lat(ii,jj) .lt.  -77.60 )  isfmask(ii,jj) =118
    if ( lon(ii,jj) .ge.  -28.3 .and. lon(ii,jj) .lt.  -27.5 .and. lat(ii,jj) .ge.  -76.50 .and. lat(ii,jj) .lt.  -76.00 )  isfmask(ii,jj) =119
    if ( lon(ii,jj) .ge.  -27.0 .and. lon(ii,jj) .lt.  -26.4 .and. lat(ii,jj) .ge.  -76.21 .and. lat(ii,jj) .lt.  -76.05 )  isfmask(ii,jj) =120
    ! East Ant.
    if ( lon(ii,jj) .ge.   39.4 .and. lon(ii,jj) .lt.   39.9 .and. lat(ii,jj) .ge.  -69.75 .and. lat(ii,jj) .lt.  -69.50 )  isfmask(ii,jj) =121
    if ( lon(ii,jj) .ge.   40.6 .and. lon(ii,jj) .lt.   41.0 .and. lat(ii,jj) .ge.  -69.00 .and. lat(ii,jj) .lt.  -68.50 )  isfmask(ii,jj) =122
    if ( lon(ii,jj) .ge.   42.1 .and. lon(ii,jj) .lt.   42.8 .and. lat(ii,jj) .ge.  -68.60 .and. lat(ii,jj) .lt.  -68.10 )  isfmask(ii,jj) =123
    if ( lon(ii,jj) .ge.   44.2 .and. lon(ii,jj) .lt.   44.9 .and. lat(ii,jj) .ge.  -68.10 .and. lat(ii,jj) .lt.  -67.70 )  isfmask(ii,jj) =124
    if ( lon(ii,jj) .ge.   44.9 .and. lon(ii,jj) .lt.   45.3 .and. lat(ii,jj) .ge.  -67.90 .and. lat(ii,jj) .lt.  -67.50 )  isfmask(ii,jj) =125
    if ( lon(ii,jj) .ge.   45.3 .and. lon(ii,jj) .lt.   46.0 .and. lat(ii,jj) .ge.  -68.00 .and. lat(ii,jj) .lt.  -67.60 )  isfmask(ii,jj) =126
    if ( lon(ii,jj) .ge.   46.0 .and. lon(ii,jj) .lt.   46.8 .and. lat(ii,jj) .ge.  -67.80 .and. lat(ii,jj) .lt.  -67.40 )  isfmask(ii,jj) =127
    if ( lon(ii,jj) .ge.   46.8 .and. lon(ii,jj) .lt.   47.85.and. lat(ii,jj) .ge.  -67.60 .and. lat(ii,jj) .lt.  -67.10 )  isfmask(ii,jj) =128
    if ( lon(ii,jj) .ge.   48.5 .and. lon(ii,jj) .lt.   49.6 .and. lat(ii,jj) .ge.  -67.40 .and. lat(ii,jj) .lt.  -67.10 )  isfmask(ii,jj) =129
    if ( lon(ii,jj) .ge.   49.6 .and. lon(ii,jj) .lt.   50.5 .and. lat(ii,jj) .ge.  -67.40 .and. lat(ii,jj) .lt.  -67.10 )  isfmask(ii,jj) =130
    if ( lon(ii,jj) .ge.   48.9 .and. lon(ii,jj) .lt.   49.7 .and. lat(ii,jj) .ge.  -67.10 .and. lat(ii,jj) .lt.  -66.90 )  isfmask(ii,jj) =131
    if ( lon(ii,jj) .ge.   56.0 .and. lon(ii,jj) .lt.   56.8 .and. lat(ii,jj) .ge.  -66.50 .and. lat(ii,jj) .lt.  -66.30 )  isfmask(ii,jj) =132
    if ( lon(ii,jj) .ge.   57.9 .and. lon(ii,jj) .lt.   58.4 .and. lat(ii,jj) .ge.  -67.40 .and. lat(ii,jj) .lt.  -66.90 )  isfmask(ii,jj) =133
    if ( lon(ii,jj) .ge.   58.5 .and. lon(ii,jj) .lt.   58.95.and. lat(ii,jj) .ge.  -67.40 .and. lat(ii,jj) .lt.  -67.20 )  isfmask(ii,jj) =134
    if ( lon(ii,jj) .ge.   58.95.and. lon(ii,jj) .lt.   60.0 .and. lat(ii,jj) .ge.  -67.70 .and. lat(ii,jj) .lt.  -67.20 )  isfmask(ii,jj) =135
    if ( lon(ii,jj) .ge.   60.7 .and. lon(ii,jj) .lt.   61.0 .and. lat(ii,jj) .ge.  -67.60 .and. lat(ii,jj) .lt.  -67.30 )  isfmask(ii,jj) =136
    if ( lon(ii,jj) .ge.   61.0 .and. lon(ii,jj) .lt.   61.6 .and. lat(ii,jj) .ge.  -67.70 .and. lat(ii,jj) .lt.  -67.30 )  isfmask(ii,jj) =137
    if ( lon(ii,jj) .ge.   76.3 .and. lon(ii,jj) .lt.   76.55.and. lat(ii,jj) .ge.  -69.60 .and. lat(ii,jj) .lt.  -69.30 )  isfmask(ii,jj) =138
    if ( lon(ii,jj) .ge.   76.55.and. lon(ii,jj) .lt.   76.86.and. lat(ii,jj) .ge.  -69.45 .and. lat(ii,jj) .lt.  -69.10 )  isfmask(ii,jj) =139
    if ( lon(ii,jj) .ge.   76.86.and. lon(ii,jj) .lt.   77.0 .and. lat(ii,jj) .ge.  -69.45 .and. lat(ii,jj) .lt.  -69.10 )  isfmask(ii,jj) =140
    if ( lon(ii,jj) .ge.   77.3 .and. lon(ii,jj) .lt.   77.9 .and. lat(ii,jj) .ge.  -69.35 .and. lat(ii,jj) .lt.  -69.03 )  isfmask(ii,jj) =141
    if ( lon(ii,jj) .ge.   77.8 .and. lon(ii,jj) .lt.   78.2 .and. lat(ii,jj) .ge.  -69.03 .and. lat(ii,jj) .lt.  -68.95 )  isfmask(ii,jj) =142
    if ( lon(ii,jj) .ge.   77.5 .and. lon(ii,jj) .lt.   78.8 .and. lat(ii,jj) .ge.  -68.85 .and. lat(ii,jj) .lt.  -68.50 )  isfmask(ii,jj) =143
    if ( lon(ii,jj) .ge.   89.8 .and. lon(ii,jj) .lt.   91.0 .and. lat(ii,jj) .ge.  -67.00 .and. lat(ii,jj) .lt.  -66.40 )  isfmask(ii,jj) =144
    if ( lon(ii,jj) .ge.   93.3 .and. lon(ii,jj) .lt.   94.5 .and. lat(ii,jj) .ge.  -67.00 .and. lat(ii,jj) .lt.  -66.40 )  isfmask(ii,jj) =145
    if ( lon(ii,jj) .ge.  106.2 .and. lon(ii,jj) .lt.  106.5 .and. lat(ii,jj) .ge.  -66.50 .and. lat(ii,jj) .lt.  -66.30 )  isfmask(ii,jj) =146
    if ( lon(ii,jj) .ge.  106.5 .and. lon(ii,jj) .lt.  106.9 .and. lat(ii,jj) .ge.  -66.55 .and. lat(ii,jj) .lt.  -66.35 )  isfmask(ii,jj) =147
    if ( lon(ii,jj) .ge.  107.4 .and. lon(ii,jj) .lt.  108.2 .and. lat(ii,jj) .ge.  -66.90 .and. lat(ii,jj) .lt.  -66.45 )  isfmask(ii,jj) =148
    if ( lon(ii,jj) .ge.  108.2 .and. lon(ii,jj) .lt.  108.5 .and. lat(ii,jj) .ge.  -66.80 .and. lat(ii,jj) .lt.  -66.65 )  isfmask(ii,jj) =149
    if ( lon(ii,jj) .ge.  123.5 .and. lon(ii,jj) .lt.  125.4 .and. lat(ii,jj) .ge.  -66.90 .and. lat(ii,jj) .lt.  -66.30 )  isfmask(ii,jj) =150
    if ( lon(ii,jj) .ge.  125.8 .and. lon(ii,jj) .lt.  126.48.and. lat(ii,jj) .ge.  -66.55 .and. lat(ii,jj) .lt.  -66.10 )  isfmask(ii,jj) =151
    if ( lon(ii,jj) .ge.  126.48.and. lon(ii,jj) .lt.  127.0 .and. lat(ii,jj) .ge.  -66.45 .and. lat(ii,jj) .lt.  -66.35 )  isfmask(ii,jj) =152
    if ( lon(ii,jj) .ge.  128.05.and. lon(ii,jj) .lt.  128.2 .and. lat(ii,jj) .ge.  -67.20 .and. lat(ii,jj) .lt.  -66.80 )  isfmask(ii,jj) =153
    if ( lon(ii,jj) .ge.  128.2 .and. lon(ii,jj) .lt.  128.8 .and. lat(ii,jj) .ge.  -67.20 .and. lat(ii,jj) .lt.  -66.80 )  isfmask(ii,jj) =154
    if ( lon(ii,jj) .ge.  128.8 .and. lon(ii,jj) .lt.  129.7 .and. lat(ii,jj) .ge.  -67.20 .and. lat(ii,jj) .lt.  -66.85 )  isfmask(ii,jj) =155
    if ( lon(ii,jj) .ge.  129.3 .and. lon(ii,jj) .lt.  129.9 .and. lat(ii,jj) .ge.  -66.85 .and. lat(ii,jj) .lt.  -66.55 )  isfmask(ii,jj) =156
    if ( lon(ii,jj) .ge.  129.7 .and. lon(ii,jj) .lt.  131.0 .and. lat(ii,jj) .ge.  -66.40 .and. lat(ii,jj) .lt.  -66.05 )  isfmask(ii,jj) =157
    !if ( lon(ii,jj) .ge.  133.7 .and. lon(ii,jj) .lt.  135.4 .and. lat(ii,jj) .ge.  -66.70 .and. lat(ii,jj) .lt.  -65.80 )  isfmask(ii,jj) =158 ! Dibble !!
    if ( lon(ii,jj) .ge.  135.4 .and. lon(ii,jj) .lt.  136.2 .and. lat(ii,jj) .ge.  -66.60 .and. lat(ii,jj) .lt.  -65.90 )  isfmask(ii,jj) =159
    if ( lon(ii,jj) .ge.  136.2 .and. lon(ii,jj) .lt.  137.0 .and. lat(ii,jj) .ge.  -66.60 .and. lat(ii,jj) .lt.  -66.25 )  isfmask(ii,jj) =160
    if ( lon(ii,jj) .ge.  137.5 .and. lon(ii,jj) .lt.  138.0 .and. lat(ii,jj) .ge.  -66.50 .and. lat(ii,jj) .lt.  -66.25 )  isfmask(ii,jj) =161
    if ( lon(ii,jj) .ge.  138.0 .and. lon(ii,jj) .lt.  138.4 .and. lat(ii,jj) .ge.  -66.70 .and. lat(ii,jj) .lt.  -66.40 )  isfmask(ii,jj) =162
    if ( lon(ii,jj) .ge.  138.4 .and. lon(ii,jj) .lt.  138.9 .and. lat(ii,jj) .ge.  -66.70 .and. lat(ii,jj) .lt.  -66.40 )  isfmask(ii,jj) =163
    if ( lon(ii,jj) .ge.  139.2 .and. lon(ii,jj) .lt.  139.7 .and. lat(ii,jj) .ge.  -66.70 .and. lat(ii,jj) .lt.  -66.40 )  isfmask(ii,jj) =164
    if ( lon(ii,jj) .ge.  139.7 .and. lon(ii,jj) .lt.  140.3 .and. lat(ii,jj) .ge.  -66.85 .and. lat(ii,jj) .lt.  -66.55 )  isfmask(ii,jj) =165
    if ( lon(ii,jj) .ge.  140.9 .and. lon(ii,jj) .lt.  141.6 .and. lat(ii,jj) .ge.  -66.95 .and. lat(ii,jj) .lt.  -66.65 )  isfmask(ii,jj) =166
    if ( lon(ii,jj) .ge.  143.5 .and. lon(ii,jj) .lt.  144.2 .and. lat(ii,jj) .ge.  -67.25 .and. lat(ii,jj) .lt.  -66.90 )  isfmask(ii,jj) =167
    if ( lon(ii,jj) .ge.  145.4 .and. lon(ii,jj) .lt.  146.2 .and. lat(ii,jj) .ge.  -67.65 .and. lat(ii,jj) .lt.  -67.48 )  isfmask(ii,jj) =168
    if ( lon(ii,jj) .ge.  145.4 .and. lon(ii,jj) .lt.  145.9 .and. lat(ii,jj) .ge.  -67.80 .and. lat(ii,jj) .lt.  -67.65 )  isfmask(ii,jj) =168
    if ( lon(ii,jj) .ge.  145.9 .and. lon(ii,jj) .lt.  146.45.and. lat(ii,jj) .ge.  -67.86 .and. lat(ii,jj) .lt.  -67.65 )  isfmask(ii,jj) =169
    if ( lon(ii,jj) .ge.  146.45.and. lon(ii,jj) .lt.  146.75.and. lat(ii,jj) .ge.  -67.98 .and. lat(ii,jj) .lt.  -67.80 )  isfmask(ii,jj) =170
    if ( lon(ii,jj) .ge.  146.15.and. lon(ii,jj) .lt.  146.45.and. lat(ii,jj) .ge.  -67.98 .and. lat(ii,jj) .lt.  -67.86 )  isfmask(ii,jj) =170
    if ( lon(ii,jj) .ge.  149.1 .and. lon(ii,jj) .lt.  150.1 .and. lat(ii,jj) .ge.  -68.60 .and. lat(ii,jj) .lt.  -68.20 )  isfmask(ii,jj) =171
    if ( lon(ii,jj) .ge.  150.1 .and. lon(ii,jj) .lt.  150.9 .and. lat(ii,jj) .ge.  -68.60 .and. lat(ii,jj) .lt.  -68.20 )  isfmask(ii,jj) =172
    if ( lon(ii,jj) .ge.  154.1 .and. lon(ii,jj) .lt.  155.3 .and. lat(ii,jj) .ge.  -68.90 .and. lat(ii,jj) .lt.  -68.50 )  isfmask(ii,jj) =173
    if ( lon(ii,jj) .ge.  155.3 .and. lon(ii,jj) .lt.  155.7 .and. lat(ii,jj) .ge.  -69.20 .and. lat(ii,jj) .lt.  -68.90 )  isfmask(ii,jj) =174
    if ( lon(ii,jj) .ge.  155.7 .and. lon(ii,jj) .lt.  156.65.and. lat(ii,jj) .ge.  -69.30 .and. lat(ii,jj) .lt.  -68.90 )  isfmask(ii,jj) =175
    if ( lon(ii,jj) .ge.  156.65.and. lon(ii,jj) .lt.  157.12.and. lat(ii,jj) .ge.  -69.30 .and. lat(ii,jj) .lt.  -68.90 )  isfmask(ii,jj) =176
    if ( lon(ii,jj) .ge.  157.12.and. lon(ii,jj) .lt.  157.7 .and. lat(ii,jj) .ge.  -69.45 .and. lat(ii,jj) .lt.  -69.00 )  isfmask(ii,jj) =177
    if ( lon(ii,jj) .ge.  158.2 .and. lon(ii,jj) .lt.  158.7 .and. lat(ii,jj) .ge.  -69.50 .and. lat(ii,jj) .lt.  -69.15 )  isfmask(ii,jj) =178
    if ( lon(ii,jj) .ge.  158.7 .and. lon(ii,jj) .lt.  159.3 .and. lat(ii,jj) .ge.  -69.60 .and. lat(ii,jj) .lt.  -69.25 )  isfmask(ii,jj) =179
    if ( lon(ii,jj) .ge.  159.3 .and. lon(ii,jj) .lt.  160.2 .and. lat(ii,jj) .ge.  -69.70 .and. lat(ii,jj) .lt.  -69.40 )  isfmask(ii,jj) =180
    if ( lon(ii,jj) .ge.  160.0 .and. lon(ii,jj) .lt.  160.8 .and. lat(ii,jj) .ge.  -70.20 .and. lat(ii,jj) .lt.  -69.70 )  isfmask(ii,jj) =181
    if ( lon(ii,jj) .ge.  162.03.and. lon(ii,jj) .lt.  162.4 .and. lat(ii,jj) .ge.  -70.45 .and. lat(ii,jj) .lt.  -70.15 )  isfmask(ii,jj) =182
    if ( lon(ii,jj) .ge.  162.6 .and. lon(ii,jj) .lt.  163.0 .and. lat(ii,jj) .ge.  -70.55 .and. lat(ii,jj) .lt.  -70.25 )  isfmask(ii,jj) =183
    if ( lon(ii,jj) .ge.  163.0 .and. lon(ii,jj) .lt.  163.5 .and. lat(ii,jj) .ge.  -70.80 .and. lat(ii,jj) .lt.  -70.50 )  isfmask(ii,jj) =184
    if ( lon(ii,jj) .ge.  165.7 .and. lon(ii,jj) .lt.  166.3 .and. lat(ii,jj) .ge.  -70.85 .and. lat(ii,jj) .lt.  -70.45 )  isfmask(ii,jj) =185
    if ( lon(ii,jj) .ge.  167.0 .and. lon(ii,jj) .lt.  168.0 .and. lat(ii,jj) .ge.  -71.07 .and. lat(ii,jj) .lt.  -70.85 )  isfmask(ii,jj) =186
    if ( lon(ii,jj) .ge.  167.6 .and. lon(ii,jj) .lt.  168.6 .and. lat(ii,jj) .ge.  -71.25 .and. lat(ii,jj) .lt.  -71.07 )  isfmask(ii,jj) =187
    if ( lon(ii,jj) .ge.  170.6 .and. lon(ii,jj) .lt.  171.0 .and. lat(ii,jj) .ge.  -71.82 .and. lat(ii,jj) .lt.  -71.60 )  isfmask(ii,jj) =188
    if ( lon(ii,jj) .ge.  170.17.and. lon(ii,jj) .lt.  170.6 .and. lat(ii,jj) .ge.  -72.03 .and. lat(ii,jj) .lt.  -71.78 )  isfmask(ii,jj) =189
    if ( lon(ii,jj) .ge.  170.0 .and. lon(ii,jj) .lt.  170.17.and. lat(ii,jj) .ge.  -72.08 .and. lat(ii,jj) .lt.  -71.98 )  isfmask(ii,jj) =190
    if ( lon(ii,jj) .ge.  169.5 .and. lon(ii,jj) .lt.  170.0 .and. lat(ii,jj) .ge.  -72.24 .and. lat(ii,jj) .lt.  -72.05 )  isfmask(ii,jj) =191
    if ( lon(ii,jj) .ge.  169.6 .and. lon(ii,jj) .lt.  170.35.and. lat(ii,jj) .ge.  -72.45 .and. lat(ii,jj) .lt.  -72.24 )  isfmask(ii,jj) =192
    if ( lon(ii,jj) .ge.  169.0 .and. lon(ii,jj) .lt.  170.35.and. lat(ii,jj) .ge.  -72.85 .and. lat(ii,jj) .lt.  -72.45 )  isfmask(ii,jj) =193
    if ( lon(ii,jj) .ge.  169.2 .and. lon(ii,jj) .lt.  169.5 .and. lat(ii,jj) .ge.  -73.30 .and. lat(ii,jj) .lt.  -72.93 )  isfmask(ii,jj) =194
    if ( lon(ii,jj) .ge.  167.23.and. lon(ii,jj) .lt.  167.4 .and. lat(ii,jj) .ge.  -73.65 .and. lat(ii,jj) .lt.  -73.52 )  isfmask(ii,jj) =195
    if ( lon(ii,jj) .ge.  167.05.and. lon(ii,jj) .lt.  167.23.and. lat(ii,jj) .ge.  -73.61 .and. lat(ii,jj) .lt.  -73.45 )  isfmask(ii,jj) =195
    if ( lon(ii,jj) .ge.  166.6 .and. lon(ii,jj) .lt.  167.05.and. lat(ii,jj) .ge.  -73.64 .and. lat(ii,jj) .lt.  -73.45 )  isfmask(ii,jj) =196
    if ( lon(ii,jj) .ge.  166.72.and. lon(ii,jj) .lt.  167.08.and. lat(ii,jj) .ge.  -73.68 .and. lat(ii,jj) .lt.  -73.64 )  isfmask(ii,jj) =196
    if ( lon(ii,jj) .ge.  166.0 .and. lon(ii,jj) .lt.  166.6 .and. lat(ii,jj) .ge.  -73.64 .and. lat(ii,jj) .lt.  -73.45 )  isfmask(ii,jj) =197
    if ( lon(ii,jj) .ge.  166.0 .and. lon(ii,jj) .lt.  166.72.and. lat(ii,jj) .ge.  -73.68 .and. lat(ii,jj) .lt.  -73.64 )  isfmask(ii,jj) =197
    if ( lon(ii,jj) .ge.  166.0 .and. lon(ii,jj) .lt.  167.0 .and. lat(ii,jj) .ge.  -73.72 .and. lat(ii,jj) .lt.  -73.68 )  isfmask(ii,jj) =197
    if ( lon(ii,jj) .ge.  166.25.and. lon(ii,jj) .lt.  167.0 .and. lat(ii,jj) .ge.  -73.90 .and. lat(ii,jj) .lt.  -73.72 )  isfmask(ii,jj) =197
    if ( lon(ii,jj) .ge.  166.1 .and. lon(ii,jj) .lt.  166.25.and. lat(ii,jj) .ge.  -73.82 .and. lat(ii,jj) .lt.  -73.72 )  isfmask(ii,jj) =198
    if ( lon(ii,jj) .ge.  165.8 .and. lon(ii,jj) .lt.  166.25.and. lat(ii,jj) .ge.  -73.93 .and. lat(ii,jj) .lt.  -73.82 )  isfmask(ii,jj) =199
    if ( lon(ii,jj) .ge.  165.5 .and. lon(ii,jj) .lt.  165.8 .and. lat(ii,jj) .ge.  -73.88 .and. lat(ii,jj) .lt.  -73.78 )  isfmask(ii,jj) =199
    if ( lon(ii,jj) .ge.  164.6 .and. lon(ii,jj) .lt.  165.15.and. lat(ii,jj) .ge.  -74.25 .and. lat(ii,jj) .lt.  -73.98 )  isfmask(ii,jj) =200
    if ( lon(ii,jj) .ge.  164.9 .and. lon(ii,jj) .lt.  165.4 .and. lat(ii,jj) .ge.  -74.65 .and. lat(ii,jj) .lt.  -74.45 )  isfmask(ii,jj) =201
    if ( lon(ii,jj) .ge.  164.1 .and. lon(ii,jj) .lt.  164.7 .and. lat(ii,jj) .ge.  -74.80 .and. lat(ii,jj) .lt.  -74.45 )  isfmask(ii,jj) =202
    if ( lon(ii,jj) .ge.  162.45.and. lon(ii,jj) .lt.  163.1 .and. lat(ii,jj) .ge.  -75.68 .and. lat(ii,jj) .lt.  -75.45 )  isfmask(ii,jj) =203
    if ( lon(ii,jj) .ge.  162.0 .and. lon(ii,jj) .lt.  162.45.and. lat(ii,jj) .ge.  -75.68 .and. lat(ii,jj) .lt.  -75.50 )  isfmask(ii,jj) =203
    if ( lon(ii,jj) .ge.  162.3 .and. lon(ii,jj) .lt.  163.1 .and. lat(ii,jj) .ge.  -75.80 .and. lat(ii,jj) .lt.  -75.68 )  isfmask(ii,jj) =204
    if ( lon(ii,jj) .ge.  162.3 .and. lon(ii,jj) .lt.  163.1 .and. lat(ii,jj) .ge.  -76.00 .and. lat(ii,jj) .lt.  -75.80 )  isfmask(ii,jj) =205
    if ( lon(ii,jj) .ge.  162.2 .and. lon(ii,jj) .lt.  162.6 .and. lat(ii,jj) .ge.  -76.11 .and. lat(ii,jj) .lt.  -76.00 )  isfmask(ii,jj) =206
    if ( lon(ii,jj) .ge.  161.8 .and. lon(ii,jj) .lt.  163.3 .and. lat(ii,jj) .ge.  -76.40 .and. lat(ii,jj) .lt.  -76.11 )  isfmask(ii,jj) =208
    if ( lon(ii,jj) .ge.  161.8 .and. lon(ii,jj) .lt.  163.3 .and. lat(ii,jj) .ge.  -76.80 .and. lat(ii,jj) .lt.  -76.40 )  isfmask(ii,jj) =209
    ! 
    if ( lon(ii,jj) .ge. -158.6 .and. lon(ii,jj) .lt. -157.4 .and. lat(ii,jj) .ge.  -77.65 .and. lat(ii,jj) .lt.  -77.40 )  isfmask(ii,jj) =210
    if ( lon(ii,jj) .ge. -156.0 .and. lon(ii,jj) .lt. -155.1 .and. lat(ii,jj) .ge.  -77.40 .and. lat(ii,jj) .lt.  -76.70 )  isfmask(ii,jj) =211
    if ( lon(ii,jj) .ge. -139.2 .and. lon(ii,jj) .lt. -138.4 .and. lat(ii,jj) .ge.  -75.30 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) =212
    if ( lon(ii,jj) .ge. -138.0 .and. lon(ii,jj) .lt. -137.04.and. lat(ii,jj) .ge.  -75.30 .and. lat(ii,jj) .lt.  -74.80 )  isfmask(ii,jj) =213
    if ( lon(ii,jj) .ge. -137.04.and. lon(ii,jj) .lt. -136.5 .and. lat(ii,jj) .ge.  -75.30 .and. lat(ii,jj) .lt.  -75.05 )  isfmask(ii,jj) =213
    !if ( lon(ii,jj) .ge. -137.04.and. lon(ii,jj) .lt. -136.7 .and. lat(ii,jj) .ge.  -75.05 .and. lat(ii,jj) .lt.  -74.70 )  isfmask(ii,jj) =214
    !if ( lon(ii,jj) .ge. -136.73.and. lon(ii,jj) .lt. -136.3 .and. lat(ii,jj) .ge.  -74.80 .and. lat(ii,jj) .lt.  -74.70 )  isfmask(ii,jj) =215
    !if ( lon(ii,jj) .ge. -136.3 .and. lon(ii,jj) .lt. -135.4 .and. lat(ii,jj) .ge.  -74.85 .and. lat(ii,jj) .lt.  -74.60 )  isfmask(ii,jj) =216

    ! remove a lake from the ice shelf mask :
    if ( lon(ii,jj) .ge.   90.0 .and. lon(ii,jj) .lt.  120.0 .and. lat(ii,jj) .ge.  -80.00 .and. lat(ii,jj) .lt.  -75.00 )  isfmask(ii,jj) = 0

    ! if no named ice shelves, define by sector :
    if ( isfmask(ii,jj) .eq. 2 ) then
      if ( lon(ii,jj) .ge.  -20.0 .and. lon(ii,jj) .lt.   45.0 ) isfmask(ii,jj) = 3
      if ( lon(ii,jj) .ge.   45.0 .and. lon(ii,jj) .lt.   90.0 ) isfmask(ii,jj) = 4
      if ( lon(ii,jj) .ge.   90.0 .and. lon(ii,jj) .lt.  165.0 ) isfmask(ii,jj) = 5
      if ( lon(ii,jj) .ge.  165.0 .and. lon(ii,jj) .le.  180.0 ) isfmask(ii,jj) = 6 !NB: do not exclude 180.00 !!
      if ( lon(ii,jj) .ge. -180.0 .and. lon(ii,jj) .lt. -140.0 ) isfmask(ii,jj) = 7
      if ( lon(ii,jj) .ge. -140.0 .and. lon(ii,jj) .lt.  -90.0 ) isfmask(ii,jj) = 8
      if ( lon(ii,jj) .ge.  -90.0 .and. lon(ii,jj) .lt.  -66.0 ) isfmask(ii,jj) = 9
      ! the remaining are in the Weddel Sea => isfmask=2
    endif


  elseif ( open_ocean_mask(ii,jj) .ge. 5.d-1 ) then

    isfmask(ii,jj) = 1  ! ocean with no ice

  else

    isfmask(ii,jj) = 0  ! grounded ice or rock emerging from the ice

  endif

enddo
enddo

Nisf = MAXVAL( isfmask )
ALLOCATE( name_isf(Nisf), name_reg(Nisf) )
ALLOCATE( melt_isf(Nisf), err_melt_isf(Nisf), area_isf(Nisf) )
name_isf(:)  = 'NOT USED           '
name_reg(:)  = 'NOT USED           '
melt_isf(:)     = 0.d0
err_melt_isf(:) = 0.d0
area_isf(:)     = 1.d-9
strlen=LEN(name_isf)

! Area of our ice shelves :
ALLOCATE( our_isf_area(Nisf) ) 
our_isf_area(:) = 0.d0
do ii=2,mx-1
do jj=2,my-1
  do kisf=1,Nisf
    if ( isfmask(ii,jj) .eq. kisf ) our_isf_area(kisf) = our_isf_area(kisf) + mesh_area 
  enddo
enddo
enddo

! NB: melt values below from Rignot 2013 unless specified otherwise

! ICE SHELF NAMES                      ! REGION NAMES                         ! ICE SHELF MELT (Gt/yr) ! MELT UNCERTAINTY           ! ICE-SHELF AREA IN RIGNOT (km^2)
name_isf(  2) = 'Weddell small isf  '; name_reg(  2) = 'Weddell            '; melt_isf(  2) =   0.d0 ; err_melt_isf(  2) =   0.d0 ; area_isf(  2) =      1.d-9
name_isf(  3) = 'DrMaudLd small isf '; name_reg(  3) = 'Droning Maud Land  '; melt_isf(  3) =   0.d0 ; err_melt_isf(  3) =   0.d0 ; area_isf(  3) =      1.d-9
name_isf(  4) = 'East Ant1 small isf'; name_reg(  4) = 'East 1             '; melt_isf(  4) =   0.d0 ; err_melt_isf(  4) =   0.d0 ; area_isf(  4) =      1.d-9
name_isf(  5) = 'East Ant2 small isf'; name_reg(  5) = 'East 2             '; melt_isf(  5) =   0.d0 ; err_melt_isf(  5) =   0.d0 ; area_isf(  5) =      1.d-9
name_isf(  6) = 'West Ross small isf'; name_reg(  6) = 'Western Ross       '; melt_isf(  6) =   0.d0 ; err_melt_isf(  6) =   0.d0 ; area_isf(  6) =      1.d-9
name_isf(  7) = 'East Ross small isf'; name_reg(  7) = 'Eastern Ross       '; melt_isf(  7) =   0.d0 ; err_melt_isf(  7) =   0.d0 ; area_isf(  7) =      1.d-9
name_isf(  8) = 'Amundsen small isf '; name_reg(  8) = 'Amundsen           '; melt_isf(  8) =   0.d0 ; err_melt_isf(  8) =   0.d0 ; area_isf(  8) =      1.d-9
name_isf(  9) = 'Bellingsh small isf'; name_reg(  9) = 'Bellingshausen     '; melt_isf(  9) =   0.d0 ; err_melt_isf(  9) =   0.d0 ; area_isf(  9) =      1.d-9
name_isf( 26) = 'Larsen B           '; name_reg( 26) = 'Weddell            '; melt_isf( 26) =  12.2  ; err_melt_isf( 26) =  14.0  ; area_isf( 26) =   6755.d0
name_isf( 42) = 'Larsen C           '; name_reg( 42) = 'Weddell            '; melt_isf( 42) =  20.7  ; err_melt_isf( 42) =  67.0  ; area_isf( 42) =  46465.d0
name_isf( 55) = 'Larsen D           '; name_reg( 55) = 'Weddell            '; melt_isf( 55) =   1.4  ; err_melt_isf( 55) =  14.0  ; area_isf( 55) =  22548.d0
name_isf( 68) = 'Larsen E           '; name_reg( 68) = 'Weddell            '; melt_isf( 68) =   1.4  ; err_melt_isf( 68) =   1.0  ; area_isf( 68) =   1184.d0
name_isf( 60) = 'Larsen F           '; name_reg( 60) = 'Weddell            '; melt_isf( 60) =   1.2  ; err_melt_isf( 60) =   0.4  ; area_isf( 60) =    828.d0
name_isf( 27) = 'Larsen G           '; name_reg( 27) = 'Weddell            '; melt_isf( 27) =   0.3  ; err_melt_isf( 27) =   0.2  ; area_isf( 27) =    412.d0
name_isf( 21) = 'Ronne              '; name_reg( 21) = 'Weddell            '; melt_isf( 21) = 113.5  ; err_melt_isf( 21) =  35.0  ; area_isf( 21) = 338887.d0
name_isf( 11) = 'Filchner           '; name_reg( 11) = 'Weddell            '; melt_isf( 11) =  41.9  ; err_melt_isf( 11) =  10.0  ; area_isf( 11) = 104253.d0
name_isf( 69) = 'Stancomb Brunt     '; name_reg( 69) = 'Dronning Maud Land '; melt_isf( 69) =   1.0  ; err_melt_isf( 69) =   7.0  ; area_isf( 69) =  36894.d0
name_isf( 43) = 'Riiser-Larsen      '; name_reg( 43) = 'Dronning Maud Land '; melt_isf( 43) =   8.7  ; err_melt_isf( 43) =   9.0  ; area_isf( 43) =  43450.d0
name_isf( 28) = 'Quar               '; name_reg( 28) = 'Dronning Maud Land '; melt_isf( 28) =   1.4  ; err_melt_isf( 28) =   0.5  ; area_isf( 28) =   2156.d0
name_isf( 12) = 'Ekstrom            '; name_reg( 12) = 'Dronning Maud Land '; melt_isf( 12) =   4.3  ; err_melt_isf( 12) =   2.0  ; area_isf( 12) =   6872.d0
name_isf( 57) = 'Atka               '; name_reg( 57) = 'Dronning Maud Land '; melt_isf( 57) =  -0.5  ; err_melt_isf( 57) =   0.4  ; area_isf( 57) =   1969.d0
name_isf( 70) = 'Jelbart            '; name_reg( 70) = 'Dronning Maud Land '; melt_isf( 70) =  -1.0  ; err_melt_isf( 70) =   3.0  ; area_isf( 70) =  10844.d0
name_isf( 44) = 'Fimbul             '; name_reg( 44) = 'Dronning Maud Land '; melt_isf( 44) =  23.5  ; err_melt_isf( 44) =   9.0  ; area_isf( 44) =  40843.d0
name_isf( 29) = 'Vigrid             '; name_reg( 29) = 'Dronning Maud Land '; melt_isf( 29) =   3.2  ; err_melt_isf( 29) =   0.7  ; area_isf( 29) =   2089.d0
name_isf( 13) = 'Nivl               '; name_reg( 13) = 'Dronning Maud Land '; melt_isf( 13) =   3.9  ; err_melt_isf( 13) =   2.0  ; area_isf( 13) =   7285.d0
name_isf( 58) = 'Lazarev            '; name_reg( 58) = 'Dronning Maud Land '; melt_isf( 58) =   6.3  ; err_melt_isf( 58) =   2.0  ; area_isf( 58) =   8519.d0
name_isf( 71) = 'Borchgrevink       '; name_reg( 71) = 'Dronning Maud Land '; melt_isf( 71) =   7.5  ; err_melt_isf( 71) =   6.0  ; area_isf( 71) =  21580.d0
name_isf( 45) = 'Roi Baudoin        '; name_reg( 45) = 'Dronning Maud Land '; melt_isf( 45) =  14.1  ; err_melt_isf( 45) =  12.0  ; area_isf( 45) =  32952.d0
name_isf( 30) = 'Prince Harald      '; name_reg( 30) = 'Dronning Maud Land '; melt_isf( 30) =  -2.0  ; err_melt_isf( 30) =   3.0  ; area_isf( 30) =   5392.d0
name_isf( 14) = 'Shirase            '; name_reg( 14) = 'Dronning Maud Land '; melt_isf( 14) =   5.7  ; err_melt_isf( 14) =   1.0  ; area_isf( 14) =    821.d0
name_isf( 59) = 'Rayner/Thyer       '; name_reg( 59) = 'East 1             '; melt_isf( 59) =   6.7  ; err_melt_isf( 59) =   1.0  ; area_isf( 59) =    641.d0
name_isf( 72) = 'Edward VIII        '; name_reg( 72) = 'East 1             '; melt_isf( 72) =   4.2  ; err_melt_isf( 72) =   0.8  ; area_isf( 72) =    411.d0
name_isf( 46) = 'Wilma/Robert/Downer'; name_reg( 46) = 'East 1             '; melt_isf( 46) =  10.d0 ; err_melt_isf( 46) =   0.6  ; area_isf( 46) =    858.d0
name_isf( 31) = 'Amery              '; name_reg( 31) = 'East 1             '; melt_isf( 31) =  35.5  ; err_melt_isf( 31) =  23.0  ; area_isf( 31) =  60654.d0
name_isf( 15) = 'Publications       '; name_reg( 15) = 'East 1             '; melt_isf( 15) =   1.5  ; err_melt_isf( 15) =   2.0  ; area_isf( 15) =   1551.d0
name_isf( 61) = 'West               '; name_reg( 61) = 'East 1             '; melt_isf( 61) =  27.2  ; err_melt_isf( 61) =  10.d0 ; area_isf( 61) =  15666.d0 
name_isf( 73) = 'Shackleton         '; name_reg( 73) = 'East 2             '; melt_isf( 73) =  72.6  ; err_melt_isf( 73) =  15.0  ; area_isf( 73) =  26080.d0
name_isf( 47) = 'Tracy Tremenchus   '; name_reg( 47) = 'East 2             '; melt_isf( 47) =   3.0  ; err_melt_isf( 47) =   2.0  ; area_isf( 47) =   2845.d0
name_isf( 32) = 'Conger/Glenzer     '; name_reg( 32) = 'East 2             '; melt_isf( 32) =   3.6  ; err_melt_isf( 32) =   1.0  ; area_isf( 32) =   1547.d0
name_isf( 16) = 'Vincennes          '; name_reg( 16) = 'East 2             '; melt_isf( 16) =   5.0  ; err_melt_isf( 16) =   2.0  ; area_isf( 16) =    935.d0
name_isf( 48) = 'Totten             '; name_reg( 48) = 'East 2             '; melt_isf( 48) =  63.2  ; err_melt_isf( 48) =   4.0  ; area_isf( 48) =   6032.d0
name_isf( 33) = 'Moscow Univ.       '; name_reg( 33) = 'East 2             '; melt_isf( 33) =  27.4  ; err_melt_isf( 33) =   4.0  ; area_isf( 33) =   5798.d0
name_isf( 17) = 'Holmes             '; name_reg( 17) = 'East 2             '; melt_isf( 17) =   6.7  ; err_melt_isf( 17) =   4.0  ; area_isf( 17) =   1921.d0 
name_isf( 62) = 'Dibble             '; name_reg( 62) = 'East 2             '; melt_isf( 62) =   8.1  ; err_melt_isf( 62) =   1.0  ; area_isf( 62) =   1482.d0
name_isf( 49) = 'Mertz              '; name_reg( 49) = 'East 2             '; melt_isf( 49) =   7.9  ; err_melt_isf( 49) =   3.0  ; area_isf( 49) =   5522.d0
name_isf( 34) = 'Ninnis             '; name_reg( 34) = 'East 2             '; melt_isf( 34) =   2.2  ; err_melt_isf( 34) =   3.0  ; area_isf( 34) =   1899.d0 
name_isf( 18) = 'Cook               '; name_reg( 18) = 'East 2             '; melt_isf( 18) =   4.6  ; err_melt_isf( 18) =   5.0  ; area_isf( 18) =   3462.d0
name_isf( 63) = 'Rennick            '; name_reg( 63) = 'East 2             '; melt_isf( 63) =   7.0  ; err_melt_isf( 63) =   1.0  ; area_isf( 63) =   3273.d0
name_isf( 74) = 'Lillie             '; name_reg( 74) = 'East 2             '; melt_isf( 74) =   3.4  ; err_melt_isf( 74) =   0.3  ; area_isf( 74) =    770.d0 
name_isf( 50) = 'Mariner            '; name_reg( 50) = 'Western Ross       '; melt_isf( 50) =   2.4  ; err_melt_isf( 50) =   0.6  ; area_isf( 50) =   2705.d0
name_isf( 35) = 'Aviator            '; name_reg( 35) = 'Western Ross       '; melt_isf( 35) =   1.4  ; err_melt_isf( 35) =   0.2  ; area_isf( 35) =    785.d0
name_isf( 19) = 'Nansen             '; name_reg( 19) = 'Western Ross       '; melt_isf( 19) =   1.1  ; err_melt_isf( 19) =   0.6  ; area_isf( 19) =   1985.d0
name_isf( 64) = 'Drygalski          '; name_reg( 64) = 'Western Ross       '; melt_isf( 64) =   7.6  ; err_melt_isf( 64) =   1.0  ; area_isf( 64) =   2338.d0
name_isf( 10) = 'Ross               '; name_reg( 10) = 'Ross               '; melt_isf( 10) =  47.7  ; err_melt_isf( 10) =  35.0  ; area_isf( 10) = 500809.d0  ! sum of absolute uncertainties for East & West
name_isf( 36) = 'Withrow            '; name_reg( 36) = 'Eastern Ross       '; melt_isf( 36) =   0.3  ; err_melt_isf( 36) =   0.4  ; area_isf( 36) =    632.d0
name_isf( 20) = 'Swinburne          '; name_reg( 20) = 'Eastern Ross       '; melt_isf( 20) =   3.8  ; err_melt_isf( 20) =   0.5  ; area_isf( 20) =    900.d0
name_isf( 65) = 'Sulzberger         '; name_reg( 65) = 'Eastern Ross       '; melt_isf( 65) =  18.2  ; err_melt_isf( 65) =   3.0  ; area_isf( 65) =  12333.d0 
name_isf( 51) = 'Nickerson          '; name_reg( 51) = 'Eastern Ross       '; melt_isf( 51) =   4.2  ; err_melt_isf( 51) =   2.0  ; area_isf( 51) =   6495.d0 
name_isf( 37) = 'Land               '; name_reg( 37) = 'Eastern Ross       '; melt_isf( 37) =   3.8  ; err_melt_isf( 37) =   1.0  ; area_isf( 37) =    640.d0
name_isf( 22) = 'Getz               '; name_reg( 22) = 'Amundsen           '; melt_isf( 22) = 144.9  ; err_melt_isf( 22) =  14.0  ; area_isf( 22) =  34018.d0
name_isf( 52) = 'Dotson             '; name_reg( 52) = 'Amundsen           '; melt_isf( 52) =  45.2  ; err_melt_isf( 52) =   4.0  ; area_isf( 52) =   5803.d0
name_isf( 38) = 'Crosson            '; name_reg( 38) = 'Amundsen           '; melt_isf( 38) =  38.5  ; err_melt_isf( 38) =   4.0  ; area_isf( 38) =   3229.d0
name_isf( 23) = 'Thwaites           '; name_reg( 23) = 'Amundsen           '; melt_isf( 23) =  97.5  ; err_melt_isf( 23) =   7.0  ; area_isf( 23) =   5499.d0
name_isf( 66) = 'Pine Island        '; name_reg( 66) = 'Amundsen           '; melt_isf( 66) = 101.2  ; err_melt_isf( 66) =   8.0  ; area_isf( 66) =   6249.d0
name_isf( 53) = 'Cosgrove           '; name_reg( 53) = 'Amundsen           '; melt_isf( 53) =   8.5  ; err_melt_isf( 53) =   2.0  ; area_isf( 53) =   3033.d0
name_isf( 39) = 'Abbot              '; name_reg( 39) = 'Bellingshausen     '; melt_isf( 39) =  51.8  ; err_melt_isf( 39) =  19.0  ; area_isf( 39) =  29688.d0
name_isf( 24) = 'Venable            '; name_reg( 24) = 'Bellingshausen     '; melt_isf( 24) =  19.4  ; err_melt_isf( 24) =   2.0  ; area_isf( 24) =   3194.d0
name_isf( 67) = 'Ferrigno           '; name_reg( 67) = 'Bellingshausen     '; melt_isf( 67) =   5.1  ; err_melt_isf( 67) =   2.0  ; area_isf( 67) =    117.d0
name_isf( 75) = 'Wilkins            '; name_reg( 75) = 'Bellingshausen     '; melt_isf( 75) =  18.4  ; err_melt_isf( 75) =  17.0  ; area_isf( 75) =  12866.d0
name_isf( 54) = 'Bach               '; name_reg( 54) = 'Bellingshausen     '; melt_isf( 54) =  10.4  ; err_melt_isf( 54) =   1.0  ; area_isf( 54) =   4579.d0
name_isf( 40) = 'Stange             '; name_reg( 40) = 'Bellingshausen     '; melt_isf( 40) =  28.0  ; err_melt_isf( 40) =   6.0  ; area_isf( 40) =   8027.d0
name_isf( 25) = 'George VI          '; name_reg( 25) = 'Bellingshausen     '; melt_isf( 25) =  89.0  ; err_melt_isf( 25) =  17.0  ; area_isf( 25) =  23434.d0
name_isf( 41) = 'Wordie             '; name_reg( 41) = 'Bellingshausen     '; melt_isf( 41) =   6.5  ; err_melt_isf( 41) =   3.0  ; area_isf( 41) =    277.d0
!
name_isf( 77) = '                   '; name_reg( 77) = '                   '; melt_isf( 77) =   0.d0 ; err_melt_isf( 77) =   0.d0 ; area_isf( 77) = 1.d-9
name_isf( 78) = '                   '; name_reg( 78) = '                   '; melt_isf( 78) =   0.d0 ; err_melt_isf( 78) =   0.d0 ; area_isf( 78) = 1.d-9
name_isf( 79) = '                   '; name_reg( 79) = '                   '; melt_isf( 79) =   0.d0 ; err_melt_isf( 79) =   0.d0 ; area_isf( 79) = 1.d-9
name_isf( 80) = '                   '; name_reg( 80) = '                   '; melt_isf( 80) =   0.d0 ; err_melt_isf( 80) =   0.d0 ; area_isf( 80) = 1.d-9
name_isf( 81) = '                   '; name_reg( 81) = '                   '; melt_isf( 81) =   0.d0 ; err_melt_isf( 81) =   0.d0 ; area_isf( 81) = 1.d-9
name_isf( 82) = '                   '; name_reg( 82) = '                   '; melt_isf( 82) =   0.d0 ; err_melt_isf( 82) =   0.d0 ; area_isf( 82) = 1.d-9
name_isf( 83) = '                   '; name_reg( 83) = '                   '; melt_isf( 83) =   0.d0 ; err_melt_isf( 83) =   0.d0 ; area_isf( 83) = 1.d-9
name_isf( 84) = '                   '; name_reg( 84) = '                   '; melt_isf( 84) =   0.d0 ; err_melt_isf( 84) =   0.d0 ; area_isf( 84) = 1.d-9
name_isf( 85) = '                   '; name_reg( 85) = '                   '; melt_isf( 85) =   0.d0 ; err_melt_isf( 85) =   0.d0 ; area_isf( 85) = 1.d-9
name_isf( 86) = '                   '; name_reg( 86) = '                   '; melt_isf( 86) =   0.d0 ; err_melt_isf( 86) =   0.d0 ; area_isf( 86) = 1.d-9
name_isf( 87) = '                   '; name_reg( 87) = '                   '; melt_isf( 87) =   0.d0 ; err_melt_isf( 87) =   0.d0 ; area_isf( 87) = 1.d-9
name_isf( 88) = '                   '; name_reg( 88) = '                   '; melt_isf( 88) =   0.d0 ; err_melt_isf( 88) =   0.d0 ; area_isf( 88) = 1.d-9
name_isf( 89) = '                   '; name_reg( 89) = '                   '; melt_isf( 89) =   0.d0 ; err_melt_isf( 89) =   0.d0 ; area_isf( 89) = 1.d-9
name_isf( 90) = '                   '; name_reg( 90) = '                   '; melt_isf( 90) =   0.d0 ; err_melt_isf( 90) =   0.d0 ; area_isf( 90) = 1.d-9
name_isf( 91) = '                   '; name_reg( 91) = '                   '; melt_isf( 91) =   0.d0 ; err_melt_isf( 91) =   0.d0 ; area_isf( 91) = 1.d-9
name_isf( 92) = '                   '; name_reg( 92) = '                   '; melt_isf( 92) =   0.d0 ; err_melt_isf( 92) =   0.d0 ; area_isf( 92) = 1.d-9
name_isf( 93) = '                   '; name_reg( 93) = '                   '; melt_isf( 93) =   0.d0 ; err_melt_isf( 93) =   0.d0 ; area_isf( 93) = 1.d-9
name_isf( 94) = '                   '; name_reg( 94) = '                   '; melt_isf( 94) =   0.d0 ; err_melt_isf( 94) =   0.d0 ; area_isf( 94) = 1.d-9
name_isf( 95) = '                   '; name_reg( 95) = '                   '; melt_isf( 95) =   0.d0 ; err_melt_isf( 95) =   0.d0 ; area_isf( 95) = 1.d-9
name_isf( 96) = '                   '; name_reg( 96) = '                   '; melt_isf( 96) =   0.d0 ; err_melt_isf( 96) =   0.d0 ; area_isf( 96) = 1.d-9
name_isf( 97) = '                   '; name_reg( 97) = '                   '; melt_isf( 97) =   0.d0 ; err_melt_isf( 97) =   0.d0 ; area_isf( 97) = 1.d-9
name_isf( 98) = '                   '; name_reg( 98) = '                   '; melt_isf( 98) =   0.d0 ; err_melt_isf( 98) =   0.d0 ; area_isf( 98) = 1.d-9
name_isf( 99) = '                   '; name_reg( 99) = '                   '; melt_isf( 99) =   0.d0 ; err_melt_isf( 99) =   0.d0 ; area_isf( 99) = 1.d-9
name_isf(100) = '                   '; name_reg(100) = '                   '; melt_isf(100) =   0.d0 ; err_melt_isf(100) =   0.d0 ; area_isf(100) = 1.d-9
name_isf(101) = '                   '; name_reg(101) = '                   '; melt_isf(101) =   0.d0 ; err_melt_isf(101) =   0.d0 ; area_isf(101) = 1.d-9
name_isf(102) = '                   '; name_reg(102) = '                   '; melt_isf(102) =   0.d0 ; err_melt_isf(102) =   0.d0 ; area_isf(102) = 1.d-9
name_isf(103) = '                   '; name_reg(103) = '                   '; melt_isf(103) =   0.d0 ; err_melt_isf(103) =   0.d0 ; area_isf(103) = 1.d-9
name_isf(104) = '                   '; name_reg(104) = '                   '; melt_isf(104) =   0.d0 ; err_melt_isf(104) =   0.d0 ; area_isf(104) = 1.d-9
name_isf(105) = '                   '; name_reg(105) = '                   '; melt_isf(105) =   0.d0 ; err_melt_isf(105) =   0.d0 ; area_isf(105) = 1.d-9
name_isf(106) = '                   '; name_reg(106) = '                   '; melt_isf(106) =   0.d0 ; err_melt_isf(106) =   0.d0 ; area_isf(106) = 1.d-9
name_isf(107) = '                   '; name_reg(107) = '                   '; melt_isf(107) =   0.d0 ; err_melt_isf(107) =   0.d0 ; area_isf(107) = 1.d-9
name_isf(108) = '                   '; name_reg(108) = '                   '; melt_isf(108) =   0.d0 ; err_melt_isf(108) =   0.d0 ; area_isf(108) = 1.d-9
name_isf(109) = '                   '; name_reg(109) = '                   '; melt_isf(109) =   0.d0 ; err_melt_isf(109) =   0.d0 ; area_isf(109) = 1.d-9
name_isf(110) = '                   '; name_reg(110) = '                   '; melt_isf(110) =   0.d0 ; err_melt_isf(110) =   0.d0 ; area_isf(110) = 1.d-9
name_isf(111) = '                   '; name_reg(111) = '                   '; melt_isf(111) =   0.d0 ; err_melt_isf(111) =   0.d0 ; area_isf(111) = 1.d-9
name_isf(112) = '                   '; name_reg(112) = '                   '; melt_isf(112) =   0.d0 ; err_melt_isf(112) =   0.d0 ; area_isf(112) = 1.d-9
name_isf(113) = '                   '; name_reg(113) = '                   '; melt_isf(113) =   0.d0 ; err_melt_isf(113) =   0.d0 ; area_isf(113) = 1.d-9
name_isf(114) = '                   '; name_reg(114) = '                   '; melt_isf(114) =   0.d0 ; err_melt_isf(114) =   0.d0 ; area_isf(114) = 1.d-9
name_isf(115) = '                   '; name_reg(115) = '                   '; melt_isf(115) =   0.d0 ; err_melt_isf(115) =   0.d0 ; area_isf(115) = 1.d-9
name_isf(116) = '                   '; name_reg(116) = '                   '; melt_isf(116) =   0.d0 ; err_melt_isf(116) =   0.d0 ; area_isf(116) = 1.d-9
name_isf(117) = '                   '; name_reg(117) = '                   '; melt_isf(117) =   0.d0 ; err_melt_isf(117) =   0.d0 ; area_isf(117) = 1.d-9
name_isf(118) = '                   '; name_reg(118) = '                   '; melt_isf(118) =   0.d0 ; err_melt_isf(118) =   0.d0 ; area_isf(118) = 1.d-9
name_isf(119) = '                   '; name_reg(119) = '                   '; melt_isf(119) =   0.d0 ; err_melt_isf(119) =   0.d0 ; area_isf(119) = 1.d-9
name_isf(120) = '                   '; name_reg(120) = '                   '; melt_isf(120) =   0.d0 ; err_melt_isf(120) =   0.d0 ; area_isf(120) = 1.d-9
name_isf(121) = '                   '; name_reg(121) = '                   '; melt_isf(121) =   0.d0 ; err_melt_isf(121) =   0.d0 ; area_isf(121) = 1.d-9
name_isf(122) = '                   '; name_reg(122) = '                   '; melt_isf(122) =   0.d0 ; err_melt_isf(122) =   0.d0 ; area_isf(122) = 1.d-9
name_isf(123) = '                   '; name_reg(123) = '                   '; melt_isf(123) =   0.d0 ; err_melt_isf(123) =   0.d0 ; area_isf(123) = 1.d-9
name_isf(124) = '                   '; name_reg(124) = '                   '; melt_isf(124) =   0.d0 ; err_melt_isf(124) =   0.d0 ; area_isf(124) = 1.d-9
name_isf(125) = '                   '; name_reg(125) = '                   '; melt_isf(125) =   0.d0 ; err_melt_isf(125) =   0.d0 ; area_isf(125) = 1.d-9
name_isf(126) = '                   '; name_reg(126) = '                   '; melt_isf(126) =   0.d0 ; err_melt_isf(126) =   0.d0 ; area_isf(126) = 1.d-9
name_isf(127) = '                   '; name_reg(127) = '                   '; melt_isf(127) =   0.d0 ; err_melt_isf(127) =   0.d0 ; area_isf(127) = 1.d-9
name_isf(128) = '                   '; name_reg(128) = '                   '; melt_isf(128) =   0.d0 ; err_melt_isf(128) =   0.d0 ; area_isf(128) = 1.d-9
name_isf(129) = '                   '; name_reg(129) = '                   '; melt_isf(129) =   0.d0 ; err_melt_isf(129) =   0.d0 ; area_isf(129) = 1.d-9
name_isf(130) = '                   '; name_reg(130) = '                   '; melt_isf(130) =   0.d0 ; err_melt_isf(130) =   0.d0 ; area_isf(130) = 1.d-9
name_isf(131) = '                   '; name_reg(131) = '                   '; melt_isf(131) =   0.d0 ; err_melt_isf(131) =   0.d0 ; area_isf(131) = 1.d-9
name_isf(132) = '                   '; name_reg(132) = '                   '; melt_isf(132) =   0.d0 ; err_melt_isf(132) =   0.d0 ; area_isf(132) = 1.d-9
name_isf(133) = '                   '; name_reg(133) = '                   '; melt_isf(133) =   0.d0 ; err_melt_isf(133) =   0.d0 ; area_isf(133) = 1.d-9
name_isf(134) = '                   '; name_reg(134) = '                   '; melt_isf(134) =   0.d0 ; err_melt_isf(134) =   0.d0 ; area_isf(134) = 1.d-9
name_isf(135) = '                   '; name_reg(135) = '                   '; melt_isf(135) =   0.d0 ; err_melt_isf(135) =   0.d0 ; area_isf(135) = 1.d-9
name_isf(136) = '                   '; name_reg(136) = '                   '; melt_isf(136) =   0.d0 ; err_melt_isf(136) =   0.d0 ; area_isf(136) = 1.d-9
name_isf(137) = '                   '; name_reg(137) = '                   '; melt_isf(137) =   0.d0 ; err_melt_isf(137) =   0.d0 ; area_isf(137) = 1.d-9
name_isf(138) = '                   '; name_reg(138) = '                   '; melt_isf(138) =   0.d0 ; err_melt_isf(138) =   0.d0 ; area_isf(138) = 1.d-9
name_isf(139) = '                   '; name_reg(139) = '                   '; melt_isf(139) =   0.d0 ; err_melt_isf(139) =   0.d0 ; area_isf(139) = 1.d-9
name_isf(140) = '                   '; name_reg(140) = '                   '; melt_isf(140) =   0.d0 ; err_melt_isf(140) =   0.d0 ; area_isf(140) = 1.d-9
name_isf(141) = '                   '; name_reg(141) = '                   '; melt_isf(141) =   0.d0 ; err_melt_isf(141) =   0.d0 ; area_isf(141) = 1.d-9
name_isf(142) = '                   '; name_reg(142) = '                   '; melt_isf(142) =   0.d0 ; err_melt_isf(142) =   0.d0 ; area_isf(142) = 1.d-9
name_isf(143) = '                   '; name_reg(143) = '                   '; melt_isf(143) =   0.d0 ; err_melt_isf(143) =   0.d0 ; area_isf(143) = 1.d-9
name_isf(144) = '                   '; name_reg(144) = '                   '; melt_isf(144) =   0.d0 ; err_melt_isf(144) =   0.d0 ; area_isf(144) = 1.d-9
name_isf(145) = '                   '; name_reg(145) = '                   '; melt_isf(145) =   0.d0 ; err_melt_isf(145) =   0.d0 ; area_isf(145) = 1.d-9
name_isf(146) = '                   '; name_reg(146) = '                   '; melt_isf(146) =   0.d0 ; err_melt_isf(146) =   0.d0 ; area_isf(146) = 1.d-9
name_isf(147) = '                   '; name_reg(147) = '                   '; melt_isf(147) =   0.d0 ; err_melt_isf(147) =   0.d0 ; area_isf(147) = 1.d-9
name_isf(148) = '                   '; name_reg(148) = '                   '; melt_isf(148) =   0.d0 ; err_melt_isf(148) =   0.d0 ; area_isf(148) = 1.d-9
name_isf(149) = '                   '; name_reg(149) = '                   '; melt_isf(149) =   0.d0 ; err_melt_isf(149) =   0.d0 ; area_isf(149) = 1.d-9
name_isf(150) = '                   '; name_reg(150) = '                   '; melt_isf(150) =   0.d0 ; err_melt_isf(150) =   0.d0 ; area_isf(150) = 1.d-9
name_isf(151) = '                   '; name_reg(151) = '                   '; melt_isf(151) =   0.d0 ; err_melt_isf(151) =   0.d0 ; area_isf(151) = 1.d-9
name_isf(152) = '                   '; name_reg(152) = '                   '; melt_isf(152) =   0.d0 ; err_melt_isf(152) =   0.d0 ; area_isf(152) = 1.d-9
name_isf(153) = '                   '; name_reg(153) = '                   '; melt_isf(153) =   0.d0 ; err_melt_isf(153) =   0.d0 ; area_isf(153) = 1.d-9
name_isf(154) = '                   '; name_reg(154) = '                   '; melt_isf(154) =   0.d0 ; err_melt_isf(154) =   0.d0 ; area_isf(154) = 1.d-9
name_isf(155) = '                   '; name_reg(155) = '                   '; melt_isf(155) =   0.d0 ; err_melt_isf(155) =   0.d0 ; area_isf(155) = 1.d-9
name_isf(156) = '                   '; name_reg(156) = '                   '; melt_isf(156) =   0.d0 ; err_melt_isf(156) =   0.d0 ; area_isf(156) = 1.d-9
name_isf(157) = '                   '; name_reg(157) = '                   '; melt_isf(157) =   0.d0 ; err_melt_isf(157) =   0.d0 ; area_isf(157) = 1.d-9
name_isf(158) = '                   '; name_reg(158) = '                   '; melt_isf(158) =   0.d0 ; err_melt_isf(158) =   0.d0 ; area_isf(158) = 1.d-9
name_isf(159) = '                   '; name_reg(159) = '                   '; melt_isf(159) =   0.d0 ; err_melt_isf(159) =   0.d0 ; area_isf(159) = 1.d-9
name_isf(160) = '                   '; name_reg(160) = '                   '; melt_isf(160) =   0.d0 ; err_melt_isf(160) =   0.d0 ; area_isf(160) = 1.d-9
name_isf(161) = '                   '; name_reg(161) = '                   '; melt_isf(161) =   0.d0 ; err_melt_isf(161) =   0.d0 ; area_isf(161) = 1.d-9
name_isf(162) = '                   '; name_reg(162) = '                   '; melt_isf(162) =   0.d0 ; err_melt_isf(162) =   0.d0 ; area_isf(162) = 1.d-9
name_isf(163) = '                   '; name_reg(163) = '                   '; melt_isf(163) =   0.d0 ; err_melt_isf(163) =   0.d0 ; area_isf(163) = 1.d-9
name_isf(164) = '                   '; name_reg(164) = '                   '; melt_isf(164) =   0.d0 ; err_melt_isf(164) =   0.d0 ; area_isf(164) = 1.d-9
name_isf(165) = '                   '; name_reg(165) = '                   '; melt_isf(165) =   0.d0 ; err_melt_isf(165) =   0.d0 ; area_isf(165) = 1.d-9
name_isf(166) = '                   '; name_reg(166) = '                   '; melt_isf(166) =   0.d0 ; err_melt_isf(166) =   0.d0 ; area_isf(166) = 1.d-9
name_isf(167) = '                   '; name_reg(167) = '                   '; melt_isf(167) =   0.d0 ; err_melt_isf(167) =   0.d0 ; area_isf(167) = 1.d-9
name_isf(168) = '                   '; name_reg(168) = '                   '; melt_isf(168) =   0.d0 ; err_melt_isf(168) =   0.d0 ; area_isf(168) = 1.d-9
name_isf(169) = '                   '; name_reg(169) = '                   '; melt_isf(169) =   0.d0 ; err_melt_isf(169) =   0.d0 ; area_isf(169) = 1.d-9
name_isf(170) = '                   '; name_reg(170) = '                   '; melt_isf(170) =   0.d0 ; err_melt_isf(170) =   0.d0 ; area_isf(170) = 1.d-9
name_isf(171) = '                   '; name_reg(171) = '                   '; melt_isf(171) =   0.d0 ; err_melt_isf(171) =   0.d0 ; area_isf(171) = 1.d-9
name_isf(172) = '                   '; name_reg(172) = '                   '; melt_isf(172) =   0.d0 ; err_melt_isf(172) =   0.d0 ; area_isf(172) = 1.d-9
name_isf(173) = '                   '; name_reg(173) = '                   '; melt_isf(173) =   0.d0 ; err_melt_isf(173) =   0.d0 ; area_isf(173) = 1.d-9
name_isf(174) = '                   '; name_reg(174) = '                   '; melt_isf(174) =   0.d0 ; err_melt_isf(174) =   0.d0 ; area_isf(174) = 1.d-9
name_isf(175) = '                   '; name_reg(175) = '                   '; melt_isf(175) =   0.d0 ; err_melt_isf(175) =   0.d0 ; area_isf(175) = 1.d-9
name_isf(176) = '                   '; name_reg(176) = '                   '; melt_isf(176) =   0.d0 ; err_melt_isf(176) =   0.d0 ; area_isf(176) = 1.d-9
name_isf(177) = '                   '; name_reg(177) = '                   '; melt_isf(177) =   0.d0 ; err_melt_isf(177) =   0.d0 ; area_isf(177) = 1.d-9
name_isf(178) = '                   '; name_reg(178) = '                   '; melt_isf(178) =   0.d0 ; err_melt_isf(178) =   0.d0 ; area_isf(178) = 1.d-9
name_isf(179) = '                   '; name_reg(179) = '                   '; melt_isf(179) =   0.d0 ; err_melt_isf(179) =   0.d0 ; area_isf(179) = 1.d-9
name_isf(180) = '                   '; name_reg(180) = '                   '; melt_isf(180) =   0.d0 ; err_melt_isf(180) =   0.d0 ; area_isf(180) = 1.d-9
name_isf(181) = '                   '; name_reg(181) = '                   '; melt_isf(181) =   0.d0 ; err_melt_isf(181) =   0.d0 ; area_isf(181) = 1.d-9
name_isf(182) = '                   '; name_reg(182) = '                   '; melt_isf(182) =   0.d0 ; err_melt_isf(182) =   0.d0 ; area_isf(182) = 1.d-9
name_isf(183) = '                   '; name_reg(183) = '                   '; melt_isf(183) =   0.d0 ; err_melt_isf(183) =   0.d0 ; area_isf(183) = 1.d-9
name_isf(184) = '                   '; name_reg(184) = '                   '; melt_isf(184) =   0.d0 ; err_melt_isf(184) =   0.d0 ; area_isf(184) = 1.d-9
name_isf(185) = '                   '; name_reg(185) = '                   '; melt_isf(185) =   0.d0 ; err_melt_isf(185) =   0.d0 ; area_isf(185) = 1.d-9
name_isf(186) = '                   '; name_reg(186) = '                   '; melt_isf(186) =   0.d0 ; err_melt_isf(186) =   0.d0 ; area_isf(186) = 1.d-9
name_isf(187) = '                   '; name_reg(187) = '                   '; melt_isf(187) =   0.d0 ; err_melt_isf(187) =   0.d0 ; area_isf(187) = 1.d-9
name_isf(188) = '                   '; name_reg(188) = '                   '; melt_isf(188) =   0.d0 ; err_melt_isf(188) =   0.d0 ; area_isf(188) = 1.d-9
name_isf(189) = '                   '; name_reg(189) = '                   '; melt_isf(189) =   0.d0 ; err_melt_isf(189) =   0.d0 ; area_isf(189) = 1.d-9
name_isf(190) = '                   '; name_reg(190) = '                   '; melt_isf(190) =   0.d0 ; err_melt_isf(190) =   0.d0 ; area_isf(190) = 1.d-9
name_isf(191) = '                   '; name_reg(191) = '                   '; melt_isf(191) =   0.d0 ; err_melt_isf(191) =   0.d0 ; area_isf(191) = 1.d-9
name_isf(192) = '                   '; name_reg(192) = '                   '; melt_isf(192) =   0.d0 ; err_melt_isf(192) =   0.d0 ; area_isf(192) = 1.d-9
name_isf(193) = '                   '; name_reg(193) = '                   '; melt_isf(193) =   0.d0 ; err_melt_isf(193) =   0.d0 ; area_isf(193) = 1.d-9
name_isf(194) = '                   '; name_reg(194) = '                   '; melt_isf(194) =   0.d0 ; err_melt_isf(194) =   0.d0 ; area_isf(194) = 1.d-9
name_isf(195) = '                   '; name_reg(195) = '                   '; melt_isf(195) =   0.d0 ; err_melt_isf(195) =   0.d0 ; area_isf(195) = 1.d-9
name_isf(196) = '                   '; name_reg(196) = '                   '; melt_isf(196) =   0.d0 ; err_melt_isf(196) =   0.d0 ; area_isf(196) = 1.d-9
name_isf(197) = '                   '; name_reg(197) = '                   '; melt_isf(197) =   0.d0 ; err_melt_isf(197) =   0.d0 ; area_isf(197) = 1.d-9
name_isf(198) = '                   '; name_reg(198) = '                   '; melt_isf(198) =   0.d0 ; err_melt_isf(198) =   0.d0 ; area_isf(198) = 1.d-9
name_isf(199) = '                   '; name_reg(199) = '                   '; melt_isf(199) =   0.d0 ; err_melt_isf(199) =   0.d0 ; area_isf(199) = 1.d-9
name_isf(200) = '                   '; name_reg(200) = '                   '; melt_isf(200) =   0.d0 ; err_melt_isf(200) =   0.d0 ; area_isf(200) = 1.d-9
name_isf(201) = '                   '; name_reg(201) = '                   '; melt_isf(201) =   0.d0 ; err_melt_isf(201) =   0.d0 ; area_isf(201) = 1.d-9
name_isf(202) = '                   '; name_reg(202) = '                   '; melt_isf(202) =   0.d0 ; err_melt_isf(202) =   0.d0 ; area_isf(202) = 1.d-9
name_isf(203) = '                   '; name_reg(203) = '                   '; melt_isf(203) =   0.d0 ; err_melt_isf(203) =   0.d0 ; area_isf(203) = 1.d-9
name_isf(204) = '                   '; name_reg(204) = '                   '; melt_isf(204) =   0.d0 ; err_melt_isf(204) =   0.d0 ; area_isf(204) = 1.d-9
name_isf(205) = '                   '; name_reg(205) = '                   '; melt_isf(205) =   0.d0 ; err_melt_isf(205) =   0.d0 ; area_isf(205) = 1.d-9
name_isf(206) = '                   '; name_reg(206) = '                   '; melt_isf(206) =   0.d0 ; err_melt_isf(206) =   0.d0 ; area_isf(206) = 1.d-9
name_isf(207) = '                   '; name_reg(207) = '                   '; melt_isf(207) =   0.d0 ; err_melt_isf(207) =   0.d0 ; area_isf(207) = 1.d-9
name_isf(208) = '                   '; name_reg(208) = '                   '; melt_isf(208) =   0.d0 ; err_melt_isf(208) =   0.d0 ; area_isf(208) = 1.d-9
name_isf(209) = '                   '; name_reg(209) = '                   '; melt_isf(209) =   0.d0 ; err_melt_isf(209) =   0.d0 ; area_isf(209) = 1.d-9
name_isf(210) = '                   '; name_reg(210) = '                   '; melt_isf(210) =   0.d0 ; err_melt_isf(210) =   0.d0 ; area_isf(210) = 1.d-9
name_isf(211) = '                   '; name_reg(211) = '                   '; melt_isf(211) =   0.d0 ; err_melt_isf(211) =   0.d0 ; area_isf(211) = 1.d-9
name_isf(212) = '                   '; name_reg(212) = '                   '; melt_isf(212) =   0.d0 ; err_melt_isf(212) =   0.d0 ; area_isf(212) = 1.d-9
name_isf(213) = '                   '; name_reg(213) = '                   '; melt_isf(213) =   0.d0 ; err_melt_isf(213) =   0.d0 ; area_isf(213) = 1.d-9

!---------------------------------------------------------------------------------
! Correction to account for different areas in Rignot et al. and in our dataset :

melt_isf(:) = melt_isf(:) * our_isf_area(:) / area_isf(:)

do kisf=10,75
  write(*,*) '~~~~ ', kisf, our_isf_area(kisf) / area_isf(kisf)
enddo

!-----------------------------------------------
! Attribute special values for pinning points :
!
! Pinning points are defined as grounded points fully surrounded
! by open ocean and ice shelves.
!
! NB : islands can be counted as tipping points.

ALLOCATE( tmppin(mx,my), rmppin(mx,my) )

do kisf=2,Nisf

   write(*,*) 'processing ', kisf

   llisf=.false.

   ! Bounds for indices :
   itmpmin= 1.d12; imin = 1
   itmpmax=-1.d12; imax = mx
   jtmpmin= 1.d12; jmin = 1
   jtmpmax=-1.d12; jmax = my
   do ii=1,mx
   do jj=1,my
     if ( isfmask(ii,jj) .eq. kisf ) then
       llisf=.true.
       if ( x(ii) .lt. itmpmin ) then
         itmpmin=x(ii)
         imin=ii
       endif
       if ( x(ii) .gt. itmpmax ) then
         itmpmax=x(ii)
         imax=ii
       endif
       if ( y(jj) .lt. jtmpmin ) then
         jtmpmin=y(jj)
         jmin=jj
       endif
       if ( y(jj) .gt. jtmpmax ) then
         jtmpmax=y(jj)
         jmax=jj
       endif
     endif
   enddo
   enddo
   imin=MAX(1,imin-1)
   imax=MIN(mx,imax+1)
   jmin=MAX(1,jmin-1)
   jmax=MIN(my,jmax+1)

   if ( llisf ) then
     tmppin(:,:)=0
     do ii=imin,imax
       if ( isfmask(ii,jmin) .eq. 0 ) tmppin(ii,jmin)= 1 
       if ( isfmask(ii,jmax) .eq. 0 ) tmppin(ii,jmax)= 1 
     enddo
     do rs=1,MAX(imax-imin,jmax-jmin)
        rmppin(:,:)=tmppin(:,:)
        do ii=imin,imax
        do jj=jmin+1,jmax-1
          if ( isfmask(ii,jj).ne.kisf .and. isfmask(ii,jj).ne. 1) then  ! this ice shelf or open ocean
            if ( tmppin(MIN(ii+1,imax),jj).eq.1 ) rmppin(ii,jj)=1
            if ( tmppin(MAX(ii-1,imin),jj).eq.1 ) rmppin(ii,jj)=1  
            if ( tmppin(ii,MIN(jj+1,jmax)).eq.1 ) rmppin(ii,jj)=1  
            if ( tmppin(ii,MAX(jj-1,jmin)).eq.1 ) rmppin(ii,jj)=1  
          endif
        enddo
        enddo
        tmppin(:,:)=rmppin(:,:)
     enddo
     ! Attribut negative ice shelf ID in case of pinning point
     do ii=imin,imax
     do jj=jmin,jmax
       if ( isfmask(ii,jj).eq.0 .and. tmppin(ii,jj).eq.0 ) then
         isfmask(ii,jj) = -kisf
       endif
     enddo
     enddo
     !! MANUAL MODIFICATIONS :
     if ( kisf .eq. 22 ) then ! Getz
       do ii=imin,imax
       do jj=jmin,MIN(my,jmax+40)
         if ( lon(ii,jj) .ge. -130.0 .and. lon(ii,jj) .lt. -120.0 .and. lat(ii,jj) .ge. -74.6 .and. isfmask(ii,jj) .eq. 0 ) isfmask(ii,jj) = -kisf
       enddo
       enddo
     endif
   endif

enddo ! kisf

DEALLOCATE( tmppin, rmppin )

!---------------------------------------                      
! Calculate nearest ocean points

ALLOCATE( GLINE(mx,my), FRONT(mx,my) )

ALLOCATE( front_bot_dep_max(Nisf), front_bot_dep_avg(Nisf), front_ice_dep_min(Nisf) )
ALLOCATE( front_ice_dep_avg(Nisf), front_min_lon(Nisf), front_max_lon(Nisf), front_min_lat(Nisf), front_max_lat(Nisf) )

GLINE(:,:) = 0
FRONT(:,:) = 0

front_bot_dep_max(:) = 0.d0
front_bot_dep_avg(:) = 0.d0
front_ice_dep_min(:) = 1.d5
front_ice_dep_avg(:) = 0.d0
front_min_lon(:) =  180.d0
front_max_lon(:) = -180.d0
front_min_lat(:) =   90.d0
front_max_lat(:) =  -90.d0

do kisf=2,Nisf

   rn_ice = 0.d0

   do ii=1,mx
   do jj=1,my
 
      iip1=MIN(ii+1,mx)
      iim1=MAX(ii-1,1)
      jjp1=MIN(jj+1,my)
      jjm1=MAX(jj-1,1)

      if ( isfmask(ii,jj) .eq. kisf                                            &
      &    .and. (      isfmask(iip1,jj) .eq. 0 .or. isfmask(iim1,jj) .eq. 0   &
      &            .or. isfmask(ii,jjp1) .eq. 0 .or. isfmask(ii,jjm1) .eq. 0 ) ) then
        GLINE(ii,jj) = kisf
      endif

      if ( isfmask(ii,jj) .eq. kisf                                            &
      &    .and. (      isfmask(iip1,jj) .eq. 1 .or. isfmask(iim1,jj) .eq. 1   &
      &            .or. isfmask(ii,jjp1) .eq. 1 .or. isfmask(ii,jjm1) .eq. 1 ) ) then
        FRONT(ii,jj) = kisf
        front_bot_dep_max(kisf) = MAX( front_bot_dep_max(kisf),  bedrock_topography(ii,jj)*(-1.0) )
        front_bot_dep_avg(kisf) =      front_bot_dep_avg(kisf) - bedrock_topography(ii,jj)
        front_ice_dep_min(kisf) = MIN( front_ice_dep_min(kisf),  ice_base_topography(ii,jj)*(-1.0) )
        front_ice_dep_avg(kisf) =      front_ice_dep_avg(kisf) - ice_base_topography(ii,jj)
        if ( kisf .ne. 10 ) then
          front_min_lon(kisf) = MIN( front_min_lon(kisf), lon(ii,jj) )
          front_max_lon(kisf) = MAX( front_max_lon(kisf), lon(ii,jj) )
          front_min_lat(kisf) = MIN( front_min_lat(kisf), lat(ii,jj) )
          front_max_lat(kisf) = MAX( front_max_lat(kisf), lat(ii,jj) )
        else ! Ross
          if ( lon(ii,jj) .lt. 0.d0 ) then
            tmplon = 3.6d2 + lon(ii,jj)
            front_min_lon(kisf) = MIN( front_min_lon(kisf), tmplon  )
            front_max_lon(kisf) = MAX( front_max_lon(kisf), tmplon  )
          else
            tmplon = lon(ii,jj)
            front_min_lon(kisf) = MIN( front_min_lon(kisf), tmplon  )
            front_max_lon(kisf) = MAX( front_max_lon(kisf), tmplon  )
          endif
          front_min_lat(kisf) = MIN( front_min_lat(kisf), lat(ii,jj) )
          front_max_lat(kisf) = MAX( front_max_lat(kisf), lat(ii,jj) )
        endif
        rn_ice = rn_ice + 1.d0
      endif

   enddo
   enddo

   front_bot_dep_avg(kisf) = front_bot_dep_avg(kisf) / rn_ice
   front_ice_dep_avg(kisf) = front_ice_dep_avg(kisf) / rn_ice

enddo

! Ross and the date line :
if ( front_min_lon(10) .gt. 180.d0 ) front_min_lon(10) = front_min_lon(10) - 360.d0
if ( front_max_lon(10) .gt. 180.d0 ) front_max_lon(10) = front_max_lon(10) - 360.d0


DEALLOCATE( bedrock_topography, ice_base_topography )

write(*,*) 'A FEW CHECKS :  '
write(*,*) '~~~~~~~~~~~~~~  '
write(*,*) '                '
kisf=66 ! Pine Island
write(*,*) TRIM(name_isf(kisf)), ' (', TRIM(name_reg(kisf)), ')'
write(*,*) '  front_bot_dep_max = ', front_bot_dep_max(kisf)
write(*,*) '  front_bot_dep_avg = ', front_bot_dep_avg(kisf)
write(*,*) '  front_ice_dep_min = ', front_ice_dep_min(kisf)
write(*,*) '  front_ice_dep_avg = ', front_ice_dep_avg(kisf)
write(*,*) '  front_min_lon = ', front_min_lon(kisf)
write(*,*) '  front_max_lon = ', front_max_lon(kisf)
write(*,*) '  front_min_lat = ', front_min_lat(kisf)
write(*,*) '  front_max_lat = ', front_max_lat(kisf)
kisf=21 ! Ronne
write(*,*) TRIM(name_isf(kisf)), ' (', TRIM(name_reg(kisf)), ')'
write(*,*) '  front_bot_dep_max = ', front_bot_dep_max(kisf)
write(*,*) '  front_bot_dep_avg = ', front_bot_dep_avg(kisf)
write(*,*) '  front_ice_dep_min = ', front_ice_dep_min(kisf)
write(*,*) '  front_ice_dep_avg = ', front_ice_dep_avg(kisf)
write(*,*) '  front_min_lon = ', front_min_lon(kisf)
write(*,*) '  front_max_lon = ', front_max_lon(kisf)
write(*,*) '  front_min_lat = ', front_min_lat(kisf)
write(*,*) '  front_max_lat = ', front_max_lat(kisf)
kisf=10 ! Ross
write(*,*) TRIM(name_isf(kisf)), ' (', TRIM(name_reg(kisf)), ')'
write(*,*) '  front_bot_dep_max = ', front_bot_dep_max(kisf)
write(*,*) '  front_bot_dep_avg = ', front_bot_dep_avg(kisf)
write(*,*) '  front_ice_dep_min = ', front_ice_dep_min(kisf)
write(*,*) '  front_ice_dep_avg = ', front_ice_dep_avg(kisf)
write(*,*) '  front_min_lon = ', front_min_lon(kisf)
write(*,*) '  front_max_lon = ', front_max_lon(kisf)
write(*,*) '  front_min_lat = ', front_min_lat(kisf)
write(*,*) '  front_max_lat = ', front_max_lat(kisf)
kisf=31 ! Amery
write(*,*) TRIM(name_isf(kisf)), ' (', TRIM(name_reg(kisf)), ')'
write(*,*) '  front_bot_dep_max = ', front_bot_dep_max(kisf)
write(*,*) '  front_bot_dep_avg = ', front_bot_dep_avg(kisf)
write(*,*) '  front_ice_dep_min = ', front_ice_dep_min(kisf)
write(*,*) '  front_ice_dep_avg = ', front_ice_dep_avg(kisf)
write(*,*) '  front_min_lon = ', front_min_lon(kisf)
write(*,*) '  front_max_lon = ', front_max_lon(kisf)
write(*,*) '  front_min_lat = ', front_min_lat(kisf)
write(*,*) '  front_max_lat = ', front_max_lat(kisf)
write(*,*) ' '

!---------------------------------------                      
! Writing new netcdf file :                                   

write(*,*) 'Creating ', TRIM(file_msk)

status = NF90_CREATE(TRIM(file_msk),NF90_NOCLOBBER,fidM)
call erreur(status,.TRUE.,'create mask file')

status = NF90_DEF_DIM(fidM,"x",mx,dimID_x)
call erreur(status,.TRUE.,"def_dimID_londim")
status = NF90_DEF_DIM(fidM,"y",my,dimID_y)
call erreur(status,.TRUE.,"def_dimID_latdim")
status = NF90_DEF_DIM(fidM,"Nisf",Nisf,dimID_Nisf)
call erreur(status,.TRUE.,"def_dimID_Nisf")
status = NF90_DEF_DIM(fidM,"StrLen",strlen,dimID_strlen)
call erreur(status,.TRUE.,"def_dimID_strlen")

status = NF90_DEF_VAR(fidM,"x",NF90_DOUBLE,(/dimID_x/),x_ID)
call erreur(status,.TRUE.,"def_var_lon_ID")
status = NF90_DEF_VAR(fidM,"y",NF90_DOUBLE,(/dimID_y/),y_ID)
call erreur(status,.TRUE.,"def_var_lat_ID")
status = NF90_DEF_VAR(fidM,"lon",NF90_DOUBLE,(/dimID_x,dimID_y/),lon_ID)
call erreur(status,.TRUE.,"def_var_lon_ID")
status = NF90_DEF_VAR(fidM,"lat",NF90_DOUBLE,(/dimID_x,dimID_y/),lat_ID)
call erreur(status,.TRUE.,"def_var_lat_ID")
status = NF90_DEF_VAR(fidM,"name_isf",NF90_CHAR,(/dimID_strlen,dimID_Nisf/),name_isf_ID)
call erreur(status,.TRUE.,"def_var_name_isf_ID")
status = NF90_DEF_VAR(fidM,"name_reg",NF90_CHAR,(/dimID_strlen,dimID_Nisf/),name_reg_ID)
call erreur(status,.TRUE.,"def_var_name_reg_ID")
status = NF90_DEF_VAR(fidM,"isfmask",NF90_SHORT,(/dimID_x,dimID_y/),isfmask_ID)
call erreur(status,.TRUE.,"def_var_isfmask_ID")
status = NF90_DEF_VAR(fidM,"GL_mask",NF90_SHORT,(/dimID_x,dimID_y/),GL_mask_ID)
call erreur(status,.TRUE.,"def_var_GLmask_ID")
status = NF90_DEF_VAR(fidM,"IF_mask",NF90_SHORT,(/dimID_x,dimID_y/),IF_mask_ID)
call erreur(status,.TRUE.,"def_var_IF_mask_ID")
status = NF90_DEF_VAR(fidM,"front_bot_dep_max",NF90_DOUBLE,(/dimID_Nisf/),front_bot_dep_max_ID)
call erreur(status,.TRUE.,"def_var_front_bot_dep_max_ID")
status = NF90_DEF_VAR(fidM,"front_bot_dep_avg",NF90_DOUBLE,(/dimID_Nisf/),front_bot_dep_avg_ID)
call erreur(status,.TRUE.,"def_var_front_bot_dep_avg_ID")
status = NF90_DEF_VAR(fidM,"front_ice_dep_min",NF90_DOUBLE,(/dimID_Nisf/),front_ice_dep_min_ID)
call erreur(status,.TRUE.,"def_var_front_ice_dep_min_ID")
status = NF90_DEF_VAR(fidM,"front_ice_dep_avg",NF90_DOUBLE,(/dimID_Nisf/),front_ice_dep_avg_ID)
call erreur(status,.TRUE.,"def_var_front_ice_dep_avg_ID")
status = NF90_DEF_VAR(fidM,"front_min_lon",NF90_DOUBLE,(/dimID_Nisf/),front_min_lon_ID)
call erreur(status,.TRUE.,"def_var_front_min_lon_ID")
status = NF90_DEF_VAR(fidM,"front_max_lon",NF90_DOUBLE,(/dimID_Nisf/),front_max_lon_ID)
call erreur(status,.TRUE.,"def_var_front_max_lon_ID")
status = NF90_DEF_VAR(fidM,"front_min_lat",NF90_DOUBLE,(/dimID_Nisf/),front_min_lat_ID)
call erreur(status,.TRUE.,"def_var_front_min_lat_ID")
status = NF90_DEF_VAR(fidM,"front_max_lat",NF90_DOUBLE,(/dimID_Nisf/),front_max_lat_ID)
call erreur(status,.TRUE.,"def_var_front_max_lat_ID")
status = NF90_DEF_VAR(fidM,"melt_isf",NF90_DOUBLE,(/dimID_Nisf/),melt_isf_ID)
call erreur(status,.TRUE.,"def_var_melt_isf_ID")
status = NF90_DEF_VAR(fidM,"err_melt_isf",NF90_DOUBLE,(/dimID_Nisf/),err_melt_isf_ID)
call erreur(status,.TRUE.,"def_var_err_melt_isf_ID")
status = NF90_DEF_VAR(fidM,"area_isf",NF90_DOUBLE,(/dimID_Nisf/),area_isf_ID)
call erreur(status,.TRUE.,"def_var_area_isf_ID")

status = NF90_PUT_ATT(fidM,lon_ID,"standard_name","longitude")
call erreur(status,.TRUE.,"put_att_lon_ID")
status = NF90_PUT_ATT(fidM,lon_ID,"units","degrees_east")
call erreur(status,.TRUE.,"put_att_lon_ID")
status = NF90_PUT_ATT(fidM,lon_ID,"long_name","longitude coordinate")
call erreur(status,.TRUE.,"put_att_lon_ID")
status = NF90_PUT_ATT(fidM,lon_ID,"name","lon")
call erreur(status,.TRUE.,"put_att_lon_ID")
status = NF90_PUT_ATT(fidM,lat_ID,"standard_name","latitude")
call erreur(status,.TRUE.,"put_att_lat_ID")
status = NF90_PUT_ATT(fidM,lat_ID,"units","degrees_north")
call erreur(status,.TRUE.,"put_att_lat_ID")
status = NF90_PUT_ATT(fidM,lat_ID,"long_name","latitude coordinate")
call erreur(status,.TRUE.,"put_att_lat_ID")
status = NF90_PUT_ATT(fidM,lat_ID,"name","lat")
call erreur(status,.TRUE.,"put_att_lat_ID")
status = NF90_PUT_ATT(fidM,isfmask_ID,"coordinates","lat lon")
call erreur(status,.TRUE.,"put_att_isfmask_ID")
status = NF90_PUT_ATT(fidM,isfmask_ID,"units","-")
call erreur(status,.TRUE.,"put_att_isfmask_ID")
status = NF90_PUT_ATT(fidM,isfmask_ID,"long_name","ice shelf mask (0 for grounded, 1 for ocean, isf ID for ice shelves)")
call erreur(status,.TRUE.,"put_att_isfmask_ID")
status = NF90_PUT_ATT(fidM,isfmask_ID,"title","isfmask")
call erreur(status,.TRUE.,"put_att_isfmask_ID")
status = NF90_PUT_ATT(fidM,name_isf_ID,"name","name_isf")
call erreur(status,.TRUE.,"put_att_name_isf_ID")
status = NF90_PUT_ATT(fidM,name_isf_ID,"name","Ice shelf name")
call erreur(status,.TRUE.,"put_att_name_isf_ID")
status = NF90_PUT_ATT(fidM,name_reg_ID,"name","name_reg")
call erreur(status,.TRUE.,"put_att_name_reg_ID")
status = NF90_PUT_ATT(fidM,name_reg_ID,"name","Region")
call erreur(status,.TRUE.,"put_att_name_reg_ID")
status = NF90_PUT_ATT(fidM,GL_mask_ID,"coordinates","lat lon")
call erreur(status,.TRUE.,"put_att_GL_mask_ID")
status = NF90_PUT_ATT(fidM,GL_mask_ID,"units","-")
call erreur(status,.TRUE.,"put_att_GL_mask_ID")
status = NF90_PUT_ATT(fidM,GL_mask_ID,"long_name","grounding zone mask (isf ID for grounding line, 0 elsewhere)")
call erreur(status,.TRUE.,"put_att_GL_mask_ID")
status = NF90_PUT_ATT(fidM,GL_mask_ID,"title","GL mask")
call erreur(status,.TRUE.,"put_att_GL_mask_ID")
status = NF90_PUT_ATT(fidM,IF_mask_ID,"coordinates","lat lon")
call erreur(status,.TRUE.,"put_att_IF_mask_ID")
status = NF90_PUT_ATT(fidM,IF_mask_ID,"units","-")
call erreur(status,.TRUE.,"put_att_IF_mask_ID")
status = NF90_PUT_ATT(fidM,IF_mask_ID,"long_name","ice-shelf front mask (isf ID for ice-shelf front, 0 elsewhere)")
call erreur(status,.TRUE.,"put_att_IF_mask_ID")
status = NF90_PUT_ATT(fidM,IF_mask_ID,"title","IF mask")
call erreur(status,.TRUE.,"put_att_IF_mask_ID")
status = NF90_PUT_ATT(fidM,front_min_lon_ID,"units","-")
call erreur(status,.TRUE.,"put_att_front_min_lon_ID")
status = NF90_PUT_ATT(fidM,front_min_lon_ID,"long_name","minimum longitude of the ice shelf front")
call erreur(status,.TRUE.,"put_att_front_min_lon_ID")
status = NF90_PUT_ATT(fidM,front_min_lon_ID,"title","front min lon")
call erreur(status,.TRUE.,"put_att_front_min_lon_ID")
status = NF90_PUT_ATT(fidM,front_max_lon_ID,"units","-")
call erreur(status,.TRUE.,"put_att_front_max_lon_ID")
status = NF90_PUT_ATT(fidM,front_max_lon_ID,"long_name","maximum longitude of the ice shelf front")
call erreur(status,.TRUE.,"put_att_front_max_lon_ID")
status = NF90_PUT_ATT(fidM,front_max_lon_ID,"title","front max lon")
call erreur(status,.TRUE.,"put_att_front_max_lon_ID")
status = NF90_PUT_ATT(fidM,front_min_lat_ID,"units","-")
call erreur(status,.TRUE.,"put_att_front_min_lat_ID")
status = NF90_PUT_ATT(fidM,front_min_lat_ID,"latg_name","minimum latitude of the ice shelf front")
call erreur(status,.TRUE.,"put_att_front_min_lat_ID")
status = NF90_PUT_ATT(fidM,front_min_lat_ID,"title","front min lat")
call erreur(status,.TRUE.,"put_att_front_min_lat_ID")
status = NF90_PUT_ATT(fidM,front_max_lat_ID,"units","-")
call erreur(status,.TRUE.,"put_att_front_max_lat_ID")
status = NF90_PUT_ATT(fidM,front_max_lat_ID,"latg_name","maximum latitude of the ice shelf front")
call erreur(status,.TRUE.,"put_att_front_max_lat_ID")
status = NF90_PUT_ATT(fidM,front_max_lat_ID,"title","front max lat")
call erreur(status,.TRUE.,"put_att_front_max_lat_ID")
status = NF90_PUT_ATT(fidM,melt_isf_ID,"units","Gt/yr")
call erreur(status,.TRUE.,"put_att_melt_isf_ID")
status = NF90_PUT_ATT(fidM,melt_isf_ID,"long_name","Total ice shelf melting flux (Rignot 2013 + Depoorter 2013 for small ice shelves)")
call erreur(status,.TRUE.,"put_att_melt_isf_ID")
status = NF90_PUT_ATT(fidM,melt_isf_ID,"title","melt")
call erreur(status,.TRUE.,"put_att_melt_isf_ID")
status = NF90_PUT_ATT(fidM,err_melt_isf_ID,"units","Gt/yr")
call erreur(status,.TRUE.,"put_att_err_melt_isf_ID")
status = NF90_PUT_ATT(fidM,err_melt_isf_ID,"long_name","Uncertainty on ice shelf melting flux")
call erreur(status,.TRUE.,"put_att_err_melt_isf_ID")
status = NF90_PUT_ATT(fidM,err_melt_isf_ID,"title","melt uncertainty")
call erreur(status,.TRUE.,"put_att_err_melt_isf_ID")
status = NF90_PUT_ATT(fidM,area_isf_ID,"units","km2")
call erreur(status,.TRUE.,"put_att_area_isf_ID")
status = NF90_PUT_ATT(fidM,area_isf_ID,"long_name","Ice shelf area in our topography data")
call erreur(status,.TRUE.,"put_att_area_isf_ID")
status = NF90_PUT_ATT(fidM,area_isf_ID,"title","Ice shelf area")
call erreur(status,.TRUE.,"put_att_area_isf_ID")

status = NF90_PUT_ATT(fidM,NF90_GLOBAL,"history","Created using identify_ice_shelf_mask_BEDMAP2.f90")
call erreur(status,.TRUE.,"put_att_GLOBAL_ID")

status = NF90_ENDDEF(fidM)                   
call erreur(status,.TRUE.,"fin_definition") 

status = NF90_PUT_VAR(fidM,lon_ID,lon)
call erreur(status,.TRUE.,"var_lon_ID")
status = NF90_PUT_VAR(fidM,lat_ID,lat)
call erreur(status,.TRUE.,"var_lat_ID")
status = NF90_PUT_VAR(fidM,x_ID,x)
call erreur(status,.TRUE.,"var_x_ID")
status = NF90_PUT_VAR(fidM,y_ID,y)
call erreur(status,.TRUE.,"var_y_ID")
status = NF90_PUT_VAR(fidM,name_isf_ID,name_isf)
call erreur(status,.TRUE.,"var_name_isf_ID")
status = NF90_PUT_VAR(fidM,name_reg_ID,name_reg)
call erreur(status,.TRUE.,"var_name_reg_ID")
status = NF90_PUT_VAR(fidM,isfmask_ID,isfmask)
call erreur(status,.TRUE.,"var_isfmask_ID")
status = NF90_PUT_VAR(fidM,GL_mask_ID,GLINE)
call erreur(status,.TRUE.,"var_GLmask_ID")
status = NF90_PUT_VAR(fidM,IF_mask_ID,FRONT)
call erreur(status,.TRUE.,"var_IF_mask_ID")
status = NF90_PUT_VAR(fidM,front_bot_dep_max_ID,front_bot_dep_max)
call erreur(status,.TRUE.,"var_front_bot_dep_max_ID")
status = NF90_PUT_VAR(fidM,front_bot_dep_avg_ID,front_bot_dep_avg)
call erreur(status,.TRUE.,"var_front_bot_dep_avg_ID")
status = NF90_PUT_VAR(fidM,front_ice_dep_min_ID,front_ice_dep_min)
call erreur(status,.TRUE.,"var_front_ice_dep_min_ID")
status = NF90_PUT_VAR(fidM,front_ice_dep_avg_ID,front_ice_dep_avg)
call erreur(status,.TRUE.,"var_front_ice_dep_avg_ID")
status = NF90_PUT_VAR(fidM,front_min_lon_ID,front_min_lon)
call erreur(status,.TRUE.,"var_front_min_lon_ID")
status = NF90_PUT_VAR(fidM,front_max_lon_ID,front_max_lon)
call erreur(status,.TRUE.,"var_front_max_lon_ID")
status = NF90_PUT_VAR(fidM,front_min_lat_ID,front_min_lat)
call erreur(status,.TRUE.,"var_front_min_lat_ID")
status = NF90_PUT_VAR(fidM,front_max_lat_ID,front_max_lat)
call erreur(status,.TRUE.,"var_front_max_lat_ID")
status = NF90_PUT_VAR(fidM,melt_isf_ID,melt_isf)
call erreur(status,.TRUE.,"var_melt_isf_ID")
status = NF90_PUT_VAR(fidM,err_melt_isf_ID,err_melt_isf)
call erreur(status,.TRUE.,"var_err_melt_isf_ID")
status = NF90_PUT_VAR(fidM,area_isf_ID,our_isf_area)
call erreur(status,.TRUE.,"var_area_isf_ID")

status = NF90_CLOSE(fidM)                    
call erreur(status,.TRUE.,"final")         

end program modif

!!================================================================
!SUBROUTINE stereo_loc_xy_to_latlon(px,py,x_in,y_in,lat_out,lon_out,RT,ex,lat_true,lon_posy)
!!---------------------------------------------------------------------------------------------
!! Writen from matlab's function [phi,lambda]=polarstereo_inv(x,y,a,e,phi_c,lambda_0)
!! Purpose: transforms lat/lon data to map coordinates for a polar stereographic system
!INTEGER,                  INTENT(IN)  :: px, py            !- matrix dimensions
!REAL*8, DIMENSION(px,py), INTENT(IN)  :: x_in, y_in        !- input x,y coordinates
!REAL*8, DIMENSION(px,py), INTENT(OUT) :: lat_out, lon_out  !- outut longitude and latitude (degrees)
!REAL*8, OPTIONAL,         INTENT(IN)  :: RT,             & !- Earth radius (m)                     [default = WGS84]
!&                                        ex,             & !- Earth's misshapenness (excentricity) [default = WGS84]
!&                                        lat_true,       & !- latitude of true scale in degrees    [default = 71 degS]
!&                                        lon_posy          !- meridian in degrees along the positive Y axis of the map [default = 0 deg]
!REAL*8, DIMENSION(px,py)              :: rho, t, chi, x, y
!REAL*8                                :: t_c, m_c, a, e, lat_c, lon_0
!INTEGER                               :: pm
!
!!- Optional arguments :
!if ( present(RT)       ) then; a    = RT       ; else; a     = 6378137.0        ; endif
!if ( present(ex)       ) then; e    = ex       ; else; e     =       0.08181919 ; endif
!if ( present(lat_true) ) then; lat_c= lat_true ; else; lat_c =     -71.0        ; endif
!if ( present(lon_posy) ) then; lon_0= lon_posy ; else; lon_0 =       0.0        ; endif
!
!!- input checking :
!if ( a .lt. 6.e6 .or. a .gt. 7.e6 ) then
!  write(*,*) '~!@#$%^* error : Earth radius is wrong >>>>>> stop !!'
!  stop
!endif
!!-
!if ( e .lt. 0.0 .or. e .gt. 1.0 ) then
!  write(*,*) '~!@#$%^* error : Earth excentricity is wrong >>>>>> stop !!' 
!  stop
!endif
!!-
!if ( lat_c .lt. -90.0 .or. lat_c .gt. 90.0 ) then
!  write(*,*) '~!@#$%^* error : latitude of true scale is wrong >>>>>> stop !!'
!  stop
!endif
!!-
!if ( lon_0 .lt. -180.0 .or. lon_0 .gt. 360.0 ) then
!  write(*,*) '~!@#$%^* error : meridian of positive Y axis is wrong >>>>>> stop !!'
!  stop
!endif
!
!!- convert to radians :
!lat_c    = deg2rad * lat_c
!lon_0    = deg2rad * lon_0
!
!!- if the standard parallel is in S.Hemi., switch signs.
!if ( lat_c .lt. 0.0 ) then 
!    pm     = -1    ! plus or minus, north lat. or south
!    lat_c  = -lat_c
!    lon_0  = -lon_0
!    x(:,:) = -x_in(:,:)
!    y(:,:) = -y_in(:,:)
!else
!    pm     = 1
!endif
!
!!- See Snyder for details.
!t_c      = tan(pi/4-lat_c   /2) / ( (1-e*sin(lat_c   )) / (1+e*sin(lat_c   )) )**(e/2)
!m_c      = cos(lat_c) / sqrt( 1 - e**2 * (sin(lat_c))**2 )
!rho(:,:) = sqrt(x(:,:)**2+y(:,:)**2)
!t(:,:)   = rho(:,:) * t_c / ( a * m_c )
!
!!- iterate to find phi, with a threshold of pi*1e-8
!!  phi_alt=pi/2 - 2 * atan(t); %guess for phi
!!  phiold=phi_alt+10; %+10 to make sure it executes the while loop
!!  while any(abs(phi_alt(:)-phiold(:)) > pi*1e-8)
!!     phiold=phi_alt;
!!     phi_alt=pi/2 - 2*atan(t.*((1-e*sin(phi_alt))./(1+e*sin(phi_alt))).^(e/2));
!!     %add a break in case it doesn't converge?
!!  end
!
!chi(:,:) = 0.5*pi - 2 * atan(t(:,:))  !- find lat with a series instead of iterating.
!
!lat_out(:,:) = chi(:,:) + ( (1./2.) * e**2 + (5./24.) * e**4 + ( 1./ 12.) * e**6 + (  13./   360.) * e**8 ) * sin(2*chi(:,:)) &
!&                       + (                  (7./48.) * e**4 + (29./240.) * e**6 + ( 811./ 11520.) * e**8 ) * sin(4*chi(:,:)) &
!&                       + (                                  + ( 7./120.) * e**6 + (  81./  1120.) * e**8 ) * sin(6*chi(:,:)) &
!&                       + (                                                      + (4279./161280.) * e**8 ) * sin(8*chi(:,:))
!lon_out(:,:) = lon_0 + atan2(x,-y)
!
!!- correct the signs and phasing :
!lat_out(:,:) = pm * lat_out(:,:)
!lon_out(:,:) = pm * lon_out(:,:)
!lon_out(:,:) = mod(lon_out(:,:)+pi,2*pi)-pi !- want longitude in the range -pi to pi
!
!!- convert back to degrees :
!lat_out(:,:) = rad2deg * lat_out(:,:)
!lon_out(:,:) = rad2deg * lon_out(:,:)
!
!END SUBROUTINE stereo_loc_xy_to_latlon

!================================================================

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
