module polarstereo

IMPLICIT NONE

PRIVATE

PUBLIC :: stereo_loc_latlon_to_xy, stereo_loc_xy_to_latlon, stereo_vec_lonlat_to_xy

REAL*8,PARAMETER :: pi      = 3.1415926535898,        &
&                   deg2rad = 0.01745329251994333333, & ! pi/180
&                   rad2deg = 57.29577951308219756114   ! 180/pi

CONTAINS

!=========================================================================================
! N. Jourdain, LGGE, Grenoble, France
! Dec. 2014
!
! Contains :
! ~~~~~~~~   - stereo_loc_latlon_to_xy : location from (lat,lon) to (x,y)
!            - stereo_loc_xy_to_latlon : location from (x,y) to (lat,lon)
!            - stereo_vec_lonlat_to_xy : vector from (lat,lon) proj to (x,y) proj
!
! some standard info :
! ~~~~~~~~~~~~~~~~~~  
!  WGS84 - radius: 6378137.0 eccentricity: 0.08181919
!     in Matlab: axes2ecc(6378137.0, 6356752.3142)
!  Hughes ellipsoid - radius: 6378.273 km eccentricity: 0.081816153
!     Used for SSM/I  http://nsidc.org/data/polar_stereo/ps_grids.html
!  International ellipsoid (following Snyder) - radius: 6378388.0 eccentricity: 0.0819919
!
! history : two first subroutines originally written by Andy Bliss as matlab scripts
! ~~~~~~~
!
! Equations from: Map Projections - A Working manual - by J.P. Snyder. 1987 
! http://kartoweb.itc.nl/geometrics/Publications/Map%20Projections%20-%20A%20Working%20manual%20-%20by%20J.P.%20Snyder.pdf
! See the section on Polar Stereographic, with a south polar aspect and known phi_c not at the pole.
!
!=========================================================================================

SUBROUTINE stereo_loc_latlon_to_xy(mx,my,lat_in,lon_in,x_out,y_out,RT,ex,lat_true,lon_posy)
!---------------------------------------------------------------------------------------------
! Writen from matlab's function [x,y]=polarstereo_fwd(phi,lambda,a,e,phi_c,lambda_0)
! Purpose: transforms lat/lon data to map coordinates for a polar stereographic system
!---------------------------------------------------------------------------------------------
INTEGER,                  INTENT(IN)  :: mx, my            !- matrix dimensions
REAL*8, DIMENSION(mx,my), INTENT(IN)  :: lat_in, lon_in    !- input longitude and latitude (degrees)
REAL*8, DIMENSION(mx,my), INTENT(OUT) :: x_out, y_out      !- output x,y coordinates
REAL*8, OPTIONAL,         INTENT(IN)  :: RT,             & !- Earth radius (m)                     [default = WGS84]
&                                        ex,             & !- Earth's misshapenness (excentricity) [default = WGS84]
&                                        lat_true,       & !- latitude of true scale in degrees    [default = 71 degS]
&                                        lon_posy          !- meridian in degrees along the positive Y axis of the map [default = 0 deg]
REAL*8, DIMENSION(mx,my)              :: rho, t, lat, lon
REAL*8                                :: t_c, m_c, a, e, lat_c, lon_0
INTEGER                               :: pm

!- Optional arguments :
if ( present(RT)       ) then; a    = RT       ; else; a     = 6378137.0        ; endif
if ( present(ex)       ) then; e    = ex       ; else; e     =       0.08181919 ; endif
if ( present(lat_true) ) then; lat_c= lat_true ; else; lat_c =     -71.0        ; endif
if ( present(lon_posy) ) then; lon_0= lon_posy ; else; lon_0 =       0.0        ; endif

!- input checking :
if ( a .lt. 6.e6 .or. a .gt. 7.e6 ) then
  write(*,*) '~!@#$%^* error : Earth radius is wrong >>>>>> stop !!'
  stop
endif
!-
if ( e .lt. 0.0 .or. e .gt. 1.0 ) then
  write(*,*) '~!@#$%^* error : Earth excentricity is wrong >>>>>> stop !!' 
  stop
endif
!-
if ( lat_c .lt. -90.0 .or. lat_c .gt. 90.0 ) then
  write(*,*) '~!@#$%^* error : latitude of true scale is wrong >>>>>> stop !!'
  stop
endif
!-
if ( lon_0 .lt. -180.0 .or. lon_0 .gt. 360.0 ) then
  write(*,*) '~!@#$%^* error : meridian of positive Y axis is wrong >>>>>> stop !!'
  stop
endif

!- convert to radians :
lat(:,:) = deg2rad * lat_in(:,:)
lat_c    = deg2rad * lat_c
lon(:,:) = deg2rad * lon_in(:,:)
lon_0    = deg2rad * lon_0

!- if the standard parallel is in S.Hemi., switch signs.
if ( lat_c .lt. 0.0 ) then 
    pm       = -1    ! plus or minus, north lat. or south
    lat(:,:) = -lat(:,:)
    lat_c    = -lat_c
    lon(:,:) = -lon(:,:)
    lon_0    = -lon_0
else
    pm       = 1
endif

!- See Snyder for details.
t(:,:)   = tan(pi/4-lat(:,:)/2) / ( (1-e*sin(lat(:,:))) / (1+e*sin(lat(:,:))) )**(e/2)
t_c      = tan(pi/4-lat_c   /2) / ( (1-e*sin(lat_c   )) / (1+e*sin(lat_c   )) )**(e/2)
m_c      = cos(lat_c) / sqrt( 1 - e**2 * (sin(lat_c))**2 )
rho(:,:) = a * m_c * t(:,:) / t_c  !- true scale at lat lat_c

x_out(:,:) =  pm * rho(:,:) * sin(lon(:,:)-lon_0)
y_out(:,:) = -pm * rho(:,:) * cos(lon(:,:)-lon_0)

END SUBROUTINE stereo_loc_latlon_to_xy


SUBROUTINE stereo_loc_xy_to_latlon(mx,my,x_in,y_in,lat_out,lon_out,RT,ex,lat_true,lon_posy)
!---------------------------------------------------------------------------------------------
! Writen from matlab's function [phi,lambda]=polarstereo_inv(x,y,a,e,phi_c,lambda_0)
! Purpose: transforms lat/lon data to map coordinates for a polar stereographic system
!---------------------------------------------------------------------------------------------
INTEGER,                  INTENT(IN)  :: mx, my            !- matrix dimensions
REAL*8, DIMENSION(mx,my), INTENT(IN)  :: x_in, y_in        !- input x,y coordinates
REAL*8, DIMENSION(mx,my), INTENT(OUT) :: lat_out, lon_out  !- outut longitude and latitude (degrees)
REAL*8, OPTIONAL,         INTENT(IN)  :: RT,             & !- Earth radius (m)                     [default = WGS84]
&                                        ex,             & !- Earth's misshapenness (excentricity) [default = WGS84]
&                                        lat_true,       & !- latitude of true scale in degrees    [default = 71 degS]
&                                        lon_posy          !- meridian in degrees along the positive Y axis of the map [default = 0 deg]
REAL*8, DIMENSION(mx,my)              :: rho, t, chi, x, y
REAL*8                                :: t_c, m_c, a, e, lat_c, lon_0
INTEGER                               :: pm

!- Optional arguments :
if ( present(RT)       ) then; a    = RT       ; else; a     = 6378137.0        ; endif
if ( present(ex)       ) then; e    = ex       ; else; e     =       0.08181919 ; endif
if ( present(lat_true) ) then; lat_c= lat_true ; else; lat_c =     -71.0        ; endif
if ( present(lon_posy) ) then; lon_0= lon_posy ; else; lon_0 =       0.0        ; endif

!- input checking :
if ( a .lt. 6.e6 .or. a .gt. 7.e6 ) then
  write(*,*) '~!@#$%^* error : Earth radius is wrong >>>>>> stop !!'
  stop
endif
!-
if ( e .lt. 0.0 .or. e .gt. 1.0 ) then
  write(*,*) '~!@#$%^* error : Earth excentricity is wrong >>>>>> stop !!' 
  stop
endif
!-
if ( lat_c .lt. -90.0 .or. lat_c .gt. 90.0 ) then
  write(*,*) '~!@#$%^* error : latitude of true scale is wrong >>>>>> stop !!'
  stop
endif
!-
if ( lon_0 .lt. -180.0 .or. lon_0 .gt. 360.0 ) then
  write(*,*) '~!@#$%^* error : meridian of positive Y axis is wrong >>>>>> stop !!'
  stop
endif

!- convert to radians :
lat_c    = deg2rad * lat_c
lon_0    = deg2rad * lon_0

!- if the standard parallel is in S.Hemi., switch signs.
if ( lat_c .lt. 0.0 ) then 
    pm     = -1    ! plus or minus, north lat. or south
    lat_c  = -lat_c
    lon_0  = -lon_0
    x(:,:) = -x_in(:,:)
    y(:,:) = -y_in(:,:)
else
    pm     = 1
endif

!- See Snyder for details.
t_c      = tan(pi/4-lat_c   /2) / ( (1-e*sin(lat_c   )) / (1+e*sin(lat_c   )) )**(e/2)
m_c      = cos(lat_c) / sqrt( 1 - e**2 * (sin(lat_c))**2 )
rho(:,:) = sqrt(x(:,:)**2+y(:,:)**2)
t(:,:)   = rho(:,:) * t_c / ( a * m_c )

!- iterate to find phi, with a threshold of pi*1e-8
!  phi_alt=pi/2 - 2 * atan(t); %guess for phi
!  phiold=phi_alt+10; %+10 to make sure it executes the while loop
!  while any(abs(phi_alt(:)-phiold(:)) > pi*1e-8)
!     phiold=phi_alt;
!     phi_alt=pi/2 - 2*atan(t.*((1-e*sin(phi_alt))./(1+e*sin(phi_alt))).^(e/2));
!     %add a break in case it doesn't converge?
!  end

chi(:,:) = 0.5*pi - 2 * atan(t(:,:))  !- find lat with a series instead of iterating.

lat_out(:,:) = chi(:,:) + ( (1./2.) * e**2 + (5./24.) * e**4 + ( 1./ 12.) * e**6 + (  13./   360.) * e**8 ) * sin(2*chi(:,:)) &
&                       + (                  (7./48.) * e**4 + (29./240.) * e**6 + ( 811./ 11520.) * e**8 ) * sin(4*chi(:,:)) &
&                       + (                                  + ( 7./120.) * e**6 + (  81./  1120.) * e**8 ) * sin(6*chi(:,:)) &
&                       + (                                                      + (4279./161280.) * e**8 ) * sin(8*chi(:,:))
lon_out(:,:) = lon_0 + atan2(x,-y)

!- correct the signs and phasing :
lat_out(:,:) = pm * lat_out(:,:)
lon_out(:,:) = pm * lon_out(:,:)
lon_out(:,:) = mod(lon_out(:,:)+pi,2*pi)-pi !- want longitude in the range -pi to pi

!- convert back to degrees :
lat_out(:,:) = rad2deg * lat_out(:,:)
lon_out(:,:) = rad2deg * lon_out(:,:)

END SUBROUTINE stereo_loc_xy_to_latlon


SUBROUTINE stereo_vec_lonlat_to_xy(mx,my,lon_in,uzo_in,vme_in,ux_out,vy_out,lon_posy)
!---------------------------------------------------------------------------------------------
! Purpose: Express (zonal,meridional) vector componants as (x,y) componants on the stereographic map
!---------------------------------------------------------------------------------------------
INTEGER,                  INTENT(IN)  :: mx, my            !- matrix dimensions
REAL*8, DIMENSION(mx,my), INTENT(IN)  :: lon_in            !- input longitude (deg East) 
REAL*8, DIMENSION(mx,my), INTENT(IN)  :: uzo_in, vme_in    !- vector on input lon,lat grid
REAL*8, DIMENSION(mx,my), INTENT(OUT) :: ux_out, vy_out    !- vector on output x,y grid
REAL*8, OPTIONAL,         INTENT(IN)  :: lon_posy          !- meridian in degrees along the positive Y axis of the map [default = 0 deg]
REAL*8, DIMENSION(mx,my)              :: lon
REAL*8                                :: lon_0

!- Optional arguments :
if ( present(lon_posy) ) then; lon_0= lon_posy ; else; lon_0 =       0.0        ; endif

!- input checking :
if ( lon_0 .lt. -180.0 .or. lon_0 .gt. 360.0 ) then
  write(*,*) '~!@#$%^* error : meridian of positive Y axis is wrong >>>>>> stop !!'
  stop
endif

!- convert to radians :
lon(:,:) = deg2rad * lon_in(:,:)
lon_0    = deg2rad * lon_0

ux_out(:,:) =   uzo_in(:,:) * cos(lon(:,:)-lon_0) + vme_in(:,:) * sin(lon(:,:)-lon_0)
vy_out(:,:) = - uzo_in(:,:) * sin(lon(:,:)-lon_0) + vme_in(:,:) * cos(lon(:,:)-lon_0)

END SUBROUTINE stereo_vec_lonlat_to_xy

!=======================================================
end module polarstereo
