program plotgrids

   use map_utils
 
   implicit none

   ! Parameters
   integer, parameter :: MAX_DOMAINS = 21

   ! Variables
   integer :: iproj_type, n_domains, io_form_output, dyn_opt
   integer :: i, j, max_dom, funit, io_form_geogrid
   integer :: interval_seconds

   integer, dimension(MAX_DOMAINS) :: parent_grid_ratio, parent_id, ixdim, jydim
   integer, dimension(MAX_DOMAINS) :: i_parent_start, j_parent_start, &
                        s_we, e_we, s_sn, e_sn, &
                        start_year, start_month, start_day, start_hour, &
                        end_year,   end_month,   end_day,   end_hour

   real :: known_lat, known_lon, stand_lon, truelat1, truelat2, known_x, known_y, &
           dxkm, dykm, phi, lambda, ref_lat, ref_lon, ref_x, ref_y
   real :: dx, dy
   real :: ri, rj, rlats, rlons, rlate, rlone
   real :: polat , rot
   real :: rparent_gridpts
   real :: xa,xb,ya,yb,xxa,xxy,yya,yyb
   real :: xs, xe, ys, ye
   integer :: jproj, jgrid, jlts, iusout, idot, ier
   integer :: ltype , idom

   real, dimension(MAX_DOMAINS) :: parent_ll_x, parent_ll_y, parent_ur_x, parent_ur_y, &
                                   dom_cen_lat, dom_cen_lon

   character (len=128) :: geog_data_path, opt_output_from_geogrid_path, opt_geogrid_tbl_path
   character (len=128), dimension(MAX_DOMAINS) :: geog_data_res 
   character (len=128) :: map_proj,ncep_proc_path, ncep_proc_prefix, &
                          ncep_proc_domain_type
   character (len=128), dimension(MAX_DOMAINS) :: start_date, end_date
   character (len=3) :: wrf_core
   character (len=1) :: gridtype

   logical :: do_tiled_output
   integer :: debug_level
   logical :: is_used
   logical :: ncep_processing, do_gwd, just_last, use_igbp, ncep_proc_grib2

   type (proj_info) :: map_projection

   namelist /share/ wrf_core, max_dom, start_date, end_date, &
                     start_year, end_year, start_month, end_month, &
                     start_day, end_day, start_hour, end_hour, &
                     interval_seconds, &
                     io_form_geogrid, opt_output_from_geogrid_path, debug_level
   namelist /geogrid/ parent_id, parent_grid_ratio, &
                      i_parent_start, j_parent_start, s_we, e_we, s_sn, e_sn, &
                      map_proj, ref_x, ref_y, ref_lat, ref_lon, dom_cen_lat, &
                      dom_cen_lon, &
                      truelat1, truelat2, stand_lon, dx, dy, &
                      geog_data_res, geog_data_path, opt_geogrid_tbl_path, &
                      ncep_processing, ncep_proc_path, ncep_proc_prefix, &
                      ncep_proc_domain_type, do_gwd, just_last, use_igbp, &
                      ncep_proc_grib2

  
   ! Set defaults for namelist variables
   debug_level = 0
   io_form_geogrid = 2
   wrf_core = 'ARW'
   max_dom = 1
   geog_data_path = 'NOT_SPECIFIED'
   ref_x = NAN
   ref_y = NAN
   ref_lat = NAN
   ref_lon = NAN
   dx = 10000.
   dy = 10000.
   map_proj = 'Lambert'
   truelat1 = NAN
   truelat2 = NAN
   stand_lon = NAN
   do i=1,MAX_DOMAINS
      geog_data_res(i) = 'default'
      parent_id(i) = 1
      parent_grid_ratio(i) = INVALID
      s_we(i) = 1
      e_we(i) = INVALID
      s_sn(i) = 1
      e_sn(i) = INVALID
      start_year(i) = 0
      start_month(i) = 0
      start_day(i) = 0
      start_hour(i) = 0
      end_year(i) = 0
      end_month(i) = 0
      end_day(i) = 0
      end_hour(i) = 0
      start_date(i) = '0000-00-00_00:00:00'
      end_date(i) = '0000-00-00_00:00:00'
   end do
   opt_output_from_geogrid_path = './'
   opt_geogrid_tbl_path = 'geogrid/'
   interval_seconds = INVALID
   
   ! Read parameters from Fortran namelist
   do funit=10,100
      inquire(unit=funit, opened=is_used)
      if (.not. is_used) exit
   end do
   open(funit,file='namelist.nps',status='old',form='formatted',err=1000)
   read(funit,share)
   read(funit,geogrid)
   close(funit)

   dxkm = dx
   dykm = dy

   known_lat = ref_lat
   known_lon = ref_lon
   known_x = ref_x
   known_y = ref_y

   ! Convert wrf_core to uppercase letters
   do i=1,3
      if (ichar(wrf_core(i:i)) >= 97) wrf_core(i:i) = char(ichar(wrf_core(i:i))-32)
   end do

   ! Before doing anything else, we must have a valid grid type 
   gridtype = ' '
   if (wrf_core == 'ARW') then
      gridtype = 'C'
      dyn_opt = 2
   else if (wrf_core == 'NMM') then
      gridtype = 'E'
      dyn_opt = 4
      CALL plot_e_grid ( e_we(1)-1 , e_sn(1)-1 , ref_lat , -1. * ref_lon , dy , dx ) 
      stop
   else if (wrf_core == 'NMB') then
      gridtype = 'B'
      dyn_opt = 4
      CALL plot_b_grid ( e_we , e_sn , ref_lat , -1. * ref_lon , dy , dx , &
                        debug_level,dom_cen_lat, dom_cen_lon, parent_grid_ratio, &
                        parent_id, max_dom   ) 
      stop
   end if

   if (gridtype /= 'C' .and. gridtype /= 'E') then
      write(6,*) 'A valid wrf_core must be specified in the namelist. '// &
                 'Currently, only "ARW" and "NMM" are supported.'
      stop
   end if

   if (max_dom > MAX_DOMAINS) then
      write(6,*) 'In namelist, max_dom must be <= ',MAX_DOMAINS,'. To run with more'// &
                ' than ',MAX_DOMAINS,' domains, increase the MAX_DOMAINS parameter.'
      stop
   end if

   ! Every domain must have a valid parent id
   do i=2,max_dom
      if (parent_id(i) <= 0 .or. parent_id(i) >= i) then
         write(6,*) 'In namelist, the parent_id of domain ',i,' must be in '// &
                   'the range 1 to ',i-1
          stop
      end if
   end do

   ! Convert map_proj to uppercase letters
   do i=1,len(map_proj)
      if (ichar(map_proj(i:i)) >= 97) map_proj(i:i) = char(ichar(map_proj(i:i))-32)
   end do

   ! Assign parameters to module variables
   if ((index(map_proj, 'LAMBERT') /= 0) .and. &
       (len_trim(map_proj) == len('LAMBERT'))) then
      iproj_type = PROJ_LC 
      rot=truelat1
      polat=truelat2
      jproj = 3

   else if ((index(map_proj, 'MERCATOR') /= 0) .and. &
            (len_trim(map_proj) == len('MERCATOR'))) then
      iproj_type = PROJ_MERC 
      rot=0.
      polat=0.
      jproj = 9

   else if ((index(map_proj, 'POLAR') /= 0) .and. &
            (len_trim(map_proj) == len('POLAR'))) then
      iproj_type = PROJ_PS 
      rot=0.
      polat=SIGN(90., ref_lat)
      jproj = 1

   else if ((index(map_proj, 'ROTATED_LL') /= 0) .and. &
            (len_trim(map_proj) == len('ROTATED_LL'))) then
      iproj_type = PROJ_ROTLL 

   else
         write(6,*) 'In namelist, invalid map_proj specified. Valid '// &
                    'projections are "lambert", "mercator", "polar", '// &
                    'and "rotated_ll".'
   end if

   n_domains = max_dom

   do i=1,n_domains
      ixdim(i) = e_we(i) - s_we(i) + 1
      jydim(i) = e_sn(i) - s_sn(i) + 1
   end do

   if (gridtype == 'E') then
      phi = dykm*real(jydim(1)-1)/2.
      lambda = dxkm*real(ixdim(1)-1)
   end if

   ! If the user hasn't supplied a known_x and known_y, assume the center of domain 1
   if (known_x == NAN) known_x = ixdim(1) / 2.
   if (known_y == NAN) known_y = jydim(1) / 2.

   ! Checks specific to C grid
   if (gridtype == 'C') then

      ! C grid does not support the rotated lat/lon projection
      if (iproj_type == PROJ_ROTLL) then
         write(6,*) 'Rotated lat/lon projection is not supported for the ARW core. '// &
                    'Valid projecitons are "lambert", "mercator", and "polar".'
         stop
      end if

      ! Check that nests have an acceptable number of grid points in each dimension
      do i=2,n_domains
         rparent_gridpts = real(ixdim(i)-1)/real(parent_grid_ratio(i))
         if (floor(rparent_gridpts) /= ceiling(rparent_gridpts)) then
            write(6,*) 'For nest ',i,' (e_we-s_we+1) must be one greater than an '// &
                       'integer multiple of the parent_grid_ratio.'
            stop
         end if
         rparent_gridpts = real(jydim(i)-1)/real(parent_grid_ratio(i))
         if (floor(rparent_gridpts) /= ceiling(rparent_gridpts)) then
            write(6,*) 'For nest ',i,' (e_sn-s_sn+1) must be one greater than an '// &
                       'integer multiple of the parent_grid_ratio.'
            stop
         end if
      end do
   end if

   do i=1,n_domains
      parent_ll_x(i) = real(i_parent_start(i))
      parent_ll_y(i) = real(j_parent_start(i))
      parent_ur_x(i) = real(i_parent_start(i))+real(ixdim(i))/real(parent_grid_ratio(i))-1.
      parent_ur_y(i) = real(j_parent_start(i))+real(jydim(i))/real(parent_grid_ratio(i))-1.
   end do

   call map_init(map_projection)

   call map_set(iproj_type, map_projection, &
                lat1=known_lat, &
                lon1=known_lon, &
                knowni=known_x, &
                knownj=known_y, &
                dx=dx, &
                stdlon=stand_lon, &
                truelat1=truelat1, &
                truelat2=truelat2, &
                ixdim=ixdim(1), &
                jydim=jydim(1))

   call ij_to_latlon(map_projection, 1., 1., rlats, rlons)
   call ij_to_latlon(map_projection, real(e_we(1)), real(e_sn(1)) , rlate, rlone)

   call opngks

   ! Set some colors
   call gscr(1, 0, 1.00, 1.00, 1.00)
   call gscr(1, 1, 0.00, 0.00, 0.00)

   ! Do not grind them with details
   jgrid=10
   jlts=-2
   iusout=1
   idot=0

   call supmap(jproj,polat,stand_lon,rot,&
               rlats,rlons,rlate,rlone, &
               jlts,jgrid,iusout,idot,ier) 

   call setusv('LW',1000)
   call perim(e_we(1)-1,1,e_sn(1)-1,1)
   call getset(xa,xb,ya,yb,xxa,xxy,yya,yyb,ltype)
   call set   (xa,xb,ya,yb, &
         1.,real(e_we(1)),1.,real(e_sn(1)),ltype)

   do idom = 2 , max_dom
      call getxy ( xs, xe, ys, ye, &
                   idom , max_dom , &
                   e_we , e_sn , &
                   parent_id , parent_grid_ratio , &
                   i_parent_start , j_parent_start )
      call line ( xs , ys , xe , ys )
      call line ( xe , ys , xe , ye )
      call line ( xe , ye , xs , ye )
      call line ( xs , ye , xs , ys )
   end do

   call frame

   call clsgks

   stop

1000 write(6,*) 'Error opening namelist.wps'
   stop
  
end program plotgrids

subroutine getxy ( xs, xe, ys, ye, &
                   dom_id , num_domains , &
                   e_we , e_sn , &
                   parent_id , parent_grid_ratio , &
                   i_parent_start , j_parent_start )

   implicit none

   integer , intent(in) :: dom_id
   integer , intent(in) :: num_domains
   integer , intent(in) , dimension(num_domains):: e_we , e_sn , &
                                                   parent_id , parent_grid_ratio , &
                                                   i_parent_start , j_parent_start
   real , intent(out) :: xs, xe, ys, ye


   !  local vars

   integer :: idom

   xs = 0.
   xe = e_we(dom_id) -1
   ys = 0.
   ye = e_sn(dom_id) -1

   idom = dom_id
   compute_xy : DO

      xs = (i_parent_start(idom) + xs  -1 ) / &    
           real(parent_grid_ratio(parent_id(idom)))
      xe = xe / real(parent_grid_ratio(idom))

      ys = (j_parent_start(idom) + ys  -1 ) / &    
           real(parent_grid_ratio(parent_id(idom)))
      ye = ye / real(parent_grid_ratio(idom))

      idom = parent_id(idom)
      if ( idom .EQ. 1 ) then
         exit compute_xy
      end if

   END DO compute_xy

   xs = xs + 1
   xe = xs + xe
   ys = ys + 1
   ye = ys + ye

end subroutine getxy

SUBROUTINE plot_b_grid ( im , jm , rlat0d , rlon0d , dphd , dlmd , debug_level, &
                         dom_cen_lat, dom_cen_lon, parent_grid_ratio, parent_id, max_dom   ) 

!  This routine generates a gmeta file of the area covered by an Arakawa b-grid.
!  We assume that NCAR Graphics has not been called yet (and will be closed
!  upon exit).  The required input fields are as from the WPS namelist file.

   IMPLICIT NONE

   !  Input map parameters for E/B grid.

   INTEGER , INTENT(IN) :: im(21) , &     ! number of H points in odd rows
                           jm(21) , &     ! number of rows
                           max_dom , &        !  total number of domains
                           debug_level         !  debug

   REAL , INTENT(IN)    :: rlat0d , & ! latitude of grid center (degrees)
                           rlon0d     ! longitude of grid center (degrees, times -1)

   REAL , INTENT(IN)    :: dphd , &   ! angular distance between rows (degrees)
                           dlmd       ! angular distance between adjacent H and V points (degrees)

   REAL, INTENT(IN) :: dom_cen_lat(21),dom_cen_lon(21)
   INTEGER, INTENT(IN) :: parent_grid_ratio(21),parent_id(21)

   !  Some local vars

   REAL :: rlat1d , &
           rlon1d

   INTEGER :: ngpwe , &
              ngpsn , &
              ilowl , &
              jlowl
     
   INTEGER :: imt , imtjm
   REAL :: latlft,lonlft,latrgt,lonrgt

   imt=2*im(1)-1
   imtjm=imt*jm(1)

   rlat1d=rlat0d
   rlon1d=rlon0d
   ngpwe=im(1)
   ngpsn=jm(1)
!        write(0,*) 'ngpwe , ngpsn: ' , ngpwe , ngpsn

   !  Get lat and lon of left and right points.

   CALL corners_bgrid ( rlat1d,rlon1d,im(1),jm(1),rlat0d,rlon0d,dphd,dlmd,&
                  ngpwe,ngpsn,ilowl,jlowl,latlft,lonlft,latrgt,lonrgt)

   !  With corner points, make map background.

   CALL mapbkg_bgrid ( im,jm,ilowl,jlowl,ngpwe,ngpsn,&
                       rlat0d,rlon0d,latlft,lonlft,latrgt,lonrgt,&
                       dlmd,dphd,debug_level, &
                       dom_cen_lat, dom_cen_lon, parent_grid_ratio, parent_id, max_dom ) 

END SUBROUTINE plot_b_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  E GRID MAP INFO BELOW     !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE plot_e_grid ( im , jm , rlat0d , rlon0d , dphd , dlmd ) 

!  This routine generates a gmeta file of the area covered by an Arakawa e-grid.
!  We assume that NCAR Graphics has not been called yet (and will be closed
!  upon exit).  The required input fields are as from the WPS namelist file.

   IMPLICIT NONE

!  15 April 2005 NCEP/EMC 
!    The Code and some instructions are provided by Tom BLACK to
!    NCAR/DTC Meral Demirtas			

!  4 May 2005  NCAR/DTC Meral DEMIRTAS  
!    - An include file (plot_inc) is added to get
!    Domain size: IM,JM
!    Central latitute and longnitute: RLAT0D,RLON0D
!    Horizontal resolution: DPHD, DLMD

!  Feb 2007 NCAR/MMM
!    Turn into f90
!    Add implicit none
!    Remove non-mapping portions
!    Make part of WPS domain plotting utility


   !  Input map parameters for E grid.

   INTEGER , INTENT(IN) :: im , &     ! number of H points in odd rows
                           jm         ! number of rows

   REAL , INTENT(IN)    :: rlat0d , & ! latitude of grid center (degrees)
                           rlon0d     ! longitude of grid center (degrees, times -1)

   REAL , INTENT(IN)    :: dphd , &   ! angular distance between rows (degrees)
                           dlmd       ! angular distance between adjacent H and V points (degrees)

   !  Some local vars

   REAL :: rlat1d , &
           rlon1d

   INTEGER :: ngpwe , &
              ngpsn , &
              ilowl , &
              jlowl
     
   INTEGER :: imt , imtjm
   REAL :: latlft,lonlft,latrgt,lonrgt

   imt=2*im-1
   imtjm=imt*jm
   rlat1d=rlat0d
   rlon1d=rlon0d
   ngpwe=2*im-1
   ngpsn=jm

   !  Get lat and lon of left and right points.

   CALL corners ( rlat1d,rlon1d,im,jm,rlat0d,rlon0d,dphd,dlmd,&
                  ngpwe,ngpsn,ilowl,jlowl,latlft,lonlft,latrgt,lonrgt)

   !  With corner points, make map background.

   CALL mapbkg_egrid ( imt,jm,ilowl,jlowl,ngpwe,ngpsn,&
                       rlat0d,rlon0d,latlft,lonlft,latrgt,lonrgt,&
                       dlmd,dphd)

END SUBROUTINE plot_e_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE corners_bgrid ( glatd,glond,im,jm,tph0d,tlm0d,dphd,dlmd,&
                     ngpwe,ngpsn,ilowl,jlowl,glatl,glonl,glatr,glonr)

   IMPLICIT NONE

   REAL , INTENT(IN) :: glatd,glond,tph0d,tlm0d,dphd,dlmd,&
                        glatl,glonl,glatr,glonr

   INTEGER , INTENT(IN) :: im,jm,ngpwe,ngpsn
   INTEGER , INTENT(OUT) :: ilowl,jlowl

   !  Local vars

   REAL , PARAMETER :: d2r = 1.74532925E-2 , r2D = 1./D2R

   REAL :: glat , glon , dph , dlm , tph0 , tlm0
   REAL :: x , y , z , tlat , tlon , tlat1 , tlat2 , tlon1 , tlon2
   REAL :: row , col
   REAL :: dlm1 , dlm2 , d1 , d2 , d3 , d4 , dmin

   INTEGER :: jmt , ii , jj , iuppr , juppr, imt
   INTEGER :: nrow , ncol

   jmt = jm/2+1
   imt = im/2+1

   !  Convert from geodetic to transformed coordinates (degrees).

   glat = glatd * d2r
   glon = glond * d2r
   dph = dphd * d2r
   dlm = dlmd * d2r
   tph0 = tph0d * d2r
   tlm0 = tlm0d * d2r

   x = COS(tph0) * COS(glat) * COS(glon-tlm0)+SIN(tph0) * SIN(glat)
   y = -COS(glat) * SIN(glon-tlm0)
   z = COS(tph0) * SIN(glat)-SIN(tph0) * COS(glat) * COS(glon-tlm0)
   tlat = r2d * ATAN(z/SQRT(x*x + y*y))
   tlon = r2d * ATAN(y/x)

   !  Find the real (non-integer) row and column of the input location on 
   !  the filled e-grid.

!        write(0,*) 'tlat, tlon: ', tlat, tlon

   row = tlat/dphd+jmt
   col = tlon/dlmd+imt
   nrow = INT(row)
   ncol = INT(col)
!        write(0,*) 'ncol, nrow: ', ncol, nrow

   tlat = tlat * d2r
   tlon = tlon * d2r

!               E2     E3
! 
! 
!                  X
!               E1     E4

   tlat1 = (nrow-jmt) * dph
   tlat2 = tlat1+dph
   tlon1 = (ncol-imt) * dlm
   tlon2 = tlon1+dlm

!        write(0,*) 'tlat1, tlat2: ', tlat1, tlat2
!        write(0,*) 'tlon1, tlon2: ', tlon1, tlon2

   dlm1 = tlon-tlon1
   dlm2 = tlon-tlon2

   d1 = ACOS(COS(tlat) * COS(tlat1) * COS(dlm1)+SIN(tlat) * SIN(tlat1))
   d2 = ACOS(COS(tlat) * COS(tlat2) * COS(dlm1)+SIN(tlat) * SIN(tlat2))
   d3 = ACOS(COS(tlat) * COS(tlat2) * COS(dlm2)+SIN(tlat) * SIN(tlat2))
   d4 = ACOS(COS(tlat) * COS(tlat1) * COS(dlm2)+SIN(tlat) * SIN(tlat1))

   dmin = MIN(d1,d2,d3,d4)

   IF      ( ABS(dmin-d1) .LT. 1.e-6 ) THEN
     ii = ncol
     jj = nrow
   ELSE IF ( ABS(dmin-d2) .LT. 1.e-6 ) THEN
     ii = ncol
     jj = nrow+1
   ELSE IF ( ABS(dmin-d3) .LT. 1.e-6 ) THEN
     ii = ncol+1
     jj = nrow+1
   ELSE IF ( ABS(dmin-d4) .LT. 1.e-6 ) THEN
     ii = ncol+1
     jj = nrow
   END IF

!      write(0,*) 'ii,jj final: ', ii, jj

   !  Now find the i and j of the lower left corner of the desired grid 
   !  region and of the upper right.

   ilowl = ii-ngpwe/2
   jlowl = jj-ngpsn/2
   iuppr = ii+ngpwe/2
   juppr = jj+ngpsn/2

!        write(0,*) 'ilowl, jlowl: ', ilowl, jlowl
!        write(0,*) 'iuppr, juppr: ', iuppr,juppr

        iuppr=im
        juppr=jm

   !  Find their geodetic coordinates.

   CALL b2t2g(ilowl,jlowl,im,jm,tph0d,tlm0d,dphd,dlmd,glatl,glonl)
!        write(0,*) 'glatl,glonl: ', glatl,glonl
   CALL b2t2g(iuppr,juppr,im,jm,tph0d,tlm0d,dphd,dlmd,glatr,glonr)
!        write(0,*) 'glatr,glonr: ', glatr,glonr

END SUBROUTINE corners_bgrid

SUBROUTINE corners ( glatd,glond,im,jm,tph0d,tlm0d,dphd,dlmd,&
                     ngpwe,ngpsn,ilowl,jlowl,glatl,glonl,glatr,glonr)

   IMPLICIT NONE

   REAL , INTENT(IN) :: glatd,glond,tph0d,tlm0d,dphd,dlmd,&
                        glatl,glonl,glatr,glonr

   INTEGER , INTENT(IN) :: im,jm,ngpwe,ngpsn
   INTEGER , INTENT(OUT) :: ilowl,jlowl

   !  Local vars

   REAL , PARAMETER :: d2r = 1.74532925E-2 , r2D = 1./D2R

   REAL :: glat , glon , dph , dlm , tph0 , tlm0
   REAL :: x , y , z , tlat , tlon , tlat1 , tlat2 , tlon1 , tlon2
   REAL :: row , col
   REAL :: dlm1 , dlm2 , d1 , d2 , d3 , d4 , dmin

   INTEGER :: jmt , ii , jj , iuppr , juppr
   INTEGER :: nrow , ncol

   jmt = jm/2+1

   !  Convert from geodetic to transformed coordinates (degrees).

   glat = glatd * d2r
   glon = glond * d2r
   dph = dphd * d2r
   dlm = dlmd * d2r
   tph0 = tph0d * d2r
   tlm0 = tlm0d * d2r

   x = COS(tph0) * COS(glat) * COS(glon-tlm0)+SIN(tph0) * SIN(glat)
   y = -COS(glat) * SIN(glon-tlm0)
   z = COS(tph0) * SIN(glat)-SIN(tph0) * COS(glat) * COS(glon-tlm0)
   tlat = r2d * ATAN(z/SQRT(x*x + y*y))
   tlon = r2d * ATAN(y/x)

   !  Find the real (non-integer) row and column of the input location on 
   !  the filled e-grid.

   row = tlat/dphd+jmt
   col = tlon/dlmd+im
   nrow = INT(row)
   ncol = INT(col)
   tlat = tlat * d2r
   tlon = tlon * d2r

!               E2     E3
! 
! 
!                  X
!               E1     E4

   tlat1 = (nrow-jmt) * dph
   tlat2 = tlat1+dph
   tlon1 = (ncol-im) * dlm
   tlon2 = tlon1+dlm

   dlm1 = tlon-tlon1
   dlm2 = tlon-tlon2

   d1 = ACOS(COS(tlat) * COS(tlat1) * COS(dlm1)+SIN(tlat) * SIN(tlat1))
   d2 = ACOS(COS(tlat) * COS(tlat2) * COS(dlm1)+SIN(tlat) * SIN(tlat2))
   d3 = ACOS(COS(tlat) * COS(tlat2) * COS(dlm2)+SIN(tlat) * SIN(tlat2))
   d4 = ACOS(COS(tlat) * COS(tlat1) * COS(dlm2)+SIN(tlat) * SIN(tlat1))

   dmin = MIN(d1,d2,d3,d4)

   IF      ( ABS(dmin-d1) .LT. 1.e-6 ) THEN
     ii = ncol
     jj = nrow
   ELSE IF ( ABS(dmin-d2) .LT. 1.e-6 ) THEN
     ii = ncol
     jj = nrow+1
   ELSE IF ( ABS(dmin-d3) .LT. 1.e-6 ) THEN
     ii = ncol+1
     jj = nrow+1
   ELSE IF ( ABS(dmin-d4) .LT. 1.e-6 ) THEN
     ii = ncol+1
     jj = nrow
   END IF

   !  Now find the i and j of the lower left corner of the desired grid 
   !  region and of the upper right.

   ilowl = ii-ngpwe/2
   jlowl = jj-ngpsn/2
   iuppr = ii+ngpwe/2
   juppr = jj+ngpsn/2

   !  Find their geodetic coordinates.

   CALL e2t2g(ilowl,jlowl,im,jm,tph0d,tlm0d,dphd,dlmd,glatl,glonl)
   CALL e2t2g(iuppr,juppr,im,jm,tph0d,tlm0d,dphd,dlmd,glatr,glonr)

END SUBROUTINE corners
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE mapbkg_bgrid ( im,jm,ilowl,jlowl,ngpwe,ngpsn,&
                          rlat0d,rlon0d,glatl,glonl,glatr,glonr,&
                          dlmd,dphd, debug_level, &
                          dom_cen_lat, dom_cen_lon, parent_grid_ratio, parent_id, max_dom ) 

!  IMPLICIT NONE

   !  Some local vars

   real, dimension(21) :: dom_cen_lat, dom_cen_lon, dx, dy
   integer, dimension(21) :: parent_grid_ratio, parent_id
   INTEGER :: max_dom, I, im(21),jm(21),debug_level
   CHARACTER (LEN=97) :: string
   REAL:: wbd_guess, sbd_guess, wbd_loc, sbd_loc
   REAL:: xs, xe, ys, ye
   REAL:: wbd(21),sbd(21)

   !  Yet more center lon messing around, hoo boy.

!  clonx=180.-rlon0d
   clonx=-rlon0d

   !  Open up NCAR Graphics

   CALL opngks
   CALL gopwk(8,9,3)
   CALL gsclip(0)

   !  Make the background white, and the foreground black.

   CALL gscr ( 1 , 0 , 1., 1., 1. )
   CALL gscr ( 1 , 1 , 0., 0., 0. )

   !  Line width twice default thickness.

   CALL setusv('LW',1800)

   !  Make map outline a solid line, not dots.

   CALL mapsti('MV',8)
   CALL mapsti('DO',0)

   !  Map outlines are political and states.

   CALL mapstc('OU','PS')

   !  Cylindrical equidistant.

   CALL maproj('CE',rlat0d,clonx,0.)

   !  Specify corner points.

   CALL mapset('CO',glatl,glonl,glatr,glonr)

   !  Lat lon lines every 5 degrees.

   CALL mapsti('GR',5)

   !  Map takes up this much real estate.

   CALL mappos( 0.03 , 0.97 , 0.03 , 0.97 )

   !  Initialize and draw map.

   CALL mapint
!   CALL mapdrw

   !  Add approx grid point tick marks

!   CALL setusv('LW',1000)
!   CALL perim(im,1,jm,1)


!!! these x y values are in degrees of latitude longitude

        if (debug_level > 0) then
        write(0,*) 'start new stuff'
        write(0,*) 'max_dom: ', max_dom
        endif

        dx(1)=dlmd
        dy(1)=dphd

        if (debug_level > 0) then
        write(0,*) 'dx(1), dy(1): ', dx(1), dy(1)
        write(0,*) 'im(1), jm(1): ', im(1), jm(1)
        endif

        wbd(1)=-((im(1)-1)/2.)*dlmd
        sbd(1)=-((jm(1)-1)/2.)*dphd

        if (debug_level > 0) then
        write(0,*) 'im(1), wbd(1): ', im(1), wbd(1) 
        write(0,*) 'jm(1), sbd(1): ', jm(1), sbd(1) 
        write(0,*) 'dom_cen_lat/lon(1): ', dom_cen_lat(1), dom_cen_lon(1)
        write(0,*) 'parent_id,ratio(i): ', parent_id(i), parent_grid_ratio(i)
        endif

        do I=2,max_dom

        if (debug_level > 0) then
        write(0,*) ' '
        write(0,*) '-----------------------------------'
        write(0,*) ' '
        write(0,*) 'dom_cen_lat(i): ', dom_cen_lat(i)
        write(0,*) 'dom_cen_lon(i): ', dom_cen_lon(i)
        endif

        call tll (dom_cen_lon(i),dom_cen_lat(i),xcenter, ycenter, dom_cen_lat(1),dom_cen_lon(1))             

        if (debug_level > 0) then
        write(0,*) 'im,jm(i): ', im(i),jm(i)
        write(0,*) 'xcenter, ycenter: ', xcenter, ycenter
!        write(0,*) 'dlmd, dphd: ', dlmd, dphd
        endif

         wbd_guess = xcenter - &
                     (im(i)-1)*0.5*dx(parent_id(i))/parent_grid_ratio(i)
         sbd_guess = ycenter - &
                     (jm(i)-1)*0.5*dy(parent_id(i))/parent_grid_ratio(i)

        if (debug_level > 0) then
        write(0,*) 'wbd_guess, sbd_guess: ', wbd_guess, sbd_guess

        write(0,*) 'i_parent_start parens: ', &
                   ( wbd_guess - wbd(parent_id(i)))/dx(parent_id(i)) 
        write(0,*) 'j_parent_start parens: ', &
                   ( sbd_guess - sbd(parent_id(i)))/dy(parent_id(i))
        endif

            i_parent_start = nint (1 + ( wbd_guess - wbd(parent_id(i)))/dx(parent_id(i)) )                      
            j_parent_start = nint (1 + ( sbd_guess - sbd(parent_id(i)))/dy(parent_id(i)) )                      

        if (debug_level > 0) then
        write(0,*) 'wbd(parent_id(i)), sbd(parent_id(i)): ', &
                    wbd(parent_id(i)), sbd(parent_id(i))

        write(0,*) 'dx,dy(parent_id(i)): ', dx(parent_id(i)), dy(parent_id(i))
        write(0,*) 'i/j_parent_start: ', i_parent_start, j_parent_start
        endif

            i_parent_end   = i_parent_start + (im(i)-1)/parent_grid_ratio(i)
            j_parent_end   = j_parent_start + (jm(i)-1)/parent_grid_ratio(i)

        if (debug_level > 0) then
        write(0,*) 'i/j_parent_end: ', i_parent_end, j_parent_end
        endif

         wbd(i)= wbd(parent_id(i)) + (i_parent_start-1)*dx(parent_id(i))
         sbd(i)= sbd(parent_id(i)) + (j_parent_start-1)*dy(parent_id(i))

        dx(i)=dx(parent_id(i))/parent_grid_ratio(i)
        dy(i)=dy(parent_id(i))/parent_grid_ratio(i)

        if (debug_level > 0) then
        write(0,*) 'dx(i), dy(i): ', dx(i), dy(i)
        endif

        wbd_loc=-((im(i)-1)/2)*dx(i)
        sbd_loc=-((jm(i)-1)/2)*dy(i)

        if (debug_level > 0) then
        write(0,*) 'wbd_loc, sbd_loc: ', wbd_loc, sbd_loc
        endif

        xs=wbd(i)
        xe=wbd(i)-2*wbd_loc
        ys=sbd(i)
        ye=sbd(i)-2*sbd_loc

        if (debug_level > 0) then
        write(0,*) 'xs, xe: ', wbd(i), wbd(i)-2*wbd_loc
        write(0,*) 'ys, ye: ', sbd(i), sbd(i)-2*sbd_loc
        write(0,*) '+ + + + + + + +'

        endif

      call line( xs , ys , xe , ys )
      call line( xe , ys , xe , ye )
      call line( xe , ye , xs , ye )
      call line( xs , ye , xs , ys )

        enddo

   CALL mapdrw

   !  Put on a quicky description.

   WRITE ( string , FMT = '("B-GRID E_WE = ",I4,", E_SN = ",I4 , &
                            &", DX = ",F6.4,", DY = ",F6.4 , &
                            &", REF_LAT = ",F8.3,", REF_LON = ",F8.3)') &
                            im(1),jm(1),dlmd,dphd,rlat0d,-1.*rlon0d
   CALL getset(xa,xb,ya,yb,xxa,xxy,yya,yyb,ltype)
!  CALL set   (xa,xb,ya,yb,0.,1.,0.,1.,lytpe)
!  CALL gsclip(0)
!  CALL pchiqu ( 0.5 , -5. , string , 8. , 0., 0. )
   CALL pchiqu ( xxa , yya - (yyb-yya)/20., string , 8. , 0., -1. )
   CALL frame

   !  Close workstation and NCAR Grpahics.

   CALL gclwk(8)
   CALL clsgks

END SUBROUTINE mapbkg_bgrid

SUBROUTINE mapbkg_egrid ( imt,jm,ilowl,jlowl,ngpwe,ngpsn,&
                          rlat0d,rlon0d,glatl,glonl,glatr,glonr,&
                          dlmd,dphd)

!  IMPLICIT NONE

   !  Some local vars

   CHARACTER (LEN=97) :: string

   !  Yet more center lon messing around, hoo boy.

!  clonx=180.-rlon0d
   clonx=-rlon0d

   !  Open up NCAR Graphics

   CALL opngks
   CALL gopwk(8,9,3)
   CALL gsclip(0)

   !  Make the background white, and the foreground black.

   CALL gscr ( 1 , 0 , 1., 1., 1. )
   CALL gscr ( 1 , 1 , 0., 0., 0. )

   !  Line width twice default thickness.

   CALL setusv('LW',2000)

   !  Make map outline a solid line, not dots.

   CALL mapsti('MV',8)
   CALL mapsti('DO',0)

   !  Map outlines are political and states.

   CALL mapstc('OU','PS')

   !  Cylindrical equidistant.

   CALL maproj('CE',rlat0d,clonx,0.)

   !  Specify corner points.

   CALL mapset('CO',glatl,glonl,glatr,glonr)

   !  Lat lon lines every 5 degrees.

   CALL mapsti('GR',5)

   !  Map takes up this much real estate.

   CALL mappos( 0.05 , 0.95 , 0.05 , 0.95 )

   !  Initialize and draw map.

   CALL mapint
   CALL mapdrw

   !  Add approx grid point tick marks

   CALL setusv('LW',1000)
   CALL perim(((imt+3)/2)-1,1,jm-1,1)

   !  Put on a quicky description.

   WRITE ( string , FMT = '("E-GRID E_WE = ",I4,", E_SN = ",I4 , &
                            &", DX = ",F6.4,", DY = ",F6.4 , &
                            &", REF_LAT = ",F8.3,", REF_LON = ",F8.3)') &
                            (imt+3)/2,jm+1,dlmd,dphd,rlat0d,-1.*rlon0d
   CALL getset(xa,xb,ya,yb,xxa,xxy,yya,yyb,ltype)
!  CALL set   (xa,xb,ya,yb,0.,1.,0.,1.,lytpe)
!  CALL gsclip(0)
!  CALL pchiqu ( 0.5 , -5. , string , 8. , 0., 0. )
   CALL pchiqu ( xxa , yya - (yyb-yya)/20., string , 8. , 0., -1. )
   CALL frame

   !  Close workstation and NCAR Grpahics.

   CALL gclwk(8)
   CALL clsgks

END SUBROUTINE mapbkg_egrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE e2t2g ( ncol,nrow,im,jm,tph0d,tlm0d,dphd,dlmd,glatd,glond)

!  IMPLICIT NONE

   REAL , PARAMETER :: D2R=1.74532925E-2 , R2D=1./D2R

   DPH=DPHD*D2R
   DLM=DLMD*D2R
   TPH0=TPH0D*D2R
   TLM0=TLM0D*D2R

!***  FIND THE TRANSFORMED LAT (POSITIVE NORTH) AND LON (POSITIVE EAST)

   TLATD=(NROW-(JM+1)/2)*DPHD
   TLOND=(NCOL-IM)*DLMD

!***  NOW CONVERT TO GEODETIC LAT (POSITIVE NORTH) AND LON (POSITIVE EAST)

   TLATR=TLATD*D2R
   TLONR=TLOND*D2R
   ARG1=SIN(TLATR)*COS(TPH0)+COS(TLATR)*SIN(TPH0)*COS(TLONR)
   GLATR=ASIN(ARG1)
   GLATD=GLATR*R2D
   ARG2=COS(TLATR)*COS(TLONR)/(COS(GLATR)*COS(TPH0))- & 
        TAN(GLATR)*TAN(TPH0)
   IF(ABS(ARG2).GT.1.)ARG2=ABS(ARG2)/ARG2
   FCTR=1.
   IF(TLOND.GT.0.)FCTR=-1.
   GLOND=TLM0D+FCTR*ACOS(ARG2)*R2D
   GLOND=-GLOND

END SUBROUTINE e2t2g

SUBROUTINE b2t2g ( ncol,nrow,im,jm,tph0d,tlm0d,dphd,dlmd,glatd,glond)

!  IMPLICIT NONE

   REAL , PARAMETER :: D2R=1.74532925E-2 , R2D=1./D2R

   DPH=DPHD*D2R
   DLM=DLMD*D2R
   TPH0=TPH0D*D2R
   TLM0=TLM0D*D2R

!***  FIND THE TRANSFORMED LAT (POSITIVE NORTH) AND LON (POSITIVE EAST)

   TLATD=(NROW-(JM+1)/2)*DPHD
   TLOND=(NCOL-(IM+1)/2)*DLMD

!***  NOW CONVERT TO GEODETIC LAT (POSITIVE NORTH) AND LON (POSITIVE EAST)

   TLATR=TLATD*D2R
   TLONR=TLOND*D2R
   ARG1=SIN(TLATR)*COS(TPH0)+COS(TLATR)*SIN(TPH0)*COS(TLONR)
   GLATR=ASIN(ARG1)
   GLATD=GLATR*R2D
   ARG2=COS(TLATR)*COS(TLONR)/(COS(GLATR)*COS(TPH0))- & 
        TAN(GLATR)*TAN(TPH0)
   IF(ABS(ARG2).GT.1.)ARG2=ABS(ARG2)/ARG2
   FCTR=1.
   IF(TLOND.GT.0.)FCTR=-1.
   GLOND=TLM0D+FCTR*ACOS(ARG2)*R2D
   GLOND=-GLOND

END SUBROUTINE b2t2g

   subroutine tll(almd,aphd,tlmd,tphd,tph0d,tlm0d)
!-------------------------------------------------------------------------------
      real, intent(in) :: almd, aphd
      real, intent(out) :: tlmd, tphd
      real, intent(in) :: tph0d, tlm0d
!-------------------------------------------------------------------------------
      real, parameter :: pi=3.141592654
      real, parameter :: dtr=pi/180.0
!
      real :: tph0, ctph0, stph0, relm, srlm, crlm
      real :: aph, sph, cph, cc, anum, denom
!-------------------------------------------------------------------------------
!
      if (tlm0d==0.0.and.tph0d==0.0) then
      tlmd=almd
      tphd=aphd
      else

      tph0=tph0d*dtr
      ctph0=cos(tph0)
      stph0=sin(tph0)
!
      relm=(almd-tlm0d)*dtr
      srlm=sin(relm)
      crlm=cos(relm)
      aph=aphd*dtr
      sph=sin(aph)
      cph=cos(aph)
      cc=cph*crlm
      anum=cph*srlm
      denom=ctph0*cc+stph0*sph
!
      tlmd=atan2(anum,denom)/dtr
      tphd=asin(ctph0*sph-stph0*cc)/dtr

      end if
!
      return
!
   end subroutine tll

