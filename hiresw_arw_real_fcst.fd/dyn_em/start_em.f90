!-------------------------------------------------------------------

   SUBROUTINE start_domain_em ( grid, allowed_to_read  &
! Actual arguments generated from Registry
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/em_dummy_new_args.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
,u_b,u_bt,v_b,v_bt,w_b,w_bt,ph_b,ph_bt,t_b,t_bt,mu_b,mu_bt,moist,moist_b,moist_bt,chem,scalar,scalar_b,scalar_bt,ozmixm, &
aerosolc_1,aerosolc_2,fdda3d,fdda2d &
!ENDOFREGISTRYGENERATEDINCLUDE
!
)

   USE module_domain
   USE module_dm
!  USE module_io_domain
   USE module_state_description
   USE module_model_constants
   USE module_bc
   USE module_bc_em
!  USE module_timing
   USE module_configure
   USE module_tiles

   USE module_physics_init

   USE module_dm

!!debug
!USE module_compute_geop

   USE module_model_constants
   IMPLICIT NONE
   !  Input data.
   TYPE (domain)          :: grid

   LOGICAL , INTENT(IN)   :: allowed_to_read

   !  Definitions of dummy arguments to this routine (generated from Registry).
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/em_dummy_new_decl.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: u_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: u_bt
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: v_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: v_bt
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: w_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: w_bt
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: ph_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: ph_bt
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: t_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: t_bt
real      ,DIMENSION(max(grid%ed31,grid%ed33),1,grid%spec_bdy_width,4)           :: mu_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),1,grid%spec_bdy_width,4)           :: mu_bt
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4,num_moist)           :: moist_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4,num_moist)           :: moist_bt
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_chem)           :: chem
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_scalar)           :: scalar
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4,num_scalar)           :: scalar_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4,num_scalar)           :: scalar_bt
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%levsiz,grid%sm33:grid%em33,num_ozmixm)           :: ozmixm
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%paerlev,grid%sm33:grid%em33,num_aerosolc)           :: aerosolc_1
real      ,DIMENSION(grid%sm31:grid%em31,1:grid%paerlev,grid%sm33:grid%em33,num_aerosolc)           :: aerosolc_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_fdda3d)           :: fdda3d
real      ,DIMENSION(grid%sm31:grid%em31,1:1,grid%sm33:grid%em33,num_fdda2d)           :: fdda2d
!ENDOFREGISTRYGENERATEDINCLUDE

   !  Structure that contains run-time configuration (namelist) data for domain
   TYPE (grid_config_rec_type)              :: config_flags

   !  Local data
   INTEGER                             ::                       &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  ips, ipe, jps, jpe, kps, kpe, &
                                  its, ite, jts, jte, kts, kte, &
                                  ij,i,j,k,ii,jj,kk,loop,error,l

   INTEGER     :: i_m

   REAL        :: p00, t00, a, p_surf, pd_surf

   REAL :: qvf1, qvf2, qvf
   REAL :: MPDT
   REAL :: spongeweight
   LOGICAL :: first_trip_for_this_domain, start_of_simulation

   REAL :: lat1 , lat2 , lat3 , lat4
   REAL :: lon1 , lon2 , lon3 , lon4
   INTEGER :: num_points_lat_lon , iloc , jloc
   CHARACTER (LEN=132) :: message

! Needed by some comm layers, grid%e.g. RSL. If needed, nmm_data_calls.inc is
! generated from the registry.  The definition of REGISTER_I1 allows
! I1 data to be communicated in this routine if necessary.
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/em_data_calls.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
!ENDOFREGISTRYGENERATEDINCLUDE
   CALL get_ijk_from_grid (  grid ,                   &
                             ids, ide, jds, jde, kds, kde,    &
                             ims, ime, jms, jme, kms, kme,    &
                             ips, ipe, jps, jpe, kps, kpe    )

   kts = kps ; kte = kpe     ! note that tile is entire patch
   its = ips ; ite = ipe     ! note that tile is entire patch
   jts = jps ; jte = jpe    ! note that tile is entire patch

   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )

   IF ( ( MOD (ide-ids,config_flags%parent_grid_ratio) .NE. 0 ) .OR. &
        ( MOD (jde-jds,config_flags%parent_grid_ratio) .NE. 0 ) ) THEN
      WRITE(message, FMT='("Nested dimensions are illegal for domain ",I2,":  Both &
         &MOD(",I4,"-",I1,",",I2,") and MOD(",I4,"-",I1,",",I2,") must = 0" )') &
         grid%id,ide,ids,config_flags%parent_grid_ratio,jde,jds,config_flags%parent_grid_ratio
      CALL wrf_error_fatal3 ( "start_em.b" , 97 ,  message ) 
   END IF

! here we check to see if the boundary conditions are set properly

   CALL boundary_condition_check( config_flags, bdyzone, error, grid%id )

!kludge - need to stop CG from resetting precip and phys tendencies to zero
!         when we are in here due to a nest being spawned, we want to still
!         recompute the base state, but that is about it
   !  This is temporary and will need to be changed when grid%itimestep is removed.

   IF ( grid%itimestep .EQ. 0 ) THEN
      first_trip_for_this_domain = .TRUE.
   ELSE
      first_trip_for_this_domain = .FALSE.
   END IF

   IF ( .not. ( config_flags%restart .or. grid%moved ) ) THEN 
       grid%itimestep=0
   ENDIF

   IF ( config_flags%restart .or. grid%moved ) THEN 
       first_trip_for_this_domain = .TRUE.
   ENDIF

   IF (config_flags%specified) THEN
!
! Arrays for specified boundary conditions
! wig: Add a combined exponential+linear weight on the mother boundaries
!      following code changes by Ruby Leung. For the nested grid, there
!      appears to be some problems when a sponge is used. The points where
!      processors meet have problematic values.

     DO loop = grid%spec_zone + 1, grid%spec_zone + grid%relax_zone
       grid%fcx(loop) = 0.1 / grid%dt * (grid%spec_zone + grid%relax_zone - loop) / (grid%relax_zone - 1)
       grid%gcx(loop) = 1.0 / grid%dt / 50. * (grid%spec_zone + grid%relax_zone - loop) / (grid%relax_zone - 1)
!      spongeweight=exp(-(loop-2)/3.)
!      grid%fcx(loop) = grid%fcx(loop)*spongeweight
!      grid%gcx(loop) = grid%gcx(loop)*spongeweight
     ENDDO

   ELSE IF (config_flags%nested) THEN
!
! Arrays for specified boundary conditions

     DO loop = grid%spec_zone + 1, grid%spec_zone + grid%relax_zone
       grid%fcx(loop) = 0.1 / grid%dt * (grid%spec_zone + grid%relax_zone - loop) / (grid%relax_zone - 1)
       grid%gcx(loop) = 1.0 / grid%dt / 50. * (grid%spec_zone + grid%relax_zone - loop) / (grid%relax_zone - 1)
!      spongeweight=exp(-(loop-2)/3.)
!      grid%fcx(loop) = grid%fcx(loop)*spongeweight
!      grid%gcx(loop) = grid%gcx(loop)*spongeweight
!      grid%fcx(loop) = 0.
!      grid%gcx(loop) = 0.
     ENDDO

     grid%dtbc = 0.

   ENDIF

!     IF ( grid%id .NE. 1 ) THEN
   
      !  Every time a domain starts or every time a domain moves, this routine is called.  We want
      !  the center (middle) lat/lon of the grid for the metacode.  The lat/lon values are
      !  defined at mass points.  Depending on the even/odd points in the SN and WE directions,
      !  we end up with the middle point as either 1 point or an average of either 2 or 4 points.
      !  Add to this, the need to make sure that we are on the correct patch to retrieve the
      !  value of the lat/lon, AND that the lat/lons (for an average) may not all be on the same
      !  patch.  Once we find the correct value for lat lon, we need to keep it around on all patches, 
      !  which is where the wrf_dm_min_real calls come in.
   
      IF      ( ( MOD(ide,2) .EQ. 0 ) .AND. ( MOD(jde,2) .EQ. 0 ) ) THEN
         num_points_lat_lon = 1
         iloc = ide/2
         jloc = jde/2
         IF ( ( ips .LE. iloc ) .AND. ( ipe .GE. iloc ) .AND. &
              ( jps .LE. jloc ) .AND. ( jpe .GE. jloc ) ) THEN
            lat1 = grid%xlat (iloc,jloc)
            lon1 = grid%xlong(iloc,jloc)
         ELSE
            lat1 = 99999.
            lon1 = 99999.
         END IF
         lat1 = wrf_dm_min_real ( lat1 )
         lon1 = wrf_dm_min_real ( lon1 )
         CALL nl_set_cen_lat ( grid%id , lat1 )
         CALL nl_set_cen_lon ( grid%id , lon1 )
      ELSE IF ( ( MOD(ide,2) .NE. 0 ) .AND. ( MOD(jde,2) .EQ. 0 ) ) THEN
         num_points_lat_lon = 2
         iloc = (ide-1)/2
         jloc =  jde   /2
         IF ( ( ips .LE. iloc ) .AND. ( ipe .GE. iloc ) .AND. &
              ( jps .LE. jloc ) .AND. ( jpe .GE. jloc ) ) THEN
            lat1 = grid%xlat (iloc,jloc)
            lon1 = grid%xlong(iloc,jloc)
         ELSE
            lat1 = 99999.
            lon1 = 99999.
         END IF
         lat1 = wrf_dm_min_real ( lat1 )
         lon1 = wrf_dm_min_real ( lon1 )
   
         iloc = (ide+1)/2
         jloc =  jde   /2
         IF ( ( ips .LE. iloc ) .AND. ( ipe .GE. iloc ) .AND. &
              ( jps .LE. jloc ) .AND. ( jpe .GE. jloc ) ) THEN
            lat2 = grid%xlat (iloc,jloc)
            lon2 = grid%xlong(iloc,jloc)
         ELSE
            lat2 = 99999.
            lon2 = 99999.
         END IF
         lat2 = wrf_dm_min_real ( lat2 )
         lon2 = wrf_dm_min_real ( lon2 )
   
         CALL nl_set_cen_lat ( grid%id , ( lat1 + lat2 ) * 0.50 )
         CALL nl_set_cen_lon ( grid%id , ( lon1 + lon2 ) * 0.50 )
      ELSE IF ( ( MOD(ide,2) .EQ. 0 ) .AND. ( MOD(jde,2) .NE. 0 ) ) THEN
         num_points_lat_lon = 2
         iloc =  ide   /2
         jloc = (jde-1)/2
         IF ( ( ips .LE. iloc ) .AND. ( ipe .GE. iloc ) .AND. &
              ( jps .LE. jloc ) .AND. ( jpe .GE. jloc ) ) THEN
            lat1 = grid%xlat (iloc,jloc)
            lon1 = grid%xlong(iloc,jloc)
         ELSE
            lat1 = 99999.
            lon1 = 99999.
         END IF
         lat1 = wrf_dm_min_real ( lat1 )
         lon1 = wrf_dm_min_real ( lon1 )
   
         iloc =  ide   /2
         jloc = (jde+1)/2
         IF ( ( ips .LE. iloc ) .AND. ( ipe .GE. iloc ) .AND. &
              ( jps .LE. jloc ) .AND. ( jpe .GE. jloc ) ) THEN
            lat2 = grid%xlat (iloc,jloc)
            lon2 = grid%xlong(iloc,jloc)
         ELSE
            lat2 = 99999.
            lon2 = 99999.
         END IF
         lat2 = wrf_dm_min_real ( lat2 )
         lon2 = wrf_dm_min_real ( lon2 )
   
         CALL nl_set_cen_lat ( grid%id , ( lat1 + lat2 ) * 0.50 )
         CALL nl_set_cen_lon ( grid%id , ( lon1 + lon2 ) * 0.50 )
      ELSE IF ( ( MOD(ide,2) .NE. 0 ) .AND. ( MOD(jde,2) .NE. 0 ) ) THEN
         num_points_lat_lon = 4
         iloc = (ide-1)/2
         jloc = (jde-1)/2
         IF ( ( ips .LE. iloc ) .AND. ( ipe .GE. iloc ) .AND. &
              ( jps .LE. jloc ) .AND. ( jpe .GE. jloc ) ) THEN
            lat1 = grid%xlat (iloc,jloc)
            lon1 = grid%xlong(iloc,jloc)
         ELSE
            lat1 = 99999.
            lon1 = 99999.
         END IF
         lat1 = wrf_dm_min_real ( lat1 )
         lon1 = wrf_dm_min_real ( lon1 )
   
         iloc = (ide+1)/2
         jloc = (jde-1)/2
         IF ( ( ips .LE. iloc ) .AND. ( ipe .GE. iloc ) .AND. &
              ( jps .LE. jloc ) .AND. ( jpe .GE. jloc ) ) THEN
            lat2 = grid%xlat (iloc,jloc)
            lon2 = grid%xlong(iloc,jloc)
         ELSE
            lat2 = 99999.
            lon2 = 99999.
         END IF
         lat2 = wrf_dm_min_real ( lat2 )
         lon2 = wrf_dm_min_real ( lon2 )
   
         iloc = (ide-1)/2
         jloc = (jde+1)/2
         IF ( ( ips .LE. iloc ) .AND. ( ipe .GE. iloc ) .AND. &
              ( jps .LE. jloc ) .AND. ( jpe .GE. jloc ) ) THEN
            lat3 = grid%xlat (iloc,jloc)
            lon3 = grid%xlong(iloc,jloc)
         ELSE
            lat3 = 99999.
            lon3 = 99999.
         END IF
         lat3 = wrf_dm_min_real ( lat3 )
         lon3 = wrf_dm_min_real ( lon3 )
   
         iloc = (ide+1)/2
         jloc = (jde+1)/2
         IF ( ( ips .LE. iloc ) .AND. ( ipe .GE. iloc ) .AND. &
              ( jps .LE. jloc ) .AND. ( jpe .GE. jloc ) ) THEN
            lat4 = grid%xlat (iloc,jloc)
            lon4 = grid%xlong(iloc,jloc)
         ELSE
            lat4 = 99999.
            lon4 = 99999.
         END IF
         lat4 = wrf_dm_min_real ( lat4 )
         lon4 = wrf_dm_min_real ( lon4 )
   
         CALL nl_set_cen_lat ( grid%id , ( lat1 + lat2 + lat3 + lat4 ) * 0.25 )
         CALL nl_set_cen_lon ( grid%id , ( lon1 + lon2 + lon3 + lon4 ) * 0.25 )
      END IF
!  END IF

   IF ( .NOT. config_flags%restart .AND. &
        (( config_flags%input_from_hires ) .OR. ( config_flags%input_from_file ))) THEN

      IF ( config_flags%map_proj .EQ. 0 ) THEN
         CALL wrf_error_fatal3 ( "start_em.b" , 307 ,  'start_domain: Idealized case cannot have a separate nested input file' )
      END IF 

      CALL nl_get_base_pres  ( 1 , p00 )
      CALL nl_get_base_temp  ( 1 , t00 )
      CALL nl_get_base_lapse ( 1 , a   )

      !  Base state potential temperature and inverse density (alpha = 1/rho) from
      !  the half eta levels and the base-profile surface pressure.  Compute 1/rho
      !  from equation of state.  The potential temperature is a perturbation from t0.

      DO j = jts, MIN(jte,jde-1)
         DO i = its, MIN(ite,ide-1)

            !  Base state pressure is a function of eta level and terrain, only, plus
            !  the hand full of constants: p00 (sea level pressure, Pa), t00 (sea level
            !  temperature, K), and A (temperature difference, from 1000 mb to 300 mb, K).

            p_surf = p00 * EXP ( -t00/a + ( (t00/a)**2 - 2.*g*grid%ht(i,j)/a/r_d ) **0.5 )

            DO k = 1, kte-1
               grid%em_pb(i,k,j) = grid%em_znu(k)*(p_surf - grid%p_top) + grid%p_top
               grid%em_t_init(i,k,j) = (t00 + A*LOG(grid%em_pb(i,k,j)/p00))*(p00/grid%em_pb(i,k,j))**(r_d/cp) - t0
               grid%em_alb(i,k,j) = (r_d/p1000mb)*(grid%em_t_init(i,k,j)+t0)*(grid%em_pb(i,k,j)/p1000mb)**cvpm
            END DO

            !  Base state mu is defined as base state surface pressure minus grid%p_top

            grid%em_mub(i,j) = p_surf - grid%p_top

            !  Integrate base geopotential, starting at terrain elevation.  This assures that
            !  the base state is in exact hydrostatic balance with respect to the model equations.
            !  This field is on full levels.

            grid%em_phb(i,1,j) = grid%ht(i,j) * g
            DO k  = 2,kte
               grid%em_phb(i,k,j) = grid%em_phb(i,k-1,j) - grid%em_dnw(k-1)*grid%em_mub(i,j)*grid%em_alb(i,k-1,j)
            END DO
         END DO
      END DO

   ENDIF

   IF(.not.config_flags%restart)THEN

!  if this is for a nested domain, the defined/interpolated fields are the _2

     IF ( first_trip_for_this_domain ) THEN

! data that is expected to be zero must be explicitly initialized as such
     grid%h_diabatic = 0.

       DO j = jts,min(jte,jde-1)
       DO k = kts,kte-1
       DO i = its, min(ite,ide-1)
           IF ( grid%imask_nostag(i,j) .EQ. 1 ) THEN
             grid%em_t_1(i,k,j)=grid%em_t_2(i,k,j)
           ENDIF
       ENDDO
       ENDDO
       ENDDO
  
       DO j = jts,min(jte,jde-1)
       DO i = its, min(ite,ide-1)
           IF ( grid%imask_nostag(i,j) .EQ. 1 ) THEN
             grid%em_mu_1(i,j)=grid%em_mu_2(i,j)
           ENDIF
       ENDDO
       ENDDO
     END IF

!  reconstitute base-state fields

    IF(config_flags%max_dom .EQ. 1)THEN
! with single domain, grid%em_t_init from wrfinput is OK to use
     DO j = jts,min(jte,jde-1)
     DO k = kts,kte-1
     DO i = its, min(ite,ide-1)
       IF ( grid%imask_nostag(i,j) .EQ. 1 ) THEN
         grid%em_pb(i,k,j) = grid%em_znu(k)*grid%em_mub(i,j)+grid%p_top
         grid%em_alb(i,k,j) = (r_d/p1000mb)*(grid%em_t_init(i,k,j)+t0)*(grid%em_pb(i,k,j)/p1000mb)**cvpm
       ENDIF
     ENDDO
     ENDDO
     ENDDO
    ELSE
! with nests, grid%em_t_init generally needs recomputations (since it is not interpolated)
     DO j = jts,min(jte,jde-1)
     DO k = kts,kte-1
     DO i = its, min(ite,ide-1)
       IF ( grid%imask_nostag(i,j) .EQ. 1 ) THEN
         grid%em_pb(i,k,j) = grid%em_znu(k)*grid%em_mub(i,j)+grid%p_top
         grid%em_alb(i,k,j) = -grid%em_rdnw(k)*(grid%em_phb(i,k+1,j)-grid%em_phb(i,k,j))/grid%em_mub(i,j)
         grid%em_t_init(i,k,j) = grid%em_alb(i,k,j)*(p1000mb/r_d)/((grid%em_pb(i,k,j)/p1000mb)**cvpm) - t0
       ENDIF
     ENDDO
     ENDDO
     ENDDO
    ENDIF

     DO j = jts,min(jte,jde-1)

       k = kte-1
       DO i = its, min(ite,ide-1)
         IF ( grid%imask_nostag(i,j) .EQ. 1 ) THEN
           qvf1 = 0.5*(moist(i,k,j,P_QV)+moist(i,k,j,P_QV))
           qvf2 = 1./(1.+qvf1)
           qvf1 = qvf1*qvf2
           grid%em_p(i,k,j) = - 0.5*(grid%em_mu_1(i,j)+qvf1*grid%em_mub(i,j))/grid%em_rdnw(k)/qvf2
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%em_alt(i,k,j) = (r_d/p1000mb)*(grid%em_t_1(i,k,j)+t0)*qvf*(((grid%em_p(i,k,j)+grid%em_pb(i,k,j))/p1000mb)**cvpm)
           grid%em_al(i,k,j) = grid%em_alt(i,k,j) - grid%em_alb(i,k,j)
         ENDIF
       ENDDO

       DO k = kte-2, 1, -1
       DO i = its, min(ite,ide-1)
         IF ( grid%imask_nostag(i,j) .EQ. 1 ) THEN
           qvf1 = 0.5*(moist(i,k,j,P_QV)+moist(i,k+1,j,P_QV))
           qvf2 = 1./(1.+qvf1)
           qvf1 = qvf1*qvf2
           grid%em_p(i,k,j) = grid%em_p(i,k+1,j) - (grid%em_mu_1(i,j) + qvf1*grid%em_mub(i,j))/qvf2/grid%em_rdn(k+1)
           qvf = 1. + rvovrd*moist(i,k,j,P_QV)
           grid%em_alt(i,k,j) = (r_d/p1000mb)*(grid%em_t_1(i,k,j)+t0)*qvf* &
                        (((grid%em_p(i,k,j)+grid%em_pb(i,k,j))/p1000mb)**cvpm)
           grid%em_al(i,k,j) = grid%em_alt(i,k,j) - grid%em_alb(i,k,j)
         ENDIF
       ENDDO
       ENDDO

     ENDDO

   ENDIF

   IF ( ( grid%id .NE. 1 ) .AND. .NOT. ( config_flags%restart ) .AND. &
       ( ( config_flags%input_from_hires ) .OR. ( config_flags%input_from_file ) ) ) THEN
      DO j = jts, MIN(jte,jde-1)
         DO i = its, MIN(ite,ide-1)
            grid%em_mu_2(i,j) = grid%em_mu_2(i,j) + grid%em_al(i,1,j) / ( grid%em_alt(i,1,j) * grid%em_alb(i,1,j) ) * &
                                    g * ( grid%ht(i,j) - grid%ht_fine(i,j) )
         END DO
       END DO
       DO j = jts,min(jte,jde-1)
       DO i = its, min(ite,ide-1)
          grid%em_mu_1(i,j)=grid%em_mu_2(i,j)
       ENDDO
       ENDDO

   END IF

   IF ( first_trip_for_this_domain ) THEN

   CALL wrf_debug ( 100 , 'module_start: start_domain_rk: Before call to phy_init' )

! namelist MPDT does not exist yet, so set it here
! MPDT is the call frequency for microphysics in minutes (0 means every step)
   MPDT = 0.

! set GMT outside of phy_init because phy_init may not be called on this
! process if, for example, it is a moving nest and if this part of the domain is not
! being initialized (not the leading edge).
   CALL domain_setgmtetc( grid, start_of_simulation )

   CALL set_tiles ( grid , grid%imask_nostag, ims, ime, jms, jme, ips, ipe, jps, jpe )

! Phy_init is not necessarily thread-safe; do not multi-thread this loop.
! The tiling is to handle the fact that we may be masking off part of the computation.
   DO ij = 1, grid%num_tiles

     CALL phy_init (  grid%id , config_flags, grid%DT, grid%RESTART, grid%em_znw, grid%em_znu,   &
                      grid%p_top, grid%tsk, grid%RADT,grid%BLDT,grid%CUDT, MPDT, &
                      grid%rthcuten, grid%rqvcuten, grid%rqrcuten,           &
                      grid%rqccuten, grid%rqscuten, grid%rqicuten,           &
                      grid%rublten,grid%rvblten,grid%rthblten,               &
                      grid%rqvblten,grid%rqcblten,grid%rqiblten,             &
                      grid%rthraten,grid%rthratenlw,grid%rthratensw,         &
                      grid%stepbl,grid%stepra,grid%stepcu,                   &
                      grid%w0avg, grid%rainnc, grid%rainc, grid%raincv, grid%rainncv,  &
                      grid%nca,grid%swrad_scat,                    &
                      grid%cldefi,grid%lowlyr,                          &
                      grid%mass_flux,                              &
                      grid%rthften, grid%rqvften,                       &
                      grid%cldfra,grid%glw,grid%gsw,grid%emiss,grid%lu_index,          &
                      grid%landuse_ISICE, grid%landuse_LUCATS,            &
                      grid%landuse_LUSEAS, grid%landuse_ISN,              &
                      grid%lu_state,                                      &
                      grid%xlat,grid%xlong,grid%albedo,grid%albbck,grid%GMT,grid%JULYR,grid%JULDAY,     &
                      grid%levsiz, num_ozmixm, num_aerosolc, grid%paerlev,  &
                      grid%tmn,grid%xland,grid%znt,grid%z0,grid%ust,grid%mol,grid%pblh,grid%tke_myj,     &
                      grid%exch_h,grid%thc,grid%snowc,grid%mavail,grid%hfx,grid%qfx,grid%rainbl, &
                      grid%tslb,grid%zs,grid%dzs,config_flags%num_soil_layers,grid%warm_rain,  &
                      grid%adv_moist_cond,                         &
                      grid%apr_gr,grid%apr_w,grid%apr_mc,grid%apr_st,grid%apr_as,      &
                      grid%apr_capma,grid%apr_capme,grid%apr_capmi,          &
                      grid%xice,grid%vegfra,grid%snow,grid%canwat,grid%smstav,         &
                      grid%smstot, grid%sfcrunoff,grid%udrunoff,grid%grdflx,grid%acsnow,      &
                      grid%acsnom,grid%ivgtyp,grid%isltyp, grid%sfcevp,grid%smois,     &
                      grid%sh2o, grid%snowh, grid%smfr3d,                    &
                      grid%DX,grid%DY,grid%f_ice_phy,grid%f_rain_phy,grid%f_rimef_phy, &
                      grid%mp_restart_state,grid%tbpvs_state,grid%tbpvs0_state,&
                      allowed_to_read, grid%moved, start_of_simulation,               &
                      ids, ide, jds, jde, kds, kde,           &
                      ims, ime, jms, jme, kms, kme,           &
                      grid%i_start(ij), grid%i_end(ij), grid%j_start(ij), grid%j_end(ij), kts, kte, &
                      ozmixm,grid%pin,                             &     ! Optional
                      grid%m_ps_1,grid%m_ps_2,grid%m_hybi,aerosolc_1,aerosolc_2,&  ! Optional
                      grid%rundgdten,grid%rvndgdten,grid%rthndgdten,         &     ! Optional
                      grid%rqvndgdten,grid%rmundgdten,                  &     ! Optional
                      grid%FGDT,grid%stepfg,                        &     ! Optional
                      grid%DZR, grid%DZB, grid%DZG,                          & !Optional urban
                      grid%TR_URB2D,grid%TB_URB2D,grid%TG_URB2D,grid%TC_URB2D,    & !Optional urban
                      grid%QC_URB2D, grid%XXXR_URB2D,grid%XXXB_URB2D,        & !Optional urban
                      grid%XXXG_URB2D, grid%XXXC_URB2D,                 & !Optional urban
                      grid%TRL_URB3D, grid%TBL_URB3D, grid%TGL_URB3D,        & !Optional urban
                      grid%SH_URB2D, grid%LH_URB2D, grid%G_URB2D, grid%RN_URB2D,  & !Optional urban
                      grid%TS_URB2D, grid%FRC_URB2D, grid%UTYPE_URB2D,      & !Optional urban
                      itimestep=grid%itimestep, fdob=grid%fdob             &
                      )

   ENDDO


 
   CALL wrf_debug ( 100 , 'module_start: start_domain_rk: After call to phy_init' )


   END IF


!
!

 !  set physical boundary conditions for all initialized variables

!-----------------------------------------------------------------------
!  Stencils for patch communications  (WCS, 29 June 2001)
!  Note:  the size of this halo exchange reflects the 
!         fact that we are carrying the uncoupled variables 
!         as state variables in the mass coordinate model, as
!         opposed to the coupled variables as in the height
!         coordinate model.
!
!                           * * * * *
!         *        * * *    * * * * *
!       * + *      * + *    * * + * * 
!         *        * * *    * * * * *
!                           * * * * *
!
!j  grid%em_u_1                          x
!j  grid%em_u_2                          x
!j  grid%em_v_1                          x
!j  grid%em_v_2                          x
!j  grid%em_w_1                          x
!j  grid%em_w_2                          x
!j  grid%em_t_1                          x
!j  grid%em_t_2                          x
!j grid%em_ph_1                          x
!j grid%em_ph_2                          x
!
!j  grid%em_t_init                       x
!
!j  grid%em_phb   x
!j  grid%em_ph0   x
!j  grid%em_php   x
!j   grid%em_pb   x
!j   grid%em_al   x
!j  grid%em_alt   x
!j  grid%em_alb   x
!
!  the following are 2D (xy) variables
!
!j grid%em_mu_1                          x
!j grid%em_mu_2                          x
!j grid%em_mub   x
!j grid%em_mu0   x
!j grid%ht    x
!j grid%msft  x
!j grid%msfu  x
!j grid%msfv  x
!j grid%sina  x
!j grid%cosa  x
!j grid%e     x
!j grid%f     x
!
!  4D variables
!
! moist                        x
!  chem                        x
!scalar                        x

!--------------------------------------------------------------

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_INIT_1.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_EM_INIT_1.inc')
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3, &
     8, 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_1, 3, 4, 0, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_2, 3, 4, 0, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_1, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_2, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_1, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_2, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_1, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_2, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_1, 3, 4, 0, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_2, 3, 4, 0, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_1, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_2, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_1, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_2, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_1, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_2, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3 , &
     8, 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_1, 3, 4, 1, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_2, 3, 4, 1, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_1, 3, 4, 1, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_2, 3, 4, 1, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_INIT_2.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_EM_INIT_2.inc')
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3, &
     6, 2, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_1, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_2, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_1, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_2, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_1, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_2, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ww, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_phb, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_1, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_2, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_1, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_2, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_1, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_2, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ww, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_phb, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3 , &
     6, 2, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_1, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_2, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ww, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_phb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_1, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_2, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ww, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_phb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_INIT_3.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_EM_INIT_3.inc')
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3, &
     7, 2, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph0, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_php, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_init, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mub, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu0, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_p, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_al, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alt, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alb, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph0, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_php, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_init, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mub, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu0, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_p, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_al, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alt, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alb, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3 , &
     7, 2, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph0, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_php, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_init, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mub, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu0, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_p, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_al, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alt, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph0, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_php, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_init, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mub, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu0, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_p, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_al, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alt, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_INIT_4.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_EM_INIT_4.inc')
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3, &
     2, 11, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_pb, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%h_diabatic, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%msft, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfu, 3, 4, 0, 0, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfv, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%f, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%e, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%sina, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%cosa, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%ht, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_potevp, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_snopcx, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_soiltb, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_pb, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%h_diabatic, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%msft, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfu, 3, 4, 0, 1, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfv, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%f, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%e, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%sina, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%cosa, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%ht, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_potevp, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_snopcx, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_soiltb, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3 , &
     2, 11, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_pb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%h_diabatic, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%msft, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfu, 3, 4, 1, 0, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfv, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%f, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%e, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%sina, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%cosa, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%ht, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_potevp, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_snopcx, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_soiltb, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_pb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%h_diabatic, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%msft, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfu, 3, 4, 1, 1, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfv, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%f, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%e, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%sina, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%cosa, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%ht, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_potevp, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_snopcx, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_soiltb, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_INIT_5.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_EM_INIT_5.inc')
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3, &
     0  &
   + num_moist   &
   + num_chem   &
   + num_scalar   &
     , 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK ( local_communicator,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK ( local_communicator,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 CALL RSL_LITE_PACK ( local_communicator,scalar ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK ( local_communicator,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK ( local_communicator,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 CALL RSL_LITE_PACK ( local_communicator,scalar ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3 , &
     0  &
   + num_moist   &
   + num_chem   &
   + num_scalar   &
     , 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK ( local_communicator,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK ( local_communicator,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 CALL RSL_LITE_PACK ( local_communicator,scalar ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK ( local_communicator,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK ( local_communicator,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 CALL RSL_LITE_PACK ( local_communicator,scalar ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/PERIOD_BDY_EM_INIT.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/PERIOD_BDY_EM_INIT.inc')
IF ( config_flags%periodic_x ) THEN
CALL RSL_LITE_INIT_PERIOD ( local_communicator_periodic, 3 , &
     18, 12, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, kpe    )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_u_1, 3, 4, 1, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_u_2, 3, 4, 1, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_v_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_v_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_w_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_w_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_init, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_phb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph0, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_php, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_pb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_al, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_alt, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_alb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu_1, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu_2, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mub, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu0, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%ht, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msft, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msfu, 3, 4, 1, 0, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msfv, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%sina, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%cosa, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%e, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%f, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
   CALL RSL_LITE_EXCH_PERIOD_X ( local_communicator_periodic , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_u_1, 3, 4, 1, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_u_2, 3, 4, 1, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_v_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_v_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_w_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_w_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_init, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_phb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph0, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_php, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_pb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_al, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_alt, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_alb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu_1, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu_2, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mub, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu0, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%ht, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msft, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msfu, 3, 4, 1, 1, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msfv, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%sina, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%cosa, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%e, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%f, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
END IF
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/PERIOD_BDY_EM_MOIST.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/PERIOD_BDY_EM_MOIST.inc')
IF ( config_flags%periodic_x ) THEN
CALL RSL_LITE_INIT_PERIOD ( local_communicator_periodic, 3 , &
     0  &
   + num_moist   &
     , 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, kpe    )
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
   CALL RSL_LITE_EXCH_PERIOD_X ( local_communicator_periodic , mytask, ntasks, ntasks_x, ntasks_y )
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
END IF
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/PERIOD_BDY_EM_CHEM.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/PERIOD_BDY_EM_CHEM.inc')
IF ( config_flags%periodic_x ) THEN
CALL RSL_LITE_INIT_PERIOD ( local_communicator_periodic, 3 , &
     0  &
   + num_chem   &
     , 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, kpe    )
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
   CALL RSL_LITE_EXCH_PERIOD_X ( local_communicator_periodic , mytask, ntasks, ntasks_x, ntasks_y )
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
END IF
!ENDOFREGISTRYGENERATEDINCLUDE


   CALL set_physical_bc3d( grid%em_u_1 , 'U' , config_flags ,                  &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )
   CALL set_physical_bc3d( grid%em_u_2 , 'U' , config_flags ,                  &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )

   CALL set_physical_bc3d( grid%em_v_1 , 'V' , config_flags ,                  &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )
   CALL set_physical_bc3d( grid%em_v_2 , 'V' , config_flags ,                  &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )

! set kinematic condition for w 

   CALL set_physical_bc2d( grid%ht , 'r' , config_flags ,                &
                           ids , ide , jds , jde , &
                           ims , ime , jms , jme , &
                           its , ite , jts , jte , &
                           its , ite , jts , jte   )

   IF ( .not. config_flags%restart ) THEN 
   CALL set_w_surface(  config_flags,                                      &
                        grid%em_w_1, grid%ht, grid%em_u_1, grid%em_v_1, grid%cf1, &
                        grid%cf2, grid%cf3, grid%rdx, grid%rdy, grid%msft,        &
                        ids, ide, jds, jde, kds, kde,                      &
                        ips, ipe, jps, jpe, kps, kpe,                      &
                        its, ite, jts, jte, kts, kte,                      &
                        ims, ime, jms, jme, kms, kme                      )
   CALL set_w_surface(  config_flags,                                      &
                        grid%em_w_2, grid%ht, grid%em_u_2, grid%em_v_2, grid%cf1, &
                        grid%cf2, grid%cf3, grid%rdx, grid%rdy, grid%msft,        &
                        ids, ide, jds, jde, kds, kde,                      &
                        ips, ipe, jps, jpe, kps, kpe,                      &
                        its, ite, jts, jte, kts, kte,                      &
                        ims, ime, jms, jme, kms, kme                      )
   END IF

! finished setting kinematic condition for w at the surface

   CALL set_physical_bc3d( grid%em_w_1 , 'W' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )
   CALL set_physical_bc3d( grid%em_w_2 , 'W' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )

   CALL set_physical_bc3d( grid%em_ph_1 , 'W' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )

   CALL set_physical_bc3d( grid%em_ph_2 , 'W' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )

   CALL set_physical_bc3d( grid%em_t_1 , 't' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )

   CALL set_physical_bc3d( grid%em_t_2 , 't' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )

   CALL set_physical_bc2d( grid%em_mu_1, 't' , config_flags ,   &
                           ids , ide , jds , jde ,  &
                           ims , ime , jms , jme ,  &
                           its , ite , jts , jte ,  &
                           its , ite , jts , jte   )
   CALL set_physical_bc2d( grid%em_mu_2, 't' , config_flags ,   &
                           ids , ide , jds , jde ,  &
                           ims , ime , jms , jme ,  &
                           its , ite , jts , jte ,  &
                           its , ite , jts , jte   )
   CALL set_physical_bc2d( grid%em_mub , 't' , config_flags ,   &
                           ids , ide , jds , jde ,  &
                           ims , ime , jms , jme ,  &
                           its , ite , jts , jte ,  &
                           its , ite , jts , jte   )
   CALL set_physical_bc2d( grid%em_mu0 , 't' , config_flags ,   &
                           ids , ide , jds , jde ,  &
                           ims , ime , jms , jme ,  &
                           its , ite , jts , jte ,  &
                           its , ite , jts , jte   )


   CALL set_physical_bc3d( grid%em_phb , 'W' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )
   CALL set_physical_bc3d( grid%em_ph0 , 'W' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )
   CALL set_physical_bc3d( grid%em_php , 'W' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )

   CALL set_physical_bc3d( grid%em_pb , 't' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )
   CALL set_physical_bc3d( grid%em_al , 't' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )
   CALL set_physical_bc3d( grid%em_alt , 't' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )
   CALL set_physical_bc3d( grid%em_alb , 't' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )
   CALL set_physical_bc3d(grid%em_t_init, 't' , config_flags ,                 &
                         ids , ide , jds , jde , kds , kde ,        &
                         ims , ime , jms , jme , kms , kme ,        &
                         its , ite , jts , jte , kts , kte ,        &
                         its , ite , jts , jte , kts , kte )

   IF (num_moist > 0) THEN

! use of (:,:,:,loop) not efficient on DEC, but (ims,kms,jms,loop) not portable to SGI/Cray

      loop_3d_m   : DO loop = 1 , num_moist
         CALL set_physical_bc3d( moist(:,:,:,loop) , 'r' , config_flags ,                 &
                                 ids , ide , jds , jde , kds , kde ,        &
                                 ims , ime , jms , jme , kms , kme ,        &
                                 its , ite , jts , jte , kts , kte ,        &
                                 its , ite , jts , jte , kts , kte )
      END DO loop_3d_m

   ENDIF


   IF (num_chem >= PARAM_FIRST_SCALAR ) THEN

! use of (:,:,:,loop) not efficient on DEC, but (ims,kms,jms,loop) not portable to SGI/Cray

      loop_3d_c   : DO loop = PARAM_FIRST_SCALAR , num_chem
         CALL set_physical_bc3d( chem(:,:,:,loop) , 'r' , config_flags ,                 &
                                 ids , ide , jds , jde , kds , kde ,        &
                                 ims , ime , jms , jme , kms , kme ,        &
                                 its , ite , jts , jte , kts , kte ,        &
                                 its , ite , jts , jte , kts , kte )
      END DO loop_3d_c

   ENDIF

   CALL set_physical_bc2d( grid%msft , 'r' , config_flags ,              &
                         ids , ide , jds , jde , &
                         ims , ime , jms , jme , &
                         its , ite , jts , jte , &
                         its , ite , jts , jte   )
   CALL set_physical_bc2d( grid%msfu , 'x' , config_flags ,              &
                         ids , ide , jds , jde , &
                         ims , ime , jms , jme , &
                         its , ite , jts , jte , &
                         its , ite , jts , jte   ) 
   CALL set_physical_bc2d( grid%msfv , 'y' , config_flags ,              &
                         ids , ide , jds , jde , &
                         ims , ime , jms , jme , &
                         its , ite , jts , jte , &
                         its , ite , jts , jte   ) 
   CALL set_physical_bc2d( grid%sina , 'r' , config_flags ,              &
                         ids , ide , jds , jde , &
                         ims , ime , jms , jme , &
                         its , ite , jts , jte , &
                         its , ite , jts , jte   ) 
   CALL set_physical_bc2d( grid%cosa , 'r' , config_flags ,              &
                         ids , ide , jds , jde , &
                         ims , ime , jms , jme , &
                         its , ite , jts , jte , &
                         its , ite , jts , jte   ) 
   CALL set_physical_bc2d( grid%e , 'r' , config_flags ,                 &
                         ids , ide , jds , jde , &
                         ims , ime , jms , jme , &
                         its , ite , jts , jte , &
                         its , ite , jts , jte   ) 
   CALL set_physical_bc2d( grid%f , 'r' , config_flags ,                 &
                         ids , ide , jds , jde , &
                         ims , ime , jms , jme , &
                         its , ite , jts , jte , &
                         its , ite , jts , jte   ) 

!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_INIT_1.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_EM_INIT_1.inc')
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3, &
     8, 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_1, 3, 4, 0, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_2, 3, 4, 0, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_1, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_2, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_1, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_2, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_1, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_2, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_1, 3, 4, 0, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_2, 3, 4, 0, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_1, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_2, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_1, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_2, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_1, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_2, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3 , &
     8, 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_1, 3, 4, 1, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_2, 3, 4, 1, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_1, 3, 4, 1, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_u_2, 3, 4, 1, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_v_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_w_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_INIT_2.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_EM_INIT_2.inc')
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3, &
     6, 2, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_1, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_2, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_1, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_2, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_1, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_2, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ww, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_phb, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_1, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_2, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_1, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_2, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_1, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_2, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ww, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_phb, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3 , &
     6, 2, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_1, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_2, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ww, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_phb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_1, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu_2, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_tke_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ww, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_phb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_INIT_3.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_EM_INIT_3.inc')
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3, &
     7, 2, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph0, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_php, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_init, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mub, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu0, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_p, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_al, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alt, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alb, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph0, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_php, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_init, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mub, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu0, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_p, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_al, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alt, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alb, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3 , &
     7, 2, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph0, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_php, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_init, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mub, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu0, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_p, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_al, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alt, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_ph0, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_php, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_t_init, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mub, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_mu0, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_p, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_al, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alt, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%em_alb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_INIT_4.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_EM_INIT_4.inc')
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3, &
     2, 11, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_pb, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%h_diabatic, 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%msft, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfu, 3, 4, 0, 0, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfv, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%f, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%e, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%sina, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%cosa, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%ht, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_potevp, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_snopcx, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_soiltb, 3, 4, 0, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_pb, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%h_diabatic, 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%msft, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfu, 3, 4, 0, 1, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfv, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%f, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%e, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%sina, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%cosa, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%ht, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_potevp, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_snopcx, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_soiltb, 3, 4, 0, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3 , &
     2, 11, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
CALL RSL_LITE_PACK ( local_communicator, grid%em_pb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%h_diabatic, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%msft, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfu, 3, 4, 1, 0, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfv, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%f, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%e, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%sina, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%cosa, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%ht, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_potevp, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_snopcx, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_soiltb, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK ( local_communicator, grid%em_pb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%h_diabatic, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK ( local_communicator, grid%msft, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfu, 3, 4, 1, 1, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%msfv, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%f, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%e, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%sina, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%cosa, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%ht, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_potevp, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_snopcx, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK ( local_communicator, grid%em_soiltb, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/HALO_EM_INIT_5.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/HALO_EM_INIT_5.inc')
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3, &
     0  &
   + num_moist   &
   + num_chem   &
   + num_scalar   &
     , 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK ( local_communicator,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK ( local_communicator,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 CALL RSL_LITE_PACK ( local_communicator,scalar ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
   CALL RSL_LITE_EXCH_Y ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK ( local_communicator,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK ( local_communicator,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 CALL RSL_LITE_PACK ( local_communicator,scalar ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 0, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
CALL RSL_LITE_INIT_EXCH ( local_communicator, 3 , &
     0  &
   + num_moist   &
   + num_chem   &
   + num_scalar   &
     , 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, MAX(1,1&
,kpe &
))
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK ( local_communicator,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK ( local_communicator,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 CALL RSL_LITE_PACK ( local_communicator,scalar ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
   CALL RSL_LITE_EXCH_X ( local_communicator , mytask, ntasks, ntasks_x, ntasks_y )
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK ( local_communicator,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK ( local_communicator,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
DO itrace = PARAM_FIRST_SCALAR, num_scalar
 CALL RSL_LITE_PACK ( local_communicator,scalar ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/PERIOD_BDY_EM_INIT.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/PERIOD_BDY_EM_INIT.inc')
IF ( config_flags%periodic_x ) THEN
CALL RSL_LITE_INIT_PERIOD ( local_communicator_periodic, 3 , &
     18, 12, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, kpe    )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_u_1, 3, 4, 1, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_u_2, 3, 4, 1, 0, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_v_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_v_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_w_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_w_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph_1, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph_2, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_init, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_phb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph0, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_php, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_pb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_al, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_alt, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_alb, 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu_1, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu_2, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mub, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu0, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%ht, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msft, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msfu, 3, 4, 1, 0, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msfv, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%sina, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%cosa, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%e, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%f, 3, 4, 1, 0, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
   CALL RSL_LITE_EXCH_PERIOD_X ( local_communicator_periodic , mytask, ntasks, ntasks_x, ntasks_y )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_u_1, 3, 4, 1, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_u_2, 3, 4, 1, 1, 'XZY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_v_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_v_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_w_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_w_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph_1, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph_2, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_t_init, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_phb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_ph0, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_php, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_pb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_al, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_alt, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_alb, 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu_1, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu_2, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mub, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%em_mu0, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%ht, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msft, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msfu, 3, 4, 1, 1, 'XY', 1, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%msfv, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%sina, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%cosa, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%e, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic, grid%f, 3, 4, 1, 1, 'XY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, 1  , 1  ,             &
ims, ime, jms, jme, 1  , 1  ,             &
ips, ipe, jps, jpe, 1  , 1                )
END IF
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/PERIOD_BDY_EM_MOIST.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/PERIOD_BDY_EM_MOIST.inc')
IF ( config_flags%periodic_x ) THEN
CALL RSL_LITE_INIT_PERIOD ( local_communicator_periodic, 3 , &
     0  &
   + num_moist   &
     , 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, kpe    )
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
   CALL RSL_LITE_EXCH_PERIOD_X ( local_communicator_periodic , mytask, ntasks, ntasks_x, ntasks_y )
DO itrace = PARAM_FIRST_SCALAR, num_moist
 CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic,moist ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
END IF
!ENDOFREGISTRYGENERATEDINCLUDE
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/PERIOD_BDY_EM_CHEM.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL wrf_debug(2,'calling inc/PERIOD_BDY_EM_CHEM.inc')
IF ( config_flags%periodic_x ) THEN
CALL RSL_LITE_INIT_PERIOD ( local_communicator_periodic, 3 , &
     0  &
   + num_chem   &
     , 0, 4, &
     0, 0, 4, &
     0, 0, 8, &
      0,  0, 4, &
      mytask, ntasks, ntasks_x, ntasks_y,   &
      ips, ipe, jps, jpe, kps, kpe    )
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 0, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
   CALL RSL_LITE_EXCH_PERIOD_X ( local_communicator_periodic , mytask, ntasks, ntasks_x, ntasks_y )
DO itrace = PARAM_FIRST_SCALAR, num_chem
 CALL RSL_LITE_PACK_PERIOD_X ( local_communicator_periodic,chem ( grid%sm31,grid%sm32,grid%sm33,itrace), 3, 4, 1, 1, 'XZY', 0, &
mytask, ntasks, ntasks_x, ntasks_y,       &
ids, ide, jds, jde, kds, kde,             &
ims, ime, jms, jme, kms, kme,             &
ips, ipe, jps, jpe, kps, kpe              )
ENDDO
END IF
!ENDOFREGISTRYGENERATEDINCLUDE

     CALL wrf_debug ( 100 , 'module_start: start_domain_rk: Returning' )

     RETURN

   END SUBROUTINE start_domain_em

