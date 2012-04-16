MODULE module_si_io_nmm

   USE module_optional_input

   IMPLICIT NONE


















  integer, parameter, public  :: i_byte  = selected_int_kind(1)      
  integer, parameter, public  :: i_short = selected_int_kind(4)      
  integer, parameter, public  :: i_long  = selected_int_kind(8)      
  integer, parameter, private :: llong_t = selected_int_kind(16)     
  integer, parameter, public  :: i_llong = max( llong_t, i_long )


  integer, parameter, public :: num_bytes_for_i_byte  = 1
  integer, parameter, public :: num_bytes_for_i_short = 2
  integer, parameter, public :: num_bytes_for_i_long  = 4
  integer, parameter, public :: num_bytes_for_i_llong = 8


  integer, parameter, private :: num_i_kinds = 4
  integer, parameter, dimension( num_i_kinds ), private :: integer_types = (/ &
       i_byte, i_short, i_long,  i_llong  /)
  integer, parameter, dimension( num_i_kinds ), private :: integer_byte_sizes = (/ &
       num_bytes_for_i_byte, num_bytes_for_i_short, &
       num_bytes_for_i_long, num_bytes_for_i_llong  /)



  integer, parameter, private :: default_integer = 2  
                                                      
                                                      
                                                      
  integer, parameter, public  :: i_kind = integer_types( default_integer )
  integer, parameter, public  :: num_bytes_for_i_kind = &
       integer_byte_sizes( default_integer )





  integer, parameter, public  :: r_single = selected_real_kind(6)  
  integer, parameter, public  :: r_double = selected_real_kind(15) 
  integer, parameter, private :: quad_t   = selected_real_kind(20) 
  integer, parameter, public  :: r_quad   = max( quad_t, r_double )


  integer, parameter, public :: num_bytes_for_r_single = 4
  integer, parameter, public :: num_bytes_for_r_double = 8
  integer, parameter, public :: num_bytes_for_r_quad   = 16


  integer, parameter, private :: num_r_kinds = 3
  integer, parameter, dimension( num_r_kinds ), private :: real_kinds = (/ &
       r_single, r_double, r_quad    /)
  integer, parameter, dimension( num_r_kinds ), private :: real_byte_sizes = (/ &
       num_bytes_for_r_single, num_bytes_for_r_double, &
       num_bytes_for_r_quad    /)



  integer, parameter, private :: default_real = 2  
                                                   


      

      REAL , DIMENSION(:,:,:) , ALLOCATABLE :: u_input , v_input , &
                                               q_input , t_input

      

      REAL , DIMENSION(:,:,:) , ALLOCATABLE :: landuse_frac_input , &
                                               soil_top_cat_input , &
                                               soil_bot_cat_input

      REAL, ALLOCATABLE:: htm_in(:,:,:),vtm_in(:,:,:)

      

      REAL , DIMENSION(:,:)   , ALLOCATABLE,save :: soilt010_input , soilt040_input , &
                                               soilt100_input , soilt200_input , &
                                               soilm010_input , soilm040_input , &
                                               soilm100_input , soilm200_input , &
                                               psfc_in,pmsl

      REAL , DIMENSION(:,:)   , ALLOCATABLE :: lat_wind, lon_wind 

      REAL , DIMENSION(:)     , ALLOCATABLE :: DETA_in, AETA_in, ETAX_in
      REAL , DIMENSION(:)     , ALLOCATABLE :: DETA1_in, AETA1_in, ETA1_in
      REAL , DIMENSION(:)     , ALLOCATABLE :: DETA2_in, AETA2_in, ETA2_in, DFL_in

      REAL , DIMENSION(:,:,:), ALLOCATABLE,save :: st_inputx , sm_inputx, sw_inputx

      

      REAL,DIMENSION(:,:),ALLOCATABLE :: dum2d
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: idum2d
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: dum3d

      LOGICAL , SAVE :: first_time_in = .TRUE.

      INTEGER :: flag_soilt010 , flag_soilt100 , flag_soilt200 , &
        	 flag_soilm010 , flag_soilm100 , flag_soilm200




      INTEGER, PARAMETER          :: var_maxdims = 5
      INTEGER, PARAMETER          :: max_staggers_xy_new = 4
      INTEGER, PARAMETER          :: max_staggers_xy_old = 3
      INTEGER, PARAMETER          :: max_staggers_z = 2
      INTEGER, PARAMETER          :: max_standard_lats = 4
      INTEGER, PARAMETER          :: max_standard_lons = 4  
      INTEGER, PARAMETER          :: max_fg_variables = 200
      INTEGER, PARAMETER          :: max_vertical_levels = 2000























      TYPE wrf_var_metadata
	 CHARACTER (LEN=8)         :: name 
	 CHARACTER (LEN=16)        :: units
	 CHARACTER (LEN=80)        :: description
	 INTEGER                   :: domain_id
	 INTEGER                   :: ndim
	 INTEGER                   :: dim_val (var_maxdims)
	 CHARACTER(LEN=4)          :: dim_desc (var_maxdims)
	 INTEGER                   :: start_index(var_maxdims)
	 INTEGER                   :: stop_index(var_maxdims)
	 INTEGER                   :: h_stagger_index
	 INTEGER                   :: v_stagger_index
	 CHARACTER(LEN=8)          :: array_order
	 CHARACTER(LEN=4)          :: field_type
	 CHARACTER(LEN=8)          :: field_source_prog
	 CHARACTER(LEN=80)         :: source_desc
	 CHARACTER(LEN=8)          :: field_time_type
	 INTEGER                   :: vt_date_start
	 REAL                      :: vt_time_start
	 INTEGER                   :: vt_date_stop
	 REAL                      :: vt_time_stop
      END TYPE wrf_var_metadata

      TYPE(wrf_var_metadata)  :: var_meta , var_info

      TYPE wrf_domain_metadata
	 INTEGER                   :: id
	 INTEGER                   :: parent_id
	 CHARACTER(LEN=8)          :: dyn_init_src
	 CHARACTER(LEN=8)          :: static_init_src 
	 INTEGER                   :: vt_date
	 REAL                      :: vt_time
	 INTEGER                   :: origin_parent_x
	 INTEGER                   :: origin_parent_y
	 INTEGER                   :: ratio_to_parent
	 REAL                      :: delta_x
	 REAL                      :: delta_y
	 REAL                      :: top_level
	 INTEGER                   :: origin_parent_z
	 REAL                      :: corner_lats_new(4,max_staggers_xy_new)
	 REAL                      :: corner_lons_new(4,max_staggers_xy_new)
	 REAL                      :: corner_lats_old(4,max_staggers_xy_old)
	 REAL                      :: corner_lons_old(4,max_staggers_xy_old)
	 INTEGER                   :: xdim
	 INTEGER                   :: ydim
	 INTEGER                   :: zdim
      END TYPE wrf_domain_metadata
      TYPE(wrf_domain_metadata) :: dom_meta

      TYPE wrf_global_metadata
	 CHARACTER(LEN=80)         :: simulation_name
	 CHARACTER(LEN=80)         :: user_desc
	 INTEGER                   :: si_version
	 INTEGER                   :: analysis_version  
	 INTEGER                   :: wrf_version
	 INTEGER                   :: post_version
	 CHARACTER(LEN=32)         :: map_projection
	 REAL                      :: moad_known_lat
	 REAL                      :: moad_known_lon
	 CHARACTER(LEN=8)          :: moad_known_loc
	 REAL                      :: moad_stand_lats(max_standard_lats)
	 REAL                      :: moad_stand_lons(max_standard_lons)
	 REAL                      :: moad_delta_x
	 REAL                      :: moad_delta_y
	 CHARACTER(LEN=4)          :: horiz_stagger_type
	 INTEGER                   :: num_stagger_xy
	 REAL                      :: stagger_dir_x_new(max_staggers_xy_new)
	 REAL                      :: stagger_dir_y_new(max_staggers_xy_new)
	 REAL                      :: stagger_dir_x_old(max_staggers_xy_old)
	 REAL                      :: stagger_dir_y_old(max_staggers_xy_old)
	 INTEGER                   :: num_stagger_z    
	 REAL                      :: stagger_dir_z(max_staggers_z)
	 CHARACTER(LEN=8)          :: vertical_coord
	 INTEGER                   :: num_domains
	 INTEGER                   :: init_date
	 REAL                      :: init_time
	 INTEGER                   :: end_date
	 REAL                      :: end_time
	 CHARACTER(LEN=4)          :: lu_source
	 INTEGER                   :: lu_water
	 INTEGER                   :: lu_ice  
      END TYPE wrf_global_metadata
      TYPE(wrf_global_metadata)   :: global_meta

CONTAINS

   SUBROUTINE read_si ( grid, file_date_string )

      USE module_soil_pre
      USE module_domain

      IMPLICIT NONE

      TYPE(domain) , INTENT(INOUT)  :: grid
      CHARACTER (LEN=19) , INTENT(IN) :: file_date_string

      INTEGER :: ids,ide,jds,jde,kds,kde           &
                ,ims,ime,jms,jme,kms,kme           &
                ,its,ite,jts,jte,kts,kte

      INTEGER :: i , j , k , loop, IMAX, JMAX

      REAL :: dummy

      CHARACTER (LEN= 8) :: dummy_char

      INTEGER :: ok , map_proj , ok_open
      REAL :: pt
      INTEGER :: num_veg_cat , num_soil_top_cat , num_soil_bot_cat

      SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
            kds = grid%sd31 ; kde = grid%ed31 ;
            ids = grid%sd32 ; ide = grid%ed32 ;
            jds = grid%sd33 ; jde = grid%ed33 ;

            kms = grid%sm31 ; kme = grid%em31 ;
            ims = grid%sm32 ; ime = grid%em32 ;
            jms = grid%sm33 ; jme = grid%em33 ;

            kts = grid%sp31 ; kte = grid%ep31 ; 
            its = grid%sp32 ; ite = grid%ep32 ; 
            jts = grid%sp33 ; jte = grid%ep33 ; 

         CASE ( DATA_ORDER_XYZ )
            ids = grid%sd31 ; ide = grid%ed31 ;
            jds = grid%sd32 ; jde = grid%ed32 ;
            kds = grid%sd33 ; kde = grid%ed33 ;

            ims = grid%sm31 ; ime = grid%em31 ;
            jms = grid%sm32 ; jme = grid%em32 ;
            kms = grid%sm33 ; kme = grid%em33 ;

            its = grid%sp31 ; ite = grid%ep31 ; 
            jts = grid%sp32 ; jte = grid%ep32 ; 
            kts = grid%sp33 ; kte = grid%ep33 ; 

         CASE ( DATA_ORDER_XZY )
            ids = grid%sd31 ; ide = grid%ed31 ;
            kds = grid%sd32 ; kde = grid%ed32 ;
            jds = grid%sd33 ; jde = grid%ed33 ;

            ims = grid%sm31 ; ime = grid%em31 ;
            kms = grid%sm32 ; kme = grid%em32 ;
            jms = grid%sm33 ; jme = grid%em33 ;

            its = grid%sp31 ; ite = grid%ep31 ; 
            kts = grid%sp32 ; kte = grid%ep32 ; 
            jts = grid%sp33 ; jte = grid%ep33 ; 

      END SELECT

      

      write(0,*) 'dum3d I allocs: ', ids,ide-1
      write(0,*) 'dum3d J allocs: ', jds,jde-1
      write(0,*) 'dum3d K allocs: ', kds,kde-1

      flag_st000010 = 0
      flag_st010040 = 0
      flag_st040100 = 0
      flag_st100200 = 0
      flag_sm000010 = 0 
      flag_sm010040 = 0
      flag_sm040100 = 0
      flag_sm100200 = 0
      flag_st010200 = 0
      flag_sm010200 = 0

      flag_soilt010 = 0
      flag_soilt040 = 0
      flag_soilt100 = 0
      flag_soilt200 = 0 
      flag_soilm010 = 0 
      flag_soilm040 = 0
      flag_soilm100 = 0
      flag_soilm200 = 0

      flag_sst      = 0
      flag_toposoil = 0

      

      num_st_levels_input = 0
      num_sm_levels_input = 0
      st_levels_input = -1
      sm_levels_input = -1

      

        write(6,*) 'enter read_si...first_time_in:: ', first_time_in

      IF ( first_time_in ) THEN

         CLOSE(12)
         OPEN ( FILE   = 'real_input_nm.global.metadata' , &
                UNIT   = 12                              , &
                STATUS = 'OLD'                           , &
                ACCESS = 'SEQUENTIAL'                    , &
                FORM   = 'UNFORMATTED'                   , &
                IOSTAT = ok_open                           )

         IF ( ok_open .NE. 0 ) THEN
            PRINT '(A)','You asked for WRF SI data, but no real_input_nm.global.metadata file exists.'
            STOP 'No_real_input_nm.global.metadata_exists'
         END IF

         READ(12) global_meta%simulation_name, global_meta%user_desc, &
                  global_meta%si_version, global_meta%analysis_version, &
                  global_meta%wrf_version, global_meta%post_version
   
         REWIND (12)

         IF      ( global_meta%si_version .EQ. 1 ) THEN
            READ(12) global_meta%simulation_name, global_meta%user_desc, &
                     global_meta%si_version, global_meta%analysis_version, &
                     global_meta%wrf_version, global_meta%post_version, &
                     global_meta%map_projection, global_meta%moad_known_lat, &
                     global_meta%moad_known_lon, global_meta%moad_known_loc, &
                     global_meta%moad_stand_lats, global_meta%moad_stand_lons, &
                     global_meta%moad_delta_x, global_meta%moad_delta_y, &
                     global_meta%horiz_stagger_type, global_meta%num_stagger_xy, &
                     global_meta%stagger_dir_x_old, global_meta%stagger_dir_y_old, &
                     global_meta%num_stagger_z, global_meta%stagger_dir_z, &
                     global_meta%vertical_coord, global_meta%num_domains, &
                     global_meta%init_date, global_meta%init_time, &
                     global_meta%end_date, global_meta%end_time
         ELSE IF ( global_meta%si_version .EQ. 2 ) THEN
            READ(12) global_meta%simulation_name, global_meta%user_desc, &
                     global_meta%si_version, global_meta%analysis_version, &
                     global_meta%wrf_version, global_meta%post_version, &
                     global_meta%map_projection, global_meta%moad_known_lat, &
                     global_meta%moad_known_lon, global_meta%moad_known_loc, &
                     global_meta%moad_stand_lats, global_meta%moad_stand_lons, &
                     global_meta%moad_delta_x, global_meta%moad_delta_y, &
                     global_meta%horiz_stagger_type, global_meta%num_stagger_xy, &
                     global_meta%stagger_dir_x_new, global_meta%stagger_dir_y_new, &
                     global_meta%num_stagger_z, global_meta%stagger_dir_z, &
                     global_meta%vertical_coord, global_meta%num_domains, &
                     global_meta%init_date, global_meta%init_time, &
                     global_meta%end_date, global_meta%end_time , &
                     global_meta%lu_source, global_meta%lu_water, global_meta%lu_ice
         END IF
         CLOSE (12)
   
         print *,'GLOBAL METADATA'
         print *,'global_meta%simulation_name', global_meta%simulation_name
         print *,'global_meta%user_desc', global_meta%user_desc
         print *,'global_meta%user_desc', global_meta%user_desc
         print *,'global_meta%si_version', global_meta%si_version
         print *,'global_meta%analysis_version', global_meta%analysis_version
         print *,'global_meta%wrf_version', global_meta%wrf_version
         print *,'global_meta%post_version', global_meta%post_version
         print *,'global_meta%map_projection', global_meta%map_projection
         print *,'global_meta%moad_known_lat', global_meta%moad_known_lat
         print *,'global_meta%moad_known_lon', global_meta%moad_known_lon
         print *,'global_meta%moad_known_loc', global_meta%moad_known_loc
         print *,'global_meta%moad_stand_lats', global_meta%moad_stand_lats
         print *,'global_meta%moad_stand_lons', global_meta%moad_stand_lons
         print *,'global_meta%moad_delta_x', global_meta%moad_delta_x
         print *,'global_meta%moad_delta_y', global_meta%moad_delta_y
         print *,'global_meta%horiz_stagger_type', global_meta%horiz_stagger_type
         print *,'global_meta%num_stagger_xy', global_meta%num_stagger_xy
         IF      ( global_meta%si_version .EQ. 1 ) THEN
            print *,'global_meta%stagger_dir_x', global_meta%stagger_dir_x_old
            print *,'global_meta%stagger_dir_y', global_meta%stagger_dir_y_old
         ELSE IF ( global_meta%si_version .EQ. 2 ) THEN
            print *,'global_meta%stagger_dir_x', global_meta%stagger_dir_x_new
            print *,'global_meta%stagger_dir_y', global_meta%stagger_dir_y_new
         END IF
         print *,'global_meta%num_stagger_z', global_meta%num_stagger_z
         print *,'global_meta%stagger_dir_z', global_meta%stagger_dir_z
         print *,'global_meta%vertical_coord', global_meta%vertical_coord
         print *,'global_meta%num_domains', global_meta%num_domains
         print *,'global_meta%init_date', global_meta%init_date
         print *,'global_meta%init_time', global_meta%init_time
         print *,'global_meta%end_date', global_meta%end_date
         print *,'global_meta%end_time', global_meta%end_time
         IF ( global_meta%si_version .EQ. 2 ) THEN
            print *,'global_meta%lu_source', global_meta%lu_source
            print *,'global_meta%lu_water', global_meta%lu_water
            print *,'global_meta%lu_ice', global_meta%lu_ice
         END IF
         print *,' '

         

        IF (.NOT. ALLOCATED (DETA_in)) ALLOCATE(DETA_in(kds:kde-1))
        IF (.NOT. ALLOCATED (AETA_in)) ALLOCATE(AETA_in(kds:kde-1))
        IF (.NOT. ALLOCATED (ETAX_in)) ALLOCATE(ETAX_in(kds:kde))

        IF (.NOT. ALLOCATED (DETA1_in)) ALLOCATE(DETA1_in(kds:kde-1))
        IF (.NOT. ALLOCATED (AETA1_in)) ALLOCATE(AETA1_in(kds:kde-1))
        IF (.NOT. ALLOCATED (ETA1_in))  ALLOCATE(ETA1_in(kds:kde))

        IF (.NOT. ALLOCATED (DETA2_in)) ALLOCATE(DETA2_in(kds:kde-1))
        IF (.NOT. ALLOCATED (AETA2_in)) ALLOCATE(AETA2_in(kds:kde-1))
        IF (.NOT. ALLOCATED (ETA2_in)) ALLOCATE(ETA2_in(kds:kde))

        IF (.NOT. ALLOCATED (DFL_in)) ALLOCATE(DFL_in(kds:kde))

         

        IF (.NOT. ALLOCATED (u_input)  ) ALLOCATE ( u_input(its:ite,jts:jte,kts:kte) )
        IF (.NOT. ALLOCATED (v_input)  ) ALLOCATE ( v_input(its:ite,jts:jte,kts:kte) )
        IF (.NOT. ALLOCATED (q_input)  ) ALLOCATE ( q_input(its:ite,jts:jte,kts:kte) )
        IF (.NOT. ALLOCATED (t_input)  ) ALLOCATE ( t_input(its:ite,jts:jte,kts:kte) )
        IF (.NOT. ALLOCATED (htm_in)  ) ALLOCATE ( htm_in(its:ite,jts:jte,kts:kte) )
        IF (.NOT. ALLOCATED (vtm_in)  ) ALLOCATE ( vtm_in(its:ite,jts:jte,kts:kte) )

        

        IF (.NOT. ALLOCATED (pmsl)              ) ALLOCATE ( pmsl(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (psfc_in)           ) ALLOCATE ( psfc_in(its:ite,jts:jte) )

        

        

        IF (.NOT. ALLOCATED (st_inputx)) ALLOCATE (st_inputx(its:ite,jts:jte,num_st_levels_alloc))
        IF (.NOT. ALLOCATED (sm_inputx)) ALLOCATE (sm_inputx(its:ite,jts:jte,num_st_levels_alloc))
        IF (.NOT. ALLOCATED (sw_inputx)) ALLOCATE (sw_inputx(its:ite,jts:jte,num_st_levels_alloc))

        IF (.NOT. ALLOCATED (soilt010_input)    ) ALLOCATE ( soilt010_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilt040_input)    ) ALLOCATE ( soilt040_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilt100_input)    ) ALLOCATE ( soilt100_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilt200_input)    ) ALLOCATE ( soilt200_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilm010_input)    ) ALLOCATE ( soilm010_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilm040_input)    ) ALLOCATE ( soilm040_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilm100_input)    ) ALLOCATE ( soilm100_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilm200_input)    ) ALLOCATE ( soilm200_input(its:ite,jts:jte) )

        IF (.NOT. ALLOCATED (lat_wind)          ) ALLOCATE (lat_wind(its:ite,jts:jte))
        IF (.NOT. ALLOCATED (lon_wind)          ) ALLOCATE (lon_wind(its:ite,jts:jte))

        
        IF (.NOT. ALLOCATED (dum2d)             ) ALLOCATE (dum2d(IDS:IDE-1,JDS:JDE-1))
        IF (.NOT. ALLOCATED (idum2d)            ) ALLOCATE (idum2d(IDS:IDE-1,JDS:JDE-1))
        IF (.NOT. ALLOCATED (dum3d)             ) ALLOCATE (dum3d(IDS:IDE-1,JDS:JDE-1,KDS:KDE-1))


      END IF

      CLOSE(13)

      write(6,*) 'file_date_string: ', file_date_string
      write(6,*) 'opening real_input_nm.d01.'//file_date_string//' as unit 13'
      OPEN ( FILE   = 'real_input_nm.d01.'//file_date_string , &
             UNIT   = 13                                     , &
             STATUS = 'OLD'                                  , &
             ACCESS = 'SEQUENTIAL'                           , &
             FORM   = 'UNFORMATTED'                            )

      IF      ( global_meta%si_version .EQ. 1 ) THEN
         READ (13) dom_meta%id,dom_meta%parent_id,dom_meta%dyn_init_src,&
                   dom_meta%static_init_src, dom_meta%vt_date, dom_meta%vt_time, &
                   dom_meta%origin_parent_x, dom_meta%origin_parent_y, &
                   dom_meta%ratio_to_parent, dom_meta%delta_x, dom_meta%delta_y, &
                   dom_meta%top_level, dom_meta%origin_parent_z, &
                   dom_meta%corner_lats_old, dom_meta%corner_lons_old, dom_meta%xdim, &
                   dom_meta%ydim, dom_meta%zdim
      ELSE IF ( global_meta%si_version .EQ. 2 ) THEN
         READ (13) dom_meta%id,dom_meta%parent_id,dom_meta%dyn_init_src,&
                   dom_meta%static_init_src, dom_meta%vt_date, dom_meta%vt_time, &
                   dom_meta%origin_parent_x, dom_meta%origin_parent_y, &
                   dom_meta%ratio_to_parent, dom_meta%delta_x, dom_meta%delta_y, &
                   dom_meta%top_level, dom_meta%origin_parent_z, &
                   dom_meta%corner_lats_new, dom_meta%corner_lons_new, dom_meta%xdim, &
                   dom_meta%ydim, dom_meta%zdim
      END IF

      print *,'DOMAIN METADATA'
      print *,'dom_meta%id=', dom_meta%id
      print *,'dom_meta%parent_id=', dom_meta%parent_id
      print *,'dom_meta%dyn_init_src=', dom_meta%dyn_init_src
      print *,'dom_meta%static_init_src=', dom_meta%static_init_src
      print *,'dom_meta%vt_date=', dom_meta%vt_date
      print *,'dom_meta%vt_time=', dom_meta%vt_time
      print *,'dom_meta%origin_parent_x=', dom_meta%origin_parent_x
      print *,'dom_meta%origin_parent_y=', dom_meta%origin_parent_y
      print *,'dom_meta%ratio_to_parent=', dom_meta%ratio_to_parent
      print *,'dom_meta%delta_x=', dom_meta%delta_x
      print *,'dom_meta%delta_y=', dom_meta%delta_y
      print *,'dom_meta%top_level=', dom_meta%top_level
      print *,'dom_meta%origin_parent_z=', dom_meta%origin_parent_z
      IF      ( global_meta%si_version .EQ. 1 ) THEN
         print *,'dom_meta%corner_lats=', dom_meta%corner_lats_old
         print *,'dom_meta%corner_lons=', dom_meta%corner_lons_old
      ELSE IF ( global_meta%si_version .EQ. 2 ) THEN
         print *,'dom_meta%corner_lats=', dom_meta%corner_lats_new
         print *,'dom_meta%corner_lons=', dom_meta%corner_lons_new
      END IF
      print *,'dom_meta%xdim=', dom_meta%xdim
      print *,'dom_meta%ydim=', dom_meta%ydim
      print *,'dom_meta%zdim=', dom_meta%zdim
      print *,' '

      
    




      IF (  abs(dom_meta%xdim - (ide-1)) .gt. 1 &
       .OR. abs(dom_meta%ydim - (jde-1)) .gt. 1 &
       .OR. abs(dom_meta%zdim - (kde-1)) .gt. 1) THEN
         PRINT '(A)','Namelist does not match the input data.'
         PRINT '(A,3I5,A)','Namelist dimensions =',ide-1,jde-1,kde-1,'.'
         PRINT '(A,3I5,A)','Input data dimensions =',dom_meta%xdim,dom_meta%ydim,dom_meta%zdim,'.'
         STOP 'Wrong_data_size'
      END IF

      

      IF        ( global_meta%si_version .EQ. 1 ) THEN
         CALL nl_set_cen_lat ( grid%id , ( dom_meta%corner_lats_old(1,1) + dom_meta%corner_lats_old(2,1) +        &
                                        dom_meta%corner_lats_old(3,1) + dom_meta%corner_lats_old(4,1) ) * 0.25 ) 
      ELSE IF ( ( global_meta%si_version .EQ. 2 ) .AND. ( global_meta%moad_known_loc(1:6) .EQ. 'CENTER' ) ) THEN
         CALL nl_set_cen_lat ( grid%id , global_meta%moad_known_lat )
      ELSE IF   ( global_meta%si_version .EQ. 2 ) THEN
         CALL nl_set_cen_lat ( grid%id , ( dom_meta%corner_lats_new(1,1) + dom_meta%corner_lats_new(2,1) +        &
                                        dom_meta%corner_lats_new(3,1) + dom_meta%corner_lats_new(4,1) ) * 0.25 ) 
      END IF




      CALL nl_set_cen_lon ( grid%id , global_meta%moad_stand_lons(1) )

      write(6,*) 'set_cen_lat... global_meta%moad_stand_lats(1): ', global_meta%moad_stand_lats(1)
      CALL nl_set_cen_lat ( grid%id , global_meta%moad_stand_lats(1) )

      CALL nl_set_truelat1 ( grid%id , global_meta%moad_stand_lats(1) )
      CALL nl_set_truelat2 ( grid%id , global_meta%moad_stand_lats(2) )

      pt = dom_meta%top_level

      IF      ( global_meta%map_projection(1:17) .EQ. 'LAMBERT CONFORMAL'   ) THEN
         map_proj = 1
      ELSE IF ( global_meta%map_projection(1:19) .EQ. 'POLAR STEREOGRAPHIC' ) THEN
         map_proj = 2
      ELSE IF ( global_meta%map_projection(1: 8) .EQ. 'MERCATOR'            ) THEN
         map_proj = 3
      ELSE IF ( global_meta%map_projection(1:14) .EQ. 'ROTATED LATLON' ) THEN
         map_proj = 203 
      ELSE
         PRINT '(A,A,A)','Undefined map projection: ',TRIM(global_meta%map_projection(1:20)),'.'
         STOP 'Undefined_map_proj_si'
      END IF
      CALL nl_set_map_proj ( grid%id , map_proj ) 
     
      write(0,*) 'global_meta%si_version: ', global_meta%si_version
      write(0,*) 'global_meta%lu_source: ', global_meta%lu_source
      write(0,*) 'global_meta%lu_water: ', global_meta%lu_water
      IF      ( global_meta%si_version .EQ. 1 ) THEN
         CALL nl_set_mminlu (grid%id, 'USGS' )
         CALL nl_set_iswater (grid%id, 16 )
      ELSE IF ( global_meta%si_version .EQ. 2 ) THEN
         CALL nl_set_mminlu ( grid%id, global_meta%lu_source )
         CALL nl_set_iswater (grid%id, global_meta%lu_water )
         CALL nl_set_isice (grid%id, global_meta%lu_ice )
      END IF

      CALL nl_set_gmt (grid%id, dom_meta%vt_time / 3600. )
      CALL nl_set_julyr (grid%id, dom_meta%vt_date / 1000 )
      CALL nl_set_julday (grid%id, dom_meta%vt_date - ( dom_meta%vt_date / 1000 ) * 1000 )

      write(6,*) 'start reading from unit 13'
      read_all_the_data : DO

         READ (13,IOSTAT=OK) var_info%name, var_info%units, &
                             var_info%description, var_info%domain_id, var_info%ndim, &
                             var_info%dim_val, var_info%dim_desc, var_info%start_index, &
                             var_info%stop_index, var_info%h_stagger_index, var_info%v_stagger_index,&
                             var_info%array_order, var_info%field_type, var_info%field_source_prog, &
                             var_info%source_desc, var_info%field_time_type, var_info%vt_date_start, &
                             var_info%vt_time_start, var_info%vt_date_stop, var_info%vt_time_stop

         IF ( OK .NE. 0 ) THEN
            PRINT '(A,A,A)','End of file found for real_input_nm.d01.',file_date_string,'.'
            EXIT read_all_the_data
         END IF


         PRINT '(A,A)','var_info%name=', var_info%name 




















        JMAX=min(JDE-1,JTE)
        IMAX=min(IDE-1,ITE)
         

         write(0,*)' read_si var_info%name=',var_info%name(1:8)

         IF      ( var_info%name(1:8) .EQ. 'T       ' ) THEN
            READ (13) dum3d
            do k=kts,kte-1
            do j=jts,JMAX
            do i=its,IMAX
              t_input(i,j,k)=dum3d(i,j,k)
            enddo
            enddo
            enddo

         ELSE IF      ( var_info%name(1:8) .EQ. 'U       ' ) THEN
            READ (13) dum3d
            do k=kts,kte-1
            do j=jts,JMAX
            do i=its,IMAX
              u_input(i,j,k)=dum3d(i,j,k)
            enddo
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'V       ' ) THEN
            READ (13) dum3d
            do k=kts,kte-1
            do j=jts,JMAX
            do i=its,IMAX
              v_input(i,j,k)=dum3d(i,j,k)
            enddo
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'Q      ' ) THEN
            READ (13) dum3d
            do k=kts,kte-1
            do j=jts,JMAX
            do i=its,IMAX
              q_input(i,j,k)=dum3d(i,j,k)
            enddo
            enddo
            enddo

         

         ELSE IF ( var_info%name(1:8) .EQ. 'LANDUSEF' ) THEN
            IF ( ( first_time_in ) .AND. ( .NOT. ALLOCATED ( landuse_frac_input) ) ) THEN
               ALLOCATE (landuse_frac_input(its:ite,jts:jte,var_info%dim_val(3)) )
            END IF
            READ (13) (((dum3d(i,j,k),i=ids,ide-1),j=jds,jde-1),k=1,var_info%dim_val(3))
            do k=1,var_info%dim_val(3)
            do j=jts,JMAX
            do i=its,IMAX
              landuse_frac_input(i,j,k)=dum3d(i,j,k)
            enddo
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'SOILCTOP' ) THEN
            IF ( ( first_time_in ) .AND. ( .NOT. ALLOCATED ( soil_top_cat_input) ) ) THEN
               ALLOCATE (soil_top_cat_input(its:ite,jts:jte,var_info%dim_val(3)) )
            END IF
            READ (13) (((dum3d(i,j,k),i=ids,ide-1),j=jds,jde-1),k=1,var_info%dim_val(3))
            do k=1,var_info%dim_val(3)
            do j=jts,JMAX
            do i=its,IMAX
              soil_top_cat_input(i,j,k)=dum3d(i,j,k)
            enddo
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'SOILCBOT' ) THEN
            IF ( ( first_time_in ) .AND. ( .NOT. ALLOCATED ( soil_bot_cat_input) ) ) THEN
               ALLOCATE (soil_bot_cat_input(its:ite,jts:jte,var_info%dim_val(3)) )
            END IF
            READ (13) (((dum3d(i,j,k),i=ids,ide-1),j=jds,jde-1),k=1,var_info%dim_val(3))
            do k=1,var_info%dim_val(3)
            do j=jts,JMAX
            do i=its,IMAX
              soil_bot_cat_input(i,j,k)=dum3d(i,j,k)
            enddo
            enddo
            enddo

         

         ELSE IF ( var_info%name(1:8) .EQ. 'PD      ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%pd(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'PSFC    ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              psfc_in(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'PMSL    ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              pmsl(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'PDTOP   ' ) THEN
            READ (13) grid%pdtop

         ELSE IF ( var_info%name(1:8) .EQ. 'PT      ' ) THEN
            READ (13) grid%pt

         

        ELSE IF ( var_info%name(1:8) .eq. 'GLAT    ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%glat(i,j)=dum2d(i,j)
            enddo
            enddo
        ELSE IF ( var_info%name(1:8) .eq. 'GLON    ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%glon(i,j)=dum2d(i,j)
            enddo
            enddo
        ELSE IF ( var_info%name(1:8) .eq. 'LAT_V   ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              lat_wind(i,j)=dum2d(i,j)
            enddo
            enddo
        ELSE IF ( var_info%name(1:8) .eq. 'LON_V   ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              lon_wind(i,j)=dum2d(i,j)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'ST000010' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%st000010(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_st000010 = 1
            num_st_levels_input = num_st_levels_input + 1
            st_levels_input(num_st_levels_input) = char2int2(var_info%name(3:8))
            do j=jts,JMAX
            do i=its,IMAX
              st_inputx(I,J,num_st_levels_input + 1) = grid%st000010(i,j)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'ST010040' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%st010040(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_st010040 = 1
            num_st_levels_input = num_st_levels_input + 1
            st_levels_input(num_st_levels_input) = char2int2(var_info%name(3:8))
            do j=jts,JMAX
            do i=its,IMAX
              st_inputx(I,J,num_st_levels_input + 1) = grid%st010040(i,j)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'ST040100' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%st040100(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_st040100 = 1
            num_st_levels_input = num_st_levels_input + 1
            st_levels_input(num_st_levels_input) = char2int2(var_info%name(3:8))
            do j=jts,JMAX
            do i=its,IMAX
              st_inputx(I,J,num_st_levels_input + 1) = grid%st040100(i,j)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'ST100200' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%st100200(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_st100200 = 1
            num_st_levels_input = num_st_levels_input + 1
            st_levels_input(num_st_levels_input) = char2int2(var_info%name(3:8))
            do j=jts,JMAX
            do i=its,IMAX
              st_inputx(I,J,num_st_levels_input + 1) = grid%st100200(i,j)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'ST010200' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%st010200(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_st010200 = 1
            num_st_levels_input = num_st_levels_input + 1
            st_levels_input(num_st_levels_input) = char2int2(var_info%name(3:8))
            do j=jts,JMAX
            do i=its,IMAX
              st_inputx(I,J,num_st_levels_input + 1) = grid%st010200(i,j)
            enddo
            enddo

        ELSE IF ( var_info%name(1:8) .EQ. 'SM000010' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%sm000010(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_sm000010 = 1
            num_sm_levels_input = num_sm_levels_input + 1
            sm_levels_input(num_sm_levels_input) = char2int2(var_info%name(3:8))
            do j=jts,JMAX
            do i=its,IMAX
              sm_inputx(I,J,num_sm_levels_input + 1) = grid%sm000010(i,j)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'SM010040' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%sm010040(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_sm010040 = 1
            num_sm_levels_input = num_sm_levels_input + 1
            sm_levels_input(num_sm_levels_input) = char2int2(var_info%name(3:8))
            do j=jts,JMAX
            do i=its,IMAX
              sm_inputx(I,J,num_sm_levels_input + 1) = grid%sm010040(i,j)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'SM040100' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%sm040100(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_sm040100 = 1
            num_sm_levels_input = num_sm_levels_input + 1
            sm_levels_input(num_sm_levels_input) = char2int2(var_info%name(3:8))
            do j=jts,JMAX
            do i=its,IMAX
              sm_inputx(I,J,num_sm_levels_input + 1) = grid%sm040100(i,j)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'SM100200' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%sm100200(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_sm100200 = 1
            num_sm_levels_input = num_sm_levels_input + 1
            sm_levels_input(num_sm_levels_input) = char2int2(var_info%name(3:8))
            do j=jts,JMAX
            do i=its,IMAX
              sm_inputx(I,J,num_sm_levels_input + 1) = grid%sm100200(i,j)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'SM010200' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%sm010200(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_sm010200 = 1
            num_sm_levels_input = num_sm_levels_input + 1
            sm_levels_input(num_sm_levels_input) = char2int2(var_info%name(3:8))
            do j=jts,JMAX
            do i=its,IMAX
               sm_inputx(I,J,num_sm_levels_input + 1) = grid%sm010200(i,j)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'SOILT010' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              soilt010_input(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_soilt010 = 1
            num_st_levels_input = num_st_levels_input + 1
            st_levels_input(num_st_levels_input) = char2int1(var_info%name(6:8))

            do j=jts,JMAX
            do i=its,IMAX
              st_inputx(I,J,num_st_levels_input + 1) = soilt010_input(I,J)
            enddo
            enddo
            write(6,*) 'num_st_levels_input=',num_st_levels_input
         ELSE IF ( var_info%name(1:8) .EQ. 'SOILT040' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              soilt040_input(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_soilt040 = 1
            num_st_levels_input = num_st_levels_input + 1
            st_levels_input(num_st_levels_input) = char2int1(var_info%name(6:8))

            do j=jts,JMAX
            do i=its,IMAX
              st_inputx(I,J,num_st_levels_input + 1) = soilt040_input(I,J)
            enddo
            enddo
            write(6,*) 'num_st_levels_input=',num_st_levels_input
         ELSE IF ( var_info%name(1:8) .EQ. 'SOILT100' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              soilt100_input(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_soilt100 = 1
            num_st_levels_input = num_st_levels_input + 1
            st_levels_input(num_st_levels_input) = char2int1(var_info%name(6:8))

            do j=jts,JMAX
            do i=its,IMAX
              st_inputx(I,J,num_st_levels_input + 1) = soilt100_input(I,J)
            enddo
            enddo
            write(6,*) 'num_st_levels_input=',num_st_levels_input
        ELSE IF ( var_info%name(1:8) .EQ. 'SOILT200' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              soilt200_input(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_soilt200 = 1
            num_st_levels_input = num_st_levels_input + 1
            st_levels_input(num_st_levels_input) = char2int1(var_info%name(6:8))

            do j=jts,JMAX
            do i=its,IMAX
              st_inputx(I,J,num_st_levels_input + 1) = soilt200_input(I,J)
            enddo
            enddo
            write(6,*) 'num_st_levels_input=',num_st_levels_input
         ELSE IF ( var_info%name(1:8) .EQ. 'SOILM010' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              soilm010_input(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_soilm010 = 1
            num_sm_levels_input = num_sm_levels_input + 1
            sm_levels_input(num_sm_levels_input) = char2int1(var_info%name(6:8))

            do j=jts,JMAX
            do i=its,IMAX
              sm_inputx(I,J,num_sm_levels_input + 1) = soilm010_input(I,J)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'SOILM040' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              soilm040_input(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_soilm040 = 1
            num_sm_levels_input = num_sm_levels_input + 1
            sm_levels_input(num_sm_levels_input) = char2int1(var_info%name(6:8))

            do j=jts,JMAX
            do i=its,IMAX
              sm_inputx(I,J,num_sm_levels_input + 1) = soilm040_input(I,J)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'SOILM100' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              soilm100_input(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_soilm100 = 1
            num_sm_levels_input = num_sm_levels_input + 1
            sm_levels_input(num_sm_levels_input) = char2int1(var_info%name(6:8))

            do j=jts,JMAX
            do i=its,IMAX
              sm_inputx(I,J,num_sm_levels_input + 1) = soilm100_input(I,J)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'SOILM200' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              soilm200_input(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_soilm200 = 1
            num_sm_levels_input = num_sm_levels_input + 1
            sm_levels_input(num_sm_levels_input) = char2int1(var_info%name(6:8))

            do j=jts,JMAX
            do i=its,IMAX
              sm_inputx(I,J,num_sm_levels_input + 1) = soilm200_input(I,J)
            enddo
            enddo

         ELSE IF ( var_info%name(1:8) .EQ. 'SEAICE  ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%xice(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'WEASD   ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%weasd(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'CANWAT  ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%canwat(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'LANDMASK' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%landmask(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'SKINTEMP' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%nmm_tsk(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'TGROUND ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
             grid%tg(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'SOILTB  ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
             grid%soiltb(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'SST     ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
               grid%sst(i,j)=dum2d(i,j)
            enddo
            enddo
            flag_sst = 1
         ELSE IF ( var_info%name(1:8) .EQ. 'GREENFRC' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%vegfrc(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'ISLOPE  ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%islope(i,j)=nint(dum2d(i,j))
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'GREENMAX' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%greenmax(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'GREENMIN' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
               grid%greenmin(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'FIS     ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%fis(i,j)=dum2d(i,j)
            enddo
            enddo
        ELSE IF ( var_info%name(1:8) .EQ. 'Z0      ' ) THEN

            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%z0(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'CMC     ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%cmc(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'HTM     ' ) THEN
            READ (13) dum3d
            do k=kts,kte-1
            do j=jts,JMAX
            do i=its,IMAX
              htm_in(i,j,k)=dum3d(i,j,k)
            enddo
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'VTM     ' ) THEN
            READ (13) dum3d
            do k=kts,kte-1
            do j=jts,JMAX
            do i=its,IMAX
              vtm_in(i,j,k)=dum3d(i,j,k)
            enddo
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'SM      ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%sm(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'ALBASE  ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%albase(i,j)=dum2d(i,j)
            enddo
            enddo
         ELSE IF ( var_info%name(1:8) .EQ. 'MXSNAL  ' ) THEN
            READ (13) dum2d
            do j=jts,JMAX
            do i=its,IMAX
              grid%mxsnal(i,j)=dum2d(i,j)
            enddo
            enddo

         

          ELSE IF ( var_info%name(1:8) .EQ. 'DETA    ' ) THEN
             READ(13) DETA_in
          ELSE IF ( var_info%name(1:8) .EQ. 'DETA1   ' ) THEN
             READ(13) DETA1_in
          ELSE IF ( var_info%name(1:8) .EQ. 'DETA2   ' ) THEN
             READ(13) DETA2_in
          ELSE IF ( var_info%name(1:8) .EQ. 'ETAX    ' ) THEN
             READ(13) ETAX_in
          ELSE IF ( var_info%name(1:8) .EQ. 'ETA1    ' ) THEN
             READ(13) ETA1_in
          ELSE IF ( var_info%name(1:8) .EQ. 'ETA2    ' ) THEN
             READ(13) ETA2_in
          ELSE IF ( var_info%name(1:8) .EQ. 'AETA    ' ) THEN
             READ(13) AETA_in
          ELSE IF ( var_info%name(1:8) .EQ. 'AETA1   ' ) THEN
             READ(13) AETA1_in
          ELSE IF ( var_info%name(1:8) .EQ. 'AETA2   ' ) THEN
             READ(13) AETA2_in
          ELSE IF ( var_info%name(1:8) .EQ. 'DFL     ' ) THEN
             READ(13) DFL_in






         

         ELSE IF ( var_info%name(1:8) .EQ. 'ZETAFULL' ) THEN
            PRINT '(A)','Oops, you put in the height data.'
            STOP 'this_is_mass_not_height'
 

         

         ELSE
print *,'------------------> skipping ', var_info%name(1:8)
            READ (13) dummy
         END IF

      END DO read_all_the_data

      CLOSE (13)

      first_time_in = .FALSE.


        sw_inputx=0.


      do k=kts,kte-1
      do j=jts,JMAX
      do i=its,IMAX
        grid%U(I,J,K)=U_input(I,J,K)
        grid%V(I,J,K)=V_input(I,J,K)
        grid%T(I,J,K)=T_input(I,J,K)
        grid%Q(I,J,K)=Q_input(I,J,K)
      enddo
      enddo
      enddo

      write(0,*) 'size sw_input: ', size(sw_input,dim=1),size(sw_input,dim=2),size(sw_input,dim=3)
      write(0,*) 'size sw_inputx: ', size(sw_inputx,dim=1),size(sw_inputx,dim=2),size(sw_inputx,dim=3)
      sw_input=0.

        write(0,*) 'maxval st_inputx(1): ', maxval(st_input(:,:,1))
        write(0,*) 'maxval st_inputx(2): ', maxval(st_input(:,:,2))
        write(0,*) 'maxval st_inputx(3): ', maxval(st_input(:,:,3))
        write(0,*) 'maxval st_inputx(4): ', maxval(st_input(:,:,4))


        do J=JTS,min(JDE-1,JTE)
         do K=1,num_st_levels_alloc
          do I=ITS,min(IDE-1,ITE)
             st_input(I,K,J)=st_inputx(I,J,K)
             sm_input(I,K,J)=sm_inputx(I,J,K)
             sw_input(I,K,J)=sw_inputx(I,J,K)
          enddo
         enddo
        enddo

        write(0,*) 'maxval st_input(1): ', maxval(st_input(:,1,:))
        write(0,*) 'maxval st_input(2): ', maxval(st_input(:,2,:))
        write(0,*) 'maxval st_input(3): ', maxval(st_input(:,3,:))
        write(0,*) 'maxval st_input(4): ', maxval(st_input(:,4,:))


         num_veg_cat      = SIZE ( grid%landusef , DIM=2 )
         num_soil_top_cat = SIZE ( grid%soilctop , DIM=2 )
         num_soil_bot_cat = SIZE ( grid%soilcbot , DIM=2 )

        do J=JTS,min(JDE-1,JTE)
         do K=1,num_soil_top_cat
          do I=ITS,min(IDE-1,ITE)
          grid%SOILCTOP(I,K,J)=soil_top_cat_input(I,J,K)
          grid%SOILCTOP_gc(I,J,K)=soil_top_cat_input(I,J,K)
          enddo
         enddo
        enddo

        do J=JTS,min(JDE-1,JTE)
         do K=1,num_soil_bot_cat
          do I=ITS,min(IDE-1,ITE)
          grid%SOILCBOT(I,K,J)=soil_bot_cat_input(I,J,K)
          grid%SOILCBOT_gc(I,J,K)=soil_bot_cat_input(I,J,K)
          enddo
         enddo
        enddo

        do J=JTS,min(JDE-1,JTE)
         do K=1,num_veg_cat
          do I=ITS,min(IDE-1,ITE)
          grid%LANDUSEF(I,K,J)=landuse_frac_input(I,J,K)
          grid%LANDUSEF_gc(I,J,K)=landuse_frac_input(I,J,K)
          enddo
         enddo
        enddo


      do K=KDS,KDE
        grid%ETAX(K)=ETAX_in(KDE-K+1)
        grid%ETA1(K)=ETA1_in(KDE-K+1)
        grid%ETA2(K)=ETA2_in(KDE-K+1)
        grid%DFL(K)=DFL_in(KDE-K+1)
      enddo

      do K=KDS,KDE-1
        grid%DETA(K)=DETA_in(KDE-K)
        grid%DETA1(K)=DETA1_in(KDE-K)
        grid%DETA2(K)=DETA2_in(KDE-K)
        grid%AETA(K)=AETA_in(KDE-K)
        grid%AETA1(K)=AETA1_in(KDE-K)
        grid%AETA2(K)=AETA2_in(KDE-K)
      enddo

   END SUBROUTINE read_si








   SUBROUTINE read_wps ( grid, filename, file_date_string, num_metgrid_levels )

      USE module_soil_pre
      USE module_domain

      INCLUDE "mpif.h"

      IMPLICIT NONE

      TYPE(domain) , INTENT(INOUT)  :: grid
      CHARACTER (LEN=19) , INTENT(IN) :: file_date_string
      CHARACTER (LEN=19)              :: VarName
      CHARACTER (LEN=150)             :: chartemp
      CHARACTER (*) , INTENT(IN) :: filename

      INTEGER :: ids,ide,jds,jde,kds,kde           &
                ,ims,ime,jms,jme,kms,kme           &
                ,its,ite,jts,jte,kts,kte

      INTEGER :: i , j , k , loop, IMAX, JMAX
      INTEGER :: DATAHANDLE, num_metgrid_levels
      INTEGER :: Sysdepinfo, Status, N
      INTEGER :: istatus,ioutcount,iret,index,ierr

      integer :: nrecs,iunit, L,hor_size


      character*132, allocatable,save :: datestr_all(:)
      character*132, allocatable,save :: varname_all(:)
      integer, allocatable,save       :: domainend_all(:,:)
      integer, allocatable,save       :: start_block(:)
      integer, allocatable,save       :: end_block(:)
      integer, allocatable,save       :: start_byte(:)
      integer, allocatable,save       :: end_byte(:)
      integer(kind=i_llong), allocatable, save           :: file_offset(:)


      REAL :: dummy,tmp,garb
      REAL, ALLOCATABLE,save:: dumdata(:,:,:)

      CHARACTER (LEN= 8) :: dummy_char

      INTEGER :: ok , map_proj , ok_open, igarb
      REAL :: pt
      INTEGER :: num_veg_cat , num_soil_top_cat , num_soil_bot_cat

      SELECT CASE ( model_data_order )
         CASE ( DATA_ORDER_ZXY )
            kds = grid%sd31 ; kde = grid%ed31 ;
            ids = grid%sd32 ; ide = grid%ed32 ;
            jds = grid%sd33 ; jde = grid%ed33 ;

            kms = grid%sm31 ; kme = grid%em31 ;
            ims = grid%sm32 ; ime = grid%em32 ;
            jms = grid%sm33 ; jme = grid%em33 ;

            kts = grid%sp31 ; kte = grid%ep31 ; 
            its = grid%sp32 ; ite = grid%ep32 ; 
            jts = grid%sp33 ; jte = grid%ep33 ; 

         CASE ( DATA_ORDER_XYZ )
            ids = grid%sd31 ; ide = grid%ed31 ;
            jds = grid%sd32 ; jde = grid%ed32 ;
            kds = grid%sd33 ; kde = grid%ed33 ;

            ims = grid%sm31 ; ime = grid%em31 ;
            jms = grid%sm32 ; jme = grid%em32 ;
            kms = grid%sm33 ; kme = grid%em33 ;

            its = grid%sp31 ; ite = grid%ep31 ; 
            jts = grid%sp32 ; jte = grid%ep32 ; 
            kts = grid%sp33 ; kte = grid%ep33 ; 

         CASE ( DATA_ORDER_XZY )
            ids = grid%sd31 ; ide = grid%ed31 ;
            kds = grid%sd32 ; kde = grid%ed32 ;
            jds = grid%sd33 ; jde = grid%ed33 ;

            ims = grid%sm31 ; ime = grid%em31 ;
            kms = grid%sm32 ; kme = grid%em32 ;
            jms = grid%sm33 ; jme = grid%em33 ;

            its = grid%sp31 ; ite = grid%ep31 ; 
            kts = grid%sp32 ; kte = grid%ep32 ; 
            jts = grid%sp33 ; jte = grid%ep33 ; 

      END SELECT

      


      flag_st000010 = 0
      flag_st010040 = 0
      flag_st040100 = 0
      flag_st100200 = 0
      flag_sm000010 = 0
      flag_sm010040 = 0
      flag_sm040100 = 0
      flag_sm100200 = 0
      flag_st010200 = 0
      flag_sm010200 = 0

      flag_soilt010 = 0
      flag_soilt040 = 0
      flag_soilt100 = 0
      flag_soilt200 = 0
      flag_soilm010 = 0
      flag_soilm040 = 0
      flag_soilm100 = 0
      flag_soilm200 = 0

      flag_sst      = 0
      flag_toposoil = 0

      

      num_st_levels_input = 0
      num_sm_levels_input = 0
      st_levels_input = -1
      sm_levels_input = -1


         CALL nl_set_iswater (grid%id, 16 )
         CALL nl_set_isice (grid%id, 24 )


      




        IF (.NOT. ALLOCATED (pmsl)              ) ALLOCATE ( pmsl(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (psfc_in)           ) ALLOCATE ( psfc_in(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (st_inputx)) ALLOCATE (st_inputx(its:ite,jts:jte,num_st_levels_alloc))
        IF (.NOT. ALLOCATED (sm_inputx)) ALLOCATE (sm_inputx(its:ite,jts:jte,num_st_levels_alloc))
        IF (.NOT. ALLOCATED (sw_inputx)) ALLOCATE (sw_inputx(its:ite,jts:jte,num_st_levels_alloc))
        IF (.NOT. ALLOCATED (soilt010_input)    ) ALLOCATE ( soilt010_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilt040_input)    ) ALLOCATE ( soilt040_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilt100_input)    ) ALLOCATE ( soilt100_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilt200_input)    ) ALLOCATE ( soilt200_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilm010_input)    ) ALLOCATE ( soilm010_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilm040_input)    ) ALLOCATE ( soilm040_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilm100_input)    ) ALLOCATE ( soilm100_input(its:ite,jts:jte) )
        IF (.NOT. ALLOCATED (soilm200_input)    ) ALLOCATE ( soilm200_input(its:ite,jts:jte) )




        




      iunit=33
      call count_recs_wrf_binary_file(iunit, trim(fileName), nrecs)
        write(0,*) 'nrecs: ', nrecs

	if (.not. allocated(varname_all)) then
      allocate (datestr_all(nrecs))
      allocate (varname_all(nrecs))
      allocate (domainend_all(3,nrecs))
      allocate (start_block(nrecs))
      allocate (end_block(nrecs))
      allocate (start_byte(nrecs))
      allocate (end_byte(nrecs))
      allocate (file_offset(nrecs))
	endif

      call inventory_wrf_binary_file(iunit, trim(filename), nrecs,  &
                      datestr_all,varname_all,domainend_all,   &
                      start_block,end_block,start_byte,end_byte,file_offset)

        do N=1,NRECS
        write(0,*) 'N,varname_all(N): ',N, varname_all(N)
        enddo

      call mpi_file_open(mpi_comm_world, trim(filename),     &
                         mpi_mode_rdonly,mpi_info_null, iunit, ierr)
      if (ierr /= 0) then
       call wrf_error_fatal3("",1573,&
"Error opening file with mpi io")
      end if

      VarName='CEN_LAT'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
      if (iret /= 0) then
        print*,VarName," not found in file"
      else

        call mpi_file_read_at(iunit,file_offset(index)+5*4,      &
                              garb,1,mpi_real4,                  &
                              mpi_status_ignore, ierr)

        if (ierr /= 0) then
          print*,"Error reading ", VarName," using MPIIO"
        else
          print*,VarName, ' from MPIIO READ= ',garb
          CALL nl_set_cen_lat ( grid%id , garb )
          write(0,*) 'cenlat= ', garb
        end if
      end if

      VarName='CEN_LON'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
        call mpi_file_read_at(iunit,file_offset(index)+5*4,      &
                              garb,1,mpi_real4,                  &
                              mpi_status_ignore, ierr)
          CALL nl_set_cen_lon ( grid%id , garb )
          CALL nl_set_stand_lon ( grid%id , garb )

      VarName='TRUELAT1'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
        call mpi_file_read_at(iunit,file_offset(index)+5*4,      &
                              garb,1,mpi_real4,                  &
                              mpi_status_ignore, ierr)
          CALL nl_set_truelat1 ( grid%id , garb )

      VarName='TRUELAT2'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
        call mpi_file_read_at(iunit,file_offset(index)+5*4,      &
                              garb,1,mpi_real4,                  &
                              mpi_status_ignore, ierr)
          CALL nl_set_truelat2 ( grid%id , garb )

      VarName='MAP_PROJ'
      call retrieve_index(index,VarName,varname_all,nrecs,iret)
        call mpi_file_read_at(iunit,file_offset(index)+5*4,     &
                              igarb,1,mpi_integer4,             &
                              mpi_status_ignore, ierr)

          CALL  nl_set_map_proj( grid%id, igarb)






    hor_size=(IDE-IDS)*(JDE-JDS)
        write(0,*) 'hor_size: ', hor_size
        write(0,*) 'IDE, JDE: ', IDE, JDE

    varName='PRES'
	if (.not. allocated(dumdata)) then
    allocate(dumdata(IDS:IDE-1,JDS:JDE-1,num_metgrid_levels))
	endif

     CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
     CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                          dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                          mpi_status_ignore, ierr)

        write(0,*) 'ierr from mpi_file_read_at for PRES: ', ierr



        write(6,*) 'post first read_from_wps_int'


       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%p_gc(I,J,K)=dumdata(I,J,K)
        if (I .eq. 1 .and. J .eq. 1) then

        endif
        if (I .eq. min(ITE,IDE-1) .and. J .eq. min(JTE,JDE-1)) then

        endif
         ENDDO
        ENDDO
       ENDDO








    varName='GHT'


    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                             mpi_status_ignore, ierr)


       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%ght_gc(I,J,K)=dumdata(I,J,K)





        if (I .eq. 1 .and. J .eq. 1) then

        endif

        if (I .eq. min(ITE,IDE-1) .and. J .eq. min(JTE,JDE-1)) then

        endif

         ENDDO
        ENDDO
       ENDDO





    varName='VEGCAT'


    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%vegcat(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='SOIL_CAT'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%input_soil_cat(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='CANWAT'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)


        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%canwat(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='SNOW'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%snow(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='SKINTEMP'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%tsk_gc(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='SOILHGT'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%toposoil(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO





    varName='SEAICE'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%xice_gc(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='ST100200'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%st100200(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='ST040100'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%st040100(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='ST010040'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)


        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%st010040(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='ST000010'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)


        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%st000010(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='SM100200'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)


        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%sm100200(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='SM040100'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)


        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%sm040100(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='SM010040'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)


        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%sm010040(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='SM000010'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)


        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%sm000010(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='PSFC'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)


        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%psfc(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='RH'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                             mpi_status_ignore, ierr)


       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%rh_gc(I,J,K)=dumdata(I,J,K)
         ENDDO
        ENDDO
       ENDDO

    varName='VV'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                             mpi_status_ignore, ierr)
       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%v_gc(I,J,K)=dumdata(I,J,K)
         ENDDO
        ENDDO
       ENDDO

    varName='UU'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                             mpi_status_ignore, ierr)

       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%u_gc(I,J,K)=dumdata(I,J,K)
         ENDDO
        ENDDO
       ENDDO

    varName='TT'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)

        if (IRET .eq. 0) then
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                             mpi_status_ignore, ierr)

       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%t_gc(I,J,K)=dumdata(I,J,K)
         ENDDO
        ENDDO
       ENDDO
        endif

    varName='RWMR'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
        if (iret .ne. 0) then





        write(0,*) 'skipping cloud reads'
        goto 979
        endif

    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                             mpi_status_ignore, ierr)
       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)

         ENDDO
        ENDDO
       ENDDO


    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                             mpi_status_ignore, ierr)
       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)

         ENDDO
        ENDDO
       ENDDO

    varName='CLWMR'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                             mpi_status_ignore, ierr)
       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)

         ENDDO
        ENDDO
       ENDDO

    varName='CICE'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                             mpi_status_ignore, ierr)
       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)

         ENDDO
        ENDDO
       ENDDO

    varName='FRIMEF'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_metgrid_levels,mpi_real4,             &
                             mpi_status_ignore, ierr)
       DO K=1,num_metgrid_levels
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)

         ENDDO
        ENDDO
       ENDDO

  979   continue




    varName='SLOPECAT'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%slopecat(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='SNOALB'
    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)


        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%snoalb(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

         num_veg_cat      = SIZE ( grid%landusef_gc , DIM=3 )
         write(0,*) 'num_veg_cat: ', num_veg_cat
         num_soil_top_cat = SIZE ( grid%soilctop_gc , DIM=3 )
         num_soil_bot_cat = SIZE ( grid%soilcbot_gc , DIM=3 )

    varName='GREENFRAC'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*12,mpi_real4,             &
                             mpi_status_ignore, ierr)

       DO K=1,12
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%greenfrac_gc(I,J,K)=dumdata(I,J,K)
         ENDDO
        ENDDO
       ENDDO

    varName='ALBEDO12M'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*12,mpi_real4,             &
                             mpi_status_ignore, ierr)

       DO K=1,12
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%albedo12m_gc(I,J,K)=dumdata(I,J,K)
         ENDDO
        ENDDO
       ENDDO

    varName='SOILCBOT'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_soil_bot_cat,mpi_real4,             &
                             mpi_status_ignore, ierr)

       DO K=1,num_soil_bot_cat
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%soilcbot_gc(I,J,K)=dumdata(I,J,K)
         ENDDO
        ENDDO
       ENDDO

    varName='SOILCAT'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%soilcat(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

        write(0,*) 'veg_cat and soil_cat sizes:::: ', num_veg_cat , num_soil_top_cat

    varName='SOILCTOP'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_soil_top_cat,mpi_real4,             &
                             mpi_status_ignore, ierr)
       DO K=1,num_soil_top_cat
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%soilctop_gc(I,J,K)=dumdata(I,J,K)
         ENDDO
        ENDDO
       ENDDO

    varName='SOILTEMP'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%tmn_gc(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='HGT_V'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%htv_gc(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='HGT_M'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%ht_gc(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='LU_INDEX'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%lu_index(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='LANDUSEF'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size*num_veg_cat,mpi_real4,             &
                             mpi_status_ignore, ierr)

       DO K=1,num_veg_cat
        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%landusef_gc(I,J,K)=dumdata(I,J,K)

        if (I .eq. 50 .and. J .eq. 50) then
        write(0,*) 'I,J,landusef_gc:: ', I,J,grid%landusef_gc(I,J,K)
        endif

         ENDDO
        ENDDO
       ENDDO

    varName='LANDMASK'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%landmask(I,J)=dumdata(I,J,1)
        if (I .eq. 50 .and. J .eq. 50) then
        write(0,*) 'I,J,landmask:: ', I,J,grid%landmask(I,J)
        endif
         ENDDO
        ENDDO







    varName='XLONG_V'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%vlon_gc(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='XLAT_V'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%vlat_gc(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='XLONG_M'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%hlon_gc(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

    varName='XLAT_M'

    CALL retrieve_index(index,VarName,varname_all,nrecs,iret)
    CALL mpi_file_read_at(iunit,file_offset(index+1),     &
                             dumdata,hor_size,mpi_real4,             &
                             mpi_status_ignore, ierr)

        DO J=JTS,min(JTE,JDE-1)
         DO I=ITS,min(ITE,IDE-1)
           grid%hlat_gc(I,J)=dumdata(I,J,1)
         ENDDO
        ENDDO

      call mpi_file_close(mpi_comm_world, ierr)

       varName='ST000010'
       flag_st000010 = 1
       num_st_levels_input = num_st_levels_input + 1
       st_levels_input(num_st_levels_input) = char2int2(varName(3:8))
       DO J=JTS,min(JTE,JDE-1)
        DO I=ITS,min(ITE,IDE-1)
           st_inputx(I,J,num_st_levels_input + 1) = grid%st000010(i,j)
        ENDDO
       ENDDO

       varName='ST010040'
       flag_st010040 = 1
       num_st_levels_input = num_st_levels_input + 1
       st_levels_input(num_st_levels_input) = char2int2(varName(3:8))
       DO J=JTS,min(JTE,JDE-1)
        DO I=ITS,min(ITE,IDE-1)
            st_inputx(I,J,num_st_levels_input + 1) = grid%st010040(i,j)
        ENDDO
       ENDDO

       varName='ST040100'
       flag_st040100 = 1
       num_st_levels_input = num_st_levels_input + 1
       st_levels_input(num_st_levels_input) = char2int2(varName(3:8))
       DO J=JTS,min(JTE,JDE-1)
        DO I=ITS,min(ITE,IDE-1)
              st_inputx(I,J,num_st_levels_input + 1) = grid%st040100(i,j)
        ENDDO
       ENDDO

       varName='ST100200'
       flag_st100200 = 1
       num_st_levels_input = num_st_levels_input + 1
       st_levels_input(num_st_levels_input) = char2int2(varName(3:8))
       DO J=JTS,min(JTE,JDE-1)
        DO I=ITS,min(ITE,IDE-1)
              st_inputx(I,J,num_st_levels_input + 1) = grid%st100200(i,j)
        ENDDO
       ENDDO

       varName='SM000010'
       flag_sm000010 = 1
       num_sm_levels_input = num_sm_levels_input + 1
       sm_levels_input(num_sm_levels_input) = char2int2(varName(3:8))
       DO J=JTS,min(JTE,JDE-1)
        DO I=ITS,min(ITE,IDE-1)
           sm_inputx(I,J,num_sm_levels_input + 1) = grid%sm000010(i,j)
        ENDDO
       ENDDO

       varName='SM010040'
       flag_sm010040 = 1
       num_sm_levels_input = num_sm_levels_input + 1
       sm_levels_input(num_sm_levels_input) = char2int2(varName(3:8))
       DO J=JTS,min(JTE,JDE-1)
        DO I=ITS,min(ITE,IDE-1)
            sm_inputx(I,J,num_sm_levels_input + 1) = grid%sm010040(i,j)
        ENDDO
       ENDDO

       varName='SM040100'
       flag_sm040100 = 1
       num_sm_levels_input = num_sm_levels_input + 1
       sm_levels_input(num_sm_levels_input) = char2int2(varName(3:8))
       DO J=JTS,min(JTE,JDE-1)
        DO I=ITS,min(ITE,IDE-1)
              sm_inputx(I,J,num_sm_levels_input + 1) = grid%sm040100(i,j)
        ENDDO
       ENDDO

       varName='SM100200'
       flag_sm100200 = 1
       num_sm_levels_input = num_sm_levels_input + 1
       sm_levels_input(num_sm_levels_input) = char2int2(varName(3:8))
       DO J=JTS,min(JTE,JDE-1)
        DO I=ITS,min(ITE,IDE-1)
              sm_inputx(I,J,num_sm_levels_input + 1) = grid%sm100200(i,j)
        ENDDO
       ENDDO

            flag_sst = 1


        sw_inputx=0.


      sw_input=0.






        write(6,*) 'to st_input...sw_input definition'


        do J=JTS,min(JDE-1,JTE)
         do K=1,num_st_levels_alloc
          do I=ITS,min(IDE-1,ITE)
             st_input(I,K,J)=st_inputx(I,J,K)
             sm_input(I,K,J)=sm_inputx(I,J,K)
             sw_input(I,K,J)=sw_inputx(I,J,K)
          enddo
         enddo
        enddo



        write(0,*) 'maxval st_input(1): ', maxval(st_input(:,1,:))
        write(0,*) 'maxval st_input(2): ', maxval(st_input(:,2,:))
        write(0,*) 'maxval st_input(3): ', maxval(st_input(:,3,:))
        write(0,*) 'maxval st_input(4): ', maxval(st_input(:,4,:))


















     end subroutine read_wps







subroutine retrieve_index(index,string,varname_all,nrecs,iret)





























  implicit none

  integer,intent(out)::iret
  integer,intent(in)::nrecs
  integer,intent(out):: index
  character(*), intent(in):: string
  character(132),intent(in)::varname_all(nrecs)

  integer i

  iret=0

  do i=1,nrecs
   if(trim(string) == trim(varname_all(i))) then
      index=i
      return
   end if
  end do

  write(6,*)' problem reading wrf nmm binary file, rec id "',trim(string),'" not found'

  iret=-1

end subroutine retrieve_index
subroutine next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)

































  implicit none

  integer(i_llong) lrecl
  integer in_unit,nreads
  integer(i_byte) buf(lrecl)
  integer(i_llong) nextbyte,locbyte,thisblock
  logical lastbuf

  integer ierr

  if(lastbuf) return

  ierr=0
  nreads=nreads+1



  thisblock = 1_i_llong + (nextbyte-1_i_llong)/lrecl

  locbyte = 1_i_llong+mod(locbyte-1_i_llong,lrecl)

  read(in_unit,rec=thisblock,iostat=ierr)buf
  lastbuf = ierr /= 0

end subroutine next_buf

subroutine inventory_wrf_binary_file(in_unit,wrf_ges_filename,nrecs, &
                                     datestr_all,varname_all,domainend_all, &
                                     start_block,end_block,start_byte,end_byte,file_offset)









































  use module_internal_header_util
  implicit none

  integer,intent(in)::in_unit,nrecs
  character(*),intent(in)::wrf_ges_filename
  character(132),intent(out)::datestr_all(nrecs),varname_all(nrecs)
  integer,intent(out)::domainend_all(3,nrecs)
  integer,intent(out)::start_block(nrecs),end_block(nrecs)
  integer,intent(out)::start_byte(nrecs),end_byte(nrecs)
  integer(i_llong),intent(out)::file_offset(nrecs)

  integer irecs
  integer(i_llong) nextbyte,locbyte,thisblock
  integer(i_byte) lenrec4(4)
  integer(i_long) lenrec,lensave
  equivalence (lenrec4(1),lenrec)
  integer(i_byte) missing4(4)
  integer(i_long) missing
  equivalence (missing,missing4(1))
  integer(i_llong),parameter:: lrecl=2**20
  integer(i_byte) buf(lrecl)
  integer i,loc_count,nreads
  logical lastbuf
  integer(i_byte) hdrbuf4(2048)
  integer(i_long) hdrbuf(512)
  equivalence(hdrbuf(1),hdrbuf4(1))
  integer,parameter:: int_field       =       530
  integer,parameter:: int_dom_ti_char =       220
  integer,parameter:: int_dom_ti_real =       140
  integer,parameter:: int_dom_ti_integer =       180
  integer hdrbufsize
  integer inttypesize
  integer datahandle,count
  character(128) element,dumstr,strdata
  integer loccode
  character(132) blanks
  integer typesize
  integer fieldtype,comm,iocomm
  integer domaindesc
  character(132) memoryorder,stagger,dimnames(3)
  integer domainstart(3),domainend(3)
  integer memorystart(3),memoryend(3)
  integer patchstart(3),patchend(3)
  character(132) datestr,varname
  real(r_single) dummy_field(1)


  integer itypesize
  integer idata(1)
  real rdata(16)

  call wrf_sizeof_integer(itypesize)
  inttypesize=itypesize

  blanks=trim(' ')

  open(in_unit,file=trim(wrf_ges_filename),access='direct',recl=lrecl)
  irecs=0
  missing=-9999
  nextbyte=0_i_llong
  locbyte=lrecl
  nreads=0
  lastbuf=.false.
  do



    do i=1,4
     nextbyte=nextbyte+1_i_llong
     locbyte=locbyte+1_i_llong
     if(locbyte > lrecl .and. lastbuf) go to 900
     if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
     lenrec4(i)=buf(locbyte)
    end do
    if(lenrec <= 0 .and. lastbuf) go to 900
    if(lenrec <= 0 .and. .not. lastbuf) go to 885
    nextbyte=nextbyte+1_i_llong
    locbyte=locbyte+1_i_llong
    if(locbyte > lrecl .and. lastbuf) go to 900
    if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)

    irecs=irecs+1
    start_block(irecs)=thisblock
    start_byte(irecs)=locbyte
    file_offset(irecs)=nextbyte-1_i_llong
    hdrbuf4(1)=buf(locbyte)
    hdrbuf4(2:4)=missing4(2:4)
    hdrbuf4(5:8)=missing4(1:4)
    datestr_all(irecs)=blanks
    varname_all(irecs)=blanks
    domainend_all(1:3,irecs)=0

    loc_count=1
    do i=2,8
       if(loc_count.ge.lenrec) exit
       loc_count=loc_count+1
       nextbyte=nextbyte+1_i_llong
       locbyte=locbyte+1_i_llong
       if(locbyte > lrecl .and. lastbuf) go to 900
       if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
       hdrbuf4(i)=buf(locbyte)
    end do

         if(lenrec==2048) write(6,*)' irecs,hdrbuf(2),int_dom_ti_char,int_field=', &
                                      irecs,hdrbuf(2),int_dom_ti_char,int_field
    if(lenrec==2048.and.(hdrbuf(2) == int_dom_ti_char .or. hdrbuf(2) == int_field &
    .or. hdrbuf(2) == int_dom_ti_real .or. hdrbuf(2) == int_dom_ti_integer)) then


       do i=9,lenrec
          loc_count=loc_count+1
          nextbyte=nextbyte+1_i_llong
          locbyte=locbyte+1_i_llong
          if(locbyte > lrecl .and. lastbuf) go to 900
          if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
          hdrbuf4(i)=buf(locbyte)
       end do

       if(hdrbuf(2) == int_dom_ti_char) then

          call int_get_ti_header_char(hdrbuf,hdrbufsize,inttypesize, &
                   datahandle,element,dumstr,strdata,loccode)
          varname_all(irecs)=trim(element)
          datestr_all(irecs)=trim(strdata)
              write(6,*)' irecs,varname,datestr = ',irecs,trim(varname_all(irecs)),trim(datestr_all(irecs))

       else if(hdrbuf(2) == int_dom_ti_real) then

          call int_get_ti_header_real(hdrbuf,hdrbufsize,inttypesize,typesize, &
                   datahandle,element,rdata,count,loccode)
          varname_all(irecs)=trim(element)

	write(6,*) 'count: ', count
              write(6,*)' irecs,varname,datestr = ',irecs,trim(varname_all(irecs)),rdata(1:count)
         
       else if(hdrbuf(2) == int_dom_ti_integer) then

          call int_get_ti_header_integer(hdrbuf,hdrbufsize,inttypesize,typesize, &
                   datahandle,element,idata,count,loccode)
          varname_all(irecs)=trim(element)

              write(6,*)' irecs,varname,datestr = ',irecs,trim(varname_all(irecs)),idata(1:count)

       else

          call int_get_write_field_header(hdrbuf,hdrbufsize,inttypesize,typesize, &
             datahandle,datestr,varname,dummy_field,fieldtype,comm,iocomm, &
             domaindesc,memoryorder,stagger,dimnames, &
             domainstart,domainend,memorystart,memoryend,patchstart,patchend)
          varname_all(irecs)=trim(varname)
          datestr_all(irecs)=trim(datestr)
          domainend_all(1:3,irecs)=domainend(1:3)
              write(6,*)' irecs,datestr,domend,varname = ', &
                  irecs,trim(datestr_all(irecs)),domainend_all(1:3,irecs),trim(varname_all(irecs))

       end if
    end if

    nextbyte=nextbyte-loc_count+lenrec
    locbyte=locbyte-loc_count+lenrec
    if(locbyte > lrecl .and. lastbuf) go to 900
    if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
    end_block(irecs)=thisblock
    end_byte(irecs)=locbyte
    lensave=lenrec
    do i=1,4
     nextbyte=nextbyte+1_i_llong
     locbyte=locbyte+1_i_llong
     if(locbyte > lrecl .and. lastbuf) go to 900
     if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
     lenrec4(i)=buf(locbyte)
    end do
    if(lenrec /= lensave) go to 890

  end do

880  continue
     write(6,*)' reached impossible place in inventory_wrf_binary_file'
     close(in_unit)
     return

885  continue
     write(6,*)' problem in inventory_wrf_binary_file, lenrec has bad value before end of file'
     write(6,*)'     lenrec =',lenrec
     close(in_unit)
     return

890  continue
     write(6,*)' problem in inventory_wrf_binary_file, beginning and ending rec len words unequal'
     write(6,*)'     begining reclen =',lensave
     write(6,*)'       ending reclen =',lenrec
     write(6,*)'               irecs =',irecs
     write(6,*)'               nrecs =',nrecs
     close(in_unit)
     return

900  continue
     write(6,*)' normal end of file reached in inventory_wrf_binary_file'
     write(6,*)'        nblocks=',thisblock
     write(6,*)'          irecs,nrecs=',irecs,nrecs
     write(6,*)'         nreads=',nreads
     close(in_unit)

end subroutine inventory_wrf_binary_file

SUBROUTINE wrf_sizeof_integer( retval )
  IMPLICIT NONE
  INTEGER retval

  retval = 4
  RETURN
END SUBROUTINE wrf_sizeof_integer

SUBROUTINE wrf_sizeof_real( retval )
  IMPLICIT NONE
  INTEGER retval

  retval = 4
  RETURN
END SUBROUTINE wrf_sizeof_real
subroutine count_recs_wrf_binary_file(in_unit,wrf_ges_filename,nrecs)








































  implicit none

  integer,intent(in)::in_unit
  character(*),intent(in)::wrf_ges_filename
  integer,intent(out)::nrecs

  integer(i_llong) nextbyte,locbyte,thisblock
  integer(i_byte) lenrec4(4)
  integer(i_long) lenrec,lensave
  equivalence (lenrec4(1),lenrec)
  integer(i_byte) missing4(4)
  integer(i_long) missing
  equivalence (missing,missing4(1))
  integer(i_llong),parameter:: lrecl=2**20
  integer(i_byte) buf(lrecl)
  integer i,loc_count,nreads
  logical lastbuf

  open(in_unit,file=trim(wrf_ges_filename),access='direct',recl=lrecl)
  nrecs=0
  missing=-9999
  nextbyte=0_i_llong
  locbyte=lrecl
  nreads=0
  lastbuf=.false.
  do



    do i=1,4
     nextbyte=nextbyte+1_i_llong
     locbyte=locbyte+1_i_llong
     if(locbyte > lrecl .and. lastbuf) go to 900
     if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
     lenrec4(i)=buf(locbyte)
    end do
    if(lenrec <= 0 .and. lastbuf) go to 900
    if(lenrec <= 0 .and. .not.lastbuf) go to 885
    nextbyte=nextbyte+1_i_llong
    locbyte=locbyte+1_i_llong
    if(locbyte > lrecl .and. lastbuf) go to 900
    if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)

    nrecs=nrecs+1

    loc_count=1
    do i=2,4
       if(loc_count.ge.lenrec) exit
       loc_count=loc_count+1
       nextbyte=nextbyte+1_i_llong
       locbyte=locbyte+1_i_llong
       if(locbyte > lrecl .and. lastbuf) go to 900
       if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
    end do
    do i=1,4
       if(loc_count.ge.lenrec) exit
       loc_count=loc_count+1
       nextbyte=nextbyte+1_i_llong
       locbyte=locbyte+1_i_llong
       if(locbyte > lrecl .and. lastbuf) go to 900
       if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
    end do
    nextbyte=nextbyte-loc_count+lenrec
    locbyte=locbyte-loc_count+lenrec
    if(locbyte > lrecl .and. lastbuf) go to 900
    if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
    lensave=lenrec
    do i=1,4
     nextbyte=nextbyte+1_i_llong
     locbyte=locbyte+1_i_llong
     if(locbyte > lrecl .and. lastbuf) go to 900
     if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
     lenrec4(i)=buf(locbyte)
    end do
    if(lenrec /= lensave) go to 890

  end do

880  continue
     write(6,*)' reached impossible place in count_recs_wrf_binary_file'
     close(in_unit)
     return

885  continue
     write(6,*)' problem in count_recs_wrf_binary_file, lenrec has bad value before end of file'
     write(6,*)'     lenrec =',lenrec
     close(in_unit)
     return

890  continue
     write(6,*)' problem in count_recs_wrf_binary_file, beginning and ending rec len words unequal'
     write(6,*)'     begining reclen =',lensave
     write(6,*)'       ending reclen =',lenrec
     close(in_unit)
     return

900  continue
     write(6,*)' normal end of file reached in count_recs_wrf_binary_file'
     write(6,*)'        nblocks=',thisblock
     write(6,*)'          nrecs=',nrecs
     write(6,*)'         nreads=',nreads
     close(in_unit)

end subroutine count_recs_wrf_binary_file

subroutine retrieve_field(in_unit,wrfges,out,start_block,end_block,start_byte,end_byte)





























 
  implicit none

  integer(i_kind),intent(in)::in_unit
  character(50),intent(in)::wrfges
  integer(i_kind),intent(in)::start_block,end_block,start_byte,end_byte
  integer(i_byte),intent(out)::out(*)

  integer(i_kind),parameter:: lrecl=2**20
  integer(i_byte) buf(lrecl)
  integer(i_kind) i,ii,k,ibegin,iend,ierr

  open(in_unit,file=trim(wrfges),access='direct',recl=lrecl)

     write(6,*)' in retrieve_field, start_block,end_block=',start_block,end_block
     write(6,*)' in retrieve_field, start_byte,end_byte=',start_byte,end_byte
  ii=0
  do k=start_block,end_block
     read(in_unit,rec=k,iostat=ierr)buf
     ibegin=1 ; iend=lrecl
     if(k == start_block) ibegin=start_byte
     if(k == end_block) iend=end_byte
     do i=ibegin,iend
        ii=ii+1
        out(ii)=buf(i)
     end do
  end do
  close(in_unit)

end subroutine retrieve_field


END MODULE module_si_io_nmm
