!STARTOFREGISTRYGENERATEDINCLUDE 'inc/em_dummy_decl.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
real                                     :: cfn
real                                     :: cfn1
real                                     :: epsts
integer                                  :: step_number
real                                     :: u_frame
real                                     :: v_frame
real                                     :: p_top
integer                                  :: stepcu
integer                                  :: stepra
integer                                  :: stepbl
logical                                  :: warm_rain
real                                     :: rdx
real                                     :: rdy
real                                     :: dts
real                                     :: dtseps
real                                     :: resm
real                                     :: zetatop
real                                     :: cf1
real                                     :: cf2
real                                     :: cf3
real                                     :: dtbc
integer                                  :: number_at_same_level
integer                                  :: itimestep
integer                                  :: oid
integer                                  :: auxhist1_oid
integer                                  :: auxhist2_oid
integer                                  :: auxhist3_oid
integer                                  :: auxhist4_oid
integer                                  :: auxhist5_oid
integer                                  :: auxinput1_oid
integer                                  :: auxinput2_oid
integer                                  :: auxinput3_oid
integer                                  :: auxinput4_oid
integer                                  :: auxinput5_oid
integer                                  :: nframes
integer                                  :: lbc_fid
logical                                  :: tiled
logical                                  :: patched
logical                                  :: input_from_file
logical                                  :: write_metadata
integer                                  :: time_step
integer                                  :: time_step_fract_num
integer                                  :: time_step_fract_den
logical                                  :: restart
integer                                  :: max_dom
integer                                  :: dyn_opt
integer                                  :: rk_ord
integer                                  :: diff_opt
integer                                  :: km_opt
integer                                  :: damp_opt
integer                                  :: isfflx
integer                                  :: ifsnow
integer                                  :: icloud
integer                                  :: num_soil_layers
integer                                  :: num_land_cat
integer                                  :: num_soil_cat
integer                                  :: spec_bdy_width
integer                                  :: spec_zone
integer                                  :: relax_zone
integer                                  :: ensdim
integer                                  :: tile_sz_x
integer                                  :: tile_sz_y
integer                                  :: numtiles
integer                                  :: debug_level
integer                                  :: irand
integer                                  :: run_days
integer                                  :: run_hours
integer                                  :: run_minutes
integer                                  :: run_seconds
integer                                  :: start_year
integer                                  :: start_month
integer                                  :: start_day
integer                                  :: start_hour
integer                                  :: start_minute
integer                                  :: start_second
integer                                  :: end_year
integer                                  :: end_month
integer                                  :: end_day
integer                                  :: end_hour
integer                                  :: end_minute
integer                                  :: end_second
integer                                  :: grid_id
integer                                  :: level
integer                                  :: s_we
integer                                  :: e_we
integer                                  :: s_sn
integer                                  :: e_sn
integer                                  :: s_vert
integer                                  :: e_vert
integer                                  :: history_interval
integer                                  :: auxhist1_interval
integer                                  :: auxhist2_interval
integer                                  :: auxhist3_interval
integer                                  :: auxhist4_interval
integer                                  :: auxhist5_interval
integer                                  :: auxinput1_interval
integer                                  :: auxinput2_interval
integer                                  :: auxinput3_interval
integer                                  :: auxinput4_interval
integer                                  :: auxinput5_interval
integer                                  :: restart_interval
integer                                  :: frames_per_outfile
integer                                  :: time_step_sound
integer                                  :: parent_id
integer                                  :: i_parent_start
integer                                  :: j_parent_start
integer                                  :: shw
integer                                  :: parent_grid_ratio
integer                                  :: parent_time_step_ratio
integer                                  :: moad_grid_ratio
integer                                  :: moad_time_step_ratio
logical                                  :: non_hydrostatic
real                                     :: dx
real                                     :: dy
real                                     :: dt
real                                     :: ztop
real                                     :: zdamp
real                                     :: dampcoef
real                                     :: smdiv
real                                     :: emdiv
real                                     :: epssm
real                                     :: khdif
real                                     :: kvdif
real                                     :: mix_cr_len
real                                     :: tke_upper_bound
real                                     :: kh_tke_upper_bound
real                                     :: kv_tke_upper_bound
real                                     :: radt
real                                     :: bldt
real                                     :: cudt
real                                     :: gsmdt
integer                                  :: julyr
integer                                  :: julday
real                                     :: gmt
logical                                  :: periodic_x
logical                                  :: symmetric_xs
logical                                  :: symmetric_xe
logical                                  :: open_xs
logical                                  :: open_xe
logical                                  :: periodic_y
logical                                  :: symmetric_ys
logical                                  :: symmetric_ye
logical                                  :: open_ys
logical                                  :: open_ye
logical                                  :: nested
logical                                  :: specified
logical                                  :: top_radiation
integer                                  :: idtad
integer                                  :: nsoil
integer                                  :: nphs
integer                                  :: ncnvc
integer                                  :: nrads
integer                                  :: nradl
integer                                  :: sigma
integer                                  :: chem_opt
integer                                  :: mp_physics
integer                                  :: ra_lw_physics
integer                                  :: ra_sw_physics
integer                                  :: bl_sfclay_physics
integer                                  :: bl_surface_physics
integer                                  :: bl_pbl_physics
integer                                  :: cu_physics
integer                                  :: h_mom_adv_order
integer                                  :: v_mom_adv_order
integer                                  :: h_sca_adv_order
integer                                  :: v_sca_adv_order
integer                                  :: io_form_input
integer                                  :: io_form_auxinput1
integer                                  :: io_form_auxinput2
integer                                  :: io_form_auxinput3
integer                                  :: io_form_auxinput4
integer                                  :: io_form_auxinput5
integer                                  :: io_form_history
integer                                  :: io_form_auxhist1
integer                                  :: io_form_auxhist2
integer                                  :: io_form_auxhist3
integer                                  :: io_form_auxhist4
integer                                  :: io_form_auxhist5
integer                                  :: io_form_restart
integer                                  :: io_form_boundary
integer                                  :: interval_seconds
integer                                  :: real_data_init_type
real                                     :: cen_lat
real                                     :: cen_lon
real                                     :: truelat1
real                                     :: truelat2
real                                     :: bdyfrq
integer                                  :: iswater
integer                                  :: isice
integer                                  :: map_proj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: u_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: u_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ru
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: v_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: v_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: w_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: w_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ww
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ph_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ph_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: phb
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: ph0
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: php
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t_init
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: tp_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: tp_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mu_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mu_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mub
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mu0
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mudf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: tke_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: tke_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: p
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: al
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: alt
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: alb
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: zx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: zy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rdz
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rdzw
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: q2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sm000010
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sm010040
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sm040100
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sm100200
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sm010200
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilm000
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilm005
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilm020
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilm040
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilm160
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilm300
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sw000010
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sw010040
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sw040100
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sw100200
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sw010200
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilw000
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilw005
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilw020
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilw040
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilw160
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilw300
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: st000010
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: st010040
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: st040100
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: st100200
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: st010200
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilt000
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilt005
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilt020
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilt040
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilt160
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: soilt300
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: landmask
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: topostdv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: toposlpx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: toposlpy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: greenmax
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: greenmin
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: albedomx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: slopecat
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: toposoil
real      ,DIMENSION(grid%sm31:grid%em31,grid%num_land_cat,grid%sm33:grid%em33)           :: landusef
real      ,DIMENSION(grid%sm31:grid%em31,grid%num_soil_cat,grid%sm33:grid%em33)           :: soilctop
real      ,DIMENSION(grid%sm31:grid%em31,grid%num_soil_cat,grid%sm33:grid%em33)           :: soilcbot
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_moist)           :: moist_2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_chem)           :: chem_1
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33,num_chem)           :: chem_2
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: u_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: u_bt
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: v_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: v_bt
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: ph_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: ph_bt
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: t_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: t_bt
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: mu_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: mu_bt
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: rqv_b
real      ,DIMENSION(max(grid%ed31,grid%ed33),grid%sd32:grid%ed32,grid%spec_bdy_width,4)           :: rqv_bt
real      ,DIMENSION(grid%sm31:grid%em31,grid%num_soil_layers,grid%sm33:grid%em33)           :: smois
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: pb
real      ,DIMENSION(grid%sm32:grid%em32)           :: fnm
real      ,DIMENSION(grid%sm32:grid%em32)           :: fnp
real      ,DIMENSION(grid%sm32:grid%em32)           :: rdnw
real      ,DIMENSION(grid%sm32:grid%em32)           :: rdn
real      ,DIMENSION(grid%sm32:grid%em32)           :: dnw
real      ,DIMENSION(grid%sm32:grid%em32)           :: dn
real      ,DIMENSION(grid%sm32:grid%em32)           :: znu
real      ,DIMENSION(grid%sm32:grid%em32)           :: znw
real      ,DIMENSION(grid%sm32:grid%em32)           :: t_base
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: z
real      ,DIMENSION(grid%spec_bdy_width)           :: fcx
real      ,DIMENSION(grid%spec_bdy_width)           :: gcx
real      ,DIMENSION(grid%sm31:grid%em31,grid%num_soil_layers,grid%sm33:grid%em33)           :: tslb
real      ,DIMENSION(grid%num_soil_layers)           :: zs
real      ,DIMENSION(grid%num_soil_layers)           :: dzs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: h_diabatic
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: msft
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: msfu
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: msfv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: f
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: e
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sina
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: cosa
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: ht
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: tsk
real      ,DIMENSION(grid%sm32:grid%em32)           :: u_base
real      ,DIMENSION(grid%sm32:grid%em32)           :: v_base
real      ,DIMENSION(grid%sm32:grid%em32)           :: qv_base
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rthcuten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rqvcuten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rqrcuten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rqccuten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rqscuten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rqicuten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: w0avg
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: rainc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: rainnc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: raincv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: rainncv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: rainbl
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: nca
integer   ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: lowlyr
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rthraten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rthratenlw
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rthratensw
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: cldfra
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: gsw
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: glw
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: xlat
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: xlong
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: albedo
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: emiss
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: lu_index
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: cldefi
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rublten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rvblten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rthblten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rqvblten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rqcblten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rqiblten
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: tmn
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: xland
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: ust
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: pblh
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: capg
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: thc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: hfx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: qfx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: snowc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mavail
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: tkesfcf
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: taucldi
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: taucldc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: defor11
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: defor22
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: defor12
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: defor33
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: defor13
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: defor23
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: xkmv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: xkmh
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: xkmhd
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: xkhv
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: xkhh
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: div
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: bn2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: th2
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: u10
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: v10
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: xice
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: smstav
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: smstot
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sfcrunoff
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: udrunoff
integer   ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: ivgtyp
integer   ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: isltyp
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: vegfra
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sfcevp
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: grdflx
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: albbck
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sfcexc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: acsnow
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: acsnom
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: snow
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: canwat
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: sst
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: weasd
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mol
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: znt
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: tke_myj
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: thz0
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: qz0
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: uz0
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: vz0
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: uz0h
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: vz0h
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: dudt
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: dvdt
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: qsfc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: akhs
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: akms
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: htop
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: hbot
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: cuppt
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: t0eta
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: q0eta
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: p0eta
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: f_ice_phy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: f_rain_phy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: f_rimef_phy
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: mass_flux
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: apr_gr
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: apr_w
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: apr_mc
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: apr_st
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: apr_as
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: apr_capma
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: apr_capme
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: apr_capmi
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%ensdim)           :: xf_ens
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33,grid%ensdim)           :: pr_ens
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rthften
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33)           :: rqvften
real      ,DIMENSION(grid%sm31:grid%em31,grid%sm33:grid%em33)           :: snowh
real      ,DIMENSION(grid%sm31:grid%em31,grid%num_soil_layers,grid%sm33:grid%em33)           :: smfr3d
#undef COPY_IN
#undef COPY_OUT
!ENDOFREGISTRYGENERATEDINCLUDE