








MODULE module_first_rk_step_part1

CONTAINS

  SUBROUTINE first_rk_step_part1 (   grid , config_flags              &
                             , moist , moist_tend               &
                             , chem  , chem_tend                &
                             , tracer, tracer_tend              &
                             , scalar , scalar_tend             &
                             , fdda3d, fdda2d                   &
                             , aerod                            &
                             , ru_tendf, rv_tendf               &
                             , rw_tendf, t_tendf                &
                             , ph_tendf, mu_tendf               &
                             , tke_tend                         &
                             , adapt_step_flag , curr_secs      &
                             , psim , psih , wspd , gz1oz0 , br , chklowq &
                             , cu_act_flag , hol , th_phy       &
                             , pi_phy , p_phy , t_phy           &
                             , dz8w , p8w , t8w                 &
                             , ids, ide, jds, jde, kds, kde     &
                             , ims, ime, jms, jme, kms, kme     &
                             , ips, ipe, jps, jpe, kps, kpe     &
                             , imsx,imex,jmsx,jmex,kmsx,kmex    &
                             , ipsx,ipex,jpsx,jpex,kpsx,kpex    &
                             , imsy,imey,jmsy,jmey,kmsy,kmey    &
                             , ipsy,ipey,jpsy,jpey,kpsy,kpey    &
                             , k_start , k_end                  &
                             , f_flux                           &
                            )
    USE module_state_description
    USE module_model_constants
    USE module_domain, ONLY : domain, domain_clock_get, get_ijk_from_subgrid
    USE module_configure, ONLY : grid_config_rec_type, model_config_rec
    USE module_radiation_driver, ONLY : pre_radiation_driver, radiation_driver
    USE module_surface_driver, ONLY : surface_driver
    USE module_cumulus_driver, ONLY : cumulus_driver
    USE module_shallowcu_driver, ONLY : shallowcu_driver
    USE module_pbl_driver, ONLY : pbl_driver
    USE module_fr_fire_driver_wrf, ONLY : fire_driver_em_step
    USE module_fddagd_driver, ONLY : fddagd_driver
    USE module_em, ONLY : init_zero_tendency
    USE module_force_scm
    USE module_convtrans_prep
    USE module_big_step_utilities_em, ONLY : phy_prep
use module_scalar_tables

    USE module_dm, ONLY : local_communicator, mytask, ntasks, ntasks_x, ntasks_y, local_communicator_periodic, wrf_dm_maxval
    USE module_comm_dm, ONLY : halo_em_phys_a_sub,halo_em_fdda_sfc_sub,halo_pwp_sub,halo_em_chem_e_3_sub, &
    halo_em_chem_e_5_sub, halo_em_hydro_noahmp_sub

    USE module_utility
    IMPLICIT NONE

    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    TYPE(WRFU_Time)                :: currentTime

    INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde,     &
                           ims, ime, jms, jme, kms, kme,     &
                           ips, ipe, jps, jpe, kps, kpe,     &
                           imsx,imex,jmsx,jmex,kmsx,kmex,    &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,    &
                           imsy,imey,jmsy,jmey,kmsy,kmey,    &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey


    LOGICAL ,INTENT(IN)                        :: adapt_step_flag
    REAL, INTENT(IN)                           :: curr_secs

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: moist
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: moist_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: chem
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: chem_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_fdda3d),INTENT(INOUT)  :: fdda3d
    REAL    ,DIMENSION(ims:ime,1:1,jms:jme,num_fdda2d),INTENT(INOUT)      :: fdda2d
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_aerod),INTENT(INOUT)   :: aerod
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: psim
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: psih
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wspd
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: gz1oz0
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: br
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: chklowq
    LOGICAL ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: cu_act_flag
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: hol

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: th_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: pi_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: p_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: dz8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: p8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t8w

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ru_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: rv_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: rw_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ph_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: tke_tend

    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: mu_tendf

    INTEGER, INTENT(IN)                           ::  k_start, k_end
    LOGICAL, INTENT(IN), OPTIONAL                 ::  f_flux


    real :: HYDRO_dt
    REAL, DIMENSION( ims:ime, jms:jme ) :: exch_temf  

    REAL, DIMENSION( ims:ime, jms:jme ) :: ht_loc, mixht
    INTEGER                             :: ij
    INTEGER  num_roof_layers
    INTEGER  num_wall_layers
    INTEGER  num_road_layers
    INTEGER  iswater
    LOGICAL  :: l_flux
    INTEGER  :: isurban
    INTEGER  rk_step
    INTEGER                         :: yr, month, day, hr, minute, sec, rc
    CHARACTER*80                    :: mesg

   INTEGER                         :: sids , side , sjds , sjde , skds , skde , &
                                      sims , sime , sjms , sjme , skms , skme , &
                                      sips , sipe , sjps , sjpe , skps , skpe

   CHARACTER (LEN=256) :: mminlu
   CHARACTER (LEN=1000) :: message





  CALL get_ijk_from_subgrid (  grid ,                   &
                            sids, side, sjds, sjde, skds, skde,    &
                            sims, sime, sjms, sjme, skms, skme,    &
                            sips, sipe, sjps, sjpe, skps, skpe    )

 
 

   l_flux=.FALSE.
   if (present(f_flux)) l_flux=f_flux

    rk_step = 1


       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles

         CALL wrf_debug ( 200 , ' call init_zero_tendency' )
         CALL init_zero_tendency ( ru_tendf, rv_tendf, rw_tendf,     &
                                   ph_tendf, t_tendf, tke_tend,      &
                                   mu_tendf,                         &
                                   moist_tend,chem_tend,scalar_tend, &
                                   tracer_tend,num_tracer,           &
                                   num_moist,num_chem,num_scalar,    &
                                   rk_step,                          &
                                   ids, ide, jds, jde, kds, kde,     &
                                   ims, ime, jms, jme, kms, kme,     &
                                   grid%i_start(ij), grid%i_end(ij), &
                                   grid%j_start(ij), grid%j_end(ij), &
                                   k_start, k_end                   )

       END DO
       !$OMP END PARALLEL DO








CALL HALO_EM_PHYS_A_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )


      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles

        CALL wrf_debug ( 200 , ' call phy_prep' )
        CALL phy_prep ( config_flags,                                    &
                        grid%mut, grid%muu, grid%muv, grid%u_2,          &
                        grid%v_2, grid%p, grid%pb, grid%alt,             &
                        grid%ph_2, grid%phb, grid%t_2, grid%tsk, moist, num_moist,   &
                        grid%rho,th_phy, p_phy, pi_phy, grid%u_phy, grid%v_phy,      &
                        p8w, t_phy, t8w, grid%z, grid%z_at_w, dz8w,      &
                        grid%p_hyd, grid%p_hyd_w, grid%dnw,              &
                        grid%fnm, grid%fnp, grid%znw, grid%p_top,        &
                        grid%rthraten,                                   &
                        grid%rthblten, grid%rublten, grid%rvblten,       &
                        grid%rqvblten, grid%rqcblten, grid%rqiblten,     &
                        grid%rucuten,  grid%rvcuten,  grid%rthcuten,     &
                        grid%rqvcuten, grid%rqccuten, grid%rqrcuten,     &
                        grid%rqicuten, grid%rqscuten,                    &
                        grid%rushten,  grid%rvshten,  grid%rthshten,     &
                        grid%rqvshten, grid%rqcshten, grid%rqrshten,     &
                        grid%rqishten, grid%rqsshten, grid%rqgshten,     &
                        grid%rthften,  grid%rqvften,                     &
                        grid%RUNDGDTEN, grid%RVNDGDTEN, grid%RTHNDGDTEN, &
                        grid%RPHNDGDTEN,grid%RQVNDGDTEN, grid%RMUNDGDTEN,&

                        grid%landmask,grid%xland,                 &

                        ids, ide, jds, jde, kds, kde,                    &
                        ims, ime, jms, jme, kms, kme,                    &
                        grid%i_start(ij), grid%i_end(ij),                &
                        grid%j_start(ij), grid%j_end(ij),                &
                        k_start, k_end                                   )
      ENDDO
      !$OMP END PARALLEL DO




     CALL domain_clock_get( grid, current_time=currentTime, &
                            current_timestr=mesg )
     CALL WRFU_TimeGet( currentTime, YY=yr, dayOfYear=day, H=hr, M=minute, S=sec, rc=rc)
         IF( rc/= WRFU_SUCCESS)THEN
         CALL wrf_error_fatal3("<stdin>",239,&
'WRFU_TimeGet failed')
         ENDIF


      CALL pre_radiation_driver ( grid, config_flags                        &
     &        ,itimestep=grid%itimestep, ra_call_offset=grid%ra_call_offset    &
     &        ,XLAT=grid%xlat, XLONG=grid%xlong, GMT=grid%gmt                  &
     &        ,julian=grid%julian, xtime=grid%xtime, RADT=grid%radt            &
     &        ,STEPRA=grid%stepra                                              &
     &        ,ht=grid%ht,dx=grid%dx,dy=grid%dy,sina=grid%sina,cosa=grid%cosa  &
     &        ,shadowmask=grid%shadowmask,slope_rad=config_flags%slope_rad     &
     &        ,topo_shading=config_flags%topo_shading                          &
     &        ,shadlen=config_flags%shadlen,ht_shad=grid%ht_shad,ht_loc=ht_loc &
     &        ,ht_shad_bxs=grid%ht_shad_bxs, ht_shad_bxe=grid%ht_shad_bxe      &
     &        ,ht_shad_bys=grid%ht_shad_bys, ht_shad_bye=grid%ht_shad_bye      &
     &        ,nested=config_flags%nested, min_ptchsz=grid%min_ptchsz          &
     &        ,spec_bdy_width=config_flags%spec_bdy_width                      &
            
     &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
     &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
     &        ,IPS=ips,IPE=ipe, JPS=jps,JPE=jpe, KPS=kps,KPE=kpe          &
     &        ,i_start=grid%i_start,i_end=min(grid%i_end, ide-1)          &
     &        ,j_start=grid%j_start,j_end=min(grid%j_end, jde-1)          &
     &        ,kts=k_start, kte=min(k_end,kde-1)                          &
     &        ,num_tiles=grid%num_tiles                                   )

      CALL wrf_debug ( 200 , ' call radiation_driver' )


CALL radiation_driver(p_top=grid%p_top,ACFRCV=grid%acfrcv,ACFRST=grid%acfrst,ALBEDO=grid%albedo,CFRACH=grid%cfrach,CFRACL=grid%cfr&
&acl,CFRACM=grid%cfracm,CUPPT=grid%cuppt,CZMEAN=grid%czmean,DT=grid%dt,DZ8W=dz8w,EMISS=grid%emiss,GLW=grid%glw,GMT=grid%gmt,GSW=gri&
&d%gsw,HBOT=grid%hbot,HTOP=grid%htop,HBOTR=grid%hbotr,HTOPR=grid%htopr,ICLOUD=config_flags%icloud,ITIMESTEP=grid%itimestep,JULDAY=g&
&rid%julday,JULIAN=grid%julian,JULYR=grid%julyr,LW_PHYSICS=config_flags%ra_lw_physics,NCFRCV=grid%ncfrcv,NCFRST=grid%ncfrst,NPHS=1,&
&o3input=config_flags%o3input,O3rad=grid%o3rad,aer_opt=config_flags%aer_opt,aerod=aerod(:,:,:,P_ocarbon:P_upperaer),swint_opt=confi&
&g_flags%swint_opt,P8W=grid%p_hyd_w,P=grid%p_hyd,PI=pi_phy,RADT=grid%radt,RA_CALL_OFFSET=grid%ra_call_offset,RHO=grid%rho,RLWTOA=gr&
&id%rlwtoa,RSWTOA=grid%rswtoa,RTHRATEN=grid%rthraten,RTHRATENLW=grid%rthratenlw,RTHRATENSW=grid%rthratensw,SNOW=grid%snow,STEPRA=gr&
&id%stepra,SWDOWN=grid%swdown,SWDOWNC=grid%swdownc,SW_PHYSICS=config_flags%ra_sw_physics,T8W=t8w,T=grid%t_phy,TAUCLDC=grid%taucldc,&
&TAUCLDI=grid%taucldi,TSK=grid%tsk,VEGFRA=grid%vegfra,WARM_RAIN=grid%warm_rain,XICE=grid%xice,XLAND=grid%xland,XLAT=grid%xlat,XLONG&
&=grid%xlong,YR=yr,ALSWVISDIR=grid%alswvisdir,ALSWVISDIF=grid%alswvisdif,ALSWNIRDIR=grid%alswnirdir,ALSWNIRDIF=grid%alswnirdif,SWVI&
&SDIR=grid%swvisdir,SWVISDIF=grid%swvisdif,SWNIRDIR=grid%swnirdir,SWNIRDIF=grid%swnirdif,SF_SURFACE_PHYSICS=config_flags%sf_surface&
&_physics,SWDDIR=grid%swddir,SWDDNI=grid%swddni,SWDDIF=grid%swddif,Gx=grid%Gx,Bx=grid%Bx,gg=grid%gg,bb=grid%bb,swdown_ref=grid%swdo&
&wn_ref,swddir_ref=grid%swddir_ref,coszen_ref=grid%coszen_ref,aer_type=config_flags%aer_type,aer_aod550_opt=config_flags%aer_aod550&
&_opt,aer_aod550_val=config_flags%aer_aod550_val,aer_angexp_opt=config_flags%aer_angexp_opt,aer_angexp_val=config_flags%aer_angexp_&
&val,aer_ssa_opt=config_flags%aer_ssa_opt,aer_ssa_val=config_flags%aer_ssa_val,aer_asy_opt=config_flags%aer_asy_opt,aer_asy_val=con&
&fig_flags%aer_asy_val,aod5502d=grid%aod5502d,angexp2d=grid%angexp2d,aerssa2d=grid%aerssa2d,aerasy2d=grid%aerasy2d,aod5503d=grid%ao&
&d5503d,DECLINX=grid%declin,SOLCONX=grid%solcon,COSZEN=grid%coszen,HRANG=grid%hrang,CEN_LAT=grid%cen_lat,Z=grid%z,ALEVSIZ=grid%alev&
&siz,no_src_types=grid%no_src_types,LEVSIZ=grid%levsiz,N_OZMIXM=num_ozmixm,N_AEROSOLC=num_aerosolc,PAERLEV=grid%paerlev,ID=grid%id,&
&CAM_ABS_DIM1=grid%cam_abs_dim1,CAM_ABS_DIM2=grid%cam_abs_dim2,CAM_ABS_FREQ_S=grid%cam_abs_freq_s,XTIME=grid%xtime,CURR_SECS=curr_s&
&ecs,ADAPT_STEP_FLAG=adapt_step_flag,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=km&
&e,i_start=grid%i_start,i_end=min(grid%i_end,ide-1),j_start=grid%j_start,j_end=min(grid%j_end,jde-1),kts=k_start,kte=min(k_end,kde-&
&1),num_tiles=grid%num_tiles,TLWDN=grid%tlwdn,TLWUP=grid%tlwup,SLWDN=grid%slwdn,SLWUP=grid%slwup,TSWDN=grid%tswdn,TSWUP=grid%tswup,&
&SSWDN=grid%sswdn,SSWUP=grid%sswup,CLDFRA=grid%cldfra,CLDFRA_MP_ALL=grid%cldfra_mp_all,LRADIUS=grid%LRADIUS,IRADIUS=grid%IRADIUS,CL&
&DFRA_DP=grid%cldfra_dp,CLDFRA_SH=grid%cldfra_sh,re_cloud=grid%re_cloud,re_ice=grid%re_ice,re_snow=grid%re_snow,has_reqc=grid%has_r&
&eqc,has_reqi=grid%has_reqi,has_reqs=grid%has_reqs,PB=grid%pb,F_ICE_PHY=grid%f_ice_phy,F_RAIN_PHY=grid%f_rain_phy,QV=moist(ims,kms,&
&jms,P_QV),F_QV=F_QV,QC=moist(ims,kms,jms,P_QC),F_QC=F_QC,QR=moist(ims,kms,jms,P_QR),F_QR=F_QR,QI=moist(ims,kms,jms,P_QI),F_QI=F_QI&
&,QS=moist(ims,kms,jms,P_QS),F_QS=F_QS,QG=moist(ims,kms,jms,P_QG),F_QG=F_QG,QNDROP=scalar(ims,kms,jms,P_QNDROP),F_QNDROP=F_QNDROP,A&
&CSWUPT=grid%acswupt,ACSWUPTC=grid%acswuptc,ACSWDNT=grid%acswdnt,ACSWDNTC=grid%acswdntc,ACSWUPB=grid%acswupb,ACSWUPBC=grid%acswupbc&
&,ACSWDNB=grid%acswdnb,ACSWDNBC=grid%acswdnbc,ACLWUPT=grid%aclwupt,ACLWUPTC=grid%aclwuptc,ACLWDNT=grid%aclwdnt,ACLWDNTC=grid%aclwdn&
&tc,ACLWUPB=grid%aclwupb,ACLWUPBC=grid%aclwupbc,ACLWDNB=grid%aclwdnb,ACLWDNBC=grid%aclwdnbc,SWUPT=grid%swupt,SWUPTC=grid%swuptc,SWD&
&NT=grid%swdnt,SWDNTC=grid%swdntc,SWUPB=grid%swupb,SWUPBC=grid%swupbc,SWDNB=grid%swdnb,SWDNBC=grid%swdnbc,LWUPT=grid%lwupt,LWUPTC=g&
&rid%lwuptc,LWDNT=grid%lwdnt,LWDNTC=grid%lwdntc,LWUPB=grid%lwupb,LWUPBC=grid%lwupbc,LWDNB=grid%lwdnb,LWDNBC=grid%lwdnbc,LWCF=grid%l&
&wcf,SWCF=grid%swcf,OLR=grid%olr,AERODM=grid%aerodm,PINA=grid%pina,AODTOT=grid%aodtot,OZMIXM=grid%ozmixm,PIN=grid%pin,M_PS_1=grid%m&
&_ps_1,M_PS_2=grid%m_ps_2,AEROSOLC_1=grid%aerosolc_1,AEROSOLC_2=grid%aerosolc_2,M_HYBI0=grid%m_hybi,ABSTOT=grid%abstot,ABSNXT=grid%&
&absnxt,EMSTOT=grid%emstot,RADTACTTIME=grid%radtacttime,ICLOUD_CU=config_flags%ICLOUD_CU,QC_CU=grid%QC_CU,QI_CU=grid%QI_CU,slope_ra&
&d=config_flags%slope_rad,topo_shading=config_flags%topo_shading,shadowmask=grid%shadowmask,ht=grid%ht,dx=grid%dx,dy=grid%dy,IS_CAM&
&MGMP_USED=grid%is_CAMMGMP_used)









  IF ( config_flags%sf_surface_physics.eq.NOAHMPSCHEME ) THEN
       IF ( config_flags%opt_run.eq.5.and.mod(grid%itimestep,grid%STEPWTD).eq.0 )  THEN






CALL HALO_EM_HYDRO_NOAHMP_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )

       ENDIF
  ENDIF




      num_roof_layers = grid%num_soil_layers 
      num_wall_layers = grid%num_soil_layers 
      num_road_layers = grid%num_soil_layers 
      CALL nl_get_iswater(grid%id, iswater)
      CALL nl_get_isurban(grid%id, isurban)
      call nl_get_mminlu(grid%id, mminlu)







CALL HALO_PWP_sub ( grid, &
  config_flags, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )


      CALL wrf_debug ( 200 , ' call surface_driver' )

      if( grid%num_nests .lt. 1 )then
          HYDRO_dt = 0
      else
          HYDRO_dt = -1
      endif

CALL surface_driver(HYDRO_dt=HYDRO_dt,sfcheadrt=grid%sfcheadrt,INFXSRT=grid%INFXSRT,soldrain=grid%soldrain,ACGRDFLX=grid%acgrdflx,&
&ACHFX=grid%achfx,ACLHF=grid%aclhf,ACSNOM=grid%acsnom,ACSNOW=grid%acsnow,AKHS=grid%akhs,AKMS=grid%akms,ALBBCK=grid%albbck,ALBEDO=gr&
&id%albedo,EMBCK=grid%embck,BR=br,CANWAT=grid%canwat,CHKLOWQ=chklowq,CT=grid%ct,DT=grid%dt,DX=grid%dx,DZ8W=dz8w,DZS=grid%dzs,FLHC=g&
&rid%flhc,FM=grid%fm,FHH=grid%fh,FLQC=grid%flqc,GLW=grid%glw,GRDFLX=grid%grdflx,GSW=grid%gsw,SWDOWN=grid%swdown,GZ1OZ0=gz1oz0,HFX=g&
&rid%hfx,HT=grid%ht,IFSNOW=config_flags%ifsnow,ISFFLX=config_flags%isfflx,FRACTIONAL_SEAICE=config_flags%fractional_seaice,SEAICE_A&
&LBEDO_OPT=config_flags%seaice_albedo_opt,SEAICE_ALBEDO_DEFAULT=config_flags%seaice_albedo_default,SEAICE_THICKNESS_OPT=config_flag&
&s%seaice_thickness_opt,SEAICE_THICKNESS_DEFAULT=config_flags%seaice_thickness_default,SEAICE_SNOWDEPTH_OPT=config_flags%seaice_sno&
&wdepth_opt,SEAICE_SNOWDEPTH_MAX=config_flags%seaice_snowdepth_max,SEAICE_SNOWDEPTH_MIN=config_flags%seaice_snowdepth_min,TICE2TSK_&
&IF2COLD=config_flags%tice2tsk_if2cold,IFNDALBSI=grid%ifndalbsi,IFNDICEDEPTH=grid%ifndicedepth,IFNDSNOWSI=grid%ifndsnowsi,ISLTYP=gr&
&id%isltyp,ITIMESTEP=grid%itimestep,JULIAN_IN=grid%julian,IVGTYP=grid%ivgtyp,LH=grid%lh,LOWLYR=grid%lowlyr,MAVAIL=grid%mavail,NUM_S&
&OIL_LAYERS=config_flags%num_soil_layers,P8W=grid%p_hyd_w,PBLH=grid%pblh,PI_PHY=pi_phy,PSFC=grid%psfc,PSHLTR=grid%pshltr,PSIH=psih,&
&BLDT=grid%bldt,CURR_SECS=curr_secs,ADAPT_STEP_FLAG=adapt_step_flag,BLDTACTTIME=grid%bldtacttime,PSIM=psim,P_PHY=grid%p_hyd,Q10=gri&
&d%q10,Q2=grid%q2,QFX=grid%qfx,QSFC=grid%qsfc,QSHLTR=grid%qshltr,QZ0=grid%qz0,RAINCV=grid%raincv,RA_LW_PHYSICS=config_flags%ra_lw_p&
&hysics,RHO=grid%rho,RMOL=grid%rmol,SFCEVP=grid%sfcevp,SFCEXC=grid%sfcexc,SFCRUNOFF=grid%sfcrunoff,SF_SFCLAY_PHYSICS=config_flags%s&
&f_sfclay_physics,SF_SURFACE_PHYSICS=config_flags%sf_surface_physics,SH2O=grid%sh2o,SHDMAX=grid%shdmax,SHDMIN=grid%shdmin,SMOIS=gri&
&d%smois,SMSTAV=grid%smstav,SMSTOT=grid%smstot,SNOALB=grid%snoalb,SNOW=grid%snow,SNOWC=grid%snowc,SNOWH=grid%snowh,SMCREL=grid%smcr&
&el,SST=grid%sst,SST_INPUT=grid%sst_input,SST_UPDATE=grid%sst_update,SSTSK=grid%sstsk,DTW=grid%dtw,SST_SKIN=grid%sst_skin,SCM_FORCE&
&_SKINTEMP=grid%scm_force_skintemp,SCM_FORCE_FLUX=grid%scm_force_flux,STEPBL=grid%stepbl,TH10=grid%th10,TH2=grid%th2,THZ0=grid%thz0&
&,TH_PHY=th_phy,TKE_PBL=grid%tke_pbl,TMN=grid%tmn,TSHLTR=grid%tshltr,TSK=grid%tsk,TYR=grid%tyr,TYRA=grid%tyra,TDLY=grid%tdly,TLAG=g&
&rid%tlag,LAGDAY=config_flags%lagday,NYEAR=grid%nyear,NDAY=grid%nday,TMN_UPDATE=grid%tmn_update,YR=yr,TSLB=grid%tslb,T_PHY=t_phy,U1&
&0=grid%u10,URATX=grid%uratx,VRATX=grid%vratx,TRATX=grid%tratx,UDRUNOFF=grid%udrunoff,UST=grid%ust,UZ0=grid%uz0,U_FRAME=grid%u_fram&
&e,U_PHY=grid%u_phy,V10=grid%v10,UOCE=grid%uoce,VOCE=grid%voce,VEGFRA=grid%vegfra,VZ0=grid%vz0,V_FRAME=grid%v_frame,V_PHY=grid%v_ph&
&y,WARM_RAIN=grid%warm_rain,WSPD=wspd,XICE=grid%xice,XLAND=grid%xland,MAX_EDOM=grid%num_ext_model_couple_dom,CPLMASK=grid%cplmask,Z&
&0=grid%z0,Z=grid%z,ZNT=grid%znt,ZS=grid%zs,ALBSI=grid%albsi,ICEDEPTH=grid%icedepth,SNOWSI=grid%snowsi,XICEM=grid%xicem,ISICE=grid%&
&landuse_isice,USTM=grid%ustm,CK=grid%ck,CKA=grid%cka,CD=grid%cd,CDA=grid%cda,ISFTCFLX=config_flags%isftcflx,IZ0TLND=config_flags%i&
&z0tlnd,SF_OCEAN_PHYSICS=config_flags%sf_ocean_physics,OML_HML0=config_flags%oml_hml0,OML_GAMMA=config_flags%oml_gamma,TML=grid%tml&
&,T0ML=grid%t0ml,HML=grid%hml,H0ML=grid%h0ml,HUML=grid%huml,HVML=grid%hvml,F=grid%f,TMOML=grid%TMOML,ISWATER=iswater,lakedepth2d=gr&
&id%lakedepth2d,savedtke12d=grid%savedtke12d,snowdp2d=grid%snowdp2d,h2osno2d=grid%h2osno2d,snl2d=grid%snl2d,t_grnd2d=grid%t_grnd2d,&
&t_lake3d=grid%t_lake3d,lake_icefrac3d=grid%lake_icefrac3d,z_lake3d=grid%z_lake3d,dz_lake3d=grid%dz_lake3d,t_soisno3d=grid%t_soisno&
&3d,h2osoi_ice3d=grid%h2osoi_ice3d,h2osoi_liq3d=grid%h2osoi_liq3d,h2osoi_vol3d=grid%h2osoi_vol3d,z3d=grid%z3d,dz3d=grid%dz3d,zi3d=g&
&rid%zi3d,watsat3d=grid%watsat3d,csol3d=grid%csol3d,tkmg3d=grid%tkmg3d,tkdry3d=grid%tkdry3d,tksatu3d=grid%tksatu3d,LakeModel=grid%s&
&f_lake_physics,lake_min_elev=grid%lake_min_elev,LakeMask=grid%LakeMask,NUMC=grid%numc,NUMP=grid%nump,SABV=grid%sabv,SABG=grid%sabg&
&,LWUP=grid%lwup,SNL=grid%snl,HISTORY_INTERVAL=config_flags%history_interval,SNOWDP=grid%snowdp,WTC=grid%wtc,WTP=grid%wtp,H2OSNO=gr&
&id%h2osno,T_GRND=grid%t_grnd,T_VEG=grid%t_veg,H2OCAN=grid%h2ocan,H2OCAN_COL=grid%h2ocan_col,T2M_MAX=grid%t2m_max,T2M_MIN=grid%t2m_&
&min,T2CLM=grid%t2clm,T_REF2M=grid%t_ref2m,H2OSOI_LIQ_S1=grid%h2osoi_liq_s1,H2OSOI_LIQ_S2=grid%h2osoi_liq_s2,H2OSOI_LIQ_S3=grid%h2o&
&soi_liq_s3,H2OSOI_LIQ_S4=grid%h2osoi_liq_s4,H2OSOI_LIQ_S5=grid%h2osoi_liq_s5,H2OSOI_LIQ1=grid%h2osoi_liq1,H2OSOI_LIQ2=grid%h2osoi_&
&liq2,H2OSOI_LIQ3=grid%h2osoi_liq3,H2OSOI_LIQ4=grid%h2osoi_liq4,H2OSOI_LIQ5=grid%h2osoi_liq5,H2OSOI_LIQ6=grid%h2osoi_liq6,H2OSOI_LI&
&Q7=grid%h2osoi_liq7,H2OSOI_LIQ8=grid%h2osoi_liq8,H2OSOI_LIQ9=grid%h2osoi_liq9,H2OSOI_LIQ10=grid%h2osoi_liq10,H2OSOI_ICE_S1=grid%h2&
&osoi_ice_s1,H2OSOI_ICE_S2=grid%h2osoi_ice_s2,H2OSOI_ICE_S3=grid%h2osoi_ice_s3,H2OSOI_ICE_S4=grid%h2osoi_ice_s4,H2OSOI_ICE_S5=grid%&
&h2osoi_ice_s5,H2OSOI_ICE1=grid%h2osoi_ice1,H2OSOI_ICE2=grid%h2osoi_ice2,H2OSOI_ICE3=grid%h2osoi_ice3,H2OSOI_ICE4=grid%h2osoi_ice4,&
&H2OSOI_ICE5=grid%h2osoi_ice5,H2OSOI_ICE6=grid%h2osoi_ice6,H2OSOI_ICE7=grid%h2osoi_ice7,H2OSOI_ICE8=grid%h2osoi_ice8,H2OSOI_ICE9=gr&
&id%h2osoi_ice9,H2OSOI_ICE10=grid%h2osoi_ice10,T_SOISNO_S1=grid%t_soisno_s1,T_SOISNO_S2=grid%t_soisno_s2,T_SOISNO_S3=grid%t_soisno_&
&s3,T_SOISNO_S4=grid%t_soisno_s4,T_SOISNO_S5=grid%t_soisno_s5,T_SOISNO1=grid%t_soisno1,T_SOISNO2=grid%t_soisno2,T_SOISNO3=grid%t_so&
&isno3,T_SOISNO4=grid%t_soisno4,T_SOISNO5=grid%t_soisno5,T_SOISNO6=grid%t_soisno6,T_SOISNO7=grid%t_soisno7,T_SOISNO8=grid%t_soisno8&
&,T_SOISNO9=grid%t_soisno9,T_SOISNO10=grid%t_soisno10,DZSNOW1=grid%dzsnow1,DZSNOW2=grid%dzsnow2,DZSNOW3=grid%dzsnow3,DZSNOW4=grid%d&
&zsnow4,DZSNOW5=grid%dzsnow5,SNOWRDS1=grid%snowrds1,SNOWRDS2=grid%snowrds2,SNOWRDS3=grid%snowrds3,SNOWRDS4=grid%snowrds4,SNOWRDS5=g&
&rid%snowrds5,T_LAKE1=grid%t_lake1,T_LAKE2=grid%t_lake2,T_LAKE3=grid%t_lake3,T_LAKE4=grid%t_lake4,T_LAKE5=grid%t_lake5,T_LAKE6=grid&
&%t_lake6,T_LAKE7=grid%t_lake7,T_LAKE8=grid%t_lake8,T_LAKE9=grid%t_lake9,T_LAKE10=grid%t_lake10,H2OSOI_VOL1=grid%h2osoi_vol1,H2OSOI&
&_VOL2=grid%h2osoi_vol2,H2OSOI_VOL3=grid%h2osoi_vol3,H2OSOI_VOL4=grid%h2osoi_vol4,H2OSOI_VOL5=grid%h2osoi_vol5,H2OSOI_VOL6=grid%h2o&
&soi_vol6,H2OSOI_VOL7=grid%h2osoi_vol7,H2OSOI_VOL8=grid%h2osoi_vol8,H2OSOI_VOL9=grid%h2osoi_vol9,H2OSOI_VOL10=grid%h2osoi_vol10,MAX&
&PATCH=config_flags%maxpatch,INEST=grid%id,ALBEDOsubgrid=grid%ALBEDOsubgrid,LHsubgrid=grid%LHsubgrid,HFXsubgrid=grid%HFXsubgrid,LWU&
&Psubgrid=grid%LWUPsubgrid,Q2subgrid=grid%Q2subgrid,SABVsubgrid=grid%SABVsubgrid,SABGsubgrid=grid%SABGsubgrid,NRAsubgrid=grid%NRAsu&
&bgrid,SWUPsubgrid=grid%SWUPsubgrid,LHsoi=grid%LHsoi,LHveg=grid%LHveg,LHtran=grid%LHtran,SLOPE_RAD=config_flags%slope_rad,TOPO_SHAD&
&ING=config_flags%topo_shading,SHADOWMASK=grid%shadowmask,SLOPE=grid%slope,SLP_AZI=grid%slp_azi,SWNORM=grid%swnorm,DECLIN=grid%decl&
&in,SOLCON=grid%solcon,COSZEN=grid%coszen,HRANG=grid%hrang,xlat_urb2d=grid%XLAT,NUM_ROOF_LAYERS=num_roof_layers,NUM_WALL_LAYERS=num&
&_wall_layers,NUM_ROAD_LAYERS=num_road_layers,DZR=grid%dzr,DZB=grid%dzb,DZG=grid%dzg,TR_URB2D=grid%tr_urb2d,TB_URB2D=grid%tb_urb2d,&
&TG_URB2D=grid%tg_urb2d,TC_URB2D=grid%tc_urb2d,QC_URB2D=grid%qc_urb2d,UC_URB2D=grid%uc_urb2d,XXXR_URB2D=grid%xxxr_urb2d,XXXB_URB2D=&
&grid%xxxb_urb2d,XXXG_URB2D=grid%xxxg_urb2d,XXXC_URB2D=grid%xxxc_urb2d,TRL_URB3D=grid%trl_urb3d,TBL_URB3D=grid%tbl_urb3d,TGL_URB3D=&
&grid%tgl_urb3d,SH_URB2D=grid%sh_urb2d,LH_URB2D=grid%lh_urb2d,G_URB2D=grid%g_urb2d,RN_URB2D=grid%rn_urb2d,TS_URB2D=grid%ts_urb2d,FR&
&C_URB2D=grid%frc_urb2d,UTYPE_URB2D=grid%utype_urb2d,SF_URBAN_PHYSICS=config_flags%sf_urban_physics,NUM_URBAN_LAYERS=config_flags%n&
&um_urban_layers,NUM_URBAN_HI=config_flags%num_urban_hi,TRB_URB4D=grid%trb_urb4d,TW1_URB4D=grid%tw1_urb4d,TW2_URB4D=grid%tw2_urb4d,&
&TGB_URB4D=grid%tgb_urb4d,TLEV_URB3D=grid%tlev_urb3d,QLEV_URB3D=grid%qlev_urb3d,TW1LEV_URB3D=grid%tw1lev_urb3d,TW2LEV_URB3D=grid%tw&
&2lev_urb3d,TGLEV_URB3D=grid%tglev_urb3d,TFLEV_URB3D=grid%tflev_urb3d,SF_AC_URB3D=grid%sf_ac_urb3d,LF_AC_URB3D=grid%lf_ac_urb3d,CM_&
&AC_URB3D=grid%cm_ac_urb3d,SFVENT_URB3D=grid%sfvent_urb3d,LFVENT_URB3D=grid%lfvent_urb3d,SFWIN1_URB3D=grid%sfwin1_urb3d,SFWIN2_URB3&
&D=grid%sfwin2_urb3d,SFW1_URB3D=grid%sfw1_urb3d,SFW2_URB3D=grid%sfw2_urb3d,SFR_URB3D=grid%sfr_urb3d,SFG_URB3D=grid%sfg_urb3d,LP_URB&
&2D=grid%lp_urb2d,HI_URB2D=grid%hi_urb2d,LB_URB2D=grid%lb_urb2d,HGT_URB2D=grid%hgt_urb2d,MH_URB2D=grid%mh_urb2d,STDH_URB2D=grid%std&
&h_urb2d,LF_URB2D=grid%lf_urb2d,GMT=grid%gmt,XLAT=grid%xlat,XLONG=grid%xlong,JULDAY=grid%julday,A_U_BEP=grid%a_u_bep,A_V_BEP=grid%a&
&_v_bep,A_T_BEP=grid%a_t_bep,A_Q_BEP=grid%a_q_bep,B_U_BEP=grid%b_u_bep,B_V_BEP=grid%b_v_bep,B_T_BEP=grid%b_t_bep,B_Q_BEP=grid%b_q_b&
&ep,SF_BEP=grid%sf_bep,VL_BEP=grid%vl_bep,A_E_BEP=grid%a_e_bep,B_E_BEP=grid%b_e_bep,DLG_BEP=grid%dlg_bep,DL_U_BEP=grid%dl_u_bep,CMR&
&_SFCDIF=grid%cmr_sfcdif,CHR_SFCDIF=grid%chr_sfcdif,CMC_SFCDIF=grid%cmc_sfcdif,CHC_SFCDIF=grid%chc_sfcdif,LANDUSEF=grid%landusef,SO&
&ILCTOP=grid%soilctop,SOILCBOT=grid%soilcbot,RA=grid%ra,RS=grid%rs,LAI=grid%lai,NLCAT=grid%num_land_cat,NSCAT=grid%num_soil_cat,VEG&
&F_PX=grid%vegf_px,SNOWNCV=grid%snowncv,ANAL_INTERVAL=config_flags%auxinput9_interval_s+config_flags%auxinput9_interval_m*60,PXLSM_&
&SMOIS_INIT=config_flags%pxlsm_smois_init,PXLSM_SOIL_NUDGE=config_flags%pxlsm_soil_nudge,alswvisdir=grid%alswvisdir,alswvisdif=grid&
&%alswvisdif,alswnirdir=grid%alswnirdir,alswnirdif=grid%alswnirdif,swvisdir=grid%swvisdir,swvisdif=grid%swvisdif,swnirdir=grid%swni&
&rdir,swnirdif=grid%swnirdif,ssib_br=grid%ssib_br,ssib_fm=grid%ssib_fm,ssib_fh=grid%ssib_fh,ssib_cm=grid%ssib_cm,ssibxdd=grid%ssibx&
&dd,ssib_lhf=grid%ssib_lhf,ssib_shf=grid%ssib_shf,ssib_ghf=grid%ssib_ghf,ssib_egs=grid%ssib_egs,ssib_eci=grid%ssib_eci,ssib_ect=gri&
&d%ssib_ect,ssib_egi=grid%ssib_egi,ssib_egt=grid%ssib_egt,ssib_sdn=grid%ssib_sdn,ssib_sup=grid%ssib_sup,ssib_ldn=grid%ssib_ldn,ssib&
&_lup=grid%ssib_lup,ssib_wat=grid%ssib_wat,ssib_shc=grid%ssib_shc,ssib_shg=grid%ssib_shg,ssib_lai=grid%ssib_lai,ssib_vcf=grid%ssib_&
&vcf,ssib_z00=grid%ssib_z00,ssib_veg=grid%ssib_veg,cldfra=grid%cldfra,ISNOW=grid%isnow,SWE=grid%swe,SNOWDEN=grid%snowden,SNOWDEPTH=&
&grid%snowdepth,TKAIR=grid%tkair,DZO1=grid%dzo1,WO1=grid%wo1,TSSN1=grid%tssn1,TSSNO1=grid%tssno1,BWO1=grid%bwo1,BTO1=grid%bto1,CTO1&
&=grid%cto1,FIO1=grid%fio1,FLO1=grid%flo1,BIO1=grid%bio1,BLO1=grid%blo1,HO1=grid%ho1,DZO2=grid%dzo2,WO2=grid%wo2,TSSN2=grid%tssn2,T&
&SSNO2=grid%tssno2,BWO2=grid%bwo2,BTO2=grid%bto2,CTO2=grid%cto2,FIO2=grid%fio2,FLO2=grid%flo2,BIO2=grid%bio2,BLO2=grid%blo2,HO2=gri&
&d%ho2,DZO3=grid%dzo3,WO3=grid%wo3,TSSN3=grid%tssn3,TSSNO3=grid%tssno3,BWO3=grid%bwo3,BTO3=grid%bto3,CTO3=grid%cto3,FIO3=grid%fio3,&
&FLO3=grid%flo3,BIO3=grid%bio3,BLO3=grid%blo3,HO3=grid%ho3,DZO4=grid%dzo4,WO4=grid%wo4,TSSN4=grid%tssn4,TSSNO4=grid%tssno4,BWO4=gri&
&d%bwo4,BTO4=grid%bto4,CTO4=grid%cto4,FIO4=grid%fio4,FLO4=grid%flo4,BIO4=grid%bio4,BLO4=grid%blo4,HO4=grid%ho4,RA_SW_PHYSICS=config&
&_flags%ra_sw_physics,t2_ndg_old=grid%t2_ndg_old,q2_ndg_old=grid%q2_ndg_old,t2_ndg_new=grid%t2_ndg_new,q2_ndg_new=grid%q2_ndg_new,s&
&n_ndg_old=grid%sn_ndg_old,sn_ndg_new=grid%sn_ndg_new,idveg=config_flags%dveg,iopt_crs=config_flags%opt_crs,iopt_btr=config_flags%o&
&pt_btr,iopt_run=config_flags%opt_run,iopt_sfc=config_flags%opt_sfc,iopt_frz=config_flags%opt_frz,iopt_inf=config_flags%opt_inf,iop&
&t_rad=config_flags%opt_rad,iopt_alb=config_flags%opt_alb,iopt_snf=config_flags%opt_snf,iopt_tbot=config_flags%opt_tbot,iopt_stc=co&
&nfig_flags%opt_stc,isnowxy=grid%isnowxy,tvxy=grid%tvxy,tgxy=grid%tgxy,canicexy=grid%canicexy,canliqxy=grid%canliqxy,eahxy=grid%eah&
&xy,tahxy=grid%tahxy,cmxy=grid%cmxy,chxy=grid%chxy,fwetxy=grid%fwetxy,sneqvoxy=grid%sneqvoxy,alboldxy=grid%alboldxy,qsnowxy=grid%qs&
&nowxy,wslakexy=grid%wslakexy,zwtxy=grid%zwtxy,waxy=grid%waxy,wtxy=grid%wtxy,tsnoxy=grid%tsnoxy,zsnsoxy=grid%zsnsoxy,snicexy=grid%s&
&nicexy,snliqxy=grid%snliqxy,lfmassxy=grid%lfmassxy,rtmassxy=grid%rtmassxy,stmassxy=grid%stmassxy,woodxy=grid%woodxy,stblcpxy=grid%&
&stblcpxy,fastcpxy=grid%fastcpxy,xsaixy=grid%xsaixy,taussxy=grid%taussxy,t2mvxy=grid%t2mvxy,t2mbxy=grid%t2mbxy,q2mvxy=grid%q2mvxy,q&
&2mbxy=grid%q2mbxy,tradxy=grid%tradxy,neexy=grid%neexy,gppxy=grid%gppxy,nppxy=grid%nppxy,fvegxy=grid%fvegxy,runsfxy=grid%runsfxy,ru&
&nsbxy=grid%runsbxy,ecanxy=grid%ecanxy,edirxy=grid%edirxy,etranxy=grid%etranxy,fsaxy=grid%fsaxy,firaxy=grid%firaxy,aparxy=grid%apar&
&xy,psnxy=grid%psnxy,savxy=grid%savxy,sagxy=grid%sagxy,rssunxy=grid%rssunxy,rsshaxy=grid%rsshaxy,bgapxy=grid%bgapxy,wgapxy=grid%wga&
&pxy,tgvxy=grid%tgvxy,tgbxy=grid%tgbxy,chvxy=grid%chvxy,chbxy=grid%chbxy,shgxy=grid%shgxy,shcxy=grid%shcxy,shbxy=grid%shbxy,evgxy=g&
&rid%evgxy,evbxy=grid%evbxy,ghvxy=grid%ghvxy,ghbxy=grid%ghbxy,irgxy=grid%irgxy,ircxy=grid%ircxy,irbxy=grid%irbxy,trxy=grid%trxy,evc&
&xy=grid%evcxy,chleafxy=grid%chleafxy,chucxy=grid%chucxy,chv2xy=grid%chv2xy,chb2xy=grid%chb2xy,chstarxy=grid%chstarxy,smcwtdxy=grid&
&%smcwtdxy,rechxy=grid%rechxy,deeprechxy=grid%deeprechxy,fdepthxy=grid%fdepthxy,areaxy=grid%areaxy,rivercondxy=grid%rivercondxy,riv&
&erbedxy=grid%riverbedxy,eqzwt=grid%eqzwt,pexpxy=grid%pexpxy,qrfxy=grid%qrfxy,qspringxy=grid%qspringxy,qslatxy=grid%qslatxy,qrfsxy=&
&grid%qrfsxy,qspringsxy=grid%qspringsxy,smoiseq=grid%smoiseq,wtddt=config_flags%wtddt,stepwtd=grid%stepwtd,ua_phys=config_flags%ua_&
&phys,flx4=grid%flx4,fvb=grid%fvb,fbur=grid%fbur,fgsn=grid%fgsn,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde,IMS=ims,IME=ime,JMS&
&=jms,JME=jme,KMS=kms,KME=kme,IPS=ips,IPE=ipe,JPS=jps,JPE=jpe,KPS=kps,KPE=kpe,I_START=grid%i_start,I_END=min(grid%i_end,ide-1),J_ST&
&ART=grid%j_start,J_END=min(grid%j_end,jde-1),KTS=k_start,KTE=min(k_end,kde-1),NUM_TILES=grid%num_tiles,te_temf=grid%te_temf,hd_tem&
&f=grid%hd_temf,fCor=grid%f,exch_temf=exch_temf,wm_temf=grid%wm_temf,hfx_force=grid%hfx_force,lh_force=grid%lh_force,tsk_force=grid&
&%tsk_force,hfx_force_tend=grid%hfx_force_tend,lh_force_tend=grid%lh_force_tend,tsk_force_tend=grid%tsk_force_tend,QV_CURR=moist(im&
&s,kms,jms,P_QV),F_QV=F_QV,QC_CURR=moist(ims,kms,jms,P_QC),F_QC=F_QC,QR_CURR=moist(ims,kms,jms,P_QR),F_QR=F_QR,QI_CURR=moist(ims,km&
&s,jms,P_QI),F_QI=F_QI,QS_CURR=moist(ims,kms,jms,P_QS),F_QS=F_QS,QG_CURR=moist(ims,kms,jms,P_QG),F_QG=F_QG,CAPG=grid%capg,EMISS=gri&
&d%emiss,HOL=hol,MOL=grid%mol,T2OBS=grid%t2obs,Q2OBS=grid%q2obs,RAINBL=grid%rainbl,SR=grid%sr,RAINSHV=grid%rainshv,RAINNCV=grid%rai&
&nncv,REGIME=grid%regime,T2=grid%t2,THC=grid%thc,QSG=grid%qsg,QVG=grid%qvg,QCG=grid%qcg,SOILT1=grid%soilt1,TSNAV=grid%tsnav,SMFR3D=&
&grid%smfr3d,KEEPFR3DFLAG=grid%keepfr3dflag,DEW=grid%dew,POTEVP=grid%POTEVP,SNOPCX=grid%SNOPCX,SOILTB=grid%SOILTB,MOSAIC_LU=config_&
&flags%mosaic_lu,MOSAIC_SOIL=config_flags%mosaic_soil,ISURBAN=isurban,MMINLU=TRIM(mminlu),SNOTIME=grid%SNOTIME,RDLAI2D=config_flags&
&%rdlai2d,usemonalb=config_flags%usemonalb,NOAHRES=grid%noahres,TSK_SAVE=grid%tsk_save,ch=grid%ch,tsq=grid%tsq,qsq=grid%qsq,cov=gri&
&d%cov,Sh3d=grid%sh3d,EL_PBL=grid%el_pbl,bl_mynn_cloudpdf=config_flags%bl_mynn_cloudpdf,fgdp=grid%fgdp,dfgdp=grid%dfgdp,vdfg=grid%v&
&dfg,grav_settling=config_flags%grav_settling,OM_TMP=grid%om_tmp,OM_S=grid%om_s,OM_U=grid%om_u,OM_V=grid%om_v,OM_DEPTH=grid%om_dept&
&h,OM_ML=grid%OM_ML,OM_LON=grid%om_lon,OM_LAT=grid%om_lat,okms=1,okme=config_flags%ocean_levels,rdx=grid%rdx,rdy=grid%rdy,msfu=grid&
&%msfu,msfv=grid%msfv,msft=grid%msft,XTIME=grid%xtime,OM_TINI=grid%om_tini,OM_SINI=grid%om_sini,id=grid%id,omdt=config_flags%omdt,s&
&f_surface_mosaic=config_flags%sf_surface_mosaic,mosaic_cat=config_flags%mosaic_cat,mosaic_cat_index=grid%mosaic_cat_index,landusef&
&2=grid%landusef2,TSK_mosaic=grid%TSK_mosaic,QSFC_mosaic=grid%QSFC_mosaic,TSLB_mosaic=grid%TSLB_mosaic,SMOIS_mosaic=grid%SMOIS_mosa&
&ic,SH2O_mosaic=grid%SH2O_mosaic,CANWAT_mosaic=grid%CANWAT_mosaic,SNOW_mosaic=grid%SNOW_mosaic,SNOWH_mosaic=grid%SNOWH_mosaic,SNOWC&
&_mosaic=grid%SNOWC_mosaic,ALBEDO_mosaic=grid%ALBEDO_mosaic,ALBBCK_mosaic=grid%ALBBCK_mosaic,EMISS_mosaic=grid%EMISS_mosaic,EMBCK_m&
&osaic=grid%EMBCK_mosaic,ZNT_mosaic=grid%ZNT_mosaic,Z0_mosaic=grid%Z0_mosaic,HFX_mosaic=grid%HFX_mosaic,QFX_mosaic=grid%QFX_mosaic,&
&LH_mosaic=grid%LH_mosaic,GRDFLX_mosaic=grid%GRDFLX_mosaic,SNOTIME_mosaic=grid%SNOTIME_mosaic,TR_URB2D_mosaic=grid%TR_URB2D_mosaic,&
&TB_URB2D_mosaic=grid%TB_URB2D_mosaic,TG_URB2D_mosaic=grid%TG_URB2D_mosaic,TC_URB2D_mosaic=grid%TC_URB2D_mosaic,QC_URB2D_mosaic=gri&
&d%QC_URB2D_mosaic,UC_URB2D_mosaic=grid%UC_URB2D_mosaic,TRL_URB3D_mosaic=grid%TRL_URB3D_mosaic,TBL_URB3D_mosaic=grid%TBL_URB3D_mosa&
&ic,TGL_URB3D_mosaic=grid%TGL_URB3D_mosaic,SH_URB2D_mosaic=grid%SH_URB2D_mosaic,LH_URB2D_mosaic=grid%LH_URB2D_mosaic,G_URB2D_mosaic&
&=grid%G_URB2D_mosaic,RN_URB2D_mosaic=grid%RN_URB2D_mosaic,TS_URB2D_mosaic=grid%TS_URB2D_mosaic,TS_RUL2D_mosaic=grid%TS_RUL2D_mosai&
&c)







      CALL wrf_debug ( 200 , ' call pbl_driver' )

CALL pbl_driver(AKHS=grid%akhs,AKMS=grid%akms,BL_PBL_PHYSICS=config_flags%bl_pbl_physics,WINDFARM_OPT=config_flags%windfarm_opt,po&
&wer=grid%power,BLDT=grid%bldt,CURR_SECS=curr_secs,ADAPT_STEP_FLAG=adapt_step_flag,BLDTACTTIME=grid%bldtacttime,BR=br,CHKLOWQ=chklo&
&wq,CT=grid%ct,DT=grid%dt,DX=grid%dx,DZ8W=dz8w,EXCH_H=grid%exch_h,EXCH_M=grid%exch_m,FM=grid%fm,FHH=grid%fh,F=grid%f,GRDFLX=grid%gr&
&dflx,GZ1OZ0=gz1oz0,HFX=grid%hfx,HT=grid%ht,ID=grid%id,ITIMESTEP=grid%itimestep,KPBL=grid%kpbl,LH=grid%lh,LOWLYR=grid%lowlyr,P8W=gr&
&id%p_hyd_w,PBLH=grid%pblh,PI_PHY=pi_phy,PSIH=psih,PSIM=psim,P_PHY=grid%p_hyd,QFX=grid%qfx,QSFC=grid%qsfc,QZ0=grid%qz0,MIXHT=mixht,&
&RA_LW_PHYSICS=config_flags%ra_lw_physics,RHO=grid%rho,RQCBLTEN=grid%rqcblten,RQIBLTEN=grid%rqiblten,RQVBLTEN=grid%rqvblten,RTHBLTE&
&N=grid%rthblten,RUBLTEN=grid%rublten,RVBLTEN=grid%rvblten,SNOW=grid%snow,STEPBL=grid%stepbl,THZ0=grid%thz0,TH_PHY=th_phy,TSK=grid%&
&tsk,T_PHY=grid%t_phy,UST=grid%ust,U10=grid%u10,UZ0=grid%uz0,U_FRAME=grid%u_frame,U_PHY=grid%u_phy,V10=grid%v10,VZ0=grid%vz0,V_FRAM&
&E=grid%v_frame,V_PHY=grid%v_phy,UOCE=grid%uoce,VOCE=grid%voce,T2=grid%t2,WARM_RAIN=grid%warm_rain,WSPD=wspd,XICE=grid%xice,XLAND=g&
&rid%xland,Z=grid%z,ZNT=grid%znt,CTOPO=grid%ctopo,CTOPO2=grid%ctopo2,FRC_URB2D=grid%frc_urb2d,A_U_BEP=grid%a_u_bep,A_V_BEP=grid%a_v&
&_bep,A_T_BEP=grid%a_t_bep,A_Q_BEP=grid%a_q_bep,B_U_BEP=grid%b_u_bep,B_V_BEP=grid%b_v_bep,B_T_BEP=grid%b_t_bep,B_Q_BEP=grid%b_q_bep&
&,SF_BEP=grid%sf_bep,VL_BEP=grid%vl_bep,A_E_BEP=grid%a_e_bep,B_E_BEP=grid%b_e_bep,DLG_BEP=grid%dlg_bep,DL_U_BEP=grid%dl_u_bep,SF_SF&
&CLAY_PHYSICS=config_flags%sf_sfclay_physics,SF_URBAN_PHYSICS=config_flags%sf_urban_physics,TKE_PBL=grid%tke_pbl,EL_PBL=grid%el_pbl&
&,WU_TUR=grid%wu_tur,WV_tur=grid%wv_tur,WT_tur=grid%wt_tur,WQ_tur=grid%wq_tur,EXCH_TKE=grid%exch_tke,RTHRATEN=grid%rthraten,IDS=ids&
&,IDE=ide,JDS=jds,JDE=jde,KDS=kds,KDE=kde,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme,I_START=grid%i_start,I_END=min(grid%i_end&
&,ide-1),J_START=grid%j_start,J_END=min(grid%j_end,jde-1),KTS=k_start,KTE=min(k_end,kde-1),NUM_TILES=grid%num_tiles,ZNU=grid%znu,ZN&
&W=grid%znw,MUT=grid%mut,P_TOP=grid%p_top,te_temf=grid%te_temf,kh_temf=grid%kh_temf,km_temf=grid%km_temf,shf_temf=grid%shf_temf,qf_&
&temf=grid%qf_temf,uw_temf=grid%uw_temf,vw_temf=grid%vw_temf,hd_temf=grid%hd_temf,lcl_temf=grid%lcl_temf,wupd_temf=grid%wupd_temf,m&
&f_temf=grid%mf_temf,thup_temf=grid%thup_temf,qtup_temf=grid%qtup_temf,qlup_temf=grid%qlup_temf,cf3d_temf=grid%cf3d_temf,cfm_temf=g&
&rid%cfm_temf,hct_temf=grid%hct_temf,flhc=grid%flhc,flqc=grid%flqc,exch_temf=exch_temf,QV_CURR=moist(ims,kms,jms,P_QV),F_QV=F_QV,QC&
&_CURR=moist(ims,kms,jms,P_QC),F_QC=F_QC,QR_CURR=moist(ims,kms,jms,P_QR),F_QR=F_QR,QI_CURR=moist(ims,kms,jms,P_QI),F_QI=F_QI,QS_CUR&
&R=moist(ims,kms,jms,P_QS),F_QS=F_QS,QG_CURR=moist(ims,kms,jms,P_QG),F_QG=F_QG,HOL=HOL,MOL=grid%mol,REGIME=grid%REGIME,QKE=grid%qke&
&,Sh3d=grid%sh3d,QKE_ADV=scalar(ims,kms,jms,P_qke_adv),bl_mynn_tkeadvect=config_flags%bl_mynn_tkeadvect,tsq=grid%tsq,qsq=grid%qsq,c&
&ov=grid%cov,DQKE=grid%dqke,QWT=grid%qWT,QSHEAR=grid%qSHEAR,QBUOY=grid%qBUOY,QDISS=grid%qDISS,bl_mynn_tkebudget=config_flags%bl_myn&
&n_tkebudget,bl_mynn_cloudpdf=config_flags%bl_mynn_cloudpdf,rmol=grid%rmol,ch=grid%ch,qcg=grid%qcg,grav_settling=config_flags%grav_&
&settling,vdfg=grid%vdfg,GWD_OPT=config_flags%gwd_opt,DTAUX3D=grid%dtaux3d,DTAUY3D=grid%dtauy3d,DUSFCG=grid%dusfcg,DVSFCG=grid%dvsf&
&cg,VAR2D=grid%var2d,OC12D=grid%oc12d,OA1=grid%oa1,OA2=grid%oa2,OA3=grid%oa3,OA4=grid%oa4,OL1=grid%ol1,OL2=grid%ol2,OL3=grid%ol3,OL&
&4=grid%ol4,MFSHCONV=config_flags%mfshconv,MASSFLUX_EDKF=grid%massflux_EDKF,ENTR_EDKF=grid%entr_EDKF,DETR_EDKF=grid%detr_EDKF,THL_U&
&P=grid%thl_up,THV_UP=grid%thv_up,RT_UP=grid%rt_up,RV_UP=grid%rv_up,RC_UP=grid%rc_up,U_UP=grid%u_up,V_UP=grid%v_up,FRAC_UP=grid%fra&
&c_up,RC_MF=grid%RC_MF,phb=grid%phb,XLAT_U=grid%xlat_u,XLONG_U=grid%xlong_u,Z_AT_W=grid%z_at_w,CLDFRA_OLD_MP=grid%cldfra_old_mp,CLD&
&FRA=grid%cldfra,RTHRATENLW=grid%rthratenlw,TAURESX2D=grid%tauresx2d,TAURESY2D=grid%tauresy2d,TPERT2D=grid%tpert2d,QPERT2D=grid%qpe&
&rt2d,WPERT2D=grid%wpert2d,WSEDL3D=grid%wsedl3d,TURBTYPE3D=grid%turbtype3d,SMAW3D=grid%smaw3d,QNC_CURR=scalar(ims,kms,jms,P_QNC),F_&
&QNC=f_qnc,QNI_CURR=scalar(ims,kms,jms,P_QNI),F_QNI=f_qni,RQNIBLTEN=grid%rqniblten,XLAT_V=grid%xlat_v,XLONG_V=grid%xlong_v,FNM=grid&
&%fnm,FNP=grid%fnp,IS_CAMMGMP_USED=grid%is_CAMMGMP_used,WSTAR=grid%wstar_ysu,DELTA=grid%delta_ysu,SCALAR=scalar,SCALAR_TEND=scalar_&
&tend,NUM_SCALAR=num_scalar,TRACER=tracer,TRACER_TEND=tracer_tend,NUM_TRACER=num_tracer,SCALAR_PBLMIX=config_flags%scalar_pblmix,TR&
&ACER_PBLMIX=config_flags%tracer_pblmix)
                                                         








      IF ((grid%sr_x > 0 .OR. grid%sr_y > 0) .AND. config_flags%ifire == 2) THEN


        if(config_flags%ifire.eq.2)then 
            





            
            call fire_driver_em_step (  grid , config_flags    &
            ,ids,ide, kds,kde, jds,jde                              &
            ,ims,ime, kms,kme, jms,jme                              &
            ,ips,ipe, kps,kpe, jps,jpe                              &
            ,grid%rho,grid%z_at_w,dz8w)
        endif


       ENDIF




      CALL wrf_debug ( 200 , ' call cumulus_driver' )








CALL cumulus_driver(grid,U=grid%u_phy,V=grid%v_phy,TH=th_phy,T=grid%t_phy,W=grid%w_2,P=grid%p_hyd,PI=pi_phy,RHO=grid%rho,ITIMESTEP&
&=grid%itimestep,DT=grid%dt,DX=grid%dx,CUDT=grid%cudt,CURR_SECS=curr_secs,ADAPT_STEP_FLAG=adapt_step_flag,CUDTACTTIME=grid%cudtactt&
&ime,RAINC=grid%rainc,RAINCV=grid%raincv,PRATEC=grid%pratec,NCA=grid%nca,CLDFRA_DP=grid%cldfra_dp,CLDFRA_SH=grid%cldfra_sh,QC_CU=gr&
&id%QC_CU,QI_CU=grid%QI_CU,HTOP=grid%cutop,HBOT=grid%cubot,KPBL=grid%kpbl,Z=grid%z,Z_AT_W=grid%z_at_w,MAVAIL=grid%mavail,PBLH=grid%&
&pblh,DZ8W=dz8w,P8W=grid%p_hyd_w,PSFC=grid%psfc,TSK=grid%tsk,TKE_PBL=grid%tke_pbl,UST=grid%ust,W0AVG=grid%w0avg,STEPCU=grid%stepcu,&
&CLDEFI=grid%cldefi,LOWLYR=grid%lowlyr,XLAND=grid%xland,APR_GR=grid%apr_gr,APR_W=grid%apr_w,APR_MC=grid%apr_mc,APR_ST=grid%apr_st,A&
&PR_AS=grid%apr_as,APR_CAPMA=grid%apr_capma,APR_CAPME=grid%apr_capme,APR_CAPMI=grid%apr_capmi,MASS_FLUX=grid%mass_flux,XF_ENS=grid%&
&xf_ens,PR_ENS=grid%pr_ens,HT=grid%ht,EDT_OUT=grid%edt_out,imomentum=grid%imomentum,clos_choice=grid%clos_choice,ishallow=config_fl&
&ags%ishallow,cugd_tten=grid%cugd_tten,cugd_qvten=grid%cugd_qvten,cugd_qcten=grid%cugd_qcten,cugd_ttens=grid%cugd_ttens,cugd_qvtens&
&=grid%cugd_qvtens,ENSDIM=config_flags%ensdim,MAXIENS=config_flags%maxiens,MAXENS=config_flags%maxens,MAXENS2=config_flags%maxens2,&
&MAXENS3=config_flags%maxens3,CU_ACT_FLAG=cu_act_flag,WARM_RAIN=grid%warm_rain,HFX=grid%hfx,QFX=grid%qfx,CLDFRA=grid%cldfra,CLDFRA_&
&MP_ALL=grid%cldfra_mp_all,TPERT2D=grid%tpert2d,GSW=grid%gsw,cugd_avedx=config_flags%cugd_avedx,k22_shallow=grid%k22_shallow,kbcon_&
&shallow=grid%kbcon_shallow,ktop_shallow=grid%ktop_shallow,xmb_shallow=grid%xmb_shallow,ktop_deep=grid%ktop_deep,PERIODIC_X=(config&
&_flags%polar.OR.config_flags%periodic_x),PERIODIC_Y=config_flags%periodic_y,IS_CAMMGMP_USED=grid%is_CAMMGMP_used,EVAPCDP3D=grid%ev&
&apcdp3d,ICWMRDP3D=grid%icwmrdp3d,RPRDDP3D=grid%rprddp3d,CAPE=grid%cape,ZMMU=grid%zmmu,ZMMD=grid%zmmd,ZMDT=grid%zmdt,ZMDQ=grid%zmdq&
&,DLF=grid%dlf,RLIQ=grid%rliq,PCONVB=grid%pconvb,PCONVT=grid%pconvt,EVAPTZM=grid%evaptzm,FZSNTZM=grid%fzsntzm,EVSNTZM=grid%evsntzm,&
&EVAPQZM=grid%evapqzm,ZMFLXPRC=grid%zmflxprc,ZMFLXSNW=grid%zmflxsnw,ZMNTPRPD=grid%zmntprpd,ZMNTSNPD=grid%zmntsnpd,ZMEIHEAT=grid%zme&
&iheat,CMFMC=grid%cmfmc,CMFMCDZM=grid%cmfmcdzm,PRECCDZM=grid%preccdzm,PRECZ=grid%precz,ZMMTU=grid%zmmtu,ZMMTV=grid%zmmtv,ZMUPGU=gri&
&d%zmupgu,ZMUPGD=grid%zmupgd,ZMVPGU=grid%zmvpgu,ZMVPGD=grid%zmvpgd,ZMICUU=grid%zmicuu,ZMICUD=grid%zmicud,ZMICVU=grid%zmicvu,ZMICVD=&
&grid%zmicvd,ZMDICE=grid%zmdice,ZMDLIQ=grid%zmdliq,dp3d=grid%dp3d,du3d=grid%du3d,ed3d=grid%ed3d,eu3d=grid%eu3d,md3d=grid%md3d,mu3d=&
&grid%mu3d,dsubcld2d=grid%dsubcld2d,ideep2d=grid%ideep2d,jt2d=grid%jt2d,maxg2d=grid%maxg2d,lengath2d=grid%lengath2d,pgcon=config_fl&
&ags%sas_pgcon,CU_PHYSICS=config_flags%cu_physics,BL_PBL_PHYSICS=config_flags%bl_pbl_physics,SF_SFCLAY_PHYSICS=config_flags%sf_sfcl&
&ay_physics,KFETA_TRIGGER=config_flags%kfeta_trigger,NSAS_DX_FACTOR=config_flags%nsas_dx_factor,IDS=ids,IDE=ide,JDS=jds,JDE=jde,KDS&
&=kds,KDE=kde,IMS=ims,IME=ime,JMS=jms,JME=jme,KMS=kms,KME=kme,IPS=ips,IPE=ipe,JPS=jps,JPE=jpe,KPS=kps,KPE=kpe,I_START=grid%i_start,&
&I_END=min(grid%i_end,ide-1),J_START=grid%j_start,J_END=min(grid%j_end,jde-1),KTS=k_start,KTE=min(k_end,kde-1),NUM_TILES=grid%num_t&
&iles,RQVCUTEN=grid%rqvcuten,RQCCUTEN=grid%rqccuten,RQSCUTEN=grid%rqscuten,RQICUTEN=grid%rqicuten,RQRCUTEN=grid%rqrcuten,RQCNCUTEN=&
&grid%rqcncuten,RQINCUTEN=grid%rqincuten,RQVBLTEN=grid%rqvblten,RQVFTEN=grid%rqvften,RTHRATEN=grid%rthraten,RTHBLTEN=grid%rthblten,&
&RUCUTEN=grid%rucuten,RVCUTEN=grid%rvcuten,RTHCUTEN=grid%rthcuten,RTHFTEN=grid%rthften,QV_CURR=moist(ims,kms,jms,P_QV),F_QV=F_QV,QC&
&_CURR=moist(ims,kms,jms,P_QC),F_QC=F_QC,QR_CURR=moist(ims,kms,jms,P_QR),F_QR=F_QR,QI_CURR=moist(ims,kms,jms,P_QI),F_QI=F_QI,QS_CUR&
&R=moist(ims,kms,jms,P_QS),F_QS=F_QS,QG_CURR=moist(ims,kms,jms,P_QG),F_QG=F_QG,ZNU=grid%znu,MP_PHYSICS=config_flags%mp_physics,GD_C&
&LOUD=grid%GD_CLOUD,GD_CLOUD2=grid%GD_CLOUD2,cfu1=grid%cfu1,cfd1=grid%cfd1,dfu1=grid%dfu1,efu1=grid%efu1,dfd1=grid%dfd1,efd1=grid%e&
&fd1,f_flux=l_flux)




     if(config_flags%cu_diag.eq.1)then
      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles
           call convtrans_prep(grid%gd_cloud,grid%gd_cloud2,grid%gd_cloud_a,&
     &          grid%QC_CU,grid%raincv,grid%raincv_a,grid%raincv_b,     &
     &          grid%gd_cldfr,moist,p_QV,p_QC,p_qi,T_PHY,P_PHY,num_moist,    &
     &          grid%gd_cloud2_a,grid%QI_CU,grid%convtrans_avglen_m,&
     &          adapt_step_flag,curr_secs,                &
     &          grid%itimestep,grid%dt,                                      &
     &          config_flags%cu_physics,                                     &
     &          ids,ide, jds,jde, kds,kde,                                   &
     &          ims,ime, jms,jme, kms,kme                                    &
     &          ,ITS=grid%i_start(ij),ITE=min(grid%i_end(ij), ide-1)                 &
     &          ,JTS=grid%j_start(ij),JTE=min(grid%j_end(ij), jde-1)                 &
     &          ,KTS=k_start, KTE=min(k_end,kde-1))
      ENDDO
      !$OMP END PARALLEL DO
     endif


      CALL wrf_debug ( 200 , ' call shallow_cumulus_driver' )


      CALL shallowcu_driver(                                              &
     &              IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde     &
     &             ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme     &
     &             ,IPS=ips,IPE=ipe, JPS=jps,JPE=jpe, KPS=kps,KPE=kpe     &
     &             ,I_START=grid%i_start, I_END=min(grid%i_end, ide-1)    &
     &             ,J_START=grid%j_start, J_END=min(grid%j_end, jde-1)    &
     &             ,KTS=k_start, KTE=min(k_end, kde-1)                    &
     &             ,NUM_TILES=grid%num_tiles                              &
     &             ,U=grid%u_phy, V=grid%v_phy, TH=th_phy, T=t_phy        &
     &             ,P=grid%p_hyd, PI=pi_phy, RHO=grid%rho, MOIST=moist    &
     &             ,NUM_MOIST=num_moist                                   &
     &             ,ITIMESTEP=grid%itimestep, DT=grid%dt, DX=grid%dx      &
     &             ,CUDT=grid%cudt                                        &
     &             ,CURR_SECS=curr_secs, ADAPT_STEP_FLAG=adapt_step_flag  &
     &             ,RAINSH=grid%rainsh, PRATESH=grid%pratesh, NCA=grid%nca&
     &             ,RAINSHV=grid%rainshv                                  &
     &             ,Z=grid%z, Z_AT_W=grid%z_at_w, DZ8W=dz8w               &
     &             ,MAVAIL=grid%mavail, PBLH=grid%pblh, P8W=grid%p_hyd_w  &
     &             ,TKE_PBL=grid%tke_pbl                                  &
     &             ,CLDFRA=grid%cldfra, CLDFRA_OLD=grid%cldfra_old        &
     &             ,CLDFRA_OLD_MP=grid%cldfra_old_mp                      &
     &             ,CLDFRA_CONV=grid%cldfra_conv                          &
     &             ,CLDFRASH=grid%cldfrash, HTOP=grid%htop, HBOT=grid%hbot&
     &             ,SHCU_PHYSICS=grid%shcu_physics                        &
     &             ,QV_CURR=moist(ims,kms,jms,P_QV)                       &
     &             ,QC_CURR=moist(ims,kms,jms,P_QC)                       &
     &             ,QR_CURR=moist(ims,kms,jms,P_QR)                       &
     &             ,QI_CURR=moist(ims,kms,jms,P_QI)                       &
     &             ,QS_CURR=moist(ims,kms,jms,P_QS)                       &
     &             ,QG_CURR=moist(ims,kms,jms,P_QG)                       &
     &             ,QNC_CURR=scalar(ims,kms,jms,P_QNC)                    & 
     &             ,QNI_CURR=scalar(ims,kms,jms,P_QNI)                    & 



     &             ,DLF=grid%dlf, RLIQ=grid%rliq, RLIQ2=grid%rliq2        &
     &             ,DLF2=grid%dlf2                                        & 
     &             ,CMFMC=grid%cmfmc, CMFMC2=grid%cmfmc2                  &
     &             ,CUSH=grid%cush, SNOWSH=grid%snowsh                    &
     &             ,ICWMRSH=grid%icwmrsh, RPRDSH=grid%rprdsh              &
     &             ,CBMF=grid%cbmf_cu, CMFSL=grid%cmfsl, CMFLQ=grid%cmflq &
     &             ,EVAPCSH=grid%evapcsh                                  &
     &             ,RQVSHTEN=grid%rqvshten, RQCSHTEN=grid%rqcshten        &
     &             ,RQRSHTEN=grid%rqrshten, RQISHTEN=grid%rqishten        &
     &             ,RQSSHTEN=grid%rqsshten, RQGSHTEN=grid%rqgshten        &
     &             ,RQCNSHTEN=grid%rqcnshten, RQINSHTEN=grid%rqinshten    &  
     &             ,RQVBLTEN=grid%rqvblten, RQVFTEN=grid%rqvften          &
     &             ,RUSHTEN=grid%rushten, RVSHTEN=grid%rvshten            &
     &             ,RTHSHTEN=grid%rthshten, RTHRATEN=grid%rthraten        &
     &             ,RTHBLTEN=grid%rthblten, RTHFTEN=grid%rthften          &
     &             ,F_QV=f_qv,F_QC=f_qc,F_QR=f_qr                         &
     &             ,F_QI=f_qi,F_QS=f_qs,F_QG=f_qg                         &
     &             ,HT=grid%ht                                            &
     &             ,SHFRC3D=grid%shfrc3d                                  & 
     &             ,IS_CAMMGMP_USED = grid%is_CAMMGMP_used                &
     &             ,WSTAR=grid%wstar_ysu,DELTA=grid%delta_ysu             &
     &             ,KPBL=grid%kpbl,ZNU=grid%znu                           &
     &             ,RAINCV=grid%raincv                                    &
     &                                                                    )




      CALL force_scm(itimestep=grid%itimestep,dt=grid%dt                  &
     &             ,scm_force=config_flags%scm_force                      &
     &             ,dx=config_flags%scm_force_dx                          &
     &             ,num_force_layers=grid%num_force_layers                &
     &             ,scm_th_adv=config_flags%scm_th_adv                    &
     &             ,scm_qv_adv=config_flags%scm_qv_adv                    &
     &             ,scm_ql_adv=config_flags%scm_ql_adv                    &
     &             ,scm_wind_adv=config_flags%scm_wind_adv                &
     &             ,scm_vert_adv=config_flags%scm_vert_adv                &
     &             ,scm_th_t_tend=config_flags%scm_th_t_tend              &
     &             ,scm_qv_t_tend=config_flags%scm_qv_t_tend              &
     &             ,scm_soilT_force=config_flags%scm_soilT_force          &
     &             ,scm_soilQ_force=config_flags%scm_soilQ_force          &
     &      ,scm_force_th_largescale=config_flags%scm_force_th_largescale &
     &      ,scm_force_qv_largescale=config_flags%scm_force_qv_largescale &
     &      ,scm_force_ql_largescale=config_flags%scm_force_ql_largescale &
     &      ,scm_force_wind_largescale=config_flags%scm_force_wind_largescale &
     &             ,u_base=grid%u_base,v_base=grid%v_base                 &
     &             ,z_base=grid%z_base                                    &
     &             ,z_force=grid%z_force,z_force_tend=grid%z_force_tend   &
     &             ,u_g=grid%u_g,v_g=grid%v_g                             &
     &             ,u_g_tend=grid%u_g_tend,v_g_tend=grid%v_g_tend         &
     &             ,w_subs=grid%w_subs, w_subs_tend=grid%w_subs_tend      &
     &             ,th_upstream_x=grid%th_upstream_x                      &
     &             ,th_upstream_x_tend=grid%th_upstream_x_tend            &
     &             ,th_upstream_y=grid%th_upstream_y                      &
     &             ,th_upstream_y_tend=grid%th_upstream_y_tend            &
     &             ,qv_upstream_x=grid%qv_upstream_x                      &
     &             ,qv_upstream_x_tend=grid%qv_upstream_x_tend            &
     &             ,qv_upstream_y=grid%qv_upstream_y                      &
     &             ,qv_upstream_y_tend=grid%qv_upstream_y_tend            &
     &             ,ql_upstream_x=grid%ql_upstream_x                      &
     &             ,ql_upstream_x_tend=grid%ql_upstream_x_tend            &
     &             ,ql_upstream_y=grid%ql_upstream_y                      &
     &             ,ql_upstream_y_tend=grid%ql_upstream_y_tend            &
     &             ,u_upstream_x=grid%u_upstream_x                        &
     &             ,u_upstream_x_tend=grid%u_upstream_x_tend              &
     &             ,u_upstream_y=grid%u_upstream_y                        &
     &             ,u_upstream_y_tend=grid%u_upstream_y_tend              &
     &             ,v_upstream_x=grid%v_upstream_x                        &
     &             ,v_upstream_x_tend=grid%v_upstream_x_tend              &
     &             ,v_upstream_y=grid%v_upstream_y                        &
     &             ,v_upstream_y_tend=grid%v_upstream_y_tend              &
     &             ,th_t_tend=grid%th_t_tend                              &
     &             ,qv_t_tend=grid%qv_t_tend                              &
     &             ,tau_x=grid%tau_x                                      &
     &             ,tau_x_tend=grid%tau_x_tend                            &
     &             ,tau_y=grid%tau_y                                      &
     &             ,tau_y_tend=grid%tau_y_tend                            &
     &             ,th_largescale=grid%th_largescale                      &
     &             ,th_largescale_tend=grid%th_largescale_tend            &
     &             ,qv_largescale=grid%qv_largescale                      &
     &             ,qv_largescale_tend=grid%qv_largescale_tend            &
     &             ,ql_largescale=grid%ql_largescale                      &
     &             ,ql_largescale_tend=grid%ql_largescale_tend            &
     &             ,u_largescale=grid%u_largescale                      &
     &             ,u_largescale_tend=grid%u_largescale_tend            &
     &             ,v_largescale=grid%v_largescale                      &
     &             ,v_largescale_tend=grid%v_largescale_tend            &
     &             ,tau_largescale=grid%tau_largescale                    &
     &             ,tau_largescale_tend=grid%tau_largescale_tend          &
     &             ,num_force_soil_layers=config_flags%num_force_soil_layers &
     &             ,num_soil_layers=config_flags%num_soil_layers          &
     &             ,soil_depth_force=grid%soil_depth_force                &
     &             ,zs=grid%zs                                            &
     &             ,tslb=grid%tslb,smois=grid%smois                       &
     &             ,t_soil_forcing_val=grid%t_soil_forcing_val            &
     &             ,t_soil_forcing_tend=grid%t_soil_forcing_tend          &
     &             ,q_soil_forcing_val=grid%q_soil_forcing_val            &
     &             ,q_soil_forcing_tend=grid%q_soil_forcing_tend          &
     &             ,tau_soil=grid%tau_soil                                &
     &             ,z=grid%z,z_at_w=grid%z_at_w                           &
     &             ,th=th_phy, qv=moist(ims,kms,jms,P_QV)                 &
     &             ,ql=moist(ims,kms,jms,P_QC)                            &
     &             ,u=grid%u_phy, v=grid%v_phy                                      &
     &             ,thten=grid%rthblten, qvten=grid%rqvblten              &
     &             ,qlten=grid%rqcblten                                   &
     &             ,uten=grid%rublten, vten=grid%rvblten                  &
     &             ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde     &
     &             ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme     &
     &             ,IPS=ips,IPE=ipe, JPS=jps,JPE=jpe, KPS=kps,KPE=kpe     &
     &             ,KTS=k_start, KTE=min(k_end,kde-1)                     &
     &              )







CALL HALO_EM_FDDA_SFC_sub ( grid, &
  local_communicator, &
  mytask, ntasks, ntasks_x, ntasks_y, &
  ids, ide, jds, jde, kds, kde,       &
  ims, ime, jms, jme, kms, kme,       &
  ips, ipe, jps, jpe, kps, kpe )






      CALL wrf_debug ( 200 , ' call fddagd_driver' )


      CALL fddagd_driver(itimestep=grid%itimestep,dt=grid%dt,xtime=grid%XTIME,      &
                  id=grid%id,                                                       &
                  RUNDGDTEN=grid%rundgdten,RVNDGDTEN=grid%rvndgdten,                &
                  RTHNDGDTEN=grid%rthndgdten,RPHNDGDTEN=grid%rphndgdten,            &
                  RQVNDGDTEN=grid%rqvndgdten,RMUNDGDTEN=grid%rmundgdten,            &
                  u_ndg_old=fdda3d(ims,kms,jms,P_u_ndg_old),                        &
                  v_ndg_old=fdda3d(ims,kms,jms,P_v_ndg_old),                        &
                  t_ndg_old=fdda3d(ims,kms,jms,P_t_ndg_old),                        &
                  ph_ndg_old=fdda3d(ims,kms,jms,P_ph_ndg_old),                      &
                  q_ndg_old=fdda3d(ims,kms,jms,P_q_ndg_old),                        &
                  mu_ndg_old=fdda2d(ims,1,jms,P_mu_ndg_old),                        &
                  u_ndg_new=fdda3d(ims,kms,jms,P_u_ndg_new),                        &
                  v_ndg_new=fdda3d(ims,kms,jms,P_v_ndg_new),                        &
                  t_ndg_new=fdda3d(ims,kms,jms,P_t_ndg_new),                        &
                  ph_ndg_new=fdda3d(ims,kms,jms,P_ph_ndg_new),                      &
                  q_ndg_new=fdda3d(ims,kms,jms,P_q_ndg_new),                        &
                  mu_ndg_new=fdda2d(ims,1,jms,P_mu_ndg_new),                        &
                  u3d=grid%u_2,v3d=grid%v_2,th_phy=th_phy,                          &
                  ph=grid%ph_2,rho=grid%rho,moist=moist,                            &
                  p_phy=p_phy,pi_phy=pi_phy,p8w=p8w,t_phy=grid%t_phy,               &
                  dz8w=dz8w,z=grid%z,z_at_w=grid%z_at_w,                            &
                  grid=grid,config_flags=config_flags,dx=grid%DX,n_moist=num_moist, &
                  STEPFG=grid%STEPFG,                                               &
                  pblh=grid%pblh,ht=grid%ht,REGIME=grid%regime,ZNT=grid%znt         &
                   ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde               &
                   ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme               &
                   ,I_START=grid%i_start,I_END=min(grid%i_end, ide-1)               &
                   ,J_START=grid%j_start,J_END=min(grid%j_end, jde-1)               &
                   ,KTS=k_start, KTE=min(k_end,kde-1)                               &
                   , num_tiles=grid%num_tiles,                                      &
                   u10=grid%u10, v10=grid%v10, th2=grid%th2, q2=grid%q2,            &
                   u10_ndg_old=grid%u10_ndg_old,                    &
                   v10_ndg_old=grid%v10_ndg_old,                    &
                   t2_ndg_old=grid%t2_ndg_old,                      &
                   th2_ndg_old=grid%th2_ndg_old,                    &
                   q2_ndg_old=grid%q2_ndg_old,                      &
                   rh_ndg_old=grid%rh_ndg_old,                      &
                   psl_ndg_old=grid%psl_ndg_old,                    &
                   ps_ndg_old=grid%ps_ndg_old,                      &
                   tob_ndg_old=grid%tob_ndg_old,                                    &
                   odis_ndg_old=grid%odis_ndg_old,                                  &
                   u10_ndg_new=grid%u10_ndg_new,                    &
                   v10_ndg_new=grid%v10_ndg_new,                    &
                   t2_ndg_new=grid%t2_ndg_new,                      &
                   th2_ndg_new=grid%th2_ndg_new,                    &
                   q2_ndg_new=grid%q2_ndg_new,                      &
                   rh_ndg_new=grid%rh_ndg_new,                      &
                   psl_ndg_new=grid%psl_ndg_new,                    &
                   ps_ndg_new=grid%ps_ndg_new,                      &
                   tob_ndg_new=grid%tob_ndg_new,                                    &
                   odis_ndg_new=grid%odis_ndg_new                                   &
                   ,IPS=ips,IPE=ipe, JPS=jps,JPE=jpe, KPS=kps,KPE=kpe               &
                   ,IMSX=imsx,IMEX=imex,JMSX=jmsx,JMEX=jmex,KMSX=kmsx,KMEX=kmex     &
                   ,IPSX=ipsx,IPEX=ipex,JPSX=jpsx,JPEX=jpex,KPSX=kpsx,KPEX=kpex     &
                   ,IMSY=imsy,IMEY=imey,JMSY=jmsy,JMEY=jmey,KMSY=kmsy,KMEY=kmey     &
                   ,IPSY=ipsy,IPEY=ipey,JPSY=jpsy,JPEY=jpey,KPSY=kpsy,KPEY=kpey     )



  END SUBROUTINE first_rk_step_part1

END MODULE module_first_rk_step_part1