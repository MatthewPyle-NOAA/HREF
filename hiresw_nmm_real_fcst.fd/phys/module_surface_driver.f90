!WRF:MEDIATION_LAYER:PHYSICS
!
MODULE module_surface_driver
CONTAINS

   SUBROUTINE surface_driver(                                         &
     &           acsnom,acsnow,akhs,akms,albedo,br,canwat             &
     &          ,chklowq,dt,dx,dz8w,dzs,glw                           &
     &          ,grdflx,gsw,swdown,gz1oz0,hfx,ht,ifsnow,isfflx        &
     &          ,isltyp,itimestep,ivgtyp,lowlyr,mavail,rmol           &
     &          ,num_soil_layers,p8w,pblh,pi_phy,pshltr,psih          &
     &          ,psim,p_phy,q10,q2,qfx,qsfc,qshltr,qz0                &
     &          ,raincv,rho,sfcevp,sfcexc,sfcrunoff                   &
     &          ,smois,smstav,smstot,snoalb,snow,snowc,snowh,stepbl   &
     &          ,th10,th2,thz0,th_phy,tmn,tshltr,tsk,tslb             &
     &          ,t_phy,u10,udrunoff,ust,uz0,u_frame,u_phy,v10,vegfra  &
     &          ,vz0,v_frame,v_phy,warm_rain,wspd,xice,xland,z,znt,zs &
     &          ,ct,tke_myj                                           &
     &          ,albbck,lh,sh2o,shdmax,shdmin,z0                      &
     &          ,flqc,flhc,psfc,sst,sst_update,t2,emiss                        &
     &          ,sf_sfclay_physics,sf_surface_physics,ra_lw_physics   &
            !  Optional urban 
     &          ,declin_urb,cosz_urb2d,omg_urb2d,xlat_urb2d           & !I urban
     &          ,num_roof_layers, num_wall_layers                     & !I urban
     &          ,num_road_layers, dzr, dzb, dzg                       & !I urban
     &          ,tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d         & !H urban
     &          ,uc_urb2d                                             & !H urban
     &          ,xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d          & !H urban
     &          ,trl_urb3d,tbl_urb3d,tgl_urb3d                        & !H urban
     &          ,sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d,ts_urb2d          & !H urban
     &          ,frc_urb2d, utype_urb2d                               & !H urban
     &          ,ucmcall                                              & ! urban
     &          , ids,ide,jds,jde,kds,kde                             &
     &          , ims,ime,jms,jme,kms,kme                             &
     &          , i_start,i_end,j_start,j_end,kts,kte,num_tiles       &
             !  Optional moisture tracers
     &           ,qv_curr, qc_curr, qr_curr                           &
     &           ,qi_curr, qs_curr, qg_curr                           &
             !  Optional moisture tracer flags
     &           ,f_qv,f_qc,f_qr                                      &
     &           ,f_qi,f_qs,f_qg                                      &
             !  Other optionals (more or less em specific)
     &          ,capg,hol,mol                                   &
     &          ,rainncv,rainbl,regime,thc                         &
     &          ,qsg,qvg,qcg,soilt1,tsnav                             &
     &          ,smfr3d,keepfr3dflag                                  &
             !  Other optionals (more or less nmm specific)
     &          ,potevp,snopcx,soiltb,sr                              &
             !  Optional observation nudging
     &          ,uratx,vratx,tratx                                    &
                                                                      )

   USE module_state_description, ONLY : SFCLAYSCHEME              &
                                       ,MYJSFCSCHEME              &
                                       ,GFSSFCSCHEME              &
                                       ,SLABSCHEME                &
                                       ,NMMLSMSCHEME              &
                                       ,LSMSCHEME                 &
                                       ,RUCLSMSCHEME
   USE module_model_constants
! *** add new modules of schemes here

   USE module_sf_sfclay
   USE module_sf_myjsfc
   USE module_sf_gfs
   USE module_sf_noahlsm
   USE module_sf_ruclsm
   USE module_sf_lsm_nmm

   USE module_sf_slab
!
   USE module_sf_sfcdiags
!

   !  This driver calls subroutines for the surface parameterizations.
   !
   !  surface layer: (between surface and pbl)
   !      1. sfclay
   !      2. myjsfc
   !  surface: ground temp/lsm scheme:
   !      1. slab
   !      2. Noah LSM
   !      99. NMM LSM (NMM core only)
!------------------------------------------------------------------
   IMPLICIT NONE
!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! Definitions
!-----------
! Theta      potential temperature (K)
! Qv         water vapor mixing ratio (kg/kg)
! Qc         cloud water mixing ratio (kg/kg)
! Qr         rain water mixing ratio (kg/kg)
! Qi         cloud ice mixing ratio (kg/kg)
! Qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- itimestep     number of time steps
!-- GLW           downward long wave flux at ground surface (W/m^2)
!-- GSW           net short wave flux at ground surface (W/m^2)
!-- SWDOWN        downward short wave flux at ground surface (W/m^2)
!-- EMISS         surface emissivity (between 0 and 1)
!-- TSK           surface temperature (K)
!-- TMN           soil temperature at lower boundary (K)
!-- XLAND         land mask (1 for land, 2 for water)
!-- ZNT           time-varying roughness length (m)
!-- Z0            background roughness length (m)
!-- MAVAIL        surface moisture availability (between 0 and 1)
!-- UST           u* in similarity theory (m/s)
!-- MOL           T* (similarity theory) (K)
!-- HOL           PBL height over Monin-Obukhov length
!-- PBLH          PBL height (m)
!-- CAPG          heat capacity for soil (J/K/m^3)
!-- THC           thermal inertia (Cal/cm/K/s^0.5)
!-- SNOWC         flag indicating snow coverage (1 for snow cover)
!-- HFX           net upward heat flux at the surface (W/m^2)
!-- QFX           net upward moisture flux at the surface (kg/m^2/s)
!-- LH            net upward latent heat flux at surface (W/m^2)
!-- REGIME        flag indicating PBL regime (stable, unstable, etc.)
!-- tke_myj       turbulence kinetic energy from Mellor-Yamada-Janjic (MYJ) (m^2/s^2)
!-- akhs          sfc exchange coefficient of heat/moisture from MYJ
!-- akms          sfc exchange coefficient of momentum from MYJ
!-- thz0          potential temperature at roughness length (K)
!-- uz0           u wind component at roughness length (m/s)
!-- vz0           v wind component at roughness length (m/s)
!-- qsfc          specific humidity at lower boundary (kg/kg)
!-- uratx         ratio of u over u10 (Added for obs-nudging)
!-- vratx         ratio of v over v10 (Added for obs-nudging)
!-- tratx         ratio of t over th2 (Added for obs-nudging)
!-- u10           diagnostic 10-m u component from surface layer
!-- v10           diagnostic 10-m v component from surface layer
!-- th2           diagnostic 2-m theta from surface layer and lsm
!-- t2            diagnostic 2-m temperature from surface layer and lsm
!-- q2            diagnostic 2-m mixing ratio from surface layer and lsm
!-- tshltr        diagnostic 2-m theta from MYJ
!-- th10          diagnostic 10-m theta from MYJ
!-- qshltr        diagnostic 2-m specific humidity from MYJ
!-- q10           diagnostic 10-m specific humidity from MYJ
!-- lowlyr        index of lowest model layer above ground
!-- rr            dry air density (kg/m^3)
!-- u_phy         u-velocity interpolated to theta points (m/s)
!-- v_phy         v-velocity interpolated to theta points (m/s)
!-- th_phy        potential temperature (K)
!-- moist         moisture array (4D - last index is species) (kg/kg)
!-- p_phy         pressure (Pa)
!-- pi_phy        exner function (dimensionless)
!-- pshltr        diagnostic shelter (2m) pressure from MYJ (Pa)
!-- p8w           pressure at full levels (Pa)
!-- t_phy         temperature (K)
!-- dz8w          dz between full levels (m)
!-- z             height above sea level (m)
!-- DX            horizontal space interval (m)
!-- DT            time step (second)
!-- PSFC          pressure at the surface (Pa)
!-- SST           sea-surface temperature (K)
!-- TSLB          
!-- ZS
!-- DZS
!-- num_soil_layers number of soil layer
!-- IFSNOW      ifsnow=1 for snow-cover effects
!
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!
!******************************************************************
!------------------------------------------------------------------ 

   INTEGER, INTENT(IN) ::                                             &
     &           ids,ide,jds,jde,kds,kde                              &
     &          ,ims,ime,jms,jme,kms,kme                              &
     &          ,kts,kte,num_tiles

   INTEGER, INTENT(IN) :: sf_sfclay_physics,sf_surface_physics,ra_lw_physics,sst_update

   INTEGER, INTENT(IN) :: ucmcall                                     !urban

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                       &
     &           i_start,i_end,j_start,j_end

   INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   ISLTYP
   INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   IVGTYP
   INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   LOWLYR
   INTEGER, INTENT(IN )::   IFSNOW
   INTEGER, INTENT(IN )::   ISFFLX
   INTEGER, INTENT(IN )::   ITIMESTEP
   INTEGER, INTENT(IN )::   NUM_SOIL_LAYERS
   INTEGER, INTENT(IN )::   STEPBL
   LOGICAL, INTENT(IN )::   WARM_RAIN
   REAL , INTENT(IN )::   U_FRAME
   REAL , INTENT(IN )::   V_FRAME
   REAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), INTENT(INOUT)::   SMOIS
   REAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), INTENT(INOUT)::   TSLB
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   GLW
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   GSW,SWDOWN
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   HT
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   RAINCV
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   SST
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   TMN
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   VEGFRA
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   XICE
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   XLAND
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT)::   MAVAIL
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT)::   SNOALB
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   ACSNOW
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   AKHS
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   AKMS
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   ALBEDO
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   CANWAT


   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   GRDFLX
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   HFX
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   RMOL
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   PBLH
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   Q2
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   QFX
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   QSFC
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   QZ0
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SFCRUNOFF
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SMSTAV
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SMSTOT
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SNOW
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SNOWC
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SNOWH
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   TH2
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   THZ0
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   TSK
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   UDRUNOFF
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   UST
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   UZ0
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   VZ0
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   WSPD
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   ZNT
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   BR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   CHKLOWQ
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   GZ1OZ0
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   PSHLTR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   PSIH
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   PSIM
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   Q10
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   QSHLTR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   TH10
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   TSHLTR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   U10
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   V10
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   PSFC
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)::   ACSNOM
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)::   SFCEVP
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)::   SFCEXC
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)::   FLHC
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)::   FLQC
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::   CT
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   DZ8W
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   P8W
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   PI_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   P_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   RHO
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   TH_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   T_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   U_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   V_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   Z
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) ::   TKE_MYJ
   REAL, DIMENSION(1:num_soil_layers), INTENT(IN)::   DZS
   REAL, DIMENSION(1:num_soil_layers), INTENT(IN)::   ZS
   REAL, INTENT(IN )::   DT
   REAL, INTENT(IN )::   DX

!  arguments for NCAR surface physics

   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT )::   ALBBCK  ! INOUT needed for NMM
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT )::   LH
   REAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), INTENT(INOUT)::   SH2O
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   SHDMAX
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   SHDMIN
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   Z0

!
! Optional
!

!
!  Observation nudging
!
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(OUT)::   uratx  !Added for obs-nudging
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(OUT)::   vratx  !Added for obs-nudging
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(OUT)::   tratx  !Added for obs-nudging
!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL, INTENT(IN), OPTIONAL ::                             &
                                                      f_qv      &
                                                     ,f_qc      &
                                                     ,f_qr      &
                                                     ,f_qi      &
                                                     ,f_qs      &
                                                     ,f_qg

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         OPTIONAL, INTENT(INOUT) ::                              &
                      ! optional moisture tracers
                      ! 2 time levels; if only one then use CURR
                      qv_curr, qc_curr, qr_curr                  &
                     ,qi_curr, qs_curr, qg_curr

   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   capg
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   emiss
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   hol
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   mol
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   regime
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(IN )::     rainncv
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   RAINBL
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   t2
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(IN )::     thc
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   qsg
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   qvg
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   qcg
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   soilt1
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   tsnav
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   potevp ! NMM LSM
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   snopcx ! NMM LSM
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   soiltb ! NMM LSM
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   sr ! NMM and RUC LSM
   REAL, DIMENSION( ims:ime, 1:num_soil_layers, jms:jme ), OPTIONAL, INTENT(INOUT)::   smfr3d
   REAL, DIMENSION( ims:ime, 1:num_soil_layers, jms:jme ), OPTIONAL, INTENT(INOUT)::   keepfr3dflag

!  LOCAL  VAR

   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ) ::v_phytmp
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ) ::u_phytmp

   REAL,       DIMENSION( ims:ime, jms:jme )          ::  ZOL

   REAL,       DIMENSION( ims:ime, jms:jme )          ::          &
                                                             QGH, &
                                                             CHS, &
                                                             CPM, &
                                                            CHS2, &
                                                            CQS2

   REAL    :: DTMIN,DTBL
!
   INTEGER :: i,J,K,NK,jj,ij
   LOGICAL :: radiation, myj, frpcpn
!-------------------------------------------------
! urban related variables are added to declaration
!-------------------------------------------------
     REAL, OPTIONAL, INTENT(IN) :: DECLIN_URB                                 !urban
     REAL, OPTIONAL , DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: COSZ_URB2D  !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: OMG_URB2D   !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: XLAT_URB2D  !urban
     INTEGER,  OPTIONAL, INTENT(IN) :: num_roof_layers                         !urban
     INTEGER,  OPTIONAL, INTENT(IN) :: num_wall_layers                         !urban
     INTEGER,  OPTIONAL, INTENT(IN) :: num_road_layers                         !urban
     REAL, OPTIONAL, DIMENSION(1:num_soil_layers), INTENT(IN) :: DZR          !urban
     REAL, OPTIONAL, DIMENSION(1:num_soil_layers), INTENT(IN) :: DZB          !urban
     REAL, OPTIONAL, DIMENSION(1:num_soil_layers), INTENT(IN) :: DZG          !urban

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: TR_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: TB_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: TG_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: TC_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: QC_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: UC_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: XXXR_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: XXXB_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: XXXG_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: XXXC_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), &       !urban
           INTENT(INOUT)  :: TRL_URB3D                                 !urban
     REAL, OPTIONAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), &       !urban
           INTENT(INOUT)  :: TBL_URB3D                                 !urban
     REAL, OPTIONAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), &       !urban
           INTENT(INOUT)  :: TGL_URB3D                                 !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: SH_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: LH_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: G_URB2D  !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: RN_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: TS_URB2D !urban
!
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: FRC_URB2D  !urban 
     INTEGER, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: UTYPE_URB2D  !urban 

     REAL,  DIMENSION( ims:ime, jms:jme )  :: PSIM_URB2D  !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: PSIH_URB2D  !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: GZ1OZ0_URB2D  !urban local var
!m     REAL, DIMENSION( ims:ime, jms:jme ) :: AKHS_URB2D  !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: AKMS_URB2D  !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: U10_URB2D   !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: V10_URB2D   !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: TH2_URB2D   !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: Q2_URB2D    !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: UST_URB2D  !urban local var

!------------------------------------------------------------------
   CHARACTER*256 :: message
!------------------------------------------------------------------
!

  if (sf_sfclay_physics .eq. 0) return

  v_phytmp = 0.
  u_phytmp = 0.
  ZOL = 0.
  QGH = 0.
  CHS = 0.
  CPM = 0.
  CHS2 = 0.
  DTMIN = 0.
  DTBL = 0.

! RAINBL in mm (Accumulation between PBL calls)

  IF ( PRESENT( rainncv ) .AND. PRESENT( rainbl ) ) THEN
    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij, i, j, k )
    DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         RAINBL(i,j) = RAINBL(i,j) + RAINCV(i,j) + RAINNCV(i,j) 
         RAINBL(i,j) = MAX (RAINBL(i,j), 0.0) 
      ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ELSE IF ( PRESENT( rainbl ) ) THEN
    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij, i, j, k )
    DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         RAINBL(i,j) = RAINBL(i,j) + RAINCV(i,j)
         RAINBL(i,j) = MAX (RAINBL(i,j), 0.0)
      ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ENDIF
! Update SST
  IF (sst_update .EQ. 1) THEN
    !$OMP PARALLEL DO   &
    !$OMP PRIVATE ( ij, i, j, k )
    DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
        IF(XLAND(i,j) .GT. 1.5)TSK(i,j)=SST(i,j)
      ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ENDIF

  IF (itimestep .eq. 1 .or. mod(itimestep,STEPBL) .eq. 0) THEN

  radiation = .false.
  myj = .false.
  frpcpn = .false.

  IF (ra_lw_physics .gt. 0) radiation = .true.

!---- 
! CALCULATE CONSTANT
 
     DTMIN=DT/60.
! Surface schemes need PBL time step for updates and accumulations
! Assume these schemes provide no tendencies
     DTBL=DT*STEPBL

! SAVE OLD VALUES


     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij, i, j, k )
     DO ij = 1 , num_tiles
       DO j=j_start(ij),j_end(ij)
       DO i=i_start(ij),i_end(ij)
! PSFC : in Pa
          PSFC(I,J)=p8w(I,kts,J)
! REVERSE ORDER IN THE VERTICAL DIRECTION
          DO k=kts,kte
            v_phytmp(i,k,j)=v_phy(i,k,j)+v_frame
            u_phytmp(i,k,j)=u_phy(i,k,j)+u_frame
          ENDDO
       ENDDO
       ENDDO
     ENDDO
     !$OMP END PARALLEL DO

     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij, i, j, k )
     DO ij = 1 , num_tiles
     sfclay_select: SELECT CASE(sf_sfclay_physics)

     CASE (SFCLAYSCHEME)
       CALL wrf_error_fatal3 ( "module_surface_driver.b" , 573 , 'SFCLAY cannot be used with NMM')
      CASE (MYJSFCSCHEME)
       IF (PRESENT(qv_curr)    .AND.  PRESENT(qc_curr) .AND.    &
                                                      .TRUE. ) THEN

        myj =.true.

            CALL wrf_debug(100,'in MYJSFC')
            CALL MYJSFC(itimestep,ht,dz8w,                         &
              p_phy,p8w,th_phy,t_phy,                              &
              qv_curr,qc_curr,                                      &
              u_phy,v_phy,tke_myj,                                 &
              tsk,qsfc,thz0,qz0,uz0,vz0,                           &
              lowlyr,                                              &
              xland,                                               &
              ust,znt,z0,pblh,mavail,rmol,                         &
              akhs,akms,                                           &
              chs,chs2,cqs2,hfx,qfx,lh,flhc,flqc,qgh,cpm,ct,       &
              u10,v10,t2,th2,tshltr,th10,q2,qshltr,q10,pshltr,               &
              ids,ide, jds,jde, kds,kde,                           &
              ims,ime, jms,jme, kms,kme,                           &
              i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )
       ELSE
         CALL wrf_error_fatal3 ( "module_surface_driver.b" , 597 , 'Lacking arguments for MYJSFC in surface driver')
       ENDIF

     CASE (GFSSFCSCHEME)
       IF (PRESENT(qv_curr) .AND. .TRUE. ) THEN
       CALL wrf_debug( 100, 'in GFSSFC' )
         CALL SF_GFS(u_phytmp,v_phytmp,t_phy,qv_curr,              &
               p_phy,CP,RCP,R_d,XLV,PSFC,CHS,CHS2,CQS2,CPM,        &
               ZNT,UST,PSIM,PSIH,                                  &
               XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,                     &
               QGH,QSFC,U10,V10,                                   &
               GZ1OZ0,WSPD,BR,ISFFLX,                              &
               EP_1,EP_2,KARMAN,itimestep,                         &
               ids,ide, jds,jde, kds,kde,                          &
               ims,ime, jms,jme, kms,kme,                          &
               i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )
        CALL wrf_debug(100,'in SFCDIAGS')
       ELSE
         CALL wrf_error_fatal3 ( "module_surface_driver.b" , 615 , 'Lacking arguments for SF_GFS in surface driver')
       ENDIF

     CASE DEFAULT
        
       WRITE( message , * )                                &
   'The sfclay option does not exist: sf_sfclay_physics = ', sf_sfclay_physics
       CALL wrf_error_fatal3 ( "module_surface_driver.b" , 622 ,  message )

     END SELECT sfclay_select
     ENDDO
     !$OMP END PARALLEL DO

     IF (ISFFLX.EQ.0 ) GOTO 430
     !$OMP PARALLEL DO   &
     !$OMP PRIVATE ( ij, i, j, k )
     DO ij = 1 , num_tiles

     sfc_select: SELECT CASE(sf_surface_physics)

     CASE (SLABSCHEME)

       IF (PRESENT(qv_curr)                            .AND.    &
           PRESENT(capg)        .AND.    &
                                                      .TRUE. ) THEN
           DO j=j_start(ij),j_end(ij)
           DO i=i_start(ij),i_end(ij)
!          CQS2 ACCOUNTS FOR MAVAIL FOR SFCDIAGS 2M Q
              CQS2(I,J)= CQS2(I,J)*MAVAIL(I,J)
           ENDDO
           ENDDO

        CALL wrf_debug(100,'in SLAB')
          CALL SLAB(t_phy,qv_curr,p_phy,flhc,flqc,  &
             psfc,xland,tmn,hfx,qfx,lh,tsk,qsfc,chklowq,          &
             gsw,glw,capg,thc,snowc,emiss,mavail,                 &
             dtbl,rcp,xlv,dtmin,ifsnow,                           &
             svp1,svp2,svp3,svpt0,ep_2,karman,eomeg,stbolt,       &
             tslb,zs,dzs,num_soil_layers,radiation,               &
             ids,ide, jds,jde, kds,kde,                           &
             ims,ime, jms,jme, kms,kme,                           &
             i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )

           DO j=j_start(ij),j_end(ij)
           DO i=i_start(ij),i_end(ij)
              SFCEVP(I,J)= SFCEVP(I,J) + QFX(I,J)*DTBL
           ENDDO
           ENDDO

        CALL wrf_debug(100,'in SFCDIAGS')
          CALL SFCDIAGS(hfx,qfx,tsk,qsfc,chs2,cqs2,t2,th2,q2,      &
                     psfc,cp,r_d,rcp,                              &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
             i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )

       ELSE
         CALL wrf_error_fatal3 ( "module_surface_driver.b" , 672 , 'Lacking arguments for SLAB in surface driver')
       ENDIF

     CASE (NMMLSMSCHEME)
       IF (PRESENT(qv_curr)    .AND.  PRESENT(rainbl)  .AND.    &
           PRESENT(potevp)     .AND.  PRESENT(snopcx)  .AND.    &
           PRESENT(soiltb)     .AND.  PRESENT(sr)      .AND.    &
                                                      .TRUE. ) THEN
           CALL wrf_debug(100,'in NMM LSM')
           CALL nmmlsm(dz8w,qv_curr,p8w,rho,                    &
                t_phy,th_phy,tsk,chs,                           &
                hfx,qfx,qgh,swdown,glw,lh,rmol,                 &
                smstav,smstot,sfcrunoff,                        &
                udrunoff,ivgtyp,isltyp,vegfra,sfcevp,potevp,    &
                grdflx,sfcexc,acsnow,acsnom,snopcx,             &
                albbck,tmn,xland,xice,qz0,                      &
                th2,q2,snowc,cqs2,qsfc,soiltb,chklowq,rainbl,   &
                num_soil_layers,dtbl,dzs,itimestep,             &
                smois,tslb,snow,canwat,cpm,rcp,sr,              &    !tslb
                albedo,snoalb,sh2o,snowh,                       &
                ids,ide, jds,jde, kds,kde,                      &
                ims,ime, jms,jme, kms,kme,                      &
                i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte )
          CALL wrf_debug(100,'back from NMM LSM')
       ELSE
         CALL wrf_error_fatal3 ( "module_surface_driver.b" , 698 , 'Lacking arguments for NMMLSM in surface driver')
       ENDIF

     CASE (LSMSCHEME)

       IF (PRESENT(qv_curr)    .AND.  PRESENT(rainbl)        .AND.    &
!          PRESENT(emiss)      .AND.  PRESENT(t2)            .AND.    &
!          PRESENT(declin_urb) .AND.  PRESENT(cosz_urb2d)    .AND.    &
!          PRESENT(omg_urb2d)  .AND. PRESENT( xlat_urb2d)    .AND.    &       
!          PRESENT(dzr)       .AND.    & 
!          PRESENT( dzb)            .AND. PRESENT(dzg)       .AND.    &
!          PRESENT(tr_urb2d) .AND. PRESENT(tb_urb2d)         .AND.    &
!          PRESENT(tg_urb2d) .AND. PRESENT(tc_urb2d) .AND.            &
!          PRESENT(qc_urb2d) .AND. PRESENT(uc_urb2d) .AND.            & 
!          PRESENT(xxxr_urb2d) .AND. PRESENT(xxxb_urb2d) .AND.        & 
!          PRESENT(xxxg_urb2d) .AND.                                  &
!          PRESENT(xxxc_urb2d) .AND. PRESENT(trl_urb3d) .AND.         &
!          PRESENT(tbl_urb3d)   .AND. PRESENT(tgl_urb3d)  .AND.       &         
!          PRESENT(sh_urb2d) .AND. PRESENT(lh_urb2d)  .AND.           &
!          PRESENT(g_urb2d)   .AND. PRESENT(rn_urb2d) .AND.           &
!          PRESENT(ts_urb2d)                          .AND.           & 
!          PRESENT(frc_urb2d) .AND. PRESENT(utype_urb2d)   .AND.      &          
                                                      .TRUE. ) THEN
!------------------------------------------------------------------
         CALL wrf_debug(100,'in NOAH LSM')
           CALL lsm(dz8w,qv_curr,p8w,t_phy,tsk,                 &
                hfx,qfx,lh,grdflx,qgh,gsw,glw,smstav,smstot,    &
                sfcrunoff,udrunoff,ivgtyp,isltyp,vegfra,        &
                albedo,albbck,znt,z0, tmn,xland,xice, emiss,    &
                snowc,qsfc,rainbl,                              & 
                num_soil_layers,dtbl,dzs,itimestep,             &
                smois,tslb,snow,canwat,                         &
                chs, chs2, cqs2, cpm,rcp,                       &    
                sh2o,snowh,                                     & !h  
                u_phy,v_phy,                                    & !I
                snoalb,shdmin,shdmax,                           & !i
                acsnom,acsnow,                                  & !o 
                ids,ide, jds,jde, kds,kde,                      &
                ims,ime, jms,jme, kms,kme,                      &
                i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte,    &
                ucmcall                                         &
!Optional urban
                ,tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,  & !H urban
                uc_urb2d,                                       & !H urban
                xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d,    & !H urban
                trl_urb3d,tbl_urb3d,tgl_urb3d,                  & !H urban
                sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d,ts_urb2d,    & !H urban
                psim_urb2d,psih_urb2d,u10_urb2d,v10_urb2d,      & !O urban
                GZ1OZ0_urb2d, AKMS_URB2D,                       & !O urban
                th2_urb2d,q2_urb2d,ust_urb2d,                   & !O urban
                declin_urb,cosz_urb2d,omg_urb2d,                & !I urban
                xlat_urb2d,                                     & !I urban
                num_roof_layers, num_wall_layers,               & !I urban
                num_road_layers, DZR, DZB, DZG,                 & !I urban
                FRC_URB2D, UTYPE_URB2D                          & ! urban
                )


           DO j=j_start(ij),j_end(ij)
           DO i=i_start(ij),i_end(ij)
              CHKLOWQ(I,J)= 1.0
              SFCEVP(I,J)= SFCEVP(I,J) + QFX(I,J)*DTBL
           ENDDO
           ENDDO
         
          CALL SFCDIAGS(HFX,QFX,TSK,QSFC,CHS2,CQS2,T2,TH2,Q2,      &
                     PSFC,CP,R_d,RCP,                              &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
             i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )

!urban
     IF(UCMCALL.eq.1) THEN
       DO j=j_start(ij),j_end(ij)                             !urban
         DO i=i_start(ij),i_end(ij)                           !urban
          IF( IVGTYP(I,J) == 1 .or. IVGTYP(I,J) == 31 .or. &  !urban
              IVGTYP(I,J) == 32 .or. IVGTYP(I,J) == 33 ) THEN !urban
!             TH2(I,J)  = TH2_URB2D(I,J)                       !urban
!             T2(I,J)   = TH2_URB2D(I,J)/(1.E5/PSFC(I,J))**RCP !urban
!m             T2(I,J)   = TH2_URB2D(I,J)                       !urban
             T2(I,J)   = FRC_URB2D(i,j)*TH2_URB2D(I,J) + (1-FRC_URB2D(i,j))*T2(I,J) !urban
             TH2(I,J) = T2(I,J)*(1.E5/PSFC(I,J))**RCP                               !urban
!m             Q2(I,J)   = Q2_URB2D(I,J)                                            !urban
             Q2(I,J)   = FRC_URB2D(i,j)*Q2_URB2D(I,J) +(1-FRC_URB2D(i,j))* Q2(I,J)  !urban
             U10(I,J)  = U10_URB2D(I,J)                       !urban
             V10(I,J)  = V10_URB2D(I,J)                       !urban
             PSIM(I,J) = PSIM_URB2D(I,J)                      !urban
             PSIH(I,J) = PSIH_URB2D(I,J)                      !urban
             GZ1OZ0(I,J) = GZ1OZ0_URB2D(I,J)                  !urban
!m             AKHS(I,J) = AKHS_URB2D(I,J)                    !urban
             AKHS(I,J) = CHS(I,J)                             !urban
             AKMS(I,J) = AKMS_URB2D(I,J)                      !urban
           END IF                                             !urban
         ENDDO                                                !urban
       ENDDO                                                  !urban
     ENDIF
!------------------------------------------------------------------

       ELSE
         CALL wrf_error_fatal3 ( "module_surface_driver.b" , 798 , 'Lacking arguments for LSM in surface driver')
       ENDIF

     CASE (RUCLSMSCHEME)
       IF (PRESENT(qv_curr)    .AND.  PRESENT(qc_curr) .AND.    &
!           PRESENT(emiss)      .AND.  PRESENT(t2)      .AND.    &
           PRESENT(qsg)        .AND.  PRESENT(qvg)     .AND.    &
           PRESENT(qcg)        .AND.  PRESENT(soilt1)  .AND.    &
           PRESENT(tsnav)      .AND.  PRESENT(smfr3d)  .AND.    &
           PRESENT(keepfr3dflag) .AND. PRESENT(rainbl) .AND.    &
                                                      .TRUE. ) THEN

           IF( PRESENT(sr) ) THEN
               frpcpn=.true.
           ELSE
               SR = 1.
           ENDIF

           CALL wrf_debug(100,'in RUC LSM')
           CALL LSMRUC(dtbl,itimestep,num_soil_layers,          &
                zs,rainbl,snow,snowh,snowc,sr,frpcpn,           &
                dz8w,p8w,t_phy,qv_curr,qc_curr,rho,             & !p8w in [pa]
                glw,gsw,emiss,chklowq,                          &
                flqc,flhc,mavail,canwat,vegfra,albedo,znt,      &
                snoalb, albbck,                                 &   !new
                qsfc,qsg,qvg,qcg,soilt1,tsnav,                  &
                tmn,ivgtyp,isltyp,xland,xice,                   &
                cp,g,xlv,stbolt,                                &
                smois,smstav,smstot,tslb,tsk,hfx,qfx,lh,        &
                sfcrunoff,udrunoff,sfcexc,                      &
                sfcevp,grdflx,acsnow,                           &
                smfr3d,keepfr3dflag,                            &
                myj,                                            &
                ids,ide, jds,jde, kds,kde,                      &
                ims,ime, jms,jme, kms,kme,                      &
                i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )

         IF(.not. MYJ) then

          CALL SFCDIAGS(HFX,QFX,TSK,QVG,CHS2,CQS2,T2,TH2,Q2,      &
                     PSFC,CP,R_d,RCP,                              &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
             i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )
         ENDIF
 

       ELSE
         CALL wrf_error_fatal3 ( "module_surface_driver.b" , 846 , 'Lacking arguments for RUCLSM in surface driver')
       ENDIF

     CASE DEFAULT

       WRITE( message , * ) &
        'The surface option does not exist: sf_surface_physics = ', sf_surface_physics
       CALL wrf_error_fatal3 ( "module_surface_driver.b" , 853 ,  message )

     END SELECT sfc_select
     ENDDO
     !$OMP END PARALLEL DO

 430 CONTINUE


! Reset RAINBL in mm (Accumulation between PBL calls)

     IF ( PRESENT( rainbl ) ) THEN
       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij, i, j, k )
       DO ij = 1 , num_tiles
         DO j=j_start(ij),j_end(ij)
         DO i=i_start(ij),i_end(ij)
            RAINBL(i,j) = 0.
         ENDDO
         ENDDO
       ENDDO
       !$OMP END PARALLEL DO
     ENDIF

   ENDIF

   END SUBROUTINE surface_driver

END MODULE module_surface_driver
