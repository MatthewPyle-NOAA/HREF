!WRF:MEDIATION_LAYER:PHYSICS
!

MODULE module_fddagd_driver
CONTAINS

!------------------------------------------------------------------
   SUBROUTINE fddagd_driver(itimestep,dt,xtime,                   &
                  id,  &
                  RUNDGDTEN,RVNDGDTEN,RTHNDGDTEN,                 &
                  RQVNDGDTEN,RMUNDGDTEN,                          &
                  u_ndg_old,v_ndg_old,t_ndg_old,q_ndg_old,mu_ndg_old,       &
                  u_ndg_new,v_ndg_new,t_ndg_new,q_ndg_new,mu_ndg_new,       &
                  u3d,v3d,th_phy,rho,moist,                       &
                  p_phy,pi_phy,p8w,t_phy,dz8w,z,z_at_w,           &
                  config_flags,DX,n_moist,                        &
                  STEPFG,                                         &
                  pblh,ht,                                        &
                  ids,ide, jds,jde, kds,kde,                      &
                  ims,ime, jms,jme, kms,kme,                      &
                  i_start,i_end, j_start,j_end, kts,kte, num_tiles)
!------------------------------------------------------------------
   USE module_configure
   USE module_state_description
   USE module_model_constants

! *** add new modules of schemes here

   USE module_fdda_psufddagd
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
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
!-- RUNDGDTEN       U tendency due to 
!                 FDDA analysis nudging (m/s^2)
!-- RVNDGDTEN       V tendency due to 
!                 FDDA analysis nudging (m/s^2)
!-- RTHNDGDTEN      Theta tendency due to 
!                 FDDA analysis nudging (K/s)
!-- RQVNDGDTEN      Qv tendency due to 
!                 FDDA analysis nudging (kg/kg/s)
!-- RMUNDGDTEN      mu tendency due to 
!                 FDDA analysis nudging (Pa/s)
!-- itimestep     number of time steps
!-- u3d           u-velocity staggered on u points (m/s)
!-- v3d           v-velocity staggered on v points (m/s)
!-- th_phy        potential temperature (K)
!-- moist         moisture array (4D - last index is species) (kg/kg)
!-- p_phy         pressure (Pa)
!-- pi_phy        exner function (dimensionless)
!-- p8w           pressure at full levels (Pa)
!-- t_phy         temperature (K)
!-- dz8w          dz between full levels (m)
!-- z             height above sea level (m)
!-- config_flags
!-- DX            horizontal space interval (m)
!-- DT            time step (second)
!-- n_moist       number of moisture species
!-- STEPFG        number of timesteps per FDDA re-calculation
!-- KPBL          k-index of PBL top
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
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!
!******************************************************************
!------------------------------------------------------------------ 
   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags
!

   INTEGER , INTENT(IN)         ::     id

   INTEGER,    INTENT(IN   )    ::     ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       kts,kte, num_tiles,        &
                                       n_moist           

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                   &
  &                                    i_start,i_end,j_start,j_end

   INTEGER,    INTENT(IN   )    ::     itimestep,STEPFG
!
   REAL,       INTENT(IN   )    ::     DT,DX,XTIME


!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(IN   )    ::                         p_phy, &
                                                          pi_phy, &
                                                             p8w, &
                                                             rho, &
                                                           t_phy, &
                                                             u3d, &
                                                             v3d, &
                                                            dz8w, &
                                                               z, &
                                                          z_at_w, &
                                                          th_phy
!
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, n_moist ),         &
         INTENT(IN ) ::                                    moist
!
!
!
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(INOUT)    ::                       RUNDGDTEN, &
                                                         RVNDGDTEN, &
                                                        RTHNDGDTEN, &
                                                        RQVNDGDTEN

   REAL,       DIMENSION( ims:ime,  jms:jme ),            &
               INTENT(INOUT)    ::                      RMUNDGDTEN

   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ),            &
               INTENT(INOUT)    ::                       u_ndg_old, &
                                                         v_ndg_old, &
                                                         t_ndg_old, &
                                                         q_ndg_old, &
                                                         u_ndg_new, &
                                                         v_ndg_new, &
                                                         t_ndg_new, &
                                                         q_ndg_new
   REAL,       DIMENSION( ims:ime,  jms:jme ),            &
               INTENT(INOUT)    ::                       mu_ndg_old, &
                                                         mu_ndg_new

!
   REAL,    DIMENSION( ims:ime , jms:jme ),     &
               INTENT(IN   ) ::           pblh, &
                                            ht

!  LOCAL  VAR

!
   INTEGER :: i,J,K,NK,jj,ij

!------------------------------------------------------------------
!
!
   END SUBROUTINE fddagd_driver
END MODULE module_fddagd_driver