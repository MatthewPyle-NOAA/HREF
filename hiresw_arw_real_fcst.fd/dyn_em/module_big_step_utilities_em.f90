!WRF:MODEL_LAYER:DYNAMICS
!



MODULE module_big_step_utilities_em

   USE module_domain
   USE module_model_constants
   USE module_state_description
   USE module_configure
   USE module_wrf_error

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE calc_mu_uv ( config_flags,                 &
                        mu, mub, muu, muv,            &
                        ids, ide, jds, jde, kds, kde, &
                        ims, ime, jms, jme, kms, kme, &
                        its, ite, jts, jte, kts, kte )

   IMPLICIT NONE
   
   ! Input data

   TYPE(grid_config_rec_type   ) ,   INTENT(IN   ) :: config_flags

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte 

   REAL, DIMENSION( ims:ime , jms:jme ) , INTENT(  OUT) :: muu, muv
   REAL, DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mu, mub

   !  local stuff

   INTEGER :: i, j, itf, jtf, im, jm

!<DESCRIPTION>
!
!  calc_mu_uv calculates the full column dry-air mass at the staggered
!  horizontal velocity points (u,v) and places the results in muu and muv.
!  This routine uses the reference state (mub) and perturbation state (mu)
!
!</DESCRIPTION>


      itf=ite
      jtf=MIN(jte,jde-1)

      IF      ( ( its .NE. ids ) .AND. ( ite .NE. ide ) ) THEN
         DO j=jts,jtf
         DO i=its,itf
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j)+mub(i,j)+mub(i-1,j))
         ENDDO
         ENDDO
      ELSE IF ( ( its .EQ. ids ) .AND. ( ite .NE. ide ) ) THEN
         DO j=jts,jtf
         DO i=its+1,itf
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j)+mub(i,j)+mub(i-1,j))
         ENDDO
         ENDDO
         i=its
         im = its
         if(config_flags%periodic_x) im = its-1
         DO j=jts,jtf
!            muu(i,j) =      mu(i,j)          +mub(i,j)
!  fix for periodic b.c., 13 march 2004, wcs
            muu(i,j) = 0.5*(mu(i,j)+mu(im,j)+mub(i,j)+mub(im,j))
         ENDDO
      ELSE IF ( ( its .NE. ids ) .AND. ( ite .EQ. ide ) ) THEN
         DO j=jts,jtf
         DO i=its,itf-1
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j)+mub(i,j)+mub(i-1,j))
         ENDDO
         ENDDO
         i=ite
         im = ite-1
         if(config_flags%periodic_x) im = ite
         DO j=jts,jtf
!            muu(i,j) =      mu(i-1,j)        +mub(i-1,j)
!  fix for periodic b.c., 13 march 2004, wcs
            muu(i,j) = 0.5*(mu(i-1,j)+mu(im,j)+mub(i-1,j)+mub(im,j))
         ENDDO
      ELSE IF ( ( its .EQ. ids ) .AND. ( ite .EQ. ide ) ) THEN
         DO j=jts,jtf
         DO i=its+1,itf-1
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j)+mub(i,j)+mub(i-1,j))
         ENDDO
         ENDDO
         i=its
         im = its
         if(config_flags%periodic_x) im = its-1
         DO j=jts,jtf
!            muu(i,j) =      mu(i,j)          +mub(i,j)
!  fix for periodic b.c., 13 march 2004, wcs
            muu(i,j) = 0.5*(mu(i,j)+mu(im,j)+mub(i,j)+mub(im,j))
         ENDDO
         i=ite
         im = ite-1
         if(config_flags%periodic_x) im = ite
         DO j=jts,jtf
!            muu(i,j) =      mu(i-1,j)        +mub(i-1,j)
!  fix for periodic b.c., 13 march 2004, wcs
            muu(i,j) = 0.5*(mu(i-1,j)+mu(im,j)+mub(i-1,j)+mub(im,j))
         ENDDO
      END IF

      itf=MIN(ite,ide-1)
      jtf=jte

      IF      ( ( jts .NE. jds ) .AND. ( jte .NE. jde ) ) THEN
         DO j=jts,jtf
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1)+mub(i,j)+mub(i,j-1))
         ENDDO
         ENDDO
      ELSE IF ( ( jts .EQ. jds ) .AND. ( jte .NE. jde ) ) THEN
         DO j=jts+1,jtf
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1)+mub(i,j)+mub(i,j-1))
         ENDDO
         ENDDO
         j=jts
         jm = jts
         if(config_flags%periodic_y) jm = jts-1
         DO i=its,itf
!             muv(i,j) =      mu(i,j)          +mub(i,j)
!  fix for periodic b.c., 13 march 2004, wcs
             muv(i,j) = 0.5*(mu(i,j)+mu(i,jm)+mub(i,j)+mub(i,jm))
         ENDDO
      ELSE IF ( ( jts .NE. jds ) .AND. ( jte .EQ. jde ) ) THEN
         DO j=jts,jtf-1
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1)+mub(i,j)+mub(i,j-1))
         ENDDO
         ENDDO
         j=jte
         jm = jte-1
         if(config_flags%periodic_y) jm = jte
         DO i=its,itf
             muv(i,j) =      mu(i,j-1)        +mub(i,j-1)
!  fix for periodic b.c., 13 march 2004, wcs
             muv(i,j) = 0.5*(mu(i,j-1)+mu(i,jm)+mub(i,j-1)+mub(i,jm))
         ENDDO
      ELSE IF ( ( jts .EQ. jds ) .AND. ( jte .EQ. jde ) ) THEN
         DO j=jts+1,jtf-1
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1)+mub(i,j)+mub(i,j-1))
         ENDDO
         ENDDO
         j=jts
         jm = jts
         if(config_flags%periodic_y) jm = jts-1
         DO i=its,itf
!             muv(i,j) =      mu(i,j)          +mub(i,j)
!  fix for periodic b.c., 13 march 2004, wcs
             muv(i,j) = 0.5*(mu(i,j)+mu(i,jm)+mub(i,j)+mub(i,jm))
         ENDDO
         j=jte
         jm = jte-1
         if(config_flags%periodic_y) jm = jte
         DO i=its,itf
!             muv(i,j) =      mu(i,j-1)        +mub(i,j-1)
!  fix for periodic b.c., 13 march 2004, wcs
             muv(i,j) = 0.5*(mu(i,j-1)+mu(i,jm)+mub(i,j-1)+mub(i,jm))
         ENDDO
      END IF

END SUBROUTINE calc_mu_uv

!-------------------------------------------------------------------------------

SUBROUTINE calc_mu_uv_1 ( config_flags,                 &
                          mu, muu, muv,                 &
                          ids, ide, jds, jde, kds, kde, &
                          ims, ime, jms, jme, kms, kme, &
                          its, ite, jts, jte, kts, kte )

   IMPLICIT NONE
   
   ! Input data

   TYPE(grid_config_rec_type   ) ,   INTENT(IN   ) :: config_flags

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte 

   REAL, DIMENSION( ims:ime , jms:jme ) , INTENT(  OUT) :: muu, muv
   REAL, DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mu

   !  local stuff

   INTEGER :: i, j, itf, jtf, im, jm

!<DESCRIPTION>
!
!  calc_mu_uv calculates the full column dry-air mass at the staggered
!  horizontal velocity points (u,v) and places the results in muu and muv.
!  This routine uses the full state (mu)
!
!</DESCRIPTION>
   
      itf=ite
      jtf=MIN(jte,jde-1)

      IF      ( ( its .NE. ids ) .AND. ( ite .NE. ide ) ) THEN
         DO j=jts,jtf
         DO i=its,itf
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j))
         ENDDO
         ENDDO
      ELSE IF ( ( its .EQ. ids ) .AND. ( ite .NE. ide ) ) THEN
         DO j=jts,jtf
         DO i=its+1,itf
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j))
         ENDDO
         ENDDO
         i=its
         im = its
         if(config_flags%periodic_x) im = its-1
         DO j=jts,jtf
            muu(i,j) = 0.5*(mu(i,j)+mu(im,j))
         ENDDO
      ELSE IF ( ( its .NE. ids ) .AND. ( ite .EQ. ide ) ) THEN
         DO j=jts,jtf
         DO i=its,itf-1
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j))
         ENDDO
         ENDDO
         i=ite
         im = ite-1
         if(config_flags%periodic_x) im = ite
         DO j=jts,jtf
            muu(i,j) = 0.5*(mu(i-1,j)+mu(im,j))
         ENDDO
      ELSE IF ( ( its .EQ. ids ) .AND. ( ite .EQ. ide ) ) THEN
         DO j=jts,jtf
         DO i=its+1,itf-1
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j))
         ENDDO
         ENDDO
         i=its
         im = its
         if(config_flags%periodic_x) im = its-1
         DO j=jts,jtf
            muu(i,j) = 0.5*(mu(i,j)+mu(im,j))
         ENDDO
         i=ite
         im = ite-1
         if(config_flags%periodic_x) im = ite
         DO j=jts,jtf
            muu(i,j) = 0.5*(mu(i-1,j)+mu(im,j))
         ENDDO
      END IF

      itf=MIN(ite,ide-1)
      jtf=jte

      IF      ( ( jts .NE. jds ) .AND. ( jte .NE. jde ) ) THEN
         DO j=jts,jtf
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1))
         ENDDO
         ENDDO
      ELSE IF ( ( jts .EQ. jds ) .AND. ( jte .NE. jde ) ) THEN
         DO j=jts+1,jtf
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1))
         ENDDO
         ENDDO
         j=jts
         jm = jts
         if(config_flags%periodic_y) jm = jts-1
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,jm))
         ENDDO
      ELSE IF ( ( jts .NE. jds ) .AND. ( jte .EQ. jde ) ) THEN
         DO j=jts,jtf-1
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1))
         ENDDO
         ENDDO
         j=jte
         jm = jte-1
         if(config_flags%periodic_y) jm = jte
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j-1)+mu(i,jm))
         ENDDO
      ELSE IF ( ( jts .EQ. jds ) .AND. ( jte .EQ. jde ) ) THEN
         DO j=jts+1,jtf-1
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1))
         ENDDO
         ENDDO
         j=jts
         jm = jts
         if(config_flags%periodic_y) jm = jts-1
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,jm))
         ENDDO
         j=jte
         jm = jte-1
         if(config_flags%periodic_y) jm = jte
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j-1)+mu(i,jm))
         ENDDO
      END IF

END SUBROUTINE calc_mu_uv_1

!-------------------------------------------------------------------------------

SUBROUTINE couple_momentum ( muu, ru, u, msfu,              &
                             muv, rv, v, msfv,              &
                             mut, rw, w, msft,              &
                             ids, ide, jds, jde, kds, kde,  &
                             ims, ime, jms, jme, kms, kme,  &
                             its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE

   ! Input data

   INTEGER ,             INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                          ims, ime, jms, jme, kms, kme, &
                                          its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(  OUT) :: ru, rv, rw

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: muu, muv, mut
   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: msfu, msfv, msft
   
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u, v, w
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   
!<DESCRIPTION>
!
! couple_momentum couples the velocities to the full column mass and
! the map factors.
!
!</DESCRIPTION>

   ktf=MIN(kte,kde-1)
   
      itf=ite
      jtf=MIN(jte,jde-1)

      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
         ru(i,k,j)=u(i,k,j)*muu(i,j)/msfu(i,j)
      ENDDO
      ENDDO
      ENDDO

      itf=MIN(ite,ide-1)
      jtf=jte

      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
           rv(i,k,j)=v(i,k,j)*muv(i,j)/msfv(i,j)
      ENDDO
      ENDDO
      ENDDO

      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)

      DO j=jts,jtf
      DO k=kts,kte
      DO i=its,itf
         rw(i,k,j)=w(i,k,j)*mut(i,j)/msft(i,j)
      ENDDO
      ENDDO
      ENDDO

END SUBROUTINE couple_momentum

!-------------------------------------------------------------------

SUBROUTINE calc_mu_staggered ( mu, mub, muu, muv,            &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte )

   IMPLICIT NONE
   
   ! Input data

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte 

   REAL, DIMENSION( ims:ime , jms:jme ) , INTENT(  OUT) :: muu, muv
   REAL, DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mu, mub

   !  local stuff

   INTEGER :: i, j, itf, jtf

!<DESCRIPTION>
!
! calc_mu_staggered calculates the full dry air mass at the staggered
! velocity points (u,v).
!
!</DESCRIPTION>
   
      itf=ite
      jtf=MIN(jte,jde-1)

      IF      ( ( its .NE. ids ) .AND. ( ite .NE. ide ) ) THEN
         DO j=jts,jtf
         DO i=its,itf
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j)+mub(i,j)+mub(i-1,j))
         ENDDO
         ENDDO
      ELSE IF ( ( its .EQ. ids ) .AND. ( ite .NE. ide ) ) THEN
         DO j=jts,jtf
         DO i=its+1,itf
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j)+mub(i,j)+mub(i-1,j))
         ENDDO
         ENDDO
         i=its
         DO j=jts,jtf
            muu(i,j) =      mu(i,j)          +mub(i,j)
         ENDDO
      ELSE IF ( ( its .NE. ids ) .AND. ( ite .EQ. ide ) ) THEN
         DO j=jts,jtf
         DO i=its,itf-1
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j)+mub(i,j)+mub(i-1,j))
         ENDDO
         ENDDO
         i=ite
         DO j=jts,jtf
            muu(i,j) =      mu(i-1,j)        +mub(i-1,j)
         ENDDO
      ELSE IF ( ( its .EQ. ids ) .AND. ( ite .EQ. ide ) ) THEN
         DO j=jts,jtf
         DO i=its+1,itf-1
            muu(i,j) = 0.5*(mu(i,j)+mu(i-1,j)+mub(i,j)+mub(i-1,j))
         ENDDO
         ENDDO
         i=its
         DO j=jts,jtf
            muu(i,j) =      mu(i,j)          +mub(i,j)
         ENDDO
         i=ite
         DO j=jts,jtf
            muu(i,j) =      mu(i-1,j)        +mub(i-1,j)
         ENDDO
      END IF

      itf=MIN(ite,ide-1)
      jtf=jte

      IF      ( ( jts .NE. jds ) .AND. ( jte .NE. jde ) ) THEN
         DO j=jts,jtf
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1)+mub(i,j)+mub(i,j-1))
         ENDDO
         ENDDO
      ELSE IF ( ( jts .EQ. jds ) .AND. ( jte .NE. jde ) ) THEN
         DO j=jts+1,jtf
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1)+mub(i,j)+mub(i,j-1))
         ENDDO
         ENDDO
         j=jts
         DO i=its,itf
             muv(i,j) =      mu(i,j)          +mub(i,j)
         ENDDO
      ELSE IF ( ( jts .NE. jds ) .AND. ( jte .EQ. jde ) ) THEN
         DO j=jts,jtf-1
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1)+mub(i,j)+mub(i,j-1))
         ENDDO
         ENDDO
         j=jte
         DO i=its,itf
             muv(i,j) =      mu(i,j-1)        +mub(i,j-1)
         ENDDO
      ELSE IF ( ( jts .EQ. jds ) .AND. ( jte .EQ. jde ) ) THEN
         DO j=jts+1,jtf-1
         DO i=its,itf
             muv(i,j) = 0.5*(mu(i,j)+mu(i,j-1)+mub(i,j)+mub(i,j-1))
         ENDDO
         ENDDO
         j=jts
         DO i=its,itf
             muv(i,j) =      mu(i,j)          +mub(i,j)
         ENDDO
         j=jte
         DO i=its,itf
             muv(i,j) =      mu(i,j-1)        +mub(i,j-1)
         ENDDO
      END IF

END SUBROUTINE calc_mu_staggered

!-------------------------------------------------------------------------------

SUBROUTINE couple ( mu, mub, rfield, field, name, &
                    msf,                          &
                    ids, ide, jds, jde, kds, kde, &
                    ims, ime, jms, jme, kms, kme, &
                    its, ite, jts, jte, kts, kte )

   IMPLICIT NONE

   ! Input data

   INTEGER ,             INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                          ims, ime, jms, jme, kms, kme, &
                                          its, ite, jts, jte, kts, kte

   CHARACTER(LEN=1) ,     INTENT(IN   ) :: name

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(  OUT) :: rfield

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mu, mub, msf
   
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: field
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   REAL , DIMENSION(ims:ime,jms:jme) :: muu , muv

!<DESCRIPTION>
!
! subroutine couple couples the input variable with the dry-air 
! column mass (mu).  
!
!</DESCRIPTION>

   
   ktf=MIN(kte,kde-1)
   
   IF (name .EQ. 'u')THEN

      CALL calc_mu_staggered ( mu, mub, muu, muv,            &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte )

      itf=ite
      jtf=MIN(jte,jde-1)

      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
         rfield(i,k,j)=field(i,k,j)*muu(i,j)/msf(i,j)
      ENDDO
      ENDDO
      ENDDO

   ELSE IF (name .EQ. 'v')THEN

      CALL calc_mu_staggered ( mu, mub, muu, muv,            &
                               ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kts, kte )

      itf=ite
      itf=MIN(ite,ide-1)
      jtf=jte

      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
           rfield(i,k,j)=field(i,k,j)*muv(i,j)/msf(i,j)
      ENDDO
      ENDDO
      ENDDO

   ELSE IF (name .EQ. 'w')THEN
      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
      DO j=jts,jtf
      DO k=kts,kte
      DO i=its,itf
         rfield(i,k,j)=field(i,k,j)*(mu(i,j)+mub(i,j))/msf(i,j)
      ENDDO
      ENDDO
      ENDDO

   ELSE IF (name .EQ. 'h')THEN
      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
      DO j=jts,jtf
      DO k=kts,kte
      DO i=its,itf
         rfield(i,k,j)=field(i,k,j)*(mu(i,j)+mub(i,j))
      ENDDO
      ENDDO
      ENDDO

   ELSE 
      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
         rfield(i,k,j)=field(i,k,j)*(mu(i,j)+mub(i,j))
      ENDDO
      ENDDO
      ENDDO
   
   ENDIF

END SUBROUTINE couple

!-----------------------------------------------------------------------

SUBROUTINE calc_ww ( mu, ru, rv, ww,               &
                     rdx, rdy, msft, dnw,          &
                     ids, ide, jds, jde, kds, kde, &
                     ims, ime, jms, jme, kms, kme, &
                     its, ite, jts, jte, kts, kte )

   IMPLICIT NONE

   ! Input data


   INTEGER ,                                   INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                                                ims, ime, jms, jme, kms, kme, &
                                                                its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: ru, rv
   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mu, msft
   REAL , DIMENSION( kms:kme ) , INTENT(IN   ) :: dnw
   
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT  ) :: ww
   REAL , INTENT(IN   )  :: rdx, rdy
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   REAL , DIMENSION( its:ite ) :: dmdt

!<DESCRIPTION>
!
!  calc_ww calculates omega using the mass-coupled velocities mu*u, mu*v.
!  The algorithm integrates the continuity equation through the column
!  followed by a diagnosis of omega.
!
!</DESCRIPTION>


    jtf=MIN(jte,jde-1)
    ktf=MIN(kte,kde-1)  
    itf=MIN(ite,ide-1)

      DO j=jts,jtf

        DO i=its,ite
          dmdt(i) = 0.
          ww(i,1,j) = 0.
          ww(i,kte,j) = 0.
        ENDDO

!!        DO k=kts,ktf+1

        DO k=kts,ktf
        DO i=its,itf

          dmdt(i) = dmdt(i) + dnw(k)* ( rdx*(ru(i+1,k,j)-ru(i,k,j))  &
                                       +rdy*(rv(i,k,j+1)-rv(i,k,j))   )

        ENDDO
        ENDDO

!               DO K=2,NZ-1
!                  ww(K,I)=ww(K-1,I)-DNW(K-1)*
!     &                  (DMDT+RDX*( xmu(i  )*u(K,I  ) 
!     &                             -xmu(im1)*u(k,im1)) )
!               END DO

        DO k=2,ktf
        DO i=its,itf

           ww(i,k,j)=ww(i,k-1,j)                                       &
                        - dnw(k-1)* ( dmdt(i)                          &
                                     +rdx*(ru(i+1,k-1,j)-ru(i,k-1,j))  &
                                     +rdy*(rv(i,k-1,j+1)-rv(i,k-1,j)) )
        ENDDO
        ENDDO
     ENDDO

END SUBROUTINE calc_ww


!-------------------------------------------------------------------------------

SUBROUTINE calc_ww_cp ( u, v, mup, mub, ww,              &
                        rdx, rdy, msft, msfu, msfv, dnw, &
                        ids, ide, jds, jde, kds, kde,    &
                        ims, ime, jms, jme, kms, kme,    &
                        its, ite, jts, jte, kts, kte    )

   IMPLICIT NONE

   ! Input data


   INTEGER ,    INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: u, v
   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mup, mub, &
                                                            msft, msfu, msfv
   REAL , DIMENSION( kms:kme ) , INTENT(IN   ) :: dnw
   
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT  ) :: ww
   REAL , INTENT(IN   )  :: rdx, rdy
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   REAL , DIMENSION( its:ite ) :: dmdt
   REAL , DIMENSION( its:ite, kts:kte ) :: divv
   REAL , DIMENSION( its:ite+1, jts:jte+1 ) :: muu, muv

!<DESCRIPTION>
!
!  calc_ww calculates omega using the velocities (u,v) and the dry-air 
!  column mass (mup+mub).
!  The algorithm integrates the continuity equation through the column
!  followed by a diagnosis of omega.
!
!</DESCRIPTION>

!<DESCRIPTION>
!
!  calc_ww_cp calculates omega using the velocities (u,v) and the 
!  column mass mu.
!
!</DESCRIPTION>

    jtf=MIN(jte,jde-1)
    ktf=MIN(kte,kde-1)  
    itf=MIN(ite,ide-1)

!  mu coupled with the appropriate map factor

      DO j=jts,jtf
      DO i=its,min(ite+1,ide)
        muu(i,j) = 0.5*(mup(i,j)+mub(i,j)+mup(i-1,j)+mub(i-1,j))/msfu(i,j)
      ENDDO
      ENDDO

      DO j=jts,min(jte+1,jde)
      DO i=its,itf
        muv(i,j) = 0.5*(mup(i,j)+mub(i,j)+mup(i,j-1)+mub(i,j-1))/msfv(i,j)
      ENDDO
      ENDDO

      DO j=jts,jtf

        DO i=its,ite
          dmdt(i) = 0.
          ww(i,1,j) = 0.
          ww(i,kte,j) = 0.
        ENDDO

        DO k=kts,ktf
        DO i=its,itf

          divv(i,k) = msft(i,j)*dnw(k)*( rdx*(muu(i+1,j)*u(i+1,k,j)-muu(i,j)*u(i,k,j))  &
                                        +rdy*(muv(i,j+1)*v(i,k,j+1)-muv(i,j)*v(i,k,j))   )

!          dmdt(i) = dmdt(i) + dnw(k)* ( rdx*(ru(i+1,k,j)-ru(i,k,j))  &
!                                       +rdy*(rv(i,k,j+1)-rv(i,k,j))   )

          dmdt(i) = dmdt(i) + divv(i,k)


        ENDDO
        ENDDO

        DO k=2,ktf
        DO i=its,itf

!           ww(i,k,j)=ww(i,k-1,j)                                       &
!                        - dnw(k-1)* ( dmdt(i)                          &
!                                     +rdx*(ru(i+1,k-1,j)-ru(i,k-1,j))  &
!                                     +rdy*(rv(i,k-1,j+1)-rv(i,k-1,j)) )

           ww(i,k,j)=ww(i,k-1,j) - dnw(k-1)*dmdt(i) - divv(i,k-1)

        ENDDO
        ENDDO
     ENDDO


END SUBROUTINE calc_ww_cp


!-------------------------------------------------------------------------------
 
SUBROUTINE calc_cq ( moist, cqu, cqv, cqw, n_moist, &
                     ids, ide, jds, jde, kds, kde,  &
                     ims, ime, jms, jme, kms, kme,  &
                     its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE
   
   ! Input data

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte 

   INTEGER ,          INTENT(IN   ) :: n_moist
   

   REAL, DIMENSION( ims:ime, kms:kme , jms:jme , n_moist ), INTENT(IN   ) :: moist
                                              
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(  OUT) :: cqu, cqv, cqw

   ! Local stuff

   REAL :: qtot
   
   INTEGER :: i, j, k, itf, jtf, ktf, ispe

!<DESCRIPTION>
!
!  calc_cq calculates moist coefficients for the momentum equations.
!
!</DESCRIPTION>

      itf=ite
      jtf=MIN(jte,jde-1)
      ktf=MIN(kte,kde-1)

      IF(  n_moist >= PARAM_FIRST_SCALAR ) THEN

        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
          qtot = 0.
!DEC$ loop count(3)
          DO ispe=PARAM_FIRST_SCALAR,n_moist
            qtot = qtot + moist(i,k,j,ispe) + moist(i-1,k,j,ispe)
          ENDDO
!           qtot = 0.5*( moist(i  ,k,j,1)+moist(i  ,k,j,2)+moist(i  ,k,j,3)+  &
!     &                  moist(i-1,k,j,1)+moist(i-1,k,j,2)+moist(i-1,k,j,3) )
!           cqu(i,k,j) = 1./(1.+qtot)
           cqu(i,k,j) = 1./(1.+0.5*qtot)
        ENDDO
        ENDDO
        ENDDO

        itf=MIN(ite,ide-1)
        jtf=jte

        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
          qtot = 0.
!DEC$ loop count(3)
          DO ispe=PARAM_FIRST_SCALAR,n_moist
            qtot = qtot + moist(i,k,j,ispe) + moist(i,k,j-1,ispe)
          ENDDO
!           qtot = 0.5*( moist(i,k,j  ,1)+moist(i,k,j  ,2)+moist(i,k,j  ,3)+  &
!     &                  moist(i,k,j-1,1)+moist(i,k,j-1,2)+moist(i,k,j-1,3) )
!           cqv(i,k,j) = 1./(1.+qtot)
           cqv(i,k,j) = 1./(1.+0.5*qtot)
        ENDDO
        ENDDO
        ENDDO

        itf=MIN(ite,ide-1)
        jtf=MIN(jte,jde-1)
        DO j=jts,jtf
        DO k=kts+1,ktf
        DO i=its,itf
          qtot = 0.
!DEC$ loop count(3)
          DO ispe=PARAM_FIRST_SCALAR,n_moist
            qtot = qtot + moist(i,k,j,ispe) + moist(i,k-1,j,ispe)
          ENDDO
!           qtot = 0.5*( moist(i,k  ,j,1)+moist(i,k  ,j,2)+moist(i,k-1,j,3)+  &
!     &                  moist(i,k-1,j,1)+moist(i,k-1,j,2)+moist(i,k  ,j,3) )
!           cqw(i,k,j) = qtot
           cqw(i,k,j) = 0.5*qtot
        ENDDO
        ENDDO
        ENDDO

      ELSE

        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           cqu(i,k,j) = 1.
        ENDDO
        ENDDO
        ENDDO

        itf=MIN(ite,ide-1)
        jtf=jte

        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
           cqv(i,k,j) = 1.
        ENDDO
        ENDDO
        ENDDO

        itf=MIN(ite,ide-1)
        jtf=MIN(jte,jde-1)
        DO j=jts,jtf
        DO k=kts+1,ktf
        DO i=its,itf
           cqw(i,k,j) = 0.
        ENDDO
        ENDDO
        ENDDO

      END IF

END SUBROUTINE calc_cq

!----------------------------------------------------------------------

SUBROUTINE calc_alt ( alt, al, alb,                  &
                      ids, ide, jds, jde, kds, kde,  &
                      ims, ime, jms, jme, kms, kme,  &
                      its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE
   
   ! Input data

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte 

   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: alb, al
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(  OUT) :: alt

   ! Local stuff

   INTEGER :: i, j, k, itf, jtf, ktf

!<DESCRIPTION>
!
! calc_alt computes the full inverse density
!
!</DESCRIPTION>

      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
      ktf=MIN(kte,kde-1)

      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
        alt(i,k,j) = al(i,k,j)+alb(i,k,j)
      ENDDO
      ENDDO
      ENDDO


END SUBROUTINE calc_alt

!----------------------------------------------------------------------

SUBROUTINE calc_p_rho_phi ( moist, n_moist,                &
                            al, alb, mu, muts, ph, p, pb,  &
                            t, p0, t0, znu, dnw, rdnw,     &
                            rdn, non_hydrostatic,          &
                            ids, ide, jds, jde, kds, kde,  &
                            ims, ime, jms, jme, kms, kme,  &
                            its, ite, jts, jte, kts, kte  )

  IMPLICIT NONE
   
   ! Input data

  LOGICAL ,          INTENT(IN   ) :: non_hydrostatic

  INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte 

  INTEGER ,          INTENT(IN   ) :: n_moist

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: alb,  &
                                                                   pb,   &
                                                                   t

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme, n_moist ), INTENT(IN   ) :: moist

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(  OUT) :: al, p

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(INOUT) :: ph

  REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN   ) :: mu, muts

  REAL, DIMENSION( kms:kme ), INTENT(IN   ) :: znu, dnw, rdnw, rdn

  REAL,   INTENT(IN   ) :: t0, p0

  ! Local stuff

  INTEGER :: i, j, k, itf, jtf, ktf, ispe
  REAL    :: qvf, qtot, qf1, qf2
  REAL, DIMENSION( its:ite) :: temp,cpovcv_v


!<DESCRIPTION>
!
! For the nonhydrostatic option, calc_p_rho_phi calculates the
! diagnostic quantities pressure and (inverse) density from the
! prognostic variables using the equation of state.
!
! For the hydrostatic option, calc_p_rho_phi calculates the
! diagnostic quantities (inverse) density and geopotential from the
! prognostic variables using the equation of state and the hydrostatic 
! equation.
!
!</DESCRIPTION>

  itf=MIN(ite,ide-1)
  jtf=MIN(jte,jde-1)
  ktf=MIN(kte,kde-1)

  cpovcv_v = cpovcv

  IF (non_hydrostatic) THEN

      IF (n_moist >= PARAM_FIRST_SCALAR ) THEN  

        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
          qvf = 1.+rvovrd*moist(i,k,j,P_QV)
          al(i,k,j)=-1./muts(i,j)*(alb(i,k,j)*mu(i,j)  &
                     +rdnw(k)*(ph(i,k+1,j)-ph(i,k,j)))
          temp(i)=(r_d*(t0+t(i,k,j))*qvf)/                 &
                        (p0*(al(i,k,j)+alb(i,k,j)))
        ENDDO
! use vector version from libmassv or from compat lib in frame/libmassv.F
        CALL vspow  ( p(its,k,j), temp(its), cpovcv_v(its), itf-its+1 )
        DO i=its,itf
           p(i,k,j)= p(i,k,j)*p0-pb(i,k,j)
        ENDDO
        ENDDO
        ENDDO

      ELSE

        DO j=jts,jtf
        DO k=kts,ktf
        DO i=its,itf
          al(i,k,j)=-1./muts(i,j)*(alb(i,k,j)*mu(i,j)  &
                     +rdnw(k)*(ph(i,k+1,j)-ph(i,k,j)))
          p(i,k,j)=p0*( (r_d*(t0+t(i,k,j)))/                     &
                        (p0*(al(i,k,j)+alb(i,k,j))) )**cpovcv  &
                           -pb(i,k,j)
        ENDDO
        ENDDO
        ENDDO

      END IF

   ELSE

!  hydrostatic pressure, al, and ph1 calc; WCS, 5 sept 2001


      IF (n_moist >= PARAM_FIRST_SCALAR ) THEN  

        DO j=jts,jtf

          k=ktf          ! top layer
          DO i=its,itf

            qtot = 0.
            DO ispe=PARAM_FIRST_SCALAR,n_moist
              qtot = qtot + moist(i,k,j,ispe)
            ENDDO
            qf2 = 1./(1.+qtot)
            qf1 = qtot*qf2

            p(i,k,j) = - 0.5*(mu(i,j)+qf1*muts(i,j))/rdnw(k)/qf2
            qvf = 1.+rvovrd*moist(i,k,j,P_QV)
            al(i,k,j) = (r_d/p1000mb)*(t(i,k,j)+t0)*qvf* &
                (((p(i,k,j)+pb(i,k,j))/p1000mb)**cvpm) - alb(i,k,j)

          ENDDO

          DO k=ktf-1,kts,-1  ! remaining layers, integrate down
            DO i=its,itf

            qtot = 0.
            DO ispe=PARAM_FIRST_SCALAR,n_moist
              qtot = qtot + 0.5*(  moist(i,k  ,j,ispe) + moist(i,k+1,j,ispe) )
            ENDDO
            qf2 = 1./(1.+qtot)
            qf1 = qtot*qf2

            p(i,k,j) = p(i,k+1,j) - (mu(i,j) + qf1*muts(i,j))/qf2/rdn(k+1)
            qvf = 1.+rvovrd*moist(i,k,j,P_QV)
            al(i,k,j) = (r_d/p1000mb)*(t(i,k,j)+t0)*qvf* &
                        (((p(i,k,j)+pb(i,k,j))/p1000mb)**cvpm) - alb(i,k,j)
            ENDDO
          ENDDO

          DO k=2,ktf+1  ! integrate hydrostatic equation for geopotential
            DO i=its,itf

!              ph(i,k,j) = ph(i,k-1,j) - (1./rdnw(k-1))*(       &
!                           (muts(i,j)+mu(i,j))*al(i,k-1,j)+    &
!                            mu(i,j)*alb(i,k-1,j)  )
              ph(i,k,j) = ph(i,k-1,j) - (dnw(k-1))*(           &
                           (muts(i,j))*al(i,k-1,j)+            &
                            mu(i,j)*alb(i,k-1,j)  )
                                                   

            ENDDO
          ENDDO

        ENDDO

      ELSE

        DO j=jts,jtf

          k=ktf          ! top layer
          DO i=its,itf

            qtot = 0.
            qf2 = 1./(1.+qtot)
            qf1 = qtot*qf2

            p(i,k,j) = - 0.5*(mu(i,j)+qf1*muts(i,j))/rdnw(k)/qf2
            qvf = 1.
            al(i,k,j) = (r_d/p1000mb)*(t(i,k,j)+t0)*qvf* &
                (((p(i,k,j)+pb(i,k,j))/p1000mb)**cvpm) - alb(i,k,j)

          ENDDO

          DO k=ktf-1,kts,-1  ! remaining layers, integrate down
            DO i=its,itf

            qtot = 0.
            qf2 = 1./(1.+qtot)
            qf1 = qtot*qf2

            p(i,k,j) = p(i,k+1,j) - (mu(i,j) + qf1*muts(i,j))/qf2/rdn(k+1)
            qvf = 1.
            al(i,k,j) = (r_d/p1000mb)*(t(i,k,j)+t0)*qvf* &
                        (((p(i,k,j)+pb(i,k,j))/p1000mb)**cvpm) - alb(i,k,j)
            ENDDO
          ENDDO

          DO k=2,ktf+1  ! integrate hydrostatic equation for geopotential
            DO i=its,itf

!              ph(i,k,j) = ph(i,k-1,j) - (1./rdnw(k-1))*(       &
!                           (muts(i,j)+mu(i,j))*al(i,k-1,j)+    &
!                            mu(i,j)*alb(i,k-1,j)  )
              ph(i,k,j) = ph(i,k-1,j) - (dnw(k-1))*(           &
                           (muts(i,j))*al(i,k-1,j)+            &
                            mu(i,j)*alb(i,k-1,j)  )
                                                   

            ENDDO
          ENDDO

        ENDDO

     END IF

   END IF

END SUBROUTINE calc_p_rho_phi

!----------------------------------------------------------------------

SUBROUTINE calc_php ( php, ph, phb,                  &
                      ids, ide, jds, jde, kds, kde,  &
                      ims, ime, jms, jme, kms, kme,  &
                      its, ite, jts, jte, kts, kte  )

   IMPLICIT NONE
   
   ! Input data

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte 

   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(IN   ) :: phb, ph
   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(  OUT) :: php

   ! Local stuff

   INTEGER :: i, j, k, itf, jtf, ktf

!<DESCRIPTION>
!
!  calc_php calculates the full geopotential from the reference state
!  geopotential and the perturbation geopotential (phb_ph).
!
!</DESCRIPTION>

      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
      ktf=MIN(kte,kde-1)

      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
        php(i,k,j) = 0.5*(phb(i,k,j)+phb(i,k+1,j)+ph(i,k,j)+ph(i,k+1,j))
      ENDDO
      ENDDO
      ENDDO

END SUBROUTINE calc_php

!-------------------------------------------------------------------------------

SUBROUTINE diagnose_w( ph_tend, ph_new, ph_old, w, mu, dt,  &
                       u, v, ht,                            &
                       cf1, cf2, cf3, rdx, rdy, msft,       &
                       ids, ide, jds, jde, kds, kde,        &
                       ims, ime, jms, jme, kms, kme,        &
                       its, ite, jts, jte, kts, kte        )

   IMPLICIT NONE

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte 

   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(IN   ) ::   ph_tend, &
                                                                     ph_new,  &
                                                                     ph_old,  &
                                                                     u,       &
                                                                     v


   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(  OUT) :: w

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: mu, ht, msft

   REAL, INTENT(IN   ) :: dt, cf1, cf2, cf3, rdx, rdy

   INTEGER :: i, j, k, itf, jtf

   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)

!<DESCRIPTION>
!
! diagnose_w diagnoses the vertical velocity from the geopoential equation.
! Used with the hydrostatic option.
!
!</DESCRIPTION>

   DO j = jts, jtf

!  lower b.c. on w

     DO i = its, itf
         w(i,1,j)=  msft(i,j)*(                              &
                  .5*rdy*(                                   &
                           (ht(i,j+1)-ht(i,j  ))             &
          *(cf1*v(i,1,j+1)+cf2*v(i,2,j+1)+cf3*v(i,3,j+1))    &
                          +(ht(i,j  )-ht(i,j-1))             &
          *(cf1*v(i,1,j  )+cf2*v(i,2,j  )+cf3*v(i,3,j  ))  ) &
                 +.5*rdx*(                                   &
                           (ht(i+1,j)-ht(i,j  ))             &
          *(cf1*u(i+1,1,j)+cf2*u(i+1,2,j)+cf3*u(i+1,3,j))    &
                          +(ht(i,j  )-ht(i-1,j))             &
          *(cf1*u(i  ,1,j)+cf2*u(i  ,2,j)+cf3*u(i  ,3,j))  ) &
                                                            )
     ENDDO

!  use geopotential equation to diagnose w

     DO k = 2, kte
     DO i = its, itf
       w(i,k,j) =  msft(i,j)*(  (ph_new(i,k,j)-ph_old(i,k,j))/dt       &
                               - ph_tend(i,k,j)/mu(i,j)        )/g 

     ENDDO
     ENDDO

   ENDDO

END SUBROUTINE diagnose_w

!-------------------------------------------------------------------------------

SUBROUTINE rhs_ph( ph_tend, u, v, ww,               &
                   ph, ph_old, phb, w,              &
                   mut, muu, muv,                   &
                   fnm, fnp,                        &
                   rdnw, cfn, cfn1, rdx, rdy, msft, &
                   non_hydrostatic,                 &
                   config_flags,                    &
                   ids, ide, jds, jde, kds, kde,    &
                   ims, ime, jms, jme, kms, kme,    &
                   its, ite, jts, jte, kts, kte    )
   IMPLICIT NONE

   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte 

   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(IN   ) ::        &
                                                                     u,   &
                                                                     v,   &
                                                                     ww,  &
                                                                     ph,  &
                                                                     ph_old, &
                                                                     phb, & 
                                                                    w

! pjj/cray
!  REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(  OUT) :: ph_tend
   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(INOUT) :: ph_tend

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: muu, muv, mut, msft

   REAL, DIMENSION( kms:kme ), INTENT(IN   ) :: rdnw, fnm, fnp

   REAL,  INTENT(IN   ) :: cfn, cfn1, rdx, rdy

   LOGICAL,  INTENT(IN   )  ::  non_hydrostatic

   ! Local stuff

   INTEGER :: i, j, k, itf, jtf, ktf, kz, i_start, j_start
   REAL    :: ur, ul, ub, vr, vl, vb
   REAL, DIMENSION(its:ite,kts:kte) :: wdwn

   INTEGER :: advective_order

   LOGICAL :: specified

!<DESCRIPTION>
!
! rhs_ph calculates the large-timestep tendency terms for the geopotential
! equation.  These terms include the advection and "gw".  The geopotential
! equation is cast in advective form, so we dont use the flux form advection
! algorithms here.
!
!</DESCRIPTION>

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   advective_order = config_flags%h_sca_adv_order 
!   advective_order = 2  !  original configuration (pre Oct 2001)

   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)
   ktf=MIN(kte,kde-1)

! advective form for the geopotential equation

   DO j = jts, jtf

     DO k = 2, kte
     DO i = its, itf
          wdwn(i,k) = .5*(ww(i,k,j)+ww(i,k-1,j))*rdnw(k-1)               &
                        *(ph(i,k,j)-ph(i,k-1,j)+phb(i,k,j)-phb(i,k-1,j))
     ENDDO
     ENDDO

     DO k = 2, kte-1
     DO i = its, itf
           ph_tend(i,k,j) = ph_tend(i,k,j)                           &
                             - (fnm(k)*wdwn(i,k+1)+fnp(k)*wdwn(i,k))
     ENDDO
     ENDDO

   ENDDO

   IF (non_hydrostatic) THEN  ! add in "gw" term.
   DO j = jts, jtf            ! in hydrostatic mode, "gw" will be diagnosed
                              ! after the timestep to give us "w"
     DO i = its, itf
        ph_tend(i,kde,j) = 0.
     ENDDO

     DO k = 2, kte
     DO i = its, itf
        ph_tend(i,k,j) = ph_tend(i,k,j) + mut(i,j)*g*w(i,k,j)/msft(i,j)
     ENDDO
     ENDDO

   ENDDO

   END IF

   IF (advective_order <= 2) THEN

!  y (v) advection

   i_start = its
   j_start = jts
   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)

   IF ( (config_flags%open_ys) .and. jts == jds ) j_start = jts+1
   IF ( (config_flags%open_ye) .and. jte == jde ) jtf = jtf-1

   DO j = j_start, jtf

     DO k = 2, kte-1
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .25*rdy*                       &
                 ( muv(i,j+1)*(v(i,k,j+1)+v(i,k-1,j+1))*               &
                  (phb(i,k,j+1)-phb(i,k,j  )+ph(i,k,j+1)-ph(i,k,j  ))  &
                  +muv(i,j  )*(v(i,k,j  )+v(i,k-1,j  ))*               &
                  (phb(i,k,j  )-phb(i,k,j-1)+ph(i,k,j  )-ph(i,k,j-1)) )
     ENDDO
     ENDDO

     k = kte
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .5*rdy*                         &
                  ( muv(i,j+1)*(cfn*v(i,k-1,j+1)+cfn1*v(i,k-2,j+1))*    &
                   (phb(i,k,j+1)-phb(i,k,j  )+ph(i,k,j+1)-ph(i,k,j  ))  &
                   +muv(i,j  )*(cfn*v(i,k-1,j  )+cfn1*v(i,k-2,j  ))*    &
                   (phb(i,k,j  )-phb(i,k,j-1)+ph(i,k,j  )-ph(i,k,j-1)) )
     ENDDO

   ENDDO

!  x (u) advection

   i_start = its
   j_start = jts
   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)

   IF ( (config_flags%open_xs) .and. its == ids ) i_start = its+1
   IF ( (config_flags%open_xe) .and. ite == ide ) itf = itf-1

   DO j = j_start, jtf

     DO k = 2, kte-1
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .25*rdx*                        &
                 ( muu(i+1,j)*(u(i+1,k,j)+u(i+1,k-1,j))*                &
                  (phb(i+1,k,j)-phb(i  ,k,j)+ph(i+1,k,j)-ph(i  ,k,j))   &
                  +muu(i  ,j)*(u(i  ,k,j)+u(i  ,k-1,j))*                &
                  (phb(i  ,k,j)-phb(i-1,k,j)+ph(i  ,k,j)-ph(i-1,k,j)) )
     ENDDO
     ENDDO
 
     k = kte
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .5*rdx*                         &
                  ( muu(i+1,j)*(cfn*u(i+1,k-1,j)+cfn1*u(i+1,k-2,j))*    &
                   (phb(i+1,k,j)-phb(i  ,k,j)+ph(i+1,k,j)-ph(i  ,k,j))  &
                   +muu(i  ,j)*(cfn*u(i  ,k-1,j)+cfn1*u(i  ,k-2,j))*    &
                   (phb(i  ,k,j)-phb(i-1,k,j)+ph(i  ,k,j)-ph(i-1,k,j)) )
     ENDDO

   ENDDO

   ELSE IF (advective_order <= 4) THEN

!  y (v) advection

   i_start = its
   j_start = jts
   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)

   IF ( (config_flags%open_ys) .and. jts == jds ) j_start = jts+1
   IF ( (config_flags%open_ye) .and. jte == jde ) jtf = jtf-1

   DO j = j_start, jtf

     DO k = 2, kte-1
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .25*rdy*           (          &
                 ( muv(i,j+1)*(v(i,k,j+1)+v(i,k-1,j+1))               &
                  +muv(i,j  )*(v(i,k,j  )+v(i,k-1,j  )) )* (1./12.)*( &
                    8.*(ph(i,k,j+1)-ph(i,k,j-1))                      &
                      -(ph(i,k,j+2)-ph(i,k,j-2))                      &
                   +8.*(phb(i,k,j+1)-phb(i,k,j-1))                    &
                      -(phb(i,k,j+2)-phb(i,k,j-2))  )   )                


     ENDDO
     ENDDO

     k = kte
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .5*rdy*           (                      &
                 ( muv(i,j+1)*(cfn*v(i,k-1,j+1)+cfn1*v(i,k-2,j+1))               &
                  +muv(i,j  )*(cfn*v(i,k-1,j  )+cfn1*v(i,k-2,j  )) )* (1./12.)*( &
                    8.*(ph(i,k,j+1)-ph(i,k,j-1))                                 &
                      -(ph(i,k,j+2)-ph(i,k,j-2))                                 &
                   +8.*(phb(i,k,j+1)-phb(i,k,j-1))                               &
                      -(phb(i,k,j+2)-phb(i,k,j-2))  )   )                

     ENDDO

   ENDDO


!  x (u) advection

   i_start = its
   j_start = jts
   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)

   IF ( (config_flags%open_xs) .and. its == ids ) i_start = its+1
   IF ( (config_flags%open_xe) .and. ite == ide ) itf = itf-1

   DO j = j_start, jtf

     DO k = 2, kte-1
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .25*rdx*(                     &
                 ( muu(i+1,j)*(u(i+1,k,j)+u(i+1,k-1,j))               &
                  +muu(i,j  )*(u(i,k,j  )+u(i,k-1,j  )) )* (1./12.)*( &
                    8.*(ph(i+1,k,j)-ph(i-1,k,j))                      &
                      -(ph(i+2,k,j)-ph(i-2,k,j))                      &
                   +8.*(phb(i+1,k,j)-phb(i-1,k,j))                    &
                      -(phb(i+2,k,j)-phb(i-2,k,j))  )   )                
     ENDDO
     ENDDO
 
     k = kte
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .5*rdx*(                                 &
                 ( muu(i+1,j)*(cfn*u(i+1,k-1,j)+cfn1*u(i+1,k-2,j))               &
                  +muu(i,j  )*(cfn*u(i  ,k-1,j)+cfn1*u(i,k-2,j)) )* (1./12.)*(   &
                    8.*(ph(i+1,k,j)-ph(i-1,k,j))                                 &
                      -(ph(i+2,k,j)-ph(i-2,k,j))                                 &
                   +8.*(phb(i+1,k,j)-phb(i-1,k,j))                               &
                      -(phb(i+2,k,j)-phb(i-2,k,j))  )     )
     ENDDO

   ENDDO

   ELSE IF (advective_order <= 6) THEN

!  y (v) advection

   i_start = its
   j_start = jts
   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)

!   IF ( (config_flags%open_ys) .and. jts == jds ) j_start = jts+1
!   IF ( (config_flags%open_ye) .and. jte == jde ) jtf = jtf-1

   IF (config_flags%open_ys .or. specified ) j_start = max(jts,jds+2)
   IF (config_flags%open_ye .or. specified ) jtf     = min(jtf,jde-3)

   DO j = j_start, jtf

     DO k = 2, kte-1
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .25*rdy* (                    &
                 ( muv(i,j+1)*(v(i,k,j+1)+v(i,k-1,j+1))               &
                  +muv(i,j  )*(v(i,k,j  )+v(i,k-1,j  )) )* (1./60.)*( &
                   45.*(ph(i,k,j+1)-ph(i,k,j-1))                      &
                   -9.*(ph(i,k,j+2)-ph(i,k,j-2))                      &
                      +(ph(i,k,j+3)-ph(i,k,j-3))                      &
                  +45.*(phb(i,k,j+1)-phb(i,k,j-1))                    &
                   -9.*(phb(i,k,j+2)-phb(i,k,j-2))                    &
                      +(phb(i,k,j+3)-phb(i,k,j-3))  )   )                


     ENDDO
     ENDDO

     k = kte
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .5*rdy* (                                &
                 ( muv(i,j+1)*(cfn*v(i,k-1,j+1)+cfn1*v(i,k-2,j+1))               &
                  +muv(i,j  )*(cfn*v(i,k-1,j  )+cfn1*v(i,k-2,j  )) )* (1./60.)*( &
                   45.*(ph(i,k,j+1)-ph(i,k,j-1))                                 &
                   -9.*(ph(i,k,j+2)-ph(i,k,j-2))                                 &
                      +(ph(i,k,j+3)-ph(i,k,j-3))                                 &
                  +45.*(phb(i,k,j+1)-phb(i,k,j-1))                               &
                   -9.*(phb(i,k,j+2)-phb(i,k,j-2))                               &
                      +(phb(i,k,j+3)-phb(i,k,j-3))  )   )                

     ENDDO

   ENDDO


!  pick up near boundary rows using 4th order stencil 
!  (open bc copy only goes out to jds-1 and jde, hence 4rth is ok but 6th is too big)

   IF ( (config_flags%open_ys) .and. jts <= jds+1 )  THEN

     j = jds+1
     DO k = 2, kte-1
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .25*rdy* (                    &
                 ( muv(i,j+1)*(v(i,k,j+1)+v(i,k-1,j+1))               &
                  +muv(i,j  )*(v(i,k,j  )+v(i,k-1,j  )) )* (1./12.)*( &
                    8.*(ph(i,k,j+1)-ph(i,k,j-1))                      &
                      -(ph(i,k,j+2)-ph(i,k,j-2))                      &
                   +8.*(phb(i,k,j+1)-phb(i,k,j-1))                    &
                      -(phb(i,k,j+2)-phb(i,k,j-2))  )   )                


     ENDDO
     ENDDO

     k = kte
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .5*rdy* (                                &
                 ( muv(i,j+1)*(cfn*v(i,k-1,j+1)+cfn1*v(i,k-2,j+1))               &
                  +muv(i,j  )*(cfn*v(i,k-1,j  )+cfn1*v(i,k-2,j  )) )* (1./12.)*( &
                    8.*(ph(i,k,j+1)-ph(i,k,j-1))                                 &
                      -(ph(i,k,j+2)-ph(i,k,j-2))                                 &
                   +8.*(phb(i,k,j+1)-phb(i,k,j-1))                               &
                      -(phb(i,k,j+2)-phb(i,k,j-2))  )   )                

     ENDDO

   END IF

   IF ( (config_flags%open_ye) .and. jte >= jde-2 )  THEN

     j = jde-2
     DO k = 2, kte-1
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .25*rdy* (                    &
                 ( muv(i,j+1)*(v(i,k,j+1)+v(i,k-1,j+1))               &
                  +muv(i,j  )*(v(i,k,j  )+v(i,k-1,j  )) )* (1./12.)*( &
                    8.*(ph(i,k,j+1)-ph(i,k,j-1))                      &
                      -(ph(i,k,j+2)-ph(i,k,j-2))                      &
                   +8.*(phb(i,k,j+1)-phb(i,k,j-1))                    &
                      -(phb(i,k,j+2)-phb(i,k,j-2))  )   )                


     ENDDO
     ENDDO

     k = kte
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .5*rdy* (                                &
                 ( muv(i,j+1)*(cfn*v(i,k-1,j+1)+cfn1*v(i,k-2,j+1))               &
                  +muv(i,j  )*(cfn*v(i,k-1,j  )+cfn1*v(i,k-2,j  )) )* (1./12.)*( &
                    8.*(ph(i,k,j+1)-ph(i,k,j-1))                                 &
                      -(ph(i,k,j+2)-ph(i,k,j-2))                                 &
                   +8.*(phb(i,k,j+1)-phb(i,k,j-1))                               &
                      -(phb(i,k,j+2)-phb(i,k,j-2))  )   )                

     ENDDO

   END IF

!  x (u) advection

   i_start = its
   j_start = jts
   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)

   IF (config_flags%open_xs .or. specified ) i_start = max(its,ids+2)
   IF (config_flags%open_xe .or. specified ) itf     = min(itf,ide-3)
   IF ( config_flags%periodic_x ) i_start = its
   IF ( config_flags%periodic_x ) itf=MIN(ite,ide-1)

   DO j = j_start, jtf

     DO k = 2, kte-1
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .25*rdx*(                     &
                 ( muu(i+1,j)*(u(i+1,k,j)+u(i+1,k-1,j))               &
                  +muu(i,j  )*(u(i,k,j  )+u(i,k-1,j  )) )* (1./60.)*( &
                   45.*(ph(i+1,k,j)-ph(i-1,k,j))                      &
                   -9.*(ph(i+2,k,j)-ph(i-2,k,j))                      &
                      +(ph(i+3,k,j)-ph(i-3,k,j))                      &
                  +45.*(phb(i+1,k,j)-phb(i-1,k,j))                    &
                   -9.*(phb(i+2,k,j)-phb(i-2,k,j))                    &
                      +(phb(i+3,k,j)-phb(i-3,k,j))  )   )                
     ENDDO
     ENDDO
 
     k = kte
     DO i = i_start, itf
        ph_tend(i,k,j)=ph_tend(i,k,j) - .5*rdx*(                                 &
                 ( muu(i+1,j)*(cfn*u(i+1,k-1,j)+cfn1*u(i+1,k-2,j))               &
                  +muu(i,j  )*(cfn*u(i  ,k-1,j)+cfn1*u(i,k-2,j)) )* (1./60.)*(   &
                   45.*(ph(i+1,k,j)-ph(i-1,k,j))                                 &
                   -9.*(ph(i+2,k,j)-ph(i-2,k,j))                                 &
                      +(ph(i+3,k,j)-ph(i-3,k,j))                                 &
                  +45.*(phb(i+1,k,j)-phb(i-1,k,j))                               &
                   -9.*(phb(i+2,k,j)-phb(i-2,k,j))                               &
                      +(phb(i+3,k,j)-phb(i-3,k,j))  )     )
     ENDDO

   ENDDO

   IF ( (config_flags%open_xs) .and. its <= ids+1 ) THEN
     i = ids + 1
     DO j = j_start, jtf
       DO k = 2, kte-1
        ph_tend(i,k,j)=ph_tend(i,k,j) - .25*rdx*(                     &
                 ( muu(i+1,j)*(u(i+1,k,j)+u(i+1,k-1,j))               &
                  +muu(i,j  )*(u(i,k,j  )+u(i,k-1,j  )) )* (1./12.)*( &
                    8.*(ph(i+1,k,j)-ph(i-1,k,j))                      &
                      -(ph(i+2,k,j)-ph(i-2,k,j))                      &
                   +8.*(phb(i+1,k,j)-phb(i-1,k,j))                    &
                      -(phb(i+2,k,j)-phb(i-2,k,j))  )   )                
       ENDDO
       k = kte
       ph_tend(i,k,j)=ph_tend(i,k,j) - .5*rdx*(                                 &
                ( muu(i+1,j)*(cfn*u(i+1,k-1,j)+cfn1*u(i+1,k-2,j))               &
                 +muu(i,j  )*(cfn*u(i  ,k-1,j)+cfn1*u(i,k-2,j)) )* (1./12.)*(   &
                   8.*(ph(i+1,k,j)-ph(i-1,k,j))                                 &
                     -(ph(i+2,k,j)-ph(i-2,k,j))                                 &
                  +8.*(phb(i+1,k,j)-phb(i-1,k,j))                               &
                     -(phb(i+2,k,j)-phb(i-2,k,j))  )     )

     ENDDO
   END IF

   IF ( (config_flags%open_xe) .and. ite >= ide-2 ) THEN
     i = ide-2
     DO j = j_start, jtf
       DO k = 2, kte-1
        ph_tend(i,k,j)=ph_tend(i,k,j) - .25*rdx*(                     &
                 ( muu(i+1,j)*(u(i+1,k,j)+u(i+1,k-1,j))               &
                  +muu(i,j  )*(u(i,k,j  )+u(i,k-1,j  )) )* (1./12.)*( &
                    8.*(ph(i+1,k,j)-ph(i-1,k,j))                      &
                      -(ph(i+2,k,j)-ph(i-2,k,j))                      &
                   +8.*(phb(i+1,k,j)-phb(i-1,k,j))                    &
                      -(phb(i+2,k,j)-phb(i-2,k,j))  )   )                
       ENDDO
       k = kte
       ph_tend(i,k,j)=ph_tend(i,k,j) - .5*rdx*(                                 &
                ( muu(i+1,j)*(cfn*u(i+1,k-1,j)+cfn1*u(i+1,k-2,j))               &
                 +muu(i,j  )*(cfn*u(i  ,k-1,j)+cfn1*u(i,k-2,j)) )* (1./12.)*(   &
                   8.*(ph(i+1,k,j)-ph(i-1,k,j))                                 &
                     -(ph(i+2,k,j)-ph(i-2,k,j))                                 &
                  +8.*(phb(i+1,k,j)-phb(i-1,k,j))                               &
                     -(phb(i+2,k,j)-phb(i-2,k,j))  )     )

     ENDDO
   END IF

   END IF

!  lateral open boundary conditions,
!  start with north and south (y) boundaries

   i_start = its
   itf=MIN(ite,ide-1)

   !  south

   IF ( (config_flags%open_ys) .and. jts == jds ) THEN

     j=jts

     DO k=2,kde
       kz = min(k,kde-1)
       DO i = its,itf
         vb =.5*( fnm(kz)*(v(i,kz  ,j+1)+v(i,kz  ,j  ))    &
                 +fnp(kz)*(v(i,kz-1,j+1)+v(i,kz-1,j  )) )
         vl=amin1(vb,0.)
         ph_tend(i,k,j)=ph_tend(i,k,j)-rdy*mut(i,j)*(      &
                              +vl*(ph_old(i,k,j+1)-ph_old(i,k,j)))
       ENDDO
     ENDDO

   END IF

   ! north

   IF ( (config_flags%open_ye) .and. jte == jde ) THEN

     j=jte-1

     DO k=2,kde
       kz = min(k,kde-1)
       DO i = its,itf
        vb=.5*( fnm(kz)*(v(i,kz  ,j+1)+v(i,kz  ,j))   &
               +fnp(kz)*(v(i,kz-1,j+1)+v(i,kz-1,j)) )
        vr=amax1(vb,0.)
        ph_tend(i,k,j)=ph_tend(i,k,j)-rdy*mut(i,j)*(      &
                   +vr*(ph_old(i,k,j)-ph_old(i,k,j-1)))
       ENDDO
     ENDDO

   END IF

   !  now the east and west (y) boundaries

   j_start = its
   jtf=MIN(jte,jde-1)

   !  west

   IF ( (config_flags%open_xs) .and. its == ids ) THEN

     i=its

     DO j = jts,jtf
       DO k=2,kde-1
         kz = k
         ub =.5*( fnm(kz)*(u(i+1,kz  ,j)+u(i  ,kz  ,j))     &
                 +fnp(kz)*(u(i+1,kz-1,j)+u(i  ,kz-1,j)) )
         ul=amin1(ub,0.)
         ph_tend(i,k,j)=ph_tend(i,k,j)-rdx*mut(i,j)*(       &
                              +ul*(ph_old(i+1,k,j)-ph_old(i,k,j)))
       ENDDO

         k = kde
         kz = k
         ub =.5*( fnm(kz)*(u(i+1,kz  ,j)+u(i  ,kz  ,j))     &
                 +fnp(kz)*(u(i+1,kz-1,j)+u(i  ,kz-1,j)) )
         ul=amin1(ub,0.)
         ph_tend(i,k,j)=ph_tend(i,k,j)-rdx*mut(i,j)*(       &
                              +ul*(ph_old(i+1,k,j)-ph_old(i,k,j)))
     ENDDO

   END IF

   ! east

   IF ( (config_flags%open_xe) .and. ite == ide ) THEN

     i = ite-1

     DO j = jts,jtf
       DO k=2,kde-1
        kz = k
        ub=.5*( fnm(kz)*(u(i+1,kz  ,j)+u(i,kz  ,j))  &
               +fnp(kz)*(u(i+1,kz-1,j)+u(i,kz-1,j)) )
        ur=amax1(ub,0.)
        ph_tend(i,k,j)=ph_tend(i,k,j)-rdx*mut(i,j)*( &
                   +ur*(ph_old(i,k,j)-ph_old(i-1,k,j)))
       ENDDO

        k = kde    
        kz = k-1
        ub=.5*( fnm(kz)*(u(i+1,kz  ,j)+u(i,kz  ,j))   &
               +fnp(kz)*(u(i+1,kz-1,j)+u(i,kz-1,j)) )
        ur=amax1(ub,0.)
        ph_tend(i,k,j)=ph_tend(i,k,j)-rdx*mut(i,j)*(  &
                   +ur*(ph_old(i,k,j)-ph_old(i-1,k,j)))

     ENDDO

   END IF

  END SUBROUTINE rhs_ph

!-------------------------------------------------------------------------------

SUBROUTINE horizontal_pressure_gradient( ru_tend,rv_tend,                &
                                         ph,alt,p,pb,al,php,cqu,cqv,     &
                                         muu,muv,mu,fnm,fnp,rdnw,        &
                                         cf1,cf2,cf3,rdx,rdy,msft,       &
                                         config_flags, non_hydrostatic,  &
                                         ids, ide, jds, jde, kds, kde,   &
                                         ims, ime, jms, jme, kms, kme,   &
                                         its, ite, jts, jte, kts, kte   )

   IMPLICIT NONE
   
   ! Input data


   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   LOGICAL, INTENT (IN   ) :: non_hydrostatic

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte 

   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(IN   ) ::        &
                                                                     ph,  &
                                                                     alt, &
                                                                     al,  &
                                                                     p,   &
                                                                     pb,  &
                                                                     php, &
                                                                     cqu, &
                                                                     cqv


   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(INOUT) ::           &
                                                                    ru_tend, &
                                                                    rv_tend

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: muu, muv, mu, msft

   REAL, DIMENSION( kms:kme ), INTENT(IN   ) :: rdnw, fnm, fnp

   REAL,  INTENT(IN   ) :: rdx, rdy, cf1, cf2, cf3

   INTEGER :: i,j,k, itf, jtf, ktf, i_start, j_start
   REAL, DIMENSION( ims:ime, kms:kme ) :: dpn
   REAL :: dpx, dpy

   LOGICAL :: specified

!<DESCRIPTION>
!
!  horizontal_pressure_gradient calculates the 
!  horizontal pressure gradient terms for the large-timestep tendency 
!  in the horizontal momentum equations (u,v).
!
!</DESCRIPTION>

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

! start with the north-south (y) pressure gradient

   itf=MIN(ite,ide-1)
   jtf=jte
   ktf=MIN(kte,kde-1)
   i_start = its
   j_start = jts
   IF ( (config_flags%open_ys .or. specified .or. &
           config_flags%nested ) .and. jts == jds ) j_start = jts+1
   IF ( (config_flags%open_ye .or. specified .or. &
           config_flags%nested ) .and. jte == jde ) jtf = jtf-1

   DO j = j_start, jtf

     IF ( non_hydrostatic )  THEN

        k=1

        DO i = i_start, itf
          dpn(i,k) = .5*( cf1*(p(i,k  ,j-1)+p(i,k  ,j))   &
                         +cf2*(p(i,k+1,j-1)+p(i,k+1,j))   &
                         +cf3*(p(i,k+2,j-1)+p(i,k+2,j))  )
          dpn(i,kde) = 0.
        ENDDO
               
        DO k=2,ktf
          DO i = i_start, itf
            dpn(i,k) = .5*( fnm(k)*(p(i,k  ,j-1)+p(i,k  ,j))  &
                           +fnp(k)*(p(i,k-1,j-1)+p(i,k-1,j)) )
          END DO
        END DO

        DO K=1,ktf
          DO i = i_start, itf
            dpy = .5*rdy*muv(i,j)*(                                            &
                     (ph (i,k+1,j)-ph (i,k+1,j-1) + ph(i,k,j)-ph(i,k,j-1))  &
                    +(alt(i,k  ,j)+alt(i,k  ,j-1))*(p (i,k,j)-p (i,k,j-1))  &
                    +(al (i,k  ,j)+al (i,k  ,j-1))*(pb(i,k,j)-pb(i,k,j-1)) )
            dpy = dpy + rdy*(php(i,k,j)-php(i,k,j-1))*              &
                (rdnw(k)*(dpn(i,k+1)-dpn(i,k))-.5*(mu(i,j-1)+mu(i,j)))
            rv_tend(i,k,j) = rv_tend(i,k,j)-cqv(i,k,j)*dpy
          END DO
        END DO

     ELSE

        DO K=1,ktf
          DO i = i_start, itf
            dpy = .5*rdy*muv(i,j)*(                                            &
                     (ph (i,k+1,j)-ph (i,k+1,j-1) + ph(i,k,j)-ph(i,k,j-1))  &
                    +(alt(i,k  ,j)+alt(i,k  ,j-1))*(p (i,k,j)-p (i,k,j-1))  &
                    +(al (i,k  ,j)+al (i,k  ,j-1))*(pb(i,k,j)-pb(i,k,j-1)) )
            rv_tend(i,k,j) = rv_tend(i,k,j)-cqv(i,k,j)*dpy
          END DO
        END DO

     END IF

   ENDDO

!  now the east-west (x) pressure gradient

   itf=ite
   jtf=MIN(jte,jde-1)
   ktf=MIN(kte,kde-1)
   i_start = its
   j_start = jts
   IF ( (config_flags%open_xs .or. specified .or. &
           config_flags%nested ) .and. its == ids ) i_start = its+1
   IF ( (config_flags%open_xe .or. specified .or. &
           config_flags%nested ) .and. ite == ide ) itf = itf-1
   IF ( config_flags%periodic_x ) i_start = its
   IF ( config_flags%periodic_x ) itf=ite

   DO j = j_start, jtf

     IF ( non_hydrostatic )  THEN

        k=1

        DO i = i_start, itf
          dpn(i,k) = .5*( cf1*(p(i-1,k  ,j)+p(i,k  ,j))   &
                         +cf2*(p(i-1,k+1,j)+p(i,k+1,j))   &
                         +cf3*(p(i-1,k+2,j)+p(i,k+2,j))  )
          dpn(i,kde) = 0.
        ENDDO
               
        DO k=2,ktf
          DO i = i_start, itf
            dpn(i,k) = .5*( fnm(k)*(p(i-1,k  ,j)+p(i,k  ,j))  &
                           +fnp(k)*(p(i-1,k-1,j)+p(i,k-1,j)) )
          END DO
        END DO

        DO K=1,ktf
          DO i = i_start, itf
            dpx = .5*rdx*muu(i,j)*(                                            &
                        (ph (i,k+1,j)-ph (i-1,k+1,j) + ph(i,k,j)-ph(i-1,k,j))  &
                       +(alt(i,k  ,j)+alt(i-1,k  ,j))*(p (i,k,j)-p (i-1,k,j))  &
                       +(al (i,k  ,j)+al (i-1,k  ,j))*(pb(i,k,j)-pb(i-1,k,j)) )
            dpx = dpx + rdx*(php(i,k,j)-php(i-1,k,j))*              &
                (rdnw(k)*(dpn(i,k+1)-dpn(i,k))-.5*(mu(i-1,j)+mu(i,j)))
            ru_tend(i,k,j) = ru_tend(i,k,j)-cqu(i,k,j)*dpx
          END DO
        END DO

     ELSE

        DO K=1,ktf
          DO i = i_start, itf
            dpx = .5*rdx*muu(i,j)*(                                            &
                        (ph (i,k+1,j)-ph (i-1,k+1,j) + ph(i,k,j)-ph(i-1,k,j))  &
                       +(alt(i,k  ,j)+alt(i-1,k  ,j))*(p (i,k,j)-p (i-1,k,j))  &
                       +(al (i,k  ,j)+al (i-1,k  ,j))*(pb(i,k,j)-pb(i-1,k,j)) )
            ru_tend(i,k,j) = ru_tend(i,k,j)-cqu(i,k,j)*dpx
          END DO
        END DO

     END IF

   ENDDO

END SUBROUTINE horizontal_pressure_gradient

!-------------------------------------------------------------------------------

SUBROUTINE pg_buoy_w( rw_tend, p, cqw, mu, mub,       &
                      rdnw, rdn, g, msft,             &
                      ids, ide, jds, jde, kds, kde,   &
                      ims, ime, jms, jme, kms, kme,   &
                      its, ite, jts, jte, kts, kte   )

   IMPLICIT NONE
   
   ! Input data

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte 

   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(IN   ) ::   p
   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(INOUT) ::   cqw


   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(INOUT) ::  rw_tend

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: mub, mu, msft

   REAL, DIMENSION( kms:kme ), INTENT(IN   ) :: rdnw, rdn

   REAL,  INTENT(IN   ) :: g

   INTEGER :: itf, jtf, i, j, k
   REAL    :: cq1, cq2


!<DESCRIPTION>
!
!  pg_buoy_w calculates the 
!  vertical pressure gradient and buoyancy terms for the large-timestep 
!  tendency in the vertical momentum equation.
!
!</DESCRIPTION>

!  BUOYANCY AND PRESSURE GRADIENT TERM IN W EQUATION AT TIME T

   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)

   DO j = jts,jtf

     k=kde
     DO i=its,itf
       cq1 = 1./(1.+cqw(i,k-1,j))
       cq2 = cqw(i,k-1,j)*cq1
       rw_tend(i,k,j) = rw_tend(i,k,j)+(1./msft(i,j))*g*(      &
                        cq1*2.*rdnw(k-1)*(  -p(i,k-1,j))  &
                        -mu(i,j)-cq2*mub(i,j)            )
     END DO

     DO k = 2, kde-1
     DO i = its,itf
      cq1 = 1./(1.+cqw(i,k,j))
      cq2 = cqw(i,k,j)*cq1
      cqw(i,k,j) = cq1
      rw_tend(i,k,j) = rw_tend(i,k,j)+(1./msft(i,j))*g*(      &
                       cq1*rdn(k)*(p(i,k,j)-p(i,k-1,j))  &
                       -mu(i,j)-cq2*mub(i,j)            )
     END DO
     ENDDO           


   ENDDO

END SUBROUTINE pg_buoy_w

!-------------------------------------------------------------------------------

SUBROUTINE w_damp( rw_tend, ww, w, mut, rdnw, dt,     &
                      w_damping,                      &
                      ids, ide, jds, jde, kds, kde,   &
                      ims, ime, jms, jme, kms, kme,   &
                      its, ite, jts, jte, kts, kte   )

   IMPLICIT NONE

   ! Input data

   INTEGER ,          INTENT(IN   ) :: w_damping

   INTEGER ,          INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte

   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(IN   ) ::   ww, w

   REAL, DIMENSION( ims:ime, kms:kme , jms:jme ), INTENT(INOUT) ::  rw_tend

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN   ) :: mut

   REAL, DIMENSION( kms:kme ), INTENT(IN   ) :: rdnw

   REAL, INTENT(IN)    :: dt
   REAL                :: cfl, cf_n, cf_d, maxcfl, maxdub, maxdeta

   INTEGER :: itf, jtf, i, j, k, maxi, maxj, maxk
   INTEGER :: some
   CHARACTER*512 :: temp
   CHARACTER (LEN=256) :: time_str
   CHARACTER (LEN=256) :: grid_str

!<DESCRIPTION>
!
!  w_damp computes a damping term for the vertical velocity when the
!  vertical Courant number is too large.  This was found to be preferable to 
!  decreasing the timestep or increasing the diffusion in real-data applications
!  that produced potentially-unstable large vertical velocities because of
!  unphysically large heating rates coming from the cumulus parameterization 
!  schemes run at moderately high resolutions (dx ~ O(10) km).
!
!</DESCRIPTION>

   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)

   some = 0
   maxcfl = 0.

   IF ( w_damping == 1 ) THEN
     DO j = jts,jtf

     DO k = 2, kde-1
     DO i = its,itf
! restructure to get rid of divide
        cf_n = abs(ww(i,k,j)*rdnw(k)*dt)
        cf_d = abs(mut(i,j))
        if(cf_n .gt. cf_d*w_beta )then
           cfl = abs(ww(i,k,j)/mut(i,j)*rdnw(k)*dt)
           IF ( cfl > maxcfl ) THEN
             maxcfl = cfl ; maxi = i ; maxj = j ; maxk = k 
             maxdub = w(i,k,j) ; maxdeta = -1./rdnw(k)
           ENDIF
           WRITE(temp,*)i,j,k,' cfl,w,d(eta)=',cfl,w(i,k,j),-1./rdnw(k)
           CALL wrf_debug ( 100 , TRIM(temp) )
           if ( cfl > 2. ) some = some + 1
           rw_tend(i,k,j) = rw_tend(i,k,j)-sign(1.,w(i,k,j))*w_alpha*(cfl-w_beta)*mut(i,j)
        endif
     END DO
     ENDDO
     ENDDO
   ELSE
! just print
     DO j = jts,jtf

     DO k = 2, kde-1
     DO i = its,itf
        cf_n = abs(ww(i,k,j)*rdnw(k)*dt)
        cf_d = abs(mut(i,j))
        if(cf_n .gt. cf_d*w_beta )then
           cfl = abs(ww(i,k,j)/mut(i,j)*rdnw(k)*dt)
           IF ( cfl > maxcfl ) THEN
             maxcfl = cfl ; maxi = i ; maxj = j ; maxk = k 
             maxdub = w(i,k,j) ; maxdeta = -1./rdnw(k)
           ENDIF
           WRITE(temp,*)i,j,k,' cfl,w,d(eta)=',cfl,w(i,k,j),-1./rdnw(k)
           CALL wrf_debug ( 100 , TRIM(temp) )
           if ( cfl > 2. ) some = some + 1
        endif
     END DO
     ENDDO
     ENDDO
   ENDIF
   IF ( some .GT. 0 ) THEN
     CALL get_current_time_string( time_str )
     CALL get_current_grid_name( grid_str )
     WRITE(wrf_err_message,*)some,                                            &
            ' points exceeded cfl=2 in domain '//TRIM(grid_str)//' at time '//TRIM(time_str)//' hours'
     CALL wrf_debug ( 0 , TRIM(wrf_err_message) )
     WRITE(wrf_err_message,*)'MAX AT i,j,k: ',maxi,maxj,maxk,' cfl,w,d(eta)=',maxcfl, &
                             maxdub,maxdeta
     CALL wrf_debug ( 0 , TRIM(wrf_err_message) )
   ENDIF

END SUBROUTINE w_damp

!-------------------------------------------------------------------------------

SUBROUTINE horizontal_diffusion ( name, field, tendency, mu,           &
                                  config_flags,                        &
                                  msfu, msfv, msft, khdif, xkmhd, rdx, rdy,   &
                                  ids, ide, jds, jde, kds, kde,        &
                                  ims, ime, jms, jme, kms, kme,        &
                                  its, ite, jts, jte, kts, kte        )

   IMPLICIT NONE
   
   ! Input data

   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,        INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte

   CHARACTER(LEN=1) ,                          INTENT(IN   ) :: name

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: field, xkmhd

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mu

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msfu,      &
                                                                msfv,      &
                                                                msft

   REAL ,                                      INTENT(IN   ) :: rdx,       &
                                                                rdy,       &
                                                                khdif

   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf

   INTEGER :: i_start, i_end, j_start, j_end

   REAL :: mrdx, mkrdxm, mkrdxp, &
           mrdy, mkrdym, mkrdyp
   REAL :: pr_inv

   LOGICAL :: specified

!<DESCRIPTION>
!
!  horizontal_diffusion computes the horizontal diffusion tendency
!  on model horizontal coordinate surfaces.
!
!</DESCRIPTION>

   pr_inv = 1./prandtl
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   ktf=MIN(kte,kde-1)
   
   IF (name .EQ. 'u') THEN

      i_start = its
      i_end   = ite
      j_start = jts
      j_end   = MIN(jte,jde-1)

      IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
      IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-1,ite)
      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite


      DO j = j_start, j_end
      DO k=kts,ktf
      DO i = i_start, i_end

         mkrdxm=msft(i-1,j)*mu(i-1,j)*xkmhd(i-1,k,j)*rdx
         mkrdxp=msft(i,j)*mu(i,j)*xkmhd(i,k,j)*rdx
         mrdx=msfu(i,j)*rdx
         mkrdym=0.5*(msfu(i,j)+msfu(i,j-1))*    &
                0.25*(mu(i,j)+mu(i,j-1)+mu(i-1,j-1)+mu(i-1,j))* &
                0.25*(xkmhd(i,k,j)+xkmhd(i,k,j-1)+xkmhd(i-1,k,j-1)+xkmhd(i-1,k,j))*rdy
         mkrdyp=0.5*(msfu(i,j)+msfu(i,j+1))*    &
                0.25*(mu(i,j)+mu(i,j+1)+mu(i-1,j+1)+mu(i-1,j))* &
                0.25*(xkmhd(i,k,j)+xkmhd(i,k,j+1)+xkmhd(i-1,k,j+1)+xkmhd(i-1,k,j))*rdy
         mrdy=msfu(i,j)*rdy

            tendency(i,k,j)=tendency(i,k,j)+( &
                            mrdx*(mkrdxp*(field(i+1,k,j)-field(i  ,k,j))  &
                                 -mkrdxm*(field(i  ,k,j)-field(i-1,k,j))) &
                           +mrdy*(mkrdyp*(field(i,k,j+1)-field(i,k,j  ))  &
                                 -mkrdym*(field(i,k,j  )-field(i,k,j-1))))
      ENDDO
      ENDDO
      ENDDO
   
   ELSE IF (name .EQ. 'v')THEN

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = jte

      IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
      IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

      DO j = j_start, j_end
      DO k=kts,ktf
      DO i = i_start, i_end

         mkrdxm=0.5*(msfv(i,j)+msfv(i-1,j))*    &
                0.25*(mu(i,j)+mu(i,j-1)+mu(i-1,j-1)+mu(i-1,j))* &
                0.25*(xkmhd(i,k,j)+xkmhd(i,k,j-1)+xkmhd(i-1,k,j-1)+xkmhd(i-1,k,j))*rdx
         mkrdxp=0.5*(msfv(i,j)+msfv(i+1,j))*    &
                0.25*(mu(i,j)+mu(i,j-1)+mu(i+1,j-1)+mu(i+1,j))* &
                0.25*(xkmhd(i,k,j)+xkmhd(i,k,j-1)+xkmhd(i+1,k,j-1)+xkmhd(i+1,k,j))*rdx
         mrdx=msfv(i,j)*rdx
         mkrdym=msft(i,j-1)*xkmhd(i,k,j-1)*rdy
         mkrdyp=msft(i,j)*xkmhd(i,k,j)*rdy
         mrdy=msfv(i,j)*rdy

            tendency(i,k,j)=tendency(i,k,j)+( &
                            mrdx*(mkrdxp*(field(i+1,k,j)-field(i  ,k,j))  &
                                 -mkrdxm*(field(i  ,k,j)-field(i-1,k,j))) &
                           +mrdy*(mkrdyp*(field(i,k,j+1)-field(i,k,j  ))  &
                                 -mkrdym*(field(i,k,j  )-field(i,k,j-1))))
      ENDDO
      ENDDO
      ENDDO
   
   ELSE IF (name .EQ. 'w')THEN

      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
      IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

      DO j = j_start, j_end
      DO k=kts+1,ktf
      DO i = i_start, i_end

         mkrdxm=msfu(i,j)*   &
                0.25*(mu(i,j)+mu(i-1,j)+mu(i,j)+mu(i-1,j))* &
                0.25*(xkmhd(i,k,j)+xkmhd(i-1,k,j)+xkmhd(i,k-1,j)+xkmhd(i-1,k-1,j))*rdx
         mkrdxp=msfu(i+1,j)*   &
                0.25*(mu(i+1,j)+mu(i,j)+mu(i+1,j)+mu(i,j))* &
                0.25*(xkmhd(i+1,k,j)+xkmhd(i,k,j)+xkmhd(i+1,k-1,j)+xkmhd(i,k-1,j))*rdx
         mrdx=msft(i,j)*rdx
         mkrdym=msfv(i,j)*   &
                0.25*(mu(i,j)+mu(i,j-1)+mu(i,j)+mu(i,j-1))* &
                0.25*(xkmhd(i,k,j)+xkmhd(i,k,j-1)+xkmhd(i,k-1,j)+xkmhd(i,k-1,j-1))*rdy
         mkrdyp=msfv(i,j+1)*   &
                0.25*(mu(i,j+1)+mu(i,j)+mu(i,j+1)+mu(i,j))* &
                0.25*(xkmhd(i,k,j+1)+xkmhd(i,k,j)+xkmhd(i,k-1,j+1)+xkmhd(i,k-1,j))*rdy
         mrdy=msft(i,j)*rdy

            tendency(i,k,j)=tendency(i,k,j)+( &
                            mrdx*(mkrdxp*(field(i+1,k,j)-field(i  ,k,j)) &
                                 -mkrdxm*(field(i  ,k,j)-field(i-1,k,j))) &
                           +mrdy*(mkrdyp*(field(i,k,j+1)-field(i,k,j  )) &
                                 -mkrdym*(field(i,k,j  )-field(i,k,j-1))))
      ENDDO
      ENDDO
      ENDDO
   
   ELSE


      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
      IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

      DO j = j_start, j_end
      DO k=kts,ktf
      DO i = i_start, i_end

         mkrdxm=msfu(i,j)*0.5*(xkmhd(i,k,j)+xkmhd(i-1,k,j))*0.5*(mu(i,j)+mu(i-1,j))*rdx*pr_inv
         mkrdxp=msfu(i+1,j)*0.5*(xkmhd(i+1,k,j)+xkmhd(i,k,j))*0.5*(mu(i+1,j)+mu(i,j))*rdx*pr_inv
         mrdx=msft(i,j)*rdx
         mkrdym=msfv(i,j)*0.5*(xkmhd(i,k,j)+xkmhd(i,k,j-1))*0.5*(mu(i,j)+mu(i,j-1))*rdy*pr_inv
         mkrdyp=msfv(i,j+1)*0.5*(xkmhd(i,k,j+1)+xkmhd(i,k,j))*0.5*(mu(i,j+1)+mu(i,j))*rdy*pr_inv
         mrdy=msft(i,j)*rdy

            tendency(i,k,j)=tendency(i,k,j)+( &
                            mrdx*(mkrdxp*(field(i+1,k,j)-field(i  ,k,j))  &
                                 -mkrdxm*(field(i  ,k,j)-field(i-1,k,j))) &
                           +mrdy*(mkrdyp*(field(i,k,j+1)-field(i,k,j  ))  &
                                 -mkrdym*(field(i,k,j  )-field(i,k,j-1))))
      ENDDO
      ENDDO
      ENDDO
           
   ENDIF

END SUBROUTINE horizontal_diffusion

!-----------------------------------------------------------------------------------------

SUBROUTINE horizontal_diffusion_3dmp ( name, field, tendency, mu,           &
                                       config_flags, base_3d,               &
                                       msfu, msfv, msft, khdif, xkmhd, rdx, rdy,   &
                                       ids, ide, jds, jde, kds, kde,        &
                                       ims, ime, jms, jme, kms, kme,        &
                                       its, ite, jts, jte, kts, kte        )

   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,        INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     its, ite, jts, jte, kts, kte

   CHARACTER(LEN=1) ,                          INTENT(IN   ) :: name

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: field, &
                                                                      xkmhd, &
                                                                      base_3d

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN   ) :: mu

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msfu,      &
                                                                    msfv,      &
                                                                    msft

   REAL ,                                      INTENT(IN   ) :: rdx,       &
                                                                rdy,       &
                                                                khdif

   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf

   INTEGER :: i_start, i_end, j_start, j_end

   REAL :: mrdx, mkrdxm, mkrdxp, &
           mrdy, mkrdym, mkrdyp
   REAL :: pr_inv

   LOGICAL :: specified

!<DESCRIPTION>
!
!  horizontal_diffusion_3dmp computes the horizontal diffusion tendency
!  on model horizontal coordinate surfaces.  This routine computes diffusion
!  a perturbation scalar (field-base_3d).
!
!</DESCRIPTION>

   pr_inv = 1./prandtl
   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   ktf=MIN(kte,kde-1)
   
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
      IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-2,ite)
      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-2,jte)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = MIN(ite,ide-1)

      DO j = j_start, j_end
      DO k=kts,ktf
      DO i = i_start, i_end

         mkrdxm=msfu(i,j)*0.5*(xkmhd(i,k,j)+xkmhd(i-1,k,j))*0.5*(mu(i,j)+mu(i-1,j))*rdx*pr_inv
         mkrdxp=msfu(i+1,j)*0.5*(xkmhd(i+1,k,j)+xkmhd(i,k,j))*0.5*(mu(i+1,j)+mu(i,j))*rdx*pr_inv
         mrdx=msft(i,j)*rdx
         mkrdym=msfv(i,j)*0.5*(xkmhd(i,k,j)+xkmhd(i,k,j-1))*0.5*(mu(i,j)+mu(i,j-1))*rdy*pr_inv
         mkrdyp=msfv(i,j+1)*0.5*(xkmhd(i,k,j+1)+xkmhd(i,k,j))*0.5*(mu(i,j+1)+mu(i,j))*rdy*pr_inv
         mrdy=msft(i,j)*rdy

            tendency(i,k,j)=tendency(i,k,j)+(                        &
                    mrdx*( mkrdxp*(   field(i+1,k,j)  -field(i  ,k,j)      &
                                   -base_3d(i+1,k,j)+base_3d(i  ,k,j) )    &
                          -mkrdxm*(   field(i  ,k,j)  -field(i-1,k,j)      &
                                   -base_3d(i  ,k,j)+base_3d(i-1,k,j) )  ) &
                   +mrdy*( mkrdyp*(   field(i,k,j+1)  -field(i,k,j  )      &
                                   -base_3d(i,k,j+1)+base_3d(i,k,j  ) )    &
                          -mkrdym*(   field(i,k,j  )  -field(i,k,j-1)      &
                                   -base_3d(i,k,j  )+base_3d(i,k,j-1) )  ) &
                                                                         ) 
      ENDDO
      ENDDO
      ENDDO

END SUBROUTINE horizontal_diffusion_3dmp

!-----------------------------------------------------------------------------------------

SUBROUTINE vertical_diffusion ( name, field, tendency,        &
                                config_flags,                 &
                                alt, mut, rdn, rdnw, kvdif,   &
                                ids, ide, jds, jde, kds, kde, &
                                ims, ime, jms, jme, kms, kme, &
                                its, ite, jts, jte, kts, kte )


   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,    INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

   CHARACTER(LEN=1) ,                          INTENT(IN   ) :: name

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ,                      &
                                               INTENT(IN   ) :: field,    &
                                                                alt

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: mut

   REAL , DIMENSION( kms:kme ) ,                   INTENT(IN   ) :: rdn, rdnw

   REAL ,                                      INTENT(IN   ) :: kvdif
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end

   REAL , DIMENSION(its:ite, jts:jte) :: vfluxm, vfluxp, zz
   REAL , DIMENSION(its:ite, 0:kte+1) :: vflux

   REAL :: rdz

   LOGICAL :: specified

!<DESCRIPTION>
!
!  vertical_diffusion
!  computes vertical diffusion tendency.
!
!</DESCRIPTION>

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   ktf=MIN(kte,kde-1)
   
   IF (name .EQ. 'w')THEN

   
   i_start = its
   i_end   = MIN(ite,ide-1)
   j_start = jts
   j_end   = MIN(jte,jde-1)

j_loop_w : DO j = j_start, j_end

     DO k=kts,ktf-1
       DO i = i_start, i_end
          vflux(i,k)= (kvdif/alt(i,k,j))*rdnw(k)*(field(i,k+1,j)-field(i,k,j))
       ENDDO
     ENDDO

     DO i = i_start, i_end
       vflux(i,ktf)=0.
     ENDDO

     DO k=kts+1,ktf
       DO i = i_start, i_end
            tendency(i,k,j)=tendency(i,k,j)                                         &
                              +rdn(k)*g*g/mut(i,j)/(0.5*(alt(i,k,j)+alt(i,k-1,j)))  &
                                         *(vflux(i,k)-vflux(i,k-1))
       ENDDO
     ENDDO

    ENDDO j_loop_w

   ELSE IF(name .EQ. 'm')THEN

     i_start = its
     i_end   = MIN(ite,ide-1)
     j_start = jts
     j_end   = MIN(jte,jde-1)

j_loop_s : DO j = j_start, j_end

     DO k=kts,ktf-1
       DO i = i_start, i_end
         vflux(i,k)=kvdif*rdn(k+1)/(0.5*(alt(i,k,j)+alt(i,k+1,j)))   &
                  *(field(i,k+1,j)-field(i,k,j))
       ENDDO
     ENDDO

     DO i = i_start, i_end
       vflux(i,0)=vflux(i,1)
     ENDDO

     DO i = i_start, i_end
       vflux(i,ktf)=0.
     ENDDO

     DO k=kts,ktf
       DO i = i_start, i_end
         tendency(i,k,j)=tendency(i,k,j)+g*g/mut(i,j)/alt(i,k,j)  &
                *rdnw(k)*(vflux(i,k)-vflux(i,k-1))
       ENDDO
     ENDDO

 ENDDO j_loop_s

   ENDIF

END SUBROUTINE vertical_diffusion


!-------------------------------------------------------------------------------

SUBROUTINE vertical_diffusion_mp ( field, tendency, config_flags, &
                                   base,                          &
                                   alt, mut, rdn, rdnw, kvdif,    &
                                   ids, ide, jds, jde, kds, kde,  &
                                   ims, ime, jms, jme, kms, kme,  &
                                   its, ite, jts, jte, kts, kte  )


   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,    INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ,                      &
                                               INTENT(IN   ) :: field,    &
                                                                alt

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: mut

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: rdn,  &
                                                                  rdnw, &
                                                                  base

   REAL ,                                      INTENT(IN   ) :: kvdif
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end

   REAL , DIMENSION(its:ite, 0:kte+1) :: vflux

   REAL :: rdz

   LOGICAL :: specified

!<DESCRIPTION>
!
!  vertical_diffusion_mp
!  computes vertical diffusion tendency of a perturbation variable
!  (field-base).  Note that base as a 1D (k) field.
!
!</DESCRIPTION>

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   ktf=MIN(kte,kde-1)
   
     i_start = its
     i_end   = MIN(ite,ide-1)
     j_start = jts
     j_end   = MIN(jte,jde-1)

j_loop_s : DO j = j_start, j_end

     DO k=kts,ktf-1
       DO i = i_start, i_end
         vflux(i,k)=kvdif*rdn(k+1)/(0.5*(alt(i,k,j)+alt(i,k+1,j)))   &
                    *(field(i,k+1,j)-field(i,k,j)-base(k+1)+base(k))
       ENDDO
     ENDDO

     DO i = i_start, i_end
       vflux(i,0)=vflux(i,1)
     ENDDO

     DO i = i_start, i_end
       vflux(i,ktf)=0.
     ENDDO

     DO k=kts,ktf
       DO i = i_start, i_end
         tendency(i,k,j)=tendency(i,k,j)+g*g/mut(i,j)/alt(i,k,j)  &
                *rdnw(k)*(vflux(i,k)-vflux(i,k-1))
       ENDDO
     ENDDO

 ENDDO j_loop_s

END SUBROUTINE vertical_diffusion_mp


!-------------------------------------------------------------------------------

SUBROUTINE vertical_diffusion_3dmp ( field, tendency, config_flags, &
                                     base_3d,                       &
                                     alt, mut, rdn, rdnw, kvdif,    &
                                     ids, ide, jds, jde, kds, kde,  &
                                     ims, ime, jms, jme, kms, kme,  &
                                     its, ite, jts, jte, kts, kte  )


   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,    INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ,                      &
                                               INTENT(IN   ) :: field,    &
                                                                alt,      &
                                                                base_3d

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: mut

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: rdn,  &
                                                                  rdnw

   REAL ,                                      INTENT(IN   ) :: kvdif
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end

   REAL , DIMENSION(its:ite, 0:kte+1) :: vflux

   REAL :: rdz

   LOGICAL :: specified

!<DESCRIPTION>
!
!  vertical_diffusion_3dmp
!  computes vertical diffusion tendency of a perturbation variable
!  (field-base_3d).  
!
!</DESCRIPTION>

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   ktf=MIN(kte,kde-1)
   
     i_start = its
     i_end   = MIN(ite,ide-1)
     j_start = jts
     j_end   = MIN(jte,jde-1)

j_loop_s : DO j = j_start, j_end

     DO k=kts,ktf-1
       DO i = i_start, i_end
         vflux(i,k)=kvdif*rdn(k+1)/(0.5*(alt(i,k,j)+alt(i,k+1,j)))   &
                    *(   field(i,k+1,j)  -field(i,k,j)               &
                      -base_3d(i,k+1,j)+base_3d(i,k,j) )
       ENDDO
     ENDDO

     DO i = i_start, i_end
       vflux(i,0)=vflux(i,1)
     ENDDO

     DO i = i_start, i_end
       vflux(i,ktf)=0.
     ENDDO

     DO k=kts,ktf
       DO i = i_start, i_end
         tendency(i,k,j)=tendency(i,k,j)+g*g/mut(i,j)/alt(i,k,j)  &
                *rdnw(k)*(vflux(i,k)-vflux(i,k-1))
       ENDDO
     ENDDO

 ENDDO j_loop_s

END SUBROUTINE vertical_diffusion_3dmp


!-------------------------------------------------------------------------------


SUBROUTINE vertical_diffusion_u ( field, tendency,              &
                                  config_flags, u_base,         &
                                  alt, muu, rdn, rdnw, kvdif,   &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte )


   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,    INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ,                      &
                                               INTENT(IN   ) :: field,    &
                                                                alt

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: muu

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: rdn, rdnw, u_base

   REAL ,                                      INTENT(IN   ) :: kvdif
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end

   REAL , DIMENSION(its:ite, 0:kte+1) :: vflux

   REAL :: rdz, zz

   LOGICAL :: specified

!<DESCRIPTION>
!
!  vertical_diffusion_u computes vertical diffusion tendency for 
!  the u momentum equation.  This routine assumes a constant eddy
!  viscosity kvdif.
!
!</DESCRIPTION>

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   ktf=MIN(kte,kde-1)

      i_start = its
      i_end   = ite
      j_start = jts
      j_end   = MIN(jte,jde-1)

      IF ( config_flags%open_xs .or. specified ) i_start = MAX(ids+1,its)
      IF ( config_flags%open_xe .or. specified ) i_end   = MIN(ide-1,ite)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite


j_loop_u : DO j = j_start, j_end

     DO k=kts,ktf-1
       DO i = i_start, i_end
         vflux(i,k)=kvdif*rdn(k+1)/(0.25*( alt(i  ,k  ,j)      &
                                        +alt(i-1,k  ,j)      &
                                        +alt(i  ,k+1,j)      &
                                        +alt(i-1,k+1,j) ) )  &
                             *(field(i,k+1,j)-field(i,k,j)   &
                               -u_base(k+1)   +u_base(k)  )
       ENDDO
     ENDDO

     DO i = i_start, i_end
       vflux(i,0)=vflux(i,1)
     ENDDO

     DO i = i_start, i_end
       vflux(i,ktf)=0.
     ENDDO

     DO k=kts,ktf-1
       DO i = i_start, i_end
         tendency(i,k,j)=tendency(i,k,j)+                             &
                g*g*rdnw(k)/muu(i,j)/(0.5*(alt(i-1,k,j)+alt(i,k,j)))* &
                              (vflux(i,k)-vflux(i,k-1))
       ENDDO
     ENDDO

 ENDDO j_loop_u
   
END SUBROUTINE vertical_diffusion_u

!-------------------------------------------------------------------------------


SUBROUTINE vertical_diffusion_v ( field, tendency,              &
                                  config_flags, v_base,         &
                                  alt, muv, rdn, rdnw, kvdif,   &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme, &
                                  its, ite, jts, jte, kts, kte )


   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type), INTENT(IN   ) :: config_flags

   INTEGER ,    INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                 ims, ime, jms, jme, kms, kme, &
                                 its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ,                      &
                                               INTENT(IN   ) :: field,    &
                                                                alt
   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: rdn, rdnw, v_base

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: muv

   REAL ,                                      INTENT(IN   ) :: kvdif
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf, jm1
   INTEGER :: i_start, i_end, j_start, j_end

   REAL , DIMENSION(its:ite, 0:kte+1) :: vflux

   REAL :: rdz, zz

   LOGICAL :: specified

!<DESCRIPTION>
!
!  vertical_diffusion_v computes vertical diffusion tendency for 
!  the v momentum equation.  This routine assumes a constant eddy
!  viscosity kvdif.
!
!</DESCRIPTION>

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   ktf=MIN(kte,kde-1)
   
      i_start = its
      i_end   = MIN(ite,ide-1)
      j_start = jts
      j_end   = MIN(jte,jde-1)

      IF ( config_flags%open_ys .or. specified ) j_start = MAX(jds+1,jts)
      IF ( config_flags%open_ye .or. specified ) j_end   = MIN(jde-1,jte)

j_loop_v : DO j = j_start, j_end
!     jm1 = max(j-1,1)
     jm1 = j-1

     DO k=kts,ktf-1
       DO i = i_start, i_end
         vflux(i,k)=kvdif*rdn(k+1)/(0.25*( alt(i,k  ,j  )      &
                                        +alt(i,k  ,jm1)      &
                                        +alt(i,k+1,j  )      &
                                        +alt(i,k+1,jm1) ) )  &
                             *(field(i,k+1,j)-field(i,k,j)   &
                               -v_base(k+1)   +v_base(k)  )
       ENDDO
     ENDDO

     DO i = i_start, i_end
       vflux(i,0)=vflux(i,1)
     ENDDO

     DO i = i_start, i_end
       vflux(i,ktf)=0.
     ENDDO

     DO k=kts,ktf-1
       DO i = i_start, i_end 
         tendency(i,k,j)=tendency(i,k,j)+                              &
                g*g*rdnw(k)/muv(i,j)/(0.5*(alt(i,k,jm1)+alt(i,k,j)))*  &
                              (vflux(i,k)-vflux(i,k-1))
       ENDDO
     ENDDO

 ENDDO j_loop_v
   
END SUBROUTINE vertical_diffusion_v

!***************  end new mass coordinate routines

!-------------------------------------------------------------------------------

SUBROUTINE calculate_full ( rfield, rfieldb, rfieldp,     &
                            ids, ide, jds, jde, kds, kde, &
                            ims, ime, jms, jme, kms, kme, &
                            its, ite, jts, jte, kts, kte )

   IMPLICIT NONE
   
   ! Input data
   
   INTEGER ,      INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                   ims, ime, jms, jme, kms, kme, &
                                   its, ite, jts, jte, kts, kte 
   
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: rfieldb, &
                                                                      rfieldp

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT  ) :: rfield
   
   ! Local indices.
   
   INTEGER :: i, j, k, itf, jtf, ktf
   
!<DESCRIPTION>
!
!  calculate_full
!  calculates full 3D field from pertubation and base field.
!
!</DESCRIPTION>

   itf=MIN(ite,ide-1)
   jtf=MIN(jte,jde-1)
   ktf=MIN(kte,kde-1)

   DO j=jts,jtf
   DO k=kts,ktf
   DO i=its,itf
      rfield(i,k,j)=rfieldb(i,k,j)+rfieldp(i,k,j)
   ENDDO
   ENDDO
   ENDDO

END SUBROUTINE calculate_full

!------------------------------------------------------------------------------

SUBROUTINE coriolis ( ru, rv, rw, ru_tend, rv_tend, rw_tend, &
                      config_flags,                          &
                      f, e, sina, cosa, fzm, fzp,            &
                      ids, ide, jds, jde, kds, kde,          &
                      ims, ime, jms, jme, kms, kme,          &
                      its, ite, jts, jte, kts, kte          )

   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type) ,           INTENT(IN   ) :: config_flags   

   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: ru_tend, &
                                                                rv_tend, &
                                                                rw_tend
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: ru, &
                                                                rv, &
                                                                rw

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: f,    &
                                                                    e,    &
                                                                    sina, &
                                                                    cosa

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm, &
                                                                  fzp
   
   ! Local indices.
   
   INTEGER :: i, j , k, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   
   LOGICAL :: specified

!<DESCRIPTION>
!
!  coriolis calculates the large timestep tendency terms in the 
!  u, v, and w momentum equations arise from the coriolis force.
!
!</DESCRIPTION>

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   ktf=MIN(kte,kde-1)

! coriolis for u-momentum equation

   i_start = its
   i_end   = ite
   IF ( config_flags%open_xs .or. specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. specified .or. &
        config_flags%nested) i_end   = MIN(ide-1,ite)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite

   DO j = jts, MIN(jte,jde-1)

   DO k=kts,ktf
   DO i = i_start, i_end
   
     ru_tend(i,k,j)=ru_tend(i,k,j) + 0.5*(f(i,j)+f(i-1,j)) &
       *0.25*(rv(i-1,k,j+1)+rv(i,k,j+1)+rv(i-1,k,j)+rv(i,k,j)) &
           - 0.5*(e(i,j)+e(i-1,j))*0.5*(cosa(i,j)+cosa(i-1,j)) &
       *0.25*(rw(i-1,k+1,j)+rw(i-1,k,j)+rw(i,k+1,j)+rw(i,k,j))

   ENDDO
   ENDDO

   IF ( (config_flags%open_xs) .and. (its == ids) ) THEN

     DO k=kts,ktf
   
       ru_tend(its,k,j)=ru_tend(its,k,j) + 0.5*(f(its,j)+f(its,j))   &
         *0.25*(rv(its,k,j+1)+rv(its,k,j+1)+rv(its,k,j)+rv(its,k,j)) &
             - 0.5*(e(its,j)+e(its,j))*0.5*(cosa(its,j)+cosa(its,j)) &
         *0.25*(rw(its,k+1,j)+rw(its,k,j)+rw(its,k+1,j)+rw(its,k,j))

     ENDDO

   ENDIF

   IF ( (config_flags%open_xe) .and. (ite == ide) ) THEN

     DO k=kts,ktf
   
       ru_tend(ite,k,j)=ru_tend(ite,k,j) + 0.5*(f(ite-1,j)+f(ite-1,j)) &
         *0.25*(rv(ite-1,k,j+1)+rv(ite-1,k,j+1)+rv(ite-1,k,j)+rv(ite-1,k,j)) &
             - 0.5*(e(ite-1,j)+e(ite-1,j))*0.5*(cosa(ite-1,j)+cosa(ite-1,j)) &
         *0.25*(rw(ite-1,k+1,j)+rw(ite-1,k,j)+rw(ite-1,k+1,j)+rw(ite-1,k,j))

     ENDDO

   ENDIF

   ENDDO

!  coriolis term for v-momentum equation

   j_start = jts
   j_end   = jte

   IF ( config_flags%open_ys .or. specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. specified .or. &
        config_flags%nested) j_end   = MIN(jde-1,jte)

   IF ( (config_flags%open_ys) .and. (jts == jds) ) THEN

     DO k=kts,ktf
     DO i=its,MIN(ide-1,ite)
   
        rv_tend(i,k,jts)=rv_tend(i,k,jts) - 0.5*(f(i,jts)+f(i,jts))    &
         *0.25*(ru(i,k,jts)+ru(i+1,k,jts)+ru(i,k,jts)+ru(i+1,k,jts))   &
             + 0.5*(e(i,jts)+e(i,jts))*0.5*(sina(i,jts)+sina(i,jts))   &
             *0.25*(rw(i,k+1,jts)+rw(i,k,jts)+rw(i,k+1,jts)+rw(i,k,jts)) 

     ENDDO
     ENDDO

   ENDIF

   DO j=j_start, j_end
   DO k=kts,ktf
   DO i=its,MIN(ide-1,ite)
   
      rv_tend(i,k,j)=rv_tend(i,k,j) - 0.5*(f(i,j)+f(i,j-1))    &
       *0.25*(ru(i,k,j)+ru(i+1,k,j)+ru(i,k,j-1)+ru(i+1,k,j-1)) &
           + 0.5*(e(i,j)+e(i,j-1))*0.5*(sina(i,j)+sina(i,j-1)) &
           *0.25*(rw(i,k+1,j-1)+rw(i,k,j-1)+rw(i,k+1,j)+rw(i,k,j)) 

   ENDDO
   ENDDO
   ENDDO


   IF ( (config_flags%open_ye) .and. (jte == jde) ) THEN

     DO k=kts,ktf
     DO i=its,MIN(ide-1,ite)
   
        rv_tend(i,k,jte)=rv_tend(i,k,jte) - 0.5*(f(i,jte-1)+f(i,jte-1))        &
         *0.25*(ru(i,k,jte-1)+ru(i+1,k,jte-1)+ru(i,k,jte-1)+ru(i+1,k,jte-1))   &
             + 0.5*(e(i,jte-1)+e(i,jte-1))*0.5*(sina(i,jte-1)+sina(i,jte-1))   &
             *0.25*(rw(i,k+1,jte-1)+rw(i,k,jte-1)+rw(i,k+1,jte-1)+rw(i,k,jte-1)) 

     ENDDO
     ENDDO

   ENDIF

! coriolis term for w-mometum 

   DO j=jts,MIN(jte, jde-1)
   DO k=kts+1,ktf
   DO i=its,MIN(ite, ide-1)

       rw_tend(i,k,j)=rw_tend(i,k,j) + e(i,j)*           &
          (cosa(i,j)*0.5*(fzm(k)*(ru(i,k,j)+ru(i+1,k,j)) &
          +fzp(k)*(ru(i,k-1,j)+ru(i+1,k-1,j)))           &
          -sina(i,j)*0.5*(fzm(k)*(rv(i,k,j)+rv(i,k,j+1)) & 
          +fzp(k)*(rv(i,k-1,j)+rv(i,k-1,j+1))))

   ENDDO
   ENDDO
   ENDDO

END SUBROUTINE coriolis

!------------------------------------------------------------------------------

SUBROUTINE perturbation_coriolis ( ru_in, rv_in, rw, ru_tend, rv_tend, rw_tend, &
                                   config_flags,                                &
                                   u_base, v_base, z_base,                      &
                                   muu, muv, phb, ph,                           &
                                   f, e, sina, cosa, fzm, fzp,                  &
                                   ids, ide, jds, jde, kds, kde,                &
                                   ims, ime, jms, jme, kms, kme,                &
                                   its, ite, jts, jte, kts, kte                )

   IMPLICIT NONE
   
   ! Input data
   
   TYPE(grid_config_rec_type) ,           INTENT(IN   ) :: config_flags   

   INTEGER ,                 INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                              ims, ime, jms, jme, kms, kme, &
                                              its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: ru_tend, &
                                                                rv_tend, &
                                                                rw_tend
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: ru_in, &
                                                                      rv_in, &
                                                                      rw,    &
                                                                      ph,    &
                                                                      phb


   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: f,    &
                                                                    e,    &
                                                                    sina, &
                                                                    cosa

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: muu, &
                                                                    muv
                                                                    

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm, &
                                                                  fzp

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: u_base,  &
                                                                  v_base,  &
                                                                  z_base
   
   ! Local storage

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) :: ru, &
                                                      rv

   REAL  :: z_at_u, z_at_v, wkp1, wk, wkm1

   ! Local indices.
   
   INTEGER :: i, j , k, ktf
   INTEGER :: i_start, i_end, j_start, j_end
   
   LOGICAL :: specified

!<DESCRIPTION>
!
!  perturbation_coriolis calculates the large timestep tendency terms in the 
!  u, v, and w momentum equations arise from the coriolis force.  This version
!  subtracts off the horizontal velocities from the initial sounding when
!  computing the forcing terms, hence "perturbation" coriolis.
!
!</DESCRIPTION>

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

   ktf=MIN(kte,kde-1)

! coriolis for u-momentum equation

   i_start = its
   i_end   = ite
   IF ( config_flags%open_xs .or. specified .or. &
        config_flags%nested) i_start = MAX(ids+1,its)
   IF ( config_flags%open_xe .or. specified .or. &
        config_flags%nested) i_end   = MIN(ide-1,ite)
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite

!  compute perturbation mu*v for use in u momentum equation

   DO j = jts, MIN(jte,jde-1)+1
   DO k=kts+1,ktf-1
   DO i = i_start-1, i_end
     z_at_v = 0.25*( phb(i,k,j  )+phb(i,k+1,j  )  &
                    +phb(i,k,j-1)+phb(i,k+1,j-1)  &
                    +ph(i,k,j  )+ph(i,k+1,j  )    &
                    +ph(i,k,j-1)+ph(i,k+1,j-1))/g
     wkp1 = min(1.,max(0.,z_at_v-z_base(k))/(z_base(k+1)-z_base(k)))
     wkm1 = min(1.,max(0.,z_base(k)-z_at_v)/(z_base(k)-z_base(k-1)))
     wk   = 1.-wkp1-wkm1
     rv(i,k,j) = rv_in(i,k,j) - muv(i,j)*(            &
                                  wkm1*v_base(k-1)    &
                                 +wk  *v_base(k  )    &
                                 +wkp1*v_base(k+1)   )
   ENDDO
   ENDDO
   ENDDO


!  pick up top and bottom v 

   DO j = jts, MIN(jte,jde-1)+1
   DO i = i_start-1, i_end

     k = kts
     z_at_v = 0.25*( phb(i,k,j  )+phb(i,k+1,j  )  &
                    +phb(i,k,j-1)+phb(i,k+1,j-1)  &
                    +ph(i,k,j  )+ph(i,k+1,j  )    &
                    +ph(i,k,j-1)+ph(i,k+1,j-1))/g
     wkp1 = min(1.,max(0.,z_at_v-z_base(k))/(z_base(k+1)-z_base(k)))
     wk   = 1.-wkp1
     rv(i,k,j) = rv_in(i,k,j) - muv(i,j)*(            &
                                 +wk  *v_base(k  )    &
                                 +wkp1*v_base(k+1)   )

     k = ktf
     z_at_v = 0.25*( phb(i,k,j  )+phb(i,k+1,j  )  &
                    +phb(i,k,j-1)+phb(i,k+1,j-1)  &
                    +ph(i,k,j  )+ph(i,k+1,j  )    &
                    +ph(i,k,j-1)+ph(i,k+1,j-1))/g
     wkm1 = min(1.,max(0.,z_base(k)-z_at_v)/(z_base(k)-z_base(k-1)))
     wk   = 1.-wkm1
     rv(i,k,j) = rv_in(i,k,j) - muv(i,j)*(            &
                                  wkm1*v_base(k-1)    &
                                 +wk  *v_base(k  )   )

   ENDDO
   ENDDO

!  compute coriolis forcing for u

   DO j = jts, MIN(jte,jde-1)

   DO k=kts,ktf
     DO i = i_start, i_end
       ru_tend(i,k,j)=ru_tend(i,k,j) + 0.5*(f(i,j)+f(i-1,j)) &
         *0.25*(rv(i-1,k,j+1)+rv(i,k,j+1)+rv(i-1,k,j)+rv(i,k,j)) &
             - 0.5*(e(i,j)+e(i-1,j))*0.5*(cosa(i,j)+cosa(i-1,j)) &
         *0.25*(rw(i-1,k+1,j)+rw(i-1,k,j)+rw(i,k+1,j)+rw(i,k,j))
     ENDDO
   ENDDO

   IF ( (config_flags%open_xs) .and. (its == ids) ) THEN

     DO k=kts,ktf
   
       ru_tend(its,k,j)=ru_tend(its,k,j) + 0.5*(f(its,j)+f(its,j))   &
         *0.25*(rv(its,k,j+1)+rv(its,k,j+1)+rv(its,k,j)+rv(its,k,j)) &
             - 0.5*(e(its,j)+e(its,j))*0.5*(cosa(its,j)+cosa(its,j)) &
         *0.25*(rw(its,k+1,j)+rw(its,k,j)+rw(its,k+1,j)+rw(its,k,j))

     ENDDO

   ENDIF

   IF ( (config_flags%open_xe) .and. (ite == ide) ) THEN

     DO k=kts,ktf
   
       ru_tend(ite,k,j)=ru_tend(ite,k,j) + 0.5*(f(ite-1,j)+f(ite-1,j)) &
         *0.25*(rv(ite-1,k,j+1)+rv(ite-1,k,j+1)+rv(ite-1,k,j)+rv(ite-1,k,j)) &
             - 0.5*(e(ite-1,j)+e(ite-1,j))*0.5*(cosa(ite-1,j)+cosa(ite-1,j)) &
         *0.25*(rw(ite-1,k+1,j)+rw(ite-1,k,j)+rw(ite-1,k+1,j)+rw(ite-1,k,j))

     ENDDO

   ENDIF

   ENDDO

!  coriolis term for v-momentum equation

   j_start = jts
   j_end   = jte

   IF ( config_flags%open_ys .or. specified .or. &
        config_flags%nested) j_start = MAX(jds+1,jts)
   IF ( config_flags%open_ye .or. specified .or. &
        config_flags%nested) j_end   = MIN(jde-1,jte)

!  compute perturbation mu*u for use in v momentum equation

   DO j = j_start-1,j_end
   DO k=kts+1,ktf-1
   DO i = its, MIN(ite,ide-1)+1
     z_at_u = 0.25*( phb(i  ,k,j)+phb(i  ,k+1,j)  &
                    +phb(i-1,k,j)+phb(i-1,k+1,j)  &
                    +ph(i  ,k,j)+ph(i  ,k+1,j)    &
                    +ph(i-1,k,j)+ph(i-1,k+1,j))/g
     wkp1 = min(1.,max(0.,z_at_u-z_base(k))/(z_base(k+1)-z_base(k)))
     wkm1 = min(1.,max(0.,z_base(k)-z_at_u)/(z_base(k)-z_base(k-1)))
     wk   = 1.-wkp1-wkm1
     ru(i,k,j) = ru_in(i,k,j) - muu(i,j)*(            &
                                  wkm1*u_base(k-1)    &
                                 +wk  *u_base(k  )    &
                                 +wkp1*u_base(k+1)   )
   ENDDO
   ENDDO
   ENDDO

!  pick up top and bottom u

   DO j = j_start-1,j_end
   DO i = its, MIN(ite,ide-1)+1

     k = kts
     z_at_u = 0.25*( phb(i  ,k,j)+phb(i  ,k+1,j)  &
                    +phb(i-1,k,j)+phb(i-1,k+1,j)  &
                    +ph(i  ,k,j)+ph(i  ,k+1,j)    &
                    +ph(i-1,k,j)+ph(i-1,k+1,j))/g
     wkp1 = min(1.,max(0.,z_at_u-z_base(k))/(z_base(k+1)-z_base(k)))
     wk   = 1.-wkp1
     ru(i,k,j) = ru_in(i,k,j) - muu(i,j)*(            &
                                 +wk  *u_base(k  )    &
                                 +wkp1*u_base(k+1)   )


     k = ktf
     z_at_u = 0.25*( phb(i  ,k,j)+phb(i  ,k+1,j)  &
                    +phb(i-1,k,j)+phb(i-1,k+1,j)  &
                    +ph(i  ,k,j)+ph(i  ,k+1,j)    &
                    +ph(i-1,k,j)+ph(i-1,k+1,j))/g
     wkm1 = min(1.,max(0.,z_base(k)-z_at_u)/(z_base(k)-z_base(k-1)))
     wk   = 1.-wkm1
     ru(i,k,j) = ru_in(i,k,j) - muu(i,j)*(            &
                                  wkm1*u_base(k-1)    &
                                 +wk  *u_base(k  )   )

   ENDDO
   ENDDO

!  compute coriolis forcing for v momentum equation

   IF ( (config_flags%open_ys) .and. (jts == jds) ) THEN

     DO k=kts,ktf
     DO i=its,MIN(ide-1,ite)
   
        rv_tend(i,k,jts)=rv_tend(i,k,jts) - 0.5*(f(i,jts)+f(i,jts))    &
         *0.25*(ru(i,k,jts)+ru(i+1,k,jts)+ru(i,k,jts)+ru(i+1,k,jts))   &
             + 0.5*(e(i,jts)+e(i,jts))*0.5*(sina(i,jts)+sina(i,jts))   &
             *0.25*(rw(i,k+1,jts)+rw(i,k,jts)+rw(i,k+1,jts)+rw(i,k,jts)) 

     ENDDO
     ENDDO

   ENDIF

   DO j=j_start, j_end
   DO k=kts,ktf
   DO i=its,MIN(ide-1,ite)
   
      rv_tend(i,k,j)=rv_tend(i,k,j) - 0.5*(f(i,j)+f(i,j-1))    &
       *0.25*(ru(i,k,j)+ru(i+1,k,j)+ru(i,k,j-1)+ru(i+1,k,j-1)) &
           + 0.5*(e(i,j)+e(i,j-1))*0.5*(sina(i,j)+sina(i,j-1)) &
           *0.25*(rw(i,k+1,j-1)+rw(i,k,j-1)+rw(i,k+1,j)+rw(i,k,j)) 

   ENDDO
   ENDDO
   ENDDO


   IF ( (config_flags%open_ye) .and. (jte == jde) ) THEN

     DO k=kts,ktf
     DO i=its,MIN(ide-1,ite)
   
        rv_tend(i,k,jte)=rv_tend(i,k,jte) - 0.5*(f(i,jte-1)+f(i,jte-1))        &
         *0.25*(ru(i,k,jte-1)+ru(i+1,k,jte-1)+ru(i,k,jte-1)+ru(i+1,k,jte-1))   &
             + 0.5*(e(i,jte-1)+e(i,jte-1))*0.5*(sina(i,jte-1)+sina(i,jte-1))   &
             *0.25*(rw(i,k+1,jte-1)+rw(i,k,jte-1)+rw(i,k+1,jte-1)+rw(i,k,jte-1)) 

     ENDDO
     ENDDO

   ENDIF

! coriolis term for w-mometum 

   DO j=jts,MIN(jte, jde-1)
   DO k=kts+1,ktf
   DO i=its,MIN(ite, ide-1)

       rw_tend(i,k,j)=rw_tend(i,k,j) + e(i,j)*           &
          (cosa(i,j)*0.5*(fzm(k)*(ru(i,k,j)+ru(i+1,k,j)) &
          +fzp(k)*(ru(i,k-1,j)+ru(i+1,k-1,j)))           &
          -sina(i,j)*0.5*(fzm(k)*(rv(i,k,j)+rv(i,k,j+1)) & 
          +fzp(k)*(rv(i,k-1,j)+rv(i,k-1,j+1))))

   ENDDO
   ENDDO
   ENDDO

END SUBROUTINE perturbation_coriolis

!------------------------------------------------------------------------------

SUBROUTINE curvature ( ru, rv, rw, u, v, w, ru_tend, rv_tend, rw_tend, &
                        config_flags,                                       &
                        msfu, msfv, fzm, fzp, rdx, rdy,                 &
                        ids, ide, jds, jde, kds, kde,                   &
                        ims, ime, jms, jme, kms, kme,                   &
                        its, ite, jts, jte, kts, kte                   )


   IMPLICIT NONE
   
   ! Input data

   TYPE(grid_config_rec_type) ,           INTENT(IN   ) :: config_flags   

   INTEGER ,                  INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                               ims, ime, jms, jme, kms, kme, &
                                               its, ite, jts, jte, kts, kte
   
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ,                     &
                                               INTENT(INOUT) :: ru_tend, &
                                                                rv_tend, &
                                                                rw_tend

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ,                     &
                                               INTENT(IN   ) :: ru,      &
                                                                rv,      &
                                                                rw,      &
                                                                u,       &
                                                                v,       &
                                                                w

   REAL , DIMENSION( ims:ime , jms:jme ) ,         INTENT(IN   ) :: msfu,    &
                                                                msfv

   REAL , DIMENSION( kms:kme ) ,                 INTENT(IN   ) :: fzm,     &
                                                                fzp

   REAL ,                                      INTENT(IN   ) :: rdx,     &
                                                                rdy
   
   ! Local data
   
!   INTEGER :: i, j, k, itf, jtf, ktf, kp1, im, ip, jm, jp
   INTEGER :: i, j, k, itf, jtf, ktf
   INTEGER :: i_start, i_end, j_start, j_end
!   INTEGER :: irmin, irmax, jrmin, jrmax

   REAL , DIMENSION( its-1:ite , kts:kte, jts-1:jte ) :: vxgm

   LOGICAL :: specified

!<DESCRIPTION>
!
!  curvature calculates the large timestep tendency terms in the 
!  u, v, and w momentum equations arise from the curvature terms.  
!
!</DESCRIPTION>

   specified = .false.
   if(config_flags%specified .or. config_flags%nested) specified = .true.

      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
      ktf=MIN(kte,kde-1)

!   irmin = ims
!   irmax = ime
!   jrmin = jms
!   jrmax = jme
!   IF ( config_flags%open_xs ) irmin = ids
!   IF ( config_flags%open_xe ) irmax = ide-1
!   IF ( config_flags%open_ys ) jrmin = jds
!   IF ( config_flags%open_ye ) jrmax = jde-1
   
! Define v cross grad m at scalar points - vxgm(i,j)

   i_start = its-1
   i_end   = ite
   j_start = jts-1
   j_end   = jte

   IF ( ( config_flags%open_xs .or. specified .or. &
        config_flags%nested) .and. (its == ids) ) i_start = its
   IF ( ( config_flags%open_xe .or. specified .or. &
        config_flags%nested) .and. (ite == ide) ) i_end   = ite-1
   IF ( ( config_flags%open_ys .or. specified .or. &
        config_flags%nested) .and. (jts == jds) ) j_start = jts
   IF ( ( config_flags%open_ye .or. specified .or. &
        config_flags%nested) .and. (jte == jde) ) j_end   = jte-1
      IF ( config_flags%periodic_x ) i_start = its-1
      IF ( config_flags%periodic_x ) i_end = ite

   DO j=j_start, j_end
   DO k=kts,ktf
   DO i=i_start, i_end
      vxgm(i,k,j)=0.5*(u(i,k,j)+u(i+1,k,j))*(msfv(i,j+1)-msfv(i,j))*rdy - &
                  0.5*(v(i,k,j)+v(i,k,j+1))*(msfu(i+1,j)-msfu(i,j))*rdx
   ENDDO
   ENDDO
   ENDDO

!  Pick up the boundary rows for open (radiation) lateral b.c.
!  Rather crude at present, we are assuming there is no
!    variation in this term at the boundary.

   IF ( ( config_flags%open_xs .or. (specified .AND. .NOT. config_flags%periodic_x) .or. &
        config_flags%nested) .and. (its == ids) ) THEN

     DO j = jts-1, jte
     DO k = kts, ktf
       vxgm(its-1,k,j) =  vxgm(its,k,j)
     ENDDO
     ENDDO

   ENDIF

   IF ( ( config_flags%open_xe .or. (specified .AND. .NOT. config_flags%periodic_x) .or. &
        config_flags%nested) .and. (ite == ide) ) THEN

     DO j = jts-1, jte
     DO k = kts, ktf
       vxgm(ite,k,j) =  vxgm(ite-1,k,j)
     ENDDO
     ENDDO

   ENDIF

   IF ( ( config_flags%open_ys .or. specified .or. &
        config_flags%nested) .and. (jts == jds) ) THEN

     DO k = kts, ktf
     DO i = its-1, ite
       vxgm(i,k,jts-1) =  vxgm(i,k,jts)
     ENDDO
     ENDDO

   ENDIF

   IF ( ( config_flags%open_ye .or. specified .or. &
        config_flags%nested) .and. (jte == jde) ) THEN

     DO k = kts, ktf
     DO i = its-1, ite
       vxgm(i,k,jte) =  vxgm(i,k,jte-1)
     ENDDO
     ENDDO

   ENDIF

!  curvature term for u momentum eqn.

   i_start = its
   IF ( config_flags%open_xs .or. specified .or. &
        config_flags%nested) i_start = MAX ( ids+1 , its )
   IF ( config_flags%open_xe .or. specified .or. &
        config_flags%nested) i_end   = MIN ( ide-1 , ite )
      IF ( config_flags%periodic_x ) i_start = its
      IF ( config_flags%periodic_x ) i_end = ite

   DO j=jts,MIN(jde-1,jte)
   DO k=kts,ktf
   DO i=i_start,i_end

      ru_tend(i,k,j)=ru_tend(i,k,j) + 0.5*(vxgm(i,k,j)+vxgm(i-1,k,j)) &
              *0.25*(rv(i-1,k,j+1)+rv(i,k,j+1)+rv(i-1,k,j)+rv(i,k,j)) &
               - u(i,k,j)*reradius &
              *0.25*(rw(i-1,k+1,j)+rw(i-1,k,j)+rw(i,k+1,j)+rw(i,k,j))

   ENDDO
   ENDDO
   ENDDO

!  curvature term for v momentum eqn.

   j_start = jts
   IF ( config_flags%open_ys .or. specified .or. &
        config_flags%nested) j_start = MAX ( jds+1 , jts )
   IF ( config_flags%open_ye .or. specified .or. &
        config_flags%nested) j_end   = MIN ( jde-1 , jte ) 

   DO j=j_start,j_end
   DO k=kts,ktf
   DO i=its,MIN(ite,ide-1)

      rv_tend(i,k,j)=rv_tend(i,k,j) - 0.5*(vxgm(i,k,j)+vxgm(i,k,j-1)) &
              *0.25*(ru(i,k,j)+ru(i+1,k,j)+ru(i,k,j-1)+ru(i+1,k,j-1)) &
                    + v(i,k,j)*reradius                               &
              *0.25*(rw(i,k+1,j-1)+rw(i,k,j-1)+rw(i,k+1,j)+rw(i,k,j))

   ENDDO
   ENDDO
   ENDDO

!  curvature term for vertical momentum eqn.

   DO j=jts,MIN(jte,jde-1)
   DO k=MAX(2,kts),ktf
   DO i=its,MIN(ite,ide-1)

      rw_tend(i,k,j)=rw_tend(i,k,j) + reradius*                              &
    (0.5*(fzm(k)*(ru(i,k,j)+ru(i+1,k,j))+fzp(k)*(ru(i,k-1,j)+ru(i+1,k-1,j))) &
    *0.5*(fzm(k)*( u(i,k,j) +u(i+1,k,j))+fzp(k)*( u(i,k-1,j) +u(i+1,k-1,j)))     &
    +0.5*(fzm(k)*(rv(i,k,j)+rv(i,k,j+1))+fzp(k)*(rv(i,k-1,j)+rv(i,k-1,j+1))) &
    *0.5*(fzm(k)*( v(i,k,j) +v(i,k,j+1))+fzp(k)*( v(i,k-1,j) +v(i,k-1,j+1))))

   ENDDO
   ENDDO
   ENDDO

END SUBROUTINE curvature

!------------------------------------------------------------------------------

SUBROUTINE decouple ( rr, rfield, field, name, config_flags, &
                      fzm, fzp,                          &
                      ids, ide, jds, jde, kds, kde,      &
                      ims, ime, jms, jme, kms, kme,      &
                      its, ite, jts, jte, kts, kte      )

   IMPLICIT NONE

   ! Input data

   TYPE(grid_config_rec_type) ,           INTENT(IN   ) :: config_flags   

   INTEGER ,                                   INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                                                ims, ime, jms, jme, kms, kme, &
                                                                its, ite, jts, jte, kts, kte

   CHARACTER(LEN=1) ,                          INTENT(IN   ) :: name

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: rfield

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN   ) :: rr
   
   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT  ) :: field
   
   REAL , DIMENSION( kms:kme ) , INTENT(IN   ) :: fzm, fzp
   
   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf
   
!<DESCRIPTION>
!
!  decouple decouples a variable from the column dry-air mass.
!
!</DESCRIPTION>

   ktf=MIN(kte,kde-1)
   
   IF (name .EQ. 'u')THEN
      itf=ite
      jtf=MIN(jte,jde-1)

      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
         field(i,k,j)=rfield(i,k,j)/(0.5*(rr(i,k,j)+rr(i-1,k,j)))
      ENDDO
      ENDDO
      ENDDO

   ELSE IF (name .EQ. 'v')THEN
      itf=MIN(ite,ide-1)
      jtf=jte

      DO j=jts,jtf
      DO k=kts,ktf
        DO i=its,itf
             field(i,k,j)=rfield(i,k,j)/(0.5*(rr(i,k,j)+rr(i,k,j-1)))
        ENDDO
      ENDDO
      ENDDO

   ELSE IF (name .EQ. 'w')THEN
      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
      DO j=jts,jtf
      DO k=kts+1,ktf
      DO i=its,itf
         field(i,k,j)=rfield(i,k,j)/(fzm(k)*rr(i,k,j)+fzp(k)*rr(i,k-1,j))
      ENDDO
      ENDDO
      ENDDO

      DO j=jts,jtf
      DO i=its,itf
        field(i,kte,j) = 0.
      ENDDO
      ENDDO

   ELSE 
      itf=MIN(ite,ide-1)
      jtf=MIN(jte,jde-1)
   ! For theta we will decouple tb and tp and add them to give t afterwards
      DO j=jts,jtf
      DO k=kts,ktf
      DO i=its,itf
         field(i,k,j)=rfield(i,k,j)/rr(i,k,j)
      ENDDO
      ENDDO
      ENDDO
   
   ENDIF

END SUBROUTINE decouple

!-------------------------------------------------------------------------------


SUBROUTINE zero_tend ( tendency,                     &
                       ids, ide, jds, jde, kds, kde, &
                       ims, ime, jms, jme, kms, kme, &
                       its, ite, jts, jte, kts, kte )


   IMPLICIT NONE
   
   ! Input data
   
   INTEGER ,                                   INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                                                                ims, ime, jms, jme, kms, kme, &
                                                                its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(INOUT) :: tendency

   ! Local data
   
   INTEGER :: i, j, k, itf, jtf, ktf

!<DESCRIPTION>
!
!  zero_tend sets the input tendency array to zero.
!
!</DESCRIPTION>

      DO j = jts, jte
      DO k = kts, kte
      DO i = its, ite
        tendency(i,k,j) = 0.
      ENDDO
      ENDDO
      ENDDO

      END SUBROUTINE zero_tend

!======================================================================
!   physics prep routines
!======================================================================

   SUBROUTINE phy_prep ( config_flags,                                &  ! input
                         mu, muu, muv, u, v, p, pb, alt, ph,          &  ! input
                         phb, t, tsk, moist, n_moist,                 &  ! input
                         mu_3d, rho, th_phy, p_phy , pi_phy ,         &  ! output
                         u_phy, v_phy, p8w, t_phy, t8w,               &  ! output
                         z, z_at_w, dz8w,                             &  ! output
                         fzm, fzp,                                    &  ! params
                         RTHRATEN,                                    &
                         RTHBLTEN, RUBLTEN, RVBLTEN,                  &
                         RQVBLTEN, RQCBLTEN, RQIBLTEN,                &
                         RTHCUTEN, RQVCUTEN, RQCCUTEN,                &
                         RQRCUTEN, RQICUTEN, RQSCUTEN,                &
                         RTHFTEN,  RQVFTEN,                           &
                         RUNDGDTEN, RVNDGDTEN, RTHNDGDTEN,            &
                         RQVNDGDTEN, RMUNDGDTEN,                      &
                         ids, ide, jds, jde, kds, kde,                &
                         ims, ime, jms, jme, kms, kme,                &
                         its, ite, jts, jte, kts, kte                )
!----------------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------------

   TYPE(grid_config_rec_type) ,     INTENT(IN   ) :: config_flags

   INTEGER ,        INTENT(IN   ) ::   ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       its, ite, jts, jte, kts, kte
   INTEGER ,          INTENT(IN   ) :: n_moist

   REAL, DIMENSION( ims:ime, kms:kme , jms:jme , n_moist ), INTENT(IN) :: moist


   REAL , DIMENSION( ims:ime, jms:jme ), INTENT(IN   )   ::     TSK, mu, muu, muv

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ,                 &
          INTENT(  OUT)                                  ::   u_phy, &
                                                              v_phy, &
                                                             pi_phy, &
                                                              p_phy, &
                                                                p8w, &
                                                              t_phy, &
                                                             th_phy, &
                                                                t8w, &
                                                              mu_3d, &
                                                                rho, &
                                                                  z, &
                                                               dz8w, &
                                                              z_at_w 

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) ,                 &
          INTENT(IN   )                                  ::      pb, &
                                                                  p, &
                                                                  u, &
                                                                  v, &
                                                                alt, &
                                                                 ph, &
                                                                phb, &
                                                                  t


   REAL , DIMENSION( kms:kme ) ,           INTENT(IN   ) ::     fzm,   &
                                                                fzp

   REAL,  DIMENSION( ims:ime , kms:kme, jms:jme ),                   &
          INTENT(INOUT)   ::                               RTHRATEN  

   REAL,  DIMENSION( ims:ime , kms:kme, jms:jme ),                   &
          INTENT(INOUT)   ::                               RTHCUTEN, &
                                                           RQVCUTEN, &
                                                           RQCCUTEN, &
                                                           RQRCUTEN, &
                                                           RQICUTEN, &
                                                           RQSCUTEN

   REAL,  DIMENSION( ims:ime, kms:kme, jms:jme )                   , &
          INTENT(INOUT)   ::                                RUBLTEN, &
                                                            RVBLTEN, &
                                                           RTHBLTEN, &
                                                           RQVBLTEN, &
                                                           RQCBLTEN, &
                                                           RQIBLTEN

   REAL,  DIMENSION( ims:ime, kms:kme, jms:jme )                   , &
          INTENT(INOUT)   ::                                RTHFTEN, &
                                                            RQVFTEN

   REAL,  DIMENSION( ims:ime, kms:kme, jms:jme )                   , &
          INTENT(INOUT)   ::                                RUNDGDTEN, &
                                                            RVNDGDTEN, &
                                                           RTHNDGDTEN, &
                                                           RQVNDGDTEN, &
                                                           RMUNDGDTEN

   INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end, i_startu, j_startv
   INTEGER :: i, j, k
   REAL    :: w1, w2, z0, z1, z2

!-----------------------------------------------------------------------

!<DESCRIPTION>
!
!  phys_prep calculates a number of diagnostic quantities needed by
!  the physics routines.  It also decouples the physics tendencies from
!  the column dry-air mass (the physics routines expect to see/update the
!  uncoupled tendencies).
!
!</DESCRIPTION>

!  set up loop bounds for this grids boundary conditions

    i_start = its
    i_end   = min( ite,ide-1 )
    j_start = jts
    j_end   = min( jte,jde-1 )

    k_start = kts
    k_end = min( kte, kde-1 )

!  compute thermodynamics and velocities at pressure points

    do j = j_start,j_end
    do k = k_start, k_end
    do i = i_start, i_end

      th_phy(i,k,j) = t(i,k,j) + t0
      p_phy(i,k,j) = p(i,k,j) + pb(i,k,j)
      pi_phy(i,k,j) = (p_phy(i,k,j)/p1000mb)**rcp
      t_phy(i,k,j) = th_phy(i,k,j)*pi_phy(i,k,j)
      rho(i,k,j) = 1./alt(i,k,j)*(1.+moist(i,k,j,P_QV))
      mu_3d(i,k,j) = mu(i,j)
      u_phy(i,k,j) = 0.5*(u(i,k,j)+u(i+1,k,j))
      v_phy(i,k,j) = 0.5*(v(i,k,j)+v(i,k,j+1))

    enddo
    enddo
    enddo

!  compute z at w points

    do j = j_start,j_end
    do k = k_start, kte
    do i = i_start, i_end
      z_at_w(i,k,j) = (phb(i,k,j)+ph(i,k,j))/g
    enddo
    enddo
    enddo

    do j = j_start,j_end
    do k = k_start, kte-1
    do i = i_start, i_end
      dz8w(i,k,j) = z_at_w(i,k+1,j)-z_at_w(i,k,j)
    enddo
    enddo
    enddo

    do j = j_start,j_end
    do i = i_start, i_end
      dz8w(i,kte,j) = 0.
    enddo
    enddo

!  compute z at p points (average of z at w points)

    do j = j_start,j_end
    do k = k_start, k_end
    do i = i_start, i_end
      z(i,k,j) = 0.5*(z_at_w(i,k,j) +z_at_w(i,k+1,j) )
    enddo
    enddo
    enddo

!  interp t and p at w points

    do j = j_start,j_end
    do k = 2, k_end
    do i = i_start, i_end
      p8w(i,k,j) = fzm(k)*p_phy(i,k,j)+fzp(k)*p_phy(i,k-1,j)
      t8w(i,k,j) = fzm(k)*t_phy(i,k,j)+fzp(k)*t_phy(i,k-1,j)
    enddo
    enddo
    enddo

!  extrapolate p and t to surface and top.
!  well use an extrapolation in z for now

    do j = j_start,j_end
    do i = i_start, i_end

! bottom

      z0 = z_at_w(i,1,j)
      z1 = z(i,1,j)
      z2 = z(i,2,j)
      w1 = (z0 - z2)/(z1 - z2)
      w2 = 1. - w1
      p8w(i,1,j) = w1*p_phy(i,1,j)+w2*p_phy(i,2,j)
      t8w(i,1,j) = w1*t_phy(i,1,j)+w2*t_phy(i,2,j)

! top

      z0 = z_at_w(i,kte,j)
      z1 = z(i,k_end,j)
      z2 = z(i,k_end-1,j)
      w1 = (z0 - z2)/(z1 - z2)
      w2 = 1. - w1

!      p8w(i,kde,j) = w1*p_phy(i,kde-1,j)+w2*p_phy(i,kde-2,j)
!!!  bug fix      extrapolate ln(p) so p is positive definite
      p8w(i,kde,j) = exp(w1*log(p_phy(i,kde-1,j))+w2*log(p_phy(i,kde-2,j)))
      t8w(i,kde,j) = w1*t_phy(i,kde-1,j)+w2*t_phy(i,kde-2,j)

    enddo
    enddo

! decouple all physics tendencies

   IF (config_flags%ra_lw_physics .gt. 0 .or. config_flags%ra_sw_physics .gt. 0) THEN

      DO J=j_start,j_end
      DO K=k_start,k_end
      DO I=i_start,i_end
         RTHRATEN(I,K,J)=RTHRATEN(I,K,J)/mu(I,J)
      ENDDO
      ENDDO
      ENDDO

   ENDIF

   IF (config_flags%cu_physics .gt. 0) THEN

      DO J=j_start,j_end
      DO I=i_start,i_end
      DO K=k_start,k_end
         RTHCUTEN(I,K,J)=RTHCUTEN(I,K,J)/mu(I,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QV .ge. PARAM_FIRST_SCALAR)THEN
         DO J=j_start,j_end
         DO I=i_start,i_end
         DO K=k_start,k_end
            RQVCUTEN(I,K,J)=RQVCUTEN(I,K,J)/mu(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QC .ge. PARAM_FIRST_SCALAR)THEN
         DO J=j_start,j_end
         DO I=i_start,i_end
         DO K=k_start,k_end
            RQCCUTEN(I,K,J)=RQCCUTEN(I,K,J)/mu(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QR .ge. PARAM_FIRST_SCALAR)THEN
         DO J=j_start,j_end
         DO I=i_start,i_end
         DO K=k_start,k_end
            RQRCUTEN(I,K,J)=RQRCUTEN(I,K,J)/mu(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR)THEN
         DO J=j_start,j_end
         DO I=i_start,i_end
         DO K=k_start,k_end
            RQICUTEN(I,K,J)=RQICUTEN(I,K,J)/mu(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF(P_QS .ge. PARAM_FIRST_SCALAR)THEN
         DO J=j_start,j_end
         DO I=i_start,i_end
         DO K=k_start,k_end
            RQSCUTEN(I,K,J)=RQSCUTEN(I,K,J)/mu(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

   ENDIF

   IF (config_flags%bl_pbl_physics .gt. 0) THEN

      DO J=j_start,j_end
      DO K=k_start,k_end
      DO I=i_start,i_end
         RUBLTEN(I,K,J) =RUBLTEN(I,K,J)/mu(I,J)
         RVBLTEN(I,K,J) =RVBLTEN(I,K,J)/mu(I,J)
         RTHBLTEN(I,K,J)=RTHBLTEN(I,K,J)/mu(I,J)
      ENDDO
      ENDDO
      ENDDO

      IF (P_QV .ge. PARAM_FIRST_SCALAR) THEN
         DO J=j_start,j_end
         DO K=k_start,k_end
         DO I=i_start,i_end
            RQVBLTEN(I,K,J)=RQVBLTEN(I,K,J)/mu(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QC .ge. PARAM_FIRST_SCALAR) THEN
         DO J=j_start,j_end
         DO K=k_start,k_end
         DO I=i_start,i_end
           RQCBLTEN(I,K,J)=RQCBLTEN(I,K,J)/mu(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (P_QI .ge. PARAM_FIRST_SCALAR) THEN
         DO J=j_start,j_end
         DO K=k_start,k_end
         DO I=i_start,i_end
            RQIBLTEN(I,K,J)=RQIBLTEN(I,K,J)/mu(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

    ENDIF

!  decouple advective forcing required by Grell-Devenyi scheme

   if ( config_flags%cu_physics == GDSCHEME ) then

      DO J=j_start,j_end
      DO I=i_start,i_end
         DO K=k_start,k_end
            RTHFTEN(I,K,J)=RTHFTEN(I,K,J)/mu(I,J)
         ENDDO
      ENDDO
      ENDDO

      IF (P_QV .ge. PARAM_FIRST_SCALAR)THEN
         DO J=j_start,j_end
         DO I=i_start,i_end
            DO K=k_start,k_end
               RQVFTEN(I,K,J)=RQVFTEN(I,K,J)/mu(I,J)
            ENDDO
         ENDDO
         ENDDO
      ENDIF

   END IF

! fdda
! note fdda u and v tendencies are staggered, also only interior points have muu/muv,
!   so only decouple those

   IF (config_flags%grid_fdda .gt. 0) THEN

      i_startu=MAX(its,ids+1)
      j_startv=MAX(jts,jds+1)

      DO J=j_start,j_end
      DO K=k_start,k_end
      DO I=i_startu,i_end
         RUNDGDTEN(I,K,J) =RUNDGDTEN(I,K,J)/muu(I,J)
      ENDDO
      ENDDO
      ENDDO
      DO J=j_startv,j_end
      DO K=k_start,k_end
      DO I=i_start,i_end
         RVNDGDTEN(I,K,J) =RVNDGDTEN(I,K,J)/muv(I,J)
      ENDDO
      ENDDO
      ENDDO
      DO J=j_start,j_end
      DO K=k_start,k_end
      DO I=i_start,i_end
         RTHNDGDTEN(I,K,J)=RTHNDGDTEN(I,K,J)/mu(I,J)
!        RMUNDGDTEN(I,J) - no coupling
      ENDDO
      ENDDO
      ENDDO
      IF (P_QV .ge. PARAM_FIRST_SCALAR) THEN
         DO J=j_start,j_end
         DO K=k_start,k_end
         DO I=i_start,i_end
            RQVNDGDTEN(I,K,J)=RQVNDGDTEN(I,K,J)/mu(I,J)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

    ENDIF

END SUBROUTINE phy_prep

!------------------------------------------------------------

   SUBROUTINE moist_physics_prep_em( t_new, t_old, t0, rho, al, alb, &
                                     p, p8w, p0, pb, ph, phb,        &
                                     th_phy, pii, pf,                &
                                     z, z_at_w, dz8w,                &
                                     dt,h_diabatic,                  &
                                     config_flags,fzm, fzp,          &
                                     ids,ide, jds,jde, kds,kde,      &
                                     ims,ime, jms,jme, kms,kme,      &
                                     its,ite, jts,jte, kts,kte      )

   IMPLICIT NONE

! Here we construct full fields
! needed by the microphysics

   TYPE(grid_config_rec_type),    INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
   INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
   INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

   REAL, INTENT(IN   )  ::  dt

   REAL, DIMENSION( ims:ime , kms:kme, jms:jme ),        &
         INTENT(IN   ) ::                           al,  &
                                                    alb, &
                                                    p,   &
                                                    pb,  &
                                                    ph,  &
                                                    phb


   REAL , DIMENSION( kms:kme ) ,           INTENT(IN   ) ::   fzm, &
                                                              fzp

   REAL, DIMENSION( ims:ime , kms:kme, jms:jme ),       &
         INTENT(  OUT) ::                         rho,  &
                                               th_phy,  &
                                                  pii,  &
                                                  pf,   &
                                                    z,  &
                                               z_at_w,  &
                                                 dz8w,  &
                                                  p8w
! pjj/cray
!                                                 p8w,  &
!                                          h_diabatic

   REAL, DIMENSION( ims:ime , kms:kme, jms:jme ),       &
         INTENT(INOUT) ::                         h_diabatic

   REAL, DIMENSION( ims:ime , kms:kme, jms:jme ),        &
         INTENT(INOUT) ::                         t_new, &
                                                  t_old

   REAL, INTENT(IN   ) :: t0, p0
   REAL                :: z0,z1,z2,w1,w2

   INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
   INTEGER :: i, j, k

!--------------------------------------------------------------------

!<DESCRIPTION>
!
!  moist_phys_prep_em calculates a number of diagnostic quantities needed by
!  the microphysics routines.
!
!</DESCRIPTION>

!  set up loop bounds for this grids boundary conditions

    i_start = its    
    i_end   = min( ite,ide-1 )
    j_start = jts    
    j_end   = min( jte,jde-1 )

    k_start = kts
    k_end = min( kte, kde-1 )

     DO j = j_start, j_end
     DO k = k_start, kte
     DO i = i_start, i_end
       z_at_w(i,k,j) = (ph(i,k,j)+phb(i,k,j))/g
     ENDDO
     ENDDO
     ENDDO

    do j = j_start,j_end
    do k = k_start, kte-1
    do i = i_start, i_end
      dz8w(i,k,j) = z_at_w(i,k+1,j)-z_at_w(i,k,j)
    enddo
    enddo
    enddo

    do j = j_start,j_end
    do i = i_start, i_end
      dz8w(i,kte,j) = 0.
    enddo
    enddo


           !  compute full pii, rho, and z at the new time-level
           !  (needed for physics).
           !  convert perturbation theta to full theta (th_phy)
           !  use h_diabatic to temporarily save pre-microphysics full theta

     DO j = j_start, j_end
     DO k = k_start, k_end
     DO i = i_start, i_end

       th_phy(i,k,j) = t_new(i,k,j) + t0
       h_diabatic(i,k,j) = th_phy(i,k,j)
       rho(i,k,j)  = 1./(al(i,k,j)+alb(i,k,j))
       pii(i,k,j) = ((p(i,k,j)+pb(i,k,j))/p0)**rcp
       z(i,k,j) = 0.5*(z_at_w(i,k,j) +z_at_w(i,k+1,j) )
       pf(i,k,j) = p(i,k,j)+pb(i,k,j)

     ENDDO
     ENDDO
     ENDDO

!  interp t and p at w points

    do j = j_start,j_end
    do k = 2, k_end
    do i = i_start, i_end
      p8w(i,k,j) = fzm(k)*pf(i,k,j)+fzp(k)*pf(i,k-1,j)
    enddo
    enddo
    enddo

!  extrapolate p and t to surface and top.
!  well use an extrapolation in z for now

    do j = j_start,j_end
    do i = i_start, i_end

! bottom

      z0 = z_at_w(i,1,j)
      z1 = z(i,1,j)
      z2 = z(i,2,j)
      w1 = (z0 - z2)/(z1 - z2)
      w2 = 1. - w1
      p8w(i,1,j) = w1*pf(i,1,j)+w2*pf(i,2,j)

! top

      z0 = z_at_w(i,kte,j)
      z1 = z(i,k_end,j)
      z2 = z(i,k_end-1,j)
      w1 = (z0 - z2)/(z1 - z2)
      w2 = 1. - w1
!      p8w(i,kde,j) = w1*pf(i,kde-1,j)+w2*pf(i,kde-2,j)
      p8w(i,kde,j) = exp(w1*log(pf(i,kde-1,j))+w2*log(pf(i,kde-2,j)))

    enddo
    enddo

   END SUBROUTINE moist_physics_prep_em

!------------------------------------------------------------------------------

   SUBROUTINE moist_physics_finish_em( t_new, t_old, t0, mut,     &
                                       th_phy, h_diabatic, dt,    &
                                       config_flags,              &
                                       ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte )

   IMPLICIT NONE

! Here we construct full fields
! needed by the microphysics

   TYPE(grid_config_rec_type),    INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
   INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
   INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime , kms:kme, jms:jme ),        &
         INTENT(INOUT) ::                         t_new, &
                                                  t_old, &
                                                 th_phy, &
                                                  h_diabatic

   REAL, DIMENSION( ims:ime , jms:jme ),  INTENT(INOUT) ::  mut


   REAL, INTENT(IN   ) :: t0, dt

   INTEGER :: i_start, i_end, j_start, j_end, k_start, k_end
   INTEGER :: i, j, k

!--------------------------------------------------------------------

!<DESCRIPTION>
!
!  moist_phys_finish_em resets theta to its perturbation value and
!  computes and stores the microphysics diabatic heating term.
!
!</DESCRIPTION>

!  set up loop bounds for this grids boundary conditions


    i_start = its    
    i_end   = min( ite,ide-1 )
    j_start = jts    
    j_end   = min( jte,jde-1 )

    k_start = kts
    k_end = min( kte, kde-1 )

!  add microphysics theta diff to perturbation theta, set h_diabatic

     DO j = j_start, j_end
     DO k = k_start, k_end
     DO i = i_start, i_end

       t_new(i,k,j) = t_new(i,k,j) + (th_phy(i,k,j)-h_diabatic(i,k,j))
       h_diabatic(i,k,j) = (th_phy(i,k,j)-h_diabatic(i,k,j))/dt
!       h_diabatic(i,k,j) = 0.

     ENDDO
     ENDDO
     ENDDO

   END SUBROUTINE moist_physics_finish_em

!----------------------------------------------------------------


   SUBROUTINE init_module_big_step
   END SUBROUTINE init_module_big_step

SUBROUTINE set_tend ( field, field_adv_tend, msf,       &
                      ids, ide, jds, jde, kds, kde,     &
                      ims, ime, jms, jme, kms, kme,     &
                      its, ite, jts, jte, kts, kte       )

   IMPLICIT NONE

   ! Input data

   INTEGER ,  INTENT(IN   ) :: ids, ide, jds, jde, kds, kde, &
                               ims, ime, jms, jme, kms, kme, &
                               its, ite, jts, jte, kts, kte

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(OUT) :: field

   REAL , DIMENSION( ims:ime , kms:kme , jms:jme ) , INTENT(IN)  :: field_adv_tend

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN)  :: msf

   ! Local data

   INTEGER :: i, j, k, itf, jtf, ktf

!<DESCRIPTION>
!
!  set_tend copies the advective tendency array into the tendency array.
!
!</DESCRIPTION>

      jtf = MIN(jte,jde-1)
      ktf = MIN(kte,kde-1)
      itf = MIN(ite,ide-1)
      DO j = jts, jtf
      DO k = kts, ktf
      DO i = its, itf
         field(i,k,j) = field_adv_tend(i,k,j)*msf(i,j)
      ENDDO
      ENDDO
      ENDDO

END SUBROUTINE set_tend

!------------------------------------------------------------------------------

    SUBROUTINE rk_rayleigh_damp( ru_tendf, rv_tendf,              &
                                 rw_tendf, t_tendf,               &
                                 u, v, w, t, t_init,              &
                                 mut, muu, muv, ph, phb,          &
                                 u_base, v_base, t_base, z_base,  &
                                 dampcoef, zdamp,                 &
                                 ids, ide, jds, jde, kds, kde,    &
                                 ims, ime, jms, jme, kms, kme,    &
                                 its, ite, jts, jte, kts, kte   )

! History:     Apr 2005  Modifications by George Bryan, NCAR:
!                  - Generalized the code in a way that allows for
!                    simulations with steep terrain.
!
!              Jul 2004  Modifications by George Bryan, NCAR:
!                  - Modified the code to use u_base, v_base, and t_base
!                    arrays for the background state.  Removed the hard-wired
!                    base-state values.
!                  - Modified the code to use dampcoef, zdamp, and damp_opt,
!                    i.e., the upper-level damper variables in namelist.input.
!                    Removed the hard-wired variables in the older version.
!                    This damper is used when damp_opt = 2.
!                  - Modified the code to account for the movement of the
!                    model surfaces with time.  The code now obtains a base-
!                    state value by interpolation using the "_base" arrays.

!              Nov 2003  Bug fix by Jason Knievel, NCAR

!              Aug 2003  Meridional dimension, some comments, and
!                        changes in layout of the code added by
!                        Jason Knievel, NCAR

!              Jul 2003  Original code by Bill Skamarock, NCAR

! Purpose:     This routine applies Rayleigh damping to a layer at top
!              of the model domain.

!-----------------------------------------------------------------------
! Begin declarations.

    IMPLICIT NONE

    INTEGER, INTENT( IN )  &
    :: ids, ide, jds, jde, kds, kde,  &
       ims, ime, jms, jme, kms, kme,  &
       its, ite, jts, jte, kts, kte

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( INOUT )  &
    :: ru_tendf, rv_tendf, rw_tendf, t_tendf

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )  &
    :: u, v, w, t, t_init, ph, phb

    REAL, DIMENSION( ims:ime, jms:jme ),  INTENT( IN )  &
    :: mut, muu, muv

    REAL, DIMENSION( kms:kme ) ,  INTENT(IN   )  &
    :: u_base, v_base, t_base, z_base

    REAL, INTENT(IN   )   &
    :: dampcoef, zdamp

! Local variables.

    INTEGER  &
    :: i_start, i_end, j_start, j_end, k_start, k_end, i, j, k, ktf, k1, k2

    REAL  &
    :: pii, dcoef, z, ztop

    REAL :: wkp1, wk, wkm1

    REAL, DIMENSION( kms:kme ) :: z00, u00, v00, t00

! End declarations.
!-----------------------------------------------------------------------

    pii = 2.0 * asin(1.0)

    ktf = MIN( kte,   kde-1 )

!-----------------------------------------------------------------------
! Adjust u to base state.

    DO j = jts, MIN( jte, jde-1 )
    DO i = its, MIN( ite, ide   )

      ! Get height at top of model
      ztop = 0.5*( phb(i  ,kde,j)+phb(i-1,kde,j)   &
                  +ph(i  ,kde,j)+ph(i-1,kde,j) )/g

      ! Find bottom of damping layer
      k1 = ktf
      z = ztop
      DO WHILE( z >= (ztop-zdamp) )
        z = 0.25*( phb(i  ,k1,j)+phb(i  ,k1+1,j)  &
                  +phb(i-1,k1,j)+phb(i-1,k1+1,j)  &
                  +ph(i  ,k1,j)+ph(i  ,k1+1,j)    &
                  +ph(i-1,k1,j)+ph(i-1,k1+1,j))/g
        z00(k1) = z
        k1 = k1 - 1
      ENDDO
      k1 = k1 + 2

      ! Get reference state at model levels
      DO k = k1, ktf
        k2 = ktf
        DO WHILE( z_base(k2) .gt. z00(k) )
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt.ktf)then
          u00(k) = u_base(k2) + ( u_base(k2) - u_base(k2-1) )   &
                              * (     z00(k) - z_base(k2)   )   &
                              / ( z_base(k2) - z_base(k2-1) )
        else
          u00(k) = u_base(k2) + ( u_base(k2+1) - u_base(k2) )   &
                              * (       z00(k) - z_base(k2) )   &
                              / ( z_base(k2+1) - z_base(k2) )
        endif
      ENDDO

      ! Apply the Rayleigh damper
      DO k = k1, ktf
        dcoef = 1.0 - MIN( 1.0, ( ztop - z00(k) ) / zdamp )
        dcoef = (SIN( 0.5 * pii * dcoef ) )**2
        ru_tendf(i,k,j) = ru_tendf(i,k,j) -                    &
                          muu(i,j) * ( dcoef * dampcoef ) *    &
                          ( u(i,k,j) - u00(k) )
      END DO

    END DO
    END DO

! End adjustment of u.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Adjust v to base state.

    DO j = jts, MIN( jte, jde   )
    DO i = its, MIN( ite, ide-1 )

      ! Get height at top of model
      ztop = 0.5*( phb(i,kde,j  )+phb(i,kde,j-1)   &
                  +ph(i,kde,j  )+ph(i,kde,j-1) )/g

      ! Find bottom of damping layer
      k1 = ktf
      z = ztop
      DO WHILE( z >= (ztop-zdamp) )
        z = 0.25*( phb(i,k1,j  )+phb(i,k1+1,j  )  &
                  +phb(i,k1,j-1)+phb(i,k1+1,j-1)  &
                  +ph(i,k1,j  )+ph(i,k1+1,j  )    &
                  +ph(i,k1,j-1)+ph(i,k1+1,j-1))/g
        z00(k1) = z
        k1 = k1 - 1
      ENDDO
      k1 = k1 + 2

      ! Get reference state at model levels
      DO k = k1, ktf
        k2 = ktf
        DO WHILE( z_base(k2) .gt. z00(k) )
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt.ktf)then
          v00(k) = v_base(k2) + ( v_base(k2) - v_base(k2-1) )   &
                              * (     z00(k) - z_base(k2)   )   &
                              / ( z_base(k2) - z_base(k2-1) )
        else
          v00(k) = v_base(k2) + ( v_base(k2+1) - v_base(k2) )   &
                              * (       z00(k) - z_base(k2) )   &
                              / ( z_base(k2+1) - z_base(k2) )
        endif
      ENDDO

      ! Apply the Rayleigh damper
      DO k = k1, ktf
        dcoef = 1.0 - MIN( 1.0, ( ztop - z00(k) ) / zdamp )
        dcoef = (SIN( 0.5 * pii * dcoef ) )**2
        rv_tendf(i,k,j) = rv_tendf(i,k,j) -                    &
                          muv(i,j) * ( dcoef * dampcoef ) *    &
                          ( v(i,k,j) - v00(k) )
      END DO

    END DO
    END DO

! End adjustment of v.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Adjust w to base state.

    DO j = jts, MIN( jte,   jde-1 )
    DO i = its, MIN( ite,   ide-1 )
      ztop = ( phb(i,kde,j) + ph(i,kde,j) ) / g
      DO k = kts, MIN( kte,   kde   )
        z = ( phb(i,k,j) + ph(i,k,j) ) / g
        IF ( z >= (ztop-zdamp) ) THEN
          dcoef = 1.0 - MIN( 1.0, ( ztop - z ) / zdamp )
          dcoef = ( SIN( 0.5 * pii * dcoef ) )**2
          rw_tendf(i,k,j) = rw_tendf(i,k,j) -  &
                            mut(i,j) * ( dcoef * dampcoef ) * w(i,k,j)
        END IF
      END DO
    END DO
    END DO

! End adjustment of w.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Adjust potential temperature to base state.

    DO j = jts, MIN( jte,   jde-1 )
    DO i = its, MIN( ite,   ide-1 )

      ! Get height at top of model
      ztop = ( phb(i,kde,j) + ph(i,kde,j) ) / g

      ! Find bottom of damping layer
      k1 = ktf
      z = ztop
      DO WHILE( z >= (ztop-zdamp) )
        z = 0.5 * ( phb(i,k1,j) + phb(i,k1+1,j) +  &
                     ph(i,k1,j) +  ph(i,k1+1,j) ) / g
        z00(k1) = z
        k1 = k1 - 1
      ENDDO
      k1 = k1 + 2

      ! Get reference state at model levels
      DO k = k1, ktf
        k2 = ktf
        DO WHILE( z_base(k2) .gt. z00(k) )
          k2 = k2 - 1
        ENDDO
        if(k2+1.gt.ktf)then
          t00(k) = t_base(k2) + ( t_base(k2) - t_base(k2-1) )   &
                              * (     z00(k) - z_base(k2)   )   &
                              / ( z_base(k2) - z_base(k2-1) )
        else
          t00(k) = t_base(k2) + ( t_base(k2+1) - t_base(k2) )   &
                              * (       z00(k) - z_base(k2) )   &
                              / ( z_base(k2+1) - z_base(k2) )
        endif
      ENDDO

      ! Apply the Rayleigh damper
      DO k = k1, ktf
        dcoef = 1.0 - MIN( 1.0, ( ztop - z00(k) ) / zdamp )
        dcoef = (SIN( 0.5 * pii * dcoef ) )**2
        t_tendf(i,k,j) = t_tendf(i,k,j) -                      &
                         mut(i,j) * ( dcoef * dampcoef )  *    &
                         ( t(i,k,j) - t00(k) )
      END DO

    END DO
    END DO

! End adjustment of potential temperature.
!-----------------------------------------------------------------------

    END SUBROUTINE rk_rayleigh_damp

!==============================================================================
!==============================================================================
                                                                                
      SUBROUTINE sixth_order_diffusion( name, field, tendency, mu, dt,  &
                                        config_flags,                   &
                                        diff_6th_opt, diff_6th_factor,  &
                                        ids, ide, jds, jde, kds, kde,   &
                                        ims, ime, jms, jme, kms, kme,   &
                                        its, ite, jts, jte, kts, kte )
                                                                                
! History:       14 Nov 2006   Name of variable changed by Jason Knievel
!                07 Jun 2006   Revised and generalized by Jason Knievel  
!                25 Apr 2005   Original code by Jason Knievel, NCAR
                                                                                
! Purpose:       Apply 6th-order, monotonic (flux-limited), numerical
!                diffusion to 3-d velocity and to scalars.
                                                                                
! References:    Ming Xue (MWR Aug 2000)
!                Durran ("Numerical Methods for Wave Equations..." 1999)
!                George Bryan (personal communication)
 
!------------------------------------------------------------------------------
! Begin: Declarations.

    IMPLICIT NONE

    INTEGER, INTENT(IN)  &
    :: ids, ide, jds, jde, kds, kde,   &
       ims, ime, jms, jme, kms, kme,   &
       its, ite, jts, jte, kts, kte
 
    TYPE(grid_config_rec_type), INTENT(IN)  &
    :: config_flags
 
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT)  &
    :: tendency
 
    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN)  &
    :: field
 
    REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN)  &
    :: mu
 
    REAL, INTENT(IN)  &
    :: dt

    REAL, INTENT(IN)  &
    :: diff_6th_factor

    INTEGER, INTENT(IN)  &
    :: diff_6th_opt

    CHARACTER(LEN=1) , INTENT(IN)  &
    :: name

    INTEGER  &
    :: i, j, k,         &
       i_start, i_end,  &
       j_start, j_end,  &
       k_start, k_end,  &
       ktf
 
    REAL  &
    :: dflux_x_p0, dflux_y_p0,  &
       dflux_x_p1, dflux_y_p1,  &
       tendency_x, tendency_y,  &
       mu_avg_p0, mu_avg_p1,    &
       diff_6th_coef

    LOGICAL  &
    :: specified
 
! End: Declarations.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin: Translate the diffusion factor into a diffusion coefficient.  See
! Durrans text, section 2.4.3, then adjust for sixth-order diffusion (not
! fourth) and for diffusion in two dimensions (not one).  For reference, a
! factor of 1.0 would mean complete diffusion of a 2dx wave in one time step,
! although application of the flux limiter reduces somewhat the effects of
! diffusion for a given coefficient.

    diff_6th_coef = diff_6th_factor * 0.015625 / ( 2.0 * dt )  

! End: Translate diffusion factor.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin: Assign limits of spatial loops depending on variable to be diffused.
! The halo regions are already filled with values by the time this subroutine
! is called, which allows the stencil to extend beyond the domains edges.

    ktf = MIN( kte, kde-1 )

    IF ( name .EQ. 'u' ) THEN

      i_start = its
      i_end   = ite
      j_start = jts
      j_end   = MIN(jde-1,jte)
      k_start = kts
      k_end   = ktf

    ELSE IF ( name .EQ. 'v' ) THEN
 
      i_start = its
      i_end   = MIN(ide-1,ite)
      j_start = jts
      j_end   = jte
      k_start = kts
      k_end   = ktf
 
    ELSE IF ( name .EQ. 'w' ) THEN

      i_start = its
      i_end   = MIN(ide-1,ite)
      j_start = jts
      j_end   = MIN(jde-1,jte)
      k_start = kts+1
      k_end   = ktf

    ELSE

      i_start = its
      i_end   = MIN(ide-1,ite)
      j_start = jts
      j_end   = MIN(jde-1,jte)
      k_start = kts
      k_end   = ktf
 
    ENDIF

! End: Assignment of limits of spatial loops.
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Begin: Loop across spatial dimensions.

    DO j = j_start, j_end
    DO k = k_start, k_end
    DO i = i_start, i_end

!------------------------------------------------------------------------------
! Begin: Diffusion in x (i index).
 
! Calculate the diffusive flux in x direction (from Xues eq. 3).
 
      dflux_x_p0 = (  10.0 * ( field(i,  k,j) - field(i-1,k,j) )    &
                     - 5.0 * ( field(i+1,k,j) - field(i-2,k,j) )    &
                     +       ( field(i+2,k,j) - field(i-3,k,j) ) )
 
      dflux_x_p1 = (  10.0 * ( field(i+1,k,j) - field(i  ,k,j) )    &
                     - 5.0 * ( field(i+2,k,j) - field(i-1,k,j) )    &
                     +       ( field(i+3,k,j) - field(i-2,k,j) ) )
 
! If requested in the namelist (diff_6th_opt=2), prohibit up-gradient diffusion
! (variation on Xues eq. 10).

      IF ( diff_6th_opt .EQ. 2 ) THEN
 
        IF ( dflux_x_p0 * ( field(i  ,k,j)-field(i-1,k,j) ) .LE. 0.0 ) THEN
          dflux_x_p0 = 0.0
        END IF
 
        IF ( dflux_x_p1 * ( field(i+1,k,j)-field(i  ,k,j) ) .LE. 0.0 ) THEN
          dflux_x_p1 = 0.0
        END IF

      END IF

! Apply 6th-order diffusion in x direction.
 
      IF      ( name .EQ. 'u' ) THEN
        mu_avg_p0 = mu(i-1,j)
        mu_avg_p1 = mu(i  ,j)
      ELSE IF ( name .EQ. 'v' ) THEN
        mu_avg_p0 = 0.25 * (       &
                    mu(i-1,j-1) +  &
                    mu(i  ,j-1) +  &
                    mu(i-1,j  ) +  &
                    mu(i  ,j  ) )
        mu_avg_p1 = 0.25 * (       &
                    mu(i  ,j-1) +  &
                    mu(i+1,j-1) +  &
                    mu(i  ,j  ) +  &
                    mu(i+1,j  ) )
      ELSE
        mu_avg_p0 = 0.5 * (        &
                    mu(i-1,j) +    &
                    mu(i  ,j) )
        mu_avg_p1 = 0.5 * (        &
                    mu(i  ,j) +    &
                    mu(i+1,j) )
      END IF
 
      tendency_x = diff_6th_coef *  &
                 ( ( mu_avg_p1 * dflux_x_p1 ) - ( mu_avg_p0 * dflux_x_p0 ) )
 
! End: Diffusion in x.
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
! Begin: Diffusion in y (j index).
 
! Calculate the diffusive flux in y direction (from Xues eq. 3).
 
      dflux_y_p0 = (  10.0 * ( field(i,k,j  ) - field(i,k,j-1) )    &
                     - 5.0 * ( field(i,k,j+1) - field(i,k,j-2) )    &
                     +       ( field(i,k,j+2) - field(i,k,j-3) ) )
 
      dflux_y_p1 = (  10.0 * ( field(i,k,j+1) - field(i,k,j  ) )    &
                     - 5.0 * ( field(i,k,j+2) - field(i,k,j-1) )    &
                     +       ( field(i,k,j+3) - field(i,k,j-2) ) )
 
! If requested in the namelist (diff_6th_opt=2), prohibit up-gradient diffusion
! (variation on Xues eq. 10).

      IF ( diff_6th_opt .EQ. 2 ) THEN
 
        IF ( dflux_y_p0 * ( field(i,k,j  )-field(i,k,j-1) ) .LE. 0.0 ) THEN
          dflux_y_p0 = 0.0
        END IF
 
        IF ( dflux_y_p1 * ( field(i,k,j+1)-field(i,k,j  ) ) .LE. 0.0 ) THEN
          dflux_y_p1 = 0.0
        END IF

      END IF
 
! Apply 6th-order diffusion in y direction.
 
      IF      ( name .EQ. 'u' ) THEN
        mu_avg_p0 = 0.25 * (       &
                    mu(i-1,j-1) +  &
                    mu(i  ,j-1) +  &
                    mu(i-1,j  ) +  &
                    mu(i  ,j  ) )
        mu_avg_p1 = 0.25 * (       &
                    mu(i-1,j  ) +  &
                    mu(i  ,j  ) +  &
                    mu(i-1,j+1) +  &
                    mu(i  ,j+1) )
      ELSE IF ( name .EQ. 'v' ) THEN
        mu_avg_p0 = mu(i,j-1)
        mu_avg_p1 = mu(i,j  )
      ELSE
        mu_avg_p0 = 0.5 * (      &
                    mu(i,j-1) +  &
                    mu(i,j  ) )
        mu_avg_p1 = 0.5 * (      &
                    mu(i,j  ) +  &
                    mu(i,j+1) )
      END IF
 
      tendency_y = diff_6th_coef *  &
                 ( ( mu_avg_p1 * dflux_y_p1 ) - ( mu_avg_p0 * dflux_y_p0 ) )
 
! End: Diffusion in y.
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
! Begin: Combine diffusion in x and y.
     
      tendency(i,k,j) = tendency(i,k,j) + tendency_x + tendency_y
 
! End: Combine diffusion in x and y.
!------------------------------------------------------------------------------

    ENDDO
    ENDDO
    ENDDO

! End: Loop across spatial dimensions.
!------------------------------------------------------------------------------
 
    END SUBROUTINE sixth_order_diffusion
 
!==============================================================================
!==============================================================================

END MODULE module_big_step_utilities_em
