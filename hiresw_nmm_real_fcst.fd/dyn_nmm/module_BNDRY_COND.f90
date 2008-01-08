!-----------------------------------------------------------------------
!
!NCEP_MESO:MODEL_LAYER: BOUNDARY CONDITION UPDATES
!
!-----------------------------------------------------------------------
!
! these define the various loop range variables
! that were defined in module_MPP. Defined as macros
! here to allow thread-safety/tile callability


! these define the various loop range variables
! that were defined in module_MPP. Defined as macros
! here to allow thread-safety/tile callability





!
!-----------------------------------------------------------------------
!
      MODULE MODULE_BNDRY_COND
!
!-----------------------------------------------------------------------
      USE MODULE_STATE_DESCRIPTION
      USE MODULE_MODEL_CONSTANTS
!-----------------------------------------------------------------------
      INCLUDE "mpif.h"
!-----------------------------------------------------------------------
      REAL :: D06666=0.06666666
!-----------------------------------------------------------------------
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE BOCOH(GRIDID,NTSD,DT0,NEST,NBC,NBOCO,LAST_TIME,TSPH    & 
     &                ,LB,ETA1,ETA2,PDTOP,PT,RES                        &
     &                ,PD_B,T_B,Q_B,U_B,V_B,Q2_B,CWM_B                  &
     &                ,PD_BT,T_BT,Q_BT,U_BT,V_BT,Q2_BT,CWM_BT           &
     &                ,PD,T,Q,Q2,CWM,PINT                               &
     &                ,MOIST,N_MOIST,SCALAR,N_SCALAR                    &
     &                ,IJDS,IJDE,SPEC_BDY_WIDTH,Z                       &  ! min/max(id,jd)
     &                ,IHE,IHW,IVE,IVW                                  &
     &                ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &                ,IMS,IME,JMS,JME,KMS,KME                          &
     &                ,ITS,ITE,JTS,JTE,KTS,KTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    BOCOH       UPDATE MASS POINTS ON BOUNDARY
!   PRGRMMR: JANJIC          ORG: W/NP22     DATE: 94-03-08
!     
! ABSTRACT:
!     TEMPERATURE, SPECIFIC HUMIDITY, AND SURFACE PRESSURE
!     ARE UPDATED ON THE DOMAIN BOUNDARY BY APPLYING THE
!     PRE-COMPUTED TENDENCIES AT EACH TIME STEP.
!     
! PROGRAM HISTORY LOG:
!   87-~??  MESINGER   - ORIGINATOR
!   95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D in HORIZONTAL
!   96-12-13  BLACK      - FINAL MODIFICATION FOR NESTED RUNS
!   98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!   00-01-06  BLACK      - MODIFIED FOR JANJIC NONHYDROSTATIC CODE
!   00-09-14  BLACK      - MODIFIED FOR DIRECT ACCESS READ
!   01-03-12  BLACK      - CONVERTED TO WRF STRUCTURE
!   02-08-29  MICHALAKES - CHANGED II=I-MY_IS_GLB+1 TO II=I
!                          ADDED CONDITIONAL COMPILATION AROUND MPI
!                          CONVERT INDEXING FROM LOCAL TO GLOBAL
!   02-09-06  WOLFE      - MORE CONVERSION TO GLOBAL INDEXING 
!   04-11-18  BLACK      - THREADED
!   05-12-19  BLACK      - CONVERTED FROM IKJ TO IJK
!   06-06-02  GOPAL      - MODIFICATIONS FOR NESTING
!     
! USAGE: CALL BOCOH FROM SUBROUTINE SOLVE_NMM
!   INPUT ARGUMENT LIST:
!
!     NOTE THAT IDE AND JDE INSIDE ROUTINE SHOULD BE PASSED IN
!     AS WHAT WRF CONSIDERS THE UNSTAGGERED GRID DIMENSIONS; THAT
!     IS, 1 LESS THAN THE IDE AND JDE SET BY WRF FRAMEWORK, JM
!  
!   OUTPUT ARGUMENT LIST: 
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!  
!     UNIQUE: NONE
!  
!     LIBRARY: NONE
!  
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM 
!$$$  
!***********************************************************************
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
      LOGICAL,INTENT(IN) :: NEST
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
      INTEGER,INTENT(IN) :: IJDS,IJDE,SPEC_BDY_WIDTH
      INTEGER,INTENT(IN) :: N_MOIST, N_SCALAR
!
      INTEGER,DIMENSION(JMS:JME),INTENT(IN) :: IHE,IHW,IVE,IVW
!
      INTEGER,INTENT(IN) :: GRIDID
      INTEGER,INTENT(IN) :: LB,NBC,NTSD
      LOGICAL,INTENT(IN) :: LAST_TIME
      INTEGER,INTENT(INOUT) :: NBOCO
!
      REAL,INTENT(IN) :: DT0,PDTOP,PT,TSPH
!
      REAL,DIMENSION(KMS:KME),INTENT(IN) :: ETA1,ETA2
!
      REAL,DIMENSION(IJDS:IJDE,1,SPEC_BDY_WIDTH,4)                      &
     &                                     ,INTENT(INOUT) :: PD_B,PD_BT
!
      REAL,DIMENSION(IJDS:IJDE,KMS:KME,SPEC_BDY_WIDTH,4)                &
     &                                ,INTENT(INOUT) :: CWM_B,Q_B,Q2_B  &
     &                                                 ,T_B,U_B,V_B 
      REAL,DIMENSION(IJDS:IJDE,KMS:KME,SPEC_BDY_WIDTH,4)                &
     &                             ,INTENT(INOUT) :: CWM_BT,Q_BT,Q2_BT  &
     &                                              ,T_BT,U_BT,V_BT 
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: RES
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: PD
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: CWM      &
     &                                                        ,PINT,Q   &
     &                                                        ,Q2,T,Z
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,NUM_MOIST)                 &
     &                                           ,INTENT(INOUT) :: MOIST
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,NUM_SCALAR)                &
     &                                          ,INTENT(INOUT) :: SCALAR
!
!
!-----------------------------------------------------------------------
!
!***  LOCAL VARIABLES
!
      INTEGER :: BF,I,IB,IBDY,II,IIM,IM,IRTN,ISIZ1,ISIZ2                &
     &          ,J,JB,JJ,JJM,JM,K,KK,N,NN,NREC,NUMGAS,NV,REC
      INTEGER :: MY_IS_GLB,MY_JS_GLB,MY_IE_GLB,MY_JE_GLB  
      INTEGER :: I_M,ILPAD1,IRPAD1,JBPAD1,JTPAD1
!
      REAL :: BCHR,CONVFAC,CWK,DT,PLYR,RRI
!
      LOGICAL :: E_BDY,W_BDY,N_BDY,S_BDY
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IM=IDE-IDS+1
      JM=JDE-JDS+1
      IIM=IM
      JJM=JM
!
      ISIZ1=2*LB
      ISIZ2=2*LB*(KME-KMS)
!
      W_BDY=(ITS==IDS)
      E_BDY=(ITE==IDE)
      S_BDY=(JTS==JDS)
      N_BDY=(JTE==JDE)
!
      ILPAD1=1
      IF(W_BDY)ILPAD1=0
      IRPAD1=1
      IF(E_BDY)IRPAD1=0
      JBPAD1=1
      IF(S_BDY)JBPAD1=0
      JTPAD1=1
      IF(N_BDY)JTPAD1=0
!
      MY_IS_GLB=ITS
      MY_IE_GLB=ITE
      MY_JS_GLB=JTS
      MY_JE_GLB=JTE
!
      DT=DT0
!
!-----------------------------------------------------------------------
!***  SOUTH AND NORTH BOUNDARIES
!-----------------------------------------------------------------------
!
!***  USE IBDY=1 FOR SOUTH; 2 FOR NORTH
!
      DO IBDY=1,2 
!
!***  MAKE SURE THE PROCESSOR HAS THIS BOUNDARY.
!
        IF((S_BDY.AND.IBDY==1).OR.(N_BDY.AND.IBDY==2))THEN
!
          IF(IBDY==1)THEN
            BF=P_YSB     ! Which boundary (YSB=the boundary where Y is at its start)
            JB=1         ! Which cell in from boundary
            JJ=1         ! Which cell in the domain
          ELSE
            BF=P_YEB     ! Which boundary (YEB=the boundary where Y is at its end)
            JB=1         ! Which cell in from boundary
            JJ=JJM       ! Which cell in the domain
          ENDIF
!
          DO I=MAX(ITS-1,IDS),MIN(ITE+1,IDE)
            PD_B(I,1,JB,BF)=PD_B(I,1,JB,BF)+PD_BT(I,1,JB,BF)*DT
            PD(I,JJ)=PD_B(I,1,JB,BF)
          ENDDO
!
!$omp parallel do                                                       &
!$omp& private(i,k)
          DO K=KTS,KTE
            DO I=MAX(ITS-1,IDS),MIN(ITE+1,IDE)
              T_B(I,K,JB,BF)=T_B(I,K,JB,BF)+T_BT(I,K,JB,BF)*DT
              Q_B(I,K,JB,BF)=Q_B(I,K,JB,BF)+Q_BT(I,K,JB,BF)*DT
              Q2_B(I,K,JB,BF)=Q2_B(I,K,JB,BF)+Q2_BT(I,K,JB,BF)*DT
              CWM_B(I,K,JB,BF)=CWM_B(I,K,JB,BF)+CWM_BT(I,K,JB,BF)*DT
              T(I,JJ,K)=T_B(I,K,JB,BF)
              Q(I,JJ,K)=Q_B(I,K,JB,BF)
              Q2(I,JJ,K)=Q2_B(I,K,JB,BF)
              CWM(I,JJ,K)=CWM_B(I,K,JB,BF)
              PINT(I,JJ,K)=ETA1(K)*PDTOP                                &
     &                    +ETA2(K)*PD(I,JJ)*RES(I,JJ)+PT
            ENDDO
          ENDDO
!
          DO I_M=1,N_MOIST
            IF(I_M==P_QV)THEN
!$omp parallel do                                                       &
!$omp& private(i,k)
              DO K=KTS,KTE
              DO I=MAX(ITS-1,IDS),MIN(ITE+1,IDE)
                MOIST(I,JJ,K,I_M)=Q(I,JJ,K)/(1.-Q(I,JJ,K))
              ENDDO
              ENDDO
            ELSE
!$omp parallel do                                                       &
!$omp& private(i,k)
              DO K=KTS,KTE
              DO I=MAX(ITS-1,IDS),MIN(ITE+1,IDE)
                MOIST(I,JJ,K,I_M)=0.
              ENDDO
              ENDDO
            ENDIF
          ENDDO
          DO I_M=2,N_SCALAR
!$omp parallel do                                                       &
!$omp& private(i,k)
            DO K=KTS,KTE
            DO I=MAX(ITS-1,IDS),MIN(ITE+1,IDE)
              SCALAR(I,JJ,K,I_M)=0.
            ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!***  WEST AND EAST BOUNDARIES
!-----------------------------------------------------------------------
!
!***  USE IBDY=1 FOR WEST; 2 FOR EAST. 
!
      DO IBDY=1,2 
!
!***  MAKE SURE THE PROCESSOR HAS THIS BOUNDARY.
!
        IF((W_BDY.AND.IBDY==1).OR.(E_BDY.AND.IBDY==2))THEN
          IF(IBDY==1)THEN
            BF=P_XSB     ! Which boundary (XSB=the boundary where X is at its start)
            IB=1         ! Which cell in from boundary 
            II=1         ! Which cell in the domain
          ELSE
            BF=P_XEB     ! Which boundary (XEB=the boundary where X is at its end)
            IB=1         ! Which cell in from boundary
            II=IIM       ! Which cell in the domain
          ENDIF
!
          DO J=MAX(JTS-1,JDS+3-1),MIN(JTE+1,JDE-2)
            IF(MOD(J,2)==1)THEN
              PD_B(J,1,IB,BF)=PD_B(J,1,IB,BF)+PD_BT(J,1,IB,BF)*DT
              PD(II,J)=PD_B(J,1,IB,BF)
            ENDIF
          ENDDO
!
!$omp parallel do                                                       &
!$omp& private(j,k)
          DO K=KTS,KTE
            DO J=MAX(JTS-1,JDS+3-1),MIN(JTE+1,JDE-2)
!
              IF(MOD(J,2)==1)THEN
                T_B(J,K,IB,BF)=T_B(J,K,IB,BF)+T_BT(J,K,IB,BF)*DT
                Q_B(J,K,IB,BF)=Q_B(J,K,IB,BF)+Q_BT(J,K,IB,BF)*DT
                Q2_B(J,K,IB,BF)=Q2_B(J,K,IB,BF)+Q2_BT(J,K,IB,BF)*DT
                CWM_B(J,K,IB,BF)=CWM_B(J,K,IB,BF)+CWM_BT(J,K,IB,BF)*DT
                T(II,J,K)=T_B(J,K,IB,BF)
                Q(II,J,K)=Q_B(J,K,IB,BF)
                Q2(II,J,K)=Q2_B(J,K,IB,BF)
                CWM(II,J,K)=CWM_B(J,K,IB,BF)
                PINT(II,J,K)=ETA1(K)*PDTOP                              &
     &                      +ETA2(K)*PD(II,J)*RES(II,J)+PT
              ENDIF
!
            ENDDO
          ENDDO
!
          DO I_M=1,N_MOIST
            IF(I_M==P_QV)THEN
!$omp parallel do                                                       &
!$omp& private(j,k)
              DO K=KTS,KTE
              DO J=MAX(JTS-1,JDS+3-1),MIN(JTE+1,JDE-2)
                IF(MOD(J,2)==1)THEN
                  MOIST(II,J,K,I_M)=Q(II,J,K)/(1.-Q(II,J,K))
                ENDIF
              ENDDO
              ENDDO
!
            ELSE
!$omp parallel do                                                       &
!$omp& private(j,k)
              DO K=KTS,KTE
              DO J=MAX(JTS-1,JDS+3-1),MIN(JTE+1,JDE-2)
                IF(MOD(J,2)==1)THEN
                  MOIST(II,J,K,I_M)=0.
                ENDIF
              ENDDO
              ENDDO
!
            ENDIF
          ENDDO
!
          DO I_M=2,N_SCALAR
!$omp parallel do                                                       &
!$omp& private(j,k)
            DO K=KTS,KTE
            DO J=MAX(JTS-1,JDS+3-1),MIN(JTE+1,JDE-2)
              IF(MOD(J,2)==1)THEN
                SCALAR(II,J,K,I_M)=0.
              ENDIF
            ENDDO
            ENDDO
          ENDDO
!
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!***  SPACE INTERPOLATION OF PD THEN REMAINING MASS VARIABLES
!***  AT INNER BOUNDARY
!-----------------------------------------------------------------------
!
!***  ONE ROW NORTH OF SOUTHERN BOUNDARY
!
      IF(S_BDY)THEN
        DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))
          PD(I,2)=0.25*(PD(I,1)+PD(I+1,1)+PD(I,3)+PD(I+1,3))
        ENDDO
      ENDIF
!
!***  ONE ROW SOUTH OF NORTHERN BOUNDARY
!
      IF(N_BDY)THEN
        DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))
          CWK=PD(I,JJM-1)
          PD(I,JJM-1)=0.25*(PD(I,JJM-2)+PD(I+1,JJM-2)                   &
     &                     +PD(I,JJM)+PD(I+1,JJM))
!
!***  NESTING TEST
!
          IF(I<=IDE-1.AND.ABS(CWK-PD(I,JJM-1))>=300.)THEN
            WRITE(0,*)'PSEUDO HYDROSTATIC IMBALANCE AT THE NORTHERN BOUNDARY AT',I,JJM-1,'GRID #',GRIDID
            WRITE(0,*)'             ',CWK/100.
            WRITE(0,*)PD(I,JJM)/100.,'               ',PD(I+1,JJM)/100.
            WRITE(0,*)'             ',PD(I,JJM-1)/100.
            WRITE(0,*)PD(I,JJM-2)/100.,'             ',PD(I+1,JJM-2)/100.
            WRITE(0,*)
          ENDIF

        ENDDO
      ENDIF
!
!***  ONE ROW EAST OF WESTERN BOUNDARY
!
      IF(W_BDY)THEN
        DO J=4,JM-3,2
!
          IF(W_BDY.AND.J>=MY_JS_GLB-JBPAD1                              &
     &            .AND.J<=MY_JE_GLB+JTPAD1)THEN
            CWK=PD(1,J)
            JJ=J
            PD(1,JJ)=0.25*(PD(1,JJ-1)+PD(2,JJ-1)+PD(1,JJ+1)+PD(2,JJ+1))
!
!***  NESTING TEST
!
             IF(ABS(CWK-PD(1,JJ))>300.)THEN
              WRITE(0,*)'PSEUDO HYDROSTATIC IMBALANCE AT THE WESTERN BOUNDARY AT',J,1,'GRID #',GRIDID
              WRITE(0,*)'             ',CWK/100.
              WRITE(0,*)PD(1,JJ+1)/100.,'               ',PD(2,JJ+1)/100.
              WRITE(0,*)'             ',PD(1,JJ)/100.
              WRITE(0,*)PD(1,JJ-1)/100,'               ',PD(2,JJ-1)/100.
              WRITE(0,*)
            ENDIF

          ENDIF
!
        ENDDO
      ENDIF
!
!***  ONE ROW WEST OF EASTERN BOUNDARY
!
      IF(E_BDY)THEN
        DO J=4,JM-3,2
!
          IF(E_BDY.AND.J>=MY_JS_GLB-JBPAD1                              &
     &            .AND.J<=MY_JE_GLB+JTPAD1)THEN
            JJ=J
            PD(IIM-1,JJ)=0.25*(PD(IIM-1,JJ-1)+PD(IIM,JJ-1)              &
     &                        +PD(IIM-1,JJ+1)+PD(IIM,JJ+1))
          ENDIF
!
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(i,j,jj,k)
      DO 200 K=KTS,KTE
!
!-----------------------------------------------------------------------
!
!***  ONE ROW NORTH OF SOUTHERN BOUNDARY
!
      IF(S_BDY)THEN
        DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))
          T(I,2,K)=(T(I,1,K)+T(I+1,1,K)+T(I,3,K)+T(I+1,3,K))*0.25
          Q(I,2,K)=(Q(I,1,K)+Q(I+1,1,K)+Q(I,3,K)+Q(I+1,3,K))*0.25
          Q2(I,2,K)=(Q2(I,1,K)+Q2(I+1,1,K)+Q2(I,3,K)+Q2(I+1,3,K))*0.25
          CWM(I,2,K)=(CWM(I,1,K)+CWM(I+1,1,K)+CWM(I,3,K)+CWM(I+1,3,K))  &
     &               *0.25
          PINT(I,2,K)=ETA1(K)*PDTOP+ETA2(K)*PD(I,2)*RES(I,2)+PT
        ENDDO
!
        DO I_M=1,N_MOIST
          IF(I_M==P_QV)THEN
            DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))
              MOIST(I,2,K,I_M)=Q(I,2,K)/(1.-Q(I,2,K))
            ENDDO
          ELSE
            DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))
              MOIST(I,2,K,I_M)=(MOIST(I,1,K,I_M)                        &
     &                         +MOIST(I+1,1,K,I_M)                      &
     &                         +MOIST(I,3,K,I_M)                        &
     &                         +MOIST(I+1,3,K,I_M))*0.25
            ENDDO
          ENDIF
        ENDDO
!
        DO I_M=2,N_SCALAR
          DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))
            SCALAR(I,2,K,I_M)=(SCALAR(I,1,K,I_M)                        &
     &                        +SCALAR(I+1,1,K,I_M)                      &
     &                        +SCALAR(I,3,K,I_M)                        &
     &                        +SCALAR(I+1,3,K,I_M))*0.25
          ENDDO
        ENDDO
!
      ENDIF
!
!***  ONE ROW SOUTH OF NORTHERN BOUNDARY
!
      IF(N_BDY)THEN
        DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))
          T(I,JJM-1,K)=(T(I,JJM-2,K)+T(I+1,JJM-2,K)                     &
     &                 +T(I,JJM,K)+T(I+1,JJM,K))                        &
     &                 *0.25
          Q(I,JJM-1,K)=(Q(I,JJM-2,K)+Q(I+1,JJM-2,K)                     &
     &                 +Q(I,JJM,K)+Q(I+1,JJM,K))                        &
     &                 *0.25
          Q2(I,JJM-1,K)=(Q2(I,JJM-2,K)+Q2(I+1,JJM-2,K)                  &
     &                  +Q2(I,JJM,K)+Q2(I+1,JJM,K))                     &
     &                  *0.25
          CWM(I,JJM-1,K)=(CWM(I,JJM-2,K)+CWM(I+1,JJM-2,K)               &
     &                   +CWM(I,JJM,K)+CWM(I+1,JJM,K))                  &
     &                   *0.25
          PINT(I,JJM-1,K)=ETA1(K)*PDTOP                                 &
     &                   +ETA2(K)*PD(I,JJM-1)*RES(I,JJM-1)+PT
        ENDDO
!
        DO I_M=1,N_MOIST
          IF(I_M==P_QV)THEN
            DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))
              MOIST(I,JJM-1,K,I_M)=Q(I,JJM-1,K)/(1.-Q(I,JJM-1,K))
            ENDDO
          ELSE
            DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))
              MOIST(I,JJM-1,K,I_M)=(MOIST(I,JJM-2,K,I_M)                &
     &                             +MOIST(I+1,JJM-2,K,I_M)              &
     &                             +MOIST(I,JJM,K,I_M)                  &
     &                             +MOIST(I+1,JJM,K,I_M))*0.25
            ENDDO

          ENDIF
        ENDDO
!
        DO I_M=2,N_SCALAR
          DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))
            SCALAR(I,JJM-1,K,I_M)=(SCALAR(I,JJM-2,K,I_M)                &
     &                            +SCALAR(I+1,JJM-2,K,I_M)              &
     &                            +SCALAR(I,JJM,K,I_M)                  &
     &                            +SCALAR(I+1,JJM,K,I_M))*0.25
          ENDDO
        ENDDO
!
      ENDIF
!
!***  ONE ROW EAST OF WESTERN BOUNDARY
!
      IF(W_BDY)THEN
        DO J=4,JM-3,2
!
          IF(W_BDY.AND.J>=MY_JS_GLB-JBPAD1                              &
     &            .AND.J<=MY_JE_GLB+JTPAD1)THEN
            JJ=J
            T(1,JJ,K)=(T(1,JJ-1,K)+T(2,JJ-1,K)                          &
     &                +T(1,JJ+1,K)+T(2,JJ+1,K))                         &
     &                *0.25
            Q(1,JJ,K)=(Q(1,JJ-1,K)+Q(2,JJ-1,K)                          &
     &                +Q(1,JJ+1,K)+Q(2,JJ+1,K))                         &
     &                *0.25
            Q2(1,JJ,K)=(Q2(1,JJ-1,K)+Q2(2,JJ-1,K)                       &
     &                 +Q2(1,JJ+1,K)+Q2(2,JJ+1,K))                      &
     &                 *0.25
            CWM(1,JJ,K)=(CWM(1,JJ-1,K)+CWM(2,JJ-1,K)                    &
     &                  +CWM(1,JJ+1,K)+CWM(2,JJ+1,K))                   &
     &                  *0.25
            PINT(1,JJ,K)=ETA1(K)*PDTOP                                  &
     &                  +ETA2(K)*PD(1,JJ)*RES(1,JJ)+PT
!
            DO I_M=1,N_MOIST
              IF(I_M==P_QV)THEN
                MOIST(1,JJ,K,I_M)=Q(1,JJ,K)/(1.-Q(1,JJ,K))     
              ELSE  
                MOIST(1,JJ,K,I_M)=(MOIST(1,JJ-1,K,I_M)                  &
     &                            +MOIST(2,JJ-1,K,I_M)                  &
     &                            +MOIST(1,JJ+1,K,I_M)                  &
     &                            +MOIST(2,JJ+1,K,I_M))*0.25
              ENDIF
            ENDDO    
!
            DO I_M=2,N_SCALAR
              SCALAR(1,JJ,K,I_M)=(SCALAR(1,JJ-1,K,I_M)                  &
     &                           +SCALAR(2,JJ-1,K,I_M)                  &
     &                           +SCALAR(1,JJ+1,K,I_M)                  &
     &                           +SCALAR(2,JJ+1,K,I_M))*0.25
            ENDDO
!
          ENDIF
!
        ENDDO
!
      ENDIF
!
!***  ONE ROW WEST OF EASTERN BOUNDARY
!
      IF(E_BDY)THEN
        DO J=4,JM-3,2
!
          IF(E_BDY.AND.J>=MY_JS_GLB-JBPAD1                              &
     &            .AND.J<=MY_JE_GLB+JTPAD1)THEN
            JJ=J
            T(IIM-1,JJ,K)=(T(IIM-1,JJ-1,K)+T(IIM,JJ-1,K)                &
     &                    +T(IIM-1,JJ+1,K)+T(IIM,JJ+1,K))               &
     &                    *0.25
            Q(IIM-1,JJ,K)=(Q(IIM-1,JJ-1,K)+Q(IIM,JJ-1,K)                &
     &                    +Q(IIM-1,JJ+1,K)+Q(IIM,JJ+1,K))               &
     &                    *0.25
            Q2(IIM-1,JJ,K)=(Q2(IIM-1,JJ-1,K)+Q2(IIM,JJ-1,K)             &
     &                     +Q2(IIM-1,JJ+1,K)+Q2(IIM,JJ+1,K))            &
     &                     *0.25
            CWM(IIM-1,JJ,K)=(CWM(IIM-1,JJ-1,K)+CWM(IIM,JJ-1,K)          &
     &                      +CWM(IIM-1,JJ+1,K)+CWM(IIM,JJ+1,K))         &
     &                      *0.25
            PINT(IIM-1,JJ,K)=ETA1(K)*PDTOP                              &
     &                      +ETA2(K)*PD(IIM-1,JJ)*RES(IIM-1,JJ)+PT
!
            DO I_M=1,N_MOIST
              IF(I_M==P_QV)THEN
                MOIST(IIM-1,JJ,K,I_M)=Q(IIM-1,JJ,K)/(1.-Q(IIM-1,JJ,K))
              ELSE
                MOIST(IIM-1,JJ,K,I_M)=(MOIST(IIM-1,JJ-1,K,I_M)                   &
     &                                +MOIST(IIM,JJ-1,K,I_M)                     &
     &                                +MOIST(IIM-1,JJ+1,K,I_M)                   &
     &                                +MOIST(IIM,JJ+1,K,I_M))*0.25
                ENDIF
              ENDDO
!
              DO I_M=2,N_SCALAR
                SCALAR(IIM-1,JJ,K,I_M)=(SCALAR(IIM-1,JJ-1,K,I_M)                    &
     &                                 +SCALAR(IIM,JJ-1,K,I_M)                      &
     &                                 +SCALAR(IIM-1,JJ+1,K,I_M)                    &
     &                                 +SCALAR(IIM,JJ+1,K,I_M))*0.25
              ENDDO
!
          ENDIF
!
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!
  200 CONTINUE
!
!-----------------------------------------------------------------------
      END SUBROUTINE BOCOH
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      SUBROUTINE BOCOV(GRIDID,NTSD,DT,LB,U_B,V_B,U_BT,V_BT              &  
     &                ,U,V                                              &
     &                ,IJDS,IJDE,SPEC_BDY_WIDTH                         &  ! min/max(id,jd)
     &                ,IHE,IHW,IVE,IVW                                  &
     &                ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &                ,IMS,IME,JMS,JME,KMS,KME                          &
     &                ,ITS,ITE,JTS,JTE,KTS,KTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    BOCOV       UPDATE WIND POINTS ON BOUNDARY
!   PRGRMMR: JANJIC          ORG: W/NP22     DATE: 94-03-08
!     
! ABSTRACT:
!     U AND V COMPONENTS OF THE WIND ARE UPDATED ON THE
!     DOMAIN BOUNDARY BY APPLYING THE PRE-COMPUTED
!     TENDENCIES AT EACH TIME STEP.  AN EXTRAPOLATION FROM
!     INSIDE THE DOMAIN IS USED FOR THE COMPONENT TANGENTIAL
!     TO THE BOUNDARY IF THE NORMAL COMPONENT IS OUTWARD.
!     
! PROGRAM HISTORY LOG:
!   87-~??  MESINGER   - ORIGINATOR
!   95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!   98-10-30  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!   01-03-13  BLACK      - CONVERTED TO WRF STRUCTURE
!   02-09-06  WOLFE      - MORE CONVERSION TO GLOBAL INDEXING 
!   04-11-23  BLACK      - THREADED
!   05-12-19  BLACK      - CONVERTED FROM IKJ TO IJK
!   06-06-02  GOPAL      - MODIFICATIONS FOR NESTING
!     
! USAGE: CALL BOCOH FROM SUBROUTINE SOLVE_NMM
!   INPUT ARGUMENT LIST:
!
!     NOTE THAT IDE AND JDE INSIDE ROUTINE SHOULD BE PASSED IN
!     AS WHAT WRF CONSIDERS THE UNSTAGGERED GRID DIMENSIONS; THAT
!     IS, 1 LESS THAN THE IDE AND JDE SET BY WRF FRAMEWORK, JM
!  
!   OUTPUT ARGUMENT LIST: 
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!  
!     UNIQUE: NONE
!  
!     LIBRARY: NONE
!  
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM 
!$$$  
!***********************************************************************
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
      INTEGER,INTENT(IN) :: IJDS,IJDE,SPEC_BDY_WIDTH
!
      INTEGER,DIMENSION(JMS:JME),INTENT(IN) :: IHE,IHW,IVE,IVW
!
      INTEGER,INTENT(IN) :: GRIDID
      INTEGER,INTENT(IN) :: LB,NTSD
!
      REAL,INTENT(IN) :: DT
!
      REAL,DIMENSION(IJDS:IJDE,KMS:KME,SPEC_BDY_WIDTH,4),INTENT(INOUT)  &
     &                                         :: U_B,V_B,U_BT,V_BT
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: U,V
!-----------------------------------------------------------------------
!
!***  LOCAL VARIABLES
!
      INTEGER :: I,II,IIM,IM,J,JJ,JJM,JM,K,N
      INTEGER :: MY_IS_GLB, MY_JS_GLB,MY_IE_GLB,MY_JE_GLB  
      INTEGER :: IBDY,BF,JB,IB
      INTEGER :: ILPAD1,IRPAD1,JBPAD1,JTPAD1
      LOGICAL :: E_BDY,W_BDY,N_BDY,S_BDY
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  TIME INTERPOLATION OF U AND V AT THE OUTER BOUNDARY
!-----------------------------------------------------------------------
!
      IM=IDE-IDS+1
      JM=JDE-JDS+1
      IIM=IM
      JJM=JM
!
      W_BDY=(ITS==IDS)
      E_BDY=(ITE==IDE)
      S_BDY=(JTS==JDS)
      N_BDY=(JTE==JDE)
!
      ILPAD1=1
      IF(ITS==IDS)ILPAD1=0
      IRPAD1=1
      IF(ITE==IDE)ILPAD1=0
      JBPAD1=1
      IF(JTS==JDS)JBPAD1=0
      JTPAD1=1
      IF(JTE==JDE)JTPAD1=0
!
      MY_IS_GLB=ITS
      MY_IE_GLB=ITE
      MY_JS_GLB=JTS
      MY_JE_GLB=JTE
!
!-----------------------------------------------------------------------
!***  SOUTH AND NORTH BOUNDARIES
!***  USE IBDY=1 FOR SOUTH; 2 FOR NORTH.
!-----------------------------------------------------------------------
!
      DO IBDY=1,2  
!
!***  MAKE SURE THE PROCESSOR HAS THIS BOUNDARY.
!
        IF((S_BDY.AND.IBDY==1).OR.(N_BDY.AND.IBDY==2))THEN
!
          IF(IBDY==1)THEN 
            BF=P_YSB     ! Which boundary (YSB=the boundary where Y is at its start)
            JB=1         ! Which cell in from Boundary 
            JJ=1         ! Which cell in the Domain
          ELSE
            BF=P_YEB     ! Which boundary (YEB=the boundary where Y is at its end)
            JB=1         ! Which cell in from Boundary
            JJ=JJM       ! Which cell in the Domain
          ENDIF
!
!$omp parallel do                                                       &
!$omp& private(i,k)
          DO K=KTS,KTE
            DO I=MAX(ITS-1,IDS),MIN(ITE+1,IDE)
              U_B(I,K,JB,BF)=U_B(I,K,JB,BF)+U_BT(I,K,JB,BF)*DT
              V_B(I,K,JB,BF)=V_B(I,K,JB,BF)+V_BT(I,K,JB,BF)*DT
              U(I,JJ,K)=U_B(I,K,JB,BF)
              V(I,JJ,K)=V_B(I,K,JB,BF)
            ENDDO
          ENDDO
!
        ENDIF
      ENDDO

!
!-----------------------------------------------------------------------
!***  WEST AND EAST BOUNDARIES
!***  USE IBDY=1 FOR WEST; 2 FOR EAST.
!-----------------------------------------------------------------------
!
      DO IBDY=1,2    
!
!***  MAKE SURE THE PROCESSOR HAS THIS BOUNDARY.
!
        IF((W_BDY.AND.IBDY==1).OR.(E_BDY.AND.IBDY==2))THEN
!
          IF(IBDY==1)THEN 
            BF=P_XSB     ! Which boundary (YSB=the boundary where Y is at its start)
            IB=1         ! Which cell in from boundary
            II=1         ! Which cell in the domain
          ELSE
            BF=P_XEB     ! Which boundary (YEB=the boundary where Y is at its end)
            IB=1         ! Which cell in from boundary
            II=IIM       ! Which cell in the domain
          ENDIF
!
!$omp parallel do                                                       &
!$omp& private(j,k)
          DO K=KTS,KTE
            DO J=MAX(JTS-1,JDS+2-1),MIN(JTE+1,JDE-1)
              IF(MOD(J,2)==0)THEN
                U_B(J,K,IB,BF)=U_B(J,K,IB,BF)+U_BT(J,K,IB,BF)*DT
                V_B(J,K,IB,BF)=V_B(J,K,IB,BF)+V_BT(J,K,IB,BF)*DT
                U(II,J,K)=U_B(J,K,IB,BF)
                V(II,J,K)=V_B(J,K,IB,BF)
              ENDIF
            ENDDO
          ENDDO
!
        ENDIF
      ENDDO

!
!-----------------------------------------------------------------------
!***  EXTRAPOLATION OF TANGENTIAL VELOCITY AT OUTFLOW POINTS
!***  BASED ON SOME DISCUSSIONS WITH ZAVISA, AND MY EXPERIMENTS
!***  ON GRAVITY PULSE FOR NESTED DOMAIN.
!-----------------------------------------------------------------------
!
      IF(GRIDID/=1)GO TO 201
!
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(i,j,jj,k)
      DO 200 K=KTS,KTE
!
!-----------------------------------------------------------------------
!
!***  SOUTHERN BOUNDARY
!
      IF(S_BDY)THEN
        DO I=max(ids+(1),its-(1)),min(ide-(2),ite+(1))
          IF(V(I,1,K)<0.)U(I,1,K)=2.*U(I,3,K)-U(I,5,K)
        ENDDO
      ENDIF
!
!***  NORTHERN BOUNDARY
!
      IF(N_BDY)THEN
        DO I=max(ids+(1),its-(1)),min(ide-(2),ite+(1))
          IF(V(I,JJM,K)>0.)                                             &
     &        U(I,JJM,K)=2.*U(I,JJM-2,K)-U(I,JJM-4,K)
        ENDDO
      ENDIF
!
!***  WESTERN BOUNDARY
!
      DO J=4,JM-3,2
        IF(W_BDY)THEN
!
          IF(W_BDY.AND.J>=MY_JS_GLB-JBPAD1                              &
     &            .AND.J<=MY_JE_GLB+JTPAD1)THEN
            JJ=J
            IF(U(1,JJ,K)<0.)                                            &
     &          V(1,JJ,K)=2.*V(2,JJ,K)-V(3,JJ,K)
          ENDIF
!
        ENDIF
      ENDDO
!
!***  EASTERN BOUNDARY
!
      DO J=4,JM-3,2
        IF(E_BDY)THEN
!
          IF(E_BDY.AND.J>=MY_JS_GLB-JBPAD1                              &
     &            .AND.J<=MY_JE_GLB+JTPAD1)THEN
            JJ=J
            IF(U(IIM,JJ,K)>0.)                                          &
     &          V(IIM,JJ,K)=2.*V(IIM-1,JJ,K)-V(IIM-2,JJ,K)
          ENDIF
!
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!
  200 CONTINUE

  201 CONTINUE
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  SPACE INTERPOLATION OF U AND V AT THE INNER BOUNDARY
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(i,j,jj,k)
      DO 300 K=KTS,KTE
!
!-----------------------------------------------------------------------
!
!***  SOUTHWEST CORNER
!
      IF(S_BDY.AND.W_BDY)THEN
        U(2,2,K)=D06666*(4.*(U(1,1,K)+U(2,1,K)+U(2,3,K))                &
     &                     + U(1,2,K)+U(1,4,K)+U(2,4,K))
        V(2,2,K)=D06666*(4.*(V(1,1,K)+V(2,1,K)+V(2,3,K))                &
     &                      +V(1,2,K)+V(1,4,K)+V(2,4,K))
      ENDIF
!
!***  SOUTHEAST CORNER
!
      IF(S_BDY.AND.E_BDY)THEN
        U(IIM-1,2,K)=D06666*(4.*(U(IIM-2,1,K)+U(IIM-1,1,K)              &
     &                          +U(IIM-2,3,K))                          &
     &                          +U(IIM,2,K)+U(IIM,4,K)+U(IIM-1,4,K))
        V(IIM-1,2,K)=D06666*(4.*(V(IIM-2,1,K)+V(IIM-1,1,K)              &
     &                          +V(IIM-2,3,K))                          &
     &                          +V(IIM,2,K)+V(IIM,4,K)+V(IIM-1,4,K))
      ENDIF
!
!***  NORTHWEST CORNER
!
      IF(N_BDY.AND.W_BDY)THEN
        U(2,JJM-1,K)=D06666*(4.*(U(1,JJM,K)+U(2,JJM,K)+U(2,JJM-2,K))    &
     &                          +U(1,JJM-1,K)+U(1,JJM-3,K)              &
     &                          +U(2,JJM-3,K))
        V(2,JJM-1,K)=D06666*(4.*(V(1,JJM,K)+V(2,JJM,K)+V(2,JJM-2,K))    &
     &                          +V(1,JJM-1,K)+V(1,JJM-3,K)              &
     &                          +V(2,JJM-3,K))
      ENDIF
!
!***  NORTHEAST CORNER
!
      IF(N_BDY.AND.E_BDY)THEN
        U(IIM-1,JJM-1,K)=                                               &
     &    D06666*(4.*(U(IIM-2,JJM,K)+U(IIM-1,JJM,K)+U(IIM-2,JJM-2,K))   &
     &               +U(IIM,JJM-1,K)+U(IIM,JJM-3,K)+U(IIM-1,JJM-3,K))
        V(IIM-1,JJM-1,K)=                                               &
     &    D06666*(4.*(V(IIM-2,JJM,K)+V(IIM-1,JJM,K)+V(IIM-2,JJM-2,K))   &
     &               +V(IIM,JJM-1,K)+V(IIM,JJM-3,K)+V(IIM-1,JJM-3,K))
      ENDIF
!
!-----------------------------------------------------------------------
!***  SPACE INTERPOLATION OF U AND V AT THE INNER BOUNDARY
!-----------------------------------------------------------------------
!
!***  ONE ROW NORTH OF SOUTHERN BOUNDARY
!
      IF(S_BDY)THEN
        DO I=max(ids+(2),its-(0)),min(ide-(2),ite+(0))
          U(I,2,K)=(U(I-1,1,K)+U(I,1,K)+U(I-1,3,K)+U(I,3,K))*0.25
          V(I,2,K)=(V(I-1,1,K)+V(I,1,K)+V(I-1,3,K)+V(I,3,K))*0.25
        ENDDO
      ENDIF
!
!***  ONE ROW SOUTH OF NORTHERN BOUNDARY
!
      IF(N_BDY)THEN
        DO I=max(ids+(2),its-(0)),min(ide-(2),ite+(0))
          U(I,JJM-1,K)=(U(I-1,JJM-2,K)+U(I,JJM-2,K)                     &
     &                 +U(I-1,JJM,K)+U(I,JJM,K))*0.25
          V(I,JJM-1,K)=(V(I-1,JJM-2,K)+V(I,JJM-2,K)                     &
     &                 +V(I-1,JJM,K)+V(I,JJM,K))*0.25
        ENDDO
      ENDIF
!
!***  ONE ROW EAST OF WESTERN BOUNDARY
!
      DO J=3,JM-2,2
        IF(W_BDY)THEN
          IF(W_BDY.AND.J>=MY_JS_GLB-JBPAD1                              &
     &            .AND.J<=MY_JE_GLB+JTPAD1)THEN
            JJ=J
            U(1,JJ,K)=(U(1,JJ-1,K)+U(2,JJ-1,K)                          &
     &                +U(1,JJ+1,K)+U(2,JJ+1,K))*0.25
            V(1,JJ,K)=(V(1,JJ-1,K)+V(2,JJ-1,K)                          &
     &                +V(1,JJ+1,K)+V(2,JJ+1,K))*0.25
          ENDIF
        ENDIF
      ENDDO
!
!***  ONE ROW WEST OF EASTERN BOUNDARY
!
      IF(E_BDY)THEN
        DO J=3,JM-2,2
          IF(E_BDY.AND.J>=MY_JS_GLB-JBPAD1                              &
     &            .AND.J<=MY_JE_GLB+JTPAD1)THEN
            JJ=J
            U(IIM-1,JJ,K)=0.25*(U(IIM-1,JJ-1,K)+U(IIM,JJ-1,K)           &
     &                         +U(IIM-1,JJ+1,K)+U(IIM,JJ+1,K))
            V(IIM-1,JJ,K)=0.25*(V(IIM-1,JJ-1,K)+V(IIM,JJ-1,K)           &
     &                         +V(IIM-1,JJ+1,K)+V(IIM,JJ+1,K))
          ENDIF
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!
  300 CONTINUE
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE BOCOV
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_BNDRY_COND
!
!-----------------------------------------------------------------------