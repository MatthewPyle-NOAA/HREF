!----------------------------------------------------------------------
!#define BIT_FOR_BIT
!----------------------------------------------------------------------
! these define the various loop range variables
! that were defined in module_MPP. Defined as macros
! here to allow thread-safety/tile callability


! these define the various loop range variables
! that were defined in module_MPP. Defined as macros
! here to allow thread-safety/tile callability





!----------------------------------------------------------------------
!
!NCEP_MESO:MODEL_LAYER: HORIZONTAL AND VERTICAL ADVECTION
!
!----------------------------------------------------------------------
!
      MODULE MODULE_ADVECTION
!
!----------------------------------------------------------------------
      USE MODULE_MODEL_CONSTANTS
      USE MODULE_EXT_INTERNAL
!----------------------------------------------------------------------
      INCLUDE "mpif.h"
!----------------------------------------------------------------------
!
      REAL,PARAMETER :: FF2=-0.64813,FF3=0.24520,FF4=-0.12189
      REAL,PARAMETER :: FFC=1.533,FBC=1.-FFC
      REAL :: CONSERVE_MIN=0.9,CONSERVE_MAX=1.1
!
!----------------------------------------------------------------------
!***  CRANK-NICHOLSON OFF-CENTER WEIGHTS FOR CURRENT AND FUTURE
!***  TIME LEVELS.
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: WGT1=0.90
      REAL,PARAMETER :: WGT2=2.-WGT1
!
!***  FOR CRANK_NICHOLSON CHECK ONLY.
!
      INTEGER :: ITEST=47,JTEST=70
      REAL :: ADTP,ADUP,ADVP,TTLO,TTUP,TULO,TUUP,TVLO,TVUP
!
!----------------------------------------------------------------------
      CONTAINS
!
!*********************************************************************** 
      SUBROUTINE ADVE(NTSD,DT,DETA1,DETA2,PDTOP                         &
     &               ,CURV,F,FAD,F4D,EM_LOC,EMT_LOC,EN,ENT,DX,DY        &
     &               ,HBM2,VBM2                                         &
     &               ,T,U,V,PDSLO,TOLD,UOLD,VOLD                        &
     &               ,PETDT,UPSTRM                                      &
     &               ,FEW,FNS,FNE,FSE                                   &
     &               ,ADT,ADU,ADV                                       &
     &               ,N_IUP_H,N_IUP_V                                   &
     &               ,N_IUP_ADH,N_IUP_ADV                               &
     &               ,IUP_H,IUP_V,IUP_ADH,IUP_ADV                       &
     &               ,IHE,IHW,IVE,IVW                                   &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &               ,IMS,IME,JMS,JME,KMS,KME                           &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    ADVE        HORIZONTAL AND VERTICAL ADVECTION
!   PRGRMMR: JANJIC          ORG: W/NP22     DATE: 93-10-28       
!     
! ABSTRACT:
!     ADVE CALCULATES THE CONTRIBUTION OF THE HORIZONTAL AND VERTICAL
!     ADVECTION TO THE TENDENCIES OF TEMPERATURE AND WIND AND THEN
!     UPDATES THOSE VARIABLES.
!     THE JANJIC ADVECTION SCHEME FOR THE ARAKAWA E GRID IS USED
!     FOR ALL VARIABLES INSIDE THE FIFTH ROW.  AN UPSTREAM SCHEME
!     IS USED ON ALL VARIABLES IN THE THIRD, FOURTH, AND FIFTH
!     OUTERMOST ROWS.  THE ADAMS-BASHFORTH TIME SCHEME IS USED.
!     
! PROGRAM HISTORY LOG:
!   87-06-??  JANJIC       - ORIGINATOR
!   95-03-25  BLACK        - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!   96-03-28  BLACK        - ADDED EXTERNAL EDGE
!   98-10-30  BLACK        - MODIFIED FOR DISTRIBUTED MEMORY
!   99-07-    JANJIC       - CONVERTED TO ADAMS-BASHFORTH SCHEME
!                            COMBINING HORIZONTAL AND VERTICAL ADVECTION
!   02-02-04  BLACK        - ADDED VERTICAL CFL CHECK
!   02-02-05  BLACK        - CONVERTED TO WRF FORMAT
!   02-08-29  MICHALAKES   - CONDITIONAL COMPILATION OF MPI
!                            CONVERT TO GLOBAL INDEXING
!   02-09-06  WOLFE        - MORE CONVERSION TO GLOBAL INDEXING
!   04-05-29  JANJIC,BLACK - CRANK-NICHOLSON VERTICAL ADVECTION
!   04-11-23  BLACK        - THREADED 
!   05-12-14  BLACK        - CONVERTED FROM IKJ TO IJK
!     
! USAGE: CALL ADVE FROM SUBROUTINE SOLVE_NMM
!   INPUT ARGUMENT LIST:
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
!   MACHINE : IBM SP
!$$$  
!***********************************************************************
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER, DIMENSION(JMS:JME),INTENT(IN) :: IHE,IHW,IVE,IVW         &
                                               ,N_IUP_H,N_IUP_V         &
     &                                         ,N_IUP_ADH,N_IUP_ADV
!
      INTEGER, DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: IUP_H,IUP_V     &
     &                                                 ,IUP_ADH,IUP_ADV
!
      INTEGER,INTENT(IN) :: NTSD
!
      REAL,INTENT(IN) :: DT,DY,EN,ENT,F4D,PDTOP
!
      REAL,DIMENSION(2600),INTENT(IN) :: EM_LOC,EMT_LOC
!
      REAL,DIMENSION(KMS:KME),INTENT(IN) :: DETA1,DETA2
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CURV,DX,F,FAD,HBM2  &
     &                                             ,PDSLO,VBM2
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: ADT,ADU,ADV
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: PETDT
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: T,TOLD   &
     &                                                        ,U,UOLD   &
     &                                                        ,V,VOLD
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(OUT) :: FEW,FNE    &
     &                                                      ,FNS,FSE
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      LOGICAL :: UPSTRM
!
      INTEGER :: I,IEND,IFP,IFQ,II,IPQ,ISP,ISQ,ISTART                   &
     &          ,IUP_ADH_J,IVH,IVL                                      &
     &          ,J,J1,JA,JAK,JEND,JGLOBAL,JJ,JKNT,JP2,JSTART            &
     &          ,K,KNTI_ADH,KSTART,KSTOP                                &
     &          ,N,N_IUPH_J,N_IUPADH_J,N_IUPADV_J
!
      INTEGER :: MY_IS_GLB,MY_IE_GLB,MY_JS_GLB,MY_JE_GLB
!
      INTEGER,DIMENSION(ITS-5:ITE+5,JTS-5:JTE+5) :: ISPA,ISQA
!
      REAL :: ADPDX,ADPDY,ARRAY3_X,CFL,CFT,CFU,CFV,CMT,CMU,CMV          &
     &       ,DTE,DTQ,F0,F1,F2,F3,FEWP,FNEP,FNSP,FPP,FSEP,HM            &
     &       ,PDOP,PDOPU,PDOPV,PP                                       &
     &       ,PVVLO,PVVLOU,PVVLOV,PVVUP,PVVUPU,PVVUPV                   &
     &       ,QP,RDP,RDPU,RDPV                                          &
     &       ,TEMPA,TEMPB,TTA,TTB,UDY                                   &
     &       ,VDX,VM,VVLO,VVLOU,VVLOV,VVUP,VVUPU,VVUPV
!
      REAL,DIMENSION(ITS-5:ITE+5,JTS-5:JTE+5) :: ARRAY0,ARRAY1          &
     &                                          ,ARRAY2,ARRAY3          &
     &                                          ,DPDE,RDPD,RDPDX,RDPDY  &
     &                                          ,TEW,TNE,TNS,TSE,TST    &
     &                                          ,UNE,UNED,UEW,UNS,USE   &
     &                                          ,USED,UST               &
     &                                          ,VEW,VNE,VNS,VSE        &
     &                                          ,VST
!
      REAL,DIMENSION(ITS-5:ITE+5,JTS-5:JTE+5,KTS:KTE) :: VAD_TEND_T     &
     &                                                  ,VAD_TEND_U     &
     &                                                  ,VAD_TEND_V
!
      REAL,DIMENSION(KTS:KTE) :: CRT,CRU,CRV,DETA1_PDTOP                &
     &                          ,RCMT,RCMU,RCMV,RSTT,RSTU,RSTV          &
     &                          ,T_K,TN,U_K,UN,V_K,VN
!
!-----------------------------------------------------------------------
!***********************************************************************
!
!                         DPDE      -----  3
!                          |                      J Increasing
!                          |
!                          |                            ^
!                         FNS       -----  2            |
!                          |                            |
!                          |                            |
!                          |                            |
!                         VNS       -----  1            |
!                          |
!                          |
!                          |
!                         ADV       -----  0  ------> Current J
!                          |
!                          |
!                          |
!                         VNS       ----- -1
!                          |
!                          |
!                          |
!                         FNS       ----- -2
!                          |
!                          |
!                          |
!                         DPDE      ----- -3
!
!***********************************************************************
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      DTQ=DT*0.25
      DTE=DT*(0.5*0.25)
!
!-----------------------------------------------------------------------
!***
!***  PRECOMPUTE DETA1 TIMES PDTOP.
!***
!-----------------------------------------------------------------------
!
      DO K=KTS,KTE
        DETA1_PDTOP(K)=DETA1(K)*PDTOP
      ENDDO
!
!-----------------------------------------------------------------------
!***
!***  INITIALIZE SOME WORKING ARRAYS TO ZERO
!***
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  COMPUTE VERTICAL ADVECTION TENDENCIES USING CRANK-NICHOLSON.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  FIRST THE TEMPERATURE
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(cft,cfu,cfv,cmt,cmu,cmv,crt,cru,crv,i,k,lmhk,lmvk        &
!$omp&        ,pdop,pdopu,pdopv,pvvlo,pvvlou,pvvlov,pvvup,pvvupu,pvvupv &
!$omp&        ,rcmt,rcmu,rcmv,rdp,rdpu,rdpv,rstt,rstu,rstv,t_k,tn       &
!$omp&        ,u_k,un,v_k,vn,vvlo,vvlou,vvlov,vvup,vvupu,vvupv)
!!$omp& private(adtp,adup,advp,ttlo,ttup,tulo,tuup,tvlo,tvup)
!-----------------------------------------------------------------------
!
      main_vertical: DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
!
!-----------------------------------------------------------------------
!
        iloop_for_t: DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
!
!-----------------------------------------------------------------------
!***  EXTRACT T FROM THE COLUMN
!-----------------------------------------------------------------------
!
          DO K=KTS,KTE
            T_K(K)=T(I,J,K)
          ENDDO
!
!-----------------------------------------------------------------------
!
          PDOP=PDSLO(I,J)
          PVVLO=PETDT(I,J,KTE-1)*DTQ
          VVLO=PVVLO/(DETA1_PDTOP(KTE)+DETA2(KTE)*PDOP)
          CMT=-VVLO*WGT2+1.
          RCMT(KTE)=1./CMT
          CRT(KTE)=VVLO*WGT2
          RSTT(KTE)=-VVLO*WGT1*(T_K(KTE-1)-T_K(KTE))+T_K(KTE)
!
!-----------------------------------------------------------------------
!
          DO K=KTE-1,KTS+1,-1
            RDP=1./(DETA1_PDTOP(K)+DETA2(K)*PDOP)
            PVVUP=PVVLO
            PVVLO=PETDT(I,J,K-1)*DTQ
            VVUP=PVVUP*RDP
            VVLO=PVVLO*RDP
            CFT=-VVUP*WGT2*RCMT(K+1)
            CMT=-CRT(K+1)*CFT+((VVUP-VVLO)*WGT2+1.)
            RCMT(K)=1./CMT
            CRT(K)=VVLO*WGT2
            RSTT(K)=-RSTT(K+1)*CFT+T_K(K)                               &
       &            -(T_K(K)-T_K(K+1))*VVUP*WGT1                        &
       &            -(T_K(K-1)-T_K(K))*VVLO*WGT1
          ENDDO
!
!-----------------------------------------------------------------------
!
          PVVUP=PVVLO
          VVUP=PVVUP/(DETA1_PDTOP(KTS)+DETA2(KTS)*PDOP)
          CFT=-VVUP*WGT2*RCMT(KTS+1)
          CMT=-CRT(KTS+1)*CFT+VVUP*WGT2+1.
          CRT(KTS)=0.
          RSTT(KTS)=-(T_K(KTS)-T_K(KTS+1))*VVUP*WGT1                    &
      &              -RSTT(KTS+1)*CFT+T_K(KTS)
          TN(KTS)=RSTT(KTS)/CMT
          VAD_TEND_T(I,J,KTS)=TN(KTS)-T_K(KTS)
!
          DO K=KTS+1,KTE
            TN(K)=(-CRT(K)*TN(K-1)+RSTT(K))*RCMT(K)
            VAD_TEND_T(I,J,K)=TN(K)-T_K(K)
          ENDDO
!
!-----------------------------------------------------------------------
!***  The following section is only for checking the implicit solution
!***  using back-substitution.  Remove this section otherwise.
!-----------------------------------------------------------------------
!       if(ntsd<=10.or.ntsd>=6000)then
!       IF(I==ITEST.AND.J==JTEST)THEN
!!
!         PVVLO=PETDT(I,J,KTE-1)*DT*0.25
!         VVLO=PVVLO/(DETA1_PDTOP(KTE)+DETA2(KTE)*PDOP)
!         TTLO=VVLO*(T(I,J,KTE-1)-T(I,J,KTE)                            &
!    &              +TN(KTE-1)-TN(KTE))
!         ADTP=TTLO+TN(KTE)-T(I,J,KTE)
!         WRITE(0,*) NTSD=,NTSD, I=,ITEST, J=,JTEST, K=,KTE     &
!    &,              ADTP=,ADTP
!         WRITE(0,*) T=,T(I,J,KTE), TN=,TN(KTE)                     &
!    &,                VAD_TEND_T=,VAD_TEND_T(I,J,KTE)
!         WRITE(0,*) 
!!
!         DO K=KTE-1,KTS+1,-1
!           RDP=1./(DETA1_PDTOP(K)+DETA2(K)*PDOP)
!           PVVUP=PVVLO
!           PVVLO=PETDT(I,J,K-1)*DT*0.25
!           VVUP=PVVUP*RDP
!           VVLO=PVVLO*RDP
!           TTUP=VVUP*(T(I,J,K)-T(I,J,K+1)+TN(K)-TN(K+1))
!           TTLO=VVLO*(T(I,J,K-1)-T(I,J,K)+TN(K-1)-TN(K))
!           ADTP=TTLO+TTUP+TN(K)-T(I,J,K)
!           WRITE(0,*) NTSD=,NTSD, I=,I, J=,J, K=,K             &
!    &,                ADTP=,ADTP
!           WRITE(0,*) T=,T(I,J,K), TN=,TN(K)                       &
!    &,                VAD_TEND_T=,VAD_TEND_T(I,J,K)
!           WRITE(0,*) 
!         ENDDO
!!
!         PVVUP=PVVLO
!         VVUP=PVVUP/(DETA1_PDTOP(KTS)+DETA2(KTS)*PDOP)
!         TTUP=VVUP*(T(I,J,KTS)-T(I,J,KTS+1)+TN(KTS)-TN(KTS+1))
!         ADTP=TTUP+TN(KTS)-T(I,J,KTS)
!         WRITE(0,*) NTSD=,NTSD, I=,I, J=,J, K=,KTS             &
!    &,              ADTP=,ADTP
!         WRITE(0,*) T=,T(I,J,KTS), TN=,TN(KTS)                     &
!    &,              VAD_TEND_T=,VAD_TEND_T(I,J,KTS)
!         WRITE(0,*) 
!       ENDIF
!       endif
!
!-----------------------------------------------------------------------
!***  End of check.
!-----------------------------------------------------------------------
!
        ENDDO iloop_for_t
!
!-----------------------------------------------------------------------
!
!***  NOW VERTICAL ADVECTION OF WIND COMPONENTS
!
!-----------------------------------------------------------------------
!
        iloop_for_uv:  DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
!
!-----------------------------------------------------------------------
!***  EXTRACT U AND V FROM THE COLUMN
!-----------------------------------------------------------------------
!
          DO K=KTS,KTE
            U_K(K)=U(I,J,K)
            V_K(K)=V(I,J,K)
          ENDDO
!
!-----------------------------------------------------------------------
!
          PDOPU=(PDSLO(I+IVW(J),J)+PDSLO(I+IVE(J),J))*0.5
          PDOPV=(PDSLO(I,J-1)+PDSLO(I,J+1))*0.5
          PVVLOU=(PETDT(I+IVW(J),J,KTE-1)+PETDT(I+IVE(J),J,KTE-1))*DTE
          PVVLOV=(PETDT(I,J-1,KTE-1)+PETDT(I,J+1,KTE-1))*DTE
          VVLOU=PVVLOU/(DETA1_PDTOP(KTE)+DETA2(KTE)*PDOPU)
          VVLOV=PVVLOV/(DETA1_PDTOP(KTE)+DETA2(KTE)*PDOPV)
          CMU=-VVLOU*WGT2+1.
          CMV=-VVLOV*WGT2+1.
          RCMU(KTE)=1./CMU
          RCMV(KTE)=1./CMV
          CRU(KTE)=VVLOU*WGT2
          CRV(KTE)=VVLOV*WGT2
          RSTU(KTE)=-VVLOU*WGT1*(U_K(KTE-1)-U_K(KTE))+U_K(KTE)
          RSTV(KTE)=-VVLOV*WGT1*(V_K(KTE-1)-V_K(KTE))+V_K(KTE)
!
!-----------------------------------------------------------------------
!
          DO K=KTE-1,KTS+1,-1
            RDPU=1./(DETA1_PDTOP(K)+DETA2(K)*PDOPU)
            RDPV=1./(DETA1_PDTOP(K)+DETA2(K)*PDOPV)
            PVVUPU=PVVLOU
            PVVUPV=PVVLOV
            PVVLOU=(PETDT(I+IVW(J),J,K-1)+PETDT(I+IVE(J),J,K-1))*DTE
            PVVLOV=(PETDT(I,J-1,K-1)+PETDT(I,J+1,K-1))*DTE
            VVUPU=PVVUPU*RDPU
            VVUPV=PVVUPV*RDPV
            VVLOU=PVVLOU*RDPU
            VVLOV=PVVLOV*RDPV
            CFU=-VVUPU*WGT2*RCMU(K+1)
            CFV=-VVUPV*WGT2*RCMV(K+1)
            CMU=-CRU(K+1)*CFU+(VVUPU-VVLOU)*WGT2+1.
            CMV=-CRV(K+1)*CFV+(VVUPV-VVLOV)*WGT2+1.
            RCMU(K)=1./CMU
            RCMV(K)=1./CMV
            CRU(K)=VVLOU*WGT2
            CRV(K)=VVLOV*WGT2
            RSTU(K)=-RSTU(K+1)*CFU+U_K(K)                               &
     &              -(U_K(K)-U_K(K+1))*VVUPU*WGT1                       &
     &              -(U_K(K-1)-U_K(K))*VVLOU*WGT1
            RSTV(K)=-RSTV(K+1)*CFV+V_K(K)                               &
     &              -(V_K(K)-V_K(K+1))*VVUPV*WGT1                       &
     &              -(V_K(K-1)-V_K(K))*VVLOV*WGT1
          ENDDO
!
!-----------------------------------------------------------------------
!
          RDPU=1./(DETA1_PDTOP(KTS)+DETA2(KTS)*PDOPU)
          RDPV=1./(DETA1_PDTOP(KTS)+DETA2(KTS)*PDOPV)
          PVVUPU=PVVLOU
          PVVUPV=PVVLOV
          VVUPU=PVVUPU*RDPU
          VVUPV=PVVUPV*RDPV
          CFU=-VVUPU*WGT2*RCMU(KTS+1)
          CFV=-VVUPV*WGT2*RCMV(KTS+1)
          CMU=-CRU(KTS+1)*CFU+VVUPU*WGT2+1.
          CMV=-CRV(KTS+1)*CFV+VVUPV*WGT2+1.
          CRU(KTS)=0.
          CRV(KTS)=0.
          RSTU(KTS)=-(U_K(KTS)-U_K(KTS+1))*VVUPU*WGT1                   &
       &               -RSTU(KTS+1)*CFU+U_K(KTS)
          RSTV(KTS)=-(V_K(KTS)-V_K(KTS+1))*VVUPV*WGT1                   &
       &               -RSTV(KTS+1)*CFV+V_K(KTS)
          UN(KTS)=RSTU(KTS)/CMU
          VN(KTS)=RSTV(KTS)/CMV
          VAD_TEND_U(I,J,KTS)=UN(KTS)-U_K(KTS)
          VAD_TEND_V(I,J,KTS)=VN(KTS)-V_K(KTS)
!
          DO K=KTS+1,KTE
            UN(K)=(-CRU(K)*UN(K-1)+RSTU(K))*RCMU(K)
            VN(K)=(-CRV(K)*VN(K-1)+RSTV(K))*RCMV(K)
            VAD_TEND_U(I,J,K)=UN(K)-U_K(K)
            VAD_TEND_V(I,J,K)=VN(K)-V_K(K)
          ENDDO
!
!-----------------------------------------------------------------------
!***  The following section is only for checking the implicit solution
!***  using back-substitution.  Remove this section otherwise.
!-----------------------------------------------------------------------
!
!       if(ntsd<=10.or.ntsd>=6000)then
!       IF(I==ITEST.AND.J==JTEST)THEN
!!
!         PDOPU=(PDSLO(I+IVW(J),J)+PDSLO(I+IVE(J),J))*0.5
!         PDOPV=(PDSLO(I,J-1)+PDSLO(I,J+1))*0.5
!         PVVLOU=(PETDT(I+IVW(J),J,KTE-1)                               &
!    &           +PETDT(I+IVE(J),J,KTE-1))*DTE
!         PVVLOV=(PETDT(I,J-1,KTE-1)                                    &
!    &           +PETDT(I,J+1,KTE-1))*DTE
!         VVLOU=PVVLOU/(DETA1_PDTOP(KTE)+DETA2(KTE)*PDOPU)
!         VVLOV=PVVLOV/(DETA1_PDTOP(KTE)+DETA2(KTE)*PDOPV)
!         TULO=VVLOU*(U(I,J,KTE-1)-U(I,J,KTE)+UN(KTE-1)-UN(KTE))
!         TVLO=VVLOV*(V(I,J,KTE-1)-V(I,J,KTE)+VN(KTE-1)-VN(KTE))
!         ADUP=TULO+UN(KTE)-U(I,J,KTE)
!         ADVP=TVLO+VN(KTE)-V(I,J,KTE)
!         WRITE(0,*) NTSD=,NTSD, I=,I, J=,J, K=,KTE             &
!    &,              ADUP=,ADUP, ADVP=,ADVP
!         WRITE(0,*) U=,U(I,J,KTE), UN=,UN(KTE)                     &
!    &,              VAD_TEND_U=,VAD_TEND_U(I,KTE)                    &
!    &,              V=,V(I,J,KTE), VN=,VN(KTE)                     &
!    &,              VAD_TEND_V=,VAD_TEND_V(I,KTE)
!         WRITE(0,*) 
!!
!         DO K=KTE-1,KTS+1,-1
!           RDPU=1./(DETA1_PDTOP(K)+DETA2(K)*PDOPU)
!           RDPV=1./(DETA1_PDTOP(K)+DETA2(K)*PDOPV)
!           PVVUPU=PVVLOU
!           PVVUPV=PVVLOV
!           PVVLOU=(PETDT(I+IVW(J),J,K-1)                               &
!    &            +PETDT(I+IVE(J),J,K-1))*DTE
!           PVVLOV=(PETDT(I,J-1,K-1)+PETDT(I,J+1,K-1))*DTE
!           VVUPU=PVVUPU*RDPU
!           VVUPV=PVVUPV*RDPV
!           VVLOU=PVVLOU*RDPU
!           VVLOV=PVVLOV*RDPV
!           TUUP=VVUPU*(U(I,J,K)-U(I,J,K+1)+UN(K)-UN(K+1))
!           TVUP=VVUPV*(V(I,J,K)-V(I,J,K+1)+VN(K)-VN(K+1))
!           TULO=VVLOU*(U(I,J,K-1)-U(I,J,K)+UN(K-1)-UN(K))
!           TVLO=VVLOV*(V(I,J,K-1)-V(I,J,K)+VN(K-1)-VN(K))
!           ADUP=TUUP+TULO+UN(K)-U(I,J,K)
!           ADVP=TVUP+TVLO+VN(K)-V(I,J,K)
!           WRITE(0,*) NTSD=,NTSD, I=,ITEST, J=,JTEST, K=,K     &
!    &,                ADUP=,ADUP, ADVP=,ADVP
!           WRITE(0,*) U=,U(I,J,K), UN=,UN(K)                       &
!    &,                VAD_TEND_U=,VAD_TEND_U(I,K)                    &
!    &,                V=,V(I,J,K), VN=,VN(K)                       &
!    &,                VAD_TEND_V=,VAD_TEND_V(I,K)
!           WRITE(0,*) 
!         ENDDO
!!
!         PVVUPU=PVVLOU
!         PVVUPV=PVVLOV
!         VVUPU=PVVUPU/(DETA1_PDTOP(KTS)+DETA2(KTS)*PDOPU)
!         VVUPV=PVVUPV/(DETA1_PDTOP(KTS)+DETA2(KTS)*PDOPV)
!         TUUP=VVUPU*(U(I,J,KTS)-U(I,J,KTS+1)+UN(KTS)-UN(KTS+1))
!         TVUP=VVUPV*(V(I,J,KTS)-V(I,J,KTS+1)+VN(KTS)-VN(KTS+1))
!         ADUP=TUUP+UN(KTS)-U(I,J,KTS)
!         ADVP=TVUP+VN(KTS)-V(I,J,KTS)
!         WRITE(0,*) NTSD=,NTSD, I=,ITEST, J=,JTEST, K=,KTS     &
!    &,              ADUP=,ADUP, ADVP=,ADVP
!         WRITE(0,*) U=,U(I,J,KTS), UN=,UN(KTS)                     &
!    &,              VAD_TEND_U=,VAD_TEND_U(I,KTS)                    &
!    &,              V=,V(I,J,KTS), VN=,VN(KTS)                     &
!    &,              VAD_TEND_V=,VAD_TEND_V(I,KTS)
!         WRITE(0,*) 
!       ENDIF
!     endif
!
!-----------------------------------------------------------------------
!***  End of check.
!-----------------------------------------------------------------------
!
        ENDDO iloop_for_uv
!
!-----------------------------------------------------------------------
!
      ENDDO main_vertical
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  COMPUTE HORIZONTAL ADVECTION TENDENCIES.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(adpdx,adpdy,adt,adu,adv,array0,array1,array2,array3      &
!$omp&        ,array3_x,dpde,f0,f1,f2,f3,fewp,fnep,fnsp,fpp,fsep,hm     &
!$omp&        ,i,ifp,ifq,ii,ipq,isp,ispa,isq,isqa,iup_adh_j,j,k         &
!$omp&        ,knti_adh,n_iupadh_j,n_iupadv_j,n_iuph_j,pp,qp            &
!$omp&        ,rdpd,rdpdx,rdpdy,tew,tne,tns,tse,tst,tta,ttb             &
!$omp&        ,uew,udy,une,uned,uns,use,used,ust                        &
!$omp&        ,vdx,vew,vm,vne,vns,vse,vst)
!-----------------------------------------------------------------------
!
      main_horizontal: DO K=KTS,KTE
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(0),jts-(4)),min(jde-(0),jte+(4))
        DO I=max(ids+(0),its-(4)),min(ide-(0),ite+(4))
          DPDE(I,J)=DETA1_PDTOP(K)+DETA2(K)*PDSLO(I,J)
          RDPD(I,J)=1./DPDE(I,J)
          TST(I,J)=T(I,J,K)*FFC+TOLD(I,J,K)*FBC
          UST(I,J)=U(I,J,K)*FFC+UOLD(I,J,K)*FBC
          VST(I,J)=V(I,J,K)*FFC+VOLD(I,J,K)*FBC
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!***  MASS FLUXES AND MASS POINT ADVECTION COMPONENTS
!***  THE NS AND EW FLUXES IN THE FOLLOWING LOOP ARE ON V POINTS
!***  FOR T.
!-----------------------------------------------------------------------
!
        DO J=max(jds+(1),jts-(3)),min(jde-(1),jte+(3))
        DO I=max(ids+(0),its-(3)),min(ide-(0),ite+(3))
!
          ADPDX=DPDE(I+IVW(J),J)+DPDE(I+IVE(J),J)
          ADPDY=DPDE(I,J-1)+DPDE(I,J+1)
          RDPDX(I,J)=1./ADPDX
          RDPDY(I,J)=1./ADPDY
!
          UDY=U(I,J,K)*DY
          VDX=V(I,J,K)*DX(I,J)
!
          FEWP=UDY*ADPDX
          FNSP=VDX*ADPDY
!
          FEW(I,J,K)=FEWP
          FNS(I,J,K)=FNSP
!
          TEW(I,J)=FEWP*(TST(I+IVE(J),J)-TST(I+IVW(J),J))
          TNS(I,J)=FNSP*(TST(I,J+1)-TST(I,J-1))
!
          UNED(I,J)=UDY+VDX
          USED(I,J)=UDY-VDX
!
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!***  DIAGONAL FLUXES AND DIAGONALLY AVERAGED WIND
!***  THE NE AND SE FLUXES ARE ASSOCIATED WITH H POINTS
!***  (ACTUALLY JUST TO THE NE AND SE OF EACH H POINT).
!-----------------------------------------------------------------------
!
        DO J=max(jds+(1),jts-(2)),min(jde-(2),jte+(2))
        DO I=max(ids+(0),its-(2)),min(ide-(0),ite+(2))
          FNEP=(UNED(I+IHE(J),J)+UNED(I       ,J+1))                    &
     &        *(DPDE(I       ,J)+DPDE(I+IHE(J),J+1))
          FNE(I,J,K)=FNEP
          TNE(I,J)=FNEP*(TST(I+IHE(J),J+1)-TST(I,J))
        ENDDO
        ENDDO
!
        DO J=max(jds+(2),jts-(2)),min(jde-(1),jte+(2))
        DO I=max(ids+(0),its-(2)),min(ide-(0),ite+(2))
          FSEP=(USED(I+IHE(J),J)+USED(I       ,J-1))                    &
     &        *(DPDE(I       ,J)+DPDE(I+IHE(J),J-1))
          FSE(I,J,K)=FSEP
          TSE(I,J)=FSEP*(TST(I+IHE(J),J-1)-TST(I,J))
!
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!***  HORIZONTAL T ADVECTION TENDENCY ADT IS ON H POINTS OF COURSE.
!-----------------------------------------------------------------------
!
        DO J=max(jds+(5),jts-(0)),min(jde-(5),jte+(0))
        DO I=max(ids+(2),its-(0)),min(ide-(2),ite+(0))
          ADT(I,J)=(TEW(I+IHW(J),J)+TEW(I+IHE(J),J)                     &
     &             +TNS(I,J-1)+TNS(I,J+1)                               &
     &             +TNE(I+IHW(J),J-1)+TNE(I,J)                          &
     &             +TSE(I,J)+TSE(I+IHW(J),J+1))                         &
     &             *RDPD(I,J)*FAD(I,J)
        ENDDO
        ENDDO
!
!
!-----------------------------------------------------------------------
!***  CALCULATION OF MOMENTUM ADVECTION COMPONENTS.
!-----------------------------------------------------------------------
!
        DO J=max(jds+(4),jts-(1)),min(jde-(4),jte+(1))
        DO I=max(ids+(0),its-(1)),min(ide-(0),ite+(1))
!
!-----------------------------------------------------------------------
!***  THE NS AND EW FLUXES ARE ON H POINTS FOR U AND V.
!-----------------------------------------------------------------------
!
          UEW(I,J)=(FEW(I+IHW(J),J,K)+FEW(I+IHE(J),J,K))                &
     &            *(UST(I+IHE(J),J)-UST(I+IHW(J),J))
          UNS(I,J)=(FNS(I+IHW(J),J,K)+FNS(I+IHE(J),J,K))                &
     &            *(UST(I,J+1)-UST(I,J-1))
          VEW(I,J)=(FEW(I,J-1,K)+FEW(I,J+1,K))                          &
     &            *(VST(I+IHE(J),J)-VST(I+IHW(J),J))
          VNS(I,J)=(FNS(I,J-1,K)+FNS(I,J+1,K))                          &
     &            *(VST(I,J+1)-VST(I,J-1))
!
!-----------------------------------------------------------------------
!***  THE FOLLOWING NE AND SE FLUXES ARE TIED TO V POINTS AND ARE
!***  LOCATED JUST TO THE NE AND SE OF THE GIVEN I,J.
!-----------------------------------------------------------------------
!
          UNE(I,J)=(FNE(I+IVW(J),J,K)+FNE(I+IVE(J),J,K))                &
     &            *(UST(I+IVE(J),J+1)-UST(I,J))
          USE(I,J)=(FSE(I+IVW(J),J,K)+FSE(I+IVE(J),J,K))                &
     &            *(UST(I+IVE(J),J-1)-UST(I,J))
          VNE(I,J)=(FNE(I,J-1,K)+FNE(I,J+1,K))                          &
     &            *(VST(I+IVE(J),J+1)-VST(I,J))
          VSE(I,J)=(FSE(I,J-1,K)+FSE(I,J+1,K))                          &
     &            *(VST(I+IVE(J),J-1)-VST(I,J))
!
!-----------------------------------------------------------------------
!
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!***  COMPUTE THE ADVECTION TENDENCIES FOR U AND V.
!***  THE AD ARRAYS ARE ON THE VELOCITY POINTS.
!-----------------------------------------------------------------------
!
        DO J=max(jds+(5),jts-(0)),min(jde-(5),jte+(0))
        DO I=max(ids+(2),its-(0)),min(ide-(2),ite+(0))
          ADU(I,J)=(UEW(I+IVW(J),J)+UEW(I+IVE(J),J)                     &
     &             +UNS(I,J-1)+UNS(I,J+1)                               &
     &             +UNE(I+IVW(J),J-1)+UNE(I,J)                          &
     &             +USE(I,J)+USE(I+IVW(J),J+1))                         &
     &             *RDPDX(I,J)*FAD(I+IVW(J),J)
!
          ADV(I,J)=(VEW(I+IVW(J),J)+VEW(I+IVE(J),J)                     &
     &             +VNS(I,J-1)+VNS(I,J+1)                               &
     &             +VNE(I+IVW(J),J-1)+VNE(I,J)                          &
     &             +VSE(I,J)+VSE(I+IVW(J),J+1))                         &
     &             *RDPDY(I,J)*FAD(I+IVW(J),J)
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
!***  END OF JANJIC HORIZONTAL ADVECTION
!
!-----------------------------------------------------------------------
!
!***  UPSTREAM ADVECTION OF T
!
!-----------------------------------------------------------------------
!
        upstream: IF(UPSTRM)THEN
!
!-----------------------------------------------------------------------
!***
!***  COMPUTE UPSTREAM COMPUTATIONS ON THIS TASKS ROWS.
!***
!-----------------------------------------------------------------------
!
          jloop_upstream: DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
!
            N_IUPH_J=N_IUP_H(J)   ! See explanation in START_DOMAIN_NMM
            DO II=0,N_IUPH_J-1
!
              I=IUP_H(IMS+II,J)
              TTA=EMT_LOC(J)*(UST(I,J-1)+UST(I+IHW(J),J)                &
     &                       +UST(I+IHE(J),J)+UST(I,J+1))
              TTB=ENT       *(VST(I,J-1)+VST(I+IHW(J),J)                &
     &                       +VST(I+IHE(J),J)+VST(I,J+1))
              PP=-TTA-TTB
              QP= TTA-TTB
!
              IF(PP<0.)THEN
                ISPA(I,J)=-1
              ELSE
                ISPA(I,J)= 1
              ENDIF
!
              IF(QP<0.)THEN
                ISQA(I,J)=-1
              ELSE
                ISQA(I,J)= 1
              ENDIF
!
              PP=ABS(PP)
              QP=ABS(QP)
              ARRAY3_X=PP*QP
              ARRAY0(I,J)=ARRAY3_X-PP-QP
              ARRAY1(I,J)=PP-ARRAY3_X
              ARRAY2(I,J)=QP-ARRAY3_X
              ARRAY3(I,J)=ARRAY3_X
            ENDDO
!
!-----------------------------------------------------------------------
!
            N_IUPADH_J=N_IUP_ADH(J)
            KNTI_ADH=1
            IUP_ADH_J=IUP_ADH(IMS,J)
!
            iloop_T: DO II=0,N_IUPH_J-1
!
              I=IUP_H(IMS+II,J)
!
              ISP=ISPA(I,J)
              ISQ=ISQA(I,J)
              IFP=(ISP-1)/2
              IFQ=(-ISQ-1)/2
              IPQ=(ISP-ISQ)/2
!
!-----------------------------------------------------------------------
!
              IF(I==IUP_ADH_J)THEN  ! Upstream advection T tendencies
!
                ISP=ISPA(I,J)
                ISQ=ISQA(I,J)
                IFP=(ISP-1)/2
                IFQ=(-ISQ-1)/2
                IPQ=(ISP-ISQ)/2
!
                F0=ARRAY0(I,J)
                F1=ARRAY1(I,J)
                F2=ARRAY2(I,J)
                F3=ARRAY3(I,J)
!
                ADT(I,J)=F0*T(I,J,K)                                    &
     &                  +F1*T(I+IHE(J)+IFP,J+ISP,K)                     &
     &                  +F2*T(I+IHE(J)+IFQ,J+ISQ,K)                     &
                        +F3*T(I+IPQ,J+ISP+ISQ,K)
!
!-----------------------------------------------------------------------
!
                IF(KNTI_ADH<N_IUPADH_J)THEN
                  IUP_ADH_J=IUP_ADH(IMS+KNTI_ADH,J)
                  KNTI_ADH=KNTI_ADH+1
                ENDIF
!
              ENDIF  ! End of upstream advection T tendency IF block
!
            ENDDO iloop_T
!
!-----------------------------------------------------------------------
!
!***  UPSTREAM ADVECTION OF VELOCITY COMPONENTS
!
!-----------------------------------------------------------------------
!
            N_IUPADV_J=N_IUP_ADV(J)
!
            DO II=0,N_IUPADV_J-1
              I=IUP_ADV(IMS+II,J)
!
              TTA=EM_LOC(J)*UST(I,J)
              TTB=EN       *VST(I,J)
              PP=-TTA-TTB
              QP=TTA-TTB
!
              IF(PP<0.)THEN
                ISP=-1
              ELSE
                ISP= 1
              ENDIF
!
              IF(QP<0.)THEN
                ISQ=-1
              ELSE
                ISQ= 1
              ENDIF
!
              IFP=(ISP-1)/2
              IFQ=(-ISQ-1)/2
              IPQ=(ISP-ISQ)/2
              PP=ABS(PP)
              QP=ABS(QP)
              F3=PP*QP
              F0=F3-PP-QP
              F1=PP-F3
              F2=QP-F3
!
              ADU(I,J)=F0*U(I,J,K)                                      &
     &                +F1*U(I+IVE(J)+IFP,J+ISP,K)                       &
     &                +F2*U(I+IVE(J)+IFQ,J+ISQ,K)                       &
     &                +F3*U(I+IPQ,J+ISP+ISQ,K)
!
              ADV(I,J)=F0*V(I,J,K)                                      &
     &                +F1*V(I+IVE(J)+IFP,J+ISP,K)                       &
     &                +F2*V(I+IVE(J)+IFQ,J+ISQ,K)                       &
     &                +F3*V(I+IPQ,J+ISP+ISQ,K)
!
            ENDDO
!
          ENDDO jloop_upstream
!
!-----------------------------------------------------------------------
!
        ENDIF upstream
!
!-----------------------------------------------------------------------
!
!***  END OF HORIZONTAL ADVECTION
!
!-----------------------------------------------------------------------
!
!***  NOW SUM THE VERTICAL AND HORIZONTAL TENDENCIES,
!***  CURVATURE AND CORIOLIS TERMS.
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
          HM=HBM2(I,J)
          VM=VBM2(I,J)
          ADT(I,J)=(VAD_TEND_T(I,J,K)+2.*ADT(I,J))*HM
!
          FPP=CURV(I,J)*2.*UST(I,J)+F(I,J)*2.
          ADU(I,J)=(VAD_TEND_U(I,J,K)+2.*ADU(I,J)+VST(I,J)*FPP)*VM
          ADV(I,J)=(VAD_TEND_V(I,J,K)+2.*ADV(I,J)-UST(I,J)*FPP)*VM
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!***  SAVE THE OLD VALUES FOR TIMESTEPPING
!-----------------------------------------------------------------------
!
        DO J=max(jds+(0),jts-(4)),min(jde-(0),jte+(4))
        DO I=max(ids+(0),its-(4)),min(ide-(0),ite+(4))
          TOLD(I,J,K)=T(I,J,K)
          UOLD(I,J,K)=U(I,J,K)
          VOLD(I,J,K)=V(I,J,K)
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!***  FINALLY UPDATE THE PROGNOSTIC VARIABLES
!-----------------------------------------------------------------------
!
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
          T(I,J,K)=ADT(I,J)+T(I,J,K)
          U(I,J,K)=ADU(I,J)+U(I,J,K)
          V(I,J,K)=ADV(I,J)+V(I,J,K)
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO main_horizontal
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ADVE
!
!-----------------------------------------------------------------------
!
!***********************************************************************
      SUBROUTINE VAD2(NTSD,DT,IDTAD,DX,DY                               &
     &               ,AETA1,AETA2,DETA1,DETA2,PDSL,PDTOP,HBM2           &
     &               ,Q,Q2,CWM,PETDT                                    &
     &               ,N_IUP_H,N_IUP_V                                   &
     &               ,N_IUP_ADH,N_IUP_ADV                               &
     &               ,IUP_H,IUP_V,IUP_ADH,IUP_ADV                       &
     &               ,IHE,IHW,IVE,IVW                                   &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &               ,IMS,IME,JMS,JME,KMS,KME                           &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    VAD2        VERTICAL ADVECTION OF H2O SUBSTANCE AND TKE
!   PRGRMMR: JANJIC          ORG: W/NP22     DATE: 96-07-19
!
! ABSTRACT:
!     VAD2 CALCULATES THE CONTRIBUTION OF THE VERTICAL ADVECTION
!     TO THE TENDENCIES OF WATER SUBSTANCE AND TKE AND THEN UPDATES
!     THOSE VARIABLES.  AN ANTI-FILTERING TECHNIQUE IS USED.
!
! PROGRAM HISTORY LOG:
!   96-07-19  JANJIC   - ORIGINATOR
!   98-11-02  BLACK    - MODIFIED FOR DISTRIBUTED MEMORY
!   99-03-17  TUCCILLO - INCORPORATED MPI_ALLREDUCE FOR GLOBAL SUM
!   02-02-06  BLACK    - CONVERTED TO WRF FORMAT
!   02-09-06  WOLFE    - MORE CONVERSION TO GLOBAL INDEXING
!   04-11-23  BLACK    - THREADED
!
! USAGE: CALL VAD2 FROM SUBROUTINE SOLVE_NMM
!   INPUT ARGUMENT LIST:
!
!   OUTPUT ARGUMENT LIST
!
!   OUTPUT FILES:
!       NONE
!   SUBPROGRAMS CALLED:
!
!     UNIQUE: NONE
!
!     LIBRARY: NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM SP
!$$$
!***********************************************************************
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
                           ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,DIMENSION(JMS:JME),INTENT(IN) :: IHE,IHW,IVE,IVW
      INTEGER,DIMENSION(JMS:JME),INTENT(IN) :: N_IUP_H,N_IUP_V          &
     &                                        ,N_IUP_ADH,N_IUP_ADV
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: IUP_H,IUP_V      &
     &                                                ,IUP_ADH,IUP_ADV
!
      INTEGER,INTENT(IN) :: IDTAD,NTSD
!
      REAL,INTENT(IN) :: DT,DY,PDTOP
!
      REAL,DIMENSION(KMS:KME),INTENT(IN) :: AETA1,AETA2,DETA1,DETA2
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: DX,HBM2,PDSL
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: PETDT
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: CWM,Q,Q2
!
!----------------------------------------------------------------------
!***  LOCAL VARIABLES
!----------------------------------------------------------------------
!
      REAL,PARAMETER :: FF1=0.500
!
      LOGICAL,SAVE :: TRADITIONAL=.TRUE.
!
      INTEGER :: I,IRECV,J,JFP,JFQ,K,LAP,LLAP
!
      INTEGER,DIMENSION(KTS:KTE) :: LA
!
      REAL*8 :: ADDT,AFRP,D2PQE,D2PQQ,D2PQW,DEP,DETAP,DPDN,DPUP,DQP     &
     &       ,DWP,E00,E4P,EP,EP0,HADDT,HBM2IJ                           &
     &       ,Q00,Q4P,QP,QP0                                            &
     &       ,RFACEK,RFACQK,RFACWK,RFC,RR                               &
     &       ,SUMNE,SUMNQ,SUMNW,SUMPE,SUMPQ,SUMPW                       &
     &       ,W00,W4P,WP,WP0

      REAL :: SFACEK,SFACQK,SFACWK
!
      REAL,DIMENSION(KTS:KTE) :: AFR,DEL,DQL,DWL,E3,E4,PETDTK           &
     &                          ,RFACE,RFACQ,RFACW,Q3,Q4,W3,W4
!
!***********************************************************************
!-----------------------------------------------------------------------
!
      ADDT=REAL(IDTAD)*DT
!
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(afr,afrp,d2pqe,d2pqq,d2pqw,del,dep,detap,dpdn,dpup       &
!$omp&        ,dql,dqp,dwl,dwp,e00,e3,e4,e4p,ep,ep0,haddt,i,j,k         &
!$omp&        ,la,lap,llap,petdtk,q00,q3,q4,q4p,qp,qp0,rfacek,rfacqk    &
!$omp&        ,rfacwk,rfc,rr,sumne,sumnq,sumnw,sumpe,sumpq,sumpw        &
!$omp&        ,w00,w3,w4,w4p,wp,wp0,sfacek,sfacqk,sfacwk)
!-----------------------------------------------------------------------
!
      main_integration : DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
!
!-----------------------------------------------------------------------
!
        main_iloop: DO I=max(ids+(1),its-(1)),min(ide-(1),ite+(1))
!
!-----------------------------------------------------------------------
!
          E3(KTE)=Q2(I,J,KTE)*0.5
!
          DO K=KTE-1,KTS,-1
            E3(K)=MAX((Q2(I,J,K+1)+Q2(I,J,K))*0.5,EPSQ2)
          ENDDO
!
          DO K=KTS,KTE
            Q3(K)=MAX(Q(I,J,K),EPSQ)
            W3(K)=MAX(CWM(I,J,K),CLIMIT)
            E4(K)=E3(K)
            Q4(K)=Q3(K)
            W4(K)=W3(K)
          ENDDO
!
          IF(TRADITIONAL)THEN
            PETDTK(KTE)=PETDT(I,J,KTE-1)*0.5
!
            DO K=KTE-1,KTS+1,-1
              PETDTK(K)=(PETDT(I,J,K)+PETDT(I,J,K-1))*0.5
            ENDDO
!
            PETDTK(KTS)=PETDT(I,J,KTS)*0.5
!
          ELSE   
!
!-----------------------------------------------------------------------
!***	PERFORM HORIZONTAL AVERAGING OF VERTICAL VELOCITY
!-----------------------------------------------------------------------
!
            PETDTK(KTE)=(PETDT(I+IHW(J-1),J-1,KTE-1)                    &
     &                  +PETDT(I+IHE(J-1),J-1,KTE-1)                    &
     &                  +PETDT(I+IHW(J+1),J+1,KTE-1)                    &
     &                  +PETDT(I+IHE(J+1),J+1,KTE-1)                    &
     &                  +PETDT(I,J,KTE-1)*4.        )*0.0625
!
            DO K=KTE-1,KTS+1,-1
              PETDTK(K)=(PETDT(I+IHW(J-1),J-1,K-1)                      &
                        +PETDT(I+IHE(J-1),J-1,K-1)                      &
     &                  +PETDT(I+IHW(J+1),J+1,K-1)                      &
     &                  +PETDT(I+IHE(J+1),J+1,K-1)                      &
     &                  +PETDT(I+IHW(J-1),J-1,K  )                      &
     &                  +PETDT(I+IHE(J-1),J-1,K  )                      &
     &                  +PETDT(I+IHW(J+1),J+1,K  )                      &
     &                  +PETDT(I+IHE(J+1),J+1,K  )                      &
     &                  +(PETDT(I,J,K-1)+PETDT(I,J,K))*4.               &
     &                                                   )*0.0625
            ENDDO
!
            PETDTK(KTS)=(PETDT(I+IHW(J-1),J-1,KTS)                      &
     &                  +PETDT(I+IHE(J-1),J-1,KTS)                      &
     &                  +PETDT(I+IHW(J+1),J+1,KTS)                      &
     &                  +PETDT(I+IHE(J+1),J+1,KTS)                      &
     &                  +PETDT(I,J,KTS)*4.        )*0.0625

          ENDIF
!
!-----------------------------------------------------------------------
!
          HADDT=-ADDT*HBM2(I,J)
!
          DO K=KTE,KTS,-1
            RR=PETDTK(K)*HADDT
!
            IF(RR<0.)THEN
              LAP=1
            ELSE
              LAP=-1
            ENDIF
!
            LA(K)=LAP
            LLAP=K+LAP
!
            IF(LLAP>KTS-1.AND.LLAP<KTE+1)THEN                             !zjmod
              RR=ABS(RR/((AETA1(LLAP)-AETA1(K))*PDTOP                   &
     &                  +(AETA2(LLAP)-AETA2(K))*PDSL(I,J)))
              IF(RR>0.9)RR=0.9
!
              AFR(K)=(((FF4*RR+FF3)*RR+FF2)*RR+FF1)*RR
              DQP=(Q3(LLAP)-Q3(K))*RR
              DWP=(W3(LLAP)-W3(K))*RR
              DEP=(E3(LLAP)-E3(K))*RR
              DQL(K)=DQP
              DWL(K)=DWP
              DEL(K)=DEP
            ELSE
              RR=0.
              AFR(K)=0.
              DQL(K)=0.
              DWL(K)=0.
              DEL(K)=0.
            ENDIF
          ENDDO
!
!-----------------------------------------------------------------------
!
          IF(LA(KTE-1)>0)THEN
            RFC=(DETA1(KTE-1)*PDTOP+DETA2(KTE-1)*PDSL(I,J))             &
     &         /(DETA1(KTE  )*PDTOP+DETA2(KTE  )*PDSL(I,J))
            DQL(KTE)=-DQL(KTE-1)*RFC
            DWL(KTE)=-DWL(KTE-1)*RFC
            DEL(KTE)=-DEL(KTE-1)*RFC
          ENDIF
!
          IF(LA(KTS+1)<0)THEN
            RFC=(DETA1(KTS+1)*PDTOP+DETA2(KTS+1)*PDSL(I,J))             &
     &         /(DETA1(KTS  )*PDTOP+DETA2(KTS  )*PDSL(I,J))
            DQL(KTS)=-DQL(KTS+1)*RFC
            DWL(KTS)=-DWL(KTS+1)*RFC
            DEL(KTS)=-DEL(KTS+1)*RFC
          ENDIF
!
          DO K=KTS,KTE
            Q4(K)=Q3(K)+DQL(K)
            W4(K)=W3(K)+DWL(K)
            E4(K)=E3(K)+DEL(K)
          ENDDO
!
!-----------------------------------------------------------------------
!***  ANTI-FILTERING STEP
!-----------------------------------------------------------------------
!
          SUMPQ=0.
          SUMNQ=0.
          SUMPW=0.
          SUMNW=0.
          SUMPE=0.
          SUMNE=0.
!
          antifiltering_limiters: DO K=KTE-1,KTS+1,-1
!
            DETAP=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
!
            Q4P=Q4(K)
            W4P=W4(K)
            E4P=E4(K)
!
            LAP=LA(K)
!
            DPDN=(AETA1(K+LAP)-AETA1(K))*PDTOP                          &
     &          +(AETA2(K+LAP)-AETA2(K))*PDSL(I,J)
            DPUP=(AETA1(K)-AETA1(K-LAP))*PDTOP                          &
     &          +(AETA2(K)-AETA2(K-LAP))*PDSL(I,J)
!
            AFRP=2.*AFR(K)*DPDN*DPDN/(DPDN+DPUP)
            D2PQQ=((Q4(K+LAP)-Q4P)/DPDN                                 &
     &            -(Q4P-Q4(K-LAP))/DPUP)*AFRP
            D2PQW=((W4(K+LAP)-W4P)/DPDN                                 &
     &            -(W4P-W4(K-LAP))/DPUP)*AFRP
            D2PQE=((E4(K+LAP)-E4P)/DPDN                                 &
     &            -(E4P-E4(K-LAP))/DPUP)*AFRP
!
            QP=Q4P-D2PQQ
            WP=W4P-D2PQW
            EP=E4P-D2PQE
!
            Q00=Q3(K)
            QP0=Q3(K+LAP)
!
            W00=W3(K)
            WP0=W3(K+LAP)
!
            E00=E3(K)
            EP0=E3(K+LAP)
!
            QP=MAX(QP,MIN(Q00,QP0))
            QP=MIN(QP,MAX(Q00,QP0))
            WP=MAX(WP,MIN(W00,WP0))
            WP=MIN(WP,MAX(W00,WP0))
            EP=MAX(EP,MIN(E00,EP0))
            EP=MIN(EP,MAX(E00,EP0))
!
            DQP=QP-Q00
            DWP=WP-W00
            DEP=EP-E00
!
            DQL(K)=DQP
            DWL(K)=DWP
            DEL(K)=DEP
!
          ENDDO antifiltering_limiters
!
!-----------------------------------------------------------------------
!
          IF(LA(KTE-1)>0)THEN
            RFC=(DETA1(KTE-1)*PDTOP+DETA2(KTE-1)*PDSL(I,J))             &
     &         /(DETA1(KTE  )*PDTOP+DETA2(KTE  )*PDSL(I,J))
            DQL(KTE)=-DQL(KTE-1)*RFC+DQL(KTE)
            DWL(KTE)=-DWL(KTE-1)*RFC+DWL(KTE)
            DEL(KTE)=-DEL(KTE-1)*RFC+DEL(KTE)
          ENDIF
!
          IF(LA(KTS+1)<0)THEN
            RFC=(DETA1(KTS+1)*PDTOP+DETA2(KTS+1)*PDSL(I,J))             &
     &         /(DETA1(KTS  )*PDTOP+DETA2(KTS  )*PDSL(I,J))
            DQL(KTS)=-DQL(KTS+1)*RFC+DQL(KTS)
            DWL(KTS)=-DWL(KTS+1)*RFC+DWL(KTS)
            DEL(KTS)=-DEL(KTS+1)*RFC+DEL(KTS)
          ENDIF
!
          DO K=KTS,KTE
            DETAP=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
            DQP=DQL(K)*DETAP 
            DWP=DWL(K)*DETAP
            DEP=DEL(K)*DETAP
!
            IF(DQP>0.)THEN
              SUMPQ=SUMPQ+DQP
            ELSE
              SUMNQ=SUMNQ+DQP
            ENDIF
            IF(DWP>0.)THEN
              SUMPW=SUMPW+DWP
            ELSE
              SUMNW=SUMNW+DWP
            ENDIF
            IF(DEP>0.)THEN
              SUMPE=SUMPE+DEP
            ELSE
              SUMNE=SUMNE+DEP
            ENDIF
          ENDDO
!
!-----------------------------------------------------------------------
!***  FIRST MOMENT CONSERVING FACTOR
!-----------------------------------------------------------------------
!
          IF(SUMPQ>1.E-9)THEN
            SFACQK=-SUMNQ/SUMPQ
          ELSE
            SFACQK=1.
          ENDIF
!
          IF(SUMPW>1.E-9)THEN
            SFACWK=-SUMNW/SUMPW
          ELSE
            SFACWK=1.
          ENDIF
!
          IF(SUMPE>1.E-9)THEN
            SFACEK=-SUMNE/SUMPE
          ELSE
            SFACEK=1.
          ENDIF
!
          IF(SFACQK<CONSERVE_MIN.OR.SFACQK>CONSERVE_MAX)SFACQK=1.
          IF(SFACWK<CONSERVE_MIN.OR.SFACWK>CONSERVE_MAX)SFACWK=1.
          IF(SFACEK<CONSERVE_MIN.OR.SFACEK>CONSERVE_MAX)SFACEK=1.
!
          RFACQK=1./SFACQK
          RFACWK=1./SFACWK
          RFACEK=1./SFACEK
!
!-----------------------------------------------------------------------
!***  IMPOSE CONSERVATION ON ANTI-FILTERING
!-----------------------------------------------------------------------
!
          DO K=KTE,KTS,-1
!
            DQP=DQL(K)
            IF(SFACQK>=1.)THEN
              IF(DQP<0.)DQP=DQP*RFACQK
            ELSE
              IF(DQP>0.)DQP=DQP*SFACQK
            ENDIF
            Q(I,J,K)=Q3(K)+DQP
!
            DWP=DWL(K)
            IF(SFACWK>=1.)THEN
              IF(DWP<0.)DWP=DWP*RFACWK
            ELSE
              IF(DWP>0.)DWP=DWP*SFACWK
            ENDIF  
            CWM(I,J,K)=W3(K)+DWP
!
            DEP=DEL(K)
            IF(SFACEK>=1.)THEN
              IF(DEP<0.)DEP=DEP*RFACEK
            ELSE
              IF(DEP>0.)DWP=DWP*SFACEK
            ENDIF  
            E3(K)=E3(K)+DEP
!
          ENDDO
!
!-----------------------------------------------------------------------
!
          HBM2IJ=HBM2(I,J)
          Q2(I,J,KTE)=MAX(E3(KTE)+E3(KTE)-EPSQ2,EPSQ2)*HBM2IJ           &
     &               +Q2(I,J,KTE)*(1.-HBM2IJ)
          DO K=KTE-1,KTS+1,-1
            Q2(I,J,K)=MAX(E3(K)+E3(K)-Q2(I,J,K+1),EPSQ2)*HBM2IJ         &
     &               +Q2(I,J,K)*(1.-HBM2IJ)
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDDO main_iloop 
!
!-----------------------------------------------------------------------
!
      ENDDO main_integration
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE VAD2
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      SUBROUTINE HAD2(                                                  &
     &                domdesc ,                                         &
     &                NTSD,DT,IDTAD,DX,DY                               &
     &               ,AETA1,AETA2,DETA1,DETA2,PDSL,PDTOP                &
     &               ,HBM2,HBM3                                         &
     &               ,Q,Q2,CWM,U,V,Z,HYDRO                              &
     &               ,N_IUP_H,N_IUP_V                                   &
     &               ,N_IUP_ADH,N_IUP_ADV                               &
     &               ,IUP_H,IUP_V,IUP_ADH,IUP_ADV                       &
     &               ,IHE,IHW,IVE,IVW                                   &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &               ,IMS,IME,JMS,JME,KMS,KME                           &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    HAD2        HORIZONTAL ADVECTION OF H2O AND TKE
!   PRGRMMR: JANJIC          ORG: W/NP22     DATE: 96-07-19
!
! ABSTRACT:
!     HAD2 CALCULATES THE CONTRIBUTION OF THE HORIZONTAL ADVECTION
!     TO THE TENDENCIES OF WATER SUBSTANCE AND TKE AND THEN
!     UPDATES THOSE VARIABLES.  AN ANTI-FILTERING TECHNIQUE IS USED.
!
! PROGRAM HISTORY LOG:
!   96-07-19  JANJIC   - ORIGINATOR
!   98-11-02  BLACK    - MODIFIED FOR DISTRIBUTED MEMORY
!   99-03-17  TUCCILLO - INCORPORATED MPI_ALLREDUCE FOR GLOBAL SUM
!   02-02-06  BLACK    - CONVERTED TO WRF FORMAT
!   02-09-06  WOLFE    - MORE CONVERSION TO GLOBAL INDEXING
!   03-05-23  JANJIC   - ADDED SLOPE FACTOR
!   04-11-23  BLACK    - THREADED
!   05-12-14  BLACK    - CONVERTED FROM IKJ TO IJK
!
! USAGE: CALL HAD2 FROM SUBROUTINE SOLVE_NMM
!   INPUT ARGUMENT LIST:
!
!   OUTPUT ARGUMENT LIST
!
!   OUTPUT FILES:
!       NONE
!   SUBPROGRAMS CALLED:
!
!     UNIQUE: NONE
!
!     LIBRARY: NONE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM SP
!$$$
!***********************************************************************
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,DIMENSION(JMS:JME),INTENT(IN) :: IHE,IHW,IVE,IVW
      INTEGER,DIMENSION(JMS:JME),INTENT(IN) :: N_IUP_H,N_IUP_V          &
     &                                        ,N_IUP_ADH,N_IUP_ADV
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: IUP_H,IUP_V      &
     &                                                ,IUP_ADH,IUP_ADV
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IDTAD,NTSD
!
      REAL,INTENT(IN) :: DT,DY,PDTOP
!
      REAL,DIMENSION(KMS:KME),INTENT(IN) :: AETA1,AETA2,DETA1,DETA2
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: DX,HBM2,HBM3,PDSL
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: U,V,Z
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: CWM,Q,Q2
!
      LOGICAL,INTENT(IN) :: HYDRO
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: FF1=0.530
!
      INTEGER :: DOMDESC
!
!
      LOGICAL :: BOT,TOP
!
      INTEGER :: I,IRECV,J,JFP,JFQ,K,LAP,LLAP,MPI_COMM_COMP
!
      INTEGER,DIMENSION(ITS-5:ITE+5,JTS-5:JTE+5,KTS:KTE) :: IFPA,IFPF   &
     &                                                     ,IFQA,IFQF   &
     &                                                     ,JFPA,JFPF   &
     &                                                     ,JFQA,JFQF
!
      REAL :: ADDT,AFRP,CRIT,D2PQE,D2PQQ,D2PQW,DEP,DESTIJ,DQP,DQSTIJ    &
     &       ,DVOLP,DWP,DWSTIJ,DZA,DZB,E00,E0Q,E1X,E2IJ,E4P,ENH,EP,EP0  &
     &       ,ESTIJ,FPQ,HAFP,HAFQ,HBM2IJ,HM,PP,PPQ00,Q00,Q0Q            &
     &       ,Q1IJ,Q4P,QP,QP0,QSTIJ,RDY,RFACE,RFACQ,RFACW,RFC           &
     &       ,RFEIJ,RFQIJ,RFWIJ,RR,SLOPAC,SPP,SQP,SSA,SSB,SUMNE,SUMNQ   &
     &       ,SUMNW,SUMPE,SUMPQ,SUMPW,TTA,TTB,W00,W0Q,W1IJ,W4P,WP,WP0   &
     &       ,WSTIJ
!
      DOUBLE PRECISION,DIMENSION(6,KTS:KTE) :: GSUMS,XSUMS
!
      REAL,DIMENSION(KTS:KTE) :: AFR,DEL,DQL,DWL,E3,E4                  &
     &                          ,Q3,Q4,W3,W4
!
      REAL,DIMENSION(ITS-5:ITE+5,JTS-5:JTE+5) :: DARE,EMH
!
      REAL,DIMENSION(ITS-5:ITE+5,JTS-5:JTE+5,KTS:KTE) :: AFP,AFQ,DEST   &
     &                                                  ,DQST,DVOL,DWST &
     &                                                  ,E1,E2,Q1,W1
!-----------------------------------------------------------------------
      integer :: nunit,ier
      save nunit
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RDY=1./DY
      SLOPAC=SLOPHT*SQRT(2.)*0.5*50.
      CRIT=SLOPAC*REAL(IDTAD)*DT*RDY*1000.
!
      ADDT=REAL(IDTAD)*DT
      ENH=ADDT/(08.*DY)
!
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=max(jds+(0),jts-(3)),min(jde-(0),jte+(3))
      DO I=max(ids+(0),its-(2)),min(ide-(0),ite+(2))
        EMH (I,J)=ADDT/(08.*DX(I,J))
        DARE(I,J)=HBM3(I,J)*DX(I,J)*DY
        E1(I,J,KTE)=MAX(Q2(I,J,KTE)*0.5,EPSQ2)
        E2(I,J,KTE)=E1(I,J,KTE)
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(dza,dzb,e1x,fpq,hm,i,j,jfp,jfq,k,pp,qp,ssa,ssb,spp,sqp   &
!$omp&        ,tta,ttb)
!-----------------------------------------------------------------------
!
      vertical_1: DO K=KTS,KTE
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(0),jts-(3)),min(jde-(0),jte+(3))
        DO I=max(ids+(0),its-(2)),min(ide-(0),ite+(2))
          DVOL(I,J,K)=DARE(I,J)*(DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J))
          Q  (I,J,K)=MAX(Q  (I,J,K),EPSQ)
          CWM(I,J,K)=MAX(CWM(I,J,K),CLIMIT)
          Q1 (I,J,K)=Q  (I,J,K)
          W1 (I,J,K)=CWM(I,J,K)
        ENDDO
        ENDDO
!
        IF(K<KTE)THEN
          DO J=max(jds+(0),jts-(3)),min(jde-(0),jte+(3))
          DO I=max(ids+(0),its-(2)),min(ide-(0),ite+(2))
            E1X=(Q2(I,J,K+1)+Q2(I,J,K))*0.5
            E1(I,J,K)=MAX(E1X,EPSQ2)
            E2(I,J,K)=E1(I,J,K)
          ENDDO
          ENDDO
        ENDIF
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(2),jts-(1)),min(jde-(2),jte+(1))
        DO I=max(ids+(1),its-(1)),min(ide-(1),ite+(1))
!
          HM=HBM2(I,J)
          TTA=(U(I,J-1,K)+U(I+IHW(J),J,K)+U(I+IHE(J),J,K)+U(I,J+1,K))   &
     &        *EMH(I,J)*HM
          TTB=(V(I,J-1,K)+V(I+IHW(J),J,K)+V(I+IHE(J),J,K)+V(I,J+1,K))   &
     &        *ENH*HBM2(I,J)
!
          SPP=-TTA-TTB
          SQP= TTA-TTB
!
          IF(SPP<0.)THEN
            JFP=-1
          ELSE
            JFP=1
          ENDIF
          IF(SQP<0.)THEN
            JFQ=-1
          ELSE
            JFQ=1
          ENDIF
!
          IFPA(I,J,K)=IHE(J)+I+( JFP-1)/2
          IFQA(I,J,K)=IHE(J)+I+(-JFQ-1)/2
!
          JFPA(I,J,K)=J+JFP
          JFQA(I,J,K)=J+JFQ
!
          IFPF(I,J,K)=IHE(J)+I+(-JFP-1)/2
          IFQF(I,J,K)=IHE(J)+I+( JFQ-1)/2
!
          JFPF(I,J,K)=J-JFP
          JFQF(I,J,K)=J-JFQ
      if(i==111.and.j==438.and.k==1)then
      endif
!
!-----------------------------------------------------------------------
          IF(.NOT.HYDRO)THEN ! z currently not available for hydro=.true.
            DZA=(Z(IFPA(I,J,K),JFPA(I,J,K),K)-Z(I,J,K))*RDY
            DZB=(Z(IFQA(I,J,K),JFQA(I,J,K),K)-Z(I,J,K))*RDY
!
            IF(ABS(DZA)>SLOPAC)THEN
              SSA=DZA*SPP
              IF(SSA>CRIT)THEN
                SPP=0. !spp*.1
              ENDIF
            ENDIF
!
            IF(ABS(DZB)>SLOPAC)THEN
              SSB=DZB*SQP
              IF(SSB>CRIT)THEN
                SQP=0. !sqp*.1
              ENDIF
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!
          FPQ=SPP*SQP*0.25
          PP=ABS(SPP)
          QP=ABS(SQP)
!
          AFP(I,J,K)=(((FF4*PP+FF3)*PP+FF2)*PP+FF1)*PP
          AFQ(I,J,K)=(((FF4*QP+FF3)*QP+FF2)*QP+FF1)*QP
!
          Q1(I,J,K)=(Q  (IFPA(I,J,K),JFPA(I,J,K),K)-Q  (I,J,K))*PP        &
       &           +(Q  (IFQA(I,J,K),JFQA(I,J,K),K)-Q  (I,J,K))*QP        &
       &           +(Q  (I,J-2,K)+Q  (I,J+2,K)                            &
       &            -Q  (I-1,J,K)-Q  (I+1,J,K))*FPQ                       &
       &           +Q(I,J,K)
!
          W1(I,J,K)=(CWM(IFPA(I,J,K),JFPA(I,J,K),K)-CWM(I,J,K))*PP        &
       &           +(CWM(IFQA(I,J,K),JFQA(I,J,K),K)-CWM(I,J,K))*QP        &
       &           +(CWM(I,J-2,K)+CWM(I,J+2,K)                            &
       &            -CWM(I-1,J,K)-CWM(I+1,J,K))*FPQ                       &
       &           +CWM(I,J,K)
!
          E2(I,J,K)=(E1 (IFPA(I,J,K),JFPA(I,J,K),K)-E1 (I,J,K))*PP        &
       &           +(E1 (IFQA(I,J,K),JFQA(I,J,K),K)-E1 (I,J,K))*QP        &
       &           +(E1 (I,J-2,K)+E1 (I,J+2,K)                            &
       &            -E1 (I-1,J,K)-E1 (I+1,J,K))*FPQ                       &
       &           +E1(I,J,K)
!
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO vertical_1
!
!-----------------------------------------------------------------------
!***  ANTI-FILTERING STEP
!-----------------------------------------------------------------------
!
      DO K=KTS,KTE
        XSUMS(1,K)=0.
        XSUMS(2,K)=0.
        XSUMS(3,K)=0.
        XSUMS(4,K)=0.
        XSUMS(5,K)=0.
        XSUMS(6,K)=0.
      ENDDO
!-----------------------------------------------------------------------
!
!***  ANTI-FILTERING LIMITERS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(d2pqe,d2pqq,d2pqw,destij,dqstij,dvolp,dwstij             &
!$omp&        ,e00,e0q,e2ij,ep0,estij,hafp,hafq,i,j,k                   &
!$omp&        ,q00,q0q,q1ij,qp0,qstij,w00,w0q,w1ij,wp0,wstij)
!-----------------------------------------------------------------------
!
      vertical_2: DO K=KTS,KTE
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
!
          DVOLP=DVOL(I,J,K)
          Q1IJ =Q1(I,J,K)
          W1IJ =W1(I,J,K)
          E2IJ =E2(I,J,K)
!
          HAFP=AFP(I,J,K)
          HAFQ=AFQ(I,J,K)
!
          D2PQQ=(Q1(IFPA(I,J,K),JFPA(I,J,K),K)-Q1IJ                     &
     &          -Q1IJ+Q1(IFPF(I,J,K),JFPF(I,J,K),K))                    &
     &          *HAFP                                                   &
     &         +(Q1(IFQA(I,J,K),JFQA(I,J,K),K)-Q1IJ                     &
     &          -Q1IJ+Q1(IFQF(I,J,K),JFQF(I,J,K),K))                    &
     &          *HAFQ
!
          D2PQW=(W1(IFPA(I,J,K),JFPA(I,J,K),K)-W1IJ                     &
     &          -W1IJ+W1(IFPF(I,J,K),JFPF(I,J,K),K))                    &
     &          *HAFP                                                   &
     &         +(W1(IFQA(I,J,K),JFQA(I,J,K),K)-W1IJ                     &
     &          -W1IJ+W1(IFQF(I,J,K),JFQF(I,J,K),K))                    &
     &          *HAFQ
!
          D2PQE=(E2(IFPA(I,J,K),JFPA(I,J,K),K)-E2IJ                     &
     &          -E2IJ+E2(IFPF(I,J,K),JFPF(I,J,K),K))                    &
     &          *HAFP                                                   &
     &         +(E2(IFQA(I,J,K),JFQA(I,J,K),K)-E2IJ                     &
     &          -E2IJ+E2(IFQF(I,J,K),JFQF(I,J,K),K))                    &
     &          *HAFQ
!
          QSTIJ=Q1IJ-D2PQQ
          WSTIJ=W1IJ-D2PQW
          ESTIJ=E2IJ-D2PQE
!
          Q00=Q  (I          ,J          ,K)
          QP0=Q  (IFPA(I,J,K),JFPA(I,J,K),K)
          Q0Q=Q  (IFQA(I,J,K),JFQA(I,J,K),K)
!
          W00=CWM(I          ,J          ,K)
          WP0=CWM(IFPA(I,J,K),JFPA(I,J,K),K)
          W0Q=CWM(IFQA(I,J,K),JFQA(I,J,K),K)
!
          E00=E1 (I          ,J          ,K)
          EP0=E1 (IFPA(I,J,K),JFPA(I,J,K),K)
          E0Q=E1 (IFQA(I,J,K),JFQA(I,J,K),K)
!
          QSTIJ=MAX(QSTIJ,MIN(Q00,QP0,Q0Q))
          QSTIJ=MIN(QSTIJ,MAX(Q00,QP0,Q0Q))
          WSTIJ=MAX(WSTIJ,MIN(W00,WP0,W0Q))
          WSTIJ=MIN(WSTIJ,MAX(W00,WP0,W0Q))
          ESTIJ=MAX(ESTIJ,MIN(E00,EP0,E0Q))
          ESTIJ=MIN(ESTIJ,MAX(E00,EP0,E0Q))
!
          DQSTIJ=QSTIJ-Q(I,J,K)
          DWSTIJ=WSTIJ-CWM(I,J,K)
          DESTIJ=ESTIJ-E1(I,J,K)
!
          DQST(I,J,K)=DQSTIJ
          DWST(I,J,K)=DWSTIJ
          DEST(I,J,K)=DESTIJ
!
          DQSTIJ=DQSTIJ*DVOLP
          DWSTIJ=DWSTIJ*DVOLP
          DESTIJ=DESTIJ*DVOLP
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
          IF(DQSTIJ>0.)THEN
            XSUMS(1,K)=XSUMS(1,K)+DQSTIJ
          ELSE
            XSUMS(2,K)=XSUMS(2,K)+DQSTIJ
          ENDIF
!
          IF(DWSTIJ>0.)THEN
            XSUMS(3,K)=XSUMS(3,K)+DWSTIJ
          ELSE
            XSUMS(4,K)=XSUMS(4,K)+DWSTIJ
          ENDIF
!
          IF(DESTIJ>0.)THEN
            XSUMS(5,K)=XSUMS(5,K)+DESTIJ
          ELSE
            XSUMS(6,K)=XSUMS(6,K)+DESTIJ
          ENDIF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO vertical_2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  GLOBAL REDUCTION
!-----------------------------------------------------------------------
!
      CALL WRF_GET_DM_COMMUNICATOR(MPI_COMM_COMP)
      CALL MPI_ALLREDUCE(XSUMS,GSUMS,6*(KTE-KTS+1)                      &
     &                  ,MPI_DOUBLE_PRECISION,MPI_SUM                   &
     &                  ,MPI_COMM_COMP,IRECV)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  END OF GLOBAL REDUCTION
!-----------------------------------------------------------------------
!
!     if(mype==0)then
!!!     if(ntsd==0)then
!!!       call int_get_fresh_handle(nunit)
!!!       close(nunit)
!         nunit=56
!!!       open(unit=nunit,file=gsums,form=unformatted,iostat=ier)
!!!     endif
!     endif
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(destij,dqstij,dwstij,i,j,k,rface,rfacq,rfacw             &
!$omp&        ,rfeij,rfqij,rfwij,sumne,sumnq,sumnw,sumpe,sumpq,sumpw)
!-----------------------------------------------------------------------
!
      vertical_3: DO K=KTS,KTE
!
!-----------------------------------------------------------------------
!       if(mype==0)then
!         write(nunit)(gsums(i,k),i=1,6)
!       endif
!!!     read(nunit)(gsums(i,k),i=1,6)
!-----------------------------------------------------------------------
!
        SUMPQ=GSUMS(1,K)
        SUMNQ=GSUMS(2,K)
        SUMPW=GSUMS(3,K)
        SUMNW=GSUMS(4,K)
        SUMPE=GSUMS(5,K)
        SUMNE=GSUMS(6,K)
!
!-----------------------------------------------------------------------
!***  FIRST MOMENT CONSERVING FACTOR
!-----------------------------------------------------------------------
!
        IF(SUMPQ>1.)THEN
          RFACQ=-SUMNQ/SUMPQ
        ELSE
          RFACQ=1.
        ENDIF
!
        IF(SUMPW>1.)THEN
          RFACW=-SUMNW/SUMPW
        ELSE
          RFACW=1.
        ENDIF
!
        IF(SUMPE>1.)THEN
          RFACE=-SUMNE/SUMPE
        ELSE
          RFACE=1.
        ENDIF
!
        IF(RFACQ<CONSERVE_MIN.OR.RFACQ>CONSERVE_MAX)RFACQ=1.
        IF(RFACW<CONSERVE_MIN.OR.RFACW>CONSERVE_MAX)RFACW=1.
        IF(RFACE<CONSERVE_MIN.OR.RFACE>CONSERVE_MAX)RFACE=1.
!
!-----------------------------------------------------------------------
!       if(mype==0.and.ntsd==181)close(nunit)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  IMPOSE CONSERVATION ON ANTI-FILTERING
!-----------------------------------------------------------------------
!
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IF(RFACQ<1.)THEN
            DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
              DQSTIJ=DQST(I,J,K)
              RFQIJ=HBM2(I,J)*(RFACQ-1.)+1.
              IF(DQSTIJ>=0.)DQSTIJ=DQSTIJ*RFQIJ
              Q(I,J,K)=Q(I,J,K)+DQSTIJ
            ENDDO
          ELSE
            DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
              DQSTIJ=DQST(I,J,K)
              RFQIJ=HBM2(I,J)*(RFACQ-1.)+1.
              IF(DQSTIJ<0.)DQSTIJ=DQSTIJ/RFQIJ
              Q(I,J,K)=Q(I,J,K)+DQSTIJ
            ENDDO
          ENDIF
        ENDDO
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IF(RFACW<1.)THEN
            DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
              DWSTIJ=DWST(I,J,K)
              RFWIJ=HBM2(I,J)*(RFACW-1.)+1.
              IF(DWSTIJ>=0.)DWSTIJ=DWSTIJ*RFWIJ
              CWM(I,J,K)=CWM(I,J,K)+DWSTIJ
            ENDDO
          ELSE
            DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
              DWSTIJ=DWST(I,J,K)
              RFWIJ=HBM2(I,J)*(RFACW-1.)+1.
              IF(DWSTIJ<0.)DWSTIJ=DWSTIJ/RFWIJ
              CWM(I,J,K)=CWM(I,J,K)+DWSTIJ
            ENDDO
          ENDIF
        ENDDO
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IF(RFACE<1.)THEN
            DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
              DESTIJ=DEST(I,J,K)
              RFEIJ=HBM2(I,J)*(RFACE-1.)+1.
              IF(DESTIJ>=0.)DESTIJ=DESTIJ*RFEIJ
              E1(I,J,K)=E1(I,J,K)+DESTIJ
            ENDDO
          ELSE
            DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
              DESTIJ=DEST(I,J,K)
              RFEIJ=HBM2(I,J)*(RFACE-1.)+1.
              IF(DESTIJ<0.)DESTIJ=DESTIJ/RFEIJ
              E1(I,J,K)=E1(I,J,K)+DESTIJ
            ENDDO
          ENDIF
        ENDDO
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
        DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
          Q  (I,J,K)=MAX(Q  (I,J,K),EPSQ)
          CWM(I,J,K)=MAX(CWM(I,J,K),CLIMIT)
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO vertical_3
!
!-----------------------------------------------------------------------
!
!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
      DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
        Q2(I,J,KTE)=MAX(E1(I,J,KTE)+E1(I,J,KTE)-EPSQ2,EPSQ2)
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      DO K=KTE-1,KTS+1,-1
!$omp parallel do                                                       &
!$omp& private(i,j)
        DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
        DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
          IF(K>KTS)THEN
            Q2(I,J,K)=MAX(E1(I,J,K)+E1(I,J,K)-Q2(I,J,K+1),EPSQ2)
          ELSE
            Q2(I,J,K)=Q2(I,J,K+1)
          ENDIF
        ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!
      END SUBROUTINE HAD2
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!!!! THE FOLLOWING  _SCAL ROUTINES ARE NOT IN IJK YET  !!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
! New routines added by Georg Grell to handle advection more like ARW
! core.  Instead of VAD2/HAD2 that advect TKE, specific humidity, and
! condensed water species all in one routine, we call VAD2/HAD2_SCAL
! with multidimensioned arrays to advect each variable.  For purposes
! here, solve_nmm.F calls this routine once for TKE, then again for
! all the species held in the moist array (qv, qc, qi, qr, qs, qg),
! then call again for number concentrations held in scalar array (qni).
! The dummy argument lstart is the starting index of the multidimensioned
! array for starting the advection since the 1st index of moist and
! scalar are actually empty placeholders (and the 2nd element is vapor,
! then qc, etc.)  When calling with single 3D array (like TKE), just
! set NUM_SCAL=1 and lstart=1.  The variable to advect is called SCAL
! herein.
!***********************************************************************
      SUBROUTINE VAD2_SCAL(NTSD,DT,IDTAD,DX,DY                          &
     &               ,AETA1,AETA2,DETA1,DETA2,PDSL,PDTOP                &
     &               ,HBM2                                              &
     &               ,SCAL,PETDT                                        &
     &               ,N_IUP_H,N_IUP_V                                   &
     &               ,N_IUP_ADH,N_IUP_ADV                               &
     &               ,IUP_H,IUP_V,IUP_ADH,IUP_ADV                       &
     &               ,IHE,IHW,IVE,IVW                                   &
     &               ,NUM_SCAL,LSTART                                   &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &               ,IMS,IME,JMS,JME,KMS,KME                           &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    VAD2_SCAL   VERTICAL ADVECTION OF SCALARS
!
!   PRGRMMR: JANJIC          ORG: W/NP22     DATE: 96-07-19
!            GRELL,PECKHAM   ORG: NOAA/FSL   DATE: 05-02-03
!     
! ABSTRACT:          
!     VAD2_SCAL CALCULATES THE CONTRIBUTION OF THE VERTICAL ADVECTION   
!     TO THE TENDENCIES OF SCALAR SUBSTANCES AND THEN UPDATES           
!     THOSE VARIABLES.  AN ANTI-FILTERING TECHNIQUE IS USED.            
!    
! PROGRAM HISTORY LOG:
!   96-07-19  JANJIC           - ORIGINATOR
!   05-02-03  GRELL,PECKHAM    - MODIFIED FOR SCALARS                   
!    
! USAGE: CALL VAD2_SCAL FROM SUBROUTINE SOLVE_NMM                       
!   INPUT ARGUMENT LIST:
!
!   OUTPUT ARGUMENT LIST
!                
!   OUTPUT FILES:
!       NONE
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
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
                           ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: LSTART,NUM_SCAL
!
      INTEGER,DIMENSION(JMS:JME),INTENT(IN) :: IHE,IHW,IVE,IVW
      INTEGER,DIMENSION(JMS:JME),INTENT(IN) :: N_IUP_H,N_IUP_V          &
     &                                        ,N_IUP_ADH,N_IUP_ADV
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: IUP_H,IUP_V      &
     &                                                ,IUP_ADH,IUP_ADV
!
      INTEGER,INTENT(IN) :: IDTAD,NTSD
!
      REAL,INTENT(IN) :: DT,DY,PDTOP
!
      REAL,DIMENSION(KMS:KME),INTENT(IN) :: AETA1,AETA2,DETA1,DETA2
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: DX,HBM2,PDSL
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: PETDT
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,1:NUM_SCAL)               &
                                                 ,INTENT(INOUT) :: SCAL
!
!----------------------------------------------------------------------
!***  LOCAL VARIABLES
!----------------------------------------------------------------------
!
      REAL,PARAMETER :: FF1=0.500
!
      LOGICAL,SAVE :: TRADITIONAL=.TRUE.
!
      INTEGER :: I,IRECV,J,JFP,JFQ,K,L,LAP,LLAP
!
      INTEGER,DIMENSION(KTS:KTE) :: LA
!
      REAL*8 :: ADDT,AFRP,D2PQQ,DETAP,DPDN,DPUP,DQP                     &
     &         ,HADDT,HBM2IJ                                            &
     &         ,Q00,Q4P,QP,QP0                                          &
     &         ,RFACQK,RFC,RR                                           &
     &         ,SUMNQ,SUMPQ
!
      REAL :: SFACQK
!
      REAL,DIMENSION(KTS:KTE) :: AFR,DEL,DQL,DWL,E3,E4,PETDTK           &
     &                          ,RFACE,RFACQ,RFACW,Q3,Q4,W3,W4
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      ADDT=REAL(IDTAD)*DT
!
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(afr,afrp,d2pqq,detap,dpdn,dpup                           &
!$omp&        ,dql,dqp,haddt,i,j,k                                      &
!$omp&        ,la,lap,llap,petdtk,q00,q3,q4,q4p,qp,qp0,rfacqk           &
!$omp&        ,rfc,rr,sfacqk,sumnq,sumpq)
!-----------------------------------------------------------------------
!
      scalar_loop: DO L=LSTART,NUM_SCAL
!
      main_integration: DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
!
!-----------------------------------------------------------------------
!
        main_iloop: DO I=max(ids+(1),its-(1)),min(ide-(1),ite+(1))
!
!-----------------------------------------------------------------------
!
          DO K=KTS,KTE
            Q3(K)=SCAL(I,J,K,L)
            Q4(K)=Q3(K)
          ENDDO
!
          IF(TRADITIONAL)THEN
            PETDTK(KTE)=PETDT(I,J,KTE-1)*0.5
!
            DO K=KTE-1,KTS+1,-1
              PETDTK(K)=(PETDT(I,J,K)+PETDT(I,J,K-1))*0.5
            ENDDO
!
            PETDTK(KTS)=PETDT(I,J,KTS)*0.5
!
          ELSE
!
!-----------------------------------------------------------------------
!***    PERFORM HORIZONTAL AVERAGING OF VERTICAL VELOCITY
!-----------------------------------------------------------------------
!
            PETDTK(KTE)=(PETDT(I+IHW(J-1),J-1,KTE-1)                    &
     &                  +PETDT(I+IHE(J-1),J-1,KTE-1)                    &
     &                  +PETDT(I+IHW(J+1),J+1,KTE-1)                    &
     &                  +PETDT(I+IHE(J+1),J+1,KTE-1)                    &
     &                  +PETDT(I,J,KTE-1)*4.        )*0.0625
!
            DO K=KTE-1,KTS+1,-1
              PETDTK(K)=(PETDT(I+IHW(J-1),J-1,K-1)                      &
                        +PETDT(I+IHE(J-1),J-1,K-1)                      &
     &                  +PETDT(I+IHW(J+1),J+1,K-1)                      &
     &                  +PETDT(I+IHE(J+1),J+1,K-1)                      &
     &                  +PETDT(I+IHW(J-1),J-1,K  )                      &
     &                  +PETDT(I+IHE(J-1),J-1,K  )                      &
     &                  +PETDT(I+IHW(J+1),J+1,K  )                      &
     &                  +PETDT(I+IHE(J+1),J+1,K  )                      &
     &                  +(PETDT(I,J,K-1)+PETDT(I,J,K))*4.               &
     &                                                   )*0.0625
            ENDDO
!
            PETDTK(KTS)=(PETDT(I+IHW(J-1),J-1,KTS)                      &
     &                  +PETDT(I+IHE(J-1),J-1,KTS)                      &
     &                  +PETDT(I+IHW(J+1),J+1,KTS)                      &
     &                  +PETDT(I+IHE(J+1),J+1,KTS)                      &
     &                  +PETDT(I,J,KTS)*4.        )*0.0625
 
          ENDIF
!
!-----------------------------------------------------------------------
!
          HADDT=-ADDT*HBM2(I,J)
!
          DO K=KTE,KTS,-1
            RR=PETDTK(K)*HADDT
!
            IF(RR<0.)THEN
              LAP=1
            ELSE
              LAP=-1
            ENDIF
!
            LA(K)=LAP
            LLAP=K+LAP
!
            IF(LLAP>KTS-1.AND.LLAP<KTE+1)THEN
              RR=ABS(RR/((AETA1(LLAP)-AETA1(K))*PDTOP                   &
     &                  +(AETA2(LLAP)-AETA2(K))*PDSL(I,J)))
              IF(RR>0.9)RR=0.9
!
              AFR(K)=(((FF4*RR+FF3)*RR+FF2)*RR+FF1)*RR
              DQP=(Q3(LLAP)-Q3(K))*RR
              DQL(K)=DQP
            ELSE
              RR=0.
              AFR(K)=0.
              DQL(K)=0.
            ENDIF
          ENDDO
!
!-----------------------------------------------------------------------
!
          IF(LA(KTE-1)>0)THEN
            RFC=(DETA1(KTE-1)*PDTOP+DETA2(KTE-1)*PDSL(I,J))             &
     &         /(DETA1(KTE  )*PDTOP+DETA2(KTE  )*PDSL(I,J))
            DQL(KTE)=-DQL(KTE-1)*RFC
          ENDIF
!
          IF(LA(KTS+1)<0)THEN
            RFC=(DETA1(KTS+1)*PDTOP+DETA2(KTS+1)*PDSL(I,J))           &
     &         /(DETA1(KTS  )*PDTOP+DETA2(KTS  )*PDSL(I,J))
            DQL(KTS)=-DQL(KTS+1)*RFC
          ENDIF
!
          DO K=KTS,KTE
            Q4(K)=Q3(K)+DQL(K)
          ENDDO
!
!-----------------------------------------------------------------------
!***  ANTI-FILTERING STEP
!-----------------------------------------------------------------------
!
          SUMPQ=0.
          SUMNQ=0.
!
!***  ANTI-FILTERING LIMITERS
!
          antifilter: DO K=KTE-1,KTS+1,-1
!
            DETAP=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
!
            Q4P=Q4(K)
!
            LAP=LA(K)
!
            DPDN=(AETA1(K+LAP)-AETA1(K))*PDTOP                          &
     &          +(AETA2(K+LAP)-AETA2(K))*PDSL(I,J)
            DPUP=(AETA1(K)-AETA1(K-LAP))*PDTOP                          &
     &          +(AETA2(K)-AETA2(K-LAP))*PDSL(I,J)
!
            AFRP=2.*AFR(K)*DPDN*DPDN/(DPDN+DPUP)
            D2PQQ=((Q4(K+LAP)-Q4P)/DPDN                                 &
     &            -(Q4P-Q4(K-LAP))/DPUP)*AFRP
!
            QP=Q4P-D2PQQ
!
            Q00=Q3(K)
            QP0=Q3(K+LAP)
!
            QP=MAX(QP,MIN(Q00,QP0))
            QP=MIN(QP,MAX(Q00,QP0))
!
            DQP=QP-Q00
!
            DQL(K)=DQP
!
          ENDDO antifilter
!
!-----------------------------------------------------------------------
!
          IF(LA(KTE-1)>0)THEN
            RFC=(DETA1(KTE-1)*PDTOP+DETA2(KTE-1)*PDSL(I,J))             &
     &         /(DETA1(KTE  )*PDTOP+DETA2(KTE  )*PDSL(I,J))
            DQL(KTE)=-DQL(KTE-1)*RFC+DQL(KTE)
          ENDIF
!
          IF(LA(KTS+1)<0)THEN
            RFC=(DETA1(KTS+1)*PDTOP+DETA2(KTS+1)*PDSL(I,J))             &
     &         /(DETA1(KTS  )*PDTOP+DETA2(KTS  )*PDSL(I,J))
            DQL(KTS)=-DQL(KTS+1)*RFC+DQL(KTS)
          ENDIF
!
          DO K=KTS,KTE
            DETAP=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
            DQP=DQL(K)*DETAP
!
            IF(DQP>0.)THEN
              SUMPQ=SUMPQ+DQP
            ELSE
              SUMNQ=SUMNQ+DQP
            ENDIF
          ENDDO
!
!-----------------------------------------------------------------------
!***  FIRST MOMENT CONSERVING FACTOR
!-----------------------------------------------------------------------
!
          IF(SUMPQ>1.E-9)THEN
            SFACQK=-SUMNQ/SUMPQ
          ELSE
            SFACQK=1.
          ENDIF
!
          IF(SFACQK<CONSERVE_MIN.OR.SFACQK>CONSERVE_MAX)SFACQK=1.
!
          RFACQK=1./SFACQK
!
!-----------------------------------------------------------------------
!***  IMPOSE CONSERVATION ON ANTI-FILTERING
!-----------------------------------------------------------------------
!
          DO K=KTE,KTS,-1
            DQP=DQL(K)
            IF(SFACQK>=1.)THEN
              IF(DQP<0.)DQP=DQP*RFACQK
            ELSE
              IF(DQP>0.)DQP=DQP*SFACQK
            ENDIF
            SCAL(I,J,K,L)=Q3(K)+DQP
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDDO main_iloop
!
!-----------------------------------------------------------------------
!
      ENDDO main_integration
!
!-----------------------------------------------------------------------
!
      ENDDO scalar_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE VAD2_SCAL
!
!-----------------------------------------------------------------------
!         
!***********************************************************************
      SUBROUTINE HAD2_SCAL(                                             &
     &                DOMDESC ,                                         &
     &                NTSD,DT,IDTAD,DX,DY                               &
     &               ,AETA1,AETA2,DETA1,DETA2,PDSL,PDTOP                &
     &               ,HBM2,HBM3                                         &
     &               ,SCAL,U,V,Z,HYDRO                                  &
     &               ,N_IUP_H,N_IUP_V                                   &
     &               ,N_IUP_ADH,N_IUP_ADV                               &
     &               ,IUP_H,IUP_V,IUP_ADH,IUP_ADV                       &
     &               ,IHE,IHW,IVE,IVW                                   &
     &               ,NUM_SCAL,LSTART                                   &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &               ,IMS,IME,JMS,JME,KMS,KME                           &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE)
!***********************************************************************
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
! SUBPROGRAM:    HAD2_SCAL   HORIZONTAL ADVECTION OF SCALAR
!   PRGRMMR: JANJIC          ORG: W/NP22     DATE: 96-07-19
!            GRELL,PECKHAM   ORG: NOAA/FSL   DATE: 05-02-03
!
! ABSTRACT:
!     HAD2_SCAL CALCULATES THE CONTRIBUTION OF THE HORIZONTAL ADVECTION
!     TO THE TENDENCIES OF SCALAR SUBSTANCES AND THEN
!     UPDATES THOSE VARIABLES.  AN ANTI-FILTERING TECHNIQUE IS USED.
!
! PROGRAM HISTORY LOG:
!   96-07-19  JANJIC           - ORIGINATOR
!   05-01-03  GRELL,PECKKHAM   - MODIFIED FOR SCALAR
!
! USAGE: CALL HAD2_SCAL FROM SUBROUTINE SOLVE_NMM
!   INPUT ARGUMENT LIST:
!
!   OUTPUT ARGUMENT LIST
!
!   OUTPUT FILES:
!       NONE
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
!    
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,DIMENSION(JMS:JME),INTENT(IN) :: IHE,IHW,IVE,IVW
      INTEGER,DIMENSION(JMS:JME),INTENT(IN) :: N_IUP_H,N_IUP_V          &
     &                                        ,N_IUP_ADH,N_IUP_ADV
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: IUP_H,IUP_V      &
     &                                                ,IUP_ADH,IUP_ADV
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IDTAD,LSTART,NTSD,NUM_SCAL
!
      REAL,INTENT(IN) :: DT,DY,PDTOP
!
      REAL,DIMENSION(KMS:KME),INTENT(IN) :: AETA1,AETA2,DETA1,DETA2
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: DX,HBM2,HBM3,PDSL
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: U,V,Z
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,1:NUM_SCAL)                &
                                                  ,INTENT(INOUT) :: SCAL
!
      LOGICAL,INTENT(IN) :: HYDRO
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: FF1=0.530
!
      INTEGER :: DOMDESC
!
!
      LOGICAL :: BOT,TOP
!
      INTEGER :: I,IRECV,J,JFP,JFQ,K,L,LAP,LLAP,MPI_COMM_COMP
!
      INTEGER,DIMENSION(ITS-5:ITE+5,JTS-5:JTE+5,KTS:KTE) :: IFPA,IFPF   &
     &                                                     ,IFQA,IFQF   &
     &                                                     ,JFPA,JFPF   &
     &                                                     ,JFQA,JFQF
!
      REAL :: ADDT,AFRP,CRIT,D2PQE,D2PQQ,D2PQW,DEP,DESTIJ,DQP,DQSTIJ    &
     &       ,DVOLP,DWP,DWSTIJ,DZA,DZB,E00,E0Q,E1X,E2IJ,E4P,ENH,EP,EP0  &
     &       ,ESTIJ,FPQ,HAFP,HAFQ,HBM2IJ,HM,PP,PPQ00,Q00,Q0Q            &
     &       ,Q1IJ,Q4P,QP,QP0,QSTIJ,RDY,RFACE,RFACQ,RFACW,RFC           &
     &       ,RFEIJ,RFQIJ,RFWIJ,RR,SLOPAC,SPP,SQP,SSA,SSB,SUMNQ,SUMPQ   &
     &       ,TTA,TTB,W00,W0Q,W1IJ,W4P,WP,WP0,WSTIJ
!
      DOUBLE PRECISION,DIMENSION(2,KTS:KTE) :: GSUMS,XSUMS
!
      REAL,DIMENSION(KTS:KTE) :: AFR,DEL,DQL,DWL,Q3,Q4
!
      REAL,DIMENSION(ITS-5:ITE+5,JTS-5:JTE+5) :: DARE,EMH
!
      REAL,DIMENSION(ITS-5:ITE+5,JTS-5:JTE+5,KTS:KTE) :: AFP,AFQ,DEST   &
     &                                                  ,DQST,DVOL,DWST &
     &                                                  ,Q1
!
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME) :: Q
!
!-----------------------------------------------------------------------
      integer :: nunit,ier
      save nunit
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RDY=1./DY
      SLOPAC=SLOPHT*SQRT(2.)*0.5*50.
      CRIT=SLOPAC*REAL(IDTAD)*DT*RDY*1000.
!
      ADDT=REAL(IDTAD)*DT
      ENH=ADDT/(08.*DY)
!
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=max(jds+(0),jts-(3)),min(jde-(0),jte+(3))
      DO I=max(ids+(0),its-(2)),min(ide-(0),ite+(2))
        EMH (I,J)=ADDT/(08.*DX(I,J))
        DARE(I,J)=HBM3(I,J)*DX(I,J)*DY
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
!
      scalar_loop: DO L=LSTART,NUM_SCAL
!
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(dza,dzb,e1x,fpq,hm,i,j,jfp,jfq,k,pp,qp,ssa,ssb,spp,sqp   &
!$omp&        ,tta,ttb)
!-----------------------------------------------------------------------
!
      vertical_1: DO K=KTS,KTE
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(0),jts-(3)),min(jde-(0),jte+(3))
        DO I=max(ids+(0),its-(2)),min(ide-(0),ite+(2))
          DVOL(I,J,K)=DARE(I,J)*(DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J))
          Q (I,J,K)=SCAL(I,J,K,L)
          Q1(I,J,K)=Q(I,J,K)
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(2),jts-(1)),min(jde-(2),jte+(1))
        DO I=max(ids+(1),its-(1)),min(ide-(1),ite+(1))
!
          HM=HBM2(I,J)
          TTA=(U(I,J-1,K)+U(I+IHW(J),J,K)+U(I+IHE(J),J,K)+U(I,J+1,K))   &
     &        *EMH(I,J)*HM
          TTB=(V(I,J-1,K)+V(I+IHW(J),J,K)+V(I+IHE(J),J,K)+V(I,J+1,K))   &
     &        *ENH*HBM2(I,J)
!
          SPP=-TTA-TTB
          SQP= TTA-TTB
!
          IF(SPP<0.)THEN
            JFP=-1
          ELSE
            JFP=1
          ENDIF
          IF(SQP<0.)THEN
            JFQ=-1
          ELSE
            JFQ=1
          ENDIF
!
          IFPA(I,J,K)=IHE(J)+I+( JFP-1)/2
          IFQA(I,J,K)=IHE(J)+I+(-JFQ-1)/2
!
          JFPA(I,J,K)=J+JFP
          JFQA(I,J,K)=J+JFQ
!
          IFPF(I,J,K)=IHE(J)+I+(-JFP-1)/2
          IFQF(I,J,K)=IHE(J)+I+( JFQ-1)/2
!
          JFPF(I,J,K)=J-JFP
          JFQF(I,J,K)=J-JFQ
!
!-----------------------------------------------------------------------
          IF(.NOT.HYDRO)THEN ! z currently not available for hydro=.true.
            DZA=(Z(IFPA(I,J,K),JFPA(I,J,K),K)-Z(I,J,K))*RDY
            DZB=(Z(IFQA(I,J,K),JFQA(I,J,K),K)-Z(I,J,K))*RDY
!
            IF(ABS(DZA)>SLOPAC)THEN
              SSA=DZA*SPP
              IF(SSA>CRIT)THEN
                SPP=0. !spp*.1
              ENDIF
            ENDIF
!
            IF(ABS(DZB)>SLOPAC)THEN
              SSB=DZB*SQP
              IF(SSB>CRIT)THEN
                SQP=0. !sqp*.1
              ENDIF
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!
          FPQ=SPP*SQP*0.25
          PP=ABS(SPP)
          QP=ABS(SQP)
!
          AFP(I,J,K)=(((FF4*PP+FF3)*PP+FF2)*PP+FF1)*PP
          AFQ(I,J,K)=(((FF4*QP+FF3)*QP+FF2)*QP+FF1)*QP
!
          Q1(I,J,K)=(Q  (IFPA(I,J,K),JFPA(I,J,K),K)-Q  (I,J,K))*PP        &
       &           +(Q  (IFQA(I,J,K),JFQA(I,J,K),K)-Q  (I,J,K))*QP        &
       &           +(Q  (I,J-2,K)+Q  (I,J+2,K)                            &
       &            -Q  (I-1,J,K)-Q  (I+1,J,K))*FPQ                       &
       &           +Q(I,J,K)
!
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO vertical_1
!
!-----------------------------------------------------------------------
!***  ANTI-FILTERING STEP
!-----------------------------------------------------------------------
!
      DO K=KTS,KTE
        XSUMS(1,K)=0.
        XSUMS(2,K)=0.
      ENDDO
!
!-----------------------------------------------------------------------
!
!***  ANTI-FILTERING LIMITERS
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(d2pqe,d2pqq,d2pqw,destij,dqstij,dvolp,dwstij             &
!$omp&        ,e00,e0q,ep0,estij,hafp,hafq,i,j,k                        &
!$omp&        ,q00,q0q,q1ij,qp0,qstij,w00,w0q,wp0,wstij)
!-----------------------------------------------------------------------
!
      vertical_2: DO K=KTS,KTE
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
!
          DVOLP=DVOL(I,J,K)
          Q1IJ =Q1(I,J,K)
!
          HAFP=AFP(I,J,K)
          HAFQ=AFQ(I,J,K)
!
          D2PQQ=(Q1(IFPA(I,J,K),JFPA(I,J,K),K)-Q1IJ                     &
     &          -Q1IJ+Q1(IFPF(I,J,K),JFPF(I,J,K),K))                    &
     &          *HAFP                                                   &
     &         +(Q1(IFQA(I,J,K),JFQA(I,J,K),K)-Q1IJ                     &
     &          -Q1IJ+Q1(IFQF(I,J,K),JFQF(I,J,K),K))                    &
     &          *HAFQ
!
          QSTIJ=Q1IJ-D2PQQ
!
          Q00=Q  (I          ,J          ,K)
          QP0=Q  (IFPA(I,J,K),JFPA(I,J,K),K)
          Q0Q=Q  (IFQA(I,J,K),JFQA(I,J,K),K)
!
          QSTIJ=MAX(QSTIJ,MIN(Q00,QP0,Q0Q))
          QSTIJ=MIN(QSTIJ,MAX(Q00,QP0,Q0Q))
!
          DQSTIJ=QSTIJ-Q(I,J,K)
!
          DQST(I,J,K)=DQSTIJ
!
          DQSTIJ=DQSTIJ*DVOLP
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
          IF(DQSTIJ>0.)THEN
            XSUMS(1,K)=XSUMS(1,K)+DQSTIJ
          ELSE
            XSUMS(2,K)=XSUMS(2,K)+DQSTIJ
          ENDIF
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO vertical_2
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  GLOBAL REDUCTION
!-----------------------------------------------------------------------
!
      CALL WRF_GET_DM_COMMUNICATOR(MPI_COMM_COMP)
      CALL MPI_ALLREDUCE(XSUMS,GSUMS,2*(KTE-KTS+1)                      &
     &                  ,MPI_DOUBLE_PRECISION,MPI_SUM                   &
     &                  ,MPI_COMM_COMP,IRECV)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  END OF GLOBAL REDUCTION
!-----------------------------------------------------------------------
!
!     if(mype==0)then
!!!     if(ntsd==0)then
!!!       call int_get_fresh_handle(nunit)
!!!       close(nunit)
!         nunit=56
!!!       open(unit=nunit,file=gsums,form=unformatted,iostat=ier)
!!!     endif
!     endif
!-----------------------------------------------------------------------
!$omp parallel do                                                       &
!$omp& private(destij,dqstij,dwstij,i,j,k,rface,rfacq,rfacw             &
!$omp&        ,rfeij,rfqij,rfwij,sumne,sumnq,sumnw,sumpe,sumpq,sumpw)
!-----------------------------------------------------------------------
!
      vertical_3: DO K=KTS,KTE
!
!-----------------------------------------------------------------------
!       if(mype==0)then
!         write(nunit)(gsums(i,k),i=1,6)
!       endif
!!!     read(nunit)(gsums(i,k),i=1,6)
!-----------------------------------------------------------------------
!
        SUMPQ=GSUMS(1,K)
        SUMNQ=GSUMS(2,K)
!
!-----------------------------------------------------------------------
!***  FIRST MOMENT CONSERVING FACTOR
!-----------------------------------------------------------------------
!
        IF(SUMPQ>1.)THEN
          RFACQ=-SUMNQ/SUMPQ
        ELSE
          RFACQ=1.
        ENDIF
!
        IF(RFACQ<CONSERVE_MIN.OR.RFACQ>CONSERVE_MAX)RFACQ=1.
!
!-----------------------------------------------------------------------
!       if(mype==0.and.ntsd==181)close(nunit)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  IMPOSE CONSERVATION ON ANTI-FILTERING
!-----------------------------------------------------------------------
!
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IF(RFACQ<1.)THEN
            DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
              DQSTIJ=DQST(I,J,K)
              RFQIJ=HBM2(I,J)*(RFACQ-1.)+1.
              IF(DQSTIJ>=0.)DQSTIJ=DQSTIJ*RFQIJ
              Q(I,J,K)=Q(I,J,K)+DQSTIJ
            ENDDO
          ELSE
            DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
              DQSTIJ=DQST(I,J,K)
              RFQIJ=HBM2(I,J)*(RFACQ-1.)+1.
              IF(DQSTIJ<0.)DQSTIJ=DQSTIJ/RFQIJ
              Q(I,J,K)=Q(I,J,K)+DQSTIJ
            ENDDO
          ENDIF
        ENDDO
!
!-----------------------------------------------------------------------
!
        DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
        DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
          SCAL(I,J,K,L)=Q(I,J,K)
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDDO vertical_3
!
!-----------------------------------------------------------------------
!
      ENDDO scalar_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE HAD2_SCAL
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END MODULE MODULE_ADVECTION
!
!-----------------------------------------------------------------------
