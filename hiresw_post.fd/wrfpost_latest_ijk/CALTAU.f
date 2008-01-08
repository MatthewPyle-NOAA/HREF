      SUBROUTINE CALTAU(TAUX,TAUY)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    CALTAU      COMPUTE U AND V WIND STRESSES
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-09-01
C     
C ABSTRACT:  THIS ROUTINE COMPUTES SURFACE LAYER U AND V
C   WIND COMPONENT STRESSES USING K THEORY AS PRESENTED
C   IN SECTION 8.4 OF "NUMBERICAL PREDICTION AND DYNAMIC
C   METEOROLOGY" BY HALTINER AND WILLIAMS (1980, JOHN WILEY
C   & SONS).
C   .     
C     
C PROGRAM HISTORY LOG:
C   93-09-01  RUSS TREADON
C   98-06-11  T BLACK - CONVERSION FROM 1-D TO 2-D
C   00-01-04  JIM TUCCILLO - MPI VERSION
C   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID OUTPUT
C   02-01-15  MIKE BALDWIN - WRF VERSION, OUTPUT IS ON MASS-POINTS
C   05-02-23  H CHUANG - COMPUTE STRESS FOR NMM ON WIND POINTS
C   05-07-07  BINBIN ZHOU - ADD RSM STRESS for A GRID     
C USAGE:    CALL CALTAU(TAUX,TAUY)
C   INPUT ARGUMENT LIST:
C     NONE     
C
C   OUTPUT ARGUMENT LIST: 
C     TAUX     - SUFACE LAYER U COMPONENT WIND STRESS.
C     TAUY     - SUFACE LAYER V COMPONENT WIND STRESS.
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       CLMAX
C       MIXLEN
C
C     LIBRARY:
C       COMMON   - 
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY C-90
C$$$  
C     
C
      use vrbls3d
      use vrbls2d
      use masks
C
C     INCLUDE/SET PARAMETERS.
C     
!      INCLUDE "parmeta"
!      INCLUDE "parmout"
      INCLUDE "params"
C
C     INCLUDE COMMON BLOCKS.
C
      INCLUDE "CTLBLK.comm"
C     
C     DECLARE VARIABLES.
      INTEGER KK(4),IVE(JM),IVW(JM)
      REAL TAUX(IM,JM),TAUY(IM,JM)
      REAL EL0(IM,JM)
      REAL, ALLOCATABLE :: EL(:,:,:)
      REAL EGRIDU(IM,JM),EGRIDV(IM,JM),EGRID4(IM,JM),EGRID5(IM,JM)
      REAL UZ0H(IM,JM),VZ0H(IM,JM)
      CHARACTER*1 AGRID
C     
C     
C********************************************************************
C     START CALTAU HERE.
C    
      ALLOCATE (EL(IM,JSTA_2L:JEND_2U,LM))
C
C     COMPUTE MASTER LENGTH SCALE.
C
!      CALL CLMAX(EL0,EGRIDU,EGRIDV,EGRID4,EGRID5)
!      CALL MIXLEN(EL0,EL)
C     
C     INITIALIZE OUTPUT AND WORK ARRAY TO ZERO.
C     
      DO J=JSTA,JEND
      DO I=1,IM
        EGRIDU(I,J) = D00
        EGRIDV(I,J) = D00
        TAUX(I,J)   = SPVAL
        TAUY(I,J)   = SPVAL
      ENDDO
      ENDDO
C     
C     COMPUTE SURFACE LAYER U AND V WIND STRESSES.
C
C     ASSUME THAT U AND V HAVE UPDATED HALOS
C
      IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
       CALL CLMAX(EL0,EGRIDU,EGRIDV,EGRID4,EGRID5)
       CALL MIXLEN(EL0,EL)

       DO J=JSTA,JEND
       DO I=1,IM
C
        LMHK = NINT(LMH(I,J))
C
C       COMPUTE THICKNESS OF LAYER AT MASS POINT.
C
        DZ  = D50*(ZINT(I,J,LMHK)-ZINT(I,J,LMHK+1))
        DZ  = DZ-Z0(I,J)
        RDZ = 1./DZ
C
C        COMPUTE REPRESENTATIVE AIR DENSITY.
C
        PSFC = PMID(I,J,LMHK)
        TV   = (H1+D608*Q(I,J,LMHK))*T(I,J,LMHK)
        RHO  = PSFC/(RD*TV)
C     
C        COMPUTE A MEAN MASS POINT WIND IN THE 
C        FIRST ATMOSPHERIC ETA LAYER.
C
        ULMH = UH(I,J,LMHK)
        VLMH = VH(I,J,LMHK)
C
C       COMPUTE WIND SHEAR COMPONENTS ACROSS LAYER.
C
        DELUDZ = (ULMH-UZ0(I,J))*RDZ
        DELVDZ = (VLMH-VZ0(I,J))*RDZ
C     
C       COMPUTE U (EGRIDU) AND V (EGRIDV) WIND STRESSES.
C
        ELSQR     = EL(I,J,LMHK)*EL(I,J,LMHK)
        TAUX(I,J) = RHO*ELSQR*DELUDZ*DELUDZ
        TAUY(I,J) = RHO*ELSQR*DELVDZ*DELVDZ

C
       END DO
       END DO
      ELSE IF(MODELNAME .EQ. 'NMM')THEN

       DO J=JSTA_M,JEND_M
        IVE(J)=MOD(J,2)
        IVW(J)=IVE(J)-1
       ENDDO
 
       DO J=JSTA_M,JEND_M
       DO I=2,IM-1
C
        LMHK = NINT(LMH(I,J)) 
        IE=I+IVE(J)
        IW=I+IVW(J)
        ZINT1=(ZINT(IW,J,LMHK)+ZINT(IE,J,LMHK)
     &   +ZINT(I,J+1,LMHK)+ZINT(I,J-1,LMHK))*D25
        ZINT2=(ZINT(IW,J,LMHK+1)+ZINT(IE,J,LMHK+1)
     &   +ZINT(I,J+1,LMHK+1)+ZINT(I,J-1,LMHK+1))*D25
        DZ  = D50*(ZINT1-ZINT2)       
        Z0V=(Z0(IW,J)+Z0(IE,J)+Z0(I,J+1)+Z0(I,J-1))*D25
        DZ  = DZ-Z0V
        RDZ = 1./DZ
C
C        COMPUTE REPRESENTATIVE AIR DENSITY.
C
        PSFC = (PMID(IW,J,LMHK)+PMID(IE,J,LMHK)
     &   +PMID(I,J+1,LMHK)+PMID(I,J-1,LMHK))*D25 
        TVV = (T(IW,J,LMHK)+T(IE,J,LMHK)
     &   +T(I,J+1,LMHK)+T(I,J-1,LMHK))*D25
        QVV = (Q(IW,J,LMHK)+Q(IE,J,LMHK)
     &   +Q(I,J+1,LMHK)+Q(I,J-1,LMHK))*D25
        TV   = (H1+D608*QVV)*TVV
        RHO  = PSFC/(RD*TV) 

C       COMPUTE WIND SHEAR COMPONENTS ACROSS LAYER.
C
        DELUDZ = (U(I,J,LMHK)-UZ0(I,J))*RDZ
        DELVDZ = (V(I,J,LMHK)-VZ0(I,J))*RDZ 

C       COMPUTE U (EGRIDU) AND V (EGRIDV) WIND STRESSES.
C                                       
        ELV1=(EL_MYJ(IW,J,LMHK)+EL_MYJ(IE,J,LMHK)
     &   +EL_MYJ(I,J+1,LMHK)+EL_MYJ(I,J-1,LMHK))*D25
        ELV2=(EL_MYJ(IW,J,LMHK-1)+EL_MYJ(IE,J,LMHK-1)
     &   +EL_MYJ(I,J+1,LMHK-1)+EL_MYJ(I,J-1,LMHK-1))*D25
        ELV=(ELV1+ELV2)/2.0  ! EL is defined at the bottom of layer
        ELSQR       =ELV*ELV
        TAUX(I,J)=RHO*ELSQR*DELUDZ*DELUDZ 
        TAUY(I,J)=RHO*ELSQR*DELVDZ*DELVDZ
	ii=im/2
	jj=(jsta+jend)/2
!        if(i.eq.ii.and.j.eq.jj)print*,'sample tau'
!     &	,RHO,ELSQR,DELUDZ,DELVDZ  
       END DO
       END DO

      END IF
C     
      DEALLOCATE(EL)
C     END OF ROUTINE.
C     
      RETURN
      END
