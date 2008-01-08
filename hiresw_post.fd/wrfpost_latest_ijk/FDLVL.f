      SUBROUTINE FDLVL(ITYPE,TFD,UFD,VFD)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    FDLVL       COMPUTES FD LEVEL T, U, V
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-22       
C     
C ABSTRACT:
C     THIS ROUTINE COMPUTES TEMPERATURE, U WIND COMPONENT,
C     AND V WIND COMPONENT ON THE NFD=6 FD LEVELS.  THE
C     HEIGHT OF THESE LEVELS (IN METERS) IS GIVEN IN THE 
C     DATA STATEMENT BELOW.  THE ALGORITHM PROCEEDS AS 
C     FOLLOWS. (AGL IN PARENTHESES)
C     
C     AT EACH MASS POINT MOVE UP VERTICALLY FROM THE LM-TH (LOWEST
C     ATMOSPHERIC) ETA LAYER.  FIND THE ETA LAYERS WHOSE 
C     HEIGHT (ABOVE GROUND) BOUNDS THE TARGET FD LEVEL HEIGHT.
C     VERTICALLY INTERPOLATE TO GET TEMPERATURE AT THIS FD
C     LEVEL.  AVERAGE THE FOUR SURROUNDING WINDS
C     TO GET A MASS POINT WIND.  VERTICALLY INTERPOLATE THESE
C     MASS POINT WINDS TO THE TARGET FD LEVEL.  CONTINUE THIS
C     PROCESS UNTIL ALL NFD=6 FD LEVELS HAVE BEEN PROCESSED.
C     MOVE ON TO THE NEXT MASS POINT.  
C     
C     AVERAGING THE FOUR ABOVE GROUND WINDS TO THE MASS POINT
C     WAS FOUND TO SMOOTH THE FIELD AND REDUCE THE OCCURRENCE
C     OF POINT PEAK WINDS FAR IN EXCESS OF THE WINDS AT 
C     ADJACENT POINTS.  MASS POINT VALUES ARE RETURNED.
C   .     
C     
C PROGRAM HISTORY LOG:
C   92-12-22  RUSS TREADON
C   93-11-23  RUSS TREADON - CORRECTED ROUTINE TO COMPUTE
C             FD LEVELS WITH REPECT TO MEAN SEA LEVEL.
C   94-01-04  MICHAEL BALDWIN - INCLUDE OPTIONS FOR COMPUTING
C                               EITHER AGL OR MSL
C   98-06-15  T BLACK - CONVERSION FROM 1-D TO 2-D
C   00-01-04  JIM TUCCILLO - MPI VERSION            
C   02-01-15  MIKE BALDWIN - WRF VERSION
C     
C USAGE:    CALL FDLVL(ITYPE,TFD,UFD,VFD)
C   INPUT ARGUMENT LIST:
C     ITYPE    - FLAG THAT DETERMINES WHETHER MSL (1) OR AGL (2)
C                   LEVELS ARE USED.
C
C   OUTPUT ARGUMENT LIST: 
C     TFD      - TEMPERATURE (K) ON FD LEVELS.
C     UFD      - U WIND (M/S) ON FD LEVELS.
C     VFD      - V WIND (M/S) ON FD LEVELS.
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C
C     LIBRARY:
C       COMMON   - 
C                  LOOPS
C                  MASKS
C                  OPTIONS
C                  INDX
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
C     
C     SET NUMBER OF FD LEVELS.
      PARAMETER (NFD=11)  ! adding 6000 m
C     
C     INCLUDE PARAMETERS.
!      INCLUDE "parmeta"
      INCLUDE "params"
C
C     INCLUDE COMMON BLOCKS.
      INCLUDE "CTLBLK.comm"
C     
C     DECLARE VARIABLES
C     
      INTEGER LVL(NFD),LHL(NFD)
      INTEGER IVE(JM),IVW(JM),JVN,JVS
      REAL DZABV(NFD), HTFD(NFD),DZABH(NFD)
      REAL TFD(IM,JM,NFD),UFD(IM,JM,NFD)
      REAL VFD(IM,JM,NFD)
      LOGICAL DONEH, DONEV
C
C     SET FD LEVEL HEIGHTS IN METERS.
      DATA HTFD  / 305.E0,457.E0,610.E0,914.E0,1524.E0,1829.E0,
     X     2134.E0,2743.E0,3658.E0,4572.E0,6000.E0/
C     
C****************************************************************
C     START FDLVL HERE
C     
C     INITIALIZE ARRAYS.
C     
!$omp  parallel do
      DO 10 IFD = 1,NFD
      DO J=JSTA,JEND
      DO I=1,IM
         TFD(I,J,IFD)    = SPVAL
         UFD(I,J,IFD)    = SPVAL
         VFD(I,J,IFD)    = SPVAL
      ENDDO
      ENDDO
 10   CONTINUE

      IF(MODELNAME .EQ. 'NMM')THEN
        JVN=1
        JVS=-1
        do J=JSTA,JEND
         IVE(J)=MOD(J,2)
         IVW(J)=IVE(J)-1
        enddo
      END IF

C
C     MSL FD LEVELS
C
      IF (ITYPE.EQ.1) THEN
	write(6,*) 'computing above MSL'
C     
C     LOOP OVER HORIZONTAL GRID.
C
      DO 50 J=JSTA_M,JEND_M
      DO 50 I=2,IM-1
         HTSFC = FIS(I,J)*GI
         LLMH  = NINT(LMH(I,J))
         IFD = 1
C     
C        LOCATE VERTICAL INDICES OF T,U,V, LEVEL JUST 
C        ABOVE EACH FD LEVEL.
C
        DO 22 IFD = 1, NFD
	 DONEH=.FALSE.
	 DONEV=.FALSE.
          DO 20 L = LM,1,-1
	    HTT = ZMID(I,J,L)
	    IF(MODELNAME .EQ. 'NMM')THEN
	     IE=I+IVE(J)
             IW=I+IVW(J)
             JN=J+JVN
             JS=J+JVS
             HTTUV = 0.25*(ZMID(IW,J,L)
     1          +ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))
            ELSE
	     HTTUV = HTT
	    END IF
	    
	    IF (.NOT. DONEH .AND. HTT.GT.HTFD(IFD)) THEN
               LHL(IFD)   = L
               DZABH(IFD) = HTT-HTFD(IFD)
	       DONEH=.TRUE.
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
               IF(HTSFC.GT.HTFD(IFD)) THEN
!mp
                LHL(IFD)=LM+1  ! CHUANG: changed to lm+1
!mp
	       ENDIF
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
!               IFD        = IFD + 1
!               IF (IFD.GT.NFD) GOTO 30
	    END IF   
	     
            IF (.NOT. DONEV .AND. HTTUV.GT.HTFD(IFD)) THEN
               LVL(IFD)   = L
               DZABV(IFD) = HTTUV-HTFD(IFD)
	       DONEV=.TRUE.
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
               IF(HTSFC.GT.HTFD(IFD)) THEN
!mp
                LVL(IFD)=LM+1  ! CHUANG: changed to lm+1
!mp
	       ENDIF
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
!               IFD        = IFD + 1
!               IF (IFD.GT.NFD) GOTO 30
            ENDIF
	    
	    IF(DONEH .AND. DONEV)GO TO 22 
 20       CONTINUE
 22     CONTINUE   	
C     
C        COMPUTE T, U, AND V AT FD LEVELS.
C
 30      CONTINUE
C
         DO 40 IFD = 1,NFD
	 
	    L = LHL(IFD)
            IF (L.LT.LM) THEN
               DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
               RDZ  = 1./DZ
               DELT = T(I,J,L)-T(I,J,L+1)
               TFD(I,J,IFD) = T(I,J,L) - DELT*RDZ*DZABH(IFD)
            ELSEIF (L.EQ.LM) THEN
               TFD(I,J,IFD) = T(I,J,L)
            ENDIF
	    
            L = LVL(IFD)
            IF (L.LT.LM) THEN
	      IF(MODELNAME .EQ. 'NMM')THEN
	       IE=I+IVE(J)
               IW=I+IVW(J)
               JN=J+JVN
               JS=J+JVS
               Z1 = 0.25*(ZMID(IW,J,L)
     1           +ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))
               Z2 = 0.25*(ZMID(IW,J,L+1)
     1           +ZMID(IE,J,L+1)+ZMID(I,JN,L+1)+ZMID(I,JS,L+1))
               DZ = Z1-Z2
              ELSE
               DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
	      END IF 
               RDZ  = 1./DZ
               DELU = UH(I,J,L) - UH(I,J,L+1)
               DELV = VH(I,J,L) - VH(I,J,L+1)
               UFD(I,J,IFD) = UH(I,J,L) - DELU*RDZ*DZABV(IFD)
               VFD(I,J,IFD) = VH(I,J,L) - DELV*RDZ*DZABV(IFD)
            ELSEIF (L.EQ.LM) THEN
               UFD(I,J,IFD)=UH(I,J,L)
               VFD(I,J,IFD)=VH(I,J,L)
            ENDIF
 40      CONTINUE
C     
C     COMPUTE FD LEVEL T, U, AND V AT NEXT K.
C
 50   CONTINUE
C     END OF MSL FD LEVELS
      ELSE
	write(6,*) 'computing above AGL'
C
C     AGL FD LEVELS 
C
C     
C     LOOP OVER HORIZONTAL GRID.
C     
      DO 250 J=JSTA_M,JEND_M
      DO 250 I=2,IM-1
         HTSFC = FIS(I,J)*GI
         LLMH  = NINT(LMH(I,J))
!         IFD   = 1
C     
C        LOCATE VERTICAL INDICES OF T,U,V, LEVEL JUST 
C        ABOVE EACH FD LEVEL.
C
        DO 222 IFD = 1, NFD
	 DONEH=.FALSE.
         DONEV=.FALSE.
         DO 220 L = LLMH,1,-1
	    HTABH = ZMID(I,J,L)-HTSFC
	    IF(MODELNAME .EQ. 'NMM')THEN
             IE=I+IVE(J)
             IW=I+IVW(J)
             JN=J+JVN
             JS=J+JVS
             HTABV = 0.25*(ZMID(IW,J,L)
     1          +ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))-HTSFC  
            ELSE
	     HTABV = HTABH
	    END IF
	    
	    IF (.NOT. DONEH .AND. HTABH.GT.HTFD(IFD)) THEN
               LHL(IFD)   = L
               DZABH(IFD) = HTABH-HTFD(IFD)
	       DONEH=.TRUE.
!               IFD        = IFD + 1
!               IF (IFD.GT.NFD) GOTO 230
            ENDIF
	     
            IF (.NOT. DONEV .AND. HTABV.GT.HTFD(IFD)) THEN
               LVL(IFD)   = L
               DZABV(IFD) = HTABV-HTFD(IFD)
	       DONEV=.TRUE.
!               IFD        = IFD + 1
!               IF (IFD.GT.NFD) GOTO 230
            ENDIF
            IF(DONEH .AND. DONEV)GO TO 222	    
 220        CONTINUE
C     
C        COMPUTE T, U, AND V AT FD LEVELS.
C
 222     CONTINUE
C
         DO 240 IFD = 1,NFD
	    L = LHL(IFD)
            IF (L.LT.LM) THEN
               DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
               RDZ  = 1./DZ
               DELT = T(I,J,L)-T(I,J,L+1)
               TFD(I,J,IFD) = T(I,J,L) - DELT*RDZ*DZABH(IFD)
            ELSE
               TFD(I,J,IFD) = T(I,J,L)
            ENDIF
	    
            L = LVL(IFD)
            IF (L.LT.LM) THEN
	      IF(MODELNAME .EQ. 'NMM')THEN
               IE=I+IVE(J)
               IW=I+IVW(J)
               JN=J+JVN
               JS=J+JVS
               Z1 = 0.25*(ZMID(IW,J,L)
     1           +ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))
               Z2 = 0.25*(ZMID(IW,J,L+1)
     1           +ZMID(IE,J,L+1)+ZMID(I,JN,L+1)+ZMID(I,JS,L+1))
               DZ = Z1-Z2
              ELSE
               DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
              END IF
              RDZ  = 1./DZ
              DELU = UH(I,J,L)-UH(I,J,L+1)
              DELV = VH(I,J,L)-VH(I,J,L+1)
              UFD(I,J,IFD) = UH(I,J,L) - DELU*RDZ*DZABV(IFD)
              VFD(I,J,IFD) = VH(I,J,L) - DELV*RDZ*DZABV(IFD)
            ELSE
              UFD(I,J,IFD) = UH(I,J,L)
              VFD(I,J,IFD) = VH(I,J,L)
            ENDIF
 240     CONTINUE
C     
C     COMPUTE FD LEVEL T, U, AND V AT NEXT K.
C
 250  CONTINUE
C     END OF AGL FD LEVELS
      ENDIF
C
C     END OF ROUTINE.
C
      RETURN
      END