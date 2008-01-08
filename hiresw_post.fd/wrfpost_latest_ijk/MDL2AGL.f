      SUBROUTINE MDL2AGL
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    MDL2P       VERT INTRP OF MODEL LVLS TO AGL HEIGHT
C   PRGRMMR: CHUANG           ORG: W/NP22     DATE: 05-05-23       
C     
C ABSTRACT:
C     FOR MOST APPLICATIONS THIS ROUTINE IS THE WORKHORSE
C     OF THE POST PROCESSOR.  IN A NUTSHELL IT INTERPOLATES
C     DATA FROM MODEL TO AGL HEIGHT SURFACES. 
C   .     
C     
C PROGRAM HISTORY LOG:
C   05-09-20  H CHUANG AND B ZHOU - ADD WIND DIFFERENCES OVER 2000 FT
C     
C USAGE:    CALL MDL2P
C   INPUT ARGUMENT LIST:
C
C   OUTPUT ARGUMENT LIST: 
C     NONE       
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       SCLFLD   - SCALE ARRAY ELEMENTS BY CONSTANT.
C       CALPOT   - COMPUTE POTENTIAL TEMPERATURE.
C       CALRH    - COMPUTE RELATIVE HUMIDITY.
C       CALDWP   - COMPUTE DEWPOINT TEMPERATURE.
C       BOUND    - BOUND ARRAY ELEMENTS BETWEEN LOWER AND UPPER LIMITS.
C       CALMCVG  - COMPUTE MOISTURE CONVERGENCE.
C       CALVOR   - COMPUTE ABSOLUTE VORTICITY.
C       CALSTRM  - COMPUTE GEOSTROPHIC STREAMFUNCTION.
C
C     LIBRARY:
C       COMMON   - CTLBLK
C                  RQSTFLD
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN 90
C     MACHINE : IBM SP
C$$$  
C
C
      use vrbls3d
      use vrbls2d
      use masks
C     
C     INCLUDE MODEL DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
C     GAMMA AND RGAMOG ARE USED IN THE EXTRAPOLATION OF VIRTUAL
C     TEMPERATURES BEYOND THE UPPER OF LOWER LIMITS OF DATA.
C     
!      INCLUDE "parmeta"
!      INCLUDE "parmout"
      INCLUDE "params"
C
!      PARAMETER (IM_JM=IM*JM,LMP1=LM+1)
      PARAMETER (GAMMA=6.5E-3,RGAMOG=RD*GAMMA/G)
      PARAMETER (LAGL=2,LAGL2=1)
C
C     INCLUDE COMMON BLOCKS.
      INCLUDE "CTLBLK.comm"
      INCLUDE "RQSTFLD.comm"
C     
C     DECLARE VARIABLES.
C     
      LOGICAL RUN,FIRST,RESTRT,SIGMA
      LOGICAL IOOMG,IOALL
!      REAL FSL(IM,JM),TSL(IM,JM),QSL(IM,JM)
!      REAL OSL(IM,JM),USL(IM,JM),VSL(IM,JM)
!      REAL Q2SL(IM,JM)
      REAL UAGL(IM,JM),VAGL(IM,JM)
      REAL GRID1(IM,JM),GRID2(IM,JM)
C
      INTEGER NL1X(IM,JM), IHE(JM),IHW(JM)
!      REAL TPRS(IM,JSTA_2L:JEND_2U,LSM)
!     + ,QPRS(IM,JSTA_2L:JEND_2U,LSM),FPRS(IM,JSTA_2L:JEND_2U,LSM)
C
!
!--- Definition of the following 2D (horizontal) dummy variables
!
!  C1D   - total condensate
!  QW1   - cloud water mixing ratio
!  QI1   - cloud ice mixing ratio
!  QR1   - rain mixing ratio
!  QS1   - snow mixing ratio
!  DBZ1  - radar reflectivity
!  DBZR1 - radar reflectivity from rain
!  DBZI1 - radar reflectivity from ice (snow + graupel + sleet)
!  DBZC1 - radar reflectivity from parameterized convection (bogused)
!
!      REAL C1D(IM,JM),QW1(IM,JM),QI1(IM,JM),QR1(IM,JM)
!     &,    QS1(IM,JM) ,DBZ1(IM,JM)
      REAL DBZ1(IM,JM),DBZR1(IM,JM),DBZI1(IM,JM),DBZC1(IM,JM)
     &,    ZAGL(LAGL),ZAGL2(LAGL2)
C
C      COMMON/JIMA/NL1X(IM,JM),ALPETUX(IM,JM),ALPET2X(IM,JM)
C
C     
C******************************************************************************
C
C     START MDL2P. 
C     
C     SET TOTAL NUMBER OF POINTS ON OUTPUT GRID.
C
C---------------------------------------------------------------
      ZAGL(1)=4000.
      ZAGL(2)=1000.
      ZAGL2(1)=609.6  ! 2000 ft
C
C     *** PART I ***
C
C     VERTICAL INTERPOLATION OF EVERYTHING ELSE.  EXECUTE ONLY
C     IF THERE'S SOMETHING WE WANT.
C
      IF (IGET(253).GT.0 .OR. IGET(279).GT.0 .OR. IGET(280).GT.0 .OR.
     &    IGET(281).GT.0 ) THEN
C
C---------------------------------------------------------------------
C***
C***  BECAUSE SIGMA LAYERS DO NOT GO UNDERGROUND,  DO ALL
C***  INTERPOLATION ABOVE GROUND NOW.
C***
C
        DO 310 LP=1,LAGL
         IF(LVLS(LP,IGET(253)).GT.0 .OR.
     &      LVLS(LP,IGET(279)).GT.0 .OR.
     &      LVLS(LP,IGET(280)).GT.0 .OR.
     &      LVLS(LP,IGET(281)).GT.0 ) THEN 
C
          jj=float(jsta+jend)/2.0
          ii=float(im)/3.0
          DO J=JSTA_2L,JEND_2U
          DO I=1,IM

C
	   DBZ1(I,J)=SPVAL
	   DBZR1(I,J)=SPVAL
	   DBZI1(I,J)=SPVAL
	   DBZC1(I,J)=SPVAL
C
C***  LOCATE VERTICAL INDEX OF MODEL MIDLAYER JUST BELOW
C***  THE AGL LEVEL TO WHICH WE ARE INTERPOLATING.
C
           LLMH=NINT(LMH(I,J))
           NL1X(I,J)=LLMH+1
           DO L=LLMH,2,-1
            ZDUM=ZMID(I,J,L)-ZINT(I,J,LLMH+1)
            IF(ZDUM.GE.ZAGL(LP))THEN
             NL1X(I,J)=L+1
	     GO TO 30
            ENDIF
           ENDDO
   30      CONTINUE	   
C
C  IF THE AGL LEVEL IS BELOW THE LOWEST MODEL MIDLAYER
C  BUT STILL ABOVE THE LOWEST MODEL BOTTOM INTERFACE,
C  WE WILL NOT CONSIDER IT UNDERGROUND AND THE INTERPOLATION
C  WILL EXTRAPOLATE TO THAT POINT
C
           IF(NL1X(I,J).EQ.(LLMH+1) .AND. ZAGL(LP).GT.0.)THEN
            NL1X(I,J)=LM
           ENDIF
C
c        if(NL1X(I,J).EQ.LMP1)print*,'Debug: NL1X=LMP1 AT '
c     1 ,i,j,lp
         ENDDO
         ENDDO
C
!mptest        IF(NHOLD.EQ.0)GO TO 310
C
!$omp  parallel do
!$omp& private(nn,i,j,ll,fact,qsat,rhl)
chc        DO 220 NN=1,NHOLD
chc        I=IHOLD(NN)
chc        J=JHOLD(NN)
c        DO 220 J=JSTA,JEND
         DO 220 J=JSTA_2L,JEND_2U
         DO 220 I=1,IM
          LL=NL1X(I,J)
C---------------------------------------------------------------------
C***  VERTICAL INTERPOLATION OF GEOPOTENTIAL, TEMPERATURE, SPECIFIC
C***  HUMIDITY, CLOUD WATER/ICE, OMEGA, WINDS, AND TKE.
C---------------------------------------------------------------------
C
CHC        IF(NL1X(I,J).LE.LM)THEN
          LLMH = NINT(LMH(I,J))
          IF(NL1X(I,J).LE.LLMH)THEN
C
C---------------------------------------------------------------------
C          INTERPOLATE LINEARLY IN LOG(P)
C***  EXTRAPOLATE ABOVE THE TOPMOST MIDLAYER OF THE MODEL
C***  INTERPOLATION BETWEEN NORMAL LOWER AND UPPER BOUNDS
C***  EXTRAPOLATE BELOW LOWEST MODEL MIDLAYER (BUT STILL ABOVE GROUND)
C---------------------------------------------------------------------
C
!          FACT=(ALSL(LP)-ALOG(PMID(I,J,LL)))/
!     &         (ALOG(PMID(I,J,LL))-ALOG(PMID(I,J,LL-1)))
           ZDUM=ZAGL(LP)+ZINT(I,J,NINT(LMH(I,J))+1)
           FACT=(ZDUM-ZMID(I,J,LL))
     &        /(ZMID(I,J,LL)-ZMID(I,J,LL-1))
C	  
	 DBZ1(I,J)=DBZ(I,J,LL)+(DBZ(I,J,LL)-DBZ(I,J,LL-1))*FACT
	 DBZR1(I,J)=DBZR(I,J,LL)+(DBZR(I,J,LL)-DBZR(I,J,LL-1))*FACT
	 DBZI1(I,J)=DBZI(I,J,LL)+(DBZI(I,J,LL)-DBZI(I,J,LL-1))*FACT
	 DBZC1(I,J)=DBZC(I,J,LL)+(DBZC(I,J,LL)-DBZC(I,J,LL-1))*FACT
c           IF(I.eq.ii.and.j.eq.jj)print*,'Debug AGL RADAR REF',
c     &     i,j,ll,zagl(lp),ZINT(I,J,NINT(LMH(I,J))+1)
c     &      ,ZMID(I,J,LL-1),ZMID(I,J,LL)
c     &     ,DBZ(I,J,LL-1),DBZ(I,J,LL),DBZ1(I,J)
c     &     ,DBZR(I,J,LL-1),DBZR(I,J,LL),DBZR1(I,J)
c     &     ,DBZI(I,J,LL-1),DBZI(I,J,LL),DBZI1(I,J)
c     &     ,DBZC(I,J,LL-1),DBZC(I,J,LL),DBZC1(I,J)
	   DBZ1(I,J)=AMAX1(DBZ1(I,J),DBZmin)
	   DBZR1(I,J)=AMAX1(DBZR1(I,J),DBZmin)
	   DBZI1(I,J)=AMAX1(DBZI1(I,J),DBZmin)
	   DBZC1(I,J)=AMAX1(DBZC1(I,J),DBZmin)
C
C FOR UNDERGROUND AGL LEVELS, ASSUME TEMPERATURE TO CHANGE 
C ADIABATICLY, RH TO BE THE SAME AS THE AVERAGE OF THE 2ND AND 3RD
C LAYERS FROM THE GOUND, WIND TO BE THE SAME AS THE LOWEST LEVEL ABOVE
C GOUND
          ELSE
	   DBZ1(I,J)=DBZmin
	   DBZR1(I,J)=DBZmin
	   DBZI1(I,J)=DBZmin
	   DBZC1(I,J)=DBZmin
          END IF
  220    CONTINUE
C
C     
C---------------------------------------------------------------------
C        *** PART II ***
C---------------------------------------------------------------------
C
C        OUTPUT SELECTED FIELDS.
C
C---------------------------------------------------------------------
C
C
!---  Radar Reflectivity
          IF((IGET(253).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=DBZ1(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
             ID(11) = NINT(ZAGL(LP))
             CALL GRIBIT(IGET(253),LP,GRID1,IM,JM)
          END IF    
!---  Radar reflectivity from rain
          IF((IGET(279).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=DBZR1(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
             ID(11) = NINT(ZAGL(LP))
             CALL GRIBIT(IGET(279),LP,GRID1,IM,JM)
          END IF    
!---  Radar reflectivity from all ice habits (snow + graupel + sleet, etc.)
          IF((IGET(280).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=DBZI1(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
             ID(11) = NINT(ZAGL(LP))
             CALL GRIBIT(IGET(280),LP,GRID1,IM,JM)
          END IF    
!---  Radar reflectivity from parameterized convection
          IF((IGET(281).GT.0) )THEN
             DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=DBZC1(I,J)
             ENDDO
             ENDDO
             ID(1:25)=0
             ID(02)=129
             ID(11) = NINT(ZAGL(LP))
             CALL GRIBIT(IGET(281),LP,GRID1,IM,JM)
          END IF    
C          
         ENDIF ! FOR LEVEL
C     
C***  END OF MAIN VERTICAL LOOP
C     
  310   CONTINUE
C***  ENDIF FOR IF TEST SEEING IF WE WANT ANY OTHER VARIABLES
C
      ENDIF
C
      IF((IGET(259).GT.0) )THEN
C
C---------------------------------------------------------------------
C***
C***  BECAUSE SIGMA LAYERS DO NOT GO UNDERGROUND,  DO ALL
C***  INTERPOLATION ABOVE GROUND NOW.
C***
C
        DO 320 LP=1,LAGL2
         IF(LVLS(LP,IGET(259)).GT.0)THEN 
C
          jj=float(jsta+jend)/2.0
          ii=float(im)/3.0
          DO J=JSTA_2L,JEND_2U
          DO I=1,IM

C
	   UAGL(I,J)=SPVAL
	   VAGL(I,J)=SPVAL
C
C***  LOCATE VERTICAL INDEX OF MODEL MIDLAYER JUST BELOW
C***  THE AGL LEVEL TO WHICH WE ARE INTERPOLATING.
C
           LLMH=NINT(LMH(I,J))
           NL1X(I,J)=LLMH+1
           DO L=LLMH,2,-1
            ZDUM=ZMID(I,J,L)-ZINT(I,J,LLMH+1)
            IF(ZDUM.GE.ZAGL2(LP))THEN
             NL1X(I,J)=L+1
	     GO TO 40
            ENDIF
           ENDDO
   40      CONTINUE	   
C
C  IF THE AGL LEVEL IS BELOW THE LOWEST MODEL MIDLAYER
C  BUT STILL ABOVE THE LOWEST MODEL BOTTOM INTERFACE,
C  WE WILL NOT CONSIDER IT UNDERGROUND AND THE INTERPOLATION
C  WILL EXTRAPOLATE TO THAT POINT
C
           IF(NL1X(I,J).EQ.(LLMH+1) .AND. ZAGL2(LP).GT.0.)THEN
            NL1X(I,J)=LM
           ENDIF
C
c        if(NL1X(I,J).EQ.LMP1)print*,'Debug: NL1X=LMP1 AT '
c     1 ,i,j,lp
         ENDDO
         ENDDO
C
!mptest        IF(NHOLD.EQ.0)GO TO 310
C
!$omp  parallel do
!$omp& private(nn,i,j,ll,fact,qsat,rhl)
chc        DO 220 NN=1,NHOLD
chc        I=IHOLD(NN)
chc        J=JHOLD(NN)
c        DO 220 J=JSTA,JEND
         DO J=JSTA_2L,JEND_2U
          IF(MODELNAME .EQ. 'NCAR' .OR.MODELNAME.EQ.'RSM')THEN
           IHW(J)=-1
           IHE(J)=1 
          ELSE
           IHW(J)=-MOD(J,2)
           IHE(J)=IHW(J)+1
          END IF	
         ENDDO 
	 
         DO 230 J=JSTA_M,JEND_M
         DO 230 I=2,IM-1
          LL=NL1X(I,J)
C---------------------------------------------------------------------
C***  VERTICAL INTERPOLATION OF GEOPOTENTIAL, TEMPERATURE, SPECIFIC
C***  HUMIDITY, CLOUD WATER/ICE, OMEGA, WINDS, AND TKE.
C---------------------------------------------------------------------
C
CHC        IF(NL1X(I,J).LE.LM)THEN
          LLMH = NINT(LMH(I,J))
          IF(NL1X(I,J).LE.LLMH)THEN
C
C---------------------------------------------------------------------
C          INTERPOLATE LINEARLY IN LOG(P)
C***  EXTRAPOLATE ABOVE THE TOPMOST MIDLAYER OF THE MODEL
C***  INTERPOLATION BETWEEN NORMAL LOWER AND UPPER BOUNDS
C***  EXTRAPOLATE BELOW LOWEST MODEL MIDLAYER (BUT STILL ABOVE GROUND)
C---------------------------------------------------------------------
C
!          FACT=(ALSL(LP)-ALOG(PMID(I,J,LL)))/
!     &         (ALOG(PMID(I,J,LL))-ALOG(PMID(I,J,LL-1)))
           ZDUM=ZAGL2(LP)+ZINT(I,J,NINT(LMH(I,J))+1)
           FACT=(ZDUM-ZMID(I,J,LL))
     &        /(ZMID(I,J,LL)-ZMID(I,J,LL-1))
C	  
           IF(MODELNAME .EQ. 'NCAR' .OR.MODELNAME.EQ.'RSM')THEN
	    UAGLU=UH(I,J,LL-1)
            UAGLL=UH(I,J,LL)
     
            VAGLU=VH(I,J,LL-1)
            VAGLL=VH(I,J,LL)
	   ELSE
            UAGLU=(UH(I+IHE(J),J,LL-1)
     &	      +UH(I+IHW(J),J,LL-1)+
     &	       UH(I,J-1,LL-1)+UH(I,J+1,LL-1))/4.0
            UAGLL=(UH(I+IHE(J),J,LL)
     &	      +UH(I+IHW(J),J,LL)+
     &	       UH(I,J-1,LL)+UH(I,J+1,LL))/4.0
     
            VAGLU=(VH(I+IHE(J),J,LL-1)
     &	      +VH(I+IHW(J),J,LL-1)+
     &	       VH(I,J-1,LL-1)+VH(I,J+1,LL-1))/4.0
            VAGLL=(VH(I+IHE(J),J,LL)
     &	      +VH(I+IHW(J),J,LL)+
     &	       VH(I,J-1,LL)+VH(I,J+1,LL))/4.0
           END IF
           UAGL(I,J)=UAGLL+(UAGLL-UAGLU)*FACT
	   VAGL(I,J)=VAGLL+(VAGLL-VAGLU)*FACT
     
C
C FOR UNDERGROUND AGL LEVELS, ASSUME TEMPERATURE TO CHANGE 
C ADIABATICLY, RH TO BE THE SAME AS THE AVERAGE OF THE 2ND AND 3RD
C LAYERS FROM THE GOUND, WIND TO BE THE SAME AS THE LOWEST LEVEL ABOVE
C GOUND
          ELSE
	   IF(MODELNAME .EQ. 'NCAR' .OR.MODELNAME.EQ.'RSM')THEN
            UAGL(I,J)=UH(I,J,NINT(LMV(I,J)))  
	    VAGL(I,J)=VH(I,J,NINT(LMV(I,J)))
	   ELSE
	    UAGL(I,J)=(UH(I+IHE(J),J,NINT(LMV(I+IHE(J),J)))
     &	      +UH(I+IHW(J),J,NINT(LMV(I+IHW(J),J)))+
     &	       UH(I,J-1,NINT(LMV(I,J-1)))+UH(I,J+1,NINT(LMV(I,J+1))))/4.0   
	    VAGL(I,J)=(VH(I+IHE(J),J,NINT(LMV(I+IHE(J),J)))
     &	      +VH(I+IHW(J),J,NINT(LMV(I+IHW(J),J)))+
     &	       VH(I,J-1,NINT(LMV(I,J-1)))+VH(I,J+1,NINT(LMV(I,J+1))))/4.0
           END IF
          END IF
  230    CONTINUE
C
C     
C---------------------------------------------------------------------
C        *** PART II ***
C---------------------------------------------------------------------
C
C        OUTPUT SELECTED FIELDS.
C
C---------------------------------------------------------------------
C
C
C---  Wind Shear (wind speed difference in knots between sfc and 2000 ft)

	     DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J)=SQRT((UAGL(I,J)-U10(I,J))**2+
     +	       (VAGL(I,J)-V10(I,J))**2)*1.943*ZAGL2(LP)/(ZAGL2(LP)-10.) 
             ENDDO
             ENDDO
             ID(1:25)=0
	     ID(10) = NINT(ZAGL2(LP))
             ID(11) = 0
             CALL GRIBIT(IGET(259),LP,GRID1,IM,JM)
C          
         ENDIF ! FOR LEVEL
C     
C***  END OF MAIN VERTICAL LOOP
C     
  320   CONTINUE
C***  ENDIF FOR IF TEST SEEING IF WE WANT ANY OTHER VARIABLES
C
      ENDIF

C
C     END OF ROUTINE.
C
      RETURN
      END
