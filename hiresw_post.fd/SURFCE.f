      SUBROUTINE SURFCE
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    SURFCE      POST SURFACE BASED FIELDS
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-21       
C     
C ABSTRACT:
C     THIS ROUTINE POSTS SURFACE BASED FIELDS.
C   .     
C     
C PROGRAM HISTORY LOG:
C   92-12-21  RUSS TREADON
C   94-08-04  MICHAEL BALDWIN - ADDED OUTPUT OF SFC FLUXES OF
C                               SENS AND LATENT HEAT AND THETA AT Z0
C   94-11-04  MICHAEL BALDWIN - ADDED INSTANTANEOUS PRECIP TYPE
C   96-03-19  MICHAEL BALDWIN - CHANGE SOIL PARAMETERS
C   96-09-25  MICHAEL BALDWIN - ADDED SNOW RATIO FROM EXPLICIT SCHEME
C   96-10-17  MICHAEL BALDWIN - CHANGED SFCEVP,POTEVP TO ACCUM.  TOOK
C                               OUT -PTRACE FOR ACSNOW,SSROFF,BGROFF.
C   97-04-23  MICHAEL BALDWIN - TOOK OUT -PTRACE FOR ALL PRECIP FIELDS
C   98-06-12  T BLACK         - CONVERSION FROM 1-D TO 2-D
C   98-07-17  MIKE BALDWIN - REMOVED LABL84
C   98-08-18  MIKE BALDWIN - COMPUTE RH OVER ICE
C   98-12-22  MIKE BALDWIN - BACK OUT RH OVER ICE
C   00-01-04  JIM TUCCILLO - MPI VERSION
C   01-10-22  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
C   02-06-11  MIKE BALDWIN - WRF VERSION ASSUMING ALL ACCUM VARS
C                            HAVE BUCKETS THAT FILL FROM T=00H ON
C   02-08-28  H CHUANG - COMPUTE FIELDS AT SHELTER LEVELS FOR WRF
C   04-12-09  H CHUANG - ADD ADDITIONAL LSM FIELDS
C   05-07-07  BINBIN ZHOU - ADD RSM MODEL
C   05-08-24  GEOFF MANIKIN - ADDED DOMINANT PRECIP TYPE
C     
C USAGE:    CALL SURFCE
C   INPUT ARGUMENT LIST:
C
C   OUTPUT ARGUMENT LIST: 
C     
C   OUTPUT FILES:
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       BOUND    - ENFORCE LOWER AND UPPER LIMITS ON ARRAY ELEMENTS.
C       DEWPOINT - COMPUTE DEWPOINT TEMPERATURE.
C       CALDRG   - COMPUTE SURFACE LAYER DRAG COEFFICENT
C       CALTAU   - COMPUTE SURFACE LAYER U AND V WIND STRESSES.
C
C     LIBRARY:
C       COMMON   - CTLBLK
C                  RQSTFLD
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : CRAY C-90
C$$$  
C
C     
C     INCLUDE GRID DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
C
      use vrbls3d   
      use vrbls2d   
      use soil
      use masks
C     
!      INCLUDE "parmeta"
!      INCLUDE "parmout"
      INCLUDE "params"
!      INCLUDE "parmsoil"
C
C     INCLUDE COMMON BLOCKS.
!      COMMON /LEAF/ XLAI
      INCLUDE "CTLBLK.comm"
      INCLUDE "RQSTFLD.comm"
!      INCLUDE "SOILDEPTH.comm"
C     
C     IN NGM SUBROUTINE OUTPUT WE FIND THE FOLLOWING COMMENT.
C     "IF THE FOLLOWING THRESHOLD VALUES ARE CHANGED, CONTACT
C     TDL/SYNOPTIC-SCALE TECHNIQUES BRANCH (PAUL DALLAVALLE
C     AND JOHN JENSENIUS).  THEY MAY BE USING IT IN ONE OF 
C     THEIR PACKING CODES."  THE THRESHOLD VALUE IS 0.01 INCH
C     OR 2.54E-4 METER.  PRECIPITATION VALUES LESS THAN THIS
C     THRESHOLD ARE SET TO MINUS ONE TIMES THIS THRESHOLD.
      PARAMETER (PTRACE = 0.000254E0)
C     
C     SET CELCIUS TO KELVIN AND SECOND TO HOUR CONVERSION.
      PARAMETER (C2K    = 273.15)
      PARAMETER (NALG    = 5)
      PARAMETER (SEC2HR = 1./3600.)
      PARAMETER (PTHRESH=0.000004)
C     
C     DECLARE VARIABLES.
C     
      LOGICAL RUN,FIRST,RESTRT,SIGMA
      INTEGER IWX1(IM,JM),NROOTS(IM,JM),IWX4(IM,JM),IWX5(IM,JM)
      REAL IWX2(IM,JM),IWX3(IM,JM)
      REAL PSFC(IM,JM),TSFC(IM,JM),QSFC(IM,JM),RHSFC(IM,JM)
      REAL ZSFC(IM,JM),THSFC(IM,JM),DWPSFC(IM,JM),EVP(IM,JM)
!      REAL ANCPRC(IM,JM),P1D(IM,JM),T1D(IM,JM),Q1D(IM,JM)
      REAL P1D(IM,JM),T1D(IM,JM),Q1D(IM,JM),ZWET(IM,JM)
      REAL EGRID1(IM,JM),EGRID2(IM,JM),UA(IM,JM),VA(IM,JM)
      REAL GRID1(IM,JM),GRID2(IM,JM),IW(IM,JM),IWM1
!      REAL SLEET(IM,JM),RAIN(IM,JM),FREEZR(IM,JM),SNOW(IM,JM)
      REAL SLEET(IM,JM,NALG),RAIN(IM,JM,NALG),
     *     FREEZR(IM,JM,NALG),SNOW(IM,JM,NALG)
      REAL SLEET1(IM,JM),RAIN1(IM,JM),FREEZR1(IM,JM)
     *    ,SNOW1(IM,JM)
      REAL SLEET2(IM,JM),RAIN2(IM,JM),FREEZR2(IM,JM)
     *    ,SNOW2(IM,JM)
      REAL SLEET3(IM,JM),RAIN3(IM,JM),FREEZR3(IM,JM)
     *    ,SNOW3(IM,JM)
      REAL SLEET4(IM,JM),RAIN4(IM,JM),FREEZR4(IM,JM)
     *    ,SNOW4(IM,JM)
      REAL SLEET5(IM,JM),RAIN5(IM,JM),FREEZR5(IM,JM)
     *    ,SNOW5(IM,JM)
      REAL DOMS(IM,JM),DOMR(IM,JM),DOMIP(IM,JM),DOMZR(IM,JM)

      REAL ECAN(IM,JM),EDIR(IM,JM),ETRANS(IM,JM),ESNOW(IM,JM)
     &,SMCDRY(IM,JM),SMCMAX(IM,JM)
      REAL RSMIN(IM,JM),SMCWLT(IM,JM),SMCREF(IM,JM)
     & ,RCS(IM,JM),RCQ(IM,JM),RCT(IM,JM),RCSOIL(IM,JM)
     & ,GC(IM,JM)     
C     
C****************************************************************************
C
C     START SURFCE.
C
C     
C***  BLOCK 1.  SURFACE BASED FIELDS.
C
C     IF ANY OF THE FOLLOWING "SURFACE" FIELDS ARE REQUESTED,
C     WE NEED TO COMPUTE THE FIELDS FIRST.
C     
      IF ( (IGET(024).GT.0).OR.(IGET(025).GT.0).OR.
     X     (IGET(026).GT.0).OR.(IGET(027).GT.0).OR.
     X     (IGET(028).GT.0).OR.(IGET(029).GT.0).OR.
     X     (IGET(154).GT.0).OR.
     X     (IGET(034).GT.0).OR.(IGET(076).GT.0) ) THEN
C     
         DO J=JSTA,JEND
         DO I=1,IM
C
C           SCALE ARRAY FIS BY GI TO GET SURFACE HEIGHT.
            ZSFC(I,J)=FIS(I,J)*GI
C
C           SURFACE PRESSURE.
            PSFC(I,J)=PINT(I,J,NINT(LMH(I,J))+1)
C     
C           SURFACE (SKIN) POTENTIAL TEMPERATURE AND TEMPERATURE.
            THSFC(I,J)=THS(I,J)
            TSFC(I,J) =THSFC(I,J)*(PSFC(I,J)/P1000)**CAPA 
C     
C           SURFACE SPECIFIC HUMIDITY, RELATIVE HUMIDITY,
C           AND DEWPOINT.  ADJUST SPECIFIC HUMIDITY IF
C           RELATIVE HUMIDITY EXCEEDS 0.1 OR 1.0.
C
            QSFC(I,J)=QS(I,J)
            QSFC(I,J)=AMAX1(H1M12,QSFC(I,J))
            TSFCK    =TSFC(I,J)
C     
            QSAT=PQ0/PSFC(I,J)
     1          *EXP(A2*(TSFCK-A3)/(TSFCK-A4))
            RHSFC(I,J)=QSFC(I,J)/QSAT

            IF (RHSFC(I,J).GT.H1 ) RHSFC(I,J) = H1
            IF (RHSFC(I,J).LT.D00) RHSFC(I,J) = D01
            QSFC(I,J)  = RHSFC(I,J)*QSAT
            EVP(I,J)   = PSFC(I,J)*QSFC(I,J)/(EPS+ONEPS*QSFC(I,J))
            EVP(I,J)   = EVP(I,J)*D001
C     
Cmp           ACCUMULATED NON-CONVECTIVE PRECIP.
Cmp            IF(IGET(034).GT.0)THEN
Cmp              IF(LVLS(1,IGET(034)).GT.0)THEN

C           ACCUMULATED PRECIP (convective + non-convective)
            IF(IGET(087).GT.0)THEN
              IF(LVLS(1,IGET(087)).GT.0)THEN
C	write(6,*) 'acprec, ancprc, cuprec: ', ANCPRC(I,J)+CUPREC(I,J),
C     +		ANCPRC(I,J),CUPREC(I,J)
!                 ACPREC(I,J)=ANCPRC(I,J)+CUPREC(I,J)
              ENDIF
            ENDIF

         ENDDO
         ENDDO
C     
C        INTERPOLATE/OUTPUT REQUESTED SURFACE FIELDS.
C     
C        SURFACE PRESSURE.
         IF (IGET(024).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PSFC(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(024),LVLS(1,IGET(024)),
     X           GRID1,IM,JM)
         ENDIF
C     
C        SURFACE HEIGHT.
         IF (IGET(025).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ZSFC(I,J)
            ENDDO
            ENDDO
!            CALL BOUND(GRID1,D00,H99999)
            ID(1:25) = 0
            CALL GRIBIT(IGET(025),LVLS(1,IGET(025)),
     X           GRID1,IM,JM)
         ENDIF
C     
C        SURFACE (SKIN) TEMPERATURE.
         IF (IGET(026).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=TSFC(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(026),LVLS(1,IGET(026)),
     X           GRID1,IM,JM)
         ENDIF
C     
C        SURFACE (SKIN) POTENTIAL TEMPERATURE.
         IF (IGET(027).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=THSFC(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(027),LVLS(1,IGET(027)),
     X           GRID1,IM,JM)
         ENDIF
C     
C        SURFACE SPECIFIC HUMIDITY.
         IF (IGET(028).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=QSFC(I,J)
            ENDDO
            ENDDO
            CALL BOUND(GRID1,H1M12,H99999)
            ID(1:25) = 0
            CALL GRIBIT(IGET(028),LVLS(1,IGET(028)),
     X           GRID1,IM,JM)
         ENDIF
C     
C        SURFACE DEWPOINT TEMPERATURE.
         IF (IGET(029).GT.0) THEN
            CALL DEWPOINT(EVP,DWPSFC)
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=DWPSFC(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(029),LVLS(1,IGET(029)),
     X           GRID1,IM,JM)
         ENDIF
C     
C        SURFACE RELATIVE HUMIDITY.
         IF (IGET(076).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=RHSFC(I,J)*100.
            ENDDO
            ENDDO
            CALL BOUND(GRID1,H1,H100)
            ID(1:25) = 0
            CALL GRIBIT(IGET(076),LVLS(1,IGET(076)),
     X           GRID1,IM,JM)
         ENDIF
C     
      ENDIF
C
C     ADDITIONAL SURFACE-SOIL LEVEL FIELDS.
C

      DO L=1,NSOIL
C     SOIL TEMPERATURE.
      IF (IGET(116).GT.0) THEN
        IF (LVLS(L,IGET(116)).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=STC(I,J,L)
            ENDDO
            ENDDO
         ID(1:25) = 0
          DTOP=0.
          DO LS=1,L-1
           DTOP=DTOP+SLDPTH(LS)
          ENDDO
          DBOT=DTOP+SLDPTH(L)
            ID(10) = NINT(DTOP*100.)
            ID(11) = NINT(DBOT*100.)
         CALL GRIBIT(IGET(116),L,GRID1,IM,JM)
        ENDIF
      ENDIF
C
C     SOIL MOISTURE.
      IF (IGET(117).GT.0) THEN
        IF (LVLS(L,IGET(117)).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SMC(I,J,L)
            ENDDO
            ENDDO
         ID(1:25) = 0
          DTOP=0.
          DO LS=1,L-1
           DTOP=DTOP+SLDPTH(LS)
          ENDDO
          DBOT=DTOP+SLDPTH(L)
            ID(10) = NINT(DTOP*100.)
            ID(11) = NINT(DBOT*100.)
         CALL GRIBIT(IGET(117),L,GRID1,IM,JM)
        ENDIF
      ENDIF
c      ENDDO
C     ADD LIQUID SOIL MOISTURE
      IF (IGET(225).GT.0) THEN
        IF (LVLS(L,IGET(225)).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SH2O(I,J,L)
            ENDDO
            ENDDO
         ID(1:25) = 0
          DTOP=0.
          DO LS=1,L-1
           DTOP=DTOP+SLDPTH(LS)
          ENDDO
          DBOT=DTOP+SLDPTH(L)
            ID(10) = NINT(DTOP*100.)
            ID(11) = NINT(DBOT*100.)
            ID(02) = 130
         CALL GRIBIT(IGET(225),L,GRID1,IM,JM)
        ENDIF
      ENDIF
! END OF NSOIL LOOP
      ENDDO
C
C     BOTTOM SOIL TEMPERATURE.
      IF (IGET(115).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SOILTB(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         ISVALUE     = 300
         ID(11) = ISVALUE
         CALL GRIBIT(IGET(115),LVLS(1,IGET(115)),
     X        GRID1,IM,JM)
      ENDIF
C
C     SOIL MOISTURE AVAILABILITY
      IF (IGET(171).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SMSTAV(I,J)*100.
            ENDDO
            ENDDO
         ID(1:25) = 0
         ID(10) = 0
         ID(11) = 100
         CALL GRIBIT(IGET(171),LVLS(1,IGET(171)),
     X        GRID1,IM,JM)
      ENDIF
C
C     TOTAL SOIL MOISTURE
      IF (IGET(036).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
	     IF(SM(I,J).GT.SMALL .AND. SICE(I,J).LT.SMALL)THEN
	      GRID1(I,J)=1000.0  ! TEMPORY FIX TO MAKE SURE SMSTOT=1 FOR WATER
	     ELSE  
              GRID1(I,J)=SMSTOT(I,J)
	     END IF 
            ENDDO
            ENDDO
         ID(1:25) = 0
         ID(10) = 0
         ID(11) = 200
         CALL GRIBIT(IGET(036),LVLS(1,IGET(036)),
     X        GRID1,IM,JM)
      ENDIF
C
C     PLANT CANOPY SURFACE WATER.
      IF ( IGET(118).GT.0 ) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=CMC(I,J)*1000.
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(118),LVLS(1,IGET(118)),
     X        GRID1,IM,JM)
      ENDIF
C
C     SNOW WATER EQUIVALENT.
      IF ( IGET(119).GT.0 ) THEN
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=SNO(I,J)*1000.
             GRID1(I,J)=SNO(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(119),LVLS(1,IGET(119)),
     X        GRID1,IM,JM)
      ENDIF
C
C     PERCENT SNOW COVER.
      IF ( IGET(120).GT.0 ) THEN
         DO J=JSTA,JEND
         DO I=1,IM
!             GRID1(I,J)=PCTSNO(I,J)
           SNEQV=SNO(I,J)
           IVEG=IVGTYP(I,J)
           IF(IVEG.EQ.0)IVEG=7
           CALL SNFRAC (SNEQV,IVEG,SNCOVR)
           GRID1(I,J)=SNCOVR*100.
         ENDDO
         ENDDO
         CALL BOUND(GRID1,D00,H100)
         ID(1:25) = 0
         CALL GRIBIT(IGET(120),LVLS(1,IGET(120)),
     X        GRID1,IM,JM)
      ENDIF
! ADD SNOW DEPTH
      IF ( IGET(224).GT.0 ) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SI(I,J)/1000.  ! SI comes out of WRF in mm
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(224),LVLS(1,IGET(224)),
     X        GRID1,IM,JM)
      ENDIF      
! ADD POTENTIAL EVAPORATION
      IF ( IGET(242).GT.0 ) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=POTEVP(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(242),LVLS(1,IGET(242)),
     X        GRID1,IM,JM)
      ENDIF

! ADD EC,EDIR,ETRANS,ESNOW,SMCDRY,SMCMAX
! ONLY OUTPUT NEW LSM FIELDS FOR NMM AND ARW BECAUSE RSM USES OLD SOIL TYPES
      IF (MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'NMM')THEN
      IF ( IGET(228).GT.0 .OR. IGET(229).GT.0
     + .OR.IGET(230).GT.0 .OR. IGET(231).GT.0
     + .OR.IGET(232).GT.0 .OR. IGET(233).GT.0) THEN
        DO J=JSTA,JEND
         DO I=1,IM
C ----------------------------------------------------------------------
!          IF(QWBS(I,J).gt.0.001)print*,'NONZERO QWBS',i,j,QWBS(I,J)
!          IF(abs(SM(I,J)-0.).lt.1.0E-5)THEN
          IF( (abs(SM(I,J)-0.)   .lt. 1.0E-5) .AND.
     &        (abs(SICE(I,J)-0.) .lt. 1.0E-5) ) THEN
           CALL ETCALC(QWBS(I,J),POTEVP(I,J),SNO(I,J),VEGFRC(I,J)
     &  ,  ISLTYP(I,J),SH2O(I,J,1:1),CMC(I,J)
     &  ,  ECAN(I,J),EDIR(I,J),ETRANS(I,J),ESNOW(I,J),SMCDRY(I,J)
     &  ,  SMCMAX(I,J) )
          ELSE
           ECAN(I,J)=0.
           EDIR(I,J)=0.
           ETRANS(I,J)=0.
           ESNOW(I,J)=0.
           SMCDRY(I,J)=0.
           SMCMAX(I,J)=0.
          END IF
         END DO
        END DO

        IF ( IGET(228).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ECAN(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(228),LVLS(1,IGET(228)),
     X        GRID1,IM,JM)
        ENDIF	

        IF ( IGET(229).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=EDIR(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(229),LVLS(1,IGET(229)),
     X        GRID1,IM,JM)
        ENDIF

        IF ( IGET(230).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ETRANS(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(230),LVLS(1,IGET(230)),
     X        GRID1,IM,JM)
        ENDIF
	

        IF ( IGET(231).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ESNOW(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
	 ID(02)= 130
         CALL GRIBIT(IGET(231),LVLS(1,IGET(231)),
     X        GRID1,IM,JM)
        ENDIF	

        IF ( IGET(232).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SMCDRY(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
	 ID(02)= 130
         CALL GRIBIT(IGET(232),LVLS(1,IGET(232)),
     X        GRID1,IM,JM)
        ENDIF

        IF ( IGET(233).GT.0 )THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SMCMAX(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
	 ID(02)= 130
         CALL GRIBIT(IGET(233),LVLS(1,IGET(233)),
     X        GRID1,IM,JM)
        ENDIF
	
      ENDIF
      END IF  ! endif for ncar and nmm options
C
C     
C
C***  BLOCK 2.  SHELTER (2M) LEVEL FIELDS.
C     
C     COMPUTE/POST SHELTER LEVEL FIELDS.
C     
      IF ( (IGET(106).GT.0).OR.(IGET(112).GT.0).OR.
     X     (IGET(113).GT.0).OR.(IGET(114).GT.0).OR.
     X     (IGET(138).GT.0) ) THEN
C
CHC  COMPUTE SHELTER PRESSURE BECAUSE IT WAS NOT OUTPUT FROM WRF       
        IF(MODELNAME .EQ. 'NCAR' .OR. MODELNAME.EQ.'RSM')THEN
         DO J=JSTA,JEND
         DO I=1,IM
          TLOW=T(I,J,NINT(LMH(I,J)))
          PSHLTR(I,J)=PSFC(I,J)*EXP(-0.068283/TLOW)
         END DO
         END DO 
	END IF 
C
C        SHELTER LEVEL TEMPERATURE
         IF (IGET(106).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
c             GRID1(I,J)=TSHLTR(I,J)
CHC CONVERT FROM THETA TO T 
             GRID1(I,J)=TSHLTR(I,J)*(PSHLTR(I,J)*1.E-5)**CAPA
             IF(GRID1(I,J).LT.200)PRINT*,'ABNORMAL 2MT ',i,j,
     +       TSHLTR(I,J),PSHLTR(I,J)
!             TSHLTR(I,J)=GRID1(I,J) 
            ENDDO
            ENDDO
            ID(1:25) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(106),LVLS(1,IGET(106)),
     X           GRID1,IM,JM)
         ENDIF
C
C        SHELTER LEVEL SPECIFIC HUMIDITY.
         IF (IGET(112).GT.0) THEN       
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=QSHLTR(I,J)
            ENDDO
            ENDDO
            CALL BOUND (GRID1,H1M12,H99999)
            ID(1:25) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(112),LVLS(1,IGET(112)),
     X           GRID1,IM,JM)
         ENDIF
C     
C        SHELTER LEVEL DEWPOINT.
         IF (IGET(113).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
              EVP(I,J)=PSHLTR(I,J)*QSHLTR(I,J)/(EPS+ONEPS*QSHLTR(I,J))
	      EVP(I,J)=EVP(I,J)*D001
            ENDDO
            ENDDO
            CALL DEWPOINT(EVP,EGRID1)
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=EGRID1(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(113),LVLS(1,IGET(113)),
     X           GRID1,IM,JM)
         ENDIF
C     
C        SHELTER LEVEL RELATIVE HUMIDITY.
         IF (IGET(114).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             P1D(I,J)=PSHLTR(I,J)
             T1D(I,J)=TSHLTR(I,J)*(PSHLTR(I,J)*1.E-5)**CAPA
             Q1D(I,J)=QSHLTR(I,J)
            ENDDO
            ENDDO
!            CALL CALRH(PSHLTR,TSHLTR,QSHLTR,EGRID1)
            CALL CALRH(P1D,T1D,Q1D,EGRID1)
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=EGRID1(I,J)*100.
            ENDDO
            ENDDO
            CALL BOUND(GRID1,H1,H100)
            ID(1:25) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(114),LVLS(1,IGET(114)),
     X           GRID1,IM,JM)
         ENDIF
C     
C        SHELTER LEVEL PRESSURE.
         IF (IGET(138).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PSHLTR(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ISVALUE = 2
            ID(10) = MOD(ISVALUE/256,256)
            ID(11) = MOD(ISVALUE,256)
            CALL GRIBIT(IGET(138),LVLS(1,IGET(138)),
     X           GRID1,IM,JM)
         ENDIF
C
      ENDIF
C
C
C     BLOCK 3.  ANEMOMETER LEVEL (10M) WINDS, THETA, AND Q.
C
      IF ( (IGET(064).GT.0).OR.(IGET(065).GT.0) ) THEN
         ID(1:25) = 0
         ISVALUE = 10
         ID(10) = MOD(ISVALUE/256,256)
         ID(11) = MOD(ISVALUE,256)
C
C        ANEMOMETER LEVEL U WIND AND/OR V WIND.
         IF ((IGET(064).GT.0).OR.(IGET(065).GT.0)) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=U10(I,J)
             GRID2(I,J)=V10(I,J)
            ENDDO
            ENDDO
            IF (IGET(064).GT.0) CALL GRIBIT(IGET(064),
     X           LVLS(1,IGET(064)),GRID1,IM,JM)
            IF (IGET(065).GT.0) CALL GRIBIT(IGET(065),
     X           LVLS(1,IGET(065)),GRID2,IM,JM)
         ENDIF
      ENDIF
C
C        ANEMOMETER LEVEL (10 M) POTENTIAL TEMPERATURE.
C   NOT A OUTPUT FROM WRF
      IF (IGET(158).GT.0) THEN
         ID(1:25) = 0
         ISVALUE = 10
         ID(10) = MOD(ISVALUE/256,256)
         ID(11) = MOD(ISVALUE,256)
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=TH10(I,J)
            ENDDO
            ENDDO
         CALL GRIBIT(IGET(158),
     X        LVLS(1,IGET(158)),GRID1,IM,JM)
       ENDIF
C
C        ANEMOMETER LEVEL (10 M) SPECIFIC HUMIDITY.
C
      IF (IGET(159).GT.0) THEN
         ID(1:25) = 0
         ISVALUE = 10
         ID(10) = MOD(ISVALUE/256,256)
         ID(11) = MOD(ISVALUE,256)
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=Q10(I,J)
            ENDDO
            ENDDO
         CALL GRIBIT(IGET(159),
     X        LVLS(1,IGET(159)),GRID1,IM,JM)
       ENDIF
C
C
C
C***  BLOCK 4.  PRECIPITATION RELATED FIELDS.
!MEB 6/17/02  ASSUMING THAT ALL ACCUMULATED FIELDS NEVER EMPTY
!             THEIR BUCKETS.  THIS IS THE EASIEST WAY TO DEAL WITH
!             ACCUMULATED FIELDS.  SHORTER TIME ACCUMULATIONS CAN
!             BE COMPUTED AFTER THE FACT IN A SEPARATE CODE ONCE
!             THE POST HAS FINISHED.  I HAVE LEFT IN THE OLD
!             ETAPOST CODE FOR COMPUTING THE BEGINNING TIME OF
!             THE ACCUMULATION PERIOD IF THIS IS CHANGED BACK
!             TO A 12H OR 3H BUCKET.  I AM NOT SURE WHAT
!             TO DO WITH THE TIME AVERAGED FIELDS, SO
!             LEAVING THAT UNCHANGED.
C     
C     SNOW FRACTION FROM EXPLICIT CLOUD SCHEME.  LABELLED AS
C      'PROB OF FROZEN PRECIP' IN GRIB, 
C      DIDN'T KNOW WHAT ELSE TO CALL IT
      IF (IGET(172).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
              IF (PREC(I,J) .LE. PTHRESH) THEN
                GRID1(I,J)=-50.
              ELSE
                GRID1(I,J)=SR(I,J)*100.
              ENDIF
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(172),LVLS(1,IGET(172)),
     X        GRID1,IM,JM)
      ENDIF
C     INSTANTANEOUS PRECIPITATION RATE.
      IF (IGET(167).GT.0) THEN
!MEB need to get physics DT
         RDTPHS=1./(DT * NPHS) 
!MEB need to get physics DT
            DO J=JSTA,JEND
            DO I=1,IM
             IF(MODELNAME .EQ. 'NCAR')THEN
              GRID1(I,J)=PREC(I,J)/DT*1000.
             ELSE IF (MODELNAME .EQ. 'NMM')THEN
              GRID1(I,J)=PREC(I,J)*RDTPHS*1000.
             ELSE IF (MODELNAME .EQ. 'RSM') THEN    !Add by Binbin 
              GRID1(I,J)=PREC(I,J)
             END IF
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(167),LVLS(1,IGET(167)),
     X        GRID1,IM,JM)
      ENDIF
C
C     INSTANTANEOUS CONVECTIVE PRECIPITATION RATE.
C     SUBSTITUTE WITH CUPPT IN WRF FOR NOW
      IF (IGET(249).GT.0) THEN
         RDTPHS=1000./DTQ2     !--- 1000 kg/m**3, density of liquid water
!         RDTPHS=1000./(TRDLW*3600.)
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J)=CPRATE(I,J)*RDTPHS
!           GRID1(I,J)=CUPPT(I,J)*RDTPHS
         ENDDO
         ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(249),LVLS(1,IGET(249)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     ACCUMULATED TOTAL PRECIPITATION.
      IF (IGET(087).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ACPREC(I,J)*1000.
            ENDDO
            ENDDO
         ID(1:25) = 0
         ITPREC     = INT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
	 IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
         ID(18)     = 0
         ID(19)     = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)     = 4
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITPREC
         ELSE
          ID(18) = IFHR-IFINCR
	  IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
!	 IF(IFMIN .GE. 1 .AND. ID(19) .GT. 256)THEN
!	  IF(ITPREC.EQ.3)ID(17)=10
!	  IF(ITPREC.EQ.6)ID(17)=11
!	  IF(ITPREC.EQ.12)ID(17)=12
!	 END IF 
         IF (ID(18).LT.0) ID(18) = 0
!	write(6,*) 'call gribit...total precip'
         CALL GRIBIT(IGET(087),LVLS(1,IGET(087)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     ACCUMULATED CONVECTIVE PRECIPITATION.
      IF (IGET(033).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=CUPREC(I,J)*1000.
            ENDDO
            ENDDO
         ID(1:25) = 0
         ITPREC     = INT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
         ID(18)     = 0
         ID(19)     = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)     = 4
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITPREC
         ELSE
          ID(18) = IFHR-IFINCR
          IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
!	write(6,*) 'call gribit...convective precip'
         CALL GRIBIT(IGET(033),LVLS(1,IGET(033)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     ACCUMULATED GRID-SCALE PRECIPITATION.
      IF (IGET(034).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=ANCPRC(I,J)*1000.
            ENDDO
            ENDDO
         ID(1:25) = 0
         ITPREC     = INT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
         ID(18)     = 0
         ID(19)     = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)     = 4
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITPREC
         ELSE
          ID(18) = IFHR-IFINCR
          IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
!	write(6,*) 'call gribit...grid-scale precip'
         CALL GRIBIT(IGET(034),LVLS(1,IGET(034)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     ACCUMULATED LAND SURFACE PRECIPITATION.
      IF (IGET(256).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=LSPA(I,J)*1000.
            ENDDO
            ENDDO
         ID(1:25) = 0
         ITPREC     = INT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
         ID(18)     = 0
         ID(19)     = IFHR
	 IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20)     = 4
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITPREC
         ELSE
          ID(18) = IFHR-IFINCR
          IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
         ENDIF
         IF (ID(18).LT.0) ID(18) = 0
         ID(02)= 130
         CALL GRIBIT(IGET(256),LVLS(1,IGET(256)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     ACCUMULATED SNOWFALL.
         IF (IGET(035).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=ACSNOW(I,J)*1000.
             GRID1(I,J)=ACSNOW(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = INT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(035),LVLS(1,IGET(035)),
     X           GRID1,IM,JM)
         ENDIF
C     
C     ACCUMULATED SNOW MELT.
         IF (IGET(121).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=ACSNOM(I,J)*1000.
             GRID1(I,J)=ACSNOM(I,J)	     
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = INT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(121),LVLS(1,IGET(121)),
     X           GRID1,IM,JM)
         ENDIF
C     
C     ACCUMULATED STORM SURFACE RUNOFF.
         IF (IGET(122).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=SSROFF(I,J)*1000.
             GRID1(I,J)=SSROFF(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = INT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(122),LVLS(1,IGET(122)),
     X           GRID1,IM,JM)
         ENDIF
C     
C     ACCUMULATED BASEFLOW-GROUNDWATER RUNOFF.
         IF (IGET(123).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
!             GRID1(I,J)=BGROFF(I,J)*1000.
             GRID1(I,J)=BGROFF(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = INT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
         IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
             IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(123),LVLS(1,IGET(123)),
     X           GRID1,IM,JM)
         ENDIF
C     
C     INSTANTANEOUS PRECIPITATION TYPE.
         IF (IGET(160).GT.0 .OR.(IGET(247).GT.0)) THEN

          CALL CALWXT(T,Q,PMID,PINT,HTM,LMH,PREC,ZINT,IWX1
     X       ,ZWET)
          IF (IGET(160).GT.0) THEN 
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=IWX1(I,J)
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW1(I,J)   = ISNO*1.0
              SLEET1(I,J)  = IIP*1.0
              FREEZR1(I,J) = IZR*1.0
              RAIN1(I,J)   = IRAIN*1.0
	      SNOW(I,J,1)    = SNOW1(I,J)
              SLEET(I,J,1)   = SLEET1(I,J)
              FREEZR(I,J,1) = FREEZR1(I,J)
              RAIN(I,J,1)    = RAIN1(I,J)
            ENDDO
            ENDDO
          ENDIF
C     
C     LOWEST WET BULB ZERO HEIGHT
           IF (IGET(247).GT.0) THEN
             DO J=JSTA,JEND
             DO I=1,IM
              GRID1(I,J)=ZWET(I,J)
             ENDDO
             ENDDO
             ID(1:25) = 0
             CALL GRIBIT(IGET(247),LVLS(1,IGET(247)),
     X           GRID1,IM,JM)
           ENDIF

C     DOMINANT PRECIPITATION TYPE
CGSM  IF DOMINANT PRECIP TYPE IS REQUESTED, 4 MORE ALGORITHMS
CGSM    WILL BE CALLED.  THE TALLIES ARE THEN SUMMED IN
CGSM    CALWXT_DOMINANT

           IF (IGET(160).GT.0) THEN   
C  RAMER ALGORITHM
            CALL CALWXT_RAMER(T,Q,PMID,PINT,LMH,PREC,IWX2)
               
C     DECOMPOSE IWX2 ARRAY
C
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=NINT(IWX2(I,J))
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW2(I,J)   = ISNO*1.0
              SLEET2(I,J)  = IIP*1.0
              FREEZR2(I,J) = IZR*1.0
              RAIN2(I,J)   = IRAIN*1.0
              SNOW(I,J,2)    = SNOW2(I,J)
              SLEET(I,J,2)   = SLEET2(I,J)
              FREEZR(I,J,2) = FREEZR2(I,J)
              RAIN(I,J,2)    = RAIN2(I,J)
            ENDDO
            ENDDO

C BOURGOUIN ALGORITHM
            CALL CALWXT_BOURG(T,Q,PMID,PINT,LMH,PREC,ZINT,IWX3)

C     DECOMPOSE IWX3 ARRAY
C
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=NINT(IWX3(I,J))
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW3(I,J)   = ISNO*1.0
              SLEET3(I,J)  = IIP*1.0
              FREEZR3(I,J) = IZR*1.0
              RAIN3(I,J)   = IRAIN*1.0
              SNOW(I,J,3)    = SNOW3(I,J)
              SLEET(I,J,3)   = SLEET3(I,J)
              FREEZR(I,J,3) = FREEZR3(I,J)
              RAIN(I,J,3)    = RAIN3(I,J)
            ENDDO
            ENDDO

C REVISED NCEP ALGORITHM
            CALL CALWXT_REVISED(T,Q,PMID,PINT,HTM,LMH,PREC,ZINT,
     X          IWX4)
C     DECOMPOSE IWX2 ARRAY
C
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=IWX4(I,J)
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW4(I,J)   = ISNO*1.0
              SLEET4(I,J)  = IIP*1.0
              FREEZR4(I,J) = IZR*1.0
              RAIN4(I,J)   = IRAIN*1.0
              SNOW(I,J,4)    = SNOW4(I,J)
              SLEET(I,J,4)   = SLEET4(I,J)
              FREEZR(I,J,4) = FREEZR4(I,J)
              RAIN(I,J,4)    = RAIN4(I,J)
            ENDDO
            ENDDO
              
C EXPLICIT ALGORITHM (UNDER 18 NOT ADMITTED WITHOUT PARENT 
C     OR GUARDIAN)
 
            CALL CALWXT_EXPLICIT(LMH,THS,PMID,PREC,SR,F_RimeF,IWX5)
C     DECOMPOSE IWX2 ARRAY
C
            DO J=JSTA,JEND
            DO I=1,IM
              IWX=IWX5(I,J)
              ISNO=MOD(IWX,2)
              IIP=MOD(IWX,4)/2
              IZR=MOD(IWX,8)/4
              IRAIN=IWX/8
              SNOW5(I,J)   = ISNO*1.0
              SLEET5(I,J)  = IIP*1.0
              FREEZR5(I,J) = IZR*1.0
              RAIN5(I,J)   = IRAIN*1.0
              SNOW(I,J,5)    = SNOW5(I,J)
              SLEET(I,J,5)   = SLEET5(I,J)
              FREEZR(I,J,5) = FREEZR5(I,J)
              RAIN(I,J,5)    = RAIN5(I,J)
            ENDDO
            ENDDO
               
           CALL CALWXT_DOMINANT(PREC,RAIN,FREEZR,SLEET,SNOW,
     X         DOMR,DOMZR,DOMIP,DOMS)
           ID(1:25) = 0
C     SNOW.
            ID(8) = 143 
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=DOMS(I,J)
            ENDDO
            ENDDO
            CALL GRIBIT(IGET(160),LVLS(1,IGET(160)),
     X           GRID1,IM,JM)
C     ICE PELLETS.
            ID(8) = 142 
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=DOMIP(I,J)
            ENDDO
            ENDDO
            CALL GRIBIT(IGET(160),LVLS(1,IGET(160)),
     X           GRID1,IM,JM)
C     FREEZING RAIN.
            ID(8) = 141 
            DO J=JSTA,JEND
            DO I=1,IM
!             if (DOMZR(I,J) .EQ. 1) THEN
!               PSFC(I,J)=PINT(I,J,NINT(LMH(I,J))+1)
!               print *, 'aha ', I, J, PSFC(I,J)
!               print *, FREEZR(I,J,1), FREEZR(I,J,2),
!     *  FREEZR(I,J,3), FREEZR(I,J,4), FREEZR(I,J,5)
!             endif
             GRID1(I,J)=DOMZR(I,J)
            ENDDO
            ENDDO
            CALL GRIBIT(IGET(160),LVLS(1,IGET(160)),
     X           GRID1,IM,JM)
C     RAIN.
            ID(8) = 140 
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=DOMR(I,J)
            ENDDO
            ENDDO
	    CALL GRIBIT(IGET(160),LVLS(1,IGET(160)),
     X           GRID1,IM,JM)
        ENDIF
      ENDIF
C     
C
C
C***  BLOCK 5.  SURFACE EXCHANGE FIELDS.
C     
C     TIME AVERAGED SURFACE LATENT HEAT FLUX.
         IF (IGET(042).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE  
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=-1.*SFCLHX(I,J)*RRNUM !change the sign to conform with Grib
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
	    IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(042),LVLS(1,IGET(042)),
     X           GRID1,IM,JM)
          END IF 
         ENDIF
C
C     TIME AVERAGED SURFACE SENSIBLE HEAT FLUX.
         IF (IGET(043).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J) = -1.* SFCSHX(I,J)*RRNUM !change the sign to conform with Grib
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
	    IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(043),LVLS(1,IGET(043)),
     X           GRID1,IM,JM)
         ENDIF
C     
C     TIME AVERAGED SUB-SURFACE SENSIBLE HEAT FLUX.
         IF (IGET(135).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J) = SUBSHX(I,J)*RRNUM
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
            IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(135),LVLS(1,IGET(135)),
     X           GRID1,IM,JM)
         ENDIF
C     
C     TIME AVERAGED SNOW PHASE CHANGE HEAT FLUX.
         IF (IGET(136).GT.0) THEN
          IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J) = SNOPCX(I,J)*RRNUM
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
            IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(136),LVLS(1,IGET(136)),
     X           GRID1,IM,JM)
         ENDIF
C     
C     TIME AVERAGED SURFACE MOMENTUM FLUX.
         IF (IGET(046).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
            IF(ASRFC.GT.0.)THEN
              RRNUM=1./ASRFC
            ELSE
              RRNUM=0.
            ENDIF
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J) = SFCUVX(I,J)*RRNUM
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITSRFC     = INT(TSRFC)
            IF(ITSRFC .ne. 0) then
             IFINCR     = MOD(IFHR,ITSRFC)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITSRFC*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 3
            IF (IFINCR.EQ.0) THEN
               ID(18) = IFHR-ITSRFC
            ELSE
               ID(18) = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(046),LVLS(1,IGET(046)),
     X           GRID1,IM,JM)
         ENDIF
C     
C     ACCUMULATED SURFACE EVAPORATION
         IF (IGET(047).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SFCEVP(I,J)*1000.
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = INT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
	 IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
	     IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(047),LVLS(1,IGET(047)),
     X           GRID1,IM,JM)
         ENDIF
C     
C     ACCUMULATED POTENTIAL EVAPORATION
         IF (IGET(137).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=POTEVP(I,J)*1000.
            ENDDO
            ENDDO
            ID(1:25) = 0
            ITPREC     = INT(TPREC)
!mp
	if (ITPREC .ne. 0) then
         IFINCR     = MOD(IFHR,ITPREC)
	 IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITPREC*60)
	else
	 IFINCR     = 0
	endif
!mp
            ID(18)     = 0
            ID(19)     = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)     = 4
            IF (IFINCR.EQ.0) THEN
             ID(18) = IFHR-ITPREC
            ELSE
             ID(18) = IFHR-IFINCR
	     IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
            CALL GRIBIT(IGET(137),LVLS(1,IGET(137)),
     X           GRID1,IM,JM)
         ENDIF
C     
C     ROUGHNESS LENGTH.
      IF (IGET(044).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=Z0(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(044),LVLS(1,IGET(044)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     FRICTION VELOCITY.
      IF (IGET(045).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=USTAR(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(045),LVLS(1,IGET(045)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     SURFACE DRAG COEFFICIENT.
      IF (IGET(132).GT.0) THEN
         CALL CALDRG(EGRID1)
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=EGRID1(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(132),LVLS(1,IGET(132)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     SURFACE U AND/OR V COMPONENT WIND STRESS
      IF ( (IGET(133).GT.0) .OR. (IGET(134).GT.0) ) THEN
         CALL CALTAU(EGRID1,EGRID2)
C     
C        SURFACE U COMPONENT WIND STRESS.
         IF (IGET(133).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=EGRID1(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(133),LVLS(1,IGET(133)),
     X           GRID1,IM,JM)
         ENDIF
C     
C        SURFACE V COMPONENT WIND STRESS
         IF (IGET(134).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=EGRID2(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(134),LVLS(1,IGET(134)),
     X           GRID1,IM,JM)
         ENDIF
      ENDIF
C     
C     INSTANTANEOUS SENSIBLE HEAT FLUX
      IF (IGET(154).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=-1.*TWBS(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(154),LVLS(1,IGET(154)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     INSTANTANEOUS LATENT HEAT FLUX
      IF (IGET(155).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=-1.*QWBS(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(155),LVLS(1,IGET(155)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     SURFACE EXCHANGE COEFF
      IF (IGET(169).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=SFCEXC(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(169),LVLS(1,IGET(169)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     GREEN VEG FRACTION
      IF (IGET(170).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=VEGFRC(I,J)*100.
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(170),LVLS(1,IGET(170)),
     X        GRID1,IM,JM)
      ENDIF
C     
C     INSTANTANEOUS GROUND HEAT FLUX
      IF (IGET(152).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=GRNFLX(I,J)
            ENDDO
            ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(152),LVLS(1,IGET(152)),
     X        GRID1,IM,JM)
      ENDIF
!    VEGETATION TYPE
      IF (IGET(218).GT.0) THEN
         DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = FLOAT(IVGTYP(I,J))
           ENDDO
         ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(218),LVLS(1,IGET(218)),
     X        GRID1,IM,JM)                                                          
      ENDIF
!
!    SOIL TYPE
      IF (IGET(219).GT.0) THEN
         DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = FLOAT(ISLTYP(I,J))
           ENDDO
         ENDDO
         ID(1:25) = 0
         CALL GRIBIT(IGET(219),LVLS(1,IGET(219)),
     X        GRID1,IM,JM)                                                          
      ENDIF
!    SLOPE TYPE
      IF (IGET(223).GT.0) THEN
         DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = FLOAT(ISLOPE(I,J))                                  
           ENDDO
         ENDDO
         ID(1:25) = 0
         ID(02)= 130
         CALL GRIBIT(IGET(223),LVLS(1,IGET(223)),
     X        GRID1,IM,JM)
                                                                                
      ENDIF
!      print*,'starting computing canopy conductance'
!
! CANOPY CONDUCTANCE
! ONLY OUTPUT NEW LSM FIELDS FOR NMM AND ARW BECAUSE RSM USES OLD SOIL TYPES
      IF (MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'NMM')THEN
      IF (IGET(220).GT.0 .OR. IGET(234).GT.0
     & .OR. IGET(235).GT.0 .OR. IGET(236).GT.0
     & .OR. IGET(237).GT.0 .OR. IGET(238).GT.0
     & .OR. IGET(239).GT.0 .OR. IGET(240).GT.0
     & .OR. IGET(241).GT.0 .OR. IGET(254).GT.0 ) THEN
      print*,'starting computing canopy conductance'     
         DO J=JSTA,JEND
           DO I=1,IM
!             IF(abs(SM(I,J)-0.).lt.1.0E-5)THEN
             IF( (abs(SM(I,J)-0.)   .lt. 1.0E-5) .AND.
     &           (abs(SICE(I,J)-0.) .lt. 1.0E-5) ) THEN
              IF(CZMEAN(I,J).GT.1.E-6) THEN
               FACTRS=CZEN(I,J)/CZMEAN(I,J)
              ELSE
               FACTRS=0.0
              ENDIF
!              SOLAR=HBM2(I,J)*RSWIN(I,J)*FACTRS
              LLMH=NINT(LMH(I,J))
	      SOLAR=RSWIN(I,J)*FACTRS
              SFCTMP=T(I,J,LLMH)
              SFCQ=Q(I,J,LLMH)
              SFCPRS=PINT(I,J,LLMH+1)
!              IF(IVGTYP(I,J).EQ.0)PRINT*,'IVGTYP ZERO AT ',I,J
!     &        ,SM(I,J)
              IVG=IVGTYP(I,J)
!              IF(IVGTYP(I,J).EQ.0)IVG=7
!              CALL CANRES(SOLAR,SFCTMP,SFCQ,SFCPRS
!     &        ,SMC(I,J,1:NSOIL),GC(I,J),RC,IVG,ISLTYP(I,J))
!
              CALL CANRES(SOLAR,SFCTMP,SFCQ,SFCPRS
     &        ,SH2O(I,J,1:NSOIL),GC(I,J),RC,IVG,ISLTYP(I,J)
     &        ,RSMIN(I,J),NROOTS(I,J),SMCWLT(I,J),SMCREF(I,J)
     &        ,RCS(I,J),RCQ(I,J),RCT(I,J),RCSOIL(I,J),SLDPTH)
               IF(abs(SMCWLT(I,J)-0.5).lt.1.e-5)print*,
     &       'LARGE SMCWLT',i,j,SM(I,J),ISLTYP(I,J),SMCWLT(I,J)
             ELSE
              GC(I,J)=0.
              RSMIN(I,J)=0.
              NROOTS(I,J)=0
              SMCWLT(I,J)=0.
              SMCREF(I,J)=0.
              RCS(I,J)=0.
              RCQ(I,J)=0.
              RCT(I,J)=0.
              RCSOIL(I,J)=0.
             END IF
           ENDDO
         ENDDO
	 
         IF (IGET(220).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = GC(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(220),LVLS(1,IGET(220)),
     X        GRID1,IM,JM)                                                          
         ENDIF	 	     

         IF (IGET(234).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RSMIN(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(234),LVLS(1,IGET(234)),
     X        GRID1,IM,JM)                                                          
         ENDIF	
	 
         IF (IGET(235).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = FLOAT(NROOTS(I,J))
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(235),LVLS(1,IGET(235)),
     X        GRID1,IM,JM)                                                          
         ENDIF	

         IF (IGET(236).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = SMCWLT(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(236),LVLS(1,IGET(236)),
     X        GRID1,IM,JM)                                                          
         ENDIF	

         IF (IGET(237).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = SMCREF(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(237),LVLS(1,IGET(237)),
     X        GRID1,IM,JM)                                                          
         ENDIF	

         IF (IGET(238).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RCS(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(238),LVLS(1,IGET(238)),
     X        GRID1,IM,JM)                                                          
         ENDIF	

         IF (IGET(239).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RCT(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(239),LVLS(1,IGET(239)),
     X        GRID1,IM,JM)                                                          
         ENDIF	

         IF (IGET(240).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RCQ(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(240),LVLS(1,IGET(240)),
     X        GRID1,IM,JM)                                                          
         ENDIF	
         
         IF (IGET(241).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RCSOIL(I,J)
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(241),LVLS(1,IGET(241)),
     X        GRID1,IM,JM)                                                          
         ENDIF	
	 
	 print*,'outputting leaf area index= ',XLAI
         IF (IGET(254).GT.0 )THEN
          DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = XLAI
           ENDDO
          ENDDO
          ID(1:25) = 0
	  ID(02)= 130
          CALL GRIBIT(IGET(254),LVLS(1,IGET(254)),
     X        GRID1,IM,JM)                                                          
         ENDIF

      ENDIF
      END IF
C     
C     END OF ROUTINE
C     
C     
C       MODEL TOP REQUESTED BY CMAQ
      IF (IGET(282).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PT
            ENDDO
            ENDDO
            ID(1:25) = 0
            CALL GRIBIT(IGET(282),LVLS(1,IGET(282)),
     X           GRID1,IM,JM)
      ENDIF
C     
C       PRESSURE THICKNESS REQUESTED BY CMAQ
      IF (IGET(283).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PDTOP
            ENDDO
            ENDDO
            ID(1:25) = 0
	    II=IM/2
	    JJ=(JSTA+JEND)/2
	    DO L=1,LM
	     IF(PINT(II,JJ,L).GT.PDTOP)EXIT
	    END DO
	    PRINT*,'hybrid boundary ',L
	    ID(10)=1
	    ID(11)=L  
            CALL GRIBIT(IGET(283),LVLS(1,IGET(283)),
     X           GRID1,IM,JM)
      ENDIF
C      
C       SIGMA PRESSURE THICKNESS REQUESTED BY CMAQ
      IF (IGET(273).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
             GRID1(I,J)=PD(I,J)
            ENDDO
            ENDDO
            ID(1:25) = 0
	    ID(10)=L+1
	    ID(11)=LM
            CALL GRIBIT(IGET(273),LVLS(1,IGET(273)),
     X           GRID1,IM,JM)
      ENDIF
              
      RETURN
      END