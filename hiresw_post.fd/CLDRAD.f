      SUBROUTINE CLDRAD
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C SUBPROGRAM:    CLDRAD       POST SNDING/CLOUD/RADTN FIELDS
C   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-08-30       
C     
C ABSTRACT:  THIS ROUTINE COMPUTES/POSTS SOUNDING, CLOUD 
C   RELATED, AND RADIATION FIELDS.  UNDER THE HEADING OF 
C   SOUNDING FIELDS FALL THE THREE ETA MODEL LIFTED INDICES,
C   CAPE, CIN, AND TOTAL COLUMN PRECIPITABLE WATER.
C
C   THE THREE ETA MODEL LIFTED INDICES DIFFER ONLY IN THE
C   DEFINITION OF THE PARCEL TO LIFT.  ONE LIFTS PARCELS FROM
C   THE LOWEST ABOVE GROUND ETA LAYER.  ANOTHER LIFTS MEAN 
C   PARCELS FROM ANY OF NBND BOUNDARY LAYERS (SEE SUBROUTINE
C   BNDLYR).  THE FINAL TYPE OF LIFTED INDEX IS A BEST LIFTED
C   INDEX BASED ON THE NBND BOUNDARY LAYER LIFTED INDICES.
C
C   TWO TYPES OF CAPE/CIN ARE AVAILABLE.  ONE IS BASED ON PARCELS
C   IN THE LOWEST ETA LAYER ABOVE GROUND.  THE OTHER IS BASED 
C   ON A LAYER MEAN PARCEL IN THE N-TH BOUNDARY LAYER ABOVE 
C   THE GROUND.  SEE SUBROUTINE CALCAPE FOR DETAILS.
C
C   THE CLOUD FRACTION AND LIQUID CLOUD WATER FIELDS ARE DIRECTLY
C   FROM THE MODEL WITH MINIMAL POST PROCESSING.  THE LIQUID 
C   CLOUD WATER, 3-D CLOUD FRACTION, AND TEMPERATURE TENDENCIES
C   DUE TO PRECIPITATION ARE NOT POSTED IN THIS ROUTINE.  SEE
C   SUBROUTINE ETAFLD FOR THESE FIELDS.  LIFTING CONDENSATION
C   LEVEL HEIGHT AND PRESSURE ARE COMPUTED AND POSTED IN
C   SUBROUTINE MISCLN.  
C
C   THE RADIATION FIELDS POSTED BY THIS ROUTINE ARE THOSE COMPUTED
C   DIRECTLY IN THE MODEL.
C     
C PROGRAM HISTORY LOG:
C   93-08-30  RUSS TREADON
C   94-08-04  MICHAEL BALDWIN - ADDED OUTPUT OF INSTANTANEOUS SFC
C                               FLUXES OF NET SW AND LW DOWN RADIATION
C   97-04-25  MICHAEL BALDWIN - FIX PDS FOR PRECIPITABLE WATER
C   97-04-29  GEOFF MANIKIN - MOVED CLOUD TOP TEMPS CALCULATION
C                               TO THIS SUBROUTINE.  CHANGED METHOD
C                               OF DETERMINING WHERE CLOUD BASE AND
C                               TOP ARE FOUND AND ADDED HEIGHT OPTION
C                               FOR TOP AND BASE.
C   98-04-29  GEOFF MANIKIN - CHANGED VALUE FOR CLOUD BASE/TOP PRESSURES
C                               AND HEIGHTS FROM SPVAL TO -500
C   98-06-15  T BLACK       - CONVERSION FROM 1-D TO 2-D
C   98-07-17  MIKE BALDWIN  - REMOVED LABL84
C   00-01-04  JIM TUCCILLO  - MPI VERSION
C   00-02-22  GEOFF MANIKIN - CHANGED VALUE FOR CLOUD BASE/TOP PRESSURES
C                               AND HEIGHTS FROM SPVAL TO -500 (WAS NOT IN
C                               PREVIOUS IBM VERSION)
C   01-10-22  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
C   02-01-15  MIKE BALDWIN - WRF VERSION
C   05-01-06  H CHUANG - ADD VARIOUS CLOUD FIELDS
C   05-07-07  BINBIN ZHOU - ADD RSM MODEL
C   05-08-30  BINBIN ZHOU - ADD CEILING and FLIGHT CONDITION RESTRICTION
C
C     
C USAGE:    CALL CLDRAD
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
C       NONE
C     LIBRARY:
C       COMMON   - RQSTFLD
C                  CTLBLK
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : IBM SP
C$$$  
C
      use vrbls3d
      use vrbls2d
      use masks
C
C     INCLUDE GRID DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
C     
!      INCLUDE "parmeta"
!      INCLUDE "parmout"
      INCLUDE "params"
C
C     INCLUDE COMMON BLOCKS.
      INCLUDE "RQSTFLD.comm"
      INCLUDE "CTLBLK.comm"
C     
C     SET CELSIUS TO KELVIN CONVERSION.
      PARAMETER (C2K=273.15)
C     
C     DECLARE VARIABLES.
C     
      LOGICAL RUN,FIRST,RESTRT,SIGMA
      LOGICAL NEED(IM,JM)
      INTEGER L1D(IM,JM)
      INTEGER IBOTT(IM,JM),IBOTCu(IM,JM),IBOTDCu(IM,JM)
     & ,IBOTSCu(IM,JM),IBOTGr(IM,JM), ITOPT(IM,JM)
     &,ITOPCu(IM,JM),ITOPDCu(IM,JM),ITOPSCu(IM,JM)
     &,ITOPGr(IM,JM)
      REAL EGRID1(IM,JM),EGRID2(IM,JM),EGRID3(IM,JM)
      REAL GRID1(IM,JM),GRID2(IM,JM),CLDP(IM,JM),
     &        CLDZ(IM,JM),CLDT(IM,JM)
      
C     B ZHOU: For aviation:
      REAL  TCLD(IM,JM), CEILING(IM,JM), FLTCND(IM,JM)
C     
C
C*************************************************************************
C     START CLDRAD HERE.
C     
C***  BLOCK 1.  SOUNDING DERIVED FIELDS.
C     
C     ETA SURFACE TO 500MB LIFTED INDEX.  TO BE CONSISTENT WITH THE
C     LFM AND NGM POSTING WE ADD 273.15 TO THE LIFTED INDEX
C
C     THE BEST (SIX LAYER) AND BOUNDARY LAYER LIFTED INDICES ARE
C     COMPUTED AND POSTED IN SUBROUTINE MISCLN.
C
      IF (IGET(030).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           EGRID1(I,J) = SPVAL
         ENDDO
         ENDDO
C
         CALL OTLIFT(EGRID1)
C
         DO J=JSTA,JEND
         DO I=1,IM
           IF(EGRID1(I,J).LT.SPVAL) GRID1(I,J)=EGRID1(I,J) +TFRZ 
         ENDDO
         ENDDO
C
         ID(1:25)=0
         ID(10)  =50
         ID(11)  =100
         CALL GRIBIT(IGET(030),LVLS(1,IGET(030)),GRID1,IM,JM)
      ENDIF
C
C     SOUNDING DERIVED AREA INTEGRATED ENERGIES - CAPE AND CIN.
C       THIS IS THE SFC-BASED CAPE/CIN (lowest 70 mb searched)
      IF ((IGET(032).GT.0).OR.(IGET(107).GT.0)) THEN
         IF ( (LVLS(1,IGET(032)).GT.0) .OR. 
     X        (LVLS(1,IGET(107)).GT.0) ) THEN
            ITYPE = 1
	    DPBND=10.E2
            CALL CALCAPE(ITYPE,DPBND,P1D,T1D,Q1D,L1D,EGRID1,EGRID2,
     X           EGRID3)
C
C           CONVECTIVE AVAILABLE POTENTIAL ENERGY.
            IF (IGET(032).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = EGRID1(I,J)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               ID(1:25)=0
               CALL GRIBIT(IGET(032),LVLS(1,IGET(032)),GRID1,IM,JM)
            ENDIF
C
C           CONVECTIVE INHIBITION.
            IF (IGET(107).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = -1.*EGRID2(I,J)
               ENDDO
               ENDDO
               CALL BOUND(GRID1,D00,H99999)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = -1.*GRID1(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(107),LVLS(1,IGET(107)),GRID1,IM,JM)
            ENDIF
         ENDIF
      ENDIF
C
C     TOTAL COLUMN PRECIPITABLE WATER (SPECIFIC HUMIDITY).
      IF (IGET(080).GT.0) THEN
         CALL CALPW(GRID1,1)
         ID(1:25)=0
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(080),LVLS(1,IGET(080)),GRID1,IM,JM)
      ENDIF
C     
C     TOTAL COLUMN CLOUD WATER
      IF (IGET(200).GT.0) THEN
         CALL CALPW(GRID1,2)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(200),LVLS(1,IGET(200)),GRID1,IM,JM)
      ENDIF
C
C     TOTAL COLUMN CLOUD ICE
      IF (IGET(201).GT.0) THEN
         CALL CALPW(GRID1,3)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(201),LVLS(1,IGET(201)),GRID1,IM,JM)
      ENDIF
C
C     TOTAL COLUMN RAIN 
      IF (IGET(202).GT.0) THEN
         CALL CALPW(GRID1,4)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(202),LVLS(1,IGET(202)),GRID1,IM,JM)
      ENDIF
C
C     TOTAL COLUMN SNOW 
      IF (IGET(203).GT.0) THEN
         CALL CALPW(GRID1,5)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(203),LVLS(1,IGET(203)),GRID1,IM,JM)
      ENDIF
C
C     TOTAL COLUMN CONDENSATE 
      IF (IGET(204).GT.0) THEN
         CALL CALPW(GRID1,6)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(204),LVLS(1,IGET(204)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN SUPERCOOLED (<0C) LIQUID WATER 
      IF (IGET(285).GT.0) THEN
         CALL CALPW(GRID1,7)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(285),LVLS(1,IGET(285)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN MELTING (>0C) ICE
      IF (IGET(286).GT.0) THEN
         CALL CALPW(GRID1,8)
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL BOUND(GRID1,D00,H99999)
         CALL GRIBIT(IGET(286),LVLS(1,IGET(286)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN SHORT WAVE T TENDENCY
      IF (IGET(290).GT.0) THEN
         CALL CALPW(GRID1,9)
         ID(1:25)=0
         CALL GRIBIT(IGET(290),LVLS(1,IGET(290)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN LONG WAVE T TENDENCY
      IF (IGET(291).GT.0) THEN
         CALL CALPW(GRID1,10)
         ID(1:25)=0
         CALL GRIBIT(IGET(291),LVLS(1,IGET(291)),GRID1,IM,JM)
      ENDIF            
!
!     TOTAL COLUMN GRID SCALE LATENT HEATING (TIME AVE)
      IF (IGET(292).GT.0) THEN
         CALL CALPW(GRID1,11)
	 IF(AVRAIN.GT.0.)THEN
           RRNUM=1./AVRAIN
         ELSE
           RRNUM=0.
         ENDIF
!$omp  parallel do
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J)=GRID1(I,J)*RRNUM
         ENDDO
         ENDDO
         ID(1:25)=0
	 ITHEAT     = INT(THEAT)
         IF (ITHEAT .NE. 0) THEN
          IFINCR     = MOD(IFHR,ITHEAT)
         ELSE
          IFINCR=0
         END IF
         ID(19) = IFHR
         IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20) = 3
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITHEAT
         ELSE
          ID(18) = IFHR-IFINCR
         ENDIF
         IF(IFMIN .GE. 1)ID(18)=ID(18)*60
         IF (ID(18).LT.0) ID(18) = 0
         CALL GRIBIT(IGET(292),LVLS(1,IGET(292)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN CONVECTIVE LATENT HEATING (TIME AVE)
      IF (IGET(293).GT.0) THEN
         CALL CALPW(GRID1,12)
	 IF(AVRAIN.GT.0.)THEN
           RRNUM=1./AVCNVC
         ELSE
           RRNUM=0.
         ENDIF
!$omp  parallel do
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J)=GRID1(I,J)*RRNUM
         ENDDO
         ENDDO
         ID(1:25)=0
	 ITHEAT     = INT(THEAT)
         IF (ITHEAT .NE. 0) THEN
          IFINCR     = MOD(IFHR,ITHEAT)
         ELSE
          IFINCR=0
         END IF
         ID(19) = IFHR
         IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
         ID(20) = 3
         IF (IFINCR.EQ.0) THEN
          ID(18) = IFHR-ITHEAT
         ELSE
          ID(18) = IFHR-IFINCR
         ENDIF
         IF(IFMIN .GE. 1)ID(18)=ID(18)*60
         IF (ID(18).LT.0) ID(18) = 0
         CALL GRIBIT(IGET(293),LVLS(1,IGET(293)),GRID1,IM,JM)
      ENDIF
!
!     TOTAL COLUMN moisture convergence
      IF (IGET(295).GT.0) THEN
         CALL CALPW(GRID1,13)
         ID(1:25)=0
         CALL GRIBIT(IGET(295),LVLS(1,IGET(295)),GRID1,IM,JM)
      ENDIF
!
!     BOTTOM AND/OR TOP OF SUPERCOOLED (<0C) LIQUID WATER LAYER
      IF (IGET(287).GT.0 .OR. IGET(288).GT.0) THEN
         DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J)=-5000.
               GRID2(I,J)=-5000.
!-- Search for the base first, then look for the top if supercooled liquid exists
               LBOT=0
               LM=NINT(LMH(I,J))
               DO L=LM,1,-1
                  QCLD=QQW(I,J,L)+QQR(I,J,L)
                  IF (QCLD.GE.QCLDmin .AND. T(I,J,L).LT.TFRZ) THEN
                     LBOT=L
                     EXIT
                  ENDIF
               ENDDO    !--- End L loop
               IF (LBOT .GT. 0) THEN
  !-- Supercooled liquid exists, so get top & bottom heights.  In this case,
  !   be conservative and select the lower interface height at the bottom of the
  !   layer and the top interface height at the top of the layer.
                  GRID1(I,J)=ZINT(I,J,LBOT+1)
                  DO L=1,LM
                     QCLD=QQW(I,J,L)+QQR(I,J,L)
                     IF (QCLD.GE.QCLDmin .AND. T(I,J,L).LT.TFRZ) THEN
                        LTOP=L
                        EXIT
                     ENDIF
                  ENDDO    !--- End L loop
                  LTOP=MIN(LBOT,LTOP)
                  GRID2(I,J)=ZINT(I,J,LTOP)
               ENDIF    !--- End IF (LBOT .GT. 0)
            ENDDO       !--- End I loop
         ENDDO          !--- End J loop
         IF (IGET(287).GT.0) THEN
            ID(1:25)=0
            CALL GRIBIT(IGET(287),LVLS(1,IGET(287)),GRID1,IM,JM)
         ENDIF
         IF (IGET(288).GT.0) THEN
            DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J)=GRID2(I,J)
            ENDDO
            ENDDO
            ID(1:25)=0
            CALL GRIBIT(IGET(288),LVLS(1,IGET(288)),GRID1,IM,JM)
         ENDIF
      ENDIF
C
C
C     Convective cloud efficiency parameter used in convection ranges
C     from 0.2 (EFIMN in cuparm in model) to 1.0   (Ferrier, Feb '02) 
      IF (IGET(197).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           GRID1(I,J) = CLDEFI(I,J)
         ENDDO
         ENDDO
         ID(1:25)=0
         ID(02)=129      !--- Parameter Table 129, PDS Octet 4 = 129)
         CALL GRIBIT(IGET(197),LVLS(1,IGET(197)),GRID1,IM,JM)
      ENDIF
C   
C
C
C***  BLOCK 2.  2-D CLOUD FIELDS.
C
C     LOW CLOUD FRACTION.
      IF (IGET(037).GT.0) THEN	  
        DO J=JSTA,JEND
        DO I=1,IM
          GRID1(I,J) = CFRACL(I,J)*100.
        ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(037),LVLS(1,IGET(037)),GRID1,IM,JM)
      ENDIF
C     
C     MIDDLE CLOUD FRACTION.
      IF (IGET(038).GT.0) THEN
        DO J=JSTA,JEND
        DO I=1,IM
           GRID1(I,J) = CFRACM(I,J)*100.
        ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(038),LVLS(1,IGET(038)),GRID1,IM,JM)
      ENDIF
C     
C     HIGH CLOUD FRACTION.
      IF (IGET(039).GT.0) THEN
        DO J=JSTA,JEND
        DO I=1,IM
           GRID1(I,J) = CFRACH(I,J)*100.
        ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(039),LVLS(1,IGET(039)),GRID1,IM,JM)
      ENDIF
C     
C     TOTAL CLOUD FRACTION (INSTANTANEOUS).
      IF ((IGET(161).GT.0) .OR. (IGET(260).GT.0)) THEN
         IF(MODELNAME .EQ. 'NCAR')THEN
          DO J=JSTA,JEND
          DO I=1,IM
            EGRID1(I,J)=CLDFRA(I,J)
          ENDDO
          ENDDO
         ELSE IF (MODELNAME.EQ.'NMM'.OR.MODELNAME.EQ.'RSM')THEN
          DO J=JSTA,JEND
          DO I=1,IM
!           EGRID1(I,J)=AMAX1(CFRACL(I,J),
!     1                 AMAX1(CFRACM(I,J),CFRACH(I,J)))
            EGRID1(I,J)=1.-(1.-CFRACL(I,J))*(1.-CFRACM(I,J))*
     &                 (1.-CFRACH(I,J))
          ENDDO
          ENDDO
         END IF
         DO J=JSTA,JEND
         DO I=1,IM
            GRID1(I,J) = EGRID1(I,J)*100.
	    TCLD(I,J)  = EGRID1(I,J)*100.         !B ZHOU, PASSED to CALCEILING
         ENDDO
         ENDDO
         IF (IGET(161).GT.0) THEN
            ID(1:25)=0
            CALL GRIBIT(IGET(161),LVLS(1,IGET(161)),GRID1,IM,JM)
         ENDIF
      ENDIF
C
C     TIME AVERAGED TOTAL CLOUD FRACTION.
         IF (IGET(144).GT.0) THEN
           IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	   ELSE IF (MODELNAME .EQ. 'NMM')THEN
            DO J=JSTA,JEND
            DO I=1,IM
!               RSUM = NCFRST(I,J)+NCFRCV(I,J)
!               IF (RSUM.GT.0.0) THEN
!                  EGRID1(I,J)=(ACFRST(I,J)+ACFRCV(I,J))/RSUM
!               ELSE
!                  EGRID1(I,J) = D00
!               ENDIF
!ADDED BRAD'S MODIFICATION
               RSUM = D00
               IF (NCFRST(I,J) .GT. 0) RSUM=ACFRST(I,J)/NCFRST(I,J)
               IF (NCFRCV(I,J) .GT. 0) 
     &            RSUM=MAX(RSUM, ACFRCV(I,J)/NCFRCV(I,J))
               EGRID1(I,J) = RSUM
            ENDDO
            ENDDO
C
            DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J) = EGRID1(I,J)*100.
            ENDDO
            ENDDO
	   END IF 
          IF(MODELNAME.EQ.'NMM')THEN
           ID(1:25)= 0
           ITCLOD     = INT(TCLOD)
           IF(ITCLOD .ne. 0) then
            IFINCR     = MOD(IFHR,ITCLOD)
            IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
           ELSE
            IFINCR     = 0
           endif

           ID(19)  = IFHR
	   IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN  !USE MIN FOR OFF-HR FORECAST
           ID(20)  = 3
           IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITCLOD
           ELSE
               ID(18)  = IFHR-IFINCR
               IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
           ENDIF
           IF (ID(18).LT.0) ID(18) = 0
          ENDIF
           CALL GRIBIT(IGET(144),LVLS(1,IGET(144)),GRID1,IM,JM)
         ENDIF
C
C     TIME AVERAGED STRATIFORM CLOUD FRACTION.
         IF (IGET(139).GT.0) THEN
           IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	   ELSE IF (MODELNAME .EQ. 'NMM')THEN	 
            DO J=JSTA,JEND
            DO I=1,IM
               IF (NCFRST(I,J).GT.0.0) THEN
                  EGRID1(I,J) = ACFRST(I,J)/NCFRST(I,J)
               ELSE
                  EGRID1(I,J) = D00
               ENDIF
            ENDDO
            ENDDO
C
            DO J=JSTA,JEND
            DO I=1,IM
              GRID1(I,J) = EGRID1(I,J)*100.
            ENDDO
            ENDDO
	   END IF 
          IF(MODELNAME.EQ.'NMM')THEN
           ID(1:25)=0
           ITCLOD     = INT(TCLOD)
	   IF(ITCLOD .ne. 0) then
            IFINCR     = MOD(IFHR,ITCLOD)
	    IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	   ELSE
	    IFINCR     = 0
           endif 
           ID(19)  = IFHR
	   IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
           ID(20)  = 3
           IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITCLOD
           ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
           ENDIF
           IF (ID(18).LT.0) ID(18) = 0
          ENDIF
           CALL GRIBIT(IGET(139),LVLS(1,IGET(139)),GRID1,IM,JM)
         ENDIF
C    
C     TIME AVERAGED CONVECTIVE CLOUD FRACTION.
         IF (IGET(143).GT.0) THEN
           IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	   ELSE IF (MODELNAME .EQ. 'NMM')THEN	 
            DO J=JSTA,JEND
            DO I=1,IM
               IF (NCFRCV(I,J).GT.0.0) THEN
                  EGRID1(I,J) = ACFRCV(I,J)/NCFRCV(I,J)
               ELSE
                  EGRID1(I,J) = D00
               ENDIF
            ENDDO
            ENDDO
C
            DO J=JSTA,JEND
            DO I=1,IM
               GRID1(I,J) = EGRID1(I,J)*100.
            ENDDO
            ENDDO
	   END IF
          IF(MODELNAME.EQ.'NMM')THEN 
           ID(1:25)=0
           ITCLOD     = INT(TCLOD)
	   IF(ITCLOD .ne. 0) then
            IFINCR     = MOD(IFHR,ITCLOD)
	    IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITCLOD*60)
	   ELSE
	    IFINCR     = 0
           endif 
           ID(19)  = IFHR
	   IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
           ID(20)  = 3
           IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITCLOD
           ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
           ENDIF
           IF (ID(18).LT.0) ID(18) = 0
          ENDIF
           CALL GRIBIT(IGET(143),LVLS(1,IGET(143)),GRID1,IM,JM)
         ENDIF
C    
C     CLOUD BASE AND TOP FIELDS 
      IF((IGET(148).GT.0) .OR. (IGET(149).GT.0) .OR.
     &    (IGET(168).GT.0) .OR. (IGET(178).GT.0) .OR.
     &    (IGET(179).GT.0) .OR. (IGET(194).GT.0) .OR.
     &    (IGET(195).GT.0) .OR. (IGET(260).GT.0) .OR. 
     &    (IGET(275).GT.0))  THEN
  !
  !--- Calculate grid-scale cloud base & top arrays (Ferrier, Feb '02)
  !
  !--- Rain is not part of cloud, only cloud water + cloud ice + snow
  !
        DO J=JSTA,JEND
          DO I=1,IM
    !
    !--- Various convective cloud base & cloud top levels
    !
            IBOTCu(I,J)=NINT(HBOT(I,J))
            IBOTDCu(I,J)=NINT(HBOTD(I,J))
            IBOTSCu(I,J)=NINT(HBOTS(I,J))
            ITOPCu(I,J)=NINT(HTOP(I,J))
            ITOPDCu(I,J)=NINT(HTOPD(I,J))
            ITOPSCu(I,J)=NINT(HTOPS(I,J))
            IF (IBOTCu(I,J)-ITOPCu(I,J) .LE. 1) THEN
              IBOTCu(I,J)=0
              ITOPCu(I,J)=100
            ENDIF
            IF (IBOTDCu(I,J)-ITOPDCu(I,J) .LE. 1) THEN
              IBOTDCu(I,J)=0
              ITOPDCu(I,J)=100
            ENDIF
            IF (IBOTSCu(I,J)-ITOPSCu(I,J) .LE. 1) THEN
              IBOTSCu(I,J)=0
              ITOPSCu(I,J)=100
            ENDIF
    !
    !--- Grid-scale cloud base & cloud top levels 
    !
    !--- Grid-scale cloud occurs when the mixing ratio exceeds QCLDmin
    !
            IBOTGr(I,J)=0
            DO L=NINT(LMH(I,J)),1,-1
              QCLD=QQW(I,J,L)+QQI(I,J,L)+QQS(I,J,L)
              IF (QCLD .GE. QCLDmin) THEN
                IBOTGr(I,J)=L
                EXIT
              ENDIF
            ENDDO    !--- End L loop
            ITOPGr(I,J)=100
            DO L=1,NINT(LMH(I,J))
              QCLD=QQW(I,J,L)+QQI(I,J,L)+QQS(I,J,L)
              IF (QCLD .GE. QCLDmin) THEN
                ITOPGr(I,J)=L
                EXIT
              ENDIF
            ENDDO    !--- End L loop
    !
    !--- Combined (convective & grid-scale) cloud base & cloud top levels 
            IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	     IBOTT(I,J)=IBOTGr(I,J)
	     ITOPT(I,J)=ITOPGr(I,J)
	    ELSE
             IBOTT(I,J)=MAX(IBOTGr(I,J), IBOTCu(I,J))
             ITOPT(I,J)=MIN(ITOPGr(I,J), ITOPCu(I,J))
	    END IF 
          ENDDO      !--- End I loop
        ENDDO        !--- End J loop
      ENDIF          !--- End IF tests 
!
!-------------------------------------------------
!-----------  VARIOUS CLOUD BASE FIELDS ----------
!-------------------------------------------------
!
!--- "TOTAL" CLOUD BASE FIELDS (convective + grid-scale;  Ferrier, Feb '02)
!
      IF ((IGET(148).GT.0) .OR. (IGET(178).GT.0)
     1     .OR. (IGET(260).GT.0) ) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            IBOT=IBOTT(I,J)
            IF (IBOT .LE. 0) THEN
              CLDP(I,J) = -50000.
              CLDZ(I,J) = -5000.
            ELSE IF (IBOT .LE. NINT(LMH(I,J))) THEN
              CLDP(I,J) = PMID(I,J,IBOT)
              IF (IBOT .EQ. LM) THEN
                CLDZ(I,J) = ZINT(I,J,LM)
              ELSE
                CLDZ(I,J) = HTM(I,J,IBOT+1)*T(I,J,IBOT+1)
     1                     *(Q(I,J,IBOT+1)*D608+H1)*ROG*
     2                     (LOG(PINT(I,J,IBOT+1))-LOG(CLDP(I,J)))
     3                     +ZINT(I,J,IBOT+1)
              ENDIF     !--- End IF (IBOT .EQ. LM) ...
            ENDIF       !--- End IF (IBOT .LE. 0) ...
          ENDDO         !--- End DO I loop
        ENDDO           !--- End DO J loop
C   CLOUD BOTTOM PRESSURE
         IF (IGET(148).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = CLDP(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(148),LVLS(1,IGET(148)),GRID1,IM,JM)
         ENDIF 
C    CLOUD BOTTOM HEIGHT
         IF (IGET(178).GT.0) THEN
  !--- Parameter was set to 148 in operational code  (Ferrier, Feb '02)
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = CLDZ(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(178),LVLS(1,IGET(178)),GRID1,IM,JM)
         ENDIF
      ENDIF
      
C    B. ZHOU: CEILING
        IF (IGET(260).GT.0) THEN                                                                                                          
            CALL CALCEILING(CLDZ,TCLD,CEILING)                                                                                   
            DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = CEILING(I,J)
             ENDDO
            ENDDO
            ID(1:25)=0
            CALL GRIBIT(IGET(260),LVLS(1,IGET(260)),GRID1,IM,JM)
         ENDIF
                                                                                                          
C    B. ZHOU: FLIGHT CONDITION RESTRICTION
        IF (IGET(261).GT.0) THEN
            CALL CALFLTCND(CEILING,FLTCND)
            DO J=JSTA,JEND
             DO I=1,IM
               GRID1(I,J) = FLTCND(I,J)
             ENDDO
            ENDDO
            ID(1:25)=0
            CALL GRIBIT(IGET(261),LVLS(1,IGET(261)),GRID1,IM,JM)
         ENDIF
!
!---  Convective cloud base pressures (deep & shallow; Ferrier, Feb '02)
!
      IF (IGET(188) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            IBOT=IBOTCu(I,J)
            IF (IBOT.GT.0 .AND. IBOT.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,IBOT)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(188),LVLS(1,IGET(188)),GRID1,IM,JM)
       ENDIF
C
!---  Deep convective cloud base pressures  (Ferrier, Feb '02)
!
      IF (IGET(192) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            IBOT=IBOTDCu(I,J)
            IF (IBOT.GT.0 .AND. IBOT.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,IBOT)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(192),LVLS(1,IGET(192)),GRID1,IM,JM)
       ENDIF 
!---  Shallow convective cloud base pressures   (Ferrier, Feb '02)
!
      IF (IGET(190) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            IBOT=IBOTSCu(I,J)  
            IF (IBOT.GT.0 .AND. IBOT.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,IBOT)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(190),LVLS(1,IGET(190)),GRID1,IM,JM)
       ENDIF
  !---  Base of grid-scale cloudiness   (Ferrier, Feb '02)
  !
      IF (IGET(194) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            IBOT=IBOTGr(I,J)
            IF (IBOT.GT.0 .AND. IBOT.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,IBOT)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(194),LVLS(1,IGET(194)),GRID1,IM,JM)
       ENDIF
!
!------------------------------------------------
!-----------  VARIOUS CLOUD TOP FIELDS ----------
!------------------------------------------------
!
!--- "TOTAL" CLOUD TOP FIELDS (convective + grid-scale;  Ferrier, Feb '02)
!
      IF ((IGET(149).GT.0) .OR. (IGET(179).GT.0) .OR.
     X    (IGET(168).GT.0) .OR. (IGET(275).GT.0)) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            ITOP=ITOPT(I,J)
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              CLDP(I,J) = PMID(I,J,ITOP)
              CLDT(I,J) = T(I,J,ITOP)
              IF (ITOP .EQ. LM) THEN
                CLDZ(I,J) = ZINT(I,J,LM)
              ELSE
                CLDZ(I,J) = HTM(I,J,ITOP+1)*T(I,J,ITOP+1)
     1                    *(Q(I,J,ITOP+1)*D608+H1)*ROG*
     2                     (LOG(PINT(I,J,ITOP+1))-LOG(CLDP(I,J)))
     3                    +ZINT(I,J,ITOP+1)
              ENDIF    !--- End IF (ITOP .EQ. LM) ...
            ELSE
              CLDP(I,J) = -50000.
              CLDZ(I,J) = -5000.
              CLDT(I,J) = -500.
            ENDIF      !--- End IF (ITOP.GT.0 .AND. ITOP.LE.LMH(I,J)) ...
          ENDDO        !--- End DO I loop
        ENDDO          !--- End DO J loop
C
C   CLOUD TOP PRESSURE
C
         IF (IGET(149).GT.0) THEN
              DO J=JSTA,JEND
              DO I=1,IM
                 GRID1(I,J) = CLDP(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(149),LVLS(1,IGET(149)),GRID1,IM,JM)
         ENDIF
C   CLOUD TOP HEIGHT
C
          IF (IGET(179).GT.0) THEN
              DO J=JSTA,JEND
              DO I=1,IM
                 GRID1(I,J) = CLDZ(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(179),LVLS(1,IGET(179)),GRID1,IM,JM)
         ENDIF
C
C   CLOUD TOP TEMPS
C
          IF (IGET(168).GT.0) THEN 
              DO J=JSTA,JEND
              DO I=1,IM
                 GRID1(I,J) = CLDT(I,J)
               ENDDO
               ENDDO
               ID(1:25)=0
               CALL GRIBIT(IGET(168),LVLS(1,IGET(168)),GRID1,IM,JM)
         ENDIF 
      ENDIF
C
Chuang  CLOUD TOP BRIGHTNESS TEMPERATURE
          IF (IGET(275).GT.0) THEN
             num_thick=0   ! for debug
             DO J=JSTA,JEND
             DO I=1,IM
               opdepth=0.
               llmh=nint(lmh(i,j))
               do k=1,llmh
                 dp=pint(i,j,k+1)-pint(i,j,k)
                 opdepth=opdepth+( abscoef*qqw(i,j,k)+
     &                   abscoefi*( qqi(i,j,k)+qqs(i,j,k) ) )*dp
                 if (opdepth > 1.) exit
               enddo
               if (opdepth > 1.) num_thick=num_thick+1   ! for debug
               k=min(k,llmh)
	     GRID1(I,J)=T(i,j,k)
             ENDDO
             ENDDO
      print *,'num_points, num_thick = ',(jend-jsta+1)*im,num_thick
!!              k=0
!! 20           opdepthu=opdepthd
!!              k=k+1
!!!              if(k.eq.1) then
!!!               dp=pint(i,j,itop+k)-pmid(i,j,itop)
!!!               opdepthd=opdepthu+(abscoef*(0.75*qqw(i,j,itop)+
!!!     &                  0.25*qqw(i,j,itop+1))+abscoefi*
!!!     &                  (0.75*qqi(i,j,itop)+0.25*qqi(i,j,itop+1)))
!!!     &                        *dp/g
!!!              else
!!               dp=pint(i,j,k+1)-pint(i,j,k)
!!               opdepthd=opdepthu+(abscoef*qqw(i,j,k)+
!!     &                        abscoefi*qqi(i,j,k))*dp
!!!              end if
!!	      
!!              lmhh=nint(lmh(i,j))
!!              if (opdepthd.lt.1..and. k.lt.lmhh) then
!!               goto 20
!!              elseif (opdepthd.lt.1..and. k.eq.lmhh) then
!!	       GRID1(I,J)=T(i,j,lmhh )
!!!               prsctt=pmid(i,j,lmhh)
!!              else
!!!	       GRID1(I,J)=T(i,j,k) 
!!               if(k.eq.1)then
!!	         GRID1(I,J)=T(i,j,k)
!!	       else if(k.eq.lmhh)then
!!	         GRID1(I,J)=T(i,j,k)
!!	       else 	 	 
!!                 fac=(1.-opdepthu)/(opdepthd-opdepthu)
!!	         GRID1(I,J)=(T(i,j,k)+T(i,j,k-1))/2.0+
!!     &             (T(i,j,k+1)-T(i,j,k-1))/2.0*fac 
!!               end if    	       
!!!               prsctt=pf(i,j,k-1)+fac*(pf(i,j,k)-pf(i,j,k-1))
!!!               prsctt=min(prs(i,j,mkzh),max(prs(i,j,1),prsctt))
!!              endif
!!!              do 30 k=2,mkzh
!!!              if (prsctt.ge.prs(i,j,k-1).and.prsctt.le.prs(i,j,k)) then
!!!               fac=(prsctt-prs(i,j,k-1))/(prs(i,j,k)-prs(i,j,k-1))
!!!               ctt(i,j)=tmk(i,j,k-1)+
!!!     &            fac*(tmk(i,j,k)-tmk(i,j,k-1))-celkel
!!!               goto 40
!!!              endif
!!!   30       continue
!!!   40       continue 
!!             END DO
!!	     END DO 
            ID(1:25)=0
!	    ID(02)=129    ! Parameter Table 129
            CALL GRIBIT(IGET(275),LVLS(1,IGET(275)),GRID1,IM,JM)
         ENDIF

!
!---  Convective cloud top pressures (deep & shallow; Ferrier, Feb '02)
!
      IF (IGET(189) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            ITOP=ITOPCu(I,J) 
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,ITOP)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(189),LVLS(1,IGET(189)),GRID1,IM,JM)
      END IF
!
!---  Deep convective cloud top pressures   (Ferrier, Feb '02)
!
      IF (IGET(193) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            ITOP=ITOPDCu(I,J)
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,ITOP)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(193),LVLS(1,IGET(193)),GRID1,IM,JM) 
      END IF
!---  Shallow convective cloud top pressures  (Ferrier, Feb '02)
!
      IF (IGET(191) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            ITOP=ITOPSCu(I,J)
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,ITOP)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(191),LVLS(1,IGET(191)),GRID1,IM,JM)
      END IF
!
!---  Top of grid-scale cloudiness  (Ferrier, Feb '02)
!
      IF (IGET(195) .GT. 0) THEN
        DO J=JSTA,JEND
          DO I=1,IM
            ITOP=ITOPGr(I,J)
            IF (ITOP.GT.0 .AND. ITOP.LE.NINT(LMH(I,J))) THEN
              GRID1(I,J) = PMID(I,J,ITOP)
            ELSE
              GRID1(I,J) = -50000.
            ENDIF
          ENDDO
        ENDDO
        ID(1:25)=0
        CALL GRIBIT(IGET(195),LVLS(1,IGET(195)),GRID1,IM,JM)
      END IF
!
!--- Convective cloud fractions from modified Slingo (1987)
!
      IF (IGET(196) .GT. 0) THEN
          DO J=JSTA,JEND
          DO I=1,IM
            GRID1(I,J)=CNVCFR(I,J)
          ENDDO
          ENDDO
          ID(1:25)=0
           CALL GRIBIT(IGET(196),LVLS(1,IGET(196)),GRID1,IM,JM)
      END IF
C***  BLOCK 3.  RADIATION FIELDS.
C     
C
C     TIME AVERAGED SURFACE SHORT WAVE INCOMING RADIATION.
         IF (IGET(126).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE  
c          print*,'ARDSW in CLDRAD=',ARDSW 
           IF(ARDSW.GT.0.)THEN
             RRNUM=1./ARDSW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
              GRID1(I,J) = ASWIN(I,J)*RRNUM
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDSW     = INT(TRDSW)
	    IF(ITRDSW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDSW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	    ELSE
	     IFINCR     = 0
            endif 	    
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDSW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF 
          CALL GRIBIT(IGET(126),LVLS(1,IGET(126)),GRID1,IM,JM)
         ENDIF
C
C     TIME AVERAGED SURFACE LONG WAVE INCOMING RADIATION.
         IF (IGET(127).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
           IF(ARDLW.GT.0.)THEN
             RRNUM=1./ARDLW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = ALWIN(I,J)*RRNUM
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDLW     = INT(TRDLW)
	    IF(ITRDLW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDLW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDLW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDLW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(127),LVLS(1,IGET(127)),GRID1,IM,JM)
         ENDIF
C
C     TIME AVERAGED SURFACE SHORT WAVE OUTGOING RADIATION.
         IF (IGET(128).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
           IF(ARDSW.GT.0.)THEN
             RRNUM=1./ARDSW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = -1.0*ASWOUT(I,J)*RRNUM
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDSW     = INT(TRDSW)
	    IF(ITRDSW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDSW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDSW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(128),LVLS(1,IGET(128)),GRID1,IM,JM)
         ENDIF
C
C     TIME AVERAGED SURFACE LONG WAVE OUTGOING RADIATION.
         IF (IGET(129).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
           IF(ARDLW.GT.0.)THEN
             RRNUM=1./ARDLW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = -1.0*ALWOUT(I,J)*RRNUM
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDLW     = INT(TRDLW)
	    IF(ITRDLW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDLW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDLW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDLW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(129),LVLS(1,IGET(129)),GRID1,IM,JM)
         ENDIF
C
C     TIME AVERAGED TOP OF THE ATMOSPHERE SHORT WAVE RADIATION.
         IF (IGET(130).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
           IF(ARDSW.GT.0.)THEN
             RRNUM=1./ARDSW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = ASWTOA(I,J)*RRNUM
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDSW     = INT(TRDSW)
	    IF(ITRDSW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDSW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDSW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDSW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(130),LVLS(1,IGET(130)),GRID1,IM,JM)
         ENDIF
C
C     TIME AVERAGED TOP OF THE ATMOSPHERE LONG WAVE RADIATION.
         IF (IGET(131).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	    GRID1=SPVAL
	    ID(1:25)=0
	  ELSE
           IF(ARDLW.GT.0.)THEN
             RRNUM=1./ARDLW
           ELSE
             RRNUM=0.
           ENDIF
           DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = ALWTOA(I,J)*RRNUM
           ENDDO
           ENDDO
            ID(1:25)=0
            ITRDLW     = INT(TRDLW)
            IF(ITRDLW .ne. 0) then
             IFINCR     = MOD(IFHR,ITRDLW)
	     IF(IFMIN .GE. 1)IFINCR= MOD(IFHR*60+IFMIN,ITRDLW*60)
	    ELSE
	     IFINCR     = 0
            endif
            ID(19)  = IFHR
	    IF(IFMIN .GE. 1)ID(19)=IFHR*60+IFMIN
            ID(20)  = 3
            IF (IFINCR.EQ.0) THEN
               ID(18)  = IFHR-ITRDLW
            ELSE
               ID(18)  = IFHR-IFINCR
	       IF(IFMIN .GE. 1)ID(18)=IFHR*60+IFMIN-IFINCR
            ENDIF
            IF (ID(18).LT.0) ID(18) = 0
	  END IF  
          CALL GRIBIT(IGET(131),LVLS(1,IGET(131)),GRID1,IM,JM)
         ENDIF
C
C     CURRENT TOP OF THE ATMOSPHERE LONG WAVE RADIATION.
         IF (IGET(274).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	   GRID1=SPVAL
	   ID(1:25)=0
	  ELSE
           DO J=JSTA,JEND
           DO I=1,IM
             GRID1(I,J) = RLWTOA(I,J)
           ENDDO
           ENDDO
           ID(1:25)=0
	  END IF  
          CALL GRIBIT(IGET(274),LVLS(1,IGET(274)),GRID1,IM,JM)
         ENDIF
C
C     CLOUD TOP BRIGHTNESS TEMPERATURE FROM TOA OUTGOING LW.
         IF (IGET(265).GT.0) THEN
	  IF(MODELNAME .EQ. 'NCAR'.OR.MODELNAME.EQ.'RSM')THEN
	   GRID1=SPVAL
	  ELSE
           DO J=JSTA,JEND
           DO I=1,IM
             IF(RLWTOA(I,J) .LT. SPVAL)
     +         GRID1(I,J) = (RLWTOA(I,J)*STBOL)**0.25
           ENDDO
           ENDDO
	  END IF  
	  ID(1:25)=0
	  ID(02)=129    ! Parameter Table 129
          CALL GRIBIT(IGET(265),LVLS(1,IGET(265)),GRID1,IM,JM)
         ENDIF
C     
C     CURRENT INCOMING SW RADIATION AT THE SURFACE.
      IF (IGET(156).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           IF(CZMEAN(I,J).GT.1.E-6) THEN
             FACTRS=CZEN(I,J)/CZMEAN(I,J)
           ELSE
             FACTRS=0.0
           ENDIF
           GRID1(I,J)=RSWIN(I,J)*FACTRS
         ENDDO
         ENDDO
C
         ID(1:25)=0
         CALL GRIBIT(IGET(156),LVLS(1,IGET(156)),GRID1,IM,JM)
      ENDIF
C     
C     CURRENT INCOMING LW RADIATION AT THE SURFACE.
      IF (IGET(157).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
          IF(MODELNAME.eq."RSM") THEN      !add by Binbin: RSM has direct RLWIN output
           GRID1(I,J)=RLWIN(I,J)
          ELSE
           IF(SIGT4(I,J).GT.0.0) THEN
             LLMH=NINT(LMH(I,J))
             TLMH=T(I,J,LLMH)
             FACTRL=5.67E-8*TLMH*TLMH*TLMH*TLMH/SIGT4(I,J)
           ELSE
             FACTRL=0.0
           ENDIF
           GRID1(I,J)=RLWIN(I,J)*FACTRL
          ENDIF
         ENDDO
         ENDDO
C
         ID(1:25)=0
         CALL GRIBIT(IGET(157),LVLS(1,IGET(157)),GRID1,IM,JM)
      ENDIF
C     
C     CURRENT OUTGOING SW RADIATION AT THE SURFACE.
      IF (IGET(141).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
           IF(CZMEAN(I,J).GT.1.E-6) THEN
             FACTRS=CZEN(I,J)/CZMEAN(I,J)
           ELSE
             FACTRS=0.0
           ENDIF
           GRID1(I,J)=RSWOUT(I,J)*FACTRS
         ENDDO
         ENDDO
C
         ID(1:25)=0
         CALL GRIBIT(IGET(141),LVLS(1,IGET(141)),GRID1,IM,JM)
      ENDIF
C     
C     CURRENT OUTGOING LW RADIATION AT THE SURFACE.
      IF (IGET(142).GT.0) THEN
               DO J=JSTA,JEND
               DO I=1,IM
                 GRID1(I,J) = RADOT(I,J)
               ENDDO
               ENDDO
         ID(1:25)=0
         CALL GRIBIT(IGET(142),LVLS(1,IGET(142)),GRID1,IM,JM)
      ENDIF
C     
C     CURRENT (instantaneous) INCOMING CLEARSKY SW RADIATION AT THE SURFACE.
      IF (IGET(262).GT.0) THEN
         DO J=JSTA,JEND
         DO I=1,IM
	   IF(CZMEAN(I,J).GT.1.E-6) THEN
             FACTRS=CZEN(I,J)/CZMEAN(I,J)
           ELSE
             FACTRS=0.0
           ENDIF
           GRID1(I,J) = RSWINC(I,J)*FACTRS
         ENDDO
	 ENDDO
         ID(1:25)=0
         CALL GRIBIT(IGET(262),LVLS(1,IGET(262)),GRID1,IM,JM)
      ENDIF

C
C     END OF ROUTINE.
C
      RETURN
      END