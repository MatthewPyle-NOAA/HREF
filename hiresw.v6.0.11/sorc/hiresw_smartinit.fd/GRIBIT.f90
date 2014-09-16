      SUBROUTINE GRIBIT(ID,RITEHD,GRID,GDIN,LUNOUT,DECI)
      use grddef
      use constants
    INTERFACE
      SUBROUTINE GET_BITS(IBM,SGDS,LEN,MG,G,ISCALE,GROUND,GMIN,GMAX,NBIT)
      DIMENSION MG(LEN),G(LEN),GROUND(LEN)
      END SUBROUTINE GET_BITS
    END INTERFACE

!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    GRIBIT      POST FIELDS IN GRIB1
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 93-06-18       
!     
! ABSTRACT:
!     THIS ROUTINE POSTS THE DATA IN THE PASSED ARRAY GRID
!     TO THE OUTPUT FILE IN GRIB1 FORMAT.
!     
! PROGRAM HISTORY LOG:
!   93-06-18  RUSS TREADON
!   93-11-23  RUSS TREADON - REMOVED CODE GENERATING GRIB INDEX FILE.
!   98-07-17  MIKE BALDWIN - REMOVED LABL84, NOW USING ID
!   02-06-17  MIKE BALDWIN - WRF VERSION
!   05-12-05  H CHUANG - ADD CAPABILITY TO OUTPUT OFF-HOUR FORECAST WHICH HAS
!               NO INPACTS ON ON-HOUR FORECAST
!  12-11-25   J. MCQUEEN  - F90 version
!     
! USAGE:    CALL GRIBIT(IFLD,ILVL,GRID,IMOUT,JMOUT)
!   INPUT ARGUMENT LIST:
!     IFLD     - FIELD ID TAG.
!     ILVL     - INTEGER TAG FOR LEVEL OF FIELD.
!     GRID     - FIELD TO BE POSTED IN GRIB.
!     IMOUT    - FIRST DIMENSION OF OUTPUT GRID.
!     JMOUT    - SECOND DIMENSION OF OUTPUT GRID.

!   OUTPUT ARGUMENT LIST: 
!     
!   OUTPUT FILES:
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!     GETENV   - SUBROUTINE TO GET VALUE OF ENVIRONMENT VARIABLE.
!     MINMAX   - DETERMINES MIN/MAX VALUES IN AN ARRAY.
!     WRYTE    - WRITE DATA OUT BY BYTES.
!     GET_BITS   - COMPUTE NUMBER OF BITS 
!     VARIOUS W3LIB ROUTINES
!     LIBRARY:
!       COMMON   - CTLBLK
!                  RQSTFLD
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN90
!     MACHINE : WCOSS    
!$$$  
!     
!     INCLUDE GRID DIMENSIONS.  SET/DERIVE PARAMETERS.

!     
!     GRIB1 PARAMETERS.
!        MNBIT  = MINIMUM NUMBER OF BITS TO USE IN PACKING.
!        MXBIT  = MAXIMUM NUMBER OF BITS TO USE IN PACKING.
!        LENPDS = LENGTH OF GRIB1 PDS.
!        LENGDS = LENGTH OF GRIB1 GDS.
!     
      PARAMETER (MXBIT=16,LENPDS=28,LENGDS=32)
      PARAMETER (SMALL=1.E-6)
      TYPE (GINFO),INTENT(IN) :: GDIN

      LOGICAL RITEHD
      CHARACTER*1  CDOT,CCOL
      CHARACTER*4  RESTHR
      CHARACTER  PROJ*6,DATSET*8
      CHARACTER*10 DESCR2,DESCR3
      CHARACTER*28 PDS
      CHARACTER*50 ENVAR
      CHARACTER*80 FNAME
      INTEGER IBDSFL(9),IGDS(18),ID(25)
      INTEGER ICENT,IYY,IMM,IDD,IHRST,DATE,OGRD

      CHARACTER*1, ALLOCATABLE    :: KBUF(:)
      REAL,        ALLOCATABLE    :: HOLDGRID(:,:)
      REAL,        INTENT(INOUT)  :: GRID(:,:)
      INTEGER,     ALLOCATABLE    :: IBMAP(:,:),IGRD(:,:)

      LOGICAL  NEWFILE
      INTEGER IH(5),LUNOUT
!     
!     SET DEFAULT GRIB1 PARAMETERS.  
!     PARAMETERS MNBIT, MXBIT, IBX, AND NBIT ARE USED 
!     IN THE CALL TO GET_BITS.
!        IBX    = DESIRED BINARY PRECISION.
!        NBIT   = NUMBER OF BITS TO USE IN PACKING DATA.
      DATA IBX,NBIT/ 0,12 /
!     
!*****************************************************************************
!     START GRIBIT HERE.

!     ALL TASKS MUST CALL COLLECT BUT ONLY TASK 0 CAN EXECUTE THE REMAINDER 
!      OF GRIBIT
      IM=GDIN%IMAX;JM=GDIN%JMAX;DATE=GDIN%DATE;IFHR=GDIN%FHR
      OGRD=GDIN%OGRD
      IFMIN=0
      IBUFLEN=30+LENPDS+LENGDS+IM*JM*(MXBIT+2)/8
      ALLOCATE(KBUF(IBUFLEN),STAT=kret)
      ALLOCATE(IGRD(IM,JM),IBMAP(IM,JM),STAT=kret)
      ALLOCATE(HOLDGRID(IM,JM),STAT=kret)

      IF (LUNOUT.EQ.70) THEN
!      OUTPUT FILE NAME LINKED TO GRID DOMAIN (CS,HI, PR...)
       DATSET='MESO'//TRIM(GDIN%REGION)
      ELSE IF (LUNOUT.EQ.71) THEN 
       DATSET='MAXMIN'
      ELSE
       DATSET='UNKOWN' 
      ENDIF
      IGRD=0.

      NEWFILE = .FALSE.
!     SET NUMBER OF OUTPUT GRID POINTS.
      IJOUT = IM*JM
!     
!     PREPARE GRIB PDS
!     
!     SET ARRAY ID VALUES TO GENERATE GRIB1 PDS.  
!        ID(1)  = NUMBER OF BYTES IN PRODUCT DEFINITION SECTION (PDS)
!        ID(2)  = PARAMETER TABLE VERSION NUMBER
!        ID(3)  = IDENTIFICATION OF ORIGINATING CENTER
!        ID(4)  = MODEL IDENTIFICATION (ALLOCATED BY ORIGINATING CENTER)
!        ID(5)  = GRID IDENTIFICATION
!        ID(6)  = 0 IF NO GDS SECTION, 1 IF GDS SECTION IS INCLUDED
!        ID(7)  = 0 IF NO BMS SECTION, 1 IF BMS SECTION IS INCLUDED
!        ID(8)  = INDICATOR OF PARAMETER AND UNITS (TABLE 2)
!        ID(9)  = INDICATOR OF TYPE OF LEVEL       (TABLE 3)
!        ID(10) = VALUE 1 OF LEVEL (=0 FOR 1-100,102,103,105,107,
!          109,111,113,115,117,119,125,160,200,201 LEVEL IS IN ID WORD 11)
!        ID(11) = VALUE 2 OF LEVEL
!        ID(12) = YEAR OF CENTURY
!        ID(13) = MONTH OF YEAR
!        ID(14) = DAY OF MONTH
!        ID(15) = HOUR OF DAY
!        ID(16) = MINUTE OF HOUR   (IN MOST CASES SET TO 0)
!        ID(17) = FCST TIME UNIT
!        ID(18) = P1 PERIOD OF TIME
!        ID(19) = P2 PERIOD OF TIME
!        ID(20) = TIME RANGE INDICATOR
!        ID(21) = NUMBER INCLUDED IN AVERAGE
!        ID(22) = NUMBER MISSING FROM AVERAGES
!        ID(23) = CENTURY
!        ID(24) = RESERVED - SET TO 0
!        ID(25) = SCALING POWER OF 10

!        PREPARE DATE PART OF GRIB PDS RECORD.

      iyy=int(date/1000000)
      icent=(iyy-1)/100 + 1
      imm=int(int(mod(date,1000000)/100)/100)
      idd=int(mod(date,10000)/100)
      ihrst=mod(date,100)

         AYEAR0     = IYY
         AMNTH0     = IMM
         ADAY0      = IDD
         AGMT0      = IHRST
         ID(01)     = 28
         IF (ID(2) .NE. 129.AND.ID(2).NE.130) ID(2)=2
         ID(03)     = 7
         ID(12)     = MOD(IYY,100)
         ID(13)     = IMM
         ID(14)     = IDD
         ID(15)     = IHRST
         ID(16)     = 0
         ID(17)     = 1

!    ASSUMING ID(18-20), (P1, P2, TIME RANGE INDICATOR) 
!    ARE PASSED IN CORRECTLY IF NOT AN INSTANTANEOUS FIELD
!   
         IF (ID(20).EQ.0) THEN
          ID(18)     = IFHR 
          ID(19)     = 0
         ENDIF
! CHUANG: TO OUTPUT OFF-HOUR FORECAST, I USED MIN INSTEAD OF HOUR AS FORECAST UNIT
! ALSO, SINCE ONLT TIME RANGE TYPE 10 USES 2 BYTES TO STORE TIME, MODIFICATION WAS
! MADE TO USE TYPE 10 AS TIME RANGE INDICATOE WHEN FORECST MINS ARE LARGER THAN 254,	
! WHICH MEANS ALL THE ACCUMULATED AND TIME-AVERAGED QUANTITY ARE VERIFIED AT ONE TIME
! INSTEAD OF AT A TIME RANGE. 
	 IF(IFMIN .GE. 1)THEN
	   ID(17)     = 0
	   TOTMIN=IFHR*60+IFMIN
	   IF(TOTMIN .LE. 256)THEN  	     
	     IF (ID(20).EQ.0)ID(18)=IFHR*60+IFMIN
           ELSE
	     ID(20)=10
	     ID(18)=IFHR*60+IFMIN     
	   END IF 
	  END IF

         ID(21)     = 0
         ID(22)     = 0
         ID(23)     = ICENT
         ID(24)     = 0

!     
!        SET OUTPUT GRID TYPE.  WE ASSUME KGYTPE HOLDS THE GRIB
!        ID FOR THE OUTPUT GRID.  

!         KGTYP = KGTYPE
          KGTYP = OGRD
!     
!        SET GRID TYPE ID(5)
!        GENERATING PROGRAM ID(4)

         ID(4) = 89 
         ID(5) = KGTYP

!        ID(6) =0 IF NO GDS SECTION, =1 IF GDS INCLUDED, 
!                 ALWAYS INCLUDE GDS

         ID(6) = 1
!     
!        SET DATA TYPE ID(8) AND SURFACE ID(9).


!     END OF GRIB PDS LABEL PREPARATION.

!     
!     SET DECIMAL SCALING (IDECI) FROM LIST IN INCLUDE FILE 
!     RQSTFLD.  A CALL TO GET_BITS WILL COMPUTE THE NUMBER OF
!     BITS NECESSARY TO PACK THE DATA BASED ON THE RANGE OF 
!     THE FIELD.  THE FIELD IS SCALED TO THIS PRECISION AND
!     RETURNED FOR PACKING BY THE GRIB PACKER.
!     
      IBITM = 0    
      IBM = 0
      SGDG  = DECI 
!     set bitmap
       HOLDGRID=GRID
       IBMAP=0       !2-d array
       where(abs(grid-spval).gt.small)ibmap=1
       ibitm=count(abs(grid-spval).gt.small)
!       print *,'IBITM',ibitm,'IJOUT',IJOUT

!        ID(7) =0 IF NO BMS SECTION, =1 IF BMS INCLUDED

      IF (IBITM.EQ.IJOUT) THEN
        ID(7) = 0
        IBM = 0
      ELSE
        ID(7) = 1
        IBM = 1
      ENDIF
      CALL GET_BITS(IBM,SGDG,IJOUT,IBMAP,GRID,IDECI,GRID,GMIN,GMAX,NBIT)

!     ID(25) = SCALING POWER OF 10
      ID(25) = IDECI
!     
!     GENERATE COMPLETE GRIB1 MESSAGE USING W3FI72.
!        ITYPE  = 0 SPECIFIES REAL DATA TO BE PACKED.
!        IGRD   = DUMMY ARRAY FOR INTEGER DATA.
!        IBITL  = NBIT TELLS W3FI72 TO PACK DATA USING NBIT BITS.
!        IPFLAG = 0 IS PDS INFORMATION IN USER ARRAY ID.
!                 1 IS PDS (GENERATED ABOVE BY W3FP12).
!        ID     = (DUMMY) ARRAY FOR USER DEFINED PDS.
!        IGFLAG = 0 TELLS W3FI72 TO MAKE GDS USING IGRID.
!                 1 IS GDS GENERATED BY USER IN ARRAY IGDS
!        IGRID  = GRIB1 GRID TYPE (TABLE B OF ON388).
!        IGDS   = ARRAY FOR USER DEFINED GDS.
!        ICOMP  = 0 FOR EARTH ORIENTED WINDS,
!                 1 FOR GRID ORIENTED WINDS.
!        IBFLAG = 0 TELLS W3FI72 TO MAKE BIT MAP FROM USER
!                 SUPPLIED DATA.
!        IBMASK = ARRAY CONTAINING USER DEFINED BIT MAP.
!        IBLEN  = LENGTH OF ARRAY IBMASK.
!        IBDSFL = ARRAY CONTAINING TABLE 11 (ON388) FLAG INFORMATION.
!        NPTS   = LENGTH OF ARRAY GRID OR IGRD.  MUST AGREE WITH IBLEN.
!     
!     INTIALIZE VARIABLES.

      ITYPE  = 0

      IBITL  = MIN(NBIT,MXBIT)

      IPFLAG = 0

!MEB  IGFLAG = 0
      IGFLAG = 1  ! set to 1 so that IGDS is defined here instead of w3lib
      IGRID  = ID(5)
      IF (IGRID.EQ.26) IGRID=6
      IGDS(1:18) = 0
      CALL W3FI71(OGRD,IGDS,IERR)
!       LAMBERT CONFORMAL:
!           IGDS( 1) = NUMBER OF VERTICAL COORDINATES
!           IGDS( 2) = PV, PL OR 255
!           IGDS( 3) = DATA REPRESENTATION TYPE (CODE TABLE 6)
!           IGDS( 4) = NO. OF POINTS ALONG X-AXIS
!           IGDS( 5) = NO. OF POINTS ALONG Y-AXIS
!           IGDS( 6) = LATITUDE OF ORIGIN (SOUTH -IVE)
!           IGDS( 7) = LONGITUTE OF ORIGIN (WEST -IVE)
!           IGDS( 8) = RESOLUTION FLAG (CODE TABLE 7)
!           IGDS( 9) = LONGITUDE OF MERIDIAN PARALLEL TO Y-AXIS
!           IGDS(10) = X-DIRECTION GRID LENGTH (INCREMENT)
!           IGDS(11) = Y-DIRECTION GRID LENGTH (INCREMENT)
!           IGDS(12) = PROJECTION CENTER FLAG (0=NORTH POLE ON PLANE,
!                                              1=SOUTH POLE ON PLANE,
!           IGDS(13) = SCANNING MODE FLAGS (CODE TABLE 8)
!           IGDS(14) = NOT USED
!           IGDS(15) = FIRST LATITUDE FROM THE POLE AT WHICH THE
!                      SECANT CONE CUTS THE SPERICAL EARTH
!           IGDS(16) = SECOND LATITUDE ...
!           IGDS(17) = LATITUDE OF SOUTH POLE (MILLIDEGREES)
!           IGDS(18) = LONGITUDE OF SOUTH POLE (MILLIDEGREES)

      ICOMP  = 1
      IF (INDEX(PROJ,'LOLA').NE.0) ICOMP = 0
      IBFLAG = 0
      IBLEN  = IJOUT
      IBDSFL(1:9) = 0

      CALL W3FI72(ITYPE,GRID,IGRD,IBITL,IPFLAG,ID,PDS, IGFLAG,IGRID,IGDS,   &
                  ICOMP,IBFLAG,IBMAP,IBLEN,IBDSFL,NPTS,KBUF,ITOT,IER)
      GRID=HOLDGRID
!     
!     EXPLICITLY SET BYTE 12 OF KBUF (BYTE 4 OF THE PDS)
!     TO 2.  THIS WILL REFER ALL QUANTITIES TO PARAMETER
!     TABLE VERSION 2 OF WHICH TABLE VERSION 1 IS A SUBSET.
!     THIS IS NEEDED BECAUSE THE W3 ROUTINES HARDWIRE THIS
!     VALUE TO 1 YET SOME OF THE OUTPUT VARIABLES ARE ONLY 
!     DEFINED IN VERSION 2 OF THE PARAMETER TABLE.

!--- Comment out; BYTE 4 (PDS Octet 4) = 2 or 129 (see ON388, Table 2)
!!      KBUF(12)=CHAR(2)

      IF (IER.NE.0) THEN
         WRITE(6,*)'GRIBIT:  W3FI72 ERROR DID NOT POST THIS FIELD'
         RETURN
      ENDIF
!     
!     ON FIRST ENTRY MAKE OUTPUT DIRECTORY.  SET SWITCH (RITEHD)
!     TO FALSE FOR SUBSEQUENT ENTRIES.
      IF (RITEHD) THEN

!        PUT FORECAST HOUR INTO DIR PREFIX FOR GRIB FILE.
         IHR = IFHR

!        GET FULL PATH FOR OUTPUT FILE FROM ENVIRONMENT VARIABLE
!        COMSP WHICH IS SET IN THE SCRIPT RUNNING THE MODEL.

!        CONSTRUCT FULL PATH-FILENAME FOR OUTPUT FILE
         CALL GETENV('COMSP',ENVAR)
         RESTHR='tm00'
         DESCR3=''
         CCOL='';CDOT='.'
	 IF(IFMIN .GE. 1) THEN
           CCOL=':'
	   WRITE(DESCR3,'(I2.2)') IFMIN
         ENDIF

!        CONSTRUCT FULL PATH-FILENAME FOR OUTPUT FILE

         WRITE(DESCR2,'(I2.2)') IHR
         IF (IHR.GE.100) WRITE(DESCR2,'(I3.3)') IHR

         FNAME = TRIM(ENVAR)//TRIM(DATSET)//TRIM(DESCR2)//TRIM(CCOL)//TRIM(DESCR3)//CDOT//RESTHR 

!        ASSIGN AND OPEN UNIT FOR GRIB DATA FILE.
         CLOSE(LUNOUT)
         CALL BAOPEN(LUNOUT,FNAME,IER)
         IF (IER.NE.0) WRITE(6,*) &
           'GRIBIT:  BAOPEN ERROR FOR GRIB DATA FILE.  IER=',IER
         WRITE(6,*)'GRIBIT:  OPENED ',LUNOUT,' FOR GRIB DATA ',FNAME
!     
!        SET OPEN-UNIT FLAGS TO FALSE.
         RITEHD = .FALSE.
         NEWFILE = .TRUE.
      ENDIF

!     WRITE GRIB1 MESSAGE TO OUTPUT FILE.
      CALL WRYTE(LUNOUT,ITOT,KBUF)

!     WRITE DIAGNOSTI! MESSAGE.
!        ID(8)  = INDICATOR OF PARAMETER AND UNITS (TABLE 2)
!        ID(9)  = INDICATOR OF TYPE OF LEVEL       (TABLE 3)
!        ID(10) = VALUE 1 OF LEVEL  (0 FOR 1-100,102,103,105,107
!              111,160   LEVEL IS IN ID WORD 11)
!        ID(11) = VALUE 2 OF LEVEL
       WRITE(6,1050) ID(8),ID(9),ID(10),ID(18),ID(19),MINVAL(GRID),MAXVAL(GRID)
 1050  FORMAT('GRIBIT:  ',5I5,2G10.3)

!     END OF ROUTINE.

      RETURN
      END
!        IGDS VARIES DEPENDING ON GRID REPRESENTATION TYPE.

!       LAT/LON GRID:
!           IGDS( 1) = NUMBER OF VERTICAL COORDINATES
!           IGDS( 2) = PV, PL OR 255
!           IGDS( 3) = DATA REPRESENTATION TYPE (CODE TABLE 6)
!           IGDS( 4) = NO. OF POINTS ALONG A LATITUDE
!           IGDS( 5) = NO. OF POINTS ALONG A LONGITUDE MERIDIAN
!           IGDS( 6) = LATITUDE OF ORIGIN (SOUTH - IVE)
!           IGDS( 7) = LONGITUDE OF ORIGIN (WEST -IVE)
!           IGDS( 8) = RESOLUTION FLAG (CODE TABLE 7)
!           IGDS( 9) = LATITUDE OF EXTREME POINT (SOUTH - IVE)
!           IGDS(10) = LONGITUDE OF EXTREME POINT (WEST - IVE)
!           IGDS(11) = LATITUDE INCREMENT
!           IGDS(12) = LONGITUDE INCREMENT
!           IGDS(13) = SCANNING MODE FLAGS (CODE TABLE 8)
!           IGDS(14) = ... THROUGH ...
!           IGDS(18) =   ... NOT USED FOR THIS GRID
!           IGDS(19) - IGDS(91) FOR GRIDS 37-44, NUMBER OF POINTS
!                      IN EACH OF 73 ROWS.

!       GAUSSIAN GRID:
!           IGDS( 1) = ... THROUGH ...
!           IGDS(10) =   ... SAME AS LAT/LON GRID
!           IGDS(11) = NUMBER OF LATITUDE LINES BETWEEN A POLE
!                      AND THE EQUATOR
!           IGDS(12) = LONGITUDE INCREMENT
!           IGDS(13) = SCANNING MODE FLAGS (CODE TABLE 8)
!           IGDS(14) = ... THROUGH ...
!           IGDS(18) =   ... NOT USED FOR THIS GRID

!       SPHERICAL HARMONICS:
!           IGDS( 1) = NUMBER OF VERTICAL COORDINATES
!           IGDS( 2) = PV, PL OR 255
!           IGDS( 3) = DATA REPRESENTATION TYPE (CODE TABLE 6)
!           IGDS( 4) = J - PENTAGONAL RESOLUTION PARAMETER
!           IGDS( 5) = K - PENTAGONAL RESOLUTION PARAMETER
!           IGDS( 6) = M - PENTAGONAL RESOLUTION PARAMETER
!           IGDS( 7) = REPRESENTATION TYPE (CODE TABLE 9)
!           IGDS( 8) = REPRESENTATION MODE (CODE TABLE 10)
!           IGDS( 9) = ... THROUGH ...
!           IGDS(18) =   ... NOT USED FOR THIS GRID

!       POLAR STEREOGRAPHIC:
!           IGDS( 1) = NUMBER OF VERTICAL COORDINATES
!           IGDS( 2) = PV, PL OR 255
!           IGDS( 3) = DATA REPRESENTATION TYPE (CODE TABLE 6)
!           IGDS( 4) = NO. OF POINTS ALONG X-AXIS
!           IGDS( 5) = NO. OF POINTS ALONG Y-AXIS
!           IGDS( 6) = LATITUDE OF ORIGIN (SOUTH -IVE)
!           IGDS( 7) = LONGITUTE OF ORIGIN (WEST -IVE)
!           IGDS( 8) = RESOLUTION FLAG (CODE TABLE 7)
!           IGDS( 9) = LONGITUDE OF MERIDIAN PARALLEL TO Y-AXIS
!           IGDS(10) = X-DIRECTION GRID LENGTH (INCREMENT)
!           IGDS(11) = Y-DIRECTION GRID LENGTH (INCREMENT)
!           IGDS(12) = PROJECTION CENTER FLAG (0=NORTH POLE ON PLANE,
!                                              1=SOUTH POLE ON PLANE,
!           IGDS(13) = SCANNING MODE FLAGS (CODE TABLE 8)
!           IGDS(14) = ... THROUGH ...
!           IGDS(18) =   .. NOT USED FOR THIS GRID

!       MERCATOR:
!           IGDS( 1) = ... THROUGH ...
!           IGDS(12) =   ... SAME AS LAT/LON GRID
!           IGDS(13) = LATITUDE AT WHICH PROJECTION CYLINDER
!                        INTERSECTS EARTH
!           IGDS(14) = SCANNING MODE FLAGS
!           IGDS(15) = ... THROUGH ...
!           IGDS(18) =   .. NOT USED FOR THIS GRID

!       LAMBERT CONFORMAL:
!           IGDS( 1) = NUMBER OF VERTICAL COORDINATES
!           IGDS( 2) = PV, PL OR 255
!           IGDS( 3) = DATA REPRESENTATION TYPE (CODE TABLE 6)
!           IGDS( 4) = NO. OF POINTS ALONG X-AXIS
!           IGDS( 5) = NO. OF POINTS ALONG Y-AXIS
!           IGDS( 6) = LATITUDE OF ORIGIN (SOUTH -IVE)
!           IGDS( 7) = LONGITUTE OF ORIGIN (WEST -IVE)
!           IGDS( 8) = RESOLUTION FLAG (CODE TABLE 7)
!           IGDS( 9) = LONGITUDE OF MERIDIAN PARALLEL TO Y-AXIS
!           IGDS(10) = X-DIRECTION GRID LENGTH (INCREMENT)
!           IGDS(11) = Y-DIRECTION GRID LENGTH (INCREMENT)
!           IGDS(12) = PROJECTION CENTER FLAG (0=NORTH POLE ON PLANE,
!                                              1=SOUTH POLE ON PLANE,
!           IGDS(13) = SCANNING MODE FLAGS (CODE TABLE 8)
!           IGDS(14) = NOT USED
!           IGDS(15) = FIRST LATITUDE FROM THE POLE AT WHICH THE
!                      SECANT CONE CUTS THE SPERICAL EARTH
!           IGDS(16) = SECOND LATITUDE ...
!           IGDS(17) = LATITUDE OF SOUTH POLE (MILLIDEGREES)
!           IGDS(18) = LONGITUDE OF SOUTH POLE (MILLIDEGREES)

!       ARAKAWA SEMI-STAGGERED E-GRID ON ROTATED LAT/LON GRID
!           IGDS( 1) = NUMBER OF VERTICAL COORDINATES
!           IGDS( 2) = PV, PL OR 255
!           IGDS( 3) = DATA REPRESENTATION TYPE (CODE TABLE 6) [201]
!           IGDS( 4) = NI  - TOTAL NUMBER OF ACTUAL DATA POINTS
!                            INCLUDED ON GRID
!           IGDS( 5) = NJ  - DUMMY SECOND DIMENSION; SET=1
!           IGDS( 6) = LA1 - LATITUDE  OF FIRST GRID POINT
!           IGDS( 7) = LO1 - LONGITUDE OF FIRST GRID POINT
!           IGDS( 8) = RESOLUTION AND COMPONENT FLAG (CODE TABLE 7)
!           IGDS( 9) = LA2 - NUMBER OF MASS POINTS ALONG
!                            SOUTHERNMOST ROW OF GRID
!           IGDS(10) = LO2 - NUMBER OF ROWS IN EACH COLUMN
!           IGDS(11) = DI  - LONGITUDINAL DIRECTION INCREMENT
!           IGDS(12) = DJ  - LATITUDINAL  DIRECTION INCREMENT
!           IGDS(13) = SCANNING MODE FLAGS (CODE TABLE 8)
!           IGDS(14) = ... THROUGH ...
!           IGDS(18) = ... NOT USED FOR THIS GRID (SET TO ZERO)

!       ARAKAWA FILLED E-GRID ON ROTATED LAT/LON GRID
!           IGDS( 1) = NUMBER OF VERTICAL COORDINATES
!           IGDS( 2) = PV, PL OR 255
!           IGDS( 3) = DATA REPRESENTATION TYPE (CODE TABLE 6) [202]
!           IGDS( 4) = NI  - TOTAL NUMBER OF ACTUAL DATA POINTS
!                            INCLUDED ON GRID
!           IGDS( 5) = NJ  - DUMMY SECOND DIMENTION; SET=1
!           IGDS( 6) = LA1 - LATITUDE LATITUDE OF FIRST GRID POINT
!           IGDS( 7) = LO1 - LONGITUDE OF FIRST GRID POINT
!           IGDS( 8) = RESOLUTION AND COMPONENT FLAG (CODE TABLE 7)
!           IGDS( 9) = LA2 - NUMBER OF (ZONAL) POINTS IN EACH ROW
!           IGDS(10) = LO2 - NUMBER OF (MERIDIONAL) POINTS IN EACH
!                            COLUMN
!           IGDS(11) = DI  - LONGITUDINAL DIRECTION INCREMENT
!           IGDS(12) = DJ  - LATITUDINAL  DIRECTION INCREMENT
!           IGDS(13) = SCANNING MODE FLAGS (CODE TABLE 8)
!           IGDS(14) = ... THROUGH ...
!           IGDS(18) = ... NOT USED FOR THIS GRID

!       ARAKAWA STAGGERED E-GRID ON ROTATED LAT/LON GRID
!           IGDS( 1) = NUMBER OF VERTICAL COORDINATES
!           IGDS( 2) = PV, PL OR 255
!           IGDS( 3) = DATA REPRESENTATION TYPE (CODE TABLE 6) [203]
!           IGDS( 4) = NI  - NUMBER OF DATA POINTS IN EACH ROW
!           IGDS( 5) = NJ  - NUMBER OF ROWS
!           IGDS( 6) = LA1 - LATITUDE OF FIRST GRID POINT
!           IGDS( 7) = LO1 - LONGITUDE OF FIRST GRID POINT
!           IGDS( 8) = RESOLUTION AND COMPONENT FLAG (CODE TABLE 7)
!           IGDS( 9) = LA2 - CENTRAL LATITUDE
!           IGDS(10) = LO2 - CENTRAL LONGTITUDE
!           IGDS(11) = DI  - LONGITUDINAL DIRECTION INCREMENT
!           IGDS(12) = DJ  - LATITUDINAL  DIRECTION INCREMENT
!           IGDS(13) = SCANNING MODE FLAGS (CODE TABLE 8)
!           IGDS(14) = ... THROUGH ...
!           IGDS(18) = ... NOT USED FOR THIS GRID
