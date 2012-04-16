      SUBROUTINE READ_SORT_CTL (LCNTRL,MAXG,MAXF,K5PDS,MDLID,KGRIDA,
     &  TYPE,SCAL,NAMES,KGRID,IG1,KSMPR,KSMPO,ITILS,JTILS,
     &  NFILS,NGRDS,NFLDS,JPDS5,JPDS6,JPDS7,JPDS16,IOUTUN,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C   SUBPROGRAM: READ_SORT_CTL
C   PRGMMR: BALDWIN          ORG: NP22        DATE: 98-08-11  
C
C ABSTRACT: READ_SORT_CTL PROCESSES THE PRDGEN CONTROL FILE.
C
C PROGRAM HISTORY LOG:
C   98-08-11  BALDWIN     ORIGINATOR
C   99-03-25  BALDWIN     MODIFY TO ADD NUMBER OF TILES TO SMOOTHING
C                          PARAMETERS (10'S DIGIT)
C
C USAGE:  CALL READ_SORT_CTL (LCNTRL,MAXG,MAXF,K5PDS,MDLID,KGRIDA,
C    &  TYPE,SCAL,NAMES,KGRID,IG1,KSMPR,KSMPO,ITILS,JTILS,
C    &  NFILS,NGRDS,NFLDS,JPDS5,JPDS6,JPDS7,JPDS16,IOUTUN,IRET)
C
C   INPUT:
C         LCNTRL            INTEGER - UNIT NUMBER OF CONTROL FILE
C         MAXG              INTEGER - MAXIMUM NUMBER OF OUTPUT FILES 
C              ("GRIDS")
C         MAXF              INTEGER - MAXIMUM NUMBER OF GRIB PARAMETERS 
C              ("FIELDS")             IN CONTROL FILE
C
C   OUTPUT:
C         K5PDS(5,MAXF)     INTEGER - PDS ELEMENTS 9-12, 21 FOR EACH FIELD
C         MDLID(MAXG,MAXF)  INTEGER - MODEL NUMBER FOR EACH GRID/FIELD
C         KGRIDA(MAXG,MAXF) INTEGER - GRID NUMBER FOR EACH GRID/FIELD
C         TYPE(MAXG,MAXF)   CHAR*4  - WMO HEADER TYPE FOR EACH GRID/FIELD
C                                      H = INTERNATIONAL     
C                                      A = AWIPS              
C                                      X = NO HEADER           
C         SCAL(MAXG,MAXF)   REAL    - PACKING PRECISION FOR EACH GRID/FIELD
C         NAMES(MAXG,MAXF)  CHAR*16 - OUTPUT NAME FOR EACH GRID/FIELD
C         KGRID(MAXG)       INTEGER - SORTED LIST OF UNIQUE GRID NUMBERS
C         IG1               INTEGER - NUMBER OF UNIQUE OUTPUT GRIDS
C         KSMPR(MAXG,MAXF)  INTEGER - NUMBER OF PRE-INTERPOLATION SMOOTHING
C                                     PASSES FOR EACH GRID/FIELD
C         KSMPO(MAXG,MAXF)  INTEGER - NUMBER OF POST-INTERPOLATION SMOOTHING
C                                     PASSES FOR EACH GRID/FIELD
C         ITILS(MAXG,MAXF)  INTEGER - NUMBER OF TILES IN X DIRECTION
C         JTILS(MAXG,MAXF)  INTEGER - NUMBER OF TILES IN Y DIRECTION
C         NFILS(MAXF)       INTEGER - NUMBER OF OUTPUT FILES THAT EACH
C                                     FIELD WILL BE WRITTEN TO
C         NGRDS(MAXF)       INTEGER - NUMBER OF OUTPUT GRIDS THAT EACH
C                                     FIELD WILL BE INTERPOLATED TO
C                                     (COULD BE LESS THAN NFILS SINCE
C                                      CERTAIN GRIDS ARE WRITTEN TO MORE
C                                      THAN ONE FILE)
C         NFLDS             INTEGER - NUMBER OF FIELDS READ
C         JPDS5(MAXF)       INTEGER - PDS OCTET 9 OF EACH FIELD
C         JPDS6(MAXF)       INTEGER - PDS OCTET 10 OF EACH FIELD
C         JPDS7(MAXF)       INTEGER - PDS OCTET 11-12 OF EACH FIELD
C         JPDS16(MAXF)      INTEGER - PDS OCTET 21 OF EACH FIELD
C         IOUTUN(MAXG,MAXF) INTEGER - OUTPUT UNIT NUMBER OF EACH GRID/FIELD
C         IRET              INTEGER - RETURN CODE
C
C   ENVIRONMENT VARIABLES: (OPTIONAL)
C           COMSP   - PATH PREFIX OF OUTPUT FILE NAMES
C           fhr     - FORECAST HOUR APPENDED TO FILE NAME
C           tmmark  - TIME MARK (tm00) APPENDED TO FILE NAME
C
C   SUBPROGRAMS CALLED:
C     UNIQUE:
C       CTL_RDR - CONTROL FILE READER
C
C   RETURN CODES:
C     IRET =   0 - NORMAL EXIT
C             -1 - INVALID FILE ENTRY
C             -2 - MISPLACED GRID/TYPE INFO 
C             -3 - INVALID GRID/TYPE ENTRY
C             -4 - INVALID GRID #  
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE : CRAY J-916
C
C$$$
C    
C    
      INTEGER KGRIDA(MAXG,MAXF),KGRID(MAXG)
      INTEGER JPDS5(MAXF),JPDS6(MAXF),
     &        JPDS7(MAXF),JPDS16(MAXF)
      INTEGER K5PDS(5,MAXF)
      INTEGER NGRDS(MAXF),NFILS(MAXF),MDLID(MAXG,MAXF)
      INTEGER IOUTUN(MAXG,MAXF)
      INTEGER KSMPO(MAXG,MAXF),KSMPR(MAXG,MAXF)
      INTEGER ITILS(MAXG,MAXF),JTILS(MAXG,MAXF)

      LOGICAL MATCH

      REAL SCAL(MAXG,MAXF)

      CHARACTER NAMES(MAXG,MAXF)*16,FNAME1(MAXG)*16
      CHARACTER FSTHR*16,RESTHR*16,ENVAR*60,DATSET*16,FNAME*80
      CHARACTER TYPE(MAXG,MAXF)*4
C
C
C
      IRET=0
C
C     READ CNTRL FILE.
C     
      DO N=1,MAXF
       NFILS(N)=0
       NGRDS(N)=0
      ENDDO

      CALL CTL_RDR  (LCNTRL, MAXG, K5PDS, MDLID, KGRIDA, TYPE,
     &     SCAL, NAMES, KSMPR, KSMPO, NFILS, NFLDS, IRET0)
      DO I=1,MAXG
      DO J=1,MAXF
       ITILS(I,J)=KSMPR(I,J)/10
       JTILS(I,J)=KSMPO(I,J)/10
       KSMPR(I,J)=MOD(KSMPR(I,J),10)
       KSMPO(I,J)=MOD(KSMPO(I,J),10)
      ENDDO
      ENDDO
      IRET=IRET0

      IF (IRET.NE.0) RETURN
C           
C
C

       DO N=1,NFLDS
        JPDS5(N)=K5PDS(1,N)
        JPDS6(N)=K5PDS(2,N)
        JPDS7(N)=K5PDS(3,N)*256+K5PDS(4,N)
        JPDS16(N)=K5PDS(5,N)
       ENDDO

C
C  MAKE GDS FOR OUTPUT GRIDS
C   ALSO KEEP TRACK OF NUMBER OF UNIQUE GRIDS, THIS WILL
C   PROBABLY NOT BE THE SAME AS THE NUMBER OF OUTPUT FILES
C
      DO M=1,MAXG
       KGRID(M)=0
      ENDDO
      IG1=0
      DO N=1,NFLDS
       NG1=NFILS(N)
       DO M=1,NG1
C
        JJ=0
        DO L=1,IG1
        IF (KGRIDA(M,N).EQ.KGRID(L)) JJ=JJ+1
        ENDDO
        IF (JJ.EQ.0) THEN
          IG1=IG1+1
          KGRID(IG1)=KGRIDA(M,N)
        ENDIF
C
        MATCH=.FALSE.
        DO L=1,M-1
        IF (KGRIDA(M,N).EQ.KGRIDA(L,N)) MATCH=.TRUE.
        ENDDO
        IF (.NOT.MATCH) NGRDS(N)=NGRDS(N)+1
C
       ENDDO
      ENDDO
C
C  FIGURE OUT UNIQUE FILE NAMES AND ASSIGN UNIT NUMBERS TO EACH FIELD
C
      INAM1=0
      DO N=1,NFLDS
       NG1=NFILS(N)
       DO M=1,NG1
        MATCH=.FALSE.
        DO L=1,INAM1
        IF (NAMES(M,N).EQ.FNAME1(L)) THEN
          MATCH=.TRUE.
          IOUTUN(M,N)=45+L
        ENDIF
        ENDDO
C
C  NEW FILE NAME
C
        IF (.NOT.MATCH) THEN
          INAM1=INAM1+1
          FNAME1(INAM1)=NAMES(M,N)
          IOUTUN(M,N)=45+INAM1
        ENDIF
       ENDDO
      ENDDO
C
C  OPEN OUTPUT FILES
C
C     
C        GET FULL PATH FOR OUTPUT FILE FROM ENVIRONMENT VARIABLE
C        COMSP WHICH IS SET IN THE SCRIPT RUNNING THE MODEL.
C        RESTHR GETS APPENDED TO THE NAME
C        FSTHR IS THE FORECAST HOUR (THIS IS REQUIRED)
C     
C        CONSTRUCT FULL PATH-FILENAME FOR OUTPUT FILE
         ENVAR = ' '
         RESTHR = ' '
         FSTHR = ' '
         CALL GETENV('COMSP',ENVAR)
         CALL GETENV('tmmark',RESTHR)
         CALL GETENV('fhr',FSTHR)
         KENV = INDEX(ENVAR,' ') -1
         IF (KENV.LE.0) KENV = LEN(ENVAR)
         KTHR = INDEX(RESTHR,' ') -1
         IF (KTHR.LE.0) KTHR = LEN(RESTHR)
         KFHR = INDEX(FSTHR,' ') -1
         IF (KFHR.LE.0) KFHR = LEN(FSTHR)

         DO L=1,INAM1
          LUNOUT=45+L
          DATSET=FNAME1(L)
          KDAT = INDEX(DATSET,' ') -1
          IF (KDAT.LE.0) KDAT = LEN(DATSET)
C     
C         CONSTRUCT FULL PATH-FILENAME FOR OUTPUT FILE
          IF (ENVAR(1:4).EQ.'    ') THEN
           IF (RESTHR(1:4).EQ.'    ') THEN
            FNAME = DATSET(1:KDAT) // FSTHR(1:KFHR)
           ELSE
            FNAME = DATSET(1:KDAT) // FSTHR(1:KFHR) // '.' // RESTHR
           ENDIF
          ELSE
           IF (RESTHR(1:4).EQ.'    ') THEN
            FNAME = ENVAR(1:KENV) // DATSET(1:KDAT) // FSTHR(1:KFHR)
           ELSE
            FNAME = ENVAR(1:KENV) // DATSET(1:KDAT) // FSTHR(1:KFHR)
     &              //'.'// RESTHR
           ENDIF
          ENDIF
C
C           OPEN UNIT FOR GRIB DATA FILE.
C
          CALL BAOPENW(LUNOUT,FNAME,IER)
          IF (IER.NE.0) THEN
            WRITE(6,*) '  ERROR = ',IER,' OPENING OUTPUT GRIB DATA ',
     X        'FILE  ',FNAME
          ELSE
            WRITE(6,*)' OPENED ',LUNOUT,' FOR GRIB DATA ',FNAME
          ENDIF
         ENDDO
C
C
C  SORT GRIDS BY GRID NUMBER
C
      DO M=1,IG1-1
       DO MM=M+1,IG1
        IF (KGRID(M).GT.KGRID(MM)) THEN
         IGRID1=KGRID(M)
         KGRID(M)=KGRID(MM)
         KGRID(MM)=IGRID1
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END
