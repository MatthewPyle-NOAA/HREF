      SUBROUTINE SERVER
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .
C SUBPROGRAM:    SERVER      PERFORMS IO TO DISK
C   PRGRMMR: TUCCILLO        ORG: IBM
C
C ABSTRACT:
C     THIS ROUTINE RECEIVES DATA FROM TASK 0 OF MPI_COMM_INTER,
C     THE FIRST TASK PERFORMING THE POST_PROCESSING, AND WRITES
C     THE DATA TO DISK
C   .
C
C PROGRAM HISTORY LOG:
C   01-06-15  TUCCILLO - ORIGINAL
C
C USAGE:    CALL SERVER
C   INPUT ARGUMENT LIST:
C     NONE
C
C   OUTPUT ARGUMENT LIST:
C     NONE
C
C   OUTPUT FILES:
C     WRITES TO FILE FNAME
C
C   SUBPROGRAMS CALLED:
C       MPI_RECV
C       BAOPEN
C       BACIO
C     UTILITIES:
C       NONE
C     LIBRARY:
C       COMMON   - CTLBLK.comm
C
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN
C     MACHINE : IBM RS/6000 SP
C$$$
      INCLUDE "CTLBLK.comm"
      INCLUDE 'mpif.h'
      LOGICAL :: DONE, NEWFILE
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR, COUNT, LUN
      CHARACTER*80 :: FNAME
      CHARACTER*1, ALLOCATABLE :: BUF(:)
C
C---------------------------------------------------------------------  
C
C     THIS CODE IS EXPECTING THE FOLLOWING MESSAGE STRUCTURE
C
C     VARIABLE     TYPE           DESCRIPTION     TAG
C=====================================================
C     DONE         LOGICAL        ARE WE DONE?    1
C     NEWFILE      LOGICAL        OPEN THE FILE?  2
C     LUN          INTEGER        FORTRAN UNIT #  3
C     FNAME        CHARACTER*80   FILE NAME       4
C     BUF          CHARACTER*1(*) BURF RECORD     5
C
C---------------------------------------------------------------------
C
      PRINT *, ' STARTING UP IO SERVER ...'
666   CONTINUE
C
C     THE FIRST MESSAGE IS A LOGICAL TO TELL US WHETHER WE ARE
C     FINISHED OR NOT
C
      CALL MPI_RECV(DONE,1,MPI_LOGICAL,
     *              0,1,MPI_COMM_INTER,STATUS,IERR)
C
      IF ( DONE ) THEN
         PRINT *, ' SHUTTING DOWN IO SERVER ...'
         RETURN   !    RETURNING TO MAIN
      END IF
C
C     DO WE NEED TO OPEN THE FILE ?
C
      CALL MPI_RECV(NEWFILE,1,MPI_LOGICAL,
     *              0,2,MPI_COMM_INTER,STATUS,IERR)
C
C     FORTRAN UNIT NUMBER
C
      CALL MPI_RECV(LUN,1,MPI_INTEGER,
     *              0,3,MPI_COMM_INTER,STATUS,IERR)
C
C     FILENAME
C
      CALL MPI_RECV(FNAME,80,MPI_CHARACTER,
     *              0,4,MPI_COMM_INTER,STATUS,IERR)
C
C     OPEN THE FILE, IF NECESSARY
C
      IF ( NEWFILE ) THEN
          CLOSE(LUN)
          CALL BAOPEN(LUN,FNAME,IER)
          PRINT *, ' FILE ',FNAME,' OPENED AS UNIT ',LUN
      END IF
C
C     DETERMINE THE SIZE OF THE BUFR RECORD AND ALLOCATE A BUFFER FOR IT
C
      CALL MPI_PROBE(0,5,MPI_COMM_INTER,STATUS,IERR)
      CALL MPI_GET_COUNT(STATUS,MPI_CHARACTER,COUNT,IERR)
      ALLOCATE( BUF( COUNT ) )
C   
C     FINALLY, GET THE BUFR RECORD
C
      CALL MPI_RECV(BUF,COUNT,MPI_CHARACTER,
     *              0,5,MPI_COMM_INTER,STATUS,IERR)
C
C     OUT TO DISK WE GO ...
C
      CALL WRYTE(LUN,COUNT,BUF)
      DEALLOCATE(BUF)
      GOTO 666
      END
