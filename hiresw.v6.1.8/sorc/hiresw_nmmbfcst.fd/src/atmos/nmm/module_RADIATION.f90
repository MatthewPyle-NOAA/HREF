

!-----------------------------------------------------------------------
!
      MODULE MODULE_RADIATION
!
!-----------------------------------------------------------------------
!
!***  THE RADIATION DRIVERS AND PACKAGES
!
!---------------------
!--- Modifications ---
!---------------------
! 2010-04-02 Vasic - Removed WFR driver
!-----------------------------------------------------------------------
!
      USE MODULE_INCLUDE
      USE MODULE_MY_DOMAIN_SPECS
      USE MODULE_RA_GFDL,ONLY   : GFDL,CAL_MON_DAY,ZENITH
      USE MODULE_RA_RRTM,ONLY   : RRTM
      USE MODULE_CONSTANTS,ONLY : CAPPA,CP,EPSQ,G,P608,PI,R_D,STBOLT
!
      USE MODULE_DM_PARALLEL,ONLY : LOOPLIMITS
!
      USE MODULE_CONTROL,ONLY : NMMB_FINALIZE

      use module_radiation_driver_nmmb,  only : radupdate_nmmb
      use machine, only : kind_phys

!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: RADIATION
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE RADIATION PACKAGE OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  Shortwave
!
      INTEGER,PARAMETER  :: GFDLSWSCHEME=99                             &  !<--- (GFDL)
                           ,SWRADSCHEME=1                               &  !<--- (Dudhia, WRF)
                           ,GSFCSWSCHEME=2                              &  !<--- (Goddard, WRF)
                           ,RRTMSWSCHEME=3                                 !<--- (RRTM)
!
!***  Longwave
!
      INTEGER,PARAMETER  :: GFDLLWSCHEME=99                             &  !<--- (GFDL)
                           ,RRTMLWSCHEME=3                                 !<--- (RRTM)
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE RADIATION(ITIMESTEP,DT,JULDAY,JULYR,XTIME,JULIAN       &
     &                    ,IHRST,NPHS,GLAT,GLON                         &
     &                    ,NRADS,NRADL                                  &
     &                    ,DSG2,SGML2,SG2,PDSG1,PSGML1,PSG1,PT,PD       &
     &                    ,T,Q                                          &
     &                    ,THS,ALBEDO                                   &
     &                    ,QV,QC,QR,QI,QS,QG,NI                         &
     &                    ,F_QV,F_QC,F_QR,F_QI,F_QS,F_QG,F_NI           &
     &                    ,NUM_WATER                                    &
     &                    ,SM,CLDFRA                                    &
     &                    ,RLWTT,RSWTT                                  &
     &                    ,RLWIN,RSWIN                                  &
     &                    ,RSWINC,RSWOUT                                &
     &                    ,RLWTOA,RSWTOA                                &
     &                    ,CZMEAN,SIGT4                                 &
     &                    ,CFRACL,CFRACM,CFRACH                         &
     &                    ,ACFRST,NCFRST                                &
     &                    ,ACFRCV,NCFRCV                                &
     &                    ,CUPPT,SNOW                                   &
     &                    ,HTOP,HBOT                                    &
     &                    ,SHORTWAVE,LONGWAVE                           &
!
     &                    ,DT_INT,JDAT                                  &
     &                    ,CW,O3                                        &
     &                    ,F_ICE,F_RAIN                                 &
     &                    ,F_RIMEF                                      &
     &                    ,SI,TSKIN                                     &
     &                    ,Z0,SICE                                      &
     &                    ,MXSNAL,SGM                                   &
     &                    ,STDH,OMGALF                                  &
     &                    ,SNOWC                                        &
     &                    ,LM)
!-----------------------------------------------------------------------
!***  NOTE ***
! RLWIN  - downward longwave at the surface (=GLW)
! RSWIN  - downward shortwave at the surface (=XXX)
! RSWINC - CLEAR-SKY downward shortwave at the surface (=SWDOWN, new for AQ)
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    RADIATION   RADIATION OUTER DRIVER
!   PRGRMMR: BLACK           ORG: W/NP22     DATE: 2002-06-04       
!     
! ABSTRACT:
!     RADIATION SERVES AS THE INTERFACE BETWEEN THE NMMB PHYSICS COMPONENT
!     AND THE WRF RADIATION DRIVER.
!     
! PROGRAM HISTORY LOG:
!   02-06-04  BLACK      - ORIGINATOR
!   02-09-09  WOLFE      - CONVERTING TO GLOBAL INDEXING
!   04-11-18  BLACK      - THREADED
!   06-07-20  BLACK      - INCORPORATED INTO NMMB PHYSICS COMPONENT
!   08-08     JANJIC     - Synchronize WATER array and Q.
!   08-11-23  janjic     - general hybrid coordinate
!     
! USAGE: CALL RADIATION FROM PHY_RUN
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM 
!$$$  
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER,INTENT(IN) :: LM,DT_INT                                   &
                           ,IHRST,ITIMESTEP,JULDAY,JULYR                &
                           ,NPHS,NRADL,NRADS,NUM_WATER
!
      INTEGER,INTENT(IN) :: JDAT(8)
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: NCFRCV,NCFRST
!
      REAL,INTENT(IN) :: DT,JULIAN,PT,XTIME
!
      REAL,DIMENSION(1:LM),INTENT(IN) :: DSG2,PDSG1,PSGML1,SGML2
!
      REAL,DIMENSION(1:LM+1),INTENT(IN) :: PSG1,SG2
!
      REAL,DIMENSION(LM+1),INTENT(IN) :: SGM
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: ALBEDO,CUPPT        &
                                                   ,GLAT,GLON           &
                                                   ,PD,SM               &
                                                   ,SNOW,SNOWC,THS,SI   &
                                                   ,TSKIN,Z0,SICE       &
                                                   ,MXSNAL,STDH
!
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(IN) :: Q,T,CW,O3        &
                                                      ,F_ICE,F_RAIN     &
                                                      ,F_RIMEF,OMGALF
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACFRCV,ACFRST    &
                                                      ,RLWIN,RLWTOA     &
                                                      ,RSWIN,RSWOUT     &
                                                      ,HBOT,HTOP        &
                                                      ,RSWINC,RSWTOA
!
      REAL,DIMENSION(IMS:IME,JMS:JME,LM),INTENT(INOUT) :: RLWTT,RSWTT
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: CFRACH,CFRACL      &
                                                    ,CFRACM,CZMEAN      &
                                                    ,SIGT4
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(INOUT) :: QV,QR,QG
      REAL,DIMENSION(:,:,:),POINTER,INTENT(INOUT) :: QC,QI,QS,NI

!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(OUT) :: CLDFRA
!
      LOGICAL,INTENT(IN) :: F_QV,F_QC,F_QR,F_QI,F_QS,F_QG,F_NI
!
      CHARACTER(99),INTENT(IN) :: LONGWAVE,SHORTWAVE
!
!---------------------
!***  Local Variables
!---------------------
!
!.......................................................................
      INTEGER :: IQS,IQE,JQS,JQE   ! Same as ITS,ITE,JTS,JTE - Changed in looplimits
      INTEGER :: I_S,I_E,J_S,J_E   ! Also represent ITS,ITE,JTS,JTE





!.......................................................................
      INTEGER :: I,II,J,IJ,JDAY,JMONTH                                  &
                ,K,KMNTH,N,NRAD
!
      INTEGER :: LW_PHYSICS=0,SW_PHYSICS=0
!
      INTEGER,DIMENSION(3) :: IDAT
      INTEGER,DIMENSION(12) :: MONTH=(/31,28,31,30,31,30,31,31          &
     &                                ,30,31,30,31/)
!
      REAL :: DAYI,GMT,HOUR,PDSL,PLYR,RADT,TIMES,TDUM
!
      real (kind=kind_phys) :: SLAG, SDEC, CDEC, SOLCON, DTSW, DTX
!
      REAL,DIMENSION(1:LM) :: QL
!
      REAL,DIMENSION(IMS:IME,JMS:JME) :: GSW                            &
     &                                  ,TOT,TSFC,XLAND                 &
     &                                  ,GLW,SWDOWN,SWDOWNC,CZEN        &
     &                                  ,CUPPTR
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM+1) :: PHINT
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:LM) :: PI3D                      &
                                             ,THRATEN,THRATENLW         &
                                             ,THRATENSW
!
      LOGICAL :: GFDL_LW, GFDL_SW, LSSWR

      INTEGER :: jj, ip  ! used for 2D threading around RRTM


      integer(4) :: ic1, crate1, cmax1
      integer(4) :: ic2, crate2, cmax2

      call system_clock(count=ic1, count_rate=crate1, count_max=cmax1)
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!*****
!***** NOTE: THIS IS HARDWIRED FOR CALLS TO LONGWAVE AND SHORTWAVE
!*****       AT EQUAL INTERVALS
!*****
!-----------------------------------------------------------------------
!
      NRAD=NRADS
      RADT=DT*NRADS/60.
!
!-----------------------------------------------------------------------
!***  NOTE:  THE NMMB HAS IJK STORAGE WITH LAYER 1 AT THE TOP.
!***         THE WRF PHYSICS DRIVERS HAVE IKJ STORAGE WITH LAYER 1 
!***         AT THE BOTTOM.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                    &
!$omp private (j,i,k,pdsl,plyr,ql)
!.......................................................................
      DO J=JTS,JTE
      DO I=ITS,ITE
!
        PDSL=PD(I,J)
        XLAND(I,J)=SM(I,J)+1.
!
!-----------------------------------------------------------------------
!***  FILL THE SINGLE-COLUMN INPUT
!-----------------------------------------------------------------------
!
        DO K=1,LM
!
          PLYR=SGML2(K)*PDSL+PSGML1(K)
!
          QL(K)=AMAX1(Q(I,J,K),EPSQ)
!
          PHINT(I,J,K)=SG2(K)*PD(I,J)+PSG1(K)
          PI3D(I,J,K)=(PLYR*1.E-5)**CAPPA
!
          THRATEN(I,J,K)=0.
          THRATENLW(I,J,K)=0.
          THRATENSW(I,J,K)=0.
        ENDDO
!
        PHINT(I,J,LM+1)=SG2(LM+1)*PD(I,J)+PSG1(LM+1)
!
!-----------------------------------------------------------------------
!
        TSFC(I,J)=THS(I,J)*(PHINT(I,J,LM+1)*1.E-5)**CAPPA
!
      ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do 
!.......................................................................
!
      GMT=REAL(IHRST)
!
!.......................................................................
!$omp parallel do                                                    &
!$omp private (k,j,i)
!.......................................................................
        DO K=1,LM
        DO J=JMS,JME
        DO I=IMS,IME
          CLDFRA(I,J,K)=0.
        ENDDO
        ENDDO
        ENDDO
!
!.......................................................................
!$omp parallel do                                                    &
!$omp private (j,i)
!.......................................................................
      DO J=JMS,JME
        DO I=IMS,IME
          CFRACH(I,J)=0.
          CFRACL(I,J)=0.
          CFRACM(I,J)=0.
          CZMEAN(I,J)=0.
          SIGT4(I,J)=0.
          SWDOWN(I,J)=0.    ! TOTAL (clear+cloudy sky) shortwave down at the surface
          SWDOWNC(I,J)=0.   ! CLEAR SKY shortwave down at the surface
          GSW(I,J)=0.       ! Net (down - up) total (clear+cloudy sky) shortwave at the surface
          GLW(I,J)=0.       ! Total longwave down at the surface
          CUPPTR(I,J)=CUPPT(I,J)   ! Temporary array set to zero in radiation
        ENDDO
      ENDDO
!.......................................................................
!$omp end parallel do 
!.......................................................................
!
!-----------------------------------------------------------------------
!***  SYNCHRONIZE MIXING RATIO IN WATER ARRAY WITH SPECIFIC HUMIDITY.
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k)
!.......................................................................
      DO K=1,LM                                            
        DO J=JMS,JME                                      
          DO I=IMS,IME                                   
            QV(I,J,K)=Q(I,J,K)/(1.-Q(I,J,K))    
          ENDDO                                        
        ENDDO                                         
      ENDDO                                          
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!
!***  CALL THE INNER DRIVER.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!***  A PRIMARY MODIFICATION TO THE WRF DRIVER IS THE SPECIFICATION
!***  OF THE PACKAGES IN THE SELECT_CASE BLOCKS BEING CHANGED FROM
!***  INTEGERS (LISTED IN THE PHYSICS SECTION OF THE WRF REGISTRY)
!***  TO CHARACTERS (AS DEFINED IN THE ESMF CONFIG FILE).
!
!-----------------------------------------------------------------------
!***  TRANSLATE THE RADIATION OPTIONS IN THE CONFIG FILE TO THEIR
!***  ANALOGS IN THE WRF REGISTRY SO THAT THE WRF RADIATION DRIVER
!***  REMAINS UNTOUCHED.
!-----------------------------------------------------------------------
!
      SELECT CASE (TRIM(SHORTWAVE))
        CASE ('gfdl')
          SW_PHYSICS=99
        CASE ('dudh')
          SW_PHYSICS=1
        CASE ('gsfc')
          SW_PHYSICS=2
        CASE ('rrtm')
          SW_PHYSICS=3
        CASE DEFAULT
          WRITE(0,*)' User selected SHORTWAVE=',TRIM(SHORTWAVE)
          WRITE(0,*)' Improper selection of SW scheme in RADIATION'
          CALL NMMB_FINALIZE
      END SELECT

      SELECT CASE (TRIM(LONGWAVE))
        CASE ('gfdl')
          LW_PHYSICS=99
        CASE ('rrtm')
          LW_PHYSICS=3
        CASE DEFAULT
          WRITE(0,*)' User selected LONGWAVE=',TRIM(LONGWAVE)
          WRITE(0,*)' Improper selection of LW scheme in RADIATION'
          CALL NMMB_FINALIZE
      END SELECT


!==========================================================================
! Put "call radupdate_nmmb" here for threading safe
!==========================================================================
!---- for forcast purpose IDAT=JDAT

       DTSW =NRADS*DT                  ! [s]
       LSSWR=MOD(ITIMESTEP,NRADS)==0

      IF (SW_PHYSICS .EQ. 3 .or. LW_PHYSICS .EQ. 3) THEN
         DTX =DT
         call radupdate_nmmb                                          &
!  ---   inputs:
     &      ( JDAT, JDAT, DTSW, DTX, LSSWR, MYPE,                     &
!  ---   outputs:
     &        SLAG, SDEC, CDEC, SOLCON                                &
     &      )
      ENDIF

!==========================================================================
!==========================================================================


!
!-----------------------------------------------------------------------
!     CALL RADIATION_DRIVER
!-----------------------------------------------------------------------

   IF (LW_PHYSICS .EQ. 0 .AND. SW_PHYSICS .EQ. 0)         RETURN

   IF (ITIMESTEP .EQ. 1 .OR. MOD(ITIMESTEP,NRAD) .EQ. 0) THEN
     GFDL_LW = .FALSE.
     GFDL_SW = .FALSE.

!---------------

     IQS = ITS_B1
     IQE = ITE_B1
!.......................................................................







     JQS = JTS_B1
     JQE = JTE_B1

!-----------------------------------------------------------------------
!***  Initialize Data
!-----------------------------------------------------------------------
!
     DO J=JQS,JQE
     DO I=IQS,IQE
        GSW(I,J)=0.
        GLW(I,J)=0.
        SWDOWN(I,J)=0.
     ENDDO
     ENDDO
!
     DO K=1,LM
     DO J=JQS,JQE
     DO I=IQS,IQE
         THRATEN(I,J,K)=0.
     ENDDO
     ENDDO
     ENDDO
!
!-----------------------------------------------------------------------
!
     lwrad_gfdl_select: SELECT CASE(lw_physics)
!
!-----------------------------------------------------------------------

        CASE (GFDLLWSCHEME)

!-- Do nothing, since cloud fractions (with partial cloudiness effects) 
!-- are defined in GFDL LW/SW schemes and do not need to be initialized.

        CASE (RRTMLWSCHEME)

!-- Do nothing, since cloud fractions is calculated in RRTM

        CASE DEFAULT

          CALL CAL_CLDFRA(CLDFRA,                               &
                          QC,QI,F_QC,F_QI,                      &
                          IDS,IDE, JDS,JDE, 1,LM+1,             &
                          IMS,IME, JMS,JME, 1,LM+1,             &
                          IQS,IQE, JQS,JQE, 1,LM  )
!-----------------------------------------------------------------------

     END SELECT lwrad_gfdl_select    

!-----------------------------------------------------------------------
!
     lwrad_select: SELECT CASE(lw_physics)
!
!-----------------------------------------------------------------------
        CASE (RRTMLWSCHEME)


!
! The purpose of this logic is to divide the domain into tiles that
! are 8 elements in I and 1 element in J, giving potentially
! many more (greater concurrency) and smaller sized (better cache
! locality) tiles that may also be the width of the vector unit
! depending on the value of 8 (defined via CPP).  Dynamic
! thread scheduling is specified to help with load imbalance in 
! RRTM radiation.  The outer loop, chunk_loop_rrtm, is over the
! total number of tiles in the 2D subdomain. For each tile, ip,
! the J index (jj) by dividing the tile index by the number of tiles 
! in a row. That index is checked to make sure it falls within the
! extent of the subdomain in J, then the starting I index (ii) 
! of the tile is computed by taking the integer modulus of the tile
! index and the number of tiles in a row.  Finally, to avoid having
! more than one J-row at the start or end of the J-extent (which can
! happen because we're skipping the first and last J-row), which check
! to make sure that the start and end of this J-iteration (J_S and J_E)
! are identical (that is, we're only doing one row in each tile).
!
!$OMP DO PRIVATE (ip,ii,jj,i_s,i_e,j_s,j_e) schedule(dynamic)
          chunk_loop_rrtm:                                          &
          DO ip=1,((1+(ite-its+1)/8)*8)*(jte-jts+1) &
                                                         ,8
            jj=jts+(ip-1)/((1+(ite-its+1)/8)*8)
            j_in_range_rrtm:                                        &
            IF ((jj.ge.jts.and.jj.le.jte)         .AND.             &
                 ((JDS+1).LE.jj .AND. jj.LE.(JDE-1))) THEN
                ii=its+mod((ip-1),((1+(ite-its+1)/8)*       &
                                                       8))
              I_S = MAX(MAX(ii,ITS),IDS)
              I_E = MIN(MIN(ii+8-1,ITE),IDE)
              J_S = jj
              J_E = jj






              IF ( I_S .LE. I_E ) THEN
                CALL RRTM(ITIMESTEP,DT,JDAT                         &
                   ,NPHS,GLAT,GLON                                  &
                   ,NRADS,NRADL                                     &
                   ,DSG2,SGML2,PDSG1,PSGML1                         &
                   ,PT,PD                                           &
                   ,T,Q,CW,O3                                       &
                   ,ALBEDO                                          &
                   ,F_ICE,F_RAIN                                    &
                 ,QV,QC,QI,QS,QR,QG,NI                              &
                 ,F_QV,F_QC,F_QI,F_QS,F_QR,F_QG,F_NI                &
                 ,NUM_WATER                                         &
                   ,SM,CLDFRA                                       &
                   ,RLWTT,RSWTT                                     &
                   ,RLWIN,RSWIN                                     &
                   ,RSWINC,RSWOUT                                   &
                   ,RLWTOA,RSWTOA                                   &
                   ,CZMEAN,SIGT4                                    &
                   ,CFRACL,CFRACM,CFRACH                            &
                   ,ACFRST,NCFRST                                   &
                   ,ACFRCV,NCFRCV                                   &
                   ,CUPPT,SNOWC,SI                                  & ! was SNOW
                   ,HTOP,HBOT                                       &
                   ,TSKIN,Z0,SICE,F_RIMEF,MXSNAL,SGM,STDH,OMGALF    &
                   ,IMS,IME,JMS,JME                                 &
                   ,I_S,I_E,J_S,J_E                                 &
                   ,LM                                              &
                   ,SOLCON                                          &
                   ,MYPE )
              ENDIF

            ENDIF j_in_range_rrtm
          ENDDO chunk_loop_rrtm
!$OMP END DO



        CASE (GFDLLWSCHEME)

                 gfdl_lw  = .true.
                 CALL GFDL(                                         &
                  DT=dt,XLAND=xland                                 &
                 ,PHINT=phint,T=t                                   &
                 ,Q=Q                                               &
                 ,QW=QC                                             &
                 ,QI=QI                                             &
                 ,QS=QS                                             &
                 ,F_QC=F_QC,F_QI=F_QI,F_QS=F_QS                     &
                 ,TSK2D=tsfc,GLW=GLW,RSWIN=SWDOWN,GSW=GSW           &
                 ,RSWINC=SWDOWNC,CLDFRA=CLDFRA,PI3D=PI3D            &
                 ,GLAT=glat,GLON=glon,HTOP=htop,HBOT=hbot           &
                 ,ALBEDO=albedo,CUPPT=cupptr                        &
                 ,SNOW=snow,G=g,GMT=gmt                             &
                 ,NSTEPRA=nrad,NPHS=nphs,ITIMESTEP=itimestep        &
                 ,XTIME=xtime,JULIAN=julian                         &
                 ,JULYR=julyr,JULDAY=julday                         &
                 ,GFDL_LW=gfdl_lw,GFDL_SW=gfdl_sw                   &
                 ,CFRACL=cfracl,CFRACM=cfracm,CFRACH=cfrach         &
                 ,ACFRST=acfrst,NCFRST=ncfrst                       &
                 ,ACFRCV=acfrcv,NCFRCV=ncfrcv                       &
                 ,RSWTOA=rswtoa,RLWTOA=rlwtoa,CZMEAN=czmean         &
                 ,THRATEN=thraten,THRATENLW=thratenlw               &
                 ,THRATENSW=thratensw                               &
                 ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=1,KDE=lm+1  &     
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=1,KME=lm+1  &
                 ,ITS=iqs,ITE=iqe, JTS=jqs,JTE=jqe, KTS=1,KTE=lm    &
                                                                    )

        CASE DEFAULT
  
             WRITE(0,*)'The longwave option does not exist: lw_physics = ', lw_physics
             CALL NMMB_FINALIZE

!-----------------------------------------------------------------------
           
     END SELECT lwrad_select    

!-----------------------------------------------------------------------
!
     swrad_select: SELECT CASE(sw_physics)
!
!-----------------------------------------------------------------------

        CASE (SWRADSCHEME)
!!!          CALL SWRAD()

        CASE (GSFCSWSCHEME)
!!!          CALL GSFCSWRAD()

        CASE (RRTMSWSCHEME)

!-- Already called complete RRTM SW/LW scheme in LW part of driver
!!!          CALL RRTM()

        CASE (GFDLSWSCHEME)

                 gfdl_sw = .true.
                 CALL GFDL(                                         &
                  DT=dt,XLAND=xland                                 &
                 ,PHINT=phint,T=t                                   &
                 ,Q=Q                                               &
                 ,QW=QC                                             &
                 ,QI=QI                                             &
                 ,QS=QS                                             &
                 ,F_QC=F_QC,F_QI=F_QI,F_QS=F_QS                     &
                 ,TSK2D=tsfc,GLW=GLW,RSWIN=SWDOWN,GSW=GSW           &
                 ,RSWINC=SWDOWNC,CLDFRA=CLDFRA,PI3D=PI3D            &
                 ,GLAT=glat,GLON=glon,HTOP=htop,HBOT=hbot           &
                 ,ALBEDO=albedo,CUPPT=cupptr                        &
                 ,SNOW=snow,G=g,GMT=gmt                             &
                 ,NSTEPRA=nrad,NPHS=nphs,ITIMESTEP=itimestep        &
                 ,XTIME=xtime,JULIAN=julian                         &
                 ,JULYR=julyr,JULDAY=julday                         &
                 ,GFDL_LW=gfdl_lw,GFDL_SW=gfdl_sw                   &
                 ,CFRACL=cfracl,CFRACM=cfracm,CFRACH=cfrach         &
                 ,ACFRST=acfrst,NCFRST=ncfrst                       &
                 ,ACFRCV=acfrcv,NCFRCV=ncfrcv                       &
                 ,RSWTOA=rswtoa,RLWTOA=rlwtoa,CZMEAN=czmean         &
                 ,THRATEN=thraten,THRATENLW=thratenlw               &
                 ,THRATENSW=thratensw                               &
                 ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=1,KDE=lm+1  &
                 ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=1,KME=lm+1  &
                 ,ITS=iqs,ITE=iqe, JTS=jqs,JTE=jqe, KTS=1,KTE=lm    &
                                                                    )

        CASE DEFAULT

             WRITE(0,*)'The shortwave option does not exist: sw_physics = ', sw_physics
             CALL NMMB_FINALIZE

!-----------------------------------------------------------------------

     END SELECT swrad_select    

!-----------------------------------------------------------------------
!
!.......................................................................



!.......................................................................
!


!     CALL NMMB_FINALIZE

   ENDIF

!-----------------------------------------------------------------------
!
        IF(TRIM(SHORTWAVE)=='rrtm')THEN
      call system_clock(count=ic2, count_rate=crate2, count_max=cmax2)
      write(0,*)'RADIATION: ',ic2-ic1
!--- RRTM already calculated variables below
          RETURN
        ENDIF
!
!-----------------------------------------------------------------------
!
!***  UPDATE FLUXES AND TEMPERATURE TENDENCIES.
!
!-----------------------------------------------------------------------
!***  SHORTWAVE
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      IF(MOD(ITIMESTEP,NRADS)==0)THEN
!-----------------------------------------------------------------------
!
        IF(TRIM(SHORTWAVE)/='gfdl')THEN
!
!-----------------------------------------------------------------------
!***  COMPUTE CZMEAN FOR NON-GFDL SHORTWAVE
!-----------------------------------------------------------------------
!
          DO J=JMS,JME
          DO I=IMS,IME
            CZMEAN(I,J)=0.
            TOT(I,J)=0.
          ENDDO
          ENDDO
!
          CALL CAL_MON_DAY(JULDAY,JULYR,JMONTH,JDAY)
          IDAT(1)=JMONTH
          IDAT(2)=JDAY
          IDAT(3)=JULYR
!
          DO II=0,NRADS,NPHS
            TIMES=ITIMESTEP*DT+II*DT
            CALL ZENITH(TIMES,DAYI,HOUR,IDAT,IHRST,GLON,GLAT,CZEN       &
     &                 ,ITS,ITE,JTS,JTE                                 &
     &                 ,IDS,IDE,JDS,JDE,1,LM+1                          &
     &                 ,IMS,IME,JMS,JME,1,LM+1                          &
     &                 ,ITS,ITE,JTS,JTE,1,LM)
            DO J=JTS,JTE
            DO I=ITS,ITE
              IF(CZEN(I,J)>0.)THEN
                CZMEAN(I,J)=CZMEAN(I,J)+CZEN(I,J)
                TOT(I,J)=TOT(I,J)+1.
              ENDIF
            ENDDO
            ENDDO
!
          ENDDO
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            IF(TOT(I,J)>0.)CZMEAN(I,J)=CZMEAN(I,J)/TOT(I,J)
          ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!***  COMPUTE TOTAL SFC SHORTWAVE DOWN FOR NON-GFDL SCHEMES
!-----------------------------------------------------------------------
!
          DO J=JTS_B1,JTE_B1
          DO I=ITS_B1,ITE_B1
!
            SWDOWN(I,J)=GSW(I,J)/(1.-ALBEDO(I,J))  
!--- No value currently available for clear-sky solar fluxes from
!    non GFDL schemes, though it's needed for air quality forecasts.
!    For the time being, set to the total downward solar fluxes.
            SWDOWNC(I,J)=SWDOWN(I,J)
!
          ENDDO
          ENDDO
!
        ENDIF   !End non-GFDL/non-RRTM block
!-----------------------------------------------------------------------
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k)
!.......................................................................
        DO J=JTS_B1,JTE_B1
          DO I=ITS_B1,ITE_B1
!
            RSWIN(I,J)=SWDOWN(I,J)
            RSWINC(I,J)=SWDOWNC(I,J)
            RSWOUT(I,J)=SWDOWN(I,J)-GSW(I,J)
!
            DO K=1,LM
              RSWTT(I,J,K)=THRATENSW(I,J,K)*PI3D(I,J,K)
            ENDDO
!
          ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  LONGWAVE
!-----------------------------------------------------------------------
!
      IF(MOD(ITIMESTEP,NRADL)==0)THEN
!
!.......................................................................
!$omp parallel do                                                       &
!$omp& private(i,j,k,tdum)
!.......................................................................
        DO J=JTS_B1,JTE_B1
          DO I=ITS_B1,ITE_B1
!
            TDUM=T(I,J,LM)
            SIGT4(I,J)=STBOLT*TDUM*TDUM*TDUM*TDUM
!
            DO K=1,LM
              RLWTT(I,J,K)=THRATENLW(I,J,K)*PI3D(I,J,K)
            ENDDO
!
            RLWIN(I,J)=GLW(I,J)
!
          ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
      ENDIF

      call system_clock(count=ic2, count_rate=crate2, count_max=cmax2)
      write(0,*)'RADIATION: ',ic2-ic1

!
!-----------------------------------------------------------------------
!
      END SUBROUTINE RADIATION
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
   SUBROUTINE radconst(XTIME,DECLIN,SOLCON,JULIAN)
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------

! !ARGUMENTS:
   REAL, INTENT(IN   )      ::       XTIME,JULIAN
   REAL, INTENT(OUT  )      ::       DECLIN,SOLCON
   REAL, PARAMETER          ::       DEGRAD=3.1415926/180.
   REAL                     ::       OBECL,SINOB,SXLONG,ARG,  &
                                     DECDEG,DJUL,RJUL,ECCFAC
! ---- local variables -----
   REAL                     ::       DPD=360./365.
!
! !DESCRIPTION:
! Compute terms used in radiation physics 
! for short wave radiation

   DECLIN=0.
   SOLCON=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.
        
   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)
        
!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:
        
   IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)
   IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)

   SXLONG=SXLONG*DEGRAD
   ARG=SINOB*SIN(SXLONG)
   DECLIN=ASIN(ARG)
   DECDEG=DECLIN/DEGRAD
!----SOLAR CONSTANT ECCENTRICITY FACTOR (PALTRIDGE AND PLATT 1976)
   DJUL=JULIAN*360./365.
   RJUL=DJUL*DEGRAD
   ECCFAC=1.000110+0.034221*COS(RJUL)+0.001280*SIN(RJUL)+0.000719*  &
          COS(2*RJUL)+0.000077*SIN(2*RJUL)
   SOLCON=1370.*ECCFAC
   
   END SUBROUTINE radconst

!---------------------------------------------------------------------
   SUBROUTINE cal_cldfra(CLDFRA,QC,QI,F_QC,F_QI,                     &
          ids,ide, jds,jde, kds,kde,                                 &
          ims,ime, jms,jme, kms,kme,                                 &
          its,ite, jts,jte, kts,kte                                  )
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------
   INTEGER,  INTENT(IN   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

!
   REAL, DIMENSION( ims:ime, jms:jme, kts:kte ), INTENT(OUT  ) ::    &
                                                             CLDFRA

   REAL, DIMENSION( ims:ime, jms:jme, kts:kte ), INTENT(IN   ) ::    &
                                                                 QI, &
                                                                 QC

   LOGICAL,INTENT(IN) :: F_QC,F_QI

   REAL thresh
   INTEGER:: i,j,k
! !DESCRIPTION:
! Compute cloud fraction from input ice and cloud water fields
! if provided.
!
! Whether QI or QC is active or not is determined from the logical
! switches f_qi and f_qc. They are passed in to the routine
! to enable testing to see if QI and QC represent active fields.
! 
!---------------------------------------------------------------------
     thresh=1.0e-6

     IF ( f_qi .AND. f_qc ) THEN
!.......................................................................
!$omp parallel do private(k,j,i)
!.......................................................................
        DO k = kts,kte
        DO j = jts,jte
        DO i = its,ite
           IF ( QC(i,j,k)+QI(I,j,k) .gt. thresh) THEN
              CLDFRA(i,j,k)=1.
           ELSE
              CLDFRA(i,j,k)=0.
           ENDIF
        ENDDO
        ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
     ELSE IF ( f_qc ) THEN
!.......................................................................
!$omp parallel do private(k,j,i)
!.......................................................................
        DO k = kts,kte
        DO j = jts,jte
        DO i = its,ite
           IF ( QC(i,j,k) .gt. thresh) THEN
              CLDFRA(i,j,k)=1.
           ELSE
              CLDFRA(i,j,k)=0.
           ENDIF
        ENDDO
        ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
!
     ELSE 
!
!.......................................................................
!$omp parallel do private(k,j,i)
!.......................................................................
        DO k = kts,kte
        DO j = jts,jte
        DO i = its,ite
           CLDFRA(i,j,k)=0.
        ENDDO
        ENDDO
        ENDDO
!.......................................................................
!$omp end parallel do
!.......................................................................
     ENDIF

   END SUBROUTINE cal_cldfra
!---------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!---------------------------------------------------------------------
!
      END MODULE MODULE_RADIATION
!
!-----------------------------------------------------------------------
