



















      MODULE MODULE_PHYSICS_CALLS


      USE MODULE_DOMAIN
      USE MODULE_DM
      USE MODULE_CONFIGURE
      USE MODULE_TILES
      USE MODULE_STATE_DESCRIPTION,ONLY : P_QV,P_QC,P_QR,P_QI,P_QS,P_QG,P_QNI,P_QNR
      USE MODULE_MODEL_CONSTANTS
      USE MODULE_RA_GFDLETA,ONLY : CAL_MON_DAY,ZENITH
      USE MODULE_RADIATION_DRIVER
      USE MODULE_SF_MYJSFC
      USE MODULE_SURFACE_DRIVER
      USE MODULE_PBL_DRIVER
      USE MODULE_GWD
      USE MODULE_CU_BMJ
      USE MODULE_CUMULUS_DRIVER
      USE MODULE_MP_ETANEW
      USE MODULE_MICROPHYSICS_DRIVER
      USE MODULE_MICROPHYSICS_ZERO_OUT


      CONTAINS



      SUBROUTINE RADIATION(NTSD,DT,JULDAY,JULYR,XTIME,JULIAN            &
     &                    ,IHRST,NPHS,GLAT,GLON                         &
     &                    ,NRADS,NRADL                                  &
     &                    ,DETA1,DETA2,AETA1,AETA2,ETA1,ETA2,PDTOP,PT   &
     &                    ,PD,RES,PINT,T,Q,MOIST,THS,ALBEDO,EPSR        &
     &                    ,F_ICE,F_RAIN                                 &
     &                    ,SM,HBM2,CLDFRA,N_MOIST,RESTRT                &
     &                    ,RLWTT,RSWTT,RLWIN,RSWIN,RSWINC,RSWOUT        &
     &                    ,RLWTOA,RSWTOA,CZMEAN                         &
     &                    ,CFRACL,CFRACM,CFRACH,SIGT4                   &
     &                    ,ACFRST,NCFRST,ACFRCV,NCFRCV                  &
     &                    ,CUPPT,VEGFRC,SNOW,HTOP,HBOT                  &
     &                    ,Z,SICE,NUM_AEROSOLC,NUM_OZMIXM               &
     &                    ,GRID,CONFIG_FLAGS                            &
     &                    ,RTHRATEN                                     &
     &                    ,IDS,IDE,JDS,JDE,KDS,KDE                      &
     &                    ,IMS,IME,JMS,JME,KMS,KME                      &
     &                    ,ITS,ITE,JTS,JTE,KTS,KTE)




























      IMPLICIT NONE



      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,IHRST,JULDAY,JULYR                          &
     &                     ,N_MOIST,NPHS,NRADL,NRADS,NTSD               &
     &                     ,NUM_AEROSOLC,NUM_OZMIXM

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: NCFRCV,NCFRST

      REAL,INTENT(IN) :: DT,PDTOP,PT,XTIME,JULIAN

      REAL,DIMENSION(KMS:KME-1),INTENT(IN) :: AETA1,AETA2,DETA1,DETA2

      REAL,DIMENSION(KMS:KME),INTENT(IN) :: ETA1,ETA2

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: ALBEDO              &
     &                                             ,EPSR,GLAT,GLON      &
     &                                             ,HBM2                &
     &                                             ,PD,RES,SICE,SM      &
     &                                             ,SNOW,THS,VEGFRC

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: CUPPT


      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: Q,T,Z

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(IN) :: F_ICE       &   
     &                                                     ,F_RAIN

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) :: RTHRATEN     

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,N_MOIST)                   &
                                                 ,INTENT(INOUT) :: MOIST

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACFRCV,ACFRST    &
     &                                                ,HBOT,HTOP        &
     &                                                ,RLWIN,RLWTOA     &
     &                                                ,RSWIN,RSWOUT     &
     &                                                ,RSWINC,RSWTOA

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: PINT     &
     &                                                        ,RLWTT    &
     &                                                        ,RSWTT

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: CFRACH,CFRACL    &
     &                                                ,CFRACM,CZMEAN    &
     &                                                ,SIGT4


      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: CLDFRA

      LOGICAL,INTENT(IN) :: RESTRT

      TYPE(DOMAIN),TARGET :: GRID

      TYPE(GRID_CONFIG_REC_TYPE),INTENT(IN) :: CONFIG_FLAGS






      INTEGER :: I,ICLOUD,IENDX,II,ISTAT,J,JDAY,JMONTH,K,KMNTH,N,NRAD

      INTEGER,DIMENSION(3) :: IDAT
      INTEGER,DIMENSION(12) :: MONTH=(/31,28,31,30,31,30,31,31          &
     &                                ,30,31,30,31/)

      REAL :: CAPA,DAYI,DPL,FICE,FRAIN,GMT,HOUR,PLYR,PSFC               &
     &       ,QI,QR,QW,RADT,TIMES,WC,TDUM

      REAL,DIMENSION(KMS:KME-1) :: QL,TL

      REAL,DIMENSION(IMS:IME,JMS:JME) :: CUPPTR,CZEN,HBOTR,HTOPR        &
     &                                  ,PDSL,REXNSFC,SWNETDN           &
     &                                  ,TOT,TOTLWDN,TOTSWDN,TOTSWDNC   &
     &                                  ,TSFC,XLAND,XLAT,XLON


      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: CLFR,DZ                &   
     &                                          ,P8W,P_PHY,PI_PHY       &
     &                                          ,RR,T8W                 &
     &                                          ,THRATENLW,THRATENSW    &
     &                                          ,TH_PHY,T_PHY,Z_PHY

      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: MOIST_TRANS

      LOGICAL :: WARM_RAIN










      NRAD=NRADS
      RADT=DT*NRADS/60.



      ALLOCATE(MOIST_TRANS(IMS:IME,KMS:KME,JMS:JME,N_MOIST),STAT=ISTAT)



      CAPA=R_D/CP



!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
      DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))

        PDSL(I,J)=PD(I,J)*RES(I,J)
        P8W(I,KTE+1,J)=PT
        XLAT(I,J)=GLAT(I,J)/DEGRAD
        XLON(I,J)=GLON(I,J)/DEGRAD
        XLAND(I,J)=SM(I,J)+1.
        PSFC=PD(I,J)+PDTOP+PT
        REXNSFC(I,J)=(PSFC*1.E-5)**CAPA
        TSFC(I,J)=THS(I,J)*REXNSFC(I,J)
        T8W(I,KTS,J)=TSFC(I,J)
        P8W(I,KTS,J)=ETA1(KTS)*PDTOP+ETA2(KTS)*PDSL(I,J)+PT
        Z_PHY(I,KTS,J)=Z(I,J,KTS)
      ENDDO
      ENDDO





!$omp parallel do                                                       &
!$omp& private(dpl,i,j,k,plyr,ql,qr,tl)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
      DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
        DO K=KTS,KTE
          DPL=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
          QL(K)=MAX(Q(I,J,K),EPSQ)
          PLYR=AETA1(K)*PDTOP+AETA2(K)*PDSL(I,J)+PT
          TL(K)=T(I,J,K)

          RR(I,K,J)=PLYR/(R_D*TL(K)*(1.+P608*QL(K)))
          T_PHY(I,K,J)=TL(K)
          TH_PHY(I,K,J)=TL(K)*(1.E5/PLYR)**CAPA
          P8W(I,K+1,J)=ETA1(K+1)*PDTOP+ETA2(K+1)*PDSL(I,J)+PT
          P_PHY(I,K,J)=PLYR
          PI_PHY(I,K,J)=(PLYR*1.E-5)**CAPA
          DZ(I,K,J)=TL(K)*(P608*QL(K)+1.)*R_D                           &
     &                 *(P8W(I,K,J)-P8W(I,K+1,J))                       &
     &                 /(P_PHY(I,K,J)*G)

          RTHRATEN(I,K,J)=0.
          THRATENLW(I,K,J)=0.
          THRATENSW(I,K,J)=0.



        ENDDO

        DO K=KTS+1,KTE
          T8W(I,K,J)=0.5*(TL(K-1)+TL(K))
        ENDDO
        T8W(I,KTE+1,J)=-1.E20

      ENDDO
      ENDDO

      ICLOUD=999

      GMT=REAL(IHRST)

!$omp parallel do                                                       &
!$omp& private(i,j,k)
      DO K=KMS,KME
        DO J=JMS,JME
        DO I=IMS,IME
          CLDFRA(I,J,K)=0.
        ENDDO
        ENDDO
      ENDDO

!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=JMS,JME
        DO I=IMS,IME
          CFRACH(I,J)=0.
          CFRACL(I,J)=0.
          CFRACM(I,J)=0.
          CZMEAN(I,J)=0.
          SIGT4(I,J)=0.
          TOTSWDN(I,J)=0.   
          TOTSWDNC(I,J)=0.  
          SWNETDN(I,J)=0.   
          TOTLWDN(I,J)=0.   
          CUPPTR(I,J)=CUPPT(I,J)   







        ENDDO
      ENDDO








      DO N=1,N_MOIST
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO K=KMS,KME
        DO J=JMS,JME
        DO I=IMS,IME
          MOIST_TRANS(I,K,J,N)=MOIST(I,J,K,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO







      CALL SET_TILES(GRID,IDS+1,IDE-1,JDS+2,JDE-2,ITS,ITE,JTS,JTE)

CALL RADIATION_DRIVER(IDS=IDS,IDE=IDE,JDS=JDS,JDE=JDE,KDS=KDS,KDE=KDE,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=KMS,KME=KME,I_START=GRID&
&%I_START,I_END=GRID%I_END,J_START=GRID%J_START,J_END=GRID%J_END,KTS=KTS,KTE=KTE,NUM_TILES=GRID%NUM_TILES,ITIMESTEP=NTSD,DT=DT,RTHR&
&ATENLW=THRATENLW,RTHRATENSW=THRATENSW,RTHRATEN=RTHRATEN,GLW=TOTLWDN,GSW=SWNETDN,SWDOWN=TOTSWDN,XLAT=XLAT,XLONG=XLON,ALBEDO=ALBEDO,&
&EMISS=EPSR,XICE=SICE,XLAND=XLAND,Z=Z,TSK=TSFC,N_AEROSOLC=NUM_AEROSOLC,PAERLEV=GRID%PAERLEV,CAM_ABS_DIM1=GRID%CAM_ABS_DIM1,CAM_ABS_&
&DIM2=GRID%CAM_ABS_DIM2,CAM_ABS_FREQ_S=GRID%CAM_ABS_FREQ_S,LEVSIZ=GRID%LEVSIZ,N_OZMIXM=NUM_OZMIXM,HTOP=HTOP,HBOT=HBOT,CUPPT=CUPPTR,&
&HTOPR=HTOPR,HBOTR=HBOTR,VEGFRA=VEGFRC,SNOW=SNOW,RHO=RR,P8W=P8W,P=P_PHY,PI=PI_PHY,DZ8W=DZ,T=T_PHY,T8W=T8W,GMT=GMT,JULDAY=JULDAY,JUL&
&YR=JULYR,NPHS=NPHS,JULIAN=JULIAN,XTIME=XTIME,LW_PHYSICS=CONFIG_FLAGS%RA_LW_PHYSICS,SW_PHYSICS=CONFIG_FLAGS%RA_SW_PHYSICS,RADT=RADT&
&,RA_CALL_OFFSET=GRID%RA_CALL_OFFSET,STEPRA=NRAD,ICLOUD=ICLOUD,WARM_RAIN=WARM_RAIN,SWDOWNC=TOTSWDNC,CLDFRA=CLFR,RSWTOA=RSWTOA,RLWTO&
&A=RLWTOA,CZMEAN=CZMEAN,CFRACL=CFRACL,CFRACM=CFRACM,CFRACH=CFRACH,ACFRST=ACFRST,NCFRST=NCFRST,ACFRCV=ACFRCV,NCFRCV=NCFRCV,F_ICE_PHY&
&=F_ICE,F_RAIN_PHY=F_RAIN,QV=MOIST_TRANS(IMS,KMS,JMS,P_QV),F_QV=F_QV,QC=MOIST_TRANS(IMS,KMS,JMS,P_QC),F_QC=F_QC,QR=MOIST_TRANS(IMS,&
&KMS,JMS,P_QR),F_QR=F_QR,QI=MOIST_TRANS(IMS,KMS,JMS,P_QI),F_QI=F_QI,QS=MOIST_TRANS(IMS,KMS,JMS,P_QS),F_QS=F_QS,QG=MOIST_TRANS(IMS,K&
&MS,JMS,P_QG),F_QG=F_QG)











      nrads_block: IF(MOD(NTSD,NRADS)==0)THEN


        IF(CONFIG_FLAGS%RA_SW_PHYSICS/=GFDLSWSCHEME)THEN





!$omp parallel do                                                       &
!$omp& private(i,j)
          DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
          DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
            CZMEAN(I,J)=0.
            TOT(I,J)=0.
          ENDDO
          ENDDO

          CALL CAL_MON_DAY(JULDAY,JULYR,JMONTH,JDAY)
          IDAT(1)=JMONTH
          IDAT(2)=JDAY
          IDAT(3)=JULYR

          DO II=0,NRADS,NPHS
            TIMES=NTSD*DT+II*DT
            CALL ZENITH(TIMES,DAYI,HOUR,IDAT,IHRST,GLON,GLAT,CZEN       &
     &                 ,max(ids+(0),its-(0)),min(ide-(0),ite+(0)),max(jds+(0),jts-(0)),min(jde-(0),jte+(0))                             &
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                         &
     &                 ,IMS,IME,JMS,JME,KMS,KME                         &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE)

!$omp parallel do                                                       &
!$omp& private(i,j)
            DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
            DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
              IF(CZEN(I,J)>0.)THEN
                CZMEAN(I,J)=CZMEAN(I,J)+CZEN(I,J)
                TOT(I,J)=TOT(I,J)+1.
              ENDIF
            ENDDO
            ENDDO

          ENDDO

!$omp parallel do                                                       &
!$omp& private(i,j)
          DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
          DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
            IF(TOT(I,J)>0.)CZMEAN(I,J)=CZMEAN(I,J)/TOT(I,J)
          ENDDO
          ENDDO





!$omp parallel do                                                       &
!$omp& private(i,j)
          DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))

            IF(HBM2(I,J)>0.5)THEN
              TOTSWDN(I,J)=SWNETDN(I,J)/(1.-ALBEDO(I,J))





              TOTSWDNC(I,J)=TOTSWDN(I,J)
            ENDIF

          ENDDO
          ENDDO

        ENDIF   


!$omp parallel do                                                       &
!$omp& private(i,iendx,j)
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IENDX=min(ide-(1),ite+(0))
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=max(ids+(1),its-(0)),IENDX

            RSWIN(I,J)=TOTSWDN(I,J)
            RSWINC(I,J)=TOTSWDNC(I,J)
            RSWOUT(I,J)=TOTSWDN(I,J)-SWNETDN(I,J)

          ENDDO
        ENDDO

!$omp parallel do                                                       &
!$omp& private(i,iendx,j,k)
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IENDX=min(ide-(1),ite+(0))
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=max(ids+(1),its-(0)),IENDX
            DO K=KTS,KTE
              RSWTT(I,J,K)=THRATENSW(I,K,J)*PI_PHY(I,K,J)
            ENDDO

          ENDDO
        ENDDO

      ENDIF nrads_block





      nradl_block: IF(MOD(NTSD,NRADL)==0)THEN

!$omp parallel do                                                       &
!$omp& private(i,iendx,j)
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IENDX=min(ide-(1),ite+(0))
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=max(ids+(1),its-(0)),IENDX

            IF(HBM2(I,J)>0.5)THEN
              TDUM=T(I,J,KTS)
              SIGT4(I,J)=STBOLT*TDUM*TDUM*TDUM*TDUM
              RLWIN(I,J)=TOTLWDN(I,J)
            ENDIF

          ENDDO
        ENDDO

!$omp parallel do                                                       &
!$omp& private(i,iendx,j,k)
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IENDX=min(ide-(1),ite+(0))
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1

          DO K=KTS,KTE
          DO I=max(ids+(1),its-(0)),IENDX
            IF(HBM2(I,J)>0.5)THEN
                RLWTT(I,J,K)=THRATENLW(I,K,J)*PI_PHY(I,K,J)
            ENDIF
          ENDDO
          ENDDO

        ENDDO

      ENDIF nradl_block





!$omp parallel do                                                       &
!$omp& private(i,iendx,j,k)
      DO K=KTS,KTE
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IENDX=min(ide-(1),ite+(0))
          IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
          DO I=max(ids+(1),its-(0)),IENDX
            CLDFRA(I,J,K)=CLFR(I,K,J)
          ENDDO
        ENDDO
      ENDDO






!$omp parallel do                                                       &
!$omp& private(i,iendx,j)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        IENDX=min(ide-(1),ite+(0))
        IF(MOD(J,2)==0.AND.ITE==IDE)IENDX=IENDX-1
        DO I=max(ids+(1),its-(0)),IENDX
          HBOT(I,J)=HBOTR(I,J)
          HTOP(I,J)=HTOPR(I,J)
          CUPPT(I,J)=CUPPTR(I,J)
        ENDDO
      ENDDO





      DO J=JTS,JTE
      DO I=ITS,ITE
        IF(HBM2(I,J)<0.5)THEN
          ACFRST(I,J)=0.
          ACFRCV(I,J)=0.
          CFRACL(I,J)=0.
          CFRACM(I,J)=0.
          CFRACH(I,J)=0.
          RSWTOA(I,J)=0.
          RLWTOA(I,J)=0.
        ENDIF
      ENDDO
      ENDDO






      DO N=2,N_MOIST
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO J=JMS,JME
        DO K=KMS,KME
        DO I=IMS,IME
          MOIST(I,J,K,N)=MOIST_TRANS(I,K,J,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(MOIST_TRANS,STAT=ISTAT)



      END SUBROUTINE RADIATION




      SUBROUTINE TURBL(NTSD,DT,NPHS,RESTRT                              &
     &                ,N_MOIST,NSOIL,SLDPTH,DZSOIL                      &
     &                ,DETA1,DETA2,AETA1,AETA2,ETA1,ETA2,PDTOP,PT       &
     &                ,SM,HBM2,VBM2,DX_ARRAY,DFRLG                      &
     &                ,CZEN,CZMEAN,SIGT4,RLWIN,RSWIN,RADOT              &

     &                ,PD,RES,PINT,T,Q,CWM,F_ICE,F_RAIN,SR              &
     &                ,Q2,U,V,THS,TSFC,SST,PREC,SNO                     &
     &                ,FIS,Z0,Z0BASE,USTAR,MIXHT,PBLH,LPBL,EL_MYJ       &
     &                ,MOIST,RMOL,MOL                                   &
     &                ,EXCH_H,EXCH_M,F,AKHS,AKMS,AKHS_OUT,AKMS_OUT      &
     &                ,THZ0,QZ0,UZ0,VZ0,QS,MAVAIL                       &
     &                ,STC,SMC,CMC,SMSTAV,SMSTOT,SSROFF,BGROFF          &
     &                ,IVGTYP,ISLTYP,VEGFRC,SHDMIN,SHDMAX,GRNFLX        &
     &                ,SNOTIME                                          &
     &                ,SFCEXC,ACSNOW,ACSNOM,SNOPCX,SICE,TG,SOILTB       &
     &                ,ALBASE,MXSNAL,ALBEDO,SH2O,SI,EPSR,EMBCK          &
     &                ,U10,V10,TH10,Q10,TSHLTR,QSHLTR,PSHLTR            &
     &                ,T2,QSG,QVG,QCG,SOILT1,TSNAV,SMFR3D,KEEPFR3DFLAG  &
     &                ,TWBS,QWBS,TAUX,TAUY,SFCSHX,SFCLHX,SFCEVP                   &
     &                ,POTEVP,POTFLX,SUBSHX                             &
     &                ,APHTIM,ARDSW,ARDLW,ASRFC                         &
     &                ,RSWOUT,RSWTOA,RLWTOA                             &
     &                ,ASWIN,ASWOUT,ASWTOA,ALWIN,ALWOUT,ALWTOA          &
     &                ,UZ0H,VZ0H,DUDT,DVDT,UGWDsfc,VGWDsfc,SFENTH              & 
     &                ,RTHBLTEN,RQVBLTEN                                &
     &                ,PCPFLG,DDATA                                     & 
     &                ,HSTDV,HCNVX,HASYW,HASYS,HASYSW,HASYNW,HLENW      & 
     &                ,HLENS,HLENSW,HLENNW,HANGL,HANIS,HSLOP,HZMAX      & 
     &                ,CROT,SROT                                        & 
     &                ,GRID,CONFIG_FLAGS                                &
     &                ,IHE,IHW,IVE,IVW                                  &
     &                ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &                ,IMS,IME,JMS,JME,KMS,KME                          &
     &                ,ITS,ITE,JTS,JTE,KTS,KTE)


































      IMPLICIT NONE



      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,N_MOIST,NPHS,NSOIL,NTSD

      INTEGER, DIMENSION(JMS:JME),INTENT(IN) :: IHE,IHW,IVE,IVW

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ISLTYP,IVGTYP

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: LPBL

      REAL,INTENT(IN) :: DT,PDTOP,PT

      REAL,INTENT(IN) :: SFENTH
      REAL,INTENT(INOUT) :: APHTIM,ARDSW,ARDLW,ASRFC

      REAL,DIMENSION(KMS:KME-1),INTENT(IN) :: AETA1,AETA2,DETA1,DETA2

      REAL,DIMENSION(KMS:KME),INTENT(IN) :: DFRLG,ETA1,ETA2

      REAL,DIMENSION(NSOIL),INTENT(IN) :: DZSOIL,SLDPTH

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CZEN,CZMEAN         &
     &                                             ,DX_ARRAY            &
     &                                             ,F,FIS,HBM2          &
     &                                             ,PD,RES              &
     &                                             ,RLWIN,RLWTOA        &
     &                                             ,RSWIN,RSWOUT,RSWTOA &
     &                                             ,SHDMIN,SHDMAX       &

     &                                             ,SIGT4               &
     &                ,HSTDV,HCNVX,HASYW,HASYS,HASYSW,HASYNW,HLENW      & 
     &                ,HLENS,HLENSW,HLENNW,HANGL,HANIS,HSLOP,HZMAX      & 
     &                ,CROT,SROT                                        & 
     &                                             ,VBM2,VEGFRC
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: SST

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: SM,EPSR,SR       & 
                                                      ,TG,SICE          &
                                                      ,EMBCK
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ALBASE,MXSNAL

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACSNOM,ACSNOW    &
     &                                                ,SNOTIME          &
     &                                                ,AKHS,AKMS        &
     &                                                ,ALBEDO           &
     &                                                ,BGROFF,CMC       &
     &                                                ,MAVAIL,MOL       &
     &                                                ,MIXHT            &
     &                                                ,PBLH,POTEVP      &
     &                                                ,POTFLX,PREC      &
     &                                                ,QCG,QS,QSG       &
     &                                                ,QVG,QZ0          &
     &                                                ,RMOL             &
     &                                                ,SFCEVP           &
     &                                                ,SFCLHX,SFCSHX    &
     &                                                ,SI,SMSTOT        &
     &                                                ,SNO,SNOPCX       &
     &                                                ,SOILT1           &
     &                                                ,SSROFF,SUBSHX    &
     &                                                ,T2,THS,THZ0      &
     &                                                ,TSFC,TSNAV       &
     &                                                ,USTAR,UZ0,UZ0H   &
     &                                                ,VZ0,VZ0H         &
     &                                                ,Z0,Z0BASE

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: AKHS_OUT,AKMS_OUT  &
     &                                              ,ALWIN,ALWOUT       &
     &                                              ,ALWTOA,ASWIN       &
     &                                              ,ASWOUT,ASWTOA      &
     &                                              ,PSHLTR,Q10,QSHLTR  &
     &                                              ,TH10,TSHLTR        &
     &                                              ,U10,V10            & 
     &                                              ,UGWDsfc,VGWDsfc      

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: GRNFLX,QWBS,RADOT  &
                                                    ,SFCEXC,SMSTAV      &
                                                    ,SOILTB,TWBS


      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: taux, tauy
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: PINT

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: CWM      &
     &                                                        ,DUDT     &
     &                                                        ,DVDT     &
     &                                                        ,Q,Q2     &
     &                                                        ,T,U,V

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) :: F_ICE    &   
     &                                                        ,F_RAIN   &
     &                                                        ,RQVBLTEN &
     &                                                        ,RTHBLTEN

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(OUT) :: EL_MYJ     &   
     &                                                      ,EXCH_H     &
     &                                                      ,EXCH_M

      REAL,DIMENSION(IMS:IME,NSOIL,JMS:JME),INTENT(INOUT) :: KEEPFR3DFLAG & 
     &                                                      ,SH2O,SMC     &
     &                                                      ,SMFR3D,STC

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,N_MOIST)                   &
     &                                           ,INTENT(INOUT) :: MOIST

      LOGICAL,INTENT(IN) :: RESTRT

      TYPE(DOMAIN),TARGET :: GRID

      TYPE(GRID_CONFIG_REC_TYPE),INTENT(IN) :: CONFIG_FLAGS


      LOGICAL,INTENT(IN) :: PCPFLG
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: DDATA






      INTEGER :: I,I_M,IDUMMY,IEND,ISFFLX,ISTAT,ISTR,J,K,KOUNT_ALL      &
     &          ,LENGTH_ROW,LLIJ,LLYR,N,SST_UPDATE,SF_URBAN_PHYSICS,NUM_URBAN_LAYERS

      INTEGER,DIMENSION(IMS:IME,JMS:JME) :: KPBL,LOWLYR

      REAL :: TRESH=0.95

      REAL :: ALTITUDE,CWML,DQDT,DTDT,DTPHS,DX,DZHALF,FACTR,FACTRL      &
     &       ,G_INV,PLYR,PSFC,QI,QL,QOLD,QR,QW,RATIOMX,RDTPHS      &
     &       ,ROG,RWMSK,SDEPTH,SNO_FACTR,TL,TLMH,TLMH4,TNEW,TSFC2       &
     &       ,U_FRAME,V_FRAME,XLVRW

      REAL :: APES,CAPA,CKLQ,EXNER,FACTOR,FFS,PQ0X,Q2SAT,QFC1,QLOWX     &
     &       ,RLIVWV,THBOT,DPL

      REAL,DIMENSION(IMS:IME,JMS:JME) :: BR,CHKLOWQ,CT,CWMLOW,ELFLX     &
     &                                  ,EXNSFC,FACTRS,FLHC,FLQC,GZ1OZ0 &
     &                                  ,ONE,PDSL,PLM,PSFC_OUT,PSIH     &
     &                                  ,PSIM,Q2X,QLOW,RAIN,RAINBL      &
     &                                  ,RLW_DN_SFC,RSW_NET_SFC         &
     &                                  ,RSW_DN_SFC                     &
     &                                  ,SFCEVPX,SFCZ,SNOW,SNOWC,SNOWH  &
     &                                  ,TH2X,THLOW,TLOW,VGFRCK         &
     &                                  ,WSPD,XLAND,REGIME,HOL

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: DUDT_PHY,DVDT_PHY,DZ   &
     &                                          ,P_PHY,P8W,PI_PHY       &
     &                                          ,RQCBLTEN,RQIBLTEN      &

     &                                          ,RR,DELP                & 
     &                                          ,T_PHY,TH_PHY,TKE       &
     &                                          ,DUDT_GWD,DVDT_GWD      & 
     &                                          ,U_PHY,V_PHY,Z

      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: MOIST_TRANS

      REAL,DIMENSION(IMS:IME,NSOIL,JMS:JME) :: ZERO_SOIL

      LOGICAL :: E_BDY,WARM_RAIN

      INTEGER :: NUM_ROOF_LAYERS,NUM_WALL_LAYERS,NUM_ROAD_LAYERS   
      INTEGER :: FRACTIONAL_SEAICE

      INTEGER :: IGS,IGE,JGS,JGE, PQ_I   
      LOGICAL :: FQ_I                    

      CHARACTER(len=255) :: message



    INTEGER  :: isurban
    CHARACTER(len=256) :: MMINLU




      ALLOCATE(MOIST_TRANS(IMS:IME,KMS:KME,JMS:JME,N_MOIST),STAT=ISTAT)

      SF_URBAN_PHYSICS=CONFIG_FLAGS%SF_URBAN_PHYSICS

      if ( config_flags%bl_pbl_physics == BOULACSCHEME ) then
         call wrf_error_fatal3("",832,&
"Cannot use BOULAC PBL with NMM")
      endif
      
      FRACTIONAL_SEAICE = CONFIG_FLAGS%FRACTIONAL_SEAICE
      IF ( FRACTIONAL_SEAICE == 1 ) THEN
         CALL wrf_error_fatal3("",838,&
"NMM cannot use FRACTIONAL_SEAICE = 1.")
      ENDIF

      DTPHS=NPHS*DT
      RDTPHS=1./DTPHS
      G_INV=1./G
      ROG=R_D*G_INV
      FACTOR=-XLV*RHOWATER/DTPHS
      CAPA=R_D/CP

      U_FRAME=0.
      V_FRAME=0.

      IDUMMY=0
      ISFFLX=1
      DX=0.
      SST_UPDATE=0

!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=JMS,JME
      DO I=IMS,IME
        UZ0H(I,J)=0.
        VZ0H(I,J)=0.
        ONE(I,J)=1.
        RMOL(I,J)=0.     
        SFCEVPX(I,J)=0.  
      ENDDO
      ENDDO

      IF(MODEL_CONFIG_REC%SF_SURFACE_PHYSICS(GRID%ID)==99)THEN
        SNO_FACTR=1.
      ELSE
        SNO_FACTR=0.001
      ENDIF

!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
      DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
        LOWLYR(I,J)=1
        VGFRCK(I,J)=100.*VEGFRC(I,J)
        SNOW(I,J)=SNO(I,J)
        SNOWH(I,J)=SI(I,J)*SNO_FACTR
        XLAND(I,J)=SM(I,J)+1.
        T2(I,J)=TSFC(I,J)
      ENDDO
      ENDDO

      IF(NTSD==0)THEN
!$omp parallel do                                                       &
!$omp& private(i,j)
        DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
        DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
          Z0BASE(I,J)=Z0(I,J)
          IF(SM(I,J)>0.5.AND.SICE(I,J)>0.5)THEN  
            SM(I,J)=0.
          ENDIF
        ENDDO
        ENDDO
      ENDIF

!$omp parallel do                                                       &
!$omp& private(i,j,k)
      DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
      DO K=KTS,KTE+1
      DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
        Z(I,K,J)=0.
        DZ(I,K,J)=0.
        EXCH_H(I,K,J)=0.
        EXCH_M(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO







!$omp parallel do                                                       &
!$omp& private(factrl,i,j,llij,tlmh)
      DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
      DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))

        PDSL(I,J)=PD(I,J)*RES(I,J)


        P8W(I,KTS,J)=PINT(I,J,KTS)
        PSFC=PINT(I,J,KTS)
        LOWLYR(I,J)=KTS     
        EXNSFC(I,J)=(1.E5/PSFC)**CAPA
        THS(I,J)=(SST(I,J)*EXNSFC(I,J))*SM(I,J)+THS(I,J)*(1.-SM(I,J))
        TSFC(I,J)=THS(I,J)/EXNSFC(I,J)
        SFCZ(I,J)=FIS(I,J)*G_INV

        IF (PCPFLG.AND.DDATA(I,J)<100.)THEN
          RAIN(I,J)=DDATA(I,J)*RHOWATER
        ELSE
          RAIN(I,J)=PREC(I,J)*RHOWATER
        ENDIF

        RAINBL(I,J)=0.
        IF(SNO(I,J)>0.)SNOWC(I,J)=1.
        LLIJ=LOWLYR(I,J)
        PLM(I,J)=(PINT(I,J,LLIJ)+PINT(I,J,LLIJ+1))*0.5
        TH2X(I,J)=T(I,J,LLIJ)*(1.E5/PLM(I,J))**CAPA
        Q2X(I,J)=Q(I,J,LLIJ)





        IF(CZMEAN(I,J)>0.)THEN
          FACTRS(I,J)=CZEN(I,J)/CZMEAN(I,J)
        ELSE
          FACTRS(I,J)=0.
        ENDIF

        IF(SIGT4(I,J)>0.)THEN
          TLMH=T(I,J,LLIJ)
          FACTRL=STBOLT*TLMH*TLMH*TLMH*TLMH/SIGT4(I,J)
        ELSE
          FACTRL=0.
        ENDIF



        RLW_DN_SFC(I,J)=RLWIN(I,J)*HBM2(I,J)*FACTRL
        RSW_NET_SFC(I,J)=(RSWIN(I,J)-RSWOUT(I,J))*HBM2(I,J)*FACTRS(I,J)



        RSW_DN_SFC(I,J)=RSWIN(I,J)*HBM2(I,J)*FACTRS(I,J)

        Z(I,KTS,J)=SFCZ(I,J)

      ENDDO
      ENDDO





!$omp parallel do                                                       &
!$omp& private(cwml,i,j,k,plyr,qi,ql,qr,qw,tl)
      DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
        DO K=KTS,KTE
        DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
          Q2(I,J,K)=MAX(Q2(I,J,K)*HBM2(I,J),EPSQ2)
          QL=MAX(Q(I,J,K),EPSQ)
          PLYR=(PINT(I,J,K)+PINT(I,J,K+1))*0.5

          TL=T(I,J,K)
          CWML=CWM(I,J,K)

          RR(I,K,J)=PLYR/(R_D*TL)
          T_PHY(I,K,J)=TL

          EXNER=(1.E5/PLYR)**CAPA
          PI_PHY(I,K,J)=1./EXNER
          TH_PHY(I,K,J)=TL*EXNER
          P8W(I,K+1,J)=PINT(I,J,K+1)

          P_PHY(I,K,J)=PLYR
          TKE(I,K,J)=0.5*Q2(I,J,K)

          RTHBLTEN(I,K,J)=0.
          RQVBLTEN(I,K,J)=0.
          RQCBLTEN(I,K,J)=0.
          RQIBLTEN(I,K,J)=0.






          DPL=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
          Z(I,K+1,J)=Z(I,K,J)+TL/PLYR*DPL*ROG*(Q(I,J,K)*P608-CWML+1.)
          DELP(I,K,J)=DPL
          DZ(I,K,J)=Z(I,K+1,J)-Z(I,K,J)
        ENDDO
      ENDDO
      ENDDO

!$omp parallel do                                                       &
!$omp& private(i,j,llyr,qlowx)
      DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
      DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
        TWBS(I,J)=0.
        QWBS(I,J)=0.
        LLYR=LOWLYR(I,J)
        THLOW(I,J)=TH_PHY(I,LLYR,J)
        TLOW(I,J)=T_PHY(I,LLYR,J)
        QLOW(I,J)=MAX(Q(I,J,LLYR),EPSQ)
        QLOWX=QLOW(I,J)/(1.-QLOW(I,J))
        QLOW(I,J)=QLOWX/(1.+QLOWX)
        CWMLOW(I,J)=CWM(I,J,LLYR)
        PBLH(I,J)=MAX(PBLH(I,J),0.)
        PBLH(I,J)=MIN(PBLH(I,J),Z(I,KTE,J))
      ENDDO
      ENDDO





!$omp parallel do                                                       &
!$omp& private(i,j,k)
      DO K=KTS,KTE
        DO J=max(jds+(1),jts-(1)),min(jde-(1),jte+(1))
          DO I=max(ids+(0),its-(1)),min(ide-(0),ite+(1))
            U_PHY(I,K,J)=(U(I+IHE(J),J,K)+U(I+IHW(J),J,K)               &
     &                   +U(I,J+1,K)+U(I,J-1,K))                        &
     &                   *0.25
            V_PHY(I,K,J)=(V(I+IHE(J),J,K)+V(I+IHW(J),J,K)               &
     &                   +V(I,J+1,K)+V(I,J-1,K))                        &
     &                   *0.25
          ENDDO
        ENDDO
      ENDDO

!$omp parallel do                                                       &
!$omp& private(i,iend,istr,j)
      DO J=max(jds+(1),jts-(1)),min(jde-(1),jte+(1))
        IF(MOD(J,2)==0)THEN
          ISTR=max(ids+(0),its-(1))
          IEND=MIN(min(ide-(0),ite+(1)),IDE-1)
        ELSE
          ISTR=MAX(max(ids+(0),its-(1)),IDS+1)
          IEND=MIN(min(ide-(0),ite+(1)),IDE-1)
        ENDIF

        DO I=ISTR,IEND
          UZ0H(I,J)=(UZ0(I+IHE(J),J)+UZ0(I+IHW(J),J)                    &
     &              +UZ0(I,J+1)+UZ0(I,J-1))*0.25

          VZ0H(I,J)=(VZ0(I+IHE(J),J)+VZ0(I+IHW(J),J)                    &
     &              +VZ0(I,J+1)+VZ0(I,J-1))*0.25

        ENDDO
      ENDDO





      DO J=JTS,JTE
      DO I=ITS,ITE
        IF(MODEL_CONFIG_REC%SF_SURFACE_PHYSICS(GRID%ID)==2.OR.          &
           MODEL_CONFIG_REC%SF_SURFACE_PHYSICS(GRID%ID)==99)THEN
          ONE(I,J)=1.
        ELSE

          ONE(I,J)=MAVAIL(I,J)
        ENDIF
      ENDDO
      ENDDO





      DO N=1,N_MOIST
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO K=KMS,KME
        DO J=JMS,JME
        DO I=IMS,IME
          MOIST_TRANS(I,K,J,N)=MOIST(I,J,K,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO





      NUM_ROOF_LAYERS=GRID%NUM_SOIL_LAYERS   
      NUM_WALL_LAYERS=GRID%NUM_SOIL_LAYERS   
      NUM_ROAD_LAYERS=GRID%NUM_SOIL_LAYERS   
      CALL nl_get_isurban(grid%id, isurban)
      call nl_get_mminlu(grid%id, mminlu)







      CALL SET_TILES(GRID,IDS,IDE-1,JDS+1,JDE-1,ITS,ITE,JTS,JTE)

CALL SURFACE_DRIVER(ACSNOM=ACSNOM,ACSNOW=ACSNOW,AKHS=AKHS,AKMS=AKMS,ALBEDO=ALBEDO,BR=BR,CANWAT=CMC,CHKLOWQ=CHKLOWQ,DT=DT,DX=DX,DZ8&
&W=DZ,DZS=DZSOIL,GLW=RLW_DN_SFC,GRDFLX=GRNFLX,GSW=RSW_NET_SFC,SWDOWN=RSW_DN_SFC,GZ1OZ0=GZ1OZ0,HFX=TWBS,HT=SFCZ,IFSNOW=IDUMMY,ISFFLX&
&=ISFFLX,FRACTIONAL_SEAICE=FRACTIONAL_SEAICE,ISLTYP=ISLTYP,ITIMESTEP=NTSD,IVGTYP=IVGTYP,LOWLYR=LOWLYR,MAVAIL=ONE,RMOL=RMOL,MOL=MOL,&
&NUM_SOIL_LAYERS=NSOIL,P8W=P8W,PBLH=PBLH,PI_PHY=PI_PHY,PSHLTR=PSHLTR,PSIH=PSIH,PSIM=PSIM,P_PHY=P_PHY,Q10=Q10,Q2=Q2X,QFX=QWBS,TAUX=T&
&AUX,TAUY=TAUY,QSFC=QS,QSHLTR=QSHLTR,QZ0=QZ0,RAINCV=RAIN,RHO=RR,SFCEVP=SFCEVPX,SFCEXC=SFCEXC,SFCRUNOFF=SSROFF,SMOIS=SMC,SMSTAV=SMST&
&AV,SMSTOT=SMSTOT,SNOALB=MXSNAL,SNOW=SNOW,SNOWC=SNOWC,SNOWH=SNOWH,STEPBL=NPHS,SST=SST,SST_UPDATE=SST_UPDATE,TH10=TH10,TH2=TH2X,T2=T&
&2,THZ0=THZ0,TH_PHY=TH_PHY,TMN=TG,TSHLTR=TSHLTR,TSK=TSFC,TSLB=STC,T_PHY=T_PHY,U10=U10,UDRUNOFF=BGROFF,UST=USTAR,UZ0=UZ0H,U_FRAME=U_&
&FRAME,U_PHY=U_PHY,V10=V10,VEGFRA=VGFRCK,VZ0=VZ0H,V_FRAME=V_FRAME,V_PHY=V_PHY,WARM_RAIN=WARM_RAIN,WSPD=WSPD,XICE=SICE,XICEM=SICE,IS&
&ICE=GRID%LANDUSE_ISICE,ISWATER=GRID%ISWATER,XLAND=XLAND,Z=Z,ZNT=Z0,ZS=SLDPTH,CT=CT,TKE_MYJ=TKE,SFENTH=SFENTH,ALBBCK=ALBASE,LH=ELFL&
&X,SH2O=SH2O,SHDMAX=SHDMAX,SHDMIN=SHDMIN,Z0=Z0BASE,FLQC=FLQC,FLHC=FLHC,PSFC=PSFC_OUT,EMISS=EPSR,EMBCK=EMBCK,SF_SFCLAY_PHYSICS=CONFI&
&G_FLAGS%SF_SFCLAY_PHYSICS,SF_SURFACE_PHYSICS=CONFIG_FLAGS%SF_SURFACE_PHYSICS,RA_LW_PHYSICS=CONFIG_FLAGS%RA_LW_PHYSICS,LAI=GRID%LAI&
&,SF_URBAN_PHYSICS=SF_URBAN_PHYSICS,NUM_URBAN_LAYERS=NUM_URBAN_LAYERS,IDS=IDS,IDE=IDE,JDS=JDS,JDE=JDE,KDS=KDS,KDE=KDE,IMS=IMS,IME=I&
&ME,JMS=JMS,JME=JME,KMS=KMS,KME=KME,I_START=GRID%I_START,I_END=GRID%I_END,J_START=GRID%J_START,J_END=GRID%J_END,KTS=KTS,KTE=KTE,NUM&
&_TILES=GRID%NUM_TILES,QV_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QV),F_QV=F_QV,QC_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QC),F_QC=F_QC,QR_CURR=MOIS&
&T_TRANS(IMS,KMS,JMS,P_QR),F_QR=F_QR,QI_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QI),F_QI=F_QI,QS_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QS),F_QS=F_Q&
&S,QG_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QG),F_QG=F_QG,RAINBL=RAINBL,LAGDAY=1,QSG=QSG,QVG=QVG,QCG=QCG,SOILT1=SOILT1,TSNAV=TSNAV,SMFR3D=&
&SMFR3D,KEEPFR3DFLAG=KEEPFR3DFLAG,POTEVP=POTEVP,SNOPCX=SNOPCX,SOILTB=SOILTB,SR=SR,NUM_ROOF_LAYERS=NUM_ROOF_LAYERS,NUM_WALL_LAYERS=N&
&UM_WALL_LAYERS,NUM_ROAD_LAYERS=NUM_ROAD_LAYERS,REGIME=REGIME,NLCAT=grid%num_land_cat,NSCAT=grid%num_soil_cat,ISURBAN=isurban,MMINL&
&U=TRIM(mminlu),SNOTIME=grid%SNOTIME,RDLAI2D=config_flags%rdlai2d,usemonalb=config_flags%usemonalb,NOAHRES=grid%noahres)







!$omp parallel do                                                       &
!$omp& private(i,j,k)
      DO J=JMS,JME
      DO K=KMS,KME
      DO I=IMS,IME
        DUDT_PHY(I,K,J)=0.
        DVDT_PHY(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO






!$omp parallel do                                                       &
!$omp& private(dzhalf,i,j)
      DO J=JTS,JTE
      DO I=ITS,ITE
        DZHALF=0.5*DZ(I,KTS,J)
        AKHS_OUT(I,J)=AKHS(I,J)*DZHALF
        AKMS_OUT(I,J)=AKMS(I,J)*DZHALF
      ENDDO
      ENDDO




      IF(CONFIG_FLAGS%MP_PHYSICS==ETAMPNEW)THEN
         FQ_I=F_QS
         PQ_I=P_QS
      ELSE
         FQ_I=F_QI
         PQ_I=P_QI
      ENDIF

CALL PBL_DRIVER(ITIMESTEP=NTSD,DT=DT,U_FRAME=U_FRAME,V_FRAME=V_FRAME,RUBLTEN=DUDT_PHY,RVBLTEN=DVDT_PHY,RTHBLTEN=RTHBLTEN,RQVBLTEN=&
&RQVBLTEN,RQCBLTEN=RQCBLTEN,RQIBLTEN=RQIBLTEN,TSK=TSFC,XLAND=XLAND,ZNT=Z0,HT=SFCZ,UST=USTAR,MIXHT=MIXHT,PBLH=PBLH,HFX=TWBS,QFX=QWBS&
&,GRDFLX=GRNFLX,U_PHY=U_PHY,V_PHY=V_PHY,TH_PHY=TH_PHY,RHO=RR,P_PHY=P_PHY,PI_PHY=PI_PHY,P8W=P8W,T_PHY=T_PHY,DZ8W=DZ,Z=Z,TKE_MYJ=TKE,&
&EL_MYJ=EL_MYJ,F=F,EXCH_H=EXCH_H,EXCH_M=EXCH_M,AKHS=AKHS,AKMS=AKMS,THZ0=THZ0,QZ0=QZ0,UZ0=UZ0H,VZ0=VZ0H,QSFC=QS,LOWLYR=LOWLYR,PSIM=P&
&SIM,PSIH=PSIH,GZ1OZ0=GZ1OZ0,U10=U10,V10=V10,WSPD=WSPD,BR=BR,CHKLOWQ=CHKLOWQ,DX=DX,STEPBL=NPHS,WARM_RAIN=WARM_RAIN,KPBL=KPBL,CT=CT,&
&LH=ELFLX,SNOW=SNOW,XICE=SICE,BL_PBL_PHYSICS=config_flags%bl_pbl_physics,RA_LW_PHYSICS=config_flags%ra_lw_physics,IDS=IDS,IDE=IDE,J&
&DS=JDS,JDE=JDE,KDS=KDS,KDE=KDE,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=KMS,KME=KME,I_START=GRID%I_START,I_END=GRID%I_END,J_START=GRID%&
&J_START,J_END=GRID%J_END,KTS=KTS,KTE=KTE,NUM_TILES=GRID%NUM_TILES,QV_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QV),F_QV=F_QV,QC_CURR=MOIST_TR&
&ANS(IMS,KMS,JMS,P_QC),F_QC=F_QC,QI_CURR=MOIST_TRANS(IMS,KMS,JMS,PQ_I),F_QI=FQ_I,QR_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QR),F_QR=F_QR,QS&
&_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QS),F_QS=F_QS,QG_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QG),F_QG=F_QG,HOL=HOL,sf_sfclay_physics=CONFIG_FLA&
&GS%SF_SFCLAY_PHYSICS,sf_urban_physics=CONFIG_FLAGS%SF_URBAN_PHYSICS)











      IF(min(ide-(0),ite+(0))==IDE)THEN
!$omp parallel do                                                       &
!$omp& private(i,j)
        DO J=JDS,JDE
        IF (J>=max(jds+(0),jts-(0)).AND.J<=min(jde-(0),jte+(0)))THEN
          TH10(min(ide-(0),ite+(0)),J)=TH10(min(ide-(0),ite+(0))-1,J)
          Q10(min(ide-(0),ite+(0)),J)=Q10(min(ide-(0),ite+(0))-1,J)
          U10(min(ide-(0),ite+(0)),J)=U10(min(ide-(0),ite+(0))-1,J)
          V10(min(ide-(0),ite+(0)),J)=V10(min(ide-(0),ite+(0))-1,J)
          TSHLTR(min(ide-(0),ite+(0)),J)=TSHLTR(min(ide-(0),ite+(0))-1,J)
          QSHLTR(min(ide-(0),ite+(0)),J)=QSHLTR(min(ide-(0),ite+(0))-1,J)
        ENDIF
        ENDDO
      ENDIF




      IF(max(jds+(0),jts-(0))==JDS)THEN
        DO J=JDS,JDS+1
        DO I=IDS,IDE
          IF (I>=max(ids+(0),its-(0)).AND.I<=min(ide-(0),ite+(0))) THEN
            TH10(I,J)=TH10(I,max(jds+(0),jts-(0))+2)
            Q10(I,J)=Q10(I,max(jds+(0),jts-(0))+2)
            U10(I,J)=U10(I,max(jds+(0),jts-(0))+2)
            V10(I,J)=V10(I,max(jds+(0),jts-(0))+2)
            TSHLTR(I,J)=TSHLTR(I,max(jds+(0),jts-(0))+2)
            QSHLTR(I,J)=QSHLTR(I,max(jds+(0),jts-(0))+2)
          ENDIF
        ENDDO
        ENDDO
      ENDIF



      IF(min(jde-(0),jte+(0))==JDE)THEN
!$omp parallel do                                                       &
!$omp& private(i,j)
        DO J=min(jde-(0),jte+(0))-1,min(jde-(0),jte+(0))
        DO I=IDS,IDE
          IF (I>=max(ids+(0),its-(0)).AND.I<=min(ide-(0),ite+(0))) THEN
            TH10(I,J)=TH10(I,min(jde-(0),jte+(0))-2)
            Q10(I,J)=Q10(I,min(jde-(0),jte+(0))-2)
            U10(I,J)=U10(I,min(jde-(0),jte+(0))-2)
            V10(I,J)=V10(I,min(jde-(0),jte+(0))-2)
            TSHLTR(I,J)=TSHLTR(I,min(jde-(0),jte+(0))-2)
            QSHLTR(I,J)=QSHLTR(I,min(jde-(0),jte+(0))-2)
          ENDIF
        ENDDO
        ENDDO
      ENDIF

      IF(CONFIG_FLAGS%SF_SFCLAY_PHYSICS==1)THEN 
!$omp parallel do                                                       &
!$omp& private(i,j)
        DO J=max(jds+(1),jts-(0)),min(jde-(1),jte+(0))
        DO I=max(ids+(0),its-(0)),min(ide-(1),ite+(0))

          IF(TSHLTR(I,J)<200..OR.TSHLTR(I,J)>350.)THEN
            WRITE(message,*)'Troublesome TSHLTR...I,J,TSHLTR,PSHLTR: '        &
                      ,I,J,TSHLTR(I,J),PSHLTR(I,J)
            CALL wrf_message(trim(message))
          ENDIF
	ENDDO
	ENDDO
      ENDIF





      IF(CONFIG_FLAGS%BL_PBL_PHYSICS/=MYJPBLSCHEME)THEN
        LENGTH_ROW=min(ide-(1),ite+(0))-max(ids+(1),its-(0))+1
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
          KPBL(I,J)=-1000
        ENDDO
        ENDDO

!$omp parallel do                                                       &
!$omp& private(altitude,i,j,k,kount_all)
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          KOUNT_ALL=0
          find_kpbl : DO K=KTS,KTE
          DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
            ALTITUDE=Z(I,K+1,J)-SFCZ(I,J)
            IF(PBLH(I,J)<=ALTITUDE.AND.KPBL(I,J)<0)THEN
              KPBL(I,J)=K
              KOUNT_ALL=KOUNT_ALL+1
            ENDIF
            IF(KOUNT_ALL==LENGTH_ROW)EXIT find_kpbl
          ENDDO
          ENDDO find_kpbl
        ENDDO
      ENDIF

      IF(MODEL_CONFIG_REC%SF_SURFACE_PHYSICS(GRID%ID)==99)THEN
        SNO_FACTR=1.
      ELSE
        SNO_FACTR=1000.
      ENDIF

!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
      DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
        SNO(I,J)=SNOW(I,J)
        SI(I,J)=SNOWH(I,J)*SNO_FACTR
        LPBL(I,J)=KTE-KPBL(I,J)+1
      ENDDO
      ENDDO





!$omp parallel do                                                       &
!$omp& private(i,j,tsfc2)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
      DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
        ASWIN (I,J)=ASWIN (I,J)+RSWIN(I,J)*HBM2(I,J)*FACTRS(I,J)
        ASWOUT(I,J)=ASWOUT(I,J)-RSWOUT(I,J)*HBM2(I,J)*FACTRS(I,J)
        ASWTOA(I,J)=ASWTOA(I,J)+RSWTOA(I,J)*HBM2(I,J)*FACTRS(I,J)
        ALWIN (I,J)=ALWIN (I,J)+RLW_DN_SFC(I,J)
        ALWOUT(I,J)=ALWOUT(I,J)-RADOT (I,J)*HBM2(I,J)
        ALWTOA(I,J)=ALWTOA(I,J)+RLWTOA(I,J)*HBM2(I,J)

        TSFC2=TSFC(I,J)*TSFC(I,J)
        RADOT(I,J)=HBM2(I,J)*EPSR(I,J)*STBOLT*TSFC2*TSFC2
        THS(I,J)=TSFC(I,J)*EXNSFC(I,J)
        PREC(I,J)=0.
      ENDDO
      ENDDO





      IGS=max(ids+(1),its-(0))
      IGE=min(ide-(1),ite+(0))
      JGS=max(jds+(2),jts-(0))
      JGE=min(jde-(2),jte+(0))












      IF (grid%gwd_opt .eq. 2) THEN
       
        CALL wrf_message("GWD usage currently may be problematic for some cases - use at own risk")

      CALL GWD_driver(U=U_PHY,V=V_PHY,T=T_PHY                           &
     &               ,Q=MOIST_TRANS(IMS,KMS,JMS,P_QV)                   &
     &               ,Z=Z,DP=DELP,PINT=P8W,PMID=P_PHY,EXNR=PI_PHY       &
     &               ,KPBL=KPBL,ITIME=NTSD                              &
     &               ,HSTDV=HSTDV,HCNVX=HCNVX,HASYW=HASYW,HASYS=HASYS   &
     &               ,HASYSW=HASYSW,HASYNW=HASYNW,HLENW=HLENW           &
     &               ,HLENS=HLENS,HLENSW=HLENSW,HLENNW=HLENNW           &
     &               ,HANGL=HANGL,HANIS=HANIS,HSLOP=HSLOP,HZMAX=HZMAX   &
     &               ,CROT=CROT,SROT=SROT                               &
     &               ,DUDT=DUDT_GWD,DVDT=DVDT_GWD                       &
     &               ,UGWDsfc=UGWDsfc,VGWDsfc=VGWDsfc                   &
     &               ,IDS=IDS,IDE=IDE,JDS=JDS,JDE=JDE,KDS=KDS,KDE=KDE   &
     &               ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=KMS,KME=KME   &
     &               ,ITS=IGS,ITE=IGE,JTS=JGS,JTE=JGE,KTS=KTS,KTE=KTE )









      DO K=KTS,KTE
      DO J=JTS,JTE
      DO I=ITS,ITE




        IF (DUDT_GWD(I,K,J) .gt. 1.6) then
        write(message,*) 'BIG DUDT_GWD:: ', I,K,J, DUDT_GWD(I,K,J)
        CALL wrf_message(message)
        DUDT_GWD(I,K,J)=1.6
        ENDIF

        IF (DUDT_GWD(I,K,J) .lt. -1.6) then
        write(message,*) 'BIG DUDT_GWD:: ', I,K,J, DUDT_GWD(I,K,J)
        CALL wrf_message(message)
        DUDT_GWD(I,K,J)=-1.6
        ENDIF

        IF (DVDT_GWD(I,K,J) .gt. 1.6) then
        write(message,*) 'BIG DVDT_GWD:: ', I,K,J, DVDT_GWD(I,K,J)
        CALL wrf_message(message)
        DVDT_GWD(I,K,J)=1.6
        ENDIF

        IF (DVDT_GWD(I,K,J) .lt. -1.6) then
        write(message,*) 'BIG DVDT_GWD:: ', I,K,J, DVDT_GWD(I,K,J)
        CALL wrf_message(message)
        DVDT_GWD(I,K,J)=-1.6
        ENDIF



        DUDT(I,J,K)=DUDT_PHY(I,K,J)+DUDT_GWD(I,K,J)
        DVDT(I,J,K)=DVDT_PHY(I,K,J)+DVDT_GWD(I,K,J)

      ENDDO
      ENDDO
      ENDDO

      ELSE  

      DO K=KTS,KTE
      DO J=JTS,JTE
      DO I=ITS,ITE
        DUDT(I,J,K)=DUDT_PHY(I,K,J)
        DVDT(I,J,K)=DVDT_PHY(I,K,J)
      ENDDO
      ENDDO
      ENDDO

      ENDIF 





      DO N=1,N_MOIST
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO J=JMS,JME
        DO K=KMS,KME
        DO I=IMS,IME
          MOIST(I,J,K,N)=MOIST_TRANS(I,K,J,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(MOIST_TRANS,STAT=ISTAT)





      E_BDY=(ITE>=IDE)

!$omp parallel do                                                       &
!$omp& private(dqdt,dtdt,i,iend,j,k,qi,qold,qr,qw,ratiomx,i_m)
      DO K=KTS,KTE
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        IEND=min(ide-(1),ite+(0))
        IF(E_BDY.AND.MOD(J,2)==0)IEND=IEND-1

        DO I=max(ids+(1),its-(0)),IEND
          DTDT=RTHBLTEN(I,K,J)*PI_PHY(I,K,J)
          DQDT=RQVBLTEN(I,K,J)         
          T(I,J,K)=T(I,J,K)+DTDT*DTPHS
          QOLD=Q(I,J,K)
          RATIOMX=QOLD/(1.-QOLD)+DQDT*DTPHS
          Q(I,J,K)=RATIOMX/(1.+RATIOMX)

          MOIST(I,J,K,P_QV)=MAX(EPSQ,(MOIST(I,J,K,P_QV)+RQVBLTEN(I,K,J)*DTPHS) )
          CWM(I,J,K)=0.

          IF(CONFIG_FLAGS%MP_PHYSICS/=ETAMPNEW)THEN
            DO I_M=1,N_MOIST
              IF(I_M==P_QC) THEN
                MOIST(I,J,K,I_M)=MAX(EPSQ,(MOIST(I,J,K,I_M)+RQCBLTEN(I,K,J)*DTPHS) )
              ELSE IF(I_M==P_QI) THEN
                MOIST(I,J,K,I_M)=MAX(EPSQ,(MOIST(I,J,K,I_M)+RQIBLTEN(I,K,J)*DTPHS) )






              ENDIF
              IF(I_M/=P_QV) CWM(I,J,K)=CWM(I,J,K)+MOIST(I,J,K,I_M)
            ENDDO
          ELSE



            QW=MAX(0.,MOIST(I,J,K,P_QC)+RQCBLTEN(I,K,J)*DTPHS )
            QI=MAX(0.,MOIST(I,J,K,P_QS)+RQIBLTEN(I,K,J)*DTPHS )  

            QR=MAX(0.,MOIST(I,J,K,P_QR) )
            MOIST(I,J,K,P_QC)=QW
            MOIST(I,J,K,P_QS)=QI
            MOIST(I,J,K,P_QR)=QR
            CWM(I,J,K)=QW+QI+QR

            IF(QI<=EPSQ)THEN
              F_ICE(I,K,J)=0.
            ELSE
              F_ICE(I,K,J)=MAX(0.,MIN(1.,QI/CWM(I,J,K)))
            ENDIF

            IF(QR<=EPSQ)THEN
              F_RAIN(I,K,J)=0.
            ELSE
              F_RAIN(I,K,J)=QR/(QW+QR)
            ENDIF

          ENDIF     

          Q2(I,J,K)=2.*TKE(I,K,J)
        ENDDO
        ENDDO

      ENDDO






!$omp parallel do                                                       &
!$omp& private(i,j,llij,xlvrw)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
      DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
        LLIJ=LOWLYR(I,J)





        TWBS(I,J)=-TWBS(I,J)
        QWBS(I,J)=-QWBS(I,J)*XLV*CHKLOWQ(I,J)








        SFCSHX(I,J)=SFCSHX(I,J)+TWBS(I,J)
        SFCLHX(I,J)=SFCLHX(I,J)+QWBS(I,J)
        XLVRW=DTPHS/(XLV*RHOWATER)
        SFCEVP(I,J)=SFCEVP(I,J)-QWBS(I,J)*XLVRW
        POTEVP(I,J)=POTEVP(I,J)-QWBS(I,J)*SM(I,J)*XLVRW
        POTFLX(I,J)=POTEVP(I,J)*FACTOR
        SUBSHX(I,J)=SUBSHX(I,J)+GRNFLX(I,J)
      ENDDO
      ENDDO





      APHTIM=APHTIM+1.
      ARDSW =ARDSW +1.
      ARDLW =ARDLW +1.
      ASRFC =ASRFC +1.


      END SUBROUTINE TURBL



      SUBROUTINE UV_H_TO_V(NTSD,DT,NPHS,UZ0H,VZ0H,UZ0,VZ0               &
     &                    ,DUDT,DVDT,U,V,HBM2,IVE,IVW                   &
     &                    ,IDS,IDE,JDS,JDE,KDS,KDE                      &
     &                    ,IMS,IME,JMS,JME,KMS,KME                      &
     &                    ,ITS,ITE,JTS,JTE,KTS,KTE)





















      IMPLICIT NONE



      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,NPHS,NTSD

      INTEGER, DIMENSION(JMS:JME),INTENT(IN) :: IVE,IVW

      REAL,INTENT(IN) :: DT

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: HBM2,UZ0H,VZ0H

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: DUDT,DVDT

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: UZ0,VZ0

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: U,V







      INTEGER :: I,IEND,J,K

      REAL :: DTPHS

      LOGICAL :: E_BDY




      DTPHS=NPHS*DT
      E_BDY=(ITE>=IDE)





!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
      DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
        UZ0(I,J)=(UZ0H(I+IVE(J),J)*HBM2(I+IVE(J),J)                     &
     &           +UZ0H(I+IVW(J),J)*HBM2(I+IVW(J),J)                     &
     &           +UZ0H(I,J+1)*HBM2(I,J+1)+UZ0H(I,J-1)*HBM2(I,J-1))*0.25
        VZ0(I,J)=(VZ0H(I+IVE(J),J)*HBM2(I+IVE(J),J)                     &
     &           +VZ0H(I+IVW(J),J)*HBM2(I+IVW(J),J)                     &
     &           +VZ0H(I,J+1)*HBM2(I,J+1)+VZ0H(I,J-1)*HBM2(I,J-1))*0.25
      ENDDO
      ENDDO





!$omp parallel do                                                       &
!$omp& private(i,iend,j,k)
      DO K=KTS,KTE
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IEND=min(ide-(1),ite+(0))
          IF(E_BDY.AND.MOD(J,2)==1)IEND=IEND-1

          DO I=max(ids+(1),its-(0)),IEND
            U(I,J,K)=(DUDT(I+IVE(J),J,K)+DUDT(I+IVW(J),J,K)             &
     &               +DUDT(I,J+1,K)+DUDT(I,J-1,K))*0.25*DTPHS           &
     &               +U(I,J,K)
            V(I,J,K)=(DVDT(I+IVE(J),J,K)+DVDT(I+IVW(J),J,K)             &
     &               +DVDT(I,J+1,K)+DVDT(I,J-1,K))*0.25*DTPHS           &
     &               +V(I,J,K)
          ENDDO
        ENDDO
      ENDDO


      END SUBROUTINE UV_H_TO_V




      SUBROUTINE CUCNVC(NTSD,DT,NCNVC,NRADS,NRADL                       &
     &                 ,GPS,RESTRT,HYDRO                                &
     &                 ,CLDEFI,N_MOIST,ENSDIM                           &
     &                 ,MOIST                                           &
     &                 ,DETA1,DETA2,AETA1,AETA2,ETA1,ETA2               &
     &                 ,F_ICE,F_RAIN                                    &

     &                 ,APR_GR,APR_W,APR_MC,TTEN,QTEN                   &
     &                 ,APR_ST,APR_AS,APR_CAPMA                         &
     &                 ,APR_CAPME          ,APR_CAPMI                   &
     &                 ,MASS_FLUX         ,XF_ENS                       &
     &                 ,PR_ENS,GSW                                      &

     &                 ,PDTOP,PT,PD,RES,PINT,T,Q,CWM,TCUCN              &
     &                 ,OMGALF,U,V,W,Z,FIS,W0AVG                        &
     &                 ,PREC,ACPREC,CUPREC,CUPPT,CPRATE                 &
     &                 ,SM,HBM2,LPBL,CNVBOT,CNVTOP                      &
     &                 ,HTOP,HBOT,HTOPD,HBOTD,HTOPS,HBOTS               &
     &                 ,RTHBLTEN,RQVBLTEN,RTHRATEN                      &
     &                 ,AVCNVC,ACUTIM,IHE,IHW                           &
     &                 ,GRID,CONFIG_FLAGS                               &
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                         &
     &                 ,IMS,IME,JMS,JME,KMS,KME                         &
     &                 ,IPS,IPE,JPS,JPE,KPS,KPE                         &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE)






















      IMPLICIT NONE



      INTEGER,INTENT(IN) :: ENSDIM                                      &
     &                     ,IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,IPS,IPE,JPS,JPE,KPS,KPE                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,N_MOIST,NCNVC,NTSD,NRADS,NRADL

      INTEGER, DIMENSION(JMS:JME),INTENT(IN) :: IHE,IHW

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: LPBL

      REAL,INTENT(IN) :: DT,GPS,PDTOP,PT

      REAL,INTENT(INOUT) :: ACUTIM,AVCNVC

      REAL,DIMENSION(KMS:KME-1),INTENT(IN) :: AETA1,AETA2,DETA1,DETA2
      REAL,DIMENSION(KMS:KME  ),INTENT(IN) :: ETA1,ETA2

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FIS,HBM2,PD,RES,SM

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACPREC,CLDEFI    &
     &                                                ,CNVBOT,CNVTOP    &
     &                                                ,CUPPT,CUPREC     &
     &                                                ,HBOT,HTOP        &
     &                                                ,HBOTD,HTOPD      &
     &                                                ,HBOTS,HTOPS      &
     &                                                ,PREC,CPRATE      &
     &                 ,APR_GR,APR_W,APR_MC                             &
     &                 ,APR_ST,APR_AS,APR_CAPMA                         &
     &                 ,APR_CAPME,APR_CAPMI                             &
     &                 ,GSW,MASS_FLUX

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(IN) :: F_ICE       &
     &                                                     ,F_RAIN

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) :: QTEN     &
     &                                                        ,RQVBLTEN &
     &                                                        ,RTHBLTEN &
     &                                                        ,RTHRATEN &
     &                                                        ,TTEN

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: CWM      &
     &                                                        ,OMGALF   &
     &                                                        ,Q,T      &
     &                                                        ,TCUCN    &
     &                                                        ,U,V      &
     &                                                        ,W,Z

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: PINT

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) :: W0AVG

      REAL,DIMENSION(IMS:IME,JMS:JME,1:ENSDIM),INTENT(INOUT) :: PR_ENS  &
     &                                                         ,XF_ENS

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,N_MOIST)                   &
     &                                           ,INTENT(INOUT) :: MOIST


      LOGICAL,INTENT(IN) :: HYDRO,RESTRT

      TYPE(DOMAIN),TARGET :: GRID

      TYPE(GRID_CONFIG_REC_TYPE),INTENT(IN) :: CONFIG_FLAGS





      INTEGER :: I,ICLDCK,IENDX,ISTAT,J,K,MNTO,N,N_TIMSTPS_OUTPUT       &
     &          ,NCUBOT,NCUTOP,NSTEP_CNV

      INTEGER,DIMENSION(IMS:IME,JMS:JME) :: KPBL,LBOT,LOWLYR,LTOP

      REAL :: CAPA,CF_HI,DPL,DQDT,DTCNVC,DTDT,FICE,FRAIN,G_INV          &
     &       ,PCPCOL,PLYR,QI,QL_K,QR,QW,RDTCNVC,TL_K,WC,WMID

      REAL,DIMENSION(KMS:KME-1) :: QL,TL

      REAL,DIMENSION(IMS:IME,JMS:JME) :: CUBOT,CUTOP,NCA,PDSL           &
     &                                  ,RAINC,SFCZ,XLAND
      REAL,DIMENSION(IMS:IME,JMS:JME) :: RAINCV

      REAL,DIMENSION(ITS:ITE,JTS:JTE) :: WMID_L

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: DZ,P8W,P_PHY,PI_PHY    &
     &                                          ,RQCCUTEN,RQRCUTEN      &
     &                                          ,RQICUTEN,RQSCUTEN      &
     &                                          ,RQVCUTEN,RR,RTHCUTEN   &
     &                                          ,T_PHY,TH_PHY           &
     &                                          ,U_PHY,V_PHY,WINT

      REAL,DIMENSION(IMS:IME,JMS:JME,ENSDIM) :: ZERO_GD

      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: MOIST_TRANS

      LOGICAL :: RESTART,WARM_RAIN
      LOGICAL,DIMENSION(IMS:IME,JMS:JME) :: CU_ACT_FLAG

      CHARACTER(LEN=255) :: message




      INTEGER :: DTEMP_CHECK=2.0
      REAL :: TCHANGE











      IF(MOD(NTSD,NRADS)==0.OR.MOD(NTSD,NRADL)==0)THEN
         DO J=JMS,JME
         DO I=IMS,IME
           HTOP(I,J)=0.
           HBOT(I,J)=REAL(KTE+1)
           CUPPT(I,J)=0.
         ENDDO
         ENDDO
      ENDIF

      IF(MOD(NTSD,NCNVC)/=0.AND.                                      &
     &   CONFIG_FLAGS%CU_PHYSICS==BMJSCHEME)RETURN
      IF(MOD(NTSD,NCNVC)/=0.AND.                                      &
     &   CONFIG_FLAGS%CU_PHYSICS==SASSCHEME)RETURN

      NSTEP_CNV=NCNVC

      RESTART=RESTRT

      IF(CONFIG_FLAGS%CU_PHYSICS==KFETASCHEME)THEN

        IF(.NOT.RESTART.AND.NTSD==0)THEN
!$omp parallel do                                                       &
!$omp& private(i,j,k)
          DO J=JTS,JTE
          DO K=KTS,KTE
          DO I=ITS,ITE
            W0AVG(I,K,J)=0.
          ENDDO
          ENDDO
          ENDDO
        ENDIF

      ENDIF





      AVCNVC=AVCNVC+1.
      ACUTIM=ACUTIM+1.

      DTCNVC=NCNVC*DT
      RDTCNVC=1./DTCNVC
      CAPA=R_D/CP
      G_INV=1./G

!$omp parallel do                                                       &
!$omp& private(I,J)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
      DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))

        PDSL(I,J)=PD(I,J)*RES(I,J)
        RAINCV(I,J)=0.
        RAINC(I,J)=0.
        P8W(I,KTS,J)=PD(I,J)+PDTOP+PT
        LOWLYR(I,J)=KTS        
        XLAND(I,J)=SM(I,J)+1.
        NCA(I,J)=0.
        SFCZ(I,J)=FIS(I,J)*G_INV

        CUTOP(I,J)=HTOP(I,J)
        CUBOT(I,J)=HBOT(I,J)






        KPBL(I,J)=KTE-LPBL(I,J)+1
      ENDDO
      ENDDO

!$omp parallel do                                                       &
!$omp& private(dpl,fice,frain,i,j,k,plyr,qi,ql,qr,qw,wc)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        DO K=KTS,KTE
        DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
          DPL=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
          QL(K)=MAX(Q(I,J,K),EPSQ)
          PLYR=AETA1(K)*PDTOP+AETA2(K)*PDSL(I,J)+PT
          TL(K)=T(I,J,K)

          RR(I,K,J)=PLYR/(R_D*TL(K)*(P608*QL(K)+1.))
          T_PHY(I,K,J)=TL(K)

          TH_PHY(I,K,J)=TL(K)*(1.E5/PLYR)**CAPA

          P8W(I,K+1,J)=ETA1(K+1)*PDTOP+ETA2(K+1)*PDSL(I,J)+PT
          P_PHY(I,K,J)=PLYR
          PI_PHY(I,K,J)=(PLYR*1.E-5)**CAPA

          RTHCUTEN(I,K,J)=0.
          RQVCUTEN(I,K,J)=0.
          RQCCUTEN(I,K,J)=0.
          RQRCUTEN(I,K,J)=0.
          RQICUTEN(I,K,J)=0.
          RQSCUTEN(I,K,J)=0.
        ENDDO

      ENDDO
      ENDDO




      IF(.NOT.HYDRO)THEN
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO K=KTS,KTE
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
          DZ(I,K,J)=Z(I,J,K+1)-Z(I,J,K)
        ENDDO
        ENDDO
        ENDDO

        IF(NTSD==0)THEN
!$omp parallel do                                                       &
!$omp& private(i,j,k)
          DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          DO K=KTS,KTE+1   
          DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
            WINT(I,K,J)=0.
          ENDDO
          ENDDO
          ENDDO

	ELSE  

         DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
           DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
             WINT(I,KTS,J)=0.
             WINT(I,KTE+1,J)=0.
           ENDDO
         ENDDO

         DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          DO K=KTS+1,KTE
           DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
             WINT(I,K,J)=0.5*(W(I,J,K)+W(I,J,K-1))
           ENDDO
          ENDDO
         ENDDO

        ENDIF
	
      ELSE   

        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
          WINT(I,KTS,J)=0.
          WINT(I,KTE+1,J)=0.
        ENDDO
        ENDDO

!$omp parallel do                                                       &
!$omp& private(i,j,k,plyr)
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
            WMID_L(I,J)=-OMGALF(I,J,KTS)*CP/(G*DT)
            PDSL(I,J)=PD(I,J)*RES(I,J)
            PLYR=AETA1(KTS)*PDTOP+AETA2(KTS)*PDSL(I,J)+PT
            DZ(I,KTS,J)=T(I,J,KTS)*(P608*Q(I,J,KTS)+1.)*R_D             &
     &                 *(P8W(I,KTS,J)-P8W(I,KTS+1,J))                   &
     &                 /(PLYR*G)
          ENDDO
        ENDDO

!$omp parallel do                                                       &
!$omp& private(i,j,k,ql_k,tl_k,wmid)
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          DO K=KTS+1,KTE
          DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
            TL_K=T_PHY(I,K,J)
            QL_K=MAX(Q(I,J,K),EPSQ)
            WMID=-OMGALF(I,J,K)*CP/(G*DT)
            WINT(I,K,J)=0.5*(WMID_L(I,J)+WMID)
            WMID_L(I,J)=WMID
            DZ(I,K,J)=TL_K*(P608*QL_K+1.)*R_D                           &
     &               *(P8W(I,K,J)-P8W(I,K+1,J))                         &
     &               /(P_PHY(I,K,J)*G)
          ENDDO
          ENDDO
        ENDDO

      ENDIF





      IF(CONFIG_FLAGS%CU_PHYSICS/=BMJSCHEME)THEN

!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO K=KTS,KTE

          DO J=max(jds+(1),jts-(1)),min(jde-(1),jte+(1))
          DO I=max(ids+(0),its-(1)),min(ide-(0),ite+(1))
            U_PHY(I,K,J)=(U(I+IHE(J),J,K)+U(I+IHW(J),J,K)               &
     &                   +U(I,J+1,K)+U(I,J-1,K))                        &
     &                   *0.25
            V_PHY(I,K,J)=(V(I+IHE(J),J,K)+V(I+IHW(J),J,K)               &
     &                   +V(I,J+1,K)+V(I,J-1,K))                        &
     &                   *0.25
          ENDDO
          ENDDO

        ENDDO

      ENDIF





      IF(.NOT.ALLOCATED(MOIST_TRANS))THEN
        ALLOCATE(MOIST_TRANS(IMS:IME,KMS:KME,JMS:JME,N_MOIST),STAT=ISTAT)
      ENDIF

      DO N=1,N_MOIST
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO K=KMS,KME
        DO J=JMS,JME
        DO I=IMS,IME
          MOIST_TRANS(I,K,J,N)=MOIST(I,J,K,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO







      CALL SET_TILES(GRID,IDS+1,IDE-1,JDS+2,JDE-2,ITS,ITE,JTS,JTE)

CALL CUMULUS_DRIVER(GRID,IDS=IDS,IDE=IDE,JDS=JDS,JDE=JDE,KDS=KDS,KDE=KDE,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=KMS,KME=KME,IPS=ips,I&
&PE=ipe,JPS=jps,JPE=jpe,KPS=kps,KPE=kpe,I_START=GRID%I_START,I_END=GRID%I_END,J_START=GRID%J_START,J_END=GRID%J_END,KTS=KTS,KTE=KTE&
&,NUM_TILES=GRID%NUM_TILES,U=U_PHY,V=V_PHY,TH=TH_PHY,T=T_PHY,W=WINT,P=P_PHY,PI=PI_PHY,RHO=RR,W0AVG=W0AVG,ITIMESTEP=NTSD,DT=DT,DX=GP&
&S,RAINC=RAINC,RAINCV=RAINCV,NCA=NCA,DZ8W=DZ,P8W=P8W,FORCET=TTEN,FORCEQ=QTEN,CLDEFI=CLDEFI,LOWLYR=LOWLYR,XLAND=XLAND,CU_ACT_FLAG=CU&
&_ACT_FLAG,WARM_RAIN=WARM_RAIN,STEPCU=NSTEP_CNV,GSW=GSW,PERIODIC_X=.FALSE.,PERIODIC_Y=.FALSE.,HTOP=CUTOP,HBOT=CUBOT,KPBL=KPBL,HT=SF&
&CZ,APR_GR=APR_GR,APR_W=APR_W,APR_MC=APR_MC,APR_ST=APR_ST,APR_AS=APR_AS,APR_CAPMA=APR_CAPMA,APR_CAPME=APR_CAPME,APR_CAPMI=APR_CAPMI&
&,MASS_FLUX=MASS_FLUX,XF_ENS=XF_ENS,PR_ENS=PR_ENS,ENSDIM=ENSDIM,MAXIENS=1,MAXENS=3,MAXENS2=3,MAXENS3=16,RTHCUTEN=RTHCUTEN,RQVCUTEN=&
&RQVCUTEN,RQCCUTEN=RQCCUTEN,RQRCUTEN=RQRCUTEN,RQICUTEN=RQICUTEN,RQSCUTEN=RQSCUTEN,RTHBLTEN=RTHBLTEN,RQVBLTEN=RQVBLTEN,RTHRATEN=RTHR&
&ATEN,CU_PHYSICS=CONFIG_FLAGS%CU_PHYSICS,QV_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QV),F_QV=F_QV,QC_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QC),F_QC&
&=F_QC,QR_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QR),F_QR=F_QR,QI_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QI),F_QI=F_QI,QS_CURR=MOIST_TRANS(IMS,KMS,&
&JMS,P_QS),F_QS=F_QS,QG_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QG),F_QG=F_QG)








      CF_HI=CONFIG_FLAGS%HISTORY_INTERVAL
      N_TIMSTPS_OUTPUT=NINT(60.*CF_HI/DT)
      MNTO=MOD(NTSD,N_TIMSTPS_OUTPUT)

      IF(MNTO>0.AND.MNTO<=NCNVC)THEN
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        IENDX=min(ide-(1),ite+(0))
        IF(MOD(J,2)==0.AND.ITE==IDE-1)IENDX=IENDX-1
        DO I=max(ids+(1),its-(0)),IENDX
          CNVBOT(I,J)=REAL(KTE+1.)
          CNVTOP(I,J)=0.
          HBOTD(I,J)=REAL(KTE+1.)
          HTOPD(I,J)=0.
          HBOTS(I,J)=REAL(KTE+1.)
          HTOPS(I,J)=0.
        ENDDO
        ENDDO
      ENDIF



!$omp parallel do                                                       &
!$omp& private(i,iendx,j,ncubot,ncutop,pcpcol)
      pcp_cloud: DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        IENDX=min(ide-(1),ite+(0))
        IF(MOD(J,2)==0.AND.ITE==IDE-1)IENDX=IENDX-1
        DO I=max(ids+(1),its-(0)),IENDX



          PCPCOL=RAINCV(I,J)*1.E-3*NSTEP_CNV
          PREC(I,J)=PREC(I,J)+PCPCOL
          ACPREC(I,J)=ACPREC(I,J)+PCPCOL
          CUPREC(I,J)=CUPREC(I,J)+PCPCOL
          CUPPT(I,J)=CUPPT(I,J)+PCPCOL
          CPRATE(I,J)=PCPCOL





          CUTOP(I,J)=MIN(CUTOP(I,J),REAL(KDE))
          CUTOP(I,J)=MAX(CUTOP(I,J),0.0)
          CUBOT(I,J)=MIN(CUBOT(I,J),REAL(KDE))
          CUBOT(I,J)=MAX(CUBOT(I,J),0.0)

          NCUTOP=NINT(CUTOP(I,J))
          NCUBOT=NINT(CUBOT(I,J))

          IF(NCUTOP>1.AND.NCUTOP<KDE)THEN
            HTOP(I,J)=MAX(CUTOP(I,J),HTOP(I,J))
            CNVTOP(I,J)=MAX(CUTOP(I,J),CNVTOP(I,J))
            IF(PCPCOL>0.)THEN
              HTOPD(I,J)=MAX(CUTOP(I,J),HTOPD(I,J))
            ELSE
              HTOPS(I,J)=MAX(CUTOP(I,J),HTOPS(I,J))
            ENDIF
          ENDIF

          IF(NCUBOT>0.AND.NCUBOT<KDE)THEN
            HBOT(I,J)=MIN(CUBOT(I,J),HBOT(I,J))
            CNVBOT(I,J)=MIN(CUBOT(I,J),CNVBOT(I,J))
            IF(PCPCOL>0.)THEN
              HBOTD(I,J)=MIN(CUBOT(I,J),HBOTD(I,J))
            ELSE
              HBOTS(I,J)=MIN(CUBOT(I,J),HBOTS(I,J))
            ENDIF
          ENDIF

        ENDDO
      ENDDO pcp_cloud





!$omp parallel do                                                       &
!$omp& private(dqdt,dtdt,i,iendx,j,k,tchange)
      DO K=KTS,KTE
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        IENDX=min(ide-(1),ite+(0))
        IF(MOD(J,2)==0.AND.ITE==IDE-1)IENDX=IENDX-1
        DO I=max(ids+(1),its-(0)),IENDX




          DQDT=RQVCUTEN(I,K,J)/(1.+MOIST_TRANS(I,K,J,P_QV))**2



          DTDT=RTHCUTEN(I,K,J)*PI_PHY(I,K,J)
          T(I,J,K)=T(I,J,K)+DTDT*DTCNVC
          Q(I,J,K)=Q(I,J,K)+DQDT*DTCNVC
          TCUCN(I,J,K)=TCUCN(I,J,K)+DTDT
          MOIST_TRANS(I,K,J,P_QV)=Q(I,J,K)/(1.-Q(I,J,K))       

          cps_select: SELECT CASE(config_flags%cu_physics)

          CASE (KFSCHEME,KFETASCHEME,GDSCHEME,SASSCHEME)
            IF(CONFIG_FLAGS%MP_PHYSICS==ETAMPNEW)THEN
              MOIST_TRANS(I,K,J,P_QS)=MAX(0.,MOIST_TRANS(I,K,J,P_QS)+RQICUTEN(I,K,J)*DTCNVC+RQSCUTEN(I,K,J)*DTCNVC)
            ELSE
              MOIST_TRANS(I,K,J,P_QI)=MAX(0.,MOIST_TRANS(I,K,J,P_QI)+RQICUTEN(I,K,J)*DTCNVC)
              MOIST_TRANS(I,K,J,P_QS)=MAX(0.,MOIST_TRANS(I,K,J,P_QS)+RQSCUTEN(I,K,J)*DTCNVC)
            ENDIF
            MOIST_TRANS(I,K,J,P_QR)=MAX(0.,MOIST_TRANS(I,K,J,P_QR)+RQRCUTEN(I,K,J)*DTCNVC)
            MOIST_TRANS(I,K,J,P_QC)=MAX(0.,MOIST_TRANS(I,K,J,P_QC)+RQCCUTEN(I,K,J)*DTCNVC)
          END SELECT cps_select

          TCHANGE=DTDT*DTCNVC
	  IF(ABS(TCHANGE)>DTEMP_CHECK)THEN
            WRITE(message,*)'BIG T CHANGE BY CONVECTION=',TCHANGE             &
                     ,' AT (',I,',',J,',',K,') FOR NTSD=',NTSD
            CALL wrf_message(trim(message))
	  ENDIF

        ENDDO
      ENDDO
      ENDDO




      DO N=1,N_MOIST
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO J=JMS,JME
        DO K=KMS,KME
        DO I=IMS,IME
          MOIST(I,J,K,N)=MOIST_TRANS(I,K,J,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO



      DEALLOCATE(MOIST_TRANS,STAT=ISTAT)



      END SUBROUTINE CUCNVC



      SUBROUTINE GSMDRIVE(NTSD,DT,NPHS,N_MOIST                          &
     &                   ,DX,DY,SM,HBM2,FIS                             &
     &                   ,DETA1,DETA2,AETA1,AETA2,ETA1,ETA2             &
     &                   ,PDTOP,PT,PD,RES,PINT,T,Q,CWM,TRAIN            &
     &                   ,MOIST,SCALAR,N_SCALAR                         &
     &                   ,F_ICE,F_RAIN,F_RIMEF,SR                       &
     &                   ,PREC,ACPREC,AVRAIN                            &
     &                   ,MP_RESTART_STATE                              &
     &                   ,TBPVS_STATE                                   &
     &                   ,TBPVS0_STATE                                  &
     &                   ,GRID,CONFIG_FLAGS                             &
     &                   ,IDS,IDE,JDS,JDE,KDS,KDE                       &
     &                   ,IMS,IME,JMS,JME,KMS,KME                       &
     &                   ,ITS,ITE,JTS,JTE,KTS,KTE)






















      IMPLICIT NONE



      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,N_MOIST,N_SCALAR,NPHS,NTSD

      REAL,INTENT(IN) :: DT,DX,DY,PDTOP,PT

      REAL,INTENT(INOUT) :: AVRAIN

      REAL,DIMENSION(KMS:KME-1),INTENT(IN) :: AETA1,AETA2,DETA1,DETA2
      REAL,DIMENSION(KMS:KME),INTENT(IN) :: ETA1,ETA2

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FIS,HBM2,PD,RES,SM

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: PINT

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: ACPREC,PREC

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(INOUT) :: CWM,Q    &
     &                                                        ,T,TRAIN

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(INOUT) :: F_ICE    &   
     &                                                        ,F_RAIN   &
     &                                                        ,F_RIMEF

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,N_MOIST)                   &
     &                                           ,INTENT(INOUT) :: MOIST
      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,N_SCALAR)                  &
     &                                          ,INTENT(INOUT) :: SCALAR



      REAL,DIMENSION(:),INTENT(INOUT) :: MP_RESTART_STATE               &
     &                                  ,TBPVS_STATE,TBPVS0_STATE

      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: SR

      TYPE(DOMAIN),TARGET :: GRID

      TYPE(GRID_CONFIG_REC_TYPE),INTENT(IN) :: CONFIG_FLAGS





      INTEGER :: I,IENDX,IJ,ISTAT,J,K,N

      INTEGER,DIMENSION(IMS:IME,JMS:JME) :: LOWLYR

      REAL :: CAPA,DPL,DTPHS,PCPCOL,PLYR,RDTPHS,RG,TNEW

      REAL,DIMENSION(KMS:KME-1) :: QL,TL

      REAL,DIMENSION(IMS:IME,JMS:JME) :: CUBOT,CUTOP,PDSL               &
     &                                  ,RAINNC,RAINNCV,XLAND

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME) :: CWM_PHY,DZ             &
     &                                          ,P8W,P_PHY,PI_PHY       &
     &                                          ,RR,T_PHY,TH_PHY

      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: MOIST_TRANS
      REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: SCALAR_TRANS

      LOGICAL :: E_BDY,F_QT,QT_PRESENT,WARM_RAIN





      ALLOCATE(MOIST_TRANS(IMS:IME,KMS:KME,JMS:JME,N_MOIST),STAT=ISTAT)
      ALLOCATE(SCALAR_TRANS(IMS:IME,KMS:KME,JMS:JME,N_SCALAR),STAT=ISTAT)





      DO N=1,N_MOIST
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO K=KMS,KME
        DO J=JMS,JME
        DO I=IMS,IME
          MOIST_TRANS(I,K,J,N)=MOIST(I,J,K,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO



      IF(CONFIG_FLAGS%MP_PHYSICS/=ETAMPNEW)THEN
        DO N=1,N_SCALAR
!$omp parallel do                                                       &
!$omp& private(i,j,k)
          DO K=KMS,KME
          DO J=JMS,JME
          DO I=IMS,IME
            SCALAR_TRANS(I,K,J,N)=SCALAR(I,J,K,N)
          ENDDO
          ENDDO
          ENDDO
        ENDDO
      ENDIF



      IF(CONFIG_FLAGS%MP_PHYSICS==ETAMPNEW)THEN
        QT_PRESENT=.TRUE.
      ELSE
        QT_PRESENT=.FALSE.
      ENDIF

      DTPHS=NPHS*DT
      RDTPHS=1./DTPHS
      CAPA=R_D/CP
      RG=1./G
      AVRAIN=AVRAIN+1.






!$omp parallel do                                                       &
!$omp& private(i,j)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
      DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))

        PDSL(I,J)=PD(I,J)*RES(I,J)
        P8W(I,KTE+1,J)=PT
        LOWLYR(I,J)=KTS        
        XLAND(I,J)=SM(I,J)+1.





        RAINNC(I,J)=0.

      ENDDO
      ENDDO





!$omp parallel do                                                       &
!$omp& private(dpl,i,j,k,plyr,ql,tl)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        DO K=KTS,KTE
        DO I=max(ids+(1),its-(0)),min(ide-(1),ite+(0))
          DPL=DETA1(K)*PDTOP+DETA2(K)*PDSL(I,J)
          QL(K)=MAX(Q(I,J,K),EPSQ)

          PLYR=(PINT(I,J,K)+PINT(I,J,K+1))*0.5
          TL(K)=T(I,J,K)

          RR(I,K,J)=PLYR/(R_D*TL(K)*(P608*QL(K)+1.))
          T_PHY(I,K,J)=TL(K)
          PI_PHY(I,K,J)=(PLYR*1.E-5)**CAPA
          TH_PHY(I,K,J)=TL(K)/PI_PHY(I,K,J)

          P8W(I,K,J)=ETA1(K)*PDTOP+ETA2(K)*PDSL(I,J)+PT
          P_PHY(I,K,J)=PLYR
          DZ(I,K,J)=DPL*RG/RR(I,K,J)
          CWM_PHY(I,K,J)=CWM(I,J,K)
        ENDDO

      ENDDO
      ENDDO






      CALL SET_TILES(GRID,IDS+1,IDE-1,JDS+2,JDE-2,ITS,ITE,JTS,JTE)

      CALL MICROPHYSICS_DRIVER(                                         &
     &                  TH=TH_PHY,RHO=RR,PI_PHY=PI_PHY,P=P_PHY          &
     &                 ,RAINNC=RAINNC,RAINNCV=RAINNCV                   &
     &                 ,DZ8W=DZ,P8W=P8W,DT=DTPHS,DX=DX,DY=DY            &
     &                 ,MP_PHYSICS=CONFIG_FLAGS%MP_PHYSICS              &
     &                 ,SPECIFIED=CONFIG_FLAGS%SPECIFIED                &
     &                        .OR.CONFIG_FLAGS%NESTED                   &
     &                 ,SPEC_ZONE=0,WARM_RAIN=WARM_RAIN                 &
     &                 ,XLAND=XLAND,ITIMESTEP=NTSD-1                    &
     &                 ,F_ICE_PHY=F_ICE,F_RAIN_PHY=F_RAIN               &
     &                 ,F_RIMEF_PHY=F_RIMEF                             &
     &                 ,LOWLYR=LOWLYR,SR=SR                             &
     &                 ,QV_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QV),F_QV=F_QV &
     &                 ,QC_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QC),F_QC=F_QC &
     &                 ,QR_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QR),F_QR=F_QR &
     &                 ,QI_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QI),F_QI=F_QI &
     &                 ,QS_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QS),F_QS=F_QS &
     &                 ,QG_CURR=MOIST_TRANS(IMS,KMS,JMS,P_QG),F_QG=F_QG &
     &                 ,QNI_CURR=SCALAR_TRANS(IMS,KMS,JMS,P_QNI),F_QNI=F_QNI  &
     &                 ,QNR_CURR=SCALAR_TRANS(IMS,KMS,JMS,P_QNR),F_QNR=F_QNR  &
     &                 ,QT_CURR=CWM_PHY,F_QT=QT_PRESENT                 &
     &                 ,MP_RESTART_STATE=MP_RESTART_STATE               &
     &                 ,TBPVS_STATE=TBPVS_STATE                         &
     &                 ,TBPVS0_STATE=TBPVS0_STATE                       &
     &                 ,IDS=IDS,IDE=IDE,JDS=JDS,JDE=JDE,KDS=KDS,KDE=KDE &
     &                 ,IMS=IMS,IME=IME,JMS=JMS,JME=JME,KMS=KMS,KME=KME &
     &                 ,I_START=GRID%I_START,I_END=GRID%I_END           &
     &                 ,J_START=GRID%J_START,J_END=GRID%J_END           &
     &                 ,KTS=KTS,KTE=KTE,NUM_TILES=GRID%NUM_TILES        &
                                                                        )

!$omp parallel do                                                       &
!$omp& private(ij)
      DO IJ=1,GRID%NUM_TILES
        CALL MICROPHYSICS_ZERO_OUT(                                     &
                     MOIST_TRANS,N_MOIST,CONFIG_FLAGS                   &
                    ,IDS,IDE,JDS,JDE,KDS,KDE                            &
                    ,IMS,IME,JMS,JME,KMS,KME                            &
                    ,GRID%I_START(IJ),GRID%I_END(IJ)                    &
                    ,GRID%J_START(IJ),GRID%J_END(IJ)                    &
                    ,KTS,KTE                                       )
      ENDDO






      E_BDY=(ITE>=IDE)




!$omp parallel do                                                       &
!$omp& private(i,iendx,j,k,tnew)
      DO K=KTS,KTE
        DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
          IENDX=min(ide-(1),ite+(0))
          IF(E_BDY.AND.MOD(J,2)==0)IENDX=IENDX-1
          DO I=max(ids+(1),its-(0)),IENDX
            TNEW=TH_PHY(I,K,J)*PI_PHY(I,K,J)
            TRAIN(I,J,K)=TRAIN(I,J,K)+(TNEW-T(I,J,K))*RDTPHS
            T(I,J,K)=TNEW
            Q(I,J,K)=MOIST_TRANS(I,K,J,P_QV)/(1.+MOIST_TRANS(I,K,J,P_QV))
            CWM(I,J,K)=CWM_PHY(I,K,J)
          ENDDO
        ENDDO
      ENDDO







!$omp parallel do                                                       &
!$omp& private(i,iendx,j,pcpcol)
      DO J=max(jds+(2),jts-(0)),min(jde-(2),jte+(0))
        IENDX=min(ide-(1),ite+(0))
        IF(E_BDY.AND.MOD(J,2)==0)IENDX=IENDX-1
        DO I=max(ids+(1),its-(0)),IENDX
          PCPCOL=RAINNCV(I,J)*1.E-3
          PREC(I,J)=PREC(I,J)+PCPCOL
          ACPREC(I,J)=ACPREC(I,J)+PCPCOL
        ENDDO
      ENDDO





      DO N=1,N_MOIST
!$omp parallel do                                                       &
!$omp& private(i,j,k)
        DO J=JMS,JME
        DO K=KMS,KME
        DO I=IMS,IME
          MOIST(I,J,K,N)=MOIST_TRANS(I,K,J,N)
        ENDDO
        ENDDO
        ENDDO
      ENDDO



      IF(CONFIG_FLAGS%MP_PHYSICS/=ETAMPNEW)THEN
        DO N=1,N_SCALAR
!$omp parallel do                                                       &
!$omp& private(i,j,k)
          DO J=JMS,JME
          DO K=KMS,KME
          DO I=IMS,IME
            SCALAR(I,J,K,N)=SCALAR_TRANS(I,K,J,N)
          ENDDO
          ENDDO
          ENDDO
        ENDDO
      ENDIF



      DEALLOCATE(MOIST_TRANS,STAT=ISTAT)
      DEALLOCATE(SCALAR_TRANS,STAT=ISTAT)



      END SUBROUTINE GSMDRIVE



      SUBROUTINE UPDATE_MOIST(MOIST,Q,CWM,F_ICE,F_RAIN,N_MOIST          &
     &                       ,IDS,IDE,JDS,JDE,KDS,KDE                   &
     &                       ,IMS,IME,JMS,JME,KMS,KME                   &
     &                       ,ITS,ITE,JTS,JTE,KTS,KTE)



      IMPLICIT NONE



      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                     &
     &                     ,IMS,IME,JMS,JME,KMS,KME                     &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE                     &
     &                     ,N_MOIST

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME),INTENT(IN) :: CWM,Q

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(IN) :: F_ICE       &   
     &                                                     ,F_RAIN

      REAL,DIMENSION(IMS:IME,JMS:JME,KMS:KME,N_MOIST),INTENT(OUT) :: MOIST





      INTEGER :: I,J,K

      REAL :: FICE,FRAIN,QI,QR,QW,WC





      DO K=KTS,KTE
      DO J=max(jds+(0),jts-(0)),min(jde-(0),jte+(0))
      DO I=max(ids+(0),its-(0)),min(ide-(0),ite+(0))
        MOIST(I,J,K,P_QV)=Q(I,J,K)/(1.-Q(I,J,K))
        WC=CWM(I,J,K)
        QI=0.
        QR=0.
        QW=0.
        FICE=F_ICE(I,K,J)
        FRAIN=F_RAIN(I,K,J)

        IF(FICE>=1.)THEN
          QI=WC
        ELSEIF(FICE<=0.)THEN
          QW=WC
        ELSE
          QI=FICE*WC
          QW=WC-QI
        ENDIF

        IF(QW>0..AND.FRAIN>0.)THEN
          IF(FRAIN>=1.)THEN
            QR=QW
            QW=0.
          ELSE
            QR=FRAIN*QW
            QW=QW-QR
          ENDIF
        ENDIF

        MOIST(I,J,K,P_QC)=QW
        MOIST(I,J,K,P_QR)=QR
        MOIST(I,J,K,P_QI)=0.
        MOIST(I,J,K,P_QS)=QI
        MOIST(I,J,K,P_QG)=0.
      ENDDO
      ENDDO
      ENDDO



      END SUBROUTINE UPDATE_MOIST





      END MODULE MODULE_PHYSICS_CALLS


