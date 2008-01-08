!-----------------------------------------------------------------------
!
      MODULE MODULE_BL_MYJPBL
!
!-----------------------------------------------------------------------
!
      USE MODULE_MODEL_CONSTANTS
!
!-----------------------------------------------------------------------
!
! REFERENCES:  Janjic (2002), NCEP Office Note 437
!              Mellor and Yamada (1982), Rev. Geophys. Space Phys.
!
! ABSTRACT:
!     MYJ UPDATES THE TURBULENT KINETIC ENERGY WITH THE PRODUCTION/
!     DISSIPATION TERM AND THE VERTICAL DIFFUSION TERM
!     (USING AN IMPLICIT FORMULATION) FROM MELLOR-YAMADA
!     LEVEL 2.5 AS EXTENDED BY JANJIC.  EXCHANGE COEFFICIENTS FOR
!     THE SURFACE AND FOR ALL LAYER INTERFACES ARE COMPUTED FROM
!     MONIN-OBUKHOV THEORY.
!     THE TURBULENT VERTICAL EXCHANGE IS THEN EXECUTED.
!
!-----------------------------------------------------------------------
!
      INTEGER :: ITRMX=5 ! Iteration count for mixing length computation
!
!     REAL,PARAMETER :: G=9.81,PI=3.1415926,R_D=287.04,R_V=461.6        &
!    &                 ,VKARMAN=0.4
      REAL,PARAMETER :: PI=3.1415926,VKARMAN=0.4
!     REAL,PARAMETER :: CP=7.*R_D/2.
      REAL,PARAMETER :: CAPA=R_D/CP
      REAL,PARAMETER :: RLIVWV=XLS/XLV,ELOCP=2.72E6/CP
      REAL,PARAMETER :: EPS1=1.E-12,EPS2=0.
      REAL,PARAMETER :: EPSL=0.32,EPSRU=1.E-7,EPSRS=1.E-7               &
     &                 ,EPSTRB=1.E-24
      REAL,PARAMETER :: EPSA=1.E-8,EPSIT=1.E-4,EPSU2=1.E-4,EPSUST=0.07  &
     &                 ,FH=1.01 
      REAL,PARAMETER :: ALPH=0.30,BETA=1./273.,EL0MAX=1000.,EL0MIN=1.   &
     &                 ,ELFC=0.23*0.5,GAM1=0.2222222222222222222        &
     &                 ,PRT=1.
      REAL,PARAMETER :: A1=0.659888514560862645                         &
     &                 ,A2x=0.6574209922667784586                       &
     &                 ,B1=11.87799326209552761                         &
     &                 ,B2=7.226971804046074028                         &
     &                 ,C1=0.000830955950095854396
      REAL,PARAMETER :: A2S=17.2693882,A3S=273.16,A4S=35.86
      REAL,PARAMETER :: ELZ0=0.,ESQ=5.0,EXCM=0.001                      &
     &                 ,FHNEU=0.8,GLKBR=10.,GLKBS=30.                   &
     &                 ,QVISC=2.1E-5,RFC=0.191,RIC=0.505,SMALL=0.35     &
     &                 ,SQPR=0.84,SQSC=0.84,SQVISC=258.2,TVISC=2.1E-5   &
     &                 ,USTC=0.7,USTR=0.225,VISC=1.5E-5                 &
     &                 ,WOLD=0.15,WWST=1.2,ZTMAX=1.,ZTFC=1.,ZTMIN=-5.
!
      REAL,PARAMETER :: SEAFC=0.98,PQ0SEA=PQ0*SEAFC
!
      REAL,PARAMETER :: BTG=BETA*G,CZIV=SMALL*GLKBS                     &
!    &                 ,EP_1=R_V/R_D-1.,ESQHF=0.5*5.0,GRRS=GLKBR/GLKBS  &
     &                 ,ESQHF=0.5*5.0,GRRS=GLKBR/GLKBS                  &
     &                 ,RB1=1./B1,RTVISC=1./TVISC,RVISC=1./VISC         &
     &                 ,ZQRZT=SQSC/SQPR
!
      REAL,PARAMETER :: ADNH= 9.*A1*A2x*A2x*(12.*A1+3.*B2)*BTG*BTG      &                  
     &                 ,ADNM=18.*A1*A1*A2x*(B2-3.*A2x)*BTG              & 
     &                 ,ANMH=-9.*A1*A2x*A2x*BTG*BTG                     &
     &                 ,ANMM=-3.*A1*A2x*(3.*A2x+3.*B2*C1+18.*A1*C1-B2)  &
     &                                *BTG                              &   
     &                 ,BDNH= 3.*A2x*(7.*A1+B2)*BTG                     &
     &                 ,BDNM= 6.*A1*A1                                  &
     &                 ,BEQH= A2x*B1*BTG+3.*A2x*(7.*A1+B2)*BTG          &
     &                 ,BEQM=-A1*B1*(1.-3.*C1)+6.*A1*A1                 &
     &                 ,BNMH=-A2x*BTG                                   &     
     &                 ,BNMM=A1*(1.-3.*C1)                              &
     &                 ,BSHH=9.*A1*A2x*A2x*BTG                          &
     &                 ,BSHM=18.*A1*A1*A2x*C1                           &
     &                 ,BSMH=-3.*A1*A2x*(3.*A2x+3.*B2*C1+12.*A1*C1-B2)  &
     &                                *BTG                              &
     &                 ,CESH=A2x                                        &
     &                 ,CESM=A1*(1.-3.*C1)                              &
     &                 ,CNV=EP_1*G/BTG                                  &
     &                 ,ELFCS=VKARMAN*BTG                               &
     &                 ,FZQ1=RTVISC*QVISC*ZQRZT                         &
     &                 ,FZQ2=RTVISC*QVISC*ZQRZT                         &
     &                 ,FZT1=RVISC *TVISC*SQPR                          &
     &                 ,FZT2=CZIV*GRRS*TVISC*SQPR                       &
     &                 ,FZU1=CZIV*VISC                                  &
     &                 ,PIHF=0.5*PI                                     &
     &                 ,RFAC=RIC/(FHNEU*RFC*RFC)                        &
     &                 ,RQVISC=1./QVISC                                 &
     &                 ,RRIC=1./RIC                                     &
     &                 ,USTFC=0.018/G                                   &
     &                 ,WNEW=1.-WOLD                                    &
     &                 ,WWST2=WWST*WWST
!
!-----------------------------------------------------------------------
!***  FREE TERM IN THE EQUILIBRIUM EQUATION FOR (L/Q)**2
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: AEQH=9.*A1*A2x*A2x*B1*BTG*BTG                   &
     &                      +9.*A1*A2x*A2x*(12.*A1+3.*B2)*BTG*BTG       &
     &                 ,AEQM=3.*A1*A2x*B1*(3.*A2x+3.*B2*C1+18.*A1*C1-B2)&
     &                      *BTG+18.*A1*A1*A2x*(B2-3.*A2x)*BTG
!
!-----------------------------------------------------------------------
!***  FORBIDDEN TURBULENCE AREA
!-----------------------------------------------------------------------
!
      REAL,PARAMETER :: REQU=-AEQH/AEQM                                 &
     &                 ,EPSGH=1.E-9,EPSGM=REQU*EPSGH
!
!-----------------------------------------------------------------------
!***  NEAR ISOTROPY FOR SHEAR TURBULENCE, WW/Q2 LOWER LIMIT
!-----------------------------------------------------------------------
! 
      REAL,PARAMETER :: UBRYL=(18.*REQU*A1*A1*A2x*B2*C1*BTG             &
     &                         +9.*A1*A2x*A2x*B2*BTG*BTG)               &
     &                        /(REQU*ADNM+ADNH)                         &
     &                 ,UBRY=(1.+EPSRS)*UBRYL,UBRY3=3.*UBRY
!
      REAL,PARAMETER :: AUBH=27.*A1*A2x*A2x*B2*BTG*BTG-ADNH*UBRY3       &
     &                 ,AUBM=54.*A1*A1*A2x*B2*C1*BTG -ADNM*UBRY3        &
     &                 ,BUBH=(9.*A1*A2x+3.*A2x*B2)*BTG-BDNH*UBRY3       &
     &                 ,BUBM=18.*A1*A1*C1           -BDNM*UBRY3         &
     &                 ,CUBR=1.                     -     UBRY3         &
     &                 ,RCUBR=1./CUBR
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!----------------------------------------------------------------------
      SUBROUTINE MYJPBL(DT,STEPBL,HT,DZ                                &
     &                 ,PMID,PINT,TH,T,EXNER,QV,CWM,U,V,RHO            &
     &                 ,TSK,QSFC,CHKLOWQ,THZ0,QZ0,UZ0,VZ0              &
     &                 ,LOWLYR,XLAND,SICE,SNOW                         &
     &                 ,TKE_MYJ,EXCH_H,USTAR,ZNT,EL_MYJ,PBLH,KPBL,CT   &
     &                 ,AKHS,AKMS,ELFLX                                &
     &                 ,RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,RQCBLTEN     &
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                 ,IMS,IME,JMS,JME,KMS,KME                        &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: STEPBL

      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: LOWLYR
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: KPBL
!
      REAL,INTENT(IN) :: DT
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: HT,SICE,SNOW       &
     &                                             ,TSK,XLAND
!
      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(IN) :: CWM,DZ     &
     &                                                     ,EXNER      &
     &                                                     ,PMID,PINT  &
     &                                                     ,QV,RHO     &
     &                                                     ,T,TH,U,V   
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: PBLH
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: AKHS,AKMS
!
      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME)                          &
     &    ,INTENT(OUT) ::                      EL_MYJ                  &
     &                                        ,RQCBLTEN,RQVBLTEN       &
     &                                        ,RTHBLTEN                &
     &                                        ,RUBLTEN,RVBLTEN        
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: CT,QSFC,QZ0     &
     &                                                ,THZ0,USTAR      &
     &                                                ,UZ0,VZ0,ZNT
!
      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME)                          &
     &    ,INTENT(INOUT) ::                    EXCH_H,TKE_MYJ
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: CHKLOWQ,ELFLX
!
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: I,J,K,KFLIP,LLOW,LMH,LMXL
!
      INTEGER,DIMENSION(ITS:ITE,JTS:JTE) :: LPBL
!
      REAL :: AKHS_DENS,AKMS_DENS,APEX,DCDT,DELTAZ,DQDT,DTDIF,DTDT     &
     &       ,DTTURBL,DUDT,DVDT,EXNSFC,PSFC,PTOP,QFC1,QLOW,QOLD        &
     &       ,RATIOMX,RDTTURBL,RG,RWMSK,SEAMASK,THNEW,THOLD,TX         &
     &       ,ULOW,VLOW,WMSK
!
      REAL,DIMENSION(KTS:KTE) :: CWMK,PK,Q2K,QK,THEK,TK,UK,VK
!
      REAL,DIMENSION(KTS:KTE-1) :: AKHK,AKMK,EL,GH,GM
!
      REAL,DIMENSION(KTS:KTE+1) :: ZHK
!
      REAL,DIMENSION(ITS:ITE,JTS:JTE) :: THSK
!
      REAL,DIMENSION(KTS:KTE,ITS:ITE) :: RHOK
!
      REAL,DIMENSION(ITS:ITE,KTS:KTE,JTS:JTE) :: APE,THE
!
      REAL,DIMENSION(ITS:ITE,KTS:KTE-1,JTS:JTE) :: AKH,AKM
!
      REAL,DIMENSION(ITS:ITE,KTS:KTE+1,JTS:JTE) :: ZINT
!
!***  Begin debugging
      REAL :: ZSL_DIAG
      INTEGER :: IMD,JMD,PRINT_DIAG
!***  End debugging
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
!***  Begin debugging
      IMD=(IMS+IME)/2
      JMD=(JMS+JME)/2
!***  End debugging
!
!***  MAKE PREPARATIONS
!
!----------------------------------------------------------------------
      DTTURBL=DT*STEPBL
      RDTTURBL=1./DTTURBL
      DTDIF=DTTURBL
      RG=1./G
!
      DO J=JTS,JTE
      DO K=KTS,KTE-1
      DO I=ITS,ITE
        AKM(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO K=KTS,KTE+1
      DO I=ITS,ITE
        ZINT(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        ZINT(I,KTE+1,J)=HT(I,J)     ! Z at bottom of lowest sigma layer
!
!!!!!!!!!
!!!!!! UNCOMMENT THESE LINES IF USING ETA COORDINATES
!!!!!!!!!
!!!!!!  ZINT(I,KTE+1,J)=1.E-4       ! Z of bottom of lowest eta layer
!!!!!!  ZHK(KTE+1)=1.E-4            ! Z of bottom of lowest eta layer
!
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO K=KTE,KTS,-1
        KFLIP=KTE+1-K
        DO I=ITS,ITE
          ZINT(I,K,J)=ZINT(I,K+1,J)+DZ(I,KFLIP,J)
          APEX=1./EXNER(I,K,J)
          APE(I,K,J)=APEX
          TX=T(I,K,J)
          THE(I,K,J)=(CWM(I,K,J)*(-ELOCP/TX)+1.)*TH(I,K,J)
        ENDDO
      ENDDO
      ENDDO
!
      EL_MYJ = 0.
!
!----------------------------------------------------------------------
      setup_integration:  DO J=JTS,JTE
!----------------------------------------------------------------------
!
        DO I=ITS,ITE
!
!***  LOWEST LAYER ABOVE GROUND MUST BE FLIPPED
!
          LMH=KTE-LOWLYR(I,J)+1
!
          PTOP=PINT(I,KTE+1,J)      ! KTE+1=KME
          PSFC=PINT(I,LOWLYR(I,J),J)
!
!***  CONVERT LAND MASK (1 FOR SEA; 0 FOR LAND)
!
          SEAMASK=XLAND(I,J)-1.
!
!***  FILL 1-D VERTICAL ARRAYS
!***  AND FLIP DIRECTION SINCE MYJ SCHEME
!***  COUNTS DOWNWARD FROM THE DOMAINS TOP
!
          DO K=KTE,KTS,-1
            KFLIP=KTE+1-K
            TK(K)=T(I,KFLIP,J)
            THEK(K)=THE(I,KFLIP,J)
            RATIOMX=QV(I,KFLIP,J)
            QK(K)=RATIOMX/(1.+RATIOMX)
            CWMK(K)=CWM(I,KFLIP,J)
            PK(K)=PMID(I,KFLIP,J)
            UK(K)=U(I,KFLIP,J)
            VK(K)=V(I,KFLIP,J)
!
!***  TKE=0.5*(q**2) ==> q**2=2.*TKE
!
            Q2K(K)=2.*TKE_MYJ(I,KFLIP,J)
!
!***  COMPUTE THE HEIGHTS OF THE LAYER INTERFACES
!
            ZHK(K)=ZINT(I,K,J)
!
          ENDDO
          ZHK(KTE+1)=HT(I,J)          ! Z at bottom of lowest sigma layer
!
!***  Begin debugging
!         IF(I==IMD.AND.J==JMD)THEN
!           PRINT_DIAG=1
!         ELSE
!           PRINT_DIAG=0
!         ENDIF
!         IF(I==227.AND.J==363)PRINT_DIAG=2
!***  End debugging
!
!----------------------------------------------------------------------
!***
!***  FIND THE MIXING LENGTH
!***
          CALL MIXLEN(LMH,UK,VK,TK,THEK,QK,CWMK                        &
     &               ,Q2K,ZHK,GM,GH,EL                                 &
     &               ,PBLH(I,J),LPBL(I,J),LMXL,CT(I,J)                 &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &               ,IMS,IME,JMS,JME,KMS,KME                          &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE)
!
!----------------------------------------------------------------------
!***
!***  SOLVE FOR THE PRODUCTION/DISSIPATION OF
!***  THE TURBULENT KINETIC ENERGY
!***
!
          CALL PRODQ2(LMH,DTTURBL,USTAR(I,J),GM,GH,EL,Q2K              &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &               ,IMS,IME,JMS,JME,KMS,KME                          &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE)
!
!----------------------------------------------------------------------
!*** THE MODEL LAYER (COUNTING UPWARD) CONTAINING THE TOP OF THE PBL
!----------------------------------------------------------------------
!
          KPBL(I,J)=KTE-LPBL(I,J)+1
!
!----------------------------------------------------------------------
!***
!***  FIND THE EXCHANGE COEFFICIENTS IN THE FREE ATMOSPHERE
!***
          CALL DIFCOF(LMH,LMXL,GM,GH,EL,TK,Q2K,ZHK,AKMK,AKHK      &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &               ,IMS,IME,JMS,JME,KMS,KME                          &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE,PRINT_DIAG)   ! debug
!
!***  COUNTING DOWNWARD FROM THE TOP, THE EXCHANGE COEFFICIENTS AKH 
!***  ARE DEFINED ON THE BOTTOMS OF THE LAYERS KTS TO KTE-1.  COUNTING 
!***  COUNTING UPWARD FROM THE BOTTOM, THOSE SAME COEFFICIENTS EXCH_H
!***  ARE DEFINED ON THE TOPS OF THE LAYERS KTS TO KTE-1.
!
          DO K=KTS,KTE-1
            KFLIP=KTE-K
            AKH(I,K,J)=AKHK(K)
            AKM(I,K,J)=AKMK(K)
            DELTAZ=0.5*(ZHK(KFLIP)-ZHK(KFLIP+2))
            EXCH_H(I,K,J)=AKHK(KFLIP)*DELTAZ
          ENDDO
!
!----------------------------------------------------------------------
!***
!***  CARRY OUT THE VERTICAL DIFFUSION OF
!***  TURBULENT KINETIC ENERGY
!***
!
          CALL VDIFQ(LMH,DTDIF,Q2K,EL,ZHK                              &
     &              ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &              ,IMS,IME,JMS,JME,KMS,KME                           &
     &              ,ITS,ITE,JTS,JTE,KTS,KTE)
!
!***  SAVE THE NEW TKE AND MIXING LENGTH.
!
          DO K=KTS,KTE
            KFLIP=KTE+1-K
            Q2K(KFLIP)=AMAX1(Q2K(KFLIP),EPSQ2)
            TKE_MYJ(I,K,J)=0.5*Q2K(KFLIP)
            IF(K<KTE)EL_MYJ(I,K,J)=EL(K)   ! EL IS NOT DEFINED AT KTE
          ENDDO
!
        ENDDO
!
!----------------------------------------------------------------------
      ENDDO setup_integration
!----------------------------------------------------------------------
!
!***  CONVERT SURFACE SENSIBLE TEMPERATURE TO POTENTIAL TEMPERATURE.
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        PSFC=PINT(I,LOWLYR(I,J),J)
        THSK(I,J)=TSK(I,J)*(1.E5/PSFC)**CAPA
      ENDDO
      ENDDO
!
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
      main_integration:  DO J=JTS,JTE
!----------------------------------------------------------------------
!
        DO I=ITS,ITE
!
!***  FILL 1-D VERTICAL ARRAYS
!***  AND FLIP DIRECTION SINCE MYJ SCHEME
!***  COUNTS DOWNWARD FROM THE DOMAINS TOP
!
          DO K=KTE,KTS,-1
            KFLIP=KTE+1-K
            THEK(K)=THE(I,KFLIP,J)
            RATIOMX=QV(I,KFLIP,J)
            QK(K)=RATIOMX/(1.+RATIOMX)
            CWMK(K)=CWM(I,KFLIP,J)
            ZHK(K)=ZINT(I,K,J)
            RHOK(K,I)=PMID(I,KFLIP,J)/(R_D*T(I,KFLIP,J)*               &
     &                                (1.+P608*QK(K)-CWMK(K)))
          ENDDO
!
!***  COUNTING DOWNWARD FROM THE TOP, THE EXCHANGE COEFFICIENTS AKH
!***  ARE DEFINED ON THE BOTTOMS OF THE LAYERS KTS TO KTE-1.  THESE COEFFICIENTS
!***  ARE ALSO MULTIPLIED BY THE DENSITY AT THE BOTTOM INTERFACE LEVEL.
!
          DO K=KTS,KTE-1
            AKHK(K)=AKH(I,K,J)*0.5*(RHOK(K,I)+RHOK(K+1,I))
          ENDDO
!
          ZHK(KTE+1)=ZINT(I,KTE+1,J)
!
          SEAMASK=XLAND(I,J)-1.
          THZ0(I,J)=(1.-SEAMASK)*THSK(I,J)+SEAMASK*THZ0(I,J)
!
          LLOW=LOWLYR(I,J)
          AKHS_DENS=AKHS(I,J)*RHOK(KTE+1-LLOW,I)
!
          IF(SEAMASK<0.5)THEN
            QFC1=XLV*CHKLOWQ(I,J)*AKHS_DENS
!
            IF(SNOW(I,J)>0..OR.SICE(I,J)>0.5)THEN
              QFC1=QFC1*RLIVWV
            ENDIF
!
            IF(QFC1>0.)THEN
              QLOW=QK(KTE+1-LLOW)
              QSFC(I,J)=QLOW+ELFLX(I,J)/QFC1
            ENDIF
!
          ELSE
            PSFC=PINT(I,LOWLYR(I,J),J)
            EXNSFC=(1.E5/PSFC)**CAPA
            QSFC(I,J)=PQ0SEA/PSFC                                      &
     &         *EXP(A2*(THSK(I,J)-A3*EXNSFC)/(THSK(I,J)-A4*EXNSFC))
          ENDIF
!
          QZ0 (I,J)=(1.-SEAMASK)*QSFC(I,J)+SEAMASK*QZ0 (I,J)
!
!***  LOWEST LAYER ABOVE GROUND MUST BE FLIPPED
!
          LMH=KTE-LOWLYR(I,J)+1
!
!----------------------------------------------------------------------
!***  CARRY OUT THE VERTICAL DIFFUSION OF
!***  TEMPERATURE AND WATER VAPOR
!----------------------------------------------------------------------
!
          CALL VDIFH(DTDIF,LMH,THZ0(I,J),QZ0(I,J)                      &
     &              ,AKHS_DENS,CHKLOWQ(I,J),CT(I,J)                    &
     &              ,THEK,QK,CWMK,AKHK,ZHK,RHOK(KTS,I)                 &
     &              ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &              ,IMS,IME,JMS,JME,KMS,KME                           &
     &              ,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!----------------------------------------------------------------------
!***
!***  COMPUTE PRIMARY VARIABLE TENDENCIES
!***
          DO K=KTS,KTE
            KFLIP=KTE+1-K
            THOLD=TH(I,K,J)
            THNEW=THEK(KFLIP)+CWMK(KFLIP)*ELOCP*APE(I,K,J)
            DTDT=(THNEW-THOLD)*RDTTURBL
            QOLD=QV(I,K,J)/(1.+QV(I,K,J))
            DQDT=(QK(KFLIP)-QOLD)*RDTTURBL
            DCDT=(CWMK(KFLIP)-CWM(I,K,J))*RDTTURBL
!
            RTHBLTEN(I,K,J)=DTDT
            RQVBLTEN(I,K,J)=DQDT/(1.-QK(KFLIP))**2
            RQCBLTEN(I,K,J)=DCDT
          ENDDO
!
!*** Begin debugging
!         IF(I==IMD.AND.J==JMD)THEN
!           PRINT_DIAG=0
!         ELSE
!           PRINT_DIAG=0
!         ENDIF
!         IF(I==227.AND.J==363)PRINT_DIAG=0
!*** End debugging
!
        PSFC=.01*PINT(I,LOWLYR(I,J),J)
        ZSL_DIAG=0.5*DZ(I,1,J)
!
!*** Begin debugging
!         IF(PRINT_DIAG==1)THEN
!
!           write(6,"(a, 2i5, 2i3, 2f8.2, f6.2, 2f8.2)") &
!           {turb4 i,j, Kpbl, Kmxl, Psfc, Zsfc, Zsl, Zpbl, Zmxl =  &
!           , i, j, KPBL(i,j), KTE-LMXL+1, PSFC, ZHK(LMH+1), ZSL_diag  &
!           , PBLH(i,j), ZHK(LMXL)-ZHK(LMH+1)
!           write(6,"(a, 2f7.2, f7.3, 3e11.4)") &
!           {turb4 tsk, thsk, qz0, q**2_0, akhs, exch_0 =  &
!           , tsk(i,j)-273.15, thsk(i,j), 1000.*qz0(i,j) &
!           , 2.*tke_myj(i,1,j), akhs(i,j), akhs(i,j)*ZSL_diag
!           write(6,"(a)") &
!           {turb5 k, Pmid, Pint_1, Tc, TH, DTH, GH, GM, EL, Q**2, Akh, EXCH_h, Dz, Dp
!           do k=kts,kte/2
!             KFLIP=KTE-K   !-- Includes the KFLIP-1 in earlier versions
!             write(6,"(a,i3, 2f8.2, 2f8.3, 3e12.4, 4e11.4, f7.2, f6.2)") &
!            {turb5 , k, .01*pmid(i,k,j),.01*pint(i,k,j), T(i,k,j)-273.15 &
!            , th(i,k,j), DTTURBL*rthblten(i,k,j), GH(KFLIP), GM(KFLIP) &
!            , el_myj(i,KFLIP,j), 2.*tke_myj(i,k+1,j), akh(i,KFLIP,j) &
!            , exch_h(i,k,j), dz(i,k,j), .01*(pint(i,k,j)-pint(i,k+1,j))
!           enddo
!
!         ELSEIF(PRINT_DIAG==2)THEN
!
!           write(6,"(a, 2i5, 2i3, 2f8.2, f6.2, 2f8.2)") &
!           }turb4 i,j, Kpbl, Kmxl, Psfc, Zsfc, Zsl, Zpbl, Zmxl =  &
!           , i, j, KPBL(i,j), KTE-LMXL+1, PSFC, ZHK(LMH+1), ZSL_diag  &
!           , PBLH(i,j), ZHK(LMXL)-ZHK(LMH+1)
!           write(6,"(a, 2f7.2, f7.3, 3e11.4)") &
!           }turb4 tsk, thsk, qz0, q**2_0, akhs, exch_0 =  &
!           , tsk(i,j)-273.15, thsk(i,j), 1000.*qz0(i,j) &
!           , 2.*tke_myj(i,1,j), akhs(i,j), akhs(i,j)*ZSL_diag
!           write(6,"(a)") &
!           }turb5 k, Pmid, Pint_1, Tc, TH, DTH, GH, GM, EL, Q**2, Akh, EXCH_h, Dz, Dp
!           do k=kts,kte/2
!             KFLIP=KTE-K   !-- Includes the KFLIP-1 in earlier versions
!             write(6,"(a,i3, 2f8.2, 2f8.3, 3e12.4, 4e11.4, f7.2, f6.2)") &
!            }turb5 , k, .01*pmid(i,k,j),.01*pint(i,k,j), T(i,k,j)-273.15 &
!            , th(i,k,j), DTTURBL*rthblten(i,k,j), GH(KFLIP), GM(KFLIP) &
!            , el_myj(i,KFLIP,j), 2.*tke_myj(i,k+1,j), akh(i,KFLIP,j) &
!            , exch_h(i,k,j), dz(i,k,j), .01*(pint(i,k,j)-pint(i,k+1,j))
!           enddo
!         ENDIF
!*** End debugging
!
!----------------------------------------------------------------------
        ENDDO
!----------------------------------------------------------------------
        DO I=ITS,ITE
!
!***  FILL 1-D VERTICAL ARRAYS
!***  AND FLIP DIRECTION SINCE MYJ SCHEME
!***  COUNTS DOWNWARD FROM THE DOMAINS TOP
!
          DO K=KTS,KTE-1
            AKMK(K)=AKM(I,K,J)
            AKMK(K)=AKMK(K)*(RHOK(K,I)+RHOK(K+1,I))*0.5
          ENDDO
!
          LLOW=LOWLYR(I,J)
          AKMS_DENS=AKMS(I,J)*RHOK(KTE+1-LLOW,I)
!
          DO K=KTE,KTS,-1
            KFLIP=KTE+1-K
            UK(K)=U(I,KFLIP,J)
            VK(K)=V(I,KFLIP,J)
            ZHK(K)=ZINT(I,K,J)
          ENDDO
          ZHK(KTE+1)=ZINT(I,KTE+1,J)
!
!----------------------------------------------------------------------
!***  CARRY OUT THE VERTICAL DIFFUSION OF
!***  VELOCITY COMPONENTS
!----------------------------------------------------------------------
!
          CALL VDIFV(LMH,DTDIF,UZ0(I,J),VZ0(I,J)                       &
     &              ,AKMS_DENS,UK,VK,AKMK,ZHK,RHOK(KTS,I)              &
     &              ,IDS,IDE,JDS,JDE,KDS,KDE                           &
     &              ,IMS,IME,JMS,JME,KMS,KME                           &
     &              ,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!
!----------------------------------------------------------------------
!***
!***  COMPUTE PRIMARY VARIABLE TENDENCIES
!***
          DO K=KTS,KTE
            KFLIP=KTE+1-K
            DUDT=(UK(KFLIP)-U(I,K,J))*RDTTURBL
            DVDT=(VK(KFLIP)-V(I,K,J))*RDTTURBL
            RUBLTEN(I,K,J)=DUDT
            RVBLTEN(I,K,J)=DVDT
          ENDDO
!
        ENDDO
!----------------------------------------------------------------------
!
      ENDDO main_integration
!
!----------------------------------------------------------------------
!
      END SUBROUTINE MYJPBL
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                          SUBROUTINE MIXLEN                            &
!----------------------------------------------------------------------
!   ******************************************************************
!   *                                                                *
!   *                   LEVEL 2.5 MIXING LENGTH                      *
!   *                                                                *
!   ******************************************************************
!
     &(LMH,U,V,T,THE,Q,CWM,Q2,Z,GM,GH,EL,PBLH,LPBL,LMXL,CT             &
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: LMH
!
      INTEGER,INTENT(OUT) :: LMXL,LPBL
!
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: CWM,Q,Q2,T,THE,U,V
!
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,INTENT(OUT) :: PBLH
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(OUT) :: EL,GH,GM
!
      REAL,INTENT(INOUT) :: CT
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K,LPBLM
!
      REAL :: A,ADEN,B,BDEN,AUBR,BUBR,BLMX,EL0,ELOQ2X,GHL,GML           &
     &       ,QOL2ST,QOL2UN,QDZL,RDZ,SQ,SREL,SZQ,TEM,THM,VKRMZ
!
      REAL,DIMENSION(KTS:KTE) :: Q1
!
      REAL,DIMENSION(KTS:KTE-1) :: DTH,ELM,REL
!
!----------------------------------------------------------------------
!**********************************************************************
!--------------FIND THE HEIGHT OF THE PBL-------------------------------
      LPBL=LMH
!
      DO K=LMH-1,1,-1
        IF(Q2(K)<=EPSQ2*FH)THEN
          LPBL=K
          GO TO 110
        ENDIF
      ENDDO
!
      LPBL=1
!
!--------------THE HEIGHT OF THE PBL------------------------------------
!
 110  PBLH=Z(LPBL)-Z(LMH+1)
!
!-----------------------------------------------------------------------
      DO K=KTS,LMH
        Q1(K)=0.
      ENDDO
!
      DO K=1,LMH-1
        DTH(K)=THE(K)-THE(K+1)
      ENDDO
!
      DO K=LMH-2,1,-1
        IF(DTH(K)>0..AND.DTH(K+1)<=0.)THEN
          DTH(K)=DTH(K)+CT
          EXIT
        ENDIF
      ENDDO
!
      CT=0.
!----------------------------------------------------------------------
      DO K=KTS,LMH-1
        RDZ=2./(Z(K)-Z(K+2))
        GML=((U(K)-U(K+1))**2+(V(K)-V(K+1))**2)*RDZ*RDZ
        GM(K)=MAX(GML,EPSGM)
!
        TEM=(T(K)+T(K+1))*0.5
        THM=(THE(K)+THE(K+1))*0.5
!
        A=THM*P608
        B=(ELOCP/TEM-1.-P608)*THM
!
        GHL=(DTH(K)*((Q(K)+Q(K+1)+CWM(K)+CWM(K+1))*(0.5*P608)+1.)      &
     &     +(Q(K)-Q(K+1)+CWM(K)-CWM(K+1))*A                            &
     &     +(CWM(K)-CWM(K+1))*B)*RDZ
!
        IF(ABS(GHL)<=EPSGH)GHL=EPSGH
        GH(K)=GHL
      ENDDO
!
!----------------------------------------------------------------------
!***  FIND MAXIMUM MIXING LENGTHS AND THE LEVEL OF THE PBL TOP
!----------------------------------------------------------------------
!
      LMXL=LMH
!
      DO K=KTS,LMH-1
        GML=GM(K)
        GHL=GH(K)
!
        IF(GHL>=EPSGH)THEN
          IF(GML/GHL<=REQU)THEN
            ELM(K)=EPSL
            LMXL=K
          ELSE
            AUBR=(AUBM*GML+AUBH*GHL)*GHL
            BUBR= BUBM*GML+BUBH*GHL
            QOL2ST=(-0.5*BUBR+SQRT(BUBR*BUBR*0.25-AUBR*CUBR))*RCUBR
            ELOQ2X=1./QOL2ST
            ELM(K)=MAX(SQRT(ELOQ2X*Q2(K)),EPSL)
          ENDIF
        ELSE
          ADEN=(ADNM*GML+ADNH*GHL)*GHL
          BDEN= BDNM*GML+BDNH*GHL
          QOL2UN=-0.5*BDEN+SQRT(BDEN*BDEN*0.25-ADEN)
          ELOQ2X=1./(QOL2UN+EPSRU)       ! repsr1/qol2un
          ELM(K)=MAX(SQRT(ELOQ2X*Q2(K)),EPSL)
        ENDIF
      ENDDO
!
      IF(ELM(LMH-1)==EPSL)LMXL=LMH
!
!----------------------------------------------------------------------
!***  THE HEIGHT OF THE MIXED LAYER
!----------------------------------------------------------------------
!
      BLMX=Z(LMXL)-Z(LMH+1)
!
!----------------------------------------------------------------------
      DO K=LPBL,LMH
        Q1(K)=SQRT(Q2(K))
      ENDDO
!----------------------------------------------------------------------
      SZQ=0.
      SQ =0.
!
      DO K=KTS,LMH-1
        QDZL=(Q1(K)+Q1(K+1))*(Z(K+1)-Z(K+2))
        SZQ=(Z(K+1)+Z(K+2)-Z(LMH+1)-Z(LMH+1))*QDZL+SZQ
        SQ=QDZL+SQ
      ENDDO
!
!----------------------------------------------------------------------
!***  COMPUTATION OF ASYMPTOTIC L IN BLACKADAR FORMULA
!----------------------------------------------------------------------
!
      EL0=MIN(ALPH*SZQ*0.5/SQ,EL0MAX)
      EL0=MAX(EL0            ,EL0MIN)
!
!----------------------------------------------------------------------
!***  ABOVE THE PBL TOP
!----------------------------------------------------------------------
!
      LPBLM=MAX(LPBL-1,1)
!
      DO K=KTS,LPBLM
        EL(K)=MIN((Z(K)-Z(K+2))*ELFC,ELM(K))
        REL(K)=EL(K)/ELM(K)
      ENDDO
!
!----------------------------------------------------------------------
!***  INSIDE THE PBL
!----------------------------------------------------------------------
!
      IF(LPBL<LMH)THEN
        DO K=LPBL,LMH-1
          VKRMZ=(Z(K+1)-Z(LMH+1))*VKARMAN
          EL(K)=MIN(VKRMZ/(VKRMZ/EL0+1.),ELM(K))
          REL(K)=EL(K)/ELM(K)
        ENDDO
      ENDIF
!
      DO K=LPBL+1,LMH-2
        SREL=MIN(((REL(K-1)+REL(K+1))*0.5+REL(K))*0.5,REL(K))
        EL(K)=MAX(SREL*ELM(K),EPSL)
      ENDDO
!
!----------------------------------------------------------------------
      END SUBROUTINE MIXLEN
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                          SUBROUTINE PRODQ2                            &
!----------------------------------------------------------------------
!   ******************************************************************
!   *                                                                *
!   *            LEVEL 2.5 Q2 PRODUCTION/DISSIPATION                 *
!   *                                                                *
!   ******************************************************************
!
     &(LMH,DTTURBL,USTAR,GM,GH,EL,Q2                                   &
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: LMH
!
      REAL,INTENT(IN) :: DTTURBL,USTAR
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(IN) :: GH,GM
      REAL,DIMENSION(KTS:KTE-1),INTENT(INOUT) :: EL
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: Q2
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K
!
      REAL :: ADEN,AEQU,ANUM,ARHS,BDEN,BEQU,BNUM,BRHS,CDEN,CRHS        &
     &       ,DLOQ1,ELOQ11,ELOQ12,ELOQ13,ELOQ21,ELOQ22,ELOQ31,ELOQ32   &
     &       ,ELOQ41,ELOQ42,ELOQ51,ELOQ52,ELOQN,EQOL2,GHL,GML          &
     &       ,RDEN1,RDEN2,RHS2,RHSP1,RHSP2,RHST2
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
      main_integration: DO K=1,LMH-1
        GML=GM(K)
        GHL=GH(K)
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE EQUILIBRIUM EQUATION
!----------------------------------------------------------------------
!
        AEQU=(AEQM*GML+AEQH*GHL)*GHL
        BEQU= BEQM*GML+BEQH*GHL
!
!----------------------------------------------------------------------
!***  EQUILIBRIUM SOLUTION FOR L/Q
!----------------------------------------------------------------------
!
        EQOL2=-0.5*BEQU+SQRT(BEQU*BEQU*0.25-AEQU)
!
!----------------------------------------------------------------------
!***  IS THERE PRODUCTION/DISSIPATION ?
!----------------------------------------------------------------------
!
        IF((GML+GHL*GHL<=EPSTRB)                                       &
     &   .OR.(GHL>=EPSGH.AND.GML/GHL<=REQU)                            &
     &   .OR.(EQOL2<=EPS2))THEN
!
!----------------------------------------------------------------------
!***  NO TURBULENCE
!----------------------------------------------------------------------
!
          Q2(K)=EPSQ2
          EL(K)=EPSL
!----------------------------------------------------------------------
!
        ELSE
!
!----------------------------------------------------------------------
!***  TURBULENCE
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE TERMS IN THE NUMERATOR
!----------------------------------------------------------------------
!
          ANUM=(ANMM*GML+ANMH*GHL)*GHL
          BNUM= BNMM*GML+BNMH*GHL
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!----------------------------------------------------------------------
!
          ADEN=(ADNM*GML+ADNH*GHL)*GHL
          BDEN= BDNM*GML+BDNH*GHL
          CDEN= 1.
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE NUMERATOR OF THE LINEARIZED EQ.
!----------------------------------------------------------------------
!
          ARHS=-(ANUM*BDEN-BNUM*ADEN)*2.
          BRHS=- ANUM*4.
          CRHS=- BNUM*2.
!
!----------------------------------------------------------------------
!***  INITIAL VALUE OF L/Q
!----------------------------------------------------------------------
!
          DLOQ1=EL(K)/SQRT(Q2(K))
!
!----------------------------------------------------------------------
!***  FIRST ITERATION FOR L/Q, RHS=0
!----------------------------------------------------------------------
!
          ELOQ21=1./EQOL2
          ELOQ11=SQRT(ELOQ21)
          ELOQ31=ELOQ21*ELOQ11
          ELOQ41=ELOQ21*ELOQ21
          ELOQ51=ELOQ21*ELOQ31
!
!----------------------------------------------------------------------
!***  1./DENOMINATOR
!----------------------------------------------------------------------
!
          RDEN1=1./(ADEN*ELOQ41+BDEN*ELOQ21+CDEN)
!
!----------------------------------------------------------------------
!***  D(RHS)/D(L/Q)
!----------------------------------------------------------------------
!
          RHSP1=(ARHS*ELOQ51+BRHS*ELOQ31+CRHS*ELOQ11)*RDEN1*RDEN1
!
!----------------------------------------------------------------------
!***  FIRST-GUESS SOLUTION
!----------------------------------------------------------------------
!
          ELOQ12=ELOQ11+(DLOQ1-ELOQ11)*EXP(RHSP1*DTTURBL)
          ELOQ12=MAX(ELOQ12,EPS1)
!
!----------------------------------------------------------------------
!***  SECOND ITERATION FOR L/Q
!----------------------------------------------------------------------
!
          ELOQ22=ELOQ12*ELOQ12
          ELOQ32=ELOQ22*ELOQ12
          ELOQ42=ELOQ22*ELOQ22
          ELOQ52=ELOQ22*ELOQ32
!
!----------------------------------------------------------------------
!***  1./DENOMINATOR
!----------------------------------------------------------------------
!
          RDEN2=1./(ADEN*ELOQ42+BDEN*ELOQ22+CDEN)
          RHS2 =-(ANUM*ELOQ42+BNUM*ELOQ22)*RDEN2+RB1
          RHSP2= (ARHS*ELOQ52+BRHS*ELOQ32+CRHS*ELOQ12)*RDEN2*RDEN2
          RHST2=RHS2/RHSP2
!
!----------------------------------------------------------------------
!***  CORRECTED SOLUTION
!----------------------------------------------------------------------
!
          ELOQ13=ELOQ12-RHST2+(RHST2+DLOQ1-ELOQ12)*EXP(RHSP2*DTTURBL)
          ELOQ13=AMAX1(ELOQ13,EPS1)
!
!----------------------------------------------------------------------
!***  TWO ITERATIONS IS ENOUGH IN MOST CASES ...
!----------------------------------------------------------------------
!
          ELOQN=ELOQ13
!
          IF(ELOQN>EPS1)THEN
            Q2(K)=EL(K)*EL(K)/(ELOQN*ELOQN)
            Q2(K)=AMAX1(Q2(K),EPSQ2)
            IF(Q2(K)==EPSQ2)THEN
              EL(K)=EPSL
            ENDIF
          ELSE
            Q2(K)=EPSQ2
            EL(K)=EPSL
          ENDIF
!
!----------------------------------------------------------------------
!***  END OF TURBULENT BRANCH
!----------------------------------------------------------------------
!
        ENDIF
!----------------------------------------------------------------------
!***  END OF PRODUCTION/DISSIPATION LOOP
!----------------------------------------------------------------------
!
      ENDDO main_integration
!
!----------------------------------------------------------------------
!***  LOWER BOUNDARY CONDITION FOR Q2
!----------------------------------------------------------------------
!
      Q2(LMH)=AMAX1(B1**(2./3.)*USTAR*USTAR,EPSQ2)
!----------------------------------------------------------------------
!
      END SUBROUTINE PRODQ2
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                           SUBROUTINE DIFCOF                           &
!   ******************************************************************
!   *                                                                *
!   *                LEVEL 2.5 DIFFUSION COEFFICIENTS                *
!   *                                                                *
!   ******************************************************************
     &(LMH,LMXL,GM,GH,EL,T,Q2,Z,AKM,AKH                                &
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE,PRINT_DIAG)   ! debug
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: LMH,LMXL
!
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: Q2,T
      REAL,DIMENSION(KTS:KTE-1),INTENT(IN) :: EL,GH,GM
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(OUT) :: AKH,AKM
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K,KINV
!
      REAL :: ADEN,AKMIN,BDEN,BESH,BESM,CDEN,D2T,ELL,ELOQ2,ELOQ4,ELQDZ &
     &       ,ESH,ESM,GHL,GML,Q1L,RDEN,RDZ
!
!*** Begin debugging
      INTEGER,INTENT(IN) :: PRINT_DIAG
!     REAL :: D2Tmin
!*** End debugging
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
      DO K=1,LMH-1
        ELL=EL(K)
!
        ELOQ2=ELL*ELL/Q2(K)
        ELOQ4=ELOQ2*ELOQ2
!
        GML=GM(K)
        GHL=GH(K)
!
!----------------------------------------------------------------------
!***  COEFFICIENTS OF THE TERMS IN THE DENOMINATOR
!----------------------------------------------------------------------
!
        ADEN=(ADNM*GML+ADNH*GHL)*GHL
        BDEN= BDNM*GML+BDNH*GHL
        CDEN= 1.
!
!----------------------------------------------------------------------
!***  COEFFICIENTS FOR THE SM DETERMINANT
!----------------------------------------------------------------------
!
        BESM=BSMH*GHL
!
!----------------------------------------------------------------------
!***  COEFFICIENTS FOR THE SH DETERMINANT
!----------------------------------------------------------------------
!
        BESH=BSHM*GML+BSHH*GHL
!
!----------------------------------------------------------------------
!***  1./DENOMINATOR
!----------------------------------------------------------------------
!
        RDEN=1./(ADEN*ELOQ4+BDEN*ELOQ2+CDEN)
!
!----------------------------------------------------------------------
!***  SM AND SH
!----------------------------------------------------------------------
!
        ESM=(BESM*ELOQ2+CESM)*RDEN
        ESH=(BESH*ELOQ2+CESH)*RDEN
!
!----------------------------------------------------------------------
!***  DIFFUSION COEFFICIENTS
!----------------------------------------------------------------------
!
        RDZ=2./(Z(K)-Z(K+2))
        Q1L=SQRT(Q2(K))
        ELQDZ=ELL*Q1L*RDZ
        AKM(K)=ELQDZ*ESM
        AKH(K)=ELQDZ*ESH
!----------------------------------------------------------------------
      ENDDO
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
!***  INVERSIONS
!----------------------------------------------------------------------
!
!     IF(LMXL==LMH)THEN
!       KINV=LMH
!       D2Tmin=0.
!
!       DO K=LMH/2,LMH-1
!         D2T=T(K-1)-2.*T(K)+T(K+1)
!         IF(D2T<D2Tmin)THEN
!           D2Tmin=D2T
!           IF(D2T<0)KINV=K
!         ENDIF
!       ENDDO
!
!       IF(KINV<LMH)THEN
!         DO K=KINV-1,LMH-1
!           RDZ=2./(Z(K)-Z(K+2))
!           AKMIN=0.5*RDZ
!           AKM(K)=MAX(AKM(K),AKMIN)
!           AKH(K)=MAX(AKH(K),AKMIN)
!         ENDDO
!
!*** Begin debugging
!         IF(PRINT_DIAG>0)THEN
!           write(6,"(a,3i3)") {turb1 lmxl,lmh,kinv=,lmxl,lmh,kinv
!           write(6,"(a,3i3)") }turb1 lmxl,lmh,kinv=,lmxl,lmh,kinv
!           IF(PRINT_DIAG==1)THEN
!             write(6,"(a)") &
!               {turb3 k, t, d2t, rdz, z(k), z(k+2), akmin, akh 
!           ELSE
!             write(6,"(a)") &
!               }turb3 k, t, d2t, rdz, z(k), z(k+2), akmin, akh 
!           ENDIF
!           DO K=LMH-1,KINV-1,-1
!             D2T=T(K-1)-2.*T(K)+T(K+1)
!             RDZ=2./(Z(K)-Z(K+2))
!             AKMIN=0.5*RDZ
!             IF(PRINT_DIAG==1)THEN
!               write(6,"(a,i3,f8.3,2e12.5,2f9.2,2e12.5)") {turb3  &
!               ,k,t(k)-273.15,d2t,rdz,z(k),z(k+2),akmin,akh(k)
!             ELSE
!               write(6,"(a,i3,f8.3,2e12.5,2f9.2,2e12.5)") }turb3  &
!               ,k,t(k)-273.15,d2t,rdz,z(k),z(k+2),akmin,akh(k)
!             ENDIF
!           ENDDO
!         ENDIF     !- IF (print_diag > 0) THEN
!       ENDIF       !- IF(KINV<LMH)THEN
!*** End debugging
!
!     ENDIF
!----------------------------------------------------------------------
!
      END SUBROUTINE DIFCOF
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
                           SUBROUTINE VDIFQ                            &
!   ******************************************************************
!   *                                                                *
!   *               VERTICAL DIFFUSION OF Q2 (TKE)                   *
!   *                                                                *
!   ******************************************************************
     &(LMH,DTDIF,Q2,EL,Z                                               &
     &,IDS,IDE,JDS,JDE,KDS,KDE                                         &
     &,IMS,IME,JMS,JME,KMS,KME                                         &
     &,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: LMH
!
      REAL,INTENT(IN) :: DTDIF
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(IN) :: EL
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: Q2
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K
!
      REAL :: ADEN,AKQS,BDEN,BESH,BESM,CDEN,CF,DTOZS,ELL,ELOQ2,ELOQ4   &
     &       ,ELQDZ,ESH,ESM,ESQHF,GHL,GML,Q1L,RDEN,RDZ
!
      REAL,DIMENSION(KTS:KTE-2) :: AKQ,CM,CR,DTOZ,RSQ2
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!***
!***  VERTICAL TURBULENT DIFFUSION
!***
!----------------------------------------------------------------------
      ESQHF=0.5*ESQ
!
      DO K=KTS,LMH-2
        DTOZ(K)=(DTDIF+DTDIF)/(Z(K)-Z(K+2))
        AKQ(K)=SQRT((Q2(K)+Q2(K+1))*0.5)*(EL(K)+EL(K+1))*ESQHF         &
     &        /(Z(K+1)-Z(K+2))
        CR(K)=-DTOZ(K)*AKQ(K)
      ENDDO
!
      CM(1)=DTOZ(1)*AKQ(1)+1.
      RSQ2(1)=Q2(1)
!
      DO K=KTS+1,LMH-2
        CF=-DTOZ(K)*AKQ(K-1)/CM(K-1)
        CM(K)=-CR(K-1)*CF+(AKQ(K-1)+AKQ(K))*DTOZ(K)+1.
        RSQ2(K)=-RSQ2(K-1)*CF+Q2(K)
      ENDDO
!
      DTOZS=(DTDIF+DTDIF)/(Z(LMH-1)-Z(LMH+1))
      AKQS=SQRT((Q2(LMH-1)+Q2(LMH))*0.5)*(EL(LMH-1)+ELZ0)*ESQHF        &
     &    /(Z(LMH)-Z(LMH+1))
!
      CF=-DTOZS*AKQ(LMH-2)/CM(LMH-2)
!
      Q2(LMH-1)=(DTOZS*AKQS*Q2(LMH)-RSQ2(LMH-2)*CF+Q2(LMH-1))          &
     &        /((AKQ(LMH-2)+AKQS)*DTOZS-CR(LMH-2)*CF+1.)
!
      DO K=LMH-2,KTS,-1
        Q2(K)=(-CR(K)*Q2(K+1)+RSQ2(K))/CM(K)
      ENDDO
!----------------------------------------------------------------------
!
      END SUBROUTINE VDIFQ
!
!----------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!---------------------------------------------------------------------
      SUBROUTINE VDIFH(DTDIF,LMH,THZ0,QZ0,RKHS,CHKLOWQ,CT             &
     &                ,THE,Q,CWM,RKH,Z,RHO                            &
     &                ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                ,IMS,IME,JMS,JME,KMS,KME                        &
     &                ,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!     ***************************************************************
!     *                                                             *
!     *         VERTICAL DIFFUSION OF MASS VARIABLES                *
!     *                                                             *
!     ***************************************************************
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
!---------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE,I,J
!
      INTEGER,INTENT(IN) :: LMH
!
      REAL,INTENT(IN) :: CHKLOWQ,CT,DTDIF,QZ0,RKHS,THZ0
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(IN) :: RKH
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: RHO
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: CWM,Q,THE
!
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K
!
      REAL :: CF,CMB,CMCB,CMQB,CMTB,CTHF,DTOZL,DTOZS                   &
     &       ,RCML,RKHH,RKQS,RSCB,RSQB,RSTB
!
      REAL,DIMENSION(KTS:KTE-1) :: CM,CR,DTOZ,RKCT,RSC,RSQ,RST
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      CTHF=0.5*CT
!
      DO K=KTS,LMH-1
        DTOZ(K)=DTDIF/(Z(K)-Z(K+1))
        CR(K)=-DTOZ(K)*RKH(K)
        RKCT(K)=RKH(K)*(Z(K)-Z(K+2))*CTHF
      ENDDO
!
      CM(KTS)=DTOZ(KTS)*RKH(KTS)+RHO(KTS)
!----------------------------------------------------------------------
      RST(KTS)=-RKCT(KTS)*DTOZ(KTS)                                    &
     &         +THE(KTS)*RHO(KTS)
      RSQ(KTS)=Q(KTS)  *RHO(KTS)
      RSC(KTS)=CWM(KTS)*RHO(KTS)
!----------------------------------------------------------------------
      DO K=KTS+1,LMH-1
        DTOZL=DTOZ(K)
        CF=-DTOZL*RKH(K-1)/CM(K-1)
        CM(K)=-CR(K-1)*CF+(RKH(K-1)+RKH(K))*DTOZL+RHO(K)
        RST(K)=-RST(K-1)*CF+(RKCT(K-1)-RKCT(K))*DTOZL+THE(K)*RHO(K)
        RSQ(K)=-RSQ(K-1)*CF+Q(K)  *RHO(K)
        RSC(K)=-RSC(K-1)*CF+CWM(K)*RHO(K)
      ENDDO
!
      DTOZS=DTDIF/(Z(LMH)-Z(LMH+1))
      RKHH=RKH(LMH-1)
!
      CF=-DTOZS*RKHH/CM(LMH-1)
      RKQS=RKHS*CHKLOWQ
!
      CMB=CR(LMH-1)*CF
      CMTB=-CMB+(RKHH+RKHS)*DTOZS+RHO(LMH)
      CMQB=-CMB+(RKHH+RKQS)*DTOZS+RHO(LMH)
      CMCB=-CMB+(RKHH     )*DTOZS+RHO(LMH)
!
      RSTB=-RST(LMH-1)*CF+RKCT(LMH-1)*DTOZS+THE(LMH)*RHO(LMH)
      RSQB=-RSQ(LMH-1)*CF+Q(LMH)  *RHO(LMH)
      RSCB=-RSC(LMH-1)*CF+CWM(LMH)*RHO(LMH)
!----------------------------------------------------------------------
      THE(LMH)=(DTOZS*RKHS*THZ0+RSTB)/CMTB
      Q(LMH)  =(DTOZS*RKQS*QZ0 +RSQB)/CMQB
      CWM(LMH)=(                RSCB)/CMCB
!----------------------------------------------------------------------
      DO K=LMH-1,KTS,-1
        RCML=1./CM(K)
        THE(K)=(-CR(K)*THE(K+1)+RST(K))*RCML
        Q(K)  =(-CR(K)*  Q(K+1)+RSQ(K))*RCML
        CWM(K)=(-CR(K)*CWM(K+1)+RSC(K))*RCML
      ENDDO
!----------------------------------------------------------------------
!
      END SUBROUTINE VDIFH
!
!---------------------------------------------------------------------
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!---------------------------------------------------------------------
      SUBROUTINE VDIFV(LMH,DTDIF,UZ0,VZ0,RKMS,U,V,RKM,Z,RHO           &
     &                ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                ,IMS,IME,JMS,JME,KMS,KME                        &
                      ,ITS,ITE,JTS,JTE,KTS,KTE,I,J)
!     ***************************************************************
!     *                                                             *
!     *        VERTICAL DIFFUSION OF VELOCITY COMPONENTS            *
!     *                                                             *
!     ***************************************************************
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
!---------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                   &
     &                     ,IMS,IME,JMS,JME,KMS,KME                   &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE,I,J
!
      INTEGER,INTENT(IN) :: LMH
!
      REAL,INTENT(IN) :: RKMS,DTDIF,UZ0,VZ0
!
      REAL,DIMENSION(KTS:KTE-1),INTENT(IN) :: RKM
      REAL,DIMENSION(KTS:KTE),INTENT(IN) :: RHO
      REAL,DIMENSION(KTS:KTE+1),INTENT(IN) :: Z
!
      REAL,DIMENSION(KTS:KTE),INTENT(INOUT) :: U,V
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: K
!
      REAL :: CF,DTOZAK,DTOZL,DTOZS,RCML,RCMVB,RHOK,RKMH
!
      REAL,DIMENSION(KTS:KTE-1) :: CM,CR,DTOZ,RSU,RSV
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      DO K=1,LMH-1
        DTOZ(K)=DTDIF/(Z(K)-Z(K+1))
        CR(K)=-DTOZ(K)*RKM(K)
      ENDDO
!
      RHOK=RHO(1)
      CM(1)=DTOZ(1)*RKM(1)+RHOK
      RSU(1)=U(1)*RHOK
      RSV(1)=V(1)*RHOK
!----------------------------------------------------------------------
      DO K=2,LMH-1
        DTOZL=DTOZ(K)
        CF=-DTOZL*RKM(K-1)/CM(K-1)
        RHOK=RHO(K)
        CM(K)=-CR(K-1)*CF+(RKM(K-1)+RKM(K))*DTOZL+RHOK
        RSU(K)=-RSU(K-1)*CF+U(K)*RHOK
        RSV(K)=-RSV(K-1)*CF+V(K)*RHOK
      ENDDO
!----------------------------------------------------------------------
      DTOZS=DTDIF/(Z(LMH)-Z(LMH+1))
      RKMH=RKM(LMH-1)
!
      CF=-DTOZS*RKMH/CM(LMH-1)
      RHOK=RHO(LMH)
      RCMVB=1./((RKMH+RKMS)*DTOZS-CR(LMH-1)*CF+RHOK)
      DTOZAK=DTOZS*RKMS
!----------------------------------------------------------------------
      U(LMH)=(DTOZAK*UZ0-RSU(LMH-1)*CF+U(LMH)*RHOK)*RCMVB
      V(LMH)=(DTOZAK*VZ0-RSV(LMH-1)*CF+V(LMH)*RHOK)*RCMVB
!----------------------------------------------------------------------
      DO K=LMH-1,1,-1
        RCML=1./CM(K)
        U(K)=(-CR(K)*U(K+1)+RSU(K))*RCML
        V(K)=(-CR(K)*V(K+1)+RSV(K))*RCML
      ENDDO
!----------------------------------------------------------------------
!
      END SUBROUTINE VDIFV
!
!-----------------------------------------------------------------------
!
!=======================================================================
      SUBROUTINE MYJPBLINIT(RUBLTEN,RVBLTEN,RTHBLTEN,RQVBLTEN,          &
     &                      TKE_MYJ,EXCH_H,RESTART,ALLOWED_TO_READ,     &
     &                      IDS,IDE,JDS,JDE,KDS,KDE,                    &
     &                      IMS,IME,JMS,JME,KMS,KME,                    &
     &                      ITS,ITE,JTS,JTE,KTS,KTE                 )
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      LOGICAL,INTENT(IN) :: ALLOWED_TO_READ,RESTART
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE,                    &
     &                      IMS,IME,JMS,JME,KMS,KME,                    &
     &                      ITS,ITE,JTS,JTE,KTS,KTE

      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(OUT) ::    EXCH_H, &
     &                                                         RUBLTEN, &
     &                                                         RVBLTEN, &
     &                                                        RTHBLTEN, &
     &                                                        RQVBLTEN, &
     &                                                         TKE_MYJ
      INTEGER :: I,J,K,ITF,JTF,KTF
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      JTF=MIN0(JTE,JDE-1)
      KTF=MIN0(KTE,KDE-1)
      ITF=MIN0(ITE,IDE-1)

      IF(.NOT.RESTART)THEN
        DO J=JTS,JTF
        DO K=KTS,KTF
        DO I=ITS,ITF
          TKE_MYJ(I,K,J)=EPSQ2
          RUBLTEN(I,K,J)=0.
          RVBLTEN(I,K,J)=0.
          RTHBLTEN(I,K,J)=0.
          RQVBLTEN(I,K,J)=0.
          EXCH_H(I,K,J)=0.
        ENDDO
        ENDDO
        ENDDO
      ENDIF

      END SUBROUTINE MYJPBLINIT
!-----------------------------------------------------------------------
!
      END MODULE MODULE_BL_MYJPBL
!
!-----------------------------------------------------------------------