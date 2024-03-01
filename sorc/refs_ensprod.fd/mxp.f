
C===================================================== PELGEV.FOR
      SUBROUTINE PELGEV(XMOM,PARA,*)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  PARAMETER ESTIMATION VIA L-MOMENTS FOR THE GENERALIZED EXTREME-VALUE
C  DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  XMOM   * INPUT* ARRAY OF LENGTH 3. CONTAINS THE L-MOMENTS LAMBDA-1,
C                  LAMBDA-2, TAU-3.
C  PARA   *OUTPUT* ARRAY OF LENGTH 3. ON EXIT, CONTAINS THE PARAMETERS
C                  IN THE ORDER XI, ALPHA, K (LOCATION, SCALE, SHAPE).
C
C  OTHER ROUTINES USED: DLGAMA
C
C  METHOD: FOR  -0.8 LE TAU3 LT 1,  K IS APPROXIMATED BY RATIONAL
C  FUNCTIONS AS IN DONALDSON (1996, COMMUN. STATIST. SIMUL. COMPUT.).
C  IF TAU3 IS OUTSIDE THIS RANGE, NEWTON-RAPHSON ITERATION IS USED.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMOM(3),PARA(3)
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/,THREE/3D0/
      DATA P8/0.8D0/,P97/0.97D0/
C
C         SMALL IS USED TO TEST WHETHER K IS EFFECTIVELY ZERO
C         EPS,MAXIT CONTROL THE TEST FOR CONVERGENCE OF N-R ITERATION
C
      DATA SMALL/1D-5/,EPS/1D-6/,MAXIT/20/
C
C         EU IS EULER'S CONSTANT
C         DL2 IS LOG(2), DL3 IS LOG(3)
C
      DATA EU/0.57721566D0/,DL2/0.69314718D0/,DL3/1.0986123D0/
C
C         COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATIONS FOR K
C
      DATA A0,A1,A2/ 0.28377530D0,-1.21096399D0,-2.50728214D0/
      DATA A3,A4   /-1.13455566D0,-0.07138022D0/
      DATA B1,B2,B3/ 2.06189696D0, 1.31912239D0, 0.25077104D0/
      DATA C1,C2,C3/ 1.59921491D0,-0.48832213D0, 0.01573152D0/
      DATA D1,D2   /-0.64363929D0, 0.08985247D0/
C
      T3=XMOM(3)
      IF(XMOM(2).LE.ZERO)GOTO 1000
      IF(DABS(T3).GE.ONE)GOTO 1000
      IF(T3.LE.ZERO)GOTO 10
C
C         RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN 0 AND 1
C
      Z=ONE-T3
      G=(-ONE+Z*(C1+Z*(C2+Z*C3)))/(ONE+Z*(D1+Z*D2))
      IF(DABS(G).LT.SMALL)GOTO 50
      GOTO 40
C
C         RATIONAL-FUNCTION APPROXIMATION FOR TAU3 BETWEEN -0.8 AND 0
C
   10 G=(A0+T3*(A1+T3*(A2+T3*(A3+T3*A4))))/(ONE+T3*(B1+T3*(B2+T3*B3)))
      IF(T3.GE.-P8)GOTO 40
C
C         NEWTON-RAPHSON ITERATION FOR TAU3 LESS THAN -0.8
C
      IF(T3.LE.-P97)G=ONE-DLOG(ONE+T3)/DL2
      T0=(T3+THREE)*HALF
      DO 20 IT=1,MAXIT
      X2=TWO**(-G)
      X3=THREE**(-G)
      XX2=ONE-X2
      XX3=ONE-X3
      T=XX3/XX2
      DERIV=(XX2*X3*DL3-XX3*X2*DL2)/(XX2*XX2)
      GOLD=G
      G=G-(T-T0)/DERIV
      IF(DABS(G-GOLD).LE.EPS*G)GOTO 30
   20 CONTINUE
      WRITE(6,7010)
   30 CONTINUE
C
C         ESTIMATE ALPHA,XI
C
   40 PARA(3)=G
c     GAM=DEXP(DLGAMA(ONE+G))   !use double precision log Gamma
      GAM=GAMMA(ONE+G)          !use Gamma
      PARA(2)=XMOM(2)*G/(GAM*(ONE-TWO**(-G)))
      PARA(1)=XMOM(1)-PARA(2)*(ONE-GAM)/G
      RETURN
C
C         ESTIMATED K EFFECTIVELY ZERO
C
   50 PARA(3)=ZERO
      PARA(2)=XMOM(2)/DL2
      PARA(1)=XMOM(1)-EU*PARA(2)
      RETURN
C
 1000 WRITE(6,7000)
      RETURN 1
C
 7000 FORMAT(' *** ERROR *** ROUTINE PELGEV : L-MOMENTS INVALID')
 7010 FORMAT(' ** WARNING ** ROUTINE PELGEV :',
     *  ' ITERATION HAS NOT CONVERGED. RESULTS MAY BE UNRELIABLE.')
      END


C===================================================== QUAGEV.FOR
C     DOUBLE PRECISION FUNCTION QUAGEV(F,PARA)
      REAL FUNCTION QUAGEV(F,PARA)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  QUANTILE FUNCTION OF THE GENERALIZED EXTREME-VALUE DISTRIBUTION
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(3)
      DATA ZERO/0D0/,ONE/1D0/
      U=PARA(1)
      A=PARA(2)
      G=PARA(3)
      IF(A.LE.ZERO)GOTO 1000
      IF(F.LE.ZERO.OR.F.GE.ONE)GOTO 10
      Y=-DLOG(-DLOG(F))
      IF(G.NE.ZERO)Y=(ONE-DEXP(-G*Y))/G
      QUAGEV=U+A*Y
      RETURN
C
   10 IF(F.EQ.ZERO.AND.G.LT.ZERO)GOTO 20
      IF(F.EQ.ONE .AND.G.GT.ZERO)GOTO 20
      WRITE(6,7000)
      QUAGEV=ZERO
      RETURN
   20 QUAGEV=U+A/G
      RETURN
C
 1000 WRITE(6,7010)
      QUAGEV=ZERO
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE QUAGEV :',
     *  ' ARGUMENT OF FUNCTION INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE QUAGEV : PARAMETERS INVALID')
      END


C===================================================== SAMLMR.FOR
      SUBROUTINE SAMLMR(X,N,XMOM,NMOM,A,B,*)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  SAMPLE L-MOMENTS OF A DATA ARRAY
C
C  PARAMETERS OF ROUTINE:
C  X      * INPUT* ARRAY OF LENGTH N. CONTAINS THE DATA, IN ASCENDING
C                  ORDER.
C  N      * INPUT* NUMBER OF DATA VALUES
C  XMOM   *OUTPUT* ARRAY OF LENGTH NMOM. ON EXIT, CONTAINS THE SAMPLE
C                  L-MOMENTS L-1, L-2, T-3, T-4, ... .
C  NMOM   * INPUT* NUMBER OF L-MOMENTS TO BE FOUND. AT MOST MAX(N,20).
C  A      * INPUT* ) PARAMETERS OF PLOTTING
C  B      * INPUT* ) POSITION (SEE BELOW)
C
C  FOR UNBIASED ESTIMATES (OF THE LAMBDA'S) SET A=B=ZERO. OTHERWISE,
C  PLOTTING-POSITION ESTIMATORS ARE USED, BASED ON THE PLOTTING POSITION
C  (J+A)/(N+B)  FOR THE J'TH SMALLEST OF N OBSERVATIONS. FOR EXAMPLE,
C  A=-0.35D0 AND B=0.0D0 YIELDS THE ESTIMATORS RECOMMENDED BY
C  HOSKING ET AL. (1985, TECHNOMETRICS) FOR THE GEV DISTRIBUTION.
C
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension X(N)
      DOUBLE PRECISION XMOM(NMOM),SUM(20)
      DATA ZERO/0D0/,ONE/1D0/

c     print*,'inside samlmr, fst=',(x(i),i=1,n)
      IF(NMOM.GT.20.OR.NMOM.GT.N)GOTO 1000
      DO 10 J=1,NMOM
   10 SUM(J)=ZERO
      IF(A.EQ.ZERO.AND.B.EQ.ZERO)GOTO 50
      IF(A.LE.-ONE.OR.A.GE.B)GOTO 1010
C
C         PLOTTING-POSITION ESTIMATES OF PWM'S
C
      DO 30 I=1,N
      PPOS=(I+A)/(N+B)
      TERM=X(I)
      SUM(1)=SUM(1)+TERM
      DO 20 J=2,NMOM
      TERM=TERM*PPOS
   20 SUM(J)=SUM(J)+TERM
   30 CONTINUE
      DO 40 J=1,NMOM
   40 SUM(J)=SUM(J)/N
      GOTO 100
C
C         UNBIASED ESTIMATES OF PWM'S
C
   50 DO 70 I=1,N
      Z=I
      TERM=X(I)
      SUM(1)=SUM(1)+TERM
      DO 60 J=2,NMOM
      Z=Z-ONE
      TERM=TERM*Z
   60 SUM(J)=SUM(J)+TERM
   70 CONTINUE
      Y=N
      Z=N
      SUM(1)=SUM(1)/Z
      DO 80 J=2,NMOM
      Y=Y-ONE
      Z=Z*Y
   80 SUM(J)=SUM(J)/Z
C
C         L-MOMENTS
C
  100 K=NMOM
      P0=ONE
      IF(NMOM-NMOM/2*2.EQ.1)P0=-ONE
      DO 120 KK=2,NMOM
      AK=K
      P0=-P0
      P=P0
      TEMP=P*SUM(1)
      DO 110 I=1,K-1
      AI=I
      P=-P*(AK+AI-ONE)*(AK-AI)/(AI*AI)
  110 TEMP=TEMP+P*SUM(I+1)
      SUM(K)=TEMP
  120 K=K-1
      XMOM(1)=SUM(1)
      IF(NMOM.EQ.1)RETURN
      XMOM(2)=SUM(2)
      IF(SUM(2).EQ.ZERO)GOTO 1020
      IF(NMOM.EQ.2)RETURN
      DO 130 K=3,NMOM
  130 XMOM(K)=SUM(K)/SUM(2)
      RETURN
C
 1000 WRITE(6,7000)
      RETURN 1
 1010 WRITE(6,7010)
      RETURN 1
 1020 WRITE(6,7020)
      RETURN 1
C
 7000 FORMAT(' *** ERROR *** ROUTINE SAMLMR : PARAMETER NMOM INVALID')
 7010 FORMAT(' *** ERROR *** ROUTINE SAMLMR :',
     *  ' PLOTTING-POSITION PARAMETERS INVALID')
 7020 FORMAT(' *** ERROR *** ROUTINE SAMLMR : ALL DATA VALUES EQUAL')
      END


C===================================================== SORT.FOR
      SUBROUTINE SORT(X,N)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  SORTS THE ARRAY X INTO ASCENDING ORDER
C
C  PARAMETERS OF ROUTINE:
C  X      *IN/OUT* ARRAY OF LENGTH N. CONTAINS THE NUMBERS TO BE SORTED.
C                  ON EXIT, CONTAINS THE SORTED NUMBERS.
C  N      * INPUT* NUMBER OF ELEMENTS TO BE SORTED
C
C  METHOD USED IS SHELL SORT WITH SEQUENCE OF INCREMENTS AS IN
C  D.F.KNUTH (1969) 'THE ART OF COMPUTER PROGRAMMING', VOL.3, P.95
C
      dimension X(N)
      IF(N.LE.1)RETURN
      J=4
      DO 10 I=1,100
      J=3*J+1
      IF(J.GE.N)GOTO 20
   10 CONTINUE
   20 CONTINUE
      M=(J/3)
      DO 60 MM=1,100
      M=M/3
      IF(M.EQ.0)RETURN
      DO 50 I=M+1,N
      TEST=X(I)
      J=I
      DO 30 JJ=1,100
      J=J-M
      IF(J.LE.0)GOTO 40
      IF(TEST.GE.X(J))GOTO 40
   30 X(J+M)=X(J)
   40 CONTINUE
   50 X(J+M)=TEST
   60 CONTINUE
      END



C************************************************
      subroutine grange(n,ld,d,dmin,dmax)
      logical ld
      dimension ld(n),d(n)
 
      dmin=1.e38
      dmax=-1.e38
      do i=1,n
        if(ld(i)) then
          dmin=min(dmin,d(i))
          dmax=max(dmax,d(i))
        endif
      enddo
 
      return
      end
 
C***************************************************
      subroutine getname (fnum,fname1,fname2)
      character*13  fname1
      character*15  fname2
      integer    fnum
 
      if (fnum.eq.1)  fname1='r_gribawips01'
      if (fnum.eq.2)  fname1='r_gribawips02'
      if (fnum.eq.3)  fname1='r_gribawips03'
      if (fnum.eq.4)  fname1='r_gribawips04'
      if (fnum.eq.5)  fname1='r_gribawips05'
      if (fnum.eq.6)  fname1='r_gribawips06'
      if (fnum.eq.7)  fname1='r_gribawips07'
      if (fnum.eq.8)  fname1='r_gribawips08'
      if (fnum.eq.9)  fname1='r_gribawips09'
      if (fnum.eq.10) fname1='r_gribawips10'
      if (fnum.eq.11) fname1='r_gribawips11'
      if (fnum.eq.12) fname1='r_gribawips12'
      if (fnum.eq.13) fname1='r_gribawips13'
      if (fnum.eq.14) fname1='r_gribawips14'
      if (fnum.eq.15) fname1='r_gribawips15'
      if (fnum.eq.16) fname1='r_gribawips16'
      if (fnum.eq.17) fname1='r_gribawips17'
      if (fnum.eq.18) fname1='r_gribawips18'
      if (fnum.eq.19) fname1='r_gribawips19'
      if (fnum.eq.20) fname1='r_gribawips20'
      if (fnum.eq.21) fname1='r_gribawips21'
      if (fnum.eq.22) fname1='r_gribawips22'
      if (fnum.eq.23) fname1='r_gribawips23'

      if (fnum.eq.1)  fname2='r_gribawips01.i'
      if (fnum.eq.2)  fname2='r_gribawips02.i'
      if (fnum.eq.3)  fname2='r_gribawips03.i'
      if (fnum.eq.4)  fname2='r_gribawips04.i'
      if (fnum.eq.5)  fname2='r_gribawips05.i'
      if (fnum.eq.6)  fname2='r_gribawips06.i'
      if (fnum.eq.7)  fname2='r_gribawips07.i'
      if (fnum.eq.8)  fname2='r_gribawips08.i'
      if (fnum.eq.9)  fname2='r_gribawips09.i'
      if (fnum.eq.10) fname2='r_gribawips10.i'
      if (fnum.eq.11) fname2='r_gribawips11.i'
      if (fnum.eq.12) fname2='r_gribawips12.i'
      if (fnum.eq.13) fname2='r_gribawips13.i'
      if (fnum.eq.14) fname2='r_gribawips14.i'
      if (fnum.eq.15) fname2='r_gribawips15.i'
      if (fnum.eq.16) fname2='r_gribawips16.i'
      if (fnum.eq.17) fname2='r_gribawips17.i'
      if (fnum.eq.18) fname2='r_gribawips18.i'
      if (fnum.eq.19) fname2='r_gribawips19.i'
      if (fnum.eq.20) fname2='r_gribawips20.i'
      if (fnum.eq.21) fname2='r_gribawips21.i'
      if (fnum.eq.22) fname2='r_gribawips22.i'
      if (fnum.eq.23) fname2='r_gribawips23.i'

      return
      end


      subroutine mean_median(n,fcsts,ave,xmed)

c...........................................................................
c  Purpose: computing ensemble mean and median for an n-member ensemble
c           at a fixed point
c
c  Input:
c	 n=the number of ensemble members
c	 fcsts=n forecasts in the ensemble
c  Output:
c	 ave=ensemble mean
c	 xmed=ensemble median
c Programer: Jun Du
c...........................................................................

      dimension fcsts(n),ens(n,3)

c 1. Initialization and ensemble mean
      ave=0.
      xmed=0.
      do 100 i=1,n
        ens(i,1)=i
        ens(i,2)=fcsts(i)
        ave=ave+fcsts(i)/float(n)
100   continue

c 2. Ordering ensemble data in an ascending manner
      call sortm(ens,n,3,2)

c 3. Finding median
      im=n/2
      if(mod(n,2).eq.1) then
       xmed=ens(im+1,2)
      else
       xmed=(ens(im,2)+ens(im+1,2))/2.
      endif

      return
      end

C**************************************
      subroutine sortm(a,n,nc,k)
      dimension a(n,nc),b(n,nc),js(n)
c
      do i1 = 1, n
      iless = 0
      imore = 0
      ieq   = 0
      aa=a(i1,k)
      do i2 = 1, n
       bb=a(i2,k)
       if ( aa.lt.bb ) iless = iless + 1
       if ( aa.gt.bb ) imore = imore + 1
       if ( aa.eq.bb ) then
          ieq   = ieq   + 1
          js(ieq) = i2
       endif
      enddo
       if ( ieq.eq.1) then
          b(imore+1,2)=aa
          b(imore+1,1)=i1
       else
        do i3 = 1, ieq
          b(imore+i3,2)=aa
          b(imore+i3,1)=js(i3)
        enddo
       endif
      enddo
      do jj= 1, n
        a(jj,3) = b(jj,1)
        a(jj,2) = b(jj,2)
      enddo
      return
      end


      subroutine probability(fst,inum,fvmin,fvmax,fvalue10
     &,fvalue25,fvalue50,fvalue75,fvalue90)

c  subroutine program    probability
c  Prgmmr: Yuejian Zhu           Org: np23          Date: 2004-09-30
c          Bo Cui                mod: wx20          Date: 2007-07-18
c
c  Change log:
c          04/03/2011, Jun Du: Modified to work for SREF
c
c This is subroutine to get 10%, 50% and 90% probability forecast               
c
c   subroutine
c              sort  ---> sorts the array x into ascending order 
c              samlmr---> sample L_moments of a dada array           
c              pelgev---> parameter estimation via L-moments for the genaralized extreme-value distribution
c              quagev---> quantile function of the generalized extreme-value diftribution                     
c
c   output 
c         favalue10    -- 10% probabilaity forecast
c         favalue90    -- 90% probabilaity forecast
c         favalue50    -- 50% probabilaity forecast
c         mode         -- mode forecast
c
c   Fortran 77 on IBMSP
c
C--------+---------+---------+---------+---------+----------+---------+--

c     implicit none
      dimension fst(inum)
      double precision fmon(3),opara(3),prob,amt
      real fvalue10,fvalue25,fvalue50,fvalue75,fvalue90,mode
c     double precision fvalue,fvalue10,fvalue50,fvalue90,mode
      integer          inum

ccc
ccc       Using L-moment ratios and GEV method
ccc
C--------+---------+---------+---------+---------+---------+---------+---------+
c     print *, 'Calculates the L-moment ratios, by prob. weighted'

c     print*,'inum=',inum
c     print*,'1/4=',int(inum/4)
c     print*,'3/4=',int(3*inum/4)
c     print*,'before sort, fst=',(fst(i),i=1,inum)
      call sort(fst,inum)
c     print*,'after sort, fst=',(fst(i),i=1,inum)
      fvmin=fst(1)
      fvmax=fst(inum)

      opara = 0.0D0

ccc for unbiased estimation, A=B=ZERO

c     print*,'before samlmr, fst=',(fst(i),i=1,inum)
      call samlmr(fst,inum,fmon,3,-0.0D0,0.0D0,*10)
c     print *, "fmon=",fmon
      call pelgev(fmon,opara,*10)
c     print *, "opara=",opara

      prob=0.1
      fvalue10=quagev(prob,opara)
      if(fvalue10.lt.fvmin) fvalue10=fvmin
c     print *, "10% value is ",fvalue10

      prob=0.25
      fvalue25=quagev(prob,opara)
c     print *, "25% value is ",fvalue25

      prob=0.5
      fvalue50=quagev(prob,opara)
c     print *, "50% value is ",fvalue50

c     mode=3*fvalue50-2*fmon(1)
c     print *, "mode1 is ",mode  

      prob=0.75
      fvalue75=quagev(prob,opara)
c     print *, "75% value is ",fvalue75

      prob=0.9
      fvalue90=quagev(prob,opara)
      if(fvalue90.gt.fvmax) fvalue90=fvmax
c     print *, "90% value is ",fvalue90

      return

10    print *, '  '
      print *, "Recalculate Probabilistic Forecast"
      write (*,'(10f8.2)') (fst(ii),ii=1,inum)

      fvalue10=0.5*(fst(1)+fst(2))
      fvalue90=0.5*(fst(inum)+fst(inum-1))
      if (fvalue90.eq.fst(inum)) then
       fvalue50=0.2*fvalue10+0.8*fvalue90
      else
       fvalue50=0.8*fvalue10+0.2*fvalue90
      endif
      fvalue25=fst(int(inum/4))
      fvalue75=fst(int(3*inum/4))

      print *, "10%,25%,50%,75% and 90% values are"
      write (*,'(5f8.2)') fvalue10,fvalue25,fvalue50,fvalue75
     &,fvalue90
      print *, '  '

      return
      end


CC On Zeus/WCOSS/EDDY no function DLGAMA(), following function is copied
Cfrom http://cpansearch.perl.org/src/AJOLMA/Statistics-Lmoments-0.03/lmoments/dlgama.f
CC
CC
      DOUBLE PRECISION FUNCTION DLGAMA(X)
C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C
C  LOGARITHM OF GAMMA FUNCTION
C
C  BASED ON ALGORITHM ACM291, COMMUN. ASSOC. COMPUT. MACH. (1966)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA SMALL,CRIT,BIG,TOOBIG/1D-7,13D0,1D9,2D36/
C
C         C0 IS 0.5*LOG(2*PI)
C         C1...C7 ARE THE COEFFTS OF THE ASYMPTOTIC EXPANSION OF DLGAMA
C
      DATA C0,C1,C2,C3,C4,C5,C6,C7/
     *   0.91893 85332 04672 742D 0,  0.83333 33333 33333 333D-1,
     *  -0.27777 77777 77777 778D-2,  0.79365 07936 50793 651D-3,
     *  -0.59523 80952 38095 238D-3,  0.84175 08417 50841 751D-3,
     *  -0.19175 26917 52691 753D-2,  0.64102 56410 25641 026D-2/
C
C         S1 IS -(EULER'S CONSTANT), S2 IS PI**2/12
C
      DATA S1/-0.57721 56649 01532 861D 0/
      DATA S2/ 0.82246 70334 24113 218D 0/
C
      DATA ZERO/0D0/,HALF/0.5D0/,ONE/1D0/,TWO/2D0/
      DLGAMA=ZERO
      IF(X.LE.ZERO)GOTO 1000
      IF(X.GT.TOOBIG)GOTO 1000
C
C         USE SMALL-X APPROXIMATION IF X IS NEAR 0, 1 OR 2
C
      IF(DABS(X-TWO).GT.SMALL)GOTO 10
      DLGAMA=DLOG(X-ONE)
      XX=X-TWO
      GOTO 20
   10 IF(DABS(X-ONE).GT.SMALL)GOTO 30
      XX=X-ONE
   20 DLGAMA=DLGAMA+XX*(S1+XX*S2)
      RETURN
   30 IF(X.GT.SMALL)GOTO 40
      DLGAMA=-DLOG(X)+S1*X
      RETURN
C
C         REDUCE TO DLGAMA(X+N) WHERE X+N.GE.CRIT
C
   40 SUM1=ZERO
      Y=X
      IF(Y.GE.CRIT)GOTO 60
      Z=ONE
   50 Z=Z*Y
      Y=Y+ONE
      IF(Y.LT.CRIT)GOTO 50
      SUM1=SUM1-DLOG(Z)
C
C         USE ASYMPTOTIC EXPANSION IF Y.GE.CRIT
C
   60 SUM1=SUM1+(Y-HALF)*DLOG(Y)-Y+C0
      SUM2=ZERO
      IF(Y.GE.BIG)GOTO 70
      Z=ONE/(Y*Y)
      SUM2=((((((C7*Z+C6)*Z+C5)*Z+C4)*Z+C3)*Z+C2)*Z+C1)/Y
   70 DLGAMA=SUM1+SUM2
      RETURN
C
 1000 WRITE(6,7000)X
      RETURN
C
 7000 FORMAT(' *** ERROR *** ROUTINE DLGAMA :',
     *  ' ARGUMENT OUT OF RANGE :',D24.16)
      END

        RECURSIVE SUBROUTINE quick_sort(list, order, mysize)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's
! Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which
! gives
! the positions of the elements in the original order.

        IMPLICIT NONE
        integer :: mysize
        REAL, DIMENSION (mysize), INTENT(IN OUT)  :: list
        INTEGER, DIMENSION (mysize), INTENT(OUT)  :: order

! Local variable
        INTEGER :: i

!        write(0,*) 'in here'
!        write(0,*) 'size(list): ', size(list)
!        write(0,*) 'size(order): ', size(order)
        DO i = 1, SIZE(list)
          order(i) = i
        END DO
!        write(0,*) 'order defined'

        CALL quick_sort_1(1, SIZE(list))

        CONTAINS

        RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
        INTEGER             :: i, j, itemp
        REAL                :: reference, temp
        INTEGER, PARAMETER  :: max_simple_sort_size = 6

        IF (right_end < left_end + max_simple_sort_size) THEN
          ! Use interchange sort for small lists
          CALL interchange_sort(left_end, right_end)

        ELSE
          ! Use partition ("quick") sort
          reference = list((left_end + right_end)/2)
          i = left_end - 1; j = right_end + 1
        
          DO
            ! Scan list from left end until element >= reference is found
            DO
              i = i + 1
              IF (list(i) >= reference) EXIT
            END DO
            ! Scan list from right end until element <= reference is found
            DO
              j = j - 1
              IF (list(j) <= reference) EXIT
            END DO


            IF (i < j) THEN
              ! Swap two out-of-order elements
              temp = list(i); list(i) = list(j); list(j) = temp
              itemp = order(i); order(i) = order(j); order(j) = itemp
            ELSE IF (i == j) THEN
              i = i + 1
              EXIT
            ELSE
              EXIT
            END IF
          END DO
        
          IF (left_end < j) CALL quick_sort_1(left_end, j)
          IF (i < right_end) CALL quick_sort_1(i, right_end)
        END IF
        
        END SUBROUTINE quick_sort_1


        SUBROUTINE interchange_sort(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
        INTEGER             :: i, j, itemp
        REAL                :: temp

        DO i = left_end, right_end - 1
          DO j = i+1, right_end
            IF (list(i) > list(j)) THEN
              temp = list(i); list(i) = list(j); list(j) = temp
              itemp = order(i); order(i) = order(j); order(j) = itemp
            END IF
          END DO
        END DO

        END SUBROUTINE interchange_sort
        END SUBROUTINE quick_sort
