C$$$me  MAIN PROGRAM DOCUMENTATION BLOCK $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C MAIN PROGRAM: SREF_ENS_GEN
C   PRGMMR: DU               ORG: NP21        DATE: 2002-01-02
C
C ABSTRACT: READS IN SREF GRIB AND OUTPUTS ENSEMBLE GRIB 
C           PRODUCTS. 
C
C PROGRAM HISTORY LOG:
c
c Old Version:
c 2001-07-10    Jun Du, initial program
c
c Generic Version: 
c
c 2005-08-01 Author: Binbin Zhou 
c
c       Re-organize the code structure and re-write the program in following features:
c
c        (1) The generic version dosn't fix GRIB# like 212 as old version. It can process
c            any GRIB#, and the domain (IM, JM) for the defined GRIB# is also flexible. 
c            As long as GRIB# is specified in 'filename', this program can automatically 
c            retrievs the domain (IM, JM, JF=IM*JM), then all arrays associated with this 
c            domain are dynamically allocated (on the fly)
c        (2) Ensemble members (IENS) are not fixed, but depends on the input file "filename'
c            in which, the all of soft linked file name are listed in order. The order and
c            number of files are changeable, the program can automatically allocate them
c        (3) The variables are changeable by using the 'variable.tbl' file, in which, two kinds of
c            variables are defined: (i) direct variables: read from GRIB files. Direct  
c            variables will also be computed for mean/spread an probability depending 
c            on the settings in each variable redcord (ii) derived variables: they are 
c            not read from GRIB files but need to do some derivation from direct variables.
c        (4) Derived variables need some algorithm to derive, Currently, only following derived 
c            variables are defined in this program:
c            (A) Thickness of two (pressure) levels: mean/spread/probability         
c            (B) 3hr, 6hr, 12hr, 24hr Accumulated precipitation: mean/spread/probability 
c            (C) 3hr, 6hr, 12hr, 24hr Accumulated snowfall: mean/spread/probability
c            (D) Precipitation type: probability
c            (E) Dominant precipitation: Mean
c            If want to process other derived product, user must provide their algorithms 
c            and then the computation code should be put into this program       
c        (5) User can request any kind of direct variables, either mean/spread or probability
c            by editing the 'variable.tbl' file (read README file for its gramma)
c 
c change log:
c 2006-01-12:Binbin Zhou
c             Add Geoff Manikin's algorithm to determine dominant precip type 
c 2006-03-01:Binbin Zhou
c             Modify wind spread computation by considering the wind directions, suggested 
c             by Jun Du            
c 2006-03-30: Jun Du: convert the unit from K to C for lifted index 
c
c 2006-04-30: Binbin Z. Add DTRA variance variables 
c
c 2006-05-15: Binbin Zhou
c             Modify to accept Beijing 2008 Olympic Game Grid#137
c
c 2006-07-21: Binbin Z. Add fog probability
c
c 2007-01-09: Binbin Z. Add freezing rain for GFS ensemble
c
c 2007-01-15: Binbin Z. Add Absolute vorticity for global ensemble since GFS has no absolute vorticity output
c 
c 2007-04-05: Binbin Z. Add High resolution WRF testing grid#137 
c
c 2007-10-09: Binbin Z. Add flight-restriction probability computation
c 
c 2009-04-07: Binbin Z. Adopted for VSREF (by adding member weights and RUC/VSREF's grid#130 ) 
c
c 2009-09-09: Binbin Z. Add convection adapted from Steve W. of GSD
C
c 2010-06-20: Binbin Z. Add precip accumulation for VRSEF version
c
c 2010-8-18:  Binbin Z. Add new fog algorithm (based on Zhou and Ferrier 2008)
c 2011-3-18:  Binbin Z. adapt for NARRE-Time Lag  (based on Zhou and Ferrier 2008)
c
c 2011-04-15: Binbin Z. Add cpat lightning products (based on David Bright's SPC program)
c 2011-04-18: Binbin Z. Add fire weather probability computation
c 2011-04-18: Binbin Z. Add LLWS code
c 2011-06-10: Binbin Z. Add cptp_dryt based on David B.'s code (but dryt_hrly_rgn3 seems not good)
c 2011-06-20: Binbin Z. Add severe thunder storm potential based on David B's code
c 2012-01-20: Binbin Z. Add missing array to deal with after copygb of hiresw WRF-east/west over CONUS 
c                       Modified to work on NSSE 
c 2012-02-21: Binbin Z. Modify accumulated precip computation method: read them directly from prcip.* files
c                       so they become direct vriables in the table 
c 
c 2013-03-28: Binbin Z. Modify to work on grib2 and structure change to reduce memory by read one variable
c                       then process it and pack its ensmeble products into output files
c
c 2013-12-21: Binbin Z. Modify  missing array from missing(jf,iens) to missing(maxvar,iens) after 
c                       hiresWRFs are all in same CONUS grid, and deal with field missing in GRIB2 file 
c                       and getGB2 issue (not stop if the field is not in the GRIB2 file
c
c 2014-03-15: Binbin Z. Add Storm-related fileds requested by Jacob Caley of EMC
c
c 2014-04-21: Binbin Z. Add Flash Flood and Intense Rain (FFaIR) summer experiment fields
c
c 2015-06-10: Matthew P. Considers the bitmap from all ingested grids,
c                        and outputs only for points where all source
c                        grids have valid data (bmap_f) .
c
c 2015-11-12: Binbin Z. Add SPC's SSEO products: Spaghetti, Max,  etc
c                       specified by 'X','S' in the product table;
c                       Also Add 10m wind and 2m RH joint prob to derived products for fire weather
c 2015-12-02: Binbin Z. Add Gaussian smoothing for probability
c
c 2015-12-09: Binbin Z. Add neighborhood max value for probability computation
c 
c 2016-02-02: Matthew P. Add Probability-Matched mean QPF
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C$$$
C

       use grib_mod
       include 'parm.inc'

      type(gribfield) :: gfld, gfld_temp, gfld_neighb_restore, gfld_ffg, 
     &                   gfld_ri, gfld_lightning
C  raw data
      real,allocatable,dimension(:,:,:)   :: rawdata_mn, rawdata_pr   !jf,iens, maxmlvl
      real,allocatable,dimension(:,:,:) :: rawdata_mn_2d
      real,allocatable,dimension(:,:,:)   :: precip                   !jf,iens, number of fcst output files
      real,allocatable,dimension(:,:,:)   :: mrk_ice                  !jf,iens, number of fcst output files
      real,allocatable,dimension(:,:,:)   :: mrk_frz                  !jf,iens, number of fcst output files
      real,allocatable,dimension(:,:,:)   :: mrk_snow                 !jf,iens, number of fcst output files

C mean
      real,allocatable,dimension(:,:) :: vrbl_mn                    !jf, maxmlvl
      real,allocatable,dimension(:,:) :: vrbl_mn_pm                 !jf, maxmlvl
      real,allocatable,dimension(:,:) :: vrbl_mn_locpm               !jf, maxmlvl
      real,allocatable,dimension(:,:) :: vrbl_mn_blend              !jf, maxmlvl
      real,allocatable,dimension(:,:) :: vrbl_mn_blendlpm           !jf, maxmlvl
      real,allocatable,dimension(:,:) :: derv_mn                    !jf, maxmlvl
      real,allocatable,dimension(:,:) :: vrbl_mn_2d, vrbl_lpm_2d
      real,allocatable,dimension(:,:) :: vrbl_pmmn_2d

C spread
      real,allocatable,dimension(:,:) :: vrbl_sp                    !jf, maxmlvl
      real,allocatable,dimension(:,:) :: derv_sp                    !jf, maxmlvl

      real,allocatable,dimension(:) :: ppt3_sp,ppt6_sp,ppt12_sp,
     *                                 ppt24_sp,s12_sp

C probability
      real,allocatable,dimension(:,:,:) :: vrbl_pr                   !jf, maxplvl, maxtlvl
      real,allocatable,dimension(:,:,:) :: derv_pr                   !jf, maxplvl, maxtlvl

C Mix,min, 10,25,50 and 90% mean products
      real,allocatable,dimension(:,:,:) :: mxp8                     !jf,maxmlvl,7

C Return Interval fixed data
      real, allocatable, dimension(:) :: return_int


C others 
       character*19, allocatable, dimension(:) :: fhead
       character*50 :: fname
       real,allocatable,dimension(:) :: p03mp01,vrbl_mn_use             !jf
       real,allocatable,dimension(:,:)  :: ptype_mn,ptype_pr,ptype_pr2 !jf,4 

       real,allocatable,dimension(:,:,:) :: derv_dtra                  !jf,maxmlvl,8 for DTRA requests
       real,allocatable, dimension(:)    :: Hsfc                       !surface height for DTRA  
 
       integer,allocatable,dimension(:,:)   :: missing                 ! to deal with missing data 
       integer,allocatable,dimension(:)     :: miss
       real,allocatable,dimension(:) ::           apoint               !iens

       integer :: patch_nx, patch_ny, ovx, ovy

C for get grib size jf=im*jm
  
       integer iyr,imon,idy,ihr
       character*50 gdss(400)
       integer IENS, gribid, kgdss(200), lengds,im,jm,km,jf

       integer  leadtime, interval, loutput

cc%%%%%%%  8. To be added if necessary ...........................
C original                                                                                                                                           
       dimension kgds(25)
       character*20 mnout,pmin,pmax,pmod,pp10,pp25,pp50,pp75,pp90
       character*20 spout,pmmnout,avgout,ffriout,locpmmnout
       character*20 prout,locavgout
       character*40 files(50),prcps(50)
     


C for variable table:
        Integer numvar, nderiv
        Character*4 vname(maxvar)
        Integer k5(maxvar), k6(maxvar), k4(maxvar)     !k4-category#, k5-paremter#, k6-surface id
        Character*1 Msignal(maxvar), Psignal(maxvar)
        Integer Mlvl(maxvar), MeanLevel(maxvar,maxmlvl),Lm
        Integer Plvl(maxvar), ProbLevel(maxvar,maxplvl),Lp
        Integer Tlvl(maxvar),Lth
        Character*1 op(maxvar)
        Real    Thrs(maxvar,maxtlvl)
        Character*5 eps       
  
c    for derived variables
        Character*4 dvname(maxvar)
        Integer dk5(maxvar), dk6(maxvar), dk4(maxvar)
        Character*1 dMsignal(maxvar), dPsignal(maxvar)
        Integer dMlvl(maxvar), dMeanLevel(maxvar,maxmlvl)
        Integer dPlvl(maxvar), dProbLevel(maxvar,maxplvl)
        Character*1 dop(maxvar)
        Integer dTlvl(maxvar)
        Real    dThrs(maxvar,maxtlvl)
        Integer MPairLevel(maxvar,maxmlvl,2)
        Integer PPairLevel(maxvar,maxplvl,2)

c   for max,min,10,25,50,90% mean products
        Character*4 qvname(maxvar)
        Integer qk5(maxvar), qk6(maxvar), qk4(maxvar)
        Character*1 qMsignal(maxvar)
        Integer qMlvl(maxvar), qMeanLevel(maxvar,maxmlvl)
        
        real  weight(30), gauss_sig                         
        character*20 filenames
        character*3 cfhr                                       
        integer est                                 !east time for convection code
        integer ifunit(50),ipunit(50)               !member, prcp member files units         
        character*2 cycle(24)
        character*7 mbrname(50)

        logical*1, allocatable:: bmap_f(:)
        integer jptyp2(4)                           !for precip type jpd2
        integer iqout(8) 
        integer III

        character Tsignal(maxvar)                   !For Gaussian smoothing
        character dTsignal(maxvar)                   !For Gaussian smoothing

        common /tbl/numvar,
     +              vname,k4,k5,k6,Mlvl,Plvl,Tlvl,
     +              MeanLevel,ProbLevel,Thrs,
     +              Msignal,Psignal,op
                                                                                                                                                                
        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop
                           
     
        common /qtbl/nmxp,
     +              qvname,qk4,qk5,qk6,qMlvl,
     +              qMeanLevel,qMsignal


       common /Tsgn/Tsignal,dTsignal 
                                                                                                                                
        data (iqout(i),i=1,8)
     +  /301,302,303,304,305,306,307,308/

        data (cycle(i),i=1,24)
     +  /'00','01','02','03','04','05','06','07','08',
     +   '09','10','11','12','13','14','15','16','17',
     +   '18','19','20','21','22','23'/

        data (jptyp2(i),i=1,4)
     +  /192,193,194,195/                     !rain, frzr, icep, snow

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c (I)   Read table file to get wanted ensemble variable information
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        numvar=0
        nderiv=0
        nmxp=0

        nunit=101
        open (nunit, file='variable.tbl', status='old')
        
        write(*,*) 'reading variable.tbl ...'                           
        call readtbl(nunit)
        write(*,*) 'read variable.tbl done'

        close(nunit)

        write(*,*) 'Check direct variable reading:'
        do i = 1, numvar
          write(*,*) vname(i),k4(i),k5(i), k6(i),Msignal(i), Mlvl(i),
     +                 (MeanLevel(i,j),j=1,Mlvl(i)),
     +    Psignal(i), Plvl(i), (ProbLevel(i,j),j=1,Plvl(i)),
     +    op(i), Tlvl(i),   (Thrs(i,j),j=1,Tlvl(i))
151      format(a4,1x,3i4,a2,i2,<Mlvl(i)>(1x,i4),a2,i2,<Plvl(i)>(1x,i4),
     +   a2, i2, <Tlvl(i)>(1x,F7.1))   
        end do

        write(*,*) 'Check derived variable reading:'
        do i = 1, nderiv
         if(dk6(i).eq.101) then
           write(*,*) dvname(i),dk4(i),dk5(i), dk6(i),dMsignal(i),
     +               dMlvl(i),(dMeanLevel(i,j),j=1,dMlvl(i)),
     +              (MPairLevel(i,j,1),j=1,dMlvl(i)),
     +              (MPairLevel(i,j,2),j=1,dMlvl(i)),
     +          dPsignal(i),dPlvl(i),(dProbLevel(i,j),j=1,dPlvl(i)),
     +              (PPairLevel(i,j,1),j=1,dPlvl(i)),
     +              (PPairLevel(i,j,2),j=1,dPlvl(i)),
     +              dop(i), dTlvl(i),   (dThrs(i,j),j=1,dTlvl(i))
         else
           write(*,*) dvname(i),dk4(i),dk5(i), dk6(i),dMsignal(i),
     +              dMlvl(i),(dMeanLevel(i,j),j=1,dMlvl(i)),
     +          dPsignal(i),dPlvl(i),(dProbLevel(i,j),j=1,dPlvl(i)),
     +             dop(i), dTlvl(i),   (dThrs(i,j),j=1,dTlvl(i))
         end if
        end do

        write(*,*) 'Check max,min,10,25,50,75,90% mean product reading:'
        do i = 1, nmxp
         write(*,*)qvname(i),qk4(i),qk5(i),qk6(i),qMsignal(i),
     +        qMlvl(i),(qMeanLevel(i,j),j=1,qMlvl(i))
        end do


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c (II): Read filename file to get # of ensemble available, grib#, vertical pressure level,
c       last forecast time, forecast outout interval and grib file heads, 
c       so that arrays can be dinamically allocatable
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       filenames="filename" 
       write(*,*) 'read ',filenames
       open(1, file=filenames, status='old')
       read(1,*) iyr,imon,idy,ihr,ifhr,gribid,KM,leadtime,interval, !ihr-cycle hour , ifhr-forecast hour
     +             loutput
       write(*,*)iyr,imon,idy,ihr,ifhr,gribid,KM,leadtime,interval,
     +             loutput

       !read(1,*) (cfhr(k),k=1,loutput) 

       IENS=0
       do irun=1,50
        if(ifhr.lt.100) then
         read(1,301,END=302) weight(irun),files(irun),mbrname(irun)
        else
         read(1,303,END=302) weight(irun),files(irun),mbrname(irun)
        end if

	write(0,*) 'files(irun)(1:4): ', files(irun)(1:4)

         if( files(irun)(1:4).ne.'href' ) then       
          mbrname(irun)=files(irun)
         end if
         write(*,301) weight(irun), files(irun),mbrname(irun)
         IENS=IENS+1
         prcps(irun)='prcip'//files(irun)(5:17)              !get prcip file names
         write(*,*) weight(irun), files(irun),prcps(irun)
       end do

       close (1)

  301  format(f7.2, 1x, a17,4x,a7)
  303  format(f7.2, 1x, a18,4x,a7)

  302  continue
      write(*,*) 'IENS=',IENS, (weight(irun),irun=1,iens)

      DO 3000 itime=ifhr,ifhr          !Here ifhr is forecast hour, e.g. .f12, '12' is ifhr

         if(itime.lt.100) then
          write(cfhr,'(i2.2)') itime
         else
          write(cfhr,'(i3.3)') itime
         end if
        print*,'Initializing variables'

       est = ihr + itime - 6                                          !Eastern time For convection computation
       if (est.lt.0) est= est + 24 


CCCC Binbin Zhou Note:

       if(gribid.eq.255) then   !For HRRR grid
         im=1799
         jm=1059
         jf=im*jm
         dx=3000.0
         dy=3000.0
       elseif (gribid .eq. 1099) then ! RRFS CONUS grid
         im=1746
         jm=1014
         jf=im*jm
         dx=3000.0
         dy=3000.0
       elseif (gribid .eq. 999) then ! AK grid
         im=825
         jm=603
         jf=im*jm
         dx=5000.0
         dy=5000.0
       elseif (GRIBID.eq.998) then ! HI 5 km grid
         im=223
         jm=170
         jf=im*jm
         dx=4500.0
         dy=4500.0
       elseif (GRIBID.eq.997) then ! PR 5 km grid
         im=340
         jm=208
         jf=im*jm
         dx=4500.0
         dy=4500.0
       else
         call makgds(gribid, kgdss, gdss, lengds, ier)
         im=kgdss(2)
         jm=kgdss(3)
         jf=kgdss(2)*kgdss(3)
         dx=1.*kgdss(8)
         dy=1.*kgdss(9)
       end if
 
       write(*,*) 'im, jm, jf,dx,dy =', im, jm, jf,dx,dy
       write(*,*) 'leadtime, interval, loutput=', 
     +             leadtime, interval, loutput

       if (.NOT.allocated(fhead)) then
        allocate(fhead(iens))   
       end if
       if (.NOT.allocated(apoint)) then
        allocate(apoint(iens)) 
       end if
       if (.NOT.allocated(missing)) then
        allocate(missing(maxvar,iens))
       end if
       if (.NOT.allocated(miss)) then
        allocate(miss(iens))
       end if

       if (.NOT.allocated(p03mp01)) then
        allocate(p03mp01(jf))
       end if


       eps=files(1)(1:4)
       missing=0 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c (III): Open all ensemble data files 
c In these lines put here, we need the date and time of the forecast,
c the name of the file.  We need to open the file.  Then we can point to
c the correct location of the forecast in the forecast file.  If this
c is an ensemble, then we also need another loop for the forecast name
c (cntl, n1, p1, etc.).  That loop should be inside the date/time loop.



       do 500 irun=1,iens

        ifunit(irun)=10+irun
        call baopenr(ifunit(irun),files(irun),ier1)
        write(*,*)'open#',irun,ifunit(irun),files(irun),'err=',ier1 

        ipunit(irun)=50+irun
        call baopenr(ipunit(irun),prcps(irun),ier2)
        write(*,*)'open#',irun,ipunit(irun),prcps(irun),'err=',ier2

500     continue 

        write(*,*) 'check numvar=',numvar

        imean=201
        isprd=202
        iprob=203
        ipmmn=204
        iavg=205
        iffri=206
        ilocpmmn=207
!        ilocavg=208

      mnout=trim(eps)//'.mean.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
      pmmnout=trim(eps)//'.pmmn.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
      locpmmnout= 
     &     trim(eps)//'.lpmm.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
      spout=trim(eps)//'.sprd.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
      prout=trim(eps)//'.prob.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
      avgout=trim(eps)//'.avrg.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
!      locavgout= 
!     &     trim(eps)//'.lavg.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
      ffriout=trim(eps)//'.ffri.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)

      call baopen (imean, mnout, iret)
      if (iret.ne.0) write(*,*) 'open ', mnout, 'err=', iret
      call baopen (isprd, spout, iret)
      if (iret.ne.0) write(*,*) 'open ', spout, 'err=', iret
      call baopen (iprob, prout, iret)
      if (iret.ne.0) write(*,*) 'open ', prout, 'err=', iret
      call baopen (ipmmn, pmmnout, iret)
      if (iret.ne.0) write(*,*) 'open ', pmmnout, 'err=', iret
      call baopen (ilocpmmn, locpmmnout, iret)
      if (iret.ne.0) write(*,*) 'open ', locpmmnout, 'err=', iret
      call baopen (iavg, avgout, iret)
      if (iret.ne.0) write(*,*) 'open ', avgout, 'err=', iret
!      call baopen (ilocavg, locavgout, iret)
!      if (iret.ne.0) write(*,*) 'open ', locavgout, 'err=', iret

! limit to only for CONUS (based on what?)
      call baopen (iffri, ffriout, iret)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c (IV) Main tasks:
c  Loop1: process direct variables
c  Loop2: process derived variables
c  Loop3: process mxp variables (min,max,mod, 10, 250 50 75 90%) 

c  Loop 1: for direct variables
c   1-0: Allocate nesessary arrays
c   1-1: In this loop, one field (nv) data are read from all members (1001 loop),
c   1-2: Compute its mean, spread, and prob,
c   1-3: Pack them into output mean, spread, prob files   
c   then, loop for next field
c   1-4: deallocate allocated arrays 
c
c      write(*,*) 'numvar before DO 2001 ', numvar
       DO 2001 nv = 1, numvar
       loop222: do III=1,1

       write(*,*) 'nv Mlvl(nv),Plvl(nv),Tlvl(nv)=',
     +    nv, Mlvl(nv),Plvl(nv),Tlvl(nv)

c  Loop  1-0: Allocate nesessary arrays
         if (.NOT.allocated(rawdata_mn)) then
           allocate (rawdata_mn(jf,iens,Mlvl(nv)))
         end if
         if (.NOT.allocated(rawdata_pr)) then
           allocate (rawdata_pr(jf,iens,plvl(nv)))
         end if
        if (.NOT.allocated(vrbl_mn)) then
           Lm=max(1,Mlvl(nv))                !Keep at least one vrbl_mn to pass it into packGB2 
           allocate (vrbl_mn(jf,Lm))         !in case Mlvl(nv)=0
           allocate (vrbl_mn_pm(jf,Lm))       
           allocate (vrbl_mn_locpm(jf,Lm))       
           allocate (vrbl_mn_blend(jf,Lm))       
           allocate (vrbl_mn_blendlpm(jf,Lm))       
         end if
        if (.NOT.allocated(vrbl_sp)) then
           Lm=max(1,Mlvl(nv))
           allocate (vrbl_sp(jf,Lm))
         end if
        if (.NOT.allocated(vrbl_pr)) then
           Lp=max(1,Plvl(nv))
           Lth=max(1,Tlvl(nv))
           allocate (vrbl_pr(jf,Lp,Lth))
        end if
        if (.not. allocated(return_int)) then
           allocate(return_int(jf))
        end if

         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Loop 1-1: Read direct variable's GRIB2 data from all members 


        loop1001: DO irun=1,iens         ! members

         jret=0
         kret=0

         write(*,*) 'For ensmeble member#', irun

         !apcp1h,apcp3hr,apcp6hr,apcp12,apcp24 have been put into 
         !every members as direct variable


          if ((k4(nv).eq.2.and.k5(nv).eq.299).or.   !These field not in hrrrgsd file
     +        (k4(nv).eq.3.and.k5(nv).eq.1).or.
     +        (k4(nv).eq.1.and.k5(nv).eq.192).or.
     +        (k4(nv).eq.6.and.k5(nv).eq.1).or.
     +        (k4(nv).eq.7.and.k5(nv).eq.197) ) then

              if(mbrname(irun)(1:4).eq.'hrrr') then ! These 8 fields are not in GSD's HRRR
                 rawdata_pr(:,irun,:)=-9999.0
                 rawdata_mn(:,irun,:)=-9999.0
                 cycle loop1001
              end if

           end if  

          !First re-set  Product template# and ID for some products 
          IF ((k4(nv).eq.1.and.k5(nv).eq.8).or.          ! APCP/SNOW
     +       (k4(nv).eq.1.and.k5(nv).eq.15).or.
     +       (k4(nv).eq.1.and.k5(nv).eq.13).or.
     +       (k4(nv).eq.1.and.k5(nv).eq.11) ) then

!      here maybe for inserting?

              if (.NOT.allocated(vrbl_mn_use)) then
                allocate (vrbl_mn_use(jf))
              endif
              if (.NOT.allocated(ptype_mn)) then
                allocate (ptype_mn(jf,4))
              end if
              if (.NOT.allocated(ptype_pr)) then
                allocate (ptype_pr(jf,4))
              end if
              if (.NOT.allocated(ptype_pr2)) then
                allocate (ptype_pr2(jf,4))
              end if
             
 
!	      if (slr_derv(1)<0.01) then
!              call preciptype (nv,ifunit,jpdtn,jf,iens,
!     +         ptype_mn,ptype_pr,ptype_pr2,slr_derv)
!              endif


!  end of inserted ptype stuff


         write(0,*) 'get APCP GRIB2 data for member ', irun

             if(mbrname(irun)(1:4).eq.'sref') then
              jpdtn=11
             else
              jpdtn=8
             end if

             jpd1=k4(nv)
	write(0,*) 'jpd1 defined(a): ', jpd1
             jpd2=k5(nv)
             jpd10=k6(nv)
             !jpd12 is determined by a specific level MeanLevel(nv,lv) to !ProbLevel(nv,lv) later on

             if (vname(nv).eq.'AP1h'.or.vname(nv).eq.'SN1h'.or. 
     &           vname(nv).eq.'FFG1') then        !AP1h,Ap3h,Ap6h, AP12,Ap24 should be hardcopy in the variable tbl
	     write(0,*) 'jpd27=1 for AP1h or SN1h'
             jpd27=1
             else if (vname(nv).eq.'AP3h'.or.vname(nv).eq.'SN3h' .or. 
     &                vname(nv).eq.'A3RI' .or. vname(nv).eq.'FFG3') then
                if(ifhr.lt.3 .or. mod(ifhr,3) .ne. 0) exit loop222
                jpd27=3
             else if (vname(nv).eq.'AP6h'.or.vname(nv).eq.'SN6h' .or.
     &                vname(nv).eq.'A6RI'.or.vname(nv).eq.'FFG6') then
                if(ifhr.lt.6) exit loop222
                jpd27=6
             else if (vname(nv).eq.'AP12'.or.vname(nv).eq.'SN12' .or. 
     &                vname(nv).eq.'A12R'.or.vname(nv).eq.'FF12') then
                if(ifhr.lt.12) exit loop222
                jpd27=12
             else if (vname(nv).eq.'AP24'.or.vname(nv).eq.'SN24' .or.
     &                vname(nv).eq.'A24R'.or.vname(nv).eq.'FF24') then
                if(ifhr.lt.24) exit loop222
                jpd27=24
             else
                write(*,*) 'Using wrong APCP name!'
             end if

           igrb2=ipunit(irun)
	
	write(0,*) 'vname(nv), igrb2: ',vname(nv), igrb2

          ELSE                   !Non-APCP/Snow

            jpd1=k4(nv)
	write(0,*) 'jpd1 defined(b): ', jpd1
            jpd2=k5(nv)
            !jpd12 is determined by a specific level MeanLevel(nv,lv) to
            !!ProbLevel(nv,lv) later on
            jpd10=k6(nv)
            jpd27=-9999

            !Non-APCP/SNOW products, SSEO products 
            if (trim(eps).eq.'sseo') then 
 
               write(*,*) 'Check sseo here ...' 

              if (vname(nv).eq.'UH3h'.or.
     +            vname(nv).eq.'RF3h'.or.
     +            vname(nv).eq.'WS3h'.or.
     +            vname(nv).eq.'US3h') then         !3hr max
                  if(ifhr.lt.3) exit loop222
                         jpd27=1
                         igrb2=ipunit(irun)
              else if (vname(nv).eq.'UH24'.or.
     +                 vname(nv).eq.'RF24'.or.
     +                 vname(nv).eq.'WS24'.or.
     +                 vname(nv).eq.'US24' ) then    !24hr max
                        if(ifhr.lt.24) exit loop222
                         jpd27=2
                         igrb2=ipunit(irun)
               else

                  igrb2=ifunit(irun)

               end if

               jpdtn=0

            else   !Non-APCP/SNOW products, Non-SSEO product      

              if(mbrname(irun)(1:4).eq.'sref'.or.
     +           mbrname(irun)(1:4).eq.'gefs') then
               jpdtn=1
              else
               jpdtn=0
              end if

              if((k4(nv).eq.2.and.k5(nv).eq.222).or.
     +         (k4(nv).eq.2.and.k5(nv).eq.223).or.
     +         (k4(nv).eq.16.and.k5(nv).eq.198).or.
     +         (k4(nv).eq.7.and.k5(nv).eq.199).or.
     +         (k4(nv).eq.2.and.k5(nv).eq.220).or.
     +         (k4(nv).eq.2.and.k5(nv).eq.221) ) then
                if(mbrname(irun)(1:4).eq.'sref') then
                 jpdtn=11
                else 
                 jpdtn=8
                end if 
              end if

              igrb2=ifunit(irun)

            end if
         
          END IF !end of APCP/SNOW AND /Non-APCP/SNOW


        loop501: do lv=1,Mlvl(nv)

          jpd12=MeanLevel(nv,lv)


!!! CAPE stuff

          if(jpd1.eq.7.and.jpd2.eq.6.and.jpd10.eq.108 .and. 
     +               (jpd12 .eq. 9000 .or. jpd12 .eq. 18000)) then
	write(0,*) 'attempting to deal with CAPE'
            jpd27=0    ! CAPE
          end if

          if(jpd1.eq.7.and.jpd2.eq.7.and.jpd10.eq.108 .and. 
     +               (jpd12 .eq. 9000 .or. jpd12 .eq. 18000)) then
	write(0,*) 'attempting to deal with CIN'
            jpd27=0    ! CIN
          end if

!!! end CAPE stuff

          if(jpd1.eq.0.and.jpd2.eq.192.and.jpd10.eq.106) then
            jpd27=10    !only 0-10 cm layer soil moisture is requested
          end if

          if(jpd1.eq.0.and.jpd2.eq.2.and.jpd10.eq.106) then
            jpd27=10    !only 0-10 cm layer soil temperature is requested
          end if

          write(0,*)'a1 - readGB2:', 
     &          igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27
          call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld, eps, jret)
          write(0,*) 'a2 - readGB ',igrb2,' for mean kret=',kret 


         if (jret .eq. 0.and.trim(eps).eq.'href') then

	write(0,*) 'in here with jf: ', jf


           if ( .not. allocated(bmap_f)) then
             allocate(bmap_f(jf))
	if (jf .ne. 37910 .and. jf .ne. 70720 .and. jf .ne. 1770444) then
             bmap_f=gfld%bmap
        else
             bmap_f=.true.
        endif
           endif

! avoid accounting for echo top bitmap (and cloud base/ceiling from HRRR) (and REFC from FV3 now)
! and FV3 soil
! and FV3 WEASD

	if (jf .ne. 37910 .and. jf .ne. 70720 .and. jf .ne. 1770444) then

           if ( jpd2.ne.197 .and. jpd2.ne.5 .and. 
     &          jpd2.ne. 192 .and. jpd2 .ne. 2 .and. 
     &          jpd2 .ne. 13 ) then

            do J=1,jf
             if ( (bmap_f(J)) .and. (.not. gfld%bmap(J))) then
              bmap_f(J)=.false.
             endif
            enddo
           endif

        endif

         endif


         if (jret.ne.0) cycle loop501

         rawdata_mn(:,irun,lv)=gfld%fld

       end do loop501    !end of read rawdata for mean 

       loop502: do lv=1,Plvl(nv)

          jpd12=probLevel(nv,lv)

        call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld,eps, kret)

	if (kret .eq. 0 .and. k5(nv) .eq. 192 .and. jpd1 .eq. 17) then
            gfld_lightning=gfld
        endif

         if (kret .ne. 0 .and. k5(nv) .eq. 220 ) then
! check for 100 hPa
	 write(0,*) 'look for 100-1000 hPa UVV'
         jpd12=100
        call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld,eps, kret)

         endif


        if (kret .eq. 0. and. trim(eps).eq.'href') then

         if ( .not. allocated(bmap_f)) then
                allocate(bmap_f(jf))
	if (jf .ne. 37910 .and. jf .ne. 70720 .and. jf .ne. 1770444) then
                bmap_f=gfld%bmap
        else
                bmap_f=.true.
        endif
         endif

	if (jf .ne. 37910 .and. jf .ne. 70720 .and. jf .ne. 1770444) then
! avoid accounting for echo top bitmap
! and FV3 soil
! and FV3 WEASD
           if ( jpd2.ne.197 .and. jpd2.ne.5 .and. 
     &          jpd2.ne. 192 .and. jpd2 .ne. 2 .and.
     &          jpd2.ne. 13  ) then
            do J=1,jf
             if ( (bmap_f(J)) .and. (.not. gfld%bmap(J))) then
              bmap_f(J)=.false.
            endif
           enddo
          endif

        endif !37910
        endif

         if (kret.ne.0) cycle loop502
          rawdata_pr(:,irun,lv)=gfld%fld
         !do igrid =1,jf
         !  if(igrid.eq.1371523) then
         !  write(*,*) irun, igrid, rawdata_pr(igrid,irun,lv)
         !  end if
         !end do

        end do loop502   ! End of reading rawdata for prob  
        
         if(jret.ne.0.or.kret.ne.0) then
           missing(nv,irun)=1
           cycle loop1001
         else
           missing(nv,irun)=0 
         end if 


         !write(*,*) 'Correct/Process some variables ...'
         
         do lv = 1, Mlvl(nv)
           do igrid=1,jf

            if(k4(nv).eq.7.and.k5(nv).eq.192) then
              rawdata_mn(igrid,irun,lv) = rawdata_mn(igrid,irun,lv)
     +                                                       - 273.15  !Jun Du: change unit from K to C for lifted index
            end if

            if(k4(nv).eq.3.and.k5(nv).eq.5.and.(k6(nv).eq.2.or.
     +        k6(nv).eq.3))                                    then    !RUC cloud base/top with negative values

              if (rawdata_mn(igrid,irun,lv).le.0.0) then               !means no cloud
                 rawdata_mn(igrid,irun,lv) = 20000.0
              end if
            end if

           end do
          end do
 
          do lv = 1, Plvl(nv)

            do igrid=1,jf
             if(k4(nv).eq.7.and.k5(nv).eq.192) then
               rawdata_pr(igrid,irun,lv) = rawdata_pr(igrid,irun,lv)
     +                                                       - 273.15  !Jun Du: change unit from K to C for lifted index
             end if

             if(k4(nv).eq.3.and.k5(nv).eq.5.and.(k6(nv).eq.2.or.
     +          k6(nv).eq.3))                 then                     !RUC cloud base/top with negative values

                if (rawdata_pr(igrid,irun,lv).le.0.0) then               !means no cloud
                    rawdata_pr(igrid,irun,lv) = 20000.0
                end if
             end if
            end do

           !Get neighborhood max value, where A,B,C,D: for different !neighborhood radius
           if (trim(Psignal(nv)).eq.'A'.or.trim(Psignal(nv)).eq.'K' .or.
     +        trim(Psignal(nv)).eq.'B'.or. trim(Psignal(nv)).eq.'L' .or.
     +        trim(Psignal(nv)).eq.'M'.or. trim(Psignal(nv)).eq.'C' .or.
     +        trim(Psignal(nv)).eq.'D' ) then
             
            write(*,*) 'Call  neighborhood_max .......'

!! this modifies rawdata_pr
             call neighborhood_max(rawdata_pr(:,irun,lv),
     +                          jf,im,jm,Psignal(nv))     
           end if

           !Get neighborhood avrage value, where H,I for different
           !neighbohood radius
           if (trim(Psignal(nv)).eq.'H'.or.
     +         trim(Psignal(nv)).eq.'I') then
                    
            write(*,*) 'Call  neighborhood_avg ....'

            !write(*,*) 'before neighbor',rawdata_pr(797267,irun,lv)
               call neighborhood_average(rawdata_pr(:,irun,lv),
     +                          jf,im,jm,Psignal(nv))     
            !write(*,*) 'after neighbor',rawdata_pr(797267,irun,lv)

           end if

          end do

         end do loop1001   !end of iens loop for getting GRIB2 data for this direct variable
       
         write(*,*) 'Get direct variable data done for', nv

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Loop 1-2:   Compute mean/spread/prob for this direct variable

         vrbl_mn=0.
         vrbl_sp=0.
         vrbl_pr=0.

        IF(trim(Msignal(nv)).eq.'M'.or.trim(Msignal(nv)).eq.'P' .or. 
     &     trim(Msignal(nv)).eq.'L') THEN !P & L for Prob-matched mean
         
         do lv = 1, Mlvl(nv)
          do igrid=1,jf

           apoint = rawdata_mn(igrid,:,lv)

	if (igrid .eq. 937093) then
	write(0,*) 'lv, apoint(937093): ',lv, apoint
	endif

           miss = missing(nv,:)

          if ((k4(nv).eq.2.and.k5(nv).eq.299).or.   !These field not in hrrrgsd file
     +        (k4(nv).eq.3.and.k5(nv).eq.1).or.
     +        (k4(nv).eq.1.and.k5(nv).eq.192).or.
     +        (k4(nv).eq.6.and.k5(nv).eq.1).or.
     +        (k4(nv).eq.7.and.k5(nv).eq.197) ) then
             do ie=1,iens
               if(mbrname(ie)(1:4).eq.'hrrr') then
C	        write(0,*) 'set miss for hrrr: ', k4(nv),k5(nv)
                miss(ie)=1
               end if
             end do
          end if


           if(k4(nv).eq.19.and.k5(nv).eq.0) then  ! Visibility: conditional mean 
             do ie =1,iens
               if(apoint(ie).lt.1.0) miss(ie)=1   !exclude spc runs's error vis 
             end do
             call get_cond_mean (apoint,iens,24056.0,amean,aspread,
     +                              miss, weight) 

           else if(k4(nv).eq.3.and.k5(nv).eq.5.and.k6(nv).eq.215)then !  Ceiling: conditional mean
             call get_cond_mean (apoint,iens,20000.0,amean,aspread,
     +                               miss,weight)

           else if(k4(nv).eq.6.and.k5(nv).eq.6.and.k6(nv).eq.1) then  !  Sfc CLOUD (fog): liquid water content
                apoint = apoint * 1000.0                              !  kg/kg -> g/kg
                call get_cond_mean_lwc(apoint,iens,0.0,amean,aspread,
     +                               miss,weight)

           else if(k4(nv).eq.3.and.k5(nv).eq.5.and.k6(nv).eq.2) then  !  Cloud base: conditional mean
                do i =1,iens
                 if(apoint(i).le.0.) apoint(i)=20000.0
                end do
                call get_cond_mean (apoint,iens,20000.0,amean,aspread,
     +                               miss,weight)

           else if(k4(nv).eq.3.and.k5(nv).eq.5.and.k6(nv).eq.3) then  !  Cloud top: conditional mean
                do i =1,iens
                 if(apoint(i).le.0.) apoint(i)=20000.0
                end do
                call get_cond_mean (apoint,iens,20000.0,amean,aspread,
     +                               miss,weight)

           else
             
             call getmean(apoint,iens,amean,aspread,miss,weight)

           end if

            vrbl_mn(igrid,lv)=amean
            vrbl_sp(igrid,lv)=aspread

            if (igrid .eq. jf) then
              write(6,*) 'igrid = jf block'
               write(6,*) 'nv,k4(nv),k5(nv): ', nv,k4(nv),k5(nv)
              
            if (trim(Msignal(nv)).eq.'P') then
               write(*,*) 'shape(rawdata_mn): ', shape(rawdata_mn)
               write(*,*) 'shape(vrbl_mn):', shape(vrbl_mn)
               write(*,*) 'shape(vrbl_mn_pm):', shape(vrbl_mn_pm)

               Lm=max(1,Mlvl(nv))                !Keep at least one vrbl_mn to pass it into packGB2 
!              allocate (vrbl_mn(jf,Lm))         !in case Mlvl(nv)=0
!              allocate (vrbl_mn_pm(jf,Lm))         !in case Mlvl(nv)=0

               !Currently, only for APCP, REFD,REFC and RETOP 
               call pmatch_mean(vname(nv),rawdata_mn,vrbl_mn,
     &             k4(nv),k5(nv),k6(nv),vrbl_mn_pm,lv,lm,jf,iens)

	       vrbl_mn_blend(:,:)=0.5*(vrbl_mn_pm(:,:)+vrbl_mn(:,:))
               write(6,*) 'maxval(vrbl_mn_blend): ', 
     &          maxval(vrbl_mn_blend)
 
              endif

	write(0,*) 'past pmatch_mean with vname(nv): ', vname(nv)

              if (trim(Msignal(nv)).eq.'L') then


	write(0,*) 'AP1h/AP3h branch'
               patch_nx=6
               patch_ny=6
               ovx=patch_nx*5
               ovy=patch_ny*5
               filt_min=0.01
               gauss_sig=2.0

               call pmatch_mean(vname(nv),rawdata_mn,vrbl_mn,
     &             k4(nv),k5(nv),k6(nv),vrbl_mn_pm,lv,lm,jf,iens)


!!! recast what is being passed in here as 2D arrays
                  
               allocate(vrbl_mn_2d(im,jm))
               allocate(vrbl_lpm_2d(im,jm))
               allocate(vrbl_pmmn_2d(im,jm))
               allocate(rawdata_mn_2d(im,jm,iens))

	       do JJ=1,jm
               do II=1,im
                 I1D=(JJ-1)*IM+II
                 vrbl_mn_2d(II,JJ)=vrbl_mn(I1D,1)
                 vrbl_pmmn_2d(II,JJ)=vrbl_mn_pm(I1D,1)
               enddo
               enddo

               do I=1,iens
	       do JJ=1,jm
               do II=1,im
                 I1D=(JJ-1)*IM+II
                 rawdata_mn_2d(II,JJ,I)=rawdata_mn(I1D,I,1)
               enddo
               enddo
               enddo

               call lpm(im,jm,iens,patch_nx,patch_ny,ovx,ovy,
     &            filt_min, gauss_sig,
     &            rawdata_mn_2d,vrbl_mn_2d,vrbl_lpm_2d,vrbl_pmmn_2d)

	       do JJ=1,jm
               do II=1,im
                 I1D=(JJ-1)*IM+II
                 vrbl_mn_locpm(I1D,1)=vrbl_lpm_2d(II,JJ)
               enddo
               enddo

            write(0,*) 'presmooth max(vrbl_mn_locpm): ', 
     &                  maxval(vrbl_mn_locpm)

             call Gsmoothing(vrbl_mn_locpm(:,1),jf,im,jm,
     +           'M','M')

	       vrbl_mn_blendlpm(:,:)=0.5*(vrbl_mn_locpm(:,:)+vrbl_mn(:,:))

	       deallocate(vrbl_mn_2d,vrbl_lpm_2d,rawdata_mn_2d)
               deallocate(vrbl_pmmn_2d)


            write(0,*) 'postsmooth max(vrbl_mn_locpm): ', 
     &                  maxval(vrbl_mn_locpm)

               endif ! lpm variable selection


              endif ! jf block
!             endif

          end do  !end of igrid loop
  
         write(*,'(a4,10f9.2)')'MEAN',(vrbl_mn(i,lv),i=jf/2,jf/2+9)
         write(*,'(a4,10f9.2)')'SPRD',(vrbl_sp(i,lv),i=jf/2,jf/2+9)

         end do  !end if Mlvl

        END IF     !End of Msignal(nv)='M' and 'P' and 'L'

        IF(trim(Msignal(nv)).eq.'X') THEN   ! Max member value

         do lv = 1, Mlvl(nv)
          do igrid=1,jf

           apoint = rawdata_mn(igrid,:,lv)
           miss = missing(nv,:)

           call getmax (apoint,iens,vmax,miss) 
           vrbl_mn(igrid,lv)=vmax

          end do
         end do

        END IF     !End of Msignal(nv)='X'


        IF (trim(Psignal(nv)).eq.'P'.or.    !P: general prob
     +      trim(Psignal(nv)).eq.'A'.or.    !A,B,C,D: prob for different neighborhood
     +      trim(Psignal(nv)).eq.'B'.or.    ! radius size
     +      trim(Psignal(nv)).eq.'C'.or. 
     +      trim(Psignal(nv)).eq.'D'.or.
     +      trim(Psignal(nv)).eq.'H'.or.
     +      trim(Psignal(nv)).eq.'K'.or.
     +      trim(Psignal(nv)).eq.'L'.or.
     +      trim(Psignal(nv)).eq.'M'.or.
     +      trim(Psignal(nv)).eq.'I') THEN


         do lv = 1, Plvl(nv)
          do lt = 1, Tlvl(nv)

             thr1 = Thrs(nv,lt)

	if (vname(nv) .eq. 'A6RI') then
	write(0,*) 'defining fname for A6RI ' , thr1
        if (thr1 .eq. 1) then
          igrb2=71
          fname="new_1y_6h.grib2.rrfs"
        endif

        if (thr1 .eq. 2) then
          igrb2=72
          fname="new_2y_6h.grib2.rrfs"
        endif

        if (thr1 .eq. 5) then
          igrb2=73
          fname="new_5y_6h.grib2.rrfs"
        endif

        if (thr1 .eq. 10) then
          igrb2=74
           fname='new_10y_6h.grib2.rrfs'
        endif

        if (thr1 .eq. 25) then
          igrb2=75
           fname='new_25y_6h.grib2.rrfs'
        endif

        if (thr1 .eq. 50) then
          igrb2=76
           fname='new_50y_6h.grib2.rrfs'
        endif

        if (thr1 .eq. 100) then
          igrb2=77
            fname='new_100y_6h.grib2.rrfs'
        endif


	write(0,*) 'trim(fname): ', trim(fname)
!               igrb2=81
        call baopenr(igrb2,trim(fname),ier1)
	write(0,*) 'ier1 from baopenr for A6RI: ', ier1
              jpdtn=-1
              jpd27=-9999
            jpd1=-9999
            jpd2=-9999
            jpd10=-9999
            jpd12=-9999
            jpd27=-9999
          write(*,*)'a3 - readGB2:',igrb2,jpdtn,
     +                jpd1,jpd2,jpd10,jpd12,jpd27

	gfld_neighb_restore=gfld

          call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld_ri, eps, jret)

	write(0,*) 'jret from readGB2: ', jret
            jpd1=1
            jpd2=8
            jpd10=1
            jpd27=6

	if (jret .eq. 0) then
	do J=1,jf
	if (gfld_ri%fld(j) .gt. 0. .and. gfld_ri%bmap(j)) then
!old	    return_int(J)=gfld_ri%fld(j)*0.0254  ! convert from "milliinches" to mm
            return_int(J)=gfld_ri%fld(j)  ! believe new data is in mm
        else
            return_int(J)=-9999.
	endif
	enddo
	endif

	gfld=gfld_neighb_restore

	call baclose(igrb2,ier1)

	write(0,*) 'min/max of 6h return_int: ', minval(return_int),
     &                maxval(return_int)


        elseif (vname(nv) .eq. 'A3RI') then
	if (thr1 .eq. 1) fname='3h_1yr.grib.href5km'
	if (thr1 .eq. 2) fname='3h_2yr.grib.href5km'
	if (thr1 .eq. 5) fname='3h_5yr.grib.href5km'
	if (thr1 .eq. 10) fname='3h_10yr.grib.href5km'
	if (thr1 .eq. 25) fname='3h_25yr.grib.href5km'
	if (thr1 .eq. 50) fname='3h_50yr.grib.href5km'
	if (thr1 .eq. 100) fname='3h_100yr.grib.href5km'
	if (thr1 .eq. 200) fname='3h_200yr.grib.href5km'
	if (thr1 .eq. 500) fname='3h_500yr.grib.href5km'
	if (thr1 .eq. 1000) fname='3h_1000yr.grib.href5km'

	write(0,*) 'trim(fname): ', trim(fname)
               igrb2=82
        call baopenr(igrb2,trim(fname),ier1)
	write(0,*) 'ier1 from baopenr for A3RI: ', ier1
              jpdtn=-1
              jpd27=-9999
            jpd1=-9999
            jpd2=-9999
            jpd10=-9999
            jpd12=-9999
            jpd27=-9999
          write(*,*)'readGB2:',igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27

	gfld_neighb_restore=gfld

          call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld_ri, eps, jret)

	write(0,*) 'jret from readGB2: ', jret
            jpd1=1
            jpd2=8
            jpd10=1
            jpd27=3

	if (jret .eq. 0) then
	do J=1,jf
	if (gfld_ri%fld(j) .gt. 0. .and. gfld_ri%bmap(j)) then
!old	    return_int(J)=gfld_ri%fld(j)*0.0254  ! convert from "milliinches" to mm
            return_int(J)=gfld_ri%fld(j)  ! believe new data is in mm
        else
            return_int(J)=-9999.
	endif
	enddo
	endif

	gfld=gfld_neighb_restore

	call baclose(igrb2,ier1)

	write(0,*) 'min/max of 3h return_int: ', minval(return_int),
     &                maxval(return_int)

        elseif (vname(nv) .eq. 'A12R' ) then
	if (thr1 .eq. 1) fname='12h_1yr.grib.href5km'
	if (thr1 .eq. 2) fname='12h_2yr.grib.href5km'
	if (thr1 .eq. 5) fname='12h_5yr.grib.href5km'
	if (thr1 .eq. 10) fname='12h_10yr.grib.href5km'
	if (thr1 .eq. 25) fname='12h_25yr.grib.href5km'
	if (thr1 .eq. 50) fname='12h_50yr.grib.href5km'
	if (thr1 .eq. 100) fname='12h_100yr.grib.href5km'
	if (thr1 .eq. 200) fname='12h_200yr.grib.href5km'
	if (thr1 .eq. 500) fname='12h_500yr.grib.href5km'
	if (thr1 .eq. 1000) fname='12h_1000yr.grib.href5km'

	write(0,*) 'trim(fname): ', trim(fname)
               igrb2=83
        call baopenr(igrb2,trim(fname),ier1)
	write(0,*) 'ier1 from baopenr: ', ier1
              jpdtn=-1
              jpd27=-9999
            jpd1=-9999
            jpd2=-9999
            jpd10=-9999
            jpd12=-9999
            jpd27=-9999
          write(*,*)'readGB2:',igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27

	gfld_neighb_restore=gfld

          call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld_ri, eps, jret)

	write(0,*) 'jret from readGB2: ', jret
            jpd1=1
            jpd2=8
            jpd10=1
            jpd27=12

	if (jret .eq. 0) then
	do J=1,jf
	if (gfld_ri%fld(j) .gt. 0. .and. gfld_ri%bmap(j)) then
!old	    return_int(J)=gfld_ri%fld(j)*0.0254  ! convert from "milliinches" to mm
            return_int(J)=gfld_ri%fld(j)  ! believe new data is in mm
        else
            return_int(J)=-9999.  ! set large so will not be exceeded
	endif
	enddo
	endif

	gfld=gfld_neighb_restore
	call baclose(igrb2,ier1)

	write(0,*) 'min/max of 12h return_int: ', minval(return_int),
     &                maxval(return_int)

        elseif (vname(nv) .eq. 'A24R' ) then

        if (thr1 .eq. 1) then 
         igrb2=84
         fname="new_1y_24h.grib2.rrfs"
        endif

        if (thr1 .eq. 2) then
         igrb2=85
         fname="new_2y_24h.grib2.rrfs"
        endif

        if (thr1 .eq. 5) then
         igrb2=86
         fname="new_5y_24h.grib2.rrfs"
        endif

        if (thr1 .eq. 10) then
         igrb2=87
         fname="new_10y_24h.grib2.rrfs"
        endif

        if (thr1 .eq. 25) then
         igrb2=88
         fname="new_25y_24h.grib2.rrfs"
        endif

        if (thr1 .eq. 50) then
         igrb2=89
         fname="new_50y_24h.grib2.rrfs"
        endif

        if (thr1 .eq. 100) then
         igrb2=90
         fname="new_100y_24h.grib2.rrfs"
        endif

	write(0,*) 'trim(fname): ', trim(fname)
!               igrb2=84
 	write(0,*) 'open as igrb2: ', igrb2
        call baopenr(igrb2,trim(fname),ier1)
	write(0,*) 'ier1 from baopenr: ', ier1
              jpdtn=-1
              jpd27=-9999
            jpd1=-9999
            jpd2=-9999
            jpd10=-9999
            jpd12=-9999
            jpd27=-9999
          write(*,*)'readGB2:',igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27

	gfld_neighb_restore=gfld

          call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld_ri, eps, jret)

	write(0,*) 'jret from readGB2 for 24hRI: ', jret
            jpd1=1
            jpd2=8
            jpd10=1
            jpd27=24

	if (jret .eq. 0) then
	do J=1,jf
	if (gfld_ri%fld(j) .gt. 0. .and. gfld_ri%bmap(j)) then
!old	    return_int(J)=gfld_ri%fld(j)*0.0254  ! convert from "milliinches" to mm
            return_int(J)=gfld_ri%fld(j)  ! believe new data is in mm
        else
            return_int(J)=-9999.  ! set large so will not be exceeded
	endif
	enddo
	endif

	gfld=gfld_neighb_restore
	call baclose(igrb2,ier1)

	write(0,*) 'min/max of 24 h return_int: ', minval(return_int),
     &                maxval(return_int)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	elseif (vname(nv) .eq. 'FFG1') then
	write(0,*) 'defining fname for FFG1 ' , thr1
	if (thr1 .eq. 1) then
!          fname="ffg1h.grib2.href5km"
          fname="href.ffg1h.5km.grib2"
        else
          write(0,*) 'ffg1h with wrong thr1 ', thr1
        endif
	write(0,*) 'trim(fname): ', trim(fname)
               igrb2=91
        call baopenr(igrb2,trim(fname),ier1)
	write(0,*) 'ier1 from baopenr for ffg1: ', ier1
              jpdtn=-1
              jpd27=-9999
            jpd1=-9999
            jpd2=-9999
            jpd10=-9999
            jpd12=-9999
            jpd27=-9999
          write(*,*)'readGB2 for ffg1:',igrb2,jpdtn,jpd1,
     +              jpd2,jpd10,jpd12,jpd27

	gfld_neighb_restore=gfld

          call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld_ffg, eps, jret)

	write(0,*) 'jret from readGB2 for ffg1: ', jret
            jpd1=1
            jpd2=194
            jpd10=1
            jpd27=1

!           if (Tsignal(nv).ne.'A') then

	if (jret .eq. 0) then
	do J=1,jf
	if (gfld_ffg%fld(j) .gt. 0.01 .and. gfld_ffg%bmap(j)) then
	    return_int(J)=gfld_ffg%fld(j)
            if (j .eq. 398683) then
            write(0,*) 'ffg1 j, return_int(J): ', j, return_int(J)
            endif
        else
            return_int(J)=-9999.
	endif
	enddo
	endif

	gfld=gfld_neighb_restore
        gfld%discipline=1   ! FFG exceedance is a hydro discipline product

	call baclose(igrb2,ier1)
        write(0,*) 'ier1 from baclose for ffg1: ', ier1

	write(0,*) 'min/max of FFG1: ', minval(return_int),
     &                maxval(return_int)





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	elseif (vname(nv) .eq. 'FFG3') then
	write(0,*) 'defining fname for FFG3 ' , thr1
!	if (thr1 .eq. 3) fname="ffg3h.grib2.href5km"
	if (thr1 .eq. 3) fname="href.ffg3h.5km.grib2"

	write(0,*) 'trim(fname): ', trim(fname)
               igrb2=92
        call baopenr(igrb2,trim(fname),ier1)
	write(0,*) 'ier1 from baopenr for ffg3: ', ier1
              jpdtn=-1
              jpd27=-9999
            jpd1=-9999
            jpd2=-9999
            jpd10=-9999
            jpd12=-9999
            jpd27=-9999
          write(*,*)'readGB2:',igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27

	gfld_neighb_restore=gfld

          call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld_ffg, eps, jret)

            jpd1=1
            jpd2=194

	write(0,*) 'jret from readGB2: ', jret
            jpd10=1
            jpd27=3

!           if (Tsignal(nv).ne.'A') then

	if (jret .eq. 0) then
	do J=1,jf
	if (gfld_ffg%fld(j) .gt. 0.01 .and. gfld_ffg%bmap(j)) then
	    return_int(J)=gfld_ffg%fld(j)
            if (j .eq. 398683) then
            write(0,*) 'ffg3 j, return_int(J): ', j, return_int(J)
            endif
        else
            return_int(J)=-9999.
	endif
	enddo
	endif

	gfld=gfld_neighb_restore

        gfld%discipline=1   ! FFG exceedance is a hydro discipline product
	call baclose(igrb2,ier1)

	write(0,*) 'min/max of FFG3: ', minval(return_int),
     &                maxval(return_int)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	elseif (vname(nv) .eq. 'FFG6') then
	write(0,*) 'defining fname for FFG6 ' , thr1
!	if (thr1 .eq. 6) fname="ffg6h.grib2.href5km"
	if (thr1 .eq. 6) fname="href.ffg6h.5km.grib2"

	write(0,*) 'trim(fname): ', trim(fname)
               igrb2=93
        call baopenr(igrb2,trim(fname),ier1)
	write(0,*) 'ier1 from baopenr: ', ier1
              jpdtn=-1
              jpd27=-9999
            jpd1=-9999
            jpd2=-9999
            jpd10=-9999
            jpd12=-9999
            jpd27=-9999
          write(*,*)'readGB2:',igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27

	gfld_neighb_restore=gfld

          call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld_ffg, eps, jret)


! should we look at bitmap of the file??

	write(0,*) 'jret from readGB2: ', jret
            jpd1=1
            jpd2=194
            jpd10=1
            jpd27=6

	if (jret .eq. 0) then
	do J=1,jf
	if (gfld_ffg%fld(j) .gt. 0.01 .and. gfld_ffg%bmap(j)) then
	    return_int(J)=gfld_ffg%fld(j)
        else
            return_int(J)=-9999.
	endif
	enddo
	endif

	gfld=gfld_neighb_restore
        gfld%discipline=1   ! FFG exceedance is a hydro discipline product

	call baclose(igrb2,ier1)

	write(0,*) 'min/max of FFG6: ', minval(return_int),
     &                maxval(return_int)

	elseif (vname(nv) .eq. 'FF12') then
	write(0,*) 'defining fname for FF12 ' , thr1
	if (thr1 .eq. 12) fname="ffg12h.grib2.href5km"

	write(0,*) 'trim(fname): ', trim(fname)
               igrb2=94
        call baopenr(igrb2,trim(fname),ier1)
	write(0,*) 'ier1 from baopenr: ', ier1
              jpdtn=-1
              jpd27=-9999
            jpd1=-9999
            jpd2=-9999
            jpd10=-9999
            jpd12=-9999
            jpd27=-9999
          write(*,*)'readGB2:',igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27

	gfld_neighb_restore=gfld

          call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld_ffg, eps, jret)

	write(0,*) 'jret from readGB2: ', jret
            jpd1=1
            jpd2=8
            jpd10=1
            jpd27=12

	if (jret .eq. 0) then
	do J=1,jf
	if (gfld_ffg%fld(j) .gt. 0.01 .and. gfld_ffg%bmap(j)) then
	    return_int(J)=gfld_ffg%fld(j)
        else
            return_int(J)=-999999.
	endif
	enddo
	endif

	gfld=gfld_neighb_restore

	call baclose(igrb2,ier1)

	write(0,*) 'min/max of FF12: ', minval(return_int),
     &                maxval(return_int)

	elseif (vname(nv) .eq. 'FF24') then
	write(0,*) 'defining fname for FF24 ' , thr1
	if (thr1 .eq. 24) fname="ffg24h.grib2.href5km"

	write(0,*) 'trim(fname): ', trim(fname)
               igrb2=95
        call baopenr(igrb2,trim(fname),ier1)
	write(0,*) 'ier1 from baopenr: ', ier1
              jpdtn=-1
              jpd27=-9999
            jpd1=-9999
            jpd2=-9999
            jpd10=-9999
            jpd12=-9999
            jpd27=-9999
          write(*,*)'readGB2:',igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27

	gfld_neighb_restore=gfld

          call readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpd27,
     +          gfld_ffg, eps, jret)

	write(0,*) 'jret from readGB2: ', jret
            jpd1=1
            jpd2=8
            jpd10=1
            jpd27=24

	if (jret .eq. 0) then
	do J=1,jf
	if (gfld_ffg%fld(j) .gt. 0.01 .and. gfld_ffg%bmap(j)) then
	    return_int(J)=gfld_ffg%fld(j)
        else
            return_int(J)=-999999.

	endif
	enddo
	endif

	gfld=gfld_neighb_restore

	call baclose(igrb2,ier1)

	write(0,*) 'min/max of FF24: ', minval(return_int),
     &                maxval(return_int)


	endif

           do igrid=1,jf

	if (vname(nv) .eq. 'SN1h' .or. vname(nv) .eq. 'SN3h' .or. 
     &      vname(nv) .eq. 'SN6h') then
            apoint = rawdata_pr(igrid,:,lv)*1. ! SLR val
       else
            apoint = rawdata_pr(igrid,:,lv) 
       endif
            miss = missing(nv,:)


          if ((k4(nv).eq.2.and.k5(nv).eq.299).or.   !These field not in hrrrgsd file
     +        (k4(nv).eq.3.and.k5(nv).eq.1).or.
     +        (k4(nv).eq.1.and.k5(nv).eq.192).or.
     +        (k4(nv).eq.6.and.k5(nv).eq.1).or.
     +        (k4(nv).eq.7.and.k5(nv).eq.197) ) then
             do ie=1,iens
               if(mbrname(ie)(1:4).eq.'hrrr') then
                miss(ie)=1
               end if
             end do
          end if


            if(k4(nv).eq.19.and.k5(nv).eq.0) then  ! Visibility: conditional mean 
             do ie =1,iens
               if(apoint(ie).lt.1.0) miss(ie)=1    !exclude spc runs's error vis 
             end do
            end if

            if(k4(nv).eq.3.and.k5(nv).eq.5.and.k6(nv).eq.2) then  !  Cloud base: conditional mean
                do i =1,iens
                 if(apoint(i).le.0.) apoint(i)=20000.0
                end do
            end if
            if(k4(nv).eq.3.and.k5(nv).eq.5.and.k6(nv).eq.3) then  !  Cloud top: conditional mean
                do i =1,iens
                 if(apoint(i).le.0.) apoint(i)=20000.0
                end do
            end if


            if(trim(op(nv)).ne.'-') then

             thr1 = Thrs(nv,lt)
             thr2 = 0.
             if (vname(nv) .ne. 'A3RI' .and. vname(nv) .ne. 'A6RI'.and.
     &           vname(nv) .ne. 'A24R' .and. vname(nv) .ne. 'FFG1'.and.
     &           vname(nv) .ne. 'FFG3' .and.
     &           vname(nv) .ne. 'FFG6' .and. vname(nv) .ne. 'FF24'.and.
     &           vname(nv) .ne. 'FF12' .and. vname(nv) .ne. 'A12R')then
             call getprob(apoint,iens,thr1,thr2,op(nv),aprob,
     +                         miss,weight)
             vrbl_pr(igrid,lv,lt)=aprob

             else

             if (igrid .eq. jf/2) then
             write(0,*) 'working on : ', vname(nv)
             endif


! thr1 needs to come from the grid point of the RI data (and FFG data...using same return_int array)
! 
! where should that data be read?  thr1 will have the year period
! over which the RI is computed

	     thr1=return_int(igrid)
             if (mod(igrid,100000) .eq. 0) then
               write(0,*) 'for igrid: ', igrid
               write(0,*) 'working using thr1: ', thr1
               write(0,*) 'maxval(apoint): ', maxval(apoint)
               write(0,*) 'op(nv): ', op(nv)

              endif
             thr2=0.
	if (thr1 .gt. 0 .and. thr1 .le. 99999.) then
             call getprob(apoint,iens,thr1,thr2,op(nv),aprob,
     +                         miss,weight)
             vrbl_pr(igrid,lv,lt)=aprob
        else
             vrbl_pr(igrid,lv,lt)=-10.
        endif

             endif

            else

             if(lt.lt.Tlvl(nv)) then
               thr1 = Thrs(nv,lt)
               thr2 = Thrs(nv,lt+1)
               call getprob(apoint,iens,thr1,thr2,op(nv),aprob,
     +                      miss,weight)
               vrbl_pr(igrid,lv,lt)=aprob
             end if

            end if

           end do ! end of igrid loop

           !Check if the probability needs to smooting? 
           if (Tsignal(nv).ne.'T') then

             !write(*,*) 'before smoothing',vrbl_pr(797267,lv,lt)

	write(0,*) 'Tsignal(nv) at gsmoothing vrbl_pr: ',
     &         nv,Tsignal(nv)
          
             call Gsmoothing(vrbl_pr(:,lv,lt),jf,im,jm,
     +           Psignal(nv),Tsignal(nv))

             !write(*,*) 'after smoothing',vrbl_pr(797267,lv,lt)

           end if 

           if(vname(nv).eq.'AP3h') then    !This var will be used in cptp nightning/thunder computation
            if(thr1.eq.0.25) then
             p03mp01(:)=vrbl_pr(:,lv,lt)
            end if
           end if

          end do !end of Tlvl
         end do !end of Plvl

        END IF ! End of Psignal(nv) = 'P,A,B,C,D,H,I' 

        IF (trim(Psignal(nv)).eq.'E'.or.        !Neighborhood fraction prob, where 
     +      trim(Psignal(nv)).eq.'F') THEN      ! E, F for different neighborhood radius size

         do lv = 1, Plvl(nv)
          do lt = 1, Tlvl(nv)
            write(*,*) 'Call  neighborhood fraction .......'
 
             if(trim(op(nv)).ne.'-') then
               thr1 = Thrs(nv,lt)
               thr2 = 0.
             else
               thr1 = Thrs(nv,lt)
               thr2 = Thrs(nv,lt+1)
             end if

             do irun=1,iens
c               write(*,*) '398918 before', rawdata_pr(398918,irun,lv)
c               write(*,*) jf,im,jm,Psignal(nv),thr1,thr2,op(nv) 
              call neighborhood_fraction(rawdata_pr(:,irun,lv),
     +         jf,im,jm,Psignal(nv),thr1,thr2,op(nv))
c               write(*,*) '398918 after', rawdata_pr(398918,irun,lv)
             end do

             

             !now rawdata_pr stores fraction at each grid 
             do igrid=1,jf
              sm=0.
              do irun=1,iens
               sm=sm+rawdata_pr(igrid,irun,lv)
              end do
               vrbl_pr(igrid,lv,lt)=100.0*sm/iens
c               if(igrid.eq.398918) then
c                write(*,*) igrid, rawdata_pr(igrid,:,lv),
c     +                     vrbl_pr(igrid,lv,lt)
c               end if
             end do

             if (Tsignal(nv).ne.'T') then
               call Gsmoothing(vrbl_pr(:,lv,lt),jf,im,jm,
     +               Psignal(nv),Tsignal(nv))
             end if

          end do
         end do

        END IF     ! End of Psignal(nv) = 'E,F'


        IF (trim(Psignal(nv)).eq.'S') THEN    !get spaghetti plots according to SPC SSEO

         do lv = 1, Plvl(nv)
          do lt = 1, Tlvl(nv)

           do igrid=1,jf

            apoint = rawdata_pr(igrid,:,lv)
            miss = missing(nv,:)
            thr1 = Thrs(nv,lt)
            thr2 = 0.
             
            call getSP(apoint,iens,thr1,thr2,op(nv),sp,
     +                         miss)
            vrbl_pr(igrid,lv,lt)=sp

           end do
         
          end do
         end do

       END IF  !! End of Psignal(nv) = 'S' 


        write(*,*) 'Ensemble computation done for direct var ', nv
c Loop 1-3:  Packing  mean/spread/prob for this direct variable

! Reset gfld%bmap with the combined version bmap_f
        if(trim(eps).eq.'href' .and. associated(gfld%bmap)) then
		gfld%bmap=bmap_f
        else
        write(0,*) 'did not set gfld%bmap to bmap_f as not associated'
        endif

!        if (trim(Msignal(nv)).eq.'M'.or.trim(Msignal(nv)).eq.'P') then
! change so PM mean field doesn't write out regular mean as well.


        if (trim(Msignal(nv)).eq.'M') then

	write(0,*) 'calling packGB2_mean for M 1,2,27 ', jpd1,jpd2, jpd27


! reset EMSL to write out as PMSL?
!	if (jpd1 .eq. 3 .and. jpd2 .eq. 192) then
!        jpd2=1
!        endif

	write(0,*) 'here with maxval(vrbl_mn(:,Lm)): ', maxval(vrbl_mn(:,Lm))
          call packGB2_mean(imean,isprd,vrbl_mn,vrbl_sp,   !jpd12 is determined inside 
     +          nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +          iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)    !gfld is used to send in other info 
          write(*,*) 'packing mean direct var for', nv

	if (jpd2 .eq. 13 .and. jpd27 .eq. 1)  then
	write(0,*) 'defining vrbl_mn_use'
	write(0,*) 'Lm is: ', Lm
	write(0,*) 'maxval(vrbl_mn(:,Lm)): ', maxval(vrbl_mn(:,Lm))
	
          jpd2loc=11
          do JJJ=1,jf
	  vrbl_mn_use(JJJ)=vrbl_mn(JJJ,Lm)
          enddo
	 
	write(0,*) 'min/max of test SNOW field ', minval(vrbl_mn_use),
     +                                       maxval(vrbl_mn_use)

!          call packGB2_mean(imean,isprd,vrbl_mn_use,vrbl_sp,   !jpd12 is determined inside 
!     +          nv,jpd1,jpd2loc,jpd10,jpd27,jf,Lm,
!     +          iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)    !gfld is used to send in other info 

        endif

	endif
        if(trim(Msignal(nv)).eq.'L') then
! LPMM
          call packGB2_mean(ilocpmmn,isprd,vrbl_mn_locpm,vrbl_sp,
     +          nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +          iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)         

! AVRG
! limit avrg to precip

!          if (vname(nv).eq.'AP1h' .or. vname(nv).eq.'AP3h' .or. 
!     &        vname(nv).eq.'AP6h' .or. vname(nv).eq.'AP24') then
!          call packGB2_mean(ilocavg,isprd,vrbl_mn_blendlpm,vrbl_sp,    
!     +          nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
!     +          iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)         
!          endif

        endif

        if(trim(Msignal(nv)).eq.'P') then
!          jpd10=10                                         !specific for probablity-match meaned QPF


! PMMN
          call packGB2_mean(ipmmn,isprd,vrbl_mn_pm,vrbl_sp,    
     +          nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +          iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)         
          write(*,*) 'packing PM mean direct var for', nv

! AVRG
! limit avrg to precip

          if (vname(nv).eq.'AP1h' .or. vname(nv).eq.'AP3h' .or. 
     &        vname(nv).eq.'AP6h' .or. vname(nv).eq.'AP24') then
          call packGB2_mean(iavg,isprd,vrbl_mn_blend,vrbl_sp,    
     +          nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +          iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)         
          endif

! MEAN (just for precip with PM type fields)

          if (vname(nv).eq.'AP1h' .or. vname(nv).eq.'AP3h' .or. 
     &        vname(nv).eq.'AP6h' .or. vname(nv).eq.'AP24') then

         write(0,*) 'make sure it gets packed like a true mean field'
         Msignal(nv)='M'
          call packGB2_mean(imean,isprd,vrbl_mn,vrbl_sp,   !jpd12 is determined inside 
     +          nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +          iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)    !gfld is used to send in other info 
          endif
         Msignal(nv)='P'

        end if

        if(trim(Psignal(nv)).eq.'P'.or.
     +     trim(Psignal(nv)).eq.'A'.or.
     +     trim(Psignal(nv)).eq.'B'.or.
     +     trim(Psignal(nv)).eq.'C'.or.
     +     trim(Psignal(nv)).eq.'D'.or.
     +     trim(Psignal(nv)).eq.'E'.or.
     +     trim(Psignal(nv)).eq.'F'.or.
     +     trim(Psignal(nv)).eq.'H'.or. 
     +     trim(Psignal(nv)).eq.'K'.or. 
     +     trim(Psignal(nv)).eq.'L'.or. 
     +     trim(Psignal(nv)).eq.'M'.or. 
     +     trim(Psignal(nv)).eq.'I'  ) then


       if (vname(nv) .eq. 'A3RI' .or. vname(nv) .eq. 'A6RI' .or. 
     &     vname(nv) .eq. 'A12R' .or. vname(nv) .eq. 'A24R' .or.
     &     vname(nv) .eq. 'FFG1' .or.
     &     vname(nv) .eq. 'FFG3' .or. vname(nv) .eq. 'FFG6' .or.
     &     vname(nv) .eq. 'FF12' .or. vname(nv) .eq. 'FF24')  then

	write(0,*) 'calling packGB2_prob for FFG/ARI'
        write(0,*) 'writing out: ', vname(nv)
        write(0,*) 'max of vrbl_pr: ', maxval(vrbl_pr)


         call packGB2_prob(iffri,vrbl_pr,             !jpd12 is determined inside
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lth,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)

        else

	write(0,*) 'calling packGB2_prob for normal'
        write(0,*) 'writing out: ', vname(nv)
        write(0,*) 'max of vrbl_pr: ', maxval(vrbl_pr)

	if (vname(nv) .eq. 'LTNG') then
            gfld=gfld_lightning
        endif

         call packGB2_prob(iprob,vrbl_pr,             !jpd12 is determined inside
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lth,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)

           write(*,*) 'packed prob direct var for', nv

        endif

        end if

        if (trim(Msignal(nv)).eq.'X') then
          !SSEO: Mean and max are using same packGB2_max code
          !    and put in same mean file, no spread file 
              if (vname(nv).eq.'UH3h'.or.
     +            vname(nv).eq.'RF3h'.or.
     +            vname(nv).eq.'WS3h'.or.
     +            vname(nv).eq.'US3h') jpd27=3
              if (vname(nv).eq.'UH24'.or.
     +            vname(nv).eq.'RF24'.or.
     +            vname(nv).eq.'WS24'.or.
     +            vname(nv).eq.'US24') jpd27=24


          call packGB2_max(imean,isprd,vrbl_mn,              !jpd12 is determined inside
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)           !gfld is used to send in other info

           write(*,*) 'packing max direct var for', nv
        end if

        if (trim(Psignal(nv)).eq.'S') then
         call packGB2_spag(iprob,vrbl_pr,             !jpd12 is determined inside
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lth,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)

           write(*,*) 'packing spaghetii  direct var for', nv
        end if

      end do loop222

c Loop 1-4: Deallocation

        deallocate (rawdata_mn)
        deallocate (rawdata_pr)
        deallocate (vrbl_mn)
        deallocate (vrbl_mn_pm)
        deallocate (vrbl_mn_locpm)
        deallocate (vrbl_mn_blend)
        deallocate (vrbl_mn_blendlpm)
        deallocate (vrbl_sp)
        deallocate (vrbl_pr)
        deallocate (return_int)

        write(*,*) 'Variable # ', nv, ' done ....'
      
2001     continue  !end of nv direct variables loop



C --------------------------------------------------
C --------------------------------------------------
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C --------------------------------------------------
C --------------------------------------------------
C
c  Loop 2: for derived variables
c    2-0: allocate necessary arrays
c    2-1: compute all derived variable by calling their correspondent subroutines 
c            Note: all input data will be read inside each subroutine 
c    2-2: pack derived varaibles
c    2-3: deallocation

       write(*,*) 'Now begin to generate derived products ......'

       DO 2002 nv = 1, nderiv

c   Loop 2-0: allocation
         gfld%discipline=0   !reset discipline# in case it was changed (i.e. fire weather)

         jpd1=dk4(nv)
	write(0,*) 'jpd1 defined(c): ', jpd1
         jpd2=dk5(nv)
         jpd10=dk6(nv)
         jpd27=-9999

         if(trim(eps).eq.'sref' ) then
           if (jpd1.eq.1.and.
     +      (jpd2.eq.8.or.jpd2.eq.11.or.jpd2.eq.15)) then
              jpdtn=11
           else
              jpdtn=1
           end if
         else
            if (jpd1.eq.1.and.
     +      (jpd2.eq.8.or.jpd2.eq.11.or.jpd2.eq.15)) then
              jpdtn=8
           else
              jpdtn=0
           end if
         end if

          Lm=max(1,dMlvl(nv))
          Lp=max(1,dPlvl(nv))
          Lth=max(1,dTlvl(nv))
          write(*,*) 'dTlvl(nv)=',dTlvl(nv),Lth

          if (.NOT.allocated(derv_mn)) then
            allocate (derv_mn(jf,Lm))         
          end if
          if (.NOT.allocated(derv_sp)) then
            allocate (derv_sp(jf,Lm))
          end if
          if (.NOT.allocated(derv_pr)) then
            allocate (derv_pr(jf,Lp,Lth))
          end if

c Loop 2-1: Computation 

cc%%%%%%% 1.mTo  see if there is thickness computation, if yes, do it
c           if (dk5(nv).eq.7.and. dk6(nv).eq.101) then
c


         
cc%%%%%%% 2. To see if there is precipitation type computation, if yes, do it
          if (dk4(nv).eq.1.and.dk5(nv).eq.19) then

              if (.NOT.allocated(ptype_mn)) then
                allocate (ptype_mn(jf,4))
              end if
              if (.NOT.allocated(ptype_pr)) then
                allocate (ptype_pr(jf,4))
              end if
              if (.NOT.allocated(ptype_pr2)) then
                allocate (ptype_pr2(jf,4))
              end if
             
             if (trim(eps).eq.'sseo') then

              call preciptype (nv,ipunit,jpdtn,jf,iens,
     +         ptype_mn,ptype_pr,ptype_pr2)

             else
 
              call preciptype (nv,ifunit,jpdtn,jf,iens,
     +         ptype_mn,ptype_pr,ptype_pr2)

             end if

              do jp=1,4
               jpd1=1
               jpd2=jptyp2(jp)
               jpd10=1
               jpd12=0
               jpd27=-999
               derv_mn(:,1)=ptype_mn(:,jp)     !for precip type, Lm,lp,Lth all are 1
               derv_sp=0.0

! ptype_pr2 has out of members defining a type
!               derv_pr(:,1,1)=ptype_pr2(:,jp)
! ptype_pr has out of all members 
               derv_pr(:,1,1)=ptype_pr(:,jp)

               gfld_temp=gfld

               if(trim(eps).eq.'href' .and. 
     +           associated(gfld_temp%bmap)) gfld_temp%bmap=bmap_f

               call packGB2_mean_derv(imean,isprd,derv_mn,
     +              derv_sp,nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +              iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld_temp)  !gfld uses what was got from previous direct variables 
      
              gfld_temp=gfld                           !some of idrtmpl() fields have been changed after packGB2_prob,

              if(trim(eps).eq.'href' .and. 
     +          associated(gfld_temp%bmap)) gfld_temp%bmap=bmap_f

               call packGB2_prob_derv(iprob,derv_pr,
     +              nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lth,
     +              iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld_temp)  !gfld uses what was got from previous direct variables

              end do

              deallocate (ptype_mn)
              deallocate (ptype_pr)
              deallocate (ptype_pr2)

              write(*,*) 'preciptype done!'
            end if   !end of preciptype computation and save

cc%%%%%%% 3. To see if there is wind speed computation, if yes, do it

	write(0,*) 'here dk4(nv) dk5(nv) dk6(nv): ', dk4(nv),
     &             dk5(nv),dk6(nv)

          if (dk4(nv).eq.2.and.dk5(nv).eq.1.and.dk6(nv).ne.108) then
            write(*,*) 'call wind'
            call wind (nv,ifunit,jpdtn,jf,im,jm,iens,Lm,Lp,Lth,
     +        derv_mn,derv_sp,derv_pr,weight,mbrname)
            write(*,*) 'Wind done'
          else if (dk4(nv).eq.2 .and. dk5(nv).eq.192 .and.
     +             dk6(nv).eq.103) then
            write(*,*) 'call bulkshear'
            call bulkshear(nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lth,
     +        derv_mn,derv_sp,derv_pr,weight,mbrname)
            write(*,*) 'Wind(shear) done'
          end if


cc%%%%%%% 4. To see if there is icing computation, if yes, do it
          if (dk4(nv).eq.19.and.dk5(nv).eq.7) then
            call get_icing (nv,ifunit,jpdtn,jf,iens,Lp,Lth,
     +              derv_pr,weight)
            write(*,*) 'Icing done'
          end if

cc%%%%%%% 5. To see if there is CAT computation, if yes, do it
          if (dk4(nv).eq.19.and.dk5(nv).eq.22) then
            call get_cat (nv,ifunit,jpdtn,jf,iens,Lp,Lth,
     +              derv_pr,im,jm,dx,dy,weight)
            dPlvl(nv)=dPlvl(nv)-1
            write(*,*) 'CAT done' 
          end if

cc%%%%%%% 6. To see if there is flight restriction  computation, if yes, do it
          if (dk4(nv).eq.19.and.dk5(nv).eq.205) then
            call  flight_res (nv,ifunit,jpdtn,jf,iens,Lp,Lth,
     +              derv_pr,weight)
             write(*,*) 'Flight restiction done'
          end if

cc%%%%%%% 7. To see if there is Hains index for fire weather computation, if yes, do it
          if (dk4(nv).eq.4.and.dk5(nv).eq.2) then
            dTlvl(nv)=dTlvl(nv)-1
            Lth=Lth-1                     !since dop='-'
            call fire_weather (nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lth,
     +            derv_mn,derv_sp,derv_pr,weight)
            
            gfld%discipline=2      !Fireweather discipline = 2, used for packing
            write(*,*) 'Hains Index done'
! Reset gfld%bmap with the combined version bmap_f
        if(trim(eps).eq.'href' .and. associated(gfld%bmap)) 
     +          gfld%bmap=bmap_f
          end if

cc%%%%%%% 7-1. To see if there is Fosberg index fire weather computation, if yes, do it
          if (dk4(nv).eq.4.and.dk5(nv).eq.4) then
             call  get_fosberg(nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lth,
     +         derv_mn,derv_sp,derv_pr,weight)

             gfld%discipline=2      !Fireweather discipline = 2, used for packing
             write(*,*) 'Fosberg Index done'
          end if

cc%%%%%%% 8. To see if there is ceiling computation, if yes, do it

          if (dk4(nv).eq.3.and.dk5(nv).eq.5.and.dk6(nv).eq.215) then
            write(0,*) 'call getceiling'
            call  getceiling (nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lth,
     +             derv_mn,derv_sp,derv_pr,weight)

            write(0,*) 'getceiling done'
          end if


cc%%%%%%% 9. To see if there is fog  computation, if yes, do it

          if(dk4(nv).eq.6.and.dk5(nv).eq.193.and.dk6(nv).eq.103.
     +                             and.itime.ge.2) then

          call new_fog(nv,ifunit,ipunit,jpdtn,jf,im,jm,dx,dy,interval,
     +      iens,Lm,Lp,Lth,derv_mn,derv_sp,derv_pr,weight)
          write(*,*) 'new_fog done'

         end if

cc%%%%%%% 10. To see if there is thickness computation, if yes, do it
          if(dk4(nv).eq.3.and.dk5(nv).eq.5.and.
     +                             dk6(nv).eq.101) then
              
            call thickness (nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lth,
     +             derv_mn,derv_sp,derv_pr,weight)

              write(*,*) 'Thickness done'
          
          end if


cc%%%%%%% 11. To see if there is LLWS computation, if yes, do it
        if(dk4(nv).eq.2.and.dk5(nv).eq.192.and.
     +                             dk6(nv).eq.1) then


         !write(*,*) 'call llws ', eps       !??? Can not write to std !output   

         call llws (nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lth,eps,
     +          derv_mn,derv_sp,derv_pr,weight)

         write(*,*) 'LLWS done'

        end if

cc%%%%%%% 11. To see if there is CONVECTION  computation, if yes, do it
        if(dk4(nv).eq.1.and.dk5(nv).eq.196.and.
     +                             dk6(nv).eq.200) then
         call getconv(nv,ipunit,jpdtn,jf,im,jm,est,iens,Lm,Lp,Lth,
     +      gribid,derv_pr,weight)

           write(*,'(a9,10f9.2)')'CONV PROB',
     +           (derv_pr(i,1,1),i=10001,10010)
            write(*,*) 'getconv done' 

        end if


cc%%%%%%% 12. To see if there is Lighning, if yes, do it
c         missing() array is not considered, so current NSSE won't do this
         if(dk4(nv).eq.17.and.dk5(nv).eq.192.and.
     +                             dk6(nv).eq.1) then
         call get_cptp_severe(nv,cycle(ihr+1),cfhr,ifunit,p03mp01,
     +       im,jm,km,jpdtn,jf,iens,Lp,Lth,derv_pr,weight,eps)
           write(*,*) 'CPTP done '
         end if


cc%%%%%%% 13. To see if there is 850-300mb mean wind, if yes, do it
         if(dk4(nv).eq.2.and.dk5(nv).eq.1.and.
     +                             dk6(nv).eq.108) then

           call meanwind(nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lth,eps,
     +          derv_pr,weight,mbrname)

            write(*,*) '850-300mb mean-wind done'
         end if

cc%%%%%%% 14. To see if there is 10m wind and 2m RH joint prob, if yes, do it
          if(dk4(nv).eq.19.and.dk5(nv).eq.235) then
            dTlvl(nv)=dTlvl(nv)-1
            Lth = Lth -1 
            call get_wind_rh_joint_prob(nv,ifunit,jpdtn,
     +          jf,iens,Lp,Lth,derv_pr,weight)

            write(*,*) 'Wind-RH joint prob done'
          end if

cc%%%%%%% 15. To see if there is snowfall for HREF
            if(dk4(nv).eq.1.and.dk5(nv).eq.13.and.
     +       trim(eps).eq.'href') then
	write(0,*) 'call snow4href'
             call snow4href(nv,ifunit,ipunit,jf,iens,Lm,Lp,Lt,
     +         derv_mn,derv_sp,derv_pr,eps)
               write(0,*) 'Snowfall done'
            end if

cc%%%%%%% 16. To see if there is 700-500mb mean omeg, if yes, do it
         if(dk4(nv).eq.2.and.dk5(nv).eq.8.and.
     +                             dk6(nv).eq.108) then

	write(0,*) 'call meanomeg '
           call meanomeg(nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lth,eps,
     +          derv_mn,derv_sp,weight,mbrname)
            write(*,*) '700-500mb mean-omeg done'
         end if


C  Loop 2-2:  Pack Non-precip type dereved mean/spread/prob products
C
          if(dk4(nv).eq.1.and.dk5(nv).eq.19 ) then  !Precip type prob is packed already
! 
! already okay
!             goto 20021
!
          else                                       !Non-precip type
               jpd1=dk4(nv)
	write(0,*) 'jpd1 defined(d): ', jpd1
               jpd2=dk5(nv)
               jpd10=dk6(nv)
               jpd27=-9999


           !Check if the probability needs to smooting? 

	write(0,*) 'dTsignal(nv) at gsmoothing derv_pr: ',
     &         nv,dTsignal(nv)
	write(0,*) 'Psignal, dPsignal: ', Psignal(nv),dPsignal(nv)
           if (dTsignal(nv).ne.'T') then

         do lv = 1, dPlvl(nv)
          do lt = 1, dTlvl(nv)

         write(*,*) 'before smoothing of derv_pr ',derv_pr(jf/2,lv,lt)

             call Gsmoothing(derv_pr(:,lv,lt),jf,im,jm,
     +           dPsignal(nv),dTsignal(nv))

          write(*,*) 'after smoothing of derv_pr ',derv_pr(jf/2,lv,lt)

           enddo
          enddo

           end if 

! Reset gfld%bmap with the combined version bmap_f
        if(trim(eps).eq.'href' .and. 
     +          associated(gfld%bmap)) gfld%bmap=bmap_f

               call packGB2_mean_derv(imean,isprd,derv_mn,
     +              derv_sp,nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +              iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)  !borrow gfld from what was got from previous direct variables
               write(*,*) 'pack _mean_derv done for', nv

               write(0,*) 'pack _prob_derv for ', nv
     +          ,jpd1,jpd2,jpd10,jpd27,iprob,jf,Lp,Lth,
     +          iens,iyr,imon,idy,ihr,ifhr,gribid

	write(0,*) 'min/max gfld%fld: ', minval(gfld%fld),maxval(gfld%fld)
	write(0,*) 'calling packGB2_prob_derv here'

               call packGB2_prob_derv(iprob,derv_pr,
     +              nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lth,
     +              iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld)  !borrow gfld from what was got from previous direct variables
               write(*,*) 'pack _prob_derv done for', nv
          end if
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Loop 2-3: deallocation

20021    continue

        deallocate (derv_mn)
        deallocate (derv_sp)
        deallocate (derv_pr)




2002   CONTINUE  ! end of derived variables loop
        deallocate (bmap_f)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c Loop 3: Process MXP variables (min,max, 10%,20%, ...90%) 
c 3-0: Allocate array mxp8 
c 3-1: Compute each mxp variable 
c 3-3: pack all mxp variables
c 3-4: deallocate
     
       write(*,*) 'Now begin to generate MXP products ......'

       if(nmxp.gt.0) then
         pmin=trim(eps)//'.pmin.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
         pmax=trim(eps)//'.pmax.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
         pmod=trim(eps)//'.pmod.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
         pp10=trim(eps)//'.pp10.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
         pp25=trim(eps)//'.pp25.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
         pp50=trim(eps)//'.pp50.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
         pp75=trim(eps)//'.pp75.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)
         pp90=trim(eps)//'.pp90.t'//cycle(ihr+1)//'z'//'.f'//trim(cfhr)

         call baopen(301,pmin,ierr)
         call baopen(302,pmax,ierr)
         call baopen(303,pmod,ierr)
         call baopen(304,pp10,ierr)
         call baopen(305,pp25,ierr)
         call baopen(306,pp50,ierr)
         call baopen(307,pp75,ierr)
         call baopen(308,pp90,ierr)

          if(trim(eps).eq.'sref' ) then
           if (jpd1.eq.1.and.
     +      (jpd2.eq.8.or.jpd2.eq.11.or.jpd2.eq.15)) then
              jpdtn=11
           else
              jpdtn=1
           end if
          else
            if (jpd1.eq.1.and.
     +      (jpd2.eq.8.or.jpd2.eq.11.or.jpd2.eq.15)) then
              jpdtn=8
           else
              jpdtn=0
           end if
          end if


       end if

       DO 2003 nv = 1, nmxp

c  Loop 3-0: allocation 

         Lq=max(1,qMlvl(nv))
         if (.NOT.allocated(mxp8)) then
          allocate (mxp8(jf,Lq,8))
         end if
c Loop 3-1: Computation


         call get_mxp(nv,ifunit,jpdtn,jf,iens,Lq,mxp8,weight) 
         write(*,*) 'call get_mxp done for ', nv

c Loop 3-2: Packing 

         gfld_temp=gfld

         jpd1=qk4(nv)
	write(0,*) 'jpd1 defined(e): ', jpd1
         jpd2=qk5(nv)
         jpd10=qk6(nv)
         jpd27=-9999


         do kq=1,8
          call packGB2_mxp(iqout(kq),kq,mxp8,    
     +      nv,jpd1,jpd2,jpd10,jpd27,jf,Lq,
     +      iens,iyr,imon,idy,ihr,ifhr,gribid,bmap_f,gfld_temp)           !gfld_temp  is used to send in other info 
         end do


         deallocate(mxp8)
2003  CONTINUE
       

       do irun=1,iens
        call baclose(ifunit(irun),ierr)
        call baclose(ipunit(irun),ierr)
       end do

       call baclose(imean,ierr)
       call baclose(ipmmn,ierr)
       call baclose(isprd,ierr)
       call baclose(iprob,ierr)
       call baclose(iavg,ierr)

       if (nmxp.gt.0) then
         call baclose(301,ierr)
         call baclose(302,ierr)
         call baclose(303,ierr)
         call baclose(304,ierr)
         call baclose(305,ierr)
         call baclose(306,ierr)
         call baclose(307,ierr)
         call baclose(308,ierr)
       end if
 
3000    continue   ! end of itime loop

      !call gf_free(gfld) 

      write(*,*) 'Entire program completed!'

      stop
      end



