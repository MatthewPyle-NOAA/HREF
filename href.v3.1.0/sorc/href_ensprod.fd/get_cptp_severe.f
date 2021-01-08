cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutine get_severe_thunder_storm: diagmose severe thuner storm ensemble prob 
c             with David Bright method
c  This subroutine is for preparing before using calibrate2_sever subroutine
c  Author: Binbin Zhou, MArch 15, 2010
c    Modification: B. Zhou: June 16, 2011: adapted from SREF's version to NARRE's version 
c    Modification: B. Zhou: Aug 16, 2013: Rewrite for GRIB2 I/O   
c
c  Input: nv, itime,fhead,kx,ky,km,jf,iens,vrbl_pr 
c  output:  derv_pr (cptp products are saved in) 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   	subroutine get_cptp_severe(nv,cyc,fhr,ifunit,p03mp01,
     +      kx,ky,km,jpdtn,jf,iens,Lp,Lth,derv_pr,wgt,eps) 

         include 'parm.inc'

c    for derived variables
        Character*4 dvname(maxvar)
        Integer dk5(maxvar),dk6(maxvar),dk4(maxvar)
        Character*1 dMsignal(maxvar), dPsignal(maxvar)
        Integer dMlvl(maxvar), dMeanLevel(maxvar,maxmlvl)
        Integer dPlvl(maxvar), dProbLevel(maxvar,maxplvl)
        Character*1 dop(maxvar)
        Integer dTlvl(maxvar)
        Real    dThrs(maxvar,maxtlvl)
        Integer MPairLevel(maxvar,maxmlvl,2)
        Integer PPairLevel(maxvar,maxplvl,2)
        Character*19 files(50) 
        real wgt(30)
        integer ifunit(iens)
        character*2 cyc, fhr
        character*5 eps

        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop


        INTEGER, intent(IN) :: nv,kx,ky,km,jf,iens,Lp,Lth
        REAL,dimension(jf),intent(IN) :: p03mp01
        REAL,dimension(jf,Lp,Lth),intent(INOUT) :: derv_pr
        real,allocatable,dimension(:,:) :: terr_ens,hel1_ens,
     +    hel3_ens,cape_ens,dryt_ens,dapelcl_ens,cptp2_ens,
     +    dcape_ens,eshr_ens,tc500_ens

        !following are for severe storm
        real,allocatable,dimension(:) :: hicapep500   !means hicape>500 prob
        real,allocatable,dimension(:) :: hicapep1000   
        real,allocatable,dimension(:) :: hicapep2000   
        real,allocatable,dimension(:) :: hicapep3000   
        real,allocatable,dimension(:) :: hicapep250  
        real,allocatable,dimension(:) :: eshrp20
        real,allocatable,dimension(:) :: eshrp30
        real,allocatable,dimension(:) :: eshrp40
        real,allocatable,dimension(:) :: eshrp50
        real,allocatable,dimension(:) :: dcapep1000
        real,allocatable,dimension(:) :: dcapep2000
        real,allocatable,dimension(:) :: dapelclp10
        real,allocatable,dimension(:) :: tmpcpm25     !500mb temperature < -20C prob
        real,allocatable,dimension(:) :: tmpcpm20     !500mb temperature < -20C prob
        real,allocatable,dimension(:) :: tmpcpm15     !500mb temperature < -15C prob
        real,allocatable,dimension(:,:) :: layer21    !to store 21 combinations of cape/dcape etc in severe storm
        real,allocatable,dimension(:) :: svr_calib
        real svr_calib_max(21)

        !following are for ctpt 
        real,allocatable,dimension(:) :: cptpp1   !means cptp>1.0 prob       
        real,allocatable,dimension(:) :: drytp2   !means dryt>2.0 prob
        real calib_cptp(jf)

        integer id_dThrs,ndv


         allocate(terr_ens(jf,iens))
         allocate(hel1_ens(jf,iens))
         allocate(hel3_ens(jf,iens))
         allocate(cape_ens(jf,iens))
         allocate(dryt_ens(jf,iens))
         allocate(dapelcl_ens(jf,iens))
         allocate(cptp2_ens(jf,iens))
         allocate(dcape_ens(jf,iens))
         allocate(eshr_ens(jf,iens))
         allocate(tc500_ens(jf,iens))
         allocate(layer21(jf,21))
         allocate(svr_calib(jf))
         

         allocate(cptpp1(jf))
         allocate(drytp2(jf))

         allocate(hicapep500(jf))
         allocate(hicapep1000(jf))
         allocate(hicapep2000(jf))
         allocate(hicapep3000(jf))
         allocate(hicapep250(jf))
         allocate(eshrp20(jf))
         allocate(eshrp30(jf))
         allocate(eshrp40(jf))
         allocate(eshrp50(jf))
         allocate(dcapep1000(jf))
         allocate(dcapep2000(jf))
         allocate(dapelclp10(jf))
         allocate(tmpcpm25(jf))
         allocate(tmpcpm20(jf))
         allocate(tmpcpm15(jf))

        write(*,*) 'In get_cptp_severe ...'
        write(*,*) 'dThrs=',dThrs(nv,:) 
        write(*,*) nv,cyc,'  ',fhr,kx,ky,km,jf,iens,Lp,Lth
        write(*,*) ifunit 
        write(*,*) 'dk4,dk5,dk6=',dk4(nv),dk5(nv),dk6(nv)

        !step 1: get capes (hicape and dcape) of each members

          call get_hicape_dcape(ifunit,eps,
     +             kx,ky,km,jf,iens,
     +             terr_ens,hel1_ens,hel3_ens,cape_ens,
     +             dryt_ens, dapelcl_ens,cptp2_ens,
     +             dcape_ens,eshr_ens,tc500_ens)

        !step 2-1: get cptpp1 (ie ensemble prob of cptp>1)    !for now only hel3_ens and dryt_ens are are used used.      
        !note: in David code, this computation is done in getens_data.f (see line 8055), but here we move it to here 

         wsum=0.0
         do i=1,iens
           wsum=wsum+wgt(i)
         end do

         !write(*,*) (wgt(i),i=1,iens),wsum

         cptpp1=0.0
         do k=1,jf
          do i=1,iens
           if(hel3_ens(k,i).ge.1.0) cptpp1(k)=cptpp1(k)+wgt(i)   !from David code: getens_data.f, hel3(hel3_ens) or ltgpar2 is CPTP 
          end do

           cptpp1(k)=(cptpp1(k)/wsum)*100.0
         end do

         !step 2-2: get p03mp01 (ie ensemble prob of p03m>0.01 inch or 0.25 mm)         
         !          IN NARRE, since it is hourly putput, use 1-hr precipitation accumulation
c          do i=1,numvar
c           if(k5(i).eq.61) then
c            id_Thrs=0
c            do lh=1, Tlvl(i)                !search threshold 0.25 is there?
c             if (Thrs(i,lh).eq.0.1) then    !for NARRE, only 1hr apcp, so use 0.1 mm as threshold
c              id_Thrs=lh
c              goto 102
c             end if
c            end do
c102         if (id_Thrs.eq.0) then
c             write(*,*) 'p03m>0.01 prob not found'
c             stop 102
c            else
c             p03mp01(:)=vrbl_pr(:,i,1,id_Thrs)
c            end if
c           end if
c          end do



         !step 2-3: get dryt>2's prob for cptp product 
         !note: in David code, this computation is done in getens_data.f (see line 8055), but here we move it to here

         drytp2=0.0
         do k=1,jf
          do i=1,iens
           if(dryt_ens(k,i).ge.2.0) drytp2(k)=drytp2(k)+wgt(i)   !from David code: getens_data.f, hel3(hel3_ens) or ltgpar2 is CPTP
          end do
           drytp2(k)=(drytp2(k)/wsum)*100.0
         end do

         !step 3-1: get hicape>1000,2000,3000,500,250 prob

         hicapep1000=0.0
         hicapep2000=0.0
         hicapep3000=0.0
         hicapep500=0.0
         hicapep250=0.0
         do k=1,jf
          do i=1,iens
           if(cape_ens(k,i).ge.1000.)hicapep1000(k)=hicapep1000(k)
     +       +wgt(i)
           if(cape_ens(k,i).ge.2000.)hicapep2000(k)=hicapep2000(k)
     +       +wgt(i)
           if(cape_ens(k,i).ge.3000.)hicapep3000(k)=hicapep3000(k)
     +       +wgt(i)
           if(cape_ens(k,i).ge.500.)hicapep500(k)=hicapep500(k)
     +       +wgt(i)
           if(cape_ens(k,i).ge.250.)hicapep250(k)=hicapep250(k)
     +       +wgt(i)  
          end do
           hicapep1000(k)=(hicapep1000(k)/wsum)*100.0
           hicapep2000(k)=(hicapep2000(k)/wsum)*100.0
           hicapep3000(k)=(hicapep3000(k)/wsum)*100.0
           hicapep500(k)=(hicapep500(k)/wsum)*100.0
           hicapep250(k)=(hicapep250(k)/wsum)*100.0
         end do

         !step 3-2: get eshr > 20,30,40,50
          eshrp20=0.
          eshrp30=0.
          eshrp40=0.
          eshrp50=0.
         do k=1,jf
          do i=1,iens
           if(eshr_ens(k,i).ge.20.0) eshrp20(k)=eshrp20(k)+wgt(i)
           if(eshr_ens(k,i).ge.30.0) eshrp30(k)=eshrp30(k)+wgt(i)
           if(eshr_ens(k,i).ge.40.0) eshrp40(k)=eshrp40(k)+wgt(i)
           if(eshr_ens(k,i).ge.50.0) eshrp50(k)=eshrp50(k)+wgt(i)
          end do
           eshrp20(k)=(eshrp20(k)/wsum)*100.0
           eshrp30(k)=(eshrp30(k)/wsum)*100.0
           eshrp40(k)=(eshrp40(k)/wsum)*100.0
           eshrp50(k)=(eshrp50(k)/wsum)*100.0
         end do

         !step 3-3: get dcape > 1000, 2000 prob
          dcapep1000=0.
          dcapep2000=0.
         do k=1,jf
          do i=1,iens
           if(dcape_ens(k,i).ge.1000.)dcapep1000(k)=dcapep1000(k)
     +       +wgt(i)
           if(dcape_ens(k,i).ge.2000.)dcapep2000(k)=dcapep2000(k)
     +       +wgt(i)
          end do
          dcapep1000(k)=(dcapep1000(k)/wsum)*100.
          dcapep2000(k)=(dcapep2000(k)/wsum)*100.
         end do

        !step 3-4: get dapelcl > 1000 prob
           dapelclp10=0.
         do k=1,jf
          do i=1,iens
           if(dapelcl_ens(k,i).ge.1000.)dapelclp10(k)=dapelclp10(k)
     +       +wgt(i)
          end do
          dapelclp10(k)=(dapelclp10(k)/wsum)*100.
         end do

        !step 3-5: get t500 < -25C -20C, -15C  prob
         tmpcpm15=0.
         tmpcpm20=0.
         tmpcpm25=0.
         do k=1,jf
          do i=1,iens
           if(tc500_ens(k,i).ge.-900.0.and.
     +       tc500_ens(k,i).le.-15.0) tmpcpm15(k)=tmpcpm15(k)+wgt(i)
           if(tc500_ens(k,i).ge.-900.0.and.
     +       tc500_ens(k,i).le.-20.0) tmpcpm20(k)=tmpcpm20(k)+wgt(i)
           if(tc500_ens(k,i).ge.-900.0.and.
     +       tc500_ens(k,i).le.-25.0) tmpcpm25(k)=tmpcpm25(k)+wgt(i)
          end do
          tmpcpm15(k)=(tmpcpm15(k)/wsum)*100.
          tmpcpm20(k)=(tmpcpm20(k)/wsum)*100.
          tmpcpm25(k)=(tmpcpm25(k)/wsum)*100.
         end do

         !step 4: get calibrated cptp and severe ensemble prob in different cases:

         write(*,*) 'dPlvl(nv),dTlvl(nv)=',dPlvl(nv),dTlvl(nv)

         do lv=1,dPlvl(nv)
          do lh = 1, dTlvl(nv)

           if(dThrs(nv,lh).eq.1.0) then        !Lightning

             if(jf.eq.23865) then              !for grid#212
              call calibrate2_hrly_rgn3(cptpp1,p03mp01,kx,ky,jf,
     +          cyc,fhr,calib_cptp)
             else                              !Other grids than 212 
              call calibrate2(cptpp1,p03mp01,jf,calib_cptp)
             end if

             write(*,*) 'Lightning done'

           else if (dThrs(nv,lh).eq.2.0) then  !Dry lightning 
             call calibrate2_dryt(drytp2,p03mp01,jf,calib_cptp)
             write(*,*) 'Dry lightning done'

           else if (dThrs(nv,lh).eq.3.0) then    !severe thunder storm

             derv_pr(:,lv,lh) = -9999. 
             layer21 = -9999.

             call calibrate2_svr(hicapep500,eshrp30,
     +          jf,'01',svr_calib)
              layer21(:,1)=svr_calib(:)

             call calibrate2_svr(hicapep500,eshrp40,
     +          jf,'02',svr_calib)
              layer21(:,2)=svr_calib(:)

             call calibrate2_svr(hicapep1000,eshrp30,
     +          jf,'03',svr_calib)
              layer21(:,3)=svr_calib(:)

             call calibrate2_svr(hicapep1000,eshrp40,
     +          jf,'04',svr_calib)
              layer21(:,4)=svr_calib(:)

             call calibrate2_svr(hicapep2000,eshrp30,
     +          jf,'05',svr_calib)
              layer21(:,5)=svr_calib(:)

             call calibrate2_svr(hicapep2000,eshrp40,
     +          jf,'06',svr_calib)
              layer21(:,6)=svr_calib(:)

             call calibrate2_svr(hicapep3000,eshrp20,
     +          jf,'07',svr_calib)
              layer21(:,7)=svr_calib(:)

             call calibrate2_svr(hicapep3000,eshrp30,
     +          jf,'08',svr_calib)
              layer21(:,8)=svr_calib(:)

             call calibrate2_svr(hicapep3000,eshrp40,
     +          jf,'09',svr_calib)
              layer21(:,9)=svr_calib(:)

             call calibrate2_svr(hicapep250,eshrp30,
     +          jf,'10',svr_calib)
              layer21(:,10)=svr_calib(:)

             call calibrate2_svr(hicapep250,eshrp40,
     +          jf,'11',svr_calib)
              layer21(:,11)=svr_calib(:)

             call calibrate2_svr(hicapep250,eshrp50,
     +          jf,'12',svr_calib)
              layer21(:,12)=svr_calib(:)

             call calibrate2_svr(hicapep500,dcapep1000,
     +          jf,'13',svr_calib)
              layer21(:,13)=svr_calib(:)

             call calibrate2_svr(hicapep500,dcapep2000,
     +          jf,'14',svr_calib)
              layer21(:,14)=svr_calib(:)

             call calibrate2_svr(hicapep1000,dcapep1000,
     +          jf,'15',svr_calib)
              layer21(:,15)=svr_calib(:)

             call calibrate2_svr(hicapep1000,dcapep2000,
     +          jf,'16',svr_calib)
              layer21(:,16)=svr_calib(:)

             call calibrate2_svr(hicapep500,dapelclp10,
     +          jf,'17',svr_calib)
              layer21(:,17)=svr_calib(:)

             call calibrate2_svr(hicapep1000,dapelclp10,
     +          jf,'18',svr_calib)
              layer21(:,18)=svr_calib(:)

             call calibrate2_svr(hicapep500,tmpcpm15,
     +          jf,'19',svr_calib)
              layer21(:,19)=svr_calib(:)

             call calibrate2_svr(hicapep500,tmpcpm20,
     +          jf,'20',svr_calib)
              layer21(:,20)=svr_calib(:)

             call calibrate2_svr(hicapep500,tmpcpm25,
     +          jf,'21',svr_calib)
              layer21(:,21)=svr_calib(:)

1200         continue
             do i=1,jf
              svr_calib_max(:)=layer21(i,:)
              calib_cptp(i)=maxval(svr_calib_max)    !here calib_cptp is re-used array name
             end do
             write(*,*) 'Severe thunder done '

           end if 

           
           derv_pr(:,lv,lh)=calib_cptp(:) 
  
           do k=1,jf
            if (derv_pr(k,lv,lh).lt.0.0) derv_pr(k,lv,lh)=0.0
           end do


          end do
         end do

          deallocate(terr_ens)
          deallocate(hel1_ens)
          deallocate(hel3_ens)
          deallocate(cape_ens)
          deallocate(dryt_ens)
          deallocate(dapelcl_ens)
          deallocate(cptp2_ens)
          deallocate(dcape_ens)
          deallocate(eshr_ens)
          deallocate(tc500_ens)
          deallocate(cptpp1)
          deallocate(drytp2)

         deallocate(hicapep500)
         deallocate(hicapep1000)
         deallocate(hicapep2000)
         deallocate(hicapep3000)
         deallocate(hicapep250)
         deallocate(eshrp20)
         deallocate(eshrp30)
         deallocate(eshrp40)
         deallocate(eshrp50)
         deallocate(dcapep1000)
         deallocate(dcapep2000)
         deallocate(dapelclp10)
         deallocate(tmpcpm25)
         deallocate(tmpcpm20)
         deallocate(tmpcpm15)
         deallocate(layer21)
         deallocate(svr_calib)


          write(*,*) 'Done get_cptp_severe ...'

           return
           end
c
c   Subroutine get_hicape_dcape: To get 8 cape-related parameters of all members
c    to compute cptp prob (cloud phy thunder parameter) and hicape/dcape 
c    to compute severe storm prob developed by David Bright 
c    This code is cut from a part of David Bright code getens_data.f
c
c   Input: ifunit,kx,ky,km,jf,iens
c   Output: terr_ens,hel1_ens,hel3_ens,cape_ens,dryt_ens, dapelcl_ens,cptp2_ens
c
c   Author: Binbin Zhou  EMC/NCEP
c     March 15, 2010
c     August 16, 2013, B. Zhou: modification to GRIB2
c
c
c
	subroutine get_hicape_dcape(ifunit,eps,
     +             kx,ky,km,jf,iens,
     +             terr_ens,hel1_ens,hel3_ens,cape_ens,
     +             dryt_ens, dapelcl_ens,cptp2_ens,
     +             dcape_ens,eshr_ens,tc500_ens)


        use grib_mod
        type (gribfield) :: gfld

	CHARACTER	ht(37)*72, dummy*72
        real terr_ens(jf,iens),hel1_ens(jf,iens),
     +        hel3_ens(jf,iens),cape_ens(jf,iens),
     +        dryt_ens(jf,iens), dapelcl_ens(jf,iens),
     +        cptp2_ens(jf,iens),
     +        dcape_ens(jf,iens),eshr_ens(jf,iens),
     +        tc500_ens(jf,iens)
        integer  ifunit(iens)

	REAL		grid(jf), 
     &                  lclht(jf), 
     &                  ltgpar3,cptp2(jf)
        integer         
     &                  delpcpn,
     &                  imnstrt,iendmem,rht(37)
     &                  ,kx,ky
        real            
     &                  sfcp(jf),terr(jf),
     &                  sfcz(jf),grid3(jf),grid4(jf),
     &                  sfcu(jf),sfcv(jf),upru(37,jf),
     &                  sfct(jf),sfcsh(jf),sfctd(jf),
     &                  uprv(37,jf),
     &                  hel1(jf), hel3(jf),
     &                  cape(jf),
     &                  tmtm(37,jf),tdtd(37,jf),
     &                  rhrh(37,jf),
     &                  tmcape(100),mxcape(100),prcape(100),tpar,
     &                  dpar,ppar,hicape,hicin,peqlvl,
     &                  plfc,teqlvl,ltgcape,hicape2,hipres2,
     &                  hieql2,hicin2,ltgpar2,
     &                  maxcape(jf),
     &                  ucape(100),vcape(100),hishru,hishrv,
     &                  hishru2,hishrv2,effshru(jf),
     &                  effshrv(jf),
     &                  dbtpar,dbppar,dbdpar,
     &                  dbsfcwd,dbp(100),dbt(100),dbd(100),
     &                  dryt(jf),dape,dapelcl(jf),smoothd,
     &                  wgt(25)
		real 
     &                wndspd(jf),egnd,eair,qlflx,qsflx,dh,
     &  	      deltatm,snowdep(jf),consnow,tgnd(jf),
     &                consnow2,emiss,sigma,tempgnd,tcondcf,shwvwm,
     &                orgsndp,ustar,vstar,scape,sorig,dttraj,maxzlvl,
     &                mxzlvl2,maxz(jf),phalf,awcwndu,awcwndv,
     &                awcu(jf),awcv(jf),awccnt,hizlcl2,
     &                hizeql2
c
	INTEGER		level (2), luns (4), ienum, ilcnt, imem,
     &                  ips1, ips2, ilevel(50),imedian, i5start,
     &                  igdpcp, ip, ipcnt, imin, ncape, ptype,kvmax,
     &                  dbncape,inumout, header(2), iemtype,irstype
	CHARACTER	time(2)*20, garout*72, satfil*132,qfull*3,
     &                  pcptim*4, gdpcpf*72, pfld*4, qrange*1,qmode*4,
     &                  cnumin*3,custar*128, pfld1*4, pfld2*4
	LOGICAL		respnd, proces, drpflg, done,mintest,maxtest,
     &                  eoff

        CHARACTER  wrt*1  
        CHARACTER*19 files(50),fname   
        CHARACTER*5  eps

        real phail(100),zhail(100),thail(100),tdhail(100),
     +       dcape,dcapep,cape_down(jf)

         integer miss(jf) 

         !ilayers = KM
         ilayers = 29    !SREF's specific humidity only reach 300mb (1~29 levels)
         ifirst = 1 ! This is just a flag because members may be time lagged.
         ibeginc = 1 ! point to begin calculations at.
         iendinc = kx*ky ! point to end calculations at.
	 delpcpn = 3 ! hours between pcpn fields.


c the time lagged means should be last in the file...what is the starting imem number?
         imnstrt = 9999	 
         ht(1) = '1000'
         ht(2) = '975'
         ht(3) = '950'
         ht(4) = '925'
         ht(5) = '900'
         ht(6) = '875'
         ht(7) = '850'
         ht(8) = '825'
         ht(9) = '800'
         ht(10) = '775'
         ht(11) = '750'
         ht(12) = '725'
         ht(13) = '700'
         ht(14) = '675'
         ht(15) = '650'
         ht(16) = '625'
         ht(17) = '600'
         ht(18) = '575'
         ht(19) = '550'
         ht(20) = '525'
         ht(21) = '500'
         ht(22) = '475'
         ht(23) = '450'
         ht(24) = '425'
         ht(25) = '400'
         ht(26) = '375'
         ht(27) = '350'
         ht(28) = '325'
         ht(29) = '300'
         ht(30) = '275'
         ht(31) = '250'
         ht(32) = '225'
         ht(33) = '200'
         ht(34) = '175'
         ht(35) = '150'
         ht(36) = '125'
         ht(37) = '100'
c
         rht(1) = 1000
         rht(2) = 975
         rht(3) = 950
         rht(4) = 925
         rht(5) = 900
         rht(6) = 875
         rht(7) = 850
         rht(8) = 825
         rht(9) = 800
         rht(10) = 775
         rht(11) = 750
         rht(12) = 725
         rht(13) = 700
         rht(14) = 675
         rht(15) = 650
         rht(16) = 625
         rht(17) = 600
         rht(18) = 575
         rht(19) = 550
         rht(20) = 525
         rht(21) = 500
         rht(22) = 475
         rht(23) = 450
         rht(24) = 425
         rht(25) = 400
         rht(26) = 375
         rht(27) = 350
         rht(28) = 325
         rht(29) = 300
         rht(30) = 275
         rht(31) = 250
         rht(32) = 225
         rht(33) = 200
         rht(34) = 175
         rht(35) = 150
         rht(36) = 125
         rht(37) = 100	 
c

	 
         write(*,*) 'In get_hicape_dcape ....'

         write(*,*) kx,ky,km,jf,iens,eps
         write(*,*) ifunit
cc^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c if flag file not found, then highcape calculated all points.
             do k=ibeginc,iendinc
              lclht(k) = 1.0
             enddo
c Read in the flag file to run pts over US only!!!
             open(unit=72,file='hicape_model_flag_gem.out',
     &            status='old',err=1243)

 7775 read(72,1) dummy
 1     format(a)
       if(dummy(2:4).eq.'ROW') then
        backspace(unit=72)

        ikcnt = kx*ky
        do j=ky,1,-1
         ikcnt = ikcnt - kx
         if(j.ge.100) then
          read(72,*) dummy(1:8),(lclht(ikcnt+i),i=1,kx)
         else
          read(72,*) dummy(1:8),irow,(lclht(ikcnt+i),i=1,kx)
         endif
        enddo
       else
        goto 7775
       endif

       close(unit=72)
c put the flag values into lclht() array...
 1243  continue ! dumped here if no flag file
        write(*,*) 'Skip open hicape_model_flag_gem.out'
        write(*,*) 'lclht=',lclht(500)
cc^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       DO 1000 imem = 1, iens

         miss=0

         do i = 1,jf
           effshru(i) = -9999.0
           effshrv(i) = -9999.0
         enddo
          
         !!!!!!!!!!!!!!!  2D DATA  !!!!!!!!!!!!!!!!!
          !jpdtn=0
          jp27=-9999
          !sfc height
          call readGB2(ifunit(imem),jpdtn,3,5,1,0,jp27,gfld,
     +                 eps,ie) !Hsfc
           sfcz=gfld%fld
          
          !sfc pressure
           call readGB2(ifunit(imem),jpdtn,3,0,1,0,jp27,gfld,
     +                 eps,ie) !Sfc Pressure
           sfcp=gfld%fld/100.0

          !2m T (in C) 
           call readGB2(ifunit(imem),jpdtn,0,0,103,2,jp27,gfld,
     +                 eps,ie)
           sfct=gfld%fld-273.15

             do k=1,jf
              if (gfld%fld(k).eq.0.0) then  !NAMnest some T2m=0.0, missing
                miss(k)=1
              else 
                miss(k)=0
              end if               
             end do

          ! 2 meter dew pt from rsm is junk...need to calc 2 meter dewpt from
          ! specific humidity 
          !Note: In David B.'s code, he uses 'mul(spfh,pres@0%none)' (GEMpak)
          !But it seems that his result is little different from GRIB1 here
          !Thus cause different results for CIN  at some grid points

          !specific humidity (kg/kg) to derive vapor mixing ratio (g/kg)

          if (trim(eps).eq.'sref') then 
             call readGB2(ifunit(imem),jpdtn,1,0,103,2,jp27,gfld,
     +                     eps,ie)
              sfctd=gfld%fld/(1-gfld%fld/0.622)       !adjusted from specific humidity to mixing ratio
              sfctd=sfctd*1000.0                      ! --> g/kg
          else                                        !NARRE has no specific humidity
                                                      !so use RH to derive
             call readGB2(ifunit(imem),jpdtn,1,1,103,2,jp27,gfld,
     +                     eps,ie)  !RH (%)
             do k=ibeginc,iendinc

               tt=17.62*sfct(k)/(243.12 + sfct(k))  
               es = 6.112 *2.71828 ** tt            !(mb)
               RH=gfld%fld(k)/100.0
               q=0.622*es*RH/sfcp(k)                !(g/g)  specific humidity            
               sfctd(k)=1000.0*q/(1-q/0.622)        !(g/kg) vapor  mixing ratio

             end do
          end if

          ! 10 m u-wind
           call readGB2(ifunit(imem),jpdtn,2,2,103,10,jp27,gfld,
     +                    eps,ie)
           sfcu=gfld%fld

          ! 10 m v-wind
           call readGB2(ifunit(imem),jpdtn,2,3,103,10,jp27,gfld,
     +                    eps,ie)
           sfcv=gfld%fld
       
          !write(*,*) 'read 2D done' 

         !!!!!!!!!!!!!!!  3D DATA  !!!!!!!!!!!!!!!!!
           do i = 1,ilayers
              !Temperature profile
               call readGB2(ifunit(imem),jpdtn,0,0,100,rht(i),
     +           jp27,gfld,eps,ie)
               do k=ibeginc,iendinc
                tmtm(i,k)=gfld%fld(k)-273.15
               end do
 
               !write(*,*) '... T ... for ', rht(i)

            !Note: David B. uses Gempak 'mul(1000,mul(quo(relh,100),mixr(tmpc,pres)))'
            !to get specific humidity (ie RH*saturated_mixing_ratio or 1000*RH*saturated_mixing_ratio in g/kg)
            !mixr(tmpc,pres) is saturated mixing ratio. Then RH*saturated_mixing_ratio is mixing ratio
            !to represent aproximately specific humidity. 


             !Pure Specific humidity (kg/kg) profile
            if(trim(eps).eq.'sref') then
              call readGB2(ifunit(imem),jpdtn,1,0,100,rht(i),
     +              jp27,gfld,eps,ie) !Specific humidity
              do k=ibeginc,iendinc
               vr=gfld%fld(k)/(1.0-gfld%fld(k)/0.332)    !first adjusted from specific humidity to mixing ratio
               emb=(vr*rht(i))/(.622+vr)
               tdtd(i,k)=(237.3*log(emb/6.108))/(17.27-log(emb/6.108))   !Mixing ratio (g/kg)
              enddo
            else
              call readGB2(ifunit(imem),jpdtn,1,1,100,rht(i),
     +               jp27,gfld,eps,ie)  !RH(%) 
              do k=ibeginc,iendinc
                 tt=17.62*tmtm(i,k)/(243.12 + tmtm(i,k))
                 es = 6.112 *2.71828 ** tt            !(mb)
                 RH=gfld%fld(k)/100.0
                 q=0.622*es*RH/rht(i)                !(g/g)  specific humidity            
                 tdtd(i,k)=1000.0*q/(1-q/0.622)      !(g/kg) vapor  mixing ratio
              end do
            end if

            !write(*,*) '... Spedific H ... for ', rht(i)

             !U profile
              call readGB2(ifunit(imem),jpdtn,2,2,100,rht(i),
     +               jp27,gfld,eps,ie)
              do k=ibeginc,iendinc
               if(lclht(k).gt.0.5) then
                  if(rht(i).le.349) then
                   upru(i,k) = upru(27,k)
                  else
                   upru(i,k) = gfld%fld(k)
                  endif
               else
                  upru(i,k) = gfld%fld(k)
               end if
              end do

             !V profile
             call readGB2(ifunit(imem),jpdtn,2,3,100,rht(i),
     +                jp27,gfld,eps,ie)
              do k=ibeginc,iendinc
               if(lclht(k).gt.0.5) then
                  if(rht(i).le.349) then
                   uprv(i,k) = uprv(27,k)
                  else
                   uprv(i,k) = gfld%fld(k)
                  endif
               else
                  uprv(i,k) = gfld%fld(k)
               end if
              end do
 
              !write(*,*) '... U, V ... for ', rht(i)

           enddo  !end of ilayers

          do i = 1,ilayers
           do k=ibeginc,iendinc
             if(tdtd(i,k).gt.tmtm(i,k))tdtd(i,k)=tmtm(i,k)-.001
              if(i.eq.1.and.
     &           tdtd(1,k).gt.sfct(k))tdtd(1,k)=sfct(k)-.001
            end do
          end do
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         write(*,*) 'Done reading 2D & 3D data for member ',imem

c Binbin: please see here:
c Now I have everything to calculate cape.  Start at surface 
c and levels above...work your way up until 500 mb level.
c
c 8/10/03...start at 175 mb agl and end at 500 mb lvl.
c 9/28/03...start at surface and end at 501 mb abv sfc.
c
          
         do 2000 k=ibeginc,iendinc
           if(miss(k).eq.1) goto 2000 
                  if(lclht(k).gt.0.9) then
                   hicape = -9999.0
                   maxzlvl = -9999.0
                   mxzlvl2 = -9999.0
                   hizlcl2 = -9999.0
                   hizeql2 = -9999.0
                   hilfc2 = -9999.0
                   hicape2 = 0.0
                   hipres2 = -9999.0
                   hieql2 = -9999.0
                   hicin2 = 0.0
                   ltgpar2 = 0.0
		   ltgpar3 = 0.0
                   telcl = -9999.0
		   hishru2 = 0.
		   hishrv2 = 0.
                   dryt(k) = 0.
	           dapelcl(k) = 0.
                  else
                   hicape = -9999.0
                   maxzlvl = -9999.0
                   mxzlvl2 = -9999.0
                   hizlcl2 = -9999.0
                   hizeql2 = -9999.0
                   hilfc2 = -9999.0
                   hicape2 = -9999.0
                   hipres2 = -9999.0
                   hieql2 = -9999.0
                   hicin2 = -9999.0
                   ltgpar2 = -9999.0
                   ltgpar3 = -9999.0
                   telcl = -9999.0
		   hishru2 = -9999.0
		   hishrv2 = -9999.0
                   dryt(k) = -9999.0
	           dapelcl(k) = -9999.0
                   goto 9280
                  endif
c get sfc td from specific humidity
                   sfctd(k) = sfctd(k)/0.622 ! vapor pressure
                   sfctd(k) = (log(sfctd(k)) - 1.81)/
     &                  (19.8 - (log(sfctd(k))-1.81))
                  sfctd(k) = sfctd(k)*273. ! dew pt deg C
                  isfc = -9999
                  rpstop = sfcp(k) - 501.  ! need to speed things up a bit...
                   if(rpstop.lt.450.01) rpstop = 450.01
                   ipstop = ilayers - 1
           !if(imem.eq.5) write(*,*) 'sfctd rpstop =',sfctd(k),rpstop
                   do ir = 1,ilayers-1
                    if(rpstop.le.rht(ir).and.
     &                 rpstop.gt.rht(ir+1)) then
                        ipstop = ir
                        goto 6566
                    endif
                   enddo
 6566              continue
               !if (imem.eq.5) write(*,*) 'Before 6566'

               do i=1,ilayers-1
                   if(sfcp(k).gt.rht(i).and.i.eq.1) then
                  imixmax = 2
                    ncape = 1
                    isfc = 1
c call 1000 mb the sfc...it is the lowest model layer.
                    prcape(ncape) = sfcp(k) 
                    tmcape(ncape) = sfct(k)
                    mxcape(ncape) = sfctd(k)
                    ucape(ncape) = sfcu(k)
                    vcape(ncape) = sfcv(k)
                     do j=1,ilayers
                      ncape = ncape + 1
                      prcape(ncape) = rht(j)
                      tmcape(ncape) = tmtm(j,k)
                      mxcape(ncape) = tdtd(j,k)
                      ucape(ncape) = upru(j,k)
                      vcape(ncape) = uprv(j,k)
                     enddo
                     tpar = 0.
                     ppar = 0.
                     dpar = 0.
                     do imix=1,imixmax
                      tpar = tpar + tmcape(imix) 
                      dpar = dpar + mxcape(imix) 
                      ppar = ppar + prcape(imix) 
                     enddo
                     tpar = tpar/float(imixmax)
                     ppar = ppar/float(imixmax)
                     dpar = dpar/float(imixmax)
c need a full sounding for the downburst calculation.  Save it now.
                     dbncape = ncape
		     dbsfcwd = (ucape(1)**2+vcape(1)**2)**0.5
		     dbsfcwd = dbsfcwd*2.24 ! m/s to mph
		     dbtpar = tpar
		     dbppar = ppar
		     dbdpar =dpar
                     do imix=1,dbncape
		      dbt(imix) = tmcape(imix)
		      dbp(imix) = prcape(imix)
		      dbd(imix) = mxcape(imix)
		     enddo		    

        !if (k.eq.5000) write(*,*) 'CALL thermodynamics 1 ...'
       call thermodynamics(tmcape,mxcape,prcape,ncape,tpar,
     &                     dpar,ppar,hicape,hicin,peqlvl,
     &                     plcl,plfc,teqlvl,ltgcape,telcl,
     &                     maxzlvl)

        !if (k.eq.5000) write(*,*) 'CALL thermodynamics 1 done '
                     if(hicape.gt.hicape2.and.hicape.gt.0.0) then
                      hicape2 = hicape
                      mxzlvl2 = maxzlvl
                      hizlcl2 = plcl
                      hizeql2 = peqlvl
                      hipres2 = ppar
                      hieql2 = teqlvl
                      hilfc2 = plfc
                      hicin2 = hicin
                      ltgpar2=((ltgcape-100.)/100.)*(-19.0-teqlvl)
		      if(ltgcape.lt.100.0.or.teqlvl.gt.-20.0)ltgpar2=0.
                      if(ltgcape.ge.75.0.and.teqlvl.le.-17.5) then
		       ltgpar3 = 1.0
		      else
		       ltgpar3 = 0.0
		      endif
c calculate the shear from lpl to .5 peql level - the "effective" shear
                      wrt(1:1) = 'z'

                      call hishear(prcape,ucape,vcape,ppar,
     &                             peqlvl,hishru,hishrv,.5,wrt,tmcape)

                      hishru2 = hishru
		      hishrv2 = hishrv
                     endif          
                     goto 9124
                   elseif(sfcp(k).le.rht(i).and.
     &                    sfcp(k).ge.rht(i+1))then
                    imixmax = 2
                    ncape = 1
                    isfc = i+1
                    prcape(ncape) = sfcp(k) 
                    tmcape(ncape) = sfct(k)
                    mxcape(ncape) = sfctd(k)
                    ucape(ncape) = sfcu(k)
                    vcape(ncape) = sfcv(k)
                     do j=i+1,ilayers
                      ncape = ncape + 1
                      prcape(ncape) = rht(j)
                      tmcape(ncape) = tmtm(j,k)
                      mxcape(ncape) = tdtd(j,k)
                      ucape(ncape) = upru(j,k)
                      vcape(ncape) = uprv(j,k)
                     enddo
                     tpar = 0.
                     ppar = 0.
                     dpar = 0.
                     do imix=1,imixmax
                      tpar = tpar + tmcape(imix) 
                      dpar = dpar + mxcape(imix) 
                      ppar = ppar + prcape(imix) 
                     enddo
                     tpar = tpar/float(imixmax)
                     ppar = ppar/float(imixmax)
                     dpar = dpar/float(imixmax)
c need a full sounding for the downburst calculation.  Save it now.
                     dbncape = ncape
		     dbsfcwd = (ucape(1)**2+vcape(1)**2)**0.5
		     dbsfcwd = dbsfcwd*2.24 ! m/s to mph
		     dbtpar = tpar
		     dbppar = ppar
		     dbdpar =dpar
                     do imix=1,dbncape
		      dbt(imix) = tmcape(imix)
		      dbp(imix) = prcape(imix)
		      dbd(imix) = mxcape(imix)
		     enddo		  
   
       !if (k.eq.5000) write(*,*) 'CALL thermodynamics 2 ...'
       call thermodynamics(tmcape,mxcape,prcape,ncape,tpar,
     &                     dpar,ppar,hicape,hicin,peqlvl,
     &                     plcl,plfc,teqlvl,ltgcape,telcl,
     &                     maxzlvl)

       if (k.eq.5000) write(*,*) 'CALL thermodynamics 2 done '
                     if(hicape.gt.hicape2.and.hicape.gt.0.0) then
                      hicape2 = hicape
                      mxzlvl2 = maxzlvl
                      hizlcl2 = plcl
                      hizeql2 = peqlvl
                      hilfc2 = plfc
                      hipres2 = ppar
                      hieql2 = teqlvl
                      hicin2 = hicin
                      ltgpar2=((ltgcape-100.)/100.)*(-19.0-teqlvl)
		      if(ltgcape.lt.100.0.or.teqlvl.gt.-20.0)ltgpar2=0.
                      if(ltgcape.ge.75.0.and.teqlvl.le.-17.5) then
		       ltgpar3 = 1.0
		      else
		       ltgpar3 = 0.0
		      endif
c calculate the shear from lpl to .5 peql level - the "effective" shear
                      wrt(1:1) = 'z'
                      call hishear(prcape,ucape,vcape,ppar,
     &                             peqlvl,hishru,hishrv,.5,wrt,tmcape)
                      hishru2 = hishru
		      hishrv2 = hishrv
                     endif
                     goto 9124
                   endif ! this is the end of sfcp if block...
                  enddo ! enddo for the i (vertical) loop
 9124             continue
                  imixmax = 1 ! 2 is  ~ 50 mb mixed layer
                   do i=isfc,ipstop
                    ncape = 0
                     do j = i,ilayers
                      ncape = ncape + 1
                      prcape(ncape) = rht(j)
                      tmcape(ncape) = tmtm(j,k)
                      mxcape(ncape) = tdtd(j,k)
                      ucape(ncape) = upru(j,k)
                      vcape(ncape) = uprv(j,k)
                     enddo                  
                     tpar = 0.
                     ppar = 0.
                     dpar = 0.
                     do imix=1,imixmax
                      tpar = tpar + tmcape(imix) 
                      dpar = dpar + mxcape(imix) 
                      ppar = ppar + prcape(imix) 
                     enddo
                     tpar = tpar/float(imixmax)
                     ppar = ppar/float(imixmax)
                     dpar = dpar/float(imixmax)

       !if (k.eq.5000) write(*,*) 'CALL thermodynamics 3 ...'

       call thermodynamics(tmcape,mxcape,prcape,ncape,tpar,
     &                     dpar,ppar,hicape,hicin,peqlvl,
     &                     plcl,plfc,teqlvl,ltgcape,telcl,
     &                     maxzlvl)

       !if (k.eq.5000) write(*,*) 'CALL thermodynamics 3 done'

                     if(hicape.gt.hicape2.and.hicape.gt.0.0) then
                      hicape2 = hicape
                      mxzlvl2 = maxzlvl
                      hizlcl2 = plcl
                      hizeql2 = peqlvl
                      hilfc2 = plfc
                      hipres2 = ppar
                      hieql2 = teqlvl
                      hicin2 = hicin
c need a full sounding for the downburst calculation.  This was built above.  Thus,
c only need to update parcel info iff this is the mu vertical profile.
		      dbtpar = tpar
		      dbppar = ppar
		      dbdpar =dpar
                      ltgpar2=((ltgcape-100.)/100.)*(-19.0-teqlvl)
		      if(ltgcape.lt.100.0.or.teqlvl.gt.-20.0)ltgpar2=0.
                      if(ltgcape.ge.75.0.and.teqlvl.le.-17.5) then
		       ltgpar3 = 1.0
		      else
		       ltgpar3 = 0.0
		      endif
                      wrt(1:1) = 'z'


                      call hishear(prcape,ucape,vcape,ppar,
     &                             peqlvl,hishru,hishrv,.5,wrt,tmcape)

                      hishru2 = hishru
		      hishrv2 = hishrv
                     endif
                   enddo
 9280            continue ! Here if flag = 0 to skip this pt.
                 grid(k) = hicape2
                 if(nint(hipres2).ne.-9999) then
                  cape(k) = sfcp(k) - hipres2
                 else
                  cape(k) = hipres2
                 endif
c aviation additions, sept 20-21 2006, DRB.
c19dec2007                 if(mxzlvl2.gt.1000.0.and.hicape2.gt.25.0) then ! put some min threshold on the cloud...
                 if(mxzlvl2.gt.1250.0.and.hicape2.gt.25.0) then ! put some min threshold on the cloud...
                  maxz(k) = mxzlvl2*0.9 + sfcz(k) ! height msl of max parcel
                 else 
                  maxz(k) = -9999.0
                  awcu(k) = -9999.0
                  awcv(k) = -9999.0
                  goto 6711
                 endif
c determine mean wind in lower half of convective cloud for system speed...
                 phalf = hizlcl2 - ((hizlcl2-hizeql2)*0.5)
                 if(hizlcl2.lt.0.0.or.hizeql2.lt.0.0) phalf = 100000.
                 awcwndu = 0.
                 awcwndv = 0.
                 awccnt = 0.
                 awcu(k) = -9999.0
                 awcv(k) = -9999.0
                 do ia = 1,ilayers
                  if(rht(ia).le.hizlcl2.and.rht(ia).ge.phalf) then
                   awcwndu = awcwndu + upru(ia,k)
                   awcwndv = awcwndv + uprv(ia,k)
                   awccnt = awccnt + 1.0
                  endif
                  if(rht(ia).lt.phalf)goto 6710 ! too high...stop looking.
                 enddo
 6710            if(awccnt.gt.0.0) awcu(k) = (awcwndu/awccnt)*1.944 ! kts
                 if(awccnt.gt.0.0) awcv(k) = (awcwndv/awccnt)*1.944 ! kts
 6711            continue
c end of aviation additions.
                 terr(k) = hieql2
                 hel1(k) = hicin2
                 hel3(k) = ltgpar2
                 cptp2(k) = ltgpar3
		 maxcape(k) = grid(k)
		 effshru(k) = hishru2
		 effshrv(k) = hishrv2
                 
c
c Compute downburst and dry thunder potential for the most unstable vertical profile.
	    dbrain = -9999.0
	    dape = 0.
            if(grid(k).gt.10.0.and.lclht(k).gt.0.5) then ! require at least 10 j/kg for downburst calc...
	     call downburst(dbtpar,dbppar,dbdpar,dbt,dbp,dbd,
     &                       dbncape,dbsfcwd,dbrain,dape)

              dryt(k) =0.0
             if(hel3(k).ge.1.0.and.hel1(k).gt.-10.0) then ! cptp >= 1...hicin >= -10
c very little pcpn at sfc...this is the SPC definition of a dry TRW... <= 0.10	    
	      if(dbrain.ge.0.01.and.dbrain.lt.0.10) dryt(k)=1.
c no pcpn at sfc but lightning conditions met!!	      
	      if(dbrain.lt.0.01.and.dbrain.gt.-999.) dryt(k)=2.
	     endif
	     dapelcl(k) = dape
	    endif
c put values in corner points for nmap2...	    
	    if(k.eq.ibeginc) dryt(k) = 1.
	    if(k.eq.iendinc) dryt(k) = 2.
	    if(k.eq.ibeginc) dapelcl(k) = 1.
	    if(k.eq.iendinc) dapelcl(k) = 2.
c end experimental downburst and dry thunderstorm calculations...

c Compute dcape (this is for severe thunder storm 
               dcape = 0.0
               dcapep = -9999.0
              if(sfcp(k).gt.1000.0) then

c ground is below 1000 mb...build the sounding.
                ilvls = 1
                phail(ilvls) = sfcp(k)
                thail(ilvls) = sfct(k)
                tdhail(ilvls) = sfctd(k)
                do i=1,ilayers
                 ilvls = ilvls + 1
                 phail(ilvls) = rht(i)
                 thail(ilvls) = tmtm(i,k)
                 tdhail(ilvls) = tdtd(i,k)
                enddo

               else

                do i=1,ilayers-1
                 if (sfcp(k).le.rht(i).and.
     &               sfcp(k).gt.rht(i+1)) then
c found the ground...now build the vertical profile.
                  ilvls = 1
                  phail(ilvls) = sfcp(k)
                  thail(ilvls) = sfct(k)
                  tdhail(ilvls) = sfctd(k)
                   do ii=i+1,ilayers
                    ilvls = ilvls + 1
                    phail(ilvls) = rht(ii)
                    thail(ilvls) = tmtm(ii,k)
                    tdhail(ilvls) = tdtd(ii,k)
                   enddo
                  goto 8545
                 endif

                enddo
 8545           continue
               endif
               call capecalc_down(phail,thail,tdhail,ilvls,dcape,dcapep)
               cape_down(k) = dcape

2000     continue ! enddo for the k (grid point) loop
ccccccccccccccccc


        !final output are following 11 fields:
        do 2010 k=ibeginc,iendinc
         if(miss(k).eq.1) then
           terr_ens(k,imem)=-999.0
           hel1_ens(k,imem)=-999.0
           hel3_ens(k,imem)=-999.0
           cape_ens(k,imem)=-999.0
           dryt_ens(k,imem)=-999.0
           dapelcl_ens(k,imem)=-999.0
           cptp2_ens(k,imem)=-999.0
           eshr_ens(k,imem)=-999.0
           dcape_ens(k,imem)=-999.0
           tc500_ens(k,imem)=-999.0
         else
           terr_ens(k,imem)=terr(k)
           hel1_ens(k,imem)=hel1(k)
           hel3_ens(k,imem)=hel3(k)
           cape_ens(k,imem)=grid(k)
           dryt_ens(k,imem)=dryt(k)
           dapelcl_ens(k,imem)=dapelcl(k)
           cptp2_ens(k,imem)=cptp2(k)
           eshr_ens(k,imem)=sqrt(effshru(k)*effshru(k)+
     +      effshrv(k)*effshrv(k))
           dcape_ens(k,imem)=cape_down(k)        
           tc500_ens(k,imem)=tmtm(21,k)    !layer 21 is 500mb
         end if
2010    continue     

        


1000    CONTINUE  !end of member loop

        write(*,*) 'Done get_hicape ...' 
        return
	END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine thermodynamics(ttt,ddd,ppp,ilev,tpar,dpar,ppar,
     &                     cape,cin,eqlvl,plcl,lfc,teql,capem20,
     &                     telcl,maxzlvl)
c
c This code is modified homework code from Atmo/HWR 524.
c Designed to give thermodynamics properties of a parcel.
c DRB. 9/25/2003. TDM, SPC, Norman, OK, 73072.
c

c
c DAVID R. BRIGHT
c
c david.bright@noaa.gov (Phone: 670-5156)
c Homework for Atmo/Hydro 524
c Due: Oct. 2, 1998

c This program will compute basic thermodynamic properties
c for a parcel given a set of basic initial conditions.
c Equations for the solutions are derived from 524 class
c notes.  An iterative approach is used to solve for the
c LCL.

       real ttt(100),ddd(100),ppp(100),tpar,dpar,ppar,cape,
     &      cin,eqlvl, tvup(1200), pup(1200), tvenv(100),
     &      renv,eenv,envtv,tvpar,lfc,lcl,plcl,pparc,teql,
     &      tup(1200)

       real tc, pmb, qg, hm, tf, tk, ppa, pin, pkp, qkg, rkg,
     &      rg, tvc, tvk, rd, rho, r, ekpa, eskpa, emb, esmb,
     &      esmb2, rskg, rsg, tdc, tdk, tdf, tdc2, tdf2, tdk2,
     &      gt, eg, wg, test, rh, rho8, qskg7, qskg8a, qskg8b,
     &      lv, dqsdt, lrd, lrm, newp, newt, newe, newr, newrh,
     &      liftm, liftmh, tvave, newtd, newtv, newth, newthv,
     &      newq, newdq, newlv, newlrm, newpk, newek, thetae,
     &      newh, hinc, output, cond, newrs, newes, thetaw,
     &      capem20,telcl,cin2,cin3,maxplvl,zlvl,maxzlvl,
     &      tempor(1200),pup2(100),tvup2(100),tup2(100)


       integer upcnt,ilev

       upcnt = 1


       tc = tpar
       pmb = ppar
       pup(upcnt) = pmb
       td = dpar
c convert td to specific humidity (g/kg)
       qg = 1000.0*(.622/pmb)*
     &   exp((td*(19.8+log(6.1)) + 273.155*log(6.1))/(273.155 + td))

       tf = tc*1.8 +32.0
       tk = tc + 273.155
c
c PRESSURE
c ppa = pres in pascal
c pin = pres in inches
c pkp = pres in kpascal
c
       ppa = pmb*100.0
       pin = ppa/3386.0
       pkp = ppa/1000.0
c
c MOISTURE & PARTIAL PRESSURES
c qkg = spec. humidity kg/kg
c rg = mixing ratio g/kg
c rkg = mixing ratio kg/kg
c rskg = saturation mixing ratio kg/kg
c rsg  = saturation mixing ratio g/kg
c tvc = virtual temp (C)
c tvk = virtual temp (kelvin)
c ekpa   = vapor pressure (kpascals)
c eskpa  = saturation vapor pressure (kpascals)
c emb    = vapor pressure (mb)
c esmb   = sat. vapor pressure (mb)
c esmb2  = sat. vapor pressure (mb) from the formula I like
c tdc    = dew point temperature (C)
c tdk    = dew point temp (K)
c tdf    = dew point temp (F)
c tdc2   = dew point temp (C) from an iterative technique
c tdk2   = dew point temp (K) from "                    "
c tdf2   = dew point temp (F) from "                    "
c dqsdt  = rate of change of sat. specific hum. with temp.
c lv     = latent heat of vaporization at temperature tc
c
c
      qkg = qg/1000.0
      rkg = qkg/(1.0 - qkg)
      rg  = rkg*1000.0
      tvk = (273.155+tc)*(1.0 + 0.61*rkg)
      tvup(upcnt) = tvk
      tup(upcnt) = tc
      tvc = tvk - 273.155
      eskpa = 0.6108*(exp((17.27*tc)/(237.3 + tc))) ! Eqn from page 6 of class notes
      esmb  = eskpa*10.0
      emb   = ((rkg*ppa)/(rkg + 0.622))/100.0   ! From Wallace and Hobbs Eqn 2.61
      ekpa  = emb/10.0
      rh    = (emb/esmb)*100.
      rskg  = 0.622*(esmb/(pmb - esmb))
      rsg   = rskg*1000.0
      qskg7 = (qkg/rh)*100. ! Page 7 formula from notes
      qskg8a= (0.622*esmb)/(pmb - esmb + 0.621*esmb) ! Eqn 8a (exact)
      qskg8b= 0.622*(esmb/pmb)   ! Eqn 8b (approximate)
      lv    = 2.501
      dqsdt = (0.622/pkp)*((4098*eskpa)/((237.3 + tc)**2))
      esmb2=6.11*(EXP(9.081*(5.9529 - (752.61/(tc+273.155)) -
     &      (0.57*LOG(tc+273.155)))))
      tdc = (237.3*log(ekpa/0.6108))/(17.27 - log(ekpa/0.6108))
      tdf = tdc*1.8 + 32.
      tdk = tdc + 273.155
       rd = 287.0
       rho = ppa/(rd*tvk) ! eqn 4 ideal gas law
       rho8 = 3.486*pkp/(275.0 + tc) ! pg 8 hydro approx
       r = rd*(1. + 0.61*qkg)
      theta=tk*(1000./pmb)**0.286
      thetav=tvk*(1000./pmb)**0.286
      lrd = 9.81/1004.
      lrm = 9.81/(1004. + lv*1000000.*dqsdt)
      lcl = td - (.001296*td + .1963)*(tc - td)
      plcl = pmb*(((lcl+273.155)/tk)**3.4965)

      newp = plcl
       upcnt = upcnt + 1
       newt = lcl ! tlcl in degC
       newes = 0.6108*(exp((17.27*newt)/(237.3 + newt))) ! Eqn from page 6 of class notes
       newes = newes*10. ! convert to mb
       newrs = rkg
       newrh = 100.
       newh = (tc-newt)/lrd
       newe = newes
       newtd = lcl
       newth = theta
       newr = newrs
       newlv= 2.501
      thetae = newth*exp((2.501*1000000.0*newr)/
     &         (1004.*(newt+273.155)))

      if(thetae.le.telcl) then
       return
      else
       telcl = thetae
      endif

c
c Now...convert the temperature of the environment (degC) to
c a virtual temperature for later use in cape calculation.
c
       do i=1,ilev
        if(ddd(i).lt.-200..or.ttt(i).lt.-200..or.ppp(i).lt.-1.)then
         cape = -9999.
         cin = -9999.
         return
        endif
        eenv = 0.6108*(exp((17.27*ddd(i))/(237.3 + ddd(i))))
        eenv = eenv*10. ! convert to mb
        renv = 0.622*(eenv/(ppp(i) - eenv))
        tvenv(i) = (ttt(i)+273.155)*(1.0 + .61*renv)
       enddo

       pup(upcnt) = newp
       tvup(upcnt) = (newt+273.155)*(1.0 + 0.61*rkg)
       tup(upcnt) = newt
         liftm = (tc-newt)/lrd
         plcl = newp
         newh = liftm
         newtd = newt
         newtv = ((newt+273.155)*(1.0 + 0.61*newrs)) - 273.155
         tvlcl = newtv
         rlcl = newrs
         newth = theta
         newthv= (newtv+273.155)*((1000./newp)**0.286)
         newq = (0.622*newes)/(newp - newes + 0.622*newes)
         newpk = newp/10.
         newek = newe/10.
         newlv = 2.501
       rnewtk = newt + 273.155
       rlrm1= (9.81/1004.0)
       rlrm2= 1. + ((newlv*1000000.*newrs)/(287.*rnewtk))
       rlrm3= 1. + ((((newlv*1000000.)**2)*newrs)/
     &(1004.*461.*(rnewtk**2)))
       newlrm = rlrm1*rlrm2/rlrm3
        newdq = (((9.81/1004.)-newlrm)*1004.)/(newlrm*newlv*1000000.)
      newh = liftm
       rdmb = 15. ! make the increment a little bigger so sref runs faster.
 300  continue
      upcnt = upcnt + 1
       newp = newp - rdmb
       hinc = ((287.*(newtv+273.155))/9.81)*log((newp+rdmb)/newp)
       newt = newt  - (newlrm*hinc)


       newh = newh + hinc
       newe = 0.6108*(exp((17.27*newt)/(237.3 + newt))) ! Eqn from page 6 of class notes
       newe = newe*10. ! convert to mb
       newr = 0.622*(newe/(newp - newe))
       newrh = 100. ! This is a fact of saturated ascent.
       newtd = newt
       newtv = ((newt+273.155)*(1.0 + 0.61*newr)) - 273.155
       pup(upcnt) = newp
       tvup(upcnt) = newtv + 273.155
       tup(upcnt) = newt

       newth = (newt+273.155)*((1000./newp)**0.286) ! theta cons. so same as before
       newthv= (newtv+273.155)*((1000./newp)**0.286)
       newq = (0.622*newe)/(newp - newe + 0.622*newe)
       newpk = newp/10.
       newek = newe/10.
       newdq= (0.622/newpk)*((4098*newek)/((237.3 + newt)**2))
       newlv = 2.501
       rnewtk = newt + 273.155
       rlrm1= (9.81/1004.0)
       rlrm2= 1. + ((newlv*1000000.*newr)/(287.*rnewtk))
       rlrm3= 1. + ((((newlv*1000000.)**2)*newr)/
     &(1004.*461.*(rnewtk**2)))
       newlrm = rlrm1*rlrm2/rlrm3
        newdq = (((9.81/1004.)-newlrm)*1004.)/(newlrm*newlv*1000000.)
       dp = roldp-newp
      thetae = newth*exp((2.501*1000000.0*newr)/
     &         (1004.*(newt+273.155)))
 309        format(a14,6(f9.1))

       if(newp.gt.ppp(ilev)) goto 300
         ipcnt = 0
         thetav = tvup(1)*((1000./pup(1))**.286)
         theta = tup(1)*((1000./pup(1))**.286)
         do rp=pup(1),pup(2),-10.
          ipcnt = ipcnt + 1
          pup2(ipcnt) = rp
          tvup2(ipcnt) = thetav/((1000./rp)**.286)
          tup2(ipcnt) = theta/((1000./rp)**.286)
         enddo
       do i=1,upcnt
        tempor(i) = pup(i)
       enddo
       do i=1,ipcnt
        pup(i) = pup2(i)
       enddo
       do i=ipcnt+1,ipcnt+upcnt-1
        pup(i) = tempor(i-ipcnt+1)
       enddo
c virtual temp
       do i=1,upcnt
        tempor(i) = tvup(i)
       enddo
       do i=1,ipcnt
        tvup(i) = tvup2(i)
       enddo
       do i=ipcnt+1,ipcnt+upcnt-1
        tvup(i) = tempor(i-ipcnt+1)
       enddo
c temperture
       do i=1,upcnt
        tempor(i) = tup(i)
       enddo
       do i=1,ipcnt
        tup(i) = tup2(i)
       enddo
       do i=ipcnt+1,ipcnt+upcnt-1
        tup(i) = tempor(i-ipcnt+1)
       enddo
       upcnt = ipcnt+upcnt-1
c
c End dividing up the adiabatic layer into 10 mb increments...

       capem20=0.
       cape=0.
       cin =0.
       cin2=0.
       cin3=0.
       zlvl=0.
       maxplvl = -9999.0
       maxzlvl = -9999.0
       jold = 1
       eqlvl = -9999.
       lfc = -9999.
       do i=1,upcnt-1

c calculate the height agl of the parcel...
        zlvl = zlvl +
     &(14.63*(tvup(i)+tvup(i+1))*log(pup(i)/pup(i+1)))

        pparc = (pup(i) + pup(i+1))/2.0
        tvpar= (tvup(i)+tvup(i+1))/2.0
        if(pup(i).le.400.0.and.cape.le.0.001) return
        do j=jold,ilev-1
         if(ppp(j).ge.pparc.and.ppp(j+1).le.pparc) then
          envtv = tvenv(j) + (((tvenv(j+1)-tvenv(j))/
     &            (ppp(j+1)-ppp(j)))*(pparc-ppp(j)))
          jold = j ! save time...start next search here!
          if(tvpar.gt.envtv.and.pup(i).le.plcl) then
           if(lfc.lt.0.0) lfc = pup(i)
           cin = cin + cin2 ! elevated stable layer...add to cin total
           cin2 = 0. ! reset elevated cin layer
           maxplvl = -9999.0 ! reset the max parcel level
           maxzlvl = -9999.0 ! reset the max parcel level
           cape = cape+287.0*(tvpar-envtv)*log(pup(i)/pup(i+1))
           if(tvpar.le.273.155.and.tvpar.ge.253.155.and.lcl.gt.-10.0)
     &capem20=capem20+287.0*(tvpar-envtv)*log(pup(i)/pup(i+1))
            eqlvl = pup(i+1)
            teql = tup(i+1)
          else
           if(tvpar.le.envtv) then ! make sure not removing cin if superadiabatic
            if(lfc.lt.0.0)
     &cin = cin+287.0*(tvpar-envtv)*log(pup(i)/pup(i+1))
            if(lfc.gt.0.0)
     &cin2 = cin2+287.0*(tvpar-envtv)*log(pup(i)/pup(i+1))

            if(lfc.gt.0.0) cin3 = cin3 + cin2
            if(abs(cin3).ge.cape.and.maxplvl.lt.0.0) then
             maxplvl = pup(i+1)
             maxzlvl = zlvl ! height of maximum parcel AGL
            elseif(maxplvl.lt.0.0.and.i.eq.upcnt-1) then
             maxplvl = pup(i+1)
             maxzlvl = zlvl ! height of maximum parcel AGL
            endif

           endif
          endif
         endif
        enddo
       enddo
       if(cape.gt.5.0.and.teql.gt.-9900.0.and.
     &    nint(maxzlvl).eq.-9999) then
         maxzlvl = zlvl
         maxplvl = pup(upcnt)
       endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine hishear(prcape,ucape,vcape,pstart,pend,
     &                    uhishra,vhishra,pctup,wrt,tmcape)
c
c The option 'z' or 'p' is whether the weighted shear is wrt to height or pressure,
c respectively.
c
c Calculates the vector shear between pstart to 60% of pend
c
       real prcape(100),ucape(100),vcape(100),pstart,pend,
     &      uhishra,vhishra,peq2,du,dv,hishr2,pctup,tmcape(100),
     &      tmean,zmean,pmean
       integer istart, istop, istop2
       character wrt*1

       if(wrt(1:1).eq.'z'.or.wrt(1:1).eq.'Z') then
        wrt(1:1) = 'z'
       else
        wrt(1:1) = 'p'
       endif

       uhishra = -9999.0
       vhishra = -9999.0
       istart = 0
       istop = 0

       if(wrt(1:1).eq.'p') then
c calculate the vector shear from lfc thru 60% of layer to eql lvl.
       peq2 = pstart - ((pstart-pend)*pctup)

       if(prcape(1).lt.pstart) istart = 1
       do i=1,100-1
        if(prcape(i).ge.(pstart-1.0).and.
     &     prcape(i+1).le.(pstart+1.0)) istart=i
        if(prcape(i).ge.(peq2-1.0).and.
     &     prcape(i+1).le.(peq2+1.0)) istop=i+1
       enddo
c
       endif ! wrt = pressure
       if(wrt(1:1).eq.'z') then
        peq2 = pstart - ((pstart-pend)*1.0)
        do i=1,100-1
         if(prcape(i).ge.(pstart-1.0).and.
     &      prcape(i+1).le.(pstart+1.0)) istart=i
         if(prcape(i).ge.(peq2-1.0).and.
     &      prcape(i+1).le.(peq2+1.0)) istop=i+1
        enddo
        if(pctup.le.0.99) then
c estimate mean temp in this layer...temp is close enuf to virtual temp. for
c this purpose.
         tmean = 0.
         do i=istart,istop
          tmean = tmean + (tmcape(i)+273.155) ! degK
         enddo
         tmean = tmean/(float(istop-istart)+1.0)
         zmean = (287.*tmean/9.81)*log(prcape(istart)/prcape(istop))
         zmean = zmean*pctup ! this is the depth...from the lpl...to calc shear over
         tmean = 0.
         istop2 = istart+nint(float(istop-istart)*pctup)
         do i=istart,istop2
          tmean = tmean + (tmcape(i)+273.155) ! degK
         enddo
         tmean = tmean/(float(istop2-istart)+1.0)
         pmean = prcape(istart)/(exp(zmean*9.81/(287.0*tmean)))

c now find this pressure level in the sounding and adjust istop appropriatly.
        do i=1,100-1
         if(prcape(i).ge.(pmean-1.0).and.
     &     prcape(i+1).le.(pmean+1.0)) istop=i+1
        enddo
        endif
       endif ! wrt = height

       if(istart.eq.0.or.istop.eq.0) then
        write(*,*) 'Warning...could not find level to calc shear'
        write(*,*) 'Looking for pstart,peq2= ',pstart,peq2
        return
       endif

c
c Okay, now determine the vector shear (convert to kts) between
c istart and istop.
c
c GO TO A BULK SHEAR CALCULATION (12/17/2003). DRB.
c
c calculate the total vector shear...
c
c calculate the bulk shear...
c
                 du = 0.
                 dv = 0.
                 du = du + ucape(istop) - ucape(istart)
                 dv = dv + vcape(istop) - vcape(istart)
                 hishr2 = (du**2 + dv**2)**0.5
                 hishr2 = hishr2*1.944 ! m/s to kts
                 uhishra = du*1.944
                 vhishra = dv*1.944

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine downburst(tpar,ppar,dpar,ttt,ppp,ddd,ilev,
     &                    dbsfcwd,dbrain,dape)

c
c Calculates CAPE and a few other parameters.  Input is a GEMPAK
c snlist file.  Edit out the text and any mandatory levels below
c the surface.  DRB. 9/28/2003.
c
       real ttt(100),ddd(100),ppp(100),tpar,dpar,ppar,cape,
     &      cin,eqlvl,lcl,lfc,teql,ttt2(100),ddd2(100),ppp2(100),
     &      twmin,eenv,gamma,delta,wetbulb,de,der,ewet,wetold,
     &      wetbulb2,rwet,twmin2,dave,tave,pave,pcpn,tlcl,evap,
     &      evlvl,pcpeff,pcpnold,dbsfcwd,dbrain,dape
       integer ilev, numpar, ilev2, wcnt


       pcpeff = 1.0 ! precipitation efficency
       call thermody_up(ttt,ddd,ppp,ilev,tpar,dpar,ppar,
     &                     cape,cin,eqlvl,lcl,lfc,teql,pcpn,tlcl,
     &                     pcpeff)
      if(cape.le.0.0) then
        dbrain = -9999.0
        dape = 0.0
        return
      endif


       tpar = tlcl
       dpar = tlcl
       ppar = lcl
c now reverse the order of the input data "below" the parcel initiation pressure...

       ilev2 = 0
       do i=ilev,1,-1
        if(ppp(i).ge.ppar) then
         ilev2 = ilev2 + 1
         ppp2(ilev2) = ppp(i)
         ttt2(ilev2) = ttt(i)
         ddd2(ilev2) = ddd(i)
        endif
       enddo


       if(cape.gt.0.0) then
        pcpnold = pcpn
        call thermody_down(ttt2,ddd2,ppp2,ilev2,tpar,dpar,ppar,
     &                     cape,cin,eqlvl,lcl,lfc,teql,pcpn,evap,
     &                     evlvl)

        if(cape.lt.0.0) cape = 0.0
        dbsfcwd=dbsfcwd + (((2.*cape)**0.5)*2.24)
        dbrain = (pcpnold - evap)/25.4 ! rain at sfc in inches
        if(dbrain.lt.0.0) dbrain = 0.0
        dape = cape
       else
        dbrain = -9999.0
        dape = 0.0
       endif

       return
       end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine thermody_up(ttt,ddd,ppp,ilev,tpar,dpar,ppar,
     &                     cape,cin,eqlvl,plcl,lfc,teql,liq,lcl,
     &                     pcpeff)
c
c This code is modified homework code from Atmo/HWR 524.
c Designed to give thermodynamics properties of a parcel.
c DRB. 9/25/2003. TDM, SPC, Norman, OK, 73072.
c

c
c DAVID R. BRIGHT
c
c david.bright@noaa.gov (Phone: 670-5156)
c Homework for Atmo/Hydro 524
c Due: Oct. 2, 1998

c This program will compute basic thermodynamic properties
c for a parcel given a set of basic initial conditions.
c Equations for the solutions are derived from 524 class
c notes.  An iterative approach is used to solve for the
c LCL.

       real ttt(100),ddd(100),ppp(100),tpar,dpar,ppar,cape,
     &      cin,eqlvl, tvup(1200), pup(1200), tvenv(100),
     &      renv,eenv,envtv,tvpar,lfc,lcl,plcl,pparc,teql,
     &      tup(1200),qup(1200),qlost,liq,qenv(100),qenvi(1200),
     &      tempor(1200),pup2(100),tvup2(100),tup2(100)

       real tc, pmb, qg, hm, tf, tk, ppa, pin, pkp, qkg, rkg,
     &      rg, tvc, tvk, rd, rho, r, ekpa, eskpa, emb, esmb,
     &      esmb2, rskg, rsg, tdc, tdk, tdf, tdc2, tdf2, tdk2,
     &      gt, eg, wg, test, rh, rho8, qskg7, qskg8a, qskg8b,
     &      lv, dqsdt, lrd, lrm, newp, newt, newe, newr, newrh,
     &      liftm, liftmh, tvave, newtd, newtv, newth, newthv,
     &      newq, newdq, newlv, newlrm, newpk, newek, thetae,
     &      newh, hinc, output, cond, newrs, newes, thetaw,
     &      pcpeff

       integer upcnt,ilev,ipcnt
c
c first...convert the temperature of the environment (degC) to
c a virtual temperature for later use in cape calculation.
c
       do i=1,ilev
        if(ddd(i).lt.-200..or.ttt(i).lt.-200..or.ppp(i).lt.-1.)then
         cape = -9999.
         cin = -9999.
         return
        endif
        eenv = 0.6108*(exp((17.27*ddd(i))/(237.3 + ddd(i))))
        eenv = eenv*10. ! convert to mb
        renv = 0.622*(eenv/(ppp(i) - eenv))
        tvenv(i) = (ttt(i)+273.155)*(1.0 + .61*renv)
        qenv(i) = renv/(1.+renv)
       enddo

       upcnt = 1

c tc    = temperature deg cel.
c pmb   = pres mb
c qg    = specific hum. g/kg
c hm    = height of mountain in meters
c cond  = percent of pcpn condensed during ascent (100% for class)

       tc = tpar
       pmb = ppar
       pup(upcnt) = pmb
       td = dpar
c convert td to specific humidity (g/kg)
       qg = 1000.0*(.622/pmb)*
     &   exp((td*(19.8+log(6.1)) + 273.155*log(6.1))/(273.155 + td))

       tf = tc*1.8 +32.0
       tk = tc + 273.155
c
       ppa = pmb*100.0
       pin = ppa/3386.0
       pkp = ppa/1000.0
c
      qkg = qg/1000.0
      rkg = qkg/(1.0 - qkg)
      rg  = rkg*1000.0
      tvk = (273.155+tc)*(1.0 + 0.61*rkg)
      tvup(upcnt) = tvk
      tup(upcnt) = tc
      qup(upcnt) = qkg
      tvc = tvk - 273.155
      eskpa = 0.6108*(exp((17.27*tc)/(237.3 + tc))) ! Eqn from page 6 of class notes
      esmb  = eskpa*10.0
      emb   = ((rkg*ppa)/(rkg + 0.622))/100.0   ! From Wallace and Hobbs Eqn 2.61
      ekpa  = emb/10.0
      rh    = (emb/esmb)*100.
      rskg  = 0.622*(esmb/(pmb - esmb))
      rsg   = rskg*1000.0
      qskg7 = (qkg/rh)*100. ! Page 7 formula from notes
      qskg8a= (0.622*esmb)/(pmb - esmb + 0.621*esmb) ! Eqn 8a (exact)
      qskg8b= 0.622*(esmb/pmb)   ! Eqn 8b (approximate)
      lv    = 2.501
      dqsdt = (0.622/pkp)*((4098*eskpa)/((237.3 + tc)**2))
      esmb2=6.11*(EXP(9.081*(5.9529 - (752.61/(tc+273.155)) -
     &      (0.57*LOG(tc+273.155)))))
      tdc = (237.3*log(ekpa/0.6108))/(17.27 - log(ekpa/0.6108))
      tdf = tdc*1.8 + 32.
      tdk = tdc + 273.155
       rd = 287.0
       rho = ppa/(rd*tvk) ! eqn 4 ideal gas law
       rho8 = 3.486*pkp/(275.0 + tc) ! pg 8 hydro approx
       r = rd*(1. + 0.61*qkg)
      theta=tk*(1000./pmb)**0.286
      thetav=tvk*(1000./pmb)**0.286
      lrd = 9.81/1004.
      lrm = 9.81/(1004. + lv*1000000.*dqsdt)

      lcl = td - (.001296*td + .1963)*(tc - td)
      plcl = pmb*(((lcl+273.155)/tk)**3.4965)

      newp = plcl
       upcnt = upcnt + 1
       newt = lcl ! tlcl in degC
       newes = 0.6108*(exp((17.27*newt)/(237.3 + newt))) ! Eqn from page 6 of class notes
       newes = newes*10. ! convert to mb
       newrs = rkg
       newrh = 100.
       newh = (tc-newt)/lrd
       newe = newes
       newtd = lcl
       newth = theta
       newr = newrs
       newlv= 2.501
      thetae = newth*exp((2.501*1000000.0*newr)/
     &         (1004.*(newt+273.155)))
ccccccccccccc       thetaw = (newt+273.155)*((1000./newp)**.286)
       pup(upcnt) = newp
       tvup(upcnt) = (newt+273.155)*(1.0 + 0.61*rkg)
       tup(upcnt) = newt
       qup(upcnt) = newr/(1.0+newr)
         liftm = (tc-newt)/lrd
         plcl = newp
         newh = liftm
         newtd = newt
         newtv = ((newt+273.155)*(1.0 + 0.61*newrs)) - 273.155
         tvlcl = newtv
         rlcl = newrs
         newth = theta
         newthv= (newtv+273.155)*((1000./newp)**0.286)
         newq = (0.622*newes)/(newp - newes + 0.622*newes)
         newpk = newp/10.
         newek = newe/10.
         newlv = 2.501
c Use more exact eqns for sat. lapse rate and dq/dT than from text notes.
c Had too much error in text note version of dq/dT.  These are from
c Fleage and Businger text Atmo Physics, page 76.
       rnewtk = newt + 273.155
       rlrm1= (9.81/1004.0)
       rlrm2= 1. + ((newlv*1000000.*newrs)/(287.*rnewtk))
       rlrm3= 1. + ((((newlv*1000000.)**2)*newrs)/
     &(1004.*461.*(rnewtk**2)))
       newlrm = rlrm1*rlrm2/rlrm3
c now calc. newdq from this eqn.
        newdq = (((9.81/1004.)-newlrm)*1004.)/(newlrm*newlv*1000000.)

      newh = liftm
      rdmb = 10. ! This is the mb increment integral calc'd at.
 300  continue
      upcnt = upcnt + 1
       newp = newp - rdmb
       hinc = ((287.*(newtv+273.155))/9.81)*log((newp+rdmb)/newp)
       newt = newt  - (newlrm*hinc)

ccccccccc       newt = (thetaw/((1000./newp)**.286)) - 273.155

       newh = newh + hinc
c use tv at lcl to calc new p. based on hypsometric eqn and scale hgt.
cccc       newp = newp*exp((-9.81*hinc)/(287.*(tvlcl+273.15)))
c NOW NEED TO UPDATE EVERYTHING...LAPSE RATE, TV, ETC.
       newe = 0.6108*(exp((17.27*newt)/(237.3 + newt))) ! Eqn from page 6 of class notes
       newe = newe*10. ! convert to mb
       newr = 0.622*(newe/(newp - newe))
       newrh = 100. ! This is a fact of saturated ascent.
       newtd = newt
       newtv = ((newt+273.155)*(1.0 + 0.61*newr)) - 273.155
       pup(upcnt) = newp
       tvup(upcnt) = newtv + 273.155
       tup(upcnt) = newt
       qup(upcnt) = newr/(1.0 + newr)

       newth = (newt+273.155)*((1000./newp)**0.286) ! theta cons. so same as before
       newthv= (newtv+273.155)*((1000./newp)**0.286)
       newq = (0.622*newe)/(newp - newe + 0.622*newe)
       newpk = newp/10.
       newek = newe/10.
       newdq= (0.622/newpk)*((4098*newek)/((237.3 + newt)**2))
       newlv = 2.501
       rnewtk = newt + 273.155
       rlrm1= (9.81/1004.0)
       rlrm2= 1. + ((newlv*1000000.*newr)/(287.*rnewtk))
       rlrm3= 1. + ((((newlv*1000000.)**2)*newr)/
     &(1004.*461.*(rnewtk**2)))
       newlrm = rlrm1*rlrm2/rlrm3


c now calc. newdq from this eqn.
        newdq = (((9.81/1004.)-newlrm)*1004.)/(newlrm*newlv*1000000.)
       dp = roldp-newp
      thetae = newth*exp((2.501*1000000.0*newr)/
     &         (1004.*(newt+273.155)))
 309        format(a14,6(f9.1))

       if(newp.gt.ppp(ilev)) goto 300

c
c Okay...now have virtual temp of updraft in tvup() array and
c the pressure levels in pup() array.  Now...calculate the cape.
c
c Set up cin calculation below the lcl...
c
         ipcnt = 0
         thetav = tvup(1)*((1000./pup(1))**.286)
         theta = tup(1)*((1000./pup(1))**.286)
         do rp=pup(1),pup(2),-10.
          ipcnt = ipcnt + 1
          pup2(ipcnt) = rp
          tvup2(ipcnt) = thetav/((1000./rp)**.286)
          tup2(ipcnt) = theta/((1000./rp)**.286)
         enddo
c now merge the dry adiabatic data into main arrays...
c pressure
       do i=1,upcnt
        tempor(i) = pup(i)
       enddo
       do i=1,ipcnt
        pup(i) = pup2(i)
       enddo
       do i=ipcnt+1,ipcnt+upcnt-1
        pup(i) = tempor(i-ipcnt+1)
       enddo
c virtual temp
       do i=1,upcnt
        tempor(i) = tvup(i)
       enddo
       do i=1,ipcnt
        tvup(i) = tvup2(i)
       enddo
       do i=ipcnt+1,ipcnt+upcnt-1
        tvup(i) = tempor(i-ipcnt+1)
       enddo
c temperture
       do i=1,upcnt
        tempor(i) = tup(i)
       enddo
       do i=1,ipcnt
        tup(i) = tup2(i)
       enddo
       do i=ipcnt+1,ipcnt+upcnt-1
        tup(i) = tempor(i-ipcnt+1)
       enddo
c specific humidity
       do i=1,upcnt
        tempor(i) = qup(i)
       enddo
       do i=1,ipcnt
        qup(i) = tempor(1) ! q is conserved below the lcl
       enddo
       do i=ipcnt+1,ipcnt+upcnt-1
        qup(i) = tempor(i-ipcnt+1)
       enddo
       upcnt = ipcnt+upcnt-1
c
c End dividing up the adiabatic layer into 10 mb increments...
c

       cape=0.
       cin =0.
       cin2=0.
       jold = 1
       eqlvl = -9999.
       lfc = -9999.
       do i=1,upcnt-1
        pparc = (pup(i) + pup(i+1))/2.0
        tvpar= (tvup(i)+tvup(i+1))/2.0
c now interpolate environment to ascending parcel...
        do j=jold,ilev-1
         if(ppp(j).ge.pparc.and.ppp(j+1).le.pparc) then
          qenvi(i) = qenv(j) + (((qenv(j+1)-qenv(j))/
     &               (ppp(j+1)-ppp(j)))*(pparc-ppp(j)))
          envtv = tvenv(j) + (((tvenv(j+1)-tvenv(j))/
     &            (ppp(j+1)-ppp(j)))*(pparc-ppp(j)))
          jold = j ! save time...start next search here!
          if(tvpar.gt.envtv.and.pup(i).le.plcl) then
           if(lfc.lt.0.0) lfc = pup(i)
           cin = cin + cin2 ! elevated stable layer...add to cin total
           cin2 = 0. ! reset elevated cin layer
           cape = cape+287.0*(tvpar-envtv)*log(pup(i)/pup(i+1))
            eqlvl = pup(i+1)
            teql = tup(i+1)
          else
           if(lfc.lt.0.0.and.tvpar.lt.envtv)
     &cin = cin+287.0*(tvpar-envtv)*log(pup(i)/pup(i+1))
           if(lfc.gt.0.0.and.tvpar.lt.envtv)
     &cin2 = cin2+287.0*(tvpar-envtv)*log(pup(i)/pup(i+1))
          endif
         endif
        enddo
       enddo

c now...sum up the liquid condensed between lcl and el...It is just the pw in the updraft.
        liq = 0.
        do i = 1,upcnt-1
         qlost = (qup(i) + qup(i+1))*0.5 ! Ave q in updraft
         qlost = .1019*qlost*(pup(i)-pup(i+1))*100.
         qlost = qlost*pcpeff ! precipiation efficency
          if(pup(i).le.plcl.and.pup(i).ge.eqlvl.and.eqlvl.gt.0.0) then
c Start summing up the condensed water excess...this is the parcel q - the average of the layer.
           liq = liq + qlost
          endif
        enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine thermody_down(ttt,ddd,ppp,ilev,tpar,dpar,ppar,
     &                     cape,cin,eqlvl,plcl,lfc,teql,rain,liq,
     &                     evlvl)
c
c This code is modified homework code from Atmo/HWR 524.
c Designed to give thermodynamics properties of a parcel.
c DRB. 9/25/2003. TDM, SPC, Norman, OK, 73072.
c

c
c DAVID R. BRIGHT
c
c david.bright@noaa.gov (Phone: 670-5156)
c Homework for Atmo/Hydro 524
c Due: Oct. 2, 1998

c This program will compute basic thermodynamic properties
c for a parcel given a set of basic initial conditions.
c Equations for the solutions are derived from 524 class
c notes.  An iterative approach is used to solve for the
c LCL.

       real ttt(100),ddd(100),ppp(100),tpar,dpar,ppar,cape,
     &      cin,eqlvl, tvup(1200), pup(1200), tvenv(100),
     &      renv,eenv,envtv,tvpar,lfc,lcl,plcl,pparc,teql,
     &      tup(1200),rain,esat,rsat,qsenv(100),qsenvi,
     &      liq,qlost,qup(1200),evlvl,rainold,qenv(100),qenvi

       real tc, pmb, qg, hm, tf, tk, ppa, pin, pkp, qkg, rkg,
     &      rg, tvc, tvk, rd, rho, r, ekpa, eskpa, emb, esmb,
     &      esmb2, rskg, rsg, tdc, tdk, tdf, tdc2, tdf2, tdk2,
     &      gt, eg, wg, test, rh, rho8, qskg7, qskg8a, qskg8b,
     &      lv, dqsdt, lrd, lrm, newp, newt, newe, newr, newrh,
     &      liftm, liftmh, tvave, newtd, newtv, newth, newthv,
     &      newq, newdq, newlv, newlrm, newpk, newek, thetae,
     &      newh, hinc, output, cond, newrs, newes, thetaw

       integer upcnt,ilev
c
c first...convert the temperature of the environment (degC) to
c a virtual temperature for later use in cape calculation.
c
       do i=1,ilev
        if(ddd(i).lt.-200..or.ttt(i).lt.-200..or.ppp(i).lt.-1.)then
         cape = -9999.
         cin = -9999.
         return
        endif
        eenv = 0.6108*(exp((17.27*ddd(i))/(237.3 + ddd(i))))
        eenv = eenv*10. ! convert to mb
        renv = 0.622*(eenv/(ppp(i) - eenv))
        tvenv(i) = (ttt(i)+273.155)*(1.0 + .61*renv)
        esat = 0.6108*(exp((17.27*ttt(i))/(237.3 + ttt(i))))
        esat = esat*10. ! convert to mb
        rsat = 0.622*(esat/(ppp(i) - esat))
        qsenv(i) = rsat/(1.0+rsat)
        qenv(i) = renv/(1.0+renv)
       enddo
       upcnt = 1

c tc    = temperature deg cel.
c pmb   = pres mb
c qg    = specific hum. g/kg
c hm    = height of mountain in meters
c cond  = percent of pcpn condensed during ascent (100% for class)

       tc = tpar
       pmb = ppar
       pup(upcnt) = pmb
       td = dpar
       qg = 1000.0*(.622/pmb)*
     &   exp((td*(19.8+log(6.1)) + 273.155*log(6.1))/(273.155 + td))

       tf = tc*1.8 +32.0
       tk = tc + 273.155
c
c PRESSURE
c ppa = pres in pascal
c pin = pres in inches
c pkp = pres in kpascal
c
       ppa = pmb*100.0
       pin = ppa/3386.0
       pkp = ppa/1000.0
      qkg = qg/1000.0
      rkg = qkg/(1.0 - qkg)
      rg  = rkg*1000.0
      tvk = (273.155+tc)*(1.0 + 0.61*rkg)
      tvup(upcnt) = tvk
      tup(upcnt) = tc
      qup(upcnt) = qkg
      tvc = tvk - 273.155
      eskpa = 0.6108*(exp((17.27*tc)/(237.3 + tc))) ! Eqn from page 6 of class notes
      esmb  = eskpa*10.0
      emb   = ((rkg*ppa)/(rkg + 0.622))/100.0   ! From Wallace and Hobbs Eqn 2.61
      ekpa  = emb/10.0
      rh    = (emb/esmb)*100.
      rskg  = 0.622*(esmb/(pmb - esmb))
      rsg   = rskg*1000.0
      qskg7 = (qkg/rh)*100. ! Page 7 formula from notes
      qskg8a= (0.622*esmb)/(pmb - esmb + 0.621*esmb) ! Eqn 8a (exact)
      qskg8b= 0.622*(esmb/pmb)   ! Eqn 8b (approximate)
      lv    = 2.501
c to calc. dqsdt use the approx. qs = 0.622*es/p.  Then, dqs/dT = (.622/p)*des/dT
c if we assume the change in qs wrt T is done at constant total pressure. The
c formula for des/dT is given on page 6 of the class notes.
      dqsdt = (0.622/pkp)*((4098*eskpa)/((237.3 + tc)**2))
c For comparison of vapor pressure, here is a calculation from the formula
c that I normally use (from the text The Ceaseless Wind).
C CALCULATE SATURATION VAPOR PRES
      esmb2=6.11*(EXP(9.081*(5.9529 - (752.61/(tc+273.155)) -
     &      (0.57*LOG(tc+273.155)))))
c The eqn for dew pt below is based on solving eqn on pg 6 of class notes for T.
      tdc = (237.3*log(ekpa/0.6108))/(17.27 - log(ekpa/0.6108))
      tdf = tdc*1.8 + 32.
      tdk = tdc + 273.155
       rd = 287.0
       rho = ppa/(rd*tvk) ! eqn 4 ideal gas law
       rho8 = 3.486*pkp/(275.0 + tc) ! pg 8 hydro approx
       r = rd*(1. + 0.61*qkg)
c
c ADIABATIC PROCESSES
c theta = potential temp (K)
c thetav= virtual pot. temp (K)
c lrd   = dry adiabatic lapse rate (C/m)
c lrm   = saturated adiabatic lapse rate (C/m)
      theta=tk*(1000./pmb)**0.286
      thetav=tvk*(1000./pmb)**0.286
      lrd = 9.81/1004.
      lrm = 9.81/(1004. + lv*1000000.*dqsdt)
      lcl = td - (.001296*td + .1963)*(tc - td)
      plcl = pmb*(((lcl+273.155)/tk)**3.4965)

      newp = plcl
       upcnt = upcnt + 1
       newt = lcl ! tlcl in degC
       newes = 0.6108*(exp((17.27*newt)/(237.3 + newt))) ! Eqn from page 6 of class notes
       newes = newes*10. ! convert to mb
       newrs = rkg
       newrh = 100.
       newh = (tc-newt)/lrd
       newe = newes
       newtd = lcl
       newth = theta
       newr = newrs
       newlv= 2.501
      thetae = newth*exp((2.501*1000000.0*newr)/
     &         (1004.*(newt+273.155)))
       pup(upcnt) = newp
       tvup(upcnt) = (newt+273.155)*(1.0 + 0.61*rkg)
       tup(upcnt) = newt
       qup(upcnt) = newrs/(1.+newrs)
         liftm = (tc-newt)/lrd
         plcl = newp
         newh = liftm
         newtd = newt
         newtv = ((newt+273.155)*(1.0 + 0.61*newrs)) - 273.155
         tvlcl = newtv
         rlcl = newrs
         newth = theta
         newthv= (newtv+273.155)*((1000./newp)**0.286)
         newq = (0.622*newes)/(newp - newes + 0.622*newes)
         newpk = newp/10.
         newek = newe/10.
         newlv = 2.501
c Use more exact eqns for sat. lapse rate and dq/dT than from text notes.
c Had too much error in text note version of dq/dT.  These are from
c Fleage and Businger text Atmo Physics, page 76.
       rnewtk = newt + 273.155
       rlrm1= (9.81/1004.0)
       rlrm2= 1. + ((newlv*1000000.*newrs)/(287.*rnewtk))
       rlrm3= 1. + ((((newlv*1000000.)**2)*newrs)/
     &(1004.*461.*(rnewtk**2)))
       newlrm = rlrm1*rlrm2/rlrm3
c now calc. newdq from this eqn.
        newdq = (((9.81/1004.)-newlrm)*1004.)/(newlrm*newlv*1000000.)

      newh = liftm
      rdmb = -10. ! This is the mb increment integral calc'd at.
 300  continue
      upcnt = upcnt + 1
       newp = newp - rdmb
       hinc = ((287.*(newtv+273.155))/9.81)*log((newp+rdmb)/newp)
       newt = newt  - (newlrm*hinc)

ccccccccc       newt = (thetaw/((1000./newp)**.286)) - 273.155

       newh = newh + hinc
c use tv at lcl to calc new p. based on hypsometric eqn and scale hgt.
cccc       newp = newp*exp((-9.81*hinc)/(287.*(tvlcl+273.15)))
c NOW NEED TO UPDATE EVERYTHING...LAPSE RATE, TV, ETC.
       newe = 0.6108*(exp((17.27*newt)/(237.3 + newt))) ! Eqn from page 6 of class notes
       newe = newe*10. ! convert to mb
       newr = 0.622*(newe/(newp - newe))
       newrh = 100. ! This is a fact of saturated ascent.
       newtd = newt
       newtv = ((newt+273.155)*(1.0 + 0.61*newr)) - 273.155
       pup(upcnt) = newp
       tvup(upcnt) = newtv + 273.155
       tup(upcnt) = newt
       qup(upcnt) = newr/(1.+newr)
c        write(*,*) 'newt,newtv,newp= ',newt,newtv,newp

       newth = (newt+273.155)*((1000./newp)**0.286) ! theta cons. so same as before
       newthv= (newtv+273.155)*((1000./newp)**0.286)
       newq = (0.622*newe)/(newp - newe + 0.622*newe)
       newpk = newp/10.
       newek = newe/10.
       newdq= (0.622/newpk)*((4098*newek)/((237.3 + newt)**2))
       newlv = 2.501
       rnewtk = newt + 273.155
       rlrm1= (9.81/1004.0)
       rlrm2= 1. + ((newlv*1000000.*newr)/(287.*rnewtk))
       rlrm3= 1. + ((((newlv*1000000.)**2)*newr)/
     &(1004.*461.*(rnewtk**2)))
       newlrm = rlrm1*rlrm2/rlrm3


c now calc. newdq from this eqn.
        newdq = (((9.81/1004.)-newlrm)*1004.)/(newlrm*newlv*1000000.)
       dp = roldp-newp
      thetae = newth*exp((2.501*1000000.0*newr)/
     &         (1004.*(newt+273.155)))
 309        format(a14,6(f9.1))
       if(newp.lt.ppp(ilev)) goto 300

       rainold = rain
       evlvl = -9999.0
       liq = 0.
       cape=0.
       cin =0.
       cin2=0.
       jold = 1
       eqlvl = -9999.
       lfc = -9999.
       do i=1,upcnt-1
        pparc = (pup(i) + pup(i+1))/2.0
        tvpar= (tvup(i)+tvup(i+1))/2.0
c now interpolate environment to ascending parcel...
        do j=jold,ilev-1
          qsenvi = qsenv(j) + (((qsenv(j+1)-qsenv(j))/
     &             (ppp(j+1)-ppp(j)))*(pparc-ppp(j)))
          qenvi = qenv(j) + (((qenv(j+1)-qenv(j))/
     &             (ppp(j+1)-ppp(j)))*(pparc-ppp(j)))
         if(ppp(j).le.pparc.and.ppp(j+1).ge.pparc) then
          envtv = tvenv(j) + (((tvenv(j+1)-tvenv(j))/
     &            (ppp(j+1)-ppp(j)))*(pparc-ppp(j)))
          jold = j ! save time...start next search here!
           cin = cin + cin2 ! elevated stable layer...add to cin total
           cin2 = 0. ! reset elevated cin layer
           cape = cape+287.0*(tvpar-envtv)*log(pup(i)/pup(i+1))
           qlost = (qup(i) + qup(i+1))*0.5 ! Ave q in downdraft
           qlost = qlost-qenvi
           qlost = -0.1019*qlost*(pup(i)-pup(i+1))*100.
           liq = liq + qlost
           if(liq.ge.rain) then
c all the precip has evaporated...follow dry adiabatic descent now.
c            write(*,*) 'All precip has now evaporated!!!'
            evlvl = pparc
            rain = 9.99e9
            thetav = tvpar*((1000./pparc)**.286)
            do jj = i,upcnt
             tvup(jj) = thetav/((1000./pup(jj))**.286)
            enddo
           endif
         endif
        enddo
       enddo
       if(rain.gt.9.00e9) liq = rainold

      return
      end

       subroutine capecalc_down(ppp,ttt,ddd,ilev,cape,ppar)

c
c Calculates CAPE and a few other parameters.  Input is a GEMPAK
c snlist file.  Edit out the text and any mandatory levels below
c the surface.  DRB. 9/28/2003.
c
       real ttt(100),ddd(100),ppp(100),tpar,dpar,ppar,cape,
     &      cin,eqlvl,lcl,lfc,teql,ttt2(100),ddd2(100),ppp2(100),
     &      twmin,eenv,gamma,delta,wetbulb,de,der,ewet,wetold,
     &      wetbulb2,rwet,twmin2,dave,tave,pave
       integer ilev, numpar, ilev2, wcnt

c ttt,ddd,ppp are 1-d arrays containing temp, dewpt, and pres,
c respectively.  They are from the ground (=1) upward, and
c are in degC and mb.


c   find the level where wet bulb potential temp is min...
c   make sure it is within 500 mb of surface...
       twmin = 99999.9
       twmin2 = 99999.9
       do i=2,ilev-1
c
        if((ppp(i-1)-ppp(i+1)).le.75.0) then
        dave = ddd(i-1)*.333 + ddd(i)*.333 + ddd(i+1)*.333
        tave = ttt(i-1)*.333 + ttt(i)*.333 + ttt(i+1)*.333
        pave = ppp(i-1)*.333 + ppp(i)*.333 + ppp(i+1)*.333
        else
        dave = ddd(i)
        tave = ttt(i)
        pave = ppp(i)
        endif     
c
        eenv = 0.6108*(exp((17.27*dave)/(237.3 + dave)))
        eenv = eenv*10. ! convert to mb
        gamma = 6.6e-4*pave
        delta = (4098.0*eenv)/((dave+237.7)**2)
        wetbulb = ((gamma*tave)+(delta*dave))/(gamma+delta) !Tw degC
c
c Now iterate to precisely determine wet bulb temp.
        wcnt = 0
 800    continue
c calc vapor pressure at wet bulb temp
        ewet = 0.6108*(exp((17.27*wetbulb)/(237.3 + wetbulb)))
        ewet = ewet*10. ! convert to mb
        rwet = 0.622*(ewet/(pave - ewet))
        de = (0.0006355*pave*(tave-wetbulb))-(ewet-eenv)
        der= (ewet*(.0091379024 - (6106.396/(273.155+wetbulb)**2)))
     &       - (0.0006355*pave)
        wetold = wetbulb
        wetbulb = wetbulb - de/der
        wcnt = wcnt + 1
        if((abs(wetbulb-wetold)/wetbulb).gt..0001.and.(wcnt.lt.11))
     &   goto 800

        wetbulb2 = wetbulb
        wetbulb = (wetbulb+273.155)*((1000./pave)**0.286)
        wetbulb = wetbulb*exp((2.5e6*rwet)/(1004.*(273.155+wetbulb2)))
c        write(*,*) 'p,theta-w= ',ppp(i),wetbulb
c look for wet bulb min temperature at least 100 mb above ground...but no more than 500 mb
        if(wetbulb.lt.twmin.and.(ppp(1)-pave).lt.500.0.and.
     &      (ppp(1)-pave).gt.100.0) then
         twpre = pave
         twmin = wetbulb
         twmin2 = wetbulb2
        endif
       enddo

       tpar = twmin2
       dpar = twmin2
       ppar = twpre
c now reverse the order of the input data "below" the parcel initiation pressure...

       ilev2 = 0
       do i=ilev,1,-1
        if(ppp(i).ge.ppar) then
         ilev2 = ilev2 + 1
         ppp2(ilev2) = ppp(i)
         ttt2(ilev2) = ttt(i)
         ddd2(ilev2) = ddd(i)
        endif
       enddo

       call d_thermodynamics(ttt2,ddd2,ppp2,ilev2,tpar,dpar,ppar,
     &                     cape,cin,eqlvl,lcl,lfc,teql)


       return
       end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c DOWNDRAFT THERMODYNAMICS
       subroutine d_thermodynamics(ttt,ddd,ppp,ilev,tpar,dpar,ppar,
     &                     cape,cin,eqlvl,plcl,lfc,teql)
c
c This code is modified homework code from Atmo/HWR 524.
c Designed to give thermodynamics properties of a parcel.
c DRB. 9/25/2003. TDM, SPC, Norman, OK, 73072.
c

c
c DAVID R. BRIGHT
c
c david.bright@noaa.gov (Phone: 670-5156)
c Homework for Atmo/Hydro 524
c Due: Oct. 2, 1998

c This program will compute basic thermodynamic properties
c for a parcel given a set of basic initial conditions.
c Equations for the solutions are derived from 524 class
c notes.  An iterative approach is used to solve for the
c LCL.

       real ttt(100),ddd(100),ppp(100),tpar,dpar,ppar,cape,
     &      cin,eqlvl, tvup(1200), pup(1200), tvenv(100),
     &      renv,eenv,envtv,tvpar,lfc,lcl,plcl,pparc,teql,
     &      tup(1200)

       real tc, pmb, qg, hm, tf, tk, ppa, pin, pkp, qkg, rkg,
     &      rg, tvc, tvk, rd, rho, r, ekpa, eskpa, emb, esmb,
     &      esmb2, rskg, rsg, tdc, tdk, tdf, tdc2, tdf2, tdk2,
     &      gt, eg, wg, test, rh, rho8, qskg7, qskg8a, qskg8b,
     &      lv, dqsdt, lrd, lrm, newp, newt, newe, newr, newrh,
     &      liftm, liftmh, tvave, newtd, newtv, newth, newthv,
     &      newq, newdq, newlv, newlrm, newpk, newek, thetae,
     &      newh, hinc, output, cond, newrs, newes, thetaw

       integer upcnt,ilev
c
c first...convert the temperature of the environment (degC) to
c a virtual temperature for later use in cape calculation.
c

       do i=1,ilev
        if(ddd(i).lt.-200..or.ttt(i).lt.-200..or.ppp(i).lt.-1.)then
         cape = -9999.
         cin = -9999.
         return
        endif
        eenv = 0.6108*(exp((17.27*ddd(i))/(237.3 + ddd(i))))
        eenv = eenv*10. ! convert to mb
        renv = 0.622*(eenv/(ppp(i) - eenv))
        tvenv(i) = (ttt(i)+273.155)*(1.0 + .61*renv)
       enddo

       upcnt = 1

c tc    = temperature deg cel.
c pmb   = pres mb
c qg    = specific hum. g/kg
c hm    = height of mountain in meters
c cond  = percent of pcpn condensed during ascent (100% for class)

       tc = tpar
       pmb = ppar
       pup(upcnt) = pmb
       td = dpar
c convert td to specific humidity (g/kg)
       qg = 1000.0*(.622/pmb)*
     &   exp((td*(19.8+log(6.1)) + 273.155*log(6.1))/(273.155 + td))
c
c Calculate the basics...


c
c TEMPERATURE
c tf = temp deg f
c tk = temp kelvin

       tf = tc*1.8 +32.0
       tk = tc + 273.155
c
c PRESSURE
c ppa = pres in pascal
c pin = pres in inches
c pkp = pres in kpascal
c
       ppa = pmb*100.0
       pin = ppa/3386.0
       pkp = ppa/1000.0
c
c MOISTURE & PARTIAL PRESSURES
c qkg = spec. humidity kg/kg
c rg = mixing ratio g/kg
c rkg = mixing ratio kg/kg
c rskg = saturation mixing ratio kg/kg
c rsg  = saturation mixing ratio g/kg
c tvc = virtual temp (C)
c tvk = virtual temp (kelvin)
c ekpa   = vapor pressure (kpascals)
c eskpa  = saturation vapor pressure (kpascals)
c emb    = vapor pressure (mb)
c esmb   = sat. vapor pressure (mb)
c esmb2  = sat. vapor pressure (mb) from the formula I like
c tdc    = dew point temperature (C)
c tdk    = dew point temp (K)
c tdf    = dew point temp (F)
c tdc2   = dew point temp (C) from an iterative technique
c tdk2   = dew point temp (K) from "                    "
c tdf2   = dew point temp (F) from "                    "
c dqsdt  = rate of change of sat. specific hum. with temp.
c lv     = latent heat of vaporization at temperature tc
c
c
      qkg = qg/1000.0
      rkg = qkg/(1.0 - qkg)
      rg  = rkg*1000.0
      tvk = (273.155+tc)*(1.0 + 0.61*rkg)
      tvup(upcnt) = tvk
      tup(upcnt) = tc
      tvc = tvk - 273.155
      eskpa = 0.6108*(exp((17.27*tc)/(237.3 + tc))) ! Eqn from page 6 of class notes
      esmb  = eskpa*10.0
      emb   = ((rkg*ppa)/(rkg + 0.622))/100.0   ! From Wallace and Hobbs Eqn 2.61
      ekpa  = emb/10.0
      rh    = (emb/esmb)*100.
      rskg  = 0.622*(esmb/(pmb - esmb))
      rsg   = rskg*1000.0
      qskg7 = (qkg/rh)*100. ! Page 7 formula from notes
      qskg8a= (0.622*esmb)/(pmb - esmb + 0.621*esmb) ! Eqn 8a (exact)
      qskg8b= 0.622*(esmb/pmb)   ! Eqn 8b (approximate)
      lv    = 2.501
c to calc. dqsdt use the approx. qs = 0.622*es/p.  Then, dqs/dT = (.622/p)*des/dT
c if we assume the change in qs wrt T is done at constant total pressure. The
c formula for des/dT is given on page 6 of the class notes.
      dqsdt = (0.622/pkp)*((4098*eskpa)/((237.3 + tc)**2))
c For comparison of vapor pressure, here is a calculation from the formula
c that I normally use (from the text The Ceaseless Wind).
C CALCULATE SATURATION VAPOR PRES
      esmb2=6.11*(EXP(9.081*(5.9529 - (752.61/(tc+273.155)) -
     &      (0.57*LOG(tc+273.155)))))
c The eqn for dew pt below is based on solving eqn on pg 6 of class notes for T.
      tdc = (237.3*log(ekpa/0.6108))/(17.27 - log(ekpa/0.6108))
      tdf = tdc*1.8 + 32.
      tdk = tdc + 273.155
c
c BASIC STATE
c rho = density (kg/m**3) based on eqn 4 (ideal gas law)
c rho8= density (kg/m**3) based on page 8 approx.
c r = gas constant
c [rd = gas constant for dry air = 287 J kg-1 K-1]
c
       rd = 287.0
       rho = ppa/(rd*tvk) ! eqn 4 ideal gas law
       rho8 = 3.486*pkp/(275.0 + tc) ! pg 8 hydro approx
       r = rd*(1. + 0.61*qkg)
c
c ADIABATIC PROCESSES
c theta = potential temp (K)
c thetav= virtual pot. temp (K)
c lrd   = dry adiabatic lapse rate (C/m)
c lrm   = saturated adiabatic lapse rate (C/m)

      theta=tk*(1000./pmb)**0.286
      thetav=tvk*(1000./pmb)**0.286
      lrd = 9.81/1004.
      lrm = 9.81/(1004. + lv*1000000.*dqsdt)
c
c
c LCL Calculations...
c
c
c First, find the LCL.  Here is the method I will use:
c Pot. temp is conserved.  I will use an iterative process
c and calc. a new temp. every .1 mb based on this fact.
c At that new pressure and temperature calc. a new sat.
c mixing ratio.  Since mixing ratio is conserved with height,
c when the sat. mixing ratio = initial mixing ratio rh is
c 100 % and LCL is found.  Then use lapse rate to back out
c height in meters...can compare with hypsometric eqn too.
c
c newp = the pressure at lcl (mb)
c newt = the temp at lcl (C)
c newe = the vapor press = sat vapor press at lcl
c newr = the mixing ratio = sat mixing ratio at lcl
c liftm = the meters from original p to lcl based on lapse rate
c liftmh = the meters from original p to lcl based on hypsometrc eqn
c tvave = the average virt. pot. temp between base and lcl.
c newtd = dew pt at lcl
c newtv = virt temp at lcl
c newth = theta at lcl
c newthv= virt pot. temp at lcl
c newq = specific humidity
c newdq = rate of change of specific hum with temp.
c newrh = rh
c newlv = latent heat.
c newlrm = sat. ad lapse rate
c

      lcl = td - (.001296*td + .1963)*(tc - td)
      plcl = pmb*(((lcl+273.155)/tk)**3.4965)

      newp = plcl
       upcnt = upcnt + 1
       newt = lcl ! tlcl in degC
       newes = 0.6108*(exp((17.27*newt)/(237.3 + newt))) ! Eqn from page 6 of class notes
       newes = newes*10. ! convert to mb
       newrs = rkg
       newrh = 100.
       newh = (tc-newt)/lrd
       newe = newes
       newtd = lcl
       newth = theta
       newr = newrs
       newlv= 2.501
      thetae = newth*exp((2.501*1000000.0*newr)/
     &         (1004.*(newt+273.155)))
ccccccccccccc       thetaw = (newt+273.155)*((1000./newp)**.286)
       pup(upcnt) = newp
       tvup(upcnt) = (newt+273.155)*(1.0 + 0.61*rkg)
       tup(upcnt) = newt
         liftm = (tc-newt)/lrd
         plcl = newp
         newh = liftm
         newtd = newt
         newtv = ((newt+273.155)*(1.0 + 0.61*newrs)) - 273.155
         tvlcl = newtv
         rlcl = newrs
c         newth = (newt+273.155)*((1000./newp)**0.286) ! theta cons. so same as before
         newth = theta
         newthv= (newtv+273.155)*((1000./newp)**0.286)
         newq = (0.622*newes)/(newp - newes + 0.622*newes)
         newpk = newp/10.
         newek = newe/10.
         newlv = 2.501
c Use more exact eqns for sat. lapse rate and dq/dT than from text notes.
c Had too much error in text note version of dq/dT.  These are from
c Fleage and Businger text Atmo Physics, page 76.
       rnewtk = newt + 273.155
       rlrm1= (9.81/1004.0)
       rlrm2= 1. + ((newlv*1000000.*newrs)/(287.*rnewtk))
       rlrm3= 1. + ((((newlv*1000000.)**2)*newrs)/
     &(1004.*461.*(rnewtk**2)))
       newlrm = rlrm1*rlrm2/rlrm3
        newdq = (((9.81/1004.)-newlrm)*1004.)/(newlrm*newlv*1000000.)

      newh = liftm
      rdmb = -10. ! This is the mb increment integral calc'd at.
 300  continue
      upcnt = upcnt + 1
       newp = newp - rdmb
       hinc = ((287.*(newtv+273.155))/9.81)*log((newp+rdmb)/newp)
       newt = newt  - (newlrm*hinc)
ccccccccc       newt = (thetaw/((1000./newp)**.286)) - 273.155

       newh = newh + hinc
c use tv at lcl to calc new p. based on hypsometric eqn and scale hgt.
cccc       newp = newp*exp((-9.81*hinc)/(287.*(tvlcl+273.15)))
c NOW NEED TO UPDATE EVERYTHING...LAPSE RATE, TV, ETC.
       newe = 0.6108*(exp((17.27*newt)/(237.3 + newt))) ! Eqn from page 6 of class notes
       newe = newe*10. ! convert to mb
       newr = 0.622*(newe/(newp - newe))
       newrh = 100. ! This is a fact of saturated ascent.
       newtd = newt
       newtv = ((newt+273.155)*(1.0 + 0.61*newr)) - 273.155
       pup(upcnt) = newp
       tvup(upcnt) = newtv + 273.155
       tup(upcnt) = newt

       newth = (newt+273.155)*((1000./newp)**0.286) ! theta cons. so same as before
       newthv= (newtv+273.155)*((1000./newp)**0.286)
       newq = (0.622*newe)/(newp - newe + 0.622*newe)
       newpk = newp/10.
       newek = newe/10.
       newdq= (0.622/newpk)*((4098*newek)/((237.3 + newt)**2))
       newlv = 2.501
       rnewtk = newt + 273.155
       rlrm1= (9.81/1004.0)
       rlrm2= 1. + ((newlv*1000000.*newr)/(287.*rnewtk))
       rlrm3= 1. + ((((newlv*1000000.)**2)*newr)/
     &(1004.*461.*(rnewtk**2)))
       newlrm = rlrm1*rlrm2/rlrm3


c now calc. newdq from this eqn.
        newdq = (((9.81/1004.)-newlrm)*1004.)/(newlrm*newlv*1000000.)
       dp = roldp-newp
      thetae = newth*exp((2.501*1000000.0*newr)/
     &         (1004.*(newt+273.155)))
 309        format(a14,6(f9.1))

       if(newp.lt.ppp(ilev)) goto 300

c
c Okay...now have virtual temp of updraft in tvup() array and
c the pressure levels in pup() array.  Now...calculate the cape.
c
c       do i=1,upcnt
c        write(*,*) pup(i),tvup(i)
c       enddo
       cape=0.
       cin =0.
       cin2=0.
       jold = 1
       eqlvl = -9999.
       lfc = -9999.
       do i=1,upcnt-1
        pparc = (pup(i) + pup(i+1))/2.0
        tvpar= (tvup(i)+tvup(i+1))/2.0
c now interpolate environment to ascending parcel...
        do j=jold,ilev-1
         if(ppp(j).le.pparc.and.ppp(j+1).ge.pparc) then
          envtv = tvenv(j) + (((tvenv(j+1)-tvenv(j))/
     &            (ppp(j+1)-ppp(j)))*(pparc-ppp(j)))
          jold = j ! save time...start next search here!
           cin = cin + cin2 ! elevated stable layer...add to cin total
           cin2 = 0. ! reset elevated cin layer
           cape = cape+287.0*(tvpar-envtv)*log(pup(i)/pup(i+1))
         endif
        enddo
       enddo
      return
      end


