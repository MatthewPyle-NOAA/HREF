c       program calibrate2_svr
c
c Just like calibrate2.f, except reads the output of the combined sevr probability table. DRB 12/10/2004
c       parameter (igemy=129, igemx=185) ! GEMPAK grid size
c

c  calibrate2_svr: to compute calibrated severe thunderstorm prob
c                   adapted from David Bright (SPC/NCEP)
c Just like calibrate_hrly.f, except reads the output of the combined probability table.
c   Author: Binbin Zhou/IMSG
c           June 9 2011
c   Input: data1,data2,jf
c   Output grd

       subroutine calibrate2_svr(data1,data2,jf,lyr,
     +                           grd)
       real data1(jf),data2(jf) 
       real prct(jf,2),prob,pa(200),pb(200),rk,bias(200),
     &       grd(jf),dum1,dum2
       integer itcnt,jstop,jdel
       character dummy*20,date*8,fhour*4,frun*2,dateetc*72
       character*2 lyr
       character*150 combine
c  
c fill all possible forecast values...	 
       itcnt = 0
       do i=0,100,10
        do j = 0,100,10
	 itcnt = itcnt + 1
         pa(itcnt) = float(i) ! fcst1
  	 pb(itcnt) = float(j) ! fcst2
        enddo
       enddo
       write(*,*) 'Number of possible combinations... ',itcnt

 1    format(a)
c Now read them in...
c
         prct(:,1)=data1(:)
         prct(:,2)=data2(:)


c
c Now the grids are in...        
c Do a little work here.  First, bin probs to 10% intervals.  Then, change verification
c to 1's or 0's.  
c
       do i=1,jf
        do k=1,2
         if(prct(i,k).lt.-1.0) prct(i,k) = -9999.0
cbzhou         if(prct(i,k).ge.-1.0.and.prct(i,k).lt.5.0)prct(i,k)= 0.0    
         if(prct(i,k).le.0.0)prct(i,k)= 0.0                                !Bzhou modify
         if(prct(i,k).gt.0.0.and.prct(i,k).lt.5.0)prct(i,k)= 5.0           !Bzhou modify
         if(prct(i,k).ge.5.0.and.prct(i,k).lt.15.) prct(i,k)=10.0
         if(prct(i,k).ge.15.0.and.prct(i,k).lt.25.)prct(i,k)=20.0
         if(prct(i,k).ge.25.0.and.prct(i,k).lt.35.)prct(i,k)=30.0
         if(prct(i,k).ge.35.0.and.prct(i,k).lt.45.)prct(i,k)=40.0
         if(prct(i,k).ge.45.0.and.prct(i,k).lt.55.)prct(i,k)=50.0
         if(prct(i,k).ge.55.0.and.prct(i,k).lt.65.)prct(i,k)=60.0
         if(prct(i,k).ge.65.0.and.prct(i,k).lt.75.)prct(i,k)=70.0
         if(prct(i,k).ge.75.0.and.prct(i,k).lt.85.)prct(i,k)=80.0
         if(prct(i,k).ge.85.0.and.prct(i,k).lt.95.)prct(i,k)=90.0
         if(prct(i,k).ge.95.0) prct(i,k)=100.0
        enddo
       enddo


c Read the file that contains the dependency coefficient and the bias... 
c
         combine='combine_probs_svr_layer'//lyr//'.out'
         open(unit=311,file=combine,status='old')
         do ii=1,itcnt
          read(311,*) pa(ii),pb(ii)
          read(311,*) dum1
          read(311,*) dum1,dum2
          read(311,*) dum1,dum2
          if((dum1+dum2).lt.0.001e-10) then
           bias(ii) = 0.0
          else
	   bias(ii) = (dum2/(dum1+dum2))*100.
          endif
 598      format(2f10.2)
 600      format(f10.2)
 601      format(8f10.2)
         enddo
	 close (311)

c Now...bin the forecasts into all psbl combinations.  Record if a hit occurred
c during this combined period too.
c

       grd = -9999.0

       do i=1,jf
           ii = 0
       loop1001: do k=1,itcnt
         if(nint(prct(i,1)).eq.nint(pa(k)).and.
     &     nint(prct(i,2)).eq.nint(pb(k))) then
            ii = k
            exit loop1001
         endif
       enddo loop1001


cbzhou       if(prct(i,1).lt.-9990.0.or.prct(i,2).lt.-9990.0.or.
       if(prct(i,1).le.0.0.or.prct(i,2).le.0.0.or.
     &      ii.eq.0) then
         grd(i) = -9999.0
        else
         grd(i) = bias(ii)
        endif

       enddo ! i loop
 
c       write(*,*)'In calibrate2_svr 4158 4159 for ',lyr,' after'
c       write(*,*) prct(4158,1),prct(4158,2), grd(4158)
c       write(*,*) prct(4159,1),prct(4159,2), grd(4159)


       write(*,*) 'Done calibrate2_svr'
       return
       end


