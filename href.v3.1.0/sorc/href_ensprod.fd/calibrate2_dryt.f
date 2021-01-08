c  calibrate2_dryt: to compute calibrated cloud phy thunder parameter(cptp) prob
c                   adapted from David Bright (SPC/NCEP)
c Just like calibrate_hrly.f, except reads the output of the combined probability table.
c   Author: Binbin Zhou/IMSG
c           June 2 2011
c   Input: drytp2,p03m01,jf
c   Output grd



       subroutine calibrate2_dryt(drytp2,p03m01,jf,grd)

c       parameter (igemy=129, igemx=185) ! GEMPAK grid size
c

       real drytp2(jf),p03m01(jf)
       real prct(jf,2),prob,pa(200),pb(200),rk,bias(200),
     &       grd(jf),dum1,dum2

       integer itcnt,jstop,jdel
       character dummy*20,date*8,fhour*4,frun*2,dateetc*72

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
c
        prct(:,1)=drytp2(:)
        prct(:,2)=100.0-p03m01(:)

c        write(*,*)'calibrate2_dryt before: drytp2  100-p03'
c        do i=1000,2000
c         write(*,*) i, prct(i,1),prct(i,2)
c        end do
c
c Now the grids are in...        
c Do a little work here.  First, bin probs to 10% intervals.  Then, change verification
c to 1's or 0's.  
       do i=1,jf
        do k=1,2
c         if(prct(i,k).lt.-1.0) prct(i,k) = -9999.0
         if(prct(i,k).le.0.0) prct(i,k) = -9999.0
         if(prct(i,k).gt.0.0.and.prct(i,k).lt.5.0)prct(i,k)= 0.0
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
c
c Read the file that contains the dependency coefficient and the bias... 
c
         open(unit=311,file='combine_pops_dryt.out',status='old')
         
         do ii=1,itcnt
          read(311,*) pa(ii),pb(ii)
          read(311,*) dum1
          read(311,*) dum1,dum2
          read(311,*) dum1,dum2
          if((dum1+dum2).lt.0.0001) then
           bias(ii) = 0.
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
       do i=1,jf
       ii = 0
       loop1001: do k=1,itcnt
         if(nint(prct(i,1)).eq.nint(pa(k)).and.
     &     nint(prct(i,2)).eq.nint(pb(k))) then
            ii = k
            exit loop1001
         endif
       enddo loop1001
 

        if(prct(i,1).lt.-9990.0.or.prct(i,2).lt.-9990.0.or.
     &     ii.eq.0) then
	 grd(i) = -9999.0
	else
 	 grd(i) = bias(ii)
	endif
	
       enddo ! i loop

c        write(*,*)'calibrate2_dryt after: drytp2  100-p03 grd'
c        do i=1000,2000
c         write(*,*) i, prct(i,1),prct(i,2),grd(i)
c        end do

              
       write(*,*) 'Done calibrate2_dryt'
       return
       end


