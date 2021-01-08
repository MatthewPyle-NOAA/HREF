c  calibrate2: to compute calibrated cloud phy thunder parameter(cptp) prob 
c              adapted from David Bright (SPC/NCEP)
c   Author: Binbin Zhou/SAIC
c           Apr 1 2010
c   Input: cptpp1,p03m01,jf
c   Output grd  
c
       subroutine calibrate2(cptpp1,p03m01,jf, grd)
c
c Just like calibrate.f, except reads the output of the combined probability table.

c       parameter (igemy=129, igemx=185) ! GEMPAK grid size
c
       real cptpp1(jf),p03m01(jf)
       real prct(jf,2),prob,pa(200),pb(200),rk,bias(200),
     &       grd(jf),dum1,dum2

       integer itcnt,jstop,jdel
       character*150 pctfila,pctfilb
       character dummy*20,date*8,fhour*4,frun*2,dateetc*72

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

       write(*,*) 'calibrate2--->'
       write(*,*) 'cptpp1(1562),p03m01(1562),jf=',
     +   cptpp1(1562),p03m01(1562),jf
       write(*,*) 'Number of possible combinations... ',itcnt

        prct(:,1)=cptpp1(:)
        prct(:,2)=p03m01(:)

c
c Now the grids are in...        
c Do a little work here.  First, bin probs to 10% intervals.  Then, change verification
c to 1's or 0's.  
       do i=1,jf
        do k=1,2
         if(prct(i,k).lt.-1.0) prct(i,k) = -9999.0
         if(prct(i,k).ge.-1.0.and.prct(i,k).lt.5.0)prct(i,k)= 0.0
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
         open(unit=301,file='combine_pops.out',status='old')
         do ii=1,itcnt
          read(301,*) pa(ii),pb(ii)
          read(301,*) dum1
          read(301,*) dum1,dum2
          read(301,*) dum1,dum2
	  bias(ii) = (dum2/(dum1+dum2))*100.
 598      format(2f10.2)
 600      format(f10.2)
 601      format(8f10.2)
         enddo
	 close (301)

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
     &      ii.eq.0) then
	 grd(i) = -9999.0
	else
 	 grd(i) = bias(ii)
	endif
	
       enddo ! i loop
       
       write(*,*) 'adjusting probabilities:' 
       write(*,*) 'gridpoint cptpp1 p03m01 prct(k,1),prct(k,2),grd(k)'
        do k=1,jf
         if(k.eq.1562) then
          write(*,'(i10,5f7.2)')k,cptpp1(K),p03m01(k),
     +     prct(k,1),prct(k,2),grd(k)
         end if
        end do
        
       write(*,*) 'Done adjusting probabilities'
       return      
       end


