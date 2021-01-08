c       program calibrate2_hrly_rgn3
c
c Just like calibrate2_hrly, except reads the gridded output for "regional" coefficients.  DRB, 8/18/2004.
c Just like calibrate.f, except reads the output of the combined probability table.

c  calibrate2_hrly_rgn3: to compute calibrated cloud phy thunder parameter(cptp) prob
c                   adapted from David Bright (SPC/NCEP)
c   In compare to calibrate2, calibrate2_hrly_rgn3 use cycle and forecast hour depedent
c     calibration file (combine_pops_rgn3_$cyc_$fhr.out) instead of constant 'combine_pops.out' file
c
c   Author: Binbin Zhou/IMSG
c           July 6 2011
c   Input: cptpp1,p03m01,igemy,igemx,jf
c   Output grd
c

       subroutine calibrate2_hrly_rgn3 (cptpp1,p03m01,igemx,igemy,
     +    jf,cyc,fhr,grd)

c       parameter (igemy=129, igemx=185) ! GEMPAK grid size
c
       real,dimension(jf),intent(IN) :: cptpp1,p03m01
       real grd(jf)
       real prct(igemy,igemx,2),prob,pa(121),pb(121),rk,
     &       grd2(igemx,igemy),dum1,dum2
       real,allocatable,dimension(:,:) :: bin
       real,allocatable,dimension(:,:,:) :: bias
       real,allocatable,dimension(:,:,:,:) :: hit

       integer igemy,igemx,itcnt,jstop,jdel
       character*150 pctfila,pctfilb,combine_pops
       character dummy*20,date*8,fhour*4,frun*2,dateetc*72
       character*2 cyc,fhr             

c fill all possible forecast values...	 
       itcnt = 0
       do i=0,100,10
        do j = 0,100,10
	 itcnt = itcnt + 1
         pa(itcnt) = float(i) ! fcst1
  	 pb(itcnt) = float(j) ! fcst2
        enddo
       enddo

       write(*,*) ' calibration2_hrly_rgn3: cyc, fhr=', cyc, fhr
       write(*,*) 'Number of possible combinations... ',itcnt

        do j=1,igemy
         do i=1,igemx
          ij=i+igemx*(j-1)
           prct(j,i,1)=cptpp1(ij)
           prct(j,i,2)=p03m01(ij)
         end do
        end do

c Now the grids are in...        
c Do a little work here.  First, bin probs to 10% intervals.  Then, change verification
c to 1's or 0's.  

       do i=1,igemx
       do j=1,igemy
        do k=1,2
         if(prct(j,i,k).lt.-1.0) prct(j,i,k) = -9999.0
         if(prct(j,i,k).ge.-1.0.and.prct(j,i,k).lt.5.0)prct(j,i,k)= 0.0
         if(prct(j,i,k).ge.5.0.and.prct(j,i,k).lt.15.) prct(j,i,k)=10.0
         if(prct(j,i,k).ge.15.0.and.prct(j,i,k).lt.25.)prct(j,i,k)=20.0
         if(prct(j,i,k).ge.25.0.and.prct(j,i,k).lt.35.)prct(j,i,k)=30.0
         if(prct(j,i,k).ge.35.0.and.prct(j,i,k).lt.45.)prct(j,i,k)=40.0
         if(prct(j,i,k).ge.45.0.and.prct(j,i,k).lt.55.)prct(j,i,k)=50.0
         if(prct(j,i,k).ge.55.0.and.prct(j,i,k).lt.65.)prct(j,i,k)=60.0
         if(prct(j,i,k).ge.65.0.and.prct(j,i,k).lt.75.)prct(j,i,k)=70.0
         if(prct(j,i,k).ge.75.0.and.prct(j,i,k).lt.85.)prct(j,i,k)=80.0
         if(prct(j,i,k).ge.85.0.and.prct(j,i,k).lt.95.)prct(j,i,k)=90.0
         if(prct(j,i,k).ge.95.0) prct(j,i,k)=100.0
        enddo
       enddo
       enddo
c
c Read the file that contains the dependency coefficient and the bias... 
c

         allocate (bin(igemx,igemy))
         allocate (bias(121,igemx,igemy))
         allocate (hit(121,2,igemx,igemy))
 
         combine_POps='combine_pops_rgn3_'//cyc//'_f'//fhr//'.out'
         write(*,*) 'combine_pops=',combine_pops
  
         open(unit=301,file=combine_pops,status='old',
     &      form='unformatted')
          read(301) pa
	  read(301) pb
	  read(301) bin       
	  read(301) hit       ! This is the reason for blank strip in lightning plots
	  close(301)

        do ii = 1,igemx
	do jj = 1,igemy
	do kk = 1,itcnt
	  if(hit(kk,1,ii,jj).lt.-9000.0.or.
     &        hit(kk,2,ii,jj).lt.-9000.0.or.
     &        bin(ii,jj).lt.-9000.0) then
             bias(kk,ii,jj) = -9999.0
	  elseif(hit(kk,1,ii,jj).lt.0.0e-15.and.
     &           hit(kk,2,ii,jj).lt.0.0e-15) then
	   iiii = -9999
	   jjjj = -9999
c find closest grid point in the area with good data...	   
	   do jjj=jj-20,jj+20
	   loop801: do iii=ii-20,ii+20

	   if(iii.lt.1.or.jjj.lt.1.or.iii.gt.igemx.or.jjj.gt.igemy) 
     &        cycle loop801	   
	   if(hit(kk,1,iii,jjj).ge.0.0e-15.and.
     &         hit(kk,2,iii,jjj).ge.0.0e-15)then
             if(abs(ii-iii).lt.abs(ii-iiii).and.
     &           abs(jj-jjj).lt.abs(jj-jjjj)) then
	      iiii = iii
	      jjjj = jjj
	     endif	          
	   endif
	   enddo loop801
	   enddo
	   if(iiii.gt.-9000.and.jjjj.gt.-9000) then
	   dum1= (hit(kk,1,iiii,jjjj)/bin(iiii,jjjj))*100.0
	   dum2= (hit(kk,2,iiii,jjjj)/bin(iiii,jjjj))*100.0	  
	   bias(kk,ii,jj) = (dum2/(dum1+dum2))*100.
	   else
	   bias(kk,ii,jj) = -9999.0
	   endif
	  else 
	   dum1= (hit(kk,1,ii,jj)/bin(ii,jj))*100.0
	   dum2= (hit(kk,2,ii,jj)/bin(ii,jj))*100.0	  
	   bias(kk,ii,jj) = (dum2/(dum1+dum2))*100.
	  endif
c	  bias(ii) = (dum2/(dum1+dum2))*100.
	enddo
	enddo
        enddo

        deallocate(bin)
        deallocate(hit)

c Now...bin the forecasts into all psbl combinations.  Record if a hit occurred
c during this combined period too.
c


       do i=1,igemx
       do j=1,igemy
       ii = 0
       loop1001: do k=1,itcnt
         if(nint(prct(j,i,1)).eq.nint(pa(k)).and.
     &     nint(prct(j,i,2)).eq.nint(pb(k))) then
            ii = k
            exit loop1001
         endif
       enddo loop1001
 
        if(prct(j,i,1).lt.-9990.0.or.prct(j,i,2).lt.-9990.0.or.
     &      ii.eq.0) then
	 grd2(i,j) = -9999.0
	else
c 	 grd2(i,j) = bias(ii)
 	 grd2(i,j) = bias(ii,i,j)
	endif
	
       enddo ! j loop
       enddo ! i loop
 
       deallocate(bias)

       do j= 1, igemy
        do i =1, igemx
          ij = i + igemx*(j-1)
          grd(ij)=grd2(i,j)
          if (grd(ij).le.-999.0) grd(ij)=0.0
        end do
       end do  
             
       write(*,*) 'calibrate2_hrly_rgn3 DONE!'     

       return
       end


