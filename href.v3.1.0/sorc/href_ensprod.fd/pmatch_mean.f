	subroutine pmatch_mean(vname,rawdata_mn,vrbl_mn,
     &         jpd1,jpd2,jpd10,vrbl_mn_pm,lv,lm,jf,iens)

        character(len=4) :: vname
	real :: rawdata_mn(jf,iens,lm),vrbl_mn(jf,lm),vrbl_mn_pm(jf,lm)
        real, allocatable :: rawdata_mn_loc(:,:,:)

	real, allocatable :: rawdata_1d(:),vrbl_mn_hold(:,:)
        integer, allocatable :: listorderfull(:),listorder(:)
        integer :: iplace,lm,lf,jf,JJ
        integer :: ibound_max, ibound_min
        real:: amin,amax

        allocate(listorderfull(jf*iens))
        allocate(listorder(jf))
	allocate(rawdata_1d(jf*iens))
        allocate(vrbl_mn_hold(jf,lm))
        allocate(rawdata_mn_loc(jf,iens,lm))

	write(0,*) 'pmatch_mean iens, jf, lm: ', iens, jf, lm

	write(0,*) 'jpd1, jpd2: ', jpd1, jpd2

! set local copy to rawdata_mn
	rawdata_mn_loc=rawdata_mn

        do I=1,iens
        do J=1,jf
        listorder(J)=J
        listorderfull((I-1)*jf+J)=(I-1)*jf+J
	
! force reflectivity type fields to be zero
         if(jpd1.eq.16.and.(jpd2.eq.195 
     &                 .or. jpd2.eq.196
     &                 .or. jpd2.eq.198) .and. 
     &   rawdata_mn(J,I,lv) .lt. 0.) then

        rawdata_1d((I-1)*jf+J)=0.
        rawdata_mn_loc(J,I,lv)=0.
        else
        rawdata_1d((I-1)*jf+J)=rawdata_mn(J,I,lv)
        endif

        enddo
        enddo

        vrbl_mn_hold=vrbl_mn

        call quick_sort(vrbl_mn,listorder,jf)
        call quick_sort(rawdata_1d,listorderfull,iens*jf)

        write(0,*) 'min,maxval(rawdata_mn): ',minval(rawdata_mn(:,:,lv))
     &                                       ,maxval(rawdata_mn(:,:,lv))
        write(0,*) 'min,maxval(rawdata_mn_loc): ', 
     &  minval(rawdata_mn_loc(:,:,lv)),maxval(rawdata_mn_loc(:,:,lv))

        write(0,*) 'minval(vrbl_mn(:,lv)): ',minval(vrbl_mn(:,lv))
        write(0,*) 'maxval(vrbl_mn(:,lv)): ',maxval(vrbl_mn(:,lv))

        write(0,*) 'minval(rawdata_1d(:)): ', minval(rawdata_1d(:))
        write(0,*) 'maxval(rawdata_1d(:)): ', maxval(rawdata_1d(:))

        vrbl_mn_pm(:,lv)=-999.

        ibound_max=0
        ibound_min=0

      ens_loop:  do J=1,jf*(iens),iens    ! loop over full ensemble, skipping

         I=1+(J-1)/(iens)
         iplace=listorder(I)

!!!  use unsorted version if looking at iplace

         if(jpd1.eq.1.and.jpd2.eq.8.and.jpd10.eq.1.and.     !APCP
     +     vrbl_mn_hold(iplace,lv).eq.0) then
          vrbl_mn_pm(iplace,lv)=0.
          cycle ens_loop
         end if

         amin=9999. 
         amax=-9999.
         do JJ=1,iens
          if (rawdata_mn_loc(iplace,JJ,lv) .gt. amax) then
                amax=rawdata_mn_loc(iplace,JJ,lv)
          endif
          if (rawdata_mn_loc(iplace,JJ,lv) .lt. amin) then
                amin=rawdata_mn_loc(iplace,JJ,lv)
          endif
         enddo
 
!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (rawdata_1d(J) .gt. amax) then

          vrbl_mn_pm(iplace,lv)=amax 
          ibound_max=ibound_max+1

         elseif (rawdata_1d(J) .lt. amin) then

          vrbl_mn_pm(iplace,lv)=amin
          ibound_min=ibound_min+1

         else

          vrbl_mn_pm(iplace,lv)=rawdata_1d(J)

         endif

      enddo ens_loop

! restore the mean value for use in possible blending

        vrbl_mn=vrbl_mn_hold

	write(0,*) 'ibound_min: ', ibound_min
	write(0,*) 'ibound_max: ', ibound_max

	deallocate(listorderfull)
        deallocate(listorder)
        deallocate(rawdata_1d)
        deallocate(vrbl_mn_hold)
        deallocate(rawdata_mn_loc)

	end subroutine pmatch_mean
