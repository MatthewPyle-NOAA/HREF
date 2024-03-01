        subroutine read_grib1(datstr)

!       Simple program that reads in GRIB1 individual RFC files
!       and extracts needed FFG data, and writes out
!       as GRIB2 encoded.  A preliminary step before each
!       GRIB2 subpiece is interpolated (at script level) onto
!       the HREF 5 km grid 227.  Those are ingested by a second
!       program that stitches them together onto a single 
!       mosaic'd product.

        USE GRIB_MOD


! grib2
      INTEGER :: LUGB,LUGI,J,JDISC,JPDTN,JGDTN
      INTEGER,DIMENSION(:) :: JIDS(200),JPDT(200),JGDT(200)
      LOGICAL :: UNPACK
      INTEGER :: K,IRET
      TYPE(GRIBFIELD) :: GFLD, GFLD_FULL
! grib2


      integer jpds(200),jgds(200),kpds(200),kgds(200)
      integer :: nx, ny, nxny

      logical*1, allocatable:: li(:)
      real, allocatable:: ffg1(:), ffg3(:)
      real, allocatable:: ffg6(:), ffg12(:)
      real, allocatable:: ffg24(:)

      real, allocatable:: ffg1_full(:), ffg3_full(:)
      real, allocatable:: ffg6_full(:), ffg12_full(:)
      real, allocatable:: ffg24_full(:)

      character*80 fname,prefx,filnam,fnameg2
      character*80 fnameg2out
      character*8 datstr

      character(len=3) :: rfc_str(12)
      DATA rfc_str /'150','152','153','154','155','156', &
                    '157','158','159','160','161','162'/
!      fname='ffg.20190326.009.156'

!      read(5,FMT='(A)') datstr
     

   RFC_LOOP:      DO II=1,12  ! loop over RFC regions

      iunit=II+11
      iunitg2=II+23

      fname='ffg.'//datstr//'.009.'//rfc_str(II)

      write(0,*) 'computed fname as: ', trim(fname)
      fnameg2=trim(fname)//'.g2'
      iunitg2out=II+35
      fnameg2out=trim(fname)//'.g2out'


        write(0,*) 'fname, fnameg2, fnameg2out: ', & 
                   fname, fnameg2, fnameg2out

        jpds=-1
        jgds=-1

       call baopenr(iunit,fname,ierr)
!       write(0,*) 'ierr from baopenr: ', ierr
	if (ierr .ne. 0) then
	write(0,*) 'missing file to read  in read_grib1...skip'
        CYCLE RFC_LOOP
        endif
       call baopenr(iunitg2,fnameg2,ierr)
!       write(0,*) 'fnameg2 ierr from baopenr: ', ierr
       call baopenw(iunitg2out,fnameg2out,ierr)
!       write(0,*) 'fnameg2out ierr from baopenw: ', ierr

        JIDS=-9999
        JPDT=-9999
        JGDT=-9999

         call getgb2(iunitg2,0,0,-1,jids,-1,jpdt,-1,jgdt,.true., &
     &               K,gfld,iret)  

!           write(0,*) 'iret from getgb2: ', iret
!           write(0,*) 'K returned: ', K
!            write(0,*) 'gfld%ngrdpts: ', gfld%ngrdpts
!           write(0,*) 'gfld%discipline: ', gfld%discipline

!  reset ipdtmpl(24) so knows is an accumulation (not missing)
           gfld%ipdtmpl(24)=1
           gfld%discipline=1
!           write(0,*) 'gfld%ipdtnum: ', gfld%ipdtnum
!           write(0,*) 'gfld%ipdtmpl(1:20): ', gfld%ipdtmpl(1:20)

           gfld%ipdtmpl(1)=0
           gfld%ipdtmpl(2)=0

        jpds(5)=221
        call getgb(iunit,0,0,0,jpds,jgds,kf,kr,kpds,kgds,li,ffg1,irets)

        write(0,*) 'kgds(2:3): ', kgds(2:3)
         nx=kgds(2)
         ny=kgds(3)

         nxny=nx*ny

        if (.not. allocated(li)) then
        allocate(li(nxny))
        allocate (ffg1(nxny))
        allocate (ffg3(nxny))
        allocate (ffg6(nxny))
        allocate (ffg12(nxny))
        allocate (ffg24(nxny))
        endif

!        write(0,*) 'shape(ffg1): ', shape(ffg1)
        jpds(5)=221
        call getgb(iunit,0,nxny,0,jpds,jgds,kf,kr,kpds,kgds,li, &
     &             ffg1,irets)
        write(6,*) 'getgb ', fname, ' irets=', irets
         if (irets .eq. 0) then
        write(0,*) '221 - maxval(ffg1): ', maxval(ffg1)
        gfld%fld=ffg1
        gfld%ipdtmpl(27)=1
        gfld%ipdtmpl(9)=23
        call putgb2(iunitg2out,gfld,iret)
!        write(0,*) 'iret from putgb2: ', iret
         endif
    

        jpds(5)=222
        call getgb(iunit,0,nxny,0,jpds,jgds,kf,kr,kpds,kgds,li, &
     &            ffg3,irets)
!        write(6,*) 'getgb ', fname, ' irets=', irets
         if (irets .eq. 0) then
        write(0,*) '222 - maxval(ffg3): ', maxval(ffg3)
        gfld%fld=ffg3
        gfld%ipdtmpl(27)=3
        gfld%ipdtmpl(9)=21
        call putgb2(iunitg2out,gfld,iret)
         endif

        jpds(5)=223
        call getgb(iunit,0,nxny,0,jpds,jgds,kf,kr,kpds,kgds,li, &
     &             ffg6,irets)
!        write(6,*) 'getgb ', fname, ' irets=', irets
         if (irets .eq. 0) then
        write(0,*) '223 - maxval(ffg6): ', maxval(ffg6)
        gfld%fld=ffg6
        gfld%ipdtmpl(27)=6
        gfld%ipdtmpl(9)=18
        call putgb2(iunitg2out,gfld,iret)
         endif

        jpds(5)=224
        call getgb(iunit,0,nxny,0,jpds,jgds,kf,kr,kpds,kgds,li, &
     &             ffg12,irets)
!        write(6,*) 'getgb ', fname, ' irets=', irets
         if (irets .eq. 0) then
        write(0,*) '224 - maxval(pcp1): ', maxval(ffg12)
        gfld%fld=ffg12
        gfld%ipdtmpl(27)=12
        gfld%ipdtmpl(9)=12
        call putgb2(iunitg2out,gfld,iret)
         endif

        jpds(5)=225
        call getgb(iunit,0,nxny,0,jpds,jgds,kf,kr,kpds,kgds,li,  &
     &             ffg24,irets)
!        write(6,*) 'getgb ', fname, ' irets=', irets
         if (irets .eq. 0) then
        write(0,*) '225 - maxval(pcp1): ', maxval(ffg24)
        gfld%fld=ffg24
        gfld%ipdtmpl(27)=24
        gfld%ipdtmpl(9)=0
        call putgb2(iunitg2out,gfld,iret)
         endif

        call baclose(iunit,ierr)

        deallocate(li)
        deallocate(ffg1)
        deallocate(ffg3)
        deallocate(ffg6)
        deallocate(ffg12)
        deallocate(ffg24)

         enddo RFC_LOOP ! for RFC regions

        end subroutine read_grib1
