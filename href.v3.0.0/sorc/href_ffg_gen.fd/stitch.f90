        subroutine stitch(datstr)

        USE GRIB_MOD

! grib2
      INTEGER :: LUGB,LUGI,J,JDISC,JPDTN,JGDTN
      INTEGER,DIMENSION(:) :: JIDS(200),JPDT(200),JGDT(200)
      LOGICAL :: UNPACK
      INTEGER :: K,IRET
      TYPE(GRIBFIELD) :: gfld, gfld_full
! grib2

      integer jpds(200),jgds(200),kpds(200),kgds(200)
      integer :: nx, ny, nxny, FFG, KK

      character*80 fname,prefx,filnam,fnameg2
      character*80 fnameout,fnamein
      character*8 datstr, ffgstr

      character(len=3) :: rfc_str(12)
      DATA rfc_str /'150','152','153','154','155','156', &
                    '157','158','159','160','161','162'/
!      fname='ffg.20190326.009.156'

!      read(5,FMT='(A)') datstr
     
      DO FFG=221,223

        if (FFG .eq. 221) ffgstr='ffg1h'
        if (FFG .eq. 222) ffgstr='ffg3h'
        if (FFG .eq. 223) ffgstr='ffg6h'
        if (FFG .eq. 224) ffgstr='ffg12h'
        if (FFG .eq. 225) ffgstr='ffg24h'

RFC_LOOP:    DO II=1,12  ! loop over RFC regions

	write(0,*) 'start loop with FFG, II: ', FFG, II
      iunit=II+11+(FFG-221)*12

      fname='ffg.'//datstr//'.009.'//rfc_str(II)
      fnamein=trim(fname)//'.g2out.227'
      write(0,*) 'computed fnamein as: ', trim(fnamein)

        if (II .eq. 1) then
      iunitg2out=77
      fnameout='full.g227.grib2_'//trim(ffgstr)
      write(0,*) 'computed fnameout as: ', trim(fnameout)
       call baopen(iunitg2out,fnameout,ierr)
       write(0,*) 'fnameout ierr from baopen: ', ierr
        endif

        jpds=-1
        jgds=-1

       call baopenr(iunit,fnamein,ierr)
	if (ierr .ne. 0) then
       write(0,*) 'ierr from baopenr: ', ierr
	write(0,*) 'skill skip this region as have no input'
        CYCLE RFC_LOOP
	endif
        
        JIDS=-9999
        JPDT=-9999
        JGDT=-9999

         call getgb2(iunit,0,0,-1,jids,-1,jpdt,-1,jgdt,.true., &
     &               K,gfld,iret)  
        if (II .eq. 1) then
         call getgb2(iunit,0,0,-1,jids,-1,jpdt,-1,jgdt,.true., &
     &               K,gfld_full,iret)  
          gfld_full%bmap(:)=.false.
        endif

           write(0,*) 'iret from getgb2: ', iret
!           write(0,*) 'K returned: ', K
!            write(0,*) 'gfld%ngrdpts: ', gfld%ngrdpts
!           write(0,*) 'gfld%discipline: ', gfld%discipline
            write(0,*) 'gfld_full%igdtmpl(8:9): ', gfld_full%igdtmpl(8:9)

!  reset ipdtmpl(24) so knows is an accumulation (not missing)
           gfld_full%ipdtmpl(24)=1
           gfld_full%discipline=1
!           write(0,*) 'gfld%ipdtnum: ', gfld%ipdtnum
!           write(0,*) 'gfld%ipdtmpl(1:20): ', gfld%ipdtmpl(1:20)

           gfld_full%ipdtmpl(1)=0
           gfld_full%ipdtmpl(2)=0

         nxny=gfld_full%igdtmpl(8)*gfld_full%igdtmpl(9)

        write(0,*) 'nxny for 227 grid: ', nxny



!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (FFG .eq. 221) then

         jpdt(27)=1
         jpdt(9)=23
         call getgb2(iunit,0,0,-1,jids,8,jpdt,-1,jgdt,.true., &
     &               K,gfld,iret)  
         if (iret .eq. 0) then
        write(0,*) '221 - maxval(ffg): ', maxval(gfld%fld)
        gfld_full%ipdtmpl(27)=gfld%ipdtmpl(27)
        gfld_full%ipdtmpl(9)=gfld%ipdtmpl(9)

        do J=1,nxny

!tst        if (gfld%bmap(J) .and. .not. gfld_full%bmap(j)) then
        if (gfld%bmap(J)) then
        gfld_full%fld(J)=gfld%fld(J)
        gfld_full%bmap(J)=.true.
!        if (mod(J,500) .eq. 0 .or. J .eq. 888959) then
!        write(0,*) 'FFG1 defined II,J ',II, J, gfld_full%fld(j)
!        endif
        endif
        enddo

         endif ! on iret=0
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (FFG .eq. 222) then

!        write(6,*) 'getgb ', fname, ' irets=', irets
         jpdt(27)=3
         jpdt(9)=21
	write(0,*) 'call getgb2 for iunit: ', iunit
         call getgb2(iunit,0,0,-1,jids,8,jpdt,-1,jgdt,.true., &
     &               K,gfld,iret)  

         if (iret .eq. 0) then
        write(0,*) '222 - maxval(ffg): ', maxval(gfld%fld)
        gfld_full%ipdtmpl(27)=gfld%ipdtmpl(27)
        gfld_full%ipdtmpl(9)=gfld%ipdtmpl(9)

        do J=1,nxny
        if (gfld%bmap(J)) then
        gfld_full%fld(J)=gfld%fld(J)
        gfld_full%bmap(J)=.true.
!        if (mod(J,500) .eq. 0) then
!        write(0,*) 'FFG3 defined II, J ',II, J, gfld_full%fld(j)
!        endif
        endif
        enddo
         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (FFG .eq. 223) then
         jpdt(27)=6
         jpdt(9)=18
         call getgb2(iunit,0,0,-1,jids,8,jpdt,-1,jgdt,.true., &
     &               K,gfld,iret)  
         if (iret .eq. 0) then
        write(0,*) '223 - maxval(ffg): ', maxval(gfld%fld)

        gfld_full%ipdtmpl(27)=gfld%ipdtmpl(27)
        gfld_full%ipdtmpl(9)=gfld%ipdtmpl(9)

        do J=1,nxny
        if (gfld%bmap(J)) then
        gfld_full%fld(J)=gfld%fld(J)
        gfld_full%bmap(J)=.true.
        endif
        enddo

         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (FFG .eq. 224) then
         jpdt(27)=12
         jpdt(9)=12
         call getgb2(iunit,0,0,-1,jids,8,jpdt,-1,jgdt,.true., &
     &               K,gfld,iret)  
         if (iret .eq. 0) then
        write(0,*) '224 - maxval(ffg): ', maxval(gfld%fld)
        gfld_full%ipdtmpl(27)=12
        gfld_full%ipdtmpl(9)=12
        do J=1,nxny
        if (gfld%bmap(J)) then
        gfld_full%fld(J)=gfld%fld(J)
        gfld_full%bmap(J)=.true.
        endif
        enddo
         endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (FFG .eq. 225) then

         jpdt(27)=24
         jpdt(9)=0
         call getgb2(iunit,0,0,-1,jids,8,jpdt,-1,jgdt,.true., &
     &               K,gfld,iret)  
         if (irets .eq. 0) then
        write(0,*) '225 - maxval(ffg): ', maxval(gfld%fld)
        gfld_full%ipdtmpl(27)=24
        gfld_full%ipdtmpl(9)=0
        do J=1,nxny
        if (gfld%bmap(J)) then
        gfld_full%fld(J)=gfld%fld(J)
        gfld_full%bmap(J)=.true.
        endif
        enddo

         endif

        endif ! FFG test

        call baclose(iunit,ierr)

         enddo RFC_LOOP ! for RFC regions

        write(0,*) 'here with maxval(gfld_full): ', &
         maxval(gfld_full%fld)



!! fill in some missing values (mostly or entirely water points??)
        do KK=1,3

        do J=2,nxny-1

! missing one

        if (gfld_full%bmap(J-1) .and. gfld_full%bmap(J+1) .and. .not. &
            gfld_full%bmap(J)) then

         write(0,*) '1 - reset J to true: ',KK, J
         gfld_full%bmap(J)=.true.
         gfld_full%fld(J)=0.5*(gfld_full%fld(J-1)+gfld_full%fld(J+1))
         write(0,*) 'defined gfld_full%fld: ', KK,J, gfld_full%fld(J)
         endif
         enddo

! this value (1746) is grid specific
        do J=3+1746,nxny-2-1746

        if ((gfld_full%bmap(J-1746).and.gfld_full%bmap(J+1746)).and. &
            .not. gfld_full%bmap(J)) then
         write(0,*) '2 - reset J to true: ',KK, J
         gfld_full%bmap(J)=.true.
         gfld_full%fld(J)=0.5*(gfld_full%fld(J-1746)+gfld_full%fld(J+1746))
         write(0,*) '2 - defined gfld_full%fld: ', KK,J, gfld_full%fld(J)
        endif

         enddo

        do J=3,nxny-2
 
        if ( (gfld_full%bmap(J-2) .or. gfld_full%bmap(J-1)) .and. &
             (gfld_full%bmap(J+1) .or. gfld_full%bmap(J+2)) &
             .and. .not. gfld_full%bmap(J)) then             
         write(0,*) '3 - reset J to true: ', KK, J
         gfld_full%bmap(J)=.true.
         gfld_full%fld(J)=maxval(gfld_full%fld(J-2:J+2))
         write(0,*) '3 - defined gfld_full%fld: ', KK,J, gfld_full%fld(J)
         endif

         enddo

        enddo ! KK

	write(0,*) 'to putgb2 call'
	write(0,*) 'maxval: ', maxval(gfld_full%fld)
         call putgb2( iunitg2out, gfld_full, iret)
        write(0,*) 'iret from putgb2 call: ', iret
        call baclose(iunitg2out,ierr)

          END DO  ! FFG

        end subroutine stitch
