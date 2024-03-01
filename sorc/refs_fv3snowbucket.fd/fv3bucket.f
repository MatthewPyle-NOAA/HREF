	program fv3bucket

!       Special program to get 1-h and 3-h snow from FV3SAR instantaneous WEASD



        USE GRIB_MOD

        real, allocatable :: pdiff(:,:)

        integer :: ihrs1, ihrs2, reset_flag
        integer :: mm,nn,oo,m,n
	character(len=2), dimension(2):: hrs
        integer :: is_hrrr

	character(len=255):: file1,file2,file3,file4,file5,testout
	character(len=150):: dirname
	character(len=150):: filename
	character(len=1):: reflag

	read(5,FMT='(A)') dirname
	read(5,FMT='(A)') filename
        read(5,FMT='(A)') hrs(1)
        read(5,FMT='(A)') hrs(2)
        read(5,FMT='(I1)') reset_flag
        read(5,*) IM, JM
        read(5,FMT='(I)') is_hrrr

        allocate(pdiff(im,jm))

	write(0,*) 'read reset_flag: ', reset_flag

	n=index(dirname,' ')-1
	m=index(filename,' ')-1

	I=2
	
	file1= dirname(1:n)//'/'//filename(1:m)
     +					//HRS(I-1)//'.grib2'
	file2= dirname(1:n)//'/'//filename(1:m)//HRS(I)//'.grib2'

	read(HRS(I-1), '(I2)' ) ihrs1
	read(HRS(I)  , '(I2)' ) ihrs2

	interv=ihrs2-ihrs1

        write(0,*) 'ihrs1, ihrs2: ', ihrs1, ihrs2
	
	if (interv .eq. 3) then
	  testout= dirname(1:n)//'/PCP3HR'//HRS(I)//'.tm00'
	elseif (interv .eq. 1) then
	  testout= dirname(1:n)//'/PCP1HR'//HRS(I)//'.tm00'
	endif

	mm=index(file1,' ')-1
	nn=index(file2,' ')-1

	mmm=index(testout,' ')-1

	write(0,*) 'file1: ', file1(1:mm)
	write(0,*) 'file2: ', file2(1:nn)
	write(0,*) 'testout: ', testout(1:mmm)

	if (mod(ihrs1,3) .eq. 0 .and. reset_flag .eq. 1) then
	reset_flag=1
	else
	reset_flag=0
	endif

        write(0,*) 'call calc_pdiff with reset_flag: ', reset_flag

	call calc_pdiff(file1(1:mm),file2(1:nn),
     &                  TESTOUT(1:mmm),
     &                  pdiff,reset_flag,ihrs1,interv,IM,JM,is_hrrr)

        write(0,*) 'past calc_pdiff'


	END program fv3bucket

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE CALC_PDIFF(FNAME1,FNAME2,TESTOUT,SPRECIP,
     &                   reset_flag,ihrs1,interv,IM,JM,is_hrrr)
        USE GRIB_MOD
        USE pdstemplates
	character(*):: FNAME1,FNAME2,testout
	integer:: reset_flag,lencheck,ihrs1,interv,is_hrrr
	logical:: FIRST

        real:: rinc(5),sprecip(IM*JM)
        real:: asnowprecip(IM*JM)
        real:: fzprecip(IM*JM)
        integer:: idat(8),jdat(8)

        INTEGER :: IBM,LEN,ISCALE,NBIT,IBITM
        INTEGER, allocatable :: ibmap(:)
        REAL, allocatable :: GROUND(:)
        REAL :: GMIN, GMAX, SGDS

C grib2
      INTEGER :: LUGB,LUGI,J,JDISC,JPDTN,JGDTN
      INTEGER,DIMENSION(:) :: JIDS(200),JPDT(200),JGDT(200)
      logical*1, allocatable:: bmap_f(:)
      LOGICAL :: UNPACK
      INTEGER :: K,IRET
      TYPE(GRIBFIELD) :: GFLD, GFLD_QPF
C grib2
        real, allocatable :: s_later(:),s_earlier(:)
        real, allocatable :: asnow_later(:),asnow_earlier(:)
        real, allocatable :: fz_later(:),fz_earlier(:)
	

	call baopenr(11,fname1,ierr1)
	call baopenr(12,fname2,ierr2)
	call baopenw(13,testout,ierr3)

	write(0,*) 'baopened ', fname1,fname2,testout

	if ( (ierr1+ierr2+ierr3) .ne. 0) then
		write(0,*) 'bad baopen!!! ', ierr1
		write(0,*) 'bad baopen!!! ', ierr2
		write(0,*) 'bad baopen!!! ', ierr3
		STOP 99
	endif

        allocate(s_earlier(IM*JM))
        allocate(s_later(IM*JM))
        allocate(fz_earlier(IM*JM))
        allocate(fz_later(IM*JM))
        allocate(asnow_earlier(IM*JM))
        allocate(asnow_later(IM*JM))
        allocate(GROUND(IM*JM))
        allocate(ibmap(IM*JM))
	
        s_earlier=0.

        JIDS=-9999
        JPDTN=-1
        JPDT=-9999
        JGDTN=-1
        JGDT=-9999
        UNPACK=.true.

        if ( .not. allocated(bmap_f)) then
        allocate(bmap_f(IM*JM))
        endif

! -------------------------------

        if (is_hrrr .eq. 0) then

        J=0
        JIDS=-9999
        JPDTN=0
        JPDT=-9999
        JPDT(2)=13
        JGDTN=-1
        JGDT=-9999
        UNPACK=.true.

        call getgb2(11,0,J,0,JIDS,JPDTN,JPDT,JGDTN,JGDT,
     &     UNPACK,K,GFLD,IRET)

	if (IRET .ne. 0) then
	write(0,*) 'bad getgb1 (13)  earlier ', IRET
	STOP 999
	endif

        s_earlier=gfld%fld
        write(0,*) 'maxval(s_earlier): ', maxval(s_earlier)

         endif

        if (is_hrrr .eq. 1) then

        J=0
        JIDS=-9999
        JPDTN=8
        JPDT=-9999
        JPDT(2)=225
        JGDTN=-1
        JGDT=-9999
        UNPACK=.true.

        call getgb2(11,0,J,0,JIDS,JPDTN,JPDT,JGDTN,JGDT,
     &     UNPACK,K,GFLD,IRET)

	if (IRET .ne. 0) then
	write(0,*) 'bad getgb1 (225)  earlier ', IRET
	STOP 999
	endif

        fz_earlier=gfld%fld
        write(0,*) 'maxval(fz_earlier): ', maxval(fz_earlier)

         endif

! -------------------------------

        J=0
        JIDS=-9999
        JPDTN=8
        JPDT=-9999
        JPDT(2)=29
        JGDTN=-1
        JGDT=-9999
        UNPACK=.true.

        call getgb2(11,0,J,0,JIDS,JPDTN,JPDT,JGDTN,JGDT,
     &     UNPACK,K,GFLD,IRET)

	if (IRET .ne. 0) then
	write(0,*) 'bad getgb1 (29) earlier ', IRET
	STOP 999
	endif

        asnow_earlier=gfld%fld
        write(0,*) 'maxval(asnow_earlier): ', maxval(asnow_earlier)

! -------------------------------

        if (is_hrrr .eq. 0) then

        J=0
        JIDS=-9999
        JPDTN=0
        JPDT=-9999
        JPDT(2)=13
        JGDTN=-1
        JGDT=-9999

        call getgb2(12,0,0,0,JIDS,JPDTN,JPDT,JGDTN,JGDT,
     &     UNPACK,K,GFLD,IRET1)

        write(0,*) 'K: ', K

	if (IRET1 .ne. 0) then
	 write(0,*) 'bad getgb later ', IRET1
	STOP 9999
	endif

        write(0,*) 'pulled gfld%ipdtnum: ', gfld%ipdtnum
!        do J=1,29
!        write(0,*) 'ingest J,gfld%ipdtmpl(J): ', J,gfld%ipdtmpl(J)
!        enddo

        bmap_f=gfld%bmap

        write(0,*) 'set s_later to gfld%fld'
        s_later=gfld%fld
        write(0,*) 'maxval(s_later): ', maxval(s_later)

          endif

        if (is_hrrr .eq. 1) then

        J=0
        JIDS=-9999
        JPDTN=8
        JPDT=-9999
        JPDT(2)=225
        JGDTN=-1
        JGDT=-9999

        call getgb2(12,0,0,0,JIDS,JPDTN,JPDT,JGDTN,JGDT,
     &     UNPACK,K,GFLD,IRET1)

        write(0,*) 'K: ', K

	if (IRET1 .ne. 0) then
	 write(0,*) 'bad getgb later (225) ', IRET1
	STOP 9999
	endif

        write(0,*) 'pulled gfld%ipdtnum: ', gfld%ipdtnum
!        do J=1,29
!        write(0,*) 'ingest J,gfld%ipdtmpl(J): ', J,gfld%ipdtmpl(J)
!        enddo

        bmap_f=gfld%bmap

        write(0,*) 'set fz_later to gfld%fld'
        fz_later=gfld%fld
        write(0,*) 'maxval(fz_later): ', maxval(fz_later)

          endif

! -------------------------------

        J=0
        JIDS=-9999
        JPDTN=8
        JPDT=-9999
        JPDT(2)=29
        JGDTN=-1
        JGDT=-9999

        call getgb2(12,0,0,0,JIDS,JPDTN,JPDT,JGDTN,JGDT,
     &     UNPACK,K,GFLD,IRET1)

        write(0,*) 'K: ', K

	if (IRET1 .ne. 0) then
	 write(0,*) 'bad getgb later ', IRET1
	STOP 9999
	endif

        write(0,*) 'pulled gfld%ipdtnum: ', gfld%ipdtnum

        bmap_f=gfld%bmap

        write(0,*) 'set asnow_later to gfld%fld'
        asnow_later=gfld%fld
        write(0,*) 'maxval(asnow_later): ', maxval(asnow_later)

! -------------------------------

! Pull in QPF to get more proper PDS for WEASD bucket

        JIDS=-9999
        JPDTN=8
        JPDT=-9999
        JPDT(2)=8
        JGDTN=-1
        JGDT=-9999

        call getgb2(12,0,0,0,JIDS,JPDTN,JPDT,JGDTN,JGDT,
     &     UNPACK,K,GFLD_QPF,IRET1)

        allocate(gfld_qpf%bmap(im*jm))

        if (iret1 .eq. 0) then
        gfld_qpf%bmap=bmap_f
        else
        write(0,*) 'choked on QPF get: ', iret1
        endif

! ---------------------------------


	if (reset_flag .eq. 1) then
	write(0,*) 'just later value'

	do NPT=1,IM*JM

        if (is_hrrr .eq. 0) then
	sprecip(NPT)=s_later(NPT)
        endif

	asnowprecip(NPT)=asnow_later(NPT)

        if (is_hrrr .eq. 1) then
	fzprecip(NPT)=fz_later(NPT)
        endif

	enddo

	else

	write(0,*) 'take normal difference ', IM*JM


	do NPT=1,IM*JM

!tmp        if (  gfld_qpf%bmap(NPT) ) then ! make sure is a valid point
!tried        if (  .not. bmap_f(NPT) ) then ! make sure is a valid point

        if (is_hrrr .eq. 0) then
	sprecip(NPT)=max(s_later(NPT)-s_earlier(NPT),0.0)
        endif

        if (is_hrrr .eq. 1) then
	fzprecip(NPT)=max(fz_later(NPT)-fz_earlier(NPT),0.0)
        endif

	asnowprecip(NPT)=asnow_later(NPT)-asnow_earlier(NPT)

!        if (asnowprecip(NPT) .ne. 0.) then
!         write(0,*) 'NPT,asnow_later(NPT),asnow_earlier(NPT),diff: ',
!     &           NPT,asnow_later(NPT),asnow_earlier(NPT),
!     &           asnow_later(NPT)-asnow_earlier(NPT)
!         endif

!        else
!         write(0,*) 'skipping due to bmap_f at NPT: ', NPT
         
!tmp        endif

	enddo


         if (interv .eq. 1) then
                 write(0,*) '1 h values'
                 if (is_hrrr .eq. 0) then
        write(0,*) 'min,max of sprecip: ', minval(sprecip), 
     &         maxval(sprecip)
                 endif
                 if (is_hrrr .eq. 1) then
        write(0,*) 'min,max of fzprecip: ', minval(fzprecip), 
     &         maxval(fzprecip)
                 endif
        write(0,*) 'min,max of asnowprecip: ', minval(asnowprecip), 
     &         maxval(asnowprecip)
        write(0,*) 'min,max of asnow_later: ', minval(asnow_later), 
     &         maxval(asnow_later)
        write(0,*) 'min,max of asnow_earlier: ', minval(asnow_earlier), 
     &         maxval(asnow_earlier)
        endif


	endif

	deallocate(s_later)
	deallocate(s_earlier)
	deallocate(asnow_later)
	deallocate(asnow_earlier)

	write(0,*) 'shape(gfld_qpf%ipdtmpl): ', shape(gfld_qpf%ipdtmpl)

        gfld_qpf%ipdtmpl(22)=1
        gfld_qpf%ipdtmpl(27)=interv

        write(0,*) 'interval specified in 27: ', interv

! -------------------------------------------


       gfld_qpf%ipdtmpl(1)=1
       gfld_qpf%ipdtmpl(2)=13
       gfld_qpf%ipdtmpl(3)=2
       gfld_qpf%ipdtmpl(4)=0
!      gfld_qpf%ipdtmpl(5)=84
       gfld_qpf%ipdtmpl(5)=134 !for rrfs, J. Du
       gfld_qpf%ipdtmpl(6)=0 ! hours cutoff
       gfld_qpf%ipdtmpl(7)=0 ! minutes cutoff
       gfld_qpf%ipdtmpl(8)=1 ! units of hours
       gfld_qpf%ipdtmpl(9)=ihrs1 ! earlier forecast time of period?
       gfld_qpf%ipdtmpl(10)=1 ! sfc
       gfld_qpf%ipdtmpl(11)=0 ! sfc
       gfld_qpf%ipdtmpl(12)=0 ! sfc
       gfld_qpf%ipdtmpl(13)=255 ! sfc
       gfld_qpf%ipdtmpl(14)=0 ! sfc
       gfld_qpf%ipdtmpl(15)=0 ! sfc

!!! need to figure out how to do this end of period date stuff right

       rinc=0.
       rinc(2)=float(ihrs1+interv)

       write(0,*) 'rinc: ', rinc
       idat=0
       idat(1)=gfld_qpf%idsect(6)
       idat(2)=gfld_qpf%idsect(7)
       idat(3)=gfld_qpf%idsect(8) 
       idat(5)=gfld_qpf%idsect(9)


       call w3movdat(rinc,idat,jdat)

       gfld_qpf%ipdtmpl(16)=jdat(1)
       gfld_qpf%ipdtmpl(17)=jdat(2)
       gfld_qpf%ipdtmpl(18)=jdat(3)
       gfld_qpf%ipdtmpl(19)=jdat(5)
       gfld_qpf%ipdtmpl(20)=0
       gfld_qpf%ipdtmpl(21)=0

        gfld_qpf%ipdtnum=8
        gfld_qpf%ipdtmpl(22)=1
        gfld_qpf%ipdtmpl(23)=0
        gfld_qpf%ipdtmpl(24)=1 ! accum?
        gfld_qpf%ipdtmpl(25)=2 ! fcst hour increments
        gfld_qpf%ipdtmpl(26)=1 ! hours?
        gfld_qpf%ipdtmpl(27)=interv
        gfld_qpf%ipdtmpl(28)=255
        gfld_qpf%ipdtmpl(29)=0

        gfld_qpf%fld=sprecip

	write(0,*) 'use gfld_qpf%idrtmpl(1:10): ', gfld_qpf%idrtmpl(1:10)

!! use GET_BITS to compute nbits?

      IBM=0
      IBITM = 0
      gfld_qpf%idrtmpl(3)=5.0
      SGDS  = gfld%idrtmpl(3)

!     set bitmap
      DO N=1,IM*JM
        IF( gfld_qpf%bmap(N) ) THEN
             ibmap(N) = 1
             ibitm = ibitm+1
        ELSE
             ibmap(N) = 0
        ENDIF
      ENDDO

!     set bitmap
      IF (IBITM.EQ.IM*JM) THEN
        IBM = 0
      ELSE
        IBM = 1
      ENDIF
      call GET_BITS(IBM,SGDS,IM*JM,ibmap,gfld_qpf%fld,
     &                ISCALE,GROUND,GMIN,GMAX,NBIT)

      write(0,*) 'returned NBIT for WEASD as: ', NBIT
!      write(0,*) 'returned ISCALE as: ', ISCALE
!      write(0,*) 'GMAX: ', GMAX


        gfld_qpf%idrtmpl(4)=NBIT

	write(0,*) 'use gfld%idrtmpl(1:10): ', gfld%idrtmpl(1:10)

        if (is_hrrr .eq. 0) then
	call putgb2(13,GFLD_QPF,IRET)
        write(0,*) 'IRET from putgb2 for sprecip ', IRET
        endif

! add asnow piece

       gfld_qpf%ipdtmpl(2)=29
       gfld_qpf%fld=asnowprecip

!! use GET_BITS to compute nbits?

      IBM=0
      IBITM = 0
      gfld_qpf%idrtmpl(3)=5.0
      SGDS  = gfld%idrtmpl(3)

!     set bitmap
      DO N=1,IM*JM
        IF( gfld_qpf%bmap(N) ) THEN
             ibmap(N) = 1
             ibitm = ibitm+1
        ELSE
             ibmap(N) = 0
        ENDIF
      ENDDO

!     set bitmap
      IF (IBITM.EQ.IM*JM) THEN
        IBM = 0
      ELSE
        IBM = 1
      ENDIF
      call GET_BITS(IBM,SGDS,IM*JM,ibmap,gfld_qpf%fld,
     &                ISCALE,GROUND,GMIN,GMAX,NBIT)

      write(0,*) 'returned NBIT for ASNOW as: ', NBIT
!      write(0,*) 'returned ISCALE as: ', ISCALE
!      write(0,*) 'GMAX: ', GMAX


        gfld_qpf%idrtmpl(4)=NBIT

	write(0,*) 'use gfld%idrtmpl(1:10): ', gfld%idrtmpl(1:10)

	call putgb2(13,GFLD_QPF,IRET)
        write(0,*) 'IRET from putgb2 for asnowprecip ', IRET

	write(0,*) 'extremes of snow,asnow: ', 
     +		maxval(sprecip),maxval(asnowprecip)

! add frzr piece
        if (is_hrrr .eq. 1) then

       gfld_qpf%ipdtmpl(2)=225
       gfld_qpf%fld=fzprecip

!! use GET_BITS to compute nbits?

      IBM=0
      IBITM = 0
      gfld_qpf%idrtmpl(2)=0.0
      gfld_qpf%idrtmpl(3)=5.0
      SGDS  = gfld%idrtmpl(3)

!     set bitmap
      DO N=1,IM*JM
        IF( gfld_qpf%bmap(N) ) THEN
             ibmap(N) = 1
             ibitm = ibitm+1
        ELSE
             ibmap(N) = 0
        ENDIF
      ENDDO

!     set bitmap
      IF (IBITM.EQ.IM*JM) THEN
        IBM = 0
      ELSE
        IBM = 1
      ENDIF
      call GET_BITS(IBM,SGDS,IM*JM,ibmap,gfld_qpf%fld,
     &                ISCALE,GROUND,GMIN,GMAX,NBIT)

      write(0,*) 'returned NBIT for FZ as: ', NBIT
!      write(0,*) 'returned ISCALE as: ', ISCALE
!      write(0,*) 'GMAX: ', GMAX


        gfld_qpf%idrtmpl(4)=NBIT

	write(0,*) 'use gfld%idrtmpl(1:10): ', gfld%idrtmpl(1:10)

	call putgb2(13,GFLD_QPF,IRET)
        write(0,*) 'IRET from putgb2 for fzprecip ', IRET
	write(0,*) 'extremes of fzprecip ', 
     +		maxval(fzprecip)

        endif

        call baclose(13,IRET)

	end subroutine calc_pdiff
