C  raw data
       use grib_mod
       real,allocatable,dimension(:,:) :: dp3,sn3,asn3,fz3 !jf,4        
       real,allocatable,dimension(:) :: dp6,dp12,dp24 !jf         
       real,allocatable,dimension(:) :: sn6,sn12,sn24 !jf         
       real,allocatable,dimension(:) :: fz6,fz12,fz24 !jf         
       real,allocatable,dimension(:) :: asn6,asn12,asn24 !jf         
 
       integer iyr,imon,idy,ihr
       character*50 gdss(400)
       integer IENS, GRIBID, kgdss(200),lengds,im,jm,km,jf
       character*40 filehead,filename(8), output, outdone
       
       integer ff
       logical do_old, do_hrrr, do_hrrr_pre, do_fv3_pre, skip_1h
       character*3 fhr(20),fog
       character*5 domain
       integer iunit,ounit, pdt9_orig, acclength
       type(gribfield) :: gfld,gfld_save,gfld_save_snow,
     +                    gfld_save_frzr,gfld_save_asn

       data (fhr(i),i=1,20)
     + /'f03','f06','f09','f12','f15','f18','f21',
     +  'f24','f27','f30','f33','f36','f39','f42','f45','f48',
     +  'f51','f54','f57','f60'/


 !  ${RUN}.m${m}.t${cyc}z. $ff 
! do_old: .false.
! do_hrrr:  .false. 
! do_hrrr_pre: .false.
! do_fv3_pre:  .false. 
! skip_1h: .false. 
! acc_length: 1
! domain: conus
! fog: yes

       read (*,*) filehead, ff, do_old, do_hrrr, do_hrrr_pre, 
     +            do_fv3_pre,skip_1h,acclength,domain,fog

	if (domain(1:5) .eq. 'conus') then
         GRIBID=255            !namnest grid
        else if (domain(1:2) .eq. 'ak') then
         GRIBID=999            !AK hiresw grid
        else if (domain(1:2) .eq. 'hi') then
         GRIBID=998            !HI hiresw grid
        else if (domain(1:2) .eq. 'pr') then
         GRIBID=997            !PR hiresw grid
        endif


	write(0,*) 'start prcip code'
	write(0,*) 'domain: ', domain
	write(0,*) 'GRIBID: ', GRIBID

	write(0,*) 'do_hrrr: ', do_hrrr
	write(0,*) 'do_hrrr_pre: ', do_hrrr_pre
	write(0,*) 'do_fv3_pre: ', do_fv3_pre
	write(0,*) 'acclength: ', acclength
	write(0,*) 'fog: ', fog

cc     RAP has one-hour accumu precip, so only one file is used
cc     NAM has no one-hour accumu precip, so two files are needed

       jf=0

       if(GRIBID.eq.255) then   !For NARRE 13km RAP grid#130
         im=1799
         jm=1059
         jf=im*jm
       elseif (GRIBID.eq.999) then ! AK  grid
         im=1649
         jm=1105
         jf=im*jm
       elseif (GRIBID.eq.998) then ! HI grid
         im=321
         jm=225
         jf=im*jm
       elseif (GRIBID.eq.997) then ! PR grid
         im=544
         jm=310
         jf=im*jm
       else
         call makgds(GRIBID, kgdss, gdss, lengds, ier)
         im=kgdss(2)
         jm=kgdss(3)
         jf=im*jm
       end if

       write(*,*) 'jf=',jf

	IF (acclength .eq. 3) THEN

       allocate(dp3(jf,8))
       allocate(dp6(jf))
       allocate(dp12(jf))
       allocate(dp24(jf))

       allocate(sn3(jf,8))
       allocate(sn6(jf))
       allocate(sn12(jf))
       allocate(sn24(jf))

       allocate(asn3(jf,8))
       allocate(asn6(jf))
       allocate(asn12(jf))
       allocate(asn24(jf))

       allocate(fz3(jf,8))
       allocate(fz6(jf))
       allocate(fz12(jf))
       allocate(fz24(jf))

       if (ff.ge.24) then
         nfile=8
       else if (ff.lt.24.and.ff.ge.12) then
         nfile=4
       else if (ff.lt.12.and.ff.ge.6) then
         nfile=2
       else
         nfile=1 
       end if
 
       nff=ff/3
       do 1000 nf=1,nfile
        
        filename(nf)=filehead(1:14)//fhr(nff)

        iunit=20+nf

        jpdtn=8    !APCP's Product Template# is  4.8 

        call baopenr(iunit,filename(nf),ierr)
        write(*,*) 'open ', filename(nf), 'ierr=',ierr

	if (ierr .eq. 0) then

!       APCP

        jpd1=1
        jpd2=8
        jpd27=3 !3 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
	write(0,*) 'populate nf of nfile: ', nf, nfile
         dp3(:,nf)=gfld%fld(:)
         if (nf.eq.1) then 
           gfld_save=gfld
	write(0,*) 'gfld_save%ipdtmpl(9) when saved: ', gfld_save%ipdtmpl(9)
	pdt9_orig=gfld_save%ipdtmpl(9)
           do i=1,gfld_save%ipdtlen
            write(*,*) i, gfld_save%ipdtmpl(i)
           end do
         end if
        else
         write(*,*) '3h readGB2 apcp error=',ie
        end if

!       FRZR

        jpd1=1
        jpd2=225
        jpd27=3 !3 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! FRZR

        if (ie.eq.0) then
         fz3(:,nf)=gfld%fld(:)
         if (nf.eq.1) then 
           gfld_save_frzr=gfld
         end if
        else
         write(*,*) '3h readGB2 frzr error=',ie
        end if

!       ASNOW

        jpd1=1
        jpd2=29
        jpd27=3 !3 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! ASNOW

        if (ie.eq.0) then
         asn3(:,nf)=gfld%fld(:)
         if (nf.eq.1) then 
           gfld_save_asn=gfld
         end if
        else
         write(*,*) '3h readGB2 asnow error=',ie
        end if

!       SNOW

        jpd1=1
        jpd2=13
        jpd27=3 !3 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !
        if (ie.eq.0) then
	write(0,*) 'populate nf of nfile: ', nf, nfile
         sn3(:,nf)=gfld%fld(:)
         if (nf.eq.1) then 
           gfld_save_snow=gfld
         end if
        else
         write(*,*) '3h readGB2 snow error=',ie
        end if

	endif

        call baclose(iunit,ierr)
        write(*,*) 'close ', filename(nf), 'ierr=',ierr
        nff=nff-1

1000  continue

        dp6=0.0
        dp24=0.0
        dp12=0.0
        sn6=0.0
        sn24=0.0
        sn12=0.0
        asn6=0.0
        asn24=0.0
        asn12=0.0
        fz6=0.0
        fz24=0.0
        fz12=0.0


       if (ff.ge.24) then
         dp6(:)=dp3(:,1)+dp3(:,2)
         dp12(:)=dp6(:)+dp3(:,3)+dp3(:,4)
         dp24(:)=dp12(:)+dp3(:,5)+dp3(:,6)+dp3(:,7)+dp3(:,8)
         sn6(:)=sn3(:,1)+sn3(:,2)
         sn12(:)=sn6(:)+sn3(:,3)+sn3(:,4)
         sn24(:)=sn12(:)+sn3(:,5)+sn3(:,6)+sn3(:,7)+sn3(:,8)
         asn6(:)=asn3(:,1)+asn3(:,2)
         asn12(:)=asn6(:)+asn3(:,3)+asn3(:,4)
         asn24(:)=asn12(:)+asn3(:,5)+asn3(:,6)+asn3(:,7)+asn3(:,8)
         fz6(:)=fz3(:,1)+fz3(:,2)
         fz12(:)=fz6(:)+fz3(:,3)+fz3(:,4)
         fz24(:)=fz12(:)+fz3(:,5)+fz3(:,6)+fz3(:,7)+fz3(:,8)
       else if (ff.lt.24.and.ff.ge.12) then
         dp6(:)=dp3(:,1)+dp3(:,2)
         dp12(:)=dp6(:)+dp3(:,3)+dp3(:,4)
         sn6(:)=sn3(:,1)+sn3(:,2)
         sn12(:)=sn6(:)+sn3(:,3)+sn3(:,4)
         asn6(:)=asn3(:,1)+asn3(:,2)
         asn12(:)=asn6(:)+asn3(:,3)+asn3(:,4)
         fz6(:)=fz3(:,1)+fz3(:,2)
         fz12(:)=fz6(:)+fz3(:,3)+fz3(:,4)
       else if (ff.lt.12.and.ff.ge.6) then
	write(0,*) 'adding to create dp6'
	write(0,*) 'maxvals of dp3 inputs: ', 
     &          maxval(dp3(:,1)),maxval(dp3(:,2))
         dp6(:)=dp3(:,1)+dp3(:,2)
         sn6(:)=sn3(:,1)+sn3(:,2)
         asn6(:)=asn3(:,1)+asn3(:,2)
         fz6(:)=fz3(:,1)+fz3(:,2)
       end if
            
!       do i=382461,382470
!        write(*,'(i10,11f8.2)')i,(dp3(i,k),k=1,8),
!     +         dp6(i),dp12(i),dp24(i)                                   
!       end do


cccccc  Then call putgb2 to store the calculated data into a grib2 file
c
c      data structure gfld is re-used for pack data since all are same
c      only gfld%fld and gfld%ipdtmpl(27) are different

Cmp   believe gfld%ipdtmpl(9) matters as well
c

c      If a field not in a GRIB2 file, getGB2 output gfld will be crashed. 
c      so use previously saved gfld_save 

       nff=ff/3      
       output='prcip3h'//filehead(5:14)//fhr(nff)
       outdone='prcipdone'//filehead(9:14)//fhr(nff)

        ounit=50+nff
        call baopen(ounit,output,ierr)

	write(0,*) 'setting gfld to gfld_save'

          if(ff.ge.24) then

	     gfld=gfld_save
             gfld%fld(:)=dp3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_snow
             gfld%fld(:)=sn3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_frzr
             gfld%fld(:)=fz3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_asn
             gfld%fld(:)=asn3(:,1)
             gfld%ipdtmpl(27)=3
	     write(0,*) 'ASN 3 h maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save
             gfld%fld(:)=dp6(:)
             gfld%ipdtmpl(27)=6
             gfld%ipdtmpl(9)=-3 + pdt9_orig
	write(0,*) 'dp6 gfld%ipdtmpl(9) now: ', gfld%ipdtmpl(9)
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_snow
             gfld%ipdtmpl(27)=6
             gfld%ipdtmpl(9)=-3 + pdt9_orig
             gfld%fld(:)=sn6(:)
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_asn
             gfld%ipdtmpl(27)=6
             gfld%ipdtmpl(9)=-3 + pdt9_orig
             gfld%fld(:)=asn6(:)
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_frzr
             gfld%ipdtmpl(27)=6
             gfld%ipdtmpl(9)=-3 + pdt9_orig
             gfld%fld(:)=fz6(:)
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save
             gfld%fld(:)=dp12(:)
             gfld%ipdtmpl(27)=12
             gfld%ipdtmpl(9)=-9 + pdt9_orig
	write(0,*) 'gfld%ipdtmpl(9) for dp12 now: ', gfld%ipdtmpl(9)
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_snow
             gfld%fld(:)=sn12(:)
             gfld%ipdtmpl(27)=12
             gfld%ipdtmpl(9)=-9 + pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_asn
             gfld%fld(:)=asn12(:)
             gfld%ipdtmpl(27)=12
             gfld%ipdtmpl(9)=-9 + pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_frzr
             gfld%fld(:)=fz12(:)
             gfld%ipdtmpl(27)=12
             gfld%ipdtmpl(9)=-9 + pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save
             gfld%fld(:)=dp24(:)
             gfld%ipdtmpl(27)=24
             gfld%ipdtmpl(9)=-21+pdt9_orig
	write(0,*) 'dp24 gfld%ipdtmpl(9) now: ', gfld%ipdtmpl(9)
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_snow
             gfld%fld(:)=sn24(:)
             gfld%ipdtmpl(27)=24
             gfld%ipdtmpl(9)=-21+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_asn
             gfld%fld(:)=asn24(:)
             gfld%ipdtmpl(27)=24
             gfld%ipdtmpl(9)=-21+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_frzr
             gfld%fld(:)=fz24(:)
             gfld%ipdtmpl(27)=24
             gfld%ipdtmpl(9)=-21+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

          else if (ff.lt.24.and.ff.ge.12) then

	     gfld=gfld_save
             gfld%fld(:)=dp3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_snow
             gfld%fld(:)=sn3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_frzr
             gfld%fld(:)=fz3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_asn
             gfld%fld(:)=asn3(:,1)
             gfld%ipdtmpl(27)=3
	     write(0,*) 'ASN 3 h maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save
             gfld%fld(:)=dp6(:)
             gfld%ipdtmpl(27)=6
             gfld%ipdtmpl(9)=-3+pdt9_orig
	write(0,*) 'dp6 gfld%ipdtmpl(9) now: ', gfld%ipdtmpl(9)
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_snow
             gfld%fld(:)=sn6(:)
             gfld%ipdtmpl(27)=6
             gfld%ipdtmpl(9)=-3+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_asn
             gfld%fld(:)=asn6(:)
             gfld%ipdtmpl(27)=6
             gfld%ipdtmpl(9)=-3+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_frzr
             gfld%fld(:)=fz6(:)
             gfld%ipdtmpl(27)=6
             gfld%ipdtmpl(9)=-3+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save
             gfld%fld(:)=dp12(:)
             gfld%ipdtmpl(27)=12
             gfld%ipdtmpl(9)=-9+pdt9_orig
	write(0,*) 'dp12(b) gfld%ipdtmpl(9) now: ', gfld%ipdtmpl(9)
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_snow
             gfld%fld(:)=sn12(:)
             gfld%ipdtmpl(27)=12
             gfld%ipdtmpl(9)=-9+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_asn
             gfld%fld(:)=asn12(:)
             gfld%ipdtmpl(27)=12
             gfld%ipdtmpl(9)=-9+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_frzr
             gfld%fld(:)=fz12(:)
             gfld%ipdtmpl(27)=12
             gfld%ipdtmpl(9)=-9+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

          else if (ff.lt.12.and.ff.ge.6) then

	     gfld=gfld_save
             gfld%fld(:)=dp3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_snow
             gfld%fld(:)=sn3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_frzr
             gfld%fld(:)=fz3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_asn
             gfld%fld(:)=asn3(:,1)
             gfld%ipdtmpl(27)=3
	     write(0,*) 'ASN 3 h maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save
             gfld%fld(:)=dp6(:)
             gfld%ipdtmpl(27)=6       
             gfld%ipdtmpl(9)=-3+pdt9_orig
	write(0,*) 'dp6(b) gfld%ipdtmpl(9) now: ', gfld%ipdtmpl(9)
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_snow
             gfld%fld(:)=sn6(:)
             gfld%ipdtmpl(27)=6       
             gfld%ipdtmpl(9)=-3+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_asn
             gfld%fld(:)=asn6(:)
             gfld%ipdtmpl(27)=6       
             gfld%ipdtmpl(9)=-3+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_frzr
             gfld%fld(:)=fz6(:)
             gfld%ipdtmpl(27)=6       
             gfld%ipdtmpl(9)=-3+pdt9_orig
             call putgb2_wrap(ounit,gfld,ierr)

           else

	     gfld=gfld_save
             gfld%fld(:)=dp3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_snow
             gfld%fld(:)=sn3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_frzr
             gfld%fld(:)=fz3(:,1)
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

	     gfld=gfld_save_asn
             gfld%fld(:)=asn3(:,1)
	     write(0,*) 'ASN 3 h maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
             gfld%ipdtmpl(27)=3
             call putgb2_wrap(ounit,gfld,ierr)

          end if
    
        write(0,*) 'Pack 3-24h APCP done for fhr',ff

! write "done" file

	write(0,*) 'acclength=',acclength
        call baclose(ounit,ierr) 
	write(0,*) 'calling just_hrly(a) with jf: ', jf
        call just_hrly(filehead, ff, jf, do_old, skip_1h)

        ELSE !every hour with acclength=1

        if (fog .eq. 'yes' ) then
         call get_temp(filehead, ff, jf)
        endif

	if (.not. do_hrrr_pre .and. .not. do_fv3_pre ) then
        write(0,*) 'I am here doing hourly at fcst hour ',ff
	write(0,*) 'acclength=',acclength
	write(0,*) 'calling just_hrly(b) with jf: ', jf
        call just_hrly(filehead, ff, jf, do_old, skip_1h )
        endif

	if (do_hrrr) then
	write(0,*) 'calling just_hrrr_3hrly'
        call just_hrrr_3hrly(filehead, ff, jf)
        endif

	if (do_hrrr_pre) then
	write(0,*) 'calling just_hrrr_3hrly_pre'
        call just_hrrr_3hrly_pre(filehead, ff, jf)
        endif

	if (do_fv3_pre) then
	write(0,*) 'calling just_fv3_3hrly_pre'
        call just_fv3_3hrly_pre(filehead, ff, jf)
        endif

	ENDIF

      stop
      end

! -----------------------------

	subroutine just_hrly(filehead, ff, jf, do_old, skip_1h)

C  raw data
       use grib_mod
       real,allocatable,dimension(:,:) ::  dphold !jf,4        
       real,allocatable,dimension(:)   ::  dp1
       real,allocatable,dimension(:,:) ::  snhold !jf,4        
       real,allocatable,dimension(:)   ::  sn1
       real,allocatable,dimension(:,:) ::  fzhold !jf,4        
       real,allocatable,dimension(:)   ::  fz1
       real,allocatable,dimension(:,:) ::  asnhold !jf,4        
       real,allocatable,dimension(:)   ::  asn1

       integer iyr,imon,idy,ihr
       character*50 gdss(400)
       integer IENS, GRIBID, kgdss(200), lengds,im,jm,km,jf
       character*40 filehead,filename(8),output,outdone
       
       integer ff
       logical do_old, skip_1h
       character*3 fhr(60)
       integer iunit,ounit,pdt9_orig
       type(gribfield) :: gfld,gfld_save_1h,gfld_2h
       type(gribfield) :: gfld_snow,     gfld_save_1h_snow ! ,gfld_2h
       type(gribfield) :: gfld_frzr,     gfld_save_1h_frzr ! ,gfld_2h
       type(gribfield) :: gfld_asn,     gfld_save_1h_asn ! ,gfld_2h

       data (fhr(i),i=1,60)
     + /'f01','f02','f03','f04','f05','f06','f07','f08','f09',
     + 'f10','f11','f12','f13','f14','f15','f16','f17','f18',
     + 'f19','f20','f21','f22','f23','f24','f25','f26','f27',
     + 'f28','f29','f30','f31','f32','f33','f34','f35','f36',
     + 'f37','f38','f39','f40','f41','f42','f43','f44','f45',
     + 'f46','f47','f48','f49','f50','f51','f52','f53','f54',
     + 'f55','f56','f57','f58','f59','f60'/
 
!	write(0,*) 'know that filehead, ff, do_old, jf', 
!     *    trim(filehead), ff, do_old, jf


!	write(0,*) 'start prcip code(b)'

!       write(0,*) 'jf=',jf

       allocate(dphold(jf,3))
       allocate(dp1(jf))
       allocate(snhold(jf,3))
       allocate(sn1(jf))
       allocate(fzhold(jf,3))
       allocate(fz1(jf))
       allocate(asnhold(jf,3))
       allocate(asn1(jf))

!! these numbers need to change for hourly

	if (do_old .and. ff .gt. 1) then
         nfile=2
        else
         nfile=1
        endif
	write(0,*) 'do_old, nfile: ', do_old, nfile
       nff=ff/1

        dphold=0.

       do 1000 nf=1,nfile
        
        filename(nf)=filehead(1:14)//fhr(nff)
	write(0,*) 'nff, ff, nf, trim(filename(nf)): ', nff,ff,
     &                nf,trim(filename(nf))

        iunit=20+nf

c default is set to APCP
        jpdtn=8    !APCP's Product Template# is  4.8 

        call baopenr(iunit,filename(nf),ierr)
        write(0,*) 'open ', iunit, filename(nf), 'ierr=',ierr

	if (ierr .eq. 0) then


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (mod(nff,3) .eq. 0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	write(0,*) 'START MOD=0 BLOCK'

!  if nfile = 2, read current and previous

!  get 1 h accumulation if possible at 3 hourly time

        jpd1=1
        jpd2=8
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP

        if (ie.eq.0) then

!	write(0,*) 'populate nf (mod=0) of nfile: ', nf, nfile
         dp1(:)=gfld%fld(:)
         if (nf.eq.1) then 
           gfld_save_1h=gfld
!  	   write(0,*) 'gfld_save(9) when saved: ', gfld_save_1h%ipdtmpl(9)
           do i=1,gfld_save_1h%ipdtlen
!            write(0,*) 'MOD=0 ', i, gfld_save_1h%ipdtmpl(i)
           end do
         end if

        else

! get 3 h accumulation that will be used in MOD=2 block to get the 1 h total

        jpd1=1
        jpd2=8
        jpd27=3 !3 hr accumulation
!	write(0,*) 'seek 3 h accum'
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP

        if (ie.eq.0) then
!	write(0,*) 'populate 3 h accum  nf of nfile: ', nf, nfile
        dphold(:,3)=gfld%fld(:)
	endif

       endif  ! ie=0 for 1 h accum


!    SNOW

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP

        if (ie.eq.0) then

!  	 write(0,*) 'populate SNOW nf (mod=0) of nfile: ', nf, nfile
         sn1(:)=gfld%fld(:)
         if (nf.eq.1) then 
           gfld_save_1h_snow=gfld
         end if

        else

        jpd1=1
        jpd2=13
        jpd27=3 !3 hr accumulation
!	write(0,*) 'seek 3 h SNOW accum'
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP

        if (ie.eq.0) then
!	write(0,*) 'populate 3 h SNOW accum  nf of nfile: ', nf, nfile
        snhold(:,3)=gfld%fld(:)
	endif

       endif  ! ie=0 for 1 h accum


!       FRZR

        jpd1=1
        jpd2=225
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! FRZR

        if (ie.eq.0) then

!  	 write(0,*) 'populate FRZR nf (mod=0) of nfile: ', nf, nfile
         fz1(:)=gfld%fld(:)
         if (nf.eq.1) then 
           gfld_save_1h_frzr=gfld
         end if

        else

        jpd1=1
        jpd2=225
        jpd27=3 !3 hr accumulation
!	write(0,*) 'seek 3 h FRZR accum'
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! FRZR

        if (ie.eq.0) then
!	write(0,*) 'populate 3 h FRZR accum  nf of nfile: ', nf, nfile
        fzhold(:,3)=gfld%fld(:)
	endif

        endif ! ie=0 for 1 h accum

!       ASNOW

        jpd1=1
        jpd2=29
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! ASNOW

        if (ie.eq.0) then

!  	 write(0,*) 'populate ASNOW nf (mod=0) of nfile: ', nf, nfile
         asn1(:)=gfld%fld(:)
         if (nf.eq.1) then 
           gfld_save_1h_asn=gfld
         end if

        else

        jpd1=1
        jpd2=29
        jpd27=3 !3 hr accumulation
!	write(0,*) 'seek 3 h ASNOW accum'
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! ASNOW 

        if (ie.eq.0) then
!	write(0,*) 'populate 3 h ASNOW accum  nf of nfile: ', nf, nfile
        asnhold(:,3)=gfld%fld(:)
	endif

        endif ! ie=0 for 1 h accum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        elseif (mod(nff,3) .eq. 1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! all should have 1 h total

!	write(0,*) 'START MOD=1 BLOCK'


        jpd1=1
        jpd2=8
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP

        if (ie.eq.0) then

	if (nf .eq. 1) then
!	write(0,*) 'mod=1, populate nf of nfile: ', nf, nfile
         dp1(:)=gfld%fld(:)
           gfld_save_1h=gfld

           do i=1,gfld_save_1h%ipdtlen
!            write(0,*) 'MOD=1 ', i, gfld_save_1h%ipdtmpl(i)
           enddo

        else

!	write(0,*) 'in here when 1 h old from 2 h block'
          dp1(:)=dphold(:,2)-gfld%fld(:)
	   gfld%fld(:)=dp1(:)
           gfld%ipdtmpl(9)=gfld%ipdtmpl(9)+1
           gfld%ipdtmpl(19)=gfld%ipdtmpl(19)+1
           gfld_save_1h=gfld

        endif  ! nf=1

         if (nf.eq.1) then 
           gfld_save_1h=gfld
  	   write(0,*) 'gfld_save_1h(9) when saved: ', gfld_save_1h%ipdtmpl(9)
!           do i=1,gfld_save_1h%ipdtlen
!            write(0,*) i, gfld_save_1h%ipdtmpl(i)
!           end do
         end if

       endif ! ie=0


! SNOW

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP

        if (ie.eq.0) then

	if (nf .eq. 1) then
!	write(0,*) 'mod=1, SNOW populate nf of nfile: ', nf, nfile
         sn1(:)=gfld%fld(:)
           gfld_save_1h_snow=gfld

!           do i=1,gfld_save_1h%ipdtlen
!            write(0,*) 'MOD=1 ', i, gfld_save_1h%ipdtmpl(i)
!           enddo

        else

!	write(0,*) 'in here when 1 h old from 2 h block'
          sn1(:)=snhold(:,2)-gfld%fld(:)
	   gfld%fld(:)=sn1(:)
           gfld%ipdtmpl(9)=gfld%ipdtmpl(9)+1
           gfld%ipdtmpl(19)=gfld%ipdtmpl(19)+1
           gfld_save_1h_snow=gfld

        endif ! nf=1

         if (nf.eq.1) then 
           gfld_save_1h_snow=gfld
         end if


        endif ! ie=0

! FRZR

        jpd1=1
        jpd2=225
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !

        if (ie.eq.0) then

	if (nf .eq. 1) then
!	write(0,*) 'mod=1, FRZR populate nf of nfile: ', nf, nfile
         fz1(:)=gfld%fld(:)
           gfld_save_1h_frzr=gfld

!           do i=1,gfld_save_1h%ipdtlen
!            write(0,*) 'MOD=1 ', i, gfld_save_1h%ipdtmpl(i)
!           enddo

        else

!	write(0,*) 'in here when 1 h old from 2 h block'
          fz1(:)=fzhold(:,2)-gfld%fld(:)
	   gfld%fld(:)=fz1(:)
           gfld%ipdtmpl(9)=gfld%ipdtmpl(9)+1
           gfld%ipdtmpl(19)=gfld%ipdtmpl(19)+1
           gfld_save_1h_frzr=gfld

        endif ! nf=1

         if (nf.eq.1) then 
           gfld_save_1h_frzr=gfld
         end if


        endif ! ie=0

! ASNOW

        jpd1=1
        jpd2=29
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !

        if (ie.eq.0) then

	if (nf .eq. 1) then
!	write(0,*) 'mod=1, ASNOW populate nf of nfile: ', nf, nfile
         asn1(:)=gfld%fld(:)
           gfld_save_1h_asn=gfld
        else

!	write(0,*) 'in here when 1 h old from 2 h block'
          asn1(:)=asnhold(:,2)-gfld%fld(:)
	   gfld%fld(:)=asn1(:)
           gfld%ipdtmpl(9)=gfld%ipdtmpl(9)+1
           gfld%ipdtmpl(19)=gfld%ipdtmpl(19)+1
           gfld_save_1h_asn=gfld

        endif ! nf=1

         if (nf.eq.1) then 
           gfld_save_1h_asn=gfld
         end if


        endif ! ie=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	elseif (mod(nff,3) .eq. 2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	write(0,*) 'START MOD=2 BLOCK'

!  if nfile = 2, read current and previous

         jpdtn=8    !APCP's Product Template# is  4.8 
         jpd1=1
         jpd2=8
         jpd27=1 !1 hr accumulation

         call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP

         if (ie.eq.0) then
! 	  write(0,*) 'mod=2, 1h accum  populate nf of nfile: ', nf, nfile
          dp1(:)=gfld%fld(:)
          dphold(:,1)=gfld%fld(:)
!          write(0,*) 'maxval(dphold(:,2)) ', maxval(dphold(:,2))

  	  if ( maxval(dphold(:,2)) .gt. 0) then 
! 	   write(0,*) 'inside maxval(dphold(:,2) '
           dp1(:)=dphold(:,2)-dphold(:,1)
	   write(0,*) 'maxval dpholds, dp1: ', maxval(dphold(:,1)), 
     +                 maxval(dphold(:,2)), maxval(dp1)
          endif

          if (nf.eq.1) then 
           gfld_save_1h=gfld
!  	   write(0,*) 'gfld_save_1h(9) when saved: ', gfld_save_1h%ipdtmpl(9)
           do i=1,gfld_save_1h%ipdtlen
!            write(0,*) 'MOD=2 ', i, gfld_save_1h%ipdtmpl(i)
           enddo
          end if

         else ! ie not = 0

          jpd1=1
          jpd2=8
          jpd27=2 !2 hr accumulation
          call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
          if (ie.eq.0) then
!	    write(0,*) 'populate mod=2, 2 h accum nf of nfile: ', nf, nfile
            dphold(:,2)=gfld%fld(:)
!	    write(0,*) 'maxval(dphold(:,2)) ', maxval(dphold(:,2))
	
	  if (nf .eq. 2) then
	   dp1(:)=dphold(:,3)-dphold(:,2)
!	   write(0,*) 'definined dp1 from difference'
!	   write(0,*) 'maxval(dp1): ', maxval(dp1)
	   gfld%fld(:)=dp1(:)
           gfld%ipdtmpl(9)=gfld%ipdtmpl(9)+2
           gfld%ipdtmpl(19)=gfld%ipdtmpl(19)+1
           gfld_save_1h=gfld
          endif ! nf=2

          endif  ! ie check
          endif  ! ie check

! SNOW 

         jpd1=1
         jpd2=13
         jpd27=1 !1 hr accumulation

         call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
         if (ie.eq.0) then
 	  write(0,*) 'mod=2, 1h accum SN populate nf of nfile: ', nf, nfile
          sn1(:)=gfld%fld(:)
          snhold(:,1)=gfld%fld(:)
	  write(0,*) 'maxval(snhold(:,2)) ', maxval(snhold(:,2))

!!! is this safe??

	if ( maxval(snhold(:,2)) .gt. 0 .and. 
     +                 maxval(snhold(:,2)) .lt. 1.e4) then 
	write(0,*) 'inside maxval(snhold(:,2) '
            sn1(:)=snhold(:,2)-snhold(:,1)
	write(0,*) 'maxval snholds, sn1: ', maxval(snhold(:,1)), 
     +                       maxval(snhold(:,2)), maxval(sn1)
        endif

         if (nf.eq.1) then 
           gfld_save_1h_snow=gfld
!           do i=1,gfld_save_1h%ipdtlen
!            write(0,*) 'MOD=2 ', i, gfld_save_1h%ipdtmpl(i)
!          enddo
         end if

         else ! ie not = 0

          jpd1=1
          jpd2=13
          jpd27=2 !2 hr accumulation
          call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
          if (ie.eq.0) then
	    write(0,*) 'populate SN mod=2, 2 h accum nf of nfile: ',nf,nfile
            snhold(:,2)=gfld%fld(:)
            write(0,*) 'maxval(snhold(:,2)) ', maxval(snhold(:,2))
	
	if (nf .eq. 2) then
	 sn1(:)=snhold(:,3)-snhold(:,2)
	 write(0,*) 'definined sn1 from difference'
	 write(0,*) 'maxval(sn1): ', maxval(sn1)
	 gfld%fld(:)=sn1(:)
         gfld%ipdtmpl(9)=gfld%ipdtmpl(9)+2
         gfld%ipdtmpl(19)=gfld%ipdtmpl(19)+1
         gfld_save_1h_snow=gfld
	endif

         endif ! ie=0
         endif ! ie=0

! FRZR

         jpd1=1
         jpd2=225
         jpd27=1 !1 hr accumulation

         call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! FRZR
         if (ie.eq.0) then
! 	  write(0,*) 'mod=2, 1h accum FZ populate nf of nfile: ', nf, nfile
          fz1(:)=gfld%fld(:)
          fzhold(:,1)=gfld%fld(:)
!	  write(0,*) 'maxval(fzhold(:,2)) ', maxval(fzhold(:,2))

!!! is this safe??

	if ( maxval(fzhold(:,2)) .gt. 0) then 
!	write(0,*) 'inside maxval(snhold(:,2) '
            fz1(:)=fzhold(:,2)-fzhold(:,1)
!	write(0,*) 'maxval fzholds, fz1: ', maxval(fzhold(:,1)), 
!     +                       maxval(fzhold(:,2)), maxval(fz1)
        endif

         if (nf.eq.1) then 
           gfld_save_1h_frzr=gfld
!           do i=1,gfld_save_1h%ipdtlen
!            write(0,*) 'MOD=2 ', i, gfld_save_1h%ipdtmpl(i)
!          enddo
         end if

         else ! ie not = 0

          jpd1=1
          jpd2=225
          jpd27=2 !2 hr accumulation
          call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! FRZR
          if (ie.eq.0) then
!	    write(0,*) 'populate FZ mod=2, 2 h accum nf of nfile: ',nf,nfile
            fzhold(:,2)=gfld%fld(:)
!            write(0,*) 'maxval(fzhold(:,2)) ', maxval(fzhold(:,2))
	
	if (nf .eq. 2) then
	 fz1(:)=fzhold(:,3)-fzhold(:,2)
!	 write(0,*) 'definined fz1 from difference'
!	 write(0,*) 'maxval(fz1): ', maxval(fz1)
	 gfld%fld(:)=fz1(:)
         gfld%ipdtmpl(9)=gfld%ipdtmpl(9)+2
         gfld%ipdtmpl(19)=gfld%ipdtmpl(19)+1
         gfld_save_1h_frzr=gfld
	endif

         endif ! ie=0
         endif ! ie=0

! ASNOW

         jpd1=1
         jpd2=29
         jpd27=1 !1 hr accumulation

         call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! ASNOW
         if (ie.eq.0) then
! 	  write(0,*) 'mod=2, 1h accum ASN populate nf of nfile: ', nf, nfile
          asn1(:)=gfld%fld(:)
          asnhold(:,1)=gfld%fld(:)

!!! is this safe??

	if ( maxval(asnhold(:,2)) .gt. 0) then 
!	write(0,*) 'inside maxval(snhold(:,2) '
            asn1(:)=asnhold(:,2)-asnhold(:,1)
!	write(0,*) 'maxval asnholds, asn1: ', maxval(asnhold(:,1)), 
!     +                       maxval(asnhold(:,2)), maxval(asn1)
        endif

         if (nf.eq.1) then 
           gfld_save_1h_asn=gfld
!           do i=1,gfld_save_1h%ipdtlen
!            write(0,*) 'MOD=2 ', i, gfld_save_1h%ipdtmpl(i)
!          enddo
         end if

         else ! ie not = 0

          jpd1=1
          jpd2=29
          jpd27=2 !2 hr accumulation
          call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! ASNOW
          if (ie.eq.0) then
!	    write(0,*) 'populate ASNOW mod=2, 2 h accum nf of nfile: ',nf,nfile
            asnhold(:,2)=gfld%fld(:)
!            write(0,*) 'maxval(fzhold(:,2)) ', maxval(fzhold(:,2))
	
	if (nf .eq. 2) then
	 asn1(:)=asnhold(:,3)-asnhold(:,2)
	 write(0,*) 'definined asn1 from difference'
!	 write(0,*) 'maxval(asn1): ', maxval(asn1)
	 gfld%fld(:)=asn1(:)
         gfld%ipdtmpl(9)=gfld%ipdtmpl(9)+2
         gfld%ipdtmpl(19)=gfld%ipdtmpl(19)+1
         gfld_save_1h_asn=gfld
	endif

         endif ! ie=0
         endif ! ie=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

	endif ! mod(nff,3)
        endif !  ierr=0

!	write(0,*) 'to bottom and baclose ' , iunit

        call baclose(iunit,ierr)
        write(*,*) 'close ', filename(nf), 'ierr=',ierr
        nff=nff-1

1000  continue


cccccc  Then call putgb2 to store the calculated data into a grib2 file
c
c      data structure gfld is re-used for pack data since all are same
c      only gfld%fld and gfld%ipdtmpl(27) are different

Cmp   believe gfld%ipdtmpl(9) matters as well
c

c      If a field not in a GRIB2 file, getGB2 output gfld will be crashed. 
c      so use previously saved gfld_save 

       nff=ff/1      
       output='prcip'//filehead(5:14)//fhr(nff)
       outdone='prcipdone'//filehead(9:14)//fhr(nff)

!        ounit=50+nfhr
        ounit=50+nff
        call baopen(ounit,output,ierr)



	if (.not. skip_1h) then

!!        Add a 1 h total for everyone
	   gfld=gfld_save_1h
           gfld%fld(:)=dp1(:)
           gfld%ipdtmpl(27)=1
	write(0,*) 'to putgb2 for ounit: ', ounit
	write(0,*) 'minval/maxval APCP 1h (gfld%fld(:)): ', 
     &    minval(gfld%fld(:)),maxval(gfld%fld(:))
             call putgb2_wrap(ounit,gfld,ierr)

	   gfld=gfld_save_1h_snow
           gfld%fld(:)=sn1(:)
           gfld%ipdtmpl(27)=1
	write(0,*) 'SNOW 1 h maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
             call putgb2_wrap(ounit,gfld,ierr)

	   gfld=gfld_save_1h_frzr
           gfld%fld(:)=fz1(:)
           gfld%ipdtmpl(27)=1
	write(0,*) 'FRZR 1 h maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
             call putgb2_wrap(ounit,gfld,ierr)

	   gfld=gfld_save_1h_asn
           gfld%fld(:)=asn1(:)
           gfld%ipdtmpl(27)=1
	write(0,*) 'ASN 1 h maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
             call putgb2_wrap(ounit,gfld,ierr)
    
!        write(0,*) 'Pack APCP done for fhr',nfhr
        write(0,*) 'Pack APCP done for fhr',nff, ' file: ', output

        endif

! write "done" file

        call baclose(ounit,ierr) 

C	if (filehead(6:8) .eq. 'm07') then
C	open(unit=11,file=outdone)
C	write(11,*) 'done'
C	close(unit=11)
C	endif

        return
      end subroutine just_hrly

c=========================================
      subroutine readGB2(igrb2,jpdtn,jpd1,jpd2,jpd27,gfld,iret)

        use grib_mod

        type(gribfield) :: gfld 
 
        integer jids(200), jpdt(200), jgdt(200)
        integer jpd1,jpd2,jpdtn
        logical :: unpck=.true. 
   

        jids=-9999  !array define center, master/local table, year,month,day, hour, etc, -9999 wildcard to accept any
        jpdt=-9999  !array define Product, to be determined
        jgdt=-9999  !array define Grid , -9999 wildcard to accept any

        jdisc=-1    !discipline#  -1 wildcard 
        jgdtn=-1    !grid template number,    -1 wildcard 
        jskp=0      !Number of fields to be skip, 0 search from beginning
        ifile=0

        jpdt(1)=jpd1   !Category #     
        jpdt(2)=jpd2   !Product # under this category     
        jpdt(27)=jpd27
        write(*,*) jpdtn,jpd1,jpd2,jpd27

         call getgb2(igrb2,ifile,jskp,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     +        unpck, jskp1, gfld,iret)

         
        return
        end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         
	subroutine just_hrrr_3hrly (filehead, ff, jf)

C  raw data
       use grib_mod
       real,allocatable,dimension(:,:) ::  dphold !jf,4        
       real,allocatable,dimension(:)   ::  dp1
       real,allocatable,dimension(:,:) ::  snhold !jf,4        
       real,allocatable,dimension(:)   ::  sn1
       integer iyr,imon,idy,ihr
       character*50 gdss(400)
       integer IENS, GRIBID, kgdss(200), lengds,im,jm,km,jf
       character*40 filehead,filename(8), output, outdone
       
       integer ff, nfm1,nfm2
       logical do_old
       character*3 fhr(48)
       integer iunit,ounit, pdt9_orig
       type(gribfield) :: gfld,gfld_save_curr
       type(gribfield) :: gfld_snow,gfld_save_curr_snow 

       data (fhr(i),i=1,48)
     + /'f01','f02','f03','f04','f05','f06','f07','f08','f09',
     + 'f10','f11','f12','f13','f14','f15','f16','f17','f18',
     + 'f19','f20','f21','f22','f23','f24','f25','f26','f27',
     + 'f28','f29','f30','f31','f32','f33','f34','f35','f36',
     + 'f37','f38','f39','f40','f41','f42','f43','f44','f45',
     + 'f46','f47','f48'/
 
	write(0,*) 'hrrr - know that filehead, ff,  jf', 
     *    trim(filehead), ff,  jf

       allocate(dphold(jf,3))
       allocate(dp1(jf))
       allocate(snhold(jf,3))
       allocate(sn1(jf))

!! these numbers need to change for hourly

        nff=ff/1
	nfm1=nff-1
	nfm2=nff-2

        dphold=0.
        snhold=0.

C       do 1001 nf=1,nfile
        
        filename(1)=filehead(1:14)//fhr(nff)
        filename(2)=filehead(1:14)//fhr(nfm1)
        filename(3)=filehead(1:14)//fhr(nfm2)

CCCCCCCCCCCCCCCCCCCCCCCCCC

        iunit=21
        call baopenr(iunit,filename(1),ierr)
        write(0,*) 'open ', iunit, filename(1), 'ierr=',ierr
	if (ierr /= 0) STOP

        jpdtn=8    !APCP's Product Template# is  4.8 

        jpd1=1
        jpd2=8
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP

        if (ie.eq.0) then
         dphold(:,1)=gfld%fld(:)
	 pdt9_orig=gfld%ipdtmpl(9)
         write(0,*) 'gfld(9) before saved: ', gfld%ipdtmpl(9)
         gfld%ipdtmpl(9)=-2 + pdt9_orig
         gfld_save_curr=gfld

         do i=1,gfld_save_curr%ipdtlen
           write(0,*) 'MOD=0 ', i, gfld_save_curr%ipdtmpl(i)
         end do

        end if

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         snhold(:,1)=gfld%fld(:)
         gfld%ipdtmpl(9)=-2 + pdt9_orig
         gfld_save_curr_snow=gfld
        endif


CCCCCCCCCCCCCCCCCCCCCCCCC

        iunit=22
        call baopenr(iunit,filename(2),ierr)
        write(0,*) 'open ', iunit, filename(2), 'ierr=',ierr
	if (ierr /= 0) STOP

        jpd1=1
        jpd2=8
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         dphold(:,2)=gfld%fld(:)
        end if

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         snhold(:,2)=gfld%fld(:)
        endif

CCCCCCCCCCCCCCCCCCCCCCCCC

        iunit=23
        call baopenr(iunit,filename(3),ierr)
        write(0,*) 'open ', iunit, filename(3), 'ierr=',ierr
	if (ierr /= 0) STOP

        jpd1=1
        jpd2=8
	jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         dphold(:,3)=gfld%fld(:)
        end if

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         snhold(:,3)=gfld%fld(:)
        endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



cccccc  Then call putgb2 to store the calculated data into a grib2 file
c
c      data structure gfld is re-used for pack data since all are same
c      only gfld%fld and gfld%ipdtmpl(27) are different

Cmp   believe gfld%ipdtmpl(9) matters as well
c

c      If a field not in a GRIB2 file, getGB2 output gfld will be crashed. 
c      so use previously saved gfld_save 

       nff=ff/1      
       output='prcip3h'//filehead(5:14)//fhr(nff)
       outdone='prcipdone'//filehead(9:14)//fhr(nff)

!        ounit=50+nfhr
        ounit=50+nff
        call baopen(ounit,output,ierr)



!!        Add a 1 h total for everyone
	   gfld=gfld_save_curr
           gfld%fld(:)=dphold(:,1)+dphold(:,2)+dphold(:,3)
           gfld%ipdtmpl(27)=3
	write(0,*) 'to putgb2 for ounit: ', ounit
	write(0,*) 'minval/maxval APCP 3h (gfld%fld(:)): ', 
     &    minval(gfld%fld(:)),maxval(gfld%fld(:))
         do i=1,gfld_save_curr%ipdtlen
         write(0,*) 'at output of apcp3h ',i,gfld%ipdtmpl(i)
         end do
           call putgb2_wrap(ounit,gfld,ierr)

	   gfld=gfld_save_curr_snow
           gfld%fld(:)=snhold(:,1)+snhold(:,2)+snhold(:,3)
           gfld%ipdtmpl(27)=3
	write(0,*) 'SNOW maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
             call putgb2_wrap(ounit,gfld,ierr)
    
!        write(0,*) 'Pack APCP done for fhr',nfhr
        write(0,*) 'Pack APCP done for fhr',nff

! write "done" file

        call baclose(ounit,ierr) 

      stop
      end 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         
	subroutine just_hrrr_3hrly_pre(filehead, ff, jf)

C  raw data
       use grib_mod
       real,allocatable,dimension(:,:) ::  dphold !jf,4        
       real,allocatable,dimension(:)   ::  dp1
       real,allocatable,dimension(:,:) ::  snhold !jf,4        
       real,allocatable,dimension(:)   ::  sn1
       integer iyr,imon,idy,ihr
       character*50 gdss(400)
       integer IENS, GRIBID, kgdss(200), lengds,im,jm,km,jf
       character*40 filehead,filename(8), output, outdone
       
       integer ff, nfm1,nfm2
       logical do_old
       character*3 fhr(48)
       character(len=6) :: term
       integer iunit,ounit, pdt9_orig
       type(gribfield) :: gfld,gfld_save_curr
       type(gribfield) :: gfld_snow,gfld_save_curr_snow 

       data (fhr(i),i=1,48)
     + /'f01','f02','f03','f04','f05','f06','f07','f08','f09',
     + 'f10','f11','f12','f13','f14','f15','f16','f17','f18',
     + 'f19','f20','f21','f22','f23','f24','f25','f26','f27',
     + 'f28','f29','f30','f31','f32','f33','f34','f35','f36',
     + 'f37','f38','f39','f40','f41','f42','f43','f44','f45',
     + 'f46','f47','f48'/
 
	write(0,*) 'hrrr 3hrly_pre - know that filehead, ff,  jf', 
     *    trim(filehead), ff,  jf

       term=".grib2"

       allocate(dphold(jf,3))
       allocate(dp1(jf))
       allocate(snhold(jf,3))
       allocate(sn1(jf))

!! these numbers need to change for hourly

        nff=ff/1
	nfm1=nff-1
	nfm2=nff-2

        dphold=0.
        snhold=0.

C       do 1001 nf=1,nfile
        
        filename(1)=filehead(1:10)//fhr(nff)//term
        filename(2)=filehead(1:10)//fhr(nfm1)//term
        filename(3)=filehead(1:10)//fhr(nfm2)//term

	write(0,*) 'filename(1): ', filename(1)

CCCCCCCCCCCCCCCCCCCCCCCCCC

        iunit=21
        call baopenr(iunit,filename(1),ierr)
        write(0,*) 'open ', iunit, filename(1), 'ierr=',ierr
	if (ierr /= 0) STOP

        jpdtn=8    !APCP's Product Template# is  4.8 

        jpd1=1
        jpd2=8
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP

        if (ie.eq.0) then
         dphold(:,1)=gfld%fld(:)
	 pdt9_orig=gfld%ipdtmpl(9)
         write(0,*) 'gfld(9) before saved: ', gfld%ipdtmpl(9)
         gfld%ipdtmpl(9)=-2 + pdt9_orig
         gfld_save_curr=gfld

         do i=1,gfld_save_curr%ipdtlen
           write(0,*) 'MOD=0 ', i, gfld_save_curr%ipdtmpl(i)
         end do

        end if

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         snhold(:,1)=gfld%fld(:)
         gfld%ipdtmpl(9)=-2 + pdt9_orig
         gfld_save_curr_snow=gfld
        endif


CCCCCCCCCCCCCCCCCCCCCCCCC

        iunit=22
        call baopenr(iunit,filename(2),ierr)
        write(0,*) 'open ', iunit, filename(2), 'ierr=',ierr
	if (ierr /= 0) STOP

        jpd1=1
        jpd2=8
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         dphold(:,2)=gfld%fld(:)
        end if

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         snhold(:,2)=gfld%fld(:)
        endif

CCCCCCCCCCCCCCCCCCCCCCCCC

        iunit=23
        call baopenr(iunit,filename(3),ierr)
        write(0,*) 'open ', iunit, filename(3), 'ierr=',ierr
	if (ierr /= 0) STOP

        jpd1=1
        jpd2=8
	jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         dphold(:,3)=gfld%fld(:)
        end if

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         snhold(:,3)=gfld%fld(:)
        endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



cccccc  Then call putgb2 to store the calculated data into a grib2 file
c
c      data structure gfld is re-used for pack data since all are same
c      only gfld%fld and gfld%ipdtmpl(27) are different

Cmp   believe gfld%ipdtmpl(9) matters as well
c

c      If a field not in a GRIB2 file, getGB2 output gfld will be crashed. 
c      so use previously saved gfld_save 

       nff=ff/1      
       output='prcip3h.t'//filehead(7:10)//fhr(nff)//term
       outdone='prcipdone'//filehead(7:8)//'_'//fhr(nff)

!        ounit=50+nfhr
        ounit=50+nff
        call baopen(ounit,output,ierr)



!!        Add a 1 h total for everyone
	   gfld=gfld_save_curr
           gfld%fld(:)=dphold(:,1)+dphold(:,2)+dphold(:,3)
           gfld%ipdtmpl(27)=3
	write(0,*) 'to putgb2 for ounit: ', ounit
	write(0,*) 'maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
         do i=1,gfld_save_curr%ipdtlen
         write(0,*) 'at output of apcp3h ',i,gfld%ipdtmpl(i)
         end do
           call putgb2_wrap(ounit,gfld,ierr)

	   gfld=gfld_save_curr_snow
           gfld%fld(:)=snhold(:,1)+snhold(:,2)+snhold(:,3)
           gfld%ipdtmpl(27)=3
	write(0,*) 'SNOW maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
             call putgb2_wrap(ounit,gfld,ierr)
    
!        write(0,*) 'Pack APCP done for fhr',nfhr
        write(0,*) 'Pack APCP done for fhr',nff

! write "done" file

        call baclose(ounit,ierr) 

      end subroutine just_hrrr_3hrly_pre

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         
	subroutine just_fv3_3hrly_pre(filehead, ff, jf)

C  raw data
       use grib_mod
       real,allocatable,dimension(:,:) ::  dphold !jf,4        
       real,allocatable,dimension(:)   ::  dp1
       real,allocatable,dimension(:,:) ::  snhold !jf,4        
       real,allocatable,dimension(:)   ::  sn1
       real,allocatable,dimension(:,:) ::  fzhold !jf,4        
       real,allocatable,dimension(:)   ::  fz1
       integer iyr,imon,idy,ihr
       character*50 gdss(400)
       integer IENS, GRIBID, kgdss(200), lengds,im,jm,km,jf
       character*40 filehead,filename(8), output, outdone
       
       integer ff, nfm1,nfm2
       logical do_old
       character*3 fhr(60)
       character(len=6) :: term
       integer iunit,ounit, pdt9_orig
       type(gribfield) :: gfld,gfld_save_curr
       type(gribfield) :: gfld_snow,gfld_save_curr_snow 
       type(gribfield) :: gfld_frzr,gfld_save_curr_frzr

       data (fhr(i),i=1,60)
     + /'f01','f02','f03','f04','f05','f06','f07','f08','f09',
     + 'f10','f11','f12','f13','f14','f15','f16','f17','f18',
     + 'f19','f20','f21','f22','f23','f24','f25','f26','f27',
     + 'f28','f29','f30','f31','f32','f33','f34','f35','f36',
     + 'f37','f38','f39','f40','f41','f42','f43','f44','f45',
     + 'f46','f47','f48','f49','f50','f51','f52','f53','f54',
     + 'f55','f56','f57','f58','f59','f60'/
 
	write(0,*) 'fv3 3hrly_pre - know that filehead, ff, jf', 
     *    trim(filehead), ff,  jf

       term=".grib2"

       allocate(dphold(jf,3))
       allocate(dp1(jf))
       allocate(snhold(jf,3))
       allocate(sn1(jf))
       allocate(fzhold(jf,3))
       allocate(fz1(jf))

!! these numbers need to change for hourly

        nff=ff/1
	nfm1=nff-1
	nfm2=nff-2

        dphold=0.
        snhold=0.

C       do 1001 nf=1,nfile
        
        filename(1)=filehead(1:10)//fhr(nff)//term
        filename(2)=filehead(1:10)//fhr(nfm1)//term
        filename(3)=filehead(1:10)//fhr(nfm2)//term

	write(0,*) 'filename(1): ', filename(1)

CCCCCCCCCCCCCCCCCCCCCCCCCC

        iunit=21
        call baopenr(iunit,filename(1),ierr)
        write(0,*) 'open ', iunit, filename(1), 'ierr=',ierr
	if (ierr /= 0) STOP

        jpdtn=8    !APCP's Product Template# is  4.8 

        jpd1=1
        jpd2=8
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP

        if (ie.eq.0) then
         dphold(:,1)=gfld%fld(:)
	 pdt9_orig=gfld%ipdtmpl(9)
         write(0,*) 'gfld(9) before saved: ', gfld%ipdtmpl(9)
         gfld%ipdtmpl(9)=-2 + pdt9_orig
         gfld_save_curr=gfld

         do i=1,gfld_save_curr%ipdtlen
           write(0,*) 'MOD=0 ', i, gfld_save_curr%ipdtmpl(i)
         end do

        end if

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         snhold(:,1)=gfld%fld(:)
         gfld%ipdtmpl(9)=-2 + pdt9_orig
         gfld_save_curr_snow=gfld
        endif

        jpd1=1
        jpd2=225
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! FRZR
        if (ie.eq.0) then
         fzhold(:,1)=gfld%fld(:)
         gfld%ipdtmpl(9)=-2 + pdt9_orig
         gfld_save_curr_frzr=gfld
        endif


CCCCCCCCCCCCCCCCCCCCCCCCC

        iunit=22
        call baopenr(iunit,filename(2),ierr)
        write(0,*) 'open ', iunit, filename(2), 'ierr=',ierr
	if (ierr /= 0) STOP

        jpd1=1
        jpd2=8
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         dphold(:,2)=gfld%fld(:)
        end if

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! snow
        if (ie.eq.0) then
         snhold(:,2)=gfld%fld(:)
        endif

        jpd1=1
        jpd2=225
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! frzr
        if (ie.eq.0) then
         fzhold(:,2)=gfld%fld(:)
        endif

CCCCCCCCCCCCCCCCCCCCCCCCC

        iunit=23
        call baopenr(iunit,filename(3),ierr)
        write(0,*) 'open ', iunit, filename(3), 'ierr=',ierr
	if (ierr /= 0) STOP

        jpd1=1
        jpd2=8
	jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         dphold(:,3)=gfld%fld(:)
        end if

        jpd1=1
        jpd2=13
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  !Large scale APCP
        if (ie.eq.0) then
         snhold(:,3)=gfld%fld(:)
        endif

        jpd1=1
        jpd2=225
        jpd27=1 !1 hr accumulation
        call readGB2(iunit,jpdtn,jpd1,jpd2,jpd27,gfld,ie)  ! frzr
        if (ie.eq.0) then
         fzhold(:,3)=gfld%fld(:)
        endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



cccccc  Then call putgb2 to store the calculated data into a grib2 file
c
c      data structure gfld is re-used for pack data since all are same
c      only gfld%fld and gfld%ipdtmpl(27) are different

Cmp   believe gfld%ipdtmpl(9) matters as well
c

c      If a field not in a GRIB2 file, getGB2 output gfld will be crashed. 
c      so use previously saved gfld_save 

       nff=ff/1      
       output='prcip3h.t'//filehead(7:10)//fhr(nff)//term
       outdone='prcipdone'//filehead(7:8)//'_'//fhr(nff)

!        ounit=50+nfhr
        ounit=50+nff
        call baopen(ounit,output,ierr)



!!        Add a 1 h total for everyone
	   gfld=gfld_save_curr
           gfld%fld(:)=dphold(:,1)+dphold(:,2)+dphold(:,3)
           gfld%ipdtmpl(27)=3
	write(0,*) 'to putgb2 for ounit: ', ounit
	write(0,*) 'maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
         do i=1,gfld_save_curr%ipdtlen
         write(0,*) 'at output of apcp3h ',i,gfld%ipdtmpl(i)
         end do
           call putgb2_wrap(ounit,gfld,ierr)

	   gfld=gfld_save_curr_snow
           gfld%fld(:)=snhold(:,1)+snhold(:,2)+snhold(:,3)
           gfld%ipdtmpl(27)=3
	write(0,*) 'SNOW maxval(gfld%fld(:)): ', maxval(gfld%fld(:))
             call putgb2_wrap(ounit,gfld,ierr)
    
!        write(0,*) 'Pack APCP done for fhr',nfhr
        write(0,*) 'Pack APCP done for fhr',nff

! write "done" file

        call baclose(ounit,ierr) 

      end subroutine just_fv3_3hrly_pre


c=============================================

      subroutine putgb2_wrap(ounit,gfld,ierr)
       use grib_mod
       type(gribfield) :: gfld
       integer :: ounit
	real, allocatable :: grnd(:)
        real :: gmin, gmax
       integer ::  nbit

	write(0,*) 'into putgb2_wrap'
	write(0,*) 'gfld%ngrdpts: ', gfld%ngrdpts
	write(6,*) 'was: gfld%idrtmpl(2:3): ', gfld%idrtmpl(2:3)

	gfld%idrtmpl(2)=0
	gfld%idrtmpl(3)=5
	write(6,*) 'now: gfld%idrtmpl(2:3): ', gfld%idrtmpl(2:3)

	allocate(grnd(gfld%ngrdpts))

!  compute nbit
      call getbit(0,abs(gfld%idrtmpl(2)), 
     +    gfld%idrtmpl(3),gfld%ngrdpts,0,gfld%fld,
     +    grnd,gmin,gmax,nbit)

	write(6,*) 'gmin,gmax,nbit: ', gmin,gmax,nbit

         gfld%idrtmpl(4)=nbit

         call putgb2(ounit,gfld,ierr)

         deallocate(grnd)

      end subroutine putgb2_wrap

