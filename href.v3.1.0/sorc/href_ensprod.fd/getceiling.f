cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine getceiling: compute ceiling height
c     
c     Author: Binbin Zhou, Jun 1, 2006
c     Modification: 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine getceiling (nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lt,  
     +             derv_mn,derv_sp,derv_pr,wgt)
        
        use grib_mod
        include 'parm.inc'

c    for derived variables
        Character*4 dvname(maxvar)
        Integer dk5(maxvar), dk6(maxvar),dk4(maxvar)
        Character*1 dMsignal(maxvar), dPsignal(maxvar)
        Integer dMlvl(maxvar), dMeanLevel(maxvar,maxmlvl)
        Integer dPlvl(maxvar), dProbLevel(maxvar,maxplvl)
        Character*1 dop(maxvar)
        Integer dTlvl(maxvar)
        Real    dThrs(maxvar,maxtlvl)
        Integer MPairLevel(maxvar,maxmlvl,2)
        Integer PPairLevel(maxvar,maxplvl,2)


        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop


        INTEGER, intent(IN) :: nv, jf, iens
        REAL,dimension(jf,Lm),intent(INOUT) :: derv_mn
        REAL,dimension(jf,Lm),intent(INOUT) :: derv_sp
        REAL,dimension(jf,Lp,Lt),intent(INOUT) :: derv_pr

        REAL, dimension(jf,iens) :: tcld,cldb,hsfc

        real CEILapoint(iens),TCLDapoint(iens),CLDBapoint(iens),
     +       HSFCapoint(iens)
        real wgt(30)

        integer miss(iens), JJ

        Character*5 eps

        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld


        !jpdtn=0
        jp27=-9999

c        write(*,*) 'In getceiling .....'
c        write(*,*) 'nv,ifunit,jf,iens,Lp,Lt',
c     +              nv,ifunit,jf,iens,Lp,Lt

        miss=0
        loop400: do irun=1,iens

           call readGB2(ifunit(irun),jpdtn,6,1,200,0,jp27,gfld,
     +                  eps,iret) !total cloud
            if(iret.eq.0) then
             tcld(:,irun)=gfld%fld
            else
             miss(irun)=1
             cycle loop400
            end if


! --------------------

           call readGB2(ifunit(irun),jpdtn,3,5,1,0,jp27,gfld,
     +                 eps,iret)   ! surface hgt
            if(iret.eq.0) then
             hsfc(:,irun)=gfld%fld
            else
             miss(irun)=1
             cycle loop400
            end if
! --------------------

           call readGB2(ifunit(irun),jpdtn,3,5,2,0,jp27,gfld,
     +                 eps,iret)   !cloud base
            if(iret.eq.0) then            

! account for bmap

        if (jf .ne. 37910 .and. jf .ne. 70720) then

            do JJ=1,jf
            if (.not. gfld%bmap(JJ)) then
             cldb(JJ,irun)=-5000.
            else
             cldb(JJ,irun)=gfld%fld(JJ)
            endif
            enddo
         else
             cldb(:,irun)=gfld%fld
         endif

            else

           call readGB2(ifunit(irun),jpdtn,3,5,215,0,jp27,gfld,
     +                 eps,iret)   !cloud ceiling

            if(iret.eq.0) then            

	write(0,*) 'believe processing HRRR cloud ceiling'

! account for bmap

        if (jf .ne. 37910 .and. jf .ne. 70720) then

            do JJ=1,jf
            if (.not. gfld%bmap(JJ)) then
             cldb(JJ,irun)=-5000.
            else
!was             cldb(JJ,irun)=gfld%fld(JJ)+hsfc(JJ,irun)
             cldb(JJ,irun)=gfld%fld(JJ)
            endif
            enddo
         else
             cldb(:,irun)=gfld%fld+hsfc(:,irun)
         endif

            else
             miss(irun)=1
             cycle loop400
            end if


	endif


           end do loop400

        do 600 igrid = 1,jf

          TCLDapoint=tcld(igrid,:)
          CLDBapoint=cldb(igrid,:)
          HSFCapoint=hsfc(igrid,:)

	if (igrid .eq. 842084) then
	write(0,*) 'miss: ', miss
	write(0,*) 'TCLDapoint: ', TCLDapoint
	write(0,*) 'CLDBapoint: ', CLDBapoint
	write(0,*) 'HSFCapoint: ', HSFCapoint
        endif

          do i = 1, iens
            if(miss(i).eq.0) then
              if(CLDBapoint(i).lt.0.0) CEILapoint(i)=20000.0    !Dec. 30, 'le'->'lt'


              if(TCLDapoint(i).ge.50.0 .and.
     +          CLDBapoint(i).ge.0.0  ) then
                CEILapoint(i) = CLDBapoint(i) - HSFCapoint(i)
                if(CEILapoint(i).lt.0.0) CEILapoint(i)=0.0
	if (igrid .eq. 842084) then
          write(0,*) 'defined CEILapoint(i) ', i, CEILapoint(i)
	endif
              else
                CEILapoint(i)=20000.0
	if (igrid .eq. 842084) then
          write(0,*) 'defaulted CEILapoint(i) ', i, CEILapoint(i)
	endif
              end if
            end if
          end do               

          call get_cond_mean (CEILapoint, iens, 
     +      20000.0, amean, aspread,miss,wgt)

               derv_mn(igrid,1)=amean
               derv_sp(igrid,1)=aspread

       
            do 30 lv=1,dPlvl(nv)
              do lh = 1, dTlvl(nv)

                if(trim(dop(nv)).ne.'-') then
                 thr1 = dThrs(nv,lh)
                 thr2 = 0.
                 call getprob(CEILapoint,iens,
     +                thr1,thr2,dop(nv),aprob,miss,wgt)
                 derv_pr(igrid,lv,lh)=aprob
                else
                  if(lh.lt.dTlvl(nv)) then
                    thr1 = dThrs(nv,lh)
                    thr2 = dThrs(nv,lh+1)
                    call getprob(CEILapoint,iens,
     +                   thr1,thr2,dop(nv),aprob,miss,wgt)
                    derv_pr(igrid,lv,lh)=aprob
                   end if
                end if
             end do

30         continue  

600      continue

        return
        end
