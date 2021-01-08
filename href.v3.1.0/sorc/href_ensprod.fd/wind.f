cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine wind: compute wind speed at 10m or high pressure levels
c     first search for index of U, V from the table -> ID_U, ID_V, if found:
c     then search level index from MeanLevel-> L, if found:
c     then compute the wind  W = SQRT(U*U+V*V)
c     
c     Author: Binbin Zhou, Aug, 5, 2005
c     Modification: 
c      03/01/2006: Binbin ZHou: modify wind speed spread computation method
c                   by considering the directions of wind vectors
c
c      04/10/2013: Binbin Z. Modified for grib2 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine wind (nv,ifunit,jpdtn,jf,im,jm,iens,Lm,Lp,Lt,
     +        derv_mn,derv_sp,derv_pr,wgt,mbrname)

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
 
        Character*5 eps                                                                                                                                                                
        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop


        INTEGER, intent(IN) :: nv, jf, iens, im, jm
        REAL,dimension(jf,Lm),intent(INOUT) :: derv_mn
        REAL,dimension(jf,Lm),intent(INOUT) :: derv_sp
        REAL,dimension(jf,Lp,Lt),intent(INOUT) :: derv_pr

        REAL, dimension(jf,iens) :: u,v      ! temporal vars      
        real, allocatable :: windspd(:,:)

        INTEGER miss(iens)
        INTEGER missing(20,iens)
        character*7 mbrname(50)

        real apoint(iens),Uapoint(iens),Vapoint(iens),wgt(30)
        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld

        jpd10=dk6(nv)
        !jpdtn=0
        jp27=-9999

        write(*,*) 'In wind .....'
        write(0,*) 'trim(dpsignal(nv)): ', trim(dpsignal(nv))
        write(*,*) 'nv,ifunit,jf,iens,Lm,Lp,Lt,jpd10',
     +              nv,ifunit,jf,iens,Lm,Lp,Lt,jpd10


	allocate(windspd(jf,iens))

        miss=0
        missing=0
        DO 500 lv=1,dMlvl(nv)

          jpd12=dMeanLevel(nv,lv)
          
          loop400: do irun=1,iens

           if(mbrname(irun).eq.'hrrrgsd'.and.jpd12.eq.80) then  !GSD HRRR use category 1 instead of 2 for 80m U and V
             jpd1=1
           else
             jpd1=2
           end if
           


	if (jpd12 .eq. 10) then


	write(0,*) 'jpdtn, jpd1: ', jpdtn, jpd1
           call readGB2(ifunit(irun),8,2,222,103,10,jp27,
     +       gfld,eps,iret)   !UMAX mean -component

           if (iret .eq. 0) then
             write(0,*) 'found UMAX'
             u(:,irun)=gfld%fld
           else

          
	write(0,*) 'no UMAX'
           call readGB2(ifunit(irun),jpdtn,jpd1,2,jpd10,jpd12,jp27,
     +       gfld,eps,iret)   !U mean -component
            if (iret.eq.0) then
             u(:,irun)=gfld%fld
            else
             missing(lv,irun)=1
             cycle loop400
            end if
          endif

        else ! non level=10

           call readGB2(ifunit(irun),jpdtn,jpd1,2,jpd10,jpd12,jp27,
     +       gfld,eps,iret)   !U mean -component
            if (iret.eq.0) then
             u(:,irun)=gfld%fld
            else
             missing(lv,irun)=1
             cycle loop400
            end if

       endif

	if (jpd12 .eq. 10) then

           call readGB2(ifunit(irun),8,2,223,103,10,jp27,
     +       gfld,eps,iret)   !VMAX mean -component

           if (iret .eq. 0) then
	     write(0,*) 'found VMAX'
             v(:,irun)=gfld%fld
           else
	     write(0,*) 'no VMAX'
           call readGB2(ifunit(irun),jpdtn,jpd1,3,jpd10,jpd12,jp27,
     +       gfld,eps,iret)   !V mean -component
            if (iret.eq.0) then
             v(:,irun)=gfld%fld
            else
             missing(lv,irun)=1
             cycle loop400
            end if

           endif
         

         else ! not level 10
           call readGB2(ifunit(irun),jpdtn,jpd1,3,jpd10,jpd12,jp27,
     +       gfld,eps,iret)   !V mean -component
            if (iret.eq.0) then
             v(:,irun)=gfld%fld
            else
             missing(lv,irun)=1
             cycle loop400
            end if


         endif

	write(0,*) 'define windspd for irun: ', irun
           do igrid = 1,jf
             windspd(igrid,irun)=sqrt(u(igrid,irun)*u(igrid,irun)+
     +                             v(igrid,irun)*v(igrid,irun))
           enddo

	write(0,*) 'windspd(igrid/2,:) ', windspd(igrid/2,:)

           end do loop400
     
           write(*,*) 'get wind speed mean for level ',jpd12

           do igrid = 1,jf
             Uapoint=u(igrid,:)
             Vapoint=v(igrid,:)
             miss=missing(lv,:)
             call getwindmean(Uapoint,Vapoint,iens,amean,
     +            aspread,miss,wgt)
             derv_mn(igrid,lv)=amean
             derv_sp(igrid,lv)=aspread
            end do

500      CONTINUE

        miss=0
        missing=0


        DO 600 lv=1,dPlvl(nv)

          jpd12=dProbLevel(nv,lv)

          loop401: do irun=1,iens

           if(mbrname(irun).eq.'hrrrgsd'.and.jpd12.eq.80) then  !GSD HRRR use category 1 instead of 2 for 80m U and V
             jpd1=1
           else
             jpd1=2
           end if


           call readGB2(ifunit(irun),jpdtn,jpd1,2,jpd10,jpd12,jp27,
     +       gfld,eps,iret)   !U Prob -component
            if (iret.eq.0) then
             u(:,irun)=gfld%fld
            else
             missing(lv,irun)=1
             cycle loop401
            end if
           call readGB2(ifunit(irun),jpdtn,jpd1,3,jpd10,jpd12,jp27,
     +       gfld,eps,iret)   !V Prob -component
            if (iret.eq.0) then
             v(:,irun)=gfld%fld
            else
             missing(lv,irun)=1
             cycle loop401
            end if

          end do loop401      

          write(*,*) 'get wind speed prob for level', jpd12
!          write(*,*) 'missing=',missing(lv,:) 


	if ( trim(dpsignal(nv)) .eq. 'K' .or. trim(dpsignal(nv)) .eq. 'L') then

!       neighborhood max

	do irun = 1, iens

!	write(0,*) 'irun: ', irun
!       write(0,*) 'dPsignal(nv): ', dPsignal(nv)
!	write(0,*) 'shape(windspd): ', shape(windspd)
!	write(0,*) 'pre windspd(jf/2,irun): ', windspd(jf/2,irun)


!	write(0,*) 'maxval(windspd(:,irun)): ', maxval(windspd(:,irun))


!! this modifies rawdata_pr
             call neighborhood_max(windspd(:,irun),
     +                          jf,im,jm,dPsignal(nv))

!	write(0,*) 'post windspd(jf/2,irun): ', windspd(jf/2,irun)
!	write(0,*) 'post maxval(windspd(:,irun)): ', maxval(windspd(:,irun))
       enddo


!	write(0,*) 'got past derv_pr for wind neighb max'

        endif

          do lh = 1, dTlvl(nv)
           do igrid = 1,jf
             miss=missing(lv,:)
	if ( trim(dpsignal(nv)) .eq. 'K' .or. trim(dpsignal(nv)) .eq. 'L') then
             apoint=windspd(igrid,:)
!	if (igrid .eq. jf/2) then
!           write(0,*) 'defined apoint(nbrmax) as: ', apoint
!        endif

        else

             apoint=sqrt(u(igrid,:)*u(igrid,:)+v(igrid,:)*v(igrid,:))
!	if (Igrid .eq. jf/2) then
!           write(0,*) 'defined apoint(pnt) as: ', apoint
!        endif

        endif
              if(trim(dop(nv)).ne.'-') then
                     thr1 = dThrs(nv,lh)
                     thr2 = 0.
                  call getprob(apoint,iens,thr1,thr2,dop(nv),aprob,
     +                miss,wgt)
                  derv_pr(igrid,lv,lh)=aprob
              else
                 if(lh.lt.dTlvl(nv)) then
                   thr1 = dThrs(nv,lh)
                   thr2 = dThrs(nv,lh+1)
                   call getprob(apoint,iens,thr1,thr2,dop(nv),aprob,
     +                miss,wgt)
                   derv_pr(igrid,lv,lh)=aprob
                 end if
              end if

           end do
          end do

600     CONTINUE

           return
           end
