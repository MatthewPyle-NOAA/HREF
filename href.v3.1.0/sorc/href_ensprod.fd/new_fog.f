


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine new_fog: compute fog LWC according to the asymptotic analysis 
c                         of radiation fog theory but adding advection term
c     
c     Author: Binbin Zhou, Aug 18 , 2010 
c     Modification: April 23, 2013: Modified for grib2 I/O
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine new_fog(nv,ifunit,ipunit,jpdtn,jf,im,jm,dx,dy,dt,
     +     iens,Lm,Lp,Lt,derv_mn,derv_sp,derv_pr,wgt)

        use grib_mod
        include 'parm.inc'

                                               
        type (gribfield) :: gfld                                                                                                                 
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
        REAL,dimension(jf,Lm),intent(INOUT) :: derv_mn, derv_sp
        REAL,dimension(jf,Lp,Lt),intent(INOUT) ::  derv_pr
        INTEGER :: JJ
        
        real apoint(iens),FOGapoint(iens)
        real u10,v10,hsfc,up(10),vp(10),hp(10),dt2,dtp(10),
     +       tp(10),t2,rh2,rhp(10),lwc,lwc1(jf,iens)

        REAL wgt(30) 

        real qw(jf),qw2d(im,jm),qadv(jf,iens),
     +       qadv2d(im,jm),u2d(im,jm),v2d(im,jm)
        real dx,dy,ddx,ddy,q_adv

        real rhthrsh(jf)

        integer miss(iens)

         !following fields data are required to read in from all members
         REAL, dimension(jf,iens)::U10m,V10m,RH2m,CLDBS,CLDTP,HS
         REAL, dimension(jf,iens)::T2m, dT2m
         REAL, dimension(jf,iens,10)::Upr,Vpr,RHpr,Hpr
         REAL, dimension(jf,iens,10)::Tpr,dTpr

         Integer p(14),dt 
         integer,dimension(iens),intent(IN) :: ifunit,ipunit

         character*5 eps

         data (p(i),i=1,14) 
     +  /1000,975,950,925,900,875,850,825,800,775,750,725,700,675/

CCCCccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Get all paramters from current and previous files 
c   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        eps='any  '
        write(*,*) 'In new_fog .....' 
        write(*,*) 'ifunit=', ifunit
        write(*,*) 'ipunit=', ipunit
        write(*,*)  nv,jf,im,jm,dx,dy,iens,Lm,Lp,Lt
        !jpdtn=0
          if (jpdtn.eq.0) then
           jpdtnp=8    !Temperature variation in prcip file is in Template 8
          else if (jpdtn.eq.1) then
           jpdtnp=11
          else 
           jpdtnp=8
          end if

        jp27=-9999

        miss=0

      bigloop:  do irun=1,iens
          !jpdtn=0
          call readGB2(ifunit(irun),jpdtn,2,2,103,10,jp27,gfld,eps,ie) !U10m 
          if(ie.eq.0) then
           U10m(:,irun)=gfld%fld
          else
             miss(irun)=1
             cycle bigloop
          end if

          call readGB2(ifunit(irun),jpdtn,2,3,103,10,jp27,gfld,eps,ie) !V10m 
          if(ie.eq.0) then
           V10m(:,irun)=gfld%fld
          else
             miss(irun)=1
             cycle bigloop
          end if

          call readGB2(ifunit(irun),jpdtn,1,1,103,2,jp27,gfld,eps,ie) !RH2m
          if(ie.eq.0) then
           RH2m(:,irun)=gfld%fld
          else
             miss(irun)=1
             cycle bigloop
          end if

          call readGB2(ifunit(irun),jpdtn,0,0,103,2,jp27,gfld,eps,ie) !T2m(current)
          if(ie.eq.0) then
           T2m(:,irun)=gfld%fld-273.15
          else
             miss(irun)=1
             cycle bigloop
          end if

          call readGB2(ifunit(irun),jpdtn,3,5,2,0,jp27,gfld,eps,ie) !CLDBS
          if(ie.eq.0) then
! account for bmap

        if (jf .ne. 37910 .and. jf .ne. 70720) then

            do JJ=1,jf
            if (.not. gfld%bmap(JJ)) then
             CLDBS(JJ,irun)=-5000.
            else
             CLDBS(JJ,irun)=gfld%fld(JJ)
            endif
            enddo
         else
             CLDBS(:,irun)=gfld%fld
         endif

          else
             miss(irun)=1
             cycle bigloop
          end if

          call readGB2(ifunit(irun),jpdtn,3,5,3,0,jp27,gfld,eps,ie) !CLDT
          if(ie.eq.0) then
           CLDTP(:,irun)=gfld%fld
          else
             miss(irun)=1
             cycle bigloop
          end if

          call readGB2(ifunit(irun),jpdtn,3,5,1,0,jp27,gfld,eps,ie) !Hsfc
          if(ie.eq.0) then
           HS(:,irun)=gfld%fld
          else
             miss(irun)=1
             cycle bigloop
          end if

          call readGB2(ipunit(irun),jpdtnp,0,0,103,2,jp27,gfld,eps,ie)       !dT2m

          if(ie.eq.0) then
           dT2m(:,irun)=gfld%fld
          else
             miss(irun)=1
             cycle bigloop
          end if
 
           do k=1,10
            lvl=p(k)
              !jpdtn=0
              call readGB2(ifunit(irun),jpdtn,2,2,100,lvl,jp27,gfld,
     +           eps,ie) 
                Upr(:,irun,k)=gfld%fld
              call readGB2(ifunit(irun),jpdtn,2,3,100,lvl,jp27,gfld,
     +           eps,ie) 
                Vpr(:,irun,k)=gfld%fld
              call readGB2(ifunit(irun),jpdtn,1,1,100,lvl,jp27,gfld,
     +            eps,ie)  
                RHpr(:,irun,k)=gfld%fld
              call readGB2(ifunit(irun),jpdtn,3,5,100,lvl,jp27,gfld,
     +            eps,ie)      
                Hpr(:,irun,k)=gfld%fld
              call readGB2(ifunit(irun),jpdtn,0,0,100,lvl,jp27,gfld,
     +            eps,ie)      
                Tpr(:,irun,k)=gfld%fld-273.15
               
              call readGB2(ipunit(irun),jpdtnp,0,0,100,lvl,jp27,gfld,
     +            eps,ie)      
                dTpr(:,irun,k)=gfld%fld
            end do

           end do bigloop

           rhthrsh=98.0
 
           ddx=dx/1000.
           ddy=dy/1000.

           qadv=0.

        do 500 k=1,iens

            qadv2d=0.

           do igrid = 1,jf

            if(miss(k).eq.0) then

             RH=RH2m(igrid,k)/100.0
             t=T2m(igrid,k)

c             tt=7.45*t/(235.0 + t)
c             es = 6.10 *10 ** tt                        !es:saturated vapor pressure, qw: specific humiduty
C        Use WMO recommended formulation (2008)
C   Guide to Metteorological Instruments and Methods of
C   Observation (CIMO Guide):

              if( tfog.ge.0.0) then
                tt=17.62*t/(243.12 + t)     !over liquid water
              else
                tt=22.46*t/(272.62 + t)     !ice fog condition, assume no supercooled water
              end if

              es = 6.112 *2.71828 ** tt
              qw(igrid)=622.0*es*RH/1000.0               !qw=622*es*RH/P    (g/kg)
             end if ! end of if missing

           end do  !end of igrid

           do j=1,jm
             do i=1,im
              ij=(j-1)*im + i
              if(miss(k).eq.0) then
                qw2d(i,j)=qw(ij)
                u2d(i,j)=U10m(ij,k)
                v2d(i,j)=V10m(ij,k)
              else
                qw2d(i,j)=0.
                u2d(i,j)=0.
                v2d(i,j)=0.
              end if     
             end do
            end do

            do j=2,jm-1                                            !general upwind scheme
             do i=2,im-1
              ij=(j-1)*im + i
              if(miss(k).eq.0) then
              
               if(u2d(i,j).ge.0.0.and.v2d(i,j).ge.0.0) then
                  qadv2d(i,j) =
     +            -u2d(i,j)*(qw2d(i,j)-qw2d(i-1,j))/ddx
     +            -v2d(i,j)*(qw2d(i,j)-qw2d(i,j-1))/ddy
          else if(u2d(i,j).gt.0.0.and.v2d(i,j).lt.0.0) then
                 qadv2d(i,j) =
     +            -u2d(i,j)*(qw2d(i,j)-qw2d(i-1,j))/ddx
     +            -v2d(i,j)*(qw2d(i,j+1)-qw2d(i,j))/ddy
          else if(u2d(i,j).lt.0.0.and.v2d(i,j).lt.0.0) then
                  qadv2d(i,j) =
     +            -u2d(i,j)*(qw2d(i+1,j)-qw2d(i,j))/ddx
     +            -v2d(i,j)*(qw2d(i,j+1)-qw2d(i,j))/ddy
          else if(u2d(i,j).lt.0.0.and.v2d(i,j).gt.0.0) then
                  qadv2d(i,j) =
     +             -u2d(i,j)*(qw2d(i+1,j)-qw2d(i,j))/ddx
     +             -v2d(i,j)*(qw2d(i,j)-qw2d(i,j-1))/ddy
               end if

              end if
             end do
            end do

            do j=1,jm
             do i=1,im
              ij=(j-1)*im + i
              qadv(ij,k)=qadv2d(i,j)
             end do
            end do

500       continue   !end of  k~iens


        do 600 lv=1,dMlvl(nv)                       
          do igrid = 1,jf

            do i=1,iens
              if(miss(i).eq.0) then
               
                 u10= U10m(igrid,i)
                 v10= V10m(igrid,i)
                 hsfc= HS(igrid,i)
                 rh2= RH2m(igrid,i)
                 t2 = T2m(igrid,i)     !this time step
                 cldb = CLDBS(igrid,i)
                 cldt = CLDTP(igrid,i)
                 q_adv=qadv(igrid,i)
                 dt2=dT2m(igrid,i)

                 do k=1,10
                  up(k)=Upr(igrid,i,k)
                  vp(k)=Vpr(igrid,i,k)
                  hp(k)=Hpr(igrid,i,k)
                  rhp(k)=RHpr(igrid,i,k)  
                  tp(k)=TPr(igrid,i,k)     !this time step
                  dtp(k)=dTPr(igrid,i,k)
                 end do
  
                 if( t2.ge.-30.0) then
                   rh_thrsh=rhthrsh(igrid)              !warm fog saturation threshold     
                 else
                   ttw=17.62*t2/(243.12 + t2)     !over water
                   tti=22.46*t2/(272.62 + t2)     !over ice
                                                        !esi = 6.112 *2.71828 ** tti
                                                        !esw = 6.112 *2.71828 ** ttw
                   rh_thrsh = 2.71828**(tti-ttw)*100.   !rh_thrsh = esi / esw, ice fog saturation threshold
                 end if

                 lwc=0.0
                 !if(i.eq.1) then
                 !  write(*,*) igrid,i,u10,v10,t2,rh2,hsfc,
     +           !     cldt,cldb,q_adv,rh_thrsh,dt,dt2 
                 !  do k=1,14
                 !write(*,'(6f9.2)') up(k),vp(k),hp(k),rhp(k),
     +           !     tp(k),dtp(k)
                 !  end do
                 !end if

                 
                 call get_new_fog(igrid,up,vp,hp,u10,v10,hsfc,
     +             rhp,rh2,t2,tp,dt2,dtp,dt,10,cldt,cldb,q_adv,
     +             rh_thrsh,lwc)
                 FOGapoint(i)=lwc 
                 lwc1(igrid,i)=lwc
               end if
            end do  ! end if iens 

            call getmean_fog(FOGapoint,iens,amean,aspread,
     +          miss,wgt)

             derv_mn(igrid,lv)=amean
             derv_sp(igrid,lv)=aspread

          end do
600     continue   


         do 700 lv=1,dPlvl(nv)
           do lh = 1, dTlvl(nv)
             do igrid = 1,jf

               FOGapoint(:)=lwc1(igrid,:)

               if(trim(dop(nv)).ne.'-') then
                thr1 = dThrs(nv,lh)
                thr2 = 0.
                call getprob(FOGapoint,iens,thr1,thr2,dop(nv),
     +               aprob,miss,wgt)
                     derv_pr(igrid,lv,lh)=aprob
                else
                 if(lh.lt.dTlvl(nv)) then
                   thr1 = dThrs(nv,lh)
                   thr2 = dThrs(nv,lh+1)
                   call getprob(FOGapoint,iens,thr1,thr2,dop(nv),
     +                  aprob,miss,wgt)
                        derv_pr(igrid,lv,lh)=aprob
                  end if
                end if

                 end do
             end do
700        continue

           return
           end




