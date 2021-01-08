cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine LLWS: compute low level wind shear 
c     
c     Author: Binbin Zhou, Aug, 7, 2005
c     05-17-2013, B. Zhou: Upgraded to grib2 I/O
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine llws (nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lt,eps,  
     +          derv_mn,derv_sp,derv_pr,wgt)


        use grib_mod
        include 'parm.inc'
        parameter (lvl=10)     !max 10 levels is set to reduce memory
        !parameter (lvl=14)     !max 14 levels is set

c    for derived variables
        Character*4 dvname(maxvar)
        Integer dk5(maxvar), dk6(maxvar), dk4(maxvar)
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
        REAL,dimension(jf,Lm),intent(INOUT) :: derv_mn,derv_sp
        REAL,dimension(jf,Lm,Lt),intent(INOUT) ::  derv_pr

        real apoint(iens),LLWSapoint(iens),llws1(jf,iens)
        real ws,u10,v10,hsfc,up(lvl),vp(lvl),hp(lvl)            
                                                                
        INTEGER miss(iens) 
        real wgt(30)

        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld

         REAL, dimension(jf,iens)::U10m,V10m,HS
         REAL, dimension(jf,iens,lvl)::Upr,Vpr,Hpr               
         Integer p(50),lvls
         character*5 eps

         data (p(i),i=1,25)
     +  /1000,975,950,925,900,875,850,825,800,775,750,725,700,675,
     +   650,625,600,575,550,525,500,475,450,425,400/

CCCCccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   Get all paramters from current and previous files 
c   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        write(*,*) 'In LLWS .....'
        write(*,*) 'ifunit=', ifunit, eps

        lvls=10

        write(*,*)  nv,jf,iens,Lm,Lp,Lt,lvls

        !jpdtn=0
        jp27=-9999

         miss=0

         do irun=1,iens

          call readGB2(ifunit(irun),jpdtn,2,2,103,10,jp27,gfld,eps,ie) !U10m 
           U10m(:,irun)=gfld%fld
          call readGB2(ifunit(irun),jpdtn,2,3,103,10,jp27,gfld,eps,ie) !V10m 
           V10m(:,irun)=gfld%fld
          call readGB2(ifunit(irun),jpdtn,3,5,1,0,jp27,gfld,eps,ie) !Hsfc
           HS(:,irun)=gfld%fld

          !write(*,*) 'read surface data done'

           do k=1,lvls
              !write(*,*) 'Level=',p(k)
            call readGB2(ifunit(irun),jpdtn,2,2,100,p(k),jp27,gfld,
     +                   eps,ie)
                Upr(:,irun,k)=gfld%fld
            call readGB2(ifunit(irun),jpdtn,2,3,100,p(k),jp27,gfld,
     +                   eps,ie)
                Vpr(:,irun,k)=gfld%fld
            call readGB2(ifunit(irun),jpdtn,3,5,100,p(k),jp27,gfld,
     +                   eps,ie)      
                Hpr(:,irun,k)=gfld%fld
            end do
 
         end do


             do lv=1,dMlvl(nv)                       !for all  levels

              do igrid = 1,jf

               do i=1,iens
                 u10=U10m(igrid,i) 
                 v10=V10m(igrid,i) 
                 hsfc=HS(igrid,i)
                 do k=1,lvls
                  up(k)=Upr(igrid,i,k)
                  vp(k)=Vpr(igrid,i,k)
                  hp(k)=Hpr(igrid,i,k)
                 end do
  
                if (u10.eq.0.0.and.v10.eq.0.0.and.hsfc.eq.0.0) then
                  ws=0.0
                  goto 500
                end if

                 call get_llws(up,vp,hp,u10,v10,hsfc,jf,lvls,ws)

500               LLWSapoint(i)=ws
                  llws1(igrid,i)=ws

                end do 


                call getmean(LLWSapoint,iens,amean,aspread,
     +                miss,wgt)
                  derv_mn(igrid,lv)=amean
                  derv_sp(igrid,lv)=aspread

c           if(amean.gt.20.0) write(*,*) 'LLWS:',igrid,amean

                  !if(igrid.eq.10000) then
                  ! write(*,*) 'mean LLWS=',amean,aspread
                  !end if

                end do

              end do

              do lv=1,dPlvl(nv)

                do lh = 1, dTlvl(nv)
                              
                 do igrid = 1,jf

                    LLWSapoint(:)=llws1(igrid,:)

                    if(trim(dop(nv)).ne.'-') then
                     thr1 = dThrs(nv,lh)
                     thr2 = 0.  
             call getprob(LLWSapoint,iens,thr1,thr2,dop(nv),aprob,
     +               miss,wgt)
                     derv_pr(igrid,lv,lh)=aprob

c                  if(igrid.eq.10000) then
c        write(*,*) 'prob LLWS=', LLWSapoint,aprob,thr1,dop(nv)
c                  end if

                    else
                     if(lh.lt.dTlvl(nv)) then
                       thr1 = dThrs(nv,lh)
                       thr2 = dThrs(nv,lh+1)
             call getprob(LLWSapoint,iens,thr1,thr2,dop(nv),aprob,
     +                miss,wgt) 
                       derv_pr(igrid,lv,lh)=aprob
                     end if
                    end if

                 end do

                end do
              end do

c            write(*,*) 'LLWS done'

           return
           end


      subroutine get_llws (up,vp,hp,u10,v10,hsfc,jf,lv,llws)
c         up,vp : u,v on pressure levels
c            hp : height of pressure levels
c       u10,v10 : 10m u and v
c          hsfc : sfc height

        INTEGER jf,lv
        REAL up(lv),vp(lv),hp(lv)
        REAL u10,v10,hsfc,llws

        !write(*,*) 'in get_llws ...'
        !write(*,*) u10,v10,hsfc,jf,lv
        !do k=1,lv
        ! write(*,*)k, up(k),vp(k),hp(k)
        !end do

! mp - for safety, predefine kp to zero
        kp = 0
!mp

         z1 = 10. + hsfc
          if (z1 .lt. hp(1) ) then
            kp = 0
          else
            do k = 1,lv-1
             if(z1.ge.hp(k).and.z1.lt.hp(k+1)) then
               kp=k
             end if
            end do
          end if

          hz1 = hp(kp+1) - z1

          dh=0.0
! mp - for safety...was seeing undefined u2 value
         u2=u10
         v2=v10
! mp 

          if((hz1+10.0).gt.609.6) then                    !609.6m=2000ft, hz1+10 means from sfc instead form 10m
            u2=u10 + (up(kp+1)-u10)*599.6/hz1
            v2=v10 + (vp(kp+1)-v10)*599.6/hz1
            z2=hsfc+609.6
          else
            do k=kp+1,lv-1
              dh=dh+hp(k+1)-hp(k)
              if((dh+hz1+10).gt.609.6) then  !ie, reach to 2000 feet
               z2=hsfc+609.6
               rt=(z2-hp(k))/(hp(k+1)-hp(k))
               u2=up(k)+(up(k+1)-up(k))*rt
               v2=vp(k)+(vp(k+1)-vp(k))*rt
               k2=k
               goto 609
              end if
            end do
          end if

609       llws=abs(sqrt((u2-u10)**2+(v2-v10)**2))
     &              /609.9
!          llws=llws*1.943*609.9                             !unit--> knot/2000ft
          llws=llws*609.9                             !unit--> m per s / 2000ft

     
       return
       end

