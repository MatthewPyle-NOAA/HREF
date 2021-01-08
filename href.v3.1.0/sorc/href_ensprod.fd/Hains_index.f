cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine compute Hains fire weather Index
c     
c     Author: Binbin Zhou, April, 2013
c     Modification: 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine fire_weather (nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lt,
     +              derv_mn,derv_sp,derv_pr,wgt)

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
        REAL,dimension(jf,Lp,Lt),intent(INOUT) ::  derv_pr

        
        real  count, amean,aspread,aprob, flt_cnd(iens),wgt(30)

        integer miss(iens)

        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld

        real,dimension(jf,iens) :: T950,T850,Td850,Hsfc,
     +   T700, T500, Td700
        real index_hains(iens)
        character*5 eps

        !jpdtn=0
        jp27=-9999

        write(*,*) 'In Hains Index  .....'
        write(*,*) 'nv,ifunit,jf,iens,Lp,Lt',
     +              nv,ifunit,jf,iens,Lp,Lt
        write(*,*) 'dTlvl=',dTlvl
    
        miss=0

        k_loop: do k=1,iens
          jpd12=dProbLevel(nv,1)
           call readGB2(ifunit(k),jpdtn,0,0,100,925,jp27,gfld,
     +          eps,iret) !T950, hirew has no 950mb, use 925 to replace
           if(iret.eq.0) then
             T950(:,k)=gfld%fld
           else
             miss(k)=1
             cycle k_loop
           end if

           call readGB2(ifunit(k),jpdtn,0,0,100,850,jp27,gfld,
     +         eps,iret) !T850
           if(iret.eq.0) then
             T850(:,k)=gfld%fld
           else
             miss(k)=1
             cycle k_loop
           end if

           call readGB2(ifunit(k),jpdtn,0,6,100,850,jp27,gfld,
     +         eps, iret) !Td850
           if(iret.eq.0) then
             Td850(:,k)=gfld%fld
           else
             miss(k)=1
             cycle k_loop
           end if

           call readGB2(ifunit(k),jpdtn,3,5,1,0,jp27,gfld,
     +         eps, iret)     !Sfc height
           if(iret.eq.0) then
             Hsfc(:,k)=gfld%fld
           else
             miss(k)=1
             cycle k_loop
           end if

           call readGB2(ifunit(k),jpdtn,0,0,100,700,jp27,gfld,
     +         eps, iret) !T700
           if(iret.eq.0) then
             T700(:,k)=gfld%fld
           else
             miss(k)=1
             cycle k_loop
           end if

           call readGB2(ifunit(k),jpdtn,0,0,100,500,jp27,gfld,
     +         eps, iret) !T500
           if(iret.eq.0) then
             T500(:,k)=gfld%fld
           else
             miss(k)=1
             cycle k_loop
           end if

           call readGB2(ifunit(k),jpdtn,0,6,100,700,jp27,gfld,
     +         eps, iret) !Td700
           if(iret.eq.0) then
             Td700(:,k)=gfld%fld
           else
             miss(k)=1
             cycle k_loop
           end if

        end do k_loop


        do 600 igrid=1,jf

           do k=1,iens
            if(miss(k).eq.0) then

                call hains_index(t950(igrid,k),
     +                 t850(igrid,k),td850(igrid,k),
     +                 t700(igrid,k),td700(igrid,k),
     +                 t500(igrid,k),hsfc(igrid,k),indx) 
                index_hains(k)=indx*1.0

             else
                index_hains(k)=0.0
             end if
           end do

           do lv=1,dMlvl(nv)
             call getmean(index_hains,iens,amean,aspread,miss,wgt)
             derv_mn(igrid,lv) = amean
             derv_sp(igrid,lv) = aspread                    
           end do

           do lv=1,dPlvl(nv)
            do lh = 1, dTlvl(nv)
              thr1 = dThrs(nv,lh)
              thr2 = dThrs(nv,lh+1)
              call getprob(index_hains,iens,thr1,thr2,dop(nv),aprob,
     +        miss,wgt)
              derv_pr(igrid,lv,lh)=aprob 
            end do
           end do

600      continue

        return
        end



	subroutine hains_index(t950,t850,td850,t700,td700,t500,
     +             hsfc,indx)

         real t950,t850,td850,t700,td700,t500,hsfc

         if(hsfc.lt.300.0) then

            dt=t950-t850
            dd=t850-td850
            if (dt.lt.3.99) then
              K1=1
            else if (dt.ge.3.99.and.dt.lt.7.99) then
              K1=2
            else
              K1=3
            end if
            if (dd.lt.5.99) then
              K2=1
            else if (dd.ge.5.99.and.dt.lt.9.99) then
              K2=2
            else
              K2=3
            end if

          else if (hsfc.ge.300.0.and.hsfc.lt.1000.0) then
            dt=t850-t700
            dd=t850-td850 
            if (dt.lt.5.99) then
              K1=1
            else if (dt.ge.5.99.and.dt.lt.10.99) then
              K1=2
            else
              K1=3
            end if
            if (dd.lt.5.99) then
              K2=1
            else if (dd.ge.5.99.and.dt.lt.12.99) then
              K2=2
            else
              K2=3
            end if

          else
            dt=t700-t500
            dd=t700-td700
            if (dt.lt.17.99) then
              K1=1
            else if (dt.ge.17.99.and.dt.lt.21.99) then
              K1=2
            else
              K1=3
            end if
            if (dd.lt.14.99) then
              K2=1
            else if (dd.ge.14.99.and.dt.lt.20.99) then
              K2=2
            else
              K2=3
            end if
          end if

          indx=K1+k2

         return
         end
