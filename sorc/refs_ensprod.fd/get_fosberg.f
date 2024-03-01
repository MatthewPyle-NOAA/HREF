ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutine get_fosberg: Compute fosberg index mean/prob
c
c  Author: Binbin Zhou,  Nov. 18, 2015
c  Modification history:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   	subroutine  get_fosberg(nv,ifunit,jpdtn,jf,iens,Lm,Lp,Lt,
     +      derv_mn,derv_sp,derv_pr,wgt)

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

        INTEGER, intent(IN) :: nv, jf, iens 

        REAL,dimension(jf,Lm),intent(INOUT) :: derv_mn
        REAL,dimension(jf,Lm),intent(INOUT) :: derv_sp
        REAL,dimension(jf,Lp,Lt),intent(INOUT) :: derv_pr

        real,dimension(jf) :: u10,v10,rh2,t2m,q2
        real,dimension(jf,iens) :: findex

        REAL,dimension(iens)  :: apoint
        REAL                  :: amean,aspread,aprob



        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld 

        REAL wgt(30)
        INTEGER miss(iens)

        miss=0 
        mbrs=iens
        jpd27=-9999
        !jpdtn=0

        write(*,*) 'In get_fosberg', ifunit
        write(*,*) 'jf,iens,Lm,Lp,Lt=',jf,iens,Lm,Lp,Lt

        DO 101 irun=1,iens

           call readGB2(ifunit(irun),0,0,0,103,2,jpd27,gfld,eps,iret)
            if (iret.eq.0) then
             t2m(:)=gfld%fld(:)
            else
               write(*,*) 'read 0,0,103,2 error', iret
            end if


           call readGB2(ifunit(irun),0,1,1,103,2,jpd27,gfld,eps,iret)
            if (iret.eq.0) then
             rh2(:)=gfld%fld(:)
            else !RH 2m not exist, then use SPFH to compute it
              call readGB2(ifunit(irun),0,1,0,103,2,jpd27,gfld,eps,iret)
              if (iret.eq.0) then 
                q2(:)=gfld%fld(:)
                do i=1,jf
                 q=q2(i) 
                 t=t2m(i) - 273.15
                 call get_rh(q,t,rh)
                 rh2(i)=rh
                end do
               else
                 write(*,*) 'read 0,1,103,2 error', iret
              end if
            end if

           call readGB2(ifunit(irun),0,2,2,103,10,jpd27,gfld,eps,iret)
             if(iret.eq.0) then
               u10(:)=gfld%fld(:)
             else
               write(*,*) 'read 2,2,103,10 error', iret
             end if
           call readGB2(ifunit(irun),0,2,3,103,10,jpd27,gfld,eps,iret)
             if(iret.eq.0) then
              v10(:)=gfld%fld(:)
             else
               write(*,*) 'read 2,3,103,10 error', iret
             end if

           do i = 1, jf
              w10=sqrt(u10(i)*u10(i)+v10(i)*v10(i))
              t2=t2m(i)-273.15
              rh=rh2(i)
              call fosberg_index(t2,w10,rh,fidx)
              findex(i,irun)=fidx 
           end do

101     CONTINUE

c         write(*,*) 'Get data done '

           derv_mn = 0.
           derv_sp = 0.
           derv_pr = 0.

          do 200 lv=1,dMlvl(nv)                                  !dMlvl(nv) always =1 for precp type
            do 300 igrid = 1,jf

             apoint(:)=findex(igrid,:) 
             call getmean(apoint,iens,amean,aspread,miss,wgt)
             derv_mn(igrid,lv) = amean
             derv_sp(igrid,lv) = aspread

300        continue 
200      continue 

        do 400 lv=1,dPlvl(nv)     !dPlvl(nv) always is 1 for precp type
          do 500 lh =1,dTlvl(nv)
            do 600 igrid = 1, jf

              apoint(:) = findex(igrid,:)

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
     +                 miss,wgt)
                   derv_pr(igrid,lv,lh)=aprob
                 end if
              end if

600      continue
500     continue
400    continue

           return
           end


         subroutine get_rh(q,t,rh) 
          real q,t,rh
          
          if( t.ge.0.0) then
            tt=17.62*t/(243.12 + t)     !over liquid water
          else
            tt=22.46*t/(272.62 + t)     !over ice, assume no supercooled water
          end if

          p0=1013.0                     !Assume all at 1013mb
          es = 6.112 *2.71828 ** tt     !saturated pressure(mb)
          qs=0.622*es/p0                !saturated specific humidity (in kg/kg)
          rh=q/qs*100.0

          return 
          end 

          
