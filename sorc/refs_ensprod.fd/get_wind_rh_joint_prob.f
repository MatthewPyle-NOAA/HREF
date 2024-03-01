ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutine get_wind_rh_joint_prob: Compute wind, RH joint prob
c
c  Author: Binbin Zhou,  Nov. 18, 2015
c  Modification history:
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   	subroutine  get_wind_rh_joint_prob(nv,ifunit,jpdtn,
     +               jf,iens,Lp,Lt,derv_pr,wgt)

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

        REAL,dimension(jf,Lp,Lt),intent(INOUT) :: derv_pr

        real,dimension(jf) :: u10,v10,rh2,t2m,q2
        real,dimension(jf,iens) :: joint

        REAL,dimension(iens)  :: apoint
        REAL                  :: aprob,jindex



        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld 

        REAL wgt(30)
        INTEGER miss(iens)

        miss=0 
        mbrs=iens
        jpd27=-9999
        !jpdtn=0

        write(*,*) 'In get_wind_rh_joint', ifunit
        write(*,*) 'jf,iens,Lp,Lt=',jf,iens,Lp,Lt

        DO 101 irun=1,iens

           call readGB2(ifunit(irun),jpdtn,1,1,103,2,jpd27,
     +                    gfld,eps,iret)
 
           if (iret.eq.0) then

             rh2(:)=gfld%fld(:)

           else !RH 2m not exist, then use SPFH to compute it

              call readGB2(ifunit(irun),jpdtn,0,0,103,2,jpd27,
     +                     gfld,eps,iret)
              if (iret.eq.0) then
                 t2m(:)=gfld%fld(:)
              else
                 write(*,*) 'read 0,0,103,2 error', iret
                 stop 123
              end if

              call readGB2(ifunit(irun),jpdtn,1,0,103,2,jpd27,
     +                     gfld,eps,iret)
              if (iret.eq.0) then 
                q2(:)=gfld%fld(:)
              else
                 write(*,*) 'read 1,0,103,2 error', iret
                 stop 456
              end if

              do i=1,jf
                q=q2(i) 
                t=t2m(i) - 273.15
                call get_rh(q,t,rh)
                rh2(i)=rh
              end do

           end if

           call readGB2(ifunit(irun),jpdtn,2,2,103,10,jpd27,
     +                   gfld,eps,iret)
             if(iret.eq.0) then
               u10(:)=gfld%fld(:)
             else
               write(*,*) 'read 2,2,103,10 error', iret
             end if
           call readGB2(ifunit(irun),jpdtn,2,3,103,10,jpd27,
     +                  gfld,eps,iret)
             if(iret.eq.0) then
              v10(:)=gfld%fld(:)
             else
               write(*,*) 'read 2,3,103,10 error', iret
             end if

           tw10=dthrs(nv,Lt)         !Lt=1 
           trh2=dthrs(nv,Lt+1)
           do i = 1, jf
              w10=sqrt(u10(i)*u10(i)+v10(i)*v10(i))
              rh=rh2(i)
              call get_joint(w10,rh,dop(nv),tw10,trh2,jindex)
              joint(i,irun)=jindex
           end do

101     CONTINUE

        write(*,*) 'Get wind10 and RH2 data done '

        derv_pr = 0.
        do 600 igrid = 1, jf
           apoint(:) = joint(igrid,:)
           call getprob(apoint,iens,1.0,0.0,'=',aprob,miss,wgt)
           derv_pr(igrid,Lp,Lt)=aprob
600     continue

        return
        end

         

         subroutine get_rh2m(q,t,rh) 
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

          
