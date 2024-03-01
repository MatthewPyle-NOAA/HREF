cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine get_icing: compte icing
c     
c     Author: Binbin Zhou, Apr, 25, 2009
c      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine get_icing (nv,ifunit,jpdtn,jf,iens,Lp,Lt, 
     +              derv_pr,wgt)

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
        REAL,dimension(jf,Lp,Lt),intent(INOUT) :: derv_pr

        real Tapoint(iens),Rapoint(iens), Wapoint(iens),
     +              Icing(iens),wgt(30)


         REAL, dimension(jf,iens) :: T, R, W  !Temperature, TH and upward-wind
         INTEGER miss(iens),missing(20,iens)

        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld

        character*5 eps

        jpd10=dk6(nv)
        !jpdtn=0
        jp27=-9999

c        write(*,*) 'In get_icing .....'
c        write(*,*) 'nv,ifunit,jf,iens,Lp,Lt,jpd10',
c     +              nv,ifunit,jf,iens,Lp,Lt,jpd10

        miss=0
        missing=0
        do 600 lv=1,dPlvl(nv) 

          jpd12=dProbLevel(nv,lv)
          
c          write(*,*) 'jpd10 jpd12=',jpd10, jpd12

          loop601: do irun=1,iens

           call readGB2(ifunit(irun),jpdtn,0,0,jpd10,jpd12,jp27,
     +       gfld,eps,iret)   !T
            if(iret.eq.0) then 
             T(:,irun)=gfld%fld
            else
             missing(lv,irun)=1
             cycle loop601
            end if

           call readGB2(ifunit(irun),jpdtn,1,1,jpd10,jpd12,jp27,
     +       gfld,eps,iret)   !RH
            if(iret.eq.0) then
             R(:,irun)=gfld%fld
            else
             missing(lv,irun)=1
             cycle loop601
            end if

           call readGB2(ifunit(irun),jpdtn,2,8,jpd10,jpd12,jp27,
     +       gfld,eps,iret)   !W
            if(iret.eq.0) then
             W(:,irun)=gfld%fld
            else
             missing(lv,irun)=1
             cycle loop601
            end if

           end do loop601
           

           do lh = 1, dTlvl(nv)
             do igrid = 1,jf

               miss=missing(lv,:)

               do i=1,iens
                 Tapoint(i)=T(igrid,i)
                 Rapoint(i)=R(igrid,i)
                 Wapoint(i)=W(igrid,i)
                 IF(Wapoint(i).lt.0.0.AND.
     +             (Tapoint(i).LE.273.0.AND.Tapoint(i).GE.251.0)
     +            .AND. Rapoint(i).GE.70.0) THEN
                     Icing(i) = 1.0
                 ELSE
                     Icing(i) = 0.0
                 END IF
                end do

               thr1 = dThrs(nv,lh)
               thr2 = 0.
               call getprob(Icing,iens,thr1,thr2,dop(nv),aprob,
     +               miss,wgt)

               derv_pr(igrid,lv,lh)=aprob

              end do  ! end of igrid
            end do ! end of dTlvl
600      continue  ! end of dPlvl

           return
           end
