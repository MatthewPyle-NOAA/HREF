ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutine preciptype: Compute precipitation type rain/freezing rain/snow
c  based on Jun Du old version
c
c  Author: Binbin Zhou, Aug. 4, 2005 
c  Modification history:
c   01/12/2006: Geoff Manikin: Add Geoff Manikin algorithm to determine dominant precip_type 
c   04/10/2013: Binbin Z. Modified for grib2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   	subroutine  preciptype(nv,ifunit,jpdtn,jf,iens,
     +      ptype_mn,ptype_pr,ptype_pr2)

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

        REAL,dimension(jf,4),intent(INOUT) :: ptype_mn
        REAL,dimension(jf,4),intent(INOUT) :: ptype_pr,ptype_pr2

        INTEGER miss(iens),pcount
        integer,dimension(iens),intent(IN) :: ifunit
        type(gribfield) :: gfld 
        integer  jpd1(4),jpd2(4),jpd10(4),jpd12(4)
        real prcptype(jf,iens,4)                 !store 4 precip types data

        data (jpd1(i),i=1,4)
     +  /1,1,1,1/
        data (jpd2(i),i=1,4)
     +  /192,193,194,195/                     !rain, frzr, icep, snow
        data (jpd10(i),i=1,4)
     +  /1,1,1,1/
        data (jpd12(i),i=1,4)
     +  /0,0,0,0/

        write(*,*) ' In preciptype ......'
        write(*,*) nv,ifunit,jf,iens

        miss=0
        mbrs=iens
        !First get necessary data from raw ensemble members 
        jpd27=-9999
        !jpdtn=0
!        loop101: DO irun=1,iens
        loop101: DO irun=iens,1,-1
          do i=1,4
           call readGB2(ifunit(irun),jpdtn,jpd1(i),jpd2(i),
     +       jpd10(i),jpd12(i),jp27,gfld,eps,iret)
           if(iret.eq.0) then
!	write(0,*) 'irun, i, max(gfld%fld): ', 
!     &          irun, i, max(gfld%fld)
             prcptype(:,irun,i)=gfld%fld
           else
             miss(irun)=1
             mbrs=mbrs - 1
             cycle loop101
             write(0,*) 'set miss for irun: ', irun
           end if
          end do
         end do loop101

         write(0,*) 'Get ptype data done '

! eliminate these prints?

	do irun=1,iens
                  if(miss(irun).eq.0) then
          write(0,*) 'irun, min/maxval(prcptype(:,irun,1))' , 
     &   irun,minval(prcptype(:,irun,1)),maxval(prcptype(:,irun,1))
                  endif
        enddo
	do irun=1,iens
                  if(miss(irun).eq.0) then
          write(0,*) 'irun, min/maxval(prcptype(:,irun,2))' , 
     &   irun,minval(prcptype(:,irun,2)),maxval(prcptype(:,irun,2))
                  endif
        enddo
	do irun=1,iens
                  if(miss(irun).eq.0) then
          write(0,*) 'irun, min/maxval(prcptype(:,irun,3))' , 
     &   irun,minval(prcptype(:,irun,3)),maxval(prcptype(:,irun,3))
                  endif
        enddo
	do irun=1,iens
                  if(miss(irun).eq.0) then
          write(0,*) 'irun, min/maxval(prcptype(:,irun,4))' , 
     &   irun,minval(prcptype(:,irun,4)),maxval(prcptype(:,irun,4))
                  endif
        enddo

           ptype_mn = 0.
           ptype_pr = 0.
           ptype_pr2 = 0.

              
          do 200 lv=1,dMlvl(nv)                                  !dMlvl(nv) always =1 for precp type
          loop300:  do igrid = 1,jf
         
                 crain = 0.
                 cfrzr = 0.
                 cslet = 0.
                 csnow = 0.

                 do irun = 1,iens 
                  if(miss(irun).eq.0) then
                   crain = crain + prcptype(igrid,irun,1)
                   cfrzr = cfrzr + prcptype(igrid,irun,2)
                   cslet = cslet + prcptype(igrid,irun,3) 
                   csnow = csnow + prcptype(igrid,irun,4)
                  end if
                 end do  


! SLR test
!          frztyp=cslet+csnow
!           if (frztyp .gt. 0) then
!              fracsn=csnow/frztyp
!              fracip=cslet/frztyp
!              SLR_derv(igrid)=fracsn*10.+2*fracip
!             else
!              SLR_derv(igrid)=10.
!            endif
! end SLR test

cc  following is part is the code copy from Geoff Manikin dominant precip type decision
cc  importance priority order:  freezing_rain(1) > snow(2) > sleet(3) > rain(4)            

          if(crain.ge.1.0.or.cfrzr.ge.1.0.or.csnow.ge.1.0.or.
     +         cslet.ge.1.0) then

            if(csnow.gt.cslet) then

              if(csnow.gt.cfrzr) then
                  if(csnow.ge.crain) then
                   ptype_mn(igrid,4)=1.      !snow
                   cycle loop300
                  else
                   ptype_mn(igrid,1)=1.      !rain
                   cycle loop300
                  end if
              else if(cfrzr.ge.crain) then
                   ptype_mn(igrid,2)=1.      !freezing rain
                   cycle loop300
              else
                   ptype_mn(igrid,1)=1.      !rain
                   cycle loop300
              end if

            else if(cslet.gt.cfrzr) then
                   if(cslet.ge.crain) then
                     ptype_mn(igrid,3)=1.      !sleet
                     cycle loop300
                   else
                     ptype_mn(igrid,1)=1.      !rain
                     cycle loop300
                   end if

            else if(cfrzr.ge.crain) then
                   ptype_mn(igrid,2)=1.      !freezing rain
                   cycle loop300
            else
                   ptype_mn(igrid,1)=1.      !rain
                   cycle loop300
            end if

            end if


          end do loop300
200      continue 

        do 400 lv=1,dPlvl(nv)     !dPlvl(nv) always is 1 for precp type
          do 500 lh =1,dTlvl(nv)
            do 600 igrid = 1, jf
   
                 crain = 0.
                 cfrzr = 0.
                 csnow = 0.
                 cslet = 0.

                if(dThrs(nv,lh).eq.1.0) then

                 do irun=1,iens
                   if(miss(irun).eq.0) then
                   crain = crain + prcptype(igrid,irun,1)
                   cfrzr = cfrzr + prcptype(igrid,irun,2)
                   cslet = cslet + prcptype(igrid,irun,3)
                   csnow = csnow + prcptype(igrid,irun,4) 
                  end if
                 end do

! without limits, there are negative values on the individual species. Why?

         pcount=amax1(crain,0.)+amax1(cfrzr,0.)+
     &          amax1(cslet,0.)+amax1(csnow,0.)

                  if (pcount.gt.0) then
                     ptype_pr2(igrid,1)=crain*100./pcount

	if ( ptype_pr2(igrid,1) .gt. 100. ) then
	write(0,*) '1- crain, cfrzr, cslet, csnow, pcount: ', 
     &              crain, cfrzr, cslet, csnow, pcount
	endif
                     ptype_pr2(igrid,2)=cfrzr*100./pcount
	if ( ptype_pr2(igrid,2) .gt. 100. ) then
	write(0,*) '2- crain, cfrzr, cslet, csnow, pcount: ', 
     &              crain, cfrzr, cslet, csnow, pcount
	endif
                     ptype_pr2(igrid,3)=cslet*100./pcount
	if ( ptype_pr2(igrid,3) .gt. 100. ) then
	write(0,*) '3- crain, cfrzr, cslet, csnow, pcount: ', 
     &              crain, cfrzr, cslet, csnow, pcount
	endif
                     ptype_pr2(igrid,4)=csnow*100./pcount
	if ( ptype_pr2(igrid,4) .gt. 100. ) then
	write(0,*) '4- crain, cfrzr, cslet, csnow, pcount: ', 
     &              crain, cfrzr, cslet, csnow, pcount
	endif
                  else
                     ptype_pr2(igrid,1)=0.
                     ptype_pr2(igrid,2)=0.
                     ptype_pr2(igrid,3)=0.
                     ptype_pr2(igrid,4)=0.
                  end if   

                  if (mbrs.gt.0) then
                     ptype_pr(igrid,1)=crain*100./mbrs
                     ptype_pr(igrid,2)=cfrzr*100./mbrs        
                     ptype_pr(igrid,3)=cslet*100./mbrs
                     ptype_pr(igrid,4)=csnow*100./mbrs
                  else
                     ptype_pr(igrid,1)=0.
                     ptype_pr(igrid,2)=0.
                     ptype_pr(igrid,3)=0.
                     ptype_pr(igrid,4)=0.
                  end if   
                end if


             if(igrid.eq.150000) then
              write(*,*)'IN preciptype 150000 ptype_pr'
              write(*,*) (ptype_pr(igrid,kp),kp=1,4)
             end if

             
600      continue
500     continue
400    continue

         write(*,*) 'call preciptype return' 

           return
           end
