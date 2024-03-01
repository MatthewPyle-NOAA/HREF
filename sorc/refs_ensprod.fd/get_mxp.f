cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine get_mxp: compute max,min,mode, 10,25,50, and 90% mean products
c     Author: Binbin Zhou, Apr. 7, 2011, based on Jun Du code
c     May 24, 2013, Binbin Zhou: Modified to grib2 I/O
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine get_mxp (nv,ifunit,jpdtn,jf,iens,Lq,mxp8,missing,wgt) 

          use grib_mod
          include 'parm.inc'

c   for max,min,10,25,50,90% mean products
        Character*4 qvname(maxvar)
        Integer qk5(maxvar), qk6(maxvar), qk4(maxvar)
        Character*1 qMsignal(maxvar)
        Integer qMlvl(maxvar), qMeanLevel(maxvar,maxmlvl)


        common /qtbl/nmxp,
     +              qvname,qk4,qk5,qk6,qMlvl,
     +              qMeanLevel,qMsignal


        INTEGER, intent(IN) :: nv, jf, iens
        REAL,dimension(jf,Lq,8),intent(INOUT) :: mxp8
        REAL,dimension(jf,iens) :: var 
        INTEGER ifunit(iens) 

        real x(iens)
        INTEGER,dimension(jf,iens),intent(IN) :: missing
        integer miss(iens)
        real wgt(30)
        type(gribfield) :: gfld


c        write(*,*) 'In mxp'
c        write(*,*)  nv,ifunit,jf,iens,Lq 

         !jpdtn=0
         jd1=qk4(nv)
         jd2=qk5(nv)
         jd10=qk6(nv)
         jp27=-9999

         do 1000 lv=1,qMlvl(nv)                       !for all levels

           do irun=1,iens

             if(qk4(nv).eq.1.and.qk5(nv).eq.8) then
              write(*,*) 'get APCP GRIB2 data for member ', irun
             else
              jd12=qMeanLevel(nv,lv)
              call readGB2(ifunit(irun),jpdtn,jd1,jd2,jd10,jd12,
     +             jp27,gfld,eps,ie)
              var(:,irun)=gfld%fld(:)
             end if

            end do
 
            do 500 igrid = 1,jf
              x(:)=var(igrid,:) 
              call mean_median(iens,x,xm,xmed)
              call probability(x,iens,vmin,vmax,p10,p25,p50,p75,p90)
              mxp8(igrid,lv,1)=vmin
              mxp8(igrid,lv,2)=vmax
              mxp8(igrid,lv,3)=3*p50-2*xm

              if((jd1.eq.1.and.jd2.eq.0).or.(jd1.eq.1.and.jd2.eq.1)
     +           .or.(jd1.eq.1.and.jd2.eq.8)) then
c                if(kv5.eq.61.or.kv5.eq.51.or.kv5.eq.52) then

                 if (mxp8(igrid,lv,3).lt.0.0)        then    !Jun: to deal with negative mod
                     mxp8(igrid,lv,3)=p50
                 end if
               end if

               mxp8(igrid,lv,4)=p10
               mxp8(igrid,lv,5)=p25
               mxp8(igrid,lv,6)=p50
               mxp8(igrid,lv,7)=p75
               mxp8(igrid,lv,8)=p90

500        continue

1000      continue          

          return
          end
