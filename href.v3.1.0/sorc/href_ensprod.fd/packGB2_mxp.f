      subroutine packGB2_mxp(imean,imxp,mxp8,
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lq,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,gfld)

        use grib_mod
        include 'parm.inc'

        type(gribfield) :: gfld


         INTEGER,intent(IN) :: imean,imxp,
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lq,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid

        REAL,dimension(jf,Lq,8),intent(IN) :: mxp8

        INTEGER,allocatable,dimension(:) ::   ipdtmpl

        integer ml

c   for max,min,10,25,50,90% mean products
        Character*4 qvname(maxvar)
        Integer qk5(maxvar), qk6(maxvar), qk4(maxvar)
        Character*1 qMsignal(maxvar)
        Integer qMlvl(maxvar), qMeanLevel(maxvar,maxmlvl)

        common /qtbl/nmxp,
     +              qvname,qk4,qk5,qk6,qMlvl,
     +              qMeanLevel,qMsignal


c        write(*,*) 'packing imxp:',imxp

c        write(*,*) nv,imean,imxp, 
c     +     jpd1,jpd2,jpd10,jpd27,jf,Lq,
c     +     iens,iyr,imon,idy,ihr,ifhr,gribid

        DO 1000 ml=1,qMlvl(nv)

            !redefine some of gfld%idsect() array elements 
            gfld%idsect(1)=7
            gfld%idsect(2)=2
            gfld%idsect(3)=0   !experimental, see Table 1.0
            gfld%idsect(4)=1   !experimental, see Table 1.1
            gfld%idsect(5)=1   !experimental, see Table 1.2
            gfld%idsect(6)=iyr  !year
            gfld%idsect(7)=imon !mon
            gfld%idsect(8)=idy  !day
            gfld%idsect(9)=ihr  !cycle time
            gfld%idsect(10)=0
            gfld%idsect(11)=0
            gfld%idsect(12)=0
            gfld%idsect(13)=1


        !Pack direct variable mean/spread

         if (ml.eq.1) then              !this is just do once

            if (jpd1.eq.1.and.(jpd2.eq.8.or.jpd2.eq.11) ) then  
             ipdtnum=12                  !ensemble APCP mean use Template 4.12
             ipdtlen=31
            else 
             ipdtnum=2                   !ensemble NON-accum mean use Template 4.2
             ipdtlen=17
            end if

            allocate (ipdtmpl(ipdtlen))

             ipdtmpl(1)=jpd1
             ipdtmpl(2)=jpd2
             ipdtmpl(3)=4
             ipdtmpl(4)=0
             ipdtmpl(5)=132              !assigned 20161214
             ipdtmpl(6)=0
             ipdtmpl(7)=0
             ipdtmpl(8)=1
             ipdtmpl(9)=ifhr             !fcst time    
             ipdtmpl(10)=jpd10
             ipdtmpl(11)=gfld%ipdtmpl(11)
             !ipdtmpl(12)= see below 
             ipdtmpl(13)=gfld%ipdtmpl(13)
             ipdtmpl(14)=gfld%ipdtmpl(14)
             ipdtmpl(15)=gfld%ipdtmpl(15)
             ipdtmpl(16)=1               !weighted mean of all members
             ipdtmpl(17)=iens            !number of members
         
            if (jpd1.eq.1.and.(jpd2.eq.8.or.jpd2.eq.11) ) then  !Template 4.12 has extra elements than Template 4.2
              !2015121205 correction: B. zhou ...
              ihr_ifhr=ihr+ifhr
              call get_ymd(iyr,imon,idy,ihr_ifhr,kyr,kmon,kdy,khr)
              !ipdtmpl(18)=iyr   !year
              !ipdtmpl(19)=imon  !mon
              !ipdtmpl(20)=idy   !day
              ipdtmpl(18)=kyr   !year
              ipdtmpl(19)=kmon  !mon
              ipdtmpl(20)=kdy   !day
              ipdtmpl(9) =ihr+ifhr-jpd27    !overwrite for APCP: Beginning time of accumulation
              !ipdtmpl(21)=ihr+ifhr         !end time of accumulation
              ipdtmpl(21)=khr         !end time of accumulation
              ipdtmpl(22)=0   
              ipdtmpl(23)=0  
              ipdtmpl(24)=1   
              ipdtmpl(25)=0   
              ipdtmpl(26)=1   
              ipdtmpl(27)=2                 !See Table 4.11, same start time (or same cycle time) 
              ipdtmpl(28)=1  
              ipdtmpl(29)=jpd27 
              ipdtmpl(30)=1
              ipdtmpl(31)=0   
             end if
 
          end if  

          if(jpd10.eq.100) then
             ipdtmpl(12)=qMeanLevel(nv,ml)*100
          else
             ipdtmpl(12)=qMeanLevel(nv,ml)
          end if

          gfld%fld=mxp8(:,ml,imxp) 

          iret=0
          !write(*,*) 'call putgb2'
          !write(*,*) gfld%ipdtnum, gfld%ipdtlen
          !write(*,*) gfld%ipdtmpl
          !call putgb2(imean,gfld,iret)
          !write(*,*) 'call putgb2 done'

          if (gfld%idrtmpl(7).eq.255) gfld%idrtmpl(7)=0

          call Zputgb2(imean,gfld,ipdtmpl,ipdtnum,ipdtlen,iret)

          if(iret.ne.0) then
           write(*,*) 'Zputgb2 mean error:',iret
          end if
          
1000    CONTINUE   


        return
        end
