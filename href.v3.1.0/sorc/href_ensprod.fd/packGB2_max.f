      subroutine packGB2_max(imean,isprd,vrbl_mn,
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,bmap,gfld)

        use grib_mod
        include 'parm.inc'

        type(gribfield) :: gfld

C for variable table:
        Integer numvar
        Character*4 vname(maxvar)
        Integer k5(maxvar), k6(maxvar), k4(maxvar)
        Character*1 Msignal(maxvar), Psignal(maxvar)
        Integer Mlvl(maxvar), MeanLevel(maxvar,maxmlvl)
        Integer Plvl(maxvar), ProbLevel(maxvar,maxplvl)
        Integer Tlvl(maxvar)
        Character*1 op(maxvar)
        Real    Thrs(maxvar,maxtlvl)


        common /tbl/numvar,
     +              vname,k4,k5,k6,Mlvl,Plvl,Tlvl,
     +              MeanLevel,ProbLevel,Thrs,
     +              Msignal,Psignal,op


         INTEGER,intent(IN) :: imean,isprd,
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid

        REAL,dimension(jf,Lm),intent(IN) :: vrbl_mn

        INTEGER,allocatable,dimension(:) ::   ipdtmpl
        LOGICAL*1:: bmap(jf)

        integer ml

        write(*,*) 'packing direct mean/spread for nv ', nv

        write(*,*) nv, imean,isprd,
     +     jpd1,jpd2,jpd10,jpd27,jf,Lm,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid



        DO 1000 ml=1,Mlvl(nv)

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

             ipdtnum=12                 !All use Template 4.12
             ipdtlen=31

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
             ipdtmpl(16)=9               !maximum of all members
             ipdtmpl(17)=iens            !number of members

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

          if(jpd10.eq.100) then
             ipdtmpl(12)=MeanLevel(nv,ml)*100
          else
             ipdtmpl(12)=MeanLevel(nv,ml)
          end if

          gfld%fld=vrbl_mn(:,ml) 

          !!!IMPORTANT 
          gfld%ibmap=0

          iret=0
          !write(*,*) 'ipdtnum, ipdtlen=',ipdtnum, ipdtlen
          !do k=1,ipdtlen
          !  write(*,*) k, ipdtmpl(k)
          !end do

          call Zputgb2(imean,gfld,ipdtmpl,ipdtnum,ipdtlen,bmap,iret)
          if(iret.ne.0) then
           write(*,*) 'Zputgb2 mean error:',iret
          end if
          
1000    CONTINUE   


        return
        end
