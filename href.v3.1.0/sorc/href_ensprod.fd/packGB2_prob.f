c  SUBROUTINE packGB2_prob: Write probability data into GRIB2 file
c
      subroutine packGB2_prob(iprob,vrbl_pr,
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lt,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,gfld)


        use grib_mod
        include 'parm.inc'

        type(gribfield) :: gfld

C for variable table:
        Integer numvar
        Character*4 vname(maxvar)
        Integer k5(maxvar), k6(maxvar), k4(maxvar)
        Character*1 Msignal(maxvar), Psignal(maxvar)
        Integer Mlvl(maxvar), MeanLevel(maxvar,maxplvl)
        Integer Plvl(maxvar), ProbLevel(maxvar,maxplvl)
        Integer Tlvl(maxvar)
        Character*1 op(maxvar)
        Real    Thrs(maxvar,maxtlvl)


        common /tbl/numvar,
     +              vname,k4,k5,k6,Mlvl,Plvl,Tlvl,
     +              MeanLevel,ProbLevel,Thrs,
     +              Msignal,Psignal,op


         INTEGER,intent(IN) :: iprob,
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lt,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid

        REAL,dimension(jf,Lp,Lt),intent(IN) :: vrbl_pr

        INTEGER,allocatable,dimension(:) ::   ipdtmpl

        integer pl

        write(*,*) 'packing direct prob for nv ',nv
	write(0,*) 'HERE WITH IYR: ', iyr
	write(0,*) 'here with idy: ', idy
	write(0,*) 'here with ihr: ', ihr

c        write(*,*) iprob, 
c     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lt,
c     +     iens,iyr,imon,idy,ihr,ifhr,gribid

c        write(*,*) 'Plvl Tlvl=',Plvl(nv),Tlvl(nv)


        DO 2000 pl=1,Plvl(nv)


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

            !redefine data representation, otherwise some (e.g. VIS)
            !prob could not be stored in (Min=0.0, Max=0.0)
            gfld%idrtmpl(1)=0
            gfld%idrtmpl(2)=-1
            gfld%idrtmpl(3)=0
            gfld%idrtmpl(4)=8
            gfld%idrtmpl(5)=0
!	if (gfld%IDRTNUM .gt. 5) then
            gfld%idrtmpl(6)=0
!	endif
!	if (gfld%IDRTNUM .gt. 6) then
            gfld%idrtmpl(7)=255
!        endif


        !Pack direct variable prob

         if (pl.eq.1) then  !this is just do once 

            if (jpd1.eq.1.and.
     &         (jpd2.eq.8.or.jpd2.eq.11 .or. jpd2.eq.13 .or.
     &          jpd2.eq.194) ) then   ! add 194 for PPFFG
             ipdtnum=9                  !ensemble APCP prob use Template 4.12
             ipdtlen=36
            else
             ipdtnum=5                   !ensemble NON-accum prob use Template 4.5
             ipdtlen=22
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
             ipdtmpl(9)=ifhr            !Forecast time
             ipdtmpl(10)=jpd10
             ipdtmpl(11)=gfld%ipdtmpl(11)
             ipdtmpl(13)=gfld%ipdtmpl(13)
             !ipdtmpl(12)= see below
             ipdtmpl(14)=gfld%ipdtmpl(14)
             ipdtmpl(15)=gfld%ipdtmpl(15)
             ipdtmpl(16)=0               !???
             ipdtmpl(17)=iens            !number of members


             if(trim(Psignal(nv)).eq.'A'.or.
     +          trim(Psignal(nv)).eq.'B'.or.
     +          trim(Psignal(nv)).eq.'C'.or.
     +          trim(Psignal(nv)).eq.'D'.or.
     +          trim(Psignal(nv)).eq.'K'.or.
     +          trim(Psignal(nv)).eq.'E'.or.
     +          trim(Psignal(nv)).eq.'F'.or.
     +          trim(Psignal(nv)).eq.'H'.or.
     +          trim(Psignal(nv)).eq.'I'  ) ipdtmpl(3)=194     !Neighborhood probability


            if (jpd1.eq.1.and.
     &         (jpd2.eq.8.or.jpd2.eq.11 .or. jpd2.eq.13  .or.
     &          jpd2.eq.194) ) then  !Template 4.9 has extra elements than Template 4.5
              !2015121205 correction: B. zhou ...
              !ihr_ifhr=ihr+ifhr
              !call get_ymd(iyr,imon,idy,ihr_ifhr,kyr,kmon,kdy,khr)
              !ipdtmpl(23)=kyr   !year
              !ipdtmpl(24)=kmon  !mon
              !ipdtmpl(25)=kdy   !day
              !ipdtmpl(9)= ihr+ifhr-jpd27   !overwite for APCP, begin time of accumulation
              !ipdtmpl(26)=khr              !end time of accumulation

              call get_time_GB2(iyr,imon,idy,ihr,ifhr,
     +            iyr1,imon1,idy1,ihr1)

              ipdtmpl(23)=iyr1   !year
              ipdtmpl(24)=imon1  !mon
              ipdtmpl(25)=idy1   !day
              ipdtmpl(9)= ifhr-jpd27   !overwite for APCP, begin time of accumulation
              ipdtmpl(26)=ihr1         !end of fcst time of accumulation

              ipdtmpl(27)=0
              ipdtmpl(28)=0
              ipdtmpl(29)=1
              ipdtmpl(30)=0
              ipdtmpl(31)=1                !See Table 4.11, same start time (or same cycle time) 
              ipdtmpl(32)=2
              ipdtmpl(33)=1
              ipdtmpl(34)=jpd27
              ipdtmpl(35)=1
              ipdtmpl(36)=0
             end if

          end if  ! p1=1?


          if(jpd10.eq.100) then
            ipdtmpl(12)=ProbLevel(nv,pl)*100
          else
            ipdtmpl(12)=probLevel(nv,pl)
          end if

          do 500 kt=1,Tlvl(nv)
             if (op(nv).eq.'<') then
              ipdtmpl(18)=0
              ipdtmpl(19)=3
              ipdtmpl(20)=Thrs(nv,kt)*1000  !since scale factor ipdtmpl(19) is 3
              ipdtmpl(21)=3
              ipdtmpl(22)=0

             else if (op(nv).eq.'>') then
              ipdtmpl(18)=1
              ipdtmpl(19)=0
              ipdtmpl(20)=0
              ipdtmpl(21)=3
              ipdtmpl(22)=Thrs(nv,kt)*1000   !since scale factor ipdtmpl(21) is 3
	if (Psignal(nv).eq.'B'.or.Psignal(nv).eq.'C'.or.
     &      Psignal(nv).eq.'D') then
              ipdtmpl(22)=(Thrs(nv,kt)+1.)*1000   !since scale factor ipdtmpl(21) is 3
	write(*,*) 'using alt Thrs val: ', Thrs(nv,kt)+1
       endif

             else if (op(nv).eq.'=') then
              ipdtmpl(18)=2
              ipdtmpl(19)=3
              ipdtmpl(20)=Thrs(nv,kt)*1000  !since scale factor ipdtmpl(19) is 3
              ipdtmpl(21)=3
              ipdtmpl(22)=Thrs(nv,kt+1)*1000  !since scale factor ipdtmpl(21) is 3
             else if (op(nv).eq.'-') then
              ipdtmpl(18)=2
              ipdtmpl(19)=0
              ipdtmpl(20)=Thrs(nv,kt)  !since scale factor ipdtmpl(19) is 0
              ipdtmpl(21)=0
              ipdtmpl(22)=Thrs(nv,kt+1)  !since scale factor ipdtmpl(21) is 0

             end if

              gfld%fld(:)=vrbl_pr(:,pl,kt)
 
              if ( gfld%ibmap .ne. 255 ) then
               gfld%ibmap=0      !important resetting
              end if 

              if (gfld%idrtmpl(7).eq.255) gfld%idrtmpl(7)=0

              iret=0
	
!	if (jpd2 .eq. 13) then
!	write(0,*) 'switch type from WEASD to SNOD for prob'
!	ipdtmpl(2)=11
!	endif

              call Zputgb2(iprob,gfld,ipdtmpl,ipdtnum,ipdtlen,iret)

             if(iret.ne.0) then
              write(*,*) 'Zputgb2 prob  error:',iret
             end if


 500      continue

2000    CONTINUE
        return
        end
