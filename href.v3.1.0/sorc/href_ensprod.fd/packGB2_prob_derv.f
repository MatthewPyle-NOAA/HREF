      subroutine packGB2_prob_derv(iprob,derv_pr,
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lt,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid,gfld)

        use grib_mod
        include 'parm.inc'

        type(gribfield) :: gfld

C for variable table:
        integer nderiv
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


         INTEGER,intent(IN) :: iprob,
     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lt,
     +     iens,iyr,imon,idy,ihr,ifhr,gribid

        REAL,dimension(jf,Lp,Lt),intent(IN) :: derv_pr

        INTEGER,allocatable,dimension(:) ::   ipdtmpl

        integer pl

        write(*,*) 'packing derv prob for nv ',nv

c        write(*,*) iprob,
c     +     nv,jpd1,jpd2,jpd10,jpd27,jf,Lp,Lt,
c     +     iens,iyr,imon,idy,ihr,ifhr,gribid

c        write(*,*) 'dPlvl dTlvl=',dPlvl(nv),dTlvl(nv)

        if(jpd10.eq.108) then
         n_dPlvl=dPlvl(nv)-1
        else
         n_dPlvl=dPlvl(nv)
        end if

        DO 2000 pl=1,n_dPlvl

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
            gfld%idrtmpl(6)=0
            gfld%idrtmpl(7)=255


        !Pack direct variable prob

         if (pl.eq.1) then  !this is just do once 

            if (jpd1.eq.1.and.(jpd2.eq.8.or.jpd2.eq.11) ) then  !
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
            !ipdtmpl(12)= see below
             ipdtmpl(13)=gfld%ipdtmpl(13)
             ipdtmpl(14)=gfld%ipdtmpl(14)
             ipdtmpl(15)=gfld%ipdtmpl(15)
             ipdtmpl(16)=0               !???
             ipdtmpl(17)=iens            !number of members

             if(trim(dPsignal(nv)).eq.'A'.or.
     +          trim(dPsignal(nv)).eq.'B'.or.
     +          trim(dPsignal(nv)).eq.'C'.or.
     +          trim(dPsignal(nv)).eq.'D'.or.
     +          trim(dPsignal(nv)).eq.'K'.or.
     +          trim(dPsignal(nv)).eq.'L'.or.
     +          trim(dPsignal(nv)).eq.'E'.or.
     +          trim(dPsignal(nv)).eq.'F'.or.
     +          trim(dPsignal(nv)).eq.'H'.or.
     +          trim(dPsignal(nv)).eq.'I'  ) then
           write(0,*) 'setting derv prob to nbmax type 194'
           ipdtmpl(3)=194     !Neighborhood probability
             endif


            if (jpd1.eq.1.and.(jpd2.eq.8.or.jpd2.eq.11) ) then  !Template 4.9 has extra elements than Template 4.5
              !2015121205 correction: B. zhou ...
              ihr_ifhr=ihr+ifhr
              call get_ymd(iyr,imon,idy,ihr_ifhr,kyr,kmon,kdy,khr)
              ipdtmpl(23)=kyr   !year
              ipdtmpl(24)=kmon  !mon
              ipdtmpl(25)=kdy   !day
              !ipdtmpl(23)=iyr   !year
              !ipdtmpl(24)=imon  !mon 
              !ipdtmpl(25)=idy   !day
              ipdtmpl(9)= ihr+ifhr-jpd27   !overwite for APCP, begin time of accumulation
              !ipdtmpl(26)=ihr+ifhr         !end time of accumulation 
              ipdtmpl(26)=khr              !end time of accumulation 
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
          end if


          if(jpd10.eq.100) then
            ipdtmpl(12)=dProbLevel(nv,pl)*100
          else if (jpd10.eq.108 .and. jpd2 .ne. 1) then
            ipdtmpl(10)=108
            ipdtmpl(11)=0
            ipdtmpl(12)=dProbLevel(nv,pl)*100
            ipdtmpl(13)=108
            ipdtmpl(14)=0
            ipdtmpl(15)=dProbLevel(nv,pl+1)*100
	write(0,*) 'in 108 non wind speed block'
          else if (jpd10.eq.108.and.jpd1.eq.2 .and. jpd2.eq.1) then
            ipdtmpl(10)=100
            ipdtmpl(11)=0
            ipdtmpl(12)=dProbLevel(nv,pl)*100
            ipdtmpl(13)=100
            ipdtmpl(14)=0
            ipdtmpl(15)=dProbLevel(nv,pl+1)*100
	write(0,*) 'in special prob block redefining 108 GRIB packing'
        write(0,*) 'dProbLevel(nv,pl): ', dProbLevel(nv,pl)
        write(0,*) 'dProbLevel(nv,pl+1): ', dProbLevel(nv,pl+1)
          else if (jpd1.eq.3.and.jpd2.eq.5.and.jpd10.eq.101) then !thickness case
             ipdtmpl(10)=100
             ipdtmpl(11)=0
             ipdtmpl(12)=PPairLevel(nv,ml,1)*100
             ipdtmpl(13)=100
             ipdtmpl(14)=0
             ipdtmpl(15)=PPairLevel(nv,ml,2)*100

          else if (jpd2 .eq. 192 .and. jpd10 .eq. 103) then ! bulk shear
             ipdtmpl(10)=103
             ipdtmpl(11)=0
             ipdtmpl(12)=PPairLevel(nv,ml,1)*100
             ipdtmpl(13)=103
             ipdtmpl(14)=0
             ipdtmpl(15)=6000
          else
            ipdtmpl(12)=dprobLevel(nv,pl)
          end if

          do 500 kt=1,dTlvl(nv)
             if (dop(nv).eq.'<') then
              ipdtmpl(18)=0
              ipdtmpl(19)=3
              ipdtmpl(20)=dThrs(nv,kt)*1000  !since scale factor ipdtmpl(19) is 3
              ipdtmpl(21)=3
              ipdtmpl(22)=0

             else if (dop(nv).eq.'>') then
              ipdtmpl(18)=1
              ipdtmpl(19)=0
              ipdtmpl(20)=0
              ipdtmpl(21)=3
              ipdtmpl(22)=dThrs(nv,kt)*1000   !since scale factor ipdtmpl(21) is 3

             else if (dop(nv).eq.'=') then
              ipdtmpl(18)=2
              ipdtmpl(19)=3
              ipdtmpl(20)=dThrs(nv,kt)*1000  !since scale factor ipdtmpl(19) is 3
              ipdtmpl(21)=3
              ipdtmpl(22)=dThrs(nv,kt+1)*1000  !since scale factor ipdtmpl(21) is 3

             else if (dop(nv).eq.'-'.or.dop(nv).eq.'^') then
              ipdtmpl(18)=2
              ipdtmpl(19)=3
              ipdtmpl(20)=dThrs(nv,kt)*1000  !since scale factor ipdtmpl(19) is 3
              ipdtmpl(21)=3
              ipdtmpl(22)=dThrs(nv,kt+1)*1000  !since scale factor ipdtmpl(21) is 3

             end if

              gfld%fld(:)=derv_pr(:,pl,kt)

              if(gfld%ibmap .ne. 255 ) then
                gfld%ibmap = 0   !important resetting
              end if 

              if (gfld%idrtmpl(7).eq.255) gfld%idrtmpl(7)=0

              iret=0
              call Zputgb2(iprob,gfld,ipdtmpl,ipdtnum,ipdtlen,iret)

             if(iret.ne.0) then
              write(*,*) 'Zputgb2 derv prob  error:',iret
             end if


 500      continue

2000    CONTINUE
        return
        end
