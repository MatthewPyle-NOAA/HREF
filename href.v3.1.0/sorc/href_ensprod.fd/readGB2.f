      subroutine readGB2(igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpdx,
     +     gfld,eps,iret)

        use grib_mod

        type(gribfield) :: gfld
        integer jids(200), jpdt(200), jgdt(200)
        integer igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,jpdx
        logical :: unpack=.true.
        character*5 eps   

c        write(*,*)'readGB2 igrb2=',igrb2,jpdtn,jpd1,jpd2,jpd10,jpd12,
c     +          jpdx,eps 

        jids=-9999  !array define center, master/local table, year,month,day, hour, etc, -9999 wildcard to accept any
        jpdt=-9999  !array define Product, to be determined
        jgdt=-9999  !array define Grid , -9999 wildcard to accept any

        jdisc=-1    !discipline#  -1 wildcard 
        jgdtn=-1    !grid template number,    -1 wildcard 
        jskp=0      !Number of fields to be skip, 0 search from beginning

        jpdt(1)=jpd1   !Category #     
        jpdt(2)=jpd2   !Product # under this category     
        jpdt(10)=jpd10 !Product vertical ID      
        !jpdt(27)=jpdx

        if(jpd10.eq.100) then
           jpdt(12)=jpd12*100   !pressure level     
        else 
           jpdt(12)=jpd12
        end if

        if (trim(eps).eq.'sseo') then
          if(jpd1.eq.1.and.(jpd2.eq.8.or.
     +                      jpd2.eq.11) )  then
            jpdt(27)=jpdx
          else
            jpdt(15)=jpdx                !to identify max of 3hr,24hr
          end if
        else
          if (jpdtn.eq.11) then
            jpdt(30)=jpdx                  !jpdtn=11, SREF accum time use jpdt(30)
          else if (jpdtn.eq.8) then 
            jpdt(27)=jpdx 
          end if 
        end if

        if(jpd1.eq.0.and.jpd2.eq.192.and.jpd2.eq.106) then
          jpdt(15)=jpdx
        end if 

        call getgb2(igrb2,0,jskp,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt,
     +     unpack, jskp1, gfld,iret)

        if(iret.ne.0) then
         write(*,*) 'getgb2 error:',iret,' in read ',igrb2,
     +   ' for var ', jpd1,jpd2,jpd10,jpd12
        end if

        return
        end
