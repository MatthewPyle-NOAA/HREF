c        Test
c        integer yr,mm,dy,fhr,cyc,yr1,mm1,dy1,hr
c         yr=2015
c         mm=4
c         dy=5
c         fhr=12
c         cyc=6
c         call get_time_GB2(yr,mm,dy,cyc,fhr,
c     +            yr1,mm1,dy1,hr)
c         write(*,*) yr,mm,dy,cyc,fhr
c         write(*,*) yr1,mm1,dy1,hr
c         stop
c         end          


        subroutine get_time_GB2(yr,mm,dy,cyc,fhr,
     +            yr1,mm1,dy1,hr)

         integer:: yr,mm,dy,fhr,cyc,yr1,mm1,dy1,hr,t,n,h
         integer:: ndy, days_in_month

           t=cyc+fhr
           n=int(t/24) 
           h=mod (t,24)  
           yr1=yr
           mm1=mm
           if ( n .lt. 1) then
             dy1=dy
             hr=t
           else 
             dy1=dy+n
             ndy=days_in_month(yr,mm)
             hr=h 
             if (dy1.gt.ndy) then
              dy1=dy1-ndy
              mm1=mm+1
               if (mm1.gt.12) then
                 mm1=1
                 yr1=yr+1
               end if
              end if
           end if
            
          return
          end

          function days_in_month(yr,mm)
          integer days_in_month,yr,mm
          select case (mm)
           case (1,3,5,7,8,10,12)  
            days_in_month = 31
           case (2)
            days_in_month  = 28
             if ( mod(yr,4) == 0) days_in_month = 29  !!Leap year 
           case default 
             days_in_month = 30
          end select

          return
          end              
       
            
