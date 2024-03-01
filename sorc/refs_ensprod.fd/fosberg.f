c      PROGRAM TEST
c      REAL T2,W10,RH2, FWI
c      write(*,*) 'input T2,W10,RH2'
c      read(*,*) T2,W10,RH2
c      CALL  FOSBERG_INDEX(T2,W10,RH2, FWI)
c      write(*,*) T2,W10,RH2, FWI
c      stop
c      end



        SUBROUTINE FOSBERG_INDEX(T2,W10,RH2,FWI)
c
c   FOSBERG_INDEX is to calculate Fosberg Fire Weather Index (FWI,1987) based on 
c    Temperature(T2, 2m in F), U10,V10 (10m U V in m/s), and RH2 (2m in %)
c    Its output FWI ranges between 0.0 ~ 100.0
c    First program: B. ZHOU, 12/10/2013
c

       real, intent (IN) :: T2,W10,RH2
       real, intent (OUT) :: FWI
       real A, B, T, W, RH   

       T=T2 * (9./5.) + 32.        !from C -> F
       W=W10*3.6/1.6               !from m/s -> miles/hr
       RH=RH2

 
       if (RH .lt. 10.0 ) then
         A=0.03229 + 0.281073*RH - 0.000578*RH*T
       else if (RH .ge. 10.0 .and. RH .lt. 50.0) then
         A=2.22749 + 0.160107*RH - 0.01478*T
       else
         A=21.0606 + 0.005565*RH*RH - 0.00035*RH*T - 0.483199*RH
       end if

       B=1-2*(A/30)+1.5*(A/30)**2-0.5*(A/30)**3
       FWI=B*sqrt(1+W*W)/0.3002
       FWI=MIN(FWI,100.0)
       FWI=MAX(FWI,0.0) 
     
       RETURN
       END
      
