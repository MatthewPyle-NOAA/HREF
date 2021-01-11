
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  getlvl function is to get the number 'nn' from 'Mnn/Pnn/Dnn/Hnn'
c   Author: Binbin Zhou
C*  B. Zhou      03/2013 Adapted onto Zeus/wcoss Linux compiler version
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
	function getlvl(S)
 	integer getlvl, Length, temp
        character*3  S, S1
  
        Length=LEN_TRIM(S)
        S1=S(2:Length)
        call ST_NUMB(S1, temp, ier)
         getlvl = temp
        return
        end
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C This function is to get grib id  from region string x 
C such as  x = 'G211/region' or just x = 'G212', 
C Where 212 is grib id
C
C Author Binbin Zhou, 
c        Mar 20, 2005
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        function ID (x)
        character(20) x, gribid
        integer p, length, y
        p = index(x,'/')
        length = len_trim(x)
        if (p.gt.0) then
         gribid = x (2:p-1)
        else
         gribid = x (2:length)
        end if

        call st_numb (gribid, y, ier)
        ID = y

        return
        end 

       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C This function is to get numerical number from a string x, 
C such as x = '100.05', '.5', '100', '100.', etc
C Author: Binbin ZHou
C         Mar 20, 2005
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	function CharToReal(x)     
        real CharToReal
        character(*) x 
        character*20 p1, p2
        integer p, length, d1, d2, pp2
        real r1,r2       

        length = len_trim(x)
        p = index(x,'.')       !e.g.  x=.15: p=1
        if (p.eq.1) then
          p1 = '0'
          p2 = x(p+1:length)
          pp2= len_trim(p2)
        else if (p.gt.1) then   !e.g. x=12.5: p=3
          p1 = x (1:p-1)
          p2 = x(p+1:length)
          pp2= len_trim(p2)
        else
          p1 = x
          p2 = '0'
          pp2 = 0
        end if

        call st_numb (p1,d1,ier)
        call st_numb (p2,d2,ier)

        r1 = float(d2)
        r2 = 0 ! initial definition

        if(pp2.eq.1) then
          r2 = r1 /10.0
        else if (pp2.eq.2) then
          r2 = r1 /100.0
        else if (pp2.eq.3) then
          r2 = r1 /1000.0
        else if (pp2.eq.4) then
          r2 = r1 /10000.0
        else if (pp2.eq.5) then
          r2 = r1 /100000.0
        else if (pp2.eq.6) then
          r2 = r1 /1000000.0
        end if
    
        if(d1 .ge. 0) then 
          CharToReal = d1 + r2
        else
          CharToReal = -1.0 * (abs(d1) + r2)
        end if

         return  
         end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   function index_table() returns the numvar index of pair (kpds5, kpds6) 
c   last appearing in the direct variable sequence from the table file
c
c   Author: Binbin Zhou
c           Aug. 4, 2005
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	function index_table(k5,k6,n5,n6, n)
          integer k5(*), k6(*), n5, n6, n
          index_table = 0
          do i = 1, n
           if(k5(i).eq.n5.and. k6(i).eq.n6 ) then
             index_table = i
             return
           end if
          end do
         return
         end
       

          function index_table_var(var,k5,k6,n5,n6,name,n)
          integer k5(*), k6(*), n5, n6, n
          character*4 var(*), name

          index_table_var = 0
          do i = 1, n
           if(k5(i).eq.n5.and. k6(i).eq.n6
     +     .and. trim(var(i)).eq.trim(name) ) then
             index_table_var = i
             return
           end if
          end do
         return
         end

 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    function index_int_array() returns the last first location
c    of element in the array x 
c
c    Author: Binbin Zhou
c           Aug. 4, 2005
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
        function index_int_array(x, element, n)
          integer x(*), element, n
          index_int_array = 0
          do i = 1, n
            if(x(i).eq.element) then
              index_int_array = i
              return
            end if
          end do
         return
         end
          

          

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine readtbl() is to read table file in which variables and derived 
c       variables are defined
c    Original Author: Binbin Zhou
c
c    05-08-02  B. Zhou original code was built    
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine readtbl(nunit) 

        INCLUDE 'parm.inc'
 
        Integer numvar, nsubstr, nm, ier, lng, getlvl
        Character*200 oneline
        Character*20  substr(50)

c    for variables
        Character*4 vname(maxvar)
        Integer k5(maxvar), k6(maxvar), k4(maxvar)
        Character*1 Msignal(maxvar), Psignal(maxvar)
        Character*3 Mn(maxvar), Pn(maxvar), Tn(maxvar)
        Integer Mlvl(maxvar), MeanLevel(maxvar,maxmlvl)
        Integer Plvl(maxvar), ProbLevel(maxvar,maxplvl)
        Character*1 op(maxvar)
        Integer Tlvl(maxvar)
        Character*20 Cthrs
        Real      Thrs(maxvar,maxtlvl)

c    for derived variables
        Character*4 dvname(maxvar)
        Integer dk5(maxvar), dk6(maxvar),dk4(maxvar)
        Character*1 dMsignal(maxvar), dPsignal(maxvar)
        Character*3 dMn(maxvar), dPn(maxvar), dTn(maxvar)
        Integer dMlvl(maxvar), dMeanLevel(maxvar,maxmlvl)
        Integer dPlvl(maxvar), dProbLevel(maxvar,maxplvl)
        Character*1 dop(maxvar)
        Integer dTlvl(maxvar)
        Real    dThrs(maxvar,maxtlvl)
        Integer MPairLevel(maxvar,maxmlvl,2)
        Integer PPairLevel(maxvar,maxplvl,2)

c   for max,min,10,25,50,90% mean products
        Character*4 qvname(maxvar)
        Integer qk5(maxvar), qk6(maxvar),qk4(maxvar)
        Character*1 qMsignal(maxvar)
        Character*3 qMn(maxvar)
        Integer qMlvl(maxvar), qMeanLevel(maxvar,maxmlvl)

        Integer missvar(maxvar,30)         !dynamically build a missing array for direct variable missing in each member
        Character Tsignal(maxvar)          !2015-12-09: new added to account on Gaussian smoothing sigma value   
        Character dTsignal(maxvar)          !2018-10-02:  new for derived Gauss smooth  


        common /tbl/numvar,
     +              vname,k4,k5,k6,Mlvl,Plvl,Tlvl,
     +              MeanLevel,ProbLevel,Thrs,
     +              Msignal,Psignal,op 

        common /dtbl/nderiv,
     +              dvname,dk4,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop


        common /qtbl/nmxp,
     +              qvname,qk4,qk5,qk6,qMlvl,
     +              qMeanLevel,qMsignal

        common /Tsgn/Tsignal ,dTsignal

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       First, read direct variable information
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        read(nunit,*) numvar  
        n1 = 1 
        n2 = numvar

c        write(*,*) 'numvar=',numvar

        do 1000 n = n1, n2       

         read (nunit, '(A)') oneline
c         write(*,*) oneline
         call ST_CLST( oneline, ' ', ' ', 50, substr, nsubstr, ier ) !note: if too many levels, need to increase 50
         lng=len_trim(substr(1))
         vname(n)=substr(1)(1:lng)

         call ST_NUMB( substr(2), k4(n), ier)
         call ST_NUMB( substr(3), k5(n), ier)
         call ST_NUMB( substr(4), k6(n), ier)

         lng=len_trim(substr(5))
         Mn(n)=substr(5)(1:lng)

         Msignal(n)=Mn(n)(1:1)
         Mlvl(n) = getlvl(Mn(n))

         if ( Mlvl(n) .gt. 0 ) then         !Mean level > 0 
           do i = 1, Mlvl(n)
             call ST_NUMB(substr(5+i), MeanLevel(n,i), ier)
           end do
         end if

         nm = 5 + Mlvl(n)

c         write(*,*) vname(n), k4(n),k5(n),k6(n),Mn(n),Msignal(n),
c     +       Mlvl(n),(MeanLevel(n,i),i=1,Mlvl(n))

         if ( nsubstr .gt. nm) then      !there are Prob requests    
  
           lng=len_trim(substr(nm+1))
           Pn(n)=substr(nm+1)(1:lng)

           Psignal(n)=Pn(n)(1:1)
           Plvl(n) = getlvl(Pn(n))
          
           do i = 1, Plvl(n)
             call ST_NUMB(substr(nm+1+i), ProbLevel(n,i), ier)
           end do

           np=nm+1+Plvl(n)         
           OP(n)=substr(np+1)

           lng=len_trim(substr(np+2))
           Tn(n)=substr(np+2)(1:lng)

           Tsignal(n)=Tn(n)(1:1)

           Tlvl(n) = getlvl(Tn(n))

           do i = 1, Tlvl(n)
             lng=len_trim(substr(np+2+i))
             Cthrs=substr(np+2+i)(1:lng)

             Thrs(n,i)=CharToReal(Cthrs) 
           end do

c         write(*,*) Pn(n),Psignal(n),Plvl(n),
c     +    (ProbLevel(n,i),i=1,Plvl(n)),
c     +    OP(n),Tn(n),Tlvl(n),(Thrs(n,i),i=1,Tlvl(n))

         end if

1000     continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          Then, read derived variable information
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         read(nunit,*,END=3000) nderiv

c         n1 = numvar + 1
c         n2 = numvar + nderiv
 
          n1 = 1
	  n2 = nderiv
          
       do 1100 n = n1, n2

         read (nunit, '(A)') oneline

         call ST_CLST( oneline, ' ', ' ', 50, substr, nsubstr, ier )
         lng=len_trim(substr(1))
         dvname(n)=substr(1)(1:lng)

         call ST_NUMB( substr(2), dk4(n), ier)
         call ST_NUMB( substr(3), dk5(n), ier)
         call ST_NUMB( substr(4), dk6(n), ier)
                                                                                                                                                                                                             
         lng=len_trim(substr(5))
         dMn(n)=substr(5)(1:lng)

         dMsignal(n)=dMn(n)(1:1)
         dMlvl(n) = getlvl(dMn(n))
                                                                                                                                                                                                             
         if ( dMlvl(n) .gt. 0 ) then         !Mean level > 0
           do i = 1, dMlvl(n)
             call ST_NUMB(substr(5+i), dMeanLevel(n,i), ier)
           end do
         end if
              
          
         if( dk6(n).eq.101) then      !Thickness, then need to read additional
c          write(*,*) 'substr=',(substr(i),i=1,nsubstr)
c          write(*,*) 'nsubstr=',nsubstr, dMlvl(n)

          do i = 1, dMlvl(n)
            call ST_NUMB(substr(5+dMlvl(n)+i*2-1),MPairLevel(n,i,1),ier)
            call ST_NUMB(substr(5+dMlvl(n)+i*2),  MPairLevel(n,i,2),ier)
c            write(*,*)'MPairLevel=',MPairLevel(n,i,1),MPairLevel(n,i,2)
          end do


          nm = 5 + 3*dMlvl(n)  
 
         else
                                                                                                                             
          nm = 5 + dMlvl(n)
 
         end if     
                                                                                                                                                                                                   
         if ( nsubstr .gt. nm) then      !there are Prob requests
                                                                                                                                                                                                             
           lng=len_trim(substr(nm+1))
           dPn(n)=substr(nm+1)(1:lng)

           dPsignal(n)=dPn(n)(1:1)
           dPlvl(n) = getlvl(dPn(n))
                                                                                                                                                                                                             
           do i = 1, dPlvl(n)
             call ST_NUMB(substr(nm+1+i), dProbLevel(n,i), ier)
           end do
 
          if(dk6(n).eq.101 .or.dk6(n).eq.104 .or.dk6(n).eq.106   !then need to read additional
     +     .or.dk6(n).eq.112 .or.dk6(n).eq.114 .or.            !level pairs
     +         dk6(n).eq.121. or.dk6(n).eq.128 .or.
     +         dk6(n).eq.141. or.dk6(n).eq.116) then

            do i = 1, dPlvl(n)
             call ST_NUMB(substr(nm+1+dPlvl(n)+i*2-1),
     +                              PPairLevel(n,i,1),ier)
             call ST_NUMB(substr(nm+1+dPlvl(n)+i*2),  
     +                              PPairLevel(n,i,2),ier)
            end do
            np = nm + 1 + 3*dPlvl(n)
          else
            np=nm + 1 + dPlvl(n)
          end if 

           dOP(n)=substr(np+1)

           lng=len_trim(substr(np+2))
           dTn(n)=substr(np+2)(1:lng)
           dTsignal(n)=dTn(n)(1:1)
           write(0,*) 'pulled dTsignal(n): ',n, dTsignal(n)

           dTlvl(n) = getlvl(dTn(n))
                                                                                                                                                                                                             
           do i = 1, dTlvl(n)
             lng=len_trim(substr(np+2+i))
             Cthrs=substr(np+2+i)(1:lng)

             dThrs(n,i)=CharToReal(Cthrs) 
c             write(*,*) "Cthrs, dThrs(n,i)=",Cthrs, dThrs(n,i)
           end do
               
         end if
 
1100    continue  

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Then read max,max,10,25,50 and 50% mean  products request (MXP)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        read(nunit,*,END=3000) nmxp

        n1 = 1
        n2 = nmxp

        do 1200 n = n1, n2

         read (nunit, '(A)') oneline

         call ST_CLST( oneline, ' ', ' ', 50, substr, nsubstr, ier ) !note: if too many levels, need to increase 50
         lng=len_trim(substr(1))
         qvname(n)=substr(1)(1:lng)

         call ST_NUMB( substr(2), qk4(n), ier)
         call ST_NUMB( substr(3), qk5(n), ier)
         call ST_NUMB( substr(4), qk6(n), ier)

         lng=len_trim(substr(5))
         qMn(n)=substr(5)(1:lng)

         qMsignal(n)=qMn(n)(1:1)
         qMlvl(n) = getlvl(qMn(n))

         if ( qMlvl(n) .gt. 0 ) then         ! qMean level > 0
           do i = 1, qMlvl(n)
             call ST_NUMB(substr(5+i), qMeanLevel(n,i), ier)
           end do
         end if

1200    continue


  
3000    return
          end         
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c getporb: compute probability accoding to the given thresholds 
c Author: Binbin Zhou
c  Aug. 4, 2005
c
c  input: x     - data
c         thrs  - threshold
c         opr    - operation (>, < , =, - )
c
c  output: prob
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine getprob(x,n,thrs1,thrs2,opr,prob,miss,wgt)
	real x(*), wgt(*), thrs1, thrs2, prob
        integer n,miss(*)
        character*1 opr
        real count


         wsum=0.0
         do i=1,n
           if(miss(i).eq.0) then
            wsum=wsum+wgt(i)
           end if 
         end do

        count = 0.0

        do i = 1, n

         if (miss(i).eq.0) then
CMP       if (trim(opr).eq.'>') then
          if (opr(1:1).eq.'>') then
            if (x(i).ge.thrs1)count=count+wgt(i)/wsum
CMP       else if (trim(opr).eq.'<') then
          else if (opr(1:1).eq.'<') then
            if (x(i).le.thrs1) count=count+wgt(i)/wsum
CMP       else if (trim(opr).eq.'=') then
          else if (opr(1:1).eq.'=') then
            if (x(i).eq.thrs1) count=count+wgt(i)/wsum
CMP       else if (trim(opr).eq.'-') then
          else if (opr(1:1).eq.'-') then
            if (x(i).ge.thrs1.and.x(i).lt.thrs2) 
     +       count = count + wgt(i)/wsum
          end if
         end if
        end do

         prob = 100.0 * count

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   subroutine getmean: compute mean and spread of one dimension array  
c   Author: Binbin Zhou
c   Aug. 3, 2005
c   Apr. 7 2009: Zhou B. Add weight (wgt) for VSREF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine getmean (x,n,mean,spread,miss,wgt)
	 real x(*), wgt(*), mean, spread
         integer n,miss(*)
      

         wsum=0.0
         do i=1,n
           if(miss(i).eq.0) then
             wsum=wsum+wgt(i)
           end if
         end do
 
         if(wsum.eq.0.0) then         !if All members are missing
           mean = -9999.0
           spread=0.
           return
         end if

         mean = 0.
         do i=1,n
           if(miss(i).eq.0) then
            mean = mean + x(i)*(wgt(i)/wsum)
           end if
         end do

         spread = 0.
         do i = 1, n
           if(miss(i).eq.0) then
            spread = spread + (wgt(i)/wsum)*(x(i)-mean)**2
           end if
         end do

         spread = sqrt (spread )

         return
         end
                          
      
        subroutine getmean_fog (x, n,  mean, spread, miss,wgt)
         real x(*), wgt(*), mean, spread
         integer n,miss(*)


         wsum=0.0
         do i=1,n
          if(miss(i).eq.0) then          
           if(x(i).gt.0.0) then
            wsum=wsum+wgt(i)
           end if
          end if
         end do

         mean = 0.
         do i=1,n
          if(miss(i).eq.0) then
           if(wsum.gt.0.0) then
             mean = mean + x(i)*(wgt(i)/wsum)
           end if
          end if
         end do

         if(mean.eq.0.0) then
          spread = 0.0
          return
         end if


cvsref         mean = mean / n

         spread = 0.
         do i = 1, n
          if(miss(i).eq.0) then
           spread = spread + (wgt(i)/wsum)*(x(i)-mean)**2
          end if
         end do

         spread = sqrt (spread )

         return
         end

c -------------------------------
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   subroutine get_cond_mean: compute conditional mean and spread of one dimension array
c      Conditional mean is the mean under some condition, if equal to certain
c      very large value, then it is not taken into mean, but spread still
c      count it (only "less than" case)
c
c   Author: Binbin Zhou
c   Aug. 3, 2005
c   Dec. 15, 2008: Use dominant numbers as mean
c   Apr. 7, 2009: Zhou B. Add weight (wgt) for VSREF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                                                     
        subroutine get_cond_mean (x,n,alarge,mean,spread,miss,wgt)
         real x(*),wgt(30), mean, spread, alarge,count, half
         integer n,miss(*)
                                                                                                     
         mean = 0.
         count = 0.

         wsum=0.0
         count = 0.
       loop1:       do i=1,n
           if(miss(i).eq.0) then
             if(x(i).ge.alarge) cycle loop1
             wsum=wsum+wgt(i)
             count = count + 1.0
           end if
         end do loop1

         half=count/2.0
    
        loop2:       do  i=1,n
          if (miss(i).eq.0) then
           if(x(i).ge.alarge) cycle loop2 
            mean = mean + x(i)*(wgt(i)/wsum)
          end if
           end do loop2
                     
         if( count .le. half ) then           !only most of member happen   
           mean = alarge
         end if 
                                                                                            
         spread = 0.
         loop3:      do i = 1, n                         
          if(x(i).ge.alarge.or.miss(i).eq.1) cycle loop3                       !!! modified in Apr. 9 2008 
           spread = spread + (wgt(i)/wsum)*(x(i)-mean)**2
         end do loop3
            
         if (mean.eq.alarge) then
             spread=0.0              !!! for VSREF, assume this. Apr. 22, 2009
         else 
             spread = sqrt (spread )
         end if
       
         return
         end

c -------------------------------
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   subroutine get_cond_lpm_mean: compute conditional mean and spread of one dimension array
c      Conditional mean is the mean under some condition, if equal to certain
c      very large value, then it is not taken into mean, but spread still
c      count it (only "less than" case)
c
c   Author: Binbin Zhou
c   Aug. 3, 2005
c   Dec. 15, 2008: Use dominant numbers as mean
c   Apr. 7, 2009: Zhou B. Add weight (wgt) for VSREF
c   Jan 8, 2021: M. Pyle create an LPM mean version
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                                                     
        subroutine get_cond_lpm_mean (x,n,alarge,mean,spread,miss,wgt)
         real x(*),wgt(30), mean, spread, alarge,count, half
         integer n,miss(*)
                                                                                                     
         mean = 0.
         count = 0.

         wsum=0.0
         count = 0.
       loop1:       do i=1,n
           if(miss(i).eq.0) then
             if(x(i).ge.alarge) cycle loop1
             wsum=wsum+wgt(i)
             count = count + 1.0
           end if
         end do loop1

         half=count/2.0
    
        loop2:       do  i=1,n
          if (miss(i).eq.0) then
           if(x(i).ge.alarge) cycle loop2 
            mean = mean + x(i)*(wgt(i)/wsum)
          end if
           end do loop2
                     
         if( count .le. half ) then           !only most of member happen   
           mean = alarge
         end if 
                                                                                            
         spread = 0.
         loop3:      do i = 1, n                         
          if(x(i).ge.alarge.or.miss(i).eq.1) cycle loop3                       !!! modified in Apr. 9 2008 
           spread = spread + (wgt(i)/wsum)*(x(i)-mean)**2
         end do loop3
            
         if (mean.eq.alarge) then
             spread=0.0              !!! for VSREF, assume this. Apr. 22, 2009
         else 
             spread = sqrt (spread )
         end if
       
         return
         end



c -------------------------------

        subroutine get_cond_mean_test(igrid,x,n,alarge,mean,spread,
     +         miss,wgt)
         real x(*),wgt(30), mean, spread, alarge,count, half
         integer n,miss(*)

         mean = 0.
         count = 0.

         wsum=0.0
         count = 0.
          loop102:       do i=1,n
          if(miss(i).eq.0) then
           if(x(i).ge.alarge) cycle loop102
           wsum=wsum+wgt(i)
           count = count + 1.0
          end if
          end do loop102
         half=count/2.0
 
c         if(igrid.eq.737601) write(*,*)'half=',half,wsum
  
         loop100:      do i=1,n
           if(miss(i).eq.0) then
            if(x(i).ge.alarge) cycle loop100
            mean = mean + x(i)*(wgt(i)/wsum)
c            if(igrid.eq.737601) write(*,*)i,'mean=',mean
           end if
         end do loop100
 
         if( count .le. half ) then           !only most of member happen

           mean = alarge
         end if

c         if(igrid.eq.737601) write(*,*)'mean=',mean


         spread = 0.
        loop200:      do i = 1, n 
          if(x(i).ge.alarge.or.miss(i).eq.1) cycle loop200                       !!! modified in Apr. 9 2008
           spread = spread + (wgt(i)/wsum)*(x(i)-mean)**2
c           if(igrid.eq.737601) write(*,*)i,'spread=',spread
          end do loop200

         if (mean.eq.alarge) then
             spread=0.0              !!! for VSREF, assume this. Apr. 22, 2009
         else
             spread = sqrt (spread )
c             if(igrid.eq.737601) write(*,*)'spread=',spread
         end if

         return
         end

        subroutine get_cond_mean_lwc(x,n,alarge,mean,spread,miss)
         real x(*), mean, spread, alarge,count
         integer n,miss(*)
                                                                                                                                                                                       
         mean = 0.
         count = 0.
                                                                                                                                                                                       
       loop_lwc:      do i=1,n
           if(x(i).eq.alarge.and.miss(i).eq.0) cycle loop_lwc
            mean = mean + x(i)
            count = count + 1.0
        end do loop_lwc
                                                                                                                                                                                       
         if(count .gt. 0.0) then
          mean = mean / count
         else
          mean = 0.0
         end if
                    
         spread = 0.
         count = 0.
         do i = 1, n
           if(x(i).gt.alarge.and.miss(i).eq.0) then
            count = count + 1.0
            spread = spread + (x(i)-mean)**2
           end if
         end do
          
         if(count .gt. 0.0) then          
            spread = sqrt (spread / count )
         else 
            spread = 0.0
         end if
           
         return
         end

 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   subroutine getwindmean: compute mean and spread of wind vector
c   Author: Binbin Zhou
c   March 1, 2006
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                                                                                                                          
        subroutine getwindmean (u,v,n, mean, spread,miss,wgt)
         real u(*),v(*), Umean, Vmean, Uspread,Vspread,
     +             mean, spread,count,wgt(*)         
         integer n,miss(*) 
                                                                                                                                          
         Umean = 0.
         Vmean = 0.
         mean = 0.
         count=0.0

         wsum=0.0
         do i=1,n
          if(miss(i).eq.0) then
c           Umean = Umean + u(i)
c           Vmean = Vmean + v(i)
c           mean = mean + sqrt(v(i)*v(i)+u(i)*u(i))
c           count=count+1.0
           wsum=wsum+wgt(i)          
          end if
         end do
        

         do i=1,n
           if(miss(i).eq.0) then
            Umean = Umean + u(i)*(wgt(i)/wsum)
            Vmean = Vmean + v(i)*(wgt(i)/wsum)
            mean = mean + sqrt(v(i)*v(i)+u(i)*u(i))
     +        *(wgt(i)/wsum)
           end if
         end do
     
         Uspread = 0.
         Vspread = 0.
         spread = 0.
         do i = 1, n
          if(miss(i).eq.0) then
           Uspread=Uspread+(wgt(i)/wsum)*(u(i)-Umean)**2
           Vspread=Vspread+(wgt(i)/wsum)*(v(i)-Vmean)**2
          end if
         end do
 
         Uspread=sqrt(Uspread)     
         Vspread=sqrt(Vspread)     
         spread = sqrt(Uspread*Uspread + Vspread*Vspread)

         return 
         end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   subroutine getmax: compute max values among members
c   Author: Binbin Zhou
c   Nov. 12, 2015
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine getmax (x,n,vmax,miss)
         real x(*), vmax
         integer n,miss(*)

         vmax=0.
         do i=1,n
           if(miss(i).eq.0) then
            if(x(i).gt.vmax) vmax=x(i)
           end if
         end do

         return
         end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   subroutine getSP: find if exceeding certain threshold to generate
c    Spaghetti plot according to SPC SSEO ensemble SP products 
c   Author: Binbin Zhou
c   Nov. 12, 2015
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine getSP(x,n,thrs1,thrs2,opr,sp,miss)
        real x(n),  thrs1, thrs2,sp,nk(n)
        integer n,miss(n)
        character*1 opr

        nk=1.

        do i = 1, n
         if (miss(i).eq.0) then
          if (opr.eq.'>') then
            if (x(i).ge.thrs1) nk(i) = 2.
          else if (opr.eq.'<') then
            if (x(i).le.thrs1) nk(i) = 2.
          else if (opr.eq.'=') then
            if (x(i).eq.thrs1) nk(i) = 2.
          else if (opr.eq.'-') then
            if (x(i).ge.thrs1.and.x(i).lt.thrs2)
     +       nk(i) = 2.
          end if
         end if
        end do

        sp=0.0
        do i =1,n
         sp= nk(i)*10.0**(n-i) + sp
        end do

        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccc      
c  getDigit: purpose: get a number "1" or "2" from 
c   string sp (e.g. 2212111) for each ensemble member
c   "2" means member exceed the threshold, "1" not exceed
c   To build spaghetti plots for SSEO 
c
c   input: 
c          sp: a number(e.g. 2202211.0) 
c           n: ensemble size 
c  output: 
c         d(n): a nember to show if a member is exceed threshold
c         "2"--> 100
c         "1" --> 0
c
c  B. Zhou: Dec 10, 2015
c
        subroutine getDigit (sp,d,n)
         real sp, d(n)
         integer k(n)

         c = 0.
         do i=n-1,0,-1
          b=10.**i
          k(n-i)= int( (sp - c)/b )
          c = c + k(n-i) * b
         !write(*,*) i, b, k(n-i), c  
         end do
          d(:) = k(:)*1.0

         do i=1,n
          if(d(i).eq.2.0) then 
            d(i)=100.0
          else  
            d(i)=0.
          end if
         end do

        return
        end
 

        subroutine get_joint(x,y,opr,t1,t2,jindex)
         real x,y,t1,t2,jindex
         character opr

         jindex=0.0
         if (opr.eq.'>') then
           if(x.ge.t1.and.y.ge.t2) jindex=1.0
         else if (opr.eq.'<') then
           if(x.le.t1.and.y.le.t2) jindex=1.0
         else if (opr.eq.'=') then
           if(x.eq.t1.and.y.eq.t2) jindex=1.0
         else if (opr.eq.'^') then
           if(x.ge.t1.and.y.le.t2) jindex=1.0
         else
           write(*,*) opr, ' operator not set up'
         end if

         return
         end

c  Subroutine Gsmoothing 
c   Purpose smoothing a gridded field using Gaussian Kernel smooting technique 
c
c   Reference: Amy Harless et al: "A report and feature-based verification stduy
c       of the CAPS 2008 storm-scale ensemble forecast for severe convection weather"
c       AMS Conference 2012 
c   Input: 
c          A: gridded field to be smoothed
c      im,jm: X and Y dimension of field A  
c        nbr: range of smoothing (in grids) 
c          s: sigma value of Gaussian smoothing
c   Output: 
c          A: Smoothed field 
c   Programmer: 
c         2015-12-02: Binbin Zhou, NCEP/EMC  
c
c
        subroutine Gsmoothing (A,jf,im,jm,ns,gs)
         real A(jf), B(jf)
         character ns,gs
         integer s,nbr

         if(ns.eq.'A'.or.ns.eq.'E'.or.ns.eq.'H') then
           nbr=5
         else if (ns.eq.'J') then
           nbr=3
         else if (ns.eq.'K') then
           nbr=4
         else if (ns.eq.'L') then
           nbr=2
         else if (ns.eq.'M') then
           nbr=1
         else if (ns.eq.'B'.or.ns.eq.'F'.or.ns.eq.'I') then
           nbr=10
         else if (ns.eq.'C') then
           nbr=15
         else if (ns.eq.'D') then
           nbr=20
         else
           nbr=1
         end if

         if(gs.eq.'A'.or.gs.eq.'E'.or.gs.eq.'H') then
           s=5
         else if (gs.eq.'K') then
           s=3
         else if (ns.eq.'L') then
           s=2
         else if (ns.eq.'M') then
           s=1
         else if (gs.eq.'B'.or.gs.eq.'F'.or.gs.eq.'I') then
           s=10
         else if (gs.eq.'C') then
           s=15
         else if (gs.eq.'D') then
           s=20
         else
           s=1
         end if

         write(*,*) 'In Gsmoothing,nbr,s=',nbr,s 

         f1=1./(3.14*s*s)
         f2=-0.5/(s*s)
         B=0.0

         do jp = 1,jm
          do ip = 1,im

             ijp=(jp-1)*im + ip
             i1 = ip - nbr
             i2 = ip + nbr
             j1 = jp - nbr
             j2 = jp + nbr
             if ( i1.le.1 ) i1=1
             if ( j1.le.1 ) j1=1
             if ( i2.ge.im ) i2=im
             if ( j2.ge.jm ) j2=jm

             Gsum=0.
             AxG=0.
             do j = j1,j2
              do i = i1,i2
                ij = (j-1)*im + i
                !if (A(ij).gt.0.0) then
                G=f1*exp(f2*((ip-i)*(ip-i)+(jp-j)*(jp-j)))
                Gsum=Gsum+G
                AxG=AxG+G*A(ij)  
                !end if
              end do
             end do

            if(Gsum.gt.0.0) then
             B(ijp)=AxG/Gsum
            end if

           end do
          end do

          A=B             

          return
          end


      FUNCTION LEAP (YEAR) RESULT (LEAPFLAG)

      IMPLICIT NONE

      INTEGER :: YEAR
      LOGICAL :: LEAPFLAG

      LEAPFLAG = .FALSE.
      IF (MOD(YEAR,4) .EQ. 0)   LEAPFLAG = .TRUE.
      IF (MOD(YEAR,100) .EQ. 0) LEAPFLAG = .FALSE.
      IF (MOD(YEAR,400) .EQ. 0) LEAPFLAG = .TRUE.
      RETURN
      END
 
      subroutine get_ymd (y,m,d,h, yy,mm,dd,hh)
      integer y,m,d,h, yy,mm,dd,hh
      integer md
      logical leap 

      yy=y
      mm=m
      dd=d
      hh=h

       if (m.eq.1.or.m.eq.3.or.m.eq.5.or.m.eq.7.or.
     +     m.eq.8.or.m.eq.10.or.m.eq.12) then
         md=31
       else if (m.eq.4.or.m.eq.6.or.m.eq.9.or.m.eq.11) then
         md=30
       else
         if(leap(y).eq..true.) then 
           md=29
         else
           md=28
         end if
       end if 
 
       if(h.ge.24) then
         hh=mod(h,24)
         kd=int(h/24)
         dd=d+kd
         if (dd.gt.md) then
           dd=dd-md
           mm=m+1
           if (mm.gt.12) then
            mm=1
            yy=y+1
           end if
         end if
       end if

       return
       end 
       

c ------------------------------------------------

c   subroutine neighborhood_max. Purpose:
c     Get neiborhood max value within a circle with radius nbr 
c     Currently, 4 radius are set: 5,10,15,20 (grids), depending
c           on user-defined in the product table 
c   Input:
c          A: gridded field 
c      im,jm: X and Y dimension of field A
c          s: sign (A,or B, or C, or D)  
c   Output:
c          A: modified gridded field
c   Programmer:
c         2015-12-09: Binbin Zhou, NCEP/EMC
c

       subroutine neighborhood_max (A,jf,im,jm,s)
         real A(jf), Amax(jf)
         real nbr,dist 
         character s 

	 write(0,*) 'in neighborhood_max with sign s: ', s

	write(0,*) 'shape(A): ', shape(A)
        write(0,*) 'jf, im, jm: ', jf, im,jm

         if(s.eq.'A') then
           nbr=8.  ! to mimic the 40 km radius
         else if (s.eq.'K') then
           nbr=4. ! try small neighborhood for mountain snow
         else if (s.eq.'J') then
           nbr=3. ! try small neighborhood for ceil/vis
         else if (s.eq.'L') then
           nbr=2. ! try very small neighborhood for lightning
         else if (s.eq.'M') then
           nbr=1. ! try very small neighborhood for lightning
         else if (s.eq.'B') then
           nbr=10.
         else if (s.eq.'C') then
           nbr=15.
         else if (s.eq.'D') then
           nbr=20.
         else 
           nbr=1.
         end if

        write(*,*) 'In neighborhood: nbr=', nbr

         do jp = 1,jm
          do ip = 1,im

           ijp=(jp-1)*im + ip

           i1 = ip - int(nbr)
           i2 = ip + int(nbr)
           j1 = jp - int(nbr)
           j2 = jp + int(nbr)
           if ( i1.le.1 ) i1=1
           if ( j1.le.1 ) j1=1
           if ( i2.ge.im ) i2=im
           if ( j2.ge.jm ) j2=jm 
           !The searching radius-nbr circle is just within the 2nbx2nbr squre:
           ! i2-i1 x j2-j1 

           Amax(ijp)=A(ijp)
           !Search max value within the circle (radius=nbr with center grid ip,jp) 
           do j = j1,j2
            do i = i1,i2
              ij = (j-1)*im + i
              dist=sqrt((ip-i)*(ip-i)+1.*(jp-j)*(jp-j)) 
              if(A(ij).gt.Amax(ijp).and.dist.le.nbr) Amax(ijp)=A(ij)
            end do
           end do

          end do
         end do

         A=Amax

       return
       end

c ------------------------------------------------

c   subroutine neighborhood_min. Purpose:
c     Get neighborhood min value within a circle with radius nbr 
c     Currently, 4 radius are set: 5,10,15,20 (grids), depending
c           on user-defined in the product table 
c   Input:
c          A: gridded field 
c      im,jm: X and Y dimension of field A
c          s: sign (A,or B, or C, or D)  
c   Output:
c          A: modified gridded field
c   Programmer:
c         2015-12-09: Binbin Zhou, NCEP/EMC
c         2021-01-08: Matthew Pyle, NCEP/EMC
c

       subroutine neighborhood_min (A,jf,im,jm,s)
         real A(jf), Amin(jf)
         real nbr,dist 
         character s 

	 write(0,*) 'in neighborhood_min with sign s: ', s

	write(0,*) 'shape(A): ', shape(A)
        write(0,*) 'jf, im, jm: ', jf, im,jm

         if(s.eq.'A') then
           nbr=8.  ! to mimic the 40 km radius
         else if (s.eq.'J') then
           nbr=3. ! try small neighborhood for ceil/vis
         else if (s.eq.'K') then
           nbr=4. ! try small neighborhood for mountain snow
         else if (s.eq.'L') then
           nbr=2. ! try very small neighborhood for lightning
         else if (s.eq.'M') then
           nbr=1. ! try very small neighborhood for lightning
         else if (s.eq.'B') then
           nbr=10.
         else if (s.eq.'C') then
           nbr=15.
         else if (s.eq.'D') then
           nbr=20.
         else 
           nbr=1.
         end if

        write(*,*) 'In neighborhood: nbr=', nbr

         do jp = 1,jm
          do ip = 1,im

           ijp=(jp-1)*im + ip

           i1 = ip - int(nbr)
           i2 = ip + int(nbr)
           j1 = jp - int(nbr)
           j2 = jp + int(nbr)
           if ( i1.le.1 ) i1=1
           if ( j1.le.1 ) j1=1
           if ( i2.ge.im ) i2=im
           if ( j2.ge.jm ) j2=jm 
           !The searching radius-nbr circle is just within the 2nbx2nbr squre:
           ! i2-i1 x j2-j1 

           Amin(ijp)=A(ijp)
           !Search max value within the circle (radius=nbr with center grid ip,jp) 
           do j = j1,j2
            do i = i1,i2
              ij = (j-1)*im + i
              dist=sqrt((ip-i)*(ip-i)+1.*(jp-j)*(jp-j)) 
              if(A(ij).lt.Amin(ijp).and.dist.le.nbr) Amin(ijp)=A(ij)
            end do
           end do

          end do
         end do

         A=Amin

       return
       end


c ------------------------------------------------

c   subroutine neighborhood_fraction. Purpose:
c     Get neiborhood fraction  within a circle with radius nbr
c     Currently, 2 radius are set: 5,10 (grids), depending
c           on user-defined in the product table
c   Input:
c          A: gridded field
c      im,jm: X and Y dimension of field A
c          s: sign (E or F) 
c      t1,t2: thresholds
c         op: operator >, or < or = 
c
c   Output:
c          A: modified gridded field that become fraction 
c   Programmer:
c         2016-01-06: Binbin Zhou, NCEP/EMC
c

       subroutine neighborhood_fraction (A,jf,im,jm,s,t1,t2,op)
         real A(jf), Afrc(jf)
         real nbr,dist,t1,t2
         character s,op

        
         if(s.eq.'E') then
           nbr=5.
         else if (s.eq.'F') then
           nbr=10.
         else
           nbr=1.
         end if

c        write(*,*) 'In neighborhood_fraction: nbr=', nbr
c        write(*,*) jf,im,jm,s,t1,t2,op

         Afrc=0.
         do jp = 1,jm
          do ip = 1,im

           ijp=(jp-1)*im + ip

           i1 = ip - int(nbr)
           i2 = ip + int(nbr)
           j1 = jp - int(nbr)
           j2 = jp + int(nbr)
           if ( i1.le.1 ) i1=1
           if ( j1.le.1 ) j1=1
           if ( i2.ge.im ) i2=im
           if ( j2.ge.jm ) j2=jm
           !The searching radius-nbr circle is just within the 2nbx2nbr squre:
           ! i2-i1 x j2-j1 

           !get fraction within the circle (radius=nbr with center grid
           !ip,jp)
           Np=0
           do j = j1,j2
            do i = i1,i2
              ij = (j-1)*im + i
              Np=Np+1
              dist=sqrt((ip-i)*(ip-i)+1.*(jp-j)*(jp-j))
              if (dist.le.nbr.and.op.eq.'>') then
               if (A(ij).ge.t1) Afrc(ijp)=Afrc(ijp)+1.
              else if (dist.le.nbr.and.op.eq.'<') then
               if (A(ij).le.t1) Afrc(ijp)=Afrc(ijp)+1. 
              else if (dist.le.nbr.and.op.eq.'=') then
               if (A(ij).ge.t1.and.A(ij).le.t2) Afrc(ijp)=Afrc(ijp)+1.
              end if
            end do
           end do
c             if(ijp.eq.398918) write(*,*) Afrc(ijp), Np
             Afrc(ijp)=Afrc(ijp)/Np
          end do
         end do
     
         
c          write(*,*) A(398918), Afrc(398918)

          A=Afrc

       return
       end

c ------------------------------------------------

c   subroutine neighborhood_average. Purpose:
c     Get neiborhood average within a circle with radius nbr
c     Currently, 2 radius are set: 5,10 (grids), depending
c           on user-defined in the product table
c   Input:
c          A: gridded field
c      im,jm: X and Y dimension of field A
c          s: sign (E or F)
c      t1,t2: thresholds
c         op: operator >, or < or =
c
c   Output:
c          A: modified gridded field that become fraction
c   Programmer:
c         2016-01-06: Binbin Zhou, NCEP/EMC
c

       subroutine neighborhood_average (A,jf,im,jm,s)
         real A(jf), Avg(jf)
         real nbr,dist
         character s


         if(s.eq.'H') then
           nbr=5.
         else if (s.eq.'I') then
           nbr=10.
         else
           nbr=1.
         end if

c        write(*,*) 'In neighborhood_avg: nbr=', nbr
c        write(*,*) jf,im,jm,s

         Avg=0.
         do jp = 1,jm
          do ip = 1,im

           ijp=(jp-1)*im + ip

           i1 = ip - int(nbr)
           i2 = ip + int(nbr)
           j1 = jp - int(nbr)
           j2 = jp + int(nbr)
           if ( i1.le.1 ) i1=1
           if ( j1.le.1 ) j1=1
           if ( i2.ge.im ) i2=im
           if ( j2.ge.jm ) j2=jm
           !The searching radius-nbr circle is just within the 2nbx2nbr squre:
           ! i2-i1 x j2-j1 

           !get fraction within the circle (radius=nbr with center grid
           !ip,jp)
           Np=0
           do j = j1,j2
            do i = i1,i2
              ij = (j-1)*im + i
              dist=sqrt((ip-i)*(ip-i)+1.*(jp-j)*(jp-j))
              if (dist.le.nbr) then
               NP=NP+1
               Avg(ijp)=Avg(ijp) + A(ij)
              end if
            end do
           end do
             Avg(ijp)=Avg(ijp)/Np
          end do
         end do


c          write(*,*) A(398918), Avg(398918)

          A=Avg

       return
       end

