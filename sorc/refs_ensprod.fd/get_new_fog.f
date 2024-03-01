ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  SUBROUTINE get_new_fog: implemented based on Zhou and Ferrier 's aysmptotic formulation
cc      2008 on JAMC but add advection term
cc
cc     Input: grid-wide u,v component profile and surface values, temperature and RH profiles
cc            and surface values, surface height, and previous time step temperature profile
cc            and surface values (are used to compute cooling rate).
cc     Output: fog LWC
cc
cc
cc    This code is for grid or point-wide computation as following steps:
cc
cc         (1) Search for which (pressure)levels the surface height is located
cc         (2) Search for which (pressure)levels the suturated top above the surface is located
cc                and get the fog layer (saturated layer) depth, if its depth > 800m, not
cc                a fog but reture, otherwise
cc         (3) Search from the fog layer top downward to see if its bottom reach the ground
cc             if its bottom not touch the ground and > 800m, it is not a fog but return,
cc             otherwise
cc         (4) Compute averaged cooling rate within the fog layer using current and previous
cc             time step's temperature profiles and surface (in C/hr)
cc         (5) Based on temperaure at nearest level from the fog top and 2m to compute
cc             averag temperature with the fog layer, based on which parameter beta is computed
cc       (5.1) From cooling rate, fog depth and beta, the critical turb exch. ceoff Kc can
cc               be computed
cc         (6) Based on the u,v at nearest level from the fog top and 10m u,v, to compute
cc             wind speeds at the fog top and the surface
cc         (7) Based on T, wind speeds at nearest level from the fog top and 2m, Ri is computed
cc             from which stability function Fm(Ri) is computed, from which turbulent excahnge
cc             coefficient K is computed
cc       (7.1) If won't further compute LWC, using K and Kc, fog exists or not can be determined
cc              at this step (K>Kc, K<Kc?), otherwise,
cc         (8) With K, cooling rate, beta, fog boundary layer FBL can be computed
cc         (9) From the FBL, cooling rate, fog depth, LWC profile with the fog layer are computed
cc
cc    Language: Fortran 77/90
cc
cc    Origninal author: Binbin Zhou, EMC/NCEP/NOAA
cc                      August 28, 2010
cc
cc    Modifcation history:
cc        Sept 9, 2010 by B. Zhou:
cc          Add startus cloud layer checking (< 400m), if yes, use modeled cloud layer as fog layer
cc               instead of using RH to check saturated layer as fog layer. This imply the modeled cloud
cc               depth is assumed correct but just adjust its LWC distrbution
cc               If clout top < 400m and cloud bottm < 50 m above the ground,
cc               find whcih level the cloud top is located, then repeat the above steps from
cc               step (4)
cc
cc        Oct. 10, 2010 by B. Zhou:
cc          Add LWC advection term from RH or specific humidity field into the Kc and LWC
cc               computations
cc
cc
cc
cc   Description diagram:
cc
cc             kp and kp+1: the level just below and abobe the surcae
cc                    kfog: fog top level
cc
cc   (1) Clear sky case
cc
cc  ---------------------------------------------------------------------------- k=14, hp(14)
cc
cc   .....................................................                       .....
cc
cc  ------------ fog top---100%-------------- kfog (=4 in this case)------------ k=4, hp(4)
cc         ^                             ^
cc         |                             |
cc   ------|-------------- 100%----------|--- kp+1 (=3 in this case)------------ k=3, hp(3)
cc      depth_fog=z_RH-hsfc              |
cc         | ............ o-o (10m)      |
cc         |   ^  ....[2m] |            z_RH=hp(kfog)
cc      ___v___|____^___Y__|_____________|___________________________ surface level (hsfc)
cc     /       |    |                    |                           \
cc-------------|----|--------------------|---kp (=2 in this case) --------------- k=2, hp(2)                                                \
cc   /         |   z2m=2+hsfc            |                             \
cc  /          |    |                    |                              \
cc /        z10m=10+hsfc                 |                               \
cc ------------|----|--------------------|-----------------------------------------k=1, hp(1)
cc/            |    |                    |                                 \
cc ____________v____v____________________v__________________________________\______sea level
cc|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc
cc  (2) Stratus cloud case (cloud top: cldt<400m)
cc
cc  ~~~~~~~cloud top (cldt < 400m)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cc         ^                             ^
cc         |                             |
cc   -kfog-|- (level just below cldt) ---|---------------------------------------k=4. hp(4)
cc         |                             |
cc         |                             |
cc      depth_fog=cldt (i.e. z_RH-hsfc)  |
cc         |                             |
cc  -------|-----------------------------|----------------------------------------k=3, hp(3)
cc         |                             |
cc         | ............ o-o (10m)      |
cc         |   ^  ....[2m] |            z_RH=cldt+hsfc
cc      ___v___|____^___Y__|_____________|___________________________ surface level (hsfc)
cc     /       |    |                    |                           \
cc-------------|----|--------------------|---kp (=2 in this case) --------------- k=2, hp(2)                                                \
cc   /         |   z2m=2+hsfc            |                             \
cc  /          |    |                    |                              \
cc /        z10m=10+hsfc                 |                               \
cc ------------|----|--------------------|-----------------------------------------k=1, hp(1)
cc/            |    |                    |                                 \
cc ____________v____v____________________v__________________________________\______sea level
cc||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine get_new_fog(igrid,up,vp,hp,u10,v10,hsfc,
     +            rhp,rh2,t2,tp,dt2,dtp,interval,lv,cldt,cldb,
     +            adv,rhthr,lwc)
c     +              fog_depth,cooling_rate,lwc_max,dwdz,
c     +       lamda,FBL,Km,bx,cx)

         parameter (Rd = 287.0, Rv=416.0, Ro=1.2)

c                    Input
c         up,vp : u,v on pressure levels
c            hp : height of pressure levels
c       u10,v10 : 10m u and v
c          hsfc : sfc height
c            tp : current temperature on  pressure levels
c           rhp : RH on pressure levels
c            t2 : current temperature at 2m
c           dt2 : 2m temperature change
c           dtp : temperature change on pressure levels
c           rh2 : previous and current RH at 2m
c       interval: time step interval (in hours)    
c           cldb: cloud base
c           cldt: cloud top
c           adv : advection term
c         rhthr : RH threshold deponding on regions
c                    Output
c           lwc : fog LWC
c
c
c  Modification 20130304: Binbin ZHou: To save space, use temperature 
c           change and current temperature instead of store temperature
c           in all forecast hours. Also use only current LWCfield instead
c           of store LWC in all forecast hours
c           


        INTEGER lv,interval
        REAL up(lv),vp(lv),hp(lv),rhp(lv),tp(lv)
        REAL u10,v10,hsfc,lwc,rh2,t2,cooling_rate,
     +       Km,lamda_1,lamda,lwc_max,Kc,cldt,cldb,
     +       dt2,dtp(lv)

        real rhthr
        real cooling(lv+1)
        character*15 ftype
        integer maybefog

        !write(*,*) 'in get_new_fog '
c        if(igrid.eq.1) then
c            write(*,*) igrid,u10,v10,t2,rh2,hsfc,
c     +      cldt,cldb,adv,rhthr,dt,dt2
c           do k=1,14
c             write(*,'(6f9.2)') up(k),vp(k),hp(k),rhp(k),
c     +            tp(k),dtp(k)
c           end do
c         end if

        NFOG1=110591     !Wenatchee WA (EAT Airport)
        NFOG2=52101      !Sherman TX
        NFOG3=75679      !Dulles
        NFOG4=67443      !Chanute, KS
        NFOG5=103040     !BISMARK, ND
        NFOG6=83368      !JFK
        NFOG7=82440      !University park Airport, PA
        NFOG8=84178      !Chicago, IL
        NFOG9=75834      !San Fracisco, CA
        NFOG10=39464     !Austin, TX

         ftype='FOG0'
         lwc=0.0
         cooling=-9.9     !20130306: add it to fix bug

         if (igrid.eq.NFOG1
     +    .or.igrid.eq.NFOG2
     +    .or.igrid.eq.NFOG3
     +    .or.igrid.eq.NFOG4
     +    .or.igrid.eq.NFOG5
     +    .or.igrid.eq.NFOG6
     +    .or.igrid.eq.NFOG7
     +    .or.igrid.eq.NFOG8
     +    .or.igrid.eq.NFOG9
     +    .or.igrid.eq.NFOG10) then
          write(*,'(i7,a6,f6.2,a6,10f6.2,a6,f6.2)') igrid, 
     +    '  RH2=',rh2, 'RHp=',(rhp(i),i=1,10), ' rhthr=',rhthr 
         end if

c        rhthr = 98.0   !assumened saturated 
c        rhthr = 90.0   !it seems that 95% is too large to capture some fogs.

c   search sfc in which pressure layers     
         z10m = 10. + hsfc
         z2m  =  2. + hsfc

         if(cldt.le.cldb) cldt=cldb

          if (z10m .lt. hp(1) ) then
            kp = 0
          else
            do k = 1,lv-1
             if(z10m.ge.hp(k).and.z10m.lt.hp(k+1)) then
               kp = k
             end if
            end do
          end if



          if(cldt.lt.400.0.and.cldb.lt.50.0) then    !yes there is stratus cloud case
c          if(cldt.lt.800.0.and.cldb.lt.50.0) then    !yes there is stratus cloud case
            kfog=0
            ftype='SFC_STRA'
            z_RH=cldt+hsfc
            !seach which level the cloud top is below
            do k = 1, lv-1
             if(z_RH.ge.hp(k).and.z_RH.lt.hp(k+1)) then
               kfog = k+1
             end if
            end do
             
            if(rh2.ge.rhthr) then  !new add 20130306 searhc saturated layer inside stratus
               maybefog=1
               goto 1000
            else                   
              do k=1,kfog
               if(rhp(k).ge.rhthr) then
                 maybefog=1
                 goto 1000
               end if
              end do
              return
            end if  
  
           goto 1000
          end if

c   CCCCCCCCCC following is non-stratus but fog situation: CCCCCCCCCCCCCCCCCCCCC

c         hp(kp+1) : the pressure level just above the surface level 
c         hp(kfog) : the pressure level just above the fog top          

c   search RH=100% (rhthr) depth above the sfc and averaged cooling rate
c          if(cldb.lt.2000.or.cldt.lt.5000.) then
c            rhthr=99.5           !20130207: under stratus, using higher threshold 
c          end if 

          z_RH=hsfc  
          kfog=0
          do k = 1, lv-1
            if(rhp(k).ge.rhthr.and.
     +         rhp(k+1).lt.rhthr ) then
              kfog=k           !find which pressure level is saturated
              z_RH = hp(kfog)  
            end if
          end do

         if((z_RH-hsfc).gt.800.0) then     !DEEP STRATUS'
           !if (igrid.eq.NFOG2) then
           !  write(*,*)'              DEP_STRA',
     +     !    'z_RH-hsfc=',z_RH-hsfc 
           !end if  
           return                      !not a fog layer but a stratus cloud
         end if

          !now search from fog top downward to find fog bottom reach the ground?
          kbottom=0
          if(kfog.ge.2) then
            do k = kfog,1,-1
             if (rhp(k).lt.rhthr) then
              kbottom=k+1
              goto 100
             end if
            end do
           else
             kbottom=1
             goto 100
          end if

100       continue
          if(rh2.ge.rhthr) then
            touch_down=0.
          else
            touch_down=hp(kbottom)-hsfc
          end if


        if(touch_down.gt.100.0) then
         if (igrid.eq.NFOG1
     +    .or.igrid.eq.NFOG2
     +    .or.igrid.eq.NFOG3
     +    .or.igrid.eq.NFOG4
     +    .or.igrid.eq.NFOG5
     +    .or.igrid.eq.NFOG6
     +    .or.igrid.eq.NFOG7
     +    .or.igrid.eq.NFOG8
     +    .or.igrid.eq.NFOG9
     +    .or.igrid.eq.NFOG10) then

            write(*,17)'                       NO_TROUCH_DOWN STRA',
     +      ' z_RH-hsfc=',z_RH-hsfc,'touch_down=',touch_down 
17          format(a41,a11,f8.1,a11,f8.1)
          end if
            return               !not a fog layer but a stratus cloud
         end if
           if(cldt.lt.5000.) then
             ftype='BLW_STRA'
           else
             ftype='PURE_FOG'
           end if   

1000      continue

          if(z_RH.le.z2m) then
            if(rh2.ge.rhthr) then 
              depth_fog = 2.0
              cooling_rate=-dt2/interval   !in C/hr 
              tfog=t2
            else
              depth_fog = 0.0
              return
            end if 
          else
             depth_fog = z_RH - hsfc
c             cooling_rate= -dtp(kfog)- dt2/interval/2.0
c             tfog=(tp(kfog)+t2)/2.0

            !Jan 25, 2013: Binbin Zhou:
            !Some users (East region WFO, Geoff Manikin etc)       
            !complain fog prob and visibility<1km are not
            !consistent. Fogs are not there but visibility shows < 1km
            !Seems the above under-estimate the cooling rate 
            !Modified as search for  maximum cooling below the fog depth:
            
            do k=kp+1,kfog+1
c test            do k=kp+1,kfog
             cooling(k)=-dtp(k)/interval
            end do
            cooling(kfog+2)=-dt2/interval  !just stored 2m cooling there

            cooling_rate=maxval (cooling)

c            if(igrid.eq.NFOG5) then
c         write(*,14) '     crate-profile=',(cooling(m),m=1,kfog+2)
c14       format(a20,14f7.3)
c         write(*,'(a12,f7.3)') 'max-cooling=',cooling_rate
c            end if

c            cooling_rate=0.
c            nn=0
c            do k=kp+1,kfog+2
c             cooling_rate=cooling_rate+cooling(k)
c             nn=nn+1
c            end do
c            cooling_rate=cooling_rate/nn

            tfog=0.
            nn=0
            do k=kp+1,kfog
             tfog=tfog+tp(k)
             nn=nn+1
            end do
             tfog=tfog+t2
             tfog=tfog/(nn+1)

          end if

         if (igrid.eq.NFOG1
     +    .or.igrid.eq.NFOG2
     +    .or.igrid.eq.NFOG3
     +    .or.igrid.eq.NFOG4
     +    .or.igrid.eq.NFOG5
     +    .or.igrid.eq.NFOG6
     +    .or.igrid.eq.NFOG7
     +    .or.igrid.eq.NFOG8
     +    .or.igrid.eq.NFOG9
     +    .or.igrid.eq.NFOG10) then

            write(*,15) ftype,'cldtp=',cldt,'cldbs=',cldb,
     +      '  fogdpth=',depth_fog, 'crate=',cooling_rate,'adv=',adv
15         format(a38,a8,f8.1,a8,f8.1,a10,f8.1,a8,f7.3,a5,d11.4)
          end if


c          tt=7.45*tfog/(235.0 + tfog)
C        Use WMO recommended formulation (2008)
C   Guide to Metteorological Instruments and Methods of
C   Observation (CIMO Guide):

          if( tfog.ge.0.0) then
            tt=17.62*tfog/(243.12 + tfog)     !over liquid water
          else
            tt=22.46*tfog/(272.62 + tfog)     !over ice, assume no supercooled water
          end if
 
          es = 6.112 *2.71828 ** tt
          Tk=tfog + 273.17
          beta = 622.0* 2.5e6 * es / (Rv*Tk*Tk*1000)

          c_adv=beta*cooling_rate/3600.0 + adv/1000.0

          if(c_adv.lt.0.0) then
            lwc=0.0
            return
          end if

          lwc_max = sqrt(c_adv*depth_fog/0.062)


c   compute FBL (Fog Boundary Layer, delta)
          if ( depth_fog .le. 2.0 ) then
            wp=sqrt(up(kp+1)*up(kp+1)+vp(kp+1)*vp(kp+1))
            w10=sqrt(u10*u10+v10*v10)
            dwdz=(wp-w10)/(hp(kp+1)-z10m)
            dtdz=(tp(kp+1)-t2)/(hp(kp+1)-z2m)
          else
            wp=sqrt(up(kfog)*up(kfog)+vp(kfog)*vp(kfog))
            w10=sqrt(u10*u10+v10*v10)
            dwdz=(wp-w10)/(hp(kfog)-z10m)
            dtdz=(tp(kfog)-t2)/(hp(kfog)-z2m)
          end if
 
          dwdz=abs(dwdz)    
          if(dwdz.eq.0.0) dwdz=0.001      


          Ri = (9.8/Tk)*(dtdz+0.01)/dwdz/dwdz 
          if (Ri.le.0.0) then                            !Unstable case: Beljaars: "The parameterization
           z=hp(kp+1)-hsf                                !of the planetary boundary layer", May 1992, ECMWF report
           y=(1.0+z/0.02)                                !You can use other formulation
           Ch=12.0*sqrt(y)/(log(y)*log(y))
           Fm=1.0-10.0*Ri/(1.+Ch*sqrt(-Ri))
          else                                           !Stable case: Holtslag, et al. BLM 2006, vol. 118
c           Fm=1.0/(1.0+10.0*Ri)                         !take long tail format, this leads to higher turbulence
                                                         !and fog LWC becomes smaller or dispersed 
                                                         !You can use other formulation
            if(Ri.ge.0.0.and.Ri.lt.0.1) then             !Sharp tail (Ric ~ 0.25)
             Fm=(1.0-5.0*Ri)*(1.0-5.0*Ri)
            else
             Fm=(1.0/20.0/Ri)*(1.0/20.0/Ri)
            end if
          end if
         
          if (depth_fog.eq.2.0) depth_fog=30.0        !modify: 20130306 use 30m as shallow fog depth instead

          lamda_1=1./(0.4*(depth_fog+0.02)) + 1./40.0
          lamda=1.0/lamda_1
          Km=lamda*lamda*abs(dwdz)*Fm


          FBL=Km/(2.0*sqrt(0.062*c_adv*depth_fog))

          Kc=1.38*sqrt(0.062*c_adv)*                  !not used, just for verification
     +        depth_fog**1.5

         
c  at last compute LWC near the surface (at z = 1.0 m)
c         if(depth_fog.le.2.0) then
c           z=1.0
c         else
c           z=0.2*depth_fog
c           z=10.0
c         end if
 
         z=10.0

         aa=z/depth_fog
         ax=sqrt(1.0-aa)
         bx=z/FBL
         cx=2.0/(1.0+exp(bx))
         lwc=lwc_max*(ax - cx) 

         if (igrid.eq.NFOG1
     +    .or.igrid.eq.NFOG2
     +    .or.igrid.eq.NFOG3
     +    .or.igrid.eq.NFOG4
     +    .or.igrid.eq.NFOG5
     +    .or.igrid.eq.NFOG6
     +    .or.igrid.eq.NFOG7
     +    .or.igrid.eq.NFOG8
     +    .or.igrid.eq.NFOG9
     +    .or.igrid.eq.NFOG10) then

           write(*,16) ' c_adv=',c_adv,' Kc=',kc,' Km=',km,'LWC=',lwc 
16        format(a28,d11.4,a7,f8.4,a7,f8.4,a7,f8.4) 
         end if

         if (lwc.lt.1.0E-4) lwc = 0.0

       end

