      subroutine ndslfv_monoadvh (grid_gr,
     &                    global_lats_a,lonsperlat,deltim)
!
! a routine to do non-iteration semi-Lagrangain advection
! considering advection  with monotonicity in interpolation
! contact: hann-ming henry juang
! program log
! 2011 02 20 : henry juang, created for ndsl advection
! 2013 06 20 : Henry Juang correct wind direction for north-south advection
!
      use gfs_dyn_machine , only : kind_grid
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_vert_def
      use gfs_dyn_coordinate_def
      use gfs_dyn_physcons
      use gfs_dyn_mpi_def
      implicit none

      real(kind=kind_grid) grid_gr(lonf*lats_node_a_max,lotgr)
      real(kind=kind_grid) plev(lonfull,levs+1)
      integer,intent(in):: global_lats_a(latg)
      integer,intent(in):: lonsperlat(latg)
      real,   intent(in):: deltim
   
      real	uulon(lonfull,levs,latpart)
      real	vvlon(lonfull,levs,latpart)
      real	qqlon(lonfull,levs*ndslhvar,latpart)
      real	rrlon(lonfull,levs*ndslhvar,latpart)

      real	vvlat(latfull,levs,lonpart)
      real	qqlat(latfull,levs*ndslhvar,lonpart)
      real	rrlat(latfull,levs*ndslhvar,lonpart)
      real      rma, rm2a, rdt2, rkt, kappa, pi

      logical 	lprint

      integer mono,mass
      integer nlevs,nvars
      integer ilan,i,j,n,k,kk,lon,lan,lat,lons_lat,jlonf,irc
      integer kp2, kdp, kqq, ktt, kuu, kvv
      integer k2 , kp , kq , kt , ku , kv
      integer k2g, kpg, kqg, ktg, kug, kvg
!
      lprint = .false.

      if( lprint ) print *,' enter ndslfv_monoadvh '
!
      mono  = 1
      mass  = 0
!
      kuu = 1
      kvv = kuu + levs
      ktt = kvv + levs
      kdp = ktt + levs
      kp2 = kdp + levs
      kqq = kp2 + levs

      nvars = ndslhvar
      nlevs = nvars * levs
!
      rdt2 = 0.5 / deltim      
!
! =================================================================
!   prepare wind and variable in flux form with gaussina weight
! =================================================================
!
!$omp parallel do schedule(dynamic,1) private(lan)
!$omp+private(lat,lons_lat,jlonf,rma,rm2a,plev,i,k,ilan)
!$omp+private(kug,kvg,ktg,kpg,k2g,kqg)
!$omp+private(ku ,kv ,kt ,kp ,k2 ,kq )

      do lan=1,lats_node_a

        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf
        rma  = 1. / cosglat(lat) / con_rerth
        rm2a = rma / cosglat(lat)
!
! wind at time step n
        do k=1,levs
          kug=g_uu+k-1
          kvg=g_vv+k-1
          do i=1,lons_lat
            ilan=i+jlonf
            uulon(i,k,lan) = grid_gr(ilan,kug) * rm2a
            vvlon(i,k,lan) = grid_gr(ilan,kvg) * rma
          enddo
        enddo
        if( lprint ) then
          call mymaxmin(uulon(1,1,lan),lons_lat,lonfull,1,' uu1 in deg')
          call mymaxmin(uulon(1,5,lan),lons_lat,lonfull,1,' uu5 in deg')
          call mymaxmin(vvlon(1,1,lan),lons_lat,lonfull,1,' vv in deg')
        endif
!
! u v h dp p2 at n-1
        plev(:,levs+1) = 0.0
        do k=levs,1,-1
          kpg=g_dpm+k-1
          do i=1,lons_lat
            ilan=i+jlonf
            plev(i,k)=plev(i,k+1)+grid_gr(ilan,kpg)
          enddo
        enddo
!
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kp=kdp+k-1
          k2=kp2+k-1
          kug=g_uum+k-1
          kvg=g_vvm+k-1
          ktg=g_ttm+k-1
          kpg=g_dpm+k-1
          do i=1,lons_lat
            ilan=i+jlonf
            qqlon(i,ku,lan) = grid_gr(ilan,kug) 
            qqlon(i,kv,lan) = grid_gr(ilan,kvg) 
            qqlon(i,kt,lan) = grid_gr(ilan,ktg) 
            qqlon(i,kp,lan) = grid_gr(ilan,kpg)
            qqlon(i,k2,lan) = plev(i,k)+plev(i,k+1)
          enddo
        enddo
! rq at n-1
        do k=1,levh
          kq=kqq+k-1
          kqg=g_rm+k-1
          do i=1,lons_lat
            ilan=i+jlonf
            qqlon(i,kq,lan) = grid_gr(ilan,kqg)
          enddo
        enddo
! add surface pressure perturbation for removing resonance
!       rkt = con_g / ( con_rd * 300.0 )
!       do i=1,lons_lat
!         ilan=i+jlonf
!         qqlon(i,kp2,lan) = alog(plev(i,1))+grid_gr(ilan,g_gz)*rkt
!       enddo
!
! save qqlon into n+1 for later as tendency
! u v h dp at n+1
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kp=kdp+k-1
          k2=kp2+k-1
          kug=g_u+k-1
          kvg=g_v+k-1
          ktg=g_t+k-1
          kpg=g_dpn+k-1
          k2g=g_p  +k-1
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,kug) = qqlon(i,ku,lan)
            grid_gr(ilan,kvg) = qqlon(i,kv,lan)
            grid_gr(ilan,ktg) = qqlon(i,kt,lan)
            grid_gr(ilan,kpg) = qqlon(i,kp,lan)
            grid_gr(ilan,k2g) = qqlon(i,k2,lan)
          enddo
        enddo
! change h to theta
!       kappa = con_rd / con_cp
!       do k=1,levs
!         kt=ktt+k-1
!         k2=kp2+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,k2,lan)*0.005)**kappa
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)/pi
!         enddo
!       enddo
! rq  no need for tendency
!
        if( lprint ) then
        print *,' ------------------------------------------- '
        ilan=1+jlonf
        call mymaxmin(grid_gr(ilan,g_u),lons_lat,lonfull,1,' n-1 u ')
        call mymaxmin(grid_gr(ilan,g_v),lons_lat,lonfull,1,' n-1 v ')
        call mymaxmin(grid_gr(ilan,g_t),lons_lat,lonfull,1,' n-1 t ')
        call mymaxmin(grid_gr(ilan,g_rt),lons_lat,lonfull,1,' n-1 q ')
        call mymaxmin(grid_gr(ilan,g_dpn),lons_lat,lonfull,1,' n-1 dp ')
        call mymaxmin(grid_gr(ilan,g_p ),lons_lat,lonfull,1,' n-1 p2 ')
        print *,' ------------------------------------------- '
        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' red u ')
        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' red v ')
        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' red t ')
        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' red q ')
        call mymaxmin(qqlon(1,kdp,lan),lons_lat,lonfull,1,' red dp ')
        call mymaxmin(qqlon(1,kp2,lan),lons_lat,lonfull,1,' red p2 ')
        endif

        call cyclic_cell_intpx(levs,lons_lat,lonf,uulon(1,1,lan))
        call cyclic_cell_intpx(levs,lons_lat,lonf,vvlon(1,1,lan))
        call cyclic_cell_intpx(nlevs,lons_lat,lonf,qqlon(1,1,lan))

        if( lprint ) then
        print *,' ------------------------------------------- '
        call mymaxmin(qqlon(1,kuu,lan),lonfull,lonfull,1,' full u ')
        call mymaxmin(qqlon(1,kvv,lan),lonfull,lonfull,1,' full v ')
        call mymaxmin(qqlon(1,ktt,lan),lonfull,lonfull,1,' full t ')
        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' full q ')
        call mymaxmin(qqlon(1,kdp,lan),lonfull,lonfull,1,' full dp ')
        call mymaxmin(qqlon(1,kp2,lan),lonfull,lonfull,1,' full p2 ')
        endif

        rrlon(:,:,lan) = qqlon(:,:,lan)
!
! first set positive advection in east-west direction 
!
        call cyclic_cell_massadvx(lonfull,levs,nvars,deltim,               &
     &                   uulon(1,1,lan),rrlon(1,1,lan),mass)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),rrlon(1,1,lan),mono)

        if( lprint ) then
        print *,' done cyclic_massadvx with mass= ',mass
        print *,' ------------------------------------------- '
        call mymaxmin(rrlon(1,kuu,lan),lonfull,lonfull,1,' advx u ')
        call mymaxmin(rrlon(1,kvv,lan),lonfull,lonfull,1,' advx v ')
        call mymaxmin(rrlon(1,ktt,lan),lonfull,lonfull,1,' advx t ')
        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' advx q ')
        call mymaxmin(rrlon(1,kdp,lan),lonfull,lonfull,1,' advx dp ')
        call mymaxmin(rrlon(1,kp2,lan),lonfull,lonfull,1,' advx p2 ')
        print *,' done the first x adv for lan=',lan
        print *,' =========================================== '
        endif
 
      enddo

! ---------------------------------------------------------------------
! mpi para from east-west full grid to north-south full grid
! ---------------------------------------------------------------------
!
! para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat

       if( lprint ) print *,' ndslfv_advect transport from we to ns '

       call para_we2ns(vvlon,vvlat,levs,global_lats_a,latg)
       call para_we2ns(rrlon,rrlat,nlevs,global_lats_a,latg)
       call para_we2ns(qqlon,qqlat,nlevs,global_lats_a,latg)

       if( lprint ) then
       print *,' ------------ after we2ns ---------------------- '
       do lon=1,mylonlen
        print *,'  lon=',lon
        call mymaxmin(vvlat(1,1  ,lon),latfull,latfull,1,' we2ns v')
        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' we2ns r')
        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' we2ns q')
       enddo
       endif
!
! ---------------------------------------------------------------------
! -------------- in north-soutn great circle -------------------
! ---------------------------------------------------------------------

       if( lprint ) then
       print *,' ndslfv_advect adv loop in y '
       print *,' mylonlen=',mylonlen
       endif

!$omp parallel do schedule(dynamic,1) private(lon,k,j,ku,kv)
       do lon=1,mylonlen
! 
        if( lprint ) print *,' lon=',lon
! convert wind before advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,latfull/2
            rrlat(j,ku,lon) = -rrlat(j,ku,lon)
            rrlat(j,kv,lon) = -rrlat(j,kv,lon)
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo

! first set advection in north-south direction in great circle through two poles
!
        call cyclic_cell_massadvy(latfull,levs,nvars,deltim,
     &                   vvlat(1,1,lon),rrlat(1,1,lon),mass)
!       call cyclic_mono_advecty (latfull,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),rrlat(1,1,lon),mono)

        if( lprint ) then
        print *,' ------------------------------------------- '
        call mymaxmin(rrlat(1,kuu,lon),latfull,latfull,1,' advry u ')
        call mymaxmin(rrlat(1,kvv,lon),latfull,latfull,1,' advry v ')
        call mymaxmin(rrlat(1,ktt,lon),latfull,latfull,1,' advry t ')
        call mymaxmin(rrlat(1,kqq,lon),latfull,latfull,1,' advry q ')
        call mymaxmin(rrlat(1,kdp,lon),latfull,latfull,1,' advry dp ')
        call mymaxmin(rrlat(1,kp2,lon),latfull,latfull,1,' advry p2 ')
        endif
!
! second set advection in north-south direction in great circle through two poles
!
        call cyclic_cell_massadvy(latfull,levs,nvars,deltim,
     &                   vvlat(1,1,lon),qqlat(1,1,lon),mass)
!       call cyclic_mono_advecty (latfull,levs,nvars,deltim,
!    &                   vvlat(1,1,lon),qqlat(1,1,lon),mono)

        if( lprint ) then
        print *,' ------------------------------------------- '
        call mymaxmin(qqlat(1,kuu,lon),latfull,latfull,1,' advqy u ')
        call mymaxmin(qqlat(1,kvv,lon),latfull,latfull,1,' advqy v ')
        call mymaxmin(qqlat(1,ktt,lon),latfull,latfull,1,' advqy t ')
        call mymaxmin(qqlat(1,kqq,lon),latfull,latfull,1,' advqy q ')
        call mymaxmin(qqlat(1,kdp,lon),latfull,latfull,1,' advqy dp ')
        call mymaxmin(qqlat(1,kp2,lon),latfull,latfull,1,' advqy p2 ')
        print *,' done with y at lon=',lon
        endif

! convert wind back after advy
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          do j=1,latfull/2
            rrlat(j,ku,lon) = -rrlat(j,ku,lon)
            rrlat(j,kv,lon) = -rrlat(j,kv,lon)
            qqlat(j,ku,lon) = -qqlat(j,ku,lon)
            qqlat(j,kv,lon) = -qqlat(j,kv,lon)
          enddo
        enddo

       enddo
!
! ----------------------------------------------------------------------
! mpi para from north-south direction to east-west direeectory 
! ----------------------------------------------------------------------
!
! para qqlat and rrlat to qqlon and rrlon

       if( lprint ) print *,' ndslfv_advect transport from ns to we '

       call para_ns2we(rrlat,rrlon,nlevs,global_lats_a,latg)
       call para_ns2we(qqlat,qqlon,nlevs,global_lats_a,latg)

       if( lprint ) then
       print *,' ------------ after ns2we ---------------------- '
       do lan=1,lats_node_a
        print *,'  lan=',lan
        call mymaxmin(rrlon(1,kqq,lan),lonfull,lonfull,1,' ns2we r')
        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' ns2we q')
       enddo
       endif

! ---------------------------------------------------------------
! ---------------- back to east-west direction ------------------
! ---------------------------------------------------------------
!      print *,' ndslfv_advect adv loop in x for last '

!$omp parallel do schedule(dynamic,1) private(lan)
!$omp+private(lat,lons_lat,jlonf,i,k,ilan)
!$omp+private(kug,kvg,ktg,kpg,k2g,kqg)
!$omp+private(ku ,kv ,kt ,kp ,k2 ,kq )

      do lan=1,lats_node_a

        lat = global_lats_a(ipt_lats_node_a-1+lan)
        lons_lat = lonsperlat(lat)
        jlonf = (lan-1)*lonf

!
! second set advection in x for the second of the pair
!
        call cyclic_cell_massadvx(lonfull,levs,nvars,deltim,               &
     &                   uulon(1,1,lan),qqlon(1,1,lan),mass)
!       call cyclic_mono_advectx (lonfull,levs,nvars,deltim,               &
!    &                   uulon(1,1,lan),qqlon(1,1,lan),mono)

        if( lprint ) then
        print *,' ------------------------------------------- '
        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' adv x q ')
        endif

        do k=1,nlevs
          do i=1,lonfull
            qqlon(i,k,lan) = 0.5 * ( qqlon(i,k,lan) + rrlon(i,k,lan) )
          enddo
        enddo

        if( lprint ) then
        print *,' ------------------------------------------- '
        call mymaxmin(qqlon(1,kuu,lan),lonfull,lonfull,1,' do full u ')
        call mymaxmin(qqlon(1,kvv,lan),lonfull,lonfull,1,' do full v ')
        call mymaxmin(qqlon(1,ktt,lan),lonfull,lonfull,1,' do full t ')
        call mymaxmin(qqlon(1,kqq,lan),lonfull,lonfull,1,' do full q ')
        call mymaxmin(qqlon(1,kdp,lan),lonfull,lonfull,1,' do full dp')
        call mymaxmin(qqlon(1,kp2,lan),lonfull,lonfull,1,' do full p2')
        endif
!
! mass conserving interpolation from full grid to reduced grid
!
        call cyclic_cell_intpx(nlevs,lonf,lons_lat,qqlon(1,1,lan))

        if( lprint ) then
        print *,' ------------------------------------------- '
        call mymaxmin(qqlon(1,kuu,lan),lons_lat,lonfull,1,' do redu u ')
        call mymaxmin(qqlon(1,kvv,lan),lons_lat,lonfull,1,' do redu v ')
        call mymaxmin(qqlon(1,ktt,lan),lons_lat,lonfull,1,' do redu t ')
        call mymaxmin(qqlon(1,kqq,lan),lons_lat,lonfull,1,' do redu q ')
        call mymaxmin(qqlon(1,kdp,lan),lons_lat,lonfull,1,' do redu dp')
        call mymaxmin(qqlon(1,kp2,lan),lons_lat,lonfull,1,' do redu p2')
        print *,' finish horizonatal advection at lan=',lan
        endif

! change theta to h
!       do k=1,levs
!         kt=ktt+k-1
!         k2=kp2+k-1
!         do i=1,lons_lat
!           pi = (qqlon(i,k2,lan)*0.005)**kappa
!           qqlon(i,kt,lan) = qqlon(i,kt,lan)*pi
!         enddo
!       enddo
! u v h dp tendency at n
        do k=1,levs
          ku=kuu+k-1
          kv=kvv+k-1
          kt=ktt+k-1
          kp=kdp+k-1
          k2=kp2+k-1
          kug=g_u+k-1
          kvg=g_v+k-1
          ktg=g_t+k-1
          kpg=g_dpn+k-1
          k2g=g_p  +k-1
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,kug) = (qqlon(i,ku,lan)-grid_gr(ilan,kug))*rdt2
            grid_gr(ilan,kvg) = (qqlon(i,kv,lan)-grid_gr(ilan,kvg))*rdt2
            grid_gr(ilan,ktg) = (qqlon(i,kt,lan)-grid_gr(ilan,ktg))*rdt2
            grid_gr(ilan,kpg) = (qqlon(i,kp,lan)-grid_gr(ilan,kpg))*rdt2
            grid_gr(ilan,k2g) = (qqlon(i,k2,lan)-grid_gr(ilan,k2g))*rdt2
          enddo
        enddo
! rq update
        do k=1,levh
          kq=kqq+k-1
          kqg=g_rt+k-1
          do i=1,lons_lat
            ilan=i+jlonf
            grid_gr(ilan,kqg) = qqlon(i,kq,lan)
          enddo
        enddo

        if( lprint ) then
        print *,' ------------------------------------------- '
        ilan=1+jlonf
        call mymaxmin(grid_gr(ilan,g_u),lons_lat,lonfull,1,' tend u ')
        call mymaxmin(grid_gr(ilan,g_v),lons_lat,lonfull,1,' tend v ')
        call mymaxmin(grid_gr(ilan,g_t),lons_lat,lonfull,1,' tend t ')
        call mymaxmin(grid_gr(ilan,g_rt),lons_lat,lonfull,1,' tend q ')
        call mymaxmin(grid_gr(ilan,g_dpn),lons_lat,lonfull,1,'tend dp')
        do k=1,levs,10
        k2g=g_p+k-1
        print *,' k=',k
        call mymaxmin(grid_gr(ilan,k2g),lons_lat,lonfull,1,' tend p2 ')
        enddo
        print *,' finish horizonatal advection at lan=',lan
        endif
!
      enddo

! 
! ===============================
!
      return
      end