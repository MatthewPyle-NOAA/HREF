      subroutine gloopb(grid_fld,     g3d_fld,       sfc_fld,
     &                  flx_fld,      aoi_fld,       nst_fld,
     &                  lats_nodes_r, global_lats_r, lonsperlar,
     &                  tstep,        phour,         sfalb,   xlon,
     &                  swh,          swhc,          hlw,     hlwc,
!    &                  nbdsw,        nbdlw,         htrsw    htrlwb,
     &                  hprime,       slag,          sdec,    cdec,
     &                  ozplin,       jindx1,        jindx2,  ddy,
     &                  phy_f3d,      phy_f2d,       xlat,    nblck,
     &                  kdt,          restart_step)
!!
!! Code Revision:
!! Sep    2009       Shrinivas Moorthi added nst_fld
!! Oct 11 2009       Sarah Lu, grid_gr replaced by gri_fld
!! Oct 16 2009       Sarah Lu, grid_fld%tracers used
!! Nov 18 2009       Sarah Lu, rain/rainc added to gbphys call arg
!! Dec 14 2009       Sarah Lu, add g3d_fld to calling argument,
!!                             update dqdt after gbphys returns dqdt_v
!! July   2010       Shrinivas Moorthi - Updated for new physics
!! Aug    2010       Shrinivas Moorthi - Recoded 3d diagnostic arrays so that
!                              trap will not occur on call to gbphys
!! Oct 18 2010       Shrinivas Moorthi - Added fscav
!! Dec 23 2010       Sarah Lu, add lgocart to gbphys call arg
!! Apr 06 2012       Henry Juang, add idea related changes
!! Oct 18 2012       Shrinivas Moorthi - Added random number realted chages
!! Jul 26 2012       Jun Wang     pass mpi info to idea_phys
!! Nov 30 2012       Jun Wang     update idea_phys using whole band radiation
!! Dec 27 2012       Jun Wang     move co2 init to gloopb
!! Mar 29 2013       Shrinivas Moorthi - Added dtphys option
!! Apr 08 2013       Jun Wang     add idea init variables to restart file
!! Oct 21 2013       Henry Juang  compute prsi from model top
!! Oct 31 2013       Xingren Wu   add flx_fld%dusfci/dvsfci
!! Mar 26 2014       Xingren Wu   add aoi_fld
!! May 05 2014       Jun Wang     add cgwf,prslrd0 in gbphys argument list
!! jun    2014       y-t hou,  revised sw sfc spectral component fluxes
!!                             and ocean albedoes (no ice contamination) for coupled mdl 
!! Jun    2014       Xingren Wu   update net SW fluxes over the ocean
!!                                (no ice contamination)
!! Jul    2014       Xingren Wu   Add Sea/Land/Ice Mask - slimsk_cpl
!! Jul -- 2014       Shrinivas Moorthi updating to semi-lagrangian and latest gfs physics 
!!                                     and some cleanup - bugfix for random number
!! Aug 01 2014       Shrinivas Moorthi - add tracer fixer option
!! Sep 16 2014       Shrinivas Moorthi - cleanup and rearrange argumets
!! Sep 30 2014       Sarah Lu, remove fscav (the option to compute tracer
!!                             scavenging in GFS is disable)
!! Dec 21 2014       Jun Wang, add g3d_fld cnv_qc to gbphys argument,
!!                             update 3D fields after gbphys
!!
      use resol_def
      use layout1
      use gg_def
      use vert_def
      use date_def
      use namelist_physics_def
      use coordinate_def
      use module_ras , only : ras_init
      use physcons, grav => con_g , rerth => con_rerth, rk => con_rocp
      use ozne_def
      use d3d_def
      use gfs_physics_sfc_flx_mod
      use gfs_physics_nst_var_mod
      use gfs_physics_aoi_var_mod
      use gfs_physics_gridgr_mod, ONLY: Grid_Var_Data
      use gfs_physics_g3d_mod,    ONLY: G3D_Var_Data
      use mpi_def,                only: mpi_r_io_r, mpi_comm_all
      use mersenne_twister
      use idea_composition, only: prlog,pr_idea,amgm,amgms,nlev_co2,k43,
     &                            nlevc_h2o,k71,gg,prsilvl
      use efield, only: efield_init

      implicit none
      include 'mpif.h'
!
      TYPE(Grid_Var_Data)       :: grid_fld
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
      TYPE(Nst_Var_Data)        :: nst_fld
      TYPE(G3D_Var_Data)        :: g3d_fld
      TYPE(AOI_Var_Data)        :: aoi_fld

!
      integer nblck, kdt, nbdsw, nbdlw
      logical restart_step
!!
      real(kind=kind_phys)  tstep, phour,slag, sdec, cdec
      real(kind=kind_phys)  plyr(levs)
      real(kind=kind_phys), dimension(ngptc)       :: dpshc, pgr
      real(kind=kind_phys), dimension(ngptc,levs)  :: prsl, prslk
     &,                           phil, gu, gv, gt, adu, adv, adt
      real(kind=kind_phys), dimension(ngptc,levs+1) :: prsi, prsik
     &,                           phii
!!
      real (kind=kind_rad), dimension(ngptc,levs,ntrac) :: gr, adr
!!
      real (kind=kind_rad), dimension(lonr,lats_node_r) :: xlon, xlat
     &,                           sfalb
!!
      real (kind=kind_rad), dimension(ngptc,levs,nblck,lats_node_r) ::
     &                            swh, swhc, hlw, hlwc
!!
      real (kind=kind_rad)  hprime(nmtvr,lonr,lats_node_r)

!idea add by hmhj
      real (kind=kind_rad) hlwd(ngptc,levs,6)
!    &,                    htrswb(ngptc,levs,nbdsw,nblck,lats_node_r)
!    &,                    htrlwb(ngptc,levs,nbdlw,nblck,lats_node_r)
!!
      real  (kind=kind_phys)
     &     phy_f3d(ngptc,levs,num_p3d+npdf3d,nblck,lats_node_r),
     &     phy_f2d(lonr,lats_node_r,num_p2d)
!!
      real (kind=kind_phys) dtp,dtf,solhr,clstp
!!
      integer              lats_nodes_r(nodes)
      integer, dimension(latr) :: global_lats_r, lonsperlar
!
      integer              i,j,k,kk,n,l,lan,lat,ii,lonrbm,jj
     &,                    lon_dim,lons_lat,nsphys
!
!timers______________________________________________________---
 
!     real*8 rtc ,timer1,timer2
!     real(kind=kind_evod) global_times_b(latr,nodes)
 
!timers______________________________________________________---
!
      logical, parameter :: flipv = .true.
      real(kind=kind_phys), parameter :: pt01=0.01, pt00001=1.0e-5
     &,                                  thousnd=1000.0
     &,                                  cons_0=0.0,   cons_24=24.0
     &,                                  cons_99=99.0, cons_1p0d9=1.0E9
     &,                                  qmin=1.0e-10
!
! for nrl/nasa ozone production and distruction rates:(input through fixio)
! ---------------------------------------------------
      integer              jindx1(lats_node_r),jindx2(lats_node_r) !for ozone interpolaton
      real(kind=kind_phys) ozplin(latsozp,levozp,pl_coeff,timeoz)
     &,                    ddy(lats_node_r)                        !for ozone interpolaton
     &,                    ozplout(levozp,lats_node_r,pl_coeff)
!!
      real(kind=kind_phys), allocatable :: acv(:,:),acvb(:,:),acvt(:,:)
      save acv,acvb,acvt
!!
!     integer, parameter :: maxran=5000
!     integer, parameter :: maxran=3000
!     integer, parameter :: maxran=6000, maxsub=6, maxrs=maxran/maxsub
      integer, parameter :: maxran=3000, maxsub=6, maxrs=maxran/maxsub
      type (random_stat) :: stat(maxrs)
      real (kind=kind_phys), allocatable, save :: rannum_tank(:,:,:)
      real (kind=kind_phys), allocatable       :: rannum(:), wrkn(:)
      integer, allocatable                     :: indxr(:)
      integer iseed, nrc, seed0, kss, ksr, iseedl, latseed
     &,       nf0,nf1,ind,nt,indod,indev
      real(kind=kind_evod) fd2, wrk(1)

      logical first
      data first/.true./
      save    first, seed0
!
      integer nlons_v(ngptc)

      real(kind=kind_phys), dimension(ngptc) :: sinlat_v,coslat_v,rqtk
      real(kind=kind_phys), dimension(ngptc,lsoil) :: smc_v,stc_v,slc_v
      real(kind=kind_phys), dimension(ngptc,levs)  :: vvel,dtdt,dqdt_v,
     &                                                cnvqc_v

      real(kind=kind_phys)  hprime_v(ngptc,nmtvr)
!    &,                     phy_f3dv(ngptc,LEVS,num_p3d+npdf3d)
     &,                     phy_f2dv(ngptc,num_p2d)
     &,                     rannum_v(ngptc,nrcm)
     &,                     ozplout_v(ngptc,levozp,pl_coeff)

      real(kind=kind_phys) trcg(latr,ntrac,2),trcj(lats_node_r,ntrac,2)
     &,                    trcp(lonr,ntrac,2)

      real(kind=kind_phys), allocatable, save :: sumtrc(:,:), adjtrc(:)

!   variables for stochastic physics
!-----------------------------------
!     real (kind=kind_phys),dimension(ngptc,levs)        ::
!    &                      uphys,vphys,tphys,qphys,sppt_wts,shum_wts,
!    &                      skebu_wts,skebv_wts,dtdt
!     real (kind=kind_phys),allocatable,dimension(:,:,:) ::
!    &                      sppt_wt,shum_wt,skebu_wt,skebv_wt
!     real (kind=kind_phys),allocatable,dimension(:,:)   ::
!    &                      trigger_wt

      real(kind=kind_phys) dt3dt_v(ngptc,levs,6), du3dt_v(ngptc,levs,4)
     &,                    dv3dt_v(ngptc,levs,4)
     &,                    dq3dt_v(ngptc,levs,5+pl_coeff)
      real(kind=kind_phys) upd_mfv(ngptc,levs), dwn_mfv(ngptc,levs)
     &,                    det_mfv(ngptc,levs)
!    &,                    det_mfv(ngptc,levs), dkh_v(ngptc,levs)
!    &,                    rnp_v(ngptc,levs)
      real(kind=kind_phys), dimension(ngptc) :: nirbmdi, nirdfdi,
     &                        visbmdi, visdfdi, nirbmui, nirdfui,
     &                        visbmui, visdfui
     &,                       aoi_du,     aoi_dv,     aoi_dt, aoi_dq
     &,                       aoi_dlw,    aoi_dsw,    aoi_dnirbm
     &,                       aoi_dnirdf, aoi_dvisbm, aoi_dvisdf
     &,                       aoi_rain,   aoi_nlw,    aoi_nsw
     &,                       aoi_nnirbm, aoi_nnirdf, aoi_nvisbm
     &,                       aoi_nvisdf

     &,                       aoi_dusfci,  aoi_dvsfci,  aoi_dtsfci
     &,                       aoi_dqsfci,  aoi_dlwsfci, aoi_dswsfci
     &,                       aoi_dnirbmi, aoi_dnirdfi, aoi_dvisbmi
     &,                       aoi_dvisdfi, aoi_nlwsfci, aoi_nswsfci
     &,                       aoi_nnirbmi, aoi_nnirdfi, aoi_nvisbmi
     &,                       aoi_nvisdfi, aoi_t2mi,    aoi_q2mi
     &,                       aoi_u10mi,   aoi_v10mi,   aoi_tseai
     &,                       aoi_psurfi

     &,                       nst_xt,  nst_xs,   nst_xu,  nst_xv, nst_xz
     &,                       nst_zm,  nst_xtts, nst_xzts, nst_d_conv
     &,                       nst_ifd, nst_dt_cool, nst_qrain

     &,                       nst_tref, nst_z_c, nst_c_0, nst_c_d
     &,                       nst_w_0,  nst_w_d

      real(kind=kind_phys) work1, tem
      integer              njeff,lon,iblk,item, nn, nnr, nnrcm, dbgu

! idea local vars
      real(kind=kind_phys) pmod(LEVS),gg1(levs),philco2(levs),
     &                     qtrac(levs,ntrac),prsilvl1(levs+1)
      real(kind=kind_phys),allocatable :: prpa(:)
      integer info,pelat1,pelatall,lanlat1
!
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!
      if (first) then
        allocate (acv(lonr,lats_node_r))
        allocate (acvb(lonr,lats_node_r))
        allocate (acvt(lonr,lats_node_r))

        allocate (sumtrc(ntrac,3), adjtrc(ntrac))
        adjtrc = 1.0
!
        if (.not. newsas .or. cal_pre) then  ! random number needed for RAS and old SAS
          if (random_clds) then ! create random number tank
!                                 -------------------------
            seed0 = idate(1) + idate(2) + idate(3) + idate(4)

            call random_setseed(seed0)
            call random_number(wrk)
            seed0 = seed0 + nint(wrk(1)*thousnd)
            if (me == 0) print *,' seed0=',seed0,' idate=',idate,
     &                           ' wrk=',wrk
!
            if (.not. allocated(rannum_tank))
     &                allocate (rannum_tank(lonr,maxran,lats_node_r))
            if (.not. allocated(rannum)) allocate (rannum(lonr*maxrs))
            lonrbm = lonr / maxsub
            if (me == 0) write(0,*)' maxran=',maxran,' maxrs=',maxrs,
     &          'maxsub=',maxsub,' lonrbm=',lonrbm,
     &          ' lats_node_r=',lats_node_r
            do j=1,lats_node_r
              iseedl = global_lats_r(ipt_lats_node_r-1+j) + seed0
              call random_setseed(iseedl)
              call random_number(rannum)
!$omp parallel do  shared(j,lonr,lonrbm,rannum,rannum_tank)
!$omp+private(nrc,nn,i,ii,k,kk)
              do nrc=1,maxrs
                nn = (nrc-1)*lonr
                do k=1,maxsub
                  kk = k - 1
                  do i=1,lonr
                    ii = kk*lonrbm + i
                    if (ii > lonr) ii = ii - lonr
                    rannum_tank(i,nrc+kk*maxrs,j) = rannum(ii+nn)
                  enddo
                enddo
              enddo
            enddo
            if (allocated(rannum)) deallocate (rannum)
          endif
        endif
!
        if (me  ==  0) then
!         write(0,*)' seed0=',seed0,' idate=',idate,' wrk=',wrk
          if (num_p3d == 3) write(0,*)' USING Ferrier-MICROPHYSICS'
          if (num_p3d == 4) write(0,*)' USING ZHAO-MICROPHYSICS'
        endif
        if (fhour == 0.0) then
          do j=1,lats_node_r
            do i=1,lonr
              phy_f2d(i,j,num_p2d) = 0.0
            enddo
          enddo
        endif
       
        if (ras) call ras_init(levs, me)
!
!***************************************idea below*****************************************
        if ( lsidea ) then
!
!find PE with lat 1
          pelat1  = -1
          lanlat1 = -1
          findlat1pe: do lan=1,lats_node_r
            lat=global_lats_r(ipt_lats_node_r-1+lan)
            if(lat == 1) then
              pelat1  = me
              lanlat1 = lan
               print *,'pelat1=',pelat1
              exit findlat1pe
            endif
          enddo findlat1pe
          call mpi_reduce(pelat1,pelatall,1,mpi_integer,MPI_MAX,0,
     &                    mpi_comm_all,info)
          if(me == 0) then
            print *,'pelatall=',pelatall
            pelat1 = pelatall
          endif
          call mpi_bcast(pelat1,1,mpi_integer,0,mpi_comm_all,info)
          print *,'pelat1=',pelat1,'lanlat1=',lanlat1,
     &            'gen_coord_hybrid=',gen_coord_hybrid,
     &            'thermodyn_id=',thermodyn_id
!
! set plyr from lat1 pe
          if(me == pelat1) then
            do k=1,levs
              plyr(k) = grid_fld%p(1,lanlat1,k)
            enddo
            print *,' plyr in gloopb ',(plyr(k),k=1,levs)
          endif
          call mpi_bcast(plyr,levs,mpi_r_io_r,pelat1,mpi_comm_all,info)
          call idea_composition_init(levs,plyr)
!
          call idea_solar_init(levs)
          call idea_tracer_init(levs)
!h2ocin
          print *,'in gloopb,nlevc_h2o=',nlevc_h2o,'k71=',k71,
     &            'levs=',levs,levs-k71+1
          allocate(prpa(nlevc_h2o))
          prpa(1:nlevc_h2o) = 100.*pr_idea(k71:levs)
          call h2ocin(prpa,nlevc_h2o,me,mpi_r_io_r,mpi_comm_all)
          deallocate(prpa)
!ion
          call efield_init()
!o3
          call o3ini(levs)
!
!co2
          if(.not. restart_step) then
            if(me == pelat1) then
!compute phil
              do n=1,ntrac
                do k=1,levs
                  qtrac(k,n) = grid_fld%tracers(n)%flds(1,lanlat1,k)
                enddo
              enddo
!             print *,'in gloopb,ps=',grid_fld%ps(1,lanlat1)
              call getphilvl(levs,ntrac, grid_fld%ps(1,lanlat1),
     &                   grid_fld%t(1,lanlat1,1:levs),qtrac,
     &                   grid_fld%dp(1,lanlat1,1:levs),gen_coord_hybrid,
     &                   thermodyn_id,philco2,prsilvl1)

!! change prsi from cb to pascal
              prsilvl1 = prsilvl1*1000.
!             print *,'prsi=',prsi(1,103)*1000.
!!compute gravity
              call gravco2(levs,philco2,sfc_fld%oro(1,lanlat1),gg1)
            endif
            call mpi_bcast(gg1,levs,mpi_r_io_r,pelat1,mpi_comm_all,info)
!
            call mpi_bcast(prsilvl1,levs+1,mpi_r_io_r,pelat1,
     &                     mpi_comm_all,info)
!
          endif
          pmod = pr_idea*100.
          if(.not. allocated(prsilvl)) then
            allocate(prsilvl(levs+1), gg(levs))
            do k=1,levs
              prsilvl(k) = prsilvl1(k)
              gg(k)      = gg1(k)
              amgms(k)   = amgm(k)
            enddo
            prsilvl(levs+1) = prsilvl1(levs+1)
          endif

!         print *,'bf co2cin,prlog=',prlog(k43:k43+3)
          call co2cin(prlog(k43),pmod(k43),amgms(k43),gg(k43),
     &                nlev_co2,me,mpi_r_io_r,mpi_comm_all)
!
          call ideaca_init(prsilvl,levs+1)

        endif
!***************************************idea above****************************************

        first = .false.

      endif
!
      if (semilag) then
        nsphys = max(int(tstep/dtphys),1)
        dtp    = tstep / nsphys
        dtf    = dtp
      else
        nsphys = max(int((tstep+tstep)/dtphys+0.9999),1)
        dtp    = (tstep+tstep)/nsphys
        dtf    = 0.5*dtp
        if(lsfwd) dtf = dtp
      endif
!
      solhr = mod(phour+idate(1),cons_24)

! **************  Ken Campana Stuff  ********************************
!...  set switch for saving convective clouds
      if(lscca.and.lsswr) then
        clstp = 1100+min(fhswr,fhour,cons_99)  !initialize,accumulate,convert
      elseif(lscca) then
        clstp = 0100+min(fhswr,fhour,cons_99)  !accumulate,convert
      elseif(lsswr) then
        clstp = 1100                           !initialize,accumulate
      else
        clstp = 0100                           !accumulate
      endif
! **************  Ken Campana Stuff  ********************************
!
!
      if (.not. newsas .or. cal_pre) then  ! random number needed for RAS and old SAS
                                           ! and when cal_pre=.true.

        if (random_clds) then
          iseed = mod(100.0*sqrt(fhour*3600),cons_1p0d9) + 1 + seed0

!         write(0,*)' After Initialization in gloopb iseed=',iseed

          nnrcm  = nrcm*nsphys
          allocate(wrkn(nnrcm), indxr(nnrcm))
          call random_setseed(iseed)
          call random_number(wrkn)
          do nrc=1,nnrcm
            indxr(nrc) = max(1, min(nint(wrkn(nrc)*maxran)+1,maxran))
          enddo
        endif
        if (allocated(wrkn)) deallocate (wrkn)
      endif
!
! do ozone i/o and latitudinal interpolation to local gaussian lats
!
      if (ntoz > 0) then
       call ozinterpol(me,lats_node_r,lats_node_r,idate,fhour,
     &                 jindx1,jindx2,ozplin,ozplout,ddy)
      endif

!-------------------------------------------------------------------------
!  for stochastic physics
!     if (sppt > tiny(sppt)) then
!        allocate(sppt_wt(lonr,lats_node_r,levs))
!        call get_pattern_sppt(
!    &     ls_node,ls_nodes,max_ls_nodes,
!    &     lats_nodes_r,global_lats_r,lonsperlar,
!    &     plnev_r,plnod_r,
!    &     sppt_wt,nsphys*dtf)
!     endif
!     if (shum > tiny(shum)) then
!        allocate(shum_wt(lonr,lats_node_r,levs))
!        call get_pattern_shum(
!    &     ls_node,ls_nodes,max_ls_nodes,
!    &     lats_nodes_r,global_lats_r,lonsperlar,
!    &     plnev_r,plnod_r,
!    &     shum_wt,nsphys*dtf)
!     endif
!     if (strig > tiny(strig)) then
!        allocate(trigger_wt(lonr,lats_node_r))
!        call get_pattern_strig(
!    &     ls_node,ls_nodes,max_ls_nodes,
!    &     lats_nodes_r,global_lats_r,lonsperlar,
!    &     plnev_r,plnod_r,
!    &     trigger_wt,nsphys*dtf)
!     endif
!     if (skeb > tiny(skeb)) then
!        allocate(skebu_wt(lonr,lats_node_r,levs))
!        allocate(skebv_wt(lonr,lats_node_r,levs))
!        call get_pattern_skeb(trie_ls(1,1,p_ze),
!    &         trio_ls(1,1,p_ze),
!    &         trie_ls(1,1,p_di),
!    &         trio_ls(1,1,p_di),
!    &         ls_node,ls_nodes,max_ls_nodes,
!    &         lats_nodes_r,global_lats_r,lonsperlar,
!    &         epsedn,epsodn,snnp1ev,snnp1od,
!    &         plnev_r,plnod_r,plnew_r,plnow_r,
!    &         skebu_wt,
!    &     
!---------------------------main latitude loop starts---------------------------------
!
      do lan=1,lats_node_r
         lat         = global_lats_r(ipt_lats_node_r-1+lan)
         lon_dim     = lon_dims_r(lan)
         lons_lat    = lonsperlar(lat)
!        pwatp       = 0.
!        trcp(:,:,:) = 0.0
         do n=1,ntrac
           do i=1,lonr
             trcp(i,n,1) = 0.0
             trcp(i,n,2) = 0.0
           enddo
         enddo

!$omp parallel do  schedule(dynamic,1) private(lon)
!$omp+private(hprime_v,stc_v,smc_v,slc_v)
!$omp+private(nlons_v,sinlat_v,coslat_v,ozplout_v,rannum_v)
!$omp+private(prslk,prsl,prsik,prsi,phil,phii,dpshc)
!$omp+private(gu,gv,gt,gr,vvel)
!$omp+private(adt,adr,adu,adv,pgr,rqtk)
!$omp+private(phy_f2dv,dtdt)
!$omp+private(dt3dt_v,du3dt_v,dv3dt_v,dq3dt_v)
!$omp+private(upd_mfv,dwn_mfv,det_mfv)
!$omp+private(dqdt_v,cnvqc_v,hlwd)
!$omp+private(njeff,iblk,i,j,k,n,item,nn,nnr,dbgu)
!$omp+private(nirbmdi, nirdfdi,visbmdi, visdfdi, nirbmui, nirdfui)
!$omp+private(visbmui,visdfui,aoi_du,aoi_dv,aoi_dt,aoi_dq)
!$omp+private(aoi_dlw,aoi_dsw,aoi_dnirbm,aoi_dnirdf,aoi_dvisbm)
!$omp+private(aoi_dvisdf,aoi_rain, aoi_nlw,aoi_nsw,aoi_nnirbm)
!$omp+private(aoi_nnirdf,aoi_nvisbm,aoi_nvisdf,aoi_dusfci)
!$omp+private(aoi_dvsfci,aoi_dtsfci,aoi_dqsfci,aoi_dlwsfci)
!$omp+private(aoi_dswsfci,aoi_dnirbmi,aoi_dnirdfi,aoi_dvisbmi)
!$omp+private(aoi_dvisdfi,aoi_nlwsfci,aoi_nswsfci,aoi_nnirbmi)
!$omp+private(aoi_nnirdfi,aoi_nvisbmi,aoi_nvisdfi,aoi_t2mi,aoi_q2mi)
!$omp+private(aoi_u10mi,aoi_v10mi,aoi_tseai,aoi_psurfi)
!$omp+private(nst_xt,nst_xs,nst_xu,nst_xv,nst_xz,nst_zm,nst_xtts)
!$omp+private(nst_xzts,nst_d_conv,nst_ifd,nst_dt_cool,nst_qrain)
!$omp+private(nst_tref,nst_z_c,nst_c_0,nst_c_d,nst_w_0,nst_w_d)

!!$omp+private(upd_mfv,dwn_mfv,det_mfv,dkh_v,rnp_v)
!!!$omp+private(njeff,iblk,i,j,k,n,item,nn,nnr,dbgu)
!!!$omp+private(njeff,iblk,ilan,i,j,k,n,item)
!!!!$omp+private(temlon,temlat,lprnt,ipt)


!---------------------------main longitude loop starts--------------------------------
        do lon=1,lons_lat,ngptc
!!

          njeff = min(ngptc,lons_lat-lon+1)
          iblk  = (lon-1)/ngptc + 1
!
!         dbgu = 300 + lon
          dbgu = 0

!         write(dbgu,*)' GLOOPB : LON=',lon,' lons_lat=',lons_lat,
!    &' njeff=',njeff,' iblk=',iblk
!    &,' zorl=',sfc_fld%zorl(njeff-3:njeff,lan),' lan=',lan
!        write(dbgu,*)' dbgu=',dbgu,' lon=',lon,' lan=',lan
!!
!     write(dbgu,*)' lan=',lan,' pgr=',pgr(i),' i=',i,' njeff=',njeff
!     print *,' lan=',lan,' pgr=',pgr(i),' grid_gr=',grid_gr(ilan,g_ps)
!    &,' i=',i,' lan=',lan
!

          do k = 1, LEVS
            do i = 1, njeff
              item = lon+i-1
              gu(i,k)     = grid_fld%u(item,lan,k)    ! real u    
              gv(i,k)     = grid_fld%v(item,lan,k)    ! real v
              gt(i,k)     = grid_fld%t(item,lan,k)    ! temperature K
              prsl(i,k)   = grid_fld%p(item,lan,k)    ! pascal  
              vvel(i,k)   = grid_fld%dpdt(item,lan,k) ! pascal/sec
!hmhj         prsi(i,k+1) = prsi(i,k) - grid_fld%dp(item,lan,k)  !pascal
            enddo
          enddo

!     write(0,*)' mintem=',minval(gt(1:njeff,:))
!    &,' maxtem=',maxval(gt(1:njeff,:)),' lan=',lan

!hmhj prsi should be computed from model top for accuracy
          do i=1,njeff
            prsi (i,levs+1) = 0.0
            prsik(i,levs+1) = 0.0
          enddo
          do k=levs,1,-1
            do i = 1, njeff
              item = lon+i-1
              prsi(i,k) = prsi(i,k+1) + grid_fld%dp(item,lan,k)  !pascal
            enddo
          enddo
          do i=1,njeff
            pgr(i) = prsi(i,1)
          enddo
          do n=1,ntrac
            do k=1,levs
              do i=1,njeff
                gr(i,k,n)= grid_fld%tracers(n)%flds(lon+i-1,lan,k) ! kg/kg
              enddo
            enddo
          enddo

          do i=1,njeff
            phil(i,levs) = 0.0 ! will force calculation of geopotential in gbphys.
            dpshc(i)     = 0.3 * prsi(i,1)
!
            nlons_v(i)   = lons_lat
            sinlat_v(i)  = sinlat_r(lat)
            coslat_v(i)  = coslat_r(lat)
          enddo

 !        if (gen_coord_hybrid .and. thermodyn_id == 3) then
            do i=1,njeff
              prslk(i,1) = 0.0 ! forces calculation of (p/p00)**kappa in gbphys
              prsik(i,1) = 0.0 ! forces calculation of (p/p00)**kappa in gbphys
            enddo
 !        else
 !          do k = 1, levs
 !            do i = 1, njeff
 !              prslk(i,k) = (prsl(i,k)*pt00001)**rk
 !              prsik(i,k) = (prsi(i,k)*pt00001)**rk
 !            enddo
 !          enddo
 !        endif

          if (ntoz .gt. 0) then
            do j=1,pl_coeff
              do k=1,levozp
                do i=1,njeff
                  ozplout_v(i,k,j) = ozplout(k,lan,j)
                enddo
              enddo
            enddo
          endif

          do k=1,lsoil
            do i=1,njeff
              item = lon+i-1
              smc_v(i,k) = sfc_fld%smc(k,item,lan)
              stc_v(i,k) = sfc_fld%stc(k,item,lan)
              slc_v(i,k) = sfc_fld%slc(k,item,lan)
            enddo
          enddo
          do k=1,nmtvr
            do i=1,njeff
              hprime_v(i,k) = hprime(k,lon+i-1,lan)
            enddo
          enddo
!!
          do j=1,num_p2d
            do i=1,njeff
              phy_f2dv(i,j) = phy_f2d(lon+i-1,lan,j)
            enddo
          enddo

          if (ldiag3d) then
            do k=1,6
              do j=1,levs
                do i=1,njeff
                  dt3dt_v(i,j,k) = dt3dt(i,j,k,iblk,lan)
                enddo
              enddo
            enddo
            do k=1,4
              do j=1,levs
                do i=1,njeff
                  du3dt_v(i,j,k) = du3dt(i,j,k,iblk,lan)
                  dv3dt_v(i,j,k) = dv3dt(i,j,k,iblk,lan)
                enddo
              enddo
            enddo
            do k=1,5+pl_coeff
              do j=1,levs
                do i=1,njeff
                  dq3dt_v(i,j,k) = dq3dt(i,j,k,iblk,lan)
                enddo
              enddo
            enddo
            do j=1,levs
              do i=1,njeff
                upd_mfv(i,j) = upd_mf(i,j,iblk,lan)
                dwn_mfv(i,j) = dwn_mf(i,j,iblk,lan)
                det_mfv(i,j) = det_mf(i,j,iblk,lan)
              enddo
            enddo
          endif

          if (lgocart) then
            do k=1,levs
              do i=1,njeff
                dqdt_v(i,k) = 0.
                upd_mfv(i,k) = 0.
                dwn_mfv(i,k) = 0.
                det_mfv(i,k) = 0.
                cnvqc_v(i,k) = 0.
              enddo
            enddo
          endif
!

!*************************************idea below**************************
          if( lsidea ) then
!hmhj debug
!           do n=1,ntrac
!             call mymaxmin(gr(1,1,n),njeff,ngptc,levs,' gr b idea ')
!           enddo
!hmhj debug
!           call mymaxmin(gt,njeff,ngptc,levs,' gt before idea_phys ')
            call idea_phys(njeff,ngptc,levs,prsi,prsl,
     &                     gu,gv,gt,gr,ntrac,dtp,lat,
     &                     solhr,slag,sdec,cdec,sinlat_v,coslat_v,
     &                     xlon(lon,lan),xlat(lon,lan),
     &                     sfc_fld%oro(lon,lan),flx_fld%coszen(lon,lan),
     &                     swh(1,1,iblk,lan),hlw(1,1,iblk,lan),hlwd,
     &                     thermodyn_id,sfcpress_id,gen_coord_hybrid,
     &                     me,mpi_r_io_r,MPI_COMM_ALL)
!hmhj debug
!            hlwd(:,:,6) = 0.0
!            do n=1,ntrac
!              print *,' =========== gr at n=',n,' ==============='
!              call mymaxmin(gr(1,1,n),njeff,ngptc,levs,' gr a idea ')
!            enddo
!            do n=1,6
!              print *,' =========== hlwd at n=',n,' ==============='
!              call mymaxmin(hlwd(1,1,n),njeff,ngptc,levs,' hlwd a idea ')
!            enddo
!            call mymaxmin(gu,njeff,ngptc,levs,' gu after idea_phys ')
!            call mymaxmin(gv,njeff,ngptc,levs,' gv after idea_phys ')
!            call mymaxmin(gt,njeff,ngptc,levs,' gt after idea_phys ')
          endif
!*************************************idea above**************************
!
!     write(dbgu,*)' before gbphys:', njeff,ngptc,levs,lsoil,lsm,       &
!    &      ntrac,ncld,ntoz,ntcw,                                       &
!    &      nmtvr,nrcm,levozp,lonr,latr,jcap,num_p3d,num_p2d,           &
!    &      kdt,lat,me,pl_coeff,ncw,flgmin,crtrh,cdmbgwd
!    &,'sizer1=',size(phy_f3dv,1),'sizer2=',size(phy_f3dv,2)
!    &,'sizer3=',size(phy_f3dv,3)
!    &,' ccwf=',ccwf,' dlqf=',dlqf,' lsidea=',lsidea
!    &,' evpco=',evpco,' wminco=',wminco
!     write(dbgu,*)' tisfc=',sfc_fld%tisfc(1:20,lan),' lan=',lan,' lon='&
!    &,           lon
!
!hmhj debug
!         call mymaxmin(gt,njeff,ngptc,levs,' gt before gbphys ')
!
          do k=1,levs
            do i=1,njeff
              dtdt(i,k) = 0.0
            enddo
          enddo

!     For tracer fixer

          do n=1,ntrac
            if (fixtrc(n)) then
              do k=1,levs
                do i=1,njeff
                  trcp(lon+i-1,n,1) = trcp(lon+i-1,n,1)
     &                         + gr(i,k,n) * (prsi(i,k) - prsi(i,k+1))
                enddo
              enddo
            endif
          enddo
          do i=1,njeff
            rqtk(i) = 0.0
          enddo

          if (cplflx) then
            do i=1,njeff
              item = lon+i-1
              nirbmdi(i)    = aoi_fld%nirbmdi(item,lan) 
              nirdfdi(i)    = aoi_fld%nirdfdi(item,lan)
              visbmdi(i)    = aoi_fld%visbmdi(item,lan)
              visdfdi(i)    = aoi_fld%visdfdi(item,lan)
              nirbmui(i)    = aoi_fld%nirbmui(item,lan)
              nirdfui(i)    = aoi_fld%nirdfui(item,lan)
              visbmui(i)    = aoi_fld%visbmui(item,lan)
              visdfui(i)    = aoi_fld%visdfui(item,lan)
!
              aoi_du(i)     = aoi_fld%dusfc(item,lan)
              aoi_dv(i)     = aoi_fld%dvsfc(item,lan) 
              aoi_dt(i)     = aoi_fld%dtsfc(item,lan) 
              aoi_dq(i)     = aoi_fld%dqsfc(item,lan) 
              aoi_dlw(i)    = aoi_fld%dlwsfc(item,lan) 
              aoi_dsw(i)    = aoi_fld%dswsfc(item,lan) 
              aoi_dnirbm(i) = aoi_fld%dnirbm(item,lan)
              aoi_dnirdf(i) = aoi_fld%dnirdf(item,lan)
              aoi_dvisbm(i) = aoi_fld%dvisbm(item,lan)
              aoi_dvisdf(i) = aoi_fld%dvisdf(item,lan) 
              aoi_rain(i)   = aoi_fld%rain(item,lan)
              aoi_nlw(i)    = aoi_fld%nlwsfc(item,lan)
              aoi_nsw(i)    = aoi_fld%nswsfc(item,lan)
              aoi_nnirbm(i) = aoi_fld%nnirbm(item,lan)
              aoi_nnirdf(i) = aoi_fld%nnirdf(item,lan)
              aoi_nvisbm(i) = aoi_fld%nvisbm(item,lan)
              aoi_nvisdf(i) = aoi_fld%nvisdf(item,lan)
            enddo
          endif
          if (nst_fcst > 0) then
            do i=1,njeff
              item = lon+i-1
              nst_xt(i)      = nst_fld%xt(item,lan) 
              nst_xs(i)      = nst_fld%xs(item,lan) 
              nst_xu(i)      = nst_fld%xu(item,lan) 
              nst_xv(i)      = nst_fld%xv(item,lan) 
              nst_xz(i)      = nst_fld%xz(item,lan) 
              nst_zm(i)      = nst_fld%zm(item,lan) 
              nst_xtts(i)    = nst_fld%xtts(item,lan) 
              nst_xzts(i)    = nst_fld%xzts(item,lan) 
              nst_d_conv(i)  = nst_fld%d_conv(item,lan) 
              nst_ifd(i)     = nst_fld%ifd(item,lan) 
              nst_dt_cool(i) = nst_fld%dt_cool(item,lan) 
              nst_qrain(i)   = nst_fld%qrain(item,lan) 
            enddo
          endif
!
          do nn=1,nsphys             ! physics sub-steps

            if (.not. newsas .or. cal_pre) then
              if (random_clds) then
                nnr = (nn-1)*nrcm
                do j=1,nrcm
                  do i=1,njeff
                    rannum_v(i,j)=rannum_tank(lon+i-1,indxr(nnr+j),lan)
                  enddo
                enddo
              else
                do j=1,nrcm
                  do i=1,njeff
                    rannum_v(i,j) = 0.6    ! This is useful for debugging
                  enddo
                enddo
              endif
            endif

            call gbphys                                                 &
!  ---  inputs:
     &    ( njeff,ngptc,levs,lsoil,lsm,ntrac,ncld,ntoz,ntcw,            &
     &      nmtvr,nrcm,levozp,lonr,latr,jcap,num_p3d,num_p2d,npdf3d,    &
     &      kdt,lat,me,pl_coeff,nlons_v,ncw,flgmin,crtrh,cdmbgwd,       &
     &      ccwf,dlqf,ctei_rm,clstp,cgwf,prslrd0,dtp,dtf,fhour,solhr,   &
     &      slag,sdec,cdec,sinlat_v,coslat_v,pgr,gu,gv,                 &
     &      gt,gr,vvel,prsi,prsl,prslk,prsik,phii,phil,                 &
     &      rannum_v,ozplout_v,pl_pres,dpshc,                           &
     &      hprime_v, xlon(lon,lan),xlat(lon,lan),                      &
     &      sfc_fld%slope (lon,lan),    sfc_fld%shdmin(lon,lan),        &
     &      sfc_fld%shdmax(lon,lan),    sfc_fld%snoalb(lon,lan),        &
     &      sfc_fld%tg3   (lon,lan),    sfc_fld%slmsk (lon,lan),        &
     &      sfc_fld%vfrac (lon,lan),    sfc_fld%vtype (lon,lan),        &
     &      sfc_fld%stype (lon,lan),    sfc_fld%uustar(lon,lan),        &
     &      sfc_fld%oro   (lon,lan),    sfc_fld%oro_uf(lon,lan),        &
     &      flx_fld%coszen(lon,lan),                                    &
     &      flx_fld%sfcdsw(lon,lan),    flx_fld%sfcnsw(lon,lan),        &

     &      nirbmdi, nirdfdi, visbmdi, visdfdi,                         &
     &      nirbmui, nirdfui, visbmui, visdfui,                         &

     &      flx_fld%sfcdlw(lon,lan),    flx_fld%tsflw (lon,lan),        &
     &      flx_fld%sfcemis(lon,lan),   sfalb(lon,lan),                 &
     &      swh(1,1,iblk,lan),swhc(1,1,iblk,lan),                       &
     &      hlw(1,1,iblk,lan),hlwc(1,1,iblk,lan), hlwd, lsidea,         &
     &      ras,pre_rad,ldiag3d,lgocart,lssav,cplflx,                   &
     &      bkgd_vdif_m,bkgd_vdif_h,bkgd_vdif_s,psautco,prautco,evpco,  &
     &      wminco,pdfcld,shcnvcw,sup,redrag,hybedmf,dspheat,           &
     &      flipv,old_monin,cnvgwd,shal_cnv,sashal,newsas,cal_pre,      &
     &      mom4ice,mstrat,trans_trac,nst_fcst,moist_adj,               &
     &      thermodyn_id, sfcpress_id, gen_coord_hybrid,levr,adjtrc,nn, &
!  ---  input/outputs:
     &      sfc_fld%hice  (lon,lan),    sfc_fld%fice  (lon,lan),        &
     &      sfc_fld%tisfc (lon,lan),    sfc_fld%tsea  (lon,lan),        &
     &      sfc_fld%tprcp (lon,lan),    sfc_fld%cv    (lon,lan),        &
     &      sfc_fld%cvb   (lon,lan),    sfc_fld%cvt   (lon,lan),        &
     &      sfc_fld%srflag(lon,lan),    sfc_fld%snwdph(lon,lan),        &
     &      sfc_fld%sheleg(lon,lan),    sfc_fld%sncovr(lon,lan),        &
     &      sfc_fld%zorl  (lon,lan),    sfc_fld%canopy(lon,lan),        &
     &      sfc_fld%ffmm  (lon,lan),    sfc_fld%ffhh  (lon,lan),        &
     &      sfc_fld%f10m  (lon,lan),    flx_fld%srunoff(lon,lan),       &
     &      flx_fld%evbsa (lon,lan),    flx_fld%evcwa (lon,lan),        &
     &      flx_fld%snohfa(lon,lan),    flx_fld%transa(lon,lan),        &
     &      flx_fld%sbsnoa(lon,lan),    flx_fld%snowca(lon,lan),        &
     &      flx_fld%soilm (lon,lan),    flx_fld%tmpmin(lon,lan),        &
     &      flx_fld%tmpmax(lon,lan),    flx_fld%dusfc (lon,lan),        &
     &      flx_fld%dvsfc (lon,lan),    flx_fld%dtsfc (lon,lan),        &
     &      flx_fld%dqsfc (lon,lan),    flx_fld%geshem(lon,lan),        &
     &      flx_fld%gflux (lon,lan),    flx_fld%dlwsfc(lon,lan),        &
     &      flx_fld%ulwsfc(lon,lan),    flx_fld%suntim(lon,lan),        &
     &      flx_fld%runoff(lon,lan),    flx_fld%ep    (lon,lan),        &
     &      flx_fld%cldwrk(lon,lan),    flx_fld%dugwd (lon,lan),        &
     &      flx_fld%dvgwd (lon,lan),    flx_fld%psmean(lon,lan),        &
     &      flx_fld%bengsh(lon,lan),    flx_fld%spfhmin(lon,lan),       &
     &      flx_fld%spfhmax(lon,lan),                                   &
     &      flx_fld%rain(lon,lan),      flx_fld%rainc(lon,lan),         &

!    &      dt3dt(1,1,1,iblk,lan), dq3dt(1,1,1,iblk,lan),               &
!    &      du3dt(1,1,1,iblk,lan), dv3dt(1,1,1,iblk,lan), dqdt_v,       & ! added for GOCART
     &      dt3dt_v, dq3dt_v,  du3dt_v, dv3dt_v, dqdt_v,cnvqc_v,        & ! added for GOCART
     &      acv(lon,lan), acvb(lon,lan), acvt(lon,lan),                 &
!    &      slc_v, smc_v, stc_v, upd_mf(1,1,iblk,lan),                  &
!    &      dwn_mf(1,1,iblk,lan), det_mf(1,1,iblk,lan),                 &
     &      slc_v, smc_v, stc_v, upd_mfv, dwn_mfv, det_mfv,             &

     &      phy_f3d(1,1,1,iblk,lan), phy_f2dv,                          &
!    &      phy_f3dv, phy_f2dv,                                         &

     &      aoi_du,     aoi_dv, aoi_dt, aoi_dq, aoi_dlw, aoi_dsw,       &
     &      aoi_dnirbm, aoi_dnirdf, aoi_dvisbm, aoi_dvisdf, aoi_rain,   &
     &      aoi_nlw,    aoi_nsw,    aoi_nnirbm, aoi_nnirdf,             &
     &      aoi_nvisbm, aoi_nvisdf,                                     &

     &      nst_xt,      nst_xs,   nst_xu,   nst_xv,     nst_xz,        &
     &      nst_zm,      nst_xtts, nst_xzts, nst_d_conv, nst_ifd,       &
     &      nst_dt_cool, nst_Qrain,                                     &
!  ---  outputs:
     &      adt, adr, adu, adv,                                         &
     &      sfc_fld%t2m   (lon,lan),    sfc_fld%q2m   (lon,lan),        &
     &      flx_fld%u10m  (lon,lan),    flx_fld%v10m  (lon,lan),        &
     &      flx_fld%zlvl  (lon,lan),    flx_fld%psurf (lon,lan),        &
     &      flx_fld%hpbl  (lon,lan),    flx_fld%pwat  (lon,lan),        &
     &      flx_fld%t1    (lon,lan),    flx_fld%q1    (lon,lan),        &
     &      flx_fld%u1    (lon,lan),    flx_fld%v1    (lon,lan),        &
     &      flx_fld%chh   (lon,lan),    flx_fld%cmm   (lon,lan),        &
     &      flx_fld%dlwsfci(lon,lan),   flx_fld%ulwsfci(lon,lan),       &
     &      flx_fld%dswsfci(lon,lan),   flx_fld%uswsfci(lon,lan),       &
     &      flx_fld%dusfci(lon,lan),    flx_fld%dvsfci(lon,lan),        &
     &      flx_fld%dtsfci(lon,lan),    flx_fld%dqsfci(lon,lan),        &
     &      flx_fld%gfluxi(lon,lan),    flx_fld%epi   (lon,lan),        &
     &      flx_fld%smcwlt2(lon,lan),   flx_fld%smcref2(lon,lan),       &
     &      flx_fld%wet1(lon,lan),      flx_fld%sr(lon,lan),            &
     &      rqtk,                                                       &! rqtkD
     &      dtdt,                                                       &

     &      aoi_dusfci,     aoi_dvsfci,  aoi_dtsfci,   aoi_dqsfci,      &
     &      aoi_dlwsfci,    aoi_dswsfci, aoi_dnirbmi,  aoi_dnirdfi,     &
     &      aoi_dvisbmi,    aoi_dvisdfi, aoi_nlwsfci,  aoi_nswsfci,     &
     &      aoi_nnirbmi,    aoi_nnirdfi, aoi_nvisbmi,  aoi_nvisdfi,     &
     &      aoi_t2mi,       aoi_q2mi,    aoi_u10mi,    aoi_v10mi,       &
     &      aoi_tseai,      aoi_psurfi,                                 &

     &      nst_Tref, nst_z_c, nst_c_0, nst_c_d, nst_w_0, nst_w_d       &
     &      )

            if (nn < nsphys) then
              do k=1,levs
                do i=1,njeff
                  gt(i,k) = adt(i,k)
                  gu(i,k) = adu(i,k)
                  gv(i,k) = adv(i,k)
                enddo
              enddo
              do n=1,ntrac
                do k=1,levs
                  do i=1,njeff
                    gr(i,k,n) = adr(i,k,n)
                  enddo
                enddo
              enddo
            else
              do n=1,ntrac
                if (fixtrc(n)) then
                  do k=1,levs
                    do i=1,njeff
                      trcp(lon+i-1,n,2) = trcp(lon+i-1,n,2)
     &                        + adr(i,k,n) * (prsi(i,k) - prsi(i,k+1))
                    enddo
                  enddo
                endif
              enddo
              if (gg_tracers) then
                do i=1,njeff
                  item = lon+i-1
                  grid_fld%rqtk(item,lan) = rqtk(i) / pgr(i)
                enddo
!             else
!               do i=1,njeff
!                 item = lon+i-1
!                 grid_fld%rqtk(item,lan) = 0.0
!               enddo
              endif
            endif
          enddo                                ! end of nsphys loop
!-------------------------------------------------------------------------

!         if(kdt==100) then
!      print *,'in gloopb,aft gbphys,kdt=',kdt,'lat=',lat,lon,'smcwlt=',
!     &     flx_fld%smcwlt2(lon:lon+3,lan),
!     &    'loc=',minloc(flx_fld%smcwlt2(lon:lon+njeff-1,lan))
!         endif
!
!!
!hmhj debug
!         do n=1,ntrac
!         call mymaxmin(gr(1,1,n),njeff,ngptc,levs,' gr a gbphys ')
!         call mymaxmin(adr(1,1,n),njeff,ngptc,levs,' adr a gbphys ')
!         enddo
!         write(0,*)' in gloopb dusfc=',flx_fld%dusfc(:,lan),' lan=',
!    &lan,' kdt=',kdt

          do k=1,lsoil
            do i=1,njeff
              item = lon+i-1
              sfc_fld%smc(k,item,lan) = smc_v(i,k)
              sfc_fld%stc(k,item,lan) = stc_v(i,k)
              sfc_fld%slc(k,item,lan) = slc_v(i,k)
            enddo
          enddo
!
          if (ldiag3d) then
            do k=1,6
              do j=1,levs
                do i=1,njeff
                  dt3dt(i,j,k,iblk,lan) = dt3dt_v(i,j,k)
                enddo
              enddo
            enddo
            do k=1,4
              do j=1,levs
                do i=1,njeff
                  du3dt(i,j,k,iblk,lan) = du3dt_v(i,j,k)
                  dv3dt(i,j,k,iblk,lan) = dv3dt_v(i,j,k)
                enddo
              enddo
            enddo
            do k=1,5+pl_coeff
              do j=1,levs
                do i=1,njeff
                  dq3dt(i,j,k,iblk,lan) = dq3dt_v(i,j,k)
                enddo
              enddo
            enddo
            do j=1,levs
              do i=1,njeff
                upd_mf(i,j,iblk,lan) = upd_mfv(i,j)
                dwn_mf(i,j,iblk,lan) = dwn_mfv(i,j)
                det_mf(i,j,iblk,lan) = det_mfv(i,j)
              enddo
            enddo
          endif
!
          if (cplflx) then
            do i=1,njeff
              item = lon+i-1
!
              aoi_fld%dusfc(item,lan)   = aoi_du(i)
              aoi_fld%dvsfc(item,lan)   = aoi_dv(i)
              aoi_fld%dtsfc(item,lan)   = aoi_dt(i)
              aoi_fld%dqsfc(item,lan)   = aoi_dq(i)
              aoi_fld%dlwsfc(item,lan)  = aoi_dlw(i)
              aoi_fld%dswsfc(item,lan)  = aoi_dsw(i)
              aoi_fld%dnirbm(item,lan)  = aoi_dnirbm(i)
              aoi_fld%dnirdf(item,lan)  = aoi_dnirdf(i)
              aoi_fld%dvisbm(item,lan)  = aoi_dvisbm(i)
              aoi_fld%dvisdf(item,lan)  = aoi_dvisdf(i)
              aoi_fld%rain(item,lan)    = aoi_rain(i)
              aoi_fld%nlwsfc(item,lan)  = aoi_nlw(i)
              aoi_fld%nswsfc(item,lan)  = aoi_nsw(i)
              aoi_fld%nnirbm(item,lan)  = aoi_nnirbm(i)
              aoi_fld%nnirdf(item,lan)  = aoi_nnirdf(i)
              aoi_fld%nvisbm(item,lan)  = aoi_nvisbm(i)
              aoi_fld%nvisdf(item,lan)  = aoi_nvisdf(i)

              aoi_fld%dusfci(item,lan)  = aoi_dusfci(i)
              aoi_fld%dvsfci(item,lan)  = aoi_dvsfci(i)
              aoi_fld%dtsfci(item,lan)  = aoi_dtsfci(i)
              aoi_fld%dqsfci(item,lan)  = aoi_dqsfci(i)
              aoi_fld%dlwsfci(item,lan) = aoi_dlwsfci(i)
              aoi_fld%dswsfci(item,lan) = aoi_dswsfci(i)
              aoi_fld%dnirbmi(item,lan) = aoi_dnirbmi(i)
              aoi_fld%dvisbmi(item,lan) = aoi_dvisbmi(i)
              aoi_fld%dvisdfi(item,lan) = aoi_dvisdfi(i)

              aoi_fld%nlwsfci(item,lan) = aoi_nlwsfci(i)
              aoi_fld%nswsfci(item,lan) = aoi_nswsfci(i)
              aoi_fld%nirbmi(item,lan)  = aoi_nnirbmi(i)
              aoi_fld%nnirdfi(item,lan) = aoi_nnirdfi(i)
              aoi_fld%nvisbmi(item,lan) = aoi_nvisbmi(i)
              aoi_fld%nvisdfi(item,lan) = aoi_nvisdfi(i)

              aoi_fld%t2mi(item,lan)    = aoi_t2mi(i)
              aoi_fld%q2mi(item,lan)    = aoi_q2mi(i)
              aoi_fld%u10mi(item,lan)   = aoi_u10mi(i)
              aoi_fld%v10mi(item,lan)   = aoi_v10mi(i)
              aoi_fld%tseai(item,lan)   = aoi_tseai(i)
              aoi_fld%psurfi(item,lan)  = aoi_psurfi(i)

              aoi_fld%oro(item,lan)     = sfc_fld%oro(item,lan)
              aoi_fld%slimsk(item,lan)  = sfc_fld%slmsk(item,lan)
            enddo
          endif
          if (nst_fcst > 0) then
            do i=1,njeff
              item = lon+i-1
              nst_fld%xt(item,lan)      = nst_xt(i)
              nst_fld%xs(item,lan)      = nst_xs(i)
              nst_fld%xu(item,lan)      = nst_xu(i)
              nst_fld%xv(item,lan)      = nst_xv(i)
              nst_fld%xz(item,lan)      = nst_xz(i)
              nst_fld%zm(item,lan)      = nst_zm(i)
              nst_fld%xtts(item,lan)    = nst_xtts(i)
              nst_fld%xzts(item,lan)    = nst_xzts(i)
              nst_fld%d_conv(item,lan)  = nst_d_conv(i)
              nst_fld%ifd(item,lan)     = nst_ifd(i)
              nst_fld%dt_cool(item,lan) = nst_dt_cool(i) 
              nst_fld%qrain(item,lan)   = nst_qrain(i)

              nst_fld%tref(item,lan)    = nst_tref(i)
              nst_fld%z_c(item,lan)     = nst_z_c(i)
              nst_fld%c_0(item,lan)     = nst_c_0(i)
              nst_fld%c_d(item,lan)     = nst_c_d(i)
              nst_fld%w_0(item,lan)     = nst_w_0(i)
              nst_fld%w_d(item,lan)     = nst_w_d(i)
            enddo
          endif
!!
!! total moist tendency (kg/kg/s): from local to global array
!!
          if (lgocart) then
            do k=1,levs
              do i=1,njeff
                g3d_fld%dqdt(lon+i-1,lan,k) = dqdt_v(i,k) 
                g3d_fld%cnv_mfc(lon+i-1,lan,k) = 
     &              (upd_mfv(i,k) + dwn_mfv(i,k))/dtf
                g3d_fld%cnv_mfd(lon+i-1,lan,k) = det_mfv(i,k)/dtf
                g3d_fld%cnv_qc(lon+i-1,lan,k) = cnvqc_v(i,k)
              enddo        
            enddo         
          endif
!
!!
          do j=1,num_p2d
            do i=1,njeff
              phy_f2d(lon+i-1,lan,j) = phy_f2dv(i,j)
            enddo
          enddo

          do k=1,levs
            do i=1,njeff
              item = lon+i-1
              grid_fld%u(item,lan,k) = adu(i,k)            
              grid_fld%v(item,lan,k) = adv(i,k)         
              grid_fld%t(item,lan,k) = adt(i,k)
            enddo
          enddo
          do n=1,ntrac
            do k=1,levs
              do i=1,njeff
                grid_fld%tracers(n)%flds(lon+i-1,lan,k) = adr(i,k,n)
              enddo
            enddo
          enddo

!hmhj debug
!     do n=1,ntrac
!       call mymaxmin(adr(1,1,n),njeff,ngptc,levs,' adr a gbphys ')
!     enddo
!!

        enddo                                   !lon
!---------------------------main longitude loop ends--------------------------------
        tem = 0.5 / lonsperlar(lat)
        do n=1,ntrac
          if (fixtrc(n)) then
            trcj(lan,n,1) = 0.0
            trcj(lan,n,2) = 0.0
            do j=1,lons_lat
              trcj(lan,n,1) = trcj(lan,n,1) + trcp(j,n,1)
              trcj(lan,n,2) = trcj(lan,n,2) + trcp(j,n,2)
            enddo
            trcj(lan,n,1) = trcj(lan,n,1) * tem
            trcj(lan,n,2) = trcj(lan,n,2) * tem
          endif
        enddo
!
      enddo                                    !lan
!---------------------------main latitude loop ends---------------------------------
!!
      call excht(lats_nodes_r, global_lats_r, trcj, trcg)
      do n=1,ntrac
        if (fixtrc(n)) then
          sumtrc(n,1) = 0.0
          sumtrc(n,2) = 0.0
          do lat=1,latr
            sumtrc(n,1) = sumtrc(n,1)
     &                    + wgt_r(min(lat,latr-lat+1))*trcg(lat,n,1)
            sumtrc(n,2) = sumtrc(n,2)
     &                    + wgt_r(min(lat,latr-lat+1))*trcg(lat,n,2)
!     write(0,*)' kdt=',kdt,' lat=',lat,' sumtrc=',sumtrc(n,:)
!    &,' trcg=',trcg(lat,n,1),trcg(lat,n,2),' n=',n
          enddo

!     write(0,*)' kdt=',kdt,' n=',n,' sumtrc=',sumtrc(n,:)

          if (abs(sumtrc(n,1)) > 1.0e-12 .and. kdt > 1) then
!         if (abs(sumtrc(n,1)) > 1.0e-12 .and. kdt > 1
!    &                             .and. n == 2) then
            tem = (sumtrc(n,3)-sumtrc(n,1)) / sumtrc(n,1)
            if (abs(tem) < 0.1) then
              adjtrc(n) = 1 + tem
            else
              adjtrc(n) = 1.0
            endif
          else
            adjtrc(n) = 1.0
          endif
          sumtrc(n,3) = sumtrc(n,2)
!       write(0,*)' kdt=',kdt,' n=',n,' adjtrc=',adjtrc(n)
!    &,           ' sumtrc=',sumtrc(n,:)
        endif
      enddo

!
      if (allocated(indxr))    deallocate (indxr)

      return
      end