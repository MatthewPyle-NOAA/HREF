
       SUBROUTINE lpm(nx,ny,n_ens,patch_nx,patch_ny,ovx, ovy, 
     &              filt_min,gauss_sigma,        
     &              var2d_ens,var2d_enmean,     
     &              var2d_lpm,var2d_pmmn)
!
!Author: Nate Snook (28 Apr. 2017)
!
!Purpose: Calculate localized PM mean of a 2D field.  Localized PM mean
!         uses all data in a non-exclusive calculation area to obtain
!         PM mean for a smaller patch within that region of dimensions
!         (patch_nx, patch_ny).  The size of the calculation area is
!         controled by the overlap variables ovx and ovy.  
!
!         Here is a diagram of the design:
!
!         +------------------------------------------------------------+
!         |  full domain                                               |
!         |        +---------------------------+                       |
!         |        | calculation ^             |                       |
!         |        | area        |             |                       |
!         |        |            ovy            |                       |
!         |        |             |             |                       |
!         |        |             v             |                       |
!         |        |         +-------+         |                       |
!         |        |         | patch |         |                       |
!         |        |<--ovx-->|       |<--ovx-->|                       |
!         |        |         +-------+         |                       |
!         |        |             ^             |                       |
!         |        |             |             |                       |
!         |        |            ovy            |                       |
!         |        |             |             |                       |
!         |        |             v             |                       |
!         |        +---------------------------+                       |
!         |                                                            |
!         |                                                            |
!         |                                                            |
!         |                                                            |
!         +------------------------------------------------------------+
!
!         This calculation is repeated for patches in the x- and y-direction
!         such that the entire domain is covered, and then a smoother is run
!         on the resulting field to ensure a smooth 2D output.  For best 
!         results, patch_nx and patch_ny of no larger than 8 are recommended.
!
! History
!
! 6/2/2017 Fanyou Kong
!   - Code revision, bug fix, and clean up
! 6/7/2017 Fanyou Kong
!   - Add a constrain with ensemble max (var2d_enmax)
!

       IMPLICIT NONE

!Inputs
       INTEGER :: nx, ny, n_ens    ! number of x, y gridpoints, ensemble members
       INTEGER :: patch_nx, patch_ny ! size of actual patches in gridpoints (recommended value <= 8)
                            !       For HMT QPF, patch_nx = patch_ny = 6 works rather well.
       INTEGER :: ovx, ovy         ! how many points of patch "overlap" to include in calculation area.
                            !       Recommended value: 3 to 10 times patch_nx and patch_ny
                            !       For HMT QPF, a value of ovx = ovy = 30 works quite well.
       REAL :: filt_min            ! LPM is set to zero where enmean < filt_min  -- this helps to remove
                            !       areas of noise that can occur due to patch overlap.  For QPF,
                            !       filt_min = 0.1 mm for 3hr accumulated precip works well.
       REAL :: gauss_sigma         ! Parameter sigma (standard deviation) to be passed to gauss_rad
                            !       Recommended value = 2 (2 gridpoints) -- this needs some 
                            !       testing
       REAL :: var2d_ens(nx, ny, n_ens)  ! 2D (nx, ny) field for each member (n_ens)
       REAL :: var2d_enmean(nx, ny)      ! 2D (nx, ny) ensemble mean field
       REAL :: var2d_pmmn(nx, ny)         ! 2D (nx, ny) PMMN field
!Outputs
       REAL :: var2d_lpm(nx, ny)         ! 2D (nx, ny) LPM field

!Calculated and used in this code
      INTEGER :: ipatches, jpatches  ! number of patches in i, j directions
      INTEGER :: ipatch, jpatch      ! iterator for patches in i, j directions
      INTEGER :: pe_west, pe_east, pe_north, pe_south ! patch edge locations (i or j point)
      INTEGER :: ce_west, ce_east, ce_north, ce_south ! calculation area edge locations
      INTEGER :: calc_nx, calc_ny    ! The number of x, y gridpoints within the calculation area
      INTEGER :: offset_i, offset_j  ! The distance of the patch edge from the W, S edge of calc area

      REAL, ALLOCATABLE :: lpm_calc(:,:)  !The LPM over the calculation area.  Must be allocated since
                                          !we do not know ahead of time how big the calc area will be.
      REAL, ALLOCATABLE :: var0(:)   ! The values from ALL ensemble members at all points on the 2D
                               ! slice of interest crammed into a 1D array (needed for call to
                               ! pm_mean)
                               ! Since this is limited to calc area, size is unknown ahead of 
                               ! time.
      INTEGER :: ix,jy

!---------------------------------------------------------------!
!Determine whether specified patch_nx, patch_ny divide evenly into the domain
!and if they do not, return an error.
!!! Kong comment: This might not be necessary - need to loosen
!      IF (MOD((ny - 3), patch_ny) /= 0) THEN
!        write(0,*) 'possible trouble(y)..exit maybe'
!       CALL arpsstop('LPM MEAN ERROR: patches do not divide evenly into domain (y dir)', 1)
!      END IF

!      IF (MOD((nx - 3), patch_nx) /= 0) THEN
!        write(0,*) 'possible trouble(x)..exit maybe'
!    CALL arpsstop('LPM MEAN ERROR: patches do not divide evenly into domain (x dir)', 1)
!      END IF

 
!       var2d_lpm=0.

! pre fill with pmmn so pmmn gets used when lpm is not defined along boundaries.
        var2d_lpm=var2d_pmmn

	write(0,*) 'nx, ny: ', nx, ny
	write(0,*) 'patch_nx, patch_ny: ',patch_nx, patch_ny

      ipatches=(nx-3)/float(patch_nx)
      jpatches=(ny-3)/float(patch_ny)

       print *,'LPM -- ipatches,jpatches:',ipatches,jpatches

!Now loop over each patch and perform the PM mean in that patch's calculation area
!   NOTE: where calculation areas extend outside the domain boundary, they are clipped to fit
      DO ipatch=1, ipatches
         DO jpatch=1, jpatches

        !Calculate patch edge locations (i for pe_west, pe_east; j for pe_south, pe_north)
        !   NOTE: the "2 + " is to deal with the non-physical gridpoint on S, W arps domain edges
        !         (it excludes the westmost/southmost point (which has index "1" in FORTRAN)
           pe_west = 2 + ((ipatch-1) * patch_nx)
           pe_south= 2 + ((jpatch-1) * patch_ny)
           pe_east = pe_west + patch_nx
           pe_north= pe_south + patch_ny

        !Calculate nominal edges of calculation area (before clipping, if any is needed)
        ce_west = pe_west - ovx
        ce_east = pe_east + ovx
        ce_south = pe_south - ovx
        ce_north = pe_north + ovx

        !Catch edge cases:
        offset_i = ovx     !If there is no clipping, offset_i is the same as ovx
        offset_j = ovy     !If there is no clipping, offset_j is the same as ovy


!!! is the 2 here arbitrary, or tied to the ARPS thing noted above?

        IF (ce_west < 2) THEN    !West edge case
            offset_i = ovx - (2 - ce_west)
            ce_west = 2
        END IF
       
        IF (ce_east > (nx - 2)) THEN  !East edge case
            ce_east = nx - 2   !No need to change offset in this case (only cares about offset
                               !from west side of calculation area, not east side)
        END IF
       
        IF (ce_south < 2) THEN    !South edge case
            offset_j = ovy - (2 - ce_south)
            ce_south = 2
        END IF
      
        IF (ce_north > (ny - 2)) THEN    ! North edge case
            ce_north = ny - 2   !No need to change offset in this case (only cares about offset
                                !from south side of calcuation area, not north side)
        END IF

        !Define calc_nx, calc_ny:
        calc_nx = ce_east + 1 - ce_west
        calc_ny = ce_north + 1 - ce_south


        !Allocate LPM storage for calculation area and var0 using calc_nx, calc_ny:
        ALLOCATE(lpm_calc(calc_nx, calc_ny))
        ALLOCATE(var0(calc_nx*calc_ny*n_ens))

        !Calculate PM mean over calculation area:
        !First, if all member values in calculation area are zero, then LPM is zero
        IF(maxval(var2d_ens(ce_west:ce_east, ce_south:ce_north, :))
     &                                                <= 0 ) THEN
            lpm_calc = 0.0
        ELSE

!! var0 has the individual ensemble data just for the patch being worked
       var0 = RESHAPE(var2d_ens(ce_west:ce_east,ce_south:ce_north,:),
     &                     (/calc_nx*calc_ny*n_ens/) )

! var0 (as a 1D array)
! var2d_enmean (as patch of full 2D)

            !This call calculates PM mean over the calc area and stores it in lpm_calc:
            CALL pmatch_mean_loc(calc_nx,calc_ny,var0, 
     &                 var2d_enmean(ce_west:ce_east, ce_south:ce_north),
     &                 lpm_calc,ce_west,ce_east,ce_south,ce_north,n_ens)


! could this be shifting things?  It was.  Added +1 to the lpm_calc indices

        !Take the result and store the patch into the array for var2d_lpm:

!	write(*,*) 'shape(lpm_calc): ', shape(lpm_calc)
!        write(*,*) 'offset_i+1, offset_i + patch_nx+1: ', 
!     &              offset_i+1, offset_i + patch_nx+1
!        write(*,*) 'offset_j+1, offset_j+1 + patch_ny: ', 
!     &              offset_j+1, offset_j+1 + patch_ny

        if ( (offset_i + patch_nx+1) .le. calc_nx .and. 
     &       (offset_j+1 + patch_ny ) .le. calc_ny) then
        var2d_lpm(pe_west:pe_east, pe_south:pe_north) =  
     &    lpm_calc(offset_i+1:offset_i + patch_nx+1, 
     &    offset_j+1:offset_j+1 + patch_ny)  

        elseif (offset_i .lt. 1 .and. offset_j .ge. 1) then
        var2d_lpm(pe_west:pe_east, pe_south:pe_north) =  
     &    lpm_calc(offset_i+1:offset_i + patch_nx+1, 
     &    offset_j:offset_j + patch_ny)  

        elseif ( offset_j .lt. 1 .and. offset_i .ge. 1) then
        var2d_lpm(pe_west:pe_east, pe_south:pe_north) =  
     &    lpm_calc(offset_i:offset_i + patch_nx, 
     &    offset_j+1:offset_j+1 + patch_ny)  

         else
        var2d_lpm(pe_west:pe_east, pe_south:pe_north) =  
     &    lpm_calc(offset_i:offset_i + patch_nx, 
     &    offset_j:offset_j + patch_ny)  

         endif
        END IF
        
        !Deallocate values for this patch; will re-allocate again for the next one.
        DEALLOCATE(lpm_calc)
        DEALLOCATE(var0)
        END DO
       END DO

        if( filt_min > 0.0 ) then
 !Filter the LPM so that it is used only above a threshold value
          DO ix=1, nx
            DO jy=1, ny
        !IF(var2d_ens(ix, jy) < filt_min) THEN
        IF(var2d_enmean(ix, jy) < filt_min) THEN
            var2d_lpm(ix, jy) = 0.0
        END IF  
           END DO
        END DO
       endif

! same sort of constraint used in HREF PM mean
! Apply ensemble max as a ceiling (6/7/2017 - Fanyou Kong)

! do not have a var2d_enmax, so do as done in HREF
! var2d_lpm=min(var2d_lpm,var2d_enmax)

!Apply a gaussian smoother to the filtered var2d_lpm to get rid of near-grid-scale noise
! fix later with Gsmooth for HREF


         END SUBROUTINE lpm
!---------------------------------------------------------------!
