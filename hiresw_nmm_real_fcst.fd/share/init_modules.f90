!WRF:MEDIATION_LAYER
!
SUBROUTINE init_modules( phase )
 USE module_bc
 USE module_configure
 USE module_driver_constants
 USE module_model_constants
 USE module_domain
 USE module_machine
 USE module_nesting
 USE module_timing
 USE module_tiles
 USE module_io_wrf
 USE module_io
 USE module_wrf_quilt
 USE module_dm
 USE module_ext_internal

! <DESCRIPTION>
! This routine USES the modules in WRF and then calls the init routines
! they provide to perform module specific initializations at the
! beginning of a run.  Note, this is only once per run, not once per
! domain; domain specific initializations should be handled elsewhere,
! such as in <a href=start_domain.html>start_domain</a>.
! 
! Certain framework specific module initializations in this file are
! dependent on order they are called. For example, since the quilt module
! relies on internal I/O, the init routine for internal I/O must be
! called first.  In the case of 1 compiles, the quilt module
! calls MPI_INIT as part of setting up and dividing communicators between
! compute and I/O server tasks.  Therefore, it must be called prior to
! module_dm, which will <em>also</em> try to call MPI_INIT if it sees
! that MPI has not be initialized yet (implementations of module_dm
! should in fact behave this way by first calling MPI_INITIALIZED before
! they try to call MPI_INIT).  If MPI is already initialized before the
! the quilting module is called, quilting will not work.
! 
! The phase argument is used to allow other superstructures like ESMF to 
! place their initialization calls following the WRF initialization call 
! that calls MPI_INIT().  When used with ESMF, ESMF will call wrf_init() 
! which in turn will call phase 2 of this routine.  Phase 1 will be called 
! earlier.  
!
! </DESCRIPTION>

 INTEGER, INTENT(IN) :: phase    ! phase==1 means return after MPI_INIT()
                                 ! phase==2 means resume after MPI_INIT()

IF ( phase == 1 ) THEN
 CALL init_module_bc
 CALL init_module_configure
 CALL init_module_driver_constants
 CALL init_module_model_constants
 CALL init_module_domain
 CALL init_module_machine

 CALL init_module_ext_internal  !! must be called before quilt
 CALL init_module_wrf_quilt    !! this *must* be called before init_module_dm
 CALL init_module_dm
ELSE
 CALL init_module_nesting
 CALL init_module_timing
 CALL init_module_tiles
 CALL init_module_io_wrf
 CALL init_module_io

! core specific initializations -- add new cores here
 CALL init_modules_nmm
ENDIF
 
END SUBROUTINE init_modules
