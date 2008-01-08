  SUBROUTINE wrf_bdyout ( fid , grid , config_flags, switch , &
                           dryrun, ierr )
    USE module_io
    USE module_wrf_error
    USE module_io_wrf
    USE module_domain
    USE module_state_description
    USE module_configure
    USE module_utility
    IMPLICIT NONE
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101
      integer, parameter  :: WRF_FILE_OPENED_FOR_WRITE            = 102
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103
      integer, parameter  :: WRF_REAL                             = 104
      integer, parameter  :: WRF_DOUBLE                           = 105
      integer, parameter  :: WRF_FLOAT=WRF_REAL
      integer, parameter  :: WRF_INTEGER                          = 106
      integer, parameter  :: WRF_LOGICAL                          = 107
      integer, parameter  :: WRF_COMPLEX                          = 108
      integer, parameter  :: WRF_DOUBLE_COMPLEX                   = 109
      integer, parameter  :: WRF_FILE_OPENED_FOR_UPDATE           = 110
! This bit is for backwards compatibility with old variants of these flags 
! that are still being used in io_grib1 and io_phdf5.  It should be removed!  
      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102
  
!WRF Error and Warning messages (1-999)
!All i/o package-specific status codes you may want to add must be handled by your package (see below)
! WRF handles these and netCDF messages only
  integer, parameter  :: WRF_NO_ERR                  =  0       !no error
  integer, parameter  :: WRF_WARN_FILE_NF            = -1       !file not found, or incomplete
  integer, parameter  :: WRF_WARN_MD_NF              = -2       !metadata not found
  integer, parameter  :: WRF_WARN_TIME_NF            = -3       !timestamp not found
  integer, parameter  :: WRF_WARN_TIME_EOF           = -4       !no more timestamps
  integer, parameter  :: WRF_WARN_VAR_NF             = -5       !variable not found
  integer, parameter  :: WRF_WARN_VAR_EOF            = -6       !no more variables for the current time
  integer, parameter  :: WRF_WARN_TOO_MANY_FILES     = -7       !too many open files
  integer, parameter  :: WRF_WARN_TYPE_MISMATCH      = -8       !data type mismatch
  integer, parameter  :: WRF_WARN_WRITE_RONLY_FILE   = -9       !attempt to write readonly file
  integer, parameter  :: WRF_WARN_READ_WONLY_FILE    = -10      !attempt to read writeonly file
  integer, parameter  :: WRF_WARN_FILE_NOT_OPENED    = -11      !attempt to access unopened file
  integer, parameter  :: WRF_WARN_2DRYRUNS_1VARIABLE = -12      !attempt to do 2 trainings for 1 variable
  integer, parameter  :: WRF_WARN_READ_PAST_EOF      = -13      !attempt to read past EOF
  integer, parameter  :: WRF_WARN_BAD_DATA_HANDLE    = -14      !bad data handle
  integer, parameter  :: WRF_WARN_WRTLEN_NE_DRRUNLEN = -15      !write length not equal to training length
  integer, parameter  :: WRF_WARN_TOO_MANY_DIMS      = -16      !more dimensions requested than training
  integer, parameter  :: WRF_WARN_COUNT_TOO_LONG     = -17      !attempt to read more data than exists
  integer, parameter  :: WRF_WARN_DIMENSION_ERROR    = -18      !input dimension inconsistent
  integer, parameter  :: WRF_WARN_BAD_MEMORYORDER    = -19      !input MemoryOrder not recognized
  integer, parameter  :: WRF_WARN_DIMNAME_REDEFINED  = -20      !a dimension name with 2 different lengths
  integer, parameter  :: WRF_WARN_CHARSTR_GT_LENDATA = -21      !string longer than provided storage
  integer, parameter  :: WRF_WARN_NOTSUPPORTED       = -22      !function not supportable
  integer, parameter  :: WRF_WARN_NOOP               = -23      !package implements this routine as NOOP

!Fatal errors 
  integer, parameter  :: WRF_ERR_FATAL_ALLOCATION_ERROR  = -100 !allocation error
  integer, parameter  :: WRF_ERR_FATAL_DEALLOCATION_ERR  = -101 !dealloc error
  integer, parameter  :: WRF_ERR_FATAL_BAD_FILE_STATUS   = -102 !bad file status


!Package specific errors (1000+)        
!Netcdf status codes
!WRF will accept status codes of 1000+, but it is up to the package to handle
! and return the status to the user.

  integer, parameter  :: WRF_ERR_FATAL_BAD_VARIABLE_DIM  = -1004
  integer, parameter  :: WRF_ERR_FATAL_MDVAR_DIM_NOT_1D  = -1005
  integer, parameter  :: WRF_ERR_FATAL_TOO_MANY_TIMES    = -1006
  integer, parameter  :: WRF_WARN_BAD_DATA_TYPE      = -1007    !this code not in either spec?
  integer, parameter  :: WRF_WARN_FILE_NOT_COMMITTED = -1008    !this code not in either spec?
  integer, parameter  :: WRF_WARN_FILE_OPEN_FOR_READ = -1009
  integer, parameter  :: WRF_IO_NOT_INITIALIZED      = -1010
  integer, parameter  :: WRF_WARN_MD_AFTER_OPEN      = -1011
  integer, parameter  :: WRF_WARN_TOO_MANY_VARIABLES = -1012
  integer, parameter  :: WRF_WARN_DRYRUN_CLOSE       = -1013
  integer, parameter  :: WRF_WARN_DATESTR_BAD_LENGTH = -1014
  integer, parameter  :: WRF_WARN_ZERO_LENGTH_READ   = -1015
  integer, parameter  :: WRF_WARN_DATA_TYPE_NOT_FOUND = -1016
  integer, parameter  :: WRF_WARN_DATESTR_ERROR      = -1017
  integer, parameter  :: WRF_WARN_DRYRUN_READ        = -1018
  integer, parameter  :: WRF_WARN_ZERO_LENGTH_GET    = -1019
  integer, parameter  :: WRF_WARN_ZERO_LENGTH_PUT    = -1020
  integer, parameter  :: WRF_WARN_NETCDF             = -1021    
  integer, parameter  :: WRF_WARN_LENGTH_LESS_THAN_1 = -1022    
  integer, parameter  :: WRF_WARN_MORE_DATA_IN_FILE  = -1023    
  integer, parameter  :: WRF_WARN_DATE_LT_LAST_DATE  = -1024

! For HDF5 only
  integer, parameter  :: WRF_HDF5_ERR_FILE                 = -200
  integer, parameter  :: WRF_HDF5_ERR_MD                   = -201
  integer, parameter  :: WRF_HDF5_ERR_TIME                 = -202
  integer, parameter  :: WRF_HDF5_ERR_TIME_EOF             = -203
  integer, parameter  :: WRF_HDF5_ERR_MORE_DATA_IN_FILE    = -204
  integer, parameter  :: WRF_HDF5_ERR_DATE_LT_LAST_DATE    = -205
  integer, parameter  :: WRF_HDF5_ERR_TOO_MANY_FILES       = -206
  integer, parameter  :: WRF_HDF5_ERR_TYPE_MISMATCH        = -207
  integer, parameter  :: WRF_HDF5_ERR_LENGTH_LESS_THAN_1   = -208
  integer, parameter  :: WRF_HDF5_ERR_WRITE_RONLY_FILE     = -209
  integer, parameter  :: WRF_HDF5_ERR_READ_WONLY_FILE      = -210
  integer, parameter  :: WRF_HDF5_ERR_FILE_NOT_OPENED      = -211
  integer, parameter  :: WRF_HDF5_ERR_DATESTR_ERROR        = -212
  integer, parameter  :: WRF_HDF5_ERR_DRYRUN_READ          = -213
  integer, parameter  :: WRF_HDF5_ERR_ZERO_LENGTH_GET      = -214
  integer, parameter  :: WRF_HDF5_ERR_ZERO_LENGTH_PUT      = -215
  integer, parameter  :: WRF_HDF5_ERR_2DRYRUNS_1VARIABLE   = -216
  integer, parameter  :: WRF_HDF5_ERR_DATA_TYPE_NOTFOUND   = -217
  integer, parameter  :: WRF_HDF5_ERR_READ_PAST_EOF        = -218
  integer, parameter  :: WRF_HDF5_ERR_BAD_DATA_HANDLE      = -219
  integer, parameter  :: WRF_HDF5_ERR_WRTLEN_NE_DRRUNLEN   = -220
  integer, parameter  :: WRF_HDF5_ERR_DRYRUN_CLOSE         = -221
  integer, parameter  :: WRF_HDF5_ERR_DATESTR_BAD_LENGTH   = -222
  integer, parameter  :: WRF_HDF5_ERR_ZERO_LENGTH_READ     = -223
  integer, parameter  :: WRF_HDF5_ERR_TOO_MANY_DIMS        = -224
  integer, parameter  :: WRF_HDF5_ERR_TOO_MANY_VARIABLES   = -225
  integer, parameter  :: WRF_HDF5_ERR_COUNT_TOO_LONG       = -226
  integer, parameter  :: WRF_HDF5_ERR_DIMENSION_ERROR      = -227
  integer, parameter  :: WRF_HDF5_ERR_BAD_MEMORYORDER      = -228
  integer, parameter  :: WRF_HDF5_ERR_DIMNAME_REDEFINED    = -229
  integer, parameter  :: WRF_HDF5_ERR_MD_AFTER_OPEN        = -230
  integer, parameter  :: WRF_HDF5_ERR_CHARSTR_GT_LENDATA   = -231
  integer, parameter  :: WRF_HDF5_ERR_BAD_DATA_TYPE        = -232
  integer, parameter  :: WRF_HDF5_ERR_FILE_NOT_COMMITTED   = -233

  integer, parameter  :: WRF_HDF5_ERR_ALLOCATION        = -2001
  integer, parameter  :: WRF_HDF5_ERR_DEALLOCATION      = -2002
  integer, parameter  :: WRF_HDF5_ERR_BAD_FILE_STATUS   = -2003
  integer, parameter  :: WRF_HDF5_ERR_BAD_VARIABLE_DIM  = -2004
  integer, parameter  :: WRF_HDF5_ERR_MDVAR_DIM_NOT_1D  = -2005
  integer, parameter  :: WRF_HDF5_ERR_TOO_MANY_TIMES    = -2006
  integer, parameter ::  WRF_HDF5_ERR_DATA_ID_NOTFOUND  = -2007

  integer, parameter ::  WRF_HDF5_ERR_DATASPACE         = -300
  integer, parameter ::  WRF_HDF5_ERR_DATATYPE          = -301
  integer, parameter :: WRF_HDF5_ERR_PROPERTY_LIST      = -302

  integer, parameter :: WRF_HDF5_ERR_DATASET_CREATE     = -303
  integer, parameter :: WRF_HDF5_ERR_DATASET_READ       = -304
  integer, parameter :: WRF_HDF5_ERR_DATASET_WRITE      = -305
  integer, parameter :: WRF_HDF5_ERR_DATASET_OPEN       = -306
  integer, parameter :: WRF_HDF5_ERR_DATASET_GENERAL    = -307
  integer, parameter :: WRF_HDF5_ERR_GROUP              = -308

  integer, parameter :: WRF_HDF5_ERR_FILE_OPEN          = -309
  integer, parameter :: WRF_HDF5_ERR_FILE_CREATE        = -310
  integer, parameter :: WRF_HDF5_ERR_DATASET_CLOSE      = -311
  integer, parameter :: WRF_HDF5_ERR_FILE_CLOSE         = -312
  integer, parameter :: WRF_HDF5_ERR_CLOSE_GENERAL      = -313

  integer, parameter :: WRF_HDF5_ERR_ATTRIBUTE_CREATE   = -314
  integer, parameter :: WRF_HDF5_ERR_ATTRIBUTE_READ     = -315
  integer, parameter :: WRF_HDF5_ERR_ATTRIBUTE_WRITE    = -316
  integer, parameter :: WRF_HDF5_ERR_ATTRIBUTE_OPEN     = -317
  integer, parameter :: WRF_HDF5_ERR_ATTRIBUTE_GENERAL  = -318
  integer, parameter :: WRF_HDF5_ERR_ATTRIBUTE_CLOSE    = -319

  integer, parameter :: WRF_HDF5_ERR_OTHERS             = -320
  integer, parameter :: WRF_HDF5_ERR_ATTRIBUTE_OTHERS   = -321

    TYPE(domain) :: grid
    TYPE(grid_config_rec_type),  INTENT(INOUT)    :: config_flags
    INTEGER, INTENT(IN) :: fid, switch
    INTEGER, INTENT(INOUT) :: ierr
    LOGICAL, INTENT(IN) :: dryrun

    ! Local data
    INTEGER ids , ide , jds , jde , kds , kde , &
            ims , ime , jms , jme , kms , kme , &
            ips , ipe , jps , jpe , kps , kpe

    INTEGER       itrace
    INTEGER , DIMENSION(3) :: domain_start , domain_end
    INTEGER , DIMENSION(3) :: memory_start , memory_end
    INTEGER , DIMENSION(3) :: patch_start , patch_end
    INTEGER i,j
    INTEGER julyr, julday, idt, iswater , map_proj
    REAL    gmt, cen_lat, cen_lon, bdyfrq , truelat1 , truelat2
    INTEGER dyn_opt, diff_opt, km_opt, damp_opt,  &
            mp_physics, ra_lw_physics, ra_sw_physics, sf_sfclay_physics, &
            sf_surface_physics, bl_pbl_physics, cu_physics
    REAL    khdif, kvdif
    INTEGER rc

    CHARACTER*256 message
    CHARACTER*80  char_junk
    INTEGER    ibuf(1)
    REAL       rbuf(1)
    CHARACTER*40            :: next_datestr

    CALL get_ijk_from_grid (  grid ,                        &
                              ids, ide, jds, jde, kds, kde,    &
                              ims, ime, jms, jme, kms, kme,    &
                              ips, ipe, jps, jpe, kps, kpe    )

    call nl_get_dyn_opt ( 1 , dyn_opt                       )

    ! note that the string current_date comes in through use association
    ! of module_io_wrf

! generated by the registry
!STARTOFREGISTRYGENERATEDINCLUDE 'inc/wrf_bdyout.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'PD_BXS'               , &  ! Data Name 
                       grid%nmm_pd_b(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XS'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bdy_width'               , &  ! Dimname 2 
                       'one_element'               , &  ! Dimname 3 
                       'bdy Mass at I,J in the sigma domain'               , &  ! Desc  
                       'Pa'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field PD_BXS memorder XS' , & ! Debug message
1, (jde-1), 1, config_flags%spec_bdy_width, 1, 1, &
1, MAX( ide , jde ), 1, config_flags%spec_bdy_width, 1, 1, &
jps, MIN( (jde-1), jpe ), 1, config_flags%spec_bdy_width, 1, 1, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'PD_BXE'               , &  ! Data Name 
                       grid%nmm_pd_b(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XE'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bdy_width'               , &  ! Dimname 2 
                       'one_element'               , &  ! Dimname 3 
                       'bdy Mass at I,J in the sigma domain'               , &  ! Desc  
                       'Pa'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field PD_BXE memorder XE' , & ! Debug message
1, (jde-1), 1, config_flags%spec_bdy_width, 1, 1, &
1, MAX( ide , jde ), 1, config_flags%spec_bdy_width, 1, 1, &
jps, MIN( (jde-1), jpe ), 1, config_flags%spec_bdy_width, 1, 1, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'PD_BYS'               , &  ! Data Name 
                       grid%nmm_pd_b(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YS'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bdy_width'               , &  ! Dimname 2 
                       'one_element'               , &  ! Dimname 3 
                       'bdy Mass at I,J in the sigma domain'               , &  ! Desc  
                       'Pa'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field PD_BYS memorder YS' , & ! Debug message
1, (ide-1), 1, config_flags%spec_bdy_width, 1, 1, &
1, MAX( ide , jde ), 1, config_flags%spec_bdy_width, 1, 1, &
ips, MIN( (ide-1), ipe ), 1, config_flags%spec_bdy_width, 1, 1, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'PD_BYE'               , &  ! Data Name 
                       grid%nmm_pd_b(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YE'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bdy_width'               , &  ! Dimname 2 
                       'one_element'               , &  ! Dimname 3 
                       'bdy Mass at I,J in the sigma domain'               , &  ! Desc  
                       'Pa'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field PD_BYE memorder YE' , & ! Debug message
1, (ide-1), 1, config_flags%spec_bdy_width, 1, 1, &
1, MAX( ide , jde ), 1, config_flags%spec_bdy_width, 1, 1, &
ips, MIN( (ide-1), ipe ), 1, config_flags%spec_bdy_width, 1, 1, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'PD_BTXS'               , &  ! Data Name 
                       grid%nmm_pd_bt(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XS'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bdy_width'               , &  ! Dimname 2 
                       'one_element'               , &  ! Dimname 3 
                       'bdy tend Mass at I,J in the sigma domain'               , &  ! Desc  
                       '(Pa)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field PD_BTXS memorder XS' , & ! Debug message
1, (jde-1), 1, config_flags%spec_bdy_width, 1, 1, &
1, MAX( ide , jde ), 1, config_flags%spec_bdy_width, 1, 1, &
jps, MIN( (jde-1), jpe ), 1, config_flags%spec_bdy_width, 1, 1, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'PD_BTXE'               , &  ! Data Name 
                       grid%nmm_pd_bt(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XE'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bdy_width'               , &  ! Dimname 2 
                       'one_element'               , &  ! Dimname 3 
                       'bdy tend Mass at I,J in the sigma domain'               , &  ! Desc  
                       '(Pa)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field PD_BTXE memorder XE' , & ! Debug message
1, (jde-1), 1, config_flags%spec_bdy_width, 1, 1, &
1, MAX( ide , jde ), 1, config_flags%spec_bdy_width, 1, 1, &
jps, MIN( (jde-1), jpe ), 1, config_flags%spec_bdy_width, 1, 1, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'PD_BTYS'               , &  ! Data Name 
                       grid%nmm_pd_bt(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YS'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bdy_width'               , &  ! Dimname 2 
                       'one_element'               , &  ! Dimname 3 
                       'bdy tend Mass at I,J in the sigma domain'               , &  ! Desc  
                       '(Pa)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field PD_BTYS memorder YS' , & ! Debug message
1, (ide-1), 1, config_flags%spec_bdy_width, 1, 1, &
1, MAX( ide , jde ), 1, config_flags%spec_bdy_width, 1, 1, &
ips, MIN( (ide-1), ipe ), 1, config_flags%spec_bdy_width, 1, 1, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'PD_BTYE'               , &  ! Data Name 
                       grid%nmm_pd_bt(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YE'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bdy_width'               , &  ! Dimname 2 
                       'one_element'               , &  ! Dimname 3 
                       'bdy tend Mass at I,J in the sigma domain'               , &  ! Desc  
                       '(Pa)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field PD_BTYE memorder YE' , & ! Debug message
1, (ide-1), 1, config_flags%spec_bdy_width, 1, 1, &
1, MAX( ide , jde ), 1, config_flags%spec_bdy_width, 1, 1, &
ips, MIN( (ide-1), ipe ), 1, config_flags%spec_bdy_width, 1, 1, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'T_BXS'               , &  ! Data Name 
                       grid%nmm_t_b(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Sensible temperature'               , &  ! Desc  
                       'K'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field T_BXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'T_BXE'               , &  ! Data Name 
                       grid%nmm_t_b(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Sensible temperature'               , &  ! Desc  
                       'K'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field T_BXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'T_BYS'               , &  ! Data Name 
                       grid%nmm_t_b(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Sensible temperature'               , &  ! Desc  
                       'K'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field T_BYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'T_BYE'               , &  ! Data Name 
                       grid%nmm_t_b(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Sensible temperature'               , &  ! Desc  
                       'K'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field T_BYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'T_BTXS'               , &  ! Data Name 
                       grid%nmm_t_bt(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Sensible temperature'               , &  ! Desc  
                       '(K)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field T_BTXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'T_BTXE'               , &  ! Data Name 
                       grid%nmm_t_bt(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Sensible temperature'               , &  ! Desc  
                       '(K)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field T_BTXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'T_BTYS'               , &  ! Data Name 
                       grid%nmm_t_bt(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Sensible temperature'               , &  ! Desc  
                       '(K)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field T_BTYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'T_BTYE'               , &  ! Data Name 
                       grid%nmm_t_bt(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Sensible temperature'               , &  ! Desc  
                       '(K)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field T_BTYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q_BXS'               , &  ! Data Name 
                       grid%nmm_q_b(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Specific humidity'               , &  ! Desc  
                       'kg kg-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q_BXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q_BXE'               , &  ! Data Name 
                       grid%nmm_q_b(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Specific humidity'               , &  ! Desc  
                       'kg kg-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q_BXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q_BYS'               , &  ! Data Name 
                       grid%nmm_q_b(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Specific humidity'               , &  ! Desc  
                       'kg kg-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q_BYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q_BYE'               , &  ! Data Name 
                       grid%nmm_q_b(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Specific humidity'               , &  ! Desc  
                       'kg kg-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q_BYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q_BTXS'               , &  ! Data Name 
                       grid%nmm_q_bt(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Specific humidity'               , &  ! Desc  
                       '(kg kg-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q_BTXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q_BTXE'               , &  ! Data Name 
                       grid%nmm_q_bt(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Specific humidity'               , &  ! Desc  
                       '(kg kg-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q_BTXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q_BTYS'               , &  ! Data Name 
                       grid%nmm_q_bt(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Specific humidity'               , &  ! Desc  
                       '(kg kg-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q_BTYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q_BTYE'               , &  ! Data Name 
                       grid%nmm_q_bt(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Specific humidity'               , &  ! Desc  
                       '(kg kg-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q_BTYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'U_BXS'               , &  ! Data Name 
                       grid%nmm_u_b(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy U component of wind'               , &  ! Desc  
                       'm s-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field U_BXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'U_BXE'               , &  ! Data Name 
                       grid%nmm_u_b(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy U component of wind'               , &  ! Desc  
                       'm s-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field U_BXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'U_BYS'               , &  ! Data Name 
                       grid%nmm_u_b(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy U component of wind'               , &  ! Desc  
                       'm s-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field U_BYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'U_BYE'               , &  ! Data Name 
                       grid%nmm_u_b(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy U component of wind'               , &  ! Desc  
                       'm s-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field U_BYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'U_BTXS'               , &  ! Data Name 
                       grid%nmm_u_bt(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend U component of wind'               , &  ! Desc  
                       '(m s-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field U_BTXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'U_BTXE'               , &  ! Data Name 
                       grid%nmm_u_bt(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend U component of wind'               , &  ! Desc  
                       '(m s-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field U_BTXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'U_BTYS'               , &  ! Data Name 
                       grid%nmm_u_bt(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend U component of wind'               , &  ! Desc  
                       '(m s-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field U_BTYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'U_BTYE'               , &  ! Data Name 
                       grid%nmm_u_bt(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend U component of wind'               , &  ! Desc  
                       '(m s-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field U_BTYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'V_BXS'               , &  ! Data Name 
                       grid%nmm_v_b(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy V component of wind'               , &  ! Desc  
                       'm s-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field V_BXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'V_BXE'               , &  ! Data Name 
                       grid%nmm_v_b(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy V component of wind'               , &  ! Desc  
                       'm s-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field V_BXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'V_BYS'               , &  ! Data Name 
                       grid%nmm_v_b(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy V component of wind'               , &  ! Desc  
                       'm s-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field V_BYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'V_BYE'               , &  ! Data Name 
                       grid%nmm_v_b(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy V component of wind'               , &  ! Desc  
                       'm s-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field V_BYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'V_BTXS'               , &  ! Data Name 
                       grid%nmm_v_bt(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend V component of wind'               , &  ! Desc  
                       '(m s-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field V_BTXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'V_BTXE'               , &  ! Data Name 
                       grid%nmm_v_bt(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend V component of wind'               , &  ! Desc  
                       '(m s-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field V_BTXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'V_BTYS'               , &  ! Data Name 
                       grid%nmm_v_bt(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend V component of wind'               , &  ! Desc  
                       '(m s-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field V_BTYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'V_BTYE'               , &  ! Data Name 
                       grid%nmm_v_bt(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend V component of wind'               , &  ! Desc  
                       '(m s-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field V_BTYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q2_BXS'               , &  ! Data Name 
                       grid%nmm_q2_b(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy 2 * Turbulence kinetic energy'               , &  ! Desc  
                       'm2 s-2'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q2_BXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q2_BXE'               , &  ! Data Name 
                       grid%nmm_q2_b(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy 2 * Turbulence kinetic energy'               , &  ! Desc  
                       'm2 s-2'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q2_BXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q2_BYS'               , &  ! Data Name 
                       grid%nmm_q2_b(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy 2 * Turbulence kinetic energy'               , &  ! Desc  
                       'm2 s-2'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q2_BYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q2_BYE'               , &  ! Data Name 
                       grid%nmm_q2_b(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy 2 * Turbulence kinetic energy'               , &  ! Desc  
                       'm2 s-2'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q2_BYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q2_BTXS'               , &  ! Data Name 
                       grid%nmm_q2_bt(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend 2 * Turbulence kinetic energy'               , &  ! Desc  
                       '(m2 s-2)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q2_BTXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q2_BTXE'               , &  ! Data Name 
                       grid%nmm_q2_bt(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend 2 * Turbulence kinetic energy'               , &  ! Desc  
                       '(m2 s-2)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q2_BTXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q2_BTYS'               , &  ! Data Name 
                       grid%nmm_q2_bt(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend 2 * Turbulence kinetic energy'               , &  ! Desc  
                       '(m2 s-2)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q2_BTYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'Q2_BTYE'               , &  ! Data Name 
                       grid%nmm_q2_bt(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend 2 * Turbulence kinetic energy'               , &  ! Desc  
                       '(m2 s-2)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field Q2_BTYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'CWM_BXS'               , &  ! Data Name 
                       grid%nmm_cwm_b(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Total condensate'               , &  ! Desc  
                       'kg kg-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field CWM_BXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'CWM_BXE'               , &  ! Data Name 
                       grid%nmm_cwm_b(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Total condensate'               , &  ! Desc  
                       'kg kg-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field CWM_BXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'CWM_BYS'               , &  ! Data Name 
                       grid%nmm_cwm_b(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Total condensate'               , &  ! Desc  
                       'kg kg-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field CWM_BYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'CWM_BYE'               , &  ! Data Name 
                       grid%nmm_cwm_b(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy Total condensate'               , &  ! Desc  
                       'kg kg-1'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field CWM_BYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'CWM_BTXS'               , &  ! Data Name 
                       grid%nmm_cwm_bt(1,kds,1,1)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Total condensate'               , &  ! Desc  
                       '(kg kg-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field CWM_BTXS memorder XSZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'CWM_BTXE'               , &  ! Data Name 
                       grid%nmm_cwm_bt(1,kds,1,2)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'XEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Total condensate'               , &  ! Desc  
                       '(kg kg-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field CWM_BTXE memorder XEZ' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'CWM_BTYS'               , &  ! Data Name 
                       grid%nmm_cwm_bt(1,kds,1,3)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YSZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Total condensate'               , &  ! Desc  
                       '(kg kg-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field CWM_BTYS memorder YSZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
IF ( grid%dyn_opt .EQ. dyn_nmm ) THEN
CALL wrf_ext_write_field (  &
                       fid                , &  ! DataHandle 
                       current_date(1:19) , &  ! DateStr 
                       'CWM_BTYE'               , &  ! Data Name 
                       grid%nmm_cwm_bt(1,kds,1,4)     , &  ! Field 
                       WRF_FLOAT          , &  ! FieldType 
                       grid%communicator , &  ! Comm
                       grid%iocommunicator , &  ! Comm
                       grid%domdesc      , &  ! Comm
                       grid%bdy_mask     , &  ! bdy_mask
                       dryrun             , &  ! flag
                       'YEZ'               , &  ! MemoryOrder
                       ''               , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
                       'bdy tend Total condensate'               , &  ! Desc  
                       '(kg kg-1)/dt'               , &  ! Units 
'inc/wrf_bdyout.inc ext_write_field CWM_BTYE memorder YEZ' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                       ierr )
END IF
DO itrace = PARAM_FIRST_SCALAR , num_scalar
  IF (BTEST(scalar_stream_table(grid%id, itrace ) , switch )) THEN
    CALL wrf_ext_write_field (  &
          fid                             , &  ! DataHandle
          current_date(1:19)              , &  ! DateStr
          TRIM(scalar_dname_table( grid%id, itrace )) // '_BXS', & !data name
          grid%scalar_B(1,kds,1,1,itrace)  , &  ! Field
                       WRF_FLOAT             , &  ! FieldType 
          grid%communicator  , &  ! Comm
          grid%iocommunicator  , &  ! Comm
          grid%domdesc       , &  ! Comm
          grid%bdy_mask       , &  ! bdy_mask
          dryrun             , &  ! flag
          'XSZ'               , &  ! MemoryOrder
          ''                , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
          scalar_desc_table( grid%id, itrace  ), & ! Desc
          scalar_units_table( grid%id, itrace  ), & ! Units
'inc/wrf_bdyout.inc ext_write_field '//TRIM(scalar_dname_table( grid%id, itrace ))//' memorder XZY' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                         ierr )
    CALL wrf_ext_write_field (  &
          fid                             , &  ! DataHandle
          current_date(1:19)              , &  ! DateStr
          TRIM(scalar_dname_table( grid%id, itrace )) // '_BXE', & !data name
          grid%scalar_B(1,kds,1,2,itrace)  , &  ! Field
                       WRF_FLOAT             , &  ! FieldType 
          grid%communicator  , &  ! Comm
          grid%iocommunicator  , &  ! Comm
          grid%domdesc       , &  ! Comm
          grid%bdy_mask       , &  ! bdy_mask
          dryrun             , &  ! flag
          'XEZ'               , &  ! MemoryOrder
          ''                , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
          scalar_desc_table( grid%id, itrace  ), & ! Desc
          scalar_units_table( grid%id, itrace  ), & ! Units
'inc/wrf_bdyout.inc ext_write_field '//TRIM(scalar_dname_table( grid%id, itrace ))//' memorder XZY' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                         ierr )
    CALL wrf_ext_write_field (  &
          fid                             , &  ! DataHandle
          current_date(1:19)              , &  ! DateStr
          TRIM(scalar_dname_table( grid%id, itrace )) // '_BYS', & !data name
          grid%scalar_B(1,kds,1,3,itrace)  , &  ! Field
                       WRF_FLOAT             , &  ! FieldType 
          grid%communicator  , &  ! Comm
          grid%iocommunicator  , &  ! Comm
          grid%domdesc       , &  ! Comm
          grid%bdy_mask       , &  ! bdy_mask
          dryrun             , &  ! flag
          'YSZ'               , &  ! MemoryOrder
          ''                , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
          scalar_desc_table( grid%id, itrace  ), & ! Desc
          scalar_units_table( grid%id, itrace  ), & ! Units
'inc/wrf_bdyout.inc ext_write_field '//TRIM(scalar_dname_table( grid%id, itrace ))//' memorder XZY' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                         ierr )
    CALL wrf_ext_write_field (  &
          fid                             , &  ! DataHandle
          current_date(1:19)              , &  ! DateStr
          TRIM(scalar_dname_table( grid%id, itrace )) // '_BYE', & !data name
          grid%scalar_B(1,kds,1,4,itrace)  , &  ! Field
                       WRF_FLOAT             , &  ! FieldType 
          grid%communicator  , &  ! Comm
          grid%iocommunicator  , &  ! Comm
          grid%domdesc       , &  ! Comm
          grid%bdy_mask       , &  ! bdy_mask
          dryrun             , &  ! flag
          'YEZ'               , &  ! MemoryOrder
          ''                , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
          scalar_desc_table( grid%id, itrace  ), & ! Desc
          scalar_units_table( grid%id, itrace  ), & ! Units
'inc/wrf_bdyout.inc ext_write_field '//TRIM(scalar_dname_table( grid%id, itrace ))//' memorder XZY' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                         ierr )
  ENDIF
ENDDO
DO itrace = PARAM_FIRST_SCALAR , num_scalar
  IF (BTEST(scalar_stream_table(grid%id, itrace ) , switch )) THEN
    CALL wrf_ext_write_field (  &
          fid                             , &  ! DataHandle
          current_date(1:19)              , &  ! DateStr
          TRIM(scalar_dname_table( grid%id, itrace )) // '_BTXS', & !data name
          grid%scalar_BT(1,kds,1,1,itrace)  , &  ! Field
                       WRF_FLOAT             , &  ! FieldType 
          grid%communicator  , &  ! Comm
          grid%iocommunicator  , &  ! Comm
          grid%domdesc       , &  ! Comm
          grid%bdy_mask       , &  ! bdy_mask
          dryrun             , &  ! flag
          'XSZ'               , &  ! MemoryOrder
          ''                , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
          scalar_desc_table( grid%id, itrace  ), & ! Desc
          scalar_units_table( grid%id, itrace  ), & ! Units
'inc/wrf_bdyout.inc ext_write_field '//TRIM(scalar_dname_table( grid%id, itrace ))//' memorder XZY' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                         ierr )
    CALL wrf_ext_write_field (  &
          fid                             , &  ! DataHandle
          current_date(1:19)              , &  ! DateStr
          TRIM(scalar_dname_table( grid%id, itrace )) // '_BTXE', & !data name
          grid%scalar_BT(1,kds,1,2,itrace)  , &  ! Field
                       WRF_FLOAT             , &  ! FieldType 
          grid%communicator  , &  ! Comm
          grid%iocommunicator  , &  ! Comm
          grid%domdesc       , &  ! Comm
          grid%bdy_mask       , &  ! bdy_mask
          dryrun             , &  ! flag
          'XEZ'               , &  ! MemoryOrder
          ''                , &  ! Stagger
                       'south_north'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
          scalar_desc_table( grid%id, itrace  ), & ! Desc
          scalar_units_table( grid%id, itrace  ), & ! Units
'inc/wrf_bdyout.inc ext_write_field '//TRIM(scalar_dname_table( grid%id, itrace ))//' memorder XZY' , & ! Debug message
1, (jde-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
jps, MIN( (jde-1), jpe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                         ierr )
    CALL wrf_ext_write_field (  &
          fid                             , &  ! DataHandle
          current_date(1:19)              , &  ! DateStr
          TRIM(scalar_dname_table( grid%id, itrace )) // '_BTYS', & !data name
          grid%scalar_BT(1,kds,1,3,itrace)  , &  ! Field
                       WRF_FLOAT             , &  ! FieldType 
          grid%communicator  , &  ! Comm
          grid%iocommunicator  , &  ! Comm
          grid%domdesc       , &  ! Comm
          grid%bdy_mask       , &  ! bdy_mask
          dryrun             , &  ! flag
          'YSZ'               , &  ! MemoryOrder
          ''                , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
          scalar_desc_table( grid%id, itrace  ), & ! Desc
          scalar_units_table( grid%id, itrace  ), & ! Units
'inc/wrf_bdyout.inc ext_write_field '//TRIM(scalar_dname_table( grid%id, itrace ))//' memorder XZY' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                         ierr )
    CALL wrf_ext_write_field (  &
          fid                             , &  ! DataHandle
          current_date(1:19)              , &  ! DateStr
          TRIM(scalar_dname_table( grid%id, itrace )) // '_BTYE', & !data name
          grid%scalar_BT(1,kds,1,4,itrace)  , &  ! Field
                       WRF_FLOAT             , &  ! FieldType 
          grid%communicator  , &  ! Comm
          grid%iocommunicator  , &  ! Comm
          grid%domdesc       , &  ! Comm
          grid%bdy_mask       , &  ! bdy_mask
          dryrun             , &  ! flag
          'YEZ'               , &  ! MemoryOrder
          ''                , &  ! Stagger
                       'west_east'               , &  ! Dimname 1 
                       'bottom_top'               , &  ! Dimname 2 
                       'bdy_width'               , &  ! Dimname 3 
          scalar_desc_table( grid%id, itrace  ), & ! Desc
          scalar_units_table( grid%id, itrace  ), & ! Units
'inc/wrf_bdyout.inc ext_write_field '//TRIM(scalar_dname_table( grid%id, itrace ))//' memorder XZY' , & ! Debug message
1, (ide-1), kds, (kde-1), 1, config_flags%spec_bdy_width, &
1, MAX( ide , jde ), kds, kde, 1, config_flags%spec_bdy_width, &
ips, MIN( (ide-1), ipe ), kds, (kde-1), 1, config_flags%spec_bdy_width, &
                         ierr )
  ENDIF
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE

    RETURN
    END