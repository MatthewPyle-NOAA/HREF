!STARTOFREGISTRYGENERATEDINCLUDE 'inc/nmm_nest_interpdown_unpack.inc'
!
! WARNING This file is generated automatically by use_registry
! using the data base in the file named Registry.
! Do not edit.  Your changes to this file will be lost.
!
CALL rsl_lite_from_parent_info(pig,pjg,retval)
DO while ( retval .eq. 1 )
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%lu_index(pig,pjg) = xv(1)
DO itrace =  PARAM_FIRST_SCALAR, num_scalar
CALL rsl_lite_from_parent_msg((((ckde-1))-(ckds)+1)*RWORDSIZE,xv)
DO k = ckds,(ckde-1)
scalar(pig,pjg,k,itrace) = xv(k)
ENDDO
ENDDO
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%snowh(pig,pjg) = xv(1)
CALL rsl_lite_from_parent_msg(RWORDSIZE,xv)
grid%rhosn(pig,pjg) = xv(1)
CALL rsl_lite_from_parent_info(pig,pjg,retval)
ENDDO
!ENDOFREGISTRYGENERATEDINCLUDE