#!/bin/sh
#
#  This makefile is to compile ensembel production generator
#  Author: Binbin Zhou
#          Aug 1, 2005
#  OS: Linux (on WCOSS)
#


INC=  ${G2_INC4} 
LIBS=   ${IP_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4} ${SP_LIB4} ${G2_LIB4}  \
        ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}
FFLAGS = -O3 -fp-model precise  -g -traceback -convert big_endian -assume noold_ldout_format -I $(INC)  #  -ftrapuv  -check all  -fp-stack-check -check bounds

SRCS=	char_process.f util.f readGB2.f print_gfld.f Zputgb2.f \
	packGB2_mean.f packGB2_prob.f packGB2_mean_derv.f packGB2_prob_derv.f \
	preciptype.f wind.f bulkshear.f get_icing.f get_cat.f getflight.f Hains_index.f getceiling.f \
	new_fog.f get_new_fog.f thickness.f llws.f getconv.f \
	get_mxp.f mxp.f packGB2_mxp.f packGB2_max.f packGB2_spag.f \
	get_cptp_severe.f calibrate2_hrly_rgn3.f calibrate2_dryt.f  calibrate2_svr.f calibrate2.f \
	meanwind.f meanomeg.f fosberg.f get_fosberg.f get_wind_rh_joint_prob.f pmatch_mean.f \
        pmatch_mean_loc.f lpm.f \
	snow4href.f get_time_packGB2.f

OBJ=${SRCS:.f=.o}

href_ensprod: $(OBJ)
	$(FC) $(FFLAGS) -o href_ensprod sref_ens_gen.f $(OBJ)  $(LIBS)
clean:
	/bin/rm  *.o href_ensprod
copy:
	/bin/cp href_ensprod ../../exec/
#

