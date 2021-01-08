#!/bin/sh
#
#  This makefile is to compile ensembel production generator
#  Author: Binbin Zhou
#  OS: Intel WCOSS
#

# FC=     ${COMP} 
INC= ${G2_INC4}
LIBS=   ${IP_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4} ${SP_LIB4} ${G2_LIB4}  \
        ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}
# FFLAGS =        -g -traceback -O3 -ftrapuv  -check all  -fp-stack-check  -fstack-protector  -assume noold_ldout_format -g -convert big_endian -I $(INC)
FFLAGS =        -fp-model strict -g -traceback -O3 -assume noold_ldout_format -g -convert big_endian -I $(INC)

all:	clean href_get_prcip copy

href_get_prcip:
	$(FC) $(FFLAGS) -o $(@) href_get_prcip.f $(LIBS)
clean:
	/bin/rm   href_get_prcip
copy:
	/bin/cp href_get_prcip ../../exec/href_get_prcip
#
