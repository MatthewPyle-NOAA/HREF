#!/bin/sh
#
#  This makefile is to compile ensemble production generator get precip code
#  Author: Binbin Zhou
#  OS: Intel WCOSS
#  2022: Jun. Du -- added enspost_get_temp.f for fog product
#

INC= ${G2_INC4}
LIBS=   ${IP_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4} ${SP_LIB4} ${G2_LIB4}  \
        ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}
FFLAGS_DEBUG =        -g -traceback -O3 -ftrapuv  -check all  -fp-stack-check  -fstack-protector  -assume noold_ldout_format -g -convert big_endian -I $(INC)
FFLAGS =        -fp-model strict -g -traceback -O3 -assume noold_ldout_format -g -convert big_endian -I $(INC)

all:	clean enspost_get_prcip install

enspost_get_prcip:
	$(FC) $(FFLAGS) -o enspost_get_prcip enspost_get_prcip.f enspost_get_temp.f $(LIBS)
debug:
	$(FC) $(FFLAGS_DEBUG) -o enspost_get_prcip enspost_get_prcip.f enspost_get_temp.f $(LIBS)
clean:
	/bin/rm   enspost_get_prcip
install:
	/bin/cp enspost_get_prcip ../../exec/enspost_get_prcip
#
