#! /bin/sh
#
# M. Pyle 20190502
#
# This makefile adopted from href_ensprod makefile:
#
#   Author: Binbin Zhou
#          Aug 1, 2005
#  OS: Linux (on WCOSS)
#
#

INC=  ${G2_INC4}
LIBS=   ${IP_LIB4} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4} ${SP_LIB4} ${G2_LIB4}  \
        ${JASPER_LIB} ${PNG_LIB} ${Z_LIB}
FFLAGS = -O3 -fp-model precise  -g -traceback -convert big_endian -assume noold_ldout_format -I $(INC) #  -ftrapuv  -check all  -fp-stack-check

OBJS = stitch.o read_grib1.o

.SUFFIXES: .F .f .o .f90

href_ffg_gen: $(OBJS)
	 $(FC) $(FFLAGS) -o href_ffg_gen ffg_gen.f90 $(OBJS)  $(LIBS)

.f90.o:
	$(FC) $(FFLAGS) -c $<

clean:
	/bin/rm  *.o href_ffg_gen
copy:
	/bin/cp href_ffg_gen  ../../exec/

