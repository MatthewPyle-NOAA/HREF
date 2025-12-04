SHELL=/bin/sh
#
SRCS=	snowbucket.f GET_BITS.f

OBJS=	snowbucket.o  GET_BITS.o

INC = ${G2_INC4}
LIBS =  ${IP_LIB4} ${BACIO_LIB4} ${SP_LIB4} ${G2_LIB4} ${W3EMC_LIB4}  \
                       ${JASPER_LIB} ${PNG_LIB} ${Z_LIB} -ljpeg
FFLAGS_DEBUG = -O3 -fp-model precise  -g -traceback -convert big_endian -assume noold_ldout_format -I $(INC) -ftrapuv  -check all  -fp-stack-check -check bounds

FFLAGS =        -fp-model strict -g -traceback -O3 -assume noold_ldout_format -g -convert big_endian -I $(INC)

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
#
all:	 clean enspost_snowbucket

$(OBJS):       
	$(FC) -I ${INC} -c $(FFLAGS) $(@:.o=.f)

enspost_snowbucket:
	$(FC) -I ${INC} $(FFLAGS) -o $(@) snowbucket.f GET_BITS.f $(LIBS)
debug:
	$(FC) -I ${INC} $(FFLAGS_DEBUG) -o enspost_snowbucket snowbucket.f GET_BITS.f $(LIBS)
clean:
	/bin/rm   enspost_snowbucket
copy:
	/bin/cp enspost_snowbucket ../../exec/enspost_snowbucket
