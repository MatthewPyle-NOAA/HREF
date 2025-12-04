SHELL=/bin/sh
#
SRCS=	WRFBUCKETG2.f

OBJS=	WRFBUCKETG2.o 

# Tunable parameters
#
# LIBS		List of libraries
# CMD		Name of the executable
#
CMD = ./enspost_3hqpf

INC = ${G2_INC4}
LIBS =  ${IP_LIB4} ${BACIO_LIB4} ${SP_LIB4} ${G2_LIB4} ${W3EMC_LIB4}  \
	                ${JASPER_LIB} ${PNG_LIB} ${Z_LIB} -ljpeg

FFLAGS_DEBUG =        -g -traceback -O0 -ftrapuv  -check all  -fp-stack-check  -fstack-protector  -assume noold_ldout_format -g -convert big_endian -I $(INC)
FFLAGS =        -fp-model strict -g -traceback -O3 -assume noold_ldout_format -g -convert big_endian -I $(INC)


# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
#
all:	 clean $(CMD)

$(CMD):		$(OBJS)
	$(FC) -I${INC} -o $(@) $(OBJS) $(LIBS)

$(OBJS):       
	$(FC) -I${INC} -c $(FFLAGS) $(@:.o=.f)

# debug:	$(OBJS_DEBUG)
#	$(FC) -I${INC} -o enspost_3hqpf $(OBJS_DEBUG) $(LIBS)

# $(OBJS_DEBUG):
#	$(FC) -I${INC} -c $(FFLAGS_DEBUG) $(@:.o=.f)


clean:
	@rm *.o ${CMD}
delete:
	@rm *.o ${CMD} 
copy:
	/bin/cp ${CMD}  ../../exec/
