SHELL=/bin/sh
#
SRCS=	WRFBUCKETG2.f

OBJS=	WRFBUCKETG2.o 

# Tunable parameters
#
# LIBS		List of libraries
# CMD		Name of the executable
#
CMD = ./href_fv3_3hqpf

LIBS =   ${G2TMPL_LIB} ${G2_LIB4} ${JASPER_LIB} ${PNG_LIB} ${Z_LIB} ${W3NCO_LIB4} ${W3EMC_LIB4} ${BACIO_LIB4} 
# INC = ${G2_INC4} {G2TMPL_INC}
INC = ${G2_INC4}

FFLAGS = -O   

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.
#
all:	 $(CMD) delete

$(CMD):		$(OBJS)
	$(FC) -I${INC} -o $(@) $(OBJS) $(LIBS)

$(OBJS):       
	$(FC) -I${INC} -c $(FFLAGS) $(@:.o=.f)
clean:
	@rm  *.o
delete:
	@rm  *.o
