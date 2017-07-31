#*****************************************************************
#                                                                *
#     Makefile for CGMF         11/17/2015                       *
#                                                                *
#*****************************************************************

LDFLAGS	 = -lm 
CPPFLAGS = -O2 # -g -O4 -Wall -Wextra -p
CC       = mpicc
CPP      = mpic++
RM       = rm
INSTALL  = /usr/bin/install
INSTDIR  = /usr/local/bin

OBJS	= terminate.o asympt.o betamain.o betaprofile.o bstrength.o builtin.o\
	 elepop.o ensdf.o etc.o extwave.o kcksyst.o levden.o\
	 masstable.o gampop.o gamcas.o gdr.o gtrans.o intwave.o omcalc.o\
	 omsetform.o omsetparm.o output.o parpop.o popinit.o ripl2levels.o\
	 setup.o specmain.o BrosaYields.o specmc.o smatrix.o stattrans.o evapinterface.o\
	 Yields.o excinterface.o excipop.o exciparm.o excilambda.o exciton.o\
	 FissionFragments.o FissionEvents.o\
	 cgmf.o

COBJS	= mt19937ar.o

PROG	= mpiCGMF

all:		$(PROG)

$(PROG):	$(OBJS) $(COBJS) 
		$(CPP) $(OBJS) $(COBJS) $(LDFLAGS) -o $(PROG) 

clean:
		$(RM) -f $(OBJS) $(COBJS) $(PROG)

install:
		$(INSTALL) -m 775 $(PROG) $(INSTDIR)/


%.o:    %.cpp
	$(CPP)   $(CPPFLAGS) -c $<

### C routines
%.o: %.c
	$(CC)    $(CPPFLAGS) -c $<

