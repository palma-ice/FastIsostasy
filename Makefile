.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = lib/include
bindir = lib/bin

# Command-line options at make call
debug    ?= 0
parallel ?= 0 

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

FC = mpif90
#FC = gfortran 

INC_NC  = -I/opt/local/include
LIB_NC  = -L/opt/local/lib -lnetcdff -L/opt/local/lib -Wl,-headerpad_max_install_names -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk -arch x86_64 -lnetcdf -lnetcdf

REMBOROOT = rembo1
INC_REMBO = -I${REMBOROOT}/librembo/include 
LIB_REMBO = -L${REMBOROOT}/librembo/include -lrembo

#YELMOROOT = yelmo
YELMOROOT = /Users/robinson/models/EURICE/yelmo_gmd
INC_YELMO = -I${YELMOROOT}/libyelmo/include 
LIB_YELMO = -L${YELMOROOT}/libyelmo/include -lyelmo

LISROOT = /Users/robinson/apps/lis/lis
INC_LIS = -I${LISROOT}/include 
LIB_LIS = -L${LISROOT}/lib -llis

FFLAGS  = -ffree-line-length-none -I$(objdir) -J$(objdir)

ifeq ($(parallel), 1)
    # Overwrite default choices with openmp relevant choices 

    LISROOT = /Users/robinson/apps/lis/lis-omp
    INC_LIS = -I${LISROOT}/include 
    LIB_LIS = -L${LISROOT}/lib/ -llis

    FFLAGS  = -I$(objdir) -J$(objdir) -m64 -ffree-line-length-none -fomit-frame-pointer -fopenmp 

endif 

LFLAGS  = $(LIB_NC) $(LIB_LIS) $(LIB_YELMO)

DFLAGS_NODEBUG = -O2
DFLAGS_DEBUG   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
DFLAGS_PROFILE = -O2 -pg


# Determine whether to use normal flags or debugging flags
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

###############################################
##
## List of source files
##
###############################################

$(objdir)/isostasy.o: $(srcdir)/isostasy.f90 $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ncio.o: $(srcdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/nml.o: $(srcdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

###############################################
##
## Compilation of complete programs
##
###############################################

test_isostasy : $(objdir)/isostasy.o $(objdir)/ncio.o $(objdir)/nml.o
		$(FC) $(DFLAGS) $(FFLAGS) $(INC_LIS) -o $(bindir)/test_isostasy.x tests/test_isostasy.f90 \
			$(LFLAGS) $^
		@echo " "
		@echo "    test_isostasy.x is ready."
		@echo " "

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM

