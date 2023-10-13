.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = libisostasy/include
bindir = libisostasy/bin

# Command-line options at make call
debug    ?= 0
parallel ?= 0 

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

#FC = mpif90
FC = gfortran 

#INC_NC  = -I/opt/local/include
#LIB_NC  = -L/opt/local/lib -lnetcdff -L/opt/local/lib -Wl,-headerpad_max_install_names -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk -arch x86_64 -lnetcdf -lnetcdf

INC_NC  = -I/opt/local/include
LIB_NC  = -L/opt/local/lib -lnetcdff -L/opt/local/lib -Wl,-headerpad_max_install_names -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX11.3.sdk -lnetcdf -lnetcdf


#LISROOT = /Users/robinson/apps/lis/lis
#INC_LIS = -I${LISROOT}/include 
#LIB_LIS = -L${LISROOT}/lib -llis
LISROOT = 
INC_LIS = 
LIB_LIS =

INC_FFTW = -I/usr/local/include
LIB_FFTW = -L/usr/local/lib -lfftw3 -lm

FFLAGS  = -ffree-line-length-none -I$(objdir) -J$(objdir)

ifeq ($(parallel), 1)
    # Overwrite default choices with openmp relevant choices 

    LISROOT = /Users/montoya/apps/lis/lis-omp
    INC_LIS = -I${LISROOT}/include 
    LIB_LIS = -L${LISROOT}/lib/ -llis

    FFLAGS  = -I$(objdir) -J$(objdir) -m64 -ffree-line-length-none -fomit-frame-pointer -fopenmp 

endif 

LFLAGS  = $(LIB_NC) $(LIB_LIS) $(LIB_FFTW)

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


$(objdir)/isostasy_defs.o: $(srcdir)/isostasy_defs.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/isostasy.o: $(srcdir)/isostasy.f90 $(objdir)/isostasy_defs.o $(objdir)/solver_elva.o $(objdir)/solver_lv_elva.o $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solver_elva.o: $(srcdir)/solver_elva.f90 $(objdir)/isostasy_defs.o $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/solver_lv_elva.o: $(srcdir)/solver_lv_elva.f90 $(objdir)/isostasy_defs.o $(objdir)/nml.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/isostasy_benchmarks.o: $(srcdir)/isostasy_benchmarks.f90 $(objdir)/isostasy_defs.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ncio.o: $(srcdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/nml.o: $(srcdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

isostasy_libs = 		$(objdir)/isostasy_defs.o \
						$(objdir)/isostasy.o \
						$(objdir)/solver_elva.o \
						$(objdir)/solver_lv_elva.o \

###############################################
##
## Compilation of complete programs
##
###############################################

# Static library compilation
isostasy-static: $(isostasy_libs)
	ar rc $(objdir)/libisostasy.a $(isostasy_libs)
	ranlib $(objdir)/libisostasy.a
	@echo " "
	@echo "    $(objdir)/libisostasy.a is ready."
	@echo " "

test_isostasy : isostasy-static $(objdir)/isostasy_benchmarks.o $(objdir)/ncio.o $(objdir)/nml.o
		$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_isostasy.x tests/test_isostasy.f90 \
			-L${CURDIR}/libisostasy/include -lisostasy $(LFLAGS) $(objdir)/isostasy_benchmarks.o $(objdir)/ncio.o $(objdir)/nml.o
		@echo " "
		@echo "    test_isostasy.x is ready."
		@echo " "

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM

