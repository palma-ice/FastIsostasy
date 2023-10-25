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

<COMPILER_CONFIGURATION>

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
			-L${CURDIR}/libisostasy/include -lisostasy $(objdir)/isostasy_benchmarks.o $(objdir)/ncio.o $(objdir)/nml.o $(LFLAGS)
		@echo " "
		@echo "    test_isostasy.x is ready."
		@echo " "

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM

