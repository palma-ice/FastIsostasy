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

