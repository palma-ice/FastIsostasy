# Shared build configuration for FastIsostasy (dependency wiring).
#
# Loaded after the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). References FFLAGS /
# FFLAGS_OPENMP (compiler) and LIB_NC (machine or auto-detected netCDF).

# Dependency paths (serial build by default).
FESMUTILSROOT = fesm-utils
FESMUTILSLIBDIR = ${FESMUTILSROOT}/include-serial
INC_FESMUTILS = -I${FESMUTILSLIBDIR}
LIB_FESMUTILS = -L${FESMUTILSLIBDIR} -lfesmutils

FFTWROOT = fesm-utils/fftw/fftw-serial
INC_FFTW = -I${FFTWROOT}/include
LIB_FFTW = -L${FFTWROOT}/lib -lfftw3 -lm

# OpenMP build (make openmp=1): swap the serial deps for OpenMP variants and
# append the compiler's OpenMP flag (FFLAGS_OPENMP, set in the compiler fragment).
ifeq ($(openmp), 1)
    FESMUTILSLIBDIR = ${FESMUTILSROOT}/include-omp

    FFTWROOT = fesm-utils/fftw/fftw-omp
    INC_FFTW = -I${FFTWROOT}/include
    LIB_FFTW = -L${FFTWROOT}/lib -lfftw3_omp -lfftw3 -lm

    FFLAGS += $(FFLAGS_OPENMP)
endif

LFLAGS_EXTRA ?= -Wl,-zmuldefs
LFLAGS = $(LIB_NC) $(LIB_FESMUTILS) $(LIB_FFTW) $(LFLAGS_EXTRA)
