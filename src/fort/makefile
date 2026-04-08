# -------------------------------------------------------------------------
# Makefile for compiling the Fortran HEOM propagation code
#
# Usage:
#   make [target] [SWITCHES="..."]
#
# If SWITCHES are not specified, defaults to:
#   SWITCHES = -DSIA
#
# Supported compile-time switches (SWITCHES):
#   -DPrint_ADOs    Print the ADOs to file every N timesteps
#
# Example:
#   make fast SWITCHES=" -DPrint_ADOs"
#
# Targets:
#   fast            Optimized build
#   debug           Debug build with bounds and warnings
#   profile         Build with profiling enabled
#   valgrind         Build for valgrind
#   clean           Remove all compiled executables
# -------------------------------------------------------------------------


### Default switches (if not specified)
SWITCHES?= -DSIA

### Calculate the name of the executable based on switches [replace -D with _, = with nothing, and sort alphabetically]
NAME := propagation
SWITCH_SUFFIX := $(strip $(shell echo $(SWITCHES) | tr ' ' '\n' | sort | sed 's/-D//g; s/=//g' | paste -sd_ -))
ifeq ($(SWITCH_SUFFIX),)
    EXENAME := $(NAME)
else
    EXENAME := $(NAME)_$(SWITCH_SUFFIX)
endif

### Compiler Selection: ('gfortran' or 'ifx')
FC = ifx
# ifx Flags 
ifeq ($(FC), ifx)
    F90FLAGS = -g -O3 -funroll-loops -heap-arrays ${SWITCHES}
    F90DEBUGFLAGS = -g -traceback -check all -warn all -O0 -heap-arrays${SWITCHES}
	F90VALGRINDFLAGS = -g -traceback -O0 -heap-arrays ${SWITCHES} # Valgrind specific (Removes -check all so Valgrind doesn't choke)
	F90PROFILEFLAGS = -g -pg -O3 -heap-arrays ${SWITCHES}
    LIBS = -qmkl 
# gfortra Flags
else ifeq ($(FC), gfortran)
    # -fno-stack-arrays is the equivalent of -heap-arrays
    F90FLAGS = -g -O3 -funroll-loops -fno-stack-arrays -ffree-line-length-none -cpp ${SWITCHES}
    F90DEBUGFLAGS = -g -fbacktrace -fcheck=all -Wall -O0 -ffree-line-length-none -cpp ${SWITCHES}
	F90VALGRINDFLAGS = -g -fbacktrace -O0 -fno-stack-arrays -ffree-line-length-none -cpp ${SWITCHES} # Valgrind specific (Removes -check all so Valgrind doesn't choke)
	F90PROFILEFLAGS = -g -pg -O3 -fno-stack-arrays -ffree-line-length-none -cpp ${SWITCHES}
    LIBS = -L${MKLROOT}/lib/intel64 -lmkl_rt # gfortran usually links MKL via explicit paths or -lblas -llapack
endif

### Files for the compilation
FILES = complex_sparse_linalg.F90 utils.F90 shared_data.F90 \
        input_output.F90 gradient.F90 integrator.F90 main.F90 $(LIBS)

### Where to store the executables
EXEC_DIR = executables
$(shell mkdir -p $(EXEC_DIR))

fast: 
	$(FC) $(F90FLAGS) -o $(EXENAME) $(FILES)  
	mv $(EXENAME) $(EXEC_DIR)/

all:
	$(FC) $(F90FLAGS)  -o  propagation $(FILES)
	mv propagation $(EXEC_DIR)/
	wait
	$(FC) $(F90FLAGS) -DLTCorr=1 -o propagation_LowTCorr $(FILES)
	mv propagation_LowTCorr $(EXEC_DIR)/
	wait
	$(FC) $(F90FLAGS) -DPrint_ADOs -o  propagation_Print_ADOs $(FILES)
	mv propagation_Print_ADOs $(EXEC_DIR)/



clean:
	rm  $(EXEC_DIR)/*

debug:
	$(FC) $(F90DEBUGFLAGS) -o $(EXENAME) $(FILES)
	mv $(EXENAME) $(EXEC_DIR)/

# f2py:
# 	gfortran -Wall -g -fPIC -O2  -c $(FILES)  
# 	f2py -c -m propagator propagation.F90  
# 	rm propagation.o
# 	mv propagator* $(EXEC_DIR)/

valgrind:
	$(FC) $(F90VALGRINDFLAGS) -o $(EXENAME) $(FILES)
	mv $(EXENAME) $(EXEC_DIR)/

profile: 
	$(FC) $(F90PROFILEFLAGS) -o $(EXENAME) $(FILES)
	mv $(EXENAME) $(EXEC_DIR)/











