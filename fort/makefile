# -------------------------------------------------------------------------
# Makefile for compiling the Fortran HEOM propagation code
#
# Usage:
#   make [target] [SWITCHES="..."]
#
# If SWITCHES are not specified, defaults to:
#   SWITCHES = -DLowTCorr
#
# Supported compile-time switches (SWITCHES):
#   -DLowTCorr      Enable low-temperature correction via double commutator
#   -DPrint_ADOs    Print the ADOs to file every N timesteps
#   -DPrune         Prune the ADOs dynamically during propagation
#
# Example:
#   make fast SWITCHES="-DLowTCorr -DPrint_ADOs"
#
# Targets:
#   fast            Optimized build
#   debug           Debug build with bounds and warnings
#   profile         Build with profiling enabled
#   f2py            Build Python module via f2py
#   clean           Remove all compiled executables
# -------------------------------------------------------------------------


# Default switches
SWITCHES?= -DLTCorr=1 -DSIA
### Calculate the name of the executable based on switches [replace -D with _, = with nothing, and sort alphabetically]
NAME := propagation
SWITCH_SUFFIX := $(strip $(shell echo $(SWITCHES) | tr ' ' '\n' | sort | sed 's/-D//g; s/=//g' | paste -sd_ -))
ifeq ($(SWITCH_SUFFIX),)
    EXENAME := $(NAME)
else
    EXENAME := $(NAME)_$(SWITCH_SUFFIX)
endif
# Compiler and flags
FC = ifx             # Intel Fortran Compiler (can switch to gfortran if needed)
F90FLAGS = -g -O3 -funroll-loops  -heap-arrays  ${SWITCHES} 
F90DEBUGFLAGS = -g -traceback -check bounds -warn all -O0 ${SWITCHES}			# -ftrapuv -heap-arrays -Debug inline-Debug-info -O0 ### WAY TOO SLOW
FILES = utils.F90 shared_data.F90 input_output.F90 gradient.F90 propagator.F90 propagation.F90 -qmkl 

# Directories
EXEC_DIR = executables

fast: 
	$(FC) $(F90FLAGS) -o $(EXENAME) $(FILES)  
	mv $(EXENAME) executables/

all:
	$(FC) $(F90FLAGS)  -o  propagation $(FILES)
	mv propagation executables/
	wait
	$(FC) $(F90FLAGS) -DLTCorr=1 -o propagation_LowTCorr $(FILES)
	mv propagation_LowTCorr executables/
	wait
	$(FC) $(F90FLAGS) -DPrint_ADOs -o  propagation_Print_ADOs $(FILES)
	mv propagation_Print_ADOs executables/



clean:
	rm  executables/*

debug:
	$(FC) $(F90DEBUGFLAGS) -o $(EXENAME) $(FILES)
	mv $(EXENAME) executables/

f2py:
	gfortran -Wall -g -fPIC -O2  -c $(FILES)  
	f2py -c -m propagator propagation.F90  
	rm propagation.o
	mv propagator* executables/

profile: 
	ifx -g -O0 -fpp -o $(EXENAME) $(FILES)
	# ifx -g -prof -o $(EXENAME) $(FILES)
	mv $(EXENAME) executables/











