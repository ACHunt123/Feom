# Compiler and flags
FC = ifx             # Intel Fortran Compiler (can switch to gfortran if needed)
F90FLAGS = -g -O2 -heap-arrays
F90DEBUGFLAGS = -g -traceback -check all -warn all -ftrapuv -heap-arrays -debug inline-debug-info -O0

# Directories
EXEC_DIR = executables

fast: 
	$(FC) $(F90FLAGS) -o propagation propagation.f90 
	mv propagation executables/

clean:
	rm  executables/*

debug:
	$(FC) $(F90DEBUGFLAGS) -o propagation propagation.f90
	mv propagation executables/

f2py:
	gfortran -Wall -g -fPIC -O2  -c propagation.f90  
	f2py -c -m propagator propagation.f90  
	rm propagation.o
	mv propagator* executables/

default: fast f2py










