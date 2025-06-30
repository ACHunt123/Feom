# Compiler and flags
FC = ifx             # Intel Fortran Compiler (can switch to gfortran if needed)
# F90FLAGS = -g -O2 -heap-arrays
F90FLAGS = -g -O3 -funroll-loops  -heap-arrays # may need to add -heap-arrays (depends on how big the arrays are)
F90DEBUGFLAGS = -g -traceback -check all -warn all -ftrapuv -heap-arrays -debug inline-debug-info -O0
FILES =  shared_data.f90 input_output.f90 propagator.f90 propagation.f90
NAME = propagation
# Directories
EXEC_DIR = executables

fast: 
	$(FC) $(F90FLAGS) -o $(NAME) $(FILES)  
	mv $(NAME) executables/

clean:
	rm  executables/*

debug:
	$(FC) $(F90DEBUGFLAGS) -o $(NAME) $(FILES)
	mv $(NAME) executables/

f2py:
	gfortran -Wall -g -fPIC -O2  -c $(FILES)  
	f2py -c -m propagator propagation.f90  
	rm propagation.o
	mv propagator* executables/

profile: 
	ifx -g -prof -o $(NAME) $(FILES)
	mv $(NAME) executables/

default: fast f2py










