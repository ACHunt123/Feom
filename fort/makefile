default: 
	# ifx -g -traceback -check all -warn all -ftrapuv -heap-arrays -debug inline-debug-info -O0 -o propagation propagation.f90 -lfftw3 -liomp5 
	# ifx -g -O2 -heap-arrays -o propagation propagation.f90 -lfftw3 -liomp5 
	f2py -c -m propagator propagation.f90 -lfftw3 -liomp5
	# mv *o executables/
	# mv propagation executables/
	mv propagator* executables/

clean:
	rm  executables/*.o executables/*.a executables/*.so executables/*.mod




