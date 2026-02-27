cd ..
ifx -g -traceback -check bounds -warn all -O0 -o main complex_sparse_linalg.F90 input_output.F90 test_main.F90 
mv main executables/
cd -
