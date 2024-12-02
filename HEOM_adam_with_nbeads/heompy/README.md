HEOM simulation program heompy by Adam Prada
============================================
(HEOM = Hierarchical Equations of Motion)

This package was used to generate the data for the following publications:
- **Adam Prada**, Eszter S. Pós, Stuart C. Althorpe; _J. Chem. Phys. 21_ March 2023; 158 (11): 114106. doi: [10.1063/5.0138250](https://doi.org/10.1063/5.0138250)
- **Adam Přáda**, "Hierarchical equations of motion", First-year PhD report, University of Cambridge, 2019, doi: [10.5281/zenodo.7566519](https://doi.org/10.5281/zenodo.7566519)
- **Adam Přáda**, "Dissipative Matsubara Dynamics", PhD Thesis, University of Cambridge, 2023, doi: [10.17863/CAM.95949](https://doi.org/10.17863/CAM.95949)

This package uses python to prepare HEOM simulations and either Python or C++
to propagate these simulations in time. It is possible to use only the Python
version with complete functionality, but none of the speedup offered by C++.
Only simulations with fixed numbers of ADOs can be done using C++ and dynamical
pruning of small ADOs requires using Python.

*In past, this package also used the Ikeda, Scholes libheom and pyheom. This is
now commented out, but if necessary it could be made to work with not much effort.*

How to use after everything is installed
----------------------------------------

The files from the bin directory can be executed with python3 with core libraries and numpy.

In the input directory, you can find `inputgen.py`, which generates a directory
for each required calculation and places the input and mpi files in them.

The calculation is then run simply as:

    heompy input_file_name calculation_name

If the `heomc` simulation type is used, this only prepares `tmp_*` files. The
propagation in time is done by the C++ program as:

    heomc input_file_name calculation_name

which must be done in the directory containing the `tmp_*` files.


If Ikeda Scholes is to be used:
------------------------------
This package uses libheom and pyheom by Ikeda and Scholes. These should be
installed first and can be found in the supplementary information here:
https://doi.org/10.1063/5.0007327

When installing these, use the flag --user for both all, libheom, pylibheom and pyheom
    python3 setup.py install --user

After that the appropriate environment variables should be set as indicated in
the bashrc file.



Potentials:
-----------
-   0: Harmonic oscillator
-   1: Shifting harmonic oscillator as used by Tanimura
-   2: Morse oscillator
-   3: Polynomial potential
-   harmonic: Harmonic oscillator
-   anharmonic: Slightly anharmonic potential
-   quartic: Quartic potential
-   harmoquartic: Mixture of harmonic and quartic potential
-   champagne: Morse oscillator (in 2D known as champagne bottle)
-   hydoh: Hydrated OH bond
-   dw1: Double well potential from Topaler and Makri's paper
-   dw2: Double well potential from Topaler and Makri's paper

Initial states:
---------------
-   wavepacket: Gaussian wave packet
-   tanimura: Tanimurai's initial state (thermal distribution of an isolated harmonic oscillator)
-   thermal: Thermal distribution for the isolated potential
-   thermal\_aren: Thermal distribution for the renormalised potential
-   loaded: Initial state loaded from a file

Types of simulations:
---------------------
-   `highT_qm`: Quantum HEOM (high temperature limit)
-   `highT_cl`: Classical HEOM
-   `lowT_qm`: Quantum HEOM (including the Matsubara terms)
-   `lowT_cl`: Classical HEOM (including the Matsubara terms = non-physical)
-   `lowT_qm_ee`: Quantum HEOM (including the Matsubara terms) in energy eigenstates representation
-   `pyheom`: Uses the library by Ikeda and Scholes
-   `heomc`: For use with the C++ bit of the code. `The same as low_T_qm_ee`

Profiling:
----------
With cProfile you can also profile existing programs, without making any
separate profiling script. Just run program with profiler

    python -m cProfile -o profile_data.pyprof script_to_profile.py

and open profile data in kcachegrind with pyprof2calltree, whose -k switch
automatically opens data in kcachegrind (kcachegrind and valgrind have to be
installed using apt and pyprof2calltree using pip or conda)

    pyprof2calltree -i profile_data.pyprof -k

Installing the C++ part
========================

### Recommended compiler:
- icpc = Intel C++ compiler

icpc (ICC) 19.0.5.281 20190815

(- g++ also works)

g++ (Ubuntu 7.4.0-1ubuntu1~18.04.1) 7.4.0

### Required libraries:
- standard C and C++ libraries
- Armadillo

Expected directory tree
-----------------------
This can be changed in the Makefile

Armadillo installed in ${SWDIR}

The rpsimc directory path is in ${HEOMPYDIR}

MKL location is specified by the standard environment variable ${MKLROOT}

MKL:
----
version used: mkl/64/2019/0/4

On computers of the Department of Chemistry, University of Cambridge, add into bashrc:

    module load mkl/64/2019/0/4
This also sets the required environment variables. In local installations
you have to source a shell script that is in the MKL folder:

    source /opt/intel/oneapi/setvars.sh >> /dev/null

or for older installations

    source mklvars.sh intel64

Armadillo installation:
-----------------------
version used: armadillo-9.600.6

Before the installation make sure to load the MKL module.
(To make sure the cmake finds it)

Install using cmake:

    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX:PATH=${SWDIR} -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc
    make
    make install

Add into bashrc:

    export LD_LIBRARY_PATH="${SWDIR}/lib:${LD_LIBRARY_PATH}"
    export PKG_CONFIG_PATH="${SWDIR}/lib/pkgconfig:${PKG_CONFIG_PATH}"

In ${SWDIR}/include/armadillo_bits/config.hpp comment out:

    #define ARMA_USE_WRAPPER
    =>
    //#define ARMA_USE_WRAPPER

The bottom of the file should have been configured by cmake to use MKL
and should look approximately like:

    // if Armadillo was installed on this system via CMake and ARMA_USE_WRAPPER is not defined,
    // ARMA_AUX_LIBS lists the libraries required by Armadillo on this system, and
    // ARMA_AUX_INCDIRS lists the include directories required by Armadillo on this system.
    // Do not use these unless you know what you are doing.
    #define ARMA_AUX_LIBS /usr/local/shared/intel/compilers_and_libraries_2019.4.227/linux/mkl/lib/intel64/libmkl_rt.so;/usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.so;/usr/lib/x86_64-linux-gnu/libpthread.so;/usr/lib/x86_64-linux-gnu/libsz.so;/usr/lib/x86_64-linux-gnu/libz.so;/usr/lib/x86_64-linux-gnu/libdl.so;/usr/lib/x86_64-linux-gnu/libm.so
    #define ARMA_AUX_INCDIRS /usr/include/hdf5/serial

(i.e. it should contain the mkl bits)

If it loads BLAS and LAPACK, then cmake did not find your MKL installation.

Compiling and linking combination of MKL and Armadillo:
-------------------------------------------------------
This is most easily done using the provided Makefile. Note that you either have
to define the $HEOMPYDIR and $SWDIR environment variable or appropriately change
the beginning of the Makefile

Running a calculation
---------------------
Run the following command:

    mpirun -n number_of_MPI_processes heomc input_file calculation_name

Profiling using Valgrind
------------------------
    valgrind --tool=callgrind ./rpsimc inputfile calculation name number_of_threads # generates callgrind.out.526 (or a different number
    kcachegrind # opens the out file using GUI

Input files
-----------
A sample input file is in the input directory of the repository together with a
python script that can generate input files. It contains dictionaries with all
the available values for each keyword input variable

Debugging
---------
Compile with debug flags 

    -g3 -debug full

and then

    gdb
    file executable.exe
    run input
    bt
