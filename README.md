# FEOM: Fortran Hierarchical Equations of Motion

FEOM is a high-performance simulation package for open quantum system dynamics, implementing the Hierarchical Equations of Motion (HEOM). Built for speed and scalability, the core propagation engine is written in Fortran with highly optimized complex sparse matrix linear algebra, while providing a seamless Python frontend for automated environment setup, matrix generation, and execution.

## Features

* **High-Performance Fortran Core**: Optimized for both `gfortran` and Intel `ifx` compilers.
* **Advanced Integrators**: Supports standard 4th-Order Runge-Kutta (RK4) and the highly efficient Short Iterative Arnoldi (SIA) method for massive state spaces.
* **Sparse Matrix Algebra**: Custom memory-efficient Complex Compressed Sparse Row (CSR) matrix-vector multiplication with OpenMP vectorization.
* **Python Integration**: Fully integrated Python wrapper via `meson-python` to handle input generation, executable routing, and output parsing.
* **Compile-Time Switches**: Easily toggle features like Auxiliary Density Operator (ADO) printing, Low Temperature Corrections (`LTCorr`), and SIA propagation.

---

## Installation

FEOM can be installed directly as a Python package (recommended) or compiled manually for standalone Fortran execution.

### Prerequisites
* **Fortran Compiler**: `gfortran` (GNU) or `ifx` (Intel)
* **Python**: Python 3.10+ with `pip`
* **Math Libraries**: Intel MKL (for `ifx`) or BLAS/LAPACK (for `gfortran`)

### Option 1: Python Package Installation (Recommended)
Installing via `pip` utilizes the `meson` build system to automatically compile the Fortran backend and make the Python frontend available in your environment.

```bash
# Clone the repository
git clone [https://github.com/yourusername/Feom.git](https://github.com/yourusername/Feom.git)
cd Feom

# Install the package and compile the backend
pip install .
```
*Note: If installing outside a virtual environment, the executables will be safely routed to your user `site-packages`.*

### Option 2: Manual Fortran Build
If you prefer to run the Fortran executable directly or are developing the core numerical routines, you can compile the project using the provided `Makefile`.

```bash
cd src/fort

# Build the optimized executable
make fast

# Build with debug flags and bounds checking
make debug

# Build for memory profiling (Valgrind)
make valgrind
```

You can also pass custom compiler switches during a manual build:
```bash
make fast SWITCHES="-DSIA -DPrint_ADOs"
```

---

## Usage

### Running via Python
The Python wrapper handles the generation of the necessary `.inp` files (Liouvillian, initial density matrix, parameters etc) and routes the correct executable based on your requested switches. See `examples/spinboson/basic_script.py` for minimal working example.

### Running Manually (Standalone Fortran)
To run the Fortran executable directly, ensure the following input files are present in your working directory:
* `Fortparams.inp` - Simulation parameters (`ns`, `dt`, `nttot`...)
* `FortLiouvillian.inp` - The CSR formatted HEOM Liouvillian matrix
* `Fortrho.inp` - The vectorized initial state (ADOs)

Then, simply execute the binary:
```bash
./executables/propagation_SIA
```

---

## Project Structure (not finished)

```text
Feom/
├── pyproject.toml              # Meson build configuration
├── meson.build                 # Build instructions for Fortran sources
├── src/
│   ├── hierarchy...                 # Python 
│   ├── initialize...                # Python 
│   └── fort/                   # Core Fortran source files
│       ├── Makefile            # Manual compilation script
│       ├── main.F90            # Main propagation loop and I/O handling
│       ├── integrator.F90      # RK4 and Short Iterative Arnoldi logic
│       ├── gradient.F90        # HEOM gradient calculations
│       ├── complex_sparse_linalg.F90 # Custom CSR matrix operations
│       ├── input_output.F90    # File reading/writing routines
│       ├── utils.F90           # Math utilities (inner products, norms, SVD)
│       └── shared_data.F90     # Global state variables
└── examples/                   # Some example scripts
└── benchmark/                  # Some bencharking scripts
└── notes/                  # Some old notes
```

---

## Contributing
Contributions, issues, and feature requests are welcome. When modifying the Fortran core, please ensure you test your changes using the `make debug` and `make valgrind` targets to catch memory leaks or out-of-bounds errors before pushing.


## References 
[1]Liping Chen and Qiang Shi. Quantum rate dynamics for proton transfer reactions in condensed
phase: The exact hierarchical equations of motion approach. The Journal of Chemical Physics,130(13):134505, April 2009.

[2] Jie Hu, Meng Luo, Feng Jiang, Rui-Xue Xu, and YiJing Yan. Pad´e spectrum decompositions of quantum distribution functions and optimal hierarchical equations of motion construction for quantum open systems. The Journal of Chemical Physics, 134(24):244106, June 2011.

[3] Akihito Ishizaki and Yoshitaka Tanimura. Quantum Dynamics of System Strongly Coupled to Low-Temperature Colored Noise Bath: Reduced Hierarchy Equations Approach. Journal of the Physical Society of Japan, 74(12):3131–3134, December 2005.

[4] Yuji Nakatsukasa, Olivier S`ete, and Lloyd N. Trefethen. The AAA Algorithm for Rational Approximation. SIAM Journal on Scientific Computing, 40(3):A1494–A1522, January 2018.
