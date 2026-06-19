project('fort', 'fortran',
  version: '0.1.0',
  default_options: [
    'buildtype=release',   
    'optimization=3',      
    'debug=false',         
    'b_ndebug=true']
)

# get the python installation
py = import('python').find_installation(pure: false)

### Select Compiler
fc = meson.get_compiler('fortran')

### Robust Linear Algebra Dependency Search
mkl_dep = []
if host_machine.system() == 'darwin'
  # macOS: Use Apple's native Accelerate framework
  mkl_dep = declare_dependency(link_args: ['-framework', 'Accelerate'])
  message('System is macOS: Using Apple Accelerate framework')
else
  # Linux/Windows logic
  if fc.get_id() == 'intel-llvm' or fc.get_id() == 'intel'
    mkl_dep = declare_dependency(link_args: ['-qmkl'])
  else
    mkl_dep = dependency('mkl-rt', required: false)
    if not mkl_dep.found()
      mkl_dep = dependency('openblas', required: false)
    endif
    if not mkl_dep.found()
      mkl_dep = fc.find_library('mkl_rt', required: false)
    endif
  endif
endif

if not mkl_dep.found() and host_machine.system() != 'darwin'
  error('Linear algebra library (MKL or OpenBLAS) not found. Please install one.')
endif

### Compiler Flags
if fc.get_id() == 'intel-llvm' or fc.get_id() == 'intel'
  add_project_arguments('-funroll-loops', '-heap-arrays', language: 'fortran')
elif fc.get_id() == 'gcc'
  add_project_arguments('-funroll-loops', '-fno-stack-arrays', '-ffree-line-length-none', '-cpp', language: 'fortran')
endif

### Fortran File locations
fortran_sources = [
  'src/fort/complex_sparse_linalg.F90',
  'src/fort/utils.F90',
  'src/fort/shared_data.F90',
  'src/fort/input_output.F90',
  'src/fort/gradient.F90',
  'src/fort/integrator.F90',
  'src/fort/main.F90'
]

# Binaries location inside the installed package
bin_dir = py.get_install_dir() / 'Feom' / 'src' / 'fort' / 'executables'

### Build all targets
executable('propagation_SIA', fortran_sources, fortran_args: ['-DSIA'], dependencies: [mkl_dep], install: true, install_dir: bin_dir)
executable('propagation', fortran_sources, fortran_args: [], dependencies: [mkl_dep], install: true, install_dir: bin_dir)
executable('propagation_Print_ADOs', fortran_sources, fortran_args: ['-DPrint_ADOs'], dependencies: [mkl_dep], install: true, install_dir: bin_dir)
executable('propagation_Print_ADOs_SIA', fortran_sources, fortran_args: ['-DSIA', '-DPrint_ADOs'], dependencies: [mkl_dep], install: true, install_dir: bin_dir)

### Install Python Package Structure (Refactored)

# 1. Install specific top-level files
py.install_sources('__init__.py', subdir: 'Feom')
py.install_sources('src/__init__.py', 'src/utils.py', subdir: 'Feom/src')
py.install_sources('src/fort/__init__.py', subdir: 'Feom/src/fort')

# 2. Recursively install Python-only directories
base_install_dir = py.get_install_dir() / 'Feom' / 'src'

install_subdir('src/initialize', install_dir: base_install_dir)
install_subdir('src/hierarchy', install_dir: base_install_dir)