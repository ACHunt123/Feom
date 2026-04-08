project('fort', 'fortran',
  version: '0.1.0',
  default_options: [
    'buildtype=release',    # Forces -O3 and removes debug overhead
    'optimization=3',       # High-level optimization
    'debug=false',          # Remove debug symbols for speed
    'b_ndebug=true'         # Disables Fortran assertions (like assert/check)
  ] 
)
# get the python
py = import('python').find_installation(pure: false)


### Select Compiler and add in the flags
fc = meson.get_compiler('fortran')
mkl_dep = []
# ifx flags
if fc.get_id() == 'intel-llvm' or fc.get_id() == 'intel' 
  add_project_arguments('-g', '-funroll-loops', '-heap-arrays', language: 'fortran')
  mkl_dep = declare_dependency(link_args: ['-qmkl'])
# gfortran flags
elif fc.get_id() == 'gcc' 
  add_project_arguments('-g', '-funroll-loops', '-fno-stack-arrays', '-ffree-line-length-none', '-cpp', language: 'fortran')
  mkl_dep = declare_dependency(link_args: ['-lmkl_rt'])
endif

### File locations
fortran_sources = [
  'src/fort/complex_sparse_linalg.F90',
  'src/fort/utils.F90',
  'src/fort/shared_data.F90',
  'src/fort/input_output.F90',
  'src/fort/gradient.F90',
  'src/fort/integrator.F90',
  'src/fort/main.F90'
]
# tell Meson to put binaries inside the Feom/src/fort/executables folder
bin_dir = py.get_install_dir() / 'Feom' / 'src' / 'fort' / 'executables'


### Build  all useful Targets (Matches Python's FORT_SWITCHES logic)
executable('propagation_SIA',
  fortran_sources,
  fortran_args: ['-DSIA'],
  dependencies: [mkl_dep],
  install: true,
  install_dir: bin_dir
)
executable('propagation',
  fortran_sources,
  fortran_args: [],
  dependencies: [mkl_dep],
  install: true,
  install_dir: bin_dir
)
executable('propagation_Print_ADOs',
  fortran_sources,
  fortran_args: ['-DPrint_ADOs'],
  dependencies: [mkl_dep],
  install: true,
  install_dir: bin_dir
)
executable('propagation_Print_ADOs_SIA',
  fortran_sources,
  fortran_args: ['-DSIA', '-DPrint_ADOs'],
  dependencies: [mkl_dep],
  install: true,
  install_dir: bin_dir
)


### Install Python init file for execution of the executable
py.install_sources(
    'src/fort/__init__.py',
    subdir: 'Feom/src/fort'
)

