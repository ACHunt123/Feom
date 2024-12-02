# Beginning of a make file

# to use specific makefile: make -f mymakefile

# format is:
# target: prerequisites
# <tab> commands
#
# To continue a line on a next line \
# use backslash
#
# $< means the first prerequisite
# or also "any dependency that is more recent than its corresponding target." 
# $^ means all the prerequisits
# $@ means the target




# Set these to represent your directory structure (I use environment variables)
# ===============================================
# Path to Armadillo installation location
# i.e. the directory inside which the include and lib directories are located
ARMADILLO_LOCATION=${SWDIR}
# Path to the rpsimc directory
PREFIX=${HEOMPYDIR}
# And MKL is specified by the $MKLROOT environment variable


#******************************************************************************
# Pick one
#******************************************************************************
CXX          = icpc -Werror -Wall -fno-common #-fopenmp
CXX          = mpicxx -Werror -Wall -fno-common # -DDEBUG -g
CXX          = g++ --pedantic -Werror -Wall -Wextra -Wshadow -fno-common -Wno-unused-variable -Wno-unused-parameter -fopenmp

# Less strict g++ and mpich based compilation
CXX          = mpicxx -Wall -fno-common # -DDEBUG -g


#******************************************************************************
# MKL Flags !!!-lpthread must come after -liomp5!!!
#******************************************************************************
# Taken from:
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
#******************************************************************************
# g++ compilation flags
#MKL_COMP = -m64 -I${MKLROOT}/include
#******************************************************************************
# Sequential program
#MKL_LINK = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
#******************************************************************************
# Parallel program Intel libiomp5
#MKL_LINK = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#******************************************************************************
# Parallel program libgomp
#MKL_LINK = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
#******************************************************************************
#******************************************************************************
# icpc/icc MKL flags
#******************************************************************************
# Sequential program
#MKL_COMP = -mkl=sequential
# or 
#MKL_COMP = -I${MKLROOT}/include

#MKL_LINK = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
#******************************************************************************
# Parallel program
#MKL_COMP = -mkl=parallel
# or 
#MKL_COMP = -I${MKLROOT}/include

#MKL_LINK = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#******************************************************************************


#******************************************************************************
# Pick one
#******************************************************************************
# For g++
MKL_COMP = -I${MKLROOT}/include -m64
MKL_LINK_ARMADILLO = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -ldl -lm #-lsz -lz 
#*******************************

# For icc/icpc
MKL_COMP = -I${MKLROOT}/include
MKL_LINK_ARMADILLO = -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -ldl -lm #-lsz -lz 
#*******************************

# -lmkl_rt by default uses Intel threading and LP64
MKL_LINK = # No need to separately link MKL
#*******************************

# Compiling without MKL
MKL_COMP =
MKL_LINK_ARMADILLO =  
MKL_LINK = -lopenblas 
#*******************************

#MKL_LINK_ARMADILLO = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#MKL_LINK_ARMADILLO = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

#******************************************************************************
# Armadillo flags
#******************************************************************************
ARMADILLO_COMP = -I$(ARMADILLO_LOCATION)/include
# If using wrapper
#ARMADILLO_LINK = -L/home/ap837/sw/lib -larmadillo
# If not using wrapper and using OpenBLAS and LAPACK
#ARMADILLO_LINK = -L/home/ap837/sw/lib -lblas -llapack
# If not using wrapper and using MKL
ARMADILLO_LINK = -L$(ARMADILLO_LOCATION)/lib $(MKL_LINK_ARMADILLO)
# + missing hdf5 flags (can be found in ~/sw/include/armadillo_bits/config.hpp)

#******************************************************************************
# Compilation
#******************************************************************************
# C++ compiler
# Preprocessing flags
CPPFLAGS     = -I$(INCLUDEDIR) $(ARMADILLO_COMP) $(MKL_COMP)
# C++ compiler flags
CXXFLAGS     = -std=c++17 -O3 #-g3
# Linker flags
LDFLAGS      = # -static
# Loaded libraries
LDLIBS       = $(ARMADILLO_LINK) $(MKL_LINK)
# Target architecture
TARGET_ARCH  = #-march=coffeelake

CXX_COMPILE = $(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH)
CXX_LINK = $(CXX) $(LDFLAGS) $(TARGET_ARCH)
#******************************************************************************
# Linking
#******************************************************************************


# Directory tree
BINDIR       = $(PREFIX)/bin
INCLUDEDIR   = $(PREFIX)/include

# Program names
PROGRAM      = heomc
TESTS        = test_general test_utils test_input test_timer

# Executable names
EXEC         = $(PROGRAM)
EXEC_TESTS   = $(addsuffix .x,$(TESTS))

# Targets are inside the include directory
TARGET       = $(addprefix $(INCLUDEDIR)/,$(EXEC))
TARGET_TESTS = $(addprefix $(INCLUDEDIR)/,$(EXEC_TESTS))

# Module names
MODULES      = general timer input utils systems
# Module object files
OBJECT_LIST  = $(addsuffix .o,$(MODULES))
# Module object files inside the include directory
OBJECTS      = $(addprefix $(INCLUDEDIR)/,$(OBJECT_LIST))

#******************************************************************************
# Special targets
#******************************************************************************
# Compulsory first target i.e. the default one
all: folder_tree $(TARGET_TESTS) $(TARGET) 

# Declares explicitly that all, install and clean are phony targets
# This ensures that they run even if files of these names are present
.PHONY: all install clean dvr clinstall

install: all
	install -d $(BINDIR)
	install $(TARGET) $(BINDIR) 
	install $(TARGET_TESTS) $(BINDIR) 

clean:
	$(RM) $(INCLUDEDIR)/*.o $(TARGET) $(TARGET_TESTS)

folder_tree:
	install -d $(INCLUDEDIR)

clinstall: clean install
#******************************************************************************
# Linking executables
#******************************************************************************
$(TARGET):  $(TARGET).o $(OBJECTS) | $(INCLUDEDIR)
	$(CXX_LINK) $< $(OBJECTS) -o $@ $(LDLIBS)

$(INCLUDEDIR)/%.x: %.o $(OBJECTS) | $(INCLUDEDIR)
	$(CXX_LINK) $< $(OBJECTS) -o $@ $(LDLIBS)

#******************************************************************************
# Generic rule for compiling objects
#******************************************************************************
$(INCLUDEDIR)/%.o: %.cpp | $(INCLUDEDIR)
	$(CXX_COMPILE) -c $< -o $@

#******************************************************************************
