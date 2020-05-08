# Platform: linux
# Model:    3 Orbital Hubbard Model for Pnictides
# Method:   Monte Carlo-Mean Field

EXENAME  = main
### ------ Personal PC compilation ------------
CXX     = g++
CPPFLAGS = -std=c++11 
LDFLAGS = -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/mkl/lib/intel64/libmkl_sequential.a /opt/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
#LDFLAGS  += -L/opt/intel/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin -liomp5 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core
### ------ Newton compilation ------------
### MUST USE: module load gcc/4.8.2
#CXX = icpc  ### Or use g++ (both works!)
#CPPFLAGS = -std=c++11
#LDFLAGS = /data/apps/lapack/3.5.0/lib/liblapack.a /data/apps/lapack/3.5.0/lib/libblas.a -lgfortran
# -Wl,--start-group /opt/intel/mkl/lib/intel64/lib/intel64/libmkl_intel_ilp64.a /opt/intel/mkl/lib/intel64/lib/intel64/libmkl_sequential.a /opt/intel/mkl/lib/intel64/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
#(DO NOT USE) LDFLAGS  = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lgfortran  (Intel libraries are not found - need fixing for threaded plasma lib)

## --- turn on for production -----------
#CPPFLAGS += -Isrc -I$(HOME)/Codes/0.Codes/ -I/usr/lib/include/      #### Look inside the src folder for header files
#CPPFLAGS += -Wall -w -fopenmp ##-Werror -Wextra #### Enable warnings and treat warnings as errors
#CPPFLAGS += -DNDEBUG #### This disables debugging
#CPPFLAGS += -O3 #### Optimization level here
#STRIP_COMMAND = strip #### "strips off" all the debugging lines of executable

# --- turn on for debugging -----------
CPPFLAGS += -Isrc -I$(HOME)/Codes/0.Codes/ -I/usr/lib/include/
CPPFLAGS += -Wall -Wextra -Wno-sign-compare -w -fopenmp #-Werror #### This enables warnings with extra debugging
CPPFLAGS += -g3 #### link gdb to file system of program
CPPFLAGS += -O0 #### Reduce compilation time and make debugging produce the expected results.
STRIP_COMMAND = true #### Keeps lines in the executable for debugging


$(EXENAME): clean main.o 
	$(CXX) $(CPPFLAGS) -o $(EXENAME)  main.o $(LDFLAGS) 
	$(STRIP_COMMAND) $(EXENAME)

all: $(EXENAME)
	 
clean:
	rm -f $(EXENAME) *.o

######## End of Makefile ########

