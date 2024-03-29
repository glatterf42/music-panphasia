##############################################################################
### compile time configuration options
FFTW3		= yes
MULTITHREADFFTW	= yes
SINGLEPRECISION	= no
HAVEHDF5        = yes
HAVEPANPHASIA	= yes
PANPHASIA_HOME  = ./panphasia
HAVEBOXLIB	= no
BOXLIB_HOME     = ${HOME}/nyx_tot_sterben/BoxLib

##############################################################################
### compiler and path settings
CC      = g++ ## icpc
FF	= gfortran-11# needed only for panphasia, or ## ifort
LINKER	= g++-11 ## ifort # need to link with ifort if using intel
OPT     = -Wall -Wno-unknown-pragmas -O3 -mtune=native #-fsanitize=address -fno-omit-frame-pointer
CFLAGS  =  
LFLAGS  = -lgsl -lgslcblas -ltirpc #-fsanitize=thread  #-fsanitize=address -fno-omit-frame-pointer
FFLAGS  = -ffixed-line-length-132 -O3 -fimplicit-none -g #-fsanitize=address -fno-omit-frame-pointer## use for gfortran
#FFLAGS = -extend_source -O3 -fimplicit-none -g ## use for ifort
CPATHS  = -I. -I$(HOME)/local/include -I/opt/local/include -I/usr/local/include -I/usr/include/hdf5/serial -I/usr/include/tirpc
LPATHS  = -L$(HOME)/local/lib -L/opt/local/lib -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial 


##############################################################################
# if you have FFTW 2.1.5 or 3.x with multi-thread support, you can enable the 
# option MULTITHREADFFTW
ifeq ($(strip $(MULTITHREADFFTW)), yes)
  ifeq ($(CC), icpc)
    CFLAGS += -openmp
    LFLAGS += -openmp
  else
    CFLAGS += -fopenmp
    LFLAGS += -fopenmp
  endif
  ifeq ($(strip $(FFTW3)),yes)
	ifeq ($(strip $(SINGLEPRECISION)), yes)
		LFLAGS  +=  -lfftw3f_threads
	else
		LFLAGS  +=  -lfftw3_threads
	endif
  else
    ifeq ($(strip $(SINGLEPRECISION)), yes)
      LFLAGS  += -lsrfftw_threads -lsfftw_threads
    else
      LFLAGS  += -ldrfftw_threads -ldfftw_threads
    endif
  endif
else
  CFLAGS  += -DSINGLETHREAD_FFTW
endif

ifeq ($(strip $(FFTW3)),yes)
  CFLAGS += -DFFTW3
endif

##############################################################################
# this section makes sure that the correct FFTW libraries are linked
ifeq ($(strip $(SINGLEPRECISION)), yes)
  CFLAGS  += -DSINGLE_PRECISION
  ifeq ($(FFTW3),yes)
    LFLAGS += -lfftw3f
  else
    LFLAGS  += -lsrfftw -lsfftw
  endif
else
  ifeq ($(strip $(FFTW3)),yes)
    LFLAGS += -lfftw3
  else
    LFLAGS  += -ldrfftw -ldfftw
  endif
endif

##############################################################################
#if you have HDF5 installed, you can also enable the following options
ifeq ($(strip $(HAVEHDF5)), yes)
  OPT += -DH5_USE_16_API -DHAVE_HDF5
  LFLAGS += -lhdf5
endif

##############################################################################
CFLAGS += $(OPT)
TARGET  = MUSIC
OBJS    = output.o transfer_function.o Numerics.o defaults.o constraints.o random.o\
		convolution_kernel.o region_generator.o densities.o cosmology.o poisson.o\
		densities.o cosmology.o poisson.o log.o main.o \
		$(patsubst plugins/%.cc,plugins/%.o,$(wildcard plugins/*.cc))

##############################################################################
# stuff for PANPHASIA
#FFLAGS += $(OPT)
ifeq ($(strip $(HAVEPANPHASIA)), yes)
  CFLAGS += -DHAVE_PANPHASIA
  OBJS += generic_lecuyer.o panphasia_routines.o
  ## if linking with g++ 
  LFLAGS += -lgfortran
  ## if linking with ifort
  #LFLAGS += -nofor-main -lstdc++
endif

##############################################################################
# stuff for BoxLib
BLOBJS = ""
ifeq ($(strip $(HAVEBOXLIB)), yes)
  IN_MUSIC = YES
  TOP = ${PWD}/plugins/nyx_plugin
  CCbla := $(CC)
  include plugins/nyx_plugin/Make.ic
  CC  := $(CCbla)
  CPATHS += $(INCLUDE_LOCATIONS)
  LPATHS += -L$(objEXETempDir)
  BLOBJS = $(foreach obj,$(objForExecs),plugins/boxlib_stuff/$(obj))
#
endif

##############################################################################
all: $(OBJS) $(TARGET) Makefile
#	cd plugins/boxlib_stuff; make

bla:
	echo $(OBJS)

ifeq ($(strip $(HAVEPANPHASIA)), yes)
generic_lecuyer.o: $(PANPHASIA_HOME)/generic_lecuyer.f90 # $(ALLDEP)
	$(FF) $(FFLAGS) -c $< 

panphasia_routines.o: $(PANPHASIA_HOME)/panphasia_routines.f $(PANPHASIA_HOME)/generic_lecuyer.f90 generic_lecuyer.o # $(ALLDEP)
	$(FF) $(FFLAGS) -c $< 
endif

ifeq ($(strip $(HAVEBOXLIB)), yes)
$(TARGET): $(OBJS) plugins/nyx_plugin/*.cpp
	cd plugins/nyx_plugin; make BOXLIB_HOME=$(BOXLIB_HOME) FFTW3=$(FFTW3) SINGLE=$(SINGLEPRECISION)
	$(CC) $(LPATHS) -o $@ $^ $(LFLAGS) $(BLOBJS) -lifcore
else
$(TARGET): $(OBJS)
	$(LINKER) $(LPATHS) -o $@ $^ $(LFLAGS)
endif

%.o: %.cc *.hh Makefile 
	$(CC) $(CFLAGS) $(CPATHS) -c $< -o $@

clean:
	rm -rf $(OBJS)
ifeq ($(strip $(HAVEBOXLIB)), yes)
	oldpath=`pwd`
	cd plugins/nyx_plugin; make realclean BOXLIB_HOME=$(BOXLIB_HOME)
endif
	cd $(oldpath)

