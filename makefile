# programming environment
COMPILER     := mpic++
INCLUDE      := -L/Users/alessandrocasalino/Software/class_public -I/Users/alessandrocasalino/Software/class_public/include -I/Users/alessandrocasalino/Software/class_public/source -L/usr/local/Cellar/hdf5/1.12.0/lib -I/usr/local/Cellar/hdf5/1.12.0/include -I/Users/alessandrocasalino/Software/LATfield2# add the path to LATfield2 and other libraries (if necessary)
LIB          := -lmpi -lfftw3 -lm -lhdf5 -lgsl -lgslcblas -lclass
HPXCXXLIB    := -lhealpix_cxx -lcfitsio

# target and source
EXEC         := gevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5

# optional compiler settings (LATfield2)
#DLATFIELD2   += -DH5_HAVE_PARALLEL
#DLATFIELD2   += -DEXTERNAL_IO # enables I/O server (use with care)
#DLATFIELD2   += -DSINGLE      # switches to single precision, use LIB -lfftw3f

# optional compiler settings (gevolution)
DGEVOLUTION  := -DPHINONLINEAR
DGEVOLUTION  += -DBENCHMARK
DGEVOLUTION  += -DEXACT_OUTPUT_REDSHIFTS
#DGEVOLUTION  += -DVELOCITY      # enables velocity field utilities
#DGEVOLUTION  += -DCOLORTERMINAL
#DGEVOLUTION  += -DCHECK_B
DGEVOLUTION  += -DHAVE_CLASS    # requires LIB -lclass
#DGEVOLUTION  += -DHAVE_HEALPIX  # requires LIB -lchealpix

# MG module
DGEVOLUTION  += -DMG		# enable MG module
DGEVOLUTION  += -DMGVERBOSE	# enable MG verbose logs

# further compiler options
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)

lccat: lccat.cpp
	$(COMPILER) $< -o $@ $(OPT) $(DGEVOLUTION) $(INCLUDE)

lcmap: lcmap.cpp
	$(COMPILER) $< -o $@ $(OPT) -fopenmp $(DGEVOLUTION) $(INCLUDE) $(LIB) $(HPXCXXLIB)

clean:
	-rm -f $(EXEC) lccat lcmap
