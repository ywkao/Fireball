######################################################################
# Project: TMVA - Toolkit for Multivariate Data Analysis             #
# Code   : Example usage for TMVA analyses                           #
###################################################################### 

MAKEFLAGS = --no-print-directory -r -s

#INCLUDE = $(shell root-config --cflags)
CXXFLAGS = $(shell root-config --cflags) -Wno-write-strings -D_FILE_OFFSET_BITS=64 -DDROP_CGAL -I../ -I../external -I../external/tcl
DELPHES_LIBS = $(shell root-config --libs) -lEG -L../ -lDelphes
LIBS    = $(shell root-config --libs) -lTMVA -lMLP -lTreePlayer -lMinuit

BINS = ana_delphes_bg

all: $(BINS)

clean:
	rm -f $(BINS)

#all:
$(BINS): % : %.cc 
	@echo "Building $@ ... "
	#$(CXX) $(CCFLAGS) $< $(INCLUDE) $(LIBS) -g -o $@
	g++ -c ana_delphes_bg.cc $(CXXFLAGS)
	g++ -o ana_delphes_bg ana_delphes_bg.o $(LDFLAGS) $(DELPHES_LIBS) $(LIBS)
	@echo "Done"


