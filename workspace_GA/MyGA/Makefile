######################################################################
# Project: TMVA - Toolkit for Multivariate Data Analysis             #
# Code   : Example usage for TMVA analyses                           #
###################################################################### 

MAKEFLAGS = --no-print-directory -r -s

INCLUDE = $(shell root-config --cflags)
LIBS    = $(shell root-config --libs) -lTMVA -lMLP -lTreePlayer -lMinuit

BINS = MyGA 

all: $(BINS)

clean:
	rm -f $(BINS)

$(BINS): % : %.cc 
	@echo "Building $@ ... "
	$(CXX) $(CCFLAGS) $< $(INCLUDE) $(LIBS) -g -o $@
	@echo "Done"


