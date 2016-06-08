#!/bin/bash
SIG=2
Execution() {
	DIR=`pwd`/$2
	if [ ! -d ${DIR} ]; then mkdir ${DIR}; fi
	root -l -b -q $1
	#root -l $1
}


#Execution ana_yield.C+ output
Execution ana_draw.C+ output

#root -l ana_uncertainty.C+

#Execution ana_draw_v1.C+ output_v1
#Execution ana_draw.C\(${SIG}\)+ output
