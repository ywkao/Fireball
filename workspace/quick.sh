#!/bin/bash

outputPath='/afs/cern.ch/user/y/ykao/work/fireball/03output'
analysisPath='/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace'
METHOD='NoGA'
SUFFIX='NoGA'

make
if [ $?==0 ];then #why this does NOT work when failing to make??
	cd ${outputPath}
	./reset.sh ${SUFFIX}
	cd ${analysisPath}
	
	for ENE in 1 #2
	do
		./ana_delphes ${METHOD} ${ENE} ${SUFFIX} >> LOG_${ENE}TeV
	done
	
	cd ${outputPath}
	./prepare.sh ${SUFFIX}
fi
#vi LOG_1TeV
#vi LOG_2TeV
