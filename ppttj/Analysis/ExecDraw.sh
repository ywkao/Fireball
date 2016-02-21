#!/bin/bash
#TagName=$1 #unweighted/pythia/delphes
#STATUS=$2
Execution() {
	DIR=`pwd`/$3
	if [ ! -d ${DIR} ]; then mkdir ${DIR}; fi
	#root -l $1
	root -l -b -q $1
	mv ${2}.pdf ${DIR}/${2}.pdf
	if [ -f ${2}.png ];  then mv ${2}.png ${DIR}/${2}.png; fi
	if [ -f ${2}.root ]; then mv ${2}.root ${DIR}/${2}.root; fi
	#if [ $? == 0 ]; then root -l ${DIR}/${2}.root; fi
	#if [ $? == 0 ]; then evince ${DIR}/${2}.pdf & fi
}

#Execution Draw_ihtmin_GenHT.C\(\"${TagName}\",${STATUS}\) Comparison_ihtmin_${TagName}_HT_status${STATUS} hist


for TagName in "unweighted" "pythia" "delphes"; do
	for STATUS in -1 1 2 3; do
		if [ ${STATUS} == -1 ]; then Execution Draw_ihtmin_GenHT.C\(\"${TagName}\",${STATUS}\) Comparison_ihtmin_${TagName}_HT_status_m1 hist
		else Execution Draw_ihtmin_GenHT.C\(\"${TagName}\",${STATUS}\) Comparison_ihtmin_${TagName}_HT_status${STATUS} hist
		fi
	done
done
