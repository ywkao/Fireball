#!/bin/bash
Nl=$1
Nj=$2
Pl=$3
Pj=$4
LPT=$5
HT=$6
MET=$7
ST=$8
DIR="/afs/cern.ch/user/y/ykao/work/fireball/01source"
PREFIX="simulation_delphes"
DIR_ppvv="/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/ppvv/Events"
FILE_ppvv="tag_1_delphes_events.root"

if [ $# == 8 ]; then # only one input parameter
	NAME="skimmed_${Nl}_${Nj}_${Pl}_${Pj}_${LPT}_${HT}_${ST}_${MET}"
elif [ $# == 1 ]; then
	NAME=$1
else
	echo "[ERROR] Please specify either dirName or selectionCuts."
	exit 1
fi

function Skimming(){
	make
	if [ $? == 0 ]; then
		mkdir skimmed
		for (( i=0; i<10; i++))
		do
			./ana_skimming ${Nl} ${Nj} ${Pl} ${Pj} ${LPT} ${HT} ${ST} ${MET} ${DIR_ppvv}/run_1M_${i}/${FILE_ppvv}  ppvv_run_1M_${i}.root
		done
		./ana_chain
		#./ana_skimming ${Nl} ${Nj} ${Pl} ${Pj} ${LPT} ${HT} ${ST} ${MET} ${DIR}/${PREFIX}_ppvv_jet_matching.root   ${PREFIX}_ppvv_skimmed.root
		./ana_skimming ${Nl} ${Nj} ${Pl} ${Pj} ${LPT} ${HT} ${ST} ${MET} ${DIR}/${PREFIX}_ppvvv_jet_matching.root  ${PREFIX}_ppvvv_skimmed.root
		./ana_skimming ${Nl} ${Nj} ${Pl} ${Pj} ${LPT} ${HT} ${ST} ${MET} ${DIR}/${PREFIX}_ppvvvv_jet_matching.root ${PREFIX}_ppvvvv_skimmed.root
		./ana_skimming ${Nl} ${Nj} ${Pl} ${Pj} ${LPT} ${HT} ${ST} ${MET} ${DIR}/${PREFIX}_pptt_jet_matching.root   ${PREFIX}_pptt_skimmed.root
		./ana_skimming ${Nl} ${Nj} ${Pl} ${Pj} ${LPT} ${HT} ${ST} ${MET} ${DIR}/${PREFIX}_ppvtt_jet_matching.root  ${PREFIX}_ppttv_skimmed.root
		./ana_skimming ${Nl} ${Nj} ${Pl} ${Pj} ${LPT} ${HT} ${ST} ${MET} ${DIR}/${PREFIX}_ppvvtt_jet_matching.root ${PREFIX}_ppttvv_skimmed.root
		./ana_skimming ${Nl} ${Nj} ${Pl} ${Pj} ${LPT} ${HT} ${ST} ${MET} ${DIR}/${PREFIX}_ppvvvtt.root	  		   ${PREFIX}_ppttvvv_skimmed.root
		./ana_skimming ${Nl} ${Nj} ${Pl} ${Pj} ${LPT} ${HT} ${ST} ${MET} ${DIR}/${PREFIX}_pptttt.root 			   ${PREFIX}_pptttt_skimmed.root
		if [ -d ${NAME} ]; then
			rm -rf ./${NAME}
		fi
		mv skimmed ./${NAME}
		echo "[INFO] Skimming Completed!"
	else
		echo "[ERROR] Compilation failed!"
	fi
}

function Calculation(){
	if [ -d ${NAME} ]; then
		sed -i '1c const char  *PREFIX  = "/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD/'${NAME}'/simulation_delphes";' path.h
		root -l -q ana_draw_ABCD.C+
		#vi -R INFOLOG
		echo ${NAME} >> LOG
		echo "[NumJet vs. NumLep]" >> LOG
		cat INFOLOG | grep "A/B" | tail -n 5 | head -n 1 >> LOG
		cat INFOLOG | grep "D/C" | tail -n 5 | head -n 1 >> LOG
		echo "[HT vs. MET]" >> LOG
		cat INFOLOG | grep "A/B" | tail -n 5 | tail -n 1 >> LOG
		cat INFOLOG | grep "D/C" | tail -n 5 | tail -n 1 >> LOG
		echo " " >> LOG
		mv INFOLOG ./${NAME}
	else
		echo "[ERROR] ${NAME} does not exist."
	fi
}

if [ $# == 1 ]; then # only one input parameter
	Calculation
else
	Skimming
	Calculation
fi

#NAME="skimmed_2_5_30_40_30_40_70_30"
#COMMAND=$1
#if [ $COMMAND == "-p" ]; then
#else
#	echo "Perform ABCD method only"
#fi

#Nl=2
#Nj=6
#Pl=30
#Pj=40
#LPT=30
#HT=70
#MET=30
#ST=130
