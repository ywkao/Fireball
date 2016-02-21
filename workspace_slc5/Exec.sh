#!/bin/bash
DIR="~/work/MG5_aMC_v2_2_3/Delphes/workspace/skimmed"
PRESELECTION=1 # 0 = w/o pre-selection; 1 = w/ pre-selection; 3 = mg5 pre-selection
Execution() {
	echo "##### `date +%F\ %H:%M:%S` #####" | tee -a INFOLOG
	./ana_delphes $1 $2 ${PRESELECTION} | tee -a tmplog
}

sed -i "1c const char *savingPath = \"`pwd`/skimmed\";" savingPath.h
make clean; make
echo "" > tmplog
Execution pptt 1
Execution ppvv 1
#Execution fireball_bp_1TeV 0
#Execution fireball_bp_2TeV 0
#Execution ppvtt 1
#Execution ppvvv 1
#Execution ppvvtt 0
#Execution ppvvvv 0
#Execution pptttt 0
#Execution ppvvvtt 0

##mv INFOLOG skimmed
#cat INFOLOG | grep Cross-section | tail -n1
#if [ $? == 0 ]; then echo "[INFO] Completed!"; fi
echo "===== `date +%F\ %H:%M:%S` =====" | tee -a INFOLOG
echo "[INFO] Completed!"

