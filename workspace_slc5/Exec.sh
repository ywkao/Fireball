#!/bin/bash
PRESELECTION=3 # 0 = w/o pre-selection; 1 = w/ pre-selection; 3 = mg5 pre-selection
Execution() {
	echo "##### `date +%F\ %H:%M:%S` #####" | tee -a INFOLOG
	./ana_delphes $1 $2 ${PRESELECTION} | tee -a tmplog
}

sed -i "1c const char *savingPath = \"`pwd`/skimmed\";" savingPath.h
cat savingPath.h

if [ ! -d skimmed ]; then mkdir skimmed; fi
make clean; make
echo "" > tmplog
Execution fireball_1TeV 0
Execution fireball_2TeV 0
Execution fireball_1.1TeV 0
Execution fireball_1.2TeV 0
Execution fireball_1.3TeV 0
Execution fireball_1.4TeV 0
Execution fireball_1.5TeV 0
Execution fireball_1.6TeV 0
Execution fireball_1.7TeV 0
Execution fireball_1.8TeV 0
Execution fireball_1.9TeV 0
#Execution pptt 1
#Execution ppvv 1
#Execution fireball_bp_1TeV 0
#Execution fireball_bp_2TeV 0
#Execution ppvtt 1
#Execution ppvvv 1
#Execution ppvvtt 0
#Execution ppvvvv 0
#Execution pptttt 0
#Execution ppvvvtt 0

echo "===== `date +%F\ %H:%M:%S` =====" | tee -a INFOLOG
echo "[INFO] Completed!"
top

