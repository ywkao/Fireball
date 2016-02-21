#!/bin/bash
TREE="Delphes" #LHEF/LHEF/Delphes
TAGNAME="gen" #unweighted/pythia/gen
FILENAME="tag_1_delphes_events.root" #unweighted_events.root/tag_1_pythia_lhe_events.root/tag_1_delphes_events.root
SOURCE="/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/pptt_Go02/Events"
TARGET="/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace/skimmed_iht"

make
if [ ! -d ${TARGET} ]; then mkdir ${TARGET}; fi
for NUM in 500 1000 1500 2000 2500; do
	./ana_check ${TREE} ${TAGNAME} "${SOURCE}/run_ht_min_10k_${NUM}/${FILENAME}" "${TARGET}/result_pptt_${TAGNAME}_ht_min_${NUM}.root"
done

echo "[INFO] Completed!"
