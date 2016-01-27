#!/bin/bash
TAR="/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace"
SOURCE="/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/pptt_Go02/Events"
#mv ${TAR}/skimmed/result_pptt.root ${TAR}/skimmed/result_pptt_keep.root
for NUM in 500 1000 1500 2000 2500; do
	REPLACE="108c chain->Add(\"${SOURCE}/run_ht_min_${NUM}/tag_1_delphes_events.root\");"
	sed -i "${REPLACE}" ana_delphes.h
	./Exec.sh
	mv ${TAR}/skimmed/result_pptt.root ${TAR}/skimmed/result_pptt_check_ht_min_${NUM}.root
done
#mv ${TAR}/skimmed/result_pptt_keep.root ${TAR}/skimmed/result_pptt.root


#mv ${TAR}/skimmed/result_pptt.root ${TAR}/skimmed/result_pptt_keep.root
#./Exec.sh
#mv ${TAR}/skimmed/result_pptt.root ${TAR}/skimmed/result_pptt_check_21.root
#sed -i '108c         chain->Add("/raid2/w/ykao/simulation/pptt_31/Events/run_100k_31_0/tag_1_delphes_events.root");' ana_delphes.h
#if [ $? == 0 ]; then 
#	./Exec.sh
#	mv ${TAR}/skimmed/result_pptt.root ${TAR}/skimmed/result_pptt_check_31.root
#fi
#mv ${TAR}/skimmed/result_pptt_keep.root ${TAR}/skimmed/result_pptt.root

#./Exec.sh
#mv ana_delphes ana_delphes_bin200
#cp -ar ${TAR}/skimmed ${TAR}/skimmed_bin200

#mv ana_delphes_LPT250 ana_delphes
#./Exec.sh
#mv ana_delphes ana_delphes_LPT250
#cp -ar ${TAR}/skimmed ${TAR}/skimmed_LPT250
#
#mv ana_delphes_ST2400 ana_delphes
#./Exec.sh
#mv ana_delphes ana_delphes_ST2400
#cp -ar ${TAR}/skimmed ${TAR}/skimmed_ST2400

echo "[INFO] Completed!"

