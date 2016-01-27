#!/bin/bash
#SOURCE="/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/pptt_test/Events"
SOURCE="/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/pptt_Go02/Events"
TARGET="/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace/skimmed"

make
for NUM in 500 1000 1500 2000 2500; do
	./ana_genParticle "${SOURCE}/run_ht_min_${NUM}/tag_1_delphes_events.root" "${TARGET}/result_pptt_gencheck_ht_min_${NUM}.root"
done

echo "[INFO] Completed!"

#./ana_genParticle "${SOURCE}/run_ht_min_2500/tag_2_delphes_events.root" "${TARGET}/result_pptt_gencheck_ht_min_2500.root"
#./ana_genParticle "${SOURCE}/run_ht_min_2000/tag_4_delphes_events.root" "${TARGET}/result_pptt_gencheck_ht_min_2000.root"
#./ana_genParticle "${SOURCE}/run_ht_min_1500/tag_1_delphes_events.root" "${TARGET}/result_pptt_gencheck_ht_min_1500.root"
#./ana_genParticle "${SOURCE}/run_ht_min_1000/tag_1_delphes_events.root" "${TARGET}/result_pptt_gencheck_ht_min_1000.root"
#./ana_genParticle "${SOURCE}/run_ht_min_500/tag_1_delphes_events.root"  "${TARGET}/result_pptt_gencheck_ht_min_500.root"
