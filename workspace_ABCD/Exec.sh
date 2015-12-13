#!/bin/bash
make
if [ $? == 0 ]; then
	./ana_skimming /afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvv_jet_matching.root simulation_delphes_ppvv_skimmed.root
	./ana_skimming /afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvv_jet_matching.root simulation_delphes_ppvvv_skimmed.root
	./ana_skimming /afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvvv_jet_matching.root simulation_delphes_ppvvvv_skimmed.root
	./ana_skimming /afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_pptt_jet_matching.root simulation_delphes_pptt_skimmed.root
	./ana_skimming /afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvtt_jet_matching.root simulation_delphes_ppttv_skimmed.root
	./ana_skimming /afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvtt_jet_matching.root simulation_delphes_ppttvv_skimmed.root
	./ana_skimming /afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_ppvvvtt.root simulation_delphes_ppttvvv_skimmed.root
	./ana_skimming /afs/cern.ch/user/y/ykao/work/fireball/01source/simulation_delphes_pptttt.root simulation_delphes_pptttt_skimmed.root
	echo "Completed!"
else
	echo "Compilation failed!"
fi
