#!/bin/bash
make clean
make
#./ana_delphes fireball_bp_1TeV 0
#./ana_delphes fireball_bp_2TeV 0
./ana_delphes pptt 1
#./ana_delphes ppvv 1
#./ana_delphes ppvtt 1
#./ana_delphes ppvvv 1
#./ana_delphes ppvvtt 0
#./ana_delphes ppvvvv 0
#./ana_delphes pptttt 0
#./ana_delphes ppvvvtt 0
##mv INFOLOG skimmed
#cat INFOLOG | grep Cross-section | tail -n1
#if [ $? == 0 ]; then echo "[INFO] Completed!"; fi
echo "[INFO] Completed!"

