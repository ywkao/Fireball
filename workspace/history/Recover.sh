#!/bin/bash
mv ana_delphes.cc test.cc
mv ana_delphes.cc.keep ana_delphes.cc
make clean
./Exec.sh
echo "[INFO] Completed!"
