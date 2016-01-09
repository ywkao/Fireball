#!/bin/bash
mv ana_delphes.cc ana_delphes.cc.keep
mv test.cc ana_delphes.cc
make clean
make
./ana_delphes
echo "[INFO] Completed!"
