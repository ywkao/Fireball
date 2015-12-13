#!/bin/bash
make
./ana_delphes_bg 1 "" >ToBeProcessed/result.log
cat ToBeProcessed/result.log | tail -n2
./Stock.sh TEST
