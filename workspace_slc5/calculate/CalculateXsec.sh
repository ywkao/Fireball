#!/bin/bash
RAID2="/raid2/w/ykao/simulation"
MG5="/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation"

grep results ${RAID2}/pptt_[45]?/crossx.html | grep -v " 18" | grep -v ">18" | awk '{printf("%5.4f %4.3f\n",$5,$8)}' > log_xsec
grep results ${MG5}/pptt_52/crossx.html | grep -v " 18" | grep -v ">18" | awk '{printf("%5.4f %4.3f\n",$4,$7)}' >> log_xsec

#grep results /raid2/w/ykao/simulation/ppvv_11/crossx.html | grep -v "0.17" | awk '{printf("%5.4f %4.3f\n",$4,$7)}' > log_xsec 
#grep results /raid2/w/ykao/simulation/ppvv_12/crossx.html | grep -v "0.17" | awk '{printf("%5.4f %4.3f\n",$4,$7)}' >> log_xsec 

NUM=`nl -n ln -w 1 log_xsec | tail -n 1 | awk '{print $1}'`

g++ -o CalculateXsec CalculateXsec.cc
./CalculateXsec ${NUM} # need to input # of data
