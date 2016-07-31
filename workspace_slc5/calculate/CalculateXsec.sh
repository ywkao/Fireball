#!/bin/bash
RAID2="/raid2/w/ykao/simulation"
MG5="/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation"

Calculate(){
    NUM=`nl -n ln -w 1 log_xsec | tail -n 1 | awk '{print $1}'`
    g++ -o CalculateXsec CalculateXsec.cc
    ./CalculateXsec ${NUM} $1 # need to input # of data
}

grep results ${RAID2}/pptt_[45]?/crossx.html | grep -v " 18" | grep -v ">18" | awk '{printf("%f %f\n",$5,$8)}' > log_xsec
grep results ${MG5}/pptt_52/crossx.html | grep -v " 18" | grep -v ">18" | awk '{printf("%f %f\n",$4,$7)}' >> log_xsec
Calculate "pptt"

grep results /raid2/w/ykao/simulation/ppvv_11/crossx.html | grep -v "0.17" | awk '{printf("%f %f\n",$4,$7)}' > log_xsec 
grep results /raid2/w/ykao/simulation/ppvv_12/crossx.html | grep -v "0.17" | awk '{printf("%f %f\n",$4,$7)}' >> log_xsec 
Calculate "ppvv"

grep results /raid1/w/ykao/simulation/ppvtt/crossx.html | grep "> 1" | grep "100k" | awk '{printf("%f %f\n",$4,$7)}' > log_xsec
grep results /raid1/w/ykao/simulation/ppvtt_01/crossx.html | grep "> 1" | awk '{printf("%f %f\n",$4,$7)}' >> log_xsec
Calculate "ppvtt"

grep results /raid1/w/ykao/simulation/ppvvv/crossx.html | grep "> 0.3" | awk '{printf("%f %f\n",$4,$7)}' > log_xsec
grep results /raid1/w/ykao/simulation/ppvvv_01/crossx.html | grep "> 0.3" | awk '{printf("%f %f\n",$4,$7)}' >> log_xsec
Calculate "ppvvv"

grep result /raid1/w/ykao/simulation/ppvvtt_01/crossx.html | awk '{printf("%f %f\n",$4,$7)}' > log_xsec
Calculate "ppvvtt"
grep result /raid1/w/ykao/simulation/pptttt_01/crossx.html | awk '{printf("%f %f\n",$4,$7)}' > log_xsec
Calculate "pptttt"
grep result /afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/simulation/ppvvvv_01/crossx.html | awk '{printf("%f %f\n",$4,$7)}' > log_xsec
Calculate "ppvvvv"

#cat log_xsec;
