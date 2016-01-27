#!/bin/bash
Dir=`pwd`
TARGET=/afs/cern.ch/user/y/ykao/work/Fireball/workspace
echo ${Dir}
cp -ar ${Dir}/*.cc ${Dir}/*.h ${Dir}/*.C ${Dir}/*.sh ${TARGET}
ls -lhrt ${TARGET}
