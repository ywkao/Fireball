#!/bin/bash
Dir=`pwd`
TARGET=/afs/cern.ch/user/y/ykao/work/Fireball/workspace_slc5
echo ${Dir}
cp -ar ${Dir}/*.cc ${Dir}/*.h ${Dir}/*.C ${Dir}/*.sh calculate ${TARGET}
echo "ls -lhrt ${TARGET}"
ls -lhrt ${TARGET}
