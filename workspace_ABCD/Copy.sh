#!/bin/bash
Dir1=/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace
Dir2=/afs/cern.ch/user/y/ykao/work/MG5_aMC_v2_2_3/Delphes/workspace_ABCD
cp -ar ${Dir1} ${Dir2} ~/work/fireball/
cp -ar ${Dir1} ${Dir2} ~/work/Fireball/
ls -lhrt ~/work/Fireball/workspace/
ls -lhrt ~/work/Fireball/workspace_ABCD/
