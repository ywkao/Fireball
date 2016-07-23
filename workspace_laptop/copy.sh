#!/bin/bash
DIR=~/Desktop/backup/workspace_`date +%y%m%d`
if [ ! -d ${DIR} ]; then mkdir ${DIR}; fi
cp -aru *.C *.h *.sh multiplicity significance systematics upperlimit README ${DIR}
echo "`rl ${DIR}`"
ls -lhrt ~/Desktop/backup
ls -lhrt ${DIR}
echo "Finish!"
