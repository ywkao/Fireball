#!/bin/bash
STATUS=-1
TARGET="unweighted" #unweighted/pythia/delphes
if [ ${TARGET} == "unweighted" ]; then FILE="unweighted_events.root";
elif [ ${TARGET} == "pythia" ]; then FILE="tag_1_pythia_lhe_events.root";
elif [ ${TARGET} == "delphes" ]; then FILE="tag_1_delphes_events.root";
else echo "Something wrong"; exit 1; fi

Execution(){
	if [ ${TARGET} == "unweighted" ]; then
		LINE=`grep -n "STATUS=" ana_${TARGET}.C| awk 'FS=":" {printf("%d\n",$1)}'`
		sed -i "${LINE}c const int STATUS=${STATUS};" ana_${TARGET}.C

		LINE=`grep -n "ihtmin=" ana_${TARGET}.h| awk 'FS=":" {printf("%d\n",$1)}'`
		sed -i "${LINE}c const int ihtmin=$1;" ana_${TARGET}.h

		LINE=`grep -n "INPUT=" ana_${TARGET}.h| awk 'FS=":" {printf("%d\n",$1)}'`
		if [ $1 == 2500 ]; then sed -i "${LINE}c const char* INPUT=Form(\"../Events/run_ihtmin_40k_%d/${FILE}\",ihtmin);" ana_${TARGET}.h;
		else sed -i "${LINE}c const char* INPUT=Form(\"../Events/run_ihtmin_20k_%d/${FILE}\",ihtmin);" ana_${TARGET}.h; fi
	else
		for NUM in 500 1000 1500 2000 2500; do
			LINE=`grep -n "STATUS=" ana_${TARGET}_${NUM}.C| awk 'FS=":" {printf("%d\n",$1)}'`
			sed -i "${LINE}c const int STATUS=${STATUS};" ana_${TARGET}_${NUM}.C
		done
	fi

	root -l -b -q Ana.C\(\"${TARGET}\"\)
}

#Execution 0
Execution 500
Execution 1000
Execution 1500
Execution 2000
Execution 2500
ls -lhrt skimmed_iht
./ExecDraw.sh ${TARGET} ${STATUS}

STATUS="m1"
if [ ! -d skimmed_${TARGET}_status_${STATUS} ]; then mkdir skimmed_${TARGET}_status_${STATUS}; fi
cp -a skimmed_iht/*${TARGET}* skimmed_${TARGET}_status_${STATUS} 

