#!/bin/bash
if [[ $1 != "" ]]; then
    for file in signalStrength_*; do
        cat signalStrength.h | head -n10 > log;
        cat $file >> log; cat log > $file;
    done
    
    for N in `ls | grep signalStrength_ | awk 'BEGIN{FS="_"}{printf("%d ",$2)}'`; do 
        echo "int Nobs = $N;" >> signalStrength_$N.h;
    done
fi

NAME=UpperLimitCrossSection
for N in `ls | grep signalStrength_ | awk 'BEGIN{FS="_"}{printf("%d ",$2)}'`; do 
    file=signalStrength_${N}.h
    sed -i "s/signalStrength.h/${file}/g" ExclusionXsec.C;
    root -l -b -q ExclusionXsec.C+; wait;
    mv ${NAME}.png ${NAME}_${N}.png;
    sed -i "s/${file}/signalStrength.h/g" ExclusionXsec.C;
    rm ExclusionXsec_C.so;
done

