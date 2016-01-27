#!/bin/bash
rm -rf ?TeV
root -l -b -q ana_draw.C+
cp -ar output 1TeV
sed -i "`grep -n "SIG=1;" ana_draw.C | awk 'FS=":"{printf"%d\n",$1}'`c 	SIG=2; XERR=1;" ana_draw.C
grep "SIG=2;" ana_draw.C

root -l -b -q ana_draw.C+
cp -ar output 2TeV
sed -i "`grep -n "SIG=2;" ana_draw.C | awk 'FS=":"{printf("%d\n",$1)}'`c 	SIG=1; XERR=1;" ana_draw.C

echo "Done!!!"
