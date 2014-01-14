#!/bin/bash
FILENAMEFIRST="EnergyDistributionT_"
FILENAMEEND=".gnu"
FILENAME=""
for((k=0;k<1;k++))
do
echo "set term gif animate delay 40 enhanced;"
echo "set output \"EnergyNetworkOut.gif\";"
for((i=0;i<500;i=i+100))
do
FILENAME=$FILENAMEFIRST$i$FILENAMEEND
if [ -f $FILENAME ]
then
 	echo load \"$FILENAME\"
fi
done
done | gnuplot -persist
echo "Finished!"
