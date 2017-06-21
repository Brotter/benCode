#!/bin/bash

################################################
#
#  Ben Rotter - September 2016
#
#  This code should run on the four ANITA clusters
#  at UH to run whatever code faster
#
#################################################

startTime=`date +%m.%d.%y_%Hh`
sharedDir="/home/brotter/nfsShared/results/templateSearch/"${startTime}
mkdir ${sharedDir}
mkdir ${sharedDir}"/log"

if [ `hostname | cut -d"." -f1` == "anitaI" ]; then
    startCore=0
elif [ `hostname | cut -d"." -f1` == "anitaII" ]; then
    startCore=64
elif [ `hostname | cut -d"." -f1` == "anitaIII" ]; then
    startCore=128
elif [ `hostname | cut -d"." -f1` == "anitaIV" ]; then
    startCore=192
else
    echo "The server isn't an anita cluster server, so you shouldn't use this script"
    exit
fi

numEntries=307151

for localCore in `seq 0 63`; do
    absoluteCore=$((localCore+startCore))
    startEntry=$((numEntries*(absoluteCore)))
    stopEntry=$((numEntries*(abosluteCore+1)))
    nice -n 10 ./templateSearch ${sharedDir}/${absoluteCore} ${startEntry} ${stopEntry} 1> ${sharedDir}/log/${absoluteCore}.log 2>&1 &
done
    