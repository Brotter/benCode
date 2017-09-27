#!/bin/bash

################################################
#
#  Ben Rotter - September 2016
#
#  Same as runClusters, but only on one and for wais events
#  
#
#################################################

startTime=`date +%m.%d.%y_%Hh`
sharedDir="/home/brotter/nfsShared/results/templateSearch/"${startTime}"_wais"
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



for localCore in `seq 0 63`; do
    absoluteCore=$((localCore+startCore))
    startEntry=$((numEntries*(absoluteCore)))
    stopEntry=$((numEntries* (absoluteCore+1)))
    ./templateSearch --wais ${sharedDir}/${absoluteCore} ${numCores} ${localCore} 1> ${sharedDir}/log/${localCore}.log 2>&1 &
done
    