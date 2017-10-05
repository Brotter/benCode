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
sharedDir="/home/brotter/nfsShared/results/cluster/"${startTime}
mkdir ${sharedDir}
mkdir ${sharedDir}"/log"

numCores=256
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

for localCore in `seq ${startCore} $((startCore+63))`; do
    echo "Starting "${localCore}
    root cluster.C\(${numCores},${localCore},\"${sharedDir}\"\)  1> ${sharedDir}/log/${localCore}.log 2>&1 &
done
