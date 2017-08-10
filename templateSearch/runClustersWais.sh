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

numCores=32

for localCore in `seq 0 $((numCores-1))`; do
    absoluteCore=$((localCore+startCore))
    startEntry=$((numEntries*(absoluteCore)))
    stopEntry=$((numEntries* (absoluteCore+1)))
    ./templateSearch --wais ${sharedDir}/${absoluteCore} ${numCores} ${localCore} 1> ${sharedDir}/log/${localCore}.log 2>&1 &
done
    