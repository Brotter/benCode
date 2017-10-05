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

numCores=64
for localCore in `seq 0 $((numCores-1))`; do
    root macros/cluster.C\(${numCores},${localCore},\"${sharedDir}\"\)  1> ${sharedDir}/log/${localCore}.log 2>&1 &
done
