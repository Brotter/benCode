#!/bin/bash

################################################
#
#  Ben Rotter - November 2017
#
#  This code should run on the four ANITA clusters
#  at UH to run whatever code faster
#
#################################################

startTime=`date +%m.%d.%y_%Hh%Mm`
echo ${startTime}
sharedDir="instStokes/"
mkdir ${sharedDir}
mkdir ${sharedDir}"/log"


numCores=64
for localCore in `seq 0 64`; do
    echo "Starting "${localCore}" to "${sharedDir}
    root -b newStokes.C\(\"waisEvents.root\",64,${localCore},false,\"${sharedDir}\"\)  1> ${sharedDir}/log/${localCore}.log 2>&1 &
done
