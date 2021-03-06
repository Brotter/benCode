#!/bin/bash

################################################
#
#  Ben Rotter - September 2016
#
#  This code should run on the four ANITA clusters
#  at UH to run whatever code faster
#
#################################################

startTime=`date +%m.%d.%y_%Hh%Mm`
sharedDir="clusteringOutputs/"${startTime}
mkdir ${sharedDir}
mkdir ${sharedDir}"/log"

numCores=64
for localCore in `seq 0 64`; do
    echo "Starting "${localCore}" to "${sharedDir}
    root -b cluster.C\(\"clusterEvents\",\"${1}\",${numCores},${localCore},\"${sharedDir}\"\)  1> ${sharedDir}/log/${localCore}.log 2>&1 &
done
