#!/bin/bash

################################################
#
#  Ben Rotter - September 2016
#
#  This code should run on the four ANITA clusters
#  at UH to run whatever code faster
#
#################################################

mkdir log

numCores=64
for localCore in `seq 0 63`; do
    echo "Starting "${localCore}
    root -b makeSNRFile.C\(${numCores},${localCore}\)  1> log/${localCore}.log 2>&1 &
done
