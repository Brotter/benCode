#!/bin/bash

################################################
#
#  Ben Rotter - September 2016
#
#  This code should run on the four ANITA clusters
#  at UH to run whatever code faster
#
#################################################


sharedDir=${ANITA3_RESULTSDIR}"/templateSearch/09.27.17_19h/reCalcKeyValues/"
mkdir ${sharedDir}
mkdir ${sharedDir}"/log"

numCores=64
for localCore in `seq 0 64`; do
    echo "Starting "${localCore}
    ./reCalcKeyValues ${numCores} ${localCore} 1> ${sharedDir}/log/${localCore}.log 2>&1 &
done
