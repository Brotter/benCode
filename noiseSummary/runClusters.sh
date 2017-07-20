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
sharedDir="/home/brotter/nfsShared/results/noiseSummary/"${startTime}
mkdir ${sharedDir}
mkdir ${sharedDir}"/log"

if [ `hostname | cut -d"." -f1` == "anitaI" ]; then
    startRun=130; stopRun=207
elif [ `hostname | cut -d"." -f1` == "anitaII" ]; then
    startRun=208; stopRun=285
elif [ `hostname | cut -d"." -f1` == "anitaIII" ]; then
    startCore=286; stopRun=362
elif [ `hostname | cut -d"." -f1` == "anitaIV" ]; then
    startCore=362; stopRun=440
else
    echo "The server isn't an anita cluster server, so you shouldn't use this script"	

    exit
fi

for run in `seq ${startRun} ${endRun}`
    nice -n 10 ./noiseSummary ${run} ${sharedDir}/${run} 1> ${sharedDir}/log/${run}.log 2>&1 &
done
    