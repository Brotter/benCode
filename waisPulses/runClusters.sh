#!/bin/bash

################################################
#
#  Ben Rotter - September 2016
#
#  This code should run on the four ANITA clusters
#  at UH to run the interferometry code faster
#
#################################################



numWaisPulses=118799

#There are 4 servers with 64 cores each, so 256 total
numEvents=$((numWaisPulses/256))

if [ `hostname | cut -d"." -f1` == "anitai" ]; then
    startSeq = 0
    stopSeq = 63
elif [ `hostname | cut -d"." -f1` == "anitaii" ]; then
    startSeq = 64
    stopSeq = 127
elif [ `hostname | cut -d"." -f1` == "anitaiii" ]; then
    startSeq = 128
    stopSeq = 191
elif [ `hostname | cut -d"." -f1` == "anitaiv" ]; then
    startSeq = 192
    stopSeq = 255
fi

for seq in `seq ${startSeq} ${stopSeq}`; do
    ./waisNewCorrelator 330 360 $((numEvents*startSeq)) $((numEvents*(startSeq+1))) ${seq} & >> log/${seq}.log
done
    