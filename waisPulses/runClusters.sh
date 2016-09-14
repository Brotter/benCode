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

if [ `hostname | cut -d"." -f1` == "anitaI" ]; then
    startSeq=0
    stopSeq=63
elif [ `hostname | cut -d"." -f1` == "anitaII" ]; then
    startSeq=64
    stopSeq=127
elif [ `hostname | cut -d"." -f1` == "anitaIII" ]; then
    startSeq=128
    stopSeq=191
elif [ `hostname | cut -d"." -f1` == "anitaIV" ]; then
    startSeq=192
    stopSeq=255
else
    echo "The server isn't an anita cluster server, so you shouldn't use this script"
    exit
fi

for i in `seq ${startSeq} ${stopSeq}`; do
    nice ./waisNewCorrelator 330 360 $((numEvents*i)) $((numEvents*(i+1))) ${i} 1> log/${i}.log 2>&1 &
done
    