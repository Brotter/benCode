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
saveDir=${HOME}"/anita16/benCode/templateSearch/09.27.17_19h"
mkdir ${saveDir}
mkdir ${saveDir}"/log"

#19 dead runs broken up into servers
if [ `hostname | cut -d"." -f1` == "anitaI" ]; then
    cores="4877 89 114 128"
elif [ `hostname | cut -d"." -f1` == "anitaII" ]; then
    cores="139 144 157 170 189"
elif [ `hostname | cut -d"." -f1` == "anitaIII" ]; then
    cores="191 199 204 220 224"
elif [ `hostname | cut -d"." -f1` == "anitaIV" ]; then
    cores="226 228 251 255"
else
    echo "The server isn't an anita cluster server, so you shouldn't use this script"	
    cores="226 227 228 251 255"
#    exit
fi

numEntries=307151 #total num of events / 256 cores


#5 runs per server, 64 cores?  lets just do 12 cores per run (normally ~54 hours goes down to like 4.5)
numCores=12
for reCore in ${cores}; do
    reCoreStart=$((numEntries*(reCore)))
    reCoreStop=$((numEntries*(reCore+1)))
    for localCore in `seq 1 ${numCores}`; do
	localCoreStart=$((reCoreStart+(numEntries/12)*(localCore-1)))
	localCoreStop=$((reCoreStart+(numEntries/12)*(localCore)))
	if [ $localCore -eq ${numCores} ]; then localCoreStop=${reCoreStop};fi
	preName=${reCore}_${localCore}
	echo "reCore:"${reCore}" localCoreStart:"${localCoreStart}" localCoreStop:"${localCoreStop}
	nice -n 11 ./templateSearch --entry ${saveDir}/${preName} ${localCoreStart} ${localCoreStop} 1> ${saveDir}/log/${preName}.log 2>&1 &
    done
done
    