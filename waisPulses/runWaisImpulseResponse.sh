#!/bin/bash

for i in `seq 0 64`; do
    ./waisImpulseResponse $((i*1856)) $(((i+1)*1856)) waisImpulseResponses/${i}.root >> waisImpulseResponses/${i}.log &
done
