#/bin/bash

#
#
# One of the last runs I did in September had a bunch of dead runs, this one throws them onto some servers so it will finish tonight
#
#
#



ssh anitai   'cd anita16/benCode/templateSearch/; git pull; make clean; make; ./runClusters.sh' 
ssh anitaii  'cd anita16/benCode/templateSearch/; git pull; make clean; make; ./runClusters.sh' 
ssh anitaiii 'cd anita16/benCode/templateSearch/; git pull; make clean; make; ./runClusters.sh' 
ssh anitaiv  'cd anita16/benCode/templateSearch/; git pull; make clean; make; ./runClusters.sh' 


