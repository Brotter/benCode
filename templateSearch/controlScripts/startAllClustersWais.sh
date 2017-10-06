#/bin/bash

#
#
# Just an attempt to get all the servers to update and start with one command, requires the other commands I've made
#
#
#



ssh anitai   'cd anita16/benCode/templateSearch/; git pull; make clean; make; ./runClustersWais.sh' 
ssh anitaii  'cd anita16/benCode/templateSearch/; git pull; make clean; make; ./runClustersWais.sh' 
ssh anitaiii 'cd anita16/benCode/templateSearch/; git pull; make clean; make; ./runClustersWais.sh' 
ssh anitaiv  'cd anita16/benCode/templateSearch/; git pull; make clean; make; ./runClustersWais.sh' 


