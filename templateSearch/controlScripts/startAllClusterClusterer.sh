#/bin/bash

#
#
# the clustering thing has started to take forever, so lets get all four servers to chew on it
#
#
#



#ssh anitai   'cd anita16/benCode/templateSearch/macros/; git pull; ./runClusterCluster.sh' 
ssh anitaii  'cd anita16/benCode/templateSearch/macros/; git pull; ./runClusterCluster.sh' 
ssh anitaiii 'cd anita16/benCode/templateSearch/macros/; git pull; ./runClusterCluster.sh' 
ssh anitaiv  'cd anita16/benCode/templateSearch/macros/; git pull; ./runClusterCluster.sh' 


