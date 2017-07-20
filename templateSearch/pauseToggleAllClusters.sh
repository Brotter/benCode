#/bin/bash

#
#
# This will either pause or unpause the simulation, putting the ROOT files into a stable state
#
#
#


ssh anitai   'killall -SIGINT templateSearch'
ssh anitaii  'killall -SIGINT templateSearch'
ssh anitaiii 'killall -SIGINT templateSearch'
ssh anitaiv  'killall -SIGINT templateSearch'


