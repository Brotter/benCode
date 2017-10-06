#/bin/bash

#
#
# So that I don't have to log into each server and update it individually.
#
# Be careful because this could totally fail
#


ssh anitai 'cd ~/anita16/anitaBuildTool; ./updateComponents.sh 0; make; make install' &
ssh anitaii 'cd ~/anita16/anitaBuildTool; ./updateComponents.sh 0; make; make install' &
ssh anitaiii 'cd ~/anita16/anitaBuildTool; ./updateComponents.sh 0; make; make install' &
ssh anitaiv 'cd ~/anita16/anitaBuildTool; ./updateComponents.sh 0; make; make install' &


