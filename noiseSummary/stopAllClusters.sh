#/bin/bash

#
#
# For when I mess up and have to stop everything
#
#
#


ssh anitai   'killall -SIGTERM noiseSummary'
ssh anitaii  'killall -SIGTERM noiseSummary'
ssh anitaiii 'killall -SIGTERM noiseSummary'
ssh anitaiv  'killall -SIGTERM noiseSummary'


