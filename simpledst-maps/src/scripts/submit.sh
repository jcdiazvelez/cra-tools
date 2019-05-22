#!/bin/sh

# A simple sh script for submitting one job to the npx cluster.
# Simply log into npx.icecube.wisc.edu, run any command as you normally 
#  would but prepend "./submit-npx.sh" and the command will be run on the
#  npx cluster. For example:
#
# ./submit-npx.sh root -l -b -q macro.C
#
# This script will create directories to store your execution script, log files,
#  errors, and std output, so you need write permission in the local directory.

# This script creates a script to be executed and another script to submit it.
# The execution script must be available *at time of job execution!*, which may
#  not be until much later and so it's stored in a directory 'npx-execs'.
# You may occasionally want to 'rm -rf npx-*' directories if they get big.
# The submission script "2sub.sub" can be discarded immediately.

# This method of passing your job into a bash script as arguments may fail
#  if you have quotes or other special characters

  # Creating output directories
  mkdir -p npx-execs npx-logs npx-out npx-error

  # Creating execution script, do not delete until job has started!
  echo "#!/bin/bash" > npx-execs/npx-$$.sh
  echo "date" >> npx-execs/npx-$$.sh
  echo "hostname" >> npx-execs/npx-$$.sh
  echo "cd `pwd`" >> npx-execs/npx-$$.sh
  #echo "eval `/data/user/juancarlos/simpledst-maps/build/scripts/setup.sh`" >> npx-execs/npx-$$.sh
  echo "$@" >> npx-execs/npx-$$.sh
  echo "date" >> npx-execs/npx-$$.sh


  chmod +x npx-execs/npx-$$.sh

  # Creating condor submission script (ClassAd)
  echo "Universe  = vanilla" > 2sub.sub
  echo "Executable = npx-execs/npx-$$.sh" >> 2sub.sub
  echo "Log = npx-logs/npx-$$.log" >> 2sub.sub
  echo "Output = npx-out/npx-$$.out" >> 2sub.sub
  echo "Error = npx-error/npx-$$.error" >> 2sub.sub
  echo "Notification = NEVER" >> 2sub.sub 
  echo "Queue" >> 2sub.sub
  condor_submit 2sub.sub


