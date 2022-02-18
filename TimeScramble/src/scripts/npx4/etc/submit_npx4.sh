#!/bin/sh

# A simple sh script for submitting one job to the npx4 cluster.
# Simply log into npx4.icecube.wisc.edu, run any command as you normally 
#  would but prepend "./submit-npx4.sh" and the command will be run on the
#  npx4 cluster. For example:
#
# ./submit-npx4.sh root -l -b -q macro.C
#
# This script will create directories to store your execution script, log files,
#  errors, and std output, so you need write permission in the local directory.

# This script creates a script to be executed and another script to submit it.
# The execution script must be available *at time of job execution!*, which may
#  not be until much later and so it's stored in a directory 'npx4-execs'.
# You may occasionally want to 'rm -rf npx4-*' directories if they get big.
# The submission script "2sub.sub" can be discarded immediately.

# This method of passing your job into a bash script as arguments may fail
#  if you have quotes or other special characters

  # Creating output directories
  mkdir -p npx4-execs npx4-logs npx4-out npx4-error

  # Creating execution script, do not delete until job has started!
  echo "#!/bin/bash" > npx4-execs/npx4-$$.sh
  echo "date" >> npx4-execs/npx4-$$.sh
  echo "hostname" >> npx4-execs/npx4-$$.sh
  echo "cd `pwd`" >> npx4-execs/npx4-$$.sh
  echo "eval `/cvmfs/icecube.opensciencegrid.org/setup.sh`" >> npx4-execs/npx4-$$.sh
  echo "$@" >> npx4-execs/npx4-$$.sh
  echo "date" >> npx4-execs/npx4-$$.sh


  chmod +x npx4-execs/npx4-$$.sh

  # Creating condor submission script (ClassAd)
  echo "Universe  = vanilla" > 2sub.sub
  echo "Executable = npx4-execs/npx4-$$.sh" >> 2sub.sub
  echo "Log = npx4-logs/npx4-$$.log" >> 2sub.sub
  echo "Output = npx4-out/npx4-$$.out" >> 2sub.sub
  echo "Error = npx4-error/npx4-$$.error" >> 2sub.sub
  echo "Notification = NEVER" >> 2sub.sub 
  echo "Queue" >> 2sub.sub
  condor_submit 2sub.sub


