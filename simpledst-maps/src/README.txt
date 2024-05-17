INSTRUCTIONS TO PRODUCE SKYMAPS

Step 0: Editing the code.
* Edit @USER_DIR@ in scripts/make-local-maps-IT.py.in, scripts/merge-reco.py.in, and the four shell scripts in run_burnsample/*/
to the user-appropriate directories. 
  - fpath should be a directory containing only the burnsample root files
  - outdir should be unique for each tier, and will contain the local and final skymaps.
  - the cmd line should hold the location of your make-local-maps file

* Edit @CVMFS_ROOTBASE@ in make-local-maps-IT.py.in and merge-reco.py.in.
  - default='/cvmfs/icecube/opensciencegrid.org/py3-v4.2.1',

* Specify which reconstruction to use in SimpleDST.cc.
  - For Tier 1, we use ShowerPlane. So, on line XX, comment out 'reco = Laputop' and uncomment 'reco = ShowerPlane.'
  - For Tiers 2, 3 and 4, we use Laputop. Do the inverse of the above.

* Rebuild the project with the following commands. Make sure you are in Cobalt and have run CVMFS (version py3-v4.2.1)
  - cd ~/cra-tools-gb-edits/simpledst-maps
  - mkdir build
  - cd build
  - cmake ../src
  - make

To Run:
Step 1: Producing Local Skymaps
* SSH into the condor machine and run CVMFS (version py3-v4.2.1). 

* Go into run_burnsample. For the tier you want to run, go into that directory and run the shell script. Make sure
that the output directory is empty of any existing folders, as the code does not overwrite any existing files.

Step 2: Combining Local Skymaps 
* SSH into cobalt and run CVMFS (version py3-v4.2.1).

* Go into run_bursample and into the directory for the tier you want to process.

* For each tier, run the below command (edit the user directory to the output directory for each tier):
 - python ../../scripts/merge-reco.py.in -c ITpass2 -o /data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output/@USER_DIR@ --sd --reco --submit_dir ./submit
