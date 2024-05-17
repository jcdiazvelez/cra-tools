#!/bin/bash

function run_makemaps {
    year=$1
    options=$2

    # Run the script to produce the dag files.
    command="python ../scripts/make-local-maps-IT.py.in -c ITpass2 -o /data/ana/CosmicRay/Anisotropy/IceTop/ITpass2/output/outputg/tier2_90bs --sd $options"
    $command

    # Submit the dag file for each year to conodr. This will create 360 fits files for each day in this tier.
    commandtwo="condor_submit_dag ./submit_$year/ITpass2.make-local-maps.dag"
    $commandtwo
}

# Define the options for each year
declare -A options_per_year
options_per_year[2011]="-d 20110514 20120512 -i --smin 5 --smax 10 --submit_dir ./submit_2011"
options_per_year[2012]="-d 20120513 20130430 -i --smin 5 --smax 9 --submit_dir ./submit_2012"
options_per_year[2013]="-d 20130501 20140506 -i --smin 5 --smax 9 --submit_dir ./submit_2013"
options_per_year[2014]="-d 20140507 20150515 -i --smin 5 --smax 8 --submit_dir ./submit_2014"


# Iterate over 10 years
for year in {2011..2014}; do
    options=${options_per_year[$year]}

    # For each year, run make-local-maps-IT.py.in with the correct nStation boundaries.
    run_makemaps $year "$options"

done
