#!/bin/bash

function run_makemaps {
    year=$1
    options=$2

    # Run the command to produce the dag scripts.
    command="python ../../scripts/make-local-maps-IT.py.in -c ITpass2 --method anti -o /data/user/gagrawal/anti/tier1 $options"
    $command

    # Submit the dag files for each year to condor. This will produce 360 fits files for each day in this tier.
    commandtwo="condor_submit_dag ./submit_$year/ITpass2.make-local-maps.dag"
    $commandtwo
}

#What dates should we use?
# Define the options for each year
declare -A options_per_year
options_per_year[2011]="-d 20110513 20120512 -i --smin 3 --smax 5 --submit_dir ./submit_2011"
options_per_year[2012]="-d 20120513 20130512 -i --smin 3 --smax 5 --submit_dir ./submit_2012"
options_per_year[2013]="-d 20130513 20140512 -i --smin 3 --smax 5 --submit_dir ./submit_2013"
options_per_year[2014]="-d 20140513 20150512 -i --smin 3 --smax 5 --submit_dir ./submit_2014"

# Iterate over 10 years
for year in {2011..2014}; do
    options=${options_per_year[$year]}

    # Run make-local-maps-IT.py.in for each year with the correct nStation boundaries.
    run_makemaps $year "$options"

done

