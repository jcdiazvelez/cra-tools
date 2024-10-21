#!/bin/bash

function run_makemaps {
    year=$1
    options=$2

    # Run the script to produce the dag files.
    command="python ../../scripts/make-local-maps-IT.py.in -c ITpass2 --method anti -o /data/user/gagrawal/anti/tier3 $options"
    $command

    # Submit the dag file for each year to condor. This will produce 360 fits files for each day in the tier.
    commandtwo="condor_submit_dag ./submit_$year/ITpass2.make-local-maps.dag"
    $commandtwo
}

#First date - 20110513, last date, 20220803, how to balance between run years, when dates overlap?
# Define the options for each year
declare -A options_per_year
options_per_year[2011]="-d 20110513 20120512 -i --smin 10 --smax 17 --submit_dir ./submit_2011"
options_per_year[2012]="-d 20120513 20130512 -i --smin 9 --smax 16 --submit_dir ./submit_2012"
options_per_year[2013]="-d 20130513 20140512 -i --smin 9 --smax 16 --submit_dir ./submit_2013"
options_per_year[2014]="-d 20140513 20150512 -i --smin 8 --smax 15 --submit_dir ./submit_2014"
options_per_year[2015]="-d 20150513 20160512 -i --smin 8 --smax 15 --submit_dir ./submit_2015"
options_per_year[2016]="-d 20160513 20170512 -i --smin 7 --smax 14 --submit_dir ./submit_2016"
options_per_year[2017]="-d 20170513 20180512 -i --smin 7 --smax 14 --submit_dir ./submit_2017"
options_per_year[2018]="-d 20180513 20190512 -i --smin 6 --smax 13 --submit_dir ./submit_2018"
options_per_year[2019]="-d 20190513 20200512 -i --smin 6 --smax 13 --submit_dir ./submit_2019"
options_per_year[2020]="-d 20200513 20210512 -i --smin 5 --smax 12 --submit_dir ./submit_2020"
options_per_year[2021]="-d 20210513 20220512 -i --smin 5 --smax 12 --submit_dir ./submit_2021"


# Iterate over 10 years
for year in {2011..2021}; do
    options=${options_per_year[$year]}

    # Run make-local-maps-IT.py.in for each year with the correct nStation boundaries.
    run_makemaps $year "$options"

done

