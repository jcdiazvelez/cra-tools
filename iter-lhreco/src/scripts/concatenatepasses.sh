#!/bin/bash

eval `$HAWCSOFT/setup.sh --externals`
eval `$AERIE/install/bin/hawc-config --env-sh`

# $1 chunk number

chunk=$1

dataset="v202-sp-rc4-cuts"
basefolder="/data/scratch/backup/fiorino/hawc/local/$dataset"
folder="$basefolder/chunk${chunk}/"

if ls $folder/*pass*gz 1> /dev/null 2>&1; then

    if [ ! -d $folder/redundant ]; then
        mkdir $folder/redundant
    fi

    for b in {0..9};
    do
        for t in {000..359};
        do
           passfiles="$folder/*bin${b}*${t}_of_360_*pass*gz"

           if ls $passfiles 1> /dev/null 2>&1; then
               aerie-apps-combine-local-maps $passfiles -o $folder/hawclocal-v202-sp-rc4-cuts-chunk00${chunk}_bin${b}_N512_${t}_of_360.fits.gz
               mv $folder/$passfiles $folder/redundant/.
           fi

           echo "Analysis Bin $b. Time bin $t of 360 CONCATENTATED"

        done
    done
fi

echo "Finished"
# batch submit
#for chunk in {0003..0209}; do sbatch --mem=300mb --time=32:00:00 -o $MAPPROD_DIR/scripts/concatenate-hawclocal-v202-sp-rc4-cuts_chunk${chunk}.out concatenatepasses.sh $chunk; done
