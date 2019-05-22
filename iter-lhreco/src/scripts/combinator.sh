#!/bin/bash
eval `$HAWCSOFT/setup.sh --externals`
eval `$AERIE/install/bin/hawc-config --env-sh`

#$1  analysis bin
#$2  time bin start
#$3  time bin stop

dataset="v202-sp-rc4-cuts"
basefolder="/data/scratch/backup/fiorino/hawc/local/${dataset}/"

for t in `seq -f "%03g" $2 $3`; 
do 
  if [ ! -e ${basefolder}/combined/hawclocal-${dataset}_bin${1}_N512_${t}_of_360.fits.gz ]; then
    aerie-apps-combine-local-maps ${basefolder}/chunk*/hawclocal-${dataset}-chunk*_bin${1}_N512_${t}_of_360.fits.gz -o ${basefolder}/combined/hawclocal-${dataset}_bin${1}_N512_${t}_of_360.fits.gz
  else
    echo "${basefolder}/combined/hawclocal-${dataset}_bin${1}_N512_${t}_of_360.fits.gz exists"
  fi
done 

# to run 
#for b in {0..9}; do for i in {0..330..30}; do j=$((i + 29)); `nohup sh combinator.sh $b $i $j > nohup_bin${b}_t${i}.out &`; done; done

# to submit to cluster
#for b in {1..9}; do for i in {0..330..30}; do j=$((i + 29)); sbatch --mem=300mb --time=32:00:00 -o $MAPPROD_DIR/scripts/combinator-hawclocal-v202-sp-rc4-cuts_bin${b}_t${i}.out combinator.sh $b $i $j; done; done

# to check finished maps
#for b in {0..9}; do for i in {000..359}; do if [ -e hawclocal-v202-sp-rc4-cuts_bin${b}_N512_${i}_of_360.fits.gz ]; then foo="foo"; else echo "hawclocal-v202-sp-rc4-cuts_bin${b}_N512_${i}_of_360.fits.gz";fi; done; done
