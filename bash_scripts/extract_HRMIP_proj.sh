#!/bin/bash

# extract-era-data.sh
#
# Usage:    extract-era-data.sh <YYYYMMDD>
#
day=$1
YYYY=$(echo $day | cut -c1-4)

OUTPUT_DIR=/home/users/pcos/exs

VAR_ID=tas

files=/badc/cmip6/data/CMIP6/HighResMIP/EC-Earth-Consortium/EC-Earth3P-HR/highres-future/r1i1p2f1/Amon/tas/gr/latest/tas_Amon_EC-Earth3P-HR_highres-future_r1i1p2f1_gr_$YYYY\01-$YYYY\12.nc 

module load jaspy

for INPUT_FILE in $files; do

    echo "[INFO] Subsetting: $INPUT_FILE"
    fname=$(basename $INPUT_FILE)
    OUTPUT_FILE=$OUTPUT_DIR/$fname
    cdo selname,$VAR_ID $INPUT_FILE $OUTPUT_FILE

done

