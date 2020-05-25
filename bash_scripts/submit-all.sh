#!/bin/bash

# submit-all.sh
#
# Usage:    submit-all.sh
#

EXTRACTOR=$PWD/bash_scripts/extract_HRMIP_proj.sh
OUTPUTS_DIR=/home/users/pcos/exs

queue=short-serial
for i in {2015..2050}; do

    day=$i
    echo "[INFO] Submitting job to LOTUS for date: $day"
    bsub -q $queue -W 00:05 -We 00:01 -o $OUTPUTS_DIR/${day}.%J.out -e $OUTPUTS_DIR/${day}.%J.err $EXTRACTOR $day 

done
