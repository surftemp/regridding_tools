#!/bin/bash
#
# Run reslicing jobs for a range of years
#
# Usage: ./launch.sh STARTYEAR ENDYEAR
#
STARTYEAR=$1
ENDYEAR=$2

for ((YEAR=$STARTYEAR;YEAR<=$ENDYEAR;YEAR+=1))
do
    echo running reslice job for $YEAR
    python reslice.py --input-sst /data/input/CDR_v3/Analysis/L4/v3.0.1/ --input-climatology /data/input/CDR_v3/Climatology/L4/v3.0.1/ --output-folder /data/input/CDR_v3/reslice --year $YEAR
done
