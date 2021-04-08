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
    python reslice.py --input-cci /data2/input/dust/CDR_v2/Analysis/L4/v2.1/ --input-c3s /data2/input/dust/ICDR_v2/Analysis/L4/v2.0/ --input-climatology /data2/input/dust/CDR_v2/Climatology/L4/v2.1/ --output-folder /data2/input/dust/reslice --year $YEAR
done
