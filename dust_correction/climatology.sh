#!/bin/bash

DOY="1"
while [ $DOY -lt "366" ];
do
    echo $DOY
    python climatology.py --doy $DOY --window 2 /data2/input/dust/CDR_v2/Analysis/L4/v2.1 /data2/input/dust/CDR_v2/Climatology/L4/v2.1/
    ((DOY=DOY+1))
done
