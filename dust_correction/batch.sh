#!/bin/bash

year="1981"
while [ $year -lt "2017" ];
do
    echo $year
    python correct.py --data-dir `pwd` --input-dir /data/input/CDR_v2/Analysis/L4/v2.1 --input-year $year --output-dir /data2/input/dust/CDR_v2/Analysis/L4/v2.1
    ((year=year+1))
done

year="2017"
while [ $year -lt "2020" ];
do
    echo $year
    python correct.py --data-dir `pwd` --input-dir /data/input/ICDR_v2/Analysis/L4/v2.0  --input-year $year --output-dir /data2/input/dust/ICDR_v2/Analysis/L4/v2.0
    ((year=year+1))
done
