#!/bin/bash

ITERATION=$1

while read group;
do
    while read -r run condition;
    do
        for (( i=1; i<=4; ++i ))
        do
            swif2 add-job -workflow solid_pmt_test -name ${run:0:8}_${run:9:4}_Ch0${i}_Iter0${ITERATION} -account halla -partition production -os general -cores 1 -disk 1GB -ram 512MB -time 24hrs -stdout /farm_out/boyu/solid_pmt_test/${run:0:8}_${run:9:4}_Ch0${i}_Iter0${ITERATION}.out -stderr /farm_out/boyu/solid_pmt_test/${run:0:8}_${run:9:4}_Ch0${i}_Iter0${ITERATION}.err sh /work/halld2/home/boyu/solid_pmt_test/run_fit.sh ${run:0:8} ${run:9:4} ${i} ${ITERATION}
            # sh /work/halld2/home/boyu/solid_pmt_test/run_fit.sh ${run:0:8} ${run:9:4} ${i} ${ITERATION}
        done
    done < results/${group}/runs.txt
done < list.txt

swif2 run -workflow solid_pmt_test