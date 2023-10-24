#!/bin/bash

MODE=$1

while read -r run channel;
do
    if [ ${channel} == "all" ]
    then
        for (( i=1; i<=5; ++i ))
        do
            if [ ${MODE} == "farm" ]
            then
                swif2 add-job -workflow solid_pmt_test -name ${run:0:8}_${run:9:4}_Ch0${i} -account halla -partition production -os general -cores 1 -disk 1GB -ram 2GB -time 4hrs -stdout /farm_out/boyu/solid_pmt_test/${run:0:8}_${run:9:4}_Ch0${i}.out -stderr /farm_out/boyu/solid_pmt_test/${run:0:8}_${run:9:4}_Ch0${i}.err sh /work/halld2/home/boyu/solid_pmt_test/run_fit.sh ${run:0:8} ${run:9:4} ${i}
            elif [ ${MODE} == "local" ]
            then
                sh /work/halld2/home/boyu/solid_pmt_test/run_fit.sh ${run:0:8} ${run:9:4} ${i}
            fi
        done
    else
        if [ ${MODE} == "farm" ]
        then
            swif2 add-job -workflow solid_pmt_test -name ${run:0:8}_${run:9:4}_Ch0${channel} -account halla -partition production -os general -cores 1 -disk 1GB -ram 2GB -time 4hrs -stdout /farm_out/boyu/solid_pmt_test/${run:0:8}_${run:9:4}_Ch0${channel}.out -stderr /farm_out/boyu/solid_pmt_test/${run:0:8}_${run:9:4}_Ch0${channel}.err sh /work/halld2/home/boyu/solid_pmt_test/run_fit.sh ${run:0:8} ${run:9:4} ${channel}
        elif [ ${MODE} == "local" ]
        then
            sh /work/halld2/home/boyu/solid_pmt_test/run_fit.sh ${run:0:8} ${run:9:4} ${channel}
        fi
    fi
    
done < list_fit.txt

swif2 run -workflow solid_pmt_test