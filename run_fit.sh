#!/bin/bash

DATE=$1
TIME=$2
READOUT=$3
ITERATION=$4

source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh
gxenv $HALLD_VERSIONS/version.xml

root -b -q -l /work/halld2/home/boyu/solid_pmt_test/pmt_fit_adc.cc'('${DATE}', '${TIME}', '${READOUT}', '${ITERATION}')'