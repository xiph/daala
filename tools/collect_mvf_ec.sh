#!/bin/sh

set -e

VIDEOS=($*)
RATES=("1 3 5 7 11 16 25 37 55 81 122 181 270")
VNAMES=()

for vid in ${VIDEOS[@]}; do
    vname=$(basename ${vid})
    vname=${vname/%_*/}
    VNAMES+=(${vname})
done

parallel -j4 'export OD_EC_ACCT_SUFFIX={1/.}-{2} ; examples/encoder_example -v {2} {1} -o /dev/null' ::: ${VIDEOS[@]} ::: ${RATES[@]}
