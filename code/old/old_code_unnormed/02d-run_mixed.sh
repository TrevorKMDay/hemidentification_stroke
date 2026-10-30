#!/bin/bash

mkdir -p results/mixed/

inputs=$(find input/mixed/ -name "*.pickle")

# Both predict hemi alone (1 LD); or
# Class: HC LH/RH + LHS RH (2 LDs)
for i in hemi class ; do

    for input in ${inputs} ; do

        lpct=$(echo "${input}" | grep -o "Lpct-0.[0-9]")

        oname="results/mixed/HC_group-all_method-${lpct}_outcome-${i}_results.csv"

        python run_model.py                         \
            --drop_cols     PatientStatus           \
            --output_name   "${oname}"              \
            "${i}" "${input}" "within" "lda"

    done


done
