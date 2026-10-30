#!/bin/bash

mkdir -p results/all/

input=input/group-all_parc-brainnetome_loocv_rZ.pickle

# Both predict hemi alone (1 LD); or
# Class: HC LH/RH + LHS RH (2 LDs)
for i in hemi class ; do

    for k in sv nn lda ; do

        oname="results/all/HC_group-all_method-${k}_outcome-${i}_results.csv"

        python run_model.py                     \
            --drop_cols     PatientStatus       \
            --output_name   ${oname}            \
            "${i}" ${input} "within" "${k}"

    done

done

