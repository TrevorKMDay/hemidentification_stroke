#!/bin/bash

home=/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/test_ln_cutoff
mkdir -p results/

hi_ver=0.4.3
file_conflict=compare
outcome=is_lhln_lesioned

for i in $(seq -w 0 5 50) ; do

    input=input-lhs_lesion-0.${i}.pickle
    oname=input-lhs_lesion-0.${i}_results.csv

    docker run --rm \
        -v ${home}/inputs:/input/                   \
        -v ${home}/results/:/results                \
        hemidentification:${hi_ver}  -v             \
            --file_conflict ${file_conflict}        \
            --keep_groups                           \
                status        train_and_test --     \
            loocv                                   \
                /input/"${input}"                   \
                /results/"${oname}"                 \
                ${outcome}

    # break

done
