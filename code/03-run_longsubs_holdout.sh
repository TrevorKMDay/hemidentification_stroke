#!/bin/bash

# Set up

home=/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/code
input_all=group-all_parc-brainnetome_normed.pickle
input_testonly=group-testonly_parc-brainnetome_normed.pickle
input_allrows=group-allrows_parc-brainnetome_normed.pickle

hi_ver=0.4.3

# Get subs

testonly=input/group-testonly_parc-brainnetome_normed.csv

long_subs=$(grep "test_long" ${testonly} | cut -d',' -f1 | uniq)

echo "Found the following longitudinal ppts:"
echo "${long_subs}" | tr '\n' ' '

mkdir -p {results,predicted}/holdouts/

# Run loop

fpman=/Users/tkmd/MyDrive/Projects/apply_hemiconnectome/bin/fingerprinting_manager.py
loadd=/Users/tkmd/MyDrive/Projects/apply_hemiconnectome/bin/load_data.py

for outcome in hemi class ; do

    for s in ${long_subs} ;  do

        # Remove prefix for filenames
        s2=${s//sub-/}

        echo "${outcome}"
        model_file="input-hc_holdout-${s2}_outcome-${outcome}_model.pickle"

        docker run --rm \
            -v ${home}/input:/input/        \
            -v ${home}/results:/results     \
            hemidentification:${hi_ver}     \
                --verbose                   \
                --file_conflict compare     \
                -k sub "!${s}"           -- \
                train                       \
                    /input/${input_allrows}     \
                    "/results/holdouts/${model_file}"     \
                    ${outcome}

        test_out="input-testonly_heldout-${s2}_outcome-${outcome}_predictedon-all_results.csv"

        docker run --rm \
            -v ${home}/input:/input/            \
            -v ${home}/results:/results         \
            -v ${home}/predicted:/predicted/    \
            hemidentification:${hi_ver}         \
                --file_conflict compare         \
                --verbose                       \
                test                            \
                    "/results/holdouts/${model_file}"    \
                    /input/${input_allrows}    \
                    "/predicted/holdouts/${test_out}"

        # break

    done


    docker run --rm \
        -v ${home}/input:/input/            \
        -v ${home}/results:/results         \
        -v ${home}/predicted:/predicted/    \
        hemidentification:${hi_ver}         \
            --file_conflict compare         \
            --verbose                       \
            test                            \
                "/results/input-all_outcome-${outcome}_model.pickle"    \
                /input/${input_allrows}    \
                "/predicted/input-all_outcome-${outcome}_predictedon-allrows_results.csv"

done