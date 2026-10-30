#!/bin/bash

home=/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/code
models=$(find results/ -name "input-all_*_model.pickle" | sort)

echo "Found models:"
echo "${models}" | tr ' ' '\n'
echo

# Update me to change docker version
hi_ver=0.4.3

input_all=group-all_parc-brainnetome_normed.pickle

hconly_hemi_model=input-hc_outcome-hemi_model.pickle
class_model=input-all_outcome-class_model.pickle

# For testing
# fp=/Users/tkmd/MyDrive/Projects/apply_hemiconnectome/bin/fingerprinting_manager.py

# Predict LHS hemisphere status from models trained only on HCs

oname="input-hc_outcome-hemi_predictedon-lhs_results.csv"

docker run --rm                             \
    -v ${home}/input:/input                 \
    -v ${home}/results:/results             \
    -v ${home}/predicted:/predicted         \
    hemidentification:${hi_ver}  -v         \
        --file_conflict compare             \
        --keep_groups PatientStatus LHS --  \
        test                                \
            "/results/${hconly_hemi_model}" \
            /input/${input_all}             \
            "/predicted/${oname}"

# Loop over the models predicting class, hemi, and patient status

for model in ${models} ; do

    m=$(basename "${model}")
    oname=${m//_model.pickle/_predictedon-all_results.csv}

    echo "${m}"
    echo "${oname}"

    docker run --rm \
        -v ${home}/input:/input         \
        -v ${home}/results:/results     \
        -v ${home}/predicted:/predicted \
        hemidentification:${hi_ver}     \
            --file_conflict compare     \
            test                        \
                "/results/${m}"         \
                /input/${input_all}    \
                "/predicted/${oname}"

done

# TEST ONLY

# Run the held-out test-only data on the regular data

input_testonly=group-testonly_parc-brainnetome_normed.pickle

docker run --rm \
    -v ${home}/input:/input/        \
    -v ${home}/results:/results     \
    -v ${home}/predicted:/predicted \
    hemidentification:${hi_ver}     \
        --file_conflict compare     \
        test                        \
            /results/${class_model} \
            /input/${input_testonly}     \
            /predicted/input-testonly_outcome-class_predictedon-all_results.csv