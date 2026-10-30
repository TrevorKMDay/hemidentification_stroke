#!/bin/bash

home=/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/code

input_all=group-all_parc-brainnetome_normed.pickle

hi_ver=0.4.3

# Predict hemi from HC data ====================================================

outcome="hemi"

echo "${outcome}"
oname="input-hc_outcome-hemi_model.pickle"

docker run --rm \
    -v ${home}/input:/input/        \
    -v ${home}/results:/results     \
    hemidentification:${hi_ver}     \
        --file_conflict compare     \
        -k PatientStatus HC         \
           status train_and_test -- \
        train                       \
            /input/${input_all}     \
            /results/${oname}       \
            ${outcome}


# Predict all outcomes from everyone ===========================================

# Both predict hemi alone (1 LD); or
# Class: HC LH/RH + LHS RH (2 LDs)
for outcome in PatientStatus hemi class ; do

    # This is the self-to-self classification

    echo "${outcome}"
    oname="input-all_outcome-${outcome}_model.pickle"

    docker run --rm \
        -v ${home}/input:/input/        \
        -v ${home}/results:/results     \
        hemidentification:${hi_ver}     \
            --file_conflict compare     \
            -k status train_and_test -- \
            train                       \
                /input/${input_all}     \
                /results/${oname}       \
                ${outcome}

done

class_model=input-all_outcome-class_model.pickle

docker run --rm \
    -v ${home}/input:/input/        \
    -v ${home}/results:/results     \
    hemidentification:${hi_ver}     \
        --file_conflict compare     \
        test                        \
            /results/${class_model} \
            /input/${input_all}     \
            /results/input-all_model-class.all_outcome.csv
