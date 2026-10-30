#!/bin/bash

home=/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/code
mkdir -p results/loocv/

input_all=group-all_parc-brainnetome_normed.pickle

hi_ver=0.4.3

# for container testing
# fpman=/Users/tkmd/MyDrive/Projects/apply_hemiconnectome/bin/fingerprinting_manager.py
# loadd=/Users/tkmd/MyDrive/Projects/apply_hemiconnectome/bin/load_data.py

file_conflict=strict

# Hemisphere =====

outcome="hemi"
oname=input-HC_outcome-hemi_loocv.csv

# Base LOOCV on good data
docker run --rm \
    -v ${home}/input:/input/                    \
    -v ${home}/results/loocv:/results           \
    hemidentification:${hi_ver}  -v             \
        --file_conflict ${file_conflict}                 \
        --keep_groups                           \
            PatientStatus HC                    \
            status        train_and_test --     \
        loocv                                   \
            /input/${input_all}                 \
            /results/${oname}                   \
            ${outcome}

# exit÷

oname=input-all_outcome-hemi_loocv.csv

docker run \
    -v ${home}/input:/input/            \
    -v ${home}/results/loocv:/results   \
    hemidentification:${hi_ver}         \
        --file_conflict ${file_conflict}         \
        --keep_groups                   \
            status train_and_test --    \
        loocv                           \
            /input/${input_all}         \
            /results/${oname}           \
            ${outcome}

# Patient Status =====

outcome=class
oname=input-all_outcome-class_loocv.csv

docker run \
    -v ${home}/input:/input/            \
    -v ${home}/results/loocv:/results   \
    hemidentification:${hi_ver}         \
        --file_conflict ${file_conflict}         \
        --keep_groups                   \
            status train_and_test --    \
        loocv                           \
            /input/${input_all}         \
            /results/${oname}           \
            ${outcome}

outcome=PatientStatus
oname=input-all_outcome-patient_loocv.csv

docker run \
    -v ${home}/input:/input/            \
    -v ${home}/results/loocv:/results   \
    hemidentification:${hi_ver}         \
        --file_conflict ${file_conflict}         \
        --keep_groups                   \
            status train_and_test --    \
        loocv                           \
            /input/${input_all}         \
            /results/${oname}           \
            ${outcome}

oname=input-allRHonly_outcome-patient_loocv.csv

docker run \
    -v ${home}/input:/input/            \
    -v ${home}/results/loocv:/results   \
    hemidentification:${hi_ver}         \
        --file_conflict ${file_conflict}         \
        --keep_groups                   \
            status train_and_test       \
            hemi   R              --    \
        loocv                           \
            /input/${input_all}         \
            /results/${oname}           \
            ${outcome}

# Stroke =====

outcome=logLV_g
oname=input-lhs_outcome-lLVg_loocv.csv

docker run \
    -v ${home}/input:/input/            \
    -v ${home}/results/loocv:/results   \
    hemidentification:${hi_ver}         \
        --file_conflict ${file_conflict}         \
        --keep_groups                   \
            status train_and_test       \
            PatientStatus LHS  --       \
        loocv                           \
            /input/${input_all}         \
            /results/${oname}           \
            ${outcome}

outcome=WAB_value_g
oname=input-lhs_outcome-WABg_loocv.csv

docker run \
    -v ${home}/input:/input/            \
    -v ${home}/results/loocv:/results   \
    hemidentification:${hi_ver}         \
        --file_conflict ${file_conflict}         \
        --keep_groups                   \
            status train_and_test       \
            PatientStatus LHS  --       \
        loocv                           \
            /input/${input_all}         \
            /results/${oname}           \
            ${outcome}

# Residualized WAB =====

outcome=residual_wab_group
oname=input-lhs_outcome-rWABg_loocv.csv

docker run \
    -v ${home}/input:/input/            \
    -v ${home}/results/loocv:/results   \
    hemidentification:${hi_ver}         \
        --file_conflict ${file_conflict}         \
        --keep_groups                   \
            status train_and_test       \
            PatientStatus LHS  --       \
        loocv                           \
            /input/${input_all}         \
            /results/${oname}           \
            ${outcome}

# LN lesion WAB =====

outcome=is_ln_lesioned
oname=input-lhs_outcome-LNyn_loocv.csv

docker run \
    -v ${home}/input:/input/            \
    -v ${home}/results/loocv:/results   \
    hemidentification:${hi_ver}         \
        --file_conflict ${file_conflict}         \
        --keep_groups                   \
            status train_and_test       \
            PatientStatus LHS  --       \
        loocv                           \
            /input/${input_all}         \
            /results/${oname}           \
            ${outcome}

# LN lesion WAB =====

outcome=aphasia
oname=input-lhs_outcome-aphasia_loocv.csv

docker run \
    -v ${home}/input:/input/            \
    -v ${home}/results/loocv:/results   \
    hemidentification:${hi_ver}         \
        --file_conflict ${file_conflict}         \
        --keep_groups                   \
            status train_and_test       \
            PatientStatus LHS  --       \
        loocv                           \
            /input/${input_all}         \
            /results/${oname}           \
            ${outcome}