#!/bin/bash

home=/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/code

input_name=group-all_parc-brainnetome_normed.pickle

mkdir -p bs_results/

n=10000

docker run \
    -v ${home}/input:/input         \
    -v ${home}/bs_results:/results  \
    hemidentification:latest        \
        null -y                     \
            /input/${input_name}    \
            /results/group-all_outcome-class_bs-${n}.csv \
            class                   \
            -n ${n}

for i in bs_results/*.csv ; do
    zip "${i//csv/zip}" "${i}"
done

rm bs_results/*.csv