#!/bin/bash

input=group-all_parc-brainnetome_normed.pickle
mkdir -p threefold/

# Train models ====

docker run \
    -v "$(pwd)"/input:/input                            \
    -v "$(pwd)"/threefold:/output                       \
    hemidentification:latest                            \
    train                                               \
        /input/group-all_parc-brainnetome_normed.pickle \
        /output/folds-12_model.pickle \
        class    \
        --keep_groups fold 1 2

docker run \
    -v "$(pwd)"/input:/input                            \
    -v "$(pwd)"/threefold:/output                       \
    hemidentification:latest                            \
    train                                               \
        /input/group-all_parc-brainnetome_normed.pickle \
        /output/folds-13_model.pickle \
        class    \
        --keep_groups fold 1 3

docker run \
    -v "$(pwd)"/input:/input                            \
    -v "$(pwd)"/threefold:/output                       \
    hemidentification:latest                            \
    train                                               \
        /input/group-all_parc-brainnetome_normed.pickle \
        /output/folds-23_model.pickle \
        class    \
        --keep_groups fold 2 3

# Test models ====

docker run \
    -v "$(pwd)"/input:/input                            \
    -v "$(pwd)"/threefold:/output                       \
    hemidentification:latest                            \
    test                                               \
        /output/folds-12_model.pickle                   \
        /input/group-all_parc-brainnetome_normed.pickle \
        /output/fold-12.csv \
        --keep_groups fold 1 2 3

docker run \
    -v "$(pwd)"/input:/input                            \
    -v "$(pwd)"/threefold:/output                       \
    hemidentification:latest                            \
    test                                               \
        /output/folds-13_model.pickle                   \
        /input/group-all_parc-brainnetome_normed.pickle \
        /output/fold-13.csv \
        --keep_groups fold 1 2 3

docker run \
    -v "$(pwd)"/input:/input                            \
    -v "$(pwd)"/threefold:/output                       \
    hemidentification:latest                            \
    test                                               \
        /output/folds-23_model.pickle                   \
        /input/group-all_parc-brainnetome_normed.pickle \
        /output/fold-23.csv \
        --keep_groups fold 1 2 3