#!/bin/bash

mkdir -p results/

# input=input/group-HC_parc-brainnetome_rZ.pickle
# outcome="hemi"

# for j in A B C ; do

#     for k in sv nn lda ; do

#         echo "Hold out: ${j}"

#         python run_model.py                                                  \
#             --drop_cols     PatientStatus                                    \
#             --output_name   "results/HC_group-${j}_method-${k}_results.csv"  \
#             "${outcome}" ${input} ${j} "${k}"

#     done

# done

input=input/group-all_parc-brainnetome_loocv_rZ.pickle
outcome="PatientStatus"

for j in $(seq -w 1 159) ; do

    for k in sv nn lda ; do

        echo "Hold out: ${j}"

        python run_model.py   -v                                                                    \
            --output_name   "results/loocv/HC_group-${j}_method-${k}_outcome-patient_results.csv"   \
            "${outcome}" ${input} "${j}" "${k}"

        break 2

    done

done