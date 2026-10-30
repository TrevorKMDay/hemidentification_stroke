#!/bin/bash

mkdir -p results/loocv/

input=input/group-HC_parc-brainnetome_loocv_rZ.pickle
outcome="hemi"

# for j in $(seq -w 63) ; do

#     for k in sv nn lda ; do

#         echo "Hold out: ${j}"
#         oname=results/loocv/HC_group-${j}_method-${k}_loocv_results.csv

#         # Background these jobs for a 3x speed-up
#         python run_model.py  -v                     \
#             --drop_cols     PatientStatus           \
#             --output_name   "${oname}"              \
#             "${outcome}" ${input} "${j}" "${k}"

#         # break 2

#     done

# done

# wait

# for j in $(seq -w 159) ; do

#     for k in sv nn lda ; do

#         echo "Hold out: ${j}"
#         oname=results/loocv/everybody_group-${j}_method-${k}_loocv_results.csv

#         python run_model.py  -v                                     \
#             --drop_cols     PatientStatus                           \
#             --output_name   "${oname}"                              \
#             class input/group-all_parc-brainnetome_loocv_rZ.pickle  \
#             "${j}" "${k}"

#         # break 2

#     done

#     # break

# done

# wait

# for j in $(seq -w 159) ; do

#     for k in sv nn lda ; do

#         echo "Hold out: ${j}"
#         oname=results/loocv/everybody_group-${j}_method-${k}_outcome-patient_hemi-R_loocv_results.csv

#         # Using only RHs, try to predict patient status for cross-validation
#         python run_model.py  -v                                     \
#             --hemis         RH                                      \
#             --output_name   "${oname}"                              \
#             PatientStatus input/group-all_parc-brainnetome_loocv_rZ.pickle  \
#             "${j}" "${k}"

#         # break 2

#     done

#     # break

# done

# wait


input_normed=group-HC_parc-brainnetome_normed.pickle
outcome="hemi"

# for j in $(seq -w 159) ; do

#     k=lda

#     echo "Hold out: ${j}"
#     oname=results/loocv/patient_status/everybodynormed_group-${j}_method-${k}_outcome-patient_hemi-R_loocv_results.csv

#     # Using only RHs, try to predict patient status for cross-validation
#     python run_model.py  -v                                     \
#         --hemis         RH                                      \
#         --output_name   "${oname}"                              \
#         PatientStatus input/group-all_parc-brainnetome_loocv_rZ.pickle  \
#         "${j}" "${k}"

# done

