# LH Stroke Analysis

I began by running the un-normed analyses. However, we revised the
hemiconnectome paper to primarily rely on the normed data, so most of the
analyses here come from "normed data." Those analyses were moved to
`old_code_unnormed/`, hence the strange numbering and lettering here.

- **HC**: healthy controls
- **LHS**: left-hemisphere stroke participants

"Class" refers to the three-way class: `HC LH`, `HC RH`, and `LHS RH`.

## Data preparation

- Two individuals without a sex value are dropped for fold balancing.
- One scan from a pilot participant is removed

## Preliminary analyses

The preliminary analyses use LOOCV to establish that the results are robust,
see `02a-run_normed_loocv.sh`.
The statistical analysis after the prediction is done in R,
`05a-preliminary_analyses_normed.R`.

### Hemisphere identification

1. **HCs only**: Uses data `group-all_parc-brainnetome_normed.pickle` and
    creates `input-HC_outcome-hemi_loocv.csv`. The LHS are dropped in the
    Docker call.

2. **All ppts**: Uses data `group-all_parc-brainnetome_normed.pickle` and
    creates `input-all_outcome-hemi_loocv.csv`.

### Patient Status

Patient status is done three ways. All predictions use the full dataset:
`group-all_parc-brainnetome_normed.pickle`

1. **Class**: Predicts membership in HC LH, HC RH, or LHS RH. Creates
    `input-all_outcome-class_loocv.csv`.

1. **PatientStatus**: Predicts membership in {HC, LHS}, regardless of
    hemisphere. Creates  `input-all_outcome-patient_loocv.csv`. This is
    not presently included in the paper.

1. **PatientStatus from RHs only**: Predicts membership in {HC, LHS}, but only
    uses RHs to create the file:
    `input-allRHonly_outcome-patient_loocv.csv`.

### Aphasia characteristics

The aphasia characteristics are log lesion volume, WAB, and residualized WAB
(against lesion) size. These use the full dataset cited above and create the
files: `input-lhs_outcome-{lLV,WAB,rWAB}g_loocv.csv` (`g` for grouped).

## Primary analyses

In these models, models are trained on all participants and then the data is
transformed into LDs using those weights. Training is done in
`02c-run_normed_all.sh` and creates the following files:

    # All participants:
    results/input-all_outcome-class_model.pickle
    results/input-all_outcome-hemi_model.pickle
    results/input-all_outcome-PatientStatus_model.pickle

    # Healthy controls only:
    results/input-hc_outcome-hemi_model.pickle

Next, the transformation is done in `04-predict_normed.sh` and those results
are saved to `predicted/`. These have the same format, but the field
`predictedon-` describes which set of data the weights were applied to.

    input-all_outcome-class_predictedon-all_results.csv
    input-all_outcome-hemi_predictedon-all_results.csv
    input-all_outcome-PatientStatus_predictedon-all_results.csv
    input-hc_outcome-hemi_predictedon-all_results.csv