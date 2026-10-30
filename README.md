# Hemisphere Fingerprinting, Stroke

## Code

### Creating inputs

`01-convert_data.R` ingests the HC/LHS data and exports the input files to
`code/input` for use by the modeling programs.

`01a-create_LRavg_data.R` creates fake data that mixes LH/RH data proportionally
to show that increasing amounts of LH information does make the prediction
trend toward an LH label. These live in `code/input/mixed/` and "*mixed*" is
the label for these data.

### Running models

The `02x-` prefix scripts run input data, placing models trained on those data
in `code/results/*.pickle`. Probably should be given like a `_model` suffix for
clarity.
These are created in order of need, and so aren't ordered
in any particular principled way.

- `xvalid`: Runs the three-fold cross-validation.

- `loocv`: Runs the leave-one-out cross-validation models

- `allin`: Create the models trained on everyone, i.e., w/o cross-validation.

- `mixed`: Create the models trained on data that is mixed LH/RH for validation

- `normed`: Train the models on data that has been normed within feature.

### Predicting data

TO DO: Re-do things top-to-bottom with new predict_model.py script.
