import argparse as ap
import pandas as pd
import pickle
import pprint as pp
import itertools
import os

# Get the arguments =====

parser = ap.ArgumentParser()

parser.add_argument("--data", "-d", nargs="+", required=True,
                    help="List of data inputs.")

parser.add_argument("--data_names", "-D", nargs="+", required=False,
                    help="List of labels to apply to the data")

parser.add_argument("--models", "-m", nargs="+", required=True,
                    help="List of models to apply to the data.")

parser.add_argument("--model_names", "-M", nargs="+", required=False,
                    help="List of models to apply to the data.")

parser.add_argument("--output", "-o", default="predicted")

args = parser.parse_args()

# Load the data ====

data = [pd.read_csv(x) for x in args.data]

assert len(args.data) == len(args.data_names), \
    f"The length of the models {len(args.data)} and the number of names" + \
    f"{len(args.data_names)} must match!"

data_and_names = zip(data, args.data_names)

# Load the models ====

models = list()

for i, f in enumerate(args.models):
    with open(f, 'rb') as f:
        models.append(pickle.load(f))

# Name the models ===

if args.model_names is not None:
    assert len(args.models) == len(args.model_names), \
        f"The length of the models ({len(args.models)}) and the number of" + \
        f"model names ({len(args.model_names)}) must match!"

model_names = args.model_names if args.model_names is not None else \
    list(range(0, len(args.models)))

models_and_names = zip(models, model_names)

# Do the preidction ===

# output_names = [f"{args.output}/input-{d}_model-{m}_predicted.csv"
#                 for d, m in itertools.product(args.data_names, model_names)]

for d, m in itertools.product(data_and_names, models_and_names):

    data, data_name = d
    model, model_name = m

    print(f"Running model '{model_name}' on input data '{data_name}'")
    # print(data[data.isna().any(axis=1)])
    # continue

    print("Classes in model:")
    print(model.classes_)

    output_name = f"{args.output}/input-{data_name}_model-{model_name}" + \
        "_predicted.csv"

    # Skip already existing files
    if os.path.exists(output_name):
        print(f"Not overwriting {output_name}!\n")
        continue

    # Extract columns to keep in output and then keep only the columns that
    # the model expects
    features = model.feature_names_in_
    input_id = data.drop(features, axis=1)
    input_data = data[features]

    # Actually run the model =====

    # This is the column of input_data values (i.e. LH/RH)
    fit_pred = model.predict(input_data)

    # This is a dataframe with the LD1, 2, 3 weights depending on which was
    # selected
    fit_xfrm = pd.DataFrame(model.transform(input_data))
    fit_xfrm.columns = ["LD" + str(x + 1) for x in
                        range(len(fit_xfrm.columns))]

    # Spit everything out
    output = pd.concat([input_id,
                        pd.DataFrame(fit_pred, columns=["predicted"]),
                        fit_xfrm],
                        axis=1)

    # print(output)

    os.makedirs("predicted/", exist_ok=True)
    output.to_csv(output_name, index=False)

    print(f"Wrote results to {output_name}.\n")