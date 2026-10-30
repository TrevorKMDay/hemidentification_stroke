import sklearn
import pickle as pkl
import pandas as pd

from pathlib import Path

from sklearn import metrics
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.feature_selection import VarianceThreshold, SelectKBest, f_classif

home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/")

data_file = home / "code" / "input" / \
    "group-all_parc-brainnetome_loocv_rZ.pickle"

output_dir = home / "code" / "modeling" / "dim_reduce"

# Load the data ====

print(f"Loading data from {data_file}")
with open(data_file, "rb") as f:
    data = pd.DataFrame(pkl.load(f))

# How to reduce columns

def sel_kbest(k, fts, y, verbose=False):

    print(f"k: {k}")

    selector = SelectKBest(f_classif, k=k)
    x_new = selector.fit_transform(fts, y=y)

    if verbose:

        feature_mask = selector.get_support()
        selected_names = fts.columns[feature_mask]

        names_ordered = [val for _, val in
                        sorted(zip(selector.scores_, fts.columns))]

        scores_ordered = sorted(selector.scores_)

        for name, score in zip(names_ordered, scores_ordered):
            print(f"Feature: {name:<20} Score: {score:.2f}")

        print(f"\nChosen Features: {list(selected_names)}")

    model = clf.fit(x_new, y)
    acc = model.score(x_new, y)

    print([k, acc])

    print(model)

    return(acc)

# Five-fold

for test_group in ["A", "B", "C", "D", "E"]:

    test1 = data[data["group"] == test_group]
    train1 = data[data["group"] != test_group]

    id = train1[["sub", "hemi"]]
    pred_col = train1["hemi"].tolist()

    print(f"Testing on {test_group}, training on the rest ...")

    # Training data
    features = train1.drop(["sub", "gender", "age", "age_group", "hemi",
                            "handedness", "handedness2", "class", "class2",
                            "group", "EHI"],
                           axis=1, errors="ignore")


    # Set up the LDA
    clf = LinearDiscriminantAnalysis()

    m1 = clf.fit(features, pred_col)

    # Subset the feature columns
    test_x = test1[m1.feature_names_in_]
    test_y = test1[["sub", "hemi"]]

    test_result = m1.predict(test_x)
    # mcc = metrics.matthews_corrcoef(test_y, test_result)

    print("Starting k-best selection")

    selector1 = SelectKBest(f_classif, k=100)
    selector1.fit(features, y=pred_col)
    scores_df = pd.DataFrame({"feature": features.columns,
                              "score": selector1.scores_})

    scores_df.to_csv(output_dir / f"group-{test_group}_selectorscores.csv")

    for k in [2**x for x in range(10)]:

        print(f"\nk={k}")

        selector = SelectKBest(f_classif, k=k)
        x_new = selector.fit_transform(features, y=pred_col)

        # What columns were selected
        x_cols = [col for col, x in
                  zip(features.columns, selector.get_support())
                  if x]

        # Create new test data
        test_new = test1[x_cols]

        m2 = LinearDiscriminantAnalysis()
        m2 = clf.fit(x_new, pred_col)

        acc = m2.score(x_new, pred_col)
        print(f"Accuracy: {round(acc, 3)}")

        # Do the out-of-group prediction
        test_y["predicted"] = m2.predict(test_new)

        test_y.to_csv(output_dir / f"group-{test_group}_k-{k}_results.csv")
        # print(test_y)

        # break