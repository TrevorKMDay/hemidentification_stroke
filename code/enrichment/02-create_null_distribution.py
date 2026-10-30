import pickle as pkl

from run_model_kfold import create_null_distribution
from pathlib import Path

home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/",
            "enrichment")

hc_file = home / ".." / "input"  / "group-all_parc-brainnetome_normed.pickle"

with open(hc_file, "rb") as f:
    hc0 = pkl.load(f)

input = ("hclhs", hc0)

# Get columns that don't start with cx_
non_pred_cols = [x for x in hc0.columns if "cx_" not in x]

print([x for x in hc0.columns if "TE" in x])

# create_null_distribution(input,
#                          id_cols=non_pred_cols,
#                          cols_to_shuffle=["class"],
#                          output_dir=home / "bs_results",
#                          n=5000)

