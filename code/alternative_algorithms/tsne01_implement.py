import pickle as pkl
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

home = Path("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke")

# Load data

input_f = home / "code" / "input" / "group-all_parc-brainnetome_normed.pickle"

with open(input_f, 'rb') as f:
    input = pkl.load(f)

# First step is a PCA before TSNE

cx_matrix = input.filter(regex = r'^cx_', axis=1)
id_matrix = input.filter(regex = r'^(?!cx_)\w+$', axis=1)

pca1 = PCA()
pca_results = pca1.fit(cx_matrix)

# Quick scree plot
fig, ax = plt.subplots()
ax.plot(np.arange(1, len(pca_results.explained_variance_ratio_) + 1), 
        pca_results.explained_variance_ratio_, 
        marker='o')
plt.show()

pca_xform = pca_results.transform(cx_matrix)

# Run the TSNE

tsne1 = TSNE(n_components=2)

tsne_results = pd.DataFrame(tsne1.fit_transform(pca_xform))
tsne_results.columns = [f"tsne_{x}" for x in tsne_results.columns]

plot_tsne = pd.concat([id_matrix, tsne_results], axis=1)

plot_tsne.to_csv(home / "code" / "alternative_algorithms" / "tsne_results_n2.csv")

tsne2 = TSNE(n_components=3)

tsne_results2 = pd.DataFrame(tsne2.fit_transform(pca_xform))
tsne_results2.columns = [f"tsne_{x}" for x in tsne_results2.columns]

plot_tsne2 = pd.concat([id_matrix, tsne_results2], axis=1)

plot_tsne2.to_csv(home / "code" / "alternative_algorithms" / "tsne_results_n3.csv")
