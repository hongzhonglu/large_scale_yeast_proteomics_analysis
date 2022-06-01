import pandas as pnd
from pathlib import Path
from time import time
import os
import pandas as pd
import umap
import numpy as np
import matplotlib.pyplot as plt
from geometricus import GeometricusEmbedding
from geometricus import MomentInvariants, SplitType


# only focus on interesing subsystems
sce_kegg_pathway = pd.read_excel("data/sce_kegg_pathway.xlsx")
# take EMP, TCA, amino acids synthesis and secondary metabolites synthesis
subsystem = ["Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)", "Biosynthesis of amino acids", "Biosynthesis of secondary metabolites"]

g1 = sce_kegg_pathway[sce_kegg_pathway["name"]==subsystem[0]]["id"].tolist()
g20 = sce_kegg_pathway[sce_kegg_pathway["name"].isin(subsystem[0:2])]["id"].tolist()
g2 = list(set(g20)-set(g1))
g30 = sce_kegg_pathway[sce_kegg_pathway["name"].isin(subsystem[0:3])]["id"].tolist()
g3 = list(set(g30)-set(g1)-set(g2))
g40 = sce_kegg_pathway[sce_kegg_pathway["name"].isin(subsystem[0:4])]["id"].tolist()
g4 = list(set(g40)-set(g1)-set(g2)-set(g3))


X_names = g1 + g2 + g3 + g4
#X_names = os.listdir("/Users/xluhon/Documents/alphafold_test/")
X_names = [x for x in X_names if x !=".DS_Store"]
invariants_kmer = []
invariants_radius = []
start_time = time()
for i, key in enumerate(X_names):
    pdb_dir = "/Users/xluhon/Documents/alphafold_pdb/" + key
    if i >= 0 and i % 50 == 0:
        print(f"{i} proteins in {(time() - start_time):.2f} seconds")
    invariants_kmer.append(MomentInvariants.from_pdb_file(pdb_dir, split_type=SplitType.KMER, split_size=16))
    invariants_radius.append(MomentInvariants.from_pdb_file(pdb_dir, split_type=SplitType.RADIUS, split_size=10))


kmer_embedder = GeometricusEmbedding.from_invariants(invariants_kmer, resolution=4)
radius_embedder = GeometricusEmbedding.from_invariants(invariants_radius, resolution=4)
print(f"Generated embeddings in {(time() - start_time):.2f} seconds")


# data dimension reduction
reducer = umap.UMAP(metric="cosine", n_components=2)
#reduced = reducer.fit_transform(np.hstack((kmer_embedder.embedding, radius_embedder.embedding)))
reduced = reducer.fit_transform(kmer_embedder.embedding)
indices1 = [i for i,x in enumerate(X_names) if x in g1]
indices2 = [i for i,x in enumerate(X_names) if x in g2]
indices3 = [i for i,x in enumerate(X_names) if x in g3]
indices4 = [i for i,x in enumerate(X_names) if x in g4]


# plot
plt.figure()
plt.scatter(reduced[indices1, 0],
            reduced[indices1, 1],
            label="g1", edgecolor="black", linewidth=0.1, alpha=0.8)
plt.scatter(reduced[indices2, 0],
            reduced[indices2, 1],
            label="g2", edgecolor="red", linewidth=0.1, alpha=0.8)
plt.scatter(reduced[indices3, 0],
            reduced[indices3, 1],
            label="g3", edgecolor="green", linewidth=0.1, alpha=0.8)
plt.scatter(reduced[indices4, 0],
            reduced[indices4, 1],
            label="g4", edgecolor="blue", linewidth=0.1, alpha=0.8)
plt.legend()











g1 = sce_kegg_pathway[sce_kegg_pathway["name"]==subsystem[0]]["id"].tolist()
X_names = g1
#X_names = os.listdir("/Users/xluhon/Documents/alphafold_test/")
X_names = [x for x in X_names if x !=".DS_Store"]
invariants_kmer = []
invariants_radius = []
start_time = time()
for i, key in enumerate(X_names):
    pdb_dir = "/Users/xluhon/Documents/alphafold_pdb/" + key
    if i >= 0 and i % 50 == 0:
        print(f"{i} proteins in {(time() - start_time):.2f} seconds")
    invariants_kmer.append(MomentInvariants.from_pdb_file(pdb_dir, split_type=SplitType.KMER, split_size=16))
    invariants_radius.append(MomentInvariants.from_pdb_file(pdb_dir, split_type=SplitType.RADIUS, split_size=10))


kmer_embedder = GeometricusEmbedding.from_invariants(invariants_kmer, resolution=4)
radius_embedder = GeometricusEmbedding.from_invariants(invariants_radius, resolution=4)
print(f"Generated embeddings in {(time() - start_time):.2f} seconds")


# data dimension reduction
reducer = umap.UMAP(metric="cosine", n_components=2)
#reduced = reducer.fit_transform(np.hstack((kmer_embedder.embedding, radius_embedder.embedding)))
reduced = reducer.fit_transform(kmer_embedder.embedding)
indices1 = [i for i,x in enumerate(X_names) if x in g1]


# plot
plt.figure()
plt.scatter(reduced[indices1, 0],
            reduced[indices1, 1],
            label="g1", edgecolor="black", linewidth=0.1, alpha=0.8)
plt.legend()
