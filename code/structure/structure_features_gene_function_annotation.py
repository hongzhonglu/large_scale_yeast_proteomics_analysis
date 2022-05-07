import pandas as pnd
from pathlib import Path
from time import time
import os
import pandas as pd
from geometricus import MomentInvariants, SplitType

invariants_kmer = []
invariants_radius = []

# a small test
glucose = "YDL245C or YDL247W or YDR342C or YDR343C or YDR345C or YDR536W or YEL069C or YFL011W or YHR092C or YHR094C or YHR096C or YJL214W or YJL219W or YJR158W or YJR160C or YLR081W or YMR011W or YNR072W or YOL156W or YDR387C"
glucose_list = glucose.split(" or ")
structure_quality_all = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")
structure_select = structure_quality_all[structure_quality_all["gene"].isin(glucose_list)]["id"].tolist()
X_names = structure_select

start_time = time()
for i, key in enumerate(X_names):
    pdb_dir = "/Users/xluhon/Documents/alphafold_pdb/" + key
    if i >= 0 and i % 50 == 0:
        print(f"{i} proteins in {(time() - start_time):.2f} seconds")
    invariants_kmer.append(MomentInvariants.from_pdb_file(pdb_dir, split_type=SplitType.KMER, split_size=16))
    invariants_radius.append(MomentInvariants.from_pdb_file(pdb_dir, split_type=SplitType.RADIUS, split_size=10))



import umap
import numpy as np
import matplotlib.pyplot as plt
from geometricus import GeometricusEmbedding

start_time = time()
kmer_embedder = GeometricusEmbedding.from_invariants(invariants_kmer, resolution=4)
radius_embedder = GeometricusEmbedding.from_invariants(invariants_radius, resolution=4)
print(f"Generated embeddings in {(time() - start_time):.2f} seconds")

all = kmer_embedder.embedding
df = pd.DataFrame(all)

samples = df.values
from scipy.cluster.hierarchy import linkage, dendrogram
mergings = linkage(samples)
dendrogram(mergings, leaf_rotation=0,leaf_font_size=10)
plt.title("Dendrograms")

