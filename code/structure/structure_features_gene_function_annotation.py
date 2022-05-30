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

structure_select0 = structure_quality_all[structure_quality_all["gene"].isin(glucose_list)]
structure_select0 = structure_select0[["id","gene"]]
structure_select0["number"] = list(range(0,20))
structure_select0.to_excel("result/glucose_transportor_classification.xlsx")
structure_select = structure_select0["id"].tolist()

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
dendrogram(mergings, leaf_rotation=0, leaf_font_size=10)
plt.title("Dendrograms")




# PCA analysis
from sklearn.decomposition import PCA
pca_test = PCA(n_components=2)
principalComponents_breast = pca_test.fit_transform(df)
principal_breast_Df = pd.DataFrame(data = principalComponents_breast, columns = ['x1', 'x2'])


import seaborn as sns
pd_null = principal_breast_Df
# plot
x0='x1'
y0='x2'
plt.figure()
sns.set_style('darkgrid')
plt.scatter(principal_breast_Df["x1"], principal_breast_Df["x2"])
plt.xlabel('principal component 1', fontsize=15)
plt.ylabel('principal component 2', fontsize=15)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
# annotate the dataset points
xs=pd_null[x0].tolist()
ys=pd_null[y0].tolist()
tlab=list(pd_null.index)
for x, y, lab in zip(xs, ys, tlab):
    plt.annotate(lab,  # this is the text (put lab here to use tlab as string)
                 (x, y),  # this is the point to label
                 textcoords="offset points",  # how to position the text
                 xytext=(0, 10),  # distance from text to points (x,y)
                 ha='center',
                 fontsize=5)
plt.show()









# plot
reducer = umap.UMAP(metric="cosine", n_components=2)
#reduced = reducer.fit_transform(np.hstack((kmer_embedder.embedding, radius_embedder.embedding)))
reduced = reducer.fit_transform(kmer_embedder.embedding)
indices1 = [i for i,x in enumerate(X_names)]
plt.scatter(reduced[indices1, 0],
            reduced[indices1, 1],
            label="g1", edgecolor="black", linewidth=0.1, alpha=0.8)



