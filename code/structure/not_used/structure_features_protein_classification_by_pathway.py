import pandas as pd
from time import time
import umap
import matplotlib.pyplot as plt
from geometricus import GeometricusEmbedding
from geometricus import MomentInvariants, SplitType
from src.protein_process import *

# only focus on interesing subsystems
sce_kegg_pathway = pd.read_excel("data/sce_kegg_pathway.xlsx")
# take EMP, TCA, amino acids synthesis and secondary metabolites synthesis
subsystem = ["Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)", "Biosynthesis of amino acids", "Biosynthesis of secondary metabolites"]

g1 = sce_kegg_pathway[sce_kegg_pathway["name"]==subsystem[0]]["id"].tolist()
X_names = g1
invariants_kmer = []
invariants_radius = []
start_time = time()
for i, key in enumerate(X_names):
    pdb_dir = "/Users/xluhon/Documents/alphafold_pdb/" + key
    if i >= 0 and i % 50 == 0:
        print(f"{i} proteins in {(time() - start_time):.2f} seconds")
    invariants_kmer.append(MomentInvariants.from_pdb_file(pdb_dir, split_type=SplitType.KMER, split_size=16))
    invariants_radius.append(MomentInvariants.from_pdb_file(pdb_dir, split_type=SplitType.RADIUS, split_size=10))


kmer_embedder = GeometricusEmbedding.from_invariants(invariants_kmer, resolution=0.5)
radius_embedder = GeometricusEmbedding.from_invariants(invariants_radius, resolution=4)
print(f"Generated embeddings in {(time() - start_time):.2f} seconds")
df = kmer_embedder.embedding



# plot
reducer = umap.UMAP(metric="cosine", n_components=3)
#reduced = reducer.fit_transform(np.hstack((kmer_embedder.embedding, radius_embedder.embedding)))
reduced = reducer.fit_transform(kmer_embedder.embedding)
indices1 = [i for i,x in enumerate(X_names) if x in g1]
plt.figure()
plt.scatter(reduced[indices1, 0],
            reduced[indices1, 1],
            label="g1", edgecolor="black", linewidth=0.1, alpha=0.8)
plt.legend()

