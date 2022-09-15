from time import time
import os
import pandas as pd
import umap
import matplotlib.pyplot as plt
from geometricus import GeometricusEmbedding
from geometricus import MomentInvariants, SplitType
from src.protein_process import *
import seaborn as sns

# only focus on interesting subsystems
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


# check the subsystem classification based on SGD database
sce_family = pd.read_excel("data/sce_protein_family.xlsx")
structure_info = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")
sce_family["id"] = singleMapping(structure_info["id"], structure_info["gene"], sce_family["gene"])
sce_family =sce_family[sce_family['Cross-reference (SUPFAM)'].notna()]
sce_family0 = sce_family[sce_family["Cross-reference (SUPFAM)"].str.contains("SSF56112")]
g1_family = sce_family0["id"].tolist()


# X_names = g1 + g2 + g3 + g4
X_names = os.listdir("/Users/xluhon/Documents/alphafold_pdb/")
X_names = [x for x in X_names if x !=".DS_Store"]
X_names = [x for x in X_names if "model" in x]

# remove outlier datapoints
outlier1 = pd.read_excel("result/cluster_analysis_all_proteins_outlier_group1.xlsx")
outlier2 = pd.read_excel("result/cluster_analysis_all_proteins_outlier_group2.xlsx")
gene_remove = outlier1["id"].tolist() + outlier2["id"].tolist()

X_names = [x for x in X_names if x not in gene_remove]

# score >= 70
structure_info = structure_info[structure_info["score"] >= 70]
X_names = [x for x in X_names if x in structure_info["id"].tolist()]

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
radius_embedder = GeometricusEmbedding.from_invariants(invariants_radius, resolution=6)
print(f"Generated embeddings in {(time() - start_time):.2f} seconds")


# data dimension reduction
reducer = umap.UMAP(metric="cosine", n_components=2)
#reduced = reducer.fit_transform(kmer_embedder.embedding)
reduced = reducer.fit_transform(np.hstack((kmer_embedder.embedding, radius_embedder.embedding)))



# plot as a whole
class_names = ["test"]
colors = ["red"]

# plot
# classification based on families
indices = list(range(0,len(X_names)))
indices1 = [i for i,x in enumerate(X_names) if x in g1_family]
plt.figure()
plt.scatter(reduced[indices, 0],
            reduced[indices, 1],
            label=class_names[0], facecolors='none', edgecolor="black", linewidth=0.1, alpha=1)
plt.scatter(reduced[indices1, 0],
            reduced[indices1, 1],
            label="g1", edgecolor="red", linewidth=0.1, alpha=0.8)
plt.xlabel('Principal component 1', fontsize=15)
plt.ylabel('Principal component 2', fontsize=15)
plt.savefig("result/figure_cluster3_map_gene_family_info.pdf", bbox_inches='tight')



# plot
# based on proteins with or without function annotation
s1 = open('/Users/xluhon/Documents/data_for_structure_align/protein_seq_for_alphafold/UP000002311_559292.fasta').readlines()
s2 = [x for x in s1 if '>' in x]

# uncharacterized protein
s2_c = [x for x in s2 if 'Uncharacterized' in x]
s2_other = list(set(s2)-set(s2_c))

# remove the detailed function annotation
s2_c = [x.split(' ')[0] for x in s2_c]
s2_other = [x.split(' ')[0] for x in s2_other]

# check the name as the protein stucture ID
s2_c = [x.split('|')[1] for x in s2_c]
s2_other = [x.split('|')[1] for x in s2_other]

s2_c = ['AF-' + x + '-F1-model_v1.pdb' for x in s2_c]
s2_other = ['AF-' + x + '-F1-model_v1.pdb' for x in s2_other]



indices = list(range(0,len(X_names)))
indices1 = [i for i,x in enumerate(X_names) if x in s2_c]
plt.figure()
sns.set_style("white")
plt.scatter(reduced[indices, 0],
            reduced[indices, 1],
            label=class_names[0], facecolors='none', edgecolor="black", linewidth=0.1, alpha=1)
plt.scatter(reduced[indices1, 0],
            reduced[indices1, 1],
            label="g1", edgecolor="red", linewidth=0.1, alpha=0.8)
plt.xlabel('Principal component 1', fontsize=15)
plt.ylabel('Principal component 2', fontsize=15)
plt.savefig("result/figure_cluster3_map_uncharacteried_protein.pdf", bbox_inches='tight')






# plot
# classification based on subsystem definition
indices1 = [i for i,x in enumerate(X_names) if x in g1]
indices4 = [i for i,x in enumerate(X_names) if x in g4]
plt.figure()
plt.scatter(reduced[indices, 0],
            reduced[indices, 1],
            label=class_names[0], facecolors='none', edgecolor="black", linewidth=0.1, alpha=1)
plt.scatter(reduced[indices1, 0],
            reduced[indices1, 1],
            label="g1", edgecolor="red", linewidth=0.1, alpha=0.8)
plt.scatter(reduced[indices4, 0],
            reduced[indices4, 1],
            label="g4", edgecolor="blue", linewidth=0.1, alpha=0.8)
plt.xlabel('Principal component 1', fontsize=15)
plt.ylabel('Principal component 2', fontsize=15)
plt.savefig("result/figure_cluster3_map_classification.pdf", bbox_inches='tight')




# data output
data_analysis = pd.DataFrame(reduced)
data_analysis["id"] = X_names
data_analysis.columns = ["x1", "x2", "id"]
data_analysis.to_excel("result/feature_for_all_protein_structure.xlsx")




# initial data analysis for some clusters
data_check = data_analysis[data_analysis["x2"] >=2]
# input the detailed ID information
structure_info = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")
data_check["score"] = singleMapping(structure_info["score"],structure_info["id"],data_check["id"])
data_check["gene"] = singleMapping(structure_info["gene"],structure_info["id"],data_check["id"])
gene01= ",".join(data_check["gene"].to_list())
print(gene01)
data_check.to_excel("result/cluster_analysis_all_proteins_outlier_group1.xlsx")


data_check = data_analysis[data_analysis["x1"] >=10]
# input the detailed ID information
structure_info = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")
data_check["score"] = singleMapping(structure_info["score"],structure_info["id"],data_check["id"])
data_check["gene"] = singleMapping(structure_info["gene"],structure_info["id"],data_check["id"])
gene01= ",".join(data_check["gene"].to_list())
print(gene01)
data_check.to_excel("result/cluster_analysis_all_proteins_outlier_group2.xlsx")


data_check = data_analysis[data_analysis["x1"] >= 7]
data_check = data_check[data_check["x1"] <= 8]
data_check = data_check[data_check["x2"] >= 7]

# input the detailed ID information
structure_info = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")
data_check["score"] = singleMapping(structure_info["score"],structure_info["id"],data_check["id"])
data_check["gene"] = singleMapping(structure_info["gene"],structure_info["id"],data_check["id"])
gene01= ",".join(data_check["gene"].to_list())
print(gene01)
data_check.to_excel("result/cluster_analysis.xlsx")
