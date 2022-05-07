import pandas as pnd
from pathlib import Path
from time import time
import os
import pandas as pd



"""
url = "https://raw.githubusercontent.com/TurtleTools/geometricus/master/example_data/MAPK_KLIFS.tsv"
mapk_df = pnd.read_csv(url, sep="\t")

print(",".join(mapk_df["PDB"]))

mapk_pdb_id_to_class = {}
for pdb_id, chain, class_name in list(zip(mapk_df["PDB"], mapk_df["CHAIN"], mapk_df["CLASS"])):
    mapk_pdb_id_to_class[(pdb_id, chain)] = class_name
len(mapk_pdb_id_to_class)

X_names = list(mapk_pdb_id_to_class.keys())
class_mapping = {"JNK": 0, "Erk": 1, "p38": 2}
y = [class_mapping[mapk_pdb_id_to_class[k]] for k in X_names]
"""


"""
import prody as pd

start_time = time()
pdbs = []
for i, (pdb_id, chain) in enumerate(X_names):
    if i > 0 and i % 50 == 0:
        print(f"{i} proteins fetched in {(time() - start_time):.2f} seconds")
    pdbs.append(pd.parsePDB(pdb_id, chain=chain))
"""


from geometricus import MomentInvariants, SplitType
X_names = os.listdir("/Users/xluhon/Documents/alphafold_pdb/")
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



import umap
import numpy as np
import matplotlib.pyplot as plt
from geometricus import GeometricusEmbedding

start_time = time()
kmer_embedder = GeometricusEmbedding.from_invariants(invariants_kmer, resolution=4)
radius_embedder = GeometricusEmbedding.from_invariants(invariants_radius, resolution=4)
print(f"Generated embeddings in {(time() - start_time):.2f} seconds")

"""
reducer = umap.UMAP(metric="cosine", n_components=4)
reduced = reducer.fit_transform(np.hstack((kmer_embedder.embedding, radius_embedder.embedding)))

class_names = ["test"]
colors = ["red"]
indices = list(range(0,len(X_names)))

plt.figure()
plt.scatter(reduced[indices, 0],
            reduced[indices, 1],
            label=class_names[0], edgecolor="black", linewidth=0.1, alpha=0.8)
plt.axis("off")
plt.legend()
"""

all = kmer_embedder.embedding
df = pd.DataFrame(all)

samples = df.values
from scipy.cluster.hierarchy import linkage, dendrogram
mergings = linkage(samples)
dendrogram(mergings, leaf_rotation=0,leaf_font_size=10)
plt.title("Dendrograms")


"""
all = kmer_embedder.embedding
df = pd.DataFrame(all)
df.index = X_names
df.to_csv("data/structure_feature_kmer_4.txt",  sep="\t")

all2 = radius_embedder.embedding
df2 = pd.DataFrame(all2)
df2.index = X_names
df2.to_csv("data/structure_feature_radius_4.txt", sep="\t")
"""


# read the orginal datasets
d1 = pd.read_csv("data/structure_feature_kmer_4.txt", sep="\t")









# read the orginal datasets
d2 = pd.read_csv("data/structure_feature_radius_4.txt", sep="\t")
new_columns = ["Unnamed: 0"] + ["para_" + str(x) for x in range(0,d2.shape[1]-1)]
d2.columns = new_columns

df = pd.concat([d1, d2], axis=1)
df.index = df.iloc[:,0]
all_columns = list(df.columns)
all_columns = [x for x in all_columns if x !="Unnamed: 0"]
df = df[all_columns]



## for classification of core and variable genes
pro_type = pd.read_excel("data/sce_protein_with_core_variable_type.xlsx")
pro_type = pro_type[["id", "gene_type"]]



# machine learning
from sklearn.model_selection import train_test_split


class1 = pro_type[pro_type["gene_type"]=="core_gene"]
class2 = pro_type[pro_type["gene_type"]=="Variable"]

test_under = pd.concat([class1, class2], axis=0)

X_names = test_under["id"].tolist()
y = test_under["gene_type"].tolist()



X_train_names, X_test_names, y_train, y_test = train_test_split(X_names, y, test_size=0.30)


df0 = df.iloc[:,1:]


X_train = df0[df0.index.isin(X_train_names)]
X_test = df0[df0.index.isin(X_test_names)]

# get the columns with all elements equal to zero in X_train
all_columns = list(X_train.columns)
all_columns = [x for x in all_columns if x !="Unnamed: 0"]

columns_used = []
for x in all_columns:
    s0 = X_train[x].tolist()
    s1 = sum(s0)
    if s1 >= 1:
        columns_used.append(x)

X_train1 = X_train[columns_used]
X_test1 = X_test[columns_used]


X_train1.to_numpy()
X_test1.to_numpy()


"""
# svc too small accuracy
from sklearn.svm import SVC
from sklearn.metrics import classification_report
clf = SVC(class_weight='balanced', probability=True)
clf.fit(X_train1, y_train)
y_pred = clf.predict(X_test1)
class_names = ["core_gene","Variable"]
print(classification_report(y_test, y_pred, class_names))
"""


# random forest
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
clf = RandomForestClassifier(max_depth=2, random_state=20)
clf.fit(X_train1, y_train)
y_pred = clf.predict(X_test1)
class_names = ["core_gene","Variable"]
print(classification_report(y_test, y_pred, class_names))



# decision tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report
clf = DecisionTreeClassifier(random_state=42, max_depth=3) # 42, 3
clf.fit(X_train1, y_train)
y_pred = clf.predict(X_test1)
class_names = ["core_gene","Variable"]
print(classification_report(y_test, y_pred, class_names))




"""
# TF prediction
pro_type = pd.read_excel("data/sce_protein_with_TF_classification.xlsx")
pro_type = pro_type[["id", "TF", "volume_per_kda"]]
pro_type.columns = ["id", "gene_type", "volume_per_kda"]


from sklearn.model_selection import train_test_split
from src.protein_process import *

X_names = pro_type["id"].tolist()
y = pro_type["gene_type"].tolist()
X_train_names, X_test_names, y_train, y_test = train_test_split(X_names, y, test_size=0.30)

df.index = df.iloc[:,0]

# add new information
df["volume_per_kda"] = singleMapping(pro_type["volume_per_kda"],pro_type["id"],df.iloc[:,0])
df0 = df.iloc[:,1:]


X_train = df0[df0.index.isin(X_train_names)]
X_test = df0[df0.index.isin(X_test_names)]

# get the columns with all elements equal to zero in X_train
all_columns = list(X_train.columns)
all_columns = [x for x in all_columns if x !="volume_per_kda"]
columns_used = []
for x in all_columns:
    s0 = X_train[x].tolist()
    s1 = sum(s0)
    if s1 >= 1:
        columns_used.append(x)

X_train1 = X_train[columns_used + ["volume_per_kda"]]
X_test1 = X_test[columns_used + ["volume_per_kda"]]


X_train1.to_numpy()
X_test1.to_numpy()


# decision tree
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import classification_report
clf = DecisionTreeClassifier(random_state=42, max_depth=3) # 42, 3
clf.fit(X_train1, y_train)
y_pred = clf.predict(X_test1)
class_names = ["No","Yes"]
print(classification_report(y_test, y_pred, class_names))
"""