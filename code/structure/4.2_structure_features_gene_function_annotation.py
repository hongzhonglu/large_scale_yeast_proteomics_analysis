from time import time
import pandas as pd
from geometricus import MomentInvariants, SplitType
import matplotlib.pyplot as plt
from geometricus import GeometricusEmbedding
from sklearn.decomposition import PCA
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram




# a small test
# check whether the structure feature could help to characterize the functions of protein homologs
glucose = "YDL245C or YDL247W or YDR342C or YDR343C or YDR345C or YDR536W or YEL069C or YFL011W or YHR092C or YHR094C or YHR096C or YJL214W or YJL219W or YJR158W or YJR160C or YLR081W or YMR011W or YNR072W or YOL156W or YDR387C"
glucose_list = glucose.split(" or ")
structure_quality_all = pd.read_excel("result/alphafold_quality_with_gene_ID.xlsx")

structure_select0 = structure_quality_all[structure_quality_all["gene"].isin(glucose_list)]
structure_select0 = structure_select0[["id","gene"]]
structure_select0["number"] = list(range(0,20))
structure_select0.to_excel("result/glucose_transportor_classification.xlsx")
structure_select = structure_select0["id"].tolist()
X_names = structure_select


# feature extraction
invariants_kmer = []
invariants_radius = []
start_time = time()
for i, key in enumerate(X_names):
    pdb_dir = "/Users/xluhon/Documents/alphafold_pdb/" + key
    if i >= 0 and i % 50 == 0:
        print(f"{i} proteins in {(time() - start_time):.2f} seconds")
    invariants_kmer.append(MomentInvariants.from_pdb_file(pdb_dir, split_type=SplitType.KMER, split_size=16)) # split size = 16
    invariants_radius.append(MomentInvariants.from_pdb_file(pdb_dir, split_type=SplitType.RADIUS, split_size=10))


kmer_embedder = GeometricusEmbedding.from_invariants(invariants_kmer, resolution=4)
#radius_embedder = GeometricusEmbedding.from_invariants(invariants_radius, resolution=4)
all = kmer_embedder.embedding
#all = radius_embedder.embedding
df = pd.DataFrame(all)
samples = df.values


# Multilevel cluster analysis
mergings = linkage(samples)
dendrogram(mergings, leaf_rotation=0, leaf_font_size=10)
plt.title("Dendrograms")




# PCA analysis
pca_test = PCA(n_components=2)
principalComponents = pca_test.fit_transform(df)
principal_Df = pd.DataFrame(data=principalComponents, columns=['x1', 'x2'])
pd_null = principal_Df
# plot
x0='x1'
y0='x2'
plt.figure()
sns.set_style('darkgrid')
plt.scatter(principal_Df["x1"], principal_Df["x2"])
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










# get the important features
from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
np.random.seed(0)


# PCA feature analysis
pca_test = PCA(n_components=2)
model = PCA(n_components=2).fit(all)
X_pc = model.transform(all)

# number of components
n_pcs= model.components_.shape[0]

# get the index of the most important feature on EACH component
# LIST COMPREHENSION HERE
most_important = [np.abs(model.components_[i]).argmax() for i in range(n_pcs)] # a function

initial_feature_names = [i for i in range(0,model.components_.shape[1])]

# get the names
most_important_names = [initial_feature_names[most_important[i]] for i in range(n_pcs)]

# LIST COMPREHENSION HERE AGAIN
dic = {'PC{}'.format(i+1): most_important_names[i] for i in range(n_pcs)}

# build the dataframe
df = pd.DataFrame(dic.items())
df.columns = ["component","index"]


# mapping to the shapemer_keys
index_detail = df["index"][df["component"]=="PC1"][0]
shapemer = kmer_embedder.shapemer_keys[index_detail] # 242 is the feature IDs
residue_indices_train = kmer_embedder.map_shapemer_to_residues(shapemer)
print("Shape-mer:", shapemer, "Number of proteins with shape-mer:", len(residue_indices_train))
print()
print("Residue indices per protein (for 10 proteins):")
for i, key in enumerate(residue_indices_train):
    print(key, residue_indices_train[key])

