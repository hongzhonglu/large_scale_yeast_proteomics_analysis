# Ref to this nice tutorial:
# http://rasbt.github.io/biopandas/tutorials/Working_with_PDB_Structures_in_DataFrames/
# volume calculation
"java -jar /Users/xluhon/Documents/ProteinVolume_1.3/ProteinVolume_1.3.jar /Users/xluhon/Documents/alphafold_pdb"


from Bio.PDB import *
import os    ##for directory
import numpy as np
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from biopandas.pdb import PandasPdb

pdbfile = "/Users/xluhon/Documents/alphafold_pdb/"
all_pdb = os.listdir(pdbfile)
ppdb = PandasPdb()
PITT_score = []
for i in all_pdb:
    print(i)
    pdb_in = pdbfile + i
    ppdb.read_pdb(pdb_in)
    ss = ppdb.df['ATOM']
    mainchain = ss[(ss['atom_name'] == 'CA')]
    bfact_mc_avg = mainchain['b_factor'].mean()
    PITT_score.append(bfact_mc_avg)

data_merge = pd.DataFrame({"id":all_pdb, "score":PITT_score})
data_merge_low_quality = data_merge[data_merge["score"] < 75]
data_merge_high_quality = data_merge[data_merge["score"] >= 75]



# based on quality, estimate which enzyme from Yeast8 need re-modelling
data_merge = pd.read_excel("result/alphafold_quality.xlsx")
data_merge["id_update"] = data_merge["id"].str.replace("AF-", "").str.replace("-F1-model_v1.pdb","")
# get the gene id based on uniprot id
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
# need a function here
def multiMapping (description, item1, item2, dataframe=True, sep=";", removeDuplicates=True):
    """get multiple description of from item1 for item2 based on mapping"""
    #description = w
    #item1 = v
    #item2 = testData
    #used for the list data
    if dataframe:
        description = description.tolist()
        item1 = item1.tolist()
        item2 = item2.tolist()
    else:
        pass
    result = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index0 = [description[index] for index in range(len(item1)) if item1[index] == item2[i]]
            if removeDuplicates:
                index1 = pd.unique(index0).tolist()
            else:
                index1 = index0
            result[i] = sep.join(str(e) for e in index1) #string cat
        else:
            result[i] = None
    return result

data_merge["gene"] = multiMapping(description=id_mapping["GeneName"], item1=id_mapping["Entry"], item2=data_merge["id_update"])
data_merge.to_excel("result/alphafold_quality_with_gene_ID.xlsx")


# input model information
# remove column of "#"
yeast_gem = pd.read_excel("/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xlsx", sheet_name="GENES")

# merge the structure information with model
yeast_gem["structure_id"]= multiMapping(description=data_merge["id"], item1=data_merge["gene"], item2=yeast_gem["NAME"])
yeast_gem["average_score"]= multiMapping(description=data_merge["score"], item1=data_merge["gene"], item2=yeast_gem["NAME"])

# score analysis
yeast_gem["average_score"] = yeast_gem["average_score"].astype(float)

#ax = yeast_gem["average_score"].plot.kde()
ax = yeast_gem["average_score"].plot.hist(bins=12, alpha=0.5)
ax.set_title("pLDDT average score")
ax.set_xlabel("Average score")
ax.set_ylabel("Density")

score_list = yeast_gem["average_score"].tolist()
score_high =[x for x in score_list if x >= 75]

# save the data
yeast_gem0 = yeast_gem[['NAME','SHORT NAME', 'structure_id','average_score']]
yeast_gem0.to_excel("result/yeast_gem_with_structure_id_and_score.xlsx")


# find position of specific protein
# data_merge["id_update"][data_merge["id_update"] == "P38427"].index[0]










