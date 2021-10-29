# Based on yeast-GEM 8.5
#

"java -jar /Users/xluhon/Documents/ProteinVolume_1.3/ProteinVolume_1.3.jar /Users/xluhon/Documents/alphafold_pdb"


from Bio.PDB import *
import os    ##for directory
import numpy as np
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from biopandas.pdb import PandasPdb

# input model information
# remove column of "#"
yeast_gem_gene = pd.read_excel("/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xlsx", sheet_name="GENES")
yeast_gem = pd.read_excel("/Users/xluhon/Documents/GitHub/yeast-GEM/model/yeast-GEM.xlsx", sheet_name="RXNS")
# extract rxn and gene
yeast_gem0 = yeast_gem[["ID","GENE ASSOCIATION"]]
yeast_gem0.columns = ["ID", "GPR"]



