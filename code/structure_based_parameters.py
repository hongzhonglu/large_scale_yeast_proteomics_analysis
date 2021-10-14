# Ref to this nice tutorial:
# http://rasbt.github.io/biopandas/tutorials/Working_with_PDB_Structures_in_DataFrames/


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

# volume calculation
"java -jar /Users/xluhon/Documents/ProteinVolume_1.3/ProteinVolume_1.3.jar /Users/xluhon/Documents/alphafold_pdb"

# meeting the following error
# reading hydrogens is turned on, but couldn't find any hydrogens in AF-A0A023PZB3-F1-model_v1!




