# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

import os
import pandas as pd

protein_list = pd.read_excel('/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/data/protein_list.xlsx')

infile = '/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/Simulation using ITASSER/'
# choose the proteins with pdb
protein_list0 = protein_list[protein_list['Simulated_ITASSER']=='YES']

# then we parse the structure parameters from each pdb file
protein_list0['model_parameter'] = [None]*len(protein_list0['geneID'])
for i, x in protein_list0.iterrows():
    print(i, x['geneID'])
    s0 = x['geneID']
    parameter_file = infile + s0 +'/cscore.txt'
    parameter = open(parameter_file,'r').readlines()
    index1 = [i for i,x in enumerate(parameter) if '#model' in x or 'Name' in x]
    index2 = [i for i,x in enumerate(parameter) if 'Model1' in x]
    para_id = parameter[index1[0]].strip('\n').split('  ')
    para_inf = parameter[index2[0]].strip('\n').split('  ')
    para_id0 = [x.strip(' ') for x in para_id if x not in ['',' ']]
    para_id0 = [x.split(':') for x in para_id0]
    para_id0 = sum(para_id0, [])
    para_inf0 = [x.strip(' ') for x in para_inf if x not in ['',' ']]
    para_dict = dict()
    for m, n in zip(para_id0,para_inf0):
        para_dict[m] = n
    protein_list0['model_parameter'][i] = para_dict

# next we extract the mode1.pdb from each file and put it into a new file
# then we parse the structure parameters from each pdb file
import shutil

protein_list0['pdbid'] = [None]*len(protein_list0['geneID'])
for i, x in protein_list0.iterrows():
    print(i, x['geneID'])
    s0 = x['geneID']
    pdb_file = infile + s0 +'/model1.pdb'
    second_file = '/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/pdb_itasser_structure_mapping/' + s0 + '@model1.pdb'
    shutil.copy(pdb_file, second_file)
    protein_list0['pdbid'][i] = s0 + '@model1'

# then we check the residues number in each pdb file and compared it with the original protein length
from Bio.PDB.PDBParser import PDBParser
protein_list0['protein_length'] = [None]*len(protein_list0['geneID'])
for i, x in protein_list0.iterrows():
    print(i, x['pdbid'])
    s0 = x['pdbid']
    in_file = '/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/pdb_itasser_structure_mapping/' + s0 + '.pdb'
    p = PDBParser()
    structure = p.get_structure(s0, in_file)
    model = structure[0]
    chainID0 = []
    for chain in model:
        chainID0.append(chain.get_id())
    one_chain = model[chainID0[0]]
    protein_list0['protein_length'][i] = len(one_chain)

# after check we can find that the in the itasser structure file, the residue length is equal to the protein length
# here to make it comparable with the swiss model, we add two new columns of the pdb annotation information
protein_list0['sstart2'] = [1]*len(protein_list0['geneID'])
protein_list0['send2'] = protein_list0['protein_length']
protein_list0['locus'] = protein_list0['geneID']
protein_list0.to_excel('/Users/xluhon/Documents/GitHub/De-nevo-protein-3D-structure-yeast/data/pdb_itasser.xlsx')
# lastly we will calculate the residue distance matrix