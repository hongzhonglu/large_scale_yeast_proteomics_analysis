# The script is used to prepare sh file for the us-align calculation

import os
import shutil


# for the yeast
# pdb dir
dir0 = "/Users/xluhon/Documents/data_for_structure_align/UP000000559_237561_CANAL_v3/"
all_file = os.listdir(dir0)

# make new dir
os.mkdir("/Users/xluhon/Documents/" + "alphafold_pdb_v3")

# copy files
output = "/Users/xluhon/Documents/" + "alphafold_pdb_v3"
for x in all_file:
    print(x)
    if ".pdb.gz" in x:
        try:
            shutil.copy(dir0 + str(x), output)
        except:
            pass

# uncompress the gz files

"""
cd /Users/xluhon/Documents/alphafold_pdb_v3
gunzip -k *.gz
"""

# further remove file in .gz format
"""
rm *.pdb.gz
"""




# classify the proteins
s1 = open('/Users/xluhon/Documents/data_for_structure_align/protein_seq_for_alphafold/UP000002311_559292.fasta').readlines()
s2 = [x for x in s1 if '>' in x]


# uncharacterized protein
s2_c = [x for x in s2 if 'Uncharacterized' in x]
s2_c2 = [x for x in s2 if 'Putative uncharacterized' in x]
s2_c = s2_c + s2_c2

s2_other = list(set(s2)-set(s2_c))

# remove the detailed function annotation
s2_c = [x.split(' ')[0] for x in s2_c]
s2_other = [x.split(' ')[0] for x in s2_other]

# check the name as the protein stucture ID
s2_c = [x.split('|')[1] for x in s2_c]
s2_other = [x.split('|')[1] for x in s2_other]


s2_c = ['AF-' + x + '-F1-model_v3.pdb' for x in s2_c]
s2_other = ['AF-' + x + '-F1-model_v3.pdb' for x in s2_other]


# from here, all code run in mac
# preprocess txt file
with open("data/organism_result_sf.sh", "w") as outfile:
    for i in s2_c:
        for j in s2_other:
            s_out = i +'@@' + j
            commond_line = './USalign ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + i + ' ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + j +\
                           ' >> ' + '/home/yeast/Documents/tm_out/' +s_out + '.txt' +'\n'
            print(commond_line)
            outfile.write(commond_line)





# orthlogFinder process
import pandas as pd
dir0 = '/Users/xluhon/Documents/data_for_structure_align/protein_seq_for_alphafold/OrthoFinder/Results_Sep02/Orthogroups/Orthogroups.xlsx'
ortholog_relation = pd.read_excel(dir0)
ortholog_relation0 = ortholog_relation.dropna()

#process:
with open("data/organism_result_two_strains.sh", "w") as outfile:
    for i, x in ortholog_relation0.iterrows():
        print(i, x)
        y1 = x[1].split(', ')
        y2 = x[2].split(', ')

        # check the name as the protein stucture ID
        y1 = [x.split('|')[1] for x in y1]
        y2 = [x.split('|')[1] for x in y2]

        y1 = ['AF-' + x + '-F1-model_v3.pdb' for x in y1]
        y2 = ['AF-' + x + '-F1-model_v3.pdb' for x in y2]

        # try to save the file
        for i in y1:
            for j in y2:
                s_out = i + '@@' + j
                commond_line = './USalign ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_canal/' + i + ' ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + j + \
                               ' >> ' + '/home/yeast/Documents/tm_out_two_strains/' + s_out + '.txt' + '\n'
                print(commond_line)
                outfile.write(commond_line)