# similarity calculation based on TM score
# this script is used to generate the sh file


import os
import pandas as pd
import shutil
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt

# only focus on interesting subsystems
sce_kegg_pathway = pd.read_excel("data/sce_kegg_pathway.xlsx")
# take EMP, TCA, amino acids synthesis and secondary metabolites synthesis
all_subsystem = sce_kegg_pathway["name"].value_counts()
subsystem = ["Glycolysis / Gluconeogenesis",  "Citrate cycle (TCA cycle)", "Oxidative phosphorylation", "Pentose phosphate pathway"]
gene_select = sce_kegg_pathway[sce_kegg_pathway["name"].isin(subsystem)]


sce_kegg_pathway['id'] = sce_kegg_pathway['id'].str.replace("_v1",'_v3')
ref_EMP = sce_kegg_pathway[sce_kegg_pathway["name"]=="Glycolysis / Gluconeogenesis"]["id"].tolist()
ref_tca = sce_kegg_pathway[sce_kegg_pathway["name"]=="Citrate cycle (TCA cycle)"]["id"].tolist()
ref_op = sce_kegg_pathway[sce_kegg_pathway["name"]=="Oxidative phosphorylation"]["id"].tolist()
ref_ppp = sce_kegg_pathway[sce_kegg_pathway["name"]=="Pentose phosphate pathway"]["id"].tolist()


# automatically generate the sh file for the calculation

# example one: for one pathway
path1 = "emp"
path2 = "emp"
sh_name = path1 + "_vs_" + path2
id_list = ref_EMP
combine0 = list(combinations(id_list,2))
with open("data/" + sh_name + ".sh", "w") as outfile:
    outfile.write('mkdir' + ' ' + '/home/yeast/Documents/' + sh_name + '\n')
    for x in combine0:
        print(x)
        i = x[0]
        j = x[1]
        s_out = i + '@@' + j
        commond_line = './USalign ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + i + ' ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + j + \
                       ' >> ' + '/home/yeast/Documents/' + sh_name + '/' + s_out + '.txt' + '\n'
        print(commond_line)
        outfile.write(commond_line)




# example two: for two pathways
path1 = "emp"
path2 = "tca"
sh_name = path1 + "_vs_" + path2
with open("data/" + sh_name + ".sh", "w") as outfile:
    outfile.write('mkdir' + ' ' + '/home/yeast/Documents/' + sh_name + '\n')
    for i in ref_EMP:
        for j in ref_tca:
            s_out = i + '@@' + j
            commond_line = './USalign ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + i + ' ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + j + \
                           ' >> ' + '/home/yeast/Documents/' + sh_name + '/' + s_out + '.txt' + '\n'
            print(commond_line)
            outfile.write(commond_line)




# example two: for two pathways
path1 = "emp"
path2 = "ppp"
sh_name = path1 + "_vs_" + path2
with open("data/" + sh_name + ".sh", "w") as outfile:
    outfile.write('mkdir' + ' ' + '/home/yeast/Documents/' + sh_name + '\n')
    for i in ref_EMP:
        for j in ref_ppp:
            s_out = i + '@@' + j
            commond_line = './USalign ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + i + ' ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + j + \
                           ' >> ' + '/home/yeast/Documents/' + sh_name + '/' + s_out + '.txt' + '\n'
            print(commond_line)
            outfile.write(commond_line)




# example two: for two pathways
path1 = "emp"
path2 = "op"
sh_name = path1 + "_vs_" + path2
with open("data/" + sh_name + ".sh", "w") as outfile:
    outfile.write('mkdir' + ' ' + '/home/yeast/Documents/' + sh_name + '\n')
    for i in ref_EMP:
        for j in ref_op:
            s_out = i + '@@' + j
            commond_line = './USalign ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + i + ' ' + '/home/yeast/data_for_structure_align/alphafold_pdb_v3_yeast/' + j + \
                           ' >> ' + '/home/yeast/Documents/' + sh_name + '/' + s_out + '.txt' + '\n'
            print(commond_line)
            outfile.write(commond_line)
