# similarity calculation based on TM score

import os
import pandas as pd
import shutil
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt

# run the code in sh files
# input the yeast-GEMs
GEM_pro_id = pd.read_excel("result/yeast_gem_with_structure_id_and_score.xlsx")
pdb = GEM_pro_id["structure_id"].tolist()

# general calculation steps
id_list = pdb
out_file = "/Users/xluhon/Documents/TMscore_result/"
os.mkdir(out_file)
combine0 = list(combinations(id_list,2))
template = "./Tmscore /Users/xluhon/Documents/alphafold_pdb/structure1.pdb /Users/xluhon/Documents/alphafold_pdb/structure2.pdb > "+ "/Users/xluhon/Documents/TMscore_result/" + "output.txt"
with open("data/organism_result_new.sh", "w") as outfile:
    for i,xx in enumerate(combine0):
        print(i,xx)
        ss1 = xx[0]
        ss2 = xx[1]
        ss10 = ss1.replace(".pdb", "")
        ss20 = ss2.replace(".pdb", "")
        out_name = ss10 + "@" + ss20
        file1 = template.replace("output", out_name).replace("structure1.pdb", ss1).replace("structure2.pdb", ss2)
        file1 = file1 + "\n"
        out_name1 = ss20 + "@" + ss10
        file2 = template.replace("output", out_name1).replace("structure1.pdb", ss2).replace("structure2.pdb", ss1)
        file2 = file2 + "\n"
        outfile.write(file1)
        outfile.write(file2)


# 11:34

# only focus on interesing subsystems
sce_kegg_pathway = pd.read_excel("data/sce_kegg_pathway.xlsx")
# take EMP, TCA, amino acids synthesis and secondary metabolites synthesis
subsystem = ["Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)", "Biosynthesis of secondary metabolites", "Biosynthesis of amino acids", "Oxidative phosphorylation", "Pentose phosphate pathway"]
gene_select = sce_kegg_pathway[sce_kegg_pathway["name"].isin(subsystem)]

ref_EMP = sce_kegg_pathway[sce_kegg_pathway["name"]=="Glycolysis / Gluconeogenesis"]["id"].tolist()
ref_tca = sce_kegg_pathway[sce_kegg_pathway["name"]=="Citrate cycle (TCA cycle)"]["id"].tolist()
ref_sm = sce_kegg_pathway[sce_kegg_pathway["name"]=="Biosynthesis of secondary metabolites"]["id"].tolist()
ref_aa = sce_kegg_pathway[sce_kegg_pathway["name"]=="Biosynthesis of amino acids"]["id"].tolist()

"""
# analyze the volume per kda based on subsystems
sce_kegg_pathway_part = sce_kegg_pathway[sce_kegg_pathway["name"].isin(subsystem)]
sce_kegg_pathway_part["volume_per_kda"] = 1000*sce_kegg_pathway_part["pro_volume"]/sce_kegg_pathway_part["MW"]
sns.catplot(x="name", y="volume_per_kda", kind="box", data=sce_kegg_pathway_part)
plt.xticks(rotation=90)
plt.savefig('result/figure/volume_per_kda_for_each_subsystems.pdf', bbox_inches='tight')
"""


# example one: for one list
# id_list = ref_EMP
# combine0 = list(combinations(id_list,2))


# combinations from two subsystems
combine0 = []
for x in ref_EMP:
    for y in ref_sm:
        s1 = [x, y]
        combine0.append(s1)

# general calculation steps
out_file = "/Users/xluhon/Documents/TMscore_result_emp_sm/"
out_sh_file = "/Users/xluhon/Documents/organism_result_emp_sm.sh"
os.mkdir(out_file)
template = "./Tmscore /Users/xluhon/Documents/alphafold_pdb/structure1.pdb /Users/xluhon/Documents/alphafold_pdb/structure2.pdb > "+ out_file + "output.txt"
with open(out_sh_file, "w") as outfile:
    for i,xx in enumerate(combine0):
        print(i,xx)
        ss1 = xx[0]
        ss2 = xx[1]
        ss10 = ss1.replace(".pdb", "")
        ss20 = ss2.replace(".pdb", "")
        out_name = ss10 + "@" + ss20
        file1 = template.replace("output", out_name).replace("structure1.pdb", ss1).replace("structure2.pdb", ss2)
        file1 = file1 + "\n"
        out_name1 = ss20 + "@" + ss10
        file2 = template.replace("output", out_name1).replace("structure1.pdb", ss2).replace("structure2.pdb", ss1)
        file2 = file2 + "\n"
        outfile.write(file1)
        outfile.write(file2)

os.system("chmod u+x " + out_sh_file)


# collect and analyze the datasets
out_file = "/Users/xluhon/Documents/TMscore_result_emp/"

def getAlltmscoreFromDir(dir0):
    all_file = os.listdir(dir0)
    tm_score = []
    for x in all_file:
        print(x)
        input0 = dir0 + x
        ss = open(input0).readlines()
        value0 = ss[16]
        value = value0.split(" ")
        value1 = float(value[5])
        tm_score.append(value1)
    return tm_score

tm_score_emp = getAlltmscoreFromDir(dir0="/Users/xluhon/Documents/TMscore_result_emp/")
tm_score_emp_aa = getAlltmscoreFromDir(dir0="/Users/xluhon/Documents/TMscore_result_emp_aa/")
tm_score_emp_tca = getAlltmscoreFromDir(dir0="/Users/xluhon/Documents/TMscore_result_emp_tca/")
tm_score_emp_sm = getAlltmscoreFromDir(dir0="/Users/xluhon/Documents/TMscore_result_emp_sm/")




from scipy.stats import ttest_ind
ttest_ind(tm_score_emp, tm_score_emp_aa)
ttest_ind(tm_score_emp, tm_score_emp_tca)



pro_g1 = pd.DataFrame({"TM-score": tm_score_emp})
pro_g1["group"] = "EMP"
pro_g2 = pd.DataFrame({"TM-score": tm_score_emp_tca})
pro_g2["group"] = "EMP_TCA"
pro_g3 = pd.DataFrame({"TM-score": tm_score_emp_aa})
pro_g3["group"] = "EMP_AA"
pro_g4 = pd.DataFrame({"TM-score": tm_score_emp_sm})
pro_g4["group"] = "EMP_SM"

# combine all pandas
pro_c = pd.concat([pro_g2, pro_g3, pro_g4], axis=0)
sns.catplot(x="group", y="TM-score", order=["EMP_TCA", "EMP_AA", "EMP_SM"], kind="box", data=pro_c)
plt.ylim(0,0.2)


