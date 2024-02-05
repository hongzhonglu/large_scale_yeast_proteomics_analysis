# this script is to transform the unit of proteomics datasets from mmol/gDW or g/gDW into molecular/cell
# based on the common existing absolute protein copy per cell due to fact that some protein copy may be higher conserved across samples
# which may be used a ruler to calibrate the protein abundance.
# 2021-11-16


import pandas as pd
import numpy as np
import statistics
# import the absolute protein copy per cell from different source
protein_copy = pd.read_excel("data/proteomics/protein_copy_combine.xlsx") # unit is molecular/cell


protein_copy_sgd = pd.read_csv("data/proteomics/sce_protein_abundance_sgd.tsv", sep='\t')
# rank the protein based on SD/Median value
protein_copy_sgd["sd/median"] = protein_copy_sgd['Abundance_SD'] / protein_copy_sgd['Abundance_median']
protein_copy_sgd = protein_copy_sgd.sort_values(by=['sd/median'], ascending=False)
protein_copy_sgd = protein_copy_sgd[protein_copy_sgd["sd/median"] > 0]
protein_copy_sgd = protein_copy_sgd[protein_copy_sgd["sd/median"] <= 0.25]
gene_ruler_potential = protein_copy_sgd["Systematic_name"].tolist()



# though the variance is quite small for these proteins, but some proteins may not be detected in most samples
# filter 2
column_all = protein_copy.columns
column_all = [x for x in column_all if "cell_system_2018" not in x]
protein_copy2 = protein_copy[column_all]
protein_copy_for_rulers = protein_copy2[protein_copy2["gene"].isin(gene_ruler_potential)]

# cycle for each row
cov_percentage = []
sd_per_mean = []
for i, x in protein_copy_for_rulers.iterrows():
    print(i)
    s = list(x)
    s = s[1:]
    total_num = len(s)
    s = [float(x) for x in s]
    s1 = [x for x in s if ~np.isnan(x)]
    cov_percentage.append(len(s1)/len(s))
    if len(s1) > 1:
        average = sum(s1)/len(s1)
        sd = statistics.stdev(s1)
        sd_ave = sd/average
        sd_per_mean.append(sd_ave)
    else:
        sd_per_mean.append(None)

protein_copy_for_rulers["cov_percentage"] = cov_percentage
protein_copy_for_rulers["sd_per_mean"] = sd_per_mean

protein_copy_for_rulers1 = protein_copy_for_rulers[protein_copy_for_rulers["sd_per_mean"] <=0.25]
protein_copy_for_rulers1 = protein_copy_for_rulers1[protein_copy_for_rulers1["cov_percentage"] >=0.5]



# one example to scale protein abundance based on protein copy ruler
# only check the  protein copy of YBL050W as the deviation is quite small across samples
copy_per_cell_YBL050W = 11448
abundance_YBL050W = 6.89555E-07 # mmol/gDW
mass_per_cell_YBL050W = 11448/6.022e+20 # mmol per cell
cell_number = abundance_YBL050W/mass_per_cell_YBL050W
# then build a coefficient
protein_copy = abundance_YBL050W/cell_number*6.022e+20 # equal to 11448
coefficient1 = 6.022e+20/cell_number
protein_copy = abundance_YBL050W * coefficient1 # equal to 11448



omics_combine = pd.read_excel("data/proteomics/omics_measured_combine.xlsx")
sampleID = omics_combine.columns
sampleID = [x for x in sampleID if "all_gene" not in x]
total_copy = []
for x in sampleID:
    sample_id = x
    gene_ref = "YBL050W"
    omics_combine1 = omics_combine[["all_gene", sample_id]]
    omics_combine1_ref = omics_combine1[omics_combine1["all_gene"] == gene_ref]
    abundance_YBL050W = omics_combine1_ref[sample_id].tolist()[0]
    mass_per_cell_YBL050W = 11448 / 6.022e+20  # mmol per cell
    cell_number = abundance_YBL050W / mass_per_cell_YBL050W
    # then build a coefficient
    protein_copy = abundance_YBL050W / cell_number * 6.022e+20
    coefficient_new = 6.022e+20 / cell_number
    curation_coefficient = 2.316173078784551 # the is another scaling coefficient from carl datasets asume that total protein copy is smaller than 1e+8
    omics_combine1["copy_per_cell"] = omics_combine1[sample_id] * coefficient_new / curation_coefficient
    omics_combine1 = omics_combine1.dropna()
    value0 = sum(omics_combine1["copy_per_cell"].tolist()) / (1e+8)
    total_copy.append(value0)
