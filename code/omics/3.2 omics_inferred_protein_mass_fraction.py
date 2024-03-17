# Note:
# all these analysis is based on the protein abundance in the unit of mmol/DCW

# import self function
from src.protein_process import *

# data preprocess
# input the protein abundance data in the unit of mmol/g DCW
omics_combine_auto = pd.read_excel("data/proteomics/omics_measured_combine_with_more_samples.xlsx")

# test
omics_combine_auto.columns = omics_combine_auto.columns.str.replace('all_gene', 'gene')
# change the unit from mmol/gDCW into g/gDCW
# some general datasets
# Get the molecular weight data using the data from SGD with more genes
mw = pd.read_csv("data/sce_protein_weight.tsv", sep="\t")
mw = mw[["locus","proteins_molecular_weight"]]
mw.columns = ["gene name", "MW"]
mw["MW_Kda"] = mw["MW"]/1000
# ID mapping between uniprot ID and gene locus IDs
id_mapping = pd.read_excel("data/uniprotGeneID_mapping.xlsx")

omics_combine_input = omics_combine_auto
omics_combine_input["MW_Kda"] = singleMapping(mw["MW_Kda"], mw["gene name"], omics_combine_input["gene"])
omics_combine_input = omics_combine_input[~omics_combine_input["MW_Kda"].isna()]

all_colum = omics_combine_input.columns
all_colum1 = [x for x in all_colum if x !='MW_Kda']
all_colum2 = [x for x in all_colum1 if x !='gene']

omics_combine_input2 = omics_combine_input[all_colum2]

# change mmol/gDW into g/gDW:
for x in all_colum2:
    omics_combine_input2[x] = omics_combine_input[x]*omics_combine_input["MW_Kda"]

omics_combine_input2['gene'] = omics_combine_input['gene']

# test the above code
out = ProMassRatio_Organelle(protein_abundance=omics_combine_input2, compartment_type="organelle")
out.to_excel("data/proteomics/ProMassRatio_across_compartment.xlsx")



