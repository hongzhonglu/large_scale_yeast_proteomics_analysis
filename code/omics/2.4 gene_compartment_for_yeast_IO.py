# first, get the compartment annotation from uniprot
# if no, then find the ortholog from sce. transfer sce annotation to IO.
# lastly, for the remaining protein, using Deep learning to predict the compartment


import os

# import self function
from src.protein_process import *

# Part 1
# Initially check how many genes could find compartment


# pro_abundance = pd.read_csv("data/proteomics/sce_protein_abundance_sgd.tsv", sep='\t')

# input data from cell system, 2018
pro_abundance = pd.read_excel("data/nature_chemical_biology_datatset_2024/41589_2024_1571_MOESM3_ESM_only_IO.xlsx", sheet_name="Table 10c. abs_prot_IO_SD108")
#pro_abundance = pd.read_excel("data/nature_chemical_biology_datatset_2024/41589_2024_1571_MOESM3_ESM_only_IO.xlsx", sheet_name="Table 11b. rel_prot_IO")

IO_gene_list = pro_abundance[["Entry"]]

IO_annotation = pd.read_excel("data/nature_chemical_biology_datatset_2024/uniprotkb_proteome_UP000029867_2024_03_18.xlsx")
IO_annotation0 = IO_annotation[["Entry","Gene Ontology (cellular component)"]]

IO_gene_list0 = pd.merge(IO_gene_list,IO_annotation0,left_on="Entry", right_on="Entry", how="left")

IO_gene_no_location = IO_gene_list0[IO_gene_list0["Gene Ontology (cellular component)"].isnull()]
IO_gene_with_location = IO_gene_list0[~IO_gene_list0["Gene Ontology (cellular component)"].isnull()]
IO_gene_with_location.columns = ["gene", "compartment"]
IO_gene_with_location00 = splitAndCombine(IO_gene_with_location["compartment"],IO_gene_with_location["gene"],sep0=";")
IO_gene_with_location00.columns = ["gene", "compartment"]
cc_annotation = IO_gene_with_location00["compartment"].tolist()
cc_annotation = [x.split("[")[0] for x in cc_annotation]
IO_gene_with_location00["compartment"] = cc_annotation
IO_gene_with_location00["compartment"] = IO_gene_with_location00["compartment"].str.strip()
IO_gene_with_location00["source"] = "uniprot"







# check protein with no compartment annotation
# firstly input the gene ortholog relation between sce and io
gene_ortholog = pd.read_excel("data/nature_chemical_biology_datatset_2024/Orthogroups_between_sce_and_IO.xlsx")
gene_ortholog.columns = ["OG", "sce", "IO"]
gene_ortholog0 = splitAndCombine(gene_ortholog['IO'], gene_ortholog['sce'], sep0=",")
gene_ortholog0 = gene_ortholog0[['V2','V1']]
gene_ortholog1 = gene_ortholog0.dropna()
gene_ortholog1.columns =['IO_gene','sce_gene']
gene_ortholog1 = gene_ortholog1[gene_ortholog1["IO_gene"] != "NA"]
gene_special = gene_ortholog1["IO_gene"].tolist()
gene_special0 = [x.split("|")[1] for x in gene_special]
gene_ortholog1["IO_gene2"] = gene_special0

IO_gene_no_location["location_from_sce"] = multiMapping(gene_ortholog1["sce_gene"],gene_ortholog1["IO_gene2"],IO_gene_no_location["Entry"])
IO_gene_no_location_g1 = IO_gene_no_location[IO_gene_no_location["location_from_sce"].isnull()]
IO_gene_no_location_g2 = IO_gene_no_location[~IO_gene_no_location["location_from_sce"].isnull()]
# calculate the mapped genes of IO IN SCE
ii = IO_gene_no_location_g2["Entry"].tolist()
ss = IO_gene_no_location_g2["location_from_sce"].tolist()
ii0 = [x for x, y in zip(ii,ss) if len(y.split(",")) <= 2]
ii1 = [x for x, y in zip(ii,ss) if len(y.split(",")) > 2]

# get gene with sce ortholog genes (no more than two)
IO_gene_with_sce_ortholog = IO_gene_no_location_g2[IO_gene_no_location_g2["Entry"].isin(ii0)]
IO_gene_with_sce_ortholog = IO_gene_with_sce_ortholog[["Entry","location_from_sce"]]
IO_gene_with_sce_ortholog0 = splitAndCombine(IO_gene_with_sce_ortholog["location_from_sce"],IO_gene_with_sce_ortholog["Entry"], sep0=",")
IO_gene_with_sce_ortholog0.columns = ["gene_IO","gene_sce"]
IO_gene_with_sce_ortholog0["gene_sce"] = IO_gene_with_sce_ortholog0["gene_sce"].str.strip()









# get gene without location from above way
gene_with_no_location = ii1 + IO_gene_no_location_g1["Entry"].tolist()
# then for these 1046 genes, we can annotate them based on deep learning
from Bio import SeqIO
infile3 = 'data/nature_chemical_biology_datatset_2024/uniprotkb_proteome_UP000029867_2024_03_18.fasta'
select_sequences = []
for record in SeqIO.parse(infile3, "fasta"):
    print(record.id)
    ss0 = record.id
    ss1 = ss0.split("|")[1]
    if ss1 in gene_with_no_location:
        record.id = ss1
        select_sequences.append(record)
SeqIO.write(select_sequences, "data/nature_chemical_biology_datatset_2024/IO_select_seq_for_location_annotation.fasta", "fasta")
# then split the seq into 6
SeqIO.write(select_sequences[0:200], "data/nature_chemical_biology_datatset_2024/IO_select_seq_for_location_annotation1.fasta", "fasta")
SeqIO.write(select_sequences[200:400], "data/nature_chemical_biology_datatset_2024/IO_select_seq_for_location_annotation2.fasta", "fasta")
SeqIO.write(select_sequences[400:600], "data/nature_chemical_biology_datatset_2024/IO_select_seq_for_location_annotation3.fasta", "fasta")
SeqIO.write(select_sequences[600:800], "data/nature_chemical_biology_datatset_2024/IO_select_seq_for_location_annotation4.fasta", "fasta")
SeqIO.write(select_sequences[800:1000], "data/nature_chemical_biology_datatset_2024/IO_select_seq_for_location_annotation5.fasta", "fasta")
SeqIO.write(select_sequences[1000:], "data/nature_chemical_biology_datatset_2024/IO_select_seq_for_location_annotation6.fasta", "fasta")







# summarize the new compartment annotation
all_file = os.listdir("data/nature_chemical_biology_datatset_2024/compartment_annotation")
all_file = [x for x in all_file if x !=".DS_Store"]
def getCompartmentFromDL(file_name):
    c11 = pd.read_table(
        "data/nature_chemical_biology_datatset_2024/compartment_annotation/" + file_name + "/sub_cellular_prediction.txt",
        header=None)
    c11 = c11.iloc[:, [0, 11]]
    c11.columns = ["gene", "compartment"]
    c12 = pd.read_table(
        "data/nature_chemical_biology_datatset_2024/compartment_annotation/" + file_name + "/sub_organellar_prediction.txt",
        header=None)
    c12 = c12.iloc[:, [0, 46]]
    c12.columns = ["gene", "compartment"]
    c1_combine = pd.concat([c11, c12])
    c1_combine0 = splitAndCombine(c1_combine["compartment"], c1_combine["gene"], sep0="|", moveDuplicate=True)
    c1_combine0 = c1_combine0[c1_combine0["V2"] != ""]
    c1_combine0.columns = ["gene", "compartment"]
    c1_combine0["compartment"] = c1_combine0["compartment"].str.replace("prediction:", "").str.strip()
    gene_list = c1_combine0["gene"].tolist()
    gene_list = [x.split(" ")[0] for x in gene_list]
    c1_combine0["gene"] = gene_list
    c1_combine0["gene"] = c1_combine0["gene"].str.replace(">", "").str.strip()
    c1_combine00 = splitAndCombine(c1_combine0["compartment"], c1_combine0["gene"], sep0=",", moveDuplicate=True)
    c1_combine00.columns = ["gene", "compartment"]
    c1_combine00["compartment"] = c1_combine00["compartment"].str.strip().str.lower()
    return c1_combine00

out0 = getCompartmentFromDL(all_file[0])

out1 = getCompartmentFromDL(all_file[1])

out2 = getCompartmentFromDL(all_file[2])

out3 = getCompartmentFromDL(all_file[3])

out4 = getCompartmentFromDL(all_file[4])

out5 = getCompartmentFromDL(all_file[5])

new_compartment = pd.concat([out0, out1, out2, out3, out4, out5])
new_compartment["source"] = "MULocDeep"




#get sce gene annotation
compartment_dict20 = getCompartmentGeneList(type="all") # Remove some compartmental annotation only with computational evidence (keep experimental evidence)
compartment_dict20_update =gene_location_curation_sce(compartment_dict20)
# change it as a list
sce_gene_list = []
for x, y in compartment_dict20_update.items():
    print(x, y)
    yy = [y0 +"@" + x for y0 in y]
    sce_gene_list = sce_gene_list + yy
sce_gene_compartment_corrected = pd.DataFrame({"gene":sce_gene_list})
sce_gene_compartment_corrected[['gene', 'comparment']] = sce_gene_compartment_corrected['gene'].str.split('@', n=1, expand=True)
# transfer sce gene annotation to IO gene
IO_gene_with_sce_ortholog00 = pd.merge(IO_gene_with_sce_ortholog0,sce_gene_compartment_corrected,left_on="gene_sce",right_on="gene",how="outer")
IO_gene_with_sce_ortholog00 = IO_gene_with_sce_ortholog00[["gene_IO","comparment"]]
IO_gene_with_sce_ortholog01 = IO_gene_with_sce_ortholog00.drop_duplicates(keep='first')
IO_gene_with_sce_ortholog02 = IO_gene_with_sce_ortholog01.dropna()
IO_gene_with_sce_ortholog02["source"] = "ortholog_transfer"
IO_gene_with_sce_ortholog02.columns = ["gene","compartment","source"]







# combine all the compartment data together
IO_compartment = pd.concat([IO_gene_with_location00,new_compartment,IO_gene_with_sce_ortholog02])

# unify the name
IO_compartment["compartment"] = IO_compartment["compartment"].str.lower()
IO_compartment["compartment"] = IO_compartment["compartment"].str.replace("peroxisome membrane","peroxisomal membrane")
IO_compartment["compartment"] = IO_compartment["compartment"].str.replace("mitochondrion outer membrane","mitochondrial outer membrane")
IO_compartment["compartment"] = IO_compartment["compartment"].str.replace("mitochondrion inner membrane","mitochondrial inner membrane")
IO_compartment["compartment"] = IO_compartment["compartment"].str.replace("mitochondrion matrix","mitochondrial matrix")
IO_compartment["compartment"] = IO_compartment["compartment"].str.replace("golgi_apparatus","golgi apparatus")
IO_compartment["compartment"] = IO_compartment["compartment"].str.replace("golgi apparatus membrane","golgi membrane")
IO_compartment["compartment"] = IO_compartment["compartment"].str.replace("extracellular space","extracellular region")
IO_compartment["compartment"] = IO_compartment["compartment"].str.replace("cell membrane","plasma membrane")
IO_compartment["compartment"] = IO_compartment["compartment"].str.replace("nucleus membrane","nuclear membrane")
IO_compartment["compartment"] = IO_compartment["compartment"].str.replace("nucleus speckle","nuclear speckle")

one_more = IO_compartment["compartment"].tolist()
one_more0 = []
for xx in one_more:
    if xx == 'endoplasmic':
        xx1 = 'endoplasmic reticulum'
        one_more0.append(xx1)
    else:
        one_more0.append(xx)

IO_compartment["compartment"] = one_more0

# double check
IO_compartment1 = IO_compartment[IO_compartment["source"] !="MULocDeep"]
IO_compartment2 = IO_compartment[IO_compartment["source"].str.contains("MULocDeep")]
list(set(IO_compartment2["compartment"])-set(IO_compartment1["compartment"]))

# remove cytoplasm
# note: the protein compartment annotation from sce could affect the result of IO in compartment annotation
IO_compartment = IO_compartment[IO_compartment["compartment"]!="cytoplasm"]
IO_compartment.to_excel("data/nature_chemical_biology_datatset_2024/IO_gene_compartment.xlsx")



