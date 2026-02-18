# note:
# this part can be rewritten as a function
import pandas as pd
import os
# import self function
from src.protein_process import *

compartment_all00 = getCompartmentGeneList(type="all")
ss0 = compartment_all00['nucleoplasm']




# reanalyze the data set
def getCompartmentGeneList(type="all"):
    """
    This function to build a compartment dict, with which we can get the gene list from the compartment name

    :param filter:
    :return:
    """

    # Input the datasets from paxDB
    # compartment = pd.read_csv("data/protein_location_sce.tsv", sep='\t')
    # extract compartment
    # compartment.columns = ['DBID', 'Systematic_name', 'Organism', 'Standard_name', 'Gene_name', 'GO_Qualifier', 'GO_Identifier', 'GO_Name', 'GO_Namespace', 'Ontology_Description', 'Annot_Type']
    # input gene annotation in Uniprot
    sce_gene = pd.read_excel("data/uniprotGeneID_mapping.xlsx")

    # using the updated version in 2026
    # Input the datasets from paxDB
    compartment = pd.read_csv("data/alliancemine_results_2026-01-03T10-33-37.tsv", sep='\t')
    compartment = compartment.iloc[:,0:11]
    # extract compartment
    compartment.columns = ['DBID', 'Systematic_name',  'Standard_name', 'Feature_type', 'Gene_qualifier','GO_Identifier', 'GO_Name', 'GO_Namespace', 'Annot_Type', 'GO_Qualifier', 'Organism']

    compartment = compartment[compartment['Systematic_name'].isin(sce_gene['GeneName'])]
    compartment1 = compartment[compartment["GO_Namespace"] == "cellular_component"]

    # filter out compartment with "complex" or "subunit"
    compartment2 = compartment1[~compartment1["GO_Name"].str.contains("complex")]
    compartment2 = compartment2[~compartment1["DBID"].str.contains("CPX-")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("subunit")]
    # firstly remove some general cellular component
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("snRNP")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("spindle")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("actin")]
    compartment2 = compartment2[~compartment2["GO_Name"].str.contains("cellular_component")]

    # more filter
    compartment2 = compartment2[compartment2['GO_Qualifier']=='located_in']  # only select "located"


    # analyze the annotation type
    annotation_type = compartment2["Annot_Type"].tolist()
    print("Annotation type:")
    print(list(set(annotation_type)))
    ratios = pd.Series(annotation_type).value_counts(normalize=True).round(4)
    print(ratios)

    # add annotation
    go_evidence_codes = {
    'IDA': 'Inferred from Direct Assay',
    'IEA': 'Inferred from Electronic Annotation',
    'HDA': 'Inferred from High Throughput Direct Assay',
    'IBA': 'Inferred from Biological aspect of Ancestor',
    'TAS': 'Traceable Author Statement',
    'ISA': 'Inferred from Sequence Alignment',
    'ND': 'No biological Data available',  # no any evidence
    'IPI': 'Inferred from Physical Interaction',
    'NAS': 'Non-traceable Author Statement',
    'ISS': 'Inferred from Sequence or structural Similarity',
    'IMP': 'Inferred from Mutant Phenotype',
    'IC': 'Inferred by Curator',
    'EXP': 'Inferred from Experiment',
    'IGI': 'Inferred from Genetic Interaction',
    'ISM': 'Inferred from Sequence Model',
    'ISO': 'Inferred from Sequence Orthology'}

    exp = [
    'IDA',  # Inferred from Direct Assay (直接实验测定)
    'HDA',  # Inferred from High Throughput Direct Assay (高通量直接实验)
    'IPI',  # Inferred from Physical Interaction (物理相互作用实验)
    'IMP',  # Inferred from Mutant Phenotype (突变体表型实验)
    'IGI',  # Inferred from Genetic Interaction (遗传相互作用实验)
    'IC',   # Inferred by Curator (人工curator推断)
    'EXP',  # Inferred from Experiment (实验证据伞状代码)
    'TAS',  # Traceable Author Statement (可追溯作者声明，已弃用)
    'NAS',  # Non-traceable Author Statement (不可追溯作者声明，已弃用)
     ]


    comput = [
    'IEA',  # Inferred from Electronic Annotation (纯电子/自动注释)
    'IBA',  # Inferred from Biological aspect of Ancestor (系统发育祖先推断)
    'ISA',  # Inferred from Sequence Alignment (序列比对)
    'ISS',  # Inferred from Sequence or structural Similarity (序列/结构相似性，广义)
    'ISM',  # Inferred from Sequence Model (序列模型)
    'ISO'   # Inferred from Sequence Orthology (正交同源性)
    ]

    # here if we remove "computational"
    # it seems that this method is wrong, as a protein could be located in multiple compartment. Thus it removes too much information
    compartment_with_evidence = compartment2[compartment2["Annot_Type"].isin(exp)]
    compartment_with_computation = compartment2[compartment2["Annot_Type"].isin(comput)]

    len(set(compartment_with_evidence["Systematic_name"].tolist()))
    len(set(compartment_with_computation["Systematic_name"].tolist()))

    # in one procedure, if a protein has no compartment annotation from manual and high-throughput, then the computational is used!
    # compartment_addition = compartment_with_computation[~compartment_with_computation["Systematic_name"].isin(compartment_with_evidence["Systematic_name"])]
    # compartment_combine = pd.concat([compartment_with_evidence, compartment_addition])

    # just for the test
    # compartment_combine_test = compartment_combine[compartment_combine['GO_Name'] == "nuclear membrane"]

    # build the dict
    compartment_dict_all = {}
    for i, x in compartment2.iterrows():
        print(i, x)
        if x['GO_Name'] in compartment_dict_all.keys():
            compartment_dict_all[x['GO_Name']] = list(set(compartment_dict_all[x['GO_Name']] + [x["Systematic_name"]]))
        else:
            compartment_dict_all[x['GO_Name']] = list(set([x["Systematic_name"]]))
    # filter
    compartment_dict_all0 = {}
    for key in compartment_dict_all.keys():
        print(key)
        value = compartment_dict_all[key]
        if len(value) >= 6:
            compartment_dict_all0[key] = value
        else:
            pass

    # for compartment annotation from manual evidence
    compartment_dict2 = {}
    for i, x in compartment_with_evidence.iterrows():
        print(i, x)
        if x['GO_Name'] in compartment_dict2.keys():
            compartment_dict2[x['GO_Name']] = list(set(compartment_dict2[x['GO_Name']] + [x["Systematic_name"]]))
        else:
            compartment_dict2[x['GO_Name']] = list(set([x["Systematic_name"]]))
    # filter
    compartment_dict20 = {}
    for key in compartment_dict2.keys():
        print(key)
        value = compartment_dict2[key]
        if len(value) >= 6:
            compartment_dict20[key] = value
        else:
            pass

    # for compartment annotation from computational evidence
    compartment_dict3 = {}
    for i, x in compartment_with_computation.iterrows():
        print(i, x)
        if x['GO_Name'] in compartment_dict3.keys():
            compartment_dict3[x['GO_Name']] = list(set(compartment_dict3[x['GO_Name']] + [x["Systematic_name"]]))
        else:
            compartment_dict3[x['GO_Name']] = list(set([x["Systematic_name"]]))
    # filter
    compartment_dict30 = {}
    for key in compartment_dict3.keys():
        print(key)
        value = compartment_dict3[key]
        if len(value) >= 6:
            compartment_dict30[key] = value
        else:
            pass
    if type == "all":
        return compartment_dict_all0
    elif type == "manual":
        return compartment_dict20
    else:
        return compartment_dict30

def gene_location_curation_sce(organelle0):
    # use some manually checked gene compartment definion
    # if the manual curated gene number for one compartment is larger, nealy equal to computational, then use the manual curation
    # otherwise using the computation prediction???
    # input the annotation from sgd
    # organelle0 = getCompartmentGeneList(type="all") # this is just for the test
    # organelle0 = compartment_all0
    gene_plasma_membrane = pd.read_excel("data/sce_compartment_curation/2026_curated/plasma_annotations.xlsx")
    gene_cell_wall = pd.read_excel("data/sce_compartment_curation/2026_curated/fungal_type_cell_wall_annotations.xlsx")
    gene_fungal_type_vacuole_membrane = pd.read_excel("data/sce_compartment_curation/2026_curated/fungal_type_vacuole_membrane_annotations.xlsx")
    gene_nucleolus = pd.read_excel("data/sce_compartment_curation/2026_curated/nucleolus_annotations.xlsx")
    gene_cytoplasm = pd.read_excel("data/sce_compartment_curation/2026_curated/cytoplasm_annotations.xlsx")
    gene_cytosol = pd.read_excel("data/sce_compartment_curation/2026_curated/cytosol_annotations.xlsx")
    gene_nucleus = pd.read_excel("data/sce_compartment_curation/2026_curated/nucleus_annotations.xlsx")
    gene_peroxisome = pd.read_excel("data/sce_compartment_curation/2026_curated/peroxisome_annotations.xlsx")

    gene_endoplasmic_reticulum = pd.read_excel("data/sce_compartment_curation/2026_curated/endoplasmic_reticulum_annotations.xlsx")
    gene_endoplasmic_reticulum_membrane = pd.read_excel("data/sce_compartment_curation/2026_curated/endoplasmic_reticulum_membrane_annotations.xlsx")
    gene_Golgi_apparatus = pd.read_excel("data/sce_compartment_curation/2026_curated/Golgi_apparatus_annotations.xlsx")
    gene_Golgi_membrane_annotations = pd.read_excel("data/sce_compartment_curation/2026_curated/Golgi_membrane_annotations.xlsx")

    ribo_gene_double_check = 'YMR242C, YOR312C, YDL081C, YJR094W-A, YCR031C, YPL198W, YHR010W, YPR043W, YGL031C, YHL015W, YIL069C, YKL006W, YGL030W, YBR191W, YOL040C, YOL121C, YJL177W, YJL189W, YLR048W, YJR145C, YGL135W, YEL054C, YGL123W, YLR185W, YDR382W, YFL034C-A, YDR447C, YGR034W, YBL087C, YPL220W, YER056C-A, YPL081W, YLR340W, YDR418W, YKL156W, YLR388W, YGL147C, YBR031W, YPL249C-A, YLR325C, YLR406C, YDR450W, YKL180W, YBR048W, YLR249W, YDL184C, YLR167W, YLR264W, YHR203C, YIL148W, YDL133C-A, YPR102C, YOL039W, YML024W, YGL076C, YNL069C, YDL075W, YOR293W, YOR063W, YNL178W, YGR085C, YML063W, YMR142C, YKR094C, YHL033C, YLR075W, YDL061C, YDR471W, YJL190C, YMR194W, YER117W, YNL096C, YFR032C-A, YLR441C, YPR132W, YBL072C, YGR118W, YLR367W, YBL027W, YBR084C-A, YGL103W, YLR029C, YLR287C-A, YMR116C, YPL090C, YDR500C, YPL143W, YPL131W, YIL133C, YJR123W, YNL302C, YMR143W, YDL130W, YML026C, YIL052C, YLR344W, YNL067W, YOR167C, YLR333C, YHR141C, YOR369C, YLR061W, YGR027C, YGR148C, YDL083C, YOR182C, YDR025W, YBR181C, YOL127W, YPL079W, YNL162W, YBR189W, YLL045C, YNL301C, YER131W, YHL001W, YER074W, YDL191W, YOR096W, YHR021C, YDL082W, YML073C, YMR121C, YJL136C, YJL191W, YBL092W, YDR012W, YFR031C-A, YKR057W, YLR448W, YER102W, YGR214W, YMR230W, YOL120C, YGL189C, YOR234C, YIL018W, YDL136W, YDR064W'
    ribo_gene_double_check = ribo_gene_double_check.split(', ')

    # mitochondrion specific
    gene_mitochondrion = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrion_annotations.xlsx")
    gene_m_Outer_membrane = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrial_outer_membrane_annotations.xlsx")
    gene_m_Inner_membrane = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrial_inner_membrane_annotations.xlsx")
    gene_m_OI_space = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrial_intermembrane_space_annotations.xlsx")
    gene_m_matrix = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrial_matrix_annotations.xlsx")
    gene_m_unassigned = pd.read_excel("data/sce_compartment_curation/2026_curated/mitochondrial_unassigned.xlsx")

    organelle1 = organelle0.copy()
    organelle1['mitochondrion_unassigned'] = gene_m_unassigned['gene'].tolist()

    organelle0_update = {}
    for y in organelle1.keys():
        print(y)
        if y == "plasma membrane":
            genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
        elif y == "fungal-type cell wall":
            genes_select = gene_cell_wall["gene"].tolist()
        elif y == "fungal-type vacuole membrane":
            genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
            genes_select = [x for x in genes_select if x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
        elif y == "endosome":
            genes_select = organelle1[y]
            genes_select = [x for x in genes_select if x not in ["YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
        #elif y == "nucleolus":
        #    genes_select = gene_nucleolus["gene"].tolist()  # for the test
        elif y == "cytoplasm":
            genes_select = gene_cytoplasm["gene"].tolist()  # for the test
        elif y == "cytosol":
            genes_select = gene_cytosol["gene"].tolist()  # for the test
        elif y == "mitochondrion":
            genes_select = gene_mitochondrion["gene"].tolist()  # for the test
        elif y == "mitochondrial outer membrane":
            genes_select = gene_m_Outer_membrane["gene"].tolist()  # for the test
        elif y == "mitochondrial inner membrane":
            genes_select = gene_m_Inner_membrane["gene"].tolist()  # for the test
        elif y == "mitochondrial intermembrane space":
            genes_select = gene_m_OI_space["gene"].tolist()  # for the test
        elif y == "mitochondrial matrix":
            genes_select = gene_m_matrix["gene"].tolist()  # for the test
        elif y == "nucleus":
            genes_select = gene_nucleus["gene"].tolist()  # for the test
        elif y == "peroxisome":
            genes_select = gene_peroxisome["gene"].tolist()  # for the test
        #elif y == "endoplasmic reticulum":
        #    genes_select = gene_endoplasmic_reticulum["gene"].tolist()  # for the test
        elif y == "endoplasmic reticulum membrane":
            genes_select = gene_endoplasmic_reticulum_membrane["gene"].tolist()  # for the test
        elif y == "Golgi apparatus":
            genes_select = gene_Golgi_apparatus["gene"].tolist()  # for the test
        elif y == "Golgi membrane":
            genes_select = gene_Golgi_membrane_annotations["gene"].tolist()  # for the test
        elif y == "ribosome":
            genes_select = ribo_gene_double_check
        else:
            genes_select = organelle1[y]
        organelle0_update[y] = list(filter(lambda x: str(x) != 'nan', genes_select))
       # remove duplicates
        for key in organelle0_update:
            organelle0_update[key] = list(dict.fromkeys(organelle0_update[key]))

    return organelle0_update


mass_fraction_final = pd.read_excel("data/proteomics/mass_fraction_combine.xlsx")

compartment_all0 = getCompartmentGeneList(type="all")
compartment_all = gene_location_curation_sce(organelle0=compartment_all0)
compartment_all['nucleoplasm'] = ss0 # nucleoplasm_proteins, this gene list is from SGD annotation in 2026!




compartment_list = list(compartment_all.keys())
compartment_out = ','.join(compartment_list)

organelles_main = [
    "mitochondrion", "nucleus", "endoplasmic reticulum", "Golgi apparatus",
    "fungal-type vacuole", "peroxisome", "endosome", "lipid droplet", "fungal-type cell wall",
    "plasma membrane", "P-body", "cytoplasmic stress granule", "spindle pole body", "ribosome", "cytosol","extracellular region", "nucleolus", "mitochondrial inner membrane"]

# "cytosol","extracellular region" 属于区室，而非细胞器
compartment1 = {k: v for k, v in compartment_all.items() if k in organelles_main}
# remove duplicates
compartment1 = {key: list(set(value)) for key, value in compartment1.items()}



# test for single compartment under different filteration strategy
compartment_single ={}
cc = 'cytosol'
cc_refine = 'cytosol' + '_refine'
compartment_single[cc] = compartment_all0[cc]
compartment_single[cc_refine] = compartment_all[cc]


def ProMassRatio_Organelle(protein_abundance, compartment=compartment1):
    """
    This function is used to calculate the organelle protein aboslute abundance as a whole
    :param protein_abundance:
    :param compartment_type:
    :return:
    """
    # test
    # sample ID information
    Sample_ID_select = list(protein_abundance.columns)
    Sample_ID_select = [x for x in Sample_ID_select if x != "gene"]
    # use some manually checked gene compartment definion
    # gene_plasma_membrane = pd.read_excel("data/gene_belong_plasma_membrane_annotations.xlsx")
    # all_compartment = ['fungal-type vacuole membrane']
    # gene_fungal_type_vacuole_membrane = pd.read_excel("data/gene_belong_fungal_type_vacuole_membrane_annotations.xlsx")
    # creat a dataframe to save the result

    all_compartment = list(compartment.keys())
    result1 = pd.DataFrame({"compartment": all_compartment})
    # run the cycle

    for col0 in Sample_ID_select:
        print(col0)
        value1 = []
        for y in all_compartment:
            print(y)
            # test
            # y = "plasma membrane"
            # col0 = "Glucose_phase_rep1(g/gDW)"
            pro_abundance = protein_abundance[['gene', col0]]
            pro_abundance.columns = ['gene', 'g/gDW']

            '''if y == "plasma membrane":
                genes_select = gene_plasma_membrane["gene"].tolist()  # for the test
            elif y == "fungal-type vacuole membrane":
                genes_select = gene_fungal_type_vacuole_membrane["gene"].tolist()  # for the test
                genes_select = [x for x in genes_select if
                                x not in ["YAL005C", "YLL024C"]]  # remove two genes for fungal type vacuole membrane
            elif y == "endosome":
                genes_select = compartment[y]
                genes_select = [x for x in genes_select if x not in ["YKR039W"]]  # remove one gene from endosome as this gene belongs to different compartments, also result in dramatic change in organelle protein volume.
            else:
                genes_select = compartment[y]'''
            genes_select = compartment[y]

            # get the sum
            pro_abundance.fillna(0, axis=1, inplace=True)
            pro_select = pro_abundance[pro_abundance['gene'].isin(genes_select)]
            sum_all = sum(pro_abundance['g/gDW'])
            sum_select = sum(pro_select['g/gDW'])
            ratio = sum_select/sum_all
            value1.append(ratio)
        result1[col0] = value1
    return result1

mass_fraction_final = pd.read_excel("data/proteomics/mass_fraction_combine.xlsx")
#mass_fraction_final = mass_fraction_final.iloc[:,0:9]

#out00 = ProMassRatio_Organelle(protein_abundance=mass_fraction_final, compartment=compartment_single)
out00 = ProMassRatio_Organelle(protein_abundance=mass_fraction_final, compartment=compartment1)
out00.to_excel("data/proteomics/main_organelle_fraction_test.xlsx")


out11 = ProMassRatio_Organelle(protein_abundance=mass_fraction_final, compartment=compartment_all)

len(compartment_all['endoplasmic reticulum'])
out11.to_excel("data/proteomics/all_organelle_fraction_test.xlsx")


len(list(set(ss0)-set(genes_select)))



# 逐行计算统计量
df = out00
df_stats = pd.DataFrame()
# 计算均值
df_stats['mean'] = df.mean(axis=1)
# 计算标准差
df_stats['std'] = df.std(axis=1)
# 计算分位数
df_stats['q25'] = df.quantile(0.25, axis=1)
df_stats['q50'] = df.quantile(0.5, axis=1)  # 中位数
df_stats['q75'] = df.quantile(0.75, axis=1)
# 计算最小值和最大值
df_stats['min'] = df.min(axis=1)
df_stats['max'] = df.max(axis=1)
# 计算范围（极差）
df_stats['range'] = df.max(axis=1) - df.min(axis=1)
df_stats.index = df['compartment'].tolist()














'''
compartment_main_c = gene_location_curation_sce(compartment_main)
compartment_main_c = {k: v for k, v in compartment_main_c.items() if k in organelles_main}


# quality check
plasma = compartment_main_c['plasma membrane']
mass_fraction_final_plasma = mass_fraction_final[mass_fraction_final['gene'].isin(plasma)]
column_sums1 = mass_fraction_final_plasma.sum()
column_sums2 = mass_fraction_final.sum()



from collections import defaultdict
def count_keys_per_value_optimized(input_dict):
    # 自动初始化值为 list 类型
    result = defaultdict(list)

    for key, value_list in input_dict.items():
        for item in value_list:
            result[item].append(key)

    return dict(result)  # 转回普通字典返回

output = count_keys_per_value_optimized(compartment_main_c)
print(output)











# new way
Sample_ID_select = list(mass_fraction_final.columns)
Sample_ID_select = [x for x in Sample_ID_select if x !="gene"]
output0 = {}
for col0 in Sample_ID_select:
    print(col0)
    value1 = []
    pro_abundance = mass_fraction_final[['gene', col0]]
    pro_abundance.columns = ['gene', 'g/gDW']
    pro_abundance = pro_abundance.dropna()

    gene = pro_abundance['gene'].tolist()
    gene_with_compartment = list(output.keys())
    gene_refine = list(set(gene) & set(gene_with_compartment))

    pro_abundance00 = pro_abundance[pro_abundance['gene'].isin(gene_refine)]
    gene_refine00 = pro_abundance00['gene'].tolist()

    pro_local = []
    for x in gene_refine00:
        print(list(output[x]))
        pro_local.append(list(output[x]))

    # 示例数据（实际替换为您的数据集，如从SGD或Ho et al.下载）
    data = {
        'protein': pro_abundance00['gene'].tolist(),
        'abundance': pro_abundance00['g/gDW'].tolist(),  # 每细胞拷贝数
        'localizations': pro_local
    }
    df = pd.DataFrame(data)

    # 您的细胞器列表
    organelles_main = [
        "mitochondrion", "nucleus", "endoplasmic reticulum", "Golgi apparatus",
        "fungal-type vacuole", "peroxisome", "endosome", "lipid droplet", "fungal-type cell wall",
        "plasma membrane", "P-body", "cytoplasmic stress granule", "ribosome", "cytosol",
        "cytoskeleton", "extracellular region"
    ]

    # Fractional分配计算
    allocation = {org: 0.0 for org in organelles_main}  # 初始化
    total_abundance = df['abundance'].sum()  # 总资源（守恒）

    for _, row in df.iterrows():
        locs = [org for org in row['localizations'] if org in organelles_main]  # 只计您的列表中
        if locs:  # 避免未注解蛋白
            fraction = row['abundance'] / len(locs)  # 均匀分配
            for loc in locs:
                allocation[loc] += fraction

    # 计算比例
    proportions = {org: (allocation[org] / total_abundance) * 100 for org in organelles_main}

    # 输出（排序显示）
    proportions_sorted = dict(sorted(proportions.items(), key=lambda x: x[1], reverse=True))
    for org, pct in proportions_sorted.items():
        print(f"{org}: {pct:.2f}%")
    output0[col0] = proportions_sorted



output_df = pd.DataFrame(output0)
output_df = output_df / 100
output_df.to_excel("data/sce_compartment_curation/mass_fraction_new_way.xlsx")
'''



