# note:
# this part can be rewritten as a function

import os
import pandas as pd

# import self function
from src.protein_process import *




"""
# Part 1
# Initially check how many genes could find compartment
compartment = getCompartmentGeneList(type="all")  # based on the automatic way
#compartment = gene_location_curation_sce(organelle0=compartment)  # based on the SGD manual curation
all_compartment = list(compartment.keys())
mitochondrion_related = ['mitochondrion',
'mitochondrial outer membrane',
'mitochondrial intermembrane space',
'mitochondrial matrix',
'mitochondrial inner membrane']
compartment_m = {x:y for x, y in compartment.items() if x in mitochondrion_related}
"""


# reanalyze the data set
gene_mitochondrion_sgd = pd.read_excel("data/sce_compartment_curation/mitochondrion_annotations_manual_v2.xlsx") # this is based sgd annotation
mitochondrion_gene = gene_mitochondrion_sgd['gene'].tolist()

# combine
mitochondrial_inner_membrane = pd.read_excel("data/sce_compartment_curation/mitochondrial_inner_membrane_annotations_manual_v2.xlsx")

# combine
mitochondrial_outer_membrane = pd.read_excel("data/sce_compartment_curation/mitochondrial_outer_membrane_annotations_manual_v2.xlsx")

# combine
mitochondrial_matrix = pd.read_excel("data/sce_compartment_curation/mitochondrial_matrix_annotations_manual_v2.xlsx")

# combine
mitochondrial_intermembrane_space_annotations = pd.read_excel("data/sce_compartment_curation/mitochondrial_intermembrane_space_annotations_manual_v2.xlsx")




# input the NC paper
gene_mitochondrion_nc = pd.read_excel("data/sce_compartment_curation/mitochondrial_suborganelle.xlsx")
gene_mitochondrion_nc = gene_mitochondrion_nc.iloc[0:987:]
gene_m_Outer_membrane = gene_mitochondrion_nc[gene_mitochondrion_nc['Outer membrane'] == "X"]
gene_m_Inner_membrane = gene_mitochondrion_nc[gene_mitochondrion_nc['Inner membrane'] == "X"]
gene_m_OI_space = gene_mitochondrion_nc[gene_mitochondrion_nc['Inter-membrane space'] == "X"]
gene_m_matrix = gene_mitochondrion_nc[gene_mitochondrion_nc['Matrix'] == "X"]



gene_m_Outer_membrane = list(set(gene_m_Outer_membrane['gene'].tolist()) & set(mitochondrion_gene))
gene_m_Inner_membrane = list(set(gene_m_Inner_membrane['gene'].tolist()) & set(mitochondrion_gene))
gene_m_OI_space = list(set(gene_m_OI_space['gene'].tolist()) & set(mitochondrion_gene))
gene_m_matrix = list(set(gene_m_matrix['gene'].tolist()) & set(mitochondrion_gene))



# further combine different data together:
gene_m_Outer_membrane_v3 = list(set(gene_m_Outer_membrane + mitochondrial_outer_membrane['gene'].tolist()))
gene_m_Inner_membrane_v3 = list(set(gene_m_Inner_membrane + mitochondrial_inner_membrane['gene'].tolist()))
gene_m_OI_space_v3 = list(set(gene_m_OI_space + mitochondrial_intermembrane_space_annotations['gene'].tolist()))
gene_m_matrix_v3 = list(set(gene_m_matrix + mitochondrial_matrix['gene'].tolist()))
gene_assigned = list(set(gene_m_Outer_membrane_v3 + gene_m_Inner_membrane_v3 + gene_m_OI_space_v3 + gene_m_matrix_v3)) # 742 proteins
gene_unassigned = list(set(mitochondrion_gene)-set(gene_assigned)) # 454 proteins



gene_m_Outer_membrane_v3 = pd.DataFrame({"gene": gene_m_Outer_membrane_v3})
gene_m_Inner_membrane_v3 = pd.DataFrame({"gene": gene_m_Inner_membrane_v3})
gene_m_OI_space_v3 = pd.DataFrame({"gene": gene_m_OI_space_v3})
gene_m_matrix_v3 = pd.DataFrame({"gene": gene_m_matrix_v3})
gene_m_unassigned_v3 = pd.DataFrame({"gene": gene_unassigned})


gene_m_Outer_membrane_v3.to_excel("data/sce_compartment_curation/mitochondrial_outer_membrane_annotations_manual_v3.xlsx")
gene_m_Inner_membrane_v3.to_excel("data/sce_compartment_curation/mitochondrial_inner_membrane_annotations_manual_v3.xlsx")
gene_m_OI_space_v3.to_excel("data/sce_compartment_curation/mitochondrial_intermembrane_space_annotations_manual_v3.xlsx")
gene_m_matrix_v3.to_excel("data/sce_compartment_curation/mitochondrial_matrix_annotations_manual_v3.xlsx")
gene_m_unassigned_v3.to_excel("data/sce_compartment_curation/mitochondrial_unassigned_manual_v3.xlsx")


# input mitochrondrial protein list from  Absolute yeast mitochondrial proteome quantification reveals trade-off between biosynthesis and energy generation during diauxic shift
# input the data from PNAS
gene_mitochondrion_pnas = pd.read_excel("data/sce_compartment_curation/pnas.1918216117.sd01 VIP Francesca.xlsx", header=1)
sce_gene = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
sce_gene_select =sce_gene[sce_gene['Entry'].isin(gene_mitochondrion_pnas['Uniprot Accession'])]
gene_mitochondrion_pnas_F = sce_gene_select

# input mitochrondrial protein list from  jianye xia
gene_mitochondrion_xia = pd.read_excel("data/sce_compartment_curation/41467_2022_30513_MOESM5_ESM xia.xlsx", header=0)
gene_mitochondrion_xia.columns =['gene','short_name']
gene_mitochondrion_xia.to_excel("data/sce_compartment_curation/mitochondrion_xia.xlsx")



# input mitochrondrial protein list from  Definition of a High-Confidence Mitochondrial Proteome at Quantitative Scale
gene_mitochondrion_marcel = pd.read_excel("data/sce_compartment_curation/mmc4_cell_reports_2017_marcel.xlsx", header=0)
gene_mitochondrion_marcel = gene_mitochondrion_marcel.iloc[:,0:2]
gene_mitochondrion_marcel.columns =['gene','Entry']




compartment1 = {}
compartment1['m_nc'] = gene_mitochondrion_nc['gene'].tolist()
compartment1['m_sgd'] = gene_mitochondrion_sgd['gene'].tolist()
compartment1['m_xia'] = gene_mitochondrion_xia['gene'].tolist()
compartment1['m_pnas'] = gene_mitochondrion_pnas_F['GeneName'].tolist()
compartment1['m_marcel'] = gene_mitochondrion_marcel['gene'].tolist()




compartment_all00 = getCompartmentGeneList(type="all")
compartment_all = gene_location_curation_sce(organelle0=compartment_all00)
compartment_list = list(compartment_all.keys())
compartment_out = ','.join(compartment_list)

organelles_main = [
    "mitochondrion", "nucleus", "endoplasmic reticulum", "Golgi apparatus",
    "fungal-type vacuole", "peroxisome", "endosome", "lipid droplet", "fungal-type cell wall",
    "plasma membrane", "P-body", "cytoplasmic stress granule", "ribosome","cytosol",
    "cytoskeleton","extracellular region"]


compartment_all123 = {k: v for k, v in compartment_all.items() if k == "mitochondrion"}  # mitochondrion
compartment_all_other = {k: v for k, v in compartment_all.items() if "itochondri" not in k}
compartment_all_f = {**compartment_all123, **compartment_all_other}
from collections import defaultdict
def count_keys_per_value_optimized(input_dict):
    # 自动初始化值为 list 类型
    result = defaultdict(list)

    for key, value_list in input_dict.items():
        for item in value_list:
            result[item].append(key)

    return dict(result)  # 转回普通字典返回

output = count_keys_per_value_optimized(compartment_all_f)
print(output)

output_f = {k: v for k, v in output.items() if 'mitochondrion' in v} #mitochondrion
# check which protein exist only in mitochondrion
output_f2 = {k: v for k, v in output_f.items() if len(v) <=2}



compartment1['m_sgd_core'] = list(output_f2.keys())



# check the latest annotation
# 文件名为附件中的名称
filename = 'data/sce_compartment_curation/2026/plasma_membrane_annotations_2026.txt'
df = pd.read_csv(filename, sep='\t', skiprows=8, header=0)
df00 = df[df['Evidence']=='BSR']
print(set(df['Evidence'].tolist())) # 其中包括 ISM
df00 = df[df['Evidence'].str.contains('ISS')]

# only select "located"
df = df[df['Qualifier']=='located in']
df = df[~df['Systematic Name/Complex Accession'].str.contains('CPX-')]
gene_sgd_2026 = df.iloc[:,0:2]
gene_sgd_2026.columns =['short_name','gene']

compartment1['plasma_sgd_2026'] = list(set(gene_sgd_2026['gene'].tolist()))
','.join(list(set(gene_sgd_2026['gene'].tolist())))
#根据注释plasma可能包括液泡膜、糖酵解以及胞外蛋白等，所以目前很不准。
membrane = 'YDR055W, YPR192W, YKL046C, YFL026W, YJL170C, YOR273C, YHL016C, YNL294C, YNL173C, YML132W, YBR078W, YNR002C, YER060W, YNL323W, YOR381W, YPR149W, YDR032C, YOR030W, YBR054W, YLR092W, YBR043C, YIR019C, YML047C, YDR459C, YIL140W, YJL171C, YGR031C-A, YGR224W, YDL135C, YKR105C, YKL209C, YLL052C, YLR025W, YOL122C, YPR124W, YCR024C-A, YLR081W, YIR006C, YMR031C, YOL158C, YJL212C, YHR092C, YGR213C, YBR068C, YGL053W, YDR033W, YNL098C, YLR120C, YHR048W, YER123W, YML013W, YGL077C, YGL186C, YDL035C, YCR037C, YJL062W, YBR294W, YDR276C, YBL069W, YOL002C, YOR348C, YKR106W, YNL194C, YDR522C, YKL178C, YOR161C, YLR121C, YMR068W, YDL012C, YJL129C, YGL051W, YGL008C, YHR135C, YDR384C, YOR008C, YBR295W, YBR296C, YPR194C, YHR094C, YNL291C, YOR390W, YGR138C, YGR060W, YJR152W, YKL220C, YJL093C, YLR194C, YER185W, YDR040C, YER118C, YOR306C, YAL030W, YJL145W, YMR034C, YIL118W, YHL040C, YBR132C, YLR096W, YDR160W, YCL073C, YAR033W, YOR047C, YCL027W, YLL043W, YKR039W, YJR054W, YGL208W, YJR066W, YGR014W, YCR010C, YMR058W, YIL105C, YIL088C, YKL094W, YOR317W, YLR332W, YLR138W, YPR201W, YDL194W, YPR075C, YDL138W, YGR191W, YCR021C, YOL103W, YGR023W, YBL042C, YPR198W, YKR050W, YMR011W, YMR319C, YOR328W, YMR307W, YOR153W, YLR353W, YGR121C, YBR069C, YNR047W, YGR266W, YNL271C, YNL283C, YJR040W, YDL222C, YOL019W, YLL010C, YLR342W, YOR327C, YNL065W, YBR021W, YJL100W, YDR420W, YOR071C, YOR378W, YER020W, YML116W, YPL279C, YPL092W, YLR214W, YCL025C, YDL161W, YDR343C, YOL105C, YDR463W, YJR086W, YNL047C, YPL232W, YFR029W, YGR241C, YDR342C, YHR161C, YGL255W, YDR039C, YER056C, YJL058C, YOL009C, YOL078W, YLR310C, YGR055W, YLR019W, YGR217W, YPL265W, YER060W-A, YHR186C, YDL019C, YML052W, YDR345C, YJL156C, YAR031W, YNL154C, YOL152W, YDR210W, YIL147C, YHL044W, YPL036W, YPR165W, YPR032W, YKR093W, YOR011W, YKL126W, YNL142W, YNL275W, YLR237W, YBL029C-A, YMR017W, YOL020W, YDR011W, YIL121W, YLR130C, YGR041W, YCR075C, YGR198W, YOL109W, YOR104W, YNL192W, YDR536W, YMR183C, YOR171C, YOL130W, YGR260W, YOR086C, YPR156C, YCR098C, YLR373C, YIL120W, YIR038C, YKR055W, YDR034W-B, YDR208W, YDR090C, YER120W, YNL257C, YOR018W, YBR008C, YBR129C, YGR152C, YNR060W, YER145C, YKL217W, YJL198W, YLR413W, YPL058C, YEL063C, YHR096C, YGL114W, YLL028W, YNL180C, YML123C, YLR414C, YMR008C, YMR063W, YDR122W, YCR017C, YIL047C, YKL051W, YER143W, YBR016W, YHR005C, YNL268W, YOR101W, YMR238W, YGR197C, YPR159W, YDR093W, YDR038C, YLR229C, YGR281W, YLR020C, YNL231C, YFL005W, YKL203C, YLL061W, YDR497C, YPL274W, YCR028C, YGL084C, YCR004C, YER166W, YHR073W, YML125C, YFL050C, YLR411W, YGR065C, YMR215W, YBR086C'
membrane0 = membrane.split(', ')
glycolysis = 'YJR009C, YGR192C, YDR050C, YJL052W, YGR254W, YCR012W, YCL040W, YBR196C, YAL038W, YHR174W'
glycolysis0 = glycolysis.split(', ')

compartment1['plasma_sgd_2026_refine'] = membrane0
compartment1['plasma_sgd_2026_refine_reduce_glycolysis'] = list(set(list(set(gene_sgd_2026['gene'].tolist())))-set(glycolysis0))
## 结论 之前plasma中包括10个糖酵解蛋白，去掉之后mass fraction降低50%。





# check the latest annotation
# 文件名为附件中的名称
filename = 'data/sce_compartment_curation/2026/mitochondrion_annotations_2026.txt'
df = pd.read_csv(filename, sep='\t', skiprows=8, header=0)
df00 = df[df['Evidence']=='BSR']
print(set(df['Evidence'].tolist())) # 其中包括 ISM
df00 = df[df['Evidence'].str.contains('ISM')]

# only select "located"
df = df[df['Qualifier']=='located in']
df = df[~df['Systematic Name/Complex Accession'].str.contains('CPX-')]
gene_sgd_2026 = df.iloc[:,0:2]
gene_sgd_2026.columns =['short_name','gene']

compartment1['mitochondrion_sgd_2026'] = list(set(gene_sgd_2026['gene'].tolist()))
','.join(list(set(gene_sgd_2026['gene'].tolist())))
# enrichment analysis, find 27 genes related to EMP, only 8 genes mainly exist in mitochrondria
glycolysis0 = non_mitochondrial_proteins = [
    'YFR053C',  # HXK1
    'YDR050C',  # TPI1
    'YJL052W',  # TKL1
    'YGR240C',  # (标准名称通常为TKL1，同上？原列表中YJL052W和YGR240C均为TKL相关)
    'YCR012W',  # PGK1
    'YKL060C',  # FBA1
    'YBR196C',  # PGI1
    'YMR205C',  # PFK2
    'YJR009C',  # TDH2
    'YKL152C',  # GPM1
    'YAL054C',  # ACS1
    'YGR254W',  # ENO1
    'YHR174W',  # ENO2
    'YGR192C',  # TDH3
    'YPL061W',  # (通常为ALD6或相关)
    'YER073W',  # ALD3
    'YGL253W',  # HXK2
    'YOR347C',  # PYK2
    'YDL168W'   # SFA1
]

compartment1['mitochondrion_sgd_2026_reduce_glycolysis'] = list(set(list(set(gene_sgd_2026['gene'].tolist())))-set(glycolysis0))
','.join(list(set(list(set(gene_sgd_2026['gene'].tolist())))-set(glycolysis0)))

# 建议从线粒体质量计算中剔除的蛋白列表
to_exclude_proteins = [
    # --- 明确的非线粒体污染物 (主要定位为胞质、细胞核或细胞壁) ---
    "YLL024C",  # SSA2: 高丰度胞质 Hsp70，极易粘附在线粒体表面
    "YBR072W",  # HSP26: 胞质小热休克蛋白
    "YLR043C",  # TRX1: 胞质硫氧还蛋白
    "YGL256W",  # ADH4: 胞质乙醇脱氢酶 (注意：不同于线粒体的 ADH3)
    "YAL035W",  # FUN12: 胞质翻译起始因子
    "YAL019W",  # FUN30: 细胞核染色质重塑因子
    "YAL001C",  # TFC3: 细胞核转录因子
    "YAL011W",  # SWC3: 细胞核 SWR1 复合物亚基
    "YLR390W-A",  # CCW14: 细胞壁糖蛋白 (明确的胞外组分)
    "YOR151C",  # YHB1: 胞质氧化氮双加氧酶

    # --- 数据库注释噪音 (Dubious ORFs: 大多不编码功能蛋白，或为重叠序列) ---
    "YNCQ0020W", "YNCQ0024C", "YNCL0009C", "YNCM0036C", "YNCQ0009W",
    "YNCK0007W", "YNCF0014C", "YNCE0014W", "YNCC0009C", "YNCD0029C",
    "YNCQ0015W", "YNCI0007C", "YNCE0022C", "YNCM0036C", "YNCQ0001W",
    "YNCH0005W", "YNCQ0004W", "YNCG0041W", "YNCQ0018W", "YNCQ0012W",
    "YNCD0016C", "YNCQ0013W", "YNCQ0027W", "YNCQ0016W", "YNCQ0019W",
    "YNCG0007C", "YNCC0010C", "YNCQ0003W", "YNCE0009C", "YNCQ0007W",
    "YNCQ0005W", "YNCQ0011W", "YNCP0010W", "YNCD0009W", "YNCQ0017W",
    "YNCQ0008W", "YNCQ0026W", "YNCQ0025W", "YNCJ0018W", "YNCQ0022W",
    "YNCQ0014W", "YNCD0014C", "YNCQ0021W", "YNCQ0023W", "YNCG0005W",
    "YNCQ0010W"
]
compartment1['mitochondrion_sgd_2026_reduce_glycolysis_refine'] = list(set(compartment1['mitochondrion_sgd_2026_reduce_glycolysis'])-set(to_exclude_proteins))





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

out00 = ProMassRatio_Organelle(protein_abundance=mass_fraction_final)
out00.to_excel("data/sce_compartment_curation/mitochondrion_fraction_under_different_input.xlsx")


