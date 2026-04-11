# note:
# this part can be rewritten as a function

import os
import pandas as pd

# import self function
from src.protein_process import *


# input gene annotation in Uniprot
sce_gene = pd.read_excel("data/uniprotGeneID_mapping.xlsx")



# reanalyze the data set
gene_mitochondrion_sgd = pd.read_excel("data/sce_compartment_curation/mitochondrion_annotations_manual_v2.xlsx") # this is based sgd annotation



# input the NC paper
gene_mitochondrion_nc = pd.read_excel("data/sce_compartment_curation/mitochondrial_suborganelle.xlsx")
gene_mitochondrion_nc = gene_mitochondrion_nc.iloc[0:987:]


# input mitochrondrial protein list from  Absolute yeast mitochondrial proteome quantification reveals trade-off between biosynthesis and energy generation during diauxic shift
# input the data from PNAS
gene_mitochondrion_pnas = pd.read_excel("data/sce_compartment_curation/pnas.1918216117.sd01 VIP Francesca.xlsx", header=1)
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
output_f = {k: v for k, v in output.items() if 'mitochondrion' in v} #mitochondrion
# check which protein exist only in mitochondrion
output_f2 = {k: v for k, v in output_f.items() if len(v) <= 1}
compartment1['m_sgd_core'] = list(output_f2.keys())





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
# remove the protein with IEA evidence. IEA (Electronic Annotation)	最低	仅基于序列相似性预测，未经验证。
df = df[~df['Evidence'].isin(['IEA'])]


gene_sgd_2026 = df.iloc[:,0:2]
gene_sgd_2026.columns =['short_name','gene']

compartment1['mitochondrion_sgd_2026'] = list(set(gene_sgd_2026['gene'].tolist()))
','.join(list(set(gene_sgd_2026['gene'].tolist())))




# enrichment analysis, find 27 genes related to EMP, only 8 genes mainly exist in mitochondria.
# how to find another 9 proteins not in mitochondria.
non_mitochondrial_proteins = [
    'YFR053C',  # HXK1, glycolysis
    'YDR050C',  # TPI1, glycolysis
    'YJL052W',  # TKL1, glycolysis
    'YGR240C',  # (标准名称通常为TKL1，同上？原列表中YJL052W和YGR240C均为TKL相关), glycolysis
    'YCR012W',  # PGK1, glycolysis
    'YKL060C',  # FBA1, glycolysis
    'YBR196C',  # PGI1, glycolysis
    'YMR205C',  # PFK2, glycolysis
    'YJR009C',  # TDH2, glycolysis
    'YKL152C',  # GPM1, glycolysis
    'YAL054C',  # ACS1, glycolysis
    'YGR254W',  # ENO1, glycolysis
    'YHR174W',  # ENO2, glycolysis
    'YGR192C',  # TDH3, glycolysis
    'YPL061W',  # (通常为ALD6或相关), glycolysis
    'YGL253W',  # HXK2, glycolysis
    'YOR347C',  # PYK2, glycolysis
    'YDL168W',   # SFA1, glycolysis
    "YLR390W-A",   # CCW14: 细胞壁糖蛋白 (明确的胞外组分)   ok
    "YLL024C",  # SSA2: 高丰度胞质 Hsp70，极易粘附在线粒体表面
    "YBR072W",  # HSP26: 胞质小热休克蛋白
    "YLR043C",  # TRX1: 胞质硫氧还蛋白
    "YAL035W",  # FUN12: 胞质翻译起始因子
    "YAL019W",  # FUN30: 细胞核染色质重塑因子
    "YAL001C",  # TFC3: 细胞核转录因子
    "YAL011W",  # SWC3: 细胞核 SWR1 复合物亚基
    "YOR151C",  # YHB1: 胞质氧化氮双加氧酶
    "YGL008C",  # (PMA1): 质膜 H+-ATPase， 丰度极高
    "YDR342C",  # 葡萄糖转运
    "YDR343C",  # 葡萄糖转运
    "YDR233C"   # Reticulon protein; involved in nuclear pore assembly and maintenance of tubular ER morphology;
]

# 这些蛋白的线粒体注解多为 HDA（高通量），而非 IDA/IMP 等直接证据支持的主要定位。
compartment1['mitochondrion_sgd_2026_reduce_glycolysis_refine'] = list(set(list(set(gene_sgd_2026['gene'].tolist())))-set(non_mitochondrial_proteins))
compartment1['mitochondrion_intersection'] = list(set(compartment1['mitochondrion_sgd_2026_reduce_glycolysis_refine']) & set(compartment1['m_marcel']))


mass_fraction_final = pd.read_excel("data/proteomics/mass_fraction_combine.xlsx")
out00 = ProMassRatio_Organelle(protein_abundance=mass_fraction_final, compartment=compartment1)
out00.to_excel("data/sce_compartment_curation/mitochondrion_fraction_under_different_input.xlsx")