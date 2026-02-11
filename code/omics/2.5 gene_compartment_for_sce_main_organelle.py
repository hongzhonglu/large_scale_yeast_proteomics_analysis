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


### output ###
df_update = df[df['Systematic Name/Complex Accession'].isin(compartment1['mitochondrion_sgd_2026_reduce_glycolysis_refine'])]   # update the name
df123 = df_update[df_update['Systematic Name/Complex Accession'].isin(sce_gene['GeneName'])]
print("quality check")
print(list(set(df['Systematic Name/Complex Accession'].tolist())-set(sce_gene['GeneName'].tolist())))
df_renamed = df123.rename(columns={'Systematic Name/Complex Accession': 'gene'})
df_renamed.to_excel("data/sce_compartment_curation/2026_curated/mitochondrion_annotations.xlsx") # update the name









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
glycolysis = 'YJR009C, YGR192C, YDR050C, YJL052W, YGR254W, YCR012W, YCL040W, YBR196C, YAL038W, YHR174W'
glycolysis0 = glycolysis.split(', ')

compartment1['plasma_sgd_2026_refine_reduce_glycolysis'] = list(set(list(set(gene_sgd_2026['gene'].tolist())))-set(glycolysis0))
## 结论 之前plasma中包括10个糖酵解蛋白，去掉之后mass fraction降低50%。

','.join(compartment1['plasma_sgd_2026_refine_reduce_glycolysis'])

# go term -- membrane, further filtering based on membrane annotation
membrane = 'YDR055W, YPR192W, YKL046C, YFL026W, YJL170C, YOR273C, YHL016C, YNL294C, YNL173C, YML132W, YBR078W, YNR002C, YER060W, YNL323W, YOR381W, YPR149W, YDR032C, YOR030W, YBR054W, YLR092W, YBR043C, YIR019C, YML047C, YDR459C, YIL140W, YJL171C, YGR031C-A, YGR224W, YDL135C, YKR105C, YKL209C, YLL052C, YLR025W, YOL122C, YPR124W, YCR024C-A, YLR081W, YIR006C, YMR031C, YOL158C, YJL212C, YHR092C, YGR213C, YBR068C, YGL053W, YDR033W, YNL098C, YLR120C, YHR048W, YER123W, YML013W, YGL077C, YGL186C, YDL035C, YCR037C, YJL062W, YBR294W, YDR276C, YBL069W, YOL002C, YOR348C, YKR106W, YNL194C, YDR522C, YKL178C, YOR161C, YLR121C, YMR068W, YDL012C, YJL129C, YGL051W, YGL008C, YHR135C, YDR384C, YOR008C, YBR295W, YBR296C, YPR194C, YHR094C, YNL291C, YOR390W, YGR138C, YGR060W, YJR152W, YKL220C, YJL093C, YLR194C, YER185W, YDR040C, YER118C, YOR306C, YAL030W, YJL145W, YMR034C, YIL118W, YHL040C, YBR132C, YLR096W, YDR160W, YCL073C, YAR033W, YOR047C, YCL027W, YLL043W, YKR039W, YJR054W, YGL208W, YJR066W, YGR014W, YCR010C, YMR058W, YIL105C, YIL088C, YKL094W, YOR317W, YLR332W, YLR138W, YPR201W, YDL194W, YPR075C, YDL138W, YGR191W, YCR021C, YOL103W, YGR023W, YBL042C, YPR198W, YKR050W, YMR011W, YMR319C, YOR328W, YMR307W, YOR153W, YLR353W, YGR121C, YBR069C, YNR047W, YGR266W, YNL271C, YNL283C, YJR040W, YDL222C, YOL019W, YLL010C, YLR342W, YOR327C, YNL065W, YBR021W, YJL100W, YDR420W, YOR071C, YOR378W, YER020W, YML116W, YPL279C, YPL092W, YLR214W, YCL025C, YDL161W, YDR343C, YOL105C, YDR463W, YJR086W, YNL047C, YPL232W, YFR029W, YGR241C, YDR342C, YHR161C, YGL255W, YDR039C, YER056C, YJL058C, YOL009C, YOL078W, YLR310C, YGR055W, YLR019W, YGR217W, YPL265W, YER060W-A, YHR186C, YDL019C, YML052W, YDR345C, YJL156C, YAR031W, YNL154C, YOL152W, YDR210W, YIL147C, YHL044W, YPL036W, YPR165W, YPR032W, YKR093W, YOR011W, YKL126W, YNL142W, YNL275W, YLR237W, YBL029C-A, YMR017W, YOL020W, YDR011W, YIL121W, YLR130C, YGR041W, YCR075C, YGR198W, YOL109W, YOR104W, YNL192W, YDR536W, YMR183C, YOR171C, YOL130W, YGR260W, YOR086C, YPR156C, YCR098C, YLR373C, YIL120W, YIR038C, YKR055W, YDR034W-B, YDR208W, YDR090C, YER120W, YNL257C, YOR018W, YBR008C, YBR129C, YGR152C, YNR060W, YER145C, YKL217W, YJL198W, YLR413W, YPL058C, YEL063C, YHR096C, YGL114W, YLL028W, YNL180C, YML123C, YLR414C, YMR008C, YMR063W, YDR122W, YCR017C, YIL047C, YKL051W, YER143W, YBR016W, YHR005C, YNL268W, YOR101W, YMR238W, YGR197C, YPR159W, YDR093W, YDR038C, YLR229C, YGR281W, YLR020C, YNL231C, YFL005W, YKL203C, YLL061W, YDR497C, YPL274W, YCR028C, YGL084C, YCR004C, YER166W, YHR073W, YML125C, YFL050C, YLR411W, YGR065C, YMR215W, YBR086C'
membrane_gene_list = membrane.split(', ')

### output ###
df_update = df[df['Systematic Name/Complex Accession'].isin(compartment1['plasma_sgd_2026_refine_reduce_glycolysis'])]   # update the name
df1234 = df_update[df_update['Systematic Name/Complex Accession'].isin(sce_gene['GeneName'])]
df123 = df1234[df1234['Systematic Name/Complex Accession'].isin(membrane_gene_list)]

print("quality check")
print(list(set(df['Systematic Name/Complex Accession'].tolist())-set(sce_gene['GeneName'].tolist())))
df_renamed = df123.rename(columns={'Systematic Name/Complex Accession': 'gene'})
df_renamed.to_excel("data/sce_compartment_curation/2026_curated/plasma_annotations.xlsx") # update the name









# check the latest annotation
# 文件名为附件中的名称
filename = 'data/sce_compartment_curation/2026/cytosol_annotations_2026.txt'
df = pd.read_csv(filename, sep='\t', skiprows=8, header=0)
df00 = df[df['Evidence']=='BSR']
print(set(df['Evidence'].tolist())) # 其中包括 ISM
df00 = df[df['Evidence'].str.contains('ISS')]

# only select "located"
df = df[df['Qualifier']=='located in']
df = df[~df['Systematic Name/Complex Accession'].str.contains('CPX-')]
df = df[~df['Systematic Name/Complex Accession'].str.contains('YNC')]

gene_sgd_2026 = df.iloc[:,0:2]
gene_sgd_2026.columns =['short_name','gene']
compartment1['cytosol_sgd_2026'] = list(set(gene_sgd_2026['gene'].tolist()))
','.join(list(set(gene_sgd_2026['gene'].tolist())))

#cytosol_double_check = list(set(compartment1['cytosol_sgd_2026']) & set(ribo_gene_double_check))
#根据这个对比，cytosol中多了核糖体蛋白 YER117W
# 根据查阅，cytosol中少了一个关键酶YCR012W
supplementary_cytosol_proteins = [
    "YOL086C", # ADH1: 乙醇脱氢酶1，酵母丰度最高的蛋白之一
    "YAL038W", # PYK1 (CDC19): 丙酮酸激酶，糖酵解最后一步的关键质量贡献者
    "YDR050C", # TPI1: 磷酸丙糖异构酶，极高丰度代谢酶
    "YBR196C", # PGI1: 磷酸葡萄糖异构酶，糖酵解途径核心成员
    "YBR118W", # TEF2: 翻译延伸因子 EF-1 alpha，与 TEF1 共同占据巨大质量比
    "YCR012W" ] #3-磷酸甘油酸激酶
compartment1['cytosol'] = list(set(gene_sgd_2026['gene'].tolist())) + supplementary_cytosol_proteins
compartment1['cytosol'] = [x for x in compartment1['cytosol'] if x !='YER117W']

### output ###
df_update = df[df['Systematic Name/Complex Accession'].isin(compartment1['cytosol'])]   # update the name
df123 = df_update[df_update['Systematic Name/Complex Accession'].isin(sce_gene['GeneName'])]
print("quality check")
print(list(set(df['Systematic Name/Complex Accession'].tolist())-set(sce_gene['GeneName'].tolist())))
df_renamed = df123.rename(columns={'Systematic Name/Complex Accession': 'gene'})
df_renamed.to_excel("data/sce_compartment_curation/2026_curated/cytosol_annotations.xlsx") # update the name






# check the latest annotation
# 文件名为附件中的名称
filename = 'data/sce_compartment_curation/2026/lipid_droplet_annotations.txt'
df = pd.read_csv(filename, sep='\t', skiprows=8, header=0)
df00 = df[df['Evidence']=='BSR']
print(set(df['Evidence'].tolist())) # 其中包括 ISM
df00 = df[df['Evidence'].str.contains('ISS')]

# only select "located"
df = df[df['Qualifier']=='located in']
df = df[~df['Systematic Name/Complex Accession'].str.contains('CPX-')]
gene_sgd_2026 = df.iloc[:,0:2]
gene_sgd_2026.columns =['short_name','gene']
compartment1['lipid_droplet'] = list(set(gene_sgd_2026['gene'].tolist()))
','.join(list(set(gene_sgd_2026['gene'].tolist())))

# 不在脂滴（Lipid Droplet）中的酿酒酵母蛋白列表
non_lipid_droplet_proteins = [
    "YJR009C",  # 虽然注释有在脂滴 这个是糖酵解，应该从脂滴中移除
    "YJL052W",  # ASH1: 虽然注释有在脂滴 这个是糖酵解，应该从脂滴中移除
    "YGR192C",  # TDH3: 虽然注释有在脂滴 这个是糖酵解，应该从脂滴中移除
    "YDR196C",  # 线粒体 (辅酶A生物合成)，有脂滴注释，参与多个细胞器
]
compartment1['lipid_droplet'] = list(set(list(set(gene_sgd_2026['gene'].tolist())))-set(non_lipid_droplet_proteins))

### output ###
df_update = df[df['Systematic Name/Complex Accession'].isin(compartment1['lipid_droplet'])]   # update the name
df123 = df_update[df_update['Systematic Name/Complex Accession'].isin(sce_gene['GeneName'])]
print("quality check")
print(list(set(df['Systematic Name/Complex Accession'].tolist())-set(sce_gene['GeneName'].tolist())))
df_renamed = df123.rename(columns={'Systematic Name/Complex Accession': 'gene'})
df_renamed.to_excel("data/sce_compartment_curation/2026_curated/lipid_droplet_annotations.xlsx") # update the name







# 文件名为附件中的名称
filename = 'data/sce_compartment_curation/2026/fungaltype_vacuole_annotations.txt' # update the name
df = pd.read_csv(filename, sep='\t', skiprows=8, header=0)
df00 = df[df['Evidence']=='BSR']
print(set(df['Evidence'].tolist())) # 其中包括 ISM
df00 = df[df['Evidence'].str.contains('ISS')]

# only select "located"
df = df[df['Qualifier']=='located in']
df = df[~df['Systematic Name/Complex Accession'].str.contains('CPX-')]
gene_sgd_2026 = df.iloc[:,0:2]
gene_sgd_2026.columns =['short_name','gene']
compartment1['fungal_type_vacuole'] = list(set(gene_sgd_2026['gene'].tolist()))
','.join(list(set(gene_sgd_2026['gene'].tolist())))

### output ###
df_update = df[df['Systematic Name/Complex Accession'].isin(compartment1['fungal_type_vacuole'])]   # update the name
df123 = df_update[df_update['Systematic Name/Complex Accession'].isin(sce_gene['GeneName'])]
print("quality check")
print(list(set(df['Systematic Name/Complex Accession'].tolist())-set(sce_gene['GeneName'].tolist())))
df_renamed = df123.rename(columns={'Systematic Name/Complex Accession': 'gene'})
df_renamed.to_excel("data/sce_compartment_curation/2026_curated/fungal_type_vacuole_annotations.xlsx") # update the name





# 文件名为附件中的名称
filename = 'data/sce_compartment_curation/2026/peroxisome_annotations.txt' # update the name
df = pd.read_csv(filename, sep='\t', skiprows=8, header=0)
df00 = df[df['Evidence']=='BSR']
print(set(df['Evidence'].tolist())) # 其中包括 ISM
df00 = df[df['Evidence'].str.contains('ISS')]

# only select "located"
df = df[df['Qualifier']=='located in']
df = df[~df['Systematic Name/Complex Accession'].str.contains('CPX-')]
gene_sgd_2026 = df.iloc[:,0:2]
gene_sgd_2026.columns =['short_name','gene']
compartment1['peroxisome'] = list(set(gene_sgd_2026['gene'].tolist()))
','.join(list(set(gene_sgd_2026['gene'].tolist())))

### output ###
df_update = df[df['Systematic Name/Complex Accession'].isin(compartment1['peroxisome'])]   # update the name
df123 = df_update[df_update['Systematic Name/Complex Accession'].isin(sce_gene['GeneName'])]
print("quality check")
print(list(set(df['Systematic Name/Complex Accession'].tolist())-set(sce_gene['GeneName'].tolist())))
df_renamed = df123.rename(columns={'Systematic Name/Complex Accession': 'gene'})
df_renamed.to_excel("data/sce_compartment_curation/2026_curated/peroxisome_annotations.xlsx") # update the name








# 文件名为附件中的名称
filename = 'data/sce_compartment_curation/2026/nucleolus_annotations.txt'
df = pd.read_csv(filename, sep='\t', skiprows=8, header=0)
df00 = df[df['Evidence']=='BSR']
print(set(df['Evidence'].tolist())) # 其中包括 ISM
df00 = df[df['Evidence'].str.contains('ISS')]

# only select "located"
df = df[df['Qualifier']=='located in']
df = df[~df['Systematic Name/Complex Accession'].str.contains('CPX-')]
gene_sgd_2026 = df.iloc[:,0:2]
gene_sgd_2026.columns =['short_name','gene']
compartment1['nucleolus'] = list(set(gene_sgd_2026['gene'].tolist()))
','.join(list(set(gene_sgd_2026['gene'].tolist())))

### output ###
df_update = df[df['Systematic Name/Complex Accession'].isin(compartment1['nucleolus'])]   # update the name
df123 = df_update[df_update['Systematic Name/Complex Accession'].isin(sce_gene['GeneName'])]
print("quality check")
print(list(set(df['Systematic Name/Complex Accession'].tolist())-set(sce_gene['GeneName'].tolist())))
df_renamed = df123.rename(columns={'Systematic Name/Complex Accession': 'gene'})
df_renamed.to_excel("data/sce_compartment_curation/2026_curated/nucleolus_annotations.xlsx") # update the name



#compartment1['m_latest'] = gene_list300['gene'].tolist() # from part 2.7
#compartment1['m_latest2'] = gene_double_check # from part 2.7

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










# other batch process
# combine
def polish_annotaiton(df):
    sce_gene = pd.read_excel("data/uniprotGeneID_mapping.xlsx")
    df = df[df['Qualifier']=='located in'] # only select "located"
    df = df[~df['Systematic Name/Complex Accession'].str.contains('CPX-')]
    # remove the protein with IEA evidence. IEA (Electronic Annotation)	最低	仅基于序列相似性预测，未经验证。
    df = df[~df['Evidence'].isin(['IEA'])]
    df = df[df['Systematic Name/Complex Accession'].isin(sce_gene['GeneName'])]
    df_renamed = df.rename(columns={'Systematic Name/Complex Accession': 'gene'})
    return df_renamed

fungaltype_cell_wall = pd.read_csv("data/sce_compartment_curation/2026/fungaltype_cell_wall_annotations.txt",sep='\t', skiprows=8, header=0)
fungaltype_cell_wall = polish_annotaiton(fungaltype_cell_wall)
fungaltype_cell_wall.to_excel('data/sce_compartment_curation/2026_curated/fungal_type_cell_wall_annotations.xlsx')


fungal_type_vacuole_membrane = pd.read_csv("data/sce_compartment_curation/2026/fungaltype_vacuole_membrane_annotations.txt",sep='\t', skiprows=8, header=0)
fungal_type_vacuole_membrane = polish_annotaiton(fungal_type_vacuole_membrane)
fungal_type_vacuole_membrane.to_excel('data/sce_compartment_curation/2026_curated/fungal_type_vacuole_membrane_annotations.xlsx')


cytoplasm = pd.read_csv("data/sce_compartment_curation/2026/cytoplasm_annotations.txt",sep='\t', skiprows=8, header=0)
cytoplasm = polish_annotaiton(cytoplasm)
cytoplasm.to_excel('data/sce_compartment_curation/2026_curated/cytoplasm_annotations.xlsx')







nucleus = pd.read_csv("data/sce_compartment_curation/2026/nucleus_annotations.txt",sep='\t', skiprows=8, header=0)
nucleus = polish_annotaiton(nucleus)
nucleus.to_excel('data/sce_compartment_curation/2026_curated/nucleus_annotations.xlsx')
gene_test = list(set(nucleus['gene'].tolist()))
','.join(gene_test)

# more check
nucleus_gene_analysis = pd.read_excel("data/sce_compartment_curation/2026/nucleus_gene_function_analysis.xlsx", sheet_name='double_check')
nucleus_metabolic_gene = nucleus_gene_analysis['gene'].tolist()
nucleus_m = ', '.join(nucleus_metabolic_gene)
nucleus_m = nucleus_m.split(', ')
nucleus_m = list(set(nucleus_m))
','.join(nucleus_m)
# 114个酵母蛋白的最终精确定位分类

# ==========================================
# 组 1: 细胞核相关蛋白 (Nucleus Group)
# 包含：转录因子、常驻核蛋白、核内组装蛋白、核/质双定位且在核内有功能的蛋白
# ==========================================
nucleus_proteins = [
    "YOL108C", # INO1: 您指出的核外周/转录相关定位
    "YOL052C", # SEM1: 核孔复合物关联
    "YJR090C", # GRR1: 争议6蛋白之一 (核内降解底物)
    "YLR153C", # ACS2: 核内乙酰CoA合成
    "YGL253W", # HXK2: 葡萄糖感应/核内转录调控
    "YGR233C", # PHO81: 磷酸盐信号传导 (核内)
    "YGR180C", # RNR4: RNR小亚基 (核常驻)
    "YMR121C", # RPL14B: 核糖体组装 (核仁)
    "YOL032W", # RPP2B: 核糖体组装 (核仁)
    "YNL123W", # NMA111: 核内蛋白酶
    "YJR049C", # TRM1: tRNA甲基转移酶 (核内)
    "YDR123C", # INO2: 转录因子
    "YGL037C", # PNC1: 争议6蛋白之一 (压力下入核)
    "YHL020C", # OPI1: 争议6蛋白之一 (转录抑制因子)
    "YMR226C", # TMA17/TMA22: 争议6蛋白之一 (核内组装相关)
    "YKL069W", # MET18: 核内Fe-S簇组装
    "YOR025W", # HST3: 组蛋白去乙酰化
    "YER117W", # RNR3: 核糖核苷酸还原酶 (核内)
    "YDL131W", # LYS21: 兼职核内染色质调控
    "YMR009W", # ADI1: MTA循环 (核内)
    "YGL103W", # RPL28: 核糖体组装 (核仁)
    "YDL182W", # LYS20: 兼职核内染色质调控
    "YJL026W", # RNR2: RNR小亚基 (核常驻)
    "YNL267W", # PIK1: 磷脂酰肌醇4-激酶 (核内池)
    "YOR047C", # STD1: 转录调节
    "YGR010W", # NMA2: NAD+合成 (核)
    "YLR260W", # LCB5: 鞘氨醇激酶 (核/质双定位)
    "YKL038W", # RGT1: 转录因子
    "YBL105C", # PKC1: 争议6蛋白之一 (穿梭激酶)
    "YMR104C", # YPK2: 争议6蛋白之一 (信号穿梭)
    "YOR209C", # NPT1: NAD+合成 (核)
    "YDR315C", # IPK1: 肌醇激酶 (核)
    "YLR328W", # NMA1: NAD+合成 (核)
    "YPL015C"  # HST1: 组蛋白去乙酰化
]

# ==========================================
# 组 2: 非细胞核蛋白 (Non-Nucleus Group)
# 包含：严格定位在细胞质、线粒体、内质网、液泡的蛋白
# ==========================================
non_nucleus_proteins = [
    "YPL111W", "YHR163W", "YDR191W", "YIL083C", "YER003C", "YDR115W",
    "YKR031C", "YIL148W", "YML126C", "YPL273W", "YDR111C", "YPL117C",
    "YGL144C", "YGR043C", "YAL062W", "YPR035W", "YFR047C", "YDL103C",
    "YNR034W", "YGR202C", "YEL058W", "YFL017C", "YGR205W", "YDR173C",
    "YGL026C", "YEL038W", "YHR074W", "YBR117C", "YAL012W", "YMR220W",
    "YJR130C", "YHR137W", "YPR069C", "YDR208W", "YER070W", "YCL047C",
    "YGL184C", "YDR158W", "YCL026C-A", "YHR043C", "YGR007W", "YOR375C",
    "YDR035W", "YMR208W", "YDR354W", "YDL205C", "YPL028W", "YHR046C",
    "YGL040C", "YPL091W", "YJR148W", "YNL045W", "YOR095C", "YBR249C",
    "YNL267W", "YLR354C", "YKR069W", "YOR323C", "YGR277C", "YER163C",
    "YNL036W", "YDR400W", "YJR009C", "YLR118C", "YML082W", "YCR053W",
    "YDR047W", "YCR036W", "YDR009W", "YPR060C", "YIL020C", "YIL145C",
    "YGR248W", "YGR208W", "YDR018C", "YNL012W", "YPR118W", "YGR192C",
    "YDR531W", "YJR139C", "YMR278W"
]

# 校验
print(f"细胞核蛋白组: {len(set(nucleus_proteins))} 个")
print(f"非细胞核蛋白组: {len(set(non_nucleus_proteins))} 个")
','.join(non_nucleus_proteins)
nucleus_refine = nucleus[~nucleus['gene'].isin(non_nucleus_proteins)]
nucleus_refine.to_excel('data/sce_compartment_curation/2026_curated/nucleus_annotations.xlsx')



# other
endosome = pd.read_csv("data/sce_compartment_curation/2026/endosome_annotations.txt",sep='\t', skiprows=8, header=0)
endosome = polish_annotaiton(endosome)
endosome.to_excel('data/sce_compartment_curation/2026_curated/endosome_annotations.xlsx')


# other
endoplasmic_reticulum = pd.read_csv("data/sce_compartment_curation/2026/endoplasmic_reticulum_annotations.txt",sep='\t', skiprows=8, header=0)
endoplasmic_reticulum = polish_annotaiton(endoplasmic_reticulum)
endoplasmic_reticulum.to_excel('data/sce_compartment_curation/2026_curated/endoplasmic_reticulum_annotations.xlsx')

# other
endoplasmic_reticulum_membrane = pd.read_csv("data/sce_compartment_curation/2026/endoplasmic_reticulum_membrane_annotations.txt",sep='\t', skiprows=8, header=0)
endoplasmic_reticulum_membrane = polish_annotaiton(endoplasmic_reticulum_membrane)
endoplasmic_reticulum_membrane.to_excel('data/sce_compartment_curation/2026_curated/endoplasmic_reticulum_membrane_annotations.xlsx')

# other
Golgi_apparatus = pd.read_csv("data/sce_compartment_curation/2026/Golgi_apparatus_annotations.txt",sep='\t', skiprows=8, header=0)
Golgi_apparatus = polish_annotaiton(Golgi_apparatus)
Golgi_apparatus.to_excel('data/sce_compartment_curation/2026_curated/Golgi_apparatus_annotations.xlsx')

# other
Golgi_membrane = pd.read_csv("data/sce_compartment_curation/2026/Golgi_membrane_annotations.txt",sep='\t', skiprows=8, header=0)
Golgi_membrane = polish_annotaiton(Golgi_membrane)
Golgi_membrane.to_excel('data/sce_compartment_curation/2026_curated/Golgi_membrane_annotations.xlsx')