# note:
# this part can be rewritten as a function
# we evaluate the reliablity of protein location annotation based on several criteria.
# 1. go evidence ?  EXP, IDA, IMP, IGI, IPI 实验证据;   HDA 高通量实验
# 2. how many compartments one protein was assigned to?
import pandas as pd
import numpy as np


organelle_hierarchy = {
    "mitochondrion": [
        "extrinsic component of mitochondrial inner membrane",
        "mitochondria-associated endoplasmic reticulum membrane contact site",
        "mitochondrial crista junction",
        "mitochondrial envelope",
        "mitochondrial inner membrane",
        "mitochondrial intermembrane space",
        "mitochondrial matrix",
        "mitochondrial membrane",
        "mitochondrial nucleoid",
        "mitochondrial outer membrane",
        "mitochondrial ribosome",
        "mitochondrion",
        "vacuole-mitochondrion membrane contact site"
    ],

    "nucleus": [
        "chromatin",
        "chromosome",
        "chromosome, centromeric region",
        "chromosome, telomeric region",
        "condensed chromosome, centromeric region",
        "condensed nuclear chromosome",
        "euchromatin",
        "kinetochore",
        "nuclear chromosome",
        "nuclear envelope",
        "nuclear inner membrane",
        "nuclear lumen",
        "nuclear membrane",
        "nuclear outer membrane",
        "nuclear periphery",
        "nuclear replication fork",
        "nucleolus",
        "nucleoplasm",
        "nucleus",
        "nucleus-vacuole junction",
        "perinuclear endoplasmic reticulum",
        "perinuclear region of cytoplasm",
        "rDNA heterochromatin"
    ],

    "endoplasmic reticulum": [
        "COPII-coated ER to Golgi transport vesicle",
        "cortical endoplasmic reticulum",
        "endoplasmic reticulum",
        "endoplasmic reticulum lumen",
        "endoplasmic reticulum membrane",
        "endoplasmic reticulum tubular network",
        "ER to Golgi transport vesicle membrane",
        "mitochondria-associated endoplasmic reticulum membrane contact site",
        "nuclear outer membrane-endoplasmic reticulum membrane network",
        "perinuclear endoplasmic reticulum",
        "rough endoplasmic reticulum membrane"
    ],

    "Golgi apparatus": [
        "cis-Golgi network",
        "Golgi apparatus",
        "Golgi cis cisterna",
        "Golgi medial cisterna",
        "Golgi membrane",
        "Golgi trans cisterna",
        "Golgi-associated vesicle",
        "trans-Golgi network"
    ],

    "fungal-type vacuole": [
        "fungal-type vacuole",
        "fungal-type vacuole lumen",
        "fungal-type vacuole membrane",
        "nucleus-vacuole junction",
        "vacuole",
        "vacuole-isolation membrane contact site",
        "vacuolar lumen",
        "vacuolar membrane"
    ],

    "peroxisome": [
        "peroxisomal matrix",
        "peroxisomal membrane",
        "peroxisome"
    ],

    "endosome": [
        "early endosome",
        "early endosome membrane",
        "endosome",
        "endosome membrane",
        "late endosome",
        "late endosome membrane"
    ],

    "lipid droplet": [
        "lipid droplet"
    ],

    "plasma membrane": [
        "cell cortex",
        "cell periphery",
        "cell tip",
        "mating projection",
        "mating projection base",
        "mating projection tip",
        "membrane",
        "membrane raft",
        "plasma membrane",
        "site of polarized growth"
    ],

    "P-body": [
        "P-body"
    ],

    "cytoplasmic stress granule": [
        "cytoplasmic stress granule"
    ],

    "spindle pole body": [
        "spindle pole body"
    ],

    "ribosome": [
        "cytosolic ribosome",
        #"mitochondrial ribosome",
        "ribosome"
    ],

    "cytosol": [
        #"cytoskeleton",
        "cytosol" #,
        #"microtubule"
    ],

    "extracellular region": [
        "ascospore wall",
        "cell wall-bounded periplasmic space",
        "extracellular region",
        "fungal-type cell wall",
        "periplasmic space",
        "prospore membrane",
        "spore wall"
    ],

    "fungal-type cell wall": [
        "fungal-type cell wall"
    ]
}


def classify_localization_with_strict_propagation(df, organelle_hierarchy):
    # 1. 权重定义
    experimental_codes = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI'}
    evidence_weights = {
        'EXP': 10, 'IDA': 10, 'IPI': 10, 'IMP': 10, 'IGI': 10,
        'HDA': 7, 'TAS': 5, 'IC': 5,
        'IBA': 3, 'ISO': 3, 'NAS': 3, 'ISS': 3, 'ISA': 3, 'ISM': 3,
        'IEA': 1, 'ND': 0
    }

    # 预处理
    working_df = df[['gene', 'organelle', 'evidence']].copy()
    working_df['weight'] = working_df['evidence'].map(evidence_weights).fillna(1)

    # 2. 映射族群（Group 即为主细胞器名称）
    def get_group(loc):
        loc_lower = str(loc).lower()
        for group_name, list_of_locs in organelle_hierarchy.items():
            if loc_lower in [l.lower() for l in list_of_locs]:
                return group_name
        return "other_structures"

    working_df['group'] = working_df['organelle'].apply(get_group)

    # 3. 基础得分计算 (保留特异性惩罚)
    gene_effective_loc_count = working_df.groupby('gene')['group'].nunique().to_dict()

    scored_df = working_df.groupby(['gene', 'organelle', 'group']).agg(
        max_ev_weight=('weight', 'max'),
        ev_count=('weight', 'count'),
        has_experimental=('evidence', lambda x: any(e in experimental_codes for e in x))
    ).reset_index()

    # 基础分逻辑
    scored_df['raw_score'] = scored_df['max_ev_weight'] + np.minimum((scored_df['ev_count'] - 1) * 0.2, 2.0)

    # 族群层级奖励 (同一族内有多个描述则加分)
    group_member_count = scored_df.groupby(['gene', 'group'])['organelle'].transform('count')
    scored_df.loc[group_member_count > 1, 'raw_score'] += 1.5

    # 计算可靠性分数
    scored_df['final_score'] = scored_df.apply(
        lambda x: x['raw_score'] / np.sqrt(gene_effective_loc_count[x['gene']]), axis=1
    )

    s_min, s_max = scored_df['final_score'].min(), scored_df['final_score'].max()
    scored_df['reliability_score'] = ((scored_df['final_score'] - s_min) / (s_max - s_min) * 10).round(2)

    # 4. 初步分类逻辑
    def initial_assignment(group_df):
        exp_group_count = group_df[group_df['has_experimental']]['group'].nunique()
        max_score = group_df['reliability_score'].max()
        results = []
        for idx, row in group_df.iterrows():
            loc_type = "Secondary/Transient"
            note = []
            if row['has_experimental']:
                if exp_group_count == 1:
                    loc_type = "Primary"
                    note.append("Confirmed by unique experimental evidence group.")
                else:
                    note.append(f"Experimental evidence in {exp_group_count} distinct groups.")
            elif row['reliability_score'] == max_score and row['reliability_score'] > 4:
                loc_type = "Primary"
                note.append("Dominant location by score.")
            elif row['reliability_score'] < 2:
                loc_type = "Unclear"
                note.append("Low reliability.")
            else:
                note.append("Secondary support.")
            results.append({'Localization_Type': loc_type, 'Note': " ".join(note)})
        return pd.DataFrame(results)

    meta = scored_df.groupby('gene', group_keys=False).apply(initial_assignment)
    scored_df = pd.concat([scored_df, meta.reset_index(drop=True)], axis=1)

    # 5. 核心修正：精准向上追溯 (Strict Propagation)
    def strict_propagate(gene_df):
        # 1. 识别该基因下哪些“组(主细胞器)”已经拥有 Primary 级别的子结构
        # 注意：排除掉主细胞器本身，只看真正的亚细胞器证据
        primary_groups = set()
        for _, row in gene_df.iterrows():
            # 如果该行被判定为 Primary，且它的名字不等于 Group 键名（即它是亚结构）
            if row['Localization_Type'] == 'Primary' and row['organelle'].lower() != row['group'].lower():
                primary_groups.add(row['group'])

        # 2. 遍历该基因的所有行，仅提升“主细胞器”条目的状态
        for idx, row in gene_df.iterrows():
            # 逻辑：如果该行是对应 Primary 亚结构的主细胞器，且目前还不是 Primary
            if row['organelle'].lower() == row['group'].lower() and row['group'] in primary_groups:
                if gene_df.at[idx, 'Localization_Type'] != 'Primary':
                    gene_df.at[idx, 'Localization_Type'] = 'Primary'
                    gene_df.at[idx, 'Note'] += " (Promoted: Primary sub-structure evidence found in this group)"

        return gene_df

    final_result = scored_df.groupby('gene', group_keys=False).apply(strict_propagate)

    return final_result[['gene', 'organelle', 'reliability_score', 'Localization_Type', 'Note']].sort_values(
        by=['gene', 'reliability_score'], ascending=[True, False]
    )

# 使用方法：
filename = 'data/sce_compartment_curation/2026/alliancemine_results_2026-01-03T10-33-37.xlsx'
df = pd.read_excel(filename)
final_scores = classify_localization_with_strict_propagation (df, organelle_hierarchy=organelle_hierarchy)
final_scores_primary = final_scores[final_scores['Localization_Type']=='Primary']
print(final_scores.head(20))

gene = final_scores[final_scores['gene']=='YGR240C']


########## 核糖体蛋白重新校正
## 确保消除计算置信区间时cytosol出现的错误趋势。目前sgd数据库注释会把核糖体蛋白归为细胞质。
ribosome_gene = pd.read_excel("data/sce_compartment_curation/2026/ribosome_uniprot.xlsx")
ribo_genes = ribosome_gene['gene'].tolist()
ribo_genes = pd.Series(ribo_genes).dropna().tolist()
','.join(ribo_genes)
# 将ribo_genes 包括所有基因进行富集分析，得到注释为ribosome的基因，然后添加YMR242C, YOR312C，得到下述完整列表。
ribo_gene_double_check = 'YMR242C, YOR312C, YDL081C, YJR094W-A, YCR031C, YPL198W, YHR010W, YPR043W, YGL031C, YHL015W, YIL069C, YKL006W, YGL030W, YBR191W, YOL040C, YOL121C, YJL177W, YJL189W, YLR048W, YJR145C, YGL135W, YEL054C, YGL123W, YLR185W, YDR382W, YFL034C-A, YDR447C, YGR034W, YBL087C, YPL220W, YER056C-A, YPL081W, YLR340W, YDR418W, YKL156W, YLR388W, YGL147C, YBR031W, YPL249C-A, YLR325C, YLR406C, YDR450W, YKL180W, YBR048W, YLR249W, YDL184C, YLR167W, YLR264W, YHR203C, YIL148W, YDL133C-A, YPR102C, YOL039W, YML024W, YGL076C, YNL069C, YDL075W, YOR293W, YOR063W, YNL178W, YGR085C, YML063W, YMR142C, YKR094C, YHL033C, YLR075W, YDL061C, YDR471W, YJL190C, YMR194W, YER117W, YNL096C, YFR032C-A, YLR441C, YPR132W, YBL072C, YGR118W, YLR367W, YBL027W, YBR084C-A, YGL103W, YLR029C, YLR287C-A, YMR116C, YPL090C, YDR500C, YPL143W, YPL131W, YIL133C, YJR123W, YNL302C, YMR143W, YDL130W, YML026C, YIL052C, YLR344W, YNL067W, YOR167C, YLR333C, YHR141C, YOR369C, YLR061W, YGR027C, YGR148C, YDL083C, YOR182C, YDR025W, YBR181C, YOL127W, YPL079W, YNL162W, YBR189W, YLL045C, YNL301C, YER131W, YHL001W, YER074W, YDL191W, YOR096W, YHR021C, YDL082W, YML073C, YMR121C, YJL136C, YJL191W, YBL092W, YDR012W, YFR031C-A, YKR057W, YLR448W, YER102W, YGR214W, YMR230W, YOL120C, YGL189C, YOR234C, YIL018W, YDL136W, YDR064W'
ribo_gene_double_check = ribo_gene_double_check.split(', ')



target_organelle = 'cytosol'
# 定义一个掩码（Mask）
is_ribo_protein = final_scores['gene'].isin(ribo_gene_double_check)
is_cytosol = final_scores['organelle'] == target_organelle
final_scores.loc[is_ribo_protein & is_cytosol, 'organelle'] = 'ribosome'
########## 核糖体蛋白重新校正
final_scores.to_excel("data/sce_compartment_curation/organelle_score.xlsx")


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

gene_list = final_scores[final_scores['gene'].isin(non_mitochondrial_proteins)]
gene_list = gene_list[gene_list['Localization_Type']=='Primary']



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
gene_list2 = final_scores[final_scores['gene'].isin(gene_sgd_2026['gene'])]
gene_list3 = gene_list2[gene_list2['organelle']=='mitochondrion']
gene_list300 = gene_list3[gene_list3['Localization_Type']=='Primary']
gene_list4 = gene_list3[gene_list3['Localization_Type']!='Primary']
gene_list5 = gene_list4[gene_list4['Note'].str.contains('Experimental evidence')]
len(set(gene_list5['gene'].tolist()))

gene_double_check = list(set(gene_list3['gene'].tolist())-set(gene_list5['gene'].tolist()))

# check manual and automatic method
common_one = list(set(non_mitochondrial_proteins) & set(gene_list4['gene'].tolist()))
common_one2 = list(set(non_mitochondrial_proteins) & set(gene_list5['gene'].tolist()))




# 计算细胞器质量百分比
import pandas as pd
import numpy as np


def analyze_full_hierarchy_mass(score_path, mass_path, organelle_hierarchy):
    """
    计算主细胞器及子细胞器的蛋白质质量百分比及 95% CI
    """
    # 1. 加载数据 (本地运行时请确保安装了 openpyxl: pip install openpyxl)
    print("正在加载 Excel 数据...")
    # 如果是本地 Excel 文件，请取消注释下面两行
    score_df = pd.read_excel(score_path)
    mass_df = pd.read_excel(mass_path)

    # 为了演示，此处读取你上传的 CSV 版本
    #score_df = pd.read_csv(score_path)
    #mass_df = pd.read_csv(mass_path)

    # 2. 建立主亚结构映射
    rev_map = {sub.lower(): main for main, subs in organelle_hierarchy.items() for sub in subs}
    score_df['main_group'] = score_df['organelle'].str.lower().map(rev_map).fillna('Other')

    # 3. 计算子细胞器概率 P_sub
    # 计算每个蛋白在所有位置的总分
    gene_totals = score_df.groupby('gene')['reliability_score'].transform('sum')
    # 概率 P = 该子位置得分 / 该蛋白总分
    score_df['p_loc'] = np.where(gene_totals > 0, score_df['reliability_score'] / gene_totals, 0)

    # 4. 识别 275 个数据列
    data_cols = [c for c in mass_df.columns if c != 'gene']
    print(f"检测到 {len(data_cols)} 个数据点。正在启动高精细度计算...")

    # 5. 合并概率分布与质量数据
    # 我们保留 organelle 这一列作为子细胞器标识
    merged = score_df[['gene', 'main_group', 'organelle', 'p_loc']].merge(mass_df, on='gene', how='inner')

    final_results = []

    # 6. 遍历 275 列
    for col in data_cols:
        # 期望质量 (Signal): Mass_i * P_{i, sub}
        merged['weighted_m'] = merged[col] * merged['p_loc']
        # 方差 (Uncertainty): Mass_i^2 * P_{i, sub} * (1 - P_{i, sub})
        merged['var_m'] = (merged[col] ** 2) * merged['p_loc'] * (1 - merged['p_loc'])

        # 按主细胞器和子细胞器同时聚合
        summary = merged.groupby(['main_group', 'organelle']).agg(
            total_mass=('weighted_m', 'sum'),
            total_var=('var_m', 'sum')
        ).reset_index()

        # 7. 计算统计指标 (百分比)
        summary['std_err'] = np.sqrt(summary['total_var'])
        summary['Mass_Percentage'] = (summary['total_mass'] * 100).round(6)
        summary['Lower_95_CI'] = ((summary['total_mass'] - 1.96 * summary['std_err']) * 100).round(6)
        summary['Upper_95_CI'] = ((summary['total_mass'] + 1.96 * summary['std_err']) * 100).round(6)
        summary['Data_Point'] = col

        final_results.append(summary)

    # 8. 合并并导出
    full_output = pd.concat(final_results, ignore_index=True)

    # 调整列顺序，方便阅读
    cols_order = ['Data_Point', 'main_group', 'organelle', 'Mass_Percentage', 'Lower_95_CI', 'Upper_95_CI']
    full_output = full_output[cols_order]

    # 保存为 Excel (本地运行)
    full_output.to_excel("SubOrganelle_Mass_Fraction_Analysis.xlsx", index=False)
    print("分析完成！结果已准备就绪。")
    return full_output

# 运行 (请确保文件名正确)
# results = analyze_organelle_mass_full('organelle_score.xlsx', 'mass_fraction_combine.xlsx', organelle_hierarchy)

# 2. 实例化并运行
# 请确保文件名与你本地或上传的文件名一致
score_file = "data/sce_compartment_curation/organelle_score.xlsx"
mass_fraction_file = "data/proteomics/mass_fraction_combine.xlsx"
results = analyze_full_hierarchy_mass(score_file, mass_fraction_file, organelle_hierarchy)
#full_output.to_excel("SubOrganelle_Mass_Fraction_Analysis.xlsx", index=False)
# 查看前几个结果
print(results.head())

