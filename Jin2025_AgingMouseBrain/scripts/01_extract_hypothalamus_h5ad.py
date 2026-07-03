"""
============================================================================
Phase 1: 从全脑h5ad提取全部下丘脑细胞
============================================================================
提取两个集合:
  1. anatomical_division (52,696 cells) - 聚类转录相似性分配
  2. region_of_interest (26,105 cells) - 解剖取材标签

添加完整注释: cell_type_major, age_group, sex_group, neuronal_subclass
============================================================================
"""
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import os, sys, gc
from pathlib import Path

# ===========================================================================
# 配置
# ===========================================================================
H5AD_LOG2 = "d:/decrepitude mouse hypothamulas/Zeng-Aging-Mouse-10Xv3-log2.h5ad"
BARCODE_ANAT = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_all_barcodes.txt"
BARCODE_DISS = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypothalamus_all_cells_dissection.csv"
CELLMETA_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/cell_metadata.csv"
CLUSTER_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/cluster.csv"
DONOR_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/donor.csv"
CROSSMAP_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/taxonomy/cell_cross_mapping_annotations.csv"
OUT_DIR = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus"

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print(" Phase 1: 提取下丘脑全部细胞")
print("=" * 60)

# ===========================================================================
# 1. 加载barcode列表
# ===========================================================================
print("\n[1] 加载barcode列表...")

with open(BARCODE_ANAT) as f:
    barcodes_anat = set(line.strip() for line in f if line.strip())
print(f"  anatomical barcodes: {len(barcodes_anat):,}")

barcodes_diss = set(pd.read_csv(BARCODE_DISS).iloc[:, 0].dropna().values)
print(f"  dissection barcodes: {len(barcodes_diss):,}")

# ===========================================================================
# 2. 加载全脑h5ad (backed mode)
# ===========================================================================
print("\n[2] 加载全脑h5ad (backed mode)...")
adata = sc.read_h5ad(H5AD_LOG2, backed='r')
print(f"  总细胞: {adata.n_obs:,}")
print(f"  总基因: {adata.n_vars:,}")

# 匹配barcodes
all_bc = set(adata.obs_names)
matched_anat = [b for b in barcodes_anat if b in all_bc]
matched_diss = [b for b in barcodes_diss if b in all_bc]
print(f"  anatomical 匹配: {len(matched_anat):,} / {len(barcodes_anat):,}")
print(f"  dissection 匹配: {len(matched_diss):,} / {len(barcodes_diss):,}")

# ===========================================================================
# 3. 加载并合并metadata
# ===========================================================================
print("\n[3] 加载metadata...")

# 3a. cell_metadata
cell_meta = pd.read_csv(CELLMETA_CSV)
cell_meta.index = cell_meta['cell_label']
print(f"  cell_metadata: {len(cell_meta):,} rows")

# 3b. cluster annotations
cluster_info = pd.read_csv(CLUSTER_CSV)
print(f"  cluster.csv: {len(cluster_info)} clusters")

# 3c. donor info
donor_info = pd.read_csv(DONOR_CSV)
donor_age_map = dict(zip(donor_info['donor_label'], donor_info['donor_age_category']))
donor_sex_map = dict(zip(donor_info['donor_label'], donor_info['donor_sex']))
print(f"  donor.csv: {len(donor_info)} donors")

# 3d. cross-mapping (has class/subclass/supertype)
try:
    cross_map = pd.read_csv(CROSSMAP_CSV)
    cross_map.index = cross_map['cell_label']
    has_crossmap = True
    print(f"  cross_mapping: {len(cross_map):,} rows")
except:
    has_crossmap = False
    print("  cross_mapping: NOT FOUND, will use cluster_name parsing")

# ===========================================================================
# 4. 定义cell_type_major分类函数
# ===========================================================================
def assign_cell_type_major(cluster_name):
    """根据cluster_name模式分配主要细胞类型"""
    if not isinstance(cluster_name, str):
        return "Unclassified"
    cn = cluster_name
    # 胶质/非神经元
    if "Microglia" in cn: return "Microglia"
    if "Astro" in cn and "ependymal" not in cn.lower(): return "Astrocyte"
    if "Tanycyte" in cn: return "Tanycyte"
    if "Ependymal" in cn or "ependymal" in cn: return "Ependymal"
    if "NFOL" in cn: return "NFOL"
    if "MFOL" in cn: return "MFOL"
    if "MOL" in cn: return "MOL"
    if "OPC" in cn: return "OPC"
    if "COP" in cn: return "COP"
    # 血管
    if "Endo" in cn and "Endo-" in cn: return "Endothelial"
    if "SMC" in cn or "SMC-" in cn: return "SMC"
    if "VLMC" in cn: return "VLMC"
    if "Peri" in cn or "Peri-" in cn: return "Pericyte"
    # 免疫
    if "BAM" in cn: return "BAM"
    if "_DC" in cn or "DC_" in cn: return "Dendritic_Cell"
    if "ABC" in cn: return "ABC"
    if "T cells" in cn or "Tcell" in cn: return "T_cell"
    # 神经元 - 按神经递质
    if "Glut" in cn and "Gaba" not in cn and "GABA" not in cn:
        return "Glutamatergic_Neuron"
    if "Gaba" in cn or "GABA" in cn:
        if "Dopa" in cn: return "GABAergic_Neuron"  # GABA-Dopa归入GABA
        if "Chol" in cn: return "GABAergic_Neuron"
        if "Hist" in cn: return "GABAergic_Neuron"
        return "GABAergic_Neuron"
    if "Dopa" in cn: return "Dopaminergic_Neuron"
    if "Sero" in cn: return "Serotonergic_Neuron"
    if "Chol" in cn: return "Cholinergic_Neuron"
    if "Hist" in cn: return "Histaminergic_Neuron"
    if "NN" in cn and any(k in cn for k in ["Glut", "Gaba"]):
        return "Neuron_Other"
    return "Unclassified"

def assign_neuronal_subclass(cluster_name):
    """识别关键下丘脑核团神经元亚类"""
    if not isinstance(cluster_name, str):
        return "Non_Neuronal"
    cn = cluster_name
    # 检查是否为神经元
    is_neuron = any(k in cn for k in ["Glut", "Gaba", "GABA", "Dopa", "Sero", "Chol", "Hist"])
    if not is_neuron:
        return "Non_Neuronal"
    # ARH - Arcuate nucleus
    if "ARH" in cn: return "ARH_Neuron"
    # DMH - Dorsomedial nucleus
    if "DMH" in cn: return "DMH_Neuron"
    # PVH/PVpo - Paraventricular nucleus
    if "PVH" in cn or "PVpo" in cn or "PVp" in cn: return "PVH_Neuron"
    # VMH - Ventromedial nucleus
    if "VMH" in cn: return "VMH_Neuron"
    # SCH - Suprachiasmatic nucleus
    if "SCH" in cn: return "SCH_Neuron"
    # LHA - Lateral hypothalamic area
    if "LHA" in cn: return "LHA_Neuron"
    # MPN/MPO - Medial preoptic
    if "MPN" in cn or "MPO" in cn: return "MPO_Neuron"
    # TU - Tuberal nucleus
    if "TU" in cn: return "TU_Neuron"
    # MM - Mammillary body
    if "MM" in cn or "MM-" in cn: return "MM_Neuron"
    # Other neurons
    return "Other_Hypothalamic_Neuron"

# ===========================================================================
# 5. 提取和注释函数
# ===========================================================================
def extract_and_annotate(matched_barcodes, label_name):
    """提取细胞子集并添加完整注释"""
    print(f"\n  [{label_name}] 提取 {len(matched_barcodes):,} cells...")

    # 子集化
    subset = adata[matched_barcodes].to_memory()
    print(f"    子集: {subset.n_obs} cells x {subset.n_vars} genes")

    # 合并cell_metadata
    common = [c for c in subset.obs_names if c in cell_meta.index]
    meta_sub = cell_meta.loc[common]
    print(f"    metadata匹配: {len(common)} cells")

    # 添加基本列 (处理可能缺失的列)
    for col in ['cluster_alias', 'region_of_interest_label',
                'anatomical_division_label', 'donor_label',
                'donor_age', 'donor_sex', 'donor_age_category']:
        if col in meta_sub.columns:
            subset.obs[col] = meta_sub[col]

    # 添加cluster_name (从cluster.csv)
    if 'cluster_alias' in subset.obs.columns:
        cluster_names = dict(zip(cluster_info['cluster_alias'],
                                 cluster_info['cluster_name']))
        subset.obs['cluster_name'] = subset.obs['cluster_alias'].map(cluster_names)

    # 添加cell_type_major
    if 'cluster_name' in subset.obs.columns:
        subset.obs['cell_type_major'] = subset.obs['cluster_name'].apply(
            assign_cell_type_major)
    elif has_crossmap:
        common_cm = [c for c in subset.obs_names if c in cross_map.index]
        subset.obs.loc[common_cm, 'cell_type_major'] = \
            cross_map.loc[common_cm, 'class_name']

    # 添加neuronal_subclass
    if 'cluster_name' in subset.obs.columns:
        subset.obs['neuronal_subclass'] = subset.obs['cluster_name'].apply(
            assign_neuronal_subclass)

    # ---- 添加age_group ----
    subset.obs['age_group'] = 'Unknown'
    if 'donor_label' in subset.obs.columns:
        subset.obs['donor_age_category_donor'] = subset.obs['donor_label'].map(donor_age_map)
        # 按照donor.csv中的age_category分组
        age_cat = subset.obs['donor_age_category_donor']
        subset.obs['age_group'] = age_cat.apply(
            lambda x: 'Aged(18m)' if x == 'aged' else ('Adult(2m)' if x == 'adult' else 'Unknown'))
    elif 'donor_age_category' in subset.obs.columns:
        subset.obs['age_group'] = subset.obs['donor_age_category'].apply(
            lambda x: 'Aged(18m)' if x == 'aged' else ('Adult(2m)' if x == 'adult' else 'Unknown'))

    n_adult = (subset.obs['age_group'] == 'Adult(2m)').sum()
    n_aged = (subset.obs['age_group'] == 'Aged(18m)').sum()
    n_unk = (subset.obs['age_group'] == 'Unknown').sum()
    print(f"    年龄: Adult={n_adult}, Aged={n_aged}, Unknown={n_unk}")

    # ---- 添加sex_group ----
    subset.obs['sex_group'] = 'Unknown'
    if 'donor_sex' in subset.obs.columns:
        sex = subset.obs['donor_sex'].str.strip('"').str.strip()
        subset.obs['sex_group'] = sex.apply(
            lambda x: 'Female' if x == 'F' else ('Male' if x == 'M' else 'Unknown'))

    n_f = (subset.obs['sex_group'] == 'Female').sum()
    n_m = (subset.obs['sex_group'] == 'Male').sum()
    print(f"    性别: Female={n_f}, Male={n_m}")

    # ---- 创建组合分组 ----
    subset.obs['group'] = subset.obs['age_group'] + '_' + subset.obs['sex_group']

    # ---- 统计cell_type_major分布 ----
    print(f"    细胞类型分布:")
    for ct, cnt in subset.obs['cell_type_major'].value_counts().items():
        print(f"      {ct}: {cnt}")

    return subset

# 运行提取
print("\n[4] 提取两个label集合...")
hypo_anat = extract_and_annotate(matched_anat, "anatomical_52696")
hypo_diss = extract_and_annotate(matched_diss, "dissection_26105")

# ===========================================================================
# 6. 保存
# ===========================================================================
print("\n[5] 保存h5ad文件...")

f_anat = os.path.join(OUT_DIR, "hypo_anatomical_52696.h5ad")
hypo_anat.write(f_anat)
size_mb = os.path.getsize(f_anat) / 1024 / 1024
print(f"  {f_anat}")
print(f"    大小: {size_mb:.0f} MB")
print(f"    细胞: {hypo_anat.n_obs:,}")
print(f"    基因: {hypo_anat.n_vars:,}")

f_diss = os.path.join(OUT_DIR, "hypo_dissection_26105.h5ad")
hypo_diss.write(f_diss)
size_mb = os.path.getsize(f_diss) / 1024 / 1024
print(f"  {f_diss}")
print(f"    大小: {size_mb:.0f} MB")
print(f"    细胞: {hypo_diss.n_obs:,}")
print(f"    基因: {hypo_diss.n_vars:,}")

# ===========================================================================
# 7. 保存metadata CSV (方便R读取)
# ===========================================================================
print("\n[6] 保存metadata CSV...")
for df, name in [(hypo_anat, "anatomical"), (hypo_diss, "dissection")]:
    meta_out = df.obs.copy()
    meta_out['cell_label'] = meta_out.index
    csv_path = os.path.join(OUT_DIR, f"hypo_{name}_metadata.csv")
    meta_out.to_csv(csv_path, index=False)
    print(f"  {csv_path} ({len(meta_out)} cells)")

# ===========================================================================
# 8. 摘要
# ===========================================================================
print("\n" + "=" * 60)
print(" Phase 1 提取完成!")
print("=" * 60)
print(f"  anatomical: {hypo_anat.n_obs:,} cells, {len(hypo_anat.obs['cell_type_major'].unique())} cell types")
print(f"  dissection: {hypo_diss.n_obs:,} cells, {len(hypo_diss.obs['cell_type_major'].unique())} cell types")
print(f"\n  Raw counts下载链接:")
print(f"  https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-raw.h5ad")
print(f"  下载后放入: d:/decrepitude mouse hypothamulas/")
