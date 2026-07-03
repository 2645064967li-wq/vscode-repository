"""
============================================================================
Phase 4: 为CellChat准备数据
从raw counts h5ad提取下丘脑细胞 → 分配cell type → 按年龄分组 → 下采样 → 导出
============================================================================
"""
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import os, gc, warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1

# ===========================================================================
# 配置
# ===========================================================================
RAW_H5AD = "d:/vscode/Jin2025_AgingMouseBrain/Zeng-Aging-Mouse-10Xv3-raw.h5ad"
BARCODE_FILE = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_all_barcodes.txt"
CELLMETA_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/cell_metadata.csv"
CLUSTER_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/cluster.csv"
DONOR_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/donor.csv"
OUT_DIR = "d:/vscode/Jin2025_AgingMouseBrain/data/cellchat"
MAX_CELLS_PER_TYPE = 1500  # 每种cell type每年龄组最多采样数
RANDOM_SEED = 42

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print(" Phase 4: 为CellChat准备数据")
print("=" * 60)

# ===========================================================================
# 1. 加载hypothalamus barcode列表
# ===========================================================================
print("\n[1] 加载下丘脑barcode列表...")
with open(BARCODE_FILE) as f:
    hypo_barcodes = [line.strip() for line in f if line.strip()]
print(f"  下丘脑barcodes: {len(hypo_barcodes):,}")

# ===========================================================================
# 2. 加载metadata
# ===========================================================================
print("\n[2] 加载metadata...")

# 2a. cluster info → cluster_name
cluster_info = pd.read_csv(CLUSTER_CSV)
cluster_to_name = dict(zip(cluster_info['cluster_alias'], cluster_info['cluster_name']))
print(f"  cluster.csv: {len(cluster_info)} clusters")

# 2b. donor info → age & sex
donor_info = pd.read_csv(DONOR_CSV)
donor_age_map = dict(zip(donor_info['donor_label'], donor_info['donor_age_category']))
donor_sex_map = dict(zip(donor_info['donor_label'], donor_info['donor_sex']))
print(f"  donor.csv: {len(donor_info)} donors")
# Count age groups
print(f"    adult donors: {sum(v == 'adult' for v in donor_age_map.values())}")
print(f"    aged donors:  {sum(v == 'aged' for v in donor_age_map.values())}")

# ===========================================================================
# 3. 从cell_metadata中提取下丘脑细胞的注释
# ===========================================================================
print("\n[3] 提取下丘脑细胞metadata...")

# cell_metadata很大 (400MB)，分块读取只提取需要的列
needed_cols = ['cell_label', 'cluster_alias', 'donor_label', 'donor_age_category']
print(f"  读取cell_metadata (筛选{len(hypo_barcodes):,}个下丘脑细胞)...")

# 创建下丘脑barcode集合用于快速查找
hypo_set = set(hypo_barcodes)

# 分块读取csv
chunks = []
for chunk in pd.read_csv(CELLMETA_CSV, usecols=needed_cols, chunksize=200000):
    mask = chunk['cell_label'].isin(hypo_set)
    chunks.append(chunk[mask])
    if len(chunks) % 5 == 0:
        total = sum(len(c) for c in chunks)
        print(f"    已处理... 匹配 {total:,} cells")

hypo_meta = pd.concat(chunks, ignore_index=True)
print(f"  匹配的hypothalamus细胞: {len(hypo_meta):,}")

# 确保barcodes一致
hypo_barcodes_matched = hypo_meta['cell_label'].tolist()
print(f"  最终barcode列表: {len(hypo_barcodes_matched):,}")

# ===========================================================================
# 4. 分配cell_type_major
# ===========================================================================
print("\n[4] 分配cell type...")

def assign_cell_type_major(cluster_name):
    """根据cluster_name分配主要细胞类型"""
    if not isinstance(cluster_name, str):
        return "Unclassified"
    cn = cluster_name

    # === 胶质细胞 ===
    if "Microglia" in cn:
        return "Microglia"
    if "Astro" in cn and "Ependymal" not in cn.lower():
        return "Astrocyte"
    if "Tanycyte" in cn:
        return "Tanycyte"
    if "Ependymal" in cn or "ependymal" in cn:
        return "Ependymal"

    # === 少突胶质细胞系 ===
    if "_COP" in cn or cn.startswith("COP") or "COP_" in cn:
        return "COP"
    if "_NFOL" in cn or cn.startswith("NFOL") or "NFOL_" in cn:
        return "NFOL"
    if "_MFOL" in cn or cn.startswith("MFOL") or "MFOL_" in cn:
        return "MFOL"
    if "_MOL" in cn or cn.startswith("MOL") or "MOL_" in cn:
        return "MOL"
    if "_OPC" in cn or cn.startswith("OPC") or "OPC_" in cn:
        return "OPC"
    if "Oligo" in cn:
        return "Oligo_Other"

    # === 血管细胞 ===
    if "Endo" in cn and ("Endo-" in cn or "Endo_" in cn):
        return "Endothelial"
    if "SMC" in cn or "SMC-" in cn:
        return "SMC"
    if "VLMC" in cn:
        return "VLMC"
    if "Peri" in cn:
        return "Pericyte"

    # === 免疫细胞 ===
    if "BAM" in cn:
        return "BAM"
    if "_DC" in cn or "DC_" in cn:
        return "Dendritic_Cell"
    if "ABC" in cn:
        return "ABC"
    if "T cells" in cn or "Tcell" in cn or "T_cell" in cn:
        return "T_cell"

    # === 神经元 ===
    if "Glut" in cn and "Gaba" not in cn and "GABA" not in cn:
        return "Glutamatergic_Neuron"
    if "Gaba" in cn or "GABA" in cn:
        return "GABAergic_Neuron"
    if "Dopa" in cn:
        return "Dopaminergic_Neuron"
    if "Sero" in cn:
        return "Serotonergic_Neuron"
    if "Chol" in cn:
        return "Cholinergic_Neuron"
    if "Hist" in cn:
        return "Histaminergic_Neuron"

    # === 其他神经元 ===
    if "NN" in cn and any(k in cn for k in ["Glut", "Gaba", "Dopa", "Sero", "Chol"]):
        return "Neuron_Other"

    return "Unclassified"

# 添加cluster_name
hypo_meta['cluster_name'] = hypo_meta['cluster_alias'].map(cluster_to_name)
n_missing = hypo_meta['cluster_name'].isna().sum()
print(f"  cluster_name缺失: {n_missing} / {len(hypo_meta)}")

# 分配cell_type_major
hypo_meta['cell_type_major'] = hypo_meta['cluster_name'].apply(assign_cell_type_major)

# ===========================================================================
# 5. 分配age_group
# ===========================================================================
print("\n[5] 分配年龄分组...")

# 优先使用cell_metadata的donor_age_category，缺失的用donor.csv补充
hypo_meta['age_group'] = hypo_meta['donor_age_category'].apply(
    lambda x: 'Aged(18m)' if x == 'aged' else ('Adult(2m)' if x == 'adult' else 'Unknown')
)

# 用donor.csv补充
missing_age = hypo_meta['age_group'] == 'Unknown'
if missing_age.sum() > 0:
    print(f"  补充{donor_age_map}中缺失年龄的 {missing_age.sum()} cells")
    hypo_meta.loc[missing_age, 'age_group'] = (
        hypo_meta.loc[missing_age, 'donor_label']
        .map(donor_age_map)
        .apply(lambda x: 'Aged(18m)' if x == 'aged' else ('Adult(2m)' if x == 'adult' else 'Unknown'))
    )

print("  年龄分布:")
for ag, cnt in hypo_meta['age_group'].value_counts().items():
    print(f"    {ag}: {cnt:,}")

print("\n  Cell type分布 (by age):")
ct_summary = hypo_meta.groupby(['cell_type_major', 'age_group']).size().unstack(fill_value=0)
print(ct_summary.to_string())

# ===========================================================================
# 6. 从raw h5ad提取counts + 下采样
# ===========================================================================
print(f"\n[6] 从raw h5ad提取raw counts并下采样 (max {MAX_CELLS_PER_TYPE} cells/type/age)...")

# 加载raw h5ad (backed mode)
print("  加载raw h5ad (backed mode)...")
raw_adata = anndata.read_h5ad(RAW_H5AD, backed='r')
print(f"  总细胞: {raw_adata.n_obs:,}, 总基因: {raw_adata.n_vars:,}")

# 按年龄组分别处理
np.random.seed(RANDOM_SEED)

for age_group in ['Adult(2m)', 'Aged(18m)']:
    print(f"\n  --- {age_group} ---")

    # 筛选该年龄组的细胞
    age_cells = hypo_meta[hypo_meta['age_group'] == age_group]
    print(f"    该年龄组总细胞: {len(age_cells):,}")

    # 按cell type下采样
    sampled_barcodes = []
    sampled_meta_list = []

    for ct in sorted(age_cells['cell_type_major'].unique()):
        ct_cells = age_cells[age_cells['cell_type_major'] == ct]
        n_total = len(ct_cells)
        n_sample = min(n_total, MAX_CELLS_PER_TYPE)

        if n_sample < 10:  # 太少的cell type跳过
            continue

        if n_total <= n_sample:
            sampled = ct_cells
        else:
            sampled = ct_cells.sample(n=n_sample, random_state=RANDOM_SEED)

        sampled_barcodes.extend(sampled['cell_label'].tolist())
        sampled_meta_list.append(sampled)
        print(f"      {ct}: {n_total:,} → {n_sample:,}")

    # 从raw h5ad提取这些细胞的counts
    valid_bc = [bc for bc in sampled_barcodes if bc in raw_adata.obs_names]
    print(f"    有效barcodes: {len(valid_bc):,} / {len(sampled_barcodes):,}")

    # 提取子集到内存
    print(f"    提取counts矩阵...")
    subset = raw_adata[valid_bc].to_memory()
    print(f"    提取完成: {subset.n_obs} cells x {subset.n_vars} genes")

    # 获取counts矩阵 (raw h5ad的X应该是raw counts)
    if sparse.issparse(subset.X):
        count_matrix = subset.X.toarray()
    else:
        count_matrix = subset.X

    # 构建cell metadata
    sampled_meta_combined = pd.concat(sampled_meta_list, ignore_index=True)
    sampled_meta_combined = sampled_meta_combined[sampled_meta_combined['cell_label'].isin(valid_bc)]
    sampled_meta_combined = sampled_meta_combined.set_index('cell_label').loc[valid_bc].reset_index()

    # 只保留需要的列
    meta_out = sampled_meta_combined[['cell_label', 'cell_type_major', 'age_group',
                                        'cluster_alias', 'cluster_name']].copy()
    meta_out.columns = ['cell_label', 'cell_type', 'age_group', 'cluster_alias', 'cluster_name']

    # =========================================================================
    # 7. 保存
    # =========================================================================
    age_label = 'young' if age_group == 'Adult(2m)' else 'aged'
    print(f"\n    保存 {age_label} 数据...")

    # 保存为压缩CSV (counts: genes × cells)
    count_df = pd.DataFrame(
        count_matrix.T,  # genes × cells
        index=subset.var_names,
        columns=valid_bc
    )
    count_path = os.path.join(OUT_DIR, f"{age_label}_raw_counts.csv.gz")
    count_df.to_csv(count_path, compression='gzip')
    print(f"    Counts: {count_path}")
    print(f"      shape: {count_df.shape}")

    # 保存metadata
    meta_path = os.path.join(OUT_DIR, f"{age_label}_metadata.csv")
    meta_out.to_csv(meta_path, index=False)
    print(f"    Metadata: {meta_path}")
    print(f"      shape: {meta_out.shape}")

    # 清理内存
    del subset, count_matrix, count_df
    gc.collect()

# ===========================================================================
# 8. 汇总
# ===========================================================================
print("\n" + "=" * 60)
print(" Phase 4 数据准备完成!")
print("=" * 60)

young_meta = pd.read_csv(os.path.join(OUT_DIR, "young_metadata.csv"))
aged_meta = pd.read_csv(os.path.join(OUT_DIR, "aged_metadata.csv"))

print(f"\n  Young(2m): {len(young_meta):,} cells, "
      f"{young_meta['cell_type'].nunique()} cell types")
print(f"  Aged(18m):  {len(aged_meta):,} cells, "
      f"{aged_meta['cell_type'].nunique()} cell types")

# 确认共同的cell types
common_ct = set(young_meta['cell_type'].unique()) & set(aged_meta['cell_type'].unique())
print(f"\n  共同的cell types: {len(common_ct)}")
for ct in sorted(common_ct):
    y_cnt = (young_meta['cell_type'] == ct).sum()
    a_cnt = (aged_meta['cell_type'] == ct).sum()
    print(f"    {ct:30s}  Young={y_cnt:,}  Aged={a_cnt:,}")

# 验证Microglia存在
if 'Microglia' in common_ct:
    print(f"\n  ✓ Microglia存在于两组数据中")
else:
    print(f"\n  ⚠️ WARNING: Microglia未在两组中同时出现!")

print(f"\n  输出目录: {OUT_DIR}/")
print(f"  文件:")
for f in sorted(os.listdir(OUT_DIR)):
    size_mb = os.path.getsize(os.path.join(OUT_DIR, f)) / 1024 / 1024
    print(f"    {f} ({size_mb:.1f} MB)")
