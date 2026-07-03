"""
============================================================================
Phase 4 v2: 为CellChat准备数据 (使用现有log2 hypothalamus h5ad)
更轻量 - 直接从已注释的hypothalamus h5ad导出，避免读取13GB raw h5ad
============================================================================
"""
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
HYPO_H5AD = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_anatomical_52696.h5ad"
OUT_DIR = "d:/vscode/Jin2025_AgingMouseBrain/data/cellchat"
MAX_CELLS_PER_TYPE = 1500
RANDOM_SEED = 42

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print(" Phase 4 v2: 从log2 hypothalamus h5ad准备CellChat数据")
print("=" * 60)

# ===========================================================================
# 1. 加载hypothalamus h5ad
# ===========================================================================
print("\n[1] 加载hypothalamus h5ad...")
adata = sc.read_h5ad(HYPO_H5AD)
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print(f"  obs columns: {list(adata.obs.columns)}")

# ===========================================================================
# 2. 检查cell_type_major和age_group
# ===========================================================================
print("\n[2] 检查cell type和age分布...")

if 'cell_type_major' not in adata.obs.columns:
    print("  ERROR: cell_type_major not found in adata.obs!")
    exit(1)

print("  Cell types:")
for ct, cnt in adata.obs['cell_type_major'].value_counts().items():
    print(f"    {ct}: {cnt:,}")

print("\n  Age distribution:")
for ag, cnt in adata.obs['age_group'].value_counts().items():
    print(f"    {ag}: {cnt:,}")

# ===========================================================================
# 3. 转换log2数据为raw-like counts
# ===========================================================================
print("\n[3] 数据转换...")

# h5ad中存储的是log2(CPM+1) 或 log2(raw+1)
# CellChat需要raw counts，但log数据也可以运行
# 将log2转回近似raw: raw ≈ 2^log2 - 1 (会丢失精度但对CellChat足够)
if sparse.issparse(adata.X):
    data_matrix = adata.X.toarray()
else:
    data_matrix = adata.X.copy()
print(f"  Data matrix shape: {data_matrix.shape}")

# 检查数据范围判断是否为log2
max_val = data_matrix.max()
min_val = data_matrix.min()
print(f"  Value range: [{min_val:.2f}, {max_val:.2f}]")

# 将log2转回线性空间
# Input: log2(x + 1), Output: x ≈ 2^vals - 1
print("  Converting log2 → linear (approx raw counts)...")
data_linear = np.exp2(data_matrix) - 1
data_linear = np.maximum(data_linear, 0)  # ensure non-negative
data_linear = np.round(data_linear).astype(np.int32)  # round to integers
print(f"  Linear value range: [{data_linear.min()}, {data_linear.max()}]")

# ===========================================================================
# 4. 按age分组并下采样
# ===========================================================================
print(f"\n[4] 下采样 (max {MAX_CELLS_PER_TYPE} cells/type/age)...")
np.random.seed(RANDOM_SEED)

for age_group in ['Adult(2m)', 'Aged(18m)']:
    print(f"\n  --- {age_group} ---")

    # 筛选该年龄组
    age_mask = adata.obs['age_group'] == age_group
    age_adata = adata[age_mask].copy()
    print(f"    该年龄组总细胞: {age_adata.n_obs:,}")

    # 按cell type下采样
    sampled_barcodes = []
    sampled_meta_list = []

    for ct in sorted(age_adata.obs['cell_type_major'].unique()):
        ct_mask = age_adata.obs['cell_type_major'] == ct
        ct_cells = age_adata[ct_mask]
        n_total = ct_cells.n_obs
        n_sample = min(n_total, MAX_CELLS_PER_TYPE)

        if n_sample < 10:
            continue

        ct_indices = np.where(ct_mask)[0]
        if n_total <= n_sample:
            sampled_idx = ct_indices
        else:
            sampled_idx = np.random.choice(ct_indices, size=n_sample, replace=False)

        sampled_barcodes.extend(age_adata.obs_names[sampled_idx].tolist())

        meta_subset = age_adata.obs.iloc[sampled_idx][['cell_type_major', 'age_group']].copy()
        meta_subset['cell_label'] = meta_subset.index
        sampled_meta_list.append(meta_subset)

        print(f"      {ct}: {n_total:,} → {n_sample:,}")

    # 提取counts
    valid_bc = sampled_barcodes
    bc_to_idx = {bc: i for i, bc in enumerate(adata.obs_names)}
    valid_idx = [bc_to_idx[bc] for bc in valid_bc if bc in bc_to_idx]
    valid_bc = [bc for bc in valid_bc if bc in bc_to_idx]
    print(f"    最终细胞数: {len(valid_bc):,}")

    # 提取子矩阵
    subset_data = data_linear[valid_idx, :]

    # 构建counts DataFrame
    count_df = pd.DataFrame(
        subset_data.T,
        index=adata.var_names,
        columns=valid_bc
    )
    count_df.index.name = 'gene'

    age_label = 'young' if age_group == 'Adult(2m)' else 'aged'

    # 保存counts
    count_path = os.path.join(OUT_DIR, f"{age_label}_raw_counts.csv.gz")
    count_df.to_csv(count_path, compression='gzip')
    print(f"    Counts saved: {count_path} ({count_df.shape[0]} genes x {count_df.shape[1]} cells)")

    # 保存metadata
    meta_out = pd.concat(sampled_meta_list, ignore_index=True)
    meta_out = meta_out[meta_out['cell_label'].isin(valid_bc)]
    meta_out.columns = ['cell_type', 'age_group', 'cell_label']
    meta_out = meta_out[['cell_label', 'cell_type', 'age_group']]

    meta_path = os.path.join(OUT_DIR, f"{age_label}_metadata.csv")
    meta_out.to_csv(meta_path, index=False)
    print(f"    Metadata saved: {meta_path} ({len(meta_out)} cells)")

    del count_df, subset_data
    gc.collect()

# ===========================================================================
# 5. 汇总
# ===========================================================================
print("\n" + "=" * 60)
print(" Phase 4 v2 完成!")
print("=" * 60)

for label in ['young', 'aged']:
    meta = pd.read_csv(os.path.join(OUT_DIR, f"{label}_metadata.csv"))
    print(f"\n  {label}: {len(meta):,} cells, {meta['cell_type'].nunique()} cell types")
    for ct, cnt in meta['cell_type'].value_counts().items():
        print(f"    {ct}: {cnt}")

# 检查Microglia
young_meta = pd.read_csv(os.path.join(OUT_DIR, "young_metadata.csv"))
aged_meta = pd.read_csv(os.path.join(OUT_DIR, "aged_metadata.csv"))
has_mg = 'Microglia' in young_meta['cell_type'].values and 'Microglia' in aged_meta['cell_type'].values
print(f"\n  Microglia available: {has_mg}")

print(f"\n  Files in {OUT_DIR}:")
for f in sorted(os.listdir(OUT_DIR)):
    sz = os.path.getsize(os.path.join(OUT_DIR, f)) / 1024 / 1024
    print(f"    {f} ({sz:.1f} MB)")
