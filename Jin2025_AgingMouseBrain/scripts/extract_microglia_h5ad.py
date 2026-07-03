"""
============================================================================
从全脑h5ad中提取下丘脑Microglia细胞的表达矩阵
保存为小h5ad文件供R分析
============================================================================
用法:
  python scripts/extract_microglia_h5ad.py
  # 或
  python3 scripts/extract_microglia_h5ad.py
============================================================================
"""
import sys
import os

# 先检查anndata是否安装
try:
    import anndata
    import scanpy as sc
    import pandas as pd
    import numpy as np
except ImportError as e:
    print(f"需要安装: {e}")
    print("运行: pip install anndata scanpy pandas numpy")
    sys.exit(1)

# 路径
H5AD_FILE = "d:/decrepitude mouse hypothamulas/Zeng-Aging-Mouse-10Xv3-log2.h5ad"
BARCODE_FILE = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_microglia_barcodes.txt"
OUT_FILE = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_microglia_expression.h5ad"
CELLMETA_FILE = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/cell_metadata.csv"

print("=" * 60)
print(" 提取下丘脑Microglia表达矩阵")
print("=" * 60)
print()

# 读入barcodes
with open(BARCODE_FILE) as f:
    target_barcodes = set(line.strip() for line in f if line.strip())
print(f"目标barcodes: {len(target_barcodes)}")

# 读入表达矩阵 (backed mode, memory efficient)
print(f"\n读取 h5ad (backed模式, 仅加载metadata)...")
adata = sc.read_h5ad(H5AD_FILE, backed='r')
print(f"总细胞: {adata.n_obs:,}")
print(f"总基因: {adata.n_vars:,}")

# 匹配barcodes
all_barcodes = adata.obs_names.tolist()
matched = [b for b in target_barcodes if b in all_barcodes]
print(f"\n匹配到的barcodes: {len(matched)} / {len(target_barcodes)}")

if len(matched) == 0:
    print("\n⚠️ barcode格式不匹配!")
    print(f"  前5个target: {list(target_barcodes)[:5]}")
    print(f"  前5个actual: {all_barcodes[:5]}")

    # Try fuzzy match - strip suffixes
    target_stems = {b.split('-')[0] for b in target_barcodes}
    for b in all_barcodes[:20]:
        stem = b.split('-')[0] if '-' in b else b
        if stem in target_stems:
            print(f"  匹配: {b} (stem={stem})")

    # Try matching all target barcodes more carefully
    matched = [b for b in target_barcodes if b in set(all_barcodes)]
    print(f"\n  精确匹配重试: {len(matched)}")

if len(matched) == 0:
    # Use cluster-based filtering from the h5ad metadata
    print("\n使用cluster筛选 Microglia (840-844) + hypothalamus...")
    if 'cluster_alias' in adata.obs.columns:
        micro_mask = (adata.obs['cluster_alias'].astype(int) >= 840) & \
                     (adata.obs['cluster_alias'].astype(int) <= 844)
        if 'anatomical_division_label' in adata.obs.columns:
            hypo_mask = adata.obs['anatomical_division_label'] == 'HY - HY'
            mask = micro_mask & hypo_mask
        else:
            mask = micro_mask
        matched_cells = adata.obs_names[mask].tolist()
        print(f"  通过cluster+region筛选: {len(matched_cells)} cells")
    else:
        print("  ✗ cluster信息不可用")
        sys.exit(1)
else:
    matched_cells = matched

# 提取子集
print(f"\n提取 {len(matched_cells)} 个Microglia细胞...")
subset_adata = adata[matched_cells].to_memory()

print(f"  子集大小: {subset_adata.n_obs} cells × {subset_adata.n_vars} genes")

# 添加额外的metadata
print("\n添加年龄/性别/聚类注释...")
cell_meta = pd.read_csv(CELLMETA_FILE)
cell_meta.index = cell_meta['cell_label']

# 筛选对应细胞
common = [c for c in subset_adata.obs_names if c in cell_meta.index]
print(f"  匹配metadata: {len(common)} cells")

if len(common) > 0:
    meta_sub = cell_meta.loc[common]
    subset_adata.obs['donor_age_category'] = meta_sub['donor_age_category']
    subset_adata.obs['donor_sex'] = meta_sub['donor_sex']
    subset_adata.obs['cluster_alias'] = meta_sub['cluster_alias']
    subset_adata.obs['region_of_interest_label'] = meta_sub['region_of_interest_label']
    subset_adata.obs['anatomical_division_label'] = meta_sub['anatomical_division_label']

    # 添加年龄分组
    subset_adata.obs['age_group'] = meta_sub['donor_age_category'].apply(
        lambda x: 'Aged(18m)' if x == 'aged' else 'Adult(2m)'
    )

    print(f"\n  年龄分布:")
    print(subset_adata.obs['age_group'].value_counts())

# 保存
print(f"\n保存到: {OUT_FILE}")
subset_adata.write(OUT_FILE)
print(f"  文件大小: {os.path.getsize(OUT_FILE) / 1024 / 1024:.1f} MB")

print("\n✅ 完成! 可在R中加载此文件进行分析")
print(f"   R: zellkonverter::readH5AD('{OUT_FILE}')")
