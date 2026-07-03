"""
============================================================================
Phase 2-3: Python scanpy 下丘脑全分析
QC + 降维 + 分性别per-cell-type DEG
============================================================================
"""
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
import os, warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1

RESULT_DIR = "d:/vscode/Jin2025_AgingMouseBrain/results"
os.makedirs(f"{RESULT_DIR}/DEGs", exist_ok=True)
os.makedirs(f"{RESULT_DIR}/figures", exist_ok=True)

# ===========================================================================
# 1. 加载数据
# ===========================================================================
print("=" * 60)
print(" Phase 2-3: 下丘脑衰老分析")
print("=" * 60)

H5AD = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_anatomical_52696.h5ad"
adata = sc.read_h5ad(H5AD)
print(f"\nLoaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# ===========================================================================
# 2. 汇总统计
# ===========================================================================
print("\n" + "=" * 60)
print(" 细胞类型 × 年龄 × 性别 汇总")
print("=" * 60)

summary = adata.obs.groupby(['cell_type_major','age_group','sex_group']).size().unstack(fill_value=0)
print(summary.to_string())

summary.to_csv(f"{RESULT_DIR}/hypo_cell_type_summary.csv")

# ===========================================================================
# 3. QC + 降维 (对所有细胞)
# ===========================================================================
print("\n[Step 1] 降维...")
sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat')
print(f"  HVGs: {adata.var.highly_variable.sum()}")

# PCA
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)
# UMAP
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)
sc.tl.umap(adata, min_dist=0.3)

# 保存UMAP坐标
umap_coords = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP1','UMAP2'])
umap_coords.to_csv(f"{RESULT_DIR}/umap_coordinates.csv")

print("  [OK] PCA + UMAP done")

# ===========================================================================
# 4. Per-cell-type DEG: Aged vs Adult (分性别)
# ===========================================================================
print("\n[Step 2] Per-cell-type DEG分析...")

DEG_MIN_CELLS = 10
cell_types = [ct for ct in adata.obs['cell_type_major'].unique()
              if (adata.obs['cell_type_major'] == ct).sum() >= DEG_MIN_CELLS * 2]

all_degs = {}
de_summary = []

for ct in cell_types:
    ct_mask = adata.obs['cell_type_major'] == ct
    ct_adata = adata[ct_mask].copy()

    print(f"\n  --- {ct}: {ct_adata.n_obs} cells ---")

    for sex in ['Female', 'Male']:
        sex_mask = ct_adata.obs['sex_group'] == sex
        sub = ct_adata[sex_mask]

        n_aged = (sub.obs['age_group'] == 'Aged(18m)').sum()
        n_adult = (sub.obs['age_group'] == 'Adult(2m)').sum()

        if n_aged < DEG_MIN_CELLS or n_adult < DEG_MIN_CELLS:
            print(f"    {sex}: Adult={n_adult}, Aged={n_aged} → 跳过")
            continue

        print(f"    {sex}: Adult={n_adult}, Aged={n_aged}")

        # 设置分组
        sub.obs['comparison'] = sub.obs['age_group'].astype(str)
        sc.tl.rank_genes_groups(
            sub, groupby='comparison',
            groups=['Aged(18m)'], reference='Adult(2m)',
            method='wilcoxon', n_genes=sub.n_vars,
            key_added='aged_vs_adult'
        )

        # 提取结果
        result = sc.get.rank_genes_groups_df(sub, group='Aged(18m)', key='aged_vs_adult')
        result = result.rename(columns={
            'names': 'gene', 'scores': 'z_score',
            'logfoldchanges': 'avg_log2FC', 'pvals': 'p_val',
            'pvals_adj': 'p_val_adj'
        })

        n_sig = (result['p_val_adj'] < 0.05).sum()
        n_up = ((result['p_val_adj'] < 0.05) & (result['avg_log2FC'] > 0)).sum()
        n_down = ((result['p_val_adj'] < 0.05) & (result['avg_log2FC'] < 0)).sum()
        print(f"      DEGs: {n_sig} sig (UP {n_up} DOWN {n_down})")

        # 保存
        key = f"{ct}_{sex}"
        result.to_csv(f"{RESULT_DIR}/DEGs/{key}_DEG.csv", index=False)
        all_degs[key] = result

        de_summary.append({
            'cell_type': ct, 'sex': sex,
            'n_cells': sub.n_obs, 'n_adult': n_adult, 'n_aged': n_aged,
            'n_sig': n_sig, 'n_up': n_up, 'n_down': n_down
        })

# ===========================================================================
# 5. 关键神经元亚类DEG
# ===========================================================================
print("\n[Step 3] 神经元亚类DEG...")

if 'neuronal_subclass' in adata.obs.columns:
    key_subclasses = ['ARH_Neuron', 'DMH_Neuron', 'PVH_Neuron',
                      'VMH_Neuron', 'SCH_Neuron', 'LHA_Neuron']

    for sc_name in key_subclasses:
        if sc_name not in adata.obs['neuronal_subclass'].values:
            continue

        sc_mask = adata.obs['neuronal_subclass'] == sc_name
        sc_adata = adata[sc_mask]
        print(f"\n  {sc_name}: {sc_adata.n_obs} cells")

        if sc_adata.n_obs < 30:
            print(f"    → 跳过 (too few cells)")
            continue

        for sex in ['Female', 'Male']:
            sex_mask = sc_adata.obs['sex_group'] == sex
            sub = sc_adata[sex_mask]
            n_aged = (sub.obs['age_group'] == 'Aged(18m)').sum()
            n_adult = (sub.obs['age_group'] == 'Adult(2m)').sum()

            if n_aged < 5 or n_adult < 5:
                continue

            sc.tl.rank_genes_groups(
                sub, groupby='age_group',
                groups=['Aged(18m)'], reference='Adult(2m)',
                method='wilcoxon', n_genes=sub.n_vars,
                key_added='aged_vs_adult'
            )
            result = sc.get.rank_genes_groups_df(sub, group='Aged(18m)', key='aged_vs_adult')
            result = result.rename(columns={
                'names': 'gene', 'scores': 'z_score',
                'logfoldchanges': 'avg_log2FC', 'pvals': 'p_val',
                'pvals_adj': 'p_val_adj'
            })
            n_sig = (result['p_val_adj'] < 0.05).sum()
            print(f"    {sex}: {n_sig} DEGs (Adult={n_adult}, Aged={n_aged})")
            result.to_csv(f"{RESULT_DIR}/DEGs/{sc_name}_{sex}_DEG.csv", index=False)
            all_degs[f"{sc_name}_{sex}"] = result

# ===========================================================================
# 6. DEG汇总表
# ===========================================================================
de_summary_df = pd.DataFrame(de_summary)
de_summary_df = de_summary_df.sort_values('n_sig', ascending=False)
de_summary_df.to_csv(f"{RESULT_DIR}/DEGs/DEG_summary_all.csv", index=False)

print(f"\n[Step 4] DEG汇总:")

# Top aging-sensitive cell types
print("\n  衰老最敏感的细胞类型 (按DEG数量):")
print(de_summary_df[['cell_type','sex','n_cells','n_sig','n_up','n_down']].head(20).to_string())

# ===========================================================================
# 7. 关键基因检查 (AD + 炎症 + 小胶质激活)
# ===========================================================================
print("\n[Step 5] 关键基因表达检查...")

ad_genes = ['Apoe', 'Trem2', 'Tyrobp', 'Cd33', 'Clu', 'Cd68', 'Cd74',
            'Aif1', 'Cx3cr1', 'Tmem119', 'P2ry12', 'Cst7', 'Lpl', 'Spp1',
            'Il1b', 'Tnf', 'Ccl2', 'C1qa', 'C1qb', 'C1qc',
            'App', 'Bace1', 'Psen1', 'Adam10']

# 在Microglia中检查
if 'Microglia' in adata.obs['cell_type_major'].values:
    mg = adata[adata.obs['cell_type_major'] == 'Microglia']
    mg_aged = mg[mg.obs['age_group'] == 'Aged(18m)']
    mg_adult = mg[mg.obs['age_group'] == 'Adult(2m)']

    gene_check = []
    for gene in ad_genes:
        # Find gene index (try exact match, then case-insensitive)
        g_idx = None
        if gene in mg.var_names:
            g_idx = gene
        elif gene.capitalize() in mg.var_names:
            g_idx = gene.capitalize()
        elif gene.upper() in mg.var_names:
            g_idx = gene.upper()

        if g_idx:
            aged_expr = mg_aged[:, g_idx].X.toarray().mean()
            adult_expr = mg_adult[:, g_idx].X.toarray().mean()
            fc = aged_expr / (adult_expr + 1e-6)
            # Wilcoxon test
            aged_vals = mg_aged[:, g_idx].X.toarray().flatten()
            adult_vals = mg_adult[:, g_idx].X.toarray().flatten()
            try:
                _, pval = stats.mannwhitneyu(aged_vals, adult_vals)
            except:
                pval = 1.0
            gene_check.append({
                'gene': g_idx, 'aged_mean': aged_expr,
                'adult_mean': adult_expr, 'log2FC': np.log2(fc),
                'p_value': pval
            })

    gene_df = pd.DataFrame(gene_check)
    gene_df = gene_df.sort_values('log2FC', ascending=False)
    gene_df.to_csv(f"{RESULT_DIR}/microglia_AD_genes_check.csv", index=False)

    print("\n  Microglia中AD/炎症基因表达变化:")
    for _, row in gene_df.iterrows():
        sig = '*' if row['p_value'] < 0.05 else ' '
        print(f"    {row['gene']:15s} | log2FC={row['log2FC']:+6.2f} | p={row['p_value']:.2e} {sig}")

# ===========================================================================
# 8. 保存
# ===========================================================================
adata.write("d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_anatomical_analyzed.h5ad")
print(f"\n[OK] 已保存: hypo_anatomical_analyzed.h5ad")
print(f"\n{'='*60}")
print(f" Phase 2-3 完成!")
print(f" 总DEG分析: {len(all_degs)} 组")
print(f" 结果目录: {RESULT_DIR}/DEGs/")
print(f"{'='*60}")
