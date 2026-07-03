"""
============================================================================
Phase 4 v5: Raw counts extraction - scanpy backed + batch extraction
============================================================================
"""
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import os, gc, warnings
warnings.filterwarnings('ignore')

RAW_H5AD = "d:/vscode/Jin2025_AgingMouseBrain/Zeng-Aging-Mouse-10Xv3-raw.h5ad"
BARCODE_FILE = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_all_barcodes.txt"
CELLMETA_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/cell_metadata.csv"
CLUSTER_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/cluster.csv"
DONOR_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/donor.csv"
OUT_DIR = "d:/vscode/Jin2025_AgingMouseBrain/data/cellchat"
MAX_CELLS_PER_TYPE = 1500
BATCH_SIZE = 500  # Extract 500 cells at a time
RANDOM_SEED = 42

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print(" Phase 4 v5: Batch extraction from raw h5ad")
print("=" * 60)

# ===========================================================================
# 1-3. Load metadata, assign cell types, sample (same as v4)
# ===========================================================================
print("\n[1] Loading barcodes & metadata...")
with open(BARCODE_FILE) as f:
    hypo_set = set(line.strip() for line in f if line.strip())
print(f"  {len(hypo_set):,} hypothalamic barcodes")

cluster_info = pd.read_csv(CLUSTER_CSV)
cluster_to_name = dict(zip(cluster_info['cluster_alias'], cluster_info['cluster_name']))
donor_info = pd.read_csv(DONOR_CSV)
donor_age_map = dict(zip(donor_info['donor_label'], donor_info['donor_age_category']))

print("  Loading cell_metadata (chunked)...")
needed_cols = ['cell_label', 'cluster_alias', 'donor_label', 'donor_age_category']
chunks = []
for chunk in pd.read_csv(CELLMETA_CSV, usecols=needed_cols, chunksize=200000):
    m = chunk['cell_label'].isin(hypo_set)
    if m.any(): chunks.append(chunk[m])
hypo_meta = pd.concat(chunks, ignore_index=True)
print(f"  Matched: {len(hypo_meta):,}")

print("\n[2] Assigning cell types...")
def assign_ct(cn):
    if not isinstance(cn, str): return "Unclassified"
    cn_s = cn
    if "Microglia" in cn_s: return "Microglia"
    if "Astro" in cn_s and "Ependymal" not in cn_s.lower(): return "Astrocyte"
    if "Tanycyte" in cn_s: return "Tanycyte"
    if "Ependymal" in cn_s or "ependymal" in cn_s: return "Ependymal"
    if "COP" in cn_s: return "COP"
    if "NFOL" in cn_s: return "NFOL"
    if "MFOL" in cn_s: return "MFOL"
    if "MOL" in cn_s: return "MOL"
    if "OPC" in cn_s: return "OPC"
    if "Endo" in cn_s: return "Endothelial"
    if "SMC" in cn_s: return "SMC"
    if "VLMC" in cn_s: return "VLMC"
    if "Peri" in cn_s: return "Pericyte"
    if "BAM" in cn_s: return "BAM"
    if "T cells" in cn_s or "Tcell" in cn_s: return "T_cell"
    if "Glut" in cn_s and "Gaba" not in cn_s and "GABA" not in cn_s: return "Glutamatergic_Neuron"
    if "Gaba" in cn_s or "GABA" in cn_s: return "GABAergic_Neuron"
    return "Unclassified"

hypo_meta['cluster_name'] = hypo_meta['cluster_alias'].map(cluster_to_name)
hypo_meta['cell_type_major'] = hypo_meta['cluster_name'].apply(assign_ct)
hypo_meta['age_group'] = hypo_meta['donor_age_category'].apply(
    lambda x: 'Aged(18m)' if x == 'aged' else ('Adult(2m)' if x == 'adult' else 'Unknown')
)
m = hypo_meta['age_group'] == 'Unknown'
if m.sum() > 0:
    hypo_meta.loc[m, 'age_group'] = (
        hypo_meta.loc[m, 'donor_label'].map(donor_age_map)
        .apply(lambda x: 'Aged(18m)' if x == 'aged' else ('Adult(2m)' if x == 'adult' else 'Unknown'))
    )

# ===========================================================================
# 4. Sample cells
# ===========================================================================
print(f"\n[3] Sampling (max {MAX_CELLS_PER_TYPE} per type per age)...")
np.random.seed(RANDOM_SEED)

# Build plan: {age_group: {cell_type: [barcodes]}}
plan = {}
for age in ['Adult(2m)', 'Aged(18m)']:
    plan[age] = {}
    age_cells = hypo_meta[hypo_meta['age_group'] == age]
    for ct in sorted(age_cells['cell_type_major'].unique()):
        ct_cells = age_cells[age_cells['cell_type_major'] == ct]
        n = min(len(ct_cells), MAX_CELLS_PER_TYPE)
        if n < 10: continue
        sampled = ct_cells if len(ct_cells) <= n else ct_cells.sample(n=n, random_state=RANDOM_SEED)
        plan[age][ct] = sampled['cell_label'].tolist()
        print(f"    {age} {ct}: {len(ct_cells):,} -> {n:,}")

# ===========================================================================
# 5. Extract from raw h5ad in batches
# ===========================================================================
print("\n[4] Extracting raw counts (batch mode)...")
raw = anndata.read_h5ad(RAW_H5AD, backed='r')
print(f"  Raw h5ad: {raw.n_obs:,} cells x {raw.n_vars:,} genes")

# Gene symbols
gene_symbols = raw.var['gene_symbol'].astype(str).values
# Handle duplicates
dup_mask = pd.Series(gene_symbols).duplicated(keep=False)
if dup_mask.any():
    print(f"  Fixing {dup_mask.sum()} duplicate gene symbols...")
    gene_symbols[dup_mask] = [f"{g}_{i}" for i, g in enumerate(gene_symbols[dup_mask])]
print(f"  Genes: {len(gene_symbols)}")

for age_group in ['Adult(2m)', 'Aged(18m)']:
    label = 'young' if age_group == 'Adult(2m)' else 'aged'
    print(f"\n  --- {age_group} ({label}) ---")

    all_rows = []
    all_bcs = []
    all_cts = []
    n_done = 0

    for ct, bcs in plan[age_group].items():
        # Extract in batches
        for i in range(0, len(bcs), BATCH_SIZE):
            batch_bc = bcs[i:i+BATCH_SIZE]
            # Filter to barcodes that exist in raw h5ad
            valid = [bc for bc in batch_bc if bc in raw.obs_names]
            if not valid: continue

            # Extract subset
            subset = raw[valid].to_memory()
            if sparse.issparse(subset.X):
                batch_data = subset.X.toarray()
            else:
                batch_data = subset.X

            all_rows.append(batch_data.astype(np.float32))
            all_bcs.extend(valid)
            all_cts.extend([ct] * len(valid))

            n_done += len(valid)
            pct = n_done * 100 // sum(len(v) for v in plan[age_group].values())
            print(f"      Extracted {n_done:,} cells ({pct}%)...")
            gc.collect()

    # Combine all batches
    print(f"    Combining {len(all_rows)} batches...")
    full_mat = np.vstack(all_rows)
    del all_rows
    gc.collect()

    print(f"    Final matrix: {full_mat.shape}")

    # Save as CSV (genes x cells)
    count_df = pd.DataFrame(full_mat.T, index=gene_symbols, columns=all_bcs)
    count_df.index.name = 'gene'

    # Handle any remaining duplicate index
    if count_df.index.duplicated().any():
        dup_count = count_df.index.duplicated().sum()
        print(f"    WARNING: {dup_count} duplicate gene symbols in output, making unique...")
        idx = count_df.index.tolist()
        seen = {}
        new_idx = []
        for g in idx:
            if g in seen:
                seen[g] += 1
                new_idx.append(f"{g}_{seen[g]}")
            else:
                seen[g] = 0
                new_idx.append(g)
        count_df.index = new_idx

    count_path = os.path.join(OUT_DIR, f"{label}_raw_counts.csv.gz")
    count_df.to_csv(count_path, compression='gzip')
    sz = os.path.getsize(count_path) / 1024 / 1024
    print(f"    Saved: {count_path} ({sz:.1f} MB)")

    # Metadata
    meta = pd.DataFrame({'cell_label': all_bcs, 'cell_type': all_cts, 'age_group': age_group})
    meta_path = os.path.join(OUT_DIR, f"{label}_metadata.csv")
    meta.to_csv(meta_path, index=False)
    print(f"    Metadata: {meta_path} ({len(meta)} cells)")

    del full_mat, count_df
    gc.collect()

raw.file.close()  # Close backed file

# ===========================================================================
# 6. Summary
# ===========================================================================
print("\n" + "=" * 60)
print(" Phase 4 v5 complete!")
print("=" * 60)
for label in ['young', 'aged']:
    meta = pd.read_csv(os.path.join(OUT_DIR, f"{label}_metadata.csv"))
    print(f"\n{label}: {len(meta):,} cells, {meta['cell_type'].nunique()} types")
    for ct, cnt in meta['cell_type'].value_counts().items():
        print(f"  {ct}: {cnt}")
ym = pd.read_csv(os.path.join(OUT_DIR, "young_metadata.csv"))
am = pd.read_csv(os.path.join(OUT_DIR, "aged_metadata.csv"))
print(f"\nMicroglia: Young={('Microglia' in ym['cell_type'].values)}, Aged={('Microglia' in am['cell_type'].values)}")
