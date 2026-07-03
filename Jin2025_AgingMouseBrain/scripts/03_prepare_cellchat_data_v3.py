"""
============================================================================
Phase 4 v3: 为CellChat准备数据 (使用现有log2 hypothalamus h5ad)
内存友好版 - 使用sparse matrix，分cell type处理
============================================================================
"""
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import os, gc, warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1

HYPO_H5AD = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_anatomical_52696.h5ad"
OUT_DIR = "d:/vscode/Jin2025_AgingMouseBrain/data/cellchat"
MAX_CELLS_PER_TYPE = 1500
RANDOM_SEED = 42

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print(" Phase 4 v3: CellChat data prep (sparse, memory-efficient)")
print("=" * 60)

# ===========================================================================
# 1. Load
# ===========================================================================
print("\n[1] Loading hypothalamus h5ad...")
adata = sc.read_h5ad(HYPO_H5AD)
print(f"  Cells: {adata.n_obs:,}, Genes: {adata.n_vars:,}")

print("\n[2] Cell types:")
for ct, cnt in adata.obs['cell_type_major'].value_counts().items():
    print(f"  {ct}: {cnt:,}")

print("\n[3] Age distribution:")
for ag, cnt in adata.obs['age_group'].value_counts().items():
    print(f"  {ag}: {cnt:,}")

# ===========================================================================
# 2. Process each age group separately
# ===========================================================================
# Data is log2(x+1), convert to linear space BUT keep as float32 sparse
# CellChat will accept float data; exact integers not required

X_sparse = adata.X
if not sparse.issparse(X_sparse):
    X_sparse = sparse.csr_matrix(X_sparse)
print(f"\n[4] Sparse matrix shape: {X_sparse.shape}, type: {type(X_sparse).__name__}")

# Convert log2 to linear (sparse operation)
print("  Converting log2 -> linear (sparse)...")
gc.collect()

# For each age group
np.random.seed(RANDOM_SEED)

for age_group in ['Adult(2m)', 'Aged(18m)']:
    print(f"\n--- {age_group} ---")
    age_mask = adata.obs['age_group'] == age_group
    n_age = age_mask.sum()
    print(f"  Total cells: {n_age:,}")

    # Sample per cell type - collect indices
    all_sampled_idx = []
    cell_types_list = []

    for ct in sorted(adata.obs['cell_type_major'].unique()):
        ct_mask = (adata.obs['cell_type_major'] == ct) & age_mask
        ct_indices = np.where(ct_mask)[0]
        n_ct = len(ct_indices)

        if n_ct < 10:
            continue

        n_sample = min(n_ct, MAX_CELLS_PER_TYPE)
        if n_ct <= n_sample:
            sampled = ct_indices
        else:
            sampled = np.random.choice(ct_indices, size=n_sample, replace=False)

        all_sampled_idx.extend(sampled)
        cell_types_list.extend([ct] * len(sampled))
        print(f"    {ct}: {n_ct:,} -> {len(sampled):,}")

    # Sort indices for efficient sparse extraction
    all_sampled_idx = sorted(all_sampled_idx)
    n_total = len(all_sampled_idx)
    print(f"  Total sampled: {n_total:,}")

    # Extract sparse rows (batch-wise to avoid memory spikes)
    BATCH_SIZE = 2000
    batch_rows = []
    batch_barcodes = []

    for i in range(0, n_total, BATCH_SIZE):
        batch_idx = all_sampled_idx[i:i+BATCH_SIZE]
        batch = X_sparse[batch_idx].toarray()  # small batch -> dense
        # Convert log2(x+1) -> linear: 2^x - 1
        batch_linear = np.exp2(batch) - 1
        batch_linear = np.maximum(batch_linear, 0)
        batch_rows.append(sparse.csr_matrix(batch_linear))
        batch_barcodes.extend(adata.obs_names[batch_idx].tolist())

        if (i // BATCH_SIZE) % 10 == 0:
            gc.collect()
            pct = min(100, (i + BATCH_SIZE) * 100 // n_total)
            print(f"    Processing... {pct}%")

    # Combine batches
    print(f"  Combining {len(batch_rows)} batches...")
    subset_sparse = sparse.vstack(batch_rows, format='csr')
    del batch_rows
    gc.collect()

    print(f"  Subset matrix: {subset_sparse.shape}")

    # Use gene symbols (not Ensembl IDs) for CellChat
    gene_symbols = adata.var['gene_symbol'].astype(str).copy()
    # Handle duplicates: use unique suffix
    mask = gene_symbols.duplicated(keep=False)
    if mask.any():
        print(f"    Handling {mask.sum()} duplicate gene symbols...")
        # For duplicates, append a number to make unique
        dup_symbols = gene_symbols[mask]
        counts = {}
        new_names = []
        for s in dup_symbols:
            counts[s] = counts.get(s, 0) + 1
            new_names.append(f"{s}_{counts[s]}")
        gene_symbols[mask] = new_names
    # Ensure all unique
    assert not gene_symbols.duplicated().any(), "Duplicate gene symbols remain!"

    count_df = pd.DataFrame(
        subset_sparse.toarray(),
        index=batch_barcodes,
        columns=gene_symbols.values
    ).T  # genes x cells
    count_df.index.name = 'gene'

    age_label = 'young' if age_group == 'Adult(2m)' else 'aged'

    count_path = os.path.join(OUT_DIR, f"{age_label}_raw_counts.csv.gz")
    count_df.to_csv(count_path, compression='gzip')
    sz_mb = os.path.getsize(count_path) / 1024 / 1024
    print(f"  Counts: {count_path} ({sz_mb:.1f} MB)")

    # Save metadata
    meta = pd.DataFrame({
        'cell_label': batch_barcodes,
        'cell_type': cell_types_list,
        'age_group': age_group
    })
    meta_path = os.path.join(OUT_DIR, f"{age_label}_metadata.csv")
    meta.to_csv(meta_path, index=False)
    print(f"  Metadata: {meta_path} ({len(meta)} cells)")

    del count_df, subset_sparse
    gc.collect()

# ===========================================================================
# 3. Summary
# ===========================================================================
print("\n" + "=" * 60)
print(" Done!")
print("=" * 60)
for label in ['young', 'aged']:
    meta = pd.read_csv(os.path.join(OUT_DIR, f"{label}_metadata.csv"))
    print(f"\n{label}: {len(meta):,} cells, {meta['cell_type'].nunique()} cell types")
    for ct, cnt in meta['cell_type'].value_counts().items():
        print(f"  {ct}: {cnt}")

has_mg = ('Microglia' in pd.read_csv(os.path.join(OUT_DIR, "young_metadata.csv"))['cell_type'].values and
          'Microglia' in pd.read_csv(os.path.join(OUT_DIR, "aged_metadata.csv"))['cell_type'].values)
print(f"\nMicroglia present in both: {has_mg}")
