"""
============================================================================
Phase 4 v4: 从raw counts h5ad高效提取下丘脑细胞raw counts
使用h5py直接读取sparse matrix，分批处理避免内存溢出
============================================================================
"""
import h5py
import numpy as np
from scipy import sparse
import pandas as pd
import os, gc, sys
import warnings
warnings.filterwarnings('ignore')

RAW_H5AD = "d:/vscode/Jin2025_AgingMouseBrain/Zeng-Aging-Mouse-10Xv3-raw.h5ad"
BARCODE_FILE = "d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_all_barcodes.txt"
CELLMETA_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/cell_metadata.csv"
CLUSTER_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/cluster.csv"
DONOR_CSV = "d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata/donor.csv"
OUT_DIR = "d:/vscode/Jin2025_AgingMouseBrain/data/cellchat"
MAX_CELLS_PER_TYPE = 1500
RANDOM_SEED = 42

os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print(" Phase 4 v4: Extract RAW counts from h5ad via h5py")
print("=" * 60)

# ===========================================================================
# 1. Load barcodes & metadata
# ===========================================================================
print("\n[1] Loading hypothalamus barcodes...")
with open(BARCODE_FILE) as f:
    hypo_barcodes = [line.strip() for line in f if line.strip()]
hypo_set = set(hypo_barcodes)
print(f"  {len(hypo_barcodes):,} barcodes")

print("\n[2] Loading metadata...")
cluster_info = pd.read_csv(CLUSTER_CSV)
cluster_to_name = dict(zip(cluster_info['cluster_alias'], cluster_info['cluster_name']))

donor_info = pd.read_csv(DONOR_CSV)
donor_age_map = dict(zip(donor_info['donor_label'], donor_info['donor_age_category']))

# Load cell metadata - filter to hypothalamus cells
print("  Loading cell_metadata (chunked)...")
needed_cols = ['cell_label', 'cluster_alias', 'donor_label', 'donor_age_category']
chunks = []
for chunk in pd.read_csv(CELLMETA_CSV, usecols=needed_cols, chunksize=200000):
    mask = chunk['cell_label'].isin(hypo_set)
    if mask.any():
        chunks.append(chunk[mask])
    gc.collect()
hypo_meta = pd.concat(chunks, ignore_index=True)
print(f"  Matched: {len(hypo_meta):,} cells")

# ===========================================================================
# 3. Assign cell types & age
# ===========================================================================
print("\n[3] Assigning cell types & age...")

def assign_cell_type_major(cluster_name):
    if not isinstance(cluster_name, str):
        return "Unclassified"
    cn = cluster_name
    if "Microglia" in cn: return "Microglia"
    if "Astro" in cn and "Ependymal" not in cn.lower(): return "Astrocyte"
    if "Tanycyte" in cn: return "Tanycyte"
    if "Ependymal" in cn or "ependymal" in cn: return "Ependymal"
    if "COP" in cn: return "COP"
    if "NFOL" in cn: return "NFOL"
    if "MFOL" in cn: return "MFOL"
    if "MOL" in cn: return "MOL"
    if "OPC" in cn: return "OPC"
    if "Endo" in cn and ("Endo-" in cn or "Endo_" in cn): return "Endothelial"
    if "SMC" in cn: return "SMC"
    if "VLMC" in cn: return "VLMC"
    if "Peri" in cn: return "Pericyte"
    if "BAM" in cn: return "BAM"
    if "ABC" in cn: return "ABC"
    if "T cells" in cn or "Tcell" in cn: return "T_cell"
    if "Glut" in cn and "Gaba" not in cn and "GABA" not in cn: return "Glutamatergic_Neuron"
    if "Gaba" in cn or "GABA" in cn: return "GABAergic_Neuron"
    if "Dopa" in cn: return "Dopaminergic_Neuron"
    if "Sero" in cn: return "Serotonergic_Neuron"
    if "Chol" in cn: return "Cholinergic_Neuron"
    if "Hist" in cn: return "Histaminergic_Neuron"
    return "Unclassified"

hypo_meta['cluster_name'] = hypo_meta['cluster_alias'].map(cluster_to_name)
hypo_meta['cell_type_major'] = hypo_meta['cluster_name'].apply(assign_cell_type_major)
hypo_meta['age_group'] = hypo_meta['donor_age_category'].apply(
    lambda x: 'Aged(18m)' if x == 'aged' else ('Adult(2m)' if x == 'adult' else 'Unknown'))
# Fill missing age from donor.csv
missing = hypo_meta['age_group'] == 'Unknown'
if missing.sum() > 0:
    hypo_meta.loc[missing, 'age_group'] = (
        hypo_meta.loc[missing, 'donor_label'].map(donor_age_map)
        .apply(lambda x: 'Aged(18m)' if x == 'aged' else ('Adult(2m)' if x == 'adult' else 'Unknown'))
    )

print("  Cell type distribution:")
for ct, cnt in hypo_meta['cell_type_major'].value_counts().items():
    print(f"    {ct}: {cnt:,}")

# ===========================================================================
# 4. Sample cells
# ===========================================================================
print(f"\n[4] Sampling (max {MAX_CELLS_PER_TYPE} per type per age)...")
np.random.seed(RANDOM_SEED)

all_sampled_bc = []
all_sampled_ct = []
all_sampled_age = []

for age_group in ['Adult(2m)', 'Aged(18m)']:
    age_cells = hypo_meta[hypo_meta['age_group'] == age_group]
    for ct in sorted(age_cells['cell_type_major'].unique()):
        ct_cells = age_cells[age_cells['cell_type_major'] == ct]
        n = min(len(ct_cells), MAX_CELLS_PER_TYPE)
        if n < 10:
            continue
        if len(ct_cells) <= n:
            sampled = ct_cells
        else:
            sampled = ct_cells.sample(n=n, random_state=RANDOM_SEED)
        all_sampled_bc.extend(sampled['cell_label'].tolist())
        all_sampled_ct.extend([ct] * len(sampled))
        all_sampled_age.extend([age_group] * len(sampled))
        print(f"    {age_group} {ct}: {len(ct_cells):,} -> {len(sampled):,}")

n_total = len(all_sampled_bc)
sampled_set = set(all_sampled_bc)
print(f"\n  Total sampled: {n_total:,}")

# ===========================================================================
# 5. Read raw h5ad and extract counts
# ===========================================================================
print("\n[5] Reading raw h5ad and extracting counts...")
f = h5py.File(RAW_H5AD, 'r')

# Read obs names (cell barcodes) - stored in cell_label
obs_index = list(f['obs']['cell_label'][:])
# Convert bytes to str if needed
if isinstance(obs_index[0], bytes):
    obs_index = [x.decode('utf-8') for x in obs_index]
print(f"  h5ad obs: {len(obs_index):,} cells")

# Build mapping: barcode -> index in h5ad
bc_to_idx = {bc: i for i, bc in enumerate(obs_index)}
valid_bc = [bc for bc in all_sampled_bc if bc in bc_to_idx]
valid_idx = [bc_to_idx[bc] for bc in valid_bc]
print(f"  Valid in h5ad: {len(valid_bc):,} / {n_total:,}")

# Read var names (gene symbols, categorical)
gene_symbols_cat = f['var']['gene_symbol']
gene_codes = gene_symbols_cat['codes'][:]
gene_categories = gene_symbols_cat['categories'][:]
# Decode categories if bytes
if isinstance(gene_categories[0], bytes):
    gene_categories = [x.decode('utf-8') for x in gene_categories]
var_names = [gene_categories[c] for c in gene_codes]
print(f"  Genes: {len(var_names):,}")
print(f"  Gene name examples: {var_names[:5]}")

# Read X matrix (sparse CSR format expected)
x_group = f['X']
x_shape = tuple(x_group.attrs['shape'])
print(f"  X shape: {x_shape}")

# For CSR format: read indptr, indices, data separately (DO NOT load all data)
x_indptr = x_group['indptr'][:]
x_indices = x_group['indices'][:]
x_data = x_group['data'][:]
n_rows, n_cols = x_shape
print(f"  X shape: {n_rows} x {n_cols}")
print(f"  nnz: {len(x_data):,} (from {len(x_data)} data elements)")

# ===========================================================================
# 6. Extract rows using CSR indexing (row-by-row, memory efficient)
# ===========================================================================
print("\n[6] Extracting per age group (CSR row indexing)...")

# Match valid_bc back to age groups
bc_to_age = dict(zip(all_sampled_bc, all_sampled_age))
bc_to_ct = dict(zip(all_sampled_bc, all_sampled_ct))

def extract_row(row_idx):
    """Extract a single row from CSR sparse matrix as dense array"""
    start = x_indptr[row_idx]
    end = x_indptr[row_idx + 1]
    row = np.zeros(n_cols, dtype=np.float32)
    row[x_indices[start:end]] = x_data[start:end]
    return row

for age_group in ['Adult(2m)', 'Aged(18m)']:
    label = 'young' if age_group == 'Adult(2m)' else 'aged'
    print(f"\n  --- {age_group} ({label}) ---")

    # Filter indices for this age group
    age_valid = [(bc, idx) for bc, idx in zip(valid_bc, valid_idx) if bc_to_age.get(bc) == age_group]
    age_bc = [x[0] for x in age_valid]
    age_idx = [x[1] for x in age_valid]
    print(f"    Cells: {len(age_idx):,}")

    # Extract rows one at a time (much less memory)
    all_rows = []
    for i, row_idx in enumerate(age_idx):
        all_rows.append(extract_row(row_idx))
        if (i + 1) % 1000 == 0:
            pct = (i + 1) * 100 // len(age_idx)
            print(f"      Extracting... {pct}%")
            gc.collect()

    full_mat = np.vstack(all_rows)
    del all_rows
    gc.collect()

    print(f"    Matrix: {full_mat.shape}")

    # Build DataFrame (genes x cells)
    count_df = pd.DataFrame(full_mat.T, index=var_names, columns=batch_bcs)
    count_df.index.name = 'gene'

    count_path = os.path.join(OUT_DIR, f"{label}_raw_counts.csv.gz")
    count_df.to_csv(count_path, compression='gzip')
    sz = os.path.getsize(count_path) / 1024 / 1024
    print(f"    Saved: {count_path} ({sz:.1f} MB)")

    # Metadata
    meta = pd.DataFrame({
        'cell_label': batch_bcs,
        'cell_type': [bc_to_ct.get(bc, 'Unknown') for bc in batch_bcs],
        'age_group': age_group
    })
    meta_path = os.path.join(OUT_DIR, f"{label}_metadata.csv")
    meta.to_csv(meta_path, index=False)
    print(f"    Metadata: {meta_path} ({len(meta)} cells)")

    del full_mat, count_df
    gc.collect()

f.close()

# ===========================================================================
# 7. Summary
# ===========================================================================
print("\n" + "=" * 60)
print(" Phase 4 v4 complete!")
print("=" * 60)

for label in ['young', 'aged']:
    meta = pd.read_csv(os.path.join(OUT_DIR, f"{label}_metadata.csv"))
    print(f"\n{label}: {len(meta):,} cells, {meta['cell_type'].nunique()} types")
    for ct, cnt in meta['cell_type'].value_counts().items():
        print(f"  {ct}: {cnt}")

# Check Microglia
ym = pd.read_csv(os.path.join(OUT_DIR, "young_metadata.csv"))
am = pd.read_csv(os.path.join(OUT_DIR, "aged_metadata.csv"))
print(f"\nMicroglia: Young={('Microglia' in ym['cell_type'].values)}, Aged={('Microglia' in am['cell_type'].values)}")
