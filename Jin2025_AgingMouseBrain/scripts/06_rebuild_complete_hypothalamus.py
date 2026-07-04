"""Rebuild and audit the complete HY subset from Allen raw counts.

This script deliberately keeps Allen's cluster assignments. It only derives
coarse labels from the official cluster names, exports donor-level summaries,
and creates a neuron-only embedding for GABA/glutamate annotation QC.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import h5py
from scipy import sparse
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors


def coarse_cell_type(cluster_name: object) -> str:
    name = "" if pd.isna(cluster_name) else str(cluster_name)
    low = name.lower()
    if "microglia" in low:
        return "Microglia"
    if "bam" in low:
        return "BAM"
    if "tanycyte" in low:
        return "Tanycyte"
    if "ependymal" in low:
        return "Ependymal"
    if "astro" in low:
        return "Astrocyte"
    if "endo" in low:
        return "Endothelial"
    if "peri" in low:
        return "Pericyte"
    if "vlmc" in low:
        return "VLMC"
    if "smc" in low:
        return "SMC"
    if "nfol" in low:
        return "NFOL"
    if "mfol" in low:
        return "MFOL"
    if "mol" in low:
        return "MOL"
    if "opc" in low:
        return "OPC"
    if "cop" in low:
        return "COP"
    if "gaba" in low:
        return "GABAergic_Neuron"
    if "glut" in low:
        return "Glutamatergic_Neuron"
    if "dopa" in low:
        return "Dopaminergic_Neuron"
    if "sero" in low:
        return "Serotonergic_Neuron"
    if "chol" in low:
        return "Cholinergic_Neuron"
    if "t cell" in low or "tcell" in low:
        return "T_cell"
    if "_dc" in low or "dc_" in low:
        return "Dendritic_Cell"
    return "Unclassified"


def dense_vector(adata: ad.AnnData, gene: str) -> np.ndarray:
    x = adata[:, gene].X
    return np.asarray(x.toarray() if sparse.issparse(x) else x).ravel()


def read_selected_csr(path: Path, selected_rows: np.ndarray) -> sparse.csr_matrix:
    """Read selected CSR rows without AnnData's large temporary int64 array."""
    with h5py.File(path, "r") as handle:
        group = handle["X"]
        shape = tuple(int(x) for x in group.attrs["shape"])
        source_indptr = group["indptr"][:]
        lengths = source_indptr[selected_rows + 1] - source_indptr[selected_rows]
        target_indptr = np.empty(len(selected_rows) + 1, dtype=np.int64)
        target_indptr[0] = 0
        np.cumsum(lengths, out=target_indptr[1:])
        nnz = int(target_indptr[-1])
        target_indices = np.empty(nnz, dtype=group["indices"].dtype)
        target_data = np.empty(nnz, dtype=group["data"].dtype)
        boundaries = np.flatnonzero(np.diff(selected_rows) != 1) + 1
        for block in np.split(np.arange(len(selected_rows)), boundaries):
            first, last = int(block[0]), int(block[-1])
            source_start = int(source_indptr[selected_rows[first]])
            source_end = int(source_indptr[selected_rows[last] + 1])
            target_start = int(target_indptr[first])
            target_end = int(target_indptr[last + 1])
            target_indices[target_start:target_end] = group["indices"][source_start:source_end]
            target_data[target_start:target_end] = group["data"][source_start:source_end]
    return sparse.csr_matrix((target_data, target_indices, target_indptr), shape=(len(selected_rows), shape[1]))

def marker_audit(neurons: ad.AnnData, out_csv: Path) -> pd.DataFrame:
    genes = [g for g in ["Slc17a6", "Slc17a7", "Slc32a1", "Gad1", "Gad2"] if g in neurons.var_names]
    rows = []
    for cell_type in ["GABAergic_Neuron", "Glutamatergic_Neuron"]:
        mask = neurons.obs["cell_type_major"].astype(str).eq(cell_type).to_numpy()
        for gene in genes:
            values = dense_vector(neurons, gene)[mask]
            rows.append({
                "cell_type": cell_type,
                "gene": gene,
                "n_cells": int(mask.sum()),
                "mean_log1p": float(values.mean()),
                "fraction_detected": float((values > 0).mean()),
            })
    result = pd.DataFrame(rows)
    result.to_csv(out_csv, index=False)
    return result


def umap_mixing_audit(neurons: ad.AnnData, out_csv: Path) -> pd.DataFrame:
    xy = neurons.obsm["X_umap"]
    labels = neurons.obs["cell_type_major"].astype(str).to_numpy()
    indices = NearestNeighbors(n_neighbors=16).fit(xy).kneighbors(return_distance=False)[:, 1:]
    opposite = (labels[indices] != labels[:, None]).mean(axis=1)
    neurons.obs["opposite_label_fraction_15nn"] = opposite
    result = (
        neurons.obs.assign(cell_type=labels)
        .groupby("cell_type", observed=True)["opposite_label_fraction_15nn"]
        .agg(n_cells="size", mean="mean", median="median", q90=lambda x: x.quantile(0.9))
        .reset_index()
    )
    result["fraction_gt_0.5"] = [float((opposite[labels == x] > 0.5).mean()) for x in result["cell_type"]]
    result["umap_silhouette"] = silhouette_score(xy, labels, sample_size=min(20000, len(labels)), random_state=0)
    result.to_csv(out_csv, index=False)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=Path("Zeng-Aging-Mouse-10Xv3-raw.h5ad"))
    parser.add_argument("--metadata", type=Path, default=Path("data/single_cell/metadata/cell_metadata.csv"))
    parser.add_argument("--clusters", type=Path, default=Path("data/single_cell/metadata/cluster.csv"))
    parser.add_argument("--out", type=Path, default=Path("data/hypothalamus_complete"))
    parser.add_argument("--results", type=Path, default=Path("results/reanalysis_2026"))
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.results.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.metadata, low_memory=False).set_index("cell_label", drop=False)
    hy_meta = metadata.loc[metadata["anatomical_division_label"].eq("HY")].copy()
    clusters = pd.read_csv(args.clusters)
    cluster_map = dict(zip(clusters["cluster_alias"].astype(str), clusters["cluster_name"]))
    hy_meta["cluster_name"] = hy_meta["cluster_alias"].astype(str).map(cluster_map)
    hy_meta["cell_type_major"] = hy_meta["cluster_name"].map(coarse_cell_type)
    hy_meta["age_group"] = hy_meta["donor_age_category"].map({"adult": "Adult(2m)", "aged": "Aged(18m)"})
    hy_meta["sex_group"] = hy_meta["donor_sex"].map({"F": "Female", "M": "Male"})

    source = sc.read_h5ad(args.raw, backed="r")
    keep = source.obs_names.isin(hy_meta.index)
    if int(keep.sum()) != len(hy_meta):
        raise RuntimeError(f"Raw/metadata mismatch: raw matched {keep.sum():,}, metadata HY {len(hy_meta):,}")
    selected_rows = np.flatnonzero(keep)
    selected_obs = source.obs.iloc[selected_rows].copy()
    selected_var = source.var.copy()
    source.file.close()
    matrix = read_selected_csr(args.raw, selected_rows)
    hypo = ad.AnnData(X=matrix, obs=selected_obs, var=selected_var)
    common_columns = [c for c in hy_meta.columns if c not in hypo.obs.columns]
    hypo.obs = hypo.obs.join(hy_meta[common_columns], how="left")
    for c in ["cluster_name", "cell_type_major", "age_group", "sex_group"]:
        hypo.obs[c] = hy_meta.loc[hypo.obs_names, c].to_numpy()
    hypo.obs_names_make_unique()
    hypo.var_names = hypo.var["gene_symbol"].astype(str).to_numpy()
    hypo.var_names_make_unique()
    hypo.layers["counts"] = hypo.X.copy()
    hypo.write_h5ad(args.out / "hypothalamus_HY_raw_82431.h5ad", compression="gzip")

    composition = (
        hypo.obs.groupby(["donor_label", "age_group", "sex_group", "cell_type_major"], observed=True)
        .size().rename("n_cells").reset_index()
    )
    composition.to_csv(args.results / "donor_celltype_counts.csv", index=False)

    sc.pp.normalize_total(hypo, target_sum=1e4)
    sc.pp.log1p(hypo)
    neuron_mask = hypo.obs["cell_type_major"].isin(["GABAergic_Neuron", "Glutamatergic_Neuron"])
    neurons = hypo[neuron_mask].copy()
    sc.pp.highly_variable_genes(neurons, n_top_genes=3000, flavor="seurat", batch_key="library_label")
    sc.tl.pca(neurons, n_comps=50, use_highly_variable=True, random_state=0)
    sc.pp.neighbors(neurons, n_neighbors=30, n_pcs=30)
    sc.tl.umap(neurons, min_dist=0.3, random_state=0)

    marker_audit(neurons, args.results / "neuron_marker_audit.csv")
    mixing = umap_mixing_audit(neurons, args.results / "neuron_umap_mixing.csv")
    neurons.obs.to_csv(args.results / "neuron_metadata_with_umap_mixing.csv")

    sc.settings.set_figure_params(dpi=150, frameon=False)
    fig = sc.pl.umap(
        neurons,
        color=["cell_type_major", "Slc17a6", "Slc32a1", "Gad1", "Gad2", "donor_label"],
        ncols=2,
        return_fig=True,
        show=False,
    )
    fig.savefig(args.results / "neuron_umap_marker_validation.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    neurons.write_h5ad(args.out / "hypothalamus_neurons_log1p_umap.h5ad", compression="gzip")

    report = {
        "n_hypothalamus": int(hypo.n_obs),
        "n_genes": int(hypo.n_vars),
        "n_neurons": int(neurons.n_obs),
        "cell_type_counts": hypo.obs["cell_type_major"].value_counts().to_dict(),
        "mixing": mixing.to_dict(orient="records"),
    }
    (args.results / "reanalysis_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
