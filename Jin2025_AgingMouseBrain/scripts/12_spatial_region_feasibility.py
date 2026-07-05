"""Audit MERFISH microglia coverage in ARH (ARC), PVH (PVN), and DMH.

This script uses metadata only. It intentionally does not require the ~47 GB
imputed expression matrix. Large CSV files are streamed in chunks.
"""

from pathlib import Path
import json
import re

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPATIAL = ROOT / "data" / "spatial"
META = SPATIAL / "MERFISH_638850" / "metadata"
TAXONOMY = ROOT / "data" / "single_cell" / "metadata" / "cluster.csv"
OUT = ROOT / "results" / "spatial_region_feasibility"
OUT.mkdir(parents=True, exist_ok=True)


def descendants(terms: pd.DataFrame, root_label: str) -> set[str]:
    """Return a label and descendants, traversing MBA identifiers."""
    root_ids = set(terms.loc[terms["label"].eq(root_label), "identifier"].dropna())
    if not root_ids:
        raise ValueError(f"Missing parcellation root: {root_label}")
    keep_ids = root_ids
    while True:
        child_ids = set(
            terms.loc[terms["parent_identifier"].isin(keep_ids), "identifier"].dropna()
        )
        updated = keep_ids | child_ids
        if updated == keep_ids:
            return set(terms.loc[terms["identifier"].isin(keep_ids), "label"].dropna())
        keep_ids = updated


terms = pd.read_csv(SPATIAL / "parcellation_term.csv")
terms["parcellation_index"] = pd.to_numeric(terms["graph_order"], errors="coerce")

roots = {
    "ARC_ARH": "AllenCCF-Ontology-2017-223",
    "PVN_PVH": "AllenCCF-Ontology-2017-38",
    "DMH": "AllenCCF-Ontology-2017-830",
}

index_to_region: dict[int, str] = {}
region_terms = []
for region, root in roots.items():
    labels = descendants(terms, root)
    selected = terms.loc[
        terms["label"].isin(labels) & terms["parcellation_index"].notna()
    ].copy()
    # Integer graph_order rows are the Allen CCF annotation entries. Decimal
    # rows are ABC visualization substructures and are not used by CCF CSV.
    selected = selected.loc[selected["parcellation_index"] % 1 == 0]
    for idx in selected["parcellation_index"].astype(int):
        if idx in index_to_region and index_to_region[idx] != region:
            raise ValueError(f"Overlapping CCF index {idx}: {index_to_region[idx]} vs {region}")
        index_to_region[idx] = region
    selected["analysis_region"] = region
    region_terms.append(selected)

region_terms_df = pd.concat(region_terms, ignore_index=True)
region_terms_df[
    ["analysis_region", "graph_order", "label", "acronym", "name", "parent_identifier"]
].to_csv(OUT / "region_definition.csv", index=False)

# MERFISH cluster_alias uses the WMB taxonomy, whereas the ageing single-cell
# object uses study-specific aliases. Resolve WMB aliases through cross-mapping.
cross_path = ROOT / "data" / "single_cell" / "taxonomy" / "cell_cross_mapping_annotations.csv"
mapped = []
for chunk in pd.read_csv(
    cross_path,
    usecols=["wmb_cluster_alias", "wmb_cluster_name", "wmb_subclass_name"],
    chunksize=500_000,
):
    mapped.append(
        chunk.loc[
            chunk["wmb_subclass_name"].str.contains("Microglia", case=False, na=False),
            ["wmb_cluster_alias", "wmb_cluster_name", "wmb_subclass_name"],
        ]
    )
microglia_clusters = pd.concat(mapped).drop_duplicates().copy()
microglia_aliases = set(microglia_clusters["wmb_cluster_alias"].astype(int))
microglia_clusters.to_csv(OUT / "microglia_cluster_aliases.csv", index=False)

# MERFISH and CCF labels are matched exactly; suffixes identify distinct cells.
microglia_labels: dict[str, tuple[int, float]] = {}
metadata_counts = {"all_cells": 0, "microglia_cells": 0}
donors, sexes, genotypes = set(), set(), set()
for chunk in pd.read_csv(
    META / "cell_metadata.csv",
    usecols=[
        "cell_label", "cluster_alias", "average_correlation_score",
        "donor_label", "donor_sex", "donor_genotype",
    ],
    chunksize=500_000,
):
    metadata_counts["all_cells"] += len(chunk)
    donors.update(chunk["donor_label"].dropna().astype(str).unique())
    sexes.update(chunk["donor_sex"].dropna().astype(str).unique())
    genotypes.update(chunk["donor_genotype"].dropna().astype(str).unique())
    m = chunk.loc[chunk["cluster_alias"].isin(microglia_aliases)].copy()
    metadata_counts["microglia_cells"] += len(m)
    normalized = m["cell_label"].astype(str)
    for label, alias, score in zip(
        normalized, m["cluster_alias"].astype(int), m["average_correlation_score"]
    ):
        microglia_labels[label] = (alias, float(score))

selected_chunks = []
ccf_counts = {"all_rows": 0, "microglia_matched": 0}
for chunk in pd.read_csv(
    META / "ccf_coordinates.csv",
    usecols=["cell_label", "x", "y", "z", "parcellation_index"],
    chunksize=500_000,
):
    ccf_counts["all_rows"] += len(chunk)
    chunk["cell_label"] = chunk["cell_label"].astype(str)
    m = chunk.loc[chunk["cell_label"].isin(microglia_labels)].copy()
    ccf_counts["microglia_matched"] += len(m)
    m["analysis_region"] = m["parcellation_index"].map(index_to_region)
    m = m.loc[m["analysis_region"].notna()].copy()
    if not m.empty:
        info = m["cell_label"].map(microglia_labels)
        m["cluster_alias"] = info.map(lambda x: x[0])
        m["average_correlation_score"] = info.map(lambda x: x[1])
        selected_chunks.append(m)

selected = pd.concat(selected_chunks, ignore_index=True) if selected_chunks else pd.DataFrame()
selected.to_csv(OUT / "merfish_microglia_ARH_PVH_DMH.csv", index=False)

if selected.empty:
    counts = pd.DataFrame(columns=["analysis_region", "cluster_alias", "n_cells"])
else:
    counts = (
        selected.groupby(["analysis_region", "cluster_alias"], observed=True)
        .size().rename("n_cells").reset_index()
    )
counts.to_csv(OUT / "microglia_counts_by_region_cluster.csv", index=False)

summary = {
    "region_definition": {
        region: sorted(int(i) for i, value in index_to_region.items() if value == region)
        for region in roots
    },
    "microglia_aliases": sorted(microglia_aliases),
    "metadata": metadata_counts,
    "ccf": ccf_counts,
    "donors": sorted(donors),
    "sexes": sorted(sexes),
    "genotypes": sorted(genotypes),
    "selected_region_microglia": (
        selected["analysis_region"].value_counts().sort_index().astype(int).to_dict()
        if not selected.empty else {}
    ),
    "expression_matrix_present": any((META.parent / "expression").glob("*.h5ad")),
}
(OUT / "feasibility_summary.json").write_text(
    json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
)

print(json.dumps(summary, indent=2, ensure_ascii=False))
