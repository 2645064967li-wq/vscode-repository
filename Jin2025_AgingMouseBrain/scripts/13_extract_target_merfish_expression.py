"""Extract ARC/PVH/DMH microglia from downloaded MERFISH raw section files."""
from pathlib import Path
import sys
import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/spatial/MERFISH_638850/expression/target_sections"
REGIONS = ROOT / "results/spatial_region_feasibility/merfish_microglia_ARH_PVH_DMH.csv"
OUTPUT = ROOT / "results/spatial_region_feasibility/merfish_microglia_ARH_PVH_DMH_raw.h5ad"
EXPECTED_BYTES = {
    "18": 81777248, "19": 132028203, "24": 120794735, "25": 153782735,
    "26": 174444290, "27": 117366374, "30": 195688673, "31": 189594214,
    "32": 195821005, "33": 174277850, "35": 164027374, "36": 236094292,
    "37": 161791758, "38": 260613863, "39": 122542266, "40": 194487147,
}

regions = pd.read_csv(REGIONS, dtype={"cell_label": "string"}).set_index("cell_label", drop=False)
pieces, problems = [], []
for section, expected_size in EXPECTED_BYTES.items():
    path = INPUT / f"C57BL6J-638850.{section}-raw.h5ad"
    if not path.exists():
        problems.append(f"missing: {path.name}")
        continue
    if path.stat().st_size != expected_size:
        problems.append(f"size mismatch: {path.name} {path.stat().st_size} != {expected_size}")
        continue
    a = ad.read_h5ad(path)
    labels = pd.Index(a.obs_names.astype(str), dtype="string")
    wanted = labels.isin(regions.index)
    if wanted.any():
        sub = a[wanted, :].copy()
        sub.obs["cell_label"] = sub.obs_names.astype(str)
        meta = regions.loc[sub.obs["cell_label"].astype("string")].copy()
        meta.index = sub.obs_names
        for column in ["analysis_region", "x", "y", "z", "parcellation_index",
                       "cluster_alias", "average_correlation_score"]:
            sub.obs[column] = meta[column].to_numpy()
        sub.obs["source_section"] = section
        pieces.append(sub)
    a.file.close()

if problems:
    print("\n".join(problems), file=sys.stderr)
    raise SystemExit(2)
if not pieces:
    raise RuntimeError("No target cells found")
combined = ad.concat(pieces, join="inner", merge="same", index_unique=None)
combined.obs_names_make_unique()
found = set(combined.obs["cell_label"].astype(str))
missing = sorted(set(regions.index.astype(str)) - found)
if missing:
    pd.Series(missing, name="cell_label").to_csv(
        OUTPUT.with_name("missing_target_merfish_cells.csv"), index=False)
    raise RuntimeError(f"Only found {len(found)}/{len(regions)} target cells")
combined.uns["region_definition"] = {
    "ARC_ARH": "Allen CCF ARH (ARC)",
    "PVN_PVH": "Allen CCF PVH and descendants (PVN)",
    "DMH": "Allen CCF DMH and descendants",
}
combined.uns["evidence_scope"] = (
    "Spatial reference is one adult male C57BL/6J mouse; it does not directly test ageing.")
combined.write_h5ad(OUTPUT, compression="gzip")
print(f"wrote {OUTPUT}")
print(f"shape: {combined.n_obs} cells x {combined.n_vars} genes")
print(combined.obs["analysis_region"].value_counts().sort_index().to_string())