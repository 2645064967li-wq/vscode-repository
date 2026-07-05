# Analysis handoff: ageing mouse hypothalamus

Last updated: 2026-07-05

## Scope

Jin et al. (Nature, 2025) mouse ageing atlas. The current analysis focuses on
young (2-month) versus aged (18-month) hypothalamus, especially microglia and
tanycytes.

## ARC/PVN/DMH spatial-mapping feasibility (2026-07-05)

- Allen CCF definitions: ARC = ARH (index 733); PVN = PVH plus descendants
  (720-730); DMH = DMH plus anterior/posterior/ventral parts (739-742).
- MERFISH taxonomy must use WMB cluster_alias 5283 for Microglia; ageing
  study aliases 840-844 are not valid for the MERFISH metadata.
- Exact cell_label matching is required. Labels with and without -1 are
  distinct cells and must not be collapsed.
- Adult-male MERFISH microglia counts: ARC 68, PVN 143, DMH 29 (240 total).
- Reference limitations: one adult male donor and a 550-gene measured panel.
  Use MERFISH for spatial label mapping only; test age/state/metabolism in the
  multi-donor full-transcriptome snRNA-seq object.
- Reproducible audit: scripts/12_spatial_region_feasibility.py; outputs:
  results/spatial_region_feasibility/.
- Sixteen target raw.h5ad sections (~2.49 GB) finished downloading to
  data/spatial/MERFISH_638850/expression/target_sections/. Files 19 and 26
  were redownloaded after detecting concurrent-write corruption.
- Extracted spatial object:
  results/spatial_region_feasibility/merfish_microglia_ARH_PVH_DMH_raw.h5ad
  (240 cells x 550 measured genes).

## ARC/PVN/DMH microglia label transfer (2026-07-05)

- Source ageing microglia object:
  data/hypothalamus_complete/microglia_HY_raw_3057.h5ad
  (3,057 cells x 32,285 genes).
- Transfer used 500 genes shared between measured MERFISH and snRNA-seq,
  per-modality log-normalization and gene-wise standardization, SelectKBest
  (100 genes), and class-balanced multinomial logistic regression.
- Validation:
  random stratified balanced accuracy 0.775, macro-F1 0.771;
  leave-one-section-out balanced accuracy 0.788, macro-F1 0.778.
- High-confidence rule: maximum probability >= 0.60 and probability margin
  over second class >= 0.20.
- Final high-confidence predictions: ARC/ARH 241, DMH 178, PVN/PVH 2,341;
  ambiguous 297.
- Outputs:
  results/spatial_region_mapping/microglia_3057_spatial_region_predictions.h5ad
  results/spatial_region_mapping/microglia_region_predictions.csv
  results/spatial_region_mapping/classifier_metrics.json
- Scripts:
  scripts/12_spatial_region_feasibility.py
  scripts/13_extract_target_merfish_expression.py
  scripts/15_transfer_microglia_regions.py

## Critical replication limitation

- Adult/young microglia are concentrated in two usable donors:
  female donor 538745 has 255 cells and male donor 549927 has 193 cells.
  Donors 618176 and 619870 have only 10 and 2 cells.
- Aged microglia are available from several donors, but this imbalance means
  region-by-age inference is not adequately replicated in the young group.
- Do not treat cells as independent replicates. Do not present cell-level
  Wilcoxon p-values as donor-level ageing evidence.
- Region-specific age DE, subtype shifts and metabolism must be labelled
  exploratory. Report donor-level values, effect sizes and sensitivity
  analyses; avoid strong significance claims.

## Next work for another AI agent

1. Recluster the 2,760 high-confidence predicted microglia using the full
   snRNA-seq transcriptome. Annotate homeostatic, DAM-like, interferon,
   inflammatory, phagocytic/lysosomal, lipid-associated and proliferative
   states using multiple markers.
2. Tabulate state proportions by donor, predicted region, age and sex.
   Preserve ambiguous cells as an excluded sensitivity group.
3. Perform exploratory donor-pseudobulk region and age comparisons only where
   both groups have sufficient cells. Report effect sizes and raw donor plots.
4. Score metabolic programs: glycolysis, OXPHOS, TCA, fatty-acid oxidation,
   cholesterol/lipid handling, mTOR/AMPK, lysosome/phagosome, ROS/glutathione,
   complement, interferon, TGFB and APOE-TREM2.
5. Regional CellChat must be restricted to cell populations with anatomical
   support: directly region-labelled ARC/ARH, PVH and DMH neurons plus
   high-confidence predicted microglia. Other glia must not be assigned to a
   nucleus without an explicit spatial mapping.
6. Run CellChat separately by sex and age when cell counts permit, using the
   existing CellChat 1.6.1 environment and parameters: CellChatDB.mouse,
   triMean, nboot 100, population.size FALSE, min.cells 30. Treat comparisons
   as exploratory because young donor replication is inadequate.
7. Prioritize CX3CL1-CX3CR1, CSF1/IL34-CSF1R, TGFB-TGFBR, complement,
   MIF, SPP1, GRN/PSAP, GAS6/PROS1-AXL/MERTK, CCL/CXCL and APOE-TREM2.

## Microglia state and metabolism analysis completed (2026-07-05)

- Scanpy repeatedly spent more than 10 minutes in Windows low-level
  initialization/PCA without producing output. The completed fallback uses
  sparse log-normalization, dispersion-ranked 2,000 HVGs, randomized
  TruncatedSVD (30 components), and KMeans selected across k=4..10 by
  silhouette. This avoids the unstable Scanpy/Numba path.
- Selected k=4, silhouette=0.102. This low value indicates a continuous
  microglial state spectrum rather than sharply separated subtypes.
- Cluster-level module annotation produced:
  Inflammatory-like 1,673 cells; Homeostatic-like 1,197;
  DAM/lipid-like 187. These are exploratory state-like labels, not definitive
  cell types. Interferon, phagolysosomal and proliferative modules remain as
  continuous scores even when they were not the top label of a cluster.
- Full scored object:
  results/microglia_region_analysis/microglia_3057_region_state_scored.h5ad
- Key tables:
  results/microglia_region_analysis/cluster_state_scores.csv
  results/microglia_region_analysis/state_counts_by_donor_region.csv
  results/microglia_region_analysis/module_scores_by_donor_region.csv
  results/microglia_region_analysis/microglia_cell_metadata_scored.csv
  results/microglia_region_analysis/gene_sets_used.json
- PCA figures:
  results/microglia_region_analysis/figures/
- Reproducible successful script:
  scripts/16b_microglia_lightweight_analysis.py
- scripts/16_microglia_region_state_metabolism.py is the Scanpy version that
  timed out on this machine and should not be used unless the environment is
  fixed.

## Regional neuron-microglia CellChat completed (2026-07-05)

- Scope is intentionally restricted to directly region-labelled GABAergic and
  glutamatergic neurons plus high-confidence spatially predicted microglia.
  Other glia were not assigned to nuclei.
- Inputs are balanced independently for every cell type using the minimum
  count across age x sex; each age contains equal female and male cells.
- Cells per age object:
  ARC 1,242 (28 microglia, 852 GABA, 362 Glut);
  PVN 1,526 (272 microglia, 950 GABA, 304 Glut);
  DMH 1,456 (32 microglia, 976 GABA, 448 Glut).
- CellChat 1.6.1, CellChatDB.mouse, triMean, nboot=100,
  population.size=FALSE, min.cells=10, seed=20260705.
- Run summaries:
  ARC young 258 LR / 43 pathways / weight 13.224;
  ARC aged 295 / 44 / 15.635.
  PVN young 385 / 44 / 12.676;
  PVN aged 415 / 51 / 9.876.
  DMH young 203 / 35 / 12.020;
  DMH aged 260 / 46 / 12.945.
- Broad pattern: ageing increased LR count in all three regions; total network
  strength increased in ARC and DMH but decreased in PVN.
- Largest microglia-related exploratory shifts include:
  DMH gain of microglial Ccl3-Ccr5 and stronger Tgfb1-Tgfbr1/Tgfbr2;
  ARC gain of microglial Ccl3/Ccl4-Ccr5 but loss of neuronal Tnr-Itga9/Itgb1
  and young ARC neuronal Csf1/Il34-Csf1r;
  aged ARC gain of neuronal Gas6-Mertk.
- Region-exclusive detection examples (not proof of biological exclusivity):
  young ARC Csf1/Il34-Csf1r; young DMH Gas6-Mertk and neuronal Tgfb2-TGFBR;
  aged PVN microglial Csf1-Csf1r; aged DMH neuronal Tgfb2/Gdf11-TGFBR and
  Gas6-Axl; aged ARC microglial Tgfb1-Acvr1b/Tgfbr2.
- All CellChat objects and tables:
  results/cellchat_region/ARC_ARH/
  results/cellchat_region/PVN_PVH/
  results/cellchat_region/DMH/
  results/cellchat_region/comparisons/
- Key comparison files:
  all_run_summaries.csv
  all_regions_age_lr_delta.csv
  all_regions_age_pathway_delta.csv
  microglia_priority_pathways_age_delta.csv
  microglia_region_specificity.csv
  microglia_region_exclusive_interactions.csv
- Scripts:
  scripts/17_prepare_region_cellchat.py
  scripts/18_run_region_cellchat.R
  scripts/19_compare_region_cellchat.py
- Critical interpretation: ARC and DMH have only 28 and 32 microglia per age
  object and young biological replication remains inadequate. CellChat
  differences are exploratory hypotheses, not donor-level statistical proof.

## Additional metabolism summaries

- Donor-level descriptive table:
  results/microglia_region_analysis/module_scores_age_descriptive.csv
- Exploratory aged-minus-young module effects:
  results/microglia_region_analysis/module_scores_age_effects_exploratory.csv
- These tables intentionally provide donor counts and effect directions
  without cell-level significance claims.

## Completed work

- Downloaded and organized Allen Brain Cell Atlas metadata and raw/log2 data.
- Extracted hypothalamic subsets and major cell-type annotations.
- Performed age- and sex-stratified differential-expression analyses.
- Performed GO/KEGG enrichment for key cell types.
- Prepared raw-count matrices for CellChat.
- Built separate young and aged CellChat objects.
- Compared global young/aged communication networks.
- Characterized microglial outgoing and incoming communication.
- Extracted bidirectional tanycyte-microglia ligand-receptor changes.

## Current CellChat objects

These large files are intentionally excluded from Git and must remain local:

- `results/cellchat/cellchat_young.rds`
- `results/cellchat/cellchat_aged.rds`

They can be regenerated with:

- `scripts/03_prepare_raw_counts_v5.py`
- `scripts/03_cellchat_analysis_v2.R`

## Latest outputs committed to Git

- `results/cellchat/differential_net_stats.csv`
- `results/cellchat/figures/`
- `results/cellchat/microglia/tanycyte_microglia_LR_changes.csv`
- `results/cellchat/microglia/Tany_to_MG_LR_changes.csv`
- `results/cellchat/microglia/MG_to_Tany_LR_changes.csv`

## Current biological observations

Age-associated tanycyte-microglia communication includes increased
`CSF1-CSF1R`, `CX3CL1-CX3CR1`, PTN-related, PSAP-related and `GRN-SORT1`
signals. Reduced or lost interactions include `PROS1-AXL`, `GAS6-AXL` and
`IL1A-IL1R1-IL1RAP`. These are computational predictions and require careful
statistical and biological validation before being presented as mechanisms.

## Recommended next steps

1. Correct and simplify `scripts/04_tanycyte_microglia_LR.R`; verify that the
   function uses the sender/receiver indices passed as arguments rather than
   global young-object indices.
2. Add explicit significance and minimum-probability filters before ranking
   ligand-receptor changes; avoid interpreting enormous fold changes caused by
   a zero denominator as quantitative effects.
3. Compare pathway-level information flow, not only individual LR pairs.
4. Validate key ligands and receptors against differential expression and the
   fraction of expressing cells in each age group.
5. Produce publication-quality figures and a concise methods/results report.

## Resume prompt

When switching assistants, use:

> Read `Jin2025_AgingMouseBrain/HANDOFF.md`, inspect the current scripts and
> results, and continue from the recommended next steps without rerunning
> completed heavy preprocessing unless necessary.


## Complete hypothalamus CellChat update (2026-07-04)

A replacement CellChat analysis was completed from the corrected 82,431-cell HY object. Four sex-stratified, cell-type-balanced objects were generated under `results/cellchat_complete/`. Each contains 9,218 cells and 11 common cell types, using CellChat 1.6.1, CellChatDB.mouse, triMean, 100 bootstraps and population.size=FALSE.

The primary Tanycyte→Microglia result is reduced GDF11–TGFBR1/ACVR2A signaling in both sexes. TGFB2–TGFBR1/TGFBR2 decreases in females but is stable in males. Young Microglia have insufficient donor replication for formal donor-level inference, so CellChat findings remain exploratory.

Read `results/cellchat_complete/ANALYSIS_REPORT.md` before continuing. Do not use the older `results/cellchat/cellchat_young.rds` and `cellchat_aged.rds` for final interpretation.
