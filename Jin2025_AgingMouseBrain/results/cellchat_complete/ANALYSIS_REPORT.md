# Complete hypothalamus CellChat analysis

## Scope and methods

- Input: complete HY raw-count object, 82,431 cells × 32,285 genes.
- CellChat: v1.6.1 with `CellChatDB.mouse` (2,019 interactions; Secreted Signaling, ECM-Receptor and Cell-Cell Contact).
- Design: Young Female, Aged Female, Young Male and Aged Male analyzed separately.
- Common cell types required at least 30 cells in every stratum. Eleven types passed: Astrocyte, Endothelial, Ependymal, GABAergic Neuron, Glutamatergic Neuron, Microglia, MOL, OPC, Pericyte, SMC and Tanycyte.
- Each cell type was downsampled to the minimum count across the four strata. Each CellChat object therefore contains exactly 9,218 cells with identical cell-type composition.
- Input was library-size normalized to 10,000 counts per cell and log1p transformed.
- Parameters: `type="triMean"`, `nboot=100`, `population.size=FALSE`, `min.cells=30`, seed 20260704.
- CellChat probabilities and permutation P values are cell-level exploratory evidence, not donor-level inferential statistics.

## Global network summary

| Group | Significant LR records | Pathways | Total communication weight |
|---|---:|---:|---:|
| Young Female | 5,134 | 79 | 701.16 |
| Aged Female | 5,544 | 85 | 701.63 |
| Young Male | 5,703 | 80 | 768.67 |
| Aged Male | 5,474 | 81 | 695.18 |

Female ageing increased the number of detected LR records while leaving aggregate weight nearly unchanged. Male ageing reduced both LR records and total communication weight. Pooled-sex interpretation is therefore inappropriate.

## Tanycyte–Microglia TGFB/TGFBR findings

All detected TGFB-family interactions were Tanycyte → Microglia; no significant Microglia → Tanycyte TGFB-family interaction was detected.

| Sex | Interaction | Young | Aged | Change |
|---|---|---:|---:|---:|
| Female | TGFB2–TGFBR1/TGFBR2 | 0.1231 | 0.0948 | -0.0283 |
| Female | TGFB2–ACVR1/TGFBR1 | 0.1000 | 0.0547 | -0.0453 |
| Female | GDF11–TGFBR1/ACVR2A | 0.1156 | 0.0411 | -0.0744 |
| Male | TGFB2–TGFBR1/TGFBR2 | 0.1041 | 0.1052 | +0.0011 |
| Male | TGFB2–ACVR1/TGFBR1 | 0.0715 | 0.0608 | -0.0107 |
| Male | GDF11–TGFBR1/ACVR2A | 0.1326 | 0.0452 | -0.0874 |
| Male | GDF11–TGFBR1/ACVR2B | 0.1246 | 0 | lost |

The balanced analysis does not support a general ageing-associated enhancement of TGFBR signaling. It supports reduced GDF11-related signaling in both sexes, reduced TGFB2 signaling in females, and receptor-complex-specific stability or mild reduction in males.

## Donor-level expression validation

- Tanycyte `Gdf11` expression decreased with age in both sexes.
- Male Tanycyte `Tgfb2` pseudobulk CPM decreased; female median `Tgfb2` CPM was approximately stable.
- Microglial `Tgfbr1`, `Acvr1`, `Acvr2a` and `Acvr2b` generally showed lower median expression in aged animals; male `Tgfbr2` was an exception.
- Critical limitation: after requiring at least 20 Microglia, each young sex has only one informative donor. Therefore the Microglia age comparison cannot support donor-level statistical significance and must remain exploratory.

## Main outputs

- Four CellChat RDS objects: one under each group directory.
- `comparisons/lr_age_difference_*.csv`: aligned LR-level age differences.
- `comparisons/pathway_age_difference_*.csv`: pathway-level differences.
- `comparisons/tanycyte_microglia_TGFB_family_age_difference.csv`: focused TGFB-family table.
- `comparisons/tgfb_donor_expression_validation.csv`: donor-level expression and detection fractions.
- `comparisons/global_network_age_sex_summary.pdf`: global network summary.
- `comparisons/tanycyte_microglia_TGFB_age_difference.pdf`: focused communication plot.
- `comparisons/tgfb_donor_expression_validation.pdf`: donor expression validation.

## Interpretation boundary

These results prioritize candidate interactions. They do not demonstrate physical ligand–receptor binding, signaling activity or causality. Spatial colocalization and independent biological replication are required before mechanistic claims.
