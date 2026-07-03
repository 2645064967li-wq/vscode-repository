# Analysis handoff: ageing mouse hypothalamus

Last updated: 2026-07-03

## Scope

Jin et al. (Nature, 2025) mouse ageing atlas. The current analysis focuses on
young (2-month) versus aged (18-month) hypothalamus, especially microglia and
tanycytes.

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

