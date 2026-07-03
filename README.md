# Single-cell transcriptomics analysis workspace

This repository contains two independent analysis projects. Raw sequencing
matrices and serialized R/Python objects are intentionally excluded from Git;
each project README documents how to obtain or regenerate them.

## Projects

### 1. GSE206861

Analysis of hypothalamic responses to adamantinomatous craniopharyngioma cyst
fluid, with emphasis on microglial activation, differential expression and
Monocle3 trajectories.

- Project documentation: [`GSE206861/README.md`](GSE206861/README.md)
- Reproducible scripts: [`GSE206861/scripts/`](GSE206861/scripts/)
- Compact result tables and figures: [`GSE206861/results/`](GSE206861/results/)

### 2. Jin2025_AgingMouseBrain

Analysis of the Jin et al. (Nature, 2025) healthy-ageing mouse brain atlas,
focused on the hypothalamus, microglia, tanycytes and age-dependent CellChat
communication.

- Project documentation: [`Jin2025_AgingMouseBrain/README.md`](Jin2025_AgingMouseBrain/README.md)
- Current handoff/checkpoint: [`Jin2025_AgingMouseBrain/HANDOFF.md`](Jin2025_AgingMouseBrain/HANDOFF.md)
- Reproducible scripts: [`Jin2025_AgingMouseBrain/scripts/`](Jin2025_AgingMouseBrain/scripts/)
- Compact result tables and figures: [`Jin2025_AgingMouseBrain/results/`](Jin2025_AgingMouseBrain/results/)

## Data policy

Files such as `.h5ad`, `.h5seurat`, `.rds`, raw 10x matrices and downloaded
archives are not committed because they are large and can be regenerated from
the public source datasets. This repository is intended to preserve analysis
logic, provenance and interpretable outputs.

