# =============================================================================
# Phase 2 备选方案: 用小h5ad直接跑Python处理，然后R加载结果
# =============================================================================
# 策略: 用Python anndata加载2GB下丘脑h5ad → 直接做QC/降维 → 导出RDS
# 绕开SeuratDisk兼容性问题
# =============================================================================

# Step 1: 用reticulate调Python
if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
library(reticulate)

# 使用已安装的Python
use_python("C:/Users/57265/AppData/Local/Programs/Python/Python312/python.exe", required = TRUE)

cat("============================================\n")
cat(" Refer: 由于R SeuratDisk兼容性问题,   \n")
cat(" 请在终端手动运行以下Python命令来预处理: \n")
cat("============================================\n\n")

cat('
# 在Git Bash中运行:
/c/Users/57265/AppData/Local/Programs/Python/Python312/python.exe -c "
import scanpy as sc
import pandas as pd
import numpy as np

# 加载下丘脑h5ad
print(\"Loading h5ad...\")
adata = sc.read_h5ad(\"d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_anatomical_52696.h5ad\")
print(f\"Cells: {adata.n_obs}, Genes: {adata.n_vars}\")

# 基础QC
print(\"\\nCell type distribution:\")
print(adata.obs[\"cell_type_major\"].value_counts())

print(\"\\nAge distribution:\")
print(adata.obs[\"age_group\"].value_counts())

print(\"\\nSex distribution:\")
print(adata.obs[\"sex_group\"].value_counts())

# 按细胞类型分组统计
print(\"\\nPer cell-type age breakdown:\")
for ct in adata.obs[\"cell_type_major\"].unique():
    sub = adata[adata.obs[\"cell_type_major\"] == ct]
    adult = (sub.obs[\"age_group\"] == \"Adult(2m)\").sum()
    aged = (sub.obs[\"age_group\"] == \"Aged(18m)\").sum()
    print(f\"  {ct}: total={sub.n_obs}, Adult={adult}, Aged={aged}\")

# 保存summary
summary = adata.obs.groupby([\"cell_type_major\", \"age_group\", \"sex_group\"]).size().reset_index(name=\"count\")
summary.to_csv(\"d:/vscode/Jin2025_AgingMouseBrain/results/hypo_anatomical_summary.csv\", index=False)
print(\"\\nSummary saved.\")
"
')