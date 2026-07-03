# GSE206861 单细胞转录组分析

## 📚 文献信息

| 项目 | 内容 |
|------|------|
| **标题** | Adamantinomatous craniopharyngioma cyst fluid can trigger inflammatory activation of microglia to damage the hypothalamic neurons by inducing the production of β-amyloid |
| **PMID** | [35525962](https://pubmed.ncbi.nlm.nih.gov/35525962/) |
| **期刊** | Journal of Neuroinflammation, 2022 |
| **GEO** | [GSE206861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206861) |

## 🧬 数据集概况

| 样本 | 描述 | 物种 |
|------|------|------|
| GSM6265811 | 小鼠下丘脑 + 人 ACP 囊液注射 | Mouse |
| GSM6265812 | 小鼠下丘脑 + PBS 注射（对照） | Mouse |
| GSM6265813 | 人 ACP 肿瘤组织 | Human |

数据格式：10x Genomics（barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz）

## 🚀 快速开始

### 1. 下载数据

**方法 A: 浏览器下载（推荐，速度最快）**

点击以下链接下载并放入 `GSE206861/` 目录：
```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206861/suppl/GSE206861_RAW.tar
```

**方法 B: R 代码下载**
```r
download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206861/suppl/GSE206861_RAW.tar",
  "GSE206861_RAW.tar",
  mode = "wb"
)
```

**方法 C: 命令行下载**
```bash
cd GSE206861/
curl -L -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206861/suppl/GSE206861_RAW.tar"
```

> ⚠️ 文件大小约 243MB，国内直连 NCBI 可能较慢，建议用浏览器或下载工具。

### 2. 安装 R 包

在 R/RStudio 中运行：
```r
source("scripts/01_install_packages.R")
```

**系统要求**：R >= 4.1.0，内存 >= 16GB（推荐 32GB）

### 3. 解压数据

```r
source("scripts/00_extract_data.R")
```

### 4. 运行分析

**一键运行**：
```r
source("scripts/run_all.R")
```

**分步运行**（推荐，可以逐步检查结果）：
```r
source("scripts/02_load_data_QC.R")              # 读入 + 质控
source("scripts/03_normalization_integration.R")  # 标准化 + 整合 + 聚类
source("scripts/04_cell_annotation.R")             # 细胞注释
source("scripts/05_differential_expression_enrichment.R")  # 差异 + 富集
source("scripts/06_advanced_analysis.R")           # 拟时序 + 细胞通讯
```

## 📂 项目结构

```
GSE206861/
├── GSE206861_RAW.tar            # 原始压缩包（需下载）
├── data/
│   ├── GSM6265811_Mouse_Cystic-fluid/  # 囊液注射组
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── GSM6265812_Mouse_Sham/          # PBS 对照组
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   ├── GSM6265813_Human_ACP/           # 人类肿瘤
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   └── matrix.mtx.gz
│   └── processed/                      # 中间处理对象 (*.rds)
├── scripts/
│   ├── 00_extract_data.R               # 数据解压整理
│   ├── 01_install_packages.R           # R 包安装
│   ├── 02_load_data_QC.R               # 读入 + 质控
│   ├── 03_normalization_integration.R  # 标准化 + 整合 + 聚类
│   ├── 04_cell_annotation.R            # 细胞注释
│   ├── 05_differential_expression_enrichment.R  # 差异 + 富集
│   ├── 06_advanced_analysis.R          # 拟时序 + 细胞通讯
│   └── run_all.R                       # 主运行脚本
├── results/                            # 结果图表
└── README.md                           # 本文件
```

## 📊 分析流程概览

```
数据下载 → 解压整理 → 质控过滤 → SCTransform 标准化
    → Harmony 整合（去批次）→ PCA/UMAP 降维 → 聚类
    → SingleR 自动注释 + Marker 手动验证
    → 组间差异分析（CF vs Sham）
    → GO/KEGG 富集分析 + GSEA
    → Monocle3 拟时序（小胶质细胞激活轨迹）
    → CellChat 细胞通讯（CD74-APP 互作验证）
```

## 🔑 可调参数

| 参数 | 脚本 | 说明 |
|------|------|------|
| `qc_mouse$nFeature_high` | 02_load_data_QC.R | 小鼠基因数上限 |
| `qc_mouse$percent.mt_max` | 02_load_data_QC.R | 线粒体比例阈值 |
| `CLUSTER_PARAMS$mouse_resolution` | run_all.R | 聚类分辨率 |
| `DE_PARAMS$logfc_threshold` | run_all.R | 差异基因 FC 阈值 |
| `run_advanced` | run_all.R | 是否运行进阶分析 |

## ⚠️ 注意事项

1. 首次运行 `01_install_packages.R` 可能需要 30-60 分钟安装所有依赖
2. 如果内存不足（<16GB），建议减小分析范围或使用服务器
3. 小鼠和人类样本分开分析，不做跨物种整合
4. 关键验证基因：*Cd68, Cd74, Il1b, Tnf* （小胶质细胞炎症）；*Npy, Fgfr2, Sst* （神经元损伤）
