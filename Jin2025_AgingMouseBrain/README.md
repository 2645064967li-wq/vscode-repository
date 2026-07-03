# Jin et al. 2025 - 小鼠健康衰老的全脑细胞类型特异性转录组特征

## 📚 文献信息

| 项目 | 内容 |
|------|------|
| **标题** | Brain-wide cell-type-specific transcriptomic signatures of healthy ageing in mice |
| **作者** | Kelly Jin, Zizhen Yao, Cindy T. J. van Velthoven, Eitan S. Kaplan, et al. (Hongkui Zeng lab) |
| **期刊** | *Nature*, Volume 638, Issue 8049, Pages 182–196 |
| **发表日期** | 2025年1月1日 (在线) / 2025年2月6日 (出版) |
| **DOI** | [10.1038/s41586-024-08350-8](https://doi.org/10.1038/s41586-024-08350-8) |
| **PMID** | 39743592 |
| **数据平台** | [Allen Brain Cell (ABC) Atlas](https://alleninstitute.github.io/abc_atlas_access/) |
| **许可协议** | CC BY 4.0 (单细胞) / CC BY-NC 4.0 |

## 🧬 数据集概况

### 单细胞转录组 (scRNA-seq): Zeng-Aging-Mouse-10Xv3

| 项目 | 内容 |
|------|------|
| **技术平台** | 10x Genomics v3 (single-nucleus RNA-seq) |
| **细胞数量** | ~1,200,000 高质量单细胞转录组 |
| **年龄组** | 年轻成年 (2个月) 和 老年 (18个月) |
| **性别** | 雄性和雌性 |
| **脑区覆盖** | 前脑 (forebrain)、中脑 (midbrain)、后脑 (hindbrain) |
| **聚类结果** | 847 个细胞聚类, 至少 14 个年龄偏倚聚类 |
| **差异基因** | 2,449 个年龄相关差异表达基因 (age-DE genes) |
| **S3路径** | `s3://allen-brain-cell-atlas/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/` |
| **Metadata** | `s3://allen-brain-cell-atlas/metadata/Zeng-Aging-Mouse-10Xv3/20250131/` |
| **Taxonomy** | `s3://allen-brain-cell-atlas/metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/` |

### 空间转录组 (MERFISH): MERFISH-C57BL6J-638850 (Allen Institute)

| 项目 | 内容 |
|------|------|
| **技术平台** | MERFISH (Multiplexed Error-Robust FISH) |
| **基因面板** | 500 个基因 (原始测量) + 8,460 个基因 (imputed) |
| **细胞数量** | ~4,000,000 细胞 |
| **CCF注册** | Allen Mouse Brain Common Coordinate Framework v3 |
| **S3路径** | `s3://allen-brain-cell-atlas/expression_matrices/MERFISH-C57BL6J-638850-imputed/20240831/` |
| **CCF坐标** | `s3://allen-brain-cell-atlas/metadata/MERFISH-C57BL6J-638850-CCF/20231215/` |
| **Metadata** | `s3://allen-brain-cell-atlas/metadata/MERFISH-C57BL6J-638850/20241115/` |

### 空间转录组 (MERFISH): Zhuang-ABCA-1~4 (Zhuang Lab)

| 项目 | 内容 |
|------|------|
| **技术平台** | MERFISH |
| **基因面板** | ~1,100 个基因 |
| **细胞数量** | ~9,000,000 细胞 (4只动物合计) |
| **S3路径** | `s3://allen-brain-cell-atlas/expression_matrices/Zhuang-ABCA-{1-4}/20230830/` |

## 🎯 分析重点：下丘脑 (Hypothalamus)

该研究发现**下丘脑是衰老的核心枢纽** (hypothalamic ageing hub)：
- 第三脑室周围的细胞类型对衰老最敏感
- 包括：tanycytes（伸展细胞）、ependymal cells（室管膜细胞）、以及弓状核(ARC)、背内侧核(DMN)、室旁核(PVN)中的特定神经元
- 这些细胞表现出神经元功能下降和免疫反应增强的双重特征

### 下丘脑相关的CCF结构ID

用于从全脑数据中提取下丘脑区域:

```
CCFv3 Structure IDs for Hypothalamus (参考):
- HY (Hypothalamus): 1097
- 子区域包括: ARH, DMH, LHA, MPN, PH, PVH, SCH, SO, VMH, etc.
```

## 🚀 快速开始

### 方案A: 使用 abc_atlas_access Python包 (推荐)

```bash
# 安装
pip install "abc_atlas_access[notebooks] @ git+https://github.com/alleninstitute/abc_atlas_access.git"

# 参考官方教程
# https://alleninstitute.github.io/abc_atlas_access/notebooks/Zeng_Aging_Mouse_10x_snRNASeq_tutorial.html
```

### 方案B: 直接下载 (使用curl/wget)

```bash
# 第一步: 下载metadata (推荐先运行，文件较小 ~1-2GB)
bash scripts/01_download_metadata.sh

# 第二步: 下载表达矩阵 (文件较大 ~30-50GB+)
bash scripts/02_download_expression_matrices.sh
```

### 方案C: R语言下载

```r
# 在R中运行
source("scripts/01_download_metadata.R")
```

## 📂 项目结构

```
Jin2025_AgingMouseBrain/
├── README.md                           # 本文件
├── data/
│   ├── single_cell/
│   │   ├── metadata/                   # scRNA-seq 元数据
│   │   │   ├── cell_metadata.csv       # 细胞注释信息
│   │   │   ├── cluster.csv             # 聚类信息
│   │   │   ├── donor.csv               # 供体信息
│   │   │   ├── library.csv             # 文库信息
│   │   │   └── aging_degenes.csv       # 年龄差异基因
│   │   └── expression/                 # 表达矩阵 (h5ad)
│   │       ├── Zeng-Aging-Mouse-10Xv3-log2.h5ad
│   │       └── Zeng-Aging-Mouse-10Xv3-raw.h5ad
│   └── spatial/
│       ├── MERFISH_638850/
│       │   ├── metadata/               # MERFISH 细胞元数据
│       │   │   ├── cell_metadata.csv
│       │   │   ├── ccf_coordinates.csv # CCFv3空间坐标
│       │   │   └── gene.csv
│       │   └── expression/             # 表达矩阵 (h5ad)
│       │       └── C57BL6J-638850-imputed-log2.h5ad
│       └── Zhuang_ABCA/                # Zhuang Lab MERFISH (可选)
│           └── metadata/
├── scripts/
│   ├── 01_download_metadata.sh         # 下载metadata
│   ├── 02_download_expression_matrices.sh  # 下载表达矩阵
│   ├── 01_download_metadata.R          # R版本下载脚本
│   └── filter_hypothalamus.R           # 下丘脑数据提取
└── results/                            # 分析结果
```

## 📊 数据下载量估计

| 数据类型 | 文件 | 大小(估计) | 必需? |
|----------|------|------------|-------|
| scRNA-seq metadata | cell_metadata.csv 等 | ~300 MB | ✅ 是 |
| scRNA-seq taxonomy | cluster_mapping 等 | ~200 MB | ✅ 是 |
| scRNA-seq log2 | Zeng-Aging-Mouse-10Xv3-log2.h5ad | ~15 GB | ✅ 推荐 |
| scRNA-seq raw | Zeng-Aging-Mouse-10Xv3-raw.h5ad | ~16 GB | ⚠️ 可选 |
| MERFISH CCF | ccf_coordinates.csv | ~600 MB | ✅ 是 |
| MERFISH metadata | cell_metadata.csv 等 | ~500 MB | ✅ 是 |
| MERFISH imputed | C57BL6J-638850-imputed-log2.h5ad | ~47 GB | ⚠️ 可选 |
| MERFISH sections | 各切片表达矩阵 | ~20 GB | ❌ 按需 |
| **总计(最小)** | | ~2 GB | |
| **总计(完整)** | | ~80-100 GB | |

## 🔬 关键分析思路

1. **下丘脑单细胞图谱构建**: 从全脑 ~120万细胞中筛选下丘脑区域
2. **衰老差异分析**: 比较年轻(2m) vs 老年(18m)小鼠下丘脑各细胞类型
3. **空间分布验证**: 用MERFISH空间数据验证关键基因在下丘脑的空间表达模式
4. **与AD关联**: 将衰老相关基因与已知AD风险基因取交集
5. **细胞通讯**: 分析下丘脑衰老过程中细胞间配体-受体互作变化

## ⚠️ 注意事项

1. **下载量巨大**: 完整数据集约 80-100 GB，请确保有足够磁盘空间
2. **网络建议**: 由于数据托管在AWS S3 us-west-2，国内下载可能需要代理/VPN
3. **内存需求**: 全脑数据很大，建议 >= 32GB RAM；下丘脑子集后16GB可处理
4. **Python依赖**: 处理h5ad文件需要 `anndata` + `scanpy` (或R中的 `Seurat` v5 + `SeuratDisk`)
5. **CCF注册**: 空间数据基于Allen CCFv3坐标系统，需要用Allen Brain Atlas API进行结构映射
