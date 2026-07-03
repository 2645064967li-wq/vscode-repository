"""
=============================================================================
使用 abc_atlas_access 下载 Jin et al. 2025 衰老小鼠全脑数据
=============================================================================
安装: pip install "abc_atlas_access[notebooks] @ git+https://github.com/alleninstitute/abc_atlas_access.git"

使用方法:
    python scripts/download_with_abc_atlas.py

参考教程:
    https://alleninstitute.github.io/abc_atlas_access/notebooks/Zeng_Aging_Mouse_10x_snRNASeq_tutorial.html
=============================================================================
"""
import pandas as pd
import numpy as np
from pathlib import Path
import os

# ============================================================================
# 配置
# ============================================================================
PROJECT_DIR = Path(__file__).parent.parent if '__file__' in dir() else Path.cwd()
DOWNLOAD_BASE = PROJECT_DIR / "data" / "abc_atlas"

# ============================================================================
# 方法1: 使用 AbcProjectCache (推荐)
# ============================================================================
def download_with_cache():
    """使用abc_atlas_access包下载数据"""
    try:
        from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
    except ImportError:
        print("请先安装 abc_atlas_access:")
        print("  pip install \"abc_atlas_access[notebooks] @ git+https://github.com/alleninstitute/abc_atlas_access.git\"")
        return

    # 初始化缓存
    DOWNLOAD_BASE.mkdir(parents=True, exist_ok=True)
    abc_cache = AbcProjectCache.from_cache_dir(DOWNLOAD_BASE)

    # 加载最新manifest
    abc_cache.load_manifest('releases/20250531/manifest.json')

    print("可用数据目录:")
    for d in abc_cache.list_directories():
        print(f"  {d}")

    # ---- 1. 下载单细胞 metadata ----
    print("\n[1] 下载 Zeng-Aging-Mouse-10Xv3 metadata...")

    # 细胞metadata
    cell_meta = abc_cache.get_metadata_dataframe(
        directory='Zeng-Aging-Mouse-10Xv3',
        file_name='cell_metadata',
        dtype={'cell_label': str}
    )
    print(f"  细胞数: {len(cell_meta):,}")
    print(f"  列名: {list(cell_meta.columns)}")

    # 聚类信息
    cluster_info = abc_cache.get_metadata_dataframe(
        directory='Zeng-Aging-Mouse-10Xv3',
        file_name='cluster'
    )

    # 年龄差异基因
    degenes = abc_cache.get_metadata_dataframe(
        directory='Zeng-Aging-Mouse-WMB-taxonomy',
        file_name='aging_degenes'
    )
    print(f"  年龄差异基因: {len(degenes):,}")

    # 细胞聚类映射
    cell_cluster_mapping = abc_cache.get_metadata_dataframe(
        directory='Zeng-Aging-Mouse-WMB-taxonomy',
        file_name='cell_cluster_mapping_annotations'
    )

    # ---- 2. 筛选下丘脑细胞 ----
    print("\n[2] 筛选下丘脑细胞...")

    # 查找脑区相关列
    region_cols = [c for c in cell_meta.columns
                   if any(kw in c.lower() for kw in
                          ['parcellation', 'region', 'structure', 'anatomy', 'area', 'location'])]
    print(f"  脑区相关列: {region_cols}")

    # 尝试筛选
    hypo_mask = pd.Series(False, index=cell_meta.index)
    for col in region_cols:
        col_mask = cell_meta[col].astype(str).str.contains(
            'hypothalamus|hypothal|Arcuate|ARH|DMH|LHA|MPN|PVH|VMH|SCH|tanycyte|ependymal',
            case=False, na=False
        )
        if col_mask.sum() > 0:
            print(f"  '{col}'中: {col_mask.sum():,} 个下丘脑相关细胞")
            hypo_mask = hypo_mask | col_mask

    if hypo_mask.sum() > 0:
        hypo_cells = cell_meta[hypo_mask]
        print(f"\n  下丘脑细胞总数: {len(hypo_cells):,}")

        # 保存细胞标签
        hypo_dir = PROJECT_DIR / "data" / "hypothalamus"
        hypo_dir.mkdir(parents=True, exist_ok=True)
        hypo_cells.to_csv(hypo_dir / "hypothalamus_cell_metadata_python.csv")
        hypo_cells['cell_label'].to_csv(
            hypo_dir / "hypothalamus_cell_labels_python.csv", index=False
        )
        print(f"  已保存到: {hypo_dir}")
    else:
        print("  未直接找到下丘脑标记，请检查cell_metadata列")
        print(f"\n  所有列名: {list(cell_meta.columns)}")
        # 显示前几行供检查
        print("\n  前5行预览:")
        print(cell_meta.head())

    # ---- 3. 下载 MERFISH metadata ----
    print("\n[3] 下载 MERFISH metadata...")

    # CCF坐标
    ccf = abc_cache.get_metadata_dataframe(
        directory='MERFISH-C57BL6J-638850-CCF',
        file_name='ccf_coordinates',
        dtype={'cell_label': str}
    )
    print(f"  MERFISH细胞数: {len(ccf):,}")

    # 筛选下丘脑
    if 'parcellation_label' in ccf.columns:
        hypo_ccf = ccf[ccf['parcellation_label'].astype(str).str.contains(
            'hypothalamus|hypothal', case=False, na=False
        )]
        print(f"  下丘脑MERFISH细胞: {len(hypo_ccf):,}")

        if len(hypo_ccf) > 0:
            hypo_ccf.to_csv(hypo_dir / "hypothalamus_MERFISH_coordinates_python.csv")

    # ---- 4. 下载表达矩阵 (可选) ----
    print("\n[4] 表达矩阵路径 (可通过abc_atlas_access直接访问):")

    from abc_atlas_access.abc_atlas_cache.anndata_utils import get_gene_data

    # 获取特定基因的下丘脑表达
    # 示例: 获取AD相关基因
    ad_genes = ['App', 'Bace1', 'Psen1', 'Apoe', 'Trem2', 'Tyrobp',
                'Cd33', 'Clu', 'Tnf', 'Il1b', 'Npy', 'Agrp', 'Pomc']

    # 转为ENSEMBL ID
    gene_meta = abc_cache.get_metadata_dataframe(
        directory='Zeng-Aging-Mouse-10Xv3',
        file_name='gene'
    )

    print(f"\n  基因metadata列: {list(gene_meta.columns)}")

    print("\n========================================")
    print(" 数据下载完成!")
    print(f" 数据目录: {DOWNLOAD_BASE}")
    print("========================================")


# ============================================================================
# 方法2: 直接通过S3 URL下载 (无需abc_atlas_access包)
# ============================================================================
def download_direct_http():
    """使用requests直接从S3下载 (适用于无Python包的情况)"""
    import requests

    BASE_URL = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"

    files_to_download = {
        # 单细胞 metadata
        "single_cell/metadata/cell_metadata.csv":
            f"{BASE_URL}/metadata/Zeng-Aging-Mouse-10Xv3/20250131/cell_metadata.csv",
        "single_cell/metadata/cluster.csv":
            f"{BASE_URL}/metadata/Zeng-Aging-Mouse-10Xv3/20250131/cluster.csv",
        "single_cell/metadata/donor.csv":
            f"{BASE_URL}/metadata/Zeng-Aging-Mouse-10Xv3/20250131/donor.csv",
        # Taxonomy
        "single_cell/taxonomy/aging_degenes.csv":
            f"{BASE_URL}/metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/aging_degenes.csv",
        "single_cell/taxonomy/cell_cluster_mapping_annotations.csv":
            f"{BASE_URL}/metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/cell_cluster_mapping_annotations.csv",
        # MERFISH
        "spatial/MERFISH_638850/metadata/cell_metadata.csv":
            f"{BASE_URL}/metadata/MERFISH-C57BL6J-638850/20241115/cell_metadata.csv",
        "spatial/MERFISH_638850/metadata/ccf_coordinates.csv":
            f"{BASE_URL}/metadata/MERFISH-C57BL6J-638850-CCF/20231215/ccf_coordinates.csv",
    }

    data_dir = PROJECT_DIR / "data"

    for local_path, url in files_to_download.items():
        dest = data_dir / local_path
        dest.parent.mkdir(parents=True, exist_ok=True)

        print(f"下载: {local_path}")
        response = requests.get(url, stream=True)
        response.raise_for_status()

        with open(dest, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)


if __name__ == '__main__':
    download_with_cache()
