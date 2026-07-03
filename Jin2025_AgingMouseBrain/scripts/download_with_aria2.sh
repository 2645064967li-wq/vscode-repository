#!/bin/bash
# =============================================================================
# 使用 aria2c 多线程下载 (速度比 curl 快 10-50 倍)
# =============================================================================
# 首先安装 aria2:
#   Git Bash: 从 https://github.com/aria2/aria2/releases 下载 aria2-*-win-64bit.zip
#             解压后将 aria2c.exe 放到 C:\Windows\System32 或 PATH 中
#   或使用 Chocolatey: choco install aria2
#
# 如果有代理, 设置:
#   export https_proxy=http://127.0.0.1:7890
# =============================================================================

set -e

BASE_URL="https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
DATA_DIR="d:/vscode/Jin2025_AgingMouseBrain/data"

echo "============================================"
echo " aria2c 高速下载 Jin et al. 2025 数据集"
echo "============================================"
echo ""

# 检查 aria2c 是否可用
if ! command -v aria2c &> /dev/null; then
    echo "❌ aria2c 未找到!"
    echo ""
    echo "安装方法:"
    echo "  1. 下载 Windows 版本: https://github.com/aria2/aria2/releases"
    echo "  2. 解压 aria2c.exe 到 PATH 目录"
    echo "  3. 或在 Git Bash 中使用 Chocolatey: choco install aria2"
    echo ""
    echo "备选: 直接用浏览器+IDM下载, 链接见 download_urls.txt"
    exit 1
fi

mkdir -p "${DATA_DIR}/single_cell/metadata/views"
mkdir -p "${DATA_DIR}/single_cell/taxonomy"
mkdir -p "${DATA_DIR}/single_cell/expression"
mkdir -p "${DATA_DIR}/spatial/MERFISH_638850/metadata/views"
mkdir -p "${DATA_DIR}/spatial/MERFISH_638850/expression"

# ============================================================================
# Phase 1: Metadata (小文件, ~1.6 GB)
# ============================================================================
echo "[Phase 1/2] 下载 Metadata (~1.6 GB)..."
echo ""

# 创建临时URL列表文件
URL_LIST=$(mktemp)
trap "rm -f ${URL_LIST}" EXIT

cat > "${URL_LIST}" << 'EOF'
# 单细胞 metadata (v20250131)
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-10Xv3/20250131/cell_metadata.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata
  out=cell_metadata.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-10Xv3/20250131/cell_annotation_colors.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata
  out=cell_annotation_colors.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-10Xv3/20250131/cell_cluster_annotations.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata
  out=cell_cluster_annotations.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-10Xv3/20250131/cluster.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata
  out=cluster.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-10Xv3/20250131/donor.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata
  out=donor.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-10Xv3/20250131/library.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata
  out=library.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-10Xv3/20250131/value_sets.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/metadata
  out=value_sets.csv
# Taxonomy
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/aging_degenes.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/taxonomy
  out=aging_degenes.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/cell_cluster_mapping_annotations.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/taxonomy
  out=cell_cluster_mapping_annotations.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/cell_cross_mapping_annotations.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/taxonomy
  out=cell_cross_mapping_annotations.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/cluster_mapping.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/taxonomy
  out=cluster_mapping.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130/cluster_mapping_pivot.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/single_cell/taxonomy
  out=cluster_mapping_pivot.csv
# MERFISH metadata
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/MERFISH-C57BL6J-638850/20241115/cell_metadata.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/spatial/MERFISH_638850/metadata
  out=cell_metadata.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/MERFISH-C57BL6J-638850-CCF/20231215/ccf_coordinates.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/spatial/MERFISH_638850/metadata
  out=ccf_coordinates.csv
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/metadata/MERFISH-C57BL6J-638850/20241115/gene.csv
  dir=d:/vscode/Jin2025_AgingMouseBrain/data/spatial/MERFISH_638850/metadata
  out=gene.csv
EOF

# 16线程并行下载metadata文件
aria2c -x 16 -s 16 -j 10 --continue=true \
    --summary-interval=30 \
    -i "${URL_LIST}"

echo ""
echo "✅ Phase 1 完成! Metadata 下载完毕"
echo ""

# ============================================================================
# Phase 2: 表达矩阵 (大文件)
# ============================================================================
echo "[Phase 2/2] 下载表达矩阵..."
echo ""

echo "选择要下载的文件:"
echo "  1) 仅 scRNA-seq log2 (~15 GB) [推荐, 最常用]"
echo "  2) scRNA-seq log2 + raw (~31 GB)"
echo "  3) 全部包括 MERFISH imputed (~78 GB)"
echo "  4) 跳过"
echo ""
read -p "请选择 (1/2/3/4): " choice

case $choice in
    1)
        echo "下载 scRNA-seq log2 表达矩阵 (~15 GB)..."
        aria2c -x 16 -s 16 --continue=true --summary-interval=60 \
            -d "${DATA_DIR}/single_cell/expression" \
            -o "Zeng-Aging-Mouse-10Xv3-log2.h5ad" \
            "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-log2.h5ad"
        ;;
    2)
        echo "下载 scRNA-seq log2 + raw (~31 GB)..."
        aria2c -x 16 -s 16 --continue=true --summary-interval=60 \
            -d "${DATA_DIR}/single_cell/expression" \
            "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-log2.h5ad" \
            "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-raw.h5ad"
        ;;
    3)
        echo "下载全部表达矩阵 (~78 GB)..."
        aria2c -x 16 -s 16 --continue=true --summary-interval=60 \
            -d "${DATA_DIR}/single_cell/expression" \
            "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-log2.h5ad" \
            "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-raw.h5ad"
        aria2c -x 16 -s 16 --continue=true --summary-interval=60 \
            -d "${DATA_DIR}/spatial/MERFISH_638850/expression" \
            "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/MERFISH-C57BL6J-638850-imputed/20240831/C57BL6J-638850-imputed-log2.h5ad"
        ;;
    4)
        echo "跳过表达矩阵下载"
        ;;
    *)
        echo "无效选择, 跳过"
        ;;
esac

echo ""
echo "============================================"
echo " 下载完成!"
echo "============================================"
du -sh "${DATA_DIR}"
