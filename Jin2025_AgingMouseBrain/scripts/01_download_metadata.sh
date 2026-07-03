#!/bin/bash
# =============================================================================
# 下载 Jin et al. 2025 衰老小鼠全脑数据 - Metadata
# =============================================================================
# 使用说明:
#   1. 在 Git Bash 中运行: bash scripts/01_download_metadata.sh
#   2. 全部metadata约 1-2 GB, 下载时间取决于网络
#   3. 下载完成后会在 data/ 目录下生成对应的metadata文件
# =============================================================================

set -e

BASE_URL="https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="${PROJECT_DIR}/data"

echo "============================================"
echo " Jin et al. 2025 - Metadata 下载"
echo " 存储位置: ${DATA_DIR}"
echo "============================================"
echo ""

# ============================================================================
# 1. 单细胞转录组 Metadata (Zeng-Aging-Mouse-10Xv3) - 20250131 最新版
# ============================================================================
SC_META_DIR="${DATA_DIR}/single_cell/metadata"
SC_META_URL="${BASE_URL}/metadata/Zeng-Aging-Mouse-10Xv3/20250131"

mkdir -p "${SC_META_DIR}"

echo "[1/3] 下载单细胞 metadata (Zeng-Aging-Mouse-10Xv3)..."
echo "      版本: 20250131"

SC_FILES=(
    "cell_metadata.csv"
    "cell_annotation_colors.csv"
    "cell_cluster_annotations.csv"
    "cluster.csv"
    "donor.csv"
    "library.csv"
    "value_sets.csv"
    "views/example_genes_all_cells_expression.csv"
)

for file in "${SC_FILES[@]}"; do
    file_dir="$(dirname "${file}")"
    if [ "${file_dir}" != "." ]; then
        mkdir -p "${SC_META_DIR}/${file_dir}"
    fi

    echo "  -> 下载 ${file}..."
    curl -L --retry 3 --retry-delay 5 \
        -o "${SC_META_DIR}/${file}" \
        "${SC_META_URL}/${file}" \
        --progress-bar
done

echo "  ✓ 单细胞 metadata 完成!"
echo ""

# ============================================================================
# 2. 单细胞 Taxonomy Metadata (Zeng-Aging-Mouse-WMB-taxonomy)
# ============================================================================
TAX_META_DIR="${DATA_DIR}/single_cell/taxonomy"
TAX_URL="${BASE_URL}/metadata/Zeng-Aging-Mouse-WMB-taxonomy/20241130"

mkdir -p "${TAX_META_DIR}"

echo "[2/3] 下载 Taxonomy metadata (WMB-taxonomy)..."
echo "      版本: 20241130"

TAX_FILES=(
    "aging_degenes.csv"
    "cell_cluster_mapping_annotations.csv"
    "cell_cross_mapping_annotations.csv"
    "cluster_mapping.csv"
    "cluster_mapping_pivot.csv"
)

for file in "${TAX_FILES[@]}"; do
    echo "  -> 下载 ${file}..."
    curl -L --retry 3 --retry-delay 5 \
        -o "${TAX_META_DIR}/${file}" \
        "${TAX_URL}/${file}" \
        --progress-bar
done

echo "  ✓ Taxonomy metadata 完成!"
echo ""

# ============================================================================
# 3. MERFISH 空间转录组 Metadata
# ============================================================================
MERFISH_META_DIR="${DATA_DIR}/spatial/MERFISH_638850/metadata"
MERFISH_URL="${BASE_URL}/metadata/MERFISH-C57BL6J-638850/20241115"
MERFISH_CCF_URL="${BASE_URL}/metadata/MERFISH-C57BL6J-638850-CCF/20231215"
MERFISH_GENE_URL="${BASE_URL}/metadata/MERFISH-C57BL6J-638850-imputed/20240831"

mkdir -p "${MERFISH_META_DIR}"

echo "[3/3] 下载 MERFISH 空间转录组 metadata..."
echo "      版本: 20241115 (cell metadata), 20231215 (CCF)"

MERFISH_FILES=(
    "cell_metadata.csv"
    "gene.csv"
    "views/cell_metadata_with_cluster_annotation.csv"
)

for file in "${MERFISH_FILES[@]}"; do
    file_dir="$(dirname "${file}")"
    if [ "${file_dir}" != "." ]; then
        mkdir -p "${MERFISH_META_DIR}/${file_dir}"
    fi

    echo "  -> 下载 ${file}..."
    curl -L --retry 3 --retry-delay 5 \
        -o "${MERFISH_META_DIR}/${file}" \
        "${MERFISH_URL}/${file}" \
        --progress-bar
done

# CCF coordinates
echo "  -> 下载 ccf_coordinates.csv (CCF空间坐标)..."
curl -L --retry 3 --retry-delay 5 \
    -o "${MERFISH_META_DIR}/ccf_coordinates.csv" \
    "${MERFISH_CCF_URL}/ccf_coordinates.csv" \
    --progress-bar

# Imputed gene list
echo "  -> 下载 imputed gene.csv..."
curl -L --retry 3 --retry-delay 5 \
    -o "${MERFISH_META_DIR}/gene_imputed.csv" \
    "${MERFISH_GENE_URL}/gene.csv" \
    --progress-bar

echo "  ✓ MERFISH metadata 完成!"
echo ""

# ============================================================================
# 显示下载文件摘要
# ============================================================================
echo "============================================"
echo " 下载完成! 文件列表:"
echo "============================================"
echo ""
echo "单细胞转录组 metadata:"
du -sh "${SC_META_DIR}" 2>/dev/null || true
echo ""
echo "Taxonomy metadata:"
du -sh "${TAX_META_DIR}" 2>/dev/null || true
echo ""
echo "MERFISH metadata:"
du -sh "${MERFISH_META_DIR}" 2>/dev/null || true
echo ""
echo "下一步:"
echo "  bash scripts/02_download_expression_matrices.sh  # 下载表达矩阵(较大)"
echo "  Rscript scripts/filter_hypothalamus.R             # 提取下丘脑数据"
