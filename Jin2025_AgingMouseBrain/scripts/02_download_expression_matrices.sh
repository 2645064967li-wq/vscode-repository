#!/bin/bash
# =============================================================================
# 下载 Jin et al. 2025 衰老小鼠全脑数据 - 表达矩阵
# =============================================================================
# 使用说明:
#   1. 确保已先运行 01_download_metadata.sh
#   2. 在 Git Bash 中运行: bash scripts/02_download_expression_matrices.sh
#   3. 下载量约 80 GB, 需要足够磁盘空间
#   4. 可以选择性下载 (通过设置下面的标志)
# =============================================================================

set -e

BASE_URL="https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="${PROJECT_DIR}/data"

# ============================================================================
# 下载开关 (根据需要调整)
# ============================================================================
DOWNLOAD_SC_RAW=false      # scRNA-seq raw counts (约16 GB, 可选)
DOWNLOAD_SC_LOG2=true      # scRNA-seq log2 normalized (约15 GB, 推荐)
DOWNLOAD_MERFISH_IMPUTED=true  # MERFISH imputed genes (约47 GB)
DOWNLOAD_MERFISH_SECTIONS=false  # MERFISH 各切片原始数据 (约20 GB)

echo "============================================"
echo " Jin et al. 2025 - 表达矩阵下载"
echo " 存储位置: ${DATA_DIR}"
echo "============================================"
echo ""
echo "下载配置:"
echo "  scRNA-seq log2:      ${DOWNLOAD_SC_LOG2}"
echo "  scRNA-seq raw:       ${DOWNLOAD_SC_RAW}"
echo "  MERFISH imputed:     ${DOWNLOAD_MERFISH_IMPUTED}"
echo "  MERFISH sections:    ${DOWNLOAD_MERFISH_SECTIONS}"
echo ""

# 计算预估大小
EST_SIZE=0
[ "${DOWNLOAD_SC_LOG2}" = true ] && EST_SIZE=$((EST_SIZE + 15))
[ "${DOWNLOAD_SC_RAW}" = true ] && EST_SIZE=$((EST_SIZE + 16))
[ "${DOWNLOAD_MERFISH_IMPUTED}" = true ] && EST_SIZE=$((EST_SIZE + 47))
[ "${DOWNLOAD_MERFISH_SECTIONS}" = true ] && EST_SIZE=$((EST_SIZE + 20))

echo "预估下载量: ~${EST_SIZE} GB"
echo ""
read -p "确认开始下载? (y/n): " -n 1 -r
echo ""
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "取消下载"
    exit 0
fi

# ============================================================================
# 1. 单细胞转录组 表达矩阵
# ============================================================================
SC_EXPR_DIR="${DATA_DIR}/single_cell/expression"
SC_EXPR_URL="${BASE_URL}/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130"
mkdir -p "${SC_EXPR_DIR}"

if [ "${DOWNLOAD_SC_LOG2}" = true ]; then
    echo ""
    echo "[1/4] 下载 scRNA-seq log2 表达矩阵..."
    echo "      文件: Zeng-Aging-Mouse-10Xv3-log2.h5ad"
    echo "      估计大小: ~15 GB"
    echo "      预计时间: 30分钟 - 2小时 (取决于网速)"
    echo ""
    curl -L --retry 5 --retry-delay 30 \
        -o "${SC_EXPR_DIR}/Zeng-Aging-Mouse-10Xv3-log2.h5ad" \
        "${SC_EXPR_URL}/Zeng-Aging-Mouse-10Xv3-log2.h5ad" \
        --progress-bar
    echo "  ✓ scRNA-seq log2 完成!"
fi

if [ "${DOWNLOAD_SC_RAW}" = true ]; then
    echo ""
    echo "[2/4] 下载 scRNA-seq raw 表达矩阵..."
    echo "      文件: Zeng-Aging-Mouse-10Xv3-raw.h5ad"
    echo "      估计大小: ~16 GB"
    echo ""
    curl -L --retry 5 --retry-delay 30 \
        -o "${SC_EXPR_DIR}/Zeng-Aging-Mouse-10Xv3-raw.h5ad" \
        "${SC_EXPR_URL}/Zeng-Aging-Mouse-10Xv3-raw.h5ad" \
        --progress-bar
    echo "  ✓ scRNA-seq raw 完成!"
fi

# ============================================================================
# 2. MERFISH 空间转录组 表达矩阵
# ============================================================================
MERFISH_EXPR_DIR="${DATA_DIR}/spatial/MERFISH_638850/expression"
mkdir -p "${MERFISH_EXPR_DIR}"

if [ "${DOWNLOAD_MERFISH_IMPUTED}" = true ]; then
    echo ""
    echo "[3/4] 下载 MERFISH imputed gene 表达矩阵..."
    echo "      文件: C57BL6J-638850-imputed-log2.h5ad"
    echo "      估计大小: ~47 GB (这是最大的单个文件!)"
    echo "      预计时间: 2-6小时 (取决于网速)"
    echo "      ⚠️  建议在稳定的网络环境下进行"
    echo ""
    echo "      基因数: 8,460 (imputed)"
    echo "      细胞数: ~4,000,000"
    echo ""

    MERFISH_IMPUTED_URL="${BASE_URL}/expression_matrices/MERFISH-C57BL6J-638850-imputed/20240831"

    # 使用断点续传模式 (如果支持)
    curl -L --retry 10 --retry-delay 60 \
        -o "${MERFISH_EXPR_DIR}/C57BL6J-638850-imputed-log2.h5ad" \
        "${MERFISH_IMPUTED_URL}/C57BL6J-638850-imputed-log2.h5ad" \
        --progress-bar || {
        echo "  ⚠️ MERFISH imputed 下载失败 (可能网络中断)"
        echo "  重新运行此脚本将从零开始下载 (S3不支持断点续传)"
    }
    echo "  ✓ MERFISH imputed 完成!"
fi

if [ "${DOWNLOAD_MERFISH_SECTIONS}" = true ]; then
    echo ""
    echo "[4/4] 下载 MERFISH 各切片原始数据..."
    echo "      共约15个切片, 每个约0.5-1.5 GB"

    MERFISH_SECTION_URL="${BASE_URL}/expression_matrices/MERFISH-C57BL6J-638850-sections/20230630"

    # 切片 01-15 (注意: 切片07缺失)
    for i in $(seq -w 1 15); do
        [ "$i" = "07" ] && continue  # 切片07不存在

        echo "  -> 下载 切片 ${i}..."
        curl -L --retry 3 --retry-delay 10 \
            -o "${MERFISH_EXPR_DIR}/C57BL6J-638850.${i}-log2.h5ad" \
            "${MERFISH_SECTION_URL}/C57BL6J-638850.${i}-log2.h5ad" \
            --progress-bar 2>/dev/null &
    done
    wait
    echo "  ✓ MERFISH sections 完成!"
fi

# ============================================================================
# 总结
# ============================================================================
echo ""
echo "============================================"
echo " 下载完成! 数据目录大小:"
echo "============================================"
du -sh "${DATA_DIR}" 2>/dev/null || true
echo ""
echo "各子目录:"
for dir in "${DATA_DIR}/single_cell/expression" "${DATA_DIR}/single_cell/metadata" "${DATA_DIR}/spatial/MERFISH_638850/expression" "${DATA_DIR}/spatial/MERFISH_638850/metadata"; do
    if [ -d "${dir}" ]; then
        echo "  $(du -sh ${dir} 2>/dev/null | cut -f1)  ${dir#$DATA_DIR/}"
    fi
done
echo ""
echo "下一步:"
echo "  Rscript scripts/filter_hypothalamus.R  # 提取下丘脑数据"
echo "  或在 R 中手动分析:"
echo "    library(Seurat)"
echo "    library(anndata)"
echo "    sce <- read_h5ad('data/single_cell/expression/Zeng-Aging-Mouse-10Xv3-log2.h5ad')"
