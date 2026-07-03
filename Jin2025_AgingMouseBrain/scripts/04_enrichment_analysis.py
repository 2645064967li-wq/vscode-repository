"""
============================================================================
Phase 4-6: GO/KEGG/GSEA 富集分析
对每个细胞类型的DEG结果进行通路富集
============================================================================
"""
import pandas as pd
import numpy as np
import os, glob, json, urllib.request, urllib.parse

RESULT_DIR = "d:/vscode/Jin2025_AgingMouseBrain/results"
DEG_DIR = f"{RESULT_DIR}/DEGs"
ENRICH_DIR = f"{RESULT_DIR}/enrichment"
os.makedirs(ENRICH_DIR, exist_ok=True)

# ===========================================================================
# 1. 使用 Enrichr API 做 GO/KEGG 富集 (无需本地安装数据库)
# ===========================================================================
ENRICHR_URL = "https://maayanlab.cloud/Enrichr"

def enrichr_enrich(gene_list, gene_set_libraries=None, max_genes=500):
    """通过Enrichr REST API进行富集分析"""
    if gene_set_libraries is None:
        gene_set_libraries = [
            'GO_Biological_Process_2023',
            'GO_Molecular_Function_2023',
            'GO_Cellular_Component_2023',
            'KEGG_2019_Mouse',
            'Reactome_2022',
            'WikiPathway_2019_Mouse'
        ]

    genes_str = '\n'.join(gene_list[:max_genes])

    # Step 1: Submit gene list
    data = {'list': (None, genes_str), 'description': (None, 'analysis')}
    encoded = urllib.parse.urlencode({'userListId': ''}).encode()

    try:
        req = urllib.request.Request(
            f"{ENRICHR_URL}/addList",
            data=urllib.parse.urlencode({'list': (None, genes_str)}).encode()
        )
        response = urllib.request.urlopen(req)
        result = json.loads(response.read())
        user_list_id = result['userListId']
    except Exception as e:
        print(f"    Enrichr API error: {e}")
        return {}

    # Step 2: Get enrichment for each library
    all_results = {}
    for lib in gene_set_libraries:
        try:
            query = urllib.parse.urlencode({
                'userListId': user_list_id,
                'backgroundType': lib
            })
            req = urllib.request.Request(f"{ENRICHR_URL}/enrich?{query}")
            response = urllib.request.urlopen(req)
            data = json.loads(response.read())

            df = pd.DataFrame(data[lib], columns=[
                'rank', 'term', 'pvalue', 'zscore', 'combined_score',
                'overlapping_genes', 'adjusted_pvalue', 'old_pvalue', 'old_adjusted_pvalue'
            ])
            if len(df) > 0:
                df['-log10_padj'] = -np.log10(df['adjusted_pvalue'].clip(lower=1e-300))
                df['library'] = lib
                all_results[lib] = df
                print(f"    {lib}: {len(df)} terms enriched")
        except Exception as e:
            print(f"    {lib}: error - {e}")

    return all_results


# ===========================================================================
# 2. 对每个细胞类型做富集分析
# ===========================================================================
print("=" * 60)
print(" Phase 4-6: GO/KEGG 富集分析")
print("=" * 60)

# 加载DEG汇总
summary = pd.read_csv(f"{DEG_DIR}/DEG_summary_all.csv")
print(f"\nDEG各组: {len(summary)}")

# 选Top aging-sensitive的细胞类型做富集
# 优先: Microglia, Tanycyte, Astrocyte, Ependymal (文献报道的衰老敏感类型)
# 外加: Glutamatergic_Neuron, GABAergic_Neuron (最大变化量)

priority_ct = ['Microglia', 'Tanycyte', 'Astrocyte', 'Ependymal',
               'Oligodendrocyte', 'OPC']

# 加载DEG文件并做富集
all_enrich_results = {}

for _, row in summary.iterrows():
    ct = row['cell_type']
    sex = row['sex']
    key = f"{ct}_{sex}"
    deg_file = f"{DEG_DIR}/{key}_DEG.csv"

    if not os.path.exists(deg_file):
        continue

    # 只对DEG>=20的组做富集
    n_sig = row['n_sig']
    if n_sig < 20:
        continue

    print(f"\n--- {key} ({n_sig} DEGs) ---")

    deg = pd.read_csv(deg_file)

    # 上调基因 (aged up)
    up_genes = deg[(deg['p_val_adj'] < 0.05) & (deg['avg_log2FC'] > 0.5)]
    up_list = up_genes.sort_values('avg_log2FC', ascending=False)['gene'].tolist()

    # 下调基因 (aged down)
    down_genes = deg[(deg['p_val_adj'] < 0.05) & (deg['avg_log2FC'] < -0.5)]
    down_list = down_genes.sort_values('avg_log2FC')['gene'].tolist()

    print(f"  UP genes (log2FC>0.5): {len(up_list)}")
    print(f"  DOWN genes (log2FC<-0.5): {len(down_list)}")

    # 上调基因富集
    if len(up_list) >= 10:
        print(f"  Running Enrichr UP...")
        up_results = enrichr_enrich(up_list, max_genes=300)
        if up_results:
            for lib, df in up_results.items():
                if len(df) > 0:
                    fname = f"{ENRICH_DIR}/{key}_UP_{lib}.csv"
                    df.to_csv(fname, index=False)
            all_enrich_results[f"{key}_UP"] = up_results

    # 下调基因富集
    if len(down_list) >= 10:
        print(f"  Running Enrichr DOWN...")
        down_results = enrichr_enrich(down_list, max_genes=300)
        if down_results:
            for lib, df in down_results.items():
                if len(df) > 0:
                    fname = f"{ENRICH_DIR}/{key}_DOWN_{lib}.csv"
                    df.to_csv(fname, index=False)
            all_enrich_results[f"{key}_DOWN"] = down_results

    # 避免API限流
    import time; time.sleep(1)

# ===========================================================================
# 3. 汇总Top富集通路
# ===========================================================================
print(f"\n{'='*60}")
print(" 生成富集汇总...")
print("=" * 60)

enrich_summary = []

for key, libs in all_enrich_results.items():
    for lib, df in libs.items():
        if 'GO_Biological_Process' in lib and len(df) > 0:
            top5 = df.nsmallest(5, 'adjusted_pvalue')
            for _, row in top5.iterrows():
                enrich_summary.append({
                    'comparison': key,
                    'term': row['term'],
                    'p.adjust': row['adjusted_pvalue'],
                    'combined_score': row['combined_score'],
                    'genes': row['overlapping_genes']
                })

if enrich_summary:
    es_df = pd.DataFrame(enrich_summary)
    es_df = es_df.sort_values('p.adjust')
    es_df.to_csv(f"{ENRICH_DIR}/enrichment_summary_top5.csv", index=False)
    print(f"\n富集汇总 (Top 5 per comparison):")
    print(es_df[['comparison', 'term', 'p.adjust']].head(30).to_string())

# ===========================================================================
# 4. Local GO enrichment using gene lists from DEGs
# ===========================================================================
print(f"\n{'='*60}")
print(" 跨细胞类型比较: 共有的衰老基因")
print("=" * 60)

# 找出多个细胞类型共有的衰老DEG
shared_genes = {}
for _, row in summary.iterrows():
    ct = row['cell_type']
    sex = row['sex']
    deg_file = f"{DEG_DIR}/{ct}_{sex}_DEG.csv"
    if not os.path.exists(deg_file):
        continue
    deg = pd.read_csv(deg_file)
    sig = deg[deg['p_val_adj'] < 0.05]
    up = set(sig[sig['avg_log2FC'] > 0]['gene'].tolist())
    down = set(sig[sig['avg_log2FC'] < 0]['gene'].tolist())
    shared_genes[ct] = {'up': up, 'down': down}

# 在所有主要细胞类型中都上调或下调的基因
main_cts = ['Microglia', 'Astrocyte', 'Tanycyte', 'Ependymal',
            'Glutamatergic_Neuron', 'GABAergic_Neuron', 'OPC', 'MOL']
main_cts = [c for c in main_cts if c in shared_genes]

if len(main_cts) >= 3:
    # Pan-aging upregulated
    pan_up = shared_genes[main_cts[0]]['up']
    for ct in main_cts[1:]:
        pan_up = pan_up & shared_genes[ct]['up']
    print(f"\n  Pan-aging UP (在所有{len(main_cts)}个类型中上调): {len(pan_up)} genes")
    if pan_up:
        print(f"    {', '.join(list(pan_up)[:20])}")

    # Pan-aging downregulated
    pan_down = shared_genes[main_cts[0]]['down']
    for ct in main_cts[1:]:
        pan_down = pan_down & shared_genes[ct]['down']
    print(f"  Pan-aging DOWN (在所有{len(main_cts)}个类型中下调): {len(pan_down)} genes")
    if pan_down:
        print(f"    {', '.join(list(pan_down)[:20])}")

    # Microglia-specific aging genes
    if 'Microglia' in shared_genes:
        mg_up = shared_genes['Microglia']['up']
        other_up = set()
        for ct in main_cts:
            if ct != 'Microglia':
                other_up = other_up | shared_genes[ct]['up']
        mg_specific = mg_up - other_up
        print(f"\n  Microglia-specific UP: {len(mg_specific)} genes")
        if mg_specific:
            print(f"    {', '.join(list(mg_specific)[:30])}")

print(f"\n{'='*60}")
print(" Phase 4-6 完成!")
print("=" * 60)
