"""
Phase 4 修复版: ENSEMBL -> Gene Symbol 转换 + Enrichr富集
"""
import scanpy as sc, pandas as pd, numpy as np, json, urllib.request, urllib.parse, os, time

RESULT_DIR = "d:/vscode/Jin2025_AgingMouseBrain/results"
DEG_DIR = f"{RESULT_DIR}/DEGs"
ENRICH_DIR = f"{RESULT_DIR}/enrichment"
os.makedirs(ENRICH_DIR, exist_ok=True)

# ---- 1. 加载基因映射 (ENSEMBL -> Symbol) ----
print("Loading gene symbol mapping...")
adata = sc.read_h5ad("d:/vscode/Jin2025_AgingMouseBrain/data/hypothalamus/hypo_anatomical_52696.h5ad")
gene_syms = adata.var['gene_symbol'].astype(str).values
ens_to_sym = dict(zip(adata.var_names, gene_syms))
print(f"  {len(ens_to_sym)} gene mappings loaded")

# ---- 2. Enrichr API ----
def enrichr_enrich(gene_symbols, libraries, max_genes=200):
    """用gene symbols做Enrichr富集"""
    genes_str = '\n'.join(gene_symbols[:max_genes])

    # Submit list
    data = urllib.parse.urlencode({'list': genes_str}).encode()
    req = urllib.request.Request('https://maayanlab.cloud/Enrichr/addList', data=data)
    try:
        resp = json.loads(urllib.request.urlopen(req).read())
        user_id = resp['userListId']
    except Exception as e:
        return {}

    results = {}
    for lib in libraries:
        time.sleep(0.5)  # rate limit
        try:
            q = urllib.parse.urlencode({'userListId': user_id, 'backgroundType': lib})
            req = urllib.request.Request(f'https://maayanlab.cloud/Enrichr/enrich?{q}')
            resp = json.loads(urllib.request.urlopen(req).read())
            if lib in resp and len(resp[lib]) > 0:
                df = pd.DataFrame(resp[lib], columns=[
                    'rank','term','pvalue','zscore','combined_score',
                    'genes','adj_pvalue','old_p','old_adj_p'])
                df['library'] = lib
                results[lib] = df
        except Exception as e:
            pass
    return results


# ---- 3. 选择关键细胞类型做富集 ----
print("\nRunning enrichment for key cell types...")

# Focus on the most important ones
key_comparisons = [
    ('Microglia', 'Female', 'Microglia_Female'),
    ('Astrocyte', 'Female', 'Astrocyte_Female'),
    ('Tanycyte', 'Male', 'Tanycyte_Male'),
    ('Tanycyte', 'Female', 'Tanycyte_Female'),
    ('Ependymal', 'Female', 'Ependymal_Female'),
    ('Glutamatergic_Neuron', 'Female', 'Glut_Female'),
    ('GABAergic_Neuron', 'Female', 'GABA_Female'),
]

libraries = ['GO_Biological_Process_2023', 'KEGG_2019_Mouse',
             'Reactome_2022', 'WikiPathway_2019_Mouse']

all_summary = []

for ct, sex, label in key_comparisons:
    deg_file = f"{DEG_DIR}/{ct}_{sex}_DEG.csv"
    if not os.path.exists(deg_file):
        print(f"  {label}: DEG file not found")
        continue

    deg = pd.read_csv(deg_file)

    # Convert ENSEMBL -> Symbol
    deg['symbol'] = deg['gene'].map(ens_to_sym)
    deg = deg.dropna(subset=['symbol'])

    # UP genes (Aged > Adult)
    up = deg[(deg['p_val_adj'] < 0.05) & (deg['avg_log2FC'] > 0.5)]
    up_syms = up.sort_values('avg_log2FC', ascending=False)['symbol'].unique().tolist()

    # DOWN genes (Aged < Adult)
    down = deg[(deg['p_val_adj'] < 0.05) & (deg['avg_log2FC'] < -0.5)]
    down_syms = down.sort_values('avg_log2FC')['symbol'].unique().tolist()

    print(f"\n{'='*50}")
    print(f" {ct} ({sex}): UP={len(up_syms)}, DOWN={len(down_syms)}")
    print(f"{'='*50}")

    for direction, gene_list in [('UP', up_syms), ('DOWN', down_syms)]:
        if len(gene_list) < 10:
            print(f"  {direction}: too few genes ({len(gene_list)}), skip")
            continue

        print(f"  {direction}: {len(gene_list)} genes -> Enrichr...")
        results = enrichr_enrich(gene_list, libraries, max_genes=min(300, len(gene_list)))

        if not results:
            print(f"    No results")
            continue

        for lib, df in results.items():
            if len(df) == 0:
                continue
            fname = f"{ENRICH_DIR}/{ct}_{sex}_{direction}_{lib}.csv"
            df.to_csv(fname, index=False)
            print(f"    {lib}: {len(df)} terms")

            # Collect top terms for summary
            top5 = df.nsmallest(5, 'adj_pvalue')
            for _, row in top5.iterrows():
                all_summary.append({
                    'cell_type': ct, 'sex': sex, 'direction': direction,
                    'library': lib, 'term': row['term'],
                    'p.adjust': row['adj_pvalue'],
                    'combined_score': row['combined_score'],
                    'genes': row['genes']
                })

    time.sleep(1)  # Rate limiting

# ---- 4. 保存汇总 ----
if all_summary:
    summ = pd.DataFrame(all_summary)
    summ = summ.sort_values('p.adjust')
    summ.to_csv(f"{ENRICH_DIR}/enrichment_summary.csv", index=False)

    print(f"\n{'='*60}")
    print(" Top Enriched Pathways Across Cell Types")
    print(f"{'='*60}")
    print(summ[['cell_type','direction','term','p.adjust']].head(25).to_string())

print(f"\nDone! Results in {ENRICH_DIR}/")
