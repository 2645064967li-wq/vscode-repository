"""Prepare sex-balanced ARC/PVH/DMH neuron-microglia CellChat inputs."""
from pathlib import Path
import gzip, shutil, tempfile, json
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
from scipy.io import mmwrite

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"data/cellchat_region"
OUT.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(20260705)
a=ad.read_h5ad(ROOT/"data/hypothalamus_complete/hypothalamus_HY_raw_82431.h5ad",backed="r")
pred=pd.read_csv(ROOT/"results/spatial_region_mapping/microglia_region_predictions.csv",index_col=0)
pred=pred.loc[pred.region_label_final!="Ambiguous"]
regions={"ARC_ARH":"ARH","PVN_PVH":"PVH","DMH":"DMH"}
obs=a.obs.copy(); obs["_idx"]=np.arange(len(obs)); obs["_group"]=""
neur=obs.cell_type_major.isin(["GABAergic_Neuron","Glutamatergic_Neuron"])
for region,pat in regions.items():
    obs.loc[neur & obs.cluster_name.astype(str).str.contains(pat,case=False,regex=False),"_region_direct"]=region
map_region=pred.region_label_final.to_dict()
micro=obs.cell_type_major.eq("Microglia")
obs.loc[micro,"_region_direct"]=obs.loc[micro].index.map(map_region)
obs.loc[micro,"_group"]="Microglia"
obs.loc[obs.cell_type_major.eq("GABAergic_Neuron"),"_group"]="GABAergic_Neuron"
obs.loc[obs.cell_type_major.eq("Glutamatergic_Neuron"),"_group"]="Glutamatergic_Neuron"
summary=[]; selected={}
for region in regions:
    z=obs.loc[obs["_region_direct"].eq(region) & obs["_group"].ne("")].copy()
    for group in ["Microglia","GABAergic_Neuron","Glutamatergic_Neuron"]:
        counts=z.loc[z._group.eq(group)].groupby(["age_group","sex_group"]).size()
        strata=[("Adult(2m)","Female"),("Adult(2m)","Male"),("Aged(18m)","Female"),("Aged(18m)","Male")]
        n=min(int(counts.get(s,0)) for s in strata)
        if n<10: raise RuntimeError(f"{region} {group} minimum stratum {n}")
        for age,sex in strata:
            pool=z.loc[z._group.eq(group)&z.age_group.eq(age)&z.sex_group.eq(sex)]
            take=rng.choice(pool._idx.to_numpy(),size=n,replace=False)
            selected.setdefault((region,age),[]).extend(take.tolist())
        summary.append({"region":region,"group":group,"per_sex_age":n,"cells_per_age":2*n})
pd.DataFrame(summary).to_csv(OUT/"balancing_summary.csv",index=False)

genes=a.var.gene_symbol.astype(str).to_numpy(); uniq,codes=np.unique(genes,return_inverse=True)
G=sparse.coo_matrix((np.ones(len(codes)),(np.arange(len(codes)),codes)),shape=(len(codes),len(uniq))).tocsr()
for (region,age),idx in selected.items():
    sub=a[np.array(idx),:].to_memory(); X=(sub.X.tocsr()@G).T.tocsr()
    meta=sub.obs.copy(); meta["cell_type_major"]=meta.cell_type_major.astype(str)
    meta.loc[meta.cell_type_major.eq("Microglia"),"cell_type_major"]="Microglia"
    meta["analysis_region"]=region
    label="young" if age.startswith("Adult") else "aged"; d=OUT/region/label; d.mkdir(parents=True,exist_ok=True)
    tmp=d/"counts.mtx"; mmwrite(tmp,X)
    with open(tmp,"rb") as fi,gzip.open(d/"counts.mtx.gz","wb") as fo: shutil.copyfileobj(fi,fo)
    tmp.unlink(); pd.Series(uniq,name="gene").to_csv(d/"genes.csv",index=False)
    meta.to_csv(d/"metadata.csv",index=False)
a.file.close()
print(pd.DataFrame(summary).to_string(index=False))