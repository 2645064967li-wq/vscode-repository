"""Lightweight full-transcriptome microglia clustering and module scoring."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results/microglia_region_analysis"; FIG=OUT/"figures"
OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(exist_ok=True)
a=ad.read_h5ad(ROOT/"data/hypothalamus_complete/microglia_HY_raw_3057.h5ad")
p=ad.read_h5ad(ROOT/"results/spatial_region_mapping/microglia_3057_spatial_region_predictions.h5ad")
po=p.obs.reindex(a.obs_names)
for c in ["predicted_region","region_probability","region_margin","region_confidence","region_label_final","prob_ARC_ARH","prob_DMH","prob_PVN_PVH"]: a.obs[c]=po[c].values
a.var_names=a.var.gene_symbol.astype(str); a.var_names_make_unique()
X=a.X.tocsr().astype(np.float32); lib=np.asarray(X.sum(1)).ravel(); Xn=sparse.diags(1e4/np.maximum(lib,1))@X; Xn.data=np.log1p(Xn.data)
mean=np.asarray(Xn.mean(0)).ravel(); mean2=np.asarray(Xn.power(2).mean(0)).ravel(); var=np.maximum(mean2-mean**2,0)
disp=var/np.maximum(mean,1e-3); hvg=np.argsort(disp)[-2000:]
svd=TruncatedSVD(n_components=30,random_state=20260705); Z=svd.fit_transform(Xn[:,hvg])
best=None
for k in range(4,11):
    lab=KMeans(k,random_state=20260705,n_init=20).fit_predict(Z)
    sil=silhouette_score(Z,lab,sample_size=min(2500,len(lab)),random_state=20260705)
    if best is None or sil>best[0]: best=(sil,k,lab)
a.obsm["X_pca"]=Z; a.obs["cluster_kmeans"]=pd.Categorical(best[2].astype(str))

states={
"Homeostatic":["P2ry12","Tmem119","Cx3cr1","Sall1","Hexb","Gpr34","Fcrls","Siglech"],
"DAM_lipid":["Apoe","Trem2","Lpl","Cst7","Ctsb","Ctsd","Tyrobp","Itgax","Spp1","Lgals3"],
"Interferon":["Ifit1","Ifit2","Ifit3","Isg15","Irf7","Stat1","Oasl2","Rsad2"],
"Inflammatory":["Il1b","Tnf","Nfkbia","Ccl2","Ccl3","Ccl4","Ptgs2","Socs3"],
"Phagolysosomal":["C1qa","C1qb","C1qc","Ctss","Lyz2","Fcgr3","Aif1","Cd68","Lamp1"],
"Proliferative":["Mki67","Top2a","Cenpf","Birc5","Tubb5","Stmn1","Pclaf"]}
metab={
"Glycolysis":["Hk1","Hk2","Pfkp","Aldoa","Gapdh","Pgk1","Pkm","Ldha","Eno1"],
"OXPHOS":["Ndufs1","Ndufs2","Ndufa9","Sdha","Uqcrc1","Cox4i1","Cox5a","Atp5f1a"],
"TCA":["Cs","Aco2","Idh3a","Ogdh","Suclg1","Sdhb","Fh","Mdh2"],
"FAO":["Cpt1a","Cpt2","Acadl","Acadm","Hadha","Hadhb","Acox1"],
"Cholesterol_lipid":["Apoe","Abca1","Abcg1","Lpl","Soat1","Npc1","Lipa","Srebf2"],
"MTOR_AMPK":["Mtor","Rptor","Rictor","Prkaa1","Prkaa2","Tsc1","Tsc2","Rheb"],
"ROS_glutathione":["Nfe2l2","Gpx1","Gpx4","Gsr","Gss","Sod1","Sod2","Cat"],
"Complement":["C1qa","C1qb","C1qc","C3","Cfb","Cfh"],
"TGFB":["Tgfb1","Tgfbr1","Tgfbr2","Smad2","Smad3","Smad4"],
"APOE_TREM2":["Apoe","Trem2","Tyrobp","Lpl","Cst7","Ctsd"]}
used={}
for name,genes in {**states,**metab}.items():
    use=[g for g in genes if g in a.var_names]; used[name]=use
    idx=a.var_names.get_indexer(use); raw=np.asarray(Xn[:,idx].mean(1)).ravel()
    a.obs["score_"+name]=(raw-raw.mean())/(raw.std()+1e-8)
sc=[f"score_{x}" for x in states]
cm=a.obs.groupby("cluster_kmeans",observed=True)[sc].mean(); assign=cm.idxmax(1).str.replace("score_","",regex=False)
a.obs["microglia_state"]=a.obs.cluster_kmeans.map(assign).astype(str)
cm.assign(assigned_state=assign).to_csv(OUT/"cluster_state_scores.csv")
a.var["highly_variable_lightweight"]=False; a.var.iloc[hvg,a.var.columns.get_loc("highly_variable_lightweight")]=True
a.uns["clustering"]={"method":"TruncatedSVD+KMeans","k":int(best[1]),"silhouette":float(best[0]),"seed":20260705}
a.var.index.name=None; a.write_h5ad(OUT/"microglia_3057_region_state_scored.h5ad",compression="gzip")
keys=["donor_label","age_group","sex_group","region_label_final","microglia_state"]
a.obs.groupby(keys,observed=True).size().rename("n_cells").reset_index().to_csv(OUT/"state_counts_by_donor_region.csv",index=False)
scores=[c for c in a.obs if c.startswith("score_")]
a.obs.groupby(["donor_label","age_group","sex_group","region_label_final"],observed=True)[scores].mean().reset_index().to_csv(OUT/"module_scores_by_donor_region.csv",index=False)
a.obs.to_csv(OUT/"microglia_cell_metadata_scored.csv"); (OUT/"gene_sets_used.json").write_text(json.dumps(used,indent=2),encoding="utf-8")
for col in ["microglia_state","region_label_final","age_group"]:
    fig,ax=plt.subplots(figsize=(7,5)); cats=pd.Categorical(a.obs[col]); q=ax.scatter(Z[:,0],Z[:,1],c=cats.codes,s=5,cmap="tab20",alpha=.7)
    ax.set(xlabel="PC1",ylabel="PC2",title=col); handles,_=q.legend_elements(); ax.legend(handles,list(cats.categories),fontsize=7,bbox_to_anchor=(1.02,1),loc="upper left"); fig.tight_layout(); fig.savefig(FIG/f"pca_{col}.png",dpi=220); plt.close(fig)
print(json.dumps({"k":best[1],"silhouette":best[0],"states":a.obs.microglia_state.value_counts().to_dict(),"regions":a.obs.region_label_final.value_counts().to_dict()},indent=2))