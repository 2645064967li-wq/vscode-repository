"""Full-transcriptome clustering, state annotation and module scoring for mapped microglia."""
from pathlib import Path
import json
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/microglia_region_analysis"; FIG=OUT/"figures"
OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(exist_ok=True)
full=ad.read_h5ad(ROOT/"data/hypothalamus_complete/microglia_HY_raw_3057.h5ad")
pred=ad.read_h5ad(ROOT/"results/spatial_region_mapping/microglia_3057_spatial_region_predictions.h5ad")
if not full.obs_names.equals(pred.obs_names):
    pred_obs=pred.obs.reindex(full.obs_names)
else:
    pred_obs=pred.obs
cols=["predicted_region","region_probability","region_margin","region_confidence",
      "region_label_final","prob_ARC_ARH","prob_DMH","prob_PVN_PVH"]
for c in cols: full.obs[c]=pred_obs[c].values
full.layers["counts"]=full.X.copy()
full.var_names=full.var["gene_symbol"].astype(str)
full.var_names_make_unique()
sc.pp.normalize_total(full,target_sum=1e4)
sc.pp.log1p(full)
full.raw=full
sc.pp.highly_variable_genes(full,n_top_genes=3000,flavor="cell_ranger",batch_key="donor_label")
hvg=full[:,full.var.highly_variable].copy()
sc.pp.scale(hvg,max_value=10)
sc.tl.pca(hvg,n_comps=40,svd_solver="randomized",random_state=20260705)
sc.pp.neighbors(hvg,n_neighbors=20,n_pcs=30,random_state=20260705)
sc.tl.umap(hvg,random_state=20260705)
sc.tl.leiden(hvg,resolution=0.5,key_added="leiden_0_5",random_state=20260705)
sc.tl.leiden(hvg,resolution=0.8,key_added="leiden_0_8",random_state=20260705)
full.obsm["X_pca"]=hvg.obsm["X_pca"]; full.obsm["X_umap"]=hvg.obsm["X_umap"]
full.obsp["connectivities"]=hvg.obsp["connectivities"]; full.obsp["distances"]=hvg.obsp["distances"]
full.obs["leiden_0_5"]=hvg.obs["leiden_0_5"].values
full.obs["leiden_0_8"]=hvg.obs["leiden_0_8"].values

state_sets={
"Homeostatic":["P2ry12","Tmem119","Cx3cr1","Sall1","Hexb","Gpr34","Fcrls","Siglech"],
"DAM_lipid":["Apoe","Trem2","Lpl","Cst7","Ctsb","Ctsd","Tyrobp","Itgax","Spp1","Lgals3"],
"Interferon":["Ifit1","Ifit2","Ifit3","Isg15","Irf7","Stat1","Oasl2","Rsad2"],
"Inflammatory":["Il1b","Tnf","Nfkbia","Ccl2","Ccl3","Ccl4","Ptgs2","Socs3"],
"Phagolysosomal":["C1qa","C1qb","C1qc","Ctss","Lyz2","Fcgr3","Aif1","Cd68","Lamp1"],
"Proliferative":["Mki67","Top2a","Cenpf","Birc5","Tubb5","Stmn1","Pclaf"]
}
metabolic_sets={
"Glycolysis":["Hk1","Hk2","Pfkp","Aldoa","Gapdh","Pgk1","Pkm","Ldha","Eno1"],
"OXPHOS":["Ndufs1","Ndufs2","Ndufa9","Sdha","Uqcrc1","Cox4i1","Cox5a","Atp5f1a"],
"TCA":["Cs","Aco2","Idh3a","Ogdh","Suclg1","Sdhb","Fh","Mdh2"],
"FAO":["Cpt1a","Cpt2","Acadl","Acadm","Hadha","Hadhb","Acox1"],
"Cholesterol_lipid":["Apoe","Abca1","Abcg1","Lpl","Soat1","Npc1","Lipa","Srebf2"],
"MTOR_AMPK":["Mtor","Rptor","Rictor","Prkaa1","Prkaa2","Tsc1","Tsc2","Rheb"],
"ROS_glutathione":["Nfe2l2","Gpx1","Gpx4","Gsr","Gss","Sod1","Sod2","Cat"],
"Complement":["C1qa","C1qb","C1qc","C3","Cfb","Cfh"],
"TGFB":["Tgfb1","Tgfbr1","Tgfbr2","Smad2","Smad3","Smad4"],
"APOE_TREM2":["Apoe","Trem2","Tyrobp","Lpl","Cst7","Ctsd"]
}
present={}
for name,genes in {**state_sets,**metabolic_sets}.items():
    use=[g for g in genes if g in full.var_names]
    present[name]=use
    if len(use)>=3: sc.tl.score_genes(full,use,score_name=f"score_{name}",random_state=20260705)

state_cols=[f"score_{x}" for x in state_sets if f"score_{x}" in full.obs]
cluster_means=full.obs.groupby("leiden_0_5",observed=True)[state_cols].mean()
cluster_state=cluster_means.idxmax(axis=1).str.replace("score_","",regex=False)
full.obs["microglia_state"]=full.obs["leiden_0_5"].map(cluster_state).astype(str)
cluster_means.assign(assigned_state=cluster_state).to_csv(OUT/"cluster_state_scores.csv")

full.write_h5ad(OUT/"microglia_3057_region_state_scored.h5ad",compression="gzip")
meta_cols=["donor_label","age_group","sex_group","region_label_final","microglia_state"]
full.obs.groupby(meta_cols,observed=True).size().rename("n_cells").reset_index().to_csv(OUT/"state_counts_by_donor_region.csv",index=False)
score_cols=[c for c in full.obs if c.startswith("score_")]
full.obs.groupby(["donor_label","age_group","sex_group","region_label_final"],observed=True)[score_cols].mean().reset_index().to_csv(OUT/"module_scores_by_donor_region.csv",index=False)
full.obs.to_csv(OUT/"microglia_cell_metadata_scored.csv")
(OUT/"gene_sets_used.json").write_text(json.dumps(present,indent=2),encoding="utf-8")

sc.settings.figdir=FIG; sc.settings.set_figure_params(dpi=180,frameon=False,figsize=(6,5))
for color in ["leiden_0_5","microglia_state","region_label_final","age_group","sex_group"]:
    sc.pl.umap(full,color=color,show=False,save=f"_{color}.png")
print(full.shape)
print(full.obs["microglia_state"].value_counts().to_string())
print(full.obs["region_label_final"].value_counts().to_string())