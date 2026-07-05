"""Transfer ARC/PVH/DMH labels from MERFISH microglia to ageing snRNA microglia."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold, LeaveOneGroupOut
from sklearn.pipeline import Pipeline

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/spatial_region_mapping"; OUT.mkdir(parents=True,exist_ok=True)
sp=ad.read_h5ad(ROOT/"results/spatial_region_feasibility/merfish_microglia_ARH_PVH_DMH_raw.h5ad")
sn=ad.read_h5ad(ROOT/"data/hypothalamus_complete/microglia_HY_raw_3057.h5ad")
sp.var_names=sp.var.gene_symbol.astype(str); sp.var_names_make_unique()
sn.var_names=sn.var.gene_symbol.astype(str); sn.var_names_make_unique()
common=sp.var_names.intersection(sn.var_names)
sp=sp[:,common].copy(); sn=sn[:,common].copy()

def lognorm_z(X):
    X=X.tocsr().astype(float) if sparse.issparse(X) else sparse.csr_matrix(X,dtype=float)
    totals=np.asarray(X.sum(1)).ravel(); scale=np.divide(1e4,totals,out=np.zeros_like(totals),where=totals>0)
    X=sparse.diags(scale)@X; X.data=np.log1p(X.data)
    A=X.toarray(); sd=A.std(0); sd[sd==0]=1
    return (A-A.mean(0))/sd
Xs=lognorm_z(sp.X); Xt=lognorm_z(sn.X)
y=sp.obs.analysis_region.astype(str).to_numpy(); groups=sp.obs.source_section.astype(str).to_numpy()
model=lambda: Pipeline([("select",SelectKBest(f_classif,k=min(100,len(common)))),
                         ("clf",LogisticRegression(max_iter=3000,class_weight="balanced",C=0.5))])
metrics={}
for name,cv in [("stratified_random",StratifiedKFold(5,shuffle=True,random_state=20260705)),
                ("leave_one_section_out",LeaveOneGroupOut())]:
    truth=[]; pred=[]
    splits=cv.split(Xs,y,groups) if name=="leave_one_section_out" else cv.split(Xs,y)
    for tr,te in splits:
        m=model(); m.fit(Xs[tr],y[tr]); pred.extend(m.predict(Xs[te])); truth.extend(y[te])
    labels=sorted(np.unique(y))
    metrics[name]={"balanced_accuracy":balanced_accuracy_score(truth,pred),
                   "macro_f1":f1_score(truth,pred,average="macro"),
                   "confusion":confusion_matrix(truth,pred,labels=labels).tolist(),
                   "labels":labels}
m=model(); m.fit(Xs,y); proba=m.predict_proba(Xt); classes=m.named_steps["clf"].classes_
order=np.argsort(proba,axis=1); best=order[:,-1]; second=order[:,-2]
mx=proba[np.arange(len(sn)),best]; margin=mx-proba[np.arange(len(sn)),second]
pred=classes[best]; confident=(mx>=0.60)&(margin>=0.20)
sn.obs["predicted_region"]=pred
sn.obs["region_probability"]=mx
sn.obs["region_margin"]=margin
sn.obs["region_confidence"]=np.where(confident,"high","ambiguous")
sn.obs["region_label_final"]=np.where(confident,pred,"Ambiguous")
for i,c in enumerate(classes): sn.obs[f"prob_{c}"]=proba[:,i]
sn.var.index.name = None
sn.write_h5ad(OUT/"microglia_3057_spatial_region_predictions.h5ad",compression="gzip")
sn.obs.to_csv(OUT/"microglia_region_predictions.csv")
metrics["n_common_genes"]=len(common); metrics["n_high_confidence"]=int(confident.sum())
metrics["prediction_counts"]=sn.obs.region_label_final.value_counts().to_dict()
(OUT/"classifier_metrics.json").write_text(json.dumps(metrics,indent=2),encoding="utf-8")
print(json.dumps(metrics,indent=2))