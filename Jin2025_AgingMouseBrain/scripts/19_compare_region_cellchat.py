from pathlib import Path
import pandas as pd, numpy as np
ROOT=Path(__file__).resolve().parents[1]; BASE=ROOT/"results/cellchat_region"; OUT=BASE/"comparisons"; OUT.mkdir(exist_ok=True)
regions=["ARC_ARH","PVN_PVH","DMH"]; summaries=[]; all_delta=[]; path_delta=[]
keys=["source","target","interaction_name","interaction_name_2","ligand","receptor","pathway_name","annotation"]
for region in regions:
    y=pd.read_csv(BASE/region/"young/lr_interactions.csv"); o=pd.read_csv(BASE/region/"aged/lr_interactions.csv")
    m=y[keys+["prob","pval"]].merge(o[keys+["prob","pval"]],on=keys,how="outer",suffixes=("_young","_aged")).fillna({"prob_young":0,"prob_aged":0,"pval_young":1,"pval_aged":1})
    m["delta_prob"]=m.prob_aged-m.prob_young; m["region"]=region
    m["microglia_involved"]=m.source.eq("Microglia")|m.target.eq("Microglia")
    m["direction"]=np.where(m.source.eq("Microglia")&m.target.ne("Microglia"),"Microglia_outgoing",np.where(m.target.eq("Microglia")&m.source.ne("Microglia"),"Microglia_incoming",np.where(m.source.eq("Microglia")&m.target.eq("Microglia"),"Microglia_autocrine","Neuron_neuron")))
    m.sort_values("delta_prob",key=abs,ascending=False).to_csv(OUT/f"{region}_age_lr_delta.csv",index=False)
    m[m.microglia_involved].sort_values("delta_prob",key=abs,ascending=False).to_csv(OUT/f"{region}_microglia_age_lr_delta.csv",index=False)
    all_delta.append(m)
    py=pd.read_csv(BASE/region/"young/pathway_interactions.csv"); po=pd.read_csv(BASE/region/"aged/pathway_interactions.csv")
    pk=["source","target","pathway_name"]; q=py[pk+["prob"]].merge(po[pk+["prob"]],on=pk,how="outer",suffixes=("_young","_aged")).fillna(0)
    q["delta_prob"]=q.prob_aged-q.prob_young;q["region"]=region;path_delta.append(q)
    summaries.extend([pd.read_csv(BASE/region/f"{age}/summary.csv") for age in ["young","aged"]])
pd.concat(summaries).to_csv(OUT/"all_run_summaries.csv",index=False)
A=pd.concat(all_delta,ignore_index=True); A.to_csv(OUT/"all_regions_age_lr_delta.csv",index=False)
P=pd.concat(path_delta,ignore_index=True);P.to_csv(OUT/"all_regions_age_pathway_delta.csv",index=False)
focus="TGFB|CX3|CSF|IL34|COMPLEMENT|MIF|SPP1|GRN|PSAP|GAS6|PROS1|AXL|MERTK|CCL|CXCL|APOE|TREM"
A[A.microglia_involved & A[["pathway_name","interaction_name_2","ligand","receptor"]].astype(str).agg(" ".join,axis=1).str.contains(focus,case=False,regex=True)].sort_values(["region","delta_prob"],key=lambda x:abs(x) if x.name=="delta_prob" else x).to_csv(OUT/"microglia_priority_pathways_age_delta.csv",index=False)
print(pd.concat(summaries).to_string(index=False))
print("\nTop microglia age shifts:")
print(A[A.microglia_involved].sort_values("delta_prob",key=abs,ascending=False)[["region","source","target","interaction_name_2","pathway_name","prob_young","prob_aged","delta_prob"]].head(20).to_string(index=False))