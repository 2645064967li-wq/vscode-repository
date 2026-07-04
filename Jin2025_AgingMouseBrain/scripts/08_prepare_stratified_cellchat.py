"""Prepare balanced sparse inputs for age-by-sex CellChat analysis."""
from pathlib import Path
import gzip, shutil
import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmwrite
ROOT=Path(r"D:/vscode/Jin2025_AgingMouseBrain")
INPUT=ROOT/"data/hypothalamus_complete/hypothalamus_HY_raw_82431.h5ad"
DB_GENES=ROOT/"data/cellchat_complete/cellchatdb_mouse_genes.txt"
OUT=ROOT/"data/cellchat_complete/stratified"
SEED=20260704; MIN_CELLS=30
GROUPS={"young_female":("Adult(2m)","Female"),"aged_female":("Aged(18m)","Female"),"young_male":("Adult(2m)","Male"),"aged_male":("Aged(18m)","Male")}
def main():
 OUT.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(SEED)
 obj=ad.read_h5ad(INPUT,backed="r"); obs=obj.obs.copy()
 counts=pd.crosstab(obs["cell_type_major"],[obs["age_group"],obs["sex_group"]]); expected=pd.MultiIndex.from_tuples(GROUPS.values()); counts=counts.reindex(columns=expected,fill_value=0)
 common=counts.index[(counts>=MIN_CELLS).all(axis=1)].tolist(); targets=counts.loc[common].min(axis=1).astype(int)
 db=[x.strip() for x in DB_GENES.read_text(encoding="utf-8").splitlines() if x.strip()]; genes=[g for g in db if g in obj.var_names]
 pd.Series(sorted(set(db)-set(genes)),name="missing_gene").to_csv(OUT/"missing_db_genes.csv",index=False)
 manifests=[]
 for label,(age,sex) in GROUPS.items():
  selected=[]
  for ct in common:
   candidates=np.flatnonzero(obs["age_group"].eq(age).to_numpy() & obs["sex_group"].eq(sex).to_numpy() & obs["cell_type_major"].eq(ct).to_numpy())
   selected.extend(rng.choice(candidates,size=int(targets[ct]),replace=False).tolist())
  selected=np.asarray(sorted(selected),dtype=int); sub=obj[selected,genes].to_memory(); mat=sub.X.T.tocoo(); gd=OUT/label; gd.mkdir(parents=True,exist_ok=True)
  plain=gd/"counts.mtx"; mmwrite(plain,mat,field="integer")
  with plain.open("rb") as src,gzip.open(gd/"counts.mtx.gz","wb") as dst: shutil.copyfileobj(src,dst)
  plain.unlink(); pd.Series(sub.var_names,name="gene").to_csv(gd/"genes.csv",index=False)
  meta=sub.obs[["cell_type_major","donor_label","age_group","sex_group","library_label"]].copy(); meta.insert(0,"cell_label",sub.obs_names); meta.to_csv(gd/"metadata.csv",index=False)
  manifests.append(meta.groupby("cell_type_major",observed=True).size().rename(label)); print(label,sub.shape,"nnz=",mat.nnz)
 pd.concat(manifests,axis=1).to_csv(OUT/"balanced_cell_counts.csv"); targets.rename("cells_per_group").to_csv(OUT/"sampling_targets.csv")
 pd.DataFrame({"parameter":["seed","min_cells","db_genes_present","cell_types"],"value":[SEED,MIN_CELLS,len(genes),";".join(common)]}).to_csv(OUT/"preparation_manifest.csv",index=False); obj.file.close()
if __name__=="__main__": main()
