from pathlib import Path
import anndata as ad, numpy as np, pandas as pd
from scipy import sparse
import matplotlib.pyplot as plt, seaborn as sns
ROOT=Path(r'D:/vscode/Jin2025_AgingMouseBrain'); H5=ROOT/'data/hypothalamus_complete/hypothalamus_HY_raw_82431.h5ad'; OUT=ROOT/'results/cellchat_complete/comparisons'; OUT.mkdir(parents=True,exist_ok=True)
GENES=['Tgfb2','Gdf11','Tgfbr1','Tgfbr2','Acvr1','Acvr2a','Acvr2b']; CTS=['Tanycyte','Microglia']
a=ad.read_h5ad(H5,backed='r'); mask=a.obs.cell_type_major.isin(CTS).to_numpy(); sub=a[mask].to_memory(); a.file.close(); totals=np.asarray(sub.X.sum(axis=1)).ravel(); rows=[]
for gene in GENES:
 x=np.asarray(sub[:,gene].X.toarray()).ravel(); log=np.log1p(x/np.maximum(totals,1)*1e4)
 tmp=sub.obs[['donor_label','age_group','sex_group','cell_type_major']].copy(); tmp['raw']=x; tmp['lognorm']=log; tmp['library']=totals
 for key,z in tmp.groupby(['donor_label','age_group','sex_group','cell_type_major'],observed=True):
  rows.append(dict(zip(['donor_label','age_group','sex_group','cell_type'],key))|{'gene':gene,'n_cells':len(z),'fraction_detected':float((z.raw>0).mean()),'mean_log1p_cp10k':float(z.lognorm.mean()),'pseudobulk_cpm':float(z.raw.sum()/max(z.library.sum(),1)*1e6),'total_gene_counts':int(z.raw.sum())})
d=pd.DataFrame(rows); d['reliable_n20']=d.n_cells>=20; d.to_csv(OUT/'tgfb_donor_expression_validation.csv',index=False)
rel=d[d.reliable_n20].copy(); summ=rel.groupby(['sex_group','cell_type','gene','age_group']).agg(n_donors=('donor_label','nunique'),median_fraction=('fraction_detected','median'),median_expression=('mean_log1p_cp10k','median'),median_pseudobulk_cpm=('pseudobulk_cpm','median')).reset_index(); summ.to_csv(OUT/'tgfb_donor_expression_summary.csv',index=False)
plot=rel[((rel.cell_type=='Tanycyte')&rel.gene.isin(['Tgfb2','Gdf11']))|((rel.cell_type=='Microglia')&rel.gene.isin(['Tgfbr1','Tgfbr2','Acvr1','Acvr2a','Acvr2b']))].copy(); plot['panel']=plot.gene+' | '+plot.cell_type
sns.set_theme(style='ticks',context='paper'); g=sns.catplot(data=plot,x='age_group',y='pseudobulk_cpm',hue='sex_group',col='panel',col_wrap=3,kind='strip',dodge=True,palette={'Female':'#CC79A7','Male':'#0072B2'},height=2.4,aspect=1); g.set_axis_labels('','Pseudobulk CPM'); g.set_titles('{col_name}'); g.figure.tight_layout(); g.figure.savefig(OUT/'tgfb_donor_expression_validation.pdf'); g.figure.savefig(OUT/'tgfb_donor_expression_validation.png',dpi=300); plt.close(g.figure)
print(summ.to_string(index=False))
