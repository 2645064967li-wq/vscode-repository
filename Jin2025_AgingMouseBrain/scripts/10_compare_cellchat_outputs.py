from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
ROOT=Path(r'D:/vscode/Jin2025_AgingMouseBrain'); BASE=ROOT/'results/cellchat_complete'; OUT=BASE/'comparisons'; OUT.mkdir(parents=True,exist_ok=True)
GROUPS=['young_female','aged_female','young_male','aged_male']
summary=pd.concat([pd.read_csv(BASE/g/'run_summary.csv') for g in GROUPS],ignore_index=True); summary.to_csv(OUT/'all_group_summary.csv',index=False)
def compare_tables(young,aged,sex,name,keys,meta_cols):
 y=pd.read_csv(BASE/young/name); a=pd.read_csv(BASE/aged/name)
 yy=y[keys+meta_cols+['prob','pval']].rename(columns={'prob':'prob_young','pval':'pval_young'})
 aa=a[keys+['prob','pval']].rename(columns={'prob':'prob_aged','pval':'pval_aged'})
 d=yy.merge(aa,on=keys,how='outer');
 for c in meta_cols:
  if c not in d: d[c]=''
 for c in ['prob_young','prob_aged']: d[c]=d[c].fillna(0.0)
 for c in ['pval_young','pval_aged']: d[c]=d[c].fillna(1.0)
 d['delta_prob']=d.prob_aged-d.prob_young; eps=1e-12; d['log2FC']=np.log2((d.prob_aged+eps)/(d.prob_young+eps)); d['sex']=sex
 d['status']=np.select([(d.prob_young==0)&(d.prob_aged>0),(d.prob_young>0)&(d.prob_aged==0)],['gained','lost'],default='shared')
 return d.sort_values('delta_prob',key=lambda x:x.abs(),ascending=False)
all_lr=[]; all_path=[]
for sex,y,a in [('Female','young_female','aged_female'),('Male','young_male','aged_male')]:
 lr=compare_tables(y,a,sex,'lr_interactions.csv',['source','target','interaction_name'],['ligand','receptor','interaction_name_2','pathway_name','annotation','evidence']); lr.to_csv(OUT/f'lr_age_difference_{sex.lower()}.csv',index=False); all_lr.append(lr)
 pa=compare_tables(y,a,sex,'pathway_interactions.csv',['source','target','pathway_name'],[]); pa.to_csv(OUT/f'pathway_age_difference_{sex.lower()}.csv',index=False); all_path.append(pa)
lrall=pd.concat(all_lr,ignore_index=True); pathall=pd.concat(all_path,ignore_index=True)
tm=lrall[((lrall.source=='Tanycyte')&(lrall.target=='Microglia'))|((lrall.source=='Microglia')&(lrall.target=='Tanycyte'))].copy(); tm.to_csv(OUT/'tanycyte_microglia_all_LR_age_difference.csv',index=False)
tg=tm[tm.interaction_name.str.contains('TGF|GDF|INH',case=False,na=False)|tm.receptor.str.contains('TGF',case=False,na=False)].copy(); tg.to_csv(OUT/'tanycyte_microglia_TGFB_family_age_difference.csv',index=False)
keys=['source','target','interaction_name']; f=all_lr[0][keys+['delta_prob','log2FC','status']].rename(columns=lambda c:c+'_female' if c not in keys else c); m=all_lr[1][keys+['delta_prob','log2FC','status']].rename(columns=lambda c:c+'_male' if c not in keys else c); c=f.merge(m,on=keys,how='outer').fillna({'delta_prob_female':0,'delta_prob_male':0}); c['direction_consistent']=np.sign(c.delta_prob_female)==np.sign(c.delta_prob_male); c.to_csv(OUT/'cross_sex_LR_consistency.csv',index=False)
flow=pathall.groupby(['sex','pathway_name'])[['prob_young','prob_aged']].sum().reset_index(); flow['delta_prob']=flow.prob_aged-flow.prob_young; flow.to_csv(OUT/'pathway_information_flow_age_difference.csv',index=False)
sns.set_theme(style='ticks',context='paper'); colors={'young_female':'#56B4E9','aged_female':'#D55E00','young_male':'#0072B2','aged_male':'#E69F00'}
fig,ax=plt.subplots(1,2,figsize=(7.2,3)); sns.barplot(summary,x='group',y='total_count',hue='group',palette=colors,legend=False,ax=ax[0]); sns.barplot(summary,x='group',y='total_weight',hue='group',palette=colors,legend=False,ax=ax[1]); ax[0].set_ylabel('Significant LR records'); ax[1].set_ylabel('Total communication weight');
for a in ax:a.tick_params(axis='x',rotation=35);sns.despine(ax=a)
fig.tight_layout(); fig.savefig(OUT/'global_network_age_sex_summary.pdf'); fig.savefig(OUT/'global_network_age_sex_summary.png',dpi=300); plt.close(fig)
if len(tg):
 tg['label']=tg.interaction_name+'\n'+tg.source+'→'+tg.target; top=tg.iloc[np.argsort(-tg.delta_prob.abs().to_numpy())[:20]]
 fig,ax=plt.subplots(figsize=(7.2,max(3,0.3*len(top)))); sns.barplot(top,y='label',x='delta_prob',hue='sex',palette={'Female':'#CC79A7','Male':'#0072B2'},ax=ax); ax.axvline(0,color='black',lw=.7); ax.set_xlabel('Aged − young communication probability'); ax.set_ylabel(''); sns.despine(); fig.tight_layout(); fig.savefig(OUT/'tanycyte_microglia_TGFB_age_difference.pdf'); fig.savefig(OUT/'tanycyte_microglia_TGFB_age_difference.png',dpi=300); plt.close(fig)
print(summary.to_string(index=False)); print('\nTGFB-family Tany/MG rows:',len(tg)); print(tg[['sex','source','target','interaction_name','prob_young','prob_aged','delta_prob','status']].to_string(index=False))
