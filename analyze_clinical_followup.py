"""Source-linked, patient-level recurrence and PD-L1 analyses in GSE243013."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
W=Path(__file__).resolve().parent; O=W/'analysis'; rng=np.random.default_rng(20260916)
meta=pd.read_csv(W/'source_review/GSE243013_primary224_analysis_metadata.csv')
score=pd.read_csv(O/'sc_composite_scores.csv')
d=meta.merge(score[['sampleID','score']+[c for c in score if c.startswith('z_')]],on='sampleID',how='left',validate='1:1')
d['male']=d.gender.eq('M').astype(int);d['LUSC']=d.histology.eq('LUSC').astype(int)
d['stage3']=np.where(d.stage_S1.fillna('').str.match('^I'),d.stage_S1.fillna('').str.match('III|IV').astype(int),np.nan)
d['score_sd']=(d.score-d.score.mean())/d.score.std()
d['TPS_ordinal']=d.PDL1_TPS_category.map({'<1%':0,'1-49%':1,'>=50%':2})
if d.TPS_ordinal.notna().sum()==0:
 print('TPS categories',d.PDL1_TPS_category.unique()); raise ValueError('Check explicit TPS coding')
d['negative_TPS_ordinal']=-d.TPS_ordinal
d.to_csv(O/'sc_clinical_followup.csv',index=False)
def cindex(t,e,s):
 t=np.asarray(t);e=np.asarray(e);s=np.asarray(s)
 valid=(t[:,None]<t[None,:])&(e[:,None]==1)
 den=valid.sum()
 return float(((s[:,None]>s[None,:])[valid].sum()+.5*(s[:,None]==s[None,:])[valid].sum())/den) if den else np.nan
rows=[]
for feature in ['score_sd','non_MPR']:
 for label,cov in [('unadjusted',[]),('clinical',['age','male','LUSC','stage3']),('response_adjusted',['age','male','LUSC','stage3','non_MPR'])]:
  if feature=='non_MPR' and label=='response_adjusted':continue
  fields=[feature]+cov
  z=d.dropna(subset=fields+['RFS_months','RFS_status']);z=z[z.RFS_months.gt(0)]
  f=PHReg(z.RFS_months,z[fields].astype(float),status=z.RFS_status,ties='efron').fit()
  ci=f.conf_int()[0]
  rows.append(dict(feature=feature,model=label,n=len(z),events=int(z.RFS_status.sum()),log_hr=f.params[0],hr=np.exp(f.params[0]),ci_low=np.exp(ci[0]),ci_high=np.exp(ci[1]),p=f.pvalues[0],score_equation_max_abs=float(np.max(np.abs(f.model.score(f.params))))))
pd.DataFrame(rows).assign(q=lambda x:multipletests(x.p,method='fdr_bh')[1]).to_csv(O/'sc_rfs_models.csv',index=False)
z=d.dropna(subset=['score_sd','RFS_months','RFS_status']);z=z[z.RFS_months.gt(0)]
cr=[]
for feat in ['score_sd','non_MPR']:
 vals=[]
 for _ in range(2000):
  b=rng.integers(0,len(z),len(z));v=z.iloc[b]; vals.append(cindex(v.RFS_months,v.RFS_status,v[feat]))
 cr.append(dict(feature=feat,n=len(z),events=int(z.RFS_status.sum()),c_index=cindex(z.RFS_months,z.RFS_status,z[feat]),ci_low=np.nanquantile(vals,.025),ci_high=np.nanquantile(vals,.975)))
pd.DataFrame(cr).to_csv(O/'sc_rfs_cindex.csv',index=False)
z=d.dropna(subset=['score','TPS_ordinal','non_MPR']);au=[]
for f in ['score','negative_TPS_ordinal']:
 ids=[np.where(z.non_MPR.to_numpy()==v)[0] for v in [0,1]];vals=[]
 for _ in range(2000):
  b=np.concatenate([rng.choice(i,len(i),True) for i in ids]);v=z.iloc[b];vals.append(roc_auc_score(v.non_MPR,v[f]))
 au.append(dict(feature=f,n=len(z),non_MPR=int(z.non_MPR.sum()),auc=roc_auc_score(z.non_MPR,z[f]),ci_low=np.quantile(vals,.025),ci_high=np.quantile(vals,.975)))
pd.DataFrame(au).to_csv(O/'sc_pdl1_benchmarks.csv',index=False)
z=d.dropna(subset=['score_sd','TPS_ordinal','non_MPR']);f=sm.GLM(z.non_MPR,sm.add_constant(z[['score_sd','TPS_ordinal']]),family=sm.families.Binomial()).fit()
ci=f.conf_int();pd.DataFrame([dict(feature=c,n=len(z),non_MPR=int(z.non_MPR.sum()),OR=np.exp(f.params[c]),ci_low=np.exp(ci.loc[c,0]),ci_high=np.exp(ci.loc[c,1]),p=f.pvalues[c]) for c in ['score_sd','TPS_ordinal']]).to_csv(O/'sc_pdl1_adjusted_models.csv',index=False)
# Retain the same response contrast with original S1 stage instead of GEO stage.
m=pd.read_csv(O/'sc_major_scores.csv');cl=meta.set_index('sampleID');st=pd.Series(np.where(cl.stage_S1.fillna('').str.match('^I'),cl.stage_S1.fillna('').str.match('III|IV').astype(int),np.nan),index=cl.index);m['stage_S1_III_IV']=m.sampleID.map(st);stage=[]
for lin,g in m[m.known_ICI & m.n_cells.ge(100)].groupby('major_cell_type'):
 for col in [x for x in g if x.endswith('_logCP10K')]:
  fields=['nonMPR','age','male','LUSC','chemo','stage_S1_III_IV'];z=g.dropna(subset=fields+[col]);fit=sm.OLS(z[col],sm.add_constant(z[fields].astype(float))).fit(cov_type='HC3');ci=fit.conf_int().loc['nonMPR']
  stage.append(dict(group=lin,module=col.replace('_logCP10K',''),n=len(z),n_nonMPR=int(z.nonMPR.sum()),estimate=fit.params.nonMPR,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=fit.pvalues.nonMPR))
pd.DataFrame(stage).assign(q=lambda x:multipletests(x.p,method='fdr_bh')[1]).to_csv(O/'sc_stage_source_sensitivity.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(cr).to_string(index=False));print(pd.DataFrame(au).to_string(index=False))
