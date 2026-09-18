"""Identical gene-set bulk proxies, with unchanged higher-nonresponse direction.
No bulk score is treated as an exact measurement of the lineage-specific score.
"""
from pathlib import Path
import ast,json
import numpy as np,pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score
import statsmodels.api as sm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests
W=Path(__file__).resolve().parent;O=W/'analysis';rng=np.random.default_rng(20260916)
tree=ast.parse((W/'remote_sources/gse243013_fullcell_module_scoring.py').read_text())
sig=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='SIGNATURES' for t in n.targets))
six=['GLYCOLYSIS','CD8_EXHAUSTION','TREG_SUPPRESSION','ADENOSINE','LACTATE_TRANSPORT','MYELOID_SUPPRESSION'];four=['CD8_EXHAUSTION','TREG_SUPPRESSION','ADENOSINE','MYELOID_SUPPRESSION']
spec=[('GSE126044','GSE126044_verified_analysis_metadata.csv'),('GSE135222','GSE135222_verified_analysis_metadata.csv'),('GSE93157','GSE93157_NSCLC35_verified_analysis_metadata.csv')]
def auc(y,s):
 y=np.asarray(y,dtype=int);s=np.asarray(s,dtype=float);ids=[np.where(y==v)[0] for v in [0,1]];v=[]
 for _ in range(2000):
  b=np.concatenate([rng.choice(i,len(i),True) for i in ids]);v.append(roc_auc_score(y[b],s[b]))
 lo,hi=np.quantile(v,[.025,.975]);degenerate=bool(lo==hi)
 return dict(n=len(y),n_nonresponse=int(y.sum()),auc=roc_auc_score(y,s),ci_low=np.nan if degenerate else lo,ci_high=np.nan if degenerate else hi,ci_status='Not estimable: degenerate empirical bootstrap under complete rank separation' if degenerate else 'Stratified patient-bootstrap percentile interval',p=mannwhitneyu(s[y==1],s[y==0],alternative='two-sided').pvalue)
def cind(t,e,s):
 t=np.asarray(t);e=np.asarray(e);s=np.asarray(s);m=(t[:,None]<t[None,:])&(e[:,None]==1)
 return ((s[:,None]>s[None,:])[m].sum()+.5*(s[:,None]==s[None,:])[m].sum())/m.sum() if m.sum() else np.nan
scores=[];coverage=[];models=[];survival=[];cindices=[];scaling=[]
for cohort,fn in spec:
 meta=pd.read_csv(W/'source_review'/fn);expr=pd.read_csv(O/f'{cohort}_expression_log.csv',index_col=0).loc[meta.expression_id]
 # Exact official symbols plus unambiguous previous-symbol mappings.
 hgnc=pd.read_csv(W/'inputs/hgnc_complete_set.txt',sep='\t',low_memory=False)
 aliases={}
 for _,r in hgnc.dropna(subset=['prev_symbol']).iterrows():
  for old in str(r.prev_symbol).split('|'):aliases.setdefault(old,set()).add(r.symbol)
 for g in set(sum(sig.values(),[])):
  if g not in expr:
   matches=[old for old,target in aliases.items() if target=={g} and old in expr]
   if len(matches)==1:expr[g]=expr[matches[0]]
 expr.to_csv(O/f'{cohort}_expression_analysis.csv')
 d=meta.set_index('expression_id');available=[]
 for name,gs in sig.items():
  measured=[g for g in gs if g in expr and np.isfinite(expr[g]).all() and expr[g].std(ddof=1)>0];ok=len(measured)>=2 and len(measured)/len(gs)>=.5
  coverage.append(dict(cohort=cohort,module=name,requested=len(gs),measured=len(measured),coverage=len(measured)/len(gs),evaluable=ok,genes=';'.join(measured),missing=';'.join(g for g in gs if g not in measured)))
  if ok:
   zg=(expr[measured]-expr[measured].mean())/expr[measured].std(ddof=1);score=zg.mean(axis=1);mu=score.mean();sd=score.std(ddof=1);d[name]=(score-mu)/sd;available.append(name)
   scaling.append(dict(cohort=cohort,module=name,score_mean=mu,score_sd=sd,n_genes=len(measured),formula='Mean gene-wise z of normalized log expression; standardized across cohort'))
  else:d[name]=np.nan
 d['six_module_proxy']=d[six].mean(axis=1,skipna=False);d['four_immune_module_proxy']=d[four].mean(axis=1,skipna=False)
 d['CD274']=((expr.CD274-expr.CD274.mean())/expr.CD274.std()) if 'CD274' in expr else np.nan
 # A favorable immune marker is oriented explicitly toward nonresponse.
 d['negative_CD274']=-d.CD274;d['negative_CYTOTOXICITY']=-d.CYTOTOXICITY
 features=['six_module_proxy','four_immune_module_proxy']+available+['negative_CD274','negative_CYTOTOXICITY']
 if cohort=='GSE126044':endpoints=[('No durable benefit','response_original_binary'),('No objective response','objective_response_binary')]
 elif cohort=='GSE135222':endpoints=[('No durable benefit','durable_clinical_benefit_binary')]
 else:endpoints=[('No objective response','objective_response_binary'),('Progressive disease','disease_control_binary')]
 for endpoint,col in endpoints:
  d[endpoint]=1-d[col]
  for feat in features:
   z=d.dropna(subset=[col,feat]);status='Estimated' if len(z) and z[col].nunique()==2 else 'Not estimable: required module coverage unavailable'
   if status!='Estimated':models.append(dict(cohort=cohort,endpoint=endpoint,feature=feat,status=status));continue
   res=auc(1-z[col],z[feat]);models.append(dict(cohort=cohort,endpoint=endpoint,feature=feat,analysis='All eligible patients',status=status,**res))
   if cohort=='GSE126044' and feat in ['six_module_proxy','four_immune_module_proxy','negative_CD274','negative_CYTOTOXICITY']:
    s=z[z.geo_sample.eq('fresh')];models.append(dict(cohort=cohort,endpoint=endpoint,feature=feat,analysis='Fresh tissue only',status='Estimated',**auc(1-s[col],s[feat])))
 if cohort=='GSE135222':
  for feat in features:
   z=d.dropna(subset=[feat,'pfs_time_days','pfs_event']);s=(z[feat]-z[feat].mean())/z[feat].std();f=PHReg(z.pfs_time_days,s.to_frame(feat),status=z.pfs_event,ties='efron').fit();ci=f.conf_int()[0]
   survival.append(dict(cohort=cohort,feature=feat,n=len(z),events=int(z.pfs_event.sum()),hr=np.exp(f.params[0]),ci_low=np.exp(ci[0]),ci_high=np.exp(ci[1]),p=f.pvalues[0]))
   vals=[]
   for _ in range(2000):
    b=rng.integers(0,len(z),len(z));v=z.iloc[b];vals.append(cind(v.pfs_time_days,v.pfs_event,v[feat]))
   cindices.append(dict(cohort=cohort,feature=feat,n=len(z),events=int(z.pfs_event.sum()),c_index=cind(z.pfs_time_days,z.pfs_event,z[feat]),ci_low=np.nanquantile(vals,.025),ci_high=np.nanquantile(vals,.975)))
 scores.append(d.reset_index())
pd.concat(scores,ignore_index=True).to_csv(O/'bulk_scores.csv',index=False)
pd.DataFrame(coverage).to_csv(O/'bulk_gene_coverage.csv',index=False);pd.DataFrame(scaling).to_csv(O/'bulk_score_scaling.csv',index=False)
r=pd.DataFrame(models);r['q']=r.groupby(['cohort','endpoint','analysis']).p.transform(lambda p:multipletests(p,method='fdr_bh')[1]);r.to_csv(O/'bulk_auc_models.csv',index=False)
r=pd.DataFrame(survival);r['q']=multipletests(r.p,method='fdr_bh')[1];r.to_csv(O/'bulk_survival_models.csv',index=False);pd.DataFrame(cindices).to_csv(O/'bulk_cindex_models.csv',index=False)
print(pd.DataFrame(models).query("feature in ['six_module_proxy','four_immune_module_proxy']").to_string(index=False));print(pd.DataFrame(cindices).head(2).to_string(index=False))
