from pathlib import Path
import json,ast,warnings
import numpy as np,pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
W=Path(__file__).resolve().parent;O=W/'analysis';O.mkdir(exist_ok=True)
meta=pd.read_csv(W/'source_review/GSE221733_AOI_clinical_map.csv')
x=pd.read_excel(W/'remote_sources/data/GSE221733/GSE221733_4301_CTA_norm.xlsx',index_col=0)
x.index=x.index.str.strip();meta.sample_title=meta.sample_title.str.strip()
d=meta.set_index('sample_title').loc[x.index].copy();d['patient']=d.patient_id_whitespace_normalized
assert not d.index.duplicated().any()
d['panck']=d.segment.eq('PanCK pos').astype(int);d['responder']=d.response.map({'Responder':1,'Non-responder':0})
d['log_nuclei']=np.log1p(d.aoinucleicount);d['log_area']=np.log1p(d.area)
d['global_median']=x.median(axis=1)
medcenter=x.sub(x.median(axis=1),axis=0)
# Read constant signature definitions without importing or executing the legacy pipeline.
tree=ast.parse((W/'remote_sources/gse243013_fullcell_module_scoring.py').read_text())
sig=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='SIGNATURES' for t in n.targets))
coverage=[];mods={}
for name,gs in sig.items():
 hit=[g for g in gs if g in x];coverage.append(dict(module=name,n_requested=len(gs),n_measured=len(hit),coverage=len(hit)/len(gs),genes=';'.join(hit),missing=';'.join(g for g in gs if g not in x)))
 if len(hit)>=2 and len(hit)/len(gs)>=.5:
  mods[name]=x[hit].mean(axis=1);mods[name+'_relative']=medcenter[hit].mean(axis=1)
d=pd.concat([d,pd.DataFrame(mods)],axis=1)
pd.DataFrame(coverage).to_csv(O/'spatial_coverage.csv',index=False)
d.reset_index().to_csv(O/'spatial_sample_scores.csv',index=False)
x.to_csv(O/'spatial_expression_deposited.csv');medcenter.to_csv(O/'spatial_expression_relative.csv')
# Preserve pairing within each original ROI, then average paired differences within patient.
def paired(expr):
 a=expr.copy();a['patient']=d.patient;a['roi_key']=d.roi_key;a['panck']=d.panck
 wide=a.set_index(['patient','roi_key','panck']).unstack('panck')
 genes=expr.columns
 delta=pd.DataFrame({g:wide[(g,1)]-wide[(g,0)] for g in genes}).dropna(how='all')
 return delta.groupby('patient').mean()
rows=[];deltas={}
for label,expr in [('Deposited',x),('Relative',medcenter)]:
 delta=paired(expr);deltas[label]=delta;delta.to_csv(O/f'spatial_paired_{label.lower()}.csv')
 for g in expr:
  v=delta[g].dropna();n=len(v);se=v.std(ddof=1)/np.sqrt(n);t=stats.ttest_1samp(v,0);ci=stats.t.ppf(.975,n-1)*se
  rows.append(dict(scale=label,gene=g,n_patients=n,estimate=v.mean(),ci_low=v.mean()-ci,ci_high=v.mean()+ci,p=t.pvalue))
gene=pd.DataFrame(rows);gene['q']=gene.groupby('scale').p.transform(lambda p:multipletests(p,method='fdr_bh')[1]);gene.to_csv(O/'spatial_gene_compartments.csv',index=False)
# Nested intercepts: ROI lies within patient. A fixed compartment effect and response interaction
# are estimated separately; continuous covariate adjustment uses recorded nuclei count.
models=[]
for m in mods:
 for adjust in [False,True]:
  z=d[d.responder.notna()].copy();z['score']=z[m]
  formula='score ~ panck * responder'+(' + log_nuclei' if adjust else '')
  try:
   with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    fit=smf.mixedlm(formula,z,groups=z.patient,re_formula='1',vc_formula={'ROI':'0 + C(roi_key)'}).fit(reml=False,method=['lbfgs','powell'],maxiter=1000,disp=False)
   for term in ['panck','responder','panck:responder']:
    ci=fit.conf_int().loc[term]
    models.append(dict(module=m,adjustment='Nuclei adjusted' if adjust else 'Unadjusted',term=term,n_AOI=len(z),n_patients=z.patient.nunique(),n_ROI=z.roi_key.nunique(),estimate=fit.params[term],ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=fit.pvalues[term],converged=bool(fit.converged),patient_variance=float(fit.cov_re.iloc[0,0]),roi_variance=float(fit.vcomp[0]),status='Estimated'))
  except Exception as e:models.append(dict(module=m,adjustment='Nuclei adjusted' if adjust else 'Unadjusted',term='all',status=str(e)))
mm=pd.DataFrame(models);mm['q']=mm.groupby(['adjustment','term','module'],dropna=False).p.transform(lambda p:p)
for (adj,term,rel),ix in mm.assign(relative=mm.module.str.endswith('_relative')).groupby(['adjustment','term','relative']).groups.items():
 p=mm.loc[ix,'p'];valid=p.notna();mm.loc[p[valid].index,'q']=multipletests(p[valid],method='fdr_bh')[1]
mm.to_csv(O/'spatial_mixed_models.csv',index=False)
# Patient means within each compartment avoid AOI pseudoreplication in response comparisons.
r=[]
for scale,expr in [('Deposited',x),('Relative',medcenter)]:
 z=expr.copy();z['patient']=d.patient;z['panck']=d.panck;z['responder']=d.responder
 z=z.groupby(['patient','panck','responder']).mean(numeric_only=True).reset_index()
 for comp,a in z.groupby('panck'):
  for g in expr:
   v=a.loc[a.responder==1,g].dropna();u=a.loc[a.responder==0,g].dropna();t=stats.ttest_ind(v,u,equal_var=False);se=np.sqrt(v.var(ddof=1)/len(v)+u.var(ddof=1)/len(u));df=(v.var(ddof=1)/len(v)+u.var(ddof=1)/len(u))**2/((v.var(ddof=1)/len(v))**2/(len(v)-1)+(u.var(ddof=1)/len(u))**2/(len(u)-1));ci=stats.t.ppf(.975,df)*se
   r.append(dict(scale=scale,panck=comp,gene=g,n_responder=len(v),n_nonresponder=len(u),estimate=v.mean()-u.mean(),ci_low=v.mean()-u.mean()-ci,ci_high=v.mean()-u.mean()+ci,p=t.pvalue))
r=pd.DataFrame(r);r['q']=r.groupby(['scale','panck']).p.transform(lambda p:multipletests(p,method='fdr_bh')[1]);r.to_csv(O/'spatial_gene_response.csv',index=False)
su=dict(AOI=len(d),patients=d.patient.nunique(),ROIs=d.roi_key.nunique(),paired_patients=len(deltas['Relative']),response_AOI=int(d.responder.notna().sum()),response_patients=d[d.responder.notna()].patient.nunique(),panck=d.segment.value_counts().to_dict(),global_median_nuclei_spearman=stats.spearmanr(d.global_median,d.log_nuclei).statistic,compartment_q05=gene.groupby('scale').q.apply(lambda p:int((p<.05).sum())).to_dict(),response_q05=r.groupby(['scale','panck']).q.apply(lambda p:int((p<.05).sum())).to_dict())
su['response_q05']={str(k):v for k,v in su['response_q05'].items()}
(O/'spatial_summary.json').write_text(json.dumps(su,indent=2),encoding='utf8')
print(json.dumps(su,indent=2));print(gene[gene.gene.isin(['LDHB','NT5E','SLC16A1','CTLA4','HK2','EPCAM','CD3D','PTPRC'])].to_string(index=False));print(mm[(mm.term=='panck:responder')&(mm.adjustment=='Nuclei adjusted')].to_string(index=False))
