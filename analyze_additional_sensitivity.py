from pathlib import Path
import numpy as np,pandas as pd,warnings
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
W=Path(__file__).resolve().parent;O=W/'analysis'
cl=pd.read_csv(W/'source_review/GSE243013_primary224_analysis_metadata.csv').set_index('sampleID')
major=pd.read_csv(O/'sc_major_scores.csv');major['center']=major.sampleID.map(cl.center_S1);cols=[c for c in major if c.endswith('_logCP10K')]
rows=[]
for lin,z in major[major.known_ICI&major.n_cells.ge(100)].groupby('major_cell_type'):
 for c in cols:
  f=smf.ols(f'{c} ~ nonMPR + age + male + LUSC + chemo + stage3 + C(center)',z).fit(cov_type='HC3');ci=f.conf_int().loc['nonMPR']
  rows.append(dict(group=lin,module=c.replace('_logCP10K',''),n=int(f.nobs),n_nonMPR=int(z.loc[f.model.data.row_labels,'nonMPR'].sum()),estimate=f.params.nonMPR,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=f.pvalues.nonMPR))
pd.DataFrame(rows).assign(q=lambda d:multipletests(d.p,method='fdr_bh')[1]).to_csv(O/'sc_center_sensitivity.csv',index=False)
t=pd.read_csv(O/'sc_tnk_composition_data.csv');t['center']=t.sampleID.map(cl.center_S1)
rr=[]
for c in cols:
 f=smf.ols(f'{c} ~ nonMPR + age + male + LUSC + chemo + stage3 + Treg_fraction + NK_fraction + terminal_fraction + log_library + C(center)',t).fit(cov_type='HC3');ci=f.conf_int().loc['nonMPR']
 rr.append(dict(module=c.replace('_logCP10K',''),n=int(f.nobs),estimate=f.params.nonMPR,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=f.pvalues.nonMPR))
pd.DataFrame(rr).assign(q=lambda d:multipletests(d.p,method='fdr_bh')[1]).to_csv(O/'sc_composition_center_sensitivity.csv',index=False)
# Spatial sensitivity: recorded nuclei and raw negative-probe background in the same model.
d=pd.read_csv(O/'spatial_sample_scores.csv');bg=pd.read_csv(W/'source_review/GSE221733_initial_probe_technical_covariates.csv');d=d.merge(bg[['geo_accession','negative_probe_mean']],on='geo_accession',validate='1:1');d['log_background']=np.log1p(d.negative_probe_mean)
d.to_csv(O/'spatial_sample_scores_technical.csv',index=False);out=[]
for mod in [c for c in d if c.endswith('_relative')]:
 z=d[d.responder.notna()].copy();z['score']=z[mod]
 with warnings.catch_warnings():
  warnings.simplefilter('ignore')
  f=smf.mixedlm('score ~ panck * responder + log_nuclei + log_background',z,groups=z.patient,re_formula='1',vc_formula={'ROI':'0 + C(roi_key)'}).fit(reml=False,method=['lbfgs','powell'],maxiter=1000,disp=False)
 for term in ['panck','responder','panck:responder']:
  ci=f.conf_int().loc[term];out.append(dict(module=mod,term=term,n_AOI=len(z),n_patients=z.patient.nunique(),n_ROI=z.roi_key.nunique(),estimate=f.params[term],ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=f.pvalues[term],converged=bool(f.converged),patient_variance=float(f.cov_re.iloc[0,0]),roi_variance=float(f.vcomp[0])))
pd.DataFrame(out).assign(q=lambda d:d.groupby('term').p.transform(lambda p:multipletests(p,method='fdr_bh')[1])).to_csv(O/'spatial_background_sensitivity.csv',index=False)
print(pd.DataFrame(rr).to_string(index=False));print(pd.DataFrame(out).query('term=="panck:responder"').to_string(index=False))
