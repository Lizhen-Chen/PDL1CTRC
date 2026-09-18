"""Reanalyse deposited per-cell normalized signature sums at the biological-sample level.
Primary inference uses one score per sample and lineage; no cell is an independent replicate.
"""
from pathlib import Path
import json, numpy as np,pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.metrics import roc_auc_score,roc_curve
W=Path(__file__).resolve().parent; I=W/'remote_sources/revision_20260916_inputs'; O=W/'analysis';O.mkdir(exist_ok=True)
SEED=20260916
sub=pd.read_csv(I/'sc_subtype_all.csv');cl=pd.read_csv(I/'sc_clinical_metadata.csv')
cl['anti_PD1']=cl['anti-PD1_therapy'].str.strip()
cl['chemotherapy']=cl['chemotherapy'].str.strip()
cl['treatment_group']=np.select([cl.anti_PD1.eq('No'),cl.anti_PD1.eq('unknowm'),cl.anti_PD1.eq('Nivolumab/Placebo'),cl.chemotherapy.eq('No')],['Chemotherapy only','Unknown','Blinded allocation','ICI with chemotherapy undocumented'],default='ICI plus chemotherapy')
cl['documented_ICI']=cl.treatment_group.isin(['ICI with chemotherapy undocumented','ICI plus chemotherapy'])
source=pd.read_csv(W/'source_review/GSE243013_sample_clinical_map.csv')
cl=cl.merge(source[['sampleID','source_patient_id','S1_matched','primary_S1_GEO_concordant_eligible','pathological_response_source_status']],on='sampleID',validate='1:1')
cl['known_ICI']=cl.primary_S1_GEO_concordant_eligible
cl['nonMPR']=cl.pathological_response.map({'non-MPR':1,'MPR':0,'pCR':0})
cl['age']=pd.to_numeric(cl.age,errors='coerce')
cl['male']=cl.gender.map({'M':1,'F':0})
cl['LUSC']=cl.cancer_type.eq('LUSC').astype(int)
cl['chemo']=cl.treatment_group.eq('ICI plus chemotherapy').astype(int)
cl['stage3']=np.where(cl.pre_treatment_staging.fillna('').str.match('^I'),cl.pre_treatment_staging.fillna('').str.match('III|IV').astype(int),np.nan)
cl.to_csv(O/'sc_sample_map.csv',index=False)
cols=[x for x in sub if x.endswith('_logCP10K')]
assert not sub.duplicated(['sampleID','sub_cell_type']).any()
def cellgroup(x):
 if x.startswith('CD4T'):return 'CD4 T'
 if x.startswith('CD8T'):return 'CD8 T'
 if x.startswith('T_gdT'):return 'Gamma delta T'
 if x.startswith('NK_'):return 'NK'
 if x.startswith('ILC'):return 'ILC3'
 if x.startswith('Mφ'):return 'Macrophage'
 if 'DC' in x:return 'Dendritic'
 if x.startswith('Neu_'):return 'Neutrophil'
 if x=='Mast cell':return 'Mast'
 if x=='Plasma_cell':return 'Plasma'
 return 'B cell'
sub['lineage']=sub.sub_cell_type.map(cellgroup)
def aggregate(x,group):
 z=x.copy()
 for c in cols+['mean_library','mean_detected']:z[c]=z[c]*z.n_cells
 z=z.groupby(['sampleID',group],as_index=False)[['n_cells']+cols+['mean_library','mean_detected']].sum()
 for c in cols+['mean_library','mean_detected']:z[c]=z[c]/z.n_cells
 return z.merge(cl,on='sampleID',validate='m:1')
fine=aggregate(sub,'lineage');major=aggregate(sub,'major_cell_type')
fine.to_csv(O/'sc_lineage_scores.csv',index=False);major.to_csv(O/'sc_major_scores.csv',index=False)
sub=sub.merge(cl,on=['sampleID','pathological_response'],validate='m:1')
sub.to_csv(O/'sc_subtype_scores.csv',index=False)
def fit_one(d,col,adjust=True):
 d=d[d.nonMPR.notna() & d[col].notna()].copy()
 base=['nonMPR']+(['age','male','LUSC','chemo','stage3'] if adjust else [])
 base=[c for c in base if d[c].nunique()>1]
 d=d.dropna(subset=base)
 if len(d)<20 or d.nonMPR.value_counts().min()<5:return None
 y=d[col].astype(float)
 X=sm.add_constant(d[base].astype(float),has_constant='add')
 if np.linalg.matrix_rank(X)<X.shape[1]:return None
 f=sm.OLS(y,X).fit(cov_type='HC3');ci=f.conf_int().loc['nonMPR']
 return dict(n=len(d),n_nonMPR=int(d.nonMPR.sum()),n_response=int((1-d.nonMPR).sum()),
    estimate=float(f.params.nonMPR),ci_low=float(ci.iloc[0]),ci_high=float(ci.iloc[1]),p=float(f.pvalues.nonMPR),
    mean_nonMPR=float(y[d.nonMPR==1].mean()),mean_response=float(y[d.nonMPR==0].mean()),
    r_squared=float(f.rsquared),model='HC3 adjusted OLS' if adjust else 'HC3 unadjusted OLS')
def contrasts(frame,group,minimum,label,filter_=None):
 d=frame[frame.known_ICI & frame.nonMPR.notna() & (frame.n_cells>=minimum)].copy()
 if filter_ is not None:d=d[filter_(d)]
 rows=[]
 for name,g in d.groupby(group):
  for col in cols:
   for adj in [False,True]:
    r=fit_one(g,col,adj)
    if r:rows.append(dict(group=name,module=col.replace('_logCP10K',''),analysis=label,**r))
 out=pd.DataFrame(rows)
 if len(out):out['q']=out.groupby('model')['p'].transform(lambda p:multipletests(p,method='fdr_bh')[1])
 return out
res=[]
res.append(contrasts(major,'major_cell_type',100,'primary_major'))
res.append(contrasts(fine,'lineage',100,'primary_resolved'))
res.append(contrasts(sub,'sub_cell_type',30,'primary_subtype'))
res.append(contrasts(major,'major_cell_type',100,'MPR_only',lambda d:d.pathological_response.ne('pCR')))
res.append(contrasts(major,'major_cell_type',100,'combination_only',lambda d:d.chemo.eq(1)))
res.append(contrasts(major,'major_cell_type',200,'minimum_200_cells'))
res.append(contrasts(major,'major_cell_type',30,'minimum_30_cells'))
geosens=major.copy();geosens['known_ICI']=geosens.documented_ICI
res.append(contrasts(geosens,'major_cell_type',100,'GEO_expanded'))
res=pd.concat(res,ignore_index=True);res.to_csv(O/'sc_all_contrasts.csv',index=False)
# Subtype abundance within its original broad lineage, including explicit zero counts.
counts=sub.pivot(index='sampleID',columns='sub_cell_type',values='n_cells').fillna(0)
totals=sub.groupby(['sampleID','major_cell_type']).n_cells.sum().unstack()
ab=[];ab_dat=[]
for name in counts:
 major_name=sub.loc[sub.sub_cell_type.eq(name),'major_cell_type'].iloc[0]
 frac=counts[name]/totals[major_name]
 df=cl.set_index('sampleID').join(frac.rename('fraction')).reset_index()
 df['n_cells']=df.sampleID.map(totals[major_name]);df['fraction_pct']=df.fraction*100
 df=df[df.known_ICI & df.nonMPR.notna() & df.n_cells.ge(100)]
 r=fit_one(df,'fraction_pct')
 if r:ab.append(dict(subtype=name,major_cell_type=major_name,**r))
 ab_dat.append(df[['sampleID','nonMPR','fraction','fraction_pct']].assign(subtype=name,major_cell_type=major_name))
ab=pd.DataFrame(ab);ab['q']=multipletests(ab.p,method='fdr_bh')[1];ab.to_csv(O/'sc_abundance_contrasts.csv',index=False)
pd.concat(ab_dat).to_csv(O/'sc_abundance_data.csv',index=False)
# Composition-adjusted T/NK scores use measured within-lineage proportions.
tnk=major[major.major_cell_type.eq('T/NK cell') & major.known_ICI & major.nonMPR.notna() & major.n_cells.ge(100)].copy()
for name,pat in [('Treg_fraction','Treg'),('NK_fraction','^NK_'),('terminal_fraction','terminal_Tex')]:
 sel=counts.columns.str.contains(pat)
 tnk[name]=tnk.sampleID.map(counts.loc[:,sel].sum(axis=1)/totals['T/NK cell'])
out=[]
for col in cols:
 for label,add in [('clinical',[]),('composition',['Treg_fraction','NK_fraction','terminal_fraction']),('composition_depth',['Treg_fraction','NK_fraction','terminal_fraction','log_library'])]:
  tnk['log_library']=np.log1p(tnk.mean_library)
  fields=['nonMPR','age','male','LUSC','chemo','stage3']+add
  d=tnk.dropna(subset=fields+[col]);f=sm.OLS(d[col],sm.add_constant(d[fields].astype(float))).fit(cov_type='HC3');ci=f.conf_int().loc['nonMPR']
  out.append(dict(module=col.replace('_logCP10K',''),adjustment=label,n=len(d),estimate=f.params.nonMPR,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=f.pvalues.nonMPR))
out=pd.DataFrame(out);out['q']=out.groupby('adjustment').p.transform(lambda p:multipletests(p,method='fdr_bh')[1]);out.to_csv(O/'sc_composition_models.csv',index=False)
tnk.to_csv(O/'sc_tnk_composition_data.csv',index=False)
# Exact six-component discovery score, same patients required for every component.
spec=[('T/NK cell','GLYCOLYSIS'),('T/NK cell','CD8_EXHAUSTION'),('T/NK cell','TREG_SUPPRESSION'),('T/NK cell','ADENOSINE'),('Myeloid cell','LACTATE_TRANSPORT'),('Myeloid cell','MYELOID_SUPPRESSION')]
eligible=major[major.known_ICI & major.nonMPR.notna() & major.n_cells.ge(100)]
d=cl.set_index('sampleID').loc[cl.loc[cl.known_ICI & cl.nonMPR.notna(),'sampleID']].copy(); names=[];scaling=[]
for lin,mod in spec:
 name=lin.split()[0]+'_'+mod;names.append(name)
 d[name]=eligible[eligible.major_cell_type.eq(lin)].set_index('sampleID')[mod+'_logCP10K']
d=d.dropna(subset=names)
for col in names:
 mu=d[col].mean();sd=d[col].std(ddof=1);d['z_'+col]=(d[col]-mu)/sd;scaling.append(dict(component=col,mean=mu,sd=sd,weight=1/6))
d['score']=d[['z_'+c for c in names]].mean(axis=1)
rng=np.random.default_rng(SEED)
def auc_summary(y,s,nboot=2000):
 y=np.asarray(y,int);s=np.asarray(s,float);ids=[np.where(y==v)[0] for v in [0,1]]
 vals=[]
 for _ in range(nboot):
  b=np.concatenate([rng.choice(i,len(i),replace=True) for i in ids]);vals.append(roc_auc_score(y[b],s[b]))
 return dict(n=len(y),events=int(y.sum()),auc=roc_auc_score(y,s),ci_low=np.quantile(vals,.025),ci_high=np.quantile(vals,.975))
au=[]
for col in ['score']+['z_'+c for c in names]:au.append(dict(feature=col,**auc_summary(d.nonMPR,d[col])))
pd.DataFrame(au).to_csv(O/'sc_auc_benchmarks.csv',index=False);d.reset_index().to_csv(O/'sc_composite_scores.csv',index=False)
pd.DataFrame(scaling).to_csv(O/'sc_score_scaling.csv',index=False)
# Paired T/NK lineage difference: independent sample is retained for the interaction contrast.
paired=[]
t=fine[fine.lineage.eq('CD8 T') & fine.n_cells.ge(100)].set_index('sampleID')
n=fine[fine.lineage.eq('NK') & fine.n_cells.ge(100)].set_index('sampleID')
for col in cols:
 z=cl.set_index('sampleID').join((t[col]-n[col]).rename('paired_difference')).reset_index()
 z=z[z.known_ICI & z.nonMPR.notna()];r=fit_one(z,'paired_difference')
 if r:paired.append(dict(module=col.replace('_logCP10K',''),**r))
p=pd.DataFrame(paired);p['q']=multipletests(p.p,method='fdr_bh')[1];p.to_csv(O/'sc_CD8_NK_interaction.csv',index=False)
summary={'all_cells':int(sub.n_cells.sum()),'primary_cells':int(sub[sub.known_ICI].n_cells.sum()),'samples':len(cl),'documented_ICI':int(cl.documented_ICI.sum()),'known_ICI':int(cl.known_ICI.sum()),'primary_samples':int((cl.known_ICI & cl.nonMPR.notna()).sum()),'response_counts':cl[cl.known_ICI].pathological_response.value_counts().to_dict(),'treatment_counts':cl.treatment_group.value_counts().to_dict(),'composite_n':len(d),'auc':au[0]}
(O/'sc_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
print(json.dumps(summary,indent=2));print(res[(res.analysis=='primary_major')&(res.model=='HC3 adjusted OLS')].sort_values('q').head(12).to_string(index=False))
print(out.to_string(index=False))
