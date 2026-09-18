from pathlib import Path
import json,requests
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
W=Path(__file__).resolve().parent;I=W/'inputs';O=W/'analysis'
old=pd.read_csv(next(I.glob('Supplementary_Table_110_*')));h=pd.read_csv(next(I.glob('Supplementary_Table_111_*')));iv=pd.read_csv(next(I.glob('Supplementary_Table_108_*')))
assert h.mr_keep.all() and not h.duplicated(['target','id.exposure','id.outcome','SNP']).any()
rows=[]
for (target,eid,oid),g in h.groupby(['target','id.exposure','id.outcome']):
 bx=g['beta.exposure'].to_numpy();by=g['beta.outcome'].to_numpy();sy=g['se.outcome'].to_numpy();w=1/sy**2
 b=np.sum(bx*by*w)/np.sum(bx**2*w);n=len(g);q=np.sum(w*(by-b*bx)**2);se=np.sqrt(1/np.sum(bx**2*w))*np.sqrt(max(1,q/(n-1))) if n>1 else sy[0]/abs(bx[0]);p=2*stats.norm.sf(abs(b/se))
 r2x=np.sum(g['beta.exposure']**2/(g['beta.exposure']**2+(g['samplesize.exposure']-2)*g['se.exposure']**2))
 r2y=np.sum(g['beta.outcome']**2/(g['beta.outcome']**2+(g['samplesize.outcome']-2)*g['se.outcome']**2))
 rows.append(dict(original_target_label=target,exposure_id=eid,outcome_id=oid,exposure=g.exposure.iloc[0],outcome=g.outcome.iloc[0],nsnp=n,method='Wald ratio' if n==1 else 'IVW multiplicative random effects',beta=b,se=se,ci_low=b-1.96*se,ci_high=b+1.96*se,p=p,Q=q if n>1 else np.nan,Q_df=n-1 if n>1 else np.nan,Q_p=stats.chi2.sf(q,n-1) if n>1 else np.nan,F_min=g.F_stat.min(),r2_exposure_approx=r2x,r2_outcome_approx=r2y,direction_r2_exposure_greater=r2x>r2y,Egger='Not estimable with fewer than 3 SNPs',MR_PRESSO='Not estimable with 1 or 2 SNPs',colocalization='Not performed; no matched molecular QTL full-locus statistics',effect_scale='Outcome trait units per exposure trait unit'))
r=pd.DataFrame(rows);r['q']=multipletests(r.p,method='fdr_bh')[1]
r.to_csv(O/'mr_complete_results.csv',index=False);iv.to_csv(O/'mr_instruments.csv',index=False);h.to_csv(O/'mr_harmonised.csv',index=False)
j=old.merge(r,left_on=['target','id.exposure','id.outcome'],right_on=['original_target_label','exposure_id','outcome_id'],validate='1:1')
su=dict(tests=len(r),nominal=int((r.p<.05).sum()),FDR_significant=int((r.q<.05).sum()),minimum_p=r.p.min(),minimum_q=r.q.min(),max_beta_discrepancy=float(np.max(np.abs(j.b-j.beta))),max_se_discrepancy=float(np.max(np.abs(j.se_x-j.se_y))),max_p_discrepancy=float(np.max(np.abs(j.pval-j.p))),instruments=len(iv),unique_snps=iv.SNP.nunique(),F_min=iv.F_stat.min(),F_max=iv.F_stat.max(),single_snp=int((r.nsnp==1).sum()),two_snp=int((r.nsnp==2).sum()))
(O/'mr_summary.json').write_text(json.dumps(su,indent=2),encoding='utf8');print(json.dumps(su,indent=2))
# Primary gene coordinates are annotations, never substituted for a molecular perturbation.
coords=[]
for gene in sorted(iv.target.unique()):
 try:
  url=f'https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json'
  rr=requests.get(url,timeout=25);rr.raise_for_status();x=rr.json();coords.append(dict(original_target_label=gene,assembly=x.get('assembly_name'),chromosome=x.get('seq_region_name'),gene_start=x.get('start'),gene_end=x.get('end'),ensembl_id=x.get('id'),source=url))
 except Exception as e:coords.append(dict(original_target_label=gene,status=type(e).__name__))
c=pd.DataFrame(coords);c.to_csv(O/'mr_gene_coordinates.csv',index=False)
if 'gene_start' in c:
 a=iv.merge(c,left_on='target',right_on='original_target_label',how='left')
 a['distance_to_gene_bp']=np.maximum(np.maximum(a.gene_start-a['pos.exposure'],a['pos.exposure']-a.gene_end),0)
 a.to_csv(O/'mr_instrument_locus_annotations.csv',index=False)
 print(a[['target','SNP','pos.exposure','gene_start','gene_end','distance_to_gene_bp']].to_string(index=False))
