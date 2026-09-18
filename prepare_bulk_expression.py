"""Normalize deposited bulk data without response-dependent orientation or fitting."""
from pathlib import Path
import gzip,json,io
import pandas as pd,numpy as np
W=Path(__file__).resolve().parent;R=W/'remote_sources/runs/cns_nature_mechanism_20260426_final2/runtime';O=W/'analysis'
def preview(path):
 with gzip.open(path,'rt') as f:
  for i in range(8):print(f.readline()[:200])
for cohort,folder,fn in [('GSE126044','ici_public_validation','GSE126044_counts.txt.gz'),('GSE135222','gap_resolution_public','GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz'),('GSE93157','ici_public_validation','GSE93157_raw_data_values.txt.gz')]:
 path=R/folder/fn
 if cohort=='GSE93157':
  path=R/folder/'GSE93157_series_matrix.txt.gz'
  with gzip.open(path,'rt') as f:lines=f.readlines()
  start=next(i for i,line in enumerate(lines) if line.startswith('!series_matrix_table_begin'))+1
  end=next(i for i,line in enumerate(lines) if line.startswith('!series_matrix_table_end'))
  df=pd.read_csv(io.StringIO(''.join(lines[start:end])),sep='\t',index_col=0)
 else:df=pd.read_csv(path,sep='\t',index_col=0)
 df=df.apply(pd.to_numeric,errors='coerce').dropna(how='all')
 if cohort=='GSE135222':
  map_=pd.read_csv(W/'inputs/hgnc_complete_set.txt',sep='\t',low_memory=False).dropna(subset=['ensembl_gene_id']).drop_duplicates('ensembl_gene_id').set_index('ensembl_gene_id').symbol
  df.index=df.index.str.split('.').str[0].map(map_);df=df[df.index.notna()].groupby(level=0).sum()
  exp=np.log1p(df)
 elif cohort=='GSE126044':
  df=df.groupby(level=0).sum();exp=np.log1p(df.div(df.sum(axis=0),axis=1)*1e6)
 else:
  exp=df.groupby(level=0).mean() # GEO author-normalized, already log2; preserve negative values.
 exp.T.rename_axis('expression_id').to_csv(O/f'{cohort}_expression_log.csv')
 print(cohort,exp.shape,exp.columns.tolist())
