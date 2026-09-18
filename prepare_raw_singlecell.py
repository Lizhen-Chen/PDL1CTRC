"""Optional raw single-cell preparation. Download the four stated GEO files into raw/GSE243013 first.
The full sparse matrix may require substantial RAM; published patient-level aggregates are bundled.
"""
from pathlib import Path
import ast,numpy as np,pandas as pd
from scipy.io import mmread
W=Path(__file__).resolve().parent;raw=W/'raw/GSE243013';out=W/'remote_sources/revision_20260916_inputs';out.mkdir(parents=True,exist_ok=True)
tree=ast.parse((W/'remote_sources/gse243013_fullcell_module_scoring.py').read_text())
sig=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign))
genes=pd.read_csv(raw/'GSE243013_genes.csv.gz').iloc[:,0].astype(str).tolist()
barcodes=pd.read_csv(raw/'GSE243013_barcodes.csv.gz').iloc[:,0].astype(str).tolist()
m=mmread(str(raw/'GSE243013_NSCLC_immune_scRNA_counts.mtx.gz')).tocsr()
if m.shape==(len(genes),len(barcodes)):m=m.T.tocsr()
assert m.shape==(len(barcodes),len(genes))
lib=np.asarray(m.sum(axis=1)).ravel();gi={g:i for i,g in enumerate(genes)}
s=pd.DataFrame({'cellID':barcodes,'nCount_RNA_stream':lib,'nFeature_RNA_stream':np.diff(m.indptr)})
for name,gs in sig.items():
 ix=[gi[g] for g in gs if g in gi];a=np.asarray(m[:,ix].sum(axis=1)).ravel();s[name+'_logCP10K']=np.log1p(10000*a/np.maximum(lib,1))
meta=pd.read_csv(raw/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz')
d=meta.merge(s,on='cellID',how='outer',validate='1:1',indicator=True);assert d['_merge'].eq('both').all()
cols=[x for x in d if x.endswith('_logCP10K')]
ag={'n_cells':('cellID','size'),**{c:(c,'mean') for c in cols},'mean_library':('nCount_RNA_stream','mean'),'mean_detected':('nFeature_RNA_stream','mean')}
for fields,name in [(['sampleID','pathological_response','major_cell_type','sub_cell_type'],'sc_subtype_all.csv'),(['sampleID','pathological_response','major_cell_type'],'sc_lineage_all.csv')]:d.groupby(fields,dropna=False).agg(**ag).reset_index().to_csv(out/name,index=False)
cols=[c for c in meta if c not in ['cellID','major_cell_type','sub_cell_type'] and meta.groupby('sampleID')[c].nunique(dropna=False).max()<=1]
meta[cols].drop_duplicates('sampleID').to_csv(out/'sc_clinical_metadata.csv',index=False)
