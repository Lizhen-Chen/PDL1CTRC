"""Publication figures from the current patient-level result tables.
All points are observed samples. Jitter changes position only, never the data.
Design patterns adapted from inspected academic-figure-skill violin/ROC assets.
"""
from pathlib import Path
import string,json,textwrap
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy import stats
from sklearn.metrics import roc_curve
from statsmodels.duration.survfunc import SurvfuncRight
W=Path(__file__).resolve().parent;A=W/'analysis';O=W/'figures';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Liberation Sans','font.size':7.5,'axes.labelsize':7.5,'axes.titlesize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':6.8,'axes.linewidth':.6,'lines.linewidth':1.1,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
C={'pCR':'#009E73','MPR':'#0072B2','non-MPR':'#D55E00','blue':'#0072B2','orange':'#E69F00','green':'#009E73','red':'#D55E00','purple':'#CC79A7','lightblue':'#56B4E9','black':'#222222'}
CM=LinearSegmentedColormap.from_list('OkabeItoDiverging',[C['blue'],'#FFFFFF',C['red']])
ORD=['pCR','MPR','non-MPR'];rng=np.random.default_rng(20260916)
MOD={'GLYCOLYSIS':'Glycolysis','LACTATE_TRANSPORT':'Lactate transport','OXPHOS':'OXPHOS','FAO':'Fatty acid oxidation','LIPID_UPTAKE':'Lipid uptake','ADENOSINE':'Adenosine','KYNURENINE':'Kynurenine','URIDINE':'Uridine','CD8_EXHAUSTION':'Exhaustion program','CYTOTOXICITY':'Cytotoxicity','TREG_SUPPRESSION':'Treg program','MYELOID_SUPPRESSION':'Myeloid program'}
major=pd.read_csv(A/'sc_major_scores.csv');sub=pd.read_csv(A/'sc_subtype_scores.csv');fine=pd.read_csv(A/'sc_lineage_scores.csv');con=pd.read_csv(A/'sc_all_contrasts.csv');comp=pd.read_csv(A/'sc_composition_models.csv');ab=pd.read_csv(A/'sc_abundance_contrasts.csv');abdata=pd.read_csv(A/'sc_abundance_data.csv');tnk=pd.read_csv(A/'sc_tnk_composition_data.csv');sc=pd.read_csv(A/'sc_composite_scores.csv')
primary=major[major.known_ICI & major.n_cells.ge(100)].copy();subp=sub[sub.known_ICI & sub.n_cells.ge(30)].copy();finep=fine[fine.known_ICI&fine.n_cells.ge(100)].copy()
registry=[]
def wrap(s,n=21):return '\n'.join('\n'.join(textwrap.wrap(part.replace('_',' '),n,break_long_words=False,break_on_hyphens=False)) for part in str(s).split('\n'))
def module(s):return MOD.get(s.replace('_relative',''),s.replace('_',' '))
def setup(n=8,height=None):
 rows=int(np.ceil(n/2));fig,ax=plt.subplots(rows,2,figsize=(183/25.4,(height or rows*58)/25.4),layout='constrained');ax=np.asarray(ax).ravel();fig.set_constrained_layout_pads(w_pad=.035,h_pad=.045,wspace=.12,hspace=.15)
 for a in ax[n:]:a.remove()
 return fig,ax[:n]
def title(ax,letter,s):
 ax.set_title(s,loc='left',pad=7,fontweight='normal');ax.text(-.19,1.08,letter,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
def save(fig,name,caption,panels,sources):
 for ext in ['pdf','svg','png','tiff']:
  kw={'dpi':600 if ext=='tiff' else 300}
  if ext=='tiff':kw['pil_kwargs']={'compression':'tiff_lzw'}
  fig.savefig(O/f'{name}.{ext}',**kw)
 # Separate compact preview is removed after visual QA.
 fig.savefig(O/f'{name}_preview.png',dpi=140)
 (O/f'{name}_legend.md').write_text(caption+'\n',encoding='utf8')
 registry.append({'figure':name,'panels':panels,'sources':sources,'caption':caption})
 plt.close(fig);print(name,flush=True)
def violin(ax,df,y,label=None,groups=ORD,x='pathological_response',colors=None):
 vals=[df.loc[df[x].eq(g),y].dropna().to_numpy() for g in groups];colors=colors or [C[g] for g in groups]
 for i,(v,col) in enumerate(zip(vals,colors)):
  if len(v)>2 and np.std(v)>0:
   z=ax.violinplot(v,positions=[i],widths=.8,showextrema=False)
   for b in z['bodies']:b.set_facecolor(col);b.set_alpha(.16);b.set_edgecolor(col);b.set_linewidth(.5)
  ax.scatter(np.full(len(v),i)+rng.uniform(-.21,.21,len(v)),v,s=5,alpha=.55,c=col,linewidths=0,rasterized=True)
  if len(v):
   q=np.quantile(v,[.25,.5,.75]);ax.plot([i,i],[q[0],q[2]],color='black',lw=1.5);ax.plot([i-.13,i+.13],[q[1],q[1]],color='black',lw=1.5)
 ax.set_xticks(range(len(groups)),[f'{g}\n(n={len(v)})' for g,v in zip(groups,vals)]);ax.set_xlim(-.55,len(groups)-.45);ax.set_ylabel(label or 'Mean log-normalized score')
def forest(ax,d,labels=None,zero=0,color=None,xlabel='Adjusted difference (95% CI)',estimate='estimate',low='ci_low',high='ci_high'):
 d=d.reset_index(drop=True);y=np.arange(len(d));colors=[color or (C['red'] if v>zero else C['blue']) for v in d[estimate]]
 for i,row in d.iterrows():ax.errorbar(row[estimate],i,xerr=[[row[estimate]-row[low]],[row[high]-row[estimate]]],fmt='o',ms=3,color=colors[i],elinewidth=.85,capsize=1.5)
 ax.set_yticks(y,[wrap(x,23) for x in (labels if labels is not None else d.index)]);ax.invert_yaxis();ax.axvline(zero,color='black',ls=':',lw=.6);ax.set_xlabel(xlabel);ax.set_ylim(len(d)-.4,-.6)
 if estimate=='auc':auc_axis(ax)
def auc_axis(ax):
 # Keep full markers and CI end caps visible at valid AUC boundaries.
 ax.set_xlim(-.06,1.06);ax.set_xticks(np.linspace(0,1,6))
def heat(ax,mat,xlabels=None,ylabels=None,vlim=None,label='Adjusted difference',signif=None,colorbar=True,nonnegative=False):
 mat=np.asarray(mat,dtype=float);vlim=vlim or max(.01,np.nanmax(np.abs(mat)));im=ax.imshow(mat,aspect='auto',cmap=LinearSegmentedColormap.from_list('OkabeSequential',['#FFFFFF',C['red']]) if nonnegative else CM,vmin=0 if nonnegative else -vlim,vmax=vlim,interpolation='nearest')
 ax.set_xticks(range(mat.shape[1]),xlabels if xlabels is not None else [],rotation=0);ax.set_yticks(range(mat.shape[0]),ylabels if ylabels is not None else []);ax.tick_params(length=0)
 if signif is not None:
  for i,j in zip(*np.where(np.asarray(signif)<.05)):ax.plot(j,i,'o',ms=1.6,color='black')
 if colorbar:
  cb=ax.figure.colorbar(im,ax=ax,orientation='horizontal',fraction=.08,pad=.12,aspect=25);cb.set_label(label,size=7);cb.ax.tick_params(labelsize=6.5,length=2)
 return im
def scatter(ax,x,y,c=None,xlabel='',ylabel='',identity=False):
 ax.scatter(x,y,s=9,alpha=.6,c=c or C['blue'],linewidths=0,rasterized=True);ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
 if identity:
  lo=min(np.nanmin(x),np.nanmin(y));hi=max(np.nanmax(x),np.nanmax(y));ax.plot([lo,hi],[lo,hi],color='black',lw=.7,ls=':')
def grouped_effect(ax,table,modules,groups,col='adjustment',group_labels=None):
 colors=[C['blue'],C['orange'],C['green'],C['purple'],C['red'],C['black']]
 for j,g in enumerate(groups):
  z=table[table[col].eq(g)].set_index('module').reindex(modules);ys=np.arange(len(modules))+(j-(len(groups)-1)/2)*.15
  ax.errorbar(z.estimate,ys,xerr=np.vstack([z.estimate-z.ci_low,z.ci_high-z.estimate]),fmt='o',ms=2.7,color=colors[j],lw=.75,label=(group_labels or groups)[j],capsize=1)
 ax.set_yticks(range(len(modules)),[wrap(module(m),20) for m in modules]);ax.invert_yaxis();ax.axvline(0,color='black',ls=':',lw=.6);ax.set_xlabel('Adjusted difference (95% CI)');ax.legend(loc='best',frameon=False,handlelength=1)
def paired_spatial(ax,gene,relative=True):
 m=pd.read_csv(A/'spatial_sample_scores.csv').rename(columns={'index':'sample_title'});e=pd.read_csv(A/('spatial_expression_relative.csv' if relative else 'spatial_expression_deposited.csv'),index_col=0);e.index=e.index.str.strip();z=m[['sample_title','patient','roi_key','panck']].copy();z['v']=z.sample_title.str.strip().map(e[gene]);p=z.pivot(index=['patient','roi_key'],columns='panck',values='v').dropna().groupby('patient').mean()
 for _,r in p.iterrows():ax.plot([0,1],r[[0,1]],color=C['black'],alpha=.24,lw=.6)
 for k,col in [(0,C['blue']),(1,C['orange'])]:ax.scatter(np.full(len(p),k),p[k],s=9,color=col,alpha=.75,linewidths=0)
 ax.set_xticks([0,1],['PanCK−','PanCK+']);ax.set_ylabel('Relative log2 expression' if relative else 'Deposited expression');ax.set_xlim(-.25,1.25)
def make_main123():
 fig,axes=setup();names=['Patient-level transcription','T/NK glycolysis','T/NK adenosine','T/NK Treg-associated program','Myeloid lactate transport','Immune subset abundance','Treg abundance within T/NK','Composition-adjusted transcription']
 for i,(a,n) in enumerate(zip(axes,names)):title(a,string.ascii_uppercase[i],n)
 cols=[x for x in sc if x.startswith('z_')];z=sc.assign(order=sc.pathological_response.map(dict(zip(ORD,range(3))))).sort_values(['order','score'])
 heat(axes[0],z[cols].T,ylabels=['T/NK glycolysis','T/NK exhaustion','T/NK Treg program','T/NK adenosine','Myeloid lactate','Myeloid program'],vlim=2.5,label='Component z score')
 axes[0].set_xticks([np.mean(np.where(z.pathological_response.eq(g))[0]) for g in ORD],ORD);axes[0].set_xlabel('212 patients, ordered within response group')
 for a,lin,mod in [(axes[1],'T/NK cell','GLYCOLYSIS'),(axes[2],'T/NK cell','ADENOSINE'),(axes[3],'T/NK cell','TREG_SUPPRESSION'),(axes[4],'Myeloid cell','LACTATE_TRANSPORT')]:violin(a,primary[primary.major_cell_type.eq(lin)],mod+'_logCP10K')
 wanted=['CD4T_Treg_FOXP3','CD4T_Treg_CCR8','CD4T_Treg_MKI67','CD8T_terminal_Tex_LAYN','NK_CD16hi_FGFBP2','CD8T_Tm_IL7R']
 z=ab.set_index('subtype').loc[wanted].reset_index();forest(axes[5],z,['FOXP3 Treg','CCR8 Treg','MKI67 Treg','LAYN terminal CD8','FGFBP2 NK','IL7R memory CD8'],xlabel='Non-MPR difference (percentage points)')
 t=tnk.copy();t['Treg_pct']=t.Treg_fraction*100;violin(axes[6],t,'Treg_pct','Treg cells (% of T/NK)')
 grouped_effect(axes[7],comp,['GLYCOLYSIS','CD8_EXHAUSTION','TREG_SUPPRESSION','ADENOSINE'],['clinical','composition_depth'],group_labels=['Clinical','Clinical + composition + depth'])
 cap='**Figure 1. Immune composition and metabolic transcription distinguish incomplete pathological response.** (A) Six standardized component scores in 212 patients, ordered by pathological response and composite score within each group. Colors saturate at ±2.5 SD for display; all values enter the analyses. (B–D) T/NK glycolysis, adenosine and Treg-associated expression scores in 221 patients. (E) Myeloid lactate-transport scores in 213 patients. (F) Adjusted differences in selected subset proportions within their broad immune lineage, in percentage points. (G) Combined Treg proportion within T/NK cells. (H) Response contrasts before and after adjustment for Treg, NK and terminally exhausted CD8 proportions and mean library size. Points in B–E and G are individual patients; violins show distributions, thick lines show medians and interquartile ranges. Blue, green and vermilion denote MPR, pCR and non-MPR. Error bars in F and H are 95% confidence intervals from patient-level ordinary least-squares models with HC3 standard errors, adjusted for age, sex, histology, stage and chemotherapy-documentation status. MPR, major pathological response; pCR, pathological complete response; Treg, regulatory T cell. Complete estimates and correction families are in Supplementary Tables S3–S5.'
 save(fig,'Figure_1',cap,8,['sc_composite_scores.csv','sc_major_scores.csv','sc_abundance_contrasts.csv','sc_composition_models.csv'])
 fig,axes=setup();names=['Lineage-specific response contrasts','Adenosine in terminal CD8 cells','Treg-associated expression in Tregs','NK exhaustion-associated expression','Paired CD8–NK contrasts','CD8 cell states and adenosine','Lactate transport in myeloid subsets','Glycolysis in B-cell subsets']
 for i,(a,n) in enumerate(zip(axes,names)):title(a,string.ascii_uppercase[i],n)
 z=con.query("analysis=='primary_resolved' and model=='HC3 adjusted OLS'");ms=['GLYCOLYSIS','ADENOSINE','CD8_EXHAUSTION','TREG_SUPPRESSION','CYTOTOXICITY','LACTATE_TRANSPORT'];gs=['CD4 T','CD8 T','NK','Macrophage','Dendritic','B cell','Plasma']
 mat=z.pivot(index='group',columns='module',values='estimate').reindex(index=gs,columns=ms);q=z.pivot(index='group',columns='module',values='q').reindex(index=gs,columns=ms)
 heat(axes[0],mat,['Glycolysis','Adenosine','Exhaustion program','Treg program','Cytotoxicity','Lactate transport'],gs,label='Adjusted non-MPR difference',signif=q);plt.setp(axes[0].get_xticklabels(),rotation=45,ha='right',rotation_mode='anchor');axes[0].tick_params(axis='x',labelsize=6.3)
 violin(axes[1],subp[subp.sub_cell_type.eq('CD8T_terminal_Tex_LAYN')],'ADENOSINE_logCP10K')
 violin(axes[2],subp[subp.sub_cell_type.eq('CD4T_Treg_FOXP3')],'TREG_SUPPRESSION_logCP10K')
 violin(axes[3],subp[subp.sub_cell_type.eq('NK_CD16hi_FGFBP2')],'CD8_EXHAUSTION_logCP10K')
 z=pd.read_csv(A/'sc_CD8_NK_interaction.csv').set_index('module').loc[ms].reset_index();forest(axes[4],z,[module(m) for m in ms],xlabel='Non-MPR effect on CD8 minus NK score')
 x=subp[subp.sub_cell_type.eq('CD8T_Tem_GZMK+GZMH+')].set_index('sampleID');y=subp[subp.sub_cell_type.eq('CD8T_terminal_Tex_LAYN')].set_index('sampleID');p=x[['ADENOSINE_logCP10K','pathological_response']].join(y.ADENOSINE_logCP10K.rename('terminal')).dropna()
 for g in ORD:
  v=p[p.pathological_response.eq(g)];axes[5].scatter(v.ADENOSINE_logCP10K,v.terminal,s=8,c=C[g],alpha=.65,linewidths=0,label=g)
 axes[5].set_xlabel('GZMK/GZMH effector-memory score');axes[5].set_ylabel('LAYN terminal CD8 score');axes[5].legend(frameon=False,loc='upper right',handletextpad=.2)
 z=con.query("analysis=='primary_subtype' and model=='HC3 adjusted OLS'")
 m=z[z.module.eq('LACTATE_TRANSPORT')&z.group.str.startswith('Mφ_')].sort_values('group');forest(axes[6],m,m.group.str.replace('Mφ_','Macrophage ',regex=False),xlabel='Adjusted difference (95% CI)')
 b=z[z.module.eq('GLYCOLYSIS')&z.group.str.match('B[nm]_')].sort_values('group');forest(axes[7],b,b.group.str.replace('_',' '),xlabel='Adjusted difference (95% CI)')
 cap='**Figure 2. Adenosine-associated transcription is concentrated in specific lymphocyte states.** (A) Adjusted response contrasts across resolved lineages and six selected programs; dots denote Benjamini–Hochberg q<0.05 across all 108 tested lineage–program combinations. (B) Adenosine-associated scores in LAYN-positive terminally exhausted CD8 cells. (C) Treg-associated scores within FOXP3-positive Tregs. (D) The exhaustion-associated gene program in FGFBP2-positive NK cells; this label describes shared inhibitory transcripts and does not reclassify NK cells as CD8 T cells. (E) Within-patient CD8-minus-NK score differences regressed on response; error bars show 95% confidence intervals. (F) Paired adenosine scores in effector-memory and terminal CD8 states; each point is a patient with at least 30 cells in each state. (G,H) Complete evaluated macrophage lactate-transport and naive/memory B-cell glycolysis contrasts. Models use patients as the unit of inference, HC3 standard errors and the clinical covariates specified in Figure 1. Lineage analyses require at least 100 cells per patient, subtype analyses at least 30; panel-specific numbers are provided in Supplementary Tables S3 and S6. In B–D and F, green, blue and vermilion denote pCR, MPR and non-MPR; black distribution summaries show the median and interquartile range. In E, G and H, vermilion and blue denote positive and negative coefficients, respectively. Distribution-panel n includes all eligible measured patients; adjusted-model n additionally excludes missing clinical covariates. Original annotation labels are retained without inferring a progenitor-exhausted state from memory-cell labels.'
 save(fig,'Figure_2',cap,8,['sc_all_contrasts.csv','sc_subtype_scores.csv','sc_CD8_NK_interaction.csv'])
 fig,axes=setup();names=['Paired tumour–stromal expression','EPCAM compartment localization','SLC16A1 compartment localization','CTLA4 compartment localization','Panel-wide compartment contrast','Selected compartment-associated genes','Compartment-associated gene programs','Compartment-by-response interactions']
 for i,(a,n) in enumerate(zip(axes,names)):title(a,string.ascii_uppercase[i],n)
 genes=['EPCAM','HK2','LDHB','SLC16A1','NT5E','CTLA4','CD3D','PTPRC'];p=pd.read_csv(A/'spatial_paired_relative.csv',index_col=0);heat(axes[0],p[genes].T,[],genes,vlim=1.4,label='PanCK+ minus PanCK− (relative log2)');axes[0].set_xlabel('21 patients, paired within ROI then averaged')
 for a,g in zip(axes[1:4],['EPCAM','SLC16A1','CTLA4']):paired_spatial(a,g)
 gene=pd.read_csv(A/'spatial_gene_compartments.csv');g=gene[gene.scale.eq('Relative')];axes[4].scatter(g.estimate,-np.log10(g.q.clip(lower=1e-16)),s=4,c=np.where(g.estimate>0,C['orange'],C['blue']),alpha=.5,linewidths=0,rasterized=True)
 for name in ['EPCAM','PTPRC','HK2','LDHB','CTLA4']:
  r=g[g.gene.eq(name)].iloc[0];axes[4].annotate(name,(r.estimate,-np.log10(max(r.q,1e-16))),xytext=(3,4),textcoords='offset points',fontsize=6.5)
 axes[4].set_xlabel('Paired relative log2 difference');axes[4].set_ylabel('−log10 q (1,812 genes)')
 t=g.set_index('gene').loc[genes].reset_index();forest(axes[5],t,t.gene,xlabel='Paired difference (95% CI)')
 mm=pd.read_csv(A/'spatial_mixed_models.csv');m=mm[(mm.adjustment=='Nuclei adjusted')&mm.module.str.endswith('_relative')]
 for a,term in [(axes[6],'panck'),(axes[7],'panck:responder')]:
  z=m[m.term.eq(term)];forest(a,z,[module(x) for x in z.module],xlabel='Relative log2 effect (95% CI)')
 cap='**Figure 3. Spatial localization separates metabolic transcripts from immune programs.** (A) Relative expression differences between paired PanCK-positive and PanCK-negative regions, first paired within each region of interest (ROI) and then averaged equally within 21 patients. (B–D) Corresponding patient-level paired expression for EPCAM, SLC16A1 and CTLA4; connecting lines link the same patient. (E) All 1,812 gene-level paired compartment contrasts, with false-discovery rates calculated over the entire tested panel. (F) Selected gene differences and 95% confidence intervals; NT5E is retained alongside the other prespecified illustrative genes. (G,H) Compartment effects among nonresponders and compartment-by-response interactions for eight evaluable programs from nested mixed-effects models, with patient and ROI random intercepts and fixed adjustment for nuclei count. These models include 67 areas of illumination from 39 patients and 44 ROIs with known response. Relative expression subtracts the within-area median of all 1,812 genes from each deposited log2 value; the deposited-scale and technical-background sensitivities are presented in Supplementary Figure S8. In B–E, orange and blue mark PanCK-positive and PanCK-negative expression or enrichment; in F–H, vermilion and blue denote positive and negative coefficients, respectively; PanCK denotes pan-cytokeratin and does not imply cell purity. Supplementary Tables S12–S15 report coverage, all tested genes, model sizes and correction families. OXPHOS, oxidative phosphorylation.'
 save(fig,'Figure_3',cap,8,['spatial_paired_relative.csv','spatial_gene_compartments.csv','spatial_mixed_models.csv'])

def km(ax,d,time,event,score,cut=None):
 d=d.dropna(subset=[time,event,score]).copy();cut=float(d[score].median()) if cut is None else cut
 for label,sel,col in [('Below median',d[score]<cut,C['blue']),('At or above median',d[score]>=cut,C['red'])]:
  v=d[sel];f=SurvfuncRight(v[time],v[event]);ax.step(np.r_[0,f.surv_times,v[time].max()],np.r_[1,f.surv_prob,f.surv_prob[-1]],where='post',color=col,label=f'{label} (n={len(v)})')
  ct=v.loc[v[event].eq(0),time];ys=np.array([1 if t<f.surv_times[0] else f.surv_prob[np.searchsorted(f.surv_times,t,side='right')-1] for t in ct]);ax.scatter(ct,ys,marker='|',color=col,s=9,linewidths=.6)
 ax.set_ylim(-.03,1.04);ax.set_xlabel('Time (months)' if 'month' in time.lower() else 'Time (days)');ax.set_ylabel('Survival probability');ax.legend(frameon=False,loc='lower left',fontsize=6.4)
def roc(ax,d,y,features,labels):
 for feat,lab,col in zip(features,labels,[C['black'],C['blue'],C['orange'],C['green']]):
  v=d.dropna(subset=[feat,y]);f,t,_=roc_curve(v[y],v[feat]);ax.plot(f,t,label=lab,c=col)
 ax.plot([0,1],[0,1],color='black',ls=':',lw=.6);ax.set(xlim=(0,1),ylim=(0,1.02),xlabel='False-positive rate',ylabel='True-positive rate');ax.legend(frameon=False,loc='lower right',fontsize=6.3)
def make_main4():
 fig,axes=setup();names=['Composite score and pathology','Apparent response discrimination','PD-L1 comparison in matched patients','Recurrence associations','Recurrence by composite score','Bulk score orientation','External response discrimination','Progression in GSE135222']
 for i,(a,n) in enumerate(zip(axes,names)):title(a,string.ascii_uppercase[i],n)
 violin(axes[0],sc,'score','Six-component mean z score');roc(axes[1],sc,'nonMPR',['score','z_T/NK_TREG_SUPPRESSION','z_T/NK_ADENOSINE'],['Six-component score','Treg-associated','Adenosine'])
 p=pd.read_csv(A/'sc_pdl1_benchmarks.csv');forest(axes[2],p,['Composite score','PD-L1 category'],zero=.5,estimate='auc',xlabel='AUC (95% CI), n=101',color=C['blue']);auc_axis(axes[2])
 p=pd.read_csv(A/'sc_rfs_models.csv');p=p[p.feature.eq('score_sd')];forest(axes[3],p,['Unadjusted','Clinical covariates','Clinical + pathology'],zero=1,estimate='hr',xlabel='Hazard ratio per SD (95% CI)',color=C['blue']);axes[3].set_xscale('log')
 f=pd.read_csv(A/'sc_clinical_followup.csv');km(axes[4],f,'RFS_months','RFS_status','score',cut=sc.score.median())
 b=pd.read_csv(A/'bulk_scores.csv');z=b[b.cohort.eq('GSE126044')].copy();z['group']=z.response_original_binary.map({0:'No benefit',1:'Durable benefit'});violin(axes[5],z,'six_module_proxy','Six-module bulk proxy',groups=['No benefit','Durable benefit'],x='group',colors=[C['red'],C['blue']])
 p=pd.read_csv(A/'bulk_auc_models.csv');p=p[p.analysis.eq('All eligible patients')&(((p.cohort.isin(['GSE126044','GSE135222']))&p.endpoint.eq('No durable benefit')&p.feature.eq('six_module_proxy'))|(p.cohort.eq('GSE93157')&p.endpoint.eq('No objective response')&p.feature.eq('four_immune_module_proxy')))]
 forest(axes[6],p,['GSE126044\nSix modules','GSE135222\nSix modules','GSE93157\nFour immune modules'],zero=.5,estimate='auc',xlabel='AUC for nonresponse (95% CI)',color=C['blue']);auc_axis(axes[6])
 km(axes[7],b[b.cohort.eq('GSE135222')],'pfs_time_days','pfs_event','six_module_proxy')
 cap='**Figure 4. The post-treatment transcriptional state associates with pathology and recurrence.** (A) Six-component scores by pathological response in 212 patients. (B) Apparent, in-sample receiver-operating-characteristic curves for non-MPR; every feature uses the same 212 patients. (C) Composite-score and oppositely oriented PD-L1-category AUCs in 101 patients with both measurements; these analyses compare association within the discovery cohort. (D) Cox associations with recurrence-free survival in 154 patients with 30 events: unadjusted, adjusted for age/sex/histology/stage, and adjusted for the same clinical variables plus pathological response. (E) Recurrence-free survival grouped at the median score of the 212-person discovery set. (F) A fixed-direction six-module bulk proxy in the 16 pretreatment GSE126044 biopsies, using the original durable-benefit classification. (G) External AUCs with 2,000 stratified patient-bootstrap percentile confidence intervals. GSE126044 and GSE135222 use no durable clinical benefit as the positive endpoint; GSE93157 uses absence of objective response and only the four measurable immune modules, explicitly a different proxy. The complete six-module score is not estimable on the GSE93157 panel. (H) Progression-free survival in GSE135222, split at that cohort’s median proxy value; 27 patients and 21 events. Higher scores retain the nonresponse/risk direction throughout. Curves show observed Kaplan–Meier estimates; vertical ticks denote censoring and median splits are descriptive, without outcome-optimized cutoffs. Points in A and F denote patients; black summary lines show medians and interquartile ranges. AUC, area under the curve; SD, standard deviation. Supplementary Tables S7–S11 give all sample/event counts, confidence intervals, component comparisons and assay coverage.'
 save(fig,'Figure_4',cap,8,['sc_auc_benchmarks.csv','sc_pdl1_benchmarks.csv','sc_rfs_models.csv','sc_clinical_followup.csv','bulk_scores.csv','bulk_auc_models.csv'])
def supp_save(fig,num,title_,body,n,sources):save(fig,f'Supplementary_Figure_{num}',f'**Supplementary Figure S{num}. {title_}.** '+body,n,sources)
def make_supplement():
 # S1: biological-sample coverage and clinical distribution.
 fig,axes=setup(6);names=['Age by pathological response','Recovered immune cells','Mean library size in T/NK cells','Broad immune-cell proportions','Collection centres','PD-L1 categories']
 for i,(a,n) in enumerate(zip(axes,names)):title(a,string.ascii_uppercase[i],n)
 cl=pd.read_csv(W/'source_review/GSE243013_primary224_analysis_metadata.csv');violin(axes[0],cl,'age','Age (years)')
 totals=major[major.known_ICI].groupby('sampleID').n_cells.sum();v=cl.copy();v['log_cells']=np.log10(v.sampleID.map(totals));violin(axes[1],v,'log_cells','log10 recovered immune cells')
 v=tnk.copy();v['loglib']=np.log10(v.mean_library);violin(axes[2],v,'loglib','log10 mean library size (UMIs)')
 m=major[major.known_ICI].pivot(index='sampleID',columns='major_cell_type',values='n_cells');m=m.div(m.sum(axis=1),axis=0)*100;m=m.reindex(cl.sort_values('pathological_response').sampleID);heat(axes[3],m.T,m.shape[0]*[''],m.columns,label='Immune-cell proportion (%)',vlim=100,nonnegative=True);axes[3].set_xlabel('224 patients')
 counts=pd.crosstab(cl.center_S1,cl.pathological_response).reindex(columns=ORD);heat(axes[4],counts,ORD,counts.index,label='Patients',vlim=counts.max().max(),nonnegative=True)
 counts=pd.crosstab(cl.PDL1_TPS_category,cl.pathological_response).reindex(index=['<1%','1-49%','>=50%','unavailable'],columns=ORD);heat(axes[5],counts,ORD,['<1%','1–49%','≥50%','Unavailable'],label='Patients',vlim=counts.max().max(),nonnegative=True)
 supp_save(fig,1,'Patient characteristics and immune-cell coverage','(A) Age. (B) Total immune cells recovered per patient. (C) Mean T/NK library size. (D) Broad lineage proportions among all recovered immune cells. (E,F) Patient counts by collection centre and PD-L1 category. The primary source-matched cohort comprises 224 independent patient/sample identifiers; 103 have non-MPR, 40 MPR and 81 pCR. All source exclusions and missing clinical fields are retained in Supplementary Table S1. D–F use a white-to-vermilion scale for nonnegative proportions and counts. Distribution summaries follow Figure 1.',6,['sc_sample_map.csv','sc_major_scores.csv','GSE243013_primary224_analysis_metadata.csv'])
 # S2: all lineages rather than the selected main-figure view.
 fig,axes=setup(6);ms=['GLYCOLYSIS','ADENOSINE','CD8_EXHAUSTION','TREG_SUPPRESSION','CYTOTOXICITY','LACTATE_TRANSPORT'];z=con.query("analysis=='primary_resolved' and model=='HC3 adjusted OLS'")
 for i,(a,m) in enumerate(zip(axes,ms)):
  title(a,string.ascii_uppercase[i],module(m));v=z[z.module.eq(m)].sort_values('group');forest(a,v,v.group)
 supp_save(fig,2,'Response contrasts across resolved immune lineages','(A–F) Glycolysis, adenosine, exhaustion-associated, Treg-associated, cytotoxicity and lactate-transport programs. All evaluated lineages are shown, rather than selecting lineages by statistical significance. Points and horizontal lines show covariate-adjusted non-MPR differences and 95% HC3 confidence intervals. Each model uses one observation per patient and requires at least 100 cells of the specified lineage. Full 12-program results and panel-specific denominators are in Supplementary Table S3.',6,['sc_all_contrasts.csv'])
 # S3: complete subtype-program matrix, separated into biological families.
 z=con.query("analysis=='primary_subtype' and model=='HC3 adjusted OLS'");groups=[sorted(x for x in z.group.unique() if x.startswith(('CD4','CD8','T_g','NK','ILC'))),sorted(x for x in z.group.unique() if x.startswith(('Mφ','DC','pDC','cDC','mDC','Neu','Mast'))),[]]
 groups[2]=sorted(set(z.group)-set(groups[0])-set(groups[1]));fig,axes=plt.subplots(3,1,figsize=(183/25.4,265/25.4),gridspec_kw={'height_ratios':[len(g) for g in groups]},layout='constrained');fig.set_constrained_layout_pads(h_pad=.07,hspace=.12)
 mods=list(MOD)
 for i,(a,g,lab) in enumerate(zip(axes,groups,['T, NK and innate lymphoid states','Myeloid states','B-cell and plasma-cell states'])):
  title(a,string.ascii_uppercase[i],lab);m=z.pivot(index='group',columns='module',values='estimate').reindex(index=g,columns=mods);q=z.pivot(index='group',columns='module',values='q').reindex(index=g,columns=mods);im=heat(a,m,[module(x) for x in mods],[x.replace('_',' ') for x in g],vlim=.65,label='Adjusted non-MPR difference',signif=q,colorbar=False);plt.setp(a.get_xticklabels(),rotation=35,ha='right',fontsize=6.5);a.tick_params(axis='y',labelsize=7)
 cb=fig.colorbar(im,ax=list(axes),orientation='horizontal',fraction=.025,pad=.045,aspect=35);cb.set_label('Adjusted non-MPR difference')
 supp_save(fig,3,'Complete subtype-level transcriptional response profiles','(A–C) All estimable combinations of the 12 specified gene programs with T/NK/innate lymphoid, myeloid, and B/plasma-cell states. Original annotation labels are retained. OXPHOS, oxidative phosphorylation. Points mark q<0.05 across the complete subtype testing family; colors saturate at ±0.65 score units for legibility, and all uncapped estimates, intervals and denominators are in Supplementary Table S3. Each patient contributes at most one mean per subtype; a minimum of 30 cells is required. No new functional cell identity is inferred from a signature name.',3,['sc_all_contrasts.csv'])
 # S4: six primary components, across cohort/threshold sensitivities.
 fig,axes=setup(6);spec=[('T/NK cell','GLYCOLYSIS'),('T/NK cell','CD8_EXHAUSTION'),('T/NK cell','TREG_SUPPRESSION'),('T/NK cell','ADENOSINE'),('Myeloid cell','LACTATE_TRANSPORT'),('Myeloid cell','MYELOID_SUPPRESSION')];labs={'primary_major':'Primary','MPR_only':'Exclude pCR','combination_only':'Recorded chemotherapy','minimum_200_cells':'At least 200 cells','minimum_30_cells':'At least 30 cells','GEO_expanded':'GEO-expanded cohort'}
 for i,(a,(lin,m)) in enumerate(zip(axes,spec)):
  title(a,string.ascii_uppercase[i],(('T/NK ' if lin=='T/NK cell' else 'Myeloid ')+module(m).lower()).replace('Myeloid myeloid','Myeloid'));z=con[con.group.eq(lin)&con.module.eq(m)&con.model.eq('HC3 adjusted OLS')&con.analysis.isin(labs)].copy();forest(a,z,z.analysis.map(labs))
 supp_save(fig,4,'Cohort and cell-count sensitivity analyses','(A–F) The six discovery components are shown for the primary source-matched cohort, the MPR-only comparison excluding pCR, cases with recorded chemotherapy, minimum-cell thresholds of 200 and 30, and the expanded GEO-annotated cohort. Every estimate is a patient-level covariate-adjusted non-MPR contrast with a 95% HC3 interval. Inclusion rules and denominators are given per model in Supplementary Table S3. The expanded cohort is a sensitivity analysis, not a second validation cohort.',6,['sc_all_contrasts.csv'])
 # S5: clinical, composition, depth and centre sensitivity for all relevant programs.
 fig,axes=setup(6);s=pd.read_csv(A/'sc_composition_center_sensitivity.csv').assign(adjustment='composition_depth_center');cc=pd.concat([comp,s],ignore_index=True)
 for i,(a,m) in enumerate(zip(axes,['GLYCOLYSIS','ADENOSINE','TREG_SUPPRESSION','CD8_EXHAUSTION','LACTATE_TRANSPORT','CYTOTOXICITY'])):
  title(a,string.ascii_uppercase[i],module(m));v=cc[cc.module.eq(m)];forest(a,v,v.adjustment.map({'clinical':'Clinical','composition':'Add cell proportions','composition_depth':'Add library size','composition_depth_center':'Add collection centre'}))
 supp_save(fig,5,'Cell composition and technical covariate adjustment','(A–F) T/NK response contrasts after sequential adjustment for clinical covariates, measured Treg/NK/terminal-CD8 proportions, mean library size, and collection centre. Dots and intervals show estimates and 95% HC3 confidence intervals. Adenosine and Treg-associated contrasts persist after the complete adjustment, whereas broad glycolytic and exhaustion-associated contrasts attenuate. These are conditional associations rather than mediation estimates. Four patients with unreported stage are excluded from clinical regression models, leaving 217 T/NK observations. Supplementary Table S5 contains the full 12-program results and false-discovery rates.',6,['sc_composition_models.csv','sc_composition_center_sensitivity.csv'])
 # S6: spatial normalization and background sensitivity.
 fig,axes=setup(6);names=['Global expression across paired regions','Raw background and deposited expression','Deposited versus relative gene contrasts','RNA abundance and nuclei count','Paired gene contrasts across scales','Background-adjusted response interactions']
 for i,(a,n) in enumerate(zip(axes,names)):title(a,string.ascii_uppercase[i],n)
 d=pd.read_csv(A/'spatial_sample_scores_technical.csv');p=d.pivot(index=['patient','roi_key'],columns='panck',values='global_median').dropna().groupby('patient').mean()
 for _,r in p.iterrows():axes[0].plot([0,1],r[[0,1]],color=C['black'],alpha=.3,lw=.6)
 axes[0].scatter(np.zeros(len(p)),p[0],s=8,c=C['blue']);axes[0].scatter(np.ones(len(p)),p[1],s=8,c=C['orange']);axes[0].set_xticks([0,1],['PanCK−','PanCK+']);axes[0].set_ylabel('Median deposited expression')
 for k,col in [(0,C['blue']),(1,C['orange'])]:
  z=d[d.panck.eq(k)];axes[1].scatter(z.log_background,z.global_median,s=9,c=col,label='PanCK+' if k else 'PanCK−');axes[3].scatter(z.log_nuclei,z.global_median,s=9,c=col)
 axes[1].set(xlabel='log(1 + negative-probe mean)',ylabel='Median deposited expression');axes[1].legend(frameon=False);axes[3].set(xlabel='log(1 + nuclei count)',ylabel='Median deposited expression')
 g=pd.read_csv(A/'spatial_gene_compartments.csv');z=g.pivot(index='gene',columns='scale',values='estimate');scatter(axes[2],z.Deposited,z.Relative,xlabel='Deposited paired difference',ylabel='Relative paired difference')
 genes=['EPCAM','HK2','LDHB','SLC16A1','NT5E','CTLA4','CD3D','PTPRC'];v=g[g.gene.isin(genes)].rename(columns={'gene':'module','scale':'adjustment'});grouped_effect(axes[4],v,genes,['Deposited','Relative'])
 mm=pd.read_csv(A/'spatial_mixed_models.csv');bg=pd.read_csv(A/'spatial_background_sensitivity.csv');v=mm[(mm.term=='panck:responder')&(mm.adjustment=='Nuclei adjusted')&mm.module.str.endswith('_relative')].copy();v['adjustment']='Nuclei';bg=bg[bg.term.eq('panck:responder')].assign(adjustment='Nuclei + background');v=pd.concat([v,bg]);grouped_effect(axes[5],v,list(bg.module),['Nuclei','Nuclei + background']);axes[5].legend(loc='upper center',bbox_to_anchor=(.5,-.25),ncol=2,frameon=False,handlelength=1)
 supp_save(fig,6,'Spatial expression scaling and technical sensitivity','(A) Patient-averaged paired differences in the median deposited expression. (B,D) Area-level relations to the mean of 75 negative probes and nuclei counts; each point is an area of illumination, shown descriptively. (C) Gene-wise patient-paired differences before and after within-area median centering. (E) Selected genes on both scales. (F) Response interaction estimates after adding raw negative-probe background to the nested patient/ROI model. Lines in E,F are 95% confidence intervals. Median centering is an explicitly defined sensitivity transformation, not the original authors’ normalization. It resolves relative abundance within an area; the deposited scale is retained in all data tables.',6,['spatial_sample_scores_technical.csv','spatial_gene_compartments.csv','spatial_background_sensitivity.csv'])
 # S7: clinical association, common denominators and observed score composition.
 fig,axes=setup(6);names=['Component correlation structure','Component AUCs in the same patients','PD-L1 strata and composite scores','PD-L1 availability and score distribution','Prognostic association and pathological response','Proportional-hazards diagnostics']
 for i,(a,n) in enumerate(zip(axes,names)):title(a,string.ascii_uppercase[i],n)
 cols=[c for c in sc if c.startswith('z_')];labels=['Glycolysis','Exhaustion','Treg','Adenosine','Lactate','Myeloid'];heat(axes[0],sc[cols].corr(method='spearman'),labels,labels,vlim=1,label='Spearman correlation');plt.setp(axes[0].get_xticklabels(),rotation=45,ha='right')
 d=pd.read_csv(A/'sc_auc_benchmarks.csv');forest(axes[1],d,['Composite']+labels,zero=.5,estimate='auc',xlabel='AUC (95% CI), n=212',color=C['blue'])
 d=pd.read_csv(A/'sc_clinical_followup.csv');violin(axes[2],d,'score','Composite score',groups=['<1%','1-49%','>=50%'],x='PDL1_TPS_category',colors=[C['blue'],C['orange'],C['green']]);d['availability']=np.where(d.TPS_ordinal.notna(),'Available','Unavailable');violin(axes[3],d,'score','Composite score',groups=['Available','Unavailable'],x='availability',colors=[C['blue'],C['orange']])
 p=pd.read_csv(A/'sc_rfs_models.csv');p=p[p.feature.eq('non_MPR')];forest(axes[4],p,['Non-MPR unadjusted','Non-MPR clinical'],zero=1,estimate='hr',xlabel='Hazard ratio (95% CI)',color=C['blue']);axes[4].set_xscale('log');axes[4].set_xticks([1,2,4,8],['1','2','4','8']);axes[4].set_xticks([],minor=True)
 ph=pd.read_csv(A/'sc_rfs_ph_assumption.csv');v=ph[(ph.feature=='score_sd')&(ph.model=='clinical')];axes[5].scatter(v.p,np.arange(len(v)),c=C['blue'],s=16);axes[5].set_yticks(range(len(v)),v.term.replace({'score_sd':'Composite score','male':'Sex','LUSC':'Histology','stage3':'Stage','GLOBAL':'Global'}));axes[5].invert_yaxis();axes[5].set(xlim=(0,1),xlabel='Schoenfeld-residual test P value')
 supp_save(fig,7,'Clinical comparisons and recurrence analyses','(A) Correlations among the six discovery components. (B) Every component and the composite on the common 212-patient subset, with 2,000 patient-bootstrap 95% AUC intervals. (C,D) Score distributions by PD-L1 category and availability; seven measurements reported only as <1% are used categorically, not replaced by zero. Points in C,D denote patients; black summary lines show medians and interquartile ranges. (E) Pathological response associations with recurrence in all 157 patients with follow-up (31 events). (F) Scaled Schoenfeld-residual diagnostics for the five-covariate clinical Cox model. The continuous-score recurrence model uses 154 patients/30 events, whereas the response-only model can include 157/31. Full estimates, source values and common-subset C-index comparisons are in Supplementary Tables S7–S8.',6,['sc_auc_benchmarks.csv','sc_clinical_followup.csv','sc_rfs_models.csv','sc_rfs_ph_assumption.csv'])
 # S8: coverage and every external feature.
 fig,axes=setup(6);names=['Measured signature coverage','GSE126044 durable benefit','GSE135222 durable benefit','GSE93157 objective response','GSE126044 specimen sensitivity','GSE135222 progression associations']
 for i,(a,n) in enumerate(zip(axes,names)):title(a,string.ascii_uppercase[i],n)
 cov=pd.read_csv(A/'bulk_gene_coverage.csv');pv=cov.pivot(index='module',columns='cohort',values='coverage').reindex(index=list(MOD));heat(axes[0],pv,pv.columns,[wrap(module(m),20) for m in pv.index],vlim=1,label='Fraction of specified genes measured',nonnegative=True);axes[0].tick_params(axis='x',labelsize=6)
 b=pd.read_csv(A/'bulk_auc_models.csv');features=['six_module_proxy','four_immune_module_proxy','GLYCOLYSIS','ADENOSINE','CD8_EXHAUSTION','TREG_SUPPRESSION','MYELOID_SUPPRESSION','negative_CD274','negative_CYTOTOXICITY'];labs={'six_module_proxy':'Six-module proxy','four_immune_module_proxy':'Four immune modules','negative_CD274':'Opposite CD274','negative_CYTOTOXICITY':'Opposite cytotoxicity'}
 for a,co,endpoint in zip(axes[1:4],['GSE126044','GSE135222','GSE93157'],['No durable benefit','No durable benefit','No objective response']):
  v=b[(b.cohort==co)&(b.endpoint==endpoint)&(b.analysis=='All eligible patients')&b.feature.isin(features)&b.auc.notna()];forest(a,v,[labs.get(f,module(f)) for f in v.feature],zero=.5,estimate='auc',xlabel='AUC for nonresponse (95% CI)',color=C['blue']);auc_axis(a)
 v=b[(b.cohort=='GSE126044')&(b.feature=='six_module_proxy')];forest(axes[4],v,[('Durable benefit' if x=='No durable benefit' else 'Objective response')+' / '+('all' if y=='All eligible patients' else 'fresh') for x,y in zip(v.endpoint,v.analysis)],zero=.5,estimate='auc',xlabel='AUC for nonresponse (95% CI)',color=C['blue']);auc_axis(axes[4])
 s=pd.read_csv(A/'bulk_survival_models.csv');v=s[s.feature.isin(features)];forest(axes[5],v,[labs.get(f,module(f)) for f in v.feature],zero=1,estimate='hr',xlabel='Hazard ratio per SD (95% CI)',color=C['blue']);axes[5].set_xscale('log');axes[5].set_xticks([.5,1,2],['0.5','1','2']);axes[5].set_xticks([],minor=True)
 supp_save(fig,8,'Assay coverage and external bulk-expression comparisons','(A) Coverage of all 12 original gene sets on each assay. GSE93157 measures one of 12 glycolysis genes and none of the three lactate-transport genes; its complete six-module proxy is therefore undefined. (B–D) Fixed-direction AUCs of the evaluable proxies and component benchmarks, with percentile confidence intervals from 2,000 stratified patient bootstraps. CD274 and cytotoxicity are explicitly reversed only to express their prespecified favorable direction as nonresponse scores. (E) The six-module proxy in all 16 versus 11 fresh GSE126044 specimens under durable-benefit and objective-response definitions. All five FFPE specimens are nonresponders, and the reverse association persists among fresh specimens. (F) Continuous-score Cox estimates for GSE135222 progression (27 patients, 21 events), with 95% confidence intervals. Points without intervals denote degenerate empirical bootstrap distributions; their confidence intervals are marked non-estimable in Supplementary Table S11. No response-derived feature weight, cutoff or direction was fitted. All data and model statuses are in Supplementary Tables S9–S11.',6,['bulk_gene_coverage.csv','bulk_auc_models.csv','bulk_survival_models.csv'])
 # S9: host metabolic context, no sparse-case logistic claims.
 fig,axes=setup(6);names=['Metabolic score and inflammation','Metabolic score and white cells','Metabolic score and NLR','CRP sensitivity analyses','White-cell sensitivity analyses','NLR sensitivity analyses']
 for i,(a,n) in enumerate(zip(axes,names)):title(a,string.ascii_uppercase[i],n)
 n=pd.read_csv(A/'nhanes_analysis_data.csv');r=pd.read_csv(A/'nhanes_models.csv')
 for a,y,lab in zip(axes[:3],['log_crp','log_wbc','log_nlr'],['log(1 + CRP, mg/L)','log(1 + WBC, 10⁹/L)','log(1 + NLR)']):scatter(a,n.metabolic_only,n[y],xlabel='Metabolic-only score (SD)',ylabel=lab)
 labs={'Primary':'Primary','Ever_smokers':'Ever smokers','Age60':'Age ≥60 years','CRP10':'CRP ≤10 mg/L','Cycle2015':'2015–2016','Cycle2017':'2017–2018'}
 for a,y in zip(axes[3:],['log_crp','log_wbc','log_nlr']):
  v=r[(r.outcome==y)&(r.exposure=='metabolic_only')];forest(a,v,v.scope.map(labs),xlabel='Survey-weighted beta per SD (95% CI)')
 supp_save(fig,9,'Metabolic-inflammatory associations in the adult population','(A–C) Individual observations among 3,045 adults aged at least 40 years with complete score components and smoking status in NHANES 2015–2018. These scatterplots are descriptive; inference uses the complex survey design. The metabolic-only score contains standardized triglyceride–glucose index and body mass index, excluding inflammatory outcomes. (D–F) Survey-weighted covariate-adjusted associations in the primary sample and predefined restrictions, with 95% Taylor-linearized confidence intervals. Primary regression uses 2,685 complete covariate records. Two-year fasting weights are divided by two and masked strata/PSUs are retained. CRP, C-reactive protein, is in mg/L; WBC, white-cell count, is in 10⁹/L; NLR, neutrophil-to-lymphocyte ratio, is dimensionless. Nine reported lung-cancer cases among 2,678 participants with known cancer status are retained descriptively, without a sparse multivariable odds-ratio model. Supplementary Tables S16–S17 give selection, missingness, non-overlap models and denominators.',6,['nhanes_analysis_data.csv','nhanes_models.csv'])
 # S10: full genetic screen, with no target ranking.
 fig,axes=setup(4);names=['All 7,224 immune-trait tests','Instrument strength','Locus-to-gene distances','Nominal P values across genetic exposures']
 for i,(a,n) in enumerate(zip(axes,names)):title(a,string.ascii_uppercase[i],n)
 m=pd.read_csv(A/'mr_complete_results.csv');p=np.sort(m.p);expected=-np.log10((np.arange(1,len(p)+1)-.5)/len(p));observed=-np.log10(p);scatter(axes[0],expected,observed,xlabel='Expected −log10 P',ylabel='Observed −log10 P',identity=True)
 ins=pd.read_csv(A/'mr_instrument_locus_annotations.csv').sort_values(['target','SNP','exposure_group']);axes[1].scatter(ins.F_stat,np.arange(len(ins)),s=12,c=C['blue']);axes[1].set_yticks(range(len(ins)),ins.target+' / '+ins.SNP);axes[1].invert_yaxis();axes[1].set_xscale('log');axes[1].set_xlabel('Instrument F statistic')
 axes[2].scatter(ins.distance_to_gene_bp/1000,np.arange(len(ins)),s=12,c=C['blue']);axes[2].set_yticks(range(len(ins)),ins.target+' / '+ins.SNP);axes[2].invert_yaxis();axes[2].set_xlabel('Distance outside annotated gene (kb)')
 for i,(target,z) in enumerate(m.groupby('original_target_label',sort=True)):
  axes[3].scatter(np.full(len(z),i)+rng.uniform(-.25,.25,len(z)),-np.log10(z.p),s=2,alpha=.22,c=C['blue'],linewidths=0,rasterized=True)
 axes[3].set_xticks(range(m.original_target_label.nunique()),sorted(m.original_target_label.unique()),rotation=45,ha='right');axes[3].set_ylabel('−log10 P');axes[3].set_xlabel('Original locus labels (alphabetical)')
 supp_save(fig,10,'Complete metabolic-trait genetic screen','(A) Quantile–quantile display of all 7,224 tests. (B) Strength of the 13 instrument–exposure records, representing 12 distinct variants. (C) Distance from each instrument to the GRCh37 annotated gene interval; zero denotes an intragenic position. Names are original locus labels, not molecularly validated target assignments. (D) All uncorrected P values by locus, ordered alphabetically. None of the tests survives Benjamini–Hochberg correction across the full family (minimum q=0.9934). Estimates are changes in standardized immune traits per unit of the relevant genetically proxied metabolic trait, not effects of therapeutic inhibition. All instruments, harmonized variants, effect estimates, two-variant heterogeneity tests and diagnostic applicability are supplied in Supplementary Tables S18–S20.',4,['mr_complete_results.csv','mr_instrument_locus_annotations.csv'])

if __name__=='__main__':
 make_main123();make_main4();make_supplement()
 # Figure 3 points to the final sequential supplementary numbering.
 for r in registry:r['caption']=r['caption'].replace('Supplementary Figure S8.','Supplementary Figure S6.') if r['figure']=='Figure_3' else r['caption']

 # Complete standalone legends while retaining the source annotation names.
 abbreviations={
 'Figure_1':'SD, standard deviation; NK, natural killer; HC3, heteroskedasticity-consistent type 3. In F, vermilion and blue denote positive and negative coefficients, respectively; colors do not encode statistical significance.',
 'Figure_2':'Bm, memory B cell; Bn, naive B cell; NK, natural killer; Treg, regulatory T cell; MPR, major pathological response; pCR, pathological complete response; HC3, heteroskedasticity-consistent type 3. Coefficient colors do not encode statistical significance.',
 'Figure_3':'Treg, regulatory T cell. Coefficient colors do not encode statistical significance.',
 'Figure_4':'MPR, major pathological response; pCR, pathological complete response; PD-L1, programmed death-ligand 1; Treg, regulatory T cell.',
 'Supplementary_Figure_1':'UMI, unique molecular identifier; NK, natural killer; PD-L1, programmed death-ligand 1; MPR, major pathological response; pCR, pathological complete response. Points denote patients; black summaries show medians and interquartile ranges. Green, blue and vermilion denote pCR, MPR and non-MPR.',
 'Supplementary_Figure_2':'NK, natural killer; Treg, regulatory T cell; MPR, major pathological response; HC3, heteroskedasticity-consistent type 3.',
 'Supplementary_Figure_3':'Annotation key: Tem, effector-memory T cell; Tfh, follicular helper T cell; Tm, memory T cell; Tn, naive T cell; Tex, exhausted T cell; Trm, tissue-resident memory T cell; Treg, regulatory T cell; Th1, type 1 helper T cell; MAIT, mucosal-associated invariant T cell; ILC3, type 3 innate lymphoid cell; gdT, gamma-delta T cell; NK, natural killer cell; Mφ, macrophage; Neu, neutrophil; cDC1/cDC2, conventional dendritic-cell types 1/2; pDC, plasmacytoid dendritic cell; Bm/Bn, memory/naive B cell; GCB, germinal-centre B cell. The source labels mDC_LAMP3 and prf_MKI67 identify LAMP3-positive dendritic-cell and MKI67-positive clusters, respectively. Gene suffixes retain deposited marker labels.',
 'Supplementary_Figure_4':'MPR, major pathological response; pCR, pathological complete response; NK, natural killer; HC3, heteroskedasticity-consistent type 3; GEO, Gene Expression Omnibus.',
 'Supplementary_Figure_5':'Treg, regulatory T cell; NK, natural killer; HC3, heteroskedasticity-consistent type 3.',
 'Supplementary_Figure_6':'ROI, region of interest; PanCK, pan-cytokeratin; OXPHOS, oxidative phosphorylation; Treg, regulatory T cell.',
 'Supplementary_Figure_7':'AUC, area under the curve; PD-L1, programmed death-ligand 1; MPR, major pathological response; Treg, regulatory T cell.',
 'Supplementary_Figure_8':'AUC, area under the curve; SD, standard deviation; OXPHOS, oxidative phosphorylation; FFPE, formalin-fixed paraffin-embedded.',
 'Supplementary_Figure_9':'NHANES, National Health and Nutrition Examination Survey; PSU, primary sampling unit; SD, standard deviation.'}
 for r in registry:
  r['caption']+=' '+abbreviations.get(r['figure'],'')
  if r['figure'] in ['Supplementary_Figure_2','Supplementary_Figure_4','Supplementary_Figure_5','Supplementary_Figure_9']:
   r['caption']+=' Vermilion and blue denote positive and negative coefficients, respectively; colors do not encode statistical significance.'
 for r in registry:(O/(r['figure']+'_legend.md')).write_text(r['caption']+'\n',encoding='utf8')
 (O/'figure_registry.json').write_text(json.dumps(registry,indent=2),encoding='utf8')
