"""Correlate every competition control cell with five external control profiles.

Primary: Pearson correlation of log1p(CPM) over one shared gene axis, including
zeros. Spearman (average ties), source-pair axes, and positive-reference genes
are sensitivity analyses. No context identity or perturbation effect is inferred.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,shutil
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import prepare_reference_controls as prep
ROOT=Path(__file__).resolve().parent
REFERENCES=['K562','HCT116','HEK293T','H1','CD4']

def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
 return h.hexdigest()

def pearson_rows(values,profiles):
 """Rows by profiles, matching independent numpy.corrcoef including constants."""
 a=np.asarray(values,dtype=np.float64);b=np.asarray(profiles,dtype=np.float64)
 a=a-a.mean(axis=1,keepdims=True);b=b-b.mean(axis=1,keepdims=True)
 denominator=np.sqrt((a*a).sum(axis=1))[:,None]*np.sqrt((b*b).sum(axis=1))[None,:]
 return np.divide(a@b.T,denominator,out=np.full((len(a),len(b)),np.nan),where=denominator>0).clip(-1,1)

def normalized_log(raw):
 raw=np.asarray(raw,dtype=np.float64);totals=raw.sum(axis=-1,keepdims=True)
 if (totals<=0).any() or not np.isfinite(raw).all() or (raw<0).any():raise ValueError('Invalid normalization input')
 return np.log1p(raw*1e6/totals)

def run(batch_size=256):
 spec=json.loads((ROOT/'manifest.json').read_text());genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.astype(str).to_numpy()
 output=ROOT/'analysis/cell_reference_correlations';output.mkdir(parents=True,exist_ok=True)
 profiles={};records={}
 for name in REFERENCES:
  path=ROOT/f'reference/control_profiles/{name}.csv';record=json.loads(path.with_suffix('.json').read_text())
  if digest(path)!=record['profile_sha256']:raise ValueError(f'{name} profile checksum changed')
  frame=pd.read_csv(path)
  if frame.gene_name.duplicated().any():raise ValueError('Duplicate profile symbols')
  profiles[name]=frame.set_index('gene_name').control_count_sum;records[name]=record
  selected=path.with_name(name+'_selected_controls.csv')
  shutil.copyfile(selected,output/f'reference_selected_controls_{name}.csv')
  records[name]['selected_controls_sha256']=digest(selected)
 common=np.array([g for g in genes if all(g in p.index for p in profiles.values())])
 if len(common)<1000:raise ValueError('Insufficient common gene coverage')
 positions={g:i for i,g in enumerate(genes)};common_positions=np.array([positions[g] for g in common])
 raw=np.stack([profiles[n].loc[common].to_numpy() for n in REFERENCES]);logged=normalized_log(raw)
 reference_ranks=rankdata(logged,axis=1,method='average')
 positive=(raw>0).all(axis=0)
 common_set=set(common)
 axis=pd.DataFrame({'gene_name':genes,'common_to_all_five':[g in common_set for g in genes]})
 for n,p in profiles.items():axis['measured_'+n]=np.isin(genes,p.index)
 axis.to_csv(output/'gene_axis_coverage.csv',index=False)
 pd.DataFrame({'gene_name':common,**{n:raw[i] for i,n in enumerate(REFERENCES)}}).to_csv(output/'reference_control_profiles_common_genes.csv',index=False)
 pd.DataFrame(pearson_rows(logged,logged),index=REFERENCES,columns=REFERENCES).to_csv(output/'reference_profile_correlations.csv',index_label='reference')
 all_frames=[];summaries=[];aggregate=[];context_log={};inputs=[];guide_records=[];condition_summary=[]
 for context in spec['contexts']:
  path=ROOT/f'context_{context}.h5ad';obj=ad.read_h5ad(path,backed='r')
  if obj.n_obs!=spec['per_context'][context]['control_cells'] or not np.array_equal(obj.var_names.to_numpy(),genes):raise ValueError('Competition axes differ')
  if set(obj.obs.target_gene.astype(str))!={spec['control_label']} or set(obj.obs.context.astype(str))!={context}:raise ValueError('Unexpected competition labels')
  n=obj.n_obs;result=np.empty((n,len(REFERENCES)));spearman=result.copy();pair=result.copy();filtered=result.copy();libs=np.zeros(n);detected=np.zeros(n);sums=np.zeros(len(genes))
  pair_axes={r:np.array([positions[g] for g in genes if g in profiles[r].index]) for r in REFERENCES}
  pair_profiles={r:normalized_log(profiles[r].loc[genes[pair_axes[r]]].to_numpy()[None,:]) for r in REFERENCES}
  extra_names=['CD4_Rest','CD4_Stim8hr','CD4_Stim48hr'];extra_profiles=[]
  for r in extra_names:
   p=pd.read_csv(ROOT/f'reference/control_profiles/{r}.csv').set_index('gene_name').control_count_sum
   extra_profiles.append(p.loc[common].to_numpy())
  extra_profiles=normalized_log(np.stack(extra_profiles));extra=np.empty((n,3))
  for start in range(0,n,batch_size):
   stop=min(start+batch_size,n);x=sparse.csr_matrix(obj.X[start:stop],dtype=np.float64)
   if not np.isfinite(x.data).all() or (x.data<0).any() or (x.data!=np.floor(x.data)).any():raise ValueError('Invalid competition counts')
   x.sum_duplicates();x.eliminate_zeros();libs[start:stop]=np.asarray(x.sum(axis=1)).ravel();detected[start:stop]=np.diff(x.indptr)
   sums+=np.asarray(x.sum(axis=0)).ravel();counts=x[:,common_positions].toarray();log=normalized_log(counts)
   result[start:stop]=pearson_rows(log,logged)
   spearman[start:stop]=pearson_rows(rankdata(log,axis=1,method='average'),reference_ranks)
   filtered[start:stop]=pearson_rows(normalized_log(counts[:,positive]),normalized_log(raw[:,positive]))
   extra[start:stop]=pearson_rows(log,extra_profiles)
   for j,r in enumerate(REFERENCES):pair[start:stop,j]=pearson_rows(normalized_log(x[:,pair_axes[r]].toarray()),pair_profiles[r]).ravel()
  frame=pd.DataFrame({'context':context,'observed_target_gene':spec['control_label'],'cell_id':obj.obs_names.astype(str),'ntc_id':obj.obs.ntc_id.astype(str).to_numpy(),'total_counts':libs.astype(np.int64),'detected_genes':detected.astype(np.int64),'common_gene_count':len(common)})
  for j,r in enumerate(REFERENCES):
   frame['pearson_log1p_cpm_'+r]=result[:,j];frame['spearman_'+r]=spearman[:,j];frame['pearson_pairwise_axis_'+r]=pair[:,j];frame['pearson_positive_reference_genes_'+r]=filtered[:,j]
   values=result[:,j];quantiles=np.quantile(values,[.05,.25,.5,.75,.95]);winners=np.argmax(result,axis=1)
   summaries.append({'context':context,'reference':r,'cells':n,'common_genes':len(common),'mean_pearson':values.mean(),'std_pearson':values.std(),'p05_pearson':quantiles[0],'p25_pearson':quantiles[1],'median_pearson':quantiles[2],'p75_pearson':quantiles[3],'p95_pearson':quantiles[4],'median_spearman':np.median(spearman[:,j]),'top_reference_cells':int((winners==j).sum()),'top_reference_fraction':float((winners==j).mean()),'top_reference_spearman_fraction':float((np.argmax(spearman,axis=1)==j).mean()),'correlation_with_log_library':float(np.corrcoef(values,np.log1p(libs))[0,1]),'median_pairwise_axis_pearson':float(np.median(pair[:,j])),'median_positive_reference_genes_pearson':float(np.median(filtered[:,j]))})
   for guide,group in frame.groupby('ntc_id',observed=True):guide_records.append({'context':context,'ntc_id':guide,'reference':r,'cells':len(group),'mean_pearson':group['pearson_log1p_cpm_'+r].mean()})
  for j,r in enumerate(extra_names):
   frame['pearson_log1p_cpm_'+r]=extra[:,j]
   condition_summary.append({'context':context,'reference':r,'cells':n,'common_genes':len(common),'median_pearson':float(np.median(extra[:,j])),'p05_pearson':float(np.quantile(extra[:,j],.05)),'p95_pearson':float(np.quantile(extra[:,j],.95))})
  frame['top_reference_pearson']=np.array(REFERENCES)[np.argmax(result,axis=1)]
  frame['top_reference_spearman']=np.array(REFERENCES)[np.argmax(spearman,axis=1)]
  sorted_r=np.sort(result,axis=1);frame['top_two_pearson_margin']=sorted_r[:,-1]-sorted_r[:,-2]
  if not np.isfinite(frame.select_dtypes('number').to_numpy()).all():raise ValueError('Undefined per-cell statistic')
  all_frames.append(frame)
  context_log[context]=normalized_log(sums[common_positions][None,:])[0]
  bulk=pearson_rows(context_log[context][None,:],logged)[0]
  bulk_spearman=pearson_rows(rankdata(context_log[context])[None,:],reference_ranks)[0]
  for j,r in enumerate(REFERENCES):aggregate.append({'context':context,'reference':r,'common_genes':len(common),'pooled_control_pearson':bulk[j],'pooled_control_spearman':bulk_spearman[j]})
  obj.file.close();inputs.append({'file':path.name,'sha256':digest(path),'cells':n})
  print(context,'completed all',n,'cells',flush=True)
 combined=pd.concat(all_frames,ignore_index=True);combined.to_csv(output/'per_cell_correlations.csv',index=False)
 summary=pd.DataFrame(summaries);summary.to_csv(output/'context_reference_summary.csv',index=False)
 pd.DataFrame(aggregate).to_csv(output/'pooled_context_correlations.csv',index=False)
 pd.DataFrame(condition_summary).to_csv(output/'cd4_condition_summary.csv',index=False)
 pd.DataFrame(guide_records).to_csv(output/'guide_reference_correlations.csv',index=False)
 contexts=list(context_log);context_matrix=np.stack(list(context_log.values()))
 pd.DataFrame(pearson_rows(context_matrix,context_matrix),index=contexts,columns=contexts).to_csv(output/'context_profile_correlations.csv',index_label='context')
 report={'created_utc':datetime.now(timezone.utc).isoformat(),'competition_cells':len(combined),'common_genes':len(common),'positive_in_all_reference_genes':int(positive.sum()),'method':'Pearson across genes of log1p(CPM), normalized on exactly the same common gene axis; every competition control cell versus pooled non-targeting reference profiles','spearman':'Average ranks with tied zeros retained','reference_order':REFERENCES,'reference_provenance':records,'extra_CD4_condition_provenance':{r:json.loads((ROOT/f'reference/control_profiles/{r}.json').read_text()) for r in ['CD4_Rest','CD4_Stim8hr','CD4_Stim48hr']},'competition_inputs':inputs,'script_sha256':digest(Path(__file__)),'reference_script_sha256':digest(ROOT/'prepare_reference_controls.py'),'versions':{'numpy':np.__version__,'pandas':pd.__version__,'anndata':ad.__version__},'validation':{'expected_rows':sum(spec['per_context'][c]['control_cells'] for c in contexts),'unique_context_cell_ids':not combined[['context','cell_id']].duplicated().any(),'all_numeric_values_finite':bool(np.isfinite(combined.select_dtypes('number').to_numpy()).all()),'all_correlations_in_range':bool(combined.filter(regex='^(pearson|spearman)').abs().le(1+1e-12).all().all())},'limitations':['Basal correlation does not establish cell-line identity or response-transfer accuracy.','Gene symbols are exact-matched; duplicates in reference profiles are summed; unmeasured genes are not zero-imputed.','Reference controls are sampled or batch-limited; CD4 pooled profile mixes conditions.','Single-cell zeros, sequencing depth and platform differences affect similarity.','Cell rows across A/B/C are not paired observations.'],'outputs':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(output.glob('*.csv'))}}
 (output/'analysis_report.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
 plot(summary,output)
 write_results(summary,aggregate,condition_summary,report)
 print(summary.pivot(index='context',columns='reference',values='median_pearson').to_string(),flush=True)
 return report

def markdown_table(frame):
 columns=list(frame.columns);lines=['| '+' | '.join(columns)+' |','| '+' | '.join(['---']*len(columns))+' |']
 for row in frame.itertuples(index=False,name=None):
  lines.append('| '+' | '.join(f'{x:.4f}' if isinstance(x,(float,np.floating)) else str(x) for x in row)+' |')
 return '\n'.join(lines)

def write_results(summary,aggregate,conditions,report):
 median=summary.pivot(index='context',columns='reference',values='median_pearson').reindex(columns=REFERENCES).reset_index()
 winner=summary.pivot(index='context',columns='reference',values='top_reference_fraction').reindex(columns=REFERENCES)*100
 winner=winner.reset_index()
 pooled=pd.DataFrame(aggregate).pivot(index='context',columns='reference',values='pooled_control_pearson').reindex(columns=REFERENCES).reset_index()
 sensitivity=summary[['context','reference','median_spearman','median_pairwise_axis_pearson','median_positive_reference_genes_pearson']]
 source=pd.DataFrame([{'reference':n,'profile_genes':report['reference_provenance'][n]['genes'],'control_cells':report['reference_provenance'][n]['control_cells'],'scope':report['reference_provenance'][n].get('scope','K562 core-control pseudobulk groups')} for n in REFERENCES])
 text='# Control-cell reference correlation results\n\n'
 text+=f"Generated {report['created_utc']}. All **{report['competition_cells']:,} competition control cells** were analyzed: 18,400 each from A/B/C. Primary comparisons use **{report['common_genes']:,} shared gene symbols**. No predicted cells or 2026 perturbed responses enter the calculation.\n\n"
 text+='## Median per-cell Pearson correlations\n\n'+markdown_table(median)+'\n\n'
 text+='Each value summarizes 18,400 separate cell-to-reference correlations of log1p(CPM), normalized on the same common axis. Full distributions and 5th/95th percentiles are in `context_reference_summary.csv`.\n\n'
 text+='## Cells closest to each reference (%)\n\n'+markdown_table(winner)+'\n\n'
 text+='These percentages describe the largest correlation among the five candidate profiles, not known cell identities or probabilities. The per-cell CSV retains the margin between first and second place and Spearman-based choices.\n\n'
 text+='## Pooled context profiles\n\n'+markdown_table(pooled)+'\n\n'
 text+='Pooling counts reduces single-cell sampling zeros; pooled correlations are not averages of the per-cell correlations. A/B/C pooled profiles are also compared directly in `context_profile_correlations.csv`.\n\n'
 text+='## Sensitivity checks\n\n'+markdown_table(sensitivity)+'\n\n'
 text+=f"Positive-reference sensitivity uses {report['positive_in_all_reference_genes']:,} genes with positive pooled counts in every reference. Pair-specific gene axes differ between references and cannot provide a fair closest-reference ranking.\n\n"
 text+='### CD4 conditions, median per-cell Pearson\n\n'+markdown_table(pd.DataFrame(conditions))+'\n\n'
 text+='## Reference scope\n\n'+markdown_table(source)+'\n\n'
 text+='H1 uses 1,536 source controls across 48 batches. CD4 uses 384 source pseudobulks across 12 donor/condition/run strata, representing 34,620 cells. HCT116/HEK293T use first-batch controls. Thus this is a reproducible sampled-reference analysis, not a complete-atlas comparison.\n\n'
 text+='## Reproduction and verification\n\nSee [CELL_REFERENCE_CORRELATIONS.md](CELL_REFERENCE_CORRELATIONS.md) for source definitions, extraction, normalization, exact commands and output schemas. `analysis/cell_reference_correlations/analysis_report.json` records source identities, parameters and hashes. Run `validate_cell_reference_correlations.py` to check complete row/identifier coverage, summaries and 45 actual-cell/reference pairs against independent NumPy/SciPy calculations; its result is saved as `validation.json`.\n\n'
 text+='Basal similarity can be driven by housekeeping expression, shared zeros, sequencing depth and platform effects. It does not prove that a reference will predict CRISPRi responses accurately. Response-transfer weights require a separate held-out cell-line experiment.\n'
 (ROOT/'CELL_REFERENCE_CORRELATION_RESULTS.md').write_text(text,encoding='utf-8')

def plot(summary,output):
 fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
 colors=['#4878a8','#f58518','#54a24b']
 for i,c in enumerate(['A','B','C']):
  sub=summary[summary.context==c].set_index('reference').loc[REFERENCES];x=np.arange(5)+.24*(i-1)
  median=sub.median_pearson.to_numpy();axes[0].errorbar(x,median,yerr=np.stack([median-sub.p05_pearson,sub.p95_pearson-median]),fmt='o',capsize=3,label='Context '+c,color=colors[i])
  axes[1].bar(x,100*sub.top_reference_fraction,width=.23,label='Context '+c,color=colors[i])
 for ax in axes:ax.set_xticks(np.arange(5),REFERENCES);ax.spines[['top','right']].set_visible(False)
 axes[0].set_ylabel('Pearson r across shared genes');axes[0].set_title('All cells: median and 5–95% interval');axes[0].legend()
 axes[1].set_ylabel('Cells with highest correlation (%)');axes[1].set_title('Closest reference among these five')
 fig.suptitle('Control expression similarity — not inferred cell identities')
 fig.savefig(output/'correlation_overview.png',dpi=170);plt.close(fig)

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--batch-size',type=int,default=256);args=parser.parse_args()
 if args.batch_size<1:parser.error('batch-size must be positive')
 run(args.batch_size)
