"""Extend genuine public responses to the final panel, preserving every reused row exactly."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd,h5py
from scipy import sparse,stats
import build_submission as b
from fit_v7_response import read_cache,summarize_raw,raw_inputs
from build_cd4_model import load_subset,eligible
from fit_lineage_submission import bh_valid
from threadpoolctl import threadpool_limits
R=Path(__file__).resolve().parent;O=R/'reference/final_v91';A=R/'experiments/final_v91/analysis'
KEYS=['truth','errors','valid','q','knockdown','row_targets','row_contexts','row_families']
def main():
 assert not (A/'source_cache.npz').exists();
 cd4receipt=json.loads((A/'cd4_collection.json').read_text())
 for name,h in cd4receipt['combined'].items():assert b.sha256(O/'cd4'/name)==h,name
 old=read_cache();m,g,t=b.read_spec(R/'final');np.testing.assert_array_equal(g,old['genes']);new=np.array(sorted(set(t)-set(old['targets'])));tp={x:i for i,x in enumerate(t)};gp={x:i for i,x in enumerate(g)}
 rows=old['targets'][old['row_targets']];take=np.isin(rows,t);parts={k:[old[k][take].copy()] for k in KEYS};parts['row_targets']=[np.array([tp[x] for x in rows[take]])];nold=int(take.sum());provenance=[]
 def append(ci,family,target,y,se,q,valid,kd,n):
  assert target in tp and target in new and y.shape==se.shape==q.shape==valid.shape==(len(g),)
  valid=valid.copy();valid[gp[target]]=False
  for k,value in [('truth',np.clip(y,-3,3).astype(np.float32)),('errors',se.astype(np.float32)),('q',q.astype(np.float32)),('valid',valid),('knockdown',float(kd)),('row_targets',tp[target]),('row_contexts',ci),('row_families',family)]:parts[k].append(np.asarray([value]))
  provenance.append(dict(context_index=ci,family=family,target_gene=str(target),source_cells=int(n),modeled_genes=int(valid.sum())))
 # Match the original full GWPS pseudobulk estimator, including effective-guide selection.
 path=R/'reference/K562_gwps_raw_bulk_01.h5ad';legacy=json.loads((R/'experiments/v4_lineage_response/analysis/k562_source/reference_model.json').read_text());assert b.sha256(path)==legacy['source_file_sha256']
 with h5py.File(path) as f:
  labels=f['obs/gene_transcript'].asstr()[:];nc=f['obs/num_cells_filtered'][:].astype(float);fold=f['obs/fold_expr'][:];core=f['obs/core_control'][:].astype(bool);col=f['var/gene_name'];symbols=f[col.attrs['categories']].asstr()[:][col[:]];matrix=f['X'][:].astype(float)
 sg=np.array(sorted(set(symbols)&set(g)));sp={v:i for i,v in enumerate(sg)};cols=np.flatnonzero(np.isin(symbols,sg));mapping=sparse.csr_matrix((np.ones(len(cols)),(cols,[sp[symbols[i]] for i in cols])),shape=(len(symbols),len(sg)));matrix=mapping.T.dot(matrix.T).T;dest=np.array([gp[x] for x in sg]);names=np.array([s.split('_',2)[1] if len(s.split('_',2))==3 and s.split('_',1)[0].isdigit() else '' for s in labels]);ctrl=(names=='non-targeting')&core&(nc>=30);base=np.average(matrix[ctrl],axis=0,weights=nc[ctrl]);bcpm=base*1e6/base.sum()
 for target in new:
  rr=np.flatnonzero((names==target)&(nc>=30)&np.isfinite(fold)&(fold<=.5))
  if not len(rr):continue
  mean=np.average(matrix[rr],axis=0,weights=nc[rr]);cpm=mean*1e6/mean.sum();n=nc[rr].sum();y0=np.clip(np.log2((cpm+1)/(bcpm+1)),-3,3);se0=np.sqrt(1/np.maximum(n*mean,1)+1/np.maximum(nc[ctrl].sum()*base,1))/np.log(2);q0=bh_valid(2*stats.norm.sf(abs(y0/se0))[None,:],(bcpm>5)[None,:])[0];y=np.zeros(len(g));se=np.full(len(g),1e6);q=np.ones(len(g));valid=np.zeros(len(g),bool);y[dest]=y0;se[dest]=se0;q[dest]=q0;valid[dest]=np.isfinite(se0)&(se0>0)&(se0<1e5)&np.isfinite(old['profiles'][0,dest]);kd=float(np.clip(1-2**min(y[gp[target]],0),0,.99));append(0,'K562',target,y,se,q,valid,kd,n)
 del matrix,mapping;gc.collect();print('Final source K562 summarized',flush=True)
 # Published CD4 merged estimates, same original eligibility/missingness convention.
 meta,dd=load_subset(O/'cd4','merged');src=np.array([i for i,x in enumerate(dd['genes']) if x in gp]);dest=np.array([gp[dd['genes'][i]] for i in src]);good=eligible(meta)
 for ci,condition in enumerate(['Rest','Stim8hr','Stim48hr'],1):
  for ri in np.flatnonzero(good&meta.target_contrast_gene_name.isin(new).to_numpy()&(meta.culture_condition.to_numpy()==condition)):
   row=meta.iloc[ri];y=np.zeros(len(g));se=np.full(len(g),1e6);q=np.ones(len(g));y[dest]=np.nan_to_num(np.clip(dd['log_fc'][ri,src],-3,3));se[dest]=np.nan_to_num(dd['lfcSE'][ri,src],nan=1e6);q[dest]=np.nan_to_num(dd['adj_p_value'][ri,src],nan=1);valid=np.isfinite(se)&(se>0)&(se<1e5)&np.isfinite(old['profiles'][ci]);valid[dest]&=np.isfinite(dd['log_fc'][ri,src]);kd=float(np.clip(1-2**min(y[gp[row.target_contrast_gene_name]],0),0,.99));append(ci,'CD4',row.target_contrast_gene_name,y,se,q,valid,kd,row.n_cells_target)
 del dd;gc.collect();print('Final source CD4 summarized',flush=True)
 for ci,line in [(5,'HCT116'),(6,'HEK293T')]:
  folder=O/'orion'/line;rec=json.loads((folder/'extraction.json').read_text());assert b.sha256(folder/'selected_counts.npz')==rec['counts_sha256'];assert b.sha256(folder/'selected_cells.csv')==rec['selected_sha256'];x=sparse.load_npz(folder/'selected_counts.npz');meta=pd.read_csv(folder/'selected_cells.csv');assert b.sha256(folder/'gene_names.csv')==rec['genes_sha256'];measured=pd.read_csv(folder/'gene_names.csv').measured.to_numpy(bool)
  # The extension reuses the exact original NTC selection. Prove counts and identities agree.
  ox,om,_,_=raw_inputs(line,g);op=om[om.gene_target=='Non-Targeting'].copy();np0=meta[meta.gene_target=='Non-Targeting'].copy();look=op.reset_index().set_index('cell_id')['index'];ix=look.loc[np0.cell_id].to_numpy();assert len(np0)==len(op);assert (x[np0.index.to_numpy()]-ox[ix]).nnz==0;del ox,om
  ss=summarize_raw(x,meta,g,new,measured,2)
  for ni in np.flatnonzero(ss['n']>=30):append(ci,line,new[ni],ss['truth'][ni],ss['errors'][ni],ss['q'][ni],ss['valid'][ni],ss['knockdown'][ni],ss['n'][ni])
  del x;gc.collect();print('Final source',line,'summarized',flush=True)
 # H1: retain original training/validation/test source response convention, no repeated controls.
 hrec=json.loads((A/'H1_new_source.json').read_text());x0,meta0,measured,_=raw_inputs('H1',g);take0=meta0.gene_target.eq('Non-Targeting').to_numpy();matrices=[x0[take0]];metas=[meta0[take0].copy()];del x0;gc.collect()
 for rec in hrec['splits']:
  folder=O/'H1'/rec['split'];assert b.sha256(folder/'counts.npz')==rec['sha256'];assert b.sha256(folder/'cells.csv')==rec['cells_sha256'];assert b.sha256(folder/'genes.csv')==rec['genes_sha256'];sg=pd.read_csv(folder/'genes.csv').gene_name.to_numpy(str);lookup=np.array([gp.get(v,-1) for v in sg]);raw=sparse.load_npz(folder/'counts.npz').tocoo();ok=lookup[raw.col]>=0;xx=sparse.csr_matrix((raw.data[ok],(raw.row[ok],lookup[raw.col[ok]])),shape=(raw.shape[0],len(g)));mm=pd.read_csv(folder/'cells.csv');mm['gene_target']=mm.target_gene;mm['sample']=mm.batch.astype(str);mm['cell_id']=mm.source_cell_id;mm['total_counts']=xx.sum(1).A1;mm['partition']=0;matrices.append(xx);metas.append(mm)
 xx=sparse.vstack(matrices,format='csr');mm=pd.concat(metas,ignore_index=True);assert mm.cell_id.is_unique;ss=summarize_raw(xx,mm,g,new,measured,2)
 for ni in np.flatnonzero(ss['n']>=30):append(4,'H1',new[ni],ss['truth'][ni],ss['errors'][ni],ss['q'][ni],ss['valid'][ni],ss['knockdown'][ni],ss['n'][ni])
 result={k:np.concatenate(v,axis=0) for k,v in parts.items()};result.update(genes=g.astype(str),targets=t.astype(str),profiles=old['profiles'])
 for k in ['truth','errors','valid','q','knockdown','row_contexts','row_families']:np.testing.assert_array_equal(result[k][:nold],old[k][take])
 assert not pd.DataFrame({'c':result['row_contexts'],'t':result['row_targets']}).duplicated().any();np.savez_compressed(A/'source_cache.npz',**result);pd.DataFrame(provenance).to_csv(A/'new_source_rows.csv',index=False)
 coverage=[]
 for ti,target in enumerate(t):
  rr=np.flatnonzero(result['row_targets']==ti);families=sorted(set(result['row_families'][rr].tolist()));coverage.append(dict(target_gene=str(target),reused_validation_target=bool(target in set(old['targets'])),source_families='+'.join(families),n_source_families=len(families),source_rows=len(rr),any_output_measured=bool(result['valid'][rr].any()) if len(rr) else False))
 pd.DataFrame(coverage).to_csv(A/'target_coverage.csv',index=False);b.write_json(A/'source_provenance.json',dict(cache_sha256=b.sha256(A/'source_cache.npz'),original_source_sha256=b.sha256(R/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),reused_rows=nold,reused_rows_bitwise_unchanged=True,new_rows=len(result['truth'])-nold,script_sha256=b.sha256(Path(__file__)),rows_by_family={str(f):int((result['row_families']==f).sum()) for f in np.unique(result['row_families'])},unsupported_targets=[x['target_gene'] for x in coverage if not x['any_output_measured']],complete_treatment_truth_for_final=False));print('Final measured response cache complete',len(result['truth']),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
