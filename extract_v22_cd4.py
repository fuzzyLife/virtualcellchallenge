"""Inspect/extract checksum-verified donor CD4 raw counts; freeze reserve by metadata."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
from v22_methods import ROOT,A,write
import build_submission as b

def main():
 source=ROOT/'reference/v22_cd4/Marson2025_D1_Rest_filtered.h5ad';dest=A/'cd4_selected';dest.mkdir(exist_ok=True)
 if (dest/'extraction.json').exists():print('Existing frozen CD4 extraction retained');return
 record=json.loads((A/'cd4_download.json').read_text());assert record['md5_verified'] and source.stat().st_size==record['bytes']
 manifest,genes,targets=b.read_spec(ROOT)
 with h5py.File(source,'r') as h:print('CD4 layers',list(h.get('layers',{})),'raw',list(h.get('raw',{})),'X',dict(h['X'].attrs),flush=True)
 obj=ad.read_h5ad(source,backed='r');meta=obj.obs.copy();var=obj.var.copy();print('CD4 obs',meta.columns.tolist(),'var',var.columns.tolist(),'shape',obj.shape,flush=True)
 # Names and raw-layer choice are explicitly inspected, never guessed from dense values.
 inspection={'obs_columns':meta.columns.tolist(),'var_columns':var.columns.tolist(),'shape':obj.shape,'obs_head':meta.head(3).astype(str).to_dict(orient='records'),'var_head':var.head(3).astype(str).to_dict(orient='records'),'layers':list(obj.layers.keys())}
 write(A/'cd4_inspection.json',inspection)
 label=next((v for v in ['gene_target','target_gene','gene','perturbation','target'] if v in meta),None)
 if label is None:raise ValueError('Inspect CD4 labels before explicit mapping; no supported target column')
 names=meta[label].astype(str).to_numpy();print('CD4 labels',pd.Series(names).value_counts().head(15).to_dict(),flush=True)
 gene_col=next((v for v in ['gene_name','gene_symbol','feature_name','symbol','gene'] if v in var),None)
 symbols=var[gene_col].astype(str).to_numpy() if gene_col else var.index.to_numpy(str)
 if len(set(symbols)&set(genes))<10000:raise ValueError('CD4 requires audited gene-ID mapping')
 controls=np.isin(names,['Non-Targeting','non-targeting','NTC','control','Control','non_targeting','NT'])
 if controls.sum()<100:raise ValueError('CD4 control annotation requires explicit audit')
 take=np.flatnonzero(controls|np.isin(names,targets));layer='counts' if 'counts' in obj.layers else 'X'
 if layer!='X':
  with h5py.File(source,'r') as h:
   from anndata.io import sparse_dataset
   matrix=sparse_dataset(h['layers'][layer])[take]
 else:matrix=obj.X[take].tocsr()
 obj.file.close();matrix=sparse.csr_matrix(matrix)
 if not np.isfinite(matrix.data).all() or (matrix.data<0).any() or not np.array_equal(matrix.data,np.round(matrix.data)):raise ValueError('CD4 selected matrix is not raw integer counts')
 lookup={g:i for i,g in enumerate(genes)};cols=np.array([lookup.get(g,-1) for g in symbols]);coo=matrix.tocoo();keep=cols[coo.col]>=0;x=sparse.csr_matrix((coo.data[keep].astype(np.int32),(coo.row[keep],cols[coo.col[keep]])),shape=(len(take),len(genes)));x.sum_duplicates();x.sort_indices()
 out=meta.iloc[take].copy();out['gene_target']=names[take];out.loc[controls[take],'gene_target']='Non-Targeting';out['cell_id']=meta.index.to_numpy(str)[take]
 batch=next((v for v in ['sample','batch','batch_id','donor'] if v in out),None);out['sample']=out[batch].astype(str) if batch else 'D1_Rest';out['source_row']=take;out['total_counts']=np.asarray(x.sum(1)).ravel();out=out.reset_index(drop=True)
 assert out.cell_id.is_unique and (out.total_counts>0).all()
 sparse.save_npz(dest/'counts.npz',x);out.to_csv(dest/'cells.csv',index=False);pd.DataFrame({'gene_name':genes,'measured':np.isin(genes,symbols)}).to_csv(dest/'genes.csv',index=False)
 counts=out.gene_target.value_counts();counts.to_csv(A/'cd4_target_counts.csv');write(dest/'extraction.json',{'source':str(source),'source_md5':record['file']['checksum'],'counts_sha256':b.sha256(dest/'counts.npz'),'cells_sha256':b.sha256(dest/'cells.csv'),'genes_sha256':b.sha256(dest/'genes.csv'),'raw_layer':layer,'source_rows':len(take),'source_genes':len(symbols),'measured_competition_genes':int(np.isin(genes,symbols).sum()),'label_column':label,'symbol_column':gene_col,'study':'Marson2025 D1 Rest; existing CD4 study family, new raw-cell evaluation source','script_sha256':b.sha256(Path(__file__))})
 # Selection by cell identities/counts only, before any treatment expression metrics.
 left=[];right=[];good=[]
 for bi,batch in enumerate(sorted(out['sample'].unique())):
  ids=np.flatnonzero((out.gene_target=='Non-Targeting')&(out['sample']==batch));ids=np.random.default_rng(22300+bi).permutation(ids)
  if len(ids)>=40:left.extend(ids[::2].tolist());right.extend(ids[1::2].tolist());good.append(batch)
 eligible=[str(t) for t in sorted(targets) if int(((out.gene_target==t)&out['sample'].isin(good)).sum())>=400][:8]
 if len(eligible)<4:raise ValueError('CD4 insufficient powered reserve; inspect counts, preserve extraction and retain V16 unless a reduced-power amendment is frozen')
 cases=[]
 for ti,t in enumerate(eligible):
  ids=np.flatnonzero((out.gene_target==t)&out['sample'].isin(good));rng=np.random.default_rng(22400+ti);actual=rng.choice(ids,400,replace=False);template=rng.choice(left,400,replace=False)
  cases.append({'target':t,'seed':22400,'actual':actual.tolist(),'template':template.tolist(),'split':'confirmation'})
 write(A/'reserve_plan.json',{'family':'CD4','left':left,'right':right,'cells':cases,'selected_by':'Alphabetical target names, counts>=400, matched-control batch support only; one seed and eight targets maximum','extraction_sha256':b.sha256(dest/'extraction.json'),'parent_plan_sha256':b.sha256(A/'plan.json'),'treatment_metrics_opened':False});print('CD4 raw reserve frozen',eligible,len(left),len(right),flush=True)
if __name__=='__main__':main()
