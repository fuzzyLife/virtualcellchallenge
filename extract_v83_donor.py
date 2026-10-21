"""Extract and freeze V83 donor cases by identities/counts only, never effect outcomes."""
from pathlib import Path
import gc,json,hashlib
import h5py,numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from v57_pie import read,write,sha
import v83_context_expansion as m

def main():
 source=m.ROOT/'reference/v25_sources/Marson2025_D3_Rest_filtered.h5ad';dest=m.A/'donor';dest.mkdir(exist_ok=True);assert not (dest/'extraction.json').exists();receipt=read(m.ROOT/'experiments/v25_cross_domain/analysis/holdout_download.json');assert receipt['md5_verified'] and source.stat().st_size==receipt['bytes']
 d=m.data();genes=d['genes'];plan=read(m.A/'plan.json');targets=plan['fresh_donor']['targets'];obj=ad.read_h5ad(source,backed='r');obs=obj.obs.copy();symbols=obj.var_names.to_numpy(str);assert 'perturbation' in obs and 'lane_id' in obs;labels=obs.perturbation.astype(str).to_numpy();control=labels=='control';assert control.sum()==10000
 original_indices=np.flatnonzero(control|np.isin(labels,targets));selected=obs.iloc[original_indices].copy();selected['gene_target']=labels[original_indices];selected.loc[control[original_indices],'gene_target']='Non-Targeting';selected['cell_id']=obs.index.to_numpy(str)[original_indices];selected['source_row']=original_indices;selected['sample']=selected.lane_id.astype(str);selected=selected.reset_index(drop=True)
 left=[];right=[];good=[]
 for bi,lane in enumerate(sorted(selected['sample'].unique())):
  ix=np.flatnonzero((selected.gene_target=='Non-Targeting')&(selected['sample']==lane));ix=np.random.default_rng(83000+bi).permutation(ix)
  if len(ix)>=40:left.extend(ix[::2].tolist());right.extend(ix[1::2].tolist());good.append(lane)
 left=np.array(left);right=np.array(right);assert len(left)>=400 and len(right)>=400;cases=[];exclusions=[]
 for ti,target in enumerate(targets):
  actual=np.flatnonzero((selected.gene_target==target)&selected['sample'].isin(good));n=len(actual)
  if n<400:exclusions.append(dict(target=target,matched_cells=n,reason='Less than400 cells in lanes with disjoint NTC support'));continue
  views=[]
  for offset in [0,620000]:
   seed=83083+ti*1000+offset;rng=np.random.default_rng(seed);truth=rng.choice(actual,400,replace=False);n_per_lane=selected.iloc[truth]['sample'].value_counts();templates=[]
   for lane,nlane in n_per_lane.items():
    pool=left[selected.iloc[left]['sample'].to_numpy()==lane]
    if len(pool)<int(nlane):raise ValueError('Cannot sample400 distinct lane-matched controls; freeze a metadata-only amendment')
    templates.extend(rng.choice(pool,int(nlane),replace=False).tolist())
   template=np.array(templates);rng.shuffle(template);assert len(set(template))==400 and not set(left)&set(right)
   views.append(dict(offset=offset,seed=seed,actual=truth.tolist(),template=template.tolist()))
  cases.append(dict(target=target,available_cells=n,views=views))
 cohort=dict(cases=cases,metadata_exclusions=exclusions,left=left.tolist(),right=right.tolist(),source_md5=receipt['file']['checksum'],plan_sha256=sha(m.A/'plan.json'),selected_before_expression_metrics=True,script_sha256=sha(__file__),all_query_treatments_excluded_from_fit=True)
 write(dest/'cohort.json',cohort);print('V83 donor metadata frozen',len(cases),'targets',len(left),'fit controls',len(right),'reference',flush=True)
 # Only after the cohort is fixed, read expression. No preprocessing fit on treatments.
 x=sparse.csr_matrix(obj.X[original_indices]);obj.file.close();assert np.isfinite(x.data).all() and (x.data>=0).all() and np.array_equal(x.data,np.round(x.data)),'Donor X is not raw integer counts; inspect explicit layers instead'
 gp={str(g):i for i,g in enumerate(genes)};cols=np.array([gp.get(str(g),-1) for g in symbols]);coo=x.tocoo();keep=cols[coo.col]>=0;aligned=sparse.csr_matrix((coo.data[keep].astype(np.int32),(coo.row[keep],cols[coo.col[keep]])),shape=(len(selected),len(genes)));aligned.sum_duplicates();aligned.sort_indices();depth=np.asarray(aligned.sum(1)).ravel();assert (depth>0).all() and selected.cell_id.is_unique
 sparse.save_npz(dest/'counts.npz',aligned);selected['total_counts']=depth;selected.to_csv(dest/'cells.csv',index=False);pd.DataFrame(dict(gene_name=genes,measured=np.isin(genes,symbols))).to_csv(dest/'genes.csv',index=False)
 write(dest/'extraction.json',dict(source_file=str(source),source_md5=receipt['file']['checksum'],counts_sha256=sha(dest/'counts.npz'),cells_sha256=sha(dest/'cells.csv'),genes_sha256=sha(dest/'genes.csv'),cohort_sha256=sha(dest/'cohort.json'),raw_layer='X verified finite nonnegative integer counts',source_cells=len(selected),measured_competition_genes=int(np.isin(genes,symbols).sum()),script_sha256=sha(__file__),target_outcomes_scored=False,whole_CD4_family_excluded=True))
 print('V83 donor integer counts extracted',aligned.shape,aligned.nnz,flush=True)
if __name__=='__main__':main()
