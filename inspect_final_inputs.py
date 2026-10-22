"""Audit the official final bundle and local response coverage without relabelling validation arrays."""
from pathlib import Path
import json, sys, importlib.metadata
import numpy as np,pandas as pd,anndata as ad,pyarrow.parquet as pq,h5py
import build_submission as b
from v6_methods import moments
R=Path(__file__).resolve().parent;O=R/'experiments/final_v91';A=O/'analysis'
def main():
 A.mkdir(parents=True,exist_ok=True);m,g,t=b.read_spec(R/'final');vm,vg,vt=b.read_spec(R)
 assert m['partition']=='test' and m['contexts']==['D','E','F'] and m['cells_per_pert']==400
 np.testing.assert_array_equal(g,vg);records=[]
 for c in m['contexts']:
  x,guide,libs,_,summary,_=b.load_control(R/'final',c,g,m);depth,p,bulk,var=moments(x)
  np.savez_compressed(A/(c+'_controls.npz'),query=p*1e6,bulk=bulk,library=libs)
  summary['sha256']=b.sha256(R/'final'/f'context_{c}.h5ad');summary['query_sha256']=b.sha256(A/(c+'_controls.npz'));records.append(summary);print('Final controls verified',c,x.shape,x.nnz,flush=True)
 cover=[]
 for line in ['HCT116','HEK293T']:
  path=R/'reference/orion_pilot'/(line+'_Batch1.parquet');tab=pq.read_table(path,columns=['gene_target','pass_guide_filter','total_counts']).to_pandas();n=tab[tab.pass_guide_filter.astype(bool)&(tab.total_counts>0)].groupby('gene_target').size()
  for target in t:cover.append(dict(source=line+'_local_batch1',target=str(target),cells=int(n.get(target,0)),eligible_ge30=bool(n.get(target,0)>=30)))
 with h5py.File(R/'reference/K562_gwps_raw_bulk_01.h5ad') as f:
  labels=f['obs/gene_transcript'].asstr()[:];n=f['obs/num_cells_filtered'][:];fold=f['obs/fold_expr'][:];names=np.array([s.split('_',2)[1] if len(s.split('_',2))==3 and s.split('_',1)[0].isdigit() else '' for s in labels])
  for target in t:
   take=(names==target)&(n>=30)&np.isfinite(fold)&(fold<=.5);cover.append(dict(source='K562_full_pseudobulk',target=str(target),cells=int(n[take].sum()),eligible_ge30=bool(take.any())))
 frame=pd.DataFrame(cover);frame.to_csv(A/'local_target_coverage.csv',index=False)
 runtime=dict(python=sys.version,executable=sys.executable,packages={p:importlib.metadata.version(p) for p in ['numpy','scipy','pandas','h5py','anndata','threadpoolctl','pyarrow','vcc-cli']})
 report=dict(partition=m['partition'],panel=m['panel_id'],contexts=m['contexts'],gene_axis_matches_validation=True,n_targets=len(t),overlap_validation_targets=len(set(t)&set(vt)),new_targets=sorted(set(t)-set(vt)),inputs={str(p.relative_to(R)):dict(bytes=p.stat().st_size,sha256=b.sha256(p)) for p in sorted((R/'final').iterdir()) if p.is_file()},control_summaries=records,local_coverage=frame.groupby('source').eligible_ge30.sum().astype(int).to_dict(),runtime=runtime,no_treatment_truth_available=True,automatic_upload=False)
 b.write_json(A/'input_audit.json',report);print(json.dumps({k:report[k] for k in ['contexts','n_targets','overlap_validation_targets','local_coverage']},indent=2))
if __name__=='__main__':main()
