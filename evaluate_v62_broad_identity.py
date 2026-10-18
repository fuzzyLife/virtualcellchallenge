"""Broader genuine-count PDS diagnostic; fixed candidate, no model reselection."""
import gc
import numpy as np
from threadpoolctl import threadpool_limits
from scipy import sparse
import v62_conserved_response as m
from v62_pool_cache import install
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder
from v25_methods import raw
from audit_transfer_validation import profile,pds
install()

def main():
 m.check();candidate=read(m.A/'production_plan.json')['candidate'];d=m.data();genes=d['genes'];targets=d['targets'];records=[]
 plan=m.A/'broad_identity_plan.json'
 if not plan.exists():write(plan,dict(candidate=candidate,auditor_sha256=sha(__file__),cohort='Every target with at least 30 actual cells in the existing aligned source extraction, no score-based filtering; same V37 inference/control pools; all 300 panel genes removed from cosine features.',purpose='Broader mean-profile identity diagnostic, not another candidate-selection set. Measured population depths vary and these public studies are already reused. Predictions are desired pooled means before cell rounding, not final packaged counts.'))
 assert sha(__file__)==read(plan)['auditor_sha256']
 for ctx in contexts():
  name=ctx['family'];dest=m.A/'broad_identity'/(name+'.json');dest.parent.mkdir(exist_ok=True)
  if dest.exists():continue
  folder,_=checked_folder(name)
  with np.load(folder/'effects.npz',allow_pickle=False) as z:q=z['query'].copy();base=z['v24'].copy();bulk=z['bulk'].copy();measured=z['measured'].copy()
  effects=m.predict(q,base,(key(name),))[candidate];x,meta,mask,source=raw(name,genes);labels=meta.gene_target.to_numpy(str)
  selected=[(ti,t,np.flatnonzero(labels==t)) for ti,t in enumerate(targets) if np.count_nonzero(labels==t)>=30];axis=measured&mask&~np.isin(genes,targets);truth=np.asarray([profile(x[ids]) for ti,t,ids in selected]);ref=profile(sparse.load_npz(folder/'reference.npz'));pr={}
  for arm,eff in [('v24',base),(candidate,effects)]:
   pp=[]
   for ti,t,ids in selected:
    factor=np.exp2(eff[ti]);factor[list(genes).index(t)]=.1;prob=bulk*factor;prob/=prob.sum();pp.append(np.log1p(50000*prob))
   pred=np.asarray(pp);pr[arm]=dict(pds=pds(pred[:,axis],truth[:,axis],ref[axis]),profile_MSE=float(np.mean((pred[:,axis]-truth[:,axis])**2)))
  assert abs(pds(np.tile(ref[axis],(len(selected),1)),truth[:,axis],ref[axis])-.5)<1e-12
  r=dict(family=name,study=key(name),targets=len(selected),feature_genes=int(axis.sum()),target_names=[str(t) for ti,t,ids in selected],actual_cell_counts=[len(ids) for ti,t,ids in selected],metrics=pr,source_record=source,plan_sha256=sha(plan));write(dest,r);records.append(r);print('V62 broad identity',name,len(selected),pr,flush=True);del x;gc.collect()
 records=[read(p) for p in sorted((m.A/'broad_identity').glob('*.json'))];write(m.A/'broad_identity_results.json',dict(records=records,scope='Genuine measured pooled-count truth; intended prediction means, before final cell rounding. Public exploratory diagnostic only. No candidate reselection.'))

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
