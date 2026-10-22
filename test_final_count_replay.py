"""Check the unchanged count emitter against a full archived V90 population per context."""
from pathlib import Path
import gc
import numpy as np,h5py,pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as b
from v6_methods import moments
from v16_methods import emit
from evaluate_atlas_transfer import load_atlas_model
R=Path(__file__).resolve().parent;A=R/'experiments/final_v91/analysis'
def main():
 m,g,t=b.read_spec(R);atlas,_=load_atlas_model();records=[]
 with np.load(R/'experiments/v90_interval_diagnostic/analysis/reference_model_arrays.npz',allow_pickle=False) as z:effects=z['effects'].copy()
 with h5py.File(R/'experiments/v90_interval_diagnostic/prediction.h5ad','r') as old:
  ptr=old['X/indptr'][:]
  for ci,c in enumerate(m['contexts']):
   x,guides,library,*_=b.load_control(R,c,g,m);depth,cp,bulk,_=moments(x);ti=0;rows=b.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([20261011,ci,ti])));own=list(g).index(t[ti]);new,diag=emit(x,rows,effects[ci,ti],own,20261011+ci*100000+ti,'residual10',(depth,cp,bulk),atlas);start=(ci*300+ti)*400;lo,hi=int(ptr[start]),int(ptr[start+400]);oldcount=sparse.csr_matrix((old['X/data'][lo:hi],old['X/indices'][lo:hi],ptr[start:start+401]-lo),shape=(400,len(g)));diff=new-oldcount;l1=float(abs(diff).sum()/2/library[rows].sum());np.testing.assert_array_equal(new.sum(1).A1,library[rows]);assert l1<=1e-5,'Count runtime changes exceed the frozen diagnostic tolerance';records.append(dict(context=c,target=str(t[ti]),cells=400,bitwise_count_values_equal=bool(diff.nnz==0),runtime_changed_transcript_fraction=l1,maximum_absolute_count_change=int(abs(diff.data).max(initial=0)),**diag));print('Archived count replay',c,l1,flush=True);del x;gc.collect()
 b.write_json(A/'count_runtime_replay.json',dict(passed=True,groups=records,tolerance_changed_transcript_fraction=1e-5,scope='Three full archived populations verify numerical/count compatibility; technical replay is not biological accuracy validation.',script_sha256=b.sha256(Path(__file__))))
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
