"""Numerically identical normalized-cell caching for the V36 sampler.

The original frozen sampling/training code remains unchanged. Cache contains
only a deterministic per-cell normalization, not group labels or outcomes.
"""
import argparse,gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v35_methods as m
import v36_methods as p
import run_v35_v36 as runner
original=p.view_stats
cache={}

def cached(x,ids,axis):
 key=tuple(axis)
 if cache.get('x') is not x or cache.get('axis')!=key:
  cache.clear();libs=np.asarray(x.sum(axis=1,dtype=np.float64)).ravel();assert (libs>0).all();y=x[:,axis].toarray().astype(np.float64);y*=1e6;y/=libs[:,None];cache.update(x=x,axis=key,values=y)
 y=cache['values'][ids]
 return y.mean(0),np.maximum(y.var(0),0)

def verify():
 from scipy import sparse
 rng=np.random.default_rng(36);x=sparse.csr_matrix(rng.poisson(2,(51,87)));axis=np.array([0,4,8,20,35,57,85]);ids=np.array([3,1,9,30,24,48]);a=original(x,ids,axis);b=cached(x,ids,axis)
 for u,v in zip(a,b):np.testing.assert_array_equal(u,v)
 data=m.old.old.read_cache();x,meta,measured,rec=m.raw('HEK293T',data['genes']);x=x[:103].copy();axis=np.flatnonzero(measured)[:71];ids=rng.permutation(x.shape[0])[:32];a=original(x,ids,axis);b=cached(x,ids,axis)
 for u,v in zip(a,b):np.testing.assert_array_equal(u,v)
 cache.clear();return {'toy_bit_exact':True,'actual_raw_subset_bit_exact':True,'raw_source':'HEK293T','raw_record':rec,'rows':len(ids),'genes':len(axis)}

def install():
 p.view_stats=cached
 old=p.build_populations
 def build(*args,**kwargs):
  try:return old(*args,**kwargs)
  finally:cache.clear();gc.collect()
 p.build_populations=build

def main():
 parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['verify','train36']);parser.add_argument('--outer',choices=m.OUTERS+['all']);args=parser.parse_args();m.setup()
 if args.stage=='verify':
  tests=verify();m.write(m.P36/'normalization_cache_amendment.json',{'parent_plan_sha256':m.sha(m.A/'plan.json'),'original_methods_sha256':m.sha(m.ROOT/'v36_methods.py'),'cache_code_sha256':m.sha(m.ROOT/'run_v35_v36_cached.py'),'equivalence_checks':tests,'reason':'Before any V36 model fit/outcome, replace repeated CSR row/column normalization with one deterministic float64 cell-by-gene normalized matrix per source. Same sample rows, RNG order, means/variances, fitting targets, objectives and epochs. Original source-preparation log retained.','V36_outcomes_available':False});print('V36 cache checks passed',flush=True)
 else:
  receipt=json.loads((m.P36/'normalization_cache_amendment.json').read_text());assert receipt['cache_code_sha256']==m.sha(m.ROOT/'run_v35_v36_cached.py');install()
  for outer in [args.outer] if args.outer else m.OUTERS:p.train(outer)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()