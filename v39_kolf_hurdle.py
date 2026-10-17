"""V39: a paired data-source ablation, adding KOLF z summaries to V38.

The external study is available in every fold. Z scores remain auxiliary
predictors; there is no conversion to counts, fold changes, or DE labels.
"""
from functools import lru_cache
import numpy as np
import v38_response_hurdle as base
ORIGINAL=base.features
ROOT=base.ROOT
OUT=ROOT/'experiments/v39_kolf_support'
A=OUT/'analysis'
@lru_cache(1)
def kolf():
 rec=base.read(A/'collection.json');p=A/'response_features.npz';assert base.b.sha256(p)==rec['arrays_sha256']
 with np.load(p,allow_pickle=False) as z:
  np.testing.assert_array_equal(z['genes'],base.data()['genes']);np.testing.assert_array_equal(z['targets'],base.data()['targets'])
  return z['z'].copy(),z['saturated'].copy()
def features(desc,ti,positions):
 x,n=ORIGINAL(desc,ti,positions);z,sat=kolf();v=z[ti,positions];prior=desc['pooled'][ti,positions]
 extra=np.column_stack([v,abs(v),np.clip(v,-1,1),v*prior,np.where(np.isfinite(v),(v*prior>0).astype(float),np.nan),np.isfinite(v).astype(float),np.where(np.isfinite(v),sat[ti,positions].astype(float),np.nan)]).astype(np.float32)
 return np.column_stack([x,extra]),n
def install():
 base.OUT=OUT;base.A=A;base.features=features

