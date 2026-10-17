"""Local estimator adapter: discard training-constant/all-missing columns.

No query or validation outcomes select columns. The scientific head recipes
and all finite, variable predictor values remain unchanged.
"""
import json
from pathlib import Path
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier as C, HistGradientBoostingRegressor as R
import v38_response_hurdle as m
class ProjectedHead:
 def __init__(self,kind,**params):self.kind=kind;self.params=params
 def fit(self,x,y,sample_weight=None):
  lo=np.min(np.where(np.isfinite(x),x,np.inf),axis=0);hi=np.max(np.where(np.isfinite(x),x,-np.inf),axis=0)
  self.columns=np.flatnonzero(np.isfinite(lo)&np.isfinite(hi)&(hi>lo));assert len(self.columns)
  self.estimator=(C if self.kind=='classifier' else R)(**self.params).fit(x[:,self.columns],y,sample_weight=sample_weight)
  return self
 def predict(self,x):return self.estimator.predict(x[:,self.columns])
 def predict_proba(self,x):return self.estimator.predict_proba(x[:,self.columns])
def classifier(**params):return ProjectedHead('classifier',**params)
def regressor(**params):return ProjectedHead('regressor',**params)
def install():
 m.HistGradientBoostingClassifier=classifier;m.HistGradientBoostingRegressor=regressor
def main():
 import argparse
 from threadpoolctl import threadpool_limits
 from run_v38_response_hurdle import check
 p=argparse.ArgumentParser();p.add_argument('--outer',choices=list(m.OUTERS)+['all']);a=p.parse_args()
 check();receipt=m.A/'binning_amendment.json'
 if not receipt.exists():
  m.write(receipt,{'parent_plan_sha256':m.b.sha256(m.A/'plan.json'),'adapter_sha256':m.b.sha256(Path(__file__)),'initial_failure':'sklearn binning _find_binning_thresholds -> NumPy sliding_window_view: window shape cannot be larger than input array shape. Initial H1 training table built; no model fitted/saved or validation outcome inspected.','change':'Fit each head on columns with >1 distinct finite TRAINING value. All-missing or constant features discarded using training rows only. Missing variable features remain NaN. No scientific threshold, seed, data split, outcome or candidate changed.','known_outcomes':'Only contract tests and table sizes; no V38 fitted predictions or scores.'})
 assert m.read(receipt)['adapter_sha256']==m.b.sha256(Path(__file__))
 install()
 with threadpool_limits(limits=2):
  for outer in [a.outer] if a.outer else m.OUTERS:m.fit(outer)
if __name__=='__main__':main()

