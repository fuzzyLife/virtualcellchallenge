"""Describe final basal similarities using pairwise finite-value Spearman correlations."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from scipy.stats import spearmanr
import build_submission as b
from fit_v7_response import read_cache
R=Path(__file__).resolve().parent;A=R/'experiments/final_v91/analysis'
def main():
 d=read_cache();records=[];contexts=['K562','CD4_Rest','CD4_Stim8hr','CD4_Stim48hr','H1','HCT116','HEK293T']
 queries={}
 for c in ['D','E','F']:
  with np.load(A/(c+'_controls.npz'),allow_pickle=False) as z:queries[c]=z['query'].copy()
 for c,q in queries.items():
  for name,p in [*zip(contexts,d['profiles']),*queries.items()]:
   ok=np.isfinite(q)&np.isfinite(p);rho=float(spearmanr(q[ok],p[ok]).statistic);records.append(dict(final_context=c,reference=name,finite_gene_pairs=int(ok.sum()),spearman=rho))
 pd.DataFrame(records).to_csv(A/'basal_rank_correlations.csv',index=False);b.write_json(A/'context_diagnostic.json',dict(pairwise_finite_values_only=True,pseudocount_added=False,cell_identity_claim=False,scope='Basal compatibility is descriptive and does not validate perturbation transfer or identify an anonymous cell line.',nearest_external={c:max([v for v in records if v['final_context']==c and v['reference'] in contexts],key=lambda v:v['spearman']) for c in queries}))
if __name__=='__main__':main()
