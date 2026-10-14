"""Independently recompute V19 control covariance and selected source treatment shifts."""
from pathlib import Path
import json,gc
import numpy as np
import build_submission as b
from v19_methods import ROOT,A,FAMILIES,load_plan,load_model
from fit_v7_response import raw_inputs,read_cache

def main():
 plan=load_plan();data=read_cache();assert b.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==plan['source_cache_sha256'];contexts={c['family']:c for c in plan['contexts']};records=[]
 for held in FAMILIES+['all']:
  model=load_model(held);positions=model['positions'];projection=model['basis']/model['scale'][:,None]
  for si,family in enumerate(model['training_families']):
   assert family!=held;x,meta,measured,_=raw_inputs(str(family),data['genes']);ids=np.array(contexts[str(family)]['left']);depth=meta.total_counts.to_numpy(float);batches=meta['sample'].to_numpy(str);labels=meta.gene_target.to_numpy(str);assert np.all(labels[ids]=='Non-Targeting') and measured[positions].all()
   values=np.log1p(x[ids][:,positions].toarray().astype(float)*1e4/depth[ids,None]);z=np.log1p(depth[ids]);names=sorted(set(batches[ids]));bi=np.array([names.index(v) for v in batches[ids]]);means=np.array([values[bi==i].mean(0) for i in range(len(names))]);centers=np.array([z[bi==i].mean() for i in range(len(names))]);zd=z-centers[bi];beta=(zd@values)/max(float(zd@zd),len(ids)*1e-12);residual=values-means[bi]-zd[:,None]*beta;score=residual@projection
   cov=score.T@score/len(score)-np.outer(score.mean(0),score.mean(0));regularized=.8*cov+.2*np.diag(np.diag(cov))+.05*max(float(np.diag(cov).mean()),1e-8)*np.eye(16);np.testing.assert_allclose(regularized,model['source_covariances'][si],rtol=1e-8,atol=1e-8)
   checked=[]
   for target in ['ABCD1','HDAC8','VIM','SMARCA5']:
    ti=int(np.flatnonzero(data['targets']==target)[0]);rows=np.flatnonzero((labels==target)&np.isin(batches,names));assert len(rows)==model['source_sizes'][si,ti]
    if len(rows)<40:continue
    v=np.log1p(x[rows][:,positions].toarray().astype(float)*1e4/depth[rows,None]);index=np.array([names.index(v) for v in batches[rows]]);adjusted=v-means[index]-(np.log1p(depth[rows])-centers[index])[:,None]*beta;expected=(adjusted@projection).mean(0);np.testing.assert_allclose(expected,model['source_shifts'][si,ti],rtol=1e-8,atol=1e-8);checked.append(target)
   records.append({'held':held,'source':str(family),'control_covariance_recomputed':True,'source_shifts_recomputed':checked});del x,values;gc.collect();print('V19 independent source math',held,family,flush=True)
 b.write_json(A/'model_input_independent_audit.json',{'passed':True,'records':records,'plan_sha256':b.sha256(A/'plan.json'),'auditor_sha256':b.sha256(Path(__file__)),'frozen_source_cache_sha256_verified':True,'method':'Dense selected-gene log1p normalization, direct within-batch depth least squares, direct score cross-products and selected treatment means; distinct from fitted sparse pipeline.'})
if __name__=='__main__':
 from threadpoolctl import threadpool_limits
 with threadpool_limits(limits=2):main()
