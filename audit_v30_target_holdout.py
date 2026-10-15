"""Replay the unseen-target stress models and verify fitting/test target separation."""
import ast,json
import numpy as np
from threadpoolctl import threadpool_limits
import v30_methods as m
import build_submission as b
D=m.A/'target_holdout'

def main():
 p=json.loads((D/'plan.json').read_text());r=json.loads((D/'results.json').read_text());assert r['plan_sha256']==b.sha256(D/'plan.json');keys,e=m.embeddings();lookup={g:i for i,g in enumerate(keys)};checks=[]
 for rec in r['records']:
  name=rec['source'];variant=rec['variant'];spec=p['sources'][name];source=m.A/'sources'/f'{name}.npz';assert b.sha256(source)==spec['source_sha256'];d=np.load(source,allow_pickle=False);file=D/f'{name}_{variant}.npz';assert b.sha256(file)==rec['model_sha256'];model=np.load(file,allow_pickle=False);train=np.array(spec['train']+spec['inner']);test=np.array(spec['test']);assert not set(train)&set(test);np.testing.assert_array_equal(d['targets'][train],rec['training_targets']);np.testing.assert_array_equal(d['targets'][test],rec['test_target_names']);choice=min(rec['inner_results'],key=lambda k:(rec['inner_results'][k]['mse_vs_constant'],ast.literal_eval(k)));assert list(ast.literal_eval(choice))==rec['config'];emb=e if variant=='string' else e[np.random.default_rng(30301).permutation(len(e))];xx=emb[[lookup[g] for g in d['targets'][train]]];np.testing.assert_allclose(xx.mean(0),model['feature_mean'],atol=1e-14)
  y=d['y'][train].astype(float);w=d['n'][train].astype(float);w=w/(w+100);w/=w.mean();gp={g:i for i,g in enumerate(m.read_cache()['genes'][d['positions']])};mask=np.ones(y.shape,bool)
  for i,g in enumerate(d['targets'][train]):
   if g in gp:mask[i,gp[g]]=False
  y[~mask]=0;np.testing.assert_allclose((w[:,None]*y).sum(0)/(w[:,None]*mask).sum(0),model['bias'],atol=1e-12);rank=rec['config'][0];err=const=ma=mb=zero=0.;coordinates=0
  for lo in range(0,len(test),128):
   ids=test[lo:lo+128];names=d['targets'][ids];qx=(emb[[lookup[t] for t in names]]-model['feature_mean'])/model['feature_scale'];pred=np.clip((qx@model['weights']@model['basis'][:,:rank].T)*model['scale']+model['bias'],-3,3);truth=d['y'][ids];valid=np.ones(truth.shape,bool)
   for i,g in enumerate(names):
    if g in gp:valid[i,gp[g]]=False
   delta=(pred-truth)[valid];dd=(model['bias'][None,:]-truth)[valid];tt=truth[valid];err+=float(delta@delta);const+=float(dd@dd);ma+=float(abs(delta).sum());mb+=float(abs(dd).sum());zero+=float(tt@tt);coordinates+=int(valid.sum())
  truth=d['y'][test];valid=np.ones(truth.shape,bool)
  for i,g in enumerate(d['targets'][test]):
   if g in gp:valid[i,gp[g]]=False
  zero_values=truth[valid];zero=float(zero_values@zero_values)
  values={'mse_vs_constant':err/const,'mae_vs_constant':ma/mb,'mse_vs_zero':err/zero}
  for k,v in values.items():assert abs(v-rec['test'][k])<1e-10,(name,variant,k,v,rec['test'][k])
  checks.append({'source':name,'variant':variant,'test_targets':len(test),'test_coordinates':coordinates,'replayed':True,'disjoint_targets':True})
 m.write(D/'independent_audit.json',{'passed':True,'plan_sha256':b.sha256(D/'plan.json'),'checks':checks,'auditor_sha256':b.sha256(__import__('pathlib').Path(__file__))});print('Unseen-target split, bias and prediction replays passed.')
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
