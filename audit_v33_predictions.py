"""Replay V33 saved query predictions with independent float64 score sums."""
import json,hashlib
from functools import lru_cache
import numpy as np
from threadpoolctl import threadpool_limits
import v32_methods as m
import v32_fallback_correction as correction
import v32_runtime as runtime
import run_v33_experiment as t
from run_v30_experiment import rows_for

def main():
 t.check();runtime.install();correction.install();m.old.embeddings=lru_cache(maxsize=1)(m.old.embeddings);data=m.old.read_cache();result=json.loads((t.A/'results.json').read_text());replays=[]
 for outer in m.OUTERS:
  file=t.A/'predictions'/f'{outer}.npz';meta=json.loads(file.with_suffix('.json').read_text());assert m.sha(file)==meta['sha256'];g=m.fit_geometry([outer]);g['basis']=g['basis'][:,:8];d=m.paired_table(g);g.update(feature_mean=d['feature_mean'],feature_scale=d['feature_scale'])
  with np.load(file,allow_pickle=False) as z:saved={k:z[k].copy() for k in z.files}
  np.testing.assert_array_equal(g['axis'],saved['axis']);np.testing.assert_allclose(g['basis'],saved['basis'],atol=1e-12);np.testing.assert_allclose(d['x'][saved['fit_indices']],saved['trainx'],atol=2e-6)
  # Explicit observed-coordinate least squares checks the saved target programs.
  for i in range(min(6,len(saved['fit_indices']))):
   ri=int(saved['fit_indices'][i]);own=int(d['own'][ri]);keep=np.arange(len(g['axis']))!=own;coef=np.linalg.lstsq(g['basis'][keep],(d['y'][ri]/g['scale'])[keep],rcond=None)[0];np.testing.assert_allclose(coef,saved['trainy'][i],atol=1e-8)
  keys=[hashlib.sha256(v.tobytes()).hexdigest() for v in saved['queryx']];nva=len(saved['validation_indices'])
  for arm in ['ridge','catboost','tabicl','identity_shuffled_tabicl']:
   lookup=dict(zip(keys,saved[arm][nva:]));
   def replay(gg,x):return np.array([lookup[hashlib.sha256(np.asarray(v,np.float32).tobytes()).hexdigest()] for v in x])
   m.predict_delta=replay
   for alpha in [.5,1.]:
    vals=[];rows=rows_for(data,outer)
    for ci in np.unique(data['row_contexts'][rows]):
     pred,base,_=correction.predict(g,data['profiles'][ci],alpha);outside=np.ones(len(data['genes']),bool);outside[g['axis']]=False;np.testing.assert_array_equal(pred[:,outside],base[:,outside]);rr=rows[data['row_contexts'][rows]==ci];se=be=ae=ab=0.
     for ri in rr:
      ti=int(data['row_targets'][ri]);valid=data['valid'][ri].copy();valid[data['genes']==data['targets'][ti]]=False;a=(pred[ti,valid]-data['truth'][ri,valid]).astype(float);b=(base[ti,valid]-data['truth'][ri,valid]).astype(float);se+=a@a;be+=b@b;ae+=abs(a).sum();ab+=abs(b).sum()
     vals.append((se/be,ae/ab))
    pair=np.mean(vals,0);r=next(v for v in result['records'] if v['outer']==outer and v['arm']==arm and v['alpha']==alpha);assert abs(pair[0]-r['mse_ratio'])<2e-6;assert abs(pair[1]-r['mae_ratio'])<2e-6;replays.append({'outer':outer,'arm':arm,'alpha':alpha,'MSE_ratio':float(pair[0]),'MAE_ratio':float(pair[1])})
  print('V33 independent prediction replay',outer,'passed',flush=True)
 m.write(t.A/'prediction_replay_audit.json',{'passed':True,'results_sha256':m.sha(t.A/'results.json'),'replays':replays,'auditor_sha256':m.sha(m.ROOT/'audit_v33_predictions.py'),'scope':'All24 outer means replayed from saved query tokens/deltas, fitting features and own-masked least squares checked, exact V24 fallback outside modeled axis. Saved TabICL GPU deltas reused; not an independent implementation of foundation-model pretraining/inference.'})
if __name__=='__main__':
 with threadpool_limits(limits=2):main()