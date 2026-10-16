"""Independent NumPy replay of three held-family learned flow predictions."""
import json,gc
from functools import lru_cache
import numpy as np
from threadpoolctl import threadpool_limits
import v32_methods as m
import v34_corrected_methods as f

def replay(net,a,template,query,target):
 state=net.state_dict();w={k:v.detach().numpy().astype(float) for k,v in state.items()};center=np.log1p(np.nan_to_num(query[a['axis']])*.01).astype(np.float32);x0=(f.logspace(template,a['axis'])-center).astype(float);x=x0.copy();c=np.repeat(np.r_[f.target_feature(target),m.ctx(a,query).astype(np.float32)][None,:],len(x),0).astype(float);own=np.flatnonzero(m.old.read_cache()['genes'][a['axis']]==target);own=int(own[0]) if len(own) else -1
 for step in range(8):
  v=np.concatenate([x,x0,np.full((len(x),1),step/8),c],1)
  for layer in [0,2]:
   v=v@w[f'body.{layer}.weight'].T+w[f'body.{layer}.bias'];v=v/(1+np.exp(-np.clip(v,-50,50)))
  velocity=v@w['body.4.weight'].T+w['body.4.bias'];x+=velocity/8
  if own>=0:x[:,own]=x0[:,own]
 return np.expm1(np.clip(x+center,0,12))*100

def main():
 m.old.embeddings=lru_cache(maxsize=1)(m.old.embeddings);p=f.check();rows=[];data=m.old.read_cache()
 for key in m.OUTERS:
  ctx=next(c for c in p['contexts'] if f.family(c['family'])==key);cell=ctx['cells'][0];net,a,r=f.load_model(key);x,meta,measured,source=f.raw(ctx['family'],data['genes']);query=f.normmean(x[ctx['left']]);query[~measured]=np.nan;template=x[cell['template']];expected=f.generated(net,a,template,query,cell['target']);actual=replay(net,a,template,query,cell['target']);np.testing.assert_allclose(actual,expected,rtol=2e-5,atol=.01);rows.append({'outer':key,'context':ctx['family'],'target':cell['target'],'cells':len(cell['template']),'modeled_genes':len(a['axis']),'maximum_CPM_absolute_difference':float(abs(actual-expected).max()),'mean_CPM_absolute_difference':float(abs(actual-expected).mean())});print('V34 NumPy inference replay',key,'passed',flush=True);del x;gc.collect()
 m.write(f.A/'numpy_inference_audit.json',{'passed':True,'records':rows,'auditor_sha256':m.sha(m.ROOT/'audit_v34_inference.py'),'scope':'Independent NumPy forward/SiLU/Euler integration, own-state lock and RNA inverse transform versus actual Torch query inference. Three representative held-family cases; no independent retraining or all-case count regeneration.'})
if __name__=='__main__':
 with threadpool_limits(limits=2):main()