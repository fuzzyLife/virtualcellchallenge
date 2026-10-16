"""Audit V35/V36 fits, held-family scores and population-view identities."""
import argparse,json
import numpy as np
import torch
from threadpoolctl import threadpool_limits
import v35_methods as m
from run_v30_experiment import rows_for

def numpy_forward(net,b,gene):
 state={k:v.detach().cpu().numpy().astype(float) for k,v in net.state_dict().items()};b={k:v.detach().cpu().numpy().astype(float) for k,v in b.items()};gene=gene.detach().cpu().numpy().astype(float)
 def linear(x,name):return x@state[name+'.weight'].T+state[name+'.bias']
 def silu(x):return x/(1+np.exp(-np.clip(x,-700,700)))
 logit=b['prior']+linear(silu(linear(b['token'],'attention.0')),'attention.2').squeeze(-1);w=np.exp(logit-logit.max(1,keepdims=True));w/=w.sum(1,keepdims=True);w=w[:,:,None]*b['mask'];mass=w.sum(1);w/=np.maximum(mass[:,None,:],1e-12)
 value=b['value'];mu=(w*value).sum(1);sd=np.sqrt(np.maximum((w*(value-mu[:,None,:])**2).sum(1),1e-8));same=b['token'][:,:,-3,None]*w;sm=same.sum(1);sv=(same*value).sum(1)/np.maximum(sm,1e-12);sp=(w*b['src']).sum(1);base=b['base'];q=b['qlog']
 features=np.stack([base,mu,sd,(w*(value>0)).sum(1),(w*(value<0)).sum(1),sv,sm,q,sp,q-sp,base*mu,b['control_var'],b['dropout'],mu-base,b['moments_known']],axis=-1)
 n,g=base.shape;x=np.concatenate([features,np.broadcast_to(gene,(n,g,16)),np.broadcast_to(b['target'][:,None,:],(n,g,16)),np.broadcast_to(b['context'][:,None,:],(n,g,4))],axis=-1);out=linear(silu(linear(silu(linear(x,'decoder.0')),'decoder.2')),'decoder.4');gate=np.clip((np.expm1(q*10)-1)/5,0,1)*(mass>0)
 return base+.5*np.tanh(out[:,:,0])*gate,1.5*np.tanh(out[:,:,1])*gate

def audit(version):
 m.setup();m.check();data=m.old.old.read_cache();folder=m.A if version==35 else m.P36;result=json.loads((folder/'summary_results.json').read_text());fits=[];replays=[]
 for outer in m.OUTERS:
  net,g,rec=m.load(outer,version);assert outer not in {m.family(s) for s in rec['sources']};assert not set(rec['fitting_targets'])&set(rec['validation_targets'])
  for r in rec['training_records']:
   assert outer in r['excluded'] and r['family'] in r['excluded'];assert r['validation']==(m.target_split(r['target'])==0)
   for token in r['tokens']:assert token['source_family'] not in r['excluded']
  history=rec['history'];metric='validation_MSE_ratio' if version==35 else 'selection_loss';assert rec['selected_epoch']==min(history,key=lambda v:v[metric])['epoch']
  # Independently evaluate actual weight arrays on one real retrieval example.
  source=rec['training_records'][0]['source'];query=data['profiles'][m.CI[source]];target=rec['training_records'][0]['target'];ex=rec['training_records'][0]['excluded'];tr=m.retrieve(g,target,query,ex);axis=g['axis'][:37];base=m.base_for(query,ex);ti=int(np.flatnonzero(data['targets']==target)[0]);cols=np.arange(len(axis));a={k:np.array([v]) for k,v in tr.items() if k!='names'}
  for k in ['value','mask','src']:a[k]=a[k][:,:,cols]
  for k,v in [('base',base[ti,axis]),('qlog',np.log1p(np.nan_to_num(query[axis]))/10),('truth',np.zeros(len(axis))),('valid',np.ones(len(axis))),('control_var',np.zeros(len(axis))),('dropout',np.zeros(len(axis))),('moments_known',np.zeros(len(axis)))]:a[k]=np.array([v])
  a['gene']=np.array([m.embedding(t,16) for t in data['genes'][axis]]);b,ge=m.tensor_batch(a,[0],np.arange(len(axis)),'cpu');net=net.cpu();yp,vp,_=net(b,ge);yn,vn=numpy_forward(net,b,ge);np.testing.assert_allclose(yp.detach().numpy(),yn,atol=2e-6,rtol=2e-5);np.testing.assert_allclose(vp.detach().numpy(),vn,atol=2e-6,rtol=2e-5)
  fits.append({'outer':outer,'fit_targets':len(rec['fitting_targets']),'validation_targets':len(rec['validation_targets']),'selected_epoch':rec['selected_epoch'],'gene_outputs':len(g['axis']),'numpy_forward_max_error':float(abs(yp.detach().numpy()-yn).max())})
 for p in sorted((folder/'predictions').glob('*.npz')):
  rec=json.loads(p.with_suffix('.json').read_text());assert m.sha(p)==rec['sha256']
  with np.load(p,allow_pickle=False) as z:
   rows=z['rows'];axis=rec['axis'];outside=np.ones(len(data['genes']),bool);outside[axis]=False
   for arm in ['learned','response_shuffled','static_retrieval']:
    pred=z[arm].astype(float);base=z['base'].astype(float);np.testing.assert_array_equal(pred[:,outside],base[:,outside]);assert np.isfinite(pred).all() and abs(pred-base).max()<=.500002
    own=[int(np.flatnonzero(data['genes']==t)[0]) for t in data['targets']];assert (pred[np.arange(300),own]==0).all();se=be=ae=ab=0.
    for ri in rows:
     ti=int(data['row_targets'][ri]);mask=data['valid'][ri].copy();mask[own[ti]]=False;y=data['truth'][ri,mask].astype(float);e=pred[ti,mask]-y;b0=base[ti,mask]-y;se+=e@e;be+=b0@b0;ae+=abs(e).sum();ab+=abs(b0).sum()
    reported=next(r for r in result['records'] if r['source']==rec['source'] and r['arm']==arm);np.testing.assert_allclose([se/be,ae/ab],[reported['mse_ratio'],reported['mae_ratio']],rtol=2e-6);replays.append({'source':rec['source'],'arm':arm,'mse_ratio':se/be,'mae_ratio':ae/ab})
 if version==36:
  for p in (folder/'populations').glob('*.json'):
   rec=json.loads(p.read_text());outer=rec['outer']
   for source in rec['sources']:
    assert source['family']!=outer;ntc=set(source['control_rows'])
    for group in source['groups']:
     assert group['validation']==(m.target_split(group['target'])==0)
     v=group['views'][0];l,r=map(set,v['treated_rows']);cl,cr=map(set,v['control_rows']);assert not l&r and not cl&cr and not (l|r)&ntc and (cl|cr)<=ntc
 m.write(folder/'independent_audit.json',{'passed':True,'result_sha256':m.sha(folder/'summary_results.json'),'models':fits,'score_replays':replays,'checks':'All fitting/validation target and entire-family exclusions, retrieval labels, selected epochs, exact fallback, own exclusion, finite bounded outputs. NumPy forward independently reimplemented; all saved summary scores replayed in float64. First independent population pair per source/target checked from explicit row identities; other pairs asserted during creation.','auditor_sha256':m.sha(m.ROOT/'audit_v35_v36.py')});print('V'+str(version),'audit passed',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('version',type=int,choices=[35,36]);args=p.parse_args()
 with threadpool_limits(limits=2):audit(args.version)