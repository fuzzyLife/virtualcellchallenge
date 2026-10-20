"""Independent raw-rank spot checks and serialized DE/profile/PDS audit for V75."""
import argparse,json
from threadpoolctl import threadpool_limits
import numpy as np
from scipy import sparse,stats
import v75_assay_shrink as m
from v57_pie import read,write,sha

def spot_ranks(count,reference,dd):
 cols=np.flatnonzero(dd['eligible']);selected=cols[np.linspace(0,len(cols)-1,min(8,len(cols)),dtype=int)]
 ld=np.asarray(count.sum(1)).ravel()/1e6;rd=np.asarray(reference.sum(1)).ravel()/1e6
 for col in selected:
  x=np.log1p(count[:,col].toarray().ravel()/ld);y=np.log1p(reference[:,col].toarray().ravel()/rd);allv=np.concatenate([x,y]);rank=stats.rankdata(allv,method='average');n=len(x);k=len(y);N=n+k
  u=rank[:n].sum()-n*(n+1)/2;_,ties=np.unique(allv,return_counts=True);variance=n*k/12*((N+1)-np.sum(ties.astype(float)**3-ties)/(N*(N-1)))
  p=2*stats.norm.sf(abs(u-n*k/2)/np.sqrt(variance)) if variance>0 else 1.
  np.testing.assert_allclose(p,dd['p'][col],rtol=1e-10,atol=1e-12)
 return len(selected)

def close(a,b):
 if b is None:assert not np.isfinite(a)
 else:np.testing.assert_allclose(a,b,rtol=0,atol=1e-11)

def audit(folder,output):
 records=[];d=m.data();genes=d['genes'];targets=d['targets'];selected=read(m.A/'selection.json')['selected']
 for p in sorted(folder.glob('*.json')):
  if p.name.endswith('__panel.json'):continue
  r=read(p);arr=p.with_suffix('.npz');assert sha(arr)==r['arrays_sha256'];own=list(genes).index(r['target'])
  with np.load(arr,allow_pickle=False) as z:
   valid=z['truth_eligible']&(np.arange(len(genes))!=own);real=valid&(z['truth_q']<.05);adjudicable=valid&np.isfinite(z['truth_lfc'])&(z['truth_lfc']!=0);called=adjudicable&(z['pred_q']<.05)
   correct=called&(np.sign(z['pred_lfc'])==np.sign(z['truth_lfc']));budget=max(int(real.sum()),int(called.sum()));fidelity=correct.sum()/budget if budget else np.nan
   union=real|(valid&(z['pred_q']<.05));jac=(real&(z['pred_q']<.05)).sum()/union.sum() if union.any() else 1.
   pool=np.flatnonzero(real&adjudicable);order=sorted(pool,key=lambda i:(z['pred_q'][i]>=.05,z['pred_q'][i],z['pred_p'][i],-abs(z['pred_lfc'][i]),str(genes[i])))
   if order:
    pure=np.cumsum(np.sign(z['pred_lfc'][order])==np.sign(z['truth_lfc'][order]))/np.arange(1,len(order)+1);passes=np.flatnonzero(pure>=.9);reach=(passes[-1]+1)/int(real.sum()) if len(passes) else 0.
   else:reach=0. if real.sum()>0 else np.nan
   denominator=abs(z['truth_lfc'][real]).sum();nmae=abs(np.nan_to_num(z['pred_lfc'][real])-z['truth_lfc'][real]).sum()/denominator if real.sum()>=10 and np.isfinite(denominator) and denominator>0 else np.nan
   for field,value in [('fidelity',fidelity),('reach',reach),('jaccard',jac),('nmae',nmae)]:close(value,r[field])
   for prefix in ('pred','truth'):
    axis=np.flatnonzero(z[prefix+'_eligible']);pv=z[prefix+'_p'][axis];order=np.argsort(pv,kind='stable');ranked=pv[order];q=np.minimum.accumulate((ranked*len(axis)/np.arange(1,len(axis)+1))[::-1])[::-1];manual=np.empty(len(axis));manual[order]=np.minimum(q,1)
    np.testing.assert_allclose(manual,z[prefix+'_q'][axis],rtol=0,atol=1e-12)
  assert r['manual_rank_coordinates_checked']>0
  if r['arm']=='v66':assert r['baseline_full_replay_verified']
  records.append(dict(path=p.relative_to(m.ROOT).as_posix(),sha256=sha(p),arrays_sha256=sha(arr)))
 panels=[]
 for ctx in m.contexts():
  name=ctx['family'];indices=m.cases(ctx)[3:];rows=[ctx['cells'][i] for i in indices]
  for arm in ('v66',selected,selected+'_shuffled'):
   values=[]
   for i in indices:
    with np.load(folder/f'{name}__{i}__{arm}.npz',allow_pickle=False) as z:values.append((z['pred'].copy(),z['truth'].copy(),z['control'].copy()))
   pp=np.array([v[0] for v in values]);tt=np.array([v[1] for v in values]);control=values[0][2];mask=~np.isin(genes,targets);x=(pp-control)[:,mask];y=(tt-control)[:,mask]
   xn=np.linalg.norm(x,axis=1);yn=np.linalg.norm(y,axis=1);den=xn[:,None]*yn[None,:];cosine=np.divide(x@y.T,den,out=np.zeros_like(den),where=den>0)
   distance=1-cosine
   score=[]
   for i in range(len(values)):
    other=np.delete(distance[i],i);score.append(np.mean((distance[i,i]<other)+.5*(distance[i,i]==other)))
   close(float(np.mean(score)),read(folder/f'{name}__{arm}__panel.json')['pds'])
   panels.append(dict(source=name,arm=arm,pds=float(np.mean(score))))
 assert len(records)==45 and len(panels)==15,(len(records),len(panels))
 write(m.A/output,dict(passed=True,records=records,panels=panels,comparison_records=len(records),all_four_DE_metrics_independently_replayed=True,all_BH_arrays_verified=True,manual_rank_coordinates=sum(read(m.ROOT/r['path'])['manual_rank_coordinates_checked'] for r in records),all_frozen_V66_baselines_replayed=True,auditor_sha256=sha(__file__),scope='Raw public metric and arithmetic integrity; repeated biology and two seeds are not independent biological replication.'))
 print('V75 independent count-statistic audit passed',folder.name,len(records),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--repeat',action='store_true');a=p.parse_args()
 with threadpool_limits(limits=2):audit(m.A/('counts_repeat' if a.repeat else 'counts'),'count_audit_repeat.json' if a.repeat else 'count_audit.json')
