"""Independent rank-sum/BH/reach replay on three regenerated count panels."""
import numpy as np
from scipy import sparse,stats
from threadpoolctl import threadpool_limits
import v45_bounded_dispersion as m
from run_v35_v36 import contexts
from run_v37_ensemble import checked_folder
from evaluate_atlas_transfer import load_atlas_model
def manual(x,control):
 x=x.toarray().astype(float);c=control.toarray().astype(float);n=len(x);nc=len(c);N=n+nc
 x/=x.sum(1,keepdims=True)/1e6;c/=c.sum(1,keepdims=True)/1e6;xm=x.mean(0);cm=c.mean(0);valid=cm>5;p=np.ones(x.shape[1])
 x=np.log1p(x);c=np.log1p(c)
 for g in np.flatnonzero(valid):
  v=np.r_[x[:,g],c[:,g]];r=stats.rankdata(v,method='average');_,t=np.unique(v,return_counts=True);variance=n*nc/12*(N+1-np.sum(t.astype(float)**3-t)/(N*(N-1)))
  u=r[:n].sum()-n*(n+1)/2;p[g]=min(1.,2*stats.norm.sf(abs(u-n*nc/2)/np.sqrt(variance))) if variance>0 else 1.
 q=np.ones_like(p);idx=np.flatnonzero(valid);order=idx[np.argsort(p[idx],kind='stable')];adjust=p[order]*len(order)/np.arange(1,len(order)+1);q[order]=np.minimum(1,np.minimum.accumulate(adjust[::-1])[::-1])
 return {'p':p,'q':q,'lfc':np.log2((xm+1e-9)/(cm+1e-9)),'eligible':valid}
def reach(pred,truth,genes,own):
 valid=truth['eligible'].copy();valid[own]=False;real=valid&(truth['q']<.05);pool=np.flatnonzero(real&np.isfinite(truth['lfc'])&(truth['lfc']!=0))
 pool=sorted(pool,key=lambda i:(pred['q'][i]>=.05,pred['q'][i],pred['p'][i],-abs(pred['lfc'][i]),str(genes[i])))
 successes=0;furthest=0
 for j,i in enumerate(pool,1):
  successes+=np.sign(pred['lfc'][i])==np.sign(truth['lfc'][i])
  if successes/j>=.9:furthest=j
 return furthest/int(real.sum()) if real.any() else np.nan
def main():
 m.check();selection=m.read(m.A/'locked_finalist.json')['selected'];d=m.data();genes=d['genes'];atlas,_=load_atlas_model();records=[]
 for name in ['H1','HEK293T','CD4_Rest_D1']:
  ctx=next(c for c in contexts() if c['family']==name);folder,_=checked_folder(name);cell=ctx['cells'][0];own=int(np.flatnonzero(genes==cell['target'])[0]);ti=int(np.flatnonzero(d['targets']==cell['target'])[0]);template=sparse.load_npz(folder/'0_template.npz');ref=sparse.load_npz(folder/'reference.npz');actual=sparse.load_npz(folder/'0_actual.npz')
  with np.load(folder/'effects.npz',allow_pickle=False) as z:cp=z['cp'];bulk=z['bulk']
  with np.load(m.A/'predictions'/(name+'.npz'),allow_pickle=False) as z:effect=z[selection][0]
  libs=template.sum(1).A1.astype(np.int64);count,diag=m.emit(template,effect,own,cell['seed']+ti,(libs,cp,bulk),atlas,selection)
  pred=manual(count,ref);truth=manual(actual,ref);file=m.A/'development'/f'{name}__0__0__{selection}.npz';saved=m.read(file.with_suffix('.json'));assert m.b.sha256(file)==saved['arrays_sha256']
  with np.load(file,allow_pickle=False) as z:
   for group,v in [('pred',pred),('truth',truth)]:
    for k in ['q','lfc']:np.testing.assert_allclose(v[k],z[group+'_'+k],rtol=1e-9,atol=1e-10)
  rr=reach(pred,truth,genes,own);assert np.isfinite(rr);np.testing.assert_allclose(rr,saved['reach'],atol=1e-12)
  np.testing.assert_array_equal(count.sum(1).A1,libs);records.append({'context':name,'target':cell['target'],'all_eligible_genes_manual_rank_sum_BH':True,'reach_independently_replayed':rr,'libraries_exact':True});print('V45 raw ranks/BH/reach PASS',name,flush=True)
 m.write(m.A/'raw_score_audit.json',{'passed':True,'records':records,'script_sha256':m.b.sha256(m.ROOT/'audit_v45_raw_scores.py'),'scope':'Independent raw rank-sum tie correction, BH and reach on three regenerated populations, not private anchors.'})
if __name__=='__main__':
 with threadpool_limits(limits=2):main()


