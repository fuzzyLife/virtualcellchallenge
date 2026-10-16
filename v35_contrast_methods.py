"""Conservative contrasts that remove neural generic-response drift.

Exploratory after primary outcomes. No target treatment labels enter these
transforms; the complete300-query panel defines the target-centering constraint.
"""
import numpy as np

def center_delta(pred,base,support,cap=.5):
 d=np.asarray(pred,float)-base;mask=np.asarray(support,bool);lower=np.maximum(-cap,-3-base);upper=np.minimum(cap,3-base);active=mask.any(0)
 lo=np.min(np.where(mask,d-upper,np.inf),axis=0);hi=np.max(np.where(mask,d-lower,-np.inf),axis=0);lo[~active]=0;hi[~active]=0
 for _ in range(42):
  mid=(lo+hi)/2;z=np.clip(d-mid,lower,upper)*mask;positive=z.sum(0)>0;lo=np.where(positive,mid,lo);hi=np.where(positive,hi,mid)
 z=np.clip(d-(lo+hi)/2,lower,upper)*mask;assert np.max(abs(z.sum(0)),initial=0)<1e-7
 return z

def transform(pred,base,shuffled,genes,targets,axis,trace,query,kind):
 support=np.zeros(base.shape,bool);has=np.array([bool(r['tokens']) for r in trace]);support[:,axis]=has[:,None]&(np.nan_to_num(query[axis])>1)[None,:]
 for i,t in enumerate(targets):support[i,np.asarray(genes)==t]=False
 if kind=='target_centered25':delta=center_delta(pred,base,support)
 elif kind=='retrieval_contrast25':delta=np.clip(np.asarray(pred,float)-shuffled,-.5,.5)*support
 else:raise ValueError(kind)
 out=np.clip(base+.25*delta,-3,3).astype(np.float32)
 for i,t in enumerate(targets):out[i,np.asarray(genes)==t]=0
 np.testing.assert_array_equal(out[~support],base[~support]);assert abs(out-base).max()<=.125002
 return out