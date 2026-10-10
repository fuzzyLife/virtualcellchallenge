"""Robust source consensus: cap isolated lineage effects before empirical Bayes."""
import numpy as np,warnings
from v11_methods import response_bank
from v7_response_methods import collapse,adaptive_shrink
NAMES=['v11','robust025','robust050','student3']

def banks(data,excluded=None,reserved=()):
 original,record=response_bank(data,excluded);result={'v11':{'consensus20':original['consensus20']}}
 train=np.ones(len(data['truth']),bool) if excluded is None else data['row_families']!=excluded
 y,var,valid,kd=collapse(data,train);gain=np.clip(.9/np.maximum(kd,.35),.8,1.8);y*=gain[:,:,None];var*=gain[:,:,None]**2
 initial=np.array([2.,.25,2.,1.,1.])[:,None,None]*np.clip(.04/np.maximum(var,.0025),.1,4.)*valid
 with warnings.catch_warnings():
  warnings.simplefilter('ignore',RuntimeWarning)
  center=np.nanmedian(np.where(valid,y,np.nan),axis=0);noise=np.sqrt(np.nanmedian(np.where(valid,var,np.nan),axis=0))
 center=np.nan_to_num(center);noise=np.nan_to_num(noise,nan=1.)
 own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']]);support=valid.sum(axis=0)
 for name in NAMES[1:]:
  yy=y.copy();w=initial.copy()
  if name.startswith('robust'):
   floor=float(name[6:])/100;radius=np.maximum(floor,2*noise)
   # Do not clip a single source: no independent evidence of an outlier exists.
   yy=np.where((support>=3)[None,:],np.clip(yy,center[None,:]-radius,center[None,:]+radius),yy)
  else:
   location=center.copy()
   for _ in range(5):
    extra=np.maximum(.25**2,noise**2)
    w=initial*4/(3+(y-location[None,:])**2/(var+extra[None,:]))
    den=w.sum(axis=0);wn=np.divide(w,den[None,:],out=np.zeros_like(w),where=den[None,:]>0)
    location=(wn*y).sum(axis=0)
  den=w.sum(axis=0);wn=np.divide(w,den[None,:],out=np.zeros_like(w),where=den[None,:]>0);ok=den>0
  mean=(wn*yy).sum(axis=0)
  # Clipping source means never claims to remove measurement or lineage uncertainty.
  variance=np.maximum((wn*wn*var).sum(axis=0)+.25*(wn*(y-mean[None,:])**2).sum(axis=0),.0025)
  mean[~ok]=0;variance[~ok]=1e6
  eb,_=adaptive_shrink(mean,variance,ok)
  signed=(w*np.sign(y)).sum(axis=0);agreement=np.divide(abs(signed),den,out=np.zeros_like(mean),where=den>0)
  gate=np.clip((agreement-.5)/.5,0,1)*np.minimum(abs(mean)/np.sqrt(variance),1.)*(support>=2)
  effect=np.clip(eb+.2*gate*(mean-eb),-3,3).astype(np.float32);effect[np.arange(len(own)),own]=0
  result[name]={'consensus20':effect}
 record.update(robust_policy='Only cap isolated source values with >=3 measured families; retain original source disagreement in variance. Student-t iterates reliability without treating conditions as independent families.')
 return result,record
