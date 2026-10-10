"""Assay-aware knockdown calibration and intervention-quality source weighting."""
import numpy as np
from v11_methods import response_bank
from v7_response_methods import collapse,adaptive_shrink,FAMILIES
NAMES=['v11','unknown_kd','quality_soft','quality50']

def banks(data,excluded=None,reserved=()):
 original,record=response_bank(data,excluded);result={'v11':{'consensus20':original['consensus20']}}
 train=np.ones(len(data['truth']),bool) if excluded is None else data['row_families']!=excluded
 own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']])
 control=data['profiles'][data['row_contexts'],own[data['row_targets']]]
 known=np.isfinite(control)&(control>5)
 adjusted=dict(data);adjusted['knockdown']=data['knockdown'].copy()
 # Missing or unquantifiable on-target readout cannot mean measured failed KD.
 adjusted['knockdown'][~known]=.9
 y,var,valid,kd=collapse(adjusted,train)
 measured=np.zeros(kd.shape,bool)
 for fi,family in enumerate(FAMILIES):
  for ti in range(len(own)):
   rows=train&(data['row_families']==family)&(data['row_targets']==ti)
   # Mixed conditions only qualify if every contributing on-target readout is quantifiable.
   measured[fi,ti]=rows.any() and known[rows].all()
 gain=np.clip(.9/np.maximum(kd,.35),.8,1.8);y*=gain[:,:,None];var*=gain[:,:,None]**2
 for name in NAMES[1:]:
  mask=valid.copy();quality=np.ones(kd.shape)
  if name=='quality_soft':quality=np.where(measured,np.clip(kd/.8,0,1),1.)
  if name=='quality50':mask &= ((~measured)|(kd>=.5))[:,:,None]
  w=np.array([2.,.25,2.,1.,1.])[:,None,None]*quality[:,:,None]*np.clip(.04/np.maximum(var,.0025),.1,4.)*mask
  den=w.sum(axis=0);wn=np.divide(w,den[None,:],out=np.zeros_like(w),where=den[None,:]>0);ok=den>0
  mean=(wn*y).sum(axis=0);variance=np.maximum((wn*wn*var).sum(axis=0)+.25*(wn*(y-mean[None,:])**2).sum(axis=0),.0025)
  mean[~ok]=0;variance[~ok]=1e6
  eb,_=adaptive_shrink(mean,variance,ok)
  agreement=np.divide(abs((w*np.sign(y)).sum(axis=0)),den,out=np.zeros_like(mean),where=den>0)
  support=(mask&(quality[:,:,None]>0)).sum(axis=0)
  gate=np.clip((agreement-.5)/.5,0,1)*np.minimum(abs(mean)/np.sqrt(variance),1.)*(support>=2)
  effect=np.clip(eb+.2*gate*(mean-eb),-3,3).astype(np.float32);effect[np.arange(len(own)),own]=0
  result[name]={'consensus20':effect}
 record.update({'unknown_or_unquantifiable_training_KD_rows':int((train&~known).sum()),
  'quality_policy':'Known KD requires finite source control own-gene CPM >5. Unknown gains are neutral (KD=.9), never amplified as a failed knockdown. Soft weight=KD/.8 clipped 0..1; hard gate>=.5 for known measurements only.'})
 return result,record
