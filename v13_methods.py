"""Training-family-only low-rank perturbation response programs."""
import numpy as np
import warnings
from sklearn.utils.extmath import randomized_svd
from v11_methods import response_bank
from v7_response_methods import collapse, combine
import build_submission as build
ROOT=build.ROOT
OUT=ROOT/'experiments/v13_response_programs'
A=OUT/'analysis'
NAMES=['v11']+[f'program{rank}_blend{blend}' for rank in [8,32,64] for blend in [25,50]]

def banks(data,excluded=None,reserved=()):
 original,record=response_bank(data,excluded)
 base=original['consensus20'].copy()
 train=np.ones(len(data['truth']),bool) if excluded is None else data['row_families']!=excluded
 y,var,valid,kd=collapse(data,train)
 gain=np.clip(.9/np.maximum(kd,.35),.8,1.8)
 y*=gain[:,:,None];var*=gain[:,:,None]**2
 mean,variance,ok=combine((y,var,valid,kd),np.array([2.,.25,2.,1.,1.]))
 own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']])
 base[np.arange(len(base)),own]=0
 fitting=~np.isin(np.arange(len(base)),reserved)
 # Whitening by source uncertainty keeps high-noise genes from defining the basis.
 # The held-out family and reserved targets never define the learned programs.
 with warnings.catch_warnings():
  warnings.simplefilter("ignore",RuntimeWarning)
  scale=np.sqrt(np.maximum(np.nanmedian(np.where(ok[fitting],variance[fitting],np.nan),axis=0),.01))
 scale=np.nan_to_num(scale,nan=1.,posinf=1.)
 scale=np.clip(scale,.1,2.)
 matrix=base[fitting]/scale
 u,s,vt=randomized_svd(matrix,n_components=64,n_iter=5,random_state=2026101013)
 result={'v11':{'consensus20':base}}
 for rank in [8,32,64]:
  basis=vt[:rank]
  projection=((base/scale)@basis.T)@basis*scale
  projection[~ok]=0
  projection[np.arange(len(base)),own]=0
  for blend in [25,50]:
   effect=np.clip(base+(blend/100)*(projection-base),-3,3).astype(np.float32)
   result[f'program{rank}_blend{blend}']={'consensus20':effect}
 record.update({'program_basis_targets':np.flatnonzero(fitting).tolist(),
  'singular_values':s.tolist(),'programs':64,'whitening':'Median valid training-source uncertainty; sqrt variance clipped .1..2',
  'source_cache_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')})
 return result,record
