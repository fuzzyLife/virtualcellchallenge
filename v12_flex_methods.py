"""Matched-Flex transfer within the existing K562 family, never as a sixth lineage."""
import json
import numpy as np
from v7_response_methods import collapse,combine,adaptive_shrink
from v11_methods import consensus_gate,response_bank
import build_submission as build
ROOT=build.ROOT;D=ROOT/'reference/v12_resources/flex_source'
VARIANTS=['v11','flex_combine','flex_prefer','flex_overlap']
def banks(data,excluded=None):
 original,training=response_bank(data,excluded=excluded)
 result={'v11':{'consensus20':original['consensus20']}}
 y,var,valid,kd=collapse(data,np.ones(len(data['truth']),bool) if excluded is None else data['row_families']!=excluded)
 gain=np.clip(.9/np.maximum(kd,.35),.8,1.8);y*=gain[:,:,None];var*=gain[:,:,None]**2
 rec=json.loads((ROOT/'experiments/v12_distribution_calibrated/analysis/flex_extraction.json').read_text());assert build.sha256(D/'summaries.npz')==rec['summary_sha256']
 with np.load(D/'summaries.npz') as z:
  fy=z['truth'].copy();fv=z['errors']**2;fm=z['valid'].copy();fkd=z['knockdown'].copy()
 fg=np.clip(.9/np.maximum(fkd,.35),.8,1.8);fy*=fg[:,None];fv*=fg[:,None]**2
 for variant in VARIANTS[1:]+['flex_shuffled']:
  yy=y.copy();vv=var.copy();mm=valid.copy()
  if excluded!='K562':
   if variant=='flex_shuffled':
    order=np.random.default_rng(2026100953).permutation(len(fy));ffy=fy[order];ffv=fv[order];ffm=fm[order]
   else:ffy=fy;ffv=fv;ffm=fm.copy()
   if variant=='flex_overlap':ffm &= valid[0].any(axis=0)[None,:]
   both=valid[0]&ffm;only=(~valid[0])&ffm
   if variant=='flex_prefer':
    yy[0,ffm]=ffy[ffm];vv[0,ffm]=ffv[ffm]
   else:
    yy[0,both]=.5*(y[0,both]+ffy[both])
    # Correlated releases/platforms: retain marginal uncertainty and disagreement;
    # never divide variance by two or add an independent-family support vote.
    vv[0,both]=.5*(var[0,both]+ffv[both])+.25*(y[0,both]-ffy[both])**2
    yy[0,only]=ffy[only];vv[0,only]=ffv[only]
   mm[0]|=ffm
  mean,variance,ok=combine((yy,vv,mm,kd),np.array([2.,.25,2.,1.,1.]))
  eb,_=adaptive_shrink(mean,variance,ok);gate,_,_=consensus_gate(yy,vv,mm,kd,[2.,.25,2.,1.,1.],mean,variance)
  effect=np.clip(eb+.2*gate*(mean-eb),-3,3).astype(np.float32)
  own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']]);effect[np.arange(len(own)),own]=0
  result[variant]={'consensus20':effect}
  if excluded=='K562':np.testing.assert_array_equal(effect,result['v11']['consensus20'])
 return result,{**training,'Flex_added_as_independent_family':False,'Flex_summary_sha256':rec['summary_sha256'],'variance_policy':'Mean marginal platform variance plus platform disagreement; no independent replicate SE reduction.'}
