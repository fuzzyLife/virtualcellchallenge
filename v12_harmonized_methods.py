"""Training-only platform offset correction and concordance for the added K562 source."""
import numpy as np,json
import build_submission as build
from v7_response_methods import collapse,combine,adaptive_shrink
from v11_methods import response_bank,consensus_gate
ROOT=build.ROOT
VARIANTS=['v11','flex_centered','flex_agreement','flex_reliability']
def banks(data,excluded=None):
 original,training=response_bank(data,excluded=excluded);base=original['consensus20'];result={'v11':{'consensus20':base}}
 if excluded=='K562':return {n:{'consensus20':base.copy()} for n in VARIANTS+['centered_shuffle']},training
 y,var,valid,kd=collapse(data,np.ones(len(data['truth']),bool) if excluded is None else data['row_families']!=excluded)
 gain=np.clip(.9/np.maximum(kd,.35),.8,1.8);y*=gain[:,:,None];var*=gain[:,:,None]**2
 rec=json.loads((ROOT/'experiments/v12_distribution_calibrated/analysis/flex_extraction.json').read_text());file=ROOT/'reference/v12_resources/flex_source/summaries.npz';assert build.sha256(file)==rec['summary_sha256']
 with np.load(file) as z:fy=z['truth'].copy();fv=z['errors']**2;fm=z['valid'].copy();fkd=z['knockdown'].copy()
 fg=np.clip(.9/np.maximum(fkd,.35),.8,1.8);fy*=fg[:,None];fv*=fg[:,None]**2
 both=valid[0]&fm;paired=(fy-y[0])*both;den=both.sum(axis=0)[None,:]-both;total=paired.sum(axis=0)[None,:]-paired
 # Leave the queried target out of the gene-wise platform offset; 20-pair regularizer.
 offset=np.divide(total,den+20,out=np.zeros_like(total),where=den>=20)
 total_new=(fy*fm).sum(axis=0)[None,:]-fy*fm;den_new=fm.sum(axis=0)[None,:]-fm
 generic=np.divide(total_new,den_new+20,out=np.zeros_like(total_new),where=den_new>=20)
 centered=fy-np.where(valid[0],offset,generic)
 own=np.array([np.flatnonzero(data['genes']==t)[0] for t in data['targets']])
 for name in VARIANTS[1:]+['centered_shuffle']:
  yy=y.copy();vv=var.copy();mm=valid.copy();ffy=centered if name in ['flex_centered','centered_shuffle'] else fy;ffm=fm.copy();ffv=fv.copy()
  if name=='centered_shuffle':
   order=np.random.default_rng(2026100957).permutation(len(fy));ffy=ffy[order];ffm=ffm[order];ffv=ffv[order]
  if name=='flex_reliability':vv[0,both]+=.25*(fy[both]-y[0,both])**2
  else:
   if name=='flex_agreement':
    other,_,available=combine((y[1:],var[1:],valid[1:],kd[1:]),np.array([.25,2.,1.,1.]))
    ffm &= np.where(valid[0],np.sign(ffy)==np.sign(y[0]),available&(np.sign(ffy)==np.sign(other)))
   shared=valid[0]&ffm;added=(~valid[0])&ffm
   yy[0,shared]=.5*(y[0,shared]+ffy[shared]);vv[0,shared]=.5*(var[0,shared]+ffv[shared])+.25*(y[0,shared]-ffy[shared])**2
   yy[0,added]=ffy[added];vv[0,added]=ffv[added];mm[0]|=ffm
  mean,variance,ok=combine((yy,vv,mm,kd),np.array([2.,.25,2.,1.,1.]));eb,_=adaptive_shrink(mean,variance,ok);gate,_,_=consensus_gate(yy,vv,mm,kd,[2.,.25,2.,1.,1.],mean,variance)
  effect=np.clip(eb+.2*gate*(mean-eb),-3,3).astype(np.float32);effect[np.arange(len(own)),own]=0;result[name]={'consensus20':effect}
 return result,{**training,'platform_correction':'Leave-target-out paired offset, minimum 20 other targets, 20-pair regularizer; unpaired new genes subtract the leave-target-out Flex generic response. No held-family observations used.'}
