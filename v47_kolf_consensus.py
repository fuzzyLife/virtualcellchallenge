"""V47 corroboration-only KOLF response amplification, without z-to-LFC conversion."""
from functools import lru_cache
import numpy as np
import v38_response_hurdle as parent
ROOT=parent.ROOT;OUT=ROOT/'experiments/v47_kolf_consensus';A=OUT/'analysis'
ARMS=('kolf_consensus25','shuffled_consensus25');OUTERS=parent.OUTERS
b=parent.b;data=parent.data;read=parent.read;write=parent.write
def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,n
 return p
@lru_cache(1)
def evidence():
 import v39_kolf_hurdle as source
 z,sat=source.kolf();supported=np.any(np.isfinite(z),axis=1);ids=np.flatnonzero(supported);permutation=ids[np.random.default_rng(47047).permutation(len(ids))]
 shuffled=z.copy();ss=sat.copy();shuffled[ids]=z[permutation];ss[ids]=sat[permutation]
 return z,sat,shuffled,ss
def combine(base,z,sat):
 mask=np.isfinite(z)&~sat&(z*base>0)&(abs(base)>=.02)
 confidence=np.where(mask,np.tanh(abs(np.nan_to_num(z))),0)
 return np.clip(base*(1+.25*confidence),-3,3).astype(np.float32)
def predict(base):
 z,sat,sh,ss=evidence();d=data();out={ARMS[0]:combine(base,z,sat),ARMS[1]:combine(base,sh,ss)}
 for a,v in out.items():
  for ti,t in enumerate(d['targets']):v[ti,d['genes']==t]=0
 return out
def emit(template,effect,own,seed,prepared,atlas,arm):return parent.emit(template,effect,own,seed,prepared,atlas,False)
def freeze():
 from run_v35_v36 import contexts
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v47_kolf_consensus.py','run_v47_kolf_consensus.py','run_v47_count_grid.py','run_v47_parallel_counts.py','test_v47_kolf_consensus.py','v39_kolf_hurdle.py','v38_response_hurdle.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
 inputs=['experiments/v39_kolf_support/analysis/collection.json','experiments/v39_kolf_support/analysis/response_features.npz','experiments/v45_bounded_dispersion/analysis/selection.json','experiments/v34_residual_flow/analysis/plan.json']
 for ctx in contexts():inputs.append('experiments/v37_component_ensemble/analysis/cache/'+ctx['family']+'/receipt.json')
 write(A/'plan.json',{'dependencies':{s:b.sha256(ROOT/s) for s in deps},'inputs':{s:b.sha256(ROOT/s) for s in inputs},'hypothesis':'The V39 learned feature head did not usefully exploit KOLF. Test direct independent experimental corroboration: only amplify a nontrivial V24 direction when KOLF agrees; never create directions, flip signs or convert z units to LFC. This preserves established response identity while strengthening independently supported effects.','coefficient':.25,'confidence':'tanh(abs(KOLF NTC-z)); finite nonsaturated coordinate, same sign as V24, abs(V24)>=.02. Gain1..1.25. Own target omitted.','negative_control':'Permute KOLF target rows only among273 supported targets, preserve missing27 target rows, same confidence formula/coefficient.','selection':'Means MSE/MAE<=1.01. All5 nulls and20 generated populations for genuine and matched shuffle. Original six raw guards, null no worse, >=2 raw metrics improve and positive oriented relative gain; genuine count gain must beat matched shuffled gain. Freeze passing genuine before24 additional count cases at2 seeds; same aggregate guards and study PDS/FID/REACH within.03.','scope':'Outcome-informed independent-data hypothesis on reused public development. KOLF iPSC related to H1, not claimed blind lineage. V46 count results not yet complete at freeze. Additional24 count cases and D3 still unopened. No automatic upload.'})
 print('V47 frozen',flush=True)

