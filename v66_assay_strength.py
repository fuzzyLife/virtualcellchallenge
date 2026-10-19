"""V66: strength-matched original-assay control for V65, before count-based selection."""
from pathlib import Path
import argparse,gc,json
import numpy as np
from threadpoolctl import threadpool_limits
import v62_conserved_response as old
from v63_identity_guard import guard
from v62_pool_cache import install
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'experiments/v66_assay_strength';A=OUT/'analysis';ARMS=('assay75',)

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v66_assay_strength.py','run_v66_experiment.py','run_v65_experiment.py','v62_conserved_response.py','v63_identity_guard.py','v62_pool_cache.py','v16_methods.py']
 plan=read(ROOT/'experiments/v65_closed_response/analysis/plan.json')
 write(A/'plan.json',dict(version=66,dependencies={n:sha(ROOT/n) for n in deps},arms=ARMS,hypothesis='Strength-matched control motivated by V65 response screen: reclosure50 and closed50 lose to V63, while closed75 slightly improves profile error. Test the same three-quarter correction on the original V63 source representation, to separate correction strength from normalization. V65 counts started; only progress/runtime logs inspected, no new count summary or accuracy used to specify this candidate.',formula='Take the exact original unguarded V62 assay_center50 array, expand its correction from V24 by 1.5 (bounded +/-0.75 and final +/-3); preserve own zero; apply exactly the original V24 80% identity-margin guard. Original source profiles and uncertainties unchanged. It is strength matched to V65 up to original float32 rounding.',selection=plan['selection'],count_guards=plan['count_guards'],validation=plan['validation'],V65_mean_result_sha256=sha(ROOT/'experiments/v65_closed_response/analysis/mean_results.json'),cross_version_rule='Only V65 closed75 and this fixed V66 assay75 can be finalists. Evaluate the same canonical 30 target/context pairs and two seeds, V24/V63/shuffled controls. Prefer a candidate passing every frozen public count guard; if both pass choose lower paired equal-study generated-profile error ratio. If neither passes, any diagnostic upload must remain explicitly experimental and report all failures; prefer lower generated-profile error among mean-eligible candidates with no >3% study profile degradation and at least two raw metric gains versus V63. Do not promote additional strengths after count outcomes. No hidden-score forecast.',D3_opened=False,automatic_upload=False))
 print('V66 strength-matched control frozen',sha(A/'plan.json'),flush=True)

def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert sha(ROOT/n)==h,n
 assert sha(ROOT/'experiments/v65_closed_response/analysis/mean_results.json')==p['V65_mean_result_sha256'];return p

def predict(query,baseline,bulk,excluded=(),shuffle=False):
 install();pred=old.predict(query,baseline,excluded,shuffle)['assay_center50'];proposal=np.clip(np.asarray(baseline,float)+1.5*(pred-np.asarray(baseline,float)),-3,3).astype(np.float32)
 d=old.data();own=np.array([list(d['genes']).index(t) for t in d['targets']]);panel=np.isin(d['genes'],d['targets']);proposal[np.arange(300),own]=0;e,trace=guard(baseline,proposal,bulk,own,panel)
 return {'assay75':e},{'assay75':trace}

def prepare():
 check();from run_v35_v36 import contexts,key
 from run_v37_ensemble import checked_folder
 for ctx in contexts():
  name=ctx['family'];path=A/'predictions'/(name+'.npz');path.parent.mkdir(exist_ok=True)
  if path.exists():assert sha(path)==read(path.with_suffix('.json'))['sha256'];continue
  folder,_=checked_folder(name)
  with np.load(folder/'effects.npz',allow_pickle=False) as z:base=z['v24'].copy();q=z['query'].copy();bulk=z['bulk'].copy();measured=z['measured'].copy()
  pred,traces=predict(q,base,bulk,(key(name),));sh,shtrace=predict(q,base,bulk,(key(name),),True)
  vp=ROOT/'experiments/v63_identity_guard/analysis/predictions'/(name+'.npz');assert sha(vp)==read(vp.with_suffix('.json'))['sha256']
  with np.load(vp,allow_pickle=False) as z:v63=z['all_effects'].copy()
  # Replay original response strength to verify the comparator's mathematical path.
  unguarded=old.predict(q,base,(key(name),))['assay_center50'];own=np.array([list(old.data()['genes']).index(t) for t in old.data()['targets']]);panel=np.isin(old.data()['genes'],old.data()['targets']);oldguard,_=guard(base,unguarded,bulk,own,panel);np.testing.assert_array_equal(v63,oldguard)
  np.savez_compressed(path,query=q,bulk=bulk,measured=measured,v24=base,v63=v63,**pred,**{k+'_shuffled':v for k,v in sh.items()});write(path.with_suffix('.json'),dict(sha256=sha(path),traces=traces,shuffled_traces=shtrace,v63_exactly_replayed=True,source_receipt_sha256=sha(folder/'receipt.json'),excluded_study=key(name)));old.bank.cache_clear();gc.collect();print('V66 prepared',name,flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare']);a=p.parse_args()
 with threadpool_limits(limits=2):globals()[a.stage]()
