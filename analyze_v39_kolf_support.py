"""Audit the optional KOLF response features without changing V38 predictions."""
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
import build_submission as b
from fit_v7_response import read_cache
def main():
 a=b.ROOT/'experiments/v39_kolf_support/analysis';rec=json.loads((a/'collection.json').read_text(encoding='utf-8'));path=a/'response_features.npz';assert b.sha256(path)==rec['arrays_sha256'];d=read_cache()
 with np.load(path,allow_pickle=False) as z:
  np.testing.assert_array_equal(z['genes'],d['genes']);np.testing.assert_array_equal(z['targets'],d['targets']);response=z['z'].copy();present=z['target_present'].copy()
 permutation=np.random.default_rng(39039).permutation(np.flatnonzero(present));lookup={i:j for i,j in zip(np.flatnonzero(present),permutation)};rows=[]
 for ri,ti in enumerate(d['row_targets']):
  if not present[ti]:continue
  truth=d['truth'][ri];base=d['valid'][ri]&np.isfinite(d['profiles'][d['row_contexts'][ri]])&(d['profiles'][d['row_contexts'][ri]]>5)
  own=int(np.flatnonzero(d['genes']==d['targets'][ti])[0]);base[own]=False
  result={'row':ri,'context':int(d['row_contexts'][ri]),'target':str(d['targets'][ti]),'family':str(d['row_families'][ri])}
  for name,index in [('identity',ti),('target_shuffled',lookup[ti])]:
   mask=base&np.isfinite(response[index]);xx=response[index,mask];yy=truth[mask]
   result[name+'_pairs']=int(mask.sum());result[name+'_spearman']=float(spearmanr(xx,yy).statistic) if mask.sum()>100 and np.ptp(xx)>0 and np.ptp(yy)>0 else None
   active=mask&(abs(truth)>=.25)&(abs(truth)>=2*d['errors'][ri])&(abs(response[index])>=.12)
   result[name+'_active_pairs']=int(active.sum());result[name+'_active_sign_agreement']=float(np.mean(np.sign(response[index,active])==np.sign(truth[active]))) if active.any() else None
  rows.append(result)
 summary={}
 for f in sorted(set(r['family'] for r in rows)):
  rr=[r for r in rows if r['family']==f];summary[f]={}
  for k in ['identity_spearman','target_shuffled_spearman','identity_active_sign_agreement','target_shuffled_active_sign_agreement']:
   v=[r[k] for r in rr if r[k] is not None and np.isfinite(r[k])];summary[f][k]=float(np.mean(v)) if v else None
  summary[f]['target_rows']=len(rr)
 file=a/'source_agreement.json';assert not file.exists();b.write_json(file,{'records':rows,'summary':summary,'collection_sha256':b.sha256(a/'collection.json'),'script_sha256':b.sha256(Path(__file__)),'scope':'Exploratory matched-intervention agreement in existing public summaries. Quantized z scores are auxiliary features, not raw counts or LFC. Whole-cell held-query inference and competition performance have not been tested for this new source. V38 does not use these data.'});print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()

