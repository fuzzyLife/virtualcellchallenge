"""Cross-assay experimental phenotype features for source-specific bilinear ridge."""
import json,shutil
from functools import lru_cache
import numpy as np
import v30_methods as core
import v38_response_hurdle as parent
ROOT=core.ROOT;OUT=ROOT/'experiments/v43_lincs_response';A=OUT/'analysis'
OUTERS=core.OUTERS;ARMS=('fingerprint125','fingerprint250','shuffled250')
b=core.b;data=parent.data;read=parent.read;write=parent.write;emit=parent.emit
def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,n
 return p
def embeddings():
 path=A/'experimental_fingerprints.npz';assert b.sha256(path)==read(A/'fingerprints.json')['arrays_sha256']
 with np.load(path,allow_pickle=False) as z:return z['keys'],z['vectors'].astype(float)
def install():
 core.A=A;core.embeddings=embeddings;core.load_plan=check;core.SEED=43043
def prepare_sources():
 check();dest=A/'sources';dest.mkdir(exist_ok=True);records=[]
 for name in core.SOURCES:
  src=ROOT/'experiments/v30_structured_response/analysis/sources'/(name+'.npz');rec=read(src.with_suffix('.json'));assert b.sha256(src)==rec['sha256']
  dst=dest/src.name
  if not dst.exists():shutil.copyfile(src,dst)
  assert b.sha256(dst)==rec['sha256'];write(dst.with_suffix('.json'),rec);records.append({'source':name,'sha256':rec['sha256'],'original_receipt_sha256':b.sha256(src.with_suffix('.json'))})
 write(A/'source_reuse.json',{'sources':records,'original_numerical_sources_unchanged':True})
def fit(name,variant):
 install();core.fit_one(name,variant)
def predict(query,excluded,base):
 install();pred={};trace={}
 for arm in ARMS:
  variant='random_full' if arm=='shuffled250' else 'string_full';alpha=.125 if arm=='fingerprint125' else .25
  p,_,t=core.predict(query,excluded,variant,(32,10.,alpha),parent=base);pred[arm]=p;trace[arm]=t
 return pred,trace
def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v43_lincs_fingerprints.py','prepare_v43_fingerprints.py','unpack_v43_lincs.py','run_v43_fingerprints.py','test_v43_fingerprints.py','v30_methods.py','run_v42_count_grid.py','v38_response_hurdle.py','v25_methods.py','v24_methods.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
 inputs=['experiments/v43_lincs_response/analysis/collection.json','experiments/v43_lincs_response/analysis/unpacked.json','experiments/v43_lincs_response/analysis/fingerprints.json','experiments/v43_lincs_response/analysis/experimental_fingerprints.npz','experiments/v34_residual_flow/analysis/plan.json']
 for name in core.SOURCES:
  for ext in ['npz','json']:inputs.append('experiments/v30_structured_response/analysis/sources/'+name+'.'+ext)
 for ctx in ['H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2']:inputs.append('experiments/v37_component_ensemble/analysis/cache/'+ctx+'/receipt.json')
 write(A/'plan.json',{'dependencies':{p:b.sha256(ROOT/p) for p in deps},'inputs':{p:b.sha256(ROOT/p) for p in inputs},'hypothesis':'Measured experimental RNAi fingerprints can transfer functional perturbation identity across assays better than static STRING annotations. Fit source-only CRISPR response bases and map phenotype fingerprints into them, excluding whole held study from pooling.','rank':32,'ridge_lambda':10,'blends':[.125,.25],'correction_cap_LFC':.5,'arms':ARMS,'negative_control':'Same model and source targets, fingerprints permuted among gene keys with fixed seed43043. Not a permutation of outputs.','selection':'Mean MSE/MAE each <=1.01, genuine beats matched shuffled250. Test both genuine arms and shuffled250 on20 count development cases. Require all original six raw count/null guards, positive average gain, at least2 metrics improved. Lock highest passing genuine candidate before24 remaining count cases, seeds0/+380000. Same aggregate guards and study PDS/FID/REACH within.03. If no genuine passes, record failure and continue with a separate hypothesis; do not fabricate a better score.','scope':'Public development repeatedly used, not blind cell lineage validation or exact anchored leaderboard replication. 24 additional count cases and D3 unopened at freeze. Missing fingerprints and measured outputs abstain to V24. RNAi off-target effects remain a limitation.','published_champion':'V24 overall0.14299785521745675','private_score_known':False})
 print('V43 plan frozen',flush=True)

