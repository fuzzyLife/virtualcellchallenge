"""V50 output-gene response propensity, separating generic calibration from target identity."""
import numpy as np
import v47_kolf_consensus as source
ROOT=source.ROOT;OUT=ROOT/'experiments/v50_response_prior';A=OUT/'analysis'
ARMS=('mobility_prior','randomized_prior');OUTERS=source.OUTERS
b=source.b;data=source.data;read=source.read;write=source.write;emit=source.emit
def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,n
 return p
def propensity(z,sat):
 valid=np.isfinite(z)&~sat;score=np.where(valid,np.tanh(abs(np.nan_to_num(z))),0)
 return {'valid':valid,'score':score,'den':valid.sum(0),'up':(score*(z>0)).sum(0),'down':(score*(z<0)).sum(0)}
def averaged(base,z,sat):
 bank=propensity(z,sat);out=base.copy()
 for ti in range(len(base)):
  den=bank['den']-bank['valid'][ti];up=bank['up']-bank['score'][ti]*(z[ti]>0);down=bank['down']-bank['score'][ti]*(z[ti]<0)
  confidence=np.divide(np.where(base[ti]>0,up,down),den,out=np.zeros(base.shape[1]),where=den>0)
  confidence=np.where(abs(base[ti])>=.02,confidence,0);out[ti]=np.clip(base[ti]*(1+.25*confidence),-3,3)
 return out.astype(np.float32)
def predict(base):
 z,sat,_,_=source.evidence();out={'mobility_prior':averaged(base,z,sat),'randomized_prior':source.predict(base)['shuffled_consensus25']}
 d=data()
 for a,v in out.items():
  for ti,t in enumerate(d['targets']):v[ti,d['genes']==t]=0
 return out
def freeze():
 from run_v35_v36 import contexts
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v50_response_prior.py','run_v50_response_prior.py','run_v50_count_grid.py','run_v50_parallel_counts.py','test_v50_response_prior.py','v47_kolf_consensus.py','v39_kolf_hurdle.py','v38_response_hurdle.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
 inputs=['experiments/v47_kolf_consensus/analysis/selection.json','experiments/v47_kolf_consensus/analysis/development_results.json','experiments/v39_kolf_support/analysis/collection.json','experiments/v39_kolf_support/analysis/response_features.npz','experiments/v34_residual_flow/analysis/plan.json']
 for ctx in contexts():inputs.append('experiments/v37_component_ensemble/analysis/cache/'+ctx['family']+'/receipt.json')
 write(A/'plan.json',{'dependencies':{s:b.sha256(ROOT/s) for s in deps},'inputs':{s:b.sha256(ROOT/s) for s in inputs},'hypothesis':'V47 shuffled target rows pass ordinary count guards while true target-matched corroboration fails. Reinterpret explicitly as a generic response-amplitude prior, not useful target identity. Average independent KOLF signed response propensity across interventions per output gene, leave current target intervention out. Compare with the fixed previously observed randomized checkpoint; never attribute its gain to target-specific biology.','arms':ARMS,'prior':'Gain1+.25*mean(tanh(abs(z))*same direction), excluding missing/saturated source coordinates and own intervention row; abs(V24)>=.02. Weak/zero effects unchanged, own target zero. z is confidence only, not converted to LFC. Clipping means this is an averaged-prior approximation, not an exact conditional-expectation claim at clipping bounds.','randomized_checkpoint':'Exact fixed V47 shuffled_consensus25 seed47047 predictions and standard emitter. Its development gain1.575% already known; its reuse is explicitly outcome-informed. Fresh generated-count confirmation is still unopened.','selection':'Mean MSE/MAE<=1.01. All5 zero-response nulls,20 development populations, original six raw count guards, >=2 raw metrics improve and positive oriented average relative gain. Choose highest passing gain, with no requirement for biological target identity because this hypothesis is generic calibration. Lock once before24 additional count cases at2 seeds. Same aggregate guards and study PDS/FID/REACH within.03. No reselection after confirmation.','known_outcomes':'V38–V48 development known; V49 count tests running. This is a new interpretation after controls, not a preregistered positive biology finding. Additional24 count cases and D3 unopened.','private_score_known':False,'automatic_upload':False})
 print('V50 generic-prior hypothesis frozen',flush=True)

