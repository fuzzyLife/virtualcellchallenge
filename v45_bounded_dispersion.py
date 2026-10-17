"""V45 bounded dispersion increments after the recorded V44 tradeoff."""
import numpy as np
import v44_response_dispersion as old
ROOT=old.ROOT;OUT=ROOT/'experiments/v45_bounded_dispersion';A=OUT/'analysis'
ARMS=('residual15','response_weighted20');OUTERS=old.OUTERS
b=old.b;data=old.data;read=old.read;write=old.write
def residual_weight(effect,arm):
 if arm=='residual15':return np.full(len(effect),.15)
 if arm=='response_weighted20':return .1+.1*np.minimum(abs(effect)/.2,1)
 raise ValueError(arm)
def emit(template,effect,own,seed,prepared,atlas,arm):
 old.residual_weight=residual_weight
 return old.emit(template,effect,own,seed,prepared,atlas,arm)
def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,n
 return p
def freeze():
 from run_v35_v36 import contexts
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v45_bounded_dispersion.py','run_v45_bounded_dispersion.py','run_v45_count_grid.py','run_v45_parallel_counts.py','test_v45_bounded_dispersion.py','v44_response_dispersion.py','v38_response_hurdle.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
 inputs=['experiments/v44_response_dispersion/analysis/selection.json','experiments/v34_residual_flow/analysis/plan.json']
 for c in contexts():inputs.append('experiments/v37_component_ensemble/analysis/cache/'+c['family']+'/receipt.json')
 write(A/'plan.json',{'dependencies':{s:b.sha256(ROOT/s) for s in deps},'inputs':{s:b.sha256(ROOT/s) for s in inputs},'hypothesis':'V44 residual20 increases REACH16.8%, keeps PDS/expression/FID/NMAE guards, but JAC loss.001432 exceeds.001 guard. Larger residual35 loses CD4 fidelity. Test smaller15% empirical residual, and a response-strength-weighted10–20% increment that keeps weak/zero responses at original compact recipe. This is explicit outcome-informed refinement, not independent discovery.','arms':ARMS,'coefficients':'residual15=.15; response_weighted20=.1+.1*min(abs(predicted LFC)/.2,1), own=.1. 65% same Poisson component; deterministic expectation compensates.','selection':'All5 real split-control nulls and20 treatment cases, original V38 six raw guards/null no worse, at least2 improved metrics and positive average oriented relative gain. Lock highest passing new arm before24 additional count cases at two seeds. Same aggregate guards plus study PDS/FID/REACH within.03. All means bitwise V24.','known_outcomes':'V38–V44 public development examined. Additional24 generated-count cases and D3 remain unopened at freeze.','scope':'Reused public development and unscaled scores, not exact private anchors or blind context confirmation. No automated upload.'})
 print('V45 frozen',flush=True)

