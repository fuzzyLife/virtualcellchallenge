"""V46 distribution ensemble: retain whole cell count vectors from two noise models."""
import hashlib
import numpy as np
from scipy import sparse
import v44_response_dispersion as old
ROOT=old.ROOT;OUT=ROOT/'experiments/v46_population_mixture';A=OUT/'analysis'
ARMS=('mixture25','mixture50','mixture75');OUTERS=old.OUTERS
b=old.b;data=old.data;read=old.read;write=old.write
_CACHE={}
def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,n
 return p
def mix(a,bb,seed,fraction):
 assert a.shape==bb.shape
 rows=np.random.default_rng(seed+460046).permutation(a.shape[0])[:int(round(fraction*a.shape[0]))]
 out=a.copy();out[rows]=bb[rows];return out,rows
def emit(template,effect,own,seed,prepared,atlas,arm):
 if arm=='v24':return old.emit(template,effect,own,seed,prepared,atlas,'v24')
 assert arm in ARMS
 libs,p,bulk=prepared;libs=np.asarray(libs,np.int64)
 h=hashlib.sha256()
 for v in [template.data,template.indices,template.indptr,effect,libs,p,bulk]:h.update(np.ascontiguousarray(v).tobytes())
 key=(int(seed),own,h.hexdigest())
 if _CACHE.get('key')!=key:
  aa,_=old.emit(template,effect,own,seed,prepared,atlas,'v24');bb,_=old.emit(template,effect,own,seed,prepared,atlas,'residual20')
  _CACHE.clear();_CACHE.update(key=key,a=aa.toarray().astype(np.int32),b=bb.toarray().astype(np.int32))
 fraction={'mixture25':.25,'mixture50':.5,'mixture75':.75}[arm];out,rows=mix(_CACHE['a'],_CACHE['b'],seed,fraction)
 np.testing.assert_array_equal(out.sum(1),libs);assert out.min()>=0
 factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum();rp=(out/libs[:,None]).mean(0);rb=out.sum(0)/libs.sum()
 return sparse.csr_matrix(out),{'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'changed_transcript_fraction':float(abs(out-template.toarray()).sum()/2/libs.sum()),'empirical20_component_fraction':fraction,'component20_cells':len(rows),'whole_cell_vectors_preserved':True}
def freeze():
 from run_v35_v36 import contexts
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v46_population_mixture.py','run_v46_population_mixture.py','run_v46_count_grid.py','run_v46_parallel_counts.py','test_v46_population_mixture.py','v44_response_dispersion.py','v38_response_hurdle.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
 inputs=['experiments/v44_response_dispersion/analysis/selection.json','experiments/v45_bounded_dispersion/analysis/selection.json','experiments/v34_residual_flow/analysis/plan.json']
 for ctx in contexts():inputs.append('experiments/v37_component_ensemble/analysis/cache/'+ctx['family']+'/receipt.json')
 write(A/'plan.json',{'dependencies':{s:b.sha256(ROOT/s) for s in deps},'inputs':{s:b.sha256(ROOT/s) for s in inputs},'hypothesis':'V44 residual20 improves REACH16.8% but misses JAC guard. V45 smaller uniform increments pass loss guards but average gain is negative. Instead of averaging parameters or count values, mix whole cell vectors from V24 and residual20. This models coexistence of compact and more variable populations with heavier tails; preserves exact integer libraries and each component cell vector.','arms':ARMS,'fractions':[.25,.5,.75],'assignment':'Seed+460046 balanced nested permutation of rows, independent of expression/labels. Means/effects/controls shared; no extra final projection or measurement resampling.','selection':'All5 nulls and20 treatment cases, original six metric/null guards, >=2 metrics improve, positive mean oriented relative gain. Lock largest qualifying gain before24 additional count cases at2 seeds; same aggregate guards and study PDS/FID/REACH within.03. If no candidate, retain failure and continue separately.','known_outcomes':'Explicit outcome-informed next hypothesis after V38–V45. Additional24 count cases and D3 still unopened at freeze.','scope':'Reused public development, not private anchor scores or blind lineage validation. No automatic upload.'})
 print('V46 frozen',flush=True)

