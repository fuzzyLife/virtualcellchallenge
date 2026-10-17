"""V48: depth-stratified integer-margin projection of a population ensemble."""
import numpy as np
from scipy import sparse
import v46_population_mixture as mixture
from v26_methods import preserve_depth_margins
ROOT=mixture.ROOT;OUT=ROOT/'experiments/v48_margin_preserving';A=OUT/'analysis'
ARMS=('projected50','projected75');OUTERS=mixture.OUTERS
b=mixture.b;data=mixture.data;read=mixture.read;write=mixture.write
def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,n
 return p
def project(count,reference,seed):
 x=count.toarray().astype(np.int32) if sparse.issparse(count) else np.asarray(count,np.int32)
 ref=reference.toarray().astype(np.int32) if sparse.issparse(reference) else np.asarray(reference,np.int32)
 fixed,moved=preserve_depth_margins(x,ref,seed,bins=16);libs=ref.sum(1,dtype=np.int64)
 error=float(abs((fixed/libs[:,None]).mean(0)-(ref/libs[:,None]).mean(0)).sum())
 if error>.0005:
  from v29_pairwise_repair import repair_pairwise
  fixed,repair=repair_pairwise(fixed,ref);error=repair['mean_reference_L1']
 else:repair={'CPM_repair_molecules':0}
 np.testing.assert_array_equal(fixed.sum(0,dtype=np.int64),ref.sum(0,dtype=np.int64));np.testing.assert_array_equal(fixed.sum(1,dtype=np.int64),libs)
 assert error<=.0005,('Mean preservation failed',error)
 return sparse.csr_matrix(fixed),{'mean_reference_L1':error,'bulk_exact':True,'depth_margin_repair_molecules':moved,'CPM_repair_molecules':repair['CPM_repair_molecules']}
def emit(template,effect,own,seed,prepared,atlas,arm):
 if arm=='v24' or not np.any(effect):
  count,diag=mixture.emit(template,effect,own,seed,prepared,atlas,'v24');return count,{**diag,'zero_response_identity':not np.any(effect)}
 assert arm in ARMS
 m_arm='mixture50' if arm=='projected50' else 'mixture75';count,diag=mixture.emit(template,effect,own,seed,prepared,atlas,m_arm)
 ref=mixture._CACHE['a'];count,extra=project(count,ref,seed+480048);libs,p,bulk=prepared;libs=np.asarray(libs,np.int64)
 factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum();rp=count.multiply(1/libs[:,None]).mean(0).A1;rb=count.sum(0).A1/libs.sum()
 return count,{**diag,**extra,'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'zero_response_identity':False}
def freeze():
 from run_v35_v36 import contexts
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v48_margin_preserving.py','run_v48_margin_preserving.py','run_v48_count_grid.py','run_v48_parallel_counts.py','test_v48_margin_preserving.py','v46_population_mixture.py','v44_response_dispersion.py','v26_methods.py','v29_pairwise_repair.py','v29_swap_kernel.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
 inputs=['experiments/v46_population_mixture/analysis/selection.json','experiments/v34_residual_flow/analysis/plan.json']
 for ctx in contexts():inputs.append('experiments/v37_component_ensemble/analysis/cache/'+ctx['family']+'/receipt.json')
 write(A/'plan.json',{'dependencies':{s:b.sha256(ROOT/s) for s in deps},'inputs':{s:b.sha256(ROOT/s) for s in inputs},'hypothesis':'V46 improves coverage but random population assignment drifts pooled margins, reducing PDS/expression. Project moderate compact/residual20 mixtures to exact V24 gene totals and exact cell libraries, with16 depth-stratum margins and CPM error<=.0005. This isolates population variation; pooled PDS must equal V24.','arms':ARMS,'fractions':[.5,.75],'zero_response_boundary':'If downstream effect vector is exactly zero, return exact V24 counts. Applies to all zero-response treatments as well as control-null cases; do not add population uncertainty absent a response. Boundary added explicitly after V46 null drift.','selection':'All5 exact identity nulls and20 treatment cases with original six raw guards, >=2 metrics improve and positive relative oriented gain. Require observed pooled profiles/PDS exactly V24 and mean-reference L1<=.0005. Lock one genuine candidate before24 additional count cases at2 seeds; same aggregate guards and study PDS/FID/REACH within.03.','scope':'Outcome-informed method after V46 failure; V47 count results incomplete at freeze. Reused public biology, not private anchored scoring. Additional24 count cases and D3 still unopened. No automatic upload.'})
 print('V48 frozen',flush=True)

