"""V49 KL-balanced, bounded-rounding population ensembles."""
import numpy as np
from scipy import sparse
import v46_population_mixture as mixture
from v49_balanced_rounding import fit_expected,integerize
ROOT=mixture.ROOT;OUT=ROOT/'experiments/v49_balanced_population';A=OUT/'analysis'
ARMS=('balanced25','balanced75');OUTERS=mixture.OUTERS
b=mixture.b;data=mixture.data;read=mixture.read;write=mixture.write
def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,n
 return p
def project(count,reference,seed):
 x=count.toarray().astype(np.int32) if sparse.issparse(count) else np.asarray(count,np.int32);ref=np.asarray(reference,np.int32);libs=ref.sum(1,dtype=np.int64);order=np.argsort(libs,kind='stable')
 attempts=[]
 for bins in [32,64,128]:
  fixed=np.zeros_like(ref);moves=0;max_error=0.;max_steps=0
  for j,rows in enumerate(np.array_split(order,min(bins,len(order)))):
   initial=x[rows].astype(float)+.01*ref[rows];totals=ref[rows].sum(0,dtype=np.int64)
   e,steps=fit_expected(initial,libs[rows],totals);part,cost=integerize(e,libs[rows],totals,seed+j);fixed[rows]=part;moves+=cost;max_error=max(max_error,float(abs(part-e).max()));max_steps=max(max_steps,steps)
  error=float(abs((fixed/libs[:,None]).mean(0)-(ref/libs[:,None]).mean(0)).sum());attempts.append({'bins':bins,'mean_reference_L1':error})
  if error<=.0005:break
 assert error<=.0005,('Bounded projection misses CPM contract',attempts)
 np.testing.assert_array_equal(fixed.sum(1,dtype=np.int64),libs);np.testing.assert_array_equal(fixed.sum(0,dtype=np.int64),ref.sum(0,dtype=np.int64))
 return sparse.csr_matrix(fixed),{'mean_reference_L1':error,'bulk_exact':True,'depth_bins':bins,'projection_attempts':attempts,'binary_rounding_transfers':moves,'maximum_entry_rounding_error':max_error,'maximum_IPF_iterations':max_steps}
def emit(template,effect,own,seed,prepared,atlas,arm):
 if arm=='v24' or not np.any(effect):
  count,diag=mixture.emit(template,effect,own,seed,prepared,atlas,'v24');return count,{**diag,'zero_response_identity':not np.any(effect)}
 assert arm in ARMS
 m_arm='mixture25' if arm=='balanced25' else 'mixture75';count,diag=mixture.emit(template,effect,own,seed,prepared,atlas,m_arm);ref=mixture._CACHE['a'];count,extra=project(count,ref,seed+490049)
 libs,p,bulk=prepared;libs=np.asarray(libs,np.int64);factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum();rp=count.multiply(1/libs[:,None]).mean(0).A1;rb=count.sum(0).A1/libs.sum()
 return count,{**diag,**extra,'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'zero_response_identity':False}
def freeze():
 from run_v35_v36 import contexts
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v49_balanced_population.py','v49_balanced_rounding.py','run_v49_balanced_population.py','run_v49_count_grid.py','run_v49_parallel_counts.py','test_v49_balanced_rounding.py','test_v49_population.py','v46_population_mixture.py','v44_response_dispersion.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
 inputs=['experiments/v49_balanced_population/analysis/rounding_contract_tests.json','experiments/v48_margin_preserving/analysis/plan.json','experiments/v46_population_mixture/analysis/selection.json','experiments/v34_residual_flow/analysis/plan.json']
 for ctx in contexts():inputs.append('experiments/v37_component_ensemble/analysis/cache/'+ctx['family']+'/receipt.json')
 write(A/'plan.json',{'dependencies':{s:b.sha256(ROOT/s) for s in deps},'inputs':{s:b.sha256(ROOT/s) for s in inputs},'hypothesis':'V48 early CD4 fidelity loss suggests exact margins alone are insufficient: the previous greedy integer margin repair can concentrate corrections in individual cells. New numerical test confirms entry errors>1 count. Fit a smooth KL projection before rounding; alternating-path binary rounding guarantees each entry is floor/ceil of expectation, exact cell libraries and gene totals. Test moderate population mixtures25%/75%.','arms':ARMS,'projection':'Per-depth KL projection initialized mixture+.01*V24, maximum500 IPF steps; bounded binary rounding. Start32 depth bins, retry64/128 only if CPM L1>.0005. This numerical fallback uses predicted reference counts only, not true treatment outcomes.','selection':'All5 exact zero-response nulls and20 populations. Original six raw guards, >=2 raw improvements, positive oriented gain; every candidate pooled profile/PDS must equal V24; rounding error<=1 and CPM L1<=.0005. Lock highest passing arm before24 additional count cases at2 seeds; same guards and per-study PDS/FID/REACH within.03.','known_outcomes':'V38–V47 full development known; V48 partial development shows CD4 loss. Outcome-informed numerical correction, not blind discovery. Additional24 count cases and D3 unopened.','scope':'Reused public biology, unanchored raw metrics. No automatic upload or private score claim.'})
 print('V49 frozen',flush=True)

