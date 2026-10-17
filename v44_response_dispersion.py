"""V24 means with moderate empirical residual weight, retaining compact sampling."""
import numpy as np
from scipy import sparse
import v38_response_hurdle as parent
ROOT=parent.ROOT;OUT=ROOT/'experiments/v44_response_dispersion';A=OUT/'analysis'
ARMS=('residual20','residual35','adaptive');OUTERS=parent.OUTERS
b=parent.b;data=parent.data;read=parent.read;write=parent.write
def check():
 p=read(A/'plan.json')
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,n
 return p
def residual_weight(effect,arm):
 if arm=='residual20':return np.full(len(effect),.20)
 if arm=='residual35':return np.full(len(effect),.35)
 if arm=='adaptive':return .1+.25*.15/(.15+abs(effect))
 raise ValueError(arm)
def initial(expected,poisson,empirical,weight):
 return (.35-weight)[None,:]*expected+.65*poisson+weight[None,:]*empirical
def emit(template,effect,own,seed,prepared,atlas,arm):
 if arm=='v24':return parent.emit(template,effect,own,seed,prepared,atlas)
 libs,p,bulk=prepared;libs=np.asarray(libs,np.int64);factor=np.exp2(np.clip(effect,-3,3))
 if own is not None:factor[own]=.1
 wp=p*factor;wp/=wp.sum();wb=bulk*factor;wb/=wb.sum()
 expected=libs[:,None]*wp;emp=template.toarray().astype(float)*factor;emp*=libs[:,None]/emp.sum(1,keepdims=True)
 weight=residual_weight(effect,arm)
 if own is not None:weight[own]=.1
 q=initial(expected,np.random.default_rng(seed).poisson(expected),emp,weight);assert (q>=0).all();q/=q.sum(1,keepdims=True)
 count=atlas.dual_moment_counts(q,wp,wb,depths=libs,seed=int(seed),iterations=200)
 np.testing.assert_array_equal(count.sum(1),libs);assert count.dtype.kind in 'iu' and count.min()>=0
 rp=(count/libs[:,None]).mean(0);rb=count.sum(0)/libs.sum()
 return sparse.csr_matrix(count),{'mean_probability_L1':float(abs(rp-wp).sum()),'bulk_probability_L1':float(abs(rb-wb).sum()),'empirical_weight_min':float(weight.min()),'empirical_weight_max':float(weight.max()),'empirical_weight_mean':float(weight.mean()),'changed_transcript_fraction':float(abs(count-template.toarray()).sum()/2/libs.sum()),'response_means_unchanged':True}
def freeze():
 from run_v35_v36 import contexts
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 deps=['v44_response_dispersion.py','run_v44_response_dispersion.py','run_v44_count_grid.py','test_v44_response_dispersion.py','v38_response_hurdle.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py','evaluate_atlas_transfer.py']
 inputs=['experiments/v34_residual_flow/analysis/plan.json','experiments/v38_response_hurdle/analysis/development_results.json','experiments/v42_identity_trust/analysis/selection.json']
 for c in contexts():inputs.append('experiments/v37_component_ensemble/analysis/cache/'+c['family']+'/receipt.json')
 write(A/'plan.json',{'dependencies':{s:b.sha256(ROOT/s) for s in deps},'inputs':{s:b.sha256(ROOT/s) for s in inputs},'hypothesis':'V10->V16 benefited from adding empirical residual variation. Previously tested reduced empirical weights, learned full transport and full empirical projection failed. Test moderate increased empirical residual weights while retaining the same Poisson draw and dual moment fit, V24 means, controls and exact libraries. Strong predicted effects receive less empirical weight in adaptive arm.','arms':ARMS,'formula':'(.35-alpha_j)*expected+.65*Poisson(expected)+alpha_j*row-normalized empirical perturbation. alpha=.2,.35 or .1+.25*.15/(.15+abs(LFC)); own target retains.1. Nonnegative alpha<=.35.','selection':'Evaluate5 real control-null and20 treatment cases for all arms, six raw metrics, original V38 guards, at least2 raw improvements, positive oriented relative gain. Highest passing arm locked before24 additional count cases at2 seeds. Same aggregate guards and per-study PDS/FID/REACH within.03. If none qualify, continue separately without claiming improvement.','mean_inputs':'Exact frozen V24 effects; no treatment labels used to infer dispersion. All generated counts still re-fitted to same desired mean/pooled profile.','scope':'Reused public development biology; raw metrics not official anchor scores. Additional24 count cases/D3 unopened at freeze. Previous failures known. No automatic upload.','private_score_known':False})
 print('V44 frozen',flush=True)

