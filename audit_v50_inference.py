"""Independent complete replay of the selected production response-prior formula."""
from pathlib import Path
import numpy as np,json
import build_submission as b
ROOT=b.ROOT;A=ROOT/'experiments/v50_response_prior/analysis'
def main():
 rec=json.loads((A/'reference_model.json').read_text());assert rec['selected_arm']=='randomized_prior';manifest,genes,targets=b.read_spec(ROOT)
 parent=ROOT/'experiments/v24_cd4_response/analysis/reference_model_arrays.npz';prior=ROOT/'experiments/v39_kolf_support/analysis/response_features.npz';model=A/'reference_model_arrays.npz'
 assert b.sha256(parent)==rec['parent_V24_arrays_sha256'] and b.sha256(prior)==rec['independent_prior_arrays_sha256'] and b.sha256(model)==rec['arrays_sha256']
 with np.load(parent,allow_pickle=False) as f:base=f['effects'].copy()
 with np.load(prior,allow_pickle=False) as f:z=f['z'].copy();sat=f['saturated'].copy();np.testing.assert_array_equal(f['genes'],genes);np.testing.assert_array_equal(f['targets'],targets)
 with np.load(model,allow_pickle=False) as f:out=f['effects'].copy();np.testing.assert_array_equal(f['genes'],genes);np.testing.assert_array_equal(f['targets'],targets)
 ids=np.flatnonzero(np.any(np.isfinite(z),axis=1));mapping=np.arange(len(targets));mapping[ids]=ids[np.random.default_rng(47047).permutation(len(ids))]
 checked=0
 for ci in range(3):
  for ti,t in enumerate(targets):
   v=base[ci,ti];row=z[mapping[ti]];valid=np.isfinite(row)&~sat[mapping[ti]]&(np.sign(row)==np.sign(v))&(abs(v)>=.02)
   expected=v.copy();expected[valid]=np.clip(v[valid]*(1+.25*np.tanh(abs(row[valid]))),-3,3);expected[genes==t]=0
   np.testing.assert_array_equal(expected,out[ci,ti]);assert np.equal(np.sign(out[ci,ti]),np.sign(v)).all();assert (abs(out[ci,ti])<=1.25*abs(v)+1e-6).all();checked+=len(genes)
 b.write_json(A/'production_inference_audit.json',{'passed':True,'contexts':3,'targets_per_context':300,'coordinates_replayed':checked,'full_exact_array_comparison':True,'signs_zeros_and_gain_bounds_checked':True,'source_mapping':'Seed47047 permutation among only273 source-supported targets, independently rebuilt. No model predict/combine function used.','model_arrays_sha256':b.sha256(model),'script_sha256':b.sha256(Path(__file__))})
 print('Independent production inference PASS',checked,'coordinates',flush=True)
if __name__=='__main__':main()

