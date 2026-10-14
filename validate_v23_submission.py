"""Full V23 package/axis/identity/gene-total audit in the scientific runtime."""
from pathlib import Path
import json,tempfile,subprocess,gc
import numpy as np,pandas as pd,anndata as ad,h5py
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
from release_v23_submission import ROOT,OUT,A,load,SEED
from v23_methods import emit,load_model,query_states
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model
import build_submission as b

def block(h,ptr,start,stop):
 lo,hi=int(ptr[start]),int(ptr[stop]);return sparse.csr_matrix((h['X/data'][lo:hi],h['X/indices'][lo:hi],ptr[start:stop+1]-lo),shape=(stop-start,18533))

def main():
 manifest,genes,targets,record,effects=load();ready=json.loads((A/'submission_ready.json').read_text());assert ready['sha256']==b.sha256(OUT/'prediction.vcc');assert ready['model_sha256']==b.sha256(A/'reference_model.json');prov=pd.read_csv(A/'prediction_control_rows.csv');assert len(prov)==360000;model=load_model('all');atlas,_=load_atlas_model();prior=ROOT/'experiments/v16_control_preserving';before=json.loads((prior/'analysis/independent_validation.json').read_text());assert before['passed'] and b.sha256(prior/'prediction.h5ad')==before['prediction_sha256']
 with tempfile.TemporaryDirectory(prefix='v23_payload_check_') as folder:
  payload=Path(folder)/'prediction.h5ad';run=subprocess.run([str(ROOT/'.competition_env/Scripts/python.exe'),str(ROOT/'unpack_v23_official_payload.py'),str(OUT/'prediction.vcc'),str(payload)],capture_output=True,text=True,encoding='utf-8',errors='replace');assert run.returncode==0,run.stderr[-2000:];metadata=json.loads(run.stdout);assert metadata==ready['package_meta'];obj=ad.read_h5ad(payload,backed='r');np.testing.assert_array_equal(obj.var_names,genes);np.testing.assert_array_equal(obj.obs.target_gene.astype(str),np.tile(np.repeat(targets,400),3));np.testing.assert_array_equal(obj.obs.context.astype(str),np.repeat(manifest['contexts'],120000));assert obj.obs_names.is_unique;obj.file.close();checked=0;replayed=0;groups=0;max_mean=0.
  with h5py.File(payload,'r') as packed,h5py.File(OUT/'prediction.h5ad','r') as original,h5py.File(prior/'prediction.h5ad','r') as baseline:
   ptr=packed['X/indptr'][:];oldptr=baseline['X/indptr'][:];np.testing.assert_array_equal(ptr,original['X/indptr'][:]);assert len(ptr)==360001 and ptr[-1]==metadata['nnz']
   for start in range(0,len(original['X/data']),4000000):
    stop=min(start+4000000,len(original['X/data']))
    for name in ['data','indices']:np.testing.assert_array_equal(packed['X/'+name][start:stop],original['X/'+name][start:stop])
    checked+=stop-start
   for ci,context in enumerate(manifest['contexts']):
    obj=ad.read_h5ad(ROOT/f'context_{context}.h5ad');np.testing.assert_array_equal(obj.var_names,genes);x=obj.X.tocsr().astype(np.int32);x.sum_duplicates();x.eliminate_zeros();x.sort_indices();library=np.asarray(x.sum(1,dtype=np.int64)).ravel();guides=obj.obs.ntc_id.astype(str).to_numpy();frame=prov[prov.context==context];rows=frame.source_control_row.to_numpy(int);np.testing.assert_array_equal(frame.source_control_cell_id,obj.obs_names.to_numpy()[rows]);np.testing.assert_array_equal(frame.ntc_id,guides[rows]);np.testing.assert_array_equal(frame.original_total_counts,library[rows]);depth,p,bulk,_=moments(x);states=query_states(x,model)
    for ti,target in enumerate(targets):
     start=(ci*300+ti)*400;count=block(packed,ptr,start,start+400);ref=block(baseline,oldptr,start,start+400);group=frame.iloc[ti*400:(ti+1)*400];rows=group.source_control_row.to_numpy(int);assert group.target_gene.eq(target).all() and len(set(rows))==400;chosen=b.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([SEED,ci,ti])));np.testing.assert_array_equal(rows,chosen);counts=np.unique(guides[rows],return_counts=True)[1];assert len(counts)==46 and counts.min()==8 and counts.max()==9;assert np.isfinite(count.data).all() and (count.data>0).all() and np.array_equal(count.data,np.floor(count.data));np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),library[rows]);np.testing.assert_array_equal(count.sum(0,dtype=np.int64),ref.sum(0,dtype=np.int64));drift=float(abs(np.asarray((count.astype(float)-ref.astype(float)).multiply((1/library[rows])[:,None]).mean(0)).ravel()).sum());max_mean=max(max_mean,drift);groups+=1
     if target in ['ABCD1','ACLY','VIM','MED25']:
      own=int(np.flatnonzero(genes==target)[0]);expected,_=emit(x[rows],effects[ci,ti],own,SEED+ci*100000+ti,'mean_preserved',(depth[rows],p,bulk),atlas,model,ti,states[rows]);assert (expected!=count).nnz==0;replayed+=1
    del x;gc.collect();print('V23 full counts, source identities, exact gene totals and replay passed',context,flush=True)
  assert checked==metadata['nnz'] and groups==900 and replayed==12
  b.write_json(A/'independent_validation.json',{'passed':True,'verified_rows':360000,'verified_groups':groups,'packaged_nonzero_values_checked':checked,'all_packaged_counts_indices_pointers_axes_match':True,'all_original_control_rows_ids_guides_libraries_verified':True,'all_gene_totals_exactly_equal_V16':True,'max_mean_CPM_probability_L1_change_from_V16':max_mean,'exact_replayed_groups':replayed,'scientific_runtime':'Python3.12.10 NumPy2.5.4 pandas3.0.6 SciPy1.18.1, BLAS2','official_validation_decompression_only':True,'package_sha256':ready['sha256'],'prediction_sha256':b.sha256(OUT/'prediction.h5ad'),'model_sha256':b.sha256(A/'reference_model.json'),'validator_sha256':b.sha256(Path(__file__)),'scope':'Full format/count/identity integrity and selected deterministic regeneration; no private accuracy guarantee.'})
 print('V23 complete payload independently verified',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
