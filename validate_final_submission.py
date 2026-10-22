"""Independently compare every packaged Final V91 count, source identity and twelve regenerated groups."""
from pathlib import Path
import tempfile,subprocess,gc
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from build_final_submission import ROOT,OUT,A,load,SEED
from v57_pie import read,write,sha
from v16_methods import emit
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model

def main():
 manifest,genes,targets,model,effects=load();ready=read(A/'submission_ready.json');assert sha(OUT/'prediction.vcc')==ready['sha256'];assert (OUT/'prediction.vcc').stat().st_size==ready['bytes'];assert ready['model_sha256']==sha(A/'reference_model.json');assert np.isfinite(effects).all();own=np.array([list(genes).index(t) for t in targets]);assert (effects[:,np.arange(300),own]==0).all()
 provenance=pd.read_csv(A/'prediction_control_rows.csv');stats=pd.read_csv(A/'prediction_group_statistics.csv');assert len(provenance)==360000 and len(stats)==900;atlas,_=load_atlas_model();replays=[];scanned=0;compared=0
 with tempfile.TemporaryDirectory(prefix='final_v91_payload_check_') as tmp:
  folder=Path(tmp).resolve();assert folder.parent==Path(tempfile.gettempdir()).resolve() and folder.name.startswith('final_v91_payload_check_');payload=folder/'prediction.h5ad';result=subprocess.run([str(ROOT/'.competition_env/Scripts/python.exe'),str(ROOT/'unpack_v24_official_payload.py'),str(OUT/'prediction.vcc'),str(payload)],capture_output=True,text=True,encoding='utf-8',errors='replace');assert result.returncode==0,result.stderr[-2000:];meta=__import__('json').loads(result.stdout);assert meta==ready['package_meta'];obj=ad.read_h5ad(payload,backed='r');np.testing.assert_array_equal(obj.var_names,genes);np.testing.assert_array_equal(obj.obs.context.astype(str),np.repeat(manifest['contexts'],120000));np.testing.assert_array_equal(obj.obs.target_gene.astype(str),np.tile(np.repeat(targets,400),3));assert obj.obs_names.is_unique;obj.file.close()
  with h5py.File(payload,'r') as packed,h5py.File(OUT/'prediction.h5ad','r') as original:
   ptr=packed['X/indptr'][:];np.testing.assert_array_equal(ptr,original['X/indptr'][:]);assert ptr.shape==(360001,) and ptr[-1]==meta['nnz']
   for start in range(0,int(ptr[-1]),4000000):
    end=min(start+4000000,int(ptr[-1]))
    for name in ['data','indices']:np.testing.assert_array_equal(packed['X/'+name][start:end],original['X/'+name][start:end])
    compared+=end-start
   for ci,c in enumerate(manifest['contexts']):
    x,guides,library,*_=build.load_control(ROOT/'final',c,genes,manifest);depth,cp,bulk,_=moments(x);obj=ad.read_h5ad(ROOT/'final'/('context_'+c+'.h5ad'),backed='r');ids=obj.obs_names.to_numpy();obj.file.close();prov=provenance[provenance.context==c];rows=prov.source_control_row.to_numpy(int);np.testing.assert_array_equal(prov.source_control_cell_id.to_numpy(str),ids[rows].astype(str));np.testing.assert_array_equal(prov.ntc_id.to_numpy(str),guides[rows].astype(str));np.testing.assert_array_equal(prov.original_total_counts,library[rows]);np.testing.assert_array_equal(prov.predicted_total_counts,library[rows])
    for ti,t in enumerate(targets):
     p0=(ci*300+ti)*400;lo,hi=int(ptr[p0]),int(ptr[p0+400]);count=sparse.csr_matrix((packed['X/data'][lo:hi],packed['X/indices'][lo:hi],ptr[p0:p0+401]-lo),shape=(400,len(genes)));assert (count.data>0).all() and np.isfinite(count.data).all() and np.equal(count.data,np.floor(count.data)).all();expected_rows=build.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([SEED,ci,ti])));actual_rows=prov.iloc[ti*400:(ti+1)*400].source_control_row.to_numpy(int);np.testing.assert_array_equal(actual_rows,expected_rows);assert len(set(actual_rows))==400;np.testing.assert_array_equal(count.sum(1).A1,library[actual_rows]);assert count.sum(1).A1.max()<=1000000;scanned+=400
     if ti in [0,50,150,299]:
      regenerated,_=emit(x,actual_rows,effects[ci,ti],int(own[ti]),SEED+ci*100000+ti,'residual10',(depth,cp,bulk),atlas);regenerated.sum_duplicates();regenerated.eliminate_zeros();regenerated.sort_indices()
      for key in ['data','indices','indptr']:np.testing.assert_array_equal(getattr(count,key),getattr(regenerated,key))
      replays.append(dict(context=c,target=str(t),full_cells=400,exact_count_replay=True));del regenerated
    print('Final V91 full-payload source/provenance/replay audit',c,flush=True);del x;gc.collect()
 assert scanned==360000 and len(replays)==12;write(A/'independent_validation.json',dict(passed=True,submission_sha256=ready['sha256'],prediction_sha256=sha(OUT/'prediction.h5ad'),all_packaged_counts_and_indices_compared=compared,cells_with_original_libraries_verified=scanned,exact_original_template_identities=True,full_groups_regenerated=replays,official_package_metadata=ready['package_meta'],accuracy_claim='Final biological accuracy is blind and is not established by file validation; see the frozen final plan.',auditor_sha256=sha(__file__)));print('Final V91 complete independent submission validation passed',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
