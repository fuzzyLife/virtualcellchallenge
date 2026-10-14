"""Verify every V24 packaged count/axis, original control identity and sampled replay."""
from pathlib import Path
import json,gc,tempfile,subprocess,hashlib
import numpy as np,pandas as pd,anndata as ad,h5py
pd.options.future.infer_string=False
from scipy.stats import rankdata
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from build_v24_submission import ROOT,OUT,A,load,SEED
from v16_methods import emit
from evaluate_atlas_transfer import load_atlas_model
from v6_methods import moments

def main():
    from v24_runtime import require_science
    require_science()
    manifest,genes,targets,model,effects=load();ready=json.loads((A/'submission_ready.json').read_text());run=json.loads((A/'run_parameters.json').read_text())
    assert ready['sha256']==build.sha256(OUT/'prediction.vcc') and ready['bytes']==(OUT/'prediction.vcc').stat().st_size
    assert ready['model_sha256']==build.sha256(A/'reference_model.json')==run['model_sha256']
    assert run['generator_sha256']==build.sha256(ROOT/'build_v24_submission.py') and run['emitter_sha256']==build.sha256(ROOT/'v16_methods.py')
    with np.load(ROOT/'experiments/v16_control_preserving/analysis/reference_model_arrays.npz',allow_pickle=False) as z:parent_effects=z['effects'].copy()
    # Independent reconstruction from the frozen teacher, with no predict/gated calls.
    from fit_v7_response import read_cache
    data=read_cache()
    with np.load(A/'models/teacher_all.npz',allow_pickle=False) as z:teacher={k:z[k].copy() for k in z.files}
    assert build.sha256(A/'models/teacher_all.npz')==model['teacher_all_sha256']
    def rho(a,b):
        ok=np.isfinite(a)&np.isfinite(b);a=rankdata(a[ok],method='average');b=rankdata(b[ok],method='average');a-=a.mean();b-=b.mean();den=np.sqrt(np.dot(a,a)*np.dot(b,b));return float(np.clip(np.dot(a,b)/den,-1,1)) if den else np.nan
    for ci,context in enumerate(manifest['contexts']):
        with np.load(A/'production_inputs'/(context+'_query.npz'),allow_pickle=False) as z:query=z['query_mean'].copy();np.testing.assert_array_equal(z['genes'],genes)
        axis=~np.isin(genes,targets);broad=np.array([rho(query[axis],p[axis]) for p in teacher['profiles']]);other=np.array([rho(query[axis],data['profiles'][v,axis]) for v in [0,4,5,6]]);compatibility=float(np.clip((np.nanmax(broad)-np.nanmax(other))/.1,0,1));state=teacher['state_axis'];r=np.array([rho(query[state],p[state]) for p in teacher['profiles']]);r=np.nan_to_num(r,nan=-1);weights=np.exp((r-r.max())/.1);weights/=weights.sum();confidence=np.einsum('c,ctg->tg',weights,teacher['confidence_de']);numerator=np.einsum('c,ctg,ctg->tg',weights,teacher['confidence_de'],teacher['value']);mu=np.divide(numerator,confidence,out=np.zeros_like(numerator),where=confidence>0);gate=np.clip(query/5,0,1);own=np.array([np.flatnonzero(genes==t)[0] for t in targets]);mu=np.asarray([np.clip(v,-3,3)*gate*gate[o] for v,o in zip(mu,own)],np.float32);mu[np.arange(300),own]=0;delta=compatibility*confidence*np.clip(mu-parent_effects[ci],-.5,.5);delta[np.arange(300),own]=0;expected=np.clip(parent_effects[ci]+.5*delta,-3,3).astype(np.float32);expected[np.arange(300),own]=0;np.testing.assert_array_equal(effects[ci],expected)
        assert compatibility==model['routing'][context]['lineage_compatibility_gate'] and abs(effects[ci]-parent_effects[ci]).max()<=.250001
    change_table=pd.read_csv(A/'production_response_changes.csv')
    own=np.array([np.flatnonzero(genes==t)[0] for t in targets]);assert (effects[:,np.arange(300),own]==0).all() and (abs(effects)<=3.000001).all()
    audit=json.loads((A/'independent_experiment_audit.json').read_text());assert audit['passed']
    provenance=pd.read_csv(A/'prediction_control_rows.csv');groups=pd.read_csv(A/'prediction_group_statistics.csv');assert len(provenance)==360000 and len(groups)==900
    with tempfile.TemporaryDirectory(prefix='v24_payload_check_') as directory:
        scratch=Path(directory).resolve();assert scratch.parent==Path(tempfile.gettempdir()).resolve() and scratch.name.startswith('v24_payload_check_')
        payload=scratch/'prediction.h5ad'
        result=subprocess.run([str(ROOT/'.competition_env/Scripts/python.exe'),str(ROOT/'unpack_v24_official_payload.py'),str(OUT/'prediction.vcc'),str(payload)],capture_output=True,text=True,encoding='utf-8',errors='replace');assert result.returncode==0,result.stderr[-2000:];metadata=json.loads(result.stdout);assert metadata==ready['package_meta']
        print('V24 payload decompressed for complete independent comparison',flush=True)
        packed=ad.read_h5ad(payload,backed='r');np.testing.assert_array_equal(packed.var_names,genes)
        np.testing.assert_array_equal(packed.obs.context.astype(str),np.repeat(manifest['contexts'],120000));np.testing.assert_array_equal(packed.obs.target_gene.astype(str),np.tile(np.repeat(targets,400),3));assert packed.obs_names.is_unique;packed.file.close()
        checked=0;scanned=0;replayed=0;atlas,_=load_atlas_model()
        with h5py.File(payload,'r') as packed,h5py.File(OUT/'prediction.h5ad','r') as original:
            ptr=packed['X/indptr'][:];np.testing.assert_array_equal(ptr,original['X/indptr'][:]);assert ptr.shape==(360001,) and ptr[0]==0 and (np.diff(ptr)>0).all()
            assert packed['X'].attrs['shape'].tolist()==[360000,18533] and ptr[-1]==metadata['nnz']
            for start in range(0,len(original['X/data']),4000000):
                stop=min(start+4000000,len(original['X/data']))
                for name in ['data','indices']:np.testing.assert_array_equal(packed['X/'+name][start:stop],original['X/'+name][start:stop])
                checked+=stop-start
                if checked%400000000==0:print('V24 payload elements compared',checked,flush=True)
            for ci,context in enumerate(manifest['contexts']):
                obj=ad.read_h5ad(ROOT/f'context_{context}.h5ad');np.testing.assert_array_equal(obj.var_names,genes);x=obj.X.tocsr().astype(np.int32);x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
                library=np.asarray(x.sum(1,dtype=np.int64)).ravel();guides=obj.obs.ntc_id.astype(str).to_numpy();prov=provenance[provenance.context==context];rows=prov.source_control_row.to_numpy(int)
                np.testing.assert_array_equal(prov.source_control_cell_id,obj.obs_names.to_numpy()[rows]);np.testing.assert_array_equal(prov.ntc_id,guides[rows]);np.testing.assert_array_equal(prov.original_total_counts,library[rows])
                for ti,target in enumerate(targets):
                    group=prov.iloc[ti*400:(ti+1)*400];assert group.target_gene.eq(target).all() and group.source_control_row.nunique()==400
                    np.testing.assert_array_equal(group.prediction_cell_offset,np.arange(400));chosen=build.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([SEED,ci,ti])));np.testing.assert_array_equal(group.source_control_row,chosen);g=group.ntc_id.value_counts();assert len(g)==46 and g.min()==8 and g.max()==9
                for start in range(ci*120000,(ci+1)*120000,1000):
                    stop=start+1000;lo,hi=int(ptr[start]),int(ptr[stop]);values=packed['X/data'][lo:hi];indices=packed['X/indices'][lo:hi]
                    assert np.isfinite(values).all() and (values>0).all() and (values==np.floor(values)).all() and indices.min()>=0 and indices.max()<18533
                    libs=np.add.reduceat(values.astype(np.int64),ptr[start:stop]-lo);np.testing.assert_array_equal(libs,provenance.predicted_total_counts.iloc[start:stop]);np.testing.assert_array_equal(libs,provenance.original_total_counts.iloc[start:stop]);assert libs.min()>0 and libs.max()<=1000000;scanned+=len(libs)
                depth,p,b,_=moments(x)
                with np.load(A/'production_inputs'/(context+'_query.npz'),allow_pickle=False) as z:np.testing.assert_array_equal(z['query_mean'],p*1e6)
                strongest=change_table[change_table.context==context].sort_values('sum_absolute_LFC_correction',ascending=False).target_gene.astype(str).tolist();selected=['ABCD1','CHKB']+[t for t in strongest if t not in ['ABCD1','CHKB']][:2]
                for target in selected:
                    ti=int(np.flatnonzero(targets==target)[0]);rows=prov.iloc[ti*400:(ti+1)*400].source_control_row.to_numpy(int)
                    chosen=build.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([SEED,ci,ti])));np.testing.assert_array_equal(rows,chosen)
                    expected,_=emit(x,rows,effects[ci,ti],int(own[ti]),SEED+ci*100000+ti,'residual10',(depth,p,b),atlas)
                    start=(ci*300+ti)*400;lo,hi=int(ptr[start]),int(ptr[start+400]);actual=sparse.csr_matrix((packed['X/data'][lo:hi],packed['X/indices'][lo:hi],ptr[start:start+401]-lo),shape=(400,18533));assert (expected-actual).nnz==0;replayed+=1
                del x,obj;gc.collect();print('V24 independent counts/identities/replay passed',context,flush=True)
    assert scanned==360000 and replayed==12 and checked==metadata['nnz']
    build.write_json(A/'independent_validation.json',{'passed':True,'verified_rows':scanned,'verified_groups':900,'packaged_nonzero_values_checked':checked,'all_packaged_counts_indices_pointers_axes_match':True,'all_production_responses_independently_reconstructed':True,
        'all_original_control_rows_ids_guides_libraries_verified':True,'all_counts_nonnegative_integer_values_and_within_library_limit':True,'exact_replayed_groups':replayed,'exact_replayed_gene_counts':replayed*400*18533,
        'package_sha256':ready['sha256'],'prediction_sha256':build.sha256(OUT/'prediction.h5ad'),'provenance_sha256':build.sha256(A/'prediction_control_rows.csv'),'model_sha256':build.sha256(A/'reference_model.json'),'validator_sha256':build.sha256(Path(__file__)),
        'scope':'Full package integrity, unique original-control provenance and selected deterministic regeneration; no hidden leaderboard accuracy guarantee.'})
    print('V24 independent complete-payload audit passed',flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
