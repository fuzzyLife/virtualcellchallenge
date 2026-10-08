"""Independent V8 split, count, provenance, payload, and deterministic replay audit."""
from pathlib import Path
import argparse,json,hashlib,tarfile,tempfile,gc
import numpy as np,pandas as pd,anndata as ad,h5py,zstandard as zstd
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from build_v8_submission import load,SEED
from v8_methods import emit
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model
from fit_v7_response import read_cache,raw_inputs
ROOT=build.ROOT;OUT=ROOT/'experiments/v8_corrected_transfer';A=OUT/'analysis'


def model_audit():
    manifest,genes,targets,record,effects=load()
    data=read_cache();assert build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')==record['source_cache_sha256']
    own=np.array([np.flatnonzero(genes==t)[0] for t in targets])
    assert (effects[:,np.arange(300),own]==0).all() and np.isfinite(effects).all() and (abs(effects)<=3.000001).all()
    plan=json.loads((A/'validation_plan.json').read_text());profile=json.loads((A/'profile_selection.json').read_text())
    count=json.loads((A/'safe_count_selection.json').read_text());confirm=json.loads((A/'safe_confirmation_result.json').read_text())
    assert count['confirmation_read']==False and profile['confirmation_read']==False
    assert record['response']==count['response'] and record['generator']==count['generator']
    assert confirm['candidate']==count['candidate'] and confirm['frozen_calibration_selection_sha256']==build.sha256(A/'safe_count_selection.json')
    assignments={'calibration':set(),'confirmation':set()};groups=0
    for family,training in profile['training'].items():
        rows=np.array(training['train_rows'],dtype=int)
        assert (data['row_families'][rows]!=family).all()
    for family in ['H1','HCT116','HEK293T']:
        x,meta,measured,source=raw_inputs(family,genes)
        for context in [c for c in plan['contexts'] if c['family']==family]:
            left,right=set(context['left']),set(context['right']);assert not left&right
            assert meta.iloc[list(left|right)].gene_target.eq('Non-Targeting').all()
            for seed in plan['seeds']:
                cells=[c for c in context['cells'] if c['seed']==seed]
                assert all(c['template']==cells[0]['template'] and c['reference']==cells[0]['reference'] for c in cells)
                for cell in cells:
                    assert len(cell['template'])==len(set(cell['template']))==400
                    assert len(cell['reference'])==len(set(cell['reference']))==400
                    assert set(cell['template'])<=left and set(cell['reference'])<=right
                    assert not set(cell['actual'])&(left|right)
                    assert meta.iloc[cell['actual']].gene_target.eq(cell['target']).all()
                    assert meta.iloc[cell['actual']]['sample'].isin(context['batches']).all()
                    assert len(cell['actual'])==len(set(cell['actual'])) and 60<=len(cell['actual'])<=200
                    assert not set(cell['replicate_left'])&set(cell['replicate_right'])
                    assignments[cell['split']].add(cell['target']);groups+=1
        del x;gc.collect()
    assert not assignments['calibration']&assignments['confirmation']
    report={'passed':True,'evaluation_target_seed_groups':groups,'globally_disjoint_calibration_confirmation_targets':True,
        'whole_lineage_response_training_exclusion':True,'one_shared_reference_and_target_independent_templates_per_stratum_seed':True,
        'distinct_actual_cells_and_disjoint_control_populations':True,'own_gene_excluded_from_effects':True,
        'model_sha256':build.sha256(A/'reference_model.json'),'validation_plan_sha256':build.sha256(A/'validation_plan.json')}
    build.write_json(A/'model_validation.json',report)
    return manifest,genes,targets,record,effects


def payload_audit(genes,targets,generated_path):
    """Compare every packaged CSR value/index with the generated H5AD after float32 encoding."""
    with tempfile.TemporaryDirectory(prefix='v8_payload_check_') as directory:
        scratch=Path(directory).resolve()
        assert scratch.parent==Path(tempfile.gettempdir()).resolve() and scratch.name.startswith('v8_payload_check_')
        payload=scratch/'pred.h5ad'
        with tarfile.open(OUT/'prediction.vcc','r:') as archive:
            member=archive.extractfile('pred.h5ad.zst')
            assert member is not None
            with member, payload.open('wb') as destination:
                zstd.ZstdDecompressor().copy_stream(member,destination)
        print('Compressed payload decompressed; comparing all encoded counts and axes',flush=True)
        packed=ad.read_h5ad(payload,backed='r')
        np.testing.assert_array_equal(packed.var_names,genes)
        np.testing.assert_array_equal(packed.obs.context.astype(str).to_numpy(),np.repeat(list('ABC'),120000))
        np.testing.assert_array_equal(packed.obs.target_gene.astype(str).to_numpy(),np.tile(np.repeat(targets,400),3))
        packed.file.close()
        values=0
        with h5py.File(payload,'r') as compressed, h5py.File(generated_path,'r') as generated:
            np.testing.assert_array_equal(compressed['X/indptr'][:],generated['X/indptr'][:])
            assert compressed['X'].attrs['shape'].tolist()==[360000,18533]
            assert len(compressed['X/data'])==len(generated['X/data'])
            for start in range(0,len(generated['X/data']),4000000):
                stop=min(start+4000000,len(generated['X/data']))
                np.testing.assert_array_equal(compressed['X/data'][start:stop],generated['X/data'][start:stop])
                np.testing.assert_array_equal(compressed['X/indices'][start:stop],generated['X/indices'][start:stop])
                values+=stop-start
                if values%400000000==0:print(f'Packaged values compared: {values:,} / {len(generated["X/data"]):,}',flush=True)
        print('Compressed payload: every count, CSR index, pointer and label matches generated H5AD',flush=True)
        return values


def package_audit():
    manifest,genes,targets,record,effects=model_audit()
    ready=json.loads((A/'submission_ready.json').read_text());run=json.loads((A/'run_parameters.json').read_text())
    assert ready['sha256']==build.sha256(OUT/'prediction.vcc') and ready['bytes']==(OUT/'prediction.vcc').stat().st_size
    assert run['generator_sha256']==build.sha256(ROOT/'build_v8_submission.py')==ready['generator_sha256']
    assert run['methods_sha256']==build.sha256(ROOT/'v8_methods.py')==ready['methods_sha256']
    from vcc.vccfile import validate_vcc,read_vcc_meta
    validate_vcc(OUT/'prediction.vcc');metadata=read_vcc_meta(OUT/'prediction.vcc');assert metadata==ready['package_meta']
    with tempfile.TemporaryDirectory(prefix='v8_matrix_check_') as directory:
        scratch=Path(directory).resolve();assert scratch.parent==Path(tempfile.gettempdir()).resolve()
        local=scratch/'prediction.h5ad';original=OUT/'prediction.h5ad';before=original.stat();digest=hashlib.sha256()
        with original.open('rb') as src,local.open('wb') as dst:
            while block:=src.read(8*1024*1024):digest.update(block);dst.write(block)
        after=original.stat();assert before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns
        assert build.sha256(local)==digest.hexdigest()
        packaged=payload_audit(genes,targets,local)
        provenance=pd.read_csv(A/'prediction_control_rows.csv');groups=pd.read_csv(A/'prediction_group_statistics.csv')
        assert len(provenance)==360000 and len(groups)==900
        atlas,_=load_atlas_model();replayed=0;scanned=0
        with h5py.File(local,'r') as handle:
            ptr=handle['X/indptr'][:];assert ptr.shape==(360001,) and ptr[0]==0 and (np.diff(ptr)>0).all()
            assert ptr[-1]==len(handle['X/data'])==len(handle['X/indices'])==metadata['nnz']
            for ci,context in enumerate('ABC'):
                original_context=ad.read_h5ad(ROOT/f'context_{context}.h5ad');x=original_context.X.tocsr().astype(np.int32)
                x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
                library=np.asarray(x.sum(axis=1,dtype=np.int64)).ravel();guides=original_context.obs.ntc_id.astype(str).to_numpy()
                prov=provenance[provenance.context==context];rows=prov.source_control_row.to_numpy(int)
                np.testing.assert_array_equal(prov.source_control_cell_id,original_context.obs_names.to_numpy()[rows])
                np.testing.assert_array_equal(prov.ntc_id,guides[rows]);np.testing.assert_array_equal(prov.original_total_counts,library[rows])
                for ti,target in enumerate(targets):
                    group=prov.iloc[ti*400:(ti+1)*400]
                    assert group.target_gene.eq(target).all() and group.source_control_row.nunique()==400
                    np.testing.assert_array_equal(group.prediction_cell_offset,np.arange(400))
                    guide_counts=group.ntc_id.value_counts();assert len(guide_counts)==46 and guide_counts.min()==8 and guide_counts.max()==9
                for start in range(ci*120000,(ci+1)*120000,1000):
                    stop=start+1000;lo,hi=int(ptr[start]),int(ptr[stop]);values=handle['X/data'][lo:hi];indices=handle['X/indices'][lo:hi]
                    assert np.issubdtype(values.dtype,np.integer) and (values>0).all() and indices.min()>=0 and indices.max()<18533
                    libs=np.add.reduceat(values.astype(np.int64),ptr[start:stop]-lo)
                    np.testing.assert_array_equal(libs,provenance.predicted_total_counts.iloc[start:stop])
                    assert libs.min()>0 and libs.max()<=1000000;scanned+=len(libs)
                depths,p,b,_=moments(x);prepared=(depths,p,b,np.asarray(x.mean(axis=0)).ravel())
                for target in ['ABCD1','ACLY','VIM','SMARCA5']:
                    ti=int(np.flatnonzero(targets==target)[0]);rows=prov.iloc[ti*400:(ti+1)*400].source_control_row.to_numpy(int)
                    expected_rows=build.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([SEED,ci,ti])))
                    np.testing.assert_array_equal(rows,expected_rows)
                    own=int(np.flatnonzero(genes==target)[0]);expected,diag=emit(x,rows,effects[ci,ti],own,SEED+ci*100000+ti,record['generator'],prepared=prepared,atlas=atlas)
                    offset=(ci*300+ti)*400;lo,hi=int(ptr[offset]),int(ptr[offset+400])
                    actual=sparse.csr_matrix((handle['X/data'][lo:hi],handle['X/indices'][lo:hi],ptr[offset:offset+401]-lo),shape=(400,18533))
                    assert (expected-actual).nnz==0;replayed+=1
                del x,original_context;gc.collect();print('V8',context,'full counts, provenance and replay passed',flush=True)
        assert scanned==360000 and replayed==12
    build.write_json(A/'independent_validation.json',{'passed':True,'verified_rows':scanned,'verified_groups':900,
        'all_packaged_values_indices_pointers_and_axes_match':True,'packaged_nonzero_values_checked':packaged,
        'all_original_control_rows_ids_guides_and_libraries_verified':True,'all_counts_integer_nonnegative_and_within_depth_limit':True,
        'exact_replayed_groups':replayed,'exact_replayed_gene_counts':replayed*400*18533,
        'model_sha256':build.sha256(A/'reference_model.json'),'package_sha256':ready['sha256'],'prediction_sha256':digest.hexdigest(),
        'validator_sha256':build.sha256(Path(__file__)),'scope':'Integrity and reproducibility; no hidden leaderboard performance claim.'})
    print('V8 independent full-payload audit passed',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model-only',action='store_true');args=p.parse_args()
    with threadpool_limits(limits=2):
        if args.model_only:model_audit();print('V8 model and sampling audit passed',flush=True)
        else:package_audit()
