"""Independently verify V6 source separation, response invariants and the full submission."""
from pathlib import Path
import argparse
import json
import hashlib
import tarfile
import tempfile
import zstandard as zstd
import numpy as np
import pandas as pd
import anndata as ad
import h5py
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from evaluate_atlas_transfer import load_atlas_model
from v6_methods import moments,emit
from validate_gene300_models import strings
ROOT=build.ROOT;OUT=ROOT/'experiments/v6_fidelity_calibrated';A=OUT/'analysis'


def model_audit():
    record=json.loads((A/'reference_model.json').read_text());_,genes,targets=build.read_spec(ROOT)
    assert build.sha256(ROOT/'fit_v6_submission.py')==record['script_sha256']
    assert build.sha256(ROOT/'v6_methods.py')==record['methods_sha256']
    for name,key in [('reference_model_arrays.npz','model_arrays_sha256'),('source_generator_scores.csv','source_scores_sha256'),('source_profile_cache.npz','source_cache_sha256')]:
        assert build.sha256(A/name)==record[key]
    oldfile=ROOT/'experiments/v4_lineage_response/analysis/reference_model_arrays.npz';newfile=ROOT/'experiments/v5_atlas_calibrated/analysis/reference_model_arrays.npz'
    assert build.sha256(oldfile)==record['v4_arrays_sha256'] and build.sha256(newfile)==record['v5_arrays_sha256']
    with np.load(oldfile,allow_pickle=False) as z:old=z['effects'].astype(np.float64)
    with np.load(newfile,allow_pickle=False) as z:new=z['effects'].astype(np.float64)
    with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
        effect=z['effects'].copy();np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets)
    expected=old.copy();weight=record['setting']['blend'];assert 0<=weight<=.25
    for ci in [1,2]:
        same=(old[ci]*new[ci])>0
        expected[ci,same]=(1-weight)*old[ci,same]+weight*new[ci,same]
    np.testing.assert_array_equal(effect,expected.astype(np.float32))
    activity=pd.read_csv(A/'context_generator_ablation.csv')
    assert build.sha256(A/'context_generator_ablation.csv')==record['context_ablation_sha256']
    budget=activity.groupby(['context','candidate']).predicted_DE_calls.mean().unstack()
    for key,value in record['context_generator_guards'].items():
        assert bool((budget[key]/budget['v4']>=.9).all())==value
    assert record['context_generator_guards'].get(record['candidate'],True)
    assert record['neural_weight']==0 and not ((effect*old<0)&(old!=0)).any() and (effect[old==0]==0).all()
    identities=json.loads((A/'source_cell_identity.json').read_text())
    meta=pd.read_csv(ROOT/'reference/v5_h1/selected_cells.csv').set_index('source_cell_id');batches=sorted(meta.batch.unique());parts={b:i%2 for i,b in enumerate(batches)}
    assert len(identities)==100
    for group in identities:
        actual,template,reference=[group[k] for k in ['actual_source_ids','template_source_ids','reference_control_ids']]
        assert len(actual)==len(template)==group['cells']<=400
        assert len(set(actual))==len(actual) and len(set(template))==len(template) and len(set(reference))==len(reference)
        assert not set(template)&set(reference) and not set(actual)&(set(template)|set(reference))
        for names,label in [(actual,group['target_gene']),(template,'non-targeting'),(reference,'non-targeting')]:
            assert (meta.loc[names].target_gene==label).all()
            assert (meta.loc[names].batch.map(parts)==group['evaluation_partition']).all()
        assert group['effect_fit_partition']!=group['evaluation_partition']
    folds=record['target_outer_choices']
    tests=[]
    for fold in folds:
        assert not set(fold['train_targets'])&set(fold['test_targets']);tests.extend(fold['test_targets'])
    assert len(tests)==len(set(tests))==25
    source_scores=pd.read_csv(A/'source_generator_scores.csv');assert len(source_scores)==800
    assert not source_scores.duplicated(['target_gene','direction','seed','candidate']).any()
    neural_ids=json.loads((A/'neural_control_identity.json').read_text())
    for group in neural_ids:
        assert not set(group['template_ids'])&set(group['reference_ids'])
        assert not set(group['fit_guide_groups'])&(set(group['template_guide_groups'])|set(group['reference_guide_groups']))
        assert not set(group['template_guide_groups'])&set(group['reference_guide_groups'])
    report={'passed':True,'all_mean_effect_coordinates_reassembled':int(effect.size),'unchanged_A':True,'zero_new_coordinates_and_direction_flips':True,
        'source_evaluation_groups':100,'source_candidate_rows':800,'all_source_ID_and_batch_separation_verified':True,'target_outer_split_verified':True,
        'neural_template_weight':0,'scope':'Model/identity/provenance checks; source efficacy metrics are unanchored comparisons, not hidden response truth.',
        'validator_sha256':build.sha256(Path(__file__)),'model_sha256':build.sha256(A/'reference_model.json')}
    build.write_json(A/'model_validation.json',report)
    return record,effect,genes,targets


def payload_audit(genes,targets,generated_path):
    """Compare every packaged CSR value/index with the generated H5AD after float32 encoding."""
    with tempfile.TemporaryDirectory(prefix='v6_payload_check_') as directory:
        scratch=Path(directory).resolve()
        assert scratch.parent==Path(tempfile.gettempdir()).resolve() and scratch.name.startswith('v6_payload_check_')
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


def _package_audit(generated_path):
    record,effect,genes,targets=model_audit();ready=json.loads((A/'submission_ready.json').read_text());run=json.loads((A/'run_parameters.json').read_text())
    assert ready['official_prep_passed'] and ready['bytes']==(OUT/'prediction.vcc').stat().st_size
    assert ready['sha256']==build.sha256(OUT/'prediction.vcc')
    assert run['generator_sha256']==build.sha256(ROOT/'build_v6_submission.py')==ready['generator_sha256']
    assert run['methods_sha256']==record['methods_sha256'] and run['model_sha256']==build.sha256(A/'reference_model.json')
    from vcc.vccfile import validate_vcc,read_vcc_meta
    validate_vcc(OUT/'prediction.vcc');metadata=read_vcc_meta(OUT/'prediction.vcc');assert metadata==ready['package_meta']
    payload_values=payload_audit(genes,targets,generated_path)
    provenance=pd.read_csv(A/'prediction_control_rows.csv');stats=pd.read_csv(A/'prediction_group_statistics.csv')
    assert len(provenance)==360000 and len(stats)==900
    source=ad.read_h5ad(generated_path,backed='r')
    np.testing.assert_array_equal(source.var_names,genes)
    np.testing.assert_array_equal(source.obs.context.astype(str).to_numpy(),np.repeat(list('ABC'),120000))
    np.testing.assert_array_equal(source.obs.target_gene.astype(str).to_numpy(),np.tile(np.repeat(targets,400),3))
    source.file.close();atlas,_=load_atlas_model();replayed=0;checked_counts=0;row_count=0
    with h5py.File(generated_path,'r') as handle:
        ptr=handle['X/indptr'][:];assert ptr.shape==(360001,) and ptr[0]==0 and np.all(np.diff(ptr)>0)
        assert ptr[-1]==len(handle['X/data'])==len(handle['X/indices'])==metadata['nnz']
        for ci,context in enumerate('ABC'):
            original_path=ROOT/f'context_{context}.h5ad'
            assert build.sha256(original_path)==next(r['sha256'] for r in record['context_original_inputs'] if r['context']==context)
            original=ad.read_h5ad(original_path);x=original.X.astype(np.int32).tocsr();x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
            library=np.asarray(x.sum(axis=1,dtype=np.int64)).ravel();guides=original.obs.ntc_id.astype(str).to_numpy();ids=original.obs_names.to_numpy()
            prov=provenance[provenance.context==context];assert len(prov)==120000
            rows=prov.source_control_row.to_numpy(dtype=int);assert rows.min()>=0 and rows.max()<18400
            np.testing.assert_array_equal(prov.source_control_cell_id,ids[rows]);np.testing.assert_array_equal(prov.ntc_id,guides[rows]);np.testing.assert_array_equal(prov.original_total_counts,library[rows])
            for ti,target in enumerate(targets):
                group=prov.iloc[ti*400:(ti+1)*400]
                assert group.target_gene.eq(target).all() and group.source_control_row.nunique()==400
                np.testing.assert_array_equal(group.prediction_cell_offset,np.arange(400))
                guide_counts=group.ntc_id.value_counts();assert len(guide_counts)==46 and guide_counts.min()==8 and guide_counts.max()==9
            for start in range(ci*120000,(ci+1)*120000,1000):
                stop=start+1000;lo,hi=int(ptr[start]),int(ptr[stop]);data=handle['X/data'][lo:hi];indices=handle['X/indices'][lo:hi]
                assert np.issubdtype(data.dtype,np.integer) and (data>0).all() and indices.min()>=0 and indices.max()<18533
                libs=np.add.reduceat(data.astype(np.int64),ptr[start:stop]-lo)
                assert libs.min()>0 and libs.max()<=1000000
                np.testing.assert_array_equal(libs,provenance.predicted_total_counts.iloc[start:stop]);row_count+=len(libs)
            depth,p,b,_=moments(x);prepared=(depth,p,b,np.asarray(x.mean(axis=0)).ravel())
            for target in ['ABCD1','ACLY','VIM','SMARCA5']:
                ti=int(np.flatnonzero(targets==target)[0]);offset=(ci*300+ti)*400;rows=prov.iloc[ti*400:(ti+1)*400].source_control_row.to_numpy(dtype=int)
                expected_rows=build.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([20261007,ci,ti])))
                np.testing.assert_array_equal(rows,expected_rows)
                own=int(np.flatnonzero(genes==target)[0]);block,diagnostic=emit(x,rows,effect[ci,ti],own,20261007+ci*100000+ti,record['setting'],atlas,prepared=prepared)
                lo,hi=int(ptr[offset]),int(ptr[offset+400])
                actual=sparse.csr_matrix((handle['X/data'][lo:hi],handle['X/indices'][lo:hi],ptr[offset:offset+401]-lo),shape=(400,18533))
                assert (block-actual).nnz==0
                checked_counts+=400*18533;replayed+=1
            del x,original;print(context,'full count scan, identities and deterministic replay passed',flush=True)
    assert row_count==360000 and replayed==12
    report={'passed':True,'all_packaged_count_values_and_CSR_indices_match_H5AD':True,'packaged_nonzero_values_checked':payload_values,'all_stored_counts_and_cell_libraries_checked':True,'verified_rows':row_count,'verified_groups':900,
        'all_original_control_rows_ids_guides_and_libraries_verified':True,'exact_replayed_groups':replayed,'exact_replayed_gene_counts':checked_counts,
        'official_package_header_and_sha256_verified':True,'validator_sha256':build.sha256(Path(__file__)),
        'model_sha256':build.sha256(A/'reference_model.json'),'prediction_sha256':build.sha256(generated_path),'package_sha256':ready['sha256'],
        'scope':'Full integrity/provenance verification and exact selected count replay, not hidden leaderboard accuracy.'}
    assert report['prediction_sha256']==build.sha256(OUT/'prediction.h5ad')
    build.write_json(A/'independent_validation.json',report);print('V6 submission independent audit passed',flush=True)


def package_audit():
    # A verified local copy avoids thousands of small compressed-HDF5 reads over SMB.
    original=OUT/'prediction.h5ad';before=original.stat()
    with tempfile.TemporaryDirectory(prefix='v6_matrix_check_') as directory:
        scratch=Path(directory).resolve()
        assert scratch.parent==Path(tempfile.gettempdir()).resolve() and scratch.name.startswith('v6_matrix_check_')
        local=scratch/'prediction.h5ad';digest=hashlib.sha256()
        with original.open('rb') as source,local.open('wb') as destination:
            while block:=source.read(8*1024*1024):
                digest.update(block);destination.write(block)
        after=original.stat()
        assert before.st_size==after.st_size==local.stat().st_size and before.st_mtime_ns==after.st_mtime_ns
        assert build.sha256(local)==digest.hexdigest()
        print('Generated H5AD cached locally with identical SHA-256',flush=True)
        _package_audit(local)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--model-only',action='store_true');args=parser.parse_args()
    with threadpool_limits(limits=2):
        if args.model_only:model_audit();print('V6 model/identity audit passed',flush=True)
        else:package_audit()
