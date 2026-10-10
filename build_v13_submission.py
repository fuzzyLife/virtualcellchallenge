"""Generate V13 from the frozen selection and package with the official CLI."""
from pathlib import Path
import argparse
import gc
import json
import subprocess
import sys
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import anndata as ad
import h5py
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from evaluate_atlas_transfer import load_atlas_model
from v6_methods import moments
from v11_methods import emit
ROOT=build.ROOT;OUT=ROOT/'experiments/v13_response_programs';A=OUT/'analysis';SEED=20261011


def load():
    manifest,genes,targets=build.read_spec(ROOT)
    model=json.loads((A/'reference_model.json').read_text(encoding='utf-8'))
    for file,digest in model['dependencies'].items():
        if build.sha256(ROOT/file)!=digest:raise ValueError('Frozen dependency changed: '+file)
    for file,key in [('reference_model_arrays.npz','model_arrays_sha256'),(model['validation_plan_file'],'validation_plan_sha256'),
                     ('count_selection.json','count_selection_sha256'),('confirmation_result.json','confirmation_result_sha256')]:
        if build.sha256(A/file)!=model[key]:raise ValueError('Frozen model record changed: '+file)
    for record in model['context_original_inputs']:
        if build.sha256(ROOT/record['file'])!=record['sha256']:raise ValueError('Context inputs changed')
    with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
        np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets);effect=z['effects'].copy()
    assert effect.shape==(3,300,18533) and np.isfinite(effect).all()
    return manifest,genes,targets,model,effect


def generate(preflight=False,workers=4):
    manifest,genes,targets,model,effect=load();A.mkdir(parents=True,exist_ok=True)
    path=OUT/'prediction.h5ad'
    if not preflight and (path.exists() or (OUT/'prediction.vcc').exists()):raise FileExistsError('V13 output exists; preserve it and use a new version for changes')
    atlas,atlas_sha=load_atlas_model();records=[];provenance=[]
    # Reuse only exact unchanged V11 effects with identical original controls, seed and emitter.
    old_handle=None;old_effect=None;old_statistics=None;old_sha=None
    old_out=ROOT/'experiments/v11_consensus_coverage';old_file=old_out/'prediction.h5ad'
    if not preflight and old_file.exists():
        from build_v11_submission import load as load_v11,SEED as OLD_SEED
        _,_,_,old_model,old_effect=load_v11();assert OLD_SEED==SEED and old_model['setting']==model['setting']
        assert old_model['context_original_inputs']==model['context_original_inputs']
        prior=json.loads((old_out/'analysis/independent_validation.json').read_text(encoding='utf-8'));assert prior['passed']
        old_sha=build.sha256(old_file);assert old_sha==prior['prediction_sha256']
        old_handle=h5py.File(old_file,'r');old_ptr=old_handle['X/indptr'][:]
        old_statistics=pd.read_csv(old_out/'analysis/prediction_group_statistics.csv').set_index(['context','target_gene'])
    settings=model['setting'];n=manifest['cells_per_pert'];contexts=manifest['contexts']
    names=['ACLY','CHKB','DHX35','ZFP62','VIM','SMARCA5','S100A11','FXYD5']
    writer=None
    if not preflight:
        obs=pd.DataFrame({'target_gene':pd.Categorical(np.tile(np.repeat(targets,n),len(contexts)),categories=targets),
            'context':pd.Categorical(np.repeat(contexts,len(targets)*n),categories=contexts)},index=[f'v13_{i:06d}' for i in range(n*len(targets)*len(contexts))])
        writer=build.SparseH5Writer(path.with_suffix('.h5ad.partial'),obs,genes)
    for ci,context in enumerate(contexts):
        x,guides,library,_,_,statistics=build.load_control(ROOT,context,genes,manifest)
        ids=ad.read_h5ad(ROOT/f'context_{context}.h5ad',backed='r');cell_ids=ids.obs_names.to_numpy();ids.file.close()
        if not preflight:statistics.to_csv(A/f'gene_statistics_{context}.csv',index=False)
        baseline=moments(x);prepared=(baseline[0],baseline[1],baseline[2],np.asarray(x.mean(axis=0)).ravel())
        jobs=[(ti,t) for ti,t in enumerate(targets) if not preflight or t in names]
        def compute(job):
            ti,target=job;seed=SEED+ci*100000+ti
            rng=np.random.default_rng(np.random.SeedSequence([SEED,ci,ti]));rows=build.guide_stratified_sample(guides,library,n,rng)
            own=int(np.flatnonzero(genes==target)[0])
            reused=old_handle is not None and np.array_equal(effect[ci,ti],old_effect[ci,ti])
            if reused:
                offset=(ci*300+ti)*n;lo,hi=int(old_ptr[offset]),int(old_ptr[offset+n])
                block=sparse.csr_matrix((old_handle['X/data'][lo:hi],old_handle['X/indices'][lo:hi],old_ptr[offset:offset+n+1]-lo),shape=(n,len(genes)))
                np.testing.assert_array_equal(np.asarray(block.sum(axis=1)).ravel(),library[rows])
                previous=old_statistics.loc[(context,target)]
                diagnostic={k:float(previous[k]) for k in ['noise_amplitude','bulk_response_amplitude','mean_probability_L1','bulk_probability_L1']}
            else:block,diagnostic=emit(x,rows,effect[ci,ti],own,seed,settings,atlas=atlas,prepared=prepared)
            depth,p,b,var=moments(block)
            row={'context':context,'target_gene':target,'target_index':ti,'reused_unchanged_v11':bool(reused),'nnz':int(block.nnz),'mean_library':float(depth.mean()),
                'mean_detected_genes':float(block.getnnz(axis=1).mean()),'total_CPM_variance':float(var.sum()*1e12),
                'own_generated_mean_raw':float(block[:,own].mean()),**diagnostic}
            return block,rows,row,depth
        with threadpool_limits(limits=2):
            with ThreadPoolExecutor(max_workers=workers) as pool:
                pending={i:pool.submit(compute,jobs[i]) for i in range(min(workers,len(jobs)))}
                for index,job in enumerate(jobs):
                    block,rows,row,depth=pending.pop(index).result();records.append(row)
                    if writer is not None:
                        writer.append(block)
                        provenance.append(pd.DataFrame({'context':context,'target_gene':job[1],'prediction_cell_offset':np.arange(n),
                            'source_control_row':rows,'source_control_cell_id':cell_ids[rows],'ntc_id':guides[rows],
                            'original_total_counts':library[rows],'predicted_total_counts':depth.astype(np.int64)}))
                    next_index=index+workers
                    if next_index<len(jobs):pending[next_index]=pool.submit(compute,jobs[next_index])
                    if preflight or (job[0]+1)%10==0:build.log(f'V13 {context}: {job[0]+1}/300; candidate {model["candidate"]}')
                    del block
        del x;gc.collect()
    pd.DataFrame(records).to_csv(A/('preflight_group_statistics.csv' if preflight else 'prediction_group_statistics.csv'),index=False)
    if preflight:
        build.write_json(A/'preflight.json',{'passed':True,'groups':len(records),'model_sha256':build.sha256(A/'reference_model.json'),
            'generator_sha256':build.sha256(Path(__file__)),'methods_sha256':build.sha256(ROOT/(model['methods_module']+'.py'))})
        return
    writer.close();path.with_suffix('.h5ad.partial').replace(path)
    if old_handle is not None:old_handle.close()
    pd.concat(provenance,ignore_index=True).to_csv(A/'prediction_control_rows.csv',index=False)
    build.write_json(A/'run_parameters.json',{'created_utc':datetime.now(timezone.utc).isoformat(),'seed':SEED,'setting':settings,'response_candidate':model['candidate'],
        'knockdown_assumption':.9,'generator_sha256':build.sha256(Path(__file__)),'methods_sha256':build.sha256(ROOT/(model['methods_module']+'.py')),
        'model_sha256':build.sha256(A/'reference_model.json'),'upstream_sha256':atlas_sha,'workers':workers,'unchanged_v11_groups_reused':int(sum(r['reused_unchanged_v11'] for r in records)),'verified_v11_h5ad_sha256':old_sha,
        'simulation':'Frozen V13 selected response model, V11 dispersion setting; guide/depth-stratified unique control samples; synthetic dual-moment integer populations; own target handled separately.'})


def package():
    manifest,genes,targets=build.read_spec(ROOT)
    verified=build.validate_h5ad(OUT/'prediction.h5ad',genes,targets,manifest['contexts'],manifest['cells_per_pert'])
    build.write_json(A/'prediction_validation.json',verified)
    from package_v13_memory_safe import main
    main()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--preflight',action='store_true');parser.add_argument('--package-only',action='store_true');parser.add_argument('--workers',type=int,default=4);parser.add_argument('--wait-for-model',action='store_true');args=parser.parse_args()
    if args.wait_for_model:
        import time
        while not (A/'reference_model.json').exists():time.sleep(20)
    if not 1<=args.workers<=8:parser.error('workers must be 1..8')
    if not args.package_only:generate(args.preflight,args.workers)
    if not args.preflight:package()
