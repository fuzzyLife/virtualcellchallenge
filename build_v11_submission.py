"""Generate V11 from the frozen selection and package with the official CLI."""
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
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from evaluate_atlas_transfer import load_atlas_model
from v6_methods import moments
from v11_methods import emit
ROOT=build.ROOT;OUT=ROOT/'experiments/v11_consensus_coverage';A=OUT/'analysis';SEED=20261011


def load():
    manifest,genes,targets=build.read_spec(ROOT)
    model=json.loads((A/'reference_model.json').read_text())
    for file,digest in model['dependencies'].items():
        if build.sha256(ROOT/file)!=digest:raise ValueError('Frozen dependency changed: '+file)
    for file,key in [('reference_model_arrays.npz','model_arrays_sha256'),('validation_plan.json','validation_plan_sha256'),
                     ('count_selection.json','count_selection_sha256'),
                     ('confirmation_result.json','confirmation_result_sha256'),('decoder_neutral_manifest.json','decoder_neutral_manifest_sha256')]:
        if build.sha256(A/file)!=model[key]:raise ValueError('Frozen model record changed: '+file)
    assert build.sha256(ROOT/'experiments/v4_lineage_response/analysis/source_components.npz')==model['v4_source_components_sha256']
    for record in model['context_original_inputs']:
        if build.sha256(ROOT/record['file'])!=record['sha256']:raise ValueError('Context inputs changed')
    with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
        np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets);effect=z['effects'].copy()
    assert effect.shape==(3,300,18533) and np.isfinite(effect).all()
    return manifest,genes,targets,model,effect


def generate(preflight=False,workers=4):
    manifest,genes,targets,model,effect=load();A.mkdir(parents=True,exist_ok=True)
    path=OUT/'prediction.h5ad'
    if not preflight and (path.exists() or (OUT/'prediction.vcc').exists()):raise FileExistsError('V11 output exists; preserve it and use a new version for changes')
    atlas,atlas_sha=load_atlas_model();records=[];provenance=[]
    settings=model['setting'];n=manifest['cells_per_pert'];contexts=manifest['contexts']
    names=['ACLY','ADNP','VIM','SMARCA5','S100A11','FXYD5','ALG12','MED25']
    writer=None
    if not preflight:
        obs=pd.DataFrame({'target_gene':pd.Categorical(np.tile(np.repeat(targets,n),len(contexts)),categories=targets),
            'context':pd.Categorical(np.repeat(contexts,len(targets)*n),categories=contexts)},index=[f'v11_{i:06d}' for i in range(n*len(targets)*len(contexts))])
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
            block,diagnostic=emit(x,rows,effect[ci,ti],own,seed,settings,atlas=atlas,prepared=prepared)
            depth,p,b,var=moments(block)
            row={'context':context,'target_gene':target,'target_index':ti,'nnz':int(block.nnz),'mean_library':float(depth.mean()),
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
                    if preflight or (job[0]+1)%10==0:build.log(f'V11 {context}: {job[0]+1}/300; candidate {model["candidate"]}')
                    del block
        del x;gc.collect()
    pd.DataFrame(records).to_csv(A/('preflight_group_statistics.csv' if preflight else 'prediction_group_statistics.csv'),index=False)
    if preflight:
        build.write_json(A/'preflight.json',{'passed':True,'groups':len(records),'model_sha256':build.sha256(A/'reference_model.json'),
            'generator_sha256':build.sha256(Path(__file__)),'methods_sha256':build.sha256(ROOT/'v11_methods.py')})
        return
    writer.close();path.with_suffix('.h5ad.partial').replace(path)
    pd.concat(provenance,ignore_index=True).to_csv(A/'prediction_control_rows.csv',index=False)
    build.write_json(A/'run_parameters.json',{'created_utc':datetime.now(timezone.utc).isoformat(),'seed':SEED,'setting':settings,'response_candidate':model['candidate'],
        'knockdown_assumption':.9,'generator_sha256':build.sha256(Path(__file__)),'methods_sha256':build.sha256(ROOT/'v11_methods.py'),
        'model_sha256':build.sha256(A/'reference_model.json'),'upstream_sha256':atlas_sha,'workers':workers,
        'simulation':'Frozen V11 consensus/amplitude and dispersion setting; guide/depth-stratified unique control samples; synthetic dual-moment integer populations; own target handled separately.'})


def package():
    manifest,genes,targets=build.read_spec(ROOT)
    verified=build.validate_h5ad(OUT/'prediction.h5ad',genes,targets,manifest['contexts'],manifest['cells_per_pert'])
    build.write_json(A/'prediction_validation.json',verified)
    from package_v11_memory_safe import main
    main()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--preflight',action='store_true');parser.add_argument('--package-only',action='store_true');parser.add_argument('--workers',type=int,default=4);parser.add_argument('--wait-for-model',action='store_true');args=parser.parse_args()
    if args.wait_for_model:
        import time
        while not (A/'reference_model.json').exists():time.sleep(20)
    if not 1<=args.workers<=8:parser.error('workers must be 1..8')
    if not args.package_only:generate(args.preflight,args.workers)
    if not args.preflight:package()
