"""Generate and officially package V5 counts from individual guide/depth-balanced control cells."""
from pathlib import Path
import argparse,gc,importlib.metadata,json,subprocess,sys,time
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from evaluate_atlas_transfer import load_atlas_model
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v5_atlas_calibrated'
A=OUT/'analysis'
SEED=20261005


def desired_profiles(mean,pooled,effect,bulk_delta,own,knockdown=.9):
    mean=np.asarray(mean,dtype=float);pooled=np.asarray(pooled,dtype=float)
    if mean.shape!=pooled.shape or effect.shape!=mean.shape or bulk_delta.shape!=mean.shape:raise ValueError('Moment gene axes differ')
    if not 0<knockdown<1:raise ValueError('Invalid knockdown')
    p=mean*np.exp2(np.clip(effect,-3,3));b=np.expm1(np.maximum(np.log1p(50000*pooled)+np.clip(bulk_delta,-3,3),0))
    p[own]*=1-knockdown;b[own]*=1-knockdown
    if not np.isfinite(p).all() or not np.isfinite(b).all() or p.sum()<=0 or b.sum()<=0:raise ValueError('Invalid desired moments')
    p/=p.sum();b/=b.sum();return p,b


def project_bulk(p,b,depths):
    z=depths/depths.mean();lower=(z.min()+.01*(1-z.min()))*p;upper=(z.max()-.01*(z.max()-1))*p
    if np.ptp(depths)==0:return p.copy()
    lo,hi=0.,1.
    while np.clip(b*hi,lower,upper).sum()<1:
        hi*=2
        if not np.isfinite(hi):raise ValueError('Infeasible bulk support')
    for _ in range(70):
        middle=(lo+hi)/2
        if np.clip(b*middle,lower,upper).sum()<1:lo=middle
        else:hi=middle
    return np.clip(b*((lo+hi)/2),lower,upper)


def emit(model,x,libraries,rows,p,b,seed,iterations=200):
    depth=libraries[rows];template=x[rows].toarray().astype(float)/depth[:,None]
    projected=project_bulk(p,b,depth)
    projection=float(np.abs(projected-b).sum())
    if projection>.03:raise ValueError(f'Bulk projection exceeds documented limit: {projection}')
    with threadpool_limits(limits=2):
        counts=model.dual_moment_counts(template,p,b,depths=depth,seed=int(seed),iterations=iterations)
    np.testing.assert_array_equal(counts.sum(axis=1,dtype=np.int64),depth)
    mean=np.mean(counts/depth[:,None],axis=0);bulk=counts.sum(axis=0,dtype=np.int64)/depth.sum()
    if not np.issubdtype(counts.dtype,np.integer) or (counts<0).any():raise ValueError('Invalid generated integer counts')
    diagnostics={'requested_mean_probability_L1':float(np.abs(mean-p).sum()),'requested_pooled_probability_L1':float(np.abs(bulk-b).sum()),
                 'projected_pooled_probability_L1':float(np.abs(bulk-projected).sum()),'bulk_feasibility_projection_L1':projection,
                 'mean_detected_genes':float(np.count_nonzero(counts,axis=1).mean()),'mean_library':float(depth.mean()),
                 'CPM_total_variance':float(np.var(counts/depth[:,None],axis=0).sum()),'iterations_limit':iterations}
    if diagnostics['requested_mean_probability_L1']>.03 or diagnostics['projected_pooled_probability_L1']>.002:
        raise ValueError('Integer moment error exceeds candidate acceptance limits')
    return sparse.csr_matrix(counts),diagnostics


def inputs():
    manifest,genes,targets=build.read_spec(ROOT);record=json.loads((A/'reference_model.json').read_text())
    assert build.sha256(ROOT/'fit_v5_submission.py')==record['script_sha256']
    assert build.sha256(A/'reference_model_arrays.npz')==record['model_arrays_sha256']
    for original in record['source_original_inputs']:assert build.sha256(ROOT/original['file'])==original['sha256']
    with np.load(A/'reference_model_arrays.npz') as z:
        assert np.array_equal(z['genes'],genes) and np.array_equal(z['targets'],targets)
        effects=z['effects'].copy();bulk=z['bulk_delta'].copy()
    assert np.isfinite(effects).all() and np.isfinite(bulk).all()
    model,model_sha=load_atlas_model();return manifest,genes,targets,record,effects,bulk,model,model_sha


def generate(preflight=False,force=False,workers=4):
    manifest,genes,targets,record,effects,bulk,model,model_sha=inputs();contexts=manifest['contexts'];n=manifest['cells_per_pert']
    path=OUT/'prediction.h5ad';package=OUT/'prediction.vcc'
    if not preflight and not force and (path.exists() or package.exists()):raise FileExistsError('V5 output exists; use --force only for an intentional rebuild')
    chosen_targets=['ACLY','ADNP','VIM','SMARCA5','S100A11','ALG12','ZNF32','MED25']
    writer=None
    if not preflight:
        obs=pd.DataFrame({'target_gene':pd.Categorical(np.tile(np.repeat(targets,n),len(contexts)),categories=targets),
                          'context':pd.Categorical(np.repeat(contexts,len(targets)*n),categories=contexts)},
                          index=[f'v5_{i:06d}' for i in range(len(contexts)*len(targets)*n)])
        writer=build.SparseH5Writer(path.with_suffix('.h5ad.partial'),obs,genes)
    results=[];row_records=[]
    for ci,context in enumerate(contexts):
        x,guides,library,mean_cpm,summary,statistics=build.load_control(ROOT,context,genes,manifest)
        mean=mean_cpm/mean_cpm.sum();pooled=np.asarray(x.sum(axis=0)).ravel().astype(float);pooled/=pooled.sum()
        ids=ad.read_h5ad(ROOT/f'context_{context}.h5ad',backed='r');cell_ids=ids.obs_names.to_numpy();ids.file.close()
        if not preflight:statistics.to_csv(A/f'gene_statistics_{context}.csv',index=False)
        jobs=[(ti,target) for ti,target in enumerate(targets) if not preflight or target in chosen_targets]
        def compute(job):
            ti,target=job
            rng=np.random.default_rng(np.random.SeedSequence([SEED,ci,ti]));rows=build.guide_stratified_sample(guides,library,n,rng)
            own=int(np.flatnonzero(genes==target)[0]);p,b=desired_profiles(mean,pooled,effects[ci,ti],bulk[ci,ti],own)
            block,diagnostics=emit(model,x,library,rows,p,b,SEED+ci*100000+ti)
            result={'context':context,'target_gene':target,'target_index':ti,'nnz':int(block.nnz),
                    'own_requested_mean_ratio':float(p[own]/mean[own]) if mean[own]>0 else 0.,
                    'own_generated_mean_raw':float(block[:,own].mean()),**diagnostics}
            return block,rows,result
        # A bounded window retains at most workers completed groups and writes in panel order.
        with threadpool_limits(limits=2):
            with ThreadPoolExecutor(max_workers=workers) as pool:
                pending={i:pool.submit(compute,jobs[i]) for i in range(min(workers,len(jobs)))}
                for index,job in enumerate(jobs):
                    ti,target=job;block,rows,result=pending.pop(index).result();results.append(result)
                    if writer is not None:
                        writer.append(block)
                        row_records.append(pd.DataFrame({'context':context,'target_gene':target,'prediction_cell_offset':np.arange(n),
                            'source_control_row':rows,'source_control_cell_id':cell_ids[rows],'ntc_id':guides[rows],'original_total_counts':library[rows]}))
                    future=index+workers
                    if future<len(jobs):pending[future]=pool.submit(compute,jobs[future])
                    if preflight or (ti+1)%10==0:build.log(f'V5 {context}: {ti+1}/300; last moment L1 {result["requested_mean_probability_L1"]:.5f}, pooled {result["requested_pooled_probability_L1"]:.5f}')
                    del block
        del x;gc.collect()
    name='preflight_group_statistics.csv' if preflight else 'prediction_group_statistics.csv';pd.DataFrame(results).to_csv(A/name,index=False)
    if preflight:
        build.write_json(A/'preflight.json',{'passed':True,'groups':len(results),'source_model_sha256':build.sha256(A/'reference_model.json'),
                       'upstream_model_sha256':model_sha,'script_sha256':build.sha256(Path(__file__)),
                       'scope':'24 stress-test groups, including new H1 targets, lineage markers and neighboring-promoter targets; checks numerical generation, not hidden response accuracy.'})
        return
    writer.close();path.with_suffix('.h5ad.partial').replace(path)
    pd.concat(row_records,ignore_index=True).to_csv(A/'prediction_control_rows.csv',index=False)
    build.write_json(A/'run_parameters.json',{'created_utc':datetime.now(timezone.utc).isoformat(),'seed':SEED,'knockdown_fraction_assumption':.9,
            'simulation':'Individual control-cell templates, separate mean-CPM/pooled moments, deterministic integer rounding with exact original depths',
            'upstream_model_sha256':model_sha,'generator_sha256':build.sha256(Path(__file__)),'source_model_sha256':build.sha256(A/'reference_model.json'),
            'versions':{name:importlib.metadata.version(name) for name in ['vcc-cli','anndata','numpy','scipy','h5py','pandas']},
            'cell_pool_size':1,'parallel_group_workers':workers,'threads_for_moment_fit':2,'moment_iterations_limit':200})


def package(force=False):
    manifest,genes,targets=build.read_spec(ROOT);path=OUT/'prediction.h5ad';destination=OUT/'prediction.vcc'
    if destination.exists() and not force:raise FileExistsError('Package already exists')
    build.log('Checking all V5 counts and officially packaging')
    verified=build.validate_h5ad(path,genes,targets,manifest['contexts'],manifest['cells_per_pert']);build.write_json(A/'prediction_validation.json',verified)
    pending=destination.with_suffix('.vcc.partial')
    command=[sys.executable,'-m','vcc','--json','prep',str(path),'-g',str(ROOT/'gene_names.csv'),'--perts',str(ROOT/'pert_counts.csv'),
              '--contexts',','.join(manifest['contexts']),'--cells-per-pert',str(manifest['cells_per_pert']),'-o',str(pending)]
    if force:command.append('--force')
    result=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (A/'official_prep_stdout.json').write_text(result.stdout,encoding='utf-8');(A/'official_prep_stderr.txt').write_text(result.stderr,encoding='utf-8')
    if result.returncode:raise RuntimeError('Official preparation failed: '+result.stderr[-2000:])
    from vcc.vccfile import validate_vcc,read_vcc_meta
    validate_vcc(pending);metadata=read_vcc_meta(pending);assert metadata['nnz']==verified['nnz'];pending.replace(destination)
    build.write_json(A/'submission_ready.json',{'file':str(destination),'bytes':destination.stat().st_size,'sha256':build.sha256(destination),
          'official_prep_passed':True,'package_meta':metadata,'prediction_validation':verified,'uploaded':False,'leaderboard_score':None,
          'generator_sha256':build.sha256(Path(__file__)),'model_sha256':build.sha256(A/'reference_model.json')})
    build.log(f'READY {destination} ({destination.stat().st_size/1024**3:.2f} GiB)')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--preflight',action='store_true');parser.add_argument('--package-only',action='store_true');parser.add_argument('--force',action='store_true');parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    if not 1<=args.workers<=8:parser.error('--workers must be in [1, 8]')
    A.mkdir(parents=True,exist_ok=True)
    if not args.package_only:generate(args.preflight,args.force,args.workers)
    if not args.preflight:package(args.force)
