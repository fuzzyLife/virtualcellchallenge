"""V76 streaming counts; reviewed V24 writer/emitter and exact unchanged-block reuse."""
from pathlib import Path
import argparse,json,gc
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd,anndata as ad,h5py
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v76_response_dilution'
A=OUT/'analysis'
from v16_methods import emit
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model
SEED=20261011

def load():
    manifest,genes,targets=build.read_spec(ROOT)
    record=json.loads((A/'reference_model.json').read_text())
    assert build.sha256(A/'reference_model_arrays.npz')==record['model_arrays_sha256']
    assert build.sha256(A/'plan.json')==record['scientific_plan_sha256']
    assert build.sha256(A/'production_plan.json')==record['production_plan_sha256']
    assert build.sha256(A/'count_assessment.json')==record['count_assessment_sha256']
    for name,h in record['dependencies'].items():assert build.sha256(ROOT/name)==h,name
    for r in record['context_original_inputs']:
        assert build.sha256(ROOT/r['file'])==r['sha256']
        assert build.sha256(A/'production_inputs'/(r['context']+'_query.npz'))==r['query_sha256']
    with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
        np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets);effect=z['effects'].copy()
    assert effect.shape==(3,300,18533) and np.isfinite(effect).all()
    return manifest,genes,targets,record,effect

def generate(preflight=False,workers=4):
    from v24_runtime import require_science
    require_science()
    manifest,genes,targets,model,effect=load();atlas,atlas_sha=load_atlas_model()
    path=OUT/'prediction.h5ad';records=[];provenance=[];writer=None
    if not preflight:
        if path.exists() or (OUT/'prediction.vcc').exists():raise FileExistsError('Preserve V76 output')
        obs=pd.DataFrame({'target_gene':pd.Categorical(np.tile(np.repeat(targets,400),3),categories=targets),'context':pd.Categorical(np.repeat(manifest['contexts'],120000),categories=manifest['contexts'])},index=[f'v76_{i:06d}' for i in range(360000)])
        writer=build.SparseH5Writer(path.with_suffix('.h5ad.partial'),obs,genes)
    old_handle=None;old_effect=None;old_sha=None
    for ci,context in enumerate(manifest['contexts']):
        x,guides,library,*_=build.load_control(ROOT,context,genes,manifest);depth,p,b,_=moments(x);prepared=(depth,p,b)
        obj=ad.read_h5ad(ROOT/f'context_{context}.h5ad',backed='r');ids=obj.obs_names.to_numpy();obj.file.close()
        selected=['ACLY','ADNP','VIM','SMARCA5','S100A11','FXYD5','CHKB','MED25'];jobs=[(ti,t) for ti,t in enumerate(targets) if not preflight or t in selected]
        def compute(job):
            ti,target=job;rows=build.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([SEED,ci,ti])));own=int(np.flatnonzero(genes==target)[0]);seed=SEED+ci*100000+ti
            reused=False
            count,diag=emit(x,rows,effect[ci,ti],None,seed,model['emitter'],prepared,atlas)
            libs=np.asarray(count.sum(1)).ravel();np.testing.assert_array_equal(libs,library[rows]);assert len(set(rows))==400
            return count,rows,{'context':context,'target_gene':target,'target_index':ti,'nnz':int(count.nnz),'mean_library':float(libs.mean()),'mean_detected_genes':float(count.getnnz(1).mean()),'reused_unchanged_v24':bool(reused),**diag}
        with threadpool_limits(limits=2),ThreadPoolExecutor(max_workers=workers) as pool:
            pending={i:pool.submit(compute,jobs[i]) for i in range(min(workers,len(jobs)))}
            for i,job in enumerate(jobs):
                count,rows,r=pending.pop(i).result();records.append(r)
                if writer is not None:
                    writer.append(count);provenance.append(pd.DataFrame({'context':context,'target_gene':job[1],'prediction_cell_offset':np.arange(400),'source_control_row':rows,'source_control_cell_id':ids[rows],'ntc_id':guides[rows],'original_total_counts':library[rows],'predicted_total_counts':np.asarray(count.sum(1)).ravel()}))
                if i+workers<len(jobs):pending[i+workers]=pool.submit(compute,jobs[i+workers])
                if preflight or (i+1)%10==0:print('V76 generation',context,i+1,'/',len(jobs),flush=True)
                del count
        del x;gc.collect()
    pd.DataFrame(records).to_csv(A/('preflight_group_statistics.csv' if preflight else 'prediction_group_statistics.csv'),index=False)
    if preflight:
        build.write_json(A/'preflight.json',{'passed':True,'groups':len(records),'model_sha256':build.sha256(A/'reference_model.json'),'generator_sha256':build.sha256(Path(__file__))});return
    writer.close();path.with_suffix('.h5ad.partial').replace(path)
    if old_handle is not None:old_handle.close()
    pd.concat(provenance,ignore_index=True).to_csv(A/'prediction_control_rows.csv',index=False)
    build.write_json(A/'run_parameters.json',{'created_utc':datetime.now(timezone.utc).isoformat(),'seed':SEED,'workers':workers,'model_sha256':build.sha256(A/'reference_model.json'),'generator_sha256':build.sha256(Path(__file__)),'emitter_sha256':build.sha256(ROOT/'v16_methods.py'),'upstream_sha256':atlas_sha,'verified_v24_h5ad_sha256':old_sha,'unchanged_v24_groups_reused':int(sum(r['reused_unchanged_v24'] for r in records))})

def package():
    manifest,genes,targets,record,effect=load();del effect
    result=build.validate_h5ad(OUT/'prediction.h5ad',genes,targets,manifest['contexts'],400);build.write_json(A/'prediction_validation.json',result)
    # Official prep runs separately in the competition environment.
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preflight',action='store_true');p.add_argument('--package-only',action='store_true');p.add_argument('--workers',type=int,default=4);args=p.parse_args()
    if not 1<=args.workers<=6:p.error('workers must be 1..6')
    if not args.package_only:generate(args.preflight,args.workers)
    if not args.preflight:package()
