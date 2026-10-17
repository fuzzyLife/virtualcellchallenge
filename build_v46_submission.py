"""Build V46 only after explicit public development/confirmation assessment."""
import argparse,gc,json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from threadpoolctl import threadpool_limits
import build_submission as b
import v46_population_mixture as m
from v46_thread_cache import PerThreadCache
m._CACHE=PerThreadCache()
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model
SEED=20261011
def fit():
 m.check();selection=m.read(m.A/'release_decision.json');assert selection['release'] and selection['all_public_guards_passed']
 from build_v24_submission import load as parent
 manifest,genes,targets,parent_record,base=parent();chosen=selection['selected'];assert chosen in m.ARMS
 path=m.A/'reference_model_arrays.npz';assert not path.exists()
 np.savez_compressed(path,effects=base.copy(),genes=np.asarray(genes,str),targets=np.asarray(targets,str))
 deps=['v46_thread_cache.py','v44_response_dispersion.py','v46_population_mixture.py','build_v46_submission.py','v38_response_hurdle.py','v16_methods.py','build_submission.py','verify_vcc_payload.py']
 b.write_json(m.A/'reference_model.json',{'status':'New V24-mean response-dispersion candidate; official leaderboard performance unknown','confirmed_improvement':False,'available_public_guards_passed':selection['all_public_guards_passed'],'selected_arm':chosen,'arrays_sha256':b.sha256(path),'selection_sha256':b.sha256(m.A/'release_decision.json'),'plan_sha256':b.sha256(m.A/'plan.json'),'parent_V24_arrays_sha256':b.sha256(m.ROOT/'experiments/v24_cd4_response/analysis/reference_model_arrays.npz'),'dependencies':{s:b.sha256(m.ROOT/s) for s in deps},'mean_effects_bitwise_V24':True,'novelty':'Whole cell vectors are drawn from compact V24 and residual20 populations with balanced seeded assignment. Mean weights unchanged; count distributions differ. Per-thread cache isolation is execution-only.','uploaded':False})
def load():
 manifest,genes,targets=b.read_spec(m.ROOT);rec=m.read(m.A/'reference_model.json')
 assert b.sha256(m.A/'reference_model_arrays.npz')==rec['arrays_sha256'];assert b.sha256(m.A/'release_decision.json')==rec['selection_sha256']
 for name,h in rec['dependencies'].items():assert b.sha256(m.ROOT/name)==h,name
 with np.load(m.A/'reference_model_arrays.npz',allow_pickle=False) as z:
  np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets);effects=z['effects'].copy()
 assert effects.shape==(3,300,18533) and np.isfinite(effects).all()
 return manifest,genes,targets,rec,effects
def generate(preflight=False,workers=4):
 manifest,genes,targets,model,effects=load();atlas,atlas_sha=load_atlas_model();path=m.OUT/'prediction.h5ad';writer=None;records=[];provenance=[]
 if not preflight:
  assert not path.exists() and not (m.OUT/'prediction.vcc').exists()
  obs=pd.DataFrame({'target_gene':pd.Categorical(np.tile(np.repeat(targets,400),3),categories=targets),'context':pd.Categorical(np.repeat(manifest['contexts'],120000),categories=manifest['contexts'])},index=[f'v46_{i:06d}' for i in range(360000)])
  writer=b.SparseH5Writer(path.with_suffix('.h5ad.partial'),obs,genes)
 for ci,c in enumerate(manifest['contexts']):
  x,guides,library,*_=b.load_control(m.ROOT,c,genes,manifest);depth,p,bulk,_=moments(x)
  obj=ad.read_h5ad(m.ROOT/f'context_{c}.h5ad',backed='r');ids=obj.obs_names.to_numpy();obj.file.close()
  chosen=['ACLY','ADNP','VIM','SMARCA5','S100A11','FXYD5','CHKB','MED25'];jobs=[(ti,t) for ti,t in enumerate(targets) if not preflight or t in chosen]
  def compute(job):
   ti,t=job;rows=b.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([SEED,ci,ti])));own=int(np.flatnonzero(genes==t)[0]);seed=SEED+ci*100000+ti
   count,diag=m.emit(x[rows],effects[ci,ti],own,seed,(depth[rows],p,bulk),atlas,model['selected_arm'])
   np.testing.assert_array_equal(count.sum(1).A1,library[rows]);assert len(set(rows))==400 and count.dtype.kind in 'iu' and count.data.max(initial=0)<=1000000
   assert library[rows].max()<=1000000 and diag['mean_probability_L1']<=.03 and diag['bulk_probability_L1']<=.031
   return count,rows,{'context':c,'target_gene':t,'target_index':ti,'nnz':int(count.nnz),'mean_library':float(library[rows].mean()),'mean_detected_genes':float(count.getnnz(1).mean()),**diag}
  with threadpool_limits(limits=2),ThreadPoolExecutor(max_workers=workers) as pool:
   pending={i:pool.submit(compute,jobs[i]) for i in range(min(workers,len(jobs)))}
   for i,job in enumerate(jobs):
    count,rows,r=pending.pop(i).result();records.append(r)
    if writer is not None:
     writer.append(count);provenance.append(pd.DataFrame({'context':c,'target_gene':job[1],'prediction_cell_offset':np.arange(400),'source_control_row':rows,'source_control_cell_id':ids[rows],'ntc_id':guides[rows],'original_total_counts':library[rows],'predicted_total_counts':count.sum(1).A1}))
    if i+workers<len(jobs):pending[i+workers]=pool.submit(compute,jobs[i+workers])
    if preflight or (i+1)%10==0:print('V46 generation',c,i+1,'/',len(jobs),'nnz',r['nnz'],flush=True)
    del count
  del x;gc.collect()
 pd.DataFrame(records).to_csv(m.A/('preflight_group_statistics.csv' if preflight else 'prediction_group_statistics.csv'),index=False)
 if preflight:
  b.write_json(m.A/'preflight.json',{'passed':True,'groups':len(records),'model_sha256':b.sha256(m.A/'reference_model.json'),'generator_sha256':b.sha256(Path(__file__))});return
 writer.close();path.with_suffix('.h5ad.partial').replace(path);pd.concat(provenance,ignore_index=True).to_csv(m.A/'prediction_control_rows.csv',index=False)
 b.write_json(m.A/'run_parameters.json',{'seed':SEED,'workers':workers,'generator_sha256':b.sha256(Path(__file__)),'model_sha256':b.sha256(m.A/'reference_model.json'),'upstream_sha256':atlas_sha,'same_guide_stratified_rows_and_seed_as_V24':True,'uploaded':False})
 result=b.validate_h5ad(path,genes,targets,manifest['contexts'],400);b.write_json(m.A/'prediction_validation.json',result);print('V46 full raw submission validated',result,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['fit','preflight','generate']);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
 if not 1<=a.workers<=6:p.error('workers must be 1..6')
 with threadpool_limits(limits=2):fit() if a.stage=='fit' else generate(a.stage=='preflight',a.workers)






