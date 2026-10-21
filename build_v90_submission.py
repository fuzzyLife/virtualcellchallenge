"""V90 streaming generation using the unchanged V66 residual10 emitter and template seeds."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,gc
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from v57_pie import read,write,sha
from v16_methods import emit
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'experiments/v90_interval_diagnostic';A=OUT/'analysis';SEED=20261011

def load():
 manifest,genes,targets=build.read_spec(ROOT);r=read(A/'reference_model.json');plan=read(A/'production_plan.json');assert sha(A/'reference_model_arrays.npz')==r['model_arrays_sha256'];assert sha(A/'production_plan.json')==r['production_plan_sha256'];assert sha(A/'release_decision.json')==r['decision_sha256'];assert read(A/'release_decision.json')['experimental_release_allowed']
 for n,h in r['dependencies'].items():assert sha(ROOT/n)==h,n
 for row in r['context_original_inputs']:assert sha(A/'production_inputs'/(row['context']+'_query.npz'))==row['query_sha256']
 with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets);effects=z['effects'].copy()
 assert effects.shape==(3,300,18533) and np.isfinite(effects).all();return manifest,genes,targets,r,effects

def generate(preflight=False,workers=4):
 manifest,genes,targets,model,effects=load();atlas,up=load_atlas_model();records=[];provenance=[];writer=None;path=OUT/'prediction.h5ad'
 if not preflight:
  assert not path.exists() and not path.with_suffix('.vcc').exists() and not path.with_suffix('.h5ad.partial').exists();obs=pd.DataFrame({'target_gene':pd.Categorical(np.tile(np.repeat(targets,400),3),categories=targets),'context':pd.Categorical(np.repeat(manifest['contexts'],120000),categories=manifest['contexts'])},index=[f'v90_{i:06d}' for i in range(360000)]);writer=build.SparseH5Writer(path.with_suffix('.h5ad.partial'),obs,genes)
 for ci,c in enumerate(manifest['contexts']):
  source=ROOT/('context_'+c+'.h5ad');expected=next(x for x in model['context_original_inputs'] if x['context']==c);assert sha(source)==expected['sha256'];x,guides,library,*_=build.load_control(ROOT,c,genes,manifest);depth,cp,bulk,_=moments(x);prepared=(depth,cp,bulk)
  with np.load(A/'production_inputs'/(c+'_query.npz'),allow_pickle=False) as z:np.testing.assert_array_equal(z['library'],library);np.testing.assert_array_equal(z['query'],cp*1e6);np.testing.assert_array_equal(z['bulk'],bulk)
  obj=ad.read_h5ad(source,backed='r');ids=obj.obs_names.to_numpy();obj.file.close();jobs=[(ti,t) for ti,t in enumerate(targets) if not preflight or ti in [0,50,150,299]]
  def compute(job):
   ti,t=job;rows=build.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([SEED,ci,ti])));own=int(np.flatnonzero(genes==t)[0]);seed=SEED+ci*100000+ti;count,diag=emit(x,rows,effects[ci,ti],own,seed,'residual10',prepared,atlas);np.testing.assert_array_equal(count.sum(1).A1,library[rows]);assert len(np.unique(rows))==400;return count,rows,dict(context=c,target_gene=str(t),target_index=ti,nnz=int(count.nnz),**diag)
  with threadpool_limits(limits=2),ThreadPoolExecutor(max_workers=workers) as pool:
   pending={i:pool.submit(compute,jobs[i]) for i in range(min(workers,len(jobs)))}
   for i,job in enumerate(jobs):
    count,rows,r=pending.pop(i).result();records.append(r)
    if writer is not None:writer.append(count);provenance.append(pd.DataFrame(dict(context=c,target_gene=str(job[1]),prediction_cell_offset=np.arange(400),source_control_row=rows,source_control_cell_id=ids[rows],ntc_id=guides[rows],original_total_counts=library[rows],predicted_total_counts=count.sum(1).A1)))
    if i+workers<len(jobs):pending[i+workers]=pool.submit(compute,jobs[i+workers])
    if preflight or (i+1)%10==0:print('V90 generation',c,i+1,'/',len(jobs),flush=True)
    del count
  del x;gc.collect()
 pd.DataFrame(records).to_csv(A/('preflight_group_statistics.csv' if preflight else 'prediction_group_statistics.csv'),index=False)
 if preflight:write(A/'preflight.json',dict(passed=True,groups=len(records)));return
 writer.close();path.with_suffix('.h5ad.partial').replace(path);pd.concat(provenance,ignore_index=True).to_csv(A/'prediction_control_rows.csv',index=False);write(A/'run_parameters.json',dict(seed=SEED,workers=workers,generator_sha256=sha(__file__),emitter_sha256=sha(ROOT/'v16_methods.py'),model_sha256=sha(A/'reference_model.json'),upstream_sha256=up));write(A/'prediction_validation.json',build.validate_h5ad(path,genes,targets,manifest['contexts'],400))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--preflight',action='store_true');a=p.parse_args();assert 1<=a.workers<=8;generate(a.preflight,a.workers)
