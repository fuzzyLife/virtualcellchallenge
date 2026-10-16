"""Build the selected V35/V36 submission, preserving the actual V24 parent."""
import argparse,gc,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as b
import v35_methods as m
import v36_methods as pop
from v6_methods import moments
from v35_contrast_methods import transform
from evaluate_atlas_transfer import load_atlas_model
ROOT=m.ROOT;OUT=m.P36.parent;A=m.P36
SEED=20261011

def fit_production():
 m.setup();m.check();selection=json.loads((A/'selection.json').read_text());version=selection['model_version'];net,g,model=m.load('all',version);assert model['plan_sha256']==m.sha(m.A/'plan.json') and model['outer']=='all';assert set(model['sources'])==set(m.old.old.SOURCES)
 from build_v24_submission import load as parent_load
 manifest,genes,targets,parent,effects=parent_load();result=[];dispersion=[];records=[];variance_model=None;variance_net=None
 if selection['selected_arm']=='v35_population_dispersion':variance_net,variance_geometry,variance_model=m.load('all',36)
 for ci,c in enumerate(manifest['contexts']):
  x,guides,libs,*_=b.load_control(ROOT,c,genes,manifest);ctrl=m.control_stats(x);query=ctrl[0]*1e6;p,v,trace=m.predict(net,g,query,[],base=effects[ci],control=ctrl if version==36 else None);p=transform(p,effects[ci],m.predict(net,g,query,[],base=effects[ci],shuffle=True)[0],genes,targets,g['axis'],trace,query,selection['selected_arm']) if selection['selected_arm'] in ['target_centered25','retrieval_contrast25'] else p;v=m.predict(variance_net,variance_geometry,query,[],base=effects[ci],control=ctrl)[1] if variance_net is not None else v;outside=np.ones(len(genes),bool);outside[g['axis']]=False;np.testing.assert_array_equal(p[:,outside],effects[ci][:,outside]);assert float(abs(p-effects[ci]).max())<=.500002;result.append(p);dispersion.append(v);records.append({'context':c,'changed_targets':int((abs(p-effects[ci]).sum(1)>0).sum()),'max_LFC_correction':float(abs(p-effects[ci]).max()),'query_sha256':__import__('hashlib').sha256(query.tobytes()).hexdigest(),'control_rows':x.shape[0],'retrieval_trace':trace});del x;gc.collect();print('Production retrieval fitted',c,records[-1]['changed_targets'],flush=True)
 file=A/'reference_model_arrays.npz';assert not file.exists();np.savez_compressed(file,effects=np.array(result),dispersion=np.array(dispersion),genes=genes,targets=targets)
 deps=['v35_contrast_methods.py','v35_methods.py','v36_methods.py','build_v35_v36_submission.py','v16_methods.py','build_submission.py']
 b.write_json(A/'reference_model.json',{'status':selection['status'],'confirmed_improvement':False,'available_public_guards_passed':selection['available_public_guards_passed'],'model_version':version,'selected_arm':selection['selected_arm'],'arrays_sha256':b.sha256(file),'fit_model_sha256':model['sha256'],'dispersion_fit_model_sha256':variance_model['sha256'] if variance_model else None,'plan_sha256':b.sha256(m.A/'plan.json'),'selection_sha256':b.sha256(A/'selection.json'),'parent_V24_arrays_sha256':b.sha256(ROOT/'experiments/v24_cd4_response/analysis/reference_model_arrays.npz'),'dependencies':{n:b.sha256(ROOT/n) for n in deps},'context_diagnostics':records,'reserve_opened':False,'uploaded':False})

def load():
 manifest,genes,targets=b.read_spec(ROOT);rec=json.loads((A/'reference_model.json').read_text());assert b.sha256(A/'reference_model_arrays.npz')==rec['arrays_sha256'] and b.sha256(A/'selection.json')==rec['selection_sha256']
 for name,h in rec['dependencies'].items():assert b.sha256(ROOT/name)==h,name
 with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
  np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets);effect=z['effects'].copy();dispersion=z['dispersion'].copy()
 assert effect.shape==(3,300,18533) and np.isfinite(effect).all() and np.isfinite(dispersion).all();rec['dispersion']=dispersion
 return manifest,genes,targets,rec,effect

def generate(preflight=False,workers=4):
 manifest,genes,targets,model,effects=load();atlas,atlas_sha=load_atlas_model();path=OUT/'prediction.h5ad';writer=None;records=[];provenance=[]
 if not preflight:
  assert not path.exists() and not (OUT/'prediction.vcc').exists();obs=pd.DataFrame({'target_gene':pd.Categorical(np.tile(np.repeat(targets,400),3),categories=targets),'context':pd.Categorical(np.repeat(manifest['contexts'],120000),categories=manifest['contexts'])},index=[f'v36_{i:06d}' for i in range(360000)]);writer=b.SparseH5Writer(path.with_suffix('.h5ad.partial'),obs,genes)
 for ci,c in enumerate(manifest['contexts']):
  x,guides,library,*_=b.load_control(ROOT,c,genes,manifest);depth,p,bulk,_=moments(x);obj=ad.read_h5ad(ROOT/f'context_{c}.h5ad',backed='r');ids=obj.obs_names.to_numpy();obj.file.close();chosen=['ACLY','ADNP','VIM','SMARCA5','S100A11','FXYD5','CHKB','MED25'];jobs=[(ti,t) for ti,t in enumerate(targets) if not preflight or t in chosen]
  def compute(job):
   ti,t=job;rows=b.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([SEED,ci,ti])));own=int(np.flatnonzero(genes==t)[0]);seed=SEED+ci*100000+ti;dd=model['dispersion'][ci,ti] if model['selected_arm'] in ['v36_dispersion','v35_population_dispersion'] else None
   count,diag=pop.emit(x[rows],effects[ci,ti],own,seed,(depth[rows],p,bulk),atlas,dd);np.testing.assert_array_equal(count.sum(1).A1,library[rows]);assert len(set(rows))==400 and count.dtype.kind in 'iu' and count.data.max(initial=0)<=1000000;assert library[rows].max()<=1000000;assert diag['mean_probability_L1']<=.03 and diag['bulk_probability_L1']<=.031
   return count,rows,{'context':c,'target_gene':t,'target_index':ti,'nnz':int(count.nnz),'mean_library':float(library[rows].mean()),'mean_detected_genes':float(count.getnnz(1).mean()),**diag}
  with threadpool_limits(limits=2),ThreadPoolExecutor(max_workers=workers) as pool:
   pending={i:pool.submit(compute,jobs[i]) for i in range(min(workers,len(jobs)))}
   for i,job in enumerate(jobs):
    count,rows,r=pending.pop(i).result();records.append(r)
    if writer is not None:
     writer.append(count);provenance.append(pd.DataFrame({'context':c,'target_gene':job[1],'prediction_cell_offset':np.arange(400),'source_control_row':rows,'source_control_cell_id':ids[rows],'ntc_id':guides[rows],'original_total_counts':library[rows],'predicted_total_counts':count.sum(1).A1}))
    if i+workers<len(jobs):pending[i+workers]=pool.submit(compute,jobs[i+workers])
    if preflight or (i+1)%10==0:print('V35/V36 generation',c,i+1,'/',len(jobs),'nnz',r['nnz'],flush=True)
    del count
  del x;gc.collect()
 pd.DataFrame(records).to_csv(A/('preflight_group_statistics.csv' if preflight else 'prediction_group_statistics.csv'),index=False)
 if preflight:
  b.write_json(A/'preflight.json',{'passed':True,'groups':len(records),'max_mean_probability_L1':max(r['mean_probability_L1'] for r in records),'model_sha256':b.sha256(A/'reference_model.json'),'generator_sha256':b.sha256(Path(__file__))});return
 writer.close();path.with_suffix('.h5ad.partial').replace(path);pd.concat(provenance,ignore_index=True).to_csv(A/'prediction_control_rows.csv',index=False);b.write_json(A/'run_parameters.json',{'seed':SEED,'workers':workers,'generator_sha256':b.sha256(Path(__file__)),'model_sha256':b.sha256(A/'reference_model.json'),'upstream_sha256':atlas_sha,'same_guide_stratified_rows_and_seed_as_V24':True,'uploaded':False});validation=b.validate_h5ad(path,genes,targets,manifest['contexts'],400);b.write_json(A/'prediction_validation.json',validation);print('Raw submission validated',validation,flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['fit','preflight','generate']);p.add_argument('--workers',type=int,default=4);args=p.parse_args()
 if not 1<=args.workers<=6:p.error('workers1..6')
 with threadpool_limits(limits=2):fit_production() if args.stage=='fit' else generate(args.stage=='preflight',args.workers)