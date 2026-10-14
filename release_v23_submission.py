"""Release a passing V23 predictor, or an explicitly byte-identical V16 fallback."""
from pathlib import Path
import argparse,json,gc,shutil
from concurrent.futures import ProcessPoolExecutor
import numpy as np,pandas as pd,anndata as ad
pd.options.future.infer_string=False
from threadpoolctl import threadpool_limits
from v23_methods import ROOT,OUT,A,load_plan,load_model,query_states,emit,write
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model
import build_submission as b
SEED=20261011

def fallback():
 selection=json.loads((A/'selection.json').read_text());assert not selection['new_candidate_passes'];audit=json.loads((A/'independent_experiment_audit.json').read_text());assert audit['passed']
 old=ROOT/'experiments/v16_control_preserving';ready=json.loads((old/'analysis/submission_ready.json').read_text());validation=json.loads((old/'analysis/independent_validation.json').read_text());assert validation['passed'] and ready['official_prep_passed'] and ready['sha256']==validation['package_sha256'];source=old/'prediction.vcc';assert source.stat().st_size==ready['bytes'] and b.sha256(source)==ready['sha256']
 dest=OUT/'prediction.vcc'
 if dest.exists():raise FileExistsError('Preserve release package')
 pending=dest.with_suffix('.vcc.partial');shutil.copyfile(source,pending);assert pending.stat().st_size==ready['bytes'] and b.sha256(pending)==ready['sha256'];pending.replace(dest)
 write(A/'submission_ready.json',{'file':str(dest),'bytes':dest.stat().st_size,'sha256':ready['sha256'],'upload_ready':True,'official_prep_passed':True,'package_meta':ready['package_meta'],'prediction_validation':ready['prediction_validation'],'new_predictor_released':False,'byte_identical_to':'experiments/v16_control_preserving/prediction.vcc','published_reference_entry':'wG3GGbaTsl2dXFn7cQy1','published_reference_overall':.1416138568460195,'recommend_upload':False,'reason':'Novel V23 predictor failed frozen generated-count guards. This is the verified champion fallback; reuploading identical predictions adds no model improvement.','selection_sha256':b.sha256(A/'selection.json'),'experiment_audit_sha256':b.sha256(A/'independent_experiment_audit.json'),'prior_full_payload_validation_sha256':b.sha256(old/'analysis/independent_validation.json'),'release_script_sha256':b.sha256(Path(__file__)),'uploaded':False});print('Upload-ready V16 champion fallback copied and hash-verified; no new-model or improvement claim.',flush=True)

def fit():
 plan=load_plan();selection=json.loads((A/'selection.json').read_text());confirmation=json.loads((A/'confirmation_result.json').read_text());audit=json.loads((A/'independent_experiment_audit.json').read_text());assert selection['new_candidate_passes'] and confirmation['new_candidate_passes'] and audit['passed'];assert selection['candidate']==confirmation['candidate']=='mean_preserved';assert not (A/'reference_model.json').exists()
 manifest,genes,targets=b.read_spec(ROOT);model=load_model('all');np.testing.assert_array_equal(model['genes'],genes);np.testing.assert_array_equal(model['targets'],targets)
 old=ROOT/'experiments/v16_control_preserving/analysis';prior=json.loads((old/'reference_model.json').read_text());assert b.sha256(old/'reference_model_arrays.npz')==prior['model_arrays_sha256'];shutil.copyfile(old/'reference_model_arrays.npz',A/'reference_model_arrays.npz')
 write(A/'reference_model.json',{'model':'V23 bounded state-dependent continuous-response population model','candidate':'mean_preserved','emitter':'mean_preserved','status':'exploratory_public_algorithm_verification','available_public_guards_passed':True,'confirmed_improvement':bool(selection['complete_primary_metrics'] and confirmation['complete_primary_metrics']),'model_arrays_sha256':b.sha256(A/'reference_model_arrays.npz'),'factor_arrays_sha256':b.sha256(A/'models/all.npz'),'source_model_metadata_sha256':b.sha256(A/'models/all.json'),'selection_sha256':b.sha256(A/'selection.json'),'confirmation_sha256':b.sha256(A/'confirmation_result.json'),'audit_sha256':b.sha256(A/'independent_experiment_audit.json'),'context_original_inputs':[{'context':c,'file':f'context_{c}.h5ad','sha256':b.sha256(ROOT/f'context_{c}.h5ad')} for c in manifest['contexts']],'dependencies':{n:b.sha256(ROOT/n) for n in set(plan['dependencies'])|{'release_v23_submission.py','package_v23_submission.py','validate_v23_submission.py','unpack_v23_official_payload.py'}},'scope':'Mean responses retained from V16; public held-family and CD4 raw-cell reserve, not private score guarantee.'})

def load():
 manifest,genes,targets=b.read_spec(ROOT);record=json.loads((A/'reference_model.json').read_text())
 for name,h in record['dependencies'].items():assert b.sha256(ROOT/name)==h
 for item in record['context_original_inputs']:assert b.sha256(ROOT/item['file'])==item['sha256']
 for name,key in [('reference_model_arrays.npz','model_arrays_sha256'),('models/all.npz','factor_arrays_sha256'),('models/all.json','source_model_metadata_sha256')]:assert b.sha256(A/name)==record[key]
 with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as z:
  np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets);effect=z['effects'].copy()
 return manifest,genes,targets,record,effect

def generate(preflight,workers):
 manifest,genes,targets,record,effects=load();model=load_model('all');atlas,_=load_atlas_model();path=OUT/'prediction.h5ad';stats=[];provenance=[];writer=None
 if not preflight:
  check=json.loads((A/'preflight.json').read_text());assert check['passed'] and check['model_sha256']==b.sha256(A/'reference_model.json')
  if path.exists() or (OUT/'prediction.vcc').exists():raise FileExistsError('Preserve V23 output')
  obs=pd.DataFrame({'target_gene':pd.Categorical(np.tile(np.repeat(targets,400),3),categories=targets),'context':pd.Categorical(np.repeat(manifest['contexts'],120000),categories=manifest['contexts'])},index=[f'v23_{i:06d}' for i in range(360000)]);writer=b.SparseH5Writer(path.with_suffix('.h5ad.partial'),obs,genes)
 for ci,context in enumerate(manifest['contexts']):
  x,guides,library,*_=b.load_control(ROOT,context,genes,manifest);depth,p,bulk,_=moments(x);states=query_states(x,model);obj=ad.read_h5ad(ROOT/f'context_{context}.h5ad',backed='r');ids=obj.obs_names.to_numpy();obj.file.close();jobs=[(ti,t) for ti,t in enumerate(targets) if not preflight or t in ['ACLY','ADNP','VIM','SMARCA5','S100A11','FXYD5','CHKB','MED25']]
  def payload(job):
   ti,target=job;rows=b.guide_stratified_sample(guides,library,400,np.random.default_rng(np.random.SeedSequence([SEED,ci,ti])));own=int(np.flatnonzero(genes==target)[0]);return (x[rows],effects[ci,ti],own,SEED+ci*100000+ti,(depth[rows],p,bulk),model,ti,states[rows],rows,context,str(target))
  with ProcessPoolExecutor(max_workers=workers) as pool:
   pending={i:pool.submit(produce,payload(jobs[i])) for i in range(min(workers,len(jobs)))}
   for i,job in enumerate(jobs):
    count,rows,rec=pending.pop(i).result();stats.append(rec)
    if writer is not None:
     writer.append(count);provenance.append(pd.DataFrame({'context':context,'target_gene':job[1],'prediction_cell_offset':np.arange(400),'source_control_row':rows,'source_control_cell_id':ids[rows],'ntc_id':guides[rows],'original_total_counts':library[rows],'predicted_total_counts':np.asarray(count.sum(1)).ravel()}))
    if i+workers<len(jobs):pending[i+workers]=pool.submit(produce,payload(jobs[i+workers]))
    if preflight or (i+1)%10==0:print('V23 generation',context,i+1,'/',len(jobs),flush=True)
  del x;gc.collect()
 pd.DataFrame(stats).to_csv(A/('preflight_groups.csv' if preflight else 'prediction_group_statistics.csv'),index=False)
 if preflight:write(A/'preflight.json',{'passed':True,'groups':len(stats),'model_sha256':b.sha256(A/'reference_model.json')});return
 writer.close();path.with_suffix('.h5ad.partial').replace(path);pd.concat(provenance,ignore_index=True).to_csv(A/'prediction_control_rows.csv',index=False);write(A/'prediction_validation.json',b.validate_h5ad(path,genes,targets,manifest['contexts'],400))

_WORKER_ATLAS=None
def produce(args):
 global _WORKER_ATLAS
 template,effect,own,seed,prepared,model,ti,z,rows,context,target=args
 with threadpool_limits(limits=2):
  if _WORKER_ATLAS is None:_WORKER_ATLAS=load_atlas_model()[0]
  count,diag=emit(template,effect,own,seed,'mean_preserved',prepared,_WORKER_ATLAS,model,ti,z)
  np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),prepared[0]);factor=np.exp2(np.clip(effect,-3,3));factor[own]=.1;wanted=prepared[1]*factor;wanted/=wanted.sum();actual=np.asarray(count.astype(float).multiply((1/np.asarray(prepared[0]))[:,None]).mean(0)).ravel();diag['final_mean_probability_L1']=float(abs(actual-wanted).sum())
  return count,rows,{'context':context,'target_gene':target,'target_index':ti,'nnz':int(count.nnz),**diag}


def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['fallback','fit','preflight','build']);p.add_argument('--workers',type=int,default=8);args=p.parse_args()
 with threadpool_limits(limits=2):
  if args.stage=='fallback':fallback()
  elif args.stage=='fit':fit()
  else:
   assert 1<=args.workers<=8
   generate(args.stage=='preflight',args.workers)
if __name__=='__main__':main()
