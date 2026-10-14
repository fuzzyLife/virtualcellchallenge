"""Reproject unchanged Stack outputs into a compatible training-only program vocabulary."""
from pathlib import Path
import os,sys,json,argparse,time,hashlib,gc
ROOT=Path(__file__).resolve().parent;R=ROOT/'reference/v20_stack';A=ROOT/'experiments/v21_stack_programs/analysis';P=ROOT/'experiments/v20_stack_transfer/analysis';os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8');sys.path.insert(0,str(R/'code/src'))
import numpy as np,torch
from scipy import sparse
from stack.models.finetune import ICL_FinetunedModel
from stack_v20_mean_adapter import predict_means
from v20_methods import summarize_pair,align_counts
from fit_v7_response import raw_inputs,read_cache

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for v in iter(lambda:f.read(8*1024**2),b''):h.update(v)
 return h.hexdigest()
def write(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
def prompt(source,target,ti,plan,genes,vocab):
 for folder in [P/'inputs/sources',A/'inputs/sources']:
  side=folder/f'{source}_{target}.json'
  if side.exists():return side,json.loads(side.read_text())
 folder=A/'inputs/sources';folder.mkdir(exist_ok=True);x,meta,measured,_=raw_inputs(source,genes);ctx=next(c for c in plan['contexts'] if c['family']==source);left=np.asarray(ctx['left']);labels=meta.gene_target.to_numpy(str);batches=meta['sample'].to_numpy(str);control={v:left[batches[left]==v] for v in sorted(set(batches[left]))};rng=np.random.default_rng(np.random.SeedSequence([20020,plan['family_order'].index(source),ti]));rows=np.flatnonzero((labels==target)&np.isin(batches,list(control)));rows=rng.permutation(rows)[:400];assert len(rows)>=40;pool={k:list(rng.permutation(v)) for k,v in control.items()};paired=[]
 for row in rows:
  key=batches[row]
  if not pool[key]:pool[key]=list(rng.permutation(control[key]))
  paired.append(pool[key].pop())
 paired=np.asarray(paired);assert np.array_equal(batches[rows],batches[paired]);records={}
 for kind,ids in [('treated',rows),('control',paired)]:
  f=folder/f'{source}_{target}_{kind}.npz';sparse.save_npz(f,align_counts(x[ids],genes,vocab));records[kind]={'rows':ids.tolist(),'unique_original_cells':len(set(ids)),'sha256':sha(f),'full_library_mean':float(meta.total_counts.iloc[ids].mean())}
 rec={'family':source,'target':target,'records':records,'vocabulary_measured_mask':np.isin(vocab,genes[measured]).tolist(),'batches_match_by_row':True,'source_input_rows':len(rows)};side=folder/f'{source}_{target}.json';write(side,rec);del x;gc.collect();return side,rec

def main(stage,family):
 torch.set_num_threads(4);torch.manual_seed(20020);torch.use_deterministic_algorithms(True);plan=json.loads((A/'plan.json').read_text())
 for name,h in plan['dependencies'].items():assert sha(ROOT/name)==h
 assert sha(P/'plan.json')==plan['parent_plan_sha256']
 if stage=='confirmation':assert json.loads((A/'selection.json').read_text())['new_candidate_passes']
 assert sha(R/'bc_large_aligned.ckpt')=='f93cf6f42f36c8a85dc570d92e801c1fc3e1f45d55741e6bb250d178a6b6ad36';ck=torch.load(R/'bc_large_aligned.ckpt',map_location='cpu',weights_only=False);model=ICL_FinetunedModel(**ck['hyper_parameters']['model_config']);model.load_state_dict({k[6:]:v for k,v in ck['state_dict'].items()},strict=True);del ck;model.to('cuda:0').eval();vocab=json.loads((R/'genes.json').read_text());info=json.loads((A/'inputs'/f'{family}.json').read_text());qfile=A/'inputs'/f'{family}_query.npz';assert sha(qfile)==info['query_arrays_sha256']
 with np.load(qfile,allow_pickle=False) as z:q={k:z[k].copy() for k in z.files}
 query=sparse.load_npz(P/'inputs'/f'{family}_counts.npz');targets=q['targets'];genes=q['genes'];ctx=next(c for c in plan['contexts'] if c['family']==family);needed=sorted({int(np.flatnonzero(targets==r['target'])[0]) for r in ctx['cells'] if r['split']==stage});g=len(genes);parent_file=P/'experts'/f'{family}_{stage}.npz';original=None;original_records={}
 if parent_file.exists():
  original_side=json.loads(parent_file.with_suffix('.json').read_text());assert sha(parent_file)==original_side['sha256']
  with np.load(parent_file,allow_pickle=False) as z:original={k:z[k].copy() for k in z.files}
  original_records={r['target']:r for r in original_side['records']}
 arrays={k:np.zeros((300,g),np.float32) for k in ['lfc','shrink','raw','shuffled_lfc','shuffled_shrink']};arrays.update(supported=np.zeros((300,g),bool),program_delta=np.zeros((300,400,16),np.float32),shuffled_program_delta=np.zeros((300,400,16),np.float32),direction=np.zeros((300,16),np.float32),valid=np.zeros(300,bool),genes=genes,targets=targets)
 if original is not None:
  for k in ['lfc','shrink','raw','shuffled_lfc','shuffled_shrink','supported']:arrays[k]=original[k].copy()
 file=A/'experts'/f'{family}_{stage}.npz';file.parent.mkdir(exist_ok=True);records=[]
 if file.exists():
  side=json.loads(file.with_suffix('.json').read_text());assert side['sha256']==sha(file) and side['plan_sha256']==sha(A/'plan.json');records=side['records']
  with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
 for ti in needed:
  if arrays['valid'][ti]:continue
  source=info['source_choices'][ti];n=info['matched_prompt_sizes'][ti];row={'family':family,'target':str(targets[ti]),'source_family':source,'pairs':[]};start=time.time()
  if source:
   assert source!=family;seed=20261020+plan['family_order'].index(family)*10000+ti
   for label,index in [('real',ti),('shuffled',info['shuffled_target_indices'][ti])]:
    target=str(targets[index]);sidepath,side=prompt(source,target,index,plan,genes,vocab);inputs={}
    for kind in ['control','treated']:
     f=sidepath.with_name(sidepath.stem+'_'+kind+'.npz');assert sha(f)==side['records'][kind]['sha256'];inputs[kind]=sparse.load_npz(f)[:n].tocsr()
    mu0,theta0=predict_means(model,inputs['control'],query,seed=seed,batch_size=4);mut,thetat=predict_means(model,inputs['treated'],query,seed=seed,batch_size=4);h0=hashlib.sha256(mu0.tobytes()).hexdigest();ht=hashlib.sha256(mut.tobytes()).hexdigest()
    if stage=='development':
     prior=next(r for r in original_records[str(targets[ti])]['pairs'] if r['kind']==label);assert h0==prior['control_mean_array_sha256'] and ht==prior['treated_mean_array_sha256'],'Foundation mean replay changed'
    supported=np.asarray(side['vocabulary_measured_mask'],bool)&np.isin(vocab,genes[q['measured']]);summary=summarize_pair(mu0,mut,q['query_mean'],q,q['measured'],supported,vocab,genes,theta0,thetat);summary['shrink']*=n/(n+100);assert summary['program_coverage']>=.8
    if label=='real':
     if original is not None:
      for k in ['lfc','shrink','raw','supported']:np.testing.assert_array_equal(summary[k],original[k][ti])
     else:
      for k in ['lfc','shrink','raw','supported']:arrays[k][ti]=summary[k]
     arrays['program_delta'][ti]=summary['program_delta']*n/(n+100);arrays['direction'][ti]=summary['direction']*n/(n+100)
    else:
     if original is not None:np.testing.assert_array_equal(summary['lfc'],original['shuffled_lfc'][ti]);np.testing.assert_array_equal(summary['shrink'],original['shuffled_shrink'][ti])
     else:arrays['shuffled_lfc'][ti]=summary['lfc'];arrays['shuffled_shrink'][ti]=summary['shrink']
     arrays['shuffled_program_delta'][ti]=summary['program_delta']*n/(n+100)
    row['pairs'].append({'kind':label,'source_target':target,'source_input_path':str(sidepath.relative_to(ROOT)).replace('\\','/'),'source_input_sha256':sha(sidepath),'control_mean_array_sha256':h0,'treated_mean_array_sha256':ht,'program_coverage':summary['program_coverage'],'parent_mean_hash_matched':stage=='development'});del mu0,mut,theta0,thetat
   real=float(np.linalg.norm(arrays['program_delta'][ti]));other=float(np.linalg.norm(arrays['shuffled_program_delta'][ti]));arrays['shuffled_program_delta'][ti]*=real/max(other,1e-12);row.update(real_program_norm=real,shuffled_program_norm_before_matching=other)
  arrays['valid'][ti]=True;records.append(row);np.savez_compressed(file,**arrays);write(file.with_suffix('.json'),{'sha256':sha(file),'plan_sha256':sha(A/'plan.json'),'query_input_sha256':sha(qfile),'parent_expert_sha256':sha(parent_file) if parent_file.exists() else None,'records':records,'mean_response_arrays_unchanged':stage=='development'});print('V21 active Stack programs',family,targets[ti],round(time.time()-start,2),'seconds',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);p.add_argument('--family',required=True);args=p.parse_args();main(args.stage,args.family)
