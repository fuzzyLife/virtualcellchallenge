"""Pinned Stack one-step conditional means for source-paired, controls-only queries."""
from pathlib import Path
import os,sys,json,argparse,time,hashlib
ROOT=Path(__file__).resolve().parent;R=ROOT/'reference/v20_stack';A=ROOT/'experiments/v20_stack_transfer/analysis';os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8');sys.path.insert(0,str(R/'code/src'))
import numpy as np,torch
from scipy import sparse
from stack.models.finetune import ICL_FinetunedModel
from stack_v20_mean_adapter import predict_means
from v20_methods import summarize_pair

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for v in iter(lambda:f.read(8*1024**2),b''):h.update(v)
 return h.hexdigest()
def write(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
def main(stage,family):
 torch.set_num_threads(4);torch.manual_seed(20020);torch.use_deterministic_algorithms(True);plan=json.loads((A/'plan.json').read_text())
 for name,h in plan['dependencies'].items():assert sha(ROOT/name)==h
 for item in plan['input_files']:assert sha(ROOT/item['path'])==item['sha256']
 assert sha(A/'assets.json')==plan['assets_manifest_sha256']
 for item in json.loads((A/'assets.json').read_text())['files']:assert sha(ROOT/item['path'])==item['sha256']
 if stage=='confirmation':assert json.loads((A/'selection.json').read_text())['new_candidate_passes']
 assert sha(R/'bc_large_aligned.ckpt')==plan['checkpoint_sha256'];ck=torch.load(R/'bc_large_aligned.ckpt',map_location='cpu',weights_only=False);model=ICL_FinetunedModel(**ck['hyper_parameters']['model_config']);assert all(k.startswith('model.') for k in ck['state_dict']);model.load_state_dict({k[6:]:v for k,v in ck['state_dict'].items()},strict=True);del ck;model.to('cuda:0').eval()
 folder=A/'inputs';rec=json.loads((folder/f'{family}.json').read_text());assert sha(folder/f'{family}_query.npz')==rec['query_arrays_sha256'];assert sha(folder/f'{family}_counts.npz')==rec['query_counts_sha256']
 with np.load(folder/f'{family}_query.npz',allow_pickle=False) as z:q={k:z[k].copy() for k in z.files}
 query=sparse.load_npz(folder/f'{family}_counts.npz');vocab=json.loads((R/'genes.json').read_text());targets=q['targets'];genes=q['genes'];context=next(c for c in plan['contexts'] if c['family']==family);needed=sorted({int(np.flatnonzero(targets==r['target'])[0]) for r in context['cells'] if r['split']==stage});file=A/'experts'/f'{family}_{stage}.npz';file.parent.mkdir(exist_ok=True);g=len(genes)
 arrays={k:np.zeros((300,g),np.float32) for k in ['lfc','shrink','raw','shuffled_lfc','shuffled_shrink']};arrays.update(supported=np.zeros((300,g),bool),program_delta=np.zeros((300,400,16),np.float32),direction=np.zeros((300,16),np.float32),valid=np.zeros(300,bool),genes=genes,targets=targets)
 records=[]
 if file.exists():
  side=json.loads(file.with_suffix('.json').read_text());assert side['sha256']==sha(file) and side['plan_sha256']==sha(A/'plan.json');records=side['records']
  with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
 for ti in needed:
  if arrays['valid'][ti]:continue
  source=rec['source_choices'][ti];start=time.time();row={'family':family,'target':str(targets[ti]),'source_family':source,'input_count':rec['matched_prompt_sizes'][ti],'pairs':[]}
  if source:
   assert source!=family;n=rec['matched_prompt_sizes'][ti];seed=20261020+plan['family_order'].index(family)*10000+ti
   for label,index in [('real',ti),('shuffled',rec['shuffled_target_indices'][ti])]:
    target=str(targets[index]);sidepath=folder/'sources'/f'{source}_{target}.json';side=json.loads(sidepath.read_text());prompts={}
    for kind in ['control','treated']:
     p=folder/'sources'/f'{source}_{target}_{kind}.npz';assert sha(p)==side['records'][kind]['sha256'];prompts[kind]=sparse.load_npz(p)[:n].tocsr();assert prompts[kind].shape[0]==n
    supported=np.asarray(side['vocabulary_measured_mask'],bool)&np.isin(vocab,genes[q['measured']]);mu0,theta0=predict_means(model,prompts['control'],query,seed=seed,batch_size=4);mut,thetat=predict_means(model,prompts['treated'],query,seed=seed,batch_size=4);assert np.isfinite(mu0).all() and np.isfinite(mut).all();summary=summarize_pair(mu0,mut,q['query_mean'],q,q['measured'],supported,vocab,genes,theta0,thetat);summary['shrink']*=n/(n+100)
    if label=='real':
     for k in ['lfc','shrink','raw','supported','program_delta','direction']:arrays[k][ti]=summary[k]
     arrays['program_delta'][ti]*=n/(n+100);arrays['direction'][ti]*=n/(n+100)
    else:arrays['shuffled_lfc'][ti]=summary['lfc'];arrays['shuffled_shrink'][ti]=summary['shrink']
    row['pairs'].append({'kind':label,'source_target':target,'source_input_sha256':sha(sidepath),'seed':seed,'control_mean_array_sha256':hashlib.sha256(mu0.tobytes()).hexdigest(),'treated_mean_array_sha256':hashlib.sha256(mut.tobytes()).hexdigest(),'supported_output_genes':int(summary['supported'].sum()),'program_coverage':summary['program_coverage']});del mu0,mut,theta0,thetat
  arrays['valid'][ti]=True;row['seconds']=time.time()-start;records.append(row);np.savez_compressed(file,**arrays);write(file.with_suffix('.json'),{'sha256':sha(file),'plan_sha256':sha(A/'plan.json'),'query_input_sha256':rec['query_arrays_sha256'],'records':records,'strict_state_dict':True,'precision':'FP32','decoder':'one-step original conditional mean/dispersion, bitwise parity certified','source_prompt_holdout':True,'pretraining_overlap_audited':False});print('Stack paired means',family,targets[ti],round(row['seconds'],2),'seconds',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);p.add_argument('--family',required=True);a=p.parse_args();main(a.stage,a.family)
