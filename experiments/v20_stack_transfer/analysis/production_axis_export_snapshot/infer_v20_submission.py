"""Strict FP32 paired-control Stack means for all 900 competition groups."""
from pathlib import Path
import os,sys,json,time,hashlib
ROOT=Path(__file__).resolve().parent;R=ROOT/'reference/v20_stack';A=ROOT/'experiments/v20_stack_transfer/analysis';os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8');sys.path.insert(0,str(R/'code/src'))
import numpy as np,torch
from scipy import sparse
from stack.models.finetune import ICL_FinetunedModel
from stack_v20_mean_adapter import predict_means
from v20_methods import summarize_pair
from v20_production import verify_production
import build_submission as b

def main():
 plan=verify_production();torch.set_num_threads(4);torch.manual_seed(20020);torch.use_deterministic_algorithms(True);assets=json.loads((A/'assets.json').read_text())
 for row in assets['files']:assert b.sha256(ROOT/row['path'])==row['sha256']
 assert b.sha256(R/'bc_large_aligned.ckpt')==assets['checkpoint_sha256'];checkpoint=torch.load(R/'bc_large_aligned.ckpt',map_location='cpu',weights_only=False);model=ICL_FinetunedModel(**checkpoint['hyper_parameters']['model_config']);assert all(k.startswith('model.') for k in checkpoint['state_dict']);model.load_state_dict({k[6:]:v for k,v in checkpoint['state_dict'].items()},strict=True);del checkpoint;model.to('cuda:0').eval();folder=A/'production_inputs';out=A/'production_experts';out.mkdir(exist_ok=True);vocab=json.loads((R/'genes.json').read_text());manifest,genes,targets=b.read_spec(ROOT)
 for ci,context in enumerate(manifest['contexts']):
  rec=json.loads((folder/(context+'.json')).read_text());qfile=folder/(context+'_query.npz');matrix=folder/(context+'_counts.npz');assert b.sha256(qfile)==rec['query_arrays_sha256'] and b.sha256(matrix)==rec['query_counts_sha256'];query=sparse.load_npz(matrix)
  with np.load(qfile,allow_pickle=False) as z:q={k:z[k].copy() for k in z.files}
  measured_names=set(genes[q['measured']]);query_vocab_mask=np.array([g in measured_names for g in vocab]);file=out/(context+'.npz');arrays={'lfc':np.zeros((300,len(genes)),np.float32),'shrink':np.zeros((300,len(genes)),np.float32),'supported':np.zeros((300,len(genes)),bool),'valid':np.zeros(300,bool),'genes':genes,'targets':targets};records=[]
  if file.exists():
   side=json.loads(file.with_suffix('.json').read_text());assert side['sha256']==b.sha256(file) and side['production_plan_sha256']==b.sha256(A/'production_plan.json');records=side['records']
   with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
  for ti,target in enumerate(targets):
   if arrays['valid'][ti]:continue
   start=time.time();source=rec['source_choices'][ti];n=rec['matched_prompt_sizes'][ti];seed=20261020+ci*10000+ti;record={'context':context,'target':str(target),'source_family':source,'input_count':n,'seed':seed}
   if source:
    sidepath=folder/'sources'/(source+'_'+str(target)+'.json');side=json.loads(sidepath.read_text());prompts={}
    for kind in ['control','treated']:
     p=folder/'sources'/(source+'_'+str(target)+'_'+kind+'.npz');assert b.sha256(p)==side['records'][kind]['sha256'];prompts[kind]=sparse.load_npz(p)[:n].tocsr();assert prompts[kind].shape[0]==n
    supported=np.asarray(side['vocabulary_measured_mask'],bool)&query_vocab_mask;mu0,theta0=predict_means(model,prompts['control'],query,seed=seed,batch_size=4);mut,thetat=predict_means(model,prompts['treated'],query,seed=seed,batch_size=4);summary=summarize_pair(mu0,mut,q['query_mean'],q,q['measured'],supported,vocab,genes,theta0,thetat)
    for key in ['lfc','shrink','supported']:arrays[key][ti]=summary[key]
    arrays['shrink'][ti]*=n/(n+100);assert np.isfinite(arrays['lfc'][ti]).all() and np.isfinite(arrays['shrink'][ti]).all();record.update(source_input_sha256=b.sha256(sidepath),control_mean_array_sha256=hashlib.sha256(mu0.tobytes()).hexdigest(),treated_mean_array_sha256=hashlib.sha256(mut.tobytes()).hexdigest(),supported_output_genes=int(summary['supported'].sum()));del mu0,mut,theta0,thetat
   arrays['valid'][ti]=True;record['seconds']=time.time()-start;records.append(record)
   if (ti+1)%20==0 or ti==299:
    partial=file.with_suffix('.partial.npz');np.savez_compressed(partial,**arrays);partial.replace(file);b.write_json(file.with_suffix('.json'),{'sha256':b.sha256(file),'production_plan_sha256':b.sha256(A/'production_plan.json'),'query_input_sha256':b.sha256(qfile),'records':records,'strict_state_dict':True,'precision':'FP32','decoder':'Original-author one-step conditional means and dispersion, bitwise adapter parity certified','pretraining_overlap_audited':False});print('Production Stack means',context,ti+1,'/ 300',flush=True)
 print('All production mean experts complete',flush=True)
if __name__=='__main__':main()
