"""Strict aligned Stack GPU smoke: paired source NTC/treatment prompts and repeat."""
from pathlib import Path
import os,sys,json,time,hashlib,warnings
ROOT=Path(__file__).resolve().parent;R=ROOT/'reference/v20_stack';A=ROOT/'experiments/v20_stack_transfer/analysis'
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8');sys.path.insert(0,str(R/'code/src'))
import numpy as np,pandas as pd,torch,anndata as ad
from scipy import sparse
from stack.models.finetune import ICL_FinetunedModel
warnings.filterwarnings('ignore',message='Observation names are not unique')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for v in iter(lambda:f.read(8*1024**2),b''):h.update(v)
 return h.hexdigest()
def main(steps=1,batch_size=4):
 pd.options.future.infer_string=False
 torch.set_num_threads(4);torch.use_deterministic_algorithms(True,warn_only=True);torch.manual_seed(20020)
 assert sha(R/'bc_large_aligned.ckpt')=='f93cf6f42f36c8a85dc570d92e801c1fc3e1f45d55741e6bb250d178a6b6ad36';ck=torch.load(R/'bc_large_aligned.ckpt',map_location='cpu',weights_only=False);cfg=ck['hyper_parameters']['model_config'];assert cfg['n_genes']==15012 and cfg['n_cells']==512;model=ICL_FinetunedModel(**cfg);state={k[6:] if k.startswith('model.') else k:v for k,v in ck['state_dict'].items()};model.load_state_dict(state,strict=True);del ck,state;model.to('cuda:0').eval()
 folder=A/'smoke_inputs';rec=json.loads((folder/'inputs.json').read_text());vocab=rec['vocabulary'];inputs={}
 for name,info in rec['records'].items():
  p=folder/(name+'.npz');assert sha(p)==info['sha256'];x=sparse.load_npz(p).astype(np.float32);assert (x.data>=0).all() and np.all(x.data==np.floor(x.data));inputs[name]=ad.AnnData(x,var=pd.DataFrame(index=vocab),obs=pd.DataFrame(index=[name+'_'+str(i) for i in range(x.shape[0])]))
 outputs={};t=time.time()
 for name,prompt in [('control','source_control'),('treated','source_treated'),('repeat_control','source_control')]:
  print('Stack paired smoke',name,flush=True)
  generated=model.get_incontext_generation(inputs[prompt],inputs['query_control'],str(R/'genes.pkl'),num_steps=steps,prompt_ratio=.25,context_ratio=.4,context_ratio_min=.2,mode='mdm',batch_size=batch_size,num_workers=0,show_progress=False,random_seed=20020,filter_organism=False)
  count=generated[0] if isinstance(generated,tuple) else generated;count=sparse.csr_matrix(count);assert count.shape==(400,15012) and np.isfinite(count.data).all() and (count.data>=0).all() and np.all(count.data==np.floor(count.data));outputs[name]=count;sparse.save_npz(A/(f'smoke_{steps}step_b{batch_size}_'+name+'.npz'),count)
 assert (outputs['control']!=outputs['repeat_control']).nnz==0
 result={'passed':True,'checkpoint_sha256':sha(R/'bc_large_aligned.ckpt'),'strict_state_dict':True,'dtype':'FP32','deterministic_algorithms':'warn_only for author CUDA median values; actual full-count repeat required','model_config':cfg,'batch_size':batch_size,'num_steps':steps,'mode':'author mdm','counts_shape':[400,15012],'bitwise_control_repeat':True,'peak_gpu_bytes':torch.cuda.max_memory_allocated(),'seconds':time.time()-t,'query_treatments_used':False,'supported_output_genes':int(sum(rec['supported_mask'])),'counts_sha256':{k:sha(A/(f'smoke_{steps}step_b{batch_size}_'+k+'.npz')) for k in outputs},'code_sha256':sha(Path(__file__))}
 (A/f'smoke_{steps}step_b{batch_size}.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--steps',type=int,default=1);p.add_argument('--batch-size',type=int,default=4);args=p.parse_args();main(args.steps,args.batch_size)
