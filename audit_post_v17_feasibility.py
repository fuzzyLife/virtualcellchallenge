"""Audit candidate assets and hardware without loading model weights or pickle code."""
from pathlib import Path
import json,re,pickletools,platform,subprocess
import numpy as np,pandas as pd
import build_submission as b
ROOT=b.ROOT;A=ROOT/'analysis/post_v17_research'

def decode_unicode_list(payload):
 """Interpret only literal lists and NumPy unicode scalar descriptions; never call globals."""
 stack=[];memo={};mark=object();globals_seen=set()
 for op,arg,pos in pickletools.genops(payload):
  name=op.name
  if name in ['PROTO','FRAME']:continue
  if name=='MARK':stack.append(mark)
  elif name=='EMPTY_LIST':stack.append([])
  elif name in ['SHORT_BINUNICODE','BINUNICODE','BINUNICODE8','SHORT_BINBYTES','BINBYTES','BINBYTES8','BININT','BININT1','BININT2']:stack.append(arg)
  elif name in ['NEWTRUE','NEWFALSE','NONE']:stack.append({'NEWTRUE':True,'NEWFALSE':False,'NONE':None}[name])
  elif name=='MEMOIZE':memo[len(memo)]=stack[-1]
  elif name in ['BINGET','LONG_BINGET']:stack.append(memo[arg])
  elif name=='STACK_GLOBAL':
   function=stack.pop();module=stack.pop();ref=(module,function)
   if ref not in [('numpy','dtype'),('numpy._core.multiarray','scalar'),('numpy.core.multiarray','scalar')]:raise ValueError('Unsupported pickle global')
   globals_seen.add(ref);stack.append(('reference',ref))
  elif name in ['TUPLE','TUPLE1','TUPLE2','TUPLE3']:
   if name=='TUPLE':
    i=max(i for i,v in enumerate(stack) if v is mark);values=tuple(stack[i+1:]);del stack[i:];stack.append(values)
   else:
    n=int(name[-1]);values=tuple(stack[-n:]);del stack[-n:];stack.append(values)
  elif name=='REDUCE':
   args=stack.pop();ref=stack.pop()
   if ref==('reference',('numpy','dtype')):
    if not re.fullmatch(r'U[0-9]+',args[0]):raise ValueError('Only unicode dtypes supported')
    stack.append({'dtype':args[0],'endian':'<'})
   elif ref in [('reference',('numpy._core.multiarray','scalar')),('reference',('numpy.core.multiarray','scalar'))]:
    dtype,raw=args
    if not isinstance(dtype,dict) or not isinstance(raw,bytes):raise ValueError('Invalid literal scalar')
    value=raw.decode('utf-32-le' if dtype['endian']=='<' else 'utf-32-be');stack.append(value)
   else:raise ValueError('No function execution supported')
  elif name=='BUILD':
   state=stack.pop();dtype=stack[-1]
   if not isinstance(dtype,dict) or state[1] not in ['<','>']:raise ValueError('Unsupported dtype state')
   dtype['endian']=state[1]
  elif name=='APPENDS':
   i=max(i for i,v in enumerate(stack) if v is mark);values=stack[i+1:];del stack[i:]
   if not isinstance(stack[-1],list):raise ValueError('Expected a literal list')
   stack[-1].extend(values)
  elif name=='STOP':
   if len(stack)!=1 or not isinstance(stack[0],list) or not all(isinstance(v,str) for v in stack[0]):raise ValueError('Invalid vocabulary structure')
   return stack[0]
  else:raise ValueError('Unsupported opcode '+name)
 raise ValueError('Incomplete vocabulary')

def main():
 manifest,genes,targets=b.read_spec(ROOT);genes=np.asarray(genes,str);targets=np.asarray(targets,str);coverage=pd.DataFrame({'gene_name':genes,'panel_target':np.isin(genes,targets)});rows=[]
 values=decode_unicode_list((A/'stack_gene_vocabulary.pkl').read_bytes());assert len(set(values))==len(values)
 coverage['stack_vocabulary']=np.isin(genes,values);rows.append({'asset':'Stack published vocabulary','vocabulary_size':len(values),'output_genes':int(coverage.stack_vocabulary.sum()),'targets':int(np.isin(targets,values).sum()),'file_sha256':b.sha256(A/'stack_gene_vocabulary.pkl')})
 for name in ['gene_text','perturbation_text','esm2']:
  d=json.loads((A/(name+'_meta.json')).read_text());coverage['PIE_'+name]=np.isin(genes,d['keys']);rows.append({'asset':'PIE '+name,'vocabulary_size':len(d['keys']),'output_genes':int(np.isin(genes,d['keys']).sum()),'targets':int(np.isin(targets,d['keys']).sum()),'file_sha256':b.sha256(A/(name+'_meta.json'))})
 coverage.to_csv(A/'feature_coverage.csv',index=False);pd.DataFrame(rows).to_csv(A/'feature_coverage_summary.csv',index=False)
 runtime={'python':platform.python_version(),'platform':platform.platform()}
 try:
  import torch
  runtime.update(torch=torch.__version__,torch_cuda=torch.version.cuda,cuda_available=torch.cuda.is_available(),gpu=torch.cuda.get_device_name(0),gpu_memory_bytes=torch.cuda.get_device_properties(0).total_memory,bfloat16_native=torch.cuda.is_bf16_supported(including_emulation=False))
 except Exception as error:runtime['torch_inspection_error']=str(error)
 assets={}
 for name in ['PIE_replogle_xdataset','PIE_replogle_wdataset','Stack-Large','Stack-Large-Aligned']:
  info=json.loads((A/(name+'_model_info.json')).read_text());tree=json.loads((A/(name+'_model_files.json')).read_text());assets[name]={'revision':info['sha'],'gated':info.get('gated'),'checkpoint_files':{v['path']:v['size'] for v in tree if v['type']=='file' and v['path'].endswith('.ckpt')}}
 training_meta=json.loads((A/'PIE_training_meta.json').read_text());training={'revision':'20c9faef76fc96fdc809871340bd93499413dfe7','perturbation_names':len(training_meta['pert_to_id']),'competition_target_overlap':sorted(set(targets)&set(training_meta['pert_to_id'])),'training_output_genes':len(training_meta['genes']),'competition_output_overlap':len(set(genes)&set(training_meta['genes'])),'metadata_sha256':b.sha256(A/'PIE_training_meta.json')}
 files=json.loads((A/'PIE_sources_files.json').read_text());feature_sizes={v['path']:v.get('size') for v in files if v['type']=='file' and v['path'].endswith('embeddings.npy')}
 b.write_json(A/'feasibility_audit.json',{'date':'2026-10-13','script_sha256':b.sha256(Path(__file__)),'coverage':rows,'PIE_exact_training_revision':training,'runtime':runtime,'assets':assets,'PIE_embedding_file_bytes':feature_sizes,'models_loaded':False,'checkpoints_downloaded':False,'new_predictions_generated':False,'vocabulary_decoding':'Non-executing pickletools interpreter restricted to literal strings/lists and numpy-unicode scalar descriptions','limitations':['Feature availability does not prove predictive coverage or checkpoint inference compatibility.','PIE within-dataset training output axis has 6642 genes; extending to new genes is an extrapolation test.','PIE cross-dataset checkpoint trained on ARC VCC25 and Orion: existing H1/HCT116/HEK293T holdouts cannot establish its unseen-family performance.','Stack pretraining overlap with public benchmark cell IDs is not yet audited.']});print(json.dumps({'coverage':rows,'PIE_exact_training_revision':training,'runtime':runtime},indent=2),flush=True)
if __name__=='__main__':main()
