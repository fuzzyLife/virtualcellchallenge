"""Strict cross-dataset PIE FP32 inference with three retained heads and donor evidence."""
import argparse,sys,time
import numpy as np
from v57_pie import *
sys.path.insert(0,str(R/'code/src'))
import torch
from pie.model.pie import PieModel,ModelConfig
from pie.model.heads import readout
from pie.data.datamodule import DataStats
from pie.data.dataset import SourceLookup,collate,load_aliases
from pie.sources.contract import read_source
from build_v57_pie_evidence import blocks

def main(source):
 p=plan();torch.set_num_threads(4);torch.manual_seed(SEED);torch.use_deterministic_algorithms(True);torch.backends.cudnn.benchmark=False
 ckpath=CACHE/'best_auprc.ckpt';assert sha(ckpath)=='e60a1304d074706405939bde94c5edf9851e78070f2d584c6b5025a002ed9a7c'
 ck=torch.load(ckpath,map_location='cpu',weights_only=False)
 cfg=ModelConfig.model_validate(ck['pie']['config']['model']).model_copy(update={'inference_chunk_size':512});stats=DataStats.model_validate(ck['pie']['data_stats'])
 assert stats.evidence_key==read(CACHE/'evidence_projection/projection.json')['original_evidence_key'] and stats.train_json_sha256==p['config']['donor_split_sha256']
 paths={n:((R/'sources'/('ncbi_subset' if n=='ncbi_text' else n)) if n not in ['smiles','l1000_tas','prism_secondary','jump_morphology'] else CACHE/'sources'/n) for n in stats.source_dims}
 sources={n:read_source(v) for n,v in paths.items()};aliases=load_aliases(R/'code/src/pie/sources/curated_aliases.yaml');lookup=SourceLookup(sources,aliases);assert lookup.dims==stats.source_dims
 inp=CACHE/'inputs'/(source+'.npz');rec=read(A/'inputs'/(source+'.json'));assert sha(inp)==rec['sha256']
 with np.load(inp,allow_pickle=False) as z:data={k:z[k].copy() for k in z.files}
 positions=np.flatnonzero(data['support']);genes=data['genes'][positions].tolist();targets=data['targets'].tolist();gt=read_source(R/'sources/gene_text')
 text=torch.from_numpy(np.array(gt.embeddings[[gt.key_to_row[g] for g in genes]],dtype=np.float32));model=PieModel(cfg,stats,stats.source_dims,text)
 model.load_state_dict({k[6:]:v for k,v in ck['state_dict'].items()},strict=True);del ck;text=None;model.cuda().eval()
 folder=CACHE/'predictions';folder.mkdir(exist_ok=True);path=folder/(source+'.npz');n=len(targets);ng=len(genes)
 output={k:np.full((n,ng),np.nan,np.float32) for k in ['lfc','p_de','delta_p']};done=np.zeros(n,bool);records={}
 if path.exists():
  with np.load(path,allow_pickle=False) as z:
   done=z['done'].copy()
   for k in output:output[k]=z[k].copy()
  records=read(A/'predictions'/(source+'.json'))['records']
 needed={targets.index(t) for t in p['cohorts'][source]['targets']};needed|={p['shuffle'][i] for i in needed}
 smoke=[]
 def predict(ti,disable=False):
  target=targets[ti];ev,diag=blocks(target,source,genes)
  if disable:ev={k:np.zeros_like(v) for k,v in ev.items()}
  tokens={}
  for name in lookup.dims:
   v=lookup.tokens(name,context='VCC_UNKNOWN_'+source,perturbation=target)
   if v is not None:tokens[name]=torch.from_numpy(v.astype(np.float32))
  own=genes.index(target) if target in genes else -1
  sample={'dir_index':0,'row':-1,'context':'VCC_UNKNOWN_'+source,'perturbation':target,'source_tokens':tokens,'gene_ids':torch.arange(ng),'ctrl_mean':torch.from_numpy(data['ctrl_mean'][positions]),'evidence':{k:torch.from_numpy(v) for k,v in ev.items()},'target_gene':own,'target_gene_idx':own}
  batch=collate([sample],stats.source_dims);batch.source_tokens={k:v.cuda() for k,v in batch.source_tokens.items()};batch.source_masks={k:v.cuda() for k,v in batch.source_masks.items()};batch.dataset_ids=batch.dataset_ids.cuda()
  with torch.inference_mode():v=readout(model(batch)[0],cfg,stats.delta_p)
  values={'lfc':v.lfc.cpu().numpy()[0],'p_de':v.p_de.cpu().numpy()[0],'delta_p':v.delta_p.cpu().numpy()[0]}
  assert all(np.isfinite(v).all() for v in values.values())
  return values,{**diag,'source_tokens':{k:int(v.sum().cpu()) for k,v in batch.source_masks.items()},'peak_GPU_bytes':int(torch.cuda.max_memory_allocated())}
 for ti in sorted(needed):
  if done[ti]:continue
  begin=time.time();values,diag=predict(ti)
  if not records:
   repeated,_=predict(ti)
   for k in values:np.testing.assert_array_equal(values[k],repeated[k])
  for k,v in values.items():output[k][ti]=v
  done[ti]=True;records[str(ti)]={'target':targets[ti],**diag,'seconds':time.time()-begin}
  # Atomic resumable chunks avoid quadratic compressed rewrites per prediction.
  if done.sum()%10==0 or ti==max(needed):
   tmp=path.with_suffix('.partial.npz');np.savez_compressed(tmp,**output,done=done,positions=positions,genes=np.array(genes),targets=data['targets']);tmp.replace(path)
   write(A/'predictions'/(source+'.json'),{'sha256':sha(path),'plan_sha256':sha(A/'plan.json'),'input_sha256':sha(inp),'evidence_receipt_sha256':sha(A/'evidence_receipt.json'),'records':records,'required_targets':sorted(needed),'strict_state_dict':True,'all_three_heads_retained':True})
  print('PIE cross',source,targets[ti],int(done.sum()),'/',len(needed),round(time.time()-begin,2),'sec',flush=True)
 # No-evidence ablation on fixed metadata-chosen count-pilot targets only.
 for target in p['cohorts'][source]['count_pilot']:
  ti=targets.index(target);v,diag=predict(ti,True);smoke.append({'target':target,'seen':target in p['seen_targets'],'mean_abs_head_change':{k:float(abs(v[k]-output[k][ti]).mean()) for k in v},'donor_genes_with_evidence':records[str(ti)]['donor_genes']})
 write(A/'predictions'/(source+'_evidence_ablation.json'),{'rows':smoke,'scope':'Fixed pilot only; numerical conditioning ablation, not an independent biological performance result'})
 print('Inference complete',source,flush=True)

if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('--source',choices=SOURCES,required=True);main(q.parse_args().source)
