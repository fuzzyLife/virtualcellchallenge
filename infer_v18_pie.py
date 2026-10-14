"""Strict pretrained PIE forward adapter for novel datasets and targets only (Python 3.12)."""
from pathlib import Path
import argparse,hashlib,json,os,sys,time
ROOT=Path(__file__).resolve().parent;R=ROOT/'reference/v18_pie';A=ROOT/'experiments/v18_pie_expert/analysis'
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8');sys.path.insert(0,str(R/'code/src'))
import numpy as np,torch
from pie.model.pie import PieModel,ModelConfig
from pie.model.heads import readout
from pie.data.datamodule import DataStats
from pie.data.dataset import SourceLookup,collate
from pie.sources.contract import read_source
from pie.data.evidence import EVIDENCE_KEYS

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for v in iter(lambda:f.read(8*1024**2),b''):h.update(v)
 return h.hexdigest()
def write(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
def main(mode,family,device):
 torch.set_num_threads(4);torch.manual_seed(42);np.random.seed(42);torch.use_deterministic_algorithms(True);torch.backends.cudnn.benchmark=False
 plan=json.loads((A/'plan.json').read_text()) if mode!='smoke' else None
 assets=json.loads((A/'assets.json').read_text());assert assets['model_revision']=='055c7a2cabdaff5659121c64fa2b609c231fc5be';assert sha(R/'best_auprc.ckpt')=='b9d5a0b3943666d892a56e7a9368f7850826b3f1fe6879e63a1a590f8add7d84'
 # Hash-verified public author checkpoint. Preserve all trained modules; no fitted weights or heads replaced.
 ck=torch.load(R/'best_auprc.ckpt',map_location='cpu',weights_only=False);assert ck['pie']['format_version']==1
 cfg=ModelConfig.model_validate(ck['pie']['config']['model']);cfg=cfg.model_copy(update={'inference_chunk_size':512});stats=DataStats.model_validate(ck['pie']['data_stats']);assert stats.evidence_datasets==['replogle']
 split=R/'train.json';assert sha(split)==stats.train_json_sha256
 train=json.loads((ROOT/'analysis/post_v17_research/PIE_training_meta.json').read_text());source={name:read_source(R/'sources'/('ncbi_subset' if name=='ncbi_text' else name)) for name in stats.source_dims};lookup=SourceLookup(source,{})
 assert lookup.dims==stats.source_dims
 gt=read_source(R/'sources/gene_text')
 with np.load(A/'inputs'/f'{family}.npz',allow_pickle=False) as z:inputs={k:z[k].copy() for k in z.files}
 genes=inputs['genes'].tolist();targets=inputs['targets'].tolist();assert not set(targets)&set(train['pert_to_id']);assert len(genes)==6290 and set(genes)<=set(train['genes']);assert np.isfinite(inputs['ctrl_mean']).all()
 text=torch.from_numpy(np.array(gt.embeddings[[gt.key_to_row[g] for g in genes]],dtype=np.float32));model=PieModel(cfg,stats,stats.source_dims,text)
 state=ck['state_dict'];assert all(k.startswith('model.') for k in state);model.load_state_dict({k[6:]:v for k,v in state.items()},strict=True);del ck,state;text=None;model.to(device).eval()
 def batch(target):
  tokens={}
  for name in lookup.dims:
   v=lookup.tokens(name,context='VCC_UNKNOWN_'+family,perturbation=target)
   if v is not None:tokens[name]=torch.from_numpy(v.astype(np.float32))
  # New dataset not present in training and no target donor records -> author-derived all-zero evidence.
  ev={k:torch.zeros(len(genes),3 if k=='evidence_prov' else 4) for k in EVIDENCE_KEYS}
  own=genes.index(target) if target in genes else -1
  sample={'dir_index':0,'row':-1,'context':'VCC_UNKNOWN_'+family,'perturbation':target,'source_tokens':tokens,'gene_ids':torch.arange(len(genes)),'ctrl_mean':torch.from_numpy(inputs['ctrl_mean']),'evidence':ev,'target_gene':own,'target_gene_idx':own}
  b=collate([sample],stats.source_dims)
  b.source_tokens={k:v.to(device) for k,v in b.source_tokens.items()};b.source_masks={k:v.to(device) for k,v in b.source_masks.items()};b.dataset_ids=b.dataset_ids.to(device)
  return b
 def one(t):
  with torch.inference_mode():v=readout(model(batch(t))[0],cfg,stats.delta_p)
  out=[v.lfc.cpu().numpy()[0],v.p_de.cpu().numpy()[0],v.delta_p.cpu().numpy()[0]]
  assert all(np.isfinite(x).all() for x in out);assert ((out[1]>=0)&(out[1]<=1)).all();return out
 if mode=='smoke':
  target='VIM' if 'VIM' in targets else targets[0];t=time.time();a=one(target);z=one(target)
  for x,y in zip(a,z):np.testing.assert_array_equal(x,y)
  record={'passed':True,'family':family,'target':target,'device':device,'torch':torch.__version__,'dtype':'float32','batch_size':1,'decoder_chunk_size':512,'gene_count':len(genes),'strict_state_dict':True,'bitwise_repeat':True,'seconds_two_forwards':time.time()-t,'peak_gpu_bytes':torch.cuda.max_memory_allocated() if device.startswith('cuda') else None,'original_train_split_sha256':stats.train_json_sha256,'original_evidence_key':stats.evidence_key,'source_mask':{k:int(v.sum()) for k,v in batch(target).source_masks.items()},'lfc_range':[float(a[0].min()),float(a[0].max())],'p_de_range':[float(a[1].min()),float(a[1].max())]}
  write(A/'smoke.json',record);print(json.dumps(record,indent=2),flush=True);return
 assert plan['dependencies']['infer_v18_pie.py']==sha(Path(__file__)) and plan['assets_sha256']==sha(A/'assets.json')
 if mode=='development':
  needed={r['target'] for c in plan['contexts'] if c['family']==family for r in c['cells'] if r['split']=='development'};indices=[targets.index(t) for t in needed];indices+= [plan['shuffle'][i] for i in indices];indices=sorted(set(indices))
 elif mode=='confirmation':
  selection=json.loads((A/'selection.json').read_text());assert selection['new_candidate_passes'];needed={r['target'] for c in plan['contexts'] if c['family']==family for r in c['cells'] if r['split']=='confirmation'};indices=[targets.index(t) for t in needed];indices+= [plan['shuffle'][i] for i in indices];indices=sorted(set(indices))
 else:indices=list(range(len(targets)))
 path=A/'inference'/f'{family}_{mode}.npz';path.parent.mkdir(exist_ok=True);lfc=np.zeros((300,len(genes)),np.float32);prob=np.zeros_like(lfc);valid=np.zeros(300,bool)
 if path.exists():
  rec=json.loads(path.with_suffix('.json').read_text());assert rec['sha256']==sha(path) and rec['plan_sha256']==sha(A/'plan.json')
  with np.load(path,allow_pickle=False) as z:lfc=z['lfc'];prob=z['p_de'];valid=z['valid']
 for i in indices:
  if valid[i]:continue
  t=time.time();value=one(targets[i]);lfc[i]=value[0];prob[i]=value[1];valid[i]=True
  np.savez_compressed(path,lfc=lfc,p_de=prob,valid=valid,genes=inputs['genes'],positions=inputs['positions'],measured=inputs['measured'],targets=inputs['targets'])
  write(path.with_suffix('.json'),{'sha256':sha(path),'plan_sha256':sha(A/'plan.json'),'asset_manifest_sha256':sha(A/'assets.json'),'inputs_sha256':sha(A/'inputs'/f'{family}.npz'),'rows_complete':int(valid.sum()),'mode':mode,'precision':'FP32','source_context':'Missing, not inferred identity','original_evidence_key':stats.evidence_key,'original_train_split_sha256':stats.train_json_sha256});print('PIE',family,targets[i],round(time.time()-t,2),'seconds',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['smoke','development','confirmation','production']);p.add_argument('--family',default='H1');p.add_argument('--device',default='cuda:0');a=p.parse_args();main(a.mode,a.family,a.device)
