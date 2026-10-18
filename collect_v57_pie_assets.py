"""Fetch hash-pinned V57 PIE assets into a resumable local temporary cache.

Only four evidence label arrays are required, not original expression matrices.
Existing V18 source embeddings are reused read-only. No model inference here.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os,tempfile,time
import requests
ROOT=Path(__file__).resolve().parent
A=ROOT/'experiments/v57_pie_cross/analysis'
CACHE=Path(tempfile.gettempdir())/'vcc_v57_pie'
REV='f6da9607226139a1040b6992736118c9b5c7c313'
DATA={'arc_vcc_25':'e7dda6065316959958b87583ddf0a15f6e69c37d','jiang':'99242356dfcf58e049d867703c5743df8d1790bf','tahoe100m':'597d6ded9f9a148331388484ffb193489121c702','x_atlas_orion':'0950c8aa5ce5b33dc9d6f0c5eade3ae445ca0ffb'}

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for v in iter(lambda:f.read(16*1024**2),b''):h.update(v)
 return h.hexdigest()

def write(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')

def fetch(item):
 rel,url,size,digest=item;p=CACHE/rel;p.parent.mkdir(parents=True,exist_ok=True)
 receipt=p.with_suffix(p.suffix+'.verified.json')
 if p.exists() and receipt.exists():
  r=json.loads(receipt.read_text());assert p.stat().st_size==size and r['sha256']==digest and r['bytes']==size
  return r
 partial=p.with_suffix(p.suffix+'.partial')
 for attempt in range(6):
  try:
   offset=partial.stat().st_size if partial.exists() else 0
   headers={'Range':f'bytes={offset}-'} if offset else {}
   with requests.get(url+f'?v57_offset={offset}&attempt={attempt}',headers=headers,stream=True,timeout=(30,120)) as r:
    r.raise_for_status()
    if offset:assert r.status_code==206 and r.headers.get('Content-Range','').startswith(f'bytes {offset}-'), 'Resume not honored'
    count=offset
    with partial.open('ab' if offset else 'wb') as f:
     last=time.monotonic()
     for block in r.iter_content(4*1024**2):
      count+=len(block);assert count<=size;f.write(block)
      if time.monotonic()-last>30:print('Download',rel,round(count/1024**3,2),'/',round(size/1024**3,2),'GiB',flush=True);last=time.monotonic()
   assert partial.stat().st_size==size
   assert sha(partial)==digest, 'Downloaded SHA256 mismatch: '+rel
   partial.replace(p);rec={'path':rel,'url':url,'bytes':size,'sha256':digest};write(receipt,rec)
   print('Verified',rel,round(size/1024**2,1),'MiB',flush=True);return rec
  except Exception as exc:
   print('Retry',rel,attempt+1,str(exc)[:160],flush=True)
   if attempt==5:raise
   time.sleep(min(2**attempt,15))

def main():
 A.mkdir(parents=True,exist_ok=True);CACHE.mkdir(exist_ok=True)
 items=[]
 model=json.loads((ROOT/'analysis/v56_pie_release/xd_current.json').read_text());r=next(v for v in model['siblings'] if v['rfilename']=='best_auprc.ckpt')
 items.append(('best_auprc.ckpt',f'https://huggingface.co/arcinstitute/PIE_replogle_xdataset/resolve/{REV}/best_auprc.ckpt',r['size'],r['lfs']['sha256']))
 for name,rev in DATA.items():
  repo='arcinstitute/PIE_'+name;meta=json.loads((ROOT/f'analysis/v56_pie_release/{name}_meta.json').read_text())
  folder=CACHE/'data'/meta['dataset'];folder.mkdir(parents=True,exist_ok=True)
  (folder/'meta.json').write_bytes((ROOT/f'analysis/v56_pie_release/{name}_meta.json').read_bytes())
  url=f'https://huggingface.co/api/datasets/{repo}/tree/{rev}/preprocessed?recursive=false'
  response=requests.get(url,timeout=60);response.raise_for_status();tree={v['path'].split('/')[-1]:v for v in response.json() if v['type']=='file'}
  write(A/(meta['dataset']+'_files.json'),{'url':url,'files':list(tree.values())})
  for file in ['ctx_ids.npy','pert_ids.npy','ctrl_means.npy','fold_changes.npy','fdr.npy','delta_p.npy','tested.npy']:
   rec=tree[file];assert rec.get('lfs',{}).get('oid',meta['array_sha256'][file])==meta['array_sha256'][file]
   items.append((f'data/{meta["dataset"]}/{file}',f'https://huggingface.co/datasets/{repo}/resolve/{rev}/preprocessed/{file}',rec['size'],meta['array_sha256'][file]))
 files=json.loads((ROOT/'analysis/post_v17_research/PIE_sources_files.json').read_text());files={v['path']:v for v in files if v['type']=='file'}
 sr='cb1aaa4e7655605bdc70a9bd77bbd62016b8c7d7'
 for name in ['smiles','l1000_tas','prism_secondary','jump_morphology']:
  for file in ['meta.json','embeddings.npy']:
   rel=name+'/'+file;rec=files[rel];digest=rec.get('lfs',{}).get('oid')
   if digest is None:
    response=requests.get(f'https://huggingface.co/datasets/arcinstitute/PIE_sources/resolve/{sr}/{rel}',timeout=60);response.raise_for_status();assert len(response.content)==rec['size'];digest=hashlib.sha256(response.content).hexdigest()
   items.append(('sources/'+rel,f'https://huggingface.co/datasets/arcinstitute/PIE_sources/resolve/{sr}/{rel}',rec['size'],digest))
 write(A/'download_plan.json',{'model_revision':REV,'data_revisions':DATA,'source_revision':sr,'bytes':sum(v[2] for v in items),'files':[{'path':v[0],'url':v[1],'bytes':v[2],'sha256':v[3]} for v in items]})
 print('Asset plan',len(items),'files',round(sum(v[2] for v in items)/1024**3,2),'GiB',flush=True)
 with ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(fetch,items))
 write(A/'download_receipts.json',{'passed':True,'files':records,'local_cache':str(CACHE),'collector_sha256':sha(Path(__file__))})
 print('All V57 assets verified',flush=True)

if __name__=='__main__':main()
