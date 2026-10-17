"""Collect unused KOLF browser response features with pinned author blob checks.

These are quantized mean NTC-z scores in log-expression space: never raw counts,
fold changes, DE q-values, or cell-level variation training examples.
"""
import hashlib,json,time
from pathlib import Path
from collections import Counter
import requests,numpy as np
import build_submission as b
ROOT=b.ROOT
OUT=ROOT/'experiments/v39_kolf_support';A=OUT/'analysis';REF=OUT/'reference'
def main():
 A.mkdir(parents=True,exist_ok=True);REF.mkdir(exist_ok=True)
 old=ROOT/'experiments/v17_learned_heterogeneity/analysis'
 tree=json.loads((old/'kolf_repo_tree.json').read_text(encoding='utf-8'));meta=json.loads((old/'kolf_all_meta.json').read_text(encoding='utf-8'));revision=tree['sha']
 entries={e['path']:e for e in tree['tree']};_,genes,targets=b.read_spec(ROOT);records=[]
 wanted=['docs/data/panel/all/meta.json']+[f'docs/data/panel/all/chunk_{i}.bin' for i in range(meta['n_chunks'])]
 for name in wanted:
  ent=entries[name];path=REF/Path(name).name;url=f'https://raw.githubusercontent.com/y-doctor/KOLF2.1J_Perturbation_Cell_Atlas/{revision}/{name}'
  if not path.exists():
   temporary=path.with_suffix(path.suffix+'.partial');count=0
   with requests.get(url,stream=True,timeout=(30,90)) as r:
    r.raise_for_status()
    with temporary.open('wb') as f:
     for block in r.iter_content(1024*1024):
      count+=len(block);assert count<=ent['size'];f.write(block)
   assert count==ent['size'];temporary.replace(path)
  payload=path.read_bytes();assert len(payload)==ent['size'];blob=hashlib.sha1(f'blob {len(payload)}\0'.encode()+payload).hexdigest();assert blob==ent['sha'],name
  records.append({'name':name,'url':url,'bytes':len(payload),'git_blob_sha1':blob,'sha256':hashlib.sha256(payload).hexdigest()})
  print('KOLF verified',Path(name).name,len(payload),flush=True)
 assert json.loads((REF/'meta.json').read_text(encoding='utf-8'))==meta
 assert len(set(meta['perts']))==len(meta['perts'])
 counts=Counter(meta['genes']);gmap={g:i for i,g in enumerate(meta['genes']) if counts[g]==1};tmap={t:i for i,t in enumerate(meta['perts'])}
 measured=np.array([g in gmap for g in genes]);present=np.array([t in tmap for t in targets]);z=np.full((len(targets),len(genes)),np.nan,np.float32);saturated=np.zeros_like(z,bool)
 pos=np.flatnonzero(measured);sourcepos=[gmap[genes[g]] for g in pos]
 for ti,t in enumerate(targets):
  if t not in tmap:continue
  index=tmap[t];ci=index//meta['chunk_size'];local=index%meta['chunk_size'];path=REF/f'chunk_{ci}.bin';nrows=min(meta['chunk_size'],meta['n_perts']-ci*meta['chunk_size'])
  mm=np.memmap(path,dtype=np.int8,mode='r',shape=(nrows,meta['n_genes']));raw=mm[local,sourcepos].copy();z[ti,pos]=raw.astype(np.float32)*meta['scale'];saturated[ti,pos]=abs(raw.astype(np.int16))==127;del mm
 file=A/'response_features.npz';assert not file.exists();np.savez_compressed(file,z=z,measured=measured,target_present=present,saturated=saturated,genes=np.asarray(genes,str),targets=np.asarray(targets,str))
 b.write_json(A/'collection.json',{'revision':revision,'author_sources':records,'arrays_sha256':b.sha256(file),'script_sha256':b.sha256(Path(__file__)),'target_overlap':int(present.sum()),'unambiguous_output_overlap':int(measured.sum()),'finite_feature_pairs':int(np.isfinite(z).sum()),'saturated_feature_pairs':int(saturated.sum()),'source_units':'int8 quantized mean NTC z-scores in author log-expression space, scale0.04; NOT counts or LFC.','use':'Optional experimental-response feature only. No training or submission use implied. All-perturbation browser export, not strong-only filtered file.','missingness':'Missing targets or ambiguous/missing gene symbols are NaN; no invented zero response.','full_raw_189GB_file_downloaded':False,'license':'Data CC BY4.0 per pinned Figshare metadata; author code MIT.','paper':'https://www.nature.com/articles/s41587-026-03199-w'})
 print('KOLF usable feature targets',int(present.sum()),'outputs',int(measured.sum()),flush=True)
if __name__=='__main__':main()

