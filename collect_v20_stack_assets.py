"""Retrieve pinned Stack model/code for the second pretrained population pilot."""
from pathlib import Path
import json,hashlib,io,zipfile,requests,pickle
from collect_v18_pie_assets import fetch
from audit_post_v17_feasibility import decode_unicode_list
ROOT=Path(__file__).resolve().parent;R=ROOT/'reference/v20_stack';REV='cacc2e4b09435c3e536d46237d10b50f222dd144'
def main():
 R.mkdir(parents=True,exist_ok=True);code=R/'code'
 if not code.exists():
  r=requests.get('https://codeload.github.com/ArcInstitute/stack/zip/'+REV,timeout=180);r.raise_for_status();z=zipfile.ZipFile(io.BytesIO(r.content));prefix=z.namelist()[0]
  for n in z.namelist():
   if not n.startswith(prefix) or n.endswith('/'):continue
   rel=Path(n[len(prefix):]);assert '..' not in rel.parts;p=code/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(z.read(n))
 tree=json.loads((ROOT/'analysis/post_v17_research/stack_tree.json').read_text())['tree'];checked=0
 for item in tree:
  if item['type']=='blob' and item['path'].startswith('src/stack/'):
   raw=(code/item['path']).read_bytes();assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==item['sha'];checked+=1
 # Reencode previously non-executingly decoded vocabulary as a plain trusted list for author loader.
 original=ROOT/'analysis/post_v17_research/stack_gene_vocabulary.pkl';assert hashlib.sha256(original.read_bytes()).hexdigest()=='d8761dfda955b9897d3251798b72361ddd0171ef707119eaf65381bed2d85dcc';vocab=decode_unicode_list(original.read_bytes());assert len(vocab)==15012 and len(set(vocab))==15012
 (R/'genes.json').write_text(json.dumps(vocab)+'\n');(R/'genes.pkl').write_bytes(pickle.dumps(vocab,protocol=4));print('Stack author code/vocabulary verified',checked,flush=True)
 fetch('https://huggingface.co/arcinstitute/Stack-Large-Aligned/resolve/b09f085dac03d170b078a5c72f550ae93686e544/bc_large_aligned.ckpt',R/'bc_large_aligned.ckpt',2613863242,'f93cf6f42f36c8a85dc570d92e801c1fc3e1f45d55741e6bb250d178a6b6ad36')
if __name__=='__main__':main()
