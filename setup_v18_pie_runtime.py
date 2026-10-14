"""Recreate pinned PIE source and its isolated local Python 3.12 runtime."""
from pathlib import Path
import argparse,hashlib,io,json,os,subprocess,tempfile,zipfile
import requests
ROOT=Path(__file__).resolve().parent;R=ROOT/'reference/v18_pie';REV='21a526add9ff3a81b4f962d9b16bbae570052c63'
def request(url):r=requests.get(url,timeout=180);r.raise_for_status();return r.content
def main(install):
 R.mkdir(parents=True,exist_ok=True);code=R/'code'
 if not code.exists():
  z=zipfile.ZipFile(io.BytesIO(request('https://codeload.github.com/ArcInstitute/pie/zip/'+REV)));prefix=z.namelist()[0]
  for n in z.namelist():
   if not n.startswith(prefix) or n.endswith('/'):continue
   rel=Path(n[len(prefix):]);assert '..' not in rel.parts;p=code/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(z.read(n))
 tree=json.loads((ROOT/'analysis/post_v17_research/pie_tree.json').read_text())['tree']
 for item in tree:
  if item['type']=='blob' and item['path'].startswith('src/pie/'):
   raw=(code/item['path']).read_bytes();assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==item['sha']
 split=R/'train.json'
 if not split.exists():split.write_bytes(request('https://huggingface.co/datasets/arcinstitute/PIE_splits/resolve/396ab9563175ee887750c9eed7ccaea6f5fdbf50/replogle_wdataset/unseen_ctx/k562/train.json'))
 assert hashlib.sha256(split.read_bytes()).hexdigest()=='6143b4d1bdc30f09a323a5af944fe9049c77c724909fb20fdde0039b7faeddfb'
 runtime=Path(tempfile.gettempdir())/'vcc_v18_local';python=runtime/'env/Scripts/python.exe'
 if install:
  base=runtime/'python312'
  if not (base/'python.exe').exists():
   z=zipfile.ZipFile(io.BytesIO(request('https://www.nuget.org/api/v2/package/python/3.12.10')))
   for n in z.namelist():
    if not n.startswith('tools/') or n.endswith('/'):continue
    rel=Path(n[6:]);assert '..' not in rel.parts;p=base/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(z.read(n))
  if not python.exists():subprocess.run([str(base/'python.exe'),'-m','venv',str(runtime/'env')],check=True)
  subprocess.run([str(python),'-m','pip','install','--no-compile','torch==2.10.0','--index-url','https://download.pytorch.org/whl/cu126'],check=True)
  freeze=ROOT/'experiments/v18_pie_expert/analysis/inference_runtime_freeze.txt'
  if freeze.exists():
   # CUDA wheel installed above. Requirements pins all remaining direct/transitive versions.
   import tempfile as tf
   lines=[v for v in freeze.read_text(encoding='utf-8-sig').splitlines() if not v.startswith('torch==')]
   with tf.NamedTemporaryFile(mode='w',suffix='.txt',delete=False) as f:f.write('\n'.join(lines));requirements=f.name
   try:subprocess.run([str(python),'-m','pip','install','--no-compile','-r',requirements],check=True)
   finally:Path(requirements).unlink()
  else:
   import tomllib
   deps=tomllib.loads((code/'pyproject.toml').read_text())['project']['dependencies'];deps=[v for v in deps if not v.startswith('torch')];subprocess.run([str(python),'-m','pip','install','--no-compile',*deps],check=True)
 if not python.exists():raise FileNotFoundError('Run --install to create '+str(python))
 if (R/'sources/ncbi_subset/meta.json').exists():subprocess.run([str(python),str(ROOT/'test_v18_pie.py'),'--upstream'],check=True)
 else:print('Next: collect_v18_pie_assets.py, then test_v18_pie.py --upstream',flush=True)
 print('Use isolated inference Python:',python,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--install',action='store_true');main(p.parse_args().install)
