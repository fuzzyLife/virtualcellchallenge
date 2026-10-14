"""Reproduce Stack inference in a standalone local Python 3.12 environment.

This never modifies the competition scoring runtime or the live shared PIE env.
"""
from pathlib import Path
import argparse,io,json,subprocess,tempfile,zipfile,sys
ROOT=Path(__file__).resolve().parent
FREEZE=ROOT/'experiments/v20_stack_transfer/analysis/inference_runtime_freeze.txt'

def run(args):return subprocess.run([str(x) for x in args],check=True,text=True,encoding='utf-8',errors='replace')
def verify(python):
 result=subprocess.run([str(python),'-m','pip','freeze'],check=True,capture_output=True,text=True,encoding='utf-8',errors='replace')
 found=set(result.stdout.splitlines());expected=set(FREEZE.read_text(encoding='utf-8-sig').splitlines());missing=expected-found
 if missing:raise RuntimeError('Runtime pins differ: '+repr(sorted(missing)))
 run([python,'-c','import sys,torch;assert sys.version_info[:2]==(3,12);assert torch.cuda.is_available();print(sys.version);print(torch.__version__);print(torch.cuda.get_device_name(0))'])
 print('All frozen Stack package pins and CUDA verified:',python,flush=True)
def main(args):
 root=Path(args.runtime_dir).resolve();python=root/'env/Scripts/python.exe'
 if args.verify_runtime:return verify(Path(args.verify_runtime))
 if args.install:
  import requests
  base=root/'python312'
  if not (base/'python.exe').exists():
   response=requests.get('https://www.nuget.org/api/v2/package/python/3.12.10',timeout=180);response.raise_for_status();archive=zipfile.ZipFile(io.BytesIO(response.content))
   for name in archive.namelist():
    if not name.startswith('tools/') or name.endswith('/'):continue
    rel=Path(name[6:]);assert '..' not in rel.parts;out=base/rel;out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(archive.read(name))
  if not python.exists():run([base/'python.exe','-m','venv',root/'env'])
  lines=FREEZE.read_text(encoding='utf-8-sig').splitlines();torchpin=next(v for v in lines if v.startswith('torch=='))
  run([python,'-m','pip','install','--no-compile',torchpin,'--index-url','https://download.pytorch.org/whl/cu126'])
  with tempfile.NamedTemporaryFile(mode='w',suffix='.txt',delete=False) as stream:stream.write('\n'.join(v for v in lines if not v.startswith('torch=='))+'\n');requirements=Path(stream.name)
  try:run([python,'-m','pip','install','--no-compile','-r',requirements])
  finally:requirements.unlink()
 if not python.exists():raise FileNotFoundError('Use --install to create '+str(python))
 verify(python)
 print('Next: collect_v20_stack_assets.py, test_v20_mean_adapter.py, then the frozen input/inference/scoring workflow.',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--install',action='store_true');p.add_argument('--runtime-dir',default=str(Path(tempfile.gettempdir())/'vcc_v20_reproduction'));p.add_argument('--verify-runtime');main(p.parse_args())
