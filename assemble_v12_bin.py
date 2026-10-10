"""Assemble one complete range download to a separate, checksum-verified inspection copy."""
import argparse,json,hashlib,os
from pathlib import Path
ROOT=Path(__file__).resolve().parent;D=ROOT/'reference/v12_resources'
def main(name):
 r=json.loads((ROOT/'analysis/post_v11_research/data/viperturb_record.json').read_text());f=next(f for f in r['files'] if f['key']==name)
 folder=D/'viperturb'/(name+'.ranges');out=D/'inspected_bins'/name;out.parent.mkdir(exist_ok=True)
 size=8*1024**2;parts=[]
 for start in range(0,f['size'],size):
  p=folder/f'{start:012d}.bin';assert p.exists() and p.stat().st_size==min(size,f['size']-start),(name,start);parts.append(p)
 md5=hashlib.md5();sha=hashlib.sha256();pending=out.with_suffix('.assembling')
 with pending.open('wb') as h:
  for p in parts:
   b=p.read_bytes();h.write(b);md5.update(b);sha.update(b)
 assert 'md5:'+md5.hexdigest()==f['checksum'] and pending.stat().st_size==f['size']
 os.replace(pending,out);out.with_suffix('.json').write_text(json.dumps({'file':str(out.relative_to(ROOT)),'bytes':f['size'],'md5':md5.hexdigest(),'sha256':sha.hexdigest()},indent=2)+'\n')
 print('Verified assembled copy',out,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('name');main(p.parse_args().name)
