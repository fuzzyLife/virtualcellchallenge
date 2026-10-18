"""Verify reused V18 source arrays and code against their original byte receipts."""
from v57_pie import *

def main():
 old=read(ROOT/'experiments/v18_pie_expert/analysis/assets.json');checked=[]
 for rec in old['files']:
  name=rec['path']
  if not name.startswith(('reference/v18_pie/sources/','reference/v18_pie/code/src/')):continue
  p=ROOT/name;assert p.stat().st_size==rec['bytes'] and sha(p)==rec['sha256'],name
  checked.append({'path':name,'bytes':rec['bytes'],'sha256':rec['sha256']})
 write(A/'reused_asset_audit.json',{'passed':True,'prior_manifest_sha256':sha(ROOT/'experiments/v18_pie_expert/analysis/assets.json'),'checked':checked,'bytes_verified':sum(v['bytes'] for v in checked),'auditor_sha256':sha(Path(__file__))})
 print('Reused assets verified:',len(checked),'files,',sum(v['bytes'] for v in checked),'bytes',flush=True)

if __name__=='__main__':main()
