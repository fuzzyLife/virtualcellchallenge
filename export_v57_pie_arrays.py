"""Export completed, checksum-verified V57 arrays from temporary cache to project."""
import shutil
from v57_pie import *

def main():
 assert read(A/'final_audit.json')['passed']
 destination=A/'arrays';destination.mkdir(exist_ok=True);records=[]
 for folder in ['inputs','predictions','calibrators','count_statistics']:
  for p in sorted((CACHE/folder).glob('*.npz')):
   target=destination/folder/p.name;target.parent.mkdir(parents=True,exist_ok=True);digest=sha(p)
   if target.exists():assert sha(target)==digest
   else:shutil.copyfile(p,target);assert sha(target)==digest
   records.append({'file':str(target.relative_to(ROOT)).replace('\\','/'),'bytes':target.stat().st_size,'sha256':digest})
 write(A/'array_exports.json',{'passed':True,'files':records,'bytes':sum(r['bytes'] for r in records),'exporter_sha256':sha(ROOT/'export_v57_pie_arrays.py')})
 print('Exported',len(records),'arrays,',round(sum(r['bytes'] for r in records)/1024**3,2),'GiB to project',flush=True)

if __name__=='__main__':main()
