"""Verify staged scientific bytes match frozen dependencies and current evidence."""
from pathlib import Path
import json,hashlib,subprocess
ROOT=Path(__file__).resolve().parent
DIRS=['experiments/v18_pie_expert/analysis','experiments/v19_covariance_forces/analysis','experiments/v20_stack_transfer/analysis','experiments/v21_stack_programs/analysis']
def git(*args,input=None):
 return subprocess.run(['git','-c','safe.directory='+str(ROOT).replace('\\','/'),'-C',str(ROOT),*args],input=input,check=True,capture_output=True).stdout
def main():
 entries={}
 for row in git('ls-files','--stage','-z').split(b'\0'):
  if not row:continue
  metadata,path=row.split(b'\t',1);mode,oid,stage=metadata.split();assert stage==b'0';entries[path.decode()]=oid.decode()
 required={};plans={}
 for folder in DIRS:
  p=ROOT/folder/'plan.json';plan=json.loads(p.read_text());plans[str(p.relative_to(ROOT)).replace('\\','/')]=hashlib.sha256(p.read_bytes()).hexdigest()
  for name,value in plan['dependencies'].items():
   if name in required:assert required[name]==value,'Cross-plan dependency differs: '+name
   required[name]=value
 production=ROOT/'experiments/v20_stack_transfer/analysis/production_plan.json'
 if production.exists():
  plan=json.loads(production.read_text());plans[str(production.relative_to(ROOT)).replace('\\','/')]=hashlib.sha256(production.read_bytes()).hexdigest()
  for name,value in plan['dependencies'].items():
   if name in required:assert required[name]==value
   required[name]=value
 chosen=set(required)|set(plans)
 for name in entries:
  if name=='experiments/v20_stack_transfer/analysis/staged_scientific_bytes_audit.json':continue
  if any(name.startswith(folder+'/') for folder in DIRS) or (name.endswith('.py') and any(v in name for v in ['v18','v19','v20','v21'])):chosen.add(name)
 chosen=sorted(chosen);assert set(chosen)<=set(entries);output=git('cat-file','--batch',input=('\n'.join(entries[n] for n in chosen)+'\n').encode());cursor=0;checks={}
 for name in chosen:
  end=output.index(b'\n',cursor);oid,kind,size=output[cursor:end].split();assert kind==b'blob';start=end+1;raw=output[start:start+int(size)];cursor=start+int(size)+1;digest=hashlib.sha256(raw).hexdigest();assert digest==hashlib.sha256((ROOT/name).read_bytes()).hexdigest(),'Staged bytes differ: '+name
  if name in required:assert digest==required[name],'Frozen staged dependency differs: '+name
  if name in plans:assert digest==plans[name]
  checks[name]=digest
 assert cursor==len(output)
 record={'passed':True,'staged_files_checked':len(chosen),'unique_frozen_dependencies_checked':len(required),'frozen_plans_sha256':plans,'file_sha256':checks,'large_regenerable_arrays_included':False,'auditor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 p=ROOT/'experiments/v20_stack_transfer/analysis/staged_scientific_bytes_audit.json';p.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');print('Staged frozen scientific bytes verified',len(chosen),'files',flush=True)
if __name__=='__main__':main()
