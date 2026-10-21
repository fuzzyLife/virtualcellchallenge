"""Inspect the user's BOLERO clone, paper, and public model metadata without executing upstream code."""
from pathlib import Path
import hashlib,json,subprocess,urllib.request,sys
ROOT=Path(__file__).resolve().parent
A=ROOT/'experiments/v87_conditional_adapter/analysis'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def get(url):
 req=urllib.request.Request(url,headers={'User-Agent':'VCC-reproducible-research'})
 with urllib.request.urlopen(req,timeout=45) as r:return json.load(r)
def main():
 A.mkdir(parents=True,exist_ok=True);dest=A/'bolero_inspection.json';assert not dest.exists()
 repo=ROOT/'bolero';safe='safe.directory='+repo.as_posix()
 def git(*args):return subprocess.check_output(['git','-c',safe,'-C',str(repo),*args],text=True,encoding='utf-8').strip()
 hf=get('https://huggingface.co/api/models/arcinstitute/bolero-models');rev=hf['sha'];files=get('https://huggingface.co/api/models/arcinstitute/bolero-models/tree/'+rev+'?recursive=true&limit=1000');configs={}
 for name in ['BorzoiSignalMultiDS/WholeHuman+AP1-Gene/config.json','BorzoiSignalMultiDS/Borzoi10M-Gene/config.json']:
  url='https://huggingface.co/arcinstitute/bolero-models/resolve/'+rev+'/'+name;d=get(url);p=A/('bolero_'+name.split('/')[1]+'_config.json');p.write_text(json.dumps(d,indent=2),encoding='utf-8');configs[name]=dict(sha256=sha(p),source=url,keys=list(d),datasets=list(d.get('dataset_records',{})))
 record=dict(paper_sha256=sha(ROOT/'geneRegPred.pdf'),paper_bytes=(ROOT/'geneRegPred.pdf').stat().st_size,paper_url='https://www.biorxiv.org/content/10.64898/2026.10.19.760026v1',clone_head=git('rev-parse','HEAD'),clone_status=git('status','--porcelain'),recent_commits=git('log','-8','--format=%h %s').splitlines(),upstream_sources={name:sha(repo/name) for name in ['AGENTS.md','README.md','docs/installation.md','src/bolero/tl/generic/module_lora_cond.py','src/bolero/tl/model/borzoi/model_lora.py','src/bolero/tl/model/borzoi/module_output.py']},hf_revision=rev,checkpoints=[dict(path=f['path'],bytes=f.get('size'),lfs=f.get('lfs')) for f in files if f.get('path','').endswith('.pt')],configs=configs,scope='Public metadata/configuration inspection only; upstream model weights not downloaded or executed. Local PDF is authoritative over the earlier README atlas counts.',script_sha256=sha(__file__))
 dest.write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps({k:record[k] for k in ['clone_head','hf_revision','checkpoints']},indent=2))
if __name__=='__main__':main()
