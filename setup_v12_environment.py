"""Recreate the isolated CPU scorer and optional project-local R runtime."""
import argparse,subprocess,sys,json,hashlib
from pathlib import Path
import requests
from collect_v12_resources import scorer
ROOT=Path(__file__).resolve().parent;D=ROOT/'reference/v12_resources/tools'
def main(python313,with_r=False):
 scorer();env=ROOT/'.v12_eval_env';subprocess.run([python313,'-m','venv',str(env)],check=True)
 python=env/'Scripts/python.exe'
 subprocess.run([str(python),'-m','pip','install','-r',str(ROOT/'v12_official_requirements.txt'),str(ROOT/'reference/v12_resources/scorer')],check=True)
 freeze=subprocess.check_output([str(python),'-m','pip','freeze'],text=True)
 (D.parent/'runtime_freeze_recreated.txt').write_text(freeze,encoding='utf-8')
 if with_r:
  D.mkdir(parents=True,exist_ok=True);url='https://cran.r-project.org/bin/windows/base/R-4.6.1-win.exe';path=D/'R-4.6.1-win.exe'
  r=requests.get(url,timeout=120);r.raise_for_status();path.write_bytes(r.content)
  md5=hashlib.md5(r.content).hexdigest();expected=requests.get('https://cran.r-project.org/bin/windows/base/md5sum.R-4.6.1.txt',timeout=60);expected.raise_for_status();assert md5 in expected.text
  dest=D/'R';subprocess.run([str(path),'/VERYSILENT','/SUPPRESSMSGBOXES','/NORESTART','/CURRENTUSER','/NOICONS','/TASKS=',f'/DIR={dest}'],check=True)
  (D/'runtime.json').write_text(json.dumps({'url':url,'md5':md5,'install_dir':str(dest)},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--python313',required=True);p.add_argument('--with-r',action='store_true');a=p.parse_args();main(a.python313,a.with_r)
