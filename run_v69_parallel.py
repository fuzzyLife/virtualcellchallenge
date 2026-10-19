"""Run independent fixed V69 context jobs with two bounded worker processes."""
import argparse,subprocess,sys,time
from pathlib import Path
from run_v35_v36 import contexts
from v69_count_bridge import ROOT,OUT


def main(stage):
    jobs=[c['family'] for c in contexts()];active={}
    try:
        while jobs or active:
            while jobs and len(active)<2:
                name=jobs.pop(0);path=OUT/f'{stage}_{name}.log';log=path.open('w',encoding='utf-8')
                proc=subprocess.Popen([sys.executable,'-X','utf8',str(ROOT/'run_v69_cached.py'),stage,'--source',name],stdout=log,stderr=subprocess.STDOUT,cwd=ROOT)
                active[name]=(proc,log);print('V69 started',stage,name,flush=True)
            for name,(proc,log) in list(active.items()):
                if proc.poll() is not None:
                    log.close();del active[name]
                    if proc.returncode:raise RuntimeError(f'{stage} {name} failed: see {OUT / (stage+"_"+name+".log")}')
                    print('V69 completed',stage,name,flush=True)
            if jobs or active:time.sleep(2)
    finally:
        for proc,log in active.values():
            if proc.poll() is None:proc.terminate()
            log.close()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation']);a=p.parse_args();main(a.stage)
