"""Run the separate robust-consensus development stage after its screen completes."""
import json,subprocess,sys,time,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v13_response_programs/analysis/quality_stage'
def main():
 while not (A/'response_screen_result.json').exists():time.sleep(20)
 subprocess.run([sys.executable,str(ROOT/'prepare_v13_quality.py')],check=True)
 (A/'official_anchors').mkdir(exist_ok=True)
 for p in (A.parent/'official_anchors').iterdir():
  if p.is_file():shutil.copyfile(p,A/'official_anchors'/p.name)
 subprocess.run([sys.executable,str(ROOT/'evaluate_v13_quality.py'),'--stage','development','--workers','4'],check=True)
if __name__=='__main__':main()
