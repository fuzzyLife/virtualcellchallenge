"""Finalize documentation and compact scientific seal only after the full payload audit."""
import sys,time,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent;decision=ROOT/'experiments/v65_closed_response/analysis/release_selection.json'
while not decision.exists():time.sleep(3)
time.sleep(2);d=json.loads(decision.read_text());v=d['selected_version']
if v is None:print('No qualifying release; preserve experiment results.',flush=True);sys.exit(0)
folder='v65_closed_response' if v==65 else 'v66_assay_strength';a=ROOT/'experiments'/folder/'analysis'
while not (a/'release_complete.json').exists():time.sleep(3)
time.sleep(2)
for args in [['document_v65_v66.py'],['seal_v65_v66_release.py'],['seal_v65_v66_release.py','--check']]:subprocess.run([sys.executable,'-X','utf8',str(ROOT/args[0]),*args[1:]],check=True)
print('V'+str(v)+' documentation and scientific seal complete; ready for local Git commit.',flush=True)
