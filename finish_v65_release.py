"""Build and independently validate only the prospectively selected diagnostic release."""
from pathlib import Path
import json,time,subprocess,sys
ROOT=Path(__file__).resolve().parent;decision=ROOT/'experiments/v65_closed_response/analysis/release_selection.json'
while not decision.exists():time.sleep(3)
time.sleep(2);r=json.loads(decision.read_text());version=r['selected_version']
if version is None:
 print('Neither finalist meets the frozen public/experimental release rule; no unsupported production package built.',flush=True);sys.exit(0)
folder='v65_closed_response' if version==65 else 'v66_assay_strength';A=ROOT/'experiments'/folder/'analysis'
assert json.loads((A/'results_audit.json').read_text())['passed'];assert json.loads((A/'boundary_audit.json').read_text())['passed']
science=sys.executable;official=str(ROOT/'.competition_env/Scripts/python.exe')
commands=[(science,['audit_v65_count_replay.py','--version',str(version)]),(science,[f'fit_v{version}_submission.py','fit']+(['--experimental'] if r['experimental'] else [])),(science,[f'build_v{version}_submission.py','--preflight']),(science,[f'build_v{version}_submission.py','--workers','4']),(official,[f'package_v{version}_submission.py']),(science,[f'validate_v{version}_submission.py'])]
for interpreter,args in commands:
 print('Selected V'+str(version)+' release phase: '+' '.join(args),flush=True)
 subprocess.run([interpreter,'-X','utf8',str(ROOT/args[0]),*args[1:]],check=True)
ready=json.loads((A/'submission_ready.json').read_text());audit=json.loads((A/'independent_validation.json').read_text());assert audit['passed'] and audit['package_sha256']==ready['sha256']
(A/'release_complete.json').write_text(json.dumps(dict(version=version,candidate=r['selected']['candidate'],independent_integrity_passed=True,package_sha256=ready['sha256'],package_bytes=ready['bytes'],experimental=r['experimental'],uploaded=False,private_score=None),indent=2)+'\n',encoding='utf-8')
print('Selected V'+str(version)+' package and independent audit complete. Nothing uploaded.',flush=True)
