"""Complete independent serialized-result audits, sensitivity and the frozen release decision."""
import sys,time,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
for version,folder in [(65,'v65_closed_response'),(66,'v66_assay_strength')]:
 a=ROOT/'experiments'/folder/'analysis'
 while not (a/'mean_results.json').exists():time.sleep(3)
 mean=json.loads((a/'mean_results.json').read_text())
 if not mean['selected']:continue
 while not (a/'count_assessment.json').exists():time.sleep(3)
 time.sleep(2)
 for cmd in [[f'audit_v{version}_results.py'],['analyze_v65_sensitivity.py','--version',str(version)]]:
  subprocess.run([sys.executable,'-X','utf8',str(ROOT/cmd[0]),*cmd[1:]],check=True)
subprocess.run([sys.executable,'-X','utf8',str(ROOT/'select_v65_v66.py')],check=True)
