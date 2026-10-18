"""Apply the recorded V64 unique-label scoring amendment without changing predictors."""
from pathlib import Path
import json,shutil
from v57_pie import sha
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v64_response_magnitude/analysis'
path=A/'plan.json';p=json.loads(path.read_text())
if 'validation_correction' in p:
 assert p['dependencies']['audit_v64_unique_panels.py']==sha(ROOT/'audit_v64_unique_panels.py')
 print('Recorded V64 validation amendment already present; no files changed.')
else:
 assert not (A/'initial_plan.json').exists();shutil.copyfile(path,A/'initial_plan.json')
 p['validation_correction']='While V64 development was running, before confirmation: audit found44 case configurations span30 target/context pairs and34 distinct measured row sets. Existing development panels already have unique labels. For any confirmation/additional panel, retain canonical first case_index per target before PDS/metric aggregation. Original evaluator outputs remain unchanged, with separate unique-panel receipts. No predictor or parameter changes. Older V63 outcomes and mean-control diagnostics were known.'
 p['dependencies']['audit_v64_unique_panels.py']=sha(ROOT/'audit_v64_unique_panels.py');path.write_text(json.dumps(p,indent=2)+'\n',encoding='utf-8')
 print('Applied documented V64 scoring amendment; numerical predictors unchanged.')
