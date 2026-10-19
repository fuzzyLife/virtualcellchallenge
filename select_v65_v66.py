"""Apply the frozen V65/V66 comparison rule after both fixed screens complete."""
from pathlib import Path
import json
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v65_closed_response/analysis'
FOLDERS={65:'v65_closed_response',66:'v66_assay_strength'}

def main():
 assert not (A/'release_selection.json').exists();records=[]
 for version,folder in FOLDERS.items():
  a=ROOT/'experiments'/folder/'analysis';mean=read(a/'mean_results.json')
  if mean['selected'] is None:records.append(dict(version=version,eligible=False,reason='Failed frozen mean screen'));continue
  r=read(a/'count_assessment.json');assert r['complete'];arm=r['selected'];b=r['summary']['v63'];p=r['summary'][arm];improved=[]
  for metric in ('pds','fidelity','reach','jaccard','nmae','profile_square_error'):
   delta=(p[metric]-b[metric])*(-1 if metric in ('nmae','profile_square_error') else 1)
   if delta>1e-12:improved.append(metric)
  study_ratios={s:r['studies'][arm][s]['profile_square_error']/r['studies']['v63'][s]['profile_square_error'] for s in ('H1','Orion','CD4')}
  diagnostic_eligible=max(study_ratios.values())<=1.03 and len(improved)>=2
  records.append(dict(version=version,candidate=arm,eligible=r['release_supported'] or diagnostic_eligible,all_public_guards_passed=r['release_supported'],failed_guards=[k for k,v in r['checks'].items() if not v],improved_raw_metrics=improved,study_profile_error_ratios=study_ratios,profile_error_ratio=r['paired_equal_study_profile_error_ratio'],count_assessment_sha256=sha(a/'count_assessment.json'),mean_result_sha256=sha(a/'mean_results.json')))
 passing=[r for r in records if r.get('all_public_guards_passed')];pool=passing or [r for r in records if r['eligible']]
 selected=min(pool,key=lambda r:r['profile_error_ratio']) if pool else None
 write(A/'release_selection.json',dict(selected_version=selected['version'] if selected else None,selected=selected,records=records,selection_code_sha256=sha(__file__),rule_plan_sha256=sha(ROOT/'experiments/v66_assay_strength/analysis/plan.json'),experimental=not bool(passing),rule='Prefer a finalist passing all predeclared guards, then smaller equal-study generated-profile MSE ratio. If none passes, apply predeclared diagnostic bounds and >=2 raw improvements; never substitute another strength.',private_score=None,automatic_upload=False))
 print(json.dumps(dict(selected_version=selected['version'] if selected else None,records=records),indent=2),flush=True)
if __name__=='__main__':main()
