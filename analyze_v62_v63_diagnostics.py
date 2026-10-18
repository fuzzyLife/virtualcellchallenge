"""Separate stable count-profile errors from weak-reference expression ratios."""
from pathlib import Path
import csv,json
import numpy as np
from v57_pie import read,write,sha
import build_submission as b
ROOT=Path(__file__).resolve().parent

def aggregate(rows,omit=None):
 values=[]
 for study in ('H1','Orion','CD4'):
  contexts=sorted({r['family'] for r in rows if r['study']==study});ratios=[]
  for context in contexts:
   rr=[r for r in rows if r['family']==context and (r['family'],r['target'])!=omit]
   if not rr:return None
   names=sorted({r['target'] for r in rr}); num=sum(np.mean([r['candidate_MSE'] for r in rr if r['target']==t]) for t in names); den=sum(np.mean([r['baseline_MSE'] for r in rr if r['target']==t]) for t in names);ratios.append(num/den)
  values.append(float(np.mean(ratios)))
 return float(np.mean(values))

def main():
 manifest,genes,targets=b.read_spec(ROOT);own={str(t):list(genes).index(t) for t in targets};a=ROOT/'experiments/v63_identity_guard/analysis';case_rows=[];inputs={};cohort_path=ROOT/'experiments/v34_residual_flow/analysis/plan.json';cohort=read(cohort_path);case_meta={(c['family'],i):v for c in cohort['contexts'] for i,v in enumerate(c['cells'])};inputs[str(cohort_path.relative_to(ROOT))]=sha(cohort_path)
 for stage in ('development','additional','confirmation'):
  folder=a/stage
  if not folder.exists():continue
  for path in sorted(folder.glob('*__guarded_assay.json')):
   r=read(path);original=folder/path.name.replace('__guarded_assay.json','__v24.json');base=read(original);inputs[str(path.relative_to(ROOT))]=sha(path);inputs[str(original.relative_to(ROOT))]=sha(original)
   errors=[]
   for metadata,file in [(r,path.with_suffix('.npz')),(base,original.with_suffix('.npz'))]:
    assert sha(file)==metadata['arrays_sha256'];inputs[str(file.relative_to(ROOT))]=metadata['arrays_sha256']
    with np.load(file,allow_pickle=False) as z:delta=z['pred']-z['truth']
    errors.append(float(np.mean(np.delete(delta,own[r['target']])**2)))
   case_rows.append(dict(stage=stage,family=r['family'],study=r['study'],target=r['target'],case_index=r['case_index'],actual_cells=len(case_meta[(r['family'],r['case_index'])]['actual']),source_case_seed=case_meta[(r['family'],r['case_index'])]['seed'],seed_offset=r['seed_offset'],candidate_MSE=errors[0],baseline_MSE=errors[1],ratio=errors[0]/errors[1]))
 assert len(case_rows)==68,('Expected 20 development + 24*2 additional populations',len(case_rows))
 dest=a/'count_profile_cases.csv'
 with dest.open('w',encoding='utf-8',newline='') as stream:
  w=csv.DictWriter(stream,fieldnames=list(case_rows[0]));w.writeheader();w.writerows(case_rows)
 groups={}
 for source in sorted({r['family'] for r in case_rows}):
  z=[r for r in case_rows if r['family']==source];groups[source]=dict(targets=len({r['target'] for r in z}),populations=len(z),ratio_of_summed_profile_errors=sum(np.mean([r['candidate_MSE'] for r in z if r['target']==t]) for t in {r['target'] for r in z})/sum(np.mean([r['baseline_MSE'] for r in z if r['target']==t]) for t in {r['target'] for r in z}),improved_populations=sum(r['ratio']<1 for r in z),median_case_ratio=float(np.median([r['ratio'] for r in z])))
 omitted=[aggregate(case_rows,(r['family'],r['target'])) for r in case_rows];ratio=aggregate(case_rows)
 from audit_v64_unique_panels import inventory
 metadata=inventory();write(a/'count_profile_summary.json',dict(equal_study_profile_error_ratio=ratio,by_context=groups,source_case_inventory=metadata,leave_one_context_target_out_min=min(omitted),leave_one_context_target_out_max=max(omitted),predefined_case_configurations=44,unique_target_context_pairs=metadata['unique_target_context_pairs'],distinct_measured_row_sets=metadata['distinct_measured_row_sets'],count_seed_populations=68,scope='Plug-in pooled-profile MSE, own-target excluded, with paired controls and same realizations. Targets receive equal weight within each context after averaging repeated case/seed errors. Ratio to V24 prediction error avoids division by near-zero measured perturbation-to-control distance. Not official jackknife/anchor-scaled expression accuracy; seeds are not new biology.'))
 broad={}
 for version,folder in [('v62','v62_conserved_response'),('v63','v63_identity_guard')]:
  r=read(ROOT/'experiments'/folder/'analysis/broad_identity_results.json');candidate='assay_center50' if version=='v62' else 'guarded_assay';study=[]
  for s in ('H1','Orion','CD4'):
   rr=[v for v in r['records'] if v['study']==s];study.append(float(np.mean([v['metrics'][candidate]['pds']-v['metrics']['v24']['pds'] for v in rr])))
  broad[version]=dict(equal_study_raw_PDS_change=float(np.mean(study)),studies=dict(zip(('H1','Orion','CD4'),study)),target_context_pairs=sum(v['targets'] for v in r['records']),scope='Desired means against genuine measured count profiles; same fixed 300-panel gene exclusion. Target counts vary by source; this is not final integer-count leaderboard performance.')
 write(a/'broader_identity_summary.json',broad);write(a/'diagnostic_analysis_audit.json',dict(input_sha256=inputs,script_sha256=sha(__file__),outputs_sha256={p.name:sha(p) for p in [dest,a/'count_profile_summary.json',a/'broader_identity_summary.json']}))
 print(json.dumps(dict(profile=ratio,contexts=groups,broad=broad),indent=2),flush=True)

if __name__=='__main__':main()
