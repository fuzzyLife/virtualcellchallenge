"""Audit source-case replication and score panels with unique perturbation labels."""
from pathlib import Path
import argparse,json
import numpy as np
from v57_pie import read,write,sha
from v8_methods import panel_metrics
from run_v37_ensemble import aggregate
import build_submission as b
ROOT=Path(__file__).resolve().parent
METRICS=('pds','fidelity','reach','jaccard','nmae','common_expression_ratio')

def inventory():
 plan=read(ROOT/'experiments/v34_residual_flow/analysis/plan.json');records=[]
 for c in plan['contexts']:
  cells=c['cells'];records.append(dict(context=c['family'],case_configurations=len(cells),unique_target_context_pairs=len({x['target'] for x in cells}),distinct_measured_row_sets=len({(x['target'],tuple(sorted(x['actual']))) for x in cells}),actual_cell_counts=sorted({len(x['actual']) for x in cells})))
 return dict(contexts=records,case_configurations=sum(r['case_configurations'] for r in records),unique_target_context_pairs=sum(r['unique_target_context_pairs'] for r in records),distinct_measured_row_sets=sum(r['distinct_measured_row_sets'] for r in records),independent_new_studies=0,interpretation='Distinct row sets can overlap and are not biological replicates. HCT116 and HEK293T repeat exactly the same measured treated rows in their two original case configurations. HEK additional panels contain duplicate target labels unless explicitly deduplicated.')

def score(version,stage):
 folder={63:'v63_identity_guard',64:'v64_response_magnitude'}[version];a=ROOT/'experiments'/folder/'analysis';candidate='guarded_assay' if version==63 else 'magnitude';dest=a/stage
 manifest,genes,targets=b.read_spec(ROOT);records=[(p,read(p)) for p in sorted(dest.glob('*.json')) if not p.name.endswith('__panel.json')];panels=[];choices=[]
 for source in sorted({r['family'] for p,r in records}):
  for offset in sorted({r['seed_offset'] for p,r in records if r['family']==source}):
   for arm in ('v24',candidate,candidate+'_shuffled'):
    rr=sorted([(p,r) for p,r in records if r['family']==source and r['seed_offset']==offset and r['arm']==arm],key=lambda x:x[1]['case_index']);chosen={}
    for p,r in rr:chosen.setdefault(r['target'],(p,r))
    vals=[]
    for p,r in chosen.values():
     arr=p.with_suffix('.npz');assert sha(arr)==r['arrays_sha256']
     with np.load(arr,allow_pickle=False) as z:vals.append((z['pred'].copy(),z['jack'].copy(),z['truth'].copy(),z['truth_jack'].copy(),r));ctrl=z['control'].copy();cj=z['control_jack'].copy()
    pred,pj,truth,tj=(np.array([v[k] for v in vals]) for k in range(4));targetnames=[v[4]['target'] for v in vals];assert len(targetnames)==len(set(targetnames))>=2
    p=panel_metrics(pred,truth,ctrl,pj,tj,cj,genes,targetnames,targets)
    for k in ('fidelity','reach','jaccard','nmae'):
     vv=[v[4][k] for v in vals if v[4].get(k) is not None and np.isfinite(v[4][k])];p[k]=float(np.mean(vv)) if vv else None
    panels.append(dict(family=source,study=vals[0][4]['study'],seed_offset=offset,arm=arm,unique_targets=len(vals),**p));choices.append(dict(family=source,seed_offset=offset,arm=arm,case_indices=[v[4]['case_index'] for v in vals],target_names=targetnames))
 summary={arm:{k:aggregate(panels,arm,k) for k in METRICS} for arm in ('v24',candidate,candidate+'_shuffled')};write(a/(stage+'_unique_panels.json'),dict(panels=panels,summary=summary,choices=choices,source_inventory=inventory(),auditor_sha256=sha(__file__),scope='Canonical first case_index per literal target within each context/seed panel; no metric-based filtering. Preserves all original count records and original failed gates. Corrected PDS has unique candidate identities, matching competition semantics.'))
 print('Unique-label panel summary',version,stage,json.dumps(summary,indent=2),flush=True)
 return summary

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--version',type=int,choices=(63,64),default=63);p.add_argument('--stage',default='additional');args=p.parse_args();score(args.version,args.stage)
