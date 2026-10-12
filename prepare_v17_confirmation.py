"""Freeze eight power-qualified verification identities, disjoint from V17 development."""
import json
import numpy as np,pandas as pd
import build_submission as b
from v17_methods import ROOT,A

def main():
 path=A/'confirmation_plan.json'
 if path.exists():raise FileExistsError('Preserve verification plan')
 assert not (A/'selection.json').exists() and not (A/'confirmation_evaluation.json').exists()
 plan=json.loads((A/'plan.json').read_text());parentsha=b.sha256(A/'plan.json');records=[]
 for fi,c in enumerate(plan['contexts']):
  used={r['target'] for r in c['cells'] if r['split']=='development'};family=c['family'];m=pd.read_csv(ROOT/('reference/v5_h1' if family=='H1' else 'reference/v7_sources/'+family)/'selected_cells.csv');labels=m['target_gene' if family=='H1' else 'gene_target'].to_numpy(str);counts=pd.Series(labels).value_counts();n=4 if family=='HEK293T' else 2;eligible=sorted(set(counts[counts>=400].index)-used-{'Non-Targeting','non-targeting'});chosen=eligible[:n];assert len(chosen)==n and not set(chosen)&used
  c['cells']=[r for r in c['cells'] if r['split']=='development']
  for seed in plan['seeds']:
   rng=np.random.default_rng(np.random.SeedSequence([171112,fi,seed]))
   for target in chosen:c['cells'].append({'target':target,'seed':seed,'actual':rng.choice(np.flatnonzero(labels==target),400,replace=False).tolist(),'template':rng.choice(c['left'],400,replace=False).tolist(),'reference':c['right'],'split':'confirmation'})
  records.append({'family':family,'verification_targets':chosen,'development_targets':sorted(used)})
 plan['parent_plan_sha256']=parentsha;plan['verification_identity_amendment']='Original draft verification accidentally repeated three current-development target/family pairs. Replaced before candidate selection or verification by alphabetical first eligible >=400-cell targets disjoint from current development. HCT116 only CAT/MYRF qualify and have been inspected in prior versions; no blind biological independence claimed. All development inputs, scores, candidates and policies unchanged.'
 for n in ['prepare_v17_confirmation.py','run_v17_confirmation.py']:plan['dependencies'][n]=b.sha256(ROOT/n)
 b.write_json(path,plan);b.write_json(A/'confirmation_identity_audit.json',{'passed':True,'targets':records,'all_disjoint_from_current_development':True,'original_plan_unchanged_sha256':parentsha,'verification_plan_sha256':b.sha256(path),'before_candidate_selection':True,'blind_new_biology':False});print('V17 verification identities frozen',records,flush=True)
if __name__=='__main__':main()
