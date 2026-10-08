"""Post-selection sensitivity to full independent NTC reference populations.

Diagnostic only: no tuning or selection from these previously confirmed targets.
"""
import json,gc
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from v10_methods import ROOT,A,SEEDS,response_bank,emit
from v8_methods import gated,panel_metrics
from v6_methods import moments,de,direction_metrics
from fit_v7_response import read_cache,raw_inputs
from audit_transfer_validation import profile_jackknife_chunked
from evaluate_atlas_transfer import load_atlas_model
import build_submission as build

def main():
 plan=json.loads((A/'validation_plan.json').read_text());selection=json.loads((A/'count_selection.json').read_text())
 data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model()
 settings={'v4_legacy':('v4','legacy'),selection['candidate']:(selection['response'],selection['generator'])}
 records=[];panels=[]
 for family in ['H1','HCT116','HEK293T']:
  bank,_=response_bank(data,excluded=family)
  x,meta,measured,source=raw_inputs(family,genes)
  for context in [c for c in plan['contexts'] if c['family']==family]:
   population=x[context['left']];depth,p,b,_=moments(population)
   prepared=(depth,p,b,np.asarray(population.mean(axis=0)).ravel());query=p*1e6
   reference=x[context['right']];cp,cj=profile_jackknife_chunked(reference)
   cells=[c for c in context['cells'] if c['split']=='confirmation' and c['seed']==SEEDS[0]]
   actual_profiles=[];actual_jacks=[];pred_profiles={n:[] for n in settings};pred_jacks={n:[] for n in settings}
   for cell in cells:
    target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(np.flatnonzero(genes==target)[0])
    actual=x[cell['actual']];truth=de(actual,reference,genes);tp,tj=profile_jackknife_chunked(actual)
    actual_profiles.append(tp);actual_jacks.append(tj)
    templates=x[cell['template']]
    for name,(response,generator) in settings.items():
     effect=gated(bank[response][ti].copy(),query,own)
     count,diag=emit(templates,np.arange(400),effect,own,SEEDS[0]+ti,generator,prepared=prepared,atlas=atlas)
     metrics=direction_metrics(de(count,reference,genes),truth,genes,own)
     pp,pj=profile_jackknife_chunked(count);pred_profiles[name].append(pp);pred_jacks[name].append(pj)
     records.append({'family':family,'context':context['context'],'target':target,'seed':SEEDS[0],
      'candidate':name,'reference_cells':len(context['right']),'actual_cells':len(cell['actual']),**metrics})
     del count
    print('Full reference checked',context['context'],target,len(context['right']),flush=True)
   for name in settings:
    score=panel_metrics(np.array(pred_profiles[name]),np.array(actual_profiles),cp,np.array(pred_jacks[name]),np.array(actual_jacks),cj,genes,[c['target'] for c in cells],targets)
    subset=pd.DataFrame([r for r in records if r['candidate']==name and r['context']==context['context']])
    panels.append({'family':family,'context':context['context'],'candidate':name,'reference_cells':len(context['right']),**subset.mean(numeric_only=True).to_dict(),**score})
   pd.DataFrame(records).to_csv(A/'full_reference_scores.csv',index=False)
   pd.DataFrame(panels).to_csv(A/'full_reference_panels.csv',index=False)
  del x;gc.collect()
 frame=pd.DataFrame(records)
 summary=frame.groupby(['candidate','family','context'])[['fidelity','reach','jaccard','nmae']].mean().groupby(['candidate','family']).mean().groupby('candidate').mean()
 summary.reset_index().to_csv(A/'full_reference_summary.csv',index=False)
 build.write_json(A/'full_reference_manifest.json',{'passed_computation':True,'candidate':selection['candidate'],
  'selection_sha256':build.sha256(A/'count_selection.json'),'plan_sha256':build.sha256(A/'validation_plan.json'),
  'methods_sha256':build.sha256(ROOT/'v10_methods.py'),'script_sha256':build.sha256(ROOT/'check_v10_reference_size.py'),
  'scores_sha256':build.sha256(A/'full_reference_scores.csv'),'summary_sha256':build.sha256(A/'full_reference_summary.csv'),
  'scope':'Post-selection sensitivity, same confirmed target/cell draws, one seed and full independent right-side controls. Not new blind evidence; no tuning/reselection. Source actual groups remain smaller than the challenge. Reference sizes 768-2144, still below challenge controls.'})
 print(summary.to_string(),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
