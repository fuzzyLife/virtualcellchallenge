"""Correct-depth, whole-lineage-held V11 generated-count validation with shared origins."""
import argparse,json,gc
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from v11_methods import ROOT,SETTINGS,response_bank,response_for
from v12_methods import emit
from v12_flex_methods import banks as flex_banks
from v12_cpu_metrics import de,direction_metrics
A=ROOT/'experiments/v12_distribution_calibrated/analysis/flex_stage'
from fit_v7_response import read_cache,raw_inputs
from v6_methods import moments
from v8_methods import panel_metrics
from audit_transfer_validation import profile_jackknife_chunked,pds
from evaluate_atlas_transfer import load_atlas_model
import build_submission as build
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']


def balanced(frame):
 return frame.groupby(['reference_mode','candidate','family','context'])[FIELDS].mean().groupby(['reference_mode','candidate','family']).mean().groupby(['reference_mode','candidate']).mean().reset_index()


def main(stage,workers):
 plan=json.loads((A/'validation_plan.json').read_text())
 assert build.sha256(ROOT/'v12_flex_methods.py')==plan['source_methods_sha256']
 assert build.sha256(ROOT/'v12_cpu_metrics.py')==plan['cpu_metrics_sha256']
 extraction=json.loads((A.parent/'flex_extraction.json').read_text())
 assert extraction['all_300_targets_present'] and extraction['unique_controls']>=1000
 assert extraction['control_identity_correction_sha256']==build.sha256(A.parent/'flex_control_identity_correction.json')
 assert extraction['qc_plan_sha256']==build.sha256(A.parent/'flex_source_qc_plan.json')
 assert build.sha256(ROOT/'reference/v12_resources/flex_source/summaries.npz')==extraction['summary_sha256']
 if stage=='confirmation':
  frozen=json.loads((A/'count_selection.json').read_text());assert not frozen['confirmation_read']
  assert frozen['development_scores_sha256']==build.sha256(A/'development_scores.csv')
  names=['v11',frozen['candidate']]
 else:names=plan['emitters']
 settings={name:SETTINGS['consensus20_noise65'] for name in names}
 SEEDS=plan['seeds']
 data=read_cache();genes=data['genes'];targets=data['targets'];atlas,_=load_atlas_model()
 inputs={'plan':build.sha256(A/'validation_plan.json'),'methods':build.sha256(ROOT/'v12_flex_methods.py'),'emitter':build.sha256(ROOT/'v12_methods.py'),'flex_summary':build.sha256(ROOT/'reference/v12_resources/flex_source/summaries.npz'),'metrics':build.sha256(ROOT/'v12_cpu_metrics.py'),
         'script':build.sha256(ROOT/'evaluate_v12_flex.py'),'settings':settings}
 cache=A/'count_cache';cache.mkdir(exist_ok=True);records=[];panels=[]
 for family in ['H1','HCT116','HEK293T']:
  banks,training=flex_banks(data,excluded=family)
  assert family not in training['training_families']
  x,meta,measured,source=raw_inputs(family,genes)
  for context in [c for c in plan['contexts'] if c['family']==family]:
   population=x[context['left']];depth,p,b,var=moments(population);query=p*1e6
   if not any(c['split']==stage for c in context['cells']):continue
   for seed in SEEDS:
    cells=[c for c in context['cells'] if c['split']==stage and c['seed']==seed]
    template=x[cells[0]['template']];libs=np.asarray(template.sum(axis=1)).ravel();prepared=(libs,p,b,var,float(np.mean(1/depth)))
    refs={'sampled400':x[cells[0]['reference']],'full_independent':x[context['right']]}
    refprofiles={label:profile_jackknife_chunked(ref) for label,ref in refs.items()}
    assert all(c['template']==cells[0]['template'] and c['reference']==cells[0]['reference'] for c in cells)
    def evaluate(cell):
     target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(np.flatnonzero(genes==target)[0])
     file=cache/f'{context["context"]}_{ti:03d}_{seed}_{stage}.npz';side=file.with_suffix('.json')
     if file.exists() and side.exists():
      rec=json.loads(side.read_text());assert rec['inputs']==inputs and rec['sha256']==build.sha256(file)
      with np.load(file,allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
      return rec['records'],arrays
     actual=x[cell['actual']];tp,tj=profile_jackknife_chunked(actual)
     truth={label:de(actual,ref,genes) for label,ref in refs.items()};arrays={'truth':tp,'truth_jack':tj};rows=[]
     for name,setting in settings.items():
      effect=response_for(banks[name],setting,query,own,ti)
      count,diag=emit(template,np.arange(400),effect,own,seed+ti,'v11',prepared=prepared,atlas=atlas)
      np.testing.assert_array_equal(np.asarray(count.sum(axis=1)).ravel(),libs)
      pp,pj=profile_jackknife_chunked(count);arrays[name+'_pred']=pp;arrays[name+'_jack']=pj
      for label,ref in refs.items():
       metrics=direction_metrics(de(count,ref,genes),truth[label],genes,own)
       rows.append({'family':family,'context':context['context'],'target':target,'seed':seed,'stage':stage,
        'candidate':name,'reference_mode':label,'reference_cells':ref.shape[0],'actual_cells':actual.shape[0],
        'mean_detected_genes':float(count.getnnz(axis=1).mean()),**metrics,**diag})
      del count
     np.savez_compressed(file,**arrays)
     sanitized=[{k:None if isinstance(v,(float,np.floating)) and not np.isfinite(v) else v for k,v in row.items()} for row in rows]
     build.write_json(side,{'inputs':inputs,'sha256':build.sha256(file),'records':sanitized})
     print('V12',stage,context['context'],seed,target,'evaluated',flush=True)
     return sanitized,arrays
    with ThreadPoolExecutor(max_workers=workers) as executor:evaluated=list(executor.map(evaluate,cells))
    for group,arrays in evaluated:records.extend(group)
    tp=np.array([arrays['truth'] for _,arrays in evaluated]);tj=np.array([arrays['truth_jack'] for _,arrays in evaluated])
    constant=profile_jackknife_chunked(template)[0];mask=~np.isin(genes,targets)
    np.testing.assert_allclose(pds(np.broadcast_to(constant,(len(cells),len(genes)))[:,mask],tp[:,mask],refprofiles['sampled400'][0][mask]),.5,atol=1e-12)
    for label,(cp,cj) in refprofiles.items():
     for name in settings:
      pp=np.array([arrays[name+'_pred'] for _,arrays in evaluated]);pj=np.array([arrays[name+'_jack'] for _,arrays in evaluated])
      score=panel_metrics(pp,tp,cp,pj,tj,cj,genes,[c['target'] for c in cells],targets)
      subset=pd.DataFrame([r for group,_ in evaluated for r in group if r['candidate']==name and r['reference_mode']==label])
      panels.append({'family':family,'context':context['context'],'seed':seed,'stage':stage,'candidate':name,'reference_mode':label,
                    **subset.mean(numeric_only=True).to_dict(),**score})
    pd.DataFrame(records).to_csv(A/f'{stage}_scores.csv',index=False);pd.DataFrame(panels).to_csv(A/f'{stage}_panels.csv',index=False)
    balanced(pd.DataFrame(panels)).to_csv(A/f'{stage}_summary.csv',index=False)
  del x;gc.collect()
 build.write_json(A/f'{stage}_evaluation.json',{'completed':True,'inputs':inputs,'scores_sha256':build.sha256(A/f'{stage}_scores.csv'),
  'panels_sha256':build.sha256(A/f'{stage}_panels.csv'),'summary_sha256':build.sha256(A/f'{stage}_summary.csv'),
  'settings':settings,'template_depths_verified_exactly':True,'shared_reference_chance_discrimination_verified':True,
  'scope':'Corrected-metric unanchored public panels; two new seeds, whole-family response exclusion, 400 and full independent controls. Not private leaderboard anchors.'})
 print(pd.read_csv(A/f'{stage}_summary.csv').to_string(index=False),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['development','confirmation'],required=True);p.add_argument('--workers',type=int,default=4);args=p.parse_args()
 with threadpool_limits(limits=2):main(args.stage,args.workers)
