"""Prospective V12 amendment: corrected metrics, powered development, locked confirmation."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import build_submission as build
from fit_v7_response import read_cache
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis'
SEEDS=[2026100911,2026100927]
EMITTERS=['v11','gamma100','gamma25','empirical100']
def main():
 A.mkdir(parents=True,exist_ok=True)
 path=A/'validation_plan.json'
 if path.exists():raise FileExistsError(path)
 old=json.loads((ROOT/'experiments/v11_consensus_coverage/analysis/validation_plan.json').read_text())
 future=json.loads((ROOT/'analysis/post_v11_research/future_experiment_plan.json').read_text())
 reserved={g for c in future['new_count_confirmation'] for g in c['targets']}
 used=set(future['excluded_previous_count_targets']) if 'excluded_previous_count_targets' in future else set(future.get('previously_count_tested_targets',[]))
 # The exact earlier exclusion key is saved independently in the original plan.
 for c in old['contexts']:used.update(z['target'] for z in c['cells'])
 data=read_cache();genes=data['genes'];targets=data['targets'];contexts=[]
 for ci,ctx in enumerate(old['contexts']):
  family=ctx['family'];directory=ROOT/('reference/v5_h1' if family=='H1' else f'reference/v7_sources/{family}')
  meta=pd.read_csv(directory/'selected_cells.csv')
  if family=='H1':meta['sample']=meta.batch.astype(str);meta['gene_target']=meta.target_gene
  else:meta['sample']=meta['sample'].astype(str)
  counts=meta.loc[meta['sample'].isin(ctx['batches']),'gene_target'].value_counts()
  with np.load(ROOT/f'experiments/v7_multisource_adaptive/analysis/{family}_summaries.npz') as z:
   valid=z['2_valid'];q=z['2_q'];control=z['2_control']
   power=((q<.05)&valid&(control>5)).sum(axis=1)
  available=[str(t) for t in targets if counts.get(t,0)>= (100 if family=='H1' else 150)]
  strong=sorted([t for t in available if t not in reserved],key=lambda t:(-power[np.flatnonzero(targets==t)[0]],t))[:4]
  confirmation=next((c['targets'] for c in future['new_count_confirmation'] if c['context']==ctx['context']),[])
  # Additional confirmation is chosen only from already viewed summary power, never count predictions.
  fresh=[t for t in available if t not in used|reserved and power[np.flatnonzero(targets==t)[0]]>=50]
  fresh.sort(key=lambda t:hashlib.sha256(f'V12-power-confirm:{t}'.encode()).hexdigest())
  extra=fresh[:2] if family!='H1' else []
  confirmation=confirmation+extra;strong=[t for t in strong if t not in extra]
  development=sorted(set(strong)|set(ctx['development_targets'][:2]))
  development=[t for t in development if t not in reserved|set(confirmation)]
  cells=[]
  for seed in SEEDS:
   rng=np.random.default_rng(seed+101*ci)
   template=rng.choice(ctx['left'],400,replace=False).tolist();reference=rng.choice(ctx['right'],400,replace=False).tolist()
   for split,names in [('development',development),('confirmation',confirmation)]:
    for t in names:
     ids=np.flatnonzero(meta.gene_target.eq(t).to_numpy()&meta['sample'].isin(ctx['batches']).to_numpy());ti=int(np.flatnonzero(targets==t)[0])
     actual=np.random.default_rng(seed+ti+71000).choice(ids,min(400,len(ids)),replace=False).tolist()
     cells.append({'target':t,'split':split,'seed':seed,'actual':actual,'template':template,'reference':reference})
  contexts.append({**{k:v for k,v in ctx.items() if k not in ['cells','development_targets','confirmation_targets']},
   'development_targets':development,'confirmation_targets':confirmation,'additional_powered_confirmation':extra,'cells':cells,
   'summary_power':{t:int(power[np.flatnonzero(targets==t)[0]]) for t in development+confirmation}})
  print(ctx['context'],'development',development,'confirmation',confirmation,flush=True)
 plan={'version':12,'contexts':contexts,'seeds':SEEDS,'emitters':EMITTERS,
  'source_response':'Frozen V11 consensus20, whole source family excluded; response changes screened separately.',
  'metric_fixes':['Wilcoxon no continuity correction, tie correction','NMAE minimum 10 real-significant genes and finite nonzero denominator','reach=0 for nonempty real budget with no adjudicable genes'],
  'selection_rules':future['selection_rules'],'guard_tolerances':{'fidelity':.01,'reach':.01,'family_fidelity':.03,'family_reach':.03,'nmae':.01,'pds':.02},
  'expression_improve_required_both_stages':True,'confirmation_policy':'Choose using development only. Evaluate baseline and frozen candidate once on reserved count confirmation; do not retune after failure.',
  'power_amendment':'Add strong development cases chosen from previously inspected full-source normal-q summaries and cell counts. Additional confirmation uses the same metadata only; summary biology is not blind. Keep all original reserved confirmation targets, including low-power ones. H1 is development only.',
  'scope':'Public whole-family held-out tests, two new seeds, actual <=400 cells without duplication, shared 400-cell template, independent 400/full controls. Private competition anchors unavailable.',
  'parent_plan_sha256':build.sha256(ROOT/'analysis/post_v11_research/future_experiment_plan.json'),
  'source_cache_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),
  'script_sha256':build.sha256(Path(__file__))}
 build.write_json(path,plan)
if __name__=='__main__':main()
