"""400-distinct-treated-cell power guard on metadata-eligible H1 V11 cases."""
import argparse,json
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from v11_methods import ROOT,A,SEEDS,SETTINGS,response_bank,response_for,emit
from fit_v7_response import read_cache,raw_inputs
from v6_methods import moments,de,direction_metrics
import build_submission as build


def main(stage):
 output=A/f'{stage}_power_validation.json'
 if output.exists():raise FileExistsError('Power diagnostic already frozen')
 plan=json.loads((A/'validation_plan.json').read_text());data=read_cache();genes=data['genes'];targets=data['targets']
 settings=SETTINGS if stage=='development' else {name:SETTINGS[name] for name in ['v10',json.loads((A/'count_selection.json').read_text())['candidate']]}
 bank,training=response_bank(data,excluded='H1');assert 'H1' not in training['training_families']
 x,meta,measured,source=raw_inputs('H1',genes);records=[];provenance=[]
 for context in [c for c in plan['contexts'] if c['family']=='H1']:
  _,p,b,_=moments(x[context['left']]);query=p*1e6
  for seed in SEEDS:
   cells=[c for c in context['cells'] if c['seed']==seed and c['split']==stage]
   template=x[cells[0]['template']];depth=np.asarray(template.sum(axis=1)).ravel();prepared=(depth,p,b)
   refs={'sampled400':x[cells[0]['reference']],'full_independent':x[context['right']]}
   for cell in cells:
    available=np.flatnonzero(meta.gene_target.eq(cell['target']).to_numpy()&meta['sample'].isin(context['batches']).to_numpy())
    if len(available)<400:continue
    target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(np.flatnonzero(genes==target)[0])
    actual_rows=np.random.default_rng(seed+ti+44000).choice(available,400,replace=False);actual=x[actual_rows]
    truths={mode:de(actual,ref,genes) for mode,ref in refs.items()}
    provenance.append({'context':context['context'],'target':target,'seed':seed,'actual':actual_rows.tolist(),
                       'template':cell['template'],'reference':cell['reference'],'full_reference':context['right']})
    for name,setting in settings.items():
     effect=response_for(bank,setting,query,own,ti);count,_=emit(template,np.arange(400),effect,own,seed+ti,setting,prepared=prepared)
     for mode,ref in refs.items():records.append({'candidate':name,'target':target,'context':context['context'],'seed':seed,'reference_mode':mode,
      'actual_cells':400,'reference_cells':ref.shape[0],**direction_metrics(de(count,ref,genes),truths[mode],genes,own)})
     print('V11',stage,'power400',target,seed,name,'evaluated',flush=True)
 df=pd.DataFrame(records);assert len(df)>0
 df.to_csv(A/f'{stage}_power_scores.csv',index=False)
 summary=df.groupby(['reference_mode','candidate'])[['fidelity','reach','jaccard','nmae','predicted_calls','reference_calls']].mean().reset_index()
 summary.to_csv(A/f'{stage}_power_summary.csv',index=False);guards={name:True for name in settings}
 for mode in ['sampled400','full_independent']:
  values=summary[summary.reference_mode==mode].set_index('candidate');baseline=values.loc['v10']
  for name in settings:
   row=values.loc[name];guards[name]&=bool(row.fidelity>=baseline.fidelity-.02 and row.nmae<=baseline.nmae+.02)
 build.write_json(output,{'completed':True,'guards':guards,'scores_sha256':build.sha256(A/f'{stage}_power_scores.csv'),
  'summary_sha256':build.sha256(A/f'{stage}_power_summary.csv'),'methods_sha256':build.sha256(ROOT/'v11_methods.py'),
  'script_sha256':build.sha256(ROOT/'audit_v11_cell_power.py'),'plan_sha256':build.sha256(A/'validation_plan.json'),'provenance':provenance,
  'policy':'Before model selection, metadata availability justified a guard on 400 unique H1 treated cells. Fidelity within .02 and NMAE within .02 on both control-reference sizes. Development ANKZF1; confirmation SLIRP. HCT/HEK lack 400 available cells and are not duplicated. H1 remains reused development biology, not independent lineage evidence.'})
 print(summary.to_string(index=False),'Power guards:',guards,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['development','confirmation'],required=True);args=p.parse_args()
 with threadpool_limits(limits=2):main(args.stage)
