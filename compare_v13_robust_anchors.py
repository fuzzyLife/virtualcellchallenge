"""Scale development results by official public source-panel anchors, never private anchors."""
import json
from pathlib import Path
import numpy as np,pandas as pd
import build_submission as build
from fit_v7_response import read_cache,raw_inputs
from audit_transfer_validation import profile_jackknife_chunked,pds
ROOT=build.ROOT;A=ROOT/'experiments/v13_response_programs/analysis/robust_stage'
MAP={'pds':'pds_cosine','fidelity':'de_wilcoxon_direction_fidelity_yield_raw','reach':'de_wilcoxon_direction_reach_raw','nmae':'de_wilcoxon_lfc_nmae','jaccard':'de_wilcoxon_sig_jaccard','common_expression_ratio':'expr_mse_unbiased_capped_norm'}
def main():
 plan=json.loads((A/'validation_plan.json').read_text());manifest=json.loads((A/'official_anchors/manifest.json').read_text());assert manifest['completed']
 for panel in manifest['panels']:
  for file,h in panel['files_sha256'].items():assert build.sha256(A/'official_anchors'/file)==h
 data=read_cache();genes=data['genes'];targets=data['targets'];frame=pd.read_csv(A/'development_panels.csv');frame=frame[(frame.seed==plan['seeds'][0])&frame.reference_mode.eq('full_independent')].copy();rows=[]
 for family in ['H1','HCT116','HEK293T']:
  x,meta,measured,source=raw_inputs(family,genes)
  for ctx in [c for c in plan['contexts'] if c['family']==family]:
   cells=[c for c in ctx['cells'] if c['split']=='development' and c['seed']==plan['seeds'][0]];name=ctx['context'];names=[c['target'] for c in cells]
   baseline=pd.read_csv(A/'official_anchors'/f'{name}_baseline.csv').set_index('statistic').loc['mean'];rep=pd.read_csv(A/'official_anchors'/f'{name}_replicate.csv').set_index('metric').replicate
   cp=profile_jackknife_chunked(x[ctx['right']])[0];mask=~np.isin(genes,names);arrays=[]
   for cell in cells:
    ti=int(np.flatnonzero(targets==cell['target'])[0]);file=A/'count_cache'/f'{name}_{ti:03d}_{cell["seed"]}_development.npz'
    with np.load(file) as z:arrays.append({k:z[k].copy() for k in z.files})
   truth=np.array([z['truth'] for z in arrays])
   for _,row in frame[frame.context.eq(name)].iterrows():
    candidate=row.candidate;pred=np.array([z[candidate+'_pred'] for z in arrays]);row=row.copy();row['pds']=pds(pred[:,mask],truth[:,mask],cp[mask])
    record={'family':family,'context':name,'candidate':candidate}
    for local,official in MAP.items():
     b=float(baseline[official]);r=float(rep[official]);assert np.isfinite(b+r) and abs(r-b)>1e-10
     scaled=(float(row[local])-b)/(r-b)
     if local=='common_expression_ratio':scaled=np.clip(scaled,0,1)
     if local=='nmae':scaled=max(scaled,-6)
     record[local]=float(scaled);record[local+'_raw']=float(row[local]);record[local+'_baseline']=b;record[local+'_replicate']=r
    record['overall']=float(np.mean([record[k] for k in MAP]));rows.append(record)
 df=pd.DataFrame(rows);df.to_csv(A/'development_source_anchored_panels.csv',index=False)
 summary=df.groupby(['candidate','family'])[list(MAP)+['overall']].mean().groupby('candidate').mean();summary.to_csv(A/'development_source_anchored_summary.csv')
 build.write_json(A/'source_anchor_comparison.json',{'completed':True,'plan_sha256':build.sha256(A/'validation_plan.json'),'scores_sha256':build.sha256(A/'development_scores.csv'),'official_anchor_manifest_sha256':build.sha256(A/'official_anchors/manifest.json'),'script_sha256':build.sha256(Path(__file__)),'scope':'First seed, full independent controls, identical actual truth panels and panel-specific target exclusion for PDS. Source anchors only: not a hosted leaderboard forecast; small anchor denominators can amplify noise.'})
 print(summary.to_string(),flush=True)
if __name__=='__main__':main()
