"""Prospectively freeze source tests and amend development power without confirmation outcomes."""
import json,hashlib,shutil
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
import build_submission as build
from fit_v7_response import read_cache,raw_inputs
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis';S=A/'flex_stage'
def main():
 S.mkdir(exist_ok=True);path=S/'validation_plan.json';assert not path.exists()
 plan=json.loads((A/'validation_plan.json').read_text());amend=json.loads((A/'confirmation_power_amendment.json').read_text())
 for extra in amend['additional_cells']:
  ctx=next(c for c in plan['contexts'] if c['context']==extra['context']);ctx['cells'].append({k:v for k,v in extra.items() if k not in ['context','family']})
  ctx['confirmation_targets']=sorted(set(ctx['confirmation_targets'])|{extra['target']})
 meta=pd.read_csv(ROOT/'reference/v7_sources/HEK293T/selected_cells.csv');meta['sample']=meta['sample'].astype(str)
 ctx=next(c for c in plan['contexts'] if c['context']=='HEK293T_controlstratum1');target='RNF10';targets=pd.read_csv(ROOT/'pert_counts.csv').iloc[:,0].to_numpy(str);ti=int(np.flatnonzero(targets==target)[0])
 assert all(c['target']!=target for c in ctx['cells'])
 ids=np.flatnonzero(meta.gene_target.eq(target).to_numpy()&meta['sample'].isin(ctx['batches']).to_numpy());assert len(ids)==120
 for seed in plan['seeds']:
  common=next(c for c in ctx['cells'] if c['seed']==seed)
  actual=np.random.default_rng(seed+ti+71000).choice(ids,len(ids),replace=False).tolist()
  ctx['cells'].append({'target':target,'split':'development','seed':seed,'actual':actual,'template':common['template'],'reference':common['reference']})
 ctx['development_targets'].append(target)
 plan.update(emitters=['v11','flex_combine','flex_prefer','flex_overlap'],
  source_test='Within-K562 equal-platform blend, prefer-Flex with legacy fallback, and overlap-only blend; V11 count emitter fixed for all.',
  source_negative_control='Target-shuffled Flex effects on source-summary screens; K562-family exclusion must erase every Flex contribution.',
  source_size_limitation='Legacy K562 is pseudobulk, so an actual size-matched cell downsample cannot be reconstructed. Report source counts and measured dispersion; do not claim that control was run.',
  development_power_amendment='RNF10, already development in HEK stratum0, added with all 120 genuine cells in stratum1. Chosen from existing summary power and metadata because one sampled-reference panel otherwise lacks NMAE; no candidate-specific outcome used.',
  confirmation_power_amendment_sha256=build.sha256(A/'confirmation_power_amendment.json'),
  parent_count_plan_sha256=build.sha256(A/'validation_plan.json'),source_methods_sha256=build.sha256(ROOT/'v12_flex_methods.py'),cpu_metrics_sha256=build.sha256(ROOT/'v12_cpu_metrics.py'),
  script_sha256=build.sha256(Path(__file__)))
 dev={c['target'] for ctx in plan['contexts'] for c in ctx['cells'] if c['split']=='development'};confirm={c['target'] for ctx in plan['contexts'] for c in ctx['cells'] if c['split']=='confirmation'};assert not dev&confirm
 build.write_json(path,plan)
 # Four unchanged official development panels and anchors can be reused byte-for-byte.
 previous=json.loads((A/'official_panel_manifest.json').read_text());records=[];D=S/'official_panels';D.mkdir(exist_ok=True);O=S/'official_anchors';O.mkdir(exist_ok=True)
 data=read_cache();genes=data['genes']
 for r in previous['panels']:
  if r['context']!='HEK293T_controlstratum1':
   records.append(r)
   for p in (A/'official_anchors').glob(r['context']+'*'):shutil.copyfile(p,O/p.name)
   continue
  x,m,measured,source=raw_inputs('HEK293T',genes);ctx=next(c for c in plan['contexts'] if c['context']==r['context'])
  cells=[c for c in ctx['cells'] if c['split']=='development' and c['seed']==plan['seeds'][0]]
  rows=ctx['right']+sum([c['actual'] for c in cells],[]);labels=['non-targeting']*len(ctx['right'])+sum([[c['target']]*len(c['actual']) for c in cells],[])
  file=D/(r['context']+'_real.h5ad');ad.AnnData(x[rows],obs=pd.DataFrame({'target':pd.Categorical(labels)},index=[str(i) for i in rows]),var=pd.DataFrame(index=genes)).write_h5ad(file,compression='gzip')
  records.append({**r,'file':str(file.relative_to(ROOT)),'sha256':build.sha256(file),'cells':len(rows),'targets':[c['target'] for c in cells]})
 build.write_json(S/'official_panel_manifest.json',{'panels':records,'plan_sha256':build.sha256(path),'scope':'Development only; official source anchors are not private competition anchors.'})
 print('Matched-Flex test plan frozen; no confirmation counts read',flush=True)
if __name__=='__main__':main()
