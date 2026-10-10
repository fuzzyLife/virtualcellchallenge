"""Export development truth only for official public-source anchors."""
import json,gc
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from threadpoolctl import threadpool_limits
import build_submission as build
from fit_v7_response import read_cache,raw_inputs
ROOT=build.ROOT;A=ROOT/'experiments/v12_distribution_calibrated/analysis';D=A/'official_panels'
def main():
 D.mkdir(exist_ok=True);plan=json.loads((A/'validation_plan.json').read_text());data=read_cache();genes=data['genes'];records=[]
 for family in ['H1','HCT116','HEK293T']:
  x,meta,measured,source=raw_inputs(family,genes)
  for ctx in [c for c in plan['contexts'] if c['family']==family]:
   cells=[c for c in ctx['cells'] if c['split']=='development' and c['seed']==plan['seeds'][0]]
   rows=ctx['right']+sum([c['actual'] for c in cells],[])
   labels=['non-targeting']*len(ctx['right'])+sum([[c['target']]*len(c['actual']) for c in cells],[])
   file=D/(ctx['context']+'_real.h5ad')
   if not file.exists():
    a=ad.AnnData(x[rows],obs=pd.DataFrame({'target':pd.Categorical(labels)},index=[str(i) for i in rows]),var=pd.DataFrame(index=genes));a.write_h5ad(file,compression='gzip')
   records.append({'family':family,'context':ctx['context'],'file':str(file.relative_to(ROOT)),'sha256':build.sha256(file),'cells':len(rows),'targets':[c['target'] for c in cells]})
   print('Official source panel exported',ctx['context'],len(rows),flush=True)
  del x;gc.collect()
 build.write_json(A/'official_panel_manifest.json',{'panels':records,'plan_sha256':build.sha256(A/'validation_plan.json'),'script_sha256':build.sha256(Path(__file__)),
  'scope':'Development data only. Official source-panel PDS excludes actual subset target genes; screening PDS excludes all 300 challenge targets. Panels and anchors are not competition ABC anchors.'})
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
