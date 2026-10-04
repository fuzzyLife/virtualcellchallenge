"""Compare representative v2 and first-submission responses against provided controls.

These are response diagnostics, not hidden-label accuracy estimates. Target
selection deliberately spans effect strength and source coverage.
"""
import json,time
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import build_submission as build
from audit_prediction import ControlRankTest,normalized_counts
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'experiments/v2_cd4_hybrid'
def main():
 manifest,genes,targets=build.read_spec(ROOT)
 coverage=pd.read_csv(OUT/'analysis/reference_target_coverage.csv')
 with np.load(OUT/'analysis/reference_model_arrays.npz',allow_pickle=False) as z:effects=z['effects'][0].copy()
 for ti,target in enumerate(targets):effects[ti,np.flatnonzero(genes==target)[0]]=0.
 energy=np.sqrt(np.mean(effects.astype(float)**2,axis=1))
 covered=np.flatnonzero(coverage.eligible_cd4_conditions.to_numpy()>0)
 ordered=covered[np.argsort(energy[covered])]
 chosen=ordered[np.rint(np.linspace(0,len(ordered)-1,8)).astype(int)].tolist()
 fallback=np.flatnonzero((coverage.eligible_cd4_conditions.to_numpy()==0)&coverage.k562_available.to_numpy())
 chosen.extend(fallback[:2].tolist());chosen.extend(np.flatnonzero(coverage.status.to_numpy()=='target-only fallback')[:2].tolist())
 frames=[];original=ad.read_h5ad(ROOT/'prediction.h5ad',backed='r');candidate=ad.read_h5ad(OUT/'prediction.h5ad',backed='r')
 for ci,context in enumerate(manifest['contexts']):
  control=ad.read_h5ad(ROOT/f'context_{context}.h5ad')
  assert np.array_equal(control.var_names.to_numpy(),genes)
  cpm=normalized_counts(control.X);keep=np.flatnonzero(np.asarray(cpm.mean(axis=0)).ravel()>5)
  rank=ControlRankTest(cpm[:,keep]);print(f'{context}: rank distributions ready',flush=True)
  for ti in chosen:
   offset=(ci*len(targets)+ti)*manifest['cells_per_pert'];end=offset+manifest['cells_per_pert']
   for name,obj in [('first_submission',original),('v2_cd4_hybrid',candidate)]:
    view=obj[offset:end,:];assert set(view.obs.target_gene.astype(str))=={targets[ti]};assert set(view.obs.context.astype(str))=={context}
    matrix=view.X.tocsr();_,q,lfc=rank.test(normalized_counts(matrix)[:,keep]);sig=(q<.05)&(genes[keep]!=targets[ti])
    frames.append({'context':context,'target_gene':targets[ti],'variant':name,'source_status':coverage.status.iloc[ti],
     'significant_downstream_genes':int(sig.sum()),'significant_up':int((sig&(lfc>0)).sum()),'significant_down':int((sig&(lfc<0)).sum()),
     'mean_library':float(matrix.sum()/matrix.shape[0]),'selected_reference_effect_rms':float(energy[ti])})
   print(f'{context}/{targets[ti]}: '+str([r['significant_downstream_genes'] for r in frames[-2:]]),flush=True)
  del control,cpm,rank
 original.file.close();candidate.file.close()
 table=pd.DataFrame(frames);table.to_csv(OUT/'analysis/paired_response_check.csv',index=False)
 summary={'scope':'Provided-control response diagnostics only; not independent performance validation or a VCC score.',
  'selection':'Eight CD4 effect-strength quantiles, two K562 fallbacks, two target-only fallbacks; deliberately nonrandom.',
  'targets':targets[chosen].tolist(),'group_summaries':table.groupby(['context','variant']).significant_downstream_genes.agg(['count','median','mean','max']).reset_index().to_dict('records')}
 build.write_json(OUT/'analysis/paired_response_check.json',summary)
 print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()