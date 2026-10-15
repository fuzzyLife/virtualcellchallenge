"""Post-rejection mean-response gap on spent public H1 treatment cells only.

The retrospective scaling coefficient uses actual treatment outcomes and is
explicitly an oracle diagnostic: never an inference rule or submission fit.
"""
from pathlib import Path
import json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
from v29_methods import ROOT,A,load_plan,write
from v26_methods import raw
from fit_v7_response import read_cache
from v12_cpu_metrics import de
import build_submission as b

def main():
 assert not json.loads((A/'release_decision.json').read_text())['released'];plan=load_plan();ctx=next(c for c in plan['contexts'] if c['family']=='H1');data=read_cache();genes=data['genes'];targets=data['targets'];x,m,measured,_=raw('H1',genes);control=x[ctx['right']];rows=[]
 for cell in ctx['cells']:
  own=int(np.flatnonzero(genes==cell['target'])[0]);truth=de(x[cell['actual']],control,genes);file=A/'parallel/treatment_H1/H1'/f'{cell["target"]}_{cell["seed"]}_v24_counts.npz';pred=de(sparse.load_npz(file),control,genes)
  for scope,extra in [('official_own_excluded',np.arange(len(genes))!=own),('nonpanel',~np.isin(genes,targets))]:
   mask=truth['eligible']&(truth['q']<.05)&extra&np.isfinite(truth['lfc'])&(abs(truth['lfc'])>1e-8);tt=truth['lfc'][mask];pp=pred['lfc'][mask];den=abs(tt).sum();slope=float(pp@tt/max(pp@pp,1e-30));rows.append({'target':cell['target'],'seed':cell['seed'],'scope':scope,'true_DE_genes':int(mask.sum()),'median_abs_true_log2FC':float(np.median(abs(tt))),'median_abs_predicted_log2FC':float(np.median(abs(pp))),'predicted_over_true_absolute_effect_sum':float(abs(pp).sum()/max(den,1e-30)),'sign_agreement_on_true_DE':float(np.mean(np.sign(pp)==np.sign(tt))),'raw_NMAE':float(abs(pp-tt).sum()/max(den,1e-30)),'oracle_L2_amplitude_scale':slope,'oracle_scaled_raw_NMAE':float(abs(slope*pp-tt).sum()/max(den,1e-30)),'prediction_file_sha256':b.sha256(file)})
 pd.DataFrame(rows).to_csv(A/'mean_response_gap.csv',index=False);write(A/'mean_response_gap.json',{'plan_sha256':b.sha256(A/'plan.json'),'release_decision_sha256':b.sha256(A/'release_decision.json'),'table_sha256':b.sha256(A/'mean_response_gap.csv'),'diagnostician_sha256':b.sha256(Path(__file__)),'scope':'Retrospective descriptive/oracle analysis of already-scored H1 targets, not a trained model or prospective validation. Oracle slopes use actual treatment values, cannot be used in query inference and are not deployed. No guard changes, reserve opening or candidate reselection.'});print(pd.DataFrame(rows).to_string(index=False),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
