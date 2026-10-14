"""Audit legacy CD4 DE-stat units/directions against newly verified raw counts."""
from pathlib import Path
import json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits
import build_submission as b
from fit_v7_response import read_cache,summarize_raw
from v22_methods import raw
ROOT=b.ROOT;A=ROOT/'experiments/v24_cd4_response/analysis'
def main():
 A.mkdir(parents=True,exist_ok=True);d=read_cache();x,m,measured,rec=raw('CD4',d['genes']);s=summarize_raw(x,m,d['genes'],d['targets'],measured,2);np.savez_compressed(A/'D2_Rest_summary.npz',**s,genes=d['genes'],targets=d['targets']);rows=[];nonpanel=~np.isin(d['genes'],d['targets'])
 for ri in np.flatnonzero((d['row_families']=='CD4')&(d['row_contexts']==1)):
  ti=int(d['row_targets'][ri]);mask=d['valid'][ri]&s['valid'][ti]&nonpanel&(s['control'][ti]>5)&(abs(d['truth'][ri])>.03)&(abs(s['truth'][ti])>.03)
  if mask.sum()<100:continue
  a=d['truth'][ri,mask];z=s['truth'][ti,mask];own=int(np.flatnonzero(d['genes']==d['targets'][ti])[0]);rows.append({'target':str(d['targets'][ti]),'raw_cells':int(s['n'][ti]),'genes':int(mask.sum()),'spearman':float(spearmanr(a,z).statistic),'sign_agreement':float(np.mean(np.sign(a)==np.sign(z))),'raw_on_legacy_slope':float(a@z/(a@a)),'raw_own_lfc':float(s['truth'][ti,own]),'legacy_own_lfc':float(d['truth'][ri,own]),'raw_knockdown':float(s['knockdown'][ti]),'legacy_knockdown':float(d['knockdown'][ri])})
 f=pd.DataFrame(rows);f.to_csv(A/'legacy_CD4_raw_agreement.csv',index=False);print(f.to_string(index=False),flush=True);b.write_json(A/'raw_summary_provenance.json',{'source_extraction_sha256':b.sha256(ROOT/'experiments/v22_state_response/analysis/cd4_selected/extraction.json'),'summary_sha256':b.sha256(A/'D2_Rest_summary.npz'),'source_cache_sha256':b.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),'audit_script_sha256':b.sha256(Path(__file__)),'scope':'Spent donor2 raw-cell data, diagnostic of existing source response estimation; no fresh reserve score.'})
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
