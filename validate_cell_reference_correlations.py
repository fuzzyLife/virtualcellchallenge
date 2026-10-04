"""Independently validate complete output axes, summaries and real-cell spot checks."""
from pathlib import Path
import hashlib,json
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis/cell_reference_correlations'
NAMES=['K562','HCT116','HEK293T','H1','CD4']
def main():
 frame=pd.read_csv(OUT/'per_cell_correlations.csv');summary=pd.read_csv(OUT/'context_reference_summary.csv')
 report=json.loads((OUT/'analysis_report.json').read_text());genes=pd.read_csv(ROOT/'gene_names.csv').gene_name.to_numpy()
 profile=pd.read_csv(OUT/'reference_control_profiles_common_genes.csv');positions={g:i for i,g in enumerate(genes)}
 axis=np.array([positions[g] for g in profile.gene_name]);refs=profile[NAMES].to_numpy().T
 refs=np.log1p(refs/refs.sum(axis=1,keepdims=True)*1e6);checks=[]
 for context in ['A','B','C']:
  f=ad.read_h5ad(ROOT/f'context_{context}.h5ad',backed='r');sub=frame[frame.context==context]
  assert len(sub)==18400 and np.array_equal(sub.cell_id.to_numpy(),f.obs_names.to_numpy())
  assert np.array_equal(sub.ntc_id.to_numpy(),f.obs.ntc_id.astype(str).to_numpy())
  columns=['pearson_log1p_cpm_'+r for r in NAMES];values=sub[columns].to_numpy()
  assert np.array_equal(sub.top_reference_pearson.to_numpy(),np.array(NAMES)[np.argmax(values,axis=1)])
  for i,name in enumerate(NAMES):
   row=summary[(summary.context==context)&(summary.reference==name)].iloc[0]
   assert abs(np.median(values[:,i])-row.median_pearson)<1e-12
   assert int((sub.top_reference_pearson==name).sum())==row.top_reference_cells
  for index in [0,9199,18399]:
   raw=f.X[index:index+1].toarray().ravel().astype(np.float64);selected=raw[axis]
   log=np.log1p(selected/selected.sum()*1e6);record=sub.iloc[index]
   assert record.total_counts==raw.sum() and record.detected_genes==np.count_nonzero(raw)
   for j,name in enumerate(NAMES):
    pearson=float(np.corrcoef(log,refs[j])[0,1]);spearman=float(spearmanr(log,refs[j]).statistic)
    assert abs(pearson-record['pearson_log1p_cpm_'+name])<1e-12
    assert abs(spearman-record['spearman_'+name])<1e-12
    checks.append({'context':context,'cell_index_0based':index,'reference':name,'pearson':pearson,'spearman':spearman})
  f.file.close()
 assert len(frame)==55200 and not frame[['context','cell_id']].duplicated().any()
 numeric=frame.select_dtypes('number').to_numpy();assert np.isfinite(numeric).all()
 for name,record in report['outputs'].items():
  assert hashlib.sha256((OUT/name).read_bytes()).hexdigest()==record['sha256'],name
 result={'passed':True,'complete_rows':len(frame),'complete_cell_ids_guides_and_contexts_match':True,'all_numeric_values_finite':True,'all_output_csv_sha256_match':True,'all_summary_medians_and_winner_counts_match':True,'independent_original_count_spot_checks':checks,'real_cell_reference_pairs_checked':len(checks),'validator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (OUT/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
 print('Verified every output row/identifier/summary and 45 real-cell correlation pairs against independent NumPy/SciPy calculations.')
if __name__=='__main__':main()
