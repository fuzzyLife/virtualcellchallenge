"""Summarize all generated V90 moment diagnostics without claiming prediction accuracy."""
from pathlib import Path
import csv,numpy as np
from v57_pie import write,sha
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v90_interval_diagnostic/analysis'
if __name__=='__main__':
 rows=list(csv.DictReader((A/'prediction_group_statistics.csv').open(encoding='utf-8')));assert len(rows)==900;assert len({(r['context'],r['target_gene']) for r in rows})==900;summary={}
 for context in ['A','B','C']:
  cases=[r for r in rows if r['context']==context];assert len(cases)==300;out={'groups':300}
  for key in ['mean_probability_L1','bulk_probability_L1']:
   values=np.array([float(r[key]) for r in cases]);assert np.isfinite(values).all() and (values>=0).all();out[key]=dict(mean=float(values.mean()),maximum=float(values.max()),minimum=float(values.min()))
  summary[context]=out
 write(A/'all_group_moment_summary.json',dict(groups=900,contexts=summary,scope='All900 requested-versus-realized probability-moment diagnostics. This measures decoder fidelity, not correctness against hidden treatment profiles.',script_sha256=sha(__file__)));print('V90 all-group decoder summary',summary,flush=True)
