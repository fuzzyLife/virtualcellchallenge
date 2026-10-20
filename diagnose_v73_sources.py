"""Compare weak-evidence source estimators without using recipient outcomes."""
import numpy as np
from pathlib import Path
from v57_pie import read,write,sha
import v73_count_assay as m
rows=[]
for name in m.SOURCES:
 p=m.A/'sources'/(name+'.npz');old=m.ROOT/'experiments/v72_batch_offset/analysis/sources'/(name+'.npz')
 with np.load(p,allow_pickle=False) as z:y=z['effect'].copy();v=z['variance'].copy();ok=z['valid'].copy();b=z['baseline_bulk'].copy()
 with np.load(old,allow_pickle=False) as z:previous=z['effect'].copy();pv=z['variance'].copy();present=z['valid'].copy()
 np.testing.assert_array_equal(ok,present)
 used=ok.any(1);np.testing.assert_allclose(b[used].sum(1),1,rtol=1e-6)
 assert np.isfinite(y).all() and np.isfinite(v).all() and (v>0).all() and np.max(abs(y))<=3
 weak=ok&(pv>1);strong=ok&(pv<.01)
 stats={}
 for label,mask in [('weak',weak),('strong',strong)]:
  stats[label]=dict(coordinates=int(mask.sum()),previous_mean=float(previous[mask].mean()) if mask.any() else None,new_mean=float(y[mask].mean()) if mask.any() else None,previous_mean_absolute=float(abs(previous[mask]).mean()) if mask.any() else None,new_mean_absolute=float(abs(y[mask]).mean()) if mask.any() else None)
 rows.append(dict(source=name,statistics=stats,previous_sha256=sha(old),new_sha256=sha(p)))
write(m.A/'source_estimator_diagnostic.json',dict(records=rows,scope='Source estimator outputs only; no independent truth, recipient outcome or accuracy evidence.',script_sha256=sha(__file__)))
print('Source-estimator diagnostic saved; all source masks/axes and matched bulk closures agree.')
