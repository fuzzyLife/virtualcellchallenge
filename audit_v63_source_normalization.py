"""Report measured source-axis coverage; do not mutate source profiles or predictions."""
from pathlib import Path
import numpy as np
from v57_pie import sha,write
ROOT=Path(__file__).resolve().parent;source=ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'
with np.load(source,allow_pickle=False) as z:profiles=z['profiles'].copy()
names=['K562','CD4_Rest','CD4_Stim8hr','CD4_Stim48hr','H1','HCT116','HEK293T']
rows=[dict(source=n,finite_mapped_genes=int(np.isfinite(p).sum()),sum_available_mean_CPM=float(np.nansum(p)),available_fraction_of_one_million=float(np.nansum(p)/1e6)) for n,p in zip(names,profiles)]
write(ROOT/'experiments/v63_identity_guard/analysis/source_normalization_audit.json',dict(source_cache_sha256=sha(source),script_sha256=sha(__file__),profiles=rows,interpretation='Mapped source profiles contain roughly68-100% of a complete per-cell CPM library. The frozen V62/V63 code uses these available CPM values directly in its log1p response transform; it does not reclose each source profile. The initial plan word closed is inaccurate for source profiles. Query controls are aligned-axis normalized; source profiles are incomplete-axis per-cell proxies, not exact pooled profiles on the private scoring gene universe. No predictor, score or source data was changed after this audit. Reclosing an incompletely measured source axis would be a separate estimator assumption, not an established correction to missing biology.'))
print('Source normalization scope recorded; no numerical data changed.')
