"""Run V33 with audited fallback and fitted feature-scaling adapters."""
from functools import lru_cache
from threadpoolctl import threadpool_limits
import json
import v32_methods as m
import v32_fallback_correction as correction
import run_v33_fixed as run
if __name__=='__main__':
 dest=run.A/'fallback_correction_amendment.json'
 if dest.exists() and not (run.A/'fallback_correction_initial_amendment.json').exists():m.write(run.A/'fallback_correction_initial_amendment.json',json.loads(dest.read_text()))
 m.write(dest,{'parent_plan_sha256':m.sha(run.A/'plan.json'),'correction_sha256':m.sha(m.ROOT/'v32_fallback_correction.py'),'wrapper_sha256':m.sha(m.ROOT/'run_v33_corrected.py'),'adapter_sha256':m.sha(m.ROOT/'run_v33_fixed.py'),'reason':'Before V33 model scoring, retain exact V24 outside supported predictions and attach the paired table fitting feature_mean/scale to the query adapter. Initial execution stopped at missing adapter key before any model scores, fits or hyperparameter outcomes. Original frozen code retained. Same models/grid/splits/gates.'});correction.install();m.old.embeddings=lru_cache(maxsize=1)(m.old.embeddings)
 with threadpool_limits(limits=2):run.main()