"""V34 pre-outcome integration fix, preserving frozen initial source and plan."""
import argparse
from functools import lru_cache
from threadpoolctl import threadpool_limits
import v34_corrected_methods as fixed
import run_v34_experiment as run
import v32_methods as m
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['train','mean','null','treatment']);p.add_argument('--outer',choices=m.OUTERS);args=p.parse_args();m.old.embeddings=lru_cache(maxsize=1)(m.old.embeddings);run.f=fixed
 dest=fixed.A/'integration_correction_amendment.json'
 if not dest.exists():m.write(dest,{'parent_plan_sha256':m.sha(fixed.A/'plan.json'),'original_source_sha256':m.sha(m.ROOT/'v34_methods.py'),'corrected_source_sha256':m.sha(m.ROOT/'v34_corrected_methods.py'),'wrapper_sha256':m.sha(m.ROOT/'run_v34_corrected.py'),'reason':'Pre-outcome code review: an own-gene label is absent during fitting, so its unconstrained flow head must not alter the autoregressive input during generative integration. Hold its state at the control value at every step in source validation and query inference. No third-model fits/outcomes existed at this amendment. Preserves the initial source and plan.'})
 else:
  import json
  rec=json.loads(dest.read_text());assert rec['corrected_source_sha256']==m.sha(m.ROOT/'v34_corrected_methods.py') and rec['wrapper_sha256']==m.sha(m.ROOT/'run_v34_corrected.py')
 with threadpool_limits(limits=2):run.main(args.stage,args.outer)