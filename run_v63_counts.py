"""V63 reuses the frozen paired count evaluator with one guarded primary."""
import argparse
from threadpoolctl import threadpool_limits
import v62_conserved_response as parent
import v63_identity_guard as guard
import run_v62_counts as runner
from v57_pie import read,write

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation','assess-development','assess-confirmation']);p.add_argument('--source');args=p.parse_args()
 runner.m=guard
 # The original evaluator calls m.data(), and reads its selected primary. Preserve that interface explicitly.
 guard.data=parent.data
 with threadpool_limits(limits=2):
  guard.check()
  if not (guard.A/'mean_results.json').exists():write(guard.A/'mean_results.json',dict(selected='guarded_assay',selection='Only fixed primary, no outcome-based selection',parent_mean_results_sha256=parent.sha(parent.A/'mean_results.json')))
  if args.stage.startswith('assess-'):runner.assess(args.stage.split('-',1)[1])
  else:runner.counts(args.stage,args.source)
