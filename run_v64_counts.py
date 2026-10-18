"""Fixed V64 paired count tests, same audited evaluator as V62/V63."""
import argparse
from threadpoolctl import threadpool_limits
import v64_response_magnitude as m
import run_v62_counts as runner
runner.m=m
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','confirmation','assess-development','assess-confirmation']);p.add_argument('--source');a=p.parse_args()
 with threadpool_limits(limits=2):
  m.check();m.parent.check()
  if a.stage.startswith('assess-'):runner.assess(a.stage.split('-',1)[1])
  else:runner.counts(a.stage,a.source)
