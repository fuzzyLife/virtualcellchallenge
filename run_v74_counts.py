"""Thin V74 binding of the frozen V70 paired-count evaluator."""
import argparse
from threadpoolctl import threadpool_limits
import v74_assay_gain as m
import run_v70_counts as r
r.m=m
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['counts','assess']);p.add_argument('--source');a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='counts':r.counts(a.source)
  else:r.assess()
