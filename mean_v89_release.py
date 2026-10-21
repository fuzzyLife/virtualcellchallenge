"""Recompute V89 public means with the immutable V86 evaluator in a separate namespace."""
import importlib.util
from threadpoolctl import threadpool_limits
import v89_projection_release as m
from v57_pie import read,write,sha
if __name__=='__main__':
 with threadpool_limits(limits=2):
  m.check();spec=importlib.util.spec_from_file_location('v89_public_evaluation',m.ROOT/'run_v86_experiment.py');runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner);runner.m=m;runner.means();write(m.A/'mean_replay_receipt.json',dict(original_evaluator_sha256=sha(m.ROOT/'run_v86_experiment.py'),unchanged_evaluation_code=True,source_predictions_exact_V86=True,script_sha256=sha(__file__)))
