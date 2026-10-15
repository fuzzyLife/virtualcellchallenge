"""Check deterministic scientific arithmetic and Atlas module IPC, not biology."""
import unittest
import numpy as np
from run_v27_with_processes import Processes
from evaluate_atlas_transfer import load_atlas_model
class IPC(unittest.TestCase):
 def test_frozen_atlas_serial_process_parity(self):
  atlas,_=load_atlas_model();rng=np.random.default_rng(27);ref=rng.poisson(6,(20,31));libs=ref.sum(1);p=(ref/libs[:,None]).mean(0);bulk=ref.sum(0)/libs.sum();q=rng.uniform(size=ref.shape);q/=q.sum(1,keepdims=True)
  def one(seed):return atlas.dual_moment_counts(q,p,bulk,depths=libs,seed=seed,iterations=200)
  serial=[one(s) for s in [27,28,29,30]]
  with Processes() as pool:parallel=list(pool.map(one,[27,28,29,30]))
  for a,c in zip(serial,parallel):np.testing.assert_array_equal(a,c)
if __name__=='__main__':unittest.main()
