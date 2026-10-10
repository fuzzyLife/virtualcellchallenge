"""Meaningful leakage and masked-assay checks for the new response program model."""
import unittest,json
import numpy as np
from threadpoolctl import threadpool_limits
from fit_v7_response import read_cache
from v13_methods import A,ROOT
from v13_robust_methods import banks,NAMES
import build_submission as build

class ResponsePrograms(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  original=read_cache();targets=original['targets'][:80];genes=original['genes']
  cols=np.flatnonzero(np.isin(genes,targets)|np.isin(np.arange(len(genes)),np.arange(30)))
  rows=np.flatnonzero(original['row_targets']<80);data={}
  for key in ['truth','errors','q','valid']:data[key]=original[key][rows][:,cols].copy()
  for key in ['row_targets','row_families','row_contexts','knockdown']:data[key]=original[key][rows].copy()
  data.update(genes=genes[cols],targets=targets)
  cls.data=data
 def test_held_family_cannot_change_predictions(self):
  first,record=banks(self.data,excluded='HCT116',reserved=[1,2,3])
  changed={k:v.copy() for k,v in self.data.items()};held=changed['row_families']=='HCT116'
  changed['truth'][held]=2.9;changed['errors'][held]=.001;changed['valid'][held]=~changed['valid'][held];changed['knockdown'][held]=.99
  second,other=banks(changed,excluded='HCT116',reserved=[1,2,3])
  self.assertNotIn('HCT116',record['training_families'])
  for name in NAMES:np.testing.assert_array_equal(first[name]['consensus20'],second[name]['consensus20'])
  for name,bank in first.items():
   value=bank['consensus20'];self.assertTrue(np.isfinite(value).all());self.assertLessEqual(abs(value).max(),3.)
   own=np.array([np.flatnonzero(self.data['genes']==t)[0] for t in self.data['targets']]);np.testing.assert_array_equal(value[np.arange(80),own],0)
  build.write_json(A/'robust_source_checks.json',{'passed':True,'held_family_mutation_exact_invariance':True,
   'finite_bounded_own_masked_effects':True,'scope':'Actual source-cache subset: 80 targets and their output coordinates plus 30 other genes; HCT116 responses/errors/masks/KD mutated together.',
   'script_sha256':build.sha256(ROOT/'test_v13_robust.py'),'methods_sha256':build.sha256(ROOT/'v13_robust_methods.py')})
if __name__=='__main__':
 with threadpool_limits(limits=2):unittest.main()
