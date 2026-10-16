"""Check corrected transfer fallback with a low-expression synthetic query."""
import unittest,numpy as np
from unittest.mock import patch
import v32_methods as m
import v32_fallback_correction as c
class Fallback(unittest.TestCase):
 def test_unsupported_outputs_and_targets_are_exact_baseline(self):
  data={'genes':np.array(['T'+str(i) for i in range(300)]+['G'+str(i) for i in range(5)]),'targets':np.array(['T'+str(i) for i in range(300)])};base=np.full((300,305),.4,np.float32);base[np.arange(300),np.arange(300)]=0;query=np.ones(305);g={'axis':np.array([302,303]),'excluded':['H1']};mask=np.ones((300,2),bool);mask[0]=False
  with patch.object(m.old,'read_cache',return_value=data),patch.object(m,'estimate',return_value=(np.ones((300,2)),mask,[])):
   p,b,_=c.predict(g,query,.5,base);outside=np.ones(305,bool);outside[g['axis']]=False;np.testing.assert_array_equal(p[:,outside],base[:,outside]);np.testing.assert_array_equal(p[0],base[0]);self.assertGreater(np.max(abs(p[1]-base[1])),0)
if __name__=='__main__':unittest.main()