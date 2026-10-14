"""Scientific invariants for compatible signed mean corrections."""
import unittest,numpy as np
from unittest.mock import patch
from v24_methods import predict,training_data
class Means(unittest.TestCase):
 def test_withheld_rows(self):
  d={'row_families':np.array(['CD4','CD4','H1']),'row_contexts':np.array([1,3,4])}
  for k in ['truth','errors','valid','q','knockdown','row_targets']:d[k]=np.arange(3)
  a=training_data(d);np.testing.assert_array_equal(a['row_contexts'],[1,4]);self.assertIs(training_data(d,True),d)
 def test_zero_gate_and_bounded_correction(self):
  g=200;genes=np.array(['T']+[f'g{i}' for i in range(g-1)]);data={'genes':genes,'targets':np.array(['T']),'profiles':np.ones((7,g))};base=np.full((1,g),.2,np.float32);base[0,0]=0;teacher={'profiles':np.ones((2,g)),'state_axis':np.arange(1,g),'value':np.ones((2,1,g)),'confidence':np.ones((2,1,g))}
  # Production shape uses300targets; this fixture tests the full source contract.
  data['targets']=np.array(['T']*300);base=np.repeat(base,300,0);teacher['value']=np.repeat(teacher['value'],300,1);teacher['confidence']=np.repeat(teacher['confidence'],300,1)
  teacher['confidence_de']=teacher['confidence'].copy();bank={'base':base,'own':np.zeros(300,int)}
  data.update(row_families=np.array(['K562']),row_contexts=np.array([0]),row_targets=np.array([0]),truth=np.zeros((1,g)),errors=np.ones((1,g)),valid=np.ones((1,g),bool),q=np.ones((1,g)),knockdown=np.array([.9]))
  with patch('v24_methods.source_arrays',return_value=bank),patch('v24_methods.pairwise_spearman',return_value=(.8,g,'ok')):
   result,d=predict(data,teacher,np.ones(g)*20,baseline=base);np.testing.assert_array_equal(result['cdf50'],base)
  calls=[(.5,g,'ok')]*4+[(.9,g,'ok')]*2+[(.9,g,'ok')]*2
  with patch('v24_methods.source_arrays',return_value=bank),patch('v24_methods.pairwise_spearman',side_effect=calls):
   result,d=predict(data,teacher,np.ones(g)*20,baseline=base);self.assertLessEqual(abs(result['cdf50']-base).max(),.250001);self.assertTrue((result['cdf50'][:,0]==0).all())
  teacher['value'][1]=-1;teacher['confidence'][1]=0;teacher['confidence_de'][1]=0
  calls=[(.5,g,'ok')]*4+[(.9,g,'ok')]*4
  with patch('v24_methods.source_arrays',return_value=bank),patch('v24_methods.pairwise_spearman',side_effect=calls):
   result,d=predict(data,teacher,np.ones(g)*20,baseline=base);np.testing.assert_allclose(result['cdf50'][:,1:],.325,atol=1e-6)
if __name__=='__main__':unittest.main()
