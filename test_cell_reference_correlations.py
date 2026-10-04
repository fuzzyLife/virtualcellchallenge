"""Independent numerical checks for the correlation analysis."""
import unittest
import numpy as np
from scipy.stats import spearmanr,rankdata
from analyze_cell_reference_correlations import pearson_rows,normalized_log
from prepare_reference_controls import select_stratified
import pandas as pd
class CorrelationTests(unittest.TestCase):
 def test_numpy_and_scipy_oracles(self):
  cells=np.array([[0,3,0,10,2],[4,0,0,1,2],[1,1,1,1,1]],dtype=float)
  refs=np.array([[3,2,0,6,1],[0,2,7,0,2]],dtype=float)
  a=normalized_log(cells);b=normalized_log(refs);actual=pearson_rows(a,b)
  for i in range(2):
   for j in range(2):
    self.assertAlmostEqual(actual[i,j],np.corrcoef(a[i],b[j])[0,1],places=12)
    ranked=pearson_rows(rankdata(a[i])[None,:],rankdata(b[j])[None,:])[0,0]
    self.assertAlmostEqual(ranked,spearmanr(a[i],b[j]).statistic,places=12)
  self.assertTrue(np.isnan(actual[2]).all())
  np.testing.assert_allclose(normalized_log(cells*7),a,rtol=1e-14)
 def test_stratification_and_invalid_input(self):
  frame=pd.DataFrame({'batch':['a','a','b','b'],'target':['ctrl','kd','ctrl','ctrl']})
  chosen=select_stratified(frame,frame.target.eq('ctrl').to_numpy(),['batch'],1)
  self.assertEqual(chosen.source_row.tolist(),[0,2])
  with self.assertRaises(ValueError):normalized_log(np.zeros((1,3)))
  with self.assertRaises(ValueError):normalized_log(np.array([[1,-1,3]]))
class RangeReaderTests(unittest.TestCase):
 def test_range_cache_integrity_and_budget(self):
  from unittest.mock import patch
  from tempfile import TemporaryDirectory
  from pathlib import Path
  from prepare_reference_controls import RangeReader
  payload=bytes(range(256))*3000
  class Response:
   def __init__(self,lo,hi):
    self.status_code=206;self.content=payload[lo:hi+1]
    self.headers={'Content-Range':f'bytes {lo}-{hi}/{len(payload)}','ETag':'"fixed"'}
   def __enter__(self):return self
   def __exit__(self,*args):return False
  class Session:
   def mount(self,*args):pass
   def close(self):pass
   def get(self,url,headers,**kwargs):
    lo,hi=map(int,headers['Range'].split('=')[1].split('-'))
    return Response(lo,min(hi,len(payload)-1))
  with TemporaryDirectory() as directory,patch('prepare_reference_controls.requests.Session',Session):
   reader=RangeReader('https://test.invalid/counts',directory,budget=2*262144)
   reader.seek(262140);self.assertEqual(reader.read(16),payload[262140:262156]);reader.close()
   block=Path(directory)/'000000000.bin';block.write_bytes(b'corrupt')
   reader=RangeReader('https://test.invalid/counts',directory,budget=262144)
   with self.assertRaisesRegex(RuntimeError,'checksum'):reader.read(1)
   reader.close()
  with TemporaryDirectory() as directory,patch('prepare_reference_controls.requests.Session',Session):
   reader=RangeReader('https://test.invalid/counts',directory,budget=262144)
   with self.assertRaisesRegex(RuntimeError,'budget'):reader.read(262145)
   reader.close()
if __name__=='__main__':unittest.main()
