"""Execute a frozen V57 stage with immutable projection reads cached.

Only read-only mmap handles and immutable JSON are reused. Numerical functions,
training samples, RNGs and frozen code are unchanged.
"""
from functools import lru_cache
from pathlib import Path
import runpy,sys
import numpy as np
import build_v57_pie_evidence as evidence
from v58_pie import ROOT,A,CACHE,PARENT_CACHE,read,write,sha

class ReadOnlyNumpy:
 def __getattr__(self,name):return getattr(np,name)
 @lru_cache(maxsize=64)
 def cached(self,path):return np.load(path,mmap_mode='r')
 def load(self,path,*args,**kwargs):
  p=Path(path)
  if not args and kwargs=={'mmap_mode':'r'} and p.parent==PARENT_CACHE/'evidence_projection':return self.cached(str(p))
  return np.load(path,*args,**kwargs)

def main():
 script=sys.argv[1];args=sys.argv[2:];assert script in ['fit_v58_pie_decoder.py','evaluate_v58_pie_counts.py','summarize_v58_panels.py','audit_v58_pie.py','document_v58_pie.py']
 meta=read(PARENT_CACHE/'evidence_projection/projection.json');p=read(A/'plan.json');axis=meta['genes'][::997]
 targets=[p['seen_targets'][0],p['unseen_targets'][0]];before={t:evidence.blocks(t,'CD4_Rest_D2',axis)[0] for t in targets}
 evidence.np=ReadOnlyNumpy();evidence.read=lru_cache(maxsize=None)(evidence.read)
 for t in targets:
  after=evidence.blocks(t,'CD4_Rest_D2',axis)[0]
  for k in after:np.testing.assert_array_equal(before[t][k],after[k])
 print('Read-only projection cache parity passed:',script,flush=True)
 write(A/'projection_cache_execution.json',{'passed':True,'wrapper_sha256':sha(Path(__file__)),'frozen_methods_unchanged':True,'parity_targets':targets,'parity_genes':len(axis),'change':'Reuse read-only mmap handles and immutable JSON; no array values or calculations changed'})
 sys.argv=[str(ROOT/script),*args];runpy.run_path(str(ROOT/script),run_name='__main__')

if __name__=='__main__':main()
