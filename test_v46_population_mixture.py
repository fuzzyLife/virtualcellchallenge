"""Distribution mixtures preserve cells, integer counts and exact row libraries."""
import numpy as np
import v46_population_mixture as m
def main():
 a=np.array([[4,0,2],[0,4,2],[2,3,1],[1,1,4]],np.int32);b=np.array([[1,1,4],[2,1,3],[5,0,1],[0,3,3]],np.int32)
 sets=[]
 for f in [.25,.5,.75]:
  out,rows=m.mix(a,b,5,f);np.testing.assert_array_equal(out.sum(1),a.sum(1));assert len(rows)==round(f*len(a));assert out.dtype.kind in 'iu'
  for i in range(len(a)):np.testing.assert_array_equal(out[i],b[i] if i in rows else a[i])
  sets.append(set(rows))
 assert sets[0]<sets[1]<sets[2];out,_=m.mix(a,a,5,.5);np.testing.assert_array_equal(out,a)
 m.A.mkdir(parents=True,exist_ok=True);m.write(m.A/'contract_tests.json',{'passed':True,'checks':['exact integer row libraries','whole component cell vectors retained','nested balanced assignments','identical-component identity'],'script_sha256':m.b.sha256(m.ROOT/'test_v46_population_mixture.py')});print('V46 mixture contracts PASS',flush=True)
if __name__=='__main__':main()

