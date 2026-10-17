"""Test supported directional amplification and safe missingness."""
import numpy as np
import v47_kolf_consensus as m
def main():
 base=np.array([[.1,-.1,0.,.01,.2,.3]],np.float32);z=np.array([[.5,-.5,5.,2.,np.nan,100.]],np.float32);sat=np.array([[False,False,False,False,False,True]])
 out=m.combine(base,z,sat);assert out[0,0]>base[0,0] and out[0,1]<base[0,1];np.testing.assert_array_equal(out[0,2:],base[0,2:]);assert (abs(out)<=1.25*abs(base)+1e-7).all();np.testing.assert_array_equal(m.combine(base,-z,sat),base)
 m.A.mkdir(parents=True,exist_ok=True);m.write(m.A/'contract_tests.json',{'passed':True,'checks':['independent concordance only','zero/weak/missing/saturated exact fallback','positive bounded gains','opposed evidence abstains'],'script_sha256':m.b.sha256(m.ROOT/'test_v47_kolf_consensus.py')});print('V47 consensus contracts PASS',flush=True)
if __name__=='__main__':main()

