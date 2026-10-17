"""Nondegenerate shape-preservation test for V48 integer-margin projection."""
import numpy as np
import v48_margin_preserving as m
def main():
 rng=np.random.default_rng(48048);n=64;g=12;libs=rng.integers(800,1200,n);p=np.ones(g)/g
 ref=np.array([rng.multinomial(int(l),p) for l in libs],np.int32)
 initial=np.array([rng.multinomial(int(l),rng.dirichlet(p*20)) for l in libs],np.int32)
 out,diag=m.project(initial,ref,48048);x=out.toarray()
 np.testing.assert_array_equal(x.sum(1),libs);np.testing.assert_array_equal(x.sum(0),ref.sum(0))
 assert np.any(x!=ref),'Projection collapsed all cells to the reference'
 assert np.mean(np.var(x/libs[:,None],axis=0))>np.mean(np.var(ref/libs[:,None],axis=0))
 m.write(m.A/'nondegenerate_variation_test.json',{'passed':True,'cells':n,'genes':g,'nonreference_population_shape_retained':True,'exact_libraries_and_totals':True,'mean_reference_L1':diag['mean_reference_L1'],'script_sha256':m.b.sha256(m.ROOT/'test_v48_variation.py')})
 print('V48 nondegenerate variation test PASS',flush=True)
if __name__=='__main__':main()

