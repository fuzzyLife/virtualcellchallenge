"""Nondegenerate bounded projection preserves marginals without collapsing shape."""
import numpy as np
import v49_balanced_population as m
def main():
 rng=np.random.default_rng(49);n=64;g=30;libs=rng.integers(800,1200,n);p=np.ones(g)/g
 ref=np.array([rng.multinomial(int(l),p) for l in libs],np.int32);x=np.array([rng.multinomial(int(l),rng.dirichlet(p*100)) for l in libs],np.int32)
 out,diag=m.project(x,ref,49);xx=out.toarray();np.testing.assert_array_equal(xx.sum(1),libs);np.testing.assert_array_equal(xx.sum(0),ref.sum(0));assert np.any(xx!=ref);assert diag['maximum_entry_rounding_error']<=1 and diag['mean_reference_L1']<=.0005
 m.write(m.A/'population_contract_tests.json',{'passed':True,'nonreference_shape_retained':True,'exact_gene_and_cell_totals':True,'diagnostics':diag,'script_sha256':m.b.sha256(m.ROOT/'test_v49_population.py')});print('V49 population contracts PASS',diag['depth_bins'],flush=True)
if __name__=='__main__':main()

