"""Synthetic identification and abstention checks for the new guide-contrast mechanism."""
from unittest.mock import patch
import numpy as np
import pandas as pd
from scipy import sparse
import v85_guide_contrast as m
from v57_pie import write,sha

def main():
    beta=np.array([-.2,.3,0.]);dose=np.array([-.8,-.8]);response=dose[:,None]*beta
    effect,conf,diag=m.estimate(dose,np.array([.0001,.0001]),response,np.full((2,3),.0001))
    np.testing.assert_allclose(effect,np.log2(1-.9*beta),atol=1e-12);assert diag['strong'] and conf[0]>0 and conf[1]>0
    swapped=m.estimate(-dose,np.array([.0001,.0001]),-response,np.full((2,3),.0001))[0];np.testing.assert_array_equal(effect,swapped)
    weak=m.estimate(np.array([-.01,-.01]),np.ones(2),response,np.ones((2,3)));assert not weak[2]['strong'] and not weak[1].any()
    opposing=m.estimate(np.array([-.8,.8]),np.full(2,.001),response,np.ones((2,3)));assert not opposing[2]['strong']
    zero=m.estimate(dose,np.full(2,.001),np.zeros((2,3)),np.ones((2,3)));np.testing.assert_array_equal(zero[0],np.zeros(3))
    # Known guide contrast plus a common trans response: subtraction removes it.
    genes=np.array(['T','G1','G2']);targets=np.array(['T']);true_beta=np.array([.2,-.3]);rows=[];meta=[]
    for lane in ['L0','L1']:
        for guide,ratio in [('NTC',1.),('g1',.1),('g2',.5)]:
            mean=np.r_[100*ratio,1000*(1+true_beta*(ratio-1)+(0 if guide=='NTC' else .4))]
            for i in range(20):rows.append(mean);meta.append(dict(gene_target='Non-Targeting' if guide=='NTC' else 'T',guide_id=guide,sample=lane))
    with patch.object(m,'data',return_value=dict(genes=genes,targets=targets)):
        fit,records,stats=m.fit_groups(sparse.csr_matrix(np.array(rows)),pd.DataFrame(meta),np.ones(3,bool))
    np.testing.assert_allclose(fit['effect'][0,1:],np.log2(1-.9*true_beta),atol=.003);assert records[0]['active_outputs']==2
    used=[i for lane in records[0]['rows'] for guide in lane['halves'] for half in guide for i in half];assert len(used)==len(set(used))==80
    assert abs(fit['pooled'][0,1]-fit['effect'][0,1])>.1
    write(m.A/'synthetic_tests.json',dict(passed=True,known_dose_response_recovered=True,shared_response_cancels=True,weak_or_opposing_dose_abstains=True,guide_label_swap_invariant=True,zero_response_exact=True,cell_halves_disjoint=True,code_sha256=sha(__file__),estimator_sha256=sha(m.ROOT/'v85_guide_contrast.py')))
    print('V85 synthetic mechanism tests passed',flush=True)
if __name__=='__main__':main()
