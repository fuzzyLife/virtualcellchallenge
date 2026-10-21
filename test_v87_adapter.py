"""Numerical fixtures for V87 conditional weight generation, not a biology benchmark."""
import copy
import numpy as np
import torch
from torch import nn
import v87_conditional_adapter as m
from v57_pie import write

def main():
 torch.manual_seed(87);layer=nn.Linear(5,7);adapter=m.AdapterLinear(layer,condition_dim=6,rank=3);x=torch.randn(11,5);c=torch.randn(11,6)
 torch.testing.assert_close(adapter(x,c),layer(x),rtol=0,atol=0)
 with torch.no_grad():adapter.b[-1].weight.normal_(std=.05);adapter.b[-1].bias.normal_(std=.05)
 a,b=adapter.factors(c);expected=torch.stack([torch.nn.functional.linear(x[i],layer.weight+.1*b[i]@a[i],layer.bias) for i in range(len(x))]);torch.testing.assert_close(adapter(x,c),expected,rtol=2e-6,atol=2e-6)
 loss=adapter(x,c).square().mean();loss.backward();assert adapter.base.weight.grad is None and adapter.base.bias.grad is None;assert adapter.a[-1].weight.grad.abs().sum()>0 and adapter.b[-1].weight.grad.abs().sum()>0
 base=m.Backbone();net=m.ConditionalOperator(base);z=torch.randn(8,m.RANK);cond=torch.randn(8,48);torch.testing.assert_close(net(z,cond),base(z),rtol=0,atol=0)
 # Finite-value ranks, including zeros, are invariant to positive monotone transforms.
 p=np.array([0.,1.,4.,4.,12.]);projection=np.random.default_rng(87).normal(size=(5,16));np.testing.assert_allclose(m.context_rank(p,projection),m.context_rank(p**3,projection),atol=0,rtol=0)
 record=dict(passed=True,zero_B_exact_base=True,explicit_per_example_collapsed_weights_match=True,base_gradients_absent=True,both_factor_gradients_nonzero=True,zero_adapter_exact_backbone=True,rank_monotonic_invariance=True,scope='Numerical implementation fixtures only, not performance evidence.')
 if (m.A/'plan.json').exists():m.check()
 write(m.A/'synthetic_tests.json',record);print('V87 adapter numerical fixtures passed')
if __name__=='__main__':main()
