"""Independent V87 study exclusion, basis, checkpoint selection and matrix audits."""
import argparse,copy
import numpy as np
import torch
from threadpoolctl import threadpool_limits
import v87_conditional_adapter as m
from v57_pie import read,write,sha

def models():
 m.check();d=m.data();records=[]
 for outer in ['H1','Orion','CD4']:
  b=m.make_bank(d,outer);poison=dict(d);rows=np.array([m.GROUP[str(f)]==outer for f in d['row_families']]);contexts=np.unique(d['row_contexts'][rows])
  for name in ['truth','errors','knockdown','valid']:
   poison[name]=d[name].copy();poison[name][rows]=False if name=='valid' else 98765
  poison['profiles']=d['profiles'].copy();poison['profiles'][contexts]=123456;b2=m.make_bank(poison,outer)
  for name in ['y','present','axis','profiles','source_indices','source_names']:np.testing.assert_array_equal(b[name],b2[name])
  assert not np.isin(d['genes'][b['axis']],d['targets']).any()
  for fake in [False,True]:
   arr,nets=m.load(outer,fake);meta=read(m.A/'models'/(outer+('_shuffled' if fake else '')+'.json'));np.testing.assert_allclose(arr['basis']@arr['basis'].T,np.eye(m.RANK),rtol=0,atol=5e-5);assert not np.intersect1d(arr['training_targets'],arr['validation_targets']).size;assert len(arr['training_targets'])==240 and len(arr['validation_targets'])==60;assert all(m.GROUP[n]!=outer for n in arr['source_names'])
   for trace in [meta['base'],*meta['adapters'].values()]:
    best=min(trace['history'],key=lambda z:z['validation_loss']);assert trace['selected_epoch']==best['epoch'];assert trace['selected_validation_loss']==best['validation_loss'];assert len(trace['history'])==101
   for name,net in nets.items():
    for layer in [net.first,net.last]:
     assert all(not p.requires_grad for p in layer.base.parameters());x=torch.randn(3,layer.base.in_features);c=torch.randn(3,48);a,bb=layer.factors(c);expected=torch.stack([torch.nn.functional.linear(x[i],layer.base.weight+.1*bb[i]@a[i],layer.base.bias) for i in range(3)]);torch.testing.assert_close(layer(x,c),expected,rtol=2e-5,atol=2e-5)
   records.append(dict(outer=outer,shuffled=fake,complete_output_genes=len(arr['axis']),training_targets=240,held_targets=60,all_selection_histories_replayed=True,conditional_matrix_equations_replayed=True))
  print('V87 independent model audit',outer,flush=True)
 write(m.A/'model_audit.json',dict(passed=True,whole_recipient_study_poisoning_exact=True,records=records,code_sha256=sha(__file__)))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['models']);a=p.parse_args();torch.set_num_threads(2)
 with threadpool_limits(limits=2):models()
