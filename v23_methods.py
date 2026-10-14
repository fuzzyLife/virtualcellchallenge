"""V23 uses frozen V22 factors; adds V16-reference integer moment repair only."""
from pathlib import Path
import json,threading
import numpy as np
from scipy import sparse
import build_submission as b
import v22_methods as parent
from v23_mean_repair import preserve_reference_means
ROOT=b.ROOT;OUT=ROOT/'experiments/v23_mean_preserving';A=OUT/'analysis';FAMILIES=parent.FAMILIES
NAMES=['v16','conditional','mean_preserved','scalar_preserved','independent_preserved','shuffled_preserved']
write=parent.write;digest=parent.digest;logcpm=parent.logcpm;raw=parent.raw;query_states=parent.query_states
LOCAL=threading.local()
def load_plan():
 plan=json.loads((A/'plan.json').read_text());parent.load_plan()
 for name,h in plan['dependencies'].items():assert b.sha256(ROOT/name)==h,'Frozen V23 dependency changed: '+name
 for name,h in plan['parent_models'].items():assert b.sha256(parent.A/'models'/name)==h
 return plan
def load_model(family):
 m=parent.load_model(family);assert b.sha256(A/'models'/(family+'.npz'))==b.sha256(parent.A/'models'/(family+'.npz'));return m
def fit_model(family):return load_model(family)
def emit(template,effect,own,seed,kind,prepared,atlas,model=None,ti=None,z=None):
 key=(seed,own,ti,parent.digest(template),np.asarray(effect).tobytes(),np.asarray(prepared[1]).tobytes(),np.asarray(prepared[2]).tobytes())
 if kind=='v16':
  count,diag=parent.emit(template,effect,own,seed,'v16',prepared,atlas,model,ti,z);LOCAL.key=key;LOCAL.reference=count;return count,diag
 if kind=='conditional':return parent.emit(template,effect,own,seed,'conditional',prepared,atlas,model,ti,z)
 mapping={'mean_preserved':'conditional','scalar_preserved':'scalar','independent_preserved':'independent','shuffled_preserved':'shuffled_state'}
 assert kind in mapping
 if getattr(LOCAL,'key',None)==key:reference=LOCAL.reference
 else:
  reference,_=parent.emit(template,effect,own,seed,'v16',prepared,atlas,model,ti,z);LOCAL.key=key;LOCAL.reference=reference
 count,diag=parent.emit(template,effect,own,seed,mapping[kind],prepared,atlas,model,ti,z)
 if (count!=reference).nnz==0:return count,dict(diag,exact_V16_gene_totals=True,mean_CPM_probability_L1_after=0.)
 fixed,repair=preserve_reference_means(count.toarray(),reference.toarray(),seed=seed)
 # Repair never uses accessible or hidden true treatments; reference is generated V16.
 return sparse.csr_matrix(fixed),dict(diag,**repair)
