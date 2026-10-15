"""V29 pair-aware global transport; no change to source mean prediction."""
from pathlib import Path
import json,importlib.metadata
import numpy as np
from scipy import sparse
from v28_methods import ROOT,ARMS,FINALIST,arms_for,load_model
from v27_methods import anchored_initial,control_calibration,relative_delta,write,digest,logcpm
from v26_methods import fit_probabilities,integerize,repair_cpm
from v29_pairwise_repair import repair_pairwise
import build_submission as b
OUT=ROOT/'experiments/v29_pairwise_transport';A=OUT/'analysis'
def load_plan():
 p=json.loads((A/'plan.json').read_text())
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,'Frozen V29 dependency changed: '+n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,'Frozen V29 input changed: '+n
 for n,v in p['execution_versions'].items():assert importlib.metadata.version(n)==v,'V29 backend changed: '+n
 return p

def project(initial,reference,seed):
 ref=reference.toarray().astype(np.int32) if sparse.issparse(reference) else np.asarray(reference,np.int32);libs=ref.sum(1,dtype=np.int64);p=(ref/libs[:,None]).mean(0);bulk=ref.sum(0)/libs.sum();q=fit_probabilities(initial,libs,p,bulk);x,n=integerize(q,libs,ref.sum(0,dtype=np.int64),seed);fixed,diag=repair_pairwise(x,ref);diag['repair_path']='joint';total=int(libs.sum());budget=int(total*.05)
 # A deterministic numerical bridge is permitted only when the hard mean
 # contract is unmet. It reads predicted reference margins, never true outcomes.
 if diag['mean_reference_L1']>.0005 and budget-diag['CPM_repair_molecules']>=4:
  used=diag['CPM_repair_molecules'];left=budget-used;cap=min(int(total*.01),left//2)
  carrier,cost=repair_cpm(fixed.copy(),ref,max_fraction=cap/total)
  trial,d2=repair_pairwise(carrier,ref,max_fraction=max(0,left-cost)/total)
  if d2['mean_reference_L1']<diag['mean_reference_L1']:
   fixed=trial;diag={**d2,'repair_path':'joint_carrier_joint','CPM_repair_molecules':used+cost+d2['CPM_repair_molecules'],'CPM_repair_budget':budget,'carrier_molecules':cost}
 assert diag['CPM_repair_molecules']<=budget
 rr=ref.T@(1/libs)-fixed.T@(1/libs);diag['maximum_gene_CPM_probability_error']=float(abs(rr).max()/len(ref));diag.update({'bulk_exact':True,'global_repair_molecules':n,'depth_stratum_repair':False})
 return sparse.csr_matrix(fixed),diag

def emit(template,effect,own,p,seed,reference,arm,model,state,fit=None,override=None):
 calibrated=control_calibration(template,fit,p,seed) if arm in ['calibrated_anchor','calibrated_relative','panel_relative','target_shuffled'] else None
 delta=None if arm in ['anchor_only','calibrated_anchor'] else relative_delta(state,effect,p,model,seed,override)
 if own is not None and delta is not None:delta[:,own]=0
 initial=anchored_initial(template,effect,own,p,delta,calibrated)
 return project(initial,reference,seed+270270)
