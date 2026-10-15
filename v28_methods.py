"""V28 global-only mean repair; prospective follow-up to spent null diagnostic."""
from pathlib import Path
import json
import numpy as np
from scipy import sparse
from v27_methods import ROOT,anchored_initial,control_calibration,relative_delta,write,digest,logcpm
from v26_methods import fit_probabilities,integerize,repair_cpm,load_decoder as old_decoder
import build_submission as b
OUT=ROOT/'experiments/v28_global_transport';A=OUT/'analysis'
ARMS=['v24','anchor_only','calibrated_anchor','relative','calibrated_relative','panel_relative','target_shuffled']
FINALIST='anchor_only'
def arms_for(key):return ARMS if key=='H1' else ['v24',FINALIST]
def load_plan():
 p=json.loads((A/'plan.json').read_text())
 for n,h in p['dependencies'].items():assert b.sha256(ROOT/n)==h,'Frozen V28 dependency changed: '+n
 for n,h in p['inputs'].items():assert b.sha256(ROOT/n)==h,'Frozen V28 input changed: '+n
 return p
def load_model(key,panel=False):
 if key!='H1':return {'record':{'sha256':None,'scope':'No neural model used for global anchor'}}
 if not panel:return old_decoder('H1','conditional')
 from v27_methods import load_model as previous
 return previous('H1',True)
def project(initial,reference,seed):
 ref=reference.toarray().astype(np.int32) if sparse.issparse(reference) else np.asarray(reference,np.int32);libs=ref.sum(1,dtype=np.int64);p=(ref/libs[:,None]).mean(0);bulk=ref.sum(0)/libs.sum();q=fit_probabilities(initial,libs,p,bulk);x,n=integerize(q,libs,ref.sum(0,dtype=np.int64),seed);x,k=repair_cpm(x,ref);drift=float(abs((x/libs[:,None]).mean(0)-p).sum());return sparse.csr_matrix(x),{'mean_reference_L1':drift,'bulk_exact':True,'global_repair_molecules':n,'CPM_repair_molecules':k,'depth_stratum_repair':False}
def emit(template,effect,own,p,seed,reference,arm,model,state,fit=None,override=None):
 calibrated=control_calibration(template,fit,p,seed) if arm in ['calibrated_anchor','calibrated_relative','panel_relative','target_shuffled'] else None
 delta=None if arm in ['anchor_only','calibrated_anchor'] else relative_delta(state,effect,p,model,seed,override)
 if own is not None and delta is not None:delta[:,own]=0
 initial=anchored_initial(template,effect,own,p,delta,calibrated)
 return project(initial,reference,seed+270270)
