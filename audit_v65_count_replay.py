"""Freshly regenerate ten finalist count populations and verify reach tie ordering."""
import argparse,importlib,gc
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
from v57_pie import read,write,sha
from run_v35_v36 import contexts
from run_v37_ensemble import checked_folder
from v16_methods import emit
from v12_cpu_metrics import de
from audit_transfer_validation import profile_jackknife_chunked

def main(version):
 m=importlib.import_module('v65_closed_response' if version==65 else 'v66_assay_strength');m.check();selected=read(m.A/'mean_results.json')['selected'];d=m.old.data();genes=d['genes'];targets=d['targets'];atlas,_=__import__('evaluate_atlas_transfer').load_atlas_model();checks=[]
 for ctx in contexts():
  name=ctx['family'];folder,_=checked_folder(name);p=m.A/'predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
  with np.load(p,allow_pickle=False) as z:eff=z[selected].copy()
  with np.load(folder/'effects.npz',allow_pickle=False) as z:cp=z['cp'].copy();bulk=z['bulk'].copy()
  reference=sparse.load_npz(folder/'reference.npz');i=0;cell=ctx['cells'][i];ti=list(targets).index(cell['target']);own=list(genes).index(cell['target']);template=sparse.load_npz(folder/f'{i}_template.npz');libs=template.sum(1).A1.astype(np.int64);actual=sparse.load_npz(folder/f'{i}_actual.npz');truth=de(actual,reference,genes)
  for offset in (0,620000):
   path=m.A/'counts'/f'{name}__{i}__{offset}__{selected}.json';r=read(path);count,_=emit(template,np.arange(len(libs)),eff[ti],own,cell['seed']+ti+offset,'residual10',(libs,cp,bulk),atlas);np.testing.assert_array_equal(count.sum(1).A1,libs);pred=de(count,reference,genes);pp,pj=profile_jackknife_chunked(count)
   with np.load(path.with_suffix('.npz'),allow_pickle=False) as z:
    np.testing.assert_array_equal(pp,z['pred']);np.testing.assert_array_equal(pj,z['jack'])
    for k in ('q','lfc','eligible'):np.testing.assert_array_equal(pred[k],z['pred_'+k])
   valid=truth['eligible']&(np.arange(len(genes))!=own);real=valid&(truth['q']<.05);pool=[j for j in range(len(genes)) if real[j] and np.isfinite(truth['lfc'][j]) and truth['lfc'][j]!=0]
   pool.sort(key=lambda j:(pred['q'][j]>=.05,pred['q'][j],pred['p'][j],-abs(pred['lfc'][j]),str(genes[j])))
   best=0;correct=0
   for rank,j in enumerate(pool,1):
    correct+=int(np.sign(pred['lfc'][j])==np.sign(truth['lfc'][j]))
    if correct/rank>=.9:best=rank
   expected=best/int(real.sum()) if real.any() else None
   assert (r['reach'] is None and expected is None) or (r['reach'] is not None and expected is not None and abs(r['reach']-expected)<1e-12)
   checks.append(dict(source=name,seed_offset=offset,target=cell['target'],profile_and_DE_arrays_exactly_replayed=True,reach_verified_with_fresh_raw_p_values=True));del count;gc.collect();print('V'+str(version),'fresh count/reach replay',name,offset,flush=True)
 write(m.A/'count_replay_audit.json',dict(passed=True,checks=checks,replayed_populations=len(checks),auditor_sha256=sha(__file__),scope='Ten prespecified first-target populations: all five contexts at both technical seeds. All stored primary profile/DE arrays replay exactly; independent Python-sorted raw-p tie ordering verifies reach. Other cases get serialized fidelity/Jaccard/NMAE/PDS audits.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--version',type=int,choices=[65,66],required=True);a=p.parse_args()
 with threadpool_limits(limits=2):main(a.version)
