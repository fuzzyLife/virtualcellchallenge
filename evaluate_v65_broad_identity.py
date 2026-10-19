"""Broader response diagnostic for both fixed finalists; never used for reselection."""
import gc,json
from pathlib import Path
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key
from run_v37_ensemble import checked_folder
from v25_methods import raw
from audit_transfer_validation import profile,pds
import v62_conserved_response as old
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v65_closed_response/analysis'

def main():
 d=old.data();genes=d['genes'];targets=d['targets'];records=[]
 for ctx in contexts():
  name=ctx['family'];folder,_=checked_folder(name);pred={}
  for version,dirname in [(65,'v65_closed_response'),(66,'v66_assay_strength')]:
   a=ROOT/'experiments'/dirname/'analysis';mean=read(a/'mean_results.json');candidate=mean['selected'];assert candidate
   path=a/'predictions'/(name+'.npz');assert sha(path)==read(path.with_suffix('.json'))['sha256']
   with np.load(path,allow_pickle=False) as z:pred['v'+str(version)]=z[candidate].copy();pred['v63']=z['v63'].copy();bulk=z['bulk'].copy();measured=z['measured'].copy()
  x,meta,mask,source=raw(name,genes);labels=meta.gene_target.to_numpy(str);chosen=[(ti,t,np.flatnonzero(labels==t)) for ti,t in enumerate(targets) if np.count_nonzero(labels==t)>=30];truth=np.asarray([profile(x[ids]) for ti,t,ids in chosen]);ref=profile(sparse.load_npz(folder/'reference.npz'));axis=measured&mask&~np.isin(genes,targets);scores={}
  for arm,eff in pred.items():
   pp=[];errors=[]
   for j,(ti,t,ids) in enumerate(chosen):
    own=list(genes).index(t);factor=np.exp2(eff[ti]);factor[own]=.1;prob=bulk*factor;prob/=prob.sum();p=np.log1p(50000*prob);pp.append(p);scored=measured&mask&(np.arange(len(genes))!=own);errors.append(float(np.mean((p[scored]-truth[j,scored])**2)))
   pp=np.asarray(pp);scores[arm]=dict(pds=pds(pp[:,axis],truth[:,axis],ref[axis]),profile_MSE=float(np.mean(errors)))
  records.append(dict(source=name,study=key(name),targets=len(chosen),actual_cell_counts=[len(ids) for ti,t,ids in chosen],target_names=[str(t) for ti,t,ids in chosen],scores=scores));print('Broader finalist diagnostic',name,len(chosen),scores,flush=True);del x;gc.collect()
 summary={}
 for arm in ('v65','v66'):
  studies={s:dict(profile_ratio=float(np.mean([r['scores'][arm]['profile_MSE']/r['scores']['v63']['profile_MSE'] for r in records if r['study']==s])),pds_change=float(np.mean([r['scores'][arm]['pds']-r['scores']['v63']['pds'] for r in records if r['study']==s]))) for s in ('H1','Orion','CD4')}
  summary[arm]=dict(studies=studies,profile_ratio=float(np.mean([r['profile_ratio'] for r in studies.values()])),pds_change=float(np.mean([r['pds_change'] for r in studies.values()])))
 write(A/'broader_finalist_diagnostic.json',dict(records=records,summary=summary,evaluator_sha256=sha(__file__),scope='All available existing-extraction targets with >=30 genuine treated cells (655 target/context pairs expected), intended pooled prediction means before count emission. Profile errors exclude own and unmeasured genes; PDS excludes all panel targets. Reused public biology; diagnostic only, no reselection or hidden-score estimate.'))
 print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
