"""Challenge the locked V70 candidate on all eligible existing source responses."""
import gc,json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v70_count_teacher as m
from v25_methods import raw
from v57_pie import read,write,sha
from audit_transfer_validation import profile,pds

def main():
 m.check();plan=read(m.A/'broad_plan.json');assert sha(__file__)==plan['evaluator_sha256'];selected=read(m.A/'selection.json')['selected'];assert selected
 d=m.data();genes=d['genes'];targets=d['targets'];records=[]
 for ctx in m.contexts():
  name=ctx['family'];dest=m.A/'broad'/(name+'.json');dest.parent.mkdir(exist_ok=True)
  if dest.exists():records.append(read(dest));continue
  folder,_=m.checked_folder(name)
  with np.load(m.A/'predictions'/(name+'.npz'),allow_pickle=False) as z:
   pred={a:z[a].copy() for a in ('v66',selected,selected+'_shuffled')};bulk=z['bulk'].copy();measured=z['measured'].copy()
  x,meta,mask,source=raw(name,genes);labels=meta.gene_target.to_numpy(str)
  chosen=[(ti,t,np.flatnonzero(labels==t)) for ti,t in enumerate(targets) if np.count_nonzero(labels==t)>=30]
  truth=np.asarray([profile(x[ids]) for ti,t,ids in chosen]);ref=profile(sparse.load_npz(folder/'reference.npz'));axis=measured&mask&~np.isin(genes,targets);scores={}
  for arm,eff in pred.items():
   pp=[];errors=[]
   for j,(ti,t,ids) in enumerate(chosen):
    own=list(genes).index(t);factor=np.exp2(eff[ti]);factor[own]=.1;p=bulk*factor;p/=p.sum();pr=np.log1p(50000*p);pp.append(pr)
    scored=measured&mask&(np.arange(len(genes))!=own);errors.append(float(np.mean((pr[scored]-truth[j,scored])**2)))
   pp=np.asarray(pp);scores[arm]=dict(pds=pds(pp[:,axis],truth[:,axis],ref[axis]),profile_MSE=float(np.mean(errors)))
  r=dict(source=name,study=m.study(name),targets=len(chosen),scores=scores,target_names=[str(t) for _,t,_ in chosen],actual_cells=[len(ids) for _,_,ids in chosen],source_record=source)
  write(dest,r);records.append(r);print('V70 broader means',name,len(chosen),scores,flush=True);del x;gc.collect()
 studies={}
 for s in ('H1','Orion','CD4'):
  rr=[r for r in records if r['study']==s]
  studies[s]=dict(profile_ratio=float(np.mean([r['scores'][selected]['profile_MSE']/r['scores']['v66']['profile_MSE'] for r in rr])),pds_change=float(np.mean([r['scores'][selected]['pds']-r['scores']['v66']['pds'] for r in rr])),shuffle_ratio=float(np.mean([r['scores'][selected+'_shuffled']['profile_MSE']/r['scores']['v66']['profile_MSE'] for r in rr])))
 ratio=float(np.mean([r['profile_ratio'] for r in studies.values()]));pds_delta=float(np.mean([r['pds_change'] for r in studies.values()]))
 checks=dict(profile=ratio<1,worst_study=max(r['profile_ratio'] for r in studies.values())<=1.02,pds=pds_delta>=-.005,beats_shuffle=ratio<float(np.mean([r['shuffle_ratio'] for r in studies.values()])))
 write(m.A/'broader_assessment.json',dict(selected=selected,studies=studies,profile_ratio=ratio,pds_change=pds_delta,checks=checks,passes=all(checks.values()),target_context_pairs=sum(r['targets'] for r in records),scope='All existing eligible >=30-cell source responses; repeated public biology, intended means only. No strength reselection or final-test guarantee.'))
 print(json.dumps(dict(studies=studies,ratio=ratio,pds_change=pds_delta,checks=checks)),flush=True)

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
