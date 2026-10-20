"""Joint V77 selection on a larger PDS panel; separate held-target validation."""
import argparse,json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v77_common_response as m
from v57_pie import read,write,sha
from audit_transfer_validation import profile,pds

def profiles(effects,bulk,own):
 factor=np.exp2(effects.astype(float));factor[np.arange(len(own)),own]=.1;p=bulk[None,:]*factor;p/=p.sum(1,keepdims=True)
 return np.log1p(50000*p)

def load(name,arms):
 p=m.A/'predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
 with np.load(p,allow_pickle=False) as z:bulk=z['bulk'].copy();pred={a:z[a].copy() for a in arms}
 p=m.A/'truth_bank'/(name+'.npz');receipt=read(p.with_suffix('.json'));assert sha(p)==receipt['sha256']
 with np.load(p,allow_pickle=False) as z:bank={k:z[k].copy() for k in z.files}
 return pred,bulk,bank,receipt

def development():
 m.check();d=m.data();genes=d['genes'];targets=d['targets'];own=np.array([list(genes).index(t) for t in targets]);arms=('v66',)+m.ARMS+tuple(a+'_shuffled' for a in m.ARMS);records=[];panels=[]
 for ctx in m.contexts():
  name=ctx['family'];pred,bulk,bank,receipt=load(name,arms);lookup={str(t):i for i,t in enumerate(bank['target_names'])};rows=np.array([lookup[t] for t in receipt['development_pds_targets']]);ti=bank['target_indices'][rows];mask=bank['measured']&~np.isin(genes,targets);folder,_=m.checked_folder(name)
  for arm in arms:
   pp=profiles(pred[arm],bulk,own)
   panels.append(dict(source=name,study=m.study(name),arm=arm,targets=len(rows),pds=pds(pp[ti][:,mask],bank['truth'][rows][:,mask],bank['control'][mask])))
   for i in m.cases(ctx)[:3]:
    t=ctx['cells'][i]['target'];index=list(targets).index(t);j=own[index];truth=profile(sparse.load_npz(folder/f'{i}_actual.npz'));axis=np.arange(len(genes))!=j
    records.append(dict(source=name,study=m.study(name),target=t,arm=arm,profile_error=float(np.sum((pp[index,axis]-truth[axis])**2))))
  print('V77 larger development panel',name,len(rows),flush=True)
 studies={};summary={}
 for arm in arms[1:]:
  studies[arm]={}
  for s in ('H1','Orion','CD4'):
   mean=sum(r['profile_error'] for r in records if r['arm']==arm and r['study']==s)/sum(r['profile_error'] for r in records if r['arm']=='v66' and r['study']==s)
   change=np.mean([r['pds'] for r in panels if r['arm']==arm and r['study']==s])-np.mean([r['pds'] for r in panels if r['arm']=='v66' and r['study']==s])
   studies[arm][s]=dict(profile_ratio=float(mean),pds_change=float(change))
  summary[arm]={k:float(np.mean([v[k] for v in studies[arm].values()])) for k in ('profile_ratio','pds_change')}
 eligible=[a for a in m.ARMS if summary[a]['pds_change']>=.01 and summary[a]['profile_ratio']<=1.005 and max(v['profile_ratio'] for v in studies[a].values())<=1.02 and summary[a]['profile_ratio']<summary[a+'_shuffled']['profile_ratio'] and summary[a]['pds_change']>summary[a+'_shuffled']['pds_change']]
 selected=min(eligible,key=lambda a:(-summary[a]['pds_change'],summary[a]['profile_ratio'])) if eligible else None
 write(m.A/'development_results.json',dict(records=records,panels=panels,studies=studies,summary=summary,eligible=eligible,selected=selected,scope='PDS development55 metadata-selected targets; profile first15 fixed cases. Reused public biology.'))
 assert not (m.A/'selection.json').exists();write(m.A/'selection.json',dict(selected=selected,eligible=eligible,development_sha256=sha(m.A/'development_results.json'),rule='Rank-focused joint criteria frozen in plan; prefer PDS, tie lower profile error.'))
 print(json.dumps(dict(selected=selected,summary=summary)),flush=True)

def held():
 m.check();selected=read(m.A/'selection.json')['selected'];assert selected;d=m.data();genes=d['genes'];targets=d['targets'];own=np.array([list(genes).index(t) for t in targets]);arms=('v66',selected,selected+'_shuffled');records=[]
 for ctx in m.contexts():
  name=ctx['family'];pred,bulk,bank,receipt=load(name,arms);keep=~np.isin(bank['target_names'],receipt['development_pds_targets']);ids=bank['target_indices'][keep];truth=bank['truth'][keep];axis=bank['measured']&~np.isin(genes,targets);scores={}
  for arm in arms:
   pp=profiles(pred[arm],bulk,own)[ids];error=[]
   for j,ti in enumerate(ids):
    mask=bank['measured'].copy();mask[own[ti]]=False;error.append(float(np.mean((pp[j,mask]-truth[j,mask])**2)))
   scores[arm]=dict(pds=pds(pp[:,axis],truth[:,axis],bank['control'][axis]),profile_MSE=float(np.mean(error)))
  records.append(dict(source=name,study=m.study(name),targets=len(ids),scores=scores));print('V77 held-target point panel',name,len(ids),scores,flush=True)
 studies={}
 for s in ('H1','Orion','CD4'):
  r=[x for x in records if x['study']==s];studies[s]=dict(profile_ratio=float(np.mean([x['scores'][selected]['profile_MSE']/x['scores']['v66']['profile_MSE'] for x in r])),pds_change=float(np.mean([x['scores'][selected]['pds']-x['scores']['v66']['pds'] for x in r])),shuffle_profile_ratio=float(np.mean([x['scores'][selected+'_shuffled']['profile_MSE']/x['scores']['v66']['profile_MSE'] for x in r])),shuffle_pds_change=float(np.mean([x['scores'][selected+'_shuffled']['pds']-x['scores']['v66']['pds'] for x in r])))
 ratio=float(np.mean([v['profile_ratio'] for v in studies.values()]));delta=float(np.mean([v['pds_change'] for v in studies.values()]))
 checks=dict(profile=ratio<=1.005,worst_study=max(v['profile_ratio'] for v in studies.values())<=1.02,pds=delta>=.005,beats_shuffle=ratio<float(np.mean([v['shuffle_profile_ratio'] for v in studies.values()])) and delta>float(np.mean([v['shuffle_pds_change'] for v in studies.values()])))
 write(m.A/'held_assessment.json',dict(selected=selected,records=records,studies=studies,profile_ratio=ratio,pds_change=delta,checks=checks,passes=all(checks.values()),target_context_pairs=sum(r['targets'] for r in records),scope='Targets excluded from current development-PDS and profile selection; public biological studies and labels have been used previously. Not new independent biology or hidden-score estimate.'))
 print(json.dumps(dict(ratio=ratio,pds_change=delta,checks=checks)),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['development','held']);args=p.parse_args()
 with threadpool_limits(limits=2):(development if args.stage=='development' else held)()
