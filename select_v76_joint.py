"""Tighten V76 development selection with a joint PDS constraint before confirmation."""
import json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v76_response_dilution as m
from audit_transfer_validation import profile,pds
from v57_pie import read,write,sha

def main():
 m.check();means=read(m.A/'development_mean_results.json');original=m.A/'selection.json';assert not (m.A/'joint_selection.json').exists()
 d=m.data();genes=d['genes'];targets=d['targets'];rows=[];arms=('v66',)+m.ARMS
 for ctx in m.contexts():
  name=ctx['family'];folder,_=m.checked_folder(name);ids=m.cases(ctx)[:3];truth=np.asarray([profile(sparse.load_npz(folder/f'{i}_actual.npz')) for i in ids]);ctrl=profile(sparse.load_npz(folder/'reference.npz'));mask=~np.isin(genes,targets)
  with np.load(m.A/'predictions'/(name+'.npz'),allow_pickle=False) as z:bulk=z['bulk'].copy();pred={a:z[a].copy() for a in arms}
  for arm in arms:
   pp=[]
   for i in ids:
    target=ctx['cells'][i]['target'];ti=list(targets).index(target);own=list(genes).index(target);factor=np.exp2(pred[arm][ti].astype(float))
    if arm=='v66':factor[own]=.1
    p=bulk*factor;p/=p.sum();pp.append(np.log1p(50000*p))
   rows.append(dict(source=name,study=m.study(name),arm=arm,pds=pds(np.asarray(pp)[:,mask],truth[:,mask],ctrl[mask])))
  print('V76 joint development PDS',name,flush=True)
 changes={}
 for arm in m.ARMS:
  changes[arm]=float(np.mean([np.mean([r['pds'] for r in rows if r['study']==s and r['arm']==arm])-np.mean([r['pds'] for r in rows if r['study']==s and r['arm']=='v66']) for s in ('H1','Orion','CD4')]))
 eligible=[a for a in m.ARMS if means['summary'][a]<=.995 and max(means['study_ratios'][a].values())<=1.02 and means['summary'][a]<means['summary'][a+'_shuffled'] and changes[a]>=-.005]
 selected=min(eligible,key=lambda a:means['summary'][a]) if eligible else None
 write(m.A/'joint_selection.json',dict(selected=selected,eligible=eligible,pds_changes=changes,records=rows,development_sha256=sha(m.A/'development_mean_results.json'),criterion='All original development mean/worst-study/shuffle criteria plus PDS >= baseline-.005, tightened before new-model confirmation. Original MSE-only ranking retained separately.',scope='Repeated public development, three targets/context; broad test remains mandatory.'))
 backup=m.A/'mse_only_preselection.json';assert not backup.exists();backup.write_bytes(original.read_bytes());write(original,read(m.A/'joint_selection.json'))
 print(json.dumps(dict(selected=selected,eligible=eligible,pds_changes=changes)),flush=True)

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
