"""Build a frozen public profile bank and a metadata-selected development PDS panel."""
import hashlib,gc,json
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v77_common_response as m
from v25_methods import raw
from v57_pie import read,write,sha
from audit_transfer_validation import profile

def main():
 m.check();d=m.data();genes=d['genes'];targets=d['targets'];folder=m.A/'truth_bank';folder.mkdir(exist_ok=True)
 for ctx in m.contexts():
  name=ctx['family'];p=folder/(name+'.npz')
  if p.exists():assert sha(p)==read(p.with_suffix('.json'))['sha256'];continue
  x,meta,measured,source=raw(name,genes);label=meta.gene_target.to_numpy(str)
  chosen=[(ti,str(t),np.flatnonzero(label==t)) for ti,t in enumerate(targets) if np.count_nonzero(label==t)>=30]
  fixed={c['target'] for c in ctx['cells']};extra=[t for ti,t,ids in chosen if len(ids)>=100 and t not in fixed];extra.sort(key=lambda t:hashlib.sha256(('v77-pds:'+t).encode()).hexdigest());extra=extra[:8]
  dev=[ctx['cells'][i]['target'] for i in m.cases(ctx)[:3]]+extra
  cohort=p.with_suffix('.cohort.json');assert not cohort.exists();write(cohort,dict(source=name,development_pds_targets=dev,extra_targets=extra,all_targets=[t for _,t,_ in chosen],cell_counts=[len(ids) for _,_,ids in chosen],selection='Metadata only; cohort written before any truth profile is computed.',source_record=source))
  qfolder,_=m.checked_folder(name);control=profile(sparse.load_npz(qfolder/'reference.npz'))
  truth=np.asarray([profile(x[ids]) for ti,t,ids in chosen]);np.savez(p,genes=genes,target_names=np.asarray([t for _,t,_ in chosen],str),target_indices=np.asarray([ti for ti,_,_ in chosen]),truth=truth,control=control,measured=measured,cell_counts=np.asarray([len(ids) for _,_,ids in chosen]))
  write(p.with_suffix('.json'),dict(sha256=sha(p),cohort_sha256=sha(cohort),source=name,targets=len(chosen),development_pds_targets=dev,script_sha256=sha(__file__),source_record=source))
  print('V77 truth bank',name,len(chosen),'development PDS',len(dev),flush=True);del x;gc.collect()

if __name__=='__main__':
 with threadpool_limits(limits=2):main()
