"""Independently verify final production effects against the actual V24 parent."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
A=ROOT/'experiments/v36_population_views/analysis'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 rec=json.loads((A/'reference_model.json').read_text());assert sha(A/'reference_model_arrays.npz')==rec['arrays_sha256'];parent=ROOT/'experiments/v24_cd4_response/analysis/reference_model_arrays.npz';assert sha(parent)==rec['parent_V24_arrays_sha256'];selection=json.loads((A/'selection.json').read_text());assert selection['selected_arm']==rec['selected_arm'];rows=[]
 with np.load(A/'reference_model_arrays.npz',allow_pickle=False) as new,np.load(parent,allow_pickle=False) as old:
  for k in ['genes','targets']:np.testing.assert_array_equal(new[k],old[k])
  genes=new['genes'];targets=new['targets'];own=np.array([np.flatnonzero(genes==t)[0] for t in targets]);effects=new['effects'];base=old['effects'];assert effects.shape==(3,300,18533) and np.isfinite(effects).all() and abs(effects).max()<=3;assert (effects[:,np.arange(300),own]==0).all()
  with np.load(ROOT/'experiments/v35_response_retrieval/analysis/models/all.npz',allow_pickle=False) as model:axis=model['axis'];outside=np.ones(len(genes),bool);outside[axis]=False
  np.testing.assert_array_equal(effects[:,:,outside],base[:,:,outside]);cap=.125002 if selection['selected_arm'] in ['target_centered25','retrieval_contrast25'] else .500002;assert abs(effects-base).max()<=cap
  for ci,c in enumerate(['A','B','C']):
   d=effects[ci].astype(float)-base[ci];rows.append({'context':c,'changed_targets':int((abs(d).sum(1)>0).sum()),'changed_coordinates':int(np.count_nonzero(d)),'RMS_log2FC_correction':float(np.sqrt(np.mean(d*d))),'maximum_log2FC_correction':float(abs(d).max()),'unchanged_unmodeled_genes':int(outside.sum())})
  assert np.count_nonzero(effects-base)>0,'Do not label an identical V24 payload as a new model'
 (A/'production_effect_audit.json').write_text(json.dumps({'passed':True,'model_sha256':sha(A/'reference_model.json'),'arrays_sha256':rec['arrays_sha256'],'records':rows,'checks':'Actual V24 parent; official gene/target ordering; finite bounds; own exclusion; exact fallback outside modeled genes; genuine changed predictions.','auditor_sha256':sha(Path(__file__))},indent=2)+'\n',encoding='utf-8');print(json.dumps(rows,indent=2),flush=True)
if __name__=='__main__':main()