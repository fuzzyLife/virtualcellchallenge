"""V84 fixed public mean test and donor-only-control inference; no hidden upload."""
import argparse,gc,json
import numpy as np,pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits
import v84_intervention_operator as m
from v57_pie import read,write,sha
from audit_transfer_validation import pds
from evaluate_v77 import profiles

def prepare_public():
 from run_v35_v36 import contexts,key
 dest=m.A/'predictions';dest.mkdir(exist_ok=True)
 for ctx in contexts():
  name=ctx['family'];out=dest/(name+'.npz')
  if out.exists():assert sha(out)==read(out.with_suffix('.json'))['sha256'];continue
  p=m.ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
  with np.load(p,allow_pickle=False) as z:q=z['query'].copy();bulk=z['bulk'].copy();base=z['assay75'].copy()
  pred,tr=m.predict(q,base,bulk,(key(name),));fake,ft=m.predict(q,base,bulk,(key(name),),True);np.savez_compressed(out,query=q,bulk=bulk,v66=base,**pred,**{a+'_shuffled':v for a,v in fake.items()});write(out.with_suffix('.json'),dict(sha256=sha(out),traces=tr,shuffled_traces=ft,V66_source_sha256=sha(p),code_sha256=sha(__file__)));print('V84 predictions',name,{k:t['changed_targets'] for k,t in tr.items()},flush=True);gc.collect()

def means():
 from run_v35_v36 import contexts,key
 d=m.data();genes=d['genes'];targets=d['targets'];own=np.array([list(genes).index(t) for t in targets]);arms=('v66',)+m.ARMS+tuple(a+'_shuffled' for a in m.ARMS);records=[]
 for ctx in contexts():
  name=ctx['family'];path=m.A/'predictions'/(name+'.npz');assert sha(path)==read(path.with_suffix('.json'))['sha256'];truthfile=m.ROOT/'experiments/v77_common_response/analysis/truth_bank'/(name+'.npz');assert sha(truthfile)==read(truthfile.with_suffix('.json'))['sha256']
  with np.load(truthfile,allow_pickle=False) as z:bank={k:z[k].copy() for k in z.files}
  with np.load(path,allow_pickle=False) as z:bulk=z['bulk'].copy();pred={k:profiles(z[k],bulk,own)[bank['target_indices']] for k in arms}
  mask=bank['measured']&~np.isin(genes,targets);scores={}
  for arm in arms:
   errors=[]
   for i,ti in enumerate(bank['target_indices']):
    axis=bank['measured'].copy();axis[own[ti]]=False;errors.append(float(np.mean((pred[arm][i,axis]-bank['truth'][i,axis])**2)))
   scores[arm]=dict(profile_MSE=float(np.mean(errors)),pds=pds(pred[arm][:,mask],bank['truth'][:,mask],bank['control'][mask]))
  records.append(dict(source=name,study=key(name),target_context_pairs=len(bank['target_indices']),scores=scores));print('V84 public mean',name,{k:round(s['profile_MSE']/scores['v66']['profile_MSE'],6) for k,s in scores.items()},flush=True)
 studies={};summary={};checks={}
 for arm in arms[1:]:
  studies[arm]={}
  for group in ['H1','Orion','CD4']:
   rows=[r for r in records if r['study']==group];ratio=np.mean([r['scores'][arm]['profile_MSE']/r['scores']['v66']['profile_MSE'] for r in rows]);delta=np.mean([r['scores'][arm]['pds']-r['scores']['v66']['pds'] for r in rows]);studies[arm][group]=dict(profile_ratio=float(ratio),pds_change=float(delta))
  summary[arm]={k:float(np.mean([r[k] for r in studies[arm].values()])) for k in ['profile_ratio','pds_change']}
 for arm in m.ARMS:
  checks[arm]=dict(profile=summary[arm]['profile_ratio']<=.995,worst_study=max(r['profile_ratio'] for r in studies[arm].values())<=1.03,pds=summary[arm]['pds_change']>=-.005,beats_shuffled=summary[arm]['profile_ratio']<summary[arm+'_shuffled']['profile_ratio'])
 eligible=[a for a in m.ARMS if all(checks[a].values())];selected=min(eligible,key=lambda a:summary[a]['profile_ratio']) if eligible else None;write(m.A/'mean_results.json',dict(records=records,studies=studies,summary=summary,checks=checks,eligible=eligible,selected=selected,scope='655 reused public pairs; recipient studies excluded from legacy response estimation. New paired operator uses only Nadig/Replogle. No hosted-score forecast.'));assert not (m.A/'selection.json').exists();write(m.A/'selection.json',dict(selected=selected,primary_always_donor_tested='program25',mean_results_sha256=sha(m.A/'mean_results.json'),strict_development_passed=bool(eligible)));print(json.dumps(dict(summary=summary,selected=selected,checks=checks)),flush=True)

def donor_predict():
 from v6_methods import moments
 from v25_methods import baseline
 import v66_assay_strength as v66
 folder=m.ROOT/'experiments/v83_context_expansion/analysis/donor';receipt=read(folder/'extraction.json');assert sha(folder/'counts.npz')==receipt['counts_sha256'];cohort=read(folder/'cohort.json');assert sha(folder/'cohort.json')==receipt['cohort_sha256'];meta=pd.read_csv(folder/'cells.csv');left=np.array(cohort['left']);assert (meta.gene_target.to_numpy(str)[left]=='Non-Targeting').all();x=sparse.load_npz(folder/'counts.npz');genes=m.data()['genes'];axis=pd.read_csv(folder/'genes.csv');np.testing.assert_array_equal(axis.gene_name.to_numpy(str),genes);_,cp,bulk,_=moments(x[left]);q=cp*1e6;q[~axis.measured.to_numpy(bool)]=np.nan
 # The CD4-excluded baseline branch never loads merged CD4 response teacher.
 parent=baseline(m.data(),q,('CD4',));base=v66.predict(q,parent,bulk,('CD4',))[0]['assay75'];pred,tr=m.predict(q,base,bulk,('CD4',));fake,ft=m.predict(q,base,bulk,('CD4',),True);out=m.A/'donor_predictions.npz';assert not out.exists();np.savez_compressed(out,genes=genes,targets=m.data()['targets'],query=q,cp=cp,bulk=bulk,v24=parent,v66=base,**pred,**{a+'_shuffled':v for a,v in fake.items()});write(out.with_suffix('.json'),dict(sha256=sha(out),extraction_sha256=sha(folder/'extraction.json'),cohort_sha256=sha(folder/'cohort.json'),traces=tr,shuffled_traces=ft,all_CD4_response_evidence_excluded=True,query_fit_NTCs=len(left),query_treatment_outcomes_used=False,code_sha256=sha(__file__)));print('V84 donor responses predicted from NTC only; all CD4 response data excluded',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','means','donor']);a=p.parse_args()
 with threadpool_limits(limits=2):{'prepare':prepare_public,'means':means,'donor':donor_predict}[a.stage]()
