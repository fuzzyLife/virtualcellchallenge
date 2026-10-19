"""Independently audit V66 serialized DE sets and unique-target profile panels."""
import math,gc
import numpy as np
from v57_pie import sha,read,write
import v66_assay_strength as m
from run_v65_experiment import cases
from run_v35_v36 import contexts,key

def close(a,b):
 if b is None or not np.isfinite(b):assert a is None or not np.isfinite(a)
 else:assert abs(a-b)<1e-10,(a,b)

def main():
 m.check();assessment=read(m.A/'count_assessment.json');assert assessment['complete'];selected=assessment['selected'];d=m.old.data();genes=d['genes'];targets=d['targets'];arms=('v24','v63',selected,selected+'_shuffled');panels=0;cases_checked=0;reuse=0
 for ctx in contexts():
  name=ctx['family'];ix=cases(ctx);names=[ctx['cells'][i]['target'] for i in ix];assert len(names)==len(set(names))
  for seed in (0,620000):
   for arm in arms:
    pp=[];tt=[];errors=[]
    for i,target in zip(ix,names):
     path=m.A/'counts'/f'{name}__{i}__{seed}__{arm}.json';r=read(path);arr=path.with_suffix('.npz');assert sha(arr)==r['arrays_sha256'];assert r['family']==name and r['target']==target and r['case_index']==i and r['seed_offset']==seed and r['arm']==arm
     if 'reused_source' in r:
      src=m.ROOT/r['reused_source'];assert sha(src)==r['reused_receipt_sha256'];old=read(src);assert sha(src.with_suffix('.npz'))==old['arrays_sha256']==r['arrays_sha256'];reuse+=1
     own=list(genes).index(target)
     with np.load(arr,allow_pickle=False) as z:
      q=z['pred_q'];l=z['pred_lfc'];tq=z['truth_q'];tl=z['truth_lfc'];valid=z['truth_eligible']&(np.arange(len(genes))!=own);real=valid&(tq<.05);called=valid&np.isfinite(tl)&(tl!=0)&(q<.05);correct=called&(np.sign(l)==np.sign(tl));budget=max(real.sum(),called.sum());fid=correct.sum()/budget if budget else np.nan
      union=real|(valid&(q<.05));jac=np.sum(real&(q<.05))/union.sum() if union.any() else 1.;denom=abs(tl[real]).sum();nmae=abs(np.nan_to_num(l[real])-tl[real]).sum()/denom if real.sum()>=10 and denom>0 and np.isfinite(denom) else np.nan
      close(fid,r['fidelity']);close(jac,r['jaccard']);close(nmae,r['nmae']);assert int(called.sum())==r['predicted_calls'] and int(real.sum())==r['reference_calls']
      pp.append(z['pred'].copy());tt.append(z['truth'].copy());control=z['control'].copy();errors.append(float(np.sum((z['pred'][np.arange(len(genes))!=own]-z['truth'][np.arange(len(genes))!=own])**2)))
      for k in ('pred','truth','control','jack','truth_jack','control_jack'):assert np.isfinite(z[k]).all()
     cases_checked+=1
    p=read(m.A/'counts'/f'{name}__{seed}__{arm}__panel.json');axis=~np.isin(genes,targets);p0=np.asarray(pp)[:,axis]-control[axis];t0=np.asarray(tt)[:,axis]-control[axis];ranks=[]
    for j,row in enumerate(p0):
     distances=[]
     for truthrow in t0:
      norm=np.linalg.norm(row)*np.linalg.norm(truthrow);distances.append(1-float(row@truthrow)/norm if norm else 1.)
     dist=np.asarray(distances);rank=np.sum(dist<dist[j])+.5*(np.sum(dist==dist[j])-1);ranks.append(1-rank/(len(names)-1))
    close(float(np.mean(ranks)),p['pds']);close(float(np.mean(errors)),p['profile_square_error']);panels+=1
 assert cases_checked==240 and panels==40
 write(m.A/'results_audit.json',dict(passed=True,serialized_cases=cases_checked,unique_target_panels=panels,verified_baseline_reuse_receipts=reuse,independently_replayed_metrics=['fidelity','jaccard','nmae','predicted_calls','reference_calls','PDS','pooled_profile_square_error'],reach_replayed_from_serialized_arrays=False,reach_limitation='Legacy arrays do not store raw p values used for tie ordering; no full reach replay claimed here. Production integrity uses fresh count regeneration.',auditor_sha256=sha(__file__),count_assessment_sha256=sha(m.A/'count_assessment.json')))
 print('V66 independent serialized metric audit passed:',cases_checked,'cases,',panels,'unique-target panels;',reuse,'verified reused baselines.')
if __name__=='__main__':main()
