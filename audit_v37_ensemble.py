"""Independent replay of saved ensemble effects, DE sets and profile arithmetic."""
import json
import numpy as np
import v35_methods as m
from run_v37_ensemble import A,CACHE,check,read,write
from audit_v35_v36_counts import equal

def main():
 check();data=m.old.old.read_cache();genes=data['genes'];targets=data['targets'];effect_checks=[];de_checks=0;panels_checked=0
 for folder in sorted(CACHE.iterdir()):
  rec=read(folder/'receipt.json');assert rec['query_treatments_used'] is False
  for parent in (m.A,m.P36):
   model=read(parent/'models'/(rec['excluded_entire_family']+'.json'));assert model['outer']==rec['excluded_entire_family'];assert all(m.family(source)!=rec['excluded_entire_family'] for source in model['sources'])
  left=rec['source_control_rows'];right=rec['reference_rows'];assert len(left)==len(set(left)) and len(right)==len(set(right));assert not set(left)&set(right)
  for cell in rec['cases']:
   assert len(cell['template'])==len(set(cell['template'])) and len(cell['actual'])==len(set(cell['actual']))
   assert set(cell['template']).issubset(left);assert not set(cell['actual'])&(set(left)|set(right))
  for name,h in rec['files'].items():assert m.sha(folder/name)==h
  with np.load(folder/'effects.npz',allow_pickle=False) as z:
   b=z['v24'];r=z['retrieval'];d=r-b
   equal_error=float(np.max(abs(z['mean50']-(b+.5*d))))
   assert equal_error<1e-7
   np.testing.assert_array_equal(z['mean50'],z['population50']);np.testing.assert_array_equal(z['mean50'],z['mean50_noise_mix'])
   c=z['consensus']-b;assert np.max(abs(c))<=np.max(abs(d))+3e-7
   assert not ((c*d)<-1e-8).any()
   for arm in ('v24','mean50','consensus','center_contrast50','population50','mean50_noise_mix'):
    matrix=z[arm]
    assert np.isfinite(matrix).all() and abs(matrix).max()<=3.000001
    for i,t in enumerate(targets):assert np.all(matrix[i,genes==t]==0)
   rms={a:float(np.sqrt(np.mean((z[a]-b)**2))) for a in ('mean50','consensus','center_contrast50')}
   changed={a:int(np.sum(np.any(z[a]!=b,axis=1))) for a in rms}
  effect_checks.append({'family':rec['family'],'effects_finite_bounded_own_masked':True,'mean50_RMS_correction':rms,'changed_targets':changed})
 for splitfile in sorted((A/'null_corrected').glob('*_source.json')):
  split=read(splitfile)['split'];left=set(split['fit_controls']);right=set(split['reference']);template=set(split['template']);actual=set(split['actual'])
  assert not left&right and not left&actual and not right&actual and template.issubset(left)
  assert split==next(ctx['null'] for ctx in __import__('run_v37_ensemble').contexts() if splitfile.name==ctx['family']+'_source.json')
 for stage in ('null','development','confirmation'):
  resultpath=A/(stage+'_results.json')
  if not resultpath.exists():continue
  result=read(resultpath)
  for r in result['records']:
   arrays=A/stage/(f"{r['family']}__{r['case_index']}__{r['seed_offset']}__{r['arm']}.npz")
   if 'arrays_path' in r:arrays=A/r['arrays_path']
   assert m.sha(arrays)==r['arrays_sha256'];assert r['row_libraries_exact']
   with np.load(arrays,allow_pickle=False) as z:
    valid=z['truth_eligible'].copy();valid[genes==r['target']]=False;tq=z['truth_q'];ty=z['truth_lfc'];pq=z['pred_q'];py=z['pred_lfc'];pe=z['pred_eligible'];equal(float(np.sum(pe&(pq<.05))),r['null_DE_calls']);real=valid&(tq<.05);called=valid&(ty!=0)&np.isfinite(ty)&(pq<.05);correct=called&(np.sign(py)==np.sign(ty));budget=max(int(real.sum()),int(called.sum()));fid=correct.sum()/budget if budget else np.nan;union=real|(valid&(pq<.05));jac=(real&(pq<.05)).sum()/union.sum() if union.any() else 1.;den=abs(ty[real]).sum();nmae=abs(np.nan_to_num(py[real])-ty[real]).sum()/den if real.sum()>=10 and den>0 else np.nan
    for k,v in [('fidelity',fid),('jaccard',jac),('nmae',nmae)]:equal(v,r[k])
   de_checks+=1
  for panel in result['panels']:
   rr=sorted([r for r in result['records'] if (r['family'],r['arm'],r['seed_offset'])==(panel['family'],panel['arm'],panel['seed_offset'])],key=lambda r:r['case_index']);values=[]
   for r in rr:
    file=A/stage/(f"{r['family']}__{r['case_index']}__{r['seed_offset']}__{r['arm']}.npz")
    with np.load(file,allow_pickle=False) as z:values.append({k:z[k].copy() for k in ('pred','jack','truth','truth_jack','control','control_jack')})
   pred=np.array([v['pred'] for v in values]);truth=np.array([v['truth'] for v in values]);pj=np.array([v['jack'] for v in values]);tj=np.array([v['truth_jack'] for v in values]);control=values[0]['control'];cj=values[0]['control_jack'];nonpanel=~np.isin(genes,targets);aa=(pred-control)[:,nonpanel];bb=(truth-control)[:,nonpanel];norm=np.linalg.norm(aa,axis=1)[:,None]*np.linalg.norm(bb,axis=1)[None,:];dist=1-np.divide(aa@bb.T,norm,out=np.zeros((len(aa),len(aa))),where=norm>0);diag=np.diag(dist)[:,None];rank=(dist<diag).sum(1)+.5*((dist==diag).sum(1)-1);equal(float(np.mean(1-rank/(len(pred)-1))),panel['pds'])
   own=[np.flatnonzero(genes==r['target'])[0] for r in rr];mask=np.ones(pred.shape,bool);mask[np.arange(len(pred)),own]=False;size=mask.sum(1);w=mask/size[:,None];mu=np.divide((pred*w).sum(0),w.sum(0),out=np.zeros(len(genes)),where=w.sum(0)>0);spread=np.sum((pred-mu)**2*w);cap=np.minimum(pj.sum(1)/size,tj.sum(1)/size);rho=min(1.,spread/cap.sum()) if cap.sum()>0 else 1.;num=float(np.sum(np.sum((pred-truth)**2*mask,1)/size-rho*cap-tj.sum(1)/size));den=float(np.sum(np.sum((truth-control)**2*mask,1)/size-tj.sum(1)/size-cj.sum()/size));equal(num,panel['common_expression_numerator']);equal(den,panel['common_expression_denominator']);equal(num/den if den>0 else np.nan,panel['common_expression_ratio']);panels_checked+=1
 write(A/'independent_audit.json',{'passed':True,'effect_checks':effect_checks,'DE_set_comparisons':de_checks,'profile_panels':panels_checked,'raw_Wilcoxon_and_reach_recomputed':False,'scope':'All saved arrays/checksums replayed; effects and raw-control/treatment identity boundaries checked. Public reused biology, not independent private validation.','auditor_sha256':m.sha(m.ROOT/'audit_v37_ensemble.py')});print('Independent V37 audit passed',de_checks,panels_checked,flush=True)

if __name__=='__main__':main()
