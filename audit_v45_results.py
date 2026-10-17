"""Independent V38/V39 replay of DE sets, PDS and capped expression arithmetic."""
import argparse,json
import numpy as np
import v38_response_hurdle as m
from run_v35_v36 import contexts
from run_v37_ensemble import checked_folder
from audit_v35_v36_counts import equal
def main(version=45):
 global m
 modules={40:'v40_identity_calibration',41:'v41_identity_contrast',42:'v42_identity_trust',43:'v43_lincs_fingerprints',44:'v44_response_dispersion',45:'v45_bounded_dispersion'}
 if version in modules:
  import importlib
  m=importlib.import_module(modules[version])
 if version==39:
  import v39_kolf_hurdle as extra
  extra.install()
 from run_v38_response_hurdle import check
 (m.check() if version>=40 else check());d=m.data();genes=d['genes'];targets=d['targets'];de_checks=0;panel_checks=0
 for ctx in contexts():
  folder,r=checked_folder(ctx['family']);left=set(r['source_control_rows']);right=set(r['reference_rows']);assert not left&right
  for cell in r['cases']:assert set(cell['template']).issubset(left) and not set(cell['actual'])&(left|right)
  null=ctx['null'];fit=set(null['fit_controls']);ref=set(null['reference']);actual=set(null['actual']);assert not fit&ref and not fit&actual and not ref&actual and set(null['template']).issubset(fit)
  file=m.A/'predictions'/(ctx['family']+'.npz');receipt=m.read(file.with_suffix('.json'));assert m.b.sha256(file)==receipt['sha256']
  with np.load(file,allow_pickle=False) as z:
   for a in ('v24',)+m.ARMS:
    assert np.isfinite(z[a]).all() and abs(z[a]).max()<=3.000001
    for j,t in enumerate(receipt['case_targets']):assert z[a][j,genes==t]==0
 for stage in ['null','development','confirmation']:
  path=m.A/(stage+'_results.json')
  if not path.exists():continue
  result=m.read(path)
  for r in result['records']:
   file=m.A/stage/(f"{r['family']}__{r['case_index']}__{r['seed_offset']}__{r['arm']}.npz");assert m.b.sha256(file)==r['arrays_sha256'];assert r['row_libraries_exact']
   with np.load(file,allow_pickle=False) as z:
    valid=z['truth_eligible'].copy();valid[genes==r['target']]=False;tq=z['truth_q'];ty=z['truth_lfc'];pq=z['pred_q'];py=z['pred_lfc'];pe=z['pred_eligible']
    equal(float(np.sum(pe&(pq<.05))),r['null_DE_calls']);real=valid&(tq<.05);called=valid&(ty!=0)&np.isfinite(ty)&(pq<.05);correct=called&(np.sign(py)==np.sign(ty));budget=max(int(real.sum()),int(called.sum()));fid=correct.sum()/budget if budget else np.nan;union=real|(valid&(pq<.05));jac=(real&(pq<.05)).sum()/union.sum() if union.any() else 1.;den=abs(ty[real]).sum();nmae=abs(np.nan_to_num(py[real])-ty[real]).sum()/den if real.sum()>=10 and den>0 else np.nan
    for k,v in [('fidelity',fid),('jaccard',jac),('nmae',nmae)]:equal(v,r[k])
   de_checks+=1
  for p in result['panels']:
   rr=sorted([r for r in result['records'] if (r['family'],r['arm'],r['seed_offset'])==(p['family'],p['arm'],p['seed_offset'])],key=lambda r:r['case_index']);values=[]
   for r in rr:
    file=m.A/stage/(f"{r['family']}__{r['case_index']}__{r['seed_offset']}__{r['arm']}.npz")
    with np.load(file,allow_pickle=False) as z:values.append({k:z[k].copy() for k in ('pred','jack','truth','truth_jack','control','control_jack')})
   pred=np.array([v['pred'] for v in values]);truth=np.array([v['truth'] for v in values]);pj=np.array([v['jack'] for v in values]);tj=np.array([v['truth_jack'] for v in values]);control=values[0]['control'];cj=values[0]['control_jack'];nonpanel=~np.isin(genes,targets);aa=(pred-control)[:,nonpanel];bb=(truth-control)[:,nonpanel];norm=np.linalg.norm(aa,axis=1)[:,None]*np.linalg.norm(bb,axis=1)[None,:];dist=1-np.divide(aa@bb.T,norm,out=np.zeros((len(aa),len(aa))),where=norm>0);diag=np.diag(dist)[:,None];rank=(dist<diag).sum(1)+.5*((dist==diag).sum(1)-1);equal(float(np.mean(1-rank/(len(pred)-1))),p['pds'])
   own=[np.flatnonzero(genes==r['target'])[0] for r in rr];mask=np.ones(pred.shape,bool);mask[np.arange(len(pred)),own]=False;size=mask.sum(1);w=mask/size[:,None];mu=np.divide((pred*w).sum(0),w.sum(0),out=np.zeros(len(genes)),where=w.sum(0)>0);spread=np.sum((pred-mu)**2*w);cap=np.minimum(pj.sum(1)/size,tj.sum(1)/size);rho=min(1.,spread/cap.sum()) if cap.sum()>0 else 1.;num=float(np.sum(np.sum((pred-truth)**2*mask,1)/size-rho*cap-tj.sum(1)/size));den=float(np.sum(np.sum((truth-control)**2*mask,1)/size-tj.sum(1)/size-cj.sum()/size))
   equal(num,p['common_expression_numerator']);equal(den,p['common_expression_denominator']);equal(num/den if den>0 else np.nan,p['common_expression_ratio']);panel_checks+=1
 m.write(m.A/'independent_audit.json',{'passed':True,'version':version,'DE_set_comparisons':de_checks,'profile_panels':panel_checks,'scoring_identities_and_axis_masks_checked':True,'all_saved_arrays_hash_verified':True,'raw_Wilcoxon_and_reach_recomputed':False,'scope':'Independent set/profile arithmetic replay from sufficient statistics. Public reused biology, not private score validation.','auditor_sha256':m.b.sha256(m.ROOT/'audit_v45_results.py')});print('Independent audit PASS',version,de_checks,panel_checks,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--version',type=int,choices=[45],default=45);a=p.parse_args();main(a.version)



