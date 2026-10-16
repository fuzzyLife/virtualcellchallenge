"""Independent arithmetic replay from saved public count-test sufficient arrays."""
import json,numpy as np
import v35_methods as m

def equal(actual,reported):
 if reported is None:assert not np.isfinite(actual)
 else:np.testing.assert_allclose(actual,reported,rtol=1e-9,atol=1e-10)

def main():
 data=m.old.old.read_cache();genes=data['genes'];targets=data['targets'];checks=[]
 for stage in ['null','treatment']:
  result=json.loads((m.P36/(stage+'_results.json')).read_text());file=m.P36/(stage+'_audit_arrays.npz');assert m.sha(file)==result['audit_arrays_sha256']
  with np.load(file,allow_pickle=False) as z:
   for r in result['records']:
    prefix=r['family']+'__'+r['target']+'__';tq=z[prefix+'truth_q'];ty=z[prefix+'truth_lfc'];valid=z[prefix+'truth_eligible'].copy();own=np.flatnonzero(genes==r['target']);valid[own]=False;pq=z[prefix+r['arm']+'_q'];py=z[prefix+r['arm']+'_lfc'];eligible=z[prefix+r['arm']+'_eligible'];equal(float(np.sum((pq<.05)&eligible)),r['null_DE_calls']);real=valid&(tq<.05);adj=valid&np.isfinite(ty)&(ty!=0);called=adj&(pq<.05);correct=called&(np.sign(py)==np.sign(ty));budget=max(int(real.sum()),int(called.sum()));fid=correct.sum()/budget if budget else np.nan;union=real|(valid&(pq<.05));jac=np.sum(real&(pq<.05))/union.sum() if union.any() else 1.;den=abs(ty[real]).sum();nmae=abs(np.nan_to_num(py[real])-ty[real]).sum()/den if real.sum()>=10 and den>0 else np.nan
    for k,v in [('fidelity',fid),('jaccard',jac),('nmae',nmae)]:equal(v,r[k])
    checks.append({'stage':stage,'family':r['family'],'target':r['target'],'arm':r['arm'],'DE_set_metrics_replayed':True})
   if stage=='treatment':
    for r in result['panels']:
     pre=r['family']+'__'+r['arm']+'__';pred=z[pre+'pred'];truth=z[pre+'truth'];control=z[pre+'control'];pj=z[pre+'jack'];tj=z[pre+'truth_jack'];cj=z[pre+'control_jack'];nonpanel=~np.isin(genes,targets);aa=(pred-control)[:,nonpanel];bb=(truth-control)[:,nonpanel];cos=(aa@bb.T)/np.maximum(np.linalg.norm(aa,axis=1)[:,None]*np.linalg.norm(bb,axis=1)[None,:],1e-300);dist=1-cos;diag=np.diag(dist)[:,None];rank=(dist<diag).sum(1)+.5*((dist==diag).sum(1)-1);pds=float(np.mean(1-rank/(len(pred)-1)));equal(pds,r['pds'])
     tr=[v['target'] for v in result['records'] if v['family']==r['family'] and v['arm']==r['arm']];own=np.array([np.flatnonzero(genes==t)[0] for t in tr]);mask=np.ones(pred.shape,bool);mask[np.arange(len(pred)),own]=False;sizes=mask.sum(1);weights=mask/sizes[:,None];mean=np.divide((pred*weights).sum(0),weights.sum(0),out=np.zeros(pred.shape[1]),where=weights.sum(0)>0);spread=float(np.sum((pred-mean)**2*weights));pc=pj.sum(1)/sizes;tc=tj.sum(1)/sizes;cc=cj.sum()/sizes;cap=np.minimum(pc,tc);rho=min(1.,spread/cap.sum()) if cap.sum()>0 else 1.;num=float(np.sum(np.sum((pred-truth)**2*mask,1)/sizes-rho*cap-tc));den=float(np.sum(np.sum((truth-control)**2*mask,1)/sizes-tc-cc));equal(num,r['common_expression_numerator']);equal(den,r['common_expression_denominator']);equal(num/den if den>0 else np.nan,r['common_expression_ratio'])
 m.write(m.P36/'independent_count_audit.json',{'passed':True,'comparisons':len(checks),'records':checks,'scope':'Reimplemented all100 saved DE-set FID/JAC/NMAE/null-call summaries and20 treatment PDS/expression panels. Raw Wilcoxon p-values and reach tie order were not independently regenerated; integer/libraries asserted during generation. These are reused public diagnostics, not private VCC anchors.','auditor_sha256':m.sha(m.ROOT/'audit_v35_v36_counts.py')});print('Independent count audit passed',len(checks),flush=True)
if __name__=='__main__':main()