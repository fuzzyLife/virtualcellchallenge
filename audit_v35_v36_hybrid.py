"""Independent DE-set/profile arithmetic audit of the exploratory hybrid arm."""
import json,numpy as np
import v35_methods as m
from run_v35_v36_hybrid import A
from audit_v35_v36_counts import equal

def main():
 result=json.loads((A/'results.json').read_text());file=A/'audit_arrays.npz';assert m.sha(file)==result['audit_arrays_sha256'];data=m.old.old.read_cache();genes=data['genes'];targets=data['targets'];records=[]
 with np.load(file,allow_pickle=False) as z:
  for r in result['records']:
   pre=r['family']+'__'+r['target']+'__';valid=z[pre+'eligible'].copy();valid[genes==r['target']]=False;tq=z[pre+'truth_q'];ty=z[pre+'truth_lfc'];pq=z[pre+'pred_q'];py=z[pre+'pred_lfc'];real=valid&(tq<.05);called=valid&(ty!=0)&np.isfinite(ty)&(pq<.05);correct=called&(np.sign(py)==np.sign(ty));budget=max(int(real.sum()),int(called.sum()));union=real|(valid&(pq<.05));fid=correct.sum()/budget if budget else np.nan;jac=(real&(pq<.05)).sum()/union.sum() if union.any() else 1.;den=abs(ty[real]).sum();nmae=abs(py[real]-ty[real]).sum()/den if real.sum()>=10 and den>0 else np.nan
   for k,v in [('fidelity',fid),('jaccard',jac),('nmae',nmae)]:equal(v,r[k])
   records.append({'family':r['family'],'target':r['target'],'passed':True})
  for r in result['panels']:
   pre=r['family']+'__';pred=z[pre+'pred'];truth=z[pre+'truth'];ctrl=z[pre+'control'];pj=z[pre+'jack'];tj=z[pre+'truth_jack'];cj=z[pre+'control_jack'];cols=~np.isin(genes,targets);a=(pred-ctrl)[:,cols];b=(truth-ctrl)[:,cols];norm=np.linalg.norm(a,axis=1)[:,None]*np.linalg.norm(b,axis=1)[None,:];dist=1-np.divide(a@b.T,norm,out=np.zeros((len(a),len(a))),where=norm>0);diag=np.diag(dist)[:,None];rank=(dist<diag).sum(1)+.5*((dist==diag).sum(1)-1);equal(float(np.mean(1-rank/(len(a)-1))),r['pds']);names=[v['target'] for v in result['records'] if v['family']==r['family']];own=[np.flatnonzero(genes==t)[0] for t in names];mask=np.ones(pred.shape,bool);mask[np.arange(len(pred)),own]=False;size=mask.sum(1);w=mask/size[:,None];mu=np.divide((pred*w).sum(0),w.sum(0),out=np.zeros(len(genes)),where=w.sum(0)>0);cap=np.minimum(pj.sum(1)/size,tj.sum(1)/size);spread=np.sum((pred-mu)**2*w);rho=min(1.,spread/cap.sum()) if cap.sum()>0 else 1.;num=np.sum(np.sum((pred-truth)**2*mask,1)/size-rho*cap-tj.sum(1)/size);den=np.sum(np.sum((truth-ctrl)**2*mask,1)/size-tj.sum(1)/size-cj.sum()/size);equal(num,r['common_expression_numerator']);equal(den,r['common_expression_denominator']);equal(num/den if den>0 else np.nan,r['common_expression_ratio'])
 assert len(records)==20 and len(result['panels'])==5
 m.write(A/'independent_audit.json',{'passed':True,'results_sha256':m.sha(A/'results.json'),'records':records,'checks':'All20 DE-set FID/JAC/NMAE summaries and5 PDS/expression panels independently replayed from saved arrays. Raw Wilcoxon/reach p-value tie order inherited, not independently regenerated. Exploratory reused cases, no blind confirmation.','auditor_sha256':m.sha(m.ROOT/'audit_v35_v36_hybrid.py')});print('Hybrid independent audit passed',flush=True)
if __name__=='__main__':main()