"""Independent saved-array and delete-target expression audit for V55."""
from pathlib import Path
import json,hashlib
import numpy as np
import run_v55_mean50 as run
ROOT=Path(__file__).resolve().parent;A=run.A

def expression(p,t,c,pj,tj,cj,own):
 # Independent masked loop accumulation of the frozen scorer's whole-gene jackknife convention.
 n,g=p.shape;weights=np.ones((n,g),float)
 weights[np.arange(n),own]=0;weights/=weights.sum(1)[:,None]
 totals=weights.sum(0);center=np.divide(np.sum(p*weights,0),totals,out=np.zeros(g),where=totals>0)
 spread=np.sum(weights*np.square(p-center));den=g-1
 cap=np.minimum(np.sum(pj,1)/den,np.sum(tj,1)/den)
 rho=min(1.,float(spread/np.sum(cap))) if np.sum(cap)>0 else 1.
 numerator=0.;denominator=0.
 for i in range(n):
  mask=np.ones(g,bool);mask[own[i]]=False
  numerator+=np.sum(np.square(p[i,mask]-t[i,mask]))/den-rho*cap[i]-np.sum(tj[i])/den
  denominator+=np.sum(np.square(t[i,mask]-c[mask]))/den-np.sum(tj[i])/den-np.sum(cj)/den
 return {'common_expression_numerator':float(numerator),'common_expression_denominator':float(denominator),'common_expression_ratio':float(numerator/denominator) if denominator>0 else None}

def main():
 run.check();d=run.read(A/'evaluation.json');_,genes,targets=run.b.read_spec(ROOT);metrics_checked=0;profiles={};max_error=0
 for r in d['records']:
  file=A/'counts'/f"{r['family']}__{r['case_index']}__{r['seed_offset']}__{r['arm']}.npz"
  assert run.sha(file)==r['arrays_sha256']
  with np.load(file,allow_pickle=False) as z:
   own=int(np.flatnonzero(genes==r['target'])[0]);valid=z['truth_eligible'].copy();valid[own]=False
   real=valid&(z['truth_q']<.05);adjudicable=valid&np.isfinite(z['truth_lfc'])&(z['truth_lfc']!=0);called=adjudicable&(z['pred_q']<.05)
   budget=max(int(real.sum()),int(called.sum()));good=called&(np.sign(z['pred_lfc'])==np.sign(z['truth_lfc']))
   fid=float(good.sum()/budget) if budget else None
   union=real|(valid&(z['pred_q']<.05));jac=float((real&(z['pred_q']<.05)).sum()/union.sum()) if union.any() else 1.
   eligible=real&np.isfinite(z['truth_lfc']);denom=abs(z['truth_lfc'][eligible]).sum()
   nmae=float(abs(np.nan_to_num(z['pred_lfc'][eligible],nan=0.,posinf=0.,neginf=0.)-z['truth_lfc'][eligible]).sum()/denom) if eligible.sum()>=10 and denom>0 else None
   for k,v in [('fidelity',fid),('jaccard',jac),('nmae',nmae)]:
    assert (v is None and r[k] is None) or (v is not None and r[k] is not None and abs(v-r[k])<2e-12),(r['family'],k)
    metrics_checked+=1
   profiles[(r['family'],r['seed_offset'],r['arm'],r['case_index'])]={k:z[k].copy() for k in ['pred','jack','truth','truth_jack','control','control_jack']}
 # Verify all original corrected null outcomes are identical for the unchanged emitter.
 nulls=run.read(ROOT/'experiments/v37_component_ensemble/analysis/null_results.json')['records'];nullchecks=[]
 for family in sorted({r['family'] for r in nulls}):
  pair=[next(r for r in nulls if r['family']==family and r['arm']==a) for a in run.ARMS]
  paths=[ROOT/'experiments/v37_component_ensemble/analysis'/r['arrays_path'] for r in pair]
  for r,p in zip(pair,paths):assert run.sha(p)==r['arrays_sha256']
  with np.load(paths[0],allow_pickle=False) as a,np.load(paths[1],allow_pickle=False) as b:
   for k in a.files:np.testing.assert_array_equal(a[k],b[k])
  assert pair[0]['null_DE_calls']==pair[1]['null_DE_calls'];nullchecks.append({'family':family,'calls':pair[0]['null_DE_calls'],'arrays_exact':True})
 def panels(exclude=None):
  results=[]
  for family in sorted({r['family'] for r in d['records']}):
   for offset in run.OFFSETS:
    for arm in run.ARMS:
     rr=[r for r in d['records'] if r['family']==family and r['seed_offset']==offset and r['arm']==arm and (family,r['target'])!=exclude]
     vals=[profiles[(family,offset,arm,r['case_index'])] for r in rr]
     p,t,pj,tj=(np.array([v[k] for v in vals]) for k in ('pred','truth','jack','truth_jack'))
     own=np.array([np.flatnonzero(genes==r['target'])[0] for r in rr])
     score=expression(p,t,vals[0]['control'],pj,tj,vals[0]['control_jack'],own)
     results.append({'family':family,'study':run.key(family),'seed_offset':offset,'arm':arm,**score})
  return results
 original=panels()
 for row in original:
  old=next(p for p in d['panels'] if all(p[k]==row[k] for k in ['family','seed_offset','arm']))
  for k in ['common_expression_numerator','common_expression_denominator','common_expression_ratio']:
   if row[k] is None:assert old[k] is None
   else:
    err=abs(row[k]-old[k]);max_error=max(max_error,err);assert err<max(2e-11,abs(row[k])*1e-11),(row['family'],k,err)
 influence=[]
 for family,target in sorted({(r['family'],r['target']) for r in d['records']}):
  rows=panels((family,target));a=run.aggregate(rows,'mean50','common_expression_ratio');b=run.aggregate(rows,'v24','common_expression_ratio')
  influence.append({'family':family,'target':target,'expression_ratio_relative_change':a/b-1})
 assert len(influence)==30
 record={'passed':True,'case_metrics_independently_replayed':metrics_checked,'expression_panels_independently_replayed':len(original),'maximum_expression_absolute_error':max_error,'null_array_identity':nullchecks,'expression_leave_context_target_out':influence,'expression_improves_every_leave_target_out':all(r['expression_ratio_relative_change']<0 for r in influence),'expression_ratio_warning':'HCT116 tiny denominator dominates raw ratios; no private anchor scaling or reliable future forecast claimed.','evaluation_sha256':run.sha(A/'evaluation.json'),'auditor_sha256':run.sha(Path(__file__))}
 run.write(A/'independent_evaluation_audit.json',record)
 print('Independent V55 audit PASS:',metrics_checked,'DE metrics,20 expression panels,5 unchanged nulls,30 target-removal evaluations.',flush=True)
 print('Expression improves in every leave-target-out:',record['expression_improves_every_leave_target_out'],flush=True)
if __name__=='__main__':main()
