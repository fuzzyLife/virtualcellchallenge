"""Post-development diagnostic: does a scalar magnitude control explain mean gains?"""
from pathlib import Path
import numpy as np
from v57_pie import read,write,sha
import build_submission as b
from audit_transfer_validation import pds
from v63_identity_guard import probabilities
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v63_identity_guard/analysis'

def main():
 assert (A/'additional_assessment.json').exists(),'Finish all fixed raw comparisons first'
 manifest,genes,targets=b.read_spec(ROOT);ownall=np.array([list(genes).index(t) for t in targets]);panel=np.isin(genes,targets);rec=[]
 plan=A/'amplitude_unique_control_plan.json'
 if not plan.exists():write(plan,dict(script_sha256=sha(__file__),timing='Post-development diagnostic, artifact estimator already fixed and fitted; not eligible for deployment or reselection.',definition='Multiply V24 LFC by the ratio of guarded/V24 nonpanel LFC vector norms per target, clipped0..2; preserve signs and zero coordinates. This matches LFC magnitude, not exactly the nonlinear pooled-profile norm. Same own-target factor. Only desired means, no new generated cells.'))
 assert read(plan)['script_sha256']==sha(__file__)
 for source in ('H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2'):
  with np.load(A/'predictions'/(source+'.npz'),allow_pickle=False) as z:base=z['all_baseline'].copy();new=z['all_effects'].copy()
  n0=np.linalg.norm(base[:,~panel],axis=1);n1=np.linalg.norm(new[:,~panel],axis=1);scale=np.clip(np.divide(n1,n0,out=np.ones(len(n0)),where=n0>0),0,2);amp=np.clip(base*scale[:,None],-3,3);amp[np.arange(300),ownall]=0
  with np.load(ROOT/'experiments/v37_component_ensemble/analysis/cache'/source/'effects.npz',allow_pickle=False) as z:bulk=z['bulk'].copy()
  rows=sorted(list((A/'development').glob(source+'__*__0__v24.json'))+list((A/'additional').glob(source+'__*__0__v24.json')))
  canonical={}
  for f in sorted(rows,key=lambda f:read(f)['case_index']):canonical.setdefault(read(f)['target'],f)
  rows=list(canonical.values())
  truth=[];own=[];names=[];indices=[];ctrl=None
  for f in rows:
   r=read(f);file=f.with_suffix('.npz');assert sha(file)==r['arrays_sha256']
   with np.load(file,allow_pickle=False) as z:truth.append(z['truth'].copy());ctrl=z['control'].copy()
   indices.append(list(targets).index(r['target']));own.append(list(genes).index(r['target']));names.append(r['target'])
  assert len(rows) in (4,6,8),(source,len(rows));truth=np.asarray(truth);ix=np.array(indices);own=np.array(own);mask=np.ones(truth.shape,bool);mask[np.arange(len(own)),own]=False;metrics={}
  for name,e in [('v24',base),('guarded_assay',new),('amplitude_only',amp)]:
   prediction=np.log1p(50000*probabilities(e,bulk,ownall)[ix]);metrics[name]=dict(profile_MSE=float(np.sum((prediction-truth)**2*mask)/mask.sum()),pds=pds(prediction[:,~panel],truth[:,~panel],ctrl[~panel]))
  rec.append(dict(family=source,study='CD4' if source.startswith('CD4') else 'Orion' if source in ('HCT116','HEK293T') else 'H1',targets=names,metrics=metrics,amplitude_scales=scale[ix].tolist()));print('Amplitude control',source,metrics,flush=True)
 def agg(arm,metric):
  return float(np.mean([np.mean([r['metrics'][arm][metric] for r in rec if r['study']==s]) for s in ('H1','Orion','CD4')]))
 summary={arm:{k:agg(arm,k) for k in ('profile_MSE','pds')} for arm in ('v24','guarded_assay','amplitude_only')}
 write(A/'amplitude_unique_control_results.json',dict(records=rec,summary=summary,script_sha256=sha(__file__),scope='Desired-mean control on 30 unique target/context combinations from previously reused measured cases; nonlinear assay norms not exactly matched. No candidate promotion or fitted-model change.'))
 print('Amplitude diagnostic summary',summary,flush=True)

if __name__=='__main__':main()
