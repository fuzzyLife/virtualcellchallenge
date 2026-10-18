"""Fixed V61 joint50 production candidate; explicit experimental status on failed guards."""
import argparse,gc
import numpy as np
from threadpoolctl import threadpool_limits
import v61_pooled_state_response as v
m=v.m;A=v.A;OUT=A.parent

def kernel(base,q,donors):
 genes,targets=m.axes();axis=~np.isin(genes,targets);rhos={s:m.spearman(q['profile'][axis],d['profile'][axis]) for s,d in donors.items()};result=base.copy();support_counts=[]
 for ti,target in enumerate(targets):
  own=list(genes).index(target);num=np.zeros((4,len(genes)));den=np.zeros_like(num);ln=np.zeros(4);ld=0.
  for s,d in donors.items():
   if d['n'][ti]<=0:continue
   count=sum(m.GROUP[t]==m.GROUP[s] for t in donors);weight=float(d['n'][ti]/(d['n'][ti]+64)*np.exp(2*rhos[s])/count);mask=d['valid'][ti];num+=weight*d['within'][ti]*mask;den+=weight*mask;ln+=weight*d['tilt'][ti];ld+=weight
  within=np.divide(num,den,out=np.zeros_like(num),where=den>0);tilt=ln/max(ld,1e-12);mask=den.any(0)&q['measured'];mask[own]=False
  new=m.logratio(m.reconstruct(q['control_state_counts'],q['control_state_means'],tilt,within),q['profile']);result[ti,mask]=np.clip(base[ti,mask]+.5*np.clip(new[mask]-base[ti,mask],-1,1),-3,3);result[ti,own]=0;support_counts.append(int(mask.sum()))
 return result,dict(source_rank_correlations=rhos,supported_coordinates=support_counts)

def check_kernel():
 v.activate();checks=[]
 for source in m.SOURCES:
  outer=m.GROUP[source];q=m.mechanism(outer,source);q={k:q[k] for k in ['profile','measured','control_state_counts','control_state_means']};donors={s:m.mechanism(outer,s) for s in m.SOURCES if m.GROUP[s]!=outer};p=m.CACHE/'predictions'/(source+'.npz');rec=m.read(m.A/'predictions'/(source+'.json'));assert m.sha(p)==rec['sha256'];saved=m.load(p);actual,_=kernel(saved['v24'],q,donors);np.testing.assert_array_equal(actual,saved['joint50']);checks.append(dict(source=source,bitwise_equal=True,coordinates=int(actual.size)));print('V61 independent production kernel parity',source,flush=True)
 m.write(A/'production_kernel_parity.json',dict(passed=True,checks=checks,kernel_sha256=m.sha(__file__),plan_sha256=m.sha(A/'plan.json'),D3_opened=False))

def freeze():
 assert not (A/'production_plan.json').exists()
 scripts=['fit_v61_submission.py','build_v61_submission.py','package_v61_submission.py','validate_v61_submission.py','run_v61_release.py','verify_v61_effects.py']
 m.write(A/'production_plan.json',dict(candidate='joint50',version=61,scientific_plan_sha256=m.sha(A/'plan.json'),candidate_selection='Fixed primary V61 joint50, no post-hoc diagnostic promotion or blending grid. V59/V60 are reported comparisons.',experimental_policy='The user requested a new release. If public gates fail, --experimental may prepare this fixed candidate for a clearly labeled diagnostic upload, preserving V24 and reporting every failed gate. It must not be described as a validated improvement.',dependencies={s:m.sha(m.ROOT/s) for s in scripts},contexts=['A','B','C'],targets=300,cells_per_target_context=400,genes=18533,template_seed=20261011,emitter='residual10',D3_opened=False,automatically_upload=False))
 print('V61 fixed production candidate frozen',m.sha(A/'production_plan.json'),flush=True)

def fit(experimental=False):
 v.activate();p=m.read(A/'production_plan.json')
 for s,h in p['dependencies'].items():assert m.sha(m.ROOT/s)==h,s
 verdict=m.read(A/'count_assessment.json');assert verdict['complete'];assert experimental or verdict['release_supported'],'Failed public gates require explicit experimental artifact mode'
 assert m.read(A/'production_kernel_parity.json')['passed']
 from build_v24_submission import load as champion
 from v6_methods import moments
 import build_submission as b
 manifest,genes,targets,oldrecord,base=champion();np.testing.assert_array_equal(genes,m.axes()[0]);np.testing.assert_array_equal(targets,m.axes()[1])
 donors={s:m.mechanism('all',s) for s in m.SOURCES};model=m.state_model('all');effects=[];routing={};controls=[]
 assert not (A/'reference_model.json').exists(),'Preserve V61 fitted model'
 for ci,context in enumerate(manifest['contexts']):
  x,guides,library,*_=b.load_control(m.ROOT,context,genes,manifest);xx=x.astype(float).multiply((1e6/library)[:,None]).tocsr();meta=dict(rows=np.arange(x.shape[0]),label=np.repeat('Non-Targeting',x.shape[0]),lane=np.repeat(context,x.shape[0]),measured=np.ones(len(genes),bool),depth=np.asarray(library))
  state,trace=m.assign(xx,meta,model);nc,cm,profile=m.group_summary(xx,state,np.arange(xx.shape[0]));q=dict(control_state_counts=nc,control_state_means=cm,profile=profile,measured=np.ones(len(genes),bool));effect,r=kernel(base[ci],q,donors);effects.append(effect);routing[context]=r
  path=A/'production_inputs'/(context+'_query.npz');digest=m.save(path,**q,baseline=base[ci],query_mean=moments(x)[1]*1e6,genes=np.asarray(genes,str),state=state,library=np.asarray(library),guides=np.asarray(guides,str));controls.append(dict(context=context,file=f'context_{context}.h5ad',sha256=m.sha(m.ROOT/f'context_{context}.h5ad'),query_sha256=digest,control_cells=int(x.shape[0]),state_counts=nc.tolist(),transform=trace))
  print('V61 query fitted',context,'changed coordinates',int(np.count_nonzero(effect!=base[ci])),flush=True);del x,xx;gc.collect()
 effects=np.asarray(effects,np.float32);path=A/'reference_model_arrays.npz';digest=m.save(path,genes=np.asarray(genes,str),targets=np.asarray(targets,str),effects=effects,parent_effects=base)
 records=[]
 for ci,c in enumerate(manifest['contexts']):
  for ti,t in enumerate(targets):records.append(dict(context=c,target_gene=str(t),sum_absolute_LFC_correction=float(abs(effects[ci,ti]-base[ci,ti]).sum()),changed_genes=int(np.count_nonzero(effects[ci,ti]!=base[ci,ti]))))
 import pandas as pd
 pd.DataFrame(records).to_csv(A/'production_response_changes.csv',index=False)
 m.write(A/'reference_model.json',dict(model='V61 pooled state-abundance and conditional response transfer',candidate='joint50',status='public_supported_experimental_candidate' if verdict['release_supported'] else 'experimental_failed_public_guards',available_public_guards_passed=verdict['release_supported'],confirmed_improvement=False,emitter='residual10',model_arrays_sha256=digest,parent_V24_arrays_sha256=m.sha(m.ROOT/'experiments/v24_cd4_response/analysis/reference_model_arrays.npz'),parent_V24_model_sha256=m.sha(m.ROOT/'experiments/v24_cd4_response/analysis/reference_model.json'),scientific_plan_sha256=m.sha(A/'plan.json'),production_plan_sha256=m.sha(A/'production_plan.json'),count_assessment_sha256=m.sha(A/'count_assessment.json'),kernel_parity_sha256=m.sha(A/'production_kernel_parity.json'),routing=routing,context_original_inputs=controls,training_sources=list(donors),dependencies=p['dependencies'],D3_opened=False,uploaded=False))
 print('V61 fixed production model saved; public support:',verdict['release_supported'],flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','check','fit']);p.add_argument('--experimental',action='store_true');a=p.parse_args()
 with threadpool_limits(limits=2):
  if a.stage=='freeze':freeze()
  elif a.stage=='check':check_kernel()
  else:fit(a.experimental)
