"""Independently check V25 family exclusions, covariance math, identities and count replays."""
from pathlib import Path
import gc,json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from threadpoolctl import threadpool_limits
import build_submission as b
from v25_methods import ROOT,A,FAMILIES,NAMES,CONTEXTS,load_plan,load_model,raw,predict,emit,digest,write
from fit_v7_response import read_cache
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model
from assess_v25_experiment import compare

def main():
 from v24_runtime import require_science
 require_science();plan=load_plan();old=json.loads((ROOT/'experiments/v24_cd4_response/analysis/plan.json').read_text());assert b.sha256(ROOT/'assess_v24_experiment.py')==old['dependencies']['assess_v24_experiment.py']
 data=read_cache();genes=data['genes'];targets=data['targets'];models={f:load_model(f) for f in FAMILIES+['all']};covariances=0
 for f,model in models.items():
  rec=json.loads((A/'models'/f'{f}.json').read_text());contexts=list(model['source_contexts']);expected=[name for name,_ in CONTEXTS if (('CD4' if name.startswith('CD4') else name)!=f)] if f!='all' else [name for name,_ in CONTEXTS]
  assert contexts==expected==rec['training_contexts'];assert rec['excluded_families']==([] if f=='all' else [f]);pos=model['positions'];assert len(pos)==len(set(pos))==4096 and not np.isin(genes[pos],targets).any();assert np.isfinite(model['scale']).all() and model['scale'].min()>0
  np.testing.assert_allclose(model['basis'].T@model['basis'],np.eye(64),atol=1e-10);np.testing.assert_array_equal(model['genes'],genes);np.testing.assert_array_equal(model['targets'],targets)
  for key in ['source_shifts','source_weights','source_covariances']:assert np.isfinite(model[key]).all()
  assert (model['source_weights']>=0).all() and (model['source_weights'].sum(0)>0).any()
 for name,_ in CONTEXTS:
  x,m,measured,_=raw(name,genes);ids=np.flatnonzero(m.gene_target.to_numpy(str)=='Non-Targeting');depth=m.total_counts.to_numpy(float)[ids];assert (depth>0).all() and (np.asarray(x[ids].sum(1)).ravel()<=depth).all();batch=m['sample'].to_numpy(str)[ids];ld=np.log1p(depth.astype(float));n=len(ids)
  for f,model in models.items():
   contexts=list(model['source_contexts'])
   if name not in contexts:continue
   si=contexts.index(name);rec=json.loads((A/'models'/f'{f}.json').read_text());np.testing.assert_array_equal(rec['source_records'][si]['control_rows'],ids);pos=model['positions'];assert measured[pos].all()
   # Independently recompute dense log-depth and within-batch nuisance removal.
   y=x[ids][:,pos].toarray().astype(float)*(1e4/depth)[:,None];np.log1p(y,out=y);zd=ld.copy()
   for lane in set(batch):
    keep=batch==lane;y[keep]-=y[keep].mean(0);zd[keep]-=zd[keep].mean()
   beta=(y.T@zd)/max(float(zd@zd),n*1e-12);y-=zd[:,None]*beta;scores=(y/model['scale'])@model['basis'];scores-=scores.mean(0);emp=scores.T@scores/n;d=np.maximum(np.diag(emp),1e-8);expected=.75*emp+.25*np.diag(d)+.1*d.mean()*np.eye(64)
   np.testing.assert_allclose(model['source_covariances'][si],expected,rtol=1e-9,atol=1e-10);assert np.linalg.eigvalsh(expected).min()>0;covariances+=1;del y,scores
  del x;gc.collect();print('V25 independent covariance audit',name,flush=True)
 selected=json.loads((A/'selection.json').read_text());atlas,_=load_atlas_model();cases=caches=replays=nulls=0;stages=[('development',plan['contexts'],selected)]
 if (A/'confirmation_result.json').exists():
  opened=json.loads((A/'reserve_execution_started.json').read_text());assert selected['new_candidate_passes'] and opened['selection_sha256']==b.sha256(A/'selection.json') and opened['reserve_plan_sha256']==b.sha256(A/'reserve_plan.json');assert json.loads((A/'cd4_raw_independent_audit.json').read_text())['passed'];stages.append(('confirmation',[json.loads((A/'reserve_plan.json').read_text())],json.loads((A/'confirmation_result.json').read_text())))
 for stage,contexts,decision in stages:
  table=pd.read_csv(A/f'{stage}_scores.csv');panels=pd.read_csv(A/f'{stage}_panels.csv');evaluation=json.loads((A/f'{stage}_evaluation.json').read_text());assert evaluation['complete'] and evaluation['scores_sha256']==b.sha256(A/f'{stage}_scores.csv') and evaluation['panels_sha256']==b.sha256(A/f'{stage}_panels.csv');assert len(table)==len(NAMES)*sum(len(c['cells']) for c in contexts)
  for saved in decision['comparisons']:assert saved==json.loads(json.dumps(compare(panels,saved['candidate'],plan['raw_guards']),default=lambda v:bool(v)))
  for ctx in contexts:
   family=ctx['family'];key='CD4' if family.startswith('CD4') else family;assert ctx['excluded']==[key];x,m,measured,_=raw(family,genes);labels=m.gene_target.to_numpy(str);assert not set(ctx['left'])&set(ctx['right']);assert (labels[ctx['left']]=='Non-Targeting').all() and (labels[ctx['right']]=='Non-Targeting').all();p,bulk=moments(x[ctx['left']])[1:3];effects,_=predict(data,models[key],x,m,ctx['left'],measured,ctx['excluded']);own=np.array([np.flatnonzero(genes==t)[0] for t in targets]);assert all((v[np.arange(300),own]==0).all() and np.isfinite(v).all() for v in effects.values())
   for name in ['force25','bures25','diagonal25','shuffled25','source_program25']:assert abs(effects[name]-effects['v24']).max()<=.250001
   for name in ['force50','shuffled50']:assert abs(effects[name]-effects['v24']).max()<=.500001
   for c in ctx['cells']:
    n=200 if family=='CD4_Rest_D1' else 400;assert len(set(c['actual']))==len(set(c['template']))==n;assert (labels[c['actual']]==c['target']).all() and set(c['template'])<=set(ctx['left']);assert not set(c['actual'])&(set(ctx['left'])|set(ctx['right']));side=A/'parallel'/f'{family}_{c["seed"]}_{stage}'/'count_cache'/f'{family}_{c["target"]}_{c["seed"]}.json';rec=json.loads(side.read_text());assert rec['sha256']==b.sha256(side.with_suffix('.npz')) and rec['inputs']['plan']==b.sha256(A/'plan.json') and rec['inputs']['geometry_model']==b.sha256(A/'models'/f'{key}.npz');assert set(r['candidate'] for r in rec['records'])==set(NAMES)
    for row in rec['records']:
     found=table[(table.family==family)&(table.target==c['target'])&(table.seed==c['seed'])&(table.candidate==row['candidate'])];assert len(found)==1 and found.iloc[0].count_sha256==row['count_sha256']
    cases+=1;caches+=1
   c=ctx['cells'][0];ti=int(np.flatnonzero(targets==c['target'])[0]);template=x[c['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64)
   for name in ['v24','force25','bures25']:
    count,_=emit(template,effects[name][ti],int(own[ti]),c['seed']+ti,(libs,p,bulk),atlas);row=table[(table.family==family)&(table.target==c['target'])&(table.seed==c['seed'])&(table.candidate==name)].iloc[0];assert digest(count)==row.count_sha256;np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs);replays+=1
   from v16_methods import emit as champion_emit
   zero=np.zeros(len(genes));null,_=emit(template,zero,None,c['seed'],(libs,p,bulk),atlas);expected,_=champion_emit(template,np.arange(len(libs)),zero,None,c['seed'],'residual10',(libs,p,bulk),atlas);assert (null!=expected).nnz==0;nulls+=1;del x;gc.collect();print('V25 independent identity/count replay audit',stage,family,flush=True)
 write(A/'independent_experiment_audit.json',{'passed':True,'model_count':len(models),'independently_recomputed_source_covariances':covariances,'frozen_cases':cases,'count_populations':cases*len(NAMES),'cache_hashes_checked':caches,'full_axis_count_replays':replays,'no_effect_emitter_V16_identity_cases':nulls,'whole_held_family_excluded_from_means_and_geometry':True,'transitive_assessor_dependency_verified_against_V24_plan':True,'plan_sha256':b.sha256(A/'plan.json'),'selection_sha256':b.sha256(A/'selection.json'),'auditor_sha256':b.sha256(Path(__file__)),'limitations':'Covariance recomputation verifies implementation, not causal identification. Source programs use linearized means. Null emitter identity does not establish null-DE calibration. Development reuses public biology, CD4 D1 has200 actual cells, HCT116 NMAE missing.'});print('V25 independent audit passed',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
