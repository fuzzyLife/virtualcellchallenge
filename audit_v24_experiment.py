"""Audit V24 frozen source/condition exclusions, all identities/caches and exact replays."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from threadpoolctl import threadpool_limits
from v24_methods import ROOT,A,load_plan,make_teacher,raw,predict,emit,training_data,write,digest,NAMES
from fit_v7_response import read_cache
from v6_methods import moments
from evaluate_atlas_transfer import load_atlas_model
from assess_v24_experiment import compare
from build_cd4_model import load_subset
import build_submission as b

def main():
 from v24_runtime import require_science
 require_science()
 plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];teacher=make_teacher(False);tr=training_data(data);assert not ((tr['row_families']=='CD4')&(tr['row_contexts']==3)).any();record=json.loads((A/'models/teacher_development.json').read_text());assert record['training_conditions']==['Rest','Stim8hr'] and record['withheld_response_condition']=='Stim48hr';meta,source=load_subset(ROOT/'reference/cd4_probe','merged')
 for name in ['merged','CE0006864_CE0008162','CE0008678_CE0010866']:
  _,values=load_subset(ROOT/'reference/cd4_probe',name);q=values['adj_p_value'];finite=q[np.isfinite(q)];assert finite.min()>=0 and finite.max()<=1
 assert abs(teacher['value']).max()<=3.000001
 for r in record['source_rows']:assert str(meta.iloc[r['source_row']].culture_condition)==r['condition'] and r['condition']!='Stim48hr'
 assert not np.isin(genes[teacher['state_axis']],targets).any();assert np.isfinite(teacher['value']).all() and (teacher['confidence_de']>=0).all() and (teacher['confidence_de']<=1).all();assert (teacher['value'][:,np.arange(300),teacher['own']]==0).all()
 selected=json.loads((A/'selection.json').read_text());confirmation=json.loads((A/'confirmation_result.json').read_text());opened=json.loads((A/'reserve_execution_started.json').read_text());assert opened['selection_sha256']==b.sha256(A/'selection.json');assert len(confirmation['comparisons'])==1 and confirmation['comparisons'][0]['candidate']==selected['candidate']
 raw_audit=json.loads((A/'cd4_raw_independent_audit.json').read_text());assert raw_audit['passed'];reserve=json.loads((A/'reserve_plan.json').read_text());assert opened['reserve_plan_sha256']==b.sha256(A/'reserve_plan.json');atlas,_=load_atlas_model();cases=0;caches=0;replays=0;nulls=0;counts={}
 for stage,contexts in [('development',plan['contexts']),('confirmation',[reserve])]:
  table=pd.read_csv(A/f'{stage}_scores.csv');panels=pd.read_csv(A/f'{stage}_panels.csv');evaluation=json.loads((A/f'{stage}_evaluation.json').read_text());assert b.sha256(A/f'{stage}_scores.csv')==evaluation['scores_sha256'] and b.sha256(A/f'{stage}_panels.csv')==evaluation['panels_sha256'];assert len(table)==len(NAMES)*sum(len(c['cells']) for c in contexts)
  comparisons=selected['comparisons'] if stage=='development' else confirmation['comparisons']
  for saved in comparisons:
   expected=compare(panels,saved['candidate'],plan['raw_guards']);assert saved==json.loads(json.dumps(expected,default=lambda v:bool(v)))
  for ctx in contexts:
   x,m,measured,_=raw(ctx['family'],genes);labels=m.gene_target.to_numpy(str);assert not set(ctx['left'])&set(ctx['right']);assert (labels[ctx['left']]=='Non-Targeting').all() and (labels[ctx['right']]=='Non-Targeting').all();p0=moments(x[ctx['left']]);p,bulk=p0[1:3];query=p*1e6;query[~measured]=np.nan;effect,diag=predict(data,teacher,query,ctx['excluded']);assert abs(effect['cdf50_de']-effect['v16']).max()<=.250001;assert (effect['cdf50_de'][np.arange(300),teacher['own']]==0).all()
   for c in ctx['cells']:
    n=200 if ctx['family']=='CD4_Rest_D1' else 400;assert len(set(c['actual']))==len(set(c['template']))==n;assert (labels[c['actual']]==c['target']).all() and set(c['template'])<=set(ctx['left']);assert not set(c['actual'])&(set(ctx['left'])|set(ctx['right']));folder=A/'parallel'/f'{ctx["family"]}_{c["seed"]}_{stage}'/'count_cache';side=folder/f'{ctx["family"]}_{c["target"]}_{c["seed"]}.json';rec=json.loads(side.read_text());assert rec['sha256']==b.sha256(side.with_suffix('.npz')) and rec['inputs']['plan']==b.sha256(A/'plan.json');assert rec['inputs']['teacher']==b.sha256(A/'models/teacher_development.npz');assert set(r['candidate'] for r in rec['records'])==set(NAMES)
    for r in rec['records']:
     found=table[(table.family==ctx['family'])&(table.target==c['target'])&(table.seed==c['seed'])&(table.candidate==r['candidate'])];assert len(found)==1 and found.iloc[0].count_sha256==r['count_sha256']
    cases+=1;caches+=1
   c=ctx['cells'][0];ti=int(np.flatnonzero(targets==c['target'])[0]);template=x[c['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64)
   for name in ['v16','cdf50_de']:
    count,_=emit(template,effect[name][ti],int(teacher['own'][ti]),c['seed']+ti,(libs,p,bulk),atlas);r=table[(table.family==ctx['family'])&(table.target==c['target'])&(table.seed==c['seed'])&(table.candidate==name)].iloc[0];assert digest(count)==r.count_sha256;np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs)
   replays+=1
   # The unchanged emitter's no-target/no-effect contract is identity to itself,
   # not literal identity to real control counts or calibrated null p-values.
   from v16_methods import emit as champion_emit
   null,_=emit(template,np.zeros(len(genes)),None,c['seed'],(libs,p,bulk),atlas);expected,_=champion_emit(template,np.arange(len(libs)),np.zeros(len(genes)),None,c['seed'],'residual10',(libs,p,bulk),atlas);assert (null!=expected).nnz==0;nulls+=1;counts[ctx['family']]=len(ctx['cells']);del x;gc.collect();print('V24 independent case/replay audit',stage,ctx['family'],flush=True)
 write(A/'independent_experiment_audit.json',{'passed':True,'frozen_case_count':cases,'count_populations':cases*len(NAMES),'cache_hashes_checked':caches,'full_axis_replayed_cases':replays,'full_axis_two_arm_replays':replays*2,'no_effect_emitter_V16_identity_cases':nulls,'condition_response_exclusion_verified':True,'teacher_source_rows_and_conditions_verified':len(record['source_rows']),'independent_raw_input_audit_sha256':b.sha256(A/'cd4_raw_independent_audit.json'),'selected_candidate':selected['candidate'],'reserve_release_guards_pass':bool(confirmation['new_candidate_passes']),'plan_sha256':b.sha256(A/'plan.json'),'selection_sha256':b.sha256(A/'selection.json'),'confirmation_sha256':b.sha256(A/'confirmation_result.json'),'auditor_sha256':b.sha256(Path(__file__)),'limitations':'Rest development overlaps public summaries; condition reserve is existing study and upstream nuisance preprocessing may share information. Development HCT116 NMAE missing. Baseline refit excludes48hr; published production V16 arrays are separately retained. Null identity does not prove null-DE calibration.'});print('V24 independent experiment audit passed',cases,'cases',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
