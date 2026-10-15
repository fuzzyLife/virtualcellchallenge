"""Independently audit V26 source identities, fit-only moments and gated count replays."""
from pathlib import Path
import argparse,json,gc
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
from v26_methods import ROOT,A,FAMILIES,KINDS,NAMES,load_plan,family_key,raw,load_decoder,geometry,features,state_values,emit,digest,write,moments
from evaluate_atlas_transfer import load_atlas_model
from v12_cpu_metrics import de
import build_submission as b

def main(inputs_only=False):
 from v24_runtime import require_science
 require_science();plan=load_plan();sources=rows_checked=models_checked=cases=replays=0;trained=[];means_checked=0
 for key in FAMILIES:
  dataset=A/'training'/f'{key}.npz'
  if not dataset.exists():continue
  trained.append(key);rec=json.loads(dataset.with_suffix('.json').read_text());assert rec['sha256']==b.sha256(dataset) and rec['plan_sha256']==b.sha256(A/'plan.json');matrix=sparse.load_npz(dataset.with_name(f'{key}_counts.npz'));assert rec['counts_sha256']==b.sha256(dataset.with_name(f'{key}_counts.npz'))
  with np.load(dataset,allow_pickle=False) as z:d={n:z[n].copy() for n in z.files}
  assert matrix.shape==(rec['rows'],18533) and (matrix.data>=0).all() and np.array_equal(matrix.data,np.rint(matrix.data));assert key not in [family_key(f) for f in rec['training_contexts']];np.testing.assert_allclose(d['state_basis'].T@d['state_basis'],np.eye(8),atol=1e-10);assert (d['state_scale']>0).all();assert not np.isin(d['genes'][d['state_axis']],d['targets']).any();assert not d['output_mask'][np.isin(d['genes'],d['targets'])].any();offset=0;group_records=rec['group_records']
  for original in rec['original_inputs']:
   family=original['context'];assert family_key(family)!=key;x,meta,measured,_=raw(family,d['genes']);ids=np.array(original['selected_original_rows']);fit=np.array(original['fit_flags'],bool);labels=meta.gene_target.to_numpy(str)[ids];batches=meta['sample'].to_numpy(str)[ids];assert len(set(ids))==len(ids);np.testing.assert_array_equal(batches,original['batch_labels']);np.testing.assert_array_equal(d['fit'][offset:offset+len(ids)],fit);libs=np.asarray(matrix[offset:offset+len(ids)].sum(1)).ravel();np.testing.assert_array_equal(libs,original['projected_libraries']);assert np.all(libs<=np.array(original['original_libraries']));take=np.random.default_rng(26262).choice(len(ids),min(32,len(ids)),replace=False);assert (matrix[offset+take]!=x[ids[take]]).nnz==0;rows_checked+=len(take)
   ctrl={}
   for lane in sorted(set(batches)):
    fitctrl=np.flatnonzero((batches==lane)&(labels=='Non-Targeting')&fit);assert len(fitctrl)>=8;block=matrix[offset+fitctrl];depth=np.asarray(block.sum(1)).ravel();ctrl[lane]=np.asarray(block.multiply((1/depth)[:,None]).mean(0)).ravel()
   for group in np.unique(d['groups'][offset:offset+len(ids)]):
    local=np.flatnonzero(d['groups'][offset:offset+len(ids)]==group);assert (labels[local]==group_records[group]['label']).all();chosen=local[fit[local]];assert len(chosen)==group_records[group]['training_cells'];block=matrix[offset+chosen];depth=np.asarray(block.sum(1)).ravel();p=np.asarray(block.multiply((1/depth)[:,None]).mean(0)).ravel();c=np.mean([ctrl[lane] for lane in batches[chosen]],axis=0);smooth=p+.001*c+1e-8*measured;smooth/=smooth.sum();np.testing.assert_allclose(d['base'][group],smooth,rtol=1e-6,atol=1e-10);np.testing.assert_allclose(d['control'][group],c,rtol=1e-6,atol=1e-10);means_checked+=1
   shuffle=d['shuffled'][offset:offset+len(ids)]-offset;assert shuffle.min()>=0 and shuffle.max()<len(ids) and len(set(shuffle))==len(ids);np.testing.assert_array_equal(fit,fit[shuffle]);np.testing.assert_array_equal(batches,batches[shuffle]);q=features(matrix[offset+take],d['base'][d['groups'][offset+take]],geometry(key)['positions']);np.testing.assert_array_equal(q,d['features'][offset+take]);mini=x[ids[take]];depth=np.asarray(mini.sum(1)).ravel();rr=np.log1p(1e4*mini[:,d['state_axis']].multiply((1/depth)[:,None]).toarray())-np.array([np.log1p(1e4*ctrl[lane][d['state_axis']]) for lane in batches[take]]);expected_state=np.clip(rr.astype(np.float32)/d['state_scale']@d['state_basis'],-4,4).astype(np.float32);np.testing.assert_allclose(expected_state,d['state'][offset+take],atol=2e-5,rtol=2e-5);sources+=1;offset+=len(ids);del x;gc.collect();print('V26 independent original-count/group audit',key,family,flush=True)
  assert offset==matrix.shape[0]
  for kind in KINDS:
   if inputs_only or not (A/'models'/f'{key}_{kind}.json').exists():continue
   model=load_decoder(key,kind);assert model['record']['excluded_family']==key and model['record']['training_contexts']==rec['training_contexts'];assert model['record']['training_record_sha256']==b.sha256(dataset.with_suffix('.json'));assert len(model['record']['history'])==12;assert all(np.isfinite(h['training_elbo']) and np.isfinite(h['validation_elbo']) for h in model['record']['history']);assert all(np.isfinite(w).all() for w in model['weights'].values());models_checked+=1
  del matrix,d;gc.collect()
 atlas=None if inputs_only else load_atlas_model()[0]
 for key in ([] if inputs_only else trained):
  result=A/f'null_{key}_result.json'
  if not result.exists():continue
  evaluation=json.loads((A/f'null_{key}_evaluation.json').read_text());table=pd.read_csv(A/f'null_{key}_scores.csv');assert evaluation['scores_sha256']==b.sha256(A/f'null_{key}_scores.csv');assert evaluation['panels_sha256']==b.sha256(A/f'null_{key}_panels.csv');decision=json.loads(result.read_text());expected={};conditional=table[table.candidate=='conditional'];expected['mean_contract']=bool((conditional.mean_reference_L1<=plan['mean_contract']['maximum_mean_CPM_probability_L1']).all());policy=plan['null_policy']
  for _,r in conditional.iterrows():
   expected[str(r.family)+'_DE_calls']=bool(r.null_DE_calls<=max(policy['false_call_floor'],r.real_null_DE_calls+policy['real_null_call_allowance']));expected[str(r.family)+'_detection']=bool(r.detection_MAE<=max(policy['detection_MAE_floor'],policy['detection_real_multiplier']*r.real_detection_MAE+policy['detection_real_allowance']));expected[str(r.family)+'_variance']=bool(r.log_variance_MAE<=max(policy['log_variance_MAE_floor'],policy['log_variance_real_multiplier']*r.real_log_variance_MAE+policy['log_variance_real_allowance']))
  assert decision['checks']==expected and bool(decision['passed'])==all(expected.values())
  for ctx in [c for c in plan['contexts'] if family_key(c['family'])==key]:
   case=ctx['null'];x,meta,measured,_=raw(ctx['family'],geometry(key)['genes']);labels=meta.gene_target.to_numpy(str);actual=set(case['actual']);reference=set(case['reference']);fit=set(case['fit_controls']);assert len(actual)==len(reference)==len(set(case['template']))==400;assert not actual&reference and not actual&fit and not reference&fit;assert set(case['template'])<=fit;assert (labels[list(actual|reference|fit)]=='Non-Targeting').all();side=A/'parallel'/f'null_{key}'/ctx['family']/f'Non-Targeting_{case["seed"]}.json';cache=json.loads(side.read_text());assert cache['sha256']==b.sha256(side.with_suffix('.npz'));assert set(r['candidate'] for r in cache['records'])==set(NAMES)
   for r in cache['records']:
    row=table[(table.family==ctx['family'])&(table.candidate==r['candidate'])];assert len(row)==1 and row.iloc[0].count_sha256==r['count_sha256']
   model=load_decoder(key,'conditional');ids=case['fit_controls'];p,bulk=moments(x[ids])[1:3];state=state_values(x,meta,ids,model['state_axis'],model['state_basis'],model['state_scale']);lookup={int(v):i for i,v in enumerate(ids)};state=state[[lookup[int(v)] for v in case['template']]];template=x[case['template']];libs=np.asarray(template.sum(1)).ravel().astype(np.int64);zero=np.zeros(x.shape[1]);base,_=emit(template,zero,None,case['seed'],'v24',(libs,p,bulk),atlas,model,state)
   for name,count in [('v24',base),('conditional',emit(template,zero,None,case['seed'],'conditional',(libs,p,bulk),atlas,model,state,reference=base)[0])]:
    row=table[(table.family==ctx['family'])&(table.candidate==name)].iloc[0];assert digest(count)==row.count_sha256;np.testing.assert_array_equal(np.asarray(count.sum(1)).ravel(),libs);np.testing.assert_array_equal(np.asarray(count.sum(0)).ravel(),np.asarray(base.sum(0)).ravel());drift=float(abs(np.asarray(count.multiply((1/libs)[:,None]).mean(0)).ravel()-np.asarray(base.multiply((1/libs)[:,None]).mean(0)).ravel()).sum());assert abs(drift-row.mean_reference_L1)<1e-10;stats=de(count,x[case['reference']],geometry(key)['genes']);calls=int(((stats['q']<.05)&stats['eligible']).sum());assert calls==int(row.null_DE_calls);replays+=1
   cases+=1;del x;gc.collect();print('V26 independent null replay audit',ctx['family'],flush=True)
 if (A/'calibration_stop.json').exists():
  stop=json.loads((A/'calibration_stop.json').read_text());assert not json.loads((A/f'null_{stop["failed_family"]}_result.json').read_text())['passed'];assert stop['calibration_sha256']==b.sha256(A/f'null_{stop["failed_family"]}_result.json');assert not (A/'development_evaluation.json').exists() and not (A/'reserve_execution_started.json').exists()
 write(A/('independent_training_input_audit.json' if inputs_only else 'independent_audit.json'),{'passed':True,'training_folds':trained,'source_count_contexts_checked':sources,'original_full_axis_rows_checked':rows_checked,'training_only_group_means_checked':means_checked,'decoder_models_checked':models_checked,'null_cases_checked':cases,'full_axis_count_replays':replays,'frozen_family_exclusions_verified':True,'frozen_null_decisions_recomputed':True,'plan_sha256':b.sha256(A/'plan.json'),'auditor_sha256':b.sha256(Path(__file__)),'scope':'This audits the training inputs, saved fits and executed null gates. It does not substitute for a treatment/generalization test. Conditional states remain observational proxies.'});print('V26 independent audit passed',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--inputs-only',action='store_true');args=p.parse_args()
 with threadpool_limits(limits=2):main(args.inputs_only)
