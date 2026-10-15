"""Independent V27 source, paired-population and full-axis replay audit."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
from v27_methods import ROOT,A,ARMS,FINALIST,load_plan,load_model,emit,write,digest
from v26_methods import raw,moments,baseline,state_values,emit as old_emit,family_key
from fit_v7_response import read_cache
from evaluate_atlas_transfer import load_atlas_model
from v12_cpu_metrics import de
import build_submission as b

def main(inputs_only=False):
 plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];ns=nr=nm=nc=0;replays=[];fitted=[]
 for key in ['H1','HCT116','HEK293T','CD4']:
  ds=A/'training'/f'{key}.npz'
  if not ds.exists():continue
  fitted.append(key);rec=json.loads(ds.with_suffix('.json').read_text());matrix=sparse.load_npz(ds.with_name(f'{key}_counts.npz'))
  assert rec['sha256']==b.sha256(ds) and rec['counts_sha256']==b.sha256(ds.with_name(f'{key}_counts.npz'));assert rec['plan_sha256']==b.sha256(A/'plan.json')
  with np.load(ds,allow_pickle=False) as z:d={n:z[n].copy() for n in z.files}
  assert not np.isin(d['genes'][d['state_axis']],targets).any();assert not d['output_mask'][np.isin(genes,targets)].any();assert not d['measured_mask'][~d['output_mask']&~np.isin(genes,targets)].any()
  assert matrix.shape[1]==18533 and np.array_equal(matrix.data,np.rint(matrix.data)) and (matrix.data>=0).all();np.testing.assert_array_equal(d['genes'],genes)
  if key=='H1':
   # Check all shared source-training arrays against independently audited V26.
   previous=ROOT/'experiments/v26_count_decoder/analysis/training/H1.npz'
   with np.load(previous,allow_pickle=False) as old:
    for name in old.files:np.testing.assert_array_equal(d[name],old[name])
  offset=0;common=np.ones(len(genes),bool)
  for rr in rec['original_inputs']:
   assert family_key(rr['context'])!=key;x,m,measured,_=raw(rr['context'],genes);common&=measured;ids=np.array(rr['selected_original_rows']);fit=np.array(rr['fit_flags'],bool);labels=m.gene_target.to_numpy(str)[ids];batches=m['sample'].to_numpy(str)[ids];n=len(ids);assert len(set(ids))==n
   np.testing.assert_array_equal(d['fit'][offset:offset+n],fit);np.testing.assert_array_equal(batches,rr['batch_labels']);take=np.random.default_rng(27272).choice(n,min(32,n),replace=False);assert (matrix[offset+take]!=x[ids[take]]).nnz==0;nr+=len(take)
   ctrl={}
   for lane in sorted(set(batches)):
    j=np.flatnonzero((batches==lane)&(labels=='Non-Targeting')&fit);assert len(j)>=8;v=matrix[offset+j];ll=v.sum(1).A1;ctrl[lane]=np.asarray(v.multiply((1/ll)[:,None]).mean(0)).ravel()
   for group in np.unique(d['groups'][offset:offset+n]):
    j=np.flatnonzero((d['groups'][offset:offset+n]==group)&fit);v=matrix[offset+j];ll=v.sum(1).A1;p=np.asarray(v.multiply((1/ll)[:,None]).mean(0)).ravel();c=np.mean([ctrl[a] for a in batches[j]],axis=0);s=p+.001*c+1e-8*measured;s/=s.sum();np.testing.assert_allclose(s,d['base'][group],rtol=1e-6,atol=1e-10);nm+=1
   offset+=n;ns+=1;del x;gc.collect()
  np.testing.assert_array_equal(common,d['measured_mask']);assert offset==matrix.shape[0]
  panel=load_model(key,True);np.testing.assert_array_equal(panel['output_mask'],common);assert panel['record']['excluded_family']==key and panel['record']['training_contexts']==rec['training_contexts'];assert len(panel['record']['history'])==12;assert all(np.isfinite(w).all() for w in panel['weights'].values());del matrix,d;gc.collect()
 for key in ([] if inputs_only else fitted):
  for stage in ['null','treatment']:
   path=A/f'{stage}_{key}_evaluation.json'
   if not path.exists():continue
   ev=json.loads(path.read_text());assert ev['complete'];table=pd.read_csv(A/f'{stage}_{key}_scores.csv');assert ev['scores_sha256']==b.sha256(A/f'{stage}_{key}_scores.csv') and ev['panels_sha256']==b.sha256(A/f'{stage}_{key}_panels.csv');assert set(table.candidate)==set(ARMS)
   for ctx in [c for c in plan['contexts'] if family_key(c['family'])==key]:
    cells=ctx['null_cases'] if stage=='null' else ctx['cells'];x,m,measured,_=raw(ctx['family'],genes)
    for cell in cells:
     ids=cell['fit_controls'] if stage=='null' else ctx['left'];refids=cell['reference'] if stage=='null' else ctx['right'];fitids=cell['calibration_fit'] if stage=='null' else [i for i in ids if i not in set(cell['template'])];sets=[set(cell['actual']),set(refids),set(cell['template']),set(fitids)];assert all(not sets[i]&sets[j] for i in range(4) for j in range(i));assert len(sets[2])==400 and len(sets[3])>=200;assert (m.gene_target.to_numpy(str)[list(sets[1]|sets[2]|sets[3])]=='Non-Targeting').all()
     folder=A/'parallel'/f'{stage}_{key}'/ctx['family'];side=folder/f'{cell["target"]}_{cell["seed"]}.json';cache=json.loads(side.read_text());assert cache['sha256']==b.sha256(side.with_suffix('.npz'));assert cache['calibration_fit_rows']==fitids
     for arm in ['v24',FINALIST,'target_shuffled']:
      f=folder/f'{cell["target"]}_{cell["seed"]}_{arm}_counts.npz';count=sparse.load_npz(f);r=table[(table.family==ctx['family'])&(table.target==cell['target'])&(table.seed==cell['seed'])&(table.candidate==arm)].iloc[0];assert digest(count)==r.count_sha256;assert count.shape==(400,18533) and np.issubdtype(count.dtype,np.integer) and (count.data>=0).all();np.testing.assert_array_equal(count.sum(1).A1,x[cell['template']].sum(1).A1)
      if arm=='v24':base=count
      else:
       np.testing.assert_array_equal(count.sum(0).A1,base.sum(0).A1);libs=count.sum(1).A1;drift=float(abs(np.asarray(count.multiply((1/libs)[:,None]).mean(0)).ravel()-np.asarray(base.multiply((1/libs)[:,None]).mean(0)).ravel()).sum());assert abs(drift-r.mean_reference_L1)<1e-10
     nc+=1
    # Replay one null and one treatment population independently, full gene axis.
    if key=='H1' and ctx['family']=='H1':
     cell=cells[0];null=stage=='null';ids=cell['fit_controls'] if null else ctx['left'];fitids=cell['calibration_fit'] if null else [i for i in ids if i not in set(cell['template'])];p,bulk=moments(x[ids])[1:3];query=p*1e6;query[~measured]=np.nan;effect=np.zeros(len(genes)) if null else baseline(data,query,ctx['excluded'])[np.flatnonzero(targets==cell['target'])[0]];own=None if null else int(np.flatnonzero(genes==cell['target'])[0]);ti=0 if null else int(np.flatnonzero(targets==cell['target'])[0]);seed=cell['seed']+ti;masked=load_model(key);panel=load_model(key,True);state=state_values(x,m,ids,masked['state_axis'],masked['state_basis'],masked['state_scale']);lookup={int(v):i for i,v in enumerate(ids)};state=state[[lookup[int(v)] for v in cell['template']]];template=x[cell['template']];libs=template.sum(1).A1;atlas,_=load_atlas_model();base,_=old_emit(template,effect,own,seed,'v24',(libs,p,bulk),atlas,masked,state);count,diag=emit(template,effect,own,p,seed,base,FINALIST,panel,state,x[fitids]);r=table[(table.family==ctx['family'])&(table.target==cell['target'])&(table.seed==cell['seed'])&(table.candidate==FINALIST)].iloc[0];assert digest(count)==r.count_sha256;dr=de(count,x[cell['reference'] if null else ctx['right']],genes);assert int(((dr['q']<.05)&dr['eligible']).sum())==r.null_DE_calls;replays.append({'stage':stage,'family':key,'target':cell['target'],'count_sha256':digest(count)});print('V27 full-axis replay passed',stage,flush=True)
    del x;gc.collect()
 sampling=pd.read_csv(A/'sampling_scores.csv');old=pd.read_csv(ROOT/'experiments/v26_count_decoder/analysis/null_H1_scores.csv');assert sampling.loc[sampling.stage=='CPM_repaired','count_sha256'].iloc[0]==old.loc[old.candidate=='empirical','count_sha256'].iloc[0]
 decision=json.loads((A/'release_decision.json').read_text()) if (A/'release_decision.json').exists() else None
 if decision:
  assert not decision['released'] and not decision['new_package_generated'] and not decision['reserve_opened'];assert not (A.parent/'prediction.vcc').exists()
  for stage in ['null','treatment']:
   r=json.loads((A/f'{stage}_{decision["failed_family"]}_result.json').read_text());assert r['passed']==all(r['checks'].values());assert decision['result_hashes'][stage]==b.sha256(A/f'{stage}_{decision["failed_family"]}_result.json')
 write(A/('independent_training_input_audit.json' if inputs_only else 'independent_audit.json'),{'passed':True,'training_folds':fitted,'original_source_contexts_checked':ns,'original_full_axis_rows_checked':nr,'fit_only_group_means_recomputed':nm,'paired_cases_checked':nc,'full_axis_replays':replays,'V26_empirical_digest_reproduced':True,'H1_training_arrays_equal_immutable_V26':True,'plan_sha256':b.sha256(A/'plan.json'),'auditor_sha256':b.sha256(Path(__file__)),'scope':'Source exclusions/counts/fit-only means, paired cell identities, saved integer populations and two final-emitter replays. No independent new biological study or leaderboard improvement is asserted.'})
 print('V27 independent audit passed',flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--inputs-only',action='store_true');args=p.parse_args()
 with threadpool_limits(limits=2):main(args.inputs_only)
