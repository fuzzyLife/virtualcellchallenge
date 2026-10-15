"""Independent V29 count/margin/identity audit and deterministic full-axis replays."""
from pathlib import Path
import json,gc
import numpy as np,pandas as pd
pd.options.future.infer_string=False
from scipy import sparse
from threadpoolctl import threadpool_limits
from v29_methods import ROOT,A,FINALIST,arms_for,load_plan,emit,write,digest
from v26_methods import raw,moments,baseline,emit as old_emit,family_key
from fit_v7_response import read_cache
from evaluate_atlas_transfer import load_atlas_model
from v12_cpu_metrics import de
import build_submission as b

def main():
 plan=load_plan();data=read_cache();genes=data['genes'];targets=data['targets'];nc=0;replays=[];checked=[];rules=plan['null_policy'];np.testing.assert_equal(len(genes),18533)
 oldaudit=ROOT/'experiments/v27_control_anchored/analysis/independent_training_input_audit.json';old=json.loads(oldaudit.read_text());assert old['passed'] and old['original_full_axis_rows_checked']==128 and old['fit_only_group_means_recomputed']==662;assert old['plan_sha256']==b.sha256(ROOT/'experiments/v27_control_anchored/analysis/plan.json');assert old['auditor_sha256']==b.sha256(ROOT/'audit_v27_experiment.py')
 for key in ['H1','HCT116','HEK293T','CD4']:
  for stage in ['null','treatment']:
   ep=A/f'{stage}_{key}_evaluation.json'
   if not ep.exists():continue
   ev=json.loads(ep.read_text());table=pd.read_csv(A/f'{stage}_{key}_scores.csv');assert ev['complete'];assert ev['scores_sha256']==b.sha256(A/f'{stage}_{key}_scores.csv') and ev['panels_sha256']==b.sha256(A/f'{stage}_{key}_panels.csv');assert set(table.candidate)==set(arms_for(key));checked.append(stage+'_'+key)
   if stage=='null':
    expected={'mean_contract':bool((table[table.candidate==FINALIST].mean_reference_L1<=.0005).all())}
    for _,r in table[table.candidate==FINALIST].iterrows():
     pr=str(r.family)+'_'+str(int(r.seed));expected[pr+'_DE_calls']=bool(r.null_DE_calls<=max(rules['false_call_floor'],r.real_null_DE_calls+rules['real_null_call_allowance']));expected[pr+'_detection']=bool(r.detection_MAE<=max(rules['detection_MAE_floor'],rules['detection_real_multiplier']*r.real_detection_MAE+rules['detection_real_allowance']));expected[pr+'_variance']=bool(r.log_variance_MAE<=max(rules['log_variance_MAE_floor'],rules['log_variance_real_multiplier']*r.real_log_variance_MAE+rules['log_variance_real_allowance']))
    rr=A/f'{stage}_{key}_result.json'
    if rr.exists():
     result=json.loads(rr.read_text());assert result['checks']==expected and result['passed']==all(expected.values())
   for ctx in [c for c in plan['contexts'] if family_key(c['family'])==key]:
    x,m,measured,_=raw(ctx['family'],genes);cells=ctx['null_cases'] if stage=='null' else ctx['cells'];folder=A/'parallel'/f'{stage}_{key}'/ctx['family']
    for cell in cells:
     null=stage=='null';ids=cell['fit_controls'] if null else ctx['left'];refids=cell['reference'] if null else ctx['right'];fitids=cell['calibration_fit'] if null else [i for i in ids if i not in set(cell['template'])];sets=[set(cell['actual']),set(refids),set(cell['template']),set(fitids)];assert all(not sets[i]&sets[j] for i in range(4) for j in range(i));assert len(sets[0])==len(cell['actual']) and len(sets[2])==400 and len(sets[3])>=200;assert (m.gene_target.to_numpy(str)[cell['actual']]==('Non-Targeting' if stage=='null' else cell['target'])).all();assert (m.gene_target.to_numpy(str)[list(sets[1]|sets[2]|sets[3])]=='Non-Targeting').all();cache=json.loads((folder/f'{cell["target"]}_{cell["seed"]}.json').read_text());assert cache['sha256']==b.sha256(folder/f'{cell["target"]}_{cell["seed"]}.npz') and cache['calibration_fit_rows']==fitids
     for arm in ['v24',FINALIST]:
      count=sparse.load_npz(folder/f'{cell["target"]}_{cell["seed"]}_{arm}_counts.npz');r=table[(table.family==ctx['family'])&(table.target==cell['target'])&(table.seed==cell['seed'])&(table.candidate==arm)].iloc[0];assert digest(count)==r.count_sha256;assert count.shape==(400,18533) and np.issubdtype(count.dtype,np.integer) and (count.data>=0).all();libs=x[cell['template']].sum(1).A1;np.testing.assert_array_equal(count.sum(1).A1,libs)
      if arm=='v24':base=count
      else:
       np.testing.assert_array_equal(count.sum(0).A1,base.sum(0).A1);drift=float(abs(np.asarray(count.multiply((1/libs)[:,None]).mean(0)).ravel()-np.asarray(base.multiply((1/libs)[:,None]).mean(0)).ravel()).sum());assert abs(drift-r.mean_reference_L1)<1e-10
     if stage=='null' and ctx['family']=='CD4_Stim48_D2' and cell['seed']==270672:
      probe=json.loads((A/'compiled_probe.json').read_text());assert digest(count)==probe['count_sha256']
     nc+=1
    # H1 first null and first treatment, plus the first later-family null, are
    # rebuilt from original NTC inputs and frozen means rather than saved counts.
    if (key=='H1' and ctx['family']=='H1') or (stage=='null' and key=='HCT116'):
     cell=cells[0];null=stage=='null';ids=cell['fit_controls'] if null else ctx['left'];refids=cell['reference'] if null else ctx['right'];p,bulk=moments(x[ids])[1:3];query=p*1e6;query[~measured]=np.nan;ti=0 if null else int(np.flatnonzero(targets==cell['target'])[0]);effect=np.zeros(len(genes)) if null else baseline(data,query,ctx['excluded'])[ti];own=None if null else int(np.flatnonzero(genes==cell['target'])[0]);seed=cell['seed']+ti;template=x[cell['template']];libs=template.sum(1).A1;atlas,_=load_atlas_model();base,_=old_emit(template,effect,own,seed,'v24',(libs,p,bulk),atlas);count,diag=emit(template,effect,own,p,seed,base,FINALIST,None,None);r=table[(table.family==ctx['family'])&(table.target==cell['target'])&(table.seed==cell['seed'])&(table.candidate==FINALIST)].iloc[0];assert digest(count)==r.count_sha256;d=de(count,x[refids],genes);assert int(((d['q']<.05)&d['eligible']).sum())==r.null_DE_calls
     if key=='H1' and null:
      from unittest.mock import patch
      from v29_pairwise_reference import repair_pairwise as python_reference
      with patch('v29_methods.repair_pairwise',python_reference):
       replay,_=emit(template,effect,own,p,seed,base,FINALIST,None,None)
      assert digest(replay)==digest(count),'Compiled/full Python reference differs'
     replays.append({'family':key,'stage':stage,'target':cell['target'],'count_sha256':digest(count)});print('V29 full-axis replay passed',key,stage,flush=True)
    del x;gc.collect()
 decision=A/'release_decision.json'
 if decision.exists():
  r=json.loads(decision.read_text());assert not r['released'] and not r['new_package_generated'] and not r['reserve_opened'];assert not (A.parent/'prediction.vcc').exists() and not (A/'reserve_execution_started.json').exists()
 write(A/'independent_audit.json',{'passed':True,'evaluations_checked':checked,'paired_cases_checked':nc,'full_axis_replays':replays,'null_checks_independently_recomputed':True,'H1_full_axis_compiled_Python_reference_parity':True,'immutable_H1_training_audit_sha256':b.sha256(oldaudit),'pretrained_panel_sha256':b.sha256(ROOT/'experiments/v27_control_anchored/analysis/models/H1_panel.npz'),'plan_sha256':b.sha256(A/'plan.json'),'auditor_sha256':b.sha256(Path(__file__)),'scope':'Frozen inputs, disjoint actual/reference/template/fit NTCs, saved full-axis integer populations/margins, null criteria and replay counts. Source neural training reused unchanged from the128-row/662-fit-mean audited H1 exclusion. Executed scope only, not a claim of unseen-study or leaderboard improvement.'});print('V29 independent audit passed',flush=True)
if __name__=='__main__':
 with threadpool_limits(limits=2):main()
