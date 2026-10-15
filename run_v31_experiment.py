"""Freeze V31 guide-aware means, calibrate only in D2, then test public development."""
from pathlib import Path
import argparse,json,sys,importlib.metadata
import numpy as np
from threadpoolctl import threadpool_limits
import v31_methods as m
import build_submission as b
DEPS=['v31_methods.py','run_v31_experiment.py','fit_v7_response.py','v25_methods.py','v24_methods.py','v16_methods.py','v15_ensemble_methods.py','v7_response_methods.py','v8_methods.py','build_submission.py']

def freeze():
 assert not (m.A/'plan.json').exists();old=json.loads((m.ROOT/'experiments/v25_cross_domain/analysis/plan.json').read_text());contexts=[]
 for c in old['contexts']:
  if c['family']=='CD4_Stim48_D2':continue
  d=dict(c);d['excluded']=['CD4'] if d['family'].startswith('CD4') else ['HCT116','HEK293T'] if d['family'] in ['HCT116','HEK293T'] else ['H1'];contexts.append(d)
 inputs=['experiments/v7_multisource_adaptive/analysis/source_cache.npz','experiments/v24_cd4_response/analysis/models/teacher_all.npz']
 for folder in ['experiments/v22_state_response/analysis/cd4_selected','experiments/v24_cd4_response/analysis/holdout','experiments/v22_state_response/analysis/cd4_D1_selected']:
  inputs.extend(folder+'/'+f for f in ['counts.npz','cells.csv','genes.csv','extraction.json'])
 m.write(m.A/'plan.json',{'version':31,'dependencies':{p:b.sha256(m.ROOT/p) for p in DEPS},'inputs':{p:b.sha256(m.ROOT/p) for p in inputs},'contexts':contexts,'arms':m.ARMS,'sources':m.CONDITIONS,'alpha':[.125,.25,.5],'calibration':'Choose alpha using whole-guide-held-out D2Rest/Stim only; baseline excludes all CD4 summary effects. Shared lane controls are permitted and their reuse is explicit; this is within-study calibration. Never choose alpha on D1/H1/Orion development outcomes.','model':'Empirical per-guide CPM mean/variance with matched lane NTCs. Source KD prior0.9, gainclip0.8..1.8. Random-effects guide pooling retains heterogeneity and bounded inverse variance. Single-guide estimates allowed with sampling uncertainty. Pooled control recomputes actual all-guide raw mean. Guide shuffles preserve target/lane group sizes.','validation':'H1 excluded; HCT116+HEK293T excluded together. CD4D1 donor transfer uses raw D2 teachers but excludes all CD4 summary evidence from baseline; this is a donor holdout, NOT a whole-CD4 or unseen-lineage holdout. Source D2 stimulated cases are never validation cases. These public development cells are previously spent; no new biological blind claim.','mean_metrics':'Candidate and baseline mean CPM from query controls times2**effect; own transcript factor0.1 for both, then total normalization. Compare clipped+/-3 smoothed0.1 log2FC vs independent reference controls and measured actual populations, excluding own output. Approximate mean screen, not count/Wilcoxon or anchored leaderboard score.','mean_gates':'Guide MSE ratio<=0.99 and MAE<=0.995 across equally weighted H1/Orion/CD4D1. CD4D1 MSE<=0.95; all family MSE<=1.02. Must beat actual pooled, within-target/lane guide-shuffled, and target-shuffled controls on BOTH MSE and MAE. Mean CPM squared-error ratio<=1.01.','release':'No count test for a failed mean candidate. A passing candidate requires separately frozen400-cell DE/count guards and once-only fresh D3 donor validation before package. No upload. Conditional transport/foundation adaptation remains contingent on useful mean response evidence.','runtime':{'python':sys.version,'versions':{p:importlib.metadata.version(p) for p in ['numpy','scipy','pandas','anndata','h5py','threadpoolctl']}}});print('V31 plan frozen',b.sha256(m.A/'plan.json'),flush=True)

def normmean(x):
 lib=x.sum(1).A1;assert (lib>0).all();return np.asarray(x.astype(float).multiply((1e6/lib)[:,None]).mean(0)).ravel()

def mean_screen():
 p=m.load_plan();rows=[];genes=m.data()['genes'];targets=m.data()['targets'];teacher,rec=m.load_teacher()
 for ctx in p['contexts']:
  family=ctx['family'];x,meta,measured,source=m.raw(family,genes);query=normmean(x[ctx['left']]);query[~measured]=np.nan;reference=normmean(x[ctx['right']]);effects={};routing={};parent=None
  for arm in m.ARMS:effects[arm],parent,routing[arm]=m.predict(query,ctx['excluded'],arm)
  effects['v24_reference']=parent
  for cell in ctx['cells']:
   target=cell['target'];ti=int(np.flatnonzero(targets==target)[0]);own=int(np.flatnonzero(genes==target)[0]);actual=normmean(x[cell['actual']]);valid=measured&(reference>5);valid[own]=False;truth=np.clip(np.log2((actual+.1)/(reference+.1)),-3,3);values={}
   for arm,effect in effects.items():
    predicted=np.nan_to_num(query)*np.exp2(effect[ti]);predicted[own]*=.1;predicted*=1e6/predicted.sum();lfc=np.clip(np.log2((predicted+.1)/(reference+.1)),-3,3);delta=lfc[valid]-truth[valid];values[arm]={'squared':float(delta@delta),'absolute':float(abs(delta).sum()),'cpm_squared':float(((predicted[valid]-actual[valid])**2).sum())}
   base=values['v24_reference']
   for arm in m.ARMS:rows.append({'family':family,'study':'Orion' if family in ['HCT116','HEK293T'] else 'CD4_D1' if family.startswith('CD4') else 'H1','target':str(target),'seed':cell['seed'],'actual_cells':len(cell['actual']),'arm':arm,'mse_ratio':values[arm]['squared']/max(base['squared'],1e-12),'mae_ratio':values[arm]['absolute']/max(base['absolute'],1e-12),'cpm_mse_ratio':values[arm]['cpm_squared']/max(base['cpm_squared'],1e-12)})
  print('V31 mean context',family,'cases',len(ctx['cells']),'gate',routing['guide']['compatibility_gate'],flush=True)
 study={};summary={}
 for arm in m.ARMS:
  study[arm]={}
  for key in ['H1','Orion','CD4_D1']:
   rr=[r for r in rows if r['study']==key and r['arm']==arm];assert rr;study[arm][key]={k:float(np.mean([r[k] for r in rr])) for k in ['mse_ratio','mae_ratio','cpm_mse_ratio']}
  summary[arm]={k:float(np.mean([r[k] for r in study[arm].values()])) for k in ['mse_ratio','mae_ratio','cpm_mse_ratio']}
 q=summary['guide'];checks={'mse_gain':q['mse_ratio']<=.99,'mae_gain':q['mae_ratio']<=.995,'donor_mse_gain':study['guide']['CD4_D1']['mse_ratio']<=.95,'worst_family':max(r['mse_ratio'] for r in study['guide'].values())<=1.02,'expression_noninferiority':q['cpm_mse_ratio']<=1.01}
 for control in ['pooled','guide_shuffled','target_shuffled']:checks['beats_'+control]=all(q[k]<summary[control][k] for k in ['mse_ratio','mae_ratio'])
 m.write(m.A/'mean_screen.json',{'plan_sha256':b.sha256(m.A/'plan.json'),'teacher_sha256':rec['sha256'],'alpha':rec['chosen_alpha'],'records':rows,'routing':routing,'study_summary':study,'summary':summary,'checks':checks,'new_candidate_passes':all(checks.values()),'scope':'Approximate mean-composition screen only; no generated cells, Wilcoxon metrics, fresh reserve or package.'});print(json.dumps({'summary':summary,'studies':study,'checks':checks},indent=2),flush=True)

def main():
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','prepare','summary']);a=p.parse_args()
 with threadpool_limits(limits=2):{'freeze':freeze,'prepare':m.prepare,'summary':mean_screen}[a.stage]()
if __name__=='__main__':main()
