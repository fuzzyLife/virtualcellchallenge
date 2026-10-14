"""Freeze V24 conditional mean tests before their count outcomes."""
from pathlib import Path
import json
import numpy as np,pandas as pd
pd.options.future.infer_string=False
import build_submission as b
from v24_methods import ROOT,A,NAMES,write

def main():
 A.mkdir(parents=True,exist_ok=True);dest=A/'plan.json'
 if dest.exists():raise FileExistsError('Preserve V24 frozen plan')
 old=json.loads((ROOT/'experiments/v22_state_response/analysis/plan.json').read_text());contexts=[dict(family=c['family'],left=c['left'],right=c['right'],cells=c['cells'],excluded=[c['family']],tier='400-cell family holdout') for c in old['contexts']]
 folder=ROOT/'experiments/v22_state_response/analysis/cd4_D1_selected';m=pd.read_csv(folder/'cells.csv');left=[];right=[];good=[]
 for i,lane in enumerate(sorted(m.lane_id.unique())):
  ids=np.flatnonzero((m.gene_target=='Non-Targeting')&(m.lane_id==lane));ids=np.random.default_rng(24011+i).permutation(ids)
  if len(ids)>=40:left.extend(ids[::2].tolist());right.extend(ids[1::2].tolist());good.append(lane)
 target_names=sorted(t for t,n in m[m.lane_id.isin(good)].gene_target.value_counts().items() if t!='Non-Targeting' and n>=200)[:8];cells=[]
 for i,target in enumerate(target_names):
  rng=np.random.default_rng(24021+i);ids=np.flatnonzero((m.gene_target==target)&m.lane_id.isin(good));cells.append({'target':target,'seed':24021,'actual':rng.choice(ids,200,replace=False).tolist(),'template':rng.choice(left,200,replace=False).tolist()})
 assert len(cells)==8;contexts.append({'family':'CD4_Rest_D1','left':left,'right':right,'cells':cells,'excluded':[],'tier':'200-cell training-like donor calibration; public Rest summary biology overlaps'})
 deps=['v24_methods.py','prepare_v24_experiment.py','run_v24_experiment.py','assess_v24_experiment.py','test_v24_means.py','v16_methods.py','v15_ensemble_methods.py','v7_response_methods.py','v12_cpu_metrics.py','v12_methods.py','v8_methods.py','v6_methods.py','build_cd4_model.py','analyze_target_rank_correlations.py','fit_v7_response.py','build_submission.py','evaluate_atlas_transfer.py'];sources=['experiments/v7_multisource_adaptive/analysis/source_cache.npz','reference/cd4_probe/merged.npz','reference/cd4_probe/merged.csv','reference/cd4_probe/CE0006864_CE0008162.npz','reference/cd4_probe/CE0006864_CE0008162.csv','reference/cd4_probe/CE0008678_CE0010866.npz','reference/cd4_probe/CE0008678_CE0010866.csv']
 plan={'version':24,'contexts':contexts,'names':NAMES,'eligible_finalists':['cdf25','cdf50','cdf50_de'],'raw_guards':old['raw_guards'],'dependencies':{n:b.sha256(ROOT/n) for n in deps},'source_inputs':{n:b.sha256(ROOT/n) for n in sources},'validation_training_conditions':['Rest','Stim8hr'],'withheld_response_condition':'Stim48hr','mean_protocol':'Remove all legacy CD4 Stim48hr response rows from base-bank fitting and teacher during development/reserve. Within-state donor support and credible null evidence, all-transcriptome pairwise-valid rank compatibility versus other source controls, top512 training-control rank-change genes for within-CD4 state weights. No added1 in rank correlations. Mean correction alpha .25/.5 times compatibility and confidence, clipped +/- .5LFC before alpha. Own target fixed. V16 residual10 emitter unchanged.',
 'development':'28 genuine400-cell cases in three previously studied source families; eight200-cell Rest D1 cases are calibration and overlap public Rest summary biology. No claim of blind donor validation. Refit V16 baseline excludes48hr responses too; it is not the exact published V16 package. Family balanced, missing NMAE explicit.',
 'reserve':'Untouched D2 Stim48hr raw counts, eight alphabetically selected competition targets>=400, disjoint lane-stratified NTC halves. Condition response rows excluded from every fitting/selection stage. Existing study; upstream DE-stat preprocessing may share nuisance estimation, so this is exploratory condition transfer, not a fully independent study or whole-lineage blind test.',
 'release':'Old available six raw guards, >=3 improvements including DE; worst-family fidelity/reach guards. Eligible cdf25/cdf50/cdf50_de must beat their exactly alpha-matched target-shuffled correction on >=2 fidelity/reach/Jaccard metrics. No distribution-shape gate for a mean expert: this phase tests response directions, not learned population shape. The DE-only arm excludes null-equivalence corrections; alpha-matched shuffles retain target-wise amplitude differences and the separate amplitude comparator diagnoses generic scaling. Freeze one winner by equal available-metric improvement ranks; no selection/retuning after once-only reserve. Real full-axis null exactV16. Failed new experts retain champion.',
 'production':'After passing reserve, standard all-source refit may use all three public CD4 conditions. Start from exact V16 published ABC effect arrays; only compatible profiles and supported coordinates receive bounded corrections. Rank correlation implies compatibility, not literal context identity. Full new package counts/axes/provenance independently audited; no automatic upload.'}
 write(dest,plan);print('V24 plan frozen',b.sha256(dest),'CDF calibration',target_names,flush=True)
if __name__=='__main__':main()
