"""V89 release test of the unsubmitted, previously audited V86 interval projection."""
from pathlib import Path
import importlib.util,shutil,argparse
import numpy as np
from threadpoolctl import threadpool_limits
from v57_pie import read,write,sha
import v86_joint_posterior as old
import v88_mixture_quantiles as fast
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v89_v86_projection/analysis';ARMS=('interval25','shared25','interval25_shuffled');data=old.data

def check():
 p=read(A/'plan.json')
 for n,h in p['code_sha256'].items():assert sha(ROOT/n)==h,n
 for n,h in p['V86_records_sha256'].items():assert sha(ROOT/n)==h,n
 return p

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists()
 names=['v89_projection_release.py','run_v89_counts.py','audit_v89_counts.py','v86_joint_posterior.py','v88_mixture_quantiles.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py'];records=['experiments/v86_joint_posterior/analysis/plan.json','experiments/v86_joint_posterior/analysis/model_audit.json','experiments/v86_joint_posterior/analysis/mean_results.json','experiments/v86_joint_posterior/analysis/count_assessment_0.json']
 write(A/'plan.json',dict(version=89,candidate='interval25',hypothesis='Release-test the previously unsubmitted V86 projection after V88 converged mixture quantiles failed donor PDS and JAC. No new response method is claimed. Original V86 seed0 improved5of6 metrics with unchanged PDS; V89 explicitly reuses that known result and prospectively tests candidate seed620000.',method='Exact V86 component kernels,100fixed EM updates and normal-moment interval25 projection,original quarter-strength bounded correction and80%V66 guard. The finite-iteration fit and normal interval are retained as an empirically conservative regularizer, not claimed converged or calibrated. V88 numerical optimization improved likelihood, not guaranteed transfer accuracy. Fast Cholesky is used solely for source likelihood arithmetic after parity checks.',public='Reuse all655 V86 public predictions and recompute the3arm summary. Whole study exclusions inherited and verified. Broad bound profile<=1.01,worststudy<=1.03,PDSloss<=.01. Strong support still requires original .5%profile gain and superiority to both controls; known V86 fails those. Nothing is called independently validated.',donor='Reuse all56 audited V86 seed0 records without new cells. Generate42 new primary/shared/shuffled populations for original14targets at offset620000; reuse14 identical V84 baseline populations. Same old donor/cohort/controls/libraries. Two seeds share biology; previous donor results motivated selection.',experimental_release='Sole candidate interval25, never promote controls. Both seeds must have profile<=1.005baseline,PDSloss<=.005,FIDloss<=.002,JAC loss<=.0005,NMAE<=1.01baseline and at least2of6 endpoints improve. Broad bounds also required. No threshold changed after second-seed result. Other failures and weak shuffled controls disclosed. No private gain forecast.',production='Only after release decision permits: fit same finite-iteration prior on all source families in a separate module namespace/output folder; retain original V66 A/B/C baseline effects,control templates,libraries andseeds; unchanged residual10countgenerator. Generate360000cells,18533genes,official prep and independent full payload/source/library checks and12fullgroup replays. No automatic upload or push.',code_sha256={n:sha(ROOT/n) for n in names},V86_records_sha256={n:sha(ROOT/n) for n in records},known_V88_seed0='PDS.56043956 versusV66.56593407;JAC.15151987 versus.15241833;release bounds fail. This motivated choosing the already known V86 candidate, not claiming fresh discovery.'))
 print('V89 follow-up protocol frozen',flush=True)

def adopt():
 check();dest=A/'predictions';dest.mkdir(exist_ok=True);items=[]
 for src in sorted((old.A/'predictions').glob('*.npz')):
  r=read(src.with_suffix('.json'));assert sha(src)==r['sha256'];out=dest/src.name;assert not out.exists();shutil.copyfile(src,out);assert sha(out)==r['sha256'];write(out.with_suffix('.json'),dict(r,reused_from_V86=str(src.relative_to(ROOT)),scientific_source_unchanged=True));items.append(dict(file=out.name,sha256=sha(out)))
 assert len(items)==6
 for name in ['model_audit.json','synthetic_tests.json']:
  assert read(old.A/name)['passed'];write(A/name,dict(read(old.A/name),reused_from_V86=True))
 counts=A/'counts/0';counts.mkdir(parents=True,exist_ok=True)
 for src in (old.A/'counts/0').glob('*'):
  if src.suffix not in ['.json','.npz']:continue
  out=counts/src.name;assert not out.exists();shutil.copyfile(src,out);assert sha(src)==sha(out)
 write(A/'adoption_receipt.json',dict(predictions=items,seed0_count_files_reused=112,no_new_seed0_cells=True,original_V86_files_unchanged=True));print('V89 adopted verified V86 predictions and seed0 results',flush=True)

_IMPL=None
def implementation():
 global _IMPL
 if _IMPL is None:
  check();spec=importlib.util.spec_from_file_location('v89_isolated_response',ROOT/'v86_joint_posterior.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);module.A=A;module.check=check;module.component_likelihood=fast.gaussian;_IMPL=module
 return _IMPL

def fit(outer,shuffle=False):assert outer=='production' and not shuffle;implementation().fit('production',False)
def predict(query,base,bulk,outer,controls=True):
 assert outer=='production' and not controls
 # Only the original primary is deployed; all original outputs are computed by the immutable method.
 pred,tr=implementation().predict(query,base,bulk,outer);return {'interval25':pred['interval25']},{'interval25':tr['interval25']}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','adopt','fit-production']);a=p.parse_args()
 with threadpool_limits(limits=2):{'freeze':freeze,'adopt':adopt,'fit-production':lambda:fit('production')}[a.stage]()
