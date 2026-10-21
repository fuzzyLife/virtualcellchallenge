"""V90 explicitly post-hoc diagnostic release of the tested V86 interval25 configuration."""
from pathlib import Path
import argparse,importlib.util
import numpy as np
from threadpoolctl import threadpool_limits
import v86_joint_posterior as old
from v88_mixture_quantiles import gaussian
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v90_interval_diagnostic/analysis'

def check():
 p=read(A/'plan.json')
 for n,h in p['code_sha256'].items():assert sha(ROOT/n)==h,n
 for n,h in p['known_evidence_sha256'].items():assert sha(ROOT/n)==h,n
 return p

def freeze():
 A.mkdir(parents=True,exist_ok=True);assert not (A/'plan.json').exists();prior=ROOT/'experiments/v89_v86_projection/analysis';d=read(prior/'release_decision.json');assert not d['experimental_release_allowed'];seeds={}
 for offset in [0,620000]:
  r=read(prior/f'count_assessment_{offset}.json');b=r['summary']['v66'];v=r['summary']['interval25'];assert v['profile_MSE']<b['profile_MSE'] and v['fidelity']>b['fidelity'] and v['nmae']<b['nmae'];assert v['pds']==b['pds'];seeds[str(offset)]=dict(profile_ratio=v['profile_MSE']/b['profile_MSE'],PDS_change=v['pds']-b['pds'],FID_change=v['fidelity']-b['fidelity'],JAC_change=v['jaccard']-b['jaccard'],reach_change=v['reach']-b['reach'],NMAE_change=v['nmae']-b['nmae'])
 for name in ['count_audit_0.json','count_audit_620000.json']:assert read(prior/name)['passed']
 evidence=['experiments/v86_joint_posterior/analysis/plan.json','experiments/v86_joint_posterior/analysis/model_audit.json','experiments/v89_v86_projection/analysis/mean_results.json','experiments/v89_v86_projection/analysis/count_assessment_0.json','experiments/v89_v86_projection/analysis/count_assessment_620000.json','experiments/v89_v86_projection/analysis/count_audit_0.json','experiments/v89_v86_projection/analysis/count_audit_620000.json','experiments/v89_v86_projection/analysis/release_decision.json','experiments/v88_mixture_quantiles/analysis/release_decision.json']
 names=['v90_diagnostic_release.py','fit_v90_submission.py','build_v90_submission.py','package_v90_submission.py','validate_v90_submission.py','audit_v90_production.py','v86_joint_posterior.py','v88_mixture_quantiles.py','v16_methods.py','v63_identity_guard.py']
 write(A/'plan.json',dict(version=90,candidate='interval25',status='POST_HOC_DIAGNOSTIC_NOT_VALIDATED_IMPROVEMENT',purpose='The user explicitly requests another leaderboard test. Release the unchanged V86 interval25 method on its stable profile/FID/NMAE improvements and unchanged donor PDS, despite known overlap tradeoffs. This is not a new model, prospective validation, passed V88/V89 research gate or a forecast of hidden-score improvement. Original failed protocols remain intact.',selection='Post hoc choice using both known donor seeds and broad public means; no independent accuracy claim. No new candidate/strength or count decoder is selected. Public evidence does not reliably rank hidden-context predictions.',method='Original V86 finite100-update covariance-mixture source prior and normal-moment interval25 projection with boundedquarter correction and80%V66 identity guard. Fit all five measured source families for production in an isolated module namespace. Fast float64 Cholesky only, after actual-source parity checks. Normal-moment intervals and finite optimization remain approximations.',known_tradeoffs=seeds,strong_support=False,known_V89_failed_checks=d['seeds'],integrity='Keep V66 A/B/C baseline,originalcontroltemplates/libraries/seeds,unchangedresidual10generator. Only complete official prep plus independent full payload/source/library and12group regeneration checks establish file readiness. These do not establish biological accuracy.',automatic_upload=False,code_sha256={n:sha(ROOT/n) for n in names},known_evidence_sha256={n:sha(ROOT/n) for n in evidence}))
 write(A/'release_decision.json',dict(candidate='interval25',strong_support=False,experimental_release_allowed=True,diagnostic_only=True,post_hoc_selection=True,prior_V88_V89_gates_failed=True,all_known_tradeoffs=seeds,known_best_V66_overall=.15088218754984492,predicted_leaderboard_score=None,reason='Stable gains in three donor metrics and unchanged PDS warrant a user-requested leaderboard diagnostic; JAC and weak shuffled controls prevent a validated-superiority claim.',automatic_upload=False));print('V90 diagnostic release frozen; original failed gates retained',flush=True)

_IMPL=None
def implementation():
 global _IMPL
 if _IMPL is None:
  check();spec=importlib.util.spec_from_file_location('v90_isolated_response',ROOT/'v86_joint_posterior.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);module.A=A;module.check=check;module.component_likelihood=gaussian;_IMPL=module
 return _IMPL

def fit(outer,shuffle=False):assert outer=='production';implementation().fit(outer,shuffle)
def predict(query,base,bulk,outer,controls=True):
 assert outer=='production' and not controls;pred,tr=implementation().predict(query,base,bulk,outer);return {'interval25':pred['interval25']},{'interval25':tr['interval25']}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','fit-production']);p.add_argument('--shuffled',action='store_true');a=p.parse_args()
 with threadpool_limits(limits=2):freeze() if a.stage=='freeze' else fit('production',a.shuffled)
