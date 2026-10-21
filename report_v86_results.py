"""Write the completed V86 experiment report without selecting a post-hoc winner."""
from importlib.metadata import version
import platform
import numpy as np
import v86_joint_posterior as m
from v57_pie import read,write,sha

def fmt(x,n=6):return 'unavailable' if x is None else f'{x:.{n}f}'

def main():
    done=read(m.A/'experiment_complete.json');means=read(m.A/'mean_results.json');count=read(m.A/'count_assessment_0.json');audit=read(m.A/'model_audit.json');assert done['completed'] and audit['passed'] and read(m.A/'count_audit_0.json')['passed']
    lines=['# V86 results: multivariate Bayesian response transfer','',
     f"**Original qualification passed: {done['strict_support']}.** The fixed primary is `interval25`; no other arm is promoted after seeing outcomes. No upload package was generated. V66/V82 remains the best published project result at **0.1508821875498449**.",'',
     '[Protocol and commands](V86_STRATEGY.md) | [Guide-contrast experiment](V85_GUIDE_CONTRAST.md)','',
     '## Different mechanisms actually tested','',
     'V85 used differences between guides as a source of on-target response evidence. Its corrected negative controls behaved properly, but every source target failed at least one fixed dose-evidence requirement, so it supplied exactly zero corrections. The original sparse-log failure, source-only correction and full first-stage summaries are preserved.','',
     'V86 learns a mixture of response-sharing covariances and computes a conditional posterior for an unseen context. It uses measurement masks and source uncertainties, excludes the entire recipient study, retains query-specific variation and tests a shared-only covariance comparator. This is an independent mash-inspired model, not a reproduction of mashr or a published competition solution.','',
     '## Synthetic and model checks','',
     f"Synthetic held-context error / copying one source: **{read(m.A/'synthetic_tests.json')['synthetic_unseen_context_MSE_ratio_to_copy']:.6f}**. This validates the implementation under its simulated assumptions; it is not real-data accuracy.",'',
     'Intervals use the conditional mixture mean and variance with a normal approximation; they are not exact mixture quantiles or established frequentist coverage. Mixture weights are plug-in estimates.','',
     'Six priors were fitted: genuine and independently target-shifted source evidence for H1, Orion and CD4 outer exclusions. Each uses 20,000 source target/output coordinates. These correlated coordinates are working-likelihood units, not independent biological observations.','',
     '| Outer exclusion | Shifted labels | EM iterations | Final maximum weight change |','| --- | --- | ---: | ---: |']
    for r in audit['models']:lines.append(f"| {r['outer']} | {r['shuffled']} | {r['EM_iterations']} | {r['final_weight_change']:.8f} |")
    lines+=['','All six fits reached the fixed 100-update limit before the 1e-7 weight-change stopping tolerance. This is a bounded pilot, not a fully converged fit. Complete source optimization is a remaining model-development step; no extra iterations or component settings were selected from recipient scores. Model metadata retain full weights and histories.','',
     'Independent checks replay Gaussian likelihoods against SciPy, conditional posteriors by direct solves, covariance positive semidefiniteness and mixture moments. Poisoning all excluded-study responses, masks, uncertainties, doses and profiles leaves every allowed source-bank field unchanged. All final prediction files also pass independent 80% identity-margin replay.','',
     '## All 655 reused public pairs','',
     '| Arm | Equal-study profile error / V66 | PDS change |','| --- | ---: | ---: |']
    for arm,r in means['summary'].items():lines.append(f"| {arm} | {r['profile_ratio']:.9f} | {r['pds_change']:+.6f} |")
    lines+=['','| Primary development check | Passed |','| --- | --- |']
    for key,val in means['checks'].items():lines.append(f'| {key} | {val} |')
    lines+=['','H1, Orion and CD4 receive equal study weight; within Orion/CD4, contexts are averaged. All original source exclusions remain enforced, including every CD4 state and donor for CD4 queries. The biology and public outcomes have been reused, so this is development evidence. Complete context-level results are retained in `analysis/mean_results.json`.','',
     '## D3 final-count diagnostic','',
     'The original metadata-fixed 14 targets use 400 distinct actual cells, 400 distinct NTC templates, original template library totals and a separate 4,994-cell NTC reference. D3 was examined in V84 and is not blind. Current V86 fitting/conditioning excludes all CD4 responses; recipient treatment values never enter its predictor.','',
     '| Arm | Profile error | PDS | Fidelity | Reach | Jaccard | NMAE |','| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for arm,r in count['summary'].items():lines.append('| '+arm+' | '+' | '.join(fmt(r[k]) for k in ['profile_MSE','pds','fidelity','reach','jaccard','nmae'])+' |')
    b=count['summary']['v66'];v=count['summary']['interval25']
    lines+=['',f"Primary profile error change: **{100*(v['profile_MSE']/b['profile_MSE']-1):+.4f}%**; PDS change **{v['pds']-b['pds']:+.6f}**.",'',
     '| Primary count check | Passed |','| --- | --- |']
    for key,val in count['checks'].items():lines.append(f'| {key} | {val} |')
    improved=[k for k in ['profile_MSE','pds','fidelity','reach','jaccard','nmae'] if (v[k]<b[k] if k in ['profile_MSE','nmae'] else v[k]>b[k])]
    lines+=['',f"The fixed primary improves {len(improved)} of the six raw point estimates: {', '.join(improved)}. Discrimination is unchanged. These gains are small: the shared-only model also improves expression, and the target-shuffled model has still lower expression error. The primary therefore fails the shuffled-expression requirement despite its favorable DE point estimates.",'']
    lines+=['','Raw public statistics above are not hosted anchor-scaled scores and are not averaged into an invented competition overall. These targets and this donor were previously used. A single technical seed is not independent biological replication.','',
     f"Second seed executed: **{done['second_seed_run']}**. It is conditional on passing both original development and first count requirements. No failed thresholds were relaxed.",'',
     '## Baseline reuse and numerical integrity','',
     'The test emitted **42 new populations** and reused **14 previously audited V66 populations**, for 56 compared records. Exact equality checks cover baseline response, CPM, bulk proportions, axes, cohort and seed; unchanged emitter/metric code is checked. Current truth/reference profiles are replayed before copying baseline evidence. The reused baseline is not claimed as newly generated data.','',
     'Emission checks 336 new raw tied-rank coordinates; 112 baseline coordinates were checked previously. Independent arithmetic replay verifies BH, directional scores, profile errors and 14-target discrimination for every one of the 56 records. Integer nonnegative counts and exact original row library totals are asserted.','',
     '## Decision and reproducibility','',
     'The fixed primary must pass development, count noninferiority and gain requirements, beat both controls, and pass a separately confirmed second seed before a release recommendation. Novelty and synthetic success do not satisfy those requirements. No new upload file or better hosted score is claimed.','',
     'Run `finish_v86_experiment.py` in a fresh analysis checkout with the committed plan and required prior data. Individual stages and numerical audits are in `V86_STRATEGY.md`. `report_v86_results.py` regenerates this report from the completed, hashed receipts. Existing scientific releases and their arrays remain unchanged.']
    (m.ROOT/'V86_RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    write(m.A/'release_decision.json',dict(primary='interval25',strict_support=done['strict_support'],development_passed=means['passes'],first_count_passed=count['passes'],second_seed_run=done['second_seed_run'],new_upload_generated=False,model_audit_sha256=sha(m.A/'model_audit.json'),count_audit_sha256=sha(m.A/'count_audit_0.json'),report_script_sha256=sha(__file__),existing_best_published=.1508821875498449,hosted_score_forecast=None))
    write(m.A/'runtime_versions.json',dict(python=platform.python_version(),packages={p:version(p) for p in ['numpy','scipy','pandas','scikit-learn','h5py','anndata','threadpoolctl']}))
    print('V86 final report generated, strict support',done['strict_support'],flush=True)
if __name__=='__main__':main()
