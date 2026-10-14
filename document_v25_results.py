"""Reproduce the V25 negative-result report and retain the published champion."""
from pathlib import Path
import json,sys,importlib.metadata
import pandas as pd
pd.options.future.infer_string=False
import build_submission as b
from v25_methods import ROOT,A,load_plan,write

def main():
 plan=load_plan();selection=json.loads((A/'selection.json').read_text());audit=json.loads((A/'independent_experiment_audit.json').read_text());assert audit['passed'] and not selection['new_candidate_passes'] and selection['candidate']=='v24';assert not (A/'reserve_execution_started.json').exists() and not (A/'confirmation_result.json').exists();f=pd.read_csv(A/'development_summary.csv').set_index('candidate');old=json.loads((ROOT/'experiments/v24_cd4_response/analysis/leaderboard_result.json').read_text());deferred=json.loads((A/'holdout_download_deferred.json').read_text(encoding='utf-8-sig'));assert not deferred['download_complete'] and not deferred['reserve_opened'];assert not (A.parent/'prediction.vcc').exists();readme=ROOT/'README.md';text=readme.read_text(encoding='utf-8');section='''<!-- v25-cross-domain:begin -->

## V25 cross-domain geometry tests (2026-10-14)

Verified the supplied conceptual posters against primary mathematical literature and tested intervention-aligned response geometry: bounded force/Bures transport, split-control agreement and a Fisher-proxy trust region. All three finalists failed the frozen development criteria across44 public cases/396 generated populations. The independent family-exclusion,20 source-covariance recomputations,44 cache checks and15 full-axis count replays passed. V24 remains the best published submission at0.14299785521745675; no V25 package is released and the donor reserve remains unopened. The new donor download was stopped after rejection, with partial ranges retained; no complete-file MD5 or raw-row audit is claimed. Metadata also identifies low-power CRISPRa fibroblast/RPE1 files, which were not added as CRISPRi labels.

Methods, literature, limits and commands: [V25_STRATEGY.md](V25_STRATEGY.md). Numeric results and evidence: [V25_RESULTS.md](V25_RESULTS.md). Artifact decision: [V25_SUBMISSION.md](V25_SUBMISSION.md). Root scripts and small evidence are versioned; prior numerical code and packages are preserved.

<!-- v25-cross-domain:end -->
'''
 marker='<!-- v25-cross-domain:begin -->';end='<!-- v25-cross-domain:end -->'
 if marker in text:text=text[:text.index(marker)]+section.strip()+text[text.index(end)+len(end):]
 else:
  lines=text.splitlines(keepends=True);text=lines[0]+'\n'+section+'\n'+''.join(lines[1:])
 readme.write_text(text,encoding='utf-8');rows=[]
 for name in ['v24','force25','force50','bures25','diagonal25','shuffled25','shuffled50','amplitude25','source_program25']:
  v=f.loc[name];rows.append('|'+name+'|'+'|'.join(f'{v[k]:.6f}' for k in ['pds','fidelity','reach','jaccard','nmae'])+'|')
 failures=[]
 for c in selection['comparisons']:failures.append('- **'+c['candidate']+'**: '+', '.join(k for k,v in c['checks'].items() if not v)+'.')
 report=f'''# V25 measured results: intervention-aligned geometry

**Retain V24.** All three predeclared finalists fail development. No production package was generated, no upload was attempted, and the new-donor reserve was not scored. This is a negative result for these bounded predictive models, not a proof that every geometry-based model fails.

The frozen plan SHA256 is `{b.sha256(A/'plan.json')}`. Four held-family fits and one all-source fit use64 training-intervention directions on4096 shared output genes. Development emits396 populations:44 public cases times9 arms, including eight200-cell CD4RestD1 cases and36 actual400-cell cases. Baseline is the V24 algorithm under the same source-family exclusions, not the published ABC package. The biology is reused development; two CD4 conditions share a study/family.

## Family-balanced raw metrics

These are local raw metrics, not the hosted anchor-scaled leaderboard components. Higher is better for the first four; lower NMAE is better. HCT116 NMAE is unavailable, so the NMAE column averages available families only.

|Arm|PDS|Fidelity|Reach|Jaccard|NMAE|
|---|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

Expression comparisons use paired per-panel ratios because simple pooled ratios are dominated by small control denominators. The mean raw expression ratio in development_summary.csv is not a hosted score. The primary expression guard passes, but none of the finalists improves three available primary metrics. Shape Wasserstein diagnostics are exploratory; this model changes response means with an unchanged residual10 emitter.

## Frozen release failures

{chr(10).join(failures)}

Force25 gains0.002008 fidelity but loses0.013475 reach and0.010000 PDS; its worst-family reach guard fails. Force50 gains0.020272 reach, while fidelity drops0.002351, PDS drops0.006667 and NMAE increases0.004629. Bures25 loses0.036298 reach and fails its shuffle/amplitude/program comparisons. The source-program control improves NMAE but loses discrimination and reach; it is a diagnostic arm, not an eligible alternative selected after seeing outcomes.

All original thresholds were retained. The conservative release_guard_amendment.json was written before opening outcome tables and can only reject further candidates. Because no finalist passes the original guards, its extra alpha0.5/Bures-specific diagnostics were not triggered. No fallback, post-hoc blend, threshold relaxation or re-selection was performed. A future experiment must freeze a new plan and retain this negative result.

## Independent verification

- Three mathematical/mechanism unit tests pass, including equal stationary covariance/Fisher geometry with different intervention responses.
- All frozen direct code/source hashes and the transitive V24 assessor hash pass.
- All five model axes, family exclusions,64-dimensional orthonormal bases, finite source shifts/weights and bounded corrections pass.
- {audit['independently_recomputed_source_covariances']} source covariances are independently recomputed from dense logCP10k counts with within-batch centering and log-depth removal.
- All {audit['frozen_cases']} case identities, actual/template uniqueness, disjoint NTC halves, cache hashes and saved count digests pass.
- {audit['full_axis_count_replays']} full18533-gene count replays match exactly, with exact template libraries; {audit['no_effect_emitter_V16_identity_cases']} no-effect checks match the unchanged V16 emitter. This is not literal real-control identity or null-DE calibration.
- Recomputed release comparisons match the saved decision exactly.

The first auditor attempt incorrectly equated projected HCT116 row sums with full-source library sizes. It stopped; the auditor was corrected to use verified full-source depths, with projected sums bounded above by those depths. No frozen predictor, fitted array, test population or threshold changed. The final audit passes. Source summaries and raw CD4 covariance controls are matched at condition/study-family level, not proven exact donor aliases; whole-CD4 exclusion prevents those related sources entering CD4 validation.

## New source evidence and reserve status

Author-release metadata finds79 fibroblast and34 RPE1 competition-target overlaps, with no overlapping target reaching200 or400 retained cells. Literal symbols were checked, so low overlap/power is not a namespace conversion issue. These are CRISPRa screens and were not treated as CRISPRi labels. They remain possible representation sources requiring separate cross-modality validation. Relevant primary sources and the poster-claim assessment are in V25_STRATEGY.md.

Marson D3Rest metadata has37 overlapping targets,14 with>=400 cells. The new raw file download was stopped after development rejected every finalist: {int(deferred['partial_bytes']):,} partial bytes in{deferred['range_files']} range files remain resumable. The full file has **not** been MD5-verified, extracted or independently row-audited; its reserve is unopened. Extraction/audit scripts are prepared, not claimed executed. No query treatment-expression metrics were used to choose this source.

## What this changes in our reasoning

Control co-fluctuation structure is insufficient as the only context-transfer rule, even when response directions come from interventions. A bigger correction trades reach against direction/amplitude error. A source-program projection can reduce fold-change error without improving discrimination, so it cannot be promoted on one favorable metric. Future tests should identify which source intervention programs transfer with context and assay, and separate mean-response error from generated-cell DE calibration. The supplied universal-isomorphism claim supplies no verified shortcut around those measurements.

## Reproduce and inspect

Use the pinned Python3.12 scientific runtime and requirements-v23-science.txt with two BLAS threads. V25_STRATEGY.md gives from-scratch commands. Existing plans, models and decision files intentionally refuse overwrite. For the saved run:

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
set "VCC_SCI_PY=C:\\Users\\animeshs\\AppData\\Local\\Temp\\vcc_v20_stack\\env\\Scripts\\python.exe"
"%VCC_SCI_PY%" test_v25_geometry.py
"%VCC_SCI_PY%" audit_v25_experiment.py
"%VCC_SCI_PY%" document_v25_results.py
```

Numeric tables and hash evidence are in experiments/v25_cross_domain/analysis. No V25 host score is claimed. Current published champion:{old['entry_id']}, overall{old['overall']}.
'''
 (ROOT/'V25_RESULTS.md').write_text(report,encoding='utf-8');(ROOT/'V25_SUBMISSION.md').write_text(f'''# V25 artifact decision

**No new V25 submission is released.** All three geometry finalists fail the frozen development criteria. The independent implementation/split/count audit passes; statistical usefulness does not. The donor reserve remains unopened. See V25_RESULTS.md and V25_STRATEGY.md.

Keep the already published V24 entry `{old['entry_id']}`, overall **{old['overall']}**. Its local package is `experiments/v24_cd4_response/prediction.vcc`, SHA256 `{old['package_sha256']}`. Re-uploading it would duplicate the published experiment.

Check its status in Windows Command Prompt:

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc status {old['entry_id']}
```

No numerical code/model/package from an earlier submission was changed. New root scripts and small evidence document the failed experiment; large model/count arrays and partial raw-source ranges are ignored by Git. No upload or Git push is performed.
''',encoding='utf-8');write(A/'release_decision.json',{'released':False,'uploaded':False,'new_package_generated':False,'reason':'All eligible geometry candidates fail frozen development release guards','retained_version':24,'retained_entry_id':old['entry_id'],'retained_overall':old['overall'],'retained_package_sha256':old['package_sha256'],'reserve_opened':False,'reserve_download_complete':False,'plan_sha256':b.sha256(A/'plan.json'),'selection_sha256':b.sha256(A/'selection.json'),'independent_audit_sha256':b.sha256(A/'independent_experiment_audit.json'),'poster_hashes_sha256':b.sha256(A/'poster_inputs.json'),'source_power_sha256':b.sha256(A/'source_power.json'),'guard_amendment_sha256':b.sha256(A/'release_guard_amendment.json'),'documenter_sha256':b.sha256(Path(__file__))});write(A/'runtime.json',{'python':sys.version,'executable':sys.executable,'versions':{p:importlib.metadata.version(p) for p in ['numpy','scipy','pandas','h5py','anndata','scikit-learn','threadpoolctl','requests']},'BLAS_threads':2,'required_runtime_guard':'v24_runtime.require_science','inherited_runtime_sha256':b.sha256(ROOT/'experiments/v24_cd4_response/analysis/runtime.json'),'plan_sha256':b.sha256(A/'plan.json')});print('V25 results documented; V24 retained, reserve unopened',flush=True)
if __name__=='__main__':main()
