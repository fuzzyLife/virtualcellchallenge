"""Render actual V59-V61 outcomes and the fixed experimental upload instructions."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from v57_pie import ROOT,read,write,sha

FOLDERS={59:ROOT/'experiments/v59_state_composition/analysis',60:ROOT/'experiments/v60_depth_adjusted_states/analysis',61:ROOT/'experiments/v61_pooled_state_response/analysis'}

def table(frame):
 cols=list(frame.columns);lines=['| '+' | '.join(map(str,cols))+' |','|'+'|'.join(['---']*len(cols))+'|']
 for row in frame.itertuples(index=False,name=None):
  values=[]
  for x in row:
   if pd.isna(x):values.append('missing')
   elif isinstance(x,(float,np.floating)):values.append(f'{x:.6f}')
   else:values.append(str(x))
  lines.append('| '+' | '.join(values)+' |')
 return '\n'.join(lines)

def main():
 A=FOLDERS[61];counts=read(A/'count_assessment.json');response=read(A/'response_assessment.json');run=read(A/'run_parameters.json');audit=read(A/'final_audit.json');assert counts['complete'] and audit['passed']
 rows=[]
 for version,folder in FOLDERS.items():
  r=read(folder/'response_assessment.json')
  for study,ratio in r['joint50_MSE_ratio_by_study'].items():rows.append(dict(version=version,study=study,primary_over_V24_MSE=ratio,beats_shuffled_and_pooled=r['beats_shuffled_and_pooled_each_study'][study],response_gate=r['passed']))
 history=pd.DataFrame(rows);history.to_csv(A/'version_comparison.csv',index=False)
 rf=pd.read_csv(A/'response_summary.csv');cf=pd.read_csv(A/'count_cases.csv');summary=cf.groupby(['source','arm'])[['fidelity','reach','jaccard','nmae','profile_MSE','predicted_calls','reference_calls']].mean().reset_index();oracle=pd.read_csv(A/'oracle_summary.csv');panels=pd.read_csv(A/'panel_diagnostics.csv');panel_summary=panels.groupby(['source','arm'])[['unique_targets','raw_panel_PDS','common_expression_numerator','common_expression_denominator']].mean().reset_index();panel_summary=panel_summary[panel_summary.arm.isin(['v24','joint50'])];null=pd.DataFrame([read(p) for p in sorted((A/'counts').glob('*_null.json'))])[['source','residual10_DE','real_control_DE']]
 fig,axs=plt.subplots(1,2,figsize=(13,5));p=history.pivot(index='study',columns='version',values='primary_over_V24_MSE');((p-1)*100).plot.bar(ax=axs[0],rot=0);axs[0].axhline(0,color='black',linewidth=1);axs[0].set_ylabel('Mean response MSE change vs V24 (%)');axs[0].set_title('Whole-study transfer: lower is better')
 o=oracle.set_index('source');gain=(1-o[['composition','within','joint','direct_mean']].div(o['zero'],axis=0))*100;gain.plot.bar(ax=axs[1],rot=25);axs[1].axhline(0,color='black',linewidth=1);axs[1].set_ylabel('Reconstruction improvement over zero (%)');axs[1].set_title('V61 independent population-view diagnostic');axs[1].legend(fontsize=8)
 fig.tight_layout();fig.savefig(A/'release_diagnostics.png',dpi=160);plt.close(fig)
 quality=pd.read_csv(FOLDERS[59]/'state_quality.csv').merge(pd.read_csv(FOLDERS[60]/'state_quality.csv'),on=['source','PC'],suffixes=('_v59','_v60'))
 partial=[]
 for path in sorted((FOLDERS[59]/'counts').glob('*.json')):
  r=read(path)
  if not path.stem.endswith('_null'):partial.extend(r['records'])
 pframe=pd.DataFrame(partial)
 review=f'''# State-response research: V59, V60 and V61

Three separately frozen versions tested the requested state-abundance/within-state decomposition. V24 remains the best published project submission (overall 0.142997855). These are repeatedly examined public datasets; source cells and generation seeds are not independent biological studies. D3 remains unopened.

## What changed and why

- **V59:** four shared statistical RNA states on 512 stable, nonpanel genes, experimental abundance ratios and conditional responses. It reuses exact V52 disjoint raw views. Before scoring, CPM closure was corrected; before transfer, absent source genes were masked for pairwise-valid basal Spearman routing. Original code and both correction receipts are preserved.
- **V60:** the NTC-only state audit found CD4 coordinates strongly associated with library depth/detection. A separately frozen variant regressed nonpanel depth/detection covariates using NTCs only. This can remove biological size/complexity too; it is not proof of technical confounding. The representation change did not pass matched-control response tests.
- **V61:** pool matched state sufficient statistics across lanes before nonlinear estimation/shrinkage. The old tiny-lane oracle diagnostic averaged unstable small-population errors and is superseded by population-level reconstruction; originals remain available. V61 retains the exact V60 state bases, source exclusions, bounds, primary/shuffle definitions and count tests. Empty individual lane-state controls borrow their lane's overall control mean before pooled support checks; this imputation is explicit.

All three deployable predictors exclude the recipient's complete study from response fitting and representation selection. The independent auditor makes recipient treatment fields unavailable and requires bitwise-identical predictions. Matching target evidence in other studies is legitimate cross-context retrieval, not unseen-target generalization.

## Whole-study response results

{table(history)}

V59 completed {len(pframe)} count comparisons on H1/HCT116 before stopping after its failed response gate; its remaining planned count cases were not run. V60 stopped after the failed response gate and has no count result claimed. V61 completed all {counts['comparisons']} planned count comparisons to characterize the fixed release candidate, including the separate state-cell diagnostic. No partial pilot is represented as a completed 264-comparison experiment.

## NTC state-quality check

{table(quality[['source','PC','library_spearman_v59','library_spearman_v60','detected_genes_spearman_v59','detected_genes_spearman_v60']])}

The state figures use lane-centered control coordinates; visual overlap cannot establish biological equivalence. The population-oracle correction is an openly post-hoc diagnostic and never enters fitted release features or candidate selection. See `evaluate_state_population_diagnostic.py --version 60` for the corrected comparison on the original lane-wise estimator.

[Detailed V61 results](V61_RESULTS.md) · [Fixed release and upload instructions](V61_SUBMISSION.md) · [V59 protocol](V59_STATE_COMPOSITION_STRATEGY.md) · [V60 protocol](V60_DEPTH_ADJUSTED_STRATEGY.md) · [V61 protocol](V61_POOLED_STATE_STRATEGY.md).
'''
 (ROOT/'V59_V61_STATE_REVIEW.md').write_text(review,encoding='utf-8')
 verdict='Public tests support the fixed primary under the recorded gates; this is still not a private-score guarantee.' if counts['release_supported'] else 'The fixed primary failed the recorded public support gates. This release is experimental and is not recommended as a validated replacement for V24.'
 result=f'''# V61 pooled state-response results

{verdict}

## Response and population reconstruction

All {response['cases']} eligible target-context cases were evaluated at five arms, using complete recipient-study exclusion. State genes and depth/detection covariates exclude all 300 panel targets. Raw matched views and immutable source checksums are retained. D3 remains unopened.

{table(rf[['study','source','arm','MSE','MAE']])}

Response gate: **{response['passed']}**. Equal-study primary/V24 MSE ratio: **{response['equal_study_mean_ratio']:.6f}**. A lower number beats V24. Primary joint50 must additionally beat matched shuffled and pooled controls; these requirements are not relaxed after outcomes.

The corrected oracle estimates a mechanism on one disjoint matched population and reconstructs the other before computing error. It uses treatment labels solely as a diagnostic. Direct mean and zero response are controls, not uploadable oracles. Lower MSE is better.

{table(oracle)}

Within-dataset reconstruction can improve while cross-study prediction worsens. That difference limits the response information transferable to anonymous competition contexts. It does not prove that these models identify cell survival, proliferation or observed individual counterfactuals.

## All integer-count comparisons

All {counts['comparisons']} comparisons completed. Each target-seed population is common across arms; 400 generated cells and up to 400 distinct measured treatment cells are compared against independent real controls. HCT116/HEK293T are one Orion study and both CD4 states are one study. Seeds repeat biology. The five mean arms use identical templates and the unchanged residual10 emitter. The state_cells50 diagnostic changes both means and cell distributions; missing mechanism coordinates in that diagnostic can dilute the V24 component, so it is not a clean variance-only test or an eligible finalist. No diagnostic was promoted after outcomes.

{table(summary)}

Fidelity/reach/Jaccard are raw higher-is-better metrics. NMAE and profile MSE are errors; missing NMAE stays missing. Profile MSE excludes the directly targeted gene; nonpanel MSE and corrected error numerators are saved separately. These are not official anchor-scaled scores. No private overall/rank forecast is made.

Recorded count/release checks:

```json
{json.dumps(counts['checks'],indent=2)}
```

## Limited common-control discrimination panels

{table(panel_summary)}

These panels contain only four or eight unique targets per source/seed. All 300 panel genes are excluded from PDS coordinates; one common NTC profile is verified bitwise. Constant predictions and the balanced cyclic label-shift control both score exactly .5. These are limited diagnostics, not official 300-target PDS or a private overall forecast. The HCT116 expression denominator is unresolved/negative; signed numerator and denominator are retained separately and unstable ratios are not averaged into a release score.

## Zero-intervention comparison

{table(null)}

The state-cell zero-intervention contract returns real unedited NTC rows. Null calibration does not establish treatment prediction. Own-target knockdown is applied once by the existing count emitter and excluded from direction scoring.

![Measured release diagnostics](experiments/v61_pooled_state_response/analysis/release_diagnostics.png)

## Production scope

The release keeps the unchanged residual10 count emitter. Of 900 target/context groups, {run['unchanged_v24_groups_reused']} reuse validated V24 counts exactly; the remaining {900-run['unchanged_v24_groups_reused']} were regenerated from the fixed V61 effects. The unsuccessful state_cells50 generator is excluded. All 360,000 cells use original guide-stratified control templates with exact integer libraries.

## Integrity and release scope

The production kernel is independently implemented and compared bitwise against all public primary predictions. Every production effect is additionally reconstructed through the original transfer implementation. Counts retain exact original library sizes and unique guide-stratified templates. Full package verification compares every stored value/index/pointer and all axes, verifies provenance and replays 12 complete groups.

The user's requested fixed experimental candidate is V61 joint50, frozen before its count outcomes. `--experimental` is explicit in its release command if gates fail. It preserves V24, reports failure and makes no claim of superior private performance. No leaderboard upload or Git push is performed.

[Protocol and commands](V61_POOLED_STATE_STRATEGY.md) · [Upload artifact](V61_SUBMISSION.md) · [Three-version review](V59_V61_STATE_REVIEW.md).

Scientific plan SHA256: `{sha(A/'plan.json')}`. Production plan SHA256: `{sha(A/'production_plan.json')}`.
'''
 (ROOT/'V61_RESULTS.md').write_text(result,encoding='utf-8')
 readyfile=A/'submission_ready.json';validationfile=A/'independent_validation.json'
 if readyfile.exists() and validationfile.exists():
  ready=read(readyfile);valid=read(validationfile);assert valid['passed'] and ready['sha256']==valid['package_sha256'];filename=Path(ready['file']);assert sha(filename)==ready['sha256']
  submission=f'''# V61 experimental submission

{verdict}

Ready file: `experiments/v61_pooled_state_response/prediction.vcc`.

- Size: {ready['bytes']:,} bytes.
- SHA256: `{ready['sha256']}`.
- Shape: 360,000 cells × 18,533 genes; 300 targets × 400 cells × A/B/C.
- Official VCC prep: passed.
- Complete packaged value/index/pointer and axis comparison: passed ({valid['packaged_nonzero_values_checked']:,} nonzeros).
- Original control provenance/library checks: all 360,000 cells.
- Exact count replay: {valid['exact_replayed_groups']} full groups.
- Private leaderboard result: not submitted; unknown.

Windows CMD:

```bat
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v61_pooled_state_response\\prediction.vcc -m "Pooled state-response transfer v61 experimental" --wait
```

This uploads an experimental diagnostic candidate. V24 remains available at `experiments/v24_cd4_response/prediction.vcc`; the new release does not overwrite it. [Measured tests](V61_RESULTS.md) and [reproduction](V61_POOLED_STATE_STRATEGY.md).
'''
 else:submission=f'# V61 release status\n\n{verdict}\n\nThe production artifact has not yet passed its final package audit. Do not upload an incomplete package.\n'
 (ROOT/'V61_SUBMISSION.md').write_text(submission,encoding='utf-8')
 readme=ROOT/'README.md';old=readme.read_text(encoding='utf-8-sig');begin='<!-- v61-state-release:begin -->';end='<!-- v61-state-release:end -->'
 if begin in old:old=old[:old.index(begin)]+old[old.index(end)+len(end):].lstrip('\n')
 block=f'''{begin}
## V61: pooled state-response experimental release

Three state-composition variants tested; V61 completed {response['cases']} public response cases and {counts['comparisons']} paired count comparisons. Public release support: **{counts['release_supported']}**. {verdict} D3 remains unopened.

[Upload file and command](V61_SUBMISSION.md) · [Measured V61 results](V61_RESULTS.md) · [V59–V61 review](V59_V61_STATE_REVIEW.md) · [Reproduce](V61_POOLED_STATE_STRATEGY.md). Root scripts retain the original methods, pre-score corrections, negative results and explicit experimental candidate lock. No automatic upload.
{end}
'''
 readme.write_text(block+'\n'+old,encoding='utf-8');write(A/'documentation.json',dict(reports={name:sha(ROOT/name) for name in ['V59_V61_STATE_REVIEW.md','V61_RESULTS.md','V61_SUBMISSION.md']},plot_sha256=sha(A/'release_diagnostics.png'),writer_sha256=sha(Path(__file__)),release_supported=counts['release_supported']))
 print(verdict,flush=True)
if __name__=='__main__':main()
