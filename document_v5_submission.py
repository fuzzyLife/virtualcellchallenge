"""Document the verified V5 submission, method, evidence, source provenance and upload command."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import build_submission as build
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'experiments/v5_atlas_calibrated';A=OUT/'analysis'


def table(frame):
    return '\n'.join(['| '+' | '.join(frame.columns)+' |','| '+' | '.join(['---']*len(frame.columns))+' |']+
        ['| '+' | '.join(str(v).replace('|','\\|') for v in row)+' |' for row in frame.itertuples(index=False,name=None)])


def main():
    ready=json.loads((A/'submission_ready.json').read_text());model=json.loads((A/'reference_model.json').read_text())
    validation=json.loads((A/'independent_validation.json').read_text());model_validation=json.loads((A/'model_validation.json').read_text())
    preflight=json.loads((A/'preflight.json').read_text());promoter=json.loads((A/'promoter_selection.json').read_text())
    bulk=json.loads((A/'bulk_selection.json').read_text());h1=json.loads((A/'h1_lineage_selection.json').read_text())
    assert ready['official_prep_passed'] and validation['passed'] and model_validation['passed'] and preflight['passed']
    assert validation['source_model_sha256']==build.sha256(A/'reference_model.json')
    assert model_validation['source_model_sha256']==build.sha256(A/'reference_model.json')
    assert model_validation['script_sha256']==validation['script_sha256']==build.sha256(ROOT/'validate_v5.py')
    assert preflight['source_model_sha256']==build.sha256(A/'reference_model.json')
    assert preflight['script_sha256']==build.sha256(ROOT/'build_v5_submission.py')
    assert validation['generated_package_checks']['package_sha256']==ready['sha256']
    assert ready['bytes']==(OUT/'prediction.vcc').stat().st_size
    assert 'Ran 5 tests' in (A/'test_results.txt').read_text(encoding='utf-8-sig') and '\nOK' in (A/'test_results.txt').read_text(encoding='utf-8-sig')
    qc=pd.read_csv(A/'h1_target_QC.csv');stats=pd.read_csv(A/'prediction_group_statistics.csv');changes=pd.read_csv(A/'promoter_effect_changes.csv')
    extraction=json.loads((ROOT/'reference/v5_h1/extraction.json').read_text());build.write_json(A/'source_extraction.json',extraction)
    pd.read_csv(ROOT/'reference/v5_h1/target_coverage.csv').to_csv(A/'h1_source_coverage.csv',index=False)
    old_targets=set(pd.read_csv(ROOT/'reference/h1_response/selected_cells.csv').target_gene)-{'non-targeting'}
    new=qc[~qc.target_gene.isin(old_targets)].copy();new['direct_response_accepted']=new.direct_response_accepted.map({True:'yes',False:'no'})
    new=new[['target_gene','selected_cells','own_log2fc','direct_response_accepted']];new.own_log2fc=new.own_log2fc.map(lambda v:f'{v:.3f}')
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    axes[0,0].bar(['V4 available','V5 available'],[13,model['usable_direct_H1_targets']],color=['#6b7280','#167d8d'])
    axes[0,0].set(ylabel='QC-accepted direct H1 targets',title='Direct public-source response coverage')
    axes[0,1].bar(['Existing donor/K562','V5 conditional prior'],[promoter['v4_donor_or_k562_baseline_loss'],promoter['held_target_loss']],color=['#6b7280','#167d8d'])
    axes[0,1].set(ylabel='Held-target neighboring-gene log2FC MSE',title='CD4 donor test: repeated source suppression required')
    for c,color in zip('ABC',['#167d8d','#b47726','#815ca1']):
        group=stats[stats.context==c]
        axes[1,0].hist(group.requested_mean_probability_L1,bins=30,alpha=.5,label=c,color=color)
        axes[1,1].hist(group.projected_pooled_probability_L1,bins=30,alpha=.5,label=c,color=color)
    axes[1,0].set(xlabel='L1 to requested mean probability',ylabel='Target groups',title='All 900 generated groups: mean-moment residual')
    axes[1,1].set(xlabel='L1 to feasible pooled probability',ylabel='Target groups',title='All 900 generated groups: pooled-moment residual')
    axes[1,0].legend(title='Context');axes[1,1].legend(title='Context')
    fig.suptitle('V5 validation diagnostics: no hidden leaderboard score',fontsize=14)
    fig.savefig(A/'v5_validation.png',dpi=180);fig.savefig(A/'v5_validation.pdf');plt.close(fig)
    command='''cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v5_atlas_calibrated\\prediction.vcc -m "Atlas-calibrated lineage and dual moments v5" --wait'''
    result=f'''# V5: verified submission ready for upload

Created for the 2026 validation panel `vcc2026-val-1`. The candidate is **generated, officially packaged and independently verified**, and has not been uploaded or scored. V4 remains published with overall -0.05067391016092112; local V5 diagnostics do not establish a leaderboard improvement.

| Check | Result |
| --- | --- |
| Submission file | `experiments/v5_atlas_calibrated/prediction.vcc` |
| Size | {ready['bytes']:,} bytes ({ready['bytes']/1024**3:.2f} GiB) |
| Shape | 360,000 cells by 18,533 genes |
| Groups | All 900 context/target combinations, 400 cells each |
| Gene order | Exact official axis |
| Counts | Every stored value is a positive integer; zeros remain sparse |
| Library sizes | All 360,000 match the selected original control-cell libraries exactly |
| Sampling | Each target has 400 unique real controls with all 46 guides balanced |
| Public H1 labels | {model['public_H1_targets']} overlapping targets; {model['usable_direct_H1_targets']} accepted for direct responses |
| Promoter corrections | {len(changes)} A-only target/neighbor effects with repeated source suppression |
| Numerical preflight | {preflight['groups']} groups passed |
| Focused tests | Five passed |
| Official packaging | Passed; VCC CLI {ready['package_meta']['cli_version']} |
| Independent source verification | {validation['original_public_H1_rows_independently_checked']} original source rows checked; all 1,536 reused selected controls exactly matched in each extra split |
| Independent source moments | {validation['source_contrast_values_recalculated']:,} values recalculated |
| Independent final effects | {validation['all_mean_and_bulk_effects_reassembled']:,} mean/pooled effects reassembled |
| Deterministic full-count replay | {validation['generated_package_checks']['exact_generated_groups_replayed']} complete groups checked |

Package SHA-256:

```text
{ready['sha256']}
```

## Submit using Windows Command Prompt

Reuse the existing login:

```cmd
{command}
```

The command uploads a new competition entry and waits for scoring. It has not been run during generation. Preserve the resulting entry ID, panel/anchor set and all six scaled/raw metrics to compare with V4. No token needs to be included in the command or committed.

## What changed

1. **Public H1 coverage:** Add released 2025 validation/test perturbation observations to the original training subset. The repeated control identifiers are checked against original counts and kept only once. Global batch partitions provide disjoint control IDs for calibration. Direct estimates require at least 40 selected target cells, at least 30% mean-expression knockdown, source control mean above 5 CPM and a negative own-target effect larger than 1.96 estimated SE. Every response gene also needs measured support and source control expression above 5 CPM.
2. **Separate expression moments:** Fit H1 mean-CPM log2FC and library-weighted pooled `log1p(50000*p)` changes independently. Pooled uncertainty retains count/depth covariance. Source batches and held targets select amplitudes. CD4 centering remains zero and its V4 response estimator is retained.
3. **Preserved cell variation:** Start each 400-cell group from individual guide/depth-balanced original controls. Fit requested mean and pooled profiles to these individual templates, allowing only the documented feasible-moment projection, then round to raw integer counts while preserving each exact original library size. No four-cell averaging is used.
4. **Conditional promoter correction:** Apply the distance-dependent prior only to A, only when a neighbor has significant negative responses in at least two qualified CD4 conditions, and preserve existing stronger negative effects. The held-target/disjoint-donor test uses {promoter['targets']} targets/{promoter['observations']} observations. Relative neighboring-gene error is {promoter['relative_v4_neighbor_loss']:.3f} of the existing donor/K562 estimator. These are selected source observations, not hidden A labels.

The 12 additional source targets and their QC outcome are:

{table(new)}

Full input sampling, uncertainties, selected coefficients, transfer rules and reproduction are documented in [V5_MODEL_STRATEGY.md](V5_MODEL_STRATEGY.md). Root scripts and all compact validation records are version-controlled. Large counts/effect caches and H5AD/VCC files remain local inputs/outputs.

![Generation and source-validation diagnostics](experiments/v5_atlas_calibrated/analysis/v5_validation.png)
'''
    (ROOT/'V5_SUBMISSION.md').write_text(result,encoding='utf-8')
    strategy=f'''# V5 model strategy and reproduction

## Source inheritance and isolation

V5 is written to `experiments/v5_atlas_calibrated/`. Its parent is the frozen V4 model and source-component arrays, whose hashes are checked. Original competition control hashes match V4. Earlier packages are preserved. The upstream numerical solver is imported from the reviewed AtlasShift clone after its SHA-256 is checked against the pinned review manifest. The clone retains its MIT license. There is no dependency upgrade or execution of the upstream full-atlas preparation pipeline.

The reviewed repositories, commits, initial control reconstruction experiment and limitations are documented in [UPSTREAM_REPOSITORY_REVIEW.md](UPSTREAM_REPOSITORY_REVIEW.md). This candidate implements the supported ideas with explicit source-quality filters and checks.

## Bounded public H1 extraction

`extract_v5_h1.py` verifies the existing 2025 training subset, reads its original cell identifiers, and extracts target/control rows from generation-pinned public validation/test URLs. Its default limit is 4 GiB of used HTTP ranges per added split, including cached blocks. It chooses the largest predefined cells-per-target/batch limit fitting the planned ranges, from 16, 8, 4, 2 or 1. Rows are selected in source order; this is a bounded sampling approximation, not the complete atlas.

All 1,536 selected training controls are matched by identifier in both extra splits and compared across the entire gene axis. Their count equality is asserted before keeping one copy in the combined training set. Combined cell IDs must be unique. Selected response observations have finite, nonnegative integer raw counts and positive total counts. The measured H1 gene axis is retained separately; unmeasured official genes are supplied by inherited response sources, not imputed source measurements.

`source_extraction.json` and `h1_source_coverage.csv` retain source URL generations/ETags, counts, byte limits and original extraction hashes. Cached HTTP ranges have SHA sidecars. Source arrays are in ignored `reference/v5_h1/`.

## Source contrasts and target/batch validation

`fit_v5_submission.py` uses the same globally sorted batch partition for every public split. Left/right calibration controls are therefore disjoint by cell identifier. It forms both partition-specific and full-source target contrasts using controls matched to represented batches.

The mean-CPM contrast is `log2((target_mean_CPM + 1)/(matched_control_mean_CPM + 1))`; its SE includes target and control mean variances. Pooled probabilities are sums of raw gene counts divided by summed libraries. For their ratio uncertainty, retain covariance through the residual `count_gene - probability_gene * cell_library`; compute variance of this residual, divide by sample size and squared mean library. The pooled contrast is `log1p(50000 * target_probability) - log1p(50000 * matched_control_probability)`. Control batch contributions are weighted by target library mass. Both contrast spaces are clipped to +/-3. Source truth significant sets use actual two-sided Mann-Whitney/Wilcoxon tests of cell CPM with within-target BH; selection calls use the stated 400-cell normal approximation.

The reused V4 H1 selector has three outer target folds/two inner folds and selects specific/prior coefficients while guarding zero-relative response MSE at 1.02. Its selected mean-CPM coefficients are `{json.dumps(model['h1_settings'],sort_keys=True)}`. Pooled gain candidates are 0.3, 0.5, 0.7 and 1.0, with five outer target folds. Gain selection must not exceed the training error of the conservative 0.7 reference. Adoption of a changed gain additionally requires no nested error increase and no composite-proxy decline. The selected pooled gain is **{model['bulk_gain']}**. This conservative reference is not a historical separate-bulk V4 model.

The mean-response choice has a measurable tradeoff. Against the conservative gain-0.7/no-prior comparison on these same public-source folds, nested response MSE is {h1['nested_procedure']['loss']:.6f} versus {h1['v3_fixed_baseline']['loss']:.6f}; direction fidelity proxy is {h1['nested_procedure']['direction_fidelity_proxy']:.3f} versus {h1['v3_fixed_baseline']['direction_fidelity_proxy']:.3f}, and Jaccard proxy is {h1['nested_procedure']['jaccard_proxy']:.3f} versus {h1['v3_fixed_baseline']['jaccard_proxy']:.3f}. Selection favors the combined direction/overlap proxy while retaining the zero-relative MSE guard. This is a source-domain comparison against a conservative estimator, not a hidden V4/V5 score comparison.

Direct-response QC and the other-target shared H1 prior use only qualified targets. Each target is excluded from its own shared prior. The H1 lineage selector operates on source contrasts before the final strict direct-target QC; its diagnostics are source-domain proxies and do not prove anonymous-context accuracy. Complete folds, IDs and candidate scores are saved.

## Context policies and promoter calibration

A keeps V4's matched-lineage CD4 response with unsupported-coordinate K562 fallback. Its pooled effect is the corresponding change in the context's pooled log profile. CD4 gain/centering settings remain `{json.dumps(model['cd4_settings_unchanged'],sort_keys=True)}`; no across-target subtraction is added.

B/C replace supported response coordinates with measured H1 mean and pooled effects. Other coordinates retain `{json.dumps(model['cross_lineage_settings'],sort_keys=True)}` CD4/K562 transfer coefficients and an other-target H1 shared prior. Pooled transfer from CD4/K562 uses the H1-calibrated gain as an explicit cross-lineage assumption. Source effects are moderated by the unchanged gene and target basal-expression gates `min(1, context_mean_CPM/5)`. Own-target downstream effects are excluded and handled by the separately stated 90% knockdown assumption.

`calibrate_v5_promoters.py` reopens the two disjoint CD4 donor subsets, permits neighboring-KD flags for this narrow calibration, and keeps the remaining on-target/multi-guide/minimum-cell/no-distal/no-low-expression QC. Five outer and three inner target folds select strength from 0, 0.25, 0.5 or 1. The baseline is each source donor's V4-style soft CD4 estimator when eligible, with independent K562 fallback otherwise. The donor-specific prior uses only other training targets. It is applied only if that source donor has significant negative effects in at least two represented conditions; held-out truth-donor results are not used to set a prediction's support gate.

Adoption requires at least ten observed targets, at least 10% nested MSE improvement and at least 70% both-donor negative signs. The accepted fitted close-neighbor log2FC is **{promoter['fitted_close_neighbor_log2fc']:.6f}** (remaining fraction approximately {2**promoter['fitted_close_neighbor_log2fc']:.3f} before profile normalization). A logarithmic distance ramp is constant through 500 bp and reaches zero correction at 5 kb. Application in A requires repeated significant suppression in full qualified CD4 source conditions and retains stronger existing negative effects. It changed {len(changes)} coordinates. Gene-boundary TSS positions are approximate promoter positions, not guide locations; this is not a universal suppression rule or a B/C correction.

## Integer count generation

`build_v5_submission.py` uses seed 20261005. Each context/target has a separate deterministic seed, 400 unique original control rows and balanced source guides/depth quantiles. Templates are individual cell-normalized profiles, with no pooling. Desired mean probability is proportional to the context mean probability times the predicted fold response. Desired pooled probability is obtained by shifting the context pooled `log1p(50000*p)` profile. The own target is multiplied by 0.1 in both profiles before renormalization.

A depth-dependent feasible moment projection bounds the pooled profile while leaving the requested mean fixed. Projection L1 must be at most 0.03. The reviewed iterative solver is capped at 200 iterations; its continuous fitting tolerance/acceptance checks remain active. Integer emission uses seeded systematic rounding and column repair while retaining every original cell depth. The resulting mean probability L1 must be at most 0.03 and pooled probability L1 against the feasible projection at most 0.002. All 900 groups passed. Maximum integer mean residual: {validation['generated_package_checks']['maximum_mean_probability_L1']:.6f}; maximum feasible pooled residual: {validation['generated_package_checks']['maximum_projected_pooled_probability_L1']:.6f}.

Generation uses a bounded window of four parallel groups by default, with two BLAS threads per fit, and writes groups in official panel order. `--workers 1` provides serial generation with the same seeds/algorithm. Tests verify parallel deterministic output. CSV records retain every original control ID, row, guide and depth plus all group residuals/projections/detection/variance diagnostics. The `CPM_total_variance` diagnostic column stores the sum of variances of cell expression probabilities (counts divided by library), in probability-squared units; multiply by 10^12 to express it as total CPM variance.

## Reproduction from Command Prompt

Prerequisites: `.competition_env`, original competition inputs, unchanged V4 source arrays/model, original CD4 donor/merged caches and H1 training subset, and the reviewed upstream clone. The exact source/clone hashes are in the parent review and V5 manifests. No model evaluation uses hidden 2026 responses or scoreboard tuning.

```cmd
.\\.competition_env\\Scripts\\python.exe extract_v5_h1.py --budget-mib 4096
.\\.competition_env\\Scripts\\python.exe calibrate_v5_promoters.py
.\\.competition_env\\Scripts\\python.exe fit_v5_submission.py
.\\.competition_env\\Scripts\\python.exe -m unittest test_v5_submission -v > experiments\\v5_atlas_calibrated\\analysis\\test_results.txt 2>&1
.\\.competition_env\\Scripts\\python.exe validate_v5.py --model-only
.\\.competition_env\\Scripts\\python.exe build_v5_submission.py --preflight --workers 4
.\\.competition_env\\Scripts\\python.exe build_v5_submission.py --workers 4
.\\.competition_env\\Scripts\\python.exe validate_v5.py
.\\.competition_env\\Scripts\\python.exe document_v5_submission.py
```

The build refuses to overwrite existing final predictions unless `--force` is explicitly supplied. `--package-only` packages an existing H5AD; it does not submit. Cached fitted contrasts require matching script/input hashes. Extracted source reuse requires checksums and pinned object identity; `extractor_code_versions.json` preserves the original and accelerated extractor producers, which use identical row-selection rules and byte limits. If the fitter is intentionally changed, preserve or regenerate its old contrast cache with its provenance rather than bypassing checks. The official VCC 0.2.2 preparation path validates raw counts and packages the file; the private upstream VCC 0.2.0 adapter is not used.

## Verification and limits

Five focused tests cover library-ratio covariance, distinct profile/knockdown semantics, held-target prior exclusion, deterministic integer/parallel generation and infeasible constraints. Independent validation spot-checks newly extracted rows against checksum-verified offline original HDF5 ranges, checks all duplicated selected control counts, recalculates every source mean/pooled effect and uncertainty, and reassembles all {validation['all_mean_and_bulk_effects_reassembled']:,} final mean/pooled effects. It also checks target/batch separation, every generated raw count, all 360,000 exact cell depths and original control records, the completed package hash and twelve complete deterministic groups.

All generated group diagnostics and official preparation output are retained. Count distributions are tied to supplied controls: perturbation-dependent cell-depth changes are not modeled. H1 source sampling, selected target panels and inferred lineage matches limit generalization. Neutral profile reconstruction and source-target validation cannot establish the hidden response score; compare the eventual published V5 result with V4 on the same panel/anchors. [V5_SUBMISSION.md](V5_SUBMISSION.md) records the verified artifact and upload command.
'''
    (ROOT/'V5_MODEL_STRATEGY.md').write_text(strategy,encoding='utf-8')
    section=f'''## V5 atlas-calibrated submission

A new candidate is generated and verified at `experiments/v5_atlas_calibrated/prediction.vcc` ({ready['bytes']:,} bytes; {ready['bytes']/1024**3:.2f} GiB). It covers all 900 groups, 360,000 cells and the exact 18,533 genes. Official VCC packaging and independent source/model/count verification passed. It has not been uploaded or scored.

V5 expands public H1 coverage to {model['public_H1_targets']} source targets ({model['usable_direct_H1_targets']} QC-accepted direct responses), estimates mean/pooled effects separately, generates from individual control-cell templates while preserving all original cell depths, and applies {len(changes)} source-supported promoter-neighbor corrections only in A. CD4 across-target subtraction remains zero. Five tests and 24 numerical preflight groups passed. [V5_SUBMISSION.md](V5_SUBMISSION.md) gives the upload command and hash; [V5_MODEL_STRATEGY.md](V5_MODEL_STRATEGY.md) gives every step, selected settings, verification and limitations.

Root scripts are `extract_v5_h1.py`, `calibrate_v5_promoters.py`, `fit_v5_submission.py`, `build_v5_submission.py`, `test_v5_submission.py`, `validate_v5.py` and `document_v5_submission.py`. Small evidence is tracked in the candidate analysis folder; source counts/caches and large predictions are regenerated locally. V4 remains the best published project result until V5 is scored.

'''
    path=ROOT/'README.md';readme=path.read_text(encoding='utf-8');heading='## V5 atlas-calibrated submission\n';anchor='## Upstream repository review and validated improvement experiments\n'
    if heading in readme:
        begin=readme.index(heading);end=readme.index(anchor,begin);readme=readme[:begin]+readme[end:]
    assert anchor in readme;readme=readme.replace(anchor,section+anchor,1);path.write_text(readme,encoding='utf-8')
    code=['extract_v5_h1.py','calibrate_v5_promoters.py','fit_v5_submission.py','build_v5_submission.py','test_v5_submission.py','validate_v5.py','document_v5_submission.py']
    artifacts=['reference_model.json','h1_contrasts.json','h1_target_QC.csv','h1_lineage_selection.json','bulk_selection.json','promoter_selection.json',
               'promoter_effect_changes.csv','source_extraction.json','extractor_code_versions.json','preflight.json','run_parameters.json','prediction_group_statistics.csv',
               'prediction_control_rows.csv','model_validation.json','independent_validation.json','submission_ready.json','test_results.txt','official_prep_stdout.json','v5_validation.png','v5_validation.pdf']
    build.write_json(A/'documentation_manifest.json',{'root_code_sha256':{p:build.sha256(ROOT/p) for p in code},
        'artifact_sha256':{p:build.sha256(A/p) for p in artifacts},
        'document_sha256':{p:build.sha256(ROOT/p) for p in ['V5_SUBMISSION.md','V5_MODEL_STRATEGY.md']},'scope':'Completed generation and verification records; no upload or hidden-score claim.'})
    print('Documented verified V5 package, methods, README and artifact provenance.',flush=True)


if __name__=='__main__':main()
