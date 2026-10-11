"""Render V14's measured outcomes and provenance without modifying frozen experiments."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import build_submission as build
ROOT=build.ROOT
A=ROOT/'experiments/v14_context_interactions/analysis'
FIELDS=['pds','fidelity','reach','jaccard','nmae','common_expression_ratio']


def table(frame):
    def value(x):
        if isinstance(x,(float,np.floating)):return 'NA' if not np.isfinite(x) else f'{x:.6g}'
        return str(x)
    return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |']+['| '+' | '.join(value(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)])


def main():
    selection=json.loads((A/'selection.json').read_text());audit=json.loads((A/'independent_audit.json').read_text())
    assert audit['passed'] and not selection['promote'],'This report describes the completed rejected V14 experiment'
    for stage in ['development','oracle']:assert json.loads((A/f'{stage}_evaluation.json').read_text())['completed']
    assert selection['development_sha256']==build.sha256(A/'development_summary.csv')
    summary=pd.read_csv(A/'development_summary.csv');panels=pd.read_csv(A/'development_panels.csv')
    primary=panels[panels.reference_mode=='full_independent']
    family=primary.groupby(['family','candidate'])[FIELDS].mean().reset_index()
    key=family[family.candidate.isin(['v13',selection['candidate']])]
    oracle=pd.read_csv(A/'oracle_panels.csv')
    oracle_table=oracle[['family','candidate','nmae','jaccard','predicted_calls','reference_calls']]
    denominator=key[['family','candidate','common_expression_ratio']].merge(
        primary.groupby(['family','candidate'])[['common_expression_numerator','common_expression_denominator']].mean().reset_index(),on=['family','candidate'])
    power=primary[primary.candidate=='v13'][['family','seed','reference_calls','nmae_eligible','nmae']]
    screen=pd.read_csv(A/'screen_summary.csv')
    failures=[k for k,v in selection['checks'].items() if not v]
    deltas=pd.DataFrame([{'metric':k,'context2_minus_v13':v} for k,v in selection['changes'].items()])
    text=f'''# V14 experiment results (2026-10-11)

**Decision: reject the context-interaction candidate; retain V13.** This research experiment does not create or upload a prediction package. V13 remains the best recorded project submission, overall 0.137883134. Public raw scores below are not private anchor-scaled leaderboard scores.

## What was tested

V14 retains the exact held-family V13 pooled response and adds a bounded, target/gene-specific interaction with two control-expression rank components. The context feature basis and all response slopes exclude the entire evaluated family. CD4's conditions remain one family. No panel target genes enter the feature basis, missing measurements abstain, and own-target effects remain excluded. Penalties 2 and 8 were screened; penalty 2 advanced. Fixed shuffled context labels and simple 1.1x amplification serve as controls. See [V14_STRATEGY.md](V14_STRATEGY.md) and root `v14_context_methods.py`.

Actual-count development compares four models across nine targets in H1, HCT116 and HEK293T, two paired seeds, 400 distinct real treatments and 400 predictions per case. Full independent controls are primary; sampled400 references are secondary. That is **72 synthetic populations**, 144 scored population/reference combinations, and 48 profile-panel rows. Seeds are repeated draws of the same nine biological contrasts, not independent lineages. Models use public training families only; held-family treatment data enter evaluation, never fitting.

Separately, nine oracle cases split 400 distinct treatment cells into 200 calibration and 200 truth cells. Three arms produce **27 population cases: 18 synthetic populations and nine independent real-treatment replicas**. Combined with development, there are 99 population cases, including 90 generated populations. Oracle calibration/truth splits are disjoint within each case, but reuse development targets and can overlap development identities across experiments; they are not independent confirmation experiments.

## Response screen

The best interaction proxy barely exceeds V13. Both shuffled labels and amplification score higher, which does not support interpreting the proxy gain as learned context biology. K562's held-out feature coverage is only 41.5%, below the predeclared 70% threshold; it abstains and cannot demonstrate an interaction benefit.

{table(screen)}

These are previously inspected public summary contrasts, not a blind test or official count metrics. Zero/generic controls can have misleadingly high set-overlap proxies for weak/empty call sets.

## Completed 400-cell count test

Higher is better for PDS, fidelity, reach and Jaccard; lower is better for NMAE and expression ratio. Scores first average within source family, then across families with available values. **NMAE is unavailable for all HCT116 primary panels**, so its displayed aggregate covers H1/HEK293T only. Missing primary metrics block promotion rather than being filled with zero.

{table(summary)}

{table(deltas)}

Only reach improves. The candidate fails these frozen checks: **{', '.join(failures)}**. Primary discrimination drops 0.111111; fidelity drops 0.005150; Jaccard drops 0.000158; NMAE increases 0.001906. The expression ratio worsens beyond the prospectively declared exploratory 1% tolerance. The unanchored utility also loses to V13 and shuffled labels. This utility is dominated by HCT116's tiny expression denominator and must not be read as an estimated competition score. Independent failures in discrimination, metric breadth and missing-metric completeness support retaining V13 without relying on that utility.

### Family results

{table(key)}

### Metric power and denominator diagnosis

The baseline's mean real DE-call counts, eligible fraction of the three target cases, and NMAE for each primary panel are:

{table(power)}

Expression numerator and denominator are retained separately. HCT116 has an extremely small control-distance denominator; enormous ratios are unstable and must not be used as a headline biological error ranking. Family ratios are calculated on common-origin panels rather than averaged per-target ratios.

{table(denominator)}

No reserved V14 algorithm-verification outcomes were opened after rejection. Earlier history already contains some of those biological contrasts, so they cannot be described as blind biological confirmation in any future run.

## Treatment-mean versus population diagnostic

The oracle uses **known treated-calibration depths in every arm**, including the V13 comparator. This deliberately isolates response/population limitations; it is not deployable baseline inference. Measured treatment mean-CPM and pooled proportions constrain compact emission; actual replicas use the original 200 calibration cells. All are evaluated against the disjoint 200 truth cells using the same independent controls and scorer.

{table(oracle_table)}

For H1, supplying measured treatment moments reduces NMAE from **0.8910 to 0.2401**, close to the independent real replica's **0.2363**. But compact oracle cells call roughly **1,643 DE genes** against **61** mean truth calls, with Jaccard **0.0260**. Independent real treatment replicas call **60.3**, with Jaccard **0.2963**. The means help response accuracy while the synthetic distribution still fails to reproduce the DE set.

HEK293T shows the same distribution warning: oracle NMAE improves from **0.5698 to 0.1669** on its eligible contrast, but its mean DE calls increase to **1,231** versus **12.3** in truth. Jaccard falls from **0.0440 to 0.00760**, compared with **0.6695** for real replicas. Two of its three truth cases do not meet NMAE's ten-DE-gene minimum. HCT116 averages only one truth DE gene at 200 cells and has no eligible NMAE case; its ratios and discrete overlap estimates have little power.

This supports **two limitations**, not a proof that distribution alone explains our leaderboard plateau. Transferred means remain imperfect. Compact populations also behave very differently from independent real replicas even when given measured moments. DE calls beyond the truth set are not all proven false positives: the 200-cell truth has limited power, moments themselves are estimated from 200 cells, and residual batch composition can matter. Fitting moments does not fix zero frequencies, heterogeneous conditional means or tails. Expression error includes jackknife terms, so measured moments do not guarantee a perfect expression ratio.

## What this changes about the next experiment

Do not further tune this two-component interaction against the same few families. Five assay-confounded families are weak supervision for per-target conditional effects, and the shuffled-label control provides no persuasive evidence that this version learns useful biology. Additional independent interventional contexts, rather than more observational A/B/C control reconstruction, would improve identifiability.

Prioritize a **distribution diagnostic before another submission**. Use cross-fitted real treatment marginals conditioned on library-depth/state strata as a diagnostic oracle; preserve the joint ranks of empirical control cells while changing only source-supported per-gene quantiles. Compare that with the same means/pooled proportions and compact emitter, including a no-perturbation arm and shuffled target assignments. Test both 200 and 400 genuine cells. Require calibrated control-null DE calls, acceptable expression numerators, and improved fidelity/overlap across adequately powered held-family treatments. Oracle-treated distributions remain diagnostic, not legal inputs for hidden competition treatments.

Only if that diagnostic closes the oracle gap should we fit a deployable cross-context transport of **treated-versus-control quantile changes**, with gene-specific dropout/variance changes, context support checks and shrinkage to V13. This is more specific than globally restoring variance: uniform gamma/full-variance approaches already failed V12. Keep means, distribution changes and source selection as separate frozen ablations; validate before packaging. These are proposed next tests, not completed improvements or a promised leaderboard gain.

## Checks, provenance and limitations

Four scientific tests passed: a known context sign reversal, invalid/missing measurements not voting zero, poor-coverage abstention, and held-family label poisoning not changing predictions. The independent audit passed **{audit['cache_files_checked']} cached profile files**, all planned original cell identities and disjoint splits, and **{sum(r['sampled_target_gene_pairs'] for r in audit['ridge_coordinate_replays'])} independently recomputed ridge coordinates**. It replays scalar weighted normal equations, verifies coefficient centers/own-target exclusion, confirms family-balanced aggregations/case counts, and verifies frozen V13 numerical dependencies remain unchanged. It does not independently rerun every DE test or establish hidden biological correctness.

The original plan's freshness inventory missed V7 `cell_split_*.json` history. The separate `metadata_correction.json`, produced before numerical selection/count outcomes were inspected, supersedes that annotation. No identities, seeds, numerical methods or candidates changed. Eight of nine development target names were previously count-tested. Of four reserved target names, only MYRF was globally unused; CAT, DCAF16 and PSMB9 were not. Original freezes are retained to preserve the audit trail.

Controls are disjoint and split within batches, but treatment truth pools source batches at family level. Assay/lineage and residual batch composition remain confounded. Public small panels are noisy; no official private anchor set is available. No new biological dataset was acquired for V14. The 1% expression tolerance is an exploratory guard, not statistical noninferiority.

### Reproduction and files

Run the experiment commands in [V14_STRATEGY.md](V14_STRATEGY.md). With completed hash-verified caches, assessment is preserved once written; it refuses overwriting. Then run:

```cmd
.\\.competition_env\\Scripts\\python.exe audit_v14_experiment.py
.\\.competition_env\\Scripts\\python.exe document_v14_experiment.py
```

Root scripts include preparation, interaction fitting, execution, assessment, history audit, independent audit, documentation and focused tests. `dispatch_v14_oracle.py` only schedules the unchanged oracle evaluator after each model cache finishes; `oracle_dispatch.json` pins this wrapper. Large model/profile caches are ignored by Git. Plans, source/model sidecars, scores, selection, audit and report hashes are retained under `experiments/v14_context_interactions/analysis/`. `report_manifest.json` pins this document's generator and input/output bytes.
'''
    output=ROOT/'V14_RESULTS.md';output.write_text(text,encoding='utf-8',newline='\n')
    marker='<!-- v14-research:begin -->';end='<!-- v14-research:end -->'
    block=f'''{marker}

## V14 context interactions and treatment-population diagnostic

Tested conserved V13 responses plus context interactions, with shuffled-context and amplitude controls. The 72-population, 400-cell development comparison rejected the interaction model: only reach improves, while PDS drops from 0.8333 to 0.7222 and other metrics worsen. Missing HCT116 NMAE and tiny expression denominators limit local power. V13 remains the best recorded submission; no V14 package was created and reserved algorithm verification was not opened.

Nine disjoint 200/200-treatment oracle splits separate response and population limitations. Measured H1 moments reduce NMAE from 0.8910 to 0.2401, but compact oracle cells still have Jaccard 0.0260 versus 0.2963 for real treatment replicas. Next tests should isolate conditional distributions rather than retune the same response head. Four focused tests and the independent split/hash/coefficient/aggregation audit passed. Methods: [V14_STRATEGY.md](V14_STRATEGY.md). Outcomes, limitations and commands: [V14_RESULTS.md](V14_RESULTS.md). Root scripts and small evidence are versioned; V13 numerical dependencies remain unchanged.

{end}'''
    readme=ROOT/'README.md';content=readme.read_text(encoding='utf-8')
    if marker in content:
        a=content.index(marker);b=content.index(end,a)+len(end);content=content[:a]+block+content[b:]
    else:content=content.replace('<!-- project-retrospective:end -->','<!-- project-retrospective:end -->\n\n'+block,1)
    readme.write_text(content,encoding='utf-8',newline='\n')
    files=['plan.json','execution_freeze.json','assessment_freeze.json','metadata_correction.json','screen_summary.csv','screen_result.json','development_scores.csv','development_panels.csv','development_summary.csv','development_evaluation.json','oracle_scores.csv','oracle_panels.csv','oracle_summary.csv','oracle_evaluation.json','selection.json','independent_audit.json']
    build.write_json(A/'report_manifest.json',{'inputs':{str((A/n).relative_to(ROOT)):build.sha256(A/n) for n in files},
        'generator_sha256':build.sha256(Path(__file__)),'output_sha256':build.sha256(output),'report':'V14_RESULTS.md','decision':'retain_v13_no_v14_submission',
        'scope':'Public raw exploratory measurements; no private-score prediction; metadata correction supersedes original freshness annotation.'})
    print('Wrote V14_RESULTS.md, README section, and report_manifest.json')


if __name__=='__main__':main()
