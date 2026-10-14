"""Reproduce the post-V20 evidence review; no fitting, inference or upload.

Only new review artifacts are written. --check verifies the deterministic outputs
against the current input files instead of writing. Python standard library only.
"""
from pathlib import Path
import argparse, csv, hashlib, io, json, math

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'analysis/post_v20_review'
INPUTS = {}
METRICS = ('pds', 'mse', 'nmae', 'fid', 'reach', 'jac')


def read(relative):
    data = (ROOT / relative).read_bytes()
    INPUTS[relative] = hashlib.sha256(data).hexdigest()
    return data.decode('utf-8-sig')


def load(relative):
    return json.loads(read(relative))


def csv_text(rows):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def markdown(rows, columns):
    lines = ['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---'] * len(columns)) + ' |']
    for row in rows:
        values = []
        for key in columns:
            value = row[key]
            values.append(f'{value:.6f}' if isinstance(value, float) else str(value))
        lines.append('| ' + ' | '.join(values) + ' |')
    return '\n'.join(lines)


NEW = [
    ('V14', 'Two-component control-rank interactions; shuffled-context and amplitude controls; 200/200 treatment oracle', 'Interaction rejected; oracle identifies mean and population deficiencies', 'rejected + diagnostic', 'V14_RESULTS.md;v14_context_methods.py'),
    ('V15', 'Nested V13/Student/source-specialist ensemble; training-family reliability; depth-stratified source quantiles', 'All four deployable count candidates fail; reserve unopened', 'rejected', 'V15_RESULTS.md;V15_STRATEGY.md;v15_ensemble_methods.py;v15_population_methods.py'),
    ('V16', 'Response recalibration, molecule transport, empirical residual emission', 'Exact V13 effects + residual10 released; best published result', 'published', 'V16_RESULTS.md;V16_LEADERBOARD_RESULTS.md;v16_methods.py'),
    ('V17', 'Cross-fitted target-wise scalar response-strength CV; generic/shuffled controls', 'Learned variation fails release guards; no basal-state-dependent response function fitted', 'rejected', 'V17_RESULTS.md;V17_STRATEGY.md;v17_methods.py'),
    ('V18', 'PIE pretrained mean/confidence experts; bounded blends and matched controls', 'Tested checkpoint/blends fail; not an exhaustive test of all PIE training variants', 'rejected', 'V18_RESULTS.md;V18_STRATEGY.md;v18_methods.py'),
    ('V19', 'Control covariance operator on 16 programs; bounded transported response', 'Fails broad improvement and amplitude comparison', 'rejected', 'V19_RESULTS.md;V19_STRATEGY.md;v19_methods.py'),
    ('V20', 'Paired Stack treatment/control calibration; bounded mean-response expert', 'Local development/reserve pass; published score loses to V16', 'published regression', 'V20_RESULTS.md;V20_LEADERBOARD_RESULTS.md;V20_STRATEGY.md;v20_methods.py'),
    ('V21', 'Vocabulary-compatible active Stack programs, matched shuffle and occupancy reweighting', 'Joint candidate fails; reserve and production unopened', 'rejected', 'V21_RESULTS.md;V21_STRATEGY.md;v21_methods.py'),
]

PUBLISHED = {
    'First': 'Conservative public K562 response transfer',
    'V2': 'Donor-calibrated CD4 hybrid',
    'V4': 'Lineage-guided sources and response coverage',
    'V5': 'Atlas/moment calibration',
    'V6': 'Fidelity-guarded response calibration',
    'V7': 'V4 effects retained; different sampling',
    'V9': 'Activation-preserving profile calibration',
    'V10': 'Pooled interventions, dose/shrinkage and compact populations',
    'V11': 'Consensus coverage and dispersion calibration',
    'V13': 'Assay-aware unknown-knockdown calibration',
    'V16': 'V13 responses with 10% empirical residual emission',
    'V20': 'Bounded Stack mean-response ensemble',
}

BODY = r'''
## 1. What the task requires

The inputs for anonymous contexts A/B/C are unperturbed control cells. The challenge asks for 400 predicted cells for each of 300 CRISPRi targets in each context, on the fixed 18,533-gene axis: 360,000 cells in one raw-count package. The 300 genes are perturbation identities, not 300 measured post-treatment inputs available for a hidden cell. Learning other genes from those 300 in controls addresses association, whereas this task requires intervention transfer to unmeasured contexts.

The evaluation and FAQ audits remain relevant: use raw nonnegative integer counts, exact axes/counts, and the official scorer. The directly perturbed gene does not contribute to the score. All six CLI component scores, including the fields named MSE and NMAE, are higher-is-better after scaling; underlying raw errors are lower-is-better. Historical rank changes cannot compare model quality as participants submit. Local unanchored metrics cannot be averaged and called the official score.

## 2. How the project progressed

The table above includes every recorded publication. The detailed inventory below includes unsubmitted and rejected experiments, which matter just as much when choosing what to do next.

1. **Controls, correlations and clustering.** Full-transcriptome and 300-target rank analyses established that the contexts differ and measured their cell dispersion. External-control mapping informed source hypotheses. Marker genes and proximity do not establish exact tissue identity or prove that a source transfers CRISPR responses well. This exploratory work was useful but did not supply hidden treatment labels.
2. **Initial transfer through V4.** K562 transfer, then donor-aware CD4 and broader lineage/source policies substantially improved performance. Several aspects changed together; the leaderboard cannot attribute their gains to one isolated mechanism. V3's NNLS mixture was built without a recorded publication.
3. **V5-V9: calibration and early validation limitations.** Better averages or five improved components did not guarantee a better overall score: fidelity losses dominated V5/V9. V7's accepted package retained V4 effects after neural/ridge alternatives failed. V8 repaired common-control discrimination checks but was not published. Correct formatting and a passing data audit establish integrity, not biological accuracy.
4. **The 300-gene reconstruction networks.** We trained ridge models and 24 neural configurations across contexts/objectives: shallow, deep, bottleneck and residual networks. Selected raw-count networks reduced error roughly 26-32% against depth baselines. Normalized improvements over ridge were modest. The raw objective also changed, so this does not isolate architecture superiority. These are control reconstruction results, not evidence of causal CRISPR prediction.
5. **Cross-modal/embedding pilot.** Tahoe/PCA response heads and generated-cell blends were tested and rejected. The advertised proteomic checkpoint was absent. We did not train a full cross-modal foundation model, so that broad family is not disproved.
6. **V10: the major recent breakthrough.** Pooling measured interventions across sources, correcting dose, shrinking after pooling and changing count emission lifted the previous best from -0.050674 to 0.116754. Response estimation and emission changed together; a full causal attribution is unavailable. The practical lesson is that measured perturbation evidence and population construction mattered more than adding an architecture to control reconstruction.
7. **V11-V13: smaller improvements.** Consensus coverage/noise calibration improved V11; V13's deployed change handled unquantifiable knockdown. Low-rank response programs, robust alternatives, extra K562 platform data and pathway priors did not produce accepted replacements. Student-t pooling was later revisited in V15; it is not an untouched lead.
8. **V12/V14-V15: distribution, interactions and ensembles.** Empirical/gamma-Poisson mixtures, more noise, thinning/immigration and response-linked variation were tried. Null realism sometimes improved while treatment metrics deteriorated. V14's target-specific rank interactions failed shuffled/context controls. V15 already tested nested expert weighting, source reliability and marginal distribution transfer; simple majority voting is not a new experiment.
9. **V16: the current champion.** Keep V13 average response estimates; combine 0.25 expected counts, 0.65 Poisson counts and 0.10 transformed empirical controls, followed by integer/moment fitting. Libraries are exact; fitted mean/bulk moments have tolerances. This small change improved the published result. It is not a fully learned treated-cell distribution.
10. **V17: learned response-strength variation.** Public treatment halves estimated excess variance along a response direction. A positive scalar gain varied all response genes together. The candidate failed the learned-versus-shuffled guard. Its scalar gain was sampled independently of the template cell's basal state; it did not learn a conditional state-response function.
11. **V18-V19: external knowledge and covariance.** PIE mean/confidence blends and a 16-program covariance transport operator failed their frozen criteria. These are valid negative results for the tested implementations, not proof that all knowledge models or transport methods fail. Control covariance alone does not identify an intervention response.
12. **V20-V21: pretrained conditional population information.** Stack mean corrections passed public algorithm verification but lost on the hosted score. Original low-vocabulary-coverage joint arms abstained; V21 repaired that and actively tested conditional program changes and occupancy reweighting. The active joint mechanism still failed matched shuffles/mean-only comparisons. V21 changed a limited empirical residual inside a mostly fixed emitter, not the entire learned population distribution.

## 3. What the evidence says our gaps are

**Transfer validation is weaker than it looks.** Most raw-cell benchmarking uses H1, HCT116 and HEK293T. More targets, seeds or thousands of generated populations do not add independent cell lineages. Earlier public biology has been repeatedly inspected. Some HCT116 NMAE panels are ineligible; checkpoint study overlap is not fully audited. V20's reserve pass followed by a hosted loss is direct evidence that this benchmark did not predict the desired transfer. An untouched algorithm case is useful, but it is not automatically independent biological validation.

**Both response estimates and population construction are deficient.** V14 already ran nine disjoint 200-calibration/200-truth oracle cases across H1/HCT116/HEK293T. In H1, giving the compact generator measured treatment moments reduced NMAE from 0.8910 to 0.2401, yet Jaccard was only 0.0260 versus 0.2963 for independent actual-treatment replicas. HEK293T had a similar warning; HCT116 was severely underpowered. These cases reused development biology, had estimated moments and limited truth power. They motivate a stronger diagnostic, not a claim that every extra synthetic DE call is a proven false positive.

**The model still mixes state, efficacy and response.** V13/V16 largely share downstream response directions across contexts, with expression gates changing amplitude/support. The earlier complete-array audit found no opposite-sign effects between contexts in that baseline. V20 did introduce context corrections and sign changes, so that limitation must not be attributed unchanged to V20. Distinct control clusters alone do not tell us how the same intervention changes each cluster.

**The population mechanism remains restrictive.** V16 preserves a little real-control variation; V17 added a target-level scalar; V21 added limited projected residuals. We have not fitted the complete proposed model that jointly separates source cell state, continuous perturbation strength, shared treatment programs and target-specific state-dependent changes, then transfers the conditional population distribution. This was suggested in earlier research notes but has not been fully executed. Repeating generic variance inflation would repeat a failed idea.

**The private cause is not identified.** V20's pre-emission absolute effect mass did not shrink uniformly; it increased slightly in B/C. The score loss therefore does not justify simply increasing dose. Aggregate hosted metrics cannot identify the failing hidden genes, contexts or responder fractions. Neither the project files nor papers reveal what every leading team does.

## 4. Next tests, in order

### Test A: expand the existing oracle into a factorial error diagnosis

Freeze accessible target/study panels before scoring. Cross two mean sources (V16 versus measured treatment-calibration moments) with two relative population shapes (V16 versus independently estimated treatment residual structure), plus an independent real-treatment replica. Train and truth cells must be disjoint; controls, library handling and evaluation references must be paired. Exclude the direct target from diagnostic learning/readouts where appropriate. Treat measured residual transplantation as a diagnostic oracle, not a deployable predictor.

Where possible use genuine 400-cell calibration and 400-cell truth groups: this requires at least 800 distinct treated cells per condition, plus sufficient independent controls. The calibration treatment draw can also serve as the real-replica comparator independent of truth; a third replica independent of calibration requires another disjoint draw. If donor-stratified groups cannot support that, report the smaller power explicitly; do not duplicate cells to manufacture a 400-cell truth panel. Freeze library/moment tolerances and record changes introduced by fitting. If shape and mean cannot be isolated exactly, quantify the remaining mismatch rather than calling them fixed.

**Decision:** estimate whether improved means, improved shape, or their combination gives the largest reproducible gain in DE reach/overlap without breaking expression/fidelity. This selects the next model component. The V14 oracle was an initial version of this question; this larger factorial study is a proposed extension.

### Test B: state-dependent continuous-response pilot

Start small and interpretable, using raw source treatment cells and guide/batch metadata where available. Fit a control-derived state encoder using training controls only. Learn continuous response strength and a few coordinated downstream programs from treated cells, excluding the 300 panel genes from the diagnostic feature axis. Separate a shared treatment/stress component from target-specific fingerprints and their interactions with state; handle guide and assay differences explicitly. Do not assume every CRISPRi cell is either a perfect knockout or a true nonresponder.

One modelling hypothesis is a baseline state plus a shared treatment component plus `response_strength * (target_program + target_by_state_interaction)`. This is a statistical transfer model, not identification of an individual cell's unobserved counterfactual. Whole studies/families must be excluded from fitted effects and latent parameters during evaluation. A cell-cycle shift may be a real treatment effect; avoid matching it away or regressing it out merely because it separates populations.

For A/B/C we observe only control states. Response-strength distributions and state interactions must be predicted from public evidence with explicit uncertainty; we cannot measure private responder fractions from those controls. Weakly supported contexts/targets should return V16. First isolate the new joint structure with V16 aggregate mean targets; allow learned mean changes only in a separately frozen second stage.

Freeze these comparison arms before scoring: exact V16; V17-style scalar strength; state-independent learned programs; state-dependent programs; the same model with state labels shuffled within batch; and an amplitude-matched control. Use compatible noise/count constraints across arms. A model must beat matched shuffles and simpler arms on actual generated-count DE metrics and null tests, not only fit training projections or spread the scatter plot.

### Test C: strengthen context coverage before trusting a release

The most immediately useful data lead is donor-resolved raw CD4 resting/stimulated counts, currently represented mostly by summary effects in our bank. The CIPHER manifest was already collected. Verify raw layers, donor/guide IDs, target overlap and cell power before download/fitting; deduplicate the study against our current CD4 evidence. This may add cell-level validation, not a new biological study. Exclude all related CD4 data from model fitting when scoring it.

iPSC/KOLF raw data could provide a genuinely different context, but current browser exports are z-score summaries and most audited donor/target groups lack 400 cells. Do not use summaries as counts or pooled donors as independent donor tests. Download a feasible supported raw subset only after checking access, cell power and axes. K562 Flex is another assay in an existing lineage, not a substitute for context diversity.

The PS Jurkat dataset is not an easy drop-in source: its genome-scale experiment has a targeted 374-gene readout and a median 13 guides per cell. Multi-guide modelling and limited output coverage prevent calling it a ready full-transcriptome single-target validation set. It can inform response-strength methodology if handled accordingly.

### Test D: freeze the release rule and a genuinely new check

Use the unchanged official count scorer, separate raw errors from scaled scores, and report metric eligibility, absolute expression numerators/denominators, DE call counts and per-family/target errors. Estimate paired uncertainty at the biological target/study level; seeds are repeated draws. Small numbers of independent contexts limit confidence even with many cells. Null, shuffle and amplitude controls remain mandatory for interpreting the proposed pilot.

Select model/hyperparameters on development only. Evaluate a single locked candidate on a new study/donor/family reserve. Audit pretraining overlap if a pretrained model is used. Report missing metrics instead of imputation; underpowered panels cannot establish six-metric validation. Use the same-panel published scores for historical comparisons, without inventing unknown private anchors. Do not relax failed historical policies retrospectively. Any new guard/tolerance must be written before its outcomes and justified against measurement variability.

Keep V16 as champion until a tested candidate earns a release. A gain is possible, but none of these proposals establishes a higher leaderboard score yet.

## 5. Literature: what is useful and what has not been done

[Song et al., perturbation-response score (PS)](https://www.nature.com/articles/s41556-025-01626-9) supports continuous partial-response estimation rather than assuming binary knockout. We already cited it in V17, but implemented only a simpler cross-fitted scalar variance estimator. Its genome-scale Jurkat study uses high-MOI multi-guide cells and a targeted readout; it is not equivalent to our required single-target full-transcriptome prediction problem.

[RNA fingerprinting author implementation](https://github.com/satijalab/rna-fingerprinting) estimates response fingerprints after accounting for baseline factors. Our pinned author code was collected but never fitted. It is a potential source-effect extraction tool; mapping an observed treated query to a dictionary does not predict an unobserved treatment from controls. A separately validated state/efficacy transfer model is still needed.

[Systema](https://www.nature.com/articles/s41587-025-02777-8) documents how systematic differences between controls and treatments can make generic predictions look strong under common metrics. We have already used generic/shuffled diagnostics. The proposed extension is explicit decomposition during learning, with those controls retained, not a claim that a new metric replaces the competition scorer. This paper studies unseen-perturbation benchmarks; its lesson is relevant but its task is not identical to unseen-context transfer.

These are established research ideas and testable adaptations, not a claim of a novel published method or a proven explanation of other teams' ranks. Fresh checks on 2026-10-14 used primary abstracts/method excerpts and author documentation; some publisher/PMC pages blocked full-text access. No new raw training data, weights, model fitting or submission was produced during this review.

## 6. Reproduce this review

Run from the project directory in Windows CMD:

```bat
.\.competition_env\Scripts\python.exe review_post_v20_experiments.py
.\.competition_env\Scripts\python.exe review_post_v20_experiments.py --check
```

The standard-library script reads all 12 recorded publication receipts, checks the six-component arithmetic and panel/anchor consistency, merges the earlier scope inventory with V14-V21, hashes the cited local evidence, and produces this document plus small CSV/JSON evidence under `analysis/post_v20_review/`. `--check` verifies deterministic output bytes. Prior numerical code, frozen plans and submission artifacts are not changed. This validates the review's traceability and arithmetic, not the proposed models' performance.
'''


def build():
    old = list(csv.DictReader(io.StringIO(read('analysis/project_retrospective/published_history.csv'))))
    records = [(int(r['version']), r['label'], r['record']) for r in old]
    records += [(16, 'V16', 'experiments/v16_control_preserving/analysis/leaderboard_result.json'),
                (20, 'V20', 'experiments/v20_stack_transfer/analysis/leaderboard_result.json')]
    assert len(records) == 12
    history, anchor = [], None
    best = -math.inf
    previous = None
    for version, label, relative in sorted(records):
        record = load(relative)
        scores = record['metrics']
        assert record['final_status'] == 'published'
        assert math.isclose(sum(scores[k] for k in METRICS)/6, record['overall'], abs_tol=1e-12)
        current = (record['partition'], record['panel'], record['anchor_set'])
        if anchor is None:
            anchor = current
        assert current == anchor, ('Noncomparable panel/anchors', label)
        total = record['overall']
        history.append(dict(version=version, label=label, entry_id=record['entry_id'], overall=total,
                            change_from_previous_upload='' if previous is None else total-previous,
                            change_from_previous_best='' if best == -math.inf else total-best,
                            new_project_best=total > best, method=PUBLISHED[label], record=relative,
                            **{k:scores[k] for k in METRICS}))
        best, previous = max(best, total), total
    inventory = list(csv.DictReader(io.StringIO(read('analysis/project_retrospective/experiment_inventory.csv'))))
    for row in inventory:
        if row['area'] == 'V13 robust':
            row['outcome'] += '; nested robust ensemble revisited and rejected in V15'
            row['evidence_status'] = 'historical lead; follow-up rejected'
    inventory += [dict(zip(('area','tested_scope','outcome','evidence_status','evidence'), row)) for row in NEW]
    evidence = set()
    for row in inventory:
        for name in row['evidence'].split(';'):
            read(name)
            evidence.add(name)
    for name in ('EVALUATION_CRITERIA_AUDIT.md','FAQ_AUDIT.md','POST_V11_RESEARCH.md',
                 'POST_V17_RESEARCH.md','v20_production.py','v13_quality_methods.py',
                 'v8_methods.py','analysis/post_v11_research/data/fingerprint_readme.md'):
        read(name)
        evidence.add(name)
    winner = max(history, key=lambda row: row['overall'])
    assert winner['label'] == 'V16'
    v16, v20 = [next(row for row in history if row['label'] == label) for label in ('V16','V20')]
    delta = v20['overall']-v16['overall']
    reach = (v20['reach']-v16['reach'])/6
    comparison = [dict(metric=key, v16=v16[key], v20=v20[key], change=v20[key]-v16[key]) for key in ('overall',)+METRICS]
    report = f'''# Experiment review and next tests after V20

Review date: 2026-10-14. Evidence snapshot, not new model results.

We have tried several distinct modelling ideas. The bottleneck is not simply a lack of originality: V10/V11 produced large gains, while later algorithms mostly reuse the same public intervention families and limited validation biology. The strongest next direction is to diagnose mean versus population error, then test state-dependent continuous responses with stronger context validation.

**V16 remains the best recorded submission: {v16['overall']:.9f}.** V20 scored {v20['overall']:.9f}, a change of {delta:+.9f} ({delta/v16['overall']:.2%}). Its reach component accounts arithmetically for {reach/delta:.1%} of the net decrease; this is not a causal explanation.

## Published history

{markdown(history, ['label','overall','change_from_previous_best','method'])}

A negative change from previous best means that upload did not beat the champion then. These are comparable scaled scores on the same panel/anchors, not historical ranks. The first upload has no predecessor.

{markdown(comparison, ['metric','v16','v20','change'])}
''' + BODY + '\n## Detailed scope inventory\n\n' + markdown(inventory, ['area','tested_scope','outcome','evidence_status']) + '\n'
    read(Path(__file__).name)
    manifest = dict(review_date='2026-10-14', purpose='Historical evidence and proposed tests; no new fitted model',
                    published_submissions=len(history), best_label=winner['label'], best_score=winner['overall'],
                    common_partition_panel_anchors=anchor, v20_change_from_v16=delta,
                    reach_share_of_net_score_decline=reach/delta,
                    local_evidence_files=len(evidence), input_sha256=dict(sorted(INPUTS.items())))
    return {
        ROOT/'POST_V20_RETROSPECTIVE.md': report,
        OUT/'published_history.csv': csv_text(history),
        OUT/'experiment_inventory.csv': csv_text(inventory),
        OUT/'v16_v20_comparison.csv': csv_text(comparison),
        OUT/'review_manifest.json': json.dumps(manifest, indent=2, allow_nan=False)+'\n',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    outputs = build()
    for path, text in outputs.items():
        data = text.encode('utf-8')
        if args.check:
            assert path.read_bytes() == data, 'Review artifact differs: '+str(path)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
    print(('Verified' if args.check else 'Wrote'), len(outputs), 'review artifacts; 12 publications, V16 retained. No model or submission changed.')


if __name__ == '__main__':
    main()
