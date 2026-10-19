"""Document the fixed V68 alternative population model and its complete diagnostics."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import v68_population_balancing as m
from v57_pie import read,write,sha


def main():
    m.check();r=read(m.A/'assessment.json');complete=read(m.A/'primary_null_assessment.json');assert read(m.A/'results_audit.json')['passed']
    assert not (r['release_supported'] and complete['passed']),'Passing primary must proceed to confirmation before release documentation.'
    studies=('H1','Orion','CD4');arms=('v63','balanced_real','shuffled_weights','uniform_real')
    metrics=('pds','fidelity','reach','jaccard','nmae','profile_square_error');rows=[]
    for metric in metrics:
        row=dict(metric=metric,**{arm:r['summary'][arm][metric] for arm in arms})
        row['change_vs_v63']=row['balanced_real']-row['v63'] if row['balanced_real'] is not None and row['v63'] is not None else None;rows.append(row)
    with (m.A/'component_comparison.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fits=[read(p) for p in (m.A/'fits').glob('*.json')];assert len(fits)==10
    nconverged=sum(f['converged'] for f in fits);scope={}
    for study in studies:
        scope[study]=r['studies']['balanced_real'][study]['profile_square_error']/r['studies']['v63'][study]['profile_square_error']
    labels={'v63':'V63','balanced_real':'Balanced real cells','shuffled_weights':'Shuffled weights','uniform_real':'Uniform real cells'}
    colors={'v63':'#284b63','balanced_real':'#278676','shuffled_weights':'#c27b36','uniform_real':'#8a8a8a'}
    fig,axes=plt.subplots(1,3,figsize=(15,4.7),layout='constrained')
    for j,arm in enumerate(arms[1:]):
        vals=[r['studies'][arm][s]['profile_square_error']/r['studies']['v63'][s]['profile_square_error'] for s in studies]
        xx=np.arange(3)+(j-1)*.17;axes[0].vlines(xx,1,vals,color=colors[arm],alpha=.6);axes[0].scatter(xx,vals,color=colors[arm],label=labels[arm])
    axes[0].axhline(1,color='black',lw=.8);axes[0].set_xticks(range(3),studies);axes[0].set_title('Generated profile error / V63');axes[0].set_ylabel('Lower is better; V63 = 1');axes[0].legend(fontsize=7)
    for j,arm in enumerate(arms):
        vals=[r['studies'][arm][s]['fidelity'] for s in studies]
        axes[1].bar(np.arange(3)+(j-1.5)*.18,vals,width=.17,color=colors[arm],label=labels[arm])
    axes[1].set_xticks(range(3),studies);axes[1].set_title('Raw DE direction fidelity');axes[1].set_ylim(bottom=0);axes[1].set_ylabel('Higher is better; unanchored')
    nulls=r['nulls'];full={v['source']:v for v in complete['records']}
    for j,(field,label,color) in enumerate([('v63_calls','V63','#284b63'),('natural_calls','Uniform real','#8a8a8a'),('balanced_real_calls','Fitted real','#278676')]):
        vals=[full[v['source']][field] if field=='balanced_real_calls' else v[field] for v in nulls]
        xx=np.arange(5)+(j-1)*.2;axes[2].scatter(xx,vals,color=color,label=label,s=25);axes[2].vlines(xx,0,vals,color=color,alpha=.6)
    axes[2].set_yscale('symlog',linthresh=1);axes[2].set_xticks(range(5),['H1','HCT116','HEK293T','CD4 rest','CD4 stim'],rotation=30);axes[2].set_title('Unperturbed-control DE calls');axes[2].set_ylabel('Lower is better; symlog scale');axes[2].legend(fontsize=7)
    for ax in axes:ax.grid(axis='y',alpha=.2)
    fig.suptitle('V68: real-cell population selection versus synthetic counts')
    figure=m.A/'population_tradeoff.png';fig.savefig(figure,dpi=160);plt.close(fig)
    def value(x):return 'undefined' if x is None else f'{x:.6f}'
    table=['| Raw diagnostic | V63 | Balanced real | Shuffled weights | Uniform real |','| --- | ---: | ---: | ---: | ---: |']
    table += ['| '+row['metric']+' | '+' | '.join(value(row[a]) for a in arms)+' |' for row in rows]
    nt=['| Context | V63 null calls | Uniform-real null | Complete fitted null |','| --- | ---: | ---: | ---: |']
    nt += [f"| {v['source']} | {v['v63_calls']} | {v['natural_calls']} | {full[v['source']]['balanced_real_calls']} |" for v in nulls]
    fails=[k for k,v in r['checks'].items() if not v]
    report=f'''# V68 results: preserve real cells, test population selection

The fixed real-cell population model **did not qualify for release**. Paired equal-study generated-profile error / V63 is **{r['paired_equal_study_profile_ratio']:.6f}** ({r['paired_equal_study_profile_ratio']-1:+.2%}). Its raw direction fidelity is **{value(r['summary']['balanced_real']['fidelity'])}**, versus V63 **{value(r['summary']['v63']['fidelity'])}**. The complete fitted null guard passes: **{complete['passed']}**. No V68 submission was generated. V66 remains the highest published project score, 0.150882188.

## What was tested

This is an alternative population generator, not another response blend. The optimizer selects weights on individual real query-control cells to approach V63's desired bulk expression, with capped inclusion probabilities and fixed depth-group quotas. Pivotal sampling selects distinct original cells. Every non-own-gene count remains intact; only the own gene receives an independent 0.1 binomial thinning before fitting/sampling. There is no synthetic Poisson population, mean refit or per-gene count repair.

The desired means are shared with V63; **realized means are not forced to match**. Controls may not contain the states required to express the predicted treatment response. This experiment tests the whole alternative decoder, not a pure variance change with equal realized means. See the [frozen method and limitations](V68_STRATEGY.md).

## Treatment outcomes

{chr(10).join(table)}

PDS, fidelity, reach and overlap are raw higher-better diagnostics. NMAE and profile-square error are lower-better errors. The paired profile ratio averages three study ratios; it differs from the ratio of the displayed mean absolute errors. Undefined case metrics are omitted within their context, and a metric remains undefined if an entire study is unresolved. These values are not the leaderboard's scaled accuracies, and no overall competition score is inferred. Ten target/context pairs are tested across five contexts and three studies. Four arms produce 40 comparison records, including ten hash-verified reused V63 outcomes. Actual treatment populations/references are identical across arms. Predicted groups have 200 or400 cells, matching the saved comparator; resting CD4 uses200.

Per-study profile ratios: {', '.join(f'{s} {scope[s]:.6f}' for s in studies)}. Fixed checks failed: {', '.join(fails)}. No parameter, control arm or mixture was promoted after outcomes. Two-target discrimination panels are coarse and unstable; their ranks cannot establish a hosted gain.

In H1, the two fixed targets HDAC8 and NFE2L1 have 157 and137 reference DE calls; balanced-real predictions make zero downstream calls in each. Restoring natural baseline variation did not reproduce their measured induced transcriptional responses. This is direct evidence for these tested cases, not proof that all occupancy models fail.

## Control-null calibration

{chr(10).join(nt)}

Uniform sampling tests the raw real-cell sampler. The additional fitted zero-response test exercises the complete balancing pipeline. The [separate addendum](V68_NULL_AUDIT.md) was frozen after some treatment outcomes but before any null outcome, because review found that uniform sampling alone did not test optimization. It adds a requirement without altering the original candidate or guards. These are fifteen null-comparison records across five contexts, including the original ten and five additional fitted-null records; they are not fifteen independent studies.

![Population and null diagnostics](experiments/v68_population_balancing/analysis/population_tradeoff.png)

## What we learned and what remains open

Technical calibration and biological response are separate problems. Selecting intact baseline cells can preserve authentic dispersion, but cannot invent perturbed states absent from the fit-control population. This test does not support replacing V63's induced-response generator with the tested occupancy-only model. It also does not justify a hybrid mixture or stronger weight fit without a new protocol. The relevant next distinction is between baseline variability and variation created by a perturbation; broadening natural variation alone is insufficient in these cases. This experiment does not isolate whether the loss comes from the control-state convex hull, fixed regularization, imperfect upstream response means or finite-sample power; an oracle diagnostic would be needed to separate them.

V67 also failed: its separate network mechanism reduced broad mean-profile error by only0.0018% and lost the target confirmation check. Together, these experiments provide two explicit negative results rather than an unsupported new upload. They do not establish which method will improve the hidden panel. D3 remains unopened; the public biology here has been used before.

## Verification and reproduction

Eight V68 numerical tests pass, including6,000 sampling simulations, objective finite-difference checks and recovery of a known population shift. All five control frames match frozen V63 query/control bulk profiles. Control/reference/treatment partitions are disjoint. All ten weight fits are feasible; **{nconverged}/10** meet the fixed1e-8 dual-gap stopping rule. Maximum final dual gap is **{max(f['duality_gap'] for f in fits):.3g}**. Effective sample sizes range **{min(f['effective_sample_size'] for f in fits):.1f} to {max(f['effective_sample_size'] for f in fits):.1f}**, always at least the number of selected cells. The remaining fits reached the preset iteration limit, so this trial is not evidence that an exactly optimized occupancy model must fail. The reported dual gaps bound suboptimality of the convex surrogate, not uncertainty in the discrete DE scores. No solver or regularizer was changed after outcomes.

The independent audit reconstructs all30 new treatment populations from original control rows, verifies non-own count preservation, selected-row uniqueness and depth quotas, independently reproduces full profiles and count-derived fold changes, checks BH adjustments, and replays four DE metrics including raw-p-value tie ordering for reach. It manually verifies16 Wilcoxon rank coordinates per population. Cached V63 baseline receipts remain bound to their original hashes; old raw-p reach ties are not independently replayed from caches that did not save p-values. Null call totals are reconciled; no independent full regeneration of the old synthetic-null generator is claimed.

Root commands: [V68_STRATEGY.md](V68_STRATEGY.md), [additional null commands](V68_NULL_AUDIT.md), `audit_v68_results.py`, `document_v68_results.py`, and `seal_v68_experiment.py --check`. The original protocol was committed as `78fb103`; the stricter pre-null addendum as `96f5665`. Scientific scripts and compact records are committed locally; large arrays are ignored and regenerable. Prior scientific release seals remain unchanged. No upload or Git push occurs automatically.
'''
    (m.ROOT/'V68_RESULTS.md').write_text(report,encoding='utf-8')
    write(m.A/'final_decision.json',dict(release_supported=False,development_complete=True,
          treatment_comparison_records=40,new_treatment_populations=30,reused_V63_comparisons=10,
          original_null_records=10,additional_complete_pipeline_null_records=5,
          failure_checks=fails,complete_fitted_null_passed=complete['passed'],
          submission_generated=False,D3_opened=False,private_score=None,
          interpretation='No validated improvement. Keep V66 highest published and V63 near-tied reference; do not create a renamed duplicate or claim the natural-control variation models induced treatment variation.'))
    path=m.ROOT/'README.md';text=path.read_text(encoding='utf-8-sig');begin='<!-- v68-real-population:begin -->';end='<!-- v68-real-population:end -->'
    block=f'''{begin}
## V67/V68: two different mechanisms tested; neither released

V67 tested an intervention-feedback network on655 public target/context pairs: only **0.0018%** lower profile error, with failed target confirmation. V68 replaced synthetic counts with constrained selection of intact real control cells (own-gene thinning only). On10 paired target/context cases, generated profile error / V63 is **{r['paired_equal_study_profile_ratio']:.6f}**; treatment release guards fail. Complete fitted null guard passes: **{complete['passed']}**. Both experiments are documented with fixed controls and independent audits; reused public biology is not blind validation. **No V67/V68 submission generated. V66 remains highest published:0.150882188.**

[Network results](V67_RESULTS.md) · [Real-population results](V68_RESULTS.md) · [V68 protocol](V68_STRATEGY.md) · [Complete null addendum](V68_NULL_AUDIT.md). Earlier status blocks below are historical. Previous source models and artifacts are preserved.
{end}

'''
    if begin in text:
        start=text.index(begin);finish=text.index(end,start)+len(end);text=text[:start]+text[finish:].lstrip('\n')
    path.write_text(block+text,encoding='utf-8');print('Documented V68 completed diagnostics; no release.',flush=True)


if __name__=='__main__':main()


