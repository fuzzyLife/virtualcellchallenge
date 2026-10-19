"""Create the reproducible, same-panel publication review and scientific trajectory figure."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v65_closed_response/analysis'

def main():
 r=json.loads((A/'published_history.json').read_text());rows=r['records'];by={x['version']:x for x in rows}
 table=['| Submission | Main change | Overall | Change vs prior best | New best? |','| --- | --- | ---: | ---: | :---: |']
 for x in rows:
  delta='initial' if x['change_from_previous_best'] is None else f"{x['change_from_previous_best']:+.6f}"
  table.append(f"| V{x['version']} | {x['method']} | {x['overall']:.6f} | {delta} | {'yes' if x['new_project_best'] else ''} |")
 best=[];current=-1e9
 for x in rows:current=max(current,x['overall']);best.append(current)
 fig,axes=plt.subplots(1,2,figsize=(13.5,4.7),gridspec_kw={'width_ratios':[1.2,1]})
 pos=list(range(len(rows)));labels=[f"V{x['version']}" for x in rows]
 axes[0].plot(pos,[x['overall'] for x in rows],color='#6b7c93',marker='o',label='Published submission');axes[0].step(pos,best,where='post',color='#008b7a',label='Best so far')
 axes[0].axhline(0,color='#dddddd',linewidth=1);axes[0].set_xticks(pos,labels,rotation=60);axes[0].set_ylabel('Official scaled overall score');axes[0].set_title('All 18 published submissions, same panel');axes[0].legend(frameon=False,fontsize=8)
 metrics=('pds','mse','nmae','fid','reach','jac');deltas=[by[63][k]-by[24][k] for k in metrics]
 axes[1].barh(metrics,[v/6 for v in deltas],color=['#008b7a' if v>=0 else '#c64d48' for v in deltas]);axes[1].axvline(0,color='#777777',linewidth=.8);axes[1].invert_yaxis();axes[1].set_xlabel('Contribution to overall change (score change / 6)');axes[1].set_title('V63 improvement over V24: +0.007425')
 fig.suptitle('Scores establish our trajectory; historical ranks depend on other teams',fontsize=12);fig.tight_layout();fig.savefig(A/'published_trajectory.png',dpi=160);plt.close(fig)
 text='''# Submission history and the next experiment

The project has **18 recorded published submissions**. Every receipt uses validation panel `vcc2026-val-1` and the same anchor set. The table compares each submission with the best earlier project score. Historical ranks are affected by other teams and dates; score improvement is the reliable comparison available here. All six CLI components are scaled accuracies, with higher values better.

'''+ '\n'.join(table)+'''

## What has actually helped

- V2 and V4 improve the initial weak transfer through donor/context calibration and response coverage, but overall remains negative.
- V10 changes pooled intervention response and population transport together. Its +0.167428 gain over the previous best V4 comes mostly from fidelity (+0.110253 contribution to overall) and PDS (+0.052860). These are arithmetic contributions, not an isolated causal test of the components.
- V11 adds consensus coverage/dispersion: +0.020250 over V10, improving five components while losing expression accuracy.
- V13, V16 and V24 add smaller gains. V24 improves all six components over V16; its context-matched CD4 response is a useful conservative reference.
- V63 gains +0.007425 over V24, driven by expression (+0.008484 contribution), with further PDS, NMAE and overlap gains. Fidelity and reach decline. Its changed response representation is worth investigating while keeping the established generator as a controlled comparator.

## What has not earned a better published score

V5/V6/V7/V9 do not beat V4. Stacking V20, retrieval/population-view correction V36, empirical response prior V50 and fixed V24/V36 blending V55 do not beat V24. State-response transfer V61 regresses to 0.116763. This does not prove that every ensemble or neural model fails; these particular changes did not transfer successfully to the hosted panel. Combining the best scalar metric from different submissions cannot construct a valid population prediction.

Many version numbers are local experiments with no published score. Their evidence is kept separate:

| Experiment family | Local finding |
| --- | --- |
| V12/V14/V15/V17-V23 | Response calibrations, specialist/ensemble and population variants; see their individual result reports. Published descendants V13/V16/V20/V24 provide the hosted evidence. |
| V25-V29 | Intervention geometry, neural decoder and control-preserving transport; no candidate passed its release checks. Better control-null behavior alone did not ensure treatment fidelity. |
| V30-V34 | STRING/construct features, guide-aware correction, paired transfer, tabular foundation model and flow model; failed mean screens or shuffled-control comparisons. |
| V35-V37 | Retrieval/residual training, population views and component ensembles; V36 was published but lost to V24; V37 blends failed local guards. |
| V38-V44 | Activity/direction heads, KOLF data, contrastive identity, LINCS fingerprints and increased residual dispersion; failed relevant gates or target-shuffled controls. |
| V45-V50 | Residual reduction, population mixtures, corroboration, exact-margin projection/rounding and response priors. Numerical improvements did not establish better biology; published V50 lost to V24. |
| V52-V55 | Joint mean/rank neural prediction, feature correction and decoder/blend diagnosis. Joint heads did not provide reliable transfer; published V55 lost to V24. |
| V56-V58 | PIE evaluation with genuine saved training-response evidence and separate heads. Mean gains did not beat adequate controls or produce a validated count-decoder gain. |
| V59-V61 | State abundance/within-state transfer and depth adjustment. Public transfer checks failed; published V61 regressed. |
| V62-V64 | Assay/common-response transformation, identity guard and scalar magnitude control. V63 improved privately despite failed public gates; V64 failed count checks and was not uploaded. |

The 300-gene control-expression reconstruction networks model observational coexpression, which does not establish causal CRISPR response transfer. The repeatedly used public cohorts also provide weak validation: 44 case configurations reduce to 30 unique target/context pairs and 34 distinct measured cell sets. Seeds and overlapping cell sets do not add independent biological studies. Local PDS often ranks candidates differently from the hosted panel; zero-shot baseline geometry guards protect our prior predictions, not ground truth.

## V65 decision

Test available-axis source normalization and normalized source perturbation profiles before assay-space transfer. Keep the V24 count generator and geometry guard, compare with **V63** as the current best reference, and include V24 and a shuffled-target diagnostic. Select only with the frozen desired-mean rule before opening new counts, then report the complete locked count comparison. This targets a concrete units/composition assumption revealed by the successful V63 estimator rather than adding another unvalidated predictor to a blend.

[Protocol and commands](V65_STRATEGY.md). The new candidate's measured results and final packaging decision are documented separately; no private V65 score is inferred.

![Same-panel publication trajectory and V63 component contributions](experiments/v65_closed_response/analysis/published_trajectory.png)

The history loader asserts score arithmetic and records every source receipt hash in `published_history.json`; `published_history.csv` contains the full component data and historical rank snapshots. Reproduce with `review_v65_history.py` followed by `document_v65_history.py`. Earlier reports remain historical snapshots and are not overwritten.
'''
 (ROOT/'V65_SUBMISSION_HISTORY.md').write_text(text,encoding='utf-8');print('Documented all 18 published submissions and saved the trajectory plot.')
if __name__=='__main__':main()
