"""Document the completed, independently verified V6 submission and reproduction."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import build_submission as build
ROOT=build.ROOT;OUT=ROOT/'experiments/v6_fidelity_calibrated';A=OUT/'analysis'


def main():
    model=json.loads((A/'reference_model.json').read_text());ready=json.loads((A/'submission_ready.json').read_text());audit=json.loads((A/'independent_validation.json').read_text())
    neural=json.loads((A/'neural_template_selection.json').read_text());assert ready['official_prep_passed'] and audit['passed']
    assert audit['package_sha256']==ready['sha256'] and audit['model_sha256']==build.sha256(A/'reference_model.json')
    summary=pd.read_csv(A/'primary_source_candidate_summary.csv');gate=pd.read_csv(A/'neural_control_gate_summary.csv');context=pd.read_csv(A/'context_generator_ablation.csv')
    source=pd.read_csv(A/'source_generator_scores.csv')
    source=source[source.target_gene.isin(model['primary_source_targets'])]
    keys=['target_gene','direction','seed']
    paired=source[source.candidate==model['candidate']].merge(source[source.candidate=='v4'],on=keys,suffixes=('_v6','_v4'),validate='one_to_one')
    metrics=['fidelity','jaccard','nmae','reach']
    for key in metrics:paired[key+'_difference']=paired[key+'_v6']-paired[key+'_v4']
    target_delta=paired.groupby('target_gene')[[k+'_difference' for k in metrics]].mean()
    target_delta.to_csv(A/'source_target_diagnostics.csv')
    rng=np.random.default_rng(20261007)
    sampled=target_delta.to_numpy()[rng.integers(0,len(target_delta),(10000,len(target_delta)))].mean(axis=1)
    uncertainty=[]
    for i,key in enumerate(metrics):
        values=target_delta[key+'_difference']
        uncertainty.append({'metric':key,'mean_v6_minus_v4':float(values.mean()),'target_bootstrap_low95':float(np.quantile(sampled[:,i],.025)),
            'target_bootstrap_high95':float(np.quantile(sampled[:,i],.975)),'targets':len(target_delta),'target_improvements':int((values<0 if key=='nmae' else values>0).sum())})
    pd.DataFrame(uncertainty).to_csv(A/'paired_source_uncertainty.csv',index=False)
    build.write_json(A/'source_uncertainty_scope.json',{'seed':20261007,'bootstrap_samples':10000,'resampling_unit':'Target; paired batch directions and emission seeds remain together.',
        'scope':'Exploratory percentile intervals on nine previously inspected public H1 targets; not competition confidence intervals or evidence of anonymous context transfer.'})
    fig,axes=plt.subplots(1,3,figsize=(15,4.7),layout='constrained')
    labels=summary.candidate.tolist();pos=np.arange(len(labels));colors=['#197e8a' if x==model['candidate'] else '#8795a7' for x in labels]
    axes[0].bar(pos,summary.fidelity,color=colors);axes[0].set_xticks(pos,labels,rotation=65,ha='right');axes[0].set(title='H1: actual rank-test direction yield',ylabel='Unanchored source fidelity')
    axes[1].bar(pos,summary.expression_ratio,color=colors);axes[1].set_xticks(pos,labels,rotation=65,ha='right');axes[1].set(title='H1: corrected expression ratio',ylabel='Source error ratio (lower is better)')
    for i,method in enumerate(['observed','moment_only','neural_10','neural_25']):
        row=gate[gate.method==method].set_index('context').loc[list('ABC')]
        axes[2].bar(np.arange(3)+(i-1.5)*.2,row.variance_ratio_to_reference,width=.2,label=method)
    axes[2].axhspan(.85,1.15,alpha=.1,color='green');axes[2].axhline(1,color='black',lw=.7);axes[2].set_xticks(np.arange(3),list('ABC'));axes[2].set(title='Disjoint control-guide variance',ylabel='CPM variance / reference variance');axes[2].legend(fontsize=9,loc='upper center',bbox_to_anchor=(.5,-.08),ncol=2)
    fig.savefig(A/'v6_validation.png',dpi=170);fig.savefig(A/'v6_validation.pdf');plt.close(fig)
    grouped=context.groupby(['context','candidate']).predicted_DE_calls.mean().unstack()
    fig,ax=plt.subplots(figsize=(9,4.4),layout='constrained');grouped.plot.bar(ax=ax)
    ax.set(title='Paired actual-context generator ablations: response activity',ylabel='Mean predicted DE calls (six selected targets)',xlabel='Context');ax.tick_params(axis='x',rotation=0);ax.legend(fontsize=8)
    fig.savefig(A/'v6_context_ablation.png',dpi=170);fig.savefig(A/'v6_context_ablation.pdf');plt.close(fig)
    rows=['| Candidate | Fidelity | Reach | Jaccard | LFC NMAE | Expression ratio | PDS |', '|---|---:|---:|---:|---:|---:|---:|']
    for row in summary.itertuples():rows.append(f'| {row.candidate} | {row.fidelity:.4f} | {row.reach:.4f} | {row.jaccard:.4f} | {row.nmae:.4f} | {row.expression_ratio:.4f} | {row.pds:.4f} |')
    settings=json.dumps(model['setting'],sort_keys=True)
    text=f'''# V6: fidelity-guarded submission ready for upload

Created, officially packaged and independently verified at `experiments/v6_fidelity_calibrated/prediction.vcc`. It has **not been uploaded or scored**. V4 remains the best published project submission; no local benchmark establishes a better hidden score.

| Check | Result |
|---|---|
| Shape | 360,000 cells by 18,533 genes |
| Groups | 900 target/context combinations, 400 cells each |
| Stored counts | {ready['prediction_validation']['nnz']:,} positive integer values |
| Maximum cell library | {ready['prediction_validation']['maximum_cell_counts']:,} counts |
| Package size | {ready['bytes']:,} bytes ({ready['bytes']/1024**3:.2f} GiB) |
| Frozen candidate | `{model['candidate']}` |
| Settings | `{settings}` |
| Neural template weight | 0 |
| Source perturbations | 25 targets, two disjoint batch directions, two generation seeds, eight configurations |
| Primary setting-selection targets | {len(model['primary_source_targets'])}; at least 200 distinct evaluation cells in both partitions |
| Counts and axes | Every stored count and all 360,000 libraries checked; exact official gene/label axes |
| Control provenance | Every original row, cell ID, guide and library checked; 400 distinct controls and 8-9 cells per guide in each group |
| Exact count replay | {audit['exact_replayed_groups']} complete groups; {audit['exact_replayed_gene_counts']:,} gene counts |
| Focused tests | Seven passed |
| Official packaging | VCC prep passed; every compressed payload count, CSR index/pointer, axis and label matches generated H5AD |

SHA-256:

```text
{ready['sha256']}
```

## What was selected

A retains all V4 mean-response effects. B/C can only blend existing effects that agree in sign with V5, with a maximum blend of 25%; V4 zeros and conflicting directions remain unchanged. The frozen candidate actually uses the settings above. No A promoter correction or neural knockdown counterfactual was added.

The benchmark separates V4-like source effects, full dual-moment generation, partial stochastic count calibration, V5-like source gain/pooled calibration and small response blends. Source blends change the held-out H1 contrast gain from 0.7 toward 1.5. The production B/C blend additionally combines the frozen lineage estimators and applies the sign/zero mask; its complete 900-response policy is not reproduced by that source-only comparison. Full V5 response replacement is a diagnostic and is ineligible for production. The authoritative production record is `reference_model.json`; `source_selection.json` retains the earlier unrestricted source diagnostic.

Partial calibration, when selected, first generates V4 binomial-thinned/depth-scaled Poisson-added counts. It mixes the emitted mean and pooled probabilities with response-implied desired probabilities at the selected weight, then fits/rounds those profiles while preserving every stochastic cell's resulting library size. Full dual generation, when selected, starts from original individual controls. Mean and pooled targets always come from the same response estimator. Zero calibration retains the V4 generator.

All groups use separately seeded draws of 400 unique real A/B/C controls balanced over the 46 guides and depth. The seed is 20261007. Own-target knockdown remains the explicit 90% assumption, separate from downstream effects; scoring excludes that gene but its library-normalization consequences still matter.

## Measured source validation

Source cells are the pinned expanded public H1 cache, with all original IDs and counts verified. Its 1,536 controls form two batch partitions, each split into 384 template controls and 384 disjoint reference controls; source evaluation groups therefore have at most 384 cells. Estimation and evaluation batches are opposite global partitions. Evaluation-template controls and DE-reference controls are disjoint IDs within the evaluation partition. Every reference group uses distinct actual perturbed cells, capped at 400 and at the template-pool size; smaller groups are not expanded by duplication. Hyperparameters are selected using other targets in each of five outer folds.

The primary targets are: {', '.join(model['primary_source_targets'])}. All smaller groups remain in `source_generator_scores.csv` and `source_candidate_summary.csv` as separately scoped diagnostics.

{chr(10).join(rows)}

These are raw, unanchored source measurements. DE uses actual two-sided asymptotic Wilcoxon tests on per-cell CPM, a control-only >5 CPM gene filter and within-group BH; the own target is removed from reported metrics. Fidelity uses the call/reference budget, reach uses the confidence ordering and 90% prefix purity, and NMAE/Jaccard use the reference significant sets. PDS excludes all official panel-target genes and ranks only the available source-target panel. Expression uses the delete-one jackknife/capped-noise ratio on the source panel. There are no competition baseline/replicate anchors here. Source PDS with fewer than two targets is undefined.

Paired primary-target differences and 10,000 target-cluster bootstrap intervals are recorded in `paired_source_uncertainty.csv`; target-level outcomes are in `source_target_diagnostics.csv`. Both emission seeds and both batch directions stay within the resampled target. These are exploratory intervals on nine previously inspected public H1 targets, rather than uncertainty estimates for the competition score.

Nested production selection metrics: `{json.dumps(model['nested_procedure'],sort_keys=True)}`.

Nested unchanged-source comparator: `{json.dumps(model['nested_v4_comparator'],sort_keys=True)}`.

Full dual fitting reduced response-call budgets markedly in the paired context audit. Production generators must therefore retain at least 90% of the V4 call budget in every context on the six selected targets. This is a conservative activity guard, not proof that these calls are biologically correct.

Adoption requires competitive held-out fidelity and overlap (maximum 0.005 raw decline) and no corrected-expression-error increase. The within-training comparison uses fidelity + overlap + 0.5*(1-expression ratio), a fixed source proxy rather than official anchor scaling. Outer results are weighted by the number of primary targets. The table's full nine-target expression ratio differs from the weighted fold ratio because the jackknife correction depends on each panel's prediction spread. If nested adoption is unsupported, the frozen model falls back to V4 response/generation. All settings were locked before generating this package.

The H1 sources and earlier source results have been inspected in previous work. This is a comparative ablation on those data, not a new independent dataset. Source calibration cannot establish anonymous A/B/C efficacy. There are no cell-level CD4 donor references in the summary release to validate this complete generator in A.

## Neural template decision

The selected train-guide-only neural models were tested on disjoint validation-guide templates and test-guide reference controls, twice per context. Training-guide moments calibrate the candidates. The test guides were inspected previously, so these are comparative control diagnostics. 25% neural weights fail the distribution gate in all three contexts; 10% fails it in C. Small neural weights can improve pooled-control error while reducing variance, illustrating why reconstruction MSE alone does not justify them for DE. A/B 10% pass the control gate but have no measured perturbation-level benefit. Production therefore retains empirical residual variability with neural weight zero.

## Actual-context generator ablation

The four response/generator combinations and two partial corrections were also generated for six selected targets per context with identical original-control draws. Actual Wilcoxon/BH calls and dispersion are retained in `context_generator_ablation.csv`. This tests emitted response activity and numerical behavior against supplied controls; it has no hidden response truth and cannot attribute V5's fidelity loss to a component.

![Source and neural validation](experiments/v6_fidelity_calibrated/analysis/v6_validation.png)

![Paired context response activity](experiments/v6_fidelity_calibrated/analysis/v6_context_ablation.png)

## Reproduce from Windows Command Prompt

Requires the unchanged competition controls, V4/V5 arrays/source caches, public H1 selected counts, neural evaluation weights and pinned reviewed upstream clone. No new source download or hidden response labels are used.

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe calibrate_v6_submission.py
.\\.competition_env\\Scripts\\python.exe evaluate_v6_neural_templates.py
.\\.competition_env\\Scripts\\python.exe audit_v6_generators.py
.\\.competition_env\\Scripts\\python.exe audit_v6_generators.py --partial-only
.\\.competition_env\\Scripts\\python.exe -m unittest test_v6_submission -v
.\\.competition_env\\Scripts\\python.exe fit_v6_submission.py
.\\.competition_env\\Scripts\\python.exe validate_v6_submission.py --model-only
.\\.competition_env\\Scripts\\python.exe build_v6_submission.py --preflight --workers 4
.\\.competition_env\\Scripts\\python.exe build_v6_submission.py --workers 4
.\\.competition_env\\Scripts\\python.exe validate_v6_submission.py
.\\.competition_env\\Scripts\\python.exe document_v6_submission.py
```

Generated outputs refuse overwrite. Use a fresh experiment version for subsequent changes. Code, compact source/selection records, all source/control identities, verification records and figures are versioned; numerical NPZ caches and large H5AD/VCC files remain local.

## Submit

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v6_fidelity_calibrated\\prediction.vcc -m "Fidelity-guarded lineage response v6" --wait
```

This uploads a new entry. Preserve its ID and six metrics to compare against V4 and V5 on the same panel/anchor set. Original planning and frozen V4/V5 effect comparisons remain in [V6_MODEL_STRATEGY.md](V6_MODEL_STRATEGY.md).
'''
    (ROOT/'V6_SUBMISSION.md').write_text(text,encoding='utf-8')
    path=ROOT/'README.md';readme=path.read_text(encoding='utf-8');heading='## V6 fidelity-guarded submission\n';anchor='## Next-submission strategy after V5 and neural modeling\n'
    if heading in readme:
        start=readme.index(heading);stop=readme.index('\n## ',start+3);readme=readme[:start]+readme[stop+1:]
    entry=heading+f'\nCreated and independently verified `experiments/v6_fidelity_calibrated/prediction.vcc` ({ready["bytes"]:,} bytes), with all 900 groups, 360,000 cells and 18,533 genes. Frozen candidate: `{model["candidate"]}`, settings `{settings}`. A retains V4 effects; B/C updates are bounded and cannot reverse V4 signs or create previously zero effects. Neural templates were evaluated and withheld after distribution/perturbation-evidence checks.\n\n'
    entry+='Seven tests, disjoint source/control identity checks, all-coordinate model replay, the full count/library scan, every packaged count/index/pointer, all source-row provenance checks and 12 exact full-group replays passed. Official VCC packaging passed. [V6_SUBMISSION.md](V6_SUBMISSION.md) records actual source benchmarks, limitations, settings, hashes and Command Prompt reproduction/upload commands. The package has not been uploaded or scored; V4 remains the best published project result.\n\n'
    latest='The newest upload-ready candidate is V6, documented in [V6_SUBMISSION.md](V6_SUBMISSION.md). It has not been uploaded or scored.\n\n'
    if latest not in readme:
        paragraph_end=readme.index('\n\n',readme.index('\n\n')+2)
        readme=readme[:paragraph_end+2]+latest+readme[paragraph_end+2:]
    assert anchor in readme;path.write_text(readme.replace(anchor,entry+anchor,1),encoding='utf-8')
    code=['v6_methods.py','calibrate_v6_submission.py','evaluate_v6_neural_templates.py','audit_v6_generators.py','fit_v6_submission.py','build_v6_submission.py','test_v6_submission.py','validate_v6_submission.py','document_v6_submission.py']
    build.write_json(A/'documentation_manifest.json',{'root_code_sha256':{p:build.sha256(ROOT/p) for p in code},'document_sha256':build.sha256(ROOT/'V6_SUBMISSION.md'),
        'artifact_sha256':{p.name:build.sha256(p) for p in A.iterdir() if p.is_file() and p.suffix in ['.json','.csv','.txt','.png','.pdf'] and p.name!='documentation_manifest.json'},
        'scope':'Completed local generation and verification. No upload or new leaderboard score.'})
    print('V6 submission, measured validation, figures and README documented.',flush=True)


if __name__=='__main__':main()
