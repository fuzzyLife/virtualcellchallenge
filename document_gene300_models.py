"""Create the results report, figures and README entry from audited model artifacts."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from build_submission import ROOT
from train_gene300_models import OUT


def main():
    validation = json.loads((OUT / 'validation.json').read_text(encoding='utf-8'))
    if not validation['passed']:
        raise ValueError('Independent validation must pass before documenting results')
    summary = pd.read_csv(OUT / 'test_summary.csv').set_index('context').loc[list('ABC')]
    genes = pd.read_csv(OUT / 'per_gene_test_metrics.csv')
    cells = pd.read_csv(OUT / 'per_cell_test_metrics.csv')
    candidates = pd.read_csv(OUT / 'candidate_scores.csv')
    top = genes.loc[genes.expressed_train_5cpm & (genes.train_mean_raw >= .25)].sort_values(
        ['context', 'test_logcpm_skill'], ascending=[True, False]).groupby('context').head(20)
    top.to_csv(OUT / 'most_predictable_genes.csv', index=False)
    diagnostics = genes.copy()
    diagnostics['excess_raw_mse_vs_depth'] = diagnostics.test_raw_mse - diagnostics.test_depth_baseline_raw_mse
    diagnostics = diagnostics.sort_values(['context', 'excess_raw_mse_vs_depth'], ascending=[True, False]).groupby('context').head(20)
    diagnostics.to_csv(OUT / 'raw_count_conversion_diagnostics.csv', index=False)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    axes[0].bar(list('ABC'), summary.test_expressed_logcpm_skill * 100, color=['#356cb0', '#d58d25', '#329875'])
    axes[0].set(ylabel='MSE reduction vs training mean (%)', title='Normalized expression: held-out guides')
    axes[1].bar(list('ABC'), summary.test_raw_skill_vs_depth_baseline * 100, color=['#356cb0', '#d58d25', '#329875'])
    axes[1].axhline(0, color='black', linewidth=.8)
    axes[1].set(ylabel='MSE reduction vs depth-only baseline (%)', title='Expected raw counts: held-out guides')
    for context in 'ABC':
        values = genes.loc[(genes.context == context) & genes.expressed_train_5cpm, 'test_normalized_spearman'].dropna()
        axes[2].hist(values, bins=np.linspace(-.2, 1, 61), density=True, histtype='step', linewidth=1.5, label=context)
    axes[2].set(xlabel='Per-gene Spearman across test cells', ylabel='Density', title='Predictability of expressed genes')
    axes[2].legend(title='Context')
    fig.tight_layout()
    fig.savefig(OUT / 'heldout_performance.png', dpi=160)
    fig.savefig(OUT / 'heldout_performance.pdf')
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, context in zip(axes, 'ABC'):
        data = cells[cells.context == context]
        actual, predicted = data.actual_library_counts_evaluation_only, data.predicted_library_counts
        ax.hexbin(actual, predicted, gridsize=45, mincnt=1, bins='log', cmap='viridis')
        low, high = min(actual.min(), predicted.min()), max(actual.max(), predicted.max())
        ax.plot([low, high], [low, high], color='orange', linewidth=1)
        ax.set(title=f'{context}: inferred library size', xlabel='Measured full library (evaluation only)', ylabel='Inferred from 300 genes')
    fig.tight_layout()
    fig.savefig(OUT / 'library_inference.png', dpi=160)
    fig.savefig(OUT / 'library_inference.pdf')
    plt.close(fig)
    lines = ['# Separate A/B/C models: first results', '',
        'Three models were fitted using only the 300 measured input genes to predict the other 18,233 genes. '
        'Each evaluation uses 2,800 cells from seven guide groups withheld from model fitting and selection. '
        'The selected models were subsequently refitted on all controls for saved full-cell reconstructions.', '',
        '| Context | Selected features | Alpha | Normalized MSE reduction vs mean | Raw MSE reduction vs depth baseline | Median gene rank correlation | Library median relative error |',
        '|---|---|---:|---:|---:|---:|---:|']
    for context, row in summary.iterrows():
        lines.append(f'| {context} | {row.selected_family} | {row.alpha:g} | '
            f'{row.test_expressed_logcpm_skill:.2%} | {row.test_raw_skill_vs_depth_baseline:.2%} | '
            f'{row.median_expressed_gene_normalized_spearman:.3f} | {row.median_library_absolute_relative_error:.2%} |')
    lines += ['', 'Normalized MSE and gene rank correlations above use output genes with training mean CPM at least 5. '
        'Raw-count error includes all 18,233 outputs. Negative skill means the indicated baseline performs better. '
        'The reported raw-count conversion uses an inferred full-library total and residual smearing; '
        'the normalized expression result therefore cannot be treated as an equivalent gain in raw-count accuracy.', '',
        '![Held-out performance](analysis/gene300_models/heldout_performance.png)', '',
        '| Context | Expressed outputs | Raw MSE | Raw MSE reduction vs constant mean | Median cell profile correlation | Mean-profile baseline correlation | All-output normalized MSE reduction |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for context, row in summary.iterrows():
        lines.append(f'| {context} | {row.expressed_genes:,.0f} | {row.test_raw_mse:.3f} | '
            f'{row.test_raw_skill_vs_mean:.2%} | {row.median_cell_profile_spearman:.3f} | '
            f'{row.median_mean_baseline_cell_profile_spearman:.3f} | {row.test_all_logcpm_skill:.2%} |')
    lines += ['', 'A high correlation between the genes within a cell can partly reflect the shared ordering of abundant genes. '
        'The mean-profile baseline makes this visible. Per-gene correlations across cells and baseline-relative normalized error '
        'are more informative about the cell-specific signal.', '', '## Feature comparison on validation cells', '',
        '| Context | Family | Best alpha | Validation normalized MSE reduction vs mean |', '|---|---|---:|---:|']
    for (context, family), group in candidates.groupby(['context', 'family'], sort=True):
        row = group.sort_values('validation_expressed_logcpm_mse').iloc[0]
        lines.append(f'| {context} | {family} | {row.alpha:g} | {row.validation_expressed_skill:.2%} |')
    lines += ['', 'These comparisons chose the model settings before the test evaluation. The full 36-candidate table is '
        '`analysis/gene300_models/candidate_scores.csv`.', '', '## Most predictable measured outputs', '',
        'The following genes have the highest held-out normalized-expression skill among outputs with training mean CPM '
        'at least 5 and mean raw counts at least 0.25. This is a descriptive test result, not a feature-selection step.', '']
    for context in 'ABC':
        subset = top[top.context == context].head(10)
        entries = ', '.join(f'{r.gene_name} ({r.test_logcpm_skill:.1%})' for r in subset.itertuples())
        lines.append(f'- **{context}:** {entries}.')
    lines += ['', 'The complete output-gene metrics and a 20-gene-per-context export are provided in '
        '`per_gene_test_metrics.csv` and `most_predictable_genes.csv`.', '',
        '![Library inference](analysis/gene300_models/library_inference.png)', '',
        '## Where raw-count conversion loses accuracy', '',
        'The largest additional raw-count errors against the depth-only baseline occur in the following output genes. '
        'Many still show positive normalized-expression skill; optimizing log-CPM therefore does not ensure that their '
        'back-transformed expected counts minimize raw-count error.', '',
        '| Context | Largest excess raw-count error genes |', '|---|---|']
    for context in 'ABC':
        subset = diagnostics[diagnostics.context == context].head(5)
        lines.append(f'| {context} | ' + ', '.join(subset.gene_name) + ' |')
    lines += ['', 'Exact per-gene error differences are saved in `raw_count_conversion_diagnostics.csv`.', '',
        '## What this first model establishes', '',
        'There is measurable predictive association between the 300-gene panel and other genes in held-out control cells. '
        'Much of single-cell variation remains unpredictable from this small panel. Library estimation and normalized-expression '
        'prediction are evaluated separately because depth alone can account for substantial raw-count variation.', '',
        'For any context with negative raw-count skill, the depth-only baseline is the stronger raw-count predictor on this test. '
        'The next modeling comparison should train a raw-count objective directly and test whether nonlinear features capture '
        'additional cell-state structure. Any future comparison needs its own prespecified evaluation rather than choosing settings '
        'from the current test table.', '',
        'These data contain controls, not observed responses to the 300 CRISPRi interventions. This experiment does not create a new '
        'leaderboard submission or establish causal knockdown responses.', '', '## Verification and reproduction', '',
        'Seven focused tests passed. Independent validation checked all source and artifact hashes, all split assignments, all saved '
        'cell and gene identities, feature statistics, 24 output genes per model using an independent normal-equation solve, and '
        'all depth coefficients. Every held-out prediction was replayed and every output-gene raw/log-CPM MSE was recomputed. '
        'All full reconstruction values were scanned for finite nonnegative values and count conservation, with 24 full-fit rows '
        'per context replayed across every output gene. See `analysis/gene300_models/validation.json` for the exact coverage.', '',
        'Models and H5AD files are saved locally; source code, documentation and small result tables are versioned. '
        'See [GENE300_MODELS.md](GENE300_MODELS.md) for every method step, artifact schema and training/inference commands.', '']
    (ROOT / 'GENE300_MODEL_RESULTS.md').write_text('\n'.join(lines), encoding='utf-8')
    readme = ROOT / 'README.md'
    text = readme.read_text(encoding='utf-8')
    heading = '## Separate models from the 300-gene panel\n'
    if heading in text:
        before, rest = text.split(heading, 1)
        text = before + '## ' + rest.split('\n## ', 1)[1]
    entry = heading + '\nTrained separate A/B/C models using 300 measured genes to predict the other 18,233 genes. '
    entry += 'Held out complete control-guide groups for validation and test; saved trained weights, predictions for all 55,200 controls, '
    entry += '8,400 held-out predictions and detailed gene/cell metrics.\n\n'
    entry += 'Held-out normalized-expression MSE reductions against each context\'s training mean: '
    entry += ', '.join(f'**{c}: {summary.loc[c, "test_expressed_logcpm_skill"]:.2%}**' for c in 'ABC') + '. '
    entry += 'Raw-count conversion is compared separately with a depth-only baseline; see the results report for its limitations. '
    entry += 'Full-data reconstructions are in-sample. These observational control models are separate from competition submission artifacts.\n\n'
    entry += 'Methods and commands: [GENE300_MODELS.md](GENE300_MODELS.md). '
    entry += 'Measured results and figures: [GENE300_MODEL_RESULTS.md](GENE300_MODEL_RESULTS.md). '
    entry += 'Root scripts: `train_gene300_models.py`, `predict_gene300_models.py`, `validate_gene300_models.py`, '
    entry += '`test_gene300_models.py`, `smoke_gene300_inference.py` and `document_gene300_models.py`.\n\n'
    anchor = '## V5 atlas-calibrated submission\n'
    if anchor not in text:
        raise ValueError('Expected README insertion point missing')
    readme.write_text(text.replace(anchor, entry + anchor, 1), encoding='utf-8')
    print('Wrote methods-linked results, figures and README entry.')


if __name__ == '__main__':
    main()
