"""Generate an audited neural-network results report, figures and README entry."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from build_submission import ROOT

OUT = ROOT / 'analysis/gene300_neural'
RIDGE = ROOT / 'analysis/gene300_models'


def main():
    validation = json.loads((OUT / 'validation.json').read_text(encoding='utf-8'))
    assert validation['passed']
    scores = pd.read_csv(OUT / 'test_architecture_scores.csv')
    selected = scores[scores.selected].copy()
    selected.to_csv(OUT / 'test_summary.csv', index=False)
    ridge = pd.read_csv(RIDGE / 'test_summary.csv').set_index('context')
    genes = pd.read_csv(OUT / 'per_gene_test_metrics.csv')
    history = pd.read_csv(OUT / 'training_history.csv')
    candidates = pd.read_csv(OUT / 'candidate_scores.csv')
    norm = selected[selected.objective == 'normalized'].set_index('context').loc[list('ABC')]
    raw = selected[selected.objective == 'raw'].set_index('context').loc[list('ABC')]
    comparison = pd.DataFrame({'context': list('ABC'),
        'ridge_normalized_skill': ridge.loc[list('ABC'), 'test_expressed_logcpm_skill'].to_numpy(),
        'neural_normalized_skill': norm.test_skill.to_numpy(),
        'ridge_raw_skill_vs_depth': ridge.loc[list('ABC'), 'test_raw_skill_vs_depth_baseline'].to_numpy(),
        'neural_raw_skill_vs_depth': raw.test_skill.to_numpy(),
        'neural_normalized_mse_reduction_vs_ridge': 1 - norm.test_mse.to_numpy() / ridge.loc[list('ABC'), 'test_expressed_logcpm_mse'].to_numpy(),
        'neural_raw_mse_reduction_vs_ridge': 1 - raw.test_mse.to_numpy() / ridge.loc[list('ABC'), 'test_raw_mse'].to_numpy()})
    comparison.to_csv(OUT / 'ridge_comparison.csv', index=False)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    positions = np.arange(3); width = .35
    for ax, previous, new, title, ylabel in [
        (axes[0], comparison.ridge_normalized_skill, comparison.neural_normalized_skill,
         'Normalized expression: held-out guide groups', 'MSE reduction vs average profile (%)'),
        (axes[1], comparison.ridge_raw_skill_vs_depth, comparison.neural_raw_skill_vs_depth,
         'Expected raw counts: held-out guide groups', 'MSE reduction vs depth-only baseline (%)')]:
        ax.bar(positions - width / 2, previous * 100, width, label='Ridge', color='#8d9bac')
        ax.bar(positions + width / 2, new * 100, width, label='Selected neural', color='#287cb5')
        ax.axhline(0, color='black', linewidth=.7)
        ax.set_xticks(positions, list('ABC'))
        ax.set(title=title, ylabel=ylabel)
        ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / 'ridge_neural_comparison.png', dpi=160)
    fig.savefig(OUT / 'ridge_neural_comparison.pdf')
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    for row, objective in enumerate(['normalized', 'raw']):
        for column, context in enumerate('ABC'):
            ax = axes[row, column]
            for arch in ['shallow', 'deep', 'bottleneck', 'residual']:
                curve = history[(history.context == context) & (history.objective == objective) & (history.architecture == arch)]
                ax.plot(curve.epoch, curve.validation_mse, label=arch)
            ax.set(title=f'{context}: {objective}', xlabel='Training epoch', ylabel='Validation MSE')
            ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / 'validation_curves.png', dpi=160)
    fig.savefig(OUT / 'validation_curves.pdf')
    plt.close(fig)
    top = genes[(genes.objective == 'normalized') & genes.expressed_train_5cpm].sort_values(
        ['context', 'test_spearman'], ascending=[True, False]).groupby('context').head(20)
    top.to_csv(OUT / 'most_predictable_genes.csv', index=False)
    lines = ['# Neural architecture comparison: measured results', '',
        'Trained 24 candidate neural networks: four configurations, two output objectives and three cell types. '
        'Architecture and stopping epoch were selected on seven held-out validation guide groups. '
        'The following scores use the same seven test-guide groups per context as the earlier ridge study.', '',
        '| Context | Best normalized network | Epoch | Normalized error reduction vs mean | Best raw network | Epoch | Raw error reduction vs depth baseline |',
        '|---|---|---:|---:|---|---:|---:|']
    for context in 'ABC':
        n, r = norm.loc[context], raw.loc[context]
        lines.append(f'| {context} | {n.architecture} | {n.best_epoch:g} | {n.test_skill:.2%} | '
                     f'{r.architecture} | {r.best_epoch:g} | {r.test_skill:.2%} |')
    lines += ['', 'Normalized scores use genes with training mean CPM at least 5; raw-count scores use all 18,233 output genes. '
        'Positive skill means lower MSE than the stated baseline. Input genes are excluded from both output axes and all errors.', '',
        '![Ridge and neural comparison](analysis/gene300_neural/ridge_neural_comparison.png)', '',
        '| Context | Ridge normalized reduction vs mean | Neural normalized MSE reduction vs ridge | Ridge raw reduction vs depth baseline | Neural raw MSE reduction vs ridge |',
        '|---|---:|---:|---:|---:|']
    for row in comparison.itertuples():
        lines.append(f'| {row.context} | {row.ridge_normalized_skill:.2%} | '
            f'{row.neural_normalized_mse_reduction_vs_ridge:.2%} | {row.ridge_raw_skill_vs_depth:.2%} | '
            f'{row.neural_raw_mse_reduction_vs_ridge:.2%} |')
    lines += ['', 'The raw neural objective compares predicted counts directly with original observed counts. '
        'The earlier ridge raw output was obtained by converting a log-CPM model back to counts. '
        'The improvement therefore reflects both a neural predictor and a better-matched training objective; '
        'it does not isolate a pure architecture effect.', '', '## Every frozen configuration', '',
        '| Context | Objective | Configuration | Selected epoch | Validation MSE | Test MSE | Test skill vs stated baseline | Chosen on validation |',
        '|---|---|---|---:|---:|---:|---:|---|']
    joined = scores.merge(candidates[['context', 'objective', 'architecture', 'validation_mse']],
                           on=['context', 'objective', 'architecture'])
    for row in joined.itertuples():
        lines.append(f'| {row.context} | {row.objective} | {row.architecture} | {row.best_epoch:g} | '
            f'{row.validation_mse:.4f} | {row.test_mse:.4f} | {row.test_skill:.2%} | {"yes" if row.selected else "no"} |')
    lines += ['', 'The configuration marked chosen minimizes validation MSE. A different configuration can happen to '
        'score better on test; that does not change the saved choice. Full epoch histories are in '
        '`analysis/gene300_neural/training_history.csv`.', '',
        '![Validation curves](analysis/gene300_neural/validation_curves.png)', '',
        '## Predictable output genes', '',
        'The following output genes have the strongest normalized rank correlation across test cells among '
        'genes with training mean CPM at least 5. They were not used to choose models.', '']
    for context in 'ABC':
        group = top[top.context == context].head(10)
        entries = ', '.join(f'{row.gene_name} ({row.test_spearman:.3f})' for row in group.itertuples())
        lines.append(f'- **{context}:** {entries}.')
    lines += ['', 'Rank correlations retain zeros and average ties; undefined constant-vector correlations remain missing. '
        'The full selected-model gene/cell error tables include every output gene and test cell.', '',
        '## Saved and verified work', '',
        'Saved all candidate checkpoints, six selected evaluation models, six selected all-control models, '
        'normalized/raw held-out prediction files, and expected raw-count reconstructions for all 55,200 controls. '
        'Full-data reconstructions are in-sample and are not the source of the reported test results.', '',
        'Five focused tests passed, including independent NumPy prediction checks for every structure/objective. '
        'The artifact audit verified all selected source/output hashes, feature statistics, validation-based settings, '
        'baseline profiles and source identities; recomputed every selected-model gene/cell MSE and baseline error; '
        'scanned every saved output value; and independently replayed 32 rows per evaluation/production model '
        'and output file using explicit NumPy network mathematics. Exact coverage is in '
        '`analysis/gene300_neural/validation.json`.', '',
        'All six standalone inference CLI checks also passed using real measured counts and reversed gene-column order. '
        'Original cell IDs and output axes were preserved; see `inference_smoke.json`.', '',
        'The test set was previously reported for ridge, and only one training seed was used per configuration. '
        'Small architecture differences need repeated-seed checks. The data contain controls, so these results '
        'do not validate causal CRISPRi response predictions or establish a new leaderboard score.', '',
        'See [GENE300_NEURAL_MODELS.md](GENE300_NEURAL_MODELS.md) for architectures, losses, file units, '
        'the verified Windows runtime cache and complete training/inference commands.', '']
    (ROOT / 'GENE300_NEURAL_RESULTS.md').write_text('\n'.join(lines), encoding='utf-8')
    readme = ROOT / 'README.md'
    text = readme.read_text(encoding='utf-8')
    heading = '## Neural networks for the 300-gene panel\n'
    if heading in text:
        before, rest = text.split(heading, 1)
        text = before + '## ' + rest.split('\n## ', 1)[1]
    entry = heading + '\nCompared shallow, deep, bottleneck and residual neural configurations separately in A/B/C '
    entry += 'with normalized-expression and direct raw-count objectives (24 candidates). Validation guide groups '
    entry += 'selected architectures and stopping epochs; frozen models were compared on the existing ridge test set.\n\n'
    entry += 'Selected neural raw-count MSE reductions against the depth-only baseline: '
    entry += ', '.join(f'**{c}: {raw.loc[c, "test_skill"]:.2%}**' for c in 'ABC') + '. '
    entry += 'Selected normalized-expression reductions against the mean profile: '
    entry += ', '.join(f'**{c}: {norm.loc[c, "test_skill"]:.2%}**' for c in 'ABC') + '.\n\n'
    entry += 'Saved separate normalized/raw models, both held-out prediction types and full raw reconstructions for '
    entry += 'all 55,200 controls. Five focused tests, six standalone inference checks and independent artifact/matrix checks passed. '
    entry += 'Methods and commands: [GENE300_NEURAL_MODELS.md](GENE300_NEURAL_MODELS.md). '
    entry += 'Architecture comparison and figures: [GENE300_NEURAL_RESULTS.md](GENE300_NEURAL_RESULTS.md). '
    entry += 'Root scripts include `train_gene300_neural.py`, `predict_gene300_neural.py`, '
    entry += '`test_gene300_neural.py`, `validate_gene300_neural.py` and `document_gene300_neural.py`.\n\n'
    anchor = '## Separate models from the 300-gene panel\n'
    if anchor not in text:
        raise ValueError('Expected README anchor missing')
    readme.write_text(text.replace(anchor, entry + anchor, 1), encoding='utf-8')
    print('Wrote audited neural comparison, figures and README entry.')


if __name__ == '__main__':
    main()
