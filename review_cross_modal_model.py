"""Audit local cross-modal paper/code/artifacts; no training or submission changes.

Run from the project root:
  .\\.competition_env\\Scripts\\python.exe review_cross_modal_model.py
  .\\.competition_env\\Scripts\\python.exe review_cross_modal_model.py --check-online

Large data files are inventoried by size, not read or hashed. Gene axes of small
mean-embedding files are inspected through HDF5; numeric embedding X is not read.
Only named pure benchmark helpers are compiled for synthetic regression probes.
"""
from pathlib import Path
import argparse, ast, hashlib, json, subprocess, warnings
from collections import Counter
from datetime import datetime, timezone
import h5py
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
CLONE = ROOT / 'cross_modal_foundation_model'
OUT = ROOT / 'analysis/cross_modal_review'


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def git(*args):
    result = subprocess.run(['git', '-c', 'safe.directory=' + CLONE.as_posix(),
                             '-C', str(CLONE), *args], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(result.stderr)
    return result.stdout.strip()


def benchmark_probes():
    path = CLONE / 'phase0/benchmark.py'
    tree = ast.parse(path.read_text(encoding='utf-8-sig'))
    names = {'_median_centre', '_cosine_knn_predict', '_pearson_safe',
             '_spearman_safe', 'compute_metrics'}
    helpers = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    assert len(helpers) == len(names)
    namespace = dict(np=np, warnings=warnings, pearsonr=pearsonr,
                     spearmanr=spearmanr, K_NEIGHBORS=5, TOP_N_DE=100)
    exec(compile(ast.Module(body=helpers, type_ignores=[]), str(path), 'exec'), namespace)
    centred = namespace['_median_centre'](np.array([1., np.nan, 3.]))
    prediction = namespace['_cosine_knn_predict'](
        np.array([1., 0.]), np.array([[1., 0.], [1., 0.]]),
        np.array([[10., np.nan], [10., 4.]]), k=2)
    short = namespace['compute_metrics'](np.zeros(3), np.zeros(3), np.zeros(3, bool))
    # A correctly conditioned mean for feature 2 would renormalize valid weight to 1.
    assert np.isfinite(centred).all() and centred[1] == 0
    assert np.allclose(prediction, [10., 2.])
    assert len(short) == 2
    return {
        'median_centre_input': [1, None, 3], 'median_centre_output': centred.tolist(),
        'missingness_lost_before_pairwise_complete_knn': True,
        'weighted_knn_prediction': prediction.tolist(),
        'valid_weight_renormalized_expected': [10., 4.],
        'metric_empty_mask_return_length': len(short), 'caller_expected_length': 4,
        'scope': 'Synthetic execution of exact named local helpers; not model efficacy.'}


def online_manifests():
    import requests
    endpoints = {
        'trove': 'https://huggingface.co/api/models/Tesorai/Trove1-70m',
        'author_code': 'https://api.github.com/repos/tesorai/cross_modal_foundation_model/contents/'}
    result = {'checked_at_utc': datetime.now(timezone.utc).isoformat()}
    for label, url in endpoints.items():
        try:
            response = requests.get(url, timeout=30)
            entry = {'url': url, 'http_status': response.status_code}
            if response.ok:
                data = response.json()
                if label == 'trove':
                    entry.update(commit=data.get('sha'), last_modified=data.get('lastModified'),
                                 files=[item['rfilename'] for item in data.get('siblings', [])])
                else:
                    entry['files'] = [item['name'] for item in data]
            result[label] = entry
        except requests.RequestException as exc:
            result[label] = {'url': url, 'error_type': type(exc).__name__}
    write_json(OUT / 'public_release_check.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-online', action='store_true')
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    paths = sorted(p for p in CLONE.rglob('*') if p.is_file() and '.git' not in p.parts)
    extensions = Counter(p.suffix.lower() or '(none)' for p in paths)
    checkpoints = [p.relative_to(CLONE).as_posix() for p in paths
                   if p.suffix.lower() in {'.pt', '.pth', '.ckpt', '.safetensors', '.bin'}]
    tracked = set(git('ls-files').splitlines())
    code = []
    for path in paths:
        if path.suffix != '.py':
            continue
        tree = ast.parse(path.read_text(encoding='utf-8-sig'))
        code.append({'path': path.relative_to(CLONE).as_posix(),
                     'tracked_by_clone_git': path.relative_to(CLONE).as_posix() in tracked, 'sha256': sha256(path),
                     'functions': '; '.join(f'{n.name}:{n.lineno}' for n in ast.walk(tree)
                                           if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)))})
    pd.DataFrame(code).to_csv(OUT / 'code_inventory.csv', index=False)
    genes = pd.read_csv(ROOT / 'gene_names.csv').gene_name.astype(str)
    targets = pd.read_csv(ROOT / 'pert_counts.csv').target_gene.astype(str)
    rows, mappings = [], []
    base = CLONE / 'phase1/data/phase1/benchmarks/depmap/gene-embs'
    for name in ['tx1_70m-mean-lt5gt70-bin.h5ad', 'tx1_1b-mean-lt5gt70-bin.h5ad',
                 'tx1_3b-mean-lt5gt70-bin.h5ad', 'pcloading-lt5gt70-bin.h5ad']:
        path = base / name
        if not path.exists():
            continue
        with h5py.File(path) as handle:
            axis = handle['obs/gene'].asstr()[:]
            shape = list(handle['X'].shape)
        counts = Counter(axis)
        positions = {g: i for i, g in enumerate(axis) if counts[g] == 1}
        rows.append({'file': path.relative_to(CLONE).as_posix(), 'bytes': path.stat().st_size,
                     'sha256': sha256(path), 'rows': shape[0], 'dimensions': shape[1],
                     'duplicate_gene_symbols': sum(n > 1 for n in counts.values()),
                     'competition_genes_covered': sum(g in positions for g in genes),
                     'competition_targets_covered': sum(t in positions for t in targets),
                     'model_type': 'RNA-only Tahoe-x1' if name.startswith('tx1') else 'PCA'})
        for target in targets:
            mappings.append({'embedding': name, 'target_gene': target,
                             'available_unique_symbol': target in positions,
                             'embedding_row': positions.get(target)})
    pd.DataFrame(rows).to_csv(OUT / 'embedding_coverage.csv', index=False)
    pd.DataFrame(mappings).to_csv(OUT / 'target_embedding_mapping.csv', index=False)
    results = CLONE / 'phase0/results/phase0_results.csv'
    if results.exists():
        table = pd.read_csv(results)
        summary = table.groupby(['dataset', 'baseline'])[
            ['pearson_full', 'pearson_top100', 'spearman_full', 'spearman_top100']].agg(['mean', 'std', 'count'])
        summary.to_csv(OUT / 'local_baseline_summary.csv')
    paper = CLONE / 'paper.pdf'
    reader = PdfReader(paper)
    fulltext = '\n\n'.join(f'PAGE {i+1}\n' + (page.extract_text() or '')
                           for i, page in enumerate(reader.pages))
    (OUT / 'paper_extracted.txt').write_text(fulltext, encoding='utf-8')
    write_json(OUT / 'code_probes.json', benchmark_probes())
    record = {'review_script_sha256': sha256(Path(__file__)),
              'clone_head': git('rev-parse', 'HEAD'), 'clone_remotes': git('remote', '-v'),
              'clone_tracked_status': git('status', '--porcelain', '--untracked-files=no'),
              'clone_commits': git('log', '-6', '--format=%H %ad %s', '--date=iso').splitlines(),
              'file_count_excluding_git': len(paths), 'bytes_excluding_git': sum(p.stat().st_size for p in paths),
              'extension_counts': dict(sorted(extensions.items())), 'checkpoint_candidates': checkpoints,
              'paper_sha256': sha256(paper), 'paper_pages': len(reader.pages),
              'gene_axis_sha256': sha256(ROOT / 'gene_names.csv'),
              'target_axis_sha256': sha256(ROOT / 'pert_counts.csv'),
              'competition_genes': len(genes), 'competition_targets': len(targets),
              'cached_benchmark_sha256': sha256(results) if results.exists() else None,
              'large_file_policy': 'Size inventory only; hashes restricted to code, paper, small mean embeddings and result table.'}
    write_json(OUT / 'review.json', record)
    if args.check_online:
        print(json.dumps(online_manifests(), indent=2))
    print(pd.DataFrame(rows)[['file', 'dimensions', 'competition_genes_covered', 'competition_targets_covered']].to_string(index=False))
    print('Review recorded at', OUT)


if __name__ == '__main__':
    main()
