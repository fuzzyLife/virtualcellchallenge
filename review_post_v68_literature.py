"""Pin and verify sources for POST_V68_RESEARCH.md; never execute remote code.

Default: read pinned public author files and write hashes/symbol inventories.
--check: verify those cached bytes and local evidence without network access.
--restore-cache: reacquire the exact recorded code files without changing the receipt.
"""
from __future__ import annotations
import argparse
import ast
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'analysis' / 'post_v68_research'
CACHE = OUT / 'cache'
RECEIPT = OUT / 'source_audit.json'
REPOS = {
    'whd1125/PopPert': ('248fc3e1428de96335bef1d0271d30ca1e4d5565', [
        'README.md', 'model.py', 'sampling.py', 'zinb/sampling.py']),
    'sanjukta7/sc-evae': ('2b367c2acd085478705e2934edfdcf5c6121833b', [
        'README.md', 'LICENSE', 'src/sc-evae/models/output_heads.py',
        'src/sc-evae/models/vae.py']),
    'hwxing3259/GPerturb': ('05d5a431413ff9dc443c82b994b571ecc563426d', [
        'README.md', 'GPerturb/GPerturb_model.py']),
}
LOCAL = [
    'v15_population_methods.py', 'v16_methods.py', 'v68_population_balancing.py',
    'V17_STRATEGY.md', 'V17_RESULTS.md', 'V27_RESULTS.md', 'V32_V34_RESULTS.md',
    'V52_JOINT_RESULTS.md', 'V66_LEADERBOARD_RESULTS.md', 'V67_RESULTS.md',
    'V68_RESULTS.md', 'experiments/v68_population_balancing/analysis/assessment.json',
]
PAPERS = [
    ('PopPert', 'https://arxiv.org/html/2609.01357v1', 'Full HTML and author code read; previously cited in V17.'),
    ('ExpressionVAE design space', 'https://www.biorxiv.org/content/10.64898/2026.06.15.732063v1.full',
     'Full page fetch failed; indexed primary passages, author page and decoder source inspected.'),
    ('GPerturb', 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12215016/', 'Full HTML methods and author code inspected.'),
    ('Perturbation-response score', 'https://www.nature.com/articles/s41556-025-01626-9',
     'Full HTML methods/data availability inspected; previously cited/tested in simplified form in V17.'),
    ('VCBench', 'https://arxiv.org/html/2604.27646v1', 'Full HTML framework and generalization discussion inspected.'),
]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def fetch(item):
    repo, rev, path = item
    url = f'https://raw.githubusercontent.com/{repo}/{rev}/{path}'
    req = urllib.request.Request(url, headers={'User-Agent': 'VCC-public-literature-audit'})
    with urllib.request.urlopen(req, timeout=30) as response:
        data = response.read(3_000_001)
    if len(data) > 3_000_000:
        raise ValueError(f'Unexpectedly large source: {url}')
    dest = CACHE / repo.replace('/', '__') / path
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    rec = dict(repository=repo, revision=rev, path=path, url=url,
               cache=dest.relative_to(ROOT).as_posix(), bytes=len(data), sha256=digest(data))
    if path.endswith('.py'):
        tree = ast.parse(data.decode('utf-8-sig'))
        rec['symbols'] = [dict(name=n.name, line=n.lineno, kind=type(n).__name__)
                          for n in ast.walk(tree)
                          if isinstance(n, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))]
    return rec


def verify(receipt):
    for item in receipt['local_evidence'] + receipt['remote_code']:
        path = ROOT / item.get('cache', item['path'])
        data = path.read_bytes()
        if digest(data) != item['sha256'] or len(data) != item['bytes']:
            raise ValueError(f'Evidence changed: {path}')
    d = json.loads((ROOT / 'experiments/v68_population_balancing/analysis/assessment.json').read_text())
    ratios = [d['studies']['balanced_real'][s]['profile_square_error'] /
              d['studies']['v63'][s]['profile_square_error'] for s in ('H1', 'Orion', 'CD4')]
    value = sum(ratios) / len(ratios)
    if abs(value - receipt['v68_equal_study_profile_ratio']) > 1e-12:
        raise ValueError('V68 paired study-ratio arithmetic mismatch')
    return value


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--restore-cache', action='store_true')
    args = ap.parse_args()
    if args.check or args.restore_cache:
        rec = json.loads(RECEIPT.read_text(encoding='utf-8'))
        if args.restore_cache:
            items = [(r['repository'], r['revision'], r['path']) for r in rec['remote_code']]
            with ThreadPoolExecutor(max_workers=4) as pool:
                restored = list(pool.map(fetch, items))
            if restored != rec['remote_code']:
                raise ValueError('Restored author source differs from the pinned audit')
    else:
        if RECEIPT.exists():
            raise FileExistsError('Receipt already exists; use --check to preserve the original audit.')
        items = [(r, rev, p) for r, (rev, paths) in REPOS.items() for p in paths]
        with ThreadPoolExecutor(max_workers=4) as pool:
            remote = list(pool.map(fetch, items))
        local = []
        for name in LOCAL:
            data = (ROOT / name).read_bytes()
            local.append(dict(path=name, bytes=len(data), sha256=digest(data)))
        d = json.loads((ROOT / LOCAL[-1]).read_text())
        ratio = sum(d['studies']['balanced_real'][s]['profile_square_error'] /
                    d['studies']['v63'][s]['profile_square_error'] for s in ('H1','Orion','CD4')) / 3
        rec = dict(review_date='2026-10-19', scope='Literature/code review, no model trained or submission generated',
                   papers=[dict(title=t, url=u, access=a) for t,u,a in PAPERS],
                   remote_code=remote, local_evidence=local, v68_equal_study_profile_ratio=ratio)
        OUT.mkdir(parents=True, exist_ok=True)
        RECEIPT.write_text(json.dumps(rec, indent=2) + '\n', encoding='utf-8')
    ratio = verify(rec)
    print(f"Verified {len(rec['remote_code'])} pinned author files, {len(rec['local_evidence'])} local evidence files.")
    print(f'V68 paired equal-study profile error / V63: {ratio:.12f}')
    print('Provenance/arithmetic checks only; no new model performance result.')


if __name__ == '__main__':
    main()
