"""Replay the standalone inference CLI on real measured inputs for all contexts."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile
import anndata as ad
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from build_submission import ROOT, write_json
from train_gene300_models import OUT, load_model, predict


def main():
    for context in 'ABC':
        if not (OUT / f'model_{context}.npz').exists():
            raise FileNotFoundError(f'Train {context} model before running this integration check')
    with zipfile.ZipFile(ROOT / 'vcc_2026_300_target_matrix.zip') as archive:
        frame = pd.read_csv(archive.open(archive.namelist()[0]))
    small = frame.groupby('context', sort=False).head(3).copy()
    small['cell_id'] = [f'smoke_{c}_{i}' for i, c in enumerate(small.context)]
    report = {}
    with tempfile.TemporaryDirectory(prefix='gene300_inference_') as directory:
        folder = Path(directory)
        source = folder / 'measured300.csv'
        small[list(small.columns)[::-1]].to_csv(source, index=False)
        for context in 'ABC':
            destination = folder / f'{context}.h5ad'
            command = [sys.executable, str(ROOT / 'predict_gene300_models.py'), '--context', context,
                       '--input', str(source), '--output', str(destination)]
            subprocess.run(command, text=True, capture_output=True, check=True, cwd=ROOT)
            model = load_model(OUT / f'model_{context}.npz')
            measured = small.loc[small.context == context, list(model['input_genes'])].to_numpy(dtype=float)
            with threadpool_limits(limits=8):
                expected, _, _ = predict(model, measured)
            actual = ad.read_h5ad(destination)
            np.testing.assert_allclose(actual.X, expected, rtol=1e-6, atol=1e-6)
            np.testing.assert_array_equal(actual.var_names, model['output_genes'])
            np.testing.assert_array_equal(actual.obsm['measured_300_raw_counts'], measured)
            np.testing.assert_array_equal(actual.obs_names, small.loc[small.context == context, 'cell_id'])
            report[context] = {'passed': True, 'cells': 3, 'output_genes': actual.n_vars,
                              'named_input_columns_reversed': True, 'original_input_ids_preserved': True}
            print(context, 'standalone prediction CLI passed', flush=True)
    write_json(OUT / 'inference_smoke.json', report)


if __name__ == '__main__':
    main()
