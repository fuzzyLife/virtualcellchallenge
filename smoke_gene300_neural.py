"""Test the standalone neural CLI against real measured cells, both objectives."""
from pathlib import Path
import subprocess
import sys
import tempfile
from neural_runtime import activate_local_runtime
activate_local_runtime()
import anndata as ad
import numpy as np
import pandas as pd
import torch
from threadpoolctl import threadpool_limits
from build_submission import ROOT, read_spec, write_json
from train_gene300_neural import OUT, load_network, predict_network


def main():
    _, genes, targets = read_spec(ROOT)
    positions = [np.flatnonzero(genes == gene)[0] for gene in targets]
    frames = []
    for context in 'ABC':
        source = ad.read_h5ad(ROOT / f'context_{context}.h5ad', backed='r')
        inputs = source.X[:3][:, positions].toarray()
        frame = pd.DataFrame(inputs, columns=targets)
        frame['context'] = context
        frame['cell_id'] = source.obs_names[:3]
        frames.append(frame)
        source.file.close()
    measured = pd.concat(frames, ignore_index=True)
    report = {}
    with tempfile.TemporaryDirectory(prefix='gene300_neural_smoke_') as folder:
        folder = Path(folder)
        csv = folder / 'inputs.csv'
        measured[list(measured.columns)[::-1]].to_csv(csv, index=False)
        for context in 'ABC':
            for objective in ('normalized', 'raw'):
                destination = folder / f'{context}_{objective}.h5ad'
                command = [sys.executable, str(ROOT / 'predict_gene300_neural.py'), '--context', context,
                    '--objective', objective, '--input', str(csv), '--output', str(destination)]
                subprocess.run(command, capture_output=True, text=True, check=True, cwd=ROOT)
                model, network = load_network(OUT / f'model_{objective}_{context}.npz')
                inputs = measured.loc[measured.context == context, list(model['input_genes'])].to_numpy(dtype=float)
                expected = predict_network(model, network, inputs)
                actual = ad.read_h5ad(destination)
                np.testing.assert_allclose(actual.X, expected, rtol=1e-5, atol=1e-5)
                np.testing.assert_array_equal(actual.var_names, model['output_genes'])
                np.testing.assert_array_equal(actual.obs_names, measured.loc[measured.context == context, 'cell_id'])
                np.testing.assert_array_equal(actual.obsm['measured_300_raw_counts'], inputs)
                assert actual.uns['objective'] == objective
                report[f'{context}_{objective}'] = {'passed': True, 'cells': 3, 'output_genes': 18233,
                    'reversed_named_input_columns': True, 'original_cell_ids_preserved': True}
                print(context, objective, 'standalone CLI passed', flush=True)
    write_json(OUT / 'inference_smoke.json', report)


if __name__ == '__main__':
    torch.set_num_threads(8)
    with threadpool_limits(limits=8):
        main()
