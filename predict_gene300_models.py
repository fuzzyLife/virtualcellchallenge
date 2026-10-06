"""Predict other genes from a CSV or ZIP containing the named 300-gene panel."""
import argparse
import csv
import io
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from train_gene300_models import OUT, load_model, write_predictions


def read_inputs(path, model, context):
    if path.suffix.lower() == '.zip':
        with zipfile.ZipFile(path) as archive:
            members = [m for m in archive.namelist() if not m.endswith('/')]
            if len(members) != 1 or not members[0].lower().endswith('.csv'):
                raise ValueError('ZIP must contain exactly one CSV')
            with archive.open(members[0]) as handle:
                text = handle.read().decode('utf-8-sig')
    else:
        text = path.read_text(encoding='utf-8-sig')
    header = next(csv.reader(io.StringIO(text)))
    if len(header) != len(set(header)):
        raise ValueError('Duplicate column names')
    frame = pd.read_csv(io.StringIO(text))
    missing = set(model['input_genes']) - set(frame.columns)
    if missing:
        raise ValueError(f'Missing required measured genes: {sorted(missing)}')
    unexpected = set(frame.columns) - set(model['input_genes']) - {'context', 'cell_id', 'ntc_id'}
    if unexpected:
        raise ValueError(f'Unexpected columns: {sorted(unexpected)}')
    if 'context' in frame:
        if not set(frame.context.astype(str)).issubset({'A', 'B', 'C'}):
            raise ValueError('Invalid context labels')
        frame = frame.loc[frame.context.astype(str) == context].copy()
    if not len(frame):
        raise ValueError('No rows for requested context')
    inputs = frame[list(model['input_genes'])].to_numpy(dtype=np.float64)
    if not np.isfinite(inputs).all() or (inputs < 0).any() or (inputs != np.floor(inputs)).any():
        raise ValueError('Measured genes must be finite nonnegative integer counts')
    obs = frame[[c for c in ('cell_id', 'ntc_id') if c in frame]].copy()
    obs['context'] = context
    obs.index = (frame.cell_id.astype(str).to_numpy() if 'cell_id' in frame
                 else np.array([f'{context}_input_row_{i}' for i in frame.index]))
    if not obs.index.is_unique:
        raise ValueError('Duplicate cell IDs')
    return inputs, obs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--context', required=True, choices=['A', 'B', 'C'])
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--model-dir', type=Path, default=OUT)
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f'Refusing to overwrite {args.output}')
    model = load_model(args.model_dir / f'model_{args.context}.npz')
    if str(model['context']) != args.context:
        raise ValueError('Model context mismatch')
    inputs, obs = read_inputs(args.input, model, args.context)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=args.threads):
        write_predictions(args.output, model, inputs, obs, 'new_input_all_controls_model')
    print(f'Saved {len(inputs):,} cells x {len(model["output_genes"]):,} other genes to {args.output}')
