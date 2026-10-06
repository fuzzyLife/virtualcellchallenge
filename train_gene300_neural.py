"""Compare four neural structures for 300-gene reconstruction in A/B/C.
Separate normalized-expression and directly optimized raw-count objectives.
"""
from __future__ import annotations
import argparse
import gc
from importlib.metadata import version
import json
from pathlib import Path
import shutil
import time
import anndata as ad
import numpy as np
import pandas as pd
from neural_runtime import activate_local_runtime
activate_local_runtime()
import torch
from torch import nn
from threadpoolctl import threadpool_limits
from build_submission import ROOT, load_control, read_spec, sha256, write_json
from train_gene300_models import features, standardize, split_guides, load_model as load_ridge, spearman

OUT = ROOT / 'analysis/gene300_neural'
RIDGE = ROOT / 'analysis/gene300_models'
SEED = 20261006
STRUCTURES = {
    'shallow': {'hidden': [128], 'activation': 'relu', 'dropout': .1, 'family': 'magnitude'},
    'deep': {'hidden': [256, 128], 'activation': 'silu', 'dropout': .1, 'family': 'magnitude'},
    'bottleneck': {'hidden': [128, 32], 'activation': 'silu', 'dropout': .05, 'family': 'magnitude'},
    'residual': {'hidden': [256, 128], 'activation': 'gelu', 'dropout': .1, 'family': 'hybrid'},
}


def log(message):
    print(time.strftime('%H:%M:%S') + ' ' + message, flush=True)


class GeneNetwork(nn.Module):
    def __init__(self, architecture, input_dim, outputs):
        super().__init__()
        spec = STRUCTURES[architecture]
        activation = {'relu': nn.ReLU, 'silu': nn.SiLU, 'gelu': nn.GELU}[spec['activation']]
        layers = []
        width = input_dim
        for hidden in spec['hidden']:
            layers.extend([nn.Linear(width, hidden), activation(), nn.Dropout(spec['dropout'])])
            width = hidden
        self.body = nn.Sequential(*layers)
        self.head = nn.Linear(width, outputs)
        nn.init.zeros_(self.head.weight); nn.init.zeros_(self.head.bias)
        self.branch = None
        if architecture == 'residual':
            self.branch = nn.Sequential(nn.Linear(input_dim, 32, bias=False), nn.Linear(32, outputs, bias=False))
            nn.init.zeros_(self.branch[-1].weight)

    def forward(self, x):
        y = self.head(self.body(x))
        return y if self.branch is None else y + self.branch(x)


def output_from_delta(delta, base, objective, clipped=False):
    if objective == 'raw':
        return base * torch.exp(torch.clamp(delta.float(), -3., 3.))
    y = base + delta.float()
    return torch.clamp(y, 0., np.log1p(1e6)) if clipped else y


def metadata(inputs, raw, labels, libraries, fit, context, architecture, objective, scope, seed):
    family = STRUCTURES[architecture]['family']
    _, mean, scale = standardize(features(inputs[fit], family))
    result = {'schema': np.array(1), 'context': np.array(context), 'architecture': np.array(architecture),
        'objective': np.array(objective), 'scope': np.array(scope), 'seed': np.array(seed),
        'family': np.array(family), 'feature_mean': mean, 'feature_scale': scale,
        'fit_row_indices': np.flatnonzero(fit), 'input_genes': np.array([], dtype=str), 'output_genes': np.array([], dtype=str)}
    if objective == 'raw':
        result['base_vector'] = np.asarray(raw[fit].sum(axis=0)).ravel().astype(np.float64) / inputs[fit].sum()
    else:
        result['base_vector'] = labels[fit].mean(axis=0, dtype=np.float64)
        if architecture == 'residual':
            prefix = 'evaluation_model' if scope == 'train_only_evaluation' else 'model'
            ridge = load_ridge(RIDGE / f'{prefix}_{context}.npz')
            np.testing.assert_array_equal(ridge['fit_row_indices'], result['fit_row_indices'])
            for key in ['family', 'feature_mean', 'feature_scale', 'coefficient', 'intercept']:
                result['ridge_' + key] = ridge[key]
    return result


def baseline(model, inputs):
    if str(model['objective']) == 'raw':
        return (inputs.sum(axis=1)[:, None] * model['base_vector']).astype(np.float32)
    if 'ridge_coefficient' in model:
        f = features(inputs, str(model['ridge_family']))
        x = (f - model['ridge_feature_mean']) / model['ridge_feature_scale']
        return (x @ model['ridge_coefficient'] + model['ridge_intercept']).astype(np.float32)
    return np.broadcast_to(model['base_vector'].astype(np.float32), (len(inputs), len(model['base_vector'])))


def save_network(path, model, network):
    arrays = dict(model)
    arrays.update({'weight__' + key: value.detach().cpu().numpy() for key, value in network.state_dict().items()})
    np.savez_compressed(path, **arrays)


def load_network(path, device='cpu'):
    with np.load(path, allow_pickle=False) as archive:
        model = {key: archive[key] for key in archive.files if not key.startswith('weight__')}
        weights = {key[8:]: torch.from_numpy(archive[key].copy()) for key in archive.files if key.startswith('weight__')}
    if int(model['schema']) != 1:
        raise ValueError('Unsupported neural model schema')
    network = GeneNetwork(str(model['architecture']), len(model['feature_mean']), len(model['output_genes']))
    network.load_state_dict(weights, strict=True)
    return model, network.to(device).eval()


def predict_network(model, network, inputs, batch=256):
    inputs = np.asarray(inputs, dtype=np.float64)
    x = (features(inputs, str(model['family'])) - model['feature_mean']) / model['feature_scale']
    device = next(network.parameters()).device
    result = np.empty((len(inputs), len(model['output_genes'])), dtype=np.float32)
    network.eval()
    with torch.no_grad():
        for start in range(0, len(inputs), batch):
            stop = min(start + batch, len(inputs))
            xb = torch.from_numpy(x[start:stop].astype(np.float32)).to(device)
            base = torch.from_numpy(np.array(baseline(model, inputs[start:stop]), copy=True)).to(device)
            prediction = output_from_delta(network(xb), base, str(model['objective']), clipped=True)
            result[start:stop] = prediction.cpu().numpy()
    if not np.isfinite(result).all() or (result < 0).any():
        raise ValueError('Nonfinite or negative neural prediction')
    return result


def train_network(model, inputs, labels, fit, val, expressed, args, context, architecture, objective, production_epochs=None):
    seed = int(model['seed'])
    torch.manual_seed(seed); np.random.seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    f = features(inputs, str(model['family']))
    x = torch.from_numpy(((f - model['feature_mean']) / model['feature_scale']).astype(np.float32)).to(args.device)
    y = torch.from_numpy(labels).to(args.device)
    # Baselines depend only on fitted parameters and the measured inputs, never on unseen labels.
    base = torch.from_numpy(np.array(baseline(model, inputs), copy=True)).to(args.device)
    rows = torch.from_numpy(np.flatnonzero(fit)).to(args.device)
    vrows = torch.from_numpy(np.flatnonzero(val)).to(args.device)
    mask = torch.from_numpy(expressed if objective == 'normalized' else np.ones(labels.shape[1], dtype=bool)).to(args.device)
    network = GeneNetwork(architecture, x.shape[1], labels.shape[1]).to(args.device)
    optimizer = torch.optim.AdamW(network.parameters(), lr=args.learning_rate, weight_decay=1e-4)
    scaler = torch.amp.GradScaler('cuda', enabled=args.device == 'cuda')
    loss_scale = float(np.mean(labels[fit] ** 2, dtype=np.float64)) if objective == 'raw' else 1.
    loss_scale = max(loss_scale, 1.)
    history = []

    def validation():
        network.eval()
        total = 0.
        with torch.no_grad():
            for start in range(0, len(vrows), args.batch_size):
                index = vrows[start:start + args.batch_size]
                pred = output_from_delta(network(x[index]), base[index], objective, clipped=True)
                total += ((pred[:, mask] - y[index][:, mask]) ** 2).sum().item()
        return total / (len(vrows) * int(mask.sum()))

    selecting = production_epochs is None
    if selecting:
        best = validation()
        best_epoch = 0
        best_state = {k: v.detach().cpu().clone() for k, v in network.state_dict().items()}
        history.append({'context': context, 'architecture': architecture, 'objective': objective,
                        'epoch': 0, 'train_mse': None, 'validation_mse': best, 'learning_rate': args.learning_rate})
    epochs = args.epochs if selecting else production_epochs
    for epoch in range(1, epochs + 1):
        network.train()
        permuted = rows[torch.randperm(len(rows), device=args.device)]
        loss_sum = 0.
        for start in range(0, len(permuted), args.batch_size):
            index = permuted[start:start + args.batch_size]
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type=args.device, dtype=torch.float16, enabled=args.device == 'cuda'):
                delta = network(x[index])
            pred = output_from_delta(delta, base[index], objective)
            loss = torch.mean((pred - y[index]) ** 2) / loss_scale
            if not torch.isfinite(loss):
                raise ValueError(f'Nonfinite training loss: {context}/{architecture}/{objective}')
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(network.parameters(), 5.)
            scaler.step(optimizer); scaler.update()
            loss_sum += loss.item() * loss_scale * len(index)
        if selecting:
            value = validation()
            history.append({'context': context, 'architecture': architecture, 'objective': objective,
                'epoch': epoch, 'train_mse': loss_sum / len(rows), 'validation_mse': value,
                'learning_rate': optimizer.param_groups[0]['lr']})
            if value < best - 1e-6:
                best, best_epoch = value, epoch
                best_state = {k: v.detach().cpu().clone() for k, v in network.state_dict().items()}
            log(f'{context}/{objective}/{architecture}: epoch {epoch}, val MSE {value:.5f}, best {best:.5f} (epoch {best_epoch})')
            if epoch - best_epoch >= args.patience:
                break
        elif epoch % 5 == 0 or epoch == epochs:
            log(f'{context}/{objective}/{architecture}: production refit epoch {epoch}/{epochs}')
    if selecting:
        network.load_state_dict(best_state)
        model['best_epoch'], model['validation_mse'] = np.array(best_epoch), np.array(best)
    else:
        model['best_epoch'] = np.array(production_epochs)
    model['parameters'] = np.array(sum(p.numel() for p in network.parameters()))
    model['learning_rate'] = np.array(args.learning_rate)
    model['weight_decay'] = np.array(1e-4)
    del x, y, base, optimizer, scaler
    gc.collect()
    if args.device == 'cuda':
        torch.cuda.empty_cache()
    return network.eval(), history


def write_output(path, model, network, inputs, obs, scope, prediction=None):
    if prediction is None:
        prediction = predict_network(model, network, inputs)
    metadata = obs.copy()
    metadata['prediction_scope'] = scope
    if str(model['objective']) == 'raw':
        metadata['predicted_library_counts'] = prediction.sum(axis=1, dtype=np.float64) + inputs.sum(axis=1)
    obj = ad.AnnData(prediction, obs=metadata, var=pd.DataFrame(index=pd.Index(model['output_genes'], name='gene_name')))
    obj.obsm['measured_300_raw_counts'] = inputs.astype(np.int32)
    obj.uns['input_genes'] = model['input_genes']
    obj.uns['objective'] = str(model['objective'])
    obj.uns['architecture'] = str(model['architecture'])
    obj.uns['units'] = 'expected raw counts' if str(model['objective']) == 'raw' else 'log1p(CPM)'
    obj.uns['prediction_scope'] = scope
    obj.write_h5ad(path, compression='lzf')


def evaluate_prediction(prediction, truth, reference, output_genes, context, objective, architecture, expressed):
    records = []
    per_cell, base_cell = np.zeros(len(truth)), np.zeros(len(truth))
    for start in range(0, truth.shape[1], 256):
        stop = min(start + 256, truth.shape[1])
        p, y = prediction[:, start:stop].astype(np.float64), truth[:, start:stop].astype(np.float64)
        b = reference[:, start:stop].astype(np.float64)
        se, sb = (p - y) ** 2, (b - y) ** 2
        per_cell += se.sum(axis=1); base_cell += sb.sum(axis=1)
        records.append(pd.DataFrame({'context': context, 'objective': objective, 'architecture': architecture,
            'gene_name': output_genes[start:stop], 'expressed_train_5cpm': expressed[start:stop],
            'test_mse': se.mean(axis=0), 'baseline_test_mse': sb.mean(axis=0),
            'test_spearman': spearman(p, y)}))
    table = pd.concat(records, ignore_index=True)
    relevant = table[table.expressed_train_5cpm] if objective == 'normalized' else table
    score = {'context': context, 'objective': objective, 'architecture': architecture,
        'test_mse': float(relevant.test_mse.mean()),
        'test_skill': float(1 - relevant.test_mse.sum() / relevant.baseline_test_mse.sum()),
        'median_expressed_gene_spearman': float(table.loc[table.expressed_train_5cpm, 'test_spearman'].median()),
        'test_all_genes_mse': float(table.test_mse.mean())}
    return score, table, per_cell / truth.shape[1], base_cell / truth.shape[1]


def run(args):
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    manifest, genes, targets = read_spec(ROOT)
    ipos = np.array([np.flatnonzero(genes == g)[0] for g in targets])
    opos = np.array([i for i, g in enumerate(genes) if g not in set(targets)])
    output_genes = genes[opos].astype(str)
    records, history, scores, gene_tables, cell_tables = [], [], [], [], []
    provenance = {'seed': SEED, 'structures': STRUCTURES, 'epochs': args.epochs, 'patience': args.patience,
        'batch_size': args.batch_size, 'learning_rate': args.learning_rate,
        'device': args.device, 'cuda_runtime': torch.version.cuda,
        'gpu': torch.cuda.get_device_name(0) if args.device == 'cuda' else None,
        'versions': {p: version(p) for p in ('torch', 'numpy', 'scipy', 'pandas', 'anndata')},
        'code_sha256': sha256(Path(__file__)), 'contexts': {}, 'outputs': {},
        'test_scope': 'same previously reported guide-group test as ridge; not used in neural fit/selection'}
    old_splits = pd.read_csv(RIDGE / 'cell_splits.csv')
    for context in args.contexts:
        log(f'{context}: validating source counts and preparing both output objectives')
        raw, guides, libraries, _, _, _ = load_control(ROOT, context, genes, manifest)
        source = ad.read_h5ad(ROOT / f'context_{context}.h5ad', backed='r')
        obs = source.obs.copy(); source.file.close()
        obs['cell_id'], obs['context'] = obs.index.astype(str), context
        inputs = raw[:, ipos].toarray().astype(np.float64)
        other = raw[:, opos].tocsr()
        roles = split_guides(guides)
        np.testing.assert_array_equal(old_splits[old_splits.context == context]['split'], roles)
        train, val, test = roles == 'train', roles == 'validation', roles == 'test'
        labels = np.empty(other.shape, dtype=np.float32)
        raw_labels = np.empty_like(labels)
        mean_cpm = np.empty(other.shape[1])
        for start in range(0, other.shape[1], 256):
            stop = min(start + 256, other.shape[1])
            counts = other[:, start:stop].toarray()
            raw_labels[:, start:stop] = counts
            cpm = counts.astype(np.float64) * (1e6 / libraries)[:, None]
            labels[:, start:stop] = np.log1p(cpm)
            mean_cpm[start:stop] = cpm[train].mean(axis=0)
        expressed = mean_cpm >= 5
        context_records = []
        for objective, target_values in [('normalized', labels), ('raw', raw_labels)]:
            for number, architecture in enumerate(STRUCTURES):
                seed = SEED + ord(context) * 100 + number * 10 + (objective == 'raw')
                model = metadata(inputs, other, labels, libraries, train, context, architecture, objective, 'train_only_evaluation', seed)
                model['input_genes'], model['output_genes'] = targets.astype(str), output_genes
                log(f'{context}: train {objective}/{architecture}')
                network, curve = train_network(model, inputs, target_values, train, val, expressed, args, context, architecture, objective)
                name = f'candidate_{context}_{objective}_{architecture}.npz'
                save_network(out / name, model, network)
                record = {'context': context, 'objective': objective, 'architecture': architecture,
                    'best_epoch': int(model['best_epoch']), 'validation_mse': float(model['validation_mse']),
                    'parameters': int(model['parameters']), 'seed': seed, 'model_file': name}
                records.append(record); context_records.append(record); history.extend(curve)
                del network, model
                gc.collect()
                if args.device == 'cuda':
                    torch.cuda.empty_cache()
                pd.DataFrame(records).to_csv(out / 'candidate_scores.csv', index=False)
                pd.DataFrame(history).to_csv(out / 'training_history.csv', index=False)
        selections = {objective: min([r for r in context_records if r['objective'] == objective], key=lambda r: r['validation_mse'])
                      for objective in ('normalized', 'raw')}
        # All architectures and epochs are locked before test expressions are scored.
        write_json(out / f'selected_{context}.json', selections)
        log(f'{context}: locked validation winners; evaluating each frozen architecture once')
        for record in context_records:
            objective, architecture = record['objective'], record['architecture']
            model, network = load_network(out / record['model_file'], args.device)
            prediction = predict_network(model, network, inputs[test])
            target_values = labels if objective == 'normalized' else raw_labels
            reference = (np.broadcast_to(labels[train].mean(axis=0, dtype=np.float64), prediction.shape) if objective == 'normalized'
                         else inputs[test].sum(axis=1)[:, None] * (np.asarray(other[train].sum(axis=0)).ravel() / inputs[train].sum()))
            score, table, per_cell, base_cell = evaluate_prediction(prediction, target_values[test], reference,
                output_genes, context, objective, architecture, expressed)
            score['selected'] = architecture == selections[objective]['architecture']
            score['best_epoch'] = record['best_epoch']
            scores.append(score)
            log(f'{context}/{objective}/{architecture}: test skill {score["test_skill"]:.4f}')
            if score['selected']:
                gene_tables.append(table)
                cells = obs.iloc[np.flatnonzero(test)].copy()
                cells['objective'], cells['architecture'] = objective, architecture
                cells['source_row'], cells['test_mse'], cells['baseline_test_mse'] = np.flatnonzero(test), per_cell, base_cell
                if objective == 'raw':
                    cells['predicted_library_counts'] = prediction.sum(axis=1, dtype=np.float64) + inputs[test].sum(axis=1)
                    cells['actual_library_counts_evaluation_only'] = libraries[test]
                cell_tables.append(cells)
                shutil.copyfile(out / record['model_file'], out / f'evaluation_{objective}_{context}.npz')
                write_output(out / f'heldout_{objective}_{context}.h5ad', model, network, inputs[test],
                             obs.iloc[np.flatnonzero(test)], 'heldout_test_guides', prediction)
            del network, model, prediction, reference
            gc.collect()
        for objective in ('normalized', 'raw'):
            selected = selections[objective]
            architecture = selected['architecture']
            production = metadata(inputs, other, labels, libraries, np.ones(len(inputs), dtype=bool), context,
                architecture, objective, 'all_controls_production', selected['seed'])
            production['input_genes'], production['output_genes'] = targets.astype(str), output_genes
            log(f'{context}: refit selected {objective}/{architecture} on all controls for {selected["best_epoch"]} epochs')
            target_values = labels if objective == 'normalized' else raw_labels
            network, _ = train_network(production, inputs, target_values, np.ones(len(inputs), dtype=bool),
                val, expressed, args, context, architecture, objective, selected['best_epoch'])
            save_network(out / f'model_{objective}_{context}.npz', production, network)
            if objective == 'raw':
                log(f'{context}: saving all-cell expected raw counts')
                write_output(out / f'reconstructed_raw_{context}.h5ad', production, network, inputs, obs, 'in_sample_all_controls_refit')
            del network, production
            gc.collect()
        provenance['contexts'][context] = {'source_sha256': sha256(ROOT / f'context_{context}.h5ad'),
            'split_cells': {r: int((roles == r).sum()) for r in ('train', 'validation', 'test')},
            'selected': selections, 'input_genes': 300, 'output_genes': len(output_genes)}
        for filename in [f'{prefix}_{objective}_{context}.npz' for prefix in ('evaluation', 'model') for objective in ('normalized', 'raw')] + [
            f'heldout_{objective}_{context}.h5ad' for objective in ('normalized', 'raw')] + [f'reconstructed_raw_{context}.h5ad']:
            provenance['outputs'][filename] = {'bytes': (out / filename).stat().st_size, 'sha256': sha256(out / filename)}
        pd.DataFrame(scores).to_csv(out / 'test_architecture_scores.csv', index=False)
        pd.concat(gene_tables, ignore_index=True).to_csv(out / 'per_gene_test_metrics.csv', index=False)
        pd.concat(cell_tables, ignore_index=True).to_csv(out / 'per_cell_test_metrics.csv', index=False)
        write_json(out / 'run_manifest.json', provenance)
        del raw, other, inputs, labels, raw_labels
        gc.collect()
        if args.device == 'cuda':
            torch.cuda.empty_cache()
    log('Finished architecture comparison, held-out evaluation and selected neural models')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--contexts', nargs='+', choices=['A', 'B', 'C'], default=['A', 'B', 'C'])
    parser.add_argument('--epochs', type=int, default=40)
    parser.add_argument('--patience', type=int, default=7)
    parser.add_argument('--batch-size', type=int, default=256)
    parser.add_argument('--learning-rate', type=float, default=.001)
    parser.add_argument('--device', choices=['cuda', 'cpu'], default='cuda')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    torch.set_num_threads(args.threads)
    if args.device == 'cuda' and not torch.cuda.is_available():
        raise RuntimeError('CUDA unavailable; explicitly use --device cpu if needed')
    with threadpool_limits(limits=args.threads):
        run(args)
