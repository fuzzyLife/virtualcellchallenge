"""Extract the previously selected external controls on a finite common gene axis."""
from pathlib import Path
from datetime import datetime, timezone
import io
import json
import hashlib
import numpy as np
import pandas as pd
import h5py
import pyarrow.parquet as pq
from scipy import sparse
from analyze_target_rank_correlations import ROOT, REFERENCES, digest
from prepare_reference_controls import column

CACHE = ROOT / 'reference/external_mapping'
PROFILE_DIR = ROOT / 'reference/control_profiles'


class CachedRangeReader(io.RawIOBase):
    """Offline HDF5 reader; each consumed block must have a valid SHA256 sidecar."""
    def __init__(self, cache, expected_identity, budget=2*1024**3):
        self.cache = Path(cache)
        self.identity = json.loads((self.cache / 'identity.json').read_text())
        if self.identity != expected_identity:
            raise ValueError('Cached source identity differs from verified profile')
        self.size = self.identity['bytes']
        self.block_size = 262144
        self.position = 0
        self.blocks = {}
        self.verified = {}
        self.used = 0
        self.budget = budget

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.position
    def seek(self, offset, whence=0):
        self.position = offset if whence == 0 else self.position + offset if whence == 1 else self.size + offset
        if self.position < 0:
            raise ValueError('Negative file offset')
        return self.position

    def read(self, count=-1):
        end = min(self.size, self.position + count) if count >= 0 else self.size
        parts = []
        while self.position < end:
            block = self.position // self.block_size
            lo = block * self.block_size
            hi = min(self.size, lo + self.block_size)
            if block not in self.blocks:
                if self.used + hi - lo > self.budget:
                    raise ValueError('Offline byte budget exceeded')
                path = self.cache / f'{block:09d}.bin'
                sidecar = path.with_suffix('.sha256')
                if not path.exists() or not sidecar.exists():
                    raise FileNotFoundError('Missing CD4 source range; run prepare_reference_controls.py --source CD4: ' + str(path))
                raw = path.read_bytes()
                sha = hashlib.sha256(raw).hexdigest()
                if len(raw) != hi-lo or sha != sidecar.read_text().strip():
                    raise ValueError('Cached range checksum mismatch')
                self.blocks[block] = raw
                self.verified[str(block)] = sha
                self.used += len(raw)
            stop = min(end, hi)
            parts.append(self.blocks[block][self.position-lo:stop-lo])
            self.position = stop
        return b''.join(parts)

    def readinto(self, buffer):
        value = self.read(len(buffer))
        buffer[:len(value)] = value
        return len(value)


def collapse_to_axis(matrix, source_genes, axis):
    """Sum duplicate measured symbols; require every requested symbol to exist."""
    source_genes = np.asarray(source_genes, dtype=str)
    axis = np.asarray(axis, dtype=str)
    if len(set(axis)) != len(axis) or not set(axis).issubset(set(source_genes)):
        raise ValueError('Requested gene axis is missing or duplicated')
    positions = {gene: i for i, gene in enumerate(axis)}
    rows = np.array([i for i, gene in enumerate(source_genes) if gene in positions], dtype=int)
    columns = np.array([positions[source_genes[i]] for i in rows], dtype=int)
    projection = sparse.csr_matrix((np.ones(len(rows)), (rows, columns)), shape=(len(source_genes), len(axis)))
    return sparse.csr_matrix(matrix, dtype=np.float64) @ projection


def load_profiles():
    prior = json.loads((ROOT / 'analysis/cell_reference_correlations/analysis_report.json').read_text())
    records = {}; profiles = {}
    for name in REFERENCES + ['CD4_Rest', 'CD4_Stim8hr', 'CD4_Stim48hr']:
        path = PROFILE_DIR / f'{name}.csv'
        record = json.loads(path.with_suffix('.json').read_text())
        expected = prior['reference_provenance'].get(name) if name in REFERENCES else prior['extra_CD4_condition_provenance'][name]
        if digest(path) != record['profile_sha256'] or record['profile_sha256'] != expected['profile_sha256']:
            raise ValueError('Previously verified external profile changed: ' + name)
        frame = pd.read_csv(path)
        if frame.gene_name.duplicated().any():
            raise ValueError('Duplicate pooled symbols')
        profiles[name] = frame.set_index('gene_name').control_count_sum
        records[name] = record
        if name in REFERENCES:
            selected = PROFILE_DIR / f'{name}_selected_controls.csv'
            if digest(selected) != expected['selected_controls_sha256']:
                raise ValueError('Previously selected control rows changed: ' + name)
            records[name] = dict(record, selected_controls_sha256=digest(selected))
    genes = pd.read_csv(ROOT / 'gene_names.csv').gene_name.to_numpy()
    axis = np.array([gene for gene in genes if all(gene in profiles[name].index for name in REFERENCES)])
    if len(axis) != 7447:
        raise ValueError('Expected the verified 7,447-gene common axis')
    return axis, profiles, records


def extract():
    CACHE.mkdir(parents=True, exist_ok=True)
    axis, profiles, records = load_profiles()
    matrices = []; metadata = []; sources = {}; reader_record = None
    for name in REFERENCES:
        chosen = pd.read_csv(PROFILE_DIR / f'{name}_selected_controls.csv')
        rows = chosen.source_row.to_numpy(dtype=int)
        if not np.array_equal(rows, np.unique(rows)):
            raise ValueError('Selected source rows must be unique and ordered')
        if name == 'K562':
            path = ROOT / 'reference/K562_gwps_raw_bulk_01.h5ad'
            if digest(path) != records[name]['source_sha256']:
                raise ValueError('K562 original checksum changed')
            with h5py.File(path, 'r') as f:
                ds = f['var/gene_name']
                genes = f[ds.attrs['categories']].asstr()[:][ds[:]]
                membership = f['obs/num_cells_filtered'][rows]
                if not np.array_equal(membership, chosen.n_cells):
                    raise ValueError('K562 group membership changed')
                values = f['X'][rows, :].astype(np.float64) * membership[:, None]
                counts = collapse_to_axis(values, genes, axis)
                chosen['native_total_counts'] = values.sum(axis=1)
            sources[name] = {'file': str(path.relative_to(ROOT)), 'sha256': digest(path)}
            chosen['unit_kind'] = 'pseudobulk_reconstructed_sum'
        elif name in ['HCT116', 'HEK293T']:
            path = ROOT / f'reference/orion_pilot/{name}_Batch1.parquet'
            gene_path = ROOT / 'reference/orion_pilot/gene_metadata.parquet'
            if digest(path) != records[name]['source_sha256'] or digest(gene_path) != records[name]['gene_metadata_sha256']:
                raise ValueError('Orion original checksum changed')
            gene_meta = pq.read_table(gene_path).to_pandas()
            destination = np.full(int(gene_meta.gene_token_id.max()) + 1, -1, dtype=int)
            positions = {gene:i for i,gene in enumerate(axis)}
            for row in gene_meta.itertuples():
                if row.gene_name in positions:
                    destination[int(row.gene_token_id)] = positions[row.gene_name]
            original = pq.read_table(path, columns=['gene_target','pass_guide_filter','total_counts','gene_token_id','gene_expression'])
            raw = np.zeros((len(rows), len(axis)))
            for i, row in enumerate(rows):
                if str(original['gene_target'][int(row)].as_py()).lower() != 'non-targeting' or original['pass_guide_filter'][int(row)].as_py() <= 0:
                    raise ValueError('Unexpected Orion control annotation')
                tokens = np.array(original['gene_token_id'][int(row)].as_py(), dtype=int)
                values = np.array(original['gene_expression'][int(row)].as_py(), dtype=float)
                if not np.isfinite(values).all() or np.any(values < 0) or values.sum() != chosen.total_counts.iloc[i]:
                    raise ValueError('Invalid Orion original counts')
                dest = destination[tokens]; present = dest >= 0
                raw[i] = np.bincount(dest[present], weights=values[present], minlength=len(axis))
            counts = sparse.csr_matrix(raw)
            chosen['native_total_counts'] = chosen.total_counts
            chosen['n_cells'] = 1
            chosen['unit_kind'] = 'single_cell'
            sources[name] = {'file': str(path.relative_to(ROOT)), 'sha256': digest(path), 'metadata_sha256': digest(gene_path)}
        elif name == 'H1':
            path = ROOT / 'reference/h1_response/selected_counts.npz'
            extraction = json.loads((path.parent / 'extraction.json').read_text())
            for file, field in [('selected_counts.npz','counts_sha256'), ('selected_cells.csv','selected_cells_sha256'), ('gene_names.csv','gene_names_sha256')]:
                if digest(path.parent / file) != extraction[field]:
                    raise ValueError('H1 cached extraction checksum changed')
            meta = pd.read_csv(path.parent / 'selected_cells.csv')
            indices = meta.index[meta.target_gene == 'non-targeting'].to_numpy()
            np.testing.assert_array_equal(meta.source_row.iloc[indices], rows)
            matrix = sparse.load_npz(path)[indices]
            genes = pd.read_csv(path.parent / 'gene_names.csv').gene_name.to_numpy()
            counts = collapse_to_axis(matrix, genes, axis)
            chosen['native_total_counts'] = np.asarray(matrix.sum(axis=1)).ravel()
            chosen['n_cells'] = 1
            chosen['unit_kind'] = 'single_cell'
            sources[name] = extraction
        else:
            reader = CachedRangeReader(PROFILE_DIR / 'ranges_CD4', records[name]['remote'])
            with h5py.File(reader, 'r') as f:
                genes = column(f['var/gene_name'])
                x = f['X']; pointer = x['indptr'][:]
                values = []; indices = []; indptr = [0]
                for row in rows:
                    lo, hi = int(pointer[row]), int(pointer[row+1])
                    v = x['data'][lo:hi].astype(np.float64); j = x['indices'][lo:hi]
                    values.append(v); indices.append(j); indptr.append(indptr[-1]+len(v))
                matrix = sparse.csr_matrix((np.concatenate(values), np.concatenate(indices), np.asarray(indptr)), shape=(len(rows), len(genes)))
                counts = collapse_to_axis(matrix, genes, axis)
                chosen['native_total_counts'] = np.asarray(matrix.sum(axis=1)).ravel()
                np.testing.assert_allclose(chosen.native_total_counts, chosen.total_counts, rtol=0, atol=1)
            reader_record = {'identity': reader.identity, 'bytes_read': reader.used, 'blocks_sha256': reader.verified, 'network_bytes': 0}
            reader.close()
            chosen['unit_kind'] = 'pseudobulk_sum'
            sources[name] = reader_record
        counts.sum_duplicates(); counts.eliminate_zeros()
        if not np.isfinite(counts.data).all() or np.any(counts.data < 0):
            raise ValueError('Invalid external control counts')
        # Independently reconstructed units must recover the previously pooled profile.
        np.testing.assert_allclose(np.asarray(counts.sum(axis=0)).ravel(), profiles[name].loc[axis], rtol=2e-7, atol=1e-5)
        chosen.insert(0, 'unit_id', [f'{name}_row_{row:08d}' for row in rows])
        chosen.insert(1, 'reference', name)
        chosen['common_total_counts'] = np.asarray(counts.sum(axis=1)).ravel()
        chosen['common_detected_genes'] = np.diff(counts.indptr)
        matrices.append(counts); metadata.append(chosen)
        print(f'Extracted {name}: {len(chosen):,} {chosen.unit_kind.iloc[0]} units', flush=True)
    combined = sparse.vstack(matrices, format='csr')
    meta = pd.concat(metadata, ignore_index=True)
    if meta.unit_id.duplicated().any() or combined.shape != (4393, len(axis)):
        raise ValueError('Unexpected external unit shape')
    sparse.save_npz(CACHE / 'counts.npz', combined)
    meta.to_csv(CACHE / 'units.csv', index=False)
    pd.DataFrame({'gene_name': axis}).to_csv(CACHE / 'gene_axis.csv', index=False)
    report = {'created_utc': datetime.now(timezone.utc).isoformat(), 'units': len(meta), 'genes': len(axis),
              'source_units': meta.groupby('reference').size().to_dict(),
              'source_cell_memberships': meta.groupby('reference').n_cells.sum().to_dict(),
              'profiles': records, 'sources': sources, 'network_bytes_downloaded': 0,
              'script_sha256': digest(Path(__file__)),
              'outputs': {name: digest(CACHE / name) for name in ['counts.npz','units.csv','gene_axis.csv']}}
    (CACHE / 'extraction.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('Extracted all 4,393 external control units without network access', flush=True)
    return report


if __name__ == '__main__':
    extract()
