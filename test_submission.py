"""Focused checks for sampling, count simulation, streaming H5AD and official prep."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import build_submission as build


class SubmissionTests(unittest.TestCase):
    def test_stratification_balance_and_reproducibility(self):
        guides = np.repeat([f'ntc-{i}' for i in range(46)], 400)
        libraries = np.tile(np.arange(1, 401), 46)
        a = build.guide_stratified_sample(guides, libraries, 400, np.random.default_rng(7))
        b = build.guide_stratified_sample(guides, libraries, 400, np.random.default_rng(7))
        self.assertTrue(np.array_equal(a, b))
        self.assertEqual(len(np.unique(a)), 400)
        counts = np.unique(guides[a], return_counts=True)[1]
        self.assertEqual(counts.min(), 8)
        self.assertEqual(counts.max(), 9)
        self.assertAlmostEqual(libraries[a].mean(), 200.5, delta=8)

    def test_quantile_bin_boundaries_never_duplicate_cells(self):
        guides = np.repeat([f'non-targeting-{i}' for i in range(46)], 400)
        libraries = np.tile(np.arange(400), 46)
        for target_index in range(200):
            rng = np.random.default_rng(np.random.SeedSequence([20261003, 0, target_index]))
            rows = build.guide_stratified_sample(guides, libraries, 400, rng)
            self.assertEqual(len(np.unique(rows)), 400)

    def test_simulation_response_and_unchanged_input(self):
        x = sparse.csr_matrix(np.tile([100, 50, 20], (1000, 1)), dtype=np.int32)
        original = x.copy()
        result = build.simulate_counts(x, np.arange(1000), np.array([.1, 2., 1.]), np.random.default_rng(3))
        means = np.asarray(result.mean(axis=0)).ravel()
        np.testing.assert_allclose(means, [10, 100, 20], atol=.6)
        self.assertEqual((x != original).nnz, 0)
        self.assertTrue(np.issubdtype(result.dtype, np.integer))
        self.assertTrue((result.data > 0).all())

    def test_factors_mapping_and_unseen_target_fallback(self):
        effects = np.array([[[1., -1.], [0., 0.]], [[.5, -.5], [0., 0.]]])
        available = np.array([[True, False], [True, False]])
        basal = np.array([[20., 80.], [30., 70.]])
        model = (effects, available, basal, np.array([0, 2]),
                 {'shrinkage': .5, 'temperature': 0.}, ['s1', 's2'], {})
        genes = np.array(['G1', 'T1', 'G2', 'T2'])
        factors, _ = build.context_factors(np.array([20., 5., 80., 5.]), model,
                                            genes, np.array(['T1', 'T2']), .9)
        self.assertGreater(factors[0, 0], 1)
        self.assertLess(factors[0, 2], 1)
        self.assertAlmostEqual(factors[0, 1], .1)
        np.testing.assert_allclose(factors[1], [1., 1., 1., .1])

    def test_zero_aware_increases_activate_zeros_and_preserve_expected_mean(self):
        # Known synthetic means; this verifies generator behavior, not VCC skill.
        x = sparse.csr_matrix(np.column_stack([np.tile([0, 20], 10000),
                            np.full(20000, 100), np.zeros(20000)]), dtype=np.int32)
        mean = np.asarray(x.mean(axis=0)).ravel()
        before = x.copy()
        out = build.simulate_counts_zero_aware(x, np.arange(len(x.indptr)-1),
                        np.array([2., .5, 2.]), np.random.default_rng(119), mean, 110.)
        np.testing.assert_allclose(np.asarray(out.mean(axis=0)).ravel(), [20., 50., 0.], atol=.2)
        self.assertGreater(out[np.arange(0, 20000, 2), 0].nnz, 9000)
        self.assertEqual(out[:, 2].nnz, 0)
        self.assertEqual((x != before).nnz, 0)
        self.assertTrue(np.issubdtype(out.dtype, np.integer))
        self.assertTrue((out.data > 0).all())

    def test_cached_rank_test_matches_scipy_with_ties_and_zeros(self):
        import audit_prediction
        self.assertTrue(audit_prediction.check_rank_kernel()['passed'])

    def test_reference_suppression_and_shrinkage_are_separate_settings(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp=Path(directory);(tmp/'reference').mkdir()
            for name in ['default','candidate']:(tmp/name).mkdir()
            genes=np.array(['G0','G1','G2','G3']);targets=np.array(['G0','G3'])
            with h5py.File(tmp/'reference/K562_gwps_raw_bulk_01.h5ad','w') as f:
                f.create_dataset('X',data=np.array([[100,100,100,100],[100,100,100,100],[10,150,300,200]],dtype=float))
                obs=f.create_group('obs')
                obs.create_dataset('gene_transcript',data=np.array(['0_non-targeting_P1','1_non-targeting_P2','2_G0_P1'],dtype=object),dtype=h5py.string_dtype())
                obs.create_dataset('num_cells_filtered',data=[100.,120.,70.])
                obs.create_dataset('fold_expr',data=[1.,1.,.1])
                obs.create_dataset('core_control',data=[True,True,False])
                var=f.create_group('var')
                categories=var.create_dataset('categories',data=genes.astype(object),dtype=h5py.string_dtype())
                encoded=var.create_dataset('gene_name',data=np.arange(4))
                encoded.attrs['categories']=categories.ref
            default=build.reference_model(tmp,genes,targets,tmp/'default')
            candidate=build.reference_model(tmp,genes,targets,tmp/'candidate',shrinkage=.7,noise_z=0.)
            self.assertEqual(default[4]['shrinkage'],.3)
            self.assertEqual(default[6]['noise_z_assumption'],1.96)
            self.assertEqual(candidate[4]['shrinkage'],.7)
            self.assertEqual(candidate[6]['noise_z_assumption'],0.)
            self.assertTrue((np.abs(candidate[0])>=np.abs(default[0])-1e-12).all())
            self.assertTrue(default[1][0,0])
            self.assertFalse(default[1][0,1])

    def test_analysis_outputs_are_isolated_from_input_root(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp=Path(directory);genes=np.array(['G0','G1','G2','G3'])
            manifest={'n_genes':4,'n_constructs':2,'control_label':'non-targeting',
                      'contexts':['A','B','C'],'cells_per_pert':2,
                      'per_context':{c:{'control_cells':4,'n_ntc_ids':2} for c in ['A','B','C']}}
            (tmp/'manifest.json').write_text(json.dumps(manifest))
            pd.DataFrame({'gene_name':genes}).to_csv(tmp/'gene_names.csv',index=False)
            pd.DataFrame({'target_gene':['G0','G3']}).to_csv(tmp/'pert_counts.csv',index=False)
            for context in manifest['contexts']:
                obs=pd.DataFrame({'target_gene':['non-targeting']*4,'context':[context]*4,
                                  'ntc_id':['g1','g1','g2','g2']},index=[f'c{i}' for i in range(4)])
                ad.AnnData(sparse.csr_matrix(np.array([[1,2,3,4]]*4),dtype=np.int32),obs=obs,
                           var=pd.DataFrame(index=genes)).write_h5ad(tmp/f'context_{context}.h5ad')
            command=[sys.executable,str(Path(build.__file__)), '--root',str(tmp),
                     '--analysis-only','--output-dir','experiments/isolated']
            result=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertTrue((tmp/'experiments/isolated/analysis/control_summary.json').exists())
            self.assertFalse((tmp/'analysis').exists())
            self.assertFalse((tmp/'prediction.h5ad').exists())

    def test_changed_model_cannot_overwrite_first_submission_root(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp=Path(directory)
            result=subprocess.run([sys.executable,str(Path(build.__file__)), '--root',str(tmp),
                                   '--noise-z','0','--force'],capture_output=True,text=True)
            self.assertEqual(result.returncode,2)
            self.assertIn('require a separate --output-dir',result.stderr)
            self.assertEqual(list(tmp.iterdir()),[])

    def test_stream_roundtrip_and_official_package(self):
        # Full official gene axis and all contexts, but an explicitly synthetic
        # two-target/two-cell panel used only within this disposable test.
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            genes = np.array([f'G{i}' for i in range(18533)])
            targets = np.array(['G0', 'G1'])
            contexts = ['A', 'B', 'C']
            obs = pd.DataFrame({'target_gene': pd.Categorical(np.tile(np.repeat(targets, 2), 3)),
                                'context': pd.Categorical(np.repeat(contexts, 4))},
                               index=[f'p{i}' for i in range(12)])
            path = tmp/'synthetic.h5ad'
            writer = build.SparseH5Writer(path, obs, genes)
            block = sparse.csr_matrix((np.array([3, 5, 7, 2]),
                                       (np.array([0, 0, 1, 1]), np.array([0, 100, 1, 200]))),
                                      shape=(2, 18533), dtype=np.int32)
            for _ in range(6):
                writer.append(block)
            writer.close()
            result = build.validate_h5ad(path, genes, targets, contexts, 2)
            self.assertEqual(result['nnz'], 24)
            read = ad.read_h5ad(path)
            self.assertEqual((read.X[:2] != block).nnz, 0)
            with h5py.File(path, 'r') as f:
                self.assertEqual(f['X/indptr'].dtype, np.dtype('int64'))
            pd.DataFrame({'gene_name': genes}).to_csv(tmp/'genes.csv', index=False)
            pd.DataFrame({'target_gene': targets}).to_csv(tmp/'targets.csv', index=False)
            command = [sys.executable, '-m', 'vcc', '--json', 'prep', str(path),
                       '-g', str(tmp/'genes.csv'), '--perts', str(tmp/'targets.csv'),
                       '--cells-per-pert', '2', '-o', str(tmp/'synthetic.vcc')]
            run = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            from vcc.vccfile import validate_vcc, read_vcc_meta
            validate_vcc(tmp/'synthetic.vcc')
            meta = read_vcc_meta(tmp/'synthetic.vcc')
            self.assertEqual(meta['nnz'], 24)
            self.assertTrue(json.loads(run.stdout)['verified_targets'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
