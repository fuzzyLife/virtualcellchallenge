"""Focused archive schema and zero/tie marker checks."""
import unittest
import tempfile
import zipfile
from pathlib import Path
import numpy as np
from analyze_target_archive import read_archive, marker_statistics

class ArchiveTests(unittest.TestCase):
    def test_gene_coverage_and_duplicate_headers_rejected(self):
        genes=[f'g{i}' for i in range(300)]
        for header in [genes[:-1]+['unexpected'],genes[:-1]+[genes[0]]]:
            with tempfile.TemporaryDirectory() as folder:
                path=Path(folder)/'input.zip'
                with zipfile.ZipFile(path,'w') as archive:
                    archive.writestr('matrix.csv',','.join(['context']+header)+'\n')
                with self.assertRaises(ValueError):read_archive(path,genes)

    def test_multiple_archive_members_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'input.zip'
            with zipfile.ZipFile(path,'w') as archive:
                archive.writestr('one.csv','context\n');archive.writestr('two.csv','context\n')
            with self.assertRaises(ValueError):read_archive(path,[f'g{i}' for i in range(300)])

    def test_marker_ties_zeros_and_full_library_normalization(self):
        groups=np.repeat(['A','B','C'],2)
        counts=np.array([[10,1,0],[20,2,0],[0,1,0],[0,2,0],[0,1,0],[0,2,0]],float)
        libraries=np.tile([100,200],3)
        output=marker_statistics(counts,libraries,groups,['specific','equal_cpm','all_zero'])
        a=output[output.group=='A'].set_index('gene_name')
        self.assertEqual(a.loc['specific','minimum_pairwise_auc'],1)
        self.assertTrue(np.isinf(a.loc['specific','fold_over_highest_other_group']))
        self.assertTrue(a.loc['specific','detection_restricted'])
        self.assertTrue(a.loc['specific','panel_enriched'])
        self.assertEqual(a.loc['equal_cpm','minimum_pairwise_auc'],.5)
        self.assertEqual(a.loc['equal_cpm','fold_over_highest_other_group'],1)
        self.assertEqual(a.loc['all_zero','auc_vs_rest'],.5)
        self.assertTrue(np.isnan(a.loc['all_zero','fold_over_highest_other_group']))
        self.assertFalse(a.loc['all_zero','panel_enriched'])
        self.assertEqual(a.loc['equal_cpm','mean_cpm'],10000)

if __name__=='__main__':unittest.main()
