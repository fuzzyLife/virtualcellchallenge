"""Checks that CSV export preserves count orientation, identifiers and target flags."""
import csv,json,tempfile,unittest
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import export_competition_csv as export
class ExportTests(unittest.TestCase):
 def test_export_matches_independent_dense_transpose_and_detects_corruption(self):
  with tempfile.TemporaryDirectory() as name:
   root=Path(name);genes=['G0','G1','G2'];contexts=['A','B','C']
   manifest={'season':'2026','partition':'val','panel_id':'synthetic','contexts':contexts,'control_label':'non-targeting','n_genes':3,'n_constructs':1,'cells_per_pert':2,
    'per_context':{c:{'control_cells':4,'n_ntc_ids':2,'ground_truth_cells':6} for c in contexts}}
   (root/'manifest.json').write_text(json.dumps(manifest));pd.DataFrame({'gene_name':genes}).to_csv(root/'gene_names.csv',index=False);pd.DataFrame({'target_gene':['G1']}).to_csv(root/'pert_counts.csv',index=False)
   expected=[]
   for i,c in enumerate(contexts):
    matrix=np.array([[0,2,7],[5,0,1],[0,3,0],[4,2,0]],dtype=np.int32)*(i+1);expected.append(matrix.T)
    obs=pd.DataFrame({'context':[c]*4,'target_gene':['non-targeting']*4,'ntc_id':['g1','g1','g2','g2']},index=[f'{c}_cell{i}' for i in range(4)])
    ad.AnnData(sparse.csr_matrix(matrix),obs=obs,var=pd.DataFrame(index=genes)).write_h5ad(root/f'context_{c}.h5ad')
   cache=root/'analysis/competition_csv_cache';export.extract(root,cache);path=root/'combined.csv';export.write_csv(root,cache,path,False,2);export.verify_csv(root,cache,path)
   table=pd.read_csv(path);columns=[c for c in table if c.startswith('context_') and '_col' in c]
   np.testing.assert_array_equal(table[columns].to_numpy(),np.concatenate(expected,axis=1))
   self.assertEqual(table.gene_name.tolist(),genes);self.assertEqual(table.is_requested_knockdown_target.tolist(),[0,1,0]);self.assertEqual(table.requested_cells_per_target_context.tolist(),[0,2,0])
   self.assertTrue(all(table.observed_target_gene=='non-targeting'));self.assertTrue(all(table.has_observed_knockdown_responses==0))
   self.assertIn('cell_id=A_cell0__ntc_id=g1',columns[0]);self.assertEqual(len(columns),12)
   info=json.loads((root/'analysis/combined_csv_export.json').read_text());self.assertTrue(info['validation']['every_numeric_value_roundtripped'])
   with path.open(newline='') as f:rows=list(csv.reader(f))
   rows[1][len(info['annotation_columns'])]='1'
   with path.open('w',newline='') as f:csv.writer(f,lineterminator='\n').writerows(rows)
   with self.assertRaises(ValueError):export.verify_csv(root,cache,path)
if __name__=='__main__':unittest.main(verbosity=2)