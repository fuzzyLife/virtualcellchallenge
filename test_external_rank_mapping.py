import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
import json
import hashlib
import numpy as np
from scipy.stats import spearmanr
from cluster_full_rank_matrix import rank_unit_rows
from extract_external_mapping_controls import collapse_to_axis, CachedRangeReader
from map_external_rank_controls import frozen_projection


class ExternalRankMappingTests(unittest.TestCase):
    def test_duplicate_symbols_sum_and_missing_symbols_reject(self):
        raw=np.array([[1,2,0,9],[0,0,4,7]])
        out=collapse_to_axis(raw,['G1','G1','G2','OTHER'],['G2','G1']).toarray()
        np.testing.assert_array_equal(out,[[0,3],[4,0]])
        with self.assertRaises(ValueError):
            collapse_to_axis(raw,['G1','G1','G2','OTHER'],['G1','MISSING'])

    def test_frozen_projection_matches_spearman_without_refitting(self):
        raw=np.array([[0,0,1,3,8,20],[20,8,3,1,0,0],[0,20,1,8,3,0]])
        z=rank_unit_rows(raw); mean=z.mean(axis=0)
        components=np.linalg.qr(np.random.default_rng(7).normal(size=(6,2)))[0].T
        pcs=(z-mean)@components.T
        model={'centroids':z.copy(),'context_mean_ranks':z.copy(),
               'pca_mean':mean,'pca_components':components,
               'pc_context_means':pcs,'pc_context_cov_inverse':np.array([np.eye(2)]*3),
               'pc_context_mahal_p95':np.array([100]*3),'context_similarity_p05':np.array([-1]*3)}
        before={k:v.copy() for k,v in model.items()}
        result=frozen_projection(raw,model)
        self.assertEqual(result.nearest_context.tolist(),['A','B','C'])
        np.testing.assert_allclose(result[['PC1','PC2']],pcs,atol=1e-7)
        for i in range(3):
            for j,c in enumerate(['A','B','C']):
                self.assertAlmostEqual(result['mean_cell_spearman_'+c].iloc[i],spearmanr(raw[i],raw[j]).statistic,places=6)
        perm=frozen_projection(raw[[2,0,1]],model)
        self.assertEqual(perm.nearest_context.tolist(),['C','A','B'])
        for key in model:np.testing.assert_array_equal(model[key],before[key])
        with self.assertRaises(ValueError):frozen_projection([[0,np.nan,1,3,8,20]],model)

    def test_offline_reader_verifies_blocks_and_rejects_corruption(self):
        with TemporaryDirectory() as directory:
            path=Path(directory); raw=b'0123456789abcdef'
            identity={'url':'https://example.invalid/source','bytes':len(raw),'etag':'pinned'}
            (path/'identity.json').write_text(json.dumps(identity))
            block=path/'000000000.bin';block.write_bytes(raw)
            block.with_suffix('.sha256').write_text(hashlib.sha256(raw).hexdigest())
            reader=CachedRangeReader(path,identity)
            reader.seek(4);self.assertEqual(reader.read(5),b'45678')
            reader.seek(-4,2);self.assertEqual(reader.read(),b'cdef')
            self.assertEqual(reader.used,len(raw));reader.close()
            block.write_bytes(b'X'+raw[1:])
            with self.assertRaises(ValueError):CachedRangeReader(path,identity).read(2)
            block.unlink()
            with self.assertRaises(FileNotFoundError):CachedRangeReader(path,identity).read(2)


if __name__=='__main__':unittest.main()
