"""Verify that V90 changes predictions, preserves its exact V66 parent and respects response bounds."""
import numpy as np
from pathlib import Path
from v57_pie import read,write,sha
ROOT=Path(__file__).resolve().parent;A=ROOT/'experiments/v90_interval_diagnostic/analysis'
def main():
 record=read(A/'reference_model.json');p=A/'reference_model_arrays.npz';assert sha(p)==record['model_arrays_sha256'];source=ROOT/'experiments/v66_assay_strength/analysis/reference_model_arrays.npz';assert sha(source)==record['parent_V66_arrays_sha256']
 with np.load(p,allow_pickle=False) as z:effects=z['effects'].copy();parent=z['parent_effects'].copy();genes=z['genes'].copy();targets=z['targets'].copy()
 with np.load(source,allow_pickle=False) as z:np.testing.assert_array_equal(parent,z['effects']);np.testing.assert_array_equal(genes,z['genes']);np.testing.assert_array_equal(targets,z['targets'])
 own=np.array([list(genes).index(t) for t in targets]);assert (effects[:,np.arange(300),own]==0).all();assert np.isfinite(effects).all();delta=effects-parent;assert np.any(delta!=0);assert np.max(abs(delta))<=.125+2e-6
 manifest=__import__('build_submission').read_spec(ROOT)[0];rows=[]
 for ci,c in enumerate(manifest['contexts']):rows.append(dict(context=c,changed_target_groups=int(np.any(delta[ci]!=0,axis=1).sum()),maximum_absolute_LFC_change=float(abs(delta[ci]).max()),mean_absolute_LFC_change=float(abs(delta[ci]).mean()),changed_gene_coordinates=int(np.count_nonzero(delta[ci]))))
 write(A/'production_response_audit.json',dict(passed=True,parent_V66_effects_exact=True,not_identical_to_V66=True,own_effect_zero=True,quarter_strength_response_bounds_passed=True,records=rows,script_sha256=sha(__file__),scope='Response integrity and genuine change, not a hosted accuracy prediction.'));print('V90 production response audit',rows,flush=True)
if __name__=='__main__':main()
