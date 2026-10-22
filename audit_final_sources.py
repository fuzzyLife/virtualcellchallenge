"""Audit final-response source axes, measured masks, checksums and exact reused rows."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from fit_v7_response import read_cache
import build_submission as b
R=Path(__file__).resolve().parent;A=R/'experiments/final_v91/analysis'
def main():
 rec=json.loads((A/'source_provenance.json').read_text());assert b.sha256(A/'source_cache.npz')==rec['cache_sha256'];m,g,t=b.read_spec(R/'final');old=read_cache()
 with np.load(A/'source_cache.npz',allow_pickle=False) as z:d={k:z[k].copy() for k in z.files}
 np.testing.assert_array_equal(d['genes'],g);np.testing.assert_array_equal(d['targets'],t);np.testing.assert_array_equal(d['profiles'],old['profiles']);n=len(d['truth']);assert all(d[k].shape==(n,len(g)) for k in ['truth','errors','q','valid']);assert all(d[k].shape==(n,) for k in ['row_targets','row_contexts','row_families','knockdown']);assert ((d['row_targets']>=0)&(d['row_targets']<300)).all();assert ((d['row_contexts']>=0)&(d['row_contexts']<7)).all();assert np.isfinite(d['truth']).all() and abs(d['truth']).max()<=3;assert np.isfinite(d['errors']).all();assert (d['errors'][d['valid']]>0).all();assert not d['valid'][np.arange(n),[list(g).index(t[i]) for i in d['row_targets']]].any();assert not pd.DataFrame({'c':d['row_contexts'],'t':d['row_targets']}).duplicated().any()
 oldname=old['targets'][old['row_targets']];mask=np.isin(oldname,t);nold=rec['reused_rows'];assert nold==mask.sum()
 for k in ['truth','errors','q','valid','knockdown','row_contexts','row_families']:np.testing.assert_array_equal(d[k][:nold],old[k][mask])
 np.testing.assert_array_equal(t[d['row_targets'][:nold]],oldname[mask]);assert not np.isin(t[d['row_targets'][nold:]],old['targets']).any();assert np.isfinite(d['q']).all() and ((d['q']>=0)&(d['q']<=1)).all();assert np.all(~d['valid']|np.isfinite(d['profiles'][d['row_contexts']]))
 coverage=pd.read_csv(A/'target_coverage.csv');assert coverage.target_gene.to_list()==list(t);summary={}
 for kind,rows in [('all',coverage),('shared',coverage[coverage.reused_validation_target]),('new',coverage[~coverage.reused_validation_target])]:summary[kind]=dict(targets=len(rows),source_family_count_histogram={str(int(k)):int(v) for k,v in rows.n_source_families.value_counts().sort_index().items()},with_any_response=int(rows.any_output_measured.sum()))
 b.write_json(A/'source_audit.json',dict(passed=True,source_cache_sha256=rec['cache_sha256'],rows=n,reused_rows=nold,new_rows=n-nold,reused_labels_errors_masks_bitwise_equal=True,final_axes_exact=True,own_coordinates_masked=True,missingness_mask_retained=True,final_treatment_truth_used=False,coverage=summary,auditor_sha256=b.sha256(Path(__file__))));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
