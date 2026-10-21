"""Record the complete V85 source-only rejection and verify its exact abstention."""
import numpy as np
import v85_guide_contrast as m
from v57_pie import read,write,sha

def main():
    m.check();summary=[]
    for stage in ['teachers','linear_teachers']:
        for condition in m.CONDITIONS:
            p=m.A/stage/(condition+'.json');r=read(p);fields=['target','target_index','guide_half_cells','lanes','strong','F','dose','active_outputs','own_under_5CPM']
            summary.append(dict(stage=stage,condition=condition,full_trace_sha256=sha(p),arrays_sha256=r['sha256'],statistics_sha256=r['statistics_sha256'],source=r['source'],real=[{k:v[k] for k in fields if k in v} for v in r['records']],permuted=[{k:v[k] for k in fields if k in v} for v in r['permuted_records']],active_targets=r['active_targets'],permuted_active_targets=r['permuted_active_targets']))
    for bank in m.teachers():assert not bank['confidence'].any() and not bank['permuted_confidence'].any()
    checked=[]
    for name,path in [('H1',m.ROOT/'experiments/v66_assay_strength/analysis/predictions/H1.npz'),('CD4_D3_Rest',m.ROOT/'experiments/v84_intervention_operator/analysis/donor_predictions.npz')]:
        with np.load(path,allow_pickle=False) as z:q=z['query'].copy();bulk=z['bulk'].copy();base=z['assay75' if name=='H1' else 'v66'].copy()
        result,_=m.predict(q,base,bulk,name)
        for arm in m.ARMS:np.testing.assert_array_equal(result[arm],base)
        checked.append(dict(recipient=name,source_sha256=sha(path),all_arms_bitwise_V66=True))
    write(m.A/'source_summary.json',dict(records=summary,strict_guide_contrast_active_targets=0,source_only_negative_control_correction=True,all_corrected_nulls_abstain=True,query_checks=checked,new_populations_emitted=0,new_submission_generated=False,reason='Every corrected source target fails at least one prespecified strong-dose condition. No confidence available for any gene; all candidate corrections are exactly zero. Generating or uploading identical populations would not test a new prediction.',source_audit_sha256=sha(m.A/'source_audit.json'),script_sha256=sha(__file__)))
    print('V85 rejected at source power check; query abstention exactly verified',flush=True)
if __name__=='__main__':main()
