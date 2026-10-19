"""Independently verify that V69 adds only the declared count mechanism."""
import gc
import numpy as np
from threadpoolctl import threadpool_limits
import build_submission as b
import v69_count_bridge as m
from build_v69_submission import load
from v57_pie import read,write,sha


def main():
    manifest,genes,targets,record,model=load()
    parent=m.ROOT/'experiments/v66_assay_strength/analysis/reference_model_arrays.npz'
    with np.load(parent,allow_pickle=False) as z:
        np.testing.assert_array_equal(z['genes'],genes);np.testing.assert_array_equal(z['targets'],targets);np.testing.assert_array_equal(z['effects'],model['effects'])
    delta,trace=m.pooled(())
    np.testing.assert_array_equal(delta,model['delta']);assert trace==record['source_trace']
    source_checks=[]
    for ci,context in enumerate(manifest['contexts']):
        x,_,_,*_=b.load_control(m.ROOT,context,genes,manifest);params=m.parameters(x)
        np.testing.assert_array_equal(params['p'],model['bulk'][ci]);np.testing.assert_array_equal(params['alpha'],model['alpha'][ci]);source_checks.append(dict(context=context,cells=x.shape[0],alpha_range=[float(params['alpha'].min()),float(params['alpha'].max())]));del x;gc.collect()
    m.check()
    write(m.A/'model_audit.json',dict(passed=True,all_V66_response_means_exactly_preserved=True,all_source_dispersion_parameters_recomputed=True,all_query_NTC_parameters_recomputed=True,contexts=source_checks,parent_arrays_sha256=sha(parent),model_arrays_sha256=sha(m.A/'reference_model_arrays.npz'),verifier_sha256=sha(m.ROOT/'verify_v69_model.py')))
    print('V69 full mean/dispersion model verification passed',flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=2):main()
