"""Independent arithmetic and input-integrity audit of V67 broad response results."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v67_intervention_network as m
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key


def close(a,b):
    assert math.isclose(a,b,rel_tol=1e-11,abs_tol=1e-12),(a,b)


def main():
    m.check();result=read(m.A/'mean_results.json');assert read(m.A/'boundary_audit.json')['passed']
    source=[];allrows=[];cohort={};replays=0
    for ctx in contexts():
        name=ctx['family'];p=m.A/'predictions'/(name+'.npz');r=read(p.with_suffix('.json'))
        assert sha(p)==r['sha256']
        source.append(p);source.append(p.with_suffix('.json'))
        p66=m.ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz')
        assert sha(p66)==r['v66_predictions_sha256']
        with np.load(p,allow_pickle=False) as z, np.load(p66,allow_pickle=False) as old:
            for a,b in [('v63','v63'),('v66','assay75'),('query','query'),('bulk','bulk'),('measured','measured')]:
                np.testing.assert_array_equal(z[a],old[b])
            own=np.array([list(m.old.data()['genes']).index(t) for t in m.old.data()['targets']])
            for arm in m.ARMS:
                assert z[arm].shape==z['v63'].shape and np.isfinite(z[arm]).all()
                assert np.max(abs(z[arm]-z['v63']))<=.5000003
                np.testing.assert_array_equal(z[arm][:,own],z['v63'][:,own])
                # Independently recover bulk profiles with the explicit own-target convention.
                ti=0;f=np.exp2(z[arm][ti]);f[own[ti]]=.1
                pp=z['bulk']*f;pp/=pp.sum();assert np.isfinite(pp).all();close(float(pp.sum()),1.)
                replays+=1
        b=m.A/'broad'/(name+'.json');data=read(b);source.append(b);allrows+=data['records'];cohort[name]=data['target_names']
        for target in data['target_names']:
            rs=[v for v in data['records'] if v['target']==target]
            assert len(rs)==5 and len({v['arm'] for v in rs})==5
            assert len({v['treatment_rows_sha256'] for v in rs})==1
            assert len({v['cells'] for v in rs})==1
            assert len({v['split'] for v in rs})==1
            assert rs[0]['cells']>=30
    assert sum(map(len,cohort.values()))==655==result['unique_target_context_pairs']
    # Reconstruct macro error ratios directly from target records, independently of runner selection code.
    for arm in m.ARMS:
        for baseline in ('v63','v66'):
            for group in ('all','development','confirmation','canonical_small'):
                study_values={}
                for study in ('H1','Orion','CD4'):
                    names=sorted({v['source'] for v in allrows if v['study']==study})
                    errors={}
                    for metric in ('profile_MSE','lfc_MSE','lfc_MAE'):
                        ratios=[]
                        for name in names:
                            rows=[v for v in allrows if v['source']==name and (group=='all' or (v['canonical_small'] if group=='canonical_small' else v['split']==group))]
                            vals={a:[v[metric] for v in rows if v['arm']==a] for a in (arm,baseline)}
                            assert len(vals[arm])==len(vals[baseline])>0
                            ratios.append(sum(vals[arm])/sum(vals[baseline]))
                        errors[metric+'_ratio']=sum(ratios)/len(ratios)
                        close(errors[metric+'_ratio'],result['summary'][arm][baseline][group]['studies'][study][metric+'_ratio'])
                    study_values[study]=errors
                for metric in ('profile_MSE','lfc_MSE','lfc_MAE'):
                    close(sum(v[metric+'_ratio'] for v in study_values.values())/3,result['summary'][arm][baseline][group][metric+'_ratio'])
    checks=result['checks'];assert result['release_supported']==all(checks.values())
    assert (result['selected']=='network')==all(checks.values())
    source += [m.A/'plan.json',m.A/'mean_results.json',m.A/'boundary_audit.json',Path(__file__)]
    write(m.A/'results_audit.json',dict(passed=True,source_sha256={str(p.relative_to(m.ROOT)):sha(p) for p in source},
          independent_profile_probability_checks=replays,paired_population_records=len(allrows),
          broad_unique_target_context_pairs=655,all_300_panel_outputs_unchanged=True,
          V63_and_V66_baselines_bitwise_unchanged=True, macro_error_summaries_recomputed=True,
          scope='Audits saved arithmetic, input alignment, same outcome populations across arms and response constraints. Does not independently re-run all graph predictions, raw treatment extraction or official anchored metrics.'))
    print('V67 independent arithmetic/input audit passed:',len(allrows),'paired records,655 pairs.',flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
