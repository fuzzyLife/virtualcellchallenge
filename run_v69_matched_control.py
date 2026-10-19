"""Additional matched quarter-mixture shuffle control; no change to selection."""
import argparse,gc
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v69_count_bridge as m
import run_v69_experiment as runner
from run_v69_cached import de,profile
from v57_pie import read,write,sha
from run_v35_v36 import contexts,key
from v62_conserved_response import data


def freeze():
    path=m.A/'matched_control_plan.json';assert not path.exists()
    write(path,dict(helper_sha256=sha(m.ROOT/'run_v69_matched_control.py'),original_plan_sha256=sha(m.A/'plan.json'),
        reason='Pre-development code review: original shuffled_shape is a full bridge population. A hybrid25 candidate additionally needs the identical quarter mixture using shuffled shape evidence; otherwise the count mechanism and mixture weight are confounded.',
        formula='Reuse saved fold-specific V66 and full shuffled_shape count populations. Apply exactly the same hybrid() seed and100/400 row allocation as hybrid25.',
        stage='Recorded after some labelled oracle diagnostic outcomes, before any V69 development/confirmation outcome. No candidate rule, fitted parameter or response is changed. Additional matched comparison can only make support stricter.',D3_opened=False))


def run(stage):
    plan=read(m.A/'matched_control_plan.json');assert sha(m.ROOT/'run_v69_matched_control.py')==plan['helper_sha256']
    runner.de=de;runner.profile_jackknife_chunked=profile;d=data();genes=d['genes'];targets=d['targets'];dest=m.A/stage
    for ctx in contexts():
        name=ctx['family'];folder,reference,*_=runner.inputs(name)
        indices=runner.unique(ctx)[:2] if stage=='development' else runner.unique(ctx)[2:]
        for i in indices:
            path=dest/f'{name}__{i}__shuffled_hybrid25.json'
            if path.exists():continue
            baseline=dest/f'{name}__{i}__v66.json';shuffled=dest/f'{name}__{i}__shuffled_shape.json';a=read(baseline);b=read(shuffled)
            old=sparse.load_npz(baseline.with_suffix('.counts.npz'));new=sparse.load_npz(shuffled.with_suffix('.counts.npz'))
            assert sha(baseline.with_suffix('.counts.npz'))==a['counts_sha256'] and sha(shuffled.with_suffix('.counts.npz'))==b['counts_sha256']
            c=ctx['cells'][i];ti=list(targets).index(c['target']);own=list(genes).index(c['target']);seed=c['seed']+ti
            count=runner.hybrid(old,new,seed);actual=sparse.load_npz(folder/f'{i}_actual.npz')
            runner.score(path,count,reference,actual,genes,own,dict(stage=stage,family=name,study=key(name),case_index=i,target=c['target'],arm='shuffled_hybrid25',cells=count.shape[0],seed=seed,source_trace=b['source_trace'],template_sha256=b['template_sha256'],reference_sha256=b['reference_sha256'],actual_sha256=b['actual_sha256'],plan_sha256=b['plan_sha256'],oracle=False,additional_control_plan_sha256=sha(m.A/'matched_control_plan.json')))
        runner.panels(stage,ctx,['shuffled_hybrid25'],reference,genes,targets);gc.collect()
    if stage=='confirmation':
        summary=runner.summary(stage);choice=read(m.A/'selection.json')['selected'];control='shuffled_hybrid25' if choice=='hybrid25' else 'shuffled_shape'
        passed=choice=='mean_bridge' or summary['summary'][choice]['profile_square_error']<summary['summary'][control]['profile_square_error']
        write(m.A/'matched_control_assessment.json',dict(candidate=choice,control=control,beats_matched_control_profile=passed,summary=summary['summary'],selection_unchanged=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','development','confirmation']);a=p.parse_args()
    with threadpool_limits(limits=2):freeze() if a.stage=='freeze' else run(a.stage)
