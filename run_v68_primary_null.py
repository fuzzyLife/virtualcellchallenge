"""Additional V68 complete-pipeline null, fixed before any null outcome is opened."""
import argparse,time
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v68_population_balancing as m
from v57_pie import read,write,sha
from run_v35_v36 import contexts
from run_v68_experiment import reference_data
from v62_conserved_response import data
from v12_cpu_metrics import de


def freeze():
    m.check();path=m.A/'primary_null_plan.json';assert not path.exists()
    assert not list((m.A/'null').glob('*.json')) and not list((m.A/'primary_null').glob('*.json'))
    write(path,dict(original_plan_sha256=sha(m.A/'plan.json'),helper_sha256=sha(__file__),
          timing='After the two H1 treatment cases and while other treatment cases run, before any null outcomes are generated or opened. Source-review correction, not outcome-based model tuning.',
          reason='The original uniform-real null tests the sampler, but does not exercise balancing on a zero-response target. Add the complete fitted zero-response pipeline as a stricter additional requirement; original methods, parameters and results unchanged.',
          protocol='Five contexts,400cells, no own edit, same actual controls/reference/null seed as uniform-real null. Fit weights to the frozen real-control bulk with the same optimizer, axis, depth constraints and ridge; sample and score full Wilcoxon/BH. Require calls<=max(5,originalV63 null calls) in every context in addition to all original guards.',
          automatic_upload=False,D3_opened=False))
    print('V68 full-pipeline null addendum frozen',flush=True)


def run():
    m.check();p=read(m.A/'primary_null_plan.json');assert p['original_plan_sha256']==sha(m.A/'plan.json') and p['helper_sha256']==sha(__file__)
    d=data();genes=d['genes'];targets=d['targets'];dest=m.A/'primary_null';dest.mkdir(exist_ok=True)
    for ci,ctx in enumerate(contexts()):
        name=ctx['family'];path=dest/(name+'.json')
        if path.exists():assert sha(path.with_suffix('.npz'))==read(path)['arrays_sha256'];continue
        framepath=m.A/'frames'/name;receipt=read(framepath/'receipt.json');assert sha(framepath/'controls.npz')==receipt['controls_sha256']
        frame=sparse.load_npz(framepath/'controls.npz');folder,reference,base,bulk,query,measured=reference_data(name)
        group,quota=m.strata(frame.sum(1).A1,400);weights,diagnostic=m.fit(frame,bulk,measured&~np.isin(genes,targets),400,group,quota)
        ids=m.sample(weights,group,quota,400,680800+ci);c=frame[ids];dd=de(c,reference,genes)
        np.savez_compressed(path.with_suffix('.npz'),weights=weights,group=group,quota=quota,selected_frame_rows=ids,p=dd['p'],q=dd['q'],lfc=dd['lfc'],eligible=dd['eligible'])
        write(path,dict(source=name,cells=400,calls=int(np.sum(dd['eligible']&(dd['q']<.05))),arrays_sha256=sha(path.with_suffix('.npz')),fit=diagnostic))
        print('V68 complete-pipeline null',name,read(path)['calls'],flush=True)


def assess():
    p=read(m.A/'primary_null_plan.json');assert p['helper_sha256']==sha(__file__);rows=[]
    for ctx in contexts():
        name=ctx['family'];a=read(m.A/'null'/(name+'__v63.json'));b=read(m.A/'primary_null'/(name+'.json'))
        rows.append(dict(source=name,v63_calls=a['calls'],balanced_real_calls=b['calls'],passed=b['calls']<=max(5,a['calls'])))
    write(m.A/'primary_null_assessment.json',dict(records=rows,passed=all(r['passed'] for r in rows),plan_sha256=sha(m.A/'primary_null_plan.json')))
    print('Complete-pipeline null guards',all(r['passed'] for r in rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','run','assess']);a=p.parse_args()
    with threadpool_limits(limits=2):globals()[a.stage]()
