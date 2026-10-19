"""Additional V69 null check with representative, depth-stratified real controls.

Original prefix-template nulls are retained. This adds diagnostics, never changes
candidate selection or searches seeds.
"""
import argparse
import numpy as np
from scipy import sparse
from threadpoolctl import threadpool_limits
import v69_count_bridge as m
from v57_pie import read,write,sha
from run_v69_experiment import frame,inputs,old_emit,hybrid
from run_v35_v36 import contexts
from v68_population_balancing import strata,uniform,sample
from v12_cpu_metrics import de
from v62_conserved_response import data


def freeze():
    path=m.A/'representative_null_plan.json';assert not path.exists()
    write(path,dict(dependencies={n:sha(m.ROOT/n) for n in ['run_v69_representative_null.py','v69_count_bridge.py','v68_population_balancing.py']},
        original_plan_sha256=sha(m.A/'plan.json'),reason='Original null_template is x[left[:400]], which can overrepresent source batches. Original HCT116 identity had578 DE calls; HEK293T identity had1041. These observations motivate this extra sample-frame test, not parameter/seed selection. Preserve all original null outcomes.',
        sampler='Depth-stratified uniform pivotal sample from all query-fit controls, 400 distinct rows, eight bins, seed690800+context index. The original V68 sampler is reused; no outcomes from this additional test observed at freeze.',arms=['template','v66','bridge','hybrid25','near_zero_bridge'],selection_changed=False,D3_opened=False))


def run():
    p=read(m.A/'representative_null_plan.json')
    for n,h in p['dependencies'].items():assert sha(m.ROOT/n)==h
    genes=data()['genes'];out=m.A/'representative_null';out.mkdir(exist_ok=True)
    for ci,ctx in enumerate(contexts()):
        name=ctx['family'];folder,reference,_,q,bp,_=inputs(name);control=frame(name);params=m.parameters(control);seed=690800+ci
        group,quota=strata(control.sum(1).A1,400);ids=sample(uniform(group,quota,400),group,quota,400,seed);template=control[ids];zero=np.zeros(len(genes))
        old=old_emit(template,zero,None,seed,np.nan_to_num(q)/1e6,bp);full,_=m.bridge(template,zero,None,seed,params)
        tiny=np.full(len(genes),1e-6);tiny[::2]*=-1;near,_=m.bridge(template,tiny,None,seed,params)
        for arm,c in dict(template=template,v66=old,bridge=full,hybrid25=hybrid(old,full,seed),near_zero_bridge=near).items():
            path=out/f'{name}__{arm}.json'
            if path.exists():continue
            result=de(c,reference,genes);write(path,dict(family=name,arm=arm,calls=int(np.sum(result['eligible']&(result['q']<.05))),selected_frame_rows=ids.tolist(),seed=seed,exact_count_identity=(c-template).nnz==0,unique_control_rows=len(set(ids)),plan_sha256=sha(m.A/'representative_null_plan.json')))
            print('V69 representative null',name,arm,read(path)['calls'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','run']);a=p.parse_args()
    with threadpool_limits(limits=2):globals()[a.stage]()
