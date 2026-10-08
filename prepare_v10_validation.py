"""Reserve fresh V10 count-validation targets after the V8 decoder tradeoff."""
import json,gc,hashlib
import numpy as np
from threadpoolctl import threadpool_limits
from v10_methods import SEEDS
from v10_methods import ROOT,A
from fit_v7_response import read_cache,raw_inputs
from prepare_v8_validation import global_target_split
from v6_methods import moments
import build_submission as build


def main():
    A.mkdir(parents=True,exist_ok=True)
    if (A/'validation_plan.json').exists():raise FileExistsError('V10 plan already frozen')
    old=json.loads((ROOT/'experiments/v8_corrected_transfer/analysis/validation_plan.json').read_text())
    used={cell['target'] for c in old['contexts'] for cell in c['cells']}
    v9=json.loads((ROOT/'experiments/v9_activation_preserving/analysis/validation_plan.json').read_text())
    used|={cell['target'] for c in v9['contexts'] for cell in c['cells']}
    data=read_cache();genes=data['genes'];targets=data['targets'];contexts=[]
    for family in ['H1','HCT116','HEK293T']:
        x,meta,measured,source=raw_inputs(family,genes)
        v7=json.loads((ROOT/f'experiments/v7_multisource_adaptive/analysis/cell_split_{family}.json').read_text())
        excluded=(used|set(v7['targets'])) if family!='H1' else set()
        for prior in [c for c in old['contexts'] if c['family']==family]:
            context={k:v for k,v in prior.items() if k not in ['cells','selected_targets','old_v7_targets_reused']}
            pool=x[context['left']];query=moments(pool)[1]*1e6
            counts=meta[meta['sample'].isin(context['batches'])].gene_target.value_counts();chosen=[]
            for split in ['calibration','confirmation']:
                minimum=60 if family=='H1' else 150
                available=[str(t) for t in targets if t not in excluded and counts.get(t,0)>=minimum and query[np.flatnonzero(genes==t)[0]]>5 and global_target_split(t)==split]
                n=2 if family=='H1' else 4
                if len(available)<n:raise ValueError(f'{prior["context"]}: only {len(available)} fresh {split} targets at count threshold {minimum}')
                available.sort(key=lambda t:hashlib.sha256(f'V10:{t}'.encode()).hexdigest());chosen.extend(available[:n])
            cells=[]
            for seed in SEEDS:
                previous=prior['cells'][0]
                rng_controls=np.random.default_rng(seed+len(contexts)*101)
                previous={'template':rng_controls.choice(context['left'],400,replace=False).tolist(),'reference':rng_controls.choice(context['right'],400,replace=False).tolist()}
                for target in chosen:
                    actual=np.flatnonzero((meta.gene_target.to_numpy()==target)&meta['sample'].isin(context['batches']).to_numpy())
                    ti=int(np.flatnonzero(targets==target)[0]);rng=np.random.default_rng(seed+71000+ti)
                    actual=rng.choice(actual,min(200,len(actual)),replace=False)
                    cells.append({'target':target,'split':global_target_split(target),'seed':seed,'actual':actual.tolist(),
                                  'template':previous['template'],'reference':previous['reference']})
            context['cells']=cells;context['selected_targets']=chosen;context['all_previously_count_tested_targets_excluded']=family!='H1'
            contexts.append(context);print('Reserved V10',context['context'],chosen,flush=True)
        del x;gc.collect()
    build.write_json(A/'validation_plan.json',{'version':10,'contexts':contexts,'seeds':SEEDS,
        'v8_plan_sha256':build.sha256(ROOT/'experiments/v8_corrected_transfer/analysis/validation_plan.json'),
        'excluded_previous_v8_targets':sorted(used),'script_sha256':build.sha256(ROOT/'prepare_v10_validation.py'),
        'source_cache_sha256':old['source_cache_sha256'],
        'scope':'HCT/HEK targets exclude V7/V8/V9 count targets. H1 is reused public development because its small 25-target panel is exhausted; H1 is NOT fresh blind biology. New count-validation targets, control-defined strata and fixed independent templates; public source summaries were previously available. Not pristine blind biology. Metadata selection uses >=150 actual source cells per technical stratum (H1 >=60), own control CPM>5. No response-value selection.'})

if __name__=='__main__':
    with threadpool_limits(limits=2):main()
