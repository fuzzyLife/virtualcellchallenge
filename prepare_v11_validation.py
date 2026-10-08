"""Freeze V11 reused development and previously untested count-confirmation targets."""
import json,gc,hashlib
import numpy as np
from threadpoolctl import threadpool_limits
from v11_methods import ROOT,A,SEEDS,SETTINGS
from fit_v7_response import read_cache,raw_inputs
from v6_methods import moments
import build_submission as build


def main():
    A.mkdir(parents=True,exist_ok=True)
    if (A/'validation_plan.json').exists():raise FileExistsError('V11 plan already frozen')
    old=json.loads((ROOT/'experiments/v10_population_transport/analysis/validation_plan.json').read_text())
    used=set()
    for version,directory in [(8,'v8_corrected_transfer'),(9,'v9_activation_preserving'),(10,'v10_population_transport')]:
        prior=json.loads((ROOT/f'experiments/{directory}/analysis/validation_plan.json').read_text())
        used|={c['target'] for ctx in prior['contexts'] for c in ctx['cells']}
    data=read_cache();genes=data['genes'];targets=data['targets'];contexts=[]
    for family in ['H1','HCT116','HEK293T']:
        x,meta,measured,source=raw_inputs(family,genes)
        v7=json.loads((ROOT/f'experiments/v7_multisource_adaptive/analysis/cell_split_{family}.json').read_text())
        for previous in [c for c in old['contexts'] if c['family']==family]:
            context={k:v for k,v in previous.items() if k not in ['cells','selected_targets','all_previously_count_tested_targets_excluded']}
            query=moments(x[context['left']])[1]*1e6
            # Original V10 calibration targets are reused and explicitly development.
            development=sorted({c['target'] for c in previous['cells'] if c['split']=='calibration'})
            excluded=used|set(v7['targets']) if family!='H1' else set()
            excluded|=set(development)
            counts=meta[meta['sample'].isin(context['batches'])].gene_target.value_counts()
            minimum=60 if family=='H1' else 150
            available=[str(t) for t in targets if t not in excluded and counts.get(t,0)>=minimum and query[np.flatnonzero(genes==t)[0]]>5]
            available.sort(key=lambda t:hashlib.sha256(f'V11-confirm:{t}'.encode()).hexdigest())
            count=2 if family=='H1' else 3
            if len(available)<count:raise ValueError(f'{context["context"]}: insufficient fresh confirmation targets')
            confirmation=available[:count];cells=[]
            for seed in SEEDS:
                rng=np.random.default_rng(seed+101*len(contexts))
                template=rng.choice(context['left'],400,replace=False).tolist()
                reference=rng.choice(context['right'],400,replace=False).tolist()
                for split,names in [('development',development),('confirmation',confirmation)]:
                    for target in names:
                        ids=np.flatnonzero(meta.gene_target.eq(target).to_numpy()&meta['sample'].isin(context['batches']).to_numpy())
                        ti=int(np.flatnonzero(targets==target)[0]);actual=np.random.default_rng(seed+ti+71000).choice(ids,min(200,len(ids)),replace=False)
                        cells.append({'target':target,'split':split,'seed':seed,'actual':actual.tolist(),'template':template,'reference':reference})
            context.update(cells=cells,development_targets=development,confirmation_targets=confirmation,
                confirmation_count_targets_previously_untested=family!='H1')
            contexts.append(context);print('V11 plan',context['context'],'development',development,'confirmation',confirmation,flush=True)
        del x;gc.collect()
    development={c['target'] for ctx in contexts for c in ctx['cells'] if c['split']=='development'}
    confirmation={c['target'] for ctx in contexts for c in ctx['cells'] if c['split']=='confirmation'}
    if development&confirmation:raise ValueError('Global target overlap: '+str(development&confirmation))
    build.write_json(A/'validation_plan.json',{'version':11,'contexts':contexts,'seeds':SEEDS,'settings':SETTINGS,
        'excluded_previous_count_targets':sorted(used),'source_cache_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),
        'v10_plan_sha256':build.sha256(ROOT/'experiments/v10_population_transport/analysis/validation_plan.json'),
        'script_sha256':build.sha256(ROOT/'prepare_v11_validation.py'),'methods_sha256':build.sha256(ROOT/'v11_methods.py'),
        'scope':'V10 calibration cases reused as development. New HCT/HEK count-confirmation targets exclude V7/V8/V9/V10; H1 remains reused development biology. Source summaries previously inspected. Whole held family excluded from response fitting. Shared target-independent 400-cell templates and 400/full independent references; template depths indexed correctly. Not pristine blind validation.'})
if __name__=='__main__':
    with threadpool_limits(limits=2):main()
