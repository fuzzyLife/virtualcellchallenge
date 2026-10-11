"""Freeze V14 identities and bounded experiments before any candidate outcomes."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from v14_context_methods import A,ROOT,NAMES
from fit_v7_response import read_cache
import build_submission as build
SEEDS=[2026101411,2026101427]


def order(name,prefix):return hashlib.sha256((prefix+name).encode()).hexdigest()


def main():
    A.mkdir(parents=True,exist_ok=True);path=A/'plan.json'
    if path.exists():raise FileExistsError('Preserve frozen V14 plan')
    data=read_cache();targets=data['targets'];used=set();parents={}
    for file in sorted((ROOT/'experiments').glob('v*/analysis/**/validation_plan.json')):
        if 'v14_' in str(file):continue
        d=json.loads(file.read_text(encoding='utf-8'))
        for ctx in d.get('contexts',[]):used.update(x['target'] for x in ctx.get('cells',[]))
        parents[str(file.relative_to(ROOT))]=build.sha256(file)
    meta={};eligible={};counts={}
    for family in ['H1','HCT116','HEK293T']:
        folder=ROOT/('reference/v5_h1' if family=='H1' else f'reference/v7_sources/{family}')
        f=pd.read_csv(folder/'selected_cells.csv')
        if family=='H1':f['gene_target']=f.target_gene.replace({'non-targeting':'Non-Targeting'});f['sample']=f.batch.astype(str)
        else:f['sample']=f['sample'].astype(str)
        meta[family]=f;c=f.gene_target.value_counts();counts[family]=c
        eligible[family]=[str(t) for t in targets if c.get(t,0)>=400 and (family=='H1' or t not in used)]
    # Pick fresh confirmation from 400-cell groups before inspecting their effects.
    reserved=set();confirmation={}
    for family in ['HCT116','HEK293T']:
        available=sorted([t for t in eligible[family] if t not in reserved],key=lambda t:order(t,'v14-confirm-'+family))
        if len(available)<2:raise ValueError('Insufficient fresh 400-cell confirmation targets for '+family+': '+str(available))
        confirmation[family]=available[:2];reserved.update(available[:2])
    contexts=[];oracle=[]
    for fi,family in enumerate(['H1','HCT116','HEK293T']):
        f=meta[family];ctrl=np.flatnonzero(f.gene_target.eq('Non-Targeting').to_numpy())
        left=[];right=[]
        for sample,group in f.iloc[ctrl].groupby('sample',sort=True):
            ids=group.index.to_numpy();ids=np.random.default_rng(2026101400+int(order(str(sample),family)[:8],16)).permutation(ids)
            middle=len(ids)//2;left.extend(ids[:middle].tolist());right.extend(ids[middle:].tolist())
        assert len(left)>=400 and len(right)>=400 and not set(left)&set(right)
        dev_available=[str(t) for t in targets if counts[family].get(t,0)>=400 and t not in reserved]
        dev=sorted(dev_available,key=lambda t:(t in used,order(t,'v14-dev-'+family)))[:3]
        assert len(dev)==3
        cells=[]
        for seed in SEEDS:
            rng=np.random.default_rng(seed+fi);template=rng.choice(left,400,replace=False).tolist();reference=rng.choice(right,400,replace=False).tolist()
            for split,names in [('development',dev),('confirmation',confirmation.get(family,[]))]:
                for target in names:
                    ti=int(np.flatnonzero(targets==target)[0]);ids=np.flatnonzero(f.gene_target.eq(target).to_numpy())
                    actual=np.random.default_rng(seed+ti+fi*1000).choice(ids,400,replace=False).tolist()
                    cells.append({'target':target,'seed':seed,'split':split,'actual':actual,'template':template,'reference':reference})
        for target in dev:
            ti=int(np.flatnonzero(targets==target)[0]);ids=np.flatnonzero(f.gene_target.eq(target).to_numpy())
            chosen=np.random.default_rng(2026101401+ti+fi*1000).choice(ids,400,replace=False)
            oracle.append({'family':family,'target':target,'calibration':chosen[:200].tolist(),'truth':chosen[200:].tolist(),'seed':2026101401+ti+fi*1000})
        contexts.append({'family':family,'context':family+'_full400','left':sorted(left),'right':sorted(right),'development_targets':dev,'confirmation_targets':confirmation.get(family,[]),'cells':cells})
        print(family,'development',dev,'confirmation',confirmation.get(family,[]),flush=True)
    plan={'version':14,'seeds':SEEDS,'candidates':NAMES,'contexts':contexts,'oracle_cases':oracle,'confirmation_targets':sorted(reserved),
          'methods':{'penalties':[2.,8.],'context_dimensions':2,'interaction_blend':.5,'maximum_correction_log2':.25,'minimum_source_families':3,'feature_minimum_coverage':.7,'feature_radius_multiplier':1.5},
          'count_policy':'400 genuine unique treated cells; 400 generated predictions; full disjoint control references primary, sampled400 secondary; three development families, two fresh count-confirmation families. Batch proportions pooled at family level; residual batch confounding is diagnosed, not assumed absent.',
          'oracle_policy':'Nine targets. 200 unique held-family treated calibration cells and 200 disjoint truth cells, same 200-cell prediction size/depths. Oracle measured moments and independent actual replicate are diagnostics only, never deployable predictors or model-selection candidates.',
          'selection_policy':'Advance the higher mean family-balanced context2/context8 summary objective to 400-cell development counts, along with fixed shuffled2 and amplitude110 controls. Promote interaction only if it beats V13 and shuffled-context utility and improves at least four available primary metrics, including one of fidelity/reach/jaccard; retain worst-family guards. Freeze before one confirmation. No margin retuned from confirmation.',
          'raw_guards':{'relative_expression_worsening':.01,'fidelity_drop':.01,'reach_drop':.02,'nmae_worsening':.01,'pds_drop':.02,'family_fidelity_drop':.03,'family_reach_drop':.05},
          'count_utility':'mean(pds,fidelity,reach,jaccard,1-nmae,1-expression_ratio); unanchored public utility, not official private score. Missing primary metrics block promotion.',
          'objective_scope':'A new prospective six-metric comparison with a declared 1% expression tolerance; this tolerance is an exploratory practical bound, not an empirically established noninferiority margin. No biological superiority claim without independent confirmation.',
          'summary_scope':'Existing public summary outcomes were inspected previously; omit reserved targets from screen evaluation. Per-target slopes may use training-family responses for those targets as intended known-target context transfer. Does not test unseen-target extrapolation.',
          'development_reuse':{c['family']:[t for t in c['development_targets'] if t in used] for c in contexts},'metadata_limitation':'Only two HCT116 targets have 400 cells and were never count-tested; reserve both for confirmation. HCT116 and H1 development therefore reuse earlier biology. HEK development prioritizes fresh targets. No cell duplication or claim of wholly fresh development.',
          'previous_count_targets':sorted(used),'previous_plans':parents,'source_cache_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),
          'dependencies':{n:build.sha256(ROOT/n) for n in ['v14_context_methods.py','prepare_v14_experiment.py','v13_quality_methods.py','v11_methods.py','v8_methods.py','v7_response_methods.py','v12_cpu_metrics.py']}}
    assert not set(t for c in contexts for t in c['development_targets'])&reserved
    assert not reserved&used
    build.write_json(path,plan)
    build.write_json(A/'split_audit.json',{'passed':True,'development_targets':{c['family']:c['development_targets'] for c in contexts},'confirmation_targets':sorted(reserved),'all_count_targets_400_unique':True,'all_oracle_halves_disjoint':True,'plan_sha256':build.sha256(path)})


if __name__=='__main__':main()
