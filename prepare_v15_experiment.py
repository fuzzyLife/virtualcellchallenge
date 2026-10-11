"""Freeze V15 count identities, candidates and decision policy before outcomes."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from v15_ensemble_methods import ROOT,A,GRID
import build_submission as build

def key(t,p):return hashlib.sha256((p+t).encode()).hexdigest()
def main():
    A.mkdir(parents=True,exist_ok=True);path=A/'plan.json'
    if path.exists():raise FileExistsError('Preserve V15 prospective plan')
    previous=json.loads((ROOT/'experiments/v14_context_interactions/analysis/plan.json').read_text());history=json.loads((ROOT/'experiments/v14_context_interactions/analysis/metadata_correction.json').read_text())
    with np.load(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz',allow_pickle=False) as z:targets=z['targets'].copy()
    eligible={};metadata={}
    for c in previous['contexts']:
        f=c['family'];folder=ROOT/('reference/v5_h1' if f=='H1' else f'reference/v7_sources/{f}');m=pd.read_csv(folder/'selected_cells.csv')
        if f=='H1':m['gene_target']=m.target_gene.replace({'non-targeting':'Non-Targeting'})
        metadata[f]=m;counts=m.gene_target.value_counts();eligible[f]=[str(t) for t in targets if counts.get(t,0)>=400]
    dev={c['family']:list(c['development_targets']) for c in previous['contexts']};protected=set(t for v in dev.values() for t in v)
    reserved={'HCT116':list(previous['contexts'][1]['confirmation_targets'])}
    hct_extra=sorted(set(eligible['HCT116'])-protected-set(reserved['HCT116']),key=lambda t:key(t,'v15devHCT116'))[:1];dev['HCT116']+=hct_extra;protected.update(hct_extra)
    h1_extra=sorted(set(eligible['H1'])-protected-set(reserved['HCT116']),key=lambda t:key(t,'v15devH1'))[:1];dev['H1']+=h1_extra;protected.update(h1_extra)
    reserved['H1']=sorted(set(eligible['H1'])-protected-set(reserved['HCT116']),key=lambda t:key(t,'v15reservedH1'))[:3]
    reserved['HEK293T']=sorted(set(eligible['HEK293T'])-protected-set(sum(reserved.values(),[])),key=lambda t:key(t,'v15reservedHEK293T'))[:4]
    all_reserved=set(sum(reserved.values(),[]));extra=sorted(set(eligible['HEK293T'])-protected-all_reserved,key=lambda t:key(t,'v15devHEK293T'))[:3];dev['HEK293T']+=extra
    assert not all_reserved&set(sum(dev.values(),[]))
    contexts=[]
    for fi,old in enumerate(previous['contexts']):
        f=old['family'];m=metadata[f];cells=[];null=[]
        for seed in previous['seeds']:
            rng=np.random.default_rng(seed+fi);template=rng.choice(old['left'],400,replace=False).tolist();reference=rng.choice(old['right'],400,replace=False).tolist()
            for split,names in [('development',dev[f]),('confirmation',reserved[f])]:
                for target in names:
                    # Identical target-index construction and seeds allow direct V14 baseline replay on shared cases.
                    ti=int(np.flatnonzero(targets==target)[0])
                    ids=np.flatnonzero(m.gene_target.eq(target));actual=np.random.default_rng(seed+ti+fi*1000).choice(ids,400,replace=False).tolist()
                    cells.append({'target':target,'seed':seed,'split':split,'actual':actual,'template':template,'reference':reference})
            ids=np.random.default_rng(seed+fi+15).choice(old['left'],400,replace=False)
            null.append({'seed':seed,'template':ids[:200].tolist(),'actual':ids[200:].tolist()})
        contexts.append({'family':f,'context':f+'_v15','left':old['left'],'right':old['right'],'development_targets':dev[f],'confirmation_targets':reserved[f],'cells':cells,'null_cases':null})
    dependencies=['v15_ensemble_methods.py','v15_population_methods.py','run_v15_experiment.py','prepare_v15_experiment.py','test_v15_ensemble.py','assess_v15_experiment.py','v13_quality_methods.py','v13_robust_methods.py','v11_methods.py','v7_response_methods.py','v8_methods.py','v6_methods.py','v12_cpu_metrics.py','fit_v7_response.py','audit_transfer_validation.py','evaluate_atlas_transfer.py','build_submission.py']
    record={'version':15,'contexts':contexts,'reserved_targets':sorted(all_reserved),'seeds':previous['seeds'],'ensemble_grid':GRID,
        'response_candidates':['v13','student50','ensemble'],'count_candidates':['v13','response','v13_quantile25','response_quantile10','response_quantile25','shuffled_target25'],
        'source_cache_sha256':build.sha256(ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),'previous_plan_sha256':build.sha256(ROOT/'experiments/v14_context_interactions/analysis/plan.json'),
        'history_correction_sha256':build.sha256(ROOT/'experiments/v14_context_interactions/analysis/metadata_correction.json'),
        'dependencies':{n:build.sha256(ROOT/n) for n in dependencies},
        'raw_guards':{'pds_drop':.02,'fidelity_drop':.01,'reach_drop':.02,'jaccard_drop':.002,'nmae_worsening':.01,'relative_expression_worsening':.02,'family_fidelity_drop':.03,'family_reach_drop':.05},
        'selection':'Freeze student50 versus nested ensemble from family-balanced response proxy. Actual-count development tests isolated response/distribution and combined candidates; require >=3 available metrics improve including DE fidelity/reach/overlap; guard PDS and worst-family losses. Use equal metric-improvement ranks among passing candidates; private anchor-scaled utility cannot be estimated. Distribution finalists must beat shuffled-target control on >=2 of fidelity/reach/overlap. Keep exact V13 if none pass.',
        'missing_metrics':'NMAE is averaged only over eligible matched panels; report absent families and do not assert complete six-metric confirmation if any primary panels missing. Near-zero/negative expression denominators are reported separately; expression guard uses H1/HEK293T ratio only if denominator >1e-4, and otherwise paired absolute numerator not ratio with 2% tolerance. This pragmatic new prospective policy is not a retroactive V14 rescue or statistical noninferiority.',
        'confirmation':'One candidate frozen after development, then evaluated once on reserved algorithm-verification target names. The biological sources were previously inspected; this is not blind biological confirmation. No retuning. Package the best supported new candidate as exploratory if confirmation is incomplete, explicitly disclose all failed guards; otherwise retain V13.',
        'null_policy':'Two 200-cell disjoint-template/actual NTC cases per family; compact/quantile10/quantile25, full independent controls. Population diagnostics only, not official perturbation score. Null calls may not rise >10% versus compact for a distribution finalist.',
        'oracle_policy':'No held treatment data enter deployable fitting. Distribution extraction uses only training-source treatments; held-family arrays are excluded by construction.',
        'history_scope':'V14 development and other older public biology are reused. Target names disjoint between current development and reserved verification, not globally unseen. No statistical ranking guarantee.',
        'previously_count_tested_names':history['complete_previous_target_exclusions']}
    build.write_json(path,record)
    print('V15 development:',dev,'reserved:',reserved,flush=True)
if __name__=='__main__':main()
