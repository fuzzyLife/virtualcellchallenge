"""Descriptive leave-case sensitivity of the V50 confirmation; no revised selection."""
from pathlib import Path
import hashlib,json,math
import numpy as np
ROOT=Path(__file__).resolve().parent
A=ROOT/'experiments/v50_response_prior/analysis';OUT=ROOT/'analysis/v51_plateau'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def aggregate(records,arm,metric,exclude=None):
    rows=[r for r in records if r['arm']==arm and (r['family'],r['target'])!=exclude]
    studies=[]
    for study in sorted({r['study'] for r in rows}):
        families=[]
        for family in sorted({r['family'] for r in rows if r['study']==study}):
            seeds=[]
            for seed in sorted({r['seed_offset'] for r in rows if r['family']==family}):
                v=[r[metric] for r in rows if r['family']==family and r['seed_offset']==seed and r.get(metric) is not None and np.isfinite(r[metric])]
                if v:seeds.append(float(np.mean(v)))
            if seeds:families.append(float(np.mean(seeds)))
        if families:studies.append(float(np.mean(families)))
    return float(np.mean(studies)) if studies else None

def main():
    p=A/'confirmation_results.json';r=json.loads(p.read_text());rows=r['records']
    base=aggregate(rows,'v24','reach');alt=aggregate(rows,'randomized_prior','reach')
    reference=json.loads((A/'confirmation_summary.json').read_text())
    assert np.isclose(base,reference['v24']['reach'],rtol=0,atol=1e-12)
    assert np.isclose(alt,reference['randomized_prior']['reach'],rtol=0,atol=1e-12)
    influence=[]
    for family,target in sorted({(r['family'],r['target']) for r in rows}):
        b=aggregate(rows,'v24','reach',(family,target));a=aggregate(rows,'randomized_prior','reach',(family,target))
        cc=[x for x in rows if x['family']==family and x['target']==target]
        influence.append(dict(family=family,target=target,reference_DE=min(x['reference_calls'] for x in cc),
            full_reach_delta=alt-base,leave_case_out_reach_delta=a-b,
            removed_case_contribution_to_delta=(alt-base)-(a-b)))
    influence.sort(key=lambda x:-abs(x['removed_case_contribution_to_delta']))
    diagnosis=json.loads((OUT/'diagnosis.json').read_text())
    power=[]
    for panel in diagnosis['panel_sensitivity']:
        m,n=panel['template_cells'],panel['reference_cells']
        power.append(dict(family=panel['family'],generated_cells=m,public_reference_cells=n,official_generated_cells=400,official_control_cells=18400,
            asymptotic_z_multiplier_if_same_AUC_and_tie_fraction=math.sqrt((1/m+1/n)/(1/400+1/18400))))
    ep=[x for x in r['panels'] if x['arm']=='v24']
    hct=np.mean([x['common_expression_ratio'] for x in ep if x['family']=='HCT116'])
    total=reference['v24']['common_expression_ratio']
    result=dict(scope='Post-publication descriptive sensitivity; no revised guard, release, or prediction of official score',
        full_public_raw_reach=dict(v24=base,v50=alt,delta=alt-base),leave_case_out=influence,
        HCT116_fraction_of_equal_study_public_expression_ratio=float(hct/6/total),
        expression_ratio_aggregation='H1 weight 1/3, HCT116+HEK293T each 1/6, two CD4 contexts each 1/6',
        idealized_control_power=dict(contexts=power,limitation='First-order calculation from actual cached template/reference sizes, not a measured private effect. Equal AUC and tie fractions assumed; official group size 400, controls 18400 from manifest.'),
        known_experiment=dict(family='CD4_Rest_D1',target='MTHFSD',reference_DE=3,
            seed_0_reach_v24=0,seed_0_reach_v50=1,second_seed_jump_absent=True),
        inputs_sha256={str(x.relative_to(ROOT)):digest(x) for x in [p,A/'confirmation_summary.json',ROOT/'manifest.json',OUT/'diagnosis.json']},
        script_sha256=digest(Path(__file__)))
    OUT.mkdir(exist_ok=True)
    (OUT/'selection_sensitivity.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['inputs_sha256','leave_case_out']},indent=2));print(json.dumps(influence[:3],indent=2))
if __name__=='__main__':main()
