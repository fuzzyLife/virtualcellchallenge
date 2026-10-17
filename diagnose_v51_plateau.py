"""Post-V50 descriptive attribution. Reuses spent public outcomes; never selects a release."""
from pathlib import Path
import csv, hashlib, json, platform
import numpy as np
from scipy import sparse, stats
from threadpoolctl import threadpool_limits
from v12_cpu_metrics import de

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis/v51_plateau'
A=ROOT/'experiments/v50_response_prior/analysis'
CACHE=ROOT/'experiments/v37_component_ensemble/analysis/cache'
INPUTS={}

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(2**20),b''): h.update(block)
    return h.hexdigest()

def read(p):
    INPUTS[str(p.relative_to(ROOT))]=sha(p)
    return json.loads(p.read_text(encoding='utf-8-sig'))

def arrays(p):
    INPUTS[str(p.relative_to(ROOT))]=sha(p)
    with np.load(p,allow_pickle=False) as z:return {k:z[k] for k in z.files}

def pds(pred,truth,ctrl):
    a,b=pred-ctrl,truth-ctrl
    norm=np.linalg.norm(a,axis=1)[:,None]*np.linalg.norm(b,axis=1)[None,:]
    dist=1-np.divide(a@b.T,norm,out=np.zeros((len(a),len(b))),where=norm>0)
    diagonal=np.diag(dist)
    ranks=(dist<diagonal[:,None]).sum(1)+.5*((dist==diagonal[:,None]).sum(1)-1)
    return float(np.mean(1-ranks/(len(a)-1)))

def finite_mean(v):
    v=[float(x) for x in v if x is not None and np.isfinite(x)]
    return float(np.mean(v)) if v else None

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    previous=json.loads((OUT/'diagnosis.json').read_text()) if (OUT/'diagnosis.json').exists() else None
    published=read(A/'published_result.json')
    result={'scope':'Descriptive attribution after V50 publication; no blind validation, no release selection, D3 unopened',
        'python':platform.python_version(),'script_sha256':sha(Path(__file__)),'published_delta':published['comparison']['delta_vs_v24']}
    plan=read(ROOT/'experiments/v34_residual_flow/analysis/plan.json')
    stages={s:read(A/(s+'_results.json')) for s in ['development','confirmation']}
    null=read(A/'null_results.json')
    result['nulls']=[{k:r[k] for k in ['family','null_DE_calls','actual_DE_calls']} for r in null['records'] if r['arm']=='v24']
    source=ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'
    with np.load(source,allow_pickle=False) as z:genes=z['genes'].copy();targets=z['targets'].copy()
    # Pin the axis source's existing provenance, without repeatedly reading its large arrays.
    result['axis_source']='experiments/v7_multisource_adaptive/analysis/source_cache.npz'
    read(ROOT/'experiments/v7_multisource_adaptive/analysis/source_provenance.json')
    scored=~np.isin(genes,targets);cases=[];panels=[];raw_null=[];rng=np.random.default_rng(51051)
    for ctx in plan['contexts']:
        family=ctx['family'];f=CACHE/family
        effects_path=A/'predictions'/(family+'.npz')
        effects=arrays(effects_path)
        assert INPUTS[str(effects_path.relative_to(ROOT))]==read(effects_path.with_suffix('.json'))['sha256']
        # Read only the control probabilities, not the many unused 300-target banks.
        INPUTS[str((f/'effects.npz').relative_to(ROOT))]=sha(f/'effects.npz')
        with np.load(f/'effects.npz',allow_pickle=False) as z:cp=z['cp'].copy()
        ref_path=f/'reference.npz';INPUTS[str(ref_path.relative_to(ROOT))]=sha(ref_path)
        reference=sparse.load_npz(ref_path);depth=reference.sum(1).A1
        refcp=reference.multiply((1e6/depth)[:,None]).mean(0).A1
        records=[];vectors={a:[] for a in ['v24','randomized_prior']};truth_profiles=[];control=None
        for i,cell in enumerate(ctx['cells']):
            stage='development' if i<4 else 'confirmation';stem=f'{family}__{i}__0__'
            own=int(np.flatnonzero(genes==cell['target'])[0]);base=None
            for arm in vectors:
                r=read(A/stage/(stem+arm+'.json'));array_path=A/stage/(stem+arm+'.npz');z=arrays(array_path)
                assert INPUTS[str(array_path.relative_to(ROOT))]==r['arrays_sha256']
                vectors[arm].append(z['pred'][scored]);control=z['control'][scored]
                if arm=='v24':truth_profiles.append(z['truth'][scored]);base=z
                valid=z['truth_eligible'].astype(bool)&np.isfinite(z['truth_lfc']);valid[own]=False
                real=valid&(z['truth_q']<.05);called=valid&(z['pred_q']<.05)
                factor=np.exp2(np.clip(effects[arm][i].astype(float),-3,3));factor[own]=.1
                prob=cp*factor;prob/=prob.sum()
                intended=np.log2((prob*1e6+1e-9)/(refcp+1e-9))
                truth=z['truth_lfc'];pred=z['pred_lfc'];mask=real if real.any() else valid
                denom=float(np.sum(truth[mask]**2))
                rec=dict(family=family,study=r['study'],target=cell['target'],case_index=i,arm=arm,
                    true_DE=int(real.sum()),predicted_DE=int(called.sum()),
                    DE_precision=float((real&called).sum()/max(called.sum(),1)),
                    DE_recall=float((real&called).sum()/max(real.sum(),1)),
                    truth_DE_sign_accuracy=float(np.mean(np.sign(pred[real])==np.sign(truth[real]))) if real.any() else None,
                    intended_DE_sign_accuracy=float(np.mean(np.sign(intended[real])==np.sign(truth[real]))) if real.any() else None,
                    squared_error_vs_zero_on_mask=float(np.sum((pred[mask]-truth[mask])**2)/denom) if denom>0 else None,
                    intended_squared_error_vs_zero_on_mask=float(np.sum((intended[mask]-truth[mask])**2)/denom) if denom>0 else None,
                    emitted_minus_intended_LFC_RMSE=float(np.sqrt(np.mean((pred[valid]-intended[valid])**2))),
                    response_change_relative_L2=float(np.linalg.norm(effects[arm][i]-effects['v24'][i])/max(np.linalg.norm(effects['v24'][i]),1e-12)),
                    raw_fidelity=r['fidelity'],raw_reach=r['reach'])
                cases.append(rec)
            records.append(base)
        # Recompute PDS on all available cases and compare the original small panels.
        truth=np.array(truth_profiles)
        labels=np.array([cell['target'] for cell in ctx['cells']]);unique=np.unique(labels)
        collapse=lambda values:np.array([np.mean(values[labels==label],axis=0) for label in unique])
        row={'family':family,'available_case_populations':len(truth),'unique_targets':len(unique),'repeated_target_populations':int(len(truth)-len(unique)),'reference_cells':len(depth),'template_cells':int(sparse.load_npz(f/'0_template.npz').shape[0])}
        for arm,values in vectors.items():
            pred=np.array(values);row[arm+'_all_case_population_PDS_not_unique_target_metric']=pds(pred,truth,control)
            row[arm+'_unique_target_PDS']=pds(collapse(pred),collapse(truth),control)
            row[arm+'_development_PDS']=pds(pred[:4],truth[:4],control)
            row[arm+'_confirmation_PDS']=pds(pred[4:],truth[4:],control)
            collapsed_pred,collapsed_truth=collapse(pred),collapse(truth)
            placebo=[pds(collapsed_pred[rng.permutation(len(unique))],collapsed_truth,control) for _ in range(100)]
            row[arm+'_label_shuffle_PDS_95pct']=[float(x) for x in np.quantile(placebo,[.025,.975])]
        panels.append(row)
        # Existing cached template and reference are disjoint real NTC pools.
        template_path=f/'0_template.npz';INPUTS[str(template_path.relative_to(ROOT))]=sha(template_path)
        key_ref=str(ref_path.relative_to(ROOT));key_template=str(template_path.relative_to(ROOT))
        can_reuse=previous and all(previous['inputs_sha256'].get(k)==INPUTS[k] for k in [key_ref,key_template])
        if can_reuse:
            raw_null.append(next(x for x in previous['cached_real_control_checks'] if x['family']==family))
        else:
            template=sparse.load_npz(template_path);dd=de(template,reference,genes)
            raw_null.append(dict(family=family,real_template_vs_reference_DE=int(np.sum(dd['eligible']&(dd['q']<.05))),
                scope='Case-0 real NTC template; not the separate saved synthetic-null split'))
        print('Audited',family,'targets',len(truth),'real NTC calls',raw_null[-1]['real_template_vs_reference_DE'],flush=True)
    result.update(cases=cases,panel_sensitivity=panels,cached_real_control_checks=raw_null,
        unique_context_target_pairs=len({(c['family'],c['target']) for c in cases}),
        unique_context_case_populations=len({(c['family'],c['case_index']) for c in cases}),
        real_control_cache_policy='Reuse completed DE check only if exact reference and template SHA256 match previous run; otherwise rerun.')
    summary=[]
    for family in [c['family'] for c in plan['contexts']]:
        for arm in ['v24','randomized_prior']:
            rr=[r for r in cases if r['family']==family and r['arm']==arm]
            summary.append(dict(family=family,arm=arm,cases=len(rr),
                median_true_DE=float(np.median([r['true_DE'] for r in rr])),
                median_predicted_DE=float(np.median([r['predicted_DE'] for r in rr])),
                cases_true_DE_below_10=sum(r['true_DE']<10 for r in rr),
                **{k:finite_mean([r[k] for r in rr]) for k in ['DE_precision','DE_recall','truth_DE_sign_accuracy','intended_DE_sign_accuracy','squared_error_vs_zero_on_mask','intended_squared_error_vs_zero_on_mask','emitted_minus_intended_LFC_RMSE','response_change_relative_L2']}))
    result['context_summary']=summary
    result['expression_panels']=[{k:p[k] for k in ['family','study','seed_offset','arm','common_expression_numerator','common_expression_denominator','common_expression_ratio']} for stage in stages.values() for p in stage['panels']]
    # Official score contributions vs the actual public relative-surrogate definition.
    result['confirmation_changes']={}
    for offset in [0,380000]:
        pp=[p for p in stages['confirmation']['panels'] if p['seed_offset']==offset]
        row=[]
        for family in [c['family'] for c in plan['contexts']]:
            base=next(p for p in pp if p['family']==family and p['arm']=='v24')
            alt=next(p for p in pp if p['family']==family and p['arm']=='randomized_prior')
            row.append(dict(family=family,study=base['study'],reach_delta=None if base['reach'] is None else alt['reach']-base['reach'],
                expression_ratio_delta=None if base['common_expression_ratio'] is None else alt['common_expression_ratio']-base['common_expression_ratio']))
        result['confirmation_changes'][str(offset)]=row
    result['inputs_sha256']=INPUTS
    with (OUT/'case_diagnostics.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(cases[0]));writer.writeheader();writer.writerows(cases)
    (OUT/'diagnosis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'contexts':summary,'confirmation_changes':result['confirmation_changes'],'panels':panels},indent=2),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=2):main()
