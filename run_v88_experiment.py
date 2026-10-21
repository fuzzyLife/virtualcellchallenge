"""Freeze and evaluate V88 without selecting another strength from recipient outcomes."""
import argparse,gc
import numpy as np
from threadpoolctl import threadpool_limits
import v88_mixture_quantiles as m
from v57_pie import read,write,sha
from evaluate_v77 import profiles
from audit_transfer_validation import pds
from run_v35_v36 import contexts,key

def freeze():
    assert not (m.A/'plan.json').exists();m.A.mkdir(parents=True,exist_ok=True)
    names=['v88_mixture_quantiles.py','run_v88_experiment.py','run_v88_counts.py','test_v88_mixture.py','audit_v88_models.py','v86_joint_posterior.py','v15_ensemble_methods.py','v62_conserved_response.py','v63_identity_guard.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
    write(m.A/'plan.json',dict(version=88,primary='exact125',hypothesis='Repair V86 source-prior underoptimization and normal-moment posterior interval approximation, then test a conservative one-eighth correction to V66. The old V86 interval25 donor result improved5of6 metrics but failed broad/shuffled criteria; its scientific records remain unchanged.',prior='Same21 covariance components, source bank, masks, variance floor, rank kernel and20000deterministic fitting coordinates as V86. Strictly concave mixture MAP sum log(Lw)+0.5sum logw, solved with positive-simplex damped Newton to relative KKT<1e-8. Shared-only control independently solved. Old training samples exactly matched; source objective compared against old100iteration EM. Production fit uses all allowed sources.',interval='Use numerical generalized-inverse mixture CDF at0.1/0.9, retaining the null point mass,32bisections with12SD finite bounds. Project V66 response onto the resulting mixture interval. This is a plug-in posterior-model interval, not demonstrated biological coverage. V86 normal moments are no longer used for bounds.',arms='Sole release candidate exact125: one-eighth of LFC difference capped+/-0.5, then unchanged80%V66 identity-margin guard. Shared125 uses independently converged shared-only prior. Exact125_shuffled independently refits full prior after source-specific target shifts and uses shifted responses at inference. Controls are never promoted. No response-specific strength search.',boundaries='Whole H1,Orion(both cell lines),or CD4(all donors/states) excluded before fitting and conditioning. Query NTCs only. No output panel genes proposed; guard closure can change non-own factors. Own effects remain0; legacy emitter applies knockdown once. Existing studies and D3 are reused, not blind.',broad_test='All655 reused public pairs. Strong support requires equal-study profile ratio<=.995,worststudy<=1.03,PDSloss<=.005,and lower profile error than both controls. All fixed outputs retained. No candidate or strength selected from these outcomes.',counts='Always evaluate sole primary and both controls on original14target400cellD3cohort and both offsets0/620000, unchanged V66 residual10. Reuse previously audited identical V84 V66 baselines only after effect/queryCP/bulk/axis/cohort/code equality.42new populations/seed;14reused baselines/seed. Full raw-p/BH/direction/profile/PDS replay.',strict_count='Original V86 noninferiority limits and >=3 improved endpoints; candidate must beat both controls on profile error on both technical seeds. Both seeds share biology. Missing metrics never counted as gains.',experimental_release='User explicitly requests another upload test. The sole candidate may be packaged as EXPERIMENTAL even if strong support fails, only if broad profile<=1.01,worststudy<=1.03,PDSloss<=.01 and BOTH donor seeds have profile<=1.005baseline,PDSloss<=.005,FIDloss<=.002,JAC loss<=.0005,NMAE<=1.01baseline and>=2of6 endpoints improve. All failures disclosed; no better-score claim. Severe regressions block this release. Selection is frozen before outcomes; no postfailure alternate control/strength promotion.',production='Fit converged all-source prior, transform the immutable V66 A/B/C effects using only original controls, keep exact V66 template identities, libraries and seeds. Generate360000integer-countcells on18533orderedgenes; unchanged residual10, official prep and full payload/provenance/integer/library and repeated full-group regeneration checks. Nothing uploaded or pushed automatically.',code_sha256={n:sha(m.ROOT/n) for n in names},source_cache_sha256=sha(m.ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz'),known_results='V66 overall.15088218755; V86 public ratio.99626584,PDSchange-.00053622; V86 D3seed0 profile-.2907%,5of6gains but shuffled better. V87 primary-.151%profile,PDS-.010996,failed. These outcomes motivated this prospectively frozen follow-up.'))
    print('V88 protocol frozen',flush=True)

def prepare():
    m.check();dest=m.A/'predictions';dest.mkdir(exist_ok=True)
    for ctx in contexts():
        name=ctx['family'];p=m.ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
        with np.load(p,allow_pickle=False) as z:query=z['query'].copy();bulk=z['bulk'].copy();base=z['assay75'].copy()
        pred,tr=m.predict(query,base,bulk,key(name));out=dest/(name+'.npz');assert not out.exists();np.savez_compressed(out,query=query,bulk=bulk,v66=base,**pred);write(out.with_suffix('.json'),dict(sha256=sha(out),source_V66_sha256=sha(p),outer=key(name),trace=tr,script_sha256=sha(__file__)));print('V88 predicted',name,tr['exact125']['changed_targets'],flush=True);gc.collect()
    p=m.ROOT/'experiments/v84_intervention_operator/analysis/donor_predictions.npz';assert sha(p)==read(p.with_suffix('.json'))['sha256']
    with np.load(p,allow_pickle=False) as z:q=z['query'].copy();cp=z['cp'].copy();bulk=z['bulk'].copy();base=z['v66'].copy();genes=z['genes'].copy();targets=z['targets'].copy()
    pred,tr=m.predict(q,base,bulk,'CD4');out=dest/'CD4_D3_Rest.npz';assert not out.exists();np.savez_compressed(out,query=q,cp=cp,bulk=bulk,v66=base,genes=genes,targets=targets,**pred);write(out.with_suffix('.json'),dict(sha256=sha(out),source_V84_sha256=sha(p),outer='CD4',trace=tr,script_sha256=sha(__file__)));print('V88 predicted D3',flush=True)

def means():
    m.check();d=m.data();genes=d['genes'];own=np.array([list(genes).index(t) for t in d['targets']]);arms=('v66',)+m.ARMS;records=[]
    for ctx in contexts():
        name=ctx['family'];p=m.A/'predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
        with np.load(p,allow_pickle=False) as z:pred={a:profiles(z[a],z['bulk'],own) for a in arms}
        p=m.ROOT/'experiments/v77_common_response/analysis/truth_bank'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
        with np.load(p,allow_pickle=False) as z:bank={k:z[k].copy() for k in z.files}
        scores={};axis=bank['measured']&~np.isin(genes,d['targets'])
        for arm in arms:
            errors=[]
            for i,t in enumerate(bank['target_indices']):
                mask=bank['measured'].copy();mask[own[t]]=False;errors.append(float(np.mean((pred[arm][t,mask]-bank['truth'][i,mask])**2)))
            scores[arm]=dict(profile_MSE=float(np.mean(errors)),pds=pds(pred[arm][bank['target_indices']][:,axis],bank['truth'][:,axis],bank['control'][axis]))
        records.append(dict(source=name,study=key(name),pairs=len(bank['target_indices']),scores=scores))
    summary={};studies={}
    for arm in m.ARMS:
        studies[arm]={}
        for group in ['H1','Orion','CD4']:
            rows=[r for r in records if r['study']==group];studies[arm][group]=dict(profile_ratio=float(np.mean([r['scores'][arm]['profile_MSE']/r['scores']['v66']['profile_MSE'] for r in rows])),pds_change=float(np.mean([r['scores'][arm]['pds']-r['scores']['v66']['pds'] for r in rows])))
        summary[arm]={k:float(np.mean([v[k] for v in studies[arm].values()])) for k in ['profile_ratio','pds_change']}
    primary=summary['exact125'];checks=dict(profile=primary['profile_ratio']<=.995,worst_study=max(v['profile_ratio'] for v in studies['exact125'].values())<=1.03,pds=primary['pds_change']>=-.005,beats_shared=primary['profile_ratio']<summary['shared125']['profile_ratio'],beats_shuffle=primary['profile_ratio']<summary['exact125_shuffled']['profile_ratio'])
    write(m.A/'mean_results.json',dict(summary=summary,studies=studies,records=records,checks=checks,passes=all(checks.values()),primary='exact125',pairs=sum(r['pairs'] for r in records),scope='655 reused public target/context pairs. Equal-study means, no private score or independent biology claim.'))
    print('V88 mean results',summary,'checks',checks,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','fit','prepare','means']);p.add_argument('--outer',choices=['H1','Orion','CD4','production']);p.add_argument('--shuffled',action='store_true');a=p.parse_args()
    with threadpool_limits(limits=2):{'freeze':freeze,'fit':lambda:m.fit(a.outer,a.shuffled),'prepare':prepare,'means':means}[a.stage]()
