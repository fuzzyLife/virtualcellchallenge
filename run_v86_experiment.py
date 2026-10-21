"""Freeze and evaluate V86 without selecting another strength from recipient outcomes."""
import argparse,gc
import numpy as np
from threadpoolctl import threadpool_limits
import v86_joint_posterior as m
from v57_pie import read,write,sha
from evaluate_v77 import profiles
from audit_transfer_validation import pds
from run_v35_v36 import contexts,key

def freeze():
    assert not (m.A/'plan.json').exists();names=['v86_joint_posterior.py','run_v86_experiment.py','run_v86_counts.py','test_v86_posterior.py','v15_ensemble_methods.py','v62_conserved_response.py','v63_identity_guard.py','v66_assay_strength.py','v16_methods.py','v12_cpu_metrics.py','v8_methods.py']
    write(m.A/'plan.json',dict(version=86,primary='interval25',hypothesis='Learn a mixture of shared and context-dependent multivariate response priors, condition on each available measured source response and uncertainty, and project V66 toward a posterior interval for an unseen context.',novelty='First covariance-mixture GP conditional response model in this project. Earlier scalar pooling/shrinkage averaged effects; V84 fitted a mean response operator on different essential targets. This model learns sharing patterns across the actual panel response coordinates and propagates conditional uncertainty. Independent mash-inspired adaptation, not mashr or a published CRISPR predictor reproduction.',sources='Legacy measured K562,CD4,H1,HCT116,HEK293T responses and masks; same bounded90%-dose cube as V66. Convert to log1p(50000*proportion) response with derivative-propagated sampling variance floor0.0025. Every panel gene excluded as an output.',kernel='Fixed common expressed nonpanel rank axis (>5CPM and finite in every allowed source), centered unit rank vectors. Source kernel0.8*exp(-(1-inner_product)/length)+0.2I; length0.05,0.2,1,4 or all-ones shared base. Query ranks use finite observed values with missing embedding coordinates zero (unknown centered-rank feature, not imputed RNA); >=90% measured required. Positive-semidefinite RBF embedding; query-specific0.2 variance remains unobserved.',mixture='Null component plus four signal SDs0.025,0.1,0.4,1.0 times five kernels =21 covariance components. Fit mixture weights to20000 deterministic source target/output coordinates with>=2 valid sources, using exact small Gaussian marginal likelihood,100 EM updates or max weight change<1e-7,0.5 component pseudocount. Gene coordinates are dependent working-likelihood units, not independent biological replicates.',boundaries='Whole H1,Orion(HCT116+HEK293T),or CD4(all states/donors) excluded from every training response, covariance fitting and conditioning. Query NTCs only for embedding, original V66 baseline, and inversion. D3 is reused after V84, not blind.',arms='Primary interval25: project V66 assay response into central80% posterior interval, invert to LFC and apply quarter of bounded+/-0.5LFC correction. posterior25 uses posterior mean instead. shared25 independently fits mixture weights restricted to shared-kernel/null components and uses interval correction. interval25_shuffled refits the full prior with independently shifted source target labels and conditions on those shifted responses.',guard='Original80% V66 identity-margin guard. Own/panel effects unchanged. Missing measured source entries abstain; entirely unsupported coordinates and insufficient query rank coverage return baseline exactly.',development='All655 reused public pairs. Primary must lower equal-study profile error>=0.5%, worst-study<=1.03, PDS loss<=0.005, and beat both shared-only and refitted target-shuffled control profile error. All arms reported; no recipient-based selection or hyperparameter tuning.',counts='Always test fixed primary, shared and refitted-shuffled control on original14-target D3 cohort, seed0, unchanged residual10 and matched400-cell templates/libraries/reference. Reuse14 previously validated V84 V66 baseline populations only after exact predictor/CPM/cohort/code equality checks. Emit42 new populations, independently replay all56. No duplicated cells or treated inputs for inference.',count_gates='Profile error must not increase; PDS loss<=.005,FID loss<=.002,reach loss<=.005,JAC loss<=.0005,NMAE<=1.01baseline; >=3 endpoints improve and profile beats both controls. A second unchanged seed620000 only if development and first count gates all pass. No package otherwise.',automatic_upload=False,literature=[dict(url='https://www.nature.com/articles/s41588-018-0268-8',idea='Multivariate adaptive shrinkage learns heterogeneous effect sharing; its genomic association benchmarks do not establish VCC performance.'),dict(url='https://stephenslab.github.io/mashr/',idea='Author method documentation; software not executed.')],code_sha256={n:sha(m.ROOT/n) for n in names},source_cache_sha256=sha(m.ROOT/'experiments/v7_multisource_adaptive/analysis/source_cache.npz')))
    print('V86 protocol frozen',flush=True)

def prepare():
    m.check();dest=m.A/'predictions';dest.mkdir(exist_ok=True)
    for ctx in contexts():
        name=ctx['family'];p=m.ROOT/'experiments/v66_assay_strength/analysis/predictions'/(name+'.npz');assert sha(p)==read(p.with_suffix('.json'))['sha256']
        with np.load(p,allow_pickle=False) as z:query=z['query'].copy();bulk=z['bulk'].copy();base=z['assay75'].copy()
        pred,tr=m.predict(query,base,bulk,key(name));out=dest/(name+'.npz');assert not out.exists();np.savez_compressed(out,query=query,bulk=bulk,v66=base,**pred);write(out.with_suffix('.json'),dict(sha256=sha(out),source_V66_sha256=sha(p),outer=key(name),trace=tr,script_sha256=sha(__file__)));print('V86 predicted',name,tr['interval25']['changed_targets'],flush=True);gc.collect()
    p=m.ROOT/'experiments/v84_intervention_operator/analysis/donor_predictions.npz';assert sha(p)==read(p.with_suffix('.json'))['sha256']
    with np.load(p,allow_pickle=False) as z:q=z['query'].copy();cp=z['cp'].copy();bulk=z['bulk'].copy();base=z['v66'].copy();genes=z['genes'].copy();targets=z['targets'].copy()
    pred,tr=m.predict(q,base,bulk,'CD4');out=dest/'CD4_D3_Rest.npz';assert not out.exists();np.savez_compressed(out,query=q,cp=cp,bulk=bulk,v66=base,genes=genes,targets=targets,**pred);write(out.with_suffix('.json'),dict(sha256=sha(out),source_V84_sha256=sha(p),outer='CD4',trace=tr,script_sha256=sha(__file__)));print('V86 predicted D3',flush=True)

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
    primary=summary['interval25'];checks=dict(profile=primary['profile_ratio']<=.995,worst_study=max(v['profile_ratio'] for v in studies['interval25'].values())<=1.03,pds=primary['pds_change']>=-.005,beats_shared=primary['profile_ratio']<summary['shared25']['profile_ratio'],beats_shuffle=primary['profile_ratio']<summary['interval25_shuffled']['profile_ratio'])
    write(m.A/'mean_results.json',dict(summary=summary,studies=studies,records=records,checks=checks,passes=all(checks.values()),primary='interval25',pairs=sum(r['pairs'] for r in records),scope='655 reused public target/context pairs. Equal-study means, no private score or independent biology claim.'))
    print('V86 mean results',summary,'checks',checks,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','fit','prepare','means']);p.add_argument('--outer',choices=['H1','Orion','CD4']);p.add_argument('--shuffled',action='store_true');a=p.parse_args()
    with threadpool_limits(limits=2):{'freeze':freeze,'fit':lambda:m.fit(a.outer,a.shuffled),'prepare':prepare,'means':means}[a.stage]()
