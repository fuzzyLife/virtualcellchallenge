"""Document V67's fixed network experiment and preserve its scientific receipts."""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import v67_intervention_network as m
from v57_pie import read,write,sha


def main():
    m.check();r=read(m.A/'mean_results.json');assert read(m.A/'results_audit.json')['passed']
    base=r['summary']['network']['v63'];s=base['all'];p66=r['summary']['network']['v66']['all']
    supported=r['release_supported'];assert not supported,'Passing candidate requires count evaluation before release documentation.'
    rows=[];traces=[]
    for path in sorted((m.A/'broad').glob('*.json')):
        z=read(path)
        for arm in ('network','network_shuffled','gain_only'):
            for ref in ('v63','v66'):
                v=[a['profile_MSE'] for a in z['records'] if a['arm']==arm]
                b=[a['profile_MSE'] for a in z['records'] if a['arm']==ref]
                pp=next(a['pds'] for a in z['panels'] if a['arm']==arm and a['group']=='all')
                pb=next(a['pds'] for a in z['panels'] if a['arm']==ref and a['group']=='all')
                rows.append(dict(source=path.stem,arm=arm,baseline=ref,targets=z['targets'],profile_ratio=float(np.mean(v)/np.mean(b)),pds_change=pp-pb))
        trace=read(m.A/'predictions'/(path.stem+'.json'))['trace']
        traces.append(dict(source=path.stem,mean_absolute_correction=trace['mean_absolute_accepted_correction'],maximum_raw_correction=trace['maximum_raw_correction'],source_models=trace['sources']))
    with (m.A/'context_comparison.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    sources=[c['family'] for c in __import__('run_v35_v36').contexts()]
    colors={'network':'#246a9e','network_shuffled':'#af6b28','gain_only':'#777777'}
    labels={'network':'Network primary','network_shuffled':'Shuffled mediator activity','gain_only':'Amplitude-only control'}
    for ax,metric,scale,title in [(axes[0],'profile_ratio',100,'Intended profile error change (%)'),(axes[1],'pds_change',1,'Raw discrimination change')]:
        for j,arm in enumerate(m.ARMS):
            vals=[next(v[metric] for v in rows if v['source']==src and v['arm']==arm and v['baseline']=='v63') for src in sources]
            if metric=='profile_ratio':vals=(np.array(vals)-1)*scale
            xx=np.arange(len(sources))+(j-1)*.18
            ax.vlines(xx,0,vals,color=colors[arm],alpha=.6)
            ax.scatter(xx,vals,color=colors[arm],label=labels[arm],s=26)
        ax.axhline(0,color='black',lw=.8);ax.set_xticks(range(len(sources)),['H1','HCT116','HEK293T','CD4 rest','CD4 stim'],rotation=25)
        ax.set_title(title);ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Lower is better; V63 = 0%');axes[1].set_ylabel('Higher is better; V63 = 0')
    axes[0].legend(fontsize=8,loc='best');fig.suptitle('V67: fixed intervention-network test on 655 public target/context pairs')
    figure=m.A/'network_tradeoff.png';fig.savefig(figure,dpi=160);plt.close(fig)
    table=['| Study | Profile error change vs V63 | PDS change vs V63 |','| --- | ---: | ---: |']
    table += [f"| {study} | {(v['profile_MSE_ratio']-1):+.3%} | {v['pds_change']:+.6f} |" for study,v in s['studies'].items()]
    fails=[k for k,v in r['checks'].items() if not v]
    text=f'''# V67 results: intervention-network rewiring

The fixed network primary **did not pass its release criteria**. Across all **655 target/context pairs**, equal-study intended profile error changes **{s['profile_MSE_ratio']-1:+.4%} versus V63** and **{p66['profile_MSE_ratio']-1:+.4%} versus V66**. Raw discrimination changes **{s['pds_change']:+.7f} versus V63**. No V67 submission cells or upload file were generated; V66 remains the highest published score, with V63 the near-tied reference.

## What was actually different

The model derives a bounded interaction operator from measured responses of the 300 target genes to each other's knockdowns. It re-solves the feedback system after changing outgoing mediator availability using recipient controls, then transports the change in source response onto V63. It is neither a STRING embedding predictor nor a control-covariance transport. The operator and saturation rule are predictive assumptions, not identified causal biology. [Protocol, formulas, limitations and references](V67_STRATEGY.md).

## Broad results

{chr(10).join(table)}

The whole-target diagnostic development profile ratio is **{base['development']['profile_MSE_ratio']:.9f}**, confirmation **{base['confirmation']['profile_MSE_ratio']:.9f}**, and canonical-small target-subset ratio **{base['canonical_small']['profile_MSE_ratio']:.9f}**, all versus V63. The canonical-small subset uses the same full treatment populations as this broad test; it does not recreate the old 200/400-cell samples. Target partitions reuse previously studied public biology and are not an independent blind test.

Matched control profile ratios versus V63: **shuffled mediator availability {r['summary']['network_shuffled']['v63']['all']['profile_MSE_ratio']:.9f}**, **amplitude-only projection {r['summary']['gain_only']['v63']['all']['profile_MSE_ratio']:.9f}**. These ratios are errors, with lower better. None is a private anchored accuracy. The primary and controls were fixed before outcomes and only the primary could qualify.

Failed fixed checks: {', '.join(fails)}. All 12 check outcomes and all subgroup metrics are in `experiments/v67_intervention_network/analysis/mean_results.json`. Count generation was conditional on passing this stage and was not run. No claim is made about V67's hosted score or its DE fidelity, reach or overlap.

![Fixed network comparison](experiments/v67_intervention_network/analysis/network_tradeoff.png)

## What this experiment establishes

This particular baseline-conditioned feedback construction does not meet the improvement needed to justify another upload. It does not rule out all regulatory-network models: only 300 intervention nodes are represented; source effects are noisy, many source measurements are absent, the 90% interventions are not small perturbations, and availability of a transcript is not equivalent to regulatory activity. The stability cap and valid-output requirement also make accepted corrections small. Per-context source coverage and correction magnitudes are saved in `network_diagnostics.json`.

Do not increase feedback strength or choose a control arm using these outcomes. Any revised model requires its own protocol and validation. Improving our validation still matters: the full and canonical-small target subsets are now compared under the same bulk-profile calculation and actual treatment populations. The project still lacks convincing study-independent evidence for a large leaderboard gain; this experiment does not supply it.

## Verification and reproducibility

Seven numerical tests pass, including recovery of a known synthetic feedback network, indirect chain propagation, exact unchanged-context identity and missing-data abstention. Recipient-evidence poisoning passes for all three studies (H1, Orion and CD4), covering 25, 573 and 604 source-summary rows respectively. All 300 panel output genes remain exactly V63, with own-target knockdown applied separately only if cells were generated.

The result auditor checks bitwise equality of saved V63/V66 baseline arrays, query/bulk profiles and measured masks; all 3,275 arm/target/context records share their treatment populations; the 655-pair cohort size and all macro error ratios independently reconcile. It checks saved arithmetic and constraints, not a full independent reimplementation of every raw extraction and graph prediction. No private treatments or D3 reserve are used. Original releases and their seals remain unchanged.

Reproduce with the commands in [V67_STRATEGY.md](V67_STRATEGY.md), then run `audit_v67_results.py` and `document_v67_results.py`. Scientific source and compact receipts are committed; large arrays are regenerable and ignored. The pre-outcome protocol was committed as `d75504c`. A separate result manifest binds the protocol, audit, outputs and their inputs; README is mutable and excluded. Nothing is uploaded or pushed automatically.
'''
    (m.ROOT/'V67_RESULTS.md').write_text(text,encoding='utf-8')
    write(m.A/'network_diagnostics.json',dict(contexts=traces,scope='Correction magnitude and source support, not independent prediction accuracy.'))
    path=m.ROOT/'README.md';old=path.read_text(encoding='utf-8-sig');begin='<!-- v67-network:begin -->';end='<!-- v67-network:end -->'
    block=f'''{begin}
## V67: new intervention-network hypothesis tested; no release

Tested feedback rewiring from measured cross-target CRISPR responses on **655 target/context pairs** with whole-study exclusion, fixed target partitions, V63/V66 comparators and two matched controls. Intended profile error changes **{s['profile_MSE_ratio']-1:+.4%} vs V63**, discrimination **{s['pds_change']:+.7f}**. **Release criteria failed; no V67 submission generated.** Seven numerical tests, all three study-poisoning audits and the independent saved-result audit pass. This is reused public biology, not blind validation. V66 remains best published at 0.150882188.

[Results and failures](V67_RESULTS.md) · [Frozen strategy and commands](V67_STRATEGY.md). Prior status blocks below are historical; original release artifacts remain unchanged.
{end}

'''
    if begin in old:
        start=old.index(begin);finish=old.index(end,start)+len(end);old=old[:start]+old[finish:].lstrip('\n')
    path.write_text(block+old,encoding='utf-8')
    print('Documented failed V67 network primary; no submission generated.',flush=True)


if __name__=='__main__':main()
