"""Generate the upstream review report and scientific figures from validated evidence."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from build_submission import sha256,write_json
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'analysis/upstream_review'


def table(frame):
    columns=list(frame.columns)
    def cell(value):
        return str(value).replace('|','\\|').replace('\n',' ')
    return '\n'.join(['| '+' | '.join(columns)+' |','| '+' | '.join(['---']*len(columns))+' |']+
                     ['| '+' | '.join(cell(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)])


def main():
    validation=json.loads((OUT/'validation.json').read_text())
    assert validation['passed']
    for name,digest in validation['artifact_sha256'].items():
        assert sha256(OUT/name)==digest,name
    for name,digest in validation['root_code_sha256'].items():
        assert sha256(ROOT/name)==digest,name
    review=json.loads((OUT/'repository_review.json').read_text())
    coverage=json.loads((OUT/'h1_coverage_probe.json').read_text())
    centering=json.loads((OUT/'centering_validation.json').read_text())
    promoter=json.loads((OUT/'promoter_probe.json').read_text())
    identity=json.loads((OUT/'historical_input_identity.json').read_text())
    counts=pd.read_csv(OUT/'h1_public_target_coverage.csv')
    stats=pd.read_csv(OUT/'control_moment_diagnostics.csv')
    candidates=pd.read_csv(OUT/'centering_candidates.csv')
    candidates['gain']=candidates.setting.map(lambda s:json.loads(s)['gain'])
    candidates['centering']=candidates.setting.map(lambda s:json.loads(s)['centering'])
    gain_one=candidates[candidates.gain==1].sort_values('centering')
    split=pd.DataFrame([{'Public H1 split':s['split'],'2026 targets':s['overlapping_2026_targets'],
                         'Available cells':s['cells'],'Control cells':s['control_cells']} for s in coverage['sources']])
    extra=[]
    for row in counts[counts.new_beyond_training].itertuples():
        source=next(s for s in ['validation','test'] if getattr(row,s)>0)
        extra.append({'Target':row.target_gene,'Public 2025 split':source,'Available target cells':getattr(row,source)})
    moment_table=[]
    for r in stats.itertuples():
        method='Observed sample' if r.method=='original_control_sample' else f'Dual moments, pool {r.pool_k}'
        moment_table.append({'Context':r.context,'Generator':method,'Pooled-profile L1':f'{r.pooled_probability_L1_from_training:.6f}',
                             'Held-out control log-bulk MSE':f'{r.held_out_control_log_bulk_MSE:.8f}',
                             'Variance / held-out controls':f'{r.CPM_variance_ratio_to_heldout_controls:.3f}'})
    colors=['#6b7280','#167d8d','#b47726']
    fig,axes=plt.subplots(2,2,figsize=(12,8.6),layout='constrained')
    ax=axes[0,0]
    ax.bar(split['Public H1 split'],split['2026 targets'],color=['#167d8d','#b47726','#815ca1'])
    for i,n in enumerate(split['2026 targets']):ax.text(i,n+.4,str(n),ha='center')
    ax.set(ylim=(0,16),ylabel='Targets overlapping the 2026 panel',title='Public H1 metadata: 13 existing + 12 additional targets')
    ax=axes[0,1]
    ax.plot(gain_one.centering,gain_one.relative_zero_loss,marker='o',color='#167d8d')
    ax.axhline(1,color='#6b7280',linestyle='--',linewidth=1)
    ax.set(xlabel='Across-target subtraction strength (gain = 1)',ylabel='Donor prediction MSE / zero-effect MSE',
           title='CD4: centering did not pass nested selection',xticks=[0,.5,1])
    x=np.arange(3);width=.24
    groups=[stats[stats.method=='original_control_sample'],stats[(stats.method!='original_control_sample')&(stats.pool_k==1)],
            stats[(stats.method!='original_control_sample')&(stats.pool_k==4)]]
    labels=['Observed 400-cell sample','Dual moments, one-cell template','Dual moments, four-cell pooling']
    for index,(group,label,color) in enumerate(zip(groups,labels,colors)):
        group=group.set_index('context').loc[list('ABC')]
        axes[1,0].bar(x+(index-1)*width,group.held_out_control_log_bulk_MSE,width,label=label,color=color)
        axes[1,1].bar(x+(index-1)*width,group.CPM_variance_ratio_to_heldout_controls,width,label=label,color=color)
    axes[1,0].set(xticks=x,xticklabels=list('ABC'),yscale='log',ylim=(5e-5,.002),ylabel='Log-bulk MSE on disjoint controls',
                  title='Neutral control reconstruction: smaller profile drift')
    axes[1,1].axhline(1,color='black',linestyle='--',linewidth=1)
    axes[1,1].set(xticks=x,xticklabels=list('ABC'),ylabel='Total CPM variance / held-out controls',
                  title='One-cell templates retain more dispersion')
    handles,legend_labels=axes[1,0].get_legend_handles_labels()
    fig.legend(handles,legend_labels,loc='outside lower center',ncol=3,fontsize=9)
    fig.suptitle('Upstream review: local validation, not a leaderboard prediction',fontsize=14)
    fig.savefig(OUT/'upstream_review.png',dpi=180)
    fig.savefig(OUT/'upstream_review.pdf')
    plt.close(fig)
    repos=review['repositories']
    repo_table=pd.DataFrame([{'Repository':r['repository'],'Pinned HEAD':r['head'],'Main commits':r['commits_main'],
                            'All reachable commits':r['commits_all_refs']} for r in repos])
    center_table=pd.DataFrame([{'Centering':f'{r.centering:g}','Relative MSE':f'{r.relative_zero_loss:.6f}',
                                'Selection proxy':f'{r.selection_score:.6f}','Fidelity proxy':f'{r.direction_fidelity_proxy:.6f}',
                                'Jaccard proxy':f'{r.jaccard_proxy:.6f}'} for r in gain_one.itertuples()])
    report=f'''# Review of AtlasShift and Arc Virtual Cell Atlas

Review date: 2026-10-05. This report combines the requests to examine both user-cloned repositories and learn from their commits. Code, historical methods, configuration and tutorial cells were read before running the small numerical experiments below. Neither repository was modified. The five original competition input hashes match AtlasShift's historical input manifest; that does not establish identical hidden evaluation references or anchors.

## Conclusion and next submission strategy

The strongest additions to test are **12 more directly observed public H1 targets** and **separate calibration of mean CPM and pooled-count profiles while retaining one-cell templates**. Promoter-neighbor effects warrant an independently calibrated prior. The held-target donor experiment did not support adding across-target subtraction to our CD4 estimator. None of these local diagnostics proves a higher hidden leaderboard score, and this review does not create or upload a new submission.

V4 remains the current published candidate: entry `w5criCkgkCWj1QOgP08E`, overall **-0.05067391016092112**. Its raw expression error rose from 1.014294 for V2 to 1.247733 even though the scaled expression component stayed at zero. Better profile calibration addresses a measured weakness alongside response coverage. Existing V4 source code and effect arrays are unchanged; the effect-cache hash was checked afresh, and package size was checked against its readiness record. The 3.1 GiB package was not rehashed in this review.

## Commit history and implementation lessons

{table(repo_table)}

The origin main hashes were also checked with read-only `git ls-remote`. Complete reachable commit metadata, paths and executable-Python classifications are in `analysis/upstream_review/*_commits.csv`. AST comparison ignores docstrings; notebook code and documentation were reviewed separately, so a count of Python-file changes does not describe every notebook update.

### AtlasShift

Source: [Virtual-Cell-Challenge-2026](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026).

Only two of its 14 commits change executable Python. Initial commit `8257ada256112e1e7374241fcc1f744a575089ad` contains the methods, reproducibility guide, configuration, tests and input hashes. Commit `490dcaac2c824b92a7997b7e8d2970b1690c3570` reorganizes the implementation into five scripts and removes much of that documentation. The remaining changes adjust documentation and names/docstrings. The historical methods were recovered using `git show`; this history does not provide 38 successive trained model experiments.

The author reports historical v38 score **0.1545618019323822**, rank 82 of 839 on September 10. These are author-reported historical results, not a current rank or an independently reproduced score. A documented v39 result is slightly lower. Matching original input hashes makes the method relevant, but unspecified hidden anchor versions prevent a strict score comparison with our current V4.

Useful mechanisms in the reviewed implementation:

- Source coverage spans K562, HCT116, HEK293T, the complete released H1 target panel, and CD4 conditions. Default source-family weights are 2, 1, 1, 2 and 0.5 respectively. CD4 conditions form one averaged family. Its full K562 source has broader gene coverage than our 7,681-gene pseudobulk overlap. Our small HCT116/HEK293T pilot batches do not provide comparable response training support.
- The full source-target panel is built before source centering and output selection. Across-target source means exclude the target's own-gene coordinate. Subtracting a shared effect can prevent repeatedly transferring generic responses; it can also remove real biology shared by targets. Our empirical test below evaluates that tradeoff for CD4 donor transfer.
- Mean cell-normalized CPM and library-weighted pooled counts are modeled separately. The published defaults use response amplitudes 0.6 and 0.3, with clipping. The count generator iteratively reconciles both requested profiles, projects incompatible moments, rounds to integer counts and preserves each cell's requested library size. Defaults are hypotheses to calibrate for our estimator, not established optimum values for A/B/C.
- Its default template averages four controls with nearby library sizes. Our paired experiment shows that this can greatly reduce cell variance. A one-cell template is a more promising adaptation for the next experiment.
- GENCODE gene-boundary TSS pairs within 5 kb produce neighboring-gene suppression caps, with remaining expression 0.15 at the shortest distances and a distance-dependent relaxation. These are heuristic caps; renormalizing profiles after applying them means they are not strict final expression ceilings.
- `compact.py` changes sparse layout/cell order while checking count equality. This could reduce upload size without improving numerical predictions. `pack.py` uses a private adapter requiring VCC 0.2.0, while our working environment has 0.2.2. Retain our official preparation path.

We imported the pinned, reviewed numerical `model.py` for focused tests; no upstream source was copied. Its MIT license is retained in the original clone. The large preparation pipeline was not run: full-source downloads/scans are substantial, and the existing environment was not upgraded to its dependency versions.

### Arc Virtual Cell Atlas

Source: [arc-virtual-cell-atlas](https://github.com/ArcInstitute/arc-virtual-cell-atlas).

The checked main branch has 59 commits; all reachable refs contain 62. Its current files provide data descriptions, access examples and five Python/R notebooks, with no trained perturbation predictor. The two historical executable-Python changes concern a subsequently removed TileDB-SOMA access helper. The useful recent changes concern release pinning, metadata, access methods and tutorials.

The challenge data here are the released **2025 H1** dataset: 150 training, 50 validation and 100 test target labels. The old three-metric description is not the six-component 2026 scoring contract. Tahoe-100M is a drug-treatment resource, and scBaseCount is observational expression data: neither supplies interchangeable CRISPRi response labels. They can support lineage/control characterization if sample provenance and treatment status are checked. Study-level disease annotations can include healthy samples. Unique-UMI `X`, EM counts and Velocyto layers require distinct treatment rather than being combined indiscriminately.

Remote branch `origin/docs/fix-egress-guidance`, commit `8fc7a9e368f0fa7098858019ba41a492c8062c69`, corrects older access/billing guidance, but was **not merged into the checked main**. Our bounded H1 HTTPS metadata reads worked anonymously without a billing project. This verifies those URLs only; it is not a claim about all atlas access or billing arrangements. We did not execute example whole-species reads or download complete atlases.

## 1. Public H1 response coverage: 13 to 25 potential targets

Only observation labels, cell identifiers and gene axes were read through bounded HTTP ranges: the H1 expression matrix was not accessed. Each split has 18,080 gene names, with 18,077 shared with the 2026 axis. Coverage counts indicate available cells before response QC.

{table(split)}

The three public splits offer **25 distinct 2026 targets**, including these **12 beyond the training-only source used in V4**:

{table(pd.DataFrame(extra))}

All 38,176 control-cell identifiers in training also appear in validation and test. Count equality was not checked. Future source validation must partition controls globally by identifier and avoid treating the public splits as independent control replicates.

These are publicly released 2025 labels, not hidden 2026 evaluation data. Before estimating effects, extract the relevant sparse expression rows, map measured symbols, check guide/batch metadata, verify on-target knockdown and minimum usable cell counts, and use matched source controls. A shared gene axis does not justify filling unmeasured genes with zero. URLs, object generations, sizes, ETags and observed coverage are recorded in `h1_coverage_probe.json` and the 900-row/300-row coverage tables.

## 2. Across-target centering: not supported for our CD4 donor estimator

The existing disjoint CD4 donor contrasts supply **{centering['targets']} targets, {centering['contrasts']} target-condition contrasts and {centering['genes']:,} genes**. We tested gain {{0.5, 0.75, 1}} and centering {{0, 0.5, 1}}. For each condition, the source center is the unweighted average of valid other training-target effects. Target conditions stay together; the target's own-gene coordinate and invalid values are excluded. The response remains softly shrunk by `0.04 / (0.04 + SE²)` and clipped to +/-3.

Five outer target folds and three inner folds prevent a held target from contributing to its centering fit or setting selection. Selection maximizes the existing direction/overlap/error proxy among candidates whose inner prediction error does not exceed the fixed V4 baseline. Adoption additionally requires no outer loss increase and a better composite proxy. All five folds select **gain 1, centering 0**.

At gain 1, the independently aggregated outer diagnostics are:

{table(center_table)}

These are donor/normal-approximation 400-cell proxies, not published FID/Jaccard scores or the official hidden-reference expression metric. Full subtraction slightly increases direction yield but worsens overlap and squared error. Keep the present CD4 centering policy. This does not rule out centering for a different cross-lineage estimator or a broader source panel.

## 3. Paired count-generation experiment: separate moments, preserve dispersion

For each context, split the 18,400 original controls into **9,200 training and 9,200 distinct held-out controls**, with a recorded seed. Estimate both the equal-cell mean CPM probability and the library-weighted pooled probability using training controls. Generate 400 neutral cells. The baseline and one-cell dual-moment reconstruction use the **same guide/depth-stratified 400 donors**; the four-cell version averages 1,600 depth-sorted controls.

All generated groups contain nonnegative integer counts and preserve their requested cell depths. The control log-bulk diagnostic is the mean squared difference between `log1p(50000 * pooled_probability)` profiles of generated and held-out controls; it is not the leaderboard metric. Variance is the sum of per-gene cell-normalized variances relative to held-out controls.

{table(pd.DataFrame(moment_table))}

One-cell dual moments reduce pooled-profile L1 drift from about 0.03 to below 0.0005, with variance ratios approximately 0.97-1.10. Four-cell pooling retains only approximately 0.30-0.38 of held-out variance and increases gene detection. This is a neutral control reconstruction test: there are no perturbation responses or own-target knockdown, so response performance and DE-call calibration remain untested. Both CPM-moment residuals and detection/depth diagnostics are retained in the CSV.

![Validated source coverage, centering and count-generation diagnostics](analysis/upstream_review/upstream_review.png)

## 4. Promoter neighbors: relevant evidence, no universal knockdown rule

The 58,953,903-byte GENCODE v47 annotation is size/SHA-256 checked. Independently parsing its unambiguous gene-boundary TSS positions reproduces **{promoter['pairs']} measured neighbor pairs across {promoter['targets_with_neighbors']} targets** within 5 kb. The reconstructed upstream CSV also matches its historical SHA-256 exactly.

Cached CD4 data yield **{promoter['cd4_observations']} qualifying target-condition observations**, with **{promoter['significant_negative']} significant negative and {promoter['significant_positive']} significant positive neighboring-gene responses** at q < 0.05. Repeated conditions are not independent genes. The descriptive filter requires at least 100 target cells, significant on-target knockdown, multiple-guide estimates, and no distal-off-target/low-target-expression flags; it deliberately retains the publisher's neighboring-gene-KD flag so collateral behavior can be examined. This does not alter the stricter V4 training filter.

Gene-boundary TSS is not the experimentally targeted guide position or the context's active promoter. The evidence supports testing a distance-, expression- and source-confidence-dependent neighboring-gene prior on held-out targets; it does not justify forcing every neighbor to 15% expression. Complete target/neighbor/condition/distance/q-value records are saved.

## Proposed sequence for the next candidate

1. Extract the 12 added public H1 targets and their matched controls with a recorded download budget. Check source counts and knockdown, then extend V4's direct-H1 coverage only for usable targets. Keep global control-ID splits.
2. Estimate mean-CPM and pooled-profile perturbation effects separately. Tune amplitudes on held-out source targets. Test one-cell dual-moment generation against the existing generator with identical effect arrays and donor IDs; assess expression error, response direction/overlap, detection and per-gene variance. Expand beyond neutral tests before generating a complete package.
3. Test promoter-neighbor corrections separately under target-held-out validation. Keep CD4 across-target subtraction at zero unless a different source experiment supports it. Obtain substantial HCT116/HEK293T perturbation data before assigning them response weights.
4. Generate and validate a candidate through the existing VCC 0.2.2 official preparation path only after these ablations support it. Record all changes and input hashes. A future hidden score remains unknown until a submitted candidate is scored.

## Reproduction and evidence

Run from the project root in **Command Prompt**. Existing V4 inputs/caches must be present: the three original H5ADs, `gene_names.csv`, `pert_counts.csv`, donor inputs referenced by `fit_lineage_submission.donor_inputs`, `reference/cd4_probe/merged.csv` / `.npz`, and the V4 reference model/package. The original clones must be present at the pinned commits above; the moment test imports `Virtual-Cell-Challenge-2026/model.py` after checking its hash. No additional dependency installation was needed in the existing `.competition_env`.

```cmd
.\\.competition_env\\Scripts\\python.exe review_upstream_repositories.py
.\\.competition_env\\Scripts\\python.exe probe_atlas_sources.py
.\\.competition_env\\Scripts\\python.exe evaluate_atlas_transfer.py
.\\.competition_env\\Scripts\\python.exe -m unittest test_upstream_review -v > analysis\\upstream_review\\test_results.txt 2>&1
.\\.competition_env\\Scripts\\python.exe validate_upstream_review.py
.\\.competition_env\\Scripts\\python.exe document_upstream_review.py
```

The probe budgets 128 MiB of H1 ranges per split, reads metadata only and downloads the pinned approximately 59 MB annotation if absent. HTTP caches reside in ignored `reference/upstream_review/`, with the existing H1 training cache reused. Cached reruns reduce network use. `--metadata-only` skips annotation/neighbor analysis; `--centering-only` and `--moments-only` isolate the numerical experiments. Full validation/documentation requires all outputs.

Four focused tests passed: held-target centering exclusion, missing/own-gene handling, deterministic integer moment matching with exact depths, and infeasible moment rejection. The independent validator checks clone cleanliness/pins, reviewed file and root-script hashes, target-fold disjointness, all 900 coverage records, source identities, original GTF reconstruction, historical input hashes, all nine moment groups and V4 model preservation. It checks stored numerical summaries and focused replay, not hidden scores or every source-expression value. `validation.json` records its scope and hashes; `documentation_manifest.json` records this report/figure/code provenance.

The root scripts and compact evidence under `analysis/upstream_review/` are version-controlled. User clones, HTTP caches and original large inputs remain local research inputs. Do not push full atlases or prediction packages into Git.
'''
    report_path=ROOT/'UPSTREAM_REPOSITORY_REVIEW.md';report_path.write_text(report,encoding='utf-8')
    section='''## Upstream repository review and validated improvement experiments

Both cloned repositories were reviewed through their reachable commit history, historical methods and current code/tutorials. [UPSTREAM_REPOSITORY_REVIEW.md](UPSTREAM_REPOSITORY_REVIEW.md) records pinned commits, source differences, tests, scientific figures and exact reproduction commands. AtlasShift has 14 commits but only two executable-Python-changing commits; Arc Atlas has 59 main commits (62 across refs) and supplies data/access tutorials rather than a trained predictor.

Public 2025 H1 metadata extends potential direct-response coverage from 13 to **25** of our 300 targets, with 12 additional targets. All 38,176 control identifiers recur across the released splits, so validation must partition controls by identifier. Response counts/knockdown still require extraction and QC. Paired neutral-control tests show that separate mean-CPM/pooled-profile count matching reduces profile drift while one-cell templates retain variance near the controls; four-cell pooling reduces variance to roughly 30-38%. Five held-target CD4 folds select zero across-target subtraction. Independently reconstructed promoter pairs cover 74 targets, with descriptive neighboring-gene suppression evidence that requires held-out calibration.

Root scripts: `review_upstream_repositories.py`, `probe_atlas_sources.py`, `evaluate_atlas_transfer.py`, `test_upstream_review.py`, `validate_upstream_review.py`, `document_upstream_review.py`. Four focused tests and the independent validator passed; evidence is under `analysis/upstream_review/`. The next candidate should first evaluate additional H1 responses and one-cell dual-moment calibration. This review preserves V4 and does not generate or upload a new prediction package.

'''
    readme_path=ROOT/'README.md';readme=readme_path.read_text(encoding='utf-8-sig')
    heading='## Upstream repository review and validated improvement experiments\n'
    anchor='## V4 lineage-guided response coverage candidate\n'
    if heading in readme:
        start=readme.index(heading);end=readme.index(anchor,start);readme=readme[:start]+readme[end:]
    assert anchor in readme
    readme=readme.replace(anchor,section+anchor,1)
    readme_path.write_text(readme,encoding='utf-8')
    write_json(OUT/'documentation_manifest.json',{'script_sha256':sha256(Path(__file__)),
        'validated_evidence_sha256':sha256(OUT/'validation.json'),
        'files_sha256':{str(path.relative_to(ROOT)).replace('\\','/'):sha256(path) for path in [report_path,OUT/'upstream_review.png',OUT/'upstream_review.pdf']},
        'scope':'Report and scientific figures generated from the hash-checked validated evidence; README section inserted before the V4 regeneration boundary.'})
    print('Wrote report, README section, PNG/PDF and documentation manifest.',flush=True)


if __name__=='__main__':main()
