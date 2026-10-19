"""Write V69 results after scientific and complete-payload audits finish."""
from pathlib import Path
import json
import numpy as np,pandas as pd
from v57_pie import read,write,sha
import v69_count_bridge as m


def table(headers,rows):
    def fmt(x):
        if x is None:return 'unavailable'
        if isinstance(x,float):return f'{x:.6f}'
        return str(x)
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(fmt(x) for x in row)+' |' for row in rows])


def main():
    ready=read(m.A/'submission_ready.json');audit=read(m.A/'independent_validation.json');assert audit['passed'] and audit['package_sha256']==ready['sha256']
    dev=read(m.A/'selection.json');conf=read(m.A/'confirmation_assessment.json');matched=read(m.A/'matched_control_assessment.json');diag=read(m.A/'diagnostic_assessment.json');choice=dev['selected'];b=conf['summary']['v66'];v=conf['summary'][choice]
    ratio=float(np.mean(list(conf['ratios'][choice].values())));failed=[k for k,z in conf['checks'].items() if not z]
    metrics=['pds','fidelity','reach','jaccard','nmae','profile_square_error']
    source=[]
    for name in m.SOURCES:
        r=read(m.A/'source_shapes'/(name+'.json'));source.append([name,len(r['targets']),sum(v['supported_coordinates'] for v in r['targets'])])
    nulls=[]
    for name in ('H1','HCT116','HEK293T','CD4_Rest_D1','CD4_Stim48_D2'):
        vals=[read(m.A/'representative_null'/f'{name}__{a}.json')['calls'] for a in ('template','v66','bridge','hybrid25','near_zero_bridge')]
        prefix=read(m.A/'null'/f'{name}__bridge.json')['calls'];nulls.append([name,*vals,prefix])
    corrected=[]
    for stage in ('development','confirmation'):
        arms=['v66',*m.ARMS] if stage=='development' else ['v66',choice,'shuffled_shape','shuffled_hybrid25']
        for a in arms:
            rows=[read(p) for p in (m.A/stage).glob('*__panel.json') if read(p)['arm']==a]
            studies=[]
            for group in ('H1','Orion','CD4'):
                vv=[r['common_expression_numerator']/r['cases'] for r in rows if r['study']==group];studies.append(float(np.mean(vv)))
            corrected.append([stage,a,*studies,float(np.mean(studies))])
    groups=pd.read_csv(m.A/'prediction_group_statistics.csv');provenance=pd.read_csv(m.A/'prediction_control_rows.csv');drift=provenance.predicted_total_counts-provenance.original_total_counts
    raw=read(m.A/'results_audit.json')
    text=f'''# V69 count-bridge release

Upload artifact: `experiments/v69_count_bridge/prediction.vcc`. Candidate: **{choice}**. Official preparation and complete-payload validation passed. Public confirmation support: **{conf['release_supported']}**. Overall release support including null checks: **{ready['fixed_public_development_and_confirmation_guards_passed']}**. This file is **{'experimental' if ready['experimental'] else 'public-supported, not privately confirmed'}**; no improved hosted score is claimed. V66 remains the highest published project score, 0.150882188.

## Change and decision

The V66 response means are preserved as model inputs. An independently implemented Gamma-Poisson count bridge changes raw integer counts through posterior intensities, binomial removal, Poisson addition and beta/gamma dispersion coupling. It is a working probabilistic model, not a claim that all raw cells are NB distributed. Source target/gene dispersion changes use split-view sign support, minimum molecule/cell counts and study-balanced shrinkage. All panel genes abstain from dispersion fitting. Shuffled evidence is an inference ablation, not a separately trained shuffled model. See [frozen protocol](V69_STRATEGY.md).

For hybrid25, exactly100 of400 cell positions receive bridge counts;300 retain the exact parent V66 rows. The mixture does not inherit the pure bridge's null identity. The mean_bridge/shape_bridge candidates use the full bridge population. No exact original library-size repair follows the bridge; final count validity and normalization are measured independently. Production samples only the selected100 bridge rows, using the same row-independent probability law with a documented new random-stream realization. Earlier code/metadata and an unused partial output are preserved; fitted arrays and selection rules are unchanged. Fresh preflight and deterministic replay verify the deployed sampler.

Development-eligible arms: {dev['eligible']}. Predeclared fallback used: **{dev['experimental_fallback']}**. Selection is locked before the remaining target cases are scored. Failed confirmation checks: {', '.join(failed) or 'none'}. Matched-shuffle profile comparison: **{matched['beats_matched_control_profile']}** against **{matched['control']}**. Unresolved metrics: {conf.get('unresolved_metrics',[])}. Confirmation cannot change the selected candidate. A separately frozen pre-development correction adds a shuffled quarter-mixture control, because the original full shuffled bridge was not a matched control for hybrid25. It changes no candidate or parameter. The confirmation assessor also preserves undefined DE measurements as unavailable and failed guards rather than attempting numeric comparisons with null.

Candidate paired equal-study profile error / V66: **{ratio:.6f}** ({(ratio-1)*100:+.2f}%). These are reused-public raw results, not the leaderboard's anchored accuracies. Raw error and noise-corrected expression are both shown, because realistic sampling changes prediction variance.

## Confirmation

{table(['Diagnostic','Fold-specific V66',choice],[[k,b[k],v[k]] for k in metrics])}

Higher PDS/FID/REACH/JAC is favorable; NMAE and raw squared profile error are lower-better errors. Missing metrics remain unavailable. Study profile ratios: {conf['ratios'][choice]}. Different cases/target panel sizes prevent direct comparison to historical raw PDS numbers or hosted scores. Target confirmation reuses the existing public cohorts and is not independent-study validation.

## Mean versus decoder diagnostic

{table(['Arm','Raw profile error','FID','REACH','JAC'],[[a,*[s[k] for k in ('profile_square_error','fidelity','reach','jaccard')]] for a,s in diag['summary'].items()])}

Diagnostic representation checks: {diag['checks']}. Ten target/context pairs use distinct treatment fitting/scoring halves (100 or200 cells per synthetic prediction). Oracle means and shapes use the fitting half only; the scoring half is separate. The measured fitting-half mean is noisy and not a perfect population oracle. The bridge diagnostic additionally clips bulk ratios to the production +/-3 log2 range and cannot activate a gene with zero fitted control probability; the old-oracle input is not constrained identically. These limitations can contribute to reconstruction error and prevent attributing every difference to distribution family. Old-oracle generation receives both measured CPM/bulk means; bridge-oracle generation uses fitted bulk means and marginal shape. Realized means differ, so this is not a strictly orthogonal shape-only comparison. The real_fit arm is the separate treatment fitting half itself, not duplicated scoring cells. These labelled-data arms cannot be used for A/B/C inference or production selection. The source-trained deployable candidates use no query treatment labels. This diagnostic fails the raw-profile representation guard; that failure is retained. It does not justify claiming a generally superior decoder.

## Count-null checks and sampling correction

{table(['Context','Representative template','V66','Full bridge','Hybrid25','Near-zero bridge','Original prefix bridge'],nulls)}

The original frozen null used the first400 query-fit controls. It can overrepresent source batches; exact count identity consequently preserved578/1041 DE calls in HCT116/HEK293T. That observation prompted a separately frozen [representative-null addendum](experiments/v69_count_bridge/analysis/representative_null_plan.json): fixed depth-stratified sampling of400 distinct query-fit controls, with no seed search. Original outcomes are retained. This addendum changes no fitting or candidate-selection rule. All arms in the additional comparison use the same selected rows. Zero calls in these comparisons do not establish a universal zero false-positive rate. The hybrid retains old-generator artifacts even when the full bridge passes.

## Noise-corrected expression

{table(['Stage','Arm','H1 signed numerator/case','Orion','CD4','Equal study'],corrected)}

These are signed, jackknife-corrected expression numerators averaged per case within context, then equally over contexts/studies. Lower is favorable; a negative estimate may arise from noise correction. Denominators/rho remain in the original panel receipts, and unresolved denominators are not floored into a claimed overall score. This descriptive analysis does not override the frozen raw-profile selection rule.

## Source evidence and limitations

{table(['Raw source','Eligible target conditions','Supported target/gene coordinates'],source)}

HCT116/HEK293T share the Orion study; both CD4 conditions share the CD4 study. Entire recipient studies are excluded from both mean and dispersion evidence in public tests. A source view split uses alternating cells, not independent guides or biological replicates. The dispersion moments are not a Bayesian posterior or NB maximum-likelihood fit. Depth conditioning, clonal/batch effects and sparse counts can bias them. The new coupling preserves empirical control dependence only approximately after shape changes. No new biological dataset or untouched D3 reserve was opened. More seeds or reused public cases do not solve the validation-transfer gap observed with V66.

## Artifact and verification

- Rows/genes/groups: **360,000 /18,533 /900**,400 cells per target/context.
- Payload size: **{ready['bytes']:,} bytes**. Nonzero entries: **{ready['package_meta']['nnz']:,}**.
- SHA-256: `{ready['sha256']}`.
- Independent payload audit: **{audit['verified_rows']:,} rows**, all stored counts/indices/pointers, all control identities, **{audit['exact_replayed_groups']}** deterministic regenerated groups; **{audit['unchanged_parent_rows_verified']:,}** unchanged parent rows checked.
- Public result audit: **{raw['count_populations']}** count populations and **{raw['manual_rank_coordinates']}** manually reproduced tied-rank coordinates, independent BH/LFC/profile checks.
- Maximum final library: **{int(groups.maximum_library.max()):,}**. Cells with changed total count: **{int((drift!=0).sum()):,}**. Library changes are allowed; exact original totals are not claimed.
- Six probability/count numerical tests and the quarter-population/cache equivalence checks pass. Package validity is not evidence of hidden accuracy.

The separate official environment runs unchanged `vcc prep`, following lossless float32/int64 staging. Every value is compared during staging, and the decompressed official payload is compared to the generated H5AD. Source arrays, model arrays and count populations are ignored regenerable outputs; compact receipts and root scripts are committed locally. No automatic upload or Git push occurs.

## Submit from Windows Command Prompt

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v69_count_bridge\\prediction.vcc -m "V66 probability-coupled count mixture v69 experimental" --wait
```

Keep the resulting entry ID. If waiting is interrupted after upload, use `vcc status ENTRY_ID --wait` rather than creating a duplicate submission. This release has no hosted result yet.
'''
    (m.ROOT/'V69_RESULTS.md').write_text(text,encoding='utf-8')
    p=m.ROOT/'README.md';old=p.read_text(encoding='utf-8-sig');assert '<!-- v69-release:begin -->' not in old
    block=f'''<!-- v69-release:begin -->
## V69 ready: probabilistic count bridge, experimental

Created and independently validated a new 360,000-cell upload artifact using **{choice}**, frozen V66 response means and source-learned count dispersion. Confirmation profile error / V66 **{ratio:.6f}**; strict public support **{ready['fixed_public_development_and_confirmation_guards_passed']}**. The pure bridge preserves zero response, while the population mixture retains baseline null artifacts. Both original prefix-template and corrected representative-control nulls are retained. **No improved hosted score is claimed; V66 remains highest published at 0.150882188.**

[Results, limitations and upload command](V69_RESULTS.md) · [Frozen protocol](V69_STRATEGY.md) · [Independent audit](audit_v69_release.py). Older blocks below are historical.
<!-- v69-release:end -->

'''
    p.write_text(block+old,encoding='utf-8')
    print('V69 results and README written',flush=True)

if __name__=='__main__':main()
