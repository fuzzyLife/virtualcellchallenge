"""Render V37 results directly from saved evidence; do not promote failed blends."""
from run_v37_ensemble import ROOT,A,ARMS,read

NAMES={'pds':'PDS','fidelity':'FID','reach':'REACH','jaccard':'JAC','nmae':'NMAE error','common_expression_ratio':'Expression error'}

def main():
 selection=read(A/'selection.json');audit=read(A/'independent_audit.json');assert audit['passed'];means=read(A/'mean_results.json')['summary'];confirmed=read(A/'confirmation_decision.json') if (A/'confirmation_decision.json').exists() else None
 release=read(A/'ensemble_release.json');selected=release['selected_arm']
 ext=read(A/'recalibrated/selection.json');extaudit=read(A/'recalibrated/independent_audit.json');assert extaudit['passed']
 if release['namespace']!='primary':
  directory=A/release['namespace'];confirmed=read(directory/'confirmation_decision.json') if (directory/'confirmation_decision.json').exists() else None
 noise=read(A/'noise_recalibrated/selection.json') if (A/'noise_recalibrated/selection.json').exists() else read(A/'noise_recalibrated/null_gate.json')
 noiseaudit=read(A/'noise_recalibrated/independent_audit.json');assert noiseaudit['passed']
 promoted=confirmed is not None and confirmed['selected']==selected and selected!='v24';preferred=selected if promoted else 'v24'
 means.update(read(A/'recalibrated/mean_results.json')['summary']);means.update(read(A/'noise_recalibrated/mean_results.json')['summary'])
 rows=[]
 for arm,s in {**selection['summaries'],**{k:v for k,v in ext['summaries'].items() if k!='v24'},**{k:v for k,v in noise.get('summaries',{}).items() if k!='v24'}}.items():
  rows.append('| '+arm+' | '+' | '.join('unresolved' if s[k] is None else f'{s[k]:.9f}' for k in NAMES)+' |')
 checks=[]
 for arm,c in {**selection['checks'],**ext['checks'],**noise['checks']}.items():checks.append('| '+arm+' | '+(', '.join(k for k,v in c.items() if not v) or 'none')+' |')
 text='# V37 component ensemble results\n\n'
 text+='**Preferred model: '+preferred.upper()+'.** '
 if selected=='v24':text+='None of the five primary blends or either corrective calibrated candidate passed the applicable fixed guards. No new upload package was created; the original V24 package remains preferred. Additional count confirmation was not opened because there was no qualifying candidate.\n\n'
 elif not promoted:text+='Development selected '+selected+', but confirmation did not establish a passing replacement. V24 remains preferred.\n\n'
 else:text+='Development and additional count confirmation passed for '+preferred+'. This is public diagnostic evidence, not a proven leaderboard gain.\n\n'
 text+='## Development: unscaled public count metrics\n\nPDS/FID/REACH/JAC increase when favorable; NMAE and expression error decrease. Equal contexts within study, then equal studies. Undefined context NMAE is omitted, never zero; entirely unresolved studies fail.\n\n'
 text+='| Arm | '+' | '.join(NAMES.values())+' |\n|---|'+'---:|'*len(NAMES)+'\n'+'\n'.join(rows)+'\n\n'
 text+='## Frozen guard failures\n\n| Arm | Failed checks |\n|---|---|\n'+'\n'.join(checks)+'\n\n'
 nulls=read(A/'null_results.json')['records'];rn=read(A/'recalibrated/null_results.json')['records'];nn=read(A/'noise_recalibrated/null_results.json')['records']
 text+='## Corrected control-null DE calls\n\n| Context | Actual independent controls | V24 | Noise mix | Recalibrated row mixture | Recalibrated noise mix |\n|---|---:|---:|---:|---:|---:|\n'
 for family in sorted({r['family'] for r in nulls}):
  q={r['arm']:r for r in nulls if r['family']==family};q2={r['arm']:r for r in rn if r['family']==family};q3={r['arm']:r for r in nn if r['family']==family}
  text+=f"| {family} | {q['v24']['actual_null_DE_calls']} | {q['v24']['null_DE_calls']} | {q['mean50_noise_mix']['null_DE_calls']} | {q2['recalibrated_population50']['null_DE_calls']} | {q3['recalibrated_noise_mix']['null_DE_calls']} |\n"
 text+='\nAll groups are unperturbed. Hundreds of synthetic DE calls versus very few actual-control calls remain an unresolved calibration problem; reduced counts do not imply correct treatment recovery. Large displayed expression-error ratios reflect small reference denominators in some contexts and are not official scaled leaderboard scores.\n\n'
 text+='## Actual-cell mean errors\n\n| Arm | Relative squared LFC error change | Relative absolute LFC error change | Relative CPM squared error change |\n|---|---:|---:|---:|\n'
 for arm,s in means.items():text+=f"| {arm} | {(s['mse_ratio']-1)*100:+.4f}% | {(s['mae_ratio']-1)*100:+.4f}% | {(s['cpm_ratio']-1)*100:+.4f}% |\n"
 text+='\n## Evidence and limitations\n\n'
 text+=f"Five boundary tests passed (including exactly-once knockdown during recalibration). Independently replayed {audit['DE_set_comparisons']+extaudit['DE_set_comparisons']+noiseaudit['DE_set_comparisons']} DE-set comparisons and {audit['profile_panels']+extaudit['profile_panels']+noiseaudit['profile_panels']} population-profile panels, with source-control/reference/treatment separation, array hashes, effect bounds and own-target masking checked. Generation asserted every row library exactly preserved. Raw Wilcoxon and reach ordering remain inherited from the frozen scorer.\n\n"
 text+='A sixth corrective mixture was frozen separately after initial H1/HCT116 panels exposed broken moment matching after row selection. It applies the original dual-moment calibration after the same balanced mixture. This is outcome-known testing with unchanged thresholds. Exact effect aliases and all25 repeated V24 baseline count diagnostics were independently verified; it is not an isolated proof that model ensembling, rather than recalibration, causes any gain. The uncalibrated row-mixture mean-error row uses its half-logFC target as a proxy; full count metrics measure the emitted population.\n\n'
 text+='A seventh final corrective candidate combines the fixed noise mixture with the same final dual-moment calibration, without row mixing. It was frozen after the other six outcomes were known, with identical thresholds and an early control-null gate. If that gate fails, treatment scoring and confirmation are not opened for this candidate. All earlier results remain unchanged.\n\n'
 text+='Initial prefix-based null sampling was corrected to the original frozen three-way control split after mean and initial null outcomes were known. Original null records and results are preserved; the amendment records this sampling-only change. Treatment counts, blend coefficients and guards were unchanged. Zero-response arms with identical emitters share the same numeric null arrays; the noise blend was generated separately.\n\n'
 text+='The public mean/source biology and prior private V24/V36 scores were already known. Development reuses 20 count cases; the remaining 24 raw cases were reserved for count-level confirmation only if a blend passed. There is no fresh lineage, no new independent biological experiment, and no evidence of statistical significance from these small diagnostic differences. D3 remains unopened. No automatic upload or Git push.\n\n'
 text+='Safe outcome-known replay without overwriting originals: `reproduce_v37.py`. Full formulas, guards and commands: [V37_ENSEMBLE_STRATEGY.md](V37_ENSEMBLE_STRATEGY.md). Saved plan, caches, means, nulls, count diagnostics, selection and independent audit are in `experiments/v37_component_ensemble/analysis/`. Root scripts: `v37_ensemble.py`, `run_v37_ensemble.py`, `test_v37_ensemble.py`, `audit_v37_ensemble.py`, `document_v37_results.py`. Large arrays and count files remain outside Git.\n'
 (ROOT/'V37_ENSEMBLE_RESULTS.md').write_text(text,encoding='utf-8',newline='\n')
 start='<!-- v37-ensemble:begin -->';end='<!-- v37-ensemble:end -->';status='None of seven tested combinations passed the applicable fixed guards; V24 remains preferred and no duplicate submission was created.' if selected=='v24' else ('Selected '+selected+'; '+('additional count confirmation passed.' if promoted else 'V24 remains preferred pending or following confirmation.'))
 entry=f'{start}\n\n## V37 component ensembles\n\nTested five fixed primary combinations and two separately frozen corrective calibrations of V24, retrieval-contrast and population-trained responses, including sign consensus, centered corrections, population mixing and empirical-noise mixing. {status} Five boundary tests passed, plus independent arithmetic and source-boundary replay. [Results](V37_ENSEMBLE_RESULTS.md), [strategy and reproducible commands](V37_ENSEMBLE_STRATEGY.md). Original model and experiment receipts remain unchanged; D3 remains sealed.\n\n{end}\n'
 p=ROOT/'README.md';s=p.read_text(encoding='utf-8')
 if start in s:
  i=s.index(start);j=s.index(end,i)+len(end);s=s[:i]+entry+s[j:]
 else:s=s.replace('# Virtual Cell Challenge 2026\n','# Virtual Cell Challenge 2026\n\n'+entry,1)
 p.write_text(s,encoding='utf-8');print('V37 results and README updated',flush=True)

if __name__=='__main__':main()
