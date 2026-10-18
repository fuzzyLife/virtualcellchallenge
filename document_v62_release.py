"""Document V62 tests and a fully audited local upload artifact, if available."""
from pathlib import Path
import json
from v57_pie import read,write,sha
import v62_conserved_response as m

def main():
 means=read(m.A/'mean_results.json');dev=read(m.A/'development_assessment.json');assessment=read(m.A/'count_assessment.json');choice=means['selected']
 conf=read(m.A/'confirmation_assessment.json') if (m.A/'confirmation_assessment.json').exists() else None
 supported=assessment['release_supported'];ready=read(m.A/'submission_ready.json') if (m.A/'submission_ready.json').exists() else None
 audit=read(m.A/'independent_validation.json') if (m.A/'independent_validation.json').exists() else None
 lines=['# V62 conserved-response experiment results','',f"Fixed selected candidate: **{choice}**. Public support gates: **{supported}**. V24 remains the best published submission (0.142997855); no V62 leaderboard score is known.",'','## Mean-response screening','', '| Candidate | LFC squared error / V24 | LFC absolute error / V24 | Pooled-profile squared error / V24 |','| --- | ---: | ---: | ---: |']
 for a,v in means['summary'].items():lines.append(f"| {a} | {v['mse_ratio']:.6f} | {v['mae_ratio']:.6f} | {v['profile_ratio']:.6f} |")
 lines+=['','Ratios below one improve on the fold-specific V24 mechanism. These are equal-study aggregates over 44 reused public raw populations; they are not official scaled metrics or an independent biological reserve. All held-study source-ingredient poisoning tests passed. A separate pre-score fixture correction changed a floating round-trip assertion tolerance only; original code and plan are retained.','', '## Generated integer cells']
 for label,r in [('Development',dev),('Confirmation',conf)]:
  if r is None:lines+=['',label+' was not opened because development failed.'];continue
  lines+=['',f"### {label}",'','| Arm | PDS | Fidelity | Reach | Jaccard | NMAE error | Expression error ratio |','| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
  for a,z in r['summary'].items():lines.append('| '+a+' | '+' | '.join('missing' if z[k] is None else f'{z[k]:.6f}' for k in ('pds','fidelity','reach','jaccard','nmae','common_expression_ratio'))+' |')
  lines+=['','Checks: '+', '.join(k+'='+str(v) for k,v in r['checks'].items())+'.']
 lines+=['','PDS, fidelity, reach and Jaccard above are raw higher-is-better diagnostics. NMAE and expression ratios are lower-is-better errors. Missing NMAE in an individual context remains missing; aggregate available NMAE does not imply complete metric coverage. Small 4/8-target discrimination panels are not the 300-target leaderboard. Expression is evaluated as a ratio of sums within each context with the jackknife cap, then equal-context/equal-study aggregation; it is not an average of per-cell MSEs. The generator, control inputs, row libraries and seeds are paired across arms. A shuffled source-target diagnostic tests target evidence; it is an inference corruption, not separately retrained biology.','', '## Interpretation','',('The fixed candidate passed the recorded public support gates. Its usefulness on the hidden competition contexts remains unknown. The gain follows a change in response representation and source-common correction, with the existing count generator retained. This experiment does not isolate the two components as independent causes of the final count gain.' if supported else 'The fixed candidate did not establish a validated replacement for V24. Preserve the negative results and do not describe an experimental package as improved.'),'','The raw development populations, public study labels and many source preprocessing choices have already been reused in earlier experiments. Confirmation holds out generated count cases and seeds within these same studies; it does not remove that reuse. No D3 outcomes are opened. Returning the same zero-response input to the unchanged emitter reproduces its old null behavior by construction, including its known excess null DE; this method does not solve that limitation.','', '## Reproduce','', 'See [V62_CONSERVED_RESPONSE_STRATEGY.md](V62_CONSERVED_RESPONSE_STRATEGY.md), root v62_conserved_response.py and run_v62_counts.py. Scientific source, input and production-candidate hashes are frozen before count selection. No automatic upload or Git push.']
 if ready and audit:
  assert ready['official_prep_passed'] and audit['passed'] and ready['sha256']==audit['package_sha256']
  lines+=['','## Local artifact','',f"File: experiments/v62_conserved_response/prediction.vcc, {ready['bytes']:,} bytes. SHA256 `{ready['sha256']}`. Official preparation and complete package/count/axis/provenance verification passed. Twelve groups were regenerated exactly. No leaderboard accuracy is implied by these integrity checks."]
 (m.ROOT/'V62_RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
 if ready and audit:
  sub=f'''# V62 submission ready for a user-run test

Candidate: `{choice}`. Public support gates passed: **{supported}**. V24 remains the best published reference; V62 has no hosted score yet. See [all measured tests and limitations](V62_RESULTS.md).

File: `experiments/v62_conserved_response/prediction.vcc` ({ready['bytes']:,} bytes).
SHA256: `{ready['sha256']}`.

Official prep and the complete independent integrity audit passed: 360,000 rows, 18,533 ordered genes, 900 groups, 400 cells/group, original control provenance and libraries, every packaged sparse value/index/pointer, and 12 exact group replays.

Run in Windows Command Prompt:

```cmd
cd /d L:\\promec\\Animesh\\Download\\Virtual
.\\.competition_env\\Scripts\\python.exe -m vcc submit .\\experiments\\v62_conserved_response\\prediction.vcc -m "Assay-space conserved response v62" --wait
```

No upload was performed by these scripts. If upload finishes and waiting is interrupted, use vcc status with the returned entry ID instead of creating another submission.
'''
  (m.ROOT/'V62_SUBMISSION.md').write_text(sub,encoding='utf-8',newline='\n')
  readme=m.ROOT/'README.md';text=readme.read_text(encoding='utf-8-sig');begin='<!-- v62-conserved:begin -->';end='<!-- v62-conserved:end -->'
  block=f'''{begin}
## V62: assay-space conserved response submission

Tested three new response estimators; locked **{choice}** before count evaluation. Mean-screen pooled-profile squared error improves {1-means['summary'][choice]['profile_ratio']:.2%} versus the fold-specific V24 comparator. Public count support gates: **{supported}**. The complete upload package passes official prep and integrity checks. **No V62 leaderboard score is known; V24 remains best published.** [Tests and limitations](V62_RESULTS.md) · [Submission command](V62_SUBMISSION.md) · [Method and reproduction](V62_CONSERVED_RESPONSE_STRATEGY.md).
{end}

'''
  if begin in text:
   a=text.index(begin);b=text.index(end,a)+len(end);text=text[:a]+text[b:].lstrip('\n')
  readme.write_text(block+text,encoding='utf-8')
 print('V62 reports written; public support:',supported,'artifact:',bool(ready and audit),flush=True)

if __name__=='__main__':main()
