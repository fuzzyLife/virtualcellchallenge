"""Write exact completed V79-V81 numerical results; never infer hosted scores."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
names={79:'v79_count_innovations',80:'v80_hurdle_controls',81:'v81_residual_temperature'}
for v,name in names.items():
 a=ROOT/'experiments'/name/'analysis';result=json.loads((a/'development_0_assessment.json').read_text());audit=json.loads((a/'development_0_independent_audit.json').read_text());assert audit['passed'] and result['selected'] is None
 suffix='_opposite' if v==81 else '_shuffled';rows=[]
 for arm,r in result['summary'].items():
  rows.append(f"|{arm}|{r['fidelity']:.9f}|{r['jaccard']:.9f}|{r['reach']:.9f}|{r['nmae_available_mean']:.9f}|{r['fallback_fraction']:.3f}|")
 negative='opposite noise expansion' if v==81 else 'expression-bin gene/target-shuffled shape evidence'
 scope='V81 additionally prespecified a control-only false-DE guard; it was conditional on development success and was not run. A lower-noise proposal is a calibration hypothesis, not a claim of biological realism.' if v==81 else 'Query-fit controls were used to estimate shape; query-treated and held-out reference-control values did not fit the model.'
 text=f'''# V{v} results: original development criteria fail; no new-model release

The protocol and executed code are recorded in`experiments/{name}/analysis/plan.json` and`execution_protocol.json`. Ten reused public development target/context populations (first two distinct targets in each of five contexts) compare V66 with each genuine arm and its {negative} control. Study means are balanced equally across H1, Orion and CD4. All raw profiles and discrimination are exactly V66; those invariants do not protect expression sampling corrections or DE means.

|Arm|Raw FID|Raw JAC|Raw REACH|Available-study raw NMAE|Fallback fraction|
|---|---:|---:|---:|---:|---:|
'''+ '\n'.join(rows)+f'''

No arm passes all original checks, so selection is null. No remaining-target confirmation, second seed, production population or upload file was generated for V{v}. These are unscaled public metrics, not the hosted anchor-scaled values. All public biology was reused; the tests do not establish a private leaderboard score.

The independent audit verifies all {audit['comparison_records']} saved baseline/candidate/control comparisons and {audit['manual_rank_coordinates']} raw tied-rank coordinates, replays all four DE statistics and all BH arrays, recomputes PDS and verifies profile identity. Every fresh V66 baseline matched the frozen earlier cache. Fallback populations are original V66 and are not evidence of successful new shape calibration.

{scope}

Reproduce the saved audit using the pinned scientific environment:

```cmd
python -X utf8 audit_shape_trials.py --version {v} --stage development
```

The full per-arm Boolean checks and per-study values are preserved in`development_0_assessment.json`; the rejection is not a post hoc choice of a favorable metric. See`MARGIN_VARIANCE_STRATEGY.md` for mechanisms, original thresholds, count margins, dataset exclusions, runtime and limitations.
'''
 (ROOT/f'V{v}_RESULTS.md').write_text(text,encoding='utf-8')
print('Completed V79-V81 result reports')
